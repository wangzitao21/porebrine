from pathlib import Path
import re

import numpy as np

from . import ROOT
from .case_input import COMPONENTS, MINERALS, read_inputs
from .chemistry import Chemistry
from .mesh import head, make_case


def literal(value):
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (int, float, np.number)):
        return format(value, '.17g')
    if isinstance(value, (list, tuple, np.ndarray)):
        return '('+' '.join(literal(v) for v in value)+')'
    return str(value)


def generate(input_directory, work, ranks=None):
    config, solids, initial = read_inputs(input_directory)
    work = Path(work)
    if work.exists() and any(work.iterdir()):
        raise ValueError('OpenFOAM work directory must be empty')
    m, p, sim, solver = (config[k] for k in ('mesh', 'physics', 'simulation', 'solver'))
    h, depth = m['spacing_m'], m['thickness_m']
    eps = 1-solids.sum(axis=0)
    field, put = make_case(work, eps, spacing=h, end=sim['end_time_s'], dt=sim['time_step_s'],
                          ports=m['ny'], thickness=depth, solver_profile=solver['linear_profile'],
                          component_scheme=solver['component_scheme'])
    # All side faces start as walls. Configured, grid-aligned segments become ports.
    mesh = work/'system/blockMeshDict'
    mesh.write_text(re.sub(r'boundary \(.*?mergePatchPairs',
        'boundary (walls {type wall;faces ((0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2));} '
        'frontAndBack {type empty;faces ((0 3 2 1) (4 5 6 7));}); mergePatchPairs',
        mesh.read_text(), flags=re.S))
    actions, patches = [], []
    length, height = m['nx']*h, m['ny']*h
    for port in config['boundaries']:
        lo, hi = port['range_m']
        boxes = {
            'west': (-h*.1, lo, -depth, h*.1, hi, 2*depth),
            'east': (length-h*.1, lo, -depth, length+h*.1, hi, 2*depth),
            'south': (lo, -h*.1, -depth, hi, h*.1, 2*depth),
            'north': (lo, height-h*.1, -depth, hi, height+h*.1, 2*depth),
        }
        box, name = boxes[port['side']], port['name']
        actions.extend([
            f'{{name {name};type faceSet;action new;source patchToFace;patches (walls);}}',
            f'{{name {name};type faceSet;action subset;source boxToFace;box {literal(box[:3])} {literal(box[3:])};}}',
        ])
        patches.append(f'{{name {name};patchInfo {{type patch;}} constructFrom set;set {name};}}')
    put('system/topoSetDict', head('topoSetDict')+'actions (\n'+'\n'.join(actions)+'\n);\n')
    put('system/createPatchDict', head('createPatchDict')+'pointSync false;patches (\n'+'\n'.join(patches)+'\n);\n')
    port_bc = ' '.join(port['name']+' {type zeroGradient;}' for port in config['boundaries'])
    scalar_bc = 'walls {type zeroGradient;} '+port_bc
    vector_bc = 'walls {type fixedValue;value uniform (0 0 0);} '+port_bc
    field('epsilon', eps, scalar_bc)
    field('p', initial.get('pressure_Pa', 0.), scalar_bc, '[1 -1 -2 0 0 0 0]')
    if 'velocity_m_s' in initial:
        velocities = np.column_stack([initial['velocity_m_s'].reshape(-1, 2), np.zeros(eps.size)])
        values = f'nonuniform List<vector>\n{len(velocities)}\n(\n'+'\n'.join(literal(row) for row in velocities)+'\n)'
    else:
        values = 'uniform (0 0 0)'
    field('U', values, vector_bc, '[0 1 -1 0 0 0 0]', 'volVectorField')
    for j in range(7):
        field(f'mineral{j}', solids[j], scalar_bc)
    field('inert', solids[7], scalar_bc)
    reservoir = []
    with Chemistry(p['temperature_C'], p['background_pressure_atm']) as chemistry:
        cache = {}
        def fractions(totals):
            key = tuple(float(v) for v in totals)
            if key not in cache:
                cache[key] = chemistry.fractions(key)[0]
            return cache[key]
        resident = fractions([config['initial_solution_mol_kg_water'][n] for n in COMPONENTS])
        if set(COMPONENTS).issubset(initial):
            totals = np.stack([initial[n].ravel() for n in COMPONENTS], axis=1)
            unique, inverse = np.unique(totals, axis=0, return_inverse=True)
            values = np.array([fractions(row) for row in unique])[inverse].reshape(*eps.shape, 7)
            for j in range(7):
                field(f'Y{j}', values[:, :, j], scalar_bc)
        else:
            for j in range(7):
                field(f'Y{j}', resident[j], scalar_bc)
        for port in config['boundaries']:
            reservoir.append(fractions([port['solution_mol_kg_water'][n] for n in COMPONENTS]))
    r = config['reactions']
    layer = r['product_layer']
    controls = {
        'temperatureC': p['temperature_C'], 'absolutePressureAtm': p['background_pressure_atm'],
        'waterPropertyPolicy': 'phreeqcCompatible', 'dx': h, 'dy': h,
        'epsilonMin': p['minimum_porosity'], 'viscosity': p['viscosity_Pa_s'],
        'gravity': p['gravity_m_s2'], 'pressureReferencePoint': [0, 0, 0],
        'chemistryPressureMode': p['chemistry_pressure_mode'], 'densityPressureMode': p['density_pressure_mode'],
        'inletChemistry': 'fixedFractions', 'inletFractions': reservoir[0],
        'rateCoefficients': [r['rate_coefficients_mol_m2_s'][name] for name in MINERALS],
        'mineralInterfaceModel': r['interface_model'], 'surfaceBlocking': r['surface_blocking'],
        'nucleationSeed': r['nucleation_seed'], 'resolvedNucleusRadius': .35*h,
        'productLayerModel': layer['model'], 'productLayerPorosity': layer['porosity'],
        'productLayerTortuosity': layer['tortuosity'], 'productLayerWaterDiffusivity': layer['water_diffusivity_m2_s'],
        'productLayerAffinityConcentration': layer['affinity_concentration_mol_m3'],
        'flowTolerance': solver['flow_tolerance'], 'couplingTolerance': solver['coupling_tolerance'],
        'conservationTolerance': solver['conservation_tolerance'],
        'kineticTolerance': solver['kinetic_relative_tolerance'], 'kineticAbsoluteTolerance': solver['kinetic_absolute_tolerance'],
        'chemistryThreads': solver['chemistry_threads'], 'analyticChemicalJacobian': True,
        'kineticMaximumAttempts': 100000, 'anchorThreshold': .01, 'anchorStrength': 1e12,
        'anchorTransitionWidth': .01, 'poreClosureWidth': 1e-8, 'phaseDirectionWidth': 1e-9,
        'maxFlowIterations': 2000, 'maxCouplingIterations': 200, 'flowCorrectors': 5,
        'couplingAcceleration': 'anderson', 'andersonDepth': 5, 'andersonDamping': .8,
        'reuseTransportMatrix': True, 'consistentMomentum': True, 'conservativeFlowClosure': True,
        'pressureRelaxation': 1, 'transportRoundoffTolerance': 1e-12,
        'adaptiveTimeStep': sim['adaptive_time_step'], 'recoverFailedTimeSteps': True,
        'maximumTimeStepRejections': 12,
        'minDeltaT': sim.get('minimum_time_step_s', min(1e-6, sim['time_step_s'])),
        'maxDeltaT': sim.get('maximum_time_step_s', sim['time_step_s']),
        'frameInterval': config['output']['frame_interval_s'],
    }
    controls.update(config.get('advanced', {}))
    text = head('reactiveProperties')+f'database "{ROOT}/src/chemistry/pitzer.dat";\n'
    text += f'diffusivity [0 2 -1 0 0 0 0] {literal(p["diffusivity_m2_s"])};\n'
    text += '\n'.join(f'{name} {literal(value)};' for name, value in controls.items())+'\n'
    text += 'hydraulicBoundaries {\n'
    for port, incoming in zip(config['boundaries'], reservoir):
        unit = {'pressure': 'value_Pa', 'head': 'value_m', 'flowRate': 'value_m3_s'}[port['type']]
        values = {'type': port['type'], 'value': port[unit], 'inletChemistry': 'fixedFractions',
                  'inletFractions': incoming, 'minimumPorosity': port.get('minimum_porosity', 0),
                  'velocityBoundary': port.get('velocity_boundary', 'zeroGradient'),
                  'concentrationBoundary': port.get('concentration_boundary', 'advectiveReservoir')}
        if 'schedule' in port:
            values['schedule'] = port['schedule']
        text += port['name']+' {'+' '.join(f'{name} {literal(value)};' for name, value in values.items())+'}\n'
    text += '}\nheterogeneousNucleation {\n'
    for j, law in enumerate(r['nucleation']):
        values = {'product': law['product'], 'substrate': law['substrate'], 'prefactor': law['prefactor_m2_s'],
                  'contactAngle': law['contact_angle_deg'], 'embryoRadius': law['embryo_radius_m']}
        if 'effective_barrier' in law:
            values.update(nucleationModel='effectiveBarrier', barrier=law['effective_barrier'])
        else:
            values.update(nucleationModel='classical', interfacialEnergyJm2=law['interfacial_energy_J_m2'])
        text += f'law{j} '+' {'+' '.join(f'{name} {literal(value)};' for name, value in values.items())+'}\n'
    text += '}\n'
    put('constant/reactiveProperties', text)
    put('system/decomposeParDict', head('decomposeParDict')+f'numberOfSubdomains {ranks or solver["mpi_ranks"]};method scotch;\n')
    return config
