import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from . import ROOT
from .dual_geometry import read_geometry_header

COMPONENTS = ('Na', 'K', 'Ca', 'Mg', 'Cl', 'S', 'C')
MINERALS = ('Halite', 'Sylvite', 'Carnallite', 'Polyhalite', 'Gypsum', 'Bischofite', 'Calcite')
INPUT_FILES = ('config.json', 'geometry.npz', 'initial_conditions.npz', 'dual_geometry.bin')


def default_config(nx=16, ny=8, spacing=1e-5):
    solution = dict.fromkeys(COMPONENTS, 0.)
    solution.update(Na=.01, Cl=.01)
    return {
        'schema_version': 1,
        'title': 'Reactive transport on a regular 2D grid',
        'mesh': {'nx': nx, 'ny': ny, 'spacing_m': spacing, 'thickness_m': spacing},
        'physics': {'temperature_C': 25., 'background_pressure_atm': 1.,
                    'viscosity_Pa_s': .0018, 'diffusivity_m2_s': 1e-9,
                    'gravity_m_s2': [0., 0., 0.], 'minimum_porosity': 1e-5,
                    'chemistry_pressure_mode': 'uniform', 'density_pressure_mode': 'thermodynamic'},
        'initial_solution_mol_kg_water': solution,
        'boundaries': [
            {'name': 'inlet', 'side': 'west', 'range_m': [0., ny*spacing],
             'type': 'pressure', 'value_Pa': .01, 'solution_mol_kg_water': solution.copy()},
            {'name': 'outlet', 'side': 'east', 'range_m': [0., ny*spacing],
             'type': 'pressure', 'value_Pa': 0., 'solution_mol_kg_water': solution.copy()},
        ],
        'reactions': {'rate_coefficients_mol_m2_s': dict.fromkeys(MINERALS, 0.),
                      'interface_model': 'phaseVector', 'surface_blocking': True,
                      'nucleation_seed': 2024, 'nucleation': [],
                      'product_layer': {'model': 'impermeable', 'porosity': .35,
                                        'tortuosity': 2., 'water_diffusivity_m2_s': 1e-9,
                                        'affinity_concentration_mol_m3': 1000.}},
        'simulation': {'end_time_s': .01, 'time_step_s': .002, 'adaptive_time_step': False},
        'solver': {'linear_profile': 'balanced', 'component_scheme': 'upwind',
                   'flow_tolerance': 1e-5, 'coupling_tolerance': 1e-8,
                   'conservation_tolerance': 1e-7, 'kinetic_relative_tolerance': 2e-8,
                   'kinetic_absolute_tolerance': 1e-11, 'chemistry_threads': 1,
                   'mpi_ranks': 1, 'timeout_s': 100.},
        'output': {'frame_interval_s': .005, 'make_figures': True},
    }


def fail(message):
    raise ValueError(message)


def keys(value, allowed, location, required=()):
    if not isinstance(value, dict):
        fail(f'{location} must be an object')
    unknown = set(value)-set(allowed)
    missing = set(required)-set(value)
    if unknown:
        fail(f'Unknown {location} keys: {sorted(unknown)}')
    if missing:
        fail(f'Missing {location} keys: {sorted(missing)}')


def number(value, location, minimum=None, maximum=None, positive=False, integer=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        fail(f'{location} must be a finite number')
    if integer and (not isinstance(value, int) or value < 1):
        fail(f'{location} must be a positive integer')
    if positive and value <= 0 or minimum is not None and value < minimum or maximum is not None and value > maximum:
        fail(f'{location} is outside its permitted range')


def boolean(value, location):
    if not isinstance(value, bool):
        fail(f'{location} must be true or false')


def choice(value, choices, location):
    if not isinstance(value, str) or value not in choices:
        fail(f'{location} must be one of {choices}')


def composition(value, names, location):
    keys(value, names, location, names)
    for name in names:
        number(value[name], f'{location}.{name}', minimum=0)


def schedule(value, location, width=2):
    if not isinstance(value, list) or not value:
        fail(f'{location} must be a nonempty list')
    previous = -math.inf
    for row in value:
        if not isinstance(row, list) or len(row) != width:
            fail(f'{location} rows need {width} values')
        for entry in row:
            number(entry, location)
        if row[0] < 0 or row[0] <= previous:
            fail(f'{location} times must be nonnegative and strictly increasing')
        previous = row[0]


# Optional solver controls are typed and cannot override the scientific schema.
ADVANCED = {
    'anchorThreshold': (float, 0), 'anchorStrength': (float, 0), 'anchorTransitionWidth': (float, 0),
    'poreClosureWidth': (float, 0), 'phaseDirectionWidth': (float, 0),
    'maxFlowIterations': (int, 6), 'maxCouplingIterations': (int, 6), 'flowCorrectors': (int, 1),
    'kineticMaximumAttempts': (int, 1), 'kineticIntegrator': (int, 0),
    'analyticChemicalJacobian': (bool, None), 'reuseTransportMatrix': (bool, None),
    'consistentMomentum': (bool, None), 'conservativeFlowClosure': (bool, None),
    'couplingAcceleration': (('anderson', 'aitken', 'none'), None),
    'andersonDepth': (int, 1), 'andersonDamping': (float, 0), 'pressureRelaxation': (float, 0),
    'linearizePressureEOS': (bool, None), 'cacheBoundaryProperties': (bool, None),
    'adaptiveFlowCorrectors': (bool, None), 'detectCouplingStagnation': (bool, None),
    'traceCoupling': (bool, None), 'accelerateInitialFlow': (bool, None),
    'accelerateFaceFlux': (bool, None), 'faceMomentumInterpolation': (('linear', 'harmonic'), None),
    'transportRoundoffTolerance': (float, 0), 'densityPressureTolerance': (float, 0),
    'minimumCouplingRelaxation': (float, 0), 'maximumCouplingRelaxation': (float, 0),
    'initializePressure': (bool, None), 'solidOwnerSurfacePolicy': (('fluidSide', 'reject', 'volumeAveraged'), None),
    'nucleusResolutionModel': (('plic', 'population'), None),
    'interfaceNormalModel': (('central', 'heightFit'), None),
    'surfaceRoutingHops': (int, 1),
    'liquidControlVolumeFraction': (float, 0),
    'isolatedLiquidControlVolumes': (bool, None),
    'dualMineralGeometry': (bool, None), 'dualIgnoreShielding': (bool, None),
    'dualGeometryInterval': (float, 0), 'dualMinimumLiquidFraction': (float, 0),
    'dualWaterNeckFraction': (float, 0),
    'timeRelativeTolerance': (float, 0), 'timePorosityTolerance': (float, 0),
    'timeInventoryTolerance': (float, 0), 'timeWaterTolerance': (float, 0),
    'timeVelocityTolerance': (float, 0), 'timeSurfaceTolerance': (float, 0),
    'maximumTimeStepRejections': (int, 0), 'recoverFailedTimeSteps': (bool, None),
}


def validate_config(config):
    base = default_config()
    keys(config, set(base)|{'advanced'}, 'config', set(base))
    if type(config['schema_version']) is not int or config['schema_version'] != 1:
        fail('Unsupported schema_version; expected 1')
    if not isinstance(config['title'], str) or not config['title'].strip():
        fail('title must be a nonempty string')
    m = config['mesh']
    keys(m, base['mesh'], 'mesh', base['mesh'])
    for name in ('nx', 'ny'):
        number(m[name], f'mesh.{name}', minimum=2, integer=True)
    for name in ('spacing_m', 'thickness_m'):
        number(m[name], f'mesh.{name}', positive=True)
    p = config['physics']
    keys(p, base['physics'], 'physics', base['physics'])
    number(p['temperature_C'], 'temperature_C', minimum=-45)
    for name in ('background_pressure_atm', 'viscosity_Pa_s', 'diffusivity_m2_s'):
        number(p[name], name, positive=True)
    number(p['minimum_porosity'], 'minimum_porosity', positive=True, maximum=.1)
    if not isinstance(p['gravity_m_s2'], list) or len(p['gravity_m_s2']) != 3:
        fail('gravity_m_s2 must contain three components')
    for value in p['gravity_m_s2']:
        number(value, 'gravity_m_s2')
    if p['gravity_m_s2'][2] != 0:
        fail('A 2D XY case requires zero out-of-plane gravity')
    choice(p['chemistry_pressure_mode'], ('uniform', 'flow'), 'chemistry_pressure_mode')
    choice(p['density_pressure_mode'], ('local', 'thermodynamic'), 'density_pressure_mode')
    composition(config['initial_solution_mol_kg_water'], COMPONENTS, 'initial_solution_mol_kg_water')
    ports = config['boundaries']
    if not isinstance(ports, list) or not ports:
        fail('At least one hydraulic boundary is required')
    names, intervals = set(), {side: [] for side in ('west', 'east', 'south', 'north')}
    pressure_reference = False
    for port in ports:
        required = ('name', 'side', 'range_m', 'type', 'solution_mol_kg_water')
        keys(port, set(required)|{'value_Pa', 'value_m', 'value_m3_s', 'schedule', 'minimum_porosity',
                                 'velocity_boundary', 'concentration_boundary'}, 'boundary', required)
        name = port['name']
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', name) or name in names or name in ('walls', 'frontAndBack'):
            fail('Boundary names must be unique identifiers; walls/frontAndBack are reserved')
        names.add(name)
        side = port['side']
        choice(side, intervals, 'boundary.side')
        count = m['ny'] if side in ('west', 'east') else m['nx']
        bounds = port['range_m']
        if not isinstance(bounds, list) or len(bounds) != 2:
            fail('boundary.range_m must contain start and end coordinates')
        for bound in bounds:
            number(bound, 'boundary.range_m', minimum=0, maximum=count*m['spacing_m']*(1+1e-12))
            if abs(bound/m['spacing_m']-round(bound/m['spacing_m'])) > 1e-8:
                fail('Boundary endpoints must align with grid edges')
        if bounds[0] >= bounds[1] or any(bounds[0] < end-1e-14 and bounds[1] > start+1e-14 for start, end in intervals[side]):
            fail('Boundary segments must have positive length and must not overlap')
        intervals[side].append(bounds)
        choice(port['type'], ('pressure', 'head', 'flowRate'), 'boundary.type')
        value_key = {'pressure': 'value_Pa', 'head': 'value_m', 'flowRate': 'value_m3_s'}[port['type']]
        if set(port)&{'value_Pa', 'value_m', 'value_m3_s'} != {value_key}:
            fail(f'{port["type"]} requires only {value_key}')
        number(port[value_key], value_key)
        pressure_reference |= port['type'] in ('pressure', 'head')
        composition(port['solution_mol_kg_water'], COMPONENTS, f'boundary.{name}.solution_mol_kg_water')
        if 'schedule' in port:
            schedule(port['schedule'], f'boundary.{name}.schedule')
        number(port.get('minimum_porosity', 0), 'boundary.minimum_porosity', minimum=0, maximum=1-1e-12)
        choice(port.get('velocity_boundary', 'zeroGradient'), ('zeroGradient', 'pressureInletOutletVelocity', 'pressureNormalInletOutletVelocity'), 'velocity_boundary')
        choice(port.get('concentration_boundary', 'advectiveReservoir'), ('advectiveReservoir', 'inletOutlet', 'fixedValue'), 'concentration_boundary')
    if not pressure_reference:
        fail('At least one pressure/head boundary is required to fix the pressure datum')
    r = config['reactions']
    keys(r, set(base['reactions'])|{'initial_coating'}, 'reactions', base['reactions'])
    composition(r['rate_coefficients_mol_m2_s'], MINERALS, 'rate_coefficients_mol_m2_s')
    choice(r['interface_model'], ('phaseVector', 'heterogeneous'), 'interface_model')
    boolean(r['surface_blocking'], 'surface_blocking')
    number(r['nucleation_seed'], 'nucleation_seed', minimum=0)
    if not isinstance(r['nucleation_seed'], int):
        fail('nucleation_seed must be an integer')
    layer = r['product_layer']
    keys(layer, base['reactions']['product_layer'], 'product_layer', base['reactions']['product_layer'])
    choice(layer['model'], ('impermeable', 'porous', 'crossGrid'), 'product_layer.model')
    number(layer['porosity'], 'product_layer.porosity', minimum=0, maximum=1-1e-12)
    number(layer['tortuosity'], 'product_layer.tortuosity', minimum=1)
    for name in ('water_diffusivity_m2_s', 'affinity_concentration_mol_m3'):
        number(layer[name], name, positive=True)
    if not isinstance(r['nucleation'], list):
        fail('nucleation must be a list')
    for law in r['nucleation']:
        required = ('product', 'substrate', 'prefactor_m2_s', 'contact_angle_deg', 'embryo_radius_m')
        keys(law, set(required)|{'effective_barrier', 'interfacial_energy_J_m2'}, 'nucleation law', required)
        choice(law['product'], MINERALS, 'nucleation.product')
        choice(law['substrate'], MINERALS, 'nucleation.substrate')
        if law['product'] == law['substrate']:
            fail('A nucleation product must differ from its substrate')
        for name in ('prefactor_m2_s', 'embryo_radius_m'):
            number(law[name], name, positive=True)
        number(law['contact_angle_deg'], 'contact_angle_deg', minimum=0, maximum=180)
        if len(set(law)&{'effective_barrier', 'interfacial_energy_J_m2'}) != 1:
            fail('Choose exactly one nucleation barrier model')
        for name in ('effective_barrier', 'interfacial_energy_J_m2'):
            if name in law:
                number(law[name], name, minimum=0, positive=name == 'interfacial_energy_J_m2')
    if r['interface_model'] == 'phaseVector' and ((r['nucleation'] and not config.get('advanced', {}).get('dualMineralGeometry', False)) or layer['model'] != 'impermeable'):
        fail('Nucleation and product layers require the heterogeneous interface model')
    if 'initial_coating' in r:
        coat = r['initial_coating']
        required = ('substrate', 'morphology', 'volume_m3', 'volume_shares', 'particle_radius_m', 'contact_angle_deg')
        keys(coat, required, 'initial_coating', required)
        choice(coat['substrate'], MINERALS, 'initial_coating.substrate')
        choice(coat['morphology'], ('particles', 'film'), 'initial_coating.morphology')
        composition(coat['volume_shares'], MINERALS, 'initial_coating.volume_shares')
        if abs(sum(coat['volume_shares'].values())-1) > 1e-12 or coat['volume_shares'][coat['substrate']] != 0:
            fail('Initial coating volume shares must sum to one and exclude the substrate')
        for name in ('volume_m3', 'particle_radius_m'):
            number(coat[name], 'initial_coating.'+name, positive=True)
        number(coat['contact_angle_deg'], 'initial_coating.contact_angle_deg', positive=True, maximum=180-1e-10)
        if r['interface_model'] != 'heterogeneous' or config.get('advanced', {}).get('nucleusResolutionModel') != 'population':
            fail('Initial coating requires heterogeneous interfaces and population representation')
    sim = config['simulation']
    keys(sim, set(base['simulation'])|{'minimum_time_step_s', 'maximum_time_step_s'}, 'simulation', base['simulation'])
    for name in ('end_time_s', 'time_step_s'):
        number(sim[name], name, positive=True)
    boolean(sim['adaptive_time_step'], 'adaptive_time_step')
    minimum = sim.get('minimum_time_step_s', min(1e-6, sim['time_step_s']))
    maximum = sim.get('maximum_time_step_s', sim['time_step_s'])
    number(minimum, 'minimum_time_step_s', positive=True)
    number(maximum, 'maximum_time_step_s', positive=True)
    if not minimum <= sim['time_step_s'] <= maximum:
        fail('Time step must lie between its minimum and maximum')
    s = config['solver']
    keys(s, base['solver'], 'solver', base['solver'])
    choice(s['linear_profile'], ('strict', 'balanced'), 'linear_profile')
    choice(s['component_scheme'], ('upwind', 'limitedLinear'), 'component_scheme')
    for name in ('flow_tolerance', 'coupling_tolerance', 'conservation_tolerance',
                 'kinetic_relative_tolerance', 'kinetic_absolute_tolerance', 'timeout_s'):
        number(s[name], name, positive=True)
    number(s['conservation_tolerance'], 'conservation_tolerance', maximum=1e-7)
    for name in ('chemistry_threads', 'mpi_ranks'):
        number(s[name], name, integer=True)
    if s['mpi_ranks'] > m['nx']*m['ny']:
        fail('mpi_ranks cannot exceed the number of cells')
    o = config['output']
    keys(o, base['output'], 'output', base['output'])
    number(o['frame_interval_s'], 'frame_interval_s', positive=True)
    boolean(o['make_figures'], 'make_figures')
    advanced = config.get('advanced', {})
    keys(advanced, ADVANCED, 'advanced')
    for name, value in advanced.items():
        kind, lower = ADVANCED[name]
        if kind is bool:
            boolean(value, name)
        elif isinstance(kind, tuple):
            choice(value, kind, name)
        else:
            number(value, name, minimum=lower)
            if kind is int and not isinstance(value, int):
                fail(f'{name} must be an integer')
    if advanced.get('maxCouplingIterations', 200) > advanced.get('maxFlowIterations', 2000):
        fail('maxCouplingIterations cannot exceed maxFlowIterations')
    if advanced.get('anchorTransitionWidth', .01) > 2*advanced.get('anchorThreshold', .01):
        fail('anchorTransitionWidth exceeds twice anchorThreshold')
    if advanced.get('surfaceRoutingHops', 2) > 4:
        fail('surfaceRoutingHops cannot exceed four Cartesian links')
    if advanced.get('liquidControlVolumeFraction', 0) > 1:
        fail('liquidControlVolumeFraction cannot exceed one')
    if advanced.get('dualMineralGeometry', False):
        if r['interface_model'] != 'phaseVector' or s['mpi_ranks'] != 1 or sim['adaptive_time_step']:
            fail('dualMineralGeometry requires phaseVector, one MPI rank and no step-doubling controller')
        if advanced.get('dualGeometryInterval', .01) <= 0:
            fail('dualGeometryInterval must be positive')
        if not 0 < advanced.get('dualMinimumLiquidFraction', .5) <= 1:
            fail('dualMinimumLiquidFraction must lie in (0,1]')
        if not 0 <= advanced.get('dualWaterNeckFraction', .1) <= 1:
            fail('dualWaterNeckFraction must lie in [0,1]')
        for law in r['nucleation']:
            if law['product'] not in MINERALS[:3] or law['substrate'] not in MINERALS[:3] or 'effective_barrier' not in law:
                fail('dualMineralGeometry supports effective-barrier nucleation of the first three minerals')
    return config


def read_json(path):
    def unique(pairs):
        result = {}
        for name, value in pairs:
            if name in result:
                fail(f'Duplicate JSON key: {name}')
            result[name] = value
        return result
    return json.loads(Path(path).read_text(encoding='utf-8-sig'), object_pairs_hook=unique)


def read_inputs(directory):
    directory = Path(directory)
    if directory.is_symlink():
        fail('input_data must not be a symlink')
    directory = directory.resolve()
    if directory.name != 'input_data':
        fail('Expected an input_data directory')
    unknown = {p.name for p in directory.iterdir()}-set(INPUT_FILES)
    if unknown:
        fail(f'Unrecognized input files: {sorted(unknown)}')
    for name in INPUT_FILES[:2]:
        if not (directory/name).is_file() or (directory/name).is_symlink():
            fail(f'Missing or linked input: {name}')
    config = validate_config(read_json(directory/'config.json'))
    m = config['mesh']
    dual = directory/'dual_geometry.bin'
    if config.get('advanced', {}).get('dualMineralGeometry', False):
        if not dual.is_file() or dual.is_symlink():
            fail('Dual mineral geometry requires input_data/dual_geometry.bin')
        header = read_geometry_header(dual)
        if (header.nx, header.ny) != (m['nx'], m['ny']) or not math.isclose(header.spacing, m['spacing_m'], rel_tol=1e-12) or not math.isclose(header.depth, m['thickness_m'], rel_tol=1e-12) or header.epsilon != config['physics']['minimum_porosity']:
            fail('Fine geometry must match the main mesh and minimum porosity')
    elif dual.exists():
        fail('Enable dualMineralGeometry when supplying dual_geometry.bin')
    shape = (m['ny'], m['nx'])
    with np.load(directory/'geometry.npz', allow_pickle=False) as geometry:
        if set(geometry.files) != {'solid_fractions'}:
            fail('geometry.npz requires only solid_fractions in fixed mineral order plus inert')
        source = geometry['solid_fractions']
        if source.dtype.kind not in 'fiu':
            fail('Solid fractions must use real numeric arrays')
        solid = source.astype(float)
    if solid.shape != (8, *shape) or not np.isfinite(solid).all() or np.any(solid < 0):
        fail(f'solid_fractions must be finite and nonnegative with shape {(8, *shape)}')
    epsilon = 1-solid.sum(axis=0)
    if np.any(epsilon < config['physics']['minimum_porosity']-1e-12):
        fail('Solid fractions exceed the minimum-porosity bound')
    fields = {}
    initial = directory/'initial_conditions.npz'
    if initial.exists():
        if initial.is_symlink():
            fail('initial_conditions.npz must not be a symlink')
        allowed = set(COMPONENTS)|{'pressure_Pa', 'velocity_m_s'}
        with np.load(initial, allow_pickle=False) as data:
            if not data.files or set(data.files)-allowed:
                fail('Unknown or empty initial_conditions.npz fields')
            if any(data[name].dtype.kind not in 'fiu' for name in data.files):
                fail('Initial fields must use real numeric arrays')
            fields = {name: data[name].astype(float) for name in data.files}
        if set(fields)&set(COMPONENTS) and not set(COMPONENTS).issubset(fields):
            fail('Spatial solution fields must provide all seven component molalities')
        for name, values in fields.items():
            expected = (*shape, 2) if name == 'velocity_m_s' else shape
            if values.shape != expected or not np.isfinite(values).all():
                fail(f'{name} must be finite with shape {expected}')
            if name in COMPONENTS and np.any(values < 0):
                fail('Initial component molalities cannot be negative')
    return config, solid, fields


def input_hashes(directory):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(directory).iterdir()) if p.is_file()}


def case_path(reference, create=False):
    path = Path(reference).expanduser()
    if not path.is_absolute() and len(path.parts) == 1:
        path = ROOT/'cases'/path
    path = path.resolve()
    base = (ROOT/'cases').resolve()
    if not path.is_relative_to(base) or path == base or path.relative_to(base).parts[0] == 'trash':
        fail('A case must be under cases/ and outside trash/')
    if set(path.relative_to(base).parts)&{'input_data', 'output', '.work'}:
        fail('A case cannot be created inside another case\'s input or output folders')
    if create and path.exists():
        fail(f'Case already exists; choose a new name: {path}')
    return path


def result_path(reference, run_name=None):
    path = Path(reference).expanduser()
    if not path.exists() and len(path.parts) == 1:
        path = ROOT/'cases'/path
    path = path.resolve()
    if (path/'input_data').is_dir():
        output = path/'output'
        if run_name:
            if not re.fullmatch(r'run_[0-9]{4,}', run_name):
                fail('Run name must be run_ followed by at least four digits')
            path = output/run_name
        else:
            runs = sorted((p for p in output.glob('run_[0-9]*') if p.is_dir() and re.fullmatch(r'run_[0-9]{4,}', p.name)), key=lambda p: int(p.name[4:]))
            complete = [p for p in runs if (p/'logs/status.json').is_file()
                        and json.loads((p/'logs/status.json').read_text())['state'] == 'complete']
            if not complete:
                fail('No completed run exists for this case')
            path = complete[-1]
    if not (path/'data/summary.csv').is_file():
        fail(f'No completed result data: {path}')
    return path
