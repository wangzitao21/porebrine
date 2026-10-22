"""Fixed H/h carnallite study; no mesh or time-refinement runs."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import ctypes
import hashlib
import json
import math
import re
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from lib.case_input import COMPONENTS, default_config
from lib.chemistry import Chemistry
from lib.foam_input import generate
from lib.native_run import execute
from lib.fields import history, metrics
from lib.dual_geometry import read_geometry_header
from lib.study import (COMPARISON_CASES, COMPARISON_SETTINGS, DEFAULT_MODES,
                       DUAL_MINERALS, LABELS, STUDY, flatten, saturation, save_rows)

MODES = tuple(LABELS)


def geometry_runtime():
    """Record the libraries actually linked by the executable."""
    output = subprocess.check_output(['ldd', str(ROOT/'bin/saltLakeFoam')], text=True)
    libraries = set()
    for line in output.splitlines():
        if 'libgeos' in line and '=>' in line:
            location = line.split('=>', 1)[1].strip().split()[0]
            if location == 'not':
                raise RuntimeError('GEOS runtime is missing; rebuild the native solver')
            libraries.add(Path(location).resolve())
    api = next((p for p in libraries if p.name.startswith('libgeos_c.')), None)
    if api is None:
        raise RuntimeError('The native solver needs its linked GEOS runtime')
    library = ctypes.CDLL(str(api))
    library.GEOSversion.restype = ctypes.c_char_p
    return library.GEOSversion().decode(), sorted(libraries)


def run(mode, end, dt, tag, inlet_flow=None, nucleation_prefactor=1e11, salt_rate_scale=1.):
    geos_version, geos_libraries = geometry_runtime()
    target=STUDY/(mode+tag)
    target.mkdir(parents=True,exist_ok=False)
    (target/'logs').mkdir()
    subprocess.run([str(ROOT/'bin/makeDualGeometry'),str(target/'initial_geometry.bin'),
                    'seeded' if mode in ('frozen_seeds','coated') else 'mixed' if mode=='mixed' else 'bare'],check=True)
    header=read_geometry_header(target/'initial_geometry.bin')
    coarse=np.fromfile(str(target/'initial_geometry.bin')+'.coarse',dtype='<f8').reshape(header.ny,header.nx,3)
    solids=np.zeros((8,header.ny,header.nx))
    solids[:3]=np.moveaxis(coarse,-1,0)
    with Chemistry() as chemistry:
        molality=saturation(chemistry)
    config=default_config(header.nx,header.ny,header.spacing)
    config['title']='Dual geometry carnallite '+mode
    config['mesh']['thickness_m']=header.depth
    config['physics']['minimum_porosity']=header.epsilon
    solution=dict.fromkeys(COMPONENTS,0.)
    solution.update(Na=molality,Cl=molality)
    config['initial_solution_mol_kg_water']=solution
    for port in config['boundaries']:
        port['solution_mol_kg_water']=solution.copy()
    config['boundaries'][0]['value_Pa']=.02
    if inlet_flow is not None:
        config['boundaries'][0].pop('value_Pa')
        config['boundaries'][0].update(type='flowRate',value_m3_s=inlet_flow)
    config['reactions']['interface_model']='phaseVector'
    config['reactions']['rate_coefficients_mol_m2_s']['Carnallite']=.05
    if mode in ('unblocked','coupled','mixed','coated'):
        config['reactions']['rate_coefficients_mol_m2_s'].update(Halite=.02*salt_rate_scale,Sylvite=.01*salt_rate_scale)
        config['reactions']['nucleation']=[dict(product=m,substrate=substrate,prefactor_m2_s=nucleation_prefactor,
                                               contact_angle_deg=60.,embryo_radius_m=1e-7,effective_barrier=.02)
                                           for substrate in DUAL_MINERALS
                                           for m in DUAL_MINERALS if m != substrate]
    config['simulation'].update(end_time_s=end,time_step_s=dt,minimum_time_step_s=1e-7,
                                maximum_time_step_s=dt,adaptive_time_step=False)
    config['solver'].update(chemistry_threads=4,mpi_ranks=1,timeout_s=7200.,coupling_tolerance=1e-8)
    config['output'].update(frame_interval_s=min(.1,end),make_figures=False)
    config['advanced']={'dualMineralGeometry':True,'dualIgnoreShielding':mode=='unblocked',
                        'dualGeometryInterval':dt,'dualMinimumLiquidFraction':.5,'maxFlowIterations':8000}
    save_rows(target/'parameters.csv',[*flatten(config),{'parameter':'runtime.geos_version','value':geos_version},
                                     {'parameter':'geometry.curvature_regularisation','value':False},
                                     {'parameter':'geometry.round_corner_segments_per_quadrant','value':4},
                                     {'parameter':'geometry.volume_fraction_tolerance','value':1e-9},
                                     {'parameter':'geometry.volume_fraction_update_quantum','value':1e-10},
                                     {'parameter':'geometry.minimum_wet_fragment_area_h2','value':1e-10},
                                     {'parameter':'geometry.simultaneous_contact_distance_h','value':1e-3},
                                     {'parameter':'geometry.deferred_volume_conservation','value':True},
                                     {'parameter':'geometry.representation','value':'polygon rings with holes'},
                                     {'parameter':'geometry.initial_seed_count','value':24 if mode in ('frozen_seeds','coated') else 0},
                                     {'parameter':'geometry.nucleation_process','value':'area-weighted hazard with reproducible exponential event thresholds'},
                                     {'parameter':'geometry.max_regions_per_fine_cell','value':4096},
                                     {'parameter':'geometry.surface_contact_lookup','value':'bounding-box tree with exact segment predicate'},
                                     {'parameter':'chemistry.interface_reservoir','value':'local water regions with implicit reaction and finite throat exchange'},
                                     {'parameter':'chemistry.local_solver_coordinates','value':'pore volume and mineral partition in retained water'},
                                     {'parameter':'closed_water.model','value':'finite connected inventory with internal displacement, diffusion, conservative split and reopening'},
                                     {'parameter':'closed_water.new_space_mapping','value':'actual swept support, normalised per connected final liquid piece'}])
    source=[ROOT/'bin/saltLakeFoam',ROOT/'bin/makeDualGeometry',ROOT/'bin/libSaltChemistry.so',Path(__file__),
            ROOT/'scripts/lib/study.py',ROOT/'scripts/lib/dual_geometry.py',ROOT/'scripts/build_geometry.sh',
            *sorted((ROOT/'src').rglob('*.hpp')),*sorted((ROOT/'src/flow').glob('*.H')),
            ROOT/'src/flow/saltLakeFoam.C',ROOT/'src/chemistry/pitzer.dat',
            *geos_libraries]
    save_rows(target/'program_hashes.csv',({'path':p.relative_to(ROOT).as_posix() if p.is_relative_to(ROOT) else 'external/'+p.name,
                                          'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in source))
    with tempfile.TemporaryDirectory(prefix='porebrine-dual-') as temp:
        inputs=Path(temp)/'input_data'
        inputs.mkdir()
        (inputs/'config.json').write_text(json.dumps(config))
        np.savez_compressed(inputs/'geometry.npz',solid_fractions=solids)
        shutil.copy2(target/'initial_geometry.bin',inputs/'dual_geometry.bin')
        generate(inputs,target/'foam',1)
        control=target/'foam/system/controlDict'
        control.write_text(re.sub(r'writeInterval\s+[^;]+;',f'writeInterval {min(.1,end):g};',control.read_text()))
    result=execute(target/'foam',target/'logs',1,4,7200,bind_cores=False)
    save_rows(target/'runtime.csv',({'metric':k,'value':v} for k,v in result.items()))
    metric=metrics(target/'foam')
    if metric['max_component_balance_error']>1e-7:
        raise RuntimeError('Conservation failed')
    h=history(target/'foam')
    print(mode,'dissolved_pct',100*(1-h['min_Carnallite'][-1]/h['min_Carnallite'][0]),
          'balance',metric['max_component_balance_error'],'wall_s',metric['solver_wall_s'],flush=True)


def main():
    parser = argparse.ArgumentParser(
        description='Run fixed-grid mineral studies. --comparison runs the three viewer cases and exports their HTML.')
    parser.add_argument('--comparison', action='store_true',
                        help='Use the shared two-rate and single-mineral comparison cases')
    parser.add_argument('--modes', nargs='+', choices=MODES)
    parser.add_argument('--end', type=float)
    parser.add_argument('--dt', type=float, default=COMPARISON_SETTINGS['dt'])
    parser.add_argument('--tag', default='')
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--inlet-flow', type=float)
    parser.add_argument('--nucleation-prefactor', type=float,
                        help='Physical nucleation prefactor in m^-2 s^-1 for every mineral pair')
    parser.add_argument('--salt-rate-scale', type=float, default=1.,
                        help='Kinetic coefficient multiplier for Halite and Sylvite')
    args = parser.parse_args()
    if args.modes and len(set(args.modes)) != len(args.modes):
        parser.error('modes must be unique')
    if '/' in args.tag or '\\' in args.tag:
        parser.error('tag must not contain a path separator')
    if args.comparison and (args.modes or args.salt_rate_scale != 1.):
        parser.error('--comparison defines its own modes and rate multipliers')
    defaults = COMPARISON_SETTINGS if args.comparison else {
        'end': 1., 'inlet_flow': None, 'nucleation_prefactor': 1e11}
    for key in ('end', 'inlet_flow', 'nucleation_prefactor'):
        if getattr(args, key) is None:
            setattr(args, key, defaults[key])
    if args.jobs < 1:
        parser.error('jobs must be positive')
    for key in ('end', 'dt', 'salt_rate_scale'):
        if not math.isfinite(getattr(args, key)) or getattr(args, key) <= 0:
            parser.error(f'{key} must be finite and positive')
    if not math.isfinite(args.nucleation_prefactor) or args.nucleation_prefactor < 0:
        parser.error('nucleation-prefactor must be finite and nonnegative')
    if args.inlet_flow is not None and (not math.isfinite(args.inlet_flow) or args.inlet_flow <= 0):
        parser.error('inlet-flow must be finite and positive')
    if args.comparison:
        cases = [(c.mode, c.tag + args.tag, c.rate_scale) for c in COMPARISON_CASES]
    else:
        cases = [(mode, args.tag, args.salt_rate_scale) for mode in (args.modes or DEFAULT_MODES)]
    existing = [str(STUDY/(mode+tag)) for mode, tag, _ in cases if (STUDY/(mode+tag)).exists()]
    if existing:
        parser.error('Existing runs are preserved; choose a new --tag: ' + ', '.join(existing))
    missing = [name for name in ('saltLakeFoam', 'makeDualGeometry', 'libSaltChemistry.so')
               if not (ROOT/'bin'/name).is_file()]
    if missing:
        parser.error('Build the native programs first: ' + ', '.join(missing))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        jobs = [pool.submit(run, mode, args.end, args.dt, tag, args.inlet_flow,
                            args.nucleation_prefactor, scale) for mode, tag, scale in cases]
        for job in jobs:
            job.result()
    if args.comparison:
        from export_reactive_viewer import export
        export(STUDY, STUDY/f'三矿物反应运移_交互结果{args.tag}.html', tag=args.tag)


if __name__ == '__main__':
    main()
