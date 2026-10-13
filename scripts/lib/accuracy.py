import json
from pathlib import Path
import tempfile

import numpy as np

from . import ROOT
from .case_input import COMPONENTS, MINERALS, read_json
from .fields import final_state, metrics


def compare(candidate, reference, global_limit=.001, local_limit=.01):
    if candidate.shape != reference.shape or not np.allclose(candidate[:, :2], reference[:, :2], rtol=0, atol=1e-15):
        raise ValueError('Comparisons require the same mesh and cell coordinates')
    rows = []
    def add(name, error, limit, normalization):
        rows.append({'metric': name, 'error': float(error), 'limit': limit,
                     'passed': bool(error <= limit), 'normalization': normalization})
    columns = [('porosity', 2), ('density', 6), ('water', 7)]
    columns += [(f'aqueous_{name}', 8+j) for j, name in enumerate(COMPONENTS)]
    columns += [(f'mineral_{name}', 15+j) for j, name in enumerate(MINERALS)]
    for name, column in columns:
        a, b = candidate[:, column], reference[:, column]
        delta = np.abs(a-b)
        add(name+'.global', delta.sum()/max(np.abs(b).sum(), 1e-30), global_limit, 'L1 difference / reference L1')
        add(name+'.local', delta.max()/max(np.abs(b).max(), 1e-30), local_limit, 'maximum difference / reference peak')
    a, b = candidate[:, 3:5], reference[:, 3:5]
    add('velocity.global', np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-30), global_limit, 'vector L2 difference / reference L2')
    add('velocity.local', np.linalg.norm(a-b, axis=1).max()/max(np.linalg.norm(b, axis=1).max(), 1e-30), local_limit, 'maximum vector difference / reference peak speed')
    a = candidate[:, 5]-candidate[:, 5].mean()
    b = reference[:, 5]-reference[:, 5].mean()
    add('pressure_shape.global', np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-30), global_limit, 'datum-removed L2 difference / reference L2')
    add('pressure_shape.local', np.max(np.abs(a-b))/max(np.max(np.abs(b)), 1e-30), local_limit, 'datum-removed maximum difference / reference peak')
    return rows


def time_refinement(run_directory):
    from prepare import create
    from run import run_case, write_json
    run_directory = Path(run_directory)
    source = run_directory/'record/input_data'
    config = read_json(source/'config.json')
    if config['simulation']['adaptive_time_step']:
        raise ValueError('Nominal time-step refinement requires adaptive_time_step=false')
    candidate = final_state(run_directory)
    report = {'passed': False, 'candidate_final_state_sha256': None,
              'reference_time_step_s': config['simulation']['time_step_s']/2,
              'global_limit': .001, 'local_limit': .01,
              'scope': 'one time-step refinement on the same grid; no spatial convergence claim'}
    from .export import sha256
    report['candidate_final_state_sha256'] = sha256(run_directory/'data/final_state.f64')
    try:
        with tempfile.TemporaryDirectory(prefix='.accuracy-check-', dir=ROOT/'cases') as temporary:
            parent = Path(temporary).resolve()
            if parent.parent != (ROOT/'cases').resolve() or not parent.name.startswith('.accuracy-check-'):
                raise ValueError('Unsafe accuracy-check directory')
            case, _ = create(parent/'reference', template=source)
            config['simulation']['time_step_s'] /= 2
            if 'maximum_time_step_s' in config['simulation']:
                config['simulation']['maximum_time_step_s'] /= 2
            if 'minimum_time_step_s' in config['simulation']:
                config['simulation']['minimum_time_step_s'] = min(config['simulation']['minimum_time_step_s'], config['simulation']['time_step_s'])
            write_json(case/'input_data/config.json', config)
            result, validation = run_case(case, figures=False)
            rows = compare(candidate, final_state(result))
            report.update(passed=all(row['passed'] for row in rows), errors=rows,
                          reference_validation=validation,
                          maximum_limit_ratio=max(row['error']/row['limit'] for row in rows))
    except BaseException as error:
        report.update(reason=str(error), error_type=type(error).__name__)
        write_json(run_directory/'logs/time_accuracy.json', report)
        raise
    write_json(run_directory/'logs/time_accuracy.json', report)
    if not report['passed']:
        raise ValueError('Time-refinement comparison failed; see logs/time_accuracy.json')
    return report
