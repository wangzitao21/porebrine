import json
from pathlib import Path

import numpy as np

from .case_input import read_inputs
from .fields import COLUMNS, final_state, frames, history, metrics, read_frame


def validate_result(run_directory):
    run_directory = Path(run_directory)
    config, _, _ = read_inputs(run_directory/'record/input_data')
    summary = metrics(run_directory)
    n = config['mesh']['nx']*config['mesh']['ny']
    target = config['simulation']['end_time_s']
    end = summary['physical_time_s']
    if not np.isfinite(list(summary.values())).all():
        raise ValueError('Nonfinite solver summary')
    if abs(end-target) > max(1e-10, abs(target)*1e-10):
        raise ValueError('Simulation did not reach its requested final time')
    balance = summary['max_component_balance_error']
    limit = config['solver']['conservation_tolerance']
    if balance < 0 or balance > limit:
        raise ValueError(f'Component conservation error {balance:g} exceeds {limit:g}')
    state = final_state(run_directory)
    if len(state) != n or len(np.unique(state[:, :2], axis=0)) != n:
        raise ValueError('Incomplete or duplicate final-state cells')
    if np.any(state[:, 2] < config['physics']['minimum_porosity']-1e-10) or np.any(state[:, 2] > 1+1e-10):
        raise ValueError('Final porosity is outside the configured bounds')
    if np.any(state[:, 6] <= 0) or np.any(state[:, 7:] < -1e-12):
        raise ValueError('Invalid density, water, or component/mineral inventories')
    stamps = frames(run_directory)
    times = np.array([float(stamp) for stamp in stamps])
    if len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0) or abs(times[-1]-target) > 1e-10:
        raise ValueError('Snapshot times are incomplete or not strictly increasing')
    initial, final = (read_frame(run_directory, stamp) for stamp in (stamps[0], stamps[-1]))
    if len(initial) != n or len(final) != n or not np.array_equal(initial[:, :2], final[:, :2]):
        raise ValueError('Snapshot cell coordinates changed')
    for values in (initial, final):
        phases = values[:, [COLUMNS.index(name) for name in ['epsilon', *('Halite', 'Sylvite', 'Carnallite', 'Polyhalite', 'Gypsum', 'Bischofite', 'Calcite'), 'inert']]]
        if np.min(phases) < -1e-6 or np.max(np.abs(phases.sum(axis=1)-1)) > 2e-6:
            raise ValueError('Mineral and pore volume fractions do not close')
    h = history(run_directory)
    if abs(h['time_s'][0]) > 1e-12 or abs(h['time_s'][-1]-target) > 1e-10 or np.any(np.diff(h['time_s']) <= 0):
        raise ValueError('History does not span the full run')
    return {'passed': True, 'cells': n, 'physical_time_s': end, 'frames': len(times),
            'max_component_balance_error': balance, 'conservation_limit': limit,
            'solver_wall_s': summary['solver_wall_s'],
            'snapshot_precision': 'float32', 'final_state_precision': 'float64'}
