import argparse
import contextlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import traceback

from lib.case_input import case_path, input_hashes, read_inputs
from lib.export import export
from lib.foam_input import generate
from lib.native_run import execute
from lib.validation import validate_result


def write_json(path, value):
    temporary = Path(path).with_suffix('.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    temporary.replace(path)


@contextlib.contextmanager
def case_lock(output):
    if os.name != 'posix':
        raise RuntimeError('Run through WSL Debian using scripts/wsl.ps1 or scripts/run.sh')
    import fcntl
    with (output/'.lock').open('a') as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('This case is already being prepared or run') from None
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def allocate_run(output):
    numbers = [int(p.name[4:]) for p in output.iterdir() if re.fullmatch(r'run_[0-9]{4,}', p.name)]
    path = output/f'run_{max(numbers, default=0)+1:04d}'
    path.mkdir()
    for name in ('data', 'figures', 'logs', 'record'):
        (path/name).mkdir()
    return path


def clean_work(work, run_directory):
    resolved = work.resolve()
    if work.is_symlink() or resolved.parent != run_directory.resolve() or work.name != '.work':
        raise ValueError('Unsafe work-directory cleanup')
    shutil.rmtree(work)


def run_case(reference, ranks=None, timeout=None, keep_work=False, figures=None):
    case = case_path(reference)
    extra = {p.name for p in case.iterdir()}-{'input_data', 'output'}
    if extra:
        raise ValueError(f'Case root may contain only input_data and output; found {sorted(extra)}')
    config, _, _ = read_inputs(case/'input_data')
    output = case/'output'
    if output.is_symlink() or output.resolve().parent != case:
        raise ValueError('output must be a real directory inside this case')
    output.mkdir(exist_ok=True)
    ranks = config['solver']['mpi_ranks'] if ranks is None else ranks
    timeout = config['solver']['timeout_s'] if timeout is None else timeout
    if isinstance(ranks, bool) or not isinstance(ranks, int) or not 1 <= ranks <= config['mesh']['nx']*config['mesh']['ny']:
        raise ValueError('MPI ranks must be positive and cannot exceed the cell count')
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Timeout must be finite and positive')
    with case_lock(output):
        started = time.monotonic()
        run_directory = allocate_run(output)
        work = run_directory/'.work'
        status = {'state': 'running', 'complete': False, 'case': str(case), 'run': run_directory.name,
                  'ranks': ranks, 'timeout_s': timeout}
        status_path = run_directory/'logs/status.json'
        write_json(status_path, status)
        try:
            before = input_hashes(case/'input_data')
            snapshot = run_directory/'record/input_data'
            shutil.copytree(case/'input_data', snapshot)
            if input_hashes(snapshot) != before or input_hashes(case/'input_data') != before:
                raise ValueError('Input files changed while their run snapshot was being created')
            write_json(run_directory/'record/input_hashes.json', before)
            generate(snapshot, work, ranks)
            remaining = timeout-(time.monotonic()-started)
            if remaining <= 0:
                raise subprocess.TimeoutExpired('input preparation', timeout)
            runtime = execute(work, run_directory/'logs', ranks, config['solver']['chemistry_threads'], remaining)
            export(work, run_directory)
            validation = validate_result(run_directory)
            write_json(run_directory/'logs/validation.json', validation)
            requested_figures = config['output']['make_figures'] if figures is None else figures
            if requested_figures:
                from lib.plotting import plot_result
                plot_result(run_directory)
            if input_hashes(case/'input_data') != before:
                status['input_changed_during_run'] = True
            if not keep_work:
                clean_work(work, run_directory)
            runtime.update(total_s=time.monotonic()-started, keep_work=keep_work)
            write_json(run_directory/'logs/runtime.json', runtime)
            status.update(state='complete', complete=True, elapsed_s=runtime['total_s'])
            write_json(status_path, status)
        except BaseException as error:
            state = 'timed_out' if isinstance(error, subprocess.TimeoutExpired) else 'cancelled' if isinstance(error, KeyboardInterrupt) else 'failed'
            status.update(state=state, reason=str(error), elapsed_s=time.monotonic()-started,
                          work_retained=work.exists())
            write_json(status_path, status)
            (run_directory/'logs/error.log').write_text(traceback.format_exc(), encoding='utf-8')
            raise
    return run_directory, validation


def main():
    parser = argparse.ArgumentParser(description='Run scientific inputs; each run creates a new output/run_XXXX directory.')
    parser.add_argument('case', help='Case name or path under cases/')
    parser.add_argument('--ranks', type=int, help='Override configured MPI ranks')
    parser.add_argument('--timeout', type=float, help='Wall-clock limit for input preparation, meshing, decomposition, and solving')
    parser.add_argument('--keep-work', action='store_true', help='Retain OpenFOAM files for diagnosis or advanced restart work')
    parser.add_argument('--no-figures', action='store_true', help='Skip the optional overview image')
    args = parser.parse_args()
    try:
        path, validation = run_case(args.case, args.ranks, args.timeout, args.keep_work,
                                    False if args.no_figures else None)
    except (ValueError, OSError, RuntimeError, subprocess.TimeoutExpired) as error:
        parser.exit(1, f'{error}\n')
    print(f'Completed: {path}\nCells: {validation["cells"]}; solver: {validation["solver_wall_s"]:.3f} s; conservation error: {validation["max_component_balance_error"]:.3g}')


if __name__ == '__main__':
    main()
