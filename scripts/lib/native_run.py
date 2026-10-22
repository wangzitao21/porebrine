import hashlib
import os
from pathlib import Path
import signal
import subprocess
import time

from . import ROOT


def command(args, work, log, timeout):
    with Path(log).open('w') as out:
        process = subprocess.Popen(args, cwd=work, stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise
    if code:
        tail = '\n'.join(Path(log).read_text(errors='replace').splitlines()[-35:])
        raise RuntimeError(f'{Path(args[0]).name} failed; see {log}\n{tail}')


def mesh_fingerprint(work):
    digest = hashlib.sha256()
    for name in ('points', 'faces', 'owner', 'neighbour', 'boundary'):
        digest.update(name.encode())
        digest.update((Path(work)/'constant/polyMesh'/name).read_bytes())
    return digest.hexdigest()


def execute(work, logs, ranks, threads, timeout, executable=None, *, bind_cores=True):
    work, logs = Path(work).resolve(), Path(logs).resolve()
    executable = Path(executable or ROOT/'bin/saltLakeFoam').resolve()
    if not executable.is_file():
        raise ValueError('Solver is missing; build the project first')
    started = time.monotonic()
    def run(args, name):
        remaining = timeout-(time.monotonic()-started)
        if remaining <= 0:
            raise subprocess.TimeoutExpired('case preparation and solver', timeout)
        command(args, work, logs/name, remaining)
    run(['blockMesh'], 'mesh.log')
    if (work/'system/topoSetDict').is_file():
        run(['topoSet'], 'topology.log')
        run(['createPatch', '-overwrite'], 'patches.log')
    run(['checkMesh'], 'mesh_check.log')
    if 'Mesh OK.' not in (logs/'mesh_check.log').read_text():
        raise RuntimeError('Mesh quality check did not pass')
    if ranks == 1:
        run([str(executable)], 'solver.log')
    else:
        run(['decomposePar'], 'decomposition.log')
        cores = set()
        for topology in Path('/sys/devices/system/cpu').glob('cpu[0-9]*/topology'):
            try:
                cores.add(((topology/'physical_package_id').read_text().strip(),
                           (topology/'core_id').read_text().strip()))
            except OSError:
                pass
        binding = ['--bind-to', 'core', '--map-by', f'slot:PE={threads}'] if bind_cores and cores and ranks*threads <= len(cores) else ['--bind-to', 'none']
        run(['mpirun', *binding, '-np', str(ranks), str(executable), '-parallel'], 'solver.log')
    return {'preparation_and_solver_s': time.monotonic()-started,
            'ranks': ranks, 'chemistry_threads_per_rank': threads, 'mesh_sha256': mesh_fingerprint(work)}
