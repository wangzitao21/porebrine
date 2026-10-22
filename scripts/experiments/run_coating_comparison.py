"""Controlled physical coating states on three regular meshes, using saltLakeFoam."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
from pathlib import Path
import subprocess
import sys
import numpy as np
from carnallite_coating import ROOT, save_rows, interface_rows
from plot_carnallite_coating import parameters
sys.path.insert(0, str(ROOT/'scripts'))
from lib.fields import history, metrics, frames, read_frame, COLUMNS

STUDY = ROOT/'cases/carnallite_coating/output'
MODES = ('baseline', 'particles', 'dense', 'porous')


def summarize(names):
    rows = []
    solver_hashes = set()
    for name in names:
        run = STUDY/name
        with (run/'program_hashes.csv').open() as stream:
            hashes = {r['path']: r['sha256'] for r in csv.DictReader(stream)}
        solver_hashes.add(hashes['bin/saltLakeFoam'])
        cfg, hist, metric = parameters(run), history(run/'foam'), metrics(run/'foam')
        surface = np.atleast_1d(np.genfromtxt(run/'interface_history.csv', delimiter=',', names=True, deletechars=''))
        mode = name.split('_')[1]
        initial = hist['min_Carnallite'][0]
        volume = lambda i: hist['min_Halite'][i]*27.1e-6+hist['min_Sylvite'][i]*37.5e-6
        row = dict(run=name, morphology=mode, dx_um=float(cfg['mesh.spacing_m'])*1e6,
                   nx=int(cfg['mesh.nx']), ny=int(cfg['mesh.ny']), dt_s=float(cfg['simulation.time_step_s']),
                   end_s=float(hist['time_s'][-1]), initial_Carnallite_mol=initial,
                   initial_product_volume_m3=volume(0), final_product_volume_m3=volume(-1),
                   final_Halite_volume_m3=hist['min_Halite'][-1]*27.1e-6,
                   final_Sylvite_volume_m3=hist['min_Sylvite'][-1]*37.5e-6,
                   dissolved_Carnallite_pct=100*(1-hist['min_Carnallite'][-1]/initial),
                   initial_coverage_pct=100*surface['coverage_Carnallite_area_weighted'][0],
                   final_coverage_pct=100*surface['coverage_Carnallite_area_weighted'][-1],
                   initial_support_m2=surface['carnallite_support_m2'][0],
                   peak_layer_thickness_um=1e6*surface['layer_thickness_max_m'].max(),
                   precipitated_Halite_mol=hist['prec_Halite'][-1], precipitated_Sylvite_mol=hist['prec_Sylvite'][-1],
                   balance_error=metric['max_component_balance_error'], wall_s=metric['solver_wall_s'])
        row.update(minimum_accepted_dt_s=metric['minimum_accepted_dt_s'],
                   maximum_accepted_dt_s=metric['maximum_accepted_dt_s'],
                   rejected_global_steps=int(metric['rejected_global_steps']))
        for stamp in frames(run/'foam'):
            f = read_frame(run/'foam', stamp)
            phase = f[:, [COLUMNS.index(k) for k in ('Carnallite','Halite','Sylvite','epsilon')]]
            if not np.isfinite(f).all() or phase.min() < -1e-10 or np.max(np.abs(phase.sum(axis=1)-1)) > 2e-7:
                raise RuntimeError(f'{name}: invalid phase inventory at {stamp}')
        if row['balance_error'] > 1e-7:
            raise RuntimeError(f'{name}: failed element/water conservation')
        if mode != 'baseline':
            np.testing.assert_allclose(volume(0), 2*np.pi*40e-6*20e-6*2e-7, rtol=2e-12)
        rows.append(row)
    if len(solver_hashes)>1:
        raise RuntimeError('Controlled comparison mixes different solver binaries')
    if rows and any(r['end_s']!=rows[0]['end_s'] for r in rows):
        raise RuntimeError('Controlled comparison mixes different final times')
    save_rows(STUDY/'coating_comparison.csv', rows)
    errors = []
    for mode in MODES:
        selected = sorted([r for r in rows if r['morphology']==mode and '_half' not in r['run']], key=lambda r:-r['dx_um'])
        for coarse, fine in zip(selected[:-1], selected[1:]):
            hc, hf = (history(STUDY/r['run']/'foam') for r in (coarse, fine))
            fc = 100*(1-hc['min_Carnallite']/hc['min_Carnallite'][0])
            ff = 100*(1-hf['min_Carnallite']/hf['min_Carnallite'][0])
            errors.append(dict(morphology=mode, coarse_dx_um=coarse['dx_um'], fine_dx_um=fine['dx_um'],
                               final_dissolution_difference_pp=fine['dissolved_Carnallite_pct']-coarse['dissolved_Carnallite_pct'],
                               max_dissolution_curve_difference_pp=float(np.max(np.abs(fc-np.interp(hc['time_s'],hf['time_s'],ff)))),
                               final_product_volume_difference_pct=100*(fine['final_product_volume_m3']/coarse['final_product_volume_m3']-1) if coarse['final_product_volume_m3'] else 0))
    save_rows(STUDY/'coating_grid_dependence.csv', errors)
    if rows:
        np.testing.assert_allclose([r['initial_Carnallite_mol'] for r in rows], rows[0]['initial_Carnallite_mol'], rtol=2e-12)
    for r in rows:
        print(f"{r['run']}: C dissolved={r['dissolved_Carnallite_pct']:.6f}% initial cover={r['initial_coverage_pct']:.4f}% balance={r['balance_error']:.3g}", flush=True)
    return rows, errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--refinements', type=int, nargs='+', default=[1,2,4])
    parser.add_argument('--modes', nargs='+', choices=MODES, default=list(MODES))
    parser.add_argument('--dt', type=float, default=.01)
    parser.add_argument('--end', type=float, default=1.5)
    parser.add_argument('--suffix', default='')
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    names = [f'coating_{mode}_r{ref}{args.suffix}' for ref in args.refinements for mode in args.modes]
    if not args.summarize_only:
        def execute_one(mode, ref):
            name=f'coating_{mode}_r{ref}{args.suffix}'
            target=STUDY/name
            if target.exists():
                cfg = parameters(target)
                if float(cfg['simulation.end_time_s']) != args.end or float(cfg['simulation.time_step_s']) != args.dt:
                    raise RuntimeError('Existing run has different controls: '+name)
                with (target/'program_hashes.csv').open() as stream:
                    hashes={r['path']:r['sha256'] for r in csv.DictReader(stream)}
                if hashes['bin/saltLakeFoam'] != hashlib.sha256((ROOT/'bin/saltLakeFoam').read_bytes()).hexdigest():
                    raise RuntimeError('Existing run used a different solver: '+name)
                metrics(target/'foam')
                return name
            cmd=[sys.executable, str(Path(__file__).with_name('carnallite_coating.py')),
                 '--name',name,'--coating-morphology',mode,'--refinement',str(ref),'--dt',str(args.dt),
                 '--end',str(args.end),'--ranks','4','--threads','1','--frame-interval','.05']
            with (STUDY/(name+'.log')).open('w') as out:
                process=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT)
            if process.returncode:
                raise RuntimeError((STUDY/(name+'.log')).read_text()[-6000:])
            return name
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            pending=[pool.submit(execute_one, mode, ref) for ref in args.refinements for mode in args.modes]
            for job in as_completed(pending):
                print('COMPLETED',job.result(),flush=True)
    summarize(names)


if __name__=='__main__':
    main()
