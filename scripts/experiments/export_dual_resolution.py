"""Self-contained animation of actual solver frames and fine mineral polygons."""
from pathlib import Path
import argparse
import base64
import csv
import json
import zlib
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from lib.dual_geometry import read_geometry
from lib.study import LABELS, DEFAULT_MODES, reactive_water
from lib.fields import history, metrics, frames, read_frame, COLUMNS

def export(study,tag,baseline_tag=None,modes=DEFAULT_MODES):
    data={'runs':[],'times':None,'extent':[320,200],'mesh':[64,40],
          'parameters':{'H_um':5,'h_um':1.25,'depth_um':20,'dt_s':.01,'NaCl_mol_kg':6.1292275648}}
    solver_hash=None
    particles={}
    particle_file=study/'particle_sizes.csv'
    if particle_file.exists():
        with particle_file.open() as f:
            for row in csv.DictReader(f):
                if int(row['mineral'])==1:
                    particles.setdefault((row['case'],float(row['time_s'])),[]).append([float(row['x_um']),float(row['y_um']),float(row['equivalent_diameter_um'])])
    for mode in modes:
        run=study/(mode+(baseline_tag if mode=='carnallite_only' and baseline_tag is not None else tag));metric=metrics(run/'foam');hist=history(run/'foam');stamps=frames(run/'foam')
        if any(abs(hist['min_'+m][0])>1e-30 for m in ('Halite','Sylvite')):
            raise RuntimeError('The new-surface viewer requires zero initial product seeds; use the retained historical HTML for seeded runs')
        with (run/'program_hashes.csv').open() as f:hashes={r['path']:r['sha256'] for r in csv.DictReader(f)}
        if solver_hash is None:solver_hash=hashes['bin/saltLakeFoam']
        if solver_hash!=hashes['bin/saltLakeFoam']:raise RuntimeError('Mixed solver versions')
        if metric['max_component_balance_error']>1e-7:raise RuntimeError('Failed mass conservation')
        times=[float(t) for t in stamps]
        if abs(times[-1]-metric['physical_time_s'])>1e-10:raise RuntimeError('Missing final frame')
        if data['times'] is None:data['times']=times
        elif times!=data['times']:raise RuntimeError('Case times do not match')
        record={'mode':mode,'metrics':metric,'frames':[],'curve':[]}
        for stamp in stamps:
            cfg,polys,phases,owners=read_geometry(run/f'foam/output/dual_geometry/{stamp}.bin')
            frame=read_frame(run/'foam',stamp)
            surface=np.atleast_1d(np.genfromtxt(run/f'foam/output/dual_geometry/{stamp}.csv',delimiter=',',names=True))
            wet=surface['Carnallite_area_m2'];area=surface['Carnallite_reaction_area_m2'] if 'Carnallite_reaction_area_m2' in surface.dtype.names else wet;covered=surface['Carnallite_covered_m2'];weights=area/area.sum() if area.sum() else area
            receivers=surface['liquid_receiver'].astype(int) if 'liquid_receiver' in surface.dtype.names else np.arange(len(surface))
            receivers=np.where(receivers>=0,receivers,np.arange(len(surface)))
            t=float(stamp);hi=int(np.argmin(np.abs(hist['time_s']-t)))
            curve={'time':t,'dissolved':float(100*(1-hist['min_Carnallite'][hi]/hist['min_Carnallite'][0])),
                   'area':float(area.sum()*1e9),'coverage':float(100*covered.sum()/(wet.sum()+covered.sum())) if (wet.sum()+covered.sum()) else 0,
                   'K':float(weights@frame[receivers,COLUMNS.index('mK')]),
                   'affinity':float(weights@(1-10**frame[receivers,COLUMNS.index('SICarnallite')])),
                   'halite':float((hist['min_Halite'][hi]-hist['min_Halite'][0])*1e12),
                   'sylvite':float((hist['min_Sylvite'][hi]-hist['min_Sylvite'][0])*1e12),
                   'halite_precip':float(hist['prec_Halite'][hi]*1e12),
                   'sylvite_precip':float(hist['prec_Sylvite'][hi]*1e12),
                   'nuclei_H':int(surface['nuclei_Halite'].sum()) if 'nuclei_Halite' in surface.dtype.names else 0,
                   'nuclei_S':int(surface['nuclei_Sylvite'].sum()) if 'nuclei_Sylvite' in surface.dtype.names else 0,
                   'closed_ppm':float(surface['closed_water_fraction'].mean()*1e6) if 'closed_water_fraction' in surface.dtype.names else 0.}
            local=reactive_water(run,stamp)
            if local is not None:curve.update(K=local['mK'],affinity=local['affinity'])
            geometry=[[] for _ in range(3)]
            for p,m in zip(polys,phases):
                if m>=0:geometry[m].append(np.round(p,7).ravel().tolist())
            record['frames'].append({'polygons':geometry,
                                     'K':np.round(frame[:,COLUMNS.index('mK')],7).tolist(),
                                     'Mg':np.round(frame[:,COLUMNS.index('mMg')],7).tolist(),
                                     'sylvite_markers':particles.get((run.name,t),[])})
            record['curve'].append(curve)
        data['runs'].append(record)
    data['solver_sha256']=solver_hash
    feedback=study/'reaction_feedback.csv'
    if feedback.exists():
        with feedback.open() as f:records=list(csv.DictReader(f))
        selected={}
        for mode in modes:
            suffix=mode+(baseline_tag if mode=='carnallite_only' and baseline_tag is not None else tag)
            match=next((r for r in records if r['case'].endswith('/'+suffix)),None)
            if match:selected[mode]={k:float(v) for k,v in match.items() if k!='case'}
        if 'carnallite_only' in selected and 'coupled' in selected:data['feedback']=selected
    data['Kmax']=max(max(f['K']) for r in data['runs'] for f in r['frames'])
    payload=base64.b64encode(zlib.compress(json.dumps(data,separators=(',',':'),allow_nan=False).encode(),9)).decode()
    template=Path(__file__).with_name('dual_resolution_viewer.html').read_text(encoding='utf-8')
    if template.count('__DATA__')!=1:raise RuntimeError('Animation payload placeholder mismatch')
    output=study/'双分辨率_新表面成核与运移.html'
    output.write_text(template.replace('__DATA__',payload),encoding='utf-8')
    print(output,output.stat().st_size,'bytes',len(data['times']),'stored frames per case')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tag',default='');p.add_argument('--baseline-tag');p.add_argument('--modes',nargs='+',choices=LABELS,default=DEFAULT_MODES);args=p.parse_args()
    export(ROOT/'cases/carnallite_coating/output/dual_resolution',args.tag,args.baseline_tag,args.modes)
