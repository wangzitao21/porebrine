"""Offline four-morphology viewer, populated only from native solver outputs."""
from pathlib import Path
import argparse
import base64
import csv
import gzip
import json
import sys
import numpy as np
from carnallite_coating import ROOT
from plot_carnallite_coating import parameters
from run_coating_comparison import summarize, STUDY, MODES
sys.path.insert(0,str(ROOT/'scripts'))
from lib.fields import COLUMNS, read_frame, frames, history, partitions, spatial_order

CHANNELS = ['Carnallite','Halite','Sylvite','epsilon','coverage_Carnallite',
            'layer_thickness_m_Carnallite','layer_transmission_Carnallite',
            'bare_m2_Carnallite','exposed_m2_Carnallite','exposed_m2_Halite','exposed_m2_Sylvite']
LABELS = dict(baseline='只有光卤石溶解', particles='分散颗粒', dense='致密薄层', porous='多孔薄层')


def dataset(name):
    run=STUDY/name; cfg=parameters(run); work=run/'foam'; hist=history(work)
    nx,ny=int(cfg['mesh.nx']),int(cfg['mesh.ny']); packed=[]; surface_curves=[]
    stamps=frames(work)
    for stamp in stamps:
        field=read_frame(work,stamp)
        surface_names=['x_m','y_m']+[k for k in CHANNELS if k not in COLUMNS]
        pieces=[]
        for part in partitions(work):
            with (part/'output/interfaces'/f'{stamp}.csv').open() as stream:
                header=stream.readline().strip().split(',')
                pieces.append(np.loadtxt(stream,delimiter=',',usecols=[header.index(k) for k in surface_names],ndmin=2))
        surf=np.concatenate(pieces);surf=surf[spatial_order(surf)]
        np.testing.assert_allclose(field[:,:2],surf[:,:2],rtol=1e-6,atol=1e-12)
        values=np.column_stack([field[:,COLUMNS.index(k)] if k in COLUMNS else surf[:,surface_names.index(k)] for k in CHANNELS])
        if not np.isfinite(values).all():
            raise RuntimeError(name+': nonfinite visualization fields')
        packed.append(values.astype('<f4'))
        area=values[:,7].sum()
        surface_curves.append([float(stamp), values[:,8].sum(), area,
                               (values[:,4]*values[:,7]).sum()/area if area else 0,
                               values[:,9].sum()+values[:,10].sum()])
    packed=np.asarray(packed,dtype='<f4')
    initial=hist['min_Carnallite'][0]; initial_area=surface_curves[0][2]
    curves=[]
    for t, effective, support, coverage, product_area in surface_curves:
        at=lambda column:float(np.interp(t,hist['time_s'],hist[column]))
        curves.append([t,100*at('min_Carnallite')/initial,
                       100*at('min_Halite')*27.1/(initial*173.7),100*at('min_Sylvite')*37.5/(initial*173.7),
                       100*coverage,100*effective/initial_area,100*product_area/initial_area,
                       at('prec_Halite'),at('prec_Sylvite')])
    mode=name.split('_')[1]
    return dict(id=name,mode=mode,label=LABELS[mode],nx=nx,ny=ny,dx=float(cfg['mesh.spacing_m'])*1e6,
                times=[float(s) for s in stamps],channels=CHANNELS,curves=curves,
                data=base64.b64encode(gzip.compress(packed.tobytes(),compresslevel=6,mtime=0)).decode('ascii'),
                values=int(packed.size),maxima=packed.max(axis=(0,1)).tolist())


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--suffix',default='_pocket')
    parser.add_argument('--refinements',nargs='+',type=int,default=[1,2,4])
    parser.add_argument('--validation',nargs='*',default=[])
    # The pre-coated control is a separate experiment, not the dynamic
    # fresh-front demonstration served by carnallite_experiment.html.
    parser.add_argument('--output',type=Path,default=STUDY/'carnallite_precoating_comparison.html')
    args=parser.parse_args()
    names=[f'coating_{m}_r{r}{args.suffix}' for r in args.refinements for m in MODES]
    rows,errors=summarize(names+args.validation)
    sets=[dataset(n) for n in names]
    # The strongest controlled comparison: initial volume fractions must be
    # identical cell by cell within every coated group, not just in total.
    for ref in args.refinements:
        initial=[read_frame(STUDY/f'coating_{m}_r{ref}{args.suffix}'/'foam','0') for m in MODES]
        for f in initial[2:]:
            np.testing.assert_array_equal(f[:,[COLUMNS.index(k) for k in ('Carnallite','Halite','Sylvite','epsilon')]],
                                          initial[1][:,[COLUMNS.index(k) for k in ('Carnallite','Halite','Sylvite','epsilon')]])
        np.testing.assert_array_equal(initial[0][:,COLUMNS.index('Carnallite')],initial[1][:,COLUMNS.index('Carnallite')])
    with (STUDY/'coating_limits.csv').open() as stream:
        limits=[{k:v if k=='morphology' else float(v) for k,v in r.items()} for r in csv.DictReader(stream)]
    payload=dict(cases=sets,summary=rows,grid=errors,limits=limits,labels=LABELS,
                 maxima=np.max([s['maxima'] for s in sets],axis=0).tolist(),
                 checks=(STUDY/'coating_tests.log').read_text()+'\n'+(STUDY/'coating_routing_checks.log').read_text())
    template=Path(__file__).with_name('coating_animation.html').read_text(encoding='utf-8')
    encoded=json.dumps(payload,ensure_ascii=False,allow_nan=False,separators=(',',':')).replace('</','<\\/')
    output=args.output
    output.write_text(template.replace('__SIMULATION_DATA__',encoded),encoding='utf-8')
    print(output,output.stat().st_size,'bytes',len(sets),'cases',flush=True)


if __name__=='__main__':
    main()
