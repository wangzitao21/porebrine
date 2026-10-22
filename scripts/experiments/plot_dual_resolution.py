"""Inspect explicit fine polygons and compare recorded fixed-grid runs."""
from pathlib import Path
import argparse
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path as PlotPath
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from lib.fields import history,metrics,frames,read_frame,COLUMNS
from lib.dual_geometry import read_geometry
from lib.study import LABELS, DEFAULT_MODES, reactive_water, save_rows

COLORS={-1:'#f5fafc',0:'#c6a76c',1:'#416b99',2:'#ca6a52'}
LINES={'carnallite_only':'#555555','frozen_seeds':'#aaaaaa','unblocked':'#286ca3','coupled':'#bf553b','mixed':'#387b64'}

def draw(ax,path,zoom=None,grid=False):
    cfg,polys,phases,owners=read_geometry(path)
    if zoom:
        xmin,xmax,ymin,ymax=zoom
        keep=[i for i,p in enumerate(polys) if p[:,0].max()>=xmin and p[:,0].min()<=xmax and p[:,1].max()>=ymin and p[:,1].min()<=ymax]
    else:keep=np.flatnonzero(phases>=0)
    ax.set_facecolor(COLORS[-1])
    for phase in (-1,0,1,2):
        selected=[polys[i] for i in keep if phases[i]==phase]
        if not selected:continue
        points=np.concatenate([np.vstack([p,p[0]]) for p in selected])
        codes=np.concatenate([[PlotPath.MOVETO]+[PlotPath.LINETO]*(len(p)-1)+[PlotPath.CLOSEPOLY] for p in selected])
        ax.add_patch(PathPatch(PlotPath(points,codes),facecolor=COLORS[phase],edgecolor='none',rasterized=True))
    nx,ny,r,H,*_=cfg
    ax.set(xlim=zoom[:2] if zoom else (110,210),ylim=zoom[2:] if zoom else (50,150),aspect='equal',xlabel='x (µm)',ylabel='y (µm)')
    if grid:
        for x in np.arange(0,nx*H*1e6+1e-6,H*1e6):ax.axvline(x,color='#777777',lw=.4,alpha=.4)
        for y in np.arange(0,ny*H*1e6+1e-6,H*1e6):ax.axhline(y,color='#777777',lw=.4,alpha=.4)

def diagnostic(path,output,cell=None):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
    fig,ax=plt.subplots(figsize=(6,6),layout='constrained')
    cfg,_,_,_=read_geometry(path)
    zoom=None
    if cell is not None:
        nx,ny,r,H,*_=cfg;x=(cell%nx)*H*1e6;y=(cell//nx)*H*1e6
        zoom=(x-2,x+H*1e6+2,y-2,y+H*1e6+2)
    draw(ax,path,zoom,True)
    ax.legend(handles=[Patch(facecolor=COLORS[m],label=n) for m,n in [(2,'Carnallite'),(1,'Sylvite'),(0,'Halite')]],loc='upper left')
    fig.savefig(output,dpi=200);plt.close(fig)

def comparison(study,tag,baseline_tag=None,modes=DEFAULT_MODES):
    rows=[];curves={};checks=[];solver_hash=None
    for mode in modes:
        run=study/(mode+(baseline_tag if mode=='carnallite_only' and baseline_tag is not None else tag));h=history(run/'foam');m=metrics(run/'foam')
        with (run/'parameters.csv').open() as f:parameters={r['parameter']:r['value'] for r in csv.DictReader(f)}
        with (run/'program_hashes.csv').open() as f:hashes={r['path']:r['sha256'] for r in csv.DictReader(f)}
        current=hashes['bin/saltLakeFoam']
        if solver_hash is None:solver_hash=current
        if current!=solver_hash:raise RuntimeError('Cases were run with different solver versions')
        if abs(h['time_s'][-1]-float(parameters['simulation.end_time_s']))>1e-10:raise RuntimeError('Incomplete run '+mode)
        if m['max_component_balance_error']>1e-7:raise RuntimeError('Conservation failure '+mode)
        ts=frames(run/'foam');surface=[]
        for stamp in ts:
            a=np.atleast_1d(np.genfromtxt(run/f'foam/output/dual_geometry/{stamp}.csv',delimiter=',',names=True))
            f=read_frame(run/'foam',stamp);wet=a['Carnallite_area_m2'];A=a['Carnallite_reaction_area_m2'] if 'Carnallite_reaction_area_m2' in a.dtype.names else wet;potential=wet+a['Carnallite_covered_m2'];w=A/A.sum() if A.sum() else A
            phase=f[:,[COLUMNS.index(x) for x in ('Halite','Sylvite','Carnallite')]]
            closed=a['closed_water_fraction'] if 'closed_water_fraction' in a.dtype.names else 0
            closure=np.max(np.abs(phase.sum(axis=1)+f[:,COLUMNS.index('epsilon')]+closed-1))
            if not np.isfinite(f).all() or phase.min() < -1e-9 or closure>2e-7:raise RuntimeError('Invalid phase state '+mode+' '+stamp)
            cfg,polys,labels,owners=read_geometry(run/f'foam/output/dual_geometry/{stamp}.bin')
            fine=np.zeros_like(phase);factor=(1-cfg[5])/(cfg[3]*1e6)**2
            for polygon,label,owner in zip(polys,labels,owners):
                if label<0:continue
                relative=polygon-polygon[0];following=np.roll(relative,-1,axis=0)
                area=np.sum(relative[:,0]*following[:,1]-relative[:,1]*following[:,0])/2
                fine[owner,label]+=area*factor
            mapping_error=np.max(np.abs(fine-phase))
            if mapping_error>6e-8:raise RuntimeError('Stored geometry and chemistry disagree '+mode+' '+stamp)
            checks.append({'mode':mode,'time_s':float(stamp),'phase_closure_error':closure,'fine_to_main_fraction_error':mapping_error})
            receivers=a['liquid_receiver'].astype(int) if 'liquid_receiver' in a.dtype.names else np.arange(len(a))
            receivers=np.where(receivers>=0,receivers,np.arange(len(a)))
            row={'mode':mode,'time_s':float(stamp),'Carnallite_area_m2':A.sum(),
                 'Carnallite_covered_m2':a['Carnallite_covered_m2'].sum(),
                 'Carnallite_coverage':1-wet.sum()/potential.sum() if potential.sum() else 0,
                 'surface_weighted_K_mol_kg':np.dot(w,f[receivers,COLUMNS.index('mK')]),
                 'surface_weighted_Carnallite_affinity':np.dot(w,1-10**f[receivers,COLUMNS.index('SICarnallite')]),
                 'maximum_liquid_regions':a['liquid_regions'].max()}
            local=reactive_water(run,stamp)
            if local is not None:
                row['surface_weighted_K_mol_kg']=local['mK']
                row['surface_weighted_Carnallite_affinity']=local['affinity']
            for name in ('Halite','Sylvite'):row[name+'_area_m2']=a[name+'_area_m2'].sum()
            surface.append(row)
        save_rows(run/'surface_history.csv',surface)
        curves[mode]=(h,surface)
        row={'mode':mode,'time_s':h['time_s'][-1],'Carnallite_dissolved_pct':100*(1-h['min_Carnallite'][-1]/h['min_Carnallite'][0]),
             'initial_Carnallite_mol':h['min_Carnallite'][0], 'final_Carnallite_mol':h['min_Carnallite'][-1],
             'coverage_initial':surface[0]['Carnallite_coverage'],'coverage_final':surface[-1]['Carnallite_coverage'],
             'area_final_m2':surface[-1]['Carnallite_area_m2'],'balance_error':m['max_component_balance_error'],
             'solver_wall_s':m['solver_wall_s'],'geometry_s':m['dual_geometry_s'],
             'chemistry_calls':m['dual_chemistry_calls'],'geometry_retries':m['dual_geometry_retries'],
             'geometry_volume_error':m['dual_max_volume_fraction_error']}
        row['max_isolated_liquid_volume_fraction']=m.get('dual_max_isolated_liquid_fraction',0)
        for name in ('Halite','Sylvite'):
            row['initial_'+name+'_mol']=h['min_'+name][0];row['final_'+name+'_mol']=h['min_'+name][-1]
            row['net_'+name+'_mol']=h['min_'+name][-1]-h['min_'+name][0]
            row['precipitated_'+name+'_mol']=h['prec_'+name][-1]
            row['new_nuclei_'+name]=m.get('dual_nuclei_'+name,0)
        rows.append(row)
    save_rows(study/('comparison'+tag+'.csv'),rows)
    save_rows(study/('validation'+tag+'.csv'),checks)
    np.testing.assert_allclose([r['initial_Carnallite_mol'] for r in rows],rows[0]['initial_Carnallite_mol'],rtol=2e-12)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for mode,(hist,surface) in curves.items():
        t=np.array([r['time_s'] for r in surface]);col=LINES[mode]
        axes[0,0].plot(hist['time_s'],100*(1-hist['min_Carnallite']/hist['min_Carnallite'][0]),label=LABELS[mode],color=col)
        axes[0,1].plot(t,[r['Carnallite_area_m2']*1e9 for r in surface],color=col)
        axes[1,0].plot(t,[r['surface_weighted_K_mol_kg'] for r in surface],color=col)
        axes[1,1].plot(t,[r['surface_weighted_Carnallite_affinity'] for r in surface],color=col)
    for ax,ylabel in zip(axes.flat,('Carnallite dissolved (%)','Reactive carnallite area (10⁻⁹ m²)','K in interface liquid receivers (mol/kg water)','Carnallite affinity (1 − Ω)')):
        ax.set(xlabel='Time (s)',ylabel=ylabel);ax.grid(alpha=.2)
    axes[0,0].legend(frameon=False,fontsize=8)
    fig.savefig(study/('comparison'+tag+'.png'),dpi=200);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    middle='unblocked' if 'unblocked' in modes else 'carnallite_only'
    for ax,(mode,label) in zip(axes,[('coupled','Initial geometry'),(middle,'Final: '+LABELS[middle]),('coupled','Final: with shielding')]):
        run=study/(mode+tag);stamp='0' if label.startswith('Initial') else frames(run/'foam')[-1]
        draw(ax,run/f'foam/output/dual_geometry/{stamp}.bin');ax.set_title(label,fontsize=10)
    fig.legend(handles=[Patch(facecolor=COLORS[m],label=n) for m,n in [(2,'Carnallite'),(1,'Sylvite'),(0,'Halite')]],loc='outside lower center',ncol=3,frameon=False)
    fig.savefig(study/('geometry'+tag+'.png'),dpi=220);plt.close(fig)
    for r in rows:print(r)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inspect',type=Path);p.add_argument('--output',type=Path);p.add_argument('--cell',type=int);p.add_argument('--tag',default='');p.add_argument('--baseline-tag');p.add_argument('--modes',nargs='+',choices=LABELS,default=DEFAULT_MODES)
    a=p.parse_args()
    if a.inspect:diagnostic(a.inspect,a.output,a.cell)
    else:comparison(ROOT/'cases/carnallite_coating/output/dual_resolution',a.tag,a.baseline_tag,a.modes)
