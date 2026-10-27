"""Generate source-backed static scientific figures and coordinate animations."""
import io
import argparse
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from ase.io import read, write
from public_data import replica_coordinates
from thermal import differences, HARTREE_EV
from workflow import ROOT, save_json, sha

BLUE='#246B99'; GOLD='#C48624'; PINK='#A52D61'; COPPER='#BC7445'


def emit(fig,folder,name):
    fig.savefig(folder/(name+'.png'),dpi=220,facecolor='white',bbox_inches='tight')
    fig.savefig(folder/(name+'.pdf'),facecolor='white',bbox_inches='tight')
    plt.close(fig)


def scene(ax,cu,oxygen,side=False):
    if side:
        direction=np.array([-0.8660254,0.5,0])
        origin=np.array([7.9950585,3.08590407,0])
        xy=np.column_stack([(cu-origin)@direction,cu[:,2]])
        ox=np.array([(oxygen-origin)@direction,oxygen[2]])
        mask=(xy[:,0]>-1.7)&(xy[:,0]<3.1)&(cu[:,2]>26)
        ax.scatter(*xy[mask].T,s=420,color=COPPER,edgecolors='#8D4B28',alpha=.65)
        ax.scatter(*ox,s=90,color=PINK,edgecolors='white',zorder=6)
        ax.set(xlim=(-1.0,2.5),ylim=(26.05,28.05),xlabel='Along hcp → fcc (Å)',ylabel='z (Å)')
    else:
        top=cu[48:]
        ax.scatter(*top[:,:2].T,s=650,color=COPPER,edgecolors='#8D4B28',alpha=.65)
        ax.scatter(*oxygen[:2],s=95,color=PINK,edgecolors='white',zorder=6)
        ax.set(xlim=(4.6,10.1),ylim=(1.4,5.7),xlabel='x (Å)',ylabel='y (Å)')
        ax.set_aspect('equal')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='analysis/final-results-20261027');parser.add_argument('--skip-animations',action='store_true');args=parser.parse_args()
    out=ROOT/args.output;out.mkdir(parents=True,exist_ok=False)
    figs=ROOT/'experiment-report/figures'
    animation=ROOT/'output/animations';animation.mkdir(parents=True,exist_ok=True)
    if not args.skip_animations:
        for name in ['oxygen-migration.gif','migration-mode.gif']:
            assert not (animation/name).exists()
    band_path=ROOT/'analysis/course-cineb-5763/audit.json'
    band=json.loads(band_path.read_text());assert band['accepted_path']
    freq={s:json.loads((ROOT/f'analysis/course-{s}-frequency-{j}/audit.json').read_text())
          for s,j in [('hcp',5761),('fcc',5762),('ts',5764)]}
    assert freq['hcp']['accepted_minimum'] and freq['fcc']['accepted_minimum'] and freq['ts']['accepted_saddle']
    # Use one homogeneous set of final-band electronic energies. Frequency
    # baselines are an independent cold-start check, not mixed into the path.
    ee=dict(zip(['hcp','ts','fcc'],np.array(band['energies_Eh'])[[0,2,4]]))
    modes={s:r['frequencies_cm-1'] for s,r in freq.items()}
    thermo=differences(ee,modes)
    raw=ROOT/'data/calculations/cineb'
    xyz=replica_coordinates((raw/'neb-1.restart').read_text())
    frames=[read(ROOT/'analysis/course-hcp-5759/final.xyz') for _ in range(5)]
    for a,p in zip(frames,xyz):a.positions=p
    cell=frames[0].cell.array;cu=xyz[0,:64];oxy=xyz[:,64]
    s=np.array(band['O_path_arc_length_A']);es=np.array(band['electronic_relative_eV'])
    result={'band_audit':str(band_path.relative_to(ROOT)),'band_sha256':sha(band_path),
        'electronic_energies_Eh':ee,'frequencies_cm-1':modes,'thermochemistry':thermo,
        'source_energy_policy':'All electronic relative energies use the final converged band evaluations; independent frequency baseline differences are recorded separately.',
        'frequency_baseline_energy_offsets_eV':{k:(freq[k]['base_energy_Eh']-ee[k])*HARTREE_EV for k in freq},
        'stationary_structures':{},'course_example_relative_eV':{'hcp':0,'ts':(-296.651725+296.664375)*HARTREE_EV,'fcc':(-296.666948+296.664375)*HARTREE_EV},
        'course_comparison_scope':'CuO.pdf page26 tabulated energies and E/eV diagram; no vibrational corrections supplied there. Compared to electronic energy only.',
        'animation_scope':'Migration uses linear interpolation between the five computed images, no extra electronic points. Imaginary-mode animation is a visualized eigenvector, not a real oscillation.'}
    for key,i in [('hcp',0),('ts',2),('fcc',4)]:
        distances=frames[i].get_distances(64,range(64),mic=True)
        nearest=np.argsort(distances)[:3]
        result['stationary_structures'][key]={'O_A':oxy[i].tolist(),'nearest_Cu_one_based':(nearest+1).tolist(),
            'nearest_O_Cu_distances_A':distances[nearest].tolist(),
            'O_height_above_three_nearest_top_Cu_A':float(oxy[i,2]-cu[nearest,2].mean()),
            'max_O_residual_force_Eh_bohr':freq[key]['max_base_force_Eh_bohr']}
        write(out/(key+'.cif'),frames[i],format='cif')
    write(out/'five-image-path.xyz',frames,format='extxyz')
    write(out/'five-image-path.pdb',frames,format='proteindatabank')
    save_json(out/'results.json',result)
    plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','axes.spines.top':False,
        'axes.spines.right':False,'axes.titleweight':'normal','figure.facecolor':'white'})

    # Actual top-view structures and the same atom mapping in every panel.
    fig,axes=plt.subplots(1,3,figsize=(11.4,3.6),constrained_layout=True)
    for ax,(key,i,title) in zip(axes,[('hcp',0,'hcp minimum'),('ts',2,'Bridge transition state'),('fcc',4,'fcc minimum')]):
        scene(ax,cu,oxy[i]);ax.set_title(title)
        for idx in [55,56,59,60]:
            ax.text(cu[idx,0],cu[idx,1],str(idx+1),ha='center',va='center',fontsize=9)
        ax.plot(oxy[:,0],oxy[:,1],':',color='#53636D',lw=1.2)
        ax.text(.04,.95,f'ΔE = {es[i]:+.4f} eV',transform=ax.transAxes,va='top',fontsize=10)
    emit(fig,figs,'stationary-structures')

    # Discrete path energies; no invented smooth minimum-energy curve.
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.1),constrained_layout=True)
    initial=np.array(band['band_history_table'])[0,1:6]
    initial_s=np.linspace(0,np.linalg.norm(oxy[-1]-oxy[0]),5)
    axes[0].plot(initial_s,(initial-initial[0])*HARTREE_EV,'o--',color='#7B858D',label='Initial interpolated band',ms=5)
    axes[0].plot(s,es,'o-',color=BLUE,label='Converged CI-NEB',lw=1.8)
    axes[0].axhline(0,color='#BEC6CC',lw=.8)
    axes[0].annotate(f'{es[2]:.4f} eV',xy=(s[2],es[2]),xytext=(s[2]+.18,es[2]+.07),
        arrowprops={'arrowstyle':'-','color':BLUE},color=BLUE)
    axes[0].set(xlabel='Cumulative O displacement (Å)',ylabel='Electronic energy relative to hcp (eV)',title='(a) Electronic path')
    axes[0].legend(fontsize=9,loc='upper right');axes[0].set_ylim(-.1,.58)
    scene(axes[1],cu,oxy[2],True)
    direction=np.array([-0.8660254,.5,0]);u=(oxy-oxy[0])@direction
    axes[1].plot(u,oxy[:,2],'o-',color=BLUE,lw=1.8)
    for i,label in [(0,'hcp'),(2,'TS'),(4,'fcc')]:axes[1].annotate(label,(u[i],oxy[i,2]),xytext=(0,10),textcoords='offset points',ha='center')
    axes[1].set_title('(b) Oxygen height along the path')
    emit(fig,figs,'cineb-path')

    fig,ax=plt.subplots(figsize=(8.7,4.5),constrained_layout=True)
    for values,col,label,offset,style in [(thermo['electronic_relative_eV'],BLUE,'Electronic E',-.08,'-'),
        (thermo['local_harmonic_relative_eV'],GOLD,'Local E + Fvib (298.15 K)',.08,'--')]:
        yy=[values[k] for k in ['hcp','ts','fcc']]
        xx=np.array([0,1.5,3.0])+offset
        ax.plot(xx,yy,style,color=col,lw=1.3,label=label)
        for x,y in zip(xx,yy):ax.plot([x-.27,x+.27],[y,y],color=col,lw=2.7,ls=style)
    ax.set(xticks=[0,1.5,3],xticklabels=['hcp','Bridge TS','fcc'],ylabel='Relative energy (eV)',ylim=(-.16,.42),xlim=(-.5,3.5))
    ax.text(1.5,.354,f'E: {thermo["electronic_forward_barrier_eV"]:.4f}\nE + Fvib: {thermo["local_harmonic_forward_barrier_eV"]:.4f}',ha='center',fontsize=10)
    ax.text(3,-.14,f'E: {es[-1]:.4f}\nE + Fvib: {thermo["local_harmonic_relative_eV"]["fcc"]:.4f}',ha='center',fontsize=10)
    ax.axhline(0,color='#D2D7DB',lw=.8,zorder=-1);ax.legend(loc='upper left',fontsize=10)
    emit(fig,figs,'energy-ladder')

    fig,axes=plt.subplots(1,3,figsize=(11.4,3.6),constrained_layout=True)
    for i,ax in enumerate(axes):
        scene(ax,cu,oxy[2],side=i==2)
        v=np.array(freq['ts']['normal_modes_O'][i])
        if i==2:
            center=np.array([(oxy[2]-np.array([*oxy[0,:2],0]))@direction,oxy[2,2]])
            vector=np.array([v@direction,v[2]])*.42
        else:center=oxy[2,:2];vector=v[:2]*.42
        ax.annotate('',xy=center+vector,xytext=center-vector,arrowprops={'arrowstyle':'<->','color':BLUE,'lw':2})
        value=modes['ts'][i]
        ax.set_title(f'{abs(value):.2f}'+('i' if value<0 else '')+' cm⁻¹')
    emit(fig,figs,'ts-modes')

    text=(raw/'calc.out').read_text()
    fig,axes=plt.subplots(1,2,figsize=(10.4,3.5),constrained_layout=True)
    for ax,field,col,unit in zip(axes,['MAX FORCE','MAX DISPLACEMENT'],[BLUE,GOLD],['Eh / bohr','bohr']):
        vals=np.array(re.findall(field+r'\s*=\s*([\d.Ee+-]+)\s*\[\s*([\d.Ee+-]+)',text),float)
        ax.semilogy(range(1,len(vals)+1),vals[:,0],'o-',color=col,label='Whole-band maximum')
        ax.axhline(vals[-1,1],ls='--',color='#505F6A',label='Course threshold')
        ax.set(xlabel='Band iteration',ylabel=unit,title=field.title());ax.legend(fontsize=9)
    emit(fig,figs,'band-convergence')

    def gif(kind):
        rendered=[];coords=[]
        for t in range(48):
            if kind=='migration':
                f=(1-np.cos(2*np.pi*t/48))/2*4
                i=min(int(f),3);weight=f-i
                pos=oxy[i]*(1-weight)+oxy[i+1]*weight
                label='Five-image path: geometric interpolation'
            else:
                v=np.array(freq['ts']['normal_modes_O'][0]);v/=np.linalg.norm(v)
                pos=oxy[2]+.32*np.sin(2*np.pi*t/48)*v
                label='187.29i cm⁻¹: unstable migration eigenvector'
            coords.append(pos.tolist())
            fig,ax=plt.subplots(1,2,figsize=(8,3.8),constrained_layout=True)
            scene(ax[0],cu,pos);scene(ax[1],cu,pos,True)
            ax[0].plot(oxy[:,0],oxy[:,1],':',color=BLUE,lw=1)
            fig.suptitle(label,fontsize=11)
            stream=io.BytesIO();fig.savefig(stream,format='png',dpi=110);plt.close(fig)
            stream.seek(0);rendered.append(Image.open(stream).convert('RGB'))
        name='oxygen-migration' if kind=='migration' else 'migration-mode'
        rendered[0].save(animation/(name+'.gif'),save_all=True,append_images=rendered[1:],duration=90,loop=0)
        for idx in [0,6,12,18,24]:rendered[idx].save(out/f'{name}-frame-{idx:02}.png')
        save_json(out/(name+'-coordinates.json'),{'frames':coords,'Cu_moves':False,'frame_duration_ms':90,'amplitude_A':.32 if kind!='migration' else None})
    if not args.skip_animations:
        gif('migration');gif('mode')
    save_json(out/'artifacts.json',{'figures':{p.name:sha(p) for p in figs.glob('*') if p.stem in ['stationary-structures','cineb-path','energy-ladder','ts-modes','band-convergence']},
        'animations':{p.name:sha(p) for p in animation.glob('*.gif')}})
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
