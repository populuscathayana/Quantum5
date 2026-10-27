"""Offline public-data verification and optional figure/report rebuild. No CP2K calls."""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from public_data import binary_matrix, signed_frequencies, displaced_hessian, replica_coordinates, restart_geometry, MASS_AU, BOHR_A
from thermal import differences

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads((ROOT/path).read_text())
def write(path,value):path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--figures',action='store_true');p.add_argument('--animations',action='store_true')
    p.add_argument('--tex',action='store_true');args=p.parse_args()
    out=ROOT/'reproduced'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
    out.mkdir(parents=True,exist_ok=False)
    manifest=load('data/manifest.json')
    check={'public_file_hashes':all(sha(ROOT/r['path'])==r['sha256'] for r in manifest['files'])}
    data=ROOT/'data/calculations';band=data/'cineb';text=(band/'calc.out').read_text()
    xyz=replica_coordinates((band/'neb-1.restart').read_text());frames=read(band/'neb-FINAL-1.xyz',index=':')
    ee=[];scfs=[]
    for i in range(1,6):
        t=(band/f'neb-BAND{i}.out').read_text()
        values=list(map(float,re.findall(r'ENERGY\| Total FORCE_EVAL.*?energy \[hartree\]\s+([-\d.]+)',t)))
        residuals=list(map(float,re.findall(r'outer SCF iter =\s*\d+ RMS gradient =\s*([\d.Ee+-]+) energy =\s*[-\d.]+\s*\n\s*outer SCF loop converged',t)))
        ee.append(values[-1]);scfs.append(len(values)==len(residuals)>0 and max(residuals)<=1e-6 and '[ABORT]' not in t)
    source=read(ROOT/'data/structures/hcp.xyz');end=read(ROOT/'data/structures/fcc.xyz')
    check.update({'band_normal_end':'PROGRAM ENDED AT' in text and 'Band task converged' in text,
        'all_band_SCFs_converged':all(scfs),'band_shape':xyz.shape==(5,65,3),
        'band_fixed_Cu':float(np.max(np.abs(xyz[:,:64]-source.positions[:64])))<1e-8,
        'band_endpoint_geometries':np.max(np.abs(xyz[0]-source.positions))<1e-8 and np.max(np.abs(xyz[-1]-end.positions))<1e-8,
        'band_final_export':np.max(np.abs(xyz-np.array([f.positions for f in frames])))<1e-9,
        'band_peak_interior':int(np.argmax(ee))==2 and ee[2]>max(ee[0],ee[-1])})
    for field in ['RMS DISPLACEMENT','MAX DISPLACEMENT','RMS FORCE','MAX FORCE']:
        v,threshold,flag=re.findall(field+r'\s*=\s*([\d.Ee+-]+)\s*\[\s*([\d.Ee+-]+)\s*\]\s*\(\s*(YES|NO)\)',text)[-1]
        check['band_'+field]=flag=='YES' and float(v)<=float(threshold)
    modes={};frequencies={};freq_results={}
    tangent=xyz[3,64]-xyz[1,64];tangent/=np.linalg.norm(tangent)
    for site,index in [('hcp',0),('ts',2),('fcc',4)]:
        d=data/(site+'-frequency');t=(d/'calc.out').read_text();m=(d/'frequency-VIBRATIONS-1.mol').read_text()
        h,_,_=displaced_hessian(t);stored=binary_matrix(d/'frequency-Hessian-1.hess')
        f=np.fromstring(re.search(r'VIB\|Frequency \(cm\^-1\)([^\n]+)',t)[1],sep=' ')
        geom=restart_geometry((d/'calc.inp').read_text())
        vectors=np.array([np.array([line.split() for line in block.strip().splitlines()],float)
            for block in re.split(r'vibration\s+\d+\s*\n',m.split('[FR-NORM-COORD]')[1],flags=re.I)[1:]])
        modecoords=np.array([row.split()[1:] for row in re.search(r'\[FR-COORD\](.*?)\[FR-NORM-COORD\]',m,re.S)[1].strip().splitlines()],float)*BOHR_A
        base=re.search(r'(?:Minimum Structure|Transition State).*?Total energy:\s+([-\d.]+).*?VIB\|\s+O\s+([^\n]+)',t,re.S)
        force=np.fromstring(base[2],sep=' ')
        replicas=list(d.glob('frequency-r-*.out'));scf_ok=[]
        for r in replicas:
            rt=r.read_text();ne=len(re.findall(r'ENERGY\| Total FORCE_EVAL',rt))
            rr=re.findall(r'outer SCF iter =\s*\d+ RMS gradient =\s*([\d.Ee+-]+) energy =\s*[-\d.]+\s*\n\s*outer SCF loop converged',rt)
            scf_ok.append(ne==len(rr)>0 and max(map(float,rr))<=1e-6 and '[ABORT]' not in rt and set(re.findall(r'Number of electrons:\s+(\d+)',rt))=={'710'})
        check.update({site+'_normal_end':'PROGRAM ENDED AT' in t and '[ABORT]' not in t,
            site+'_SCFs':len(replicas)==5 and all(scf_ok),
            site+'_composition':geom.get_chemical_symbols()==['Cu']*64+['O'],
            site+'_source_coordinates':np.max(np.abs(geom.positions-xyz[index]))<1e-8,
            site+'_Hessian_from_force_differences':np.max(np.abs(h-stored.T*MASS_AU/1e6))<5.1e-8,
            site+'_Hessian_frequencies':np.max(np.abs(signed_frequencies(stored)-f))<1e-5,
            site+'_curvature':sum(f<0)==int(site=='ts') and sum(f>0)==3-int(site=='ts'),
            site+'_mode_source_geometry':np.max(np.abs(modecoords-geom.positions))<3e-7,
            site+'_three_O_modes_only':vectors.shape==(3,65,3) and np.max(np.abs(vectors[:,:64]))==0,
            site+'_small_residual_force':max(abs(force))<4.5e-4})
        frequencies[site]=f.tolist();modes[site]=vectors[:,64].tolist()
        freq_results[site]={'frequencies_cm-1':f.tolist(),'force_Eh_bohr':force.tolist(),'normal_modes_O':vectors[:,64].tolist()}
        if site=='ts':check['unstable_mode_along_path']=abs(vectors[0,64]@tangent)/np.linalg.norm(vectors[0,64])>.9999
    for name,expected in [('surface',True),('symmetry-baseline',False)]:
        d=data/name;initial=restart_geometry((d/'calc.inp').read_text());final=restart_geometry((d/'slab-opt-1.restart').read_text())
        delta=float(np.linalg.norm(final.positions[:32]-initial.positions[:32],axis=1).max())
        check[name+'_constraint_result']=(delta<1e-8)==expected
    initial=(data/'failed-initial/calc.out').read_text()
    check['historical_failure_not_misrepresented']='[ABORT]' in initial and 'outer SCF loop FAILED' in initial
    thermo=differences(dict(zip(['hcp','ts','fcc'],np.array(ee)[[0,2,4]])),frequencies)
    reference=load('analysis/final-results-20261027/results.json')['thermochemistry']
    check['same_electronic_reference']=abs(thermo['electronic_forward_barrier_eV']-reference['electronic_forward_barrier_eV'])<1e-12
    check['same_local_free_energy']=abs(thermo['local_harmonic_forward_barrier_eV']-reference['local_harmonic_forward_barrier_eV'])<1e-12
    check={k:bool(v) for k,v in check.items()}
    result={'accepted':all(check.values()),'checks':check,'manifest_files':len(manifest['files']),
        'energies_Eh':ee,'frequencies':freq_results,'thermochemistry':thermo,
        'scope':'Offline scientific checks of public files; original private scheduler/provenance checks are retained as historical audits, not rerun here.'}
    write(out/'verification.json',result)
    assert result['accepted'],{k:v for k,v in check.items() if not v}
    if args.figures or args.animations:
        work=out/'figure-work';work.mkdir()
        for directory in ['data','analysis','scripts']:shutil.copytree(ROOT/directory,work/directory)
        (work/'experiment-report/figures').mkdir(parents=True)
        command=[sys.executable,str(work/'scripts/build_results.py'),'--output','rebuilt-results']
        if not args.animations:command.append('--skip-animations')
        proc=subprocess.run(command,cwd=work,capture_output=True,text=True)
        (out/'figure-rebuild.txt').write_text(proc.stdout+proc.stderr);proc.check_returncode()
    if args.tex:
        exe=shutil.which('xelatex');assert exe,'Install XeLaTeX with ctex and Fandol, then retry --tex'
        report=out/'experiment-report';shutil.copytree(ROOT/'experiment-report',report)
        for i in range(2):
            proc=subprocess.run([exe,'-interaction=nonstopmode','-halt-on-error','oxygen_migration.tex'],cwd=report,capture_output=True,text=True)
            (out/f'xelatex-pass{i+1}.txt').write_text(proc.stdout+proc.stderr);proc.check_returncode()
    print(json.dumps({'accepted':True,'checks':len(check),'output':str(out),'figures':args.figures or args.animations,'animations':args.animations,'tex':args.tex},ensure_ascii=False))

if __name__=='__main__':main()
