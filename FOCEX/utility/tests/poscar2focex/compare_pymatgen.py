"""Compare poscar2focex.x with pymatgen's KPathSetyawanCurtarolo: sub-case, path, conventional and
primitive lattice parameters, reduced k-points and |k| of every path point.

    python compare_pymatgen.py /path/to/poscar2focex.x DIR [case ...]
"""
import numpy as np, warnings, os, sys, subprocess, re
warnings.filterwarnings('ignore')
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.kpath import KPathSetyawanCurtarolo
exe=sys.argv[1]; root=sys.argv[2]; cases=sys.argv[3:] or sorted(os.listdir(root))
lab=lambda s: s.replace('\\Gamma','G').replace('\\Sigma','Sg').replace('_','')
def frame(p):
    a,b,c=p[:3]; al,be,ga=np.radians(p[3:]); cx=c*np.cos(be); cy=c*(np.cos(al)-np.cos(ga)*np.cos(be))/np.sin(ga)
    return np.array([[a,0,0],[b*np.cos(ga),b*np.sin(ga),0],[cx,cy,np.sqrt(c*c-cx*cx-cy*cy)]])
allok=True
for c in cases:
    d=os.path.join(root,c)
    r=subprocess.run([exe,'POSCAR','-f'],cwd=d,capture_output=True,text=True)
    if r.returncode: print(f'{c:14s} FORTRAN FAILED: {(r.stdout+r.stderr).strip()[-200:]}'); allok=False; continue
    m=re.search(r"Setyawan-Curtarolo lattice (\w+) \((\w+)\)",r.stdout); fname=m.group(2) if m else "?"
    L=[l.split('#')[0].split() for l in open(d+'/structure.params')]
    latp=np.array(list(map(float,L[0]))); T=np.array(list(map(float,L[1]))).reshape(3,3)
    C=frame(latp); P=T@C; Rp=2*np.pi*np.linalg.inv(P).T
    pl=Lattice(P).parameters
    K=[l.split() for l in open(d+'/kpbs.params').readlines()[3:]]
    fk=[(k[0],np.array(list(map(float,k[1:4])))) for k in K]
    s=Structure.from_file(d+'/POSCAR'); kp=KPathSetyawanCurtarolo(s)
    pname=kp.name; pc=kp.conventional.lattice.parameters; pp=kp.prim.lattice.parameters
    pts=[kp.kpath['path'][0][0]]
    for seg in kp.kpath['path']:
        for x in seg:
            if x!=pts[-1]: pts.append(x)
    Rpm=kp.prim.lattice.reciprocal_lattice.matrix
    seq_ok=[lab(x) for x in pts]==[k for k,_ in fk]
    conv_ok=np.allclose(sorted(pc[:3]),sorted(latp[:3]),rtol=1e-4) and np.allclose(pc,latp,atol=1e-3)
    prim_ok=np.allclose(pp,pl,atol=1e-3)
    red_ok=prim_ok and all(np.allclose(kp.kpath['kpoints'][x],v,atol=1e-6) for x,(_,v) in zip(pts,fk)) if seq_ok else False
    len_ok=seq_ok and all(abs(np.linalg.norm(kp.kpath['kpoints'][x]@Rpm)-np.linalg.norm(v@Rp))<1e-5 for x,(_,v) in zip(pts,fk))
    ok=(pname==fname) and seq_ok and conv_ok and len_ok
    allok&=ok
    print(f"{c:14s} {'OK  ' if ok else 'DIFF'} S-C {fname:6s} (pymatgen {pname:6s}) path {'=' if seq_ok else '!='}  conv {'=' if conv_ok else '!='}  prim {'=' if prim_ok else '!='}  reduced k {'=' if red_ok else '!='}  |k| {'=' if len_ok else '!='}")
    if not ok or not prim_ok:
        print(f"      conv F {np.round(latp,4)}\n      conv P {np.round(pc,4)}\n      prim F {np.round(pl,4)}\n      prim P {np.round(pp,4)}")
        if not seq_ok: print('      path F',[k for k,_ in fk],'\n      path P',[lab(x) for x in pts])
print('ALL OK' if allok else 'SOME DIFFER')
