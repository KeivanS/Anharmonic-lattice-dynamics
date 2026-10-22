"""Invariance of poscar2focex.x: every structure in DIR is rewritten 3 times with a random unimodular
change of basis, a random rotation, a random origin shift and 1e-6 A of noise; structure.params
(cell, primitive vectors and atom lines) and kpbs.params must come out identical.
Run compare_pymatgen.py (or poscar2focex.x in every DIR/<name>) first: it makes the reference.

    python invariance.py /path/to/poscar2focex.x DIR
"""
import numpy as np, os, sys, subprocess, shutil
exe=sys.argv[1]; root=sys.argv[2]; rng=np.random.default_rng(7)
def readp(f):
    L=open(f).read().splitlines(); A=np.array([[float(x) for x in L[i].split()[:3]] for i in (2,3,4)])*float(L[1])
    sp=L[5].split(); cnt=list(map(int,L[6].split())); F=np.array([[float(x) for x in L[8+i].split()[:3]] for i in range(sum(cnt))])
    return A,sp,cnt,F
def writep(f,A,sp,cnt,F):
    with open(f,'w') as h:
        h.write('t\n1.0\n'+''.join('%.12f %.12f %.12f\n'%tuple(v) for v in A)+' '.join(sp)+'\n'+' '.join(map(str,cnt))+'\nDirect\n'+''.join('%.12f %.12f %.12f\n'%tuple(v) for v in F))
def rot():
    q,_=np.linalg.qr(rng.normal(size=(3,3))); return q if np.linalg.det(q)>0 else -q
def uni():
    while True:
        M=rng.integers(-1,2,size=(3,3))
        if round(abs(np.linalg.det(M)))==1: return M
def result(d):
    L=[l.split('#')[0].split() for l in open(d+'/structure.params')]
    lat=np.array(list(map(float,L[0]))); T=np.array(list(map(float,L[1])))
    K=open(d+'/kpbs.params').read().split('\n')[2:]
    at=np.array([list(map(float,l[2:5])) for l in L[19:]]); ty=[int(l[1]) for l in L[19:]]
    return lat,T,K,(at,ty)
allok=True
for c in sorted(os.listdir(root)):
    d=os.path.join(root,c); ref=result(d)
    A,sp,cnt,F=readp(d+'/POSCAR')
    oks=[]
    for trial in range(3):
        M=uni(); R=rot(); A2=(M@A)@R.T                      # new rows, rotated
        X=F@A + rng.normal(scale=1e-6,size=F.shape)         # Cartesian (old frame)
        X=(X+rng.uniform(-2,2,3))@R.T                        # shifted origin, rotated
        F2=(X@np.linalg.inv(A2))%1.0
        t=os.path.join(root+'_inv',f'{c}_{trial}'); os.makedirs(t,exist_ok=True)
        writep(t+'/POSCAR',A2,sp,cnt,F2)
        r=subprocess.run([exe,'POSCAR','-f'],cwd=t,capture_output=True,text=True)
        if r.returncode: oks.append('FAIL'); continue
        lat,T,K,at=result(t)
        same=np.allclose(lat,ref[0],atol=1e-5) and np.allclose(T,ref[1],atol=1e-8) and K==ref[2]
        # atoms: same set relative to atom 1 (origin may differ)
        same_at=np.allclose(at[0],ref[3][0],atol=2e-5) and at[1]==ref[3][1]
        oks.append('ok' if same and same_at else ('cell ok, atoms other setting' if same else 'DIFF'))
    allok&=all(o=='ok' for o in oks)
    print(f'{c:14s}',oks)
print('ALL INVARIANT' if allok else 'NOT ALL IDENTICAL')
