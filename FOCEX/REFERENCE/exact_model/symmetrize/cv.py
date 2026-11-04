import numpy as np, sys
import scipy.linalg as sl
def _svd(M,full_matrices=False): return sl.svd(M,full_matrices=full_matrices,lapack_driver='gesvd')
d=sys.argv[1]; sizes=[tuple(map(int,s.split('x'))) for s in sys.argv[2].split(',')]  # nsnap x rows_per_snap
nfold=int(sys.argv[3]); cut=1e-10
h=np.fromfile(d+'/amat_dump.bin',dtype=np.int32,count=6); nh,nt,nr,nH,nf,nc=[int(v) for v in h]
a=np.fromfile(d+'/amat_dump.bin',dtype=np.float64,offset=24); o=0
C=a[o:o+nh*nc].reshape(nc,nh).T; o+=nh*nc; bC=a[o:o+nh]; o+=nh
A=a[o:o+nf*nc].reshape(nc,nf).T; o+=nf*nc; b=a[o:o+nf]
bCs=bC.copy()
snap=np.concatenate([np.repeat(np.arange(n)+sum(m for m,_ in sizes[:i]),r) for i,(n,r) in enumerate(sizes)])
assert len(snap)==nf, (len(snap),nf)
ns=snap.max()+1; rng=np.random.default_rng(1); fold=rng.permutation(ns)%nfold
sets={'T':np.r_[0:nt],'TR':np.r_[0:nt+nr],'TH':np.r_[0:nt,nt+nr:nh],'TRH':np.r_[0:nh],'TR0':np.r_[0:nt+nr],'TR0H':np.r_[0:nh]}
bC0=bC.copy(); bC0[nt:nt+nr]=0
def svdsolve(M,v):
    u,s,vt=_svd(M); m=s>cut*s.max(); return vt[m].T@((u[:,m].T@v)/s[m])
def exact(rows,At,bt):
    Cr=C[rows]; xp=svdsolve(Cr,bC[rows]); u,s,vt=_svd(Cr,True); r=(s>cut*s.max()).sum(); K=vt[r:].T
    return xp+K@svdsolve(At@K,bt-At@xp)
def soft(rows,At,bt): return svdsolve(np.vstack([C[rows],At]),np.r_[bC[rows],bt])
print('%s: %d snapshots, %d folds; held-out relative residual ||A x - b||/||b|| (%%), and training residual'%(d.split('/')[-1],ns,nfold))
for mode,fn in (('exact (enforce_inv 1)',exact),('soft (enforce_inv 0)',soft)):
    for name,rows in sets.items():
        bC=bC0 if '0' in name else bCs
        num=0;den=0;trn=0;trd=0
        for k in range(nfold):
            te=np.isin(snap,np.where(fold==k)[0]); x=fn(rows,A[~te],b[~te])
            num+=np.sum((A[te]@x-b[te])**2); den+=np.sum(b[te]**2)
        x=fn(rows,A,b); tr=np.linalg.norm(A@x-b)/np.linalg.norm(b)
        print('  %-22s rows %-4s held-out %7.3f %%   full-data fit %7.3f %%'%(mode,name,100*np.sqrt(num/den),100*tr))
