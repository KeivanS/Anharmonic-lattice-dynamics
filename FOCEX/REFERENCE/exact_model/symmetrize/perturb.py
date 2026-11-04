import sys, numpy as np
# perturb atom coordinates of structure.params (conventional reduced units) by up to +-amp Angstrom per axis length
src,dst,amp,seed=sys.argv[1],sys.argv[2],float(sys.argv[3]),int(sys.argv[4])
rng=np.random.default_rng(seed)
l=open(src).read().split('\n'); abc=[float(v) for v in l[0].split()[:3]]
nat=int(l[10].split()[0]); start=19
out=list(l)
for i in range(nat):
    w=l[start+i].split(); x=np.array([float(v) for v in w[2:5]])
    x+=rng.uniform(-amp,amp,3)/np.array(abc)
    out[start+i]='%s %s %15.10f %15.10f %15.10f'%(w[0],w[1],*x)
open(dst,'w').write('\n'.join(out))
