import numpy as np, sys
d=sys.argv[1]
h=np.fromfile(d+'/amat_dump.bin',dtype=np.int32,count=6)
nh,nt,nr,nH,nf,nc=[int(v) for v in h]
a=np.fromfile(d+'/amat_dump.bin',dtype=np.float64,offset=24)
o=0
ahom=a[o:o+nh*nc].reshape(nc,nh).T; o+=nh*nc
bhom=a[o:o+nh]; o+=nh
af=a[o:o+nf*nc].reshape(nc,nf).T; o+=nf*nc
bf=a[o:o+nf]
print('rows: transl',nt,'rot',nr,'huang',nH,'force',nf,'cols',nc)
# FC2 columns: info from fc2_irr.dat (columns 2..)
L=[l.split() for l in open(d+'/fc2_irr.dat')]
n1=sum(1 for _ in open(d+'/fc1_irr.dat'))
c=5.2724287802
dz=np.zeros(nc); dist=np.zeros(nc); isfc2=np.zeros(nc,bool)
for k,l in enumerate(L):
    j=n1+k; isfc2[j]=True; dist[j]=float(l[7]); dz[j]=abs(float(l[14]))*c
def kernel(C,cut=1e-10):
    u,s,vt=np.linalg.svd(C); r=(s>cut*s.max()).sum(); return vt[r:].T, vt[:r].T, s, u, r
nuse=nt if len(sys.argv)>2 and sys.argv[2]=="tr" else nt+nr
C=ahom[:nuse]; dC=bhom[:nuse]
K,R,s,u,r=kernel(C)
xp=np.linalg.lstsq(C,dC,rcond=1e-10)[0]
H=ahom[nt+nr:]; bH=bhom[nt+nr:]
labs=[]
vo={(0,0):1,(1,1):2,(2,2):3,(1,2):4,(0,2):5,(0,1):6}
ax='xyz'
for al in range(3):
  for be in range(al,3):
    for ga in range(3):
      for de in range(ga,3):
        if vo[(al,be)]<=vo[(ga,de)]: continue
        labs.append(f'[{ax[al]}{ax[be]},{ax[ga]}{ax[de]}]')
AK=af@K
y=np.linalg.lstsq(AK,bf-af@xp,rcond=1e-10)[0]
x0=xp+K@y
res=lambda x: np.linalg.norm(af@x-bf)/np.linalg.norm(bf)*100
print('kernel dim',K.shape[1],'exact rows',nuse,'fit %.3f %%'%res(x0))
r0=H@x0-bH
print('Huang violation of the transl+rot fit (eV), RHS, and contribution by |dz| of the pair:')
bins=[(-.01,.01),(.01,2.7),(2.7,5.0),(5.0,5.3),(5.3,99)]
print('%-10s %9s %9s'%('row','H x-b','b')+''.join('  dz in %4.1f-%4.1f'%b for b in bins))
for i in range(nH):
    contrib=[np.sum((H[i]*x0)[isfc2&(dz>lo)&(dz<=hi)]) for lo,hi in bins]
    print('%-10s %9.3f %9.3f'%(labs[i],r0[i],bH[i])+''.join('%17.3f'%v for v in contrib))
np.save(d+'/x0.npy',x0)
