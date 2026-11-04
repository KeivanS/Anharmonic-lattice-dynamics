import numpy as np, sys
# shift each supercell atom by delta of its primitive atom (z only, wurtzite GaN, reduced units of c)
c=5.2724287802
zo=np.array([0.5005917000,0.0005865644,0.8777501978,0.3777384344]); zs=np.array([0.5005891322,0.0005891322,0.8777443161,0.3777443161])
inp,out=sys.argv[1],sys.argv[2]
l=open(inp).read().split('\n')
L=np.array([[float(v) for v in l[i].split()[:3]] for i in (2,3,4)])*float(l[1])
n=[int(v) for v in l[6].split()]; N=sum(n); mode=l[7].strip()[0].lower()
types=sum([[t]*k for t,k in enumerate(n)],[])
pos=np.array([[float(v) for v in l[8+i].split()[:3]] for i in range(N)])
cart=pos@L if mode=='d' else pos.copy()
mx=0
for i in range(N):
    z=(cart[i,2]/c)%1.0
    cand=[0,1] if types[i]==0 else [2,3]
    d=[abs(((z-zo[k]+0.5)%1)-0.5) for k in cand]; k=cand[int(np.argmin(d))]; mx=max(mx,min(d)*c)
    cart[i,2]+= (zs[k]-zo[k])*c
print(inp,'atoms',N,'max match error (A) %.2e'%mx)
newpos=cart@np.linalg.inv(L) if mode=='d' else cart
for i in range(N): l[8+i]='  %.12f %.12f %.12f'%tuple(newpos[i])
open(out,'w').write('\n'.join(l))
