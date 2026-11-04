import numpy as np, sys
exec(open(sys.argv[2]).read().split("AK=af@K")[0])
cut=1e-10
res=lambda x: np.linalg.norm(af@x-bf)/np.linalg.norm(bf)*100
def pinv_solve(C,d):
    u,s,vt=np.linalg.svd(C,full_matrices=False); m=s>cut*s.max()
    return vt[m].T@((u[:,m].T@d)/s[m])
def cfit(C,d):
    xp=pinv_solve(C,d)
    u,s,vt=np.linalg.svd(C); r=(s>cut*s.max()).sum(); K=vt[r:].T
    AK=af@K; y=np.linalg.lstsq(AK,bf-af@xp,rcond=cut)[0]
    return xp+K@y, xp
T=ahom[:nt]; dT=bhom[:nt]
cases={'transl only':(T,dT),
 'transl+all 15 Huang (code)':(np.vstack([T,H]),np.r_[dT,bH]),
 'transl+all Huang, RHS of forbidden rows=0':(np.vstack([T,H]),np.r_[dT,np.where(np.arange(nH)>=13,bH,0)]),
 'transl+[zz,xx] only':(np.vstack([T,H[13:14]]),np.r_[dT,bH[13:14]])}
xref=None
for n,(C,d) in cases.items():
    x,xp=cfit(C,d)
    if xref is None: xref=x
    print('%-45s fit %7.3f %%  |xp|max %9.3e  max|dFC2| vs transl-only %9.3e  [zz,xx] viol %8.4f'%(n,res(x),np.abs(xp).max(),np.abs((x-xref)[isfc2]).max(),(H[13]@x-bH[13])))
