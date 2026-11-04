import numpy as np, sys
exec(open(sys.argv[2]).read().split("AK=af@K")[0])   # reuse loading/labels from huang_diag.py
print('Huang row norms (max |coef|) and RHS:')
for i in range(nH): print('%-9s %12.3e %12.3e'%(labs[i],np.abs(H[i]).max(),bH[i]))
for name,CC in [('transl',ahom[:nt]),('transl+Huang',np.vstack([ahom[:nt],H])),('transl+rot+Huang',ahom)]:
    s=np.linalg.svd(CC,compute_uv=False); 
    print(name,'rank(1e-10)',(s>1e-10*s.max()).sum(),' smallest kept s/smax:',', '.join('%.2e'%v for v in sorted(s[s>1e-10*s.max()]/s.max())[:4]),
          ' largest dropped:', '%.2e'%(s[s<=1e-10*s.max()].max()/s.max()) if (s<=1e-10*s.max()).any() else '-')
