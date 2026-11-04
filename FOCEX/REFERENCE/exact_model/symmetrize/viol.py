import sys,glob,re
d=sys.argv[1]; L=open(glob.glob(d+'/log*.dat')[0]).read().split('\n')
def block(key):
    i=[k for k,l in enumerate(L) if key in l][-1]; out=[]
    for l in L[i+1:]:
        if l.startswith('*******') or not l.strip(): break
        w=l.split()
        if 'if not=15' in l: continue
        if not w[0].isdigit(): break
        out.append([float(v) for v in w[1:]])
    return out
t=block('Violation of translational'); r=block('Violation of rotational'); h=block('Violation of Huang')
vt=max(abs(x[0]) for x in t if x); vr=max(abs(x[0]-x[1]) for x in r if len(x)>1); vh=max(abs(x[0]) for x in h if x) if h else float('nan')
m=open(d+'/mech.dat').read().split('\n')
hk=[l for l in m if l.startswith('Huang: max')]; hk=hk[0].split()[-2] if hk else '-'
n=[k for k,l in enumerate(m) if 'ELASTIC_LONGWAVE: CRYSTAL FRAME C_ij' in l][0]
c=[[float(v) for v in m[n+1+i].split()] for i in range(6)]
fit=[l for l in L if 'Relative residual' in l or 'Percent error, ||' in l]
fit=fit[-1].split()[-2] if fit else '-'
print('%-7s fit %6s | max viol: transl %8.1e rot %8.1e Huang %8.1e eV  mech Huang-stress %8s GPa | C11 %5.1f C12 %5.1f C13 %5.1f C33 %5.1f C44 %5.1f C66 %5.1f'%(d.split('/')[-1],fit,vt,vr,vh,hk,c[0][0],c[0][1],c[0][2],c[2][2],c[3][3],c[5][5]))
