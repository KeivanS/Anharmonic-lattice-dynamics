"""ref.py a rL NX NZ [n1,n2,n3]: exact ballistic channels per A^2 (normal n, default z) of the chain model E = sum_pairs L s,
L = K0 (1 - s0/rL^2)^2 (model3.py), tetragonal a x a x 2.5, one Ar atom. Independent of FOCEX:
M/A = (2 pi/V) (1/N) sum_k delta(nu - nu_k) |dnu/dk_z|/2 per band (3 degenerate bands), histogram on a fine mesh."""
import sys, numpy as np
a, rL = float(sys.argv[1]), float(sys.argv[2]); c = 2.5; m = 39.948; K0 = 0.05; cnst = 521.1098918
R = np.array([(i*a, j*a, k*c) for i in range(-4, 5) for j in range(-4, 5) for k in range(-6, 7)])
d = np.linalg.norm(R, axis=1); sel = (d > 1e-9) & (d < rL); R = R[sel]; L = K0*(1 - d[sel]**2/rL**2)**2
print('neighbours (|R|, L):', sorted(set((round(x, 3), round(y, 6)) for x, y in zip(d[sel], L))))
nrm = np.array([float(x) for x in sys.argv[5].split(",")]) if len(sys.argv) > 5 else np.array([0, 0, 1.0]); nrm = nrm/np.linalg.norm(nrm)
def nu_v(k):
    ph = k @ R.T
    w2 = (2/m)*((1 - np.cos(ph)) @ L)
    dw2 = (2/m)*((np.sin(ph)*L) @ (R @ nrm))
    nu = cnst*np.sqrt(np.maximum(w2, 0)); return nu, cnst**2*dw2/(2*np.maximum(nu, 1e-30))
nx, nz = int(sys.argv[3]), int(sys.argv[4])
gx = (np.arange(nx) + 0.5)/nx*2*np.pi/a; gz = (np.arange(nz) + 0.5)/nz*2*np.pi/c
K = np.array(np.meshgrid(gx, gx, gz, indexing='ij')).reshape(3, -1).T
nu, v = nu_v(K); numax = nu.max(); V = a*a*c
print(f'nu_max = {numax:.6f} 1/cm')
edges = np.linspace(0, numax*1.05, 211); h, _ = np.histogram(nu, edges, weights=np.abs(v)/2)
M = 3*(2*np.pi/V)*h/len(nu)/np.diff(edges)
for e, mm in zip(0.5*(edges[1:] + edges[:-1]), M): print(f'{e:12.5f} {mm:12.6f}')
