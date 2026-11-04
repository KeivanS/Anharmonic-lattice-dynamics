"""exact minimum (P=0, no forces) of the full exact model (model3, PAR) for run_MAT, by Newton on the
supercell (strain Voigt xx yy zz yz xz xy, X = Xi F^T + u); as bprime/mkrelaxed.py but with C on."""
import sys, os, json, numpy as np
EM = os.path.expanduser('~/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/REFERENCE/exact_model')
sys.path.insert(0, EM)
from model3 import Model, PAR
from model import read_structure
mat = sys.argv[1]; run = f'{EM}/run_{mat}'
info = json.load(open(f'{run}/setup.json'))
Xi, latx = np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy')
m = Model(info['rL'], info['Ra'], dict(info.get('par', PAR))).build(Xi, latx)
ps_ref = m.ps.copy(); N = len(Xi); n = 6 + 3 * N
VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
def eta_of(e):
    t = np.zeros((3, 3))
    for k, (a, b) in enumerate(VOIGT):
        t[a, b] = t[b, a] = e[k]
    return t
def conf(q):
    F = np.eye(3) + eta_of(q[:6])
    return Xi @ F.T + q[6:].reshape(N, 3), ps_ref @ F.T
def grad(q):
    X, ps = conf(q); m.ps = ps
    dd, s = m._s(X)
    dEdx = -m.forces(X); dEdps = 2 * m.g(s)[:, None] * dd
    G = dEdx.T @ Xi + dEdps.T @ ps_ref
    return np.concatenate([[G[a, b] if a == b else G[a, b] + G[b, a] for a, b in VOIGT], dEdx.ravel()])
def hess(q, h=1e-4):
    Hm = np.zeros((n, n))
    for k in range(n):
        e = np.zeros(n); e[k] = h
        Hm[:, k] = (-grad(q + 2 * e) + 8 * grad(q + e) - 8 * grad(q - e) + grad(q - 2 * e)) / (12 * h)
    return 0.5 * (Hm + Hm.T)
q = np.zeros(n)
g0 = grad(q)
print('reference: max|dE/de| (eV)', np.abs(g0[:6]).max(), ' max|force| (eV/A)', np.abs(g0[6:]).max())
for it in range(40):
    g = grad(q)
    print(it, 'max|grad|', np.abs(g).max())
    if np.abs(g).max() < 1e-11:
        break
    q = q - np.linalg.pinv(hess(q), rcond=1e-10) @ g
# primitive-cell index of each supercell atom, from structure.params of the run
sp = f'{run}/e1all/structure.params'
if not os.path.exists(sp): sp = f'{run}/structure.params'
conv, prim, atoms = read_structure(sp)
tau = []
for x in Xi:
    t = [k for k, (_, p) in enumerate(atoms) if np.allclose(np.round(np.linalg.solve(prim.T, x - p)), np.linalg.solve(prim.T, x - p), atol=1e-6)]
    tau.append(t[0])
tau = np.array(tau); u = q[6:].reshape(N, 3)
u0 = np.array([u[tau == t].mean(0) for t in range(len(atoms))]); spread = max(np.abs(u[tau == t] - u0[t]).max() for t in range(len(atoms)))
u0 -= u0[0]
print('strain (Voigt)', q[:6].tolist())
print('u0 - u0(1) per primitive atom (A)', u0.tolist(), ' spread over images', spread)
json.dump({'strain': q[:6].tolist(), 'u0': u0.tolist(), 'maxgrad': float(np.abs(grad(q)).max())}, open(f'{sys.argv[2]}', 'w'), indent=1)
