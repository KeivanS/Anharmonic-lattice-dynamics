"""model_min_c.py COULOMBRUN OUT.json: exact minimum (P=0, no forces) of the exact model + point charges
(coulomb/mkcoulomb.py truth: model3 of the base run + Ewald point charges q in scalar eps, periodic supercell).
Gradient exact: model part analytic (as model_min.py); Coulomb forces analytic (Ewald), Coulomb dE/deta by
central differences of the Ewald energy (h = 1e-5). Steps: the model's Hessian (finite differences of its
analytic gradient), iterated until the total gradient is below 1e-9 (the Coulomb Hessian is left out of the
step, which only slows the convergence). Writes strain (Voigt) and u0 - u0(1) per primitive atom as model_min.py."""
import sys, os, json, numpy as np
EM = os.path.expanduser('~/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/REFERENCE/exact_model')
sys.path.insert(0, EM)
from model3 import Model, PAR
from model import read_structure
crun, outj = sys.argv[1], sys.argv[2]
cj = json.load(open(f'{crun}/coulomb.json')); base = cj['base']
info = json.load(open(f'{base}/setup.json'))
Xi, latx = np.load(f'{base}/Xi.npy'), np.load(f'{base}/latx.npy')
# Ewald of mkcoulomb.py (same charges, eps, lambda)
sys.argv = ['mkcoulomb.py', base, '/dev/null', ','.join(str(v) for v in cj['q']), str(cj['eps'])]
exec(open(f'{EM}/coulomb/mkcoulomb.py').read().split('lam = 0.30')[0])
lam = cj['lam']
m = Model(info['rL'], info['Ra'], dict(info.get('par', PAR))).build(Xi, latx)
ps_ref = m.ps.copy(); N = len(Xi); n = 6 + 3 * N
VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
def eta_of(e):
    t = np.zeros((3, 3))
    for k, (a, b) in enumerate(VOIGT): t[a, b] = t[b, a] = e[k]
    return t
def grad_model(q):
    F = np.eye(3) + eta_of(q[:6]); X = Xi @ F.T + q[6:].reshape(N, 3); m.ps = ps_ref @ F.T
    dd, s = m._s(X); dEdx = -m.forces(X); dEdps = 2 * m.g(s)[:, None] * dd
    G = dEdx.T @ Xi + dEdps.T @ ps_ref
    return np.concatenate([[G[a, b] if a == b else G[a, b] + G[b, a] for a, b in VOIGT], dEdx.ravel()])
def grad_coul(q, h=1e-5):
    F = np.eye(3) + eta_of(q[:6]); U = q[6:].reshape(N, 3)
    X = Xi @ F.T + U; Lat = latx @ F.T
    g = np.zeros(n); g[6:] = -ewald(X, Lat, lam)[1].ravel()
    for k in range(6):            # dE/de_k at fixed u: X = Xi F^T + u with F = 1 + eta(e)
        e = np.zeros(6); e[k] = h
        Fp = F + eta_of(e); Fm = F - eta_of(e)
        g[k] = (ewald(Xi @ Fp.T + U, latx @ Fp.T, lam)[0] - ewald(Xi @ Fm.T + U, latx @ Fm.T, lam)[0]) / (2 * h)
    return g
def hess_model(q, h=1e-4):
    H = np.zeros((n, n))
    for k in range(n):
        e = np.zeros(n); e[k] = h
        H[:, k] = (-grad_model(q + 2 * e) + 8 * grad_model(q + e) - 8 * grad_model(q - e) + grad_model(q - 2 * e)) / (12 * h)
    return 0.5 * (H + H.T)
q = np.zeros(n)
Hi = None
for it in range(200):
    g = grad_model(q) + grad_coul(q)
    print(it, 'max|grad|', np.abs(g).max(), 'strain', np.round(q[:6], 7).tolist(), flush=True)
    if np.abs(g).max() < 1e-9:
        break
    if Hi is None or it % 5 == 0:
        Hi = np.linalg.pinv(hess_model(q), rcond=1e-10)
    q = q - Hi @ g
sp = f'{base}/structure.params'
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
json.dump({'strain': q[:6].tolist(), 'u0': u0.tolist(), 'maxgrad': float(np.abs(g).max())}, open(outj, 'w'), indent=1)
