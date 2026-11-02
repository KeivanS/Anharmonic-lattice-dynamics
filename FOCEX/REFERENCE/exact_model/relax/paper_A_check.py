"""paper_A_check.py MAT STATE (STATE = ref, or a model_min.py json for the exact minimum): the exact strain
Hessian of the model at a force-free state (internal coordinates relaxed; finite differences of dE/de),
unrelaxed and relaxed, against (1) FOCEX's tangent V[c_BH + stress - QGQ/V + S2] (relax_reference) and
(2) the ALATDYN paper's Sym(A) - QGQ/V with A = -1/2 sum K d d / V. Also a direct check of
E2 = -1/4 sum (eps d) Phi (eps d) against energies (the paper's step; see cluster_identity.py)."""
import sys, os, json, numpy as np
EM = os.path.expanduser('~/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/REFERENCE/exact_model')
sys.path.insert(0, EM)
from model3 import Model, PAR
from model import read_structure
mat = sys.argv[1]; run = f'{EM}/run_{mat}'
info = json.load(open(f'{run}/setup.json'))
Xi, latx = np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy')
m = Model(info['rL'], info['Ra'], dict(PAR)).build(Xi, latx)
ps_ref = m.ps.copy(); N = len(Xi)
VO = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
def eta_of(e):
    t = np.zeros((3, 3))
    for k, (a, b) in enumerate(VO):
        t[a, b] = t[b, a] = e[k]
    return t
def grad_at(X0, P0, e, u):
    F = np.eye(3) + eta_of(e)
    X = X0 @ F.T + u; m.ps = P0 @ F.T
    dd, s = m._s(X)
    dEdx = -m.forces(X); dEdps = 2 * m.g(s)[:, None] * dd
    G = dEdx.T @ X0 + dEdps.T @ P0
    return np.array([G[a, b] if a == b else G[a, b] + G[b, a] for a, b in VO]), dEdx
def relax(X0, P0, e, u):
    for it in range(30):
        F = np.eye(3) + eta_of(e); X = X0 @ F.T + u; m.ps = P0 @ F.T
        f = m.forces(X)
        if np.abs(f).max() < 1e-12: break
        Hs = m.hessian(X)
        du = np.linalg.lstsq(Hs, f.ravel(), rcond=1e-10)[0].reshape(N, 3)
        u = u + du
    return u
# state: reference cell, internal coordinates relaxed
if sys.argv[2] == 'ref':
    Fs = np.eye(3)
else:
    Fs = np.eye(3) + eta_of(json.load(open(sys.argv[2]))['strain'])
u = relax(Xi @ Fs.T, ps_ref @ Fs.T, np.zeros(6), np.zeros((N, 3)))
X0 = Xi @ Fs.T + u; P0 = ps_ref @ Fs.T
m.ps = P0
lat = latx @ Fs.T; V = abs(np.linalg.det(lat)); conv, prim, atoms = read_structure(f'{run}/e1all/structure.params')
nprim = len(atoms); Vc = V / (N // nprim)
dd, s = m._s(X0); sig = 2 * np.einsum('p,pa,pb->ab', m.g(s), dd, dd) / V
print('state: max|force|', np.abs(m.forces(X0)).max(), ' stress (GPa) diag', np.diag(sig) * 160.21766208)
# exact relaxed Hessian by central differences of G(e) with u relaxed (positions X0 F^T + u)
h = 1e-5; Hrel = np.zeros((6, 6))
for k in range(6):
    e = np.zeros(6); e[k] = h
    up = relax(X0, P0, e, np.zeros((N, 3))); um = relax(X0, P0, -e, np.zeros((N, 3)))
    Hrel[:, k] = (grad_at(X0, P0, e, up)[0] - grad_at(X0, P0, -e, um)[0]) / (2 * h)
Hrel = (Hrel + Hrel.T) / 2 / (N // nprim)        # per primitive cell
# FOCEX formula from the exact FC2 at the state
Hs = m.hessian(X0).reshape(N, 3, N, 3); il = np.linalg.inv(lat)
ip = np.linalg.inv(prim); cell0 = []
for t, (sp, p) in enumerate(atoms):
    for nn in range(N):
        fr = (Xi[nn] - p) @ ip
        if np.abs(fr - np.rint(fr)).max() < 1e-6: cell0.append(nn); break
tau_of = np.zeros(N, int)
for nn in range(N):
    for t, (sp, p) in enumerate(atoms):
        fr = (Xi[nn] - p) @ ip
        if np.abs(fr - np.rint(fr)).max() < 1e-6: tau_of[nn] = t
A = np.zeros((3, 3, 3, 3)); Q = np.zeros((nprim, 3, 3, 3)); phi = np.zeros((nprim, 3, nprim, 3))
for ti, i in enumerate(cell0):
    for j in range(N):
        dv = X0[j] - X0[i]; f = dv @ il; dv = (f - np.rint(f)) @ lat
        A += -0.5 * np.einsum('ag,b,d->abgd', Hs[i, :, j, :], dv, dv)
        Q[ti] += np.einsum('ag,b->agb', Hs[i, :, j, :], dv)       # Q(tau a; g b)
        phi[ti, :, tau_of[j], :] += Hs[i, :, j, :]
A /= Vc
CBH = A + np.einsum('bagd->abgd', A) - np.einsum('agbd->abgd', A)
dl = np.eye(3)
ST = 0.5 * (np.einsum('ab,gd->abgd', dl, sig) + np.einsum('gd,ab->abgd', dl, sig)) - np.einsum('ag,bd->abgd', dl, sig) - np.einsum('bd,ag->abgd', dl, sig)
sym = lambda T: 0.25 * (T + T.transpose(1, 0, 2, 3) + T.transpose(0, 1, 3, 2) + T.transpose(1, 0, 3, 2))
c = sym(CBH + ST)
Qs = 0.5 * (Q + Q.transpose(0, 1, 3, 2)); Qr = Qs[1:].reshape(3 * (nprim - 1), 3, 3)
G = np.linalg.inv(phi[1:, :, 1:, :].reshape(3 * (nprim - 1), 3 * (nprim - 1)))
crel = c - np.einsum('kab,kl,lcd->abcd', Qr, G, Qr) / Vc
H = np.zeros((6, 6))
for k in range(6):
    Ek = eta_of(np.eye(6)[k])
    for l in range(6):
        El = eta_of(np.eye(6)[l])
        H[k, l] = Vc * (np.einsum('ab,abcd,cd', Ek, crel, El) + np.sum(sig * (Ek @ El + El @ Ek)) / 2)
np.set_printoptions(precision=6, suppress=True, linewidth=160)
print('exact relaxed Hessian per cell (eV):'); print(Hrel)
print('FOCEX formula:'); print(H)
print('max |diff| (eV)', np.abs(H - Hrel).max(), ' in GPa', np.abs(H - Hrel).max() / Vc * 160.21766208)
Hp = (H + H.T) / 2
print('pair-symmetrized: max |diff| (eV)', np.abs(Hp - Hrel).max(), ' GPa', np.abs(Hp - Hrel).max() / Vc * 160.21766208)
print('largest diff at', np.unravel_index(np.abs(Hp - Hrel).argmax(), H.shape))
# unrelaxed: exact FD (u fixed at 0 on top of X0 F^T) vs c (no QGQ)
Hun = np.zeros((6, 6)); z = np.zeros((N, 3))
for k in range(6):
    e = np.zeros(6); e[k] = h
    Hun[:, k] = (grad_at(X0, P0, e, z)[0] - grad_at(X0, P0, -e, z)[0]) / (2 * h)
Hun = (Hun + Hun.T) / 2 / (N // nprim)
Hu = np.zeros((6, 6))
for k in range(6):
    Ek = eta_of(np.eye(6)[k])
    for l in range(6):
        El = eta_of(np.eye(6)[l])
        Hu[k, l] = Vc * (np.einsum('ab,abcd,cd', Ek, c, El) + np.sum(sig * (Ek @ El + El @ Ek)) / 2)
Hu = (Hu + Hu.T) / 2
print('unrelaxed: max |formula - exact| (eV)', np.abs(Hu - Hun).max(), np.unravel_index(np.abs(Hu - Hun).argmax(), Hu.shape))
QGQ = Hu - Hp
print('exact relaxation part (Hun - Hrel) vs formula QGQ part: max diff', np.abs((Hun - Hrel) - QGQ).max())
print('Q antisymmetric part (rotational inv. with zero forces => 0): max', np.abs(Q - Q.transpose(0, 1, 3, 2)).max())

# ---- the paper: E = 1/2 V Sym(A) ee (unrelaxed) and Sym(A) - QGQ/V (relaxed), A = -1/2 sum K d d / V
cp = sym(A); cpr = cp - np.einsum('kab,kl,lcd->abcd', Qr, G, Qr) / Vc
def hmat(cc, stressterm):
    M = np.zeros((6, 6))
    for k in range(6):
        Ek = eta_of(np.eye(6)[k])
        for l in range(6):
            El = eta_of(np.eye(6)[l])
            M[k, l] = Vc * np.einsum('ab,abcd,cd', Ek, cc, El) + (Vc * np.sum(sig * (Ek @ El + El @ Ek)) / 2 if stressterm else 0)
    return (M + M.T) / 2
for name, hu, hr in [('paper Sym(A)', hmat(cp, False), hmat(cpr, False)), ('paper Sym(A) + sigma term', hmat(cp, True), hmat(cpr, True)), ('FOCEX Born-Huang+stress', Hu, Hp)]:
    print(f'{name:28s} unrelaxed max|formula-exact| {np.abs(hu-Hun).max():.4e} eV   relaxed {np.abs(hr-Hrel).max():.4e} eV  (max|H| {np.abs(Hrel).max():.3f})')
# brackets [ab,gl] = A[a,g,b,l] (A[a,b,g,d] = -1/2 sum K^{ag} d^b d^d / Vc); Huang: [ab,gl] = [gl,ab]
Br = A.transpose(0, 2, 1, 3)
print('Huang (proper): max|[ab,gl]-[gl,ab]| (eV/A^3)', np.abs(Br - Br.transpose(2, 3, 0, 1)).max(), ' max|[ab,gl]|', np.abs(Br).max())
print('rL, Ra', info['rL'], info['Ra'], ' supercell vector lengths', np.linalg.norm(lat, axis=1))
print('exact unrelaxed Hun:'); print(Hun)
print('paper unrelaxed Vc Sym(A):'); print(hmat(cp, False))
# brute force: E(e) - E(0) along Voigt 1 and 6 (u fixed), per primitive cell
def Etot(e):
    F = np.eye(3) + eta_of(e); m.ps = P0 @ F.T
    return m.energy(X0 @ F.T) if hasattr(m, 'energy') else None
for k in (0, 5):
    hh = 1e-4; e = np.zeros(6); e[k] = hh
    if Etot(e) is not None:
        print('Voigt', k + 1, 'd2E/de2 per cell from energies', (Etot(e) - 2 * Etot(0 * e) + Etot(-e)) / hh**2 / (N // nprim))
# direct: E2 = -1/4 sum_{i in cell, j} (eps d) Phi_ij (eps d), against energies, for several eps
rng = np.random.default_rng(1)
for trial in range(3):
    ev = rng.normal(size=6) if trial else np.array([1., 1, 0, 0, 0, 0])
    eps = eta_of(ev)
    e2 = 0.
    for ti, i in enumerate(cell0):
        for j in range(N):
            dv = X0[j] - X0[i]; f = dv @ il; dv = (f - np.rint(f)) @ lat
            w = eps @ dv
            e2 += -0.25 * w @ Hs[i, :, j, :] @ w
    hh = 1e-4
    eE = (Etot(hh * ev) - 2 * Etot(0 * ev) + Etot(-hh * ev)) / hh**2 / (N // nprim) / 2
    print('eps', np.round(ev, 3), ' E2 from energies', eE, '  -1/4 sum (eps d) Phi (eps d)', e2)
# FD check of Hs on one pair block
Xp = X0.copy(); i0, j0 = cell0[1], 7; hh = 1e-5
Xp[j0, 0] += hh; fp = m.forces(Xp); Xp[j0, 0] -= 2 * hh; fm = m.forces(Xp)
m.ps = P0
print('Hs check: FD', -(fp - fm)[i0] / (2 * hh), ' Hs', Hs[i0, :, j0, 0])
