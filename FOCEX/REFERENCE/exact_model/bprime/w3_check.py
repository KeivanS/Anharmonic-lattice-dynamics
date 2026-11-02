"""w3_check.py: W3(p1,p1,p1) of the relaxed ZnO model (D = 0, C = 0) three ways:
 truth  = third derivative of the supercell energy along p1 = (eta1, u1) (5-point, exact for degree 4),
 A-form = -1/2 sum_{tau,j} dPhi^{ab}_{tau j} (S_j-S_tau)_a (S_j-S_tau)_b, dPhi = sum_k Psi_{tau j k}(S_k-S_tau)
          (what FOCEX's pressure_derivative_b computed), from the Hessian's derivative along p1,
 BH-form = the strain-strain part through the Born-Huang brackets of dPhi (+ the stress of the deformed
          crystal), plus the internal parts, i.e. d/dlam of p1^T H(lam p1) p1 with H's strain block = V c."""
import sys, os, json
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.argv = ['mkrelaxed.py', 'ZnO', '__none__']
src = open(os.path.join(here, 'mkrelaxed.py')).read().replace("if os.path.exists(d):\n    sys.exit(f'{d} exists')", "").split("os.mkdir(d)")[0]
exec(src)
GPA = 160.21766208
X0, P0, lat = Xeq.copy(), pseq.copy(), lat_eq.copy()
Vs = abs(np.linalg.det(lat)); ncell = N // 4


def conf(qq):
    Fq = np.eye(3) + eta_of(qq[:6])
    return X0 @ Fq.T + qq[6:].reshape(N, 3), P0 @ Fq.T


def energy(qq):
    X, ps = conf(qq); m.ps = ps
    return m.energy(X)


def grad(qq):
    X, ps = conf(qq); m.ps = ps
    dd, s = m._s(X)
    dEdx = -m.forces(X); dEdps = 2 * m.g(s)[:, None] * dd
    G = dEdx.T @ X0 + dEdps.T @ P0
    return np.concatenate([[G[a, b] if a == b else G[a, b] + G[b, a] for a, b in VOIGT], dEdx.ravel()])


Hm = np.zeros((n, n)); h = 1e-4
for k in range(n):
    e = np.zeros(n); e[k] = h
    Hm[:, k] = (-grad(2 * e) + 8 * grad(e) - 8 * grad(-e) + grad(-2 * e)) / (12 * h)
Hm = (Hm + Hm.T) / 2
Hee, Hue, Huu = Hm[:6, :6], Hm[6:, :6], Hm[6:, 6:]
Ga = np.linalg.pinv(Huu, rcond=1e-10)
Hrel = Hee - Hue.T @ Ga @ Hue
v = Vs * np.array([1, 1, 1, 0, 0, 0.])
eta1 = -np.linalg.solve(Hrel, v); u1 = -Ga @ Hue @ eta1; p1 = np.concatenate([eta1, u1])
B = Vs / (v @ np.linalg.solve(Hrel, v))
T1 = eta_of(eta1); V2 = Vs * (np.trace(T1) ** 2 - np.trace(T1 @ T1))
lam = 1e-3 / np.abs(eta1).max()
E = [energy(k * lam * p1) for k in (-2, -1, 1, 2)]
W3true = (E[3] - 2 * E[2] + 2 * E[1] - E[0]) / (2 * lam ** 3)
print(f'B = {B * GPA:.6f} GPa; B\' from true W3 = {B ** 2 / Vs * (W3true + 3 * V2) - 1:.6f} (reference -0.157427)')

# A-form with the exact dPhi (derivative of the supercell Hessian along p1)
def hess_at(t):
    X, ps = conf(t * p1); m.ps = ps
    return m.hessian(X).reshape(N, 3, N, 3)
dl = 1e-4 / np.abs(eta1).max()
dPhi = (hess_at(dl) - hess_at(-dl)) / (2 * dl)
S = (X0 @ T1.T) + u1.reshape(N, 3)          # displacement field per unit lambda (u1 per supercell atom)
il = np.linalg.inv(lat)
W3A = 0
for i in range(N):
    for j in range(N):
        dv = X0[j] - X0[i]; f = dv @ il; img = -np.rint(f) @ lat      # minimum image of j
        dS = (S[j] + img @ T1.T) - S[i]
        W3A += -0.5 * dS @ dPhi[i, :, j, :] @ dS
print(f'W3 true {W3true:.6e}, A-form {W3A:.6e} (eV per supercell, lambda = unit pressure units)')
print(f'B\' from A-form W3 = {B ** 2 / Vs * (W3A + 3 * V2) - 1:.6f}')

# --- B' from B(+-lam) of the crystal deformed along p1, elastic constants by Born-Huang (elastic_born_huang) ---
conv_, prim_, atoms_ = read_structure(f'{run}/structure.params'); ip_ = np.linalg.inv(prim_)
cell0 = []; sub = np.zeros(N, int)
for nn in range(N):
    for t, (sp, p) in enumerate(atoms_):
        fr = (Xi[nn] - p) @ ip_
        if np.abs(fr - np.rint(fr)).max() < 1e-6:
            sub[nn] = t
for t in range(len(atoms_)):
    cell0.append(int(np.where(sub == t)[0][0]))
nat = len(atoms_)
Phi0 = m.hessian(X0).reshape(N, 3, N, 3) if False else None
m.ps = P0; H0 = m.hessian(X0).reshape(N, 3, N, 3)
dvec = {}
for i in cell0:
    for j in range(N):
        dv = X0[j] - X0[i]; f = dv @ il; dvec[(i, j)] = (f - np.rint(f)) @ lat
U1 = u1.reshape(N, 3)


def c_relaxed(Phi, shift, sigma, V):
    """relaxed Lagrangian elastic constants (eV/A^3, 3x3x3x3) of the deformed crystal"""
    A = np.zeros((3, 3, 3, 3)); Q = np.zeros((nat, 3, 3, 3)); ph0 = np.zeros((nat, 3, nat, 3))
    for t, i in enumerate(cell0):
        for j in range(N):
            dd = dvec[(i, j)] + shift(i, j)
            A += -0.5 * np.einsum('ag,b,d->abgd', Phi[i, :, j, :], dd, dd)
            Q[t] += np.einsum('ag,d->agd', Phi[i, :, j, :], dd)
            ph0[t, :, sub[j], :] += Phi[i, :, j, :]
    A /= V; dl = np.eye(3)
    C = A + A.transpose(1, 0, 2, 3) - A.transpose(0, 2, 1, 3) \
        + 0.5 * (np.einsum('ab,gd->abgd', dl, sigma) + np.einsum('gd,ab->abgd', dl, sigma)) \
        - np.einsum('ag,bd->abgd', dl, sigma) - np.einsum('bd,ag->abgd', dl, sigma)
    C = 0.25 * (C + C.transpose(1, 0, 2, 3) + C.transpose(0, 1, 3, 2) + C.transpose(1, 0, 3, 2))
    Qs = 0.5 * (Q + Q.transpose(0, 1, 3, 2)).reshape(3 * nat, 3, 3)
    G = np.zeros((3 * nat, 3 * nat)); G[3:, 3:] = np.linalg.inv(ph0.reshape(3 * nat, 3 * nat)[3:, 3:])
    QGQ = np.einsum('sab,st,tgd->abgd', Qs, G, Qs) / V
    return C - QGQ


def B_of(C, P, V):
    """B = V/(v Ht^-1 v), Ht = V (c + S2) + P V V2: the Hessian in the stretch eta (tensor components) of a
    crystal under hydrostatic pressure P; S2 = d^2(sigma:eta^2/2), sigma = -P"""
    Mc = np.array([[C[a, b, g, d] for (g, d) in VOIGT] for (a, b) in VOIGT]) * np.outer([1, 1, 1, 2, 2, 2.], [1, 1, 1, 2, 2, 2.])
    for k in range(6):
        for l in range(6):
            Ek, El = eta_of(np.eye(6)[k]), eta_of(np.eye(6)[l])
            Mc[k, l] += np.sum(-P * np.eye(3) * (Ek @ El + El @ Ek)) / 2
    Mv = np.zeros((6, 6)); Mv[0, 1] = Mv[1, 0] = Mv[0, 2] = Mv[2, 0] = Mv[1, 2] = Mv[2, 1] = 1; Mv[3, 3] = Mv[4, 4] = Mv[5, 5] = -2
    Ht = V * Mc + P * V * Mv; vv = V * np.array([1, 1, 1, 0, 0, 0.])
    return V / (vv @ np.linalg.solve(Ht, vv))


Vc = Vs / ncell
Bs = []
for sgn in (1, -1):
    lm = sgn * 1e-3 / np.abs(eta1).max() * 0 + sgn * 1e-4      # pressure step (eV/A^3)
    Phi_l = H0 + lm * dPhi
    shift = lambda i, j, lm=lm: lm * (T1 @ dvec[(i, j)] + U1[j] - U1[i])
    V_l = Vc * np.linalg.det(np.eye(3) + lm * T1)
    C_l = c_relaxed(Phi_l, shift, -lm * np.eye(3), V_l)
    Bs.append(B_of(C_l, lm, V_l))
C0 = c_relaxed(H0, lambda i, j: 0, np.zeros((3, 3)), Vc)
print(f'B(0) from c_relaxed = {B_of(C0, 0, Vc) * GPA:.6f} GPa;  B\' = (B(+l)-B(-l))/(2 l) = {(Bs[0] - Bs[1]) / 2e-4:.6f}  (reference -0.157427)')
