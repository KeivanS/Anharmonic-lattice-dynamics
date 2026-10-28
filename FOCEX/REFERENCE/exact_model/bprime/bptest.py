"""bptest.py RUN: test the B0' formula of FOCEX/notes/generalized_huang.tex (Eq. Bp) on the exact model
(pair + 3-body, 4-body coefficient C = 0), at the model's zero-pressure equilibrium.

Truth:   relax cell (6 symmetric strains) and atoms at pressures p, V(p) -> B = -V dp/dV, B' = dB/dp.
Formula: C from the exact Hessian of E(eta, u) (5-point differences of the analytic gradient, exact for a
         quartic), relaxed by the Schur complement; eta_1 = -Ht^-1 v; p_1 = (eta_1, u_1);
         W3(p1,p1,p1) = d^3/dlam^3 E(lam p1) (5-point stencil, exact for a quartic);
         B' = B^2/V [W3 + 3 V2(eta1,eta1)] - 1.
The energy is exact under homogeneous deformation because model3 depends only on the pair vectors
d_p = X_j - X_i + ps_p, and the image vectors ps are deformed with the cell."""
import sys, json, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..'))
from model3 import Model, PAR

run = sys.argv[1]
info = json.load(open(f'{run}/setup.json'))
Xi, latx = np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy')
par = dict(PAR); par['C'] = 0.0                       # pair + 3-body only
m = Model(info['rL'], info['Ra'], par).build(Xi, latx)
ps_ref = m.ps.copy(); N = len(Xi); V_ref = abs(np.linalg.det(latx))
VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def eta_of(e):
    """symmetric 3x3 from 6 components; off-diagonal eta_ab = eta_ba = e (no factor 1/2)"""
    t = np.zeros((3, 3))
    for k, (a, b) in enumerate(VOIGT):
        t[a, b] = t[b, a] = e[k]
    return t


class State:
    """deformation of a base configuration (X0, ps0, V0): x = (1+eta) X0 + u, ps = (1+eta) ps0"""
    def __init__(self, X0, ps0, V0):
        self.X0, self.ps0, self.V0 = X0, ps0, V0

    def conf(self, q):
        F = np.eye(3) + eta_of(q[:6])
        return self.X0 @ F.T + q[6:].reshape(N, 3), self.ps0 @ F.T, F

    def energy(self, q):
        X, ps, _ = self.conf(q); m.ps = ps
        return m.energy(X)

    def grad(self, q):
        X, ps, F = self.conf(q); m.ps = ps
        d, s = m._s(X); g = m.g(s)
        dEdx = -m.forces(X)                                   # dE/dx_i
        dEdps = 2 * g[:, None] * d                            # dE/dps_p
        G = dEdx.T @ self.X0 + dEdps.T @ self.ps0             # dE/dF_ab (3x3)
        ge = np.array([G[a, b] if a == b else G[a, b] + G[b, a] for a, b in VOIGT])
        return np.concatenate([ge, dEdx.ravel()])

    def volume(self, q):
        return self.V0 * np.linalg.det(np.eye(3) + eta_of(q[:6]))

    def dvol(self, q):
        F = np.eye(3) + eta_of(q[:6]); D = np.linalg.det(F) * np.linalg.inv(F).T
        return np.concatenate([[D[a, b] if a == b else D[a, b] + D[b, a] for a, b in VOIGT], np.zeros(3 * N)]) * self.V0


def relax(st, p, q0):
    f = lambda q: st.energy(q) + p * st.volume(q)
    gfun = lambda q: st.grad(q) + p * st.dvol(q)
    r = minimize(f, q0, jac=gfun, method='BFGS', options=dict(gtol=1e-11, maxiter=20000))
    for _ in range(3):                                     # polish
        r = minimize(f, r.x, jac=gfun, method='BFGS', options=dict(gtol=1e-12, maxiter=20000))
    return r.x, np.abs(gfun(r.x)).max()


# 1. zero-pressure equilibrium of the model, then make it the base
st0 = State(Xi, ps_ref, V_ref)
q, gmax = relax(st0, 0.0, np.zeros(6 + 3 * N))
Xeq, pseq, _ = st0.conf(q); Veq = st0.volume(q)
print(f'{run}: zero-pressure equilibrium: max|grad| {gmax:.1e}, V/V_ref {Veq / V_ref:.6f}, strain {np.round(q[:6], 5)}')
st = State(Xeq, pseq, Veq); n = 6 + 3 * N; z = np.zeros(n)
m.ps = pseq
print(f'  max |F| at equilibrium {np.abs(m.forces(Xeq)).max():.1e} eV/A; stress (GPa) {np.round(m.stress(Xeq, Veq) * 160.21766, 9).tolist()}')

# 2. truth: V(p) by relaxation at small pressures (eV/A^3), B and B' from ln V(p)
def truth(dp):
    ps_ = np.array([-2, -1, 0, 1, 2]) * dp
    lv = []
    for p in ps_:
        qp, gm = relax(st, p, z)
        lv.append(np.log(st.volume(qp) / Veq))
    c = np.polyfit(ps_, lv, 4)[::-1]                       # ln V = c0 + c1 p + c2 p^2 + ...
    B = -1 / c[1]; Bp = 2 * c[2] * B ** 2
    return B, Bp
for dp in (0.004, 0.002):
    B, Bp = truth(dp)
    print(f'  truth (dp = {dp * 160.21766:.3f} GPa): B = {B * 160.21766:.6f} GPa, B\' = {Bp:.6f}')

# 3. formula: exact Hessian (5-point differences of the analytic gradient) in the strains
h = 1e-3
Hcol = []
for k in range(6):
    e = np.zeros(n); e[k] = h
    Hcol.append((-st.grad(2 * e) + 8 * st.grad(e) - 8 * st.grad(-e) + st.grad(-2 * e)) / (12 * h))
Hcol = np.array(Hcol).T                                    # n x 6: d grad / d eta_k
Hee, Hue = Hcol[:6, :], Hcol[6:, :]
m.ps = pseq; Huu = m.hessian(Xeq)
Ga = np.linalg.pinv(Huu, rcond=1e-10)                      # removes the 3 translations
Hrel = Hee - Hue.T @ Ga @ Hue                              # V0 C in the 6 strain coordinates
v = Veq * np.array([1, 1, 1, 0, 0, 0.])
eta1 = -np.linalg.solve(Hrel, v); u1 = -Ga @ Hue @ eta1
p1 = np.concatenate([eta1, u1])
B_f = Veq / (v @ np.linalg.solve(Hrel, v))
lam = 0.2 / np.abs(eta1).max() * 1e-2                      # strains of 0.2 %
E = [st.energy(k * lam * p1) for k in (-2, -1, 1, 2)]
W3 = (E[3] - 2 * E[2] + 2 * E[1] - E[0]) / (2 * lam ** 3)
T1 = eta_of(eta1); V2 = Veq * (np.trace(T1) ** 2 - np.trace(T1 @ T1))
Bp_f = B_f ** 2 / Veq * (W3 + 3 * V2) - 1
print(f'  formula: B = {B_f * 160.21766:.6f} GPa, B\' = {Bp_f:.6f}   (|u_1| max {np.abs(u1).max():.3e} A per eV/A^3)')
# without the internal relaxation in W3 (u_1 = 0 along the path), and the paper's/cubic shortcuts
E0 = [st.energy(k * lam * np.concatenate([eta1, 0 * u1])) for k in (-2, -1, 1, 2)]
W3n = (E0[3] - 2 * E0[2] + 2 * E0[1] - E0[0]) / (2 * lam ** 3)
print(f'  same with u_1 = 0 in W3: B\' = {B_f ** 2 / Veq * (W3n + 3 * V2) - 1:.6f}')
iso = np.array([1, 1, 1, 0, 0, 0.]); li = 2e-3
Ei = [st.energy(np.concatenate([k * li * iso, np.zeros(3 * N)])) for k in (-2, -1, 1, 2)]
f3 = (Ei[3] - 2 * Ei[2] + 2 * Ei[1] - Ei[0]) / (2 * li ** 3)
print(f'  isotropic strain, no relaxation: f\'\'\' = {f3:.6e} eV; 1 - f\'\'\'/(27VB) = {1 - f3 / (27 * Veq * B_f):.6f}; '
      f'paper f\'\'\'/(6*27VB) = {f3 / (6 * 27 * Veq * B_f):.6f}; code -f\'\'\'/(27VB) = {-f3 / (27 * Veq * B_f):.6f}')
