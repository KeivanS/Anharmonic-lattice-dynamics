"""bptest2.py RUN D: as bptest.py, with an extra pair term D (s - s0)^3 on the anharmonic pairs (pair
potential, cubic in s), the 4-body term off, and Newton relaxations with the full Hessian.
The energy is then of degree 6 in the positions; derivatives by 5-point differences, step-size checked."""
import sys, json, numpy as np
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..'))
from model3 import Model, PAR

run, D = sys.argv[1], float(sys.argv[2])
info = json.load(open(f'{run}/setup.json'))
Xi, latx = np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy')
par = dict(PAR); par['C'] = 0.0
m = Model(info['rL'], info['Ra'], par).build(Xi, latx)
s_ref = m._s(Xi)[1]; anh = s_ref < info['Ra'] ** 2
base_g, base_E = m.g, m.energy
def g3(s):
    out = base_g(s); out[anh] += 3 * D * (s[anh] - s_ref[anh]) ** 2
    return out
def E3(X):
    s = m._s(X)[1]
    return base_E(X) + D * np.sum((s[anh] - s_ref[anh]) ** 3)
m.g, m.energy = g3, E3                                     # forces() uses m.g
ps_ref = m.ps.copy(); N = len(Xi); V_ref = abs(np.linalg.det(latx)); n = 6 + 3 * N
VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
GPA = 160.21766


def eta_of(e):
    t = np.zeros((3, 3))
    for k, (a, b) in enumerate(VOIGT):
        t[a, b] = t[b, a] = e[k]
    return t


class State:
    def __init__(self, X0, ps0, V0):
        self.X0, self.ps0, self.V0 = X0, ps0, V0
    def conf(self, q):
        F = np.eye(3) + eta_of(q[:6])
        return self.X0 @ F.T + q[6:].reshape(N, 3), self.ps0 @ F.T
    def energy(self, q):
        X, ps = self.conf(q); m.ps = ps
        return m.energy(X)
    def grad(self, q):
        X, ps = self.conf(q); m.ps = ps
        d, s = m._s(X)
        dEdx = -m.forces(X); dEdps = 2 * m.g(s)[:, None] * d
        G = dEdx.T @ self.X0 + dEdps.T @ self.ps0
        return np.concatenate([[G[a, b] if a == b else G[a, b] + G[b, a] for a, b in VOIGT], dEdx.ravel()])
    def volume(self, q):
        return self.V0 * np.linalg.det(np.eye(3) + eta_of(q[:6]))
    def dvol(self, q):
        F = np.eye(3) + eta_of(q[:6]); Dm = np.linalg.det(F) * np.linalg.inv(F).T
        return np.concatenate([[Dm[a, b] if a == b else Dm[a, b] + Dm[b, a] for a, b in VOIGT], np.zeros(3 * N)]) * self.V0
    def hess(self, q, p=0.0, h=1e-4):
        """5-point differences of the analytic gradient of E + p V"""
        gf = lambda x: self.grad(x) + p * self.dvol(x)
        Hm = np.zeros((n, n))
        for k in range(n):
            e = np.zeros(n); e[k] = h
            Hm[:, k] = (-gf(q + 2 * e) + 8 * gf(q + e) - 8 * gf(q - e) + gf(q - 2 * e)) / (12 * h)
        return 0.5 * (Hm + Hm.T)


def newton(st, p, q):
    for it in range(30):
        gq = st.grad(q) + p * st.dvol(q)
        if np.abs(gq).max() < 1e-11:
            break
        q = q - np.linalg.pinv(st.hess(q, p), rcond=1e-10) @ gq
    return q, np.abs(st.grad(q) + p * st.dvol(q)).max()


st0 = State(Xi, ps_ref, V_ref)
q, gm = newton(st0, 0.0, np.zeros(n))
Xeq, pseq = st0.conf(q); Veq = st0.volume(q)
st = State(Xeq, pseq, Veq); z = np.zeros(n)
print(f'{run.split("/")[-1]} D={D}: equilibrium max|grad| {gm:.1e}, V/V_ref {Veq / V_ref:.6f}, strain {np.round(q[:6], 5).tolist()}')

# truth
for dp in (0.004, 0.002, 0.001):
    pp = np.array([-2, -1, 0, 1, 2]) * dp; lv = []
    for p in pp:
        qp, g_ = newton(st, p, z)
        lv.append(np.log(st.volume(qp) / Veq))
    c = np.polyfit(pp, lv, 4)[::-1]
    B = -1 / c[1]
    print(f'  truth  dp={dp * GPA:.3f} GPa: B = {B * GPA:.6f} GPa, B\' = {2 * c[2] * B ** 2:.6f}  (max|grad| {g_:.0e})')

# formula
H = st.hess(z)
Hee, Hue, Huu = H[:6, :6], H[6:, :6], H[6:, 6:]
Ga = np.linalg.pinv(Huu, rcond=1e-10)
Hrel = Hee - Hue.T @ Ga @ Hue
v = Veq * np.array([1, 1, 1, 0, 0, 0.])
eta1 = -np.linalg.solve(Hrel, v); u1 = -Ga @ Hue @ eta1; p1 = np.concatenate([eta1, u1])
B_f = Veq / (v @ np.linalg.solve(Hrel, v))
T1 = eta_of(eta1); V2 = Veq * (np.trace(T1) ** 2 - np.trace(T1 @ T1))
def W3(path, strain):
    lam = strain / np.abs(eta1).max()
    E = [st.energy(k * lam * path) for k in (-2, -1, 1, 2)]
    return (E[3] - 2 * E[2] + 2 * E[1] - E[0]) / (2 * lam ** 3)
for strain in (2e-3, 1e-3, 5e-4):
    w = W3(p1, strain)
    print(f'  formula (strain {strain:.0e}): B = {B_f * GPA:.6f} GPa, B\' = {B_f ** 2 / Veq * (w + 3 * V2) - 1:.6f}')
w0 = W3(np.concatenate([eta1, 0 * u1]), 1e-3)
print(f'  formula with u_1 = 0 in W3: B\' = {B_f ** 2 / Veq * (w0 + 3 * V2) - 1:.6f};  |u_1| max {np.abs(u1).max():.2e} A per eV/A^3')
