"""elastic_check.py MAT [relaxed]: FOCEX's elastic-constant formula against the exact model.
A_{ab,gd} = -1/(2V) sum_{tau in cell, j} Phi^{ag}_{tau j} d_b d_d (a,g displacements; what FOCEX computes,
atld_raw), from the model's supercell Hessian (minimum-image d; the FC2 range must be inside the supercell).
Exact: H = d^2E/de de for F = 1 + eta(e), eta symmetric (6 tensor components), at the state (reference, or
relaxed with D = 0, C = 0 as bptest2.py), unrelaxed (no internal relaxation). With E = E0 + V s:eL + V/2 c eL eL,
eL = eta + eta^2/2: H/V = c + [d^2 (s:eta^2/2) / de de]; so c = H/V - that term.
Compares, in Voigt (GPa): Sym(A) (FOCEX until 2026-11-01), Born-Huang C_BH = [ag,bd]+[bg,ad]-[ab,gd] with
[ab,gd] = A_{ag,bd}, C_BH - Sym(d_ag s_bd) (FOCEX's stress term), and
C_BH + Sym[ (d_ab s_gd + d_gd s_ab)/2 - d_ag s_bd - d_bd s_ag ], against c."""
import sys, os, json
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
from model3 import Model, PAR
from model import read_structure

mat = sys.argv[1]; relaxed = len(sys.argv) > 2 and sys.argv[2] == 'relaxed'
run = os.path.join(here, '..', f'run_{mat}')
info = json.load(open(f'{run}/setup.json'))
Xi, latx = np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy')
par = dict(PAR)
if relaxed:
    par['C'] = 0.0
m = Model(info['rL'], info['Ra'], par).build(Xi, latx)
ps_ref = m.ps.copy(); N = len(Xi); n = 6 + 3 * N
VO = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
GPA = 160.21766208


def eta_of(e):
    t = np.zeros((3, 3))
    for k, (a, b) in enumerate(VO):
        t[a, b] = t[b, a] = e[k]
    return t


def make(X0, P0):
    def conf(q):
        F = np.eye(3) + eta_of(q[:6])
        return X0 @ F.T + q[6:].reshape(N, 3), P0 @ F.T

    def grad(q):
        X, ps = conf(q); m.ps = ps
        dd, s = m._s(X)
        dEdx = -m.forces(X); dEdps = 2 * m.g(s)[:, None] * dd
        G = dEdx.T @ X0 + dEdps.T @ P0
        return np.concatenate([[G[a, b] if a == b else G[a, b] + G[b, a] for a, b in VO], dEdx.ravel()])
    return conf, grad


def hess(grad, h=1e-4):
    Hm = np.zeros((n, n))
    for k in range(n):
        e = np.zeros(n); e[k] = h
        Hm[:, k] = (-grad(2 * e) + 8 * grad(e) - 8 * grad(-e) + grad(-2 * e)) / (12 * h)
    return 0.5 * (Hm + Hm.T)


X0, P0, lat = Xi.copy(), ps_ref.copy(), latx.copy()
if relaxed:
    conf, grad = make(Xi, ps_ref); q = np.zeros(n)
    for it in range(30):
        g = grad(q)
        if np.abs(g).max() < 1e-11:
            break
        q = q - np.linalg.pinv(hess(lambda x: grad(q + x)), rcond=1e-10) @ g
    F = np.eye(3) + eta_of(q[:6]); X0, P0 = conf(q); lat = latx @ F.T
conf, grad = make(X0, P0)
H = hess(grad)
V = abs(np.linalg.det(lat)); Vc = V / (N // len(read_structure(f'{run}/structure.params')[2]))
m.ps = P0
d, s = m._s(X0)
sig = 2 * np.einsum('p,pa,pb->ab', m.g(s), d, d) / V            # (1/V) dE/d eta, eV/A^3
fx = np.abs(m.forces(X0)).max()
# the sigma:eta^2/2 term in the 6 tensor components
S2 = np.zeros((6, 6))
for k in range(6):
    for l in range(6):
        Ek, El = eta_of(np.eye(6)[k]), eta_of(np.eye(6)[l])
        S2[k, l] = np.sum(sig * (Ek @ El + El @ Ek)) / 2
fac = np.array([1, 1, 1, 2, 2, 2.])
c_exact = (H[:6, :6] / V - S2) / np.outer(fac, fac) * GPA

conv, prim, atoms = read_structure(f'{run}/structure.params'); ip = np.linalg.inv(prim)
cell0 = []
for t, (sp, p) in enumerate(atoms):
    for nn in range(N):
        fr = (Xi[nn] - p) @ ip
        if np.abs(fr - np.rint(fr)).max() < 1e-6:
            cell0.append(nn); break
Hs = m.hessian(X0).reshape(N, 3, N, 3); il = np.linalg.inv(lat)
A = np.zeros((3, 3, 3, 3))
for i in cell0:
    for j in range(N):
        dv = X0[j] - X0[i]; f = dv @ il; dv = (f - np.rint(f)) @ lat
        A += -0.5 * np.einsum('ag,b,d->abgd', Hs[i, :, j, :], dv, dv)
A *= GPA / Vc
sg = sig * GPA
vo = lambda T: np.array([[T[a, b, c, dd] for (c, dd) in VO] for (a, b) in VO])
sym = lambda T: 0.25 * (T + T.transpose(1, 0, 2, 3) + T.transpose(0, 1, 3, 2) + T.transpose(1, 0, 3, 2))
CBH = np.einsum('abgd->abgd', A) + np.einsum('bagd->abgd', A) - np.einsum('agbd->abgd', A)
dl = np.eye(3)
ST = np.einsum('ag,bd->abgd', dl, sg)
np.set_printoptions(precision=4, suppress=True, linewidth=160)
print(f'{mat} {"relaxed (C=0,D=0)" if relaxed else "reference"}: max|force| {fx:.1e} eV/A, stress (GPa) diag {np.round(np.diag(sg), 5)}')
print('exact c (GPa):'); print(c_exact)
STB = 0.5 * (np.einsum('ab,gd->abgd', dl, sg) + np.einsum('gd,ab->abgd', dl, sg)) \
    - (ST + np.einsum('bd,ag->abgd', dl, sg))
for name, T in (('Sym(A)', sym(A)), ('Born-Huang C_BH', sym(CBH)), ('C_BH - Sym(d_ag s_bd)', sym(CBH - ST)),
                ('C_BH + stress terms', sym(CBH + STB))):
    print(f'{name:24s}: max |C - c_exact| = {np.abs(vo(T) - c_exact).max():.2e} GPa')
if len(sys.argv) > 3:
    print('c_exact - C_BH (Voigt, GPa):'); print(c_exact - vo(sym(CBH)))
    print('Voigt of Sym(d_ag s_bd):'); print(vo(sym(ST)))
    for nm, T in (('d_ab s_gd', np.einsum('ab,gd->abgd', dl, sg)), ('d_ag s_bd+d_bd s_ag', ST + np.einsum('bd,ag->abgd', dl, sg))):
        print(nm); print(vo(sym(T)))
