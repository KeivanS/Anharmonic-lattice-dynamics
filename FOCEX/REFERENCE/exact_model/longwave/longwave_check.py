"""longwave_check.py MAT: the long-wave (sound-speed) route to the elastic constants and the stress, on the
exact model (2026-11-02). State: the stressed reference of run_MAT with the internal coordinates relaxed at
fixed lattice (forces 0, stress not 0). From the exact supercell Hessian (minimum image):
  Gamma_ik(qh) V = sum_tt' K2 - (K1 T)^T G (K1 T),  K1 = sum H (qh.d), K2 = -1/2 sum H (qh.d)^2,
  G = pseudo-inverse of K0 outside the 3 translations (the internal relaxation),
along the 13 directions of FOCEX's elastic_longwave; fit of the 36 brackets B[(ik),(jl)]
(Gamma_ik = sum_jl B[(ik),(jl)] qh_j qh_l); stress from the asymmetry
  B[ab,cd] - B[cd,ab] = d_ab s_cd - d_cd s_ab  (deviatoric part only);
C (21) from B - d_ik s_jl, against the exact relaxed c (strain Hessian with the internal coordinates relaxed,
minus the s:eta^2/2 term, as elastic_check.py)."""
import sys, os, json
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
from model3 import Model, PAR
from model import read_structure

mat = sys.argv[1]
run = os.path.join(here, '..', f'run_{mat}')
info = json.load(open(f'{run}/setup.json'))
Xi, latx = np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy')
m = Model(info['rL'], info['Ra'], dict(PAR)).build(Xi, latx)
ps_ref = m.ps.copy(); N = len(Xi); n = 6 + 3 * N
VO = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
V6 = lambda i, j: VO.index(tuple(sorted((i, j))))
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


# relax the internal coordinates at fixed lattice (Newton on the positions only, strain fixed at 0)
conf, grad = make(Xi, ps_ref)
q = np.zeros(n)
for it in range(30):
    g = grad(q)
    if np.abs(g[6:]).max() < 1e-11:
        break
    Hq = hess(lambda x: grad(q + x))
    q[6:] -= np.linalg.pinv(Hq[6:, 6:], rcond=1e-10) @ g[6:]
X0, P0 = conf(q)
conf, grad = make(X0, P0)
H = hess(grad)
V = abs(np.linalg.det(latx))
conv, prim, atoms = read_structure(f'{run}/structure.params'); nat = len(atoms); Vc = V / (N // nat)
m.ps = P0
d, s = m._s(X0)
sig = 2 * np.einsum('p,pa,pb->ab', m.g(s), d, d) / V            # eV/A^3
fx = np.abs(m.forces(X0)).max()
S2 = np.zeros((6, 6))
for k in range(6):
    for l in range(6):
        Ek, El = eta_of(np.eye(6)[k]), eta_of(np.eye(6)[l])
        S2[k, l] = np.sum(sig * (Ek @ El + El @ Ek)) / 2
fac = np.array([1, 1, 1, 2, 2, 2.])
Hee, Heu, Huu = H[:6, :6], H[:6, 6:], H[6:, 6:]
Hrel = Hee - Heu @ np.linalg.pinv(Huu, rcond=1e-10) @ Heu.T       # internal coordinates relaxed
c_exact = (Hrel / V - S2) / np.outer(fac, fac) * GPA

# primitive-cell atoms and the sublattice of every supercell atom
ip = np.linalg.inv(prim); cell0 = []; sub = np.zeros(N, int)
for nn in range(N):
    for t, (sp, p) in enumerate(atoms):
        fr = (Xi[nn] - p) @ ip
        if np.abs(fr - np.rint(fr)).max() < 1e-6:
            sub[nn] = t
            if np.abs(fr).max() < 1e-6:
                cell0.append(nn)
cell0 = sorted(cell0, key=lambda k: sub[k])
Hs = m.hessian(X0).reshape(N, 3, N, 3); il = np.linalg.inv(latx)
dirs = np.array([(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1),
                 (1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)], float)
nd = 3 * nat
rows, rhs = [], []
for qd in dirs:
    qh = qd / np.linalg.norm(qd)
    K0 = np.zeros((nd, nd)); K1 = np.zeros((nd, nd)); K2 = np.zeros((nd, nd))
    for t, i in enumerate(cell0):
        for j in range(N):
            dv = X0[j] - X0[i]; f = dv @ il; dv = (f - np.rint(f)) @ latx
            x = qh @ dv; tp = sub[j]
            K0[3*t:3*t+3, 3*tp:3*tp+3] += Hs[i, :, j, :]
            K1[3*t:3*t+3, 3*tp:3*tp+3] += Hs[i, :, j, :] * x
            K2[3*t:3*t+3, 3*tp:3*tp+3] += -0.5 * Hs[i, :, j, :] * x * x
    T = np.zeros((nd, 3))
    for t in range(nat):
        T[3*t:3*t+3, :] = np.eye(3)
    Pp = np.eye(nd) - T @ T.T / nat
    G = np.linalg.pinv(Pp @ K0 @ Pp, rcond=1e-10)
    Gam = (T.T @ K2 @ T - (K1 @ T).T @ G @ (K1 @ T)) / Vc * GPA
    Gam = (Gam + Gam.T) / 2
    for i in range(3):
        for k in range(i, 3):
            r = np.zeros(36)
            for j in range(3):
                for l in range(3):
                    r[6 * V6(i, k) + V6(j, l)] += qh[j] * qh[l]
            rows.append(r); rhs.append(Gam[i, k])
A = np.array(rows); b = np.array(rhs)
x = np.linalg.lstsq(A, b, rcond=None)[0]; B = x.reshape(6, 6)
B4 = np.array([[[[B[V6(a, b_), V6(c, d_)] for d_ in range(3)] for c in range(3)] for b_ in range(3)] for a in range(3)])
# stress (deviatoric) from the asymmetry; trace fixed to the exact one for the comparison of C below
dl = np.eye(3)
Ms, rs = [], []
for a in range(3):
    for b_ in range(3):
        for c in range(3):
            for d_ in range(3):
                r = np.zeros(6)
                for k, (p_, q_) in enumerate(VO):
                    E = np.zeros((3, 3)); E[p_, q_] = E[q_, p_] = 1
                    r[k] = dl[a, b_] * E[c, d_] - dl[c, d_] * E[a, b_]
                Ms.append(r); rs.append(B4[a, b_, c, d_] - B4[c, d_, a, b_])
Ms.append(np.array([1, 1, 1, 0, 0, 0.])); rs.append(0.0)         # trace 0: the pressure is not measured
sv = np.linalg.lstsq(np.array(Ms), np.array(rs), rcond=None)[0]
s_lw = eta_of(sv) - np.diag(np.diag(eta_of(sv))) + np.diag(sv[:3])
relres = np.linalg.norm(np.array(Ms)[:-1] @ sv - np.array(rs)[:-1]) / max(np.linalg.norm(rs), 1e-30)
sg = sig * GPA; sdev = sg - np.trace(sg) / 3 * dl
s_use = s_lw + np.trace(sg) / 3 * dl
# C from B - d_ik s_jl : B'[(ik),(jl)] = (C_ijkl + C_ilkj)/2
pairs = [(I, J) for I in range(6) for J in range(I, 6)]
pidx = {p: n_ for n_, p in enumerate(pairs)}
Ar, br = [], []
for i in range(3):
    for k in range(3):
        for j in range(3):
            for l in range(3):
                r = np.zeros(21)
                for (I, J) in ((V6(i, j), V6(k, l)), (V6(i, l), V6(k, j))):
                    r[pidx[(min(I, J), max(I, J))]] += 0.5
                Ar.append(r); br.append(B4[i, k, j, l] - dl[i, k] * s_use[j, l])
cfit = np.linalg.lstsq(np.array(Ar), np.array(br), rcond=None)[0]
crel = np.linalg.norm(np.array(Ar) @ cfit - np.array(br)) / np.linalg.norm(br)
C = np.zeros((6, 6))
for (I, J), v in zip(pairs, cfit):
    C[I, J] = C[J, I] = v
np.set_printoptions(precision=4, suppress=True, linewidth=160)
print(f'{mat}: internal coordinates relaxed at fixed lattice: max|force| {fx:.1e} eV/A')
print('exact stress (GPa):'); print(sg)
print('exact deviatoric stress (GPa):'); print(sdev)
print('long-wave deviatoric stress from the asymmetry of the 36 brackets (GPa):'); print(s_lw)
print(f'max |s_lw - s_dev| = {np.abs(s_lw - sdev).max():.2e} GPa; rel. residual of the asymmetry relation {relres:.1e}')
print(f'max |B - B^T| (GPa) = {np.abs(B - B.T).max():.4f}')
print('exact relaxed c (GPa):'); print(c_exact)
print('long-wave C from B - d_ik s_jl (GPa):'); print(C)
print(f'max |C - c_exact| = {np.abs(C - c_exact).max():.2e} GPa; rel. misfit of the 21-C fit {crel:.1e}')
C0 = np.zeros((6, 6))
x0 = np.linalg.lstsq(np.array(Ar), np.array([B4[i, k, j, l] for i in range(3) for k in range(3) for j in range(3) for l in range(3)]), rcond=None)[0]
for (I, J), v in zip(pairs, x0):
    C0[I, J] = C0[J, I] = v
print(f'without the stress term: max |C - c_exact| = {np.abs(C0 - c_exact).max():.2e} GPa')
