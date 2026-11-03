"""Brackets of the Born-charge (dipole) force constants alone, E = 0, for the user's ZnO (t_b2 inputs).
Phi_tt'(r) = Z_t T(r) Z_t'^T, T = -grad grad phi, phi = 1/(sqrt(det eps) sqrt(r.eps^-1.r)), e^2 = 14.399645 eV A.
C_tt'(q) = sum_R Phi(R + x_t' - x_t) exp(i q.(R + x_t' - x_t)) by Ewald (parameter lam);
E = 0: the bare G = 0 term (4pi/V) q q/(q eps q) is left out (as FOCEX's ewald_no_macro).
Unrelaxed brackets: sum_tt' C_tt'(t qh) = sum Phi - t^2/2 sum Phi (qh.r)^2 + ..., so
  Gamma_ik(qh) = sum_jl B[(ik),(jl)] qh_j qh_l = (q^2 coefficient) / V.
Fit the 36 B from 13 directions, then the deviatoric 'stress' of the asymmetry B[ab,cd]-B[cd,ab] = d_ab s_cd - d_cd s_ab.
Also Q(t a; i j) = sum_t' Phi^{ai} r_j (linear term): its part antisymmetric in (i j) feeds the relaxation term.
"""
import numpy as np
from math import erfc, sqrt, pi, exp
import sys

E2 = 14.399645
GPA = 160.21766
A = np.array([[3.2755440732, 0, 0], [-1.6377720366, 2.8367043786, 0], [0, 0, 5.2842160216]])  # rows
X = np.array([[-0.0000000002, 1.8911362525, 2.6452346814], [1.6377720368, 0.9455681261, 0.0030428756],
              [-0.0000000002, 1.8911362525, 4.6487052590], [1.6377720368, 0.9455681261, 2.0065974395]])
EPS = np.diag([5.757, 5.757, 6.002])
Z = [np.diag([2.702060, 2.702060, 2.841310]), np.diag([2.701390, 2.701390, 2.841480]),
     np.diag([-2.701550, -2.701550, -2.840940]), np.diag([-2.702600, -2.702600, -2.841000])]
if 'neutral' in sys.argv:                      # remove the small net charge (sum Z = -7e-4 xx)
    zm = sum(Z) / 4; Z = [z - zm for z in Z]
nat = 4
V = abs(np.linalg.det(A))
B = 2 * pi * np.linalg.inv(A).T                # rows = reciprocal vectors
ei = np.linalg.inv(EPS); w, U = np.linalg.eigh(EPS); eih = U @ np.diag(w ** -0.5) @ U.T
sdet = sqrt(np.linalg.det(EPS))


def t_real(r, lam):
    """-grad grad [erfc(lam D)/D]/sqrt(det eps), D = |eps^-1/2 r|"""
    x = eih @ r; rho = np.linalg.norm(x)
    e = erfc(lam * rho); g = 2 * lam / sqrt(pi) * exp(-(lam * rho) ** 2)
    f1 = -e / rho ** 2 - g / rho
    f2 = 2 * e / rho ** 3 + g * (2 / rho ** 2 + 2 * lam ** 2)
    h = np.outer(x, x) / rho ** 2 * (f2 - f1 / rho) + np.eye(3) * f1 / rho
    return -(eih @ h @ eih) / sdet


def kernel(q, lam, nr=6, ng=8):
    """T-kernel blocks K[t][t'](q) (3x3), without Z, E = 0"""
    K = np.zeros((nat, nat, 3, 3), complex)
    rng = range(-nr, nr + 1)
    Rs = [i * A[0] + j * A[1] + k * A[2] for i in rng for j in rng for k in rng]
    Gs = [i * B[0] + j * B[1] + k * B[2] for i in range(-ng, ng + 1) for j in range(-ng, ng + 1) for k in range(-ng, ng + 1)]
    for t in range(nat):
        for u in range(nat):
            d = X[u] - X[t]
            for R in Rs:
                r = R + d
                if np.linalg.norm(r) < 1e-8:
                    continue
                K[t, u] += t_real(r, lam) * np.exp(1j * q @ r)
            for G in Gs:
                k = q + G; kek = k @ EPS @ k
                if kek < 1e-14:
                    continue
                fac = 4 * pi / V / kek * np.exp(-kek / (4 * lam ** 2))
                if np.linalg.norm(G) < 1e-10:
                    fac -= 4 * pi / V / kek          # E = 0: bare G = 0 term removed
                K[t, u] += fac * np.outer(k, k) * np.exp(-1j * G @ d)
            if t == u:
                K[t, u] -= 4 * lam ** 3 / (3 * sqrt(pi)) * ei / sdet   # Gaussian self term (r = 0 excluded)
    return K


def ctot(q, lam):
    K = kernel(q, lam)
    return sum(E2 * Z[t] @ K[t, u] @ Z[u].T for t in range(nat) for u in range(nat))


# 1) Ewald check: independent of lam
q = np.array([0.03, 0.02, 0.05])
c1, c2 = ctot(q, 0.6), ctot(q, 0.9)
print('lam check: max|C(0.6)-C(0.9)| / max|C| =', np.abs(c1 - c2).max() / np.abs(c1).max())

# 2) q^2 coefficient along 13 directions (even part, Richardson from h, 2h)
lam = 0.7; h = 0.01
dirs = [np.eye(3)[i] for i in range(3)]
dirs += [v / np.linalg.norm(v) for v in ([1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1],
                                          [1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1])]
c0 = ctot(np.zeros(3), lam)
voigt = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (1, 2): 3, (2, 1): 3, (0, 2): 4, (2, 0): 4, (0, 1): 5, (1, 0): 5}
rows, rhs = [], []
for qh in dirs:
    e1 = (ctot(h * qh, lam) + ctot(-h * qh, lam)) / 2 - c0
    e2 = (ctot(2 * h * qh, lam) + ctot(-2 * h * qh, lam)) / 2 - c0
    coef = ((16 * e1 - e2) / 12 / h ** 2).real          # q^2 coefficient, error O(h^4)
    gam = coef / V * GPA
    for i in range(3):
        for k in range(i, 3):
            row = np.zeros(36)
            for j in range(3):
                for l in range(3):
                    row[6 * voigt[i, k] + voigt[j, l]] += qh[j] * qh[l]
            rows.append(row); rhs.append(gam[i, k])
M = np.array(rows); print('rank of the 36-bracket design =', np.linalg.matrix_rank(M), ' rows', len(rows))
bb, *_ = np.linalg.lstsq(M, np.array(rhs), rcond=None)
bb = bb.reshape(6, 6)
print('fit misfit (GPa) =', np.sqrt(np.mean((M @ bb.ravel() - rhs) ** 2)))
print('Coulomb-only brackets B[(ik),(jl)] (GPa), Voigt:\n', np.array2string(bb, precision=4, suppress_small=True))
b4 = lambda a, b, c, d: bb[voigt[a, b], voigt[c, d]]
dl = np.eye(3)
# s_dev from the asymmetry, trace 0
pairs = [(a, b) for a in range(3) for b in range(a, 3)]
R, y = [], []
for (a, b) in pairs:
    for (c, d) in pairs:
        row = np.zeros(6)
        row[voigt[c, d]] += dl[a, b]; row[voigt[a, b]] -= dl[c, d]
        R.append(row); y.append(b4(a, b, c, d) - b4(c, d, a, b))
R.append(np.array([1, 1, 1, 0, 0, 0.])); y.append(0)
s, *_ = np.linalg.lstsq(np.array(R), np.array(y), rcond=None)
print('deviatoric stress of the Coulomb brackets (GPa) xx yy zz yz xz xy =', np.round(s, 5))
print('asymmetry-relation misfit (GPa) =', np.sqrt(np.mean((np.array(R) @ s - y) ** 2)))

# 3) Q(t a; i j) from the linear (odd) term: sum_t' C_tt'(t qh) linear part = i t sum Phi (qh.r)
Qa = np.zeros((nat, 3, 3, 3))      # t, a, i(force partner index), j(r)
for j in range(3):
    qh = np.eye(3)[j]
    K1, K2 = kernel(h * qh, lam), kernel(-h * qh, lam)
    K3, K4 = kernel(2 * h * qh, lam), kernel(-2 * h * qh, lam)
    for t in range(nat):
        o1 = sum(E2 * Z[t] @ (K1[t, u] - K2[t, u]) @ Z[u].T for u in range(nat)) / (2 * h)
        o2 = sum(E2 * Z[t] @ (K3[t, u] - K4[t, u]) @ Z[u].T for u in range(nat)) / (4 * h)
        Qa[t, :, :, j] = ((4 * o1 - o2) / 3).imag
anti = max(abs(Qa[t, a, i, j] - Qa[t, a, j, i]) for t in range(nat) for a in range(3) for i in range(3) for j in range(3))
print('max |Q(t a;i j) - Q(t a;j i)| of the Coulomb FCs (eV/A) =', anti, '  max|Q| =', np.abs(Qa).max())
