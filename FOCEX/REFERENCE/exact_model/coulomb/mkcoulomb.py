"""mkcoulomb.py BASE OUT [q1,q2,...] [eps]: exact test of FOCEX born_flag 2 (2026-11-03).
Truth = the exact short-range model of BASE (run_ZnO, run_HfO2: model3, stressed reference, FC1-FC4 exact)
      + point charges q_s per species (POSCAR1 order; default 2,-2) in a scalar dielectric eps, Ewald-summed in the periodic supercell
        (nonlinear in the displacements, rotationally invariant, an energy).
OUT gets BASE's inputs with FORCEDISP1 = model + Coulomb forces (and energies), stress.params = model + Coulomb
stress (sigma = (1/V) dE/deta, as setup.json), dielectric.params: born_flag 2, eps 1, Z = +-q 1.
For scalar Z and eps FOCEX's hybrid subtraction equals the full Coulomb force minus its reference value, so the
fitted short-range FC2-FC4 must equal BASE's exact ones; FC1 shifts by -F_Coulomb(x0).
Checks printed: Ewald parameter independence of E and F, F = -dE/dx by differences."""
import sys, os, json, shutil
import numpy as np
from math import pi, sqrt
from scipy.special import erfc

E2 = 14.399645                       # e^2/(4 pi eps0), eV A
base, out = sys.argv[1], sys.argv[2]
qs = [float(v) for v in sys.argv[3].split(',')] if len(sys.argv) > 3 else [2.0, -2.0]
eps = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
Xi = np.load(f'{base}/Xi.npy'); U = np.load(f'{base}/U.npy'); L = np.load(f'{base}/latx.npy')   # rows = vectors
info = json.load(open(f'{base}/setup.json'))
# species of the supercell atoms: FORCEDISP1/POSCAR1 order; POSCAR1 counts line
P = open(f'{base}/POSCAR1').read().split('\n')
cnt = [int(v) for v in (P[5].split() if P[5].split()[0].isdigit() else P[6].split())]
Q = np.concatenate([np.full(c, qs[i]) for i, c in enumerate(cnt)])
assert len(Q) == len(Xi) and abs(Q.sum()) < 1e-12


def ewald(X, Lat, lam):
    """energy (eV) and forces (eV/A) of charges Q at X in the periodic cell Lat (rows), scalar eps"""
    V = abs(np.linalg.det(Lat)); B = 2 * pi * np.linalg.inv(Lat).T
    rc = 7.0 / lam * sqrt(eps); gc = 14.0 * lam / sqrt(eps)
    n = [int(rc * np.linalg.norm(B[i]) / (2 * pi)) + 2 for i in range(3)]
    I = np.array([(i, j, k) for i in range(-n[0], n[0] + 1) for j in range(-n[1], n[1] + 1) for k in range(-n[2], n[2] + 1)])
    Rs = I @ Lat
    d = X[None, :, :] - X[:, None, :]                       # d[i,j] = x_j - x_i
    E = 0.0; F = np.zeros_like(X)
    for R in Rs:
        r = d + R; rr = np.linalg.norm(r, axis=2)
        m = (rr > 1e-10) & (rr < rc)
        D = rr / sqrt(eps); qq = Q[:, None] * Q[None, :]
        # phi(r) = erfc(lam D)/(D sqrt(det eps)), D = r/sqrt(eps), det eps = eps^3:  = erfc(lam r/sqrt eps)/(r eps)
        with np.errstate(divide='ignore', invalid='ignore'):
            phi = np.where(m, erfc(lam * D) / (rr * eps), 0)
            dphi = np.where(m, (-erfc(lam * D) / rr ** 2 - 2 * lam / sqrt(pi * eps) * np.exp(-(lam * D) ** 2) / rr) / eps, 0)
        E += 0.5 * (qq * phi).sum()
        # force on i: -d/dx_i sum_j q_i q_j phi(|x_j - x_i + R|) = sum_j q_i q_j phi'(r) (r_vec/r)
        with np.errstate(divide='ignore', invalid='ignore'):
            fac = np.where(m, qq * dphi / rr, 0)
        F += (fac[:, :, None] * r).sum(1)
    Ig = np.array([(i, j, k) for i in range(-12, 13) for j in range(-12, 13) for k in range(-12, 13)])
    G = Ig @ B; gn = np.linalg.norm(G, axis=1); G = G[(gn > 1e-10) & (gn < gc)]
    g2 = (G ** 2).sum(1) * eps
    w = 4 * pi / V * np.exp(-g2 / (4 * lam ** 2)) / g2
    S = (Q[None, :] * np.exp(1j * X @ G.T).T).sum(1)        # S(G) = sum q e^{i G x}
    E += 0.5 * (w * abs(S) ** 2).sum()
    ph = np.exp(-1j * X @ G.T)                               # (n, nG)
    F += (Q[:, None, None] * (w[None, :, None] * G[None, :, :] * np.imag(ph * S[None, :])[:, :, None])).sum(1) * (-1)
    E -= lam / sqrt(pi) / eps ** 1.5 * (Q ** 2).sum()        # self term
    return E * E2, F * E2


def stress(X, Lat, lam, h=1e-5):
    V = abs(np.linalg.det(Lat)); s = np.zeros((3, 3))
    for a in range(3):
        for b in range(a, 3):
            M = np.zeros((3, 3)); M[a, b] = M[b, a] = h / (1 if a == b else 2)
            ep = ewald(X @ (np.eye(3) + M).T, Lat @ (np.eye(3) + M).T, lam)[0]
            em = ewald(X @ (np.eye(3) - M).T, Lat @ (np.eye(3) - M).T, lam)[0]
            s[a, b] = s[b, a] = (ep - em) / (2 * h) / V    # off-diagonal: e_ab = e_ba = h/2, dE = sigma_ab h V
    return s


lam = 0.30
E0, F0 = ewald(Xi, L, lam); E0b, F0b = ewald(Xi, L, 0.40)
print('Ewald check: |E(0.30)-E(0.40)| =', abs(E0 - E0b), ' max|F(0.30)-F(0.40)| =', np.abs(F0 - F0b).max())
x = Xi + U[0]; e, f = ewald(x, L, lam)
fd = np.zeros(6); k = 0
for i, a in ((0, 0), (5, 2), (40, 1)):
    xp = x.copy(); xp[i, a] += 1e-5; xm = x.copy(); xm[i, a] -= 1e-5
    print(f'  F = -dE/dx check atom {i} dir {a}: F {f[i, a]:.10f}  -dE/dx {-(ewald(xp, L, lam)[0] - ewald(xm, L, lam)[0]) / 2e-5:.10f}')
print('max |F_Coulomb(x0)| (eV/A) =', np.abs(F0).max(), ' sum F0 =', np.abs(F0.sum(0)).max())
sC = stress(Xi, L, lam) * 160.21766208
print('Coulomb stress (GPa):\n', np.round(sC, 6))
sM = np.array(info['stress_GPa'])

os.makedirs(out, exist_ok=False)
for fn in ('structure.params', 'latdyn.params', 'kpbs.params', 'default.params', 'POSCAR1'):
    shutil.copy(f'{base}/{fn}', out)
S = open(f'{out}/structure.params').read().split('\n')
S[4] = '1 1 1 1            # transl, rot, Huang (stress.params), enforce_inv'
open(f'{out}/structure.params', 'w').write('\n'.join(S))
open(f'{out}/stress.params', 'w').write('\n'.join(' '.join(f'{v:.12f}' for v in r) for r in sM + sC) + '\n')
natp = int(S[10].split('#')[0])
Zs = [qs[int(S[19 + t].split()[1]) - 1] for t in range(natp)]      # atom lines: number, type, coordinates
with open(f'{out}/dielectric.params', 'w') as fd:
    fd.write('2   # born_flag: exact Coulomb test (mkcoulomb.py), point charges in scalar eps\n')
    for i in range(3):
        fd.write(' '.join(f'{eps if i == j else 0:.6f}' for j in range(3)) + ('   # epsilon_inf\n' if i == 0 else '\n'))
    for t in range(natp):
        for i in range(3):
            fd.write(' '.join(f'{Zs[t] if i == j else 0:.6f}' for j in range(3)) + (f'   # Born charge of atom {t+1}\n' if i == 0 else '\n'))
src = open(f'{base}/FORCEDISP1').read().split('\n')
outl = []; k = -1; i = 0; snap = 0
while i < len(src):
    l = src[i]
    if 'POSITION' in l:
        outl.append(l); hdr = src[i + 1].split()
        x = Xi + U[snap]; e, f = ewald(x, L, lam)
        outl.append(f'  {hdr[0]}  {float(hdr[1]) + e - E0:.12f}  = t,Etot(eV)')
        for a in range(len(Xi)):
            v = src[i + 2 + a].split()
            pos = [float(t) for t in v[:3]]; fr = [float(t) + f[a, c] for c, t in enumerate(v[3:6])]
            outl.append(' '.join(f'{t:20.15f}' for t in pos) + '  ' + ' '.join(f'{t:22.15e}' for t in fr))
        i += 2 + len(Xi); snap += 1
    else:
        outl.append(l); i += 1
open(f'{out}/FORCEDISP1', 'w').write('\n'.join(outl))
json.dump(dict(base=os.path.abspath(base), q=qs, eps=eps, lam=lam, F0_coulomb=F0.tolist(),
               stress_coulomb_GPa=sC.tolist(), stress_model_GPa=sM.tolist()), open(f'{out}/coulomb.json', 'w'), indent=1)
print('snapshots', snap, ' written to', out)
