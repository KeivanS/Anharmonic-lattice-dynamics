"""mkrelaxed.py MAT DIR: FOCEX inputs in ../run_MAT/DIR at the relaxed state (P = 0, no forces) of the model of
bptest2.py with D = 0 (model3 with the 4-body term off: energy of degree 4), whose B and B' are known from full
relaxations under pressure (bp2_all.log). Same 100 displacements (U.npy); forces cut after the cubic term
(orders separated exactly, as mkcubic.py); structure.params: relaxed cell and positions, ranks 1-3, Huang on,
enforce_inv=1; stress.params = 0."""
import sys, os, json, shutil
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
from model3 import Model, PAR
from model import read_structure

mat, sub = sys.argv[1], sys.argv[2]
run = os.path.join(here, '..', f'run_{mat}')
d = os.path.join(run, sub)
if os.path.exists(d):
    sys.exit(f'{d} exists')
info = json.load(open(f'{run}/setup.json'))
Xi, latx, U = np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy'), np.load(f'{run}/U.npy')
par = dict(PAR); par['C'] = 0.0
m = Model(info['rL'], info['Ra'], par).build(Xi, latx)
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
for it in range(30):
    g = grad(q)
    if np.abs(g).max() < 1e-11:
        break
    q = q - np.linalg.pinv(hess(q), rcond=1e-10) @ g
print(f'{mat}: relaxed, max|grad| {np.abs(grad(q)).max():.1e}, strain {np.round(q[:6], 6).tolist()}')
F = np.eye(3) + eta_of(q[:6])
Xeq, pseq = conf(q); m.ps = pseq
lat_eq = latx @ F.T
print(f'  max |force| at the relaxed state {np.abs(m.forces(Xeq)).max():.1e} eV/A')

os.mkdir(d)
L = np.array([0, 1, -1, 2, -2]); Vm = np.vander(L, 4, increasing=True); maxsep = 0
with open(f'{d}/FORCEDISP1', 'w') as fd:
    for k, u in enumerate(U):
        Fs = np.array([m.forces(Xeq + l * u) for l in L]).reshape(len(L), -1)
        c, *_ = np.linalg.lstsq(Vm, Fs, rcond=None)
        maxsep = max(maxsep, np.abs(Vm @ c - Fs).max())
        Fc = (c[0] + c[1] + c[2]).reshape(-1, 3)
        fd.write('  # POSITION (Ang)     TOTAL FORCE (eV/Ang)\n')
        fd.write(f'  {k + 1}  0.0  = t,Etot(eV)\n')
        for p, f in zip(Xeq + u, Fc):
            fd.write(''.join(f' {v:22.15f}' for v in p) + ''.join(f' {v:24.15e}' for v in f) + '\n')
print(f'  forces cut after the cubic term; FC5 = 0 check {maxsep:.1e}')

# POSCAR1: relaxed supercell, same species lines
P = open(f'{run}/POSCAR1').read().split('\n')
out = P[:2] + [''.join(f' {c:22.15f}' for c in v) for v in lat_eq] + P[5:8] + \
      [''.join(f' {c:22.15f}' for c in p) for p in Xeq]
open(f'{d}/POSCAR1', 'w').write('\n'.join(out) + '\n')

# structure.params: relaxed conventional cell (6 numbers) and atoms in its reduced coordinates
conv, prim, atoms = read_structure(f'{run}/structure.params')
S = open(f'{run}/structure.params').read().split('\n')
conv_eq = conv @ F.T
a, b, c = np.linalg.norm(conv_eq, axis=1)
ang = lambda x, y: np.degrees(np.arccos(x @ y / np.linalg.norm(x) / np.linalg.norm(y)))
S[0] = f'{a:.15f} {b:.15f} {c:.15f} {ang(conv_eq[1], conv_eq[2]):.15f} {ang(conv_eq[0], conv_eq[2]):.15f} ' \
       f'{ang(conv_eq[0], conv_eq[1]):.15f}  # relaxed conventional cell (mkrelaxed.py)'
ip = np.linalg.inv(prim)
for t, (sp, p) in enumerate(atoms):
    for nn in range(N):                       # an atom of sublattice t: its relaxation u_t
        fr = (Xi[nn] - p) @ ip
        if np.abs(fr - np.rint(fr)).max() < 1e-6:
            ut = q[6 + 3 * nn:9 + 3 * nn]; break
    red = (p @ F.T + ut) @ np.linalg.inv(conv_eq)
    w = S[19 + t].split()
    S[19 + t] = f'{w[0]} {w[1]} ' + ' '.join(f'{x:.15f}' for x in red)
# snap round-off (~1e-12) to the exact values: FOCEX's symmetry analysis needs them to < 1e-10
from fractions import Fraction
w0 = S[0].split()
if abs(float(w0[0]) - float(w0[1])) < 1e-9:
    w0[1] = w0[0]
for k in (3, 4, 5):
    if abs(float(w0[k]) - round(float(w0[k]))) < 1e-6:
        w0[k] = str(round(float(w0[k])))
S[0] = ' '.join(w0[:6]) + '  # relaxed conventional cell (mkrelaxed.py)'
for t in range(len(atoms)):
    w = S[19 + t].split(); y = []
    for v in map(float, w[2:5]):
        f = Fraction(v).limit_denominator(6)
        y.append(f'{float(f):.15f}' if abs(v - float(f)) < 1e-9 else f'{v:.15f}')
    S[19 + t] = f'{w[0]} {w[1]} ' + ' '.join(y)
S[3] = '1 1 1 0 0 0 0 0   # ranks 1-3 (relaxed, cubic-truncated model)'
S[4] = '1 1 1 1            # transl, rot, Huang, enforce_inv'
open(f'{d}/structure.params', 'w').write('\n'.join(S))
open(f'{d}/stress.params', 'w').write('0 0 0\n0 0 0\n0 0 0\n')
for f in ('dielectric.params', 'latdyn.params', 'kpbs.params', 'default.params'):
    shutil.copy(f'{run}/{f}', d)
print(d)
