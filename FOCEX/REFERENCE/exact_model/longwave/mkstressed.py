"""mkstressed.py MAT DIR: FOCEX inputs at a force-free but stressed state of the exact model (2026-11-02):
the reference of run_MAT with the internal coordinates relaxed at fixed lattice (as longwave_check.py).
Same 100 displacements (U.npy), exact model forces (ranks 1-4 exact); structure.params with the relaxed
atom positions, stress.params = the exact stress there (GPa), Huang rows on, enforce_inv=1, relaxation off."""
import sys, os, json, shutil
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
from model3 import Model, PAR
from model import read_structure

mat, d = sys.argv[1], sys.argv[2]
run = os.path.join(here, '..', f'run_{mat}')
if os.path.exists(d):
    sys.exit(f'{d} exists')
info = json.load(open(f'{run}/setup.json'))
Xi, latx, U = np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy'), np.load(f'{run}/U.npy')
m = Model(info['rL'], info['Ra'], dict(PAR)).build(Xi, latx)
N = len(Xi)
X = Xi.copy()
for it in range(30):                            # Newton on the positions, lattice fixed
    f = m.forces(X).ravel()
    if np.abs(f).max() < 1e-12:
        break
    H = m.hessian(X).reshape(3 * N, 3 * N)
    X = X + (np.linalg.pinv(H, rcond=1e-10) @ f).reshape(N, 3)
fx = np.abs(m.forces(X)).max()
V = abs(np.linalg.det(latx))
d_, s = m._s(X)
sig = 2 * np.einsum('p,pa,pb->ab', m.g(s), d_, d_) / V * 160.21766208
print(f'{mat}: internal coordinates relaxed, max|force| {fx:.1e} eV/A; stress (GPa):'); print(sig)
os.mkdir(d)
with open(f'{d}/FORCEDISP1', 'w') as fd:
    for k, u in enumerate(U):
        F = m.forces(X + u)
        fd.write('  # POSITION (Ang)     TOTAL FORCE (eV/Ang)\n')
        fd.write(f'  {k + 1}  0.0  = t,Etot(eV)\n')
        for p, f in zip(X + u, F):
            fd.write(''.join(f' {v:22.15f}' for v in p) + ''.join(f' {v:24.15e}' for v in f) + '\n')
P = open(f'{run}/POSCAR1').read().split('\n')
out = P[:8] + [''.join(f' {c:22.15f}' for c in p) for p in X]
open(f'{d}/POSCAR1', 'w').write('\n'.join(out) + '\n')
conv, prim, atoms = read_structure(f'{run}/structure.params')
S = open(f'{run}/structure.params').read().split('\n')
ip = np.linalg.inv(prim); ic = np.linalg.inv(conv)
for t, (sp, p) in enumerate(atoms):
    for nn in range(N):
        fr = (Xi[nn] - p) @ ip
        if np.abs(fr - np.rint(fr)).max() < 1e-6:
            ut = X[nn] - Xi[nn]; break
    red = (p + ut) @ ic
    w = S[19 + t].split()
    S[19 + t] = f'{w[0]} {w[1]} ' + ' '.join(f'{x:.15f}' for x in red)
S[4] = '1 1 1 1            # transl, rot, Huang (stress.params), enforce_inv'
open(f'{d}/structure.params', 'w').write('\n'.join(S))
open(f'{d}/stress.params', 'w').write('\n'.join(' '.join(f'{v:.12f}' for v in r) for r in sig) + '\n')
for f in ('dielectric.params', 'latdyn.params', 'kpbs.params', 'default.params'):
    shutil.copy(f'{run}/{f}', d)
print(d)
