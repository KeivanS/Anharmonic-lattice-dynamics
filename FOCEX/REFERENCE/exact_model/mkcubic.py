"""mkcubic.py MAT DIR [NOISE SEED]: FOCEX inputs in run_MAT/DIR whose forces are the exact model's Taylor series cut
after the cubic term, F = -(FC1 + FC2 u + 1/2 FC3 u u), so that a fit of ranks 1-3 reproduces the exact
FC1..FC3. The orders are separated exactly (the forces are a cubic polynomial in u, FC5 = 0):
  F(l u) = c0 + l c1 + l^2 c2 + l^3 c3, least squares over l = 0, +-1, +-2 (the residual checks FC5 = 0).
Same snapshots (U.npy), sites (Xi.npy), supercell (latx.npy) as setup.py; structure.params: ranks 1-3,
Huang rows on with stress.params, enforce_inv=1. NOISE (eV/A): Gaussian noise of that standard deviation
added to every force component (numpy default_rng(SEED)), to test the error estimates of the fit."""
import json, os, sys, shutil
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model3 import Model

CFG = {'CaF2': (5.0, 4.2), 'ZnO': (4.72, 3.5), 'HfO2': (3.87, 3.22)}
mat, sub = sys.argv[1], sys.argv[2]
noise = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
run = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'run_{mat}')
d = os.path.join(run, sub)
if os.path.exists(d):
    sys.exit(f'{d} exists')
os.mkdir(d)
U, Xi, latx = np.load(f'{run}/U.npy'), np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy')
mod = Model(*CFG[mat]).build(Xi, latx)
L = np.array([0, 1, -1, 2, -2])
V = np.vander(L, 4, increasing=True)            # rows: 1, l, l^2, l^3
maxsep = 0
with open(f'{d}/FORCEDISP1', 'w') as fd:
    for k, u in enumerate(U):
        Fs = np.array([mod.forces(Xi + l * u) for l in L]).reshape(len(L), -1)
        c, res, *_ = np.linalg.lstsq(V, Fs, rcond=None)          # 5 points, 4 coefficients: residual = FC5 check
        maxsep = max(maxsep, np.abs(V @ c - Fs).max())
        F = (c[0] + c[1] + c[2]).reshape(-1, 3)
        if noise > 0:
            F = F + rng.normal(0.0, noise, F.shape)
        fd.write('  # POSITION (Ang)     TOTAL FORCE (eV/Ang)\n')
        fd.write(f'  {k + 1}  0.0  = t,Etot(eV)\n')
        for p, f in zip(Xi + u, F):
            fd.write(''.join(f' {v:22.15f}' for v in p) + ''.join(f' {v:24.15e}' for v in f) + '\n')
print(f'{mat}: {len(U)} snapshots, noise {noise} eV/A; max |cubic fit of F(l u) - F| over l = 0,+-1,+-2: {maxsep:.2e} eV/A (FC5 = 0 check)')
for f in ('dielectric.params', 'latdyn.params', 'kpbs.params', 'default.params', 'POSCAR1'):
    shutil.copy(f'{run}/{f}', d)
S = open(f'{run}/structure.params').read().split('\n')
S[3] = '1 1 1 0 0 0 0 0   # ranks 1-3 (cubic-truncated exact model)'
S[4] = '1 1 1 1            # transl, rot, Huang (stress.params), enforce_inv'
open(f'{d}/structure.params', 'w').write('\n'.join(S))
s = np.array(json.load(open(f'{run}/setup.json'))['stress_GPa'])
open(f'{d}/stress.params', 'w').write('\n'.join(' '.join(f'{v:.12f}' for v in r) for r in s) + '\n')
print(d)
