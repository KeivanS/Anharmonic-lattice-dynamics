"""setup.py MAT RUN: exact-model FOCEX inputs for crystal MAT (folder with structure.params, kpbs.params
from poscar2focex.x) in RUN: SC_SNAPS supercell + snapshots, model3 forces on the exact ideal sites +
the SC_SNAPS displacements, FORCEDISP1 at full precision, structure.params with ranges from the model."""
import json, os, shutil, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import read_structure, read_poscar, shells
from model3 import Model

CFG = {   # supercell (conventional units), rL (harmonic pair range), Ra (anharmonic range)
    'CaF2': dict(sc=[[2, 0, 0], [0, 2, 0], [0, 0, 2]], rL=5.0, Ra=4.2),
    'ZnO':  dict(sc=[[3, 0, 0], [0, 3, 0], [0, 0, 2]], rL=4.72, Ra=3.5),
    'HfO2': dict(sc=[[2, 0, 0], [0, 2, 0], [0, 0, 2]], rL=3.87, Ra=3.22),
}
MASS = {'Ca': 40.078, 'F': 18.9984, 'Zn': 65.409, 'O': 15.9994, 'Hf': 178.49}
mat, run = sys.argv[1], sys.argv[2]
nsnap = int(sys.argv[3]) if len(sys.argv) > 3 else 100
target = float(sys.argv[4]) if len(sys.argv) > 4 else 0.04      # rms displacement (A)
cfg = CFG[mat]; src = os.path.abspath(mat)
os.makedirs(run, exist_ok=True); os.chdir(run)
conv, prim, atoms = read_structure(f'{src}/structure.params')
SL = open(f'{src}/structure.params').read().split('\n')
names = SL[9].split('#')[0].split(); ntyp = len(names)

# shells: each model range must sit in a gap of every atom's distance list
sc = np.array(cfg['sc']); latx = sc @ conv
T = np.array([i * latx[0] + j * latx[1] + k * latx[2] for i in range(-2, 3) for j in range(-2, 3) for k in range(-2, 3)])
nt = np.linalg.norm(T, axis=1); inr = nt[nt > 1e-6].min() / 2
assert cfg['rL'] < inr, (cfg['rL'], inr)
sh = shells(prim, atoms, inr + 0.5)
alld = np.concatenate(sh)
for r in (cfg['rL'], cfg['Ra']):
    gap = np.abs(alld - r).min()
    assert gap > 0.02, f'range {r} only {gap} from a neighbour distance'
nsh2 = [int((s < cfg['rL']).sum()) for s in sh]
nsh3 = [int((s < cfg['Ra']).sum()) for s in sh]

# SC_SNAPS
cell = [SL[0].split('#')[0].strip(), SL[1].split('#')[0].strip(), '1.0 0   # scale, convcoord=0',
        str(ntyp), ' '.join(str(sum(1 for a in atoms if a[0] == t + 1)) for t in range(ntyp)),
        ' '.join(f'{MASS[n]:.4f}' for n in names) + ' ' + ' '.join('0' for _ in names) + '   # masses, charges',
        ' '.join(names)]
for l in SL[19:19 + len(atoms)]:
    cell.append(' '.join(l.split()[2:5]))
open('cell.inp', 'w').write('\n'.join(cell) + '\n')
open('supercell.inp', 'w').write('\n'.join(' '.join(map(str, r)) for r in sc) + '   # conventional units\n')
Tk = 300.0
for it in range(4):
    open('snaps.inp', 'w').write(f'300   # average frequency (1/cm)\n{Tk:g}   # temperature (K)\n{nsnap}   # snapshots\n')
    r = subprocess.run(os.path.expanduser('~/ALADYN/bin/sc_snaps.x'), capture_output=True, text=True)
    lat0, cnt, X0, head = read_poscar('poscar_000')
    U = []
    inv = np.linalg.inv(lat0)
    for k in range(1, nsnap + 1):
        x = read_poscar(f'poscar_{k:03d}')[2]
        df = (x - X0) @ inv
        U.append((df - np.rint(df)) @ lat0)
    U = np.array(U); rms = np.sqrt((U ** 2).mean())
    print(f'T={Tk:.1f} K: rms displacement {rms:.4f} A, max {np.abs(U).max():.4f}')
    if abs(rms / target - 1) < 0.1:
        break
    Tk *= (target / rms) ** 2

# exact ideal sites of the supercell, in the order of poscar_000
assert np.abs(lat0 - latx).max() < 1e-6, np.abs(lat0 - latx).max()
ip = np.linalg.inv(prim)
Xi = np.zeros_like(X0); spec = np.repeat(np.arange(1, len(cnt) + 1), cnt); dev = []
for n, x in enumerate(X0):
    best = None
    for t, (sp, p) in enumerate(atoms):
        fr = (x - p) @ ip; d = np.linalg.norm((fr - np.rint(fr)) @ prim)
        if best is None or d < best[0]:
            best = (d, t, np.rint(fr) @ prim + p)
    assert atoms[best[1]][0] == spec[n]
    dev.append(best[0]); Xi[n] = best[2]
print(f'poscar_000 vs exact sites: {max(dev):.2e} A')

mod = Model(cfg['rL'], cfg['Ra']).build(Xi, latx)
vol = abs(np.linalg.det(latx))
sig = mod.stress(Xi, vol)
F0 = mod.forces(Xi)
info = dict(mat=mat, rL=cfg['rL'], Ra=cfg['Ra'], inradius=inr, T=Tk, rms=float(rms), natom=len(Xi),
            npairs=len(mod.pi), ntrip=len(mod.trip), nquad=len(mod.quad), nsh2=nsh2, nsh3=nsh3,
            stress_eV_A3=sig.tolist(), stress_GPa=(sig * 160.21766).tolist(), maxF0=float(np.abs(F0).max()))
print(json.dumps(info))
with open('FORCEDISP1', 'w') as fd:
    for k in range(nsnap):
        x = Xi + U[k]
        E = mod.energy(x); F = mod.forces(x)
        fd.write('  # POSITION (Ang)     TOTAL FORCE (eV/Ang)\n')
        fd.write(f'  {k + 1}  {E:.12f}  = t,Etot(eV)\n')
        for p, f in zip(x, F):
            fd.write(''.join(f' {v:22.15f}' for v in p) + ''.join(f' {v:24.15e}' for v in f) + '\n')
with open('POSCAR1', 'w') as fp:
    fp.write(f'{mat} exact-model supercell\n1.0\n')
    for v in latx:
        fp.write(''.join(f' {c:22.15f}' for c in v) + '\n')
    fp.write(' '.join(names) + '\n' + ' '.join(map(str, cnt)) + '\nCartesian\n')
    for p in Xi:
        fp.write(''.join(f' {c:22.15f}' for c in p) + '\n')
np.save('U.npy', U); np.save('Xi.npy', Xi); np.save('latx.npy', latx)
shutil.copy(f'{src}/kpbs.params', '.')
shutil.copy(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'example', 'Ge', 'default.params'), '.')
open('latdyn.params', 'w').write('''8 8 8  # kmesh for Brillouin zone integration
0 0 0  1 0 0   # shift in units of mesh (should be less than 1); normal to surface for reflectivity
300 1500.0 # frequency mesh size for dos, omega_max(1/cm)
10        # width of gaussian broadening (1/cm)
.False.   # Verbosity
0 900 31  # Tmin(K), Tmax(K), number of T values
0         # 0 for quantum (BE) occupations; 1 for classical
''')
diel = ['0    # Born charges off: the model has no Coulomb term',
        '1 0 0   # dielectric matrix', '0 1 0', '0 0 1']
for t in range(len(atoms)):
    diel += [f'0 0 0   # Born charge of atom {t + 1}', '0 0 0', '0 0 0']
open('dielectric.params', 'w').write('\n'.join(diel) + '\n')
L = list(SL)
L[3] = '1 1 1 1 0 0 0 0   # ranks 1-4 (exact-model test)'
L[4] = '1 1 0 1            # transl, rot, Huang off (model is stressed), enforce_inv'
L[11] = '1                  # fc2range=1: FC2 shells from the next line'
L[12] = ' '.join(map(str, nsh2)) + f'   # rank-2 shells: r < rL={cfg["rL"]}'
L[13] = ' '.join(map(str, nsh3)) + f'   # rank-3 shells: r < Ra={cfg["Ra"]}'
L[14] = ' '.join(map(str, nsh3)) + f'   # rank-4 shells: r < Ra={cfg["Ra"]}'
open('structure.params', 'w').write('\n'.join(L))
json.dump(info, open('setup.json', 'w'), indent=1)
