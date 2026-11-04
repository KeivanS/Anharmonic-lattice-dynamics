"""compare_relaxed.py RUNDIR MINJSON MAT: FOCEX relax_reference output (log, fc1-3.dat) against the exact
model at its minimum (model_min.py): strain, u0, and every FC1/FC2/FC3 term of fcN.dat evaluated exactly
at the relaxed positions (model terms fixed at the reference; positions X_ref F^T + u0(tau))."""
import sys, os, json, re, numpy as np
from scipy.spatial import cKDTree
EM = os.path.expanduser('~/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/REFERENCE/exact_model')
sys.path.insert(0, EM)
from model3 import Model, PAR
from model import read_structure, cluster
d, mj, mat = sys.argv[1], sys.argv[2], sys.argv[3]
info = json.load(open(f'{EM}/run_{mat}/setup.json')); mn = json.load(open(mj))
conv, prim, atoms = read_structure(f'{d}/structure.params')
log = open([os.path.join(d, f) for f in os.listdir(d) if f.startswith('log') and f.endswith('.dat')][0]).read()
Pn = np.array([[float(v) for v in l.split()[:3]] for l in open(f'{d}/POSCAR_relaxed').read().splitlines()[2:5]])
Fm = np.linalg.solve(prim, Pn).T          # rows: prim @ F^T = Pn
e6 = np.array([Fm[0, 0] - 1, Fm[1, 1] - 1, Fm[2, 2] - 1, Fm[1, 2], Fm[0, 2], Fm[0, 1]])
print('F - F^T (rotation part, should be 0):', np.abs(Fm - Fm.T).max())
u0 = np.array([[float(v) for v in l.split('=')[1].split()[1:]] for l in re.findall(r'RELAX_REFERENCE: tau, u0\(tau\) \(A\) =.*', log)])
print('strain FOCEX', e6, '\n       model', np.array(mn['strain']), '\n max diff', np.abs(e6 - mn['strain']).max())
print('u0 FOCEX', u0.tolist(), '\n   model', mn['u0'], '\n max diff', np.abs(u0 - np.array(mn['u0'])).max())
for l in re.findall(r'RELAX_REFERENCE:.*', log):
    if 'iter' in l or re.match(r'RELAX_REFERENCE:\s+\d', l) or 'max' in l or 'dropped' in l or 'beyond' in l or 'FC4 terms' in l:
        print(' ', l)
# reference atompos (lat_fc_unrelaxed.dat) and exact FCs at the relaxed positions
L = open(f'{d}/lat_fc_unrelaxed.dat').read().split('\n')
k = next(i for i, l in enumerate(L) if 'Neighborshell atoms' in l)
nat = int(L[k + 1].split()[1])
apos = np.array([[float(v) for v in L[k + 2 + i].split()[1:4]] for i in range(nat)])
VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
F = np.eye(3)
for kk, (a, b) in enumerate(VOIGT):
    F[a, b] += mn['strain'][kk]
    if a != b: F[b, a] += mn['strain'][kk]
U0 = np.array(mn['u0'])
R = np.linalg.norm(apos, axis=1).max() + info['rL'] + info['Ra'] + 1
X, spc, tau, nn = cluster(prim, atoms, R)
tree = cKDTree(X)
Xr = X @ F.T + U0[tau]
models, blocks = {}, {}
def fc_block(at):
    key = tuple(at)
    if key in blocks: return blocks[key]
    if at[0] not in models:
        models[at[0]] = Model(info['rL'], info['Ra'], dict(info.get('par', PAR))).build(X, None, [at[0]])
    m = models[at[0]]; r = len(at)
    dd, s = m._s(Xr)
    sg = np.array([(m.pj == a).astype(float) - (m.pi == a).astype(float) for a in at])
    p, q, h = m.hp, m.hq, m.hv; I = np.eye(3)
    if r == 1:
        out = np.einsum('p,p,pa->a', m.g(s), sg[0], 2 * dd)
    elif r == 2:
        out = 2 * I * np.sum(m.g(s) * sg[0] * sg[1])
        w = h * sg[0][p] * sg[1][q]; nz = w != 0
        out = out + np.einsum('n,na,nb->ab', w[nz], 2 * dd[p[nz]], 2 * dd[q[nz]])
    else:
        out = np.zeros((3, 3, 3))
        for (k1, k2, k3), perm in (((0, 1, 2), 'abc'), ((0, 2, 1), 'acb'), ((1, 2, 0), 'bca')):
            w = h * sg[k1][p] * sg[k2][p] * sg[k3][q]; nz = w != 0
            v = 4 * np.einsum('n,nc->c', w[nz], dd[q[nz]])
            out = out + np.einsum(f'{perm}->abc', np.einsum('ab,c->abc', I, v))
    blocks[key] = out
    return out
dist, im = tree.query(apos); assert dist.max() < 1e-6
for r in (1, 2, 3):
    fn = f'{d}/fc{r}.dat'
    worst = big = 0.0; nt = 0; ex = None
    for l in open(fn):
        x = l.split()
        if not x or x[0].startswith('#'): continue
        iat = [int(x[2 + 2 * j]) - 1 for j in range(r)]; ixyz = [int(x[3 + 2 * j]) - 1 for j in range(r)]
        v = float(x[2 + 2 * r])
        e = fc_block([int(im[a]) for a in iat])[tuple(ixyz)]
        nt += 1; big = max(big, abs(e))
        if abs(v - e) > worst: worst = abs(v - e); ex = (iat, ixyz, v, e)
    print(f'rank {r}: {nt} terms of fc{r}.dat: max |FOCEX - exact at the minimum| {worst:.3e} (max |FC| {big:.3g})', ex)
