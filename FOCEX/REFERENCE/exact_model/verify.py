"""verify.py DIR: FOCEX run DIR (inside run_MAT) against the exact FCs of model3.
Needs FOCEX_DUMP_MAP=1 (map_dump.bin) and FOCEX_DUMP_AMAT=k (amat_dump.bin).
 1. x_exact: exact value of every irreducible FC (term listed in fcN_irr.dat, = column order)
 2. map test: every term of every group, sum_ti mat(t,ti) x_exact(col) = exact FC of that term
 3. completeness: every cluster of the model (pairs < rL; pairs, triplets, quadruplets < Ra) with a
    primitive-cell atom is present in the groups of the ranks it contributes to
 4. invariance rows: ahom x_exact - bhom (by kind); force rows: aforce x_exact - bforce
 5. FOCEX's fitted irreducible FCs (fcN_irr.dat) and all fcN.dat values against the exact ones"""
import sys, os, json, numpy as np
from itertools import product
from scipy.spatial import cKDTree
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from model3 import Model
from model import read_structure, cluster

d = sys.argv[1]
info = json.load(open(os.path.join(d, '..', 'setup.json')))
conv, prim, atoms = read_structure(os.path.join(d, 'structure.params'))
nprim = len(atoms)
L = open(os.path.join(d, 'lat_fc.dat' if os.path.exists(os.path.join(d, 'lat_fc.dat')) else 'lat0_fc.dat')).read().split('\n')
k = next(i for i, l in enumerate(L) if 'Neighborshell atoms' in l)
nat = int(L[k + 1].split()[1])
apos = np.array([[float(v) for v in L[k + 2 + i].split()[1:4]] for i in range(nat)])

# ---------------------------------------------------------------- map dump
def rd(f, n, dt):
    return np.frombuffer(f.read(n * np.dtype(dt).itemsize), dtype=dt)
maps = {}
with open(os.path.join(d, 'map_dump.bin'), 'rb') as f:
    while True:
        h = f.read(8)
        if len(h) < 8:
            break
        rnk, ngr = np.frombuffer(h, np.int32)
        G = []
        for g in range(ngr):
            nt, nti = rd(f, 2, np.int32)
            cols = rd(f, nti, np.int32)
            iat = rd(f, rnk * nt, np.int32).reshape(nt, rnk)
            ixyz = rd(f, rnk * nt, np.int32).reshape(nt, rnk)
            mat = rd(f, nt * nti, np.float64).reshape(nti, nt).T
            G.append((cols, iat, ixyz, mat))
        maps[int(rnk)] = G
ranks = sorted(maps)

# ---------------------------------------------------------------- exact FCs (vectorised blocks)
used = sorted({int(a) - 1 for r in ranks for (_, iat, _, _) in maps[r] for a in iat.ravel()})
R = np.linalg.norm(apos[used], axis=1).max() + info['rL'] + info['Ra'] + 1
X, spc, tau, nn = cluster(prim, atoms, R)
tree = cKDTree(X)
dist, im = tree.query(apos[used]); assert dist.max() < 1e-6, dist.max()
cl = {u: int(m) for u, m in zip(used, im)}            # atompos index -> cluster index
models = {}
def model_for(c0):
    if c0 not in models:
        models[c0] = Model(info['rL'], info['Ra']).build(X, None, [c0])
    return models[c0]
blocks = {}
def fc_block(at):
    """exact FC tensor 3^r for the cluster atoms at (first one the centre)"""
    key = tuple(at)
    if key in blocks:
        return blocks[key]
    m = model_for(at[0]); r = len(at)
    dd, s = m._s(X)
    sg = np.array([(m.pj == a).astype(float) - (m.pi == a).astype(float) for a in at])   # r x npairs
    if r == 1:
        out = np.einsum('p,p,pa->a', m.g(s), sg[0], 2 * dd)
    elif r > 4:
        out = np.zeros((3,) * r)
    else:
        sel = np.all(sg[:, m.hp] != 0, axis=0) | True
        p, q, h = m.hp, m.hq, m.hv
        if r == 2:
            g = m.g(s)
            out = 2 * np.eye(3) * np.sum(g * sg[0] * sg[1])
            w = h * sg[0][p] * sg[1][q]; nz = w != 0
            out = out + np.einsum('n,na,nb->ab', w[nz], 2 * dd[p[nz]], 2 * dd[q[nz]])
        elif r == 3:
            out = np.zeros((3, 3, 3)); I = np.eye(3)
            for (k1, k2, k3), perm in (((0, 1, 2), 'abc'), ((0, 2, 1), 'acb'), ((1, 2, 0), 'bca')):
                w = h * sg[k1][p] * sg[k2][p] * sg[k3][q]; nz = w != 0
                v = 4 * np.einsum('n,nc->c', w[nz], dd[q[nz]])          # 2 delta * 2 d
                T = np.einsum('ab,c->abc', I, v)                         # indices (k1,k2,k3)
                out = out + np.einsum(f'{perm}->abc', T)
        else:
            out = np.zeros((3, 3, 3, 3)); I = np.eye(3)
            for (k1, k2, k3, k4) in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)):
                w = np.sum(h * sg[k1][p] * sg[k2][p] * sg[k3][q] * sg[k4][q])
                T = 4 * w * np.einsum('ab,cd->abcd', I, I)               # indices (k1,k2,k3,k4)
                inv = np.argsort([k1, k2, k3, k4])
                out = out + np.transpose(T, inv)
    blocks[key] = out
    return out
def exact(iat, ixyz):
    at = [cl[int(a) - 1] for a in iat]
    return fc_block(at)[tuple(int(x) - 1 for x in ixyz)]

# ---------------------------------------------------------------- 1. x_exact from fcN_irr.dat
def read_fc(fn, r):
    out = []
    if not os.path.exists(fn):
        return out
    for l in open(fn):
        x = l.split()
        if not x or x[0].startswith('#'):
            continue
        out.append(([int(x[2 + 2 * j]) for j in range(r)], [int(x[3 + 2 * j]) for j in range(r)], float(x[2 + 2 * r])))
    return out
ncol = 0; colrank = []
irr = {}
for r in ranks:
    irr[r] = read_fc(os.path.join(d, f'fc{r}_irr.dat'), r)
xe, xf = [], []
for r in ranks:
    for iat, ixyz, v in irr[r]:
        xe.append(exact(iat, ixyz)); xf.append(v); colrank.append(r)
xe, xf, colrank = np.array(xe), np.array(xf), np.array(colrank)
ncols_map = max(int(c) for r in ranks for (cols, _, _, _) in maps[r] for c in cols)
fitted = len(xe) == ncols_map
if not fitted:
    # run still fitting (no fcN_irr.dat yet): x_exact per group by least squares of mat x = exact FCs of
    # all its terms; the map test below is then the residual of that fit (zero iff representable)
    xe = np.zeros(ncols_map); colrank = np.zeros(ncols_map, int)
    for r in ranks:
        for (cols, iat, ixyz, mat) in maps[r]:
            kept = cols > 0
            if kept.any():
                e = np.array([exact(iat[t], ixyz[t]) for t in range(len(iat))])
                xe[cols[kept] - 1] = np.linalg.lstsq(mat[:, kept], e, rcond=None)[0]
                colrank[cols[kept] - 1] = r
    print(f'{d}: no complete fcN_irr.dat; x_exact from the map by least squares per group')
print(f'{d}: columns {len(xe)}, {ncols_map} (map), ranks {ranks}')
np.save(os.path.join(d, 'x_exact.npy'), xe); np.save(os.path.join(d, 'colrank.npy'), colrank)

# ---------------------------------------------------------------- 2. map test
for r in ranks:
    nterm = 0; worst = 0.0; big = 0.0; bad = 0; ex = []
    for (cols, iat, ixyz, mat) in maps[r]:
        kept = cols > 0
        for t in range(len(iat)):
            rec = np.dot(mat[t, kept], xe[cols[kept] - 1]) if kept.any() else 0.0
            e = exact(iat[t], ixyz[t]) if kept.any() else 0.0
            nterm += 1; err = abs(rec - e); big = max(big, abs(e))
            if err > worst:
                worst = err
            if err > 1e-9 * max(1.0, abs(e)):
                bad += 1
                if len(ex) < 3:
                    ex.append((iat[t].tolist(), ixyz[t].tolist(), rec, e))
    print(f'  map rank {r}: {len(maps[r])} groups, {nterm} terms: max |sum mat x_exact - exact| {worst:.2e} '
          f'(max |FC| {big:.3g}), wrong {bad}', ex if ex else '')

# ---------------------------------------------------------------- 3. completeness
present = {r: set() for r in ranks}
for r in ranks:
    for (cols, iat, ixyz, mat) in maps[r]:
        if (cols > 0).any():
            for t in range(len(iat)):
                present[r].add(frozenset(int(a) - 1 for a in iat[t]))
back = {v: k for k, v in cl.items()}
missing = {r: 0 for r in ranks}; checked = {r: 0 for r in ranks}; notinlist = 0
for t0 in range(nprim):
    c0 = cl[t0]
    m = model_for(c0)
    pr = [(m.pi[p], m.pj[p], np.sqrt(m._s(X)[1][p])) for p in range(len(m.pi))]
    sets = []
    for i, j, dist in pr:
        if c0 in (i, j):
            sets.append(({i, j}, [2] if dist >= info['Ra'] else [2, 3, 4]))
    for T in m.trip:
        sets.append(({m.pi[p] for p in T} | {m.pj[p] for p in T}, [3, 4]))
    for Qd in m.quad:
        sets.append(({m.pi[p] for p in Qd} | {m.pj[p] for p in Qd}, [4]))
    for s_, rr in sets:
        if any(a not in back for a in s_):
            notinlist += 1; continue
        fs = frozenset(back[a] for a in s_)
        for r in rr:
            if r in ranks:
                checked[r] += 1
                if fs not in present[r]:
                    missing[r] += 1
print('  completeness: model clusters checked per rank', checked, 'missing from the groups', missing,
      'atoms outside atompos', notinlist)

# ---------------------------------------------------------------- 4. invariance and force rows
with open(os.path.join(d, 'amat_dump.bin'), 'rb') as f:
    nh, ntr, nro, nhu, nf, nc = np.frombuffer(f.read(24), np.int32)
    A = np.frombuffer(f.read(8 * nh * nc), np.float64).reshape(nc, nh).T
    bh = np.frombuffer(f.read(8 * nh), np.float64)
    Af = np.frombuffer(f.read(8 * nf * nc), np.float64).reshape(nc, nf).T
    bf = np.frombuffer(f.read(8 * nf), np.float64)
assert nc == len(xe)
res = A @ xe - bh; sc = np.abs(A) @ np.abs(xe) + np.abs(bh) + 1e-30
o = 0
for name, n in (('translational', ntr), ('rotational', nro), ('Huang', nhu)):
    if n:
        print(f'  {name:13s} rows {n:6d}: max |row.x_exact - b| {np.abs(res[o:o + n]).max():.2e}, '
              f'relative to max row scale {np.abs(res[o:o + n]).max() / sc[o:o + n].max():.2e}, max |b| {np.abs(bh[o:o + n]).max():.3g}')
    o += n
rf = Af @ xe - bf
print(f'  force rows {nf}: max |A x_exact - F| {np.abs(rf).max():.2e} eV/A (max |F| {np.abs(bf).max():.3f}), '
      f'rms {np.sqrt(np.mean(rf ** 2)):.2e}')

# ---------------------------------------------------------------- 5. FOCEX's fit
if not fitted:
    sys.exit()
for r in ranks:
    msk = colrank == r
    if msk.any():
        e = np.abs(xf[msk] - xe[msk]); sc_ = np.abs(xe[msk]).max()
        print(f'  fit fc{r}_irr.dat: {msk.sum()} terms, max |FOCEX - exact| {e.max():.3g} '
              f'(max |FC| {sc_:.3g}, rel {e.max() / max(sc_, 1e-30):.2g})')
for r in ranks:
    T = read_fc(os.path.join(d, f'fc{r}.dat'), r)
    if T:
        e = np.array([abs(v - exact(a, x)) for a, x, v in T])
        print(f'  fit fc{r}.dat: {len(T)} terms, max |FOCEX - exact| {e.max():.3g}')

# ---------------------------------------------------------------- 6. support (for LASSO runs)
# the LASSO selection rule of lasso_constrained.f90: |w_j x_j| > 1e-8 max |w x|, w_j = rms of column j of
# the force rows (here those of the dumped snapshots, so w is approximate). The written x is not
# sparse (off the support |w x| ~ 1e-8 max, or round-off after the refit): counting non-zeros
# cannot see a dropped term.
w = np.sqrt(np.mean(Af ** 2, 0)); w = np.maximum(w, 1e-8 * w.max())
nz_e = np.abs(xe) > 1e-9 * np.abs(xe).max()
zf = w * xf; sel = np.abs(zf) > 1e-8 * np.abs(zf).max()
for r in ranks:
    msk = colrank == r
    e = np.abs(xf[msk] - xe[msk])
    print(f'  support fc{r}: exact non-zero {nz_e[msk].sum()} of {msk.sum()}; selected by FOCEX {sel[msk].sum()}; '
          f'exact non-zero not selected {np.sum(nz_e[msk] & ~sel[msk])}; '
          f'max err on exact zeros {e[~nz_e[msk]].max() if (~nz_e[msk]).any() else 0:.3g}, on non-zeros {e[nz_e[msk]].max() if nz_e[msk].any() else 0:.3g}')
