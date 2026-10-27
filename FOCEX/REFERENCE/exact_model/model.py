"""Exact test model for FOCEX: a degree-4 polynomial potential, so its force constants
of ranks 1..4 are exact and those of rank >= 5 vanish.

  E = sum_pairs  A (rc^2 - r^2)^2                               r < rc
    + B sum_{triplets, all pairs < Rm} [f_ab f_ac + f_ab f_bc + f_ac f_bc]
    + C sum_{quadruplets, all pairs < Rm} [f_ab f_cd + f_ac f_bd + f_ad f_bc]
  f(r) = Rm^2 - r^2

Triplets give FC3/FC4 on 3 distinct atoms, quadruplets FC4 on 4 distinct atoms, pairs the
2-body terms (Phi_iij, Phi_iiij, Phi_iijj) out to rc. The lists of pairs, triplets and
quadruplets are fixed at the reference structure (cutoffs applied there only), so E is a
polynomial of the displacements everywhere: no cutoff can be crossed, any amplitude works."""
import numpy as np
from itertools import combinations
from scipy.spatial import cKDTree

A, B, C = 0.02, 0.05, 0.03


def read_structure(path):
    L = [l.split('#')[0].split('!')[0] for l in open(path).read().split('\n')]
    a, b, c, al, be, ga = [float(x) for x in L[0].split()[:6]]
    p2c = np.array([float(x) for x in L[1].split()[:9]]).reshape(3, 3)   # rows R0i in conv units
    w = L[2].split()
    scal = float(w[0]); convcoord = int(w[1]) if len(w) > 1 and w[1].lstrip('-').isdigit() else 0
    nat = int(L[10].split()[0])
    ca, cb, cg = np.cos(np.radians([al, be, ga])); sg = np.sin(np.radians(ga))
    conv = np.array([[a, 0, 0], [b * cg, b * sg, 0],
                     [c * cb, c * (ca - cb * cg) / sg, 0]])
    conv[2, 2] = np.sqrt(c * c - conv[2, 0] ** 2 - conv[2, 1] ** 2)
    conv *= scal
    prim = p2c @ conv
    atoms = []
    for l in L[19:19 + nat]:
        x = l.split()
        red = np.array([float(v) for v in x[2:5]])
        atoms.append((int(x[1]), red @ (conv if convcoord == 0 else prim)))
    return conv, prim, atoms


def read_poscar(path):
    L = open(path).read().split('\n')
    s = float(L[1].split()[0])
    lat = np.array([[float(v) for v in L[i].split()[:3]] for i in (2, 3, 4)]) * s
    i = 5
    while True:
        try:
            cnt = [int(v) for v in L[i].split()]
            if cnt:
                break
        except ValueError:
            pass
        i += 1
    head = L[:i + 1]
    i += 1
    if L[i].strip()[0] in 'sS':
        i += 1
    cart = L[i].strip()[0] in 'cCkK'
    n = sum(cnt)
    X = np.array([[float(v) for v in L[i + 1 + j].split()[:3]] for j in range(n)])
    X = X * s if cart else X @ lat
    return lat, cnt, X, head


def cluster(prim, atoms, R):
    """sites of the infinite crystal within R of the origin: positions, species, tau, n"""
    n = int(np.ceil(R / np.min(np.linalg.svd(prim)[1]))) + 2
    g = np.arange(-n, n + 1)
    N = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
    T = N @ prim
    P, S, TAU, NN = [], [], [], []
    for t, (sp, p) in enumerate(atoms):
        q = T + p
        k = np.linalg.norm(q, axis=1) < R
        P.append(q[k]); S += [sp] * int(k.sum()); TAU += [t] * int(k.sum()); NN.append(N[k])
    return np.vstack(P), np.array(S), np.array(TAU), np.vstack(NN)


def shells(prim, atoms, rmax, tol=1e-4):
    """distinct neighbour distances (< rmax) of each primitive atom"""
    pts = cluster(prim, atoms, rmax + 1 + max(np.linalg.norm(p) for _, p in atoms))[0]
    out = []
    for _, p in atoms:
        d = np.sort(np.linalg.norm(pts - p, axis=1))
        d = d[(d > tol) & (d < rmax)]
        u = [d[0]]
        for x in d[1:]:
            if x - u[-1] > tol:
                u.append(x)
        out.append(np.array(u))
    return out


class Model:
    """terms of the potential on a set of points X; images given by integer cell offsets
    (lat) for a periodic supercell, none for a finite cluster."""

    def __init__(self, rc, Rm):
        self.rc, self.Rm = rc, Rm

    def build(self, X, lat=None, centers=None):
        rc, Rm = self.rc, self.Rm
        if lat is not None:
            imgs = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
        else:
            imgs = np.zeros((1, 3), int); lat = np.eye(3)
        Y = np.vstack([X + t @ lat for t in imgs])
        idx = np.tile(np.arange(len(X)), len(imgs)); img = np.repeat(imgs, len(X), axis=0)
        tree = cKDTree(Y)
        if centers is None:
            centers = range(len(X))
        pairs = {}
        def pid(i, ni, j, nj):        # pair (i at cell ni, j at cell nj): index in the pair list
            dn = tuple(np.asarray(nj) - np.asarray(ni))
            k = pairs.get((i, j, dn))
            if k is None:
                k = pairs.get((j, i, tuple(-np.array(dn))))
            if k is None:
                k = pairs[(i, j, dn)] = len(pairs)
            return k
        nb = {}
        for i in centers:
            nb[i] = []
            for m in tree.query_ball_point(X[i], rc):
                j, nj = idx[m], tuple(img[m])
                if j == i and nj == (0, 0, 0):
                    continue
                pid(i, (0, 0, 0), j, nj)
                if np.linalg.norm(Y[m] - X[i]) < Rm:
                    nb[i].append((j, nj, Y[m]))
        tri, quad, tl, ql = set(), set(), [], []
        for i in centers:
            me = (i, (0, 0, 0), X[i])
            for a, b in combinations(nb[i], 2):
                if np.linalg.norm(a[2] - b[2]) >= Rm:
                    continue
                at = [me, a, b]
                key = canon(at)
                if key in tri:
                    continue
                tri.add(key)
                tl.append([pid(at[p][0], at[p][1], at[q][0], at[q][1]) for p, q in ((0, 1), (0, 2), (1, 2))])
            for a, b, c in combinations(nb[i], 3):
                if max(np.linalg.norm(a[2] - b[2]), np.linalg.norm(a[2] - c[2]), np.linalg.norm(b[2] - c[2])) >= Rm:
                    continue
                at = [me, a, b, c]
                key = canon(at)
                if key in quad:
                    continue
                quad.add(key)
                ql.append([pid(at[p][0], at[p][1], at[q][0], at[q][1])
                           for p, q in ((0, 1), (2, 3), (0, 2), (1, 3), (0, 3), (1, 2))])
        P = sorted(pairs.items(), key=lambda kv: kv[1])
        self.pi = np.array([k[0] for k, _ in P]); self.pj = np.array([k[1] for k, _ in P])
        self.ps = np.array([np.array(k[2]) @ lat for k, _ in P]).reshape(-1, 3)
        self.trip = np.array(tl, dtype=int).reshape(-1, 3)
        self.quad = np.array(ql, dtype=int).reshape(-1, 6)
        return self

    def _s(self, X):
        d = X[..., self.pj, :] - X[..., self.pi, :] + self.ps
        return d, np.einsum('...k,...k->...', d, d)

    def energy(self, X):
        d, s = self._s(X)
        rc2, R2 = self.rc ** 2, self.Rm ** 2
        E = A * np.sum((rc2 - s) ** 2, axis=-1)
        f = R2 - s
        if len(self.trip):
            fa, fb, fc = (f[..., self.trip[:, c]] for c in range(3))
            E = E + B * np.sum(fa * fb + fa * fc + fb * fc, axis=-1)
        if len(self.quad):
            q = [f[..., self.quad[:, c]] for c in range(6)]
            E = E + C * np.sum(q[0] * q[1] + q[2] * q[3] + q[4] * q[5], axis=-1)
        return E

    def forces(self, X):
        d, s = self._s(X)
        rc2, R2 = self.rc ** 2, self.Rm ** 2
        g = -2 * A * (rc2 - s)          # dE/ds per pair
        f = R2 - s
        for c, (u, v) in enumerate(((1, 2), (0, 2), (0, 1))):
            np.add.at(g, self.trip[:, c], -B * (f[self.trip[:, u]] + f[self.trip[:, v]]))
        for c in range(6):
            np.add.at(g, self.quad[:, c], -C * f[self.quad[:, c ^ 1]])
        gd = 2 * g[:, None] * d
        F = np.zeros_like(X)
        np.add.at(F, self.pj, -gd)
        np.add.at(F, self.pi, gd)
        return F

    def distances(self, X):
        return np.sqrt(self._s(X)[1])


def canon(at):
    """canonical key of a cluster [(atom, cell, pos)]: sorted (atom, cell - cell of the
    smallest member)"""
    items = sorted((a, tuple(n)) for a, n, _ in at)
    n0 = np.array(items[0][1])
    return tuple((a, tuple(np.array(n) - n0)) for a, n in items)
