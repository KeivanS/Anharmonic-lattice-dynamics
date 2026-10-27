"""Exact test model for FOCEX (2026-10-26, user's design): pair + 3-body + 4-body, a polynomial of
degree 4 in the positions, so the force constants of ranks 1..4 are exact and rank >= 5 vanish.
All terms are quadratic forms in the squared distances s = |x_j - x_i|^2 of their atoms:

  E = sum_{pairs r<rL}  K_p s_p                              harmonic, long range (K_p fixed at the reference)
    + sum_{pairs r<Ra}  A (s_p - s0_p)^2                     anharmonic pair (s0_p at the reference)
    + sum_{triplets}    B |AB x AC|^2                        area^2 (x4), all pairs < Ra
    + sum_{quadruplets} C [(AB.CD)^2 + (AC.BD)^2 + (AD.BC)^2] opposite edges, all pairs < Ra
    each anharmonic term minus its part linear in (s - s0) at the reference (still a quadratic form
    in s), so stress and FC1 come from the harmonic term only and scale with K0

The symmetrised face-normal product (AB x AC).(AB x AD) has no 4-atom part (sum of the face normals
of a tetrahedron is 0), hence the opposite-edge form. Term lists are fixed at the reference
structure, so E is a polynomial everywhere: E = L.s + 1/2 s^T H s, and with s quadratic in x
  FC1 = sum g s'            FC2 = sum g s'' + sum H s' s'
  FC3 = sum H (s''_ab s'_c + s''_ac s'_b + s''_bc s'_a)
  FC4 = sum H (s''_ab s''_cd + s''_ac s''_bd + s''_ad s''_bc),   g = L + H s
exactly (s''' = 0)."""
import numpy as np
from itertools import combinations
from scipy.spatial import cKDTree

PAR = dict(K0=0.05, A=0.03, B=0.02, C=0.006)
H3 = np.array([[-2, 2, 2], [2, -2, 2], [2, 2, -2]]) / 16.0               # pairs ab ac bc
_M4 = np.array([[2, -2, -2, -2, -2, 4], [0, 2, -2, -2, 4, -2], [0, 0, 2, 4, -2, -2],
                [0, 0, 0, 2, -2, -2], [0, 0, 0, 0, 2, -2], [0, 0, 0, 0, 0, 2]]) / 4.0
H4 = _M4 + _M4.T                                                         # pairs ab ac ad bc bd cd
Q3 = ((0, 1), (0, 2), (1, 2))
Q4 = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


class Model:
    def __init__(self, rL, Ra, par=PAR):
        self.rL, self.Ra, self.p = rL, Ra, dict(par)

    def build(self, X, lat=None, centers=None):
        """terms on points X; lat = supercell vectors (rows) for periodic images, None for a cluster;
        centers: only terms containing one of them (all if None)."""
        rL, Ra = self.rL, self.Ra
        if lat is not None:
            imgs = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
        else:
            imgs = np.zeros((1, 3), int); lat = np.eye(3)
        Y = np.vstack([X + t @ lat for t in imgs])
        idx = np.tile(np.arange(len(X)), len(imgs)); img = np.repeat(imgs, len(X), axis=0)
        tree = cKDTree(Y)
        centers = range(len(X)) if centers is None else centers
        pairs = {}

        def pid(i, ni, j, nj):
            dn = tuple(int(v) for v in np.asarray(nj) - np.asarray(ni))
            k = pairs.get((i, j, dn))
            if k is None:
                k = pairs.get((j, i, tuple(-v for v in dn)))
            if k is None:
                k = pairs[(i, j, dn)] = len(pairs)
            return k
        nb = {}
        for i in centers:
            nb[i] = []
            for m in tree.query_ball_point(X[i], rL):
                j, nj = idx[m], tuple(img[m])
                if j == i and nj == (0, 0, 0):
                    continue
                pid(i, (0, 0, 0), j, nj)
                if np.linalg.norm(Y[m] - X[i]) < Ra:
                    nb[i].append((j, nj, Y[m]))
        seen, tl, ql = set(), [], []
        for i in centers:
            me = (i, (0, 0, 0), X[i])
            for a, b in combinations(nb[i], 2):
                if np.linalg.norm(a[2] - b[2]) >= Ra:
                    continue
                at = [me, a, b]; key = canon(at)
                if key not in seen:
                    seen.add(key); tl.append([pid(at[p][0], at[p][1], at[q][0], at[q][1]) for p, q in Q3])
            for a, b, c in combinations(nb[i], 3):
                if max(np.linalg.norm(a[2] - b[2]), np.linalg.norm(a[2] - c[2]), np.linalg.norm(b[2] - c[2])) >= Ra:
                    continue
                at = [me, a, b, c]; key = canon(at)
                if key not in seen:
                    seen.add(key); ql.append([pid(at[p][0], at[p][1], at[q][0], at[q][1]) for p, q in Q4])
        P = sorted(pairs.items(), key=lambda kv: kv[1])
        self.pi = np.array([k[0] for k, _ in P]); self.pj = np.array([k[1] for k, _ in P])
        self.ps = np.array([np.array(k[2]) @ lat for k, _ in P]).reshape(-1, 3)
        self.trip = np.array(tl, dtype=int).reshape(-1, 3); self.quad = np.array(ql, dtype=int).reshape(-1, 6)
        s0 = self._s(X)[1]
        pr = self.p
        self.L = np.where(s0 < rL ** 2, pr['K0'] * (1 - s0 / rL ** 2) ** 2, 0.0)
        anh = s0 < Ra ** 2
        # sparse symmetric H as (p, q, h), both orders present
        P_, Q_, V_ = [np.where(anh)[0]], [np.where(anh)[0]], [np.full(anh.sum(), 2 * pr['A'])]
        for T, Hc, c in ((self.trip, H3, pr['B']), (self.quad, H4, pr['C'])):
            n = Hc.shape[0]
            for u in range(n):
                for v in range(n):
                    if Hc[u, v] != 0:
                        P_.append(T[:, u]); Q_.append(T[:, v]); V_.append(np.full(len(T), c * Hc[u, v]))
        self.hp, self.hq, self.hv = (np.concatenate(z) for z in (P_, Q_, V_))
        # centred anharmonic terms: their linear part in s at the reference is removed, so dE/ds at
        # the reference (stress, FC1) comes from the harmonic K_p s_p alone
        np.add.at(self.L, self.hp, -self.hv * s0[self.hq])
        self.const = 0.5 * np.sum(self.hv * s0[self.hp] * s0[self.hq])
        return self

    def _s(self, X):
        d = X[..., self.pj, :] - X[..., self.pi, :] + self.ps
        return d, np.einsum('...k,...k->...', d, d)

    def g(self, s):
        """dE/ds per pair (s of one configuration)"""
        out = self.L.copy()
        np.add.at(out, self.hp, self.hv * s[self.hq])
        return out

    def energy(self, X):
        s = self._s(X)[1]
        if s.ndim == 1:
            return self.L @ s + 0.5 * np.sum(self.hv * s[self.hp] * s[self.hq]) + self.const
        return s @ self.L + 0.5 * np.sum(self.hv * s[:, self.hp] * s[:, self.hq], axis=1) + self.const

    def forces(self, X):
        d, s = self._s(X)
        gd = 2 * self.g(s)[:, None] * d
        F = np.zeros_like(X)
        np.add.at(F, self.pj, -gd); np.add.at(F, self.pi, gd)
        return F

    def hessian(self, X):
        """analytic d2E/dx dx (3N x 3N) at configuration X"""
        d, s = self._s(X); g = self.g(s); N = len(X)
        Hm = np.zeros((N, 3, N, 3)); I = np.eye(3)
        for i, j, gp in zip(self.pi, self.pj, g):
            Hm[i, :, i, :] += 2 * gp * I; Hm[j, :, j, :] += 2 * gp * I
            Hm[i, :, j, :] -= 2 * gp * I; Hm[j, :, i, :] -= 2 * gp * I
        for p, q, h in zip(self.hp, self.hq, self.hv):
            dp, dq = 2 * d[p], 2 * d[q]
            for a, sa in ((self.pj[p], 1), (self.pi[p], -1)):
                for b, sb in ((self.pj[q], 1), (self.pi[q], -1)):
                    Hm[a, :, b, :] += h * sa * sb * np.outer(dp, dq)
        return Hm.reshape(3 * N, 3 * N)

    def stress(self, X, vol):
        """(1/V) dE/d(eta) for x -> (1+eta) x (all atoms, no relaxation), eV/A^3"""
        d, s = self._s(X); g = self.g(s)
        return 2 * np.einsum('p,pa,pb->ab', g, d, d) / vol

    def fc(self, X, tup):
        """analytic FC for tup = ((atom, dir 0..2), ...) at configuration X (rank 1..4; 0 beyond)"""
        d, s = self._s(X)
        r = len(tup)
        if r > 4:
            return 0.0
        atoms = {a for a, _ in tup}
        rel = np.where(np.isin(self.pi, list(atoms)) | np.isin(self.pj, list(atoms)))[0]
        def s1(p, t):
            a, al = t
            return 2 * d[p, al] * (int(self.pj[p] == a) - int(self.pi[p] == a))
        def s2(p, t, u):
            (a, al), (b, be) = t, u
            if al != be:
                return 0.0
            ends = (self.pi[p], self.pj[p])
            if a not in ends or b not in ends:
                return 0.0
            return 2.0 if a == b else -2.0
        val = 0.0
        if r <= 2:
            g = self.g(s)
            for p in rel:
                val += g[p] * (s1(p, tup[0]) if r == 1 else s2(p, tup[0], tup[1]))
            if r == 1:
                return val
        m = np.isin(self.hp, rel) & np.isin(self.hq, rel)
        for p, q, h in zip(self.hp[m], self.hq[m], self.hv[m]):
            if r == 2:
                val += h * s1(p, tup[0]) * s1(q, tup[1])
            elif r == 3:
                a, b, c = tup
                val += h * (s2(p, a, b) * s1(q, c) + s2(p, a, c) * s1(q, b) + s2(p, b, c) * s1(q, a))
            else:
                a, b, c, e = tup
                val += h * (s2(p, a, b) * s2(q, c, e) + s2(p, a, c) * s2(q, b, e) + s2(p, a, e) * s2(q, b, c))
        return val


def canon(at):
    items = sorted((a, tuple(int(v) for v in n)) for a, n, _ in at)
    n0 = np.array(items[0][1])
    return tuple((a, tuple(int(v) for v in np.array(n) - n0)) for a, n in items)
