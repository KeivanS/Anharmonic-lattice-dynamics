#!/usr/bin/env python3
"""poscar2focex.py -- POSCAR -> FOCEX structure.params, in a unique, readable cell.

    python3 poscar2focex.py POSCAR [-o structure.params] [--poscar-out POSCAR_focex] [--force]

1. Niggli reduction (Krivy & Gruber 1976, epsilon tests of Grosse-Kunstleve et al. 2004):
   the primitive cell with the shortest lattice vectors -- the same cell whatever primitive
   cell of the crystal you start from.
2. Orientation: the crystal is rotated so that the axes of its conventional cell lie along
   x, y, z (spglib's standard orientation).  A conventional cell that is itself primitive
   (P lattice) and reduced is used as it is; otherwise, of the bases of the reduced cell the
   one that reads best is taken: vector i along axis i, else perpendicular to it, shorter
   vectors first, right-handed, most positive components.  E.g. fcc (0,a,a),(a,0,a),(a,a,0);
   body-centred tetragonal (a,0,0),(0,a,0),(a/2,a/2,c/2); hexagonal (a,0,0),(a/2,a*sqrt3/2,0),c.
3. Writes
   structure.params  line 1: the conventional vectors a1 a2 a3 in Cartesian (9 numbers, A),
                             in the frame of POSCAR_focex
                     line 2: the primitive vectors in units of the conventional ones
                     line 3: 1.0 0  (scale; coordinates in conventional units)
                     atom lines: atom number, atom type, coordinates in units of the
                             CONVENTIONAL vectors
                     the lines in between: masses and names from the POSCAR, the other
                             values the defaults of FOCEX/example/Ge, marked "default, check"
   POSCAR_focex      the reduced, rotated primitive cell (same Cartesian frame as line 1):
                     build the DFT supercells (POSCARi/FORCEDISPi) from this file.

Needs numpy and spglib (pip install spglib).  The atom order of the POSCAR is kept.
"""
import argparse
import itertools
import os
import sys

import numpy as np

MASS = {
    'H': 1.0079, 'He': 4.0026, 'Li': 6.9410, 'Be': 9.0122, 'B': 10.8110, 'C': 12.0107,
    'N': 14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.1797, 'Na': 22.9898, 'Mg': 24.3050,
    'Al': 26.9815, 'Si': 28.0855, 'P': 30.9738, 'S': 32.0650, 'Cl': 35.4530, 'Ar': 39.9480,
    'K': 39.0983, 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670, 'V': 50.9415, 'Cr': 51.9961,
    'Mn': 54.9380, 'Fe': 55.8450, 'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.5460,
    'Zn': 65.4090, 'Ga': 69.7230, 'Ge': 72.6400, 'As': 74.9216, 'Se': 78.9600,
    'Br': 79.9040, 'Kr': 83.7980, 'Rb': 85.4678, 'Sr': 87.6200, 'Y': 88.9059,
    'Zr': 91.2240, 'Nb': 92.9064, 'Mo': 95.9400, 'Tc': 98.0000, 'Ru': 101.0700,
    'Rh': 102.9055, 'Pd': 106.4200, 'Ag': 107.8682, 'Cd': 112.4110, 'In': 114.8180,
    'Sn': 118.7100, 'Sb': 121.7600, 'Te': 127.6000, 'I': 126.9045, 'Xe': 131.2930,
    'Cs': 132.9055, 'Ba': 137.3270, 'La': 138.9055, 'Ce': 140.1160, 'Pr': 140.9076,
    'Nd': 144.2420, 'Pm': 145.0000, 'Sm': 150.3600, 'Eu': 151.9640, 'Gd': 157.2500,
    'Tb': 158.9254, 'Dy': 162.5000, 'Ho': 164.9303, 'Er': 167.2590, 'Tm': 168.9342,
    'Yb': 173.0400, 'Lu': 174.9670, 'Hf': 178.4900, 'Ta': 180.9479, 'W': 183.8400,
    'Re': 186.2070, 'Os': 190.2300, 'Ir': 192.2170, 'Pt': 195.0840, 'Au': 196.9666,
    'Hg': 200.5900, 'Tl': 204.3833, 'Pb': 207.2000, 'Bi': 208.9804, 'Po': 210.0000,
    'At': 210.0000, 'Rn': 220.0000, 'Fr': 223.0000, 'Ra': 226.0000, 'Ac': 227.0000,
    'Th': 232.0381, 'Pa': 231.0359, 'U': 238.0289, 'Np': 237.0000, 'Pu': 244.0000,
    'Am': 243.0000, 'Cm': 247.0000, 'Bk': 247.0000, 'Cf': 251.0000, 'Es': 252.0000,
    'Fm': 257.0000, 'Md': 258.0000, 'No': 259.0000, 'Lr': 262.0000
}


# ── POSCAR ───────────────────────────────────────────────────────────────────
def read_poscar(path):
    L = open(path).read().splitlines()
    s = float(L[1].split()[0])
    A = np.array([[float(x) for x in L[i].split()[:3]] for i in range(2, 5)])
    if s < 0:
        s = (-s / abs(np.linalg.det(A))) ** (1.0 / 3.0)
    A = A * s
    if not L[5].split() or L[5].split()[0][0].isdigit():
        sys.exit('%s: no element line (line 6); a VASP-5 POSCAR is needed' % path)
    species = L[5].split()
    counts = [int(x) for x in L[6].split()]
    i0 = 8 if L[7].strip()[:1] in 'Ss' else 7
    cart = L[i0].strip()[:1] in 'CcKk'
    X = np.array([[float(x) for x in L[i0 + 1 + k].split()[:3]] for k in range(sum(counts))])
    F = X * s @ np.linalg.inv(A) if cart else X
    types = [k for k, n in enumerate(counts) for _ in range(n)]
    return {'title': L[0].strip(), 'A': A, 'species': species, 'counts': counts, 'F': F, 'types': types}


# ── Niggli reduction (lattice vectors as rows) ───────────────────────────────
def _g(L):
    a, b, c = L
    return a @ a, b @ b, c @ c, 2 * b @ c, 2 * a @ c, 2 * a @ b


def _sgn(x, e):
    return 1 if x > e else (-1 if x < -e else 0)


def niggli(A, tol=1e-5):
    """(B, M): Niggli-reduced, right-handed lattice B = M @ A, M integer and unimodular."""
    L = np.array(A, dtype=float)
    e = tol * abs(np.linalg.det(L)) ** (2.0 / 3.0)
    for _ in range(1000):
        A_, B_, C_, xi, eta, zeta = _g(L)
        if A_ > B_ + e or (abs(A_ - B_) <= e and abs(xi) > abs(eta) + e):          # A1
            L = np.array([-L[1], -L[0], -L[2]]); continue
        if B_ > C_ + e or (abs(B_ - C_) <= e and abs(eta) > abs(zeta) + e):        # A2
            L = np.array([-L[0], -L[2], -L[1]]); continue
        s = (_sgn(xi, e), _sgn(eta, e), _sgn(zeta, e))
        want = 1 if s[0] * s[1] * s[2] == 1 else -1                               # A3 / A4
        for i, j, k in itertools.product((1, -1), repeat=3):
            v = (j * k * xi, i * k * eta, i * j * zeta)
            if (want == 1 and all(x > e for x in v)) or (want == -1 and all(x <= e for x in v)):
                L = np.array([i * L[0], j * L[1], k * L[2]]); break
        A_, B_, C_, xi, eta, zeta = _g(L)
        if abs(xi) > B_ + e or (abs(xi - B_) <= e and 2 * eta < zeta - e) or (abs(xi + B_) <= e and zeta < -e):
            L[2] = L[2] - np.sign(xi) * L[1]; continue                             # A5
        if abs(eta) > A_ + e or (abs(eta - A_) <= e and 2 * xi < zeta - e) or (abs(eta + A_) <= e and zeta < -e):
            L[2] = L[2] - np.sign(eta) * L[0]; continue                            # A6
        if abs(zeta) > A_ + e or (abs(zeta - A_) <= e and 2 * xi < eta - e) or (abs(zeta + A_) <= e and eta < -e):
            L[1] = L[1] - np.sign(zeta) * L[0]; continue                           # A7
        t = xi + eta + zeta + A_ + B_
        if t < -e or (abs(t) <= e and 2 * (A_ + eta) + zeta > e):                  # A8
            L[2] = L[2] + L[0] + L[1]; continue
        break
    else:
        sys.exit('Niggli reduction did not converge')
    if np.linalg.det(L) < 0:
        L = -L
    M = np.rint(L @ np.linalg.inv(A)).astype(int)
    if round(abs(np.linalg.det(M))) != 1:
        sys.exit('Niggli reduction gave a non-unimodular transformation')
    return M @ np.array(A, dtype=float), M


# ── orientation and the best-reading basis ───────────────────────────────────
TURNS = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))   # 180-degree turns about x, y, z


def nice_basis(B):
    """(B', N, score): the right-handed basis of the same lattice with the same vector
    lengths that reads best (see the module doc), B' = N @ B."""
    B = np.asarray(B, dtype=float)
    lens = sorted(np.linalg.norm(B, axis=1))
    tol = 1e-4 * lens[-1]
    cand = []
    for c in itertools.product(range(-2, 3), repeat=3):
        v = np.array(c) @ B
        if any(c) and any(abs(np.linalg.norm(v) - l) < tol for l in lens):
            cand.append((np.array(c), v))

    def cost(v, i):                 # along axis i: 0; off-axis and perpendicular to i: 1
        z = np.abs(v) < tol
        if z.sum() == 2:
            return 0 if not z[i] else 2
        return 1 if (z.sum() == 1 and z[i]) else 2
    best = None
    for (n1, v1), (n2, v2), (n3, v3) in itertools.permutations(cand, 3):
        N = np.array([n1, n2, n3])
        if round(np.linalg.det(N)) != 1:
            continue
        V = np.array([v1, v2, v3])
        L = np.linalg.norm(V, axis=1)
        if any(abs(x - y) > tol for x, y in zip(sorted(L), lens)):
            continue
        inv = sum(L[i] > L[j] + tol for i in range(3) for j in range(i + 1, 3))
        score = (cost(v1, 0) + cost(v2, 1) + cost(v3, 2) + 1.5 * inv, -int(np.sum(V > tol)),
                 tuple(-np.round(V, 6).ravel()))
        if best is None or score < best[0]:
            best = (score, N, V)
    return best[2], best[1], best[0][:2]


def dataset(B, F, types):
    import spglib
    for prec in (1e-3, 1e-2, 1e-1):
        ds = spglib.get_symmetry_dataset((B, F, list(types)), symprec=prec)
        if ds is not None:
            return ds, prec
    sys.exit('spglib found no symmetry dataset for this structure')


def get(ds, key):
    return ds[key] if isinstance(ds, dict) else getattr(ds, key)


def canonical(A, F, types):
    """Niggli-reduced, rotated cell: (V, F', how); V rows in Cartesian (A)."""
    B, M = niggli(A)
    F = (F @ np.linalg.inv(M)) % 1.0
    ds, _ = dataset(B, F, types)
    R = np.array(get(ds, 'std_rotation_matrix'))
    Bs = B @ R.T
    # a primitive conventional cell (P lattice) with the reduced lengths is used as it is
    P = np.array(get(ds, 'transformation_matrix'))
    Ni = np.linalg.inv(P).T
    if round(abs(np.linalg.det(P))) == 1 and np.allclose(Ni, np.rint(Ni), atol=1e-6):
        Ni = np.rint(Ni).astype(int)
        Cv = Ni @ Bs
        if np.allclose(sorted(np.linalg.norm(Cv, axis=1)), sorted(np.linalg.norm(B, axis=1)),
                       rtol=1e-6) and np.linalg.det(Cv) > 0:
            Cv[np.abs(Cv) < 1e-10] = 0.0
            F = (F @ np.linalg.inv(Ni)) % 1.0
            F[np.isclose(F, 1.0, atol=1e-10)] = 0.0
            return Cv, F, get(ds, 'international')
    best = None
    for t in TURNS:
        V, N, score = nice_basis(Bs * np.array(t))
        if best is None or score < best[0]:
            best = (score, V, N)
    V, N = best[1].copy(), best[2]
    V[np.abs(V) < 1e-10] = 0.0
    F = (F @ np.linalg.inv(N)) % 1.0
    F[np.isclose(F, 1.0, atol=1e-10)] = 0.0
    return V, F, get(ds, 'international')


def describe(L):
    n = np.linalg.norm(L, axis=1)
    ang = lambda i, j: np.degrees(np.arccos(np.clip(L[i] @ L[j] / n[i] / n[j], -1, 1)))
    return '%.4f %.4f %.4f A, %.2f %.2f %.2f deg' % (n[0], n[1], n[2], ang(1, 2), ang(0, 2), ang(0, 1))


def clean(x, nd=10):
    x = np.round(np.asarray(x, dtype=float), nd)
    x[x == 0] = 0.0                       # no -0.0
    return x


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('poscar')
    ap.add_argument('-o', '--out', default='structure.params')
    ap.add_argument('--poscar-out', default='POSCAR_focex')
    ap.add_argument('--force', action='store_true', help='overwrite existing output files')
    a = ap.parse_args()
    for f in (a.out, a.poscar_out):
        if os.path.exists(f) and not a.force:
            sys.exit('%s exists; use --force to overwrite it' % f)
    try:
        import spglib  # noqa: F401
    except ImportError:
        sys.exit('poscar2focex.py needs spglib:  pip install spglib')

    p = read_poscar(a.poscar)
    V, F, sg = canonical(p['A'], p['F'], p['types'])
    ds, prec = dataset(V, F, p['types'])
    P = np.array(get(ds, 'transformation_matrix'))      # x_conv = P x_prim (+ origin shift)
    C = np.linalg.inv(P).T @ V                           # conventional vectors (rows), same frame
    prim_in_conv = P.T                                   # row i: primitive vector i in units of C
    Xc = F @ P.T                                         # atoms in conventional units, same origin
    assert np.allclose(prim_in_conv @ C, V, atol=1e-8)
    assert np.allclose(Xc @ C, F @ V, atol=1e-8)
    n = len(F)
    nt = len(p['species'])
    masses = []
    for el in p['species']:
        if el.split('_')[0] not in MASS:
            sys.exit('unknown element %r on line 6 of the POSCAR' % el)
        masses.append(MASS[el.split('_')[0]])

    c1 = ' '.join('%.10f' % x for x in clean(C).ravel())
    c2 = ' '.join('%.12g' % x for x in clean(prim_in_conv, 12).ravel())
    per_atom = lambda v: ' '.join([str(v)] * n)
    lines = [
        '%s  # conventional vectors a1 a2 a3 in Cartesian (A); %s' % (c1, describe(C)),
        '%s  # primitive vectors R01 R02 R03 in units of the conventional vectors' % c2,
        '1.0 0              # scale for line 1; 0 = coordinates in units of the conventional vectors',
        '1 1 1 0 0 0 0 0    # default, check: include rank r=1..8 in the fit (0/1)',
        '1 1 0 0            # default, check: invariances: translational, rotational, Huang, enforcement',
        '0 300              # default, check: itemp (0 = no Boltzmann weighting), tempk (K)',
        '1 .false.          # default, check: number of POSCARi/FORCEDISPi pairs, verbosity',
        '%d                  # number of atom types' % nt,
        '%s  # mass of each atom type (amu)' % ' '.join('%.4f' % m for m in masses),
        '%s  # name of each atom type' % ' '.join(p['species']),
        '%d                  # number of atoms in the primitive cell' % n,
        '0                  # default, check: fc2range (0 = FC2 range from the supercell, 1 = next line)',
        '%s  # default, check: shells for rank 2 (one per atom)' % per_atom(11),
        '%s  # default, check: shells for rank 3' % per_atom(3),
    ] + ['%s  # default, check: shells for rank %d' % (per_atom(0), r) for r in range(4, 9)]
    for k in range(n):
        x = clean(Xc[k])
        lines.append('%d %d %s' % (k + 1, p['types'][k] + 1, ' '.join('%.10f' % v for v in x))
                     + ('  # atom number, atom type, coordinates in units of CONVENTIONAL vectors' if k == 0 else ''))
    with open(a.out, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

    out = ['%s (Niggli-reduced, axes along the conventional cell; poscar2focex.py)' % p['title'],
           '1.0'] + ['  %.10f %.10f %.10f' % tuple(clean(v)) for v in V]
    out += ['  ' + '  '.join(p['species']), '  ' + '  '.join(str(c) for c in p['counts']), 'Direct']
    out += ['  %.10f %.10f %.10f' % tuple(clean(f)) for f in F]
    with open(a.poscar_out, 'w') as fh:
        fh.write('\n'.join(out) + '\n')

    print('%s: %s, %d atoms, space group %s (symprec %g)' % (a.poscar, p['title'], n, sg, prec))
    print('  input cell        : %s' % describe(p['A']))
    print('  primitive (Niggli): %s' % describe(V))
    for v in V:
        print('                      %12.6f %12.6f %12.6f' % tuple(clean(v, 6)))
    print('  conventional      : %s' % describe(C))
    print('  wrote %s and %s (build the DFT supercells from %s)' % (a.out, a.poscar_out, a.poscar_out))


if __name__ == '__main__':
    main()
