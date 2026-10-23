"""Stand-in for vasp_std, for testing focex-agent.py without VASP.

Run in a VASP run folder (POSCAR, INCAR): writes OUTCAR (and CONTCAR) with the markers the
agent and read_outcar.x read, the forces and energy of a Lennard-Jones crystal with only
nearest neighbours:  V(r) = 4 eps [(sig/r)^12 - (sig/r)^6] * s(r), s a smooth switch from
1 at R_ON to 0 at R_C.  For fcc with a = 2^(1/6) sig sqrt(2) the first neighbours sit at the
minimum of V and the second ones (at a) are beyond R_C: the crystal is at equilibrium with
zero stress, so a relaxation only copies POSCAR to CONTCAR.

    VASP_STD=/path/to/fake_vasp  (a script that runs this file with python)
"""
import math
import os
import sys

import numpy as np

EPS, SIG = 0.0104, 3.40          # eV, A (argon)
R_ON, R_C = 4.2, 4.6             # A: switch region


def switch(r):
    """s(r) and ds/dr: 1 below R_ON, 0 above R_C, smooth (C2) in between."""
    if r <= R_ON:
        return 1.0, 0.0
    if r >= R_C:
        return 0.0, 0.0
    x = (r - R_ON) / (R_C - R_ON)
    s = 1 - 10 * x**3 + 15 * x**4 - 6 * x**5
    ds = (-30 * x**2 + 60 * x**3 - 30 * x**4) / (R_C - R_ON)
    return s, ds


def pair(r):
    """V(r) and dV/dr."""
    sr6 = (SIG / r) ** 6
    v = 4 * EPS * (sr6 * sr6 - sr6)
    dv = 4 * EPS * (-12 * sr6 * sr6 + 6 * sr6) / r
    s, ds = switch(r)
    return v * s, dv * s + v * ds


def read_poscar(path):
    L = open(path).read().splitlines()
    scale = float(L[1].split()[0])
    A = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)])
    if scale < 0:
        scale = (-scale / abs(np.linalg.det(A))) ** (1 / 3)
    A = A * scale
    species = L[5].split()
    counts = [int(x) for x in L[6].split()]
    i0 = 8 if L[7].strip()[:1] in 'Ss' else 7
    cart = L[i0].strip()[:1] in 'CcKk'
    n = sum(counts)
    X = np.array([[float(x) for x in L[i0 + 1 + k].split()[:3]] for k in range(n)])
    X = X * scale if cart else X @ A
    return L, A, species, counts, X


def energy_forces(A, X):
    """Total energy and forces (eV, eV/A) of the periodic crystal (cell rows A, positions X)."""
    n = len(X)
    inv = np.linalg.inv(A)
    E, F = 0.0, np.zeros_like(X)
    # images needed to cover R_C
    h = [int(math.ceil(R_C * np.linalg.norm(inv[:, i]))) + 1 for i in range(3)]
    shifts = [np.array([a, b, c]) @ A for a in range(-h[0], h[0] + 1)
              for b in range(-h[1], h[1] + 1) for c in range(-h[2], h[2] + 1)]
    for i in range(n):
        for j in range(n):
            d0 = X[j] - X[i]
            for s in shifts:
                d = d0 + s
                r = np.linalg.norm(d)
                if r < 1e-8 or r >= R_C:
                    continue
                v, dv = pair(r)
                E += 0.5 * v
                F[i] += dv * d / r          # force on i from j: -dV/dr_i = dV/dr * d/r
    return E, F


def main():
    incar = open('INCAR').read().upper() if os.path.isfile('INCAR') else ''
    tags = {}
    for line in incar.splitlines():
        line = line.split('!')[0].split('#')[0]
        if '=' in line:
            k, v = line.split('=', 1)
            tags[k.strip()] = v.strip()
    nsw = int(float(tags.get('NSW', '0') or 0))
    ibrion = int(float(tags.get('IBRION', '-1') or -1))
    relax = nsw > 0 and ibrion in (1, 2, 3)
    dfpt = tags.get('LEPSILON', '').startswith('.T') or tags.get('IBRION', '') in ('7', '8')
    L, A, species, counts, X = read_poscar('POSCAR')
    E, F = energy_forces(A, X)
    out = [' vasp.fake (focex-agent test: nearest-neighbour Lennard-Jones)',
           f'   number of dos      NEDOS =    301   number of ions     NIONS = {len(X):6d}',
           ' ------------------------ aborting loop because EDIFF is reached ----------------------------',
           f'  free  energy   TOTEN  = {E:20.8f} eV',
           ' POSITION                                       TOTAL-FORCE (eV/Angst)',
           ' -----------------------------------------------------------------------------------']
    out += [' %12.5f %12.5f %12.5f    %14.6f %14.6f %14.6f' % (*x, *f) for x, f in zip(X, F)]
    out += [' -----------------------------------------------------------------------------------']
    if dfpt:
        out += [' MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)',
                ' ------------------------------------------------------',
                '           1.000000     0.000000     0.000000',
                '           0.000000     1.000000     0.000000',
                '           0.000000     0.000000     1.000000',
                ' ------------------------------------------------------',
                ' BORN EFFECTIVE CHARGES (including local field effects) (in |e|, cummulative output)',
                ' -------------------------------------------------------------------------------']
        for k in range(len(X)):
            out += [f' ion {k + 1:4d}', '    1     0.00000     0.00000     0.00000',
                    '    2     0.00000     0.00000     0.00000', '    3     0.00000     0.00000     0.00000']
    if relax:
        out.append(' reached required accuracy - stopping structural energy minimisation')
    out += [' General timing and accounting informations for this job:',
            ' ========================================================', '',
            '                  Total CPU time used (sec):        0.010',
            '                            Elapsed time (sec):        0.010']
    open('OUTCAR', 'w').write('\n'.join(out) + '\n')
    open('CONTCAR', 'w').write('\n'.join(L) + '\n')
    open('OSZICAR', 'w').write(f'   1 F= {E:.8E} E0= {E:.8E}  d E =0.0\n')
    print('fake vasp: %d atoms, E = %.6f eV, max |F| = %.3e eV/A%s' %
          (len(X), E, np.abs(F).max(), ' (relaxation)' if relax else ''))


if __name__ == '__main__':
    sys.exit(main())
