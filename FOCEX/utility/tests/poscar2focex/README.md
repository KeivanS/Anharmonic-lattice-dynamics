# Tests of poscar2focex.x

`poscar2focex.x` is checked against pymatgen (`SpacegroupAnalyzer` and
`KPathSetyawanCurtarolo`), which implements the Setyawan-Curtarolo convention the program
follows. The scripts are Python only because pymatgen is the reference; the program itself
needs no library.

```bash
./run_tests.sh ~/vaspenv/bin/python ~/ALADYN/bin/poscar2focex.x
```

- `make_structures.py DIR` -- 26 structures, one or more per Setyawan-Curtarolo lattice:
  CUB (BaZrO3, FeS2 pyrite), FCC (Ge, NaCl), BCC, TET (rutile), BCT1, BCT2, ORC (GeSe),
  ORCF1, ORCF2, ORCF3 (the boundary 1/a^2 = 1/b^2 + 1/c^2), ORCI, ORCC, HEX (wurtzite, MoTe2,
  CdI2), RHL1 (Bi), RHL2, MCL, MCLC1, MCLC3, MCLC5, TRI1a, TRI1b. MCLC2 and MCLC4 are not
  covered: each needs an exact equality between lattice parameters, which pymatgen tests with
  floating-point `==`.
- `compare_pymatgen.py X DIR` -- for every structure: sub-case name, path (sequence of
  labels), conventional and primitive lattice parameters, reduced coordinates and length |k|
  of every path point. Prints `ALL OK` when all agree.
- `invariance.py X DIR` -- every structure rewritten 3 times with a random unimodular change
  of basis, a random rotation, a random origin shift and 1e-6 A of noise; `structure.params`
  (cell, primitive vectors, atom lines) and `kpbs.params` must be identical to those of the
  original POSCAR. Prints `ALL INVARIANT`.

Result on 2026-10-22 (gfortran 15.3, `-O2 -fcheck=all`): `ALL OK` (26 structures) and
`ALL INVARIANT` (78 inputs).

Also checked by hand: a skewed, rotated Ge primitive cell gives the `structure.params` of
`FOCEX/example/Ge` (cell, primitive vectors, atoms at (0,0,0) and (1/4,1/4,1/4)); FOCEX run
with it on the example data gives byte-identical `fc2.dat`, `fc3.dat`, `fc2_irr.dat`,
`fc3_irr.dat`, `lat_fc.dat` and `mech.dat` (fit error 0.036 %), and FOCEX reads the
generated FCC `kpbs.params` (X = (0,1,0), W = (1/2,1,0), K = (3/4,3/4,0), L = (1/2,1/2,1/2),
U = (1/4,1,1/4) in conventional units).
