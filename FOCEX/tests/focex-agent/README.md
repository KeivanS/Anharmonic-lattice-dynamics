# End-to-end test of focex-agent.py without VASP

`run_test.sh` runs the whole workflow (relax, POSCAR0, SCF, snapshots, supercell runs,
FORCEDISP1, FOCEX, plots) on fcc argon, with `fake_vasp.py` in place of VASP:

- `fake_vasp.py` -- energy and forces of a nearest-neighbour Lennard-Jones crystal
  (eps = 0.0104 eV, sigma = 3.40 A, smooth switch-off between 4.2 and 4.6 A). For fcc with
  a = 2^(1/6) sigma sqrt(2) = 5.39716 A the first neighbours sit at the minimum of the pair
  potential and the second ones are beyond the cutoff, so the crystal is at equilibrium with
  zero stress. It writes the OUTCAR lines focex-agent and read_outcar.x read (NIONS, TOTEN,
  POSITION / TOTAL-FORCE, "EDIFF is reached", "reached required accuracy", "General timing").
  Forces checked against -dE/dx (central differences) to 8 digits.
- `exact_phonons.py` -- the exact harmonic phonons of that model,
  D(k) = V''(r0)/m sum_R (1 - cos k.R) e e^T over the 12 neighbours, compared with FOCEX's
  `bs_freq.dat` at the same k-points.

The input POSCAR is a skewed, rotated primitive cell, so the cell reduction is tested too.

Result on 2026-10-22 (FOCEX from main 2711efd plus the mechanical0 fix below, SC_SNAPS
53afb91 built with -O2). SC_SNAPS draws new random snapshots on each run, so the last digits
change from run to run:

| step | result |
| --- | --- |
| relax, poscar0 | Fm-3m, a = 5.39716 A, cell in the FOCEX frame |
| snaps | 2x2x2 conventional, 32 atoms; T 30 K -> 101.3 K for rms 0.0198 A (target 0.02) |
| vasp_snaps, collect | 51 runs; max sum of forces 5e-6 eV/A |
| focex | fit error 0.50-0.53 %, 4896 equations / 37 unknowns, no imaginary modes |
| phonons vs exact | max difference 0.22-0.23 1/cm = 0.41-0.43 % of the highest frequency (52.9 1/cm); FOCEX slightly higher everywhere (quartic terms sampled by the 0.02 A displacements, not fitted) |
| elastic constants vs exact (C11 = 2k/a, C12 = C44 = k/a) | C11 3.069 (3.052), C12 1.538 (1.526), C44 1.537 (1.526), B 2.048 (2.035) GPa |
| plots | bvel.plt, all.plt, qha.plt: ok |

Found by this test:

- `~/ALADYN/bin/sc_snaps.x` (built with `-fcheck=all`) aborted with "Recursive call to
  nonrecursive procedure 'findif'": `d2v` called `findif(r,d1v)` and `d1v` calls `findif`.
  Fixed in SC_SNAPS (`d2v` writes the same central difference out); any flags work now.
- FOCEX stopped with exit code 2 and no message after "ENTERING MECHANICAL0" for a crystal
  with one atom per cell (no `ibz_*`, DOS or `thermo_QHA.dat`): in `mechanical0` the test
  `map(2)%keep(counteri(2,g,ti))` was before the `ti` loop, with `ti` undefined (-fcheck=all:
  index 0). Moved inside the loop (FOCEX/dynamical.f90); Ge, BZO, MoTe2, GeSe outputs are
  unchanged (they do not use mechanical0).
- focex-agent reported that FOCEX step as OK, because the fit error and `bands.dat` existed.
  `focex_once` now also requires exit code 0 (checked with a wrapper that exits 2: the step
  is FAILED with "FOCEX stopped with exit code 2").
