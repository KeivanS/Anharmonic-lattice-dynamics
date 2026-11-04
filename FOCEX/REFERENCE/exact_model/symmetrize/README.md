# Positions off the symmetry: why Huang and rotational rows failed on DFT data (2026-11-04)

Requested by K. Esfarjani (GaN and ZnO: imaginary modes and C66 -7000 GPa with Huang; ZnO fits of 59 %
with enforce_inv 1 or 2). Diagnosis, fix and tests by Claude. Every number below was printed by FOCEX or by
the scripts in this folder.

## What was wrong

- DFT-relaxed positions are a few 1e-5 A off the space group. GaN (structure.params): Ga z 0.5005917 and
  0.0005866 differ by 0.5000051 instead of 0.5, so the 6_3 screw is off by 2.7e-5 A. ZnO: 8.4e-5 A.
- `symmetry_init` accepts an operation when every atom is mapped within tol_sym = 1e-3 A. The force constants
  are then symmetric exactly (built from the irreducible ones with the exact operations).
- The Huang and rotational rows are built from the actual positions (`rr = atompos(j) - atompos(i)`). Rows
  that vanish by symmetry pick up coefficients of the size of the position error.
  - GaN (6mm): [xy,xz], [yz,xx], [yz,yy] are odd under C2z. Their coefficients were 1.4e-3 (largest row 334).
  - These rows are above the SVD cut (1e-10 relative), so they act as real constraints.
- `huang_fit.py` on the matrices of the GaN run (FOCEX_DUMP_AMAT):
  - Translational rows plus the one real Huang relation [zz,xx]=[xx,zz]: FC2 changes by at most 3e-3 eV/A^2
    from the translational-only fit, and the fit goes from 1.202 to 1.214 %.
  - All 15 rows: FC2 changes by 1.0 eV/A^2, even with the forbidden rows' right-hand sides set to 0.
- The rotational rows have the same defect. GaN, transl+rot exact (1 1 0 1): 119 %, -38920 1/cm.
- The analytic models were never affected: their positions are symmetric to 1.6e-10 A.
- stress.params was not symmetric either (GaN sxy = 0.011 GPa, sxx != syy). That gave right-hand sides to rows
  whose coefficients are round-off (7e-9). The SVD dropped them, so this did no harm, but it was inconsistent.
- The cause is not the stress or the Coulomb part of the right-hand side:
  - the right-hand side of the real row is 0.12 eV, while the fit without Huang violates that row by 6.9 eV;
  - with symmetric positions the Huang rows are consistent with the other rows to 2e-7 eV.

## Fix

- `symmetrize_positions` (force_constant_symmetries.f90), called after make_unitcell, followed by symmetry_init
  again. The positions are replaced by their group average.
  - The fractional translations are re-estimated from all atoms, not from atom 1 only.
  - It is repeated until the change is below 1e-13 A, and the centroid is kept.
  - The log line SYMMETRIZE_POSITIONS gives the deviation, the shift and the positions to put back in
    structure.params.
  - The displacements are FORCEDISP - POSCAR (read_force_position_data), so the force data are unchanged.
- `pointgroup_average2`: stress.params is averaged over the point group (Huang rows and read_stress_params).
- Z* and eps are averaged over the group (check_born_symmetry, ios.f90; 2026-11-03).
- `relax_reference`: the net force on the cell is removed in pi_prime.
  - With enforce_inv 0 the translational rows are not exact; GaN FC1 summed to -2e-5 eV/A.
  - With atom 1 held fixed, this net force left the two Ga 1.05e-6 A apart.
  - Now: soft fit 3.5e-14 A; exact fit 1.4e-16 A (5e-10 before).

## Tests

**1. Exact low-symmetry model.** HfO2 (monoclinic, 12 atoms, stressed reference), ranks 1-3, enforce_inv 1 with
transl+rot+Huang.
- Inputs: those of the earlier Coulomb tests (coulomb/README.md), i.e. `c_base` (no charges) and `r3_on` (point
  charges, all right-hand sides).
- "Perturbed": the positions in structure.params only, by up to +-5e-5 A (`perturb.py ... 5e-5 7`). POSCAR and
  FORCEDISP are unchanged, so the exact FCs are unchanged.
- Errors are relative to the exact FCs, from coulomb/compare.py.

| run | positions | FOCEX | fit | FC2 rel | FC3 rel |
|---|---|---|---|---|---|
| no charges | exact | new | 0.000 % | 6.6e-8 | 2.3e-7 |
| no charges | perturbed | 9.5 | 92.5 % | 0.91 | - |
| no charges | perturbed | new | 0.000 % | 1.7e-6 | 1.9e-5 |
| point charges | exact | new | 0.000 % | 6.6e-8 | 2.3e-7 |
| point charges | perturbed | 9.5 | 1429 % | 7.6 | - |
| point charges | perturbed | new | 0.000 % | 1.6e-6 | 2.0e-5 |

- With exact positions the FCs agree to the printed digits.
- Perturbed, new code: what remains is the part of the perturbation that respects the symmetry. It cannot be
  detected and is of order delta/a.
- Huang check (mech.dat, Huang minus stress): 1.4e-8 GPa in all new runs.

**2. DFT data**, new code against 9.5 (+ the Z* symmetrization of 2026-11-03), long-wave C (E = 0) in GPa.
- GaN, user's setting 1 1 1 0:
  - before: C11 -7081, C66 -7288, -210 1/cm;
  - after: 316.7, 96.7, 137.3, 297.9, 90.5, 110.0 (C11, C12, C13, C33, C44, C66), no imaginary modes.
- GaN 1 1 1 1: 1.94 %, stable.
- ZnO 1 1 0 1: 58.7 % and -2.2e5 1/cm, then 1.30 % and stable. enforce_inv 0, 1 and 2 agree within 4 GPa.
- The hand-symmetrized inputs (scratch hsym_*) gave the same results as the code.

**3. Is a re-expansion of the FCs about the symmetric reference needed?** No.
- Test: the GaN POSCAR supercells were shifted by the same delta (`sym_poscar.py`), so that u is measured from
  the symmetric reference.
- FC2 changed by 7e-4 eV/A^2 (4e-5 relative), FC3 by 4e-4 relative, C by at most 0.5 GPa.

**4. Cross-validation** (`cv.py`, dumps of the new code).
- Whole snapshots are held out (GaN 6 folds over 42, ZnO 7 over 21). The Python full-data fits reproduce
  FOCEX's numbers.
- Held-out relative force residual, % (exact = enforce_inv 1; soft = enforce_inv 0, in brackets):

| rows | GaN | ZnO |
|---|---|---|
| T | 1.256 | 0.804 |
| T + Huang | 1.267 | 0.809 |
| T + rot | 2.002 (1.659) | 1.451 (1.337) |
| T + rot + Huang | 2.005 (1.665) | 1.455 (1.339) |

- Huang is neutral for the force prediction. It constrains the q -> 0 brackets, which supercell forces hardly
  determine.
- The rotational rows nearly double the held-out error. With their Coulomb right-hand side set to 0:
  GaN 1.924, ZnO 1.486. Open.

**5. Regressions.** Exact ZnO + Coulomb, ranks 1-4: positions moved 9e-11 A; FCs equal to 8e-10 relative;
band frequencies identical; mech.dat differs only in the sign of zeros.

## Open

- **Rotational rows with born_flag 2 and anisotropic Z*/eps:** the right-hand sides are inconsistent at 2.5e-3
  eV of 0.22 (`max|C xpart - d|`). With Z* and eps replaced by their traces: 2.4e-7. The subtracted FC3 (trace
  model) is not the rotational partner of the tensor Born FC2 (coulomb/README.md).
- **Rotational rows raise the held-out error** (test 4).
- **The relaxed K obeys Huang less well than Phi:** 2e-2 against 4e-7 GPa in the all-exact fit. Not explained.
- **ZnO long-wave C33 is 45-49 GPa** in every setting.

## Scripts

All are scratch tools with GaN constants inside (c, the z values in sym_poscar.py).
- `huang_diag.py DIR [tr]`: Huang violation of a fit, by pair dz.
- `huang_rows.py DIR huang_diag.py`: row norms and right-hand sides, SVD of the constraint sets.
- `huang_fit.py DIR huang_diag.py`: fits with all or only the real Huang rows.
- `viol.py DIR`: invariance violations and C from a run.
- `cv.py DIR NxR[,NxR] NFOLD`: cross-validation; needs scipy (`~/myenv/bin/python`).
- `perturb.py SRC DST AMP SEED`: perturbed structure.params.
- `sym_poscar.py IN OUT`: POSCAR shifted by delta.
