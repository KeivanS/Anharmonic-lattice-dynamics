# FOCEX full test on an exact model (2026-10-26/27)

Reference for reproducing the test of FOCEX 9.0 (working tree of 2026-10-27, LASSO fixes included)
against a model whose force constants of ranks 1–4 are known exactly. Plan by K. Esfarjani; scripts
and runs by Claude. Every number below was copied from the log named next to it.

## The model (`model3.py`, `model.py`)

E = L·s + ½ sᵀHs in the squared distances s of the atoms. This is a polynomial of degree 4 in the
positions, so FC1–FC4 are exact and FC5 and higher vanish. The term lists are fixed at the
reference structure, so no cutoff is crossed at any amplitude.

| term | form | range |
|---|---|---|
| harmonic pair | K s | r < rL |
| anharmonic pair | A (s − s0)² | r < Ra |
| 3-body | B \|AB × AC\|² (area squared) | all pairs < Ra |
| 4-body | C [(AB·CD)² + (AC·BD)² + (AD·BC)²] (opposite edges) | all pairs < Ra |

- Each anharmonic term has its part linear in (s − s0) removed. So the stress and FC1 come only
  from the harmonic term, and the reference is **stressed**: no relaxation, FC1 ≠ 0 where the site
  symmetry allows it.
- PAR: K0 0.05, A 0.03, B 0.02, C 0.006.
- The planned 4-body term was the face-normal product (AB × AC)·(AB × AD). Symmetrised over the
  tetrahedron it has no 4-atom part, because the face normals of a tetrahedron sum to zero. The
  opposite-edge form replaces it.

## Structures and data (`setup.py MAT run_MAT`)

The cells come from `poscar2focex.x` (`MAT/`, `p2f.log`). Supercells and 100 snapshots come from
SC_SNAPS; T was tuned to an rms displacement of 0.04 Å (`run_MAT.setup.log`). The forces are the
model's, written at full precision (FORCEDISP1).

| structure | space group | supercell | atoms | rL / Ra (Å) | irreducible FCs (FC1–4) |
|---|---|---|---|---|---|
| CaF2 fluorite | Fm-3m | 2×2×2 conv. | 96 | 5.0 / 4.2 | 0 / 23 / 127 / 1349 = 1499 |
| ZnO wurtzite | P6₃mc | 3×3×2 | 72 | 4.72 / 3.5 | 2 / 68 / 292 / 1565 = 1927 |
| HfO2 | P2₁/c | 2×2×2 | 96 | 3.87 / 3.22 | 9 / 294 / 1263 / 6411 = 7977 |

The stress of the stressed reference is in `run_MAT/setup.json`. It goes to FOCEX as
`stress.params` (GPa; `mkrun.sh` writes it).

## Checks of the model itself

- **`fdtest.py`** (`fdtest.log`):
  - Analytic forces and Hessian against finite differences: within 3e-8, which is the FD error.
  - FC2, FC3 and FC4 at the reference against differences of the analytic Hessian: within 9e-12.
- **`huang.py`** (`huang.log`): the stressed Huang relation
  S_abcd − S_cdab = −2Ω(δ_ab σ_cd − δ_cd σ_ab) holds on the exact FC2 to 2e-13.
  - In ZnO and HfO2 the right-hand side is 0.23 and 0.43. With the opposite sign it would fail by
    0.45 and 0.86.

## FOCEX checks (`verify.py run_MAT/eE`, log in `run_MAT/eE/verify.log`)

FOCEX dumps its maps and matrices (`FOCEX_DUMP_MAP=1`, `FOCEX_DUMP_AMAT=k`). verify.py checks them
against the exact FCs, with x_exact being the exact value of every irreducible FC:

1. **Map:** every term of every group satisfies Σ mat·x_exact = the exact FC.
2. **Completeness:** every model cluster with an atom in the primitive cell is in a group.
3. **Invariance rows:** C x_exact = b_hom, with translational, rotational and stressed-Huang rows.
4. **Force rows:** A x_exact = F.
5. **Fit:** FOCEX's fcN_irr.dat and fcN.dat against the exact values.
6. **Support:** the FCs that FOCEX's own rule selects (|w x| > 1e-8 max, with w the rms of the
   column) against the exact non-zero set.
   - The written x is never exactly sparse, so counting non-zeros cannot detect a dropped FC.
     Section 6 can.

### Map, invariance and force rows (identical in every run of a structure; worst values)

| check | CaF2 | ZnO | HfO2 |
|---|---|---|---|
| map rank 1 / 2 / 3 / 4, max error | — / 2e-13 / 6e-15 / 4e-16 | 8e-15 / 3e-10 / 8e-10 / 9e-11 | 3e-15 / 4e-15 / 9e-15 / 2e-16 |
| rank-4 terms checked | 4,172,541 | 746,854 | 849,690 |
| wrong terms, clusters missing | 0, 0 | 0, 0 | 0, 0 |
| translational rows | 457: 2.4e-13 | 2276: 5.5e-10 | 4356: 1e-14 |
| rotational rows | 3370: 3.1e-13 | 9645: 2.7e-9 | 25677: 5.3e-14 |
| Huang rows (max \|b\|) | 15: 2.4e-15 (0) | 15: 7.2e-9 (0.227) | 15: 5.6e-9 (0.43) |
| force rows, max \|A x − F\| (eV/Å) | 864: 5.8e-14 | 648: 7.4e-11; 21600: 1.1e-10 | 864: 1.2e-14 |

- ZnO's map, translational, rotational and force residuals (up to 3e-9) are larger than those of
  CaF2 and HfO2 (≤ 3e-13). So are the Huang rows of the two stressed crystals, ZnO and HfO2 (7e-9
  and 6e-9, against b of 0.23 and 0.43).
  - Why is not known.
  - ZnO's cell is exact: a = 3.2496 and γ = 120 in structure.params, as `poscar2focex.x` writes
    them from `ZnO/POSCAR`.
- `e1all` used all 100 snapshots (21600 force rows dumped). The other runs dumped 3 snapshots, but
  every run fits all 100.

### Fits (max |FOCEX − exact| / max |FC| per rank; fcN_irr.dat carries 7 digits)

| run | enforce_inv | FC1 | FC2 | FC3 | FC4 | support (true non-zero: selected, lost) |
|---|---|---|---|---|---|---|
| CaF2 e0 | 0 SVD | — | 1.6e-7 | 1.2e-7 | 9e-14 | all selected, 0 lost |
| CaF2 e1 | 1 kernel | — | 1.6e-7 | 1.2e-7 | 2e-13 | all, 0 lost |
| CaF2 e2n | 2 LASSO | — | 3e-7 | 5.6e-6 | **8.5e-5** | 18/61/207 selected (true 18/61/190), 0 lost |
| CaF2 e3n | 3 LASSO+refit | — | 1.6e-7 | 1.2e-7 | 6e-12 | 18/61/190 = exactly the true set |
| ZnO e0 | 0 | 8e-8 | 9e-8 | 1.6e-7 | 1e-9 | 0 lost |
| ZnO e1all | 1 | 8e-8 | 9e-8 | 1.6e-7 | 9e-10 | 0 lost |
| ZnO e2n | 2 | 6e-6 | 7e-7 | 1e-5 | **3.6e-4** | 2/47/243/579 (true 2/47/229/553), 0 lost |
| ZnO e3n | 3 | 8e-8 | 9e-8 | 1.6e-7 | 8e-10 | 2/47/229/553 = exactly the true set |
| HfO2 e0n | 0 | pending | | | | |
| HfO2 e1 | 1 | 1e-7 | 6.5e-8 | 2.2e-7 | 5e-10 | 9/210/861/1380 = true, 0 lost |
| HfO2 e2n | 2 | pending | | | | |
| HfO2 e3n | 3 | pending | | | | |

- The force fit is 0.000 % in every run.
- enforce_inv=0, 1 and 3 reproduce the exact FCs to the digits written.
- enforce_inv=2 (LASSO without the refit) loses no true FC. Its errors (≤ 3.6e-4 relative, FC4) are
  the L1 shrinkage, which the refit removes.

Run time and memory, from `time.err`:

| run | time (h) | peak memory (GB) |
|---|---|---|
| CaF2 e0, e1 | 1.4, 1.3 | 3.6 |
| CaF2 e2n, e3n | 2.3, 1.9 | 3.6 |
| ZnO e0, e1all | 0.7, 0.6 | 2.3 |
| ZnO e2n, e3n | 3.2 each | — |
| HfO2 e1 | 9.2 | 9.2 |

## LASSO: what was wrong and what was changed (`lasso_constrained.f90`)

- **The λ grid stopped too early.** With exact data the CV error falls all the way to λ → 0, but
  the grid ended at 1e-4 λ_max.
  - With that floor ZnO enforce_inv=3 lost 76 true FCs (FC4 error 1.4); in the standalone test on
    3 snapshots the floor lost 108.
  - Fix: while the CV minimum is at the end of the grid, the grid is extended by two decades at a
    time, down to 1e-8 λ_max. The log says when this happens.
  - Each extension recomputes the CV from the start, about 2.3 times the cost. Reusing the λ
    values already done would remove that; it is not done.
- **The tie rule was too loose.** It takes the largest λ within tol·⟨b²⟩ of the CV minimum. With
  tol = 1e-10 it took a λ where 3 true FCs were still off the support; tol is now 1e-13.
- **The refit's particular solution (`shift_part`) used the wrong solver.** NR `svdcmp` of B Bᵀ did
  not converge in 18 of 242 calls on the full ZnO data (eigenvalues clustered at 0 and 1), and it
  carried on with the unconverged result. It is replaced by a Householder QR with column pivoting
  of the off-support rows.
- **Standalone check** (`lasso/`, ZnO, all 21600 rows, grid floor 1e-4 then extended):
  - With the refit: support = exactly the 831 true FCs; errors 1.8e-13 / 1.8e-10 / 6.4e-10 / 1.1e-8
    (FC1–4).
  - Without the refit: 874 selected, 0 lost, FC4 error 4.3e-3.
- **How to read FOCEX's LASSO log lines:**
  - "force constants selected" is the LASSO support at the chosen λ, before the refit. For ZnO
    enforce_inv=3 it is 1199. The refit then sets the 368 truly-zero ones in it to ≤ 2.7e-9.
  - "non-zeros … after projection" counts |x| above the absolute write-out margin. "Projection" is
    a leftover label from the old enforce_inv=1 path.
  - "interior point did not converge at some λ" (CaF2 e2n, info 2) came with duality gaps of
    1e-15 to 1e-18.

### Conclusion on LASSO

- With accurate forces and the refit (enforce_inv=3), it found exactly the true support in all
  tested structures, and the FCs to round-off.
- Without the refit (enforce_inv=2), it loses no true FC, but the FCs carry shrinkage of order
  1e-4–1e-3 (relative, FC4).
- What noisy forces do to the support selection is not measured here. The 1 % noise test of
  2026-10-24 (cubic model, before these fixes) gave an FC4 max error of 206 for SVD, 9.2 for
  enforce_inv=2 and 4.8 for enforce_inv=3, of 288.

## Code changes the test found or needed (FOCEX working tree, not yet committed)

1. **The term analysis forced Φ_ab = Φ_ba.** `unique_force_constant` was called with the same
   array, `iatom`, as input and output (force_constant_symmetries.f90). When the canonical form
   reordered the atoms, the later direction combinations of that operation used the wrong atoms.
   - The output now goes to a separate array, `iatomu`.
   - Found here on wurtzite (FC2 yz = zy).
   - The MoTe2 xy = yx block noted in CLAUDE.md fits the same mechanism. MoTe2 has not been
     rechecked against an exact model since the fix.
2. **`xrowop2` had no partial pivoting** (force_constant_symmetries.f90).
3. **The rank-1 rotational rows (FC1 torque) are dropped.** On a stressed reference they depend on
   which periodic image is used (invariance_ranks.f90 and the old routine).
4. **The Huang rows are generalised to a stressed reference**, with σ from the optional
   `stress.params` (0 if absent, echoed in the log).
   - The inhomogeneous rows are solved as x = x_p + K y (enforce_inv=1) and z = z0 + Q t (LASSO).
   - Files: setup_amatrix.f90, constrained_fit.f90, lasso_constrained.f90, focex.f90.
5. **The LASSO fixes above.**
6. **Test dumps:** `FOCEX_DUMP_MAP=1` (map_dump.bin) and `FOCEX_DUMP_AMAT=k` (amat_dump.bin,
   force rows of the first k snapshots).

## DFT regression (enforce_inv = 1, 2, 3; ~/PROJECTS/MATERIALS inputs)

- `head`: binary built from git HEAD 3be3a83.
- `old`: the working tree before the LASSO fixes.
- `new`: the working tree with the LASSO fixes.
- Phonon range = min..max of bs_freq.dat column 7 (1/cm).
- Script: `regression/run.sh MAT E TAG BIN`.

| material | run | fit (%) | phonons | LASSO kept, info | grid extensions |
|---|---|---|---|---|---|
| Ge | head e1 | 0.036 | 0.08..296.90 | | |
| Ge | new e1 | 0.036 | 0.08..296.90 | | |
| Ge | old e2 | 0.083 | 0.08..297.10 | 115 of 173, 3 | 0 |
| Ge | new e2 | 0.036 | 0.08..296.90 | 166 of 173, 3 | 2 |
| Ge | old/new e3 | 0.044 | 0.08..296.90 | 107 of 173, 0 | 0 |
| BZO | head/new e1 | 3.287 | −52.75..690.90 | | |
| BZO | old/new e2 | 3.824 | 0.08..902.90 | 101 of 213, 0 | 0 |
| BZO | old/new e3 | 3.306 | −74.03..909.80 | 145 of 213, 0 | 0 |
| MoTe2 | head e1 | 35.269 | −89.60..283.50 | | |
| MoTe2 | new e1 | 32.194 | −81.96..293.60 | | |
| MoTe2 | old e2 | 32.194 | −81.93..293.60 | 153 of 153, 3 | 0 |
| MoTe2 | new e2 | 32.194 | −81.96..293.60 | 153 of 153, 3 | 2 |
| MoTe2 | old/new e3 | 32.212 | −81.82..293.60 | 137 of 153, 0 | 0 |
| GeSe | head e1 | 6.111 | −0.37..223.40 | | |
| GeSe | new e1 | 5.755 | −0.36..223.40 | | |
| GeSe | old/new e2 | 5.779 | −0.37..222.90 | 1640 of 1718, 0 | 0 |
| GeSe | old/new e3 | 5.756 | −0.36..223.40 | 1640 of 1718, 0 | 0 |

- "old/new" means the outputs are identical (all files compared). For Ge and BZO, head e1 and new
  e1 are identical.
- MoTe2 and GeSe e1 differ from HEAD. Two changes act there: items 1 and 3 above (the MoTe2 xy
  block; the rank-1 rotational rows that fail when atoms feel forces, as in GeSe and MoTe2 Te z).
  Item 2 (pivoting) may also act. Which change gives which part of the difference has not been
  isolated.
- BZO's three methods give different spectra. Which one is right is not established.

## Reproducing

```bash
cd FOCEX && make                     # build FOCEX from the commit that holds this folder
cd REFERENCE/exact_model
python setup.py CaF2 run_CaF2        # only to regenerate data: SC_SNAPS samples anew, so the
                                     # stored run_MAT/FORCEDISP1 will not be reproduced bit for bit
python fdtest.py run_CaF2 ; python huang.py run_CaF2
FOCEX_X=~/ALADYN/bin/focex-9.0.x ./run.sh CaF2 3 1 n   # -> run_CaF2/e3n (refuses an existing dir)
python verify.py run_CaF2/e3n
```

- Python: `~/dftphonon-env/bin/python` (numpy, scipy).
- Standalone LASSO, in a scratch directory:
  ```bash
  python lasso/mkin.py <run>             # writes lasso_in.bin from the run's amat_dump.bin
  gfortran -O2 lasso/svd_extract.f90 ../../lasso_constrained.f90 lasso/drv.f90 -o drv.x
  ./drv.x 1e-4 41 1                      # epsratio, nλ, refit
  python lasso/ev.py ...
  ```
- B₀′ test of FOCEX/notes/generalized_huang.tex, Eq. (Bp), from this folder: `python bprime/bptest2.py run_ZnO -0.005`. It takes the run folder and D, the coefficient of an extra pair term D(s−s₀)³. The output for all six cases is in `bprime/bp2_all.log`; formula and full relaxation agree to 10⁻⁶ in B′.
- Not stored here: the dumps (amat_dump.bin and map_dump.bin, up to 2 GB), maps.dat, fc3.dat and
  fc4.dat, and the phonon outputs. They are regenerated by the runs.
