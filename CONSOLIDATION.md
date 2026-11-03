# FOCEX / THERMACOND consolidation -- status (branch `consolidate`, not committed)

Goal (2026-10-22): one copy of the routines the two codes share. Done in this worktree
(`~/PROJECTS/Anharmonic-lattice-dynamics-consolidate`, branch `consolidate`); nothing is
committed. Each step was checked with the regression harness before the next.

## Layout

`common/` holds the shared sources; both Makefiles compile them (`vpath %.f90 ../common`),
each into its own object folder with its own flags.

| file | content | from |
| --- | --- | --- |
| `constants.f90` `params.f90` `geometry.f90` `ios_mod.f90` `atoms_force_constants.f90` `lattice.f90` `svd_stuff.f90` `linalgb.f90` | the modules of `FOCEX/fc_modules.f90`, split one per file, unchanged | FOCEX |
| `tools.f90` | `sort`, `make_reciprocal_lattice_2pi`, `inverse_real`, `get_upper_bounds`, `delta_k`, `my_erfc` | FOCEX extratools (+ THERMACOND's `sort` loop) |
| `born.f90` | module `born` (epsilon, Born charges, NA arrays) | FOCEX ewald_born |
| `ewald.f90` | module `ewald`: Ewald grids and the dynamical-matrix sums | FOCEX ewald_born |
| `nonanal.f90` | `non_anal`, `non_analq0`, `dyn_coulomb_pure_step`, `set_dynamical_NA_3N3N`, `naq` | FOCEX aux_na |
| `dielectric.f90` | `read_dielectric` (epsilon, inverse, Born charges with the acoustic sum rule) | FOCEX ios |

FOCEX keeps its Ewald *forces* in `FOCEX/ewald_born.f90`, now module `ewald_forces`.

THERMACOND uses from `common/`: constants, params, geometry, ios, born, tools, ewald, nonanal,
dielectric. Its own modules keep the THERMACOND-only variables (`tc_params`, `tc_ios`); its
`atoms_force_constants`, `lattice`, `kpoints`, `eigen`, `om_dos`, `tetrahedron`, `svd_stuff`
are still its own (see "Next").

## Steps and checks

Harness: `tests/regression/run.sh SRC LABEL` builds both codes from a tree and runs
FOCEX on Ge, BZO, MoTe2/TEST, TEST_GeSe and THERMACOND on Ge 6^3 job 2, Ge 4^3 job 3 and BZO
job 0 (serial and `mpirun -np 3`); `compare.py A B` compares every output file. Inputs and
runs are in `~/PROJECTS/consolidate-ref`. `STRICT=1` compiles with `-ffp-contract=off` (see
step 6b). Two FOCEX outputs print uninitialised memory and are masked: `gr` on the log line
"Initial Pressure(GPa), gr,xagne" and the `atomposindx` column of `supercell_N.xyz` (never set).

| step | change | result |
| --- | --- | --- |
| 1 | FOCEX built from `common/` (fc_modules.f90 split) | all files identical |
| 2 | THERMACOND: shared `constants` (same values) | identical |
| 3 | shared `params`; THERMACOND-only inputs in `tc_params`; `use params, only :` lists | identical (log: include_fc printed 1:4) |
| 4 | shared `geometry` (THERMACOND-only routines were unused) | identical |
| 5 | shared `ios`; THERMACOND unit numbers that differ from FOCEX's (`udos`, `uibz`, `ugrun`) kept in `tc_ios` | results identical; `phonlog.dat` number formats of `write_out` |
| 6a | FOCEX: `tools`, `born`, `ewald`, `nonanal`, `dielectric` moved to `common/` | FOCEX identical (all 4 materials) |
| 6b | THERMACOND: FOCEX's dielectric reader and non-analytical term | Ge: results identical with `-ffp-contract=off` (ulp-level otherwise, fused multiply-add); BZO: now FOCEX's phonons (below) |

### Step 6b: THERMACOND phonons of polar materials

Before: THERMACOND added its own Gaussian-damped NA term (rho = |G|/4, no ASR on the charges)
to FOCEX's `fc2.dat`. For `born_flag = 2` (Ewald forces subtracted before the fit, the
setting recommended for polar materials) its BZO phonons differed from FOCEX's by up to
280 1/cm over the whole zone, so any kappa of a polar material was wrong.

Now THERMACOND calls FOCEX's `read_dielectric`, sets `coef`, the Ewald grids and the q=0 term
as `focex.f90` does (`setup_nonanalytic`), and adds `set_dynamical_NA_3N3N` and its
q-derivative (`add_nonanalytic`). BZO, FOCEX band path: |omega| agrees to <= 1.05 1/cm away
from Gamma (median 0.09); FOCEX evaluates every point shifted by 3e-4|G|, and v.dk at the
worst points is 0.4-0.7 1/cm. At Gamma itself FOCEX's shifted point gives a direction-dependent
LO-TO splitting, THERMACOND's naq the angular average. `born_flag = 6` needs FOCEX's supercell
grids: THERMACOND stops with a message.

### Also found

- FOCEX `sort` made n-1 selection passes, leaving the n-th entry unsorted when the array is
  padded (`make_grid_shell_ewald`). The shared `sort` makes n passes; FOCEX outputs unchanged
  on the 4 test materials.
- FOCEX's IBZ (`getstarp`, reduced coordinates) agrees with spglib for Ge, BZO, MoTe2 and GeSe
  (counts and star sizes); the IBZ bugs fixed in THERMACOND on 2026-10-22 were THERMACOND's.
- THERMACOND's `get_freq` took `abs(omega^2)` and both codes' hermitisation skipped the last
  column (`do j=t+1,n-1`): fixed on `main` (2026-10-22: omega^2 < -25 1/cm^2 stops, whole upper
  triangle); the merge applies the same loop to `add_nonanalytic`.

## Rebased onto FOCEX 9.1 (2026-10-29)

The worktree state on 2711efd (backup: `~/PROJECTS/consolidate-ref/backup_worktree_2026-10-29/`)
was merged with `main` at 9f8ad81 and is now uncommitted on that commit.
- `common/` modules re-split from `main`'s `fc_modules.f90` at the same boundaries (their
  concatenation equals `main`'s file); this brings `fcflag` (atoms_force_constants) and `bhom`
  (svd_stuff, `maxterms` gone).
- FOCEX Makefile: `main`'s new sources (fc34_range, fc34_maps, invariance_ranks, lasso_*, ...).
- `read_input_from_latfc` (new on `main`): `use tc_ios, only : ufc0`.
- Harness: binaries `focex-9.1.x`, `tcond-2.0.x`, `tcond-2.0-serial.x` built with `DESTDIR`;
  BZO and MoTe2 inputs replaced by the `~/PROJECTS/MATERIALS` copies (the old ones fail
  FOCEX 9.0's Born and POSCAR checks; kept in `consolidate-ref/inputs_old_2026-10-22/`).

Checks, `STRICT=1`, against 9f8ad81: FOCEX Ge, BZO, MoTe2, GeSe all files identical;
THERMACOND Ge 6^3 job 2 and 4^3 job 3, serial and 3 ranks, all result files identical.
`phonlog.dat`: dropped `RECIPROCAL_LATTICE:` headers, the added `read_dielectric` /
`SETUP_NONANALYTIC` echo, timing; BZO's unstable Gamma mode omega^2 -112360 -> -3449 1/cm^2
(step 6b). Reading of latdyn/structure/kpbs/dielectric.params, lat_fc.dat, fcN.dat: same code
as `main` apart from `use` lines; unit numbers of the two `ios` modules that coincide (31, 81,
93) are never both used in THERMACOND.

## Next (not started)

- `atoms_force_constants` and `lattice`: FOCEX's `lattice` contains `make_unitcell`, which runs
  FOCEX's whole setup (`symmetry_init`, its `force_constants_init`, supercell grids). Sharing it
  means THERMACOND setting up its lattice and symmetry the FOCEX way, from `structure.params`
  and `lat_fc.dat` without POSCARs: a design decision.
- `kpoints` / IBZ: done 2026-11-03 (user): one routine, common/ibz.f90 (`ibz_reduce`, `laue_operations`), THERMACOND's
  algorithm. FOCEX's `getstarp` version used the lattice point group: wrong for crystals of lower symmetry than
  their lattice (CdI2 P-3m1, 8x8x8: 50 IBZ points instead of 65); the agreement with spglib noted above holds
  only because Ge, BZO, MoTe2, GeSe have the full symmetry of their lattices.
- dynamical matrix (short-range part), `eigen`, `om_dos`, `tetrahedron`, `svd_stuff`.
- Step 7 (started 2026-10-30): phonons from one routine (`get_freq`, `get_frequencies`,
  `freq_vel`, `gruneisen`), so THERMACOND takes FOCEX 9.1's exact-k / Gamma-limit rule.
  `gruneisen` done 2026-11-01: `common/gruneisen.f90` (module `gruneisen_rules`) in both Makefiles;
  FOCEX's and THERMACOND's `gruneisen` are wrappers that pass their FC2/FC3 terms and dynamical
  matrix. THERMACOND's own routine removed. Checked on the exact HfO2 model, not with the harness
  (user: no Ge harness rerun); see CLAUDE.md "Mode Gruneisen parameters".
