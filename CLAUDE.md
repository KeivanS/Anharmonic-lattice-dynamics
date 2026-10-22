# ALATDYN — working notes

Fortran suite for anharmonic lattice dynamics: **FOCEX** (extracts force
constants from force–displacement data by SVD), **SCOP8**, **THERMACOND**,
**ANFOMOD**. Documentation is Sphinx, in `docs/`, published on Read the Docs.

## Ground rule

**Verify against the source and an actual run.** The published docs and the
shipped example outputs have drifted from the code more than once. Read the
`read`/`write` statements, then build and run. Build in a scratch directory —
the Makefiles drop `.o`/`.mod` files next to the sources.

## Building and running FOCEX

```bash
cd FOCEX && mkdir -p ~/BIN && make       # ~3 min; produces ~/BIN/v8.18
```

`FOCEX/Makefile` does **not** include the top-level `make.inc`: it hardcodes
`FF=gfortran`, `FLAGS=-O2 -fcheck=all`, names the binary `exe=v$(ver)` (`ver=8.18` -> `v8.18`), and
`mv`s it to `~/BIN`, which it does not create. The utility programs are separate
(`cd FOCEX/utility && make`, installs to `make.inc`'s `DESTDIR=~/ALADYN/bin`) —
that builds `read_outcar.x`, `read_qe.x`, `poscar2xyz.x` and `poscar2focex.x`.
`poscar2focex.x` (since 2026-10-22, standalone Fortran, no spglib): POSCAR -> structure.params,
POSCAR_focex, kpbs.params in the Setyawan-Curtarolo convention (pymatgen's analyzer.py
`get_conventional_standard_structure(international_monoclinic=False)` and kpath.py
`KPathSetyawanCurtarolo` were the reference). Crystal system from the crystal's own rotations;
rhombohedral primitive (2,1,1),(-1,1,1),(-1,-2,1)/3 (pymatgen's order is left-handed and mirrors
the crystal); atoms in a canonical setting (proper rotations of the cell, then origin, then input
order) so the output is independent of the input basis/orientation/origin. Verified: 26 test
structures (all S-C sub-cases but MCLC2/4) match pymatgen in sub-case, path, conventional and
primitive cells and reduced k; 78 randomly skewed/rotated/shifted inputs give identical files;
skewed rotated Ge -> FOCEX example gives byte-identical fc2/fc3/fc*_irr/lat_fc/mech.dat.
`FOCEX/utility/poscar2focex.py` (numpy + spglib, the earlier version; focex-agent.py imports it): POSCAR →
`structure.params` (Niggli-reduced, rotated into FOCEX's frame of the conventional
cell; 6-number line 1, line 2 = spglib transformation, atoms in conventional units, other lines = Ge
defaults marked "default, check") + `POSCAR_focex` (the frame of line 1). Same cell
as VASP-Flow's `modules/cell_reduce.py`. Verified 2026-10-16: generated from a
skewed, rotated Ge cell, the Ge example gives identical fc2/fc3 (78 / 173). `sc_snaps.x`
comes from the SC_SNAPS repo; the stale bundled copy was removed 2026-08-14.

FOCEX takes no arguments and no stdin. Run it in a directory holding:
`structure.params`, `dielectric.params`, `latdyn.params`, `kpbs.params`,
`default.params`, `POSCAR1`, `FORCEDISP1` (plus `POSCARi`/`FORCEDISPi` for
extra supercells).

## The Ge example — verified reference

`FOCEX/example/Ge` runs clean as of 2026-08-07 (216-atom supercell, 36
snapshots, ranks 1–3). Reproduces:

| quantity | value |
| --- | --- |
| irreducible FC2 / FC3 | 78 / 95 (173 unknowns) |
| transl / rot / Huang constraints | 14 / 84 / 15 |
| fit error `\|\|F_dft-F_fit\|\|/\|\|F_dft\|\|` | 0.036 % |
| C11, C12, C44 | 121.1, 68.0, 47.4 GPa |
| Bulk modulus (Hill), Debye T | 85.7 GPa, 313.6 K |

The `fc*.dat`, `log*.dat`, `svd-all.dat` etc. **stored** in that folder are from
an older version — different column layout. Don't use them as reference.

## Regression testing FOCEX

Build a reference binary from an untouched checkout and a second from the working
tree, run both on the same inputs, and compare **every** output file, not just
`fc*.dat`. Test inputs are in `~/PROJECTS/MATERIALS`; the three that span the
symmetry range are **BZO** (cubic perovskite, ~2 min), **MoTe2/TEST**
(hexagonal, deliberately short range, ~20 s) and **TEST_GeSe** (orthorhombic,
low symmetry, ~7 min). `NaCl` has no POSCAR/FORCEDISP and cannot run a fit.

Runs are reproducible as of 2026-08-07. Before the `finitedif_vel` fix they were
not, and any phonon-level comparison was meaningless — check that first if
outputs differ for no reason.

## focex-agent.py (FOCEX/, since 2026-10-21)

POSCAR -> VASP relax/SCF/DFPT -> SC_SNAPS -> supercell VASP runs -> FORCEDISP1 -> FOCEX ->
plots, restartable (`state.json`, `report.log`). Manual: `docs/source/runfocexagent.rst`;
example input: `FOCEX/focex-agent.example.txt`. Run it with `~/vaspenv/bin/python` (needs
numpy, spglib, pymatgen; Homebrew python3 has none of them).

- Imports VASP-Flow (`~/VASP-Flow` or `$VASP_FLOW_DIR`): `InstructionParser` for the VASP
  keywords, `VASPInputGenerator` for INCAR/KPOINTS/MAGMOM/GGA+U/KPAR-NCORE/k-path/DFPT,
  `vasp-agent.py` for POTCAR and profiles. Static INCARs are made from VASP-Flow's
  relaxation INCAR, so their tag lists go under `incar_raw['relax']`. Don't copy that code.
- Its own blocks (`INCAR_snaps ... INCAR_snaps_END`, `structure_params`, ...) are cut out
  before VASP-Flow's parser sees the file -- that parser rejects them otherwise.
- POSCAR0 is rotated so its spglib conventional cell is exactly FOCEX's/SC_SNAPS's frame
  from a,b,c,alpha,beta,gamma; `cell.inp` and `structure.params` share lines 1-2.
  Verified: Ge example gives byte-identical fc2/fc3/mech.dat to a plain FOCEX run.
- Supercells are multiples of the conventional cell. A skewed HNF supercell (Ge, 60 atoms)
  made FOCEX stop in `set_dynmat_std` with `R+tp-t not inside!` (default FC2 range reaches
  outside the supercell's WS cell) -- a FOCEX issue still open.
- Born charges are transposed from VASP's (field, displacement) to FOCEX's
  `charge(displacement, field)` (`aux_na.f90` contracts q with the 2nd index). The VASP
  convention is from memory, not checked against VASP docs: `BORN_TRANSPOSE` switches it.
- DATA mode (2026-10-23): `DATA: <folder>` (POSCARi/FORCEDISPi, optional primitive POSCAR and
  dielectric.params in the data frame; `SPECIES:` if no element line) -> `init_data`: lattice
  symmetrized (`ideal_frame`: spglib std_lattice parameters, same basis/order; MgO-Safoura data
  7.6e-7 A off), Kabsch rotation into the FOCEX frame applied to POSCARi, FORCEDISPi (pos+force),
  eps/Z (R^T T R); tasks relax..collect marked provided. Verified on MgO-Safoura: NOBORN fit
  7.197 %, 35 unknowns, 647 1/cm, B 184.0, G 117.3 GPa = direct run in the original frame.
- FOCEX's Zener ratio and mean sound speed / Debye T use tensor components as if the conventional
  axes were x,y,z: MgO in the rhombohedral orientation 806 K vs 982 K in the FOCEX frame (B, G
  invariant). Open: make vbar frame-invariant.
- Normal mode writes structure.params from the relaxed CONTCAR unsymmetrized (~1e-7): same risk as
  above (lost FC groups at 1.6e-8 non-orthogonality); not yet idealized.
- Test it without VASP: point `VASP_STD` at a fake that writes OUTCAR/CONTCAR/EIGENVAL with
  the markers `outcar_status` looks for (General timing, EDIFF is reached, reached required
  accuracy, BORN EFFECTIVE CHARGES).

## THERMACOND (consolidated with FOCEX, 2026-10-21)

Manual: `docs/source/runthermacond.rst`. Theory: npj Comput. Mater. 11, 303 (2025).

- Build: `cd THERMACOND && make` (mpif90 -> kap8) or `make serial` (gfortran +
  `mpi_stub.f90` -> kap8_serial); objects in obj_mpi/ obj_serial/; `DEBUG=1` for -g
  -fbounds-check (the old default, ~20x slower). On this Mac mpif90 needs
  `OMPI_FC=/opt/homebrew/bin/gfortran SDKROOT=$(xcrun --show-sdk-path)`.
- Run: `kap8` (serial) or `mpirun -np N kap8`; same code. Ranks > 0 work in `.mpi_rank_NNN/`
  (symlinked inputs, removed at the end) so only rank 0's files survive. Since 2026-10-22
  (degenerate averaging, below) Ge 6^3 job 2 and 4^3 job 3 give byte-identical kappa_vs_T.dat for
  1 and 3 ranks and for GCC 12.2 -O2, GCC 15.3 -O2 and -O0. Before that the result depended on the
  compiler (Ge 4^3 RTA 13.25 / 13.34 / 13.29 W/mK). `/usr/local/gfortran` (GCC 12) comes first on
  this Mac's PATH, so plain `make serial` uses a different compiler from mpif90.
- Shares structure.params and latdyn.params with FOCEX (13-line latdyn, FOCEX reads 1-7;
  template FOCEX/latdyn.params). job (line 8): 0 bands/DOS, 1 +linewidths on the band
  path, 2 +V3*delta, RTA, direct (line 9: 0 CG, 1 SVD), 3 +V3*delta, RTA, iterative.
  All nT temperatures in one run; V3*delta is T-independent: line 10 = 0 writes
  V3DELTA/k######.bin per IBZ point (header nc, ndyn, nibz, k) and reuses it.
- kappa.f90 was rewritten 2026-10-21 (old 16-branch job control incl. split-job mode is
  in git history before that). The physics routines in sub.f90 are unchanged apart
  from guards that let calculate_w3_ibz_split_sq / calculate_FGR_mpi2 be called once per
  IBZ point.
- OPEN: direct solution (job 2) disagrees with the iterative one (job 3) -- Ge 6^3, 300 K:
  RTA 29.52, iterative 31.21 (isotropic), direct 28.11/28.10/28.17 even with SVD (exact
  solve, residual 5e-14). So the IBZ collision matrix M = C*S (CollisionMatrix_mpi) or its
  RHS (cal_RHS_IBZ) is wrong; pre-existing. Also: cg_quadratic is plain CG (needs a
  symmetric matrix; M is not): residual 7e-3 after 7 iterations; its printed
  "|Ax-b|" is really funcmin. mode_kappa_d (old CG path) uses star-averaged velocities,
  mode_kappa (RTA, SVD, now all paths) does not; mode_kappa resets tauinv_eff inside its
  band loop.
- OPEN: isotope scattering (line 11) not wired into jobs 2/3 (stops with a message); the
  Gamma-point optical calculation (old job 50, or 51+52 with tetrahedra) not yet in job 0
  -- ask the user which variant.
- Physics fixes 2026-10-22 (Ge results above are after them):
  - velocities: Hellmann-Feynman (`get_freq_vg` in rest.f90; `get_freq` wraps it), incl. d/dq of
    the NA term in `nonanal`; finite differences only in `veltest.dat` (FOCEX layout). FD = HF to
    2e-5 median on non-degenerate modes, 5e-5 where the NA term changes v by up to 9 %.
  - |V3|^2 averaged over degenerate modes at q1, q2, q3 (`average_degenerate_v3`, kappa.f90;
    `degtol` = 1e-3 1/cm). Degenerate eigenvectors are compiler-dependent; Schur's lemma: each
    partner's rate is the trace / size. Ge symmetry-degenerate modes come out split by up to
    1.5e-5 1/cm, probably the 7-digit fc2.dat.
  - Wigner conductivity in the RTA (`wigner_kappa`, kappa_wigner.dat; Simoncelli, Marzari,
    Mauri 2019): velocity operator (dD)_ss'/(w_s+w_s') + atom-position-gauge term
    i(w_s^2-w_s'^2)<tau>. kappa_P (degenerate pairs, trace form) is basis-independent and
    isotropic for Ge; Ge kappa_P/RTA = 1.13, 1.03, 1.01 on 4^3, 6^3, 8^3. Gamma_s = 1/tau taken
    as in calculate_RTA_distributionfunction; its 2pi convention vs the angular W of the
    Lorentzian (matters for kappa_C only) is not verified.
  - distribution_F.dat: F per IBZ point and mode, RTA and final, every T.
  - IBZ: get_weights builds stars modulo G (one star had been split over several "IBZ" points:
    Ge fcc 4^3 11 -> 8, 6^3 26 -> 16) and takes as representative the WS-folded member with the
    largest u.k, u = (0.9, sqrt2-1, sqrt5-2): one wedge (cubic x>=y>=z>=0, hexagonal 0-30 deg).
    kstarsinfbz1 matches mod G and drops duplicate arms. getstar uses the Laue group (isgop +
    inversion), not the lattice point group, and its mod-G test uses k.r_i/2pi.
  - op_kmatrix (integer, conventional reciprocal basis) had been applied to Cartesian k, F and v
    in getstar, syop_f, syop_f2 and the collision-matrix routines: wrong for hexagonal,
    trigonal, monoclinic and triclinic cells (MoTe2 included) until 2026-10-22; now op_matrix.
    IBZ checked against spglib 2.7.0 (get_ir_reciprocal_mesh, time reversal) on 4x4x4 and 6x6x4:
    Ge 8/34, Imma test 18/42, wurtzite 12/21, P3m1 test 13/24, star sizes identical.
- Degenerate blocks (item (c), 2026-10-22): P1/P2 (|V3|^2 x tetrahedron weight) are also averaged
  over the sets at q2 and q3 (`average_degenerate_p`; sums, hence Q and RTA, unchanged), so
  F = tau (v + Delta) with Delta common to a set, and kappa_vs_T.dat reports per-mode value +
  kappa_P - per-mode RTA (RTA column = kappa_P). Test hook THERMACOND_ROTATE_DEGENERATE=<seed>
  (random unitary in each degenerate set, `rotate_degenerate` in rest.f90): Ge 4^3 job 3 block
  RTA and iterative identical to all digits for seeds 0/1/2; per-mode values moved 3 %.
- OPEN: that block iterative value is anisotropic, basis-independently: Ge 4^3 100 K
  41.035/41.027/41.126 (RTA exactly isotropic) -- kernel symmetry (tetrahedra orientation or the
  star handling of q2/q3 in update_mpi); check with the direct-vs-iterative item. The
  direct-vs-iterative numbers above predate the 2026-10-22 fixes.
- Checks (2026-10-22, rules agreed with the user): omega^2 < -25 1/cm^2 stops `get_freq_vg`
  (error stop), -25..0 -> warning and |omega^2| (same -25 as focex-agent); D -> (D+D^+)/2 in
  full (the old loop `j=t+1,ndim-1` missed the last column), warning if max|D-D^+|/max|D| > 1e-6
  (same in FOCEX `make_hermitian`); empty fc3.dat stops (READ_FC3) instead of infinite tau / NaN.
  The symmetrisation exposed the per-mode basis dependence below (Ge 4^3 RTA moved 0.5 % for a
  1e-16 change of D; kappa_P moved 1e-6).
- OPEN (seen, not fixed): THERMACOND's NA term (Gaussian-damped, rho=|G|/4, no exp(iq(ta-tb))
  phase) may not match FOCEX's NA treatment -- check on a polar material.
- k-mesh: the compiled `make_kp_reg` is the one in modules.f90 (kpoints.f90 is not in the
  Makefile): g1,g2,g3 = g01,g02,g03 = reciprocal of lat_fc.dat's primitive cell (read_lattice
  now stops if structure.params' r0i differ); even meshes Gamma-centred (shift -1/2), odd from
  Gamma; latdyn line-2 shift ignored. get_k_info assumes exactly that labelling (verified 4^3-7^3:
  k, k+G, -k, IBZ; counts = spglib 8/10/16/20; `check_mesh_labels` runs every time). Odd meshes
  stopped in make_kp_reg_tet ("is not a kpoint": corners shifted by -(g1+g2+g3)/2 unconditionally)
  until 2026-10-23.
- job 1 is slow: function_self_w3 recomputes V3 over the mesh for each (k, mode); use few
  points per segment in kpbs.params.

## How the force constants are solved for (enforce_inv)

`structure.params` line 5, 4th value:

- **0** — invariance rows stacked with the force rows, SVD on the whole system.
- **1** — constrained fit inside the kernel of the invariance matrix
  (`fit_kernel_basis`, `constrained_fit.f90`). `x = K y`, so the invariances are
  exact by construction and `K` *is* the elimination of the dependent FCs.
  Until 2026-08-07 this instead projected the unconstrained solution onto the
  kernel, which is **not** a constrained fit and gave ~3400 % residual on MoTe2
  where the constrained fit gives 35 %. Results predating `83b709b` differ.
- **2** — same, but the reduced problem is solved by LASSO (`lasso.f90`), with
  the L1 weight chosen by 5-fold cross-validation on held-out force rows. Worth
  it only when the model is underdetermined: on MoTe2 with 880 free parameters
  and 1152 equations it keeps 88 and stays physical while the SVD produces
  ±15000 cm⁻¹ phonons; on well-determined GeSe it keeps 789 of 795 and matches
  the SVD. Sparsity is in the free parameters `y`, not in `x = K y`.

`svd_set` declares its matrix `intent(inout)` and `svdcmp` returns U in it, so
**it destroys the matrix you pass**. Give it a copy.

## Things that bite

- `sc_snaps.x` built before 2026-10-22 aborts with `Recursive call to nonrecursive
  procedure 'findif'` when compiled with `-fcheck=all` (`d2v` called `findif(r,d1v)`, and
  `d1v` calls `findif`). Fixed in SC_SNAPS: `d2v` writes the same central difference out;
  identical freqs.dat/modes.dat, any flags work. Rebuilt copies in `~/BIN` and `~/ALADYN/bin`.
- The parameter files are read positionally by list-directed `read`. A missing
  line or a missing value silently swallows the next record and crashes later
  with a confusing message. `default.params` must have 36 lines; line 2 of
  `latdyn.params` must carry **six** numbers (shift *and* normal).
- `bands.dat` columns 6+ are ω², not ω (`bandgen.plt` applies `521.1*sqrt()`).
- Group velocities (bs_freq.dat, bands.dat, ibz_bands.dat) are Hellmann-Feynman since 2026-10-22:
  `finitedif_vel` returns `v1` from `get_freq`; before, it returned the finite difference of
  energy-sorted bands (BZO crossing: 0.17/8.50 km/s FD vs 0.65/9.08 HF). Both stay in
  `veltest.dat`. Regression on Ge, BZO, MoTe2, GeSe: only those velocity columns changed.
- The log line "Initial Pressure(GPa), gr,xagne" prints an uninitialised `gr` (~1e-314).
- fc*.dat carry ~7 significant digits; symmetry-equivalent FCs then differ at ~1e-7.
- `conductance.dat` is written but its band loop is `do j=1,-1`, so it is all
  zeros.
- `structure.params` line 3 is `scal [convcoord]`. It is read as text, so the
  old one-value form (`5.70 # scale`) still works and means `convcoord=0`.
  `convcoord=0`: atom positions are in units of the conventional vectors
  (line 1); otherwise in units of the primitive vectors (line 2, columns of
  `prim_to_conv`). Both are converted to Cartesian in `make_unitcell` before
  `unitcell`, which expects Cartesian input.
- `structure.params` line 1 is always the **6 numbers** `a b c alpha beta gamma`;
  FOCEX (and SC_SNAPS, THERMACOND) put a along x, b in the xy plane. A 9-number
  Cartesian line 1 existed 2026-10-15..21 (e284db1) and was removed at the user's
  request. POSCARi/FORCEDISPi/dielectric tensors are Cartesian and must be in that
  frame; a DFT cell in another orientation gives *nearly* integer "rs1 in units of
  r0i" and "supercell vector k of POSCAR ... is not integer" (e.g. AgCrP2S6_PARA, a
  along y). Keep the DFT frame with line 1 `1 1 1 90 90 90`, line 2 the Cartesian
  primitive vectors, line 3 `1.0 1`, atoms in primitive units; or let
  focex-agent.py / poscar2focex.x rotate the cell before the supercells are made.
- `POSCARi` takes **either** VASP layout since 2026-08-14: `read_supercell` in
  `ios.f90` reads the line after the third lattice vector as text and retries on
  the next one if it does not parse as `natom_type` integers, so the VASP-5
  element-name line is optional. `POSCAR1` is now `cp poscar_000 POSCAR1`; the
  old `sed '6d'` recipe still works and is still in old scripts.
- `FORCEDISPi` has two accepted block formats: "POSITION" (absolute positions
  in Å, forces in eV/Å — what `read_outcar.x`/`read_qe.x` write, preferred) and
  "vasprun" (displacements in Bohr, forces scaled by 27.2116). **Open question:**
  that factor is a Hartree, not a Rydberg, and carries no 1/a_B — flagged to the
  author, unresolved.
- `read_qe.x` takes its input filename on **stdin**: `echo pw.out | read_qe.x`.
- POSCARi species (2026-10-22): `check_input_poscar_consistency_new` first requires each atom
  to sit on a site of the species its counts line gives it (order of structure.params) --
  this fixes the origin where positions cannot (MgO sublattices); if impossible it falls back
  to species from positions with a WARNING (MoTe2 test POSCAR1: counts 64 128 but Te first).
  Both "cannot map" stops are `error stop` now. Cartesian POSCAR positions are multiplied by
  line 2 (VASP rule): MgO-ANH3-cubic's POSCAR1 (16.86 with positions in A) was invalid.
- Lattice precision: lines 1-2 with irrational ratios to 7-8 digits (1.7320508) make the
  Cartesian symmetry ops non-orthogonal (1.6e-8) and the FC symmetry analysis silently wrong
  (MgO in an orthorhombic cell: one FC3 group lost, all FC3 = 0; full precision fixes it;
  exact mechanism inside the term analysis not pinned down). symmetry_init warns above 1e-10.
- Band order (bs_freq.dat, bands.dat, bs_grun.dat; `band_sort_bs`, 2026-10-22): Hungarian
  assignment on |<e(k-1)|e(k)>|^2 after aligning degenerate subspaces (polar decomposition),
  5 x mean level spacing window as before. The old greedy version left unconnected bands at
  their initial index: GeSe had a copied row at 179 of 180 k-points. Ge/BZO/MoTe2 unchanged.
- `setup_maps`: `collect_force_constants` counts a candidate group before dropping it when all
  its terms vanish, so its test is `ngroups.gt.maxgroups` with maxgroups = exact + 1 (the
  `.ge.` test stopped MgO rank 3 with "sizing disagreement").

## Docs

```bash
pip install -r docs/requirements.txt && cd docs && make html
```

`docs/source/runfocex.rst` is the FOCEX manual — inputs line by line, output
columns, troubleshooting keyed to the code's actual error strings. Keep the
build at zero warnings.

## Related repos

- `github.com/KeivanS/SC_SNAPS`, cloned at `~/PROJECTS/SC_SNAPS` — the **only**
  home of SC-SNAPS, with a Flask GUI. The stale copy in `FOCEX/utility/` was
  removed 2026-08-14; it had diverged from this one across ~2200 lines and read
  an older `cell.inp` (line 3 now takes two values, line 6 takes masses *and*
  charges). Documented in `docs/source/runscsnaps.rst`; that repo has its own
  `CLAUDE.md` with the details.

## Where the git clone is

As of 2026-10-15 **this directory is the git clone** (`main`, tracking
`origin/main`) — commit from here. The old clone under
`~/Documents/Documents - MAE-/GitHub/` is no longer the working copy.
Current focus is the **FOCEX** folder.

`git` over HTTPS to GitHub works **only with the Bash sandbox disabled**;
`WebFetch` and sandboxed `curl` cannot reach github.com at all. Other hosts
(readthedocs.io, PyPI) are fine.
