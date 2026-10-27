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
cd FOCEX && make       # ~3 min; produces ~/ALADYN/bin/focex-9.0.x
```

Executables (2026-10-23, user's rule): all in `~/ALADYN/bin` (`DESTDIR`; the top-level `make`
passes make.inc's), named `<name>-<version>.x` (user: hyphens for all): `focex-9.0.x` (9.0 since 2026-10-23, user: invariance
rows for all ranks + the Born check; 9.1 is reserved for reading some FC ranks and fitting the rest;
8.18 and earlier: `~/BIN/v8.18`),
THERMACOND `tcond-2.0.x` / `tcond-2.0-serial.x` (was `kap8` in its
source folder). Top-level `make` (= `make all`, user's spec 2026-10-23): python-tools, focex (FOCEX +
FOCEX/utility), sc_snaps, thermacond, anfomod; single targets of those names, `thermacond-serial`,
`scop8` (not in all: SCOP8 does not compile -- `VA_math.f95:7416` calls `dotproduct(...)`, which does
not exist, geometry.f95 has only the `.dot.` operator; from 07f9f456, YuanLiang21120, 2025-04-02;
user: ignore SCOP8 for now). Verified in a scratch copy: make all exit 0, make clean exit 0.
Old utility builds were kept in `~/ALADYN/bin/old_2026-10-23/`. Tags (2026-10-23, approved): focex-8.19; the agent is
versioned separately, focex-agent-1.0, focex-agent-1.1 (thermacond task; `__version__` in focex-agent.py and diagnose.py, manual title and
its Version history section); 1.2 (2026-10-27, snap_000 ISIF=2 -> stress.params; no tag); VASP-Flow v2.6 (at ef92326) and v2.7.
`FOCEX/Makefile` does **not** include the top-level `make.inc`: it hardcodes
`FF=gfortran`, `FLAGS=-O2 -fcheck=all`, `exe=focex-$(ver).x`, `DESTDIR ?= ~/ALADYN/bin`. It has
CRLF line endings. The utility programs are separate
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

## focex-agent.py (FOCEX/AGENT/, since 2026-10-21; moved there 2026-10-23)

POSCAR -> VASP relax/SCF/DFPT -> SC_SNAPS -> supercell VASP runs -> FORCEDISP1 -> FOCEX ->
plots, restartable (`state.json`, `report.log`). Everything is in `FOCEX/AGENT/`: the script,
`diagnose.py`, `focexagent_manual.tex/.pdf` (same content as `docs/source/runfocexagent.rst`: change
both, rebuild both, Sphinx at zero warnings), `focex-agent-instructions.txt` (the template),
`tests/` (argon, fake VASP). `FOCEX/focex-agent.py` is a forwarder for old run scripts. The script
finds FOCEX/ (templates, utility/) as the parent of its folder (realpath, so symlinks work);
`make python-tools` (part of the top-level `make`) links them as `focex-agent`, `focex-diagnose` in DESTDIR.
Run it with `~/dftphonon-env/bin/python` (numpy, spglib, pymatgen; Homebrew python3 has none). That one
environment serves VASP-Flow too (made 2026-10-23 with the 50 packages of the old `~/vaspenv`, which the
ZnO/GaN/BAs run scripts still use and which is kept until the user removes it).

- Command line (2026-10-23, user's spec): `focex-agent.py [--instr FILE] [--coord FILE]
  [--profile NAME] [--force]` sets up; defaults `focex-agent-instructions.txt` (not
  instructions.txt: VASP-Flow uses that name) and `POSCAR`. Old `init -i/-s` still accepted.
  `PROJECT_DIR:` (relative to the instruction file; default the Project name) is where the run
  goes; the user then submits `run_workflow.sh` from there (both modes). `run`/`plots` subcommands
  are what the scripts call. The project's copy is `focex-agent-instructions.txt` (old projects:
  `instructions.txt`, still read).
- THERMACOND task (2026-10-23, user's design): last task `thermacond`, run iff THERMACOND_JOB
  (latdyn line 8; else line 8 of a latdyn_params block; default 0 = skipped) is 1/2/3. Copies
  lat_fc/fc2/fc3/structure/dielectric/kpbs from 05_focex to 06_thermacond, writes latdyn lines 1, 6,
  8, 9 from THERMACOND_MESH / _T / _JOB / _SOLVER, runs THERMACOND_X (newest ~/ALADYN/bin/tcond-*.x,
  MPI build) with the profile's MPI command; report: <= 7 temperatures, kappa_RTA, full, kappa_C,
  warning if full < 0.99 RTA (open direct-solver issue). diagnose.py: TC_PAT. Tested on argon (6^3,
  RTA isotropic, ~1/T; job 0 skipped; missing binary held).
- Reference stress (2026-10-27): snap_000 (undistorted supercell) runs with ISIF=2 (others ISIF=1 =
  trace only, vasp.at wiki); collect writes 04_supercells/stress.params = -(OUTCAR "in kB", XX YY ZZ
  XY YZ ZX)/10 GPa (user: VASP positive = cell wants to expand; wiki agrees; PSTRESS not subtracted,
  report warns if nonzero), copied to every FOCEX folder. fake_vasp.py prints the LJ virial stress
  (= -dE/dV to 6 digits). Argon test: stress 0, fit 0.51 %, max |w - exact| 0.2 1/cm.
- Snapshot temperature (user's rule, 2026-10-23): a `TEMPERATURE:` given is used as is (effective
  phonons at that T); only without it is the default 30 K retuned to `DISP_TARGET` (Ar: 50 K kept,
  rms 0.0138 A). `SUPERCELL:` 9 integers = supercell.inp rows, RS_i = sum_j n_ij R0_j (primitive).
- SC_SNAPS (2026-10-23, user's diagnosis): its spring model (central springs, stiffness ~1/r^4, rcut =
  2 x NN) kept only the nearest periodic image of each pair, so a 4-atom fcc cube had 3 of 12 NN
  springs, 6 zero modes, displacements ~1e6 A and `***` in poscar_NNN. get_normal_modes now sums all
  images within rcut: 3 zero modes, X: 6 TA + 3 LA, LA/TA 1.252 equal in the 4-atom and 32-atom
  cubes. Ordinary supercells change too (Ar 2x2x2 up to 1.05 THz: 2nd-neighbour images at L/2 tie):
  only the sampling of the displacements, not FOCEX's force constants. Guard kept: modes below 1e-3 x
  max are not sampled (warning); the agent stops on `***` in a snapshot.
- Modes (user's labels): `MODE: from_scratch | extract_only` (extract_only needs `DATA:`; DATA
  without MODE implies extract_only).
- `diagnose.py <PROJECT_DIR>...`: read-only; task status (done/FAILED/INTERRUPTED/not reached,
  running via squeue or ps, not started); VASP causes from VASP-Flow's `vasp_diagnose.py`
  (`error_lines`, `match_pattern`, `CAUSE_TEXT`, split out of its `classify` 2026-10-23; VASP-Flow
  diagnosis on 40 materials identical before/after), FOCEX stops by message (`FOCEX_PAT`).
  Writes diagnosis/report.txt + fix_it.sh (guard: nothing running, state.json/report.log
  unchanged; backup; INCAR tags into the failed runs' INCARs and the instruction block;
  resubmit). Held for review: unknown cause, fix already in the INCAR, new supercell. Tested:
  ZnO copy (ws_outside, held), argon with two broken snapshots (fixed and finished).
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
  outside the supercell's WS cell) -- a FOCEX issue still open. Same stop on ZnO 4x3x2 (2026-10-23):
  setup_FC2_in_supercell keeps a group if one member is inside the WS cell, set_dynmat_std stops
  if any member is outside; a supercell without the lattice's symmetry splits groups (the 8.3 A
  group: +-2a1 inside, +-2a2 outside). The agent had picked 4x3x2 over 3x3x2 by 8e-13 A in face
  distance; choose_supercell now rounds it to 1e-4 A and then takes the fewest atoms.
- Born charges are transposed from VASP's (field, displacement) to FOCEX's
  `charge(displacement, field)` (`aux_na.f90` contracts q with the 2nd index). The VASP
  convention is from memory, not checked against VASP docs: `BORN_TRANSPOSE` switches it.
- DATA mode (2026-10-23): `DATA: <folder>` (POSCARi/FORCEDISPi, optional primitive POSCAR and
  dielectric.params in the data frame; `SPECIES:` if no element line) -> `init_data`: lattice
  symmetrized (`ideal_frame`: spglib std_lattice parameters, same basis/order; MgO-Safoura data
  7.6e-7 A off), Kabsch rotation into the FOCEX frame applied to POSCARi, FORCEDISPi (pos+force),
  eps/Z (R^T T R); tasks relax..collect marked provided. Verified on MgO-Safoura: NOBORN fit
  7.197 %, 35 unknowns, 647 1/cm, B 184.0, G 117.3 GPa = direct run in the original frame.
- Crystal frame (2026-10-23, `crystal_frame`/`symmetrize_tensor4`/`report_crystal_frame` in
  dynamical.f90): symmetry axes from op_matrix(isgop) put along x,y,z (cubic 4-/2-folds;
  hex/trig/tetr high axis on z, 2-fold or shortest lattice vector on x; ortho three 2-folds;
  mono 2-fold on y); C averaged over the point group then rotated; Zener ratio and the
  100/110/111 sound speeds use it (they used structure.params axes: MgO rhombohedral orientation
  labelled crystal [110] as "100"). Distinct constants: Ge/BZO/MgO 3, MoTe2 6 (5 + C66), GeSe 9.
- Normal mode writes structure.params from the relaxed CONTCAR unsymmetrized (~1e-7): same risk as
  above (lost FC groups at 1.6e-8 non-orthogonality); not yet idealized.
- Test it without VASP: point `VASP_STD` at a fake that writes OUTCAR/CONTCAR/EIGENVAL with
  the markers `outcar_status` looks for (General timing, EDIFF is reached, reached required
  accuracy, BORN EFFECTIVE CHARGES).

## THERMACOND (consolidated with FOCEX, 2026-10-21)

Manual: `docs/source/runthermacond.rst`. Theory: npj Comput. Mater. 11, 303 (2025).

- Build: `cd THERMACOND && make` (mpif90 -> ~/ALADYN/bin/tcond-2.0.x) or `make serial` (gfortran +
  `mpi_stub.f90` -> tcond-2.0-serial.x); objects in obj_mpi/ obj_serial/; `DEBUG=1` for -g
  -fbounds-check (the old default, ~20x slower). On this Mac mpif90 needs
  `OMPI_FC=/opt/homebrew/bin/gfortran SDKROOT=$(xcrun --show-sdk-path)`.
- Run: `tcond-2.0.x` (serial) or `mpirun -np N tcond-2.0.x`; same code. Ranks > 0 work in `.mpi_rank_NNN/`
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
    as in calculate_RTA_distributionfunction: the FWHM in angular units (= 2 Sigma''), omega in
    rad/s -- verified 2026-10-23: the s=s' terms give C v^2/Gamma = RTA (kappa_P/RTA 1.01, Ge 8^3),
    and at X (Ge 8^3, 300 K) RTA 1/tau vs job-1 2 Im Sigma = 0.69/0.22/0.75 per mode (Lorentzian
    eta = 10 vs tetrahedra; no uniform 2 or 2pi). kappa_C not yet compared with Simoncelli or an
    independent code (user: later).
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
- structure.params optional (2026-10-23): `read_input_from_latfc` (rest.f90) takes the cell,
  species, masses, positions from lat_fc.dat, whose primitive-vector lines now carry the
  conventional vectors in columns 4-6 (FOCEX write_lat_fc; list-directed readers -- ANFOMOD,
  THERMACOND -- take columns 1-3). Ge 4^3 job 3 byte-identical with/without structure.params.
- job 1 is slow: function_self_w3 recomputes V3 over the mesh for each (k, mode); use few
  points per segment in kpbs.params.

## PLANNED RELEASE (user, 2026-10-23, very important): read some FC ranks, fit the others

include_fc(rank)=2 (structure.params line 4): fcN.dat of that rank is input, the include_fc=1
ranks are fitted. Force part: the read rank's columns times its known FCs go to the right-hand
side (A_fit x_fit = b - A_read x_read) and those columns leave amat and the unknowns. Invariance
part (built rank by rank: translational within a rank, rotational couples rank n and n+1, Huang
rank 2): the same column partition makes the rows coupling a read and a fitted rank inhomogeneous,
C_fit x_fit = -C_read x_read; enforce_inv=1/2 then needs x = x_p + K y with a particular solution.
Existing pieces: focex.f90 opens fcN_irr.dat and calls read_fcs (readfc), nindepfc counts only
include_fc=1 ranks; ios.f90 ~136 "fix this to be also able to read only some files and fit the
rest". To check: RHS subtraction in read_snaps_set_aforce; the invariance rows. Discuss the design
with the user before coding.

## Invariance rows: one loop over ranks (2026-10-23, invariance_ranks.f90; default)

`set_translational_inv_constraints` / `set_rotational_inv_constraints` (module invariance_ranks)
build the rows of every included rank (up to 8) in one pass over the group terms (row key = prefix
atoms + directions [+ nu]; one column function, offset + current(), skipping dropped FC2 terms).
The rank-by-rank routines are kept as `set_*_inv_constraints_old` (setup_amatrix.f90), selected by
`FOCEX_INVARIANCE=old`; `FOCEX_INVARIANCE_COMPARE=1` builds both, writes inv_rows_old/new.dat and
inv_fcs.dat, fits with the chosen one. include_fc=2 is not handled yet (marked: RHS b).
- Bug in the old rotational routine, fixed 2026-10-23 (user's diagnosis: the cnt3 offsets must
  use current()): rank-3 rows, FC2 part used res1+cnt3+ti with cnt3 counting dropped FC2 groups
  (fc2range=0 drops groups outside the supercell WS cell; groups are ordered by first atom, so
  this hits crystals with >= 2 inequivalent atoms when rank 3 is fitted): now keep check +
  res1+current(2,g,ti). Rank-2 rows: FC2 part had no keep check (a dropped term went to the
  column of the previous kept one); FC1 part summed rank-2 group sizes/keep flags: now
  current(1,g,ti).
- Verified after the fix, old vs new rows (compare mode): exact model below, 6 lattices + trigonal
  in the standard frame: same rotational row space, residuals 1e-16, exact fits. BZO with FC3
  (2 shells): same row space, fit 25.986 % (e0) / 25.151 % (e1) for both (before the fix the old
  rows differed: rank 42 each, 47 stacked). DFT regression (Ge, BZO, MoTe2, GeSe, MgO, e1): new =
  before to all digits, old-fixed within 0.1 1/cm. Translational: new contains old (old skipped
  a2<a1 etc.: MoTe2 33->43, GeSe 123->220, MgO 6->7 independent rows), fits unchanged. Argon agent
  test: max |w - exact| 0.236 1/cm.
- Exact test (pair potential A(rc^2-r^2)^2: FCs of ranks 1-4 exact; 7 crystal systems, Bravais
  cell doubled with 2 species, Huang off because the model is stressed): new rows hold on the exact
  FCs to round-off and the constrained fit reproduces the forces to 1e-11 eV/A; the old rotational
  rows before the fix: residual 3e-2..2e-1, constrained force error up to 0.6 eV/A.
- FOCEX with the "keep the DFT frame" structure.params (line 1 `1 1 1 90 90 90`, Cartesian line 2)
  fails when a symmetry operation is not an integer matrix in that cube basis: doubled trigonal
  cell (C2/m, 2-fold off the axes): unconstrained exact data fit to 4e-6 only, FC4 wrong (5.5 vs
  0); the same data in FOCEX's standard frame (agent extract_only) fit exactly. Cause not located
  (probably op matrices assumed integer in the conventional basis). Until fixed, use the standard
  frame (poscar2focex.x / the agent) for such crystals.

## Range of FC3..FC8: fc3flag..fc8flag (2026-10-24, user's design; fc34_range.f90, fc34_maps.f90)

structure.params line 12 = `fc2range [fc3flag ... fc8flag]`; absent flags = 1 = nshells only, as before
(maps.dat byte-identical to HEAD on BZO, MoTe2, MgO, Ge, GeSe). fcflag(rank)=0: a cluster is kept iff
every pair of it is a kept FC2 pair (its FC2 group has a member inside the supercell WS cell, the
test of setup_FC2_in_supercell / exclude_beyond_sc; table pair_keep built in setup_maps from
map(2) before rank 3) AND (all pairs within nshells(rank), or 2 distinct atoms, or collinear).
Filtered in the enumeration, so map(rank>=3) holds only kept groups (complete orbits) and needs no
keep/current (user asked; told). User (2026-10-26): keep the original routines unchanged; the new
version has its own names: setup_maps_fc34, size_force_constants_fc34, collect_force_constants_fc34
(fc34_maps.f90, copies + the rule); focex.f90 calls setup_maps_fc34 iff some fcflag(3:8)=0.
others.f90 = HEAD except bhom (Exact-model full test below). write_clusters -> clusters<rank>.dat (edges,
nd, type, group, tau; fcflag=0 route only). Checked: 2-body FC3 distances = kept FC2 distances
exactly (BZO, MoTe2, Ge, GeSe). Size (ranks 3,4 indep): BZO 77,25 -> 534,1722;
Ge 95,2 -> 329,862; GeSe 654,18 -> 4434,11083; MoTe2 76,7 -> 169,590; setup_maps BZO 14 s -> 8 min.
Exact-model test of the maps (2026-10-26, scratch exact/: model.py gen.py check.py): degree-4
polynomial potential (pairs to rc, triplets/quadruplets with all pairs < Rm; term lists fixed at the
reference, so FCs of ranks 1-4 exact at any amplitude, 0.03 A), forces on the materials' own
supercells; check.py compares every fcN_irr.dat term and 2000 fcN.dat terms per rank with the exact
FC (polynomial fit on a 5^m grid) and evaluates inv_rows_*.dat (FOCEX_INVARIANCE_COMPARE=1) on
them. fcflag=0, FC2 range inside the WS cell (fc2range=1; with fc2range=0 a 2-body cluster and its
image are both unknowns and only their sum is fixed): Ge, BZO, GeSe (rot rows off) fit 0.000 %,
all FCs = exact to the 7 written digits. Pre-existing failures found (original route too):
- rank-1 rotational rows (sum_i Phi_i x X_i = 0 over the primitive cell) are violated when the
  atoms feel forces at the reference (GeSe, MoTe2 Te z): only those 2 of 9957 rows fail on the
  exact FCs, the rest hold to 1e-9; fit 29 % (fc34) / 40.6 % (original route). Fixed 2026-10-26:
  rows dropped (user); see Exact-model full test.
- MoTe2 FC2 group of the in-plane Mo-Mo pair (3.499 A): maps.dat sets d2U/dy1dx13 =
  d2U/dx1dy13 (symmetric xy block); the exact model, invariant under all 24 ops of P6_3/mmc
  (spglib) to 2e-6, has xy=-13.0, yx=-5.37. Translational rows fail on the exact FCs; fit 2.8 %
  on the original route, 1.9 % fc34, rot rows off. FOCEX finds the 24 ops; error in the term
  analysis. Cause found 2026-10-26 on wurtzite (yz=zy): Exact-model full test, fix 1; MoTe2 not
  rechecked with an exact model since. (MoTe2 test input also: b=1.7320508 -> ops orthogonal to 6.5e-9 only;
  POSCAR1 cell differs from structure.params by 7.5e-4 A.)
- MgO (rhombohedral frame): HEAD with no flags and nshells(3)=6 6 also stops "Error in
  force_constants: zero equation"; with fcflag=0 rank 3 has groups with 385 independent terms.
  structure.params vs POSCAR1 consistent (cell 4e-6 A, sites 8e-5 A, species right), symmetry
  ops orthogonal to 4e-16. Cause not located. (Largest-pivot xrowop2, in the code since 2026-10-26, removed the stop
  but not the wrong groups.)

## How the force constants are solved for (enforce_inv)

`structure.params` line 5, 4th value:

- **0** — invariance rows stacked with the force rows, SVD on the whole system.
- **1** — constrained fit inside the kernel of the invariance matrix
  (`fit_kernel_basis`, `constrained_fit.f90`). `x = K y`, so the invariances are
  exact by construction and `K` *is* the elimination of the dependent FCs.
  Until 2026-08-07 this instead projected the unconstrained solution onto the
  kernel, which is **not** a constrained fit and gave ~3400 % residual on MoTe2
  where the constrained fit gives 35 %. Results predating `83b709b` differ.
- **1** with a stressed reference: x = x_p + K y (Huang rows with b /= 0, `particular_solution`).
- **2** — constrained LASSO on the irreducible FCs themselves (`lasso_constrained.f90`, called from
  `lasso_kernel_basis` in constrained_fit.f90): min 1/2N||Ax-b||^2 + lam sum w_j|x_j| (w = column
  rms), invariances exact, interior point; support |w x| > 1e-8 max; lam by 5-fold CV. **3** = 2 +
  least-squares refit on the support under the invariances. Old LASSO on the kernel coefficients:
  FOCEX_LASSO=old. History: memory lasso-work-log, FOCEX/notes/constrained_lasso.tex. The MoTe2/GeSe
  numbers once given here (88 kept, 789 of 795) were the old LASSO.

`svd_set` declares its matrix `intent(inout)` and `svdcmp` returns U in it, so
**it destroys the matrix you pass**. Give it a copy.

## Exact-model full test (2026-10-26/27; FOCEX/REFERENCE/exact_model, README.md there)

Model (user's plan): E = L.s + 1/2 s^T H s in squared distances (pair K, anharmonic pair, 3-body
area^2, 4-body opposite edges -- the face-normal product has no 4-atom part), stressed reference;
CaF2 (96 atoms), wurtzite ZnO (72), monoclinic HfO2 (96), 100 SC_SNAPS snapshots each, 1499/1927/7977
unknowns. verify.py checks map terms, completeness, invariance rows, force rows, fits and the support
by FOCEX's rule (the written x is never exactly sparse: non-zero counts cannot show a dropped FC).
All pass; e0/e1/e3 exact to the 7 written digits; e2 loses no FC, FC4 off <= 4e-4 (shrinkage).
Code changes it brought:
1. force_constant_symmetries.f90: unique_force_constant was called with iatom as input and output;
   when the canonical form reordered the atoms, the next direction combinations of the same
   operation used the wrong atoms, imposing Phi(i a, j b) = Phi(i b, j a) (wurtzite FC2 yz=zy;
   MoTe2 xy=yx same mechanism). Output now in iatomu.
2. xrowop2: partial pivoting (largest element of the column); same echelon form in exact arithmetic.
3. Rank-1 rotational rows (FC1 torque) dropped (user): image-dependent when FC1 /= 0.
4. Stressed Huang: S_abcd - S_cdab = -2 Omega (d_ab sig_cd - d_cd sig_ab); optional stress.params
   (3x3 GPa, sig = (1/V) dE/deta, tensile > 0, absent = 0, echoed in the log); rows with b /= 0 ->
   x = x_p + K y (e1) and z = z0 + Q t (LASSO); bhom in others.f90.
5. LASSO (lasso_constrained.f90): lam grid extended by 2 decades while the CV minimum is at its end,
   to 1e-8 lam_max (old floor 1e-4 lost 76 true FCs on ZnO e3); tie tolerance 1e-10 -> 1e-13 <b^2>;
   refit particular solution by pivoted QR (NR svdcmp did not converge in 18/242 calls). Each
   extension recomputes the CV (~2.3x cost; reusing the done lam not implemented).
6. Test dumps: FOCEX_DUMP_MAP=1 -> map_dump.bin; FOCEX_DUMP_AMAT=k -> amat_dump.bin (ahom, bhom,
   force rows of the first k snapshots).
DFT regression (Ge, BZO, MoTe2, GeSe, e1-e3; table in the README): Ge, BZO e1 = HEAD; MoTe2 e1
35.269 -> 32.194 %, GeSe 6.111 -> 5.755 % (changes 1-3, not separated); Ge e2 0.083 -> 0.036 %;
the rest identical. The log label "non-zeros ... after projection" (focex.f90) is left over from the
old e1 path: it counts |x| > margin.

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
- QHA (`QHA_free_energy_vs_strain`, 2026-10-23): P(V) is now -dF/dV of the same F(V): the static
  pressure P0 - B dV/V0 was missing from P (Ge: P=0 at +10 % strain, F minimum at +0.25 %), F
  lacked -P0 dV, and the phonon pressure uses the exact -dln(omega)/dlnV of the omega(V) used
  (gamma (1+e)/(1 +- 3 gamma e), over V not V0). Harmonic-only fits keep P_agne = NkT/V (user's
  term for a harmonic thermal expansion) with F_agne = -NkT ln(V/V0) in F, so it does expand
  (BZO 300 K dV/V 0.43 -> 0.56 %). temp_free.dat columns 10-15: F0,P0,F_agne,P_agne,F_ph,P_ph.
- Born tensors of equivalent atoms must sit on the right site (dielectric.params blocks follow
  the structure.params atom order). ~/PROJECTS/MATERIALS/BZO had the three O tensors cyclically
  permuted (O at (1/2,0,1/2) given -4.85 on xx; its Zr-O bond is along y): with born_flag=2 the fit
  was 25.4 %, imaginary Gamma modes (-131 1/cm), cubic degeneracy broken; with them on the bond
  directions 3.46 % (born_flag=0: 2.89 %; free linear fit of the same forces: 3.2 %). Since 9.0
  `check_born_symmetry` (ios.f90, after set_neighbor_list) tests Z(k') = R Z(k) R^T for every
  space-group op (site symmetry included) and eps = R eps R^T: > 5e-2 relative stops, and every
  reassignment of the blocks within a species is tried (BZO: block order 1 2 4 5 3, or atom order
  1 2 5 3 4 in structure.params -- the user's fix -- written to dielectric.params.symmetric);
  > 1e-4 warns (GeSe, TEST_GeSe: equivalent Ge with Z_yy 2.598 / 2.630, 1.0e-2). RbBr 1e-17.
- OPEN (BZO, 2026-10-23): with the correct Z and born_flag=2 the zone-boundary spectrum differs
  from born_flag=0 (R: top 916 vs 768 1/cm, above the Gamma LO 610) -- the NA term added back to
  the dynamical matrix may not equal the Ewald force subtracted from the data.
- `conductance.dat` is written but its band loop is `do j=1,-1`, so it is all
  zeros.
- Elastic constants (dynamical.f90, 2026-10-27; formulas in theory.rst "Elastic constants"):
  checks with WARNING banner + size (tol 1e-5): ASR of phi and pi0, sum_tau Q = 0, rank-2 rotational
  Q(tau a;b,g)-Q(tau a;g,b) = pi0_g d_ab - pi0_b d_ag (from teta, before Q is symmetrized), and
  generalized Huang [ab,gl]-[gl,ab] = d_ab s_gl - d_gl s_ab with s from stress.params (else 0, as in
  the Huang rows; tol 1e-3 max|bracket|). No torque test (sum_tau pi0 x tau depends on which images
  form the cell when pi0 /= 0; the rank-1 rotational fit rows were dropped for the same reason).
  Q is symmetrized in its strain indices before xi, sigmav, QGQ (C itself unchanged: symmetrize4
  already did it). stress.params (GPa, tensile > 0, same file as the Huang rows) replaces the FC1
  estimate sum_tau pi0 tau/V of the reference stress; sigma_res = sigma_ref - pi0 G Q/V;
  C = Sym(A-QGQ)/V - Sym(d_ag sigma_res_bd), all 4 terms (was sigma0 unrelaxed, 2 of 4 terms: shear
  correction doubled). Verified: Ge, BZO unchanged; Ge with stress.params -1 GPa: C11 +1, C44 +0.5,
  C12 same, Huang unaffected; deviatoric (0.5,0,-0.5) flags Huang by 1.0 GPa. With enforce_inv=0
  BZO/MoTe2 trip the rotational check at 1e-5, GeSe 1.6e-3 (and ASR 3e-4).
  Bravais lattices (mechanical0) take the same stress correction and Huang check (Ar -1 GPa: C11 +1,
  C44 +0.5). Helpers read_stress_params, check_huang_stress (others.f90 has its own check_huang).
  supercell_*.xyz: atom_sc(i)%cell%atomposindx (only written there) is assigned since 2026-10-27 in
  identify_atoms_in_supercell (index in atompos, 0 beyond its shells; verified: atompos(idx) = r_sc
  + constant shift for every atom of Ge, GeSe). theory.rst Eq. 8: 1/2 chi (user).
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

`docs/source/focexagent_guide.rst` (2026-10-27, user): short focex-agent guide with example inputs,
installation at its end; it and the full agent manual (runfocexagent.rst) come last, before the FAQ
(the agent was added after ALATDYN). Keep its keywords in step with the template and the full manual.

`docs/source/runfocex.rst` is the FOCEX manual — inputs line by line, output
columns, troubleshooting keyed to the code's actual error strings. Keep the
build at zero warnings.

## Related repos

- SC-SNAPS is the folder `SC_SNAPS/`, a **git subtree** of `github.com/KeivanS/SC_SNAPS`
  (added 2026-10-23 at 66267c1, full history; user's choice over a submodule so that a plain
  clone/pull brings it). Edit it only here; publish to the separate repository with
  `git subtree push --prefix=SC_SNAPS sc_snaps main` (remote `sc_snaps`; sandbox off). A change
  made in that repository directly comes in with `git subtree pull --prefix=SC_SNAPS sc_snaps
  main`. `~/PROJECTS/SC_SNAPS` is the old clone, no longer the working copy. Built by
  `make sc_snaps` (top level, BINDIR = DESTDIR with ~ expanded: its Makefile quotes it) ->
  ~/ALADYN/bin/sc_snaps.x; the agent looks there, then ~/BIN. With a Flask GUI. The stale copy in `FOCEX/utility/` was
  removed 2026-08-14; it had diverged from this one across ~2200 lines and read
  an older `cell.inp` (line 3 now takes two values, line 6 takes masses *and*
  charges). Documented in `docs/source/runscsnaps.rst`; that repo has its own
  `CLAUDE.md` with the details.

## Where the git clone is

As of 2026-10-15 **this directory is the git clone** (`main`, tracking
`origin/main`) — commit from here. The old clone under
`~/Documents/Documents - MAE-/GitHub/` is no longer the working copy.
Current focus is the **FOCEX** folder.

`origin` is **SSH** since 2026-10-27 (`git@github.com:KeivanS/Anharmonic-lattice-dynamics.git`):
`~/.ssh/id_ed25519` is in the agent and authenticates as KeivanS, so `git fetch`/`git push` need
no password and work **with the Bash sandbox on**. Over HTTPS they do not: the `osxkeychain`
helper has no github.com entry (no `.netrc`, no `GH_TOKEN` either), so a push only got
`could not read Username`, and HTTPS needed the sandbox disabled even to connect. The SC_SNAPS
clone at `~/PROJECTS/SC_SNAPS` is still on an HTTPS remote. `WebFetch` and sandboxed `curl`
cannot reach github.com at all. Other hosts (readthedocs.io, PyPI) are fine.
