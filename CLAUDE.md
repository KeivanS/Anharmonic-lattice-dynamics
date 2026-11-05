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
cd FOCEX && make       # ~3 min; produces ~/ALADYN/bin/focex-9.5.x
```

Executables (2026-10-23, user's rule): all in `~/ALADYN/bin` (`DESTDIR`; the top-level `make`
passes make.inc's), named `<name>-<version>.x` (user: hyphens for all): `focex-9.5.x` (9.5 2026-11-03: include_fc=2, born_flag 2 hybrid +
Coulomb RHS, elastic_longwave, IBZ from the crystal's Laue group (common/ibz.f90), 8-digit fc files; 9.4 2026-11-02: IR/Raman activity and irreps at Gamma, computed character tables, Gamma threshold, allplots; 9.3 2026-11-02: relaxation of an unrelaxed reference, default.params line 5; 9.2 2026-11-01: elastic constants, B0', QHA, Gruneisen,
sd -- see below; 9.1 since 2026-10-29, user: frequencies at the
exact k, Gamma limit, no stop for split FC2 groups; 9.0 2026-10-23: invariance rows for all ranks + the
Born check; reading some FC ranks and fitting the rest got 9.5 (2026-11-03; 9.3 went to the relaxation, 9.4 to the Gamma symmetry, user 2026-11-02);
8.18 and earlier: `~/BIN/v8.18`),
THERMACOND `tcond-2.0.x` / `tcond-2.0-serial.x` (was `kap8` in its
source folder). Top-level `make` (= `make all`, user's spec 2026-10-23; thermacond-serial added 2026-11-01, user):
python-tools, focex (FOCEX +
FOCEX/utility), sc_snaps, thermacond, thermacond-serial, anfomod; single targets of those names,
`scop8` (not in all: SCOP8 does not compile -- `VA_math.f95:7416` calls `dotproduct(...)`, which does
not exist, geometry.f95 has only the `.dot.` operator; from 07f9f456, YuanLiang21120, 2025-04-02;
user: ignore SCOP8 for now). Verified in a scratch copy: make all exit 0, make clean exit 0.
Old utility builds were kept in `~/ALADYN/bin/old_2026-10-23/`. Tags (2026-10-23, approved): focex-8.19; the agent is
versioned separately, focex-agent-1.0, focex-agent-1.1 (thermacond task; `__version__` in focex-agent.py and diagnose.py, manual title and
its Version history section); 1.2 (2026-10-27, snap_000 ISIF=2 -> stress.params; no tag); VASP-Flow v2.6 (at ef92326) and v2.7.
`FOCEX/Makefile` does **not** include the top-level `make.inc`: it hardcodes
`FF=gfortran`, `FLAGS=-O3 -fopenmp` (since 2026-11-04; `make DEBUG=1`: `-O2 -fcheck=all -fbacktrace -fopenmp`,
`make clean` when switching), `exe=focex-$(ver).x`, `DESTDIR ?= ~/ALADYN/bin`
(a test build in a copied tree still moves its binary there: use `make DESTDIR=<scratch>`; 2026-11-03 a test
build replaced the installed focex-9.5.x). It has
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

## Speed: OpenMP, -O3, large files (2026-11-04, user)

- Released as FOCEX 9.7 (tag focex-9.7, 2026-11-04, user approved); the "After 9.6" items are in 9.7.

- Rule (user): OpenMP loops must give the same numbers as the serial code. Every parallel loop is over
  independent outputs (each element computed by one iteration, sums kept in their old order), checked
  byte for byte with `cmp` on the outputs. A change that only reorders sums is reported with its size.
- Threads: `common/omp_threads.f90` (FOCEX and THERMACOND). `FOCEX_NUM_THREADS` / `TCOND_NUM_THREADS` if
  set; 1 under an MPI launcher (OMPI_COMM_WORLD_SIZE, PMI_SIZE, MPI_LOCALNRANKS; THERMACOND also nprocs > 1);
  `SLURM_CPUS_PER_TASK`; else cores - max(1,cores/5), i.e. leave ~20 % for other work (user: 8 of the Mac's
  10). OMP_NUM_THREADS is not used (the user's ~/.bashrc sets it to 1 for MPI codes). Printed in log/times.dat.
- FOCEX svdcmp (svd_ridge.f90): right Householder step by column blocks + OpenMP over columns, bit-identical
  (scratch driver 3000x1200, 6000x2400; HfO2r all outputs). Batching the QR rotations was tried: slower, dropped.
  `kernel_and_particular`: one SVD gives the kernel and xpart (particular_solution did a second SVD of the
  same matrix). U, V dumps (elimination.dat, matrix_svd-*.dat; 2.1 GB + 572 MB for HfO2r) only with verbose.
- -O3 vs -O2 -fcheck=all: FCs, frequencies, thermo identical (4 regression cases, HfO2r); mech.dat/log
  differ only in quantities that are 0 (1e-9). HfO2r relax: 7843 s -> 6520 s (-O3) -> 1596 s (8 threads).
  Left: building the rotational rows (735 s), translational rows (247 s), force rows, all serial.
- THERMACOND (`-fopenmp`, -O2 kept): threads in calculate_w3_ibz_split_sq (q2 loop), calculate_FGR_mpi2
  ((jj,kk); module array `tet` is `!$omp threadprivate`, copyin), CollisionMatrix_mpi (n2 loop, block offsets
  in closed form d=3((n-1)ndn+l-1)+1), update_mpi ((l1,xyz)). Ge j2_6, j3_4 serial and MPI -np 3: every
  output file byte-identical to the old build. 8 threads: j2_6 83.2 -> 12.5 s, j3_4 27.6 -> 6.3 s.
  CG (cg.f90, dense matmul on root) left serial: ~2 ms here ("direct solution" in times.dat).
- THERMACOND prints `reading 0 <garbage>` (modules.f90:853): nkdir, ndir printed before they are read. Harmless.

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

`FOCEX/example/Ge/reference_output/` (2026-10-27) holds the outputs of focex-9.0 built from
ccba45b (fc*.dat, lat_fc.dat, mech.dat, bs_freq.dat, ibz_dos.dat, thermo_QHA.dat, log); they
match the table. structure.params line 12 there is `0 1 1 1 1 1 1` (outputs identical to `0`).
User's rule for the manual: show the working path only, no history or failure lists; fix the
code/scripts instead (collection scripts in FOCEX/utility now start FORCEDISP1 afresh).

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
`tests/` (argon, fake VASP). `FOCEX/focex-agent.py` was a forwarder to it for the run scripts of
pre-2026-10-23 projects; removed 2026-10-30 (user: one copy only). A project whose
`run_workflow.sh` still names that path needs the path changed to `FOCEX/AGENT/focex-agent.py`
(or a fresh setup). The script
finds FOCEX/ (templates, utility/) as the parent of its folder (realpath, so symlinks work);
`make python-tools` (part of the top-level `make`) links them as `focex-agent`, `focex-diagnose` in DESTDIR.
Run it with `~/myenv/bin/python` (numpy, spglib, pymatgen; Homebrew python3 has none). That one
environment serves VASP-Flow too. History: `~/vaspenv` until 2026-10-23, `~/dftphonon-env` until 2026-11-03,
`~/myenv` since (user). Checked 2026-11-03: a new venv (`python3.14 -m venv`, Homebrew Python 3.14.4, made 12:22),
ase 3.29.0, flask 3.1.3, phonopy 3.5.0, spglib 2.7.0, numpy 2.4.4, pymatgen; `~/vaspenv`, `~/venvs`,
`~/dftphonon-env` and `~/anaconda3` are gone. Scripts, manuals, VASP-Flow use `~/myenv`.

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
- SLURM account (2026-10-27, user): the account is set once in the user's cluster profile
  (copy of VASP-Flow's profiles/slurm.json; rivanna.json removed from VASP-Flow); ACCOUNT: in the
  instruction file overrides it; neither -> init stops before making the project (slurm_account in
  VASP-Flow). Docs/templates use `--profile mycluster` as the placeholder, no site names.
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
- More snapshots (2026-10-29, user): SNAPSHOTS raised after `snaps` is done -> `add_snapshots` runs
  sc_snaps.x in 03_snaps/added_NNN-MMM (same cell/supercell/snaps.inp, count = difference; new
  clock seed) and copies poscar_NNN..MMM; old ones untouched; supercells..plots redone through the
  missing-folder check, only new folders run VASP. Lower SNAPSHOTS: nothing removed (report says use
  unused/). Argon 20 -> 30: old poscars/OUTCARs byte-identical, 10 VASP runs, FORCEDISP1 31, rms
  0.0203/0.0205 A, max |w - exact| 0.224 1/cm.
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
- Supercells are multiples of the conventional cell. A supercell without the lattice's symmetry
  splits FC2 groups (ZnO 4x3x2, 8.3 A group: +-2a1 inside the WS cell, +-2a2 outside);
  setup_FC2_in_supercell keeps a group if one member is inside. Until 2026-10-29
  set_dynamical_matrix_std / fourier_fc2fit_4D then stopped with `R+tp-t not inside!`; the stops
  were removed (user) after an exact harmonic test (central springs to 6.1 A on ZnO, forces summed
  over all images; scratch only): 4x3x2 with split groups gives every FC2 exact and all bands
  within print rounding; the symmetric 3x3x2 (same shortest vector 9.75 A) fits the forces
  exactly but its groups at 5.0-6.0 A lie wholly on the WS boundary, so only FC(d)+FC(d+L) is
  fixed: FCs off by 0.1, bands by 11.8 1/cm, exact only at the commensurate k. Members outside
  the WS cell are logged once (log_fc2_outside_ws, others.f90). The agent had picked 4x3x2 over
  3x3x2 by 8e-13 A in face distance; choose_supercell rounds it to 1e-4 A and takes the fewest atoms.
- Phonon frequencies (2026-10-29, user): freq_vel (dynamical.f90) gives omega at the exact k,
  eigenvectors and HF velocities at k + 3e-4 |G| uinc (degeneracies); Gamma = limit q = 1e-6 along
  the path segment (band path, get_frequencies ipath=1) or uinc (mesh); the 3 rigid translations
  there are set to 0 (log FREQ_VEL), gruneisen gives them gamma 0, QHA skips them.
  Gamma test |k| < gamma_tol = 1e-5 (phonon_rules) in both gamma_limit and freq_vel (2026-11-02, user:
  LaP path ended at |k| = 2.7e-10: freq_vel zeroed the translations, gamma_limit (then 1e-12) did not see
  Gamma -> long-wave gamma along uinc, 0.207/-0.558/1.580 next to 0.594/0.594/0.704; checked by a debug
  print, then rerun: only that row changes). A k equal to a G vector (LaP path G(1 1 1)) is NOT Gamma:
  omega -1.68/-1.68/-4.07 1/cm there (D(G) - D(0) ~1e-5 in omega^2; user: NA/Ewald term, not tested);
  user: avoid such points in the path. finitedif_vel /
  veltest.dat only with FOCEX_VELTEST=1. At exactly q=0 naq returns the angular average of the NA
  term (BZO: 3 x 497.8 instead of 462,462,620), hence the limit. Regression vs d85e3d0 (Ge, BZO,
  MoTe2, GeSe): fits identical; off Gamma omega changes by the removed shift (<= 1.25 1/cm, BZO);
  thermo changes = the old Gamma acoustic modes (0.1-1 1/cm with gamma 5..-530) removed, GeSe
  <gamma> 0.7856 -> 0.7683 (= old data without them, 0.7682); BZO's rest: one mode at the
  |gamma| > 1e3 -> 0 cutoff. THERMACOND has its own get_frequencies (rest.f90), unshifted already.
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

- Build: `cd THERMACOND && make` (mpif90 -> ~/ALADYN/bin/tcond-2.0.x) or `make serial` (mpif90's compiler +
  `mpi_stub.f90` -> tcond-2.0-serial.x); objects in obj_mpi/ obj_serial/; `DEBUG=1` for -g
  -fbounds-check (the old default, ~20x slower). On this Mac mpif90 needs
  `OMPI_FC=/opt/homebrew/bin/gfortran SDKROOT=$(xcrun --show-sdk-path)`: since 2026-11-02 the top-level
  Makefile sets and exports both on macOS when /opt/homebrew/bin/gfortran exists (unless already set), so
  `make all` works in a fresh terminal (user); `make -C THERMACOND` alone still needs them exported.
- Run: `tcond-2.0.x` (serial) or `mpirun -np N tcond-2.0.x`; same code. Ranks > 0 work in `.mpi_rank_NNN/`
  (symlinked inputs, removed at the end) so only rank 0's files survive. Since 2026-10-22
  (degenerate averaging, below) Ge 6^3 job 2 and 4^3 job 3 give byte-identical kappa_vs_T.dat for
  1 and 3 ranks and for GCC 12.2 -O2, GCC 15.3 -O2 and -O0. Before that the result depended on the
  compiler (Ge 4^3 RTA 13.25 / 13.34 / 13.29 W/mK). Since the per-mode columns were added (item (c)
  below) only the RTA and block columns (2-13) are compiler-independent: Ge 4^3 job 3, GCC 12.2 vs
  16.2 (both -O2), per-mode iterative columns 14+ differ by up to 0.2 W/mK at 100 K (degenerate
  eigenvectors); same compiler, serial = MPI byte for byte (2026-10-29). `/usr/local/gfortran` (GCC 12)
  comes first on this Mac's PATH; mpif90 uses Homebrew GCC 16.2 via OMPI_FC. Since 2026-10-29 the
  Makefile's serial FC = `mpif90 --showme:command` (gfortran if no answer), so with OMPI_FC and SDKROOT
  exported both builds use GCC 16.2 (installed serial = MPI byte for byte, Ge 4^3 job 3).
- Shares structure.params and latdyn.params with FOCEX (13-line latdyn, FOCEX reads 1-7;
  template FOCEX/latdyn.params). job (line 8): 0 bands/DOS, 1 +linewidths on the band
  path, 2 +V3*delta, RTA, direct (line 9: 0 CG, 1 SVD), 3 +V3*delta, RTA, iterative.
  All nT temperatures in one run; V3*delta is T-independent: line 10 = 0 writes
  V3DELTA/k######.bin per IBZ point (header nc, ndyn, nibz, k) and reuses it. Since 2026-10-29
  V3DELTA/input_* = copies of lat_fc/fc2/fc3/structure/dielectric; any byte differs -> k*.bin
  removed, recomputed (check_v3_inputs, rank 0). Ge 4^3: job 3 after job 2 reads all 8 points,
  kappa_vs_T.dat = fresh job 3 byte for byte (1 and 3 ranks); changed fc3 -> recomputed. Agent:
  thermacond redone when 06_thermacond/latdyn.params != tc_latdyn() (THERMACOND_* changed).
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

## Read some FC ranks, fit the others (include_fc=2; implemented 2026-11-03, user)

read_fcs (ios.f90; the old read_fcs/_3/_old and readfc removed, user), given_fcs.f90 (xfull, full_col), invariance_ranks.f90 (RHS of the rows, merge_given_fcs),
setup_amatrix.f90 (force rows: column_of / xfull; Huang rows dropped when FC2 is read). Reads fcN.dat (all terms) +
lat_fc.dat of the earlier run (atom numbers -> (tau, n)), least squares per group for the irreducible values; after
the fit the read ranks are merged and marked include_fc=1, so the output and everything downstream is a full fit.
setup_maps now calls get_dim for include_fc /= 0 (sizes of read ranks). Atom numbering: the read files' numbers are
translated through their lat_fc.dat; the output (all ranks, new lat_fc.dat) is in this run's numbering. Tests: REFERENCE/exact_model/coulomb/README.md. The original plan, kept for the record:


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

## Invariance rows: one loop over ranks (2026-10-23, invariance_ranks.f90)

`set_translational_inv_constraints` / `set_rotational_inv_constraints` (module invariance_ranks)
build the rows of every included rank (up to 8) in one pass over the group terms (row key = prefix
atoms + directions [+ nu]; one column function, offset + current(), skipping dropped FC2 terms).
The rank-by-rank routines `set_*_inv_constraints_old`, `FOCEX_INVARIANCE=old` and
`FOCEX_INVARIANCE_COMPARE` (inv_rows_*.dat, inv_fcs.dat) were removed 2026-11-01 (user: remove once
validated; validation below; in git history). include_fc=2 is not handled yet (marked: RHS b).
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
FC (polynomial fit on a 5^m grid) and evaluated inv_rows_*.dat (FOCEX_INVARIANCE_COMPARE=1, removed 2026-11-01; verify.py uses the amat dump) on
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

## Irreducible k-points: one routine for FOCEX and THERMACOND (2026-11-03, user; common/ibz.f90)

module ibz_rules: laue_operations (lattice point-group ops that are a space-group rotation or minus one) and
ibz_reduce (THERMACOND's get_weights algorithm of 2026-10-22: stars modulo G, weight = mesh points in the star,
representative = WS-folded member with the largest u.k, u = (0.9, sqrt2-1, sqrt5-2); mesh lookup by sorted reduced
coordinates, any regular mesh). FOCEX get_weights (extratools.f90) and THERMACOND get_weights (modules.f90) are
wrappers that keep their output files. Removed: FOCEX get_weights2/3/4, getstarp (lattice point group: wrong for a
crystal of lower symmetry than its lattice), THERMACOND kappa9.f90 (not built). Kept: FOCEX getstar (used by the
disabled conductance loop in dynamical.f90 and the uncalled extract_fourier).
- Checked: unit test vs spglib get_ir_reciprocal_mesh (time reversal), 9 structures (Si, wurtzite, CdI2 P-3m1,
  approximate quartz, pyrite Pa-3, P4/m, HfO2, P-1, R-3) x 4 meshes: same stars, weights, wedge, own-star; FOCEX
  regression (Ge, BZO, MoTe2, GeSe): same stars, thermo outputs within 1 in the last printed digit; THERMACOND Ge 4^3
  job 3: every output byte-identical; CdI2 exact model through FOCEX and THERMACOND (8x8x8, 6x6x4): = spglib
  (old FOCEX: 50 instead of 65 IBZ points on 8x8x8). Scratch scripts: ibztest/test.py, ex_cdi2/cmp_spg.py, cmp_tc.py.

## Ballistic conductance (conductance.dat; 2026-11-03, user; test in FOCEX/REFERENCE/exact_model/chain)
- M(om)/A = 2pi sum_lambda int d3k/(2pi)^3 delta(om-om_k) |v.n|/2 in 1/A^2 (user asked for per A^2). |v.n|/2,
  not v.n: the full-zone sum of v.n is 0 by v(-k)=-v(k). Code: tet_sum(|v.n|/2 in m/s)/(c_light*volume_r0*1d-8).
- Mesh values from the IBZ: v(kpc(l)) = rlaue(:,:,kop_ibz2fbz(l)) velocibz(:,j,mapibz(l)) (module kpoints,
  filled by get_weights through ibz_reduce's optional oparm). Replaced the getstar loop (lattice point group,
  disabled with do j=1,-1).
- Verified on exact chains: plateaus 3 and 6 channels/cell to 0.2%, continuous case 0.8% rms, normal (1,2,3) too.
- SC_SNAPS gives no transverse displacements when all its springs (r < 2 r_nn) lie along a chain: use
  Gaussian displacements for such models (setup_chain.py).
- extract_fourier (fourier_stuff.f90) still calls getstar; left as it is (user, 2026-11-03). Known problems:
  global kibz instead of the k* mesh, real fk/uk, mod-based index decoding, no fractional-translation phases.

## Mode Gruneisen parameters (2026-11-01, user's request; common/gruneisen.f90)

One routine for FOCEX and THERMACOND (module gruneisen_rules): gamma = -<e|W(k)|e>/(6 omega^2),
W = sum Psi y exp(ik.R)/sqrt(mm'), y = x + xi, xi = -Gam Tr Q (internal relaxation per unit isotropic
strain, `internal_strain_trace`, from the FC2 terms; = FOCEX's Tr xi of get_phi_zeta_Xi to 2e-15 A). Each
code passes its FC2/FC3 term lists (FOCEX from map/fcs, THERMACOND from fcs_2/fcs_3, first atom in the
primitive cell only) and its dynamical matrix (`get_dynmat`, `get_dynmat_tc`: NA term included).
- The 3 rigid translations at Gamma (eivl=0 from freq_vel) get the long-wave limit gamma = -b/(6a)
  along the direction gamma_limit uses (path segment; mesh: uinc): A = T+D2T - T+D1 G D1T, B = T+W2T -
  T+W1GD1T - T+D1GW1T + T+D1GW0GD1T. D0..D2 from the code's own D at +-h,+-2h,+-3h (h=1e-3|G|; even/odd
  parts fitted exactly), W0..W2 exact. Log: GRUNEISEN_LONGWAVE (checks |D0T|,|W0T|,|T+D1T|,|T+W1T|).
  Derivation + tests: FOCEX/notes/generalized_huang.tex Part III. Also (user, 2026-11-01) the 3 acoustic
  modes at 0 < |k| V^(1/3) < 0.05 (qrad; V^(1/3) = 2pi/gsize) take the limit along k/|k| (modes = largest
  translation_weight). Checked: THERMACOND, HfO2 path 400 pts/segment: inside = exact gamma(k) to <= 4e-4
  (print), just outside off by 5-6e-3 (the shift effect below). qdir normalized inside (IBZ Gamma came in
  as uinc unnormalized: gamma unchanged, a and b scaled).
- Bug fixed: FOCEX built W with exp(iq.(R+taup-tau)) (since fd1297c, 2025-02-17) but its eigenvectors
  come from set_dynamical_matrix_std's exp(ik.R) (R since at least 3ee84a6): gamma right only at q=0,
  off by up to 3.2 on the exact HfO2 model. THERMACOND's routine (removed) had no relaxation.
- Verified on the exact model (FOCEX/REFERENCE/exact_model): gruneisen_q0.py (limit vs direct, CaF2,
  ZnO, HfO2, 4 directions, +-relaxation: q^2 convergence, <= 8e-6 at q=1e-3; direct at q=1e-6 off by up
  to 0.72); mkcubic.py (forces cut after the cubic term, orders separated exactly) + FOCEX ranks 1-3 on
  HfO2 + gruneisen_focex_check.py: FOCEX = exact to the printed 5e-5 at all 519 path points off Gamma,
  THERMACOND to its 4 digits, Gamma limit = -1.691288 -1.344395 -0.658566 in both. bs_grun.dat columns are
  in band-connection order (compare by frequency). Other FOCEX outputs byte-identical except QHA
  (thermo_QHA/temp_free use ibz gamma: HfO2 model <gamma>(300 K) -1.149 -> -0.922).
- Eigenvectors (user, 2026-11-01: the 3e-4 shift is only a numerical device): mode_gruneisen diagonalizes D at
  the point itself and W inside each degenerate set (tol 1e-7 max omega^2), then maps onto the code's bands
  (by omega^2; in a set by overlap with the code's eigenvector). THERMACOND HfO2: = exact gamma(k) to its 4
  digits everywhere (was 9e-3 next to Gamma).
- NA term: formula verified in Python (gruneisen_q0_na.py, ZnO + rigid-ion Z, 1.5e-6 at q=1e-3); the codes'
  route for a polar crystal (D0..D2 fitted from get_dynmat with FOCEX's NA term) not yet run.
- born_flag > 2: WARNING in read_dielectric (both codes; user: use 0 or 2, 1 only for testing); manual says so.
  born_flag=1: nothing subtracted (subtract_coulomb_force returns for <= 1), NA term added to D; the
  log message said otherwise until 2026-11-01 (fixed).

## Elastic constants from sound speeds (elastic_longwave, 2026-11-02; user's request)

- FOCEX/dynamical.f90 elastic_longwave, called at the end of the mech.dat block: along 13 directions
  (axes, face and body diagonals) gruneisen_rules acoustic_tensor (common/gruneisen.f90; the D0/D1/D2
  expansion and G of gruneisen_longwave, factored into longwave_parts, output of gruneisen_longwave
  unchanged) gives A; rho A_ik = sum_jl B[(ik),(jl)] qh_j qh_l, B = the 36 relaxed brackets (rank 36; a general
  Voigt C_(ij),(kl) would have rank 31). B[ab,cd]-B[cd,ab] = d_ab s_cd - d_cd s_ab gives the deviatoric
  stress (pressure invisible: absorbed exactly by a Voigt-symmetric C); C from B - d_ik s_jl with
  s = s_dev + Tr(sigmav)/3. Assumes a force-free reference (else Q's antisymmetric part adds to the asymmetry).
- Switch ewald_no_macro (common/ewald.f90): ewald_2nd_deriv_hat keeps of its G=0 term only (exp(-x)-1)
  (Z.q)(Z.q)/(q.eps.q): E = 0 elastic constants (VASP's). Only for born_flag 2,3,6,7,8; reset after.
- Checks: exact model (REFERENCE/exact_model/longwave): python longwave_check.py (exact Hessian) s_dev 1e-13,
  C 5e-10 GPa for ZnO, HfO2, CaF2 stressed (internal coordinates relaxed); FOCEX on mkstressed.py ZnO inputs:
  s_dev and C to the printed digits; FOCEX relaxed ZnO model: C = exact strain Hessian (1e-3 GPa). User's ZnO
  data born 0: = brackets route within 0.06 GPa, other outputs byte-identical; LaP (centrosymmetric) E=0 = with NA.
- Found with it: mech.dat's bracket route has no long-range term (get_phi_zeta, residuals, mechanical use the
  fitted FC2 only): with born_flag 2 its PHI can be indefinite and C wrong (user's ZnO: C33 -58, C44 -33 GPa).
  User's ZnO born 2: the brackets imply s_zz - p = 0.91 GPa vs 0.01 from FOCEX's stress (born 0: consistent);
  cause not checked (possibly the Huang rows constrain only the short-range FC2). Open.
- Bravais lattices: longwave_parts divided by 0 (P D0 P = 0, no optical modes) -> NaN, zheev in elastic_longwave
  looped forever (argon agent test, 2026-11-02); now G = 0 for nat = 1. Argon: long-wave C = mechanical0's
  bracket C (3.074/1.540/1.537 GPa), B - B^T 2e-11. mechanical0's cubic group-velocity block and the unused
  get_cij removed (user, 2026-11-02): elastic_longwave covers them.
- Argon test (AGENT/tests/fake_vasp.py): TOTEN now after the force block, as VASP; read_outcar.f90 since
  8d633ee wrote a step only when TOTEN followed TOTAL-FORCE, so the test stopped at collect. Since 2026-11-03
  (user) read_outcar pairs the t-th TOTAL-FORCE block with the t-th TOTEN, each found from the top of the file
  (find_nth rewinds; count_word): order-independent. Checked: TOTEN before/after the forces same output; two
  real OUTCARs (4, 20 steps) byte-identical to the 8d633ee build.
- ELASTIC_MODULI table (2026-11-03, user): moduli_from_c (dynamical.f90) applies the mech.dat formulas (Voigt/Reuss/Hill
  B, G, Young, Poisson, 2C44/(C11-C12), universal anisotropy) to the bracket C and to each long-wave C, crystal frame.
  Before, the moduli came from the bracket C only. inverse_real (common/tools.f90) overwrites its first argument (no
  explicit interface: an intent(in) actual is not caught): pass a copy. Ge: bracket column = the old mech.dat lines.
- One atom per cell: the G (optical-mode) terms of A and B are skipped explicitly (gruneisen_longwave,
  acoustic_tensor; user); outputs byte-identical (argon, ZnO born 0).
- born2 ZnO asymmetry, located (2026-11-03; scratch coulomb_stress/coul_brackets.py, Ewald E=0 of the Z*/eps
  dipole FCs alone, lam-independent to 4e-11): the SR brackets satisfy Huang (6e-5 GPa), the Coulomb FCs do
  not: their unrelaxed brackets carry a deviatoric stress (-0.108,-0.108,+0.216) GPa, and their Q has
  exactly the residual-force form Q(ta;bg)-Q(ta;gb) = pi_g d_ab - pi_b d_ag (misfit 2e-9) with pi_z = -0.133
  (Zn), +0.133 (O) eV/A: the dipole FCs are those of a stressed, not force-free charge lattice. FOCEX imposes
  the invariances (rotational, Huang) on the SR part with the total pi and sigma, so SR + Coulomb violates
  them by pi_C, sigma_C. 0.216 of the 0.906 GPa is the stress; the rest (through pi_C and G) not yet
  computed. Possible fix (not decided): SR rows with RHS -pi_C (rotational) and sigma - sigma_C (Huang).
- Block-by-block check of the born_flag-2 Coulomb model (2026-11-03, user; scratch coulomb_stress/coul_force.py,
  coul_force2.py: ewaldforce/dewapot reimplemented, F1 = e^2 Z1 sum_2 Z2 (-grad phi(r12)), Ewald): F0 = +0.1826
  (Zn z) = FOCEX log. Translational: F(x+c) = F(x) to 5e-15, q=0 FC sums zero (exact by construction).
  Rotational: scalar Z and eps -> F(Rx) = R F(x) to 5e-9 and dF/deta_bg - dF/deta_gb = F_g d_ab - F_b d_ag to
  6e-10 (user's expectation holds); with ZnO's tensor Z, eps: F(Rx) - R F(x) = 0.49 of the first-order change,
  rotational form violated by 0.178 eV/A (Z, eps fixed in space). Per pair (one Zn-O pair, no images): the force
  formula's FC (Z1 Z2)_ac d_b d_c phi is not symmetric (xz 0.170, zx 0.188) and differs from the Z1 T Z2^T
  added back to D (xz = zx = 0.179); q=0 sums agree only by symmetry. So born_flag 2 subtracts a harmonic part
  /= the D_NA added back when Z is anisotropic (possibly the BZO zone-boundary OPEN item; not tested).
- BUG fixed 2026-11-03 (GaN, user's run): dfc_hyb (Phi_C - Phi_m, ewald_born.f90) was a saved array computed on the
  first subtract_coulomb_force call only; with two POSCAR/FORCEDISP pairs (24 and 72 atoms) the second supercell read
  past its end: bforce of structure 2 garbage (-3.6e5), fit error 2.98 / 3.24 / 712 / 4.95 % on identical runs
  (user's 655 %), one segfault, one hang. Now computed per call (per supercell) and freed. GaN born 2: 1.719 % on
  3 runs, byte-identical; single-supercell runs byte-identical to before. Found by FOCEX_DUMP_AMAT: only bforce rows
  of structure 2 differed between runs. -fcheck=all did not report the out-of-bounds read (not checked why).
  The exact-model tests used one supercell only, so they could not see it.
- chi (dielectric, 2026-11-03, user): modes now from D at q = 0 with ewald_no_macro (TO; born_flag 2,3,6,7,8), not
  q = 1e-6*normal (there the mode along normal was LO: GaN Im chi_xx at E1(LO) 726, eps0_xx 8.26 vs eps0_yy 10.01;
  now both 10.0065, peaks at E1(TO) 553.9 and A1(TO) 529.9, the transverse modes of the old q||x run). Acoustic modes
  left out by translation_weight, not as bands 1-3. gamma_mode_symmetry(ndn,to): to=.true. writes chi_modes.dat
  (TO modes, irreps), chi.plt labels the IR-active peaks from it. gamma_modes.dat unchanged.
- GaN instability located (2026-11-03, open): the Huang rows. Soft rows (enforce_inv 0), born 2, both sets:
  transl+rot 1.714 %, C11/C12/C33/C44/C66 = 313.6/99.5/303.7/91.6/107.1 GPa, no imaginary modes; transl+Huang
  1.213 %, C66 -51.3, -17 1/cm; all three (user's) C66 -21. Not the stress (no stress.params: same) and not the
  Coulomb part of the RHS (max 0.18 eV). Set 2 alone: Huang rows give -3425 1/cm, without them C66 83 GPa, stable.
  The stable solution violates Huang rows 14, 15 ([zz,xx]-[xx,zz] type) by -6.2 eV (RHS 0.12); mech.dat
  |[ab,gl]-[gl,ab]| 10.5 of 504 GPa. Cause not found. ZnO: same symptoms (user).
- SOLVED 2026-11-04: GaN/ZnO Huang and rotational failures were caused by positions slightly off the symmetry,
  not by ill-conditioning (the 2026-11-03 entry saying so was wrong). GaN Ga z 0.5005917/0.0005866 (6_3 screw off by
  2.7e-5 A), ZnO 8.4e-5 A: below tol_sym 1e-3, so the ops are found and the FCs symmetrized exactly, but the
  Huang/rotational rows use the actual positions: rows zero by symmetry ([xy,xz],[yz,xx],[yz,yy] in 6mm, odd under
  C2z) got coefficients 1.4e-3 (max 334) and acted as spurious constraints (dFC2 1.0 eV/A^2 even with their RHS 0;
  only the real row [zz,xx]: dFC2 3e-3, fit 1.202 -> 1.214 %). Checked from FOCEX_DUMP_AMAT in Python and by runs.
  Fix: symmetrize_positions (force_constant_symmetries.f90, called after make_unitcell, symmetry_init rerun) and
  pointgroup_average2 on stress.params (setup_amatrix, read_stress_params). GaN 1 1 1 0: C66 -7288 -> 110, stable;
  1 1 1 1: 1.94 %, stable. ZnO 1 1 0 1: 58.7 % (-2e5 1/cm) -> 1.30 %; enforce_inv 0/1/2 agree within 4 GPa.
  Exact models: positions symmetric to 1.6e-10 A, hence never affected. Still open: ZnO long-wave C33 45-49 GPa
  (exp ~210) in all enforce_inv; GaN/ZnO bracket-route C not positive definite (short-range FC2 only, born 2).
- Checks of the fix (2026-11-04): FCs about the symmetric reference (POSCAR supercells shifted by the same delta, so
  u is measured from it): FC2 changes 7e-4 eV/A^2 (4e-5 rel), FC3 4e-4 rel, C <= 0.5 GPa: no re-expansion needed.
  relax_reference K vs Phi, Huang check (GPa): 1100 Phi 9.15 -> K 9.11; 1110 2.4e-4 -> 1.2e-2; 1111 4e-7 -> 2.0e-2
  (relaxation does not repair Huang; why K loses 2e-2: unknown). Relaxation from symmetric positions leaves
  u0(atom 2)_z = 1.05e-6 A (cause unknown). Exact ZnO+Coulomb ranks 1-4 regression: FCs 8e-10 rel, freqs identical.
- enforce_inv=1, born_flag 2: max|C xpart - d| = 2.5e-3 eV (GaN): rotational rows' Coulomb RHS (|b| <= 0.22)
  inconsistent with the translational rows at ~1 %; Huang rows consistent (2e-7). Z*, eps replaced by their traces:
  2.4e-7. So the cause is the tensor anisotropy (subtracted model: harmonic part full-tensor Phi_C, anharmonic part
  trace charges). OPEN, user to decide the treatment. User first said no 9.6 until resolved;
  on 2026-11-04 the user asked to tag 9.6 now with these items open.
- 9.6 (2026-11-04, user): which positions are used where (read from the code). symmetry_init builds atompos from
  the symmetrized atompos0, so the rotational rows (invariance_ranks build_inv_rows: atompos(mu,l)), the Huang rows
  (setup_amatrix: rr = atompos(j)-atompos(i)) and relax_reference (x = atompos) use symmetric positions. The
  translational rows have no positions. Displacements stay FORCEDISP - POSCAR (read_force_position_data:
  dsp - atom_sc%equilibrium_pos) and forces as read. FC1: the irreducible FC1 obey the symmetry, but sum_tau
  Phi_tau = 0 is the rank-1 translational row, exact only with enforce_inv 1/2 (GaN enforce_inv 0: -2e-5 eV/A).
- After 9.6 (2026-11-04, user): (1) relax_reference subtracts the average force per atom from FC1 before relaxing
  (sum Pi = 0; logged before/after); fc1_unrelaxed.dat keeps the fit; pi_prime still removes the net force of the
  Phi S + ... terms (FC2-4 ASR inexact with enforce_inv 0), printed at the end. Wording corrected: the relaxation
  never removed the net force from the FCs, only ignored it. (2) Born charges on (born_flag /= 0, Z /= 0): the
  bracket-route C and everything from it (moduli, sound speeds, Debye T, B, B0') go to the log only, with a
  warning; mech.dat keeps ELASTIC_LONGWAVE and the moduli table without the bracket column. Bracket sums use the
  fitted FC2 only (get_phi_zeta_Xi: fcs). NOT changed: QHA uses rbulkmod/hbulkmod and B0' of the bracket route.
  (3) cubic_elastic call commented: b0p always overwritten by pressure_derivative_b (no early return) when FC3 is
  in, and the QHA uses 1+2 gamma without FC3; its value was wrong (ZnO -216.6 vs 3.57).
  (4) Log comments for the Huang numbers and the Qiu lines. (5) ELASTIC_LONGWAVE: pass label (E = 0 / NA) now
  before its rho v^2 < 0 warnings (they used to appear under the previous block). ZnO (user's DFT, relaxed,
  1 1 0 1): the E = 0 C is unstable along 101/111 (C33 45, C13 141); negative log sound speeds come from the
  bracket C (C44 -31.5, C66 -33.2), sound_speeds itself checked by hand (mysqrt = -sqrt|x|).
  (6) fc2range = 0 range from the largest-volume supercell: confirmed by swapping the GaN POSCAR/FORCEDISP order
  (same lgridmax, 186 kept FC2, fc2 to 1.2e-8); log line of find_WS_largest_SC printed i instead of imax (fixed).
  (7) sound_speeds: log note for rho v^2 < 0. (8) OPEN: relax_reference has no Coulomb/Born term: with born_flag 2
  it relaxes with the short-range FC2 only (gradient misses Phi_C S, strain Hessian misses the Born brackets).
  Relaxation tested only on the exact ZnO model without Coulomb; not on monoclinic or with charges.
- 2026-11-04 (afternoon): relax_reference now has Phi_C (born_flag 2): phi_c (lr_c at q = 0, E = 0, self term
  by ASR -- first version overwrote it: sum_tau phi_c 4.6, Ga pair split 1.6e-5 A; fixed, 6e-15), q_c (lr_moments),
  a_c (lr_brackets), r -> F r when deformed; Psi_C, Chi_C empty term lists. GaN strain xx/zz 1.46e-3/-2.35e-3
  -> 5.9e-6/6.8e-4. sgfract re-derived after relaxing (GaN relaxed had no gamma_modes.dat/chi_modes.dat: 6_3
  translation off by 1.7e-3 > tol_sym, also in 9.6). Supercell inconsistency messages (read_supercell volume
  check, report_unmapped in check_input_poscar_consistency_new). Monoclinic test model HfO2r (setup.py, model3
  'Rb' option), run_HfO2r(_q): results pending at this commit. Note: coulomb_lr keeps the reference positions
  after a relaxation (ready flag), used by the mech.dat Huang/rotational checks.
- Night 2026-11-04 (user asleep; asked: commit and push if diagnosed and fixed, report what/how/tests):
  exact monoclinic HfO2 (ranks 1-3, stressed, +/- point charges), structure.params positions perturbed by 5e-5 A:
  9.5 fits 92.5 % / 1429 %, FC2 wrong by 0.91 / 7.6 rel; new code 0.000 %, FC2 1.7e-6 rel (the undetectable
  symmetric part of the perturbation). relax_reference pi_prime ignores the net force (it cannot be
  relaxed; soft-fit ASR error had left 1e-6 A between equivalent Ga; the fitted FCs keep it). Cross-validation: Huang neutral for force prediction, rotational rows
  double the held-out error (open). Record: FOCEX/REFERENCE/exact_model/symmetrize/README.md. Manual entry
  "After 9.5 (not tagged)"; Makefile ver stays 9.5 (user: no 9.6 until everything is resolved).
- Z*/eps symmetrization (check_born_symmetry, ios.f90; 2026-11-03, user; uncommitted): group average after the ASR,
  epsinv recomputed. Exact ZnO+Coulomb: change 4e-16, FCs equal to 1e-9 relative, sorted frequencies identical.
- Band Gruneisen divergence (GaN, user's question): gamma = -<e|W|e>/(6 omega^2); the numerator omega^2*gamma is
  smooth along G-M (-239, -772, -1573, -2553, -3640 1/cm^2), omega^2 crosses 0 between dk 0.052 and 0.069: the spike
  (-116) is the 1/omega^2 of an unstable harmonic branch, not large anharmonicity.
- GaN (2026-11-03; cause found 2026-11-04: positions off the symmetry, see SOLVED above): enforce_inv=1 with rotational rows conflicts with the data: set 2 alone, ranks 1-2:
  born 0 13.3 %, born 2 22.4 % (translational only, born 2: 7.7 %); ranks 1-3 (FC3 2 shells): born 0 36 %, born 2
  123 %; translational only, ranks 1-3, both sets: 1.20 %. Cause not found (truncated FC2 of a polar crystal for
  born 0; for born 2 possibly the Coulomb RHS of the rotational rows, rank 3 in particular). Long-wave C of the
  born 2 fit: C11 186 < C12 228, C66 -21 GPa, small imaginary TA near Gamma along G-M.
- born_flag 2 hybrid (2026-11-03, user: "go ahead with the hybrid"; details and tests in
  REFERENCE/exact_model/coulomb/README.md): subtracted F = F_m(x) - F_m(x0) - (Phi_C - Phi_m)u, F_m = Ewald force of
  charges tr(Zi Zj^T)/3 (ewaldforce), Phi_C = Coulomb dipole FCs (= D_NA added back), Phi_m its trace-model harmonic part
  (hybrid_harmonic, d2ewapot). Invariance rows get RHS -A x_sub whenever they are included (no switch, user):
  Huang -(S_C(ab,gd)-S_C(gd,ab)) (setup_amatrix), rank-2 rotational -M:eps (invariance_ranks), S_C and M from
  coulomb_lr.f90 (small-q E=0 Ewald). mech.dat Huang/rotational checks add the same terms. Exact tests (point charges,
  ZnO ranks 1-4, HfO2 cubic ranks 1-3): exact to the printed digits; each RHS is needed. Rank-3+ rotational rows still
  homogeneous at first; since 2026-11-03 rank-3 rows get -A_full x_sub too (lr_rot3_pair / lr_rot3_onsite; checked
  against an independent Python route, README). Rank 4 for tensor eps: OPEN. Build in a clean folder: a copied build
  folder kept a stale dynamical.o once (make did not recompile).

## Elastic constants, B0', QHA (FOCEX 9.2, 2026-11-01; user's questions)

- Elastic constants (elastic_born_huang, mechanical and mechanical0): c = Sym{[ag,bd]+[bg,ad]-[ab,gd]
  + (d_ab s_gd + d_gd s_ab)/2 - d_ag s_bd - d_bd s_ag} - Sym(QGQ)/V, [ab,gd] = atld_raw(a,g,b,d). Until
  2026-11-01 Sym(A) - Sym(d_ag s_bd): off-diagonal and shear constants wrong (exact relaxed ZnO model:
  C12 171 vs 221, C13 137 vs 197, C44 108 vs 78, C66 146 vs 121 GPa; diagonal right). Checked
  (REFERENCE/exact_model/bprime/elastic_check.py vs the model Hessian): <= 6e-10 GPa relaxed ZnO, CaF2 and
  CaF2 under 14.5 GPa; 0.2 GPa on the stressed ZnO/HfO2 references (residual forces, not in the formula).
  All earlier FOCEX elastic constants (and B, G, sound speeds, Debye T) were affected.
- B0' (pressure_derivative_b, overwrites cubic_elastic's b0p): B(+-lam) of the crystal deformed along the
  relaxed unit-pressure path p1 (FC2 + lam Psi.S, deformed d, pressure), B = V/(v Ht^-1 v), Ht = V(c+S2)+PV2;
  B' = (B(+)-B(-))/2lam, lam = 1e-4 eV/A^3. The A-form W3 of the notes (p H[Psi S] p) is wrong like Sym(A).
  Verified: FOCEX on bprime/mkrelaxed.py inputs (relaxed model, D=0, C=0, cubic-truncated forces):
  ZnO B 288.001543 B' -0.157427, CaF2 587.376258 / 0, = full relaxations under pressure (bp2_all.log).
  Assumes hydrostatic reference stress, no residual forces (warnings).
- QHA (QHA_free_energy_vs_strain): Murnaghan static EOS (B0 = Reuss, B0' from FC3; no FC3: 1 + 2 gamma_Leont'ev,
  Dugdale-MacDonald, user 2026-11-02; gamma formula = Agne et al. SI Eq. S17), phonon F
  with omega(V) from gamma; V(T) = root of P (cubic), B_T = -V dP/dV there, alpha_V = (dP/dT)_V/B_T
  (identity; summed analytically), lin CTE = alpha_V/3 = dlnL/dT (user: the definition), Cp = Cv + VTB_T
  alpha_V^2. Harmonic-only (no FC3): Leont'ev gamma for all modes + P_agne = NkT/V (user's term) as before.
  Checked: alpha = numerical dlnL/dT of the table to 0.8 % (ZnO model, with FC3) and 0.7 % (FC2 only).
  Before: V(T) Boltzmann-weighted over strains (kT_eff = range/100), CTE secant eq_strain/T, B(T) =
  B0(1-strain(4+3<gamma>)), static energy quadratic in dV (B0' = -1).

## Fit errors sd(rank) (2026-11-01, user: make sure it is correct)

sigma(j) = standard error of FC j: s sqrt([K (A_K^T A_K)^-1 K^T]_jj) (regularized as svbksb), s^2 =
||Ax-b||^2/(m-p) (svd_set; it now also returns V diag(s sqrt f(w)) so fit_kernel_basis maps it by K).
sd(rank) = mean sigma of the irreducible FCs of the rank (log only). Before: enforce_inv=1 never set
sigma (sd from an unset array), enforce_inv=0 scaled by sqrt(svdcut), sd divided again by sqrt(rows).
LASSO and itemp=1: 0 (not estimated). Verified: exact ZnO model ranks 1-3 + Gaussian force noise 1e-3:
s = 0.9948e-3 / 1.0005e-3 / 0.9776e-3 (seeds 1, 2; enforce_inv=1, 0); (x-x_exact)/sigma rms FC2 1.03,
FC3 0.84-0.90 (FC3 errors 10-15 % below the estimate, cause not known), 72-76 % within 1 sigma.
Jackknife (user's idea, done once, not in the code): 5 refits each leaving out 20 of the 100 noisy ZnO
snapshots, sqrt(4/5 sum (x_i - xbar)^2) / SVD sigma: median FC2 1.04, FC3 0.91; rms ratio 0.96.

## Things that bite

- No GNU-extension intrinsics: `transform_input_structure` (`common/lattice.f90`) used
  `dcosd/dsind`, which Rivanna's gfortran does not have (`Function 'dcosd' has no IMPLICIT
  type`; the Mac's gfortran 10.3 accepts them, so it built here). Standard `cos/sin` with
  `d2r = pi/180` since 2026-11-03, with the near-zeros snapped to 0/1 so right angles stay
  exact as `dcosd` made them. No Makefile flag can help: `-fall-intrinsics` only matters when
  `-std=` has switched extensions off, and the Makefiles pass no `-std`. All 16 `common/`
  modules compile under `-std=f2008`; keep it that way.
- `sc_snaps.x` built before 2026-10-22 aborts with `Recursive call to nonrecursive
  procedure 'findif'` when compiled with `-fcheck=all` (`d2v` called `findif(r,d1v)`, and
  `d1v` calls `findif`). Fixed in SC_SNAPS: `d2v` writes the same central difference out;
  identical freqs.dat/modes.dat, any flags work. Rebuilt copies in `~/BIN` and `~/ALADYN/bin`.
- The parameter files are read positionally by list-directed `read`. A missing
  line or a missing value silently swallows the next record and crashes later
  with a confusing message. `default.params` has 4 lines (the 32 maxterms/maxtermzero/maxtermsindep/maxgroups lines were
  removed 2026-10-27: setup_maps sizes its arrays itself; an old 36-line file still works); line 2 of
  `latdyn.params` must carry **six** numbers (shift *and* normal).
- `bands.dat` columns 6+ are ω², not ω (`bandgen.plt` applies `521.1*sqrt()`). Plot scripts and `bandplot` are in
  `FOCEX/PLOTTING/` since 2026-11-02 (user); the agent reads its templates from there, `make` links bandplot.
- `FOCEX/PLOTTING/allplots` (2026-11-02, user): fills bandgen/bandgrun/bvel/all/qha/thermo/fcs/chi.plt for a
  FOCEX folder and makes the PDFs; KTICS.BS for the xtics line; caps frequency 2000 1/cm, group velocity
  20 km/s, Gruneisen [-10:10]; NAME = folder above 05_focex, else the parent folder (--name overrides). all.plt colours the
  dispersion by the band Gruneisen parameters (bs_grun.dat, same band order as bands.dat: checked row by
  row) when present. The agent's makeplots.sh runs it (per_band removed: no template uses its marker).
  Checked on the exact-model ZnO run (all 9 PDFs, pages looked at) and with capped/missing data.
- Group velocities (bs_freq.dat, bands.dat, ibz_bands.dat) are Hellmann-Feynman since 2026-10-22:
  `finitedif_vel` returns `v1` from `get_freq`; before, it returned the finite difference of
  energy-sorted bands (BZO crossing: 0.17/8.50 km/s FD vs 0.65/9.08 HF). Both stay in
  `veltest.dat`. Regression on Ge, BZO, MoTe2, GeSe: only those velocity columns changed.
- The log line "Initial Pressure(GPa), gr,xagne" prints an uninitialised `gr` (~1e-314).
- fc*.dat, fc*_irr.dat carry 8 significant digits since 2026-11-03 (g15.8, user; 7 before); symmetry-equivalent FCs
  then differ at ~1e-8. Readers are list-directed (FOCEX, THERMACOND), so the wider field is harmless.
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
- Snapshot energies (read_force_position_data, 2026-11-02, user): converted to eV per primitive cell, then
  the lowest per-cell value of FORCEDISP1 subtracted, so several supercells are on one scale (log: "Energy
  assigned (eV per primitive cell ...)" with file and cell count). Before, emin was supercell 1's lowest
  TOTAL energy, subtracted before the division. Checked: Ge example tiled 2x along a1 as POSCAR2/FORCEDISP2
  (2E): identical per-cell energies for 108 and 216 cells; one supercell: same values as before.
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
- Array sizes: no user-set limits (default.params has 4 lines); size_force_constants grows its own
  work arrays (start 64/16/64, doubling). The "maxterms too small" retry messages of force_constants
  are no longer printed (2026-10-27); the final work sizes are in the log (SIZE_FORCE_CONSTANTS:
  work bounds). A real mismatch still stops: "SETUP_MAPS: sizing disagreement".
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

`origin` is **HTTPS** (`https://github.com/KeivanS/Anharmonic-lattice-dynamics.git`; checked
2026-10-29). Authentication: `gh` is logged in as KeivanS (token in the macOS keychain, git
protocol https) and git's `osxkeychain` helper supplies it, so `git fetch`/`git push` need no
password and worked with the Bash sandbox on (push of 6d8c541). SSH does not work here: no key in
`~/.ssh` (the id_ed25519 of 2026-10-27 is gone), `ssh -T git@github.com` -> Permission denied; user
chose HTTPS + gh token (no new key). The repo, VASP-Flow and SC_SNAPS are public: pulling needs no
login anywhere, pushing from another machine needs its own `gh auth login`. The SC_SNAPS clone at
`~/PROJECTS/SC_SNAPS` is on HTTPS too. `WebFetch` and sandboxed `curl`
cannot reach github.com at all. Other hosts (readthedocs.io, PyPI) are fine.

## SCOP8: static strain and u0 -> FCs, review for FOCEX (2026-11-02, user: note problems, fix later)
Read (not run: SCOP8 does not compile). Notation to reuse in FOCEX: strain(3,3) = eta, atomic_deviation(3,tau)
= u0, S(:,i) = eta.(R_i+tau_i) + u0(type_tau(i)) (static displacement), trialfc2_value = effective FC2,
myfc2_value%phi / myfc3_value%psi / myfc4_value%chi, GradientV_utau / GradientV_eta / GradientV_cor.
- Correct (VA_math.f95 GetV_avg_And_GradientV_avg, UpdateTrialFC2 in check.f90): static part of the effective FC2
  = 2 dV/d<YY> = Phi2 + sum_m Psi S_m + 1/2 sum_mn Chi (S_n-S_i)(S_m-S_i) (+ the <YY> term); equals the Taylor
  re-expansion at x0+S (the -S_i parts vanish by ASR on m, n).
- PROBLEM 1: FC1 (residual forces) is not in <V> nor in GradientV_utau/GradientV_eta: myfc1 is used only in
  `residuals` (diagnostics). An unrelaxed reference is treated as if relaxed (missing pi0.S, i.e. the
  residual forces and the reference stress they carry).
- PROBLEM 2: `residuals` sigma0 = sum_tau pi0 tau/V: image-dependent when pi0 /= 0 (FOCEX replaced it by
  stress.params); no stress.params in SCOP8.
- TO CHECK 3: add_Pressure adds V(eta) eta:sigma; the work of a constant external stress is sigma:(change of
  shape), which differs at second order in eta (hydrostatic: P dV = P V0 (tr e + ((tr e)^2 - tr e^2)/2 + det e)
  vs P V0 (1+tr e+..) tr e). Algebra only, not tested.
- TO CHECK 4: mechanical2 uses atld0 - QGQ - sigma0(be,ga) delta(al,de): the pre-9.2 FOCEX form (Sym(A) -
  stress term) gave wrong off-diagonal/shear constants; not compared.
- Cubic "gradient w.r.t strain, term 2" (IF found ... END IF) is empty; intent unknown.

## Relaxing an unrelaxed reference with the fitted FCs (relax_reference.f90, 2026-11-02, user)
default.params line 5 = 1 (only 1 switches it on; absent/other = off). After the fit: minimize the Taylor
expansion over strain eta and u0 (S = eta(R+tau)+u0, paper/SCOP8 notation), re-expand K = Phi+Psi S+Chi SS/2,
Psi' = Psi+Chi S (S_j - S_tau forms), project on the irreducible FCs, strain lattice/positions/grids by
M = 1+eta (+u0); downstream uses them; stress.params read as 0 (params: relax_ref, ref_relaxed); unrelaxed
fit kept as fcN_unrelaxed.dat, fcN_irr_unrelaxed.dat, lat_fc_unrelaxed.dat; POSCAR_relaxed written.
- The paper's A = -1/2 sum K d d as strain Hessian is WRONG (first version gave ZnO strain -0.01831 vs
  exact -0.017837). Cause (2026-11-02, user asked why): substituting u = eta X into 1/2 sum u_i Phi_ij u_j is
  right; the step to -1/4 sum (u_i-u_j) Phi_ij (u_i-u_j) is not: it needs Phi_ij^{ab} = Phi_ij^{ba}, but only
  Phi_ij^{ab} = Phi_ji^{ba} holds. The ASR form drops the part of each 3x3 block antisymmetric in its
  Cartesian indices (0 for central pair forces, Cauchy case). Checked (REFERENCE/exact_model/relax/cluster_identity.py, 40-atom
  ZnO cluster of the exact model): E2 from energies = 1/2 u Phi u = 1061.81 / 435.42 eV, the -1/4 form
  1025.80 / 471.43 = 1/2 u Phi_sym u exactly; crystal (relax/paper_A_check.py, ZnO at the minimum, sigma=0):
  Vc Sym(A) C12-type 49.91 vs exact 53.96 eV, Born-Huang 6e-3 eV off. Same cause as the pre-9.2 elastic
  constants. Q (linear in d, row ASR only) is right. Now: stage 1 relaxes u0 in
  the reference cell (stress += int Q^T du, exact), stage 2 Newton in strain with u0 relaxed at every RK4
  point (force-free path) and PK2 stress integrated with J F^-1F^-1 c_rel F^-TF^-T, c_rel = elastic_born_huang
  (current brackets, current stress) - QGQ/V.
- elastic_born_huang lacked the (ab)<->(gd) average: at a non-hydrostatic stressed state its tensor is not
  pair-symmetric (ZnO model at 15 GPa, internal relaxed: H13 51.452 / H31 51.420, exact 51.436 = average;
  relax/paper_A_check.py, its first part). Average added in elastic_born_huang (all callers; hydrostatic: no change).
- Terms whose first atom is outside the primitive cell are skipped (as get_phi_zeta_Xi).

- Checked (exact ZnO model, stressed reference 14.9 GPa, FC1 0.032 eV/A, ranks 1-4, e1; REFERENCE/exact_model/relax/
  model_min.py, compare_relaxed.py): FOCEX vs the model's own minimum: strain 7.8e-10, u0 1.3e-9 A, FC2 at the
  minimum 2.7e-6 eV/A^2 (7 printed digits), FC3 5e-7; residual stress 2e-10 GPa, Huang 2e-8 GPa; 3 strain
  steps. Ranks 1-3 only: strain off 1e-5 (FC4 missing). Ge (line 5 absent): elastic constants identical.

## IR / Raman activity at Gamma (gamma_mode_symmetry, dynamical.f90, 2026-11-02, user)
Called after dielectric; writes gamma_modes.dat. D at exactly q = 0 (NA term direction-averaged: keeps the
full symmetry; at normal*1e-6 the LO-TO split breaks it). Operation g maps atom k to img (r' = S r + sgfract,
as check_born_symmetry); chi(g) = Tr(E+ Gam(g) E) per degenerate set; IR: (1/|G|) sum chi* Tr S, Raman:
(1/|G|) sum chi* ((Tr S)^2 + Tr S^2)/2, parity chi(i)/dim. No character tables.
- Checked: ZnO exact model (C6v): A1 + 2B1 + E1 + 2E2 (E2 Raman only, B1 silent, A1/E1 both); HfO2 exact
  model (C2h): 18 g Raman-only, 15 optical u IR-only; LaP (Oh, Born): T1u IR yes, Raman no, u, strength
  0.962 = 3 x 0.3207 of dielectric's table.
- Character tables (2026-11-02, user: computed, not looked up; modular for electrons/magnetic later):
  common/point_group.f90, module point_group_rules: pg_setup(nop, ops, ulog, ier) (multiplication table,
  classes by conjugation, characters by the class-multiplication method: class matrices made normal by
  D^-1/2 M D^1/2, D = diag(class sizes), Hermitian + anti-Hermitian parts diagonalized together with zheev),
  pg_reduce(chi per op, multiplicities), pg_write_table, pg_chi_vector / pg_chi_sym2 (IR / Raman).
  Labels: A/B (cubic: all A), E, T, g/u, '/'', serial number only when needed; classes in a fixed order
  (E, C2, C3, C4, C6, i, m, S6, S4, S3; principal-axis C2 and sigma_h first), so the numbering does not depend
  on the order of the operations (it did: ZnO gave the IR E mode as E2); E1 = chi(C2 principal) < 0 with a
  6-fold axis, T1 = chi(C4 or S4) > 0 (as the standard tables); other 1/2 (B1/B2, A1/A2 where conventions
  differ) may not be Mulliken's; complex conjugate pairs (C3, C4, C6, S4, S6, C3h, C4h, C6h, T, Th) as
  Ea/Eb, pg_partner.
  Unitary operations only (magnetic: corepresentations to add; spin: double groups to add).
  Checked: FOCEX/REFERENCE/point_groups/test_point_groups.f90, all 32 groups: order, classes, dimensions,
  orthogonality (ier 0); Oh T1u = 3 1 0 -3 -1 -1 -1 0 1 1 as the standard table.
  FOCEX gamma_modes.dat (checked): ZnO Gamma = 2A1 + 2B2 + 2E2 + 2E1 (IR+Raman E1 243, Raman-only E2 196/404,
  silent B 332/682; Mulliken calls the silent ones B1: sigma_v/sigma_d choice); HfO2 9Ag + 9Au + 9Bg + 9Bu,
  acoustic Au + 2Bu; LaP 2T1u.
- Symmetry tolerance (2026-11-02, user): tol_sym = 1e-3 A (common/params.f90, parameter) for accepting the
  operations AND finding atom images: symmetry finder (force_constant_symmetries tol), findatom_cart (was
  ncmp, i.e. default.params tolerance 1e-4 per component), check_born_symmetry and gamma_mode_symmetry (were
  1e-4 fractional); all as Cartesian distance. ncmp/tolerance unchanged elsewhere (FC equation coefficients).
  User's ZnO (Zn pair 8e-5 A off the 6_3 screw, from the relaxation): stopped with "Error0 ... atom not
  found"; now runs, fit 1.756 % (symmetrized positions: 1.757 %), Gamma frequencies within 0.3 1/cm. User: do
  not move the atoms. HfO2 exact model: fc2, mech, bs_freq, gamma_modes byte-identical. THERMACOND builds.
- gamma_mode_symmetry: acoustic set = all modes with translation_weight > translation_min (was |omega| < 1
  1/cm; ZnO with Born charges has -2.6/2.0 1/cm at exactly q = 0, then listed as IR-active optical).
