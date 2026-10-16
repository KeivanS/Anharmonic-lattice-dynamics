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
cd FOCEX && mkdir -p ~/BIN && make       # ~4 min; produces ~/BIN/v16
```

`FOCEX/Makefile` does **not** include the top-level `make.inc`: it hardcodes
`FF=gfortran`, `FLAGS=-O2 -fcheck=all`, names the binary from `exe=v16`, and
`mv`s it to `~/BIN`, which it does not create. The utility programs are separate
(`cd FOCEX/utility && make`, installs to `make.inc`'s `DESTDIR=~/ALADYN/bin`) —
that builds `read_outcar.x`, `read_qe.x` and `poscar2xyz.x` only.
`FOCEX/utility/poscar2focex.py` (numpy + spglib, not in the Makefile): POSCAR →
`structure.params` (Niggli-reduced, conventional axes along x,y,z; 9-number line 1,
line 2 = spglib transformation, atoms in conventional units, other lines = Ge
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

- `sc_snaps.x` (in the SC_SNAPS repo) aborts at run time with `Recursive call to
  nonrecursive procedure 'findif'` **when built with `-fcheck=all`** (which this
  repo's `make.inc` sets, so it bites anyone reusing those flags):
  `-fcheck=all` implies `-fcheck=recursion`, and `findif` does recurse via
  `d1v`/`d2v`. Verified: `-O2` fine, `-O2 -fcheck=all` aborts, adding
  `-frecursive` fixes it. SC_SNAPS's own Makefile uses plain `-O2`.
- The parameter files are read positionally by list-directed `read`. A missing
  line or a missing value silently swallows the next record and crashes later
  with a confusing message. `default.params` must have 36 lines; line 2 of
  `latdyn.params` must carry **six** numbers (shift *and* normal).
- `bands.dat` columns 6+ are ω², not ω (`bandgen.plt` applies `521.1*sqrt()`).
- `conductance.dat` is written but its band loop is `do j=1,-1`, so it is all
  zeros.
- `structure.params` line 3 is `scal [convcoord]`. It is read as text, so the
  old one-value form (`5.70 # scale`) still works and means `convcoord=0`.
  `convcoord=0`: atom positions are in units of the conventional vectors
  (line 1); otherwise in units of the primitive vectors (line 2, columns of
  `prim_to_conv`). Both are converted to Cartesian in `make_unitcell` before
  `unitcell`, which expects Cartesian input.
- `structure.params` line 1 is **6 numbers** `a b c alpha beta gamma` (FOCEX
  then puts a along x, b in the xy plane) **or 9 numbers**, the conventional
  vectors a1 a2 a3 in Cartesian, kept as given (since 2026-10-15; detected by
  whether 9 reals parse). POSCARi/FORCEDISPi/dielectric tensors are Cartesian
  and must be in the same frame, so a DFT cell not oriented a∥x needs the
  9-number form — otherwise "rs1 in units of r0i" comes out *nearly* integer
  and FOCEX stops with "supercell vector k of POSCAR ... is not integer"
  (e.g. AgCrP2S6_PARA, whose cell has a along y; it is 2×1×2 and runs with
  poscar0's vectors on line 1).
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
