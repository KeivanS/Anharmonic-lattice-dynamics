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
(`cd FOCEX/utility && make`, installs to `make.inc`'s `DESTDIR=~/ALADYN/bin`).

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

## Things that bite

- `sc_snaps.x` aborts at run time with `Recursive call to nonrecursive procedure
  'findif'` **when built with `-fcheck=all`** (which `make.inc` sets):
  `-fcheck=all` implies `-fcheck=recursion`, and `findif` does recurse via
  `d1v`/`d2v`. Verified: `-O2` fine, `-O2 -fcheck=all` aborts, adding
  `-frecursive` fixes it. Not changed in the Makefiles — documented only.
- The parameter files are read positionally by list-directed `read`. A missing
  line or a missing value silently swallows the next record and crashes later
  with a confusing message. `default.params` must have 36 lines; line 2 of
  `latdyn.params` must carry **six** numbers (shift *and* normal).
- `bands.dat` columns 6+ are ω², not ω (`bandgen.plt` applies `521.1*sqrt()`).
- `conductance.dat` is written but its band loop is `do j=1,-1`, so it is all
  zeros.
- `POSCARi` uses the **VASP-4** layout — no element-name line. `POSCAR1` is
  `sc_snaps.x`'s `poscar_000` with line 6 deleted.
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

- `github.com/KeivanS/SC_SNAPS`, cloned at `~/PROJECTS/SC_SNAPS` — the current
  home of SC-SNAPS, with a Flask GUI. Its `cell.inp` format has changed relative
  to the copy bundled in `FOCEX/utility/` (line 3 takes two values, line 6 takes
  masses *and* charges). Documented in `docs/source/runscsnaps.rst`; that repo
  has its own `CLAUDE.md` with the details.

## Where the git clone is

**This directory is not a git repo** — it is an unzipped snapshot. The clone to
commit from is:

```
"/Users/ke4c/Documents/Documents - MAE-/GitHub/Anharmonic-lattice-dynamics"
```

(note the space in the path). Work there, or edit here and port the files over —
as of 2026-08-07 the two trees are identical apart from deliberate edits. That
clone also holds untracked user work (`FOCEX/SNAPS/`, `FOCEX/todo`) and a local
uncommitted `.gitignore` line — leave all three alone. It was 214 commits behind
origin/main on 2026-08-07; check with `git fetch && git status -sb` before
assuming it is current.

`git` over HTTPS to GitHub works **only with the Bash sandbox disabled**;
`WebFetch` and sandboxed `curl` cannot reach github.com at all. Other hosts
(readthedocs.io, PyPI) are fine.
