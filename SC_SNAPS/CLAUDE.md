# SC-SNAPS — working notes

Generates thermally displaced supercell snapshots for phonon / force-constant
work. Two parts: `sc_snaps.f90` (the generator) and `sc-snaps-gui.py` (a local
Flask web GUI). Feeds **FOCEX** in the ALATDYN suite
(`github.com/KeivanS/Anharmonic-lattice-dynamics`) — that repo has its own
`CLAUDE.md`, and the user-facing manual for *this* code lives there as
`docs/source/runscsnaps.rst`.

**Where to edit (2026-10-23):** this code is the `SC_SNAPS/` folder of
Anharmonic-lattice-dynamics, a git subtree of `github.com/KeivanS/SC_SNAPS`. Edit it there and
publish with `git subtree push --prefix=SC_SNAPS sc_snaps main` from that repository; a change
made directly in SC_SNAPS comes back with `git subtree pull --prefix=SC_SNAPS sc_snaps main`.

## Ground rule

Verify against the source and an actual run. The input files are read
positionally by list-directed Fortran `read`; a missing value silently swallows
the next line and fails somewhere else entirely. Build in a scratch directory —
`gfortran` drops `.mod` files next to the sources.

## Build and run

```bash
make compile                 # gfortran -O2 -> $BINDIR/sc_snaps.x  (default ~/BIN)
make run                     # GUI on http://localhost:5050, opens a browser
lsof -ti :5050 | xargs kill -9   # port is hardcoded; free it before restarting
```

Deps: Flask (required for the GUI), numpy + matplotlib + scipy (only for the
plot panels — imported lazily, so the GUI runs without them and the plot
endpoints return `{"error": "Missing library..."}` with HTTP 500). `scipy.stats`
is used for the Gaussian fits in the component and velocity histograms, so it is
genuinely needed, not just named in the error string. Jmol optional.

The user's `python3` is Homebrew 3.14 and is `EXTERNALLY-MANAGED` (PEP 668):
`pip install` into it is refused. `make venv` creates `.venv` with
`--system-site-packages` (brew already provides numpy and matplotlib as the
`numpy` and `python-matplotlib` formulae) and installs the rest. The Makefile
then picks `.venv/bin/python` up automatically via a `wildcard` test.

Any flags work since 2026-10-22, `-fcheck=all` included: `d2v` used to call `findif(r,d1v)`,
which re-enters `findif` (`d1v` calls it) and aborted with `Recursive call to nonrecursive
procedure 'findif'` under `-fcheck=recursion`. `d2v` now writes the same central difference out
(same step h); freqs.dat and modes.dat are byte-identical to before.

## Input format (this is the part that has drifted)

`cell.inp` — **differs from the older copy bundled in ALATDYN's `FOCEX/utility`**:

| line | contents |
| --- | --- |
| 1 | `a b c alpha beta gamma` of the conventional cell (angles in degrees) |
| 2 | 9 numbers: primitive vectors in conventional units (commas OK as separators) |
| 3 | **`scale, convcoord`** — two values. `convcoord=0` → atom coords *and* `supercell.inp` are in conventional units; ≠0 → primitive units |
| 4 | number of atom types |
| 5 | atoms of each type in the primitive cell |
| 6 | **`ntype` masses then `ntype` charges**, one line |
| 7 | element names (2 chars) |
| 8+ | one line per atom, reduced coords; grouped by type; trailing text ignored |

`supercell.inp` — 3×3 integer matrix, rows = supercell vectors, in conventional
or primitive units per `convcoord` (the README says "primitive", which is only
half right). `snaps.inp` — avg frequency (cm⁻¹), temperature (K), n snapshots
(≤999: `character(3) csn`, and units are `70+s`).

Outputs: `poscar_000` (undistorted) … `poscar_NNN`, `snapshots.xyz`, `log.dat`,
`freqs.dat`, `modes.dat`. POSCARs are Cartesian with 6 columns (position +
velocity). **`POSCAR1` for FOCEX = `sed '6d' poscar_000`** (FOCEX wants the
VASP-4 layout with no element-name line).

## Fixed here (2026-08-07)

- `defaults.json` — the `cell` template was missing the `convcoord` value on
  line 3, so the GUI's out-of-the-box default crashed `sc_snaps.x` with
  `Bad integer for item 1 in list input`. This is the first thing a new user
  sees. Note the Python `DEFAULT_CELL` in `sc-snaps-gui.py` was already correct;
  `defaults.json` overrides it.
- `Makefile` — `make compile` did `mv -f poscar2xyz.py $(BINDIR)/iposcar2xyz.py`:
  a typo *and* destructive (removed the file from the working tree, while the
  GUI's converter default points at `poscar2xyz.py`). Now `cp -f` under the
  right name.

## Still open

- `README.md`'s `cell.inp` example uses `4.247  9` (primitive-unit coordinates)
  but then gives `0 0 0` / `0.5 0.5 0.5`, which is rocksalt in *conventional*
  units — should be `4.247 0`. The clone URL in the README is also still
  `YOUR_USERNAME/sc-snaps-gui.git`.
- The GUI has no `poscar_000`-to-`POSCAR1` button; users must remember the
  `sed '6d'` step.
- `_parse_cell_inp()` / `_parse_supercell_inp()` read `WORKDIR`, not the
  browser-selected working directory, unlike the other endpoints.

## Environment note

`WebFetch` and sandboxed `curl` cannot reach github.com from this machine, but
`git` over HTTPS **does** work when the Bash sandbox is disabled. So: read repos
from local clones, and expect `git push`/`ls-remote` to work.
