# SC-Snaps GUI

A browser-based graphical interface for the `sc_snaps.x` supercell snapshot
generator. Runs locally as a Flask web app — no internet connection required.

---

## What it does

`sc_snaps.x` generates thermally-displaced supercell snapshots for use in
phonon calculations. This GUI lets you:

- Edit the three input files (`cell.inp`, `snaps.inp`, `supercell.inp`)
  directly in the browser with format hints (a default will first be uploded)
- Save edits to disk explicitly before running
- Run `sc_snaps.x` and watch the live output stream
- Browse generated `poscar_*` files and the snapshots.xyz file
- Click on any generated poscar or snapshots.xyz to visualize the generated structures using Jmol — all with one click
- Vesta or similar visualization softwares are also possible; 
  just replace the jmol path in the configuration section

---

## Requirements

| Requirement | Notes |
|-------------|-------|
| Python ≥ 3.8 | |
| [Flask](https://flask.palletsprojects.com/) | required to start the GUI |
| NumPy, Matplotlib, SciPy | only for the analysis plots; without them the GUI still runs and the plot panels report the missing package |
| [gfortran](https://gcc.gnu.org/fortran/) | or a similar Fortran compiler, to build `sc_snaps.x` from `sc_snaps.f90` |
| [Jmol](https://jmol.sourceforge.net/) | structure visualizer (optional) |

---

## Installation & first-time setup

```bash
git clone https://github.com/KeivanS/SC_SNAPS.git
cd SC_SNAPS
make venv           # one-time: create .venv and install the Python packages
make compile        # compile sc_snaps.f90 -> $BINDIR/sc_snaps.x  (default ~/BIN)
```

`make venv` creates a local virtual environment with `--system-site-packages`,
so any NumPy/Matplotlib you already have is reused and only the missing packages
are downloaded. Once `.venv` exists the Makefile uses it automatically; nothing
needs to be activated.

You need the virtual environment on any Python marked *externally managed*
(PEP 668) — Homebrew's and most Linux distributions' — where `pip install flask`
is refused with `error: externally-managed-environment`. If your Python is not
externally managed you can skip `make venv` and simply
`pip install flask numpy matplotlib scipy`.

## To open the GUI in the browser window

```bash
make run                        # opens the GUI in the browser  http://localhost:5050
```

## If you want to reopen the GUI in the browser window, first kill the existing one using

```bash
lsof -ti :5050 | xargs kill -9  # to kill the GUI browser window in case you want to rerun it
```


```makefile
PYTHON     = .venv/bin/python   # or python3 when there is no .venv
SC_SNAPS_X = ~/BIN/sc_snaps.x
# VISUALIZER = ~/BIN/jmol           # Jmol shell script or .app bundle
# POSCAR2XYZ = ~/BIN/poscar2xyz.py  # POSCAR → XYZ converter (is provided; if you move it to ~/BIN, provide the path)
```

All paths set here become the default values shown in the browser on startup.
They can be overridden at any time in the browser without restarting.

---

## Input files

### `cell.inp` — primitive cell
```
1 1 1   90 90 90          # conventional cell: a b c  α β γ
 0 0.5 0.5, 0.5 0 0.5, 0.5 0.5 0   # primitive vectors (in conventional units)
4.247  9                  # lattice parameter scale (Å), 0 for coordinates in conventional cell units, otherwise coordinates are read in primitive cell units
2                         # number of atom types
1 1                       # number of atoms of each type
24.31  16.00              # atomic masses
Mg O                      # element names
  0 0 0                   # reduced coordinates (conventional lattice if 0 on line 3, otherwise primitive)
  0.5 0.5 0.5
```

NOTE: if you want to directly input the primitive cell information, use 1 1 1 90 90 90 on the first line and use the second line for the cartesian units of the primitive cell AND  on the thrid line use scale (=1) followed by a non-zero number like 9

### `snaps.inp` — snapshot parameters
```
400    # average phonon frequency (cm⁻¹)
30     # temperature (K) 
51     # number of snapshots  supercell_type 
```
T=30K is a reasonable value for temperature. Check the displacements (modes.dat) distribution in the GUI to make sure it is what you want.

### `supercell.inp` — supercell dimensions
```
3 0 0
0 3 0
0 0 3
```
3×3 integer matrix: supercell vectors in terms of the primitive cell vectors.
A diagonal matrix gives an n×n×n supercell, but you can also use non-diagonal elements.

---

## Output

After running, the working directory contains:

| File | Description |
|------|-------------|
| `poscar_000` | Unshifted reference supercell (POSCAR format) |
| `poscar_001` … `poscar_N` | Thermally displaced snapshots |
| `poscar_j.xyz` | XYZ-format versions (created on demand by clicking "Convert & open in Jmol") |
| `snapshots.xyz` | all the snapshots and velocities if needed 

---

## Workflow after you type "make run" and get the browser window

1. **Set working directory** — type a path (new or existing) and optionally
   click **Load existing files** to populate the text areas from disk
2. **Edit** the three input files in the text areas
3. Click **💾 Save input files** to write edits to disk
4. Click **▶ Run sc_snaps.x** — files are saved and the job starts;
   live output streams in the log box
5. After completion, generated snapshots (poscar_j) appear as grey chips:
6. Click any chip to view its content inline
7. Click **👁 Convert & open in Jmol** (POSCAR) or **👁 Open in Jmol** (XYZ)
   to launch the visualizer

---

## File structure

```
sc-snaps-gui/
├── sc-snaps-gui.py      # Flask server + single-page browser app
├── Makefile             # make compile / make run
└── README.md
```

---

## License

MIT — see [LICENSE](LICENSE) for details.
