# FOCEX — FOrce Constant EXtraction

FOCEX, part of ALATDYN (Anharmonic LATtice DYNamics), extracts the force
constants (FCs) of rank 2, 3, 4 and above — up to rank 8 — from force–displacement
data of atoms in one or several supercells, computed with DFT (or any other
method). From the FCs it then calculates the phonon spectrum and a range of
thermodynamic properties: elastic constants, total and free energy, vibrational
entropy and heat capacity, Gruneisen parameter and thermal expansion
coefficient. In the current version the volume dependence of the frequencies is
obtained, for simplicity and speed, through the Gruneisen parameter computed
from the formula proposed by Leontiev.

Other codes in the suite, such as THERMACOND, additionally compute phonon
lifetimes, scattering rates and thermal conductivity.

**Documentation:** <https://aladyn.readthedocs.io/en/latest/runfocex.html>
(source in `docs/source/runfocex.rst`) — it describes the whole workflow, every
line of every input file, and every output file.

## Installing

FOCEX needs nothing but a Fortran compiler (gfortran or Intel):

```bash
cd Anharmonic-lattice-dynamics/FOCEX
make                    # -> ~/ALADYN/bin/focex-8.19.x
```

Edit the `FF` and `FLAGS` variables at the top of the `Makefile` to change the
compiler or its options; the executable is `focex-<ver>.x` (version `ver`, now 8.19),
installed in `DESTDIR` (default `~/ALADYN/bin`; `make DESTDIR=...` for another place).

The helper programs that convert DFT output (`read_outcar.x`, `read_qe.x`,
`poscar2xyz.x`) are built separately:

```bash
mkdir -p ~/ALADYN/bin   # or the DESTDIR set in ../make.inc
cd utility
make
```

`poscar2focex.x POSCAR` (built by `make` above, no external library) writes
`structure.params`, `POSCAR_focex` (the cell to build the supercells from) and
`kpbs.params` from a POSCAR, with the cell and the band path of Setyawan and
Curtarolo (2010). `utility/poscar2focex.py` (numpy + spglib) is the earlier
version, kept for `focex-agent.py`.

The snapshot generator `sc_snaps.x` is not built here but from `../SC_SNAPS`
(`make sc_snaps` in the top folder installs it in `~/ALADYN/bin`), which also has a
graphical interface. That folder is identical to <https://github.com/KeivanS/SC_SNAPS>.

## Running

Put the seven input files in one directory and run the binary there:

```
structure.params  dielectric.params  latdyn.params  kpbs.params  default.params
POSCAR1  FORCEDISP1        (POSCAR2/FORCEDISP2 ... for additional supercells)
```

A complete, ready-to-run example (germanium, 216-atom supercell, 36 snapshots)
is in `example/Ge`; see `example/INPUT_FILE_INFO.md`.

## Plotting

The `*.plt` files in this directory are gnuplot scripts for the standard figures
— dispersion (`bandgen.plt`, `bvel.plt`, `all.plt`), force-constant decay
(`fcs.plt`, `pair.plt`), thermodynamics (`qha.plt`, `thermo.plt`), dielectric
response (`chi.plt`) and the k-meshes and Wigner–Seitz cells (`kp.plt`,
`ws.plt`). Run them in the directory holding the output files.
