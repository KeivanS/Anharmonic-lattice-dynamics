# FOCEX examples

Each subfolder is a complete, self-contained FOCEX run: it holds every input
file needed, plus the output of a previous run for comparison.

The full description of every input and output file is in the manual,
[Running FOCEX](https://aladyn.readthedocs.io/en/latest/runfocex.html)
(source: `docs/source/runfocex.rst`). This page is only a quick summary.

## Ge

Germanium in the diamond structure: a 216-atom supercell (3x3x3 of the 8-atom
conventional cubic cell, a = 5.7022565 Ang) and 36 snapshots in which all atoms
are thermally displaced. Ranks 1, 2 and 3 are fitted.

### Running it

```bash
mkdir -p ~/runs/Ge && cd ~/runs/Ge
cp <repo>/FOCEX/example/Ge/{structure.params,dielectric.params,latdyn.params,kpbs.params,default.params,POSCAR1,FORCEDISP1} .
~/BIN/v16                 # or whatever you named the FOCEX binary
```

The run is serial and takes a few minutes on one core. Expected results:

| quantity | value |
| --- | --- |
| snapshots read | 36 |
| irreducible FC2 / FC3 | 78 / 95 (173 unknowns in total) |
| translational / rotational / Huang constraints | 14 / 84 / 15 |
| fit error `||F_dft - F_fit|| / ||F_dft||` | 0.036 % |
| C11, C12, C44 | 121.1, 68.0, 47.4 GPa |
| Bulk modulus (Hill) | 85.7 GPa |
| Debye temperature | 313.6 K |

### Input files

| file | contents |
| --- | --- |
| `structure.params` | primitive cell, ranks to fit, FC ranges, invariance and fit options |
| `dielectric.params` | Born-charge flag, dielectric tensor, Born effective charges (all zero here: Ge is non-polar) |
| `latdyn.params` | k-mesh, DOS mesh and broadening, temperature range for the thermodynamics |
| `kpbs.params` | k-point path for the phonon band structure |
| `default.params` | numerical thresholds and array limits; a `0` on a line means "use the built-in default" |
| `POSCAR1` | equilibrium positions of the supercell, old (VASP-4) POSCAR format — **no element-name line** |
| `FORCEDISP1` | the 36 snapshots: positions/displacements, forces and total energy |

`FORCEDISP1` in this example uses the "vasprun" block format, in which the
header line carries the energy after `(eV):` and each atom line contains a
displacement (in Bohr) followed by a force. The alternative — and preferred —
format, produced by `utility/read_outcar.x` and `utility/read_qe.x`, starts each
block with a line containing the word `POSITION` and lists absolute positions in
Ang and forces in eV/Ang. Both are described in the manual.

### Output files

The interesting ones are `fc2.dat`/`fc2_irr.dat` and `fc3.dat`/`fc3_irr.dat`
(the full and the irreducible force constants, in eV/Ang^2 and eV/Ang^3),
`lat_fc.dat` (the structure information that goes with them — the other ALATDYN
codes need it), `bs_freq.dat` and `bands.dat` (phonon dispersion), `ibz_dos.dat`
(density of states), `mech.dat` (elastic tensor and moduli), `thermo_QHA.dat`
(quasi-harmonic thermodynamics) and the `log*.dat` file, which records
everything the code decided and should be read after every run.

> **Note.** The output files stored in `Ge/` were produced by an older version
> of the code. Their numbers, and the column layout of the `fc*.dat` files,
> differ slightly from what the current version writes. Use a fresh run as the
> reference.
