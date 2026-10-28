#!/bin/bash
# Runs VASP on every snapshot poscar_* written by sc_snaps.x and collects the positions and forces
# into FORCEDISP1 for FOCEX. Needs INCAR (single point: IBRION = -1, NSW = 0), KPOINTS, POTCAR here.
export OMP_NUM_THREADS=2  # parallelization option, to be tested/optimized by the user
rm -f FORCEDISP1
for x in poscar_*   # these are snapshot outputs of sc_snaps.x
do
cp $x POSCAR
mpirun -np 8 vasp_std
# this reads the vasp output OUTCAR and extracts the positions and forces into pos-forc.dat
read_outcar.x
# collect the force-positions from every snapshot into a single FORCEDISP1 file to be read by FOCEX
cat pos-forc.dat >> FORCEDISP1
done
