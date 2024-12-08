#!/bin/bash
export OMP_NUM_THREADS=2  # parallelization option, to be tested/optimized by the user
for x in poscar_*   # these are snapshot outputs of sc_snaps.x 
do
cp $x POSCAR
mpirun -np 8 vasp_std
# this reads the vasp output OUTCAR and extracts the positions and forces into pos-forc.dat
read_outcar.x           
# collect the force-positions from every snapshot into a single FORCEDISP1 file to be read by FOCEX
cat pos-forc.dat >> FORCEDISP1    
done
