#!/bin/bash
export OMP_NUM_THREADS=${VASP_OMP_NUM_THREADS:-1}

# run.sh for relax
# VASP binary : ~/BIN/vasp_std
# MPI ranks   : 12
# Launch      : mpirun -np 12 ~/BIN/vasp_std

set -e
HERE="$(cd "$(dirname "$0")" && pwd)"

# POTCAR is a symlink to the project-level POTCAR built by vasp-agent.py.
if [ ! -f "$HERE/POTCAR" ]; then
    echo "ERROR: POTCAR not found in $HERE"
    echo "  Re-run vasp-agent.py to rebuild it."
    exit 1
fi

# SCF: pull relaxed geometry from 01_relax/CONTCAR
if [ -f "$HERE/copy_from_relax.sh" ]; then
    bash "$HERE/copy_from_relax.sh"
fi

# Bands/DOS: pull CHGCAR + POSCAR from 02_scf
if [ -f "$HERE/copy_from_scf.sh" ]; then
    bash "$HERE/copy_from_scf.sh"
fi

# Run VASP
echo "Starting relax (~/BIN/vasp_std) at $(date)"
mpirun -np 12 ~/BIN/vasp_std > vasp.out 2>&1

# Convergence check
if grep -q "reached required accuracy" OUTCAR 2>/dev/null; then
    echo "  OK: converged at $(date)"
else
    echo "  WARNING: may not have converged -- check OUTCAR"
fi
