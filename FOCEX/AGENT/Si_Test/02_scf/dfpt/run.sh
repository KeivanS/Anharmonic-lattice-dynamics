#!/bin/bash
export OMP_NUM_THREADS=${VASP_OMP_NUM_THREADS:-1}

# run.sh — VASP DFPT: Born effective charges + dielectric tensor
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

if [ ! -f "$HERE/POTCAR" ]; then echo "ERROR: POTCAR not found"; exit 1; fi
[ -f "$HERE/copy_from_scf.sh" ] && bash "$HERE/copy_from_scf.sh"

echo "Starting Born-charge / dielectric response run at $(date)"
mpirun -np 12 ~/BIN/vasp_std > vasp.out 2>&1

grep -q "reached required accuracy" vasp.out 2>/dev/null \
    && echo "  OK: converged" \
    || echo "  WARNING: may not have converged — check vasp.out"

# Extract Born charges → BORN file for phonopy NAC correction
if [ -f "$HERE/OUTCAR" ]; then
    echo "Extracting Born charges ..."
    python3 "$HERE/extract_born.py" "$HERE/OUTCAR"
fi
echo "DFPT done. See born_charges.txt and BORN (phonopy NAC format)."
