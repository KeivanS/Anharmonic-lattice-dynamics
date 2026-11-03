#!/bin/bash
HERE="$(cd "$(dirname "$0")" && pwd)"
SCF_DIR="$HERE/.."
# Guard: 02_scf may have been converted in place for LOBSTER (ISYM=0/-1).
# Its CHGCAR is still the symmetric converged density (the conversion is a
# fixed-charge NSCF with LCHARG=.FALSE.), but its WAVECAR is a full-mesh
# symmetry-off one that this ISYM=2 step cannot use — skip it and let VASP
# regenerate the orbitals from CHGCAR.
ISYM_SCF=$(sed 's/[!#].*//' "$SCF_DIR/INCAR" 2>/dev/null \
           | awk -F= 'toupper($1) ~ /^[ \t]*ISYM[ \t]*$/ {gsub(/[ \t]/,"",$2); print $2}' \
           | tail -1)
if [ "$ISYM_SCF" = "0" ] || [ "$ISYM_SCF" = "-1" ]; then
    echo "  NOTE: 02_scf/INCAR has ISYM=$ISYM_SCF (converted for LOBSTER) — skipping its"
    echo "        WAVECAR; DFPT (ISYM=2) will start from CHGCAR instead."
    if ! sed 's/[!#].*//' "$SCF_DIR/INCAR" | grep -qiE 'ICHARG[ \t]*=[ \t]*11|LCHARG[ \t]*=[ \t]*\.FALSE\.'; then
        echo "  WARNING: 02_scf is symmetry-off and NOT a fixed-charge NSCF — its CHGCAR"
        echo "           may have been generated with symmetry off. Consider re-running"
        echo "           a symmetric (ISYM=2) SCF before DFPT."
    fi
else
    cp "$SCF_DIR/WAVECAR" "$HERE/" 2>/dev/null && echo "  WAVECAR copied from 02_scf" \
        || echo "  WARNING: WAVECAR not found in 02_scf"
fi
cp "$SCF_DIR/CHGCAR" "$HERE/" 2>/dev/null && echo "  CHGCAR copied from 02_scf" \
    || echo "  WARNING: CHGCAR not found in 02_scf"
if [ -f "$SCF_DIR/POSCAR" ] && [ -s "$SCF_DIR/POSCAR" ]; then
    if [ "$SCF_DIR/POSCAR" -nt "$HERE/POSCAR" ]; then
        cp "$SCF_DIR/POSCAR" "$HERE/POSCAR"
        echo "  POSCAR updated from 02_scf (POSCAR is newer)"
    else
        echo "  Keeping existing POSCAR (local copy is newer than 02_scf)"
    fi
else
    echo "  WARNING: 02_scf/POSCAR not found; keeping existing POSCAR"
fi
