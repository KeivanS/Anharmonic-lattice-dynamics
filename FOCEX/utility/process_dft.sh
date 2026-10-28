#!/bin/bash
# Collects FORCEDISP1 for FOCEX from finished VASP runs: every OUTCAR below the given directory
# (one snapshot per subdirectory). FORCEDISP1 is written in the directory the script is run from.
if [ $# -ne 1 ]; then
   echo "usage: process_dft.sh <directory containing the snapshot runs>"
   exit 1
fi
workdir=$(pwd)
rm -f FORCEDISP1
n=0
for outcar in $(find "$1" -name OUTCAR | sort); do
   (cd "$(dirname "$outcar")" && read_outcar.x && cat pos-forc.dat >> "$workdir/FORCEDISP1" && rm pos-forc.dat)
   n=$((n+1))
done
echo "$n OUTCAR files collected into FORCEDISP1"
