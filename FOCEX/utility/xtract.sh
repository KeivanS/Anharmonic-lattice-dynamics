#!/bin/bash
# Collects FORCEDISP1 for FOCEX from finished VASP runs stored here as OUTCAR1, OUTCAR2, ...
rm -f FORCEDISP1
for x in OUTCAR[0-9]*
do
cp $x OUTCAR
read_outcar.x
cat pos-forc.dat >> FORCEDISP1
done
