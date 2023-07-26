#!/bin/bash
# extracts forcedisp from dft output
rm OUTCAR
for x in OUTCAR*
do
cp $x OUTCAR
readoutcar.x
cat pos-forc.dat >> FORCEDISP1 
done
