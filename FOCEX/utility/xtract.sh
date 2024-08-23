#!/bin/bash
# extracts forcedisp from dft output
rm FORCEDISP1
for x in OUTCAR*
do
cp $x OUTCAR
read_outcar.x
cat pos-forc.dat >> FORCEDISP1 
done
