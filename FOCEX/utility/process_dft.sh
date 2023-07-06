#!/bin/bash
if [ -f "readoutcar.x" ];
then
rm -rf ./readoutcar.x
fi
if [ -f "readpfpwscf.x" ];
then
rm -rf ./readpfpwscf.x
fi
if [ $# -eq 0 ];
then
echo "Please provide the OUTCAR file or QE outputfile as an argument or" 
echo "directory containing VASP outcars or directory containing QE runs"
fi
if [ $# -eq 1 ];
then
make all
if [ -d $1 ];
then
echo "Directory detected: directory name, $1"
totaloutcars=`find $1 -name "OUTCAR*"|wc -l`
echo "Number of OUTCARs present: $totaloutcars"
if [ -f "FORCEDISP" ];
then
rm FORCEDISP
fi
if [ -f "./OUTCAR" ];
then
rm ./OUTCAR
fi
for (( i=1; i<=$totaloutcars; i++ )){
eachoutcar=`find $1 -name "OUTCAR*"|sed -n ${i}p`
cp $eachoutcar .
./readoutcar.x
wait $!
cat pos-forc.dat >> ./FORCEDISP
rm pos-forc.dat
}
rm ./OUTCAR
rm -rf *.x
fi
if [ -f $1 ];
then
echo "OUTCAR file provided"
./readoutcar.x OUTCAR
wait $!
mv pos-forc.dat FORCEDISP
rm -rf *.x
fi
fi
