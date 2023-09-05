#!/bin/bash
workdir=`pwd`
if [ $# -ne 1 ];
then
echo "To execute the script"
echo "./process_dft.sh vasp_file_or_directory --> for VASP output"
echo "./process_dft.sh qe_file_or_directory --> for QE output"
fi
if [ $# -eq 1 ];
then
count_outcar=`find $1 -name "OUTCAR"|wc -l`
if [ ${count_outcar} -eq 0 ];
then
echo "No OUTCAR file present or directory containing it"
fi
if [ ${count_outcar} -ne 0 ];
then
echo "${count_outcar} OUTCAR file present"
for (( i=1; i<=${count_outcar}; i++ )){
each_outcar=`find $1 -name "OUTCAR"|sed -n ${i}p`
each_outcar_dir=`find $1 -name "OUTCAR"|sed -n ${i}p|awk -F "OUTCAR" '{print $1}'`
cd ${each_outcar_dir}
~/aladyn/readoutcar.x
wait $!
cat pos-forc.dat >> ${workdir}/FORCEDISP1
rm pos-forc.dat
cd ${workdir}
}
fi
fi
