#!/bin/bash
#set -o noexec
set -o verbose
set -o xtrace
#BSUB-q short 
#BSUB-o output
#BSUB-J jobash
#BSUB-n 1 gf
#############################################################################

JOBNAME="168P.16"
nibz=64
nibz_proc=10   # First choose nibz_proc and num_proc = nibz/nibz_proc
num_proc=7  # dependent variable
alphamix=1

sed -e 's/yyy/1 /' params.phon0 > params.phon1

sed -e 's/xxx/ 0 1 0 /' params.phon1 > params.phon2

for (( i = 1; i <= $num_proc; i++ ))
do

# for two digit num_proc
if [ $i -lt 10 ] ; then
 DIRNAME="00$i"
else
   if [ $i -lt 100 ] ; then
       DIRNAME="0$i"
   else
       DIRNAME="$i"
   fi
fi

rm -r $DIRNAME
mkdir $DIRNAME
cd ./$DIRNAME

cp ../kpbs.in .
cp ../params.inp .
cp ../params.phon2 ./params.phon
cp ../params.born .
####cp ../merge.inp .
cp ../lat_fc.dat .
cp ../fc2.dat .
cp ../fc3.dat .
cp ../kap7_sy_tet .


start_ibz=`echo $nibz_proc*$i-$nibz_proc+1 | bc`
if [ $i -ne $num_proc ] ; then
 end_ibz=`echo $nibz_proc* $i  | bc`
else
 end_ibz=`echo $nibz | bc`
fi

cat > ksubset.inp << EOF
$start_ibz $end_ibz
EOF

cat > kap.sh << EOF
#!/bin/bash 
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p standard
#SBATCH -A elmgroup
#SBATCH -t 00:03:00 
#SBATCH --mem-per-cpu=30000
#SBATCH --mail-type=ALL 
###SBATCH --mail-user=sn7sb@virginia.edu 
 

srun ./kaptet
mv v33sq.dat ./v33sq.$DIRNAME.dat

EOF

chmod +x kap.sh
sbatch kap.sh

cd ..

done

