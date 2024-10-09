#!/bin/bash
#set -o noexec
#set -o verbose
#set -o xtrace
#BSUB-q short 
#BSUB-o output
#BSUB-J jobash
#BSUB-n 1 gf
#############################################################################

JOBNAME="168P.16"
nibz=$(cat ./kibz.dat)
nibz_proc=2     # First choose nibz_proc and num_proc = nibz/nibz_proc
r=$((nibz%nibz_proc))
d=$((nibz/nibz_proc))
if [ $r -ne 0 ] ; then
   s=$((1+$d))
else
   s=$d
fi

num_proc=$s   # dependent variable

sed  "s/xxx/$nibz_proc/" params.phon00 > params.phon11 
sed  "s/yyy/5 /" params.phon11 > params.phon22

for (( i = 1; i <= $num_proc; i++ ))
do

# for two digit num_proc
if [ $i -lt 10 ] ; then
   DIRNAME="k00$i"
   D="00$i"
else
   if [ $i -lt 100 ] ; then
       DIRNAME="k0$i"
       D="0$i"
   else
       DIRNAME="k$i"
       D="$i"
   fi
fi

rm -r $DIRNAME
mkdir $DIRNAME
cd ./$DIRNAME

cp ../kpbs.params .
cp ../structure.params .
cp ../params.phon22 ./latdyn.params
cp ../dielectric.params .
cp ../lat_fc.dat .
cp ../fc2.dat .
cp ../fc3.dat .
cp ../fc4.dat .
cp ../kap8 .


start_ibz=`echo $nibz_proc*$i-$nibz_proc+1 | bc`
if [ $i -ne $num_proc ] ; then
 end_ibz=`echo $nibz_proc* $i  | bc`
else
 end_ibz=`echo $nibz | bc`
fi

cat > ksubset.inp << EOF
$start_ibz $end_ibz
EOF

cat > kap1.sh << EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p standard
#SBATCH -A elmgroup_standard
#SBATCH -t 01:30:00
#SBATCH --mail-type=ALL
###SBATCH --mail-user=


srun ./kap8

cp ./Cs_sycollisionm.dat ../Cs_sycollisionm.$D.dat
##cp ./RHS_IBZ.dat ../
##cp ./FRTA.dat ../

EOF

chmod +x kap1.sh
sbatch kap1.sh


cd ..

done
