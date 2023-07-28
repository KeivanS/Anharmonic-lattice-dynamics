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

sed -e 's/xxx/ 2 1 0 /' params.phon1 > params.phon2

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

cd ./$DIRNAME

cp ../params.phon2 ./params.phon

cat > kap2.sh << EOF
#!/bin/bash 
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p standard
#SBATCH -A elmgroup
#SBATCH -t 00:03:00 
#SBATCH --mem-per-cpu=30000
#SBATCH --mail-type=ALL 
####SBATCH --mail-user=sn7sb@virginia.edu 
 

srun ./kaptet
mv v33sq_delta.dat ../v33sq_delta.$DIRNAME.dat

EOF


chmod +x kap2.sh
sbatch kap2.sh

cd ..

done

sed -e 's/xxx/ 1 1 1 /' params.phon1 > params.phon

cat > run.sh << EOF
#!/bin/bash 
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p standard
#SBATCH -A elmgroup
#SBATCH -t 00:30:00 
#SBATCH --mem-per-cpu=30000
#SBATCH --mail-type=ALL 
####SBATCH --mail-user=sn7sb@virginia.edu 
 

srun ./kaptet

EOF

#sbatch run.sh
