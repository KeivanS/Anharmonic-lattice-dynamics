#!/bin/bash 
#SBATCH -N 2
#SBATCH --ntasks-per-node=25
#SBATCH -p parallel
#SBATCH -A elmgroup_standard
#SBATCH -t 05:20:00 
####SBATCH --mem-per-cpu=9000
#SBATCH --mail-type=ALL 
####SBATCH --mail-user=sn7sb@virginia.edu 

module purge
module load gcc/11.4
module load openmpi 

srun ./kap8

