#!/bin/bash 
#SBATCH -N 2
#SBATCH --ntasks-per-node=17
#SBATCH -p parallel
#SBATCH -A elmgroup_standard
#SBATCH -t 1:20:00 
#SBATCH --mail-type=ALL 
####SBATCH --mail-user= 

module purge
module load gcc/11.4
module load openmpi 

srun ./kap8

