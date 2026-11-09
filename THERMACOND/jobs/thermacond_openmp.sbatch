#!/bin/bash
# THERMACOND without MPI (the serial build, tcond-<version>-serial.x) on one node, multithreaded with OpenMP:
# one task with all the cores of the node as threads, and one copy of the V3 arrays for all of them
# (about 24 bytes x nk x (3 x atoms)^3, e.g. 2 GB for 6 atoms on 24x24x24).
# Submit it in the folder that holds the THERMACOND inputs:   sbatch thermacond_openmp.sbatch
#
# Modules: the compiler THERMACOND was built with (./configure writes the ones loaded when it ran).
# To find module names on any cluster:  module spider gcc            (the versions there are)
#                                       module spider gcc/11.4.0     (what has to be loaded before it)
#SBATCH --job-name=thermacond
#SBATCH --partition=@PARTITION@
#SBATCH --account=@ACCOUNT@
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=@CORES@
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

@MODULES@
cd "$SLURM_SUBMIT_DIR"
export TCOND_NUM_THREADS=$SLURM_CPUS_PER_TASK      # THERMACOND's threads (it does not read OMP_NUM_THREADS)
@TCOND_SERIAL@ > thermacond.out 2>&1
