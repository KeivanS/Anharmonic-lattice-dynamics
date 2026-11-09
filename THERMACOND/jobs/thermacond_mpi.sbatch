#!/bin/bash
# THERMACOND with MPI on a SLURM cluster: one MPI rank per core, nodes x ntasks-per-node ranks.
# Submit it in the folder that holds the THERMACOND inputs:   sbatch thermacond_mpi.sbatch
#
# Memory: every rank holds its own V3 arrays, about 24 bytes x nk x (3 x atoms)^3 (nk = n1 n2 n3 of
# latdyn.params line 1; e.g. 2 GB per rank for 6 atoms on 24x24x24). If the ranks do not fit in the node's
# memory, use thermacond_openmp.sbatch (one copy shared by all threads), or fewer ranks with several threads
# each (export TCOND_NUM_THREADS=<threads per rank>, and --cpus-per-task with fewer --ntasks-per-node).
#
# Modules: the compiler and MPI THERMACOND was built with (./configure writes the ones loaded when it ran).
# To find module names on any cluster:  module spider openmpi        (the versions there are)
#                                       module spider openmpi/5.0.10 (what has to be loaded before it)
#SBATCH --job-name=thermacond
#SBATCH --partition=@PARTITION@
#SBATCH --account=@ACCOUNT@
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=@CORES@
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

@MODULES@
cd "$SLURM_SUBMIT_DIR"
srun @TCOND_MPI@ > thermacond.out 2>&1
