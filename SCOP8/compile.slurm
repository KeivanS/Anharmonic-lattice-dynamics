#!/bin/bash
#SBATCH -A elmgroup
#SBATCH --time=15:00:00   # 3 hours
#SBATCH --output=outputfile
#SBATCH --error=errorfile
#SBATCH --mem-per-cpu=6000
#SBATCH --ntasks-per-node=20
#SBATCH --partition standard  #-partition serial  ## parallel
#SBATCH -N 1 


#module load quantumespresso
#module load gnuplot

gfortran -c constants.f90
gfortran -c geometry.f95
gfortran -c force_constants.f90
gfortran -c zhegv.f

gfortran -c mods9.f90
gfortran -c MatrixDiagonalize.f90
gfortran -c modules_tetra.f90
gfortran -c extratools.f90

gfortran -c Structure_info.f95
gfortran -c DFT_force_constants.f95
gfortran -c kp_1d.f90

gfortran -c others3_nshells.f90

gfortran -c Fourier_force_constants.f95
gfortran -c Iteration_parameters.f95

gfortran -c VA_math.f95
gfortran -c Broyden.f95

gfortran -c check.f90
gfortran -c ConjugateGradient.f90
gfortran -c force_update.f90

ar rcs libAll.a *.o
gfortran -o test main.f95 libAll.a 
rm *.o *.mod

#./test
