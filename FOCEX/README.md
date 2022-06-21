#### This code, FOCEX (FOrce Constant EXtraction) included in ALADYN (Anharmonic LAttice DYNamics) employs the force constant calculation, 2nd, 3rd and 4th order to be latter used for other thermodynamical properties including phonon spectra, phonon lifetimes, scattering rates and thermal conductivity. The installation of FOCEX is simple and just requires the fortran compiler

- To install the FOCEX code, clone the repo "github_link_here" or download the ALADYN code and follow the instructions below:

cd Anharmonic-lattice-dynamics/FOCEX

- edit the Makefile. generally you do not need to change the Makefile. In case if you need to change the compiler flags adjust the FLAGS option by uncommenting or adding the available option in

FLAGS= -fcheck=all #-Wall #-fbounds-check #-fbacktrace #-ffpe-trap=invalid,zero,overflow,underflow,inexact #-inline-debug-info #-p #-pedantic #-std=f95

- From that directory make the file by

make

- This generates the binary `fc234-13`
