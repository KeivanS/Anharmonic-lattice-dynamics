#### This code, FOCEX (FOrce Constant EXtraction) included in ALATDYN (Anharmonic LATtice DYNamics) extracts the force constants (FCs) of ranks 2,3,4 .. up to 8th order for a force-displacement data of atoms in one or many supercells calculated using DFT (or any other method). From the FCs, one can calculate the phonon spectrum as well as other thermodynamic properties including elastic constants, total and free energy, vibrational entropy and heat capacity, Gruneisen parameter and thermal expansion coefficient...
Other codes in the suite, such as THERMACOND can in addition calculate phonon lifetimes, scattering rates and thermal conductivity. The installation of FOCEX is simple and just requires a fortran compiler (gfortran or intel).

- To install the FOCEX code, clone the repo "github_link_here" or download the ALATDYN code and follow the instructions below:

cd Anharmonic-lattice-dynamics/FOCEX

- edit the Makefile to select the compiling command and related flags, and the path to your BIN directory where the executable will be copied.  In case if you need to change the compiler flags adjust the FLAGS option by uncommenting or adding the available option

FLAGS= -fcheck=all 

- Type make to compile and copy the binary to your BIN folder.

- This generates the binary `fcx-$version.x` where $version is the version number that appears in the Makefile

- The provided *.plt files can be used to plot the results and analyze the data/run.
