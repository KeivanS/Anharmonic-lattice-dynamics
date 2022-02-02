## SCOP8 Readme

- A Bismuth example set of input files have been put in <input_sample> folder, notice that to run the program, input files should be put in the same folder as the code files

---

#### Code Files Explanation

- *main.f95* is the main program file, can modify according to comments on different type of calculations
- *Broyden.f95* is the root finding module that utilizes Broyden's method
- *check.f90* contains several test subroutines, including free energy landscape calculations
- *ConjugateGradient.f90* is the root finding module that utilizes conjugate gradient method, currently not in use
- *DFT_force_constants.f95* declare force constants objects, allocate them and read info from fc#.dat files
- *extratools.f90* is the legacy code contains various utility functions
- *force_constants.f90* legacy code that mostly from **FOCEX**(I think?)
- *force_update.f90* a module written in 2020 summer mostly contains updating scheme subroutines(which is wrong), the useful subroutines for now are those with ASR checking and ASR fixing
- *Fourier_force_constants.f95* is a small module that mainly used for Fourier transform force constants
- *geometry.f95* contains various constants, math operation interfaces, etc.
- *Iteration_parameters.f95* declare variational parameters, allocate and initialize them, then also read important parameters such as loop control, temperature, etc from files
- *kp_1d* legacy code that relates with k mesh, can probably be combined into other module
- *MatrixDiagonalize.f90* modified legacy code that performs matrix diagonalization
- *mods9.f90* legacy code that contains useful subroutines
- *modules_tetra.f90* a module for tetrahedron method k mesh generation
- *others3_nshells.f90* legacy code that contains useful subroutines
- *Structure_info.f95* declare atom and structure related object, allocate and initialize them, also read parameters from legacy input files lat_fc.dat, params.inp, etc
- *VA_math.f95* major module for variational approach, contains subroutines for free energy and its gradients calculation and many others

---

#### Input Files Explanation

- <fc2.dat>,<fc3.dat>,<fc4.dat> are force constants file extracted from runing **FOCEX**, values are labeled with atom number and directional indexes, stored as objects in this code
- <lat_fc.dat> contains some structure info of crystal and a complete atom list of the supercell, stored as objects in this code
- <params.born> has the born charge and electric constsant
- <params.inp> contains some structure info
- <params.phon> contains some phonon related info, such as k mesh size, etc
- <iteration_parameters.in> is <u>user defined</u> input file, what every line does is noted directly in the sample file
- <kpbs.in> is <u>user defined</u> input file. It is only for post process, gives the high symmetry k point path that user chooses to calculate phonon dispersion, etc
- <targetInitialize.dat> is both <u>user defined</u> input file and code's output file. Can be used to manually assign initial values to corresponding variational parameters. When a run finishes, it will be updated as the optimized variational parameters at this temperature.

---

#### Selected Output Files Explanation

- <convergence.dat> records the converging history of Broyden. Col 1 is the iteration number, Col 2 is the L1 norm threshold, Col 3 is the current free energy
- <Dispersion.dat> is the phonon band data. First column is the index of k point. Each following columns corresponds to a band
- <dos...dat> files start with dos are phonon DOS 
- <eigenvalues.dat> and <eigenvectors.dat> record the eigenvalues and eigenvectors of current run
- <GradientF.dat> contains all the current variational parameters value, all the corresponding free energy gradients and free energy, from first iteration to the last iteration.
- <output.txt> contains various info about the running
- <result.txt> contains few selected results such as final free energy, final volume, gruneisen, specific heat, etc at this temperature
- <targetInitialize.dat> is explained above
- <Y_square.dat> contains every yy value and the its gamma point approximation

---

#### How to Run Code on Windows

1. Put that MingW folder in C:\ (download [here](https://drive.google.com/file/d/1mdHpw7Eac_hwmtHLrHkKdj9zlLljesz8/view?usp=sharing), the newest version of MingW may not work properly) 
2. Open control panel > system advanced settings > add path of "C:\MingW", see below

<img src="img\how.jpg" style="zoom:50%;" />

3. Install codeblocks, download [here](www.codeblocks.org/downloads/)

4. Codeblocks > settings > compiler > set ==gnu fortran== as default and auto detect compilers (you may need to uncheck all the optional compiling options)
