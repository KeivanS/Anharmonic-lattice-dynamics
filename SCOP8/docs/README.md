## SCOP8 Readme

- A Bismuth example set of input files have been put in <input_sample> folder, notice that to run the program, input files should be put in the same folder as the code files

---

#### Code Files Explanation

- *test.cbp* is the code::blocks project file that contains dependency and so on, only for Windows

- *main.f95* is the main program file, can modify according to comments on different type of calculations
- *Broyden.f95* is the root finding module that utilizes Broyden's method
- *check.f90* contains several test subroutines, including free energy landscape calculations
- *ConjugateGradient.f90* is the root finding module that utilizes conjugate gradient method, currently not in use
- *DFT_force_constants.f95* declare force constants objects, allocate them and read info from fc#.dat files
- *extratools.f90* is the legacy code contains various utility functions
- *force_constants.f90* legacy code that mostly from **FOCEX**
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
  - The new and old math regarding whether to couple $\eta$ with $\lang yy\rang$ can be switched between subroutine 'GetV_avg_And_GradientV_avg2' and 'GetV_avg_And_GradientV_avg' 


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

#### Output Files Explanation

- <convergence.dat> records the converging history of Broyden. Col 1 is the iteration number, Col 2 is the L1 norm threshold, Col 3 is the current free energy
- <Dispersion.dat> is the phonon band data. First column is the index of k point. Each following columns corresponds to a band
- <dos...dat> files start with 'dos' are phonon DOS using different methods.
- <eigenvalues.dat> and <eigenvectors.dat> record the eigenvalues and eigenvectors of current run
- <GradientF.dat> contains all the current variational parameters value, all the corresponding free energy gradients and free energy, from first iteration to the last iteration.
- <output.txt> contains various info about the running
- <result.txt> contains few selected results such as final free energy, final volume, gruneisen, specific heat, etc at this temperature
- <targetInitialize.dat> is explained above, which have all the converged/optimized results for variational parameters $u_\tau$, $\eta$ and $K$, this file also serves as an input file for next run if 'inherit' option is selected
- <Y_square.dat> contains every yy value and the its gamma point approximation

---

#### How to Install/Compile Code on Windows

1. Put that MinGW folder in C:\ (download [here](https://drive.google.com/file/d/1mdHpw7Eac_hwmtHLrHkKdj9zlLljesz8/view?usp=sharing), the newest version of MingW may not work properly) 
2. Open control panel > system advanced settings > add path of "C:\MingW", see below

<img src="img\how.JPG" style="zoom:50%;" />

- For code editing:

  - Install code::blocks, download [here](www.codeblocks.org/downloads/)

  - Codeblocks > settings > compiler > set **gnu fortran** as default and auto detect compilers (this step is optional as we're not compile the code in code::blocks)

    <img src="img\how2.JPG" style="zoom:50%;" />

  - Double click the file *test.cbp* to open the project in code::blocks and modify code as wish

- For code compiling:

  - **Important preparations for MPI library installation and link**: follow [here](https://abhila.sh/writing/3/mpi_instructions.html)
  - Direct to the installation folder of MinGW, find msys -> 1.0 -> msys.bat, this is the MinGW shell terminal we will use for compile and execute the code
  - Open the terminal and direct to the SCOP8 folder on your Windows PC, then input following command to compiler, the order matters

  > gfortran -c constants.f90
  >
  > gfortran -c geometry.f95
  >
  > gfortran -c force_constants.f90
  >
  > gfortran -c zhegv.f
  >
  > gfortran -c mods9.f90
  >
  > gfortran -c MatrixDiagonalize.f90
  >
  > gfortran -c modules_tetra.f90
  >
  > gfortran -c extratools.f90
  >
  > gfortran -c Structure_info.f95
  >
  > gfortran -c DFT_force_constants.f95
  >
  > gfortran -c kp_1d.f90
  > 
  > gfortran -c others3_nshells.f90
  > 
  > gfortran -c Fourier_force_constants.f95
  > 
  > gfortran -c Iteration_parameters.f95
  > 
  > gfortran -c mpi_params.f95
  > 
  > gfortran -c Broyden.f95
  > 
  > gfortran -c VA_math.f95
  > 
  > gfortran -c force_update.f90
  > 
  > gfortran -c check.f90
  > 
  > gfortran -c ConjugateGradient.f90
  >
  > ar rcs libAll.a *.o
  >
  > gfortran -o out.exe main.f95 -IC:'\lib\mpi\include' -LC:'\lib\mpi\lib' -lmsmpi libAll.a libmsmpi.a

  - Execute use 'mpiexec', for example, if want to use 4 processes, input below

  > mpiexec -np 4 ./out.exe

---

