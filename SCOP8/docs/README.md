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

- <>
