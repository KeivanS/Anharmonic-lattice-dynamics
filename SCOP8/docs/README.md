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

#### How to Install/Compile Code on Windows(UPDATE: vscode)

1. Install Visual Studio Code (download [here](https://code.visualstudio.com/download)). This is the ultimate all-in-one code development environment for Windows. Not only it's compatible with and switch with one click among C++, Python, Java, Fortran,...etc, it can also be integrated with Git, MinGW bash, WSL for Linux subsystem, Jupyter Notebook, PDF Reader and more. I strongly recommend this as the one and only code editor for everyone.

2. Install the latest MinGW (download [here](https://sourceforge.net/projects/mingw/)). I suggest installing it on C drive. Notice that you may later need download extra package via 'mingw-get', it's located in C:\MinGW\bin\mingw-get.exe

3. Add MinGW to the system environment by:  Open control panel > system advanced settings > add path of "C:\MinGW", see below

<img src="img\how.JPG" style="zoom:50%;" />

4. Setup Fortran environment with VScode:
- Install 'Modern Fortran' extension (download [here](https://marketplace.visualstudio.com/items?itemName=fortran-lang.linter-gfortran)) or press `Ctrl+Shift+X` in VScode and search for 'Modern Fortran'. There are also some useful extensions such as 'TODO Highlight'
- Integrate msys bash(1.0) for compiling and running the code (assuming you installed MinGW in C:\MinGW):
  - Create a file `Run_MSYS.bat` in C:\MinGW\msys\1.0\, copy and paste the script below into the file
  ```
  @rem Do not use "echo off" to not affect any child calls.
  @SETLOCAL
  @SETLOCAL ENABLEEXTENSIONS

  :: Figure out where msys's root folder. If you want, you could just add the folder in the line
  :: below.
  @set MSYSROOT=
  @if "x%MSYSROOT%"=="x" @if exist "%~dp0msys.bat" @set MSYSROOT=%~dp0
  @if "x%MSYSROOT%"=="x" @if exist "%~dp0.msys-root" @set /P MSYSROOT=<%~dp0.msys-root
  @if "x%MSYSROOT%"=="x" (
  @echo Could not locate your mysys root folder.
  @set /P MSYSROOT=Location:
  )
  :: Read as MSYSROOT.trim()
  @if not "x%MSYSROOT%"=="x" (
  @for /f "tokens=* delims= " %%a in ("%MSYSROOT%") do @set MSYSROOT=%%a
  @for /f "useback tokens=*" %%a in ('%MSYSROOT%') do @set MSYSROOT=%%~a
  @if not "%MSYSROOT:~-1%"=="\" @set MSYSROOT=%MSYSROOT%\
  )
  :: Verify that root folder exists
  @if not exist "%MSYSROOT%" (
  @echo "%MSYSROOT%" is not a valid folder. Please check for .msys-root in %~dp0, or if you entered the path above, please rerun this script and select a valid folder.
  @exit /B 1
  ) else (
  @if not "%MSYSROOT%"=="%~dp0" @echo %MSYSROOT%>%~dp0.msys-root
  )

  :: Home Folder
  :: If you'd prefer the home directory set to your C:\Users\Username folder, uncomment the two lines
  :: below.
  @rem @if not exist "%HOME%" @set HOME=%HOMEDRIVE%%HOMEPATH%
  @rem @if not exist "%HOME%" @set HOME=%USERPROFILE%
  @if not exist "%HOME%" @if not "%MSYSROOT%"=="" @set HOME=%MSYSROOT%home\%USERNAME%
  @if not "x%WD%"=="x" @set WD=
  @set PLINK_PROTOCOL=ssh
  @if not exist "%WD%msys-1.0.dll" @set WD=%MSYSROOT%\bin\
  @set MSYSCON=sh.exe

  :: Default action, open msys and go to the current folder.
  @set OLDCD=%CD%
  @if not "x%OLDCD%"=="x" @set CURRCD=%CD%
  :: Get the current console ("OEM") codepage.
  @for /f %%i in ('"%MSYSROOT%bin\getcp.exe" -oem') do @set cp_oem=%%i
  :: Get the current GUI ("ANSI") codepage.
  @for /f %%i in ('"%MSYSROOT%bin\getcp.exe" -ansi') do @set cp_ansi=%%i
  :: Set the console codepage to match the GUI codepage.
  @chcp %cp_ansi% > nul < nul
  @if not "x%OLDCD%"=="x" (
  @"%MSYSROOT%bin\bash.exe" -l -i -c "cd \"$CURRCD\"; exec /bin/bash -rcfile ~/.bash_profile"
  ) else (
  @"%MSYSROOT%bin\bash.exe" -l
  )
  :: Store the error level returned by bash.
  @set ErrorLevel=%ErrorLevel%
  :: Restore the original console codepage.
  @chcp %cp_oem% > nul < nul
  :: If we had a current directory at the store of the script, go back to it.
  @if not "x%OLDCD%"=="x" chdir /D "%OLDCD%"

  :: quit script with the current error level.
  @exit /b %ErrorLevel%
  ```

  - Download `getcp.exe` [here](https://github.com/msysgit/msysgit/tree/master/mingw/bin) and put it into C:\MinGW\msys\1.0\bin

  - Open VScode and press `Ctral+Shift+P`, type 'setting' and find the 'Open User Setting(JSON)' to open the `settings.json` file, add following lines
  ```
  "terminal.external.windowsExec": "C:\\MinGW\\msys\\1.0\\bin\\sh.exe",
  ```
  also add a block in the terminal dropdown menu as shown below:
  <img src="img\json_settings.png" style="zoom:50%;" />

5. Setup for MPI in Windows follow the steps [here](https://abhila.sh/writing/3/mpi_instructions.html)

6. Compile and Run the code: 
  - Open 'Developer Command Prompt for VS 2022', direct to the SCOP8 folder on your local machine
  - Type `code .` NOTICE: this should be the only way to open VScode 
  - Press `Ctrl+Shift+\`` to open the terminal window (remember to select 'msys' in the dropdown menu as the default terminal is windows powershell) and direct to the SCOP8 folder on your Windows PC, then input following command to compiler, the order matters

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

7. Git can also be integrated into VScode for repository fetch and change commitments with 1-click.