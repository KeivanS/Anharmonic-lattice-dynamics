Installation
============

Requirement
-----------

* Fortran Compiler (GCC or Intel)

Install
-------

.. _focex-install:

FOCEX
^^^^^^

* To install the FOCEX code, clone the repository `https://github.com/KeivanS/Anharmonic-lattice-dynamics.git` or download the ALADYN code and follow the instructions
* Edit the Makefile, generally you do not need to change the Makefile. In case if you need to change the compiler flags adjust the FLAGS option inside **Anharmonic-lattice-dynamics/FOCEX/Makefile** by uncommenting or adding the available option

``FLAGS= -fcheck=all -fdec-math #-Wall #-fbounds-check #-fbacktrace #-ffpe-trap=invalid,zero,overflow,underflow,inexact #-inline-debug-info #-p #-pedantic #-std=f95``

Navigate to the FOCEX directory and simply make the executable by using ``make all``. This will create executable ``focex.x`` for FOCEX code and other two executable ``readoutcar.x`` and ``readpwscf.x``. The latter two executables are from the utility folder for creating ``FORCEDISP`` file using VASP outcar and Quantum Espresso output log file. Executables ``focex.x``, ``readoutcar.x`` and ``readpwscf.x`` will be moved in your home directory inside ``aladyn`` folder. User can export the path by ``export PATH=${PATH}:~/aladyn`` so that it become available for the current shell session and put it in the terminal shell profile such as ``.bashrc`` or ``.bash_profile`` so that it is available system wide. 



.. _scop8-install:

SCOP8
^^^^^^

* On Linux

  * To install the SCOP8 code, clone the repository `https://github.com/KeivanS/Anharmonic-lattice-dynamics.git` or download the ALADYN code and follow the instructions
  * Edit the Makefile, generally you do not need to change the Makefile. In case if you need to change the compiler flags adjust the FLAGS option inside **Anharmonic-lattice-dynamics/SCOP8/Makefile** by uncommenting or adding the available option

  ``FLAGS= #-O3 #-C # to check everything O3 #g -p #for profiling with gprof``

  * From the **Anharmonic-lattice-dynamics/SCOP8/Makefile**, simply ``make`` and this will create the binary ``main.x``

* On Windows

  * Put MingW package under the root folder `download here <https://drive.google.com/file/d/1mdHpw7Eac_hwmtHLrHkKdj9zlLljesz8/view?pli=1>`_. The newest version of MingW may not work properly 
  * Open control panel > system advanced settings > add path of "root:\\MingW", see below
  .. image:: Ref/how.jpg
   :width: 1000
  * Install IDE code\:\:blocks, `download here <https://www.codeblocks.org/downloads/>`_.
  * Set up environment in code\:\:blocks: go to menu > settings > compiler > set **gnu fortran** as default and auto detect compilers (you may need to uncheck all the optional compiling options)
  * Load project in code\:\:blocks by double clicking *test.cbp*
  * Compile and run the code by pressing F9

THERMACOND
^^^^^^^^^^^
* To compile the THERMACOND code, clone the repository https://github.com/KeivanS/Anharmonic-lattice-dynamics.git or download the ALADYN code and follow the instructions
* It is enough to run ``make`` in the main directory, but a suitable Makefile must be present in the directory. THERMACOND needs a gfortran compiler. After compilation succeeds, a binary file ``kap7_sy_tet`` will be created in the main directory. There is no need for input from the terminal. It can be invoked simply as ``./ kap7_sy_te``.
* In order to compute collision matrices parallelly, there is a bash script file ``....`` that user can use to distribute k-mesh to many cores.

