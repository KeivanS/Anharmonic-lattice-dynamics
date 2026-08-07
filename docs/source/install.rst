Installation
============

Requirements
------------

* A Fortran compiler — GNU ``gfortran`` (version 9 or newer) or Intel
  ``ifort``/``ifx``. FOCEX and ANFOMOD need nothing else.
* An MPI implementation (``mpif90``) — only for SCOP8 and THERMACOND.
* ``make``.
* Optional, for plotting the results: ``gnuplot``.

No external numerical library is required: the LAPACK routine used for the
diagonalisation (``zhegv``) is bundled with the source.

Getting the code
----------------

.. code-block:: bash

   git clone https://github.com/KeivanS/Anharmonic-lattice-dynamics.git
   cd Anharmonic-lattice-dynamics

Building everything at once
---------------------------

The top-level ``make.inc`` holds the compiler settings shared by the four codes:

.. code-block:: makefile

   DESTDIR = ~/ALADYN/bin

   # ANFOMOD, FOCEX
   FF = gfortran

   # SCOP8, THERMACOND
   FC = mpif90
   F90 = mpif90

   FLAGS = -O2 -fcheck=all -ffree-line-length-0

Edit it if you use a different compiler or want a different installation
directory, then

.. code-block:: bash

   make          # builds ANFOMOD, FOCEX, SCOP8 and THERMACOND
   make clean    # removes the object and module files

Finally add the installation directory to your ``PATH`` — put the line in your
``~/.bashrc`` or ``~/.zshrc`` so that it is set in every session:

.. code-block:: bash

   export PATH=${PATH}:~/ALADYN/bin

.. _focex-install:

FOCEX
^^^^^

.. code-block:: bash

   mkdir -p ~/BIN            # see the note below
   cd FOCEX
   make

This compiles the sources listed in ``FOCEX/Makefile`` and produces a single
executable. Two details of that Makefile are worth knowing:

* the name of the executable is set by the ``exe`` variable near the top
  (currently ``v16``, following the version number ``ver`` of the code), and
* the executable is moved to ``~/BIN``, **not** to the ``DESTDIR`` of
  ``make.inc``. Create that directory first, or change the ``mv`` line to
  install somewhere else.

``FOCEX/Makefile`` also has its own ``FF`` and ``FLAGS`` variables. The default

.. code-block:: makefile

   FF = gfortran
   FLAGS = -O2 -fcheck=all

is a good compromise between speed and safety. Useful alternatives, commented
out in the file, are ``-fbacktrace -g`` (for debugging) and
``-Wall -Wextra -fimplicit-none`` (for development).

The utility programs
^^^^^^^^^^^^^^^^^^^^

The tools that prepare the input of FOCEX and post-process DFT output are built
separately:

.. code-block:: bash

   mkdir -p ~/ALADYN/bin     # or whatever DESTDIR you set in make.inc
   cd FOCEX/utility
   make

This produces four executables and moves them to ``DESTDIR``:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Executable
     - Purpose
   * - ``sc_snaps.x``
     - Builds a supercell and generates thermally displaced snapshots
       (``poscar_000``, ``poscar_001``, …) from ``cell.inp``,
       ``supercell.inp`` and ``snaps.inp``.
   * - ``read_outcar.x``
     - Extracts positions, forces and energy from a VASP ``OUTCAR`` in the
       current directory into ``pos-forc.dat``.
   * - ``read_qe.x``
     - The same for a Quantum Espresso (PWSCF) output file, whose name is read
       from standard input: ``echo pw.out | read_qe.x``.
   * - ``poscar2xyz.x``
     - Converts a POSCAR to ``.xyz`` for visualisation.

.. warning::

   ``sc_snaps.x`` stops at run time with ``Recursive call to nonrecursive
   procedure 'findif'`` when it is built with ``-fcheck=all``, as the default
   ``make.inc`` above does: ``-fcheck=all`` implies ``-fcheck=recursion``, and
   the finite-difference helper does call itself through ``d1v``/``d2v``.
   Compiled without the checks (plain ``-O2``) it runs. To keep the checks, add
   ``-frecursive``:

   .. code-block:: bash

      gfortran -O2 -fcheck=all -frecursive -ffree-line-length-0 sc_snaps.f90 -o sc_snaps.x

   The standalone SC-SNAPS release (:doc:`runscsnaps`) builds with plain ``-O2``
   and is not affected.

.. _scop8-install:

SCOP8
^^^^^

* On Linux

  * Adjust the ``FLAGS`` option in ``SCOP8/Makefile`` if you need different
    compiler flags:

    .. code-block:: makefile

       FLAGS= #-O3 #-C # to check everything O3 #g -p #for profiling with gprof

  * Run ``make`` in ``SCOP8``; this creates the binary ``main.x``.

* On Windows

  * Put the MinGW package under the root folder,
    `download here <https://drive.google.com/file/d/1mdHpw7Eac_hwmtHLrHkKdj9zlLljesz8/view?pli=1>`_.
    The newest version of MinGW may not work properly.
  * Open *Control Panel > System advanced settings* and add the path of
    ``root:\MingW``, see below:

    .. image:: Ref/how.jpg
       :width: 1000

  * Install the IDE code\:\:blocks,
    `download here <https://www.codeblocks.org/downloads/>`_.
  * Set up the environment in code\:\:blocks: *menu > settings > compiler*, set
    **GNU Fortran** as default and auto-detect compilers (you may need to
    uncheck all the optional compiling options).
  * Load the project in code\:\:blocks by double clicking *test.cbp*.
  * Compile and run the code by pressing F9.

THERMACOND
^^^^^^^^^^

* Run ``make`` in the ``THERMACOND`` directory. THERMACOND needs a ``gfortran``
  compiler and MPI. After a successful compilation, a binary named after the
  version (e.g. ``kap8`` for version 8) is created in that directory; it takes
  no input from the terminal and is invoked simply as ``./kap8``.
* To compute the collision matrices in parallel, the bash script ``splitjob.sh``
  distributes the k-mesh over many cores.

Checking the installation
-------------------------

The quickest end-to-end test is the germanium example shipped with FOCEX; see
:ref:`focex-quickstart`. It exercises the symmetry analysis, the SVD fit and the
whole phonon/thermodynamics post-processing, and finishes in a few minutes on a
single core.
