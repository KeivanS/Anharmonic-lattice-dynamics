Download and installation
=========================

Requirements
------------

* A Fortran compiler — GNU ``gfortran`` (version 9 or newer) or Intel
  ``ifort``/``ifx``. FOCEX and ANFOMOD need nothing else.
* An MPI implementation (``mpif90``) — only for THERMACOND (``make thermacond-serial``
  builds it without MPI) and SCOP8.
* ``make``.
* Optional, for plotting the results: ``gnuplot``.

No external numerical library is required: the LAPACK routine used for the
diagonalisation (``zhegv``) is bundled with the source.

Getting the code
----------------

The source is on GitHub: https://github.com/KeivanS/Anharmonic-lattice-dynamics
(in a browser: *Code* → *Download ZIP*). From a terminal, clone the repository once (this also
brings SC-SNAPS, in ``SC_SNAPS/``):

.. code-block:: bash

   git clone https://github.com/KeivanS/Anharmonic-lattice-dynamics.git
   cd Anharmonic-lattice-dynamics
   ./configure --install    # this computer's settings (next section); ends with "Ready"
   make                        # = make all

To update it later, from inside that folder:

.. code-block:: bash

   git pull
   ./configure              # again: a pull may change what the Makefiles need
   make                        # rebuild

If you have edited files of your own, ``git pull`` keeps them; it stops, changing nothing, when the
update touches a file you edited. Then ``git stash``, ``git pull``, ``git stash pop`` puts your
edits aside, updates, and puts them back.

.. _configure:

Configuring: ``configure``
-----------------------------

Once, after cloning (and again after installing or moving compilers, MPI, VASP or Python):

.. code-block:: bash

   ./configure              # look, test, write make.local and ~/.aladyn.env
   ./configure --install    # the same, and one line in ~/.bash_profile (macOS) or ~/.bashrc
                            # (Linux; ~/.zshrc with zsh) that reads ~/.aladyn.env in every shell
   ./configure --build      # the same, then make

It looks at the computer (operating system, cores, memory, SLURM, environment modules), finds the
Fortran compiler (``gfortran``, else ``ifx`` or ``ifort``) and the MPI wrapper (``mpif90``, and the
compiler it calls), compiles and runs a small OpenMP program on two threads and an MPI program on two
ranks (on a cluster login node it only compiles), and writes

* ``make.local``: ``DESTDIR``, the Fortran compiler ``FF`` (FOCEX, FOCEX/utility, ANFOMOD, SC-SNAPS), the
  MPI wrapper ``MPIFC`` (THERMACOND), ``FLAGS``, and on macOS with Homebrew's Open MPI ``OMPI_FC`` and
  ``SDKROOT``. ``make.inc`` (in git: the defaults) reads it at its end, so its values win; it is not in
  git, so ``git pull`` never conflicts with it (the old one is kept as ``make.local.bak``). The
  top-level ``make`` passes these to every part.
* ``~/.aladyn.env``: ``PATH`` with ``DESTDIR``, ``OMPI_FC`` and ``SDKROOT`` (macOS), and the settings
  of the FOCEX agent (below), which it also checks; ``--no-agent`` skips them.

Options: ``--prefix DIR`` (``DESTDIR``), ``--fc``, ``--mpifc`` (the compilers to use), and the agent's
``--potcar``, ``--vasp``, ``--python``. On a cluster load the compiler and MPI modules first: there it
also drafts the VASP-Flow profile of the FOCEX agent, ``~/VASP-Flow/profiles/<cluster name>.json``
(``--profile-name`` for another name), from ``sinfo`` (default partition, cores per node),
``sacctmgr`` (account) and ``module avail`` (the loaded modules and the newest ``vasp`` module),
keeps an existing one, and makes it the default (``FOCEX_PROFILE`` in ``~/.aladyn.env``). Check it.

Building everything at once
---------------------------

One ``make`` in the top folder builds every part and installs the executables in
``DESTDIR``, set in the top-level ``make.inc``:

.. code-block:: makefile

   DESTDIR = ~/ALADYN/bin

   # ANFOMOD, FOCEX
   FF = gfortran

   # SCOP8, THERMACOND
   FC = mpif90
   F90 = mpif90

   FLAGS = -O2 -fcheck=all -ffree-line-length-0

Edit it if you want a different installation directory (or give ``make DESTDIR=...``);
FOCEX and THERMACOND set their own compiler options in their Makefiles. Then

.. code-block:: bash

   make                    # = make all: everything below except SCOP8
   make focex              # FOCEX (focex-9.9.x) and its utilities (read_outcar.x, ...)
   make sc_snaps           # SC-SNAPS (sc_snaps.x, poscar2xyz.py)
   make thermacond         # THERMACOND with MPI (tcond-2.2.x); thermacond-serial: without
   make anfomod            # ANFOMOD
   make python-tools       # links focex-agent, focex-diagnose, bandplot, allplots, kappaplots, shell-optimizer, aladyn-gui
   make scop8              # SCOP8, not part of "make": it does not compile at present
   make clean              # removes the object and module files

On macOS with Homebrew's Open MPI, ``configure`` puts ``OMPI_FC`` and ``SDKROOT`` into ``make.local``
(without it: ``export OMPI_FC=/opt/homebrew/bin/gfortran`` and ``export SDKROOT=$(xcrun --show-sdk-path)``
before ``make``; see :doc:`runthermacond`).
The sections below describe each part and how to build it on its own.

Finally add the installation directory to your ``PATH`` — put the line in your
``~/.bashrc`` or ``~/.zshrc`` so that it is set in every session:

.. code-block:: bash

   export PATH=${PATH}:~/ALADYN/bin

.. _focex-install:

FOCEX
^^^^^

.. code-block:: bash

   cd FOCEX
   make                      # -> ~/ALADYN/bin/focex-9.9.x

This compiles the sources listed in ``FOCEX/Makefile`` and produces a single
executable. Two details of that Makefile are worth knowing:

* the name of the executable is ``focex-<ver>.x``, with the version number ``ver`` set
  near the top (currently 9.9; 8.18 and earlier were installed as ``~/BIN/v<ver>``), and
* it is moved to ``DESTDIR``, ``~/ALADYN/bin`` by default (``make DESTDIR=...`` for another
  place; the top-level ``make`` passes the ``DESTDIR`` of ``make.inc``).

``FOCEX/Makefile`` also has its own ``FF`` and ``FLAGS`` variables. Since 2026-11-04 the
default is ``-O3 -fopenmp``: the SVDs of the fit run on OpenMP threads. ``make DEBUG=1``
gives the earlier ``-O2 -fcheck=all`` build with array-bound checks (slower; for finding
bugs). Run ``make clean`` when switching between the two. Both give the same force
constants and frequencies.

Threads (FOCEX and THERMACOND): the program chooses the number itself and writes it in the
log and in ``times.dat``. It does not use ``OMP_NUM_THREADS``.

* ``FOCEX_NUM_THREADS`` (``TCOND_NUM_THREADS``) if set;
* 1 when started by ``mpirun`` with more than one process;
* ``SLURM_CPUS_PER_TASK`` in a Slurm job (``#SBATCH --cpus-per-task=N``);
* otherwise all cores but about 20 %, so that the machine stays usable (8 of 10).

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
   * - ``read_outcar.x``
     - Extracts positions, forces and energy from a VASP ``OUTCAR`` in the
       current directory into ``pos-forc.dat``.
   * - ``read_qe.x``
     - The same for a Quantum Espresso (PWSCF) output file, whose name is read
       from standard input: ``echo pw.out | read_qe.x``.
   * - ``poscar2xyz.x``
     - Converts a POSCAR to ``.xyz`` for visualisation.
   * - ``poscar2focex.x``
     - POSCAR to ``structure.params``, ``POSCAR_focex`` and ``kpbs.params`` (below).

``poscar2focex.x POSCAR`` (options ``-o``, ``-p``, ``-k`` for the three output names,
``-s symprec`` in Angstrom, default 1e-3, ``-f`` to overwrite) finds the symmetry of the
crystal itself, with no external library: Niggli reduction, the operations of the
lattice, those of the crystal (with their fractional translations), the crystal system and
the Bravais lattice. The cell and the band-structure path follow Setyawan and Curtarolo,
Comput. Mater. Sci. 49, 299 (2010), the convention of pymatgen's ``KPathSetyawanCurtarolo``
and of VASP-Flow's k-path:

* ``structure.params``: line 1 a b c alpha beta gamma of the conventional cell, line 2 the
  primitive vectors in conventional units, line 3 ``1.0 0``, the atoms in conventional
  units; masses and names from the POSCAR, the other lines the defaults of
  ``FOCEX/example/Ge``, marked ``default, check``;
* ``POSCAR_focex``: the primitive cell in the frame FOCEX builds from line 1 (a along x,
  b in the xy plane); build the supercells from it;
* ``kpbs.params``: the Setyawan-Curtarolo path for the lattice and its sub-case (FCC,
  BCT1, ORCF2, RHL1, MCLC3, TRI1a, ...), in primitive reduced coordinates.

The output does not depend on how the input cell is written: any primitive basis,
orientation or origin gives the same files. For a crystal of lower symmetry than its
lattice (pyrite, wurtzite) the setting of the atoms is fixed by a canonical rule, and a
chiral crystal is never mirrored. Checked against pymatgen on 26 structures covering all
Setyawan-Curtarolo lattices and 24 of their 26 sub-cases (MCLC2 and MCLC4 are exact
equalities between the lattice parameters), and on 78 randomly rewritten inputs; the test
scripts and the results are in ``FOCEX/utility/tests/poscar2focex``.

``FOCEX/utility/poscar2focex.py`` (Python, numpy and spglib) is the earlier version; it is
kept because ``focex-agent.py`` uses its functions.

.. note::

   The snapshot generator ``sc_snaps.x`` is built from the ``SC_SNAPS`` folder, not here:
   ``make sc_snaps`` in the top folder (or the top-level ``make``) installs it in
   ``~/ALADYN/bin``; see :doc:`runscsnaps`. The folder is identical to the separate
   repository `github.com/KeivanS/SC_SNAPS <https://github.com/KeivanS/SC_SNAPS>`__.

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
    `download here <https://drive.google.com/file/d/1mdHpw7Eac_hwmtHLrHkKdj9zlLljesz8/view?pli=1>`__.
    The newest version of MinGW may not work properly.
  * Open *Control Panel > System advanced settings* and add the path of
    ``root:\MingW``, see below:

    .. image:: Ref/how.jpg
       :width: 1000

  * Install the IDE code\:\:blocks,
    `download here <https://www.codeblocks.org/downloads/>`__.
  * Set up the environment in code\:\:blocks: *menu > settings > compiler*, set
    **GNU Fortran** as default and auto-detect compilers (you may need to
    uncheck all the optional compiling options).
  * Load the project in code\:\:blocks by double clicking *test.cbp*.
  * Compile and run the code by pressing F9.

THERMACOND
^^^^^^^^^^

* Run ``make`` in the ``THERMACOND`` directory. THERMACOND needs a ``gfortran``
  compiler and MPI (``make serial``: no MPI). The binary is named after the version,
  ``tcond-2.2.x`` (serial: ``tcond-2.2-serial.x``), and installed in ``~/ALADYN/bin``. It takes no input from the
  terminal: ``tcond-2.2.x``, or ``mpirun -np N tcond-2.2.x`` (:doc:`runthermacond`).

.. _focex-configure:

Setting up the FOCEX agent: ``focex-configure.sh``
---------------------------------------------------

The FOCEX agent (:doc:`focexagent_guide`) also runs the DFT calculations, so it has to know where
VASP, MPI, the POTCAR library, VASP-Flow and a suitable Python are. ``configure`` (above) runs this
script and reads its file from ``~/.aladyn.env``; on its own:

.. code-block:: bash

   cd FOCEX/AGENT
   bash focex-configure.sh               # look and report only; writes ~/.focex-agent.env
   bash focex-configure.sh --install     # the same, and adds one line to your shell start-up file
                                         # that reads ~/.focex-agent.env in every terminal and job

It looks at, in this order: the computer (operating system, processor, cores, Slurm), the compilers
(``gfortran``, ``mpif90``; needed only to build FOCEX and SC-SNAPS), the FOCEX, SC-SNAPS and
``read_outcar.x`` binaries, Python with numpy, spglib and pymatgen, VASP-Flow, the POTCAR library,
``vasp_std``, MPI and ``gnuplot``. It then writes the settings to ``~/.focex-agent.env`` and checks
them by building a POTCAR with VASP-Flow's own functions and by starting ``focex-agent.py
--version``. It ends with ``Ready`` and the number of warnings. When the search finds the wrong
program or nothing, tell it where to look:

.. code-block:: text

   --potcar DIR     the POTCAR library (PBE): the folder with Si/, O/, Zn/, ... each holding a POTCAR
   --vasp FILE      the vasp_std executable
   --python FILE    the Python interpreter of the environment with numpy, spglib, pymatgen
   --env FILE       where to write the settings (default ~/.focex-agent.env)

Run it again after installing or moving any of these programs. FOCEX alone (without the agent)
needs none of this: ``make`` and the ``PATH`` line above are enough. The step-by-step installation
of the agent is in its manual, ``FOCEX/AGENT/focexagent_manual.pdf``.

Checking the installation
-------------------------

The quickest end-to-end test is the germanium example shipped with FOCEX; see
:ref:`focex-quickstart`. It exercises the symmetry analysis, the SVD fit and the
whole phonon/thermodynamics post-processing, and finishes in a few minutes on a
single core.
