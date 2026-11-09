.. _focex-agent:

FOCEX agent: from a crystal structure to phonons
================================================

Version 1.10 (2026-11-09, with FOCEX 9.9). The same manual as a PDF: ``FOCEX/AGENT/focexagent_manual.pdf`` in the repository; the agent, its
instruction file and ``diagnose.py`` are in ``FOCEX/AGENT``.

What the FOCEX agent does
-------------------------

FOCEX obtains the force constants of a crystal (harmonic, cubic, quartic, ...) by fitting them to
the forces that density-functional theory (DFT) gives on the atoms of a supercell in which every
atom has been displaced a little at random (a *snapshot*). From the force constants it computes the
phonon dispersion, density of states, group velocities, Grüneisen parameters, elastic constants,
sound speeds and the thermal expansion; THERMACOND then computes the lattice thermal conductivity
from the same files.

Making those forces takes many steps: relax the crystal, orient the cell the way FOCEX expects,
compute the Born charges of an insulator, choose a supercell, generate the snapshots, run one DFT
calculation per snapshot, collect the forces, write the five input files of FOCEX and run it. Each
step has its own input files and its own ways to fail. The FOCEX agent, ``focex-agent.py``, does
them all from one short *instruction file*, checks the result of every step, writes what it did to
``report.log``, and can be restarted at any time: finished steps are not repeated. When a run
stops, ``diagnose.py`` finds out where and why, and writes a script that repairs it
(:ref:`agent-diagnose`).

It works in one of two modes:

.. list-table::
   :header-rows: 1
   :widths: 16 40 44

   * - Mode
     - You give
     - The agent does
   * - **from_scratch**
     - a crystal structure (``POSCAR``) and the instruction file
     - VASP relaxation, SCF and Born charges; SC-SNAPS snapshots; one VASP run per snapshot;
       FOCEX; plots (:ref:`agent-scratch`)
   * - **extract_only**
     - existing supercells and their forces (``POSCAR1``, ``FORCEDISP1``, ...) and the
       instruction file
     - rotates the data into the FOCEX frame; FOCEX; plots (:ref:`agent-extract`)

.. figure:: _static/img/focexagent_flow.png
   :width: 90%
   :align: center

   Green: VASP; orange: SC-SNAPS, FOCEX and their tools; blue: your input.

What you need, and how to install it
------------------------------------

You do not need to have used VASP-Flow before: the agent uses it as a library, to write the VASP
input files and the batch scripts. This section has everything you need to set up.

The programs
~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 45 33

   * - What
     - Needed for
     - How to get it
   * - VASP 5.4 or newer
     - every DFT step (from_scratch only)
     - your licence; on a cluster usually a module
   * - POTCAR library
     - VASP: the PAW potentials, one folder per element
     - comes with VASP
   * - MPI
     - running VASP in parallel (``mpirun`` or ``srun``)
     - system or cluster module
   * - Python ≥ 3.8 with numpy, spglib, pymatgen
     - the agent and ``diagnose.py``
     - a virtual environment, below
   * - VASP-Flow
     - writes the VASP inputs and batch scripts
     - ``git clone``, below
   * - FOCEX, ``read_outcar.x``
     - the fit; reading the forces from the OUTCARs
     - this repository, ``make focex``
   * - SC-SNAPS (``sc_snaps.x``)
     - the snapshots (from_scratch only)
     - this repository, ``make sc_snaps``
   * - gnuplot
     - the plots (PDF)
     - ``brew install gnuplot``, ``apt install gnuplot``, or a module
   * - a Fortran compiler
     - building FOCEX and SC-SNAPS
     - ``gfortran`` (GCC)

Install, once
~~~~~~~~~~~~~

The same commands on a workstation and on the login node of a cluster. All the settings are made by
one command, ``./configure --install`` (step 4): nothing is edited by hand and nothing needs to be
activated afterwards. :ref:`configure` says what it does.

On a cluster, first, in this terminal only: ``module load gcc openmpi python`` (``module spider gcc``
lists the versions; ``configure`` saves the modules for every shell and job).

.. code-block:: bash

   cd ~
   # 1. the codes: FOCEX, SC-SNAPS, THERMACOND, the agent; and VASP-Flow, which the agent imports
   git clone https://github.com/KeivanS/Anharmonic-lattice-dynamics.git
   git clone https://github.com/KeivanS/VASP-Flow.git

   # 2. one Python environment for VASP-Flow and the FOCEX agent
   python3 -m venv ~/myenv
   ~/myenv/bin/pip install numpy spglib pymatgen matplotlib flask

   # 3. configure: finds and tests everything, writes all the settings; ends with "Ready"
   cd Anharmonic-lattice-dynamics
   ./configure --install                # on a cluster add: --account <your SLURM account>

   # 4. build (a few minutes) -> ~/ALADYN/bin
   make

   # 5. check, in a new terminal
   focex-agent --version                # focex-agent 1.10

The agent finds the programs at these places; a program installed elsewhere is named in the
instruction file (``FOCEX_X``, ``SC_SNAPS_X``, ``READ_OUTCAR_X``, :ref:`agent-keywords`). ``focex-agent`` and
``focex-diagnose`` in ``~/ALADYN/bin`` are launchers that start ``FOCEX/AGENT/focex-agent.py`` and
``diagnose.py`` with ``~/myenv``'s Python, so they work in any shell and batch job.

When step 3 reports a problem, the end of its report gives the command that fixes it; fix it and run
``./configure --install`` again. The usual ones:

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - message
     - what to do
   * - ``gfortran 8.5: version 9 or newer is needed``
     - cluster: ``module load gcc openmpi python``, then step 3 again; macOS: ``brew install gcc``
   * - ``~/myenv/bin/python does not start``
     - the module of its Python was not found: load it (``module spider python``), or make ``~/myenv``
       again with the Python loaded now (step 2)
   * - ``vasp_std cannot run here: N of its libraries are missing``
     - an old VASP build: on a cluster nothing to do (the VASP module is used); on a workstation rebuild
       VASP or give ``--vasp``
   * - ``account X is not one of yours`` / ``several SLURM accounts``
     - ``./configure --install --account <name>``
   * - ``module ... gives no vasp_std`` / ``misses libraries``
     - ``module spider vasp``, then ``--vasp-module <name>``
   * - ``no POTCAR library found``
     - ``--potcar /path/to/potpaw_PBE.54``
   * - ``the module gnuplot/... has no pdfcairo terminal``
     - nothing to do: the agent writes the plot scripts; run ``makeplots.sh`` (in ``05_focex``) on a
       computer whose gnuplot makes PDFs

Tell VASP-Flow where VASP, MPI and the POTCARs are
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``./configure`` does all of this; this section says what it set, for checking. **The POTCAR
library** is given by an environment variable, which ``configure`` finds and writes into
``~/.focex-agent.env`` (read by every shell and job):

.. code-block:: bash

   export VASP_POTCAR_DIR=/path/to/potpaw_PBE.54     # holds Zn/, O/, Ga/, ...

**VASP and MPI** come from a *profile*, a small JSON file in ``~/VASP-Flow/profiles/``, chosen with
``--profile`` (or ``PROFILE:`` in the instruction file):

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - profile
     - use it for
   * - (none)
     - a workstation. VASP from the variables ``VASP_STD`` (default ``~/BIN/vasp_std``) and
       ``MPI_LAUNCH`` (default ``mpirun -np``), set in ``~/.bashrc``; the number of ranks is
       ``MPI:`` in the instruction file; ``run_workflow.sh`` is started with ``bash``
   * - ``workstation``
     - the same, with the paths written in ``profiles/workstation.json`` instead
   * - ``slurm``
     - a template for a SLURM cluster: edit ``profiles/slurm.json`` (VASP command, ``modules`` to
       load, partition, nodes, tasks per node, time, account)
   * - *your own* (``mycluster`` below)
     - your cluster: ``./configure`` on the cluster drafts it (``<cluster name>.json``, from ``sinfo``,
       ``sacctmgr`` and ``module avail``, the VASP module tested; drafted again on a rerun, one written
       by hand is kept) and makes it the default (``FOCEX_PROFILE`` in ``~/.aladyn.env``); when the
       VASP module unloads the Python modules, ``mpi_cmd`` is ``~/ALADYN/bin/vasp-env srun``
       (:ref:`configure`). Or copy ``slurm.json`` to ``mycluster.json``, edit it, and use
       ``--profile mycluster``

With a SLURM profile, ``run_workflow.sh`` begins with the ``#SBATCH`` lines and the
``module load`` lines of the profile and is submitted with ``sbatch``. VASP runs on all tasks of the
job, nodes x ``ntasks_per_node`` (``NODES:`` and ``NTASKS_PER_NODE:`` in the instruction file replace
the profile's); KPAR and NCORE are chosen for that number. ``MPI: N`` fewer than the tasks is run as
``srun -n N``; more than the tasks stops the set-up.

Check the set-up
~~~~~~~~~~~~~~~~

.. code-block:: bash

   ~/myenv/bin/python ~/VASP-Flow/check_setup.py local     # a workstation;  "slurm" on a cluster

It checks Python and its packages, the POTCAR library, VASP and MPI (and on a cluster the
profile), and prints for every problem the command that fixes it. **Ready.** at the end means you
can start.

When is the Python environment needed?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is the same environment as in VASP-Flow's manual: one, ``~/myenv``, serves both (DFT with
VASP-Flow, phonons with the FOCEX agent). It never has to be activated: the launchers
``focex-agent``, ``focex-diagnose``, ``aladyn-gui``, ``allplots``, ``kappaplots`` start the scripts with it,
``~/.focex-agent.env`` puts its ``bin`` first in ``PATH``, and the scripts the agent writes
(``run_workflow.sh``, ``resubmit_failed.sh``, ``diagnosis/fix_it.sh``) contain the full path of the
interpreter. (On a cluster it needs the module of its Python, which ``~/.aladyn.env`` and the profile load.) The plots need only gnuplot; without it
the agent writes the plot scripts and you run ``makeplots.sh`` later, where gnuplot is.

Quick start
-----------

The steps of a from_scratch calculation (:ref:`agent-extract` has those of extract_only):

#. **Make a folder** for the material with its crystal structure, in VASP's POSCAR format, named
   ``POSCAR``. Any cell of the crystal will do -- primitive or not, any orientation: the agent
   reduces and orients it.
#. **Copy the instruction file** ``FOCEX/AGENT/focex-agent-instructions.txt`` into that folder and
   edit it. Only two lines are required: ``Project:`` (a name) and ``RANKS:`` (which force
   constants to fit, e.g. ``1 1 1 0 0 0 0 0`` for harmonic and cubic). Set also ``FUNCTIONAL:``
   and ``MPI:`` (the number of MPI ranks for VASP); every other line is optional and shows its
   default.

   **Where the calculations go** is also set in this file, by ``PROJECT_DIR:``: a path relative to
   the folder of the instruction file (``PROJECT_DIR: runs/ZnO_PBEsol``) or absolute
   (``PROJECT_DIR: /scratch/me/ZnO``, ``~`` allowed). Without it the project folder is a folder
   named after ``Project:`` next to the instruction file. The agent creates it, with any missing
   parent folders; if it holds a project already, the agent stops unless ``--force`` is given.
#. **Set up the project**: the agent creates the project folder with the relaxation inputs and
   ``run_workflow.sh``, and prints the command to start it. This takes seconds and runs no VASP.
#. **Start the workflow** from the project folder, and follow ``report.log``.
#. **If it stops**, run ``diagnose.py`` (:ref:`agent-diagnose`).

.. code-block:: bash

   cd ~/work/ZnO                      # holds POSCAR
   A=~/Anharmonic-lattice-dynamics/FOCEX/AGENT
   cp $A/focex-agent-instructions.txt .
   nano focex-agent-instructions.txt  # Project, RANKS, FUNCTIONAL, MPI

   focex-agent                        # on a cluster the profile of ./configure (or --profile NAME)
   cd ZnO
   sbatch run_workflow.sh             # without a profile: nohup bash run_workflow.sh &
   tail -f report.log

The options of ``focex-agent.py`` are all optional:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - option
     - meaning
   * - ``--instr FILE``
     - the instruction file, if it is not called ``focex-agent-instructions.txt``
   * - ``--coord FILE``
     - the crystal structure, if it is not called ``POSCAR`` (from_scratch)
   * - ``--profile NAME``
     - the VASP-Flow profile (previous section); default: ``PROFILE:`` of the instruction file,
       else none
   * - ``--force``
     - set up again a project folder that exists already. Its record of finished tasks is reset,
       so use it only to start over, after removing the old run folders

.. note::

   **Threads.** You do not have to do anything about OpenMP threads: the agent runs VASP with one
   thread per MPI rank (``OMP_NUM_THREADS=1``). Without that, a BLAS library built with OpenMP
   starts one thread per core in every MPI rank and the cores are oversubscribed many times (a
   96-atom ZnO supercell made no SCF step in 30 minutes on 8 ranks). Only if your VASP is built
   for hybrid MPI + OpenMP runs and you want them, set ``VASP_OMP_NUM_THREADS`` to the number of
   threads per rank before you start the workflow.

.. _agent-scratch:

Mode from_scratch: from a POSCAR to phonons
-------------------------------------------

``run_workflow.sh`` runs ``focex-agent.py run <PROJECT_DIR>``, which does the eleven tasks below in
order. Each finished task is recorded in ``state.json`` and in ``report.log``, with what was
checked; a task that fails stops the workflow, and starting ``run_workflow.sh`` again continues
from that task. A job that hits its time limit is simply submitted again.

The supercell tasks are checked again before they are skipped:

- **A removed or incomplete ``04_supercells/snap_NNN`` folder** is written again, and the
  ``supercells`` task and every task after it are redone. Only the folders without a converged
  OUTCAR run VASP.
- **To leave snapshots out on purpose**, move their folders to ``04_supercells/unused/``. They
  are then not recreated, and ``FORCEDISP1`` is collected again from the others.
- **A ``snap_NNN`` folder with no ``03_snaps/poscar_NNN``** stops the workflow.
- **To add snapshots**, raise ``SNAPSHOTS`` in the project's ``focex-agent-instructions.txt`` and
  start ``run_workflow.sh`` again. The extra snapshots are made at the same temperature and numbered
  after the existing ones, which are kept; only the new ones run VASP, then FOCEX is run again on
  all of them. A lower ``SNAPSHOTS`` removes nothing.

``relax``
   The structure you gave is reduced to its primitive cell (Niggli reduction) and rotated into
   the orientation FOCEX and SC-SNAPS use (:ref:`agent-frame`); the result is ``POSCAR.init``.
   VASP relaxes the atoms and the cell shape and volume in ``01_relax``. If the forces are not yet
   below EDIFFG after NSW ionic steps, the relaxation continues from its last structure (CONTCAR),
   up to three times.
``poscar0``
   The relaxed cell is reduced and oriented again and written as ``POSCAR0``: every later step
   uses this cell and this Cartesian frame. The report gives its space group and lattice
   parameters.
``scf``
   A self-consistent calculation on POSCAR0 in ``02_scf`` (same INCAR, no ionic steps, same
   k-mesh), which also gives the band gap. If the cell did not change in the previous task, it
   starts from the charge density and wave functions of the relaxation.
``dfpt``
   For an insulator (band gap above 0.05 eV, ``DFPT: auto``): density-functional perturbation
   theory in ``06_dfpt`` (as in VASP-Flow; ``02_scf/dfpt`` in projects made before 2026-11-07)
   gives the high-frequency dielectric constant and the Born effective charges of every atom
   (:ref:`agent-born`). ``DFPT: on`` forces it, ``off`` skips it. With the finite-field method
   (``LCALCEPS``) the field ``EFIELD_PEAD`` is set to half of VASP's limit for Zener tunnelling,
   ``e|E| |a_i| < E_g/(10 N_i)`` (E_g the SCF gap, N_i the k-points along a_i; at most VASP's default
   0.01 eV/Å), unless you set it in ``INCAR_dfpt``. The values are read from the ``OUTCAR``;
   ``extract_born.py`` also writes them to ``born_charges.txt`` and ``BORN`` (its output in
   ``extract_born.out``). The result is checked: no Zener warning from VASP, Born charges summing
   to zero (within 0.1 e per component), ε∞ between 1 and 100. After a Zener warning the run is
   repeated once with the smaller field. If a check still fails, ``report.log`` gives a WARNING and
   the charges are not used: FOCEX runs with ``born_flag`` 0, unless you give them in a
   ``dielectric_params`` block.
``snaps``
   The supercell is chosen (:ref:`agent-supercell`). SC-SNAPS (``03_snaps``)
   then writes the undistorted supercell ``poscar_000`` and 50 snapshots, each atom displaced as in
   a thermal ensemble of harmonic oscillators of frequency ``AVG_FREQ`` at temperature
   ``TEMPERATURE``. **A temperature you give is used as it is**, e.g. to sample the
   displacements of that temperature for effective phonons; the report gives the resulting rms
   displacement. Only when ``TEMPERATURE`` is not given, the default 30 K is adjusted (SC-SNAPS
   rerun at most twice) until the rms displacement is within 15 % of ``DISP_TARGET`` (0.02 Å):
   large enough for clean forces, small enough for the Taylor expansion of the energy. The
   undistorted supercell is checked to be an
   exact multiple of POSCAR0 before any VASP job is run.
``supercells``
   One folder per snapshot, ``04_supercells/snap_NNN``, with its POSCAR, the SCF INCAR without
   symmetry (``ISYM = 0``) and a k-mesh equivalent to the SCF one (:ref:`agent-kmesh`); magnetic
   moments given per atom are carried over to the supercell.
``vasp_snaps``
   The VASP runs of the snapshots, one after the other. A run that fails does not stop the others;
   the failed ones are listed in ``report.log``, and ``resubmit_failed.sh`` reruns only those. An
   existing OUTCAR is never overwritten (it is kept as ``OUTCAR.1``, ...).
``collect``
   ``read_outcar.x`` reads the positions and forces of every run into
   ``04_supercells/FORCEDISP1``; ``POSCAR1`` is the undistorted supercell. The report gives the
   largest net force on a supercell (should be ≈ 0) and the largest difference between the
   positions in OUTCAR and in POSCAR. The undistorted supercell (``snap_000``, run with
   ``ISIF = 2``; the other snapshots use ``ISIF = 1``, which gives only the trace) also gives the
   stress of the reference structure: ``stress.params`` (GPa, tensile positive, i.e. minus the
   ``in kB`` line of OUTCAR divided by 10) is written and copied to the FOCEX folders, where it
   enters the Huang rows and the elastic constants (:ref:`focex-stress-params`).
``focex``
   The five input files of FOCEX are written in ``05_focex`` (``structure.params``,
   ``latdyn.params``, ``kpbs.params``, ``dielectric.params``, ``default.params``) and FOCEX is run
   with the invariances held exactly (``enforce_inv = 1``). With Born charges a second run without
   them is made in ``05_focex/NOBORN``, for comparison; with fewer than ``LASSO_RATIO`` (2)
   equations per unknown, a third run with LASSO and a refit (``enforce_inv = 3``), which keeps only
   the force constants the data need, is made in ``05_focex/LASSO``. The report gives the fit
   error, the numbers of equations and unknowns, the lowest ω² (a warning below −25 cm⁻²:
   imaginary modes) and the highest frequency. If the DOS range of ``latdyn.params`` does not fit
   the frequencies, FOCEX is run again with one that does.
``plots``
   The gnuplot scripts ``bvel.plt``, ``all.plt``, ``qha.plt`` are filled in and run in every FOCEX
   folder (:ref:`agent-results`); ``makeplots.sh`` remakes them.
``thermacond``
   THERMACOND in ``07_thermacond``: linewidths or the thermal conductivity, when ``THERMACOND_JOB``
   is 1, 2 or 3; skipped for 0, the default (:ref:`agent-thermacond`).

Defaults
~~~~~~~~

- ``ACCURACY: Accurate``: ENCUT 600 eV, EDIFF 10⁻⁷ eV, EDIFFG −0.001 eV/Å, PREC Accurate in every
  VASP step. ``Normal``: 450 eV, 10⁻⁶, −0.003, PREC Normal -- cheaper, but the residual forces of
  the relaxed structure are larger (0.01 eV/Å in a ZnO test), which adds noise to the force
  constants. DFPT uses EDIFF 10⁻⁸.
- Relaxation and SCF k-mesh: VASP-Flow's ``coarse`` density (1000 k-points per reciprocal atom),
  ``fine`` (6000) if every element is a metal. Each subdivision is rounded to the nearest integer;
  VASP-Flow itself rounds to even numbers (multiples of 6 in plane for hexagonal cells), the agent
  does not. The FOCEX mesh (``FOCEX_MESH``) is made the same way.
- 50 snapshots. Each gives three equations per atom of the supercell: 20 snapshots of a 96-atom
  supercell give 5760 equations for the 264 unknowns of ZnO with ranks 1--3, so fewer snapshots are
  often enough (``SNAPSHOTS``).
- ``RANGE: 6 2 1 0 0 0 0``: neighbour shells kept for ranks 2, 3, 4, ...; ``FC2RANGE: 0``: the
  harmonic force constants are kept over the whole supercell.
- FOCEX k-mesh ``fine``; the band-structure path is found automatically.

.. _agent-supercell:

The supercell
~~~~~~~~~~~~~

``SUPERCELL:`` takes one of three forms:

- ``SMALL`` (the default): the multiple m1×m2×m3 of the conventional cell with the fewest atoms
  whose faces are at least 10 Å apart, at most 100 atoms;
- ``LARGE``: the same with 14 Å and at most 300 atoms. In both, if no such multiple is small
  enough, the one with the largest face distance; it keeps the symmetry of the lattice;
- **the supercell matrix itself**: nine integers, the three rows of the matrix that the agent
  writes to SC-SNAPS's ``supercell.inp``. The first three integers n11 n12 n13 give the first
  supercell vector as a combination of the primitive vectors of POSCAR0,
  R\ :sub:`S1` = n11 R\ :sub:`01` + n12 R\ :sub:`02` + n13 R\ :sub:`03`; the next three give
  R\ :sub:`S2` and the last three R\ :sub:`S3` in the same way (the agent's ``cell.inp`` has
  ``convcoord = 1``, primitive units). For example
  ``SUPERCELL: 2 0 0  0 2 0  0 0 2`` is 2×2×2 primitive cells, and for an fcc crystal
  ``SUPERCELL: -1 1 1  1 -1 1  1 1 -1`` is its conventional cube (4 primitive cells). The
  determinant must be positive.

The report gives the matrix and the number of atoms (for SMALL and LARGE also the face distance).

.. _agent-frame:

The FOCEX frame
~~~~~~~~~~~~~~~

FOCEX and SC-SNAPS build the conventional cell from a, b, c, α, β, γ with a along x and b in the
xy plane, and the supercells, forces and Born charges must be in that same Cartesian frame. The
agent rotates the cell into it at the start, so everything VASP computes is already in the right
frame. ``structure.params`` lists the conventional cell (line 1), the primitive vectors in its
units (line 2) and the atoms in conventional coordinates.

.. _agent-kmesh:

The k-mesh of the supercells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the supercell is a diagonal multiple of the primitive cell, its mesh is the SCF mesh divided by
the multiples, rounded up (6 6 4 with a 3×3×2 supercell gives 2 2 2). Otherwise (every multiple of
a centred conventional cell) it keeps the k-point spacing of the SCF mesh. ``KPOINTS_snaps`` sets
it yourself (:ref:`agent-override`).

.. _agent-extract:

Mode extract_only: FOCEX on existing data
-----------------------------------------

Use it when the supercell forces exist already -- from an earlier calculation, another code, or a
colleague -- and you want FOCEX's force constants and phonons. No VASP and no SC-SNAPS are run.

**The data folder** (named in the instruction file with ``DATA:``) holds

- ``POSCAR1`` and ``FORCEDISP1``: a supercell and the positions and forces of its snapshots, in
  FOCEX's format (what ``read_outcar.x`` or ``read_qe.x`` write); more pairs
  ``POSCAR2``/``FORCEDISP2``, ... are all used. Any orientation.
- optionally ``POSCAR``: the primitive cell, in the same orientation as the data; without it,
  spglib finds the primitive cell of ``POSCAR1``.
- optionally ``dielectric.params``: the dielectric constant and the Born charges, in the frame of
  the data, one 3×3 block per atom of that primitive cell. With more than one atom of a species the
  primitive ``POSCAR`` is required, to fix the order of the atoms.

``SPECIES: Mg O`` names the elements when the POSCARs have no element line.

**What the agent does.** At set-up it symmetrizes the lattice of the primitive cell (DFT lattice
parameters given to 6--7 digits make FOCEX's symmetry analysis lose force constants), finds the
rotation into the FOCEX frame, and rotates the supercells, the positions and forces in
``FORCEDISPi``, the dielectric constant and the Born charges with it. The report gives the rotation
and the size of the symmetrization. The workflow then has only the tasks focex and plots.

.. code-block:: bash

   cat > focex-agent-instructions.txt <<'EOF'
   Project: MgO
   MODE: extract_only
   DATA: ~/data/MgO_forces          # POSCAR1, FORCEDISP1, dielectric.params
   SPECIES: Mg O
   RANKS: 1 1 1 0 0 0 0 0
   EOF
   python ~/Anharmonic-lattice-dynamics/FOCEX/AGENT/focex-agent.py
   cd MgO && bash run_workflow.sh          # FOCEX takes minutes: no batch job needed

Checked on MgO data given in a rhombohedral orientation: the fit error (7.197 %), the 35 unknowns,
the frequencies and the elastic moduli are those of a FOCEX run in the original frame.

.. _agent-keywords:

The instruction file: every keyword
-----------------------------------

One keyword per line, ``KEYWORD: value``; ``#`` starts a comment. The VASP keywords are those of
VASP-Flow (its manual lists all of them); the common ones are below.

.. list-table::
   :header-rows: 1
   :widths: 20 21 59

   * - keyword
     - default
     - meaning
   * - **project**
     -
     -
   * - ``Project``
     - (required)
     - a name, used in titles and as the default folder
   * - ``PROJECT_DIR``
     - the Project name
     - folder of the calculations, relative to the instruction file or absolute; created by the
       agent
   * - ``MODE``
     - from_scratch
     - or extract_only
   * - ``DATA``
     - --
     - extract_only: the data folder
   * - ``SPECIES``
     - --
     - extract_only: element names, if the POSCARs have none
   * - ``PROFILE``
     - none
     - the VASP-Flow profile, as ``--profile``
   * - **force constants**
     -
     -
   * - ``RANKS``
     - (required)
     - 8 flags for ranks 1--8: ``1 1 0 0 0 0 0 0`` harmonic, ``1 1 1 0 0 0 0 0`` + cubic,
       ``1 1 1 1 0 0 0 0`` + quartic
   * - ``RANGE``
     - 6 2 1 0 0 0 0
     - neighbour shells for ranks 2--8
   * - ``FC2RANGE``
     - 0
     - 0: harmonic range = the supercell
   * - ``LASSO_RATIO``
     - 2
     - below this many equations per unknown, an extra LASSO run
   * - **VASP (VASP-Flow keywords)**
     -
     -
   * - ``FUNCTIONAL``
     - R2SCAN
     - PBE, PBEsol, LDA, AM05, R2SCAN, HSE06, VV10
   * - ``MPI``
     - 1; with a SLURM profile nodes x ntasks_per_node
     - MPI ranks for VASP (with a SLURM profile: fewer than the job's tasks runs as ``srun -n N``, more
       stops the set-up)
   * - ``ACCURACY``
     - Accurate
     - or Normal (:ref:`agent-scratch`)
   * - ``ENCUT``
     - by ACCURACY
     - plane-wave cutoff (eV)
   * - ``KMESH_DENSITY``
     - coarse or fine
     - k-mesh of the relaxation and SCF
   * - ``MAGMOM``
     - --
     - magnetic moments, e.g. ``MAGMOM: Cr=3.0``
   * - ``GGA_U``
     - auto
     - Hubbard U from VASP-Flow's table: ``auto`` applies it to d/f elements when the compound has
       O, S, Se, Te or a halogen; ``ON`` always; ``OFF`` never
   * - ``DFPT``
     - auto
     - Born charges: auto (if there is a gap), on, off
   * - ``BORN_TRANSPOSE``
     - yes
     - :ref:`agent-born`
   * - **snapshots**
     -
     -
   * - ``SUPERCELL``
     - SMALL
     - SMALL, LARGE, or the 9 integers of the supercell matrix, as in ``supercell.inp`` (rows,
       primitive units; :ref:`agent-supercell`)
   * - ``SNAPSHOTS``
     - 50
     - number of snapshots; can be raised later (:ref:`agent-scratch`)
   * - ``TEMPERATURE``
     - 30
     - K; given: used as is; not given: 30 K adjusted to DISP_TARGET
   * - ``AVG_FREQ``
     - 400
     - cm⁻¹, the frequency of the oscillators of SC-SNAPS
   * - ``DISP_TARGET``
     - 0.02
     - Å, the rms displacement wanted (only without TEMPERATURE)
   * - **FOCEX**
     -
     -
   * - ``FOCEX_MESH``
     - fine
     - k-mesh of latdyn.params: coarse, fine, a density, or 3 integers
   * - ``KPATH``
     - 0
     - 0: automatic band path; 1: give a ``kpbs_params`` block
   * - ``BS_NKPTS``
     - 40
     - points per segment of the band path
   * - **THERMACOND** (:ref:`agent-thermacond`)
     -
     -
   * - ``THERMACOND_JOB``
     - 0
     - 0 none; 1 linewidths; 2 kappa, direct; 3 kappa, iterative (latdyn line 8)
   * - ``THERMACOND_MESH``
     - the FOCEX mesh
     - n1 n2 n3 (line 1)
   * - ``THERMACOND_T``
     - 0 900 31
     - Tmin Tmax (K), number of temperatures (line 6)
   * - ``THERMACOND_SOLVER``
     - 0
     - job 2: 0 conjugate gradient, 1 SVD (line 9)
   * - ``VASP_RUN``
     - here
     - ``external``: the agent writes the inputs of each VASP stage and stops as WAITING; the runs are
       made elsewhere (``aladyn-gui``: bundle or submit, :doc:`aladyngui`) and the agent continues from
       their outputs
   * - ``THERMACOND_RUN``
     - here
     - the same for THERMACOND
   * - **programs**
     -
     -
   * - ``FOCEX_X``
     - newest ``~/ALADYN/bin/focex-*.x``
     - the FOCEX binary
   * - ``SC_SNAPS_X``
     - ``~/ALADYN/bin/sc_snaps.x``
     - SC-SNAPS
   * - ``READ_OUTCAR_X``
     - ``~/ALADYN/bin/read_outcar.x``
     - reads OUTCARs
   * - ``THERMACOND_X``
     - newest ``~/ALADYN/bin/tcond-*.x`` (MPI build; the serial one if there is none)
     - THERMACOND; the serial build (``tcond-*-serial.x``) runs as one process on the job's cores as
       threads (no ``srun``), see :doc:`runthermacond`

.. _agent-override:

Changing the files the agent writes
-----------------------------------

Every file the agent generates can be changed *from the instruction file*, with a block named after
the file: a line with its name, the contents, and a line with the name followed by ``_END``.

- **INCAR blocks are merged tag by tag**: the tags you give replace or add to those the agent
  writes; the others stay. ``INCAR`` applies to every VASP step; ``INCAR_relax``, ``INCAR_scf``,
  ``INCAR_dfpt``, ``INCAR_snaps`` to one.
- **Magnetic order and +U are the same in every step**: ``MAGMOM``, ``ISPIN``, ``NUPDOWN``,
  ``LNONCOLLINEAR``, ``LSORBIT``, ``SAXIS``, ``M_CONSTR``, ``I_CONSTRAINED_M``, ``LAMBDA``, ``LDAU``,
  ``LDAUTYPE``, ``LDAUL``, ``LDAUU``, ``LDAUJ`` and ``LMAXMIX``, given in any INCAR block (``INCAR_scf``
  too), go to the relaxation, the SCF, DFPT and the supercells; the same tag with different values in
  two blocks stops the agent. A ``MAGMOM`` (or ``M_CONSTR``) with one value, or three for
  non-collinear moments, per atom of the cell is written for the supercells with each atom's value
  taken from its site in the cell (an antiferromagnet keeps its order); a list of the supercell's
  length is used as it is, any other length stops the supercells task. ``LDAUL/U/J`` are per species
  and the supercells have the species of the cell. Give the per-atom list in the order of your
  POSCAR: the agent keeps that order. ``04_supercells/snapshots_mag.xyz`` shows the moments of every
  supercell (Jmol's xyz format, the moment as the vector of each atom, collinear ones along z; in Jmol
  ``vectors on``); the POSCARs are left as VASP reads them.
- **Every other block replaces the whole file**, so it must be complete: ``KPOINTS`` (relaxation
  and SCF), ``KPOINTS_snaps``, ``cell_inp``, ``supercell_inp``, ``snaps_inp`` (SC-SNAPS),
  ``structure_params``, ``latdyn_params``, ``kpbs_params``, ``dielectric_params``,
  ``default_params`` (FOCEX).

Simple choices -- the supercell, the number of snapshots, the FOCEX k-mesh -- have their own
keyword; use a block only for what has none.

**Example.** Wurtzite ZnO without Hubbard U (VASP-Flow would otherwise put its tabulated U on Zn,
an oxide), ``LASPH`` in every VASP step, stricter electronic convergence only in the supercell
runs, a 3×3×2 supercell, 20 snapshots, a denser FOCEX mesh, and your own band path
Γ--M--K--Γ--A, which replaces ``kpbs.params``:

.. code-block:: text

   Project: ZnO
   RANKS: 1 1 1 0 0 0 0 0
   FUNCTIONAL: PBEsol
   GGA_U: OFF
   MPI: 40
   SUPERCELL: 3 0 0  0 3 0  0 0 2       # rows: supercell vectors in primitive units
   SNAPSHOTS: 20
   FOCEX_MESH: 12 12 8
   KPATH: 1                              # the band path is given below

   INCAR                                 # every VASP step (merged tag by tag)
     LASPH = .TRUE.
   INCAR_END
   INCAR_snaps                           # the supercell runs only
     NELM = 200
     ALGO = All
   INCAR_snaps_END
   kpbs_params                           # replaces kpbs.params (FOCEX's band path)
   1        # coordinates in units of the primitive reciprocal vectors
   60       # k-points per segment
   4        # number of segments
   G   0.000000 0.000000 0.000000
   M   0.500000 0.000000 0.000000
   K   0.333333 0.333333 0.000000
   G   0.000000 0.000000 0.000000
   A   0.000000 0.000000 0.500000
   kpbs_params_END

After the set-up you can read ``01_relax/INCAR`` before anything runs.

.. _agent-born:

Born charges and the dielectric constant
----------------------------------------

In an insulator, the long-range Coulomb interaction between the displaced ions splits the
longitudinal and transverse optical modes near Γ (LO--TO splitting). FOCEX includes it with the
Born effective charges Z* and the high-frequency dielectric constant ε∞; without them the optical
branches near Γ are wrong.

- from_scratch: the agent computes them (task dfpt) when the SCF shows a band gap, writes them to
  ``dielectric.params``, and makes for comparison a FOCEX run without them in
  ``05_focex/NOBORN``. ``DFPT: off`` skips them; ``on`` computes them for any material.
- extract_only: give them in ``dielectric.params`` in the data folder.
- your own values, in either mode: a ``dielectric_params`` block (it replaces the file).

.. note::

   VASP prints Z*\ :sub:`ij` with i the direction of the field and j that of the displacement;
   FOCEX expects the transpose, which the agent writes. ``BORN_TRANSPOSE: no`` writes them as VASP
   prints them; you should not need it.

.. _agent-results:

The project folder and the results
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 37 63

   * - file or folder
     - what it is
   * - ``focex-agent-instructions.txt``
     - the copy of your instruction file that the runs read
   * - ``report.log``
     - one entry per task: OK or FAILED, minutes, what was checked -- read this first
   * - ``state.json``
     - the finished tasks and the numbers passed between them
   * - ``run_workflow.sh``
     - the job; ``resubmit_failed.sh`` reruns the failed supercell runs
   * - ``POSCAR.input``, ``POSCAR.init``, ``POSCAR0``
     - your structure; reduced and oriented; relaxed
   * - ``01_relax``, ``02_scf``
     - VASP: relaxation, SCF
   * - ``03_snaps``
     - SC-SNAPS input and the snapshots ``poscar_NNN``
   * - ``04_supercells``
     - ``snap_NNN/`` (one VASP run each), ``POSCAR1``, ``FORCEDISP1``
   * - ``05_focex``
     - FOCEX: inputs, outputs, plots; ``NOBORN/`` and ``LASSO/`` when made
   * - ``06_dfpt``
     - VASP: Born charges and ε∞ (``born_charges.txt``), as in VASP-Flow; ``02_scf/dfpt`` in
       projects made before 2026-11-07
   * - ``07_thermacond``
     - THERMACOND (``THERMACOND_JOB`` 1-3): ``kappa_vs_T.dat``, ``kappa_wigner.dat``, ...
   * - ``diagnosis``
     - written by ``diagnose.py``: ``report.txt``, ``fix_it.sh``, backups

In ``05_focex``: the force constants ``fc2.dat``, ``fc3.dat``, ... with the cell in ``lat_fc.dat``;
the dispersion ``bs_freq.dat``; the elastic constants and moduli ``mech.dat``; thermodynamics and
thermal expansion ``thermo_QHA.dat``; and three PDFs: ``<Project>.pdf`` (the dispersion coloured by
group velocity, and the group velocities against frequency), ``all-<Project>.pdf`` (Grüneisen
parameters, density of states, dispersion) and ``qha.pdf`` (thermal expansion, Grüneisen
parameter, bulk modulus, thermodynamic functions, free energy against volume). :doc:`runfocex`
describes every output file. The large VASP files (WAVECAR, CHGCAR, CHG, ...) are deleted when no
longer needed; the OUTCARs never.

.. _agent-diagnose:

When a run stops: diagnose.py
-----------------------------

.. code-block:: bash

   A=~/Anharmonic-lattice-dynamics/FOCEX/AGENT
   python $A/diagnose.py ~/work/ZnO/ZnO     # one or more project folders
   less ~/work/ZnO/ZnO/diagnosis/report.txt
   bash ~/work/ZnO/ZnO/diagnosis/fix_it.sh

``diagnose.py`` itself changes nothing. It lists every task as done, FAILED, INTERRUPTED (the job
was killed inside it) or not reached, or reports that the project is still running or not yet
started. For a failed task it finds the cause from ``report.log``, the OUTCARs, ``vasp.out``,
``focex.out`` and the SLURM logs, and writes ``diagnosis/report.txt`` (what failed, where, the
evidence, what the fix will change) and ``diagnosis/fix_it.sh``. VASP failures are recognised with
the rules of VASP-Flow's diagnosis tool; FOCEX, SC-SNAPS and ``read_outcar.x`` stops by their
messages.

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - cause
     - meaning
     - what fix_it.sh does
   * - ``timeout``
     - the job hit its time limit
     - walltime × 2, resubmit
   * - ``scf_nelm``
     - SCF not converged in NELM steps
     - ALGO = All, NELM = 200
   * - ``edddav``
     - the Davidson eigensolver failed
     - ALGO = All
   * - ``zbrent``
     - the relaxation's line search failed
     - IBRION = 1
   * - ``relax_unconv``
     - relaxation not converged
     - resubmit (continues from CONTCAR)
   * - ``oom``
     - out of memory
     - KPAR = 1
   * - ``mpi_init``, ``not_run``, ``interrupted``
     - MPI did not start; never ran; killed
     - resubmit
   * - others
     - disk quota, symmetry, an unknown stop, ...
     - held, with the error lines

New INCAR tags go into the INCARs of the failed runs and into the matching block of the project's
instruction file (``INCAR_snaps``, ...), so that later runs use them too. Before it changes
anything, ``fix_it.sh`` checks that nothing is running and that nothing changed since the
diagnosis, and copies every file it edits to ``diagnosis/backup_<date>/``; then it resubmits
(``sbatch`` with a SLURM profile, otherwise ``bash`` in the background). When the fix needs your
judgement -- a cause it does not know, a fix that was already tried, a new supercell -- the script
is *held*: it prints what to do and stops, and ``bash fix_it.sh --anyway`` runs its suggested
commands as they are. Run ``diagnose.py`` again when the resubmitted job has finished.

.. _agent-thermacond:

Thermal conductivity
--------------------

The last task, ``thermacond``, runs THERMACOND on the force constants of ``05_focex``. Whether it
runs is decided by the THERMACOND job, line 8 of ``latdyn.params``, which you set in the
instruction file:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - ``THERMACOND_JOB``
     - what THERMACOND does
   * - 0 (default)
     - nothing: the task is skipped (FOCEX has already given the bands, DOS and Grüneisen
       parameters)
   * - 1
     - three-phonon linewidths of every mode along the band path (``lw-selfenergy.dat``)
   * - 2
     - thermal conductivity: RTA and the direct solution of the BTE
   * - 3
     - thermal conductivity: RTA and the iterative solution of the BTE

.. code-block:: text

   RANKS: 1 1 1 0 0 0 0 0        # cubic force constants are needed
   THERMACOND_JOB: 2
   THERMACOND_MESH: 12 12 12     # latdyn.params line 1 for THERMACOND (default: the FOCEX mesh)
   THERMACOND_T: 100 900 9       # line 6: Tmin, Tmax (K), number of temperatures
   # THERMACOND_SOLVER: 0        # line 9, job 2: 0 conjugate gradient, 1 SVD

The task copies ``lat_fc.dat``, ``fc2.dat``, ``fc3.dat``, ``structure.params``,
``dielectric.params`` and ``kpbs.params`` from ``05_focex`` to ``07_thermacond``, writes there
FOCEX's ``latdyn.params`` with lines 1, 6, 8 and 9 from these keywords, and runs ``tcond-2.2.x``
with the MPI command of the profile; ``05_focex`` is not changed. A ``latdyn_params`` block
(:ref:`agent-override`) sets all 13 lines for both FOCEX and THERMACOND, and its line 8 is the job
unless ``THERMACOND_JOB`` is given. The mesh is the expensive choice: the FOCEX default (``fine``)
is usually far too dense, so set ``THERMACOND_MESH`` and converge kappa with it.

The report gives, at up to seven of the temperatures, kappa_RTA, the full (direct or iterative)
solution and the coherences term kappa_C of the Wigner formulation (W/mK, average of the
diagonal), and a warning if the full solution falls below the RTA, which it should not do (the
direct solver has an open issue: compare with job 3). The complete results are in
``07_thermacond``: ``kappa_vs_T.dat``, ``kappa_wigner.dat``, ``Qvalue.dat`` (scattering rates) and
the other files described in :doc:`runthermacond`. An existing project gets the task by adding
``THERMACOND_JOB`` to its ``focex-agent-instructions.txt`` and submitting ``run_workflow.sh``
again: the finished tasks are skipped. Changing ``THERMACOND_JOB``, ``_MESH``, ``_T`` or
``_SOLVER`` later and submitting again runs the task again; V3*delta is read back from
``07_thermacond/V3DELTA`` when the force constants and the mesh are the same (job 3 after job 2,
new temperatures).

Version history
---------------

1.10
   2026-11-09: installation in one command, ``./configure --install``, with a table of the problems it
   reports: on a cluster it finds and loads the modules ``~/myenv`` needs and gnuplot, writes the
   ``module load`` line into ``~/.aladyn.env``, takes the account from ``--account``, drafts the profile again
   on a rerun, tests the VASP module and, when it unloads the Python modules, runs VASP through
   ``~/ALADYN/bin/vasp-env``; a ``vasp_std`` with missing libraries is rejected; ``make python-tools`` writes
   launchers with ``~/myenv``'s Python (no activation); ``sinfo``'s ``40+`` cores no longer break the
   profile's JSON; a gnuplot without the ``pdfcairo`` terminal counts as no gnuplot (scripts written, the run
   goes on). Fixed: ``FOCEX_PROFILE`` replaced the ``--profile`` of the set-up at run time (1.8, 1.9).

1.9
   2026-11-09: the serial THERMACOND build runs as one process on threads (``TCOND_NUM_THREADS``), not with
   ``srun``, also in aladyn-gui's bundles (one task, ``--cpus-per-task``); without an MPI build the agent takes
   the serial one; ``run_workflow.sh`` says where its module lines come from; ``./configure`` writes THERMACOND
   job examples (MPI, OpenMP) to ``~/ALADYN/jobs``.

1.8
   2026-11-09: with a SLURM profile VASP runs on the job's tasks (nodes x ntasks_per_node; before, without
   an ``MPI`` line, on one core), ``MPI: N`` below them as ``srun -n N``, above them an error; ``NODES:``,
   ``NTASKS_PER_NODE:``, ``WALLTIME:`` reach the job header; ``FOCEX_PROFILE``; ``./configure`` drafts the
   profile on a cluster; the installation block follows ``./configure``.

1.7
   2026-11-09: MAGMOM, ISPIN, LDAU and the other magnetic and +U tags of any INCAR block apply to every
   step; per-atom MAGMOM / M_CONSTR carried to the supercell atoms through their sites (antiferromagnets),
   shown in ``04_supercells/snapshots_mag.xyz``; THERMACOND 2.2 plots after the thermacond task.

1.6
   2026-11-07: ``VASP_RUN: external`` and ``THERMACOND_RUN: external`` (inputs written, the agent stops as
   WAITING and continues from the outputs; used by ``aladyn-gui``); a ``thermacond_params`` block (the
   latdyn.params of ``07_thermacond``, used as it is); the SCF reads the gap also when its OUTCAR came from
   elsewhere, and reuses the relaxation's CHGCAR only when it is there.

1.5
   2026-11-07: DFPT in its own folder ``06_dfpt``, as in VASP-Flow, and THERMACOND in
   ``07_thermacond`` (older projects keep ``02_scf/dfpt`` and ``06_thermacond``); ``EFIELD_PEAD``
   from the gap and the k-mesh, checks of the Born charges and ε∞, ``born_charges.txt`` written
   (task dfpt). A Cu2O run had given ε∞ = 207 and charges summing to -3.2 e: VASP's default field
   was 2.7 times its own Zener limit.

1.4
   2026-11-07: the FOCEX fits hold the invariances exactly (``enforce_inv = 1``; was 0) and the
   extra LASSO run uses the refit (``enforce_inv = 3``; was 2), as the FOCEX FAQ recommends;
   ``diagnose.py`` reports the same version.

1.3
   2026-11-06: ``focex-configure.sh`` finds the programs, the POTCAR library, VASP and MPI and
   writes ``~/.focex-agent.env``; installation and quick start as numbered steps; the minimal
   instruction file and a table of every default and how to change it; without ``MPI`` VASP runs
   on one core; a diagram of the project folder (PDF manual).

1.2
   2026-10-27, with FOCEX 9.0: the undistorted supercell ``snap_000`` is run with ``ISIF = 2``, and
   ``collect`` writes its stress, minus the ``in kB`` line of OUTCAR divided by 10, as
   ``stress.params`` (GPa, tensile positive) for FOCEX (Huang rows and elastic constants); the test
   VASP stand-in prints the stress.

1.1
   2026-10-23: task ``thermacond``, switched on by ``THERMACOND_JOB`` (latdyn.params line 8), with
   ``THERMACOND_MESH``, ``THERMACOND_T``, ``THERMACOND_SOLVER``; a given ``TEMPERATURE`` is kept;
   ``SUPERCELL`` matrix and ``PROJECT_DIR`` documented; SC-SNAPS periodic images fixed.

1.0
   2026-10-23, with FOCEX 8.19: modes from_scratch and extract_only; instruction file
   ``focex-agent-instructions.txt``, ``--instr``, ``--coord``, ``PROJECT_DIR``; ``diagnose.py``;
   this manual.
