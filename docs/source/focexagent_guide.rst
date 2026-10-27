.. _focex-agent-guide:

FOCEX agent: user's guide
=========================

``focex-agent.py`` (``FOCEX/AGENT``, version 1.2) runs the whole chain from a crystal structure
to force constants and phonons: VASP relaxation, SCF and Born charges, SC-SNAPS snapshots, one
VASP run per snapshot, FOCEX, plots and, if asked, THERMACOND. It is driven by one short
instruction file, checks every step, logs it in ``report.log`` and can be restarted at any time.
It was added to the suite after ALATDYN; it is not needed to run FOCEX by hand. This page is the
short version; :ref:`focex-agent` is the full manual.

Two modes
---------

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - ``MODE:``
     - input
     - runs
   * - ``from_scratch`` (default)
     - a ``POSCAR`` (any cell, any orientation)
     - VASP, SC-SNAPS, FOCEX
   * - ``extract_only``
     - a folder with ``POSCAR1``, ``FORCEDISP1`` (``DATA:``)
     - FOCEX only

Example 1: from a POSCAR (from_scratch)
---------------------------------------

A folder with two files. ``POSCAR``:

.. code-block:: text

   Si (diamond)
   1.0
     0.0000  2.7155  2.7155
     2.7155  0.0000  2.7155
     2.7155  2.7155  0.0000
   Si
   2
   Direct
     0.00 0.00 0.00
     0.25 0.25 0.25

``focex-agent-instructions.txt`` (only ``Project`` and ``RANKS`` are required; the template
``FOCEX/AGENT/focex-agent-instructions.txt`` lists every keyword with its default):

.. code-block:: text

   Project: Si
   RANKS: 1 1 1 0 0 0 0 0      # harmonic and cubic force constants
   FUNCTIONAL: PBEsol
   MPI: 16                     # MPI ranks for VASP
   THERMACOND_JOB: 3           # thermal conductivity after FOCEX (iterative solution)
   THERMACOND_MESH: 12 12 12

Set up (seconds, no VASP), then start the workflow:

.. code-block:: bash

   source ~/dftphonon-env/bin/activate
   python ~/Anharmonic-lattice-dynamics/FOCEX/AGENT/focex-agent.py --profile rivanna
   cd Si && sbatch run_workflow.sh        # workstation, no profile: nohup bash run_workflow.sh &
   tail -f report.log

Other keywords often used:

.. code-block:: text

   PROJECT_DIR: /scratch/me/Si   # where the calculations go (default: folder named after Project)
   RANGE: 6 2 1 0 0 0 0          # neighbour shells for ranks 2, 3, 4, ...
   SUPERCELL: LARGE              # SMALL (default, >= 10 A), LARGE (>= 14 A), or 9 integers
   SNAPSHOTS: 30                 # default 50
   MAGMOM: Cr=3.0
   DFPT: auto                    # Born charges: auto (if there is a gap) | on | off

Any file the agent writes can be amended (INCAR, tag by tag) or replaced by a block named after it:

.. code-block:: text

   INCAR_snaps
     ALGO = Normal
   INCAR_snaps_END

Example 2: existing forces (extract_only)
-----------------------------------------

The data folder holds ``POSCAR1`` and ``FORCEDISP1`` (as written by ``read_outcar.x`` or
``read_qe.x``), optionally the primitive ``POSCAR`` and ``dielectric.params``, in any orientation;
the agent rotates everything into the FOCEX frame.

.. code-block:: text

   Project: MgO
   MODE: extract_only
   DATA: ~/data/MgO_forces        # POSCAR1, FORCEDISP1, dielectric.params
   SPECIES: Mg O                  # if the POSCARs have no element line
   RANKS: 1 1 1 0 0 0 0 0

.. code-block:: bash

   python ~/Anharmonic-lattice-dynamics/FOCEX/AGENT/focex-agent.py
   cd MgO && bash run_workflow.sh        # FOCEX only: minutes

Results
-------

``report.log`` (one line per task: OK or FAILED, what was checked), then in the project folder
``01_relax``, ``02_scf`` (and ``dfpt``), ``03_snaps``, ``04_supercells`` (``FORCEDISP1``, and
``stress.params``: the stress of the undistorted supercell, used by FOCEX for the Huang rows and
the elastic constants), ``05_focex`` (force constants ``fc*.dat``, dispersion ``bs_freq.dat``,
``mech.dat``, ``thermo_QHA.dat`` and PDF plots) and ``06_thermacond`` (``kappa_vs_T.dat``).

When a run stops
----------------

.. code-block:: bash

   python ~/Anharmonic-lattice-dynamics/FOCEX/AGENT/diagnose.py Si

It is read-only: it writes ``Si/diagnosis/report.txt`` (what failed and why) and, when it can,
``fix_it.sh``, which changes the inputs and resubmits. Starting ``run_workflow.sh`` again also
continues from the failed task; finished tasks are not repeated.

Installation
------------

Needed: VASP and its POTCAR library (from_scratch only), MPI, Python 3.8 or newer, gfortran,
gnuplot, and `VASP-Flow <https://github.com/KeivanS/VASP-Flow>`_, which the agent uses for the VASP
inputs and batch scripts.

.. code-block:: bash

   cd ~
   git clone https://github.com/KeivanS/Anharmonic-lattice-dynamics.git
   (cd Anharmonic-lattice-dynamics && make focex sc_snaps)   # -> ~/ALADYN/bin
   git clone https://github.com/KeivanS/VASP-Flow.git        # elsewhere: export VASP_FLOW_DIR=...
   python3 -m venv ~/dftphonon-env
   ~/dftphonon-env/bin/pip install numpy spglib pymatgen matplotlib
   export VASP_POTCAR_DIR=/path/to/potpaw_PBE.54             # put it in ~/.bashrc

VASP and MPI come from a VASP-Flow profile (``--profile``: ``workstation``, ``slurm``,
``rivanna``, or your own copy in ``~/VASP-Flow/profiles``); without one, ``VASP_STD`` and
``MPI_LAUNCH`` are used. ``python ~/VASP-Flow/check_setup.py local`` (``slurm`` on a cluster)
checks the installation and prints the command that fixes each problem. Details:
:ref:`focex-agent`.
