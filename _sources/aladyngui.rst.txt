.. _aladyn-gui:

ALADYN GUI: the whole calculation in a browser
==============================================

``aladyn-gui`` (``GUI/aladyn-gui.py``, version 0.1, 2026-11-07) sets up and follows a calculation from a
POSCAR to phonons and the thermal conductivity: VASP relaxation, SCF and Born charges, SC-SNAPS snapshots,
the supercell VASP runs, FOCEX and THERMACOND. Underneath it is the FOCEX agent (:doc:`runfocexagent`):
each project is an agent project folder, and the GUI writes its instruction file and runs its tasks.

The user's guide, step by step and with every button: ``GUI/aladyngui_manual.pdf`` in the repository.

Starting it
-----------

.. code-block:: bash

   make python-tools                 # once: links aladyn-gui into ~/ALADYN/bin
   cd ~/calc && aladyn-gui           # projects go in the current folder (or --workdir DIR)

and open http://localhost:5070. It needs the agent's Python environment (numpy, spglib, pymatgen,
VASP-Flow) and Flask; the POTCAR library is found as for the agent (``VASP_POTCAR_DIR``).

On a cluster, start it on the login node and reach it from your computer through an ssh tunnel:

.. code-block:: bash

   ssh -L 5070:localhost:5070 you@cluster      # then, on the cluster: cd /scratch/you && aladyn-gui
                                               # and open http://localhost:5070 on your computer

The guide
---------

The left panel lists the steps in order: Project, Relaxation, SCF and DFPT, Snapshots, Supercell DFT,
FOCEX, THERMACOND. Each is marked done, next, waiting (inputs ready, runs outstanding), failed or later.
Below it, *What to do next* says which button to press. The GUI stops before FOCEX and before THERMACOND,
so that their input files can be reviewed first.

Where each DFT stage runs
-------------------------

On the *Run* tab:

* **Continue** does the quick steps here (POSCAR0, snapshots, supercell folders, reading results that are
  back) and stops at the next DFT stage with its inputs ready. It never starts VASP itself.
* **Run here**: VASP runs on this computer (the relaxation of a small cell, for example).
* **Make bundle**: a ``.tar.gz`` in the project's ``bundles/`` folder (download link on the page). On the
  cluster:

  .. code-block:: bash

     tar xzf Cu2O_vasp_snaps_<date>.tar.gz
     cd Cu2O && ./submit_all.sh       # all supercells at once; DFPT after the SCF; collect.py at the end
     squeue -u $USER                  # wait until the jobs are gone
     ./pack_results.sh                # -> results_<...>.tar.gz

  Copy ``results_*.tar.gz`` back and use **Import results**, then **Continue**.
* **Submit here** (when ``sbatch`` exists where the GUI runs): the same scripts, submitted at once; after
  the jobs, **Continue**.

The SLURM settings (partition, account, nodes, tasks per node, walltimes, ``module load`` lines, the VASP
and THERMACOND binaries on the cluster) are on the *Run* tab and kept in ``hpc.json``. The POTCAR is built
here and goes into the bundle. KPAR and NCORE in the INCARs are adjusted to the number of ranks the job
will run with.

For the supercells the last job runs ``04_supercells/collect.py`` (Python standard library), which does what
``read_outcar.x`` does: ``FORCEDISP1``, ``POSCAR1`` and ``stress.params`` with the same checks (converged
runs, net force, positions). ``pack_results.sh`` then brings back these files and the small logs, not the
OUTCARs (``./pack_results.sh --outcars`` adds them). On the 21 Cu₂O snapshots of 2026-11-07,
``collect.py`` and ``read_outcar.x`` gave identical numbers.

With DFPT the SCF and DFPT go in one bundle; the DFPT job waits for the SCF, skips itself when the SCF
shows no gap (``DFPT: auto``), and sets ``EFIELD_PEAD`` from the gap (``set_efield.py``).

The input files
---------------

The *Params* tab shows the five FOCEX files as FOCEX will get them, each line with a comment after its
values. Lines marked ``PLACEHOLDER`` are filled in by a later step or need your values: the lattice
before the relaxation, and the Born charges before DFPT. A file you change is kept in the instruction file
(as a block) and used as it is; **Back to the generated file** undoes that. A ``structure.params`` saved
before the relaxation gets the relaxed lattice and positions when the relaxation is done, keeping your
other lines. A saved ``dielectric.params`` is used instead of the DFPT result.

The *THERMACOND* tab holds the ``latdyn.params`` of ``07_thermacond``: line 8 the job (0 none,
1 linewidths, 2 and 3 the thermal conductivity), line 1 the q mesh, line 6 the temperatures. THERMACOND
runs here (with the "MPI ranks here" of the *Run* tab), as a bundle, or submitted.

The *Setup* tab shows the whole instruction file; *Results* the PDF plots of ``05_focex`` and
``07_thermacond``; *Snapshots* the SC-SNAPS displacement histograms; *Log* the job output and
``report.log``.

Projects
--------

A new project needs a name, a POSCAR and the main settings (functional, DFPT, force-constant ranks and
shells, supercell, snapshots, temperature, THERMACOND job). If the POSCAR is already relaxed *with the same
DFT settings*, tick the box (or press **Skip: POSCAR is relaxed**); otherwise relax it first, since the
supercell runs must use the settings of the relaxation.

**New project from the open one** keeps the relaxed cell, the SCF and DFPT, and makes the snapshots again
with another temperature, supercell or number of snapshots. A second supercell (``FORCEDISP2``) is a
separate project.

Only VASP is supported for now; the code-specific parts (job lines, result files) are in
``FOCEX/AGENT/bundle.py`` (``DFT_CODES``), so that Quantum ESPRESSO or SIESTA can be added there.
