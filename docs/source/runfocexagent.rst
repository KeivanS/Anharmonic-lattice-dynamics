.. _focex-agent:

FOCEX agent: from a POSCAR to phonons
=====================================

``FOCEX/AGENT/focex-agent.py`` runs the FOCEX workflow in one of two ways:

* **full** (``-s POSCAR``): from the primitive cell of a crystal -- VASP relaxation and SCF,
  Born charges, SC-SNAPS snapshots, the supercell VASP runs, the force collection, the five
  ``*.params`` files, FOCEX and the plots;
* **FOCEX only** (``DATA: <folder>`` in the instructions): from existing ``POSCARi`` /
  ``FORCEDISPi`` (and optionally Born charges) -- the ``*.params`` files, FOCEX and the plots;
  see *Existing data* below.

Every step uses defaults that can be overridden, and every step is checked and written to
``report.log``.

Requirements
------------

* `VASP-Flow <https://github.com/KeivanS/VASP-Flow>`_ in ``~/VASP-Flow`` (or ``$VASP_FLOW_DIR``):
  the agent imports its INCAR/KPOINTS generator, GGA+U and MAGMOM rules, KPAR/NCORE rule,
  POTCAR builder, SLURM profiles and high-symmetry k-path. The POTCAR library is
  ``$VASP_POTCAR_DIR``, as for VASP-Flow.
* Python with numpy, spglib and pymatgen (VASP-Flow's environment, e.g. ``~/vaspenv/bin/python``;
  a bare Homebrew ``python3`` has none of them).
* ``sc_snaps.x`` (SC-SNAPS), ``read_outcar.x`` (``FOCEX/utility``), the FOCEX binary, and gnuplot
  for the plots.

Running it
----------

.. code-block:: bash

   cp FOCEX/AGENT/focex-agent.example.txt instructions.txt      # edit: Project and RANKS are required
   ~/vaspenv/bin/python FOCEX/AGENT/focex-agent.py init -i instructions.txt -s POSCAR --profile rivanna
   cd <Project> && sbatch run_workflow.sh                  # bash run_workflow.sh on a workstation

   # FOCEX only, on existing POSCAR1/FORCEDISP1 (no VASP, no SC-SNAPS): DATA: <folder> in instructions.txt
   ~/vaspenv/bin/python FOCEX/AGENT/focex-agent.py init -i instructions.txt
   ~/vaspenv/bin/python FOCEX/AGENT/focex-agent.py run <Project>

   ~/vaspenv/bin/python FOCEX/AGENT/focex-agent.py run <Project> --only snaps     # redo one task
   ~/vaspenv/bin/python FOCEX/AGENT/focex-agent.py run <Project> --retry-failed   # = resubmit_failed.sh
   ~/vaspenv/bin/python FOCEX/AGENT/focex-agent.py plots <Project>                # remake the plots

``run_workflow.sh`` calls ``focex-agent.py run <Project>``, which does the tasks below in order.
A finished task is recorded in ``state.json`` and skipped, so after a time limit the same script
is simply submitted again. ``focex-agent.py run <Project> --only <task>`` reruns one task.

========== ========================================================================================
task       what it does
========== ========================================================================================
relax      Niggli-reduced, oriented cell (``POSCAR.init``); VASP relaxation in ``01_relax``,
           restarted from CONTCAR up to three times if EDIFFG is not reached
poscar0    CONTCAR is Niggli-reduced and oriented again: ``POSCAR0``, the cell and Cartesian
           frame of every later step
scf        ``02_scf``: same INCAR (ISIF = 1, NSW = 0) and KPOINTS; reuses CHGCAR/WAVECAR when
           the second reduction left the basis unchanged; band gap from EIGENVAL
dfpt       ``02_scf/dfpt``, only if there is a gap (``DFPT: auto``): epsilon_inf and Born charges
snaps      ``03_snaps``: ``cell.inp``, ``supercell.inp``, ``snaps.inp``; SC-SNAPS; mean and standard
           deviation of the displacements; if the rms displacement is not within 15 % of
           ``DISP_TARGET`` the temperature is rescaled (T proportional to the mean square) and
           SC-SNAPS rerun, at most twice; the
           undistorted supercell is checked to be commensurate with POSCAR0 before any VASP run
supercells ``04_supercells/snap_NNN``: INCAR as the SCF plus ISYM = 0 and no WAVECAR/CHGCAR,
           Gamma mesh (see *Supercell k-mesh* below), MAGMOM expanded to the supercell
vasp_snaps the supercell runs, one after the other; failures are diagnosed in ``report.log``
           and ``resubmit_failed.sh`` reruns them (no OUTCAR is ever overwritten: an old one
           is kept as ``OUTCAR.1``, ``OUTCAR.2`` ...)
collect    ``read_outcar.x`` in every folder, in order: ``FORCEDISP1``; ``POSCAR1`` = ``poscar_000``;
           checks the net force and that the OUTCAR positions are those of the POSCAR
focex      ``05_focex``: the five ``*.params`` files and FOCEX. With Born charges, also a
           ``NOBORN`` run (born_flag = 0). If equations/unknowns < ``LASSO_RATIO`` (2), also a
           ``LASSO`` run (enforce_inv = 2) and a comparison. The DOS range of ``latdyn.params``
           is reset to 1.1 x the highest frequency and FOCEX rerun when it does not fit
plots      ``bvel.plt``, ``all.plt``, ``qha.plt`` with their placeholders filled (k ticks, ranges,
           one plot entry per band) and ``makeplots.sh``, in every FOCEX folder
========== ========================================================================================

For each FOCEX run the report gives the fit error, the number of equations and unknowns, the
lowest omega^2 (a warning below -25 cm^-2) and the highest frequency (a warning above
2500 cm^-1).

Defaults
--------

* ``ACCURACY: Accurate`` -- ENCUT 600 eV, EDIFF 1E-7, EDIFFG -0.001, PREC = Accurate in every
  step; ``Normal`` -- 450 eV, 1E-6, -0.003, PREC = Normal. DFPT uses EDIFF = 1E-8.
* Relaxation and SCF k-mesh: VASP-Flow ``coarse`` (1000 k-points per reciprocal atom), ``fine``
  (6000) if every element is a metal. The report warns if the SCF of a non-metal finds no gap.
* ``SUPERCELL: SMALL`` -- the m1 x m2 x m3 multiple of the conventional cell with the fewest atoms
  whose smallest face distance is at least 10 A, within 100 atoms (``LARGE``: 14 A, 300 atoms).
  Multiples of the conventional cell are used rather than the most compact Hermite normal form:
  a skewed supercell of the same size can put FC2 pairs outside its Wigner-Seitz cell, and FOCEX
  then stops in ``set_dynmat_std`` with ``R+tp-t not inside!``.
* 50 snapshots, 30 K, 400 cm^-1, target rms displacement 0.02 A (all atoms together: light atoms
  move more than heavy ones).
* ``latdyn.params`` mesh ``fine``; ``kpbs.params`` from the automatic path, in primitive units.
* ``RANKS`` has no default. ``RANGE: 6 2 1 0 0 0 0`` shells for ranks 2 to 8, ``FC2RANGE: 0``.

``structure.params`` is written as: line 1 a b c alpha beta gamma of the conventional cell,
line 2 the primitive vectors in conventional units, line 3 ``1.0 0`` (atoms in conventional
units), line 5 ``1 1 1 0``, line 7 ``1 .false.``. POSCAR0 is rotated so that its conventional
cell lies exactly where FOCEX and SC-SNAPS build it from these six numbers (a along x, b in the
xy plane), so POSCAR1, FORCEDISP1 and the Born charges are all in that frame.

Existing data (FOCEX only)
--------------------------

To fit force constants to forces that are already computed, give ``DATA: <folder>`` in the
instructions (no ``-s``): ``focex-agent.py init -i instructions.txt`` and then ``run``. The folder
holds ``POSCAR1``, ``FORCEDISP1`` (and ``POSCAR2``, ``FORCEDISP2``, ... -- all pairs are used)
in any orientation, and optionally

* ``POSCAR``: the primitive cell in the frame of the data; otherwise spglib finds it from
  ``POSCAR1`` (keeping its orientation);
* ``dielectric.params``: epsilon_inf and the Born charges in the frame of the data, one block per
  atom of that ``POSCAR`` (without a ``POSCAR`` this is accepted only with one atom per species).

``SPECIES: Mg O`` names the species when the POSCARs have no element line. ``init`` symmetrizes
the lattice of the primitive cell (spglib's idealized lattice parameters, same basis and atom
order: lattice parameters given to 6-7 digits make FOCEX's symmetry analysis lose force-constant
groups), rotates it into the FOCEX frame, rotates ``POSCARi`` and the positions and forces of
``FORCEDISPi`` with the same matrix (best fit) and epsilon and Z as R^T T R, writes ``POSCAR0`` and
``04_supercells/``, and marks relax ... collect as provided; ``run`` then does ``focex`` and
``plots`` as for a full run, with the same blocks and keywords (``RANKS``, ``RANGE``,
``FC2RANGE``, ``structure_params``, ...). The report gives the rotation and the size of the
symmetrization.

Checked on MgO data in a rhombohedral orientation (``MgO-Safoura``): without Born charges the
fit error (7.197 %), the 35 unknowns, the highest frequency (647 1/cm), the bulk and shear
moduli are those of a FOCEX run in the original frame; the Zener ratio and the mean sound speed
(Debye temperature 982 K instead of 806 K) are meaningful only in the FOCEX frame, where the
conventional axes are x, y, z.

Supercell k-mesh
----------------

If the supercell matrix (rows, primitive units) is diagonal up to a permutation, the mesh is
the primitive one divided by the multiples, rounded up: 6 10 8 with 2 x 2 x 2 gives 3 5 4,
exactly commensurate. Otherwise -- which includes every multiple of a centred conventional cell,
e.g. fcc or bcc -- the supercell mesh only keeps the k-point spacing of the primitive mesh
(n_i = round(b_i(sc) / dk), b the lengths of the reciprocal vectors, dk the geometric mean of b_i / N_i of the primitive mesh); it is
then not in general exactly commensurate. ``KPOINTS_snaps`` sets it by hand.

Project folder
--------------

=========================== ==============================================================
``instructions.txt``        copy of the input (``init`` appends the k-mesh density it chose)
``POSCAR.input``            the structure given; ``POSCAR.init`` reduced and oriented
``POSCAR0``                 relaxed cell, frame of every later step
``state.json``              finished tasks, supercell matrix, gap, epsilon, Born charges
``report.log``              one entry per task: OK/FAILED, minutes, what was checked
``run_workflow.sh``         the batch script (VASP-Flow profile); ``resubmit_failed.sh``
``01_relax``, ``02_scf``    VASP; ``02_scf/dfpt`` the DFPT run
``03_snaps``                ``cell.inp``, ``supercell.inp``, ``snaps.inp``, ``poscar_NNN``
``04_supercells``           ``snap_NNN/`` one VASP run each; ``FORCEDISP1``, ``POSCAR1``
``05_focex``                the FOCEX run and plots; ``NOBORN/``, ``LASSO/`` when made
=========================== ==============================================================

WAVECAR, CHG, CHGCAR, WFULL, TMPCAR and PROCAR are deleted from ``01_relax`` after the SCF, from
``02_scf`` after the DFPT (or when it is skipped), and from each ``snap_NNN`` after its run.
OUTCARs are never deleted.

Thermal conductivity
--------------------

THERMACOND reads ``lat_fc.dat``, ``fc2.dat`` and ``fc3.dat`` of the FOCEX run, so give at least
``RANKS: 1 1 1 0 0 0 0 0``. The agent does not run it: in ``05_focex`` set lines 8-13 of
``latdyn.params`` (job, solver, ...; see :ref:`thermacond-latdyn`) and run ``kap8`` there,
or in a copy with ``structure.params``, ``dielectric.params``, ``kpbs.params`` and the three
``*.dat`` files.

Overriding any file
-------------------

A block named after a file replaces it, e.g.::

   structure_params
   ... the whole file ...
   structure_params_END

``INCAR`` (every step) and ``INCAR_relax``, ``INCAR_scf``, ``INCAR_dfpt``, ``INCAR_snaps`` are merged
tag by tag into the generated INCARs. The other blocks are written verbatim: ``KPOINTS``
(relaxation and SCF), ``KPOINTS_snaps``, ``cell_inp``, ``supercell_inp``, ``snaps_inp``,
``structure_params``, ``latdyn_params``, ``kpbs_params`` (required with ``KPATH: 1``),
``dielectric_params``, ``default_params``. ``END_<NAME>`` is also accepted, in any case.

Born charges
------------

VASP prints Z*\ :sub:`ij` with i the field and j the displacement direction. FOCEX contracts the
wave vector with the second index of ``charge(i,j)`` (``aux_na.f90``), so the agent writes the
transpose. ``BORN_TRANSPOSE: no`` writes them as printed.
