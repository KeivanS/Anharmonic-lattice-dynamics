.. _focex-agent:

FOCEX agent: from a POSCAR to phonons
=====================================

``FOCEX/focex-agent.py`` runs the whole FOCEX workflow from the primitive cell of a crystal:
VASP relaxation and SCF, Born charges, SC-SNAPS snapshots, the supercell VASP runs, the
force collection, the five ``*.params`` files, FOCEX and the plots. Every step uses defaults
that can be overridden, and every step is checked and written to ``report.log``.

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

   cp FOCEX/focex-agent.example.txt instructions.txt      # edit: Project and RANKS are required
   ~/vaspenv/bin/python FOCEX/focex-agent.py init -i instructions.txt -s POSCAR --profile rivanna
   cd <Project> && sbatch run_workflow.sh                  # bash run_workflow.sh on a workstation

   ~/vaspenv/bin/python FOCEX/focex-agent.py run <Project> --only snaps     # redo one task
   ~/vaspenv/bin/python FOCEX/focex-agent.py run <Project> --retry-failed   # = resubmit_failed.sh
   ~/vaspenv/bin/python FOCEX/focex-agent.py plots <Project>                # remake the plots

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
