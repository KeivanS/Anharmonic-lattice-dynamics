Running SC-SNAPS
================

Overview
--------

**SC-SNAPS** builds a supercell from a primitive cell and generates a set of
snapshots in which *all* atoms are displaced at once, sampled from a classical
canonical (thermal) distribution at a chosen temperature. Those snapshots are
the configurations you send to your DFT code; the resulting forces are what
:doc:`runfocex` fits.

Displacing every atom in every snapshot, rather than one atom at a time, is what
makes the approach economical: each DFT run then yields :math:`3 N_{sc}`
independent pieces of information about the force constants instead of three.

SC-SNAPS is distributed from its own repository,
`github.com/KeivanS/SC_SNAPS <https://github.com/KeivanS/SC_SNAPS>`_, which
contains the Fortran generator, a browser-based graphical interface, a
POSCAR-to-XYZ converter and example input files. An older copy of the generator
alone is still bundled in ``FOCEX/utility``; prefer the standalone repository.

.. note::

   The two versions differ in their input format. In the current version, line 3
   of ``cell.inp`` carries **two** numbers (scale factor *and* the
   ``convcoord`` flag) and line 6 carries masses **and** charges. Input files
   written for the bundled copy will not work unchanged.

How it works
^^^^^^^^^^^^

#. The primitive cell is read and replicated to fill the supercell.
#. A model dynamical matrix is built from a simple harmonic pair interaction
   between neighbors within twice the shortest bond length, and diagonalised.
   The resulting frequencies are then rescaled so that their mean square equals
   the "average frequency" you specify — the model only has to give a
   *reasonable* set of eigenvectors and a *realistic spread* of frequencies.
#. For each snapshot, every mode :math:`\lambda` (except the three
   translational ones) is given a random classical thermal amplitude,

   .. math::
      u_i \;\propto\; \sum_{\lambda} \frac{e_{i\lambda}}{\omega_\lambda}
      \sqrt{\frac{k_B T}{m_i}} \; \xi_\lambda ,

   with :math:`\xi_\lambda` drawn from a normal distribution; a conjugate set of
   Gaussians gives the corresponding velocities.

Because the frequencies enter only through the ratio
:math:`\sqrt{T}/\omega`, the *shape* of the displacement distribution is
realistic even though the model potential is crude — you are choosing a
displacement amplitude, expressed in the physically meaningful units of a
temperature.

Installation
------------

.. code-block:: bash

   git clone https://github.com/KeivanS/SC_SNAPS.git
   cd SC_SNAPS
   make compile          # builds sc_snaps.x and installs it in ~/BIN

``make compile`` compiles ``sc_snaps.f90`` with ``gfortran -O2`` and moves the
executable to ``$(BINDIR)``, which defaults to ``~/BIN``. Override it with
``make compile BINDIR=/somewhere/else``.

For the graphical interface you additionally need:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Requirement
     - Notes
   * - Python ≥ 3.8
     -
   * - Flask
     - ``pip install flask`` — required to start the GUI at all.
   * - NumPy, Matplotlib, SciPy
     - Only for the analysis plots. Without them the GUI still runs, but the
       plot panels return "Missing library".
   * - `Jmol <https://jmol.sourceforge.net/>`_
     - Optional structure viewer. Any viewer that accepts a file argument
       (VESTA, OVITO, …) works — just point the GUI at it.

.. warning::

   If you compile ``sc_snaps.f90`` with ``-fcheck=all`` (as the ALATDYN
   ``make.inc`` does) the program aborts at run time with ``Recursive call to
   nonrecursive procedure 'findif'``: the finite-difference helper calls itself
   through ``d1v``/``d2v``, which ``-fcheck=recursion`` flags. SC-SNAPS's own
   Makefile uses plain ``-O2`` and is unaffected. If you need the checks, add
   ``-frecursive``.

Input files
-----------

All three files are read with list-directed Fortran ``read`` statements: only
the leading values of each line are parsed, and anything after them is a
comment. The number and order of lines is fixed.

``cell.inp`` — the primitive cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   1 1 1   90 90 90                 # 1  CONVENTIONAL cell: a b c alpha beta gamma
   0 0.5 0.5, 0.5 0 0.5, 0.5 0.5 0  # 2  primitive vectors in conventional units
   4.247  0                         # 3  scale factor (Ang); convcoord flag
   2                                # 4  number of atom types
   1 1                              # 5  number of atoms of each type in the primitive cell
   24.31  16.00   2 -2              # 6  the masses, then the charges
   Mg O                             # 7  element names
     0   0   0                      # 8  reduced coordinates, atom 1 (in units of conventional cell)
     0.5 0.5 0.5                    # 9  reduced coordinates, atom 2

**Line 1** — the conventional cell. ``a b c`` are multiplied by the scale factor
of line 3; the angles are in degrees.

**Line 2** — nine numbers read as three groups of three: the primitive
translation vectors expressed in units of the conventional cell vectors. Commas
are accepted as separators. Use ``1 0 0  0 1 0  0 0 1`` if the conventional cell
*is* the primitive cell.

**Line 3** — ``scale``, then ``convcoord``. **Both values are required.** The
flag decides the units of the atomic coordinates below *and* of
``supercell.inp``:

* ``convcoord = 0`` — reduced coordinates of the **conventional** cell;
* ``convcoord ≠ 0`` — reduced coordinates of the **primitive** cell.

**Line 4** — number of distinct elements.

**Line 5** — how many atoms of each type the primitive cell contains.

**Line 6** — ``ntype`` masses (uma) followed by ``ntype`` charges, all on one
line. The charges are carried through to ``snapshots.xyz`` for visualisation and
are not used in the sampling.

**Line 7** — the element names (two characters each).

**Lines 8 onward** — one line per atom, giving its three reduced coordinates in units of CONVENTIONAL vectors.
The atoms must be grouped by type, in the order declared on line 5. Trailing
text on these lines (an oxidation-state label, for instance) is ignored.

``supercell.inp`` — the supercell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   2 0 0
   0 3 0
   0 0 4

A 3×3 integer matrix whose **rows** are the supercell translation vectors, in
units of the conventional cell when ``convcoord = 0`` and of the primitive cell
otherwise. A diagonal matrix gives the usual :math:`n_1 \times n_2 \times n_3`
supercell; off-diagonal entries let you build a more compact cell of the same
volume, which is often a better use of a given atom budget. The number of atoms
generated is :math:`|\det| \times` (atoms in the primitive cell), and SC-SNAPS
stops if the volume ratio is not an integer.

``snaps.inp`` — the sampling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   400    # 1  average phonon frequency (1/cm)
   100    # 2  temperature (K)
   15     # 3  number of snapshots

**Average frequency** — the model frequencies are rescaled so that their mean
square matches this value. Set it near the middle of your material's phonon
spectrum (a Debye-like scale); 400 cm\ :sup:`-1` ≈ 12 THz is a reasonable
starting point for a hard covalent solid, less for a soft or heavy one. It fixes
the overall size of the displacements together with the temperature.

**Temperature** — the sampling temperature. Choose it at or somewhat above the
temperature at which you intend to use the force constants: large enough that
the anharmonic terms show in the forces, small enough that the truncated Taylor
expansion still applies. Check the amplitude of displacements in the histogram.
0.01-0.02 Ang is a reasonable value

**Number of snapshots** — 20–50 usually suffices; the practical ceiling is 999.
See the sizing rule in :doc:`runfocex`.

Output files
------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - File
     - Contents
   * - ``poscar_000``
     - The **undistorted** supercell, in VASP POSCAR format, Cartesian.
   * - ``poscar_001`` …
     - One thermally displaced snapshot each, same format.
   * - ``snapshots.xyz``
     - All snapshots in XYZ format, each atom line carrying element, position,
       charge and velocity — for visualisation and for the GUI's histograms.
   * - ``log.dat``
     - Cells and reciprocal cells, generated supercell coordinates, shortest
       bond length, cutoff, and the model frequencies.
   * - ``freqs.dat``, ``modes.dat``
     - The eigenvalues and eigenvectors of the model dynamical matrix.

The ``poscar_*`` files carry six columns per atom: the three Cartesian
coordinates followed by the three velocity components. VASP ignores the extra
columns.

.. important::

   ``POSCAR1`` for FOCEX is ``poscar_000`` **with the element-name line (line 6)
   deleted** — FOCEX expects the older VASP-4 layout in which the atom counts
   follow the third lattice vector directly:

   .. code-block:: bash

      sed '6d' poscar_000 > POSCAR1

   The ``poscar_0nn`` files themselves are fed to VASP unchanged.

The graphical interface
-----------------------

.. code-block:: bash

   make run     # starts the server and opens http://localhost:5050

The GUI is a local Flask application — nothing leaves your machine and no
internet connection is needed. The working directory is the one you launch it
from, and can be changed in the browser at any time without restarting.

The port (5050) is fixed. To restart after a crash, or to free the port:

.. code-block:: bash

   lsof -ti :5050 | xargs kill -9

What the interface offers
^^^^^^^^^^^^^^^^^^^^^^^^^

* **Paths.** The location of ``sc_snaps.x``, of the structure viewer and of the
  ``poscar2xyz.py`` converter, editable in the browser.
* **Input editors.** ``cell.inp``, ``snaps.inp`` and ``supercell.inp`` in three
  text areas, pre-filled with templates and annotated with format hints. *Load
  existing files* fills them from the working directory; **Save input files**
  writes them back. Editing alone does not touch the disk.
* **Run and watch.** *Run sc_snaps.x* saves the inputs, starts the job, and
  streams its output live into the log box.
* **Browse and visualise.** Generated ``poscar_*`` files and ``snapshots.xyz``
  appear as clickable chips; a click shows the file inline, and *Convert & open
  in Jmol* converts it to XYZ and launches the viewer.
* **Analysis plots.** Beyond what is needed simply to run the code, the GUI
  plots the quantities you want to see before committing to a set of expensive
  DFT runs:

  * a **thermal displacement histogram**, per element, with mean and standard
    deviation of :math:`|r - r_0|`;
  * the same resolved into **Cartesian components**;
  * a **velocity histogram**;
  * the **frequency spectrum** of the model dynamical matrix, from
    ``freqs.dat``.

.. tip::

   The displacement histogram is the fastest sanity check on your sampling. Mean
   displacements of a few hundredths of an Å are typical; if they approach a
   tenth of the bond length, the snapshots are too anharmonic for a
   Taylor-expanded force-constant model, and you should lower the temperature or
   raise the average frequency.

``defaults.json``
^^^^^^^^^^^^^^^^^

If a file named ``defaults.json`` exists in the working directory, it overrides
the interface's built-in defaults — the three executable paths and the three
input templates:

.. code-block:: json

   {
     "execpath":  "~/BIN/sc_snaps.x",
     "visualizer": "~/bin/jmol",
     "converter": "~/BIN/poscar2xyz.py",
     "cell":      "1 1 1   90 90 90\n 0 0.5 0.5, 0.5 0 0.5, 0.5 0.5 0\n 4.247 0\n ...",
     "snaps":     "400\n 300\n 10",
     "supercell": "2 0 0\n 0 2 0\n 0 0 2"
   }

An invalid or absent file is silently ignored and the built-in defaults are
used. Note that a template stored here must itself be a valid input file — in
particular, line 3 of the ``cell`` template needs both the scale factor and the
``convcoord`` flag.

``poscar2xyz.py``
^^^^^^^^^^^^^^^^^

The converter used by the *open in Jmol* buttons is also a standalone tool:

.. code-block:: bash

   python3 poscar2xyz.py poscar_001              # -> poscar_001.xyz
   python3 poscar2xyz.py poscar_*                # batch
   python3 poscar2xyz.py poscar_001 -o out.xyz   # explicit name

It handles both the modern and the old POSCAR layouts, direct and Cartesian
coordinates, negative (volume) scale factors and selective-dynamics lines.

From snapshots to force constants
---------------------------------

Once the DFT runs are done, collect the results into ``POSCAR1`` and
``FORCEDISP1`` and run FOCEX. That part of the workflow — the converters, the
file formats and the FOCEX input files — is described in :doc:`runfocex`.
