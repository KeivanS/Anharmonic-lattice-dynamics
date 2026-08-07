FAQ
===

General
-------

**Which code do I need for what?**

Start with FOCEX: every other code in the suite consumes its force constants.
Then use THERMACOND for phonon lifetimes and lattice thermal conductivity,
SCOP8 for the temperature-dependent equilibrium structure and free energy of an
anharmonic crystal, and ANFOMOD to run molecular dynamics with the extracted
force field.

**Can I use a DFT code other than VASP or Quantum Espresso?**

Yes. FOCEX only reads the two plain-text files ``POSCARi`` and ``FORCEDISPi``;
the converters for VASP and QE are a convenience, not a requirement. Any code
that gives you forces on all atoms of a supercell will do — write the results in
the format described in :ref:`focex-forcedisp` (positions in Å, forces in eV/Å).

**Does FOCEX run in parallel?**

No, FOCEX is serial. The expensive part of a study is the DFT calculation of the
snapshots, which is embarrassingly parallel: run the snapshots as independent
jobs (see ``FOCEX/utility/runall-slurm.sh``).

Preparing the input
-------------------

**How large must the supercell be?**

The range of the harmonic force constants can never exceed the largest sphere
inscribed in the Wigner-Seitz cell of the supercell. As a rule of thumb, aim for
a supercell of at least ~10 Å in each direction, or about 100 atoms, and check convergence by
repeating the fit with a larger one. FOCEX prints the cutoff it derived from
your supercell in the log file.

**How many snapshots do I need?**

Enough that the number of force components, :math:`3 N_{sc}` per snapshot,
comfortably exceeds the number of independent force constants (printed in the
log as ``total # of independent FCs of all rank``). 10–20 snapshots are usually
sufficient for a harmonic fit in a large supercell, 30–50 when cubic and quartic
terms are included. 

**At which temperature should I sample the snapshots?**

At, or somewhat above, the temperature at which you intend to use the force
constants. The displacements must be large enough for the anharmonic terms to
show up in the forces, but small enough that the Taylor expansion truncated at
your highest rank is still accurate. Room temperature or lower is a reasonable default
for the sampling temperature in ``snaps.inp``. For Harmonic force constants only, a typical 
displacement amplitude of 0.01-0.02 Ang should be good. 

**Should I displace one atom at a time instead?**

You can — FOCEX accepts any set of displaced configurations, including the
single-atom displacements of the frozen-phonon method. But thermal snapshots in
which all atoms move at once carry far more information per DFT run, which is
why SC-SNAPS generates them (:doc:`runscsnaps`).

**Which SC-SNAPS should I use — the one in FOCEX/utility or the separate repo?**

The separate repository, `github.com/KeivanS/SC_SNAPS
<https://github.com/KeivanS/SC_SNAPS>`_. It is the maintained version and adds a
graphical interface. Its ``cell.inp`` format has changed, so old input files
need editing; see :doc:`runscsnaps`.

**Do I need the Born effective charges?**

Only for polar materials, where the long-range dipole-dipole interaction makes
the force constants decay slowly and produces the LO-TO splitting. For non-polar
crystals (Si, Ge) or for metals set ``born_flag = 0`` and leave the charge blocks at
zero. For polar crystals compute :math:`\epsilon^\infty` and :math:`Z^*` with
your DFT code (DFPT), put them in ``dielectric.params`` and use
``born_flag = 2``.

Running and interpreting
------------------------

**How do I know whether the fit is good?**

Check, in this order: the relative force error at the end of ``svd-all.dat``,
the condition number of the fit matrix listed just below the singular values,
the invariance violations reported in the log, and finally the physics — the
acoustic branches must go to zero linearly at :math:`\Gamma`, and no branch may
be imaginary. See "Checking the quality of the fit" in :doc:`runfocex`.

**My acoustic branches are imaginary near** :math:`\Gamma` **.**

The acoustic sum rule is not satisfied. Set ``itrans = 1`` and if needed ``irot = 1; ihuang = 1`` in
``structure.params`` and, if the violations persist, ``enforce_inv = 1`` so that
the invariances are imposed exactly by projection rather than in the
least-squares sense.

**The optical frequencies are wrong for my polar material.**

The long-range electrostatics is missing or incomplete: provide
:math:`\epsilon^\infty` and the Born charges in ``dielectric.params`` and set
``born_flag = 2``. The LDA or GGA underestimation of the bandgap can also lead to 
a larger dielectric constant and a lower LO-TO splitting. In some materials with transition
metals, adding a Hubbard U term may improve the accuracy, otherwise use of hybrid functionals
such as HSE will also improve accuracy of :math:`\epsilon^\infty` .

**FOCEX stops immediately with an end-of-file or "bad real number" message.**

One of the parameter files has the wrong number of lines or values. These files
are read positionally: ``default.params`` must have 36 lines and line 2 of
``latdyn.params`` must contain six numbers (mesh shift *and* surface normal).
Copy the templates from the ``FOCEX`` directory rather than reusing files from
an older version of the code.

**Which output files do the other ALATDYN codes need?**

``lat_fc.dat`` together with ``fc1.dat``, ``fc2.dat``, ``fc3.dat`` and, if you
fitted it, ``fc4.dat``.

**Why is the log file named** ``log1dfB00_00_03_0_tr00.dat`` **?**

The name encodes the settings of the run — the number of data files, the FC2
range mode, the Born flag, the shells used for ranks 2/3/4 and the invariances
imposed — so that runs with different options do not overwrite each other. The
decoding is given in :ref:`focex-logname`.

Compilation
-----------

**``sc_snaps.x`` compiles but stops with "Recursive call to nonrecursive
procedure 'findif'".**

It was built with ``-fcheck=all``, which turns on ``-fcheck=recursion``, and the
finite-difference helper genuinely recurses. Build it with plain ``-O2``, or
keep the checks and add ``-frecursive``; see :ref:`focex-install`.

**``make`` fails on the last line with "No such file or directory".**

The Makefile moves the executable to a directory that must already exist:
``~/BIN`` for FOCEX, the ``DESTDIR`` of ``make.inc`` for the utility programs.
Create it with ``mkdir -p`` and re-run ``make``.
