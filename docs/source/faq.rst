FAQ
===

Questions that came up while using FOCEX, with the answers we found. The main manual describes
the input and the output; the advice and the explanations of the subtler points are here.

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

The SVDs of the fit run on OpenMP threads (since 9.7); FOCEX chooses the number itself (about
80 % of the cores) and ``FOCEX_NUM_THREADS`` overrides it. The rest is serial. The expensive part
of a study is the DFT calculation of the snapshots: run them as independent jobs
(``FOCEX/utility/runall-slurm.sh``).

Installation
------------

**``make`` fails on the last line with "No such file or directory".**

The Makefile moves the executable to a directory that must already exist: the
``DESTDIR`` of ``make.inc`` (``~/ALADYN/bin``) for the utility programs. FOCEX and
THERMACOND create theirs. Create it with ``mkdir -p`` and re-run ``make``.

**Where are the programs after** ``make`` **?**

In ``~/ALADYN/bin``, named with their version (``focex-9.9.x``, ``tcond-2.1.x``). ``make
python-tools`` links the Python tools there: ``focex-agent``, ``allplots``, ``bandplot`` and
``shell-optimizer``. Job scripts that call an older ``focex-9.x.x`` must be updated by hand.

**After a** ``git pull`` **the build fails or the results look like the old version.**

Run ``make clean && make`` once when the compiler flags changed (9.7 introduced OpenMP and
``-O3``), and check with ``git describe --tags`` which version you have.

Preparing the input
-------------------

**How large must the supercell be?**

The range of the harmonic force constants can never exceed the largest sphere
inscribed in the Wigner-Seitz cell of the supercell. As a rule of thumb, aim for
a supercell of at least ~10 Å in each direction, or about 100 atoms, and check convergence by
repeating the fit with a larger one. FOCEX prints the cutoff it derived from
your supercell in the log file. A mode at a wavevector that the supercell does not contain (for
example the R point of a cubic perovskite in a 3×3×3 supercell) is only interpolated: to decide
whether it is stable, add a supercell that contains it (2×2×2) or compute it directly.

**How many snapshots do I need?**

Enough that the number of force components, :math:`3 N_{sc}` per snapshot,
comfortably exceeds the number of independent force constants (printed in the
log as ``total # of independent FCs of all rank``). 10–20 snapshots are usually
sufficient for a harmonic fit in a large supercell, 30–50 when cubic and quartic
terms are included. With single-atom displacements every inequivalent atom must be displaced at
least once, or its force constants are not determined at all.

**At which temperature should I sample the snapshots?**

At, or somewhat above, the temperature at which you intend to use the force
constants. The displacements must be large enough for the anharmonic terms to
show up in the forces, but small enough that the Taylor expansion truncated at
your highest rank is still accurate. Room temperature or lower is a reasonable default
for the sampling temperature in ``snaps.inp``. For harmonic force constants only, a typical
displacement amplitude of 0.01–0.02 Å is good. Displacements of 0.01 Å determine FC2 very well
but long-range FC3 and FC4 hardly at all.

**Should I displace one atom at a time instead?**

You can — FOCEX accepts any set of displaced configurations, including the
single-atom displacements of the frozen-phonon method. But thermal snapshots in
which all atoms move at once carry far more information per DFT run, which is
why SC-SNAPS generates them (:doc:`runscsnaps`).

**All my snapshots must come from the same DFT settings — how can I check that?**

Every snapshot must use the same functional, pseudopotentials, cutoff, k-mesh and convergence
criteria, and the same reference structure. A quick test: in a snapshot where one or two atoms
are displaced, atoms far from them (say more than 9 Å) must feel the same forces as in the
undisplaced reference. In a MoTe\ :sub:`2` data set three of seven snapshots failed this test
(far forces ~0 instead of the reference's ±4.9×10\ :sup:`-3` eV/Å, everywhere else agreement to
5×10\ :sup:`-5`); without them the FC2 fit error went from 12.7 % to 3.2 %. Remove such snapshots
(and make sure every inequivalent atom is still displaced somewhere).

**Do I need the Born effective charges?**

Only for polar materials, where the long-range dipole-dipole interaction makes
the force constants decay slowly and produces the LO-TO splitting. For non-polar
crystals (Si, Ge) or for metals set ``born_flag = 0`` and leave the charge blocks at
zero. For polar crystals compute :math:`\epsilon^\infty` and :math:`Z^*` with
your DFT code (DFPT), put them in ``dielectric.params`` and use
``born_flag = 2``. The DFPT run must be consistent with the supercell runs: same functional,
pseudopotentials and cutoff (ENCUT), and a k-mesh commensurate with the supercell ones. The fit
only sees the combination :math:`Z^*/\sqrt{\epsilon}`, so inconsistent values leave a long-range
error in the forces that the short-range force constants cannot absorb.

Choosing the ranges and the ranks
---------------------------------

**How do I choose the number of neighbour shells for each atom?**

Look at the shells of each atom (``pairs.dat``, or the shell table in the log: radius, number of
neighbours, gap to the previous shell) and put the cutoff where the gap between consecutive shells
is large. Then check convergence: increase the range until the fit error, the phonons and the
elastic constants stop changing. ``shell-optimizer DIR`` does this for you: it reads the shells,
proposes cutoffs at the large gaps, runs FOCEX for increasing FC2 and then FC3 ranges, and reports
the smallest FC2 range that is converged, stable and has a small fit error, and the largest FC3 range
whose added force constants are still determined by the data, together with the lines to put in
``structure.params``. If no range is stable with your invariances, it tells you which one to switch off. ``shell-optimizer --dry DIR`` only lists the shells and the candidates.

**Should I use the default FC2 range** (``fc2range = 0``) **?**

It takes every pair inside the supercell. That is often too many: in layered crystals it includes
many weak interlayer pairs that the data barely determine, and the neighbour search may not reach
them (the run stops with ``atom not found ... increase the cutoff length in 3rd line of
default.params``). Pairs that sit exactly on the boundary of the supercell can only be determined
in part; FOCEX now limits the damage (see "cut moved" below), but a range chosen from the shell
gaps is safer.

**Which ranks should I include?**

FC2 for phonons and elastic constants; FC3 for Grüneisen parameters, thermal expansion and thermal
conductivity; FC4 only when the displacements are large enough to see it and the data determine it
(check its error bars: in BaZrO\ :sub:`3` with 0.01–0.02 Å displacements only 38 of 399 FC4 were
larger than twice their error, and adding them made the fit and the bulk modulus worse). The range
of FC3 is usually much shorter than that of FC2.

Invariances and the fitting method
----------------------------------

.. _faq-enforce-inv:

**Which** ``enforce_inv`` **should I use?**

``1``. The invariances then hold exactly and every force constant gets an error bar. With ``0``
they hold only approximately, and how well depends on the data (check the violation tables in the
log). ``3`` (LASSO with a refit) is for selecting force constants when you include many of them and
the forces are accurate; ``2`` only as a diagnostic. In MgO with two snapshots, ``1`` fitted the
forces 25 times better than ``3`` (0.17 % against 4.3 %) with the same phonons: LASSO removed force
constants the data support.

**What do** ``enforce_inv = 2`` **and** ``3`` **actually do?**

For each value :math:`\lambda` of the weight of the :math:`L_1` term, FOCEX computes the LASSO
solution (many force constants exactly zero) and, on the force constants it kept, a "refit": an
ordinary least-squares fit with all the others held at zero and the invariances still exact. It
then chooses :math:`\lambda` (by cross-validation, or by the noise test of ``default.params`` line
6) and returns the LASSO solution (``2``) or the refit (``3``) at that :math:`\lambda`. The LASSO
values are shrunk towards zero; the refit removes that shrinkage. Flowcharts of both are in the
technical notes ``FOCEX/notes/focex_9x_algorithms.pdf``.

**Should the fit error go down when I add an invariance?**

No — on the data being fitted it can only stay the same or go up, because the invariances restrict
the fit. What they should do is make the model more physical and its predictions on data not used
in the fit better. So compare the fit error with and without the invariance: a small increase is
normal; a large increase means the model is not consistent with that invariance, and the ranks,
the ranges, the stress or the Born charges need to be reconsidered. Examples from MoTe\ :sub:`2`:

* Huang on, with an unrelaxed reference and no ``stress.params`` (Huang then assumes zero stress):
  fit 3.2 % and the interlayer shear mode at −65 cm\ :sup:`-1`; Huang off: 1.6 % and 26 cm\ :sup:`-1`.
* rotational invariance on, with FC3 to 2 shells and FC2 to 10 shells (the rotational relations tie
  FC2 to FC3 and cannot hold when FC3 is cut much shorter): fit 2.8 %, the breathing mode at
  −29 cm\ :sup:`-1` and a negative C\ :sub:`33`; off: 0.84 %, 37 cm\ :sup:`-1`, all stable. With FC3 to
  3 shells it was stable again with it on.

**When should I switch Huang** (``ihuang``) **or rotational** (``irot``) **invariance off?**

Huang: when the reference is not relaxed and you do not have its stress from the DFT run (give
it in ``stress.params`` if you have it). Rotational: when the FC3 range is much shorter than the FC2
range and switching it off removes instabilities or lowers the fit error a lot.
``shell-optimizer`` reruns its chosen model with each of them off and tells you if one is a problem.

**My acoustic branches are imaginary near** :math:`\Gamma` **.**

The acoustic sum rule is not satisfied. Set ``itrans = 1`` and use ``enforce_inv = 1``, which holds
the invariances exactly. If only a few points very close to :math:`\Gamma` are slightly negative
with Born charges on, compare with ``born_flag = 0``: it is then a long-range (q→0) effect rather
than an instability.

Judging the result
------------------

**How do I know whether the fit is good?**

* the fit error (``Relative residual`` or ``Percent error`` in the log): a few percent at most;
* the error bars: each force constant with its standard error is in the log, the mean per rank
  (``sd(rank)``) at its end; force constants smaller than twice their error are not determined;
* the invariance violations in the log (with ``enforce_inv = 1`` at round-off level);
* the phonons: no imaginary frequencies (other than where expected), acoustic branches linear at
  :math:`\Gamma`;
* the elastic constants: positive, symmetric, and for a hexagonal crystal
  :math:`C_{66} = (C_{11}-C_{12})/2`;
* the Grüneisen parameters: of order 1; very large values belong to modes with very low frequency
  (:math:`\gamma \propto 1/\omega^2`) — for soft modes look at :math:`d\omega^2/d\ln V` instead;
* predictions on snapshots that were not used in the fit, if you can spare a few: the force error
  there should be close to the fit error.

**Different ranges give very different results for one mode.**

Then the data do not determine that mode. Example: BaZrO\ :sub:`3` in a 3×3×3 supercell gave the
R-point tilt mode at +49 cm\ :sup:`-1` with one setting and −34 cm\ :sup:`-1` with another; R is not in
that supercell. Add data that contain the mode.

**The fit residual got worse when I turned on LASSO** (``enforce_inv = 2, 3``).

Expected: the residual is measured on the snapshots the model was fitted to, and LASSO uses fewer
force constants. Judge the result on the phonons and on the error bars. If it is much worse (as in
the MgO example above), use ``enforce_inv = 1``.

**The optical frequencies are wrong for my polar material.**

The long-range electrostatics is missing or incomplete: provide
:math:`\epsilon^\infty` and the Born charges in ``dielectric.params`` and set
``born_flag = 2``, consistent with the supercell runs (see above). The LDA or GGA underestimation
of the bandgap can also lead to a larger dielectric constant and a lower LO-TO splitting. In some
materials with transition metals, adding a Hubbard U term may improve the accuracy, otherwise use
of hybrid functionals such as HSE will also improve accuracy of :math:`\epsilon^\infty`.

Messages during the run
-----------------------

**"SVD: singular values next to the cut ...; cut moved to the gap at ..."**

Some combinations of force constants are (almost) not determined by the data, and their singular
values sat next to the SVD cutoff, where they would have been amplified. FOCEX moved the cutoff into
the gap above them. Nothing to do; it is the cure for force constants of 10\ :sup:`2`–10\ :sup:`5`
eV/Å\ :sup:`2` that appeared before 9.9 (MgO with two displacements and the default range).

**"SVD: measured relative noise ... refit with the cut ..." and then "refit kept" or "refit rejected"**

The fit error is larger than the force noise of ``default.params`` line 2 assumes. FOCEX tried a
larger cutoff and kept it only if the fit stayed within the noise. "Rejected ... model error" means
the misfit comes from the model (too short a range, missing ranks, an inconsistent invariance), not
from noise: improve the model.

**"WARNING: the invariance rows with a right-hand side cannot all hold"**

The invariance relations that carry a right-hand side (from ``stress.params``, from ranks read with
``include_fc = 2``, or from the Coulomb terms with ``born_flag = 2``) contradict each other or the
other relations. Check that the read force constants, the stress and the Born charges belong to this
structure and these settings.

**"atom not found ... increase the cutoff length in 3rd line of default.params"**

The requested range (often the default FC2 range of a large supercell) reaches beyond the
neighbour list. Choose a shorter range, or increase line 3 of ``default.params``.

**"CHECK_BORN_SYMMETRY: Born charges or dielectric tensor break the crystal symmetry"** and
**POSCAR atoms that are "nearest to the site ... at ... A from it"**

The Born charges are on the wrong atoms or in another frame, or POSCAR and ``structure.params`` do
not describe the same lattice, origin or atom order. The messages give the details; see
:doc:`runfocex` ("Error messages").

**FOCEX stops immediately with an end-of-file or "bad real number" message.**

One of the parameter files has the wrong number of lines or values. These files
are read positionally: ``default.params`` has 4 lines (optional 5th to 7th) and line 2 of
``latdyn.params`` must contain six numbers (mesh shift *and* surface normal).
Start from the templates in the ``FOCEX`` directory.

Output files
------------

**Which output files do the other ALATDYN codes need?**

``lat_fc.dat`` together with ``fc1.dat``, ``fc2.dat``, ``fc3.dat`` and, if you
fitted it, ``fc4.dat``.

**Why is the log file named** ``log1dfB00_00_03_0_tr00.dat`` **?**

The name encodes the settings of the run — the number of data files, the FC2
range mode, the Born flag, the shells used for ranks 2/3/4 and the invariances
imposed — so that runs with different options do not overwrite each other. The
decoding is given in :ref:`focex-logname`.
