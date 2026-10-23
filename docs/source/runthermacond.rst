Running THERMACOND
==================

THERMACOND computes phonon properties and the lattice thermal conductivity from the
harmonic and cubic force constants fitted by FOCEX, by solving the linearized phonon
Boltzmann transport equation on the irreducible wedge of the Brillouin zone, with
tetrahedron integration of the energy-conserving delta functions. The theory and the
algorithms are described in S. Nayeb Sadeghi, S. Lee and K. Esfarjani,
*npj Comput. Mater.* **11**, 303 (2025), https://doi.org/10.1038/s41524-025-01673-8.

THERMACOND reads the same parameter files as FOCEX, so a FOCEX run folder is also a
THERMACOND run folder.

Building
--------

.. code-block:: bash

   cd THERMACOND
   make            # parallel build with mpif90  -> ~/ALADYN/bin/tcond-2.0.x
   make serial     # serial build, no MPI library -> ~/ALADYN/bin/tcond-2.0-serial.x
   make DEBUG=1    # either of the above with -g -fbounds-check instead of -O2

Objects and module files go to ``obj_mpi/`` or ``obj_serial/``. On macOS with
Homebrew's Open MPI, ``mpif90`` may pick an older ``gfortran`` than the one Open MPI was
built with (``Cannot read module file ... created by a different version of GNU
Fortran``), and the link may fail with ``library not found for -lSystem``. Both are fixed
by

.. code-block:: bash

   export OMPI_FC=/opt/homebrew/bin/gfortran
   export SDKROOT=$(xcrun --show-sdk-path)

Running
-------

Serial and parallel runs are the same program; the k-points of the irreducible zone are
shared among the processes:

.. code-block:: bash

   tcond-2.0.x                 # one process (or tcond-2.0-serial.x)
   mpirun -np 8 tcond-2.0.x    # eight processes

A serial run is the reference for a parallel one: both give the same output files. The
results do not depend on the compiler or the optimisation level either (checked with GCC 12
and 15, ``-O2`` and ``-O0``). Build ``make serial`` with the same compiler as ``mpif90``
(``make serial FC=...``) if the two are to be compared file by file.
Rank 0 writes all output in the run folder; the other processes work in
``.mpi_rank_NNN/``, which are removed at the end.

Input files
-----------

``structure.params``, ``dielectric.params``, ``kpbs.params``
    The FOCEX files (:ref:`focex-inputs`). From ``structure.params`` THERMACOND reads the
    lattice (lines 1–3, including ``convcoord`` on line 3), the ranks (line 4), the
    atom types, masses and names, and the atomic positions.

``lat_fc.dat``, ``fc2.dat``, ``fc3.dat``
    Written by FOCEX. THERMACOND reads ``lat_fc.dat`` by its section headers and takes
    the number of force constants from the lines present in ``fc2.dat`` and
    ``fc3.dat`` (FOCEX omits those below its ``margin``).

.. _thermacond-latdyn:

``latdyn.params``
    13 lines; FOCEX reads the first 7.

    .. code-block:: text

       14 14 14      # 1  k-mesh n1 n2 n3
       0 0 0  0 0 1  # 2  mesh shift (fractions of a step), surface normal
       300 350       # 3  number of DOS points, omega_max (1/cm)
       10 30         # 4  Gaussian DOS width, eta = Im(omega) in the self-energy (1/cm)
       .false.       # 5  verbosity
       0 900 31      # 6  Tmin, Tmax (K), number of temperatures (equally spaced)
       0             # 7  0 quantum (Bose-Einstein), 1 classical occupations
       2             # 8  job (below)
       0             # 9  direct solver: 0 conjugate gradient, 1 SVD
       0             # 10 V3*delta: 0 on disk in V3DELTA/ (reused), 1 in memory
       0             # 11 isotope scattering: 0 no, 1 natural abundances
       0 0 0 0       # 12 calc_cross, q_cross(3)
       100000        # 13 sample length (microns)

    Each line is checked: a missing line or value stops the program with the line
    number and what was expected.

    ========= ===============================================================================
    job       what is computed
    ========= ===============================================================================
    0         bands and Gruneisen parameters along ``kpbs.params``; frequencies, velocities
              and DOS (Gaussian and tetrahedron) on the k-mesh. Done for every job.
    1         0 + three-phonon linewidths of every mode along the band path, at every
              temperature
    2         0 + V3*delta in the irreducible zone; at every temperature the RTA and the
              direct solution of the BTE (line 9)
    3         0 + V3*delta; at every temperature the RTA and the iterative solution, for
              when the collision matrix does not fit in memory
    ========= ===============================================================================

    V3*delta (the three-phonon matrix elements times the tetrahedron weights) does not
    depend on temperature and is computed once. With line 10 = 0 it is written to
    ``V3DELTA/``, one file per irreducible k-point, and a later run on the same mesh
    reads it back instead of recomputing it: a new temperature range then costs
    seconds instead of the V3 time. Delete ``V3DELTA/`` after changing the force
    constants.

    Isotope scattering uses the natural isotopic abundances of the elements named in
    ``structure.params`` (built-in table).

How it is computed
------------------

* Group velocities: Hellmann-Feynman, v = <e|dD/dq|e>/2omega, including the derivative of the
  non-analytic (Born charge) term. The finite-difference velocities are written next to them in
  ``veltest.dat`` as a check; they agree to about 1e-4 except at band crossings.
* ``structure.params`` is optional: without it the cell, species, masses and positions are
  taken from ``lat_fc.dat`` (written by FOCEX since 23 October 2026, with the conventional
  cell). Results are identical either way (Ge 4x4x4 job 3).
* k-mesh (``latdyn.params`` line 1): integer fractions of the primitive reciprocal vectors of
  ``lat_fc.dat`` (the run stops if ``structure.params`` has a different primitive cell). Even
  meshes are Gamma-centred (reduced coordinates -1/2 ... 1/2 - 1/n), odd ones start at Gamma;
  the shift on line 2 is not used. Mesh points are found everywhere by their reduced
  coordinates modulo G (``get_k_info``); this is checked at start-up for k, k+G, -k and the
  irreducible points (``CHECK_MESH_LABELS`` in ``phonlog.dat``). Odd meshes stopped in the
  tetrahedron set-up before 23 October 2026.
* Irreducible zone: the stars use the Laue group of the crystal (its point group and
  inversion); the representative of each star is the image in the first Brillouin zone
  furthest along a fixed generic direction, so the irreducible points fill one wedge
  (x >= y >= z >= 0 for a cubic lattice, 0-30 degrees for a hexagonal one). The number of
  points and their weights agree with spglib.
* Degenerate modes (frequencies within 1e-3 1/cm): their eigenvectors are defined only up to a
  rotation within the degenerate set, so ``|V3|^2`` is averaged over the set at q1, q2 and q3,
  and so is ``|V3|^2`` times the tetrahedron weight (P1, P2) at q2 and q3. Then every partner
  of a set has the same lifetime and the same collision kernel, F = tau (v + Delta) with Delta
  common to the set, and the conductivity is per-mode value + kappa_P - per-mode RTA: the
  one-velocity-per-mode term is replaced by the trace over the set of the velocity matrix
  (the degenerate block of kappa_P). This is what ``kappa_vs_T.dat`` reports; checked by
  random rotations of the degenerate eigenvectors (``THERMACOND_ROTATE_DEGENERATE=<seed>``,
  test only): Ge 4x4x4 RTA and iterative results identical to all printed digits, while the
  per-mode values moved by 3 %.
* Wigner conductivity: in the RTA, ``kappa_wigner.dat`` gives the conductivity of Simoncelli,
  Marzari and Mauri (Nat. Phys. 15, 809 (2019)): kappa_P from pairs of degenerate modes (the
  usual RTA term, in a form that does not depend on the eigenvector basis) and kappa_C from
  the coherences between different modes. kappa_C is negligible for simple crystals and
  matters for complex or strongly anharmonic ones.

Output files
------------

============================== ======================================================================
file                           contents
============================== ======================================================================
``bs_freqs.dat``               band structure: index, path length, k, frequencies (1/cm)
``bs_grun.dat``                Gruneisen parameters along the path
``dos.dat``, ``dos_tet.dat``   DOS, Gaussian and tetrahedron
``KPOINT.FBZ``, ``KPOINT.IBZ`` k-mesh and irreducible wedge with weights
``veloc.dat``                  group velocities on the mesh (Hellmann-Feynman)
``veltest.dat``                finite-difference and Hellmann-Feynman velocities, band path and mesh
``kappa_vs_T.dat``             one line per temperature: RTA (= kappa_P) and direct (job 2) or
                               iterative (job 3) conductivity tensor xx yy zz xy xz yz (W/mK),
                               degenerate modes as blocks; job 2 adds the relative residual of the
                               linear solution, job 3 the iterations; then the per-mode values
                               (RTA, direct/iterative), which depend on the degenerate basis
``kappa_RTA.dat``              conductivity per band in the RTA, for each temperature
``kappa_wigner.dat``           jobs 2, 3: one line per temperature: RTA, kappa_P, kappa_C and
                               kappa_P + kappa_C (Wigner, RTA linewidths), xx yy zz xy xz yz
``distribution_F.dat``         jobs 2, 3: distribution function F (IBZ point, mode, x y z), RTA
                               and final, for each temperature
``tau_RTA.dat``                lifetimes (total, normal, Umklapp) in the irreducible wedge
``Qvalue.dat``                 Q and the three-phonon scattering rates in the irreducible wedge
``kappa_direct.dat``           job 2: conductivity per band from the direct solution
``tau_RTAandeff.dat``          job 2: RTA and effective lifetimes
``kappa_iter.dat``             job 3: conductivity per band at each iteration
``lw-selfenergy.dat``          job 1: linewidths along the path, one line per mode and k
``lw-selfenergy2.dat``         job 1: the same, one line per k
``V3DELTA/``                   V3*delta, one file per irreducible k-point (line 10 = 0)
``phonlog.dat``                log of the run
``times.dat``                  CPU time of each step
============================== ======================================================================

Checks
------

* :math:`\omega^2 < -25` cm\ :sup:`-2` at any k-point stops the run (instability or bad fit);
  smaller negative values (round-off near :math:`\Gamma`) are replaced by :math:`|\omega^2|`
  with a warning in ``phonlog.dat``. ``focex-agent.py`` flags the same threshold.
* The dynamical matrix is replaced by :math:`(D+D^\dagger)/2`; a warning in ``phonlog.dat``
  reports :math:`\max|D-D^\dagger|/\max|D| > 10^{-6}`, which points to an error elsewhere.
* An empty ``fc3.dat`` (all cubic FCs zero in the FOCEX fit) stops the run instead of giving
  infinite lifetimes.

Example
-------

With the FOCEX Ge example run in a folder (``FOCEX/example/Ge``), THERMACOND runs in the
same folder: set lines 1, 6 and 8 of ``latdyn.params``, e.g. ``6 6 6``, ``300 300 1``
and ``2``, then ``tcond-2.0.x``. On a 6x6x6 mesh (16 irreducible points) this takes about a minute
and a half on one core (almost all of it V3*delta) and gives at 300 K kappa_P = 30.90 W/mK in
``kappa_wigner.dat`` (kappa_C = 0.017); rerunning with
``100 400 4`` on line 6 reuses ``V3DELTA/`` and takes a few seconds. Converged values need
much denser meshes (the paper uses up to 30x30x30 and extrapolates in 1/n).

Status
------

As of October 2026:

* isotope scattering (line 11 = 1) is not yet available with jobs 2 and 3;
* the optical properties at the zone centre (self-energy at q = 0, dielectric function)
  of the previous job control are not yet part of job 0;
* the per-band files (``kappa_RTA.dat``, ``kappa_iter.dat``, ``kappa_direct.dat``,
  ``distribution_F.dat``) are per mode and depend on the degenerate basis; the totals in
  ``kappa_vs_T.dat`` do not;
* the iterative result for Ge 4x4x4 is anisotropic by 0.24 % (zz high), independently of the
  basis, while the RTA is exactly isotropic: the collision kernel does not have the full
  symmetry (tetrahedron orientation or the star handling of q2, q3), being checked together
  with the next item;
* the direct solution (job 2) and the iterative one (job 3) do not agree; the direct
  (IBZ-reduced) collision matrix is being checked;
* before 22 October 2026 the k-point symmetry was wrong for hexagonal, trigonal, monoclinic
  and triclinic cells.
