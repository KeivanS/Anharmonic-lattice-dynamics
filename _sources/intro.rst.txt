Introduction
============

What is ALATDYN?
----------------

Anharmonic-LATtice-DYNamics (ALATDYN) is a suite of stand-alone codes that
compute the thermodynamic and thermal transport properties of crystalline solids
from data on their forces and potential energy as a function of atomic
positions, obtained from density functional theory or classical molecular
dynamics:

* **FOrce Constant EXtraction (FOCEX)** extracts the harmonic and anharmonic
  force constants from force-displacement data, and computes from them the
  phonon spectrum, elastic constants and quasi-harmonic thermodynamic
  properties. Its output is the input of all the other codes.
* **Self-COnsistent PHonons (SCOP8)** calculates the state of equilibrium of a
  crystal at a given temperature and pressure from the force constants of the
  high-symmetry phase, by variational minimisation of the free energy.
* **THERMAl CONDuctivity (THERMACOND)** calculates phonon lifetimes, scattering
  rates and the lattice thermal conductivity from Boltzmann transport theory.
* **ANharmonic FOrce-field MOlecular Dynamics (ANFOMOD)** runs molecular
  dynamics with the extracted force field.

**SC-SNAPS**, distributed separately, sits in front of FOCEX: it builds the
supercell and the thermally displaced snapshots that you send to your DFT code,
and comes with a browser-based graphical interface. See :doc:`runscsnaps`.

A converter, **INterface FOrce-field CONverter (INFOCON)**, which translates
force constants between the formats of FOCEX, ShengBTE, ALAMODE and PHONOPY, is
under development.

Where to start
--------------

If you are new to the suite, read :doc:`theory` for the definitions and
conventions, then :doc:`install`, and then run the germanium example of
:doc:`runfocex` — it takes a few minutes and exercises the whole FOCEX pipeline.

Citation
--------

Please cite the following articles:

* K. Esfarjani and H. T. Stokes, *Phys. Rev. B* **77**, 144112 (2008).
  [`doi:10.1103/PhysRevB.77.144112 <https://doi.org/10.1103/PhysRevB.77.144112>`__]

* K. Esfarjani, H. T. Stokes and G. Chen, *Phys. Rev. B* **84**, 085204 (2011).
  [`doi:10.1103/PhysRevB.84.085204 <https://doi.org/10.1103/PhysRevB.84.085204>`__]

* K. Esfarjani, H. Stokes, S. Nayeb Sadeghi, Y. Liang, B. Timalsina, H. Meng,
  J. Shiomi, B. Liao and R. Sun, *ALATDYN: a set of Anharmonic LATtice DYNamics
  codes to compute thermodynamic and thermal transport properties of crystalline
  solids*, *Computer Physics Communications* **312**, 109575 (2025).
  [`doi:10.1016/j.cpc.2025.109575 <https://doi.org/10.1016/j.cpc.2025.109575>`__]

Getting help
------------

Questions, bug reports and suggestions are welcome on the
`ALADYN forum <https://matsci.org/c/aladyn/57>`_ hosted at the MatSci Community
Discourse, or through the
`GitHub issue tracker <https://github.com/KeivanS/Anharmonic-lattice-dynamics/issues>`_.

Acknowledgment
--------------

This project was supported by the following organizations:

* Supercomputer center of the IMR, Tohoku University
* Financial support by the University of California Energy Institute
* Cyberinfrastructure for Sustained Scientific Innovation, National Science Foundation
* Research Computing, University of Virginia

Contact
-------

| Prof. Keivan Esfarjani <k1 AT virginia DOT edu>
| Dept of Materials Science and Engineering
| University of Virginia
| Charlottesville, VA 22904
| 434-924-8029
| https://engineering.virginia.edu/faculty/keivan-esfarjani
