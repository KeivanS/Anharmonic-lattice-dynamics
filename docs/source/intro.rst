Introduction
============

What is ALADYN?
---------------

Anharmonic-LAttice-DYNamics (ALADYN) contains 5 sets of stand-alone codes with the following functions:

* **FOrce Constant EXtraction (FOCEX)** to extract force constants from force-displacements data, the output of which can be used as input to the following codes;
* **Self-CoOnsistent Phonon (SCOP)** to calculate the state of equilibrium of a crystal at a given temperature and pressure solely from the set of force constants parameters of the high-symmetry phase. It is based on a variational free energy minimization.
* **THERMAl CONDuctivity (THERMACOND)** to calculate the phonon spectrum, group velocities and lifetimes, and lattice thermal conductivity based on Boltzmann equation theory;
* **ANharmonic FOrce-field MOlecular Dynamics (ANFOMOD)** to perform molecular dymics from the above-developed force field;
* **INterface FOrce-field CONverter (INFOCON)** to convert extracted force constants by one of the codes FOCEX, ShengBTE, ALAMODE, and PHONOPY the format needed as input for the phonon calculation by another of these codes so that these codes become intercompatible.

Citation
--------

Please cite the following articles:

* K. Esfarjani and H.T. Stokes, *Phys Rev B* **77**, 144112 (2008).
  [`Link <https://doi.org/10.1103/PhysRevB.77.144112>`__]

* K. Esfarjani, H. T. Stokes anf G. Cgen, *Phys Rev B* 84, 085204 (2011). DOI: 10.1103/PhysRevB.84.085204
  [`Link <https://doi.org/10.1103/PhysRevB.84.085204>`__]

* K. Esfarjani, H.T. Stokes, S. Nayeb Sadeghi, Y. Liang, B. Timalsina, H. Meng, J. Shiomi, B. Liao, R. Sun; ALADYN: a set of Anharmonic LAttice DYNamics codes to compute thermodynamic and thermal transport properties of crystalline solids
  [`Link <https://arxiv.org/abs/2501.02113>`__]

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
