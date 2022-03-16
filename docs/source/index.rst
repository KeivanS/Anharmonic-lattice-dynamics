Anharmonic-LAttice-DYNamics (ALADYN)
====================================

This repository contains 5 sets of stand-alone codes:
- **FOrce Constant EXtraction (FOCEX)** to extract force constants from force-displacements data, the output of which can be used as input to the following codes;
- **Self-CoOnsistent Phonon (SCOP)** to calculate the state of equilibrium of a crystal at a given temperature and pressure solely from the set of force constants parameters of the high-symmetry phase;
- **THERMAl CONDuctivity (THERMACOND)** to calculate the phonon spectrum, group velocities and lifetimes, and lattice thermal conductivity based on Boltzmann equation theory;
- **ANharmonic FOrce-field MOlecular Dynamics (ANFOMOD)** to perform molecular dymics from the above-developed force field;
- **INterface FOrce-field CONverter (INFOCON)** to convert extracted force constants by one of the codes FOCEX, ShengBTE, ALAMODE, and PHONOPY the format needed as input for the phonon calculation by another of these codes so that these codes become intercompatible.

.. note::

   This project is under active development.

Contents
--------

.. toctree::
