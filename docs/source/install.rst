Installation
============

Requirement
-----------

* Fortran Compiler (GCC or Intel)

Install
-------

.. _focex-install:

FOCEX
^^^^^^

* To install the FOCEX code, clone the repository `https://github.com/KeivanS/Anharmonic-lattice-dynamics.git` or download the ALADYN code and follow the instructions
* Edit the Makefile, generally you do not need to change the Makefile. In case if you need to change the compiler flags adjust the FLAGS option inside **Anharmonic-lattice-dynamics/FOCEX/Makefile** by uncommenting or adding the available option

``FLAGS= -fcheck=all -fdec-math #-Wall #-fbounds-check #-fbacktrace #-ffpe-trap=invalid,zero,overflow,underflow,inexact #-inline-debug-info #-p #-pedantic #-std=f95``

* From the **Anharmonic-lattice-dynamics/FOCEX/Makefile**, simply ``make`` and this will create the binary ``focex.x``



