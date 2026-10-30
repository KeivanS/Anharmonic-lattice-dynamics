Download
========

The source is on GitHub: https://github.com/KeivanS/Anharmonic-lattice-dynamics
(in a browser: *Code* → *Download ZIP*).

From a terminal, clone the repository once (this also brings SC-SNAPS, in ``SC_SNAPS/``):

.. code-block:: bash

   git clone https://github.com/KeivanS/Anharmonic-lattice-dynamics.git
   cd Anharmonic-lattice-dynamics

To update it later, from inside that folder:

.. code-block:: bash

   git pull
   make            # rebuild

(:doc:`install` describes the build.)

If you have edited files of your own, ``git pull`` keeps them; it stops, changing nothing, when the
update touches a file you edited. Then ``git stash``, ``git pull``, ``git stash pop`` puts your
edits aside, updates, and puts them back.
