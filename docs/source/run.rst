Running FOCEX
==============

.. FOrce Constant EXtraction (FOCEX)
.. ---------------------------------

.. role:: raw-math(raw)
	:format: latex html

.. only:: html
	:math:`\\require{mediawiki-texvc}`

Workflow of FOCEX
-----------------
The  workflow of the FOCEX code is shown below to give a general understanding of how the code is structured and executed as a black box. Many of the intermediate steps in this process flow may even not be necessary to actually run the code. The one which are important are highlighted on the workflow diagram and are later explained in detail. 

.. image:: ./WorkFlow-FOCEX-Website.svg
  :width: 600
  :align: center

General steps for running FOCEX to get the force constants are given below:

Description of input and output files
-------------------------------------

* **Preparing FOCEX data files**

    FOCEX code accepts DFT force displacement files prepared by any software like VASP, or Quantum Espresso. We provide utility tools to extract forces and displacements in a supercell from VASP OUTCAR or Quantum Espresso output files. All the generated displacement-force data for each atomic snapshot are concatenated into a signle called FORCEDISP1. The equilibrium atomic positions in the supercell should be in a file called POSCAR1 which has the exact format of the VASP POSCAR file. If there are other supercells, the equilibrium positions and displacement-force data should be stored in POSCAR2 and FORCEDISP2, etc... 

.. collapse:: utility tools

    Something small enough to escape casual notice.

* **Preparing input files**

    The input file(s) required to run FOCEX are ``structure.params``, ``dielectric.params``, ``kpbs.params``, POSCAR1 and FORCEDISP1.  

.. collapse:: structure.params

    something about structure.params file

.. collapse:: dielectric.params

    something about dielectric.params file

.. collapse:: POSCAR1

    something about POSCAR1

.. collapse:: FORCEDISP1

    The format of this file is : First line is a comment and should contain " # POSITION ...", second line should contain an integer followed by the total energy of the supercell snapshot in eV. Units for positions or displacements should be in Angstrom, and forces in eV/Ang. 

.. collapse:: kpbs.params

    The first line contains a flag. If 0, the kpoints are given in reduced units of the primitive vectors of the reciprocal space, else they should be in reduced coordinates of the conventional lattice vectors of the reciprocal space.
The second line contains the number of kpoints along each direction, The third line contains the number of direction paths. The following lines contain the name of the special point followed by the 3 reduced components of the special kpoint in units of primitive (if flag=0) or conventional (if flag is non-zero) vectors of the reciprocal lattice. 
	
* **output files**

	* fc1.dat and fc2_irr.dat

	* fc2.dat and fc2_irr.dat

	* fc3.dat and fc3_irr.dat

	* lat_fc.xyz 

	* log***.dat
 
Example of Running Ge
---------------------
FOrce Constant Extraction (FOCEX) is a code to extract force constants from force-displacements data, the output of which can be used as input to the following codes.
The installation of FOCEX has to be done before using it and is given in section :ref:`focex-install`. This code, FOCEX (FOrce Constant EXtraction) included in ALADYN (Anharmonic LAttice DYNamics) employs the
force constant calculation, 2nd, 3rd and 4th order to be latter used for other thermodynamical properties. The installation of FOCEX is simple and just require the
fortran compiler.

An example of **Ge** is provided inside **FOCEX/example** which contains the needed input files described above, and the FOCEX output. In this
calculation a single Ge atom in the 2x2x2 Ge supercell (64 atoms) is displaced by 4% to evaluate the forces. The force-displacement data is stored in ``FORCEDISP1`` file.
The equilibrium coordinates are in ``POSCAR1`` in the older VASP POSCAR format as below. 

.. code-block:: python

  Ge8 # this is a comment
  1.0  # scale factor
  11.5257244110         0.0000000000         0.0000000000  # supercell
  0.0000000000        11.5257244110         0.0000000000
  0.0000000000         0.0000000000        11.5257244110
  64  # number of atoms in the supercell
  Cartesian
  2.881431103         0.000000000         2.881431103
  2.881431103         0.000000000         8.644293308
  2.881431103         5.762862206         2.881431103
  2.881431103         5.762862206         8.644293308
  8.644293308         0.000000000         2.881431103
  8.644293308         0.000000000         8.644293308
  8.644293308         5.762862206         2.881431103
  8.644293308         5.762862206         8.644293308
  1.440715551         1.440715551         1.440715551
  1.440715551         1.440715551         7.203577757
  1.440715551         7.203577757         1.440715551
  1.440715551         7.203577757         7.203577757
  7.203577757         1.440715551         1.440715551
  7.203577757         1.440715551         7.203577757
  7.203577757         7.203577757         1.440715551
  7.203577757         7.203577757         7.203577757
  2.881431103         2.881431103         0.000000000
  2.881431103         2.881431103         5.762862206
  2.881431103         8.644293308         0.000000000
  2.881431103         8.644293308         5.762862206
  8.644293308         2.881431103         0.000000000
  8.644293308         2.881431103         5.762862206
  8.644293308         8.644293308         0.000000000
  8.644293308         8.644293308         5.762862206
  1.440715551         4.322146654         4.322146654
  1.440715551         4.322146654        10.085008860
  1.440715551        10.085008860         4.322146654
  1.440715551        10.085008860        10.085008860
  7.203577757         4.322146654         4.322146654
  7.203577757         4.322146654        10.085008860
  7.203577757        10.085008860         4.322146654
  7.203577757        10.085008860        10.085008860
  0.000000000         0.000000000         0.000000000
  0.000000000         0.000000000         5.762862206
  0.000000000         5.762862206         0.000000000
  0.000000000         5.762862206         5.762862206
  5.762862206         0.000000000         0.000000000
  5.762862206         0.000000000         5.762862206
  5.762862206         5.762862206         0.000000000
  5.762862206         5.762862206         5.762862206
  4.322146654         1.440715551         4.322146654
  4.322146654         1.440715551        10.085008860
  4.322146654         7.203577757         4.322146654
  4.322146654         7.203577757        10.085008860
  10.085008860         1.440715551         4.322146654
  10.085008860         1.440715551        10.085008860
  10.085008860         7.203577757         4.322146654
  10.085008860         7.203577757        10.085008860
  0.000000000         2.881431103         2.881431103
  0.000000000         2.881431103         8.644293308
  0.000000000         8.644293308         2.881431103
  0.000000000         8.644293308         8.644293308
  5.762862206         2.881431103         2.881431103
  5.762862206         2.881431103         8.644293308
  5.762862206         8.644293308         2.881431103
  5.762862206         8.644293308         8.644293308
  4.322146654         4.322146654         1.440715551
  4.322146654         4.322146654         7.203577757
  4.322146654        10.085008860         1.440715551
  4.322146654        10.085008860         7.203577757
  10.085008860         4.322146654         1.440715551
  10.085008860         4.322146654         7.203577757
  10.085008860        10.085008860         1.440715551
  10.085008860        10.085008860         7.203577757

Here, only the type of atom is not present in ``POSCAR1`` as compared to the new format of VASP POSCAR file. Similarly, ``FORCEDISP1`` is a force-displacement data format
accepted by FOCEX code and its format for example in the case of Ge is given below.

.. code-block:: python

  # POSITION (ang)     TOTAL FORCE (eV/Ang)  
     1       -289.18629538 =t, Etot(eV)     # snapshot #1
   2.8929600000000000        0.0000000000000000        2.8814299999999999      -0.11758299999999999       -0.0000000000000000       -0.0000000000000000
   2.8814299999999999        0.0000000000000000        8.6442899999999998        4.9600000000000002E-004  -0.0000000000000000       -0.0000000000000000
   2.8814299999999999        5.7628599999999999        2.8814299999999999        4.9600000000000002E-004  -0.0000000000000000       -0.0000000000000000
   2.8814299999999999        5.7628599999999999        8.6442899999999998       -4.5640000000000003E-003  -0.0000000000000000       -0.0000000000000000
   8.6442899999999998        0.0000000000000000        2.8814299999999999        3.2899999999999997E-004  -0.0000000000000000       -0.0000000000000000
   8.6442899999999998        0.0000000000000000        8.6442899999999998       -1.5699999999999999E-004  -0.0000000000000000       -0.0000000000000000
   8.6442899999999998        5.7628599999999999        2.8814299999999999       -1.5699999999999999E-004  -0.0000000000000000       -0.0000000000000000
   8.6442899999999998        5.7628599999999999        8.6442899999999998        2.6699999999999998E-004  -0.0000000000000000       -0.0000000000000000
   1.4407200000000000        1.4407200000000000        1.4407200000000000        2.8677000000000001E-002  -1.9474999999999999E-002   1.9474999999999999E-002
   1.4407200000000000        1.4407200000000000        7.2035799999999997       -4.8099999999999998E-004  -7.1400000000000001E-004   3.2800000000000000E-004
   1.4407200000000000        7.2035799999999997        1.4407200000000000       -4.8099999999999998E-004  -3.2800000000000000E-004   7.1400000000000001E-004
   1.4407200000000000        7.2035799999999997        7.2035799999999997        2.4350000000000001E-003  -4.0000000000000003E-005   4.0000000000000003E-005
   7.2035799999999997        1.4407200000000000        1.4407200000000000       -1.8400000000000000E-004  -3.3000000000000000E-004   3.3000000000000000E-004
   7.2035799999999997        1.4407200000000000        7.2035799999999997        8.3999999999999995E-005  -4.3999999999999999E-005   3.4000000000000000E-005
   7.2035799999999997        7.2035799999999997        1.4407200000000000        8.3999999999999995E-005  -3.4000000000000000E-005   4.3999999999999999E-005
   7.2035799999999997        7.2035799999999997        7.2035799999999997       -2.1599999999999999E-004  -3.4400000000000001E-004   3.4400000000000001E-004
   2.8814299999999999        2.8814299999999999        0.0000000000000000       -4.0619999999999996E-003   7.2599999999999997E-004  -7.2599999999999997E-004
   2.8814299999999999        2.8814299999999999        5.7628599999999999       -4.0780000000000000E-003  -7.3099999999999999E-004  -7.3099999999999999E-004
   2.8814299999999999        8.6442899999999998        0.0000000000000000       -4.0780000000000000E-003   7.3099999999999999E-004   7.3099999999999999E-004
   2.8814299999999999        8.6442899999999998        5.7628599999999999       -4.0619999999999996E-003  -7.2599999999999997E-004   7.2599999999999997E-004
   8.6442899999999998        2.8814299999999999        0.0000000000000000        1.0700000000000000E-004   6.0000000000000002E-005  -6.0000000000000002E-005
   8.6442899999999998        2.8814299999999999        5.7628599999999999        1.0300000000000000E-004  -5.7000000000000003E-005  -5.7000000000000003E-005
   8.6442899999999998        8.6442899999999998        0.0000000000000000        1.0300000000000000E-004   5.7000000000000003E-005   5.7000000000000003E-005
   8.6442899999999998        8.6442899999999998        5.7628599999999999        1.0700000000000000E-004  -6.0000000000000002E-005   6.0000000000000002E-005
   1.4407200000000000        4.3221499999999997        4.3221499999999997       -4.8099999999999998E-004   3.2800000000000000E-004  -7.1400000000000001E-004
   1.4407200000000000        4.3221499999999997        10.085010000000000        2.4350000000000001E-003   4.0000000000000003E-005  -4.0000000000000003E-005
   1.4407200000000000        10.085010000000000        4.3221499999999997        2.8677000000000001E-002   1.9474999999999999E-002  -1.9474999999999999E-002
   1.4407200000000000        10.085010000000000        10.085010000000000       -4.8099999999999998E-004   7.1400000000000001E-004  -3.2800000000000000E-004
   7.2035799999999997        4.3221499999999997        4.3221499999999997        8.3999999999999995E-005   3.4000000000000000E-005  -4.3999999999999999E-005
   7.2035799999999997        4.3221499999999997        10.085010000000000       -2.1599999999999999E-004   3.4400000000000001E-004  -3.4400000000000001E-004
   7.2035799999999997        10.085010000000000        4.3221499999999997       -1.8400000000000000E-004   3.3000000000000000E-004  -3.3000000000000000E-004
   7.2035799999999997        10.085010000000000        10.085010000000000        8.3999999999999995E-005   4.3999999999999999E-005  -3.4000000000000000E-005
   0.0000000000000000        0.0000000000000000        0.0000000000000000        1.6290000000000000E-003  -7.2499999999999995E-004   1.4820000000000000E-003
   0.0000000000000000        0.0000000000000000        5.7628599999999999        1.6290000000000000E-003   7.2499999999999995E-004  -1.4820000000000000E-003
   0.0000000000000000        5.7628599999999999        0.0000000000000000        4.2000000000000002E-004  -5.8000000000000000E-005  -8.8500000000000004E-004
   0.0000000000000000        5.7628599999999999        5.7628599999999999        4.2000000000000002E-004   5.8000000000000000E-005   8.8500000000000004E-004
   5.7628599999999999        0.0000000000000000        0.0000000000000000        1.6550000000000000E-003  -7.3399999999999995E-004  -1.5030000000000000E-003
   5.7628599999999999        0.0000000000000000        5.7628599999999999        1.6550000000000000E-003   7.3399999999999995E-004   1.5030000000000000E-003
   5.7628599999999999        5.7628599999999999        0.0000000000000000        4.2000000000000002E-004  -5.7000000000000003E-005   8.8000000000000003E-004
   5.7628599999999999        5.7628599999999999        5.7628599999999999        4.2000000000000002E-004   5.7000000000000003E-005  -8.8000000000000003E-004
   4.3221499999999997        1.4407200000000000        4.3221499999999997        2.8958999999999999E-002   2.0014000000000001E-002   2.0014000000000001E-002
   4.3221499999999997        1.4407200000000000        10.085010000000000       -4.8400000000000000E-004   7.0799999999999997E-004   3.3000000000000000E-004
   4.3221499999999997        7.2035799999999997        4.3221499999999997       -4.8400000000000000E-004   3.3000000000000000E-004   7.0799999999999997E-004
   4.3221499999999997        7.2035799999999997        10.085010000000000        2.4480000000000001E-003   4.0000000000000003E-005   4.0000000000000003E-005
   10.085010000000000        1.4407200000000000        4.3221499999999997       -1.8000000000000001E-004   3.2899999999999997E-004   3.2899999999999997E-004
   10.085010000000000        1.4407200000000000        10.085010000000000        7.7999999999999999E-005   3.4999999999999997E-005   2.8000000000000000E-005
   10.085010000000000        7.2035799999999997        4.3221499999999997        7.7999999999999999E-005   2.8000000000000000E-005   3.4999999999999997E-005
   10.085010000000000        7.2035799999999997        10.085010000000000       -2.1200000000000000E-004   3.4699999999999998E-004   3.4699999999999998E-004
   0.0000000000000000        2.8814299999999999        2.8814299999999999        1.6290000000000000E-003  -1.4820000000000000E-003   7.2499999999999995E-004
   0.0000000000000000        2.8814299999999999        8.6442899999999998        4.2000000000000002E-004   8.8500000000000004E-004   5.8000000000000000E-005
   0.0000000000000000        8.6442899999999998        2.8814299999999999        1.6290000000000000E-003   1.4820000000000000E-003  -7.2499999999999995E-004
   0.0000000000000000        8.6442899999999998        8.6442899999999998        4.2000000000000002E-004  -8.8500000000000004E-004  -5.8000000000000000E-005
   5.7628599999999999        2.8814299999999999        2.8814299999999999        1.6550000000000000E-003   1.5030000000000000E-003   7.3399999999999995E-004
   5.7628599999999999        2.8814299999999999        8.6442899999999998        4.2000000000000002E-004  -8.8000000000000003E-004   5.7000000000000003E-005
   5.7628599999999999        8.6442899999999998        2.8814299999999999        1.6550000000000000E-003  -1.5030000000000000E-003  -7.3399999999999995E-004
   5.7628599999999999        8.6442899999999998        8.6442899999999998        4.2000000000000002E-004   8.8000000000000003E-004  -5.7000000000000003E-005
   4.3221499999999997        4.3221499999999997        1.4407200000000000       -4.8400000000000000E-004  -3.3000000000000000E-004  -7.0799999999999997E-004
   4.3221499999999997        4.3221499999999997        7.2035799999999997        2.4480000000000001E-003  -4.0000000000000003E-005  -4.0000000000000003E-005
   4.3221499999999997        10.085010000000000        1.4407200000000000        2.8958999999999999E-002  -2.0014000000000001E-002  -2.0014000000000001E-002
   4.3221499999999997        10.085010000000000        7.2035799999999997       -4.8400000000000000E-004  -7.0799999999999997E-004  -3.3000000000000000E-004
   10.085010000000000        4.3221499999999997        1.4407200000000000        7.7999999999999999E-005  -2.8000000000000000E-005  -3.4999999999999997E-005
   10.085010000000000        4.3221499999999997        7.2035799999999997       -2.1200000000000000E-004  -3.4699999999999998E-004  -3.4699999999999998E-004
   10.085010000000000        10.085010000000000        1.4407200000000000       -1.8000000000000001E-004  -3.2899999999999997E-004  -3.2899999999999997E-004
   10.085010000000000        10.085010000000000        7.2035799999999997        7.7999999999999999E-005  -3.4999999999999997E-005  -2.8000000000000000E-005
  # POSITION (ang)     TOTAL FORCE (eV/Ang)  
     1       -289.18629538 =t, Etot(eV)     # snapshot #2
   2.8929600000000000        0.0000000000000000        2.8814299999999999      -0.11758299999999999       -0.0000000000000000       -0.0000000000000000
   2.8814299999999999        0.0000000000000000        8.6442899999999998        4.9600000000000002E-004  -0.0000000000000000       -0.0000000000000000
...


The first line in ``FORCEDISP1`` is the header and should contain the word POSITION. The second line
consists of the energy of structure in electron volt and the lines after second are positions (first three columns, :math:`x`, :math:`y` and :math:`z`) and forces
(the last three columns :math:`F_x`, :math:`F_y` and :math:`F_z` are in :math:`eV/{\\A}`) respectively. If there are many force-displacement snapshots of the structure, then the other snapshots follow these lines in the same format. There is a tool and a shell script for converting the VASP outcar or QE outputfile into the FORCEDISP1 format for the FOCEX code. It is available in **utility** folder inside FOCEX. To convert the VASP outcar into FORCEDISP1 user need to execute ``./process_dft.sh name_of_vasp_directory(s)`` or ``./process_dft.sh name_of_vasp_file`` shell script within utility folder. Here, ``name_of_vasp_directory(s)`` is the multiple directory containing VASP runs or QE runs or ``name_of_vasp_file`` is the OUTCAR file(s) for VASP or the outputfile for QE runs. This shell script will call the ``readoutcar.x`` or ``readpfpwscf.x`` taking as input the DFT output from VASP or QE. The shell script can be tailored as per your needs. As for the other inputs file of FOCEX, they are given below

``structure.params``

.. code-block:: python

  1 1 1 90 90 90          # a, b, c, alpha, beta, gamma of the conventional cell
  0 .5 .5   .5 0 .5   .5 .5 0 # reduced coordinates of primitive lattice (in this case FCC) in terms of conventional lattice (in this case cubic)
  5.7628622055            # scale factor for lattice size
  1 1 1 1                 # include FC1234, 1st, 2nd, third and fourth order harmonic force constant(s) in the fitting process. 1 is to include and 0 is to not include
  1 1 0 0                 # invariances to impose, (translational, rotational, Huang) last is enforce inv using elimination
  0 300                   # temperature and whether or not implement it (do not implement if 0,2, or ..)
  1 .true.                # number  of FORCEDISPi files, verbosity
  1                       # type of atoms
  72.64                   # masses of each type of atoms
  Ge                      # names of atoms
  2                       # number of atoms of each type in the primitive cell
  1                       # flag for setting the range of FC2s (if 0 take default; else use below)
  5 5                     # number of shells for rank 2 (harmonic) for each atom if not default
  1 1                     # number of shells for rank 3 (cubic) for each atom (there is no default value)
  1 1                     # number of shells for rank 4 (quartic) for each atom (there is no default value)
  1 1 0 0 0               # atom index, type of atom, position x, position y, position z in units of a,b,c of the conventional lattice
  2 1 0.25 0.25 0.25      # atom index, type of atom, position x, position y, position z in units of a,b,c

The fitting is done using singular value decomposition based on the requested symmetry constraints and ``POSCAR1`` and ``FORCEDISP1`` i.e. by creating the force displacement matrix. 


The next input file is ``dielectric.params``. It is required to get the phonon dispersion (and eventually thermal conductivity using ``THERMACOND``). It consists of a flag which is zero if the Born charges are to be excluded. The second, third and fourth lines contain the dielectric constant tensor values which is written as follows in the example folder inside ``FOCEX``

``dielectric.params``

.. code-block:: python

  0               # do not include Born charges in the fitting
  2.5078 0.0  0.0 # for example, dielectric constant (necessary but not used if flag is zero) 
  0.0 2.5078 0.0
  0.0 0.0 2.5078
  0 0 0           # Born charge tensor of atom 1 in order it appears in the ``structure.params``
  0 0 0 
  0 0 0
  0 0 0           # Born charge tensor of atom 2 in order it appears in the ``structure.params``
  0 0 0
  0 0 0 

``kpbs.params``

.. code-block:: python

  1  # use direct coordinates of the conventional reciprocal cell
  30 # number of kpoints along each direction
  4  # number of directions
  G 0 0 0 
  K 0.75 0.75 0
  X 1 1 0
  G 0 0 0
  L 0.5 0.5 0.5

Now, put the ``POSCAR1``, ``FORCEDISP1`` , ``structure.params``, ``dielectric.params`` and ``kpbs.params`` in same directory, simply run ``focex.x``. After successful
run ``fc2.dat``, ``fc2_irr.dat``, ``fc3.dat``, ``fc3_irr.dat``, ``fc4.dat`` and ``fc4_irr.dat`` along with other output files and log file should be available. ``fc2.dat``, ``fc2_irr.dat`` are the fitted second order force constants in eV/Ang^2. The former contains all the pairs regardless of symmetry, and the latter contains the irreducible ones from which all the rest is contructed using crystal symmetry. Likewise ``fc3.dat``, ``fc3_irr.dat`` and ``fc4.dat``, ``fc4_irr.dat`` contain the third order and fourth order full and irreducible force constants in eV/Ang^3 and eV/Ang^4 respectively. Users are advised to look for more information in the log file ``log***.dat`` to know more details about the run.
 
