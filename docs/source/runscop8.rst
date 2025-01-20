Running SCOP8
=============

.. Self-CoOnsistent Phonon (SCOP8)
.. --------------------------------------
.. role:: raw-math(raw)
  :format: latex html

.. only:: html
  :math:`\\require{mediawiki-texvc}`

Workflow of SCOP8
-----------------
The workflow of the SCOP8 code is outlined below to provide a general overview of its structure and execution as a black box. While many intermediate steps in this process may not be essential for running the code, those that are critical are highlighted in the workflow diagram and subsequently explained in detail.

.. image:: Ref/SCOP8_FLowchart.png
 :width: 300
 :align: center

Description of input and output files
-------------------------------------

* **Preparing input files**

	In order to run SCOP8, you need 9 or 10 files (depends on the choice).
	
	The first group of input file(s) required to run SCOP8 are ``default.params`` and ``dielectric.params``.
	Those are simply the same input file(s) required to run FOCEX, they should be kept consistent with the FOCEX run. Refer to 'Running FOCEX' for further details about these files.
    
	The second group of input file(s) required to run SCOP8 are ``structure.params``, ``lat_fc.dat``, ``fc2.dat``, ``fc3.dat`` and ``fc4.dat``.  
	Those are generated from running FOCEX beforehand, they are raw input and should be left as is. Refer to 'Running FOCEX' for further details about these files.

	The third group of input file(s) required to run SCOP8 are ``control.params``, ``params.phon``, ``kpbs.params`` and ``targetInitialize.dat``(depends)
	Those are input files that uniquely for SCOP8 and can be customized by users, further details will be explained below

	* **control.params**

	The control.params is main 'control file' with several options.

	.. code-block:: python

		.False. # whether to start with random value for variational parameter strain tensor and atomic displacements. 
		.True. # whether to inherit previous running result. If set to 'True', this means that this run will start with variational parameters values from the file 'targetInitialize.dat'. So be sure to have this file in current path while turning this option on.
		.False. # whether include pressure factor, the stress tensor matrix can be customized in line 12 13 & 14
		.False. # whether force ASR in K after each Broyden iteration. This option should generally be left off since acoustic sum rule should be automatically preserved in K given proper input fc2, fc3 and fc4 values.
		.False. # whether use high temperature limit in <YY> calculation when dealing with soft modes. If this is turned on, the default strategy of manually shifting all negative eigenvalues up will be skipped, instead the high temperature limit will be employed for and only for negative eigenvalues.
		54321 # only matters if set the 1st line to .True., if you want random start, this is the seed
		0.02 # only matters if set the 1st line to .True., range for random start
		0.0 0.0 0.0 # only matters if set the 1st line to .True., center for random start
		300 # max_iteration number
		-0.006 # pmix parameter here, usually -0.006
		10 # temperature in Kelvin
		25 0.0 0.0
		0.0 25 0.0 
		0.0 0.0 25 # set stress tensor, 'pressure' flag hard coded
		3 #this is the number of variational parameters you want to fix, this number can be set to negative, then it means how many of variational parameters you want to free
		1 2 3 #variational parameters indexes list that you want to be fixed/free, the total number of indexes should match the value above


	*  **params.phon**

	The ``params.phon`` is another control file but mostly deprecated. The only really meaningful line is the first line which defines the size of k point mesh.
	This file contains also parameters for DOS calculation.

	.. code-block:: python

		18 18 18   .false.   kmesh  , if true calculate only dos-bs
		1 6    # lamin,lamax
		0 0 0 0.5 0.5 0.5 2300 0.1300 0.4000    # shift_xyz, (fine mesh in IBZ) 0<sx<sy<sz
		300  5    # wmesh and wmax (in cm^-1)
		4.0  10.0     # width of gaussian broadening for DOS (in cm^_1), imaginary part of omega(lifetimes)
		1d-11      # tau0(s)  added to the inverse of the relaxation time
		.false.    # verbosity
		0 0 0      # wshift in 3 directions; it pushes the phonon frequencies up
		380 380 1    # tmin, tmax, ntemp
		1   0  0  0 # iterative_BTE,readv3,usetetra: read/write if it is=1; collect if it's 2 ; use svd
		1             # ncpu;  below is path to v33sq.dat files  
		/
		100 4    # max_iter, n_dig_accuracy, deprecated
		0.00000001  #v3_threshold, used for checking, deprecated
		0      # 1 for classical and 0 for quantum, deprecated
		0 0 0 0     # for cross setion if=1
		0      # for the 3-phonon matrix elt if=1
		1      # scale length

	* **kbps.params**

	The ``kpbs.params`` file define the q-point path for band dispersion plot. 

	.. code-block:: python

		0 # The first line can be set to 0 or 1 where 0 means conventional lattice and 1 means primitive lattice. 
		40 # The second line is the number of k points along each direction
		7 # The third line is the number of directions for the band plot.
		G 0 0.0001 0.0001 # The following lines should be input with format of special q point label(G, K, L, X, etc) followed by its 3d coordinates.
		K 0 0.75 0.75  
		X 0 1 1  
		G 0.9999 1 1    
		L 0.5 0.5 0.5 
		W 0. 0.5 1  
		X 0 0 1   
		G 0 0. 0.0001    gamma 


	* **targetInitialize.dat**

	This is the main output file that lists final values for all the variational parameters. It also serves as a customizable input file if user choose to set the 2nd line of control.params to True.

	.. code-block:: python

		0.0000000000000000        0.0000000000000000        0.0000000000000000     # the 'relaxed' tau_0, usually (0,0,0) since the center atom is fixed and used as a reference point
		7.3258075924813355E-011   7.3241067962641583E-011   7.3225737274375089E-011 # the 'relaxed' tau_1, if the lattice has more than two atoms (types), there will be more lines below
		9.7839998026831273E-004   3.2204510841042777E-012   3.2167153648262385E-012 # the 1st row of 'relaxed' strain tensor a.k.a eta_xx, eta_xy, eta_xz
		3.2204510841042777E-012   9.7839998014266562E-004   3.2126854822356125E-012 # the 2nd row of 'relaxed' strain tensor a.k.a eta_yx, eta_yy, eta_yz
		3.2167153648262385E-012   3.2126854822356125E-012   9.7839997972007462E-004 # the 3rd row of 'relaxed' strain tensor a.k.a eta_zx, eta_zy, eta_zz
		1           1           1           1           1   11.484396453112186   # this line and all the lines below corresponds to trial force constants, the columns are: index, atom_1, xyz_1, atom_2, xyz_2, K  
		2           1           1           2           2   11.484396453112325     
		3           1           1           3           3   11.484396453113003     
		4           2           2           1           1   11.484396453112071     
		5           2           2           2           2   11.484396453112852     
		6           2           2           3           3   11.484396453113003     
		7           1           2           1           1  -2.8395049908152701     
		8           2           1           1           1  -2.8395049908151790     
		9           1           2           2           2  -2.8395049908153513     
		10           2           1           2           2  -2.8395049908153425     
		11           1           2           3           3  -2.8395049908153900
		...
		23312           2         839           3           3   0.0000000000000000     
		23313           1        1160           2           2   0.0000000000000000     
		23316           2         817           2           2   0.0000000000000000     

* **Explaining output files**

	Most of the output files can be ignored, since they are for checking or logging, or simply legacy output that should be deprecated.
	Import output files will be briefly explained below.

	* **convergence.dat** 
	This file keeps a record of free energy value and L1 norm of all gradients at each iteration. So you can monitor how well the Broyden loops converge.

	.. code-block:: python

 		iteration,  free energy,                L1 norm    
			1 ,   7.7695601608371171E-003 ,  0.12192370206947493     
			2 ,   5.4742390823587731E-003 ,  0.12192370153985427     
			3 ,   1.9890803839332921E-003 ,  0.12192370154367714     
			4 ,   2.0003365953782347E-003 ,  0.12192370105296366     
			5 ,   1.0012513345730565E-005 ,  0.12192370105384057     
			6 ,   9.8580778483004610E-006 ,  0.12192370105387422     
			7 ,   9.8873791186356400E-006 ,  0.12192370104103692     
			8 ,   2.8162657844481114E-007 ,  0.12192370104041715     
			9 ,   1.8986235130117904E-007 ,  0.12192370104011624     
			10 ,   1.5464059224802259E-007 ,  0.12192370104020522     
			11 ,   1.0357871383178676E-007 ,  0.12192370104049430     
			12 ,   1.7166077234971523E-008 ,  0.12192370104056582     
			13 ,   5.7005150331804107E-009 ,  0.12192370104059483     
			14 ,   5.4702303408328802E-009 ,  0.12192370104113985     
			15 ,   6.7558572080945439E-011 ,  0.12192370104115086     

	* **Dispersion.dat**

	This file is for ploting the phonon dispersion as the file name suggests. Simily we have output files like ``dos_gauss.dat``, ``dos_tet.dat``, etc for the corresponding plot purposes.

	.. code-block:: python

	  q point index, band1, band2, band3, ... 
		1   2.9242458690853157E-003   2.9426593567020258E-003   2.9426606011758892E-003   15.200445394847394        15.200445395546224        15.200445395547092     
		2  0.26453273368343705       0.29992250241135027       0.46473813700486238        15.193704422889539        15.194830030019409        15.197208518698289     
		3  0.52727363740652622       0.59778066306805544       0.92753330396543299        15.173806732746305        15.178012413287949        15.187544762083835     
		4  0.78652244152886719       0.89163644909276552        1.3865038106147869        15.141710311252552        15.150076991747932        15.171592760596656     
		5   1.0406193052774069        1.1796769554150286        1.8398093401632922        15.098953274317608        15.111161267282284        15.149577082937030     
			...
		298   3.1829648451308912        3.2047348694808524        11.097412640591298        12.194293467387048        14.469149408840680        14.472597987374058     
		299   3.1756241291339333        3.1811019189677405        11.104582714239395        12.197751480924026        14.473360701151645        14.474226063172505     
		300   3.1731715122948594        3.1731715122953355        11.106976501294767        12.198906836802358        14.474768794358322        14.474768794359123     

	* **eigenvalues.dat**

	This file usually will print nothing unless there are negative eigenvalues a.k.a soft modes. 

	* **GradientF.dat**

	This file is the history of every variational parameter values and its free energy gradients values at every iteration, the FinalGradientF.dat is just these of the last iteration for convenience, since GradientF.dat is usually a large file with too many iterations.

	.. code-block:: python

		# the file serves as a logfile, so it's well printed and self-explainary. It looks like this below:
		current interation #:           1
		temperature=   40.000000000000000     
		F0= ( 0.12152291457977416     ,  0.0000000000000000     )
		V0= (  6.0941496242618617E-002,  0.0000000000000000     )
		free energy=  0.12192370398307743     
		=============GradientF:trial fc2====================
		largest gradient=    3.7915793395093544E-004
		||Atomic deviation u_tau(:)||
				1 x variable=   0.0000000000000000       gradient=   2.9773399986220660E-011   FIXED   
				1 y variable=   0.0000000000000000       gradient=   2.9774619876295032E-011   FIXED   
				1 z variable=   0.0000000000000000       gradient=   2.9774363665769146E-011   FIXED   
				2 x variable=   6.4873485806078854E-011  gradient=  -2.9773398395431535E-011   FREE    
				2 y variable=   6.4873383295777969E-011  gradient=  -2.9774613452670168E-011   FREE    
				2 z variable=   6.4873468950255914E-011  gradient=  -2.9774355445975634E-011   FREE    
		||Strain Tensor||
		xx variable=   1.2673966516492396E-003  gradient=  -3.7584995895641470E-004   FREE    
		xy variable=  -4.1761063898985007E-011  gradient=  -5.5646293115543233E-011   FREE    
		xz variable=  -4.1761065621560632E-011  gradient=  -5.5646710083691172E-011   FREE    
		yx variable=  -4.1761063898985007E-011  gradient=  -5.5646293115543233E-011   FREE    
		yy variable=   1.2673966546147906E-003  gradient=  -3.7584990147655333E-004   FREE    
		yz variable=  -4.1761040526752187E-011  gradient=  -5.5642905807692281E-011   FREE    
		zx variable=  -4.1761065621560632E-011  gradient=  -5.5646710083691172E-011   FREE    
		zy variable=  -4.1761040526752187E-011  gradient=  -5.5642905807692281E-011   FREE    
		zz variable=   1.2673966451871321E-003  gradient=  -3.7585008745702762E-004   FREE    
		||Force Constants||
				1 x           1 x variable=   13.135477166947316       gradient=   3.7915789332831906E-004   FREE    
				1 y           1 y variable=   13.135477166947151       gradient=   3.7915787513753685E-004   FREE    
				1 z           1 z variable=   13.135477166951206       gradient=   3.7915793392784281E-004   FREE    
				2 x           2 x variable=   13.135477166947448       gradient=   3.7915789326081750E-004   FREE    
				2 y           2 y variable=   13.135477166947071       gradient=   3.7915787517572852E-004   FREE    
				2 z           2 z variable=   13.135477166951164       gradient=   3.7915793395093544E-004   FREE    
				...
		...
		...
		current interation #:           999
		...

	* **FinalGradientF.dat**

	This is basically the last part of the previous file, for a quick look at what the final results are.

	* **output.txt** and **logfile.txt**

	These two files contain runtime info that only should be referred when try to debug the code. Otherwise ignore them.

	* **result.txt**

	This file contains most post-process results, such as final free energy, gruneisen, specific heat and elastic constants.
	It might look different as more post-process subroutines can be added later on.
	Some thermodynamics properties of interest might also be printed out in result.txt

	.. code-block:: python
		
		first Unpertubed Free Energy F0= ( 0.12152291457977416     ,  0.0000000000000000     )
		final translational vector =   0.0000000000000000        2.6988099999999999        2.6988099999999999        2.6988099999999999        0.0000000000000000        2.6988099999999999        2.6988099999999999        2.6988099999999999        0.0000000000000000     
		final F0 =  ( 0.12152325743018681     ,  0.0000000000000000     )
		final Free Energy F=F0+<V-V0> =  ( 0.12192370104115086     ,  0.0000000000000000     )
		Temperature   40.000000000000000     
		Current Volume is:   39.464183071187328     
		my calculated gruneisen: -0.94758846679605901     
		my calculated specific heat:   1.4541875815899401     
		my calculated bulk modulus:  0.66968252261197314     
		calculated beta =   -1.0806246818743969E-006
		elastic
		-7.2854364387643216       -5.3617679414071420       -5.3617679414070265        8.8636066792102938E-009   1.2069888990825881E-009   1.2068871752796343E-009
		-5.3617679414070869       -7.2854364387640818       -5.3617679414068569        1.2068554253301697E-009   8.8636579568164597E-009   1.2068605589486403E-009
		-5.3617679414071437       -5.3617679414070549       -7.2854364387640889        1.2069699917876204E-009   1.2069495447408507E-009   8.8635875559109161E-009
		-8.3631174514409110E-009  -4.7071115546290751E-009  -4.7071123636814083E-009  -14.881121710429831        1.2393144638154903E-010   1.2404539918885701E-010
		-4.7070043759087784E-009  -8.3630394947033851E-009  -4.7070017343739520E-009   1.2392663922588006E-010  -14.881121710429678        1.2398792816686518E-010
		-4.7071008543127752E-009  -4.7070844606782160E-009  -8.3632129559818321E-009   1.2406723758800231E-010   1.2405596603844226E-010  -14.881121710429220     

		compliance
		-1.3082394889446023       0.55462664638390591       0.55462664638388204       -6.8925908783530570E-010   2.6922679641009340E-010   2.6922979880380062E-010
		0.55462664638390735       -1.3082394889446218       0.55462664638385151        2.6923738137719943E-010  -6.8925938322747944E-010   2.6923313072073446E-010
		0.55462664638388670       0.55462664638387926       -1.3082394889446018        2.6922303958373513E-010   2.6923172291170923E-010  -6.8926030193366936E-010
		3.8435147375417462E-010  -7.3318969926278909E-011  -7.3318868646515789E-011 -0.24081093946118309       -2.0054972257712852E-012  -2.0073412455396363E-012
		-7.3321363602604091E-011   3.8435273817860167E-010  -7.3321694278298720E-011  -2.0054194348880938E-012 -0.24081093946118556       -2.0064112315912288E-012
		-7.3322430162109725E-011  -7.3324482369402487E-011   3.8436130858012292E-010  -2.0076946413059313E-012  -2.0075122415880799E-012 -0.24081093946119300     



Example
-------

still writing.