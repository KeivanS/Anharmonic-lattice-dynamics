!zoom id:https://virginia.zoom.us/j/4349248029
!!main program
!!mpi used for parallelisation

Program SCOP8

    use mpi
    use mpi_params

    IMPLICIT NONE

    call MPI_Init(mpi_err)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,mpi_err)
    call MPI_Comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

    ! For allgather
    if (.not. allocated(vals_displs)) allocate(vals_displs(nprocs))
    if (.not. allocated(vecs_displs)) allocate(vecs_displs(nprocs))
    if (.not. allocated(offsets)) allocate(offsets(nprocs))
    if (.not. allocated(vals_count)) allocate(vals_count(nprocs))
    if (.not. allocated(vecs_count)) allocate(vecs_count(nprocs))

    offsets=0
    vals_displs=0
    vecs_displs=0
    ! main subroutine, all-in-one
    CALL ThreeVariableRoutine

    ! In testing code this needs to be before every stop statement
    call MPI_Finalize(mpi_err) !MODIFY:comment on 11/03/2023

End Program SCOP8
!===============================================================================================
!===============================================================================================
!===============================================================================================

SUBROUTINE ThreeVariableRoutine
    USE VA_math
    USE force_update
    USE broy
    USE CG
    USE check

    IMPLICIT NONE

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~VARIABLES DECLARATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INTEGER :: direction1, direction2,atom1
    INTEGER i,j,k_number,k,l,m,eigen_number,n,ier,nv
    INTEGER t_atom1, t_atom2, t_xyz1, t_xyz2
    INTEGER mx,unit_number,unit_number2,unit_number3,unit_number4,unit_number5
    !-----for reciprocal part use-------
    INTEGER :: multiplier=4, nk,nband
    REAL(8) wmin,dk, dq
    REAL(8) err,ft, ft2, threshold
    REAL(8),ALLOCATABLE:: eigen_temp(:,:),integrate_dos(:),total_dos(:),afunc_dos(:),junk(:)
    !-----for iteration---------
    INTEGER :: method=1 !0=simple update, 1=Broyden, 2=Conjugate Gradient(not fully implemented yet)
    REAL(8) :: fret,fp,eps=1.e-10
    REAL(8),DIMENSION(:),ALLOCATABLE :: x,f,g,h,f_ex,x2,f2
    CHARACTER name*2
    !---------for test----------
    TYPE(fc2_index),DIMENSION(:),ALLOCATABLE :: ex_indiefc2_index
    REAL(8) :: temp, temp_check,temp_check2,volume,start_F, step, step_1, step_2
    INTEGER :: item
    COMPLEX(8) :: dumb
    INTEGER :: guess
    REAL(8) :: accept
    REAL(8) :: cell_vol
    INTEGER :: uio !UPDATE: FOCEX_ec
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~OPEN OUTPUT AND LOG FILES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    unit_number=33;unit_number2=34;unit_number3=35;unit_number4=36;unit_number5=37
    path_out = 'output/'
    OPEN(unit_number,file=trim(path_out)//'result.txt',status='unknown',action='write')
    OPEN(unit_number2,FILE=trim(path_out)//'output.txt')
    OPEN(unit_number3,FILE=trim(path_out)//'convergence.dat')
    OPEN(unit_number4,FILE=trim(path_out)//'the_comparison.dat')
    OPEN(unit_number5,FILE=trim(path_out)//'cputime.dat',status='unknown',position='append',action='write')
    
    uio = 45
    WRITE(*,*) 'this makes no sense'
    OPEN(uio,file=trim(path_out)//'mech.dat',status='unknown')
    
    CALL CPU_TIME(cputime)
    WRITE(unit_number2,*)'***************************************************************************'
    WRITE(unit_number5,*)'Start time: ', cputime
    cputime_1 = cputime
    cputime_i = cputime

    !~~~~~~~~~~~~~~~~~~~READ LATTICE INFORMATIONS~~~~~~~~~~~~~~~~~~~~~

    CALL read_force_constants

    CALL CPU_TIME(cputime)
    WRITE(unit_number5,*) 'Time takes for <read_force_constants>',cputime-cputime_1
    cputime_1 = cputime
    WRITE(*,*)'*************************Check Reading files*************************'
    WRITE(unit_number2,*)'*************************Check Reading files*************************'
    WRITE(*,*) '*************************Read lat_fc.dat*************************'
    WRITE(unit_number2,*) '*************************Read lat_fc.dat*************************'
    WRITE(*,*) 'atom number in primitive cell=',atom_number
    WRITE(unit_number2,*) 'atom number in primitive cell=',atom_number
    WRITE(*,*) 'total atom number=',tot_atom_number
    WRITE(unit_number2,*) 'total atom number=',tot_atom_number
    WRITE(*,*) 'tranlational vector=',trans_vec
    WRITE(unit_number2,*) 'tranlational vector=',trans_vec
    WRITE(*,*) 'NUMBER of PRIMITIVE CELLS:',cell_number
    WRITE(unit_number2,*) 'NUMBER of PRIMITIVE CELLS:',cell_number
    WRITE(*,*) '*************************End of lat_fc.dat check*************************'
    WRITE(unit_number2,*) '*************************End of lat_fc.dat check*************************'
    WRITE(*,*)
    WRITE(unit_number2,*)

    !~~~~~~~~~~~~~~~~~~~~~READ INDEPENDENT/NONINDEPENDENT FC MAP~~~~~~~~~~~~~~~~~~

    variational_parameters_size(3)=eff_fc2_terms

    ! resolution = threshold to check and fix asr in real fcs
    resolution = 1e-7 * max_fc2
    CALL fix_asr_fc2
    CALL fix_asr_fc3
    CALL fix_asr_fc4
    ! get the max distance between two input atoms
    ! legacy variable for checking purpose, may not be used
    R_0 = MAX(lengtha(trans_vec(:,1)),lengtha(trans_vec(:,2)),lengtha(trans_vec(:,3)))*maxneighbors


    WRITE(unit_number2,*)"*************************FORCE CONSTANTS INITIALIZATION*************************"
    WRITE(unit_number2,*)'number of atoms within fc2 cutoff region: ', SIZE(atoms_fc2)
    WRITE(unit_number2,*)'number of atoms within fc3 cutoff region: ', SIZE(atoms_fc3)
    WRITE(unit_number2,*)'number of atoms within fc4 cutoff region: ', SIZE(atoms_fc4)
    CALL CPU_TIME(cputime)
    WRITE(unit_number5,*) 'Time takes for ASR check and fix:',cputime-cputime_1
    cputime_1 = cputime

    !~~~~~~~~~~~~~~~~~READ K POINT MESH SIZE AND GENERATE~~~~~~~~~~~~~~~~~~

    !-------read params.phon-----------
    CALL read_params
    !------------setup k mesh----------
    nk=nc1*nc2*nc3
    k_number=nk
    wmin=-0.1;  nband=1
    dk=(length(g1)+length(g2)+length(g3))/(nc1+nc2+nc3)
WRITE(*,*) 'MARK a'
    CALL kvector_Update(nband,nk)!don't forget to turn this on for tetrahedron k mesh
WRITE(*,*) 'MARK b'

    !~~~~~~~~~~~~~~~~GET THE FOURIER CONSTANTS LIST~~~~~~~~~~~~~~~~~

    !--------read params.born----------
    CALL read_born
WRITE(*,*) 'MARK c'

    CALL Allocate_FFC(k_number)
WRITE(*,*) 'MARK d'

    ! CALL Fourier_Force_Constants_Calculation(k_number) !MODIFY: 03/12/2024 this is useless

    !~~~~~~~~~~~~~~~~~~~GET THE ITERATION && VARIATIONAL PARAMETERS~~~~~~~~~~~~~~~~~~~~~

    if (mpi_rank==0) then

    WRITE(*,*) '*************************Check of Iteration Parameters*************************'
    WRITE(unit_number2,*) '*************************Check of Iteration Parameters*************************'

    endif

    ! Probably OK for all process to read same file
WRITE(*,*) 'MARK0'
    CALL read_iteration_parameters_file
WRITE(*,*) 'MARK1'
    CALL initiate_var
WRITE(*,*) 'MARK2'

    if (mpi_rank==0) then
    WRITE(*,*) 'tolerance for convergence of f: ',tolerance2
    WRITE(unit_number2,*) 'tolerance for convergence of f: ',tolerance2
    WRITE(*,*) 'converted temperature: ',temperature*(100*h_plank*c_light)/k_b
    WRITE(unit_number2,*) 'converted temperature: ',temperature*(100*h_plank*c_light)/k_b
    WRITE(*,*) 'max iteration number: ',max_it_number
    WRITE(unit_number2,*) 'max iteration number: ',max_it_number
    WRITE(*,*)'Variationals type1: Atomic deviation in primitive cells; the number is: ',variational_parameters_size(1)
    WRITE(unit_number2,*)'Variationals type1: Atomic deviation in primitive cells; the number is: ',variational_parameters_size(1)
    WRITE(*,*)'Variationals type2: Strain; the number is: ',variational_parameters_size(2)
    WRITE(unit_number2,*)'Variationals type2: Strain; the number is: ',variational_parameters_size(2)
    WRITE(*,*)'Variationals type3: Trial Force Constants(RANK2); the number is: ',variational_parameters_size(3)
    WRITE(unit_number2,*)'Variationals type3: Trial Force Constants(RANK2); the number is: ',variational_parameters_size(3)
    WRITE(*,*)
    WRITE(unit_number2,*)
    WRITE(*,*)'I use the real fc2 for the 1st set of trial fc2 value'
    WRITE(unit_number2,*)'I use the real fc2 for the 1st set of trial fc2 value'
    WRITE(*,*)'*************************End of Iteration Parameters check*************************'
    WRITE(unit_number2,*)'****************End of Iteration Parameters check*************************'
    WRITE(*,*)
    WRITE(unit_number2,*)
    CALL CPU_TIME(cputime)
    WRITE(unit_number5,*)'Cputime after reading all the input:', cputime
    cputime_1 = cputime
    endif

WRITE(*,*) '!~~~~~~~~~~~~~~~~~INITIALIZE MATRIX PARAMETERS~~~~~~~~~~~~~~~~~~'
    !~~~~~~~~~~~~~~~~~INITIALIZE MATRIX PARAMETERS~~~~~~~~~~~~~~~~~~

    YY_flag=.False.
    eigen_number=d*atom_number
    ndyn = eigen_number
    CALL allocate_eigen(eigen_number,k_number)
    ! CALL Allocate_Gradients !NOTE: temporarily commented off for not enough memory

WRITE(*,*) '!~~~~~~~~~~~~INITIALIZE THE VARIATIONAL PARAMETERS~~~~~~~~~~~~~~~'
    !~~~~~~~~~~~~INITIALIZE THE VARIATIONAL PARAMETERS~~~~~~~~~~~~~~~
    IF(inherit) CALL target_update
    IF(rand_start) CALL test_update

    !~~~~~~~~~~~~~~~~~~~~~~~~~~UTILITY SECTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !---------------1. Quck Thermo Calculations---------------
    !! Turn off comment for a quick phonon band, gruneisen, specific heat,
    !! elastic constants calculation given input files + targetInitialize.dat
    ! CALL initiate_yy(kvector)
    ! CALL matrix_2nd_deriv(vals, vecs)
    ! CALL Get_dispersion
    ! CALL calc_gruneisen

    !---------------2. Quick Elastic Constants Calculations
    IF(atom_number.eq.1) THEN
        CALL GetElastic_final
        CALL calc_modulus
    ELSE
        ! CALL GetElastic_Wallace !my previous own 
        !UPDATE: FOCEX_ec

        ALLOCATE(gama(ndyn-3,ndyn-3))
        !MODIFY: need to calculate volume_r0 first
        CALL calculate_volume(r1, r2, r3, volume_r0)
        CALL get_phi_zeta_Xi(uio) !ndyn,atld0,gama,phi,zeta,teta,xi,qiu,uio)
        CALL residuals (uio) !ndyn,xi,zeta,phi,gama,sigma0,y0,pi0,uio)
        CALL mechanical2(elastic,uio) !ndyn,atld0,sigma0,phi,zeta,xi,qiu,gama,elastic,uio)
        CALL calc_modulus
        CLOSE(uio)
    END IF
    STOP
    !-------3. Quick Free Energy Landscape Calculations-------
    !! A utility subroutine for free energy landscape calculation
    !! with only selected variational parameters free to change by step 
    !! comment off the specific block which you want to test with
    step = 0.005
    ! CALL initiate_yy(kvector)

    !NOTE: <small_test2(step)> is a general validation test for Broyden method
    ! it compares the difference of all x(:) and f(:) from finite difference method
    ! CALL small_test2(step) 

    !NOTE: <small_test3(i, step, n)> is a specific validation test for Broyden method
    ! it compares the difference of f(i) from finite difference method with x(i) moving n steps.
    ! However, the x(i) can only be u0 or eta; not K.
    ! CALL small_test3(15,step,4)
    
    !NOTE: <small_test3_yy(i, step, n)> does the same thing as above, this time for <yy>
    ! CALL small_test3_yy(9,step,4) 

    !NOTE: <small_test(i, step, n)> calculates the free energy F and f(i)
    ! given selected variational variable x(i) move n steps
    ! CALL small_test(6,0.01d0,30) 

    !NOTE: <small_test_ex(i, j, step, n)> calculates the free energy F and f(i)
    ! given 2 selected variables x(i) and x(j) ranged in (-n*step, n*step]
    ! start_i, start_j can be designated manually in the subroutine
    ! deprecated
    ! CALL small_test_ex(4,8,step,10)

    !NOTE: <rhom_contour(i,j,step1, step2, n)> does the same thing as above
    ! with x(i) ranged in (-n*step1 + start_i, n*step1 + start_i]
    ! x(j) ranged in (-n*step2 + start_j, n*step2 + start_j]
    ! start_i, start_j can be designated in targetInitialize.dat, so turn on inherit option
    ! or assigned manually in the subroutine
    ! used for contour plot
    step_1 = 0.0001; step_2 = 0.0001
    CALL rhom_contour(4,8,step_1,step_2,40)
    STOP

    !~~~~~~~~~~MAKE AN INITIAL GUESS BASED ON FC2 DIAGONALIZATION~~~~~~~~~~~~~~~~
    !NOTE: might not used
    accept = -1d-10
    min_eival = accept

    !~~~~~~~~~~~~~~~~~~~~~~~FIRST RUN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (mpi_rank==0) then

        WRITE(*,*)'**********************Check 1st Variational Approach calculation**********************'
        WRITE(unit_number2,*)'*********************Check 1st Variational Approach calculation*********************'
        CALL CPU_TIME(cputime)
        WRITE(unit_number5,*)'Cputime before 1st Variational Approach run:', cputime
        cputime_1 = cputime
    endif

    iter_rec = 0
    CALL asr_checkfc2(trialfc2_value)!this is just to check if input fc2/inherited trialfc2 satisfies asr
    
    CALL GetF_trial_And_GradientF_trial(kvector)

    start_F = REAL(F_trial)

    if (mpi_rank==0) then

        CALL CPU_TIME(cputime)
        WRITE(unit_number5,*)'Time takes for the 1st variational approach run', cputime-cputime_1
        cputime_1 = cputime

        WRITE(*,*)'start_F', start_F
        WRITE(*,*) 'first Unpertubed Free Energy F0=',F0
        WRITE(unit_number2,*) 'first Unpertubed Free Energy F0=',F0
        WRITE(unit_number,*) 'first Unpertubed Free Energy F0=',F0
        WRITE(*,*) 'first Corresponding Potential Energy V0=',V0
        WRITE(unit_number2,*) 'first Corresponding Potential Energy V0=',V0
        WRITE(*,*) 'first Trial Free Energy=F0+V_avg-V0=',F_trial
        WRITE(unit_number2,*) 'first Trial Free Energy=F0+V_avg-V0=',F_trial
        WRITE(*,*)'*********************End of Variational Approach calculation Check***********************'
        WRITE(unit_number2,*)'******************End of Variational Approach calculation Check********************'
        WRITE(*,*)
        WRITE(unit_number2,*)
        WRITE(unit_number2,*)'ITERATION#:',0
        WRITE(unit_number2,*) &
        &'*************************Variational Parameters Check After a Broyden iteration*************************'
        WRITE(unit_number2,*) '||Atomic deviation u_tau(:)||'
        DO j=1,atom_number
            WRITE(unit_number2,*)'atom: ', j
            WRITE(unit_number2,*) atomic_deviation(:,j)
        END DO
        WRITE(unit_number2,*)
        WRITE(unit_number2,*)'||Strain Tensor||'
        DO j=1, d
            WRITE(unit_number2,*) strain(j,:)
        END DO
        WRITE(unit_number2,*)
        WRITE(unit_number2,*)'||Force Constants||'
        DO j=1,eff_fc2_terms
            WRITE(unit_number2,*)'---------------------------------'
            WRITE(unit_number2,*) 'effective fc2 number: ', j
            t_atom1 = eff_fc2_index(j)%iatom_number
            t_atom2 = eff_fc2_index(j)%jatom_number
            t_xyz1 = eff_fc2_index(j)%iatom_xyz
            t_xyz2 = eff_fc2_index(j)%jatom_xyz
            WRITE(unit_number2,*) '1st atomic index:', t_atom1
            WRITE(unit_number2,*) '1st directional index:', t_xyz1
            WRITE(unit_number2,*) '2nd atomic index:', t_atom2
            WRITE(unit_number2,*) '2nd 2directional index:', t_xyz2
            WRITE(unit_number2,*) 'effective fc2 value:', trialfc2_value(t_atom1, t_atom2)%phi(t_xyz1, t_xyz2)
            WRITE(unit_number2,*) 'corresponding (input file) fc2 value:',myfc2_value(t_atom1, t_atom2)%phi(t_xyz1, t_xyz2)
            WRITE(unit_number2,*)
        END DO

        WRITE(*,*)'****************************************************************************************************'
        WRITE(*,*)'*************************Broyden/CG Iterations*************************'
        WRITE(unit_number2,*)'*********************Broyden/CG Iterations*********************'
        WRITE(unit_number3,*)'*****iteration #, L1 norm of all gradients, free energy value*****'
        WRITE(unit_number3,*) 'iteration 0 free energy ', REAL(F_trial)

    endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~SETUP BROYDEN/CG ITERATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    iter = 0; ft=100; ft2=100; threshold=100; itmax=max_it_number
    pmix = my_pmix
    err = 1e-9
    mx=SIZE(GradientF_trial)
    IF(ALLOCATED(x)) DEALLOCATE(x)
    IF(ALLOCATED(f)) DEALLOCATE(f)
    IF(ALLOCATED(g)) DEALLOCATE(g)
    IF(ALLOCATED(h)) DEALLOCATE(h)
    ALLOCATE(x(mx),f(mx))
    ALLOCATE(g(mx),h(mx))
    !if there is any constrained parameter
    IF(.NOT.unleashed) THEN
        item = SUM(variational_parameters_size) - SIZE(fixed_params)
        IF(ALLOCATED(x2)) DEALLOCATE(x2)
        IF(ALLOCATED(f2)) DEALLOCATE(f2)
        ALLOCATE(x2(item),f2(item))
        mx=item
    END IF
    CALL allocatebroy(mx,itmax)
    fp=F_trial;f=GradientF_trial;g=-f;h=g !initialize

    !~~~~~~~~~~~~~~~~~~~~ENTER Broyden/CG LOOP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !GRAY: entry point of Broyden/CG iterations
    DO while ( iter<itmax .and. threshold>err )

          iter = iter+1
          iter_rec = iter

        if (mpi_rank==0) then
            WRITE(*,*)'ITERATION#:',iter
            WRITE(unit_number2,*)'ITERATION#:',iter
            CALL CPU_TIME(cputime)
            cputime_1 = cputime
        endif
    !---------------update old x(:) for next loop-----------------
        CALL combine_variants(x)  

        PRINT*,'process ',mpi_rank,' f(:),x(:) constructed, iter=',iter

        IF(unleashed) THEN
            IF(method.eq.0) THEN
                CALL SimpleUpdate(mx,x,f)
            ELSE IF(method.eq.1) THEN
                CALL bro90(mx,x,f,iter)
            ELSE IF(method.eq.2) THEN
                f=h
                CALL linmin(x,f,mx,fret)
                IF(2.*abs(fret-fp).le.err*(abs(fret)+abs(fp)+EPS)) THEN
                    WRITE(unit_number2,*)'!!!!! Terminated because of condition 1 !!!!!'
                    CALL decompose_variants(x) !assign x(:) back since terminated early
                    EXIT
                END IF
                fp=fret
            END IF
        ELSE
            CALL select_xy(fixed_params,x,x2)
            CALL select_xy(fixed_params,f,f2)
            IF(method.eq.0) THEN
                CALL SimpleUpdate(mx,x2,f2)
            ELSE IF(method.eq.1) THEN
                CALL bro90(mx,x2,f2,iter)
            ELSE IF(method.eq.2) THEN
                f=h
                CALL select_xy(fixed_params,f,f2)
                CALL linmin(x2,f2,mx,fret)
                IF(2.*abs(fret-fp).le.err*(abs(fret)+abs(fp)+EPS)) THEN
                    WRITE(unit_number2,*)'!!!!! Terminated because of condition 1 !!!!!'
                    CALL decompose_variants(x) !assign x(:) back since terminated early
                EXIT
                END IF
                fp=fret
            END IF
            CALL release_x(fixed_params,x2,x)
        END IF

    !---------------Calculate Gradients with new x(:) for next loop-----------------

        CALL decompose_variants(x)
        CALL GetF_trial_And_GradientF_trial(kvector) 
        f(:)=GradientF_trial(:) 

        !NOTE: check and fix ASR in trial fc2, theoretically ASR should be automatically preserved
        CALL asr_checkfc2(trialfc2_value) 

        IF(method.eq.2) THEN
            CALL update_fgh(SIZE(GradientF_trial),f,g,h)
            !NOTE: after this, f(:) no longer equal GradientF_trial(:)
        END IF

    !----------------Calculate L1 norm of gradients---------------------------------
        ft = 0d0; ft2 = 0d0
        IF(unleashed) THEN
            !threshold for all variables
            DO i=1,variational_parameters_size(1)+variational_parameters_size(2)
                ft = ft + ABS(f(i)*x(i))
            END DO
            DO i=variational_parameters_size(1)+variational_parameters_size(2)+1, SIZE(x)
                j = i - variational_parameters_size(1) - variational_parameters_size(2)
                ft = ft + ABS(f(i)*trialfc2_record(j))
            END DO
            ft = ft/mx/(temperature*100*h_plank*c_light*6.242d+18) !dimensionless
            threshold = ft
            WRITE(*,*)' norm of f = ',ft
            WRITE(unit_number2,*)' norm of f = ',ft
        ELSE
            !threshold for free variables only
            ft2=normGradients(free_var,x2,f2)
            ft2=ft2/SIZE(free_var)/(temperature*100*h_plank*c_light*6.242d+18)!dimensionless
            threshold = ft2
            WRITE(*,*)' norm of free f = ',ft2
            WRITE(unit_number2,*)' norm of free f = ',ft2
        END IF

    !~~~~~~~~~~~~~~~~Variational Parameters Check After a Broyden/CG iteration~~~~~~~~~~~~~~~~~~~
        if(mpi_rank==0) then
            WRITE(unit_number2,*) &
                &'********************Variational Parameters Check After a Broyden iteration*********************'
            
            WRITE(*,*)'***************************************************************************'
            WRITE(*,*) 'Translational Vector now is: ',trans_vec
            WRITE(unit_number2,*)'Translational Vector now is: ',trans_vec
            WRITE(*,*) 'Free Energy F=F0+<V-V0>=:', F_trial
            WRITE(unit_number2,*) 'Free Energy F=F0+<V-V0>=:', F_trial
            WRITE(*,*)
            WRITE(unit_number2,*)
            CALL CPU_TIME(cputime)
            WRITE(unit_number2,*) 'cputime at the',iter,'th iteration', cputime
            WRITE(unit_number3,*) iter,',',threshold,',',REAL(F_trial)
            CALL CPU_TIME(cputime)
            WRITE(unit_number5,*) 'Time takes for the ',iter,'th iteration:',cputime-cputime_1
        endif

    END DO !GRAY: the do while loop for the whole Broyden

    !~~~~~~~~~~~~~~~~~~~~~~~EXIT from Broyden/CG LOOP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    WRITE(*,*)'*************************End of Broyden/CG Iterations*************************'
    WRITE(unit_number2,*)'********************End of Broyden/CG Iterations*********************'
    WRITE(unit_number2,*)''
    WRITE(unit_number2,*)''

    CALL printFinalGradients(x) !final gradients value and variational parameters value in a separate file
    CALL printResults

!-------------------------------------------------------------------------------------
    if(mpi_rank==0) then
        WRITE(*,*)'final translational vector =',trans_vec
        WRITE(unit_number,*)'final translational vector =',trans_vec
        WRITE(*,*) 'final F0 = ', F0
        WRITE(unit_number,*) 'final F0 = ', F0
        WRITE(*,*) 'final Free Energy F=F0+<V-V0> = ',F_trial
        WRITE(unit_number,*) 'final Free Energy F=F0+<V-V0> = ',F_trial
        WRITE(*,*)'Temperature',temperature*(100*h_plank*c_light)/k_b
        WRITE(unit_number,*)'Temperature',temperature*(100*h_plank*c_light)/k_b
        CALL asr_checkfc2(trialfc2_value) !check final ASR

        CALL CPU_TIME(cputime)
        WRITE(unit_number5,*)'Total running time: ', cputime-cputime_i
    endif
!-------------------------------------------------------------------------------------

    !~~~~~~~~~~~~~~~~~~~~~~~POST PROCESS && CHECK~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !------------thermal expansion test---------------
    CALL R_Update
    CALL calculate_volume(r1,r2,r3,volume)
    WRITE(unit_number2,*)'Current Translational Vector 1: ',trans_vec(1,:)
    WRITE(unit_number2,*)'Current Translational Vector 2: ',trans_vec(2,:)
    WRITE(unit_number2,*)'Current Translational Vector 3: ',trans_vec(3,:)
    WRITE(unit_number2,*)'Current Volume is:',volume
    WRITE(unit_number,*)'Current Volume is:',volume

    !------------final phonon dispersion-------------------
    CALL initiate_yy(kvector) !Get new eivals+deals negative eigenvalues
    CALL GetF0_and_V0 !use new eivals for new F0 and V0

!    CALL updateK !turn off this when trialfc2 is also set free
    CALL Get_dispersion! to output dispersion data for plotting, frequencies in cm^-1

    !-----add self energy?-----
    ksubset_F = 665
    OPEN(ksubset_F,FILE='ksubset.inp')
    WRITE(ksubset_F,'(2I5)') 1, nibz
    CLOSE(ksubset_F)
!    CALL read_ksubset(ksubset,nibz)
!    CALL CPU_TIME(cputime_1)
!    ! regular self energy on coarse kpc(:,:) mesh
!    CALL mySelfEnergy_re(ksubset,temperature,'pre_35.dat')
!    CALL CPU_TIME(cputime_2)
!    WRITE(*,*) 'Time cost of self-energy part = ',cputime_2-cputime_1
!CALL Get_Dispersion_SE

    !-----final phonon DOS------ !temporarily commented off for check
!    CALL calc_dos_gauss
!    CALL calc_dos_tet

    !-----final elastic constants-----
    IF(atom_number.eq.1) THEN
        CALL GetElastic_final
    ELSE
        ! CALL GetElastic_Wallace !my previous own 
        !UPDATE: FOCEX_ec
        uio = 345
        OPEN(uio,file=trim(path_out)//'mech.dat')
        ALLOCATE(gama(ndyn-3,ndyn-3))
        !MODIFY: need to calculate volume_r0 first
        CALL calculate_volume(r1, r2, r3, volume_r0)
        CALL get_phi_zeta_Xi(uio) !ndyn,atld0,gama,phi,zeta,teta,xi,qiu,uio)

        CALL residuals (uio) !ndyn,xi,zeta,phi,gama,sigma0,y0,pi0,uio)
        CALL mechanical2(elastic,uio) !ndyn,atld0,sigma0,phi,zeta,xi,qiu,gama,elastic,uio)
        CLOSE(uio)
    END IF
    !-----final gruneisen------ 
    CALL calc_gruneisen

    

    WRITE(unit_number2,*)'Start Free Energy = ', start_F
    WRITE(unit_number2,*)'Final Free Energy = ', REAL(F_trial)
    WRITE(*,*)'...Finish'

CLOSE(unit_number)
CLOSE(unit_number2)
CLOSE(unit_number3)
CLOSE(unit_number4)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEALLOCATE EVERYTHING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL deallocate_tetra
    DEALLOCATE(x,f)
    DEALLOCATE(g,h)
    IF(.not.unleashed) DEALLOCATE(x2,f2)
!   DEALLOCATE(junk,kpc,wk,eigen_temp ,total_dos,integrate_dos,afunc_dos)


8 FORMAT(99(2x,g14.8))
9 format(a,99(1x,g12.6))
10 FORMAT(2(I6.1),2(I6.1,I2.1),4x,F9.6)
11 FORMAT(2(I2,A),3(A,F9.6))
!15 format(99(2x,g14.8))

END SUBROUTINE ThreeVariableRoutine
!===============================================================================================
SUBROUTINE dynamicSelect(limit)
    USE CG
    IMPLICIT NONE

    LOGICAL, DIMENSION(:),ALLOCATABLE :: dynamicFlag
    REAL(8), INTENT(in) :: limit
    INTEGER :: i, counter, j

    ALLOCATE(dynamicFlag(SIZE(GradientF_trial)))
    counter = 0
    DO i=1, SIZE(GradientF_trial)
        IF(ABS(GradientF_trial(i)).le.limit) THEN
            dynamicFlag(i) = .True.
            counter = counter + 1
        ELSE
            dynamicFlag(i) = .False.
        END IF
    END DO

    !reassign which parameters to fix
    IF(ALLOCATED(fixed_params)) DEALLOCATE(fixed_params)
    ALLOCATE(fixed_params(counter))
    j=1
    DO i=1, SIZE(GradientF_trial)
        IF(dynamicFlag(i)) THEN
            fixed_params(j) = i
            j = j + 1
        END IF
    END DO

    DEALLOCATE(dynamicFlag)
    WRITE(*,*)'Total Free Variational Parameter Number this iteration: ', &
                    &SIZE(GradientF_trial) - SIZE(fixed_params)
END SUBROUTINE dynamicSelect
!===============================================================================================