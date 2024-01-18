MODULE Iteration_parameters
    !!This module is for reading iteration parameters from various config files
    USE DFT_force_constants

    IMPLICIT NONE

    INTEGER :: iter_rec
    INTEGER :: seed
    INTEGER :: max_it_number,indie_fc2
    REAL(8) :: temperature,tolerance2,danger
    REAL(8) :: my_pmix !for user input
    REAL(8) :: rand_range, rand_center(3) !for user input
    INTEGER,DIMENSION(:),ALLOCATABLE :: fixed_params
    REAL(8),ALLOCATABLE :: atomic_deviation(:,:),strain(:,:) !atomic_deviation(xyz,atom_type),strain(xyz,xyz)
    REAL(8),ALLOCATABLE :: atomic_deviation_sav(:,:),strain_sav(:,:)
    TYPE(fc2_value),DIMENSION(:,:),ALLOCATABLE :: trialfc2_value,prev_trialfc2_value,trialfc2_value_initial
    TYPE(fc2_value),DIMENSION(:,:),ALLOCATABLE :: yy_value,yy_value_initial
    INTEGER,DIMENSION(:),ALLOCATABLE :: free_var
    LOGICAL :: unleashed = .false.
    LOGICAL :: inherit = .false.
    LOGICAL :: rand_start = .false.
    LOGICAL :: pressure = .True.
    LOGICAL :: forceK_asr = .False.
    LOGICAL :: highT_limit = .False.
    REAL(8) :: stress(d,d) !for pressure test
    CHARACTER path_out*29  !UPDATE: for output path

CONTAINS
!--------------------------------------------------------------------------------------------
 subroutine read_params
     !! Read phonon parameters from <params.phon>
 use io2
 use om_dos
 use phi3
 use kpoints
 use atoms_force_constants
 implicit none
! integer i

  open(uparams,file='params.phon',status='old')
  write(6,*) 'READ_PARAMS: opening params.phon'

  read(uparams,*)nc1,nc2,nc3,dos_bs_only ! coarse kmesh for big calculations
  NC(1)=nc1 ;NC(2)=nc2 ;NC(3)=nc3
  !write(ulog,*) 'READ_PARAMS: just read nib ',nc1
  read(uparams,*)lamin,lamax
  read(uparams,*)shftx,shfty,shftz ! kmesh for dos calculation
  read(uparams,*)wmesh,wmax  ! no of w mesh points and max frequency for DOS
  !write(ulog,*) 'READ_PARAMS: just read wmax ',wmax
  read(uparams,*)width,etaz  ! width of gaussian broadening for DOS,imag part of om
  !write(ulog,*) 'READ_PARAMS: just read width,etaz ',width,etaz
  read(uparams,*)tau0        ! large relaxation time: 1/tau0 to be added to gamma
  !write(ulog,*) 'READ_PARAMS: just read tau0 ',tau0
  read(uparams,*)verbose
  read(uparams,*)wshift      ! small shift added to make negative modes go away
  !write(ulog,*) 'READ_PARAMS: just read wshift ',wshift
  read(uparams,*)tmin,tmax,ntemp      ! temps in Kelvin to calculate thermal expansion and other thermal ppties
  read(uparams,*)iter,readv3,usetetra,isvd      ! if=1 then read from a file, else generate
  read(uparams,*)ncpu
  read(uparams,'(a)')v3path
  !write(ulog,*) 'READ_PARAMS: path is:',v3path
  read(uparams,*)max_iter, n_dig_acc
!, conv_error, conv_max_error, conv_diff, conv_max_diff, &
!  & conv_diff_kap, conv_max_diff_kap,conv_iter,update
   ! maximum iteration number for iterative solution, convergence criteria, update ratio
  read(uparams,*)v3_threshold        ! in 1/cm keep all v33 of norm above v3_threshold
  read(uparams,*)classical     ! if =1 then use classical distribution (kT/hw)
  read(uparams,*)cal_cross,qcros     ! if =1 then calculate the cross section at qcros (reduced U)
  read(uparams,*)threemtrx          ! if =1 then calculate the 3-phonon matrix elements
  !read(uparams,*)scalelengths       ! multiply all lattice and atomic coordinates by scalelength

  close(uparams)

  !write(ulog,*) 'READ_PARAMS: read and kpoints allocated'
  !write(ulog,3)'nc1,nc2,nc3=',nc1,nc2,nc3
  !write(ulog,*)'wmesh, wmax=',wmesh,wmax
  !write(ulog,*)'readv3=',readv3  !,writev3

!  write(ulog,*)'Tmin,Tmax(K=',tmin,tmax
!  if (classical .eq. 1) then
!     write(ulog,*)'Calculation is done in the CLASSICAL limit'
!  else
!     write(ulog,*)'Calculation is done in the QUANTUM limit'
!  endif
!  if (dos_bs_only) then
!     write(ulog,*)' DOS will only be calculated (on the coarse mesh) ',nc1,nc2,nc3
!  endif

3 format(a,6(1x,i6))
 end subroutine read_params
!----------------------------------------------------------------------------------------------------
    SUBROUTINE initiate_var
         !! initiate variational parameters to trivial start
         IMPLICIT NONE
         INTEGER :: tau1,direction1,direction2
         INTEGER :: seed

         IF(ALLOCATED(atomic_deviation)) DEALLOCATE(atomic_deviation)
         IF(ALLOCATED(strain)) DEALLOCATE(strain)
         IF(ALLOCATED(trialfc2_value)) DEALLOCATE(trialfc2_value)
         IF(ALLOCATED(prev_trialfc2_value)) DEALLOCATE(prev_trialfc2_value)
         IF(ALLOCATED(yy_value)) DEALLOCATE(yy_value)
         IF(ALLOCATED(strain_sav)) DEALLOCATE(strain_sav)
         IF(ALLOCATED(atomic_deviation_sav)) DEALLOCATE(atomic_deviation_sav)

         ALLOCATE(atomic_deviation(d,atom_number),atomic_deviation_sav(d,atom_number))
         ALLOCATE(strain(d,d),strain_sav(d,d))
         ALLOCATE(trialfc2_value(atom_number,tot_atom_number))
         ALLOCATE(prev_trialfc2_value(atom_number,tot_atom_number))
         ALLOCATE(yy_value(atom_number,tot_atom_number))

         atomic_deviation=0d0;atomic_deviation_sav=0d0

         strain=0d0;strain_sav=0d0

         trialfc2_value=myfc2_value
         prev_trialfc2_value=trialfc2_value

        !**** allocate trialfc2_value_initial to keep the first trialfc2 from input data ****
        !note: this is for if user wants to 'freeze' some selected fc2, which means
        !...corresponding trialfc2_value are frozen to its trialfc2_value_initial
        IF(ALLOCATED(trialfc2_value_initial)) DEALLOCATE(trialfc2_value_initial)
        ALLOCATE(trialfc2_value_initial(atom_number,tot_atom_number))
        trialfc2_value_initial=trialfc2_value

    END SUBROUTINE initiate_var
!----------------------------------------------------------------------------------------------------
    SUBROUTINE make_rhombohedral
        !!just for Bismuth 385K test, manually make it a rhombohedral
        IMPLICIT NONE
        REAL(8) :: x,y,utau

        x=0d0
        y=0d0
        utau = 0d0


        strain(:,1) = (/x,y,y/)
        strain(:,2) = (/y,x,y/)
        strain(:,3) = (/y,y,x/)

        atomic_deviation(:,2) = (/utau,utau,utau/)

    END SUBROUTINE make_rhombohedral
!----------------------------------------------------------------------------------------------------
    SUBROUTINE test_update
        !!Randomize eta and utau for random start, optional
        !!run initiate_var first
        IMPLICIT NONE
        INTEGER :: i,j,k,l
        INTEGER ::atom1,atom2,xyz1,xyz2
        REAL(8) :: x(3)
        !randomize strain eta
        CALL srand(seed)
        strain(:,1) = strain(:,1) + rand_range*(/rand()-0.5,rand()-0.5,rand()-0.5/)+rand_center
        strain(:,2) = strain(:,2) + rand_range*(/rand()-0.5,rand()-0.5,rand()-0.5/)+rand_center
        strain(:,3) = strain(:,3) + rand_range*(/rand()-0.5,rand()-0.5,rand()-0.5/)+rand_center

        !randomize atomic deviation utau
        CALL RANDOM_NUMBER(x)
        atomic_deviation(:,1) = atomic_deviation(:,1) + (2*x-1)*0.05
        CALL RANDOM_NUMBER(x)
        atomic_deviation(:,2) = atomic_deviation(:,2) + (2*x-1)*0.05

        !randomize effective fc2, not needed


    END SUBROUTINE test_update
!----------------------------------------------------------------------------------------------------
    SUBROUTINE target_update
    !!Initialize variational parameters from 'targetInitialize.dat' file, optional
    !!inherited should be 'true'
    !!run initiate_var first
        IMPLICIT NONE
        INTEGER :: i,idx
        INTEGER :: atom1,atom2,xyz1,xyz2
        REAL(8) :: temp
        
        !NOTE: not path output
        OPEN(21,FILE='targetInitialize.dat',STATUS='old',ACTION='read')

        !initialize utau, one vector/atom per line
        DO i=1,atom_number
            READ(21,*) atomic_deviation(:,i)
        END DO
        !initialize strain, 3x3
        DO i=1,3
            READ(21,*) strain(:,i)
        END DO
        !initialize trailfc2
        DO i=1,SIZE(myfc2_index)
            IF(myfc2_index(i)%iatom_number.gt.atom_number) CYCLE
            READ(21,*) idx,atom1,atom2,xyz1,xyz2,temp
            trialfc2_value(atom1,atom2)%phi(xyz1,xyz2) = temp
        END DO
        CLOSE(21)

    END SUBROUTINE
!----------------------------------------------------------------------------------------------------
    SUBROUTINE select_xy(choice,x_in,x_out)
    !!utility subroutine for freeze/free variational params
    !!choice is fixed_params
    !!select x feed into this iteration, used after <combine_x> before <bro90>
        IMPLICIT NONE
        INTEGER :: i,j,idx
        INTEGER,INTENT(in),DIMENSION(:) :: choice
        REAL(8),INTENT(in),DIMENSION(:) :: x_in
        REAL(8),INTENT(out), DIMENSION(:), ALLOCATABLE :: x_out

        IF(ALLOCATED(x_out)) DEALLOCATE(x_out)
        ALLOCATE(x_out(SUM(variational_parameters_size)-SIZE(choice)))
        j=1;idx=1
        DO i=1,SIZE(x_in)
            IF(i.ne.choice(j)) THEN
                x_out(idx) = x_in(i)
                idx = idx + 1
            ELSE
                j = j + 1
            END IF
        END DO
    END SUBROUTINE select_xy

    SUBROUTINE release_x(choice,x_in,x_out)
    !!utility subroutine for freeze/free variational params
    !!reassign x after this iteration, used after <bro90>/<cg> before <decompose_x>
        IMPLICIT NONE
        INTEGER :: i,j,idx
        INTEGER :: temp,counter
        INTEGER :: atom1,atom2,xyz1,xyz2
        INTEGER,INTENT(in),DIMENSION(:) :: choice
        INTEGER,ALLOCATABLE,DIMENSION(:) :: inv_choice
        REAL(8),INTENT(in),DIMENSION(:) :: x_in
        REAL(8),INTENT(out), DIMENSION(:), ALLOCATABLE :: x_out

        !expand x to cover all the variational variants
        IF(ALLOCATED(x_out)) DEALLOCATE(x_out)
        ALLOCATE(x_out(SUM(variational_parameters_size)))
        x_out = 0!default is 0, but need to also maintain the old fixed values if they were not fixed to 0

        !first preserve old fixed atomic deviation and strain tensor values
        DO i=1,variational_parameters_size(1)
            IF(MOD(i,d).eq.0) THEN
                x_out(i)=atomic_deviation(d,INT(i/d))
            ELSE
                x_out(i)=atomic_deviation(MOD(i,d),INT(i/d+1))
            END IF
        END DO

        DO i=variational_parameters_size(1)+1,variational_parameters_size(1)+variational_parameters_size(2)
            temp=i-variational_parameters_size(1)
            IF(MOD(temp,d).eq.0) THEN
                x_out(i)=strain(INT(temp/d),d)
            ELSE
                x_out(i)=strain(INT(temp/d+1),MOD(temp,d))
            END IF
        END DO

        !second, get free variable idx,then overwrite corresponding x_out(:) with x_in(:) from Broyden
        !i for x_out(:), j for choice(:), idx for x_in(:)
        ALLOCATE(inv_choice(SUM(variational_parameters_size)-SIZE(choice)))
        j=1;idx=1;
        DO i=1,SUM(variational_parameters_size)
            IF(i.ne.choice(j)) THEN
                x_out(i) = x_in(idx)
                inv_choice(idx) = i
                idx = idx + 1
            ELSE
                j = j + 1
            END IF
        END DO

        !third, just for the special case: if there is only one free diagonal strain element
        counter = 0
        temp = variational_parameters_size(1) + variational_parameters_size(2)
        DO i=1,SIZE(inv_choice)
            IF(inv_choice(i).gt.variational_parameters_size(1) .AND. inv_choice(i).le.temp) THEN
                counter = counter + 1
                j = i
            END IF
        END DO
        !make the 3 diagonal strain element equal
        IF(counter.eq.1) THEN
            temp = variational_parameters_size(1)
            x_out(temp+1) = x_in(j)
            x_out(temp+5) = x_in(j)
            x_out(temp+9) = x_in(j)
        END IF

        !fourth, preserve fixed fc2, if there is any: counter != 0
        counter = 0
        temp = variational_parameters_size(1) + variational_parameters_size(2)
        DO i=1,SIZE(choice)
            IF(choice(i).gt.temp) THEN
                counter = 1
                j = i
                EXIT
            END IF
        END DO

        IF(counter.ne.0) THEN
            DO i=j,SIZE(choice)
                temp = choice(i)- variational_parameters_size(1)-variational_parameters_size(2)

                atom1 = eff_fc2_index(temp)%iatom_number
                atom2 = eff_fc2_index(temp)%jatom_number
                xyz1 = eff_fc2_index(temp)%iatom_xyz
                xyz2 = eff_fc2_index(temp)%jatom_xyz

                x_out(choice(i))=trialfc2_value_initial(atom1,atom2)%phi(xyz1,xyz2)
                !which means for the chosen fc2, we are 'freezing' it to always be its input value
            END DO
        END IF

        IF(ALLOCATED(free_var)) DEALLOCATE(free_var)
        ALLOCATE(free_var(SIZE(inv_choice)),source=inv_choice)
        DEALLOCATE(inv_choice)

    END SUBROUTINE release_x
!----------------------------------------------------------------------------------------------------
    SUBROUTINE read_iteration_parameters_file
    !!read iteration parameters from <iteration_parameters.in>
    !!this should be called with variational_parameters_size(3) known

        IMPLICIT NONE

        INTEGER :: unit_number,n,i,j
        INTEGER,DIMENSION(:),ALLOCATABLE :: free_params

        !UPDATE: for path output initialize
        path_out = 'output/'

        unit_number = 60
!        OPEN(unit_number,file='iteration_parameters.in',status='old',action='read')
        OPEN(unit_number,file='control.params',status='old',action='read') !modified 11/10/2023

        READ(unit_number,*) rand_start !start randomly or not
        READ(unit_number,*) inherit !use results from last run to target initialize or not
        READ(unit_number,*) pressure !whether include pressure factor
        READ(unit_number,*) forceK_asr !whether force ASR in K after each Broyden
        READ(unit_number,*) highT_limit !whether use high temperature limit in <yy> calculation
!        READ(unit_number,*) tolerance2 !modified 11/10/2023
        tolerance2 = 0.0001
        READ(unit_number,*) seed !seed to generate random number
        READ(unit_number,*) rand_range !range to generate random number
        READ(unit_number,*) rand_center !center to generate random number 
        READ(unit_number,*) max_it_number
        READ(unit_number,*) my_pmix
        READ(unit_number,*) temperature

        DO i=1,3
            READ(unit_number,*) stress(i,:) !read stress tensor
        END DO
        stress = stress*1d-21/ee !unify the unit, file is in GPa

!        READ(unit_number,'(E5.0)') danger !threshold for negative eigenvalues handling, modified 11/10/2023
        danger = 1d-4
        READ(unit_number,*) n !number of fixed variational parameter
        IF(ALLOCATED(fixed_params)) DEALLOCATE(fixed_params)
        !if n variational parameters are fixed and listed
        IF(n.gt.0) THEN
            ALLOCATE(fixed_params(n))
            READ(unit_number,*) (fixed_params(i),i=1,n)
        !if all are free
        ELSEIF(n.eq.0) THEN
            unleashed = .true.
            ALLOCATE(fixed_params(1)) !this doesn't matter here, once unleashed flag sets .true.
            fixed_params = 0
        !if -n variational parameters are free and listed
        ELSE
            ALLOCATE(free_params(-n))
            READ(unit_number,*) (free_params(i),i=1,-n)
            ALLOCATE(fixed_params(SUM(variational_parameters_size)-SIZE(free_params)))

            j=1
            DO i=1,SUM(variational_parameters_size)
                IF(.NOT.ANY(free_params==i)) THEN
                    fixed_params(j) = i
                    j = j + 1
                END IF
            END DO
        ENDIF
        CLOSE(unit_number)

        temperature=temperature*k_b/100/h_plank/c_light!convert T(K) to T(1/cm)
    END SUBROUTINE read_iteration_parameters_file
!----------------------------------------------------------------------------------------------------
    SUBROUTINE drop_fc2terms
        !!get the number of independent fc2
        IMPLICIT NONE
        INTEGER :: i,temp,atom1,atom2
        temp=0
        DO i=1,fc_terms(2)
            atom1=myfc2_index(i)%iatom_number
            IF(atom1.gt.atom_number) cycle !only use the center cell atoms
            atom2=myfc2_index(i)%jatom_number
            IF(atom1.eq.atom2) cycle !use ASR to drop one term
            temp=temp+1
        END DO
        !indie_fc2=temp
        !variational_parameters_size(3) = temp !fixed the center cell
    END SUBROUTINE drop_fc2terms
!----------------------------------------------------------------------------------------------------
    SUBROUTINE asr_checkfc2(dummyFC)
        !!check asr on selected [myfc2_value] object, i.e. Phi or K
        IMPLICIT NONE
        TYPE(fc2_value),INTENT(INOUT),DIMENSION(:,:) :: dummyFC
        TYPE(fc2_value),DIMENSION(:),ALLOCATABLE :: checkfc2
        INTEGER :: atom1,atom2,direction1,direction2
        INTEGER :: i,j

        OPEN(55,FILE='asr_check_fc2.dat',STATUS='unknown',POSITION='append',ACTION='write')

        WRITE(55,*)'CHECKING ASR IN SELECTED FC2 IN ITERATION #', iter_rec
        ALLOCATE(checkfc2(atom_number))
        DO i=1,atom_number
            checkfc2(i)%phi=0
        END DO

        DO i=1,atom_number
        DO direction1=1,d
        DO direction2=1,d

        DO j=1,tot_atom_number
            checkfc2(i)%phi(direction1,direction2) = checkfc2(i)%phi(direction1,direction2)&
            &+dummyFC(i,j)%phi(direction1,direction2)
        END DO

        IF(ABS(checkfc2(i)%phi(direction1,direction2)).gt.1d-8) THEN
            WRITE(55,*)"found corresponding fc2 not satisfy ASR:"
            WRITE(55,*) checkfc2(i)%phi(direction1,direction2)
            WRITE(55,*) i,get_letter(direction1),get_letter(direction2)

            IF(forceK_asr) THEN
                dummyFC(i,i)%phi(direction1,direction2) = &
                & dummyFC(i,i)%phi(direction1,direction2) - checkfc2(i)%phi(direction1,direction2)
            END IF

        END IF

        END DO !direction2 loop
        END DO !direction1 loop
        END DO !i loop
        WRITE(55,*)'======================================='
        DEALLOCATE(checkfc2)
        CLOSE(55)
    END SUBROUTINE asr_checkfc2
!-------------------------------------------------------------------------------------------------------
END MODULE Iteration_parameters
