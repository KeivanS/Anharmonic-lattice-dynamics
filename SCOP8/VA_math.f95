Module VA_math
!!Module VA_math to get trial Free Energy and its gradients
!!it contains most of the math subroutines
    USE Iteration_parameters
    USE Fourier_force_constants
    USE MatrixDiagonalize
    USE kp
    USE tetrahedron
    USE om_dos
    USE phi3

    IMPLICIT NONE

    ! To keep a track of updated fc number
    INTEGER :: fc2_1, fc2_2, fc2_3
    INTEGER :: fc3_1, fc3_2, fc3_3
    INTEGER :: neg_lambda,neg_k !index for the most negative eigenvalue
    INTEGER,ALLOCATABLE,DIMENSION(:) :: cutoff_fixed

    COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE :: GradientV_avg !replaced by following MARK
    COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: GradientV_utau
    COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: GradientV_eta
    TYPE(fc2_value),DIMENSION(:,:),ALLOCATABLE :: GradientV_cor
    TYPE(fc2_value),DIMENSION(:,:),ALLOCATABLE :: YY_Compare

    COMPLEX(8) :: V0,F0,V_avg,F_trial,V_avg2
    COMPLEX(8),DIMENSION(:),ALLOCATABLE :: GradientF_trial

    REAL(8) :: conv_to_cart(3,3)
    REAL(8) :: cputime, cputime_1, cputime_2
    REAL(8),DIMENSION(:),ALLOCATABLE :: YY_record,trialfc2_record
    REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE :: approx_term
    LOGICAL :: YY_flag

    REAL(8) :: coefficient(3)!to test Psy<yy>(R+tau) as a coefficient for linear eta
    REAL(8) :: min_eival !record the minimum eigenvalue
    COMPLEX(8), DIMENSION(:,:,:,:),ALLOCATABLE :: phi_test,yy_test

    INTEGER :: ksubset(2), ksubset_F ! for v35 calculation

    COMPLEX(8) :: for_check

    TYPE omega_index !just for record the soft mode index
        INTEGER :: idx_q,idx_la
    END TYPE omega_index
    TYPE(omega_index),DIMENSION(:),ALLOCATABLE :: soft

    CONTAINS
!=========================================FOR INITIAL GUESS=======================================================
    SUBROUTINE initiate_guess
    !!to get a guessing start point for atomic deviation and strain
    !!based on input data
        IMPLICIT NONE
        INTEGER :: i,j,al,be
        INTEGER :: avg_lambda
        REAL(8),DIMENSION(6) ::eta, sig !stress tensor

        REAL(8) :: coeff !temperature dependent coeff to make utau around 0.1

        IF(neg_k.eq.0) RETURN !no too negative eigenvalues found, don't guess

!        coeff = temperature/k_b*100*h_plank*c_light *k_b/ee*max_fc3/max_fc2**2 !/ABS(eivals(neg_lambda,neg_k))
        coeff = max_fc2/max_fc3
!        coeff = SQRT(coeff)
        !to make this number a value between 0.2~1
!        IF(coeff.gt.1d0) THEN
!            coeff = MOD(INT(coeff),100)/100d0
!        END IF

        !initial guess for utau(:,j),linear combination for 3 degenerate
        avg_lambda = 3*INT(neg_lambda/3)
        DO j=2,atom_number
        DO al=1,3
            atomic_deviation(al,j) = REAL(eivecs(al+3*(j-1),avg_lambda+1,neg_k)) + &
            &REAL(eivecs(al+3*(j-1),avg_lambda+2,neg_k)) + REAL(eivecs(al+3*(j-1),avg_lambda+3,neg_k))
        END DO
            atomic_deviation(:,j) = atomic_deviation(:,j)*coeff
        END DO

        !initial guess for strain
        CALL GetStress(sig)

        eta = 0d0
        DO i=1,6
        DO j=1,6
            eta(i) = eta(i) + compliance(i,j)*sig(j)
        END DO
        END DO

        DO i=1,3
        DO j=1,3
            IF(i.eq.j) THEN
                strain(i,j) = -eta(i)
            ELSE
                strain(i,j) = -eta(9-i-j)
            END IF
        END DO
        END DO

        WRITE(*,*) 'atomic deviation',atomic_deviation(:,2)
        WRITE(*,*) 'strain', strain
    END SUBROUTINE initiate_guess
!---------------------------------------------------------------------------

    SUBROUTINE GetStress(sig)
    !!approximate the stress tensor, for the eivecs(tau*xyz,neg_lambda, neg_k)
        IMPLICIT NONE
        INTEGER :: i,j,al,be,ga
        INTEGER :: v1, voigt
        INTEGER :: tau_j,avg_lambda
        REAL(8) :: cell_volume,Rij_be
        REAL(8) :: eivec_diff
        REAL(8),DIMENSION(6),INTENT(OUT) :: sig


        CALL calculate_volume(r1,r2,r3,cell_volume)

        sig = 0d0
        DO al=1,3
        DO be=1,3

        v1 = voigt(al,be)

        DO i=1,atom_number
        DO j=1,tot_atom_number
        DO ga=1,3

            tau_j = every_atom(j)%type_tau
            Rij_be = every_atom(i)%R(be)+every_atom(i)%tau(be)-every_atom(j)%R(be)-every_atom(j)%tau(be)

!            sig(v1) = sig(v1) + myfc2_value(i,j)%phi(al,ga)*&
!            &(eivecs(ga+3*(i-1),neg_lambda,neg_k)-eivecs(ga+3*(tau_j-1),neg_lambda,neg_k))*Rij_be

            avg_lambda = 3*INT(neg_lambda/3)
            eivec_diff = REAL(eivecs(ga+3*(i-1),avg_lambda+1,neg_k)-eivecs(ga+3*(tau_j-1),avg_lambda+1,neg_k)) + &
                & REAL(eivecs(ga+3*(i-1),avg_lambda+2,neg_k)-eivecs(ga+3*(tau_j-1),avg_lambda+2,neg_k)) + &
                & REAL(eivecs(ga+3*(i-1),avg_lambda+3,neg_k)-eivecs(ga+3*(tau_j-1),avg_lambda+3,neg_k))
            sig(v1) = sig(v1) + myfc2_value(i,j)%phi(al,ga)*eivec_diff*Rij_be

        END DO !gamma loop
        END DO !atom j loop
        END DO !atom i loop


        END DO !beta loop
        END DO !alpha loop

        sig = sig/cell_volume

    END SUBROUTINE GetStress
!=================================================================================================================    
    SUBROUTINE Allocate_Gradients
    !! allocate free energy gradients f(:) and <V> gradients
        IMPLICIT NONE

        INTEGER :: size_GradientF_trial
        INTEGER :: i,j

        size_GradientF_trial=sum(variational_parameters_size)

        IF(ALLOCATED(GradientF_trial)) DEALLOCATE(GradientF_trial)
        ALLOCATE(GradientF_trial(size_GradientF_trial))!

        IF(ALLOCATED(GradientV_utau)) DEALLOCATE(GradientV_utau)
        ALLOCATE(GradientV_utau(d,atom_number))

        IF(ALLOCATED(GradientV_eta)) DEALLOCATE(GradientV_eta)
        ALLOCATE(GradientV_eta(d,d))

        IF(ALLOCATED(GradientV_cor)) DEALLOCATE(GradientV_cor)
        ALLOCATE(GradientV_cor(atom_number,tot_atom_number))

        IF(ALLOCATED(YY_Compare)) DEALLOCATE(YY_Compare)
        ALLOCATE(YY_Compare(atom_number,tot_atom_number))
        DO i=1,atom_number
        DO j=1,tot_atom_number
            YY_Compare(i,j)%phi=0d0
        END DO
        END DO

        IF(ALLOCATED(cutoff_fixed)) DEALLOCATE(cutoff_fixed)
        ALLOCATE(cutoff_fixed(ifc2_terms))
    END SUBROUTINE Allocate_Gradients
!=================================================================================================================
    SUBROUTINE GetF_trial_And_GradientF_trial(kvector)
    !!major overall subroutine to
    !!1.call subroutine <initiate_yy> that diagonalize Dynmat, calculate <YY>, etc
    !!2.call subroutines that calculates F0, V0 and <V>
    !!3.combine all the energy gradients into one array f(:)
    !!4.other utilities such as symetrize strain, add pressure, etc
        IMPLICIT NONE

        TYPE(vector),INTENT(INOUT)::kvector(:)

        INTEGER :: size_GradientF_trial,i,j,k,l,temp,temp1,temp2
        INTEGER :: atom1, atom2,R1,R2,tau1,tau2,direction1,direction2
        COMPLEX(8) :: temp_sum,mass1,mass2,compromise(d*atom_number,SIZE(kvector))
        REAL(8) :: negative_eivals,max_eivals,min_eivals
        INTEGER :: k_number

        !**** allocate an 1D array to record <eff_fc2> for L1 norm calculation ****
        IF(ALLOCATED(trialfc2_record)) DEALLOCATE(trialfc2_record)
        ALLOCATE(trialfc2_record(eff_fc2_terms))
        trialfc2_record = 0d0
        j=0
        DO i=1,SIZE(myfc2_index)
            atom1 = myfc2_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2 = myfc2_index(i)%jatom_number
            direction1 = myfc2_index(i)%iatom_xyz
            direction2 = myfc2_index(i)%jatom_xyz
            j = j + 1
            trialfc2_record(j) = trialfc2_value(atom1,atom2)%phi(direction1,direction2)
            !here the trialfc2_value are the new ones after Broyden
        END DO
        !below are to update f(:) using already updated x(:)

        k_number=SIZE(kvector)
        size_GradientF_trial=sum(variational_parameters_size)

        GradientF_trial=0

        !******************ENERGY CALCULATIONS************************
        !we have to diagonalize matrix here
        CALL initiate_yy(kvector)!calculate new <yy> for every time
!CALL check_degeneracy
!CALL special_check
!CALL updatePhi !just for test
!STOP
        CALL GetF0_and_V0 ! not really needed here
        CALL GetV_avg_And_GradientV_avg(kvector) ! translational invariant form that actually works
!        CALL TreatGradientV ! modify GradientV_utau GradientV_eta with a matrix
        F_trial=F0+V_avg-V0

        !*************************************************************

WRITE(*,*)'======================='
WRITE(34,*)'======================='
WRITE(*,*) 'F0=:',REAL(F0)
WRITE(34,*) 'F0=:',REAL(F0)
WRITE(*,*) 'V0=:',REAL(V0)
WRITE(34,*) 'V0=:',REAL(V0)
WRITE(*,*) 'V=:',REAL(V_avg)
WRITE(34,*) 'V=:',REAL(V_avg)
WRITE(*,*) 'F=F0+V-V0:',REAL(F_trial)
WRITE(34,*) 'F=F0+V-V0:',REAL(F_trial)
!WRITE(unitnumber2,'(f10.4)') F_trial

        !**************MAKE GradientF_trial AN ARRAY******************
        !force symmetry on strain gradients?
        CALL sym_strain
        !combine all GradientV_avg and get a 1d array GradientF_trial
        !gradients of free energy w.r.t atomic deviation
        i=0
        DO while(i<variational_parameters_size(1))
            i=i+1
            IF(MOD(i,d).eq.0) THEN
                GradientF_trial(i)=GradientV_utau(d,INT(i/d))
            ELSE
                GradientF_trial(i)=GradientV_utau(MOD(i,d),INT(i/d+1))
            END IF
        END DO !i loop

        !gradients of free energy w.r.t strain
        DO while(i<variational_parameters_size(1)+variational_parameters_size(2))
            i=i+1
            temp=i-variational_parameters_size(1)
            IF(MOD(temp,d).eq.0) THEN
                GradientF_trial(i)=GradientV_eta(INT(temp/d),d)
            ELSE
                GradientF_trial(i)=GradientV_eta(INT(temp/d+1),MOD(temp,d))
            END IF
        END DO !i loop

        !gradients of free energy w.r.t <YY>, then subtract 1/2 effective fc2(which should give 0)
        !...as a substitute of gradients w.r.t effective fc2 itself
        DO while(i<SUM(variational_parameters_size))
            i=i+1
            temp=i-variational_parameters_size(1)-variational_parameters_size(2)
            tau1=eff_fc2_index(temp)%iatom_number
            atom2=eff_fc2_index(temp)%jatom_number
            direction1=eff_fc2_index(temp)%iatom_xyz
            direction2=eff_fc2_index(temp)%jatom_xyz

            GradientF_trial(i)=GradientV_cor(tau1,atom2)%phi(direction1,direction2)-&
            &0.5*trialfc2_value(tau1,atom2)%phi(direction1,direction2)
        END DO !i loop
 WRITE(*,*)       '***************************************************************************'
        !*************************************************************

        !-----add pressure-----
        IF(pressure) THEN
            CALL add_Pressure
        END IF

7 format(a,f10.8)
    END SUBROUTINE GetF_trial_And_GradientF_trial
!================================================================================================================
    SUBROUTINE GetEigen(kvector)
    !!to Fourier transform trial force constant K into 'trialffc2_value'
    !!then diagonalize trial force constant K dynamic matrix 'trialffc2_matrix'
    !!call subroutines in MatrixDiagonalize module
        IMPLICIT NONE
         TYPE(vector),INTENT(IN) :: kvector(:)
         TYPE(ffc2_value),DIMENSION(:,:,:),ALLOCATABLE :: trialffc2_value
         COMPLEX(8),ALLOCATABLE::trialffc2_matrix(:,:),temp(:,:)
         COMPLEX(8),ALLOCATABLE :: vec_test(:,:),vec_test2(:,:),test(:,:),test2(:,:)
         REAL(8) :: limit,val_test(2),qp(d)
         INTEGER :: i,j,k,l,atom1,atom2,R1,R2,tau1,tau2,la
         INTEGER :: n,ier,nv,k_number,ndim

         !---------for test---------
         INTEGER :: test_i,test_j
         INTEGER :: test_l,test_temp
         COMPLEX(8) :: test_sum
         REAL(8) :: test_diff_r,test_diff_i
         !--------------------------

         i=0;k=0;k_number=SIZE(kvector)
         ALLOCATE(trialffc2_value(atom_number,atom_number,k_number))
         ALLOCATE(trialffc2_matrix(d*atom_number,d*atom_number),temp(d*atom_number,d*atom_number))
!****
!ALLOCATE(test(2,2),vec_test(2,2))
!****

         trialffc2_matrix=0

!initialize trialffc2_value as zero, so that we can directly calculate from arbitrarily given trialfc2_value
          Do while(k<k_number)
            k=k+1
            DO i=1,SIZE(myfc2_index)
                !get the R, tau,atom indexes of every fc2
                atom1=myfc2_index(i)%iatom_number !atom index, that's basically R*tau
                atom2=myfc2_index(i)%jatom_number
                R1=every_atom(atom1)%type_R       !cell_vec index,which is the R index
                R2=every_atom(atom2)%type_R
                tau1=every_atom(atom1)%type_tau   !atom type index, which is the tau index
                tau2=every_atom(atom2)%type_tau
                trialffc2_value(tau1,tau2,k)%FFC_D=0
            END DO !i loop
          End Do !k loop

          k=0;i=0

outer:    Do while(k<k_number)
            k=k+1
!------------------------------my old codes, without born but faster------------------------------
!            DO i=1,atom_number
!                tau1=i
!            DO j=1,tot_atom_number
!                tau2=every_atom(j)%type_tau
!                R2=every_atom(j)%type_R
!                trialffc2_value(tau1,tau2,k)%FFC_D = trialffc2_value(tau1,tau2,k)%FFC_D+&
!                &trialfc2_value(tau1,j)%phi*EXP(ci*(kvector(k).dot.cell_vec(:,R2)))/&
!                &SQRT(iatom(tau1)%mass*iatom(tau2)%mass)
!            END DO !j loop
!            END DO !i loop
!
!!we have every value for trialffc2, now need to put it into matrix form: trialffc2_matrix(d*atom_number,d*atom_number)
!!!!NOTICE: trialffc2_matrix has to be an Hermitian!
!            Do i=1,d*atom_number
!              Do j=1,i
!
!                IF(MOD(i,d).eq.0) THEN
!                    IF(MOD(j,d).eq.0) THEN
!                        trialffc2_matrix(i,j)=trialffc2_value(INT(i/d),INT(j/d),k)%FFC_D(d,d)
!                    ELSE
!                        trialffc2_matrix(i,j)=trialffc2_value(INT(i/d),INT(j/d+1),k)%FFC_D(d,MOD(j,d))
!                    END IF
!                ELSE
!                    IF(MOD(j,d).eq.0) THEN
!                        trialffc2_matrix(i,j)=trialffc2_value(INT(i/d+1),INT(j/d),k)%FFC_D(MOD(i,d),d)
!                    ELSE
!                        trialffc2_matrix(i,j)=trialffc2_value(INT(i/d+1),INT(j/d+1),k)%FFC_D(MOD(i,d),MOD(j,d))
!                    END IF
!                END IF
!
!              End Do
!            End Do
!---------------------------------------------------------------------------------------------------
!*********************added born correction to every trialffc2_value*********************
!notes: calculate dynmat & ddyn by <calculate_dynmat> then call <nonanal> to
!add born correction, then let trialffc2_matrix(i,j) equal to new dynmat(i,j)
            qp = kvector(k)%component(:)
            ndim = atom_number*d
            CALL calculate_dynmat(qp)
            !comment the if clause to turn off born term
            !problem at gamma point
            IF (k.ne.1) THEN
                CALL nonanal(qp,dynmat,ndim,ddyn)
            END IF
            DO i=1,ndim
                DO j=1, i
                    trialffc2_matrix(i,j)=dynmat(i,j)
                END DO !j loop
            END DO !i loop

!eigenvalue = 0d0 problem check POSITION 1
!OPEN(75,FILE='DynamicMat.dat',STATUS='unknown',POSITION='append',ACTION='write')
!IF(k.eq.1) THEN
!    WRITE(75,*)'current iteration # ',iter_rec
!    WRITE(75,*)'gamma point dynamic matrix:'
!    DO j=1,ndim
!        WRITE(75,5) REAL(trialffc2_matrix(:,j))
!    END DO
!END IF
!CLOSE(75)
!!***********************************************************************************


!Force Hermitian: fill  in the rest with corresponding conjg
            DO i=1,d*atom_number-1
                DO j=i+1,d*atom_number
                    trialffc2_matrix(i,j)=CONJG(trialffc2_matrix(j,i))
                END DO
            END DO

!fix rounding error or machine error
!maybe this causes the problem of eigenvalue = 0d0?
!            limit=danger
!            DO i=1,d*atom_number
!                DO j=1,d*atom_number
!                    IF(ABS(trialffc2_matrix(i,j)).lt.limit) THEN
!                        trialffc2_matrix(i,j)=0
!                    END IF
!                END DO
!            END DO

!eigenvalue = 0d0 problem check POSITION 2
!OPEN(75,FILE='DynamicMat.dat',STATUS='unknown',POSITION='append',ACTION='write')
!IF(k.eq.1) THEN
!    WRITE(75,*)'gamma point dynamic matrix AFTER:'
!    DO j=1,ndim
!        WRITE(75,5) REAL(trialffc2_matrix(:,j))
!    END DO
!END IF
!CLOSE(75)
!****************************
DO i=1,ndim
DO j=1, ndim
dynmat(i,j) = trialffc2_matrix(i,j)
dynmat_record(i,j,k) = dynmat(i,j)
END DO !
END DO !for test, no other uses
!****************************

!use module zhegv to diagonalize trialffc2_matrix for this k, and get eivecs(:,:,k), eivals(:,k)
            n=d*atom_number
            nv=d*atom_number
            ier=0

            CALL diagonalize(n,trialffc2_matrix,eivals(:,k),nv,eivecs(:,:,k),ier)

!a subroutine to transpose eivecs(atom_type*direction,lambda,k) to eivecs_t(lambda,atom_type*direction,k)
            CALL dagger_eigen(eivecs(:,:,k),eivecs_t(:,:,k))

!--------------------check consistency-----------------------
OPEN(71,FILE='dynmat_check.dat',STATUS='unknown',ACTION='write')
IF(k.eq.1) THEN
    WRITE(71,*)
    WRITE(71,*)'iteration # ',iter_rec
    WRITE(71,*)'=====dynamic matrix index i,j,at which k vector,difference of &
&<from fc2 fourier transform> vs <revert by multiply eigenvec with eigenval>'
END IF
DO test_i=1,ndim
DO test_j=1,ndim
    test_sum = CMPLX(0d0,0d0)
    DO test_l=1,ndim !sum all lambda, for any given (i,j)
        test_sum = test_sum + eivecs(test_i,test_l,k)*eivecs_t(test_l,test_j,k)*eivals(test_l,k)
    END DO
    IF(k.eq.1) THEN
        test_diff_r = ABS(REAL(dynmat(test_i,test_j))-REAL(test_sum))
        test_diff_i = ABS(AIMAG(dynmat(test_i,test_j))-AIMAG(test_sum))
        WRITE(71,9)'real part',test_i,test_j,k,test_diff_r
        WRITE(71,9)'imaginary part',test_i,test_j,k,test_diff_i
    END IF
    IF(test_diff_r.gt.1d-15) THEN
        WRITE(71,*)'real part not match'
        WRITE(71,8)test_i,test_j,k,test_diff_r
        WRITE(71,*)
    ELSEIF(test_diff_i.gt.1d-15) THEN
        WRITE(71,*)'imaginary part not match'
        WRITE(71,8)test_i,test_j,k,test_diff_i
        WRITE(71,*)
    END IF
END DO
END DO
CLOSE(71)
!------------------------------------------------------------

         End Do outer

        DEALLOCATE(trialffc2_matrix)
        DEALLOCATE(trialffc2_value)
        !CLOSE(56)

!****
!test(1,1)=1d0;test(1,2)=2d0;test(2,1)=2d0;test(2,2)=3d0
!CALL diagonalize(2,test,val_test,2,vec_test,ier)
!OPEN(74,FILE='simple_matCheck.dat',STATUS='unknown',ACTION='write')
!WRITE(74,*) 'Eigenvalues:',val_test
!WRITE(74,*) 'First Eigenvector:',vec_test(1,:)
!WRITE(74,*) 'Second Eigenvector:',vec_test(2,:)
!****
!WRITE(*,*) 'eigenvector(:,3,1): ', eivecs(:,3,1)
5 format(6(G16.7,2X))
6 format(2i5,99(1x,g10.4))
7 format((1X,F8.6,A,F8.6,A,2X))
8 format(3(I4,2X),G16.7)
9 format(A,3(I4,2X),G16.7)
    END SUBROUTINE GetEigen
!================================================================================================================
    SUBROUTINE GetF0_and_V0
    !!calculate F0 and <V0>
        IMPLICIT NONE
        INTEGER :: i,k
        F0=0
        V0=0
        DO i=1,SIZE(eivals,DIM=1) !sum for each eigenmode
            DO k=1,SIZE(eivals,DIM=2) !sum for each k
                IF(eivals(i,k).gt.0) THEN
                    !SI unit
                    F0 = F0+100*h_plank*c_light*temperature*&
                    &LOG(2*SINH(SQRT(eivals(i,k)*ee*1d20/uma)/(2*temperature*200*pi*c_light)))
                    V0 = V0+0.25*hbar*SQRT(eivals(i,k)*ee*1d20/uma)/&
                    &TANH(SQRT(eivals(i,k)*ee*1d20/uma)/(2*temperature*200*pi*c_light))
                END IF
                !V0 = V0+0.5*100*h_plank*c_light*temperature !classic
            END DO !k loop
        END DO !i loop
        !energy unit is ev
        F0=F0/SIZE(kvector)/ee
        V0=V0/SIZE(kvector)/ee
    END SUBROUTINE GetF0_and_V0
!================================================================================================================
SUBROUTINE initiate_yy(kvector)
   !!major subroutine that
   !!1.Call diagonalization subroutine
   !!2.apply shift to eigenvalues if negative
   !!3.Call <YY> calculation subroutine
   !!it runs every iteration,  'initiate' maybe improper name
    IMPLICIT NONE
    TYPE(vector),INTENT(IN)::kvector(:)

    INTEGER :: k_number
    INTEGER :: i,j,k,l
    INTEGER :: atom1,xyz1,atom2,xyz2
    INTEGER :: unitnumber,gthb
    REAL(8) :: min_eivals,negative_eivals

    INTEGER :: counter !count how many negative eigenvalues

    counter = 0

    unitnumber=40
    gthb = 41
    OPEN(unitnumber,FILE='eigenvalues.dat',STATUS='unknown',ACTION='write',POSITION='append')
    OPEN(gthb,FILE='eigenvectors.dat',STATUS='unknown',ACTION='write')
    WRITE(unitnumber,*)
    WRITE(gthb,*)
    WRITE(unitnumber,*) 'iteration # = ',iter_rec
    WRITE(gthb,*) 'iteration # = ', iter_rec

    OPEN(47,FILE='Y_square.dat',STATUS='unknown',ACTION='write')
    WRITE(47,*)"This is Iteration #",iter_rec
    WRITE(47,*)"xyz_1,atom_1,xyz_2,atom_2, gamma point correction, <YY> "

    ALLOCATE(phi_test(d,atom_number,d,tot_atom_number))

    k_number = SIZE(kvector)
    CALL GetEigen(kvector)!get eigenvalues from dynmat

    !**** fix if there is any negative eigenvalues based on input fc2 ****
    !*(old)if there are negative eigenvalues but larger than the threshold, just shift will be fine
    !*(1)if there are negative eigenvalues: drop Broyden values, manual update trialfc2 = 2*GradientVCor
    !*(2)if there still are negative eigenvalues, rollback to prev_trial_fc2(freeze eff fc2)*

    min_eivals=MINVAL(MINVAL(eivals, DIM=2))
    min_eival = min_eivals
    WRITE(unitnumber,*)'minimum eigenvalue = ', min_eival
    WRITE(*,*) 'minimum eigenvalue = ',min_eival
    neg_lambda = 0; neg_k = 0 !don't forget to initialize them to be 0 as a flag(no too negative eival)
outer:    DO k=1,SIZE(eivals,dim=2)
        DO l=1,SIZE(eivals,dim=1)

            !below is used for the initial guess generation
            !find the most negative eigenvalue's lambda and k
            IF(min_eivals.lt.-1d-5) THEN
                IF(eivals(l,k).le.min_eivals) THEN
                    neg_lambda = l
                    neg_k = k
                END IF
            END IF

            !print out negative eigenvalues record
            IF(eivals(l,k).le.0d0) THEN
                counter = counter + 1 !find one more soft mode
                negative_eivals=eivals(l,k)
                WRITE(unitnumber,*)
                WRITE(unitnumber,*)'lambda',l,'#kvector',k
                WRITE(unitnumber,*)'eigenvalue=',eivals(l,k)

            END IF
            !print out gamma point eigenvectors record
            IF(k.eq.1) THEN
                WRITE(gthb,*)'lambda',l,'#kvector gamma'
                WRITE(gthb,5) 'atom1 eigenvector (x,y,z)=',REAL(eivecs(1:3,l,1))
                WRITE(gthb,5) 'atom2 eigenvector (x,y,z)=',REAL(eivecs(4:6,l,1))
            END IF
            !very unlikely to happen, that min_eivals are exactly 0
            IF(k.eq.1 .AND. min_eivals.eq.0d0) THEN
                !this could happen to optic bands at gamma point in some weird cases
                IF(eivals(4,1).eq.0d0) eivals(4,1) = 0.0001 * eivals(4,2)
                IF(eivals(5,1).eq.0d0) eivals(5,1) = 0.0001 * eivals(5,2)
                IF(eivals(6,1).eq.0d0) eivals(6,1) = 0.0001 * eivals(6,2)
            END IF

            IF(min_eivals.lt.0d0) THEN
!@@@@@@@@@@@@@@@@@@@@@@Below are how to deal with negative eigenvalues@@@@@@@@@@@@@@@@@@@@@@@@@
            !****option.1 Simple shift method,the only method works stably****
!                eivals(l,k) = eivals(l,k) - min_eivals*1.0001
            !****option.2 The two step fix****
!                IF(min_eivals.gt.-danger) THEN
!                    eivals(l,k) = eivals(l,k) - min_eivals*1.0001
!                ELSE
!                    CALL CheckFixEigen
!                    EXIT outer
!                END IF
            !****option.3 do nothing, but modify Get_Y_square****
            END IF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        END DO
    END DO outer

!------------------------------------------NEW--------------------------------------------------
    !record all the negative eigenvalues index and shift those
    IF(ALLOCATED(soft)) DEALLOCATE(soft)
    IF(counter.ne.0) THEN !if there is any negative eigenvalue
        ALLOCATE(soft(counter))
        i=0
        DO k=1,SIZE(eivals,DIM=2)
        DO l=1,SIZE(eivals, DIM=1)
            IF(eivals(l,k).le.0d0) THEN
                i=i+1
                soft(i)%idx_la = l
                soft(i)%idx_q = k
            END IF
            eivals(l,k) = eivals(l,k) - min_eivals*1.0001
        END DO
        END DO
    END IF
!-----------------------------------------------------------------------------------------------
WRITE(*,*)'minimum eigenvalue after shift:',MINVAL(MINVAL(eivals,DIM=2))

    !**** calculate <YY> and store them in yy_value ****
    DO i=1,SIZE(myfc2_index)
        atom1 = myfc2_index(i)%iatom_number
        IF(atom1.gt.atom_number) CYCLE
        atom2 = myfc2_index(i)%jatom_number
        xyz1 = myfc2_index(i)%iatom_xyz
        xyz2 = myfc2_index(i)%jatom_xyz

        yy_value(atom1,atom2)%phi(xyz1,xyz2) = Get_Y_square(k_number,xyz1,atom1,xyz2,atom2)/k_number


!****Test the value of gamma point correction****
WRITE(47,6) get_letter(xyz1),atom1,get_letter(xyz2),atom2,&
&for_check,yy_value(atom1,atom2)%phi(xyz1,xyz2)
WRITE(47,*)
!****Test if using <yy> can recover the trial fc2****
!phi_test(xyz1,atom1,xyz2,atom2) = Get_phi_test(k_number,xyz1,atom1,xyz2,atom2)
!WRITE(47,8) get_letter(xyz1),atom1,get_letter(xyz2),atom2,phi_test(xyz1,atom1,xyz2,atom2),&
!&trialfc2_value(atom1,atom2)%phi(xyz1,xyz2),&
!&ABS(REAL(phi_test(xyz1,atom1,xyz2,atom2))-trialfc2_value(atom1,atom2)%phi(xyz1,xyz2))

    END DO

DEALLOCATE(phi_test)
!DEALLOCATE(yy_test)
    CLOSE(47)
    CLOSE(unitnumber)
    CLOSE(gthb)

5 FORMAT(a,3(3X,G16.7))
6 FORMAT(2(A2,I3),(G16.7,SP,G16.7,"i"),1(3X,G16.7))
7 FORMAT(2(A2,I3),2(3X,G16.7))
8 FORMAT(2(A2,I3),(G16.7,SP,G16.7,"i"),2(3X,G16.7))
END SUBROUTINE initiate_yy
!================================================================================================================
SUBROUTINE CheckFixEigen
!! fix if there is any negative eigenvalues based on input fc2
!!(1)if there are too negative eigenvalues: drop Broyden values, manual update trialfc2 = 2*GradientVCor*
    IMPLICIT NONE
    INTEGER :: k,l,i,j
    INTEGER :: unitnumber
    REAL(8) :: min_eivals,negative_eivals

    unitnumber=40
    !*(1)
    DO i=1,SIZE(trialfc2_value,dim=1)
    DO j=1,SIZE(trialfc2_value,dim=2)
        trialfc2_value(i,j)%phi = 2*GradientV_cor(i,j)%phi
    END DO
    END DO

    CALL GetEigen(kvector)
    min_eivals=MINVAL(MINVAL(eivals, DIM=2))
    WRITE(*,*) min_eivals
    OPEN(unitnumber,FILE='eigenvalues.dat',STATUS='unknown',ACTION='write',POSITION='append')
    WRITE(unitnumber,*) 'eigenvalues 1st check not passed, enter fixing routing (1)'
    WRITE(unitnumber,*) 'After manual update trialfc2 = 2*GradientVCor*...'
outer:    IF(min_eivals.lt.0d0) THEN
        DO k=1,SIZE(eivals,dim=2)
        DO l=1,SIZE(eivals,dim=1)
            IF(eivals(l,k).lt.0d0) THEN
                negative_eivals=eivals(l,k)
                WRITE(unitnumber,*)
                WRITE(unitnumber,*)'lambda',l,'#kvector',k
                WRITE(unitnumber,*)'eigenvalue=',eivals(l,k)
            END IF
            !if now the min eival is larger than threshold
            IF(min_eivals.lt.0d0 .AND. min_eivals.gt.-danger) THEN
                eivals(l,k) = eivals(l,k) - min_eivals*1.0001
            !if still have too negative eival,rollback trialfc2
            ELSEIF(min_eivals.le.-danger) THEN
                CALL CheckFixEigenAgain
                EXIT outer
            END IF
        END DO
        END DO
    END IF outer

    CLOSE(unitnumber)
END SUBROUTINE CheckFixEigen

SUBROUTINE CheckFixEigenAgain
!!(2)if there still are negative eigenvalues, rollback to prev_trial_fc2(freeze eff fc2)*
    IMPLICIT NONE

    INTEGER :: k,l,i,j
    INTEGER :: unitnumber
    REAL(8) :: min_eivals,negative_eivals

    unitnumber=40

    !*(2)
    trialfc2_value = prev_trialfc2_value

    CALL GetEigen(kvector)
    min_eivals=MINVAL(MINVAL(eivals, DIM=2))
    WRITE(*,*) min_eivals
    OPEN(unitnumber,FILE='eigenvalues.dat',STATUS='unknown',ACTION='write',POSITION='append')
    WRITE(unitnumber,*) 'eigenvalues 2nd check not passed, enter fixing routing (2)'
    WRITE(unitnumber,*) 'After rollback to prev_trial_fc2(freeze eff fc2)...'
    IF(min_eivals.lt.0d0) THEN
        DO k=1,SIZE(eivals,dim=2)
        DO l=1,SIZE(eivals,dim=1)
            IF(eivals(l,k).lt.0d0) THEN
                negative_eivals=eivals(l,k)
                WRITE(unitnumber,*)
                WRITE(unitnumber,*)'lambda',l,'#kvector',k
                WRITE(unitnumber,*)'eigenvalue=',eivals(l,k)
            END IF
            eivals(l,k) = eivals(l,k) - min_eivals*1.0001
        END DO
        END DO
    END IF
    CLOSE(unitnumber)

END SUBROUTINE CheckFixEigenAgain
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNCTION ifExistYY(atom1,atom2) RESULT(found)

!!if <Y(direction1,atom1)Y(direction2,atom2)> exist by translational symmetry
!!return the correspondence atom2 number, else return 0
    IMPLICIT NONE
    INTEGER,INTENT(in) :: atom1,atom2
    INTEGER :: n1,n2,n3,tau
    INTEGER :: i,j
    INTEGER :: found

    n1 = every_atom(atom2)%n1-every_atom(atom1)%n1
    n2 = every_atom(atom2)%n2-every_atom(atom1)%n2
    n3 = every_atom(atom2)%n3-every_atom(atom1)%n3
    tau = every_atom(atom2)%type_tau

    DO i=1,SIZE(every_atom)
        IF(every_atom(i)%n1.ne.n1) CYCLE
        IF(every_atom(i)%n2.ne.n2) CYCLE
        IF(every_atom(i)%n3.ne.n3) CYCLE
        IF(every_atom(i)%type_tau.ne.tau) CYCLE
        EXIT
    END DO
    !In fortran after the DO loop, i will automatically +1
    !which means in the last loop if there is still no match
    !i will be SIZE(every_atom)+1
    IF(i.eq.(SIZE(every_atom)+1)) THEN
        found = 0 !0 means no match
    ELSE
        found = i
    END IF

END FUNCTION ifExistYY
!=========================================================================================================
SUBROUTINE GetV_avg_And_GradientV_avg(kvector)
!!calculate <V> gradients w.r.t strain(:,:), atomic_deviation(:,:), yy_value(:,:)
!!translation invariant ver.
    IMPLICIT NONE

    TYPE(vector),INTENT(IN)::kvector(:)
    LOGICAL :: condition1,condition2,condition3,condition4

    INTEGER :: i=0,j=0,m=0,n=0,o=0,p=0 !atomic loop control
    INTEGER :: new_i, new_j, new_m, new_n, new_o, new_p !mapped new index for myfc_values
    INTEGER :: l=0,l2=0 !eigen loop control
    INTEGER :: k=0,k2=0 !kvector loop control
    INTEGER :: atom1, atom2
    INTEGER :: R1,R2,R3,R4,R5,R6,tau1,tau2,tau3,tau4,tau5,tau6 !atomic indexes
    INTEGER :: tau !atomic indexes for specific gradient
    INTEGER :: found1,found2,found3,found4,found5,found6!found atomic indexes
    INTEGER :: found
    INTEGER :: direction1,direction2,direction3,direction4
    INTEGER :: flag(2)
    INTEGER :: temp1, temp2,k_number

    REAL(8) :: weight(2),cputime
    REAL(8) :: nbe
    REAL(8),DIMENSION(d) :: R1_vec,R2_vec,R3_vec,R4_vec,R5_vec,R6_vec
    REAL(8),DIMENSION(d) :: tau1_vec,tau2_vec,tau3_vec,tau4_vec
    REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE :: Y_square
    !Y_square(direction1,total atom number1,direction2,total atom number2) is the thermal average of correlation dynamic terms
    !Y_square should be real by definition, but in calculations it might be complex due to eivecs
    REAL(8),DIMENSION(:,:),ALLOCATABLE :: S        !S(d,tot_atom_number) is the static displacement
    !--------------------------------------------------------------------
    COMPLEX(8) :: check
    REAL(8) :: threshold
    REAL(8) :: inspect1,inspect2,inspect3,inspect4,inspect_tot
    REAL(8) :: short
    TYPE(fc2_value),DIMENSION(:,:),ALLOCATABLE :: term1,term2,term3,term4
    REAL(8),DIMENSION(3,3) :: term5,term6,term7,term8
    REAL(8),DIMENSION(6,6,6,6) :: C4
    INTEGER :: d1,d2,d3,d4,d5,d6,d7,d8
    INTEGER :: idx1,idx2,idx3,idx4
    INTEGER :: vij,vkl,vik,vjl,vjk,vil

    REAL(8),DIMENSION(6,6) :: C2,A2,A
    INTEGER :: ntindp,rnk, voigt,g,t
    INTEGER :: counter1,counter2

    REAL(8) :: dummy

ALLOCATE(term1(atom_number,tot_atom_number),&
&term2(atom_number,tot_atom_number),&
&term3(atom_number,tot_atom_number),&
&term4(atom_number,tot_atom_number))
DO i=1,atom_number
DO j=1,tot_atom_number
    term1(i,j)%phi = 0d0
    term2(i,j)%phi = 0d0
    term3(i,j)%phi = 0d0
    term4(i,j)%phi = 0d0
END DO
END DO

term5=0d0;term6=0d0;term7=0d0;term8=0d0
C4 = 0d0
C2 = 0d0
A2 = 0d0
A = 0d0
!WRITE(34,*)'===========================CHECK Gradients==========================================='

    ALLOCATE(S(d,tot_atom_number))
    ALLOCATE(Y_square(d,tot_atom_number,d,tot_atom_number))
    Y_square=0d0
    k_number=SIZE(kvector)

!OPEN(52,FILE='check_EC2.dat',STATUS='unknown',POSITION='append',ACTION='write')
!OPEN(54,FILE='check_EC4.dat',STATUS='unknown',POSITION='append',ACTION='write')
!OPEN(79,FILE='GradientVcor.dat',STATUS='unknown',POSITION='append',ACTION='write')
!----------------------------------------------calculate S for every atom------------------------------------------
    DO i=1,tot_atom_number
        S(:,i)=(strain(:,:).dot.(every_atom(i)%R+every_atom(i)%tau))+atomic_deviation(:,every_atom(i)%type_tau)
    END DO
    i=0
    direction1=0
    direction2=0
!--------------------record YY value for threshold calculation-------------------
    IF(ALLOCATED(YY_record)) DEALLOCATE(YY_record)
    ALLOCATE(YY_record(eff_fc2_terms))
    YY_record = 0d0

WRITE(47,*) "===== ALL Y_square at this iteration ====="
j = 0
DO i=1,SIZE(myfc2_index)
    atom1=myfc2_index(i)%iatom_number
    direction1=myfc2_index(i)%iatom_xyz
    atom2=myfc2_index(i)%jatom_number
    direction2=myfc2_index(i)%jatom_xyz

    !assign the iterated yy_value to Y_square
    IF(atom1.le.atom_number) THEN
        Y_square(direction1,atom1,direction2,atom2) = yy_value(atom1,atom2)%phi(direction1,direction2)
        j = j + 1
        !record current <YY> for threshold calculation
        YY_record(j) = Y_square(direction1,atom1,direction2,atom2)
    ELSEIF(atom1.gt.atom_number .AND. atom2.le.atom_number) THEN
        Y_square(direction1,atom1,direction2,atom2) = yy_value(atom2,atom1)%phi(direction2,direction1)!symmetry
    ELSE
        Y_square(direction1,atom1,direction2,atom2) = 0d0
    END IF
WRITE(47,*) get_letter(direction1),atom1,get_letter(direction2),atom2,Y_square(direction1,atom1,direction2,atom2)
END DO
WRITE(47,*) '========================================================'

    !up to here the <YY> ranges (0tau1,R2tau2) have all been calculated

    i=0;j=0;k=0;l=0;direction1=0;direction2=0;atom1=0;atom2=0
!-------------------------------------------initialize V and GradientV---------------------------------------------------------
    V_avg=0;
    GradientV_utau=0
    GradientV_eta=0
    DO i=1,atom_number
        DO j=1,tot_atom_number
            GradientV_cor(i,j)%phi=0
        END DO
    END DO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!----reference check for elastic constants rank 2-----
!CALL check_latfc
!CALL check_read_fc2
!STOP
!CALL read_fc23
!CALL set_huang_inv_constraints
!---- test elastic constants rank 2-------
!C2 = 0d0;A2 = 0d0
!DO i=1,3
!DO j=1,3
!DO k=1,3
!DO l=1,3
!
!vij = voigt(i,j)
!vkl = voigt(k,l)
!vik = voigt(i,k)
!vjl = voigt(j,l)
!
!DO m=1,atom_number
!DO n=1,tot_atom_number
!    !----below is for A2_ikjl----
!    A2(vik,vjl) = A2(vik,vjl) - 0.125*&
!    &myfc2_value(m,n)%phi(i,k)*(every_atom(n)%R(j)+every_atom(n)%tau(j)-every_atom(m)%tau(j))*&
!    &(every_atom(n)%R(l)+every_atom(n)%tau(l)-every_atom(m)%tau(l))&
!    &-0.125*myfc2_value(m,n)%phi(j,l)*(every_atom(n)%R(i)+every_atom(n)%tau(i)-every_atom(m)%tau(i))*&
!    &(every_atom(n)%R(k)+every_atom(n)%tau(k)-every_atom(m)%tau(k))&
!    &-0.125*myfc2_value(m,n)%phi(i,l)*(every_atom(n)%R(j)+every_atom(n)%tau(j)-every_atom(m)%tau(j))*&
!    &(every_atom(n)%R(k)+every_atom(n)%tau(k)-every_atom(m)%tau(k))&
!    &-0.125*myfc2_value(m,n)%phi(j,k)*(every_atom(n)%R(i)+every_atom(n)%tau(i)-every_atom(m)%tau(i))*&
!    &(every_atom(n)%R(l)+every_atom(n)%tau(l)-every_atom(m)%tau(l))
!    !----below is for A_ijkl----
!    A(vij,vkl) = A(vij,vkl) - 0.5*&
!    &myfc2_value(m,n)%phi(i,k)*S(m,j)*S(n,l)
!    !----below is derived from brutal force C_ijkl----
!    short = myfc2_value(m,n)%phi(i,k)*(every_atom(n)%R(l)+every_atom(n)%tau(l))*&
!            &(2*every_atom(m)%tau(j)-every_atom(n)%R(j)-every_atom(n)%tau(j))&
!            &+myfc2_value(m,n)%phi(k,i)*(every_atom(n)%R(j)+every_atom(n)%tau(j))*&
!            &(2*every_atom(m)%tau(l)-every_atom(n)%R(l)-every_atom(n)%tau(l))
!    IF(i.eq.j .AND. k.eq.l) THEN
!        C2(vij,vkl) = C2(vij,vkl) + 0.25*short
!    ELSEIF((i.ne.j .AND. k.eq.l) .OR. (i.eq.j .AND. k.ne.l)) THEN
!        C2(vij,vkl) = C2(vij,vkl) + 0.5*short
!    ELSE
!        C2(vij,vkl) = C2(vij,vkl) + short
!    END IF
!END DO
!END DO
!
!END DO
!END DO
!END DO
!END DO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!---------------------------------------QUADRATIC TERMS-----------------------------------------------------------
    DO i=1,atom_number
        tau1=every_atom(i)%type_tau
        tau1_vec=every_atom(i)%tau
        R1_vec=every_atom(i)%R
    DO j=1,tot_atom_number
        tau2=every_atom(j)%type_tau
        tau2_vec=every_atom(j)%tau
        R2_vec=every_atom(j)%R

        !free energy, two terms
        V_avg=V_avg+0.5*(myfc2_value(tau1,j)%phi.dot.Y_square(:,tau1,:,j))&
            &-0.25*((S(:,j)-S(:,i)).dot.myfc2_value(tau1,j)%phi.dot.(S(:,j)-S(:,i)))

        !gradient w.r.t <YY>, one term
        GradientV_cor(tau1,j)%phi = GradientV_cor(tau1,j)%phi+0.5*myfc2_value(tau1,j)%phi

        DO direction1=1,d

        !gradient w.r.t atomic deviation
        GradientV_utau(direction1,tau1)=GradientV_utau(direction1,tau1)&
        &+(myfc2_value(tau1,j)%phi(direction1,:).dot.(S(:,j)-S(:,tau1)))

        DO direction2=1,d

        !gradient w.r.t strain

        GradientV_eta(direction1,direction2)=GradientV_eta(direction1,direction2)&
        &+(myfc2_value(tau1,j)%phi(direction1,:).dot.S(:,j))&
        &*(0.5*tau1_vec(direction2)-0.25*R2_vec(direction2)-0.25*tau2_vec(direction2))&
        &+(myfc2_value(tau1,j)%phi(:,direction1).dot.(0.5*S(:,i)-0.25*S(:,j)))&
        &*(R2_vec(direction2)+tau2_vec(direction2))&
        &+(myfc2_value(tau1,j)%phi(direction2,:).dot.S(:,j))&
        &*(0.5*tau1_vec(direction1)-0.25*R2_vec(direction1)-0.25*tau2_vec(direction1))&
        &+(myfc2_value(tau1,j)%phi(:,direction2).dot.(0.5*S(:,i)-0.25*S(:,j)))&
        &*(R2_vec(direction1)+tau2_vec(direction1))

        END DO !direction2
        END DO !direction1

    END DO !j loop *atomic
    END DO !i loop *atomic

WRITE(*,*) 'FINISH QUADRATIC TERMS CALCULATION'

check = V_avg
WRITE(*,*) 'check quadratic:', check

coefficient(:) = 0d0
!-------------------------------------------CUBIC TERMS------------------------------------------------------
    DO i=1,atom_number
        tau1=every_atom(i)%type_tau
        tau1_vec=every_atom(i)%tau
        R1_vec=every_atom(i)%R
    DO j=1,tot_atom_number
        IF(.NOT.ANY(fc3_unique_idx==j)) CYCLE !skip if j not valid
        tau2=every_atom(j)%type_tau
        tau2_vec=every_atom(j)%tau
        R2_vec=every_atom(j)%R
    DO m=1,tot_atom_number
        IF(.NOT.ANY(fc3_unique_idx==m)) CYCLE !skip if m not valid
        tau3=every_atom(m)%type_tau
        tau3_vec=every_atom(m)%tau
        R3_vec=every_atom(m)%R

        !map j,m for valid fc3 search
        new_j=find_loc(fc3_unique_idx,j)
        new_m=find_loc(fc3_unique_idx,m)
        !free energy
        V_avg=V_avg-1d0/12*(myfc3_value(tau1,new_j,new_m)%psi.dot.(S(:,m)-S(:,i)).dot. &
        &(S(:,j)-S(:,i)).dot.(S(:,j)-S(:,i)))&
        &+1d0/2*(myfc3_value(tau1,new_j,new_m)%psi.dot.(S(:,m)-S(:,i)).dot.Y_square(:,tau1,:,j))

        !gradient w.r.t <YY>
        GradientV_cor(tau1,j)%phi=GradientV_cor(tau1,j)%phi+0.5*(myfc3_value(tau1,new_j,new_m)%psi&
                                     &.dot.S(:,m))

        DO direction1=1,d
        !gradient w.r.t atomic deviation, term 1
        GradientV_utau(direction1,tau1)=GradientV_utau(direction1,tau1)&
        &+1d0/2*(myfc3_value(tau1,new_j,new_m)%psi(direction1,:,:).dot.(S(:,m)-S(:,tau1)).dot.(S(:,j)-S(:,tau1)))

        !gradient w.r.t atomic deviation, term 2
        GradientV_utau(direction1,tau3)=GradientV_utau(direction1,tau3)&
        &+1d0/2*(myfc3_value(i,new_j,new_m)%psi(:,:,direction1).dot.Y_square(:,i,:,j))

        DO direction2=1,d

        !gradient w.r.t strain, term 1

        GradientV_eta(direction1,direction2)=GradientV_eta(direction1,direction2)&
        &-1d0/12*(myfc3_value(tau1,new_j,new_m)%psi(direction1,:,:).dot.(S(:,m)-S(:,i)).dot.(S(:,j)-S(:,i)))&
        &*(R2_vec(direction2)+tau2_vec(direction2)-tau1_vec(direction2))&
        &-1d0/12*(myfc3_value(tau1,new_j,new_m)%psi(direction2,:,:).dot.(S(:,m)-S(:,i)).dot.(S(:,j)-S(:,i)))&
        &*(R2_vec(direction1)+tau2_vec(direction1)-tau1_vec(direction1))&
        &-1d0/12*(myfc3_value(tau1,new_j,new_m)%psi(:,direction1,:).dot.(S(:,m)-S(:,i)).dot.(S(:,j)-S(:,i)))&
        &*(R2_vec(direction2)+tau2_vec(direction2)-tau1_vec(direction2))&
        &-1d0/12*(myfc3_value(tau1,new_j,new_m)%psi(:,direction2,:).dot.(S(:,m)-S(:,i)).dot.(S(:,j)-S(:,i)))&
        &*(R2_vec(direction1)+tau2_vec(direction1)-tau1_vec(direction1))&
        &-1d0/12*(myfc3_value(tau1,new_j,new_m)%psi(:,:,direction1).dot.(S(:,j)-S(:,i)).dot.(S(:,j)-S(:,i)))&
        &*(R3_vec(direction2)+tau3_vec(direction2)-tau1_vec(direction2))&
        &-1d0/12*(myfc3_value(tau1,new_j,new_m)%psi(:,:,direction2).dot.(S(:,j)-S(:,i)).dot.(S(:,j)-S(:,i)))&
        &*(R3_vec(direction1)+tau3_vec(direction1)-tau1_vec(direction1))&
        &+0.5d0*(myfc3_value(tau1,new_j,new_m)%psi(:,:,direction1).dot.Y_square(:,tau1,:,j))&
        &*(R3_vec(direction2)+tau3_vec(direction2)-tau1_vec(direction2))&
        &+0.5d0*(myfc3_value(tau1,new_j,new_m)%psi(:,:,direction2).dot.Y_square(:,tau1,:,j))&
        &*(R3_vec(direction1)+tau3_vec(direction1)-tau1_vec(direction1))

!----------------------------------just for test-------------------------------
!        IF(direction1.ne.direction2) THEN
!            dummy = 0.5d0*(myfc3_value(tau1,new_j,new_m)%psi(:,:,direction1).dot.myfc2_value(tau1,j)%phi(:,:))&
!                &*(R3_vec(direction2)+tau3_vec(direction2)-tau1_vec(direction2))&
!                &+0.5d0*(myfc3_value(tau1,new_j,new_m)%psi(:,:,direction2).dot.myfc2_value(tau1,j)%phi(:,:))&
!                &*(R3_vec(direction1)+tau3_vec(direction1)-tau1_vec(direction1))
!            IF(direction1.eq.1 .AND. direction2.eq.2) THEN !eta_xy
!                coefficient(1) = coefficient(1) + dummy
!            ELSEIF(direction1.eq.1 .AND. direction2.eq.3) THEN !eta_xz
!                coefficient(2) = coefficient(2) + dummy
!            ELSEIF(direction1.eq.2 .AND. direction2.eq.3) THEN  !eta_yz
!                coefficient(3) = coefficient(3) + dummy
!            END IF
!        END IF
!-------------------------------------------------------------------------------
        !gradient w.r.t strain,term 2
        found = ifExistYY(j,m)
        IF(found.ne.0) THEN

        END IF

        END DO !direction2 loop
        END DO !direction1 loop

    END DO !m loop *atomic
    END DO !j loop *atomic
    END DO !i loop *atomic

WRITE(*,*) 'FINISH CUBIC TERMS CALCULATION'

check = V_avg - check
WRITE(*,*) 'check cubic:', check

!------------------------------------------------QUARTIC TERMS-----------------------------------------------------
    DO i=1,atom_number
        tau1=every_atom(i)%type_tau
        tau1_vec=every_atom(i)%tau
        R1_vec=every_atom(i)%R
    DO j=1,tot_atom_number
        IF(.NOT.ANY(fc4_unique_idx==j)) CYCLE !skip if j not valid
        tau2=every_atom(j)%type_tau
        tau2_vec=every_atom(j)%tau
        R2_vec=every_atom(j)%R
    DO m=1,tot_atom_number
        IF(.NOT.ANY(fc4_unique_idx==m)) CYCLE !skip if m not valid
        tau3=every_atom(m)%type_tau
        tau3_vec=every_atom(m)%tau
        R3_vec=every_atom(m)%R
    DO n=1,tot_atom_number
        IF(.NOT.ANY(fc4_unique_idx==n)) CYCLE !skip if n not valid
        tau4=every_atom(n)%type_tau
        tau4_vec=every_atom(n)%tau
        R4_vec=every_atom(n)%R

        !map j,m,n for valid fc4 search
        new_j=find_loc(fc4_unique_idx,j)
        new_m=find_loc(fc4_unique_idx,m)
        new_n=find_loc(fc4_unique_idx,n)


        !free energy, term 1&2
        V_avg = V_avg - 1d0/48*(myfc4_value(tau1,new_j,new_m,new_n)%chi.dot.(S(:,n)-S(:,i)).dot.(S(:,m)-S(:,i)).dot. &
        &(S(:,j)-S(:,i)).dot.(S(:,j)-S(:,i)))+1d0/4*(myfc4_value(tau1,new_j,new_m,new_n)%chi.dot.(S(:,n)-S(:,i)).dot. &
        &(S(:,m)-S(:,i)).dot.Y_square(:,tau1,:,j))


        !free energy, term 3
        found = ifExistYY(m,n)
        IF(found.ne.0) THEN
            V_avg = V_avg + 3d0/24*(myfc4_value(tau1,new_j,new_m,new_n)%chi.dot.Y_square(:,tau3,:,found).dot.Y_square(:,tau1,:,j))
        END IF

        !gradient w.r.t <YY>, term 1
        GradientV_cor(tau1,j)%phi=GradientV_cor(tau1,j)%phi&
        &+1d0/4*(myfc4_value(tau1,new_j,new_m,new_n)%chi.dot.(S(:,n)-S(:,i)).dot.(S(:,m)-S(:,i)))

        !gradient w.r.t <YY>, term 2
        found = ifExistYY(m,n)
        IF(found.ne.0) THEN
            GradientV_cor(tau1,j)%phi=GradientV_cor(tau1,j)%phi&
            &+3d0/24*(myfc4_value(tau1,new_j,new_m,new_n)%chi.dot.Y_square(:,tau3,:,found))

            GradientV_cor(tau3,found)%phi=GradientV_cor(tau3,found)%phi&
            &+3d0/24*(Y_square(:,i,:,j).dot.myfc4_value(tau1,new_j,new_m,new_n)%chi)
        END IF



        DO direction1=1,d

        !gradient w.r.t atomic deviation, term 1
        GradientV_utau(direction1,tau1)=GradientV_utau(direction1,tau1)&
        &+1d0/6*(myfc4_value(tau1,new_j,new_m,new_n)%chi(direction1,:,:,:).dot. &
        &(S(:,n)-S(:,tau1)).dot.(S(:,m)-S(:,tau1)).dot.(S(:,j)-S(:,tau1)))

        !gradient w.r.t atomic deviation, term 2
        GradientV_utau(direction1,tau3)=GradientV_utau(direction1,tau3)&
        &+1d0/2*(myfc4_value(i,new_j,new_m,new_n)%chi(:,:,direction1,:).dot.(S(:,n)-S(:,i)).dot.Y_square(:,i,:,j))

        DO direction2=1,d

        GradientV_eta(direction1,direction2)=GradientV_eta(direction1,direction2)&
        &-1d0/48*(myfc4_value(tau1,new_j,new_m,new_n)%chi(direction1,:,:,:).dot. &
        &(S(:,n)-S(:,i)).dot.(S(:,m)-S(:,i)).dot.(S(:,j)-S(:,i)))*(R2_vec(direction2)+tau2_vec(direction2)-tau1_vec(direction2))&
        &-1d0/48*(myfc4_value(tau1,new_j,new_m,new_n)%chi(direction2,:,:,:).dot. &
        &(S(:,n)-S(:,i)).dot.(S(:,m)-S(:,i)).dot.(S(:,j)-S(:,i)))*(R2_vec(direction1)+tau2_vec(direction1)-tau1_vec(direction1))&
        &-1d0/16*(myfc4_value(tau1,new_j,new_m,new_n)%chi(:,direction1,:,:).dot. &
        &(S(:,n)-S(:,i)).dot.(S(:,m)-S(:,i)).dot.(S(:,j)-S(:,i)))*(R2_vec(direction2)+tau2_vec(direction2)-tau1_vec(direction2))&
        &-1d0/16*(myfc4_value(tau1,new_j,new_m,new_n)%chi(:,direction2,:,:).dot. &
        &(S(:,n)-S(:,i)).dot.(S(:,m)-S(:,i)).dot.(S(:,j)-S(:,i)))*(R2_vec(direction1)+tau2_vec(direction1)-tau1_vec(direction1))&
        &+0.5*(myfc4_value(tau1,new_j,new_m,new_n)%chi(:,:,direction1,:).dot. &
        &(S(:,n)-S(:,i)).dot.Y_square(:,tau1,:,j))*(R3_vec(direction2)+tau3_vec(direction2)-tau1_vec(direction2))&
        &+0.5*(myfc4_value(tau1,new_j,new_m,new_n)%chi(:,:,direction2,:).dot. &
        &(S(:,n)-S(:,i)).dot.Y_square(:,tau1,:,j))*(R3_vec(direction1)+tau3_vec(direction1)-tau1_vec(direction1))

        !gradient w.r.t strain, term 2
        found = ifExistYY(m,n)
        IF(found.ne.0) THEN

        END IF

        END DO !direction2
        END DO !direction1


    END DO !n loop *atomic
    END DO !m loop *atomic
    END DO !j loop *atomic
    END DO !i loop *atomic

WRITE(*,*)'FINISH QUARTIC TERMS CALCULATION'

check = V_avg - check
WRITE(*,*) 'check quartic:', check

!-------------fix the gradient of eta if direction1=direction2----------
DO direction1=1,d
DO direction2=1,d
    IF(direction1.eq.direction2) THEN
        GradientV_eta(direction1,direction2) = 0.5*GradientV_eta(direction1,direction2)
    END IF
END DO
END DO

!**************Check Specific Term Block****************

!-----------------convert A into C according to Wallace------------------
!DO i=1,3
!DO j=1,3
!DO k=1,3
!DO l=1,3
!    vij = voigt(i,j)
!    vkl = voigt(k,l)
!    vil = voigt(i,l)
!    vjk = voigt(j,k)
!    vik = voigt(i,k)
!    vjl = voigt(j,l)
!
!    C2(vij,vkl) = A2(vik,vjl)+A2(vjk,vil)-A2(vij,vkl)
!
!END DO
!END DO
!END DO
!END DO
!------------------------------------------------------------------------
!WRITE(54,*)"elastic constants rank 4 check"
!WRITE(52,*)"elastic constants rank 2 check"
!WRITE(52,*)"idx1,idx2,A_hat,C,A_tilda"
!DO idx1=1,6
!DO idx2=1,6
!    IF(ABS(A2(idx1,idx2)).gt.1d-8 .OR.ABS(C2(idx1,idx2)).gt.1d-8&
!    & .OR.ABS(A(idx1,idx2)).gt.1d-8) THEN
!        WRITE(52,6) idx1,', ',idx2, ', ',&
!        &A2(idx1,idx2),C2(idx1,idx2),A(idx1,idx2)
!    END IF
!DO idx3=1,6
!DO idx4=1,6
!    IF(ABS(C4(idx1,idx2,idx3,idx4)).gt.1d-10) THEN
!        WRITE(54,*) idx1,', ',idx2, ', ',&
!        &idx3,', ',idx4,', ',&
!        &C4(idx1,idx2,idx3,idx4)
!    END IF
!
!END DO
!END DO
!END DO
!END DO
!
!WRITE(52,*) '  eigenvalues for C2 matrix'
!WRITE(52,*)' eig1 = 2*V0*(C11-C12)=',2*(C2(1,1)-C2(1,2))
!WRITE(52,*)' eig2 = 3*V0*(C11+2*C12)=',3*(C2(1,1)+2*C2(1,2))
!WRITE(52,*)' eig3 = 4*V0*C44=',4*C2(4,4)

!WRITE(54,*)"strain gradients terms check"
!WRITE(54,*)"===== xyz1, xyz2, term5, term6, term7, term8 ===== "
!DO direction1=1,3
!DO direction2=1,3
!    WRITE(54,*) get_letter(direction1),',  ',get_letter(direction2),',  ',&
!    & term5(direction1,direction2),',  ', term6(direction1,direction2),',  ',&
!    & term7(direction1,direction2),',  ', term8(direction1,direction2)
!END DO
!END DO
!WRITE(54,*)'============================================'
!DO i=1,atom_number
!DO direction1=1,d
!DO direction2=1,d
!    inspect1 = 0d0; inspect2 = 0d0; inspect3 = 0d0; inspect4 = 0d0
!    inspect_tot = 0d0
!DO j=1,tot_atom_number
!    inspect1 = inspect1 + term1(i,j)%phi(direction1,direction2)
!    inspect2 = inspect2 + term2(i,j)%phi(direction1,direction2)
!    inspect3 = inspect3 + term3(i,j)%phi(direction1,direction2)
!    inspect4 = inspect4 + term4(i,j)%phi(direction1,direction2)
!    inspect_tot = inspect_tot + GradientV_cor(i,j)%phi(direction1,direction2)
!END DO
!
!threshold = 1e-8
!
!IF(ABS(inspect1).gt.threshold) THEN
!    WRITE(54,*) "term1 sum failed ASR check: ", inspect1
!    WRITE(54,*) i,get_letter(direction1),get_letter(direction2)
!    WRITE(54,*) '-----------------------------------------'
!END IF
!
!IF(ABS(inspect2).gt.threshold) THEN
!    WRITE(54,*) "term2 sum failed ASR check: ", inspect2
!    WRITE(54,*) i,get_letter(direction1),get_letter(direction2)
!    WRITE(54,*) '-----------------------------------------'
!END IF
!
!IF(ABS(inspect3).gt.threshold) THEN
!    WRITE(54,*) "term3 sum failed ASR check: ", inspect3
!    WRITE(54,*) i,get_letter(direction1),get_letter(direction2)
!    WRITE(54,*) '-----------------------------------------'
!END IF
!IF(ABS(inspect4).gt.threshold) THEN
!    WRITE(54,*) "term4 sum failed ASR check: ", inspect4
!    WRITE(54,*) i,get_letter(direction1),get_letter(direction2)
!    WRITE(54,*) '-----------------------------------------'
!END IF
!
!IF(inspect_tot.ne.0) THEN
!    WRITE(54,*) "gradients sum failed ASR check: ", inspect_tot
!    WRITE(54,*) i,get_letter(direction1),get_letter(direction2)
!    WRITE(54,*) '-----------------------------------------'
!END IF
!
!END DO
!END DO
!END DO
!WRITE(54,*)'============================================'
!
!rnk = 2
!DO j=1, map(rnk)%ngr
!    ntindp = map(rnk)%ntind(j)
!    WRITE(54,*)'----group',j,'-----'
!    DO i=1, map(rnk)%nt(j)
!        atom1 = map(rnk)%gr(j)%iat(1,i)
!        IF(atom1.gt.atom_number) CYCLE
!        atom2 = map(rnk)%gr(j)%iat(2,i)
!        direction1 = map(rnk)%gr(j)%ixyz(1,i)
!        direction2 = map(rnk)%gr(j)%ixyz(2,i)
!        WRITE(54,'(2(a10,f12.9),f7.3)')'term1: ' ,term1(atom1,atom2)%phi(direction1,direction2),&
!        &'  term2: ',term2(atom1,atom2)%phi(direction1,direction2),&
!        &map(rnk)%gr(j)%mat(i,1:ntindp)
!        WRITE(54,*)
!        WRITE(54,'(2(a10,f12.9),f7.3)')'term3: ' ,term3(atom1,atom2)%phi(direction1,direction2),&
!        &'  term4: ',term4(atom1,atom2)%phi(direction1,direction2),&
!        &map(rnk)%gr(j)%mat(i,1:ntindp)
!        WRITE(54,*)
!        WRITE(54,*) 'gradient <YY>: ',GradientV_cor(atom1,atom2)%phi(direction1,direction2)
!        WRITE(54,*)
!    END DO
!END DO

!DO i=1,SIZE(myfc2_index)
!    atom1 = myfc2_index(i)%iatom_number
!    IF(atom1.gt.atom_number) CYCLE
!    atom2 = myfc2_index(i)%jatom_number
!    direction1 = myfc2_index(i)%iatom_xyz
!    direction2 = myfc2_index(i)%jatom_xyz
!
!    WRITE(79,7) get_letter(direction1),atom1,get_letter(direction2),atom2,&
!    &GradientV_cor(atom1,atom2)%phi(direction1,direction2),&
!    &myfc2_value(atom1,atom2)%phi(direction1,direction2)
!END DO
!
!DEALLOCATE(term1,term2,term3,term4)
!!CLOSE(52)
!!CLOSE(54)
!CLOSE(79)
!STOP
!*******************************************************

    DEALLOCATE(S)
    DEALLOCATE(Y_square)
CLOSE(47)

5 format(4i4,99(4x,g10.4))
6 format(2(i4,a),3(2x,g11.4))
7 FORMAT(2(A2,I3),2(3X,G16.7))
END SUBROUTINE GetV_avg_And_GradientV_avg
!=====================================================================================================================================
    SUBROUTINE CheckFixMatrix(negative_eivals)
    !!Fix by shift all the matrix elements, just for test purpose
        IMPLICIT NONE
        INTEGER :: k=0,l=0,i,j
        INTEGER :: atom1,atom2,tau1,tau2
        REAL(8) :: mass1,mass2
        REAL(8),INTENT(IN) :: negative_eivals

        !shift the trial fc2 matrix
        DO i=1,fc_terms(2)
            atom1=myfc2_index(i)%iatom_number
            atom2=myfc2_index(i)%jatom_number
            tau1=every_atom(atom1)%type_tau
            tau2=every_atom(atom2)%type_tau
            mass1=iatom(tau1)%mass
            mass2=iatom(tau2)%mass
            trialfc2_value(atom1,atom2)%phi=trialfc2_value(atom1,atom2)%phi-1.1*negative_eivals*SQRT(mass1*mass2)
        END DO
        !Force ASR
        DO i=1,atom_number
            trialfc2_value(i,i)%phi=0
            DO j=1,tot_atom_number
                IF(j.ne.i) THEN
                    trialfc2_value(i,i)%phi=trialfc2_value(i,i)%phi-trialfc2_value(i,j)%phi
!WRITE(*,*)'Eigenvalues=',eivals
!WRITE(*,*)'Is It Large? trialfc2(i,i)=',trialfc2_value(i,i)%phi
                END IF
            END DO
        END DO

    END SUBROUTINE CheckFixMatrix
!--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE CheckFixMatrix2(negative_eivals)
        !!directly shift the eigenvalues, just for test purpose
        IMPLICIT NONE
        INTEGER :: k=0,l=0,i,j
        REAL(8),INTENT(IN) :: negative_eivals

        DO l=1,SIZE(eivals,DIM=1)
            DO k=1,SIZE(eivals,DIM=2)
                eivals(l,k)=eivals(l,k)-negative_eivals+0.01
            END DO !k loop
        END DO !l loop
    END SUBROUTINE CheckFixMatrix2
!--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE CheckFixMatrix3(min_eivals,max_eivals)
    !!linear transformation the eigenvalues, just for test purpose
        IMPLICIT NONE
        INTEGER :: k=0,l=0,i,j
        REAL(8),INTENT(IN) :: min_eivals,max_eivals

        DO l=1,SIZE(eivals,DIM=1)
            DO k=1,SIZE(eivals,DIM=2)
                eivals(l,k)=(eivals(l,k)-min_eivals)/(max_eivals-min_eivals)*max_eivals+1E-10
            END DO !k loop
        END DO !l loop

    END SUBROUTINE CheckFixMatrix3
!----------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE CheckFixASR
    !!force ASR in trialfc2 K
        IMPLICIT NONE
        TYPE(fc2_value),DIMENSION(:),ALLOCATABLE :: check_trialfc2
        INTEGER :: i,j,k,direction1,direction2,temp=0

        !OPEN(55,FILE='CheckASR.dat',STATUS='unknown',ACTION='write')

        ALLOCATE(check_trialfc2(atom_number))
        DO i=1,atom_number
            check_trialfc2(i)%phi=0
        END DO

        DO i=1,atom_number
            DO direction1=1,d
                DO direction2=1,d
                    DO j=1,tot_atom_number
                        check_trialfc2(i)%phi(direction1,direction2)=check_trialfc2(i)%phi(direction1,direction2)+&
                        &    trialfc2_value(i,j)%phi(direction1,direction2)
                    END DO
                    IF(check_trialfc2(i)%phi(direction1,direction2).ne.0) THEN
                        temp=temp+1
                        !WRITE(55,*) 'BEFORE',trialfc2_value(i,i)%phi(direction1,direction2)
                    END IF
                END DO
            END DO
        END DO
        IF(temp.eq.0) THEN
            WRITE(*,*)'trial fc2 satisfy ASR!'
        ELSE
            WRITE(*,*) 'trial fc2 do not satisfy ASR, try fixing'
            DO i=1,atom_number
                DO direction1=1,d
                    DO direction2=1,d
                        trialfc2_value(i,i)%phi(direction1,direction2)=trialfc2_value(i,i)%phi(direction1,direction2)-&
                        &                                 check_trialfc2(i)%phi(direction1,direction2)
                        check_trialfc2(i)%phi(direction1,direction2)=0
                        !WRITE(55,*) 'AFTER',trialfc2_value(i,i)%phi(direction1,direction2)
                    END DO
                END DO
            END DO
        END IF

        !CLOSE(55)
        DEALLOCATE(check_trialfc2)
    END SUBROUTINE CheckFixASR
!======================================================================================================================================
SUBROUTINE select_atom(R1_vec,R2_vec,type_tau,condition,found)
!!given R1 and R2 and atom type, find if there is an atom at R1+R2 with same atom type, and give its atom index
    IMPLICIT NONE
    REAL(8),INTENT(IN),DIMENSION(d) :: R1_vec,R2_vec
    REAL(8),DIMENSION(d) :: R_combine,R
    INTEGER,INTENT(IN) :: type_tau
    INTEGER,INTENT(OUT) :: found
    LOGICAL,INTENT(OUT) :: condition
    INTEGER :: i,j,tau

    condition=.FALSE.
    R_combine=R1_vec+R2_vec
outer: DO i=1,tot_atom_number
            R=every_atom(i)%R
            tau=every_atom(i)%type_tau
            IF(type_tau.eq.tau) THEN
    inner:      DO j=1,SIZE(R)
                    IF(R_combine(j).ne.R(j)) THEN
                        EXIT inner
                    ELSE
                        IF(j.eq.SIZE(R)) THEN
                            found=i
                            condition=.TRUE.
                            EXIT outer
                        END IF
                    END IF
                END DO inner
            END IF
       END DO outer

END SUBROUTINE select_atom
!======================================================================================================================================
SUBROUTINE ConvertAtomicIndex(R1_vec,R2_vec,type_tau,found)
!!duplicate version of previous subroutine, only return atom index
!!Need to take tau2 index into consideration
    IMPLICIT NONE
    REAL(8),INTENT(IN), DIMENSION(d) :: R1_vec, R2_vec
    REAL(8),DIMENSION(d) :: R_combine,R
    INTEGER, INTENT(IN) :: type_tau
    INTEGER,INTENT(OUT) :: found
    INTEGER :: i,j,tau

    R_combine=R1_vec+R2_vec
outer:  DO i=1,tot_atom_number
            R=every_atom(i)%R
            tau=every_atom(i)%type_tau
            IF(type_tau.eq.tau) THEN
      inner:    DO j=1,SIZE(R)
                    IF(R_combine(j).ne.R(j)) THEN
                        exit inner
                    ELSE
                        IF(j.eq.SIZE(R)) THEN
                            found=i !found the ith atom
                            exit outer
                        END IF
                    END IF
                END DO inner
            END IF
        END DO outer
END SUBROUTINE ConvertAtomicIndex
!-------------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CheckRange(R1_vec,R2_vec,condition)
!!duplicate version of previous subroutine, only return if found or not
!!just to check if the corresponding cell is inside the farthest shell or not
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: condition
    REAL(8),INTENT(IN), DIMENSION(d) :: R1_vec, R2_vec
    REAL(8),DIMENSION(d) :: R_combine,R
    INTEGER :: i,j=0

    R_combine=R1_vec+R2_vec
WRITE(*,*) 'Check R_combine',R_combine
outer:  DO i=1, tot_atom_number
            R=every_atom(i)%R
inner:      DO j=1, SIZE(R)
                IF(R_combine(j).ne.R(j)) THEN
                    condition=.FALSE.
                    exit inner
                    ELSE
                        IF(j.eq.SIZE(R)) THEN
                            condition=.TRUE.
                            exit outer
                        END IF
                END IF
            END DO inner
        END DO outer
END SUBROUTINE CheckRange
!======================================================================================================================================
    SUBROUTINE calculate_dynmat(qp)
    !! calculate dynmat & ddyn at q-point qp(d)
    !! dynmat & ddyn are global variable declared in module [Fourier_force_constants]
        IMPLICIT NONE

        INTEGER :: i,j,k,l
        INTEGER :: atom1,atom2,R1,R2,tau1,tau2,n,nv
        REAL(8) :: limit
        REAL(8),INTENT(IN) :: qp(d)
        REAL(8) :: rr(3)
        TYPE(ffc2_value),DIMENSION(:,:),ALLOCATABLE :: dynmat_temp
        TYPE(ffc2_value),DIMENSION(:,:,:),ALLOCATABLE :: ddyn_temp

        ALLOCATE(dynmat_temp(atom_number,atom_number))
        ALLOCATE(ddyn_temp(atom_number,atom_number,d))

        !initialize ddyn_temp,dynmat_temp
        DO i=1,fc_terms(2)
            !get the R, tau,atom indexes
            atom1=myfc2_index(i)%iatom_number !atom index, that's basically R*tau
            atom2=myfc2_index(i)%jatom_number
            R1=every_atom(atom1)%type_R       !cell_vec index,which is the R index
            R2=every_atom(atom2)%type_R
            tau1=every_atom(atom1)%type_tau   !atom type index, which is the tau index
            tau2=every_atom(atom2)%type_tau
            dynmat_temp(tau1,tau2)%FFC_D=CMPLX(0d0,0d0)
            DO j=1,d
                ddyn_temp(tau1,tau2,j)%FFC_D=CMPLX(0d0,0d0)
            END DO !j loop
        END DO !i loop

        !calculate ddyn_temp,dynmat_temp
         DO i=1,atom_number
            tau1=i
            DO j=1,tot_atom_number
                tau2=every_atom(j)%type_tau
                R2=every_atom(j)%type_R
                dynmat_temp(tau1,tau2)%FFC_D = dynmat_temp(tau1,tau2)%FFC_D+&
                &trialfc2_value(tau1,j)%phi*EXP(ci*(qp.dot.cell_vec(:,R2)))/&
                &SQRT(iatom(tau1)%mass*iatom(tau2)%mass)

                !this below is for test, should be commented out
!                rr = every_atom(j)%R + every_atom(j)%tau - every_atom(i)%R - every_atom(i)%tau
!                dynmat_temp(tau1,tau2)%FFC_D = dynmat_temp(tau1,tau2)%FFC_D+&
!                &trialfc2_value(tau1,j)%phi*EXP(ci*(qp.dot.rr))/&
!                &SQRT(iatom(tau1)%mass*iatom(tau2)%mass)

                DO l=1,d
                    ddyn_temp(tau1,tau2,l)%FFC_D = ddyn_temp(tau1,tau2,l)%FFC_D+&
                    &trialfc2_value(tau1,j)%phi*EXP(ci*(qp.dot.cell_vec(:,R2)))/&
                    &SQRT(iatom(tau1)%mass*iatom(tau2)%mass)*cell_vec(l,R2)
                END DO !l loop
            END DO !j loop
        END DO !i loop
!------------------------------------------------------------------------

        !assign ddyn from ddyn_temp, dynmat from dynmat_temp; to make them into arrays
        Do i=1,d*atom_number
            Do j=1,d*atom_number
                IF(MOD(i,d).eq.0) THEN
                    IF(MOD(j,d).eq.0) THEN
                        dynmat(i,j)=dynmat_temp(INT(i/d),INT(j/d))%FFC_D(d,d)
                        ddyn(i,j,:)=ddyn_temp(INT(i/d),INT(j/d),:)%FFC_D(d,d)
                    ELSE
                        dynmat(i,j)=dynmat_temp(INT(i/d),INT(j/d)+1)%FFC_D(d,MOD(j,d))
                        ddyn(i,j,:)=ddyn_temp(INT(i/d),INT(j/d)+1,:)%FFC_D(d,MOD(j,d))
                    END IF
                ELSE
                    IF(MOD(j,d).eq.0) THEN
                        dynmat(i,j)=dynmat_temp(INT(i/d)+1,INT(j/d))%FFC_D(MOD(i,d),d)
                        ddyn(i,j,:)=ddyn_temp(INT(i/d)+1,INT(j/d),:)%FFC_D(MOD(i,d),d)
                    ELSE
                        dynmat(i,j)=dynmat_temp(INT(i/d)+1,INT(j/d)+1)%FFC_D(MOD(i,d),MOD(j,d))
                        ddyn(i,j,:)=ddyn_temp(INT(i/d)+1,INT(j/d)+1,:)%FFC_D(MOD(i,d),MOD(j,d))
                    END IF
                END IF
            End Do !j loop
        End Do !i loop

        !fix rounding error or machine error
!        limit=1E-10
!        DO i=1,d*atom_number
!            DO j=1,d*atom_number
!                IF((ABS(dynmat(i,j)).lt.limit)) THEN
!                    dynmat(i,j)=CMPLX(0d0,0d0)
!                END IF
!            END DO
!        END DO


        DEALLOCATE(dynmat_temp)
        DEALLOCATE(ddyn_temp)

    END SUBROUTINE calculate_dynmat
!==================================================================================================
SUBROUTINE get_vgr(qp,ndn,vg,eival,eivec) 
!!get the group velocity/eival/eivec for a single k point
!!for group velocity method
!! this is a refined and crucial version of original <get_freq>
!! use this otherwise it will have negative w
!! also use this because the dynmat comes from my trialfc2 which is updated every Broyden loop
    USE Fourier_force_constants
    IMPLICIT NONE

    INTEGER :: ier,nd2,al
    INTEGER :: i,j,k,l,n,nv
    INTEGER,INTENT(IN) :: ndn
    INTEGER, ALLOCATABLE :: mp(:)
    REAL(8),INTENT(IN) :: qp(d)
    REAL(8), ALLOCATABLE :: eivl(:)
    REAL(8),INTENT(OUT) :: vg(d,ndn),eival(ndn)
    COMPLEX(8),ALLOCATABLE :: eivc(:,:)
    COMPLEX(8),INTENT(OUT) :: eivec(ndn,ndn)

!since we are gonna use finite difference method, the import part is to get eigenvalues
    !use module zhegv to diagonalize trialffc2_matrix for this k, and get eivecs(:,:,k), eivals(:,k)
            CALL calculate_dynmat(qp)
            ! ----------Born Charge----------
            IF (.NOT. ALL(qp.eq.kp_gamma)) THEN
                CALL nonanal(qp,dynmat,ndyn,ddyn)
            END IF
            ! -------------------------------
            n=d*atom_number
            nv=d*atom_number
            ier=0
            !CALL diagonalize(n,dynmat,eival,nv,eivec,ier)
!-------revision Feb.12--------
ALLOCATE(eivl(ndn),eivc(ndn,ndn),mp(ndn))
CALL diagonalize(n,dynmat,eivl,nv,eivc,ier)
CALL sort(ndn,eivl,mp,ndn) !crucial
do j=1,ndn
   eival(j) = cnst*sqrt(abs(eivl(mp(j)))) ! so that all frequencies are positive
!      eival(j) = abs(eivl(mp(j))) ! so that all frequencies are positive
   do l=1,ndn
      eivec(l,j) = eivc(l,mp(j))
   enddo
enddo
!------------------------------

            do j=1,ndn
                if (eival(j).lt.-1d-6) THEN
                    write(*,*)' Eigenvalue ',j,' is negative=',eival(j)
                    STOP
                ELSEIF(eival(j).lt.1d-18) then
                    eival(j)=1d-18
                end if
            end do
    ! now calculate the group velocities from Hellman-Feynmann formula(for each component al=1,3)
    vg=0
    do al=1,d
    do l=1,ndn
      vg(al,l)=0 !initialize the velocity
      do k=1,ndn
         dynmat(k,l)=0 !initialize the dynamic matrix?
         do j=1,ndn
            dynmat(k,l)=dynmat(k,l)+ddyn(k,j,al)*eivec(j,l)
         enddo
      enddo
      do k=1,ndn
         vg(al,l)=vg(al,l)+dynmat(k,l)*conjg(eivec(k,l))
      enddo
      vg(al,l)=vg(al,l)/2/sqrt(abs(eival(l)))*cnst*1d-10*100*2*pi
    enddo
    enddo

END SUBROUTINE get_vgr
!---------------------------------------------------------------------------
SUBROUTINE finitedif_vgr(q0,ndn,vgr,evl0,evc0)
!! this is a refined and crucial version of original <finitedif_vel>
!! use this otherwise it will have segment fault
    !calculates the group velocities from finite difference
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ndn
    REAL(8),INTENT(IN) :: q0(d)
    REAL(8),INTENT(OUT) :: vgr(d,ndn),evl0(ndn)
    COMPLEX(8),INTENT(OUT) :: evc0(ndn,ndn)

    INTEGER :: i,j
    REAL(8) :: q1(d),dq,om0,om1,aux(d,ndn),aux0(d,ndn)
    REAL(8) :: evlp(ndn),evlm(ndn)
    COMPLEX(8) :: evct(ndn,ndn)


    dq=3d-4  ! this is the ideal dq (at least for si)

    CALL get_vgr(q0,ndn,aux0,evl0,evc0)

    DO i=1,d
        q1=q0;q1(i)=q0(i)+dq
        CALL get_vgr(q1,ndn,aux,evlp,evct)
        q1=q0;q1(i)=q0(i)-dq
        CALL get_vgr(q1,ndn,aux,evlm,evct)
        vgr(i,:)=(SQRT(ABS(evlp))-SQRT(ABS(evlm)))/2/dq
    END DO

END SUBROUTINE finitedif_vgr
!=====================================================================================================================================
subroutine get_frequencies(nkp,kp,dk,ndn,eival,nv,eivec,vg)
!! also outputs the group velocities from HF theorem (exact) in units of c_light
!! eival is w not w^2, it's in unit cm^(-1)
 implicit none
 integer, intent(in) :: nkp,ndn,nv ! no of wanted eivecs
 real(8), intent(in) :: kp(3,nkp),dk(nkp)
 real(8), intent(out):: eival(ndn,nkp),vg(3,ndn,nkp)
 complex(8), intent(out) :: eivec(ndn,ndn,nkp)
 integer i,j,k,l,ier,nd2,al,ll
 integer, allocatable :: mp(:) !Map for ascending band sorting
 real(8), allocatable :: eivl(:)
 real(8) absvec,om
 complex(8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:)
 real(8) mysqrt
 character ext*3

 integer, allocatable :: mp1(:,:) !Map for projection band sorting
 real(8), allocatable :: eival_tmp(:) !Temp matrices for reordering eigenvalues
 complex(8), allocatable :: eivec_tmp(:,:)
 real(8), allocatable :: vg_tmp(:,:)

! open files and write the group velocities
! do l=1,ndn
!    write(ext,'(i3.3)')l
!    open(uvel+l,file='veloc-'//ext//'.dat')
!    write(uvel+l,*)'# la,i1,i2,i3,nk,kp(nk),eival(l,mapibz(nk)),vg(l,nk),length(vg(l,nk))'
! enddo
 nd2 = min(ndn,12)
!allocate(dynmat(ndn,ndn),mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3))
 allocate(mp(ndn),eivl(ndn),eivc(ndn,ndn))
! write(uio,*)'# dynmat allocated, printing: Kpoint, eival(j),j=1,ndyn)'

 kloop: do i=1,nkp

! write(uio,*)'############ before setting up dynmat, kp number ',i
    !call finitedif_vel(kp(:,i),ndn,vg(:,:,i),eivl,eivc)
    !try my own version
    call finitedif_vgr(kp(:,i),ndn,vg(:,:,i),eivl,eivc)

    do j=1,ndn
!      eival(j,i) = sqrt(abs(eivl(j)))*cnst
       eival(j,i) = eivl(j)
       do l=1,nv
          eivec(j,l,i) = eivc(j,l)  !mp(l))
       enddo
    enddo


 enddo kloop

! Band structure projection sorting, added by Bolin
!   allocate(mp1(ndn,nkp))
!    write(*,*) "Projection Band Sorting for band structures!"
!   call band_sort_bs(nkp,ndn,kp,dk,eival,eivec,mp1)
!   allocate(eival_tmp(ndn),eivec_tmp(ndn,ndn),vg_tmp(3,ndn))
!   do i=1,nkp
!      eival_tmp=eival(:,i)
!      eivec_tmp=eivec(:,:,i)
!      vg_tmp=vg(:,:,i)
!      do j=1,ndn
!         eival(j,i)=eival_tmp(mp1(j,i))
!         eivec(:,j,i)=eivec_tmp(:,mp1(j,i))
!         vg(:,j,i)=vg_tmp(:,mp1(j,i))
!      end do
!      call write_eigenvalues(i,dk(i),kp(:,i),eival(:,i),ndn,uio)
!   end do
!   deallocate(mp1,eival_tmp,eivec_tmp,vg_tmp)

 deallocate(eivl,eivc,mp)  !,dynmat,ddyn)

 3 format(a,i5,9(1x,f19.9))
 4 format(99(1x,2(1x,g9.3),1x))
 5 format(a,i5,99(1x,2(1x,g9.3),1x))
 6 format(2i6,2x,99(1x,f9.3))
 7 format(i6,2x,3(1x,f7.3),1x,99(1x,f7.3))

 end subroutine get_frequencies
!======================================================================================================================================
    SUBROUTINE check_conj_eivec
        !!some check subroutine, not used
        IMPLICIT NONE
        INTEGER :: temp,l,k
        REAL(8) :: real_part, imaginary_part
        REAL(8) :: threshold

        threshold = 1d-8

        OPEN(70,FILE='check_conj_eivec.dat',STATUS='unknown',ACTION='WRITE')
        WRITE(70,*) '===== tau_alpha, lambda, k, eivec ====='

        DO k=1,SIZE(kvector)
        DO temp=1,atom_number*d
        DO l=1,atom_number*d
            real_part = ABS(REAL(eivecs(temp,l,k)) - REAL(eivecs_t(l,temp,k)))
            imaginary_part = ABS(AIMAG(eivecs(temp,l,k))+AIMAG(eivecs_t(l,temp,k)))

            IF(real_part.gt.threshold) THEN
                WRITE(70,*)'real part wrong:'
                WRITE(70,9) temp,l,k,eivecs(temp,l,k),eivecs_t(temp,l,k)
                WRITE(70,*)
            ELSEIF(imaginary_part.gt.threshold) THEN
                WRITE(70,*)'imaginary part wrong:'
                WRITE(70,9) temp,l,k,eivecs(temp,l,k),eivecs_t(temp,l,k)
                WRITE(70,*)
            END IF
        END DO
        END DO
        END DO
9 FORMAT(3I3,2(G16.7,SP,G16.7,"i"))
        CLOSE(70)
    END SUBROUTINE check_conj_eivec
    
    FUNCTION Get_phi_test(k_number,direction1,atom1,direction2,atom2)
    !!some check subroutine, not used
        IMPLICIT NONE
        INTEGER :: i,j,k,l,steps
        INTEGER :: R1,tau1,R2,tau2,temp1,temp2,temp3,temp4
        INTEGER, INTENT(IN) :: k_number,atom1,direction1,atom2,direction2

        REAL(8) :: qp(3),R2_vec(3)


        COMPLEX(8) :: evc0(d*atom_number,d*atom_number),Get_phi_test

        INTEGER :: cutoff,cutoff_true
        LOGICAL :: auxiliary,FLAG

        !retrieve atomic indexes
        R1=every_atom(atom1)%type_R
        tau1=every_atom(atom1)%type_tau
        R2=every_atom(atom2)%type_R
        tau2=every_atom(atom2)%type_tau
        R2_vec(:)=every_atom(atom2)%R

        !re-get needed indexes
        temp1 = d*(tau1-1)+direction1
        temp2 = d*(tau2-1)+direction2


        Get_phi_test = CMPLX(0,0)
        DO l=1,SIZE(eivals,DIM=1)
        DO k=1,SIZE(kvector)
            qp=kvector(k)%component(:)
            Get_phi_test=Get_phi_test+&
            &SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/k_number*&
!            &dynmat_record(temp1,temp2,k)*EXP(-ci*(qp.dot.cell_vec(:,R2)))
            &eivals(l,k)*(eivecs(temp1,l,k)*eivecs_t(l,temp2,k))*EXP(-ci*(qp.dot.cell_vec(:,R2)))

        END DO !k loop
        END DO !l loop

    END FUNCTION Get_phi_test
!--------------------------------------------------------------------------------------------
    FUNCTION Get_Y_square(k_number,direction1,atom1,direction2,atom2)
    !!Major subroutine that calculate <YY> for specific indexes given
    !!utilize analytical approximation for diverging terms
        IMPLICIT NONE
        INTEGER :: i,j,k,l,steps
        INTEGER :: R1,tau1,R2,tau2,temp1,temp2,temp3,temp4
        INTEGER, INTENT(IN) :: k_number,atom1,direction1,atom2,direction2
        INTEGER :: unitnumber

        REAL(8) :: om_max,term,nbe,check,limit,vgr(d,d*atom_number),evl0(d*atom_number),delta_k(d)
        REAL(8),DIMENSION(:),ALLOCATABLE :: om,array!parameter omega; Y_square only sums lambda but not q
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: arg,func,res
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: n_qlambda !Bose-Einstein distribution
        REAL(8) :: coefficientA(6)
        REAL(8) :: delta_cubic,denominator(3),term1,term2
!        REAL(8) :: Get_Y_square
        COMPLEX(8) :: Get_Y_square

        COMPLEX(8) :: evc0(d*atom_number,d*atom_number)

        INTEGER :: cutoff,cutoff_true
        LOGICAL :: auxiliary,FLAG
        LOGICAL :: match

        !retrieve atomic indexes
        R1=every_atom(atom1)%type_R
        tau1=every_atom(atom1)%type_tau
        R2=every_atom(atom2)%type_R
        tau2=every_atom(atom2)%type_tau

!        ALLOCATE(func(SIZE(eivals,DIM=1),k_number),array(k_number))
        ALLOCATE(n_qlambda(SIZE(eivals,DIM=1),SIZE(eivals,DIM=2)))
!******************************BELOW IS TETRAHEDRON METHOD*****************************
!        array=0d0
!        !get needed indexes
!        IF(direction1.eq.d) THEN
!            temp1=d*tau1
!        ELSE
!            temp1=d*(tau1-1)+direction1
!        END IF
!
!        IF(direction2.eq.d) THEN
!            temp2=d*tau2
!        ELSE
!            temp2=d*(tau2-1)+direction2
!        END IF
!
!        !initialize func(k)
!        func=0
!        DO k=1,k_number ! k point loop
!        DO l=1,3
!
!            !---------------------------CLASSICAL VERSION------------------------------
!            !func(l,k)=func(l,k)+&
!            !    &temperature/eivals(l,k)/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/&
!            !    &(ee*1d20)*100*h_plank*c_light*&
!            !    &REAL(eivecs_t(l,temp1,k)*eivecs(temp2,l,k))*COS(kvector(k).dot.cell_vec(:,R2))
!            !IMPORTANT: updated k dot updated R
!            !--------------------------------------------------------------------------
!
!            !---------------------------QUANTUM VERSION--------------------------------
!            !IF(ABS(eivals(l,k)).lt.1E-10) THEN
!            !    eivals(l,k)=ABS(eivals(l,k))
!            !END IF
!
!            n_qlambda(l,k)=nbe(SQRT(eivals(l,k))*cnst,temperature,0) !0 for quantum ver.
!
!            !func is in S.I units
!            func(l,k)=func(l,k)+&
!                &hbar/2/SQRT(eivals(l,k))/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/SQRT(ee*1d20*uma)*&
!                &(2*n_qlambda(l,k)+1)*&
!                & REAL(eivecs_t(l,temp1,k)*eivecs(temp2,l,k))*COS(kvector(k).dot.cell_vec(:,R2))*1d20 !angstrom^2
!
!        END DO !l loop for eigen index sum
!        END DO !k loop for every k point

        !initialize om(steps)
!        steps=500
!        ALLOCATE(om(steps))
!        om_max=MAXVAL(MAXVAL(eivals,DIM=2)) !get the largest omega
!        om_max=SQRT(om_max)
!
!        DO i=1,steps
!            om(i)=om_max/steps*i
!        END DO

        !prepare arg(k)
!        ALLOCATE(arg(SIZE(eivals,DIM=1),SIZE(eivals,DIM=2)))
!        DO k=1,k_number
!            DO l=1,3
!                IF(ABS(eivals(l,k)).lt.1E-10) THEN
!                    arg(l,k)=0d0
!                ELSE
!                    arg(l,k)=SQRT(SQRT(eivals(l,k)))
!                    !arg(l,k)=SQRT(eivals(l,k))
!                    !arg(l,k)=eivals(l,k)
!                END IF
!            END DO ! l loop
!        END DO ! k loop



!-------------------------find the minimum trust om(i) by group velocity. PREVIOUSLY------------------------------------
!        temp1=tet(7)%p(1)%i;temp2=tet(7)%p(2)%i
!        check=ABS(kvector(temp2)%component(1)-kvector(temp1)%component(1))
!        delta_k=(/0d0,SQRT(5d0)*check,0d0/)
!        CALL finitedif_vgr(delta_k,d*atom_number,vgr,evl0,evc0) !temp1-th corner
!        cutoff=INT(SQRT(ABS(vgr(:,3)).dot.delta_k)*steps/om_max)
!        WRITE(*,*) 'Group Velocity: the critical om is at:',cutoff
!
!        !get every F(om) after the cutoff for the 3 acoustic bands
!        ALLOCATE(res(steps,SIZE(eivals,DIM=1)))
!        DO l=1,3
!        DO i=cutoff-1,steps
!            CALL tet_sum(om(i),k_number,arg(l,:),func(l,:),res(i,l),array)
!        END DO !i loop
!        END DO !l loop
!***********************check F(omega)***************************
!IF(tau1.eq.1 .AND. direction1.eq.1 .AND. atom2.eq.1 .AND. direction2.eq.1) THEN
!OPEN(44,FILE='F_om.dat',STATUS='unknown',ACTION='write')
!check=0d0
!DO i=1,steps
!    WRITE(44,*)om(i),SUM(res(i,1:6)) !plot the F(om) for 3 acoustic bands
    !check=check+SUM(res(i,:))*om_max/steps
!END DO
!CLOSE(44)
!WRITE(*,*)'Integral of DOS=',check
!WRITE(*,*)'Maximum Frequency=',om_max*cnst
!ENDIF
!*****************************************************************
!---------------------------find the minimum trust om(i) from the graph--------------------------------------------------
!        limit=(om_max-danger)/steps
!        array=0
!        cutoff_true=1
!        DO i=cutoff,steps-1
!            array(i)=ABS(SUM(res(i+1,:))-SUM(res(i,:)))
!            IF(array(i).lt.limit) THEN
!                cutoff_true=i+1
!                !WRITE(*,*) 'From Graph: the critical om:',om(cutoff_true)
!                WRITE(*,*)'From Graph: the critical om is at:',cutoff_true
!                EXIT
!            END IF
!        END DO
!--------------------------fix the found cutoff for the rest of the iterations--------------------------------------------
!        cutoff=cutoff_true

        !tweak the diverging part
!        DO l=1,3
!            CALL tet_sum(om(cutoff),k_number,arg(l,:),func(l,:),res(cutoff,l),array)
!            coefficientA(l)=res(cutoff,l)/om(cutoff)**3 *TANH(om(cutoff)**2*cnst/2/temperature)
!        END DO
!
!        DO i=1,cutoff
!        DO l=1,3
!            res(i,l)=coefficientA(l)*om(i)**3 / TANH(om(i)**2*cnst/2/temperature)
!        END DO !l loop
!        END DO!i loop


!-----------------------------SUM UP THE ACOUSTIC MODES--------------------------------
!        term=0d0
!        DO l=1,3
!        DO i=1,steps
!            term=term+res(i,l)
!        END DO !i loop
!        END DO !l loop
!
!        Get_Y_square=term/steps*om_max*k_number !k_number=nk, which will be dividing outside this subroutine
!*******************ABOVE IS TETRAHEDRON METHOD****************************

!*******************CURRENT METHOD****************************
        !re-get needed indexes
        temp1 = d*(tau1-1)+direction1
        temp2 = d*(tau2-1)+direction2

        !initiate Get_Y_square
        Get_Y_square=CMPLX(0d0,0d0)

        !Acoustic bands
        DO l=1,3
        DO k=2,SIZE(kvector) !drop gamma point region

!-------------------------------------------NEW----------------------------------------------------
            !skip if they were soft
            match = .FALSE.
            IF(min_eival.le.0d0) THEN
           checksoft_a: DO i=1,SIZE(soft)
                IF(l.eq.soft(i)%idx_la .AND. k.eq.soft(i)%idx_q) THEN
                    match = .TRUE.
                    EXIT checksoft_a
                END IF
            END DO checksoft_a
            END IF

            IF(match) CYCLE
!--------------------------------------------------------------------------------------------------

            !****DO THIS IF CHOOSE option.1 or option.2****
            n_qlambda(l,k)=nbe(SQRT(eivals(l,k))*cnst,temperature,0)
            Get_Y_square=Get_Y_square+&
                &hbar/2/SQRT(eivals(l,k))/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/SQRT(ee*1d20*uma)*&
                &(2*n_qlambda(l,k)+1)*&
                &eivecs(temp1,l,k)*eivecs_t(l,temp2,k)*EXP(-ci*(kvector(k).dot.cell_vec(:,R2)))*1d20

            !****ONLY DO THIS IF CHOOSE option.3 WHEN DEAL WITH NEGATIVE W****
!            IF(eivals(l,k).gt.0d0) THEN
!                n_qlambda(l,k)=nbe(SQRT(eivals(l,k))*cnst,temperature,0)
!
!                Get_Y_square=Get_Y_square+&
!                &hbar/2/SQRT(eivals(l,k))/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/SQRT(ee*1d20*uma)*&
!                &(2*n_qlambda(l,k)+1)*&
!                &eivecs(temp1,l,k)*eivecs_t(l,temp2,k)*EXP(-ci*(kvector(k).dot.cell_vec(:,R2)))*1d20
!            ELSE
!                n_qlambda(l,k)=-nbe(SQRT(-eivals(l,k))*cnst,temperature,1) !keep it as it is so the negative sign in the front
!
!                Get_Y_square=Get_Y_square+&
!                &hbar/2/SQRT(-eivals(l,k))/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/SQRT(ee*1d20*uma)*&
!                &(2*n_qlambda(l,k)+1)*&
!                &eivecs(temp1,l,k)*eivecs_t(l,temp2,k)*EXP(-ci*(kvector(k).dot.cell_vec(:,R2)))*1d20
!            ENDIF




IF(isnan(REAL(Get_Y_square))) THEN
WRITE(34,*)'--------------POSITION1---------------'
WRITE(34,*)'Get_Y_square', Get_Y_square
WRITE(34,*)'lambda',l,'#kvector',k
WRITE(34,*)'eivals',eivals(l,k)
WRITE(34,*)'n_qlambda',n_qlambda(l,k) !HERE COMES THE PROBLEM
WRITE(34,*)'eivecs',eivecs(temp2,l,k),eivecs_t(l,temp1,k)
WRITE(34,*)'cell_vec',cell_vec(:,R2)
STOP
END IF
        END DO !k loop
        END DO !l loop

        !correction term for 3 acoustic bands at Gamma point
        denominator(1) = (eivals(3,1)+eivals(3,3)-2*eivals(3,2))
        denominator(2) = (eivals(2,1)+eivals(2,3)-2*eivals(2,2))
        denominator(3) = (eivals(1,1)+eivals(1,3)-2*eivals(1,2))

        delta_cubic = 3d0*volume_g/(4*pi*k_number)
        term1 = delta_cubic*(4d0/3/denominator(1)+2d0/3/denominator(2)+2d0/denominator(3))*1d-20*uma/ee
        term2 = delta_cubic*(2d0/3/denominator(1)+4d0/3/denominator(2))*1d-20*uma/ee

        IF(direction1.eq.direction2) THEN
            IF(direction1.eq.1 .OR. direction1.eq.2) THEN
                Get_Y_square = Get_Y_square + 2*pi*temperature*100*h_plank*c_light&
                        &/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/uma*term1*1d20
            ELSE
                Get_Y_square = Get_Y_square + 4*pi*temperature*100*h_plank*c_light&
                        &/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/uma*term2*1d20
            END IF
        END IF !xx,yy,zz should be equal

!-------------------------------------------NEW----------------------------------------------------
        !correction term for those soft bands that get dropped at certain q
        IF(min_eival.le.0d0) THEN !first make sure there IS any negative mode
        DO i=1,SIZE(soft)
            l=soft(i)%idx_la
            k=soft(i)%idx_q

            IF(k.eq.1) CYCLE !already fixed gamma point above

            IF(l.le.3) THEN !acoustic?
                IF(k.le. SIZE(eivals,DIM=2)-2) THEN
                    denominator(1) = (eivals(3,k)+eivals(3,k+2)-2*eivals(3,k+1))
                    denominator(2) = (eivals(2,k)+eivals(2,k+2)-2*eivals(2,k+1))
                    denominator(3) = (eivals(1,k)+eivals(1,k+2)-2*eivals(1,k+1))
                ELSE
                    denominator(1) = (eivals(3,k-2)+eivals(3,k)-2*eivals(3,k-1))
                    denominator(2) = (eivals(2,k-2)+eivals(2,k)-2*eivals(2,k-1))
                    denominator(3) = (eivals(1,k-2)+eivals(1,k)-2*eivals(1,k-1))
                END IF
            ELSE
                l=INT((l-1)/3)*3
                IF(k.le. SIZE(eivals,DIM=2)-2) THEN
                    denominator(1) = (eivals(l+3,k)+eivals(l+3,k+2)-2*eivals(l+3,k+1))
                    denominator(2) = (eivals(l+2,k)+eivals(l+2,k+2)-2*eivals(l+2,k+1))
                    denominator(3) = (eivals(l+1,k)+eivals(l+1,k+2)-2*eivals(l+1,k+1))
                ELSE
                    denominator(1) = (eivals(l+3,k-2)+eivals(l+3,k)-2*eivals(l+3,k-1))
                    denominator(2) = (eivals(l+2,k-2)+eivals(l+2,k)-2*eivals(l+2,k-1))
                    denominator(3) = (eivals(l+1,k-2)+eivals(l+1,k)-2*eivals(l+1,k-1))
                END IF
            END IF

            delta_cubic = 3d0*volume_g/(4*pi*k_number)
            term1 = delta_cubic*(4d0/3/denominator(1)+2d0/3/denominator(2)+2d0/denominator(3))*1d-20*uma/ee
            term2 = delta_cubic*(2d0/3/denominator(1)+4d0/3/denominator(2))*1d-20*uma/ee

            IF(direction1.eq.direction2) THEN
                IF(direction1.eq.1 .OR. direction1.eq.2) THEN
                    Get_Y_square = Get_Y_square + 2*pi*temperature*100*h_plank*c_light&
                            &/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/uma*term1*1d20
                ELSE
                    Get_Y_square = Get_Y_square + 4*pi*temperature*100*h_plank*c_light&
                            &/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/uma*term2*1d20
                END IF
            END IF !xx,yy,zz should be equal

        END DO !soft loop
        END IF
!--------------------------------------------------------------------------------------------------
!IF(isnan(Get_Y_square)) THEN
!WRITE(34,*)'--------------POSITION2---------------'
!WRITE(34,*)'term',term
!STOP
!END IF
!***** check the contribution of this gamma point correction term *****
for_check = (0d0,0d0)

IF(direction1.eq.direction2) THEN
    for_check = 2*pi*temperature*100*h_plank*c_light&
                    &/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/uma*term1*1d20
ELSE
    for_check = 4*pi*temperature*100*h_plank*c_light&
                    &/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/uma*term2*1d20
END IF
!**********************************************************************


!WRITE(*,*) 'Extra Term=:', 2*pi*temperature*100*h_plank*c_light&
!                               &/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/uma*term*1d20

        !Optic bands
        DO l=4,SIZE(eivals,DIM=1)
        DO k=1,SIZE(kvector)

!-------------------------------------------NEW----------------------------------------------------
            !skip if they were soft
            match = .FALSE.
            IF(min_eival.le.0d0) THEN
           checksoft_o: DO i=1,SIZE(soft)
                IF(l.eq.soft(i)%idx_la .AND. k.eq.soft(i)%idx_q) THEN
                    match = .TRUE.
                    EXIT checksoft_o
                END IF
            END DO checksoft_o
            END IF

            IF(match) CYCLE
!--------------------------------------------------------------------------------------------------


            !****DO THIS IF CHOOSE option.1 or option.2****
            n_qlambda(l,k)=nbe(SQRT(eivals(l,k))*cnst,temperature,0)
            Get_Y_square=Get_Y_square+&
                 &hbar/2/SQRT(eivals(l,k))/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/SQRT(ee*1d20*uma)*&
                 &(2*n_qlambda(l,k)+1)*&
                 &eivecs(temp1,l,k)*eivecs_t(l,temp2,k)*EXP(-ci*(kvector(k).dot.cell_vec(:,R2)))*1d20

            !****ONLY DO THIS IF CHOOSE option.3 WHEN DEAL WITH NEGATIVE W****
!            IF(eivals(l,k).gt.0d0) THEN
!                n_qlambda(l,k)=nbe(SQRT(eivals(l,k))*cnst,temperature,0)
!
!                Get_Y_square=Get_Y_square+&
!                 &hbar/2/SQRT(eivals(l,k))/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/SQRT(ee*1d20*uma)*&
!                 &(2*n_qlambda(l,k)+1)*&
!                 &eivecs(temp1,l,k)*eivecs_t(l,temp2,k)*EXP(-ci*(kvector(k).dot.cell_vec(:,R2)))*1d20
!            ELSE
!                n_qlambda(l,k)=-nbe(SQRT(-eivals(l,k))*cnst,temperature,1) !keep it as it is so the negative sign in the front
!
!                Get_Y_square=Get_Y_square+&
!                 &hbar/2/SQRT(-eivals(l,k))/SQRT(iatom(tau1)%mass*iatom(tau2)%mass)/SQRT(ee*1d20*uma)*&
!                 &(2*n_qlambda(l,k)+1)*&
!                 &eivecs(temp1,l,k)*eivecs_t(l,temp2,k)*EXP(-ci*(kvector(k).dot.cell_vec(:,R2)))*1d20
!            ENDIF



!IF(isnan(Get_Y_square)) THEN
!WRITE(34,*)'--------------POSITION3---------------'
!WRITE(34,*)'lambda',l,'#kvector',k
!WRITE(34,*)'eivals',eivals(l,k)
!WRITE(34,*)'n_qlambda',n_qlambda(l,k)
!WRITE(34,*)'eivecs',eivecs(temp2,l,k),eivecs_t(l,temp1,k)
!WRITE(34,*)'cell_vec',cell_vec(:,R2)
!STOP
!END IF

        END DO !k loop
        END DO !l loop
        DEALLOCATE(n_qlambda)
!        DEALLOCATE(res,arg,om,func,array)
    END FUNCTION Get_Y_square
!======================================================================================================================================
    SUBROUTINE Get_SKfreq(kpoint,freq)
    !!get omega^2 at a single k-point
    !!only used for dispersion plot
        IMPLICIT NONE
        INTEGER :: ndim
        INTEGER :: atom1,atom2,R1,R2,tau1,tau2,direction1,direction2
        INTEGER :: i,j,l,n,nv,ier
        INTEGER, ALLOCATABLE :: mp(:)
        INTEGER, INTENT(IN) :: kpoint
        REAL(8),INTENT(OUT) :: freq(:)
        REAL(8),ALLOCATABLE :: freq_square(:)
        REAL(8) :: limit,min_eivals,max_eivals, qp(d)
        COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: junk
        COMPLEX(8) :: dummy(d,d)

        ndim = d*atom_number ! should be

        ALLOCATE(junk(d*atom_number,d*atom_number))
        ALLOCATE(skffc2_value(atom_number,atom_number),skffc2_matrix(d*atom_number,d*atom_number))


        DO tau1=1,atom_number
        DO tau2=1,atom_number
            skffc2_value(tau1,tau2)%FFC_D=0d0
        END DO !tau2 loop
        END DO !tau1 loop


        !calculate all the dynmat value at kp_bs(:,kpoint) for tau1 tau2
        DO i=1,atom_number
            tau1=i
        DO j=1,tot_atom_number
            tau2=every_atom(j)%type_tau
            R2=every_atom(j)%type_R

            dummy=trialfc2_value(tau1,j)%phi*EXP(ci*(kp_bs(:,kpoint).dot.cell_vec(:,R2)))
            skffc2_value(tau1,tau2)%FFC_D=skffc2_value(tau1,tau2)%FFC_D+dummy/SQRT(iatom(tau1)%mass * iatom(tau2)%mass)

        END DO !j loop
        END DO !i loop


        !we have every value for skffc2_value, now need to put it into matrix form: skffc2_matrix(d*atom_number,d*atom_number)
        !the matrix has to be an Hermitian
        Do i=1,d*atom_number
            Do j=1,i

                IF(MOD(i,d).eq.0) THEN
                    IF(MOD(j,d).eq.0) THEN
                        skffc2_matrix(i,j)=skffc2_value(INT(i/d),INT(j/d))%FFC_D(d,d)
                    ELSE
                        skffc2_matrix(i,j)=skffc2_value(INT(i/d),INT(j/d+1))%FFC_D(d,MOD(j,d))
                    END IF
                ELSE
                    IF(MOD(j,d).eq.0) THEN
                        skffc2_matrix(i,j)=skffc2_value(INT(i/d+1),INT(j/d))%FFC_D(MOD(i,d),d)
                    ELSE
                        skffc2_matrix(i,j)=skffc2_value(INT(i/d+1),INT(j/d+1))%FFC_D(MOD(i,d),MOD(j,d))
                    END IF
                END IF

            End Do !j loop
        End Do!i loop

        !*********************add born correction to every skffc2_value*********************
        !notes: calculate dynmat & ddyn by <calculate_dynmat> then call <nonanal> to
        !add born correction, then let skffc2_matrix(i,j) equal to new dynmat(i,j)
        qp = kp_bs(:,kpoint)
        ndim = SIZE(dynmat, DIM=1)
        CALL calculate_dynmat(qp)
        ! problem at gamma point
        IF (kpoint.ne.1) THEN !needs revision, replace this with qp.ne.gamma point for 2nd visit
            CALL nonanal(qp,dynmat,ndim,ddyn)
        END IF
        DO i=1,ndim
            DO j=1, i
                skffc2_matrix(i,j)=dynmat(i,j)
            END DO !j loop
        END DO !i loop
        !***********************************************************************************


        !Force Hermitian: fill  in the rest with corresponding conjg
        DO i=1,d*atom_number-1
            DO j=i+1,d*atom_number
                skffc2_matrix(i,j)=CONJG(skffc2_matrix(j,i))
            END DO
        END DO


        !fix rounding error or machine error
            !limit=danger
            !DO i=1,d*atom_number
            !    DO j=1,d*atom_number
            !        IF((ABS(skffc2_matrix(i,j)).lt.limit)) THEN
            !            skffc2_matrix(i,j)=0
            !        END IF
            !    END DO
            !END DO

        !use module zhegv to diagonalize skffc2_matrix for this k, and get omega^2
            n=d*atom_number
            nv=d*atom_number
            ier=0
            ALLOCATE(freq_square(n))
            CALL diagonalize(n,skffc2_matrix,freq_square(:),nv,junk(:,:),ier)

        DEALLOCATE(skffc2_value,skffc2_matrix)
        DEALLOCATE(junk)


        !**********Fix if there is negative eigenvalues***********?????

        !min_eivals=0;max_eivals=MAXVAL(MAXVAL(eivals,DIM=2))

        !DO l=1,n
        !    IF(freq_square(l).lt.min_eivals) THEN
        !        min_eivals=freq_square(l)
        !        WRITE(*,*) 'Negative Eigenvalue in dispersion',freq_square(l)
        !    END IF
        !END DO

        !DO l=1,n
        !    freq_square(l)=(freq_square(l)-min_eivals)/(max_eivals-min_eivals)*max_eivals
        !END DO

!        DO l=1,n
!            IF(freq_square(l).lt.0d0) THEN
!                WRITE(*,*) 'Negative Eigenvalue in dispersion',freq_square(l)
!            END IF
!            IF(freq_square(l).le.0d0 .AND. freq_square(l).gt.-danger) THEN
!                freq_square(l)=freq_square(l)+danger
!            ELSEIF(freq_square(l).le.-danger) THEN
!                WRITE(*,*) 'Too Negative eigenvalues in dispersion:',freq_square(l)
!                STOP
!            ELSE
!                CYCLE
!            END IF
!        END DO

        ALLOCATE(mp(ndim))
        CALL sort(ndim,freq_square, mp,ndim) !crucial
        DO j=1,ndim
           freq(j) = cnst*sqrt(abs(freq_square(mp(j))))
        END DO
        !*********************************************************
    END SUBROUTINE Get_SKfreq
!----------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE Get_AKfreq(kp_bs,SE)
    !!get omega^2 at every k-point in kp_bs(:,k), the k path
    !!only used for dispersion plot
        IMPLICIT NONE
        INTEGER :: i,j,unit_number
        REAL(8),INTENT(IN),DIMENSION(:,:) :: kp_bs
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: freq
        LOGICAL,INTENT(in) :: SE

        ALLOCATE(freq(d*atom_number,SIZE(kp_bs,DIM=2)))

        unit_number=78
        OPEN(unit_number,FILE='Dispersion.dat',STATUS='Unknown',ACTION='write')

        DO i=1,SIZE(kp_bs,DIM=2)
            CALL Get_SKfreq(i,freq(:,i))
!            DO j=1,SIZE(freq,DIM=1)
!                IF (ABS(freq(j,i)).lt.1E-10) THEN
!                    freq_square(j,i)=0d0
!                END IF ! very small negative eigenvalue
!            END DO !j loop

            IF (SE) THEN
                freq(:,i) = freq(:,i) + REAL(nself(:,i)) + REAL(uself(:,i))
            END IF
            WRITE(unit_number,*) i,freq(:,i)/33.356
        END DO !i loop
        DEALLOCATE(freq)
        CLOSE(unit_number)
    END SUBROUTINE Get_AKfreq
!==============================================================================================
!==============================================================================================
!=====================================POST PROCESS=============================================
!==============================================================================================
!==============================================================================================
!----------------------For Output Phonon Dispersion-------------------------
    SUBROUTINE Get_Dispersion
    !! for plot dispersion
        IMPLICIT NONE
        CALL make_kp_bs2
        CALL Get_AKfreq(kp_bs,.false.)
    END SUBROUTINE Get_Dispersion
    
    SUBROUTINE Get_Dispersion_SE
    !! for plot dispersion with self energy, not used and not finished
        IMPLICIT NONE
        INTEGER nkbs,unit_number1,unit_number2,unit_number3
        INTEGER i,j,ksubset(2)
        CHARACTER(40) :: file1, file2, file3
        REAL(8),ALLOCATABLE :: junk(:,:,:)

        ! get params
        ndyn = d*atom_number
!***********************step 1 get all om2 which are on the FBZ****************************
        ! om2 and ev2 are stored in global variables eigenvals & eigenvecs

        !allocate eivalibz(ndyn,nibz),eivecibz(ndyn,ndyn,nibz),velocibz(3,ndyn,nibz) and also dkc(nibz)
        nkc = SIZE(kpc, dim=2)
        CALL allocate_eig(ndyn,nibz,nkc)
        ALLOCATE(dkc(nkc))
        dkc=0
        CALL get_frequencies(nkc,kpc,dkc,ndyn,eigenval,ndyn,eigenvec,veloc)
        DEALLOCATE(dkc)  !,eivalibz,eivecibz)
        !-----------------check om2-------------
        unit_number2 = 218
        file2 = 'Basemesh_om2.dat'
        OPEN(unit_number2, FILE=TRIM(file2),STATUS='unknown')
        DO i=1, SIZE(eigenval, dim=2)
            WRITE(unit_number2,*) i, eigenval(:,i)
        END DO
        CLOSE(unit_number2)

!***********step 2: get all om1 which are on the selected path kp_bs(:,:) ***************
        ! om1 and ev1 are stored in global variables eivals & eivecs
        CALL make_kp_bs2
        nkbs = SIZE(kp_bs, dim = 2)

        ALLOCATE(dkc(nkbs))
        dkc = 0
        IF(ALLOCATED(eivals)) THEN !important, because they ARE allocated in VA
            DEALLOCATE(eivals,eivecs,eivecs_t)
        END IF
        CALL allocate_eigen(ndyn,nkbs)
        ALLOCATE(junk(d,ndyn,nkbs))
        CALL get_frequencies(nkbs,kp_bs,dkc,ndyn,eivals,ndyn,eivecs,junk)
        DEALLOCATE(junk)
        DEALLOCATE(dkc)
        !-----------------check om1-------------
        unit_number1 = 438
        file1 = 'Dispersion_om1.dat'
        OPEN(unit_number1, FILE=TRIM(file1),STATUS='unknown')
        DO i=1, SIZE(eivals, dim=2)
            WRITE(unit_number1,*) i, eivals(:,i)
        END DO
        CLOSE(unit_number1)
!*****************************************************************************************
        ! ksubset needs to be decided, cause nibz will change after make_kp_bs2
!        CALL read_ksubset(ksubset,nibz)
        ksubset = (/1,nkbs/)

        CALL mySelfEnergy(ksubset,temperature)
        ! after that nself(:,:) and uself(:,:) are collected

        !---------------dispersion with self-energy, no need for new way----------------
!        unit_number3 = 124
!        file3 = 'Dispersion_se.dat'
!        OPEN(unit_number3, FILE=TRIM(file3),STATUS='unknown')
!        DO i=1, SIZE(eivals, dim=2)
!            WRITE(unit_number3,*) i, eivals(:,i) + REAL(uself(:,i) + nself(:,i))
!        END DO
!        CLOSE(unit_number3)

    END SUBROUTINE Get_Dispersion_SE
!================================================================================================
subroutine write_lat_fc(ngrps,ntrms)
!!Force Constants I/O, legacy code not used
 use svd_stuff
 use ios
 use force_constants_module
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
 real(8) rij
 integer i,j,ngrps(4),ntrms(4),j_sc

 write(ufco,*)' Crystal data: translation vectors of the primitive cell '
 write(ufco,9)r01
 write(ufco,9)r02
 write(ufco,9)r03
 write(ufco,*)' Crystal data: atoms in primitive cell: label,type,x,y,z,mass '
 write(ufco,*)natoms0
 do i=1,natoms0
   write(ufco,6)i,atom0(i)%name, atom0(i)%at_type,atom0(i)%equilibrium_pos,atom0(i)%mass
 enddo
 write(ufco,*)' Crystal data written ************************************'
 write(ufco,*)' Included ranks of FCs '
 write(ufco,*)include_fc(1),include_fc(2),include_fc(3),include_fc(4)
 write(ufco,*)' Number of FCs for each rank '
 write(ufco,*)ntrms (1),ntrms(2),ntrms(3),ntrms(4)
 write(ufco,*)' Number of independent FCs for each rank '
 write(ufco,*)ngrps(1),ngrps(2),ngrps(3),ngrps(4)
 write(ufco,*)' Neighborshell atoms: i,x,y,z,type_tau,n1,n2,n3 '
 write(ufco,*)natoms
 do i=1,natoms
   rij = length(atompos(:,i)-atompos(:,1))
   write(ufco,7)i,(atompos(j,i),j=1,3), iatomcell0(i),(iatomcell(j,i),j=1,3),rij
 enddo

 open(173,file='primlatt.xyz')
 write(173,*) natoms0
 write(173,*) ngrps(1),ngrps(2),ngrps(3),ngrps(4)
 do i=1,natoms0
   write(173,8)atom0(i)%name,atom0(i)%equilibrium_pos
!  write(173,7)atom0(i)%name,(atompos(j,i),j=1,3)
 enddo
 close(173)

 open(173,file='latfc.xyz')
 write(173,*) natoms
 write(173,*) ngrps(1),ngrps(2),ngrps(3),ngrps(4)
 do i=1,natoms
      call findatom_sc(iatomcell(:,i),iatomcell0(i),j_sc)
   write(173,8)atom0(iatomcell0(i))%name,(atompos(j,i),j=1,3),j_sc
 enddo
 close(173)

6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(a,3(2x,f19.10),i6)
9 format(9(2x,f19.10))
 end subroutine write_lat_fc

subroutine write_output_fc2
 use svd_stuff
 use ios
 use force_constants_module
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
 integer rnk,t,ti,i,res,j,rs !ni(4),nt(4)
 integer iat(4),ixyz(4),g,ng,term,term2,cnt2,frm
 real(8) rij,bunit,one,fcd,trace,dij
 character frmt*2,goh*48,ln*1,geh*47

 bunit = ryd/ab/ab
 one =1d0

6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(2(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
9 format(9(2x,f19.10))

! first write the crystal data
! nt=0; ni=0
! do i=1,4
!    if(map(i)%ngr.gt.0) then
!      nt(i)= sum(map(i)%nt(:))
!      ni(i)= sum(map(i)%ntind(:))
!    endif
! enddo
! call write_lat_fc(ni,nt)
 call write_lat_fc(map(:)%ntotind,map(:)%ntot)

!----------------------------------------
 res = 0
 do rnk=1,4
  if ( include_fc(rnk) .ne. 0 ) then
    frm=30+rnk
    write(ln,'(i1)')rnk
    goh='(a1,i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,2(g14.8,2x),f5.2)'
    geh='(i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g14.8,f8.4,2x,f9.5)'
    write(frmt,'(i2)')30+rnk
    write(ulog,*)' FOR RANK=',rnk,' format=',frmt
    write(*,*)' FOR RANK=',rnk,' format=',frmt
    write(ufc1-1+rnk,*)'# RANK ',rnk,' tensors :term,group,(iatom,ixyz)_2 d^nU/dx_{i,alpha}^n'
!   call get_dim(rnk,map(rnk),ndindp(rnk),ndfull(rnk)) ! no need already called

  ng=map(rnk)%ngr ! number of groups
  cnt2=0
  term = 0
  term2= 0
  do g=1,map(rnk)%ngr  ! index of a given group
       if(g.gt.1) cnt2=cnt2+map(rnk)%ntind(g-1)

    ! write in the log and fcn_fit.dat: cnt2+ti is the position of indep_fc of that rank
    do ti=1,map(rnk)%ntind(g)  ! index of independent terms in that group g
       iat(1:rnk)  = map(rnk)%gr(g)%iatind (:,ti)
       ixyz(1:rnk) = map(rnk)%gr(g)%ixyzind(:,ti)
       term = term+1
       if (rnk.eq.2) then
         rij = length(atompos(:,iat(1))-atompos(:,iat(2)))
       else
         rij = 0
       endif
       write(ulog,goh) map(rnk)%err(cnt2+ti),g,ti,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+cnt2+ti),fcs(res+cnt2+ti)/ryd*ab**rnk,rij
       write(ufit1-1+rnk,geh) ti,g,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+cnt2+ti),one,rij
    enddo

    ! write in the fcn.dat file
    do t=1,map(rnk)%nt(g)  ! index of dependent terms in that group g
       iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
       ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
       term2= term2+1
       fcd = 0
       ! must find the corresponding index of the indep term t <-> ti
       do ti=1,map(rnk)%ntind(g)
          ! this is the index of the indep FC coming in the A*FC=b matrix product
          fcd = fcd + fcs(res+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
       enddo
       if( abs(fcd) .gt. margin) then
          write(ufc1-1+rnk,geh)t,g, (iat(j),ixyz(j),j=1,rnk),fcd,one
       endif
    enddo
  enddo
!  res = res+ndindp(rnk)
  res = res+map(rnk)%ntotind
  endif
 enddo

write(ulog,*)'******* Trace for the harmonic FCs ********'
 open(456,file='trace_fc.dat')

! write the trace of FC2
 rnk=2
  iloop: do i=1,natoms0
  jloop: do j=1,natoms
     rij = length(atompos(:,i)-atompos(:,j))
     trace=0
!     rs=ndindp(1)
     rs=map(1)%ntotind
     if ( include_fc(rnk) .ne. 0 ) then
        ng=map(rnk)%ngr ! number of groups
        cnt2=0
        term2= 0
  do g=1,map(rnk)%ngr  ! index of a given group
       if(g.gt.1) cnt2=cnt2+map(rnk)%ntind(g-1)
    do t=1,map(rnk)%nt(g)  ! index of independent terms in that group g
       iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
       ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
  !    i=iat(1) ; j=iat(2)
       if (iat(1).ne.i) cycle !iloop
       if (iat(2).ne.j) cycle ! jloop
!      write(*,*)'i,j,term2,al,be=',i,j,term2,ixyz(1),ixyz(2)
       fcd = 0
       ! must find the corresponding index of the indep term t <-> ti
       do ti=1,map(rnk)%ntind(g)
          ! this is the index of the indep FC coming in the A*FC=b matrix product
          fcd = fcd + fcs(rs+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
       enddo
       dij = length(atompos(:,iat(1))-atompos(:,iat(2)))
       if (ixyz(1).eq.ixyz(2)) then
          term2= term2+1
          trace = trace+ fcd
  !       write(*,*)'al,term2,trace=',ixyz(1),term2,trace
       endif
    enddo
  enddo
       if(trace.ne.0) then
           write(ulog,8) i,j,dij,trace
           write(456,8) i,j,dij,trace
       endif
 !     if (term2.ne.3) write(456,*)'#ERROR: there are ',term2,' terms for rij=',dij

     endif
  enddo jloop
  enddo iloop

 close(456)

write(ulog,*)'***************** END OF FC Trace ******************'


!----------------------------------------
 125  format(a)

 if (res.ne.ngr) then
    write(ulog,*)'WRITE_OUTPUT: sum(nterms),ngr=',res,ngr
    write(ulog,*)'WRITE_OUTPUT: they should be equal!'
 endif
 end subroutine write_output_fc2

!!legacy code, maybe used somewhere
subroutine read_fcs_2(iunit,rank)
 use svd_stuff
 use ios
 use force_constants_module
 use atoms_force_constants
 use params
 implicit none
 integer rank,iunit,t,res,i
 character line*99

 read(iunit,'(a)') line
 t = 0
if ( rank .eq. 1) then

 res = 0
 do i=1,nterms(rank)
       read(iunit,*,err=91)t,igroup_1(t), &
& iatomterm_1(1,t),ixyzterm_1(1,t),  &
& fcs_1(igroup_1(t)),ampterm_1(t)
 enddo
 return
91 write(ulog,*)'READ_FCS: rank=',rank,' reached end of file after t=',t
 if (t.eq.0) then
    nterms(1) = 1
    iatomterm_1(1:1,1:1) = 1
    ixyzterm_1 (1:1,1:1) = 1
    igroup_1 (1:1) = 0
    ampterm_1(1:1) = 0
    include_fc(rank) = 0       ! if nothing, then exclude from fitting
 endif
!----------------------------------------

elseif ( rank .eq. 2) then

 res =  igroup_1(nterms(1))
 do i=1,nterms(rank)
       read(iunit,*,err=92)t,igroup_2(t), &
& iatomterm_2(1,t),ixyzterm_2(1,t),iatomterm_2(2,t),ixyzterm_2(2,t),  &
& fcs_2(igroup_2(t)),ampterm_2(t)
 enddo
 return
92 write(ulog,*)'READ_FCS: rank=',rank,' reached end of file after t=',t
!----------------------------------------

elseif ( rank .eq. 3) then

 res =  igroup_1(nterms(1)) + igroup_2(nterms(2))
 do i=1,nterms(rank)
       read(iunit,*,err=93)t,igroup_3(t), &
& iatomterm_3(1,t),ixyzterm_3(1,t),iatomterm_3(2,t),ixyzterm_3(2,t), &
& iatomterm_3(3,t),ixyzterm_3(3,t),  &
& fcs_3(igroup_3(t)),ampterm_3(t)
 enddo
 return
93 write(ulog,*)'READ_FCS: rank=',rank,' reached end of file after t=',t
!----------------------------------------

elseif ( rank .eq. 4) then

 res =  igroup_1(nterms(1)) + igroup_2(nterms(2)) + igroup_3(nterms(3))
 do i=1,nterms(rank)
       read(iunit,*,err=94)t,igroup_4(t), &
& iatomterm_4(1,t),ixyzterm_4(1,t),iatomterm_4(2,t),ixyzterm_4(2,t), &
& iatomterm_4(3,t),ixyzterm_4(3,t),iatomterm_4(4,t),ixyzterm_4(4,t), &
& fcs_4(igroup_4(t)),ampterm_4(t)
 enddo
 return
94 write(ulog,*)'READ_FCS: rank=',rank,' reached end of file after t=',t
else

 write(ulog,*)' READ_FCS: rank must be from 1 to 4, not ',rank

endif

write(ulog,*)'READ_FCS: error!! should not have gotten here!'

stop

 end subroutine read_fcs_2
!
!-------------------------BELOW IS FOR SELF-ENERGY PART---------------------
 SUBROUTINE mySelfEnergy(ksub,tempk)
 !!legacy code for self energy, currently not used
        IMPLICIT NONE
        INTEGER, INTENT(in) :: ksub(2) ! for pool size
        INTEGER ::nkbs,lambda,nqs,monitor,i,j,om_idx
        INTEGER :: check_k1,check_l1
        REAL(8), INTENT(in) :: tempk
        REAL(8), ALLOCATABLE :: omega(:) ! temp

        INTEGER :: om_number,step
        REAL(8) :: max_om,delta_om
        COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE :: new_nself, new_uself

        !------------ prepare input params-------------
        CALL read_fc23
        CALL updateFCS2 !maybe no use

        const33 = 1d30*(hbar/uma)**1.5d0/(sqrt(ee*1d20/uma))**1.5d0*ee/h_plank/100/c_light
        ndyn = d*atom_number
        nkbs = SIZE(kp_bs,dim=2)


        ! for test purpose
!        check_k1 = 22; check_l1 = 1
!        CALL function_self_36(check_k1,kp_bs(:,check_k1),check_l1,&
!        &eivals(check_l1,check_k1),tempk,nself(check_k1,check_l1), uself(check_k1,check_l1))
!        WRITE(*,*) "Check value nself:", nself(check_k1,check_l1)
!        WRITE(*,*) "Check value uself:", uself(check_k1,check_l1)
!STOP
        ! output for check purpose
        monitor = 421
        OPEN(monitor, FILE='self.dat', STATUS='unknown')

        ! compatible with old way, turn off for the new way
!        ALLOCATE(omega(ndyn*nkbs))
!        IF(ALLOCATED(nself)) DEALLOCATE(nself)
!        IF(ALLOCATED(uself)) DEALLOCATE(uself)
!        ALLOCATE(nself(ksub(2),ndyn),uself(ksub(2),ndyn))

        ! compatible with new way, turn off for the old way
        om_number = 100
        max_om = MAXVAL(MAXVAL(eivals,DIM=2))
        delta_om = max_om/om_number

        ALLOCATE(new_nself(ksub(2),ndyn,om_number+1))
        ALLOCATE(new_uself(ksub(2),ndyn,om_number+1))

        IF(ALLOCATED(nself)) DEALLOCATE(nself)
        IF(ALLOCATED(uself)) DEALLOCATE(uself)
        ALLOCATE(nself(ksub(2),om_number),uself(ksub(2),om_number))

        IF(ALLOCATED(omega)) DEALLOCATE(omega)
        ALLOCATE(omega(om_number))

        ! q, lambda loop for every omega

        DO nqs = ksub(1), ksub(2)
            ! old way: only calculate nself & uself on the kbps path of eivals
!        DO lambda = 1, ndyn
!            IF (ALL(kp_bs(:,nqs).eq.kp_gamma)) CYCLE ! skip if gamma point
!            om_idx = lambda + ndyn*(nqs-1)
!            omega(om_idx) = eivals(lambda, nqs) !om1
!            CALL function_self_36(nqs,kp_bs(:,nqs),lambda,omega(om_idx),tempk,nself(nqs, lambda), uself(nqs, lambda))
!            WRITE(monitor, 7) 'nqs=',nqs, 'lambda=',lambda
!            WRITE(monitor, *) 'nself=',nself(nqs,lambda)
!            WRITE(monitor, *) 'uself=',uself(nqs,lambda)
!            WRITE(monitor,*) '----------------------------------'
!        END DO !lambda loop
            ! new way: calculate nself & uself on selected omega values, another nested loop
        DO step=1, om_number
        DO lambda=1,ndyn
            omega(step) = delta_om*step !can't start from 0, because of den
            CALL function_self_36(nqs,kp_bs(:,nqs),lambda,omega(step),tempk,new_nself(nqs,lambda,step),new_uself(nqs,lambda,step))
            !sum up the lambda?
            nself(nqs,step) = nself(nqs,step) + new_nself(nqs,lambda,step)
            uself(nqs,step) = uself(nqs,step) + new_uself(nqs,lambda,step)
        END DO !lambda loop
            WRITE(monitor, 8) 'nqs=',nqs, 'omega^2=',omega(step)
            WRITE(monitor, *) 'nself=',nself(nqs,step)
            WRITE(monitor, *) 'uself=',uself(nqs,step)

        END DO !om loop

        END DO !q loop

        ! compatible output with new way
        OPEN(443,FILE='nself_contour.dat',STATUS='unknown')
        OPEN(445,FILE='uself_contour.dat',STATUS='unknown')

        DO nqs=ksub(1),ksub(2)
        DO step=1,om_number+1
            WRITE(443,*) nqs, delta_om*step,nself(nqs,step)
            WRITE(445,*) nqs, delta_om*step,uself(nqs,step)
        END DO !step
        END DO !nqs
        CLOSE(443)
        CLOSE(445)

        DEALLOCATE(omega)
        DEALLOCATE(new_nself,new_uself)
        CLOSE(monitor)
7 FORMAT(2(a,I5,2x))
8 FORMAT(1(A,I5,2X),1(A,G16.7,2X))
    END SUBROUTINE mySelfEnergy
!------------------------------------------------------------------------------------------------------------------------------
subroutine function_self_36(nqs,q,la,omega,temp,nself,uself)
!!legacy code for self energy, currently not used
!!notes: I modified this a bit
!!this version uses the 5-component array v33_squared
!! frequencies, temperatures are in cm^-1 (same unit as wmax)
!! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
!! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
!! it assumes the input momentum q is on the already-generated kibz mesh and
!! the second argument of V33(q,k2,k3) k2 covers the FBZ and similarly k3=-q-k2

 use kpoints
 use om_dos
 use phi3
 use eigen
 use geometry
!use params
!use constants
!use lattice
!use io2
 implicit none
 integer, intent(in) :: la,nqs
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer :: i,j,k,nk2,i2,j2,t,i0,j0,k0,i3,al,be,ga
 integer :: ta1, ta2, ta3
 integer iq,jq,kq,inside,ik,jk,kk,nq,jk1,jk2,jk3,n2,nk3,l2,l3
! integer, save :: cntr
 integer :: check_file
 real(8) om1,om2,om3,eta,tk,v32,k2(3),k3(3),v33sqr,epsl,qdr
 real(8) mi,mj,mk,rx2(3),rx3(3),den,curr_small
 real(8), allocatable :: kpc3(:,:),eigenval3(:,:),junk(:,:,:)
 TYPE(vector) :: rr(d)
 complex(8) omz,xx,self,eiqr,ev1,ev2,ev3
 complex(8),allocatable :: eigenvec3(:,:,:)
 complex(8) ocfunc_36
 integer grid

curr_small = 100d0
!*******

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm

! input q must be a IBZ vector, nq is its eivalibz index
! call get_k_info_cent(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)

 write(ulog,5)'SELF_35: nqs, omega,etaz=',nqs,omega,etaz
 v32 = v3_threshold*v3_threshold
 nself=0 ; uself=0;
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 !cntr = 0

! for every input q, calculate om3(l3,k3) where k3= -q-k2
! it's ON-THE-FLY so this part is inside and flushed everytime the input q changes

    ALLOCATE(dkc(nkc))
    dkc=0
    ALLOCATE(kpc3(3,nkc),eigenval3(ndyn,nkc),eigenvec3(ndyn,ndyn,nkc))
    DO i = 1, nkc
        kpc3(:,i) = -q(:) - kpc(:,i)
    END DO
    !ALLOCATE(junk(d,ndyn,nkc))
    CALL get_frequencies(nkc,kpc3,dkc,ndyn,eigenval3,ndyn,eigenvec3,veloc)
    DEALLOCATE(dkc)

check_file = 435
OPEN(check_file,FILE='k1_22.dat',STATUS='unknown')

 k2loop: do n2=1,nkc

            k2=kpc(:,n2)

            ! check if k2 or k3 is gamma or its equivalents?
            !call get_k_info(k2-shft,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
            if(n2.eq.1) cycle k2loop

            k3=-q-k2

            if (ALL(k3.eq.kp_gamma)) cycle k2loop
            if (sqrt(sum(k3*k3)).lt.1e-3) cycle k2loop! in case of rounding error
            ! in case -q-k2 brings k3 to the edge which is actually gamma point
            ! 8 possible cases
            grid = 1
            CALL check_g(k3,rr1,grid)
            CALL check_g(k3,rr2,grid)
            CALL check_g(k3,rr3,grid)

            if (grid.eq.1) cycle k2loop   ! k3 is a G-vector in which case V33=0

!notes: ndyn = eigen_number = d*atom_number is calculated outside
        do l2=1,ndyn

        do l3=1,ndyn

    !-----------------calculation part for the v35 substitute-------------
            xx = cmplx(0d0,0d0)
tloop:  do t=1,nterms(3)

           i3 = iatomterm_3(1,t)
           i0 = iatomcell0(i3)
           if (i0.ne.i3) cycle tloop

           al = ixyzterm_3(1,t)
           be = ixyzterm_3(2,t)
           ga = ixyzterm_3(3,t)
           mi = atom0(i0)%mass
           j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
           k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
           ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
           rx2(:) = atompos(:,j) - atompos(:,j0)
           rx3(:) = atompos(:,k) - atompos(:,k0)
    !  be careful: it has to be the translations R not atompos!

           eiqr = exp( ci* ((k2 .dot. rx2) + (k3 .dot. rx3)) )

            ev1 = eivecs(ta1,la,nqs)!calculated outside
            ev2 = eigenvec(ta2,l2,n2)!calculated outside
            ev3 = eigenvec3(ta3,l3,n2)!calculated inside

           xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
           & ev1 * ev2 * ev3
    !&      eivcibz(ta1,l1,n1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
    ! make sure for kp at the FBZ boundary it gives (-1) or 1


    ! n2 is already in the IBZ, n3 is not but for now we are using the full mesh
    ! and nibz=nkc;  therefore n1 is also in the full zone and is OK, otherwise
    ! for n1 and n3 it is necessary to rotate the eigenvector from the IBZ to its
    ! true position in the recirpocal space. (the op_kmatrix will do the job)
    !call naiveCheck(n1,nk2,nk3, l1,l2,l3,ta1,ta2,ta3,t, eiqr,xx, counter)
        enddo tloop
            om1 = omz !eivals(la,nqs) !calculated outside
            om2 = eigenval(l2,n2) !calculated outside
            om3 = eigenval3(l3,n2) !calculated inside

        den = sqrt(om1*om2*om3)/cnst/sqrt(cnst)

IF(den.le.curr_small) curr_small = den
!WRITE(check_file,*) '==========================================='
!WRITE(check_file,*) 'n2=',n2,'l2=',l2,'l3=',l3
!WRITE(check_file,*) 'om1=',om1
!WRITE(check_file,*) 'om2=',om2
!WRITE(check_file,*) 'om3=',om3
!WRITE(check_file,*) 'den=',den
!WRITE(check_file,*)

        if (den.ne.0) then
           xx = xx / den/2d0/sqrt(2d0) * const33

!WRITE(check_file,*) 'xx=',xx

        else   ! usually acoustic eigenvalues at k=0 are not strictly zero, but a small number/original comment
            ! but that small number is 1E-18 which is on the denominator and explodes
           stop
           ! instead of stop, just skip over that den will make it too small

        endif


        v33sqr = xx*conjg(xx) !v3sq(nqs,n2,la,l2,l3)


        if(v33sqr.lt.v32) cycle ! drop the ones that are below the v3_threshold


        self = ocfunc_36(temp,omz,om2,om3)!need a modified ocfunc(temp,omz,l2,n2,l3,nk3)
!WRITE(check_file,*) 'v33sqr=',v33sqr
!WRITE(check_file,*) 'self=',self
!WRITE(check_file,*)
!WRITE(check_file,*) 'This is added into the sum: v33sqr*self=', v33sqr*self
        inside = 1
    !---------------- figure out the inside param-----------------------------
        epsl=5d-10
        !WRITE(*,*)'Check the rr value', rr%component !should they be 0?
        !WRITE(*,*)'Check the q value', q
        rr = (/rr1,rr2,rr3/)
        DO i = 1, d
            qdr = (k3.dot.rr(i)) + epsl
            !PRINT *, 'qdr=',qdr
            ! include from -0.00005 to 0.99995 included (include 0.0 and exclude 1.0)
            if (qdr.lt.0 .or. qdr.ge.1) THEN
                inside = 0
                exit
            end if
        END DO
    !------------------------------------------------------------------------

       if(inside.eq.1) then                ! normal process
            nself = nself + v33sqr * self
                 !WRITE(*,*) 'QUICK CHECK: nself=',nself
       else                                ! umklapp process
            uself = uself + v33sqr * self
       endif
    enddo ! l3
    enddo ! l2
 enddo k2loop
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la, kpc(:,nq),omega,self, nself,uself
!
WRITE(check_file,*)
WRITE(check_file,*) 'Minimum den=', curr_small
WRITE(check_file,*) 'const33=',const33
CLOSE(check_file)

3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i4,9(2x,3(1x,f7.3)))

 end subroutine function_self_36
!--------------------------------------------------------------------------------------------------
    SUBROUTINE check_g(kpoint,rr,grid)
    !! tells if k3 is on a g-vector, grid will be 1
        IMPLICIT NONE
        REAL(8) :: g(3), temp
        REAL(8),INTENT(in) :: kpoint(3)
        TYPE(vector),INTENT(in) :: rr
        INTEGER, INTENT(inout) :: grid
!        LOGICAL, INTENT(out) :: grid

        REAL(8) aux,qdr,epsl
        integer n,i,insid
        temp = kpoint.dot.rr
        aux = temp-FLOOR(temp)

        IF(aux.gt.1E-4) THEN
            grid = 0*grid
        ELSE
            grid = grid
        END IF

    END SUBROUTINE check_g
!--------------------------------------------------------------------------------------------------
    SUBROUTINE mySelfEnergy_re(ksub, tempk,filename)
        !! v3sq(:,:,:,:,:) flattened into v33sq(:) and nq1(:),nq2(:),nq3(:)
    !! only works on kpc(:,:) mesh
        IMPLICIT NONE
        INTEGER, INTENT(in) :: ksub(2) ! for pool size
        INTEGER ::lambda,q,monitor,unit_number5,i,j
        REAL(8), INTENT(in) ::  tempk
        REAL(8), ALLOCATABLE :: omega(:) ! temp
        CHARACTER(30),INTENT(in) :: filename
        CHARACTER(30) :: file1, file2
        !------------ prepare input params-------------
        const33 = 1d30*(hbar/uma)**1.5d0/(sqrt(ee*1d20/uma))**1.5d0*ee/h_plank/100/c_light

        ndyn = d*atom_number
        !allocate eivalibz(ndyn,nibz),eivecibz(ndyn,ndyn,nibz),velocibz(3,ndyn,nibz) and also dkc(nibz)
        ALLOCATE(omega(ndyn))
        CALL allocate_eig(ndyn,nibz,nkc)

        ALLOCATE(dkc(nkc))
        dkc=0
        CALL get_frequencies(nkc,kpc,dkc,ndyn,eigenval,ndyn,eigenvec,veloc)
        DEALLOCATE(dkc)  !,eivalibz,eivecibz)

        ! irreducible Brillouin zone part is not changed
        ALLOCATE(dkc(nibz))
        dkc=0
        CALL get_frequencies(nibz,kibz,dkc,ndyn,eivalibz,ndyn,eivecibz,velocibz)
        DEALLOCATE(dkc)

        !-----------------monitor eigenvalues-------------
        unit_number5 = 438
        file1 = 'eigenvalue_' //TRIM(filename)
        OPEN(unit_number5, file=TRIM(file1),status='unknown')
        WRITE(unit_number5,*) "check const33=", const33
        DO i=1, SIZE(eigenval, dim=1)
            DO j=1, SIZE(eigenval, dim=2)
                WRITE(unit_number5,*) "nb=", i, "nk=",j, "eigenval=",eigenval(i,j)
            END DO
        END DO
        CLOSE(unit_number5)

        CALL read_fc23
        nv3split=nkc*(ksub(2)-ksub(1)+1)*ndyn**3
        CALL allocate_v33(nv3split)
        CALL calculate_w3_ibz_split_sq(ksub,nv3split,ndyn,nkc,kpc)

        ! for test purpose
!        CALL function_self_new_sq(kibz(:,2),2,eivalibz(2,2),tempk,nself(2,2), uself(2,2))
!        WRITE(*,*) "Check value nself:", nself(2,2)
!        WRITE(*,*) "Check value uself:", uself(2,2)

        ! output for check purpose
        monitor = 421
        file2 = 'self_'//TRIM(filename)
        OPEN(monitor, FILE=TRIM(file2),STATUS='unknown')

        ! q, lambda loop for every omega
        DO lambda = 1, ndyn
            DO q = ksub(1), ksub(2)
                omega(lambda) = eivalibz(lambda, q)
                CALL function_self_new_sq(kibz(:,q),lambda,omega(lambda),tempk,nself(q,lambda),uself(q,lambda))
                WRITE(monitor, 7) 'q=',q, 'lambda=',lambda
                WRITE(monitor, *) 'nself=',nself(q,lambda)
                WRITE(monitor, *) 'uself=',uself(q,lambda)
                WRITE(monitor,*) '----------------------------------'
            END DO !q loop
        END DO !lambda loop

        CLOSE(unit_number5)
        CLOSE(monitor)
7 FORMAT(2(a,I5,2x))
    END SUBROUTINE mySelfEnergy_re
!------------------------------------------------------------------
 subroutine read_ksubset (ksubset,nib)
 !!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!! this subroutine reads ksubset.inp
 use io2
 implicit none
 integer ksubset_f,tmp,nib
 integer ksubset(2)

 ksubset_f=9998
 ksubset=0
 open(ksubset_f,file='ksubset.inp', status='old')
 read(ksubset_f,*) ksubset
 close(ksubset_f)
 if (ksubset(2).lt.ksubset(1)) then
     write(ulog,*) '!!!!!!!!!!!WARNING!!!!!!!!KSUBSET=', ksubset(1), ksubset(2)
     write(ulog,*) 'ckeck your ksubset input'
     tmp=ksubset(1)
     ksubset(1)=ksubset(2)
     ksubset(2)=tmp
 endif
 if (ksubset(2) .gt. nib) then
     write(ulog,*) '!!!!!!!!!!!WARNING!!!!!!!!KSUBSET=', ksubset(1), ksubset(2)
     write(ulog,*) 'ckeck your ksubset input; k(2)>nibz; setting k(2)=nibz'
     ksubset(2)=nib
 endif
 end subroutine read_ksubset
!=========================================================================================
 subroutine matrix_kpt_sy
!!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!! print out k2 points, omega, sum(1,3) V33(k1,la1,k2,la2,k3,la3)^2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer j,nq,jq,la,jk1,jk2,jk3,nqq2,nk1,inside,q2,l2
 real(8) omq2(ndyn),sumz(ndyn)

 open(udos,file='v33square_k.dat', status='unknown')
 write(udos,*) 'i, q, omq, v33square'

 do q2=1,nibz
 sumz=0
 do l2=1,ndyn
    call get_k_info(kibz(:,q2),nc,nqq2,jk1,jk2,jk3,g1,g2,g3,inside)
    omq2(l2) = eigenval(l2,nqq2)
    !write(debug,*) 'start v3loop...q2,l2,nv3', q2, l2, nv3
    v3loop: do j=1,nv3
      !write(debug,*) 'j, la2(j), l2, nq2(j), nqq2', j, la2(j), l2, nq2(j), nqq2
      if (la2(j) .ne. l2) cycle v3loop
      if (nq2(j) .ne. nqq2) cycle v3loop
      sumz(l2) = sumz(l2) + (v33(j)*conjg(v33(j)))
      !write(debug,*) 'sumz, v33(j)', sumz, v33(j)
    enddo v3loop

    !write(debug,*) 'mapinv(q2),kibz(:,q2), wibz(q2)', mapinv(q2), kibz(:,q2), wibz(q2)
    sumz(l2)=sumz(l2)*wibz(q2)

    !write(ulog,*) 'matrix_kpt_sy, q2, nibz', q2, nibz
    !write(debug,*) 'kibz(:,q2), wibz(q2)', kibz(:,q2), wibz(q2)

  enddo
  write(udos,8) mapinv(q2), kibz(:,q2), (omq2(l2),sumz(l2),l2=1,ndyn)
  enddo

 close(udos)

8 format(1x,9(1x,g11.4))

 end subroutine matrix_kpt_sy
!==================================================================================================
 subroutine matrix_elt(q1,q2,l1,l2,l3,w33,inside)
 !!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!!for input modes (qi,li) i=1,2,3 calculates the 3-phonon matrix element w33(1,2,3) on the fly
 use kpoints
 use phi3
 use svd_stuff
 use eigen
 use atoms_force_constants
 use om_dos

 implicit none
 integer, intent(in) :: l1,l2,l3
 real(8), intent(in) :: q1(3),q2(3)
 integer, intent(out) :: inside
 complex(8), intent (out) :: w33
 integer i3,nk1,nk2,nk3,i0,al,be,ga,j,j0,k,k0,t,ta1,ta2,ta3,i1,j1,k1
 real(8) q3(3),mi,mj,mk,rx2(3),rx3(3),den
 complex(8) xx,eiqr

! write(*,*)'matrix_elt',q1,q2,l1,l2,l3
! write(*,*)'matrix_elt, const33=',const33

!***************************
! here we assume q1 is IBZ (non-shifted) but q2 is FBZ and shifted
!***************************???

! call get_k_info_cent(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
! call get_k_info_cent(q2-shft,NC,nk2,i1,j1,k1,g1,g2,g3,inside)
 call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
 call get_k_info(q2-shft,NC,nk2,i1,j1,k1,g1,g2,g3,inside)
 w33 = cmplx(0d0,0d0)
 q3 = -q2-q1
! call get_k_info_cent(q3-shft,NC,nk3,i1,j1,k1,g1,g2,g3,inside
 call get_k_info(q3-shft,NC,nk3,i1,j1,k1,g1,g2,g3,inside)
! nk1, nk2 are checked outside but nk3 could still be 1
IF(nk1.eq.1 .OR. nk2.eq.1 .OR. nk3.eq.1) THEN
    w33 = 0d0
    STOP
END IF
! write(*,*)'matrix_elt: q3=',q3,nk3,inside
 xx = cmplx(0d0,0d0)

 tloop: do t=1,nterms(3)

! write(*,*)'matrix_elt: tloop=',t,nterms(3)
       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop ! first index must be in primitive cell
       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j) - atompos(:,j0)
       rx3(:) = atompos(:,k) - atompos(:,k0)
       eiqr = exp( ci* ((q2 .dot. rx2) + (q3 .dot. rx3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
!&      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)
&      eivecibz(ta1,l1,mapibz(nk1))*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)

 enddo tloop
! den = sqrt(8*eigenval(l1,nk1)*eigenval(l2,nk2)*eigenval(l3,nk3))
 den = sqrt(8*eivalibz(l1,mapibz(nk1))*eigenval(l2,nk2)*eigenval(l3,nk3))
 if (den.ne.0) then
! the last term cnst^1.5 is new and corrects for previous cnst*sqrt(eiv)
    w33 = xx / den * const33 *sqrt(cnst*cnst*cnst)
 else
    write(*,*)'MATRIX_ELT: den=0 ',den
    stop
 endif

 end subroutine matrix_elt
!==================================================================================================
subroutine matrix_elt_full(q,q2,q3,omq,om2,om3,evq,ev2,ev3,w33)
!!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!!for input modes (qi,li) i=1,2,3 calculates the 3-phonon matrix element w33(1,2,3) on the fly
 use kpoints
 use phi3
 use svd_stuff
 use eigen
 use atoms_force_constants
 use om_dos

!use lattice
!use io2  **included in MODULE kpoints
!use params ** included in MODULE lattice
!use constants ** included in MODULE lattice
!use force_constants_module ** included in MODULE atoms_force_constants
 implicit none
! integer, intent(in) :: la
 real(8), intent(in) :: q(3),q2(3),q3(3),omq,om2,om3
 complex(8), intent(in) :: evq(ndyn),ev2(ndyn),ev3(ndyn)
! integer, intent(out) :: inside
 complex(8), intent (out) :: w33
 integer i3,nk1,nk2,nk3,al,be,ga,j1,k1,i0,j0,k0,t,ta1,ta2,ta3
 real(8) mi,mj,mk,rx2(3),rx3(3),den
 complex(8) xx,eiqr

! write(*,*)'matrix_elt',q1,q2,l1,l2,l3
! write(*,*)'matrix_elt, const33=',const33
 w33 = cmplx(0d0,0d0)
 xx = cmplx(0d0,0d0)

 tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop ! first index must be in primitive cell
       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j1 = iatomterm_3(2,t)  ;  j0 = iatomcell0(j1) ;  mj = atom0(j0)%mass
       k1 = iatomterm_3(3,t)  ;  k0 = iatomcell0(k1) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j1) - atompos(:,j0)
       rx3(:) = atompos(:,k1) - atompos(:,k0)
       eiqr = exp( ci* ((q2.dot. rx2) + (q3 .dot. rx3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
&      evq(ta1)*ev2(ta2)*ev3(ta3)
! &      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)

 enddo tloop
 den = sqrt(8*omq*om2*om3/(cnst*cnst*cnst))
 if (den.ne.0) then
    w33 = xx / den * const33
 else
    write(*,*)'MATRIX_ELT: den=0 ',den
    stop
 endif

 end subroutine matrix_elt_full
!==================================================================================================
!subroutine get_freq is in extratools.f90
!==================================================================================================
subroutine function_self_tetra(nqs,q,la,omega,temp,nself,uself)
!!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!! the tetrahedron method requires direct access to v3 for given k,la indices.
!! we therefore assume v3sq is calculated and used
!! frequencies, temperatures are in cm^-1 (same unit as wmax)
!! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
!! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
!! it assumes the 2nd momentum k2 is on the already-generated kmesh and
!! is the second argument of V33(q,k2,k3) where k2 covers the FBZ and k3=-q-k2
 use phi3
 use eigen
 use om_dos
 use tetrahedron
!use lattice
!use kpoints
!use io2
!use params
!use constants
 implicit none
 integer, intent(in) :: la,nqs
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,nk1,nk2,nk3,l1,l2,l3
 integer, save :: cntr
 real(8) om2,nb2,om3,nb3,nbe,resm,resp,eta,tk,k1(3),k3(3),q1(3),q3(3)
 complex(8) omz,self,ocfunc
 real(8), allocatable :: eivp(:),eivm(:),funcp(:),funcm(:),junkp(:),junkm(:)
 integer, allocatable:: ins(:)

 allocate(eivp(nkc),funcp(nkc),eivm(nkc),funcm(nkc),junkp(nkc),junkm(nkc),ins(nkc))

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
! call get_k_info_cent(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! jq = (nq-1)*ndyn+la  !nq+(la-1)*nkc

! omq = sqrt(abs(eigenval(la,nq))) * cnst
! write(ulog,5)'SELF_NEW: last cntr, omega,etaz=',cntr,omega,etaz
 nself=0 ; uself=0;
!cntr = 0

 do l2=1,ndyn
 do l3=1,ndyn

    do nk2=1,nkc
       k3=-q-kpc(:,nk2)
!      call get_k_info_cent(k3-shft,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
       call get_k_info(k3-shft,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
       ins(nk2)=inside
       om3 = eigenval(l3,nk3) ; nb3=nbe(om3,temp,classical)
       om2 = eigenval(l2,nk2) ; nb2=nbe(om2,temp,classical)

! define the arguments of delta(omega-om1-om3)
       eivp(nk2)=om2+om3
       eivm(nk2)=om2-om3
! define the weighting function
       funcp(nk2)=v3sq(nqs,nk2,la,l2,l3)*(1+nb2+nb3)
       funcm(nk2)=v3sq(nqs,nk2,la,l2,l3)*(nb3-nb2)
    enddo

    call tet_sum(omega,nkc,eivp,funcp,resp,junkp)  ! res= sum_k func(k)*delta(omega-eiv(k))
    call tet_sum(omega,nkc,eivm,funcm,resm,junkm)  ! res= sum_k func(k)*delta(omega-eiv(k))

    do nk2=1,nkc
!   cntr = cntr+1
       if(ins(nk2).eq.1) then                ! normal process
!        nself = nself + cmplx(0,resp)/2
         nself = nself + cmplx(0,junkp(nk2))/2
       else                                ! umklapp process
!        uself = uself + cmplx(0,res)/2
         uself = uself + cmplx(0,junkp(nk2))/2
       endif

!   cntr = cntr+1
       if(ins(nk2).eq.1) then                ! normal process
!        nself = nself + cmplx(0,resm)
         nself = nself + cmplx(0,junkm(nk2))
       else                                ! umklapp process
!        uself = uself + cmplx(0,res)
         nself = nself + cmplx(0,junkm(nk2))
       endif
    enddo

 enddo
 enddo

 nself = nself  *2*pi*pi
 uself = uself  *2*pi*pi
! self  = nself + uself

 deallocate(eivp,eivm,funcp,funcm,junkp,junkm,ins)

! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_tetra
!==================================================================================================

!--------------------------------------------------------------------------------------------------
subroutine function_self_new_sq(q,la,omega,temp,nself,uself)
!!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!! this version uses the general index in v33sq
!! frequencies, temperatures are in cm^-1 (same unit as wmax)
!! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
!! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
!! it assumes the input momentum q is on the already-generated kmesh and
!! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1,nk3,l1,l2,l3
 integer, save :: cntr
 real(8) om1,eta,tk,v32,k2(3),k3(3),q1(3) !,sh(3)
 complex(8) omz,self,ocfunc

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
! sh=(-0.5d0)*(g1+g2+g3)

! input q must be a IBZ vector
 !call get_k_info_cent(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
 jq = (nq-1)*ndyn+la  !nq+(la-1)*nkc

! om1 = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 write(ulog,5)'SELF_NEW: last cntr, omega,etaz=',cntr,omega,etaz
! omz = cmplx(omega,-etaz)   ! this is in cm^-1
 v32 = v3_threshold*v3_threshold
 nself=0 ; uself=0;
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 cntr = 0
 nv3loop: do j=1,nv3split
!  q must be the 2nd argument in v33
    if ( jq .ne. nq1(j) ) cycle nv3loop
    if ( v33sq(j) < v32) cycle nv3loop

    l3=1+mod(nq3(j)-1,ndyn) ; jk3=1+(nq3(j)-l3)/ndyn
    l2=1+mod(nq2(j)-1,ndyn) ; jk2=1+(nq2(j)-l2)/ndyn
    l1=1+mod(nq1(j)-1,ndyn) ; jk1=1+(nq1(j)-l1)/ndyn  ! this is (q,la):
    k2=kpc(:,jk2)
    k3 = -q-k2
    !call get_k_info_cent(k3-shft,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
    call get_k_info(k3-shft,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
    if(jk1.ne.nq) then
      write(ulog,3)'SELF_NEW: nq.ne.jk1 ',nq,jk1
      stop
    endif
    if(nk3.ne.jk3) then
      write(ulog,3)'SELF_NEW: nk3.ne.jk3 ',nk3,jk3
      write(ulog,3)'Check consistency between nki and those of v33.dat'
      write(ulog,3)'SELF_NEW:j,l,q,om=',j,la,q,omega
      write(ulog,4)'SELF_NEW: q,k2,k3=',q,k2,k3
      write(ulog,*)'SELF_NEW: nqi    =',nq1(j),nq2(j),nq3(j)
      write(ulog,*)'SELF_NEW:l1,l2,l3=',l1,l2,l3
      write(ulog,4)'SELF_NEW:kpc(jki)=',kpc(:,jk1),kpc(:,jk2),kpc(:,jk3)
!      write(ulog,4)'SELF_NEW:la,la2(j)=',la,la2(j)
      write(*,3)'SELF_NEW: nk3.ne.jk3 ',nk3,jk3
      write(*,3)'SELF_NEW: j,la,q,om =',j,la,q,omega
      stop
    endif

    cntr = cntr+1
    self = ocfunc(temp,omz,l2,jk2,l3,jk3)
    if(inside.eq.1) then                ! normal process
        nself = nself + v33sq(j) * self
    else                                ! umklapp process
        uself = uself + v33sq(j) * self
    endif
 enddo nv3loop
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself


3 format(a,i6,i6,3x,3(1x,f8.3),1x,g12.5,1x,9(2x,2(1x,g11.4)))
4 format(a,9(1x,g11.4))
5 format(a,i8,99(1x,g11.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_new_sq
!--------------------------------------------------------------------------------------------------
 subroutine function_self_w2_sy(q,la,omega,temp,nself,uself)
!!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!! calculates the self-energy on the fly for q-point in the generated kmesh and Lorentzian delta
!! frequencies, temperatures are in cm^-1 (same unit as wmax)
!! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
!! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
!! it assumes the input momentum q is on the already-generated kmesh and
!! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use phi3
 use eigen
 use om_dos
 !use params
 !use constants
 !use io2
 !use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1,nk3,ins2,ifbz,l1,l3
 integer, save :: cntr
 real(8) omq,om3,omk,eta,tk,term,k1(3),k2(3),q1(3),q3(3)
 complex(8) omz,self,ocfunc,oc2,s2,xx


 integer k_FBZ,q_ibz
 real(8) v33_square_sum,tempk
 complex(8) delta_tot_sum, delta1_sum, delta2_sum, delta3_sum, delta4_sum
 complex(8) delta_tot, delta1, delta2, delta3, delta4
 v33_square_sum=0; delta_tot_sum=0; delta1_sum=0; delta2_sum=0; delta3_sum=0; delta4_sum=0

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
! write(*,*) 'self_w: tempk=',tk
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! write(*,*) 'self_w: k_info: nq=',nq
 omq = eigenval(la,nq)
 eta= etaz
! write(ulog,5)'SELF_W: last cntr, omega,etaz=',cntr,omega,etaz
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0

 write(ulog,*) 'entered self energy calculation... at q and la', nq, la

 cntr = 0
 do ifbz=1,nkc
    q3=kpc(:,ifbz)
    call get_k_info(q3,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
    if (nk3.ne.ifbz) then
       write(ulog,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       write(*,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       stop
    endif
    k1 = -q-q3
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,ins2)

    do l3=1,ndyn
       om3 = eigenval(l3,nk3)
    do l1=1,ndyn
       omk = eigenval(l1,nk1)
!      call matrix_elt     (q,q3,la,l3,l1,xx,inside)
       call matrix_elt_full(q,q3,k1,omq,om3,omk,eigenvec(:,la,nq),eigenvec(:,l3,nk3),eigenvec(:,l1,nk1),xx)
       call check_inside_bz(k1,g1,g2,g3,inside)
       if(inside.ne.ins2) then
          write(*,*)'SELF_W:ERROR:inside .ne. ins2 ',inside,ins2,k1
          stop
       endif

       ! for raw data output
       call ocfunc_sy(temp,omz,l1,nk1,l3,nk3,delta1,delta2,delta3,delta4,delta_tot)

       v33_square_sum=v33_square_sum+xx*conjg(xx)
       delta1_sum=delta1_sum+delta1
       delta2_sum=delta2_sum+delta2
       delta3_sum=delta3_sum+delta3
       delta4_sum=delta4_sum+delta4
       delta_tot_sum=delta_tot_sum+delta_tot
       ! end for raw data output


       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)

       jk1 = nk1+(l1-1)*nkc
       jk3 = ifbz+(l3-1)*nkc

       self = ocfunc(temp,omz,l1,nk1,l3,nk3)
       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif
    enddo
    enddo

 enddo
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself


! Write Raw data.....
! this is to make sure k and all its stars have the same self energy
 tempk=temp*(100*h_plank*c_light)/k_b
! do k_FBZ=1, nkc
!   write(*,*) 'k_FBZ, mapibz(k_FBZ), nq', k_FBZ, mapibz(k_FBZ), nq
!   if(mapibz(k_FBZ) .eq. q_ibz) then   ! IBZ -> FBZ
!     write(*,*) 'matched--------------------------------------------------------'
!     write(self_detail,8) tempk, k_FBZ, kpc(:,k_FBZ), la, omega, &  ! ki, kpoint, lambda, omega
!& v33_square_sum, real(delta1_sum), aimag(delta1_sum), real(delta2_sum), aimag(delta2_sum) &
!& , real(delta3_sum), aimag(delta3_sum), real(delta4_sum), aimag(delta4_sum), real(delta_tot_sum), aimag(delta_tot_sum)
!   endif
! enddo

   write(self_detail,8) tempk, q, la, omega, &
& v33_square_sum, real(delta1_sum), aimag(delta1_sum), real(delta2_sum), aimag(delta2_sum) &
& , real(delta3_sum), aimag(delta3_sum), real(delta4_sum), aimag(delta4_sum), real(delta_tot_sum), aimag(delta_tot_sum)


! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
6 format(i6,6(1x,f9.4),2x,3(1x,i2))
7 format(a,4i3,99(1x,f7.3))
8 format(99(5x,g16.8))

 end subroutine function_self_w2_sy
!--------------------------------------------------------------------------------------------------
subroutine function_self_w3(q,la,omega,temp,nself,uself)
!!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!! this one is for arbitrary q
!! uses lorentzian delta with 4 terms, and calculates both real and imaginary parts
!! In this subroutine, which calculates like self_w on the fly, the second momentum
!! is restricted to be on the kpc mesh, but the first and third momenta can be anywhere
!! frequencies, temperatures are in cm^-1 (same unit as wmax)
!! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
!! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
!! the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use om_dos
 use phi3
 use eigen

 !use io2
 !use params
 !use constants
 !use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer inside,ik,l3,l2
! integer, save :: cntr
 real(8) omq,eta,term,k3(3),k2(3),eivq(ndyn),eiv2(ndyn),eiv3(ndyn),vg(3,ndyn)
 real(8) etacut,arg,delta_l,om3,om2,nb3,nb2,nbe,ires,rres
 complex(8) omz,self,xx,evq(ndyn,ndyn),ev2(ndyn,ndyn),ev3(ndyn,ndyn)

! be careful band sorting should be avoided in this case; otherwise need to
! find the nearest kpoint in kpc, and sort eivq according to it

 !call get_freq(q,ndyn,vg,eivq,evq) !replaced by my get_vgr()
 call get_vgr(q,ndyn,vg,eivq,evq)
 omq=eivq(la)
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0
 etacut = eta*3000 ! error due to non-zero eta is 1/9E6

 do ik=1,nkc
    k2(:)=kpc(:,ik)
    ev2=eigenvec(:,:,ik)
    k3 = -q-k2
    call check_inside_bz(k3,g1,g2,g3,inside)

    !call get_freq(k3,ndyn,vg,eiv3,ev3) !replaced by my get_vgr()
    call get_vgr(k3,ndyn,vg,eiv3,ev3)
    do l2=1,ndyn
       om2=eigenval(l2,ik)
       nb2=nbe(om2,temp,classical)
    do l3=1,ndyn

! first check see if the energy is conserved
       om3=eiv3(l3)
       nb3=nbe(om3,temp,classical)
       ires=0; rres=0

       arg=(omega-om3-om2)
       if (abs(arg) .lt. etacut) then
          ires=(nb3+nb2+1)*delta_l(arg,eta)
          rres=(nb3+nb2+1)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega+om3+om2)
       if (abs(arg) .lt. etacut) then
          ires=ires-(nb3+nb2+1)*delta_l(arg,eta)
          rres=rres-(nb3+nb2+1)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega-om3+om2)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb2-nb3)*delta_l(arg,eta)
          rres=rres+(nb2-nb3)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega-om2+om3)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb3-nb2)*delta_l(arg,eta)
          rres=rres+(nb3-nb2)*arg/(arg*arg+eta*eta)
       endif

       self=cmplx(rres,ires*pi)
 !     if (abs(self).lt.1d-6) cycle loop1   ! self is in cm

       call matrix_elt_full(q,k2,k3,omq,om2,om3,evq(:,la),ev2(:,l2),ev3(:,l3),xx)

       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)

       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif

 enddo
 enddo
 enddo
 nself = nself /nkc *pi
 uself = uself /nkc *pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w3

!--------------------------------------------------------------------------------------------------
 subroutine function_self_w4(q,la,omega,temp,nself,uself)
!!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!! this one is for arbitrary q
!! uses gaussian delta with 3 terms, and just calculates the imaginary part.
!! In this subroutine, which calculates like self_w on the fly, the second momentum
!! is restricted to be on the kpc mesh, but the first and third momenta can be anywhere
!! frequencies, temperatures are in cm^-1 (same unit as wmax)
!! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
!! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
!! the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use om_dos
 use phi3
 use eigen

 !use params
 !use constants
 !use io2
 !use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer inside,ik,l3,l2
! integer, save :: cntr
 real(8) omq,eta,term,k3(3),k2(3),eivq(ndyn),eiv2(ndyn),eiv3(ndyn),vg(3,ndyn)
 real(8) etacut,arg,delta_g,om3,om2,nb3,nb2,nbe,ires,rres
 complex(8) omz,self,xx,evq(ndyn,ndyn),ev2(ndyn,ndyn),ev3(ndyn,ndyn)

! be careful band sorting should be avoided in this case; otherwise need to
! find the nearest kpoint in kpc, and sort eivq according to it
 call get_freq(q,ndyn,vg,eivq,evq)
 omq=eivq(la)
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0
 etacut = eta*6 ! error due to non-zero eta is 1/9E6

 do ik=1,nkc
    k2=kpc(:,ik)
    ev2=eigenvec(:,:,ik)
    k3 = -q-k2
    call check_inside_bz(k3,g1,g2,g3,inside)

    call get_freq(k3,ndyn,vg,eiv3,ev3)

    do l2=1,ndyn
       om2=eigenval(l2,ik)
       nb2=nbe(om2,temp,classical)
    do l3=1,ndyn

! first check see if the energy is conserved
       om3=eiv3(l3)
       nb3=nbe(om3,temp,classical)
       ires=0; rres=0

       arg=(omega-om3-om2)
       if (abs(arg) .lt. etacut) then
          ires=(nb3+nb2+1)*delta_g(arg,eta)
!         rres=(nb3+nb2+1)*arg/(arg*arg+eta*eta)
       endif

! this affects the lifetime of low-frequency phonons, can be removed...
! contributes only to noremal processes, and might yield a higher power of w!
!       arg=(omega+om3+om2)
!       if (abs(arg) .lt. etacut) then
!          ires=ires-(nb3+nb2+1)*delta_g(arg,eta)
!!         rres=rres-(nb3+nb2+1)*arg/(arg*arg+eta*eta)
!       endif

       arg=(omega-om3+om2)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb2-nb3)*delta_g(arg,eta)
!         rres=rres+(nb2-nb3)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega-om2+om3)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb3-nb2)*delta_g(arg,eta)
!         rres=rres+(nb3-nb2)*arg/(arg*arg+eta*eta)
       endif

       self=cmplx(rres,ires*pi)
 !     if (abs(self).lt.1d-6) cycle loop1   ! self is in cm

       call matrix_elt_full(q,k2,k3,omq,om2,om3,evq(:,la),ev2(:,l2),ev3(:,l3),xx)

       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)

       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif

 enddo
 enddo
 enddo
 nself = nself /nkc *pi
 uself = uself /nkc *pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w4
!===========================================================
 subroutine calculate_v35(ksubibz,ndn,nib,eivlibz,eivcibz,nk,eival,eivec)
!!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!!all of the write() will give compiling time error, so commented out
!! this version is used when direct access is needed to v33 for given k_indices and bands
!! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
!! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
!! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
!! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
!! between 0 and gi; i=1,2,3 ; units are in cm^-1
!! this works for k1 in the IBZ and k2 in the FBZ on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
!! to avoid ambiguity in the case of degenerate eigenstates, and to satisfy symmetry exactly,
!! not to mention saving time, we calculate only for 0<q_i<G/2 and complete by symmetry
!
 use kpoints
 use phi3
 use eigen
 use atoms_force_constants
 use svd_stuff
 implicit none
 integer, intent(in) :: ksubibz(2),ndn,nk,nib
 integer l1,l2,l3,n2,n3,al,be,ga,i0,j0,k0,j,k,t,ta1,ta2,ta3
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,np1,np2,np3,n1 ,nibzsub
 real(8) q1(3),q2(3),q3(3),mi,mj,mk,rx2(3),rx3(3),den,eivlibz(ndn,nib),eival(ndn,nk)
 complex(8) eivec(ndn,ndn,nk),xx,eiqr,eivcibz(ndn,ndn,nib)
 character(4) cibz1,cibz2
 character(6) cnk
 character(132) filname

INTEGER :: counter = 0
! write(ulog,*)'calculating V3(k2,k3,l1,l2,l3) ##########################'
 write(cibz1,'(i4.4)') ksubibz(1)
 write(cibz2,'(i4.4)') ksubibz(2)
 write(cnk,'(i6.6)') nkc
 nibzsub=ksubibz(2)-ksubibz(1)+1

 !filname=trim(v3path)//'v35-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
 filname='v35-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
 !write(ulog,*)'opening the file ',filname

 !open(uv3,file=trim(filname),status='unknown' ,FORM='UNFORMATTED')
 open(uv3,file=filname,status='unknown')
 write(uv3,*)nkc,ndn,ksubibz(1),ksubibz(2)
!--------------------------------------------------
OPEN(207,file='monitor_v35.dat',status='unknown')

 v3sq = cmplx(0d0,0d0)
! not initializing might be better if some array elements are not properly initialized

 loop1: do n1=ksubibz(1),ksubibz(2)
    np1=n1-ksubibz(1)+1
    q1 = kibz(:,n1) ! =kpc(:,mapinv(n1))
!   call get_k_info_cent(q1,NC,nk1,i2,j2,k2,g1,g2,g3,inside)
    call get_k_info(q1,NC,nk1,i2,j2,k2,g1,g2,g3,inside)
    if(mapinv(n1).ne.nk1) then
!      write(ulog,*)'mapinv(n1),nk1,inside=',mapinv(n1),nk1,inside
!      write(ulog,*)'q1=',q1
!      write(ulog,*)'ijk1=',i2,j2,k2
      stop
    endif
    write(*,2)'nibz,nk,kibz=',n1,nk1,q1
    !write(ulog,2)'nibz,nk,kibz=',n1,nk1,q1
    if(nk1.eq.1) cycle loop1
 loop2: do n2=1,nkc   ! second argument k2 needs to be only in the IFBZ
    q2 = kpc(:,n2)
!   call get_k_info_cent(q2-shft,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    call get_k_info(q2-shft,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(n2.ne.nk2) then
!      write(ulog,*)'n2,nk2,inside=',n2,nk2,inside
!      write(ulog,*)'q2=',q2
!      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif
    if(nk2.eq.1) cycle loop2

    q3=-q1-q2
!   call get_k_info_cent(q3-shft,NC,nk3,i2,j2,k2,g1,g2,g3,inside)
    call get_k_info(q3-shft,NC,nk3,i2,j2,k2,g1,g2,g3,inside)
    if(nk3.eq.1) cycle loop2

!   write(*,2)'nk123=',nk1,nk2,nk3

 do l1=1 ,ndn
 do l2=1 ,ndn

!   np1= l1+(nk1-1)*ndn
!   np2= l2+(nk2-1)*ndn
!   if(np1.gt.np2) cycle    ! permutation symmetry can not be used that easily
! as the thwo k1 and k2 meshes differ. In case they overlap, in that overlap region one can
! use permiation symmetry which can be headache.

 do l3=1 ,ndn
!   np3= l3+(nk3-1)*ndn
!   if(np2.gt.np3) cycle

    xx = cmplx(0d0,0d0)
    tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop

       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j) - atompos(:,j0)
       rx3(:) = atompos(:,k) - atompos(:,k0)
!  be careful: it has to be the translations R not atompos!
!&      exp(ci*(kp(:,n2) .dot. atompos(:,j) + kp(:,n3) .dot. atompos(:,k))) /  &
       eiqr = exp( ci* ((q2 .dot. rx2) + (q3 .dot. rx3)) )

       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
!&      eivec(ta1,l1,nk1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
&      eivcibz(ta1,l1,n1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
! make sure for kp at the FBZ boundary it gives (-1) or 1


! n2 is already in the IBZ, n3 is not but for now we are using the full mesh
! and nibz=nkc;  therefore n1 is also in the full zone and is OK, otherwise
! for n1 and n3 it is necessary to rotate the eigenvector from the IBZ to its
! true position in the recirpocal space. (the op_kmatrix will do the job)
!call naiveCheck(n1,nk2,nk3, l1,l2,l3,ta1,ta2,ta3,t, eiqr,xx, counter)
    enddo tloop


    den = sqrt(eivlibz(l1,n1)*eival(l2,nk2)*eival(l3,nk3))/cnst/sqrt(cnst)

    if (den.ne.0) then
       xx = xx / den/2d0/sqrt(2d0) * const33

!call naiveCheck2(n1,nk2,nk3,l1,l2,l3,xx,den)
!      write(*,6)'nk_i,la_i,xx=',nk1,nk2,nk3,l1,l2,l3,xx
    else   ! usually acoustic eigenvalues at k=0 are not strictly zero, but a small number/original comment
        ! but that small number is 1E-18 which is on the denominator and explodes
!       write(ulog,*)'calculate_v3: denominator=0 !! '
!       write(ulog,*)'n1,n2,n3,l123=', n1,n2,nk3,l1,l2,l3
!       write(ulog,*)'q1=', q1
!       write(ulog,*)'q2=', q2
!       write(ulog,*)'q3=', q3
!       write(ulog,*)'eivs123=',eivlibz(l1,n1),eival(l2,n2),eival(l3,nk3)
       stop
       ! instead of stop, just skip over that den will make it too small

    endif

!notes: v3sq values are stored at this step
!notes: v3sq() is originally  allocated in <read_all_v3sq>

    v3sq(np1,nk2,l1,l2,l3)=xx*conjg(xx)
    write(uv3,*)n1,nk2,l1,l2,l3,v3sq(np1,nk2,l1,l2,l3)

 enddo !l3
 enddo !l2
 enddo !l1
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop2 !n2
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop1 !n1

 close(uv3)

2 format(a,2(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g14.8))
6 format(a,6(i5),9(2x,g14.8))

 end subroutine calculate_v35
!===========================================================
subroutine calculate_w3_ibz_split_sq(ibz_subset,nv3_split,ndn,nk,kp)
!!subroutines extracted from subs.f90
!!used for phonon lifetime calculations
!! This is a modified calculate_w3_ibz by sy.
!! It reads ibz_subset (ibz_subset(1)-starting index, ibz_subset(2)-end index in nibz)
!! This subroutine will not write header in output file. This should be done by another code for merging v33.dat

!! Below is for the original calculate_w3_ibz
!! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
!! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
!! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
!! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
!! between 0 and gi; i=1,2,3 ; units for V3 are in cm^-1
!! this works for two k's (k2,k3) in the prim cell on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
!! this subroutine calculates V3(k1,la1,k2,la2,k3,la3)  with k1=-k2-k3+G  and all ki in the primitive G-cell
!! with the restriction of k2 in the IRREDUCIBLE BZ defined by ki(3,ni) included in kp(3,nk)
!!
 use kpoints
 use phi3
 use atoms_force_constants
 use svd_stuff

 implicit none

 integer nk,ndn,l1,l2,l3,n2,n1,j,k,nibzsub
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,indx
 integer ibz_subset(2),nv3_split
 real(8) q1(3),q2(3),q3(3),eival(ndn,nk) !,sh(3)
 real(8) kp(3,nk)  !ki(3,ni),
 complex(8) eivec(ndn,ndn,nk),xx
 real cputim
 character(4) cibz1,cibz2
 character(6) cnk

! v33 = cmplx(0d0,0d0)
 v33sq=0d0
 indx=0
! sh=(-0.5d0)*(g1+g2+g3) ! this is the shift applied to get kibz
! COMMENT: ulog=30, but its declared both in io2 and ios which is confusing the compiler
 write(30,*)'entered v3 subroutine...'
 write(30,*)'V3: nc(123)=',nc
 write(30,*)'starting ibz point=',ibz_subset(1)
 write(30,*)'ending   ibz point=', ibz_subset(2)
 nibzsub=ibz_subset(2)-ibz_subset(1)+1
 write(30,*)'number of ibz subset points=', nibzsub
 write(cibz1,'(i4.4)') ibz_subset(1)
 write(cibz2,'(i4.4)') ibz_subset(2)
 write(cnk,'(i6.6)') nk

   open(uv3,file='v33-'//cnk//'-'//cibz1//'-'//cibz2//'.dat',status='unknown') !,FORM='UNFORMATTED')
   write(uv3,*)nv3_split,ndn,ibz_subset(1),ibz_subset(2)

!************************************************
! First argument q1 needs to be only in the IFBZ
!************************************************
 loop2: do n1=ibz_subset(1),ibz_subset(2)       ! calculate v33 only for subset of k points
    q1 = kibz(:,n1) !kp(:,mapinv(n2))  ! is also ki(:,n2)  ! should be
    !write(debug,*) 'q2=' q2
    !call get_k_info_cent(q1,NC,nk1,i2,j2,k2,g1,g2,g3,inside)
    call get_k_info(q1,NC,nk1,i2,j2,k2,g1,g2,g3,inside)
     if(mapinv(n1).ne.nk1) then
       write(30,*)'n1,mapinv(n1),nk1,inside=',n1,mapinv(n1),nk1,inside
       write(30,*)'q1=',q1
       write(30,*)'ijk=',i2,j2,k2
       stop
     endif

 call cpu_time(cputim)  !date_and_time(time=tim)
 write(*,1) 'W3_IBZ_SPLIT: nibz,nk1,q...',n1,nk1,q1
 write(30,1) 'entering nibz,nk1,q...',n1,nk1,q1
 write(30,'(a,f12.4)')'  TIME IS ',cputim

    if(nk1.eq.1) cycle loop2
 loop3: do n2=1 ,nk
    q2 = kp(:,n2)   ! third argument is in the whole FBZ coarse mesh
    !write(debug,*) 'q3=', q3
    !call get_k_info_cent(q2-shft,NC,nk2,i3,j3,k3,g1,g2,g3,inside)
    call get_k_info(q2-shft,NC,nk2,i3,j3,k3,g1,g2,g3,inside)
    if(n2.ne.nk2) then
      write(30,*)'n2,nk2,inside=',n2,nk2,inside
      write(30,*)'q2=',q2
      write(30,*)'ijk2=',i3,j3,k3
      stop
    endif
    if(nk2.eq.1) cycle loop3
    q3 = -q2-q1

    !call get_k_info_cent(q3-shft,NC,nk3,i1,j1,k1,g1,g2,g3,inside)
    call get_k_info(q3-shft,NC,nk3,i1,j1,k1,g1,g2,g3,inside)
    if(nk3.eq.1) cycle loop3
 ! write(ulog,2)'nk123=',nk1,nk2,nk3
  !write(ulog,*)'q2,q3,q1',q2,q3,q1
 do l1=1 ,ndn
 do l2=1 ,ndn
 do l3=1 ,ndn
    call matrix_elt(q1,q2,l1,l2,l3,xx,inside)
    indx=indx+1
!    v33(indx)=xx
    v33sq(indx)=xx*conjg(xx)

!   write(uv3,5) n1,n2,l1,l2,l3,v33sq(indx)

     nq1(indx)= (nk1-1)*ndn+l1
     nq2(indx)= (nk2-1)*ndn+l2
     nq3(indx)= (nk3-1)*ndn+l3

     !note: temporarily commented out for saving space on Rivanna
     !write(uv3,5) nq1(indx),nq2(indx),nq3(indx),v33sq(indx)
! reverse mapping is: la=mod(nq,ndn) ; nk=1+(nq-la)/ndn
!   nq2(indx)= nk2
!   nq3(indx)= nk3
!   la1(indx)= l1
!   la2(indx)= l2
!   la3(indx)= l3
    !write(debug,*),'l1,l3,l2,indx,v33',l1,l3,l2,indx,v33(indx)
!   if (mod(indx,1000).eq.1) write(*,*) 'indx,v33=',indx,xx
 enddo
 enddo
 enddo
 enddo loop3
 enddo loop2

 nv3 = indx
 write(30,*)' V33: total size of this array is =',nv3

! if( writev3.eq.1) then
! if the v33 are calculated, they are automatically written into a file
   !write(uv3,*) nv3_split
!   do j=1 ,nv3_split
!   enddo
   close(uv3)
! endif

1 format(a,2(1x,i5),9(1x,g12.5))
2 format(a,3(3x,i5),9(2x,g12.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g12.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(3(i6),2x,g14.7) ! original one doesn't match
6 format(a,6(i5),9(2x,g15.8))
7 format(6(i6),9(2x,g15.8))

 end subroutine calculate_w3_ibz_split_sq
!==========================================================================
!--------------------------------------------------------------------------
    ! check subroutine rta_sub(ksub,tempk)

!=================================All below are for checking purpose===================================================================
SUBROUTINE naiveCheck(n1,nk2,nk3,l1,l2,l3,ta1,ta2,ta3,t, eiqr,xx, counter)
!!check subroutine, not used
    USE eigen
    USE svd_stuff
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n1, nk2, nk3, l1, l2, l3,ta1,ta2,ta3,t
    INTEGER, INTENT(inout) :: counter
    COMPLEX(8),INTENT(in) :: eiqr,xx
    IF(n1.eq.2 .and. nk2.eq.3 .and.l1.eq.1 .and. l2.eq.4 .and. l3.eq.1) THEN
        counter = counter + 1
        WRITE(207,*) '==================================================='
        WRITE(207,7) 'n1=',n1,'nk2=',nk2,'nk3=',nk3
        WRITE(207,7) 'l1=', l1, 'l2=', l2, 'l3=', l3
        WRITE(207,*) 'eiqr=', eiqr
        WRITE(207,7) 'tau1=',ta1,'tau2=',ta2,'tau3=',ta3
        WRITE(207,*) 'eivcibz(ta1, l1, n1)=',eivecibz(ta1,l1,n1)
        WRITE(207,*) 'eivec(ta2,l2,nk2)=',eigenvec(ta2,l2,nk2)
        WRITE(207,*) 'eivec(ta3,l3,nk3)=',eigenvec(ta3,l3,nk3)

        WRITE(207,*)'!!!!!!!!!!!!!ei*ei*ei=',eivecibz(ta1,l1,n1)*&
        &eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)

        WRITE(207,*) 'fcs_3(t)=', fcs_3(t)
        WRITE(207,*) 'xx before den is:', xx
        WRITE(207,*) 'This is the', counter, 'th item in the sum.'
        WRITE(207,*)
    END IF
7 format(3(a,i5,3x))
END SUBROUTINE naiveCheck
!----------------------------
SUBROUTINE naiveCheck2(n1,nk2,nk3,l1,l2,l3,xx, den)
!!check subroutine, not used
    USE eigen
    USE svd_stuff
    USE constants
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n1, nk2, nk3, l1, l2, l3
    COMPLEX(8),INTENT(in) :: xx
    REAL(8),INTENT(in) :: den
    IF(n1.eq.2 .and. nk2.eq.3 .and.l1.eq.1 .and. l2.eq.4 .and. l3.eq.1) THEN
        WRITE(207,*) '***************************************************'
        WRITE(207,*) '***************************************************'
        WRITE(207,*) 'eivlibz(l1,n1)=', eivalibz(l1, n1)
        WRITE(207,*) 'eival(l2,nk2)=', eigenval(l2, nk2)
        WRITE(207,*) 'eival(l3,nk3)=', eigenval(l3, nk3)
        WRITE(207, *)'cnst=', cnst
        WRITE(207,*) ' den = ', den
        WRITE(207,*) ' xx now  is:', xx
        WRITE(207,*) '***************************************************'
        WRITE(207,*) '***************************************************'
    END IF
END SUBROUTINE naiveCheck2
!SUBROUTINE CheckYY(Y_square,flag)
!    IMPLICIT NONE
!    REAL(8),DIMENSION(d,atom_number,d,tot_atom_number),INTENT(IN) :: Y_square
!    LOGICAL, INTENT(OUT) :: flag
!
!    INTEGER :: i,atom1,atom2,direction1,direction2
!    REAL(8) :: difference
!    INTEGER :: unitnumber
!unitnumber=75
!OPEN(unitnumber,FILE='YY_Compare.dat',STATUS='unknown',POSITION='append',ACTION='write')
!
!    flag=.FALSE.
!    DO i=1,ifc2_terms
!        atom1=indiefc2_index(i)%iatom_number
!        atom2=indiefc2_index(i)%jatom_number
!        direction1=indiefc2_index(i)%iatom_xyz
!        direction2=indiefc2_index(i)%jatom_xyz
!
!        difference=ABS(Y_square(direction1,atom1,direction2,atom2)-YY_Compare(atom1,atom2)%phi(direction1,direction2))
!WRITE(unitnumber,'(4(1X,A10,1X,I3))')'direction1',direction1,'atom1',atom1,&
!& 'direction2',direction2,'atom2',atom2
!WRITE(unitnumber,'(A10,2X,F10.8)') 'YY Change:',difference
!        !IF(difference.le.(danger*1d3).AND.(i.ne.ifc2_terms)) THEN
!        !    CYCLE
!        !ELSEIF(difference.le.(danger*1d3).AND.(i.eq.ifc2_terms)) THEN
!        !    flag=.TRUE.
!        !ELSE
!        !    EXIT
!        !END IF
!
!        IF(difference.gt.(danger*1d3)) THEN
!            EXIT
!        ELSE
!            IF(i.eq.ifc2_terms) THEN
!                flag=.TRUE.
!                EXIT
!            END IF
!            CYCLE
!        END IF
!    END DO
!
!    !remember to update YY_Compare, it cannot be in the same loop for comparison
!    DO i=1,ifc2_terms
!        atom1=indiefc2_index(i)%iatom_number
!        atom2=indiefc2_index(i)%jatom_number
!        direction1=indiefc2_index(i)%iatom_xyz
!        direction2=indiefc2_index(i)%jatom_xyz
!
!        YY_Compare(atom1,atom2)%phi(direction1,direction2)=Y_square(direction1,atom1,direction2,atom2)
!    END DO
!
!CLOSE(unitnumber)
!
!END SUBROUTINE CheckYY
!------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE CheckASR
!    IMPLICIT NONE
!    INTEGER :: direction1, direction2
!    INTEGER :: unitnumber
!    INTEGER :: atom1,atom2
!    REAL(8) :: check
!
!    unitnumber=20
!    OPEN(unitnumber,FILE='CheckASR.dat',STATUS='unknown',POSITION='append',ACTION='write')
!
!    DO atom1=1,atom_number
!    DO direction1=1,d
!    DO direction2=1,d
!        check=SUM(trialfc2_value(atom1,:)%phi(direction1,direction2))
!        !IF(ABS(check).gt.1d-10) THEN
!        WRITE(unitnumber,'(3(X,I3))') atom1,direction1,direction2
!        WRITE(unitnumber,*)'ASR Check: ',check
!        !END IF
!    END DO
!    END DO
!    END DO
!
!    WRITE(unitnumber,*)
!    DO atom1=1,atom_number
!    DO atom2=1,atom_number
!    DO direction1=1,d
!    DO direction2=1,d
!
!    IF(trialfc2_value(atom1,atom2)%phi(direction1,direction2).ne.myfc2_value(atom1,atom2)%phi(direction1,direction2)) THEN
!        WRITE(unitnumber,*)'Found Difference: atom1, atom2, direction1, direction2'
!        WRITE(unitnumber,'(4I4)')atom1,atom2,direction1,direction2
!        WRITE(unitnumber,*)'Original fc2 value:',myfc2_value(atom1,atom2)%phi(direction1,direction2)
!        WRITE(unitnumber,*)'Current fc2 value:',trialfc2_value(atom1,atom2)%phi(direction1,direction2)
!        WRITE(unitnumber,*)
!    END IF
!
!    END DO
!    END DO
!    END DO
!    END DO
!    WRITE(unitnumber,*)'=============================================================='
!    CLOSE(unitnumber)
!END SUBROUTINE CheckASR
!!=========================================================================================================================
SUBROUTINE FixASR2
!!Force ASR in FC2
    IMPLICIT NONE
    INTEGER :: direction1, direction2
    INTEGER :: unitnumber
    INTEGER :: atom1
    REAL(8) :: check

    unitnumber=21
    OPEN(unitnumber,FILE='FixASR2.dat',STATUS='unknown',POSITION='append',ACTION='write')

    DO atom1=1,atom_number
    DO direction1=1,d
    DO direction2=1,d
        check=SUM(trialfc2_value(atom1,:)%phi(direction1,direction2))
        IF(ABS(check).gt.danger) THEN
            trialfc2_value(atom1,atom1)%phi(direction1,direction2) = &
            &trialfc2_value(atom1,atom1)%phi(direction1,direction2) - check
        WRITE(unitnumber,*)'ASR Fix by: ',check
        END IF
    END DO
    END DO
    END DO

    WRITE(unitnumber,*)'=============================================================='
    CLOSE(unitnumber)
END SUBROUTINE FixASR2
!!=========================================================================================================================
!SUBROUTINE CheckASR3
!    IMPLICIT NONE
!    INTEGER direction1, direction2, direction3
!    INTEGER :: unitnumber
!    INTEGER :: atom1,atom3,new
!    REAL(8) :: check
!
!    unitnumber=22
!    OPEN(unitnumber,FILE='CheckASR3.dat',STATUS='unknown',POSITION='append',ACTION='write')
!
!    DO atom1=1,atom_number
!    DO atom3=1,tot_atom_number
!        IF (.NOT.ANY(fc3_unique_idx==atom3)) THEN
!            CYCLE
!        ELSE
!            new = find_loc(fc3_unique_idx, new)
!        END IF
!    DO direction1=1,d
!    DO direction2=1,d
!    DO direction3=1,d
!        check=SUM(myfc3_value(atom1,new,:)%psi(direction1,direction2,direction3))
!        IF(ABS(check).gt.danger) THEN
!            WRITE(unitnumber,*) 'ASR not satisfied'
!            WRITE(unitnumber,'(5I3)') atom1, new, direction1, direction2, direction3
!            WRITE(unitnumber,*) check
!        END IF
!    END DO
!    END DO
!    END DO
!    END DO
!    END DO
!
!    CLOSE(unitnumber)
!END SUBROUTINE CheckASR3
!!=========================================================================================================================
!SUBROUTINE CheckASR4
!    IMPLICIT NONE
!    INTEGER :: direction1, direction2,direction3,direction4
!    INTEGER :: unitnumber
!    INTEGER :: atom1,atom3,atom4,new3,new4
!    REAL(8) :: check
!
!    unitnumber=23
!    OPEN(unitnumber,FILE='CheckASR4.dat',STATUS='unknown',POSITION='append',ACTION='write')
!
!    DO atom1=1,atom_number
!    DO atom3=1,tot_atom_number
!        IF(.NOT.ANY(fc4_unique_idx==atom3)) THEN
!            CYCLE
!        ELSE
!            new3 = find_loc(fc4_unique_idx,atom3)
!        END IF
!    DO atom4=1,tot_atom_number
!        IF(.NOT.ANY(fc4_unique_idx==atom4)) THEN
!            CYCLE
!        ELSE
!            new4 = find_loc(fc4_unique_idx, atom4)
!        END IF
!    DO direction1=1,d
!    DO direction2=1,d
!    DO direction3=1,d
!    DO direction4=1,d
!        check=SUM(myfc4_value(atom1,:,new3,new4)%chi(direction1,direction2,direction3,direction4))
!        IF(ABS(check).ne.0d0) THEN
!            WRITE(unitnumber,*) 'ASR not satisfied'
!            WRITE(unitnumber,'(3I3)') atom1, new3, new4
!            WRITE(unitnumber,'(4I3)') direction1, direction2, direction3, direction4
!            WRITE(unitnumber,*) check
!        END IF
!    END DO
!    END DO
!    END DO
!    END DO
!    END DO
!    END DO
!    END DO
!
!WRITE(unitnumber,*) '================================================================================'
!    CLOSE(unitnumber)
!
!END SUBROUTINE CheckASR4
!--------------------------------------------------------------------------------------------------
    SUBROUTINE sym_strain
      !!force symmetry on strain gradients
        IMPLICIT NONE
        INTEGER :: i,j
        DO i=1,d-1
            DO j=i+1,d
                GradientV_eta(i,j) = GradientV_eta(j,i)
                strain(i,j) = strain(j,i)
            END DO
        END DO
    END SUBROUTINE sym_strain
!--------------------------------------------------------------------------------------------------
     SUBROUTINE R_Update
     !!update translational vectors with strain eta
        IMPLICIT NONE
        INTEGER :: i=0
        !restore original
        cell_vec=cell_save
        trans_vec=trans_save

        !update
        UpdateTranslationalVector: DO i=1,d
                                       trans_vec(:,i)=trans_vec(:,i)+(strain(:,:).dot.trans_vec(:,i))
                                  END DO UpdateTranslationalVector

        r1%component=trans_vec(:,1)
        r2%component=trans_vec(:,2)
        r3%component=trans_vec(:,3)

        call make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
        call make_reciprocal_lattice(g1,g2,g3,rr1,rr2,rr3)
        call calculate_volume(r1,r2,r3,volume_r)
        call calculate_volume(g1,g2,g3,volume_g)

!        UpdateCell: DO i=1,cell_number+1
!                        cell_vec(:,i)=cell_vec(:,i)+(strain(:,:).dot.cell_vec(:,i))
!                    END DO UpdateCell
        !no need
        !UpdateAtom:DO i=1,tot_atom_number
                     !every_atom(i)%R = cell_vec(:,every_atom(i)%type_R)
                   !End DO UpdateAtom


    END SUBROUTINE R_Update
!--------------------------------------------------------------------------------------------------
    SUBROUTINE TreatGradientV
    !!adjust <V> gradients of utau and eta by a matrix
    !!this subroutine is used for the guessing start
        IMPLICIT NONE
        INTEGER :: i,j,voigt_idx1,voigt_idx2
        INTEGER :: atom1,atom2,xyz1,xyz2
        INTEGER,DIMENSION(2) :: temp_xyz
        REAL(8) :: cell_volume
        REAL(8),DIMENSION(6,6) :: S
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: temp

        !modify GradientV_utau
        ALLOCATE(temp(d,atom_number))
        temp = 0d0
        j = 0
        DO i=1,d*atom_number
            IF(ANY(i.eq.fixed_params)) THEN
               CYCLE !means that variable i is fixed
            ELSE
                j = j + 1 !found a free 'utau'
                atom1 = INT((i+d-1)/d)
                xyz1 = MOD((i+d-1),d)+1
                DO atom2=1,atom_number
                DO xyz2=1,d
                    temp(xyz1,atom1) = temp(xyz1,atom1)+inv_phiTau(atom1,atom2)%phi(xyz1,xyz2)*&
                    &GradientV_utau(xyz2,atom2)
                END DO !xyz2
                END DO !atom2
            END IF
        END DO
        GradientV_utau = temp !modify GradientV_utau
        DEALLOCATE(temp)

        !modify GradientV_eta
        S=compliance
        CALL calculate_volume(r1,r2,r3,cell_volume)
        ALLOCATE(temp(d,d))
        temp = 0d0
        DO i=d*atom_number+1,d*atom_number+9
            IF(ANY(i.eq.fixed_params)) THEN
                CYCLE !means that variable i is fixed
            ELSE
                j = j + 1 !found a free 'eta'(strain)
                voigt_idx1 = voigtMap(i-d*atom_number)
                temp_xyz = inv_voigtMap(voigt_idx1)
                DO xyz1=1,d
                DO xyz2=1,d
                    voigt_idx2 = voigtMap((xyz1-1)*d+xyz2)
                    temp(temp_xyz(1),temp_xyz(2)) = temp(temp_xyz(1),temp_xyz(2)) + &
                    & S(voigt_idx1,voigt_idx2)*GradientV_eta(xyz1,xyz2)/cell_volume
                END DO
                END DO

            END IF
        END DO

        DO xyz2=1,d-1
        DO xyz1=xyz2+1,d
            temp(xyz1,xyz2) = temp(xyz2,xyz1)
        END DO
        END DO
        GradientV_eta=temp !modify GradientV_eta
        DEALLOCATE(temp)
    END SUBROUTINE TreatGradientV
!--------------------------------------------------------------------------------------------------
    SUBROUTINE SimpleUpdate(mx,x,f)
     !!just simple update,same as broy, it updates x, not used
        IMPLICIT NONE
        INTEGER :: i,j,temp,voigt_idx1,voigt_idx2
        INTEGER :: atom1,atom2, xyz1,xyz2
        INTEGER,INTENT(in) :: mx
        REAL(8),INTENT(INOUT),DIMENSION(:) :: x
        REAL(8),INTENT(in),DIMENSION(:) :: f
        REAL(8) :: update !for the temporal dotproduct
        REAL(8) :: cell_volume !for unit cell volume
        REAL(8) :: unit_conv !to convert the unit
        !--------------------------------------
        REAL(8),DIMENSION(6,6) :: S

        !way1: hard coded
        S(:,:)=0d0
        S(1,1)=7.68;S(2,2)=S(1,1);S(3,3)=S(1,1)
        S(1,2)=-2.14;S(2,1)=S(1,2);S(3,1)=S(1,2);S(1,3)=S(1,2);S(2,3)=S(1,2);S(3,2)=S(1,2)
        S(4,4)=12.6;S(5,5)=S(4,4);S(6,6)=S(4,4)
        unit_conv = 0.16!1d-12*1d30*1d-19
        !way2: use my own calculated C2
        S=compliance
        unit_conv=1
        !--------------------------------------
        update = 0d0
        !temporary, only works for the special case
        j = 0
        DO i=1,d*atom_number
            IF(ANY(i.eq.fixed_params)) THEN
               CYCLE !means that variable i is fixed
            ELSE
                j = j + 1 !found a free 'utau'

                atom1 = INT((i+d-1)/d)
                xyz1 = MOD((i+d-1),d)+1
                update = 0d0
                DO atom2=1,atom_number
                DO xyz2=1,d
                    update = update+inv_phiTau(atom1,atom2)%phi(xyz1,xyz2)*&
                    &GradientV_utau(xyz2,atom2)
                END DO !xyz2
                END DO !atom2
                x(j) = x(j) - update
            END IF
        END DO

        CALL calculate_volume(r1,r2,r3,cell_volume)
        update = 0d0
        DO i=d*atom_number+1,d*atom_number+9
            IF(ANY(i.eq.fixed_params)) THEN
                CYCLE !means that variable i is fixed
            ELSE
                j = j + 1 !found a free 'eta'(strain)

                temp = i-d*atom_number
                voigt_idx1 = voigtMap(temp)
                update = 0d0
                DO xyz1=1,d
                DO xyz2=1,d
                    voigt_idx2 = voigtMap((xyz1-1)*d+xyz2)
                    update = update + S(voigt_idx1,voigt_idx2)*&
                    &GradientV_eta(xyz1,xyz2)/cell_volume
                END DO
                END DO

                x(j) = x(j) - update*unit_conv
            END IF
        END DO

        DO i=d*atom_number+10,d*atom_number+9+eff_fc2_terms
            IF(ANY(i.eq.fixed_params)) THEN
                CYCLE !means that variable i is fixed
            ELSE
                j = j + 1 !found a free eff_fc2
                x(j) = x(j)+2*f(j)
            END IF
        END DO

    END SUBROUTINE SimpleUpdate
!--------------------------------------------------------------------------------------------
    SUBROUTINE updateK
     !!update all K, using new strain and atomic deviation, but no <yy>
 !!used for specific testing cases, should not be used in real calculations
        IMPLICIT NONE

        INTEGER :: i
        INTEGER :: atom1, atom2, xyz1, xyz2
        INTEGER :: atom3, atom4
        INTEGER :: new_atom2, new_atom3,new_atom4

        REAL(8),DIMENSION(:,:),ALLOCATABLE :: S

        ALLOCATE(S(d,tot_atom_number))
        DO i=1,tot_atom_number
            S(:,i) = (strain(:,:).dot.(every_atom(i)%R+every_atom(i)%tau)) + atomic_deviation(:,every_atom(i)%type_tau)
        END DO

        !quadratic term Phi
        DO atom1=1,atom_number
        DO atom2=1,tot_atom_number
            trialfc2_value(atom1,atom2)%phi = myfc2_value(atom1,atom2)%phi
        END DO
        END DO

        !cubic term Psi*S
        DO atom1=1,atom_number

        DO atom2=1,tot_atom_number
            IF(.NOT.ANY(fc3_unique_idx==atom2)) CYCLE

        DO atom3=1,tot_atom_number
            IF(.NOT.ANY(fc3_unique_idx==atom3)) CYCLE

            new_atom2 = find_loc(fc3_unique_idx,atom2)
            new_atom3 = find_loc(fc3_unique_idx,atom3)

            trialfc2_value(atom1,atom2)%phi(:,:) = trialfc2_value(atom1,atom2)%phi(:,:) + &
                & (myfc3_value(atom1,new_atom2,new_atom3)%psi(:,:,:).dot.S(:,atom3))

        END DO !atom3 loop
        END DO !atom2 loop
        END DO !atom1 loop

        !quartic term 0.5*Chi*SS
        DO atom1=1,atom_number
        DO atom2=1,tot_atom_number
            IF(.NOT.ANY(fc4_unique_idx==atom2)) CYCLE !skip if j not valid
        DO atom3=1,tot_atom_number
            IF(.NOT.ANY(fc4_unique_idx==atom3)) CYCLE !skip if m not valid
        DO atom4=1,tot_atom_number
            IF(.NOT.ANY(fc4_unique_idx==atom4)) CYCLE !skip if n not valid

            new_atom2 = find_loc(fc4_unique_idx,atom2)
            new_atom3 = find_loc(fc4_unique_idx,atom3)
            new_atom4 = find_loc(fc4_unique_idx,atom4)


            trialfc2_value(atom1,atom2)%phi(:,:) = trialfc2_value(atom1,atom2)%phi(:,:) + &
                & 0.5*(myfc4_value(atom1,new_atom2,new_atom3,new_atom4)%chi(:,:,:,:)&
                &.dot.(S(:,atom4)-S(:,atom1)).dot.(S(:,atom3)-S(:,atom1)))
        END DO !atom4 loop
        END DO !atom3 loop
        END DO !atom2 loop
        END DO !atom1 loop
        trialfc2_value_initial = trialfc2_value
        DEALLOCATE (S)

    END SUBROUTINE updateK
!--------------------------------------------------------------------------------------------
    SUBROUTINE updateFCS2
    !!copy the final trialfc2 values to fcs_2
        IMPLICIT NONE
        INTEGER :: atom1, atom2, xyz1, xyz2
        INTEGER :: i

        DO i=1,SIZE(fcs_2)
            atom1 = iatomterm_2(1,i)
            atom2 = iatomterm_2(2,i)
            xyz1 = ixyzterm_2(1,i)
            xyz2 = ixyzterm_2(2,i)
            IF(atom1.le.atom_number) THEN
                fcs_2(i) = trialfc2_value(atom1,atom2)%phi(xyz1,xyz2)
            ELSEIF(atom1.gt.atom_number .AND. atom2.le.atom_number) THEN
                fcs_2(i) = trialfc2_value(atom2,atom1)%phi(xyz2,xyz1)
            ENDIF
        END DO

    END SUBROUTINE updateFCS2
!--------------------------------------------------------------------------------------------
    FUNCTION normGradients(selected,x,f) RESULT(norm)
    !!calculate the weighted L1 norm of gradients f(:)
    !!used for Broyden threshold check
        IMPLICIT NONE
        INTEGER,DIMENSION(:),INTENT(in) :: selected
        REAL(8),DIMENSION(:),INTENT(in) :: x,f
        REAL(8) :: norm
        INTEGER :: i,j,baffle

        !calculate atomic deviation & strain related gradients norm
        baffle=0
        DO i=1, SIZE(selected)
           IF(selected(i).gt.(variational_parameters_size(1)+variational_parameters_size(2))) THEN
                EXIT
           END IF
           baffle=baffle+1
        END DO

        norm = 0d0
        DO i=1,baffle
            norm=norm+ABS(f(i)*x(i))
        END DO

        !calculate <YY> related gradients norm
        DO i=baffle+1,SIZE(x)
            j = selected(i) - variational_parameters_size(1) - variational_parameters_size(2)
            norm=norm+ABS(f(i)*trialfc2_record(j))
        END DO

    END FUNCTION normGradients
!----------------------------------------------------------------------------------------------------
    SUBROUTINE printGradients(x)
    !!output all gradients and variational parameters in a file 'GradientF.dat'
    !!for every iteration
        IMPLICIT NONE
        INTEGER :: checkF, i, j
        INTEGER :: atom1,xyz1,atom2,xyz2
        REAL(8),DIMENSION(:),INTENT(in) :: x!all variational parameters
        CHARACTER(10) :: flag

        checkF = 19
        OPEN(checkF,FILE='GradientF.dat',STATUS='unknown',ACTION='write',POSITION='append')
        WRITE(checkF,*)'current interation #:',iter_rec
        WRITE(checkF,*)'temperature=',temperature*c_light*h_plank*100/k_b
        WRITE(checkF,*)'F0=', F0
        WRITE(checkF,*)'V0=', V0
        WRITE(checkF,*)'free energy=',REAL(F_trial)
        WRITE(checkF,*)'=============GradientF:trial fc2===================='

        WRITE(checkF,*)'largest gradient= ',MAXVAL(ABS(GradientF_trial))

        WRITE(checkF,*) '||Atomic deviation u_tau(:)||'
        DO i=1, variational_parameters_size(1)
            IF(ANY(fixed_params.eq.i)) THEN
                flag = '  FIXED'
            ELSE
                flag = '  FREE'
            END IF
            WRITE(checkF,*) (i-1)/3+1, get_letter(MOD(i-1,3)+1),&
            &' variable=',x(i),' gradient=',REAL(GradientF_trial(i)),flag
        END DO

        WRITE(checkF,*) '||Strain Tensor||'
        DO i=variational_parameters_size(1)+1,variational_parameters_size(1)+variational_parameters_size(2)
            IF(ANY(fixed_params.eq.i)) THEN
                flag = '  FIXED'
            ELSE
                flag = '  FREE'
            END IF
            j = i-variational_parameters_size(1)
            WRITE(checkF,*) get_letter((j-1)/3+1), get_letter(MOD(j-1,3)+1),&
            &' variable=',x(i),' gradient=',REAL(GradientF_trial(i)),flag
        END DO

        WRITE(checkF,*) '||Force Constants||'
        DO i=variational_parameters_size(1)+variational_parameters_size(2)+1,SUM(variational_parameters_size)
            IF(ANY(fixed_params.eq.i)) THEN
                flag = '  FIXED'
            ELSE
                flag = '  FREE'
            END IF

            j = i-variational_parameters_size(1)-variational_parameters_size(2)
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^old^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            atom1 = eff_fc2_index(j)%iatom_number
            xyz1 = eff_fc2_index(j)%iatom_xyz
            atom2 = eff_fc2_index(j)%jatom_number
            xyz2 = eff_fc2_index(j)%jatom_xyz
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            WRITE(checkF,*) atom1,get_letter(xyz1),atom2,get_letter(xyz2),&
            &' variable=',x(i),' gradient=',REAL(GradientF_trial(i)),flag
        END DO

        WRITE(checkF,*)
        CLOSE(checkF)
    END SUBROUTINE printGradients
!-------------------------------------------------------------------------------------------------------
    SUBROUTINE checkGradientYY
    !!check subroutine, not used
        IMPLICIT NONE
        INTEGER :: checkYY,rnk,ntindp
        INTEGER :: i,j
        INTEGER :: atom1, atom2, direction1, direction2
        TYPE(fc2_value),DIMENSION(:),ALLOCATABLE :: checkASR

        ALLOCATE(checkASR(atom_number))
        DO i=1, atom_number
            checkASR(i)%phi = 0d0
        END DO

        checkYY = 20
        OPEN(checkYY,FILE='GradientYY.dat',STATUS='unknown',ACTION='write',POSITION='append')
        WRITE(checkYY,*)'==========================================='
        rnk = 2
        DO j=1, map(rnk)%ngr
            ntindp = map(rnk)%ntind(j)
            WRITE(checkYY,*)'----group',j,'-----'
            DO i=1, map(rnk)%nt(j)
                atom1 = map(rnk)%gr(j)%iat(1,i)
                IF(atom1.gt.atom_number) CYCLE
                atom2 = map(rnk)%gr(j)%iat(2,i)
                direction1 = map(rnk)%gr(j)%ixyz(1,i)
                direction2 = map(rnk)%gr(j)%ixyz(2,i)
                WRITE(checkYY,'(f8.5,f7.3)') GradientV_cor(atom1,atom2)%phi(direction1,direction2),&
                &map(rnk)%gr(j)%mat(i,1:ntindp)
!                checkASR(atom1)%phi(direction1,direction2) = checkASR(atom1)%phi(direction1,direction2)&
!                &+GradientV_cor(atom1,atom2)%phi(direction1,direction2)
            END DO
        END DO

        DO i=1,atom_number
        DO direction1=1,d
        DO direction2=1,d
        DO j=1, tot_atom_number
            checkASR(i)%phi(direction1,direction2) = checkASR(i)%phi(direction1,direction2)&
            &+GradientV_cor(i,j)%phi(direction1,direction2)
        END DO
        IF(checkASR(i)%phi(direction1,direction2).ne.0) THEN
            WRITE(checkYY,*)'ASR not satisfied'
            WRITE(checkYY,*)i,get_letter(direction1),get_letter(direction2)
            WRITE(checkYY,*)checkASR(i)%phi(direction1,direction2)
        END IF
        END DO
        END DO
        END DO

        CLOSE(checkYY)
    END SUBROUTINE checkGradientYY
!-------------------------------------------------------------------------------------------------------
    SUBROUTINE printResults
    !!major subroutine that output the targetinitialization.dat file
    !!that file has the optimized/converged/solved variational parameters x(:)
    !!that file can also be used as an input for 'inherited' start, where user can assign start value for specific vari para
        implicit none

        INTEGER :: i,j
        INTEGER :: atom1,atom2,xyz1,xyz2

        OPEN(21,file='targetInitialize.dat',status='unknown',action='write')
        !output utau, one vector/atom per line
        DO i=1,atom_number
            WRITE(21,*) atomic_deviation(:,i)
        END DO
        !output strain,3x3
        DO i=1,3
            WRITE(21,*) strain(:,i)
        END DO
        !output all the trialfc2
        DO i=1,SIZE(myfc2_index)
            atom1 = myfc2_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2 = myfc2_index(i)%jatom_number
            xyz1 = myfc2_index(i)%iatom_xyz
            xyz2 = myfc2_index(i)%jatom_xyz
            WRITE(21,*)i,atom1,atom2,xyz1,xyz2,trialfc2_value(atom1,atom2)%phi(xyz1,xyz2)
        END DO
        CLOSE(21)
    END SUBROUTINE printResults
!===============================POST PROCESS RELATED================================
    SUBROUTINE my_set_omdos(lambda,mesh)
     !!my subroutine to set om(i), for phonon DOS calculation
    !!wmax is not read from params.phon but the actual max eivals
        IMPLICIT NONE
        INTEGER :: i
        INTEGER, INTENT(in) :: lambda, mesh

        IF(.NOT.ALLOCATED(om)) ALLOCATE(om(mesh))
        IF(.NOT.ALLOCATED(dos)) ALLOCATE(dos(0:lambda,mesh))

        ndyn2 = lambda

!        wmax = MAXVAL(MAXVAL(eigenval,DIM=2))
        wmax = MAXVAL(MAXVAL(eivals,DIM=2))
        wmax = SQRT(wmax)

        DO i=1,mesh
            om(i)=wmax*(0.0001 + (i-1)/(mesh-1.)) !these om(i)*cnst should give the regular cm^-1 phonons
        END DO

    END SUBROUTINE my_set_omdos
!-------------------------------------------------------------------------------------------------------
    SUBROUTINE calc_dos_gauss
    !! calculate phonon dos using gauss broadening, params defined in params.phon
        IMPLICIT NONE
        INTEGER :: lambda,i,j
        REAL(8),ALLOCATABLE,DIMENSION(:) :: evc
        CHARACTER :: wxt*3

        WRITE(*,*) 'running gaussian dos'

        OPEN(udos,FILE='dos_gauss.dat',STATUS='unknown')

!        CALL set_omdos(ndyn,wmesh)  ! allocates om(:) and dos(:,:) arrays of size mesh
        CALL my_set_omdos(ndyn,wmesh) !should I use my own ?

        !only IBZ
        ALLOCATE(evc(nibz))
        dos = 0d0
        DO lambda=1,ndyn
            j=0
            DO i=1,nibz
                j=j+1
!                evc(j)=eigenval(lambda,mapinv(i))!!!!!!!!!!!!!!!!
                evc(j)=SQRT(eivals(lambda,mapinv(i)))
            END DO
            CALL calculate_dos(nibz,evc,wibz,wmesh,om,dos(lambda,:))
            WRITE(wxt,'(i2.2)')lambda
            dos(0,:) = dos(0,:) + dos(lambda,:)
        END DO
        CALL write_dos
        WRITE(*,*) 'gaussian dos done'

        !all BZ
!        ALLOCATE(evc(nkc))
!        dos = 0d0
!        DO lambda=1,ndyn
!            j=0
!            DO i=1,nkc
!                j=j+1
!!                evc(j)=eigenval(lambda,mapinv(i))!!!!!!!!!!!!!!!!
!                evc(j)=SQRT(eivals(lambda,mapinv(i)))
!            END DO
!            CALL calculate_dos(nkc,evc,wk,wmesh,om,dos(lambda,:))
!            WRITE(wxt,'(i2.2)')lambda
!            dos(0,:) = dos(0,:) + dos(lambda,:)
!        END DO
!        CALL write_dos
!        WRITE(*,*) 'gaussian dos done'


        DEALLOCATE(dos)
        DEALLOCATE(evc)
        CLOSE(udos)

    END SUBROUTINE calc_dos_gauss
    
     subroutine calculate_dos(mx,eival,wkp,mesh,omega,ds)
     !! use gaussian broadening to calculate phonon DOS
         implicit none
         integer i,j,mx,mesh
         real(8) x,wkp(mx),delta,ds(mesh),omega(mesh),eival(mx)  !,cnst

        ! wkp=1d0/(nkx*nky*nkz)
        ! write(udos,*)'# wkp,width =',wkp,width

            do i=1,mesh
               ds(i) = 0
               do j=1,mx
                  x = (eival(j) - omega(i))/width *cnst   ! these are all in cm^-1
        !         x = (eival(j) - ene*ene)/width/width
                  if ( abs(x) .gt. 5 ) cycle
                  ds(i) = ds(i) + delta(x)/width*wkp(j) !/mx
               enddo
            enddo

3 format(i5,9(3x,g11.5))

     end subroutine calculate_dos
!-------------------------------------------------------------------------------------------------------
    SUBROUTINE calc_dos_tet
    !! use tetrahedron to calculate phonon DOS
        IMPLICIT NONE

        INTEGER :: lambda,i
        REAL(8) :: wmin,res
        REAL(8),ALLOCATABLE,DIMENSION(:) :: evc,lg
        REAL(8),ALLOCATABLE,DIMENSION(:) :: total_dos, integrate_dos, func

        WRITE(*,*)'running tetrahedron dos'

        OPEN(104,FILE='dos_tet.dat',STATUS='unknown')
        OPEN(105,FILE='dos2_tet.dat',STATUS='unknown')

!        CALL set_omdos(ndyn,wmesh)  ! allocates om(:) and dos(:,:) arrays of size mesh
        CALL my_set_omdos(ndyn,wmesh) !should I use my own?

        WRITE(104,*) '#la,om(i),dos(la,i),total_dos(i),sum(dos_tet(1:i))*dw'
        WRITE(105,*) '#junk,om(i),total_dos(i),integrate_dos(i),tretrados(test)'

        ALLOCATE(integrate_dos(wmesh),total_dos(wmesh),func(nkc))
        ALLOCATE(evc(nkc))

        wmin = 0d0
        total_dos = 0d0

        DO lambda=1,ndyn
!            evc=eigenval(lambda,:)!!!!!!
            evc=SQRT(eivals(lambda,:))
            CALL calc_tet(wmesh,wmin,wmax,nkc,om,kpc,evc,evc) !wmin and wmax seems not used here?
            integrate_dos=0d0
            DO i=1,wmesh
                total_dos(i)=total_dos(i)+dos_tet(i)   ! this is the sum over all bands
                integrate_dos(i)=integrate_dos(i)+sum(total_dos(1:i))*(wmax-wmin)/wmesh  ! this is the cumulative dos
                WRITE(104,3)lambda,om(i),dos_tet(i),total_dos(i),sum(dos_tet(1:i))*(wmax-wmin)/wmesh
            END DO
        END DO

        ALLOCATE(lg(nkc))
        func=1
        DO i=1,wmesh
!            evc=eigenval(1,:)!!!!!!!!!!
            evc=SQRT(eivals(1,:))
            CALL tet_sum(om(i),nkc,evc,func,res,lg)
            WRITE(105,4) om(i),',',total_dos(i),',',integrate_dos(i),',',res
        END DO

        WRITE(*,*) 'tetrahedron dos done'

        DEALLOCATE(dos)
        DEALLOCATE(evc,lg)
        DEALLOCATE(integrate_dos,total_dos,func)
        CLOSE(104)
        CLOSE(105)

3 format(i5,4(3x,g11.5))
4 format(g11.5,3(3x,a,3x,g11.5))
    END SUBROUTINE calc_dos_tet
!-------------------------------------------------------------------------------------------------------
    SUBROUTINE calc_gruneisen
    !!calculate gruneisen
        IMPLICIT NONE
        INTEGER :: i,j
        INTEGER,ALLOCATABLE :: mp(:)
        REAL(8),ALLOCATABLE :: sorted(:)
        REAL(8) :: cv, gama

        WRITE(*,*) 'running thermal calculation'

        CALL read_fc23

        !calculate mode gruneisen on the q path
        CALL make_kp_bs2
        CALL allocate_eig_bs(ndyn,nkp_bs,ndyn)
        CALL get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,ndyn,eigenvec_bs,veloc)
        OPEN(644,file='path_grun.dat',status='unknown')
        CALL gruneisen(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,eigenvec_bs,644,grun_bs)

        CLOSE(644)

        ALLOCATE(mp(ndyn),sorted(ndyn))
        mp = 0; sorted = 0d0

        OPEN(106,FILE='Dispersion_grun.dat',STATUS='unknown')
        DO i=1,SIZE(grun_bs,dim=2)
            CALL sort(ndyn,REAL(grun_bs(:,i)), mp,ndyn) !crucial
            DO j=1,ndyn
               sorted(j) = REAL(grun_bs(mp(j),i))
            END DO
            WRITE(106,*) i, sorted(:)
        END DO
        CLOSE(106)

        CALL deallocate_eig_bs
        DEALLOCATE(dk_bs)
        DEALLOCATE(mp,sorted)

        !calculate mode gruneisen on all k mesh
        CALL allocate_eig(ndyn,nibz,nkc)
        IF(ALLOCATED(dkc)) DEALLOCATE(dkc)
        ALLOCATE(dkc(nkc))
        dkc=0
        CALL get_frequencies(nkc,kpc,dkc,ndyn,eigenval,ndyn,eigenvec,veloc)
!        DEALLOCATE(dkc)

!        ALLOCATE(dkc(nibz))
!        dkc=0
!
!        CALL get_frequencies(nibz,kibz,dkc,ndyn,eivalibz,ndyn,eivecibz,velocibz)
!        DEALLOCATE(dkc)


        OPEN(645,file='all_grun.dat',status='unknown')

        CALL gruneisen(nkc,kpc,dkc,ndyn,eigenval,eigenvec,645,grun)

        CLOSE(645)

        DEALLOCATE(dkc)


        !calculate average gruneisen
!        CALL subst_eivecs(ndyn,nkc,eigenval,eigenvec,kpc,npos,mappos)

        CALL gruneisen_fc

        OPEN(646,file='thermal.dat',status='unknown')

        CALL calculate_thermal(nkc,wk,ndyn,eigenval,grun,tmin,tmax,646,veloc)

        CALL my_thermal(cv,gama)

        WRITE(33,*) 'my calculated gruneisen:',gama
        WRITE(33,*) 'my calculated specific heat:',cv

        CLOSE(646)
        DEALLOCATE(grun)

        WRITE(*,*) 'thermal calculation done'
    END SUBROUTINE calc_gruneisen
!-------------------------------------------------------------------------------------------------------
    SUBROUTINE my_thermal(cv,gama) !everthing is in standard unit
    !!calculate multiple thermal dynamical terms
        IMPLICIT NONE
        INTEGER :: la,k
        REAL(8) :: cv_nk,  x, temp, nbe, nbx
        REAL(8),INTENT(out) :: cv, gama
        REAL(8) :: volume

        CALL calculate_volume(r1,r2,r3,volume) !get volume, in angstrom^3
        temp = temperature*(100*h_plank*c_light)/k_b ! in K

        cv = 0d0; gama = 0d0
        DO k=1,nkc
        DO la=1,ndyn
            x= SQRT(ABS(eivals(la,k)))*cnst*100*c_light*h_plank/k_b/temp

            nbx=nbe(x,1d0,0) ! use quantum one?
            cv_nk = x*x*nbx*(1+nbx)
            cv = cv + cv_nk*wk(k) !weighted?

            gama = grun(la,k)*cv_nk*wk(k) !unit?
        END DO
        END DO
        gama = gama/cv
        cv = cv/atom_number*n_avog*k_b


    END SUBROUTINE my_thermal
!-------------------------------------------------------------------------------------------------------
subroutine get_k_info3(q,nkt,ex) !,i1,j1,k1,gg1,gg2,gg3,inside)
!! legacy code
!! scans the kpoints to identify the index of the input q

 implicit none
 real(8) q(3),ll
 integer inside,i1,j1,k1,nkt,mc(3),i
 logical ex

3 format(a,i6,2x,3(1x,g9.3),2x,f10.4,L3)

 nkt=0; ex=.false.
 loop: do i=1,nkc
!    if( (q(1).myeq.kpc(1,i)) .and. (q(2).myeq.kpc(2,i)) .and. (q(3).myeq.kpc(3,i)) ) then
    ll=length(q-kpc(:,i))
!   if( length(q-kpc(:,i)) .lt. 1d-4) then
!       write(*,3)'i,k,q,k-q=',i,kpc(:,i),q,ll
    if( ll .lt. 1d-4) then
       nkt=i
       ex=.true.
       exit loop
    endif
 enddo loop
! write(*,*)'nkt=',nkt
 if (nkt.eq.0) then
   write(34,3)'GET_K_INFO3: qpoint not found! nq,q,l,exist?',nkt,q,length(q),ex
!  stop
 endif
 end subroutine get_k_info3
!-------------------------------------------------------------------------------------------------------
    subroutine subst_eivecs(ndn,nk,eival,eivec,kp,nposi,map)
!! legacy code, not used
!! this subroutine keeps the eivecs(q) for q in the nposi list and replaces the
!! other eivec(q') (for q'=-q) by the conjugate of eivec(q) in order to get
!! rid of the ambiguity in the case of the degenerate case, and assure e(-q)=conjg(e(q)).

 implicit none
 integer nk,ndn,nposi,map(nk),i1,j1,k1,n,inside,mk3,n3,la,mu
 real(8) kp(3,nk) ,eival(ndn)
 complex(8) eivec(ndn,ndn,nk),aux(ndn)
 logical exists

 write(34,*)'SUBST_EIVECS: nk,nposi=',nk,nposi,'=============='
! do n=1,nk
!    write(ulog,*)'ik,map(ik)=',n,map(n)
! enddo
! do n=nposi+1,nk
 do n=1,nposi  ! these are the kpoints that have a negative in the kpc list
    n3=map(n)      ! actual index of the kpoint
!   call get_k_info(-kp(:,n3)-shft,NC,mk3,i1,j1,k1,g1,g2,g3,inside)
!write(ulog,*)'npos,ik,k ',n,n3,kp(:,n3)
    call get_k_info3(-kp(:,n3),mk3,exists)
  if (exists) then
    do la=1,ndn
       aux=eivec(:,la,n3)
!      write(ulog,*)'for kpoint,branch=',n3,la
       do mu=1,ndn
          if(abs(eivec(mu,la,mk3) - conjg(aux(mu))) .gt. 1d-4) then
          if(abs(eivec(mu,la,mk3) + conjg(aux(mu))) .gt. 1d-4) then
            write(34,4)'WILL SWITCH:mu,k,eivl,ev(-k),ev*(k) =',mu,kp(:,n3),eival(mu),eivec(mu,la,mk3) , conjg(aux(mu))
          endif
          endif
       enddo
       eivec(:,la,mk3) = conjg(aux)
    enddo
  endif
 enddo
 write(34,*)'EXITING SUBST_EIVECS =========================='

4 format(a,i4,99(1x,f9.3))

 end subroutine subst_eivecs
!-------------------------------------------------------------------------------------------------------
    SUBROUTINE gruneisen(nkp,kp,dk,ndn,eivl,eivc,ugr,grn)
!! legacy code, calculate mode gruneisen
!! takes the eigenvalues (w) and eigenvectors calculated along some
!! crystalline directions and calculates the corresponding mode
!! gruneisen parameters
        IMPLICIT NONE
        INTEGER :: ik,i0,nkp,la,al,be,ga,j,k,j0,k0,ta1,ta2,t,ugr,ndn
        REAL(8) :: mi,mj,rx3(3),rx2(3),qq(3),denom,qdotr,omk,mysqrt
        REAL(8) :: kp(3,nkp),dk(nkp),eivl(ndn,nkp)
        COMPLEX(8) :: zz,one,term
        COMPLEX(8) :: grn(ndn,nkp), eivc(ndn,ndn,nkp)

        one = cmplx(1d0,0d0)
!        OPEN(ugr,FILE='gruneisen.dat',STATUS='unknown')

        WRITE(ugr,*)'# la,nk,dk(nk),kp(:,nk),om(la,nk),gruneisen(la,nk))'

        DO ik=1,nkp
        qq(:) = kp(:,ik)

        DO la=1,ndn
        !   write(ulog,6)'Starting gruneisen parameters for kpoint & band# ',ik,la,kp(:,ik)
        !   write(ulog,*)' i,la,t,fcs_3(t),term(2),6*eival(la,i),rr3(3),grun(la,i)'
        grn(la,ik) = 0
        denom = 6 * eivl(la,ik)**2/cnst/cnst
        omk = eivl(la,ik)
        DO i0=1,natoms0
             mi = atom0(i0)%mass
           tloop: DO t=1,nterms(3)

             IF ( i0 .ne. iatomterm_3(1,t) ) cycle tloop
             al = ixyzterm_3(1,t)
             be = ixyzterm_3(2,t)
             ga = ixyzterm_3(3,t)
        !        i0 = iatomcell0(iatomterm_3(1,t))
             j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ; mj = atom0(j0)%mass
             k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k)
             ta1= al + 3*(i0-1)
             ta2= be + 3*(j0-1)

        ! rr2(:)=iatomcell(1,j)*r1+iatomcell(2,j)*r2+iatomcell(3,j)*r3
        !  be careful: it has to be the translations R not atompos!
             rx2(:) = atompos(:,j) - atompos(:,j0)  ! R
             rx3(:) = atompos(:,k)                  ! R+tau

             qdotr =  ( qq .dot. rx2)
             zz = cdexp( ci * qdotr )
        !! term = - ampterm_3(t)*fcs_3(igroup_3(t))*zz*eivc(ta1,la,i)*conjg(eivc(ta2,la,i))*rr3(ga)/sqrt(mi*mj)
             term = - fcs_3(t) * zz * eivc(ta2,la,ik)*conjg(eivc(ta1,la,ik))*rx3(ga)/sqrt(mi*mj)
             grn(la,ik) = grn(la,ik) + term
        !        write(ulog,7)i,la,t,fcs_3(t),rr2,qdotr,zz,rr3,grn(la,ik)
           ENDDO tloop
        ENDDO
        grn(la,ik) = grn(la,ik)/denom
        IF (aimag(grn(la,ik)) .gt. 1d-4) then
!           WRITE(ulog,*) 'GRUNEISEN: ... has large imaginary part!! for nk,la=',ik,la,grn(la,ik)
        !      stop
        ENDIF
        WRITE(ugr,6)' ',la,ik,dk(ik),kp(:,ik),omk,REAL(grn(la,ik))
        ENDDO
        !   write(ugr,8)ik,dk(ik),kp(:,ik),(real(grn(la,ik)),la=1,ndn)
        ENDDO

!    CLOSE(ugr)
5 format(4i7,9(1x,g10.4))
6 format(a,2i5,99(1x,g10.4))
7 format(i5,i5,i6,99(1x,g10.4))
8 format(i8,99(1x,g10.4))
! deallocate(eivl,eivc,kg,grn)

    END SUBROUTINE gruneisen
!-------------------------------------------------------------------------------------------------------
subroutine gruneisen_fc

 implicit none
 integer i0,al,be,ga,j,k,k0,t,t2,igr2,t2old,grold
 real(8) rx3(3),phi  !mi,mj,

 if( .not. allocated(grun_fc)) allocate(grun_fc(nterms(2)))
 write(34,*)'GRUNEISEN_FC: term , phi2_ij , grun'
 do igr2=1,ngroups(2)
 t2old = 0 ; grold = 0 ; grun_fc=0
 do t2=1,nterms(2)
    if( igroup_2(t2) .eq. grold ) cycle
    if( igroup_2(t2) .ne. igr2) cycle
    i0 = iatomcell0(iatomterm_2(1,t2))
    j  = iatomterm_2(2,t2)
    al = ixyzterm_2(1,t2)
    be = ixyzterm_2(2,t2)
!    phi= ampterm_2(t2)*fcs_2(igroup_2(t2))
    phi= fcs_2(t2)
! k1    grun_fc(t2) = 0
!   grn = 0
    tloop: do t=1,nterms(3)
         if ( i0 .ne. iatomterm_3(1,t) ) cycle tloop
         if ( j  .ne. iatomterm_3(2,t) ) cycle tloop
         if ( al .ne. ixyzterm_3(1,t) ) cycle tloop
         if ( be .ne. ixyzterm_3(2,t) ) cycle tloop
         ga = ixyzterm_3(3,t)
         k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k)
         rx3(:) = atompos(:,k)
         grun_fc(t2) = grun_fc(t2) - rx3(ga)  / 6 /phi*fcs_3(t)
!        grn = grn - rr3(ga)  / 6 /phi*fcs_3(t) !ampterm_3(t) *fcs_3(igroup_3(t))
!         write(ulog,6)'phi3: ', t ,phi,-ampterm_3(t)*fcs_3(igroup_3(t)) , rr3(ga),grn
!         write(ulog,6)'phi3: ', t ,phi,-fcs_3(t) , rr3(ga),grn
    enddo tloop
    grold = igr2
    write(34,7)igr2,i0,al,j,be,phi,grun_fc(t2)
 enddo
 enddo

6 format(a,i7,9(1x,g10.4))
7 format('****',i7,2(4x,'(',i3,',',i1,')'),9(1x,g10.4))
 end subroutine gruneisen_fc
!-------------------------------------------------------------------------------------------------------
    SUBROUTINE calculate_thermal(nk,wk,ndyn,eival,grn,tmn,tmx,ual,veloc)
        !!legacy code, calculate multiple thermal dynamical terms
        IMPLICIT NONE
        INTEGER nk,ndyn,b,k,ual,itemp,nat,ntmp,iter,al,be
        REAL(8) wk(nk),eival(ndyn,nk),veloc(3,ndyn,nk),nbx,cv,cv_nk,cv2,nbe,dlogv,dalpha,kapoverl(3,3)
        REAL(8) temp,x,alpha,bulk_modulus,b0,tmn,tmx,gama,a1,c11,c44,etot,free,pres,pres0,magv
        COMPLEX(8) grn(ndyn,nk)

        ! bulk modulus = a_0^2/V d^2E/da^2 evaluated at a_0 (equilibrium lattice parameter)
        ! call mechanical(b0,c11,c44)  ! this is the bulk_mod at T=0

        ! B(T,Veq)=B(T=0,Veq(T=0)) - pres0
WRITE(*,*) 'CHECK MARK2'
        nat = ndyn/3
!        OPEN(ual,FILE='thermal.dat',STATUS='unknown')
        WRITE(ual,'(a132)')'# temperature(K) ,alpha (1/K) , Cv (J/K/mol), gama , E_tot(J/mol) , &
                    &   E_free , Pressure(GPa) , P0 , Bulk_mod(GPa) , kappa/L(nm) '

        ntmp=ntemp !60
        DO itemp=1,ntmp

            temp=tmn+(tmx-tmn)*(itemp-1)**3/(ntmp-1d0+1d-8)**3   ! this temp is in Kelvin
        IF (temp.le.0) THEN
           WRITE(34,*) 'temperature not in the proper range!!', temp
           STOP
        ENDIF
        alpha=1d2; dalpha=1d9; iter=0
        DO WHILE (abs(dalpha) .gt. abs(alpha)/1000 .and. iter .lt. 50)
            iter=iter+1
            dlogv=temp*alpha
            CALL mechanical(b0,c11,c44,dlogv)  ! this is the bulk_mod at T=0
            CALL energies(nk,wk,ndyn,eival,grn,temp,etot,free,pres,pres0,cv2)
            bulk_modulus = b0 - pres0

            a1=0 ; cv=0 ; gama = 0; kapoverl=0
! in the sum over k, all quantities are per unitcell
       DO k=1,nk
       DO b=1,ndyn
!         if(eival(b,k) .lt.0) then
!            x=0
!            write(ulog,*) 'ALPHA: negative eival for band,kp# ',b,k,eival(b,k)
!            write(ulog,*) ' will use its absolute value instead!'
!      else
!         endif
! to convert eV/A^2/mass to Hz we need to take the sqrt, multiply by cnst*100*c_light
            x=(h_plank*eival(b,k)*100*c_light)/k_b/temp
            IF (x.gt.60) THEN
              cv_nk = 0
            ELSE
              nbx=nbe(x,1d0,classical)
              cv_nk = x*x*nbx*(1+nbx) !/4/sinh(x/2)/sinh(x/2)
            ENDIF
            cv = cv + cv_nk*wk(k)   ! this is in J/K units : per unit cell
            a1 = a1 + grn(b,k)*cv_nk*wk(k)
            magv=sqrt(dot_product(veloc(:,b,k),veloc(:,b,k)))
          DO al=1,3
          DO be=1,3
             kapoverl(al,be)=kapoverl(al,be)+cv_nk*wk(k)*veloc(al,b,k)*veloc(be,b,k)/(magv+1d-20)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
! multiplied by MFP(nm), it is the thermal conductivity
       kapoverl = kapoverl *k_b/volume_r*1d30*c_light*1d-9
       gama = a1 / cv
       cv = cv/nat*n_avog*k_b  ! this is per mole    )/(volume_r*1d-30)  !
       IF (abs(cv2-cv).gt.1d-5 ) then
          WRITE(34,4)'CALCULATE_THERMAL:temp, cv from energies ne cv ',iter,temp,cv2,cv
!      stop
       ENDIF
       dalpha = a1/bulk_modulus /n_avog*nat/(volume_r*1d-30) - alpha
       alpha = a1/bulk_modulus /n_avog*nat/(volume_r*1d-30)
       WRITE(*,4)'CALCULATE_THERMAL:i,T,a,da=',iter,temp,alpha,dalpha,magv,kapoverl
    ENDDO
    WRITE(ual,3)temp,alpha,cv,gama,etot,free,pres*1d-9,pres0*1d-9,bulk_modulus*1d-9,kapoverl
    ENDDO

!    CLOSE(ual)
3 format(99(2x,g10.4))
4 format(a,i3,9(2x,g11.5))

END SUBROUTINE calculate_thermal

subroutine mechanical(bulk,c11,c44,dlogv)
!!calculate bulk modulus?
 implicit none
 integer i0,al,be,j,t
 real(8) bulk,rija,rijb,c11,c44,dlogv

! write(ulog,*)'BULK_MOD: i0, al, j, be, rija,rijb,fcs_2(t),bulk'
 bulk=0; c11=0 ; c44=0
 do i0=1,natoms0
 do t=1,nterms(2)
    if ( i0 .ne. iatomterm_2(1,t) ) cycle
    al =  ixyzterm_2(1,t)
    be =  ixyzterm_2(2,t)
    j  = iatomterm_2(2,t)
    rija = atompos(al,j)-atompos(al,i0)
    rijb = atompos(be,j)-atompos(be,i0)
!   write(ulog,3)i0, al, j, be, rija,rijb,fcs_2(t),bulk
    bulk = bulk - fcs_2(t)*(1-2*grun_fc(t)*dlogv)*rija*rijb
    if (al.eq.1 .and. be.eq.1) then
       c11 = c11 - fcs_2(t)*(1-2*grun_fc(t)*dlogv)*rija*rijb
    endif
 enddo
 enddo
 bulk = bulk/volume_r/18
 c11 = c11/volume_r/2

 write(34,*)'BULK_MODULUS: in eV/A^3 ',bulk
 bulk = bulk*1d30*ee
 write(34,*)'BULK_MODULUS: in SI, Mbar ',bulk,bulk*1d-11
 c11 = c11*1d30*ee
 write(34,*)'C11 in SI, Mbar ',c11,c11*1d-11
3 format(4(i4),9(2x,f10.5))

 end subroutine mechanical

 subroutine energies(nk,wk,ndyn,eival,grn,temp,etot,free,pres,pres0,cv)
 !!calculate specific heat?
!! calculate total and free energies within QHA, at a given temperature (temp in Kelvin)
 implicit none
 integer nk,ndyn,b,k,nat
 real(8) wk(nk),eival(ndyn,nk)
 real(8) temp,x,cv_nk,cv,hw,free,etot,pres,nbe,mdedv,pres0,nbx
 complex(8) grn(ndyn,nk)

    nat = ndyn/3
    if (temp.le.0) then
       write(34,*)'temperature not in the proper range!!',temp
       stop
    endif
    etot=0 ; cv=0 ; free=0 ; pres=0
    mdedv= 0.35388/(20.8**3-20.0**3)/ab**3 ! this is -dE/dV in eV/A^3
    mdedv = mdedv*1d+30*ee ! (in J/m^3 )
    do k=1,nk
    do b=1,ndyn
       if(eival(b,k) .lt.0) then
          x=0
          write(34,*) 'ALPHA: negative eival for band,kp# ',b,k,eival(b,k)
          write(34,*) ' will use its absolute value instead!'
!      else
       endif
! to convert eV/A^2/mass to Hz we need to take the sqrt, multiply by cnst*100*c_light
       x=(h_plank*eival(b,k)*100*c_light)/k_b/temp
       if (x.gt.60) then
           cv_nk = 0
       else
           nbx=nbe(x,1d0,classical)
           cv_nk = x*x*nbx*(1+nbx) !/4/sinh(x/2)/sinh(x/2)
       endif
       hw = x*k_b*temp  ! hbar*omega in Joules
       cv  = cv + cv_nk*wk(k)   ! this is in J/K units : per unit cell
       etot= etot + hw * (nbe(x,1d0,classical) + 0.5d0) * wk(k)
       free= free + (0.5d0*hw + k_b*temp*log(1-exp(-x))) * wk(k)
       pres= pres + hw * (nbe(x,1d0,classical) + 0.5d0) * wk(k) * grn(b,k)
    enddo
    enddo
    cv = cv/nat*n_avog*k_b  ! this is per mole    )/(volume_r*1d-30)  !
    etot = etot/nat*n_avog   ! convert from Joule/cell to Joule per mole
    free = free/nat*n_avog   ! convert from Joule/cell to Joule per mole
    pres0= pres/(volume_r*1d-30) ! * 1d-8  ! 1d-8 is to convert to kbar
    pres = pres0+ mdedv
3 format(9(2x,g11.5))

 end subroutine energies
!-------------------------------------------------------------------------------------------------------
 SUBROUTINE GetElastic_final
 !!elastic constant approximation from simple formula, using trial fc2
    IMPLICIT NONE
    INTEGER :: i,j,al,be,ga,de
    INTEGER :: v1, v2, voigt
    REAL(8) :: Rij_ga,Rij_de,cell_volume
    REAL(8),DIMENSION(6,6) :: temp !call inverse matrix will destroy original matrix, so

    CALL calculate_volume(r1,r2,r3,cell_volume)

    elastic = 0d0;compliance = 0d0!don't forget to initialize as 0

    !4 nested xyz loop
    DO al=1,3
    DO be=1,3
    DO ga=1,3
    DO de=1,3

    v1 = voigt(al,be)
    v2 = voigt(ga,de)

    !atom ij sum
    DO i=1,atom_number
    DO j=1,tot_atom_number
        Rij_ga = every_atom(j)%R(ga)+every_atom(j)%tau(ga)-every_atom(i)%R(ga)-every_atom(i)%tau(ga)
        Rij_de = every_atom(j)%R(de)+every_atom(j)%tau(de)-every_atom(i)%R(de)-every_atom(i)%tau(de)
        elastic(v1,v2) = elastic(v1,v2) + trialfc2_value(i,j)%phi(al,be)*Rij_ga*Rij_de/cell_volume
    END DO !atom i loop
    END DO !atom j loop

    END DO !de loop
    END DO !ga loop
    END DO !be loop
    END DO !al loop

    !get inverse elastic
    temp = elastic
    CALL invers_r(temp, compliance,6)

    WRITE(33,*)'elastic'
    DO i=1,6
        WRITE(33,*) elastic(i,:)/cell_volume*100*SQRT(2d0)
    END DO
    WRITE(33,*)

    WRITE(33,*)'compliance'
    DO i=1,6
        WRITE(33,*) compliance(i,:)
    END DO
    WRITE(33,*)
 END SUBROUTINE GetElastic_final
!-------------------------------------------------------------------------------------------------------
 SUBROUTINE GetElastic2
  !!elastic constants calculation based on Wallace book page 82
 !!his S(miu) is my u_tau; his u_ij is my eta_ij; his eta_ij is something else
 !!problem that his u_tau and eta are relevant while mine are independent, his coefficient X should be 0
    IMPLICIT NONE

    INTEGER :: i,j,k,l,m,n
    INTEGER :: vij,vkl,vik,vjl,vjk,vil
    INTEGER :: fyetran(2)
    REAL(8),DIMENSION(6,6) :: C2,A2
    INTEGER :: voigt
    REAL(8) :: short,cell_volume

    CALL calculate_volume(r1,r2,r3,cell_volume)
    C2=0d0;A2=0d0

    OPEN(52,FILE='check_EC2.dat',STATUS='unknown',POSITION='append',ACTION='write')

    DO i=1,3
    DO j=1,3
    DO k=i,3
    DO l=j,3

    vij = voigt(i,j)
    vkl = voigt(k,l)
    vik = voigt(i,k)
    vjl = voigt(j,l)

    DO m=1,atom_number
    DO n=1,tot_atom_number
        !----below is for A2_ikjl----
        A2(vik,vjl) = A2(vik,vjl) - 0.5/cell_volume*&
        &trialfc2_value(m,n)%phi(i,k)*(every_atom(n)%R(j)+every_atom(n)%tau(j)-every_atom(m)%tau(j))*&
        &(every_atom(n)%R(l)+every_atom(n)%tau(l)-every_atom(m)%tau(l))

    END DO !n loop
    END DO !m loop


    END DO !l loop
    END DO !k loop
    END DO !j loop
    END DO !i loop
!    DO vik=1,6
!    DO vjl=1,6
!        fyetran=inv_voigtMap(vik)
!        i=fyetran(1);k=fyetran(2)
!        fyetran=inv_voigtMap(vjl)
!        j=fyetran(1);l=fyetran(2)
!        DO m=1,atom_number
!        DO n=1,tot_atom_number
!            !----below is for A2_ikjl----
!            A2(vik,vjl) = A2(vik,vjl) - 0.5/cell_volume*&
!            &myfc2_value(m,n)%phi(i,k)*(every_atom(n)%R(j)+every_atom(n)%tau(j)-every_atom(m)%tau(j))*&
!            &(every_atom(n)%R(l)+every_atom(n)%tau(l)-every_atom(m)%tau(l))
!
!        END DO !n loop
!        END DO !m loop
!
!    END DO
!    END DO

    !------------convert A into C according to Wallace------------------
    DO i=1,3
    DO k=1,3
    DO j=i,3
    DO l=k,3
        vij = voigt(i,j)
        vkl = voigt(k,l)
        vil = voigt(i,l)
        vjk = voigt(j,k)
        vik = voigt(i,k)
        vjl = voigt(j,l)

        C2(vij,vkl) = A2(vik,vjl)+A2(vjk,vil)-A2(vij,vkl)

    END DO
    END DO
    END DO
    END DO

    elastic = 0d0;compliance = 0d0
    WRITE(52,*)'==== vij,vkl,Elastic Constants C2 ===='
    DO vij=1,6
    DO vkl=1,6
        IF(ABS(C2(vij,vkl)).gt.1d-8) THEN
            WRITE(52,*) vij,vkl,C2(vij,vkl)
        END IF
        elastic(vij,vkl) = C2(vij,vkl)
    END DO
    END DO
    CALL invers_r(elastic,compliance,6)
    elastic = C2 !recover if needed
    CLOSE(52)

    WRITE(33,*)'elastic'
    DO i=1,6
        WRITE(33,*) elastic(i,:)/cell_volume*100*SQRT(2d0)
    END DO
    WRITE(33,*)

    WRITE(33,*)'compliance'
    DO i=1,6
        WRITE(33,*) compliance(i,:)
    END DO
    WRITE(33,*)
 END SUBROUTINE GetElastic2
!-------------------------------------------------------------------------------------------------------
 SUBROUTINE Extract_xy(idx,odx1,odx2) 
  !!utility subroutine, for pressure part
  !!i.e. give idx = 4(xy), output odx1 = 1(x), odx2 = 2(y). etc
    IMPLICIT NONE
    INTEGER, INTENT(in) :: idx
    INTEGER, INTENT(out) :: odx1, odx2

    odx1 = INT((idx+2)/d)
    odx2 = MOD(idx+2,d)+1

 END SUBROUTINE Extract_xy
!-------------------------------------------------------------------------------------------------------
 SUBROUTINE Infer_xyz(idx,odx1,odx2) 
  !!utility subroutine, for pressure part
  !!i.e. give idx = 2(y), output odx1 = 3(z), odx2 = 1(x). etc
    IMPLICIT NONE
    INTEGER, INTENT(in) :: idx
    INTEGER, INTENT(out) :: odx1, odx2

    SELECTCASE(idx)
        CASE(1)
            odx1 = 2; odx2 = 3
        CASE(2)
            odx1 = 3; odx2 = 1
        CASE(3)
            odx1 = 1; odx2 = 2
    ENDSELECT

 END SUBROUTINE Infer_xyz
!-------------------------------------------------------------------------------------------------------
SUBROUTINE add_Pressure
!!add pressure factor to energy and its gradients
    IMPLICIT NONE
    INTEGER :: i, temp
    INTEGER :: xyz1,xyz2
    INTEGER :: odx1,odx2
    REAL(8) :: cell_volume
    TYPE(vector) :: r1_now, r2_now, r3_now
    REAL(8) :: mix,A(3),B(3),C(3),strain_dot_stress

    !----- get current volume with current trans_vec -----
    r1_now%component = trans_vec(1,:) + (strain(1,:).dot.trans_vec(:,:))
    r2_now%component = trans_vec(2,:) + (strain(2,:).dot.trans_vec(:,:))
    r3_now%component = trans_vec(3,:) + (strain(3,:).dot.trans_vec(:,:))
    CALL calculate_volume(r1_now,r2_now,r3_now,cell_volume)

    !----- add V*eta*sigma term to free energy -----
    strain_dot_stress = 0d0
    DO xyz1=1,3
    DO xyz2=1,3
        strain_dot_stress = strain_dot_stress + strain(xyz1,xyz2)*stress(xyz1,xyz2)
    END DO
    END DO
    F_trial = F_trial + cell_volume*strain_dot_stress

    !----- add regarding terms to free energy gradients -----
    DO i=variational_parameters_size(1)+1,variational_parameters_size(1)+variational_parameters_size(2)
        temp=i-variational_parameters_size(1)
        !----- term1: V*sigma -----
        GradientF_trial(i) = GradientF_trial(i) + cell_volume*stress(INT((temp+2)/d),MOD(temp+2,d)+1)

        !----- term2: dV/deta * eta * sigma -----
        CALL Extract_xy(temp,xyz1,xyz2) !eta(xyz1,xyz2)
        CALL Infer_xyz(xyz1,odx1,odx2) ! R2.cross.R3 = R2(odx1)R3(odx2) - R2(odx2)R3(odx1)
        A = r1_now%component !(1+eta)R1
        B = r2_now%component !(1+eta)R2
        C = r3_now%component !(1+eta)R3


        mix = trans_vec(1,xyz2)*(B(odx1)*C(odx2)-B(odx2)*C(odx1)) + &
            & trans_vec(2,xyz2)*(C(odx1)*A(odx2)-C(odx2)*A(odx1)) + &
            & trans_vec(3,xyz2)*(A(odx1)*B(odx2)-A(odx2)*B(odx1))

        GradientF_trial(i) = GradientF_trial(i) + mix*strain_dot_stress

    END DO

END SUBROUTINE add_Pressure
!-------------------------------------------------------------------------------------------------------

End Module VA_math
