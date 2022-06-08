MODULE check
    USE force_update
    USE broy

    CONTAINS
!==========================================================================================================
    !this test a single variable x(i) versus F and f(i), for many steps
    !this test subroutine should be run after <Allocate_Gradients>
    SUBROUTINE small_test(i,step,many)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: i, many ! which x:f, how many steps
        INTEGER :: j
        REAL(8),INTENT(IN) :: step !incremental step value for x(i)
        REAL(8),DIMENSION(:),ALLOCATABLE :: x,f
        CHARACTER name*2

        WRITE(name,'(i2)') i
        OPEN(517,FILE='gradient_vs_x_'//name//'.dat',POSITION='append')
        OPEN(518,FILE='energy_vs_x_'//name//'.dat',POSITION='append')
        DO j=1,many
            WRITE(*,*)'==========Run#:',j,' ========='
            CALL assign_x(i,step)
            CALL GetF_trial_And_GradientF_trial(kvector)
            CALL collect_xf(x,f)
            WRITE(517,*) x(i),f(i)
            WRITE(518,*) x(i),REAL(F_trial)
            DEALLOCATE(x,f)
        END DO
        CLOSE(517)
        CLOSE(518)
        WRITE(*,*) 'Finished Testing'
    END SUBROUTINE small_test
!---------------------------------------------------------------------------------------------------------
    !this test two variable x(i) and x(j) versus F and f(i), for many steps
    !this test subroutine should be run after <Allocate_Gradients>
    SUBROUTINE small_test_ex(i,j,step,many)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: many,i,j ! how many steps,x(i),x(j)
        INTEGER :: idx,jdx
        REAL(8),INTENT(IN) :: step !incremental step value for x(i)
        REAL(8),DIMENSION(:),ALLOCATABLE :: x,f
        REAL(8) :: rollback
        CHARACTER name*8

        WRITE(name,'(i2,a,i2)') i,'_to_',j
        OPEN(517,FILE='gradient_vs_x_'//name//'.dat',POSITION='append')
        OPEN(518,FILE='energy_vs_x_'//name//'.dat',POSITION='append')

        CALL assign_x(i,-many*step)
        CALL assign_x(j,-many*step)

        !outer loop: iterate x(i)
        DO idx=1,many*2
            WRITE(*,*)'------Run# for idx:',idx,' ------'
            CALL assign_x(i,step)
            !inner loop: iterate x(j)
            rollback = 0d0
            DO jdx=1,many*2
                WRITE(*,*)'------Run# for jdx:',jdx,' ------'
                CALL assign_x(j,step)
                rollback = rollback + step

                CALL GetF_trial_And_GradientF_trial(kvector) !not accurate but more realistic
                !-------------this is too slow---------------
!                IF(i.le.15 .AND. j.le.15) THEN
!                    CALL GetF0_and_V0 !F0 and <V0>
!                    CALL GetV_avg_And_GradientV_avg(kvector)!<V>
!                ELSE
!                    CALL GetV_avg_And_GradientV_avg(kvector)!<dV/d<YY>> and <V>
!                    CALL UpdateTrialFC2 !update effective fc2
!                    CALL GetEigen(kvector) !update eivals
!                    CALL GetF0_and_V0 !F0 and <V0>
!                END IF
!                F_trial=F0+V_avg-V0
                !--------------------------------------------
                CALL collect_xf(x,f)
                WRITE(517,*) x(i),x(j),f(i),f(j)
                WRITE(518,*) x(i),x(j),REAL(F_trial)
                DEALLOCATE(x,f)
            END DO !inner loop
            CALL assign_x(j,-rollback) !rollback to starting x(j)
        END DO !outer loop

        CLOSE(517)
        CLOSE(518)

        WRITE(*,*) 'Finished Testing'
    END SUBROUTINE small_test_ex
!---------------------------------------------------------------------------------------------------------
    !after assign i, j make the rest of strain and atomic deviation to be rhom symmetry
    SUBROUTINE rhom_symmetrize
        IMPLICIT NONE
        INTEGER :: i,j

        DO i=1,3
            atomic_deviation(i,2) = atomic_deviation(1,2) !utau_2^x
        DO j=1,3
            IF(i.eq.j) THEN
                strain(i,j) = strain(1,1) !eta^xx
            ELSE
                strain(i,j) = strain(1,2) !eta^xy
            END IF
        END DO
        END DO

    END SUBROUTINE rhom_symmetrize
!---------------------------------------------------------------------------------------------------------
    !just for rhom test, so i,j values are limited to {4,7,8} which stands for utau_2^x,eta^xx, eta^xy
    SUBROUTINE rhom_contour(i,j,step1,step2,many)
        IMPLICIT NONE

        INTEGER,INTENT(IN):: many,i,j ! how many steps,x(i),x(j)
        INTEGER :: idx,jdx
        REAL(8),INTENT(IN) :: step1,step2 !incremental step value for x(i),x(j)
        REAL(8),DIMENSION(:),ALLOCATABLE :: x,f
        REAL(8) :: rollback
        CHARACTER name*8

        WRITE(name,'(i2,a,i2)') i,'_to_',j
        OPEN(517,FILE='gradient_vs_x_'//name//'.dat',POSITION='append')
        OPEN(518,FILE='energy_vs_x_'//name//'.dat',POSITION='append')

        CALL assign_x(i,-many*step1) !initialize x_i to starting point: -many*step1
        CALL assign_x(j,-many*step2) !initialize x_j to starting point: -many*step2
        CALL rhom_symmetrize

        !outer loop: iterate x(i)
        DO idx=1,many*2
            WRITE(*,*)'------Run# for idx:',idx,' ------'
            CALL assign_x(i,step1)
            !inner loop: iterate x(j)
            rollback = 0d0
            DO jdx=1,many*2
                WRITE(*,*)'------Run# for jdx:',jdx,' ------'
                CALL assign_x(j,step2)
                rollback = rollback + step2

                CALL rhom_symmetrize
!                CALL updateK !for inherited run, turn this off
                CALL GetF_trial_And_GradientF_trial(kvector) !it automatically diagonalize new K, calculate F

                CALL collect_xf(x,f)
                WRITE(517,*) x(i),x(j),f(i),f(j)
                WRITE(518,*) x(i),x(j),REAL(F_trial)
                DEALLOCATE(x,f)
            END DO !inner loop
            CALL assign_x(j,-rollback) !rollback to starting x(j)
        END DO !outer loop

        CLOSE(517)
        CLOSE(518)

        WRITE(*,*) 'Finished Testing'

    END SUBROUTINE rhom_contour
!---------------------------------------------------------------------------------------------------------
    !this test all variable x(:), f(:) from formula vs. finite difference
    !this test subroutine should be run after <initiate_yy>
    SUBROUTINE small_test2(step)
        IMPLICIT NONE
        INTEGER :: i,num_var
        INTEGER :: tau1,atom2,xyz1,xyz2,temp
        REAL(8),INTENT(in) :: step
        REAL(8) :: F1,F2,f_math,f_fd,discrepancy
        REAL(8),DIMENSION(:),ALLOCATABLE :: x_init,f_init
        REAL(8),DIMENSION(:),ALLOCATABLE :: x,f

        OPEN(519,FILE='gradients_test_all.dat')
        WRITE(519,*)'THIS TEST ALL THE GRADIENTS: step=',step
        WRITE(519,*)'==== #, x, gradient from formula, gradient from finite difference ===='
        num_var = SUM(variational_parameters_size)
        ALLOCATE(x(num_var),f(num_var))

        CALL combine_variants(x_init) !record initial x generated by <test_update>
        CALL decompose_variants(x_init)
        CALL GetF_trial_And_GradientF_trial(kvector)
        f_init = GradientF_trial !record initial f generated by first <GetF_....>
        F1 = F_trial !record initial F generated by first <GetF_....>

        DO i=1,num_var
            x = x_init !start from initial x

            !notice for strain, eta_xy and eta_yx should be simultaneously increased
            SELECTCASE(i)
                CASE(8) !xy
                    x(8) = x(8) + step
                    x(10) = x(10) + step
                CASE(10) !yx
                    x(8) = x(8) + step
                    x(10) = x(10) + step
                CASE(9) !xz
                    x(9) = x(9) + step
                    x(13) = x(13) + step
                CASE(13) !zx
                    x(9) = x(9) + step
                    x(13) = x(13) + step
                CASE(12) !yz
                    x(12) = x(12) + step
                    x(14) = x(14) + step
                CASE(14) !zy
                    x(12) = x(12) + step
                    x(14) = x(14) + step
                CASE DEFAULT
                    x(i) = x(i) + step
            ENDSELECT

            CALL decompose_variants(x)
            CALL GetF_trial_And_GradientF_trial(kvector)
            !get new f_math (in the middle),notice <yy> part
            IF(i.le.15) THEN
                f_math = GradientF_trial(i)
            ELSE
                temp = i-15
                tau1=eff_fc2_index(temp)%iatom_number
                atom2=eff_fc2_index(temp)%jatom_number
                xyz1=eff_fc2_index(temp)%iatom_xyz
                xyz2=eff_fc2_index(temp)%jatom_xyz
                f_math = GradientV_cor(tau1,atom2)%phi(xyz1,xyz2)
            ENDIF

            SELECTCASE(i)
                CASE(8) !xy
                    x(8) = x(8) + step
                    x(10) = x(10) + step
                CASE(10) !yx
                    x(8) = x(8) + step
                    x(10) = x(10) + step
                CASE(9) !xz
                    x(9) = x(9) + step
                    x(13) = x(13) + step
                CASE(13) !zx
                    x(9) = x(9) + step
                    x(13) = x(13) + step
                CASE(12) !yz
                    x(12) = x(12) + step
                    x(14) = x(14) + step
                CASE(14) !zy
                    x(12) = x(12) + step
                    x(14) = x(14) + step
                CASE DEFAULT
                    x(i) = x(i) + step
            ENDSELECT
            CALL decompose_variants(x)
            CALL GetF_trial_And_GradientF_trial(kvector) !get F2 (add 2 step to x(i))
            F2 = F_trial
            f_fd = (F2-F1)/(2*step)

            discrepancy = ABS(f_math-f_fd)
            IF(discrepancy.gt.1e-2) THEN
                WRITE(519,9) i,x_init(i),f_math,f_fd,'WRONG'
            ELSE
                WRITE(519,9)i,x_init(i),f_math,f_fd,'CORRECT'
            END IF
        END DO

        DEALLOCATE(x,f)
        CLOSE(519)
        9 format((I3,3X),3(G12.6,3X),A)
    END SUBROUTINE small_test2

    SUBROUTINE small_test3(i,step,n)
    !!this test for x(i),change step(epsilon) for n times,
    !!how the discrepancy between f_math and f_fd look
    !!however this subroutine cannot check for yy gradients
    !!used for logrithm plot
        IMPLICIT NONE
        INTEGER,INTENT(in) :: i,n
        INTEGER :: j,num_var,m
        INTEGER :: tau1,atom2,xyz1,xyz2,temp
        REAL(8),INTENT(INOUT) :: step
        REAL(8) :: F1,F2,f_math,f_fd,discrepancy
        REAL(8) :: V1,V2,v_math,v_fd
        REAL(8),DIMENSION(:),ALLOCATABLE :: x_init,f_init
        REAL(8),DIMENSION(:),ALLOCATABLE :: x,f
        CHARACTER name*3

        WRITE(name,'(i3)') i


        OPEN(520,FILE='gradients_test_x'//name//'_log.dat')
        WRITE(520,*) '==== epsilon, discrepancy, gradient from formula, finite difference ===='
        num_var = SUM(variational_parameters_size)
        ALLOCATE(x(num_var),f(num_var))

        CALL combine_variants(x_init) !record initial x generated by <test_update>
        CALL decompose_variants(x_init)
!        CALL GetF_trial_And_GradientF_trial(kvector)


!        f_init = GradientF_trial !record initial f generated by first <GetF_....>
!
!        f_math = GradientF_trial(i)
        CALL initiate_yy(kvector)
        CALL GetV_avg_And_GradientV_avg2(kvector)
        IF(i.le.variational_parameters_size(1)) THEN
            IF(MOD(i,d).eq.0) THEN
                v_math=GradientV_utau(d,INT(i/d))
            ELSE
                v_math=GradientV_utau(MOD(i,d),INT(i/d+1))
            END IF
        ELSE
            temp = i-variational_parameters_size(1)
            IF(MOD(temp,d).eq.0) THEN
                v_math=GradientV_eta(INT(temp/d),d)
            ELSE
                v_math=GradientV_eta(INT(temp/d+1),MOD(temp,d))
            END IF
        END IF

        DO j=1,n
            x = x_init !start from initial x

            !notice for strain, eta_xy and eta_yx should be simultaneously increased
            SELECTCASE(i)
                CASE(8) !xy
                    x(8) = x(8) - step
                    x(10) = x(10) - step
                CASE(10) !yx
                    x(8) = x(8) - step
                    x(10) = x(10) - step
                CASE(9) !xz
                    x(9) = x(9) - step
                    x(13) = x(13) - step
                CASE(13) !zx
                    x(9) = x(9) - step
                    x(13) = x(13) - step
                CASE(12) !yz
                    x(12) = x(12) - step
                    x(14) = x(14) - step
                CASE(14) !zy
                    x(12) = x(12) - step
                    x(14) = x(14) - step
                CASE DEFAULT
                    x(i) = x(i) - step
            ENDSELECT

            CALL decompose_variants(x)
!            CALL GetF_trial_And_GradientF_trial(kvector)
!            F1 = F_trial

            CALL GetV_avg_And_GradientV_avg2(kvector)
            V1 = REAL(V_avg)

            !notice for strain, eta_xy and eta_yx should be simultaneously increased
            SELECTCASE(i)
                CASE(8) !xy
                    x(8) = x(8) + 2*step
                    x(10) = x(10) + 2*step
                CASE(10) !yx
                    x(8) = x(8) + 2*step
                    x(10) = x(10) + 2*step
                CASE(9) !xz
                    x(9) = x(9) + 2*step
                    x(13) = x(13) + 2*step
                CASE(13) !zx
                    x(9) = x(9) + 2*step
                    x(13) = x(13) + 2*step
                CASE(12) !yz
                    x(12) = x(12) + 2*step
                    x(14) = x(14) + 2*step
                CASE(14) !zy
                    x(12) = x(12) + 2*step
                    x(14) = x(14) + 2*step
                CASE DEFAULT
                    x(i) = x(i) + 2*step
            ENDSELECT

            CALL decompose_variants(x)
!            CALL GetF_trial_And_GradientF_trial(kvector) !get F2 (add 2 step to x(i))
!            F2 = F_trial
!            f_fd = (F2-F1)/(2*step)

!            discrepancy = ABS(f_math-f_fd)
            CALL GetV_avg_And_GradientV_avg2(kvector)
            V2 = REAL(V_avg)
            v_fd = (V2-V1)/(2*step)
            discrepancy = ABS(v_math-v_fd)

            WRITE(520,*) step,',',discrepancy,v_math,v_fd
!            WRITE(520,*) LOG10(step),',',LOG10(discrepancy)

            step = step / 10
        END DO

        DEALLOCATE(x,f)
        CLOSE(520)
    END SUBROUTINE small_test3
!---------------------------------------------------------------------------------------------------------
    SUBROUTINE small_test3_yy(i,step,n)
    !!this test for a selected yy(K)
    !!change step(epsilon) for n times
    !!how the discrepancy between v_math and v_fd look
    !!used for logrithm plot
        IMPLICIT NONE
        INTEGER,INTENT(in) :: i,n
        INTEGER :: j,num_var
        INTEGER :: tau1,atom2,xyz1,xyz2
        REAL(8),INTENT(INOUT) :: step
        REAL(8) :: V1,V2,v_math,v_fd,discrepancy,yy_init

        CHARACTER name*3

        WRITE(name,'(i3)') i+15

        OPEN(504,FILE='gradients_test_x'//name//'_log.dat')
        WRITE(504,*) '==== epsilon, discrepancy, gradient from formula, finite difference ===='
        num_var = SUM(variational_parameters_size)

        tau1=eff_fc2_index(i)%iatom_number
        atom2=eff_fc2_index(i)%jatom_number
        xyz1=eff_fc2_index(i)%iatom_xyz
        xyz2=eff_fc2_index(i)%jatom_xyz

        CALL initiate_yy(kvector)
        CALL GetV_avg_And_GradientV_avg2(kvector)

        yy_init = yy_value(tau1,atom2)%phi(xyz1,xyz2)
        v_math = GradientV_cor(tau1,atom2)%phi(xyz1,xyz2)
        DO j=1,n
            yy_value(tau1,atom2)%phi(xyz1,xyz2) = yy_init - step
            CALL GetV_avg_And_GradientV_avg2(kvector)
            V1= REAL(V_avg)
            yy_value(tau1,atom2)%phi(xyz1,xyz2) = yy_init + step
            CALL GetV_avg_And_GradientV_avg2(kvector)
            V2 = REAL(V_avg)
            v_fd = (V2-V1)/(2*step)
            discrepancy = ABS(v_math-v_fd)

            WRITE(504,*) step,',',discrepancy,v_math,v_fd
            step = step/10
        END DO

        CLOSE(504)
    END SUBROUTINE small_test3_yy
!---------------------------------------------------------------------------------------------------------
    !add 'step' value to variational_parameter(i)
    SUBROUTINE assign_x(i,step)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: i !index
        REAL(8),INTENT(IN) :: step !value wants to be added to original x(i)

        INTEGER :: mx,temp
        INTEGER :: tau1, atom2, direction1, direction2
        mx = SUM(variational_parameters_size)
        IF(i.le.variational_parameters_size(1)) THEN
            IF(MOD(i,d).eq.0) THEN
                atomic_deviation(d,INT(i/d)) = atomic_deviation(d,INT(i/d)) + step
            ELSE
                atomic_deviation(MOD(i,d),INT(i/d+1)) = atomic_deviation(MOD(i,d),INT(i/d+1)) + step
            END IF
        END IF

        IF(i.gt.variational_parameters_size(1) .AND. i.le.mx-variational_parameters_size(3)) THEN
            temp=i-variational_parameters_size(1)
            IF(MOD(temp,d).eq.0) THEN
                strain(INT(temp/d),d) = strain(INT(temp/d),d) + step
            ELSE
                strain(INT(temp/d+1),MOD(temp,d)) = strain(INT(temp/d+1),MOD(temp,d)) + step
            END IF
        END IF

        IF(i.gt.mx-variational_parameters_size(3)) THEN
            temp=i-variational_parameters_size(1)-variational_parameters_size(2)

            tau1=eff_fc2_index(temp)%iatom_number
            atom2=eff_fc2_index(temp)%jatom_number
            direction1=eff_fc2_index(temp)%iatom_xyz
            direction2=eff_fc2_index(temp)%jatom_xyz

            yy_value(tau1,atom2)%phi(direction1,direction2) =&
            & yy_value(tau1,atom2)%phi(direction1,direction2) + step
        END IF
    END SUBROUTINE assign_x
!---------------------------------------------------------------------------------------------------------
    SUBROUTINE collect_xf(x,f)
        IMPLICIT NONE
        REAL(8),DIMENSION(:),ALLOCATABLE,INTENT(out) :: x,f
        INTEGER :: mx,i,temp
        INTEGER :: tau1, atom2, direction1, direction2

        mx=SUM(variational_parameters_size)
        ALLOCATE(x(mx),f(mx))

        !collect f(:)
        DO i=1,mx
            f(i) = REAL(GradientF_trial(i))
        END DO

        !collect x(:)
        DO i=1,variational_parameters_size(1) !atom deviation
            IF(MOD(i,d).eq.0) THEN
                x(i)=atomic_deviation(d,INT(i/d))
            ELSE
                x(i)=atomic_deviation(MOD(i,d),INT(i/d+1))
            END IF
        END DO

        DO i=variational_parameters_size(1)+1,mx-variational_parameters_size(3) !strain
            temp=i-variational_parameters_size(1)
            IF(MOD(temp,d).eq.0) THEN
                x(i)=strain(INT(temp/d),d)
            ELSE
                x(i)=strain(INT(temp/d+1),MOD(temp,d))
            END IF
        END DO

        DO i=mx-variational_parameters_size(3)+1,mx
            temp=i-variational_parameters_size(1)-variational_parameters_size(2)

            tau1=eff_fc2_index(temp)%iatom_number
            atom2=eff_fc2_index(temp)%jatom_number
            direction1=eff_fc2_index(temp)%iatom_xyz
            direction2=eff_fc2_index(temp)%jatom_xyz

            x(i)=yy_value(tau1,atom2)%phi(direction1,direction2)
       END DO

       !check
!       OPEN(824,FILE='collected_x.dat')
!       DO i=1,mx
!            WRITE(824,*) x(i),f(i)
!       END DO
!       CLOSE(824)

    END SUBROUTINE collect_xf
!---------------------------------------------------------------------------------------------------------
    !used after call <GetV_avg_And_GradientV_avg>
    SUBROUTINE UpdateTrialFC2
        IMPLICIT NONE
        INTEGER :: tau1, atom2, direction1, direction2

        !****turn off for actual calculation****
!        CALL set_yy_0 !for QHA test
!        CALL GetV_avg_And_GradientV_avg(kvector) !update GradientV_cor if all <yy> are made 0
        !***************************************

        DO tau1=1,atom_number
        DO atom2=1,tot_atom_number
        DO direction1=1,d
        DO direction2=1,d

            trialfc2_value(tau1,atom2)%phi(direction1,direction2)=2*GradientV_cor(tau1,atom2)%phi(direction1,direction2)

        END DO
        END DO
        END DO
        END DO

    END SUBROUTINE UpdateTrialFC2

    !set all <yy> to 0 for QHA test
    SUBROUTINE set_yy_0
        IMPLICIT NONE
        INTEGER :: tau1,atom2,direction1,direction2

        DO tau1=1,atom_number
        DO atom2=1,tot_atom_number
        DO direction1=1,d
        DO direction2=1,d

            yy_value(tau1,atom2)%phi(direction1,direction2) = 0d0

        END DO
        END DO
        END DO
        END DO

    END SUBROUTINE set_yy_0
!---------------------------------------------------------------------------------------------------------
!========================================================================================================
    ! update every atom position using info directly from variational params
    SUBROUTINE atompos_Update2
        IMPLICIT NONE
        INTEGER :: i

        !Firstly, store all the previous 'every_atom(:)' for the force update mapping part
        IF(ALLOCATED(prev_atom)) DEALLOCATE(prev_atom)
        ALLOCATE(prev_atom(tot_atom_number))
        DO i=1, tot_atom_number
            prev_atom(i) = every_atom(i)
        END DO

        !tot_atom_number will stay the same

        !Finally, store all information in newly allocated 'every_atom(:)'
        DO i=1, tot_atom_number
            every_atom(i)%label_number = i !nothing really changes for this one
            !below the trans_vec has updated before
            every_atom(i)%R = every_atom(i)%n1*trans_vec(:,1) + &
                &every_atom(i)%n2*trans_vec(:,2) + every_atom(i)%n3*trans_vec(:,3)
            !below the tau position has updated before
            every_atom(i)%tau = iatom(every_atom(i)%type_tau)%pos_tau
            !update x,y,z
            every_atom(i)%x = every_atom(i)%R(1) + every_atom(i)%tau(1)
            every_atom(i)%y = every_atom(i)%R(2) + every_atom(i)%tau(2)
            every_atom(i)%z = every_atom(i)%R(3) + every_atom(i)%tau(3)
            !n1,n2,n3,type_tau will stay the same
            atompos(1,i) = every_atom(i)%x
            atompos(2,i) = every_atom(i)%y
            atompos(3,i) = every_atom(i)%z
            iatomcell(1,i) = every_atom(i)%n1
            iatomcell(2,i) = every_atom(i)%n2
            iatomcell(3,i) = every_atom(i)%n3
            iatomcell0(i) = every_atom(i)%type_tau

            !every_atom(i)%type_R will be updated in subroutine <cellvec_Update>
        END DO
        !for checking
        OPEN(222, file="new2_lat_fc.dat",status='unknown',action='write')
        DO i=1, tot_atom_number
            WRITE(222,*) i,every_atom(i)%x, every_atom(i)%y, every_atom(i)%z
        END DO
        CLOSE(222)

    END SUBROUTINE atompos_Update2
 !---------------------------------------------------------------------------------------------------
 SUBROUTINE uni_fc_update2(rnk, atoms, xyzs)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: rnk, atoms(:), xyzs(:)
        SELECTCASE(rnk)
            CASE(2)
                CALL renew_fc2_2(atoms, xyzs)
            CASE(3)
                CALL renew_fc3_2(atoms, xyzs)
            CASE(4)
                CALL renew_fc4_2(atoms, xyzs)
        ENDSELECT

    END SUBROUTINE uni_fc_update2
!---------------------------------------------------
    SUBROUTINE renew_fc2_2(atoms, xyzs)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: atoms(:), xyzs(:)
        INTEGER :: i, j,k,l,idx
        INTEGER :: atom1,atom2,atom3,atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER :: new_atom1, new_atom2

        REAL(8) :: temp
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: deltau
        REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE :: Y_square
        INTEGER :: found,tau3
        INTEGER :: k_number
        INTEGER :: flag(2)
        REAL(8) :: weight(2)

        ![0] get all the <YY> term again
        k_number = SIZE(kvector)
        ALLOCATE(Y_square(d,tot_atom_number,d,tot_atom_number))
        Y_square=0d0
        !--------calculate independent Y_square for each k,lambda---------
        WRITE(47,*) "===== All independent Y_square ====="
        DO i=1,ifc2_terms
            atom1=indiefc2_index(i)%iatom_number
            direction1=indiefc2_index(i)%iatom_xyz
            atom2=indiefc2_index(i)%jatom_number
            direction2=indiefc2_index(i)%jatom_xyz
            Y_square(direction1,atom1,direction2,atom2)=Get_Y_square(k_number,direction1,atom1,direction2,atom2)/k_number !N=SIZE(kvector)
            WRITE(47,*) get_letter(direction1),atom1,get_letter(direction2),atom2,Y_square(direction1,atom1,direction2,atom2)

        END DO !i loop
        !--------------------generate the rest Y_square by grouping----------
        WRITE(47,*) "===== The rest of Y_square ====="
        DO i=1,d !iatom_xyz
        DO j=1,atom_number !iatom_number
        DO k=1,d !jatom_xyz
        DO l=1,tot_atom_number !jatom_number
            flag=myfc2_group(i,j,k,l)%group !get the group this fc2 belongs to
            weight=myfc2_group(i,j,k,l)%mat !get the weight factor
            !get corresponding indie fc2 labels
            IF(flag(1).ne.0 .AND. flag(2).eq.0) THEN
                !WRITE(47,*) 'Flag=: ', flag
                !WRITE(47,*)'direction1: ',i,'atom1: ',j,'direction2: ',k,'atom2: ',l
                atom1=indiefc2_index(flag(1))%iatom_number
                direction1=indiefc2_index(flag(1))%iatom_xyz
                atom2=indiefc2_index(flag(1))%jatom_number
                direction2=indiefc2_index(flag(1))%jatom_xyz
                Y_square(i,j,k,l)=Y_square(direction1,atom1,direction2,atom2)*weight(1)
            ELSEIF(flag(1).ne.0 .AND. flag(2).ne.0) THEN
                atom1=indiefc2_index(flag(1))%iatom_number
                direction1=indiefc2_index(flag(1))%iatom_xyz
                atom2=indiefc2_index(flag(1))%jatom_number
                direction2=indiefc2_index(flag(1))%jatom_xyz
                Y_square(i,j,k,l)=Y_square(direction1,atom1,direction2,atom2)*weight(1)

                atom1=indiefc2_index(flag(2))%iatom_number
                direction1=indiefc2_index(flag(2))%iatom_xyz
                atom2=indiefc2_index(flag(2))%jatom_number
                direction2=indiefc2_index(flag(2))%jatom_xyz
                Y_square(i,j,k,l)=Y_square(i,j,k,l)+&
                &Y_square(direction1,atom1,direction2,atom2)*weight(2)
            ELSE
                CYCLE
            END IF
        END DO
        END DO
        END DO
        END DO
        WRITE(47,*) '========================================================'
        !up to here the <YY> ranges (0tau1,R2tau2) have all been calculated

!WRITE(36,*)"atom1,xyz1,atom2,xyz2",atoms(1),get_letter(xyzs(1)),atoms(2),get_letter(xyzs(2))

        ![1]get all the static deviation s
        ALLOCATE(deltau(d,SIZE(every_atom)))
        DO atom1=1,SIZE(every_atom)
            deltau(:,atom1)=strain(:,:).dot.every_atom(atom1)%R+every_atom(atom1)%tau+&
            &atomic_deviation(:,every_atom(atom1)%type_tau)
        END DO
!WRITE(34,*)"check mark 1"
        ![2]fc4 term sum in the formula: 0.5*u*Chi*u
        DO idx = 1, SIZE(myfc4_index)
            atom1=myfc4_index(idx)%iatom_number
            new_atom1=atom1
            IF(new_atom1.ne.atoms(1)) CYCLE

            atom2=myfc4_index(idx)%jatom_number
            new_atom2=atom2
            IF(new_atom2.ne.atoms(2)) CYCLE

            direction1 = myfc4_index(idx)%iatom_xyz
            IF(direction1.ne.xyzs(1)) CYCLE !Low efficient
            direction2 = myfc4_index(idx)%jatom_xyz
            IF(direction2.ne.xyzs(2)) CYCLE !Low efficient

            atom3=myfc4_index(idx)%katom_number
            tau3 = every_atom(atom3)%type_tau !for <YY> matching
            atom4=myfc4_index(idx)%latom_number
            direction3 = myfc4_index(idx)%katom_xyz
            direction4 = myfc4_index(idx)%latom_xyz

            temp = 0d0
            !if corresponding fc4 term exists, add into
            IF(ANY(fc4_unique_idx==atom2).AND.ANY(fc4_unique_idx==atom3)&
                &.AND.ANY(fc4_unique_idx==atom4)) THEN
                !get corresponding idx pos
                atom2 = find_loc(fc4_unique_idx,atom2)
                atom3 = find_loc(fc4_unique_idx,atom3)
                atom4 = find_loc(fc4_unique_idx,atom4)

                temp = 0.5*myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)&
                &*deltau(direction3,atom3)*deltau(direction4,atom4)

                !remember to add <YY> term also
                found = ifExistYY(atom3,atom4)
                IF(found.ne.0) THEN
                    temp = temp&
                    &+0.5*myfc4_value(atom1,atom2,atom3,atom4)%chi(direction1,direction2,direction3,direction4)&
                    &*Y_square(direction3,tau3,direction4,found)
                END IF
                !remember to recover the atom labels
                atom2=myfc4_index(idx)%jatom_number
                atom3=myfc4_index(idx)%katom_number
                atom4=myfc4_index(idx)%latom_number
!WRITE(34,*)"atom2,atom3,atom4",atom2,atom3,atom4
!WRITE(34,*)"check term:",temp
            END IF

            myfc2_value(atoms(1),atoms(2))%phi(xyzs(1),xyzs(2)) = &
            &myfc2_value(atoms(1),atoms(2))%phi(xyzs(1),xyzs(2)) + temp
        END DO
!WRITE(34,*)"check mark 2"
        ![3]fc3 term sum in the formula: Psi*u
        DO idx = 1, SIZE(myfc3_index)
            atom1 = myfc3_index(idx)%iatom_number
            new_atom1 = atom1
            IF(new_atom1.ne.atoms(1)) CYCLE !Low efficiency

            atom2 = myfc3_index(idx)%jatom_number
            new_atom2 = atom2
            IF(new_atom2.ne.atoms(2)) CYCLE !Low efficiency

            direction1 = myfc3_index(idx)%iatom_xyz
            IF(direction1.ne.xyzs(1)) CYCLE !Low efficiency
            direction2 = myfc3_index(idx)%jatom_xyz
            IF(direction2.ne.xyzs(2)) CYCLE !Low efficiency

            atom3 = myfc3_index(idx)%katom_number
            direction3 = myfc3_index(idx)%katom_xyz

            temp = 0d0
            !if corresponding fc3 term exists, add into
            IF(ANY(fc3_unique_idx==atom2).AND.ANY(fc3_unique_idx==atom3)) THEN
                !get corresponding idx pos
                atom2 = find_loc(fc3_unique_idx,atom2)
                atom3 = find_loc(fc3_unique_idx,atom3)

                temp = myfc3_value(atom1,atom2,atom3)%psi(direction1,direction2,direction3)*deltau(direction3,atom3)
                !remember to recover the atom labels
                atom2 = myfc3_index(idx)%jatom_number
                atom3 = myfc3_index(idx)%katom_number
            END IF

            myfc2_value(atoms(1),atoms(2))%phi(xyzs(1),xyzs(2)) = &
            &myfc2_value(atoms(1),atoms(2))%phi(xyzs(1),xyzs(2)) + temp
        END DO

        DEALLOCATE(deltau)
!WRITE(34,*)'Successfully updated this fc2!'
    END SUBROUTINE renew_fc2_2
!---------------------------------------------------
    SUBROUTINE renew_fc3_2(atoms, xyzs)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: atoms(:), xyzs(:)
        INTEGER :: i, j, idx
        INTEGER :: atom1,atom2,atom3,atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER :: new_atom1, new_atom2, new_atom3

        REAL(8) :: temp
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: deltau

        ![1]get all the deviation u first
        ALLOCATE(deltau(d,SIZE(every_atom)))
        DO atom1=1,SIZE(every_atom)
            deltau(:,atom1)=strain(:,:).dot.every_atom(atom1)%R+every_atom(atom1)%tau+&
            &atomic_deviation(:,every_atom(atom1)%type_tau)
        END DO
        ![2]update myfc3_value using formula
        !FORMULA: Psi_new = Psi_old + SUM(Chi_old*deltaU)
        DO i=1, SIZE(myfc4_index)
            atom1=myfc4_index(i)%iatom_number
            new_atom1=atom1
            IF(new_atom1.ne.atoms(1)) CYCLE !Low efficiency

            atom2=myfc4_index(i)%jatom_number
            new_atom2=atom2
            IF(new_atom2.ne.atoms(2)) CYCLE !Low efficiency

            atom3=myfc4_index(i)%katom_number
            new_atom3=atom3
            IF(new_atom3.ne.atoms(3)) CYCLE !Low efficiency

            direction1 = myfc4_index(i)%iatom_xyz
            IF(direction1.ne.xyzs(1)) CYCLE !Low efficiency
            direction2 = myfc4_index(i)%jatom_xyz
            IF(direction2.ne.xyzs(2)) CYCLE !Low efficiency
            direction3 = myfc4_index(i)%katom_xyz
            IF(direction3.ne.xyzs(3)) CYCLE !Low efficiency

            atom4=myfc4_index(i)%latom_number
            direction4 = myfc4_index(i)%latom_xyz

            IF(ANY(fc3_unique_idx==new_atom2).AND.ANY(fc3_unique_idx==new_atom3)) THEN
                new_atom2 = find_loc(fc3_unique_idx,new_atom2)
                new_atom3 = find_loc(fc3_unique_idx,new_atom3)

                !need to make sure atomic index for fc4 are also valid
                IF(ANY(fc4_unique_idx==atom2).AND.ANY(fc4_unique_idx==atom3)&
                &.AND.ANY(fc4_unique_idx==atom4)) THEN
                    atom2 = find_loc(fc4_unique_idx,atom2)
                    atom3 = find_loc(fc4_unique_idx,atom3)
                    atom4 = find_loc(fc4_unique_idx,atom4)

                    myfc3_value(new_atom1,new_atom2,new_atom3)%psi(xyzs(1),xyzs(2),xyzs(3)) = &
                    & myfc3_value(new_atom1,new_atom2,new_atom3)%psi(xyzs(1),xyzs(2),xyzs(3)) + &
                    & myfc4_value(atom1, atom2, atom3, atom4)%chi(direction1, direction2, direction3, direction4) * &
                    & deltau(direction4,atom4)
                END IF
            END IF
        END DO

        DEALLOCATE(deltau)

    END SUBROUTINE renew_fc3_2
!---------------------------------------------------
    SUBROUTINE renew_fc4_2(atoms, xyzs)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: atoms(:), xyzs(:)
        !Do nothing, fc4 is already updated in the prepare_fc4

    END SUBROUTINE renew_fc4_2
!---------------------------------------------------------------------------------------------------------

    !just update fc by their value
    SUBROUTINE all_fc_update2
        IMPLICIT NONE
        INTEGER :: rnk, i, j
        INTEGER :: atom1, atom2, atom3, atom4
        INTEGER :: direction1, direction2, direction3, direction4
        INTEGER, ALLOCATABLE,DIMENSION(:) :: atoms, xyzs
        INTEGER :: new_atom1, new_atom2, new_atom3
        INTEGER,ALLOCATABLE,DIMENSION(:) :: directions
        REAL(8) :: temp

        !-------------------fc4-------------------------
        !nothing needs to be done
        !------------------fc3-------------------------
        rnk = 3
        ALLOCATE(atoms(rnk),xyzs(rnk))
        DO j=1, map(rnk)%ngr
            DO i=1, map(rnk)%nt(j)
                atoms(1) = map(rnk)%gr(j)%iat(1,i)
                IF(atoms(1).gt.atom_number) CYCLE
                atoms(2) = map(rnk)%gr(j)%iat(2,i)
                atoms(3) = map(rnk)%gr(j)%iat(3,i)

                xyzs(1) = map(rnk)%gr(j)%ixyz(1,i)
                xyzs(2) = map(rnk)%gr(j)%ixyz(2,i)
                xyzs(3) = map(rnk)%gr(j)%ixyz(3,i)

                CALL uni_fc_update2(rnk, atoms, xyzs)
            END DO
        END DO
        DEALLOCATE(atoms,xyzs)

        DO i=1, SIZE(myfc3_index)
            atom1=myfc3_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2=myfc3_index(i)%jatom_number
            atom3=myfc3_index(i)%katom_number
            direction1=myfc3_index(i)%iatom_xyz
            direction2=myfc3_index(i)%jatom_xyz
            direction3=myfc3_index(i)%katom_xyz

            IF(ANY(fc3_unique_idx==atom2).AND.ANY(fc3_unique_idx==atom3)) THEN
                atom2 = find_loc(fc3_unique_idx,atom2)
                atom3 = find_loc(fc3_unique_idx,atom3)
                myfc3_index(i)%psi_temp = &
                &myfc3_value(atom1,atom2,atom3)%psi(direction1,direction2,direction3)
            ELSE
                myfc3_index(i)%psi_temp = 0d0
            END IF
        END DO
        WRITE(34,*)"fc3 update finished!"
        !------------------fc2-------------------------
        rnk = 2
        ALLOCATE(atoms(rnk),xyzs(rnk))
        DO j=1, map(rnk)%ngr
            DO i=1, map(rnk)%nt(j)
                atoms(1) = map(rnk)%gr(j)%iat(1,i)
                IF(atoms(1).gt.atom_number) CYCLE
                atoms(2) = map(rnk)%gr(j)%iat(2,i)

                xyzs(1) = map(rnk)%gr(j)%ixyz(1,i)
                xyzs(2) = map(rnk)%gr(j)%ixyz(2,i)

                CALL uni_fc_update2(rnk, atoms, xyzs)
            END DO
        END DO
        DEALLOCATE(atoms,xyzs)

        DO i=1, SIZE(myfc2_index)
            atom1=myfc2_index(i)%iatom_number
            IF(atom1.gt.atom_number) CYCLE
            atom2=myfc2_index(i)%jatom_number
            direction1=myfc2_index(i)%iatom_xyz
            direction2=myfc2_index(i)%jatom_xyz
            myfc2_index(i)%phi_temp = myfc2_value(atom1,atom2)%phi(direction1,direction2)
        END DO

        WRITE(34,*)"fc2 update finished!"
        WRITE(*,*)"=======fc2 update checked======="
        CALL fc3_update_check
        CALL fc4_update_check
    END SUBROUTINE all_fc_update2
!==============================================================================================
!    !keep original atom list, don't call fcinit or setup_maps
!    !don't update grouping of fcs, just their values
!    SUBROUTINE half_Update
!        IMPLICIT NONE
!        CALL trans_Update
!        CALL convs_Update(conv_to_cart)
!        CALL latticeparameters_Update(conv_to_cart)
!        CALL primitivelattices_Update(conv_to_cart)
!        CALL atompos0_Update(conv_to_cart)
!        CALL atompos_Update2
!        CALL all_fc_update2
!        WRITE(34,*)'====================FINISH UPDATE ROUTINES==================='
!    END SUBROUTINE half_Update
!========================================================================================================

END MODULE check
