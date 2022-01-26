MODULE Fourier_force_constants
!!This module is for Fourier fcs calculation
    USE DFT_force_constants !this DFT... module itself uses structure_info module
    USE kp
    USE born

    IMPLICIT NONE
!---------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
    TYPE ffc2_value
        COMPLEX(8),DIMENSION(d,d):: FFC_D
    END TYPE

    TYPE ffc3_value
        COMPLEX(8),DIMENSION(d,d,d) :: FFC_C
    END TYPE

    TYPE ffc4_value
        COMPLEX(8),DIMENSION(d,d,d,d) :: FFC_B
    END TYPE
!-----------------------------------------------------------------------------------
    TYPE(ffc2_value),DIMENSION(:,:,:),ALLOCATABLE :: myffc2_value !index atom1_type, atom2_type, kvector
    TYPE(ffc3_value),DIMENSION(:,:,:,:),ALLOCATABLE :: myffc3_value !index atom1_label_number, atom2_type,atom3_type, kvector
    TYPE(ffc4_value),DIMENSION(:,:,:,:,:),ALLOCATABLE :: myffc4_value

    TYPE(ffc2_value),DIMENSION(:,:),ALLOCATABLE :: skffc2_value !dynmat value at a single k-point
    COMPLEX(8),ALLOCATABLE :: skffc2_matrix(:,:) !dynamat at a single k-point
!---------------------------for connector with legacy codes--------------------------
    COMPLEX(8), ALLOCATABLE :: dynmat(:,:),ddyn(:,:,:)
    COMPLEX(8), ALLOCATABLE :: dynmat_record(:,:,:)
    COMPLEX(8), ALLOCATABLE :: nself(:,:),uself(:,:)
    CONTAINS
!====================================================================================================
    SUBROUTINE Allocate_FFC(k_number)
    !! Allocate fourier force constant variables
        IMPLICIT NONE

        INTEGER,INTENT(IN):: k_number
        INTEGER :: tau1,tau2,tau3,tau4,atom1,atom2
        INTEGER :: k
        INTEGER :: ndn

        ALLOCATE(myffc2_value(atom_number,atom_number,k_number))

!no need to calculate fourier force constants in this method

        !ALLOCATE(myffc3_value(tot_atom_number,atom_number,atom_number,k_number))


        !ALLOCATE(myffc4_value(tot_atom_number,tot_atom_number,atom_number,atom_number,k_number))
        ndn = d*atom_number
        ALLOCATE(dynmat(ndn,ndn),ddyn(ndn,ndn,d))

        !-------------------below are not used-----------------
        ALLOCATE(dynmat_record(ndn,ndn,k_number))
        ALLOCATE(nself(k_number,ndn),uself(k_number,ndn))
    END SUBROUTINE Allocate_FFC
!==================================================================================================
     subroutine read_born
     !! this subroutine read born charge from params.born
    ! reads ! born ! effective ! charge ! and ! dielectric ! constants ! from ! para.born
         integer i,j,k
         real(8) asr

         open(uborn,file='params.born',status='old')
         allocate(zeu(3,3,natoms0)) !JS
         read(uborn,*) rho,born_flag   ! if 0 use default
         do i=1,3
            read(uborn,*)(epsil(i,j),j=1,3)
         end do
         do k=1,natoms0
            do i=1,3
               read(uborn,*)(zeu(i,j,k),j=1,3)
            end do
         end do

        ! enforce ASR to the charges
         do i=1,3
         do j=1,3
            asr=0
            do k=1,natoms0
               asr=asr+ zeu(i,j,k)
            end do
            if(abs(asr).gt.1d-9) then ! split the correction equally between the atoms
               asr=asr/natoms0
               zeu(i,j,k)=zeu(i,j,k)-asr
            endif
         end do
         end do


     end subroutine read_born
!====================================================================================================

!====================================================================================================
    SUBROUTINE Fourier_Force_Constants_Calculation(k_number)
    !! This subroutine Fourier transform FC2s, uncomment for FC3 and FC4(but they are not used) 
        IMPLICIT NONE
        INTEGER :: i,j,R1,R2,R3,R4,tau1,tau2,tau3,tau4,atom1,atom2,atom3,atom4
        INTEGER :: k
        INTEGER, INTENT(IN) :: k_number
        COMPLEX(8) :: dummy(d,d),dummy2(d,d,d)

!---------------------------------------------------------------------------------------------------
         DO tau1=1,atom_number
            DO tau2=1,atom_number
                DO k=1,k_number
                    myffc2_value(tau1,tau2,k)%FFC_D=0
                END DO
            END DO
        END DO


        DO k=1,k_number
        DO i=1,atom_number
            tau1=i
        DO j=1,tot_atom_number
            tau2=every_atom(j)%type_tau
            R2=every_atom(j)%type_R

            dummy=myfc2_value(tau1,j)%phi*EXP(ci*(kvector(k).dot.cell_vec(:,R2)))
            myffc2_value(tau1,tau2,k)%FFC_D=myffc2_value(tau1,tau2,k)%FFC_D+dummy/SQRT(iatom(tau1)%mass * iatom(tau2)%mass)
        END DO !j loop
        END DO !i loop
        END DO !k loop

!----------------------------------------------------------------------------------------------------
!        DO atom1=1,tot_atom_number
!            DO tau2=1,atom_number
!                DO tau3=1,atom_number
!                    DO k=1,k_number
!                        myffc3_value(atom1,tau2,tau3,k)%FFC_C=0
!                    END DO
!                END DO
!            END DO
!        END DO



!outer3: DO k=1, k_number
!inner3:     DO i=1,fc_terms(3)

!               atom1=myfc3_index(i)%iatom_number
!               atom2=myfc3_index(i)%jatom_number
!               atom3=myfc3_index(i)%katom_number
!               R1=every_atom(atom1)%type_R
!               R2=every_atom(atom2)%type_R
!               R3=every_atom(atom3)%type_R
!               tau1=every_atom(atom1)%type_tau
!               tau2=every_atom(atom2)%type_tau
!               tau3=every_atom(atom3)%type_tau

!               dummy2= myfc3_value(atom1,atom2,atom3)%psi*exp(ci*(kvector(k).dot.(cell_vec(:,R3)-cell_vec(:,R2))))/&
!               &       SQRT(iatom(tau2)%mass * iatom(tau3)%mass)
!               myffc3_value(atom1,tau2,tau3,k)%FFC_C = myffc3_value(atom1,tau2,tau3,k)%FFC_C + dummy2

!           END DO inner3
!       END DO outer3
!-------------------------------------------------------------------------------------------------------
!        DO atom1=1,tot_atom_number
!            DO atom2=1,tot_atom_number
!                DO tau3=1,atom_number
!                    DO tau4=1,atom_number
!                        DO k=1,k_number
!                            myffc4_value(atom1,atom2,tau3,tau4,k)%FFC_B=0
!                        END DO
!                    END DO
!                END DO
!            END DO
!        END DO


!outer4: DO k=1, k_number
!inner4:     DO i=1,fc_terms(4)

!               atom1=myfc4_index(i)%iatom_number
!               atom2=myfc4_index(i)%jatom_number
!               atom3=myfc4_index(i)%katom_number
!               atom4=myfc4_index(i)%latom_number
!               R1=every_atom(atom1)%type_R
!               R2=every_atom(atom2)%type_R
!               R3=every_atom(atom3)%type_R
!               R4=every_atom(atom4)%type_R
!               tau1=every_atom(atom1)%type_tau
!               tau2=every_atom(atom2)%type_tau
!               tau3=every_atom(atom3)%type_tau
!               tau4=every_atom(atom4)%type_tau

!               myffc4_value(atom1,atom2,tau3,tau4,k)%FFC_B = myffc4_value(atom1,atom2,tau3,tau4,k)%FFC_B + &
!               &                                       myfc4_value(atom1,atom2,atom3,atom4)%chi* &
!               &                                       exp(ci*(kvector(k).dot.(cell_vec(:,R4)-cell_vec(:,R3))))/&
!               &                                       SQRT(iatom(tau3)%mass * iatom(tau4)%mass)

!           END DO inner4
!       END DO outer4

!----------------------------------------------------------------------------------------------------------
    END SUBROUTINE Fourier_Force_Constants_Calculation

END MODULE Fourier_force_constants
