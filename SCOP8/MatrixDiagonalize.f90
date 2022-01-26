!! This module is for diagonalize dynamic matrix
MODULE MatrixDiagonalize
    USE eigen
    IMPLICIT NONE

    !with a letter 's', it covers all the eivec corresponding to every k point
    !so the last extra dimension tells which k point it is related with
    ! eivals and eivecs are equivalent with eigenval and eigenvec in module [eigen]
    COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE :: eivecs, eivecs_t
    REAL(8),DIMENSION(:,:),ALLOCATABLE :: eivals

    CONTAINS

!==========================================================
    !!Allocate eigen related variables
    !!eigen_number = d*atom_number = ndyn = nb
    SUBROUTINE allocate_eigen(eigen_number,k_number)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: eigen_number,k_number
        ALLOCATE(eivecs(eigen_number,eigen_number,k_number))!(direction*atom type,lambda,q)
        ALLOCATE(eivecs_t(eigen_number,eigen_number,k_number))!(lambda,direction*atom,q)
        ALLOCATE(eivals(eigen_number,k_number))
        return
    END SUBROUTINE allocate_eigen
!==========================================================
    !!Diagonalize dynamic matrix mat, calculate eigenvalues eival and eigenvectors eivec
    !!this subroutine will change the mat, after calling, mat will no longer have the original value
    !! n=size of mat; nv is the number of needed eigenvectors
    !!here the eival and eivec are dummy variables...
    !!... and only corresponds to 1 specified k point e.g. eivecs(:,:,1)
    SUBROUTINE diagonalize(n,mat,eival,nv,eivec,ier)
 
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: n,ier,nv  !,i,j
    COMPLEX(8),INTENT(INOUT) :: mat(n,n),eivec(n,n)
    REAL(8),INTENT(INOUT) ::  eival(n)

    ! This is used by ZHEGV
    REAL(8),ALLOCATABLE :: rwork(:)
    COMPLEX(8),ALLOCATABLE:: work(:)
    INTEGER :: lwork

    ! external ZHEEV

    ! n = size(mat(:,1))
    IF (n .ne. size(eival) ) THEN
        WRITE(*,*)' EIGCH, size inconsistency:mat,eival=',n, size(eival)
        STOP
    END IF

    lwork = max(1,2*n-1)
    ALLOCATE(work(lwork),rwork(max(1,3*n-2)))


    IF(nv.eq.0) THEN
        CALL zheev('N','U',n,mat,n,eival,work,lwork,rwork,ier)
    ELSE
        CALL zheev('V','U',n,mat,n,eival,work,lwork,rwork,ier)
        eivec = mat
    END IF
    DEALLOCATE(work,rwork)
!4 format(i5,1x,3(f6.3,1x,f6.3,4x))

    END SUBROUTINE diagonalize
!=============================================================================
    !! get the transposed conjugate eigenvector eivec
    SUBROUTINE dagger_eigen(eivec,eivec_t)
        IMPLICIT NONE
        COMPLEX(8),INTENT(INOUT),DIMENSION(:,:) :: eivec,eivec_t
        eivec_t = transpose(conjg(eivec))
    END SUBROUTINE dagger_eigen
!=============================================================================
    !! test subroutine for diagonalization
    SUBROUTINE test_diagonalization(s)
        IMPLICIT NONE
        INTEGER,INTENT(in) :: s
        INTEGER :: i,j,l
        COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: epsil,epsil_t
        REAL(8),DIMENSION(:),ALLOCATABLE :: omega
        COMPLEX(8),DIMENSION(:,:),ALLOCATABLE :: test_matrix
        REAL(8) :: test_sum,test_diff

        CALL srand(2021)

        ALLOCATE(omega(s),epsil(s,s),epsil_t(s,s),test_matrix(s,s))

        !randomly generate a Hermitian test matrix
        DO i=1,s
            DO j=1,i
                test_matrix(i,j) = CMPLX(i*j*rand(),(i+j)*rand())
            END DO
        END DO
        DO i=1,s-1
            DO j=i+1,s
                test_matrix(i,j) = CONJG(test_matrix(j,i))
            END DO
        END DO

        CALL diagonalize(s,test_matrix,omega,s,epsil,0)
        CALL dagger_eigen(epsil,epsil_t)

        OPEN(72,FILE='diagonalization_check.dat',STATUS='unknown',ACTION='write')
        WRITE(72,*)'i, j, discrepancy'
        DO i=1,s
        DO j=1,s
            test_sum = 0d0
            DO l=1,s !sum all lambda, for any given (i,j)
                test_sum = test_sum + epsil(i,l)*epsil_t(l,j)*omega(l)
            END DO
            test_diff = ABS(test_matrix(i,j)-test_sum)
            IF(test_diff.gt.1d-8) THEN
                WRITE(72,*)i,j,test_diff
            END IF
        END DO
        END DO

        CLOSE(72)
    END SUBROUTINE test_diagonalization

END MODULE
