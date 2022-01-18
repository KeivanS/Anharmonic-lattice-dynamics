!=========================================================================
      module broy
      use Iteration_parameters
      implicit none
      real(8), allocatable :: xin(:,:),ubroy(:,:),vbroy(:,:),fbroy(:,:)
      integer itmax
      integer, allocatable :: effective_index(:)
      real(8) pmix

      contains
!------------------------------------------------------------------------
        subroutine allocatebroy(n,m)
        implicit none
        integer, intent(in) :: n,m !n is size of array, m is the iteration number

        if(allocated(xin)) deallocate(xin)
        if(allocated(ubroy)) deallocate(ubroy)
        if(allocated(vbroy)) deallocate(vbroy)
        if(allocated(fbroy)) deallocate(fbroy)
        allocate ( xin(n,m),ubroy(n,m),vbroy(n,m),fbroy(n,m) )
        !call drop_fc2terms
        !allocate(effective_index(indie_fc2))
        return
        end subroutine allocatebroy
!------------------------------------------------------------------------
        subroutine deallocatebroy
        implicit none

        deallocate (xin,ubroy,vbroy,fbroy)
        !call drop_fc2terms
        !deallocate(effective_index)
        return
        end subroutine deallocatebroy
!------------------------------------------------------------------------
        subroutine combine_variants(x)
            !mx is the size of x, it equals the size of GradientF_trial, which is f
            !x is the: atomic_deviation(d,atom_number),strain(d,d),yy_value(atom_number,tot_atom_number)%phi(d,d)
            !...combined into a single rank array
            implicit none
            integer checkcombine
            integer i,j,mx,temp,temp1,temp2,tau1,atom2,direction1,direction2
            real(8),intent(out),dimension(:),allocatable::x !it's already allocated though

!******************************
checkcombine=14
OPEN(checkcombine,file='combine.dat',STATUS='unknown',ACTION='write')
!******************************
            mx=SUM(variational_parameters_size)
            allocate(x(mx))

WRITE(checkcombine,*) '===================atomic deviation====================='

            do i=1,variational_parameters_size(1) !u_tau
                IF(MOD(i,d).eq.0) THEN
                    x(i)=atomic_deviation(d,INT(i/d))
                ELSE
                    x(i)=atomic_deviation(MOD(i,d),INT(i/d+1))
                END IF
WRITE(checkcombine,*) x(i)
            end do
WRITE(checkcombine,*) '===================strain tensor============================'

            do i=variational_parameters_size(1)+1,mx-variational_parameters_size(3) !eta strain
                temp=i-variational_parameters_size(1)
                IF(MOD(temp,d).eq.0) THEN
                    x(i)=strain(INT(temp/d),d)
                ELSE
                    x(i)=strain(INT(temp/d+1),MOD(temp,d))
                END IF
WRITE(checkcombine,*) x(i)
            end do
WRITE(checkcombine,*) '=========================<YY>========================='

           do i=mx-variational_parameters_size(3)+1,mx
                temp=i-variational_parameters_size(1)-variational_parameters_size(2)

                tau1=eff_fc2_index(temp)%iatom_number
                atom2=eff_fc2_index(temp)%jatom_number
                direction1=eff_fc2_index(temp)%iatom_xyz
                direction2=eff_fc2_index(temp)%jatom_xyz

                x(i)=trialfc2_value(tau1,atom2)%phi(direction1,direction2)
WRITE(checkcombine,4) tau1,direction1,atom2,direction2,trialfc2_value(tau1,atom2)%phi(direction1,direction2)
WRITE(checkcombine,5) i,x(i)
           end do
close(checkcombine)
4 FORMAT(4I4,F10.5)
5 FORMAT(I4,F10.5)
        end subroutine combine_variants
!------------------------------------------------------------------------
        subroutine decompose_variants(x)
            !it's a inverse process of combine_variants(x)
            !after one call of bro90, the x is updated. when need to decompose and distribute x, for VA_math to get f
            implicit none
            INTEGER :: checkdecompose,flag(2), checknewfc2
            integer i,j,k,l
            integer mx,temp,temp1,temp2,tau1,atom2,direction1,direction2
            real(8) :: weight(2)
            real(8),intent(in),dimension(:) :: x
!*************
checkdecompose=15
OPEN(checkdecompose,file='decompose.dat',status='unknown',action='write')
!*************
            mx=size(x)
WRITE(checkdecompose,*) '===================atomic deviation====================='
            do i=1,variational_parameters_size(1)
                IF(MOD(i,d).eq.0) THEN
                    atomic_deviation(d,INT(i/d))=x(i)
                ELSE
                    atomic_deviation(MOD(i,d),INT(i/d+1))=x(i)
                END IF
WRITE(checkdecompose,*) x(i)
            end do

WRITE(checkdecompose,*) '===================strain tensor====================='

            do i=variational_parameters_size(1)+1,variational_parameters_size(1)+variational_parameters_size(2)
                temp=i-variational_parameters_size(1)
                IF(MOD(temp,d).eq.0) THEN
                    strain(INT(temp/d),d)=x(i)
                ELSE
                    strain(INT(temp/d+1),MOD(temp,d))=x(i)
                END IF

WRITE(checkdecompose,*) x(i)

            end do

WRITE(checkdecompose,*) '===================trial force====================='
            do i=variational_parameters_size(1)+variational_parameters_size(2)+1,mx
                temp=i-variational_parameters_size(1)-variational_parameters_size(2)

                tau1=eff_fc2_index(temp)%iatom_number
                atom2=eff_fc2_index(temp)%jatom_number
                direction1=eff_fc2_index(temp)%iatom_xyz
                direction2=eff_fc2_index(temp)%jatom_xyz

                !save the trialfc2 value from last iteration first
                prev_trialfc2_value(tau1,atom2)%phi(direction1,direction2)=&
                &trialfc2_value(tau1,atom2)%phi(direction1,direction2)
                !then renew the trialfc2 value using x(:)
                trialfc2_value(tau1,atom2)%phi(direction1,direction2)=x(i)
WRITE(checkdecompose,4) tau1,atom2,get_letter(direction1),get_letter(direction2),&
&trialfc2_value(tau1,atom2)%phi(direction1,direction2)
WRITE(checkdecompose,5) i,x(i)

            end do

close(checkdecompose)
4 FORMAT(2I4,2A2,F10.5)
5 FORMAT(I4,F10.5)

        end subroutine decompose_variants
!====================================================================================================
        subroutine decompose_variants2(x)
            !it's a inverse process of combine_variants(x)
            !after one call of bro90, the x is updated. when need to decompose and distribute x, for VA_math to get f
            implicit none
            INTEGER :: checkdecompose,flag(2), checknewfc2
            integer i,j,k,l
            integer mx,temp,temp1,temp2,tau1,atom1,atom2,direction1,direction2
            integer rnk,ntindp
            real(8) :: weight(2)
            real(8),intent(in),dimension(:) :: x
!*************
checkdecompose=15
OPEN(checkdecompose,file='decompose.dat',status='unknown',action='write')
checknewfc2=16
OPEN(checknewfc2,file='newfc2.dat',status='unknown',action='write')
!*************
            mx=size(x)
WRITE(checkdecompose,*) '===================atomic deviation====================='
            do i=1,variational_parameters_size(1)
                IF(MOD(i,d).eq.0) THEN
                    atomic_deviation(d,INT(i/d))=x(i)
                ELSE
                    atomic_deviation(MOD(i,d),INT(i/d+1))=x(i)
                END IF
WRITE(checkdecompose,*) x(i)
            end do

WRITE(checkdecompose,*) '===================strain tensor====================='

            do i=variational_parameters_size(1)+1,variational_parameters_size(1)+variational_parameters_size(2)
                temp=i-variational_parameters_size(1)
                IF(MOD(temp,d).eq.0) THEN
                    strain(INT(temp/d),d)=x(i)
                ELSE
                    strain(INT(temp/d+1),MOD(temp,d))=x(i)
                END IF
WRITE(checkdecompose,*) x(i)
            end do

WRITE(checkdecompose,*) '===================trial force====================='
            do i=variational_parameters_size(1)+variational_parameters_size(2)+1,mx
                temp=i-variational_parameters_size(1)-variational_parameters_size(2)
                tau1=indiefc2_index(temp)%iatom_number
                atom2=indiefc2_index(temp)%jatom_number
                direction1=indiefc2_index(temp)%iatom_xyz
                direction2=indiefc2_index(temp)%jatom_xyz
                trialfc2_value(tau1,atom2)%phi(direction1,direction2)=x(i)
WRITE(checkdecompose,4) tau1,atom2,get_letter(direction1),get_letter(direction2),&
&trialfc2_value(tau1,atom2)%phi(direction1,direction2)
WRITE(checkdecompose,5) i,x(i)
            end do
WRITE(checkdecompose,*) '================assign other trial force===================='
            !--------------IMPORTANT:update the rest of trialfc2_value----------------------
            !*NEW* method to distribute
            rnk = 2
            DO temp=variational_parameters_size(1)+variational_parameters_size(2)+1,mx
               j=temp-variational_parameters_size(1)-variational_parameters_size(2)
               ntindp=map(rnk)%ntind(j)
               DO i=1, map(rnk)%nt(j)
                    atom1 = map(rnk)%gr(j)%iat(1,i)
                    IF(atom1.gt.atom_number) CYCLE
                    atom2 = map(rnk)%gr(j)%iat(2,i)
                    direction1 = map(rnk)%gr(j)%ixyz(1,i)
                    direction2 = map(rnk)%gr(j)%ixyz(2,i)
                    trialfc2_value(atom1,atom2)%phi(direction1,direction2)=x(temp)
                    IF(ntindp.eq.1) THEN
                        trialfc2_value(atom1,atom2)%phi(direction1,direction2)=&
                        &trialfc2_value(atom1,atom2)%phi(direction1,direction2)&
                        &*map(rnk)%gr(j)%mat(i,1) !multiply weight factor
                    END IF

WRITE(checkdecompose,4) atom1,atom2,get_letter(direction1),get_letter(direction2),&
&trialfc2_value(atom1,atom2)%phi(direction1,direction2)
WRITE(checkdecompose,5) temp,x(temp)

               END DO
            END DO

            !this is not complete, but missing trialfc2(atom,atom)%phi(direction1,direction2)
            !however, it fill be filled up by CheckFixASR

            !force symmetry on strain(:,:)?

            !call CheckFixASR
close(checkdecompose)
close(checknewfc2)
4 FORMAT(2I4,2A2,F10.5)
5 FORMAT(I4,F10.5)

        end subroutine decompose_variants2

!------------------------------------------------------------------------
!this subroutine is to check if the trial fc2 is recovered correctly from
!subroutine <decompose_variants>
SUBROUTINE check_decompose
    IMPLICIT NONE
    INTEGER :: atom1,atom2,direction1,direction2
    INTEGER :: i,j,ntindp,rnk
    INTEGER :: checkK

    checkK = 21
    OPEN(checkK,FILE='decompose_check_K.dat',STATUS='unknown',ACTION='write',POSITION='append')
    WRITE(checkK,*)'==========================================='
    rnk = 2
    DO j=1, map(rnk)%ngr
        ntindp = map(rnk)%ntind(j)
        WRITE(checkK,*)'----group',j,'-----'
        DO i=1, map(rnk)%nt(j)
            atom1 = map(rnk)%gr(j)%iat(1,i)
            IF(atom1.gt.atom_number) CYCLE
            atom2 = map(rnk)%gr(j)%iat(2,i)
            direction1 = map(rnk)%gr(j)%ixyz(1,i)
            direction2 = map(rnk)%gr(j)%ixyz(2,i)
            WRITE(checkK,'(f8.5,f7.3)') trialfc2_value(atom1,atom2)%phi(direction1,direction2),&
            &map(rnk)%gr(j)%mat(i,1:ntindp)
        END DO
    END DO
    CLOSE(checkK)
END SUBROUTINE check_decompose
!------------------------------------------------------------------------
        !this subroutine extend f(:) used in Broyden by
        !the indiefc2 <-----> extendfc2 relation
        subroutine extend_f(arrayIn,arrayOut)
            IMPLICIT NONE
            INTEGER :: atoms(2),xyzs(2)
            TYPE(fc2_index),DIMENSION(:),ALLOCATABLE :: ex_indiefc2_index
            REAL(8),DIMENSION(:),INTENT(IN) :: arrayIn
            REAL(8),DIMENSION(:),INTENT(OUT),ALLOCATABLE :: arrayOut

            INTEGER :: i, mx, temp,idx

            CALL Get_FC2pairs
            CALL Extend_indiefc2(ex_indiefc2_index)
            mx=variational_parameters_size(1)+variational_parameters_size(2)+&
            &SIZE(ex_indiefc2_index) !size has changed to a uniform number
            ALLOCATE(arrayOut(mx))

            !strain and atomic deviation part can be copied directly
            DO i=1,variational_parameters_size(1)+variational_parameters_size(2)
                arrayOut(i) = arrayIn(i)
            END DO

            !fc2 part need to be extended
            temp = variational_parameters_size(1)+variational_parameters_size(2)
            DO i=1,SIZE(ex_indiefc2_index)
                atoms(1)=ex_indiefc2_index(i)%iatom_number
                atoms(2)=ex_indiefc2_index(i)%jatom_number
                xyzs(1)=ex_indiefc2_index(i)%iatom_xyz
                xyzs(2)=ex_indiefc2_index(i)%jatom_xyz
                !index for arrayOut
                idx = variational_parameters_size(1)+variational_parameters_size(2)+i

                IF(find_indiefc2(atoms,xyzs)) THEN
                    !index for arrayIn, shift by one if found
                    temp = temp + 1
                    !here it is not very robust, I use the rule that
                    !indie fc2 is always in correct order:
                    !atomic number label(from small to big)
                    !   directional label(xx->xy->xz->yy->yz->zz)
                    arrayOut(idx) = arrayIn(temp)

                ELSE
                    arrayOut(idx) = 0d0 !if not found, set to 0
                END IF
            END DO

            !for check
!            WRITE(34,*) '-----f(:) before extended-----'
!            DO i=1,SIZE(arrayIn)
!                WRITE(34,*)i,arrayIn(i)
!            END DO
!            WRITE(34,*)'-----f(:) after extended-----'
!            DO i=1,SIZE(arrayOut)
!                WRITE(34,*)i,arrayOut(i)
!            END DO

            DEALLOCATE(ex_indiefc2_index)
        end subroutine extend_f
!------------------------------------------------------------------------
        subroutine calf(x,f,m)
        ! this is the subroutine which calculates the array f for the input array x
        ! in the actual code, it can have any name, but the 3 arguments must be as (x,f,m)
        ! where both arrays x and f are of dimension m,
        ! the same name should be used in the calling routine (here testbro) within the dowhile loop above
        implicit none
        integer i, m
        real(8),intent(inout)::  x(m),f(m)

    ! this is an example; the force includes also a non-linear term.
    ! in VA codes, f(m)->GradientF_trial(3*tot_atom_number,3*tot_atom_number,3)
        do i=1,m
            f(i) = i**2 * x(i) + i*0.01*x(i)**3
        enddo

        end subroutine calf
!-----------------------------------------------------------------------
        subroutine bro90(mx,x,f,it)
!
! uses the broyden method to make fbroy(xin)=0 at the next few
! iterations. It starts with a linear mixing scheme with mix,
! it gives the new input value for the next iteration
! It basically solves for rhoin, to make fbroy equal to 0
! mx is the actual dimension of the x and f arrays
!
      !use broy
      implicit none
      integer i,j,l,it,jmin,mx
      real(8) x(mx),f(mx),sum,s(itmax),vtmp(mx),ftmp(mx) !s size is not good here
!
      jmin = max((it-40),1)!prev: 40 coefficient can be lowered

!      jmin = 1
      if (it.ge.itmax) then
         print*,' In BROYDEN, it=',it,' is greater than ',itmax
         return
      endif
!
! store in new arrays for iterations
      xin  (1:mx,it) = x
      fbroy(1:mx,it) = f

! it=1
      if( it.eq.1 ) then
          xin(:,it+1) = xin(:,it) + pmix * fbroy(:,it)
          write(*,*)'BRO90: iteration 1 done; return!'
          x = xin(1:mx,it+1)
          return
      endif

! sum in the V vector and create V vector
!
      ftmp(1:mx) = fbroy(1:mx,it)-fbroy(1:mx,it-1)
      sum=dot_product(ftmp,ftmp)
      vbroy(1:mx,it-1) = ftmp(1:mx)/sum

!     write(*,5)(vbroy(i,it-1),i=1,n)
 5    format(5(2x,g10.4))
!
! sum in the U vector and create U vector
      do j=jmin,it-2
         s(j) = 0.0d0
         vtmp(1:mx) = vbroy(1:mx,j)
         s(j) = dot_product(vtmp,ftmp)
      enddo

      do i=1,mx
         ubroy(i,it-1) = -pmix *(fbroy(i,it)-fbroy(i,it-1))
         do j=jmin,it-2
            ubroy(i,it-1) = ubroy(i,it-1) - s(j)*(ubroy(i,j) - (xin(i,j+1)-xin(i,j)) )
         enddo
      enddo
!
      do j=jmin,it-1
         s(j) = 0.0d0
         vtmp(1:mx) = vbroy(1:mx,j)
         ftmp(1:mx) = fbroy(1:mx,it)
         s(j) = dot_product(vtmp,ftmp)
      enddo!

      do i=1,mx
         xin(i,it+1) = xin(i,it) + pmix * fbroy(i,it)
         do j=jmin,it-1
            xin(i,it+1) = xin(i,it+1) + s(j)*(ubroy(i,j) - (xin(i,j+1)-xin(i,j)) )
         enddo
      enddo

      x = xin(1:mx,it+1)

      end  subroutine bro90
!----------------------------------------------------------------------------------------

      end module broy
!----------------------------------------------------------------------


