 module linalgb
 use constants, only : r15

  interface append_array
     module procedure append_array_1d,append_array_2d
  end interface

  interface symmetrize_res
     module procedure symmetrize_res1,symmetrize_res2
  end interface

 contains

 subroutine append_array_1d(a,b,c)
!! appends array b to the end of a and stores the result in c
 implicit none
! integer, intent(in) :: na,nb
! real(r15), intent(in):: a(na),b(nb)
! real(r15), allocatable, intent(inout):: c(:)
! real(r15) :: c(size(a)+size(b))
 real(r15), intent(in):: a(:),b(:)
 real(r15), allocatable :: c(:),aux(:)
 integer nc,na,nb

 na=size(a(:));
 nb=size(b(:));
 nc=na+nb
 if (allocated(c)) deallocate (c)
 allocate(c(nc))
 c=reshape(a,shape=(/nc/),pad=b)

 return

! could also use
! allocate(aux(nc)) ! this is to allow calls like append_array(a,b,a)
! aux(1:na)=a
! aux(na+1:na+nb)=b
! if (allocated(c)) deallocate (c)
! allocate(c(nc))
! c=aux
! deallocate(aux)

 end subroutine append_array_1d
!---------------------------------------
 subroutine append_array_2d(a,b,c)
!! appends array b to the end of a and stores the result in c
 implicit none
! integer, intent(in) :: na,nb
! real(r15), intent(in):: a(na),b(nb)
! real(r15), allocatable, intent(inout):: c(:)
! real(r15) :: c(size(a)+size(b))
 real(r15), intent(in):: a(:,:),b(:,:)
 real(r15), allocatable :: c(:,:) ,aux(:,:)
 integer nc,na,nb,col

 col=size(a,2) ! (1,:));
 na =size(a,1)   ! (:,1));
 nb =size(b,1)   ! (:,1));
 nc =na+nb

! c=reshape(transpose(a),shape=(/col,nc/),pad=transpose(b),order=(/1,2/))
! c=transpose(c)

! return

! could also use
 allocate(aux(nc,col)) ! this is to allow calls like append_array(a,b,a)
 aux(1:na,:)=a
 aux(na+1:na+nb,:)=b
 if (allocated(c)) deallocate (c)
 allocate(c(nc,col))
 c=aux
 deallocate(aux)

 end subroutine append_array_2d
!---------------------------------------
 subroutine symmetrize_res1(mat2,res)
! enforces mat2(a,b)-mat2(b,a)=res(b); symmetrize mat then mat2=sym(mat)+res/2 and sym(mat)-res/2
 use constants, only : r15
 implicit none
 real(r15), intent(inout) :: mat2(3,3) 
 real(r15), intent(in) :: res(3) 
 real(r15) mean(3,3) 
 integer n

 mean=(mat2+transpose(mat2))/2
 mat2(1,2)=mean(1,2)+res(2)/2
 mat2(2,1)=mean(1,2)-res(2)/2
 mat2(2,3)=mean(3,2)+res(3)/2
 mat2(3,2)=mean(3,2)-res(3)/2
 mat2(3,1)=mean(1,3)+res(1)/2
 mat2(1,3)=mean(1,3)-res(1)/2
 
 end subroutine symmetrize_res1
!---------------------------------------
 subroutine symmetrize_res2(mat2,res)
! enforces mat2(a,b)-mat2(b,a)=res(a,b)-res(b,a); consider mat=mat2-res ; symmetrize mat then mat2=sym(mat)+res
 use constants, only : r15
 implicit none
 real(r15), intent(inout) :: mat2(:,:) 
 real(r15), intent(in) :: res(:,:) 
 real(r15), allocatable:: mean(:,:) 
 integer n

 n=size(mat2,1)
 allocate(mean(n,n))

 mean=(mat2-res+transpose(mat2-res))/2
 mat2=mean+res
 deallocate(mean)
 end subroutine symmetrize_res2

!---------------------------------------
 subroutine symmetrize2(n,mat2)
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(inout) :: mat2(n,n) 
 real(r15), allocatable:: mean(:,:) 
 real(r15) dif 

! n=size(mat2,1)
 allocate(mean(n,n))
 dif=maxval(abs(mat2-transpose(mat2)))
! relative (user's rule 2026-11-07): warn if max|M-M^T| > 1e-6 max|M| (was 1e-4 absolute)
 if(dif.gt.1d-6*max(maxval(abs(mat2)),tiny(1d0))) write(*,*)'SYMMETRIZE2: max|M-M^T|/max|M|=',dif/maxval(abs(mat2))
 mean=(mat2+transpose(mat2))/2
 mat2=mean
 deallocate(mean)
 end subroutine symmetrize2
!---------------------------------------
 subroutine lin_solve(a,b,x)
 use constants, only : r15
 real(r15), intent(in) :: a(:,:),b(:) 
 real(r15), intent(out) :: x(:)
! integer, parameter :: n = size(A,1)
 real(r15) factor,sum
 real(r15), allocatable:: aux(:,:),bb(:)
 integer i,j,k,ier,n

 n =size(a,1)
! stupid way!!
 allocate(aux(n,n),bb(n)) 
! call xmatinv(n,a,ainv,ier)
! if(ier.eq.0)  then
!    x=matmul(ainv,b)
! else
!    write(*,*)' output of xmatinv from lin_solve has error ',ier
!    stop
! endif

 aux=a; bb=b
        ! Forward Elimination
 do i=1,n-1
     do k=i+1,n
         ! Calculate elimination factor
         factor = aux(k,i)/aux(i,i)  

         ! Modify row k to subtract a multiple of row i
         do j=i+1,n
             aux(k,j) = aux(k,j) - factor*aux(i,j)
         enddo
         bb(k) = bb(k) - factor*bb(i)
     enddo
 enddo

 ! Back Substitution
 x(n) = bb(n) / aux(n,n)
 do i=n-1,1,-1
     sum = 0
     do j=i+1,n
         sum = sum + aux(i,j)*x(j)
     enddo
     x(i) = (bb(i) - sum)/aux(i,i)
 enddo

 deallocate(aux,bb)

 end subroutine lin_solve

!---------------------------------------
 subroutine lin_solve33(a,b,x)
 use constants, only : r15
 real(r15), intent(in) :: a(3,3),b(3) 
 real(r15), intent(out) :: x(3)
 real(r15) factor,deta,d1,d2,d3,det
 real(r15) aux(3,3),bb(3)
 integer i,j,k,ier,n

 n = 3
 deta=det33(a)
 if(abs(deta).gt.1d-9) then
    aux=a ; aux(:,1)=b; x(1)=det33(aux)/deta
    aux=a ; aux(:,2)=b; x(2)=det33(aux)/deta
    aux=a ; aux(:,3)=b; x(3)=det33(aux)/deta
 else
    write(*,*)' lin_solve33 error : det=',deta
    stop
 endif

 end subroutine lin_solve33
!------------------------------------------------------------------------------
      function det33(mat) result(y)
 use constants, only : r15
      implicit none
!
!!	FIND THE DETERMINANT OF A 3X3 MATRIX MAT
!
      real(r15) ,intent(in) :: mat(3,3)
      real(r15) y

      y=mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))  &
     & -mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1))  &
     & +mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))

      end function det33
!---------------------------------------
 end module linalgb


