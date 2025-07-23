! program test_cg
! implicit none
! integer, parameter :: n=20
! real(8) a(n,n),b(n),x0(n),xmin(n),funcmin
! integer i,j
!
! !  CALL init_random_seed()
!
!    call random_number(a)
!    call random_number(b)
!    call random_number(x0)
!    do i=1,n
!       a(i,i)=a(i,i)+1+i*i
!    enddo
!    call cg_quadratic(n,a,b,x0,funcmin,100,4)
! 
! end program test_cg
!================================================
 subroutine cg_quadratic(n,a,b,vmin,funcmin,max_iter,n_dig_accur)
! this routine finds the minimum of 1/2 v.a.v - b.v wrt v for given a and b
! n_dig_accur is the requested number of digits accuracy
! initial guess v0 is CHANGED on output
 implicit none
 integer n,i,iter,n_dig_accur,max_iter
 real(8) tolerance,err,v2,vmin(n),b(n),a(n,n),num,den,funcmin
 real(8), allocatable :: conj(:),grad(:),aux(:),v0(:)

 allocate( conj(n),grad(n),aux(n),v0(n))

 tolerance=1d0/10d0**(2*n_dig_accur)
 write(*,*)'CG: tolerance is=',tolerance
 iter=0 ; err=9999
 v0=vmin
 grad=matmul(a,v0)-b
 conj=-grad
 do while ( iter.lt.max_iter )  
   iter=iter+1 
   aux=matmul(a,conj)
   num=dot_product(grad,conj)
   den=dot_product(aux,conj)
   v0=grad  ! save grad of the previous iteration in the dummy variable v0
   vmin=vmin-num/den*conj   ! this is the new minimum
   grad=grad-num/den*aux    ! now update grad

! we define the error to be (deltav/v)**2
   v2=dot_product(vmin,vmin)
   err=(num*num/den/den)*dot_product(conj,conj)/v2 ! error is (delta_v/v)^2
   if ( err .lt. tolerance) exit  ! this gives about 4,5 digits accuracy in the answer

   num=dot_product(grad,grad)
   den=dot_product(v0,v0)
   conj=-grad+num/den*conj
   write(*,3)'CG: iter,err,vmin^2,grad^2=',iter,err,v2,num
!  write(10,2)vmin,matmul(a,vmin)-b
 enddo

 v0= matmul(a,vmin)
 funcmin=0.5*dot_product(vmin,v0)-dot_product(vmin,b)
! aux=v0-b   ! this is the actual value of Ax-b
 write(*,3)'CG: loop ended; iter, error, funcmin,|Ax-b|=',iter,err,funcmin
! write(*,1)'vmin=',vmin
! write(*,1)'Av-b=',aux

1 format(a,99(1x,g10.4))
2 format(99(1x,g10.4))
3 format(a,i5,9(2x,g10.4))
 deallocate(aux,conj,grad,v0)

 end subroutine cg_quadratic
