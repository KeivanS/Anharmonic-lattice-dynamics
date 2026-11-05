!============================================================
 module broyden_mix
!! Modified Broyden method for f(x) = 0: bro90 of SCOP8/Broyden.f95, copied 2026-11-04 (user: "Broyden's method is
!! implemented in SCOP8 and you can use that subroutine"; SCOP8 itself does not compile). The algorithm is unchanged;
!! changed: no module Iteration_parameters, arrays allocated by broyden_init, no screen output.
!! Each call takes x_it and f_it = f(x_it) and returns x_{it+1} in x. The first step is x + pmix f, so f must be a
!! step direction (relax_reference passes minus the preconditioned gradient, with pmix = 1).
 use constants, only : r15
 implicit none
 real(r15), allocatable :: xin(:,:),ubroy(:,:),vbroy(:,:),fbroy(:,:)
 integer itmax
 real(r15) pmix
 contains

 subroutine broyden_init(n,m,pm)
 integer, intent(in) :: n,m   ! size of x, maximum number of iterations
 real(r15), intent(in) :: pm
 if(allocated(xin)) deallocate(xin,ubroy,vbroy,fbroy)
 allocate(xin(n,m),ubroy(n,m),vbroy(n,m),fbroy(n,m))
 xin=0; ubroy=0; vbroy=0; fbroy=0
 itmax=m; pmix=pm
 end subroutine broyden_init

 subroutine broyden_free
 if(allocated(xin)) deallocate(xin,ubroy,vbroy,fbroy)
 end subroutine broyden_free

 subroutine bro90(mx,x,f,it)
!! uses the broyden method to make fbroy(xin)=0 at the next few iterations. It starts with a linear mixing
!! scheme with pmix and gives the new input value for the next iteration. mx is the dimension of x and f.
 integer, intent(in) :: mx,it
 real(r15), intent(inout) :: x(mx)
 real(r15), intent(in) :: f(mx)
 integer i,j,jmin
 real(r15) sm,s(itmax),vtmp(mx),ftmp(mx)
 jmin = max((it-40),1)
 if (it.ge.itmax) return
 xin  (1:mx,it) = x
 fbroy(1:mx,it) = f
 if( it.eq.1 ) then
    xin(:,it+1) = xin(:,it) + pmix * fbroy(:,it)
    x = xin(1:mx,it+1)
    return
 endif
 ftmp(1:mx) = fbroy(1:mx,it)-fbroy(1:mx,it-1)
 sm=dot_product(ftmp,ftmp)
 vbroy(1:mx,it-1) = ftmp(1:mx)/sm
 do j=jmin,it-2
    vtmp(1:mx) = vbroy(1:mx,j)
    s(j) = dot_product(vtmp,ftmp)
 enddo
 do i=1,mx
    ubroy(i,it-1) = -pmix *(fbroy(i,it)-fbroy(i,it-1))
    do j=jmin,it-2
       ubroy(i,it-1) = ubroy(i,it-1) - s(j)*(ubroy(i,j) - (xin(i,j+1)-xin(i,j)) )
    enddo
 enddo
 do j=jmin,it-1
    vtmp(1:mx) = vbroy(1:mx,j)
    ftmp(1:mx) = fbroy(1:mx,it)
    s(j) = dot_product(vtmp,ftmp)
 enddo
 do i=1,mx
    xin(i,it+1) = xin(i,it) + pmix * fbroy(i,it)
    do j=jmin,it-1
       xin(i,it+1) = xin(i,it+1) + s(j)*(ubroy(i,j) - (xin(i,j+1)-xin(i,j)) )
    enddo
 enddo
 x = xin(1:mx,it+1)
 end subroutine bro90

 end module broyden_mix
