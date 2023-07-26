!==============================================================
 module constants
 implicit none
!--! specific precisions, usually same as real and double precision
integer, parameter :: r6 = selected_real_kind(6) 
integer, parameter :: r15 = selected_real_kind(15) ! should work on any compiler
 integer, parameter :: dp=8 !(for nag compilers use kd=2 for double precision; 8 for gfortran)
 complex(kind=r15) :: ci=cmplx(0d0,1d0)
! complex(8) :: ci=cmplx(0d0,1d0)
 real(kind=r15) :: pi=3.14159265358979323846d0
 real(r15) :: h_plank= 6.62606896d-34
 real(r15) :: n_avog= 6.023d23
 real(r15) :: k_b = 1.3806504d-23       ! J/K
 real(r15) :: c_light = 2.99792458d+8   ! in (m/s)
 real(r15) :: hbar = 1.054571628d-34  !h_plank/2/pi
 real(r15) :: ee = 1.60217653d-19
 real(r15) :: eps0 = 8.854187817d-12
 real(r15) :: me = 9.1093826d-31
 real(r15) :: uma = 1.66053878d-27 ! kg
 real(r15) :: cnst= 521.1098918 ! sqrt(ee*1d20/uma)/pi/200/c_light ! converts sqrt(eV/A/A/uma) to cm^-1
 real(r15) :: ryd= 27.2116  !me*ee**3/hbar**2/16/pi/pi/eps0/eps0 (in eV)
 real(r15) :: ab = 0.529177  ! hbar*hbar/me/ee/ee*4*pi*eps0
 real(r15) :: kbe=8.617343e-05  ! = value of 1 Kelvin in eV

 end module constants
!==============================================
module params
 real(8) tolerance,margin,scalelengths,alfaborn,tempk
 integer nconfigs,classical,ntemp,fdfiles,cal_cross,threemtrx,lamin,lamax,ncpu,n_dig_acc,itemp
 integer nshells(4,20)  ! up to which shell to include for each rank of FC (read from input file)
 integer include_fc(4),nsmax  ! whether to include FCs of that rank ,max# of shells looping
 real(8) rcut(4),tau0,wshift(3)
 real(8) tmin,tmax,qcros(3),svdc
 real(8) rcutoff_ewa,gcutoff_ewa,eta  ! these are for ewaldsums
 logical verbose

end module params
!============================================================
 module geometry

  type vector
    real(8) :: x,y,z
  end type vector
  type point
    real(8) :: x,y,z
  end type point
!-----------------------------
   interface operator(*)
     module procedure multiply_by_real,multiply_by_integer ; end interface
   interface operator(/)
     module procedure divide_by_scalar; end interface
   interface assignment(=)
     module procedure array_eq_vector , array_eq_point , vector_eq_array , &
&                     vector_eq_point , point_eq_vector, point_eq_array
   end interface
   interface operator(+)
     module procedure addition_avv , addition_vav, addition !,  &
! &                     addition_ava , addition_vaa
   end interface
   interface operator(-)
     module procedure subtraction,subtraction_av,subtraction_va; end interface
   interface operator(.dot.)
     module procedure dotproduct_v,dotproduct_a,  &
&                     dotproduct_av,dotproduct_va; end interface
   interface operator(.cross.)
     module procedure crossproduct_v,crossproduct_a; end interface
   interface operator(.myeqz.)
     module procedure myequal0, myequal0array, myequal0vector ; end interface
   interface operator(.myeq.)
     module procedure myequal, myequalarray, myequalvector ; end interface

  interface length
     module procedure lengthv,lengtha ; end interface

  interface bring_to_cell_c
     module procedure bring_to_cell_cv,bring_to_cell_ca
  end interface

  interface bring_to_cell_d
     module procedure bring_to_cell_dv,bring_to_cell_da
  end interface

  interface calculate_volume
     module procedure calvol_a, calvol_v
  end interface

  interface v2a
     module procedure v2a_r,v2a_a
  end interface

! interface cart_to_direct
!    module procedure cart_to_direct_aa,cart_to_direct_av,cart_to_direct_v
! end interface

! interface direct_to_cart
!    module procedure direct_to_cart_aa,direct_to_cart_av,direct_to_cart_v
! end interface

  interface reduce
     module procedure reduce_v,reduce_v2,reduce_a1,reduce_a4,reduce_a5
  end interface

  interface det
     module procedure idet,rdet
  end interface

  interface bring_to_center
     module procedure bring_to_center_a, bring_to_center_v
  end interface bring_to_center

   contains


subroutine calvol_a(r1,r2,r3,om)
implicit none
real(8) om
real(8) r1(3),r2(3),r3(3),cross12(3)

cross12 = r1 .cross. r2
om = abs(r3 .dot. cross12)

end subroutine calvol_a
!-----------------------------------
subroutine calvol_v(r1,r2,r3,om)
implicit none
real(8) om
type(vector) cross12,r1,r2,r3

cross12 = r1 .cross. r2
om = abs(r3 .dot. cross12)

end subroutine calvol_v
!------------------------------------------------------------------------------
      function rdet(mat)
      implicit none
!
!!	FIND THE DETERMINANT OF A 3X3 real MATRIX MAT
!
      real(8) rdet,mat(3,3)
      rdet=mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))  &
     &  - mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1))  &
     &  + mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))

      end function rdet
!------------------------------------------------------------------------------
      function idet(mat)
      implicit none
!
!! FIND THE DETERMINANT OF A 3X3 integer MATRIX MAT
!
      integer idet,mat(3,3)
      idet=mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))  &
     &   - mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1))  &
     &   + mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))

      end function idet
!------------------------------------------------------------------------------
   function crossproduct_v(v,w) result(cross)
     type(vector) cross !, intent(out) ::
     type(vector), intent(in) :: v,w
     cross%x = v%y*w%z-v%z*w%y
     cross%y = v%z*w%x-v%x*w%z
     cross%z = v%x*w%y-v%y*w%x
   end function crossproduct_v

   function crossproduct_a(v,w) result(cross)
     real(8) cross(3) !, intent(out) ::
     real(8), intent(in) :: v(3),w(3)
     cross(1) = v(2)*w(3)-v(3)*w(2)
     cross(2) = v(3)*w(1)-v(1)*w(3)
     cross(3) = v(1)*w(2)-v(2)*w(1)
   end function crossproduct_a
!-----------------------------------
   function dotproduct_v(v,w) result(dot)
     real(8) dot !, intent(out) ::
     type(vector), intent(in) :: v,w
     dot = w%x * v%x + w%y * v%y + w%z * v%z
   end function dotproduct_v
!-----------------------------------
   function dotproduct_a(v,w) result(dot)
     real(8) dot !, intent(out) ::
     real(8), intent(in) :: v(:),w(:)
!    dot = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
     dot = sum(v*w)
   end function dotproduct_a
!-----------------------------------
   function dotproduct_av(v,w) result(dot)
     real(8) dot !, intent(out) ::
     real(8), intent(in) :: v(3)
     type(vector), intent(in) :: w
     dot = v(1)*w%x + v(2)*w%y + v(3)*w%z
   end function dotproduct_av
!-----------------------------------
   function dotproduct_va(w,v) result(dot)
     real(8) dot !, intent(out) ::
     real(8), intent(in) :: v(3)
     type(vector), intent(in) :: w
     dot = v(1)*w%x + v(2)*w%y + v(3)*w%z
   end function dotproduct_va
!-----------------------------------
   function myequal0(v,w) result(eq)
     use params
     real(8), intent(in)::  v,w
     logical eq !, intent(out) :: eq

        if (abs(v-w) .lt. tolerance) then
           eq=.true.
        else
           eq=.false.
        endif

   end function myequal0
!-----------------------------------
   function myequal(v,w) result(eq)
     use params
     real(8), intent(in)::  v,w
     logical eq !, intent(out) :: eq
!    if (v.eq.0 ) then
!       if (w.eq.0) then
!          eq=.true.
!       else
!          eq=.false.
!       endif
!    else
!       if (abs(v-w)/abs(v) .lt. tolerance) then
        if (abs(v-w) .lt. tolerance) then
           eq=.true.
        else
           eq=.false.
        endif
!    endif
   end function myequal
!-----------------------------------
   function myequal0vector(v,w) result(eq)
     type(vector), intent(in)::  v,w
     logical eq !, intent(out) :: eq
  if ( (v%x .myeqz. w%x) .and. (v%y .myeqz. w%y) .and. (v%z .myeqz. w%z) ) then
           eq=.true.
     else
           eq=.false.
     endif
   end function myequal0vector
!-----------------------------------
   function myequalvector(v,w) result(eq)
     type(vector), intent(in)::  v,w
     logical eq !, intent(out) :: eq
     if ( (v%x .myeq. w%x) .and. (v%y .myeq. w%y) .and. (v%z .myeq. w%z) ) then
           eq=.true.
     else
           eq=.false.
     endif
   end function myequalvector
!-----------------------------------
   function myequal0array(v,w) result(eq)
     real(8), dimension(:), intent(in) ::  v,w
     logical eq !, intent(out) :: eq
     integer i,n
     i = size(v) ; n=size(w)
     if (n .ne. i) then
        print*, 'MYEQUAL0ARRAY: the input arrays are of different size ',i,n
        stop
     else
        eq = .true.
        loop: do n=1,i
           if ( .not. myequal0(v(n),w(n)) ) then
              eq = .false.
              exit loop
           endif
        enddo loop
     endif
   end function myequal0array
!-----------------------------------
   function myequalarray(v,w) result(eq)
     real(8), dimension(:), intent(in) ::  v,w
     logical eq !, intent(out) :: eq
     integer i,n
     i = size(v) ; n=size(w)
     if (n .ne. i) then
        print*, 'MYEQUALARRAY: the input arrays are of different size ',i,n
        stop
     else
        eq = .true.
        loop: do n=1,i
           if ( .not. myequal(v(n),w(n)) ) then
              eq = .false.
              exit loop
           endif
        enddo loop
     endif
   end function myequalarray
!-----------------------------------
   function addition(v,w) result(add)
     type(vector) add !, intent(out) ::
     type(vector), intent(in) :: v,w
     add%x = v%x + w%x
     add%y = v%y + w%y
     add%z = v%z + w%z
   end function addition
!-----------------------------------
   function addition_vav(v,w) result(add)
     type(vector) add !, intent(out) ::
     type(vector), intent(in) :: v
     real(8), intent(in) :: w(3)
     add%x = v%x + w(1)
     add%y = v%y + w(2)
     add%z = v%z + w(3)
   end function addition_vav
!-----------------------------------
   function addition_avv(w,v) result(add)
     type(vector) add !, intent(out) ::
     type(vector), intent(in) :: v
     real(8), intent(in) :: w(3)
     add%x = v%x + w(1)
     add%y = v%y + w(2)
     add%z = v%z + w(3)
   end function addition_avv
!-----------------------------------
   function subtraction(v,w) result(dif)
     type(vector) dif !, intent(out) ::
     type(vector), intent(in) :: v,w
     dif%x = v%x - w%x
     dif%y = v%y - w%y
     dif%z = v%z - w%z
   end function subtraction
!-----------------------------------
   function subtraction_va(v,w) result(dif)
     type(vector) dif !, intent(out) ::
     type(vector), intent(in) :: v
     real(8), intent(in) :: w(3)
     dif%x = v%x - w(1)
     dif%y = v%y - w(2)
     dif%z = v%z - w(3)
   end function subtraction_va
!-----------------------------------
   function subtraction_av(w,v) result(dif)
     type(vector) dif !, intent(out) ::
     type(vector), intent(in) :: v
     real(8), intent(in) :: w(3)
     dif%x = v%x - w(1)
     dif%y = v%y - w(2)
     dif%z = v%z - w(3)
   end function subtraction_av
!-----------------------------------
   function multiply_by_integer(s,v) result(sv)
   integer, intent(in) :: s
   type(vector), intent(in) :: v
   type(vector) sv !, intent(out)::
     sv%x=s*v%x
     sv%y=s*v%y
     sv%z=s*v%z
   end function multiply_by_integer
!-----------------------------------
   function multiply_by_real(s,v) result(sv)
   real(8), intent(in) :: s
   type(vector), intent(in) :: v
   type(vector) sv !, intent(out)::
     sv%x=s*v%x
     sv%y=s*v%y
     sv%z=s*v%z
   end function multiply_by_real
!-----------------------------------
   function divide_by_scalar(v,s) result(sv)
   real(8), intent(in) :: s
   type(vector), intent(in) :: v
   type(vector) sv !, intent(out)::
     sv%x=1/s*v%x
     sv%y=1/s*v%y
     sv%z=1/s*v%z
   end function divide_by_scalar
!-----------------------------------
   function distance(v,w) result(d)
   type(point), intent(in) :: v,w
   real(8) d !, intent(out)::
     d = sqrt( (v%x-w%x)**2 +(v%y-w%y)**2 +(v%z-w%z)**2 )
   end function distance
!-----------------------------------
    subroutine point_eq_array(p,a)
      real(8), intent(in) :: a(3)
      type(point), intent(out) :: p
      p%x = a(1)
      p%y = a(2)
      p%z = a(3)
    end subroutine point_eq_array
!-----------------------------------
    subroutine point_eq_vector(p,v)
      type(vector), intent(in) :: v
      type(point), intent(out) :: p
      p%x = v%x
      p%y = v%y
      p%z = v%z
    end subroutine point_eq_vector
!-----------------------------------
    subroutine vector_eq_point(v,p)
      type(point), intent(in) :: p
      type(vector), intent(out) :: v
      v%x = p%x
      v%y = p%y
      v%z = p%z
    end subroutine vector_eq_point
!-----------------------------------
    subroutine vector_eq_array(vv,aa)
      real(8), intent(in) :: aa(3)
      type(vector), intent(out) :: vv
      vv%x = aa(1)
      vv%y = aa(2)
      vv%z = aa(3)
    end subroutine vector_eq_array
!-----------------------------------
    subroutine array_eq_vector(a,v)
      type(vector), intent(in) :: v
      real(8), intent(out) :: a(3)
      a(1) = v%x
      a(2) = v%y
      a(3) = v%z
    end subroutine array_eq_vector
!-----------------------------------
    subroutine array_eq_point(a,p)
      type(point), intent(in) :: p
      real(8), intent(out):: a(3)
      a(1) = p%x
      a(2) = p%y
      a(3) = p%z
    end subroutine array_eq_point
!-----------------------------------
   function vect(a,b) result(v)
     type(point), intent(in) :: a,b
     type(vector) v
     v%x= b%x-a%x
     v%y= b%y-a%y
     v%z= b%z-a%z
   end function vect
!-----------------------------------
   function lengthv(v) result(l)
     real(8) l
     type(vector), intent(in) :: v
     l = sqrt(v%x*v%x+v%y*v%y+v%z*v%z)
   end function lengthv
!-----------------------------------
   function lengtha(v) result(l)
     real(8) l
     real(8), dimension(:), intent(in) :: v
     integer i
     l = 0d0
     do i=1,size(v)
        l = l + v(i)*v(i)
     enddo
     l = sqrt(l)
   end function lengtha
!-----------------------------------
 subroutine reduce_v(v,q1,q2,q3,w)
!! takes cartesian coordinates and returns direct coordinates
 implicit none
 type(vector) , intent(in) :: v,q1,q2,q3
 type(vector) , intent(out) :: w
 real(8) a1,a2,a3

 a1 = v .dot. q1
 a2 = v .dot. q2
 a3 = v .dot. q3
 w%x = a1 ; w%y = a2 ; w%z = a3

 end subroutine reduce_v
!-----------------------------------
 subroutine reduce_v1(v,q1,q2,q3,w)
!! takes cartesian coordinates and returns direct coordinates
 implicit none
 type(vector) , intent(in) :: v,q1,q2,q3
 real(8) , intent(out) :: w(3)
 real(8) a1,a2,a3

 w(1) = v .dot. q1
 w(2) = v .dot. q2
 w(3) = v .dot. q3

 end subroutine reduce_v1
!-----------------------------------
 subroutine reduce_v2(v,q1,q2,q3,w)
!! takes cartesian coordinates and returns direct coordinates
 implicit none
 type(vector) , intent(in) :: q1,q2,q3
 real(8) , intent(in) :: v(3)
 real(8) , intent(out) :: w(3)
 real(8) a1,a2,a3

 w(1) = v .dot. q1
 w(2) = v .dot. q2
 w(3) = v .dot. q3

 end subroutine reduce_v2
!-----------------------------------------
 subroutine reduce_a1(a,q1,q2,q3,w)
! takes cartesian coordinates and returns direct coordinates
 implicit none
 type(vector) , intent(in) :: q1,q2,q3
 type(vector) v,w
 real(8), intent(in) :: a(3)
 real(8) a1,a2,a3

 v%x = a(1);v%y = a(2);v%z = a(3)
 a1 = v .dot. q1
 a2 = v .dot. q2
 a3 = v .dot. q3
 w%x = a1 ; w%y = a2 ; w%z = a3

 end subroutine reduce_a1
!-----------------------------------------
 subroutine reduce_a4(a,q1,q2,q3,w)
! takes cartesian coordinates and returns direct coordinates
 implicit none
 real(8) , intent(in) :: a(3),q1(3),q2(3),q3(3)
 type(vector) w
 real(8) a1,a2,a3

 a1 = a .dot. q1
 a2 = a .dot. q2
 a3 = a .dot. q3
 w%x = a1 ; w%y = a2 ; w%z = a3

 end subroutine reduce_a4
!-----------------------------------------
 subroutine reduce_a5(a,q1,q2,q3,w)
! takes cartesian coordinates and returns direct coordinates
 implicit none
 real(8), intent(in) :: a(3),q1(3),q2(3),q3(3)
 real(8), intent(out) :: w(3)

 w(1) = a .dot. q1
 w(2) = a .dot. q2
 w(3) = a .dot. q3

 end subroutine reduce_a5
!-----------------------------------------
 function bring_to_cell_dv(v) result(w)
! takes direct coordinates and returns direct coordinates
 implicit none
 type(vector) v,w
 real(8) a1,a2,a3

 a1 = v%x ; a2 = v%y ; a3 = v%z
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
 w%x = a1 ; w%y = a2 ; w%z = a3

 end function bring_to_cell_dv
!-----------------------------------------
 function bring_to_cell_da(a) result(w)
! takes direct coordinates and returns direct coordinates
 implicit none
 type(vector) w
 real(8) a(3)

 a(1) = a(1) - floor(a(1))
 a(2) = a(2) - floor(a(2))
 a(3) = a(3) - floor(a(3))
 w%x = a(1) ; w%y = a(2) ; w%z = a(3)

 end function bring_to_cell_da
!-----------------------------------------
 subroutine bring_to_cell_ca(a,r1,r2,r3,g1,g2,g3,b)
! takes cart coordinates and returns cart coordinates within the supercell
! assumes r_i.g_j=delta_ij (no factor of 2pi included)
 implicit none
 real(8), intent(in) :: a(3)
 real(8), intent(out) :: b(3)
 type(vector), intent(in) :: g1,g2,g3,r1,r2,r3
 type(vector) v,w
 real(8) a1,a2,a3

 v%x = a(1);v%y = a(2);v%z = a(3)
! get direct coordinates first
 a1 = v .dot. g1
 a2 = v .dot. g2
 a3 = v .dot. g3
! bring into cell: between 0 and cell
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
! convert to cartesian coordinates
 w = a1*r1+a2*r2+a3*r3
 b(1)=w%x ; b(2)=w%y ; b(3)=w%z

 end subroutine bring_to_cell_ca
!-----------------------------------
 subroutine bring_to_cell_cv(v,r1,r2,r3,g1,g2,g3,w)
! takes cart coordinates and returns cart coordinates within the supercell
! assumes r_i.g_j=delta_ij (no factor of 2pi included)
 implicit none
 type(vector), intent(in) :: v,g1,g2,g3,r1,r2,r3
 type(vector), intent(out) :: w
 real(8) a1,a2,a3

! get direct coordinates first
 a1 = v .dot. g1
 a2 = v .dot. g2
 a3 = v .dot. g3
!print*,'v=',v
!print*,'g1,g2,g3=',g1,g2,g3
!print*,'a1,a2,a3=',a1,a2,a3
! bring into cell: between 0 and cell
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
!print*,'NEW a1,a2,a3=',a1,a2,a3
! convert to cartesian coordinates
 w = a1*r1+a2*r2+a3*r3

 end subroutine bring_to_cell_cv
!-----------------------------------------
 subroutine bring_to_center_v(v,r1,r2,r3,g1,g2,g3,w)
 use constants
! takes cart coordinates and returns cart coordinates within the FBZ
 implicit none
 type(vector), intent(in) :: v,g1,g2,g3,r1,r2,r3
 type(vector), intent(out) :: w
 real(8) a1,a2,a3

! get direct coordinates first
 a1 = v .dot. g1
 a2 = v .dot. g2
 a3 = v .dot. g3
 a1 = a1/2/pi ;  a2 = a2/2/pi ;  a3 = a3/2/pi
!print*,'v=',v
!print*,'g1,g2,g3=',g1,g2,g3
!print*,'a1,a2,a3=',a1,a2,a3
! bring into cell: between -cell/2 and cell/2
 a1 = a1 - floor(a1) - 0.5
 a2 = a2 - floor(a2) - 0.5
 a3 = a3 - floor(a3) - 0.5
!print*,'NEW a1,a2,a3=',a1,a2,a3
! convert to cartesian coordinates
 w = a1*r1+a2*r2+a3*r3

 end subroutine bring_to_center_v
!-----------------------------------------
 subroutine bring_to_center_a(v,r1,r2,r3,g1,g2,g3,w)
 use constants
! takes cart coordinates and returns cart coordinates within the FBZ
 implicit none
 type(vector), intent(in) :: g1,g2,g3,r1,r2,r3
 real(8), intent(in) :: v(3)
 real(8), intent(out) :: w(3)
 real(8) a1,a2,a3

! get direct coordinates first
 a1 = v .dot. g1
 a2 = v .dot. g2
 a3 = v .dot. g3
 a1 = a1/2/pi ;  a2 = a2/2/pi ;  a3 = a3/2/pi
!print*,'v=',v
!print*,'g1,g2,g3=',g1,g2,g3
!print*,'a1,a2,a3=',a1,a2,a3
! bring into cell: between -cell/2 and cell/2
 a1 = a1 - floor(a1) - 0.5
 a2 = a2 - floor(a2) - 0.5
 a3 = a3 - floor(a3) - 0.5
!print*,'NEW a1,a2,a3=',a1,a2,a3
! convert to cartesian coordinates
 w = a1*r1+a2*r2+a3*r3

 end subroutine bring_to_center_a
!-----------------------------------------
 function v2a_r(r0) result(rx)
 implicit none
 type(vector), intent(in):: r0
 real(8) rx(3) !, intent(out)::
 rx(1)=r0%x
 rx(2)=r0%y
 rx(3)=r0%z
 end function v2a_r
!-----------------------------------------
 function v2a_a(r0) result(rx)
 implicit none
 type(vector), dimension(:), intent(in):: r0
 real(8), dimension(3,size(r0)) :: rx !, intent(out)::
 integer i

 do i=1,size(r0)
    rx(1,i)=r0(i)%x
    rx(2,i)=r0(i)%y
    rx(3,i)=r0(i)%z
 enddo
! rx(2)=r0%y
! rx(3)=r0%z
 end function v2a_a
!-----------------------------------------
 function a2v(r0) result(rx)
 implicit none
 type(vector) rx
 real(8), intent(in):: r0(3)
 rx%x=r0(1)
 rx%y=r0(2)
 rx%z=r0(3)
 end function a2v

 end module geometry
!==============================================
module ios
 integer, parameter:: ulog=30,uposcar=10,uparams=11, utraj=40,ufco=20,  &
&         umap=60,umatrx=50,utimes=70,ufc1=21,ufc2=22,ufc3=23,ufc4=24,  &
&         ufc=80,ufit1=31,ufit2=32,ufit3=33,ufit4=34,ucor=93,uborn=12,  &
&         ujunk=79,uibz=80,uibs=81,ugr=82


  interface write_out
     module procedure write_outim, write_outrm, write_outiv, write_outrv, write_outv &
&    , write_outi, write_outr,write_outcv
  end interface

 contains

!-----------------------------
!  subroutine write_outrm(unit,string,n,m,var)
  subroutine write_outrm(unit,string,var)
  implicit none
  character*(*) string
  integer n,m,unit,l,i,j
  real(8), dimension(:,:), intent(in) :: var

  l=len_trim(string)
  n=size(var(:,1)) ; m=size(var(1,:))
  write(unit,*)' write_outrm called;nl,nc=',n,m
  write(unit,4)string(1:l)//' is='
  do i=1,n
     write(unit,5)(var(i,j),j=1,m)
  enddo
4 format(a,99(1x,g13.6))
5 format(199(1x,g12.5))
  end subroutine write_outrm
!-----------------------------
  subroutine write_outim(unit,string,n,m,var)
  implicit none
  character*(*) string
  integer n,m,unit,i,l,j
  integer var(n,m)

  l=len_trim(string)
  write(unit,4)string(1:l)//' is=',((var(i,j),i=1,n),j=1,m)
4 format(a,99(1x,i8))

  end subroutine write_outim
!-----------------------------
  subroutine write_outrv(unit,string,var)
  implicit none
  character*(*) string
  integer unit,l,i
  real(8), dimension(:) :: var

  l=len_trim(string)
  write(unit,4)string(1:l)//' is=',(var(i),i=1,size(var))
4 format(a,99(1x,g13.6))
  end subroutine write_outrv
!-----------------------------
  subroutine write_outcv(unit,string,var)
  implicit none
  character*(*) string
  integer unit,l,i
  complex(8), dimension(:) :: var

  l=len_trim(string)
  write(unit,4)string(1:l)//' is=',(var(i),i=1,size(var))
4 format(a,99(3x,g11.4,1x,g11.4))
  end subroutine write_outcv
!-----------------------------
  subroutine write_outiv(unit,string,var)
  implicit none
  character*(*) string
  integer unit,i,l
  integer, dimension(:) :: var

  l=len_trim(string)
  write(unit,4)string(1:l)//' is=',(var(i),i=1,size(var))
4 format(a,99(1x,i8))

end subroutine write_outiv
!-----------------------------
  subroutine write_outr(unit,string,var)
  implicit none
  character*(*) string
  integer unit,l
  real(8) var

  l=len_trim(string)
  write(unit,4)string(1:l)//' is=',var
4 format(a,99(1x,g13.6))
end subroutine write_outr
!-----------------------------
  subroutine write_outi(unit,string,var)
  implicit none
  character*(*) string
  integer n,unit,l
  integer var

  l=len_trim(string)
  write(unit,4)string(1:l)//' is=',var
4 format(a,99(1x,i8))

end subroutine write_outi
!-----------------------------
  subroutine write_outv(unit,string,var)
  use geometry
  implicit none
  character*(*) string
  integer unit,l
  type(vector) var

  l=len_trim(string)
  write(unit,4)string(1:l)//' is=',var
4 format(a,3(1x,g13.6))
end subroutine write_outv

end module ios
!===========================================================
 module atoms_force_constants
! the type atom_id0 concerns the properties of the primitive unit cell
! it s atoms, masses, neighbors, shells and force constants, whereas the
! type atomic refers to atoms in the supercell, their displacements and forces.
 use geometry
 implicit none
 integer natom_super_cell,nmax, fc2flag

!-------------------------
 type cell_id         ! id of atoms: position within cell, and cell coordinates
    integer n(3)
    integer tau
!    integer at_type
 end type
!-------------------------
 type shell
    integer no_of_neighbors  ! within that shell
    real(8) rij              ! radius of that shell
    type(cell_id) :: neighbors(296)  ! id of atoms in that shell
 end type
!-------------------------
 type atomic
!   type(tensor1) u       ! displacement from equilibrium
!   type(tensor1) force   ! force applied on it
    type(vector) equilibrium_pos
    type(cell_id) cell
    integer at_type
    real(8) mass
    real(8) charge(3,3)
 end type
!-------------------------
 type atom_id0   ! all characteristics of an atom in the primitive central cell
    integer at_type
    character name*2
    real(8) mass,charge(3,3)
    integer nshells       ! how many shells are included=actual dim(shell)
    type(vector) equilibrium_pos
    type(shell), allocatable ::  shells(:) 
!    type(tensor2r), pointer :: phi2(:)  ! refers to atom numbers neighbor of i
!    type(tensor3), pointer :: phi3(:,:)
!    type(tensor4), pointer :: phi4(:,:,:)
 end type
!-------------------------
! everything with 0 refers to the primitive cell

 integer  natom_type
 integer, allocatable:: natom(:),atom_type(:), map_rtau_sc(:,:)!label_of_primitive_atoms(:),
 real(8), allocatable:: atompos0(:,:),mas(:),force(:,:,:),displ(:,:,:),vel(:,:)
 real(8), allocatable:: cur(:,:)
 character(2), allocatable:: atname(:)
 type (atom_id0), allocatable :: atom0(:)
 type (atomic), allocatable :: atom_sc(:),atom_shl(:)

!---------------------------
! Module for saved data
!      module force_constants_module
! maximum shells of nearest neighbors (along radial direction) , and actual # of neighborshells
      integer maxneighbors,maxshell
!     parameter(maxneighbors=18 )
! maximum number of atoms out to maxneighbors
      integer maxatoms,imaxat
!     parameter(maxatoms=2800 )
! op_matrix(k,j,i), matrix for the ith point operator
      double precision op_matrix(3,3,48)
! op_kmatrix(k,j,i), matrix for the ith point operator acting on k vector
      double precision op_kmatrix(3,3,48)
      integer lattpgcount
! isgopcount, number of operators in space group
      integer isgopcount
! isgop(i), point operation in ith space group operator
      integer isgop(48)
! sgfract(j,i), jth cartesian coordinate of fractional in ith space group operator
      double precision sgfract(3,48)
! iatomop(j,i), point operator that takes jth atom into ith atom
      integer iatomop(:,:)
! atomopfract(k,j,i), kth cartesian coordinate of fractional to accompany iatomop
      double precision atomopfract(:,:,:)
! natom_prim_cell, number of atoms in the primitive unit cell
      integer natom_prim_cell
! natoms, number of atoms out to maxneighbors
      integer natoms
! atompos(j,i), jth cartesian coordinate of ith atom
!  double precision atompos(3,maxatoms)
      real(8), allocatable :: atompos(:,:)
! iatomcell(j,i), linear combination of basis vectors of the primitive lattice
! that takes us to the unit cell containing the ith atom
!  integer iatomcell(3,maxatoms)
      integer, allocatable ::  iatomcell(:,:),iatomcell0(:)
! iatomcell0(i), identity of atom in unit cell at origin equivalent to ith atom
!  integer iatomcell0(maxatoms)
! iatomneighbor(j,i), nearest neighbor shell of jth atom in primitive unit
! cell at the origin that contains the ith atom
      integer iatomneighbor(:,:)
      allocatable iatomneighbor,iatomop,atomopfract
!      end module force_constants_module
! -------------------------------------

contains

 subroutine allocate_disp_forc(n)  ! needed for the md code
    integer n
!   allocate( disp(3,n),forc(3,n),vel(3,n),cur(3,n),ncel1(n),ncel2(n),ncel3(n) )
!    allocate( disp(3,n),forc(3,n),vel(3,n),cur(3,n) )
 end subroutine allocate_disp_forc

 subroutine allocate_pos(n,m)
    integer n,m
    allocate( displ(3,n,m),force(3,n,m) )
 end subroutine allocate_pos

 subroutine allocate_mass(n)
    integer n
    allocate( mas(n) ,natom(n), atname(n) )
 end subroutine allocate_mass

 subroutine allocate_primcell(n0)
 integer n0
!print*,'n0=',n0
 allocate( atom0(n0),atompos0(3,n0), atom_type(n0) )  !label_of_primitive_atoms(n0),
 end subroutine allocate_primcell

 subroutine allocate_supercell(n)
 integer n
 allocate( atom_sc(n) )
 end subroutine allocate_supercell

!------------------------------------
!subroutine set_neighbor_list (maxshell)
!! Generates a set of neighbors to the primitive cell (up to maxshell, which should be large enough)
!! and finds which shell corresponds to the cutoff rcut set by the WS cell
!! this is needed in order to know what force constants to associate
!! to a pair ij, or ijk or ijkl
!! if a neighbor j is known then one can get the vect(i,j) between
!! their equilibrium positions
!! use atoms_force_constants , only : natom_prim_cell
!use params
!use lattice
!use ios
!use geometry
!implicit none
!! integer, parameter :: max=500
!integer, intent(in) :: maxshell
!real(8), intent(in) :: rcut
!integer i0,j,shel_count,counter,msort(maxatoms),l
!real(8) dist(maxatoms),d_old,d_min,rmax

!! allocate( atom0(1:natom_prim_cell)%shells(maxshells) )
!rmax = 0 ; dist = 1d10
!i0loop: do i0=1,natom_prim_cell
!   allocate( atom0(i0)%shells(0:maxshell) )
!   do j=1,natoms
!      call calculate_distance(i0,j,atompos,maxatoms,dist(j) )
!!      if ( iatomneighbor(i0,j) .eq. nshells(2) .and. dist(j) .gt. rmax ) then
!      if ( iatomneighbor(i0,j) .eq. nshells(2,i0) .and. dist(j) .gt. rmax ) then
!           rmax = dist(j)
!      endif
!   enddo
!   call sort(natoms,dist,msort,maxatoms)

!   shel_count = -1
!   d_old = 0
!   write(ulog,*)' ========================================================='
!   write(ulog,*)' neighborlist of atom i0 =',i0

!   jloop: do j=1,min(500,natoms)
!      d_min = dist(msort(j))  ! start with the shortest distance
!      if ( (d_min .myeq. d_old) .and. (shel_count.ne.-1) ) then   ! same shell
!         counter = counter + 1
!      else    ! new shell
!         shel_count = shel_count+1
!         counter = 1
!         if ( shel_count .gt. maxshell ) then
!            write(ulog,*)maxshell,' shells completed'
!            write(ulog,*)'shel_count=',shel_count,' exceeded it for j=',j
!            exit jloop
!         endif
!         d_old = d_min
!      endif
!      if ( counter .gt. 296) then
!         write(ulog,*) ' counter in neighbors exceeded 296 ', counter
!         stop
!      endif
!      atom0(i0)%shells(shel_count)%no_of_neighbors = counter
!      atom0(i0)%shells(shel_count)%rij = d_min
!      atom0(i0)%shells(shel_count)%neighbors(counter)%tau = iatomcell0(msort(j))
!      atom0(i0)%shells(shel_count)%neighbors(counter)%n   = iatomcell(:,msort(j))
!      if(verbose) then
!!       write(ulog,*)' ------------------------------------'
!         write(ulog,5)' count, j, of tau,n1,n2,n3,shel#,nb#,dij,rj =',j,msort(j),iatomcell0(msort(j)),   &
!                iatomcell(:,msort(j)),iatomneighbor(i0,msort(j)),counter,d_min,(atompos(l,msort(j)),l=1,3)
!      endif
!   enddo jloop

!   do shel_count = 0 , maxshell
!      write(ulog,*)'shell#, neigh#=',shel_count,atom0(i0)%shells(shel_count)%no_of_neighbors
!   enddo
!! also initialize atom0%equilibrium_pos

!! shift the position of the first atom to the  origin
!!do i=1,natom_prim_cell
!!   atompos(:,i) = atompos(:,i)-atompos(:,1)
!!    atom0(i0)%equilibrium_pos = atompos0(1,i0)*r01 + atompos0(2,i0)*r02  &
!!&                             + atompos0(3,i0)*r03
!   atom0(i0)%equilibrium_pos%x = atompos(1,i0)
!   atom0(i0)%equilibrium_pos%y = atompos(2,i0)
!   atom0(i0)%equilibrium_pos%z = atompos(3,i0)
!!enddo
!! they might not be inside the primitive cell

!enddo i0loop
! rcut(2) = rmax

! format(a,8(2x,i4))
! format(a,' (',i4,2(1x,i4),1x,')')
! format(a,1x,i4,' (',f8.4,2(1x,f8.4),1x,')')
! format(a,3(2x,i4),' (',i4,2(1x,i4),1x,')',2(1x,i3),1x,f10.7,2x,9(1x,f9.6))
! format(a,1x,f10.7,2x,9(1x,f9.6))

!end subroutine set_neighbor_list


 end module atoms_force_constants
!===========================================================
 module lattice
 use geometry
 use constants
 implicit none
 type(vector) gs1,gs2,gs3,rs1,rs2,rs3  ! rri=translation vectors of the supercell
 type(vector) r1conv,r2conv,r3conv,g1conv,g2conv,g3conv
 type(vector) r01,r02,r03,g01,g02,g03  ! tr vect of prim cell and its recip spce
 real(8) volume_r,volume_g,volume_r0,volume_g0
 real(8) lattice_parameter,latticeparameters(6),primitivelattice(3,3)
 real(8),dimension(3,3):: conv_to_cart,prim_to_conv,conv_to_prim,prim_to_cart,cart_to_prim
 real(8) box(3),boxg(3),density
 real(8) r0g(3,3)
! real(8), allocatable:: rws_weights(:),gws_weights(:)
 integer n1min,n2min,n3min,n1max,n2max,n3max !,nrgrid,nggrid
! integer nr1(3),nr2(3),nr3(3)
! real(8), allocatable:: rgrid(:,:),ggrid(:,:),xgrid(:,:)
 real(8) gws26(3,26),rws26(3,26)  ! superlattice shells defining the WS of SL

  interface make_rg
     module procedure make_rga,make_rgv
  end interface

 contains

!-----------------------------------------
 subroutine transform_input_structure
!! takes input variables latticeparameters(6),primitivelattice(3,3) to construct the primitive cell
 use ios
 implicit none
 real(8) latp(6),temp(3,3),tempi(3,3)
 integer ier

!k1
!      d2r = 4d0*datan(1d0)/180d0 ! 3.1415926535897932384626/180d0
!k1
      conv_to_cart=0
      prim_to_conv=primitivelattice
      call xmatinv(prim_to_conv,conv_to_prim,ier)
      if(ier.ne.0) then
         write(*,*)'prim_to_conv does not have an inverse!'
         stop
      endif
      latp=latticeparameters ! for short
!      latp(4:6) = latp(4:6)*d2r   ! convert angles from degree to radian
      conv_to_cart(1,1)=latp(1)
      conv_to_cart(1,2)=latp(2)*dcosd(latp(6))
      conv_to_cart(2,2)=latp(2)*dsind(latp(6))
      conv_to_cart(1,3)=latp(3)*dcosd(latp(5))
      conv_to_cart(2,3)=latp(3)*(dcosd(latp(4))   &
     &     -dcosd(latp(6))*dcosd(latp(5)))  /dsind(latp(6))
      conv_to_cart(3,3)=sqrt(latp(3)**2-conv_to_cart(1,3)**2   &
     &     -conv_to_cart(2,3)**2)
      write(*,*)'conv_to_cart has conventional vectors in its columns='
      call write_out(ulog,' conventional cell(1) =',conv_to_cart(:,1))
      call write_out(ulog,' conventional cell(2) =',conv_to_cart(:,2))
      call write_out(ulog,' conventional cell(3) =',conv_to_cart(:,3))
      r1conv=conv_to_cart(:,1)
      r2conv=conv_to_cart(:,2)
      r3conv=conv_to_cart(:,3)

      call make_reciprocal_lattice_2pi(r1conv,r2conv,r3conv,g1conv,g2conv,g3conv)
      call write_out(ulog,' conventional recip(1) =',g1conv)
      call write_out(ulog,' conventional recip(2) =',g2conv)
      call write_out(ulog,' conventional recip(3) =',g3conv)

! prim_to_cart(i,j) is the ith cartesian coordinate of the jth primitive translation vector
      prim_to_cart=matmul(conv_to_cart,prim_to_conv)

      r01 = prim_to_cart(:,1)
      r02 = prim_to_cart(:,2)
      r03 = prim_to_cart(:,3)
      call write_out(ulog,'primitive lattice vector #1= ',prim_to_cart(:,1))
      call write_out(ulog,'primitive lattice vector #2= ',prim_to_cart(:,2))
      call write_out(ulog,'primitive lattice vector #3= ',prim_to_cart(:,3))

      call make_reciprocal_lattice_2pi(r01,r02,r03,g01,g02,g03)
      call write_out(ulog,' g01= ' ,g01)
      call write_out(ulog,' g02= ' ,g02)
      call write_out(ulog,' g03= ' ,g03)

 end subroutine transform_input_structure
!-----------------------------------------
 subroutine make_unitcell
!! lattice parameters and atompos0 read from structure.params, it
!! constructs the primitive vectors r0i and reciprocal vectors g0i
!! initializes the object atom0%equilibrium_pos, mass, name, type and charge
!! finds the symmetry operations of the unit cell
!! calculates number and position of the atoms within rcut(2)=15 of the primitive cell
 use ios
 use params
 use geometry
 use atoms_force_constants
! use force_constants_module
! use svd_stuff
 use constants
 implicit none
 character line*90,name*2
 integer i,j,tt,ttyp,natoms2,n1
 real(8) mss ,x1,x2,x3
 type(vector) tau1,vv

! generate r0i and g0i rom latticeparameters, primitivelattice
  call transform_input_structure
  call calculate_volume(r01,r02,r03,volume_r0)
  call calculate_volume(g01,g02,g03,volume_g0)

  write(ulog,*)' PRIMITIVE CELL *************************'
  call write_out(ulog,' r01=',r01)
  call write_out(ulog,' r02=',r02)
  call write_out(ulog,' r03=',r03)
  call write_out(ulog,' Volume of primitive cell ',volume_r0)
! call write_out(ulog,' g01=',g01)
! call write_out(ulog,' g02=',g02)
! call write_out(ulog,' g03=',g03)
  call write_out(ulog,' Volume of reciprocal primitive cell ',volume_g0)


  box(1) = length(r01)
  box(2) = length(r02)
  box(3) = length(r03)
  boxg(1) = length(g01)
  boxg(2) = length(g02)
  boxg(3) = length(g03)
  write(ulog,3)' box lengths = ',box
  write(ulog,3)' boxg lengths = ',boxg
  write(ulog,3)' volume_r0,volume_g0 = ',volume_r0,volume_g0

  write(ulog,*)' Reduced Atomic positions, mass ====================== '
  do i=1,natom_prim_cell
     atom0(i)%equilibrium_pos=atompos0(:,i)
     write(ulog,8)atom0(i)%name,atom0(i)%at_type,atompos0(:,i),atom0(i)%mass
  enddo
  write(ulog,*)' Cartesian Atomic positions,  ====================== '
  do i=1,natom_prim_cell
     x1=atom0(i)%equilibrium_pos%x
     x2=atom0(i)%equilibrium_pos%y
     x3=atom0(i)%equilibrium_pos%z
     atom0(i)%equilibrium_pos%x = x1*r1conv%x + x2*r2conv%x + x3*r3conv%x
     atom0(i)%equilibrium_pos%y = x1*r1conv%y + x2*r2conv%y + x3*r3conv%y
     atom0(i)%equilibrium_pos%z = x1*r1conv%z + x2*r2conv%z + x3*r3conv%z
     write(ulog,3)atom0(i)%name,atom0(i)%equilibrium_pos
  enddo

! generate neighbors in atompos up to a disctance rcut(2)=15 Ang
  rcut(2)=15  ! this is enough for a supercell of length 30 Ang

! get largest # of atoms within rcut; and nsmax=rcut/shortest translation+2 
  call get_upper_bounds(r01,r02,r03,rcut(2),maxatoms,nsmax)


  write(ulog,*)'*****************************************************************'
  write(ulog,6)'UPPER_BOUNDS for rcut(2)=',rcut(2),' maxatoms,nsmax=',maxatoms,nsmax
  write(ulog,*)'*****************************************************************'

  allocate(atompos(3,maxatoms), iatomcell(3,maxatoms),iatomcell0(maxatoms))

! find symmetry operations; the above allocation is for this subroutine
! requires maxneighbors to generate atoms within these shells 
  call force_constants_init(natom_prim_cell,atom_type,atompos0)

! cat_to_prim had reciprocal lattice vectors on its lines
  call write_out(ulog,'reciprocal lattice vectors ',cart_to_prim)


  density=sum(atom0(:)%mass)/volume_r0*uma*1d27
  write(ulog,*)' density (uma/Ang^3)=',sum(atom0(:)%mass)/volume_r0
  write(ulog,*)' density (  g/ cm^3)=',density
  write(ulog,*)' COORDINATES read successfully '

3 format(a,9(2x,g13.6))
4 format(i6,3(2x,g12.5),4(3x,i4))
6 format(a,(2x,g12.5),a,4(3x,i4))
8 format(a2,i3,9(1x,g13.6))

 end subroutine make_unitcell
!-----------------------------------------
! subroutine get_components(q,n,i,j,k )  !,g1,g2,g3)
!! for a given q-vector, it finds its integer components assuming it was
!! created as: nk=0 do i=0,N(1)-1 ; do j=0,N(2)-1; do k=0,N(3)-1; nk=nk+1 ;q=i*G1/N1+j*G2/N2+k*G3/N3
!! as the ones created in make_kp_reg with zero shift
!! it works even if there's a shift less than 0.5
! implicit none
! real(8) q(3) !,g1(3),g2(3),g3(3),r1(3),r2(3),r3(3)
!! type(vector) g1,g2,g3,r1,r2,r3
! integer i,j,k,n(3)
!
!! call make_reciprocal_lattice(g1,g2,g3,r1,r2,r3)
! i = nint(1+ n(1)* (q.dot.r1)/2/pi )
! j = nint(1+ n(2)* (q.dot.r2)/2/pi )
! k = nint(1+ n(3)* (q.dot.r3)/2/pi )
!
! end subroutine get_components
!-----------------------------------------
 subroutine make_rga(x1,x2,x3,q1,q2,q3,n)   ! n(i,j)=xi.qj
! matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
! matmul(r0g,v) gives the 3 reduced coordinates of an atom in the supercell if its
! reduced coordinates in the primitive cell are given by v
 use geometry
 implicit none
 real(8), intent(in ) :: x1(3),x2(3),x3(3),q1(3),q2(3),q3(3)
 real(8), intent(out) :: n(3,3)

 n(1,1) = x1 .dot. q1
 n(1,2) = x1 .dot. q2
 n(1,3) = x1 .dot. q3
 n(2,1) = x2 .dot. q1
 n(2,2) = x2 .dot. q2
 n(2,3) = x2 .dot. q3
 n(3,1) = x3 .dot. q1
 n(3,2) = x3 .dot. q2
 n(3,3) = x3 .dot. q3

 end subroutine make_rga
!-----------------------------------------
 subroutine make_rgv(x1,x2,x3,q1,q2,q3,n)   ! n(i,j)=xi.qj
! matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
! matmul(r0g,v) gives the 3 reduced coordinates of an atom in the supercell if its
! reduced coordinates in the primitive cell are given by v
 use geometry
 implicit none
 type(vector), intent(in) :: x1,x2,x3,q1,q2,q3
 real(8), intent(out) :: n(3,3)

 n(1,1) = x1 .dot. q1
 n(1,2) = x1 .dot. q2
 n(1,3) = x1 .dot. q3
 n(2,1) = x2 .dot. q1
 n(2,2) = x2 .dot. q2
 n(2,3) = x2 .dot. q3
 n(3,1) = x3 .dot. q1
 n(3,2) = x3 .dot. q2
 n(3,3) = x3 .dot. q3

 end subroutine make_rgv
!-----------------------------------------
 subroutine set_neighbor_list (rcut2,mxshell)
!! Generates a set of neighbors to the primitive cell up to rcut2
!! and calculates the corresponding number of neighborshells
!! maxshell = actual max number of neighbor shells  
! this is needed in order to know what force constants to associate
! to a pair ij, or ijk or ijkl
! if a neighbor j is known then one can get the vect(i,j) between
! their equilibrium positions
 use atoms_force_constants !, only : natom_prim_cell, atom0,maxatoms
 use params
 use ios
 use geometry
 implicit none
 integer, intent(out) :: mxshell
 real(8), intent(in) :: rcut2
 integer i0,j,shel_count,counter,msort(maxatoms),l
 real(8) dist(maxatoms),d_old,d_min,rmax

 mxshell=maxneighbors ! do not want to go beyond 50 shells, hopefully is large enough for rcut2

! allocate( atom0(1:natom_prim_cell)%shells(mxshells) )
 rmax = 0  ! is the largest radius among the natoms neighbors   
 dist = 1d10
 i0loop: do i0=1,natom_prim_cell
    allocate( atom0(i0)%shells(0:mxshell) )
    do j=1,natoms
       call calculate_distance(i0,j,atompos,maxatoms,dist(j) )
! this dist(j) is supposed to be less than rcut(2)
!      if ( iatomneighbor(i0,j) .eq. nshells(2,i0) .and. dist(j) .gt. rmax ) then
!      if ( dist(j) .gt. rcut2 ) then
!           write(ulog,*)'SET_NEIGHBOR_LIST: dist(j)>rcut ',i0,j,dist(j),rcut2
!      endif
    enddo
!   if (rmax.lt.rcut2) then
!      write(ulog,*)'rmax=',rmax,' is less than input cutoff of ',rcut2
!      write(ulog,*)'need to increase mxshell from ',mxshell,' or decrease rcut=',rcut2
!   endif
    call sort(natoms,dist,msort,maxatoms)

    shel_count = -1
    d_old = 0
    write(ulog,*)' ========================================================='
    write(ulog,*)' neighborlist of atom i0 =',i0

    jloop: do j=1,natoms
       d_min = dist(msort(j))  ! start with the shortest distance
       if ( (d_min .myeq. d_old) .and. (shel_count.ne.-1) ) then   ! same shell
          counter = counter + 1  ! counts atoms within a given shell 
       else    ! new shell
          if (d_min.gt.rcut2) cycle
          shel_count = shel_count+1
          counter = 1
          if ( shel_count .gt. mxshell ) then
             write(ulog,*)mxshell,' shells completed'
             write(ulog,*)'shel_count=',shel_count,' exceeded it for j=',j
             write(ulog,*)'Need to increase maxneighbors from ',maxneighbors,' or decrease rcut from ',rcut2
             exit jloop
          endif
          d_old = d_min
       endif
!      if ( counter .gt. 296) then
!         write(ulog,*) ' counter in neighbors exceeded 296 ', counter
!         stop
!      endif
       atom0(i0)%shells(shel_count)%no_of_neighbors = counter
       atom0(i0)%shells(shel_count)%rij = d_min
       atom0(i0)%shells(shel_count)%neighbors(counter)%tau = iatomcell0(msort(j))
       atom0(i0)%shells(shel_count)%neighbors(counter)%n   = iatomcell(:,msort(j))
       if(verbose) then
!       write(ulog,*)' ------------------------------------'
          write(ulog,5)' count, j, of tau,n1,n2,n3,shel#,nb#,dij,rj =',j,msort(j),iatomcell0(msort(j)),   &
&                iatomcell(:,msort(j)),iatomneighbor(i0,msort(j)),counter,d_min,(atompos(l,msort(j)),l=1,3)
       endif
    enddo jloop

    mxshell=shel_count
    maxshell=mxshell  ! to save it in the module atoms_force_constants
    write(ulog,*)'SET_NEIGHBOR_LIST: largest#of shells=',mxshell,' for modified maxneighbors=',maxneighbors
    do shel_count = 0 , mxshell
       write(ulog,*)'shell#, neigh#=',shel_count,atom0(i0)%shells(shel_count)%no_of_neighbors
    enddo
! also initialize atom0%equilibrium_pos

! shift the position of the first atom to the  origin
!do i=1,natom_prim_cell
!   atompos(:,i) = atompos(:,i)-atompos(:,1)
!    atom0(i0)%equilibrium_pos = atompos0(1,i0)*r01 + atompos0(2,i0)*r02  &
!&                             + atompos0(3,i0)*r03
    atom0(i0)%equilibrium_pos%x = atompos(1,i0)
    atom0(i0)%equilibrium_pos%y = atompos(2,i0)
    atom0(i0)%equilibrium_pos%z = atompos(3,i0)
!enddo
! they might not be inside the primitive cell

 enddo i0loop

2 format(a,8(2x,i4))
3 format(a,' (',i4,2(1x,i4),1x,')')
4 format(a,1x,i4,' (',f8.4,2(1x,f8.4),1x,')')
5 format(a,3(2x,i4),' (',i4,2(1x,i4),1x,')',2(1x,i3),1x,f10.7,2x,9(1x,f9.6))
6 format(a,1x,f10.7,2x,9(1x,f9.6))

 end subroutine set_neighbor_list
!-------------------------------------

 end module lattice
!============================================================
  module svd_stuff
! stuff related to Harold's subroutines, mappings and the svd matrices

 type groupmatrix
    real(8), allocatable:: mat(:,:)
    integer, allocatable:: iat(:,:),ixyz(:,:),iatind(:,:),ixyzind(:,:)
 end type groupmatrix
 type fulldmatrix
    type(groupmatrix), allocatable :: gr(:)
    integer, allocatable:: nt(:),ntind(:)  !,igr(:)
    character(1), allocatable:: err(:)
    integer ngr,ntotind,ntot  ! number of groups, total number of independent terms and full terms
 end type fulldmatrix

 integer maxrank,itrans,irot,ihuang,enforce_inv, size_kept_fc2,ngrnk(4)
 integer, allocatable :: keep_grp2(:),nlines(:)
 integer nterms(4),maxterms(4),maxtermsindep(4),ngroups(4),maxtermzero(4),maxgroups(4)
! integer ndindp(4),ndfull(4)
 parameter(maxrank=4)
 integer, allocatable :: nterm(:),ntermsindep(:)
 integer, allocatable :: iatmtermindp(:,:,:),ixyztermindp(:,:,:)
 integer, allocatable :: iatmtrm(:,:,:),ixyztrm(:,:,:)
 real(8), allocatable :: mapmat(:,:,:)
 real(8) svdcut,radius(4)
 integer, allocatable:: iatomterm_1(:,:),ixyzterm_1(:,:),igroup_1(:),map_1(:)
 integer, allocatable:: iatomterm_2(:,:),ixyzterm_2(:,:),igroup_2(:),map_2(:)
 integer, allocatable:: iatomterm_3(:,:),ixyzterm_3(:,:),igroup_3(:),map_3(:)
 integer, allocatable:: iatomterm_4(:,:),ixyzterm_4(:,:),igroup_4(:),map_4(:)
 character(1), allocatable:: err_1(:),err_2(:),err_3(:),err_4(:)
 type(fulldmatrix) map(4)
 real(8), allocatable:: fcrnk(:,:)
 real(8), allocatable:: ampterm_1(:),fcs_1(:)
 real(8), allocatable:: ampterm_2(:),fcs_2(:),grun_fc(:)
 real(8), allocatable:: ampterm_3(:),fcs_3(:)
 real(8), allocatable:: ampterm_4(:),fcs_4(:)
 real(8), allocatable:: amat(:,:),bmat(:),sigma(:),fcs(:),ahom(:,:),kernelbasis(:,:)
 real(8), allocatable:: a11ia12(:,:),fc1(:)
 real(8), allocatable:: atransl(:,:),btransl(:),arot(:,:),brot(:),ahuang(:,:),bhuang(:), &
 &                      aforce(:,:),bforce(:),energies(:)
 real(8), allocatable:: fc_ind(:) ! This is for getting force constant by reading FC2.dat file
 integer inv_constraints,force_constraints,dim_al,dim_ac,n_indep,newdim_al,dim_hom   &
&        ,transl_constraints, rot_constraints, huang_constraints,nindepfc
contains

 subroutine set_maxterms
   maxterms(1)=100
   maxterms(2)=500
   maxterms(3)=1800
   maxterms(4)=2000
   maxtermzero(1)=500
   maxtermzero(2)=2000
   maxtermzero(3)=5000
   maxtermzero(4)=8000
   maxtermsindep(1)=10
   maxtermsindep(2)=70
   maxtermsindep(3)=150
   maxtermsindep(4)=300
   maxgroups(1)=10
   maxgroups(2)=52
   maxgroups(3)=150
   maxgroups(4)=300
 end subroutine set_maxterms

!--------------------------------------------
 subroutine remove_zeros_all(nl,nc)
! removes lines which are all zero from the amatrix
 use ios
 use params
! use svd_stuff  !for itrans,irot,ihuang
 implicit none
 integer, intent(in) :: nc
 integer, intent(inout) :: nl
! real(8), dimension(:,:), intent(inout) :: amatr
! real(8), dimension(:)  , intent(inout) :: bmatr
! real(8), intent(inout) :: amatr(nl,nc)
! real(8), intent(inout) :: bmatr(nl)
 integer i,j,nz,tra,rot,hua,tra_out,rot_out,hua_out
 real(8), allocatable:: at(:,:),bt(:),ar(:,:),br(:),ah(:,:),bh(:),aa(:,:),bb(:)
 real(8) small,suml
 integer, allocatable:: zz(:)

 small=1d-7 !10  displacements or cos(t) are never that small

 tra=0; rot=0; hua=0
 if (itrans .ne. 0) tra= transl_constraints
 if (irot   .ne. 0) rot= rot_constraints
 if (ihuang .ne. 0) hua= huang_constraints

 call remove_zeros(tra,nc,amat(1:tra,1:nc),bmat(1:tra),  &
 &  tra_out,at(1:tra_out,1:nc),bt(1:tra_out))
 call remove_zeros(rot,nc,amat(tra+1:tra+rot,1:nc),bmat(tra+1:tra+rot),  &
 &  rot_out,ar(1:rot_out,1:nc),br(1:rot_out))
 call remove_zeros(hua,nc,amat(tra+rot+1:tra+rot+hua,1:nc),bmat(tra+rot+1:tra+rot+hua), &
 &  hua_out,ah(1:hua_out,1:nc),bh(1:hua_out))

 nl=tra_out+rot_out+hua_out+force_constraints
!  nl=tra    +rot    +hua    +force_constraints
 allocate(aa(nl,nc),bb(nl))
 aa(1:tra_out,1:nc)=at
 bb(1:tra_out   )=bt
 aa(tra_out+1:tra_out+rot_out,1:nc)=ar
 bb(tra_out+1:tra_out+rot_out   )=br
 aa(tra_out+rot_out+1:tra_out+rot_out+hua_out,1:nc)=ah
 bb(tra_out+rot_out+1:tra_out+rot_out+hua_out   )=bh
 aa(tra_out+rot_out+hua_out+1:nl,1:nc)=aforce
 bb(tra_out+1:tra_out   )=bt

! call remove_zeros(rot,nc,amat(tra+1:tra+rot,1:nc),bmat(tra+1:tra+rot), &
! & rot_out,aa(tra_out+1:tra_out+rot_out,1:nc),bb(tra_out+1:tra_out+rot_out))
! call remove_zeros(hua,nc,amat(tra+rot+1:tra+rot+hua,1:nc),bmat(tra+tro+1:tra+rot+hua), &
! &  hua_out,aa(tra_out+rot_out+1:tra_out+rot_out+hua_out,1:nc), &
! & bb(tra_out+rot_out+1:tra_out+rot_out+hua_out))

 allocate(zz(nl))
 zz = 1
 nz = 0
 do i=1,nl
! get the norm1 of each line
!   suml = 0
!   do j=1,nc
!      if ( abs(amatr(i,j)) .ge. suml ) suml = abs(amatr(i,j))
!   enddo
    suml = maxval(abs(amat(i,:)))     ! sqrt(sum(amatr(i,:)*amatr(i,:))/nc)
    if (suml.lt.small) then
       if (abs(bmat(i)) .lt. small) then
          if(verbose) write(ulog,*)' line #',i,' of Amatrix and Bmatrix is zero; it will be removed!'
       else
          write(ulog,*)' Inconsistency in line i of Amatrix which is zero!',i,suml
          write(ulog,*)' It will be removed; Corresponding line in bmatr is=',bmat(i)
          write(ulog,*)' Need to increase your range of FCs'
       endif
! record the line index for later removal
       zz(i) = 0
       nz = nz + 1  ! nz counts the number of zero lines
    endif
 enddo
 write(ulog,*)' number of lines=0 in amat are=',nz
 allocate(aa(nl-nz,nc),bb(nl-nz))

 j = 0
 do i=1,nl
    if (zz(i).ne.0) then
       j = j+1
       aa(j,:) = amat(i,:)
       bb(j) = bmat(i)
    elseif(i.le.transl_constraints) then
       tra=tra - 1 !nsl_constraints=transl_constraints-1
    elseif(i.le.transl_constraints+rot_constraints) then
       rot=rot - 1 !_constraints=rot_constraints-1
    elseif(i.le.transl_constraints+rot_constraints+huang_constraints) then
       hua = hua-1 !ng_constraints=huang_constraints-1
    else
       force_constraints=force_constraints-1
    endif
 enddo
 if (j.ne.nl-nz) then
    write(ulog,*)'REMOVE ZEROS: Inconsistency! j.ne. nl.nz ',j,nl-nz
    stop
 endif

 inv_constraints = tra+rot+hua !transl_constraints+rot_constraints+huang_constraints

 write(ulog,4)'After remove_zeros: NEW tra,rot,hua,force_constaints=',tra,  &
 &       rot,hua,force_constraints

 nl=j
 deallocate(amat,bmat)
 allocate(amat(nl,nc),bmat(nl))
 amat = aa
 bmat = bb
 write(ulog,*)'REMOVE ZEROS: new number of lines in amat=',nl
 deallocate(aa,bb,zz)
 4 format(a,4(i9))

 end subroutine remove_zeros_all

!--------------------------------------------
 subroutine remove_zeros(nl,nc,matin,vectin,nlout,matout,vectout)
! removes lines which are all zero from the matrix mat and the vector v
 use ios , only : ulog
! use params
 implicit none
 integer, intent(out) :: nlout
 integer, intent(in) :: nl,nc
 real(8), dimension(nl,nc), intent(in) :: matin
 real(8), dimension(nl)   , intent(in) :: vectin
 real(8), dimension(nl,nc), intent(out) :: matout
 real(8), dimension(nl   ), intent(out) :: vectout
 integer i,j,nz
! real(8), allocatable:: aa(:,:),bb(:)
 real(8) small,suml
 integer, allocatable:: zz(:)

 small=1d-7 !10  displacements or cos(t) are never that small

 allocate(zz(nl))
 zz = 1
 nz = 0
 do i=1,nl
    suml = maxval(abs(matin(i,:)))     ! sqrt(sum(amatr(i,:)*amatr(i,:))/nc)
    if (suml.lt.small) then
       if (abs(vectin(i)) .lt. small) then
!         if(verbose) write(ulog,*)' line #',i,' of Amatrix and Bmatrix is zero; it will be removed!'
       else
          write(ulog,*)' Inconsistency in line i of Amatrix which is zero!',i,suml
          write(ulog,*)' It will be removed; Corresponding line in bmatr is=',vectin(i)
          write(ulog,*)' Need to increase your range of FCs, or include higher orders of FCs'
       endif
! record the line index for later removal
       zz(i) = 0
       nz = nz + 1  ! nz counts the number of zero lines
    endif
 enddo
 write(ulog,*)' number of lines=0 in matin are=',nz
 nlout=nl-nz
 write(ulog,4)'After removing zeros: NEW # of lines=',nlout
! if(allocated(matout)) deallocate (matout)
! if(allocated(vectout)) deallocate (vectout)
! allocate(matout(nlout,nc),vectout(nlout))
! allocate(aa(nl-nz,nc),bb(nl-nz))

 j = 0
 do i=1,nl
    if (zz(i).ne.0) then
       j = j+1
       matout(j,:) = matin(i,:)
       vectout(j ) = vectin(i)
    endif
 enddo

! deallocate(matin,vectin)
! allocate(matout(nlout,nc),vectout(nlout))
! matin = aa
! vectin= bb

 deallocate(zz)

 4 format(a,4(i9))

 end subroutine remove_zeros

  end module svd_stuff
!===========================================================
  module born
  real(8) epsil(3,3),epsinv(3,3)
!  real(8), allocatable:: zeu(:,:,:)
  real(8) rho
  integer born_flag
  end module born

!===========================================================
subroutine sort3(n1,n2,n3,mp)
implicit none
integer n1,n2,n3,mp(3)

if( n2.ge.n1) then
  if (n3.ge.n2) then
     mp(1)=n1
     mp(2)=n2
     mp(3)=n3
  elseif( n3.ge.n1) then
     mp(1)=n1
     mp(2)=n3
     mp(3)=n2
  else
     mp(1)=n3
     mp(2)=n1
     mp(3)=n2
  endif
elseif( n3.ge.n1) then
     mp(1)=n2
     mp(2)=n1
     mp(3)=n3
 elseif(n3.ge.n2) then
     mp(1)=n2
     mp(2)=n3
     mp(3)=n1
 else
     mp(1)=n3
     mp(2)=n2
     mp(3)=n1
 endif
end subroutine sort3
!===========================================================
 subroutine check_inside_fbz_old(q,g1,g2,g3,inside)
 use geometry
 use constants
 implicit none
 real(8), intent(in):: q(3)
 type(vector), intent(in):: g1,g2,g3
 integer, intent(out) :: inside
 real(8) qdotg,gg,gt(3),sm
 sm = -.000001  ! on the boundary is also counted as inside, but once only

 inside = 0
!------------------- along the diagonal 100
       qdotg = ( q .dot. g1) + sm
       gg = g1 .dot. g1
       if(abs(qdotg) .gt. gg/2) return  ! this is the Bragg condition |q.G| < G.G/2
       qdotg = ( q .dot. g2) + sm
       gg = g2 .dot. g2
       if(abs(qdotg) .gt. gg/2) return
       qdotg = ( q .dot. g3) + sm
       gg = g3 .dot. g3
       if(abs(qdotg) .gt. gg/2) return
!------------------- along the diagonal 110
       gt = g1+g2
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1-g2
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1+g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1-g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return

       gt = g3+g2
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g3-g2
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
!------------------- along the diagonal 111
       gt = g1+g2+g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1+g2-g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1-g2+g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1-g2-g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
!-------------------
 inside = 1

 end subroutine check_inside_fbz_old


! end module kpoints
!===========================================================

