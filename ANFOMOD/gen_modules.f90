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

 real(8) tolerance,margin,scalelengths,alfaborn !,temperature_k
 integer nconfigs,classical,ntemp,fdfiles,cal_cross,threemtrx,lamin,lamax,ncpu,n_dig_acc
 integer nshells(4,20)  ! up to which shell to include for each rank of FC (read from input file)
 integer include_fc(4),nsmax  ! whether to include FCs of that rank ,max# of shells looping
 real(8) rcut(4),tau0,wshift(3)
 real(8) tmin,tmax,qcros(3),svdc,lmax ! lmax is the cutoff length of FC2 limited by the supercell WS
 real(8) rcutoff_ewa,gcutoff_ewa,eta! these are for ewaldsums
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
     module procedure lengthv,lengtha,lengthi ; end interface

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

  interface cart_to_direct
     module procedure cart_to_direct_aa,cart_to_direct_av,cart_to_direct_v
  end interface

  interface direct_to_cart
     module procedure direct_to_cart_aa,direct_to_cart_av,direct_to_cart_v
  end interface

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
   function lengthi(v) result(l)
     real(8) l
     integer, intent(in) :: v(:)
     l = sqrt(real(dot_product(v,v)))
   end function lengthi
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
!-----------------------------------
 subroutine cart_to_direct_v(v,g1,g2,g3,w)
! takes cart coordinates and returns cart coordinates within the supercell
 use constants, only : pi
 implicit none
 type(vector) v,w,g1,g2,g3

 w%x = (v .dot. g1)/(2*pi)
 w%y = (v .dot. g2)/(2*pi)
 w%z = (v .dot. g3)/(2*pi)

 end subroutine cart_to_direct_v
!-----------------------------------
 subroutine cart_to_direct_av(a,g1,g2,g3,w)
! takes cart coordinates and returns cart coordinates within the supercell
 implicit none
 real(8) a(3)
 type(vector) v,w,g1,g2,g3

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call cart_to_direct_v(v,g3,g3,g3,w)

 end subroutine cart_to_direct_av
!-----------------------------------
 subroutine cart_to_direct_aa(a,g1,g2,g3,b)
! takes cart coordinates and returns cart coordinates within the supercell
 implicit none
 real(8) a(3),b(3)
 type(vector) v,w,g1,g2,g3

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call cart_to_direct_v(v,g1,g2,g3,w)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end subroutine cart_to_direct_aa
!-----------------------------------
 subroutine direct_to_cart_v(v,r1,r2,r3,w)
! takes direct coordinates and returns cart coordinates
 implicit none
 type(vector) v,w,r1,r2,r3

 w%x = v%x*r1%x + v%y*r2%x + v%z*r3%x
 w%y = v%x*r1%y + v%y*r2%y + v%z*r3%y
 w%z = v%x*r1%z + v%y*r2%z + v%z*r3%z

 end subroutine direct_to_cart_v
!-----------------------------------
 subroutine direct_to_cart_av(a,r1,r2,r3,w)
! takes direct coordinates and returns cart coordinates
 implicit none
 real(8) a(3)
 type(vector) v,w,r1,r2,r3

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call direct_to_cart_v(v,r1,r2,r3,w)

 end subroutine direct_to_cart_av
!-----------------------------------
 subroutine direct_to_cart_aa(a,r1,r2,r3,b)
! takes direct coordinates and returns cart coordinates
 implicit none
 real(8) a(3),b(3)
 type(vector) v,r1,r2,r3,w

 b=v2a(a(1)*r1+a(2)*r2+a(3)*r3)
! v%x = a(1) ; v%y = a(2) ; v%z = a(3)
! call direct_to_cart_v(v,r1,r2,r3,w)
! b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end subroutine direct_to_cart_aa
!------------------------------------
 subroutine dir2cart_g(q,g1,g2,g3,k)
! takes q in direct coordinates and outputs k in cartesian
 real(8) q(3),k(3)
 type(vector)  g1,g2,g3

 k(:) = q(1)*g1 + q(2)*g2 + q(3)*g3

 end subroutine dir2cart_g
!------------------------------------
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

 integer, parameter:: uposcar=10,uparams=11,uborn=12, &
&         ufco=20,ulog=30,utraj=40,umatrx=50,umap=60,utimes=70,ufc=80, &
&         ufc1=21,ufc2=22,ufc3=23,ufc4=24,  &
&         ufit1=31,ufit2=32,ufit3=33,ufit4=34,  &
&         ujunk=79,uibz=80,uibs=81,ugrun=82,uband=83,ucor=93

  interface write_out
     module procedure write_outim, write_outrm, write_outiv, write_outrv, &
&    write_outv , write_outi, write_outr,write_outcv,write_outrm3
  end interface

 contains

!-----------------------------
  subroutine write_outrm3(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer n1,n2,n3,unt,l,j
  real(8), dimension(:,:,:), intent(in) :: var

  l=len_trim(string)
!  n1=size(var(:,1,1)) ; n2=size(var(1,:,1)) ; n3=size(var(1,1,:))
  n1=size(var,1) ; n2=size(var,2) ; n3=size(var,3)
! write(unt,*)' write_outrm called;nl,n2,n3=',n1,n2,n3
  write(unt,4)string(1:l)//' is='
!  do j=1,n2
!     write(unt,5)var(:,j,1)
!  enddo
  do j=1,n3
     call write_out(unt,' block #j ',j)
     call write_out(unt,' ',var(:,:,j))
  enddo
4 format(a,99(1x,g11.4))
5 format(199(1x,g11.4))
  end subroutine write_outrm3
!-----------------------------
!  subroutine write_outrm(unit,string,n,m,var)
  subroutine write_outrm(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer n,m,unt,l,i,j
  real(8), dimension(:,:), intent(in) :: var

  l=len_trim(string)
  n=size(var,1) ; m=size(var,2)
! write(unit,*)' write_outrm called;nl,nc=',n,m
  write(unt,4)string(1:l)//' is='
  do i=1,n
     write(unt,5)var(i,:)
  enddo
4 format(a,99(1x,g13.6))
5 format(199(1x,f12.5))
  end subroutine write_outrm
!-----------------------------
  subroutine write_outim(unt,string,var)
  implicit none
  integer, dimension(:,:), intent(in) :: var
  integer, intent(in) :: unt
  character*(*), intent(in) :: string
  integer n,m,i,l,j

  l=len_trim(string)
  n=size(var,1) ; m=size(var,2)
  do i=1,n
     write(unt,4)string(1:l)//' is=',var(i,:)
  enddo
4 format(a,99(1x,i8))

  end subroutine write_outim
!-----------------------------
  subroutine write_outrv(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,l,i
  real(8), dimension(:) :: var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var(:)
4 format(a,99(1x,g13.6))
  end subroutine write_outrv
!-----------------------------
  subroutine write_outcv(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,l,i
  complex(8), dimension(:) :: var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var(:)
4 format(a,99(3x,g11.4,1x,g11.4))
  end subroutine write_outcv
!-----------------------------
  subroutine write_outiv(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,i,l
  integer, dimension(:) :: var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var(:)
4 format(a,99(1x,i8))

end subroutine write_outiv
!-----------------------------
  subroutine write_outr(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,l
  real(8) var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var
4 format(a,99(1x,g13.6))
end subroutine write_outr
!-----------------------------
  subroutine write_outi(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,l
  integer var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var
4 format(a,99(1x,i8))

end subroutine write_outi
!-----------------------------
  subroutine write_outv(unt,string,var)
  use geometry
  implicit none
  character*(*), intent(in) :: string
  integer unt,l
  type(vector) var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var
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
    integer atompos  ! corresponding index in the atompos list
 end type
!-------------------------
 type shell
    integer no_of_neighbors  ! within that shell
    real(8) radius              ! radius of that shell
    type(cell_id) :: neighbors(296)  ! id of atoms in that shell ! keep up to 296 neighbors
 end type
!-------------------------
 type atomic
!   type(tensor1) u       ! displacement from equilibrium
!   type(tensor1) force   ! force applied on it
    type(vector) equilibrium_pos
    type(cell_id) cell
    integer at_type
    character name*2
    real(8) mass,charge(3,3)
 end type
!-------------------------
 type atom_id0   ! all characteristics of an atom in the primitive central cell
    integer at_type,tau !  (n1,n2,n3)=0 for atoms in the primitive cell
    character name*2
    real(8) mass,charge(3,3)
    integer nshells       ! how many shells are included for that atom = actual dim(shell)
    type(vector) equilibrium_pos
    type(shell), allocatable ::  shells(:)  ! of dim atom0(i)%nshells
!    type(tensor2r), pointer :: phi2(:)  ! refers to atom numbers neighbor of i
!    type(tensor3), pointer :: phi3(:,:)
!    type(tensor4), pointer :: phi4(:,:,:)
 end type
!-------------------------
! everything with 0 refers to the primitive cell

 integer  natom_type
 integer, allocatable:: natom(:),atom_type(:), map_rtau_sc(:,:)!label_of_primitive_atoms(:),
 real(8), allocatable:: atompos0(:,:),mas(:),force(:,:,:),displ(:,:,:),vel(:,:),energy(:)
! real(8), allocatable:: cur(:,:)
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
!     integer iatomcell0(maxatoms)
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

 subroutine allocate_edf(nat,ncfg)
    integer nat,ncfg
    allocate( displ(3,nat,ncfg),force(3,nat,ncfg),energy(ncfg) )
 end subroutine allocate_edf

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
 real(8) lattice_parameter,latticeparameters(6),primitivelattice(3,3),gconv_to_cart(3,3),gprim_to_cart(3,3)
 real(8),dimension(3,3):: conv_to_cart,prim_to_conv,conv_to_prim,prim_to_cart,cart_to_prim
 real(8) box(3),boxg(3),density
 real(8) r0g(3,3)
! real(8), allocatable:: rws_weights(:),gws_weights(:)
! integer nr1(3),nr2(3),nr3(3)
! real(8), allocatable:: rgrid(:,:),ggrid(:,:),xgrid(:,:)
 real(8) gws26(3,26),rws26(3,26),invn_sc(3,3)  ! superlattice shells defining the WS of SL
 integer n_sc(3,3)

  interface make_rg
     module procedure make_rga,make_rgv
  end interface

  interface check_int
     module procedure check_d,check_3,check_a
  end interface

 contains

!-------------------------------------------
 subroutine make_r0g
! matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
! matmul(r0g,v) gives the 3 reduced coordinates of an atom in the supercell if its
! coordinates are v(1)*r01+v(2)*r02+v(3)*r03 i.e. if v is the reduced coords in primitive units
 use constants , only : pi
 implicit none

 r0g(1,1) = r01 .dot. gs1
 r0g(2,1) = r01 .dot. gs2
 r0g(3,1) = r01 .dot. gs3
 r0g(1,2) = r02 .dot. gs1
 r0g(2,2) = r02 .dot. gs2
 r0g(3,2) = r02 .dot. gs3
 r0g(1,3) = r03 .dot. gs1
 r0g(2,3) = r03 .dot. gs2
 r0g(3,3) = r03 .dot. gs3
 r0g=r0g/(2*pi)

 end subroutine make_r0g
!-------------------------------------------
 subroutine check_a(r,a,ier,g1,g2,g3)
! subroutine to check whether r is an integer multiple of (r1,r2,r3)
! whose reciprocal basis is defined by (g1,g2,g3)
! output ai are the coefficients of its linear combination on this basis
! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use constants, only : pi
 use ios
 implicit none
 type(vector) :: r,g1,g2,g3
 real(8) a(3),ep
 integer ier,n(3)
 logical is_integer

 ep = tolerance ! displacements larger than 0.001 A are not tolerated
!call write_out(ulog,'CHECK: R ',r)
 a(1) = (r .dot. g1)/(2*pi)  ! g's have a 2 pi
 a(2) = (r .dot. g2)/(2*pi)
 a(3) = (r .dot. g3)/(2*pi)
 n=floor(a)
 if(.not.is_integer(a)) then
   ier = 1
 else
   ier = 0
 endif

 end subroutine check_a
 !============================================================
 subroutine check_3(r,a1,a2,a3,ier,g1,g2,g3)
!! subroutine to check whether r is an integer multiple of (r01,r02,r03)
!! output ai are the coefficients of its linear combination on this basis
!! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use constants, only : pi
 use ios
 implicit none

 type(vector) :: r,g1,g2,g3
 real(8) a1,a2,a3,ep
 integer ier

 ep = tolerance ! 0.001 ! displacements larger than 0.001 A are not tolerated
!call write_out(ulog,'CHECK: R ',r)
 a1 = (r .dot. g1)/(2*pi)  ! g's have a 2 pi
 a2 = (r .dot. g2)/(2*pi)
 a3 = (r .dot. g3)/(2*pi)
 if(abs(a1-nint(a1)).ge.ep .or.abs(a2-nint(a2)).ge.ep   &
&   .or.abs(a3-nint(a3)).ge.ep ) then
   ier = 1
!  write(ulog,*) ' R is not a multiple of r0s , check your inputs '
!  write(ulog,3) 'ier=1;r and r.gi/2pi=',r,a1,a2,a3
!  stop
 else
!  write(ulog,*) ' R is a multiple of r0s'
!  write(ulog,3) ' n1,n2,n3  =',a1,a2,a3
   ier = 0
 endif
3 format(a,9(1x,g11.4))

 end subroutine check_3
!-------------------------------------------
 subroutine check_d(r,a1,a2,a3,ier)
! subroutine to check whether r is an integer vector
! output ai are the coefficients of its linear combination on this basis
! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use ios
 implicit none
 type(vector) :: r
 real(8) a1,a2,a3,ep
 integer ier

 ep = tolerance ! 1d-3
!call write_out(ulog,'CHECK: R ',r)
 a1 = r%x
 a2 = r%y
 a3 = r%z
 if(abs(a1-nint(a1)).ge.ep .or.abs(a2-nint(a2)).ge.ep   &
&   .or.abs(a3-nint(a3)).ge.ep ) then
   ier = 1
   write(ulog,3) 'ier=1; r is not integer ',a1,a2,a3
 else
   ier = 0
 endif
3 format(a,9(1x,g11.4))

 end subroutine check_d
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
    type(groupmatrix), allocatable :: gr(:)  ! groups of fcs
    integer, allocatable:: nt(:),ntind(:)  ! full and independent fc terms for each group
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
 real(8) svdcut !,rayon(4)
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
   maxterms(1)=10 !100
   maxterms(2)=4000 !500
   maxterms(3)=5000 !1800
   maxterms(4)=3000 !2000
   maxtermzero(1)=20 !500
   maxtermzero(2)=1000 !2000
   maxtermzero(3)=5000!5000
   maxtermzero(4)=3000 !8000
   maxtermsindep(1)=10
   maxtermsindep(2)=150
   maxtermsindep(3)=150
   maxtermsindep(4)=300
   maxgroups(1)=10
   maxgroups(2)=100
   maxgroups(3)=150
   maxgroups(4)=300
 end subroutine set_maxterms


  end module svd_stuff

!===========================================================

  module born

  real(8) epsil(3,3),epsinv(3,3)
!  real(8), allocatable:: zeu(:,:,:)
! real(8) rho
  integer born_flag

  end module born

!===========================================================

 module linalgb

  interface append_array
     module procedure append_array_1d,append_array_2d
  end interface

 contains 

 subroutine append_array_1d(a,b,c)
!! appends array b to the end of a and stores the result in c
 implicit none
! integer, intent(in) :: na,nb
! real(8), intent(in):: a(na),b(nb)
! real(8), allocatable, intent(inout):: c(:)
! real(8) :: c(size(a)+size(b))
 real(8), intent(in):: a(:),b(:)
 real(8), allocatable :: c(:),aux(:)
 integer nc,na,nb,col

 na=size(a(:));
 nb=size(b(:));
 nc=na+nb
 c=reshape(a,shape=(/nc/),pad=b)
 return

! could also use
 allocate(aux(nc)) ! this is to allow calls like append_array(a,b,a)
 aux(1:na)=a
 aux(na+1:na+nb)=b
 if (allocated(c)) deallocate (c)
 allocate(c(nc))
 c=aux
 deallocate(aux)

 end subroutine append_array_1d
!---------------------------------------
 subroutine append_array_2d(a,b,c)
!! appends array b to the end of a and stores the result in c
 implicit none
! integer, intent(in) :: na,nb
! real(8), intent(in):: a(na),b(nb)
! real(8), allocatable, intent(inout):: c(:)
! real(8) :: c(size(a)+size(b))
 real(8), intent(in):: a(:,:),b(:,:)
 real(8), allocatable :: c(:,:) ,aux(:,:)
 integer nc,na,nb,col

 col=size(a,2) ! (1,:));
 na =size(a,1)   ! (:,1));
 nb =size(b,1)   ! (:,1));
 nc =na+nb

!  c=reshape(transpose(a),shape=(/col,nc/),pad=transpose(b),order=(/2,1/))
  c=reshape(transpose(a),shape=(/col,nc/),pad=transpose(b),order=(/1,2/))
  c=transpose(c)
! c=reshape(a,shape=(/nc,col/),pad=b,order=(/1,2/))
  return

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

 end module linalgb
