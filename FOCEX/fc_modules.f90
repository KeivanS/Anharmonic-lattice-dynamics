!==============================================================
 module constants
 implicit none
!--! specific precisions, usually same as real and double precision
! integer, parameter :: r6 = selected_real_kind(6)
integer, parameter,public :: r15 = kind(1.0d0)
! integer, parameter,public :: r15 = selected_real_kind(15) ! should work on any compiler
integer, parameter,public :: c15 = 8
!integer, parameter :: c6 = selected_complex_kind(6)
! integer, parameter :: dp=8 !(for nag compilers use kd=2 for double precision; 8 for gfortran)
! complex(kind=8) :: ci=dcmplx(0d0,1d0)
 complex(r15) :: ci=cmplx(0d0,1d0)
 real(r15) :: pi=3.14159265358979323846d0
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
 use constants, only : r15
 real(r15) tolerance,margin,scalelengths,alfaborn,tempk,rcutoff
 integer nconfigs,classical,ntemp,fdfiles,cal_cross,threemtrx,lamin,lamax,ncpu,n_dig_acc,itemp
 integer nshells(4,20)  ! up to which shell to include for each rank of FC (read from input file)
 integer include_fc(4) !,nsmax  ! whether to include FCs of that rank ,max# of shells looping
! real(r15) rcut(4),tau0,wshift(3)
 real(r15) tmin,tmax,qcros(3),svdc,lmax ! lmax is the cutoff length of FC2 limited by the supercell WS
 real(r15) rcutoff_ewa,gcutoff_ewa,eta! these are for ewaldsums
 logical verbose

end module params
!============================================================
 module geometry
 use constants, only : r15

  type vector
    real(r15) :: x,y,z
  end type vector
  type point
    real(r15) :: x,y,z
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
&                     dotproduct_av,dotproduct_va
   end interface

   interface operator(.cross.)
     module procedure crossproduct_v,crossproduct_a; end interface
   interface operator(.myeqz.)
     module procedure myequal0, myequal0array, myequal0vector ; end interface

   interface operator(.myeq.)
     module procedure myequal, myequal_i,myequalarray_i, &
&         myequalarray_r, myequalvector ; end interface

  interface length
     module procedure lengthv,lengtha,lengthi ; end interface

  interface calculate_volume
     module procedure calvol_a, calvol_v
  end interface

  interface v2a
     module procedure v2a_r,v2a_a
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

  interface is_integer
     module procedure is_integer_1,is_integer_a
  end interface is_integer

   contains

 function is_integer_1(i) result(is)
!! checks if the variable is integer
 use params, only : tolerance
 implicit none
 real(r15) i
 logical is

 if (abs(i-nint(i)).lt.tolerance ) then
    is=.true.
 else
    is=.false.
 endif
 end function is_integer_1
!-----------------------------------------
 function is_integer_a(i) result(is)
!! checks if the variable is integer
 use params, only : tolerance
 implicit none
 real(r15) i(:)
 logical is
 integer j,n

 n=size(i)
 do j=1,n
    if(abs(i(j)-nint(i(j))).gt.tolerance ) then
       is=.false.
       return
    endif
 enddo
 is=.true.

 end function is_integer_a
!-----------------------------------------

subroutine calvol_a(r1,r2,r3,om)
implicit none
real(r15) om
real(r15) r1(3),r2(3),r3(3),cross12(3)

cross12 = r1 .cross. r2
om = abs(r3 .dot. cross12)

end subroutine calvol_a
!-----------------------------------
subroutine calvol_v(r1,r2,r3,om)
implicit none
real(r15) om
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
      real(r15) rdet,mat(3,3)
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
     real(r15) cross(3) !, intent(out) ::
     real(r15), intent(in) :: v(3),w(3)
     cross(1) = v(2)*w(3)-v(3)*w(2)
     cross(2) = v(3)*w(1)-v(1)*w(3)
     cross(3) = v(1)*w(2)-v(2)*w(1)
   end function crossproduct_a
!-----------------------------------
   function dotproduct_v(v,w) result(dot)
     real(r15) dot !, intent(out) ::
     type(vector), intent(in) :: v,w
     dot = w%x * v%x + w%y * v%y + w%z * v%z
   end function dotproduct_v
!-----------------------------------
   function dotproduct_a(v,w) result(dot)
     real(r15) dot !, intent(out) ::
     real(r15), intent(in) :: v(:),w(:)
!    dot = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
     dot = sum(v*w)
   end function dotproduct_a
!-----------------------------------
   function dotproduct_av(v,w) result(dot)
     real(r15) dot !, intent(out) ::
     real(r15), intent(in) :: v(3)
     type(vector), intent(in) :: w
     dot = v(1)*w%x + v(2)*w%y + v(3)*w%z
   end function dotproduct_av
!-----------------------------------
   function dotproduct_va(w,v) result(dot)
     real(r15) dot !, intent(out) ::
     real(r15), intent(in) :: v(3)
     type(vector), intent(in) :: w
     dot = v(1)*w%x + v(2)*w%y + v(3)*w%z
   end function dotproduct_va
!-----------------------------------
   function myequal0(v,w) result(eq)
     use params
     real(r15), intent(in)::  v,w
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
     real(r15), intent(in)::  v,w
     logical eq

        if (abs(v-w) .lt. tolerance) then
           eq=.true.
        else
           eq=.false.
        endif

   end function myequal
!-----------------------------------
   function myequal_i(v,w) result(eq)
     use params
     integer, intent(in)::  v,w
     logical eq

        if (abs(v-w) .lt. tolerance) then
           eq=.true.
        else
           eq=.false.
        endif

   end function myequal_i
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
     real(r15), dimension(:), intent(in) ::  v,w
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
   function myequalarray_r(v,w) result(eq)
     real(r15), dimension(:), intent(in) ::  v,w
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
   end function myequalarray_r
!-----------------------------------
   function myequalarray_i(v,w) result(eq)
     integer, dimension(:), intent(in) ::  v,w
     logical eq !, intent(out) :: eq
     integer i,n
     i = size(v) ; n=size(w)
     if (n .ne. i) then
        print*, 'MYEQUALARRAY: the input arrays are of different size ',i,n
        stop
     else
        eq = .true.
        loop: do n=1,i
           if ( .not. myequal_i(v(n),w(n)) ) then
              eq = .false.
              exit loop
           endif
        enddo loop
     endif
   end function myequalarray_i
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
     real(r15), intent(in) :: w(3)
     add%x = v%x + w(1)
     add%y = v%y + w(2)
     add%z = v%z + w(3)
   end function addition_vav
!-----------------------------------
   function addition_avv(w,v) result(add)
     type(vector) add !, intent(out) ::
     type(vector), intent(in) :: v
     real(r15), intent(in) :: w(3)
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
     real(r15), intent(in) :: w(3)
     dif%x = v%x - w(1)
     dif%y = v%y - w(2)
     dif%z = v%z - w(3)
   end function subtraction_va
!-----------------------------------
   function subtraction_av(w,v) result(dif)
     type(vector) dif !, intent(out) ::
     type(vector), intent(in) :: v
     real(r15), intent(in) :: w(3)
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
   real(r15), intent(in) :: s
   type(vector), intent(in) :: v
   type(vector) sv !, intent(out)::
     sv%x=s*v%x
     sv%y=s*v%y
     sv%z=s*v%z
   end function multiply_by_real
!-----------------------------------
   function divide_by_scalar(v,s) result(sv)
   real(r15), intent(in) :: s
   type(vector), intent(in) :: v
   type(vector) sv !, intent(out)::
     sv%x=1/s*v%x
     sv%y=1/s*v%y
     sv%z=1/s*v%z
   end function divide_by_scalar
!-----------------------------------
   function distance(v,w) result(d)
   type(point), intent(in) :: v,w
   real(r15) d !, intent(out)::
     d = sqrt( (v%x-w%x)**2 +(v%y-w%y)**2 +(v%z-w%z)**2 )
   end function distance
!-----------------------------------
    subroutine point_eq_array(p,a)
      real(r15), intent(in) :: a(3)
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
      real(r15), intent(in) :: aa(3)
      type(vector), intent(out) :: vv
      vv%x = aa(1)
      vv%y = aa(2)
      vv%z = aa(3)
    end subroutine vector_eq_array
!-----------------------------------
    subroutine array_eq_vector(a,v)
      type(vector), intent(in) :: v
      real(r15), intent(out) :: a(3)
      a(1) = v%x
      a(2) = v%y
      a(3) = v%z
    end subroutine array_eq_vector
!-----------------------------------
    subroutine array_eq_point(a,p)
      type(point), intent(in) :: p
      real(r15), intent(out):: a(3)
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
     real(r15) l
     integer, intent(in) :: v(:)
     l = sqrt(real(dot_product(v,v)))
   end function lengthi
!-----------------------------------
   function lengthv(v) result(l)
     real(r15) l
     type(vector), intent(in) :: v
     l = sqrt(v%x*v%x+v%y*v%y+v%z*v%z)
   end function lengthv
!-----------------------------------
   function lengtha(v) result(l)
     real(r15) l
     real(r15), dimension(:), intent(in) :: v
!    integer i
     l = sqrt(dot_product(v,v))
!    l = 0d0
!    do i=1,size(v)
!       l = l + v(i)*v(i)
!    enddo
!    l = sqrt(l)
   end function lengtha
!-----------------------------------
 subroutine reduce_v(v,q1,q2,q3,w)
!! takes cartesian coordinates and returns direct coordinates
 implicit none
 type(vector) , intent(in) :: v,q1,q2,q3
 type(vector) , intent(out) :: w
 real(r15) a1,a2,a3

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
 real(r15) , intent(out) :: w(3)

 w(1) = v .dot. q1
 w(2) = v .dot. q2
 w(3) = v .dot. q3

 end subroutine reduce_v1
!-----------------------------------
 subroutine reduce_v2(v,q1,q2,q3,w)
!! takes cartesian coordinates and returns direct coordinates
 implicit none
 type(vector) , intent(in) :: q1,q2,q3
 real(r15) , intent(in) :: v(3)
 real(r15) , intent(out) :: w(3)

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
 real(r15), intent(in) :: a(3)
 real(r15) a1,a2,a3

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
 real(r15) , intent(in) :: a(3),q1(3),q2(3),q3(3)
 type(vector) w
 real(r15) a1,a2,a3

 a1 = a .dot. q1
 a2 = a .dot. q2
 a3 = a .dot. q3
 w%x = a1 ; w%y = a2 ; w%z = a3

 end subroutine reduce_a4
!-----------------------------------------
 subroutine reduce_a5(a,q1,q2,q3,w)
! takes cartesian coordinates and returns direct coordinates
 implicit none
 real(r15), intent(in) :: a(3),q1(3),q2(3),q3(3)
 real(r15), intent(out) :: w(3)

 w(1) = a .dot. q1
 w(2) = a .dot. q2
 w(3) = a .dot. q3

 end subroutine reduce_a5
!-----------------------------------------
 subroutine bring_to_center_v(v,r1,r2,r3,g1,g2,g3,w)
 use constants
! takes cart coordinates and returns cart coordinates within the FBZ
 implicit none
 type(vector), intent(in) :: v,g1,g2,g3,r1,r2,r3
 type(vector), intent(out) :: w
 real(r15) a1,a2,a3

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
 real(r15), intent(in) :: v(3)
 real(r15), intent(out) :: w(3)
 real(r15) a1,a2,a3

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
 real(r15) rx(3) !, intent(out)::
 rx(1)=r0%x
 rx(2)=r0%y
 rx(3)=r0%z
 end function v2a_r
!-----------------------------------------
 function v2a_a(r0) result(rx)
 implicit none
 type(vector), dimension(:), intent(in):: r0
 real(r15), dimension(3,size(r0)) :: rx !, intent(out)::
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
 real(r15), intent(in):: r0(3)
 rx%x=r0(1)
 rx%y=r0(2)
 rx%z=r0(3)
 end function a2v

 end module geometry
!==============================================
module ios
 use constants, only : r15
 integer, parameter:: uposcar=10,uparams=11,uborn=12, &
&         ufco=20,ulog=30,utraj=40,umatrx=50,umap=60,utimes=70,ufc=80, &
&         ufc1=21,ufc2=22,ufc3=23,ufc4=24,  &
&         ufit1=31,ufit2=32,ufit3=33,ufit4=34,  &
&         ujunk=79,uibz=80,uibs=81,ugrun=82,uband=83,ucor=93,utherm=110


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
  real(r15), dimension(:,:,:), intent(in) :: var

  l=len_trim(string)
  n1=size(var,1) ; n2=size(var,2) ; n3=size(var,3)
! write(unt,*)' write_outrm called;nl,n2,n3=',n1,n2,n3
  write(unt,4)string(1:l)//' is='
  do j=1,n1
     call write_out(unt,' block #j ',j)
     call write_out(unt,' ',var(j,:,:))
  enddo
4 format(a,99(1x,g11.4))
  end subroutine write_outrm3
!-----------------------------
  subroutine write_outrm(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer n,m,unt,l,i,j
  real(r15), dimension(:,:), intent(in) :: var

  l=len_trim(string)
  n=size(var,1) ; m=size(var,2)
! write(unit,*)' write_outrm called;nl,nc=',n,m
  write(unt,*)string(1:l)//' is='
  do i=1,n
     write(unt,4)var(i,:)
  enddo
4 format(99(1x,f13.5))
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
  real(r15), dimension(:) :: var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var(:)
4 format(a,99(1x,f13.5))
  end subroutine write_outrv
!-----------------------------
  subroutine write_outcv(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,l,i
  complex(r15), dimension(:) :: var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var(:)
4 format(a,99(3x,f13.5,1x,f13.5))
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
  real(r15) var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var
4 format(a,99(1x,f13.5))
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
4 format(a,3(1x,f13.5))
end subroutine write_outv

end module ios
!===========================================================
 module atoms_force_constants
! the type atom_id0 concerns the properties of the primitive unit cell
! it s atoms, masses, neighbors, shells and force constants, whereas the
! type atomic refers to atoms in the supercell, their displacements and forces.
 use constants, only : r15
 use geometry
 implicit none
 integer natom_super_cell,nmax, fc2flag

!-------------------------
 type cell_id         ! id of atoms: position within cell, and cell coordinates
    integer n(3)
    integer tau
    integer atomposindx  ! corresponding index in the atompos list
 end type
!-------------------------
 type shell
    integer no_of_neighbors  ! within that shell
    real(r15) radius              ! radius of that shell
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
    real(r15) mass,charge(3,3)
 end type
!-------------------------
 type atom_id0   ! all characteristics of an atom in the primitive central cell
    integer at_type,tau !  (n1,n2,n3)=0 for atoms in the primitive cell
    character name*2
    real(r15) mass,charge(3,3)
    integer nshells       ! how many shells are included for that atom = actual dim(shell)
    type(vector) equilibrium_pos
    type(shell), allocatable ::  shells(:)  ! of dim atom0(i)%nshells
!   type(shell), pointer ::  shells(:)  ! of dim atom0(i)%nshells
!    type(tensor2), pointer :: phi2(:)  ! refers to atom numbers neighbor of i
!    type(tensor3), pointer :: phi3(:,:)
!    type(tensor4), pointer :: phi4(:,:,:)
 end type
!-------------------------
! everything with 0 refers to the primitive cell

 integer  natom_type
 integer, allocatable:: natom(:),atom_type(:), map_rtau_sc(:,:)!label_of_primitive_atoms(:),
 real(r15), allocatable:: atompos0(:,:),mas(:),force(:,:,:),displ(:,:,:),vel(:,:),energy(:)
 real(r15), allocatable:: cur(:,:)
 character(2), allocatable:: atname(:)
 type (atom_id0), allocatable :: atom0(:)
 type (atomic), allocatable :: atom_sc(:),atom_shl(:)

!---------------------------
! Module for saved data
!      module force_constants_module
! maximum shells of nearest neighbors (along radial direction) , and actual # of neighborshells
      integer maxneighbors,maxshells
!     parameter(maxneighbors=18 )
! maximum number of atoms out to maxneighbors
      integer maxatoms,imaxatm,imaxnei
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
      real(r15), allocatable :: atompos(:,:)
! iatomcell(j,i), linear combination of basis vectors of the primitive lattice
! that takes us to the unit cell containing the ith atom
      integer, allocatable ::  iatomcell(:,:),iatomcell0(:)
! iatomcell0(i), identity of atom in unit cell at origin equivalent to ith atom
! iatomneighbor(j,i), nearest neighbor shell of jth atom in primitive unit
! cell at the origin that contains the ith atom
      integer iatomneighbor(:,:)
      allocatable iatomneighbor,iatomop,atomopfract
!      end module force_constants_module
! -------------------------------------

contains

! subroutine allocate_disp_forc(n)  ! needed for the md code
!    integer n
!   allocate( disp(3,n),forc(3,n),vel(3,n),cur(3,n),ncel1(n),ncel2(n),ncel3(n) )
!    allocate( disp(3,n),forc(3,n),vel(3,n),cur(3,n) )
! end subroutine allocate_disp_forc

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
!subroutine set_neighbor_list (maxshell)
!  Generates a set of neighbors to the primitive cell (up to maxshell, which should be large enough)
!  and finds which shell corresponds to the cutoff rcut set by the WS cell
!  this is needed in order to know what force constants to associate
!  to a pair ij, or ijk or ijkl
!  if a neighbor j is known then one can get the vect(i,j) between
!  their equilibrium positions
! use atoms_force_constants , only : natom_prim_cell
!use params
!use lattice
!use ios
!use geometry
!implicit none
!! integer, parameter :: max=500
!integer, intent(in) :: maxshell
!real(r15), intent(in) :: rcut
!integer i0,j,shel_count,counter,msort(maxatoms),l
!real(r15) dist(maxatoms),d_old,d_min,rmax

!  allocate( atom0(1:natom_prim_cell)%shells(maxshells) )
! rmax = 0 ; dist = 1d10
! i0loop: do i0=1,natom_prim_cell
!   allocate( atom0(i0)%shells(0:maxshell) )
!   do j=1,natoms
!      call calculate_distance(i0,j,atompos,maxatoms,dist(j) )
!      if ( iatomneighbor(i0,j) .eq. nshells(2) .and. dist(j) .gt. rmax ) then
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
!      atom0(i0)%shells(shel_count)%radius = d_min
!      atom0(i0)%shells(shel_count)%neighbors(counter)%tau = iatomcell0(msort(j))
!      atom0(i0)%shells(shel_count)%neighbors(counter)%n   = iatomcell(:,msort(j))
!      if(verbose) then
!        write(ulog,*)' ------------------------------------'
!         write(ulog,5)' count, j, of tau,n1,n2,n3,shel#,nb#,dij,rj =',j,msort(j),iatomcell0(msort(j)),   &
!                iatomcell(:,msort(j)),iatomneighbor(i0,msort(j)),counter,d_min,(atompos(l,msort(j)),l=1,3)
!      endif
!   enddo jloop

!   do shel_count = 0 , maxshell
!      write(ulog,*)'shell#, neigh#=',shel_count,atom0(i0)%shells(shel_count)%no_of_neighbors
!   enddo
!  also initialize atom0%equilibrium_pos

!  shift the position of the first atom to the  origin
! do i=1,natom_prim_cell
!    atompos(:,i) = atompos(:,i)-atompos(:,1)
!     atom0(i0)%equilibrium_pos = atompos0(1,i0)*r01 + atompos0(2,i0)*r02  &
! &                             + atompos0(3,i0)*r03
!   atom0(i0)%equilibrium_pos%x = atompos(1,i0)
!   atom0(i0)%equilibrium_pos%y = atompos(2,i0)
!   atom0(i0)%equilibrium_pos%z = atompos(3,i0)
! enddo
!  they might not be inside the primitive cell

!enddo i0loop
! rcut(2) = rmax

! format(a,8(2x,i4))
! format(a,' (',i4,2(1x,i4),1x,')')
! format(a,1x,i4,' (',f8.4,2(1x,f8.4),1x,')')
! format(a,3(2x,i4),' (',i4,2(1x,i4),1x,')',2(1x,i3),1x,f10.7,2x,9(1x,f9.6))
! format(a,1x,f10.7,2x,9(1x,f9.6))

!end subroutine set_neighbor_list

!--------------------------------
     subroutine test_ntau
     use ios, only : ulog
     implicit none
     integer i,tau1,n1(3),tau2,n2(3),iatom1,iatom2,err

     err=0
     do i=1,natoms
    !     call findatom_cart(v,iatom(irank))
        call get_n_tau_r    (atompos(:,i),n1,tau1)
        call get_n_tau_r_mod(atompos(:,i),n2,tau2)
        call find_atompos(n1,tau1,iatom1)
        call find_atompos(n2,tau2,iatom2)
        if(.not.(n1.myeq.n2) .or. tau1.ne.tau2) then
           write(*   ,4)'TEST_NTAU: i,i1,i2,tau1,tau2, n1,n2=', &
&                   i,iatom1,iatom2,tau1,tau2,n1,n2
           write(ulog,4)'TEST_NTAU: i,i1,i2,tau1,tau2, n1,n2=', &
&                   i,iatom1,iatom2,tau1,tau2,n1,n2
           err=1
        endif
     enddo
     if(err.eq.0) then
        write(ulog,*) 'Both get_n_tau_r and get_n_tau_r_mod yield the same results '
     else
        write(ulog,*) ' CAREFUL: get_n_tau_r and get_n_tau_r_mod yield different results! '
     endif
 4   format(a,3(1x,i5),3x,2i3,2(4x,3(1x,i2)))

     end subroutine test_ntau

 end module atoms_force_constants
!===========================================================
 module lattice
 use geometry
 use constants
 implicit none
 type(vector) gs1,gs2,gs3,rs1,rs2,rs3  ! rri=translation vectors of the supercell
 type(vector) r1conv,r2conv,r3conv,g1conv,g2conv,g3conv
 type(vector) r01,r02,r03,g01,g02,g03  ! tr vect of prim cell and its recip spce
 real(r15) volume_r,volume_g,volume_r0,volume_g0
 real(r15) lattice_parameter,latticeparameters(6),primitivelattice(3,3)
 real(r15),dimension(3,3):: conv_to_cart,cart_to_conv,prim_to_conv,conv_to_prim,  &
  &                       prim_to_cart,cart_to_prim
 real(r15) box(3),boxg(3),density
 real(r15) r0g(3,3)
! real(r15), allocatable:: rws_weights(:),gws_weights(:)
! integer nr1(3),nr2(3),nr3(3)
! real(r15), allocatable:: rgrid(:,:),ggrid(:,:),xgrid(:,:)
 real(r15) gws26(3,26),rws26(3,26),invn_sc(3,3)  ! superlattice shells defining the WS of SL
 real(r15) g0ws26(3,26),r0ws26(3,26)
 integer n_sc(3,3)

  interface make_rg
     module procedure make_rga,make_rgv
  end interface

  interface check_int
     module procedure check3,check_a
  end interface

  interface cart_to_direct
     module procedure cart_to_direct_aa,cart_to_direct_av,cart_to_direct_v
  end interface

  interface direct_to_cart
     module procedure direct_to_cart_aa,direct_to_cart_av,direct_to_cart_v
  end interface

  interface bring_to_cell_c
     module procedure bring_to_cell_cv,bring_to_cell_ca
  end interface

  interface bring_to_cell_d
     module procedure bring_to_cell_dv,bring_to_cell_da
  end interface

 contains

!-------------------------------------------
 subroutine make_r0g
! matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
! matmul(r0g,v) gives the 3 reduced coordinates of v in units of the supercell translations 
! coordinates are v(1)*r01+v(2)*r02+v(3)*r03
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
! subroutine make_rg(x1,x2,x3,q1,q2,q3,n)
!  matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
!  matmul(r0g,v) gives the 3 reduced coordinates of an atom in the supercell if its
!  reduced coordinates in the primitive cell are given by v
! use geometry
! implicit none
! real(r15) x1(3),x2(3),x3(3),q1(3),q2(3),q3(3),n(3,3)

! n(1,1) = x1 .dot. q1
! n(1,2) = x1 .dot. q2
! n(1,3) = x1 .dot. q3
! n(2,1) = x2 .dot. q1
! n(2,2) = x2 .dot. q2
! n(2,3) = x2 .dot. q3
! n(3,1) = x3 .dot. q1
! n(3,2) = x3 .dot. q2
! n(3,3) = x3 .dot. q3

! end subroutine make_rg
! -------------------------------------------
 subroutine check_a(r,a,ier,g1,g2,g3)
!! subroutine to check whether r is an integer multiple of (r01,r02,r03)
!! output ai are the coefficients of its linear combination on this basis
!! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use constants, only : pi
 use ios
 implicit none

 type(vector) :: r,g1,g2,g3
 real(r15) a(3)
 integer ier,n(3)
! logical is_integer

! lengths smaller than tolerance are considered to be zero
!call write_out(ulog,'CHECK: R ',r)
 a(1) = (r .dot. g1)/(2*pi)  ! g's have a 2 pi
 a(2) = (r .dot. g2)/(2*pi)
 a(3) = (r .dot. g3)/(2*pi)
 n=floor(a)
 if(is_integer_a(a)) then
!  write(ulog,*) ' R is a multiple of r0s'
!  write(ulog,3) ' n1,n2,n3  =',a1,a2,a3
   ier = 0
 else
   ier = 1
!  write(ulog,*) ' R is not a multiple of r0s , check your inputs '
!  write(ulog,3) 'ier=1;r and r.gi/2pi=',r,a
!  stop
 endif
3 format(a,9(1x,g11.4))

 end subroutine check_a
!-------------------------------------------
 subroutine check3(r,a1,a2,a3,ier,g1,g2,g3)
!! subroutine to check whether r is an integer multiple of (r01,r02,r03)
!! output ai are the coefficients of its linear combination on this basis
!! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use constants, only : pi
 use ios
 implicit none

 type(vector) :: r,g1,g2,g3
 real(r15) a1,a2,a3
 integer ier

! lengths smaller than tolerance are considered to be zero
! call write_out(ulog,'CHECK: R ',r)
 a1 = (r .dot. g1)/(2*pi)  ! g's have a 2 pi
 a2 = (r .dot. g2)/(2*pi)
 a3 = (r .dot. g3)/(2*pi)
 if(abs(a1-nint(a1)).ge.tolerance .or.abs(a2-nint(a2)).ge.tolerance   &
&   .or.abs(a3-nint(a3)).ge.tolerance ) then
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

 end subroutine check3
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
 real(r15) a1,a2,a3
 integer ier

!call write_out(ulog,'CHECK: R ',r)
 a1 = r%x
 a2 = r%y
 a3 = r%z
 if(abs(a1-nint(a1)).ge.tolerance .or.abs(a2-nint(a2)).ge. tolerance   &
&   .or.abs(a3-nint(a3)).ge.tolerance ) then
   ier = 1
   write(ulog,3) 'ier=1; r is not integer ',a1,a2,a3
 else
   ier = 0
 endif
3 format(a,9(1x,g11.4))

 end subroutine check_d
!-----------------------------------------
 subroutine transform_input_structure
!! takes input variables latticeparameters(6),primitivelattice(3,3) to construct the primitive cell
 use ios
! use constants, only : pi
! use atoms_force_constants, only : xmatinv
 implicit none
 external xmatinv
 real(r15) latp(6)
 integer ier

!k1
!      d2r = 4d0*datan(1d0)/180d0 ! 3.1415926535897932384626/180d0
!k1
      conv_to_cart=0
      prim_to_conv=primitivelattice
!     call write_out(6,'Prim_to_conv=prim vectors in reduced units of conv in lines ',prim_to_conv)
!     call write_out(ulog,'Prim_to_conv=prim vectors in reduced units of conv in lines ',prim_to_conv)
      call xmatinv(3,prim_to_conv,conv_to_prim,ier)
!     call write_out(6,'conv_to_prim has primitive reciprocal vectors in units of recip conv ',conv_to_prim)
!     call write_out(ulog,'conv_to_prim has primitive reciprocal vectors in units of recip conv ',conv_to_prim)
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
      write(*,*)'conv_to_cart has conventional vectors on its columns='
      write(*,*)'prim_to_cart has primitive vectors on its columns'
      write(*,*)'cart_to_conv has conventional reciprocal vectors on its lines (no 2pi)'
      write(*,*)'cart_to_prim has primitive reciprocal vectors on its lines (no 2pi) '
      call write_out(ulog,' Conventional cell(conv_to_cart) on columns ',conv_to_cart)
! go to reciprocal space by matrix inversion
      call xmatinv(3,conv_to_cart,cart_to_conv,ier)
      call write_out(ulog,' Conv. reciprocal lattice on lines  (no 2pi)',cart_to_conv)

! the usual way
      r1conv=conv_to_cart(:,1)
      r2conv=conv_to_cart(:,2)
      r3conv=conv_to_cart(:,3)
! this has an extra factor of 2pi wrt cart_to_conv
      call make_reciprocal_lattice_2pi(r1conv,r2conv,r3conv,g1conv,g2conv,g3conv)
      call write_out(ulog,' Conventional reciprocal(1) with 2pi',g1conv)
      call write_out(ulog,' Conventional reciprocal(2) with 2pi',g2conv)
      call write_out(ulog,' Conventional reciprocal(3) with 2pi',g3conv)
!     call write_out(ulog,' Conv recip on lines/2pi ',cart_to_conv)
! prim_to_cart(i,j) is the ith cartesian coordinate of the jth primitive translation vector
      prim_to_cart=matmul(conv_to_cart,prim_to_conv)
      call write_out(ulog,' Primitive translations: prim_to_cart ( on columns ) ',prim_to_cart)

      r01 = prim_to_cart(:,1)
      r02 = prim_to_cart(:,2)
      r03 = prim_to_cart(:,3)

      call make_reciprocal_lattice_2pi(r01,r02,r03,g01,g02,g03)

!     call write_out(ulog,' g01= ' ,g01)
!     call write_out(ulog,' g02= ' ,g02)
!     call write_out(ulog,' g03= ' ,g03)

      call xmatinv(3,prim_to_cart,cart_to_prim,ier)
      if(ier.ne.0)then
        write(6,*)'Error in transform_input_structure: primitive_lattice '   &
     &       //'is singular'
        stop
      endif
      call write_out(ulog,' Primitive reciprocal: cart_to_prim (no 2pi, on lines)  ',cart_to_prim)

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
 use constants
 implicit none
 integer i
 real(r15)particle_density,volcut  ! x1,x2,x3,
 external unitcell,symmetry_init

! generate r0i and g0i rom latticeparameters, primitivelattice
  call transform_input_structure
  call calculate_volume(r01,r02,r03,volume_r0)
  call calculate_volume(g01,g02,g03,volume_g0)

  write(ulog,*)' PRIMITIVE CELL *************************'
  call write_out(ulog,' r01',r01)
  call write_out(ulog,' r02',r02)
  call write_out(ulog,' r03',r03)
  call write_out(ulog,' Volume of primitive cell ',volume_r0)
  call write_out(ulog,' g01',g01)
  call write_out(ulog,' g02',g02)
  call write_out(ulog,' g03',g03)
  call write_out(ulog,' Volume of reciprocal primitive cell ',volume_g0)
  call write_out(6,' primitive lattice (lines)',transpose(prim_to_cart))
  call write_out(6,' reciprocal lattice (lines, no 2pi) ',cart_to_prim)


  box(1) = length(r01)
  box(2) = length(r02)
  box(3) = length(r03)
  boxg(1) = length(g01)
  boxg(2) = length(g02)
  boxg(3) = length(g03)
  write(ulog,3)' Lengths of r0i = ',box
  write(ulog,3)' Lengths of g0i = ',boxg
! write(ulog,3)' volume_r0,volume_g0 = ',volume_r0,volume_g0

  write(ulog,*)' Reduced Atomic positions in conventionals, mass ====================== '
  do i=1,natom_prim_cell
     atom0(i)%equilibrium_pos=atompos0(:,i)
     write(ulog,8)atom0(i)%name,atom0(i)%at_type,atom0(i)%tau,atompos0(:,i),atom0(i)%mass
  enddo

! convert to cartesian
  write(ulog,*)' reduced Atomic positions in primitive, mass ====================== '
  do i=1,natom_prim_cell
     atompos0(:,i)=matmul(conv_to_cart,atompos0(:,i))
     write(ulog,8)' '//atom0(i)%name,atom0(i)%at_type,  &
&    atom0(i)%tau,matmul(cart_to_prim,atompos0(:,i)),atom0(i)%mass
! bring to [0 r0i[
     call unitcell(cart_to_prim,prim_to_cart,atompos0(:,i), atompos0(:,i))
     atom0(i)%equilibrium_pos = atompos0(:,i)
  enddo

  write(ulog,*)' Cartesian Atomic positions, mass ====================== '
  do i=1,natom_prim_cell
     write(ulog,8)'After unitcell '//atom0(i)%name,atom0(i)%at_type,atom0(i)%tau,atom0(i)%equilibrium_pos,atom0(i)%mass
  enddo

! generate neighbors in atompos up to a disctance rcut(2)
!  write(ulog,*)'A default cutoff range of ',rcut(2),' has been chosen for FC2s'

! get largest # of atoms within rcut; and nsmax=rcut/shortest translation+2
!  call get_upper_bounds(r01,r02,r03,rcut(2),maxatoms,nsmax)

  density=sum(atom0(:)%mass)/volume_r0*uma*1d27
  write(ulog,*)' density (uma/Ang^3)=',sum(atom0(:)%mass)/volume_r0
  write(ulog,*)' density (  g/ cm^3)=',density
  write(ulog,*)' COORDINATES read successfully '
  write(   *,*)' COORDINATES read successfully '

! hardwire maxatoms and maxneighbors; after fc_init maxatoms -> natoms, maxneighbors -> nshell
  particle_density=natom_prim_cell/volume_r0
  volcut=4*pi/3d0*rcutoff**3   ! volume of a sphare of radius rcutoff 
  maxatoms = nint(volcut * particle_density)
  write(ulog,*)'Rcutoff=',rcutoff,' suggests maxatoms=',maxatoms
  maxatoms=max(maxatoms,1500)
  write(ulog,*)' the value of maxatom will be updated to ',maxatoms

!   rcutoff=(3*maxatoms/particle_density/4/pi)**0.333333
!   write(ulog,*)'Initially, take Maxatoms=',maxatoms,' suggests Rcutoff=',rcutoff

! FCs will exist at the most up to maxshells
  maxshells=min(maxval(nshells)+2,20)  ! default for largest number of cubic shells < 20
!  everything is based on maxneighbors, maxatoms and maxshells will be updatedaccordingly
  maxneighbors=min(3*maxshells,((2*maxshells)**3)/48)

! maxshells=50
! maxneighbors=800

  write(ulog,*) 'Default maxatoms    =', maxatoms
  write(ulog,*) 'Default maxshells   =', maxshells
  write(ulog,*) 'Default maxneighbors=', maxneighbors
  write(ulog,*)'*****************************************************************'

  imaxatm=1
     if(allocated(atompos   )) deallocate(atompos   )
     if(allocated(iatomcell )) deallocate(iatomcell )
     if(allocated(iatomcell0)) deallocate(iatomcell0)
     allocate(atompos(3,maxatoms), iatomcell(3,maxatoms),iatomcell0(maxatoms))

  imaxnei=1
! we'll fix maxatoms since the radius=15 is large enough, and adjust maxneighbors
  do while(imaxnei.eq.1)
     maxneighbors=maxneighbors+100
     write(*,*)' Increasing maxeighbors to ',maxneighbors

! find symmetry operations; the above allocation is for this subroutine
! requires maxneighbors to generate atoms within these shells
     call symmetry_init(natom_prim_cell,atom_type,atompos0)
  enddo
!    ----------------        DO not want to increase maxshells
! do while(imaxshl.eq.1)
!    maxshells=maxshells+5
!    write(*,*)' updating maxshells to ',maxshells
!    call force_constants_init(natom_prim_cell,atom_type,atompos0)
! enddo
  write(*   ,*)' maxatoms, natoms=',maxatoms,natoms
  write(ulog,*)' maxatoms, natoms=',maxatoms,natoms

  call test_ntau
! cat_to_prim had reciprocal lattice vectors on its lines
! call write_out(ulog,'reciprocal lattice vectors ',cart_to_prim)


3 format(a,9(2x,g13.6))
4 format(i6,3(2x,g12.5),4(3x,i4))
6 format(a,(2x,g12.5),a,4(3x,i6))
8 format(a,2i3,9(1x,g13.6))

 end subroutine make_unitcell
!-----------------------------------------
! subroutine get_components(q,n,i,j,k )  !,g1,g2,g3)
!! for a given q-vector, it finds its integer components assuming it was
!! created as: nk=0 do i=0,N(1)-1 ; do j=0,N(2)-1; do k=0,N(3)-1; nk=nk+1 ;q=i*G1/N1+j*G2/N2+k*G3/N3
!! as the ones created in make_kp_reg with zero shift
!! it works even if there's a shift less than 0.5
! implicit none
! real(r15) q(3) !,g1(3),g2(3),g3(3),r1(3),r2(3),r3(3)
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
 real(r15), intent(in ) :: x1(3),x2(3),x3(3),q1(3),q2(3),q3(3)
 real(r15), intent(out) :: n(3,3)

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
 real(r15), intent(out) :: n(3,3)

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
 subroutine set_neighbor_list
!! Generates a set of neighbors to the primitive cell up to
!! and calculates the corresponding number of neighborshells
! this is needed in order to know what force constants to associate
! to a pair ij, or ijk or ijkl
! if a neighbor j is known then one can get the vect(i,j) between
! their equilibrium positions
 use atoms_force_constants !, only : natom_prim_cell, atom0,maxatoms
 use params
 use ios
 use geometry
 implicit none
 integer i0,j,jj,s,l,k,shel_count,minat,cnt,nbs,mxshell
 integer nb(0:maxneighbors),msort(maxatoms)  ! because dist is of size natoms
 real(r15) dist(maxatoms),d_old,d_min(0:maxneighbors),rmax,dij,d_new

 call write_out(ulog,'SET_NEIGHBOR_LIST: entering with maxneighbors ',maxneighbors)
 call write_out(ulog,'SET_NEIGHBOR_LIST: entering with maxshells    ',maxshells )
 call write_out(ulog,'SET_NEIGHBOR_LIST: entering with natoms       ',natoms)
 write(ulog,*)' the neighbor list will stop at maxshells '
! allocate( atom0(natom_prim_cell)%shells(0:mxshell)) !%neighbors(maxneighbors) )
! allocate( atom0(natom_prim_cell)%shells(0:mxshell)%neighbors(maxneighbors) )
 rmax = 0  ! is the largest radius among the natoms neighbors
 mxshell = 0

 i0loop: do i0=1,natom_prim_cell
    dist = 1d10
! calculate pair i0-j distances, and sort them
    cnt=0
    do j=1,natoms
       call calculate_distance(i0,j,atompos,natoms,dij)  ! dij=|atompos(i0)-atompos(j)|
       dist(j)=dij
    enddo
    call sort(natoms,dist,msort,maxatoms)

! count the corresponding # of shells and nb, # of atoms within that shell, for each i0
    shel_count = -1
    d_old = 1e4
    d_min=0
    nb=1
    minat=min(natoms,350)  ! about 7 shells at the most: 7^3=343
    write(ulog,*)' ========================================================='
    write(ulog,*)' first ',minat,' neighborlist of atom i0 =',i0

    jloop: do j=1,natoms !minat
       d_new = dist(msort(j))  ! start in order of increasing distances from d=0

       if ( abs(d_new-d_old).gt.tolerance ) then   ! new shell
          if(shel_count.eq.maxneighbors) exit jloop
          shel_count = shel_count+1
          nb(shel_count)=1
          d_old = d_new
          d_min(shel_count)=d_new
          write(*,4)' NEW shell , dij= ',shel_count,d_new
          if ( shel_count .gt. maxneighbors ) then
             write(   *,*)'maxneighbors=', maxneighbors,' shells completed up to distance ',d_new
             write(   *,*)'It  corresponds to ',msort(j),' th nearest atom'
             write(   *,*)'The largest # of shells is set to ',shel_count-1 ,' which should be ',maxneighbors
             write(ulog,*)'maxneighbors=', maxneighbors,' shells completed up to distance ',d_new
             write(ulog,*)'It  corresponds to ',msort(j),' th nearest atom'
             write(ulog,*)'The largest # of shells is set to ',shel_count-1 ,' which should be ',maxneighbors
             exit jloop
          endif

          iatomneighbor(i0,msort(j))=shel_count

! after maxshell iatomneighbor values are not reliable ! FIXED !

       else  ! same shell
          nb(shel_count)=nb(shel_count)+1  ! counts how many atoms within shell=shel_count
          if(j.lt.150) write(*,9)'shell# , nb# , j , dij= ',shel_count,nb(shel_count),j,d_new
          iatomneighbor(i0,msort(j))=shel_count
       endif
    enddo jloop

    shel_count=shel_count-1
    atom0(i0)%nshells=shel_count
    write(*,*)'allocating atom0(',i0,')%shells(:) to size ',shel_count
    if(allocated(atom0(i0)%shells)) deallocate(atom0(i0)%shells)
    allocate(atom0(i0)%shells(0:shel_count))  ! %neighbors dimension is fixed to 296
    cnt=0  ! counts atoms from nearest to farthest, shell by shell
    do s=0,shel_count
       atom0(i0)%shells(s)%no_of_neighbors = nb(s)
       atom0(i0)%shells(s)%radius = d_min(s)
       do l=1,nb(s)  ! l is the index of the neighbor of i0 in that shell

          cnt=cnt+1  ! counts the total # of atoms up to this shell
          jl:do j=1,natoms
       !     if(cnt.eq.msort(j)) then
             jj=msort(j)
             if(d_min(s) .myeq. length(atompos(:,jj)-atompos(:,i0))) then
! make sure it is not already assigned
               do k=1,l-1
                  if( jj.eq.atom0(i0)%shells(s)%neighbors(k)%atomposindx) cycle jl
               enddo

               atom0(i0)%shells(s)%neighbors(l)%tau = iatomcell0(jj)
               atom0(i0)%shells(s)%neighbors(l)%n   = iatomcell(:,jj)
               atom0(i0)%shells(s)%neighbors(l)%atomposindx = jj
      !        iatomneighbor(i0,msort(j))= s ! atom0(i0)%shells(s)%neighbors(l)
!    if(verbose) then
!        write(ulog,*)' ------------------------------------'
    write(ulog,5)'shl,nb#,atomj,tauj,nj,shel#,dij,Rj =',s,l,jj,iatomcell0(jj),   &
    &             iatomcell(:,jj),iatomneighbor(i0,jj),d_min(s),atompos(:,jj)
!    endif
               exit jl
             endif
          enddo jl
       enddo
!      nbs= atom0(i0)%shells(s)%no_of_neighbors
!      write(ulog,*)atom0(i0)%shells(s)%neighbors(1:nbs)
    enddo

! mxshell is the largest nearest neighbor shell number, and will replace initially guessed maxneighbors
    if(mxshell.le.shel_count) then
        mxshell=shel_count
        write(*   ,*)'i0,shel_count,mxshell=',i0,shel_count,mxshell
        write(ulog,*)'i0,shel_count,mxshell=',i0,shel_count,mxshell
    endif

    do s = 0 , shel_count
       nbs= atom0(i0)%shells(s)%no_of_neighbors
       write(ulog,7)'shell#, neigh#, rij/rij(1),rij=',s,atom0(i0)%shells(s)%no_of_neighbors   &
            ,atom0(i0)%shells(s)%radius/atom0(i0)%shells(1)%radius,atom0(i0)%shells(s)%radius
       write(ulog,*)'    nb#,jatom,tau,     n    '
       do l=1,nbs
          write(ulog,8)' --  ',l,atom0(i0)%shells(s)%neighbors(l)%atomposindx  &
&                ,atom0(i0)%shells(s)%neighbors(l)%tau  &
&                ,(atom0(i0)%shells(s)%neighbors(l)%n(j),j=1,3)
       enddo
    enddo

! also initialize atom0%equilibrium_pos
    atom0(i0)%equilibrium_pos%x = atompos(1,i0)
    atom0(i0)%equilibrium_pos%y = atompos(2,i0)
    atom0(i0)%equilibrium_pos%z = atompos(3,i0)

 enddo i0loop

 write(ulog,*)'SET_NEIGHBOR_LIST: largest#of neighbor shells=', mxshell
 maxneighbors=mxshell
 write(ulog,*)'# Maxneighbors updated to ',maxneighbors


2 format(a,8(2x,i4))
3 format(a,' (',i4,2(1x,i4),1x,')')
4 format(a,1x,i4,' (',f8.4,2(1x,f8.4),1x,')')
5 format(a,4(1x,i4),' (',i2,2(1x,i2),')',1x,i4,9(1x,f9.4))
6 format(a,1x,f10.7,2x,9(1x,f9.6))
7 format(a,1x,i3,1x,i4,2(1x,f11.5),2x,999('[',i5,',',i2,'(',3(i2),')]'))
8 format(a,i3,i5,1x,i2,' (',i2,2(1x,i2),')')
9 format(a,3i4,3x,9(1x,f9.3))

 end subroutine set_neighbor_list
!-------------------------------------
 subroutine cart_to_direct_v(v,w)
! takes cart coordinates and returns cart coordinates within the supercell
! use lattice
 use geometry
 use constants, only : pi,r15
 implicit none
 type(vector), intent(in) :: v
 type(vector), intent(out) :: w

 w%x = (v .dot. g01)/(2*pi)
 w%y = (v .dot. g02)/(2*pi)
 w%z = (v .dot. g03)/(2*pi)

 end subroutine cart_to_direct_v
!-----------------------------------
 subroutine cart_to_direct_av(a,w)
! takes cart coordinates and returns cart coordinates within the supercell
! use lattice
 use geometry
 use constants, only : r15
 implicit none
 type(vector), intent(out) :: w
 real(r15), intent(in) ::  a(3)
 type(vector) v

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call cart_to_direct_v(v,w)

 end subroutine cart_to_direct_av
!-----------------------------------
 subroutine cart_to_direct_aa(a,b)
! takes cart coordinates and returns cart coordinates within the supercell
 use geometry
! use lattice
 use constants, only : r15
 implicit none
 real(r15), intent(in) ::  a(3)
 real(r15), intent(out) ::  b(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call cart_to_direct_v(v,w)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end subroutine cart_to_direct_aa
!-----------------------------------
 subroutine direct_to_cart_v(v,w)
! takes direct coordinates and returns cart coordinates
 use geometry
! use lattice
 implicit none
 type(vector), intent(in) :: v
 type(vector), intent(out) :: w

 w%x = v%x*r01%x + v%y*r02%x + v%z*r03%x
 w%y = v%x*r01%y + v%y*r02%y + v%z*r03%y
 w%z = v%x*r01%z + v%y*r02%z + v%z*r03%z

 end subroutine direct_to_cart_v
!-----------------------------------
 subroutine direct_to_cart_av(a,w)
! takes direct coordinates and returns cart coordinates
 use geometry
! use lattice
 use constants, only : r15
 implicit none
 real(r15), intent(in) ::  a(3)
 type(vector), intent(out) :: w
 type(vector) v

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call direct_to_cart_v(v,w)

 end subroutine direct_to_cart_av
!-----------------------------------
 subroutine direct_to_cart_aa(a,b)
! takes direct coordinates and returns cart coordinates
 use geometry
! use lattice
 use constants, only : r15
 implicit none
 real(r15), intent(in) ::  a(3)
 real(r15), intent(out) ::  b(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call direct_to_cart_v(v,w)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end subroutine direct_to_cart_aa
!------------------------------------
 subroutine dir2cart_g(q,k)
! takes q in direct coordinates and outputs k in cartesian
 use geometry
! use lattice
 use constants, only : r15
 real(r15), intent(in) ::  q(3)
 real(r15), intent(out) ::  k(3)

 k(:) = q(1)*g01 + q(2)*g02 + q(3)*g03

 end subroutine dir2cart_g
!------------------------------------
 subroutine dirconv2cart(q,k)
! takes q in direct coordinates of the conventional cell and outputs k in cartesian
 use geometry
! use lattice
 use ios
 use constants, only : r15
 real(r15), intent(in) ::  q(3)
 real(r15), intent(out) ::  k(3)

 k(:) = q(1)*g1conv + q(2)*g2conv + q(3)*g3conv

! k2=matmul(cart_to_prim,q)
! call write_out(ulog,'k1=',k1)
! call write_out(ulog,'k2=',k2)

 end subroutine dirconv2cart
!-----------------------------------------
 function bring_to_cell_dv(v) result(w)
! takes direct coordinates and returns direct coordinates
 implicit none
 type(vector) v,w
 real(r15) a1,a2,a3

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
 real(r15) a(3)

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
 real(r15), intent(in) :: a(3)
 real(r15), intent(out) :: b(3)
 type(vector), intent(in) :: g1,g2,g3,r1,r2,r3
 type(vector) v,w
 real(r15) a1,a2,a3

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
 real(r15) a1,a2,a3

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

 end module lattice
!============================================================
  module svd_stuff
 use constants, only : r15
! stuff related to Harold's subroutines, mappings and the svd matrices

 type groupmatrix
    real(r15), allocatable:: mat(:,:)
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
 real(r15), allocatable :: mapmat(:,:,:)
 real(r15) svdcut !,rayon(4)
 integer, allocatable:: iatomterm_1(:,:),ixyzterm_1(:,:),igroup_1(:),map_1(:)
 integer, allocatable:: iatomterm_2(:,:),ixyzterm_2(:,:),igroup_2(:),map_2(:)
 integer, allocatable:: iatomterm_3(:,:),ixyzterm_3(:,:),igroup_3(:),map_3(:)
 integer, allocatable:: iatomterm_4(:,:),ixyzterm_4(:,:),igroup_4(:),map_4(:)
 character(1), allocatable:: err_1(:),err_2(:),err_3(:),err_4(:)
 type(fulldmatrix) map(4)
 real(r15), allocatable:: fcrnk(:,:)
 real(r15), allocatable:: ampterm_1(:),fcs_1(:)
 real(r15), allocatable:: ampterm_2(:),fcs_2(:),grun_fc(:)
 real(r15), allocatable:: ampterm_3(:),fcs_3(:)
 real(r15), allocatable:: ampterm_4(:),fcs_4(:)
 real(r15), allocatable:: amat(:,:),bmat(:),sigma(:),fcs(:),ahom(:,:),kernelbasis(:,:)
 real(r15), allocatable:: a11ia12(:,:),fc1(:)
 real(r15), allocatable:: atransl(:,:),btransl(:),arot(:,:),brot(:),ahuang(:,:),bhuang(:), &
 &                      aforce(:,:),bforce(:),energies(:)
 real(r15), allocatable:: fc_ind(:) ! This is for getting force constant by reading FC2.dat file
 integer inv_constraints,force_constraints,dim_al,dim_ac,n_indep,newdim_al,dim_hom   &
&        ,transl_constraints, rot_constraints, huang_constraints,nindepfc
contains

! subroutine set_maxterms
!   maxterms(1)=40 !100
!   maxterms(2)=4000 !500
!   maxterms(3)=5000 !1800
!   maxterms(4)=3000 !2000
!   maxtermzero(1)=20 !500
!   maxtermzero(2)=1000 !2000
!   maxtermzero(3)=5000!5000
!   maxtermzero(4)=3000 !8000
!   maxtermsindep(1)=10
!   maxtermsindep(2)=150
!   maxtermsindep(3)=150
!   maxtermsindep(4)=300
!   maxgroups(1)=10
!   maxgroups(2)=100
!   maxgroups(3)=150
!   maxgroups(4)=300
! end subroutine set_maxterms

!--------------------------------------------
 subroutine remove_zeros_all(nl,nc)
! removes lines which are all zero from the amatrix
 use ios
 use params
! use svd_stuff  !for itrans,irot,ihuang
 implicit none
 integer, intent(in) :: nc
 integer, intent(inout) :: nl
! real(r15), dimension(:,:), intent(inout) :: amatr
! real(r15), dimension(:)  , intent(inout) :: bmatr
! real(r15), intent(inout) :: amatr(nl,nc)
! real(r15), intent(inout) :: bmatr(nl)
 integer i,j,nz,tra,rot,hua,tra_out,rot_out,hua_out
 real(r15), allocatable:: aa(:,:),bb(:) !at(:,:),bt(:),ar(:,:),br(:),ah(:,:),bh(:),
 real(r15) small,suml
 integer, allocatable:: zz(:)

 small=1d-6 !10  displacements or cos(t) are never that small

 tra=0; rot=0; hua=0
 if (itrans .ne. 0) tra_out= transl_constraints
 if (irot   .ne. 0) rot_out= rot_constraints
 if (ihuang .ne. 0) hua_out= huang_constraints

! zero lines have already been removed in the corresponding subroutines!
! call remove_zeros(tra,nc,amat(1:tra,1:nc),bmat(1:tra),  &
! &  tra_out,at(1:tra_out,1:nc),bt(1:tra_out))
! call remove_zeros(rot,nc,amat(tra+1:tra+rot,1:nc),bmat(tra+1:tra+rot),  &
! &  rot_out,ar(1:rot_out,1:nc),br(1:rot_out))
! call remove_zeros(hua,nc,amat(tra+rot+1:tra+rot+hua,1:nc),bmat(tra+rot+1:tra+rot+hua), &
! &  hua_out,ah(1:hua_out,1:nc),bh(1:hua_out))

 nl=tra_out+rot_out+hua_out+force_constraints
! allocate(aa(nl,nc),bb(nl))
! aa(1:tra_out,1:nc)=at
! bb(1:tra_out   )=bt
! aa(tra_out+1:tra_out+rot_out,1:nc)=ar
! bb(tra_out+1:tra_out+rot_out   )=br
! aa(tra_out+rot_out+1:tra_out+rot_out+hua_out,1:nc)=ah
! bb(tra_out+rot_out+1:tra_out+rot_out+hua_out   )=bh
! aa(tra_out+rot_out+hua_out+1:nl,1:nc)=aforce
! bb(tra_out+1:tra_out   )=bt

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
          write(ulog,2)' Inconsistency in line i of Amatrix which is zero!',i,suml
          write(ulog,3)' It will be removed; Corresponding line in bmatr is=',bmat(i)
          write(ulog,*)' Need to increase your range of FCs'
       endif
! record the line index for later removal
       zz(i) = 0
       nz = nz + 1  ! nz counts the number of zero lines
    endif
 enddo
 write(ulog,*)' number of lines=0 in amat are=',nz
 allocate(aa(nl-nz,nc),bb(nl-nz))

 deallocate(amat,bmat)
 allocate(amat(nl-nz,nc),bmat(nl-nz))
 amat = aa
 bmat = bb
 write(ulog,*)'REMOVE ZEROS: new number of lines in amat=',nl
 deallocate(aa,bb,zz)

 return

 j = 0
 do i=1,nl
    if (zz(i).ne.0) then
       j = j+1
       aa(j,:) = amat(i,:)
       bb(j) = bmat(i)
    elseif(i.le.transl_constraints) then
       tra_out=tra_out - 1 !nsl_constraints=transl_constraints-1
    elseif(i.le.transl_constraints+rot_constraints) then
       rot_out=rot_out - 1 !_constraints=rot_constraints-1
    elseif(i.le.transl_constraints+rot_constraints+huang_constraints) then
       hua_out = hua_out-1 !ng_constraints=huang_constraints-1
    else
       force_constraints=force_constraints-1
    endif
 enddo
 if (sum(zz).ne.nl-nz) then
    write(ulog,*)'REMOVE ZEROS: Inconsistency! sum(zz).ne. nl.nz ',sum(zz),nl-nz
    stop
 endif
 if (j.ne.nl-nz) then
    write(ulog,*)'REMOVE ZEROS: Inconsistency! j.ne. nl.nz ',j,nl-nz
    stop
 endif

 inv_constraints = tra_out+rot_out+hua_out !transl_constraints+rot_constraints+huang_constraints

 write(ulog,4)'After remove_zeros: NEW tra,rot,hua,force_constaints=',tra_out,  &
 &       rot_out,hua_out,force_constraints

 nl=j
 deallocate(amat,bmat)
 allocate(amat(nl,nc),bmat(nl))
 amat = aa
 bmat = bb
 write(ulog,*)'REMOVE ZEROS: new number of lines in amat=',nl
 deallocate(aa,bb,zz)

 2 format(a,i9,g11.4)
 3 format(a,g11.4)
 4 format(a,4(i9))

 end subroutine remove_zeros_all

!--------------------------------------------
 subroutine remove_zeros(nl,nc,matin,vectin,nlout,matout,vectout)
!! removes lines which are all zero from the matrix mat and the vector v
 use ios , only : ulog
 use params, only : verbose
 implicit none
 integer, intent(out) :: nlout
 integer, intent(in) :: nl,nc
 real(r15), dimension(nl,nc), intent(in) :: matin
 real(r15), dimension(nl)   , intent(in) :: vectin
 real(r15), allocatable, dimension(:,:), intent(out) :: matout
 real(r15), allocatable, dimension(:)  , intent(out) :: vectout
 !real(r15), dimension(nl,nc), intent(out) :: matout
 !real(r15), dimension(nl   ), intent(out) :: vectout
 integer i,j,nz
 real(r15) small,suml
 integer zz(nl)

 small=1d-6 !10  displacements or cos(t) are never that small

 zz = 1
 nz = 0
 do i=1,nl
    suml = maxval(abs(matin(i,:)))     ! sqrt(sum(amatr(i,:)*amatr(i,:))/nc)
!   if(verbose) write(*,*)' line#, suml ',i,suml
    if (suml.lt.small) then
       if (abs(vectin(i)) .lt. small) then
!         if(verbose) write(ulog,*)' line #',i,' of Amatrix and Bmatrix is zero; it will be removed!'
       else
          write(ulog,3)' Inconsistency in line i of Amatrix and Bmatrix!',i,suml,vectin(i)
          write(ulog,*)' It will be removed; increase range of FCs, or include higher orders of FCs'
       endif
! record the line index for later removal
       zz(i) = 0
       nz = nz + 1  ! nz counts the number of zero lines
    endif
 enddo
 write(ulog,*)' number of lines=0 in matin are=',nz
 write(*   ,*)' number of lines=0 in matin are=',nz
 nlout=nl-nz
 write(ulog,4)'After removing zeros: NEW # of lines=',nlout
 write(*,4)'After removing zeros: NEW # of lines=',nlout
 if (sum(zz).ne.nl-nz) then
    write(ulog,*)'REMOVE ZEROS: Inconsistency! sum(zz).ne. nl.nz ',sum(zz),nl-nz
    write(*   ,*)'REMOVE ZEROS: Inconsistency! sum(zz).ne. nl.nz ',sum(zz),nl-nz
    stop
 endif

 if(allocated(matout)) deallocate (matout)
 if(allocated(vectout)) deallocate (vectout)
 allocate(matout(nlout,nc),vectout(nlout))

 write(*   ,*)'REMOVE ZEROS: matout,vectout allocated; nlout,nc=',nlout,nc

 j = 0
 do i=1,nl
    if (zz(i).ne.0) then
       j = j+1
       matout(j,:) = matin(i,:)
       vectout(j ) = vectin(i)
    endif
 enddo

! deallocate(matin,vectin)


 3 format(a,i9,2(1x,g11.4))
 4 format(a,4(i9))

 end subroutine remove_zeros

  end module svd_stuff

!===========================================================
  module born
 use constants, only : r15
  real(r15) epsil(3,3),epsinv(3,3)
!  real(r15), allocatable:: zeu(:,:,:)
! real(r15) rho
  integer born_flag
  end module born

!===========================================================
 module linalgb
 use constants, only : r15

  interface append_array
     module procedure append_array_1d,append_array_2d
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
 subroutine symmetrize_res(mat2,res)
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
 end subroutine symmetrize_res

!---------------------------------------
 subroutine symmetrize2(n,mat2)
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(inout) :: mat2(n,n) 
 real(r15), allocatable:: mean(:,:) 

! n=size(mat2,1)
 allocate(mean(n,n))
 mean=(mat2+transpose(mat2))/2
 mat2=mean
 deallocate(mean)
 end subroutine symmetrize2


 end module linalgb


