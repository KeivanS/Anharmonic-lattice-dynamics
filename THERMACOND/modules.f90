!==============================================
! Modules here
!==============================================
 module constants
 implicit none
 real(8) :: pi=3.14159265358979323846d0
 complex(8) :: ci=cmplx(0d0,1d0)
 real(8) :: h_plank= 6.62606896d-34
 real(8) :: n_avog= 6.023d23
 real(8) :: k_b = 1.3806504d-23       ! J/K:
 real(8) :: c_light = 2.99792458d+8   ! in (m/s)
 real(8) :: hbar = 1.054571628d-34  !h_plank/2/pi  !
 real(8) :: ee = 1.60217653d-19
 real(8) :: eps0 = 8.854187817d-12
 real(8) :: me = 9.1093826d-31
 real(8) :: uma = 1.66053878d-27 ! kg
! real(8) :: cnst = sqrt(ee*1d20/uma)/pi/200/c_light ! converts sqrt(eV/A/A/uma) to cm^-1
  real(8) :: cnst= 521.1098918
! real(8) :: ryd= me*ee**4/2/pi/pi/hbar**2/16/pi/pi/eps0/eps0
  real(8) :: ryd= 27.2116  
! real(8) :: ab = hbar*hbar/me/ee/ee*4*pi*eps0
  real(8) :: ab = 0.529177
! kbe=8.617343e-05 
 end module constants
!==============================================
module params
 real(8) tolerance,margin,lmicron
 integer nconfigs,classical,ntemp,fdfiles,calc_cross,threemtrx,lamin,lamax
 integer nshells(4)   ! up to which shell to include for each rank of FC
 integer include_fc(4),nsmax  ! whether to include FCs of that rank
 real(8) rcut(4),tau0,wshift(3)
 real(8) tmin,tmax,q_cross(3)
 logical verbose
 !!mine
 real(8) alphamix
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

! interface cart_to_direct
!    module procedure cart_to_direct_aa,cart_to_direct_av,cart_to_direct_v
! end interface

! interface direct_to_cart
!    module procedure direct_to_cart_aa,direct_to_cart_av,direct_to_cart_v
! end interface

  interface reduce
     module procedure reduce_v,reduce_a
  end interface

  interface bring_to_center
     module procedure bring_to_center_a, bring_to_center_v
  end interface bring_to_center

   contains

!============================================================

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
   function v2a(w) result(a)
   type(vector), intent(in) :: w
   real(8) a(3)
   a(1) = w%x
   a(2) = w%y
   a(3) = w%z
   end function v2a
!-----------------------------------
   function myequal0(v,w) result(eq)
     use params
     real(8), intent(in)::  v,w
     logical eq !, intent(out) :: eq

        if (abs(v-w) .lt. 1d-5) then
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
!  function addition_vaa(v,w) result(add)
!    real(8) add(3) !, intent(out) ::
!    type(vector), intent(in) :: v
!    real(8), intent(in) :: w(3)
!    add(1) = v%x + w(1)
!    add(2) = v%y + w(2)
!    add(3) = v%z + w(3)
!  end function addition_vaa
!-----------------------------------
!  function addition_ava(w,v) result(add)
!    real(8) add(3) !, intent(out) ::
!    type(vector), intent(in) :: v
!    real(8), intent(in) :: w(3)
!    add(1) = v%x + w(1)
!    add(2) = v%y + w(2)
!    add(3) = v%z + w(3)
!  end function addition_ava
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
! takes direct coordinates and returns direct coordinates
 implicit none
 type(vector) v,q1,q2,q3,w
 real(8) a1,a2,a3

 a1 = v .dot. q1
 a2 = v .dot. q2
 a3 = v .dot. q3
 w%x = a1 ; w%y = a2 ; w%z = a3

 end subroutine reduce_v
!-----------------------------------------
 subroutine reduce_a(a,q1,q2,q3,w)
! takes direct coordinates and returns direct coordinates
 implicit none
 type(vector) v,q1,q2,q3,w
 real(8) a1,a2,a3,a(3)

 v%x = a(1);v%y = a(2);v%z = a(3)
 a1 = v .dot. q1
 a2 = v .dot. q2
 a3 = v .dot. q3
 w%x = a1 ; w%y = a2 ; w%z = a3

 end subroutine reduce_a
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
 type(vector) v,w,g1,g2,g3,r1,r2,r3
 real(8) a1,a2,a3,a(3),b(3)

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
 type(vector) v,g1,g2,g3,w,r1,r2,r3
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
 type(vector) v,g1,g2,g3,r1,r2,r3,w
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
 type(vector) g1,g2,g3,r1,r2,r3
 real(8) a1,a2,a3,v(3),w(3)

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

 end module geometry
!==============================================
module ios
 integer, parameter:: ulog=30,uposcar=10,uparams=11, utraj=40,ufco=20,  &
&         umap=60,umatrx=50,utimes=70,ufc1=21,ufc2=22,ufc3=23,ufc4=24,   &
&         ufc=80,ufit1=31,ufit2=32,ufit3=33,ufit4=34,uborn=12

 integer, parameter::ubs=51,udos=41,ueiv=61,ukap=71,uslf=3000, &
&         ufc0=28,ualph=27,uvel=31,umod=2000,  &
&         fbz=93,bzf=95,uibz=94,ibs=96,ucors=65,uv3=85,ugrun=81,urate=88, &
&         ucof=87,ucross=89, k_kappa=9000 

 integer, parameter:: debug=9999, self_detail=9990, iself_bs=9980, ksubset_bs_inp=9970    ! sy. 

  interface write_out
     module procedure write_outiv, write_outrv, write_outim, write_outrm &
&    , write_outi, write_outr , write_outv,write_outcv,write_outcm
  end interface

 contains

!-----------------------------
  subroutine write_outcm(unit,string,var)
  implicit none
  character*(*) string
  integer n,m,unit,l,i,j
  complex(8), dimension(:,:), intent(in) ::  var(:,:)

  l=len_trim(string)
  n=size(var(:,1)) ; m=size(var(1,:))
  write(unit,4)string(1:l)//' is='
  do i=1,n
     write(unit,5) (var(i,j),j=1,m)
  enddo
4 format(a,99(1x,g12.6))
5 format(199(1x,g12.6))
  end subroutine write_outcm
!-----------------------------
  subroutine write_outrm(unit,string,var)
  implicit none
  character*(*) string
  integer n,m,unit,l,i,j
  real(8), dimension(:,:), intent(in) ::  var(:,:)

  l=len_trim(string)
  n=size(var(:,1)) ; m=size(var(1,:))
  write(unit,4)string(1:l)//' is='
  do i=1,n
     write(unit,5) (var(i,j),j=1,m)
  enddo
4 format(a,99(1x,g12.6))
5 format(199(1x,g12.6))
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
4 format(a,99(1x,g12.5))
  end subroutine write_outrv
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
4 format(a,99(1x,g12.6))
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
4 format(a,3(1x,g12.6))
end subroutine write_outv
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

end module ios
!===========================================================
 module atoms_force_constants
! the type atom_id0 concerns the properties of the primitive unit cell
! it s atoms, masses, neighbors, shells and force constants, whereas the
! type atomic refers to atoms in the supercell, their displacements and forces.
 use geometry
 implicit none
 integer natom_prim_cell,natom_super_cell,maxshells
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
    real(8) mass
 end type
!-------------------------
 type atom_id0   ! all characteristics of an atom in the primitive central cell
    integer at_type,tau
    character name*2
    real(8) mass,charge(3,3)
    integer nshells       ! how many shells are included=actual dim(shell)
    type(vector) equilibrium_pos
    type(shell), allocatable ::  shells(:)  !0:maxshells)
!    type(tensor2r), pointer :: phi2(:)  ! refers to atom numbers neighbor of i
!    type(tensor3), pointer :: phi3(:,:)
!    type(tensor4), pointer :: phi4(:,:,:)
 end type
!-------------------------
! everything with 0 refers to the primitive cell

 integer  natom_type
 integer, allocatable:: natom(:),atom_type(:) !label_of_primitive_atoms(:),
! integer, allocatable:: ncel1(:),ncel2(:),ncel3(:)
 real(8), allocatable:: atompos0(:,:),mas(:),force(:,:,:),displ(:,:,:)
 real(8), allocatable:: forc(:,:),disp(:,:),vel(:,:),cur(:,:)
 character*2, allocatable:: atname(:)
 type (atom_id0), allocatable :: atom0(:)
 type (atomic), allocatable :: atom_sc(:),atom_shl(:)

! maximum number of shells of nearest neighbors
      integer maxneighbors
!     parameter(maxneighbors=18 )
! maximum number of atoms out to maxneighbors
      integer maxatoms,imaxat
!     parameter(maxatoms=2800 )
! op_matrix(k,j,i), matrix for the ith point operator
      double precision op_matrix(3,3,48)
      double precision op_kmatrix(3,3,48)
      integer lattpgcount
! isgopcount, number of operators in space group
      integer isgopcount
! isgop(i), point operation in ith space group operator
      integer isgop(48)
! sgfract(j,i), jth cartesian coordinate of fractional in ith space group
! operator
      double precision sgfract(3,48)
! iatomop(j,i), point operator that takes jth atom into ith atom
      integer iatomop(:,:)
! atomopfract(k,j,i), kth cartesian coordinate of fractional to accompany
! iatomop
      double precision atomopfract(:,:,:)
! natoms0, number of atoms in the primitive unit cell
      integer natoms0
! natoms, number of atoms out to maxneighbors
      integer natoms
! atompos(j,i), jth cartesian coordinate of ith atom
   !  double precision atompos(3,maxatoms)
      real(8), allocatable :: atompos(:,:)
! iatomcell(j,i), linear combination of basis vectors of the primitive lattice
! that takes us to the unit cell containing the ith atom
   !  integer iatomcell(3,maxatoms)
      integer, allocatable ::  iatomcell(:,:),iatomcell0(:)
! iatomcell0(i), identity of atom in unit cell at origin equivalent to ith
! atom
   !  integer iatomcell0(maxatoms)
! iatomneighbor(j,i), nearest neighbor shell of jth atom in primitive unit
! cell at the origin that contains the ith atom
      integer iatomneighbor(:,:)
      allocatable iatomneighbor,iatomop,atomopfract
contains

 subroutine allocate_disp_forc(n)  ! needed for the md code
    integer n
!   allocate( disp(3,n),forc(3,n),vel(3,n),cur(3,n),ncel1(n),ncel2(n),ncel3(n) )
    allocate( disp(3,n),forc(3,n),vel(3,n),cur(3,n) )
 end subroutine allocate_disp_forc

 subroutine allocate_pos(n,m)
    integer n,m
!   allocate( pos(n,m),force(3,n,m) )
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
 subroutine set_neighbor_list
! this is needed in order to know what force constants to associate
! to a pair ij, or ijk or ijkl
! if a neighbor j is known then one can get the vect(i,j) between
! their equilibrium positions
 use params
! use lattice
 use ios , only : ulog
 use geometry
 implicit none
! integer, parameter :: max=500
 integer i0,j,shel_count,counter,msort(maxatoms),l
 real(8) dist(maxatoms),d_old,d_min,rmax

! allocate( atom0(1:natoms0)%shells(maxshells) )
 rmax = 0 ; dist = 1d10
 do i0=1,natoms0
    allocate( atom0(i0)%shells(0:maxshells) )
    do j=1,natoms
       call calculate_distance(i0,j,atompos,maxatoms,dist(j) )
       if ( iatomneighbor(i0,j) .eq. nshells(2) .and. dist(j) .gt. rmax ) then
            rmax = dist(j)
       endif
    enddo
    call sort(natoms,dist,msort,maxatoms)

    shel_count = -1
    d_old = 0
    write(ulog,*)' ========================================================='
    write(ulog,*)' neighborlist of atom i0 =',i0

    jloop: do j=1,min(500,natoms)
       d_min = dist(msort(j))
       if ( (d_min .myeq. d_old) .and. (shel_count.ne.-1) ) then   ! same shell
          counter = counter + 1
       else    ! new shell
          shel_count = shel_count+1
          counter = 1
          if ( shel_count .gt. maxshells ) then
             write(ulog,*)maxshells,' shells completed'
             write(ulog,*)'shel_count=',shel_count,' exceeded it for j=',j
             exit jloop
          endif
          d_old = d_min
       endif
       if ( counter .gt. 296) then
          write(ulog,*) ' counter in neighbors exceeded 296 ', counter
          stop
       endif
       atom0(i0)%shells(shel_count)%no_of_neighbors = counter
       atom0(i0)%shells(shel_count)%rij = d_min
       atom0(i0)%shells(shel_count)%neighbors(counter)%tau = iatomcell0(msort(j))
       atom0(i0)%shells(shel_count)%neighbors(counter)%n   = iatomcell(:,msort(j))
       write(ulog,*)' ------------------------------------'
       write(ulog,5)' count, j, of tau =',j,msort(j),iatomcell0(msort(j)),   &
&                                        iatomcell(:,msort(j))
!      write(ulog,3)' cell (n1,n2,n3) ='
       write(ulog,2)' shell# , neighb# =',shel_count,counter
       write(ulog,6)' neighbor distance=',d_min
       write(ulog,4)' neighborshel,xyz =',iatomneighbor(i0,msort(j)),(atompos(l,msort(j)),l=1,3)
    enddo jloop
    do shel_count = 0 , maxshells
       write(ulog,*)'shell#, neigh#=',shel_count,atom0(i0)%shells(shel_count)%no_of_neighbors
    enddo
! also initialize atom0%equilibrium_pos

! shift the position of the first atom to the  origin
!do i=1,natoms0
!   atompos(:,i) = atompos(:,i)-atompos(:,1)
!    atom0(i0)%equilibrium_pos = atompos0(1,i0)*r01 + atompos0(2,i0)*r02  &
!&                             + atompos0(3,i0)*r03
    atom0(i0)%equilibrium_pos%x = atompos(1,i0)
    atom0(i0)%equilibrium_pos%y = atompos(2,i0)
    atom0(i0)%equilibrium_pos%z = atompos(3,i0)
!enddo
! they might not be inside the primitive cell

 enddo
 rcut(2) = rmax

2 format(a,8(2x,i4))
3 format(a,' (',i4,2(1x,i4),1x,')')
4 format(a,1x,i4,' (',f8.4,2(1x,f8.4),1x,')')
5 format(a,3(2x,i4),' (',i4,2(1x,i4),1x,')')
6 format(a,2x,g16.7)

 end subroutine set_neighbor_list

 end module atoms_force_constants
!===========================================================
 module lattice
 use geometry
 use constants
 implicit none
 type(vector) r1,r2,r3,g1,g2,g3     ! translation vectors of the supercell
 type(vector) r01,r02,r03,g01,g02,g03  ! tr vect of prim cell and its recip spce
 type(vector) r1conv,r2conv,r3conv,g1conv,g2conv,g3conv
 real(8) volume_r,volume_g,lattice_parameter,latticeparameters(6),primitivelattice(3,3)
 real(8) conv_to_cart(3,3),conv_to_prim(3,3), cart_to_prim(3,3), prim_to_cart(3,3),prim_to_conv(3,3)
 real(8) box(3),boxg(3)
 real(8) r0g(3,3)
 integer n1min,n2min,n3min,n1max,n2max,n3max !,NC(3),NF(3)

  contains 
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

!function bring_to_prim_cell_cav(a) result(w)
! takes cartesian coordinates and returns cart coordinates within the primcell
!use geometry
!implicit none
!type(vector) v,w, bring_to_prim_cell_cv
!real(8) a(3)

!v%x = a(1) ; v%y = a(2) ; v%z = a(3)
!w =  bring_to_prim_cell_cv(v)

!end function bring_to_prim_cell_cav
!-----------------------------------------
!function bring_to_super_cell_cav(a) result(w)
! takes cart coordinates and returns cart coordinates within the supercell
!use geometry
!implicit none
!type(vector) v,w, bring_to_super_cell_cv
!real(8) a(3)

!v%x = a(1) ; v%y = a(2) ; v%z = a(3)
!w = bring_to_super_cell_cv(v)

!end function bring_to_super_cell_cav
!-----------------------------------------
!function bring_to_prim_cell_caa(a) result(b)
! takes cartesian coordinates and returns cart coordinates within the primcell
!use geometry
!implicit none
!type(vector) v,w, bring_to_prim_cell_cv
!real(8) a(3),b(3)

!v%x = a(1) ; v%y = a(2) ; v%z = a(3)
!w =  bring_to_prim_cell_cv(v)
!b(1) = w%x ; b(2) = w%y ; b(3) = w%z

!end function bring_to_prim_cell_caa
!-----------------------------------------
!function bring_to_super_cell_caa(a) result(b)
! takes cart coordinates and returns cart coordinates within the supercell
!use geometry
!implicit none
!type(vector) v,w, bring_to_super_cell_cv
!real(8) a(3),b(3)

!v%x = a(1) ; v%y = a(2) ; v%z = a(3)
!w = bring_to_super_cell_cv(v)
!b(1) = w%x ; b(2) = w%y ; b(3) = w%z

!end function bring_to_super_cell_caa
!-----------------------------------------
!function bring_to_prim_cell_cv(v) result(w)
! takes cartesian coordinates and returns cart coordinates within the primcell
!use geometry
!implicit none
!type(vector) v,w
!real(8) a1,a2,a3

! get direct coordinates first
!a1 = v .dot. g01
!a2 = v .dot. g02
!a3 = v .dot. g03
! bring into cell: between 0 and cell
!a1 = a1 - floor(a1)
!a2 = a2 - floor(a2)
!a3 = a3 - floor(a3)
! convert to cartesian coordinates
!w = a1*r01+a2*r02+a3*r03

!end function bring_to_prim_cell_cv
!-----------------------------------------
!function bring_to_super_cell_cv(v) result(w)
! takes cart coordinates and returns cart coordinates within the supercell
!use geometry
!implicit none
!type(vector) v,w
!real(8) a1,a2,a3

! get direct coordinates first
!a1 = v .dot. g1
!a2 = v .dot. g2
!a3 = v .dot. g3
! bring into cell: between 0 and cell
!a1 = a1 - floor(a1)
!a2 = a2 - floor(a2)
!a3 = a3 - floor(a3)
! convert to cartesian coordinates
!w = a1*r1+a2*r2+a3*r3

!end function bring_to_super_cell_cv
!-----------------------------------
 subroutine cart_to_direct_v(v,w)
! takes cart coordinates and returns cart coordinates within the supercell
 use geometry
 implicit none
 type(vector) v,w

 w%x = v .dot. g1
 w%y = v .dot. g2
 w%z = v .dot. g3

 end subroutine cart_to_direct_v
!-----------------------------------
 subroutine cart_to_direct_av(a,w)
! takes cart coordinates and returns cart coordinates within the supercell
 use geometry
 implicit none
 real(8) a(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call cart_to_direct_v(v,w)

 end subroutine cart_to_direct_av
!-----------------------------------
 subroutine cart_to_direct_aa(a,b)
! takes cart coordinates and returns cart coordinates within the supercell
 use geometry
 implicit none
 real(8) a(3),b(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call cart_to_direct_v(v,w)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end subroutine cart_to_direct_aa
!-----------------------------------
 subroutine direct_to_cart_v(v,w)
! takes direct coordinates and returns cart coordinates
 use geometry
 implicit none
 type(vector) v,w

 w%x = v%x*r1%x + v%y*r2%x + v%z*r3%x
 w%y = v%x*r1%y + v%y*r2%y + v%z*r3%y
 w%z = v%x*r1%z + v%y*r2%z + v%z*r3%z

 end subroutine direct_to_cart_v
!-----------------------------------
 subroutine direct_to_cart_av(a,w)
! takes direct coordinates and returns cart coordinates
 use geometry
 implicit none
 real(8) a(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call direct_to_cart_v(v,w)

 end subroutine direct_to_cart_av
!-----------------------------------
 subroutine direct_to_cart_aa(a,b)
! takes direct coordinates and returns cart coordinates
 use geometry
 implicit none
 real(8) a(3),b(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call direct_to_cart_v(v,w)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end subroutine direct_to_cart_aa
!------------------------------------
 subroutine dir2cart_g(q,k)
! takes q in direct coordinates and outputs k in cartesian
 use geometry
 real(8) q(3),k(3)

 k(:) = q(1)*g1 + q(2)*g2 + q(3)*g3

 end subroutine dir2cart_g
!------------------------------------

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
    integer, allocatable:: nt(:),ntind(:),igr(:)
    character*1, allocatable:: err(:)
    integer ngr
 end type fulldmatrix

 integer maxrank,maxgroups,itrans,irot,ihuang,enforce_inv
 integer nterms(4),maxterms(4),maxtermsindep(4),ngroups(4)
 integer ndindp(4),ndfull(4)
 parameter(maxrank=4)
 integer, allocatable :: nterm(:),ntermsindep(:)
 integer, allocatable :: iatomtermindep(:,:,:),ixyztermindep(:,:,:)
 integer, allocatable :: iatmtrm(:,:,:),ixyztrm(:,:,:)
 real(8), allocatable :: mapmat(:,:,:)
 real(8) svdcut,radius(4)
 integer, allocatable:: iatomterm_1(:,:),ixyzterm_1(:,:),igroup_1(:),map_1(:)
 integer, allocatable:: iatomterm_2(:,:),ixyzterm_2(:,:),igroup_2(:),map_2(:)
 integer, allocatable:: iatomterm_3(:,:),ixyzterm_3(:,:),igroup_3(:),map_3(:)
 integer, allocatable:: iatomterm_4(:,:),ixyzterm_4(:,:),igroup_4(:),map_4(:)
 character*1, allocatable:: err_1(:),err_2(:),err_3(:),err_4(:)
 type(fulldmatrix) map(4)
 real(8), allocatable:: ampterm_1(:),fcs_1(:)
 real(8), allocatable:: ampterm_2(:),fcs_2(:),grun_fc(:)
 real(8), allocatable:: ampterm_3(:),fcs_3(:)
 real(8), allocatable:: ampterm_4(:),fcs_4(:)
 real(8), allocatable:: amat(:,:),bmat(:),sigma(:),fcs(:)
 real(8), allocatable:: a11ia12(:,:),fc1(:)
 real(8), allocatable:: arot(:,:),brot(:),atransl(:,:),ahuang(:,:),aforce(:,:),bforce(:)
 integer inv_constraints,force_constraints,dim_al,dim_ac,n_indep,newdim_al    &
&        ,transl_constraints, rot_constraints, huang_constraints,ngr
contains

 subroutine set_maxterms
   maxterms(1)=15
   maxterms(2)=1200 
   maxterms(3)=10
   maxterms(4)=10 
   maxtermsindep(1)=5
   maxtermsindep(2)=30
   maxtermsindep(3)=10
   maxtermsindep(4)=10
   maxgroups=30
 end subroutine set_maxterms

  end module svd_stuff
!============================================================
  module born
  real(8) epsil(3,3)
  real(8), allocatable:: zeu(:,:,:)
  real(8) rho 
  integer born_flag
  end module born

!===================================================================
module tetrahedron

        implicit none
        type point
                real(8) :: w,c,d
                integer :: i
        end type point

        type tetra
                type(point) :: p(4)
        end type tetra

        real(8), allocatable :: calct(:,:),doswt(:,:)
        type(tetra), allocatable :: tet(:)

!        integer num1_tet, num2_tet, num3_tet, num4_tet, num5_tet

contains

subroutine allocate_tetra(n1,n2,n3,wm,nd) !num1_tet,num2_tet,num3_tet,num4_tet,num5_tet)
 implicit none
 integer, intent(in):: n1,n2,n3,wm,nd ! num1_tet, num2_tet, num3_tet, num4_tet, num5_tet


 allocate(calct(wm,nd)) !num4_tet,num5_tet))
 allocate(doswt(wm,nd)) !num4_tet,num5_tet))
 allocate(tet(n1*n2*n3*6))  !(num1_tet)*(num2_tet)*(num3_tet)*6))
end subroutine allocate_tetra

subroutine deallocate_tetra
 deallocate(calct,doswt,tet)
end subroutine deallocate_tetra

end module tetrahedron



!===========================================================
module pre_matrix_elt
complex(8), allocatable :: eiqr2(:,:), eivec2(:,:,:,:)

contains

subroutine allocate_eiqr2(num1,num2)
 allocate(eiqr2(num1,num2))
end subroutine allocate_eiqr2
subroutine allocate_eivec2(num1,num3,num4)
 allocate(eivec2(num1,3,num3,num4))
end subroutine allocate_eivec2

subroutine deallocate_eiqr2
 deallocate(eiqr2)
end subroutine deallocate_eiqr2
subroutine deallocate_eivec2
 deallocate(eivec2)
end subroutine deallocate_eivec2

end module pre_matrix_elt


!===========================================================
module mod_ksubset2
integer num_pfiles
integer, allocatable :: ksubset2(:,:)

contains

 subroutine allocate_ksubset2(num)
   allocate(ksubset2(num,2))
 end subroutine allocate_ksubset2
 subroutine deallocate_ksubset2
   deallocate(ksubset2)
 end subroutine deallocate_ksubset2

end module mod_ksubset2


!==========================================================
 module phi3_sy
 complex(8), allocatable :: v33_md(:,:,:,:,:)
 real(8), allocatable :: v33sq1(:,:,:,:,:)
! v3 : (nk1,nk2,la1,la2,la3)
! integer nv3,readv3,writev3, iter, split, calc_kappa ! sy added iter,split, calc_kappa
! real(8) v3_threshold,const33

 contains

  subroutine allocate_v33_sy(n1,n2,l1,l2,l3)
    allocate(v33_md(n1,n2,l1,l2,l3) )
  end subroutine allocate_v33_sy
  subroutine deallocate_v33_sy
    deallocate(v33_md)
  end subroutine deallocate_v33_sy
  subroutine allocate_v33_sq1(n1,n2,l1,l2,l3)
    allocate(v33sq1(n1,n2,l1,l2,l3))
  end subroutine allocate_V33_sq1
  subroutine deallocate_v33_sq1
    deallocate(v33sq1)
  end subroutine
  subroutine allocate_phi3_sy(nk,nl,nw,ni)
    integer nk,nw,nl,ni
    allocate(v33_md(ni,nk,nl,nl,nl))
  end subroutine allocate_phi3_sy

 end module phi3_sy




!===========================================================
module exactBTE2

integer ncol
real(8), allocatable:: F_RTA(:,:,:), F1(:,:,:), Fs1(:,:,:) ,F2(:,:,:), Qvalue(:,:), iselfenergy(:,:), tau(:,:,:), error(:,:,:), diff(:,:,:)
real(8), allocatable:: F1_old(:,:,:)
real(8), allocatable:: Qvalue_N(:,:), Qvalue_U(:,:), tauinv_N(:,:), tauinv_U(:,:), tauinv_tot(:,:), tauinv_eff(:,:)
real(8), allocatable:: diff_kap(:,:,:,:),kappa(:,:,:), kappa_k(:,:,:,:), kappa_RTA(:,:,:), kappa_k_RTA(:,:,:,:)
! F1,F2: (kp1,la1,xyz)    , Qvalue,Avalue,iselfenergy: (kp,la)
real(8), allocatable:: P1(:,:,:,:,:), P2(:,:,:,:,:),Piso(:,:),Pisos(:,:),giso(:)
! P1,P2: (kp1,kp2,la1,la2,la3)   don't need kp3 since it is determined by momentum conservation.
real(8), allocatable:: frequency(:,:),dist(:,:)   ! BE distribution function
!real(8), allocatable:: v33sq(:,:,:,:,:)

real(8),allocatable::DirectSolution(:),sig(:),col_matrix(:,:),RHS(:),sy_matrix(:,:),coll(:,:),F_MRHS(:),Collision_Matrix(:,:),Collision_Matrixmpi(:,:),FFF(:,:,:),kappa_new(:,:,:),kappa_old(:,:,:),e(:,:,:)!,giso(:)
real(8) errorr,ermax,svdcut

contains

!subroutine allocate_iter (ni,nk,ndn)
!implicit none
!integer, intent(in) :: ni,nk,ndn
!allocate (F_RTA(ni,ndn,3),Qvalue(ni,ndn),iselfenergy(ni,ndn),tau(ni,ndn,3)) 
!allocate (F1(nk,ndn,3),F2(nk,ndn,3),F1_old(nk,ndn,3),error(nk,ndn,3),diff(nk,ndn,3)) 
!allocate (diff_kap(nk,ndn,3,3),kappa(ndn,3,3),kappa_RTA(ndn,3,3))    ! 3 for xyz
!allocate (kappa_k(ni,ndn,3,3),kappa_k_RTA(ni,ndn,3,3))    ! 3 for xyz
!allocate (dist(nk,ndn),frequency(nk,ndn))
!allocate (Qvalue_N(ni,ndn),Qvalue_U(ni,ndn),tauinv_N(ni,ndn),tauinv_U(ni,ndn),tauinv_tot(ni,ndn),tauinv_eff(ni,ndn))
!end subroutine allocate_iter

subroutine allocate_iter (ni,nk,ndn)
implicit none
integer, intent(in) :: ni,nk,ndn
!allocate (F_RTA(ni,ndn,3),Qvalue(ni,ndn),iselfenergy(ni,ndn),tau(ni,ndn,3))
!allocate (F1(nk,ndn,3),F2(nk,ndn,3),F1_old(nk,ndn,3),error(nk,ndn,3),diff(nk,ndn,3))
!allocate (diff_kap(nk,ndn,3,3),kappa(ndn,3,3),kappa_RTA(ndn,3,3))    ! 3 for xyz
!allocate (kappa_k(ni,ndn,3,3),kappa_k_RTA(ni,ndn,3,3))    ! 3 for xyz
allocate (dist(nk,ndn),frequency(nk,ndn))
allocate (Qvalue_N(ni,ndn),Qvalue_U(ni,ndn),tauinv_N(ni,ndn),tauinv_U(ni,ndn),tauinv_tot(ni,ndn),tauinv_eff(ni,ndn))
allocate (Qvalue(ni,ndn),kappa(ndn,3,3))
end subroutine allocate_iter

subroutine allocate_FGR (n1,n2,l1,l2,l3)
implicit none
integer, intent(in) :: n1,n2,l1,l2,l3
allocate (P1(n1,n2,l1,l2,l3),P2(n1,n2,l1,l2,l3),Piso(n1,l1))      
end subroutine allocate_FGR

subroutine deallocate_iter
!deallocate(F_RTA,F1,F2,Qvalue,iselfenergy,tau,error,diff,F1_old)
!deallocate(diff_kap,kappa,kappa_k,kappa_RTA,kappa_k_RTA)
!deallocate(dist,frequency)
!deallocate(Qvalue_N,Qvalue_U,tauinv_N,tauinv_U,tauinv_tot,tauinv_eff)
end subroutine deallocate_iter

subroutine deallocate_FGR
deallocate(P1,P2)
end subroutine deallocate_FGR

end module exactBTE2


!==========================================================
 module om_dos
 integer wmesh,ndyn2
 real(8) wmax,width,etaz
 real(8), allocatable :: dos(:,:),om(:)

   contains

   subroutine set_omdos(la,mesh)
   implicit none
   integer la,mesh,i
   allocate(om(mesh),dos(la+1,mesh))
   ndyn2 = la  ! this is a copy of ndyn
   do i=1,mesh
      om(i) = wmax *(0.0001 + (i-1)/(mesh-1.))
   enddo
   end subroutine set_omdos

   subroutine write_dos
   use ios
   implicit none
   integer i,la
   real(8) sumdos

   open(udos,file='dos.dat',status='unknown')

   write(udos,*)'# om(i),i,integrated_dos,total_dos,(dos(n+1-la,i),la=1,n))'
   sumdos=0
   do i=1,wmesh-1
      sumdos=sumdos+dos(ndyn2+1,i)*wmax/wmesh
      write(udos,3)om(i),i,sumdos,(dos(ndyn2+2-la,i),la=1,min(10,ndyn2+1))
   enddo
   sumdos=sumdos-dos(ndyn2+1,1)*wmax/wmesh/2d0  ! 1/2 the weight for last point
   write(udos,3)om(wmesh),wmesh,sumdos,(dos(ndyn2+2-la,wmesh),la=1,min(10,ndyn2+1))

   close(udos)
 
3  format(g11.5,2x,i5,99(1x,g10.4))
   end subroutine write_dos

 end module om_dos
!==========================================================
 module eigen
 integer ndyn
! eigenval for the coarse mesh (nkc,kpc),
! eigenval_c for (nibz_coarse,kibz), eigenval_f for (nbz_fine,kibzf)

 real(8), allocatable :: veloc_f(:,:,:),veloc(:,:,:)
 real(8), allocatable :: eigenval_f(:,:),eigenval_bs(:,:),eigenval(:,:)
! sqr of phonon freqs in eV/A/A/m0
 complex(8), allocatable :: eigenvec_f(:,:,:),eigenvec_bs(:,:,:),eigenvec(:,:,:),grun(:,:),grun_bs(:,:)

    contains

    subroutine allocate_eig(nb,ni,nk) ! nb for band, nk for coarse mesh in FBZ
    integer nb,nk,ni
      allocate( eigenval(nb,ni),eigenvec(nb,nb,ni),grun(nb,ni)  ,veloc(3,nb,nk) )
    end subroutine allocate_eig
!---------------------------------
    subroutine allocate_eig_bs(nb,nk,nv) ! nb for band, nk for band structure mesh
    integer nb,nk,nv
     if (nv.ne.0) then
      allocate( eigenval_bs(nb,nk),eigenvec_bs(nb,nv,nk), grun_bs(nb,nk) ,veloc(3,nb,nk) )
     else
      allocate( eigenval_bs(nb,nk),eigenvec_bs(1,1,1), grun_bs(nb,nk)) 
     endif
    end subroutine allocate_eig_bs
!---------------------------------
    subroutine deallocate_eig_bs ! nb for band, nk for band structure mesh
      if(allocated(veloc)) deallocate(veloc)
      deallocate( eigenval_bs, eigenvec_bs, grun_bs ) 
    end subroutine deallocate_eig_bs
!---------------------------------
!   subroutine allocate_eig_f(nb,nk) ! nb for band, nk for mesh in IBZ
!   integer nb,nk
!     allocate( eigenval_f(nb,nk),eigenvec_f(nb,nb,nk) ) ! ,veloc_f(3,nb,nk) )
!   end subroutine allocate_eig_f
!---------------------------------
    subroutine deallocate_eig ! nb for band, nk for coarse mesh in FBZ
      deallocate( eigenval,eigenvec,grun,veloc )
    end subroutine deallocate_eig

 end module eigen
!==========================================================
 module phi3
 complex(8), allocatable :: selfN(:,:,:) ,selfU(:,:,:) !v3(:,:,:,:,:),
 real(8), allocatable :: v33sq(:),v33sq_5(:,:,:,:,:),v33s8(:,:,:,:,:),P1smpi(:,:,:,:,:),P2smpi(:,:,:,:,:),p1sm_mpi(:,:,:,:,:),p2sm_mpi(:,:,:,:,:)
 real(8), allocatable :: isempi(:,:),tsempi(:,:),ise_mpi(:,:),tse_mpi(:,:)
 integer, allocatable:: nq1(:),nq2(:),nq3(:),la1(:),la2(:),la3(:)
 integer nv3,readv3,writev3, iter, split, calc_kappa, max_iter, job, iso ! sy added iter,split, calk, max_iter
 real(8) v3_threshold,const33, conv_error , conv_max_error, conv_max_diff, conv_diff, conv_diff_kap, conv_max_diff_kap,conv_iter,update_mix   ! sy added conv_norm, conv_max, update
 integer ksub_size
 character(99) v3path

 contains

  subroutine allocate_v33sq(n)
    allocate(v33sq(n),nq1(n),nq2(n),nq3(n),la1(n),la2(n),la3(n) )
  end subroutine allocate_v33sq
  subroutine deallocate_v33sq
    deallocate(v33sq,nq1,nq2,nq3,la1,la2,la3)
  end subroutine deallocate_v33sq
  subroutine allocate_phi3(nk,nl,nw,ni)
  integer nk,nw,nl,ni
!   allocate(delta3(nk,nl,nw),gamma3(nk,nl,nw),v3(nk,nk,nl,nl,nl))
   allocate(selfN(ni,nl,nw),selfU(ni,nl,nw))
!    allocate(v3(ni,nk,nl,nl,nl))
  end subroutine allocate_phi3

 end module phi3
!==========================================================
 module kpoints
! this includes 4 kpoint meshes: one for the band structure: kp_bs(3,nkp_bs)
! one shifted MP-coarse mesh around gamma for interpolation : kpc(3,nkc)
! a mapping of this coarse mesh into the FBZ: kpc(3,mapbz(1:nbz)) 
! one fine grid in the irreducible FBZ for kappa/dos sums: kibz(3,nibz)
 use lattice   ! it is also used for all included subroutines
 use ios
 use constants
 implicit none
 integer nbs1,nbs2,nbs3,nkp_bs   ! for the band structure
 integer nc(3),nkc,nbz,nibz_coarse,nkcc  ! coarse mesh,its FBZ folded mesh,IBZ_coarse mesh
 integer nib1,nib2,nib3,nbz_fine  ! for the fine mesh in the Irreducible BZ
 integer nshl,nbkmax,npos           ! max no of neighbors=20
 integer, allocatable :: nb(:,:),nbk(:),mapbz(:) ,mapibz(:),mappos(:),mapinv(:)  ! nb_list of k
 real(8) shft(3),deltak,kcutnbf,kcutnbc   ! kcutnb for nbrlst and interpolation
 real(8), allocatable :: kpc(:,:),wk(:),wkf(:),kp_bs(:,:),gg(:,:),dk_bs(:)
 real(8), allocatable :: kbzf(:,:) ,kibz(:,:),wibz(:),kbz(:,:)
 character(LEN=3), allocatable :: kname_bs(:) 
 integer nibz
 logical dos_bs_only

 contains

  subroutine allocate_kp_bs(n)
  integer n
    allocate (kp_bs(3,n),dk_bs(n))
  end subroutine allocate_kp_bs
!-------------------------------------------
  subroutine deallocate_kp_bs
    deallocate (kp_bs,dk_bs) 
  end subroutine deallocate_kp_bs
!-------------------------------------------
 !-------------------------------------------
   subroutine make_kp_bs
 !! this subroutine sets up the kpoints for the band structure from the
 !! list written in the file kpbs.params (given in direct coordinates of the CONVENTIONAL cell)
   use ios
   use geometry
   use lattice !, only : g01,g02,g03
   use constants, only : pi
     implicit none
     integer :: i,j,nk,uio,nkdir,ubs1,units ,ndir
     real(kind=8) lk,q(3),k_conv(3),k_prim(3),k_cart(3)
     real(kind=8), allocatable :: ki(:,:),kf(:,:) ,kext_bs(:,:)
     character(LEN=6), allocatable :: dk(:)

     write(*,*)'entering MAKE_KP_BS'
     uio = 67
     ubs1 = 68
     open(uio,file='kpbs.params',status='old')
     open(ubs1,file='redtocart')

     read(uio,*) units  ! if 0 conventional, else primitive
     write(*,*)'reading ',units,nkdir,ndir
     read(uio,*) nkdir  ! number of kpoints along each direction
     write(*,*)'reading ',units,nkdir,ndir
     read(uio,*) ndir  ! number of directions for the band structure
     write(ubs1,*)'reading units,nkdir,ndir=',units,nkdir,ndir
     write(ubs1,*)'# k_name , k in primitive, k in gconv units, k in cartesian '
     nkp_bs = ndir*nkdir    ! +1 is for including the last point
     write(*,*)'reading nkp_bs= ',nkp_bs
     allocate(kp_bs(3,nkp_bs),ki(3,ndir),kf(3,ndir),dk_bs(nkp_bs),kname_bs(ndir+1),kext_bs(3,ndir+1),dk(ndir+1))

     do i=1,ndir+1
        read(uio,*) kname_bs(i)(2:2),kext_bs(:,i)   ! in direct coordinates of conventional
        kname_bs(i)(1:1)='"'
        kname_bs(i)(3:3)='"'
        write(*,'(i5,1x,a3,9(1x,f9.4))')i,kname_bs(i),kext_bs(:,i)
     enddo

     do i=1,ndir
        if(units.eq.0) then ! reading in units of conventional
           k_conv=kext_bs(:,i)
 ! convert to units of primitive (k_prim) and cartesian (k_cart)
           k_prim=matmul(transpose(prim_to_conv),kext_bs(:,i))
           k_cart=2*pi*(matmul(cart_to_prim,k_prim))
 ! this q is cartesian units
           q=kext_bs(1,i  )*g1conv+kext_bs(2,i  )*g2conv+kext_bs(3,i  )*g3conv
           if(length(q-k_cart).gt.1d-5) write(*,9)'MAKE_KP_BS: Error,units,q,k_cart',units,q,k_cart
           ki(:,i)=q       ! ki should also be in cartesian
           q=kext_bs(1,i+1)*g1conv+kext_bs(2,i+1)*g2conv+kext_bs(3,i+1)*g3conv
           kf(:,i)=q

 ! conventional reduced and primitive
           write(ubs1,9)kname_bs(i),i, k_prim,k_conv,k_cart

        else  ! reading in units of primitive
           q=kext_bs(:,i)
           ki(:,i) = q(1)*g01+q(2)*g02+q(3)*g03
           q=kext_bs(:,i+1)
           kf(:,i) = q(1)*g01+q(2)*g02+q(3)*g03
 ! convert to units of conventional
           k_prim=kext_bs(:,i)
           k_cart=2*pi*matmul(cart_to_prim,k_prim)
           k_conv=matmul(transpose(conv_to_cart),k_cart)/2/pi

 ! conventional reduced and primitive
           write(ubs1,9)kname_bs(i),i, kext_bs(:,i),k_conv,k_cart

        endif
     enddo

     i=ndir+1
     if(units.eq.0) then ! reading in units of conventional
           k_cart=kext_bs(1,i)*g1conv+kext_bs(2,i)*g2conv+kext_bs(3,i)*g3conv
           k_prim=matmul(transpose(prim_to_conv),kext_bs(:,i))
           write(ubs1,9)kname_bs(i),i, k_prim,kext_bs(:,i),k_cart
     else
           k_cart=kext_bs(1,i)*g01+kext_bs(2,i)*g02+kext_bs(3,i)*g03
           k_conv=matmul(transpose(conv_to_cart),k_cart)/2/pi
           write(ubs1,9)kname_bs(i),i, kext_bs(:,i),k_conv,k_cart
     endif

     close(uio)
     close(ubs1)

   9 format(a,i4,9(1x,f9.5))

   write(ulog,*)' Kpoints for band structure generated from kpbs.params'
   write(*,*)' Kpoints for band structure generated from kpbs.params'

   write(ulog,*) g01
   write(ulog,*) length(g01)
   write(ulog,*) g02
   write(ulog,*) length(g02)
   write(ulog,*) g03
   write(ulog,*) length(g03)

  open(ubs1,file='KPOINT.BS',status='unknown')
  write(ubs1,*)'# k in cartesian and in direct units of primitive and conventional reciprocal lattice'
  open(88,file='KTICS.BS',status='unknown')
     nk = 0
 !    dk_bs(1) = 0
 ! now for each direction set up the coordinates of the kpoints
     kext_bs(1,1)=0    ! here using kext_bs(1,:) as a dummy variable
     do i = 1,ndir !-1
        write(ubs1,7) i, matmul(prim_to_cart,ki(:,i)),ki(:,i), matmul(transpose(conv_to_prim),ki(:,i))
        q(:) = kf(:,i)-ki(:,i)  ! length of each section
        lk = length(q) !/length(g01)  ! length of each section
        do j = 1,nkdir
           nk = nk+1
          if ( nk.ge.2) then
             if (j.ge.2) then
               dk_bs(nk) = dk_bs(nk-1) + lk/(nkdir-1+1d-8)
             else
               dk_bs(nk) = dk_bs(nk-1)
             endif
           else  ! nk = 1
             dk_bs(nk) = 0 !dk_bs(nk-1)
           endif
           kp_bs(:,nk) = ki(:,i) + (j-1)*q(:)/(nkdir-1+1d-8)
 
    write(*,*)' dk_bs =', nk,dk_bs(nk)

        enddo
        kext_bs(1,i+1)=real(dk_bs(nk))   ! here using kext_bs(1,:) as a dummy variable
     enddo
 !   kext_bs(1,ndir+1)=dk_bs(nk)

     do i=1,ndir
        write(dk(i),'(f6.3)') kext_bs(1,i)
     enddo
     write(dk(i),'(f6.3)') kext_bs(1,ndir+1)-0.002

 !   write(88,*)'set xtics ( ',(kname_bs(i),dk(i),",",i=1,ndir+1),' )'
     write(88,*)'set xtics ( ',(kname_bs(i),dk(i),",",i=1,ndir),kname_bs(ndir+1),dk(ndir+1),' )'
     close(ubs1)
     close(88)

     deallocate(ki,kf,kname_bs)
 3   format(9(2x,f10.4))
 7   format(i4,99(2x,f10.4))
 5   format(a,i4,9(2x,f10.4))
 !4   format(a,(ndir+1)(a,1x,f6.4,a))

     write(ulog,*)' MAKE_KP_BS: Done! with ',nkp_bs,' points'
     write(*,*)' MAKE_KP_BS: Done! with ',nkp_bs,' points'

  end subroutine make_kp_bs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****  end subroutine make_kp_bs (erased)
!-------------------------------------------
  subroutine make_kp_bs_old
! this subroutine sets up the kpoints for the band structure from the 
! list written in the file kpbs.in (given in direct coordinates)
    implicit none
    integer :: nkdir,ndir,i,j,nk,uio
    real(8) ep,lk,q(3),q2(3)
    real(8), allocatable :: ki(:,:),kf(:,:)

    uio = 167
    open(uio,file='kpbs.in',status='old')

    read(uio,*) nkdir  ! number of kpoints along each direction 
!   if (nkdir.eq.1) then
!       write(ulog,*)'make_kp_bs2: nkdir must be equal or larger than 2'
!       stop
!   endif
! includes firstpoint not the last, except for the last direction
    read(uio,*) ndir  ! number of directions for the band structure
    nkp_bs = ndir*nkdir    ! +1 is for including the last point
    allocate(kp_bs(3,nkp_bs),ki(3,ndir),kf(3,ndir),dk_bs(nkp_bs))

    do i=1,ndir
!       read(uio,*) ki(:,i),kf(:,i)   ! in direct coordinates
       read(uio,*) q(:),q2(:)   ! in direct coordinates

!       call cart_to_direct_aa(q ,ki(:,i))
!       call cart_to_direct_aa(q2,kf(:,i))
       call dir2cart_g(q ,ki(:,i))
       call dir2cart_g(q2,kf(:,i))
    enddo 

    close(uio)

    nk = 0
!    dk_bs(1) = 0
! now for each direction set up the coordinates of the kpoints
    do i = 1,ndir
       q(:) = kf(:,i)-ki(:,i)  ! length of each section
       lk = length(q)/length(g1)  ! length of each section
       do j = 1,nkdir
          nk = nk+1
          if ( nk.ge.2) then
            if (j.ge.2) then
              dk_bs(nk) = dk_bs(nk-1) + lk/(nkdir-1+1d-8)
            else
              dk_bs(nk) = dk_bs(nk-1) 
            endif
          else  ! nk = 1
            dk_bs(nk) = 0 !dk_bs(nk-1)
          endif
          kp_bs(:,nk) = ki(:,i) + (j-1)*q(:)/(nkdir-1+1d-8)
          write(ibs,3) kp_bs(:,nk)
       enddo
    enddo

    deallocate(ki,kf)
3   format(9(2x,f10.4))

    write(ulog,*)' MAKE_KP_BS: Done! with ',nkp_bs,' points'

  end subroutine make_kp_bs_old
!-------------------------------------------
  subroutine make_kp_reg(nc,gg1,gg2,gg3,shft,kpt,wkt)
! generate a regular mesh from 0 to g_i, with eventual shift 
    use geometry
    use params
    implicit none
    integer, intent(in) :: nc(3)
    real(8), intent(in) :: shft(3) 
    type(vector), intent(in) :: gg1,gg2,gg3
    type(vector) rr1,rr2,rr3
    real(8) vg1(3),vg2(3),vg3(3),sh(3)  
    integer :: i,j,k,nk
    real(8) dkc(nc(1)*nc(2)*nc(3))
    real(8) q(3),ep(3),sx,sy,sz,kpt(3,nc(1)*nc(2)*nc(3)),wkt(nc(1)*nc(2)*nc(3))
    real(8) a1,a2,a3

    open(126,file='KPOINT.MP',status='unknown')
    write(126,*)'# i ,kp(:,i), wk(i) , kp_red'
    ep = 0d0  ! shift by a small amount so that only one boundary K remains in the FBZ after folding
    vg1=v2a(gg1)
    vg2=v2a(gg2)
    vg3=v2a(gg3)
    sh = (-0.5d0)*(gg1+gg2+gg3)
    nk = 0
    do i = 1,nc(1)
    do j = 1,nc(2)
    do k = 1,nc(3)
       nk = nk+1
!      q = ((i-1+shft(1))/nc(1))*g1 + ((j-1+shft(2))/nc(2))*g2 + ((k-1+shft(3))/nc(3))*g3 + ep
       q = (float(i-1)/nc(1))*vg1 + (float(j-1)/nc(2))*vg2 + (float(k-1)/nc(3))*vg3 + ep
       kpt(:,nk) = q(:) 
       wkt(nk)=1d0/nkc
       if (nk.eq.1) then
          dkc(nk)=length(q)
       else
          dkc(nk)=length(q-kpt(:,nk-1))
       endif
      
    enddo
    enddo
    enddo

    write(ulog,*)'KP_REG: Number of regular kpoints generated is=',nk

    if (mod(nc(1),2).eq.0 .and. mod(nc(2),2).eq.0 .and. mod(nc(3),2).eq.0 ) then
       do i=1,nk
          kpt(:,i) = kpt(:,i)+sh(:)
       enddo
    endif

    call make_reciprocal_lattice(gg1,gg2,gg3,rr1,rr2,rr3)

    do i=1,nk
       q=kpt(:,i)
       a1 = (v2a(rr1) .dot. q)
       a2 = (v2a(rr2) .dot. q)
       a3 = (v2a(rr3) .dot. q)
       write(126,2)i,q,wkt(i),a1,a2,a3
    enddo

 2  format(i7,2x,3(1x,f12.5),5x,f9.5)
 3  format(3(i3),2x,i6,2x,3(1x,f12.5),5x,f9.5)
    close(126)

  end subroutine make_kp_reg
!-------------------------------------------
  subroutine make_kp_FBZ(nk)
! generate a regular cubic mesh of kpoints with eventual shift, within FBZ
    use geometry
    implicit none
    integer :: i,j,k,nx,ny,nz,nktot,nk
    logical inside
    real(8) q(3),qdotg,ep,gmaxx,gmaxy,gmaxz
    real(8), allocatable :: kpt(:,:)

    gmaxx = max(abs(g1%x),abs(g2%x),abs(g3%x))
    gmaxy = max(abs(g1%y),abs(g2%y),abs(g3%y))
    gmaxz = max(abs(g1%z),abs(g2%z),abs(g3%z))
!    dk = (min(gmaxx,gmaxy,gmaxz))/nk
    nx = nint(gmaxx/deltak)+1
    ny = nint(gmaxy/deltak)+1
    nz = nint(gmaxz/deltak)+1
    write(ulog,*)'MAKE_KP: gmax_xyz=',gmaxx,gmaxy,gmaxz
    write(ulog,*)'MAKE_KP: del_k,nk=',deltak,nk
    write(ulog,*)'MAKE_KP: nk_xyz  =',nx,ny,nz
    nktot = (nx*ny*nz*8)
    allocate(kpt(3,nktot))

    nk = 0
    do i = -nx,nx
    do j = -ny,ny
    do k = -nz,nz
       q(1) = (i+shft(1))*deltak
       q(2) = (j+shft(2))*deltak
       q(3) = (k+shft(3))*deltak
       call check_inside_fbz(q,g1,g2,g3,inside)
       if (inside) then
          nk = nk+1
          kpt(:,nk) = q(:)
          write(ulog,4)nk,q
          write(fbz,3)q
       endif
    enddo
    enddo
    enddo
    write(ulog,*)nk,' kpoints generated in the FBZ'
!    call allocate_kp(nk)
    allocate( kbz(3,nk),wk(nk) )
    wk(:)=1d0/(nk)
    kbz(:,1:nk) = kpt(:,1:nk)
    deallocate(kpt)

3 format(9(2x,f9.4))
4 format(i6,9(2x,f9.4))

  end subroutine make_kp_FBZ
!-------------------------------------------
  subroutine map_ibz(nk,kp,mapi,ni)
  implicit none
  integer nk,mapi(nk),i,j,l,ni
  real(8) kp(3,nk)
  logical inside

  l=0 ; j=0 ; mapi=0
  do i=1,nk
     call check_inside_irred_fbz(kp(:,i),inside)
     if (inside) then
        j=j+1
        mapi(j) = i
     else
        l=l+1
        mapi(nk-l+1) = i
     endif
  enddo
  ni = j
  if(j+l .ne. nk)then
     write(*,*)' pb in map_ibz j,l,nk=',j,l,nk
  endif
  end subroutine map_ibz
!-------------------------------------------
  subroutine make_kp_MP(nx,ny,nz,sx,sy,sz,nibz,kiz,uio,scal)
! generate a regular mesh(nx,ny,nz) within the box -gi/2,gi/2 + shift(sx,sy,sz)
! and folds some of these vectors within the IFBZ, stored into another kmesh
! kibzc(3,nibz_coarse) 
! sx=sy=sz=0 gives the Monkhorst-Pack scheme
! odd nx,ny,nz contain the gamma point (using sx=sy=sz=0)
    use geometry
    implicit none
!   logical inside
    integer :: i,j,k,nx,ny,nz,ns,nkbz,l,nkfine,nibz,uio
    real(8) sx,sy,sz,q(3),kiz(3,nx*ny*nz),scal
    integer, allocatable:: mapk(:),map2(:)
    real(8), allocatable :: kl(:),junk(:,:),kp_fine(:,:),kz(:,:)

    open(123,file='KPOINT.MP',status='unknown')
    write(ulog,*)'MAKE_KP_MP: nk_xyz  =',nx,ny,nz
    write(ulog,5)'MAKE_KP_MP: shift_xyz  =',sx,sy,sz

    allocate( kp_fine(3,nx*ny*nz) )

    nkfine = 0 ; nkbz = 0
    do i = 0,nx-1
    do j = 0,ny-1
    do k = 0,nz-1
       q = dfloat(2*i+1-nx)/(2*nx)*g1 + dfloat(2*j+1-ny)/(2*ny)*g2 &
&        + dfloat(2*k+1-nz)/(2*nz)*g3 + (sx/nx)*g1+ (sy/ny)*g2+ (sz/nz)*g3
       nkfine = nkfine+1
       write(123,3) q
       kp_fine(:,nkfine) = q(:)   ! this ddefines the MP mesh
    enddo
    enddo
    enddo
    write(ulog,*)nkfine,' fine kpoints generated in the -gi/2,gi/2 FBZ parallelipiped'

! the weights must really be calculated from the weight of each point (to be corrected)
    wk(:) = 1d0/nkfine

! now fold them within FBZ
    allocate(mapk(nkfine),kz(3,nkfine))
    call fold_in_fbz(nkfine,kp_fine,nshl,gg,nkbz,kz,mapk)
    write(ulog,*)nkbz,' kpoints generated in the FBZ'
!    deallocate(gg)

! now sort them
    allocate(junk(3,nkbz))
    junk(1:3,1:nkbz) = kz(1:3,1:nkbz)
    allocate(kl(nkbz),map2(nkbz))
    do i=1,nkbz
       kl(i) = length(junk(:,i))
    enddo
    call sort(nkbz,kl,map2,nkbz)
    deallocate(kz)
    allocate(kz(3,nkbz))
    write(ulog,*)'MAKE_KP_MP: WRITING SORTED KPOINTS IN THE FBZ **************'
    do i=1,nkbz
       kz(:,i) = junk(:,map2(i))  * scal  ! expand a bit beyond the FBZ 
!       write(ulog,4)i,length(kz(:,i)),kz(:,i)
        write(uio+30,3)kz(:,i)
    enddo
    deallocate(kl,map2)
!   junk(1:3,1:nkbz) = kz(1:3,1:nkbz)

! now find the ones in the irreducible FBZ (these would also be sorted in length)
    call select_IBZ(nkbz,kz,nibz,junk)
    do i=1,nibz
       kiz(:,i) = junk(:,i)
       write(uio,3)kiz(:,i)
    enddo

    deallocate(junk,kp_fine,kz)
    close(123)

3 format(9(2x,f11.5))
4 format(i6,9(2x,f11.5))
5 format(a,9(2x,f11.5))

  end subroutine make_kp_MP
!----------------------------------------------
 subroutine fold_in_fbz(nk,kp,nshl,gg,nbz,kz,mapz)
! folds the k-mesh defined by kp_fine into the FBZ and stores the result in kz
 use lattice
 use ios
 implicit none
 integer, intent(in) :: nk,nshl
 integer, intent(out) :: nbz,mapz(nk)
 real(8), intent(in) :: kp(3,nk),gg(3,nshl)
 real(8), intent(out) :: kz(3,nk)
 integer i,ns,j
 real(8) qaux(3),dmin,dist
 logical inside

 write(ulog,*)'FOLD IN FBZ: nk,nshl=',nk,nshl
 write(ulog,*)'FOLD IN FBZ: ik,shell,i,q ====================='
 nbz=0
 do i=1,nk
 foldloop: do ns=1,nshl  ! should be at least 26
     qaux = kp(:,i)+gg(:,ns)
     call check_inside_fbz(qaux,g1,g2,g3,inside)
     if(inside) then
! make sure it is different from previous points: find smallest neighbor distance
        dmin = 1d9
        jloop: do j=1,nbz
           dist = length(qaux-kz(:,j))
           if( dist .lt. dmin) then
              dmin = dist
           endif
        enddo jloop
! if smallest distance is not 0 then accept the point
        if ( dmin .gt. 1.d-6 ) then
           nbz=nbz+1
           mapz(nbz) = i
           kz(:,nbz) = qaux(:)
!          write(ulog,3)i,ns,nbz,qaux
        endif
     endif
 enddo foldloop
 enddo
 write(ulog,*)'FOLD IN FBZ:  DONE! ==============================='

3 format(3(i6),9(2x,f8.4))

 end subroutine fold_in_fbz
!----------------------------------------------
  subroutine get_kindex(q,nkc,kpc,ik)
  use  geometry
  implicit none
  integer ik,nkc,i,j
  real(8) kpc(3,nkc),q(3)

  ik=0
  mainlp: do i=1,nkc
     do j=1,3
        if (.not. (abs(q(j)-kpc(j,i)) .myeq. 0d0) ) exit
        if (j.eq.3) then
           ik=i
           exit mainlp
        endif
     enddo
  enddo mainlp
  if (ik.eq.0) then
     write(*,*)'GET_INDEX: could not find the index of',q
     stop
  endif

  end subroutine get_kindex
!----------------------------------------------
 subroutine get_weights(nk,kp)
!! takes as input kp(3,nk) in cartesian and finds based on crystal symmetries, the
!! weights associated with each kpoint, and then normalizes them
!! nibz is the final number of kpoints stored in kibz in the irreducible BZ
!! the other important output is the mapping of the full kmesh onto the ones
!! folded in the irreducible zone : mapibz(j=1:nk) is the index of the k in the IBZ
!! corresponding to the argument j
!! for i=1,nibz mapinv(i) gives the index k of the kpoints generated from
!! n1,n2,n3 loops
 use lattice
 use constants
 use ios
 use geometry
 use params
 implicit none
 integer nk,nkibz,i,j,l,narms,kvecop(48)
 integer, allocatable :: mcor(:)
 real(8) zro,q(3),kvecstar(3,48),sw,skc(3)
 real(8) kp(3,nk) !,wk(nk)
 real(8) , allocatable :: k2(:,:),lg(:),w2(:)
 logical exists,foundit

  open(fbz,file='KPOINT.FBZ',status='unknown')
  open(uibz,file='KPOINT.IBZ',status='unknown')
 allocate(k2(3,nk),w2(nk),mapibz(nk),mapinv(nk))
 zro=0d0  ! shift
 write(ulog,*)'MAKE_KP_IBZ: generating kpoints in the irreducible FBZ '
 write(*   ,*)'MAKE_KP_IBZ: nk,wk(nk),k2(:,nk)'
 call write_out(ulog,'reduced primitivelattice=',primitivelattice)
 nibz=0 ; mapibz=0; mapinv=0; w2=1
 kploop: do i=1,nk
!    q = kp(:,i)
! kp is in cartesian coordinates, we need to convert it to reduced units:
    q(1)=(kp(:,i)  .dot. r1) /2/pi
    q(2)=(kp(:,i)  .dot. r2) /2/pi
    q(3)=(kp(:,i)  .dot. r3) /2/pi
! below the cartesian components of kp are needed
    call getstar(kp(:,i),primitivelattice,narms,kvecstar,kvecop)
!   call getstar(q      ,primitivelattice,narms,kvecstar,kvecop)
 !  write(ulog,4)'i,kp,q=',i,kp(:,i),q
! see if already exists among the previously assigned ones
    exists=.False.
    if(verbose) write(ulog,4)'list of kvecstar(l),l=1,narms for kp_red=',i,q
    lloop: do l=1,narms
    if(verbose)   write(ulog,4)'stars are:',l,kvecstar(:,l)
! set weight for the first kpoint where nibz=0
        jloop: do j=1,nibz
          if (length(k2(:,j) - kvecstar(:,l)) .lt.1d-4) then
! first bring the star in the FBZ, then compare to the existing points
             exists=.True.
             w2(j)=w2(j)+1d0
             mapibz(i)=j
!            write(ulog,4)' this kpoint turned out to exist ',j,k2(:,j)
             exit lloop
          endif
        enddo jloop
    enddo lloop

!   write(ulog,4)' kpoint, folded one  ',i,kp(:,i),exists,j,k2(:,j)
    if(exists ) then
       cycle kploop
    else
       nibz=nibz+1
!  choose the kpoint star in the first quadrant: at least  works for cubic systems     
    !  call getstar(kp(:,i),primitivelattice,narms,kvecstar,kvecop)
    !  foundit=.False.
    !  strloop: do l=1,narms
    !     q= kvecstar(:,l)
    !     if ( q(1).ge.q(2) .and. q(2).ge.q(3) .and. q(3).ge.-1d-5) then
    !        k2(:,nibz)=q
    !        foundit=.True.
    !        exit strloop
    !     endif
    !  enddo strloop

       if (nibz.eq.1) w2(nibz) = 1
       mapibz(i) = nibz
       mapinv(nibz)=i
    !  if (.not. foundit) then    
          k2(:,nibz)=kp(:,i)
          write(ulog,4)'new vector*:',nibz,w2(nibz),k2(:,nibz),length(k2(:,nibz)),q
          write(*,4)'new vector*:',nibz,w2(nibz),k2(:,nibz),length(k2(:,nibz))
    !  else
    !     write(ulog,4)'new vector :',nibz,w2(nibz),k2(:,nibz),length(k2(:,nibz))
    !     write(*,4)'new vector :',nibz,w2(nibz),k2(:,nibz),length(k2(:,nibz))
    !  endif
    endif
 ! test for FCC 
 !  q(:) = kp(:,i)
 !  if ( q(1).ge.q(2) .and. q(2).ge.q(3) .and. q(3).ge.-1d-5) then
 !     write(bzf,3)i,q
 !  endif

 enddo kploop
 write(ulog,*)'GET_WEIGHTS: generated ',nibz,' points in the irreducible FBZ'

! define kibz and wibz ---------------------------------------------
 if ( allocated(wibz)) deallocate(wibz)
 if ( allocated(kibz)) deallocate(kibz)
 allocate(wibz(nibz),kibz(3,nibz))
 wibz(1:nibz) = w2(1:nibz)
 kibz(:,1:nibz)=k2(:,1:nibz)
! k2 is in reduced units
! do j=1,nibz
!    kibz(:,j)=k2(1,j)*g1+k2(2,j)*g2+k2(3,j)*g3 
! enddo

 deallocate(k2,w2)
! kp(:,1+nibz:nk)=9d9 
! wk(1+nibz:nk)=0d0 
 sw = sum(wibz(1:nibz))
 wibz = wibz/sw

! sort and write out the kpoints and their weights-------------------------
! allocate(lg(nibz),mcor(nibz))
! do i=1,nibz
!    lg(i)=length(kibz(:,i))
! enddo
! call sort(nibz,lg,mcor,nibz)
! do i=1,nibz
!   write(ibzc,3)i,wk(mcor(i)),kibz(:,mcor(i)),lg(mcor(i))
! enddo

 if ( allocated(lg)) deallocate(lg)
 if ( allocated(mcor)) deallocate(mcor)
 allocate(lg(nk),mcor(nk))
 do i=1,nk
    lg(i)=length(kp(:,i))
 enddo
 call sort(nk,lg,mcor,nk)
 write(fbz,*)'#i,l,mapibz(l),kp(l),length(kp(l)),kibz(mapibz(l)),wibz(mapibz(l)) l=mcor(i)'
 do i=1,nk
     j=mcor(i)
!    write(fbz,3)i,kp(:,mcor(i)),length(kp(:,mcor(i))),kibz(:,mapibz(mcor(i))),wibz(mapibz(mcor(i)))
    write(fbz,2)i,j,mapibz(j),kp(:,j),length(kp(:,j)),kibz(:,mapibz(j)),wibz(mapibz(j))
 enddo
 write(uibz,*)'#i,kibz(:,i),wibz(i),kibz_reduced,length(kibz(i))'
open(988,file='kibz.dat')
write(988,*)nibz
close(988)
 do i=1,nibz
    q(1)=(kibz(:,i)  .dot. r1) /2/pi
    q(2)=(kibz(:,i)  .dot. r2) /2/pi
    q(3)=(kibz(:,i)  .dot. r3) /2/pi
    write(uibz,3)i,kibz(:,i),wibz(i),q,length(kibz(:,i))
 enddo

 deallocate(mcor,lg)
 close(fbz)
 close(uibz)

2 format(3i7,2x,3(1x,f9.4),2x,f9.4,2x,3(1x,f9.5),3x,f14.8)
3 format(i7,2x,3(1x,f9.4),2x,f14.8,2x,3(1x,f9.5),3x,f9.5)
4 format(a,i7,2x,f9.5,2x,99(1x,f9.5))
5 format(a,2x,99(1x,f9.5))

 end subroutine get_weights
!===========================================================
!     subroutine getstar(kvec,primlatt,narms,kvecstar,kvecop)
!! find atom in data base
!! arguments:
!!     kvec(i) (input), ith dimensionless component of k vector
!!     narms (output), number of star vectors associated with kvec
!!     kvecstar(3,1:narms), all the stars of kvec
!!     kvecop(i), the symmetry operation number for the star vecor i
!     use atoms_force_constants
!     implicit none 
!     
!     integer, intent(out):: narms,kvecop(48)
!     real(8), intent(in) :: kvec(3),primlatt(3,3)
!     real(8), intent(out):: kvecstar(3,48)
!     integer i,j,k,n,ncmp
!     real(8) v(3),v2(3),v3(3),kvecstarp(3,48)
!
!
!     narms=0
!      print*,'lattpgcount=',lattpgcount
!     iloop: do i=1,lattpgcount
! apply symmetry operation to k to get v=kstar      
!       ! call xvmlt(op_kmatrix(1,1,i),kvec,v,3,3,3)
!       v=matmul(op_kmatrix(:,:,i),kvec)
! find the reduced coordinates of v and store in v2        
!       call xmatmlt(v,primlatt,v2,1,3,3,1,3,1)
!
! now check if v2 - any_previous_v2 is integer (differ by a G vector)
! if so, skip; else store this v as a new star vector
!       do j=1,narms
!       ! subtract previous_v2(=kvecstarp) from v2; result is v3
!         call xvsub(v2,kvecstarp(1,j),v3,3)
!         v3=v2-kvecstarp(:,j)
!         do k=1,3
!           n=nint(v3(k))
!           if(ncmp(v3(k)-n).ne.0)exit
!           if(k.eq.3)cycle iloop  ! goto next sym_op iff v3=integer
!         enddo
!       enddo
!       narms=narms+1
!       kvecstar(1:3,narms)=v(1:3)
!       kvecstarp(1:3,narms)=v2(1:3)
!       kvecop(narms)=i
!     enddo iloop
!     end  
!===========================================================

 end module kpoints
!===========================================================
