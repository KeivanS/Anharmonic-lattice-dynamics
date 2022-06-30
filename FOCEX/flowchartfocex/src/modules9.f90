!==============================================
! Modules here
!==============================================

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
! integer nv3,readv3,writev3, iter, split, calk ! sy added iter,split, calk
! real(8) v3_threshold,const33

 contains

  subroutine allocate_v33_sy(n1,n2,l1,l2,l3)
    allocate(v33_md(n1,n2,l1,l2,l3) )
  end subroutine allocate_v33_sy
  subroutine deallocate_v33_sy
    deallocate(v33_md)
  end subroutine deallocate_v33_sy
  subroutine allocate_v33_sq(n1,n2,l1,l2,l3)
    allocate(v33sq1(n1,n2,l1,l2,l3))
  end subroutine allocate_V33_sq
  subroutine deallocate_v33_sq
    deallocate(v33sq1)
  end subroutine
  subroutine allocate_phi3_sy(nk,nl,ni)
  integer nk,nl,ni
!   allocate(delta3(nk,nl,nw),gamma3(nk,nl,nw),v3(nk,nk,nl,nl,nl))
!   allocate(selfN(ni,nl,nw),selfU(ni,nl,nw))
    allocate(v33_md(ni,nk,nl,nl,nl))
  end subroutine allocate_phi3_sy

 end module phi3_sy


!===========================================================
module exactBTE2

integer ncol,nlin
real(8), allocatable:: F_RTA(:,:,:), F1(:,:,:), F2(:,:,:), Qvalue(:,:), iselfenergy(:,:), tau(:,:,:), error(:,:,:), diff(:,:,:)
real(8), allocatable:: mval(:,:),qval(:),coll(:,:),rhs(:,:),rhsi(:,:),ne_dist(:,:),ne_dist_rta(:,:)
real(8), allocatable:: Qvalue_N(:,:), Qvalue_U(:,:) !, tauinv_N(:,:), tauinv_U(:,:), tauinv_tot(:,:), tauinv_eff(:,:,:)
real(8), allocatable:: diff_kap(:,:,:,:),kappa(:,:,:), kappa_k(:,:,:,:), kappa_RTA(:,:,:), kappa_k_RTA(:,:,:,:)
! F1,F2: (kp1,la1,xyz)    , Qvalue,Avalue,iselfenergy: (kp,la)
real(8), allocatable:: P1(:,:,:,:,:), P2(:,:,:,:,:)
! P1,P2: (kp1,kp2,la1,la2,la3)   don't need kp3 since it is determined by momentum conservation.
real(8), allocatable:: frequency(:,:) ,dist(:,:)   ! BE distribution function
!real(8), allocatable:: v33sq(:,:,:,:,:)


contains

!subroutine allocate_v33sq(nk,ndn)
!allocate(v33sq(nk,nk,ndn,ndn,ndn))
!end subroutine allocate_v33sq

subroutine allocate_iter_kap (nk,ndn)
 implicit none
 integer nk,ndn
!allocate (F_RTA(nk,ndn,3),F1(nk,ndn,3),F2(nk,ndn,3),F1_old(nk,ndn,3),Qvalue(nk,ndn), &
!&         iselfenergy(nk,ndn),tau(nk,ndn,3),error(nk,ndn,3),diff(nk,ndn,3))
 allocate (diff_kap(nk,ndn,3,3),kappa(ndn,3,3) ) !,kappa_k(nk,ndn,3,3),kappa_RTA(ndn,3,3), &
!&         kappa_k_RTA(nk,ndn,3,3))    ! 3 for xyz
!allocate (dist(nk,ndn),frequency(nk,ndn))
!allocate(tauinv_N(nk,ndn),tauinv_U(nk,ndn),tauinv_tot(nk,ndn),tauinv_eff(nk,ndn,3))
end subroutine allocate_iter_kap

subroutine allocate_iter0 (ni,nk,ndn)
allocate (Qvalue_N(ni,ndn),Qvalue_U(ni,ndn),qval(ni*ndn),mval(nk*ndn,ni*ndn))
end subroutine allocate_iter0

subroutine allocate_iter1 (ni,nk,ndn)
allocate (rhs(ndn*nk,3),rhsi(ndn*ni,3),coll(ndn*ni,ndn*ni),ne_dist(ni*ndn,3),ne_dist_rta(ni*ndn,3))
end subroutine allocate_iter1

subroutine allocate_iter2 (n1,n2,l1,l2,l3)
allocate (P1(n1,n2,l1,l2,l3),P2(n1,n2,l1,l2,l3))
end subroutine allocate_iter2

subroutine deallocate_iter_kap
!deallocate(dist,frequency)
!deallocate(F_RTA,F1,F2,Qvalue,iselfenergy,tau,error,diff,F1_old)
deallocate(diff_kap,kappa) !,kappa_k,kappa_RTA,kappa_k_RTA)
end subroutine deallocate_iter_kap

subroutine deallocate_iter0
deallocate(Qvalue_N,Qvalue_U,qval,mval)
end subroutine deallocate_iter0

subroutine deallocate_iter1
deallocate(rhs,rhsi,coll) !tauinv_N,tauinv_U,tauinv_tot,tauinv_eff)
end subroutine deallocate_iter1

subroutine deallocate_iter2
deallocate(P1,P2)
end subroutine deallocate_iter2

!subroutine deallocate_v33sq
!deallocate(v33sq)
!end subroutine deallocate_V33sq

end module exactBTE2

!==============================================================
 module constants
 implicit none
 real(8) :: pi=3.14159265358979323846d0
 complex(8) :: ci=dcmplx(0d0,1d0)
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
 real(8) tolerance,margin,scalelengths
 integer nconfigs,classical,ntemp,fdfiles,cal_cross,threemtrx,lamin,lamax,ncpu,n_dig_acc,isvd
 integer nshells(4,10)   ! up to which shell to include for each rank of FC
 integer include_fc(4)  ! whether to include FCs of that rank
 real(8) rcut(4),tau0,wshift(3)
 real(8) tmin,tmax,qcros(3),svdc
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
&         ufc=80,ufit1=31,ufit2=32,ufit3=33,ufit4=34,ucor=93,uborn=12
real(8) temperature

  interface write_out
     module procedure write_outiv, write_outrv, write_outim, write_outrm &
&    , write_outi, write_outr , write_outv
  end interface

 contains

!-----------------------------
  subroutine write_outrm(unit,string,n,m,var)
  implicit none
  character*(*) string
  integer n,m,unit,l,i,j
  real(8) var(n,m)

  l=len_trim(string)
  write(unit,4)string(1:l)//' is=',((var(i,j),i=1,n),j=1,m)
4 format(a,99(1x,g13.6))
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
  integer unit,l
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
! Module for saved data
      module force_constants_module
! maximum number of shells of nearest neighbors
      integer maxneighbors
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
      integer, allocatable :: iatomop(:,:)
! atomopfract(k,j,i), kth cartesian coordinate of fractional to accompany iatomop
      double precision, allocatable :: atomopfract(:,:,:)
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
! iatomcell0(i), identity of atom in unit cell at origin equivalent to ith atom
!  integer iatomcell0(maxatoms)
! iatomneighbor(j,i), nearest neighbor shell of jth atom in primitive unit
! cell at the origin that contains the ith atom
      integer, allocatable:: iatomneighbor(:,:)
!      allocatable iatomneighbor,iatomop,atomopfract
      end module force_constants_module
!============================================================
 module lattice
 use geometry
 use constants
 implicit none
 type(vector) r1,r2,r3,g1,g2,g3  ,rr1,rr2,rr3   ! translation vectors of the supercell
 type(vector) r01,r02,r03,g01,g02,g03  ! tr vect of prim cell and its recip spce
 real(8) volume_r,volume_g,lattice_parameter,latticeparameters(6),primitivelattice(3,3)
 real(8) box(3),boxg(3)
 real(8) r0g(3,3)
 integer n1min,n2min,n3min,n1max,n2max,n3max,NC(3),NF(3)
 integer nr1(3),nr2(3),nr3(3)

! contains
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

 end module lattice
!===========================================================
 module atoms_force_constants
! the type atom_id0 concerns the properties of the primitive unit cell
! it s atoms, masses, neighbors, shells and force constants, whereas the
! type atomic refers to atoms in the supercell, their displacements and forces.
 use geometry
 use force_constants_module
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
    integer at_type
    real(8) mass
 end type
!-------------------------
 type atom_id0   ! all characteristics of an atom in the primitive central cell
    integer at_type
    character name*2
    real(8) mass,charge
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
 character(2), allocatable:: atname(:)
 type (atom_id0), allocatable :: atom0(:)
 type (atomic), allocatable :: atom_sc(:),atom_shl(:)

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
 use force_constants_module
 use params
 use lattice
 use ios
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
!      if ( iatomneighbor(i0,j) .eq. nshells(2) .and. dist(j) .gt. rmax ) then
       if ( iatomneighbor(i0,j) .eq. nshells(2,i0) .and. dist(j) .gt. rmax ) then
            rmax = dist(j)
       endif
    enddo
    call sort(natoms,dist,msort,maxatoms)

    shel_count = -1
    d_old = 0
    write(ulog,*)' ========================================================='
    write(ulog,*)' neighborlist j, of atom i0 in atompos array =',i0

    jloop: do j=1,min(500,natoms)  ! assign the first 500 neighbors
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
       write(ulog,2)' shell# , neighb#, neighbdist =',shel_count,counter,(atompos(l,msort(j)),l=1,3),d_min
!       write(ulog,6)' neighbor distance=',d_min
!      write(ulog,4)' neighborshel,xyz =',iatomneighbor(i0,msort(j))
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

2 format(a,2(2x,i4),', (',3(1x,f8.4),') ,',3x,f10.5)
!3 format(a,' (',i4,2(1x,i4),1x,')')
!4 format(a,1x,i4,' (',f8.4,2(1x,f8.4),1x,')')
5 format(a,3(2x,i4),' (',i4,2(1x,i4),1x,')')
!6 format(a,2x,g16.7)

 end subroutine set_neighbor_list

 end module atoms_force_constants
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
    character(1), allocatable:: err(:)
    integer ngr,ntotind,ntot  ! number of groups, total number of independent terms and full terms
 end type fulldmatrix

 integer maxrank,maxgroups,itrans,irot,ihuang,enforce_inv
 integer nterms(4),maxterms(4),maxtermsindep(4),ngroups(4)
! integer ndindp(4),ndfull(4)
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
 character(1), allocatable:: err_1(:),err_2(:),err_3(:),err_4(:)
 type(fulldmatrix) map(4)
 real(8), allocatable:: ampterm_1(:),fcs_1(:)
 real(8), allocatable:: ampterm_2(:),fcs_2(:),grun_fc(:)
 real(8), allocatable:: ampterm_3(:),fcs_3(:)
 real(8), allocatable:: ampterm_4(:),fcs_4(:)
 real(8), allocatable:: amat(:,:),bmat(:),sigma(:),fcs(:),ahom(:,:),overl(:,:)
 real(8), allocatable:: a11ia12(:,:),fc1(:)
 real(8), allocatable:: arot(:,:),brot(:),atransl(:,:),ahuang(:,:),aforce(:,:),bforce(:)
 real(8), allocatable:: fc_ind(:) ! This is for getting force constant by reading FC2.dat file
 integer inv_constraints,force_constraints,dim_al,dim_ac,n_indep,newdim_al,dim_hom   &
&        ,transl_constraints, rot_constraints, huang_constraints,ngr
real(8), allocatable:: energy(:)
contains

 subroutine set_maxterms
   maxterms(1)=15
   maxterms(2)=200
   maxterms(3)=200
   maxterms(4)=600
   maxtermsindep(1)=5
   maxtermsindep(2)=30
   maxtermsindep(3)=30
   maxtermsindep(4)=60
   maxgroups=90
 end subroutine set_maxterms

  end module svd_stuff
!============================================================
 function bring_to_prim_cell_cav(a) result(w)
! takes cartesian coordinates and returns cart coordinates within the primcell
 use lattice
 use geometry
 implicit none
 type(vector) v,w, bring_to_prim_cell_cv
 real(8) a(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w =  bring_to_prim_cell_cv(v)

 end function bring_to_prim_cell_cav
!-----------------------------------------
 function bring_to_super_cell_cav(a) result(w)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
 use geometry
 implicit none
 type(vector) v,w, bring_to_super_cell_cv
 real(8) a(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w = bring_to_super_cell_cv(v)

 end function bring_to_super_cell_cav
!-----------------------------------------
 function bring_to_prim_cell_caa(a) result(b)
! takes cartesian coordinates and returns cart coordinates within the primcell
 use lattice
 use geometry
 implicit none
 type(vector) v,w, bring_to_prim_cell_cv
 real(8) a(3),b(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w =  bring_to_prim_cell_cv(v)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end function bring_to_prim_cell_caa
!-----------------------------------------
 function bring_to_super_cell_caa(a) result(b)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
 use geometry
 implicit none
 type(vector) v,w, bring_to_super_cell_cv
 real(8) a(3),b(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w = bring_to_super_cell_cv(v)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end function bring_to_super_cell_caa
!-----------------------------------------
 function bring_to_prim_cell_cv(v) result(w)
! takes cartesian coordinates and returns cart coordinates within the primcell
 use lattice
 use geometry
 implicit none
 type(vector) v,w
 real(8) a1,a2,a3

! get direct coordinates first
 a1 = v .dot. g01
 a2 = v .dot. g02
 a3 = v .dot. g03
! bring into cell: between 0 and cell
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
! convert to cartesian coordinates
 w = a1*r01+a2*r02+a3*r03

 end function bring_to_prim_cell_cv
!-----------------------------------------
 function bring_to_super_cell_cv(v) result(w)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
 use geometry
 implicit none
 type(vector) v,w
 real(8) a1,a2,a3

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

 end function bring_to_super_cell_cv
!-----------------------------------
 subroutine cart_to_direct_v(v,w)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
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
 use lattice
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
 use lattice
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
 use lattice
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
 use lattice
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
 use lattice
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
 use lattice
 real(8) q(3),k(3)

 k(:) = q(1)*g1 + q(2)*g2 + q(3)*g3

 end subroutine dir2cart_g
!===========================================================
  module born
  real(8) epsil(3,3)
  real(8), allocatable:: zeu(:,:,:)
  real(8) rho
  integer born_flag
  end module born
!===========================================================

