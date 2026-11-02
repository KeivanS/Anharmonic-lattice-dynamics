 module geometry
 use constants, only : r15

  type vector
    real(r15) :: x,y,z
  end type vector
  type point
    real(r15) :: x,y,z
  end type point
  type svector
    real(r15) :: s,v(3)
  end type svector
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

   interface operator(.myeq.)
     module procedure myequal, myequal_i,myequalarray_i, &
&         myequalarray_r, myequalvector ; end interface

   interface myeqz
     module procedure myequal0, myequal0array, myequal0vector ; end interface

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

  interface trace
     module procedure trace_r,trace_c,trace_i
  end interface trace

   contains

!-----------------------------------------
 function trace_r(mat) result(tr)
 use constants, only : r15
 implicit none
 integer i,n
 real(r15), dimension(:,:), intent(in) :: mat
 real(r15) tr 

 n=size(mat,1) !; write(*,*)'TRACE: n=',n
! write(*,*)'TRACE: n=',n,' mat=',mat
 tr=0
 do i=1,n
    tr=tr+mat(i,i)
 enddo

 end function trace_r
!-----------------------------------------
 function trace_c(mat) result(tr)
 use constants, only : r15
 implicit none
 integer i,n
 complex(r15), dimension(:,:), intent(in) :: mat
 complex(r15) tr 

 n=size(mat,1) !; write(*,*)'TRACE: n=',n
! write(*,*)'TRACE: n=',n,' mat=',mat
 tr=0
 do i=1,n
    tr=tr+mat(i,i)
 enddo

 end function trace_c
!-----------------------------------------
 function trace_i(mat) result(tr)
 use constants, only : r15
 implicit none
 integer i,n
 integer, dimension(:,:), intent(in) :: mat
 integer tr 

 n=size(mat,1); 
! write(*,*)'TRACE: n=',n,' mat=',mat
 tr=0
 do i=1,n
    tr=tr+mat(i,i)
 enddo

 end function trace_i
!-----------------------------------------
 function is_integer_1(i) result(is)
!! checks if the variable is integer
 use params, only : tolerance
 implicit none
 real(r15), intent(in):: i
 logical is

 if (abs(i-nint(i)).gt.tolerance ) then
! if (abs(i-nint(i)).gt.1d-6 ) then
    is=.false.
    return
 endif
 is=.true.

 end function is_integer_1
!-----------------------------------------
 function is_integer_a(i) result(is)
!! checks if the variable is integer
 use params, only : tolerance
 implicit none
 real(r15), intent(in):: i(:)
 logical is
 integer j,n

 n=size(i)
 do j=1,n
   if(abs(i(j)-nint(i(j))).gt.tolerance ) then
!  if(abs(i(j)-nint(i(j))).gt.1d-6 ) then
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
 if(myeqz(om)) then
    write(*,*)'CALVOL: Volume is zero! check your basis vectors!'
    write(*,'(a,9(1x,g12.5))')'r1,r2,r3=',r1,r2,r3
!   stop
 endif

end subroutine calvol_a
!-----------------------------------
subroutine calvol_v(r1,r2,r3,om)
implicit none
real(r15) om
type(vector) cross12,r1,r2,r3

cross12 = r1 .cross. r2
om = abs(r3 .dot. cross12)
 if(myeqz(om)) then
    write(*,*)'CALVOL: Volume is zero! check your basis vectors!'
    write(*,'(a,9(1x,g12.5))')'r1,r2,r3=',r1,r2,r3
!   stop
 endif

end subroutine calvol_v
!------------------------------------------------------------------------------
      function rdet(mat)
      implicit none
!
!!	FIND THE DETERMINANT OF A 3X3 real MATRIX MAT
!
      real(r15) rdet,mat(3,3)
      rdet=mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))  &
     &   + mat(1,2)*(mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3))  &
     &   + mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))

      end function rdet
!------------------------------------------------------------------------------
      function idet(mat)
      implicit none
!
!! FIND THE DETERMINANT OF A 3X3 integer MATRIX MAT
!
      integer idet,mat(3,3)
      idet=mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))  &
     &   + mat(1,2)*(mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3))  &
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
   function myequal0(v) result(eq)
     use params
     real(r15), intent(in)::  v
     logical eq !, intent(out) :: eq

        if (abs(v) .lt. tolerance) then
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
   function myequal0vector(v) result(eq)
     type(vector), intent(in)::  v
     logical eq !, intent(out) :: eq
  if ( (v%x .myeq. 0d0) .and. (v%y .myeq. 0d0) .and. (v%z .myeq. 0d0) ) then
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
   function myequal0array(v) result(eq)
     real(r15), dimension(:), intent(in) ::  v
     logical eq !, intent(out) :: eq
     integer i,n
     i = size(v) 
        eq = .true.
        loop: do n=1,i
           if ( .not. myequal0(v(n)) ) then
              eq = .false.
              exit loop
           endif
        enddo loop
     
   end function myequal0array
!-----------------------------------
   function myequalarray_r(v,w) result(eq)
     use params
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
           if ( abs(v(n)-w(n)) .gt. tolerance ) then
              eq = .false.
              exit loop
           endif
        enddo loop
     endif
   end function myequalarray_r
!-----------------------------------
   function myequalarray_i(v,w) result(eq)
     use params
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
           if ( abs(v(n)-w(n)) .gt. tolerance ) then
  !        if ( .not. myequal_i(v(n),w(n)) ) then
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
