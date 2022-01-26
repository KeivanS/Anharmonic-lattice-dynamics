!!THIS MODULE IS SEPERATED FROM modules9.f90
!!legacy code, no comment
!==============================================
module params
 real(8) tolerance,margin,alfaborn
 integer nconfigs,classical,ntemp,fdfiles,cal_cross,threemtrx,lamin,lamax,ncpu,n_dig_acc,isvd
 integer nshells(4,20)   ! up to which shell to include for each rank of FC
 integer include_fc(4)  ! whether to include FCs of that rank
 real(8) rcut(4),tau0,wshift(3)
 real(8) tmin,tmax,qcros(3),svdcut
 logical verbose
!********************************************
INTEGER,PARAMETER :: d=3 !dimension
!********************************************
end module params
!===============================================
!THIS MODULE IS SEPERATED FROM modules9.f90
!with changed type declaration, removed interface operator not used
!added/adjusted original interace operator
!============================================================
 module geometry
  !use Structure_info !now the dimension is in module [params]
  use params
  IMPLICIT NONE
  type vector
    real(8),dimension(d) :: component
  end type vector
  !type point
    !real(8),dimension(d) :: coordinate
  !end type point
!------------------------------------------------------------
    interface operator(+)
     module procedure addition_vva,addition_ava,addition_vaa!,&
     !&addition,addition_avv,,addition_vav
   end interface
!------------------------------------------------------------
    interface operator(-)
     module procedure subtraction,subtraction_av,subtraction_va;
    end interface
!------------------------------------------------------------
    interface operator(*)
     module procedure multiply_by_real,multiply_by_integer
    end interface
!------------------------------------------------------------
    interface operator(/)
     module procedure divide_by_scalar
    end interface
!------------------------------------------------------------
    interface operator(.cross.) !only works for d=3??
     module procedure crossproduct_v,crossproduct_a
    end interface
!------------------------------------------------------------
    interface operator(.myeq.)
     module procedure myequal, myequalarray, myequalvector
    end interface
!------------------------------------------------------------
    interface operator(.myeqz.)
     module procedure myequal0, myequal0array, myequal0vector
    end interface
!------------------------------------------------------------
!-------------------------------------------------------------
    interface operator(.dot.)
        module procedure dotproduct_v,dotproduct_a,  &
&                         dotproduct_av,dotproduct_va,dotproduct_sv,&
&                         dotproduct_vt2,dotproduct_vt3,dotproduct_vt4,&
&                         dotproduct_t2v,dotproduct_t3v,dotproduct_t4v,&

&                         dotproductC_a,dotproductC_sv,&
&                         dotproductC_vt2,&
&                         dotproductC_t2v,dotproductC_t3v,dotproductC_t4v,&

&                         dotproductM_a,dotproductM_t2v,dotproductM_t3v,dotproductM_t4v,&
&                         dotproductIM_a,dotproductIM_t2v,dotproductIM_t3v,dotproductIM_t4v,&

&                         doubledot,tripledot,tripledotR,quadrupledot,quadrupledotR
    end interface
!------------------------------------------------------------
    interface operator(.idot.)
        module procedure doubledot2,tripledot2,quadrupledot2
    end interface
!------------------------------------------------------------
    interface operator(.odot.)
        module procedure doubledot3
    end interface
!-------------------------------------------------------------
    interface length
     module procedure lengthv,lengtha
    end interface

contains
!============================================================
!AMBIGUOUS AND NOT NECESSARY!
    !function addition(v,w) result(add)
    ! type(vector) add !, intent(out) ::
    ! type(vector), intent(in) :: v,w
    ! add%component(1) = v%component(1) + w%component(1)
    ! add%component(2) = v%component(2) + w%component(2)
    ! add%component(3) = v%component(3) + w%component(3)
    !end function addition
!-----------------------------------
    function addition_vaa(v,w) result(add)
     real(8) add(3) !, intent(out) ::
     type(vector), intent(in) :: v
     real(8), intent(in) :: w(3)
     add(1) = v%component(1) + w(1)
     add(2) = v%component(2) + w(2)
     add(3) = v%component(3) + w(3)
    end function addition_vaa
!-----------------------------------
   function addition_ava(w,v) result(add)
     real(8) add(3) !, intent(out) ::
     type(vector), intent(in) :: v
     real(8), intent(in) :: w(3)
     add(1) = v%component(1) + w(1)
     add(2) = v%component(2) + w(2)
     add(3) = v%component(3) + w(3)
   end function addition_ava
!-----------------------------------
   function addition_vva(w,v) result(add)
     real(8) :: add(3)
     type(vector), intent(in) :: w,v
     add(1) = w%component(1)+v%component(1)
     add(2) = w%component(2)+v%component(2)
     add(3) = w%component(3)+v%component(3)
   end function addition_vva
!-----------------------------------!AMBIGUOUS AND NOT NECESSARY!
   !function addition_avv(w,v) result(add)
   !  type(vector) add !, intent(out) ::
   !  type(vector), intent(in) :: v
   !  real(8), intent(in) :: w(3)
   !  add%component(1) = v%component(1) + w(1)
   !  add%component(2) = v%component(2) + w(2)
   !  add%component(3) = v%component(3) + w(3)
   !end function addition_avv
!----------------------------------- !AMBIGUOUS AND NOT NECESSARY!
   !function addition_vav(v,w) result(add)
   !  type(vector) add !, intent(out) ::
   !  type(vector), intent(in) :: v
   !  real(8), intent(in) :: w(3)
   !  add%component(1) = v%component(1) + w(1)
   !  add%component(2) = v%component(2) + w(2)
   !  add%component(3) = v%component(3) + w(3)
   !end function addition_vav
!-----------------------------------
   function subtraction(v,w) result(dif)
     real(8) dif(3) !, intent(out) ::
     type(vector), intent(in) :: v,w
     dif(1) = v%component(1) - w%component(1)
     dif(2) = v%component(2) - w%component(2)
     dif(3) = v%component(3) - w%component(3)
   end function subtraction
!-----------------------------------
   function subtraction_va(v,w) result(dif)
     real(8) dif(3) !, intent(out) ::
     type(vector), intent(in) :: v
     real(8), intent(in) :: w(3)
     dif(1) = v%component(1) - w(1)
     dif(2) = v%component(2) - w(2)
     dif(3) = v%component(3) - w(3)
   end function subtraction_va
!-----------------------------------
   function subtraction_av(w,v) result(dif)
     real(8) dif(3) !, intent(out) ::
     type(vector), intent(in) :: v
     real(8), intent(in) :: w(3)
     dif(1) = v%component(1) - w(1)
     dif(2) = v%component(2) - w(2)
     dif(3) = v%component(3) - w(3)
   end function subtraction_av
!-----------------------------------
   function multiply_by_integer(s,v) result(sv)
   integer, intent(in) :: s
   type(vector), intent(in) :: v
   real(8) :: sv(3)
   !type(vector) sv !, intent(out)::
     sv(1)=s*v%component(1)
     sv(2)=s*v%component(2)
     sv(3)=s*v%component(3)
   end function multiply_by_integer
!-----------------------------------
   function multiply_by_real(s,v) result(sv)
   real(8), intent(in) :: s
   type(vector), intent(in) :: v
   type(vector) sv !, intent(out)::
     sv%component(1)=s*v%component(1)
     sv%component(2)=s*v%component(2)
     sv%component(3)=s*v%component(3)
   end function multiply_by_real
!-----------------------------------
   function divide_by_scalar(v,s) result(sv)
   real(8), intent(in) :: s
   type(vector), intent(in) :: v
   type(vector) sv !, intent(out)::
     sv%component(1)=1/s*v%component(1)
     sv%component(2)=1/s*v%component(2)
     sv%component(3)=1/s*v%component(3)
   end function divide_by_scalar
!-----------------------------------
   function crossproduct_v(v,w) result(cross)
     type(vector) cross !, intent(out) ::
     type(vector), intent(in) :: v,w

     cross%component(1) = v%component(2)*w%component(3)-v%component(3)*w%component(2)
     cross%component(2) = v%component(3)*w%component(1)-v%component(1)*w%component(3)
     cross%component(3) = v%component(1)*w%component(2)-v%component(2)*w%component(1)
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
     real(8) dot
     type(vector), intent(in) :: v,w
     dot = sum(v%component * w%component)
   end function dotproduct_v
!-----------------------------------
   function dotproduct_a(v,w) result(dot)
     real(8) dot
     real(8), intent(in) :: v(:),w(:)
     dot = sum(v*w)
   end function dotproduct_a
!-----------------------------------
   function dotproduct_av(v,w) result(dot)
     real(8) dot
     real(8), intent(in) :: v(d)
     type(vector), intent(in) :: w
     dot = sum(v * w%component)
   end function dotproduct_av
!-----------------------------------
   function dotproduct_va(w,v) result(dot)
     real(8) dot
     real(8), intent(in) :: v(d)
     type(vector), intent(in) :: w
     dot = sum(v * w%component)
   end function dotproduct_va
!-----------------------------------

!****************************************************************************************************************
!expand .dot. operation for tensors
function dotproduct_vt2(v,w) result(dot)
     real(8), dimension(d) :: dot
     real(8), intent(in),dimension(d) :: v
     real(8), intent(in),dimension(d,d) :: w
     integer i
     do i=1,d
        dot(i) = sum(v*w(:,i))
     end do
end function dotproduct_vt2
!---------------------------------------
function dotproduct_sv(v,w) result(dot)
    real(8),dimension(d) :: dot
    real(8),intent(in):: v
    real(8),intent(in),dimension(d) :: w
    dot=v*w
end function
!---------------------------------------
function dotproduct_t2v(v,w) result(dot)
     real(8), dimension(d) :: dot
     real(8), intent(in),dimension(d) :: w
     real(8), intent(in),dimension(d,d) :: v
     integer i
     do i=1,d
        dot(i) = sum(v(i,:)*w)
     end do
end function dotproduct_t2v
!----------------------------------------
function dotproduct_t3v(v,w) result(dot)
    real(8),dimension(d,d) :: dot
    real(8),intent(in),dimension(d) :: w
    real(8),intent(in),dimension(d,d,d) :: v
    integer i
    do i=1,d
        dot(i,:)=dotproduct_t2v(v(i,:,:),w)
    end do
end function dotproduct_t3v
!-------------------------------------------
function dotproduct_vt3(v,w) result(dot)
    real(8),dimension(d,d) :: dot
    real(8),intent(in),dimension(d,d,d) :: w
    real(8),intent(in),dimension(d) :: v
    integer i
    do i=1,d
        dot(:,i)=dotproduct_vt2(v,w(:,:,i))
    end do
end function dotproduct_vt3
!-------------------------------------------
function dotproduct_t4v(v,w) result(dot)
    real(8),dimension(d,d,d) :: dot
    real(8),intent(in),dimension(d) :: w
    real(8),intent(in),dimension(d,d,d,d) :: v
    integer i
    do i=1,d
        dot(i,:,:)=dotproduct_t3v(v(i,:,:,:),w)
    end do
end function dotproduct_t4v
!-------------------------------------------
function dotproduct_vt4(v,w) result(dot)
    real(8),dimension(d,d,d) :: dot
    real(8),intent(in),dimension(d,d,d,d) :: w
    real(8),intent(in),dimension(d) :: v
    integer i
    do i=1,d
        dot(:,:,i)=dotproduct_vt3(v,w(:,:,:,i))
    end do
end function dotproduct_vt4
!****************************************************************************************************************
!expand .dot. operation for complex number
function dotproductC_vt2(v,w) result(dot)
     complex(8), dimension(d) :: dot
     complex(8), intent(in),dimension(d) :: v
     complex(8), intent(in),dimension(d,d) :: w
     integer i
     do i=1,d
        dot(i) = sum(v*w(:,i))
     end do
end function dotproductC_vt2
!-----------------------------------------
function dotproductC_sv(v,w) result(dot)
    complex(8),dimension(d) :: dot
    complex(8),intent(in):: v
    complex(8),intent(in),dimension(d) :: w
    dot=v*w
end function dotproductC_sv
!----------------------------------------------
 function dotproductC_a(v,w) result(dot)
     complex(8) dot
     complex(8), intent(in) :: v(:),w(:)
     dot = sum(v*w)
   end function dotproductC_a
!----------------------------------------------
 function dotproductM_a(v,w) result(dot)
     complex(8) dot
     complex(8), intent(in) :: v(:)
     real(8),intent(in)::w(:)
     dot = sum(v*w)
   end function dotproductM_a
!----------------------------------------------
 function dotproductIM_a(v,w) result(dot)
     complex(8) dot
     real(8), intent(in) :: v(:)
     complex(8),intent(in)::w(:)
     dot = sum(v*w)
   end function dotproductIM_a
!----------------------------------------------
function dotproductC_t2v(v,w) result(dot)
     complex(8), dimension(d) :: dot
     complex(8), intent(in),dimension(d) :: w
     complex(8), intent(in),dimension(d,d) :: v
     integer i
     do i=1,d
        dot(i) = sum(v(i,:)*w)
     end do
end function dotproductC_t2v
!---------------------------------------------
function dotproductM_t2v(v,w) result(dot)
     complex(8), dimension(d) :: dot
     real(8), intent(in),dimension(d) :: w
     complex(8), intent(in),dimension(d,d) :: v
     integer i
     do i=1,d
        dot(i) = sum(v(i,:)*w)
     end do
end function dotproductM_t2v
!----------------------------------------------
function dotproductIM_t2v(v,w) result(dot)
     complex(8), dimension(d) :: dot
     complex(8), intent(in),dimension(d) :: w
     real(8), intent(in),dimension(d,d) :: v
     integer i
     do i=1,d
        dot(i) = sum(v(i,:)*w)
     end do
end function dotproductIM_t2v
!---------------------------------------------
function dotproductC_t3v(v,w) result(dot)
    complex(8),dimension(d,d) :: dot
    complex(8),intent(in),dimension(d) :: w
    complex(8),intent(in),dimension(d,d,d) :: v
    integer i
    do i=1,d
        dot(i,:)=dotproductC_t2v(v(i,:,:),w)
    end do
end function dotproductC_t3v
!---------------------------------------------
function dotproductM_t3v(v,w) result(dot)
    complex(8),dimension(d,d) :: dot
    real(8),intent(in),dimension(d) :: w
    complex(8),intent(in),dimension(d,d,d) :: v
    integer i
    do i=1,d
        dot(i,:)=dotproductM_t2v(v(i,:,:),w)
    end do
end function dotproductM_t3v
!---------------------------------------------
function dotproductIM_t3v(v,w) result(dot)
    complex(8),dimension(d,d) :: dot
    complex(8),intent(in),dimension(d) :: w
    real(8),intent(in),dimension(d,d,d) :: v
    integer i
    do i=1,d
        dot(i,:)=dotproductIM_t2v(v(i,:,:),w)
    end do
end function dotproductIM_t3v
!------------------------------------------
function dotproductC_t4v(v,w) result(dot)
    complex(8),dimension(d,d,d) :: dot
    complex(8),intent(in),dimension(d) :: w
    complex(8),intent(in),dimension(d,d,d,d) :: v
    integer i
    do i=1,d
        dot(i,:,:)=dotproductC_t3v(v(i,:,:,:),w)
    end do
end function dotproductC_t4v
!-----------------------------------------
function dotproductM_t4v(v,w) result(dot)
    complex(8),dimension(d,d,d) :: dot
    real(8),intent(in),dimension(d) :: w
    complex(8),intent(in),dimension(d,d,d,d) :: v
    integer i
    do i=1,d
        dot(i,:,:)=dotproductM_t3v(v(i,:,:,:),w)
    end do
end function dotproductM_t4v
!-----------------------------------------
function dotproductIM_t4v(v,w) result(dot)
    complex(8),dimension(d,d,d) :: dot
    complex(8),intent(in),dimension(d) :: w
    real(8),intent(in),dimension(d,d,d,d) :: v
    integer i
    do i=1,d
        dot(i,:,:)=dotproductIM_t3v(v(i,:,:,:),w)
    end do
end function dotproductIM_t4v
!-----------------------------------------
function doubledot(v,w) result(dot)
    real(8) :: dot
    real(8),dimension(d) :: temp
    real(8),intent(in),dimension(d,d) ::v,w
    integer i
    do i=1,d
        temp(i)=sum(v(i,:)*w(i,:))
    end do
    dot=sum(temp)
end function doubledot
!-----------------------------------------
function tripledot(v,w) result(dot)
    real(8),dimension(d) :: dot
    real(8),dimension(d,d) :: temp
    real(8),intent(in),dimension(d,d) ::w
    real(8),intent(in),dimension(d,d,d) ::v
    integer i,j
    do i=1,d
        do j=1,d
            temp(i,j)=sum(v(i,j,:)*w(j,:))
        end do
        dot(i)=sum(temp(i,:))
    end do
end function tripledot
!-----------------------------------------
function tripledotR(v,w) result(dot)
    real(8),dimension(d) :: dot
    real(8),dimension(d,d) :: temp
    real(8),intent(in),dimension(d,d) ::v
    real(8),intent(in),dimension(d,d,d) ::w
    integer i,j
    do i=1,d
        do j=1,d
            temp(j,i)=sum(v(:,j)*w(:,j,i))
        end do
        dot(i)=sum(temp(:,i))
    end do
end function tripledotR
!-----------------------------------------
function quadrupledot(v,w) result(dot)
    real(8),dimension(d,d) :: dot
    real(8),intent(in),dimension(d,d,d,d) :: v
    real(8),intent(in),dimension(d,d) :: w
    integer i,j
    do i=1,d
        do j=1,d
            dot(i,j)=doubledot(v(i,j,:,:),w)
        end do
    end do
end function quadrupledot
!-----------------------------------------
function quadrupledotR(v,w) result(dot)
    real(8),dimension(d,d) :: dot
    real(8),intent(in),dimension(d,d) :: v
    real(8),intent(in),dimension(d,d,d,d) :: w
    integer i,j
    do i=1,d
        do j=1,d
            dot(i,j)=doubledot(v,w(:,:,i,j))
        end do
    end do
end function quadrupledotR
!*****************************************************************************************************************
function doubledot2(v,w) result(idot)
    real(8) :: idot
    real(8),dimension(d) :: temp
    real(8),intent(in),dimension(d,d) :: v,w
    integer i
    do i=1,d
        temp(i)=sum(v(i,:)*w(:,i))
    end do
    idot=sum(temp)
end function doubledot2
!-------------------------------------------
function tripledot2(v,w) result(idot)
    real(8),dimension(d) :: idot
    real(8),dimension(d,d) :: temp
    real(8),intent(in),dimension(d,d) ::w
    real(8),intent(in),dimension(d,d,d) ::v
    integer i,j
    do i=1,d
        do j=1,d
            temp(i,j)=sum(v(i,j,:)*w(:,j))
        end do
        idot(i)=sum(temp(i,:))
    end do
end function tripledot2
!-----------------------------------------
function quadrupledot2(v,w) result(idot)
    real(8),dimension(d,d) :: idot
    real(8),intent(in),dimension(d,d,d,d) :: v
    real(8),intent(in),dimension(d,d) :: w
    integer i,j
    do i=1,d
        do j=1,d
            idot(i,j)=doubledot2(v(i,j,:,:),w)
        end do
    end do
end function quadrupledot2
!------------------------------------------
function doubledot3(v,w) result(odot)
    real(8),dimension(d) :: odot
    real(8),intent(in),dimension(d,d,d) :: v
    real(8),intent(in),dimension(d,d) :: w
    integer i,j
    odot=0
    do j=1,d
        do i=1,d
            odot(j)=odot(j)+sum(v(i,j,:)*w(i,:))
        end do
    end do
end function doubledot3
!-------------------------------------------
   function myequal(v,w) result(eq)
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
!-------------------------------------------
!-------------------------------------------
   function myequal0(v,w) result(eq)
     real(8), intent(in)::  v,w
     logical eq
        if (abs(v-w) .lt. 1d-5) then
           eq=.true.
        else
           eq=.false.
        endif
!    endif
   end function myequal0
!-------------------------------------------
!-------------------------------------------
   function myequalvector(v,w) result(eq)
     type(vector), intent(in)::  v,w
     logical eq !, intent(out) :: eq
     INTEGER :: i
     eq=.true.
     do i=1,d
        if((v%component(i)).myeq.(w%component(i))) then
            CYCLE
        ELSE
            eq=.false.
            EXIT
        ENDIF
     end do
     !if ( (v%x .myeq. w%x) .and. (v%y .myeq. w%y) .and. (v%z .myeq. w%z) ) then
     !      eq=.true.
     !else
     !      eq=.false.
     !endif
   end function myequalvector
!-------------------------------------------
 function myequal0vector(v,w) result(eq)
     type(vector), intent(in)::  v,w
     logical eq
     INTEGER :: i
     eq=.true.
     do i=1,d
        if((v%component(i)).myeqz.(w%component(i))) then
            CYCLE
        ELSE
            eq=.false.
            EXIT
        ENDIF
     end do
   end function myequal0vector
!-------------------------------------------
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
!-----------------------------------------
  function myequal0array(v,w) result(eq)
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
           if ( .not. myequal0(v(n),w(n)) ) then
              eq = .false.
              exit loop
           endif
        enddo loop
     endif
   end function myequal0array
!-----------------------------------------
   function lengthv(v) result(l)
     real(8) l
     type(vector), intent(in) :: v
     integer :: i

     l=0d0
     do i=1,d
        l=l+v%component(i)**2
     end do
     l= sqrt(l)
     !l = sqrt(v%x*v%x+v%y*v%y+v%z*v%z)
   end function lengthv
!-----------------------------------
   !in Fortran 2008, can use NORM2(v) directly
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
!------------------inverse matrix----------------------
 subroutine invers_r(a,b,n)
      implicit none
      integer n,imax,k,j,i,ii
      real(8) ::  amax
      real(8) a(n,n),b(n,n),atmp,btmp,amult,div

      b=0d0
      do i=1,n
         b(i,i)=1d0
      enddo

  MAIN_LOOP: do k=1,n

         amax=-1.d-31
         imax=1
         if (k .ne. n) then

           do i=k,n
             if (dabs(a(i,k)) .ge. amax) then
               amax=dabs(a(i,k))
               imax=i
             endif
           enddo
           if (imax.ne.k) then
             do j=1,n
               atmp=a(imax,j)
               a(imax,j)=a(k,j)
               a(k,j)=atmp
               btmp=b(imax,j)
               b(imax,j)=b(k,j)
               b(k,j)=btmp
             enddo
           endif

         endif

         div=a(k,k)

         if(dabs(div).gt.1.d-14) then

           do ii=1,n
             a(k,ii)=a(k,ii)/div
             b(k,ii)=b(k,ii)/div
           enddo
!          a(k,:)=a(k,:)/div
!          b(k,:)=b(k,:)/div

         endif

         do i=1,n
            amult=a(i,k)
            if (i.eq.k) cycle
            do j=1,n
               a(i,j)=a(i,j)-amult*a(k,j)
               b(i,j)=b(i,j)-amult*b(k,j)
            enddo
         enddo

   enddo MAIN_LOOP

 end subroutine invers_r
 !-----------------------------------------------------
 subroutine invers(a,b,n)
      implicit none
      integer n,imax,k,j,i,ii
      real(8) ::  amax
      complex(8) a(n,n),b(n,n),atmp,btmp,amult,div

      b=0d0
      do i=1,n
         b(i,i)=1d0
      enddo

  MAIN_LOOP: do k=1,n

         amax=-1.d-31
         imax=1
         if (k .ne. n) then

           do i=k,n
             if (cdabs(a(i,k)) .ge. amax) then
               amax=cdabs(a(i,k))
               imax=i
             endif
           enddo
           if (imax.ne.k) then
             do j=1,n
               atmp=a(imax,j)
               a(imax,j)=a(k,j)
               a(k,j)=atmp
               btmp=b(imax,j)
               b(imax,j)=b(k,j)
               b(k,j)=btmp
             enddo
           endif

         endif

         div=a(k,k)

         if(cdabs(div).gt.1.d-14) then

           do ii=1,n
             a(k,ii)=a(k,ii)/div
             b(k,ii)=b(k,ii)/div
           enddo
!          a(k,:)=a(k,:)/div
!          b(k,:)=b(k,:)/div

         endif

         do i=1,n
            amult=a(i,k)
            if (i.eq.k) cycle
            do j=1,n
               a(i,j)=a(i,j)-amult*a(k,j)
               b(i,j)=b(i,j)-amult*b(k,j)
            enddo
         enddo

   enddo MAIN_LOOP

 end subroutine invers
 !-----------------------------------------------------
 FUNCTION voigtMap(input) RESULT(output)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: input
    INTEGER :: output
    IF(input.eq.1) output=1 !xx
    IF(input.eq.2 .OR. input.eq.4) output=6 !xy
    IF(input.eq.3 .OR. input.eq.7) output=5 !xz
    IF(input.eq.5) output=2 !yy
    IF(input.eq.6 .OR. input.eq.8) output=4 !yz
    IF(input.eq.9) output=3 !zz
END FUNCTION voigtMap
 !-----------------------------------------------------

 FUNCTION inv_voigtMap(input) RESULT(output)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: input
    INTEGER,DIMENSION(2) :: output
    SELECTCASE(input)
        CASE(1)
            output(1)=1;output(2)=1 !xx
        CASE(2)
            output(1)=2;output(2)=2 !yy
        CASE(3)
            output(1)=3;output(2)=3 !zz
        CASE(4)
            output(1)=2;output(2)=3 !yz or zy
        CASE(5)
            output(1)=1;output(2)=3 !xz or zx
        CASE(6)
            output(1)=1;output(2)=2 !xy or yx
    ENDSELECT
 END FUNCTION inv_voigtMap
 !-----------------------------------------------------
 end module geometry
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
module ios
 use geometry
 integer, parameter:: ulog=30,uposcar=10, utraj=40,ufco=20,  &
&         umap=60,umatrx=50,utimes=70,ufc1=21,ufc2=22,ufc3=23,ufc4=24,   &
&         ufc=80,ufit1=31,ufit2=32,ufit3=33,ufit4=34,ucor=93,uborn=12

!contradictory with it in io2
!integer, parameter :: uparams = 11
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
4 format(a,99(1x,g12.6))
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
4 format(a,99(1x,g12.6))
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

  implicit none
  character*(*) string
  integer unit,l
  type(vector) var

  l=len_trim(string)
  write(unit,4)string(1:l)//' is=',var
4 format(a,3(1x,g12.6))
end subroutine write_outv
!-----------------------------
subroutine ustring(m,lineout,nrank,iatom,ixyz)

    implicit none
      integer i,j,k,m,n,nrank,iatom(nrank),ixyz(nrank)
      character lineout*80,xyz(3)
      data xyz/'x','y','z'/
!     write(*,*)' Entering ustring with rank=',nrank
      lineout(m+1:m+1)='d'
      m=m+1
      write(lineout(m+1:m+1),'(i1)')nrank
      m=m+1
      lineout(m+1:m+2)='U/'
      m=m+2
      do i=1,nrank
        lineout(m+1:m+2)='d'//xyz(ixyz(i))
        m=m+2
        n=iatom(i)
        if(n.lt.10)then
          write(lineout(m+1:m+1),'(i1)')n
          m=m+1
        else if(n.lt.100)then
          write(lineout(m+1:m+2),'(i2)')n
          m=m+2
        else if(n.lt.1000)then
          write(lineout(m+1:m+3),'(i3)')n
          m=m+3
        else
          write(lineout(m+1:m+4),'(i4)')n
          m=m+4
        endif
      enddo
!     write(*,*)' Exiting ustring'

end subroutine ustring
!--------------------
subroutine warn(unt)
 implicit none
 integer unt

 write(unt,*)'********************************************************************'
 write(unt,*)'|                                                                  |'
 write(unt,*)'|                                                                  |'
 write(unt,*)'|      W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!     |'
 write(unt,*)'|      W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!     |'
 write(unt,*)'|      W    W  A    A  R    R  N N  N  II  N N  N  G       !!!     |'
 write(unt,*)'|      W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   !      |'
 write(unt,*)'|      WW  WW  A    A  R   R   N   NN  II  N   NN  G    G          |'
 write(unt,*)'|      W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!     |'
 write(unt,*)'|                                                                  |'
 write(unt,*)'|   There MIGHT be FC cancellation in the sum and perhaps errors   |'
 write(unt,*)'|                                                                  |'
 write(unt,*)'********************************************************************'
 end subroutine warn

end module ios
!===========================================================
 module atoms_force_constants
! the type atom_id0 concerns the properties of the primitive unit cell
! it s atoms, masses, neighbors, shells and force constants, whereas the
! type atomic refers to atoms in the supercell, their displacements and forces.
 use ios
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
    INTEGER at_type ! I add this because it's used in subroutine identify_atoms_in_supercell
    real(8) mass
 end type
!-------------------------
 type atom_id0   ! all characteristics of an atom in the primitive central cell
    integer at_type
    character name*2
    real(8) mass,charge
    integer nshells       ! how many shells are included=actual dim(shell)
    type(vector) equilibrium_pos !BE CAREFUL WITH vector TYPE
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
! atompos0 is using reduced unit of conventional lattice vector
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
IF(.NOT.ALLOCATED(mas)) ALLOCATE(mas(n))
IF(.NOT.ALLOCATED(natom)) ALLOCATE(natom(n))
IF(.NOT.ALLOCATED(atname)) ALLOCATE(atname(n))
    !allocate( mas(n) ,natom(n), atname(n) )
 end subroutine allocate_mass

 subroutine allocate_primcell(n0)
 integer n0
!print*,'n0=',n0
IF(.NOT.ALLOCATED(atom0)) ALLOCATE(atom0(n0))
IF(.NOT.ALLOCATED(atompos0)) ALLOCATE(atompos0(3,n0))
IF(.NOT.ALLOCATED(atom_type)) ALLOCATE(atom_type(n0))
!allocate( atom0(n0),atompos0(3,n0), atom_type(n0) )  !label_of_primitive_atoms(n0),
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
!******************redudant module use declaration************************
 !use force_constants_module
 !use params
 !use lattice
 !use ios
 !use geometry
 implicit none
! integer, parameter :: max=500
 integer i0,j,shel_count,counter,msort(maxatoms),l
 real(8) dist(maxatoms),d_old,d_min,rmax

! allocate( atom0(1:natoms0)%shells(maxshells) )
 rmax = 0 ; dist = 1d10
 do i0=1,natoms0
    !MARK MODIFIED, bad modification, revised
    IF(ALLOCATED(atom0(i0)%shells)) DEALLOCATE(atom0(i0)%shells)
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
       atom0(i0)%shells(shel_count)%neighbors(counter)%tau = iatomcell0(msort(j))!atom type
       atom0(i0)%shells(shel_count)%neighbors(counter)%n   = iatomcell(:,msort(j))!n1,n2,n3
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
    atom0(i0)%equilibrium_pos%component(1) = atompos(1,i0)
    atom0(i0)%equilibrium_pos%component(2) = atompos(2,i0)
    atom0(i0)%equilibrium_pos%component(3) = atompos(3,i0)
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
!============================================================
!===========================================================
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
    integer ngr, ntotind, ntot
 end type fulldmatrix

 integer maxrank,maxgroups,itrans,irot,ihuang,enforce_inv
 integer nterms(4),maxterms(4),maxtermsindep(4),ngroups(4)
 integer ndindp(4),ndfull(4)
 parameter(maxrank=4)
 integer, allocatable :: nterm(:),ntermsindep(:)
 integer, allocatable :: iatomtermindep(:,:,:),ixyztermindep(:,:,:)
 integer, allocatable :: iatmtrm(:,:,:),ixyztrm(:,:,:)
 real(8), allocatable :: mapmat(:,:,:)
 !real(8) svdcut,radius(4) !svdcut is already declared in module [params]
 real(8) radius(4)
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
 integer inv_constraints,force_constraints,dim_al,dim_ac,n_indep,newdim_al, dim_hom    &
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

