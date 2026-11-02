module ios
 use constants, only : r15
 integer, parameter:: uposcar=10,uparams=11,uborn=12,usym=29, &
&         ufco=20,ulog=30,utraj=40,umatrx=50,umap=60,utimes=70,ufc=80, &
&         ufc1=21,ufc2=22,ufc3=23,ufc4=24,  &
&         ufit1=31,ufit2=32,ufit3=33,ufit4=34,  &
&         ujunk=79,uibz=80,uibs=81,ugrun=82,uband=83,ufbz=84,ucor=93,utherm=110,udos=111


  interface write_out
     module procedure write_outim, write_outrm, write_outiv, write_outrv, &
&    write_outv, write_outi, write_outr, write_outcv, write_outrm3,  &
&    write_outc, write_outcm 
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
  subroutine write_outcm(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer n,m,unt,l,i,j
  complex(r15), dimension(:,:), intent(in) :: var

  l=len_trim(string)
  n=size(var,1) ; m=size(var,2)
! write(unit,*)' write_outrm called;nl,nc=',n,m
  write(unt,*)string(1:l)//' is='
  do i=1,n
     write(unt,4)var(i,:)
  enddo
4 format(99(1x,f13.5))
  end subroutine write_outcm
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
  real(r15), intent(in) :: var(:)

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var(:)
4 format(a,99(1x,f13.5))
  end subroutine write_outrv
!-----------------------------
  subroutine write_outcv(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,l,i
  complex(r15), intent(in) :: var(:)

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var(:)
4 format(a,99(3x,f13.5,1x,f13.5))
  end subroutine write_outcv
!-----------------------------
  subroutine write_outiv(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,i,l
  integer, intent(in):: var(:)

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var(:)
4 format(a,99(1x,i8))

end subroutine write_outiv
!-----------------------------
  subroutine write_outc(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,l
  complex(r15), intent(in):: var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var
4 format(a,99(1x,f13.5))
end subroutine write_outc
!-----------------------------
  subroutine write_outr(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,l
  real(r15), intent(in):: var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var
4 format(a,99(1x,f13.5))
end subroutine write_outr
!-----------------------------
  subroutine write_outi(unt,string,var)
  implicit none
  character*(*), intent(in) :: string
  integer unt,l
  integer, intent(in):: var

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
  type(vector), intent(in):: var

  l=len_trim(string)
  write(unt,4)string(1:l)//' is=',var%x,var%y,var%z
4 format(a,3(1x,f13.5))
end subroutine write_outv

end module ios
!===========================================================
