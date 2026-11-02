! test_point_groups.f90: the 32 crystallographic point groups generated from their generators (Cartesian);
! pg_setup (common/point_group.f90) must give the standard order, number of classes and irrep dimensions,
! and pass its orthogonality checks (ier 0). Expected output: test_point_groups.out. Build and run:
!   gfortran -ffree-line-length-0 -o t ../../../common/constants.f90 ../../../common/point_group.f90 \
!            ../../zhegv.f test_point_groups.f90 && ./t > new.out && diff new.out test_point_groups.out
program t
use constants, only : r15
use point_group_rules
implicit none
real(r15) :: g(3,3,4), ops(3,3,48), s(3,3), c3h(3,3), c6h(3,3), c4(3,3), c3d(3,3), id(3,3), mx(3,3), mz(3,3), c2x(3,3), c2z(3,3), inv(3,3)
real(r15) pi, a
integer ng, nop, i, j, k, ier, n0
logical found, grew
character(len=6) :: names(32)
integer :: ngen(32)
character(len=200) line
pi = 4*atan(1d0)
id = 0; do i=1,3; id(i,i)=1; enddo
inv = -id
mz = id; mz(3,3) = -1
mx = id; mx(1,1) = -1
c2x = -id; c2x(1,1) = 1
c2z = -id; c2z(3,3) = 1
c4 = 0; c4(1,2) = -1; c4(2,1) = 1; c4(3,3) = 1
a = 2*pi/3; c3h = id; c3h(1,1)=cos(a); c3h(1,2)=-sin(a); c3h(2,1)=sin(a); c3h(2,2)=cos(a)
a = 2*pi/6; c6h = id; c6h(1,1)=cos(a); c6h(1,2)=-sin(a); c6h(2,1)=sin(a); c6h(2,2)=cos(a)
c3d = 0; c3d(1,3)=1; c3d(2,1)=1; c3d(3,2)=1
do k = 1, 32
  select case (k)
  case(1); names(k)='C1';  ng=0
  case(2); names(k)='Ci';  ng=1; g(:,:,1)=inv
  case(3); names(k)='C2';  ng=1; g(:,:,1)=c2z
  case(4); names(k)='Cs';  ng=1; g(:,:,1)=mz
  case(5); names(k)='C2h'; ng=2; g(:,:,1)=c2z; g(:,:,2)=inv
  case(6); names(k)='D2';  ng=2; g(:,:,1)=c2z; g(:,:,2)=c2x
  case(7); names(k)='C2v'; ng=2; g(:,:,1)=c2z; g(:,:,2)=mx
  case(8); names(k)='D2h'; ng=3; g(:,:,1)=c2z; g(:,:,2)=c2x; g(:,:,3)=inv
  case(9); names(k)='C4';  ng=1; g(:,:,1)=c4
  case(10);names(k)='S4';  ng=1; g(:,:,1)=matmul(mz,c4)
  case(11);names(k)='C4h'; ng=2; g(:,:,1)=c4; g(:,:,2)=inv
  case(12);names(k)='D4';  ng=2; g(:,:,1)=c4; g(:,:,2)=c2x
  case(13);names(k)='C4v'; ng=2; g(:,:,1)=c4; g(:,:,2)=mx
  case(14);names(k)='D2d'; ng=2; g(:,:,1)=matmul(mz,c4); g(:,:,2)=c2x
  case(15);names(k)='D4h'; ng=3; g(:,:,1)=c4; g(:,:,2)=c2x; g(:,:,3)=inv
  case(16);names(k)='C3';  ng=1; g(:,:,1)=c3h
  case(17);names(k)='S6';  ng=2; g(:,:,1)=c3h; g(:,:,2)=inv
  case(18);names(k)='D3';  ng=2; g(:,:,1)=c3h; g(:,:,2)=c2x
  case(19);names(k)='C3v'; ng=2; g(:,:,1)=c3h; g(:,:,2)=mx
  case(20);names(k)='D3d'; ng=3; g(:,:,1)=c3h; g(:,:,2)=c2x; g(:,:,3)=inv
  case(21);names(k)='C6';  ng=1; g(:,:,1)=c6h
  case(22);names(k)='C3h'; ng=2; g(:,:,1)=c3h; g(:,:,2)=mz
  case(23);names(k)='C6h'; ng=2; g(:,:,1)=c6h; g(:,:,2)=inv
  case(24);names(k)='D6';  ng=2; g(:,:,1)=c6h; g(:,:,2)=c2x
  case(25);names(k)='C6v'; ng=2; g(:,:,1)=c6h; g(:,:,2)=mx
  case(26);names(k)='D3h'; ng=3; g(:,:,1)=c3h; g(:,:,2)=mz; g(:,:,3)=c2x
  case(27);names(k)='D6h'; ng=3; g(:,:,1)=c6h; g(:,:,2)=c2x; g(:,:,3)=inv
  case(28);names(k)='T';   ng=3; g(:,:,1)=c2z; g(:,:,2)=c2x; g(:,:,3)=c3d
  case(29);names(k)='Th';  ng=4; g(:,:,1)=c2z; g(:,:,2)=c2x; g(:,:,3)=c3d; g(:,:,4)=inv
  case(30);names(k)='O';   ng=2; g(:,:,1)=c4; g(:,:,2)=c3d
  case(31);names(k)='Td';  ng=2; g(:,:,1)=matmul(mz,c4); g(:,:,2)=c3d
  case(32);names(k)='Oh';  ng=3; g(:,:,1)=c4; g(:,:,2)=c3d; g(:,:,3)=inv
  end select
  nop = 1; ops(:,:,1) = id
  do i = 1, ng; call add(g(:,:,i)); enddo
  grew = .true.
  do while (grew)
    grew = .false.; n0 = nop
    do i = 1, n0; do j = 1, n0
       call add(matmul(ops(:,:,i), ops(:,:,j)))
    enddo; enddo
    if (nop > n0) grew = .true.
  enddo
  call pg_setup(nop, ops(:,:,1:nop), 6, ier)
  write(line,'(a6,a,i3,a,i3,a,i2,a)') names(k),' order',nop,' classes',pg_ncls,' ier',ier,'  dims'
  do i = 1, pg_nirr; write(line,'(a,1x,i1)') trim(line), pg_dim(i); enddo
  write(*,'(a)') trim(line)
  write(line,'(a)') '        labels:'
  do i = 1, pg_nirr; line = trim(line)//' '//trim(pg_label(i)); enddo
  write(*,'(a)') trim(line)
  if (names(k).eq.'C6v' .or. names(k).eq.'Td' .or. names(k).eq.'Oh' .or. names(k).eq.'C3' .or. names(k).eq.'D3h' .or. names(k).eq.'C6') call pg_write_table(6,'      ')
enddo
contains
 subroutine add(m)
 real(r15), intent(in) :: m(3,3)
 integer l
 do l = 1, nop
   if (maxval(abs(ops(:,:,l)-m)) < 1d-8) return
 enddo
 nop = nop + 1; ops(:,:,nop) = m
 end subroutine add
end program t
