! generates from the coordinates of the atoms in the unit cell
! all the coordinates in a (nx,ny,nz) lattice
 module cell
 implicit none
 real(8) :: r1(3),r2(3),r3(3),g1(2),g2(3),g3(3),volume
 end module cell
!-------------------------------------------------
 module coord
 implicit none
 real(8) , allocatable :: x(:),y(:),z(:),mass(:)
 integer , allocatable :: typ(:)
 integer :: np
 character coordtype*1

 contains 
  
    subroutine allocate_coord(n)
	implicit none
	integer :: n

	allocate(x(n),y(n),z(n),typ(n),mass(n))
	
	end subroutine allocate_coord
	 
 end module coord
!-------------------------------------------------
 module params
 implicit none
 real(8) :: scale
 character :: coord_type*1 
 integer :: n1,n2,n3
 end module params
!-------------------------------------------------
 program lattice
! this program generates a lattice and atoms inside it
 use cell
 use params
 use coord
 implicit none
 
 call read_input 

 call read_coord

 call write_coord

 end
!----------------------------------
 subroutine read_input
 use cell
 use params
 use coord
 implicit none
 real(8) :: cross(3)

 open(10,file='lattice.data',status='old')
! these are the translation vectors of the primitive cell in cartesian
 read(10,*) r1
 read(10,*) r2
 read(10,*) r3
 read(10,*) scale
 read(10,*) n1,n2,n3
 
 call cross_product(r1,r2,cross)

 call dot_product(r3,cross,volume)

 volume = abs(volume)

 print*,' volume of the prim cell is=',volume

 end subroutine read_input
!----------------------------------
 subroutine read_coord
 use cell
 use params
 use coord
 implicit none
 integer :: i

! open(20,file='coordinates.data',status='old')
! these are the atoms within the primitive cell in cartesian coordinates scaled
 read(10,'(a)') coordtype

 read(10,*) np
 
 call allocate_coord(np)

 do i=1,np
    read(10,*) x(i),y(i),z(i),typ(i),mass(i)
 enddo

 close(10)

 end subroutine read_coord
!----------------------------------
 subroutine cross_product(a,b,c)
 implicit none
 real(8) :: a(3),b(3),c(3)

 c(1) = a(2)*b(3)-a(3)*b(2)
 c(2) = a(3)*b(1)-a(1)*b(3)
 c(3) = a(1)*b(2)-a(2)*b(1)

 end subroutine cross_product
!----------------------------------
 subroutine dot_product(a,b,c)
 implicit none
 real(8) :: a(3),b(3),c

 c = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

 end subroutine dot_product
!----------------------------------
 subroutine write_coord
 use cell
 use params
 use coord
 implicit none
 integer :: i,j,k,p,one

 one =1
 open(30,file='coordinates.out',status='unknown')

 write(30,1)'# np,n1,n2,n3,bx,by,bz=',np*n1*n2*n3,n1,n2,n3, &
 &         scale*n1,scale*n2,scale*n3
1 format(a,4(1x,i5),3(1x,f10.5)) 
 do i=0,n1-1
 do j=0,n2-1
 do k=0,n3-1
   do p=1,np
     if (coordtype .eq. 'c' .or. coordtype .eq. 'C' ) then 
    write(30,3) ( x(p)+i*r1(1)+j*r2(1)+k*r3(1) ) * scale , &
 &	            ( y(p)+i*r1(2)+j*r2(2)+k*r3(2) ) * scale , &
 &              ( z(p)+i*r1(3)+j*r2(3)+k*r3(3) ) * scale , &
 &              typ(p),mass(p)
     elseif (coordtype .eq. 'd' .or. coordtype .eq. 'D' ) then 
    write(30,3) ( x(p)+i )/n1  , &
 &	            ( y(p)+j )/n2  , &
 &              ( z(p)+k )/n3  , &
 &              typ(p),mass(p)
     else
        print*,' coordtype must either be c or d not ',coordtype
        pause
     endif
   enddo
 enddo
 enddo
 enddo

 close(30)
3 format(3(1x,f15.5),i9,2x,f15.6)

 end subroutine write_coord