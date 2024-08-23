!===================================================================
 subroutine make_kp_reg_tet(nx,ny,nz,sx,sy,sz,kpt,wkt)
! generate a regular mesh for tetrahedron method (same as make_kp_reg)

    use geometry
    use params
    use tetrahedron
    use lattice
    use ios
    use kpoints , only : nc,nkc

    implicit none
    integer :: i,j,k,l,nx,ny,nz,nk,nt,np,ierr,pt(8)
    real(8) q(3),ep(3),sx,sy,sz,kpt(3,nx*ny*nz),wkt(nx*ny*nz),shft(3)
    real(8) q1(3),q2(3),q3(3),q4(3),q5(3),q6(3),q7(3),q8(3)
    integer inside, i1, j1, k1, nk1

!    type(tetra), allocatable :: tet(:)
    open(126,file='KPOINT_tet.MP',status='unknown')
    write(126,*)'# i ,kp(:,i), wk(i)'
    ep = 0d0  ! shift by a small amount so that only one boundary K remains in the FBZ after folding
    shft = (-0.5d0)*(g1+g2+g3)

!    allocate(tet((nx-1)*(ny-1)*(nz-1)*6), stat=ierr)
!    if (ierr/=0) print*, "tet : Allocation failed"

    nk = 0                                                                  
    do i = 1,nx
    do j = 1,ny
    do k = 1,nz
       nk = nk+1                                                                !makes kpt(3,nk) mesh in reciprocal space
       q = ((i-1+sx)/nx)*g1 + ((j-1+sy)/ny)*g2 + ((k-1+sz)/nz)*g3 + ep
       kpt(:,nk) = q(:) 
       wkt(nk)=1d0/nkc                                                          !1 divided by number of points
      
    enddo
    enddo
    enddo
 
! assigning kpt number to each corner of tetrahedra   
    np = 0
    nt = 0
    do i = 1,nx
    do j = 1,ny
    do k = 1,nz
       np = np+1
!       if (k.eq.nz .or. j.eq.ny .or. i.eq.nx) cycle
       do l = 1,6
          q1 = ((i-1+sx)/nx)*g1 + ((j-1+sy)/ny)*g2 + ((k-1+sz)/nz)*g3 + shft
          q2 = q1 + (1d0/nx)*g1
          q3 = q1 + (1d0/ny)*g2
          q4 = q1 + (1d0/nx)*g1 + (1d0/ny)*g2
          q5 = q1 + (1d0/nz)*g3
          q6 = q1 + (1d0/nz)*g3 + (1d0/nx)*g1
          q7 = q1 + (1d0/nz)*g3 + (1d0/ny)*g2
          q8 = q1 + (1d0/nz)*g3 + (1d0/nx)*g1 + (1d0/ny)*g2

          call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
          pt(1)=nk1
          call get_k_info(q2,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
          pt(2)=nk1
          call get_k_info(q3,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
          pt(3)=nk1
          call get_k_info(q4,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
          pt(4)=nk1
          call get_k_info(q5,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
          pt(5)=nk1
          call get_k_info(q6,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
          pt(6)=nk1
          call get_k_info(q7,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
          pt(7)=nk1
          call get_k_info(q8,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
          pt(8)=nk1

!          write(*,*) 'i,j,k=',i,j,k
!          write(*,*) 'pt=',pt
!          write(*,*) 'q1=',q1
!          write(*,*) 'q2=',q2
!          write(*,*) 'q3=',q3
!          write(*,*) 'q4=',q4
!          write(*,*) 'q5=',q5
!          write(*,*) 'q6=',q6
!          write(*,*) 'q7=',q7
!          write(*,*) 'q8=',q8

! stop

          nt = nt+1
	  
          if (l.eq.1) then   			! 1,2,3,6   corner points of each tetrahedron
             tet(nt)%p(1)%i = pt(1) 
             tet(nt)%p(2)%i = pt(2)
	     tet(nt)%p(3)%i = pt(3)
	     tet(nt)%p(4)%i = pt(6)
          elseif (l.eq.2) then			! 2,3,4,6
             tet(nt)%p(1)%i = pt(2) 
             tet(nt)%p(2)%i = pt(3)
	     tet(nt)%p(3)%i = pt(4)
	     tet(nt)%p(4)%i = pt(6)
	  elseif (l.eq.3) then			! 1,3,5,6
             tet(nt)%p(1)%i = pt(1)
             tet(nt)%p(2)%i = pt(3)
	     tet(nt)%p(3)%i = pt(5)
	     tet(nt)%p(4)%i = pt(6)
          elseif (l.eq.4) then			! 3,4,6,8
             tet(nt)%p(1)%i = pt(3)
             tet(nt)%p(2)%i = pt(4)
	     tet(nt)%p(3)%i = pt(6)
	     tet(nt)%p(4)%i = pt(8)
          elseif (l.eq.5) then			! 3,5,6,7
             tet(nt)%p(1)%i = pt(3)
             tet(nt)%p(2)%i = pt(5)
	     tet(nt)%p(3)%i = pt(6)
	     tet(nt)%p(4)%i = pt(7)
          else           			! 3,6,7,8
             tet(nt)%p(1)%i = pt(3)
             tet(nt)%p(2)%i = pt(6)
	     tet(nt)%p(3)%i = pt(7)
	     tet(nt)%p(4)%i = pt(8)
          endif
      enddo
   enddo
   enddo
   enddo

    write(ulog,*)'KP_REG: Number of regular kpoints generated is=',nk

    if (mod(nx,2).eq.0 .and. mod(ny,2).eq.0 .and. mod(nz,2).eq.0 ) then
       do i=1,nk
          kpt(:,i) = kpt(:,i)+shft(:)
       enddo
    endif

    do i=1,nk
       write(126,2)i,kpt(:,i),wkt(i)
    enddo

 2  format(i7,2x,3(1x,f12.5),5x,f9.5)
 3  format(3(i3),2x,i6,2x,3(1x,f12.5),5x,f9.5)
    close(126)

  end subroutine make_kp_reg_tet
!============================================================================
subroutine eigen_tet(nx,ny,nz,la,nkp)
! calculate eigenvalues at each corner of tetrahedra
use constants
use tetrahedron
use eigen
use kpoints

implicit none
integer j,k,la,nx,ny,nz,n,nkp,w

real(8) k_leng, k_leng_min
!real(8), intent(in) :: eival(n,nkp)
!type(tetra), intent(inout) :: tet((nx-1)*(ny-1)*(nz-1)*6)

k_leng_min=1.0

do j=1,(nx)*(ny)*(nz)*6                       !iterate over all tetrahedrons
   do k=1,4                                         !iterate over four corners of tetrahedron
!      do l=1,n                                      !iterate over band
         w = tet(j)%p(k)%i                          !label for kpt(3,i)
         tet(j)%p(k)%w=cnst*sqrt(eigenval(la,w))     !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
         
!         k_leng=sqrt(kpc(1,w)**2.0+kpc(2,w)**2.0+kpc(3,w)**2.0)
!         tet(j)%p(k)%w=100*k_leng
!
!         if (k_leng < k_leng_min .and. k_leng .ne. 0.0) then
!                 k_leng_min=k_leng
!         endif

!!         if (k_leng <= 1d-4) then
!!              tet(j)%p(k)%w(l)=1d-4
!!         endif

!      enddo
   enddo
enddo

!write(*,*) 'kleng_min=',k_leng_min

end subroutine eigen_tet



!============================================
subroutine weight_tet(nx,ny,nz,e)            !this subroutine allocates the weights tet(i)%p(j)%c(k) to each corner of the tetrahedron
use tetrahedron                                     !nx,ny,nz are meshing params.  nb is #ofbands. e is energy to do calc at. tet is the array of tetrahedrons
use constants
implicit none
integer i,j,k,l,m,nb,nx,ny,nz,uns(4),zero
real :: e1,e2,e3,e4,e21,e31,e41,e32,e42,e43,en(4),en1(4),a1,a2,a3,a4,b11,b12,b21,b22,b31,b32,b41,b42
real :: c11,c12,c21,c22,c31,c32,c41,c42
real(8) :: e
!type(tetra), intent(inout) :: tet((nx-1)*(ny-1)*(nz-1)*6)

real(8) small1, small2
small1=1d-10
small2=1d-10


zero=0

!   do i=1,(nx-1)*(ny-1)*(nz-1)*6    ! allocate weights for real terms
!      do j=1,nb
!         do k=1,4
!            en(k)=tet(i)%p(k)%w(j)
!         enddo            
!         do k=1,4
!            a1 = 1      
!            do l=1,4
!               if (l.eq.k) cycle
!               a1 = a1*(en(l)-en(k))
!            enddo
!            a2 = 0
!            do l=1,4
!               if (l.eq.k) cycle
!               a2 = a2 + (e-en(l))/(en(k)-en(l))*log(abs(e-en(k)))
!            enddo
!            a3 = 0
!            do l=1,4
!               if (l.eq.k) cycle
!               a4 = 1
!               do m=1,4
!                  if (m.eq.l) cycle
!                  a4 = a4*(en(m)-en(l))
!               enddo
!               a3 = a3 + (e-en(l))**3/a4*log(abs(e-en(l)))/(en(k)-en(l))
!            enddo
!            tet(i)%p(k)%d(j) = (e-en(k))**2/a1*(1+a2)+a3
!         enddo   
!      enddo   
!   enddo   

do i=1,(nx)*(ny)*(nz)*6    ! allocate weights for imaginary terms
!   do j=1,nb
      do k=1,4
         en(k)=tet(i)%p(k)%w         !takes four corner energies and puts them in an array
         !en(k)=100.15+.1*k
      enddo
      en1=en
      call ssort(en,uns,4)              !sorts the four corner energies from highest to lowest, i.e. en(1) = highest, en(4) = lowest
!      do k=1,4                          !this do loop assigns the original corner to the newly sorted en(i) array. see below.
     if (en(1).ne.en(2) .and. en(2).ne.en(3) .and. en(3).ne.en(4)) then
        do k=1,4
           if (en(k).eq.en1(1)) then
              uns(k) = 1
           elseif (en(k).eq.en1(2)) then
              uns(k) = 2
           elseif (en(k).eq.en1(3)) then
              uns(k) = 3
           elseif (en(k).eq.en1(4)) then
              uns(k) = 4
           endif
        enddo 
      else
         do k=1,4
            if (en(k).eq.en1(1)) then
               uns(k) = 1
            elseif (en(k).eq.en1(2)) then
               uns(k) = 2
            elseif (en(k).eq.en1(3)) then
               uns(k) = 3
            elseif (en(k).eq.en1(4)) then
               uns(k) = 4
            endif
            if (k.ne.1) then
               do l=1,k-1
                  if (uns(l).eq.uns(k) .and. uns(k).eq.1) then
                     if (en(k).eq.en1(2)) then
                        uns(k) = 2
                     elseif (en(k).eq.en1(3)) then
                        uns(k) = 3
                     elseif (en(k).eq.en1(4)) then
                        uns(k) = 4
                     endif
                  elseif (uns(l).eq.uns(k) .and. uns(k).eq.2) then
                     if (en(k).eq.en1(3)) then
                        uns(k) = 3
                     elseif (en(k).eq.en1(4)) then 
                        uns(k) = 4
                     endif
                  elseif (uns(l).eq.uns(k) .and. uns(k).eq.3) then
                     uns(k) = 4
                  endif
               enddo
            endif
         enddo
      endif

 !     enddo
      e1 = en(4)
      e2 = en(3)
      e3 = en(2)
      e4 = en(1)
      e21 = e2-e1
      e31 = e3-e1
      e41 = e4-e1
      e32 = e3-e2
      e42 = e4-e2
      e43 = e4-e3

! from here sy to avoid numerical error
      if (abs(e21) < small1 .and. e21 > 0.0) then
           e21=small2
      endif
      if (abs(e21) < small1 .and. e21 < 0.0) then
           e21=-small2
      endif

      if (abs(e31) < small1 .and. e31 > 0.0) then
           e31=small2
      endif
      if (abs(e31) < small1 .and. e31 < 0.0) then
           e31=-small2
      endif

      if (abs(e41) < small1 .and. e41 > 0.0) then
           e41=small2
      endif
      if (abs(e41) < small1 .and. e41 < 0.0) then
           e41=-small2
      endif

      if (abs(e42) < small1 .and. e42 > 0.0) then
           e42=small2
      endif
      if (abs(e42) < small1 .and. e42 < 0.0) then
           e42=-small2
      endif

      if (abs(e43) < small1 .and. e43 > 0.0) then
           e43=small2
      endif
      if (abs(e43) < small1 .and. e43 < 0.0) then
           e43=-small2
      endif


! to here sy


      if (e.lt.e1 .or. e.gt.e4) then
          tet(i)%p(uns(4))%c = 0     !p(uns(4)) is the point with the lowest energy
          tet(i)%p(uns(3))%c = 0     !p(uns(3)) is the point with the second lowest energy
          tet(i)%p(uns(2))%c = 0     !p(uns(2)) is the point with the second highest energy
          tet(i)%p(uns(1))%c = 0     !p(uns(1)) is the point with the highest energy
      elseif (e1.le.e .and. e.le.e2) then
          tet(i)%p(uns(4))%c = ((e2-e)/e21+(e3-e)/e31+(e4-e)/e41)*(e-e1)**2/(e41*e31*e21)
          tet(i)%p(uns(3))%c = (e-e1)**3/(e21**2*e31*e41)
          tet(i)%p(uns(2))%c = (e-e1)**3/(e21*e31**2*e41)
          tet(i)%p(uns(1))%c = (e-e1)**3/(e21*e31*e41**2)
      elseif (e2.le.e .and. e.le.e3) then
          c11 = (e3-e)/e31**2
          c12 = (e4-e)/e41**2
          c21 = (e3-e)/e32**2
          c22 = (e4-e)/e42**2
          c31 = (e-e2)/e32**2
          c32 = (e-e1)/e31**2
          c41 = (e-e2)/e42**2
          c42 = (e-e1)/e41**2
          b11 = (e3-e)*(e-e2)/(e42*e32)+(e4-e)*(e-e1)/(e41*e42)+(e3-e)*(e-e1)/(e32*e41)
          b12 = (e4-e)*(e-e1)/(e42*e31)+(e4-e)*(e-e2)/(e42*e32)+(e3-e)*(e-e1)/(e31*e32)
          b21 = (e3-e)*(e-e2)/(e42*e31)+(e4-e)*(e-e2)/(e42*e41)+(e3-e)*(e-e1)/(e31*e41)
          b22 = (e3-e)*(e-e2)/(e32*e31)+(e4-e)*(e-e1)/(e41*e31)+(e4-e)*(e-e2)/(e32*e41)
          b31 = (e3-e)*(e-e2)/(e42*e31)+(e4-e)*(e-e2)/(e42*e41)+(e3-e)*(e-e1)/(e31*e41)
          b32 = (e3-e)*(e-e2)/(e42*e32)+(e4-e)*(e-e1)/(e41*e42)+(e3-e)*(e-e1)/(e32*e41)
          b41 = (e3-e)*(e-e2)/(e32*e31)+(e4-e)*(e-e1)/(e41*e31)+(e4-e)*(e-e2)/(e32*e41)
          b42 = (e4-e)*(e-e1)/(e42*e31)+(e4-e)*(e-e2)/(e42*e32)+(e3-e)*(e-e1)/(e31*e32)
          tet(i)%p(uns(4))%c = .5*(c11*b11+c12*b12)
          tet(i)%p(uns(3))%c = .5*(c21*b21+c22*b22)
          tet(i)%p(uns(2))%c = .5*(c31*b31+c32*b32)
          tet(i)%p(uns(1))%c = .5*(c41*b41+c42*b42)
      else
          tet(i)%p(uns(4))%c = (e4-e)**3/(e41**2*e42*e43)
          tet(i)%p(uns(3))%c = (e4-e)**3/(e41*e42**2*e43)
          tet(i)%p(uns(2))%c = (e4-e)**3/(e41*e42*e43**2)
          tet(i)%p(uns(1))%c = ((e-e3)/e43+(e-e2)/e42+(e-e1)/e41)*(e4-e)**2/(e41*e42*e43)
      endif
!   enddo
enddo

end subroutine weight_tet


!=======================================================================================
subroutine calc_tet(nx,ny,nz,ndn,kpt)        
!nx,ny,nz are meshing params. ndn is # of bands
!! calct is the result of the integration method for a function F(k)
!! doswt is the result of the integration method for a function F(k) = 1
use tetrahedron                                                             
use constants                                                               
use om_dos
use geometry
use ios
implicit none                                      

integer :: i,j,k,kk,l,nx,ny,nz,mesh,ndn
real(8) :: kfunc(nx*ny*nz),kpt(3,nx*ny*nz),kq
real(8) integrate_dos, total_dos

open(udos,file='dos_tet.dat',status='unknown')
open(udos+1,file='weightet.dat',status='unknown')

calct = 0
doswt = 0

! kfunc(i)=coefficient of the delta function in sum_k  kfunc(k)*delta(om-om(k))
do i=1,nx*ny*nz                                             
   kfunc(i) = 1  ! 1 gives the dos

!   kq = length(kpt(:,i)) 
!   kfunc(i) = 1/(exp(100*kq)-1+1d-4)
!   kfunc(i) = kq
!   kfunc(i) = 1/(kq**2+1d-4)   !kfunc(i) = 1
enddo

do i=1,wmesh   !for each energy om(i), assigns weights to corners of all tetrahedrons
  do k=1,ndn
     call eigen_tet(nx,ny,nz,k,nx*ny*nz)
     call weight_tet(nx,ny,nz,om(i))

      do j=1,(nx)*(ny)*(nz)*6  
         do l=1,4
            doswt(i,k) = doswt(i,k) + tet(j)%p(l)%c/((nx)*(ny)*(nz)*6)      !dos calculation
            calct(i,k) = calct(i,k) + tet(j)%p(l)%c/((nx)*(ny)*(nz)*6)*kfunc(tet(j)%p(l)%i) 
         enddo
      enddo 
  enddo
enddo

!do i=1,wmesh
!   write(udos,*)om(i),(calct(i,j),j=1,n)
!enddo

integrate_dos=0.0
total_dos=0.0
do i=1,wmesh
   total_dos=sum(doswt(i,:))
   integrate_dos=integrate_dos+total_dos*wmax/wmesh
   write(udos,15) om(i),i,integrate_dos,total_dos,(doswt(i,ndn-j+1),j=1,ndn)
enddo

integrate_dos=0.0
total_dos=0.0
do i=1,wmesh
   total_dos=sum(calct(i,:))
   integrate_dos=integrate_dos+total_dos*wmax/wmesh
   write(udos+1,15) om(i),i,integrate_dos,total_dos,(calct(i,ndn-j+1),j=1,ndn)
enddo

close(udos); close(udos+1)

15 format(99(1x,g10.4))


end subroutine calc_tet
!=======================================================================================
SUBROUTINE SSORT (X, IY, N)   !sorting algorithm, (insertion sort)
IMPLICIT NONE

INTEGER N
REAL X(1:N)
INTEGER IY(N)
REAL TEMP
INTEGER I, ISWAP(1), ITEMP, ISWAP1
INTRINSIC MAXLOC
DO 200 I=1,N-1
ISWAP=MAXLOC(X(I:N))
ISWAP1=ISWAP(1)+I-1
IF(ISWAP1.NE.I) THEN
TEMP=X(I)
X(I)=X(ISWAP1)
X(ISWAP1)=TEMP
ITEMP=IY(I)
IY(I)=IY(ISWAP1)
IY(ISWAP1)=ITEMP
ENDIF
200 CONTINUE
RETURN
END
!===============================================================
