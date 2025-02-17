!===================================================================
module tetrahedron

        implicit none
        type point
                real(8) :: w,c,d  
! w=argument of the delta function at that point
! c=resulting weight at that point given the coefficient of the delta function 
! d=? not used for the moment ; perhaps weights for the Hilbert transform
                integer :: i   ! index of the assigned kpoint
        end type point

        type tetra
                type(point) :: p(4)
        end type tetra

        integer nk_tet
        real(8), allocatable :: afunc_dos_tet(:),dos_tet(:),eig(:),afunc(:)
        type(tetra), allocatable :: tet(:)

contains
!-------------------------------------
 subroutine allocate_tetra(n1,n2,n3,mesh)
    implicit none
    integer n1, n2, n3 , mesh

    allocate( afunc_dos_tet(mesh),dos_tet(mesh) )
    allocate( tet(n1*n2*n3*6),eig(n1*n2*n3),afunc(n1*n2*n3) )
 end subroutine allocate_tetra
!-------------------------------------
 subroutine deallocate_tetra
    deallocate(tet,eig,afunc,afunc_dos_tet,dos_tet)
 end subroutine deallocate_tetra
!-------------------------------------
 subroutine make_kp_reg_tet
! generate a regular mesh for tetrahedron method (same as make_kp_reg)
! kpoints kpc(nkc) and their weights wk(nkc) are accessed through module kpoints,  otherwise all
! other variables needed are internal to the module tetrahedron
! nk_tet is a copy of nkc which is defined in modeule tetrahedron
! there is no shift in this first call
    use geometry
    use lattice , only : g01,g02,g03
    use kpoints

    implicit none
    integer :: i,j,k,l,nk,nt,np,ierr,pt(8),mc(3)
    real(8) q(3),ep(3),sh(3) 
    real(8) q1(3),q2(3),q3(3),q4(3),q5(3),q6(3),q7(3),q8(3)
    integer inside , i1, j1, k1, nk1
    integer inside2, i2, j2, k2, nk2

    open(126,file='KPOINT_tet.mp',status='unknown')
    write(126,*)'# i ,kp(:,i), wk(i) , kp_reduced'
!   open(127,file='tet.mp',status='unknown')
!   write(127,*)'# i ,tet(i)%p(1...4)%i'
    ep = 0d0  ! shift by a small amount so that only one boundary K remains in the FBZ after folding
    sh = shift(1)/nc(1)*v2a(g01) + shift(2)/nc(2)*v2a(g02) + shift(3)/nc(3)*v2a(g03)
    mc=nc !mc(1)=nc(1); mc(2)=nc(2); mc(3)=nc(3)

!call mypause('entering make_kp_reg_tet',shft(1))

    nkc=nc(1)*nc(2)*nc(3)
    nk_tet=nkc
    nk = 0
    do i = 1,nc(1)
    do j = 1,nc(2)
    do k = 1,nc(3)
       nk = nk+1                      !makes kpc(3,nk) mesh in reciprocal space
!      q = ((i-1+shftx)/nc(1))*g1 + ((j-1+shfty)/nc(2))*g2 + ((k-1+shftz)/nc(3))*g3 
! get_k_info works with the non-shifted kpoints
       q = ((i-1d0)/nc(1))*g01 + ((j-1d0)/nc(2))*g02 + ((k-1d0)/nc(3))*g03  
       kpc(:,nk) = q  + sh !shift works for tetra but IBZ gets much larger
       wk(nk)=1d0/nkc              !1 divided by number of points
       write(126,2)nk,kpc(:,nk), wk(nk), cart2red_g(kpc(:,nk))
    enddo
    enddo
    enddo

! assigning kpc number to each corner of tetrahedra
    np = 0
    nt = 0
    do i = 1,nc(1)
    do j = 1,nc(2)
    do k = 1,nc(3)
       np = np+1
!       if (k.eq.nc(3) .or. j.eq.nc(2) .or. i.eq.nc(1)) cycle
!  q1 = ((i-1+shftx)/nc(1))*g1 + ((j-1+shfty)/nc(2))*g2 + ((k-1+shftz)/nc(3))*g3 + shft
          q1 = kpc(:,np) - sh
          q2 = q1 + (1d0/nc(1))*g01
          q3 = q1 + (1d0/nc(2))*g02
          q4 = q1 + (1d0/nc(1))*g01 + (1d0/nc(2))*g02
          q5 = q1 + (1d0/nc(3))*g03
          q6 = q1 + (1d0/nc(3))*g03 + (1d0/nc(1))*g01
          q7 = q1 + (1d0/nc(3))*g03 + (1d0/nc(2))*g02
          q8 = q1 + (1d0/nc(3))*g03 + (1d0/nc(1))*g01 + (1d0/nc(2))*g02


! regardless of the shift included in q, we can get the kpoint indices
          call get_k_info(q1,mc,nk1,i1,j1,k1,inside)
          if (nk1.ne.np) then
             write(*,4)'i,np,nk1,q=',i,np,nk1,q1
             stop
          endif
          pt(1)=nk1
          call get_k_info(q2,mc,nk1,i1,j1,k1,inside)
!         call get_k_info2(q2,mc,nk2,i2,j2,k2,inside2) as a test of get_k_info
!         if (i1.ne.i2 .or. j1.ne.j2 .or. k1.ne.k2 .or. nk1.ne.nk2 ) then !.or. inside.ne.inside2
!            write(*,'(3(1x,f9.4))') reduce_g(q2)
!            write(*,*)'i1,i2=',i1,i2
!            write(*,*)'j1,j2=',j1,j2
!            write(*,*)'k1,k2=',k1,k2
!            write(*,*)'nk1,nk2=',nk1,nk2
!             write(*,*)'ins1,ins2=',inside,inside2
!            stop
!         endif 
          pt(2)=nk1
          call get_k_info(q3,mc,nk1,i1,j1,k1,inside)
          pt(3)=nk1
          call get_k_info(q4,mc,nk1,i1,j1,k1,inside)
          pt(4)=nk1
          call get_k_info(q5,mc,nk1,i1,j1,k1,inside)
          pt(5)=nk1
          call get_k_info(q6,mc,nk1,i1,j1,k1,inside)
          pt(6)=nk1
          call get_k_info(q7,mc,nk1,i1,j1,k1,inside)
          pt(7)=nk1
          call get_k_info(q8,mc,nk1,i1,j1,k1,inside)
          pt(8)=nk1
! could also use the more time-consuming subroutine get_k_info3(q,nk1,exists)  which loops over kpc(:,1:nkc)

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

       do l = 1,6
          nt = nt+1

          if (l.eq.1) then               ! 1,2,3,6   corner points of each tetrahedron
             tet(nt)%p(1)%i = pt(1)
             tet(nt)%p(2)%i = pt(2)
             tet(nt)%p(3)%i = pt(3)
             tet(nt)%p(4)%i = pt(6)
          elseif (l.eq.2) then            ! 2,3,4,6
             tet(nt)%p(1)%i = pt(2)
             tet(nt)%p(2)%i = pt(3)
             tet(nt)%p(3)%i = pt(4)
             tet(nt)%p(4)%i = pt(6)
          elseif (l.eq.3) then            ! 1,3,5,6
             tet(nt)%p(1)%i = pt(1)
             tet(nt)%p(2)%i = pt(3)
             tet(nt)%p(3)%i = pt(5)
             tet(nt)%p(4)%i = pt(6)
          elseif (l.eq.4) then            ! 3,4,6,8
             tet(nt)%p(1)%i = pt(3)
             tet(nt)%p(2)%i = pt(4)
             tet(nt)%p(3)%i = pt(6)
             tet(nt)%p(4)%i = pt(8)
          elseif (l.eq.5) then            ! 3,5,6,7
             tet(nt)%p(1)%i = pt(3)
             tet(nt)%p(2)%i = pt(5)
             tet(nt)%p(3)%i = pt(6)
             tet(nt)%p(4)%i = pt(7)
          else                            ! 3,6,7,8
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

!   if (mod(nc(1),2).eq.0 .and. mod(nc(2),2).eq.0 .and. mod(nc(3),2).eq.0 ) then
!      do i=1,nk
!         kpc(:,i) = kpc(:,i)+shft(:)
!      enddo
!   endif

!   do i=1,nk
!      write(126,2)i,kpc(:,i),wk(i)
!   enddo
!   do i=1,nt
!      write(127,1)i,(tet(i)%p(l)%i,l=1,4)
!   enddo

 1  format(9(1x,i6))
 2  format(i7,2x,3(1x,f10.5),5x,f9.5,3x,9(f11.5,1x))
 3  format(3(i3),2x,i6,2x,3(1x,f12.5),5x,f9.5)
 4  format(a,3(1x,i6),2x,9(1x,f12.5),5x,f9.5)
    close(126)
!   close(127)

  end subroutine make_kp_reg_tet
!-------------------------------------
 subroutine eigen_tet(nkp,eig)
! assign eigenvalues to each corner of tetrahedra

 implicit none
 integer j,k,nkp,n
 real(8), intent(in) :: eig(nkp)

 do j=1,nkp*6                 !iterate over all tetrahedra
    do k=1,4                  !iterate over four corners of tetrahedron
       n = tet(j)%p(k)%i      !label for kpt(3,i)
       tet(j)%p(k)%w=eig(n)   !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
    enddo
 enddo

 end subroutine eigen_tet
!-------------------------------------
 subroutine weight_tet(energy)
!this subroutine allocates the weights tet(i)%p(j)%c to each corner of the tetrahedron
! energy is the other argument of the delta function. tet is the array of tetrahedra
 implicit none
 integer i,j,k,l,m,nb,uns(4),zero
 real :: e1,e2,e3,e4,e21,e31,e41,e32,e42,e43,en(4),en1(4),a1,a2,a3,a4,b11,b12,b21,b22,b31,b32,b41,b42
 real :: c11,c12,c21,c22,c31,c32,c41,c42
 real(8), intent(IN) :: energy
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

 do i=1,nk_tet*6    ! allocate weights for imaginary terms
!   do j=1,nb
      do k=1,4
         en(k)=tet(i)%p(k)%w         !takes four corner energies and puts them in an array
      enddo
      en1=en
      call ssort(en,uns,4)           !sorts the four corner energies from highest(1) to lowest(4)
!      do k=1,4                 !this do loop assigns the original corner to the newly sorted en(i) array.
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


      if (energy.lt.e1 .or. energy.gt.e4) then
          tet(i)%p(uns(4))%c = 0     !p(uns(4)) is the point with the lowest energy
          tet(i)%p(uns(3))%c = 0     !p(uns(3)) is the point with the second lowest energy
          tet(i)%p(uns(2))%c = 0     !p(uns(2)) is the point with the second highest energy
          tet(i)%p(uns(1))%c = 0     !p(uns(1)) is the point with the highest energy
      elseif (e1.le.energy .and. energy.le.e2) then
          tet(i)%p(uns(4))%c = ((e2-energy)/e21+(e3-energy)/e31+(e4-energy)/e41)*(energy-e1)**2/(e41*e31*e21)
          tet(i)%p(uns(3))%c = (energy-e1)**3/(e21**2*e31*e41)
          tet(i)%p(uns(2))%c = (energy-e1)**3/(e21*e31**2*e41)
          tet(i)%p(uns(1))%c = (energy-e1)**3/(e21*e31*e41**2)
      elseif (e2.le.energy .and. energy.le.e3) then
          c11 = (e3-energy)/e31**2
          c12 = (e4-energy)/e41**2
          c21 = (e3-energy)/e32**2
          c22 = (e4-energy)/e42**2
          c31 = (energy-e2)/e32**2
          c32 = (energy-e1)/e31**2
          c41 = (energy-e2)/e42**2
          c42 = (energy-e1)/e41**2
          b11 = (e3-energy)*(energy-e2)/(e42*e32)+(e4-energy)*(energy-e1)/(e41*e42)+(e3-energy)*(energy-e1)/(e32*e41)
          b12 = (e4-energy)*(energy-e1)/(e42*e31)+(e4-energy)*(energy-e2)/(e42*e32)+(e3-energy)*(energy-e1)/(e31*e32)
          b21 = (e3-energy)*(energy-e2)/(e42*e31)+(e4-energy)*(energy-e2)/(e42*e41)+(e3-energy)*(energy-e1)/(e31*e41)
          b22 = (e3-energy)*(energy-e2)/(e32*e31)+(e4-energy)*(energy-e1)/(e41*e31)+(e4-energy)*(energy-e2)/(e32*e41)
          b31 = (e3-energy)*(energy-e2)/(e42*e31)+(e4-energy)*(energy-e2)/(e42*e41)+(e3-energy)*(energy-e1)/(e31*e41)
          b32 = (e3-energy)*(energy-e2)/(e42*e32)+(e4-energy)*(energy-e1)/(e41*e42)+(e3-energy)*(energy-e1)/(e32*e41)
          b41 = (e3-energy)*(energy-e2)/(e32*e31)+(e4-energy)*(energy-e1)/(e41*e31)+(e4-energy)*(energy-e2)/(e32*e41)
          b42 = (e4-energy)*(energy-e1)/(e42*e31)+(e4-energy)*(energy-e2)/(e42*e32)+(e3-energy)*(energy-e1)/(e31*e32)
          tet(i)%p(uns(4))%c = .5*(c11*b11+c12*b12)
          tet(i)%p(uns(3))%c = .5*(c21*b21+c22*b22)
          tet(i)%p(uns(2))%c = .5*(c31*b31+c32*b32)
          tet(i)%p(uns(1))%c = .5*(c41*b41+c42*b42)
      else
          tet(i)%p(uns(4))%c = (e4-energy)**3/(e41**2*e42*e43)
          tet(i)%p(uns(3))%c = (e4-energy)**3/(e41*e42**2*e43)
          tet(i)%p(uns(2))%c = (e4-energy)**3/(e41*e42*e43**2)
          tet(i)%p(uns(1))%c = ((energy-e3)/e43+(energy-e2)/e42+(energy-e1)/e41)*(e4-energy)**2/(e41*e42*e43)
      endif
!   enddo
 enddo

 end subroutine weight_tet
!-------------------------------------
 subroutine calc_tet(mesh,emin,emax,nk,omt,kpt,eig1,afunc)
! nk is the #of kpoints, mesh and emin,emax are energy meshing params defining omt.
! kpt and eig and afunc are the arrays of kpoints eigenvalues
! which are arguments of the delta function and the weighting function in front of delta, all of
! which are defined at the kpoints
! output: afunc_dos_tet(E) is the result of the integration method for a function afunc(k)*delta(E-eig(k))
! dos_tet is the result of the integration method for a function afunc(k) = 1
 implicit none
 integer :: i,j,la,l,nk,mesh
 real(8) :: emin,emax,omt(mesh),kpt(3,nk), eig1(nk), afunc(nk)

! initialize calct
 afunc_dos_tet = 0
 dos_tet = 0

! do i=1,nk
!    write(34,4)i,eig1(i)
! enddo

! calculate dos for each band (good if no band crossing or k-mesh very fine)

 call eigen_tet(nk,eig1)    ! assign the eigenvalues to each tetrahedron
 do i=1,mesh           !for each energy omt(i), assigns weights to corners of all tetrahedrons
    call weight_tet(omt(i))  ! k1,tet)

    do j=1,6*nk
       do l=1,4
                dos_tet(i) =       dos_tet(i) + tet(j)%p(l)%c/(nk*6)
          afunc_dos_tet(i) = afunc_dos_tet(i) + tet(j)%p(l)%c/(nk*6)*afunc(tet(j)%p(l)%i)
!            write(33,3)j,la,l,tet(j)%p(l)%i,tet(j)%p(l)%c,tet(j)%p(l)%w
       enddo
    enddo

 enddo

3 format(4(2x,i8),3x,5(1x,g10.4))
4 format(1(2x,i8),3x,55(1x,g10.4))

 end subroutine calc_tet
!-------------------------------------
! subroutine tet_sum(om,nk,arg,func,res,array)
 subroutine tet_sum(om,nk,arg,func,res)
! for given om, it calculates res=sum_k delta(om-arg(k))*func(k)
! nk is the #of kpoints, 
! arg and func are the arrays of kpoints eigenvalues, arguments of the delta 
! function and the weighting function in front of delta, all of
! which are defined at the kpoints
 implicit none
 integer,intent(in) :: nk
 real(8),intent(in) :: om, arg(nk), func(nk)  ! res=sum_k array(k)
 real(8),intent(out) :: res !,array(nk) ! res=sum_k array(k)
 integer :: j,l,k

! assign the eigenvalues to each tetrahedron
 call eigen_tet(nk,arg)    

! assign weights to corners of all tetrahedrons
 call weight_tet(om) 

! do the sum over k using the tetrahedron method
! array=0
 res = 0! this is the dos(om)*func(om)
 do j=1,6*nk ! number of possible tetrahedra
    do l=1,4 ! number of corners for each tetrahedron
       k = tet(j)%p(l)%i  ! this is the k index
!      array(k) = array(k) + tet(j)%p(l)%c/(nk*6)*func(k)
       res = res + tet(j)%p(l)%c/(nk*6)*func(k)
!      write(33,3)j,l,tet(j)%p(l)%i,tet(j)%p(l)%c,tet(j)%p(l)%w
    enddo
 enddo

3 format(3(2x,i8),3x,5(1x,g10.4))
4 format(1(2x,i8),3x,55(1x,g10.4))

 end subroutine tet_sum
!-------------------------------------
SUBROUTINE SSORT (X, IY, N)   !sorting algorithm, (insertion sort)
IMPLICIT NONE

INTEGER N
REAL X(1:N)
INTEGER IY(N)
REAL TEMP
INTEGER I, ISWAP(1), ITEMP, ISWAP1
INTRINSIC MAXLOC

 DO I=1,N-1
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
 enddo

END SUBROUTINE SSORT

end module tetrahedron
!===================================================================
