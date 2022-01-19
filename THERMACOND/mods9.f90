!==========================================================
 module io2
 integer, parameter:: ulog=30,ubs=50,uparams=13,udos=40,ueiv=60,ukap=175,uklem=71,uslf=3000, &
& ufc0=20,ufc1=21,ufc2=22,ufc3=23,ufc4=24,utimes=11,ualph=27,uvel=31,ufine=67,umod=2000,  &
& fbz=93,bzf=95,ibz=94,ibs=96,ucors=65,ueivibz=64,uv3=85,ugrun=80,urate=88,ucof=87,  &
& udist=290, ucross=89, k_kappa=9000 
 integer, parameter:: debug=9999, self_detail=9990, iself_bs=9980, ksubset_bs_inp=9970  ! sy. 
 end module io2
!==========================================================
 module om_dos
 integer wmesh,ndyn2
 real(8) wmax,width,etaz
 real(8), allocatable :: dos(:,:),om(:)

   contains

   subroutine set_omdos(la,mesh)
   implicit none
   integer la,mesh,i
   allocate(om(mesh),dos(0:la,mesh))
   ndyn2 = la  ! this is a copy of ndyn
   do i=1,mesh
      om(i) = wmax *(0.0001 + (i-1)/(mesh-1.))
   enddo
   end subroutine set_omdos

   subroutine write_dos
   use io2
   implicit none
   integer i,la
   real(8) sumdos

   write(udos,*)'# om(i),i,integrated_dos,total_dos,(dos(la,i),la=1,n))'
   sumdos=0
   do i=1,wmesh
      sumdos=sumdos+dos(0,i)*wmax/wmesh
      write(udos,3)om(i),i,sumdos,(dos(la,i),la=0,min(10,ndyn2))
   enddo
   sumdos=sumdos-dos(0,1)*wmax/wmesh/2d0
   
3  format(g11.5,2x,i5,99(1x,g10.4))
   end subroutine write_dos

 end module om_dos
!==========================================================
 module eigen
 integer ndyn,nkc2  ! this is a copy of nkc
! eigenval for the coarse mesh (nkc,kpc),

 real(8), allocatable :: velocibz(:,:,:),veloc(:,:,:)
 real(8), allocatable :: eigenval_bs(:,:),eigenval(:,:),eivalibz(:,:)
! sqr of phonon freqs in eV/A/A/m0
 complex(8), allocatable :: eigenvec_bs(:,:,:),eigenvec(:,:,:),grun(:,:),grun_bs(:,:),eivecibz(:,:,:)

    contains

    subroutine allocate_eig(nb,ni,nk) ! nb for band, nk for coarse mesh in FBZ
    integer nb,nk,ni
      allocate( eigenval(nb,nk),eigenvec(nb,nb,nk),veloc   (3,nb,nk), grun(nb,nk) )
      allocate( eivalibz(nb,ni),eivecibz(nb,nb,ni),velocibz(3,nb,ni) )
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
    subroutine deallocate_eig ! nb for band, nk for coarse mesh in FBZ
      deallocate( eigenval,eigenvec,grun,veloc,eivalibz,eivecibz,velocibz )
    end subroutine deallocate_eig

!---------------------------------
 function onedim(nb,nk)
! calculates the running index onedim
! calculates the running index onedim
! use eigen
 implicit none
 integer nk,nb,onedim
 ! onedim=(nb-1)*nkc2+nk !if outermost loop is j=1,ndyn and innermost loop is i=1,nkc2
 onedim=(nk-1)*ndyn+nb   !if innermost loop is j=1,ndyn and outermost loop is i=1,nkc2
 end function onedim
!---------------------------------
 function nkpt(onedm,nb)
! assumes the outerloop is over kpoints and the innerloop over bands
! use eigen
 implicit none
 integer nkpt,onedm,nb
 nkpt = int((onedm-1)/nb)+1
 end function nkpt
!---------------------------------
 function nband(onedm,nb)
! assumes the outerloop is over kpoints and the innerloop over bands
! use eigen
 implicit none
 integer nband,onedm,nb
 nband = mod(onedm-1,nb)+1
 end function nband

 end module eigen
!==========================================================
 module phi3
 complex(8), allocatable :: v3(:,:,:,:,:),selfN(:,:,:) ,selfU(:,:,:),v33(:)
 real(8), allocatable :: v33sq(:),v3sq(:,:,:,:,:) ! the latter used for tetrahedron
 integer, allocatable:: nq1(:),nq2(:),nq3(:),la1(:),la2(:),la3(:)
 integer nv3,nv3split,readv3,iter,  max_iter !,conv_iter ! sy added iter,split, calk, max_iter !, split, calk,writev3, k1 moved conv_iter to integers
 real(8) v3_threshold,const33 
!, conv_error , conv_max_error, conv_max_diff, conv_diff, &
! & conv_diff_kap, conv_max_diff_kap,update ! sy added conv_norm, conv_max, update
 integer ksub_size,usetetra
 character(132) v3path

 contains

  subroutine allocate_v33(n)
    allocate(v33sq(n),nq1(n),nq2(n),nq3(n) ) !,la1(n),la2(n),la3(n) )
  end subroutine allocate_v33
  subroutine deallocate_v33
    deallocate(v33sq,nq1,nq2,nq3 ) !,la1,la2,la3)
  end subroutine deallocate_v33
!  subroutine allocate_phi3(nk,nl,nw,ni)
!  integer nk,nw,nl,ni
!!   allocate(delta3(nk,nl,nw),gamma3(nk,nl,nw),v3(nk,nk,nl,nl,nl))
!!   allocate(selfN(ni,nl,nw),selfU(ni,nl,nw))
!    allocate(v3(ni,nk,nl,nl,nl))
!  end subroutine allocate_phi3

 end module phi3
!==========================================================
 module kpoints
! this includes 4 kpoint meshes: one for the band structure: kp_bs(3,nkp_bs)
! one shifted MP-coarse mesh around gamma for interpolation : kpc(3,nkc)
! a mapping of this coarse mesh into the FBZ: kpc(3,mapbz(1:nbz)) 
! one fine grid in the irreducible FBZ for kappa/dos sums: kibz(3,nibz)
 use lattice   ! it is also used for all included subroutines
 use io2
 use constants
 implicit none
 integer nbs1,nbs2,nbs3,nkp_bs   ! for the band structure
 integer nc1,nc2,nc3,nkc,nbz,nibz_coarse,nkcc  ! coarse mesh,its FBZ folded mesh,IBZ_coarse mesh
 integer nib1,nib2,nib3,nbz_fine  ! for the fine mesh in the Irreducible BZ
 integer nshl,nbkmax,npos           ! max no of neighbors=20
 integer, allocatable :: nb(:,:),nbk(:),mapbz(:) ,mapibz(:),mappos(:),mapinv(:)  ! nb_list of k
 real(8) shftx,shfty,shftz,shft(3),deltak,kcutnbf,kcutnbc   ! kcutnb for nbrlst and interpolation
 real(8), allocatable :: kpc(:,:),wk(:),wkf(:),kp_bs(:,:),gg(:,:),dk_bs(:),dkc(:)
 real(8), allocatable :: kbzf(:,:) ,kibz(:,:),wibz(:),kbz(:,:)
 integer nibz
 logical dos_bs_only

 contains

  subroutine allocate_kp_bs(n)
  integer n
    allocate (kp_bs(3,n),dk_bs(n))
  end subroutine allocate_kp_bs
!-------------------------------------------
  subroutine make_kp_bs_fcc
use lattice
    implicit none
    integer :: i,j
    real(8) ep,ax,newk,oldk,q(3)

    call allocate_kp_bs(2*nbs1+2*nbs2+2*nbs3)
    ax = lattice_parameter
    ep = 0.001
    nkp_bs = 0
! the following is for FCC
    do i = 1,nbs2    ! Gamma(0,0,0) ---> K(3pi/4,3pi/4,0)--->X(0,0,pi)
       nkp_bs = nkp_bs+1
       kp_bs(1,nkp_bs) = ep + 2*(i-1.)*pi/ax/(nbs2-1.)
       kp_bs(2,nkp_bs) = ep + 2*(i-1.)*pi/ax/(nbs2-1.)
       kp_bs(3,nkp_bs) = 0
       if (i.ge. 2) then
!         dk_bs(nkp_bs)=dk_bs(nkp_bs-1)+length(kp_bs(:,nkp_bs)-kp_bs(:,nkp_bs-1))
          q(:)=kp_bs(:,nkp_bs)-kp_bs(:,nkp_bs-1)
          newk=oldk+length(q)
       else
!         dk_bs(nkp_bs)=length(kp_bs(:,nkp_bs))
          newk=length(kp_bs(:,nkp_bs))
       endif
       dk_bs(nkp_bs)= newk
       oldk = newk
    enddo
    do i = nbs1,1,-1    ! X(pi,0,0) ---> Gamma(0,0,0)
       nkp_bs = nkp_bs+1
       kp_bs(1,nkp_bs) = ep + 2*(i-1.)*pi/ax/(nbs1-0.)
       kp_bs(2,nkp_bs) = 0
       kp_bs(3,nkp_bs) = 0
       dk_bs(nkp_bs)=dk_bs(nkp_bs-1)+  2*pi/ax/nbs1 !length(kp_bs(:,nkp_bs)-kp_bs(:,nkp_bs-1))
    enddo
    do i = 1,nbs1    ! Gamma(0,0,0) --->  L(pi,pi,pi)
       nkp_bs = nkp_bs+1
       kp_bs(1,nkp_bs) = (i-1.)*pi/ax/(nbs1-0.) + ep
       kp_bs(2,nkp_bs) = (i-1.)*pi/ax/(nbs1-0.) + ep
       kp_bs(3,nkp_bs) = (i-1.)*pi/ax/(nbs1-0.) + ep
       dk_bs(nkp_bs)=dk_bs(nkp_bs-1)+length(kp_bs(:,nkp_bs)-kp_bs(:,nkp_bs-1))
    enddo
    do i = 1,nbs2    ! L(pi,pi,pi)  ---> X(pi,0,0)
       nkp_bs = nkp_bs+1
       kp_bs(1,nkp_bs) = pi/ax + (i-1.)*pi/ax/(nbs2-0.)
       kp_bs(2,nkp_bs) = pi/ax - (i-1.)*pi/ax/(nbs2-0.)
       kp_bs(3,nkp_bs) = pi/ax - (i-1.)*pi/ax/(nbs2-0.)
       dk_bs(nkp_bs)=dk_bs(nkp_bs-1)+length(kp_bs(:,nkp_bs)-kp_bs(:,nkp_bs-1))
    enddo
    do i = 1,nbs3    ! X(pi,0,0)    ---> W(pi,pi/2,0)
       nkp_bs = nkp_bs+1
       kp_bs(1,nkp_bs) = 2*pi/ax
       kp_bs(2,nkp_bs) = ep + (i-1.)*pi/ax/(nbs3-0.)
       kp_bs(3,nkp_bs) = 0
       dk_bs(nkp_bs)=dk_bs(nkp_bs-1)+length(kp_bs(:,nkp_bs)-kp_bs(:,nkp_bs-1))
    enddo
    do i = 1,nbs3    ! W(pi,pi/2,0) ---> L(pi,pi,pi)
       nkp_bs = nkp_bs+1
       kp_bs(1,nkp_bs) = 2*pi/ax - (i-1)*pi/ax/(nbs3-0.)
       kp_bs(2,nkp_bs) = pi/ax
       kp_bs(3,nkp_bs) = ep + (i-1.)*pi/ax/(nbs3-0.)
       dk_bs(nkp_bs)=dk_bs(nkp_bs-1)+length(kp_bs(:,nkp_bs)-kp_bs(:,nkp_bs-1))
    enddo
!   do i = 1,nbs3    ! L(pi,pi,pi) ---> K(3pi/4,3pi/4,0)
!      nkp_bs = nkp_bs+1
!      kp_bs(1,nkp_bs) = pi/ax + (i-1)*pi/ax/(nbs3-0.)/2.
!      kp_bs(2,nkp_bs) = pi/ax + (i-1)*pi/ax/(nbs3-0.)/2.
!      kp_bs(3,nkp_bs) = pi/ax - (i-1)*pi/ax/(nbs3-0.)
!   enddo
!   do i = 1,nbs3    ! K(3pi/4,3pi/4,0) --> W(pi,pi/2,0)
!      nkp_bs = nkp_bs+1
!      kp_bs(1,nkp_bs) = 3*pi/ax/2. + (i-1)*pi/ax/(nbs3-0.)/2.
!      kp_bs(2,nkp_bs) = 3*pi/ax/2. - (i-1)*pi/ax/(nbs3-0.)/2.
!      kp_bs(3,nkp_bs) = 0
!   enddo

    do i=1,nkp_bs
       write(ibs,3) kp_bs(:,i)
    enddo
3   format(9(2x,f10.4))

  end subroutine make_kp_bs_fcc
!-------------------------------------------
  subroutine make_kp_bs2
! this subroutine sets up the kpoints for the band structure from the 
! list written in the file kpbs.in (given in direct coordinates)
 use io2
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

 open(ibs,file='KPOINT.BS',status='unknown')
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
 close(ibs)

    deallocate(ki,kf)
3   format(9(2x,f10.4))

    write(ulog,*)' MAKE_KP_BS: Done! with ',nkp_bs,' points'

  end subroutine make_kp_bs2
!-------------------------------------------
  subroutine make_kp_cubic(nx,ny,nz,kpg)
  use lattice
    implicit none
    integer :: i,j,k,nx,ny,nz,nk
    real(8) ax,ay,az,kpg(3,nx*ny*nz)

!    call allocate_kp(nx*ny*nz)
!    allocate( wk(nx*ny*nz) )
!   ax = boxx; ay=boxy; az=boxz
    ax = lattice_parameter
    ay = lattice_parameter
    az = lattice_parameter
    nk = 0
    do i = 1,nx
    do j = 1,ny
    do k = 1,nz
       nk = nk+1
       kpg(1,nk) = -2*pi/ax + i*4*pi/ax/(nx) - pi/2/ax/nx  ! this is shifted
       kpg(2,nk) = -2*pi/ay + j*4*pi/ay/(ny) - pi/2/ay/ny
       kpg(3,nk) = -2*pi/az + k*4*pi/az/(nz) - pi/2/az/nz
!      kpg(1,nk) = -2*pi/ax + i*4*pi/ax/(nx)
!      kpg(2,nk) = -2*pi/ay + j*4*pi/ay/(ny)
!      kpg(3,nk) = -2*pi/az + k*4*pi/az/(nz)
    enddo
    enddo
    enddo
!    wk(:)=1d0/(nx*ny*nz)

  end subroutine make_kp_cubic
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
       q(1) = (i+shftx)*deltak
       q(2) = (j+shfty)*deltak
       q(3) = (k+shftz)*deltak
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
  subroutine make_kp_reg(nx,ny,nz,sx,sy,sz,kpt,wkt)
! generate a regular mesh from 0 to g_i, with eventual shift 
    use geometry
    use params
    implicit none
    integer :: i,j,k,nx,ny,nz,nk
    real(8) q(3),ep(3),sx,sy,sz,kpt(3,nx*ny*nz),wkt(nx*ny*nz)

    open(126,file='KPOINT.MP',status='unknown')
    write(126,*)'# i ,kp(:,i), wk(i)'
    ep = 0d0  ! shift by a small amount so that only one boundary K remains in the FBZ after folding
    shft = (-0.5d0)*(g1+g2+g3)

    nk=nx*ny*nz
    allocate(dkc(nk))

    nk = 0
    do i = 1,nx
    do j = 1,ny
    do k = 1,nz
       nk = nk+1
       q = ((i-1+sx)/nx)*g1 + ((j-1+sy)/ny)*g2 + ((k-1+sz)/nz)*g3 + ep
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

  end subroutine make_kp_reg
!-------------------------------------------
  subroutine make_kp_coarse(nx,ny,nz)
! generate a regular mesh (no shift) from -g_i to g_i, throw away kpoints larger
! that some radius a bit larger than the FBZ's larger radius, the points in the
! FBZ are identified with a mapping integer array called mapbz:
! mapbz(1:nbz) is the index of the corse kpoints which are in the FBZ
! also adds a finer mesh near Gamma for better fitting there.
! input: ni, output: kpc,nkc,nbz,mapbz
    use geometry
    use params
    implicit none
    integer :: i,j,k,n4,nx,ny,nz,nk,indexn !,nkbz
    real(8) q(3),ep(3),kmax,g6
    real(8), allocatable :: kp(:,:),kz(:,:),kl(:),junk(:,:)
    integer, allocatable :: mp(:)
!   type(vector) vv
!    real(8) kp(3,nk)

    nk = (2*nx+1)*(2*ny+1)*(2*nz+1)
    allocate(kp(3,nk),kl(nk))
    open(126,file='KPOINT.COARSE',status='unknown')
    ep = 0d0  !5*tolerance  ! 0.00005
!   write(ulog,*)' KP_COARSE: i, kpc(:,i) ##############'

! shift by a small amount so that only one boundary K remains in the FBZ after folding
    nk = 0
    do i = -nx,nx
    do j = -ny,ny
    do k = -nz,nz
       nk = nk+1
       q = ((i+0d0+shftx)/nx)*g1 + ((j+0d0+shfty)/ny)*g2 + ((k+0d0+shftz)/nz)*g3 + ep
       kp(:,nk) = q(:)
       kl(nk) = length(q)
!      write(ulog,4)nk,q
!      write(126,3)q
       if( nk .ne. indexn(i,j,k,nx,ny,nz) ) then
!       if( nk .ne. (k+nz+1)+(j+ny)*(2*nz+1)+(i+nx)*(2*ny+1)*(2*nz+1) ) then
          write(ulog,*)' ERROR in the coarse labeling or indexn labeling'
          stop
       endif
    enddo
    enddo
    enddo
    write(ulog,*)' Number of raw coarse kpoints generated is=',nk
    allocate(mp(nk))

! sort according to their lengths
    call sort(nk,kl,mp,nk)

! find the radius of the sphere containing the FBZ
    kmax = 0
    do i=1,15 ! nshl do not need the third shell here
       g6=length(gg(:,i))
       if (g6 .gt. kmax) then
          kmax = g6
       endif
    enddo
    g6 = kmax*0.7
    write(ulog,*)' radius of the coarse mesh is chosen to be=',g6

! find how many kp are within this sphere
    iloop: do i=1,nk
       if (kl(mp(i)) .gt. g6) then
          nkc = i
          exit iloop
       endif
    enddo iloop
    write(ulog,*)' size of the coarse mesh is=',nkc

! store the points within the spere in junk
    allocate(junk(3,nkc))
    do i=1,nkc
       junk(:,i) = kp(:,mp(i))
    enddo

! now fold in the FBZ and find out how many kpoints are in there and their mapping
    allocate(mapbz(nkc),kz(3,nkc))
    call fold_in_fbz(nkc,junk,1,gg,nbz,kz,mapbz)

! mapbz is the mapping of a point in the FBZ to the index of that point in kpc
    write(ulog,*)'Number of coarse kpoints generated in the FBZ =',nbz
    do i=1,nbz
       write(fbz,3)junk(:,mapbz(i))
    enddo
    deallocate(kp,kz,mp)

! add a small fine mesh around the gamma point
    n4 = 4
    allocate(kz(3,n4*n4*n4),kp(3,n4*n4*n4))
    call make_kp_cubic(n4,n4,n4,kz)
    kz = kz/7  ! make it 7 times finer
    kz = kz+0.06*kmax/7   ! shift it
    j=0
    do i=1,n4**3
       g6 = length(kz(:,i))
       if (g6.lt.kmax/9) then
          j=j+1
          kp(:,j) = kz(:,i)
       endif
    enddo
    write(ulog,*)j,' mesh points added near the gamma point'

! now allocate and store the coarse mesh in this array, and the KPOINT.COARSE file
    allocate(kpc(3,nkc+j))
    do i=1,nkc
       kpc(:,i) = junk(:,i)     ! junk(:,mapbz(i)) lists the folded vectors in the FBZ
       write(126,3)kpc(:,i)
    enddo

! add the extra sphere around gamma
    do i=nkc+1,nkc+j
       kpc(:,i) = kp(:,i-nkc)
       write(126,3)kpc(:,i)
    enddo
    nkcc = nkc     ! this is the # of coarse mesh
    nkc = nkc+j    ! this is the total # of points including the small cubic mesh

    close(126)
    deallocate(kp,kl,junk,kz)


!   allocate(mapib(nbz))
!   call map_ibz(nbz,kz,mapib,nkc_ib)
!   deallocate(kz)

3 format(9(2x,f9.4))
4 format(i6,9(2x,f9.4))

  end subroutine make_kp_coarse
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
 use io2
 implicit none
 integer, intent(in) :: nk,nshl
 integer, intent(out) :: nbz,mapz(nk)
 real(8), intent(in) :: kp(3,nk),gg(3,nshl)
 real(8), intent(out) :: kz(3,nk)
 integer i,ns,j
 real(8) qaux(3),dmin,dist
 logical inside

! write(ulog,*)'FOLD IN FBZ: ik,shell,i,q ====================='
 nbz=0
 do i=1,nk
 foldloop: do ns=1,nshl
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
 subroutine get_weights(nk,kp)
! takes as input kp(3,nk) in cartesian and finds based on crystal symmetries, the
! weights associated with each kpoint, and then normalizes them
! nibz is the final number of kpoints stored in kibz in the irreducible BZ
! the other important output is the mapping of the full kmesh onto the ones
! folded in the irreducible zone : mapibz(j=1:nk) is the index of the k in the IBZ
! corresponding to the argument j
! for i=1,nibz mapinv(i) gives the index k of the kpoints generated from
! n1,n2,n3 loops
 use lattice
 use constants
 use io2
 use geometry
 use params
 implicit none
 integer nk,nkibz,i,j,l,narms,kvecop(48),aux
 integer, allocatable :: mcor(:)
 real(8) zro,q(3),kvecstar(3,48),sw,skc(3)
 real(8) kp(3,nk) ,st(3),qaux(3)
 real(8) , allocatable :: k2(:,:),lg(:),w2(:)
 logical exists,foundit

 open(ibz,file='KPOINT.IBZ',status='unknown')
 allocate(k2(3,nk),w2(nk),mapibz(nk),mapinv(nk))
 zro=0d0  ! shift
 write(ulog,*)'MAKE_KP_IBZ: generating kpoints in the irreducible FBZ '
 write(*,*)'MAKE_KP_IBZ: nk,wk(nk),k2(:,nk)'
 write(ulog,*)'primitivelattice=',primitivelattice
 nibz=0 ; mapibz=0; w2=1

! initialize mapvinv with identity so that later elements can be switched
 do i=1,nk  
    mapinv(i)=i
 enddo

! main loop to identify points in the FBZ
 kploop: do i=1,nk

    call getstar(kp(:,i),primitivelattice,narms,kvecstar,kvecop)

! if one of the stars differs from an existing kibz, then skip k
    do l=1,narms
       do j=1,nibz
          qaux= (k2(:,j) - kvecstar(:,l))
          st(1)=(qaux(:) .dot. r1) /2/pi
          st(2)=(qaux(:) .dot. r2) /2/pi
          st(3)=(qaux(:) .dot. r3) /2/pi
!         write(ulog,7)'i,lstar,jibz,(kstar-kibz)_red=',i,l,j,st
          qaux=0
          if ( st-nint(st) .myeq. qaux ) then ! this is an already irreducible point
             w2(j)=w2(j)+1
             cycle kploop
          else   
          endif
       enddo 
    enddo
! by now kstar was none of the generated kibz, so it is a new one
    nibz=nibz+1
    k2(:,nibz)=kp(:,i)

    if (nibz.eq.1) w2(nibz) = 1
    mapibz(i) = nibz
       
! here need to switch two elements of mapinv; so we add this line
    aux=mapinv(nibz)
    mapinv(nibz)=i
    mapinv(i)=aux

    write(ulog,4)'new vector*:',nibz,w2(nibz),k2(:,nibz),length(k2(:,nibz)),q
    write(*,4)'new vector*:',nibz,w2(nibz),k2(:,nibz),length(k2(:,nibz))

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
!
! this sorting and writing is donw in the main routine
! if ( allocated(lg)) deallocate(lg)
! if ( allocated(mcor)) deallocate(mcor)
! allocate(lg(nk),mcor(nk))
! do i=1,nk
!    lg(i)=length(kp(:,i))
! enddo
! call sort(nk,lg,mcor,nk)
! write(fbz,*)'#i,l,mapibz(l),kp(l),length(kp(l)),kibz(mapibz(l)),wibz(mapibz(l)) l=mcor(i)'
! do i=1,nk
!     j=mcor(i)
!!    write(fbz,3)i,kp(:,mcor(i)),length(kp(:,mcor(i))),kibz(:,mapibz(mcor(i))),wibz(mapibz(mcor(i)))
!    write(fbz,2)i,j,mapibz(j),kp(:,j),length(kp(:,j)),kibz(:,mapibz(j)),wibz(mapibz(j))
! enddo
! deallocate(mcor,lg)
! close(fbz)

 write(ibz,*)'#i,kibz(:,i),wibz(i),kibz_reduced,length(kibz(i))'
 open(345,file='mapinv.dat')
 write(345,*)'i,mapinv(i)'
 do i=1,nk !,nibz
  if(i.le.nibz)  then
    q(1)=(kibz(:,i)  .dot. r1) /2/pi
    q(2)=(kibz(:,i)  .dot. r2) /2/pi
    q(3)=(kibz(:,i)  .dot. r3) /2/pi
    write(ibz,3)i,kibz(:,i),wibz(i),q,length(kibz(:,i))
  endif
    write(345,*)i,mapinv(i)
 enddo
 close(ibz)
 close(345)

2 format(3i7,2x,3(1x,f9.4),2x,f9.4,2x,3(1x,f9.5),3x,f14.8)
3 format(i7,2x,3(1x,f9.4),2x,f14.8,2x,3(1x,f9.5),3x,f9.5)
4 format(a,i7,2x,f9.5,2x,99(1x,f9.5))
5 format(a,2x,99(1x,f9.5))
6 format(a,i7,2x,i5,2x,f9.5,2x,99(1x,f9.5))
7 format(a,i7,2x,i5,2x,i5,2x,f9.5,2x,99(1x,f9.5))

 end subroutine get_weights
!===========================================================
!----------------------------------------------
 subroutine get_weights2(nk,kp)
! takes as input kp(3,nk) in cartesian and finds based on crystal symmetries, the
! weights associated with each kpoint, and then normalizes them
! nibz is the final number of kpoints stored in kibz in the irreducible BZ
! the other important output is the mapping of the full kmesh onto the ones
! folded in the irreducible zone : mapibz(j=1:nk) is the index of the k in the IBZ
! corresponding to the argument j
! for i=1,nibz mapinv(i) gives the index k of the kpoints generated from
! n1,n2,n3 loops
 use lattice
 use constants
 use io2
 use geometry
 use params
 implicit none
 integer nk,nkibz,i,j,l,narms,kvecop(48),aux
 integer, allocatable :: mcor(:)
 real(8) zro,q(3),kvecstar(3,48),sw,skc(3)
 real(8) kp(3,nk) ,st(3),qaux(3)
 real(8) , allocatable :: k2(:,:),lg(:),w2(:)
 logical exists,foundit

 open(ibz,file='KPOINT.IBZ',status='unknown')
 allocate(k2(3,nk),w2(nk),mapibz(nk),mapinv(nk))
 zro=0d0  ! shift
 write(ulog,*)'MAKE_KP_IBZ: generating kpoints in the irreducible FBZ '
 write(*,*)'MAKE_KP_IBZ: nk,wk(nk),k2(:,nk)'
 write(ulog,*)'primitivelattice=',primitivelattice
 nibz=0 ; mapibz=0; w2=1

! initialize mapvinv with identity so that later elements can be switched
 do i=1,nk  
    mapinv(i)=i
 enddo

! main loop to identify points in the FBZ
 kploop: do i=1,nk
!    q = kp(:,i)
! kp is in cartesian coordinates, we need to convert it to reduced units:
    q(1)=(kp(:,i)  .dot. r1) /2/pi
    q(2)=(kp(:,i)  .dot. r2) /2/pi
    q(3)=(kp(:,i)  .dot. r3) /2/pi
! below the cartesian components of kp are needed
    call getstar(kp(:,i),primitivelattice,narms,kvecstar,kvecop)
  ! call getstar(q      ,primitivelattice,narms,kvecstar,kvecop)

! now shift them back by G in the same region[-G/2:G/2[ to identify them 
    do l=1,narms
       write(ulog,4)'i,kp,q=',i,kp(:,i),q,kvecstar(:,l)
       st(1)=(kvecstar(:,l) .dot. r1) /2/pi
       st(2)=(kvecstar(:,l) .dot. r2) /2/pi
       st(3)=(kvecstar(:,l) .dot. r3) /2/pi
! this is needed if the original mesh is centered at 0
 !     qaux(1)=st(1)-nint(st(1))  
 !     qaux(2)=st(2)-nint(st(1))
 !     qaux(3)=st(3)-nint(st(3))
       qaux=st
       kvecstar(:,l)=qaux(1)*g1+qaux(2)*g2+qaux(3)*g3
       write(ulog,4)'moved kvecstar in the BZ=',i,kvecstar(:,l),qaux
    enddo
       
! see if already exists among the previously assigned ones
    exists=.False.
    if(verbose) write(ulog,6)'list of kvecstar(l),l=1,narms for kp_red=',i,narms,q
    lloop: do l=1,narms
    if(verbose)   write(ulog,4)'stars are:',l,kvecstar(:,l)
! set weight for the first kpoint where nibz=0
        jloop: do j=1,nibz
          if (k2(:,j) .myeq. kvecstar(:,l)) then
! first bring the star in the FBZ, then compare to the existing points
             exists=.True.
             w2(j)=w2(j)+1d0
             mapibz(i)=j
             write(ulog,4)' this kpoint turned out to exist ',j,k2(:,j)
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
       
! here need to switch two elements of mapinv; so we add this line
       aux=mapinv(nibz)
       mapinv(nibz)=i
       mapinv(i)=aux

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
!
! this sorting and writing is donw in the main routine
! if ( allocated(lg)) deallocate(lg)
! if ( allocated(mcor)) deallocate(mcor)
! allocate(lg(nk),mcor(nk))
! do i=1,nk
!    lg(i)=length(kp(:,i))
! enddo
! call sort(nk,lg,mcor,nk)
! write(fbz,*)'#i,l,mapibz(l),kp(l),length(kp(l)),kibz(mapibz(l)),wibz(mapibz(l)) l=mcor(i)'
! do i=1,nk
!     j=mcor(i)
!!    write(fbz,3)i,kp(:,mcor(i)),length(kp(:,mcor(i))),kibz(:,mapibz(mcor(i))),wibz(mapibz(mcor(i)))
!    write(fbz,2)i,j,mapibz(j),kp(:,j),length(kp(:,j)),kibz(:,mapibz(j)),wibz(mapibz(j))
! enddo
! deallocate(mcor,lg)
! close(fbz)

 write(ibz,*)'#i,kibz(:,i),wibz(i),kibz_reduced,length(kibz(i))'
 open(345,file='mapinv.dat')
 write(345,*)'i,mapinv(i)'
 do i=1,nk !,nibz
  if(i.le.nibz)  then
    q(1)=(kibz(:,i)  .dot. r1) /2/pi
    q(2)=(kibz(:,i)  .dot. r2) /2/pi
    q(3)=(kibz(:,i)  .dot. r3) /2/pi
    write(ibz,3)i,kibz(:,i),wibz(i),q,length(kibz(:,i))
  endif
    write(345,*)i,mapinv(i)
 enddo
 close(ibz)
 close(345)

2 format(3i7,2x,3(1x,f9.4),2x,f9.4,2x,3(1x,f9.5),3x,f14.8)
3 format(i7,2x,3(1x,f9.4),2x,f14.8,2x,3(1x,f9.5),3x,f9.5)
4 format(a,i7,2x,f9.5,2x,99(1x,f9.5))
5 format(a,2x,99(1x,f9.5))
6 format(a,i7,2x,i5,2x,f9.5,2x,99(1x,f9.5))

 end subroutine get_weights2
!===========================================================
subroutine sort3(n1,n2,n3,mp)
implicit none
integer n1,n2,n3,mp(3)

if( n2.ge.n1) then
  if (n3.ge.n2) then
     mp(1)=n1
     mp(2)=n2
     mp(3)=n3
  elseif( n3.ge.n1) then
     mp(1)=n1
     mp(2)=n3
     mp(3)=n2
  else
     mp(1)=n3
     mp(2)=n1
     mp(3)=n2
  endif
elseif( n3.ge.n1) then
     mp(1)=n2
     mp(2)=n1
     mp(3)=n3
 elseif(n3.ge.n2) then
     mp(1)=n2
     mp(2)=n3
     mp(3)=n1
 else
     mp(1)=n3
     mp(2)=n2
     mp(3)=n1
 endif
end subroutine sort3
!----------------------------------------------


 end module kpoints
!===========================================================
