!===========================================================

 module kpoints
!! this includes 4 kpoint meshes: one for the band structure: kp_bs(3,nkp_bs)
! one Monkhorst-Pack mesh around gamma for interpolation : kpc(3,nkc)
! a mapping of this mesh into the FBZ: kpc(3,mapbz(1:nbz))
! one grid in the irreducible FBZ for kappa/dos sums: kibz(3,nibz)
 use lattice   ! it is also used for all included subroutines
 use ios
 use constants
 implicit none
 integer nbs1,nbs2,nbs3,nkp_bs   ! for the band structure
 integer nc(3),nkc,nbz,nkcc  ! coarse mesh,its FBZ folded mesh,IBZ_coarse mesh
 integer nib1,nib2,nib3,nbz_fine  ! for the fine mesh in the Irreducible BZ
 integer nshl,nbkmax,npos,nkmesh,nrmesh,ngibz         ! max no of neighbors=20
 integer, allocatable :: nb(:,:),nbk(:),mapbz(:) ,mapibz(:),mappos(:),mapinv(:)  ! nb_list of k
 real(r15) shift(3),normal(3),deltak,kcutnbf,kcutnbc,k_cart(3),k_conv(3),k_prim(3)  ! kcutnb for nbrlst and interpolation
 real(r15), allocatable :: kpc(:,:),wk(:),kp_bs(:,:),dk_bs(:),dkc(:),gibz(:,:),wgibz(:)
 real(r15), allocatable :: kbzf(:,:) ,kibz(:,:),wibz(:),kbz(:,:),kpfbz(:,:),kmesh(:,:),rmesh(:,:)
 integer, allocatable :: fbz2ibz(:),kop_ibz2fbz(:),fbzKstar(:)
 integer, allocatable :: ibzarms(:),ibz2fbz(:)
 character(LEN=3), allocatable :: kname_bs(:)
 integer nibz
 logical dos_bs_only

 contains

!-------------------------------------------
  subroutine allocate_kp_bs(n)
  integer n
    allocate (kp_bs(3,n),dk_bs(n))
  end subroutine allocate_kp_bs
!-------------------------------------------
  subroutine deallocate_kp_bs
    deallocate (kp_bs,dk_bs)
  end subroutine deallocate_kp_bs
!-------------------------------------------
  subroutine make_kp_bs
!! this subroutine sets up the kpoints for the band structure from the
!! list written in the file kpbs.params (given in direct coordinates of the CONVENTIONAL cell)
  use ios
  use geometry
  use lattice , only : g01,g02,g03,volume_g0
  use constants, only : pi
    implicit none
    integer :: i,j,nk,uio,nkdir,ubs,units ,ndir
    real(r15) lk,q(3),k_conv(3),k_prim(3),k_cart(3),rand(3)
    real(r15), allocatable :: ki(:,:),kf(:,:) ,kext_bs(:,:)
    character(LEN=6), allocatable :: dk(:)

! we use this small random shift to move away from high symmetry points and remove degeneracies
! but need to worry about band connections -> band_sort (in case finite_diff is used)
! if dk is small enough, but larger than dq=1d-4 below, then this should not affect bandsorting.
  call random_number(rand)
  rand=0 !(2*rand-1)*1d-4 *volume_g0**0.33333333 

    write(*,*)'entering MAKE_KP_BS'
    uio = 67
    ubs = 68
    open(uio,file='kpbs.params',status='old')
    open(ubs,file='KPOINT.BS')

    read(uio,*) units  ! if 0 conventional, else primitive
    write(*,*)'reading ',units
    read(uio,*) nkdir  ! number of kpoints along each direction
    write(*,*)'reading ',nkdir
    read(uio,*) ndir  ! number of directions for the band structure
    write(ubs,*)'# reading units,nkdir,ndir=',units,nkdir,ndir
    write(ubs,*)'# k_name , k in primitive, k in gconv units, k in cartesian '
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
          write(ubs,9)'# '//kname_bs(i),i, k_prim,k_conv,k_cart

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
          write(ubs,9)'# '//kname_bs(i),i, kext_bs(:,i),k_conv,k_cart

       endif
    enddo

    i=ndir+1
    if(units.eq.0) then ! reading in units of conventional
          k_cart=kext_bs(1,i)*g1conv+kext_bs(2,i)*g2conv+kext_bs(3,i)*g3conv
          k_prim=matmul(transpose(prim_to_conv),kext_bs(:,i))
          write(ubs,9)'# '//kname_bs(i),i, k_prim,kext_bs(:,i),k_cart
    else
          k_cart=kext_bs(1,i)*g01+kext_bs(2,i)*g02+kext_bs(3,i)*g03
          k_conv=matmul(transpose(conv_to_cart),k_cart)/2/pi
          write(ubs,9)'# '//kname_bs(i),i, kext_bs(:,i),k_conv,k_cart
    endif

    close(uio)
!   close(ubs)

  9 format(a,i4,9(1x,f9.5))

  write(ulog,*)' Kpoints for band structure generated from kpbs.params'
  write(*,*)' Kpoints for band structure generated from kpbs.params'


 open(88,file='KTICS.BS',status='unknown')
    nk = 0
!    dk_bs(1) = 0
! now for each direction set up the coordinates of the kpoints
    kext_bs(1,1)=0    ! here using kext_bs(1,:) as a dummy variable
    do i = 1,ndir !-1
       q(:) = kf(:,i)-ki(:,i)  ! length of each section
       lk = length(q)/length(g01)  ! length of each section
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
          kp_bs(:,nk) = ki(:,i) + (j-1)*q(:)/(nkdir-1+1d-8) + rand
          write(ubs,7)nk, kp_bs(:,nk)
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
    close(88)
    close(ubs)

    deallocate(ki,kf,kname_bs)
3   format(9(2x,f10.4))
7   format(i4,99(2x,f10.4))
5   format(a,i4,9(2x,f10.4))
!4   format(a,(ndir+1)(a,1x,f6.4,a))

    write(ulog,*)' MAKE_KP_BS: Done! with ',nkp_bs,' points'
    write(*,*)' MAKE_KP_BS: Done! with ',nkp_bs,' points'

  end subroutine make_kp_bs
!-------------------------------------------
  subroutine make_grid(ncs,gg1,gg2,gg3,shft,kpt,wkt)
!! generate a regular mesh integer multiples of g_i, with eventual shift
    use geometry
    use params
    use lattice, only : cart2red_g
    implicit none
    integer, intent(in) :: ncs(3)
    real(r15), intent(in) :: shft(3) 
    real(r15), intent(out):: kpt(3,ncs(1)*ncs(2)*ncs(3)),wkt(ncs(1)*ncs(2)*ncs(3))
    type(vector), intent(in) :: gg1,gg2,gg3
    integer :: i,j,k,nk
    open(126,file='grid.MP',status='unknown')
    write(126,*)'# i ,kp(:,i), wk(i), kp_red'

    nkc=ncs(1)*ncs(2)*ncs(3)

    nk = 0
    do i = 1,ncs(1)
    do j = 1,ncs(2)
    do k = 1,ncs(3)
       nk = nk+1
       kpt(:,nk) = (i-1+shft(1))*gg1 + (j-1+shft(2))*gg2 + (k-1+shft(3))*gg3 
       wkt(nk)=1d0/nkc
       write(126,2)nk,kpt(:,nk), wkt(nk), cart2red_g(kpt(:,nk))
    enddo
    enddo
    enddo

    write(ulog,*)'MAKE_GRID: Number of regular kpoints generated is=',nkc

 2  format(i7,2x,3(1x,f12.5),5x,f9.5,2x,3(1x,f10.5))
    close(126)

  end subroutine make_grid
!-------------------------------------------
  subroutine make_kp_reg(ncs,gg1,gg2,gg3,shft,kpt,wkt)
!! generate a regular mesh from 0 to g_i, with eventual shift
    use geometry
    use params
    implicit none
    integer, intent(in) :: ncs(3)
    real(r15), intent(in) :: shft(3) !,g1(3),g2(3),g3(3)
    real(r15), intent(out):: kpt(3,ncs(1)*ncs(2)*ncs(3)),wkt(ncs(1)*ncs(2)*ncs(3))
    type(vector), intent(in) :: gg1,gg2,gg3
    type(vector) rr1,rr2,rr3
    real(r15) a1,a2,a3

    integer :: i,j,k,nk
    real(r15) q(3),ep(3)
    open(126,file='KPOINT.MP',status='unknown')
    write(126,*)'# i ,kp(:,i), wk(i), kp_red'
    ep = 0d0  ! shift by a small amount so that only one boundary K remains in the FBZ after folding
!   shft = (-0.5d0)*(g1+g2+g3)

    nkc=ncs(1)*ncs(2)*ncs(3)
    allocate(dkc(nkc))

    nk = 0
    do i = 1,ncs(1)
    do j = 1,ncs(2)
    do k = 1,ncs(3)
       nk = nk+1
       q = ((i-1+shft(1))/ncs(1))*gg1 + ((j-1+shft(2))/ncs(2))*gg2 + ((k-1+shft(3))/ncs(3))*gg3 + ep
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

    call make_reciprocal_lattice_v(gg1,gg2,gg3,rr1,rr2,rr3)

    do i=1,nkc
       q=kpt(:,i)
       a1=rr1.dot.q
       a2=rr2.dot.q
       a3=rr3.dot.q
       write(126,2)i,q,wkt(i),a1,a2,a3
    enddo

 2  format(i7,2x,3(1x,f12.5),5x,f9.5,2x,3(1x,f10.5))
 3  format(3(i3),2x,i6,2x,3(1x,f12.5),5x,f9.5)
    close(126)

  end subroutine make_kp_reg
!-------------------------------------------
 subroutine get_k_info(q,N,nk,i,j,k,inside)
! for a vector q(3) in the primitive cell of the reciprocal space, defined on
! a mesh N(1),N(2),N(3), this subroutine finds the three indices of q
! and its number nk based on the triple loop
! nk=0 do i=0,N(1)-1 ; do j=0,N(2)-1; do k=0,N(3)-1; nk=nk+1 ;q=i*G1/N1+j*G2/N2+k*G3/N3
! it works even if q is outside because get_components finds i,j,k in the FBZ
 use geometry
 implicit none
 integer, intent(in):: N(3)
 real(r15), intent(in):: q(3)
 integer, intent (out) :: nk,i,j,k,inside
 integer indexg

! this was for when the kpoints were between 0 and N-1
    call get_components_g(q,N,i,j,k,inside)
    nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
 !  nk= indexg(i,j,k,-N(1)/2,N(1)/2-1,-N(2)/2,N(2)/2-1,-N(3)/2,N(3)/2-1)
    return

! this is for when the k vectors are between -N/2 and N/2-1
  if (mod(N(1),2).eq.0 .and. mod(N(2),2).eq.0 .and. mod(N(3),2).eq.0 ) then
! use this if shifted by -0.5(g1+g2+g3)
    call get_components_g_centered (q,N,i,j,k,inside)
!   write(*,*)'for even ni, ijk=',i,j,k
  else   ! do the regulat thing
    call get_components_g(q,N,i,j,k,inside)
  endif
  nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
! write(*,*)'the corresponding nk=',nk

 end subroutine get_k_info
!----------------------------------------------
 subroutine get_k_info_cent(q,N,nk,i,j,k,inside)
! for a vector q(3) in the primitive cell of the CENTERED-reciprocal space, defined on
! a mesh N(1),N(2),N(3), this subroutine finds the three indices of q
! and its number nk based on the triple loop
! nk=0 do i=0,N(1)-1 ; do j=0,N(2)-1; do k=0,N(3)-1; nk=nk+1 ;q=i*G1/N1+j*G2/N2+k*G3/N3

 use geometry
 implicit none
 integer, intent(in):: N(3)
 real(r15), intent(in):: q(3)
 integer, intent (out) :: nk,i,j,k
 integer, intent (inout) :: inside
 integer indexg

! this was for when the kpoints were between 0 and N-1
    call get_components_g_centered (q,N,i,j,k,inside)
    nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
 !  nk= indexg(i,j,k,-N(1)/2,N(1)/2-1,-N(2)/2,N(2)/2-1,-N(3)/2,N(3)/2-1)
    return

! this is for when the k vectors are between -N/2 and N/2-1
  if (mod(N(1),2).eq.0 .and. mod(N(2),2).eq.0 .and. mod(N(3),2).eq.0 ) then
! use this if shifted by -0.5(g1+g2+g3)
    call get_components_g_centered (q,N,i,j,k,inside)
!   write(*,*)'for even ni, ijk=',i,j,k
  else   ! do the regulat thing
    call get_components_g(q,N,i,j,k,inside)
  endif
  nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
! write(*,*)'the corresponding nk=',nk

 end subroutine get_k_info_cent
!----------------------------------------------
 subroutine get_k_info2(q,N,nk,i,j,k,inside)
! scans the kpoints to identify the index of the input q
! also works is outside of the FBZ cell
 use geometry
 use params, only : tolerance
 use constants, only : pi
 use ios
 use lattice, only : prim_to_cart,cart_to_prim
 implicit none
 integer, intent(in):: N(3)
 real(r15), intent(in):: q(3)
 integer, intent (out) :: nk,i,j,k,inside
 integer l,nq
 real(r15) q2(3),bring_to_prim_cell_caa

! first bring q in between 0 and g0i 
 call unitcell(transpose(prim_to_cart)/(2*pi),transpose(cart_to_prim)*2*pi,q,q2)
 write(*,4)'q =',q ,cart2red_g(q)
 write(*,4)'q2=',q2,cart2red_g(q2)
 inside=0
 nk=0; 
 loop: do l=1,nkc
    if( length(q2-kpc(:,l)) .lt. tolerance ) then
       nk=l
       inside=1
       exit loop
    endif
 enddo loop
! write(*,*)'nkt=',nkt
 if (inside.eq.1) then
 ! assumed i=1,n1; j=1,n2; k=1,n3 so that nk = k + (j-1)*n3 + (i-1)*n2*n3
    k=1+mod(nk-1,N(3))
    j=1+mod((nk-k)/N(3),N(2))
    i=1+(nk-k-(j-1)*N(3))/(N(2)*N(3))
    nq = k + (j-1)*N(3) + (i-1)*N(2)*N(3)
 !  if (nk.eq.0 .or. nq.ne.nk) then
       write(ulog,3)'GET_K_INFO2: qpoint not found! nk,nq,i,j,k,q ',nk,nq,i,j,k,cart2red_g(q),q
 !     stop
 !  endif
 else
   i=0;j=0;k=0;nk=0
 endif

3 format(a,5i4,9(1x,f10.4))
4 format(a,9(1x,f10.4))

 end subroutine get_k_info2
!----------------------------------------------
 subroutine get_k_info3(q,nkt,exists) 
! scans the kpoints to identify the index of the input q
 use geometry
 use params, only : tolerance
 use ios
 implicit none
 real(r15), intent(in):: q(3)
 integer, intent(out):: nkt
 logical, intent(out):: exists
 integer i

3 format(a,i6,2x,3(1x,g10.3),2x,f10.4,L3)

 nkt=0; exists=.false.
 loop: do i=1,nkc
    if( length(q-kpc(:,i)) .lt. tolerance ) then
       nkt=i
       exists=.true.
       exit loop
    endif
 enddo loop
! write(*,*)'nkt=',nkt
 if (nkt.eq.0) then
   write(ulog,3)'GET_K_INFO3: qpoint not found! nq,q,l,exist?',nkt,q,length(q),exists
!  stop
 endif
 end subroutine get_k_info3
!----------------------------------------------

 subroutine get_k_info_general(q,kp,nk,nq,exists) 
! scans the kpoints to identify the index of the input q
 use geometry
 use params, only : tolerance
 use ios
 implicit none
 integer, intent(in):: nk
 real(r15), intent(in):: q(3),kp(3,nk)
 integer, intent(out):: nq
 logical, intent(out):: exists
 integer i

3 format(a,i6,2x,3(1x,g10.3),2x,f10.4,L3)

 nq=0; exists=.false.
 loop: do i=1,nk
    if( length(q-kp(:,i)) .lt. tolerance ) then
       nq=i
       exists=.true.
       exit loop
    endif
 enddo loop
 if (nq.eq.0) then
   write(ulog,3)'GET_K_INFO3: qpoint not found! nq,q,l,exist?',nq,q,length(q),exists
!  stop
 endif
 end subroutine get_k_info_general
!----------------------------------------------

end module kpoints

