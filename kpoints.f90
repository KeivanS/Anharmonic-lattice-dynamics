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
 integer nc(3),nkc,nbz,nibz_coarse,nkcc  ! coarse mesh,its FBZ folded mesh,IBZ_coarse mesh
 integer nib1,nib2,nib3,nbz_fine  ! for the fine mesh in the Irreducible BZ
 integer nshl,nbkmax,npos,nkmesh,nrmesh,ngibz         ! max no of neighbors=20
 integer, allocatable :: nb(:,:),nbk(:),mapbz(:) ,mapibz(:),mappos(:),mapinv(:)  ! nb_list of k
 real(8) shift(3),deltak,kcutnbf,kcutnbc   ! kcutnb for nbrlst and interpolation
 real(8), allocatable :: kpc(:,:),wk(:),kp_bs(:,:),dk_bs(:),dkc(:),gibz(:,:),wgibz(:)
 real(8), allocatable :: kbzf(:,:) ,kibz(:,:),wibz(:),kbz(:,:),kpfbz(:,:),kmesh(:,:),rmesh(:,:)
 integer, allocatable :: fbz2ibz(:),kop_ibz2fbz(:),fbzKstar(:)
 integer, allocatable :: ibzarms(:),ibz2fbz(:)
 integer nibz
 logical dos_bs_only

 contains

  subroutine allocate_kp_bs(n)
  integer n
    allocate (kp_bs(3,n),dk_bs(n))
  end subroutine allocate_kp_bs
!-------------------------------------------
  subroutine make_kp_bs
! this subroutine sets up the kpoints for the band structure from the
! list written in the file kpbs.in (given in direct coordinates)
  use ios
  use geometry
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
    call allocate_kp_bs(nkp_bs)
    allocate(ki(3,ndir),kf(3,ndir))

    do i=1,ndir
       read(uio,*) q(:),q2(:)   ! in direct coordinates of g0i
       call dir2cart_g(q ,ki(:,i))
       call dir2cart_g(q2,kf(:,i))
    enddo

    close(uio)

 open(uibs,file='KPOINT.BS',status='unknown')
    nk = 0
    do i = 1,ndir  ! for each direction set up the coordinates of the kpoints
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
          kp_bs(:,nk) = ki(:,i) + (j-1)*q(:)/(nkdir-1+1d-8)
          write(uibs,3) kp_bs(:,nk)
       enddo
    enddo
 close(uibs)

    deallocate(ki,kf)
3   format(9(2x,f10.4))

    write(ulog,*)' MAKE_KP_BS: Done! with ',nkp_bs,' points'

  end subroutine make_kp_bs
!-------------------------------------------
  subroutine make_kp_reg(ncs,gg1,gg2,gg3,shft,kpt,wkt)
!! generate a regular mesh from 0 to g_i, with eventual shift
    use geometry
    use params
    implicit none
    integer, intent(in) :: ncs(3)
    real(8), intent(in) :: shft(3) !,g1(3),g2(3),g3(3)
    real(8), intent(out):: kpt(3,ncs(1)*ncs(2)*ncs(3)),wkt(ncs(1)*ncs(2)*ncs(3))
    type(vector), intent(in) :: gg1,gg2,gg3
    type(vector) rr1,rr2,rr3
    real(8) a1,a2,a3

    integer :: i,j,k,nk
    real(8) q(3),ep(3)
    open(126,file='KPOINT.MP',status='unknown')
    write(126,*)'# i ,kp(:,i), wk(i), kp_red'
    ep = 0d0  ! shift by a small amount so that only one boundary K remains in the FBZ after folding
!   shft = (-0.5d0)*(g1+g2+g3)

    nk=ncs(1)*ncs(2)*ncs(3)
    allocate(dkc(nk))

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

    do i=1,nk
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
!----------------------------------------------
 subroutine fold_in_fbz_new(nk,kp,gg,nbz,kz,mapz)
!! folds the k-mesh defined by kp into the FBZ and stores the result in kz
 use lattice
 use ios
 implicit none
 integer, intent(in) :: nk
 integer, intent(out) :: nbz,mapz(nk)
 real(8), intent(in) :: kp(3,nk),gg(3,26)
 real(8), intent(out) :: kz(3,nk)
 integer i,ns,j,inside
 real(8) qaux(3),dmin,dist


! write(ulog,*)'FOLD IN FBZ: ik,shell,i,q ====================='
 nbz=0
 do i=1,nk
 foldloop: do ns=1,26
     qaux = kp(:,i)+gg(:,ns)
     call check_inside_ws(qaux,gg,inside)
     if(inside.eq.1) then
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

 end subroutine fold_in_fbz_new
!----------------------------------------------
 subroutine fold_in_bz_new(q,gshel,qin)
!! folds the vector q into the FBZ and stores the result in qin
! use lattice
 use ios , only : ulog
 use geometry
 implicit none
 real(8), intent(in):: gshel(3,26)
 real(8), intent(in) :: q(3)
 real(8), intent(out) :: qin(3)
 integer i,j,k,inside
 real(8) gt(3),gmin,leng
 real(8) w(3)

 inside=0
 call check_inside_ws(q,gshel,inside)
 if(inside.eq.1)      return

 qin=q
 leng=length(q)
 do j=1,26 !while (leng.gt.gmin)
    if(length(q-gshel(:,j)).lt.leng) then
       qin=q-gshel(:,j)
       leng=length(qin)
    endif
 enddo

 call check_inside_ws(qin,gshel,inside)
 if(inside.eq.1) then
    return
 else
    write(*,4)'fold_in_bz_new: starting q=',q
    write(*,4)'fold_in_bz_new: final  qin=',qin
    write(*,*)'fold_in_bz_new: still not inside FBZ!'
    write(ulog,*)'fold_in_bz_new: still not inside FBZ!'
    stop
 endif

! call reduce(q,(gshel(:,1)),(gshel(:,2)),(gshel(:,3)),w)
! do i=1,3
!    w(i)=mod(w(i),1d0)-0.5d0
! enddo
! qin=w(1)*(gshel(:,1))+w(2)*(gshel(:,2))+w(3)*(gshel(:,3))

! do i=1,26

!!    qin = q+gshel(3,i)
!    qin = v2a(w)+gshel(3,i)
!    call check_inside_ws(qin,gshel,inside)
!    if(inside.eq.1) then
!        return
!    endif
! enddo

! qin=1d50
! write(ulog,*)'FOLD IN FBZ:  DONE! ==============================='

3 format(3(i6),9(2x,f8.4))
4 format(a,9(1x,g11.4))

 end subroutine fold_in_bz_new
!----------------------------------------------
  subroutine map_ibz2(nk,kp,mapi,ni)
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
  end subroutine map_ibz2

end module kpoints