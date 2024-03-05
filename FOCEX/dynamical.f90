!===========================================================

 subroutine set_band_Structure 
 use svd_stuff, only : map
 use kpoints
 use atoms_force_constants, only : natom_prim_cell
 use ios, only : ulog , uibs, uband, ugrun, write_out
 use eigen !, only : ndyn
 use constants, only : r15,ee
 use lattice, only : volume_r0,g0ws26
 use mech
 use fourier, only : nrgrid,rgrid,rws_weights !,ggrid,nggrid !,gws_weights
 implicit none
 integer i,j,k,uio
 real(r15) elastic(6,6) ,compliance(6,6), c1(6,6),q2(3),sf,dsf(3) ! coord(3,3) ,
 real(r15), allocatable :: foldedk(:,:)

  ndyn = 3*natom_prim_cell

  write(*,*)'Entered set_band_structure '
  write(ulog,*)'map(ntotind)=',map(:)%ntotind
!---------------------------------------------------------------
! do a band structure calculation along the symmetry directions

!  call cubic_on_gconv(coord)
  
4 format(a,3(1x,f9.4),1x,g13.6)
  do i=1,5
  do j=1,5
  do k=1,5
     q2=i*gs1+j*gs2+k*gs3 
     call structure_factor_recip(q2,nrgrid,rgrid,rws_weights,sf,dsf)
     write(6,4)'dyncoul: qred,sf(q)=',matmul(transpose(prim_to_cart),q2)/(2*pi),sf
  enddo
  enddo
  enddo

  call make_kp_bs

  call allocate_eig_bs(ndyn,nkp_bs,ndyn) !3rd index is for eivecs needed for gruneisen

! now output eigenval is directly frequencies in 1/cm
  write(*,*)' calling get_freqs'
  call enforce_asr_phi
  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,ndyn,eigenvec_bs,veloc_bs)

  write(*,*)' Writing eigenvalues'
  open(uband,file='bs_freq.dat')
  call write_all_eigenvalues(nkp_bs,dk_bs,kp_bs,eigenval_bs,veloc_bs,ndyn,uband)
  close(uband)

  uio=321
  open(uio,file='mech.dat')
  allocate(phi(ndyn,ndyn),xi(ndyn,3,3),qiu(ndyn,3,3),zeta(ndyn,ndyn,3),teta(ndyn,ndyn,3) &
&            ,y0(ndyn),pi0(ndyn),qiuv(ndyn,6)) 
  if(ndyn.gt.3) then ! only makes sense for non-bravais lattices
    allocate(gama(ndyn-3,ndyn-3))
    call get_phi_zeta_Xi(uio) !ndyn,atld0,gama,phi,zeta,teta,xi,qiu,uio)
    call residuals (uio) !ndyn,xi,zeta,phi,gama,sigma0,y0,pi0,uio)
    call mechanical(elastic,compliance,uio) !ndyn,atld0,sigma0,phi,zeta,xi,qiu,gama,elastic,uio)
!   call mechanical2(atld0,elastic,ndyn,xi,teta,phi,gama)
  else
    call mechanical0(elastic,atld0,uio)
  endif
  close(uio)
  write(*,*)' calling gruneisen'
  open(ugrun,file='bs_grun.dat')
  call gruneisen(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,eigenvec_bs,ugrun,grun_bs)
  close(ugrun)

  call deallocate_eig_bs
  call deallocate_kp_bs

!  call sound_speeds (elastic)
  call write_out(ulog,' Elastic tensor from Vgr ',elastic)

! now take kpoints within the IBZ for lattice summations
  call read_latdyn
  call make_kp_reg_ws(nc,g01,g02,g03,shift,kpt,wkt)
  allocate(foldedk(3,nkc))
  call fold_in_WS_BZ(nkc,kpt,g0ws26,foldedk)   ! also calls get_weights

  call allocate_eig_ibz(ndyn,nibz) 

! now output eigenval is directly frequencies in 1/cm
  write(*,*)' calling get_freqs'
  call get_frequencies(nibz,kibz,dk_bs,ndyn,eivalibz,ndyn,eivecibz,velocibz)
  call gruneisen(nibz,kibz,dk_bs,ndyn,eivalibz,eivecibz,ugrun,grunibz)

  open(utherm,file='thermo.dat')
  write(utherm,*)'# tempk,  etot,  free,  pres,  cv'
  do i=1,21
     tempk=(i-1)*50
     if(i.eq.1) tempk=10
     call thermo(nibz,wibz,ndyn,eivalibz,grunibz,tempk,etot,free,pres,cv)
     write(utherm,7) tempk,etot,free,pres,cv
  enddo
  close(utherm)

 end subroutine set_band_Structure
!===========================================================
 subroutine get_frequencies(nkp,kp,dk,ndn,eival,nv,eivec,vg)
! also outputs the group velocities from HF theorem (exact) in units of c_light
 use params
 use constants
 use born
 use geometry
 use ios, only : uband
 use eigen, only : mysqrt
 implicit none
 integer, intent(in) :: nkp,ndn,nv ! no of wanted eivecs
 real(r15), intent(in) :: kp(3,nkp),dk(nkp)
 real(r15), intent(out):: eival(ndn,nkp),vg(3,ndn,nkp)
 complex(r15), intent(out) :: eivec(ndn,ndn,nkp)
 integer i,j
! integer, allocatable :: mp(:) !Map for ascending band sorting
! real(r15), allocatable :: eivl(:)
! complex(kind=8), allocatable:: eivc(:,:),eivec_tmp(:,:)
 integer mp1(ndn,nkp) !,mp(ndn) !Map for projection band sorting
 real(r15) eival_tmp(ndn),vg_tmp(3,ndn)   !Temp matrices for reordering eigenvalues
 complex(r15) eivec_tmp(ndn,ndn) !,eivc(ndn,ndn)

! allocate(mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3))

! open files and write the group velocities
  open(330,file='veltest.dat')

 kloop: do i=1,nkp

    write(*,*)' Calling finitedif for k# ',i
    call finitedif_vel(kp(:,i),ndn,vg(:,:,i),eival(:,i),eivec(:,1:nv,i)) ! don't need group velocities
!   call get_freq(kp(:,i),ndn,vg(:,:,i),eival(:,i),eivec(:,1:nv,i) )

 enddo kloop


! Band structure projection sorting, added by Bolin
!  write(*,*) "Projection Band Sorting for band structures!"
!  call band_sort_bs(nkp,ndn,kp,eivec,mp1)
!  open(uband,file='bands.dat')
!  write(uband,*)'# i,dk(i),kp(:,i), freq(i),i=1,nband ; velocity_al(i),i=1,nband ; |velocity(i)| '
!  do i=1,nkp
!     eival_tmp=eival(:,i)
!     eivec_tmp=eivec(:,:,i)
!     vg_tmp=vg(:,:,i)
!     do j=1,ndn
!        eival(j,i)=eival_tmp(mp1(j,i))
!        eivec(:,j,i)=eivec_tmp(:,mp1(j,i))
!        vg(:,j,i)=vg_tmp(:,mp1(j,i))
!     end do
!     write(uband,2)i,dk(i),kp(:,i),(cnst*mysqrt(eival(j,i)),j=1,ndn)  &
!  &                   ,(vg(:,j,i),length(vg(:,j,i)),j=1,ndn)
!  end do
!  close(uband)


 2 format(i5,1x,4(1x,f8.3),999(1x,g11.4))
 3 format(a,i5,9(1x,f19.9))
 4 format(99(1x,2(1x,g10.3),1x))
 5 format(a,i5,99(1x,2(1x,g10.3),1x))
 6 format(2i6,2x,99(1x,f9.3))
 7 format(i6,2x,3(1x,f7.3),1x,99(1x,f7.3))

 end subroutine get_frequencies
!=====================================================
  subroutine finitedif_vel(q0,ndn,vgr,evl0,evc0)
! calculates the group velocities in units of c_light from finite difference
! it would give zero near band crossings, thus HF is a better way to do it
  use constants
  use lattice, only : prim_to_cart,volume_g0
  use geometry
  implicit none
  integer, intent(in) :: ndn
  real(r15),intent(in) :: q0(3)
  real(r15),intent(out) :: vgr(3,ndn),evl0(ndn)
  complex(r15),intent(out) :: evc0(ndn,ndn)
  integer i,l
  real(r15) q1(3),dq,v2(3,ndn),evlp(ndn),evlm(ndn) ,rand(3),qr(3)
  complex(r15) evct(ndn,ndn)

  call random_number(rand)
!  write(*,4)'FD: rand=',rand
  dq=1d-4 *volume_g0**0.33333333  
  rand=(2*rand-1)*dq
  qr=q0+ rand  ! add a small random # to break degeneracies

! call get_freq(q0,ndn,v2,evl0,evc0)
  call get_freq(qr,ndn,vgr,evl0,evc0)

! return

  q1=qr
  do i=1,3
     q1(i)=qr(i)+dq
     call get_freq(q1,ndn,v2,evlp,evct)  ! v2 used here as dummy
     q1(i)=qr(i)-dq
     call get_freq(q1,ndn,v2,evlm,evct)
     v2(i,:)=(evlp-evlm)/2/dq / 2/sqrt(abs(evl0)) *cnst*1d-8 *2*pi *c_light
     q1(i)=qr(i) ! back to qr
  enddo
!  write(30,6)'q,dq,om=',q0,dq,evl0
  write(330,4)'#la , v_FD ; v_PERT *** FOR q_red= ',matmul(transpose(prim_to_cart),qr)/(2*pi),' *** qcart=',qr
  do l=1,ndn
     write(330,5)'la ',l,v2(:,l), length(v2(:,l)), vgr(:,l),length(vgr(:,l))
  enddo


 4 format(a,3(1x,f7.4),a,3(1x,g11.4))
 5 format(a,i5,2(5x,f11.3,3(1x,f11.3)))
 6 format(a,3(1x,f9.3),1x,g12.5,9(1x,f9.4))

  end subroutine finitedif_vel
!===========================================================
 subroutine get_freq(kp,ndn,vg,eival,eivec)
!! given a kpoint kp, it sets up the dynamical matrix and diagonalizes it
!! eival is in eV/uma/Ang^2 ; group velocity calculated from Hellman-Feynmann theorem
 use params
 use ios
 use constants
 use geometry, only : is_integer,v2a
 use born
 use lattice, only : prim_to_cart,r01,r02,r03
 implicit none
 integer, intent(in) :: ndn
 real(r15), intent(in) :: kp(3)
 real(r15), intent(out):: eival(ndn),vg(3,ndn)
 complex(r15), intent(out) :: eivec(ndn,ndn)
 integer j,l,ier,nd2,al,mp(ndn) 
 integer ,save :: nksave
! integer, allocatable :: mp(:)
! real(r15), allocatable :: eivl(:)
 real(r15) absvec,k_prim(3),eivl(ndn),kdotr,vg2(3,ndn)
! complex(kind=8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:),temp(:,:)
 complex(r15) dynmat(ndn,ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3),temp(ndn,ndn)

 nd2 = min(ndn,12)

nksave=nksave+1
!write(*,3)'get_freq for ik,k=',nksave,kp

!   write(*,*) ' calling set_dynamical_matrix'
    call set_dynamical_matrix(kp,dynmat,ndn,ddyn)
!      do l=1,ndn
!         write(*,8)(dynmat(l,j),j=1,nd2)
!      enddo

!   write(*,*) ' calling nonanal'
    call nonanal(kp,dynmat,ndn,ddyn)
!      do l=1,ndn
!         write(*,8)(dynmat(l,j),j=1,nd2)
!      enddo

! is k a reciprocal lattice vector?
    k_prim=matmul(transpose(prim_to_cart),kp)/(2*pi)
    if (is_integer(k_prim(1)) .and. &
     &  is_integer(k_prim(2)) .and. &
     &  is_integer(k_prim(3)) ) then
       write(ulog,3)' calling check_asr for kp=',nksave,kp,k_prim
!      call check_asr(ndn,dynmat)  ! checks and enforces asr
    endif

    if (verbose) then
       write(ulog,3)' ===================================================='
       write(ulog,3)' THE DYNAMICAL MATRIX for KP# (reduced,cart)=',nksave,k_prim,kp(:)
       do l=1,ndn
          write(ulog,8)(dynmat(l,j),j=1,nd2)
       enddo
       write(ulog,3)' ===================================================='
!      do al=1,3
          write(ulog,3)' THE k-DERIVATIVE OF THE DYNAMICAL MATRIX =',al
          do l=1,ndn
             write(ulog,8)(ddyn(l,j,3 ),j=1,nd2)
          enddo
!      enddo
    endif

! store dynmat in temp
    temp=dynmat
    call diagonalize(ndn,temp,eivl,ndn,eivc,ier)

! sort eivals in ascending order
    call sort(ndn,eivl,mp,ndn)
!   write(ulog,6)'map=',mp

    if (ier.ne.0) then
      write(   *,*)' ier=',ier
      write(ulog,*)' ier=',ier, ' EIGENVALUES ARE...'
      write(ulog,4) eivl
    endif

    do j=1,ndn
       eival(j) = eivl(mp(j))  ! so that all frequencies are positive
       eivec(:,j) = eivc(:,mp(j))
    enddo

! if pure imaginary make eigenvectors real as they can be multiplied by any arbitrary number
!   do l=1,ndn
!      absvec = sum(abs(real(eivec(:,l))))
!      if (absvec .lt. 1d-3) then
!         eivec(:,l)=-dcmplx(0,1)*eivec(:,l)
!      endif
!   enddo

  if (verbose) then

    do l=1,ndn
        write(ulog,3)' GET_FREQ:-----------  BAND ',l,eival(l)
        do j=1,ndn/3
           write(ulog,5)'atom ',j,eivec(3*(j-1)+1,l),eivec(3*(j-1)+2,l),eivec(3*(j-1)+3,l)
        enddo
    enddo

! now calculate the square of frequencies based on dynmat
     temp = matmul(dynmat,eivec)
     temp = matmul(transpose(conjg(eivec)),temp)
     absvec = sum(cdabs(temp(:,:)))
     write(ulog,*)' TEST: below square of eivals should appear in the diagonal if e.D.e is right'
     do j=1,ndn
       write(ulog,9)'e.D.e=',j, temp(j,j)
       do l=1+j,ndn-1
          if (abs(temp(l,j)).gt.1d-5*absvec) then
             write(*,7)'l,j,[e.D.e](l,j)=',l,j,temp(l,j)
             stop
          endif
       enddo
     enddo
!    write(ulog,*)' three components of ddyn are '
!    do al=1,3
!    do j=1,ndn
!      write(ulog,9)'al=',al, ddyn(j,:,al)
!    enddo
!    enddo

  endif

! now calculate the group velocities from Hellman-Feynmann formula(for each component al=1,3)
! V_pert(k,l)= < e_l(k) | dD/dk | e_l(k) > / 2 om(k,l)
!   write(*,*) ' calculating group velocities'
    do al=1,3
      temp = matmul(transpose(conjg(eivec)),matmul(ddyn(:,:,al),eivec))
      call diagonalize(ndn,temp,vg2(al,:),ndn,eivc,ier)
      do l=1,ndn
!     vg(al,l)= dble( dot_product(conjg(eivec(:,l)),matmul(ddyn(:,:,al),eivec(:,l))) )
         vg(al,l)=vg2(al,l)/2/sqrt(abs(eival(l))+1d-24)*cnst*1d-10*100*2*pi*c_light
      enddo
    enddo
 if (verbose) then
    write(ulog,*)' three components of the velocity of all bands in m/s are '
    do al=1,3
      write(ulog,5)'alpha,v_alpha(lambda) (pert)=',al,vg(al,:)
    enddo
 endif

! deallocate(eivl,eivc,dynmat,mp,ddyn,temp)

 3 format(a,i5,9(1x,f11.4))
 4 format(99(1x,2(1x,g10.3),1x))
 5 format(a,i5,99(1x,2(1x,g10.3),1x))
 6 format(a,99(1x,i5))
 7 format(a,2i6,2x,99(1x,f7.3))
 8 format(99(1x,2(1x,f9.3),1x))
 9 format(a,i5,99(1x,2(1x,f9.3),1x))

 end subroutine get_freq
!==========================================================
 subroutine set_dynamical_matrix(kpt,dynmat,ndim,ddyn)
 use params
 use lattice
 use kpoints
 use atoms_force_constants
 use geometry
 use svd_stuff
 use constants, only : r15
 implicit none
 integer, intent(in) :: ndim
 complex(r15), intent(out) :: dynmat(ndim,ndim) ,ddyn(ndim,ndim,3)
 real(r15), intent(in) :: kpt(3)
 complex(r15) junk
 real(r15) mi,mj,nrm1,rr(3) !,delt(3)
 integer tau,j,taup,al,be,i3,j3,t,cnt2,ti,ired,g

4 format(a,4(1x,i5),9(2x,f9.4))
9 format(i6,99(1x,f9.3))

! delt = wshift !*wshift
 ddyn   = dcmplx(0d0,0d0)
 dynmat = dcmplx(0d0,0d0)
 do tau=1,natom_prim_cell
!   write(*,3)'tau,ri0=',tau,atom0(tau)%equilibrium_pos ,atompos(:,tau)
 do al=1,3
    i3 = al+3*(tau-1)
    mi = atom0(tau)%mass
! write(ulog,*) 'i,al,mass=',tau,al,i3,mi
    cnt2=0   ! counter of ntotind, to find position of each term in the array fc2
    gloop: do g=1,map(2)%ngr

       if(keep_grp2(g).ne.1) cycle gloop
       if(g.gt.1) then
         if (keep_grp2(g-1).eq.1) cnt2=cnt2+map(2)%ntind(g-1)
       endif
       tiloop: do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
          ired=cnt2+ti + map(1)%ntotind  ! the latter is the size of fc1
          if(cnt2+ti .gt. size_kept_fc2) then
             write(*,*)'SET_DYNMAT: size of fc2 exceeded, ired=',ired
             stop
          endif

          tloop: do t=1,map(2)%nt(g)
             if ( tau .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
             be = map(2)%gr(g)%ixyz(2,t)
             j  = map(2)%gr(g)%iat(2,t)
             taup = iatomcell0(j)    ! atom_sc(j)%cell%tau  is incorrect
             mj = atom0(taup)%mass
             j3 = be+3*(taup-1)
             rr = atompos(:,j) - atompos(:,taup) ! this is just R 
! k1 2-22-24  this is ith a different phase
!            rr = atompos(:,j) - atompos(:,tau)   ! this is R+taup-tau
! k1 2-22-24  this is ith a different phase
             junk = fcs(ired) * map(2)%gr(g)%mat(t,ti) * exp(ci*(kpt .dot. rr))/sqrt(mi*mj)
             dynmat(i3,j3) = dynmat(i3,j3) + junk
             ddyn(i3,j3,:) = ddyn(i3,j3,:) + junk*ci*rr(:)
          enddo tloop
       enddo tiloop
    enddo gloop
 enddo
 enddo

! if (verbose) then
!  write(ulog,*)'SET_DYNAMICAL_MATRIX: d(dyn)/dk is:'
!  do t=1,3
!     write(ulog,*)'=======component of v ',t
!     do j=1,ndim
!        write(ulog,9) j, ddyn(j,:,t)
!     enddo
!  enddo
! endif

3 format(a,i5,9(f8.3))
5 format(a,5i5,9(f8.3))
 nrm1 = maxval(cdabs(dynmat))  !sum(cdabs(dynmat(:,:)))/(ndim*ndim)
! make sure it is hermitian
 do t=1,ndim
    if (abs(aimag(dynmat(t,t))) .gt. 1d-5*nrm1) then !abs(real(dynmat(t,t))) ) then
       write(ulog,*)' dynmat is not hermitian on its diagonal'
       write(ulog,*)' diagonal element i=',t,dynmat(t,t)
       write(ulog,*)' setting its imaginary part to zero!'
!      stop
!   else
       dynmat(t,t) = dcmplx(real(dynmat(t,t)),0d0)
    endif
    do j=t+1,ndim-1
      if (abs(aimag(dynmat(t,j))+aimag(dynmat(j,t))) .gt. 1d-5*nrm1 ) then
        write(ulog,*)' dynmat is not hermitian in AIMAG of its off-diagonal elts'
        write(ulog,*)' off-diagonal element i,j=',t,j,aimag(dynmat(t,j)),aimag(dynmat(j,t))
        write(ulog,*)' compared to max(abs(dynmat))=',nrm1
!       stop
      elseif(abs(real(dynmat(t,j))-real(dynmat(j,t))) .gt. 1d-5*nrm1 ) then
         write(ulog,*)' dynmat is not hermitian in REAL of its off-diagonal elts'
         write(ulog,*)' off-diagonal element i,j=',t,j,real(dynmat(t,j)),real(dynmat(j,t))
         write(ulog,*)' compared to max(abs(dynmat))=',nrm1
!        stop
      else
    endif
!   enforcing it to be hermitian in anycase!!'
    mi=(aimag(dynmat(t,j))-aimag(dynmat(j,t)))/2
    mj=(real (dynmat(t,j))+real (dynmat(j,t)))/2
    dynmat(t,j) = dcmplx(mj, mi)
    dynmat(j,t) = dcmplx(mj,-mi)
    do i3=1,3  ! ddyn is anti hermitian (ddyn^*T = -ddyn)
       mi=(aimag(ddyn(t,j,i3))+aimag(ddyn(j,t,i3)))/2
       mj=(real (ddyn(t,j,i3))-real (ddyn(j,t,i3)))/2
       ddyn(t,j,i3) = dcmplx(mj, mi)
       ddyn(j,t,i3) = -dcmplx(mj,-mi)
    enddo
  enddo
 enddo


 end subroutine set_dynamical_matrix
!===========================================================
 subroutine nonanal(q,dynmat,ndim,ddyn)
!! adds the non-analytical Coulomb interaction between Born charges to
!! the dynamical matrix
 use constants
 use lattice
 use atoms_force_constants
 use born
 use geometry
 use ewald
 use ios, only : ulog,write_out

 integer, intent(in) :: ndim
 complex(r15), intent(inout) :: dynmat(ndim,ndim),ddyn(ndim,ndim,3)
 real(r15), intent(in) :: q(3)
 real(r15) ma,mb  
 integer na,nb,al,be  
 real(r15) dyn_coul(3,3) ,ddn(3,3,3) 

  if(born_flag.le.0) return ! only for BF > 0 the non-analytical term is added to the dynamical matrix 

!  dynmat=0; ddyn=0 they are already initialized in set_dynmat
  do na=1,natom_prim_cell
     ma = atom0(na)%mass
     do nb=1,natom_prim_cell
        mb = atom0(nb)%mass
        call dyn_coulomb(na,nb,q,dyn_coul,ddn)

!    call write_out(ulog,' NONANAL: dyn_coul ',dyn_coul)
        do al = 1,3
        do be = 1,3
           dynmat(al+3*(na-1),be+3*(nb-1)) = dynmat(al+3*(na-1),be+3*(nb-1))+ dyn_coul(al,be)  /sqrt(ma*mb)
           ddyn  (al+3*(na-1),be+3*(nb-1),:) = ddyn(al+3*(na-1),be+3*(nb-1),:)   + ddn(al,be,:)/sqrt(ma*mb)
        enddo
        enddo

     enddo
  enddo


 end subroutine nonanal
!===========================================================
 subroutine dyn_coulomb(tau,taup,q,dyn,ddyn)
!! calculates the q component of the Coulomb force to subtract from fourier transform of total force
! sdyn(tau,taup) = [z(tau)^al,ga.q^ga]  [z(taup)^be.de.q^de] * sf(q) &
!    &          1/eps0/(q.epsr.q)/volume_r0
! DOES NOT INCLUDE THE MASS DENOMINATOR
 use atoms_force_constants
 use geometry
 use lattice, only :  volume_r0 , prim_to_cart ,g0ws26
 use fourier, only : nrgrid,rgrid,rws_weights ,ggrid,nggrid !,gws_weights
 use born, only : epsil !,born_flag
 use constants , only : eps0,ee ,pi,r15
 use ios, only : ulog,write_out
 use params, only : verbose
 implicit none
 integer, intent(in) :: tau,taup
 real(r15), intent(in) :: q(3)
 real(r15), intent(out) :: dyn(3,3),ddyn(3,3,3)  ! element_(tau,taup)(q) of the Coulomb dynamical mat
 integer al,be,i
 real(r15) zqa(3),zqb(3),coef,mt(3,3),sf,dsf(3),qeq,dqeq(3),q2(3)


    dyn=0; ddyn=0

 coef=ee*1d10/eps0/volume_r0 ! to convert 1/ang^2 to eV/ang , includes cancellation of 4pi
! call write_out(6,' DYN_COULOMB: epsil ',epsil)
! denom=sum( (/ (epsil(i,i),i=1,3)/) )   ! this is the trace function
! call write_out(6,' DYN_COULOMB: trace(epsil) ',denom)
 q2=q
 if(length(q).lt.1d-6) then
    q2(1)=1d-6
    q2(2)=1d-6
    q2(3)=1d-6
 endif

   qeq=dot_product(q2,matmul(epsil,q2))
!write(6,*)'Dyn_coulomb: coef,q2,qeq=',coef,q2,qeq
   dqeq=matmul((epsil+transpose(epsil)),q2)

   if (fc2flag.eq.0)then
      call structure_factor_recip(q2,nrgrid,rgrid,rws_weights,sf,dsf)
   else
       sf=exp(-2*length(q2)**2*volume_r0**0.666666)
       dsf(:)=-4*q2(:)*volume_r0**0.666666*sf 
       do i=1,26
          sf=sf+exp(-2*length(q2+g0ws26(:,i))**2*volume_r0**0.666666)
          dsf=dsf-4*(q2+g0ws26(:,i))*volume_r0**0.666666+exp(-2*length(q2+g0ws26(:,i))**2*volume_r0**0.666666)
       enddo
   endif

   zqa=matmul(atom0(tau)%charge,q2)
   zqb=matmul(atom0(taup)%charge,q2)
!call write_out(6,' DYN_COULOMB: sf  ',sf )
!call write_out(6,' DYN_COULOMB: dqeq ',dqeq)
!call write_out(6,' DYN_COULOMB: zqa ',zqa)
!call write_out(6,' DYN_COULOMB: zqb ',zqb)

   do al=1,3
   do be=1,3
      dyn(al,be)=zqa(al)*zqb(be)/qeq * coef * sf
      ddyn(al,be,:) = ((atom0(tau)%charge(al,:)*zqb(be)+zqa(al)*atom0(taup)%charge(be,:))*sf +  &
 &                        zqa(al)*zqb(be)*(dsf(:)-sf*dqeq(:)/qeq))/qeq*coef
   enddo
   enddo

 if(length(q).lt.1d-6) then
! if(verbose) then
     write(6,4)' NA 3x3 block in DYN_COULOMB: tau,taup, sf,q,q_red=',tau,taup,sf,q,matmul(transpose(prim_to_cart),q)/(2*pi)
     do al=1,3
        write(6,3)(dyn(al,be),be=1,3)
     enddo  
! endif
 endif

3 format(3(1x,f9.4))
4 format(a,2i4,9(1x,f9.4))

 end subroutine dyn_coulomb
!============================================================
 subroutine gruneisen(nkp,kp,dk,ndn,eivl,eivc,ugr2,grn)
! takes the eigenvalues (w^2) and eigenvectors calculated along some
! crystalline directions and calculates the corresponding mode
! gruneisen parameters
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use eigen
 use constants
 implicit none
 integer, intent(in) :: nkp,ndn,ugr2
 real(r15), intent(in) :: kp(3,nkp),dk(nkp),eivl(ndn,nkp)
 complex(r15), intent(in) :: eivc(ndn,ndn,nkp)
 complex(r15), intent(out) :: grn(ndn,nkp)
 integer ik,tau,la,al,be,ga,j,k,taup,k0,ta1,ta2,t,cnt3,g,ti,ired
 real(r15) mi,mj,rr3(3),rr2(3),qq(3),qdotr,omk !,mysqrt
 complex(r15) zz,one,term

 one = dcmplx(1d0,0d0)
! write(ugr2,*)'# la, nk , dk(nk) ,     kp(1:3,nk)     , om(la,nk),Re(gruneisen(la,nk)), Im(gr))'
 write(ugr2,*)'# nk , dk(nk) ,     kp(1:3,nk)     , om(la,nk),Re(gruneisen(la,nk)) '
  do ik=1,nkp
     qq(:) = kp(:,ik)
  do la=1,ndn
!    write(ulog,6)'Starting gruneisen parameters for kpoint & band# ',ik,la,kp(:,ik)
!    write(ulog,*)' fcs_3(t),term(2),6*eival(la,i),rr3(3),grun(la,i)'
     grn(la,ik) = 0
     omk = mysqrt(eivl(la,ik))*cnst
     do tau=1,natom_prim_cell
        mi = atom0(tau)%mass
        cnt3=0
        gloop: do g=1,map(3)%ngr
           if(g.gt.1) then
              cnt3=cnt3+map(3)%ntind(g-1)  ! cnt3 counts all the terms in the previous groups
           endif
           tloop: do t=1,map(3)%nt(g)  !nterms(3)

              if ( tau .ne. map(3)%gr(g)%iat(1,t) ) cycle tloop
              al = map(3)%gr(g)%ixyz(1,t)  !ixyzterm_3(1,t)
              be = map(3)%gr(g)%ixyz(2,t)  !ixyzterm_3(2,t)
              ga = map(3)%gr(g)%ixyz(3,t)  !ixyzterm_3(3,t)
              j  = map(3)%gr(g)%iat(2,t)   ! this atompos index
              k  = map(3)%gr(g)%iat(3,t)
              taup = iatomcell0(j)     ! atom_sc(j)%cell%tau  is incorrect
              k0 = iatomcell0(k)     !atom_sc(k)%cell%tau is incorrect
              mj = atom0(taup)%mass
              ta1= al + 3*(tau-1)
              ta2= be + 3*(taup-1)
              rr3(:) = atompos(:,k)                  ! R"+tau"
              rr2(:) = atompos(:,j) - atompos(:,taup)  ! R consistent with dynmat which uses j-j0
! K1 2-9-24 this one is consistnt with exp(iq.(r+taup-tau))
!             rr2(:) = atompos(:,j) - atompos(:,tau)  ! R consistent with dynmat which uses j-i0
! K1 2-9-24 this one is consistnt with exp(iq.(r+taup-tau))
              qdotr =  ( qq .dot. rr2)
              zz = cdexp( ci * qdotr ) 
              do ti=1,map(3)%ntind(g)  ! index of independent terms in that group g
                 ired=cnt3+ti + map(1)%ntotind + size_kept_fc2 !map(2)%ntotind
!! term = - ampterm_3(t)*fcs_3(igroup_3(t))*zz*eivc(ta1,la,i)*conjg(eivc(ta2,la,i))*rr3(ga)/sqrt(mi*mj)
                 term = -zz * eivc(ta2,la,ik)*conjg(eivc(ta1,la,ik))/sqrt(mi*mj) &
         &             * (fcs(ired) * rr3(ga)*map(3)%gr(g)%mat(t,ti))
                 grn(la,ik) = grn(la,ik) + term

! if (ik.eq.nkp .and. map(3)%gr(g)%mat(t,ti).ne.0) then
!    write(*,9)'t,ti,la,g,tau,ired,R2,R3,zz,fc,term=',t,ti,la,g,tau,ired,rr2,rr3,zz,fcs(ired),term
! endif
              enddo
           enddo tloop
        enddo gloop
     enddo
     
     grn(la,ik) = grn(la,ik)/6/eivl(la,ik)
     if(abs(grn(la,ik)).gt.1d3) grn(la,ik)=0  ! substitute very large gama by zero

     if (aimag(grn(la,ik)) .gt. 1d-4) then
       write(ulog,*)' GRUNEISEN: ... has large imaginary part!! for nk,la=',ik,la,grn(la,ik)

     endif
!   write(ugr2,6)' ',ik,la,dk(ik),kp(:,ik),omk,grn(la,ik)
  enddo
  write(ugr2,8)ik,dk(ik),kp(:,ik),( mysqrt(eivl(la,ik))*cnst , la=1,ndn), &
&              (real(grn(la,ik)),la=1,ndn)

 enddo

5 format(4i7,9(1x,g11.4))
6 format(a,2i5,99(1x,g11.4))
7 format(i5,i5,i6,99(1x,g11.4))
8 format(i8,66(1x,f9.4),6(1x,f9.4))
9 format(a,6i4,66(1x,g11.4),6(1x,f9.4))
! deallocate(eivl,eivc,kg,grn)
 end subroutine gruneisen
!============================================================
 subroutine get_phi_zeta_Xi(uio) !ndn,atld0,gama,phi,zeta,teta,xi,qiu,uio)
!! calculates some matrices useful for later operations, in addition to the  "standard" elastic constants
 use ios , only: ulog, write_out
 use lattice, only : volume_r0 
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 use params, only : verbose
 use mech
 use eigen, only : ndyn
 use linalgb
 implicit none
 integer, intent(in) :: uio !,ndn
! real(r15), intent(out) :: gama(max(1,ndn-3),max(1,ndn-3)),phi(ndn,ndn),zeta(ndn,ndn,3), &
!&            xi(ndn,3,3),teta(ndn,ndn,3),atld0(3,3,3,3),qiu(ndn,3,3)
 integer tau,taup,al,be,ga,de,j,t,g,ti,s,cnt2,ired,la,nl,nc,nq
 real(r15)  rij(3),junk,gam(ndyn,ndyn),tm(max(1,ndyn-3),max(1,ndyn-3)),constr(3,3,3),am(3,3,3,3)
 real(r15)  matr(3,3),c1(6,6),c0(6,6),aux(ndyn,ndyn)

 write(  * ,*)' ********** ENTERING get_phi_zeta_x *************'
 write(ulog,*)' ********** ENTERING get_phi_zeta_x *************'

 zeta=0; phi=0; xi=0;  teta=0; atld0=0;constr=0
 do al=1,3
 do ga=1,3
 do tau=1,natom_prim_cell
 cnt2=0  ! cumulative number of independent terms up to group (g-1)
 gloop: do g=1,map(2)%ngr
    if(keep_grp2(g).ne.1) cycle
    if(g.gt.1) then
       if (keep_grp2(g-1).eq.1) cnt2=cnt2+map(2)%ntind(g-1)
    endif
    do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
       ired=cnt2+ti + map(1)%ntotind  ! index of the indep FC2 in the list fcs
       if(cnt2+ti .gt. size_kept_fc2) then
          write(*,*)'MECHANICAL: size of fc2 exceeded, ired=',ired
          stop
       endif
       tloop: do t=1,map(2)%nt(g)
          if ( tau .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
          if                                     ( ga .ne. map(2)%gr(g)%ixyz(2,t) ) cycle tloop
          j  = map(2)%gr(g)%iat(2,t)
          taup = iatomcell0(j)    ! atom_sc(j)%cell%tau  is incorrect
          rij = atompos(:,j) - atompos(:,tau)
          junk = fcs(ired)* map(2)%gr(g)%mat(t,ti)  ! that is phi(al,ga,tau,j)
          nl=al+3*(tau-1) ; nc=ga+3*(taup-1)    ! dynmat dimensions
          phi(nl,nc) =phi(nl,nc)+junk ! = sum_R fc2(0,tau;R,Taup)
          do la=1,3
             constr(al,ga,la)=constr(al,ga,la) + junk * rij(la)  ! = sum_R,tau,taup fc2(0,tau;R,Taup) (R+taup-tau)_la
          enddo
          do la=1,3
             teta(nl,nc,la)=teta(nl,nc,la)+junk * atompos(la,j)  ! = sum_R fc2(0,tau;R,Taup) (R+taup)_la
             zeta(nl,nc,la)=zeta(nl,nc,la)+junk * rij(la)  ! = sum_R fc2(0,tau;R,Taup) (R+taup-tau)_la
          enddo
          do be=1,3
          do de=1,3
             atld0(al,be,ga,de)=atld0(al,be,ga,de) - junk*rij(be)*rij(de)/2
          enddo
          enddo
       enddo tloop
    enddo
 enddo gloop
 enddo
 enddo
 enddo


 qiu=0
 do la=1,ndyn
    do al=1,3
    do ga=1,3
    do tau=1,natom_prim_cell
       qiu(la,al,ga)=qiu(la,al,ga)+teta(la,3*(tau-1)+al,ga)
    enddo
    enddo
    enddo
! ensure qiu is symmetric
   if(maxval(abs(qiu(la,:,:)-transpose(qiu(la,:,:)))).gt.1d-6) call symmetrize2(3,qiu(la,:,:))
 enddo

 do la =1,3
    call write_out(ulog,' CONSTR=\sum_R,tau,taup phi_ij * R_ij (eV/Ang) = 0',constr(:,:,la))
 enddo

 junk=maxval((phi-transpose(phi))*(phi-transpose(phi)))
 if (junk.gt.1d-12) then
    call write_out(   6,' PHI NOT SYMMETRIC ',phi)
    call write_out(ulog,' PHI NOT SYMMETRIC ',phi)
    stop
 endif
 if (verbose) then
    call write_out(ulog,' PHI=sum_R phi(tau,R+taup)  (eV/Ang^2)',phi)
 endif

 do la =1,3 ! for each la, zeta(:,:,la) is ANTIsymmetric wrt first 2 indices
    if (verbose) then
       call write_out(ulog,' ZETA=sum_R phi(tau,R+taup)(R+taup-tau) (eV/Ang)',zeta(:,:,la))
    endif
    aux=zeta(:,:,la)+transpose(zeta(:,:,la))
    junk=maxval(aux*aux)
    if (junk.gt.1d-12) then
       call write_out(   6,' Zeta NOT ANTISYMMETRIC:z+z^T ',aux)
       call write_out(ulog,' Zeta NOT ANTISYMMETRIC:z+z^T ',aux)
    endif
 enddo

3 format(a,2i5,9(1x,f10.4))

! to get gama, invert phi: gam.phi=1a ---------------------------------------
  gam=phi  ! force constants and the images sum_R K(tau,R+taup)
  gama=0 ! gama.phi=1
  if(ndyn.gt.3) then
     gama=gam(4:ndyn,4:ndyn)
     call inverse_real(gama,tm,ndyn-3) 
     gam(1:3,:)=0 ; gam(:,1:3)=0
     gam(4:ndyn,4:ndyn)=tm  
     gama=tm  
     junk=maxval((gam-transpose(gam))*(gam-transpose(gam)))
     if (junk.gt.1d-12) then
        call write_out(   6,' |G-G^T|^2 ',junk)
        call write_out(ulog,' |G-G^T|^2 ',junk)
        call write_out(   6,' GAMA NOT SYMMETRIC ',gam)
        call write_out(ulog,' GAMA NOT SYMMETRIC ',gam)
        stop
     endif
     if (verbose) then
        call write_out(ulog,' Gamma (should be symmetric) ',gam)
     endif
  endif

! this part calculates xi(tau,ga;al,be) = du(tau,ga)/d eta_al,be = - Gam. qiu  ! should be symmetric under al <-> be
! xi=0 
! do al=1,3
! do be=1,3
!    do tau=1,natom_prim_cell
!    do ga=1,3
!       j=ga+3*(tau-1)
!
!       do s=1,natom_prim_cell
!          nc=al+3*(s-1)
!! Should not symmetrize; strictly speaking, the formula below is correct
!          xi(j,al,be)=xi(j,al,be)-dot_product(gam(j,:),teta(:,nc,be))
!       enddo 
!
!    enddo 
!    enddo 
! enddo 
! enddo 

!  xi=-gama*qiu ; this way both xi and qiu are symmetric wrt al<->be
 constr=0
 do s=1,ndyn
    ga=mod(s-1,3)+1
 do al=1,3
 do be=1,3
    xi(s,al,be)= -dot_product(gam(s,:),qiu(:,al,be))
    constr(al,be,ga)=constr(al,be,ga)+qiu(s,al,be)
 enddo
 enddo
    call write_out(ulog,' Symmetrized xi (Ang) ', xi(s,:,:))
 enddo
 if(maxval(abs(constr)).gt.1d-5) call write_out(ulog,' sum_tau(qiu) non zero ',constr)

! do tau=1,ndyn
!     junk=maxval((xi(tau,:,:)-transpose(xi(tau,:,:)))*(xi(tau,:,:)-transpose(xi(tau,:,:))))
!     if (junk.gt.1d-12) then
!        write(   6,5)' for tau, |Xi-Xi^T|^2 ',tau,junk
!        write(ulog,5)' for tau, |Xi-Xi^T|^2 ',tau,junk
!        call write_out(   6,' Xi NOT SYMMETRIC ',xi(tau,:,:))
!        call write_out(ulog,' Xi NOT SYMMETRIC ',xi(tau,:,:))
!!       stop
!     endif
!     call symmetrize2(xi(tau,:,:))
!     call symmetrize2(qiu(tau,:,:))
!     if (verbose) then
!        call write_out(ulog,' Symmetrized xi (Ang) ', xi(tau,:,:))
!        call write_out(ulog,' Symmetrized qiu(Ang) ',qiu(tau,:,:))
!     endif 
! enddo

! calculate A from teta
! am=0
! do al=1,3
! do be=1,3
! do ga=1,3
! do de=1,3
!    do tau=1,natom_prim_cell 
!       am(al,be,ga,de)=am(al,be,ga,de)+qiu(3*(tau-1)+al,ga,de)* atompos(be,tau)
!    enddo
! enddo
! enddo
! enddo
! enddo

! am   =am   /volume_r0*ee*1d30*1d-9
 atld0=atld0/volume_r0*ee*1d30*1d-9
 call convert_to_voigt(atld0,c0)
! call convert_to_voigt(am   ,c1)
 call write_out (ulog,' Elastic Tensor (standard term; old formula) in GPa, in voigt ',c0)
! call write_out (ulog,' Elastic Tensor AM=teta*tau in GPa, in voigt ',c1)

4 format(a,99(1x,f10.4))
5 format(a,i5, 99(1x,f10.4))
 end subroutine get_phi_zeta_Xi
!============================================================
 subroutine mechanical0(elastic,atld0,uio)
!! calculates elastic constants, only in the cubic case, from sound velocities atlong 110
!! also calculates the standard term -0.5*phi_ik R_jl R_jl (translationally invariant formulation)
!! it is called for a bravais lattice natom_prim_cell=1

 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 implicit none
 integer, intent(in) :: uio
 real(r15), intent(out) :: elastic(6,6),atld0(3,3,3,3) 
 integer tau,taup,al,be,ga,de,i,j,t,g,ti,cnt2,ired,voigt,la,nl,nc,nu,g1,g2,mu,nq,s,n1,n2,ndn
 real(r15)  rij(3),junk,c11,c12,c44,rho_SI,constr(3,3,3) &
         ,c1(6,6),c2(6,6),c3(6,6),q(3)   &
 &       ,ct(3,3,3,3),halfsum,halfdif  &
 &       ,vg(3,3*natom_prim_cell),eivl(3*natom_prim_cell),eivl0(3*natom_prim_cell)
 complex(r15) eivc(3*natom_prim_cell,3*natom_prim_cell),eivc0(3*natom_prim_cell,3*natom_prim_cell)

 write(  * ,*)' ********** ENTERING MECHANICAL0 *************'
 write(uio,*)' ********** ENTERING MECHANICAL0 *************'
 write(uio,*)' Elastic constant of a cubic crystal from group velocities '
 ndn=3*natom_prim_cell
 rho_SI=sum(atom0(:)%mass)*uma/volume_r0*1d30
 write(uio,*)'density in kg/m^3 is=',rho_SI
 q=(/1d-6,1d-6,0d0/)
 call get_freq(q,ndn,vg,eivl0,eivc0)
 call write_out(uio,'Vg from perturbation q=1d-6 ',vg)
 write(uio,*)'Length(Vg_pert) is=',(length(vg(:,j)),j=1,3)

 q=(/1d-3,1d-3,0d0/)
! from FD or perturbation group velocities
 call get_freq(q,ndn,vg,eivl,eivc); eivl=eivl-eivl0
 call write_out(uio,'Vg from PERT q=1d-3(1,1,0) ',vg)
 write(uio,3)'Length(Vg_pert) for 3 acoustic modes =',(length(vg(:,j)),j=1,3)
 write(uio,3)'Corresponding rho v^2 (GPa)           =',(rho_SI * 1d-9*length(vg(:,j))**2 , j=1,3)

3 format(a,9(1x,g11.4))

! eivl is the square of phonon frequencies eivl=omega^2 -> convert to eV
      eivl=eivl/(q.dot.q)*rho_SI  *ee/uma*1d-9
      call write_out(uio,'rho_cell*w^2 /q^2 (GPa) ',eivl)

      c44    =rho_SI*(length(vg(:,1)))**2 *1d-9  ! assumes eivc(z,1)=1
      halfdif=rho_SI*(length(vg(:,2)))**2 *1d-9
      halfsum=rho_SI*(length(vg(:,3)))**2 *1d-9 - c44
      c12=halfsum-halfdif
      c11=halfsum+halfdif
      call write_out(uio,'cij from rho vg^2 if e_1=z :',(/c11,c12,c44/))

      c44    =rho_SI*(length(vg(:,2)))**2 *1d-9  ! assumes eivc(z,2)=1
      halfdif=rho_SI*(length(vg(:,1)))**2 *1d-9
      halfsum=rho_SI*(length(vg(:,3)))**2 *1d-9 - c44
      c12=halfsum-halfdif
      c11=halfsum+halfdif
      call write_out(uio,'cij from rho vg^2 if e_2=z :',(/c11,c12,c44/))

! call write_out(uio,'GF dyn_mat eivals (eV/Ang^2) ',eivl)
  call finitedif_vel(q,ndn,vg,eivl,eivc); eivl=eivl-eivl0
  call write_out(uio,'Vg from FD q=1d-3(1,1,0) ',vg)
  write(uio,*)'Length(Vg_FD) is=',(length(vg(:,j)),j=1,3)
  write(uio,*)'rho v^2 (GPa) is=',(rho_SI * 1d-9*length(vg(:,j))**2 , j=1,3)
      c44    =rho_SI*(length(vg(:,1)))**2 *1d-9  ! assumes eivc(z,1)=1
      halfdif=rho_SI*(length(vg(:,2)))**2 *1d-9
      halfsum=rho_SI*(length(vg(:,3)))**2 *1d-9 - c44
      c12=halfsum-halfdif
      c11=halfsum+halfdif
      call write_out(uio,'CUBIC_ONLY: cij from rho vg^2 if e_1=z :',(/c11,c12,c44/))

      c44    =rho_SI*(length(vg(:,2)))**2 *1d-9  ! assumes eivc(z,2)=1
      halfdif=rho_SI*(length(vg(:,1)))**2 *1d-9
      halfsum=rho_SI*(length(vg(:,3)))**2 *1d-9 - c44
      c12=halfsum-halfdif
      c11=halfsum+halfdif
      call write_out(uio,'CUBIC_ONLY: cij from rho vg^2 if e_2=z :',(/c11,c12,c44/))

 write(uio,*)' ********************************************'


! calculation of the standard Bravais formula in atld0
 atld0=0; constr=0  ! atld0 = -0.25 sum_r,tau,taup phi(tau;R+taup)*(R+taup-tau)  ------------
 do al=1,3
 do ga=1,3
 do tau=1,natom_prim_cell
 cnt2=0  ! cumulative number of independent terms up to group (g-1)
 gloop: do g=1,map(2)%ngr
    if(keep_grp2(g).ne.1) cycle
    if(g.gt.1) then
       if (keep_grp2(g-1).eq.1) cnt2=cnt2+map(2)%ntind(g-1)
    endif
    do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
       ired=cnt2+ti + map(1)%ntotind  ! index of the indep FC2 in the list fcs
       if(cnt2+ti .gt. size_kept_fc2) then
          write(*,*)'MECHANICAL: size of fc2 exceeded, ired=',ired
          stop
       endif
       tloop: do t=1,map(2)%nt(g)
          if ( tau .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
          if( map(2)%gr(g)%ixyz(2,t) .ne. ga) cycle tloop
          j  = map(2)%gr(g)%iat(2,t)
          taup = iatomcell0(j)    ! atom_sc(j)%cell%tau  is incorrect
          rij = atompos(:,j) - atompos(:,tau)
          junk = fcs(ired)* map(2)%gr(g)%mat(t,ti)  ! that is phi(al,ga,tau,j)
          nl=al+3*(tau-1) ; nc=ga+3*(taup-1)  ! dynmat dimensions
          do la=1,3
             constr(al,ga,la)=constr(al,ga,la) + junk * rij(la)  ! = sum_R,tau,taup fc2(0,tau;R,Taup) (R+taup-tau)_la
          enddo
          do be=1,3
          do de=1,3
             atld0(al,be,ga,de)=atld0(al,be,ga,de) - junk*(rij(be)*rij(de))/2
          enddo
          enddo
       enddo tloop
    enddo
 enddo gloop
 enddo
 enddo
 enddo
 atld0=atld0/volume_r0*ee*1d30*1d-9

 do la =1,3
    call write_out(uio,' CONSTR=\sum_R,tau,taup phi_ij * R_ij (eV/Ang) = 0',constr(:,:,la))
 enddo

 call convert_to_voigt(atld0,elastic)
 call write_out (uio,' Elastic Tensor (standard term; old formula) in GPa, in voigt ',elastic)

 end subroutine mechanical0
!============================================================
 subroutine residuals(uio) !ndn,xi,zeta,phi,gam,sigma0,y0,pi0,uio)
!! given phi,gam=1/phi,and xi=d_u/d_eta, calculates residual forces pi0, stresses sigma0 and displacements y0 
!! sigma0=sum pi*tau ; y0=-gama*pi ; y_eq=y0-gama*qiu*u_eq ; a'*u_eq=(qiu*gama*pi-sigma0)
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 use linalgb, only : symmetrize2,symmetrize_res
 use mech
 use eigen, only : ndyn
 implicit none
 integer, intent(in) :: uio !,ndn 
! real(r15), intent(in) :: xi(ndn,3,3),phi(ndn,ndn),zeta(ndn,ndn,3),gam(ndn-3,ndn-3)
! real(r15), intent(out) :: sigma0(3,3),y0(3*natom_prim_cell),pi0(3*natom_prim_cell)   
 integer tau,taup,al,be,ga,de,i,j,t,g,ti,cnt2,ired,voigt,la,nl,nc,nu,s,nq,delta_k
 real(r15)  rij(3),junk,constr(3,3,3),res(3,3) ,mat2(3,3) 

 write(  * ,*)' ********** ENTERING RESIDUALS *************, ndyn=',ndyn
 write(uio,*)' ********** ENTERING RESIDUALS *************, ndyn=',ndyn

! pi0 is the residual force ------------------------------------
 pi0=0
 do tau=1,natom_prim_cell
 do al=1,3
    i=al+3*(tau-1)
    cnt2=0
    do g=1,map(1)%ngr  ! ineq. terms and identify neighbors
      if(g.gt.1) cnt2 = cnt2 + map(1)%ntind(g-1)
    do t=1,map(1)%nt(g)
       if ( map(1)%gr(g)%ixyz(1,t) .ne. al ) cycle
       if ( map(1)%gr(g)%iat(1,t) .ne. tau ) cycle
       do ti=1,map(1)%ntind(g)
          ired = cnt2+ti    ! this is the corresponding index of i in ared
          pi0(i)=pi0(i)+fcs(ired) * map(1)%gr(g)%mat(t,ti)
       enddo
    enddo
    enddo
 enddo
    write(uio,3)'Residual -force(eV/Ang): tau,pi0(tau,:)=',tau,(pi0(3*(tau-1)+al),al=1,3)
 enddo

3 format(a,i4,99(1x,f11.5))

 y0=0 ! y0 = -Gama*pi correction to equilibrium positions  at eta=0----------------------
 do tau=2,natom_prim_cell
    do al=1,3
    i=3*(tau-1)+al
    do j=4,ndyn
       y0(i)=y0(i)-gama(i-3,j-3)*pi0(j)
    enddo
    enddo
    write(uio,3)' tau,y0(tau,:)(Ang)=',tau,(y0(3*(tau-1)+al),al=1,3)
 enddo

! residual Stress tensor (under no strain) ------------------------------------
 sigma0=0
 do al=1,3
 do be=1,3
    do tau=1,natom_prim_cell
         sigma0(al,be)=sigma0(al,be)+atompos(be,tau)*pi0(3*(tau-1)+al) !  &
    enddo
 enddo
 enddo
 call write_out(uio,' sigma(eta=0) before symmetrization (eV) ',sigma0)
! should be symmetric according to rotational invariance
 call symmetrize2(3,sigma0)

 qiuv=0
 do al=1,3
 do be=al,3
    sigmav(voigt(al,be))=sigma0(al,be)
    qiuv(:,voigt(al,be))=qiu(:,al,be)
 enddo
 enddo
 do la=1,6
    sigmav(la)=sigmav(la)+dot_product(qiuv(:,la),y0)
 enddo
 call write_out(uio,'Residual stress sigma(eta) (eV) ',sigmav)
 sigma0 = sigma0/volume_r0*1d30*ee*1d-9
 call write_out(uio,' sigma(eta=0) after symmetrization (GPa) ',sigma0)
! enforces mat2(a,b)-mat2(b,a)=res(a,b)-res(b,a); consider mat=mat2-res ; symmetrize mat then mat2=sym(mat)+res
! check below
 s=0
 do al=1,3
    do tau=1,natom_prim_cell
       do be=1,3
       do ga=1,3
          junk=0
          do taup=1,natom_prim_cell
            junk=junk+teta(3*(tau-1)+al,3*(taup-1)+be,ga)-teta(3*(tau-1)+al,3*(taup-1)+ga,be)
          enddo
          junk=junk + (pi0(3*(tau-1)+be)*delta_k(al,ga)-pi0(3*(tau-1)+ga)*delta_k(al,be))
          if(abs(junk).gt.1d-5) then
             write(ulog,*)'GET_PHI_XI: rot invce violation in theta:tau,al,be,ga,rot=',tau,al,be,ga,junk
             s=1
          endif
       enddo
       enddo
    enddo
 enddo

! if(s.ne.0) then 
!   write(*,*)' symmetrizing now teta '
!! Impose rotational invariance on teta: sum_taup teta(tau,al;taup,be;ga) + pi(taup,be) delta(al,ga) symm in be<->ga
!    do al=1,3
!    do tau=1,natom_prim_cell
!       la=al+3*(tau-1)
!       do be=1,3
!       do ga=1,3
!          res(be,ga)= -( pi0(3*(tau-1)+be)*delta_k(al,ga)-pi0(3*(tau-1)+ga)*delta_k(al,be) )
!          mat2(be,ga)=0
!          do taup=1,natom_prim_cell
!            mat2(be,ga)=mat2(be,ga)+teta(3*(tau-1)+al,3*(taup-1)+be,ga)-teta(3*(tau-1)+al,3*(taup-1)+ga,be)
!          enddo
!       enddo
!       enddo
!       call symmetrize_res(mat2,res)
!    enddo
!    enddo
!
! endif

 end subroutine residuals
!===================================================
 subroutine mechanical(elastic,compliance,uio) !ndn,atld1,sigma0,phi,zeta,xi,qiu,gama,elastic,uio)
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 use mech
 use eigen, only : ndyn
 implicit none
 integer, intent(in) :: uio !,ndn 
! real(r15), intent(in) :: xi(ndn,3,3),phi(ndn,ndn),zeta(ndn,ndn,3),sigma0(3,3), &
!&                         qiu(ndn,3,3),gama(ndn-3,ndn-3)
 real(r15), intent(out) :: elastic(6,6),compliance(6,6)
 integer tau,taup,al,be,ga,de,i,j,t,g,ti,cnt2,ired,la,nl,nc,nu,s,nq,delta_k,ier
 real(r15) c1(6,6),c2(6,6),c3(6,6),cq(6,6) ,atld2(3,3,3,3),atld3(3,3,3,3),qgq(3,3,3,3)


 call convert_to_voigt(atld0,c1)
 call write_out (uio,' Elastic Tensor (first term in GPa) in voigt ',c1)

 qgq=0 

 atld2=0    ! cross terms : phi*(R+tau)*Xi = zeta*delta (Xi) needs symmtrization --------------------
 do al=1,3
 do be=1,3
 do ga=1,3
 do de=1,3

    qgq(al,be,ga,de)=dot_product(qiu(4:ndyn,al,be),matmul(gama,qiu(4:ndyn,ga,de)))
!   do tau=1,natom_prim_cell
!      nl=al+3*(tau -1) ; nu=ga+3*(tau -1)  
!   do taup=1,natom_prim_cell
!      nc=ga+3*(taup-1) ; nq=al+3*(taup-1)  
!      atld2(al,be,ga,de) = atld2(al,be,ga,de) - 0.5*(xi(nq,al,be)-xi(nl,al,be))*zeta(nl,nc,de)  & 
!      &                                       - 0.5*(xi(nc,ga,de)-xi(nu,ga,de))*zeta(nl,nc,be)  
!   enddo
!   enddo
 enddo
 enddo
 enddo
 enddo
 qgq  =qgq  /volume_r0*ee*1d30*1d-9
!atld2=atld2/volume_r0*ee*1d30*1d-9
!call convert_to_voigt(atld2,c2)
!call write_out (uio,' Elastic Tensor (second term in GPa) in voigt ',c2)
 call convert_to_voigt(qgq,cq)
 call write_out (uio,' Elastic Tensor correction QGQ in voigt ',cq)


! atld3=0  ! last term  phi xi xi ; it leads to a symmetric 6x6 tensor ---------------------------
 do al=1,3
 do be=1,3
 do ga=1,3
 do de=1,3
!   do tau=1,natom_prim_cell
!      nl=al+3*(tau-1) ; nu=ga+3*(tau-1) 
!   do taup=1,natom_prim_cell
!      nc=ga+3*(taup-1) ; nq=al+3*(taup-1) 
!      atld3(al,be,ga,de) = atld3(al,be,ga,de) - 0.5*phi(nl,nc)*(xi(nq,al,be)-xi(nl,al,be)) &
!      &                                                       *(xi(nc,ga,de)-xi(nu,ga,de)) 
!   enddo 
!   enddo 
!! !!  is is + or - ?
    atld0(al,be,ga,de)=atld0(al,be,ga,de) - (sigma0(al,ga)*delta_k(be,de)+ sigma0(be,de)*delta_k(al,ga))/2
 enddo
 enddo
 enddo
 enddo
! atld3=atld3/volume_r0*ee*1d30*1d-9
! call convert_to_voigt(atld3,c3)
! call write_out (uio,' Xi^2 term in Elastic Tensor in voigt ',c3)
 
! should have c2+c3=-cq

! call write_out (uio,' Total Elastic Tensor (new formula in GPa) in voigt ',c1+c2+c3-cq)

 atld0=atld0-qgq
 call symmetrize4(3,atld0)

! these should be symmetric wr (al,be <-> ga,de)
! can(?) symmetrize wr (al<->be) and (ga<->de) 
!
! check symmetry
  do al=1,3
  do be=1,3
  do ga=1,3
  do de=1,3
    if (abs(atld0(al,be,ga,de)-atld0(ga,de,al,be)).gt.1d-4)  &
 &       write(uio,*)al,be,ga,de,atld0(al,be,ga,de),atld0(ga,de,al,be)
    if (abs(atld0(al,be,ga,de)-atld0(be,al,ga,de)).gt.1d-4)  &
 &       write(uio,*)al,be,ga,de,atld0(al,be,ga,de),atld0(be,al,ga,de)
    if (abs(atld0(al,be,ga,de)-atld0(al,be,de,ga)) .gt.1d-4)  &
 &       write(uio,*)al,be,ga,de,atld0(al,be,ga,de),atld0(al,be,de,ga)
!   if (abs(atld2(al,be,ga,de)-atld2(ga,de,al,be)).gt.1d-4)  &
!&       write(uio,*)al,be,ga,de,atld2(al,be,ga,de),atld2(ga,de,al,be)
!   if (abs(atld3(al,be,ga,de)-atld3(ga,de,al,be)).gt.1d-4)  &
!&       write(uio,*)al,be,ga,de,atld3(al,be,ga,de),atld3(ga,de,al,be)
!    if (abs(atld(al,be,ga,de)-atld(ga,de,al,be)).gt.1d-4)  &
! &       write(uio,*)al,be,ga,de,atld(al,be,ga,de),atld(ga,be,al,de)
!    ahat(al,ga,be,de)=0.5*(atld(al,be,ga,de)+atld1(al,de,ga,be))
  enddo
  enddo
  enddo
  enddo

  call convert_to_voigt(atld0,elastic)

! do al=1,3
! do be=1,3
! do ga=1,3
! do de=1,3
! !  ct(al,be,ga,de)=ahat(al,ga,be,de)+ahat(be,ga,al,de)-ahat(al,be,ga,de)
!    ct(al,be,ga,de)=atld1(al,ga,be,de)+atld2(al,ga,be,de)+atld3(al,ga,be,de)
! enddo
! enddo
! enddo
! enddo
! call convert_to_voigt(c1+c2+c3,elastic)
! elastic=c1+c2+c3
 call write_out (6,' Elastic Tensor(GPa) ',elastic)

 call xmatinv(6,elastic,compliance,ier)

 call write_out (6,' Compliance Tensor(1/GPa) ',compliance)

 sigmav  =sigmav  /volume_r0*ee*1d30*1d-9
 u0v=-matmul(compliance,sigmav)

 call write_out (6,' Total residual strain',u0v)

 y0(4:ndyn) =y0(4:ndyn)-matmul(gama,matmul(qiuv(4:ndyn,:),u0v))

 call write_out (6,' Total residual displacements (Ang)',y0)

 end subroutine mechanical
!============================================================
 subroutine mechanical2(atld0,elastic,ndn,xi,tet,phi,gam)
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 implicit none
 real(r15), intent(out) :: elastic(6,6)
 integer, intent(in) :: ndn 
 real(r15), intent(in) :: xi(ndn,3,3),phi(ndn,ndn),tet(ndn,ndn,3),gam(ndn-3,ndn-3),atld0(3,3,3,3)
 integer tau,taup,al,be,ga,de,i,j,t,g,ti,cnt2,ired,voigt,la,nl,nc,nu,s
          !,longt,transz
 real(r15)  rij(3),junk,atld1(3,3,3,3),atld2(3,3,3,3),atld3(3,3,3,3),c11,c12,c44,rho_SI,constr(3,3,3) &
 &       ,c1(6,6),c2(6,6),c3(6,6) 
! &       ,sigma0(3,3),y0(3*natom_prim_cell) ,pi0(3*natom_prim_cell)   
! &       ,sv(3*natom_prim_cell,3*natom_prim_cell),sw(3*natom_prim_cell),wcut &   ! used for svd
! &       ,temp(3*natom_prim_cell,3*natom_prim_cell)  &

 write(  * ,*)' ********** ENTERING MECHANICAL *************, ndn=',ndn
 write(ulog,*)' ********** ENTERING MECHANICAL *************, ndn=',ndn

! this inversion is with SVD
! temp=phi  ! phi will be overwritten
! call svdcmp(temp,ndn,ndn,sw,sv)
! wcut=maxval(abs(sw))*svdcut
! write(ulog,*)'w less than wcut is set to zero=',wcut
! do j=1,ndn
!    write(ulog,*)'after svd, before cutting, j,w(j)=',j,sw(j)
!    if(abs(sw(j)).lt.wcut) sw(j)=0d0
! enddo
! gam=0  ! here we construct gam, the inverse matrix of phi
! do nl=1,ndn
! do nc=1,ndn
!    do j=1,ndn
!       if(sw(j).ne.0d0) gam(nl,nc)=gam(nl,nc)+sv(nl,j)/sw(j)*temp(nc,j)
!    enddo
! enddo
! enddo

! test of inversion
! eye=0d0 ! construct the identity matrix
! do j=1,ndn
!    eye(j,j)=1d0   ! make the identity matrix
! enddo
! sv=matmul(gam,phi)
! junk = maxval(abs(sv-eye))
! if (junk .gt. 1d-8) then
!   write(ulog,*) ' SVD_MECHANICAL: gam is not inverse of phi ', junk
!   stop
! endif
! call write_out(ulog,' Gamma ',gam)

! pi0 is the residual force ------------------------------------
! pi0=0
! do tau=1,natom_prim_cell
! do al=1,3
!    i=al+3*(tau-1)
!    cnt2=0
!    do g=1,map(1)%ngr  ! ineq. terms and identify neighbors
!      if(g.gt.1) cnt2 = cnt2 + map(1)%ntind(g-1)
!    do t=1,map(1)%nt(g)
!       if ( map(1)%gr(g)%ixyz(1,t) .ne. al ) cycle
!       if ( map(1)%gr(g)%iat(1,t) .ne. tau ) cycle
!       do ti=1,map(1)%ntind(g)
!          ired = cnt2+ti    ! this is the corresponding index of i in ared
!          pi0(i)=pi0(i)+fcs(ired) * map(1)%gr(g)%mat(t,ti)
!       enddo
!    enddo
!    enddo
! enddo
!    write(ulog,*)'Residual force(eV/Ang): tau,pi0(tau,:)=',tau,(pi0(3*(tau-1)+al),al=1,3)
! enddo


! y0=0 ! y0 = -Gama*pi correction to equilibrium positions ----------------------
! do tau=2,natom_prim_cell
!    do al=1,3
!    i=3*(tau-1)+al
!    do j=4,ndn
!       y0(i)=y0(i)-gam(i-3,j-3)*pi0(j)
!    enddo
!    enddo
!    write(ulog,*)'Position correction(Ang): tau,u0(tau,:)=',tau,(y0(3*(tau-1)+al),al=1,3)
! enddo
!
! residual Stress tensor (under no strain) ------------------------------------
! sigma0=0
! do al=1,3
! do be=1,3
!    sigma0(al,be) = dot_product(pi0(:),xi(:,al,be))
!    do tau=1,natom_prim_cell
!!      if(al.ne.be) then
!         sigma0(al,be)=sigma0(al,be)+(atompos(al,tau)*pi0(be)+atompos(be,tau)*pi0(al))/2
!!      else
!!        sigma0(al,be)=sigma0(al,be)+atompos(al,i0)*pi0(be)
!!      endif
!    enddo
! enddo
! enddo
! call write_out(ulog,'Residual stress sigma(eta=0) (eV) ',sigma0)
! sigma0 = sigma0/volume_r0*1d30*ee*1d-9
! call write_out(ulog,'Residual stress sigma(eta=0) (GPa) ',sigma0)

! call convert_to_voigt(atld0,c1)
! call write_out (ulog,' Elastic Tensor (standard term; old formula) in GPa, in voigt ',c1)

! Standard Bravais term here: theta tau------------------------------------------
  atld1=0  
  do al=1,3
  do be=1,3
  do ga=1,3
  do de=1,3
     do tau=1,natom_prim_cell
     do taup=1,natom_prim_cell
        nl=al+3*(tau-1)
        nc=ga+3*(taup-1)
 
        atld1(al,be,ga,de) = atld1(al,be,ga,de) + atompos(be,tau)*tet(nl,nc,de) 
!
!       if(be.ne.al .and. de.ne.ga) then
!          atld1(al,be,ga,de) = atld1(al,be,ga,de) + atompos(ga,i0)*(tet(g2,n1,be)+tet(g2,n2,al))   & 
!          &                                       + atompos(de,i0)*(tet(g1,n1,be)+tet(g1,n2,al))
!       elseif(be.eq.al .and. de.ne.ga) then
!          atld1(al,be,ga,de) = atld1(al,be,ga,de) + (atompos(ga,i0)*tet(g2,n1,be)  &
!          &                                       +  atompos(de,i0)*tet(g1,n1,be))/2  
!       elseif(be.ne.al .and. de.eq.ga) then
!          atld1(al,be,ga,de) = atld1(al,be,ga,de) + (atompos(ga,i0)*tet(g2,n1,be)  &
!          &                                       +  atompos(ga,i0)*tet(g2,n2,al))/2  
!       elseif(be.eq.al .and. de.eq.ga) then
!          atld1(al,be,ga,de) = atld1(al,be,ga,de) + (atompos(de,i0)*tet(g1,n1,be))/4
!       endif
     enddo
     enddo
  enddo
  enddo
  enddo
  enddo
  atld1=atld1/volume_r0*ee*1d30*1d-9
  call convert_to_voigt(atld1,c1)
  call write_out (ulog,' Main Elastic Tensor (new formula ** WRONG **in GPa) in voigt ',c1)

 atld2=0    ! cross terms : phi*(R+tau)*Xi -----------------------------
 do al=1,3
 do be=1,3
 do ga=1,3
 do de=1,3
    do tau=1,natom_prim_cell
       nl=al+3*(tau-1)
    do taup=1,natom_prim_cell
       nc=ga+3*(taup-1)
       atld2(al,be,ga,de) = atld2(al,be,ga,de) + xi(nl,al,be)*tet(nl,nc,de) + &
       &                                         xi(nc,ga,de)*tet(nc,nl,be) 
!      phi(nl,nc)*atompos(be,tau) * xi(nc,ga,de) 
    enddo
    enddo
 enddo
 enddo
 enddo
 enddo
 atld2=atld2/volume_r0*ee*1d30*1d-9
 call convert_to_voigt(atld2,c2)
 call write_out (ulog,' Elastic Tensor (second term in GPa) in voigt ',c2)


 atld3=0  ! last term  phi xi xi ------------------------------------------
 do al=1,3
 do be=1,3
 do ga=1,3
 do de=1,3
    do tau=1,natom_prim_cell
       nl=al+3*(tau-1) !;g1=la+3*(j0-1)
    do taup=1,natom_prim_cell
       nc=ga+3*(taup-1) !;g2=mu+3*(i0-1)
!          atld2(al,ga,be,de) = atld2(al,ga,be,de) -0.5* phi(nl,nc)*(xi(nl,al,ga)-xi(g1,al,ga))*(xi(g2,be,de)-xi(nc,be,de))
          atld3(al,be,ga,de) = atld3(al,be,ga,de) + phi(nl,nc)*xi(nl,al,be)*xi(nc,ga,de)
    enddo 
    enddo 

 enddo
 enddo
 enddo
 enddo
 atld3=atld3/volume_r0*ee*1d30*1d-9
 call convert_to_voigt(atld3,c3)
 call write_out (ulog,' Xi^2 term in Elastic Tensor in voigt ',c3)
 
 call write_out (ulog,' Total Elastic Tensor (new formula in GPa) in voigt ',c1+c2+c3)

! for the second term, we invert sum_N fc2(tau,N+taupp) in gamma and use it to get X
! do i0=1,natom_prim_cell
! do al=1,3
!    nl=al+3*(i0-1) 
!   ex(nl,:,al)=ex(nl, , )+matmul(gam,tet(:,:,al))  ! this is the X matrix=d u0/d eta
!   tet(:,:,al)=matmul(gam,tet(:,:,al))  ! this is the X matrix=d u0/d eta
! enddo
! enddo

! check symmetry
! do al=1,3
! do be=1,3
! do ga=1,3
! do de=1,3
!   if (abs(atld(al,be,ga,de)-atld(ga,be,al,de)).gt.1d-4)  &
!&       write(ulog,*)al,be,ga,de,atld(al,be,ga,de),atld(ga,be,al,de)
!   if (abs(atld(al,be,ga,de)-atld(ga,de,al,be)).gt.1d-4)  &
!&       write(ulog,*)al,be,ga,de,atld(al,be,ga,de),atld(ga,be,al,de)
!   ahat(al,ga,be,de)=0.5*(atld(al,be,ga,de)+atld1(al,de,ga,be))
! enddo
! enddo
! enddo
! enddo

! do al=1,3
! do be=1,3
! do ga=1,3
! do de=1,3
! !  ct(al,be,ga,de)=ahat(al,ga,be,de)+ahat(be,ga,al,de)-ahat(al,be,ga,de)
!    ct(al,be,ga,de)=atld1(al,ga,be,de)+atld2(al,ga,be,de)+atld3(al,ga,be,de)
! enddo
! enddo
! enddo
! enddo
! call convert_to_voigt(c1+c2+c3,elastic)
 elastic=c1+c2+c3
! call write_out (6,' Elastic Tensor ',elastic)

 end subroutine mechanical2
!============================================================
 subroutine sound_speeds (elastic)
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants, only : r15
 implicit none
 real(r15), intent(out) :: elastic(6,6)
 integer i0,j0,al,be,ier,j,t,g,ti,cnt2,ired,mu,nu,ndn
 real(r15)  rj(3),junk,khat_red(3)
 real(r15)     eivl(3*natom_prim_cell)  ,c11,c12,c44
 complex(r15)  eivc(3*natom_prim_cell,3*natom_prim_cell)
 complex(r15) dynr2(3*natom_prim_cell,3*natom_prim_cell)

 ndn=3*natom_prim_cell
 write(ulog,*)' ********** ENTERING SOUND SPEEDS *************'
 khat_red=(g01+g02)/length(g01+g02)
 dynr2=0
! calculate the matrix sum Phi (R.khat)^2/2 for k along (110)
 do al=1,3
 do be=1,3
 do i0=1,natom_prim_cell
 cnt2=0  ! cumulative number of independent terms up to group (g-1)
 gloop: do g=1,map(2)%ngr
    if(keep_grp2(g).ne.1) cycle
    if(g.gt.1) then
       if (keep_grp2(g-1).eq.1) cnt2=cnt2+map(2)%ntind(g-1)
    endif
    do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
       ired=cnt2+ti + map(1)%ntotind  ! index of the indep FC2 in the list fcs
       if(cnt2+ti .gt. size_kept_fc2) then
          write(*,*)'SET_DYNMAT: size of fc2 exceeded, ired=',ired
          stop
       endif
       tloop: do t=1,map(2)%nt(g)
          if ( i0 .ne. map(2)%gr(g)%iat (1,t) .or. &
&              al .ne. map(2)%gr(g)%ixyz(1,t) .or. &
&              be .ne. map(2)%gr(g)%ixyz(2,t) ) cycle tloop

          j = map(2)%gr(g)%iat(2,t)
          j0 = iatomcell0(j)
          rj = atompos(:,j) - atompos(:,j0)
! need its reduced coordinates
          rj=matmul(cart_to_prim,rj) !/(2*pi)
          mu=3*(i0-1)+al
          nu=3*(j0-1)+be
          junk = fcs(ired)* map(2)%gr(g)%mat(t,ti)  ! that is phi(al,ga)
          dynr2(mu,nu)=dynr2(mu,nu)+ fcs(ired) *(dot_product(rj,khat_red))/2d0 ! rj(1)+rj(2))**2/4
       enddo tloop
    enddo
 enddo gloop
 enddo
 enddo
 enddo

! diagonalize dynr2 to get speeds of sound
    call diagonalize(ndn,dynr2,eivl,ndn,eivc,ier)

! assumes eivl is in eV
    call get_cij(eivl,eivc,ndn,c11,c12,c44)
    call write_out(ulog,'From sound_speeds:c11,c12,c44(GPa)',(/c11,c12,c44/))
    elastic=0d0
    elastic(1,1)=c11
    elastic(2,2)=c11
    elastic(3,3)=c11
    elastic(1,2)=c12
    elastic(1,3)=c12
    elastic(3,2)=c12
    elastic(2,1)=c12
    elastic(3,1)=c12
    elastic(2,3)=c12
    elastic(4,4)=c44
    elastic(5,5)=c44
    elastic(6,6)=c44

 end subroutine sound_speeds
!============================================================
 subroutine thermo(nk,wk,ndyn,eival,grn,tempk,etot,free,pres,cv)
! calculate total and free energies within QHA, at a given temperature (temp in Kelvin)
 use ios
 use params
 use lattice
 use constants
 implicit none
 integer nk,ndyn,b,k,nat
 real(r15) wk(nk),eival(ndyn,nk)
 real(r15) x,cv_nk,cv,hw,free,etot,pres,nbe,mdedv,pres0,nbx
 complex(r15) grn(ndyn,nk)

    nat = ndyn/3
    if (tempk.le.0) then
       write(ulog,*)'temperature not in the proper range!!',tempk
       stop
    endif
    etot=0 ; cv=0 ; free=0 ; pres=0
    mdedv= 0.35388/(20.8**3-20.0**3)/ab**3 ! this is -dE/dV in eV/A^3
    mdedv = mdedv*1d+30*ee ! (in J/m^3 )
    do k=1,nk
    do b=1,ndyn
       if(eival(b,k) .lt.0) then
          x=0
          write(ulog,*) 'THERMO: negative eival for band,kp# ',b,k,eival(b,k)
          write(ulog,*) ' will use its absolute value instead!'
!      else
       endif
! to convert eV/A^2/mass to Hz we need to take the sqrt, multiply by cnst*100*c_light
       x=(h_plank*sqrt(abs(eival(b,k)))*cnst*100*c_light)/k_b/tempk
       if (x.gt.60) then
           cv_nk = 0
       else
           nbx=nbe(x,1d0,classical)
           cv_nk = x*x*nbx*(1+nbx) !/4/sinh(x/2)/sinh(x/2)
       endif
       hw = x*k_b*temp  ! hbar*omega in Joules
       cv  = cv + cv_nk*wk(k)   ! this is in J/K units : per unit cell
       etot= etot + hw * (nbe(x,1d0,classical) + 0.5d0) * wk(k)
       free= free + (0.5d0*hw + k_b*temp*log(1-exp(-x))) * wk(k)
       pres= pres + hw * (nbe(x,1d0,classical) + 0.5d0) * wk(k) * real(grn(b,k))
    enddo
    enddo
    cv = cv/nat*n_avog*k_b  ! this is per mole    )/(volume_r0*1d-30)  !
    etot = etot/nat*n_avog   ! convert from Joule/cell to Joule per mole
    free = free/nat*n_avog   ! convert from Joule/cell to Joule per mole
    pres0= pres/(volume_r0*1d-30)  * 1d-9  ! 1d-9 to convert to GPa; 1d-8 is to convert to kbar
    pres = pres0 !+ mdedv
3 format(9(2x,g11.4))

 end subroutine thermo
!=======================================================
 subroutine check_mdyn(ndyn,nk,kp,eivec,eival)
! form the dynamical matrix from existing eigenvectors and eigenvalues
! and diagonalize and verify consistency :eivec(-q)=-conjg(eivec(q))
 use lattice  ! needed to access NC,g1,g2,g3
 use kpoints
 use constants, only : r15
 implicit none
 integer ndyn,nk,j,k,l,i,ier,i1,j1,k1,mi,inside
 real(r15) kp(3,nk),eival(ndyn,nk),eivl(ndyn)
 complex(r15) eivec(ndyn,ndyn,nk),dynm(ndyn,ndyn),eivc(ndyn,ndyn),d2(ndyn,ndyn)

 do i=1,nk
    dynm=0
    do j=1,ndyn
    do k=1,ndyn
      do l=1,ndyn
        dynm(j,k)=dynm(j,k)+eival(l,i)*eivec(j,l,i)*conjg(eivec(k,l,i))
      enddo
    enddo
    enddo
    call get_k_info(-kp(:,i),NC,mi,i1,j1,k1,inside)
    d2=0
    do j=1,ndyn
    do k=1,ndyn
      do l=1,ndyn
        d2(j,k)=d2(j,k)+eival(l,mi)*eivec(j,l,mi)*conjg(eivec(k,l,mi))
      enddo
    enddo
    enddo

    do j=1,ndyn
    do k=1,ndyn
       if(abs( d2(j,k)-conjg(dynm(j,k))) .gt.1d-5) then
         write(*,4)'CHECK_MDYN: j,k,d(-q),d(q)=',j,k,d2(j,k),dynm(j,k)
       endif
    enddo
    enddo
!   return

    call diagonalize(ndyn,dynm,eivl,ndyn,eivc,ier)
    do j=1,ndyn
       if(abs(eivl(j)-eival(j,i)).gt.1d-4) then
          write(*,3)'CHECK_MDYN:j,eivl,eival=',j,eivl(j),eival(j,i)
       endif
       do k=1,ndyn
         if(abs(eivc(k,j)-eivec(k,j,i)).gt.1d-4) then
         if(abs(eivc(k,j)+eivec(k,j,i)).gt.1d-4) then
            write(*,4)'CHECK_MDYN:j,k,eiv(j),eivecs(k,j)=',j,k,eivl(j),eivc(k,j),eivec(k,j,i)
         endif
         endif
       enddo
    enddo
 enddo
3 format(a,i6,9(1x,f15.6))
4 format(a,2i6,f15.7,9(1x,f11.4))

 end subroutine check_mdyn
!===========================================================
 subroutine write_eivecs(ndyn,nk,kp,eigenvec)
 use lattice
 use ios
 use geometry
 use constants
 use kpoints ! for shift
 implicit none
 integer nk,ndyn,i,la,j,i2,j2,k2,inside,l,i3,j3,k3
 real(r15) kp(3,nk),w(3)
 complex(r15) eigenvec(ndyn,ndyn,nk)

 write(ulog,*)'WRITE_EIGVECS: eivc(q)-conjg(eivc(-q)),eivc(q),eivc(-q)'
do i=1,nk
! bring kp(j) into the FBZ

   call get_k_info(kp(:,i)-shift,NC,l,i3,j3,k3,inside)
   if (i.ne.l) then
      write(ulog,*) 'i .ne. l ',i,l
      stop
   endif
!  call bring_to_cell_c(-kp(:,i),g1,g2,g3,r1,r2,r3,w)
!  w =w /2/pi
   w = -kp(:,i)+shift
   call get_k_info(w,NC,j,i2,j2,k2,inside)
!  write(ulog,5)'k,-k=',i3,j3,k3,' ** ',i2,j2,k2,kp(:,i),w

  if (j.lt.nk) then
   do la=1,ndyn
   do l=1,ndyn
      if(abs(eigenvec(l,la,i)-conjg(eigenvec(l,la,j))).gt.1d-4) then
      if(abs(eigenvec(l,la,i)+conjg(eigenvec(l,la,j))).gt.1d-4) then
          write(ulog,8)'la,l,k,-k=', la,l,i,j,eigenvec(l,la,i)-conjg(eigenvec(l,la,j)) &
      &     ,eigenvec(l,la,i),eigenvec(l,la,j)
      endif
      endif
   enddo
   enddo
  else
    write(ulog,4)'kpoint i in FBZ=',i,kp(:,i)
    write(ulog,4)'was taken to -k=',j,w
  endif
enddo
4 format(a,i5,2x,99(1x,f8.3))
5 format(a,3i5,a,3i5,2x,99(1x,f8.3))
8 format(a,4i6,99(1x,f8.3))

 end subroutine write_eivecs
!==========================================================
 subroutine enforce_asr_phi
! enforces asr on force constants instead of the dynamical matrix
 use params
 use lattice
 use kpoints
 use atoms_force_constants
 use geometry
 use svd_stuff
 use constants, only : r15
 implicit none
 real(r15) mi,mj,nrm1,rr(3),asr,junk,fcold
 integer i0,j,j0,al,be,i3,i4,j3,t,cnt2,ti,ired,g,iredsave

4 format(a,4(1x,i5),9(2x,f9.4))
9 format(i6,99(1x,f9.3))

 do i0=1,natom_prim_cell
!   write(*,3)'i0,ri0=',i0,atom0(i0)%equilibrium_pos ,atompos(:,i0)
 do al=1,3
 do be=1,3
    i3 = al+3*(i0-1)
    i4 = be+3*(i0-1) ! used for asr
    mi = atom0(i0)%mass
! write(ulog,*) 'i,al,mass=',i0,al,i3,mi
    asr = 0 
    cnt2=0  ! counter of ntotind, to find position of each term in the array fc2
    gloop: do g=1,map(2)%ngr

       if(keep_grp2(g).ne.1) cycle gloop
       if(g.gt.1) then
         if (keep_grp2(g-1).eq.1) cnt2=cnt2+map(2)%ntind(g-1)
       endif
       tiloop: do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
          ired=cnt2+ti + map(1)%ntotind  ! the latter is the size of fc1
          if(cnt2+ti .gt. size_kept_fc2) then
             write(*,*)'SET_DYNMAT: size of fc2 exceeded, ired=',ired
             stop
          endif

          tloop: do t=1,map(2)%nt(g)
             if ( i0 .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t)  &
&                .or. be .ne. map(2)%gr(g)%ixyz(2,t) ) cycle tloop
             j  = map(2)%gr(g)%iat(2,t)
             j0 = iatomcell0(j)
! this produces the diagonal term fc(i0,al;i0,be)
             if ( i0 .eq. map(2)%gr(g)%iat(2,t) ) iredsave=ired
             mj = atom0(j0)%mass
             j3 = be+3*(j0-1)
             rr = atompos(:,j) - atompos(:,j0)
             junk = fcs(ired) * map(2)%gr(g)%mat(t,ti)
             asr=asr+junk
          enddo tloop
       enddo tiloop
    enddo gloop
    fcold=fcs(iredsave)
    fcs(iredsave)=fcs(iredsave)-asr
    if(abs(asr).gt.1d-6) then
       write(*   ,5)'i0,al,be,asr_corr,oldfc,newfc=',i0,al,be,asr,fcold,fcs(iredsave)
       write(ulog,5)'i0,al,be,asr_corr,oldfc,newfc=',i0,al,be,asr,fcold,fcs(iredsave)
    endif
 enddo
 enddo
 enddo
5 format(a,3i4,9(1x,g13.6))
 end subroutine enforce_asr_phi
!===========================================================
 subroutine enforce_asr_simple(n,dyn)
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 integer nat,al,be,io,jo,lo,iat,jat
 complex(r15), intent(inout) :: dyn(n,n)
 complex(r15) sumd

 nat=n/3
 do iat=1,nat
 do al=1,3
    io=al+3*(iat-1)
    do be=1,3
       lo=be+3*(iat-1)
       sumd=0
       do jat=1,nat
          if (iat.eq.jat) cycle
          jo=be+3*(jat-1)
          sumd=sumd-dyn(io,jo)
       enddo
       write(*,4)'enforce_asr: old, new d(i,i)=',io,lo,dyn(io,lo),sumd
       dyn(io,lo)=sumd
    enddo
 enddo
 enddo
4 format(a,2i4,9(1x,g11.4))
 end subroutine enforce_asr_simple
!===========================================================
 subroutine check_asr(n,dyn)
 use constants, only : r15
 use ios, only : ulog
 implicit none
 integer, intent(in) :: n
 complex(r15), intent(inout) :: dyn(n,n)
 integer nat,al,be,io,jo,iat,jat
 complex(r15) sumd

 nat=n/3
 do iat=1,nat
 do al=1,3
    io=al+3*(iat-1)
    do be=1,3
!      write(*,*)'i,al,be=',iat,al,be
       sumd=cmplx(0d0,0d0)
       do jat=1,nat
          jo=be+3*(jat-1)
          sumd=sumd+dyn(io,jo)
       enddo
       if(abs(sumd).gt.9d-4) then
          write(*,5)'ASR_CHECK: i0,al,be,sum=',iat,al,be,sumd
          write(ulog,5)'ASR_CHECK: i0,al,be,sum=',iat,al,be,sumd
          call enforce_asr_simple(n,dyn)
!      else
!         write(*,5)'ASR_CHECK PASSED: i0,ial,be,sum=',iat,al,be,sumd
       endif
    enddo
 enddo
 enddo

5 format(a,3i3,2(1x,f10.6))

 end subroutine check_asr
!===========================================================
subroutine band_sort_bs(nkp,ndyn,kp,eivec,emap)
 use constants, only : r15
implicit none
integer, intent(in) :: nkp,ndyn
real(r15), intent(in) :: kp(3,nkp)
!real(r15), intent(in) :: eival(ndyn,nkp),dk(nkp)
complex(r15), intent(in) :: eivec(ndyn,ndyn,nkp)
integer, intent(out) :: emap(ndyn,nkp)

integer i,j,k,l
integer ntbd !Number of connections to be determined
complex(r15), allocatable :: overlap(:,:)
real(r15), allocatable :: overlap_q(:,:)

!First estimate the maximum band derivative bmax=max(dE/dk)
!bmax=0
!fsaf=10
!allocate(dk(2:nkp))
!do i=2,nkp
!   dk(i)=sqrt((kp(1,i)-kp(1,i-1))**2+(kp(2,i)-kp(2,i-1))**2+(kp(3,i)-kp(3,i-1))**2)
!   do j=1,ndyn
!      bmax = max(bmax,abs(eival(j,i)-eival(j,i-1))/dk(i))
!   end do
!end do

!Start from the first k-point, which will be used as the reference for other k points
do i=1,ndyn
   emap(i,1)=i
end do

!from the second k-point to the last, calculate the overlap matrix and do band sorting
allocate (overlap(ndyn,ndyn),overlap_q(ndyn,ndyn))

do i=2,nkp
    ntbd=ndyn
    overlap=matmul(transpose(dconjg(eivec(:,:,i-1))),eivec(:,:,i))
    do j=1,ndyn
    do k=1,ndyn
        overlap_q(j,k)=dble ( overlap(j,k)*dconjg(overlap(j,k)) )
    end do
    end do

!Step 0: if two bands with energy difference larger than fsaf*bmax, the bands are not supposed to be connected
!    write(*,*) "energy window for band sorting:", fsaf*bmax*dk(i)
!    do j=1,ndyn-1
!    do k=j+1,ndyn
!       if (abs(eival(j,i-1)-eival(k,i)) .ge. fsaf*bmax*dk(i)) then
!          overlap_q(j,k)=0
!       end if
!    end do
!    end do

!Step 1: find matrix elements larger than 0.6, which indicates a strong connectivity
    do j=1,ndyn
    do k=1,ndyn
        if (overlap_q(j,k) .ge. 0.6) then
            do l=1,ndyn
                if ( emap(l,i-1) .eq. j) then
                    emap(l,i)=k
                end if
            end do
                ntbd=ntbd-1
                overlap_q(j,:)=0  !Connected ones will not be examined again
                overlap_q(:,k)=0
        end if
    end do
    end do

    if (ntbd .ne. 0) then !If there are unresolved connections remaining
!Step 2: find the largest remaining matrix elements in each row, and connect eigenvalues connected by that.
    do j=1,ndyn
        if (maxval(overlap_q(j,:)) .ge. 0.1) then
kkloop:     do k=1,ndyn
               if (overlap_q(j,k) .eq. maxval(overlap_q(j,:))) then
                  if (overlap_q(j,k) .lt. 0.3) then
                     write(*,*) "Warning: the overlap matrix element is &
   & smaller than 0,3, there may be ambiguity for band connectivity"

 write(*,*) "k-point:",kp(:,i),"index of kp:",i,"Matrix element:",overlap_q(j,k)
                  end if
                  do l=1,ndyn
                     if ( emap(l,i-1) .eq. j) then
                        emap(l,i)=k
                     end if
                  end do
                  ntbd=ntbd-1
                  overlap_q(j,:)=0
                  overlap_q(:,k)=0
                  exit kkloop
               end if
            end do kkloop
        end if
    end do
    end if

    if (ntbd .ne. 0) then !If there are still unresolved connections remaining
!Step 3: connect the remaining pairs in a ascending order
    do j=1,ndyn
    do k=1,ndyn
        if (overlap_q(j,k) .ne. 0) then
write(*,*) "Warning: the overlap matrix element is smaller than 0.3, &
& there may be ambiguity to determine the band connectivity"
write(*,*) "k-point:",kp(:,i),"index of kp:",i,"Matrix element:",overlap_q(j,k)
           do l=1,ndyn
              if ( emap(l,i-1) .eq. j) then
                 emap(l,i)=k
              end if
           end do
        ntbd=ntbd-1
        overlap_q(j,:)=0
        overlap_q(:,k)=0
        end if
    end do
    end do
    end if

 !      write(*,*) 'map for',i, emap(:,i)

    if (ntbd .ne. 0) then
        write(*,*) 'error: band sorting failed'
        write(*,*) 'k-point:',i,kp(:,i),ntbd
    end if
end do

deallocate(overlap,overlap_q)

end subroutine band_sort_bs
!==========================================================
 subroutine diagonalize(n,mat,eival,nv,eivec,ier)
! n=size of mat; nv is the number of needed eigenvectors
 use constants, only : r15
 implicit none
 integer, intent(in) :: n,nv
 integer, intent(out) :: ier
 complex(r15), intent(in) :: mat(n,n)
 complex(r15), intent(out) :: eivec(n,n)
 real(r15), intent(out) :: eival(n)
! This is used by eigch
! real(r15), allocatable :: w(:,:)
! integer, allocatable :: lw(:)
! This is used by ZHEGV
  real(r15), allocatable :: rwork(:)
  complex(r15), allocatable :: work(:)
  integer lwork  !, zero


! n = size(mat(:,1))
 if (n .ne. size(eival) ) then
    write(   *,*)' EIGCH, size inconsistency:mat,eival=',n, size(eival)
    stop
 endif

! tol = -1.d0
! zero= 0  ! for now do not compute any eigenvectors
! ier=0
! allocate(w(n,7),lw(n))
! call eigch(mat,n,n,n,nv,tol,w,lw,eival,eivec,ier)
! deallocate(w,lw)
! return

  lwork = max(1,2*n-1)
  allocate(work(lwork),rwork(max(1,3*n-2)))
  if(nv.eq.0) then
    call zheev('N','U',n,mat,n,eival,work,lwork,rwork,ier)
  else
    call zheev('V','U',n,mat,n,eival,work,lwork,rwork,ier)
    eivec = mat
  endif
  deallocate(work,rwork)
4 format(i5,1x,3(f6.3,1x,f6.3,4x))

 end subroutine diagonalize
!==========================================================
 subroutine write_all_eigenvalues(nk,dk,kp,eival,vg,n,uio)
 use constants
 use eigen, only : mysqrt
 use geometry
 implicit none
 integer, intent(in) :: n,uio,nk
 real(r15), intent(in) :: dk(nk),kp(3,nk),eival(n,nk),vg(3,n,nk)
 integer j,k

 do k=1,nk
 do j=1,n
    write(uio,3)k,j,dk(k),kp(:,k),cnst*mysqrt(eival(j,k)),vg(:,j,k),length(vg(:,j,k))
 enddo
 enddo

2 format(i5,1x,4(1x,f8.3),999(1x,g11.4))
3 format(2i5,1x,4(1x,f8.3),999(1x,g11.4))

 end subroutine write_all_eigenvalues
!==========================================================
 subroutine calculate_dos(mx,eival,wkp,mesh,omega,ds)
! use gaussian broadening
 use ios
 use params
 use om_dos
 use constants
 implicit none
 integer i,j,mx,mesh
 real(r15) x,wkp(mx),delta,ds(mesh),omega(mesh),eival(mx)  !,cnst

! wkp=1d0/(nkx*nky*nkz)
! write(udos,*)'# wkp,width =',wkp,width
    width=(omega(2)-omega(1))/3

    do i=1,mesh
       ds(i) = 0
       do j=1,mx
          x = (eival(j) - omega(i))/width   ! these are all in cm^-1
!         x = (eival(j) - ene*ene)/width/width
          if ( abs(x) .gt. 5 ) cycle
          ds(i) = ds(i) + delta(x)/width*wkp(j) !/mx
       enddo
    enddo

3 format(i5,9(3x,g11.4))

 end subroutine calculate_dos
!==========================================================
 subroutine calculate_jdos(mx,eival,wkp,mesh,omega,udosj)
! use gaussian broadening to calculate the joint dos
 use ios
 use params
 use om_dos
 use constants
 implicit none
 integer i,j,k,mx,mesh,udosj
 real(r15) xp,xm,wkp(mx),delta_g,dsp,dsm,omega(mesh),eival(mx),iself,tmp,one
 complex(r15) omz

! wkp=1d0/(nkx*nky*nkz)
 one = 1d0
 write(udosj,*)'# width =',width
 tmp =0.07      ! the temp width =0.07/cm is equiv to 0.1 K
! then (1+n1+n2)=1 should produce dsp results

! open(udosj+1,file='jdos.discrete',status='unknown')
! do j=1,mx
! do k=j,mx
!    xp= cnst*(sqrt(abs(eival(j)))+sqrt(abs(eival(k))))
!    xm= cnst*(sqrt(abs(eival(j)))-sqrt(abs(eival(k))))
!    write(udosj+1,2) xp,one
!!    if (xm.gt.0)
!    write(udosj+1,2) abs(xm),-one
! enddo
! enddo
! close(udosj+1)

    write(udosj,*)'# i,omega(i),jdsm(i),jdsp(i), occupation-weighted jdos'
    do i=1,mesh
       dsp=0 ; dsm=0
       omz=2d0*dcmplx(omega(i),-etaz)
       do j=1,mx
       do k=1,mx
          xp= (cnst*(sqrt(abs(eival(j)))+sqrt(abs(eival(k)))) - 2*omega(i))   ! these are all in cm^-1
          xm= (cnst*(sqrt(abs(eival(j)))-sqrt(abs(eival(k)))) - 2*omega(i))   ! these are all in cm^-1
          if ( abs(xp) .lt. 6 ) then
             dsp = dsp + delta_g(xp,width)*wkp(j)*wkp(k) !* width*sqrt(2*pi)
          endif
          if ( abs(xm) .lt. 6 ) then
             dsm = dsm + delta_g(xm,width)*wkp(j)*wkp(k) !* width*sqrt(2*pi)
          endif

!          iself = aimag(oc2(tmp,omz,j,k))/pi*wkp(j)*wkp(k)*tmp
          iself = (dsp+dsm)/pi
       enddo
       enddo
       write(udosj,3) i,2*omega(i),dsm,dsp,iself
    enddo
    close(udosj)

2 format(9(3x,g11.4))
3 format(i5,9(3x,g11.4))

 end subroutine calculate_jdos
!==========================================================
 subroutine lifetime_dos(q,omq,dosp,dosm)
! use gaussian broadening to calculate the joint dos
! for given q and mode index la, calculates sum_1,2
! delta(1-q1-/+q2)delta(w_q-w1-/+w2)
 use ios
 use params
 use om_dos
 use constants
 use eigen
 use kpoints
 implicit none
 integer l,i1,j1,k1,nq1,ip,jp,kp,im,jm,km,nqp,nqm,la1,la2,inside
 real(r15) q(3),q1(3),qp(3),qm(3),dosp,dosm,om1,omp,omm,omq,delta_l

 dosp=0; dosm=0
 do l=1,nkc
    q1=kpc(:,l)
    call get_k_info_cent(q1-shift,NC,nq1,i1,j1,k1,inside)
    if(l.ne.nq1) then
      write(ulog,*)'n1,nq1,inside=',l,nq1,inside
      write(ulog,*)'q1=',q1
      write(ulog,*)'ijk1=',i1,j1,k1
      stop
    endif
    qp=-q-q1
    qm=-q+q1
    call get_k_info_cent(qp-shift,NC,nqp,ip,jp,kp,inside)
    call get_k_info_cent(qm-shift,NC,nqm,im,jm,km,inside)
    do la1=1,ndyn
       om1 = eigenval(la1,nq1)
       do la2=1,ndyn
          omp = eigenval(la2,nqp)
          omm = eigenval(la2,nqm)
          dosp=dosp+delta_l(omq-om1-omp,etaz)
          dosm=dosm+delta_l(omq-om1+omm,etaz)
       enddo
    enddo
 enddo
 end subroutine lifetime_dos
!===========================================================
 subroutine phase_space_lifetime(i1,i2,nk,eival)
 use ios
 use om_dos
 use kpoints
 use eigen
 use constants
 implicit none
 integer i,la,nk,i1,i2,udos
 real(r15) q(3),omq,dosp,dosm,eival(ndyn,nk)
 character cik1*4,cik2*4


 write(cik1,'(i4.4)') i1
 write(cik2,'(i4.4)') i2
 udos=459
 open(udos,file='lifetime_dos-'//cik1//'-'//cik2//'.dat ',status='unknown')
 write(udos,*)'# la,q,omega(q,la),jdsm(i),jdsp(i)'

    do i=i1,i2
       q=kibz(:,i)
!       call get_freq(q,ndyn,vq,eivl,eivc)
       do la=1,ndyn
          omq = eival(la,i)
          call lifetime_dos(q,omq,dosp,dosm)
          write(udos,3) la,q,omq,dosm/nkc,dosp/nkc
       enddo
    enddo
 close(udos)

3 format(i5,9(3x,g11.4))

 end subroutine phase_space_lifetime
!===========================================================
  subroutine cubic_on_gconv(coord)
! direct coords of conv cubic on present conv
  use lattice
  use ios
  use geometry
 use constants, only : r15
 implicit none
 !real(r15), intent(out) :: coord(3,3)
 real(r15) coord(3,3)

  call write_out(6,'coord before',coord)
 coord(1,1)=(g02+g03-g01).dot.r1conv
 coord(2,1)=(g02+g03-g01).dot.r2conv
 coord(3,1)=(g02+g03-g01).dot.r3conv
 coord(1,2)=(g01+g03-g02).dot.r1conv
 coord(2,2)=(g01+g03-g02).dot.r2conv
 coord(3,2)=(g01+g03-g02).dot.r3conv
 coord(1,3)=(g02+g01-g03).dot.r1conv
 coord(2,3)=(g02+g01-g03).dot.r2conv
 coord(3,3)=(g02+g01-g03).dot.r3conv
  call write_out(6,'coord after ',coord)

  end subroutine cubic_on_gconv
!============================================================
 subroutine get_cij(eivl,eivc,ndn,c11,c12,c44)
!! reads eigenvalues of L_ik = C_ij,kl qhat_j qhat_l in eV along 110 and
!! outputs the elastic constants of a CUBIC crystal
 use atoms_force_constants , only : natom_prim_cell
 use lattice, only : volume_r0
 use constants, only : ee,r15  !c_light,uma,
 use ios !, only : ulog
 implicit none
 integer, intent(in) :: ndn
 real(r15), intent(inout) :: eivl(ndn)
 complex(r15), intent(in) :: eivc(ndn,ndn)
 real(r15), intent(out) :: c11,c12,c44
 real(r15) pol,polt,junk,junt,q(3),halfdif,halfsum
 integer la,i,transz,longt

 ! if polarization is out of plane, i.e. along z, we get c44
 ! which polarization is along z
     call write_out(ulog,'mass_cell w^2 /q^2 ',eivl)
     q=(/1,1,0/)/sqrt(2d0)  ! normalized wave vector
     pol=-3; polt=-3 ! longitudinal and transverse polarizations of the mode
     do la=1,3  ! look at the 3 transverse modes la
! dentify transverse shear and in-plane modes by calculting the
! dot product of eigenvector polarizations summed over all atoms with q
        junk=0; junt=0;
        do i=1,natom_prim_cell
           junk=junk+abs(dot_product(q,eivc(3*(i-1)+1:3*(i-1)+3 ,la))) ! long q which is in xy plane
           junt=junt+abs(eivc(3*(i-1)+3,la))  ! along z
           write(*,7)'mode,iatom,projected polalong 110 and 001=',la,i,junk,junt
        enddo
! find the modes with in-plane polarizations (junk) and those along z(junt)
        if(junt.gt.polt) then  ! save the largest projections perp to q
          polt=junt
          transz=la
        endif
        if(junk.gt.pol) then  ! save the largest projections || to q
          pol=junk
          longt=la
        endif
     enddo
     write(*,*)'transverse mode in the z direction is',transz
     write(*,*)'transverse mode in the 1-10 direction is',6-longt-transz
     write(*,*)'longitudinal mode in the 110 direction is',longt
 ! rescale eivl from eV to GPA
     eivl = eivl * ee*1d21/volume_r0
     c44 = eivl(transz)
     halfsum = eivl(longt)-c44
     halfdif = eivl(6-longt-transz)
     c11=halfsum+halfdif
     c12=halfsum-halfdif
     call write_out(6   ,'c11,c12,c44(GPa)',(/c11,c12,c44/))
7 format(a,2i6,2(1x,g11.4))

 end subroutine get_cij


