 subroutine set_band_Structure !(map) !fcs,ngr) !,dyn)
!! reads irreducible 2nd order force constants, and computes the phonon spectrum
 use svd_stuff, only : map
 use kpoints
 use atoms_force_constants, only : natom_prim_cell
 use ios, only : ulog , uibs, uband, ugrun
 use eigen !, only : ndyn
!  integer, intent(in) :: ngr
!  real(8), intent(in) :: fcs(ngr)
!  complex(8), allocatable, intent(inout) :: dyn(:,:)
 real(8) elastic(6,6) , coord(3,3)

  ndyn = 3*natom_prim_cell

  write(*,*)'Entered set_band_structure '
  write(ulog,*)'map(ntotind)=',map(:)%ntotind
!---------------------------------------------------------------
! do a band structure calculation along the symmetry directions

! call cubic_on_gconv(coord)
  call write_out(6,'Direct coords of G_cubic on G_conv',coord)

  call make_kp_bs

! write(ulog,*)' KEXT_bs on gconv is '
! write(*,*)' KEXT_bs on gconv is '
! do i=1,ndir+1
!    q=0
!    do j=1,3
!!      q(j)=dot_product(coord(:,j),kext_bs(:,i))
!    do l=1,3
!       q(j)=q(j)+coord(l,j)*kext_bs(l,i)
!    enddo
!    enddo
!    q=matmul(transpose(coord),kext_bs(:,i))
!    write(ulog,4)i,q
!    write(*,4)i,q
! enddo
  4 format(i4,9(1x,f9.5))

  write(ulog,*)' Kpoints for band structure generated from kpbs.params'
  write(*,*)' Kpoints for band structure generated from kpbs.params'
  call allocate_eig_bs(ndyn,nkp_bs,ndyn) !3rd index is for eivecs needed for gruneisen

! now output eigenval is directly frequencies in 1/cm
  write(*,*)' calling get_freqs'
  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,ndyn,eigenvec_bs,veloc_bs)

  open(uband,file='bs_freq.dat')
  call write_all_eigenvalues(nkp_bs,dk_bs,kp_bs,eigenval_bs,veloc_bs,ndyn,uband)
  close(uband)

  open(ugrun,file='bs_grun.dat')
  call gruneisen(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,eigenvec_bs,ugrun,grun_bs)
  close(ugrun)

  call deallocate_eig_bs
  call deallocate_kp_bs

  call mechanical(elastic)
  call write_out(ulog,' Elastic tensor from FCs ',elastic)
! call sound_speeds(elastic)
  call write_out(ulog,' Elastic tensor from Vgr ',elastic)

 end subroutine set_band_Structure
!===========================================================
 subroutine get_frequencies(nkp,kp,dk,ndn,eival,nv,eivec,vg)
! also outputs the group velocities from HF theorem (exact) in units of c_light
 use params
 use constants
 use born
 use geometry
 use ios, only : uband
 implicit none
 integer, intent(in) :: nkp,ndn,nv ! no of wanted eivecs
 real(8), intent(in) :: kp(3,nkp),dk(nkp)
 real(8), intent(out):: eival(ndn,nkp),vg(3,ndn,nkp)
 complex(8), intent(out) :: eivec(ndn,ndn,nkp)
 integer i,j,nd2
 integer, allocatable :: mp(:) !Map for ascending band sorting
 real(8), allocatable :: eivl(:)
 complex(8), allocatable:: eivc(:,:)
 integer, allocatable :: mp1(:,:) !Map for projection band sorting
 real(8), allocatable :: eival_tmp(:) !Temp matrices for reordering eigenvalues
 complex(8), allocatable :: eivec_tmp(:,:)
 real(8), allocatable :: vg_tmp(:,:)

! open files and write the group velocities
  open(330,file='veltest.dat')
! do l=1,ndn
!    write(ext,'(i3.3)')l
!    open(uvel+l,file='veloc-'//ext//'.dat')
!    write(uvel+l,*)'# la,i1,i2,i3,nk,kp(nk),eival(l,mapibz(nk)),vg(l,nk),length(vg(l,nk))'
! enddo
 nd2 = min(ndn,12)
!allocate(dynmat(ndn,ndn),mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3))
 allocate(mp(ndn),eivl(ndn),eivc(ndn,ndn))
! write(uio,*)'# dynmat allocated, printing: Kpoint, eival(j),j=1,ndyn)'

 kloop: do i=1,nkp

! write(uio,*)'############ before setting up dynmat, kp number ',i
    call finitedif_vel(kp(:,i),ndn,vg(:,:,i),eival(:,i),eivec(:,1:nv,i)) ! don't need group velocities
!   call get_freq(kp(:,i),ndn,vg(:,:,i),eival(:,i),eivec(:,1:nv,i) )

 enddo kloop


! Band structure projection sorting, added by Bolin
   allocate(mp1(ndn,nkp))
   write(*,*) "Projection Band Sorting for band structures!"
   call band_sort_bs(nkp,ndn,kp,eivec,mp1)
   allocate(eival_tmp(ndn),eivec_tmp(ndn,ndn),vg_tmp(3,ndn))
   open(uband,file='bands.dat')
   do i=1,nkp
      eival_tmp=eival(:,i)
      eivec_tmp=eivec(:,:,i)
      vg_tmp=vg(:,:,i)
      do j=1,ndn
         eival(j,i)=eival_tmp(mp1(j,i))
         eivec(:,j,i)=eivec_tmp(:,mp1(j,i))
         vg(:,j,i)=vg_tmp(:,mp1(j,i))
      end do
      call write_eigenvalues(i,dk(i),kp(:,i),eival(:,i),vg(:,:,i),ndn,uband)
   end do
   close(uband)
   deallocate(mp1,eival_tmp,eivec_tmp,vg_tmp)
   deallocate(eivl,eivc,mp)  !,dynmat,ddyn)

 3 format(a,i5,9(1x,f19.9))
 4 format(99(1x,2(1x,g9.3),1x))
 5 format(a,i5,99(1x,2(1x,g9.3),1x))
 6 format(2i6,2x,99(1x,f9.3))
 7 format(i6,2x,3(1x,f7.3),1x,99(1x,f7.3))

 end subroutine get_frequencies
!=====================================================
  subroutine finitedif_vel(q0,ndn,vgr,evl0,evc0)
! calculates the group velocities in units of c_light from finite difference
! it would give zero near band crossings, thus HF is a better way to do it
  use constants
  use lattice, only : prim_to_cart
  use geometry
  implicit none
  integer ndn,i,l
  real(8) q0(3),vgr(3,ndn),q1(3),dq,v2(3,ndn)
  real(8) evl0(ndn),evlp(ndn),evlm(ndn) !,vpert(3,ndn)
  complex(8) evc0(ndn,ndn),evct(ndn,ndn)

  dq=1d-4  ! this is the ideal dq (at least for si)

  call get_freq(q0,ndn,v2,evl0,evc0)

! return

  do i=1,3
     q1=q0; q1(i)=q0(i)+dq
     call get_freq(q1,ndn,vgr,evlp,evct)
     q1=q0; q1(i)=q0(i)-dq
     call get_freq(q1,ndn,vgr,evlm,evct)
     vgr(i,:)=(evlp-evlm)/2/dq / 2/sqrt(abs(evl0)) *cnst*1d-8 *2*pi !*c_light
  enddo
!  write(30,6)'q,dq,om=',q0,dq,evl0
  write(330,4)'*la , v_FD ; v_PERT *** FOR q_red,qcart= ',matmul(transpose(prim_to_cart),q0)/(2*pi),' ***',q0
  do l=1,ndn
     write(330,5)'la ',l,vgr(:,l)*c_light,  &
&    length(vgr(:,l))*c_light,v2(:,l)*c_light,length(v2(:,l))*c_light
  enddo


 4 format(a,3(1x,f7.4),a,3(1x,g11.4))
 5 format(a,i5,2(5x,f11.3,3(1x,f11.3)))
 6 format(a,3(1x,f9.3),1x,g12.6,9(1x,f9.4))

  end subroutine finitedif_vel
!===========================================================
 subroutine get_freq(kp,ndn,vg,eival,eivec)
!! given a kpoint kp, it sets up the dynamical matrix and diagonalizes it
!! eival is in eV/uma/Ang^2
 use params
 use ios
 use constants
 use born
 use lattice, only : prim_to_cart
 implicit none
 integer, intent(in) :: ndn
 real(8), intent(in) :: kp(3)
 real(8), intent(out):: eival(ndn),vg(3,ndn)
 complex(8), intent(out) :: eivec(ndn,ndn)
 integer j,k,l,ier,nd2,al
 integer ,save :: nksave
 integer, allocatable :: mp(:)
 real(8), allocatable :: eivl(:)
 real(8) absvec,k_prim(3)
 complex(8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:),temp(:,:)
 real(8) khat(3)

 nd2 = min(ndn,12)
 allocate(dynmat(ndn,ndn),mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3),temp(ndn,ndn))
! write(uio,*)'# dynmat allocated, printing: Kpoint, eival(j),j=1,ndyn)'
! write(uio,*)'############ before setting up dynmat, kp number ',i

nksave=nksave+1
!write(*,3)'get_freq for ik,k=',nksave,kp

    call set_dynamical_matrix(kp,dynmat,ndn,ddyn)

!JS: call nonanalytical term for born effective change
    if (kp(1)==0.d0 .AND. kp(2)==0.d0 .AND. kp(3)==0.d0) then
        khat=kp(:)+1.0D-10
    else
        khat=kp(:)
    endif
    call nonanal(khat,dynmat,ndn,ddyn)

    if (verbose) then
       k_prim=matmul(transpose(prim_to_cart),kp)/(2*pi)
       write(ulog,3)' ===================================================='
       write(ulog,3)' THE DYNAMICAL MATRIX for KP(red,cart)=',ndn,k_prim,kp(:)
       do l=1,ndn
          write(ulog,8)(dynmat(l,j),j=1,nd2)
       enddo
       write(ulog,3)' ===================================================='
!      do al=1,3
!         write(ulog,3)' THE k-DERIVATIVE OF THE DYNAMICAL MATRIX =',al
!         do l=1,ndn
!            write(ulog,8)(ddyn(l,j,al),j=1,nd2)
!         enddo
!      enddo
    endif

! store dynmat in temp
    temp=dynmat
    call diagonalize(ndn,dynmat,eivl,ndn,eivc,ier)

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
    do l=1,ndn
       absvec = sum(abs(real(eivec(:,l))))
       if (absvec .lt. 1d-3) then
          eivec(:,l)=-cmplx(0,1)*eivec(:,l)
       endif
    enddo

  if (verbose) then

    do l=1,ndn
        write(ulog,3)' GET_FREQ:-----------  BAND ',l,eival(l)
        do j=1,ndn/3
           write(ulog,5)'atom ',j,eivec(3*(j-1)+1,l),eivec(3*(j-1)+2,l),eivec(3*(j-1)+3,l)
        enddo
    enddo

! now calculate the square of frequencies based on dynmat
!    dynmat=temp
!    temp = matmul(dynmat,eivc)
!    temp = matmul(transpose(conjg(eivc)),temp)
!    write(ulog,*)' TEST: below square of eivals should appear in the diagonal if e.D.e is right'
!    do j=1,ndn
!      write(ulog,9)'e.D.e=',j, temp(j,:)
!    enddo
!    write(ulog,*)' three components of ddyn are '
!    do al=1,3
!    do j=1,ndn
!      write(ulog,9)'al=',al, ddyn(j,:,al)
!    enddo
!    enddo

  endif

! now calculate the group velocities from Hellman-Feynmann formula(for each component al=1,3)
    vg=0
    do al=1,3

! first sandwich ddyn between eigenvectors
!      temp=matmul(ddyn(:,:,al),eivc)
!      dynmat=matmul(transpose(conjg(eivc)),temp)
!      do l=1,ndn
!         write(ulog,9)'al,d(Dyn)/dk=',al, dynmat(l,:)
!         vg(al,l)=dynmat(mp(l),mp(l))/(2*sqrt(abs(eival(l))))*cnst*1d-10*100*2*pi !*c_light
!      enddo

  
!   dynmat=matmul(ddyn(:,:,al),eivec)
!   dynmat=matmul(transpose(conjg(eivec)),dynmat)
!   vg(al,:)=0
    do l=1,ndn
      vg(al,l)=0
      do k=1,ndn
         dynmat(k,l)=0
         do j=1,ndn
            dynmat(k,l)=dynmat(k,l)+ddyn(k,j,al)*eivec(j,l)
         enddo
      enddo
      do k=1,ndn
          vg(al,l)=vg(al,l)+dynmat(k,l)*conjg(eivec(k,l))
      enddo
      vg(al,l)=dot_product(conjg(eivec(:,l)),matmul(ddyn(:,:,al),eivec(:,l)))
      vg(al,l)=vg(al,l)/2/sqrt(abs(eival(l))+1d-24)*cnst*1d-10*100*2*pi
    enddo
    enddo
 if (verbose) then
    write(ulog,*)' three components of the velocity of all bands in m/s are '
    do al=1,3
      write(ulog,5)'alpha,v_alpha(lambda) (pert)=',al,vg(al,:)*c_light
    enddo
 endif

 deallocate(eivl,eivc,dynmat,mp,ddyn,temp)

 3 format(a,i5,9(1x,f11.4))
 4 format(99(1x,2(1x,g9.3),1x))
 5 format(a,i5,99(1x,2(1x,g9.3),1x))
 6 format(a,99(1x,i5))
 7 format(i6,2x,3(1x,f7.3),1x,99(1x,f7.3))
 8 format(99(1x,2(1x,f9.3),1x))
 9 format(a,i5,99(1x,2(1x,f9.3),1x))

 end subroutine get_freq
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
 real(8), intent(in) :: kp(3,nkp),dk(nkp),eivl(ndn,nkp)
 complex(8), intent(in) :: eivc(ndn,ndn,nkp)
 complex(8), intent(out) :: grn(ndn,nkp)
 integer ik,i0,la,al,be,ga,j,k,j0,k0,ta1,ta2,t,cnt3,g,ti,ired
 real(8) mi,mj,rr3(3),rr2(3),qq(3),qdotr,omk !,mysqrt
 complex(8) zz,one,term

 one = cmplx(1d0,0d0)
! write(ugr2,*)'# la, nk , dk(nk) ,     kp(1:3,nk)     , om(la,nk),Re(gruneisen(la,nk)), Im(gr))'
 write(ugr2,*)'# nk , dk(nk) ,     kp(1:3,nk)     , om(la,nk),Re(gruneisen(la,nk)) '
  do ik=1,nkp
     qq(:) = kp(:,ik)
  do la=1,ndn
!    write(ulog,6)'Starting gruneisen parameters for kpoint & band# ',ik,la,kp(:,ik)
!    write(ulog,*)' fcs_3(t),term(2),6*eival(la,i),rr3(3),grun(la,i)'
     grn(la,ik) = 0
     omk = mysqrt(eivl(la,ik))*cnst
     do i0=1,natom_prim_cell
        mi = atom0(i0)%mass
        cnt3=0
        gloop: do g=1,map(3)%ngr
           if(g.gt.1) cnt3=cnt3+map(3)%ntind(g-1)  ! cnt3 counts all the terms in the previous groups
           tloop: do t=1,map(3)%nt(g)  !nterms(3)

              if ( i0 .ne. map(3)%gr(g)%iat(1,t) ) cycle tloop
              al = map(3)%gr(g)%ixyz(1,t)  !ixyzterm_3(1,t)
              be = map(3)%gr(g)%ixyz(2,t)  !ixyzterm_3(2,t)
              ga = map(3)%gr(g)%ixyz(3,t)  !ixyzterm_3(3,t)
              j  = map(3)%gr(g)%iat(2,t)   ! this atompos index
              k  = map(3)%gr(g)%iat(3,t)
              j0 = iatomcell0(j)     ! atom_sc(j)%cell%tau  is incorrect
              k0 = iatomcell0(k)     !atom_sc(k)%cell%tau is incorrect
              mj = atom0(j0)%mass
              ta1= al + 3*(i0-1)
              ta2= be + 3*(j0-1)
!             rr2 = atompos(:,j) - atom0(i0)%equilibrium_pos  ! Rij
              rr3(:) = atompos(:,k)                  ! R+tau
     !        rr2(:) = atompos(:,j) - atompos(:,i0)  ! R
     !        rr3(:) = atompos(:,k) - atompos(:,i0)  
     !        rr2(:) = atompos(:,j) - atompos(:,i0)  ! R  consistent with dynmat which uses j-i0
              rr2(:) = atompos(:,j) - atompos(:,j0)  ! R not consistent with dynmat which uses j-i0
              qdotr =  ( qq .dot. rr2)
              zz = cdexp( ci * qdotr )
              do ti=1,map(3)%ntind(g)  ! index of independent terms in that group g
                 ired=cnt3+ti + map(1)%ntotind + size_kept_fc2 !map(2)%ntotind
!! term = - ampterm_3(t)*fcs_3(igroup_3(t))*zz*eivc(ta1,la,i)*conjg(eivc(ta2,la,i))*rr3(ga)/sqrt(mi*mj)
                 term = -zz * eivc(ta2,la,ik)*conjg(eivc(ta1,la,ik))/sqrt(mi*mj) &
         &             * (fcs(ired) * rr3(ga)*map(3)%gr(g)%mat(t,ti))
                 grn(la,ik) = grn(la,ik) + term

! if (ik.eq.nkp .and. map(3)%gr(g)%mat(t,ti).ne.0) then
!    write(*,9)'t,ti,la,g,i0,ired,R2,R3,zz,fc,term=',t,ti,la,g,i0,ired,rr2,rr3,zz,fcs(ired),term
! endif
              enddo
           enddo tloop
        enddo gloop
     enddo

     grn(la,ik) = grn(la,ik)/6/eivl(la,ik)
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
9 format(a,6i4,66(1x,g10.4),6(1x,f9.4))
! deallocate(eivl,eivc,kg,grn)
 end subroutine gruneisen
!============================================================
 subroutine mechanical(elastic)
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 implicit none
 real(8), intent(out) :: elastic(6,6)
 integer i0,j0,al,be,ga,de,i,j,t,g,ti,cnt2,ired,voigt,la,nl,nc,nd,nu, &
&        longt,transz,ndn
 real(8)  rij(3),junk,junt,polt,pol,atld(3,3,3,3),ahat(3,3,3,3),c11,c12,c44,rho &
 &       ,ct(3,3,3,3),kapr(3*natom_prim_cell,3*natom_prim_cell),halfsum,halfdif  &
 &       ,gam(3*natom_prim_cell,3*natom_prim_cell)  &
 &       ,kap(3*natom_prim_cell,3*natom_prim_cell)  &
 &       ,temp(3*natom_prim_cell-3,3*natom_prim_cell-3),q(3)  &
 &       ,vg(3,3*natom_prim_cell),eivl(3*natom_prim_cell),eivl0(3*natom_prim_cell)
 complex(8) eivc(3*natom_prim_cell,3*natom_prim_cell),eivc0(3*natom_prim_cell,3*natom_prim_cell)
 real(8), pointer:: gama(:,:)

 write(ulog,*)' ********** ENTERING MECHANICAL *************' 
 write(ulog,*)' Elastic constant of a cubic crystal from group velocities '
 ndn=3*natom_prim_cell
 rho=sum(atom0(:)%mass)*uma/volume_r0*1d30
 write(ulog,*)'density in kg/m^3 is=',rho
 q=(/1d-9,1d-9,0d0/)
 call get_freq(q,ndn,vg,eivl0,eivc0)

 q=(/1d-3,1d-3,0d0/)
! from FD or perturbation group velocities
 call get_freq(q,ndn,vg,eivl,eivc); eivl=eivl-eivl0

! eivl is the square of phonon frequencies eivl=omega^2 -> convert to eV
      eivl=eivl/(q.dot.q)*rho*volume_r0 / (uma*1d30)
      call write_out(ulog,'mass_cell w^2 /q^2 (eV) ',eivl)
! now eivl is in eV
!     call get_cij(eivl,eivc,ndn,c11,c12,c44)
!     call write_out(ulog,'GET_CIJ:c11,c12,c44(GPa)',(/c11,c12,c44/))

 write(ulog,*)'Vg_pert is=',(length(vg(:,j))*c_light,j=1,3)
      c44    =rho*(length(vg(:,1))*c_light)**2 *1d-9  ! assumes eivc(z,1)=1
      halfdif=rho*(length(vg(:,2))*c_light)**2 *1d-9
      halfsum=rho*(length(vg(:,3))*c_light)**2 *1d-9 - c44
      c12=halfsum-halfdif
      c11=halfsum+halfdif
      call write_out(ulog,'from rho vg^2 if e_1=z :',(/c11,c12,c44/))
      c44    =rho*(length(vg(:,2))*c_light)**2 *1d-9  ! assumes eivc(z,2)=1
      halfdif=rho*(length(vg(:,1))*c_light)**2 *1d-9
      halfsum=rho*(length(vg(:,3))*c_light)**2 *1d-9 - c44
      c12=halfsum-halfdif
      c11=halfsum+halfdif
      call write_out(ulog,'from rho vg^2 if e_2=z :',(/c11,c12,c44/))
 call write_out(ulog,'GF dyn_mat eivals (eV/Ang^2) ',eivl)
 call finitedif_vel(q,ndn,vg,eivl,eivc); eivl=eivl-eivl0
 write(ulog,*)'Vg_FD   is=',(length(vg(:,j))*c_light,j=1,3)
      c44    =rho*(length(vg(:,1))*c_light)**2 *1d-9  ! assumes eivc(z,1)=1
      halfdif=rho*(length(vg(:,2))*c_light)**2 *1d-9
      halfsum=rho*(length(vg(:,3))*c_light)**2 *1d-9 - c44
      c12=halfsum-halfdif
      c11=halfsum+halfdif
      call write_out(ulog,'from rho vg^2 if e_1=z :',(/c11,c12,c44/))
      c44    =rho*(length(vg(:,2))*c_light)**2 *1d-9  ! assumes eivc(z,2)=1
      halfdif=rho*(length(vg(:,1))*c_light)**2 *1d-9
      halfsum=rho*(length(vg(:,3))*c_light)**2 *1d-9 - c44
      c12=halfsum-halfdif
      c11=halfsum+halfdif
      call write_out(ulog,'from rho vg^2 if e_2=z :',(/c11,c12,c44/))
 call write_out(ulog,'FD dyn_mat eivals (eV/Ang^2) ',eivl)
 write(ulog,*)' ********************************************' 



! Standard Bravais term here: 
 atld=0 ; kapr=0; kap=0
 do al=1,3
 do ga=1,3

 do i0=1,natom_prim_cell
 cnt2=0  ! cumulative number of inependent terms up to group (g-1)
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
          if ( i0 .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
          if( map(2)%gr(g)%ixyz(2,t) .ne. ga) cycle tloop
          j  = map(2)%gr(g)%iat(2,t)
          j0 = iatomcell0(j)    ! atom_sc(j)%cell%tau  is incorrect
          rij = atompos(:,j) - atom0(i0)%equilibrium_pos  ! Rij
          junk = fcs(ired)* map(2)%gr(g)%mat(t,ti)  ! that is phi(al,ga)
          nl=al+3*(i0-1) ; nc=ga+3*(j0-1)  ! dynmat dimensions
          kap(nl,nc) =kap(nl,nc)+junk
          do la=1,3
             nu=3*(j0-1)+la
             kapr(nl,nu)=kapr(nl,nu)+junk*rij(la)
          enddo
 do be=1,3
 do de=1,3
          atld(al,be,ga,de) = atld(al,be,ga,de) -0.25* junk*(rij(be)*rij(de)+rij(de)*rij(be))   
 enddo
 enddo
            
       enddo tloop
    enddo
 enddo gloop
 enddo

 enddo
 enddo
 atld=atld/volume_r0
 gam=kap  ! force constants and the images sum_R K(tau,R+taup)
 kap=0.5*(kap+transpose(kap))/volume_r0   ! symmetrize
 call write_out(ulog,' sum_R phi(tau,R+taup) ',kap)
 call write_out(ulog,' sum_R phi(tau,R+taup) (R+taup-tau) ',kapr)
 nd=3*natom_prim_cell
 allocate(gama(nd-3,nd-3))
 gama=gam(4:nd,4:nd)
 
 if (nd.eq.3) return  ! it's a bravais lattice and the first part is sufficient

 call inverse_real(gama,temp,nd-3)
 gam(1:3,:)=0 ; gam(:,1:3)=0
 gam(4:nd,4:nd)=temp
 
! for the second term, we invert sum_N fc2(tau,N+taupp) in gamma and use it to get X
 kapr=matmul(gam,kapr)  ! this is the X matrix=d u0/d eta

! check symmetry
 do al=1,3
 do be=1,3
 do ga=1,3
 do de=1,3
!   if (abs(atld(al,be,ga,de)-atld(ga,be,al,de)).gt.1d-4)  &
!&       write(ulog,*)al,be,ga,de,atld(al,be,ga,de),atld(ga,be,al,de)
!   if (abs(atld(al,be,ga,de)-atld(ga,de,al,be)).gt.1d-4)  &
!&       write(ulog,*)al,be,ga,de,atld(al,be,ga,de),atld(ga,be,al,de)
    ahat(al,ga,be,de)=0.5*(atld(al,be,ga,de)+atld(ga,be,al,de))
 enddo
 enddo
 enddo
 enddo
    
 do al=1,3
 do be=1,3
 do ga=1,3
 do de=1,3
 !  ct(al,be,ga,de)=ahat(al,ga,be,de)+ahat(be,ga,al,de)-ahat(al,be,ga,de)
    ct(al,be,ga,de)=(ahat(al,ga,be,de)+ahat(ga,al,be,de))/4+(ahat(al,ga,be,de)+ahat(al,ga,de,be))/4
!   write(ulog,*)al,be,ga,de,ct(al,be,ga,de)
 enddo
 enddo
 enddo
 enddo

 ct=ct*ee*1d30*1d-9

 do al=1,3
 do be=al,3
 do ga=1,3
 do de=ga,3
    i=voigt(al,be); j=voigt(ga,de)
    elastic(i,j)=ct(al,be,ga,de)
!   if (abs(elastic(i,j)) .gt. 1d-4 ) write(ulog,*)'i,i,elastic=',i,j,elastic(i,j)
 enddo
 enddo
 enddo
 enddo
   call write_out (ulog,' Elastic Tensor ',elastic)
   call write_out (6,' Elastic Tensor ',elastic)
!do t=1,map(2)%nt(g)   !nterms(2)
!   if ( i0 .ne.  map(2)%gr(g)%iat(1,t) cycle  !iatomterm_2(1,t) ) cycle

!   al = map(2)%gr(g)%ixyz(1,t) !ixyzterm_2(1,t)
!   be = map(2)%gr(g)%ixyz(2,t) !ixyzterm_2(2,t)
!   j  = map(2)%gr(g)%iat(2,t) !iatomterm_2(2,t)
!   rija = atompos(al,j)-atompos(al,i0)
!   rijb = atompos(be,j)-atompos(be,i0)
!!  write(ulog,3)i0, al, j, be, rija,rijb,fcs_2(t),bulk
!   bulk = bulk - fcs_2(t)*(1-2*grun_fc(t)*dlogv)*rija*rijb
!   if (al.eq.1 .and. be.eq.1) then
!      c11 = c11 - fcs_2(t)*(1-2*grun_fc(t)*dlogv)*rija*rijb
!   endif
!enddo
!enddo
!bulk = bulk/volume_r0/18
!c11 = c11/volume_r0/2

!write(ulog,*)'BULK_MODULUS: in eV/A^3 ',bulk
!bulk = bulk*1d30*ee
!write(ulog,*)'BULK_MODULUS: in SI, Mbar ',bulk,bulk*1d-11
!c11 = c11*1d30*ee
!write(ulog,*)'C11 in SI, Mbar ',c11,c11*1d-11
!3 format(4(i4),9(2x,f10.5))

 end subroutine mechanical
!============================================================
 subroutine sound_speeds !(elastic)
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 implicit none
 ! real(8), intent(out) :: elastic(6,6)
 integer i0,j0,al,be,ier,j,t,g,ti,cnt2,ired,mu,nu,ndn
 real(8)  rj(3),junk,khat_red(3)
 real(8)     eivl(3*natom_prim_cell),c11,c12,c44
 complex(8)  eivc(3*natom_prim_cell,3*natom_prim_cell)
 complex(8) dynr2(3*natom_prim_cell,3*natom_prim_cell)

 ndn=3*natom_prim_cell
 write(ulog,*)' ********** ENTERING SOUND SPEEDS *************' 
 khat_red=(g01+g02)/length(g01+g02)
 dynr2=0
! calculate the matrix sum Phi (R.khat)^2/2 for k along (110)
 do al=1,3
 do be=1,3
 do i0=1,natom_prim_cell
 cnt2=0  ! cumulative number of inependent terms up to group (g-1)
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
!   call get_cij(eivl,eivc,ndn,c11,c12,c44)
!   call write_out(ulog,'From sound_speeds:c11,c12,c44(GPa)',(/c11,c12,c44/))
 
 end subroutine sound_speeds
!============================================================
 subroutine thermo(nk,wk,ndyn,eival,grn,temp,etot,free,pres,pres0,cv)
! calculate total and free energies within QHA, at a given temperature (temp in Kelvin)
 use ios
 use params
 use lattice
 use constants
 implicit none
 integer nk,ndyn,b,k,nat
 real(8) wk(nk),eival(ndyn,nk)
 real(8) temp,x,cv_nk,cv,hw,free,etot,pres,nbe,mdedv,pres0,nbx
 complex(8) grn(ndyn,nk)

    nat = ndyn/3
    if (temp.le.0) then
       write(ulog,*)'temperature not in the proper range!!',temp
       stop
    endif
    etot=0 ; cv=0 ; free=0 ; pres=0
    mdedv= 0.35388/(20.8**3-20.0**3)/ab**3 ! this is -dE/dV in eV/A^3
    mdedv = mdedv*1d+30*ee ! (in J/m^3 )
    do k=1,nk
    do b=1,ndyn
       if(eival(b,k) .lt.0) then
          x=0
          write(ulog,*) 'ALPHA: negative eival for band,kp# ',b,k,eival(b,k)
          write(ulog,*) ' will use its absolute value instead!'
!      else
       endif
! to convert eV/A^2/mass to Hz we need to take the sqrt, multiply by cnst*100*c_light
       x=(h_plank*sqrt(abs(eival(b,k)))*cnst*100*c_light)/k_b/temp
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
       pres= pres + hw * (nbe(x,1d0,classical) + 0.5d0) * wk(k) * grn(b,k)
    enddo
    enddo
    cv = cv/nat*n_avog*k_b  ! this is per mole    )/(volume_r0*1d-30)  !
    etot = etot/nat*n_avog   ! convert from Joule/cell to Joule per mole
    free = free/nat*n_avog   ! convert from Joule/cell to Joule per mole
    pres0= pres/(volume_r0*1d-30) ! * 1d-8  ! 1d-8 is to convert to kbar
    pres = pres0+ mdedv
3 format(9(2x,g11.5))

 end subroutine thermo
!=======================================================
 subroutine check_mdyn(ndyn,nk,kp,eivec,eival)
! form the dynamical matrix from existing eigenvectors and eigenvalues
! and diagonalize and verify consistency :eivec(-q)=-conjg(eivec(q))
 use lattice  ! needed to access NC,g1,g2,g3
 use kpoints
 implicit none
 integer ndyn,nk,j,k,l,i,ier,i1,j1,k1,mi,inside
 real(8) kp(3,nk),eival(ndyn,nk),eivl(ndyn)
 complex(8) eivec(ndyn,ndyn,nk),dynm(ndyn,ndyn),eivc(ndyn,ndyn),d2(ndyn,ndyn)

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
 real(8) kp(3,nk),w(3)
 complex(8) eigenvec(ndyn,ndyn,nk)

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
 complex(8), intent(inout) :: dynmat(ndim,ndim),ddyn(ndim,ndim,3)
 real(8), intent(in) :: q(3)
 real(8) zag(3),zbg(3),qeq,ma,mb,rr(3),dqeq(3)
 integer na,nb,i,j,al
 real(8) eps_scale,q2,term,asr(3,3),dyn_coul(3,3) ,sf,ddn(3,3,3)

! add anyways
! if(born_flag.eq.0) return   

! rho=1d6
! eps_scale=8.8541878176D-12/1D10/ee
! gg=2*pi/(volume_r0)**0.33
! if ( gg .lt. 4*rho) then
! if (born_flag.eq.1)  rho = gg/4d0
! write(30,*)'ADOPTED VALUE of RHO = ',rho
! endif

!  dynmat=0; ddyn=0 they are already initialized in set_dynmat
  do na=1,natom_prim_cell
     ma = atom0(na)%mass
!    write(ulog,*)'na,ma=',na,ma
     do nb=1,natom_prim_cell
        mb = atom0(nb)%mass
!    write(ulog,*)'nb,mb=',nb,mb
        call dyn_coulomb(na,nb,q,dyn_coul,ddn)

!    call write_out(ulog,' NONANAL: dyn_coul ',dyn_coul)
        do i = 1,3
        do j = 1,3
           dynmat(i+3*(na-1),j+3*(nb-1)) = dynmat(i+3*(na-1),j+3*(nb-1))+ dyn_coul(i,j)/sqrt(ma*mb)
           do l=1,3
              ddyn(i+3*(na-1),j+3*(nb-1),l) = ddyn(i+3*(na-1),j+3*(nb-1),l)+ ddn(i,j,l)/sqrt(ma*mb)
           enddo
        enddo     
        enddo     

     enddo
  enddo


! term=1d0/(qeq*eps_scale)/volume_r0 *exp(-qeq/(rho*rho))
!
! do na = 1,natom_prim_cell
!    ma = atom0(na)%mass
!    zag = matmul(atom0(na)%charge,q)
!    do nb = 1,natom_prim_cell
!       mb = atom0(nb)%mass
!       zbg = matmul(atom0(nb)%charge,q)
!       rr = atompos(:,na)-atompos(:,nb)
!       do i = 1,3
!       do j = 1,3
!          dynmat(i+3*(na-1),j+3*(nb-1)) = dynmat(i+3*(na-1),j+3*(nb-1))+ &
!      &     zag(i)*zbg(j)*term/sqrt(ma*mb)
! ALSO NEED TO CORRECT THE GROUP VELOCITIES
!         do al=1,3
!            ddyn(i+3*(na-1),j+3*(nb-1),al) = ddyn(i+3*(na-1),j+3*(nb-1),al)+term/sqrt(ma*mb) *  &
! &          (atom0(na)%charge(al,i)*zbg(j)+atom0(nb)%charge(al,j)*zag(i)-dqeq(al)*zag(i)*zbg(j)/qeq)
!         enddo
!       end do
!       end do
!    end do
! end do


! check ASR
!if(length(q).le.1d-4) then
!  do na=1,natom_prim_cell
!     asr=0
!     do nb=1,natom_prim_cell
!        call dyn_coulomb(na,nb,q,dyn_coul)
!        asr=asr+dyn_coul
!     enddo     
!
!     do i=1,3
!     do j=1,3
!        if(abs(asr(i,j)).ge.1d-4)  write(ulog,*)'NONANAL:ASR broken ',i,j,asr(i,j)
!     enddo
!     enddo
!  enddo
!endif

 end subroutine nonanal
!===========================================================
 subroutine enforce_asr_simple(n,dyn)
 implicit none
 integer n,nat,i,j,io,jo,lo,iat,jat
 complex(8) dyn(n,n),sum

 nat=n/3
 do i=1,3
 do iat=1,nat
    io=i+3*(iat-1)
    do j=1,3
       sum=0
       jo=j+3*(iat-1)
       do jat=1,nat
          if (iat.eq.jat) cycle
          lo=j+3*(jat-1)
          sum=sum-dyn(io,lo)
       enddo
       dyn(io,jo)=sum
    enddo
 enddo
 enddo
 end subroutine enforce_asr_simple
!===========================================================
subroutine band_sort_bs(nkp,ndyn,kp,eivec,emap)
implicit none

integer, intent(in) :: nkp,ndyn
real(8), intent(in) :: kp(3,nkp)
!real(8), intent(in) :: eival(ndyn,nkp),dk(nkp)
complex(8), intent(in) :: eivec(ndyn,ndyn,nkp)
integer, intent(out) :: emap(ndyn,nkp)

integer i,j,k,l
integer ntbd !Number of connections to be determined


complex(8), allocatable :: overlap(:,:)
real(8), allocatable :: overlap_q(:,:)

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
        overlap_q(j,k)=overlap(j,k)*dconjg(overlap(j,k))
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
 subroutine set_dynamical_matrix(kpt,dynmat,ndim,ddyn)
 use params
 use lattice
 use kpoints
 use atoms_force_constants
 use geometry
 use svd_stuff
 implicit none
 integer, intent(in) :: ndim
 complex(8), intent(out) :: dynmat(ndim,ndim) ,ddyn(ndim,ndim,3)
 real(8), intent(in) :: kpt(3)
 complex(8) junk
 real(8) mi,mj,nrm1,rr(3) !,delt(3)
 integer i0,j,j0,al,be,i3,j3,t,cnt2,ti,ired,g

4 format(a,4(1x,i5),9(2x,f9.4))
9 format(i6,99(1x,f9.3))

! delt = wshift !*wshift
 ddyn   = cmplx(0d0,0d0)
 dynmat = cmplx(0d0,0d0)
 do i0=1,natom_prim_cell 
!   write(*,3)'i0,ri0=',i0,atom0(i0)%equilibrium_pos ,atompos(:,i0)
 do al=1,3
    i3 = al+3*(i0-1)
    mi = atom0(i0)%mass
! write(ulog,*) 'i,al,mass=',i0,al,i3,mi
    cnt2=0   ! counter of ntotind, to find position of each term in the array fc2
    gloop: do g=1,map(2)%ngr

       if(keep_grp2(g).ne.1) cycle gloop
       if(g.gt.1) cnt2=cnt2+map(2)%ntind(g-1) 
       do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
          ired=cnt2+ti + map(1)%ntotind  ! the latter is the size of fc1
          if(cnt2+ti .gt. size_kept_fc2) then
             write(*,*)'SET_DYNMAT: size of fc2 exceeded, ired=',ired
             stop
          endif

       tloop: do t=1,map(2)%nt(g)
          if(keep_grp2(g).ne.1) cycle tloop
          if ( i0 .ne. map(2)%gr(g)%iat(1,t) .or. al .ne. map(2)%gr(g)%ixyz(1,t) ) cycle tloop
          be = map(2)%gr(g)%ixyz(2,t)
          j  = map(2)%gr(g)%iat(2,t)
          j0 = iatomcell0(j)    ! atom_sc(j)%cell%tau  is incorrect
          mj = atom0(j0)%mass
          j3 = be+3*(j0-1)
!         rr = atompos(:,j) - atom0(i0)%equilibrium_pos  ! Rij
          rr = atompos(:,j) - atompos(:,j0)
             junk = fcs(ired)* exp(ci*(kpt .dot. rr))/sqrt(mi*mj)*map(2)%gr(g)%mat(t,ti)
             dynmat(i3,j3) = dynmat(i3,j3) + junk
             ddyn(i3,j3,:) = ddyn(i3,j3,:) + junk*ci*rr(:)
       enddo tloop
       enddo
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
!      stop
!   else
       dynmat(t,t) = cmplx(real(dynmat(t,t)),0d0)
    endif
  do j=t+1,ndim-1
    if (abs(aimag(dynmat(t,j))+aimag(dynmat(j,t))) .gt. 1d-5*nrm1 ) then
       write(ulog,*)' dynmat is not hermitian in AIMAG of its off-diagonal elts'
       write(ulog,*)' off-diagonal element i,j=',t,j,dynmat(t,j),dynmat(j,t)
       write(ulog,*)' comparing to max(abs(dynmat))=',nrm1
!      stop
    elseif(abs(real(dynmat(t,j))-real(dynmat(j,t))) .gt. 1d-5*nrm1 ) then
       write(ulog,*)' dynmat is not hermitian in REAL of its off-diagonal elts'
       write(ulog,*)' off-diagonal element i,j=',t,j,dynmat(t,j),dynmat(j,t)
       write(ulog,*)' comparing to max(abs(dynmat))=',nrm1
!      stop
    else
    endif
! enforcing it to be hermitian!!
    mi=(aimag(dynmat(t,j))-aimag(dynmat(j,t)))/2
    mj=(real (dynmat(t,j))+real (dynmat(j,t)))/2
    dynmat(t,j) = cmplx(mj, mi)
    dynmat(j,t) = cmplx(mj,-mi)
  enddo
 enddo

! enforce it to be real if all imaginary componenents are very small
! all = sum(abs(aimag(dynmat(:,:))))/(ndim*ndim)
! if (all .lt. sum(abs(real(dynmat(:,:))))/(ndim*ndim)*1d-10) then
!    dynmat=(dynmat+conjg(transpose(dynmat)))/2d0
! endif

 end subroutine set_dynamical_matrix
!==========================================================
 subroutine diagonalize(n,mat,eival,nv,eivec,ier)
! n=size of mat; nv is the number of needed eigenvectors
 implicit none
 integer, intent(in) :: n,nv
 integer, intent(out) :: ier
 complex(8), intent(in) :: mat(n,n)
 complex(8), intent(out) :: eivec(n,n)
 real(8), intent(out) :: eival(n)
! This is used by eigch
! real(8), allocatable :: w(:,:)
! integer, allocatable :: lw(:)
! This is used by ZHEGV
  real(8), allocatable :: rwork(:)
  complex(8), allocatable :: work(:)
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
 real(8), intent(in) :: dk(nk),kp(3,nk),eival(n,nk),vg(3,n,nk)
 integer j,k

 do k=1,nk
 do j=1,n
    write(uio,3)k,j,dk(k),kp(:,k),cnst*mysqrt(eival(j,k)),vg(:,j,k)*c_light,length(vg(:,j,k))*c_light
 enddo
 enddo

2 format(i5,1x,4(1x,f8.3),999(1x,g11.4))
3 format(2i5,1x,4(1x,f8.3),999(1x,g11.4))

 end subroutine write_all_eigenvalues
!==========================================================
 subroutine write_eigenvalues(i,dk,kp,eival,vg,n,uio)
 use constants
 use eigen, only : mysqrt
 use geometry
 implicit none 
 integer, intent(in) :: i,n,uio
 real(8), intent(in) :: dk,kp(3),eival(n),vg(3,n)
 integer j

 write(uio,2)i,dk,kp,(cnst*mysqrt(eival(j)),j=1,n) ,(vg(:,j)*c_light,length(vg(:,j))*c_light,j=1,n)

! do j=1,n
!    write(uio+1,3)i,j,dk,kp,cnst*mysqrt(eival(j))
! enddo

2 format(i5,1x,4(1x,f8.3),999(1x,g11.4))
3 format(2i5,1x,4(1x,f8.3),999(1x,g11.4))

 end subroutine write_eigenvalues
!==========================================================
 subroutine calculate_dos(mx,eival,wkp,mesh,omega,ds)
! use gaussian broadening
 use ios
 use params
 use om_dos
 use constants
 implicit none
 integer i,j,mx,mesh
 real(8) x,wkp(mx),delta,ds(mesh),omega(mesh),eival(mx)  !,cnst

! wkp=1d0/(nkx*nky*nkz)
! write(udos,*)'# wkp,width =',wkp,width

    do i=1,mesh
       ds(i) = 0
       do j=1,mx
          x = (eival(j) - omega(i))/width   ! these are all in cm^-1
!         x = (eival(j) - ene*ene)/width/width
          if ( abs(x) .gt. 5 ) cycle
          ds(i) = ds(i) + delta(x)/width*wkp(j) !/mx
       enddo
    enddo

3 format(i5,9(3x,g11.5))

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
 real(8) xp,xm,wkp(mx),delta_g,dsp,dsm,omega(mesh),eival(mx),iself,tmp,one
 complex(8) omz

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
       omz=2*cmplx(omega(i),-etaz)
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

2 format(9(3x,g11.5))
3 format(i5,9(3x,g11.5))

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
 real(8) q(3),q1(3),qp(3),qm(3),dosp,dosm,om1,omp,omm,omq,delta_l

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
 real(8) q(3),omq,dosp,dosm,eival(ndyn,nk) 
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

3 format(i5,9(3x,g10.4))

 end subroutine phase_space_lifetime
!===========================================================
  subroutine cubic_on_gconv(coord)
! direct coords of conv cubic on present conv
  use lattice
  use ios
  use geometry
 implicit none
 !real(8), intent(out) :: coord(3,3)
 real(8) coord(3,3)

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
!===========================================================
 subroutine dyn_coulomb(tau,taup,q,dyn,ddyn)
!! calculates the q component of the Coulomb force to subtract from fourier transform of total force
! sdyn(tau,taup) = [z(tau)^al,ga.q^ga]  [z(taup)^be.de.q^de] * sf(q) &
!    &          1/eps0/(q.epsr.q)/volume_r0 
! DOES NOT INCLUDE THE MASS DENOMINATOR
 use atoms_force_constants
 use lattice, only : volume_r , volume_r0
 use fourier  , only : nrgrid,nggrid,rgrid,ggrid,rws_weights,gws_weights
 use born, only : epsil ,born_flag 
 use constants , only : eps0,ee,pi
 use ios, only : ulog,write_out
 implicit none
 integer, intent(in) :: tau,taup
 real(8), intent(in) :: q(3) 
 real(8), intent(out) :: dyn(3,3),ddyn(3,3,3)  ! element_(tau,taup)(q) of the Coulomb dynamical mat
 integer al,be,i
 real(8) zqa(3),zqb(3),denom,coef,mt(3,3),sf,dsf(3),qeq,dqeq(3),trn,trd
 real(8), allocatable::ep1(:)

    dyn=0; ddyn=0
! if(born_flag.eq.0) then
!    return ! rho = 2*pi/(volume_r0)**0.33 / 4d0
! endif

 coef=ee*1d10/eps0/volume_r0 ! to convert 1/ang^2 to eV/ang , includes cancellation of 4pi
! call write_out(6,' DYN_COULOMB: epsil ',epsil)
  denom=sum( (/ (epsil(i,i),i=1,3)/) )   ! this is the trace function
! call write_out(6,' DYN_COULOMB: trace(epsil) ',denom)
 if(length(q).lt.1d-6)then

   mt=matmul(atom0(tau)%charge,atom0(taup)%charge)
! call write_out(6,' DYN_COULOMB: Z_tau Z_taup ',mt)
   dyn=0
   do al=1,3
   do be=1,3
      trn=sum( (/ (mt(i,i),i=1,3)/) ) 
      trd=sum( (/ (epsil(i,i),i=1,3)/) ) 
      dyn(al,be)=trn/trd*coef !trace(mt)/trace(epsil)*coef
      ddyn(al,be,:)= 0
   enddo
   enddo
!call write_out(ulog,' DYN_COULOMB: dyn(q=0) ',dyn)

 else  

   call structure_factor(q,nrgrid,rgrid,rws_weights,sf,dsf)
   qeq=dot_product(q,matmul(epsil,q))
   dqeq=matmul((epsil+transpose(epsil)),q)
   zqa=matmul(atom0(tau)%charge,q)
   zqb=matmul(atom0(taup)%charge,q)
!call write_out(ulog,' DYN_COULOMB: sf  ',sf )
!call write_out(ulog,' DYN_COULOMB: qeq ',qeq)

   do al=1,3
   do be=1,3
      dyn(al,be)=zqa(al)*zqb(be)/qeq * coef * sf 
      do i=1,3
         ddyn(al,be,i) = ((atom0(tau)%charge(al,i)*zqb(be)+zqa(al)*atom0(taup)%charge(be,i))*sf +  &
 &                        zqa(al)*zqb(be)*(dsf(i)-sf*dqeq(i)/qeq))/qeq*coef
      enddo
   enddo
   enddo

! call write_out(6,' DYN_COULOMB: dyn ',dyn)
 endif

 end subroutine dyn_coulomb
!============================================================
 subroutine get_cij(eivl,eivc,ndn,c11,c12,c44)
!! reads eigenvalues of L_ik = C_ij,kl qhat_j qhat_l in eV along 110 and 
!! outputs the elastic constants of a CUBIC crystal
 use atoms_force_constants , only : natom_prim_cell
 use lattice, only : volume_r0
 use constants, only : uma,c_light,ee
 use ios !, only : ulog
 implicit none
 integer, intent(in) :: ndn
 real(8), intent(inout) :: eivl(ndn)
 complex(8), intent(in) :: eivc(ndn,ndn)
 real(8), intent(out) :: c11,c12,c44
 real(8) pol,polt,junk,junt,q(3),halfdif,halfsum
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
