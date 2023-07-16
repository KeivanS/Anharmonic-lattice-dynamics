 subroutine set_dynmat(fcs,ngr,dyn)
!! reads irreducible 2nd order force constants, and computes the dynamical matrix
 use svd_stuff, only : map
 use kpoints
 use atoms_force_constants, only : natom_prim_cell
 use ios, only : ulog
  integer, intent(in) :: ngr
  real(8), intent(in) :: fcs(ngr)
  complex(8), allocatable, intent(inout) :: dyn(:,:)
  integer ndyn

  ndyn = 3*natom_prim_cell

!---------------------------------------------------------------
! do a band structure calculation along the symmetry directions

  call make_kp_bs

  write(ulog,*)' Kpoints for band structure generated from kpbs.in'
  call allocate_eig_bs(ndyn,nkp_bs,ndyn) !3rd index is for eivecs needed for gruneisen

! now output eigenval is directly frequencies in 1/cm
!  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,ndyn,eigenvec_bs,veloc)

 end subroutine set_dynmat
!===========================================================
 subroutine get_frequencies(nkp,kp,dk,ndn,eival,nv,eivec,vg)
! also outputs the group velocities from HF theorem (exact) in units of c_light
 use params
 use constants
 use born
 use geometry
 implicit none
 integer, intent(in) :: nkp,ndn,nv ! no of wanted eivecs
 real(8), intent(in) :: kp(3,nkp),dk(nkp)
 real(8), intent(out):: eival(ndn,nkp),vg(3,ndn,nkp)
 complex(8), intent(out) :: eivec(ndn,ndn,nkp)
 integer i,j,k,l,ier,nd2,al,ll
 integer, allocatable :: mp(:) !Map for ascending band sorting
 real(8), allocatable :: eivl(:)
 real(8) absvec,om
 complex(8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:)
 real(8) khat(3),mysqrt
 character ext*3

 integer, allocatable :: mp1(:,:) !Map for projection band sorting
 real(8), allocatable :: eival_tmp(:) !Temp matrices for reordering eigenvalues
 complex(8), allocatable :: eivec_tmp(:,:)
 real(8), allocatable :: vg_tmp(:,:)

! open files and write the group velocities
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
    call finitedif_vel(kp(:,i),ndn,vg(:,:,i),eivl,eivc)

!   call set_dynamical_matrix(kp(:,i),dynmat,ndn,ddyn)
!
!!JS: call nonanalytical term for born effective change
!    if (kp(1,i)==0.d0 .AND. kp(2,i)==0.d0 .AND. kp(3,i)==0.d0) then
!        khat=kp(:,i)+1.0D-10
!    else
!        khat=kp(:,i)
!    endif
!    call nonanal(khat,dynmat,ndn,ddyn)
!
!    if (verbose) then
!       write(ulog,3)' ======================================================================'
!       write(ulog,3)' THE DYNAMICAL MATRIX IS:',i,kp(:,i)
!       do l=1,ndn
!          write(ulog,4)(dynmat(l,j),j=1,nd2)
!       enddo
!    endif
!
!    call diagonalize(ndn,dynmat,eivl,nv,eivc,ier)
!
! sort eivals in ascending order
!   call sort(ndn,eivl(:),mp,ndn)

!   if (ier.ne.0 .or. verbose) then
!     write(   *,*)' ier=',ier
!     write(ulog,*)' ier=',ier, ' EIGENVALUES ARE...'
!     write(ulog,*)
!   endif

    do j=1,ndn
       eival(j,i) = eivl(j) !mp(j))
       do l=1,nv
          eivec(j,l,i) = eivc(j,l)  !mp(l))
       enddo
    enddo

! if pure imaginary make eigenvectors real as they can be multiplied by any arbitrary number
!    do l=1,nv
!       absvec = sum(abs(real(eivec(:,l,i)))**2)
!       if (absvec .lt. 1d-1) then
!          eivec(:,l,i)=cmplx(0,1)*eivec(:,l,i)   ! this does not hurt anything
!       endif
!    enddo
!
!! put the acoustic bands in the first 3 branches 1:Ta 2:Ta 3:La
!! eivl(mp(i)) is the properly sorted array
!
!!   call sort_polarizations(ndn,eival(:,i),eivec(:,:,i),mp,kp(:,i))
!    do j=1,ndn
!       eivl(j)=eival(mp(j),i)
!       do l=1,nv
!          eivc(j,l) = eivec(j,mp(l),i)
!       enddo
!    enddo
!    call write_eigenvalues(i,dk(i),kp(:,i),eivl,ndn,uio)
!
!    if (verbose) then
!      do l=1,nv
!        write(ulog,3)' -----------  BAND ',l,eivl(l)
!        ll=l !mp(l)
!        do j=1,ndn/3
!           write(ulog,5)'atom ',j,eivec(3*(j-1)+1,ll,i),eivec(3*(j-1)+2,ll,i),eivec(3*(j-1)+3,ll,i)
!!          write(ulog,5)'atom ',j,eivc(3*(j-1)+1,l),eivc(3*(j-1)+2,l),eivc(3*(j-1)+3,l)
!         enddo
!      enddo
!    endif

! now calculate the group velocities from Hellman-Feynmann formula(dynmat=dummy)
!    do al=1,3
!    do l=1,ndn
!      vg(al,l,i)=0
!      do k=1,ndn
!         dynmat(k,l)=0
!         do j=1,ndn
!            dynmat(k,l)=dynmat(k,l)+ddyn(k,j,al)*eivec(j,l,i)
!         enddo
!      enddo
!      do k=1,ndn
!         vg(al,l,i)=vg(al,l,i)+dynmat(k,l)*conjg(eivec(k,l,i))
!      enddo
!      vg(al,l,i)=vg(al,l,i)/2/sqrt(abs(eival(l,i)))*cnst*1d-10*100*2*pi
!    enddo
!    enddo
!
!    do l=1,ndn
!       om=sqrt(abs(eival(l,i)))*cnst
!       write(uvel+l,6)l,i,kp(:,i),om,vg(:,l,i)*c_light,length(vg(:,l,i))*c_light
!    enddo
!
!!   do j=1,ndn
!!      eival(j,i)=eivl(j)
!!      do l=1,nv
!!         eivec(j,l,i)=eivc(j,l)
!!      enddo
!!   enddo
!
!     if(ier .ne. 0) stop
! write(uio,*)'==============='
! write(uio,4)(eivl(j),j=1,ndn)
! write(uio,*)'==============='
! do l=1,ndn
!    write(uio,4)(eivc(l,j),j=1,ndn)
! enddo
!   call write_eigenvalues(i,dk(i),kp(:,i),eival(:,i),ndn,uio)
!   write(uio,7)i,kp(:,i),(eivl(j),j=1,ndn)

 enddo kloop

! Band structure projection sorting, added by Bolin
!   allocate(mp1(ndn,nkp))
!    write(*,*) "Projection Band Sorting for band structures!"
!   call band_sort_bs(nkp,ndn,kp,dk,eival,eivec,mp1)
!   allocate(eival_tmp(ndn),eivec_tmp(ndn,ndn),vg_tmp(3,ndn))
!   do i=1,nkp
!      eival_tmp=eival(:,i)
!      eivec_tmp=eivec(:,:,i)
!      vg_tmp=vg(:,:,i)
!      do j=1,ndn
!         eival(j,i)=eival_tmp(mp1(j,i))
!         eivec(:,j,i)=eivec_tmp(:,mp1(j,i))
!         vg(:,j,i)=vg_tmp(:,mp1(j,i))
!      end do
!      call write_eigenvalues(i,dk(i),kp(:,i),eival(:,i),ndn,uio)
!   end do
!   deallocate(mp1,eival_tmp,eivec_tmp,vg_tmp)

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
  implicit none
  integer ndn,i,j
  real(8) q0(3),vgr(3,ndn),q1(3),dq,om0,om1
  real(8) evl0(ndn),evlp(ndn),evlm(ndn)
  complex(8) evc0(ndn,ndn),evct(ndn,ndn)

  dq=3d-4  ! this is the ideal dq (at least for si)

  call get_freq(q0,ndn,vgr,evl0,evc0)

  return

  do i=1,3
     q1=q0; q1(i)=q0(i)+dq
     call get_freq(q1,ndn,vgr,evlp,evct)
     q1=q0; q1(i)=q0(i)-dq
     call get_freq(q1,ndn,vgr,evlm,evct)
     vgr(i,:)=(evlp-evlm)/2/dq / 2/sqrt(abs(evl0)) *cnst*1d-8 *2*pi !*c_light
  enddo
  write(30,6)'q,dq,om=',q0,dq,evl0
  do i=1,3
     write(30,5)'i,vgr(i,l)=',i,vgr(i,:)*c_light
  enddo
  write(30,*)'*******************************************'


 4 format(a,3(1x,f9.3),a)
 5 format(a,i5,99(1x,f9.3))
 6 format(a,3(1x,f9.3),1x,g12.6,9(1x,f9.4))

  end subroutine finitedif_vel
!===========================================================
 subroutine get_freq(kp,ndn,vg,eival,eivec)
 use params
 use ios
 use constants
 use born
 implicit none
 integer, intent(in) :: ndn
 real(8), intent(in) :: kp(3)
 real(8), intent(out):: eival(ndn)
 complex(8), intent(out) :: eivec(ndn,ndn)
 integer i,j,k,l,ier,nd2,al
 integer, allocatable :: mp(:)
 real(8), allocatable :: eivl(:)
 real(8) absvec,vg(3,ndn)
 complex(8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:),temp(:,:)
 real(8) khat(3)

 nd2 = min(ndn,12)
 allocate(dynmat(ndn,ndn),mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3),temp(ndn,ndn))
! write(uio,*)'# dynmat allocated, printing: Kpoint, eival(j),j=1,ndyn)'
! write(uio,*)'############ before setting up dynmat, kp number ',i

    call set_dynamical_matrix(kp,dynmat,ndn,ddyn)

!JS: call nonanalytical term for born effective change
    if (kp(1)==0.d0 .AND. kp(2)==0.d0 .AND. kp(3)==0.d0) then
        khat=kp(:)+1.0D-10
    else
        khat=kp(:)
    endif
    call nonanal(khat,dynmat,ndn,ddyn)

    if (verbose) then
       write(ulog,3)' ======================================================================'
       write(ulog,3)' THE DYNAMICAL MATRIX for KP=',ndn,kp(:)
       do l=1,ndn
          write(ulog,8)(dynmat(l,j),j=1,nd2)
       enddo
       write(ulog,3)' ======================================================================'
       do al=1,3
          write(ulog,3)' THE k-DERIVATIVE OF THE DYNAMICAL MATRIX =',al
          do l=1,ndn
             write(ulog,8)(ddyn(l,j,al),j=1,nd2)
          enddo
       enddo
    endif

! store dynmat in temp
    temp=dynmat
    call diagonalize(ndn,dynmat,eivl,ndn,eivc,ier)

! sort eivals in ascending order
    call sort(ndn,eivl,mp,ndn)
!   write(ulog,6)'map=',mp

    if (ier.ne.0 .or. verbose) then
      write(   *,*)' ier=',ier
      write(ulog,*)' ier=',ier, ' EIGENVALUES ARE...'
      write(ulog,4) eivl
    endif

    do j=1,ndn
       eival(j) = eivl(mp(j))  ! so that all frequencies are positive
       do l=1,ndn
          eivec(l,j) = eivc(l,mp(j))
       enddo
    enddo

! if pure imaginary make eigenvectors real as they can be multiplied by any arbitrary number
    do l=1,ndn
       absvec = sum(abs(real(eivec(:,l))))
       if (absvec .lt. 1d-3) then
          eivec(:,l)=cmplx(0,1)*eivec(:,l)
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
     dynmat=temp
     temp = matmul(dynmat,eivc)
     temp = matmul(transpose(conjg(eivc)),temp)
     write(ulog,*)' TEST: below square of eivals should appear in the diagonal if e.D.e is right'
     do j=1,ndn
       write(ulog,9)'e.D.e=',j, temp(j,:)
     enddo
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
      vg(al,l)=vg(al,l)/2/sqrt(abs(eival(l)))*cnst*1d-10*100*2*pi
    enddo
    enddo

 deallocate(eivl,eivc,dynmat,mp,ddyn,temp)

 3 format(a,i5,9(1x,f19.9))
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
! use force_constants_module
 use svd_stuff
 use eigen
 use constants
 implicit none
 integer ik,i0,nkp,la,al,be,ga,j,k,j0,k0,ta1,ta2,t,ugr2,ndn    ! ,i3,j3,k3
 real(8) mi,mj,rr3(3),rr2(3),qq(3),denom,qdotr,omk,mysqrt
 real(8) kp(3,nkp),dk(nkp),eivl(ndn,nkp)
 complex(8) zz,one,term
 complex(8) grn(ndn,nkp), eivc(ndn,ndn,nkp)

 one = cmplx(1d0,0d0)
 write(ugr2,*)'# la,nk,dk(nk),kp(:,nk),om(la,nk),gruneisen(la,nk))'
 do ik=1,nkp
    qq(:) = kp(:,ik)
 do la=1,ndn
!   write(ulog,6)'Starting gruneisen parameters for kpoint & band# ',ik,la,kp(:,ik)
!   write(ulog,*)' i,la,t,fcs_3(t),term(2),6*eival(la,i),rr3(3),grun(la,i)'
    grn(la,ik) = 0
    denom = 6 * eivl(la,ik)
    omk = mysqrt(eivl(la,ik))*cnst
    do i0=1,natom_prim_cell
         mi = atom0(i0)%mass
       tloop: do t=1,nterms(3)

         if ( i0 .ne. iatomterm_3(1,t) ) cycle tloop
         al = ixyzterm_3(1,t)
         be = ixyzterm_3(2,t)
         ga = ixyzterm_3(3,t)
!        i0 = iatomcell0(iatomterm_3(1,t))
         j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ; mj = atom0(j0)%mass
         k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k)
         ta1= al + 3*(i0-1)
         ta2= be + 3*(j0-1)

! rr2(:)=iatomcell(1,j)*r1+iatomcell(2,j)*r2+iatomcell(3,j)*r3
!  be careful: it has to be the translations R not atompos!
         rr2(:) = atompos(:,j) - atompos(:,j0)  ! R
         rr3(:) = atompos(:,k)                  ! R+tau
         qdotr =  ( qq .dot. rr2)
         zz = cdexp( ci * qdotr )
!! term = - ampterm_3(t)*fcs_3(igroup_3(t))*zz*eivc(ta1,la,i)*conjg(eivc(ta2,la,i))*rr3(ga)/sqrt(mi*mj)
         term = - fcs_3(t) * zz * eivc(ta2,la,ik)*conjg(eivc(ta1,la,ik))*rr3(ga)/sqrt(mi*mj)
         grn(la,ik) = grn(la,ik) + term
!        write(ulog,7)i,la,t,fcs_3(t),rr2,qdotr,zz,rr3,grn(la,ik)
       enddo tloop
    enddo
    grn(la,ik) = grn(la,ik)/denom
    if (aimag(grn(la,ik)) .gt. 1d-4) then
       write(ulog,*)' GRUNEISEN: ... has large imaginary part!! for nk,la=',ik,la,grn(la,ik)
!      stop
    endif
    write(ugr2,6)' ',la,ik,dk(ik),kp(:,ik),omk,real(grn(la,ik))
 enddo
!   write(ugr,8)ik,dk(ik),kp(:,ik),(real(grn(la,ik)),la=1,ndn)
 enddo
5 format(4i7,9(1x,g10.4))
6 format(a,2i5,99(1x,g10.4))
7 format(i5,i5,i6,99(1x,g10.4))
8 format(i8,99(1x,g10.4))
! deallocate(eivl,eivc,kg,grn)
 end subroutine gruneisen
!============================================================
 subroutine mechanical_old(bulk,c11,c44,dlogv)
 use ios
 use lattice
 use geometry
 use atoms_force_constants
! use force_constants_module
 use svd_stuff
 use constants
 implicit none
 integer i0,al,be,j,t
 real(8) bulk,rija,rijb,c11,c44,dlogv

! write(ulog,*)'BULK_MOD: i0, al, j, be, rija,rijb,fcs_2(t),bulk'
 bulk=0; c11=0 ; c44=0
 do i0=1,natom_prim_cell
 do t=1,nterms(2)
    if ( i0 .ne. iatomterm_2(1,t) ) cycle
    al =  ixyzterm_2(1,t)
    be =  ixyzterm_2(2,t)
    j  = iatomterm_2(2,t)
    rija = atompos(al,j)-atompos(al,i0)
    rijb = atompos(be,j)-atompos(be,i0)
!   write(ulog,3)i0, al, j, be, rija,rijb,fcs_2(t),bulk
    bulk = bulk - fcs_2(t)*(1-2*grun_fc(t)*dlogv)*rija*rijb
    if (al.eq.1 .and. be.eq.1) then
       c11 = c11 - fcs_2(t)*(1-2*grun_fc(t)*dlogv)*rija*rijb
    endif
 enddo
 enddo
 bulk = bulk/volume_r/18
 c11 = c11/volume_r/2

 write(ulog,*)'BULK_MODULUS: in eV/A^3 ',bulk
 bulk = bulk*1d30*ee
 write(ulog,*)'BULK_MODULUS: in SI, Mbar ',bulk,bulk*1d-11
 c11 = c11*1d30*ee
 write(ulog,*)'C11 in SI, Mbar ',c11,c11*1d-11
3 format(4(i4),9(2x,f10.5))

 end subroutine mechanical_old
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
    cv = cv/nat*n_avog*k_b  ! this is per mole    )/(volume_r*1d-30)  !
    etot = etot/nat*n_avog   ! convert from Joule/cell to Joule per mole
    free = free/nat*n_avog   ! convert from Joule/cell to Joule per mole
    pres0= pres/(volume_r*1d-30) ! * 1d-8  ! 1d-8 is to convert to kbar
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
    call get_k_info(-kp(:,i),NC,mi,i1,j1,k1,g01,g02,g03,inside)
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

   call get_k_info(kp(:,i)-shift,NC,l,i3,j3,k3,g01,g02,g03,inside)
   if (i.ne.l) then
      write(ulog,*) 'i .ne. l ',i,l
      stop
   endif
!  call bring_to_cell_c(-kp(:,i),g1,g2,g3,r1,r2,r3,w)
!  w =w /2/pi
   w = -kp(:,i)+shift
   call get_k_info(w,NC,j,i2,j2,k2,g01,g02,g03,inside)
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
 use constants
 use lattice
 use atoms_force_constants
 use born
 use geometry

 integer, intent(in) ::ndim
 complex(8), intent(inout) :: dynmat(ndim,ndim),ddyn(ndim,ndim,3)
 real(8), intent(in) :: q(3)
 real(8) zag(3),zbg(3),qeq,ma,mb,rr(3),dqeq(3)
 integer na,nb,i,j,al
 real(8) eps_scale,q2,gg,term

 eps_scale=8.8541878176D-12/1D10/ee
 gg=2*pi/(volume_r0)**0.33
! if ( gg .lt. 4*rho) then
 if (born_flag.eq.0)  rho = gg/4d0
! write(30,*)'ADOPTED VALUE of RHO = ',rho
! endif
 rho2=rho*rho


 q2=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
 qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+ &
&       q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+ &
&       q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
 dqeq(1) = 2*epsil(1,1)*q(1)+(epsil(1,2)+epsil(2,1))*q(2)+(epsil(1,3)+epsil(3,1))*q(3)
 dqeq(2) = 2*epsil(2,2)*q(2)+(epsil(2,3)+epsil(3,2))*q(3)+(epsil(2,1)+epsil(1,2))*q(1)
 dqeq(3) = 2*epsil(3,3)*q(3)+(epsil(3,1)+epsil(1,3))*q(1)+(epsil(3,2)+epsil(2,3))*q(2)

 term=1d0/(qeq*eps_scale)/volume_r0*exp(-q2/rho2)

 do na = 1,natom_prim_cell
   ma = atom0(na)%mass
 do nb = 1,natom_prim_cell
   mb = atom0(nb)%mass
   rr = atompos(:,na)-atompos(:,nb)
   do i=1,3
      zag(i) = q(1)*zeu(1,i,na)+q(2)*zeu(2,i,na)+q(3)*zeu(3,i,na)
      zbg(i) = q(1)*zeu(1,i,nb)+q(2)*zeu(2,i,nb)+q(3)*zeu(3,i,nb)
   end do
   do i = 1,3
   do j = 1,3
!write(*,*)
!dynmat(i+3*(na-1),j+3*(nb-1)),
!4*pi*zag(i)*zbg(j)/(qeq*eps_scale)/om0/sqrt(ma*mb)*exp(-q2/rho2),ma,mb,zag(i),zbg(j)
!4*pi*zag(i)*zbg(j)/(qeq*eps_scale)/om0/sqrt(ma*mb)*exp(-q2/rho2)
      dynmat(i+3*(na-1),j+3*(nb-1)) = dynmat(i+3*(na-1),j+3*(nb-1))+ &
      & zag(i)*zbg(j)*term/sqrt(ma*mb)
  !! ALSO NEED TO CORRECT THE GROUP VELOCITIES
     do al=1,3
        ddyn(i+3*(na-1),j+3*(nb-1),al) = ddyn(i+3*(na-1),j+3*(nb-1),al)+term/sqrt(ma*mb) *  &
 &          (zeu(al,i,na)*zbg(j)+zeu(al,j,nb)*zag(i)-dqeq(al)*zag(i)*zbg(j)/qeq)
     enddo
   end do
   end do
 end do
 end do

 end subroutine nonanal
!=========================================================
 subroutine threeph_phase_space(ndn,ni,ki,nk,kp,eival,eivec)
! this subroutine computes for each q the sum_q2 delta(w-w_q2 \pm w_q3)
! and finds the kpoints available
 use kpoints
 use ios
 use lattice
 use geometry
 use atoms_force_constants
! use force_constants_module
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ni,ndn,l1,l2,l3,n2,n3,j,k
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk3,indx,igam
 real(8) q1(3),q2(3),q3(3),eival(ndn,nk)
 real(8) ki(3,ni),kp(3,nk)
 complex(8) eivec(ndn,ndn,nk),xx

 indx=0
 write(ulog,*)'V3: nc(123)=',nc

    q1=0
    call get_kindex(q1,nk,kp,igam)
    q1 = kp(:,igam+3)

 loop3: do n3=1 ,nk
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    call get_k_info(q3,NC,nk3,i3,j3,k3,g01,g02,g03,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif

!! properly initilize q2 !!!!!!!!
    q1 = -q2-q3

    call get_k_info(q1,NC,nk1,i1,j1,k1,g01,g02,g03,inside)

  write(*,2)'nk13=',nk1,nk3

 enddo loop3
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)


2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))

 end subroutine threeph_phase_space
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
subroutine band_sort_bs(nkp,ndyn,kp,dk,eival,eivec,map)
implicit none

integer, intent(in) :: nkp,ndyn
real(8), intent(in) :: kp(3,nkp)
real(8), intent(in) :: eival(ndyn,nkp)
complex(8), intent(in) :: eivec(ndyn,ndyn,nkp)
integer, intent(out) :: map(ndyn,nkp)

integer i,j,k,l
integer ntbd !Number of connections to be determined

real(8) bmax !Maximum band derivative max(dE/dk)
real(8) dk(nkp)
real(8) fsaf !Safety factor

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
    map(i,1)=i
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
                if ( map(l,i-1) .eq. j) then
                    map(l,i)=k
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
                     if ( map(l,i-1) .eq. j) then
                        map(l,i)=k
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
              if ( map(l,i-1) .eq. j) then
                 map(l,i)=k
              end if
           end do
        ntbd=ntbd-1
        overlap_q(j,:)=0
        overlap_q(:,k)=0
        end if
    end do
    end do
    end if

        write(*,*) 'map for',i
        write(*,*) map(:,i)

    if (ntbd .ne. 0) then
        write(*,*) 'error: band sorting failed'
        write(*,*) 'k-point:',i,kp(:,i),ntbd
    end if
end do

deallocate(overlap,overlap_q)

end subroutine band_sort_bs
!==========================================================
 subroutine set_dynamical_matrix_old(kpt,dynmat,ndim,ddyn)
 use params
 use lattice
 use kpoints
 use atoms_force_constants
 use geometry
 use svd_stuff
 implicit none
 integer, intent(in) :: ndim
 complex(8), intent(out) :: dynmat(ndim,ndim),ddyn(ndim,ndim,3)
 real(8), intent(in) :: kpt(3)
 complex(8) junk
 real(8) mi,mj,all,rr(3),delt(3)
 integer i0,j,j0,al,be,i3,j3,t !,ired

4 format(a,4(1x,i5),9(2x,f9.4))
9 format(i6,99(1x,f9.3))

 delt = wshift !*wshift
 ddyn   = cmplx(0d0,0d0)
 dynmat = cmplx(0d0,0d0)
 do i0=1,natom_prim_cell
 do al=1,3
    i3 = al+3*(i0-1)
    mi = atom0(i0)%mass
! write(ulog,*) 'i,al,mass=',i0,al,i3,mi
    tloop: do t=1,nterms(2)
       if ( i0 .eq. iatomterm_2(1,t) .and. al .eq. ixyzterm_2(1,t) ) then
          be =  ixyzterm_2(2,t)
          j  = iatomterm_2(2,t)
          j0 = iatomcell0(j)
          mj = atom0(j0)%mass
          j3 = be+3*(j0-1)
          rr = atompos(:,j)-atompos(:,j0)
!          ired = igroup_2(t)
!          junk = fcs_2(ired)*ampterm_2(t)* exp(ci*(kpt .dot. rr))/sqrt(mi*mj)
          junk = fcs_2(t)* exp(ci*(kpt .dot. rr))/sqrt(mi*mj)
!  write(ulog,4) 't,j,be,mass,dyn=',t,j,be,j3,mj,junk
          dynmat(i3,j3) = dynmat(i3,j3) + junk
          ddyn(i3,j3,:) = ddyn(i3,j3,:) + junk*ci*rr(:)
!         if (be.eq.al .and. j0.eq.i0 ) then
!         if (be.eq.al .and. (length(rr-r1) .myeq. 0d0)) then
! same type and same cartesian coord but not the same atom, just take the images
!     write(ulog,5)'i0,j0,j,i3,j3,rij=',i0,j0,j,i3,j3,rr
!            dynmat(i3,j3) = dynmat(i3,j3) - delt/2
!         endif
       endif
    enddo tloop
! this leaves  gamma eivals unchanged
!   dynmat(i3,i3) = dynmat(i3,i3) + delt(al)*(1-cos(kpt .dot. r1))
! but this one shifts everything up by delt=wshift
    dynmat(i3,i3) = dynmat(i3,i3) + delt(al)
 enddo
 enddo

 if (verbose) then
  write(ulog,*)'SET_DYNAMICAL_MATRIX: d(dyn)/dk is:'
  do t=1,3
     write(ulog,*)'=======component of v ',t
     do j=1,ndim
        write(ulog,9) j, ddyn(j,:,t)
     enddo
  enddo
 endif

5 format(a,5i5,9(f8.3))
 all = sum(cdabs(dynmat(:,:)))/(ndim*ndim)
! make sure it is hermitian
 do t=1,ndim
    if (abs(aimag(dynmat(t,t))) .gt. 9d-4*abs(real(dynmat(t,t))) ) then
       write(ulog,*)' dynmat is not hermitian on its diagonal'
       write(ulog,*)' diagonal element i=',t,dynmat(t,t)
!      stop
!   else
       dynmat(t,t) = cmplx(real(dynmat(t,t)),0d0)
    endif
  do j=t+1,ndim-1
    if (abs(aimag(dynmat(t,j))+aimag(dynmat(j,t))) .gt. 9d-4*all ) then
       write(ulog,*)' dynmat is not hermitian in AIMAG of its off-diagonal elts'
       write(ulog,*)' off-diagonal element i,j=',t,j,dynmat(t,j),dynmat(j,t)
       write(ulog,*)' comparing to avg(abs(dynmat))=',all
!      stop
    elseif(abs(real(dynmat(t,j))-real(dynmat(j,t))) .gt. 9d-4*all ) then
       write(ulog,*)' dynmat is not hermitian in REAL of its off-diagonal elts'
       write(ulog,*)' off-diagonal element i,j=',t,j,dynmat(t,j),dynmat(j,t)
       write(ulog,*)' comparing to avg(abs(dynmat))=',all
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
!  do t=1,ndim
!  do j=1,ndim
!     dynmat(t,j)=(dynmat(t,j)+conjg(dynmat(t,j)))/2d0
!  enddo
!  enddo
! endif

 end subroutine set_dynamical_matrix_old
!==========================================================
 subroutine diagonalize(n,mat,eival,nv,eivec,ier)
! n=size of mat; nv is the number of needed eigenvectors
 implicit none
 integer n,ier,zero,nv  !,i,j
 complex(8) mat(n,n),eivec(n,n)
 real(8)  eival(n),tol
! This is used by eigch
! real(8), allocatable :: w(:,:)
! integer, allocatable :: lw(:)
! This is used by ZHEGV
  real(8), allocatable :: rwork(:)
  complex(8), allocatable :: work(:)
  integer lwork


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
 subroutine write_eigenvalues(i,dk,kp,eival,n,uio)
 use constants
 implicit none
 integer i,j,n,uio
 integer sig(n)
 real(8), dimension(n) :: eival
 real(8) kp(3),dk,mysqrt

! n = size(eival)
! write(ueiv,3)i-1,kp,(eival(j),j=1,n)
! do j=1,n
!    sig(j) = 1
!    if (eival(j).lt.0) then
!       eival(j) = -eival(j)
!       sig(j) = -1d0
!    endif
! enddo
! write(ueiv,3)i-1,kp,(cnst*sig(j)*sqrt(eival(j)),j=1,n)
 write(uio,3)i-1,dk,kp,(cnst*mysqrt(eival(j)),j=1,n)

3 format(i5,1x,4(1x,f8.3),999(1x,g11.5))
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
 complex(8) oc2,omz

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

          iself = aimag(oc2(tmp,omz,j,k))/pi*wkp(j)*wkp(k)*tmp
       enddo
       enddo
       write(udosj,3) i,2*omega(i),dsm,dsp,iself
    enddo
    close(udosj)

2 format(9(3x,g11.5))
3 format(i5,9(3x,g11.5))

 end subroutine calculate_jdos
!==========================================================
 subroutine lifetime_dos(q,la,omq,dosp,dosm)
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
 integer i,j,k,l,la,nq,i1,j1,k1,nq1,ip,jp,kp,im,jm,km,nqp,nqm,la1,la2,inside
 real(8) q(3),q1(3),qp(3),qm(3),dosp,dosm,om1,omp,omm,omq,delta_l
 complex(8) oc2,omz

 dosp=0; dosm=0
 do l=1,nkc
    q1=kpc(:,l)
    call get_k_info_cent(q1-shift,NC,nq1,i1,j1,k1,g01,g02,g03,inside)
    if(l.ne.nq1) then
      write(ulog,*)'n1,nq1,inside=',l,nq1,inside
      write(ulog,*)'q1=',q1
      write(ulog,*)'ijk1=',i1,j1,k1
      stop
    endif
    qp=-q-q1
    qm=-q+q1
    call get_k_info_cent(qp-shift,NC,nqp,ip,jp,kp,g01,g02,g03,inside)
    call get_k_info_cent(qm-shift,NC,nqm,im,jm,km,g01,g02,g03,inside)
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
 subroutine phase_space_lifetime(i1,i2,nk,kp,eival)
 use ios
 use om_dos
 use kpoints
 use eigen
 use constants
 implicit none
 integer i,la,nk,i1,i2,udos
 real(8) q(3),omq,dosp,dosm,kp(3,nk),eival(ndyn,nk)  !,vq(3,ndyn)
 complex(8) eivc(ndyn,ndyn)
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
          call lifetime_dos(q,la,omq,dosp,dosm)
          write(udos,3) la,q,omq,dosm/nkc,dosp/nkc
       enddo
    enddo
 close(udos)

3 format(i5,9(3x,g10.4))

 end subroutine phase_space_lifetime
!===========================================================
