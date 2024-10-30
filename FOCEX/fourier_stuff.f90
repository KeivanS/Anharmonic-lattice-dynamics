!============================================================

 module fourier
 use constants, only : r15
  integer nrgrid,nggrid,nsubgrid,nreg,n3reg(3)
  real(r15), allocatable :: rgrid(:,:),rws_weights(:),rgridreg(:,:),regweights(:)
  real(r15), allocatable :: ggrid(:,:),gws_weights(:),ggridreg(:,:)
  real(r15), allocatable :: subgrid(:,:),subgrid_weights(:)
  integer, allocatable :: i1r(:),i2r(:),i3r(:),i1g(:),i2g(:),i3g(:)

  interface fourier_r2k
     module procedure fourier_r2k_r,fourier_r2k_c,fr2k_5
  end interface
  interface fourier_k2r
     module procedure fourier_k2r_r,fourier_k2r_c,fk2r_5
  end interface

  contains
!-------------------------------------------
 subroutine fr2k_5(nr,rgrid2,fr,wei,ng,ggrid2,fk)
!! for an array of size nk defined on a mesh rmesh(3,nr), this subroutine
!! calculates its cosine Fourier transform on the kmesh, conjugate to rmesh
! a cosine transform is enough for even functions, otherwise R,-R must be symmetrized
! from r to k all weights are =1

 use geometry
 implicit none
 integer, intent(in) :: nr,ng
 real(r15), intent(in) :: fr(nr),wei(nr),rgrid2(3,nr),ggrid2(3,ng)
 real(r15), intent(out) :: fk(ng)
 integer ik,ir

 fk=0
 do ik=1,size(fk) !ng
    do ir=1,size(fr) !nr
       fk(ik)=fk(ik)+cos(ggrid2(:,ik).dot.rgrid2(:,ir)) * fr(ir) * wei(ir)
    enddo
 enddo

 end subroutine fr2k_5
!-------------------------------------------
 subroutine fk2r_5(nk,ggrid2,fk,wei,nr,rgrid2,fr)
!! for an array of size nk defined on a mesh rmesh(3,nr), this subroutine
!! calculates its cosine Fourier transform on the kmesh, conjugate to rmesh
! a cosine transform is enough for even functions, otherwise R,-R must be symmetrized
! from r to k all weights are =1

 use geometry
 implicit none
 integer, intent(in) :: nr,nk
 real(r15), intent(in) :: fk(nk),wei(nk),rgrid2(3,nr),ggrid2(3,nk)
 real(r15), intent(out) :: fr(nr)
 integer ik,ir

 fr=0
 do ir=1,size(fr) !nr
    do ik=1,size(fk) !nk
       fr(ir)=fr(ir)+cos(ggrid2(:,ik).dot.rgrid2(:,ir)) * fk(ik) * wei(ik)
    enddo
 enddo

 end subroutine fk2r_5
!-------------------------------------------
 subroutine fourier_r2k_r(fr,fk)
!! for an array of size nk defined on a mesh rmesh(3,nr), this subroutine
!! calculates its cosine Fourier transform on the kmesh, conjugate to rmesh
! a cosine transform is enough for even functions, otherwise R,-R must be symmetrized
! from r to k all weights are =1

 use geometry
 implicit none
 real(r15), intent(in) :: fr(:) !nrgrid)
 real(r15), intent(out) :: fk(:) !nggrid)
 integer ik,ir

 fk=0
 do ik=1,size(fk) !nggrid
    do ir=1,size(fr) !nrgrid
       fk(ik)=fk(ik)+cos(ggrid(:,ik).dot.rgrid(:,ir)) * fr(ir) * rws_weights(ir)
    enddo
 enddo

 end subroutine fourier_r2k_r
!-------------------------------------------
 subroutine fourier_k2r_r(fk,fr)
!! for an array of size nk defined on a mesh kmesh(3,nk), this subroutine
!! calculates its cosine Fourier transform on the rmesh, conjugate to kmesh
! there is k,-k symmetry from time reversal, so a cosine transform is enough
 use geometry
 implicit none
 real(r15), intent(in) ::  fk(:) !nggrid)
 real(r15), intent(out) :: fr(:) !nrgrid)
 integer ik,ir

 fr=0
 do ir=1,size(fr) !nrgrid
    do ik=1,size(fk) !nggrid
       fr(ir)=fr(ir)+cos(ggrid(:,ik).dot.rgrid(:,ir)) * fk(ik) * gws_weights(ik)
    enddo
 enddo

 end subroutine fourier_k2r_r
!-------------------------------------------
 subroutine fourier_r2k_c(fr,fk)
!! for an array of size nk defined on a mesh rmesh(3,nr), this subroutine
!! calculates its complex Fourier transform on the kgrid, conjugate to rgrid
! a cosine transform is enough for even functions, otherwise R,-R must be symmetrized
! from r to k all weights are =1

 use geometry
 use constants, only : ci
 implicit none
 complex(r15), intent(in) :: fr(:) !nrgrid)
 complex(r15), intent(out) :: fk(:) !nggrid)
 integer ik,ir

 fk=0
 do ik=1,size(fk) !nggrid
    do ir=1,size(fr) !nrgrid
       fk(ik)=fk(ik)+cdexp(ci*(ggrid(:,ik).dot.rgrid(:,ir))) * fr(ir) * rws_weights(ir)
    enddo
 enddo

 end subroutine fourier_r2k_c
!-------------------------------------------
 subroutine fourier_k2r_c(fk,fr)
!! for an array of size nk defined on a mesh kmesh(3,nk), this subroutine
!! calculates its cosine Fourier transform on the rmesh, conjugate to kmesh
! there is k,-k symmetry from time reversal, so a cosine transform is enough
 use geometry
 use constants, only : ci
 implicit none
 complex(r15), intent(in) ::  fk(:) !nggrid)
 complex(r15), intent(out) :: fr(:) !nrgrid)
 integer ik,ir

 fr=0
 do ir=1,size(fr) !nrgrid
    do ik=1,size(fk) !nggrid
       fr(ir)=fr(ir)+cdexp(-ci*(ggrid(:,ik).dot.rgrid(:,ir))) * fk(ik) * gws_weights(ik)
    enddo
 enddo

 end subroutine fourier_k2r_c
!-------------------------------------------
 subroutine extract_fourier(ncfg,nat,dsp,frc,nrmesh,rmesh,nkmesh,kmesh,maprtausc)
!! this subroutine takes force-displacements data X(r,tau,s), Fourier transforms to get f(k,tau,s)
!! and u(k,tau,s) for tau in the primitive cell, and s=snapshot#, and a supercell translation vector
!! k (+ its stars) .
!! Then extracts D(k;tau,tau') defined by F(k,tau,s)=-D(k,tau,tau')*u(k,tau',s) using data from all
!! snapshots "s" by direct SVD inversion and enforcing group symmetry relations on D(k) and its stars
!! then calculates phonon dipersion by diagonalization, and also calculates phi(tau,tau'+R) by inverse
!! Fourier transformation.
!! output is the force constants phi (tau,R+tau') where R belongs to the R mesh (primitive lattice)
 use atoms_force_constants !, only : natom_prim_cell,op_matrix,op_kmatrix
 use svd_stuff
 use lattice, only : primitivelattice,cart_to_prim,prim_to_cart,red2cart_g,cart2red_g
 use kpoints, only : kibz,nibz !,wibz
 use ios , only : ulog,umatrx ,ufc2,utimes,write_out
 use geometry
 use params, only : verbose
 use atoms_force_constants
 implicit none
 integer, intent(in) :: ncfg,nat,nrmesh,nkmesh,maprtausc(natom_prim_cell,nrmesh)
 real(r15), intent(in) :: dsp(3,nat,ncfg),frc(3,nat,ncfg)
 real(r15), intent(in) :: rmesh(3,nrmesh),kmesh(3,nkmesh)
 real(r15), allocatable:: ur(:,:),uk(:,:),fr(:,:),fk(:,:),amatk(:,:),bmatk(:),coefs(:,:),sigm(:)
 complex(r15), allocatable:: dynk(:,:,:) ,aux(:) ,d2(:,:),phi(:,:,:)
 integer nat0,icfg,ik,iarm,l,i,j,tau,ir,iatom,nu,nkstar,dim_l,dimdyn, &
 &       kvecop(48),narms,nsym,dim_c,n2,iop
 real(r15) error,ermax,sig,kvecstar(3,48),sk(3),symk(3,3)
 real tim

! nat0=nat/nrmesh
! dimal=3*nat0*ncfg
 nat0=natom_prim_cell
 dimdyn=(3*nat0*(3*nat0+1))/2    ! size of 1D dynamical matrix for each k
 write(*,*)'EXTRACT_FOURIER: nat0,rmesh=',nat0,rmesh
 allocate(ur(3*nat0,nrmesh),uk(3*nat0,nrmesh),fr(3*nat0,nrmesh),fk(3*nat0,nrmesh))
 allocate(dynk(3*nat0,3*nat0,nkmesh),phi(3*nat0,3*nat0,nkmesh),d2(3*nat0,3*nat0))
! &        amatk(dimal,dimac,nkmesh),bmatk(dimal,nkmesh))

! first need to group the kmesh vectors into arms of vectors within IBZ
!    call get_weights4(nkmesh,kmesh,ngibz,mapibz,gibz,wgibz)

  ibzloop: do ik=1,nibz

!    call getstar(kibz(:,ik),primitivelattice,narms,kvecstar,kvecop)
! getstar acts on reduced vectors
     call getstar(cart2red_g(kibz(:,ik)),primitivelattice,narms,kvecstar,kvecop)

! solve for all the d(k) and its stars together
     dim_c=dimdyn*narms
     nsym=dimdyn*(narms-1)  ! nunmber of symmetry relations
     dim_l=nsym + ncfg*3*nat0*narms
     allocate( coefs(dimdyn,dimdyn) )
     allocate( amatk(dim_l,dim_c),bmatk(dim_l))
     amatk=0; bmatk=0;

! symmetry part of the amatk *************************************************

     do iarm=1,narms  ! kvecop(iarm) is the label of the symmetry matrix

        iop=kvecop(iarm)
        symk= matmul(transpose(cart_to_prim),matmul(  &
&          op_kmatrix(:,:,iop),transpose(prim_to_cart)))
! first the symmetry relations between D(k) and D(S(k)): D(S(k))=D(k)SS
!       sk=matmul(op_kmatrix(:,:,iop),kibz(:,ik))   ! star of kibz
        sk=matmul(symk,kibz(:,ik))   ! star of kibz
!       if(length(sk-kvecstar(:,iarm)) .gt. 1d-4) then
        if(length(sk-red2cart_g(kvecstar(:,iarm))) .gt. 1d-4) then
            write(*,2)'PROBLEM IN SYMMETRY OPERATIONS for arm# and kibz=',iarm,kibz(:,ik)
            write(*,2)'symmetry operation # ',iop
            write(*,3)symk(1,:)
            write(*,3)symk(2,:)
            write(*,3)symk(3,:)
            write(*,3)'sk,kvecstar=',sk,kvecstar(:,iarm),red2cart_g(kvecstar(:,iarm))
            stop
         else
            write(*,1)'for arm=',iarm,' symmetry matrix # ',kvecop(iarm),' was used'
         endif

! assuming iarm=1 corresponds to kibz itself (iop_kmatrix=I)
         if(iarm.eq.1) cycle

! amatk((iarm-2)*dimdyn+1:(iarm-1)*dimdyn,(iarm-1)*dimdyn+1:iarm*dimdyn) = - Identity
         coefs=0
         do l=1,dimdyn
            coefs(l,l)= -1d0
         enddo
         amatk((iarm-2)*dimdyn+1:(iarm-1)*dimdyn,(iarm-1)*dimdyn+1:iarm*dimdyn)=coefs

! amatk((iarm-2)*dimdyn+1:(iarm-1)*dimdyn,1:dimdyn)  obtained by symmetry
!         call get_relation_dynstar(nat0,kibz(:,ik),iarm,kvecstar,iop,op_kmatrix(:,:,iop),dimdyn,coefs)
!        call get_relation_dynstar(nat0,iop,op_kmatrix(:,:,iop),dimdyn,coefs)
         call get_relation_dynstar(nat0,iop,symk,dimdyn,coefs)
         amatk((iarm-2)*dimdyn+1:(iarm-1)*dimdyn,1:dimdyn)=coefs(1:dimdyn,1:dimdyn)
         bmatk((iarm-2)*dimdyn+1:(iarm-1)*dimdyn)=0d0

     enddo

! force-displacement part of the amatk**********************

     armloop: do iarm=1,narms  ! loop over columns

! find nkstar, the index of the Kvecstar(iarm) in the kmesh
!       call find_in_mesh(kvecstar(:,iarm),nkmesh,kmesh,nkstar)
        call find_in_mesh(red2cart_g(kvecstar(:,iarm)),nkmesh,kmesh,nkstar)

        snaploop: do icfg=1,ncfg

! get the force vector -> bmatk=-F(nu,kstar,icfg)
! and displacement vector -> amatk(icfg,kstar block)=u(nu',kstar,icfg) for every kstar
! store for that snapshot the supercell data in an array and then fourier transform
           do tau=1,nat0      ! tau index
           do j=1,3
              nu=j+3*(tau-1)  ! store fk(j,tau) and uk(j,tau) in one array indexed by nu
              do ir=1,nrmesh  ! R index
                 iatom=maprtausc(tau,ir)  ! atom index in supercell of (tau,ir)
                 fr(nu,ir)=frc(j,iatom ,icfg)
                 ur(nu,ir)=dsp(j,iatom ,icfg)
              enddo

              call fourier_r2k(fr(nu,:),fk(nu,:))
              call fourier_r2k(ur(nu,:),uk(nu,:))

           enddo
           enddo

           n2=nsym+3*nat0*narms*(icfg-1)+(iarm-1)*3*nat0
           bmatk(n2+1:n2+3*nat0)=-fk(:,nkstar)

! for 1D storage of upper triangluar part of dynmat: do i=1,3N; do j=i,3N aux(l)=D(i,j)
! and for given nu=(i-1)*3N+j, j=mod(nu,3N) and i=(nu-j)/3N + 1 and d(i,j)uk(j)=-fk(i)
           do nu=1,dimdyn
              j=mod(nu,3*nat0)
              i=(nu-j)/(3*nat0)+1
              amatk(n2+i , (iarm-1)*dimdyn+nu )=uk(i,nkstar)
           enddo

           if(verbose) then
             write(umatrx,*)'***** Configuration # ',icfg,' for kibz=',ik
             do j=1,3*nat0
                write(umatrx,3)amatk(n2+j,:),bmatk(n2+j)
             enddo
           endif
        enddo snaploop
     enddo armloop

! now perform the svd for that kibz
    allocate(aux(dim_c),sigm(dim_c))
    call svd_set(dim_l,dim_c,amatk,bmatk,aux,sigm,svdcut,error,ermax,sig,'svd_k.dat')

    write(ulog,5)'FOR ik,kibz=',ik,kibz(:,ik)
    write(ulog,*)'After svd, || F_dft-F_fit || / || F_dft || (K)=',sig

  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after SVD in EXTRACT_FOURIER                IS ',tim

! transform aux to dynamical matrix for later Fourier transformation to get FCs
    do iarm=1,narms
! find the index in kmesh of all these kstars
!      call find_in_mesh(kvecstar(:,iarm),nkmesh,kmesh,nkstar)
       call find_in_mesh(red2cart_g(kvecstar(:,iarm)),nkmesh,kmesh,nkstar)
       do nu=1,dimdyn
          j=mod(nu,3*nat0)
          i=(nu-j)/(3*nat0)+1
          dynk(i,j,iarm)=     aux((iarm-1)*dimdyn+nu)
          dynk(j,i,iarm)=conjg(aux((iarm-1)*dimdyn+nu))
       enddo

       if(iarm.eq.1) then
          d2=dynk(:,:,1)
       else
! make sure dynk(iarm) for iarm.ge.2 can be obtained from dynk(1)=d2 by proper rotation
   !      s3=op_kmatrix(:,:,kvecop(iarm))
!         call rotate(d2,s3,d2rotated)
   !      do i=1,3*nat0
   !      do j=1,3*nat0
   !         if(iatomop(i,j).ne.0) then
!            if(cdabs(dynk(i,j,iarm)-d2rotated(i,j)) .gt. 1d-4 ) then
!               write(*,6)'i,j,dynk,d2rot=',i,j,dynk(i,j,iarm),d2rotated(i,j)
   !         endif
   !      enddo
   !      enddo
       endif

    enddo

  enddo ibzloop

!************************************************
!! make sure wk and dynk are properly initialized and w defined and allocated
!************************************************

! k2r: get force constants in supercell from dynamical matrix on kmesh
! make sure k and rmesh are nsc/n0 or boundary atoms repeat...
  do i=1,3*nat0
  do j=1,3*nat0
     call fourier_k2r_c(dynk(i,j,:),phi(i,j,:))
  enddo
  enddo

! write the results
  do ik=1,nrmesh
     write(ufc2,'(a,i3,a,3(1x,f9.3))')'EXTRACT_FOURIER: i=',i,'R(i)=',rmesh(:,ik)
     do i=1,3*nat0
        call write_out(ufc2,'  ',phi(i,:,ik))
     enddo
  enddo



1 format(a,i4,a,i4,a)
2 format(a,i4,999(1x,g10.3))
3 format(999(1x,g10.3))
5 format(a,i6,99(1x,f9.4))
6 format(a,2i6,99(1x,g11.4))

 end subroutine extract_fourier
!===================================================================
! subroutine rotate(d2,s3,d2r)
! The relation is: D^al,be (tau,tau',S(k))=S^al,al' S^be,be' D^al',be' (iS(tau),iS(tau'),k)
! implicit none
! complex(r15), dimension(:,:), intent(in) :: d2
! real(r15), dimension(3,3), intent(in) :: s3
! complex(r15), dimension(:,:), intent(out) :: d2r
! integer ndim,i,j,al,be,tau,taup,als,bes,taus,taups

! ndim=size(d2(1,:))

! end subroutine rotate
!===================================================================
 subroutine get_relation_dynstar(nat,iop,s3,ndim,coefs)
! Relates the "rotated" dynamical matrices D(k*) to and D(k)
! if D is stored in 1D array (1:ndim) where ndim=3N(3N+1)/2, then D(k*) = coeff*D(k)
! The relation is: D^al,be (tau,tau',S(k))=S^al,al' S^be,be' D^al',be' (iS(tau),iS(tau'),k)
! where iS is the inverse or transpose of S
! the 1D compression is: l(=1:ndim) = l=(i-1)*3N+j, for i=1:3N ; j=i:3N
! reverse transformation for given l: j=mod(l,3N) and i=(l-j)/3N + 1 and d(i,j)uk(j)=-fk(i)
!! iop=symmetry operation number (1<iop<48)
!! symk = corresponding rotation matrix transformed to cartesian basis
!! ndim = size of the uppder diagonal dynamical matrix = 3*nat*(3*nat+1)/2
 use geometry
 use atoms_force_constants
 implicit none
 integer, intent(in) :: ndim,nat,iop
 real(r15), intent(in) :: s3(3,3)
 real(r15), intent(out):: coefs(ndim,ndim)
 integer l,c,tau,taup,i,j,ic,jc,al,be,als,bes,taus,taups

 if(ndim .ne. (3*nat*(3*nat+1))/2 ) then
     write(*,*)'ERROR in get_relation_dynstar: nat, ndim=',nat,ndim
     stop
 endif

 coefs=0

 do tau=1,nat
 do al=1,3
    i=3*(tau-1)+al

 do taup=1,nat
 do be=1,3
    j=3*(taup-1)+be   ! these 4 loops scan the line l

    l=(i-1)*3*nat+j  ! line index of D(k*)
!   rttp=v2a(atom0(tau)%equilibrium_pos-atom0(taup)%equilibrium_pos)

    do c=1,ndim      ! this scans the column c
! for this column index find corresponding alpha, tau etc...

       jc=mod(c,3*nat)
       ic=(c-jc)/(3*nat) + 1

       als  = mod(ic,3)
       taus = (ic-als)/3 + 1
       bes  = mod(jc,3)
       taups= (jc-bes)/3 + 1

!      rtstps=v2a(atom0(taus)%equilibrium_pos-atom0(taups)%equilibrium_pos)

! we must make sure S(tau)=taus and S(taup)=taups ...
       if ( iatomop(tau,taus) .eq. iop .and. iatomop(taup,taups) .eq. iop) then
          coefs(l,c)=s3(al,als)*s3(be,bes)
       endif
!     if (length(rtstps-matmul(s3,rttp)) .lt. 1d-4) then
!        write(*,6)'r_taus,taups - S(r_tautaup) =', rtstps-matmul(s3,rttp)
!        matched=.true.
!     elseif (length(matmul(s3,rtstps)-rttp)) .lt. 1d-4) then
!        write(*,6)'S(r_taus,taups) - r_tautaup =', matmul(s3,rtstps)-rttp

    enddo

 enddo
 enddo
 enddo
 enddo

6 format(a,9(1x,f10.4))
7 format(a,2(1x,i4),a,2(1x,i4))

 end  subroutine get_relation_dynstar
!===================================================================
 subroutine find_in_mesh(kstar,nkmesh,kmesh,nkstar)
 use geometry
 implicit none
 integer, intent(in) :: nkmesh
 integer, intent(out):: nkstar
 real(r15), intent(in) :: kstar(3),kmesh(3,nkmesh)
 integer i

 nkstar=0
 do i=1,nkmesh
    if(length(kstar-kmesh(:,i)) .lt. 1d-4) then
       nkstar=i
       exit
    endif
 enddo

 if (nkstar.eq.0) then
    write(*,'(a,3(1x,f9.3),a)')' kvector ',kstar,' not found among kmesh'
    stop
 endif

 end subroutine find_in_mesh

 end module fourier

