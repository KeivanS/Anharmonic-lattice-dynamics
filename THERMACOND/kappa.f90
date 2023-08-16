!=========================================================
 program kappa_sy
! we define the GF as the displ-displ correlation function (pure real):
! G0(k,wn,aa)=2w_ka/betahbar(w_ka^2+wn^2-2w_ka S(k,wn;aa)) with wn=2 pi n/beta hbar (finite T)
! the self-energy is then (off diagonal in the polarization modes a,b :
!S(k,wn;ab)=18beta^2 sum_{1,2} V(1,2,-ka)V(-1,-2,kb) [G01(w1)G02(w2)d(1+2-k)]
!iwn-> w+ie ; S=-D+iT then
! Dn(k;ab)=-hbar^2 beta/(16 N)*1/sqrt(wka*wkb) sum_{1,2} V(1,2,-ka) V(-1,-2,kb)*
! 1/(w1 w2) delta(k1+k2-k) *
![(n1+n2+1)/(w1+w2-w)+(n1+n2+1)/(w1+w2+w)+(n2-n1)/(w1-w2+w)+(n2-n1)/(w1-w2-w)]
! and
!Tn(k;ab)=-pi hbar/(16 N)*1/sqrt(wka*wkb) sum_{1,2} V(1,2,-ka) V(-1,-2,kb)*
! 1/(w1 w2) delta(k1+k2-k) *
![(n1+n2+1)d(w1+w2-w)-(n1+n2+1)/(w1+w2+w)+(n2-n1)d(w1-w2-w)-(n2-n1)/(w1-w2+w)]
! this was half the inverse lifetime (its diagonal element)
! units are eV/A^2,3 for FCs and will be kept throughout except for output files
! where they will be converted into cm^-1; same holds for FC3 (in eV/A^3)
! frequencies as well as perturbation V3, tempearature and self energies are in cm^-1
! g's and kpoints have a factor of 2pi included in them so that exp(ig.r)=1
! this version calculates kappa, dos and lifetimes within the irreducible FBZ
! it also includes the full (iterative) solution of Boltzmann equation
! by Keivan Esfarjani (Nov 11th 2010)
!
! TODO LIST:
! substitute om=w_q+Delta_q+I Gamma in the argument of the self-energy see how much
! the resutls change. i.e. do the self-consistent calculation, that would 
! also get rid of the ambiguity in the eta value.
! check Im(self) against the diagonal matrix elements of the collision matrix
! write the iteration subroutine (in effect inversion matrix) and solve the
! full BTE, to compare with RTA
! Profile the code to see if it can be optimized
! Add the contribution of 4th order term in the self energy (real part)->
! affects the cross section
! Add the option of adding an extra scattering (impurity of boundary) to the
! collision matrix and RTA to get the total phonon lifetimes and thermal
! conductivity
! Add the contribution of born charges to the cubic lifetimes (zero at zone
! center, and about 10% of total at zone boundary)
! review the QHA part... requires the FC calculation for different volumes but
! we can also use the cubic FCs to calculate FCs at these volumes -> Add this
! sort modes correctly according to their polarization not order ofeigenvalues.
!
! define an etamin=etamax/50 and etamax=wmax/min(nk1,nk2,nk3) and use the 
! self-consistent scheme but keeping always eta between etamax and etamin
!
! output the momentum and energy resolved jdos like the caltech paper
!
! output the memory usage
 use io2
 use kpoints
 use om_dos
 use force_constants_module
 use atoms_force_constants
 use phi3
 use phi3_sy
 use eigen
 use params
 use constants
 use born
! use exactBTE ! sy
 use exactBTE2
 use mod_ksubset2
 use tetrahedron


 implicit none
 integer i,j,mx,la,ncm,ierr,i2,j2,k2,inside,indexg,itmp,igam,ig
 real(8), allocatable :: dkc(:),omega(:),sigp(:)
 integer, allocatable :: iw(:)
 real(8) temp,tempk,dq,qx,qy,qz,kappa_klem,kappa_q(3,3),q0(3),q1(3),om0,om1,vel0
 real(8) ommax,mass_cell,vgmax,mfp,sigma,qq(3),cross_section
 complex(8) selfw
 real(8) , allocatable :: tau_klem_inv(:,:),kap(:,:,:),vgr(:,:)
 real(8) , allocatable :: c_matrix(:,:),inv_relax_t(:),evc(:),evc1(:), evc2(:), fw(:),taud(:),evlt(:)
 complex(8) , allocatable :: nself(:),uself(:),self(:,:,:,:),evct(:,:)
! real(8), allocatable :: eivr(:,:,:)
 real cputim
 character now*10,today*8,zone*5,ext*2,wxt*3

 integer ksubset(2), nv3_split, ksubset_bs(2)! sy 
 integer mergemode, num_files                   ! sy
 real(8), allocatable :: tau_inv_n(:,:), tau_inv_u(:,:)   ! sy
 integer convergence,iter_cont ! sy

 !mine
integer safo,nss,nsl,nxyz,nn, narms,kvecop(48)
real(8) nqsi(3), fsRTA(3),kvecstar(3,48),kvecstarsd(3,48)
 integer i1,j1,k1,vs,vsdf
 integer ss,si,sa,ms,sd,kvecopsd(48)
 real(8) fs(3),qs(3), qsi(3),ve2(3), we2(3)
 real(8) kvecstars(3,48),qss(3),difvs
 integer ds, kvecops(48)
 integer ni,i3,j3,k3,nz,ni0,sz,ssz,szz
integer rs,klrs,col1,col2
integer colkil,coll1,coll2
 !integer nkc, ndyn

! integer, allocatable ::  mapibzs(:),mapinvs(:)
  !!integer nibz
 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(cputim)

  open(ufc2,file='fc2.dat'    ,status='old')
  open(ufc3,file='fc3.dat'    ,status='unknown')
  open(ulog,file='phonlog.dat',status='unknown')
  open(utimes,file='times.dat' ,status='unknown')
  open(ueiv,file='eivals.dat',status='unknown')
! open(uklem,file='klem.dat',status='unknown')
! open(fbz,file='KPOINT.FBZ',status='unknown')
! open(ibz,file='KPOINT.IBZ',status='unknown')
  open(ibs,file='KPOINT.BS',status='unknown')
  open(ucors,file='eiv-coarse.dat',status='unknown')
! open(ufine,file='eiv-fine.dat',status='unknown')
  open(ugrun,file='mode_grun.dat',status='unknown')
  open(ugrun+1,file='all_grun.dat',status='unknown')
  open(ualph,file='alpha.dat',status='unknown')
! open(urate,file='rates.dat',status='unknown')
  open(debug,file='debug.dat',status='unknown')       ! sy. temp
  open(self_detail,file='selfenergy.dat',status='unknown')   ! sy.
  open(ksubset_bs_inp,file='ksubset_bs.inp',status='unknown')   ! sy
 
  write(self_detail,*) 'T(K) nk, la, v33_square_sum, delta1_sum, delta2_sum, delta3_sum, delta4_sum, delta_tot_sum'  ! sy


 write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' Program kap7 was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(ulog,'(a,f12.4)')' STARTING TIME OF THE PROGRAM IS ',cputim
 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' Program kap7 was launched at ',today(1:4)//'/' &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(utimes,'(a,f12.4)')' At the start of the Program,        TIME IS ',cputim

2 format(a,2x,f19.12,9(2x,g11.5))
3 format(i6,99(1x,g11.5))
4 format(99(2x,2(1x,f7.3)))
5 format(99(1x,f7.3))
6 format(9(2x,g12.6))
7 format(i6,9(2x,g11.5))
8 format(2i6,99(1x,g11.5))
9 format(i6,i6,3(1x,f7.3),1x,9(1x,g10.4))
!mine
20 format(2i6,2x,99(1x,f9.3))

! this is used to convert v33 to cm^_1
  const33 = 1d30*(hbar/uma)**1.5d0/(sqrt(ee*1d20/uma))**1.5d0*ee/h_plank/100/c_light


  call read_params
  print*,' file params.phon read'

  call read_input_fit
  print*,' file params.inp read'


  allocate(zeu(3,3,natoms0)) !JS
  call read_born
  print*,' file params.born read'


  call read_lattice
 

  nkc = nc1*nc2*nc3
  dq = (volume_g /nkc)**(1d0/3)
  write(ulog,*)' nkc, dq(ang^-1)=', nkc,dq

  call read_fc23

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After readings,                     TIME IS ',cputim

  ndyn = 3*natoms0


!---------------------------------------------------------------
! do a band structure calculation along the symmetry directions

  call make_kp_bs2
  write(ulog,*)' Kpoints for band structure generated from kpbs.in'
  call allocate_eig_bs(ndyn,nkp_bs,ndyn) !3rd index is for eivecs needed for gruneisen

  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,ndyn,eigenvec_bs,ueiv,veloc)

  call gruneisen(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,eigenvec_bs,ugrun,grun_bs)

  deallocate(veloc)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After band structure,               TIME IS ',cputim


!write(debug,*) 'eigen vectors...nkp, eigenvec_bs'  ! sy
!do i=1,nkp_bs   ! sy
!write(debug,119) 'i,kp,eigenvec_bs',i,kp_bs(:,i),real(eigenvec_bs(1,1,i)), real(eigenvec_bs(4,1,i)), real(eigenvec_bs(2,1,i)), real(eigenvec_bs(5,1,i)), real(eigenvec_bs(3,1,i)), real(eigenvec_bs(6,1,i))  ! sy
!enddo ! sy
!119 format(99(2x,g10.3)) ! sy


!---------------------------------------------------------------
! calculation of eigenvalues, Gruneisen in the FBZ, and DOS using the coarse mesh

! first generate the COARSE mesh, kpc and set the weights
  nkc = nc1*nc2*nc3
  allocate(kpc(3,nkc),mappos(nkc),wk(nkc))

! for even nci, they are shifted by -(g1+g2+g3)/2
  call make_kp_reg(nc1,nc2,nc3,shftx,shfty,shftz,kpc,wk)        ! original one.
!  call make_kp_reg_tet(nc1,nc2,nc3,shftx,shfty,shftz,kpc,wk,tet) ! tetrahedron method 

  call get_negatives(nkc,kpc,mappos,npos)
!write(*,*)'nibz' nibz
! ************** BEGIN: USE IBZ MESH ************************
!mapibzs=0
!mapinvs=0
  !!!!allocate(mapibz(nc1*nc2*nc3),mapinv(nc1*nc2*nc3))  for mods2-final.f90
!mapibzs=0
!mapinvs=0
      call get_weights(nc1*nc2*nc3,kpc)
!!!!call get_weights(nc1*nc2*nc3,kpc,nibz,narms,kvecstar,kvecop)   !mods2-final.f90
    !!call get_weights(nc1*nc2*nc3,kpc,nibz,mapibz,mapinv,sd,kvecstarsd,kvecopsd)

write(*,*)'safoura-code kappa7 can get in here'

! test of group velocities by HF and finite diff

! allocate(vgr(3,ndyn))
! q0 =kpc(:,2 )
! call finitedif_vel(q0,ndyn,vgr)
! deallocate(vgr) !,evlt,evct)

! output is nibz vectors kibz of weight wibz; mapibz is their mapping to IBZ
! nibz = nkc
! allocate(wibz(nibz),kibz(3,nibz),mapibz(nkc),mapinv(nkc))
! do i=1,nkc
!    kibz(:,i) = kpc(:,i)  !mappos(i))
!    mapibz(i)=i
!    mapinv(i)=i
!    wibz(i)=wk(i)
! enddo
! ************** END:  USE IBZ MESH ************************

!---------------------------------------------------------------
 
!- for merge v33 mode by sy
! call read_merge(mergemode, num_files)           ! read mergemode, number of files, file names
! if (mergemode .eq. 1)  then
!    nv3 = nkc*nibz*ndyn**3
!    write(ulog,*) 'ENTERING MERGE MODE....'
!    call merge_v3(num_files)          ! merge splitted v3 files and make v3 for all k points in IBZ
!    stop                   ! stop this program
! endif
! - end of merge v33 mode by sy


! allocates eigenval,eigenvec and grun of size (ndyn,nibz) and veloc(3,ndyn,nkc)
!    call allocate_eig(ndyn,nibz,nkc)   ! last one (nkc) is for veloc
     call allocate_eig(ndyn,nkc,nkc)   ! last one (nkc) is for veloc
     allocate(dkc(nkc),iw(nkc))
     dkc=0
     write(ulog,*)' eigenval allocated, entering the loop for DOS calculation'


     call get_frequencies(nkc,kpc,dkc,ndyn,eigenval,ndyn,eigenvec,ucors,veloc)
!    do i=1,nkc
!       write(ulog,3)i,kpc(:,i),veloc(:,1,i)
!    enddo

! gruneisen in the volume is needed for the thermal expansion coefficient
!      call gruneisen(nkc,kpc,dkc,ndyn,eigenval,eigenvec,ugrun+1,grun)
!!    call gruneisen(nibz,kibz,dkc,ndyn,eigenval,eigenvec,ugrun,grun)
  do i=1,nkc
     iw(i)=i
  enddo
!tst call get_group_velocity(nc1,nc2,nc3,kpc,ndyn,nkc,iw,eigenval,veloc)
     deallocate(dkc,iw)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After om(k) calculation within FBZ, TIME IS ',cputim

! check to see if D(-k)=conjg(D(k))
     call check_mdyn(ndyn,nkc,kpc,eigenvec,eigenval) ! would not work with nkc.ne.nibz

     call write_eivecs(ndyn,nkc,kpc,eigenvec)

     call subst_eivecs(ndyn,nkc,eigenvec,kpc,npos,mappos)

! calculate dos and joint-dos (added over all bands)
     call set_omdos(ndyn,wmesh)  ! allocates om and dos arrays of size mesh
     open(udos,file='dos.dat',status='unknown')
     allocate(evc(nibz))
     dos = 0
     do la=1,ndyn
        j=0
        do i=1,nibz
          j=j+1
          evc(j)=eigenval(la,mapinv(i))
        enddo
        call calculate_dos(nibz,evc,wibz,wmesh,om,dos(la,:))
!       write(wxt,'(i2.2)')la
        dos(ndyn+1,:) = dos(ndyn+1,:) + dos(la,:)
     enddo
     call write_dos
     close(udos) ; deallocate(evc)
     write(*,*) 'dos1 done'

! calculate dos using tetrahedron method
     call allocate_tetra(nc1,nc2,nc3,wmesh,lamax)
     call make_kp_reg_tet(nc1,nc2,nc3,shftx,shfty,shftz,kpc,wk)
!     call eigen_tet(nc1,nc2,nc3,lamax,nkc)
     open(udos,file='dos_tet.dat',status='unknown')
     open(udos+1,file='dos2_tet.dat',status='unknown')
     call calc_tet(nc1,nc2,nc3,wmesh,wmax,lamax,udos,udos+1,kpc)
!     call dos_tet(nc1,nc2,nc3,wmesh,wmax,lamax,udos,kpc)
     ! need to delcare calct, doswt, tet
     close(udos)!; close(udos+1)
     write(*,*) 'dos2 done'

!stop

! calculate the thermal expansion coefficient using phonons and gruneisen
!    call read_grun(nkc,ndyn,kpc,grun)

!     call gruneisen_fc

!     call calculate_thermal(nkc,wk,ndyn,eigenval,grun,tmin,tmax,ualph)

     close(ualph)
     deallocate(grun)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After calculate_dos, and thermal    TIME IS ',cputim

     mx = ndyn*nibz !nkc
     allocate(evc(mx),wkf(mx))
     j=0
     do la=1,ndyn
     do i=1,nibz !nkc
         j=j+1
         evc(j)=eigenval(la,mapinv(i))
         wkf(j)=wibz(i)  !1d0/mx
     enddo
     enddo

!    call phase_space_lifetime
    write(*,*) 'phase space lifetime done'
!    call phase_space_lifetime_bs    ! sy
!    write(*,*) 'start to calculate phase_space_lifetime_bs_occ...'  ! sy
    !call jdos_bs_sy   ! sy
!    write(*,*) 'phase space lifetime calculation done...'   ! sy
    !write(ulog,*) 'phase space lifetime done...'     ! sy
!   call jdos_sy   ! sy
!    write(ulog,*) 'JDOS done...'    ! sy
  !  call matrix_kpt_sy    ! sy
   ! write(ulog,*) 'matrix_kpt done...'    ! sy

     call deallocate_eig_bs

!    open(udos,file='jdos.dat',status='unknown')
!    call calculate_jdos(mx,evc,wkf,wmesh,om,udos)
!    close(udos)
     deallocate(evc,wkf,dos)
    write(ulog,*) 'deallocate done'

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After calculate_jdos,               TIME IS ',cputim

!     if (dos_bs_only ) stop             ! sy
!    write(ulog,*) 'after stop point'    ! sy

!----------------------------------------------------
! the anharmonic V3 matrix will be calculated on the coarse mesh
! call allocate_phi3(nkc,ndyn,wmesh,nibz)
! nv3 = (nkc+1)*nkc*ndyn**3/2
  nv3 = nkc*nibz*ndyn**3

! test to see whether the eivecs at opposite kpoints are now conjugate
  call write_eivecs(ndyn,nkc,kpc,eigenvec)
! call write_eivecs(ndyn,nibz,kibz,eigenvec)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' Start V33 read or calculation,               TIME IS ',cputim

write(*,*) 'iter,split,read,write,calk', iter, split, readv3, writev3, calk

do i=1,nkc
  
    call get_k_info(kpc(:,i) ,NC,j,i1,j1,k1,g1,g2,g3,inside)
    write(*,*)'i,j,inside=',i,j,inside

enddo

do i=1,nkc
  
    qq=kpc(:,i)+g1+g2+g3
    call get_k_info(qq ,NC,j,i1,j1,k1,g1,g2,g3,inside)
    write(*,*)'i,j,inside=',i,j,inside

enddo

!do i=1,nkc
!    qq=kpc(:,i)
!    qq(1)=qq(1)+0.1
!    call get_k_info(qq ,NC,j,i1,j1,k1,g1,g2,g3,inside)
!    write(*,*)'i,j,inside=',i,j,inside
!
!enddo
!stop

if ( iter .eq. 0)  then
    
    if ( split .eq. 0 ) then  ! no split. K1's original code
        
        ! Go from phi_3 (R,tau,alpha) to phi_3 (k,lambda)
        if ( readv3 .eq. 0) then      ! calculate V33 from scratch and store
             write(ulog,*) 'MAIN: calculating V3...'
             call allocate_v33(nv3)
!            call calculate_v3_new (ndyn,nkc,kpc,eigenval,eigenvec)
!            call calculate_v3     (ndyn,nkc,kpc,eigenval,eigenvec)
             call calculate_w3_ibz (ndyn,nibz,kibz,nkc,kpc,eigenval,eigenvec)

        elseif ( readv3 .eq. 1) then      ! read V33 from file and store
             write(ulog,*) 'MAIN: reading V3 from a file...'
!            call read_v3_all_formatted  (nkc,nkc,ndyn)
!            call read_v3_all_unformatted(nkc,nkc,ndyn)
             if ( allocated(v33))  call deallocate_v33
!            call read_v3_new_formatted
             call read_v3_new_unformatted
        endif
!       if (readv3 .ge. 2) then self-energy calculation proceeds on the fly


    elseif ( split .eq. 1 )  then  ! split IBZ. 
        
        !- for calculating v33 in split mode. sy
        if ( readv3 .eq. 0 ) then   ! calculate v33_split from scratch and store
             write(*,*) 'entering v33 calculation in split mode'   ! temp
             write(ulog,*) 'MAIN: calculating v3_split...'
             call read_ksubset(ksubset)
             nv3_split=nkc*(ksubset(2)-ksubset(1)+1)*ndyn**3
             write(ulog,*) 'nkc, ksubset, nv3_split',nkc,ksubset(1),ksubset(2),nv3_split
             call allocate_v33(nv3_split)
             call calculate_w3_ibz_split (ksubset,nv3_split,ndyn,nkc,kpc,nkc,kpc,eigenval,eigenvec)
        !- end of calculating v33 in split mode. sy

        
        !- start calculating v33 on the fly in split mode. sy
        elseif (readv3 .eq. 2) then    ! calculate v33 on the fly in split k space
             write(ulog,*) 'MAIN: will v3_split on the fly...'
             call read_ksubset(ksubset)
             nv3_split=nkc*(ksubset(2)-ksubset(1)+1)*ndyn**3
             write(ulog,*) 'nkc, ksubset, nv3_split',nkc,ksubset(1),ksubset(2),nv3_split
        !- end calculating v33 on the fly in split mode. sy

        endif
    endif


elseif ( iter .eq. 1 ) then   ! for exact solution of BTE
    
    if ( readv3 .eq. -1 )  then
        call merge_Pmatrix    ! merge collision matrix. obsolete.

    elseif ( readv3 .eq. 0 ) then
        
        write(*,*) 'entering cal v33 in FBZ split...'
        call read_ksubset(ksubset)  ! read index of initial and final kpoints set in the shellscript
        nv3_split=nkc*(ksubset(2)-ksubset(1)+1)*ndyn**3
!       call allocate_v33_sq(ksubset(2)-ksubset(1)+1,nkc,ndyn,ndyn,ndyn)        
  ! k1    allocate(v33(ksubset(2)-ksubset(1)+1,nkc,ndyn,ndyn,ndyn))        
        call allocate_v33(nv3_split)
        call calculate_w3_fbz_split_sq (ksubset,nv3_split,ndyn,nkc,kpc,eigenval,eigenvec)
 !       call deallocate_V33_sq
         call cpu_time(cputim)  !date_and_time(time=tim)
         write(utimes,'(a,f12.4)')' After v33 calculation,                TIME IS ',cputim

    elseif ( readv3 .eq. 2) then     ! read v33sq.xxx.dat from the main directory and calculate v33sq*delta and write them into files.
        
        write(ulog,*) 'entering cal v33*delta...'
        call read_ksubset(ksubset) ! read ksubset
        write(ulog,*) 'read_ksubset done...'
        call calculate_v3sq_delta (ksubset,ndyn,nkc,kpc)  ! calculate v3sq*delta

        call cpu_time(cputim)  !date_and_time(time=tim)
        write(utimes,'(a,f12.4)')' After v33_delta calculation,            TIME IS ',cputim

    elseif  ( readv3 .eq. 1 .and. calk .eq. 0)  then     ! obsolete

        write(*,*) 'entering cal Pmatrix...'      

        call read_v3_new_unformatted_sq(nkc,ndyn)   ! this will allocate v33. also read nv3 (=nv3_split)
!        call allocate_P (nv3)
        write(*,*) 'read_v3 done...'

        call read_ksubset(ksubset)
        nv3_split=nkc*(ksubset(2)-ksubset(1)+1)*ndyn**3

        temp=tmin*k_b/(100*h_plank*c_light)   ! Kelvin to cm^-1. tmin in params.phon will be used for iterative solution
        call Pmatrix(ksubset,nkc,kpc,ndyn,eigenval,temp)           ! calculate collision matrix, save it to file
!        call Pmatrix2(ksubset,nkc,kpc,ndyn,eigenval,temp)           ! read v33, calculate P matrix and save it to file
!        call allocate_iter1(nkc*ndyn)
!        call allocate_iter1(nkc,ndyn)

!        call selfenergy2(ksubset,nkc,kpc,ndyn,eigenval,temp)
        
    endif

endif


call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After calculation/reading of V3,    TIME IS ',cputim

!    call check_sym_of_v3(ndyn)

!----------------------------------------------------
! do statistics on the distribution of v33 matrix elements
!   if ( readv3 .lt. 2) then
!      mx=1000
!      allocate(omega(nv3))
!      omega = abs(v33)
!      call histogram(nv3,omega,mx)
!      deallocate(omega)
!   endif 

!call cpu_time(cputim)  !date_and_time(time=tim)
!write(utimes,'(a,f12.4)')' After calculation of V33-histogram, TIME IS ',cputim

!----------------------------------------------------
! calculate the cubic matrix elements^2 at each frequency
   if (threemtrx .eq. 1 .and. readv3 .lt. 2) then

     open(udos,file='3mtrx.dat',status='unknown')
     do i=1,wmesh/2
        ommax = om(2*i)
        call three_ph_matrix(ommax,etaz,vgmax)
        write(udos,6)ommax,vgmax
     enddo
     close(udos)

   endif

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After calculation of 3-ph-matrix    TIME IS ',cputim


if (threemtrx .eq. 1 .and. readv3 .lt. 2) then
call matrix_kpt_sy  ! sy
endif

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After calculation of matrix_kpt_sy    TIME IS ',cputim

!----------------------------------------------------
! calculate the self-energy of the LA mode at (K=0) versus (omega and) T

   if( allocated(omega) ) deallocate(omega)
   if( allocated(uself) ) deallocate(uself)
   if( allocated(nself) ) deallocate(nself)
   allocate(omega(ndyn),uself(ndyn),nself(ndyn),kap(ndyn,3,3))

!----------------------------------------------------
! calculate the cross section at (K_red=qcros) versus omega and T
 if (cal_cross.eq.1) then

   open(ucross,file='cross-section.dat')
   allocate(sigp(ndyn))

   qq = qcros(1)*g1+qcros(2)*g2+qcros(3)*g3
   crloop: do itmp=1,ntemp
     tempk=tmin+(tmax-tmin)*(itmp-1)**3/(ntemp-1.0000001d0)**3   ! this temp is in Kelvin
     temp = tempk*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1
     write(ucross,*)"# la,tempk,omega,nself(la),uself(la),sigp(la)"
     do la=lamin,lamax
     do j=1,wmesh/2
        ommax=om(2*j)*1.2  ! take fewer om points to save time and go beyond band-range
!       call get_kindex(qq,nkc,kpc,igam)
!       omega(la)=sqrt(eigenval(la,igam)) * cnst
        call self_energy(qq ,la,ommax,temp,nself(la),uself(la),readv3)
        selfw = nself(la)+uself(la)
        sigp(la) = cross_section(qq,la,ommax,eigenval(la,igam),selfw)
        write(ucross,7)la,tempk,ommax,nself(la),uself(la),sigp(la)
     enddo
        write(ucross,*)" "
     enddo
   enddo crloop

   close(ucross)
   deallocate(sigp)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After calculation of cross(Gama,T)  TIME IS ',cputim

 endif


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!------------- Keivan's original code (serial calculation with RTA)
!not tested.

if  (iter .eq. 0  .and.  split .eq. 0 .and. calk .eq. 1) then  ! Keivan's original code

!----------------------------------------------------
! temperature loop for calculation of relaxation times and kappa
  open(ukap,file='kappa.dat',status='unknown')
  write(ukap,2)'# i, temperature(K) , kappa(la) (W/mK) ...    , kappa_total (W/mK) '
  mass_cell = sum( atom0(:)%mass )

! do la=1,ndyn
!    write(wxt,'(i2.2)')la
!    open(uslf+100+la,file='rself-'//wxt//'.temp')
!    open(uslf+200+la,file='iself-'//wxt//'.temp')
!    open(ukap+la,file='mode_kap_'//wxt//'.temp',status='unknown')
!    write(ukap+la,2)'# i, temperature(K) , kappa(la,i,j) (W/mK)'
! enddo
  open(uslf+1  ,file='rself.temp')
  open(uslf+500,file='iself.temp')
  open(umod,file='mode_kap.temp',status='unknown')
  write(umod,2)'# la, i, temperature(K) , kappa(la,i,j) (W/mK)'

  allocate(tau_klem_inv(ndyn,nkc))
  temploop: do itmp=1,ntemp
    allocate(evc(nibz))
    tempk=tmin+(tmax-tmin)*(itmp-1)**2/(ntemp-1.0000001d0)**2   ! this temp is in Kelvin
    temp = tempk*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1
    tau_klem_inv = 0d0
    modeloop: do la=1,ndyn
!     write(uslf+100+la,2)'#nk,tempk,omega(k,la),real(self(la)),normal,umklapp '
!     write(uslf+200+la,2)'#nk,tempk,omega(k,la),imag(self(la)),normal,umklapp,tau(ps),MFP(nm) '
      write(uslf+1  ,2)'#la,nk,tempk,omega(k,la),real(self(la)),normal,umklapp '
      write(uslf+500,2)'#la,nk,tempk,omega(k,la),imag(self(la)),normal,umklapp,tau(ps),MFP(nm) '
! SELF-ENEGY calculation on the coarse mesh
      do j=1,nibz
         call get_kindex(kibz(:,j),nkc,kpc,igam)
         if(igam.ne.mapinv(j)) then
            write(*,*)'MAIN: temp,mode,kibz=',tempk,la,j
            write(*,*)'get_kindex found igam-',igam,mapinv(j)
            stop
         endif

         omega(la) = sqrt(eigenval(la,igam)) * cnst
         call self_energy(kibz(:,j),la,omega(la),temp,nself(la),uself(la),readv3)
         evc(j)=2*aimag(nself(la)+uself(la))   ! used for inverse relaxation time
! find relaxation times of all bands over the full zone (stored in tau_klem_in)
         do k2=1,nkc
            if(mapibz(k2).eq.j) then
               if (tau_klem_inv(la,k2) .ne.0) then
                  write(ulog,*)'tau_inv already assigned!!',la,k2,tau_klem_inv(la,k2)
                  stop
               else
               tau_klem_inv(la,k2)=evc(j)
               endif
            endif
         enddo
         mfp=length(veloc(:,la,igam)) / evc(j) * 1d7   ! in nanometers
!        write(uslf+200+la,7)j,tempk,omega(la),aimag(nself(la)+uself(la)),aimag(nself(la)),aimag(uself(la)),1d10/c_light/evc(j) ,mfp
!        write(uslf+100+la,7)j,tempk,omega(la),real (nself(la)+uself(la)),real (nself(la)),real (uself(la))
         write(uslf+500,8)la,j,tempk,omega(la),aimag(nself(la)+uself(la)),aimag(nself(la)),aimag(uself(la)),1d10/c_light/evc(j) ,mfp
         write(uslf+1  ,8)la,j,tempk,omega(la),real (nself(la)+uself(la)),real (nself(la)),real (uself(la))
      enddo
      write(uslf+1  ,*)' '
      write(uslf+500,*)' '

! Thermal conductivity from the RTA to the solution of Boltzmann equation
!     call mode_thermal_conductivity(nkc,wk,ndyn,la,evc,veloc,eigenval,temp,kappa_q)

! this is from k's in the full FBZ
      call mode_thermal_conductivity(nkc,wk,tau_klem_inv(la,:),veloc(:,la,:),eigenval(la,:),temp,kappa_q)

      kap(la,:,:) = kappa_q(:,:)    ! here used as a dummy array
      omega(la)=(kap(la,1,1)+kap(la,2,2)+kap(la,3,3))/3
      write(umod,8)la,itmp,tempk,(kap(la,i,i),i=1,3),((kap(la,i,j),j=i+1,3),i=1,2)

    enddo modeloop

    write(ukap,3)itmp,tempk,(omega(la),la=1,ndyn),sum(omega)

    deallocate(evc)
  enddo temploop
  deallocate(tau_klem_inv)
  close(uslf+1  )
  close(uslf+500)
  close(umod)
  close(ukap)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After calculation of linewidth,     TIME IS ',cputim

!================================================================================
  close(ukap)

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(cputim)
write(ulog,'(a,f12.4)')' ENDING TIME OF THE PROGRAM IS ',cputim
write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' This program was ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
write(utimes,'(a,f12.4)')' At the end of the Program,          TIME IS ',cputim
write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' This program was ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

!----------------------------------------------------

 deallocate(kpc,wk)
call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After the kappa loop,               TIME IS ',cputim

 close(ulog)
 close(ueiv)
 close(ufc2)
 close(ufc3)
 close(utimes)
 close(ibz)
 close(bzf)
 close(fbz)
 close(ibs)
 close(ucors)
 close(ualph)
! close(urate)





!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!------------- Keivan's original code + split (split calculation with RTA) 
!not tested
elseif (iter .eq. 0 .and. split .eq. 1 .and. calk .eq. 1) then  ! Keivan's original code + distribute k-poitns to many cores
  
  open(ukap,file='kappa.dat',status='unknown')
  write(ukap,2)'# i, temperature(K) , kappa(la) (W/mK) ...    , kappa_total (W/mK) '
  open(k_kappa,file='k_kappa.dat',status='unknown')    ! sy
  write(k_kappa,*)'# T(K), ki, kpc(3), la, freq(cm-1), tau_inv_tot, tau_inv_normal, tau_inv_umklapp, group_vel(3), mfp(nm), cv, nBE, kxx,kyy,kzz,kxy,kxz,kyz, k_tot(W/mK)'  ! sy
  mass_cell = sum( atom0(:)%mass )

! do la=1,ndyn
!    write(wxt,'(i2.2)')la
!    open(uslf+100+la,file='rself-'//wxt//'.temp')
!    open(uslf+200+la,file='iself-'//wxt//'.temp')
!    open(ukap+la,file='mode_kap_'//wxt//'.temp',status='unknown')
!    write(ukap+la,2)'# i, temperature(K) , kappa(la,i,j) (W/mK)'
! enddo
  open(uslf+1  ,file='rself.temp')
  open(uslf+500,file='iself.temp')
  open(umod,file='mode_kap.temp',status='unknown')
  write(umod,2)'# la, i, temperature(K) , kappa(la,i,j) (W/mK)'

 ! allocate(tau_klem_inv(ndyn,nkc))
  temploop2: do itmp=1,ntemp

    call cpu_time(cputim)  !date_and_time(time=tim)      ! sy
    write(utimes,'(a,i1,a,i1,a,f15.8)')' Start Kappa calculation at',itmp,'/',ntemp,' TIME IS ',cputim   ! sy

    allocate(evc(nibz))
    allocate(evc1(nibz))     ! sy
    allocate(evc2(nibz))     ! sy
    allocate(tau_inv_n(ndyn,nkc))    ! sy
    allocate(tau_inv_u(ndyn,nkc))    ! sy
    allocate(tau_klem_inv(ndyn,nkc)) ! sy
    tau_inv_n=0; tau_inv_u=0    ! sy

    tempk=tmin+(tmax-tmin)*(itmp-1)**2/(ntemp-1.0000001d0)**2   ! this temp is in Kelvin
    temp = tempk*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1
    tau_klem_inv = 0d0
    modeloop2: do la=1,ndyn
!     write(uslf+100+la,2)'#nk,tempk,omega(k,la),real(self(la)),normal,umklapp '
!     write(uslf+200+la,2)'#nk,tempk,omega(k,la),imag(self(la)),normal,umklapp,tau(ps),MFP(nm) '
      write(uslf+1  ,2)'#la,nk,tempk,omega(k,la),real(self(la)),normal,umklapp '
      write(uslf+500,2)'#la,nk,tempk,omega(k,la),imag(self(la)),normal,umklapp,tau(ps),MFP(nm) '

! SELF-ENEGY calculation on the coarse mesh
!     do j=1,nibz
      do j=ksubset(1),ksubset(2)  ! modified by sy. for split mode
         call get_kindex(kibz(:,j),nkc,kpc,igam)
         if(igam.ne.mapinv(j)) then
            write(*,*)'MAIN: temp,mode,kibz=',tempk,la,j
            write(*,*)'get_kindex found igam-',igam,mapinv(j)
            stop
         endif

         omega(la) = sqrt(eigenval(la,igam)) * cnst
!! k1
  !      call self_energy2(j,kibz(:,j),la,omega(la),temp,nself(la),uself(la),readv3)
    call function_self_w2_sy(j,kibz(:,j),la,omega(la),temp,nself(la),uself(la))   ! sy
!! k1
         evc(j)=2*aimag(nself(la)+uself(la))   ! used for inverse relaxation time
         evc1(j)=2*aimag(nself(la))    ! sy
         evc2(j)=2*aimag(uself(la))    ! sy
         write(ulog,*) 'self energy calculation done....la, ibz',la,j    ! sy

! find relaxation times of all bands over the full zone (stored in tau_klem_in)
        do k2=1,nkc
            if(mapibz(k2).eq.j) then
               if (tau_klem_inv(la,k2) .ne.0) then
                  write(ulog,*)'tau_inv already assigned!!',la,k2,tau_klem_inv(la,k2)
                  stop
               else
               tau_klem_inv(la,k2)=evc(j)
               tau_inv_n(la,k2)=evc1(j)     ! sy
               tau_inv_u(la,k2)=evc2(j)     ! sy
               endif
            endif
         enddo

         mfp=length(veloc(:,la,igam)) / evc(j) * 1d7   ! in nanometers
!        write(uslf+200+la,7)j,tempk,omega(la),aimag(nself(la)+uself(la)),aimag(nself(la)),aimag(uself(la)),1d10/c_light/evc(j) ,mfp
!        write(uslf+100+la,7)j,tempk,omega(la),real (nself(la)+uself(la)),real (nself(la)),real (uself(la))
         write(uslf+500,8)la,j,tempk,omega(la),aimag(nself(la)+uself(la)),aimag(nself(la)),aimag(uself(la)),1d10/c_light/evc(j) ,mfp
         write(uslf+1  ,8)la,j,tempk,omega(la),real (nself(la)+uself(la)),real (nself(la)),real (uself(la))
      enddo
      write(uslf+1  ,*)' '
      write(uslf+500,*)' '

! Thermal conductivity from the RTA to the solution of Boltzmann equation
!     call mode_thermal_conductivity(nkc,wk,ndyn,la,evc,veloc,eigenval,temp,kappa_q)

! this is from k's in the full FBZ
      write(ulog,*) 'entering kappa calculation...'    ! sy
      call mode_thermal_conductivity_sy(tempk,la,nkc,kpc,wk,tau_klem_inv(la,:),tau_inv_n(la,:),tau_inv_u(la,:),veloc(:,la,:),eigenval(la,:),temp,kappa_q)
      kap(la,:,:) = kappa_q(:,:)    ! here used as a dummy array
      omega(la)=(kap(la,1,1)+kap(la,2,2)+kap(la,3,3))/3
      write(umod,8)la,itmp,tempk,(kap(la,i,i),i=1,3),((kap(la,i,j),j=i+1,3),i=1,2)

    enddo modeloop2

    write(ukap,3)itmp,tempk,(omega(la),la=1,ndyn),sum(omega)

    deallocate(evc)
    deallocate(evc1)    ! sy
    deallocate(evc2)    ! sy
    deallocate(tau_inv_n)    ! sy
    deallocate(tau_inv_u)    ! sy
    deallocate(tau_klem_inv)
 enddo temploop2
  close(uslf+1  )
  close(uslf+500)
  close(umod)
  close(ukap)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After calculation of linewidth,     TIME IS ',cputim


!================================================================================

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(cputim)
write(ulog,'(a,f12.4)')' ENDING TIME OF THE PROGRAM IS ',cputim
write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' This program was ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
write(utimes,'(a,f12.4)')' At the end of the Program,          TIME IS ',cputim
write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' This program was ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

!----------------------------------------------------

 deallocate(kpc,wk)
call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After the kappa loop,               TIME IS ',cputim

 close(ulog)
 close(ueiv)
 close(ufc2)
 close(ufc3)
 close(utimes)
 close(ibz)
 close(bzf)
 close(fbz)
 close(ibs)
 close(ucors)
 close(ualph)
! close(urate)






!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!------------- Exact solution of BTE by sy
!elseif (iter .ne.0  .and. readv3 .eq. 1 .and. calk .eq. 1) then  ! Exact solution of BTE sy
elseif (iter .eq.1  .and. readv3 .eq. 1 .and. calk .eq. 1) then  ! Exact solution of BTE sy

write(ulog,*) 'entering solving BTE iteratively'
write(ulog,*) 'convergence criteria',conv_diff_kap,conv_max_diff_kap
write(ulog,*) 'continuous iterations',conv_iter
write(ulog,*) 'update ratio',update

call allocate_iter1 (nkc,ndyn)    ! allocate memory for iterative solution
call allocate_iter2 (nkc,nkc,ndyn,ndyn,ndyn)

write(ulog,*) 'allocate_iter done'

!call read_col_matrix(nkc,ndyn)    ! read col_matrix.dat
!write(*,*) 'read_col_matrix done'
!write(*,*) 'P1(1)=', P1(1,1,1,1,1)

tempk = tmin     ! use tmin in params.phon as temperature. tmax is not being used.

call distribution(nkc,ndyn,tempk)  ! calculate BE distribution function for each state
write(ulog,*) 'distribution done'


!call mode_thermal_RTA(tempk,nkc,ndyn,kpc)  ! calculate thermal conductivity based on relaxation time from direct_RTA
!write(ulog,*) 'mode_thermal_RTA done'
!write(*,*) 'mode_thermal_RTA done'


!call cal_Q(nkc,ndyn)              ! calculate Q value
!call cal_Qm2(nkc,ndyn,kpc,tempk)        ! read v33 files and calculate P matrix on the fly
!call allocate_iter2(nkc,nkc,ndyn,ndyn,ndyn)   ! allocate collision matrix P1, P2
!P1=0.0; P2=0.0
!call cal_Qm3(nkc,ndyn,kpc,tempk)              ! calculat P1, P2 first, then calculate sum of diagonal terms in collision matrix
!call cal_Qm3_tet(nkc,ndyn,kpc,tempk)           ! same as cal_Qm3, but tetrahedron method is used.
call cal_Qm3_tet2(nkc,ndyn,kpc,tempk)           ! read v33sq_delta.xxx.dat, calculate P1, P2, Qvalue, and save them into memory


write(ulog,*) 'cal_Q done'


call RTA(nkc,kpc,ndyn,tempk,veloc,eigenval)               ! calculate RTA solution

write(ulog,*) 'RTA done'











!if(iter .eq. 2) then  ! use collision matrix 


!******* for collision(FBZ*FBZ) and only x-axis for velocity and F *******!

!allocate(col_matrix(nkc*ndyn,nkc*ndyn),RHS2(nkc*ndyn),xx2(nkc*ndyn),sig(nkc*ndyn) )
!call cal_RHS_x(nkc,kpc,ndyn,tempk,veloc,eigenval,RHS2)
 
!xopen(451,file='RHStest.dat')
!xklrs=nkc*ndyn
!xdo rs=1,klrs
!x write(401,*)rs, RHS2(rs)
!xenddo
!xclose(451) 

!call cal_collisionM_FBZ_x(nkc,kpc,ndyn,col_matrix)

!xopen(453,file='collisiontest.dat')
!xdo col1=1,klrs
!x   do col2=1,klrs
!x write(403,*)col1,col2,col_matrix(col1,col2)
!xenddo
!xenddo
!xclose(453)

!svdcut = 1d-9
!call svd_set(nkc*ndyn,nkc*ndyn,col_matrix,RHS2,xx2,sig,svdcut,errorr,ermax,sigma,'svd_output_FBZx.dat')

!*******  end (for FBZ and only x-axis)  *******!



!****** for collision(FBZ,FBZ) with all axes for velocity and F *******!

!allocate(col_matrix(nkc*ndyn*3,nkc*ndyn*3),RHS2(nkc*ndyn*3),xx2(nkc*ndyn*3),sig(nkc*ndyn*3))
!call cal_RHS(nkc,kpc,ndyn,tempk,veloc,eigenval,RHS2)
!call cal_collisionM2(nkc,kpc,ndyn,col_matrix)
!svdcut = 1d-9
!call svd_set(nkc*ndyn*3,nkc*ndyn*3,col_matrix,RHS2,xx2,sig,svdcut,errorr,ermax,sigma,'svd_output_FBZ.dat')

!******* end (for FBZ and all axes) *******!



!******* for collision matrix(IBZ*FBZ)*Symmetry matrix(FBZ*IBZ) & RHS in IBZ *******!

!*allocate(col_matrix(nibz*ndyn*3,nkc*ndyn*3),RHS2(nibz*ndyn*3),sy_matrix(nkc*ndyn*3,nibz*ndyn*3))
!call sy_matrix_IBZ_FBZ2(nibz,kibz,nkc,kpc,ndyn,sy_matrix)
!*call sy_matrix_IBZ_FBZ2(nibz,kibz,nkc,kpc,ndyn,sy_matrix)
!*call cal_RHS_IBZ(nibz,kibz,ndyn,tempk,veloc,eigenval,RHS2)
!*call cal_collisionM2_IBZ2(nkc,kpc,nibz,kibz,ndyn,col_matrix)
!*allocate(coll(nibz*ndyn*3,nibz*ndyn*3))
!*coll=MATMUL(col_matrix,sy_matrix)
!*deallocate(col_matrix,sy_matrix)

!*open(607,file='colsy.dat')

!*colkil=nibz*ndyn*3
!*do coll1=1,colkil
!*   do coll2=1,colkil
!*      write(607,*) coll1,coll2,coll(coll1,coll2)
!*   enddo
!*enddo
!*close(607)


!*allocate(xx2(nibz*ndyn*3),sig(nibz*ndyn*3))
!*svdcut = 1d-9
!*call svd_set(nibz*ndyn*3,nibz*ndyn*3,coll,RHS2,xx2,sig,svdcut,errorr,ermax,sigma,'svd_output_IBZ.dat')
!allocate(F_MRHS(nkc*ndyn*3))
!F_MRHS = MATMUL(sy_matrix,xx2)
!call mode_thermal_noniter2(nkc,ndyn,tempk,veloc,eigenval,F_MRHS) 
!*call mode_thermal_noniter22(nibz,nkc,kibz,kpc,ndyn,tempk,veloc,eigenval,xx2)
!*call write_thermal_noniter2(nkc,ndyn,tempk,veloc,eigenval,kpc)

!****** end ( C-matrix* S-matrix) *******!


!******* collision matrix and symmetry matrix in one matrix (IBZ*IBZ) & RHS in IBZ *******!

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' Before svd approach,                TIME IS ',cputim

allocate(xx2(nibz*ndyn*3),sig(nibz*ndyn*3))
allocate(col_sy(nibz*ndyn*3,nibz*ndyn*3),RHS2(nibz*ndyn*3))
call sycollisionM(nibz,kibz,ndyn,nkc,kpc,col_sy)
call cal_RHS_IBZ(nibz,kibz,ndyn,tempk,veloc,eigenval,RHS2)
!allocate(xx2(nibz*ndyn*3),sig(nibz*ndyn*3))
svdcut = 1d-9
call svd_set(nibz*ndyn*3,nibz*ndyn*3,col_sy,RHS2,xx2,sig,svdcut,errorr,ermax,sigma,'svd_output_IBZ2.dat')
call mode_thermal_noniter23(nibz,nkc,kibz,kpc,ndyn,tempk,veloc,eigenval,xx2)
call write_thermal_noniter2(nkc,ndyn,tempk,veloc,eigenval,kpc)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' after svd approach,                 TIME IS ',cputim



!allocate( col_matrix(nkc*ndyn*3,nkc*ndyn*3), RHS2(nkc*ndyn*3),xx2(nkc*ndyn*3),sig(nkc*ndyn*3) )


!call cal_RHS(nkc,kpc,ndyn,tempk,veloc,eigenval,RHS2)           
!open(401,file='RHStest.dat')
!klrs=nkc*ndyn*3
!do rs=1,klrs
!   write(401,*)rs, RHS2(rs)
!enddo
!close(401)



!call cal_collisionM(nkc,kpc,ndyn,col_matrix)

!open(403,file='collisiontest.dat')

!do col1=1,klrs
!   do col2=1,klrs
!   write(403,*)col1,col2,col_matrix(col1,col2)

!enddo
!enddo
!close(403)


!svdcut = 1d-9
!call svd_set(nkc*ndyn*3,nkc*ndyn*3,col_matrix,RHS2,xx2,sig,svdcut,errorr,ermax,sigma,'svd_output.dat')

write(ulog,*) 'SVD done'

!******* end (collision matrix and symmetry matrix in one matrix (IBZ*IBZ))*******!










!write the distribution functioni xx  here

!elseif (iter.eq.1) then  ! iterative solution here


!endif 
F1(:,:,:)=1.0*F_RTA(:,:,:) ! for iteration (initial guess)



convergence=0              ! if converged, convergence = 1
iter_cont=0                ! number of iterations while satisfying convergence criteria


!open(170,file='fRTA.dat')

  ! allocate(dkc(nibz))
  ! dkc=0

do i=1, max_iter     ! max_iter : maximum iteration allowed by params.phon

    write(ulog,*) '================= iteration',i

!   call cal_F2(nkc,kpc,ndyn)
  !!nkc=nc1*nc2*nc3 
   !call get_weights(nkc,kpc)
   !call cal_F2m2(nkc,kpc,nibz,kibz,ndyn,tempk)              ! calculate F value
   call cal_F2m23(nkc,kpc,nibz,kibz,ndyn,tempk)
   

    write(ulog,*) 'cal_F2 done'
   F1_old(:,:,:) = F1(:,:,:)         ! F1_old for convergence check.
   F1(:,:,:) = F1(:,:,:) + update*(F2(:,:,:)-F1(:,:,:))  ! update F. update is read from params.phon (0.5 is recommended for Bi)

   !allocate(dkc(nibz))
   !dkc=0
   !call get_frequencies(nibz,kibz,dkc,ndyn,eivalibz,ndyn,eivecibz,velocibz)

 !@@@  call mode_thermal_iter(nibz,nkc,ndyn,tempk,veloc,eigenval)  ! thermal conductivity calculation based on F values from above.
   call mode_thermal_iter2(nibz,nkc,kibz,kpc,ndyn,tempk,veloc,eigenval)
   write(ulog,*) 'mode_thermal_iter done...'

!   call check_conv(i,nkc,kpc,ndyn,tempk,convergence,iter_cont)
!@@@   call check_conv2(i,nkc,kpc,ndyn,tempk,convergence,iter_cont)   ! check convergence.
   call check_conv23(i,nkc,kpc,nibz,kibz,ndyn,tempk,convergence,iter_cont)

   write(ulog,*) 'check_conv done...'

   !write(ulog,*) 'i,convergence,iter_cont',i,convergence,iter_cont
   !write(ulog,*) 'conv_diff_kap,conv_max_diff_kap,conv_iter',conv_diff_kap,conv_max_diff_kap,conv_iter

   call write_thermal(i,nkc,ndyn,tempk,veloc,eigenval,kpc)         ! write output file for thermal conductivity

   if (convergence .ne. 0 .and. iter_cont .gt. conv_iter) then      ! If converged
        exit
   endif
 
enddo



!close(170)
  !deallocate(dkc)
call deallocate_iter2

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After calculation of linewidth,     TIME IS ',cputim


!================================================================================

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(cputim)
write(ulog,'(a,f12.4)')' ENDING TIME OF THE PROGRAM IS ',cputim
write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' This program was ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
write(utimes,'(a,f12.4)')' At the end of the Program,          TIME IS ',cputim
write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' This program was ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

!----------------------------------------------------

 deallocate(kpc,wk)
call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After the kappa loop,               TIME IS ',cputim

 close(ulog)
 close(ueiv)
 close(ufc2)
 close(ufc3)
 close(utimes)
 close(ibz)
 close(bzf)
 close(fbz)
 close(ibs)
 close(ucors)
 close(ualph)
! close(urate)



endif

 end program kappa_sy



!===============================================================
  subroutine get_kindex(q,nkc,kpc,ik)
  use  geometry
  implicit none
  integer ik,nkc,i,j
  real(8) kpc(3,nkc),q(3)

  ik=0
  mainlp: do i=1,nkc
     do j=1,3
        if (.not. (abs(q(j)-kpc(j,i)) .myeq. 0d0) ) exit
        if (j.eq.3) then
           ik=i
           exit mainlp
        endif
     enddo
  enddo mainlp
  if (ik.eq.0) then
     write(*,*)'GET_INDEX: could not find the index of',q
     stop
  endif

  end subroutine get_kindex
!===============================================================
  subroutine mean_sd(x,mean,sd)
  implicit none
  integer n,i
  real(8) x(:),mean,sd,sd2

  n=size(x)
  mean = sum(x)/n
  sd2  = sum(x*x)/n
  sd   = sqrt(sd2 - mean*mean)

  end subroutine mean_sd
!=====================================================
  subroutine histogram(m,x,mesh)
! calculates the histogram of the data set x (needs not be sorted)
! and writes the distribution function in the file 'histo.dat'
  implicit none
  integer, intent(in):: m,mesh
  real(8), intent(in):: x(m)
  integer i,j,cnt,unit
  real(8) xmin,xmax,dx,sume,cume,e(mesh)

  unit=123
  open(unit,file='histo.dat')
  write(unit,*)'# j,  x(j) , e(j) , e(j)/sume/dx,accumulated e '

  cnt = size(x)
  xmax= maxval(x)
  xmin= minval(x)
  dx  = (xmax-xmin)/mesh
  write(*,5)'HISTO: size,xmin,xmax,dx=',cnt,xmin,xmax,dx

  e=0
  do i=1,cnt
! j=position of the bin to which x(i) belongs
     j= int((x(i)-xmin)/dx)+1
     if(j.eq.mesh+1) j=mesh
     if (j.lt.1 .or. j.gt.mesh) then
        write(*,*)'j is out of range ',j
        stop
     endif
     write(unit,4)j,xmin+(j-0.5)*dx,e(j),e(j)/sume/dx,cume
  enddo
 4 format(i8,9(2x,g12.6))
 5 format(a,i8,9(2x,g12.6))

  close(unit)

  end subroutine histogram
!=====================================================
!! k1
! subroutine self_energy(q_ibz,k,la,omega,temp,nself,uself,iread)
  subroutine self_energy(k,la,omega,temp,nself,uself,iread)
!! k1
  implicit none
  integer iread,la,q_ibz
  complex(8) nself,uself
  real(8) temp,omega,k(3)

  if(iread .lt. 2) then
!    write(*,*)' entering self_new '
    call function_self_new(k,la,omega,temp,nself,uself)
!    call function_self_sc(k,la,omega,temp,nself,uself)
  elseif(iread.eq.2)then  ! on the fly but for k belonging to kpoint-mesh kpc
!    write(*,*)' entering self_w '
!    write(*,*)' args :k,la,om,temp=',k,la,omega,temp
!     call function_self_w2(k,la,omega,temp,nself,uself)
!! k1
 !   call function_self_w2_sy(q_ibz,k,la,omega,temp,nself,uself)   ! sy
!! k1
  elseif(iread.eq.3)then  ! on the fly but for arbitrary kpoint
     call function_self_w3(k,la,omega,temp,nself,uself)
  elseif(iread.eq.4)then  ! only imaginary part using gaussian for arbitrary kpoint
     call function_self_w4(k,la,omega,temp,nself,uself)
  endif

  end subroutine self_energy
!=====================================================
  subroutine memory_usage
  integer mb,complex_mem,real_mem
  real_mem=8
  complex_mem=16
  mb=1024*1024

!  real_mem*n/mb ! is the memory usage in Mb for a real array of size n

  end subroutine memory_usage

