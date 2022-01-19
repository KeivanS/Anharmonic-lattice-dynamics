!=========================================================
 program thermal_split
! units are eV/A^2,3 for FCs and will be kept throughout except for output files
! where they will be converted into cm^-1; same holds for FC3 (in eV/A^3)
! frequencies as well as perturbation V3, tempearature and self energies are in cm^-1
! g's and kpoints have a factor of 2pi included in them so that exp(ig.r)=1
! this version calculates the contribution of a subset of kpoints in the IBZ to kappa,
! the kpoints over which integration is done span the whole FBZ and can be shifted.
! it also includes the full (iterative) solution of Boltzmann equation
! by Keivan Esfarjani (July 11th 2014)
!
! For iterative solution kpc is only shifted by -0.5*(g1+g2+g3) and kibz is a subset of 
! the kpc mesh, v3sq is written into a file in the form v3sq(ibz,fbz,ndn,ndn,ndn) 
!
! For the RTA case, kpc may or may not be shifted from gamma point, v3sq is stored 
! in the form n1(i),n2(i),n3(i),v3sq(i) i=1,nibz*nkc*ndyn**3  and n1(i)= (nk1-1)*ndyn+l1
! where even for the kibz mesh in ni, the nk1 is the index of kpc mesh
!
! use get_k_info_cent(q) for q in IBZ, and get_k_info_cent(q-shft) for q in FBZ which has a shftxyz 
!
! TODO LIST:
! substitute om=w_q+Delta_q+I Gamma in the argument of the self-energy see how much
! the resutls change. i.e. do the self-consistent calculation, that would
! also get rid of the ambiguity in the eta value.
! check Im(self) against the diagonal matrix elements of the collision matrix (1+n1+n2 versus n1n2/nq)
! write the iteration subroutine (in effect inversion matrix) and solve the
! full BTE, to compare with RTA
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
! output the momentum and energy resolved jdos like the caltech paper
! workout all isothermal elastic and even third order elastic constants


! in parallel to the iterative solution for comparison.
! take the -k dependence on P_12^3 into consideration see if there is any change.
! think about using symmetries to reduce memory and disk space

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
 use exactBTE2
 use mod_ksubset2
 use tetrahedron
! added by k1 on July 17-2014
 use geometry
 use lattice

 implicit none
 integer i,j,mx,la,ncm,ierr,i2,j2,k2,inside,indexg,itmp,igam,jks,al
 real(8), allocatable :: omega(:),sigp(:),total_dos(:),integrate_dos(:),func(:)
 real(8) temp,tempk,dq,qx,qy,qz,kap_ovr_t(3,3),res,errsvd,ermax , ni,nbe
 real(8) ommax,mass_cell,vgmax,mfp,sigma,qq(3),cross_section,wmin,mysqrt,funcmin
 real(8) , allocatable :: kap(:,:,:),vgr(:,:),sig(:)
 real(8) , allocatable :: c_matrix(:,:),evc(:),evc1(:), evc2(:), fw(:),taud(:),evlt(:)
 complex(8) , allocatable :: nself(:),uself(:)
 real cputim
 character now*10,today*8,zone*5,ext*2,wxt*3
! additions by SY
 integer ksubset(2)
 real(8), allocatable :: lg(:)
 integer, allocatable :: mcor(:)
 integer n1,n2,l1,l2,l3

! call get_v3_indices(21600,12,1000,6,6,6,n1,n2,l1,l2,l3)
! write(*,*) n1,n2,l1,l2,l3
! stop
!==================== START OF THE CODE =======================

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(cputim)

  open(ulog,file='phonlog.dat',status='unknown')
  open(utimes,file='times.dat' ,status='unknown')
  open(ugrun,file='mode_grun.dat',status='unknown')
  open(ugrun+1,file='all_grun.dat',status='unknown')
  open(ualph,file='alpha.dat',status='unknown')
  open(udist,file='distr.dat',status='unknown')

 write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' Program THERMALsplit was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(ulog,'(a,f12.4)')' STARTING TIME OF THE PROGRAM IS ',cputim
 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' Program THERMALsplit was launched at ',today(1:4)//'/' &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(utimes,'(a,f12.4)')' At the start of the Program,        TIME IS ',cputim

2 format(a,2x,f19.12,9(2x,g11.5))
3 format(i6,99(1x,g11.5))
4 format(99(2x,2(1x,f7.3)))
5 format(99(1x,f7.3))
6 format(99(2x,g12.6))
7 format(a,99(2x,g11.5))
8 format(2i6,99(1x,g11.5))
9 format(i6,i6,3(1x,f7.3),1x,9(1x,g10.4))

! this is used to convert v33 to cm^_1
  const33 = 1d30*(hbar/uma)**1.5d0/(sqrt(ee*1d20/uma))**1.5d0*ee/h_plank/100/c_light
  tolerance=1d-4    ! to equate two real numbers


  call read_params
  print*,' file params.phon read'

  call read_input_fit
  print*,' file params.inp read'

  allocate(zeu(3,3,natoms0)) !JS
  call read_born
  print*,' file params.born read'

  call read_lattice

  nkc = nc1*nc2*nc3
  nkc2=nkc
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

! now output eigenval is directly frequencies in 1/cm
  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,ndyn,eigenvec_bs,veloc)

  call gruneisen(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,eigenvec_bs,ugrun,grun_bs)
  close(ugrun)

!*************************************************************
! can also calculate the phase space integral along kpoints
!*************************************************************

! do it at 300K  ; if negative, do not include the BE distribution functions
! tempk=300
! call jdos_sy(ndyn,nkp_bs,kp_bs,eigenval_bs,dk_bs,tempk)
!  write(ulog,*) 'phase space lifetime for band structure done...'

   call deallocate_eig_bs
   deallocate(dk_bs)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After band structure,               TIME IS ',cputim


!---------------------------------------------------------------
! calculation of eigenvalues, Gruneisen in the FBZ, and DOS

! first generate the mesh, kpc and set the weights
  nkc = nc1*nc2*nc3
  allocate(kpc(3,nkc),mappos(nkc),wk(nkc))
! k1  call allocate_tetra(nc1,nc2,nc3,wmesh,lamax)
  call allocate_tetra(nc1,nc2,nc3,wmesh)

  call make_kp_reg_tet   ! also sets wk=1/nkc

! this suboutinte produces KIBZ  and WIBZ and writes them in KPOINTS.IBZ
  call get_weights(nc1*nc2*nc3,kpc)

! this is the shift to center to produce the smallest KIBZ set 
  shft= (-0.5d0)*(g1+g2+g3) 
  do i=1,nkc
     kpc(:,i)=kpc(:,i)+shft
  enddo
  do i=1,nibz
     kibz(:,i)=kibz(:,i)+shft
  enddo

! this is the shift of the FBZ kpoints by 0.5 mesh for later integrations
! can not use shifted kpc since veloc should be the same for both sets of kpoints even in the RTA case where iter=0
! if (iter .eq. 0) then  
!    shft=(shftx/nc1)*g1+(shfty/nc2)*g2+(shftz/nc3)*g3
!    write(ulog,7)'small shift vector is:',shft
!    do i=1,nkc
!       kpc(:,i)=kpc(:,i)+shft
!    enddo
! else 
     shft=0
! endif

! sort kpoints and write into file 
  allocate(lg(nkc),mcor(nkc))
  do i=1,nkc
     lg(i)=length(kpc(:,i))
  enddo
  call sort(nkc,lg,mcor,nkc)
  open(fbz,file='KPOINT.FBZ',status='unknown')
! write(fbz,*)'#i,l,mapibz(l),kp(l),length(kp(l)),kibz(mapibz(l)),wibz(mapibz(l)) l=mcor(i)'
  write(fbz,*)'#i,l,mapibz(l),kp(l),length(kp(l)) '
  do i=1,nkc
     j=mcor(i)
     write(fbz,8)i,j,mapibz(j),kpc(:,j),length(kpc(:,j)) !,kibz(:,mapibz(j)),wibz(mapibz(j))
  enddo
  close(fbz)
  deallocate(lg,mcor)

! if mappos(i), (i=1,npos) is the index of the "positive" kpoints which have a negative
! the negative k's are kpc(:,mappos(i)), i=npos+1,nkc
  call get_negatives(nkc,kpc,mappos,npos)

!---------------------------------------------------------------
! allocates eigenval,eigenvec and grun of size (ndyn,nkc) and veloc(3,ndyn,nkc)
! as well as eivalibz and eivecibz(ndyn,nibz)
  call allocate_eig(ndyn,nibz,nkc)   ! 2nd one is for IBZ, 3rd one is for FBZ mesh
  write(ulog,*)' eigenval, eivalibz... allocated, now do dispersion and DOS calculation'
  allocate(dkc(nkc))
  dkc=0
  call get_frequencies(nkc,kpc,dkc,ndyn,eigenval,ndyn,eigenvec,veloc)
  deallocate(dkc)  !,eivalibz,eivecibz)

  allocate(dkc(nibz)) !,eivalibz(ndyn,nibz),eivecibz(ndyn,ndyn,nibz),velocibz(3,ndyn,nibz))
  dkc=0
  call get_frequencies(nibz,kibz,dkc,ndyn,eivalibz,ndyn,eivecibz,velocibz)
    
  open(uvel,file='veloc.dat')
  open(uvel+1,file='velocibz.dat')
  write(uvel,*)'# la,nk,kp(nk),om(la,nk),vg(la,nk),length(vg(la,nk))'
  do la=1,ndyn
  do i=1,nkc
write(uvel,8)la,i,kpc(:,i),eigenval(la,i),veloc(:,la,i)*c_light,length(veloc(:,la,i))*c_light
  enddo
  do i=1,nibz
write(uvel+1,8)la,i,kibz(:,i),eivalibz(la,i),velocibz(:,la,i)*c_light,length(velocibz(:,la,i))*c_light
  enddo
  enddo
  close(uvel)
  close(uvel+1)

! gruneisen in the volume is needed for the thermal expansion coefficient
  call gruneisen(nkc,kpc,dkc,ndyn,eigenval,eigenvec,ugrun+1,grun)
  close(ugrun+1)
  deallocate(dkc)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After om(k) calculation within FBZ, TIME IS ',cputim

! replace eivec(-q) by conjg(eivec(q))
  call subst_eivecs(ndyn,nkc,eigenval,eigenvec,kpc,npos,mappos)

! check to see if D(-k)=conjg(D(k))
!    call check_mdyn(ndyn,nkc,kpc,eigenvec,eigenval) ! would only work with nkc
!    call write_eivecs(ndyn,nkc,kpc,eigenvec)

  call flush()

!----------------------------------------
! calculate dos (added over all bands)
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
     write(wxt,'(i2.2)')la
     dos(0,:) = dos(0,:) + dos(la,:)
  enddo
  call write_dos
  close(udos) ; deallocate(evc)
  write(*,*) 'gaussian dos done'

! calculate dos using tetrahedron method
  open(udos,file='dos_tet.dat',status='unknown')
  open(udos+1,file='dos2_tet.dat',status='unknown')
  write(udos,*)'#la,om(i),dos(la,i),total_dos(i),sum(dos_tet(1:i))*dw'  
  write(udos+1,*)'#junk,om(i),total_dos(i),integrate_dos(i),tretrados(test)'
  allocate(integrate_dos(wmesh),total_dos(wmesh),func(nkc))
  wmin=0d0
  total_dos=0
  allocate(evc(nkc))
  do la=1,ndyn
     evc=eigenval(la,:)
     call calc_tet(wmesh,wmin,wmax,nkc,om,kpc,evc,evc)
     integrate_dos=0
     do i=1,wmesh
        total_dos(i)=total_dos(i)+dos_tet(i)   ! this is the sum over all bands
        integrate_dos(i)=integrate_dos(i)+sum(total_dos(1:i))*(wmax-wmin)/wmesh  ! this is the cumulative dos
        write(udos,3)la,om(i),dos_tet(i),total_dos(i),sum(dos_tet(1:i))*(wmax-wmin)/wmesh 
     enddo
  enddo

  allocate(lg(nkc))
  func=1
  do i=1,wmesh
     evc=eigenval(1,:)
     call tet_sum(om(i),nkc,evc,func,res,lg)
     write(udos+1,3)la,om(i),total_dos(i),integrate_dos(i),res
  enddo
  deallocate(integrate_dos,total_dos,evc,lg)
  close(udos); close(udos+1)
  write(*   ,*) 'tetrahedron dos done'
  write(ulog,*) 'tetrahedron dos done'
  deallocate(dos)

  call read_ksubset(ksubset,nibz)
  if(readv3.eq.0) then
     call phase_space_lifetime(ksubset(1),ksubset(2),nibz,kibz,eivalibz)
     write(*   ,*) 'phase space lifetime done'
     write(ulog,*) 'phase space lifetime done'
  endif

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After phase_space_lifetime,       TIME IS ',cputim

! do it at 300K  ; if negative, do not include the BE distribution functions
!     tempk=-300
!     call jdos_sy(ndyn,nkc,kpc,eigenval,tempk)
!     write(ulog,*) 'phase space lifetime for band structure done...'
!
!call cpu_time(cputim)  !date_and_time(time=tim)
!write(utimes,'(a,f12.4)')' After jdos_sy,                    TIME IS ',cputim
!
!----------------------------------------
! calculate the thermal expansion coefficient using phonons and gruneisen
!    call read_grun(nkc,ndyn,kpc,grun)

  call gruneisen_fc

! call calculate_thermal(nibz,wibz,ndyn,eivalibz,grun,tmin,tmax,ualph,velocibz)
  call calculate_thermal(nkc,wk,ndyn,eigenval,grun,tmin,tmax,ualph,veloc)

  close(ualph)
  deallocate(grun)

  write(*   ,*) 'calculate_thermal done'

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After calculate_thermal           TIME IS ',cputim


  if (dos_bs_only ) stop             ! sy

! test to see whether the eivecs at opposite kpoints are now conjugate
  call write_eivecs(ndyn,nkc,kpc,eigenvec)

!========================================================================
! the anharmonic V3 matrix will be calculated on the k-mesh

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' Start V33 read or calculation,      TIME IS ',cputim

  nv3split=nkc*(ksubset(2)-ksubset(1)+1)*ndyn**3
  write(ulog,*) 'nkc, ksubset, nv3split',nkc,ksubset(1),ksubset(2),nv3split

! calculate and write the matrix elements v33^2 in split mode from scratch, do RTA calculation
  if ( readv3 .eq. 0 ) then   

       write(*,*) 'entering v33 calculation in split mode'   
       write(ulog,*) 'MAIN: calculating v3_split...'
!      if (iter.eq.1) then  ! full iterative needs v3 with 5 arguments: v35
          allocate(v3sq(ksubset(2)-ksubset(1)+1,nkc,ndyn,ndyn,ndyn))
          call calculate_v35(ksubset,ndyn,nibz,eivalibz,eivecibz,nkc,eigenval,eigenvec) 
          write(*   ,*) 'calculate_v3 done'
          call flush()
!      else   ! for RTA, one can store v3sq in one single array
!         call allocate_v33(nv3split)
!         call calculate_w3_ibz_split_sq(ksubset,nv3split,ndyn,nkc,kpc)
!         write(*   ,*) 'calculate_w3_ibz_split done'
!      endif

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' End of V33 calculation,        TIME IS ',cputim

       tempk=tmin
       temp = tempk*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1

       call rta_sub(ksubset,tempk)  !nibz,kibz,eivalibz,velocibz,temp)

       ! for iterative case, calculate and write the collision matrix at T=temp into files
       if (iter.eq.1) then
          call calculate_sub_collision(temp,ksubset)
       endif
    
       deallocate(v3sq)

! read v33 matrix elements from split files and for iter=0 do RTA calculations for the ksubset
! cwfor iter=1 calculate and write collision matrix elements for the ksubset 
  elseif (readv3 .eq. 1) then   
       write(ulog,*) 'MAIN: reading v3_split from file...'
!      if (iter.eq.1) then ! full iterative needs v3 with 5 arguments: v35
         allocate(v3sq(ksubset(2)-ksubset(1)+1,nkc,ndyn,ndyn,ndyn))
         call read_v35(ksubset,ndyn)
         write(*   ,*) 'read_v3 done'
!      elseif(iter.eq.0) then
!         call allocate_v33(nv3split)
!         call read_v3_new_unformatted_k1(ksubset,nv3split,nkc)
!         write(*   ,*) 'read_v3_new_unformatted_k1 done'
!      endif

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' End V33 read ,        TIME IS ',cputim
 
     temploop: do i=1,ntemp
       tempk=tmin+(i-1)*(tmax-tmin)/(ntemp-1+1d-9)
       temp = tempk*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1

       ! do an RTA calculation no matter iter=0 or not
       call rta_sub(ksubset,tempk)  !nibz,kibz,eivalibz,velocibz,tmin,tmax,ntemp)
       call flush()

       ! for iterative case, calculate and write the collision matrix at T=temp into files
       if (iter.eq.1) then
          call calculate_sub_collision(temp,ksubset)
       endif
     enddo temploop

       deallocate(v3sq)

! read and collect the collision matrix parts to solve BTE and calculate kappa
  elseif (readv3 .eq. 2) then    

       if(iter.ne.1) then
          write(ulog,*) 'MAIN: readv3=2 and iter.ne.1 are not compatible, check your input' 
          stop
       endif

! dimensions of the collision matrix
       nlin=ndyn*nkc
       ncol=ndyn*nibz

! v33 is not needed, only collision matrix elts are read from different files and put together
       call read_collision(temp) 

! put everything inside mval
  !*** do i=1,nibz
  !***    mval(i,i)=mval(i,i)+qval(i)
  !*** enddo
       if( allocated(mval) ) deallocate(mval) 
       allocate(mval(nibz*ndyn,nibz*ndyn))
       mval=0
       do i=1,nibz*ndyn
          la=nband(i,ndyn)
           j=nkpt(i,ndyn)
           ni=nbe(eivalibz(la,j),temp,classical)
           mval(i,i)=qval(i)     ! that's also RTA!
       enddo

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' collision matrix built,        TIME IS ',cputim

       allocate(kappa(ndyn,3,3),kappa_rta(ndyn,3,3))
       call allocate_iter1(nibz,nkc,ndyn) ! allocates rhs,rhsi,coll,ne_dist
!      call initialize_phi(temp)
       call initialize_phi(temp,nibz,ndyn,eivalibz,velocibz,qval,ne_dist)
       write(udist,*)'###### initialized distribution ######'
       call write_dist(ncol,ne_dist,udist)
       call calculate_kappa_iter(nibz,ndyn,wibz,temp,ne_dist_rta,eivalibz,velocibz,kappa_rta) 
       open(ukap,file='kap_RTA.dat')
          call write_kappa(ndyn,kappa_rta,ukap)
       close(ukap)


       call right_hand_side(nlin,rhs,temp)    ! in veloc/temp i.e. in cm

! solve mval*ne_dost=rhs
       if (isvd.eq.1) then
   
         allocate(sig(ncol))
!        call svd_set(nlin,ncol,mval,rhs,ne_dist,sig,svdcut,errsvd,ermax)
         call svd_set(ncol,ncol,mval,rhs(1:ncol,:),ne_dist,sig,svdcut,errsvd,ermax)

         deallocate (mval,qvalue_N,qvalue_U,qval)

call cpu_time(cputim)  
write(utimes,'(a,f12.4)')' SVD done,        TIME IS ',cputim

       else

         call rearrange_into_ibz(ncol,nlin,mval,qval,rhs,coll,rhsi)  ! coll is in 1/cm

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' collision put into reduced form,        TIME IS ',cputim

! initialize using RTA distribution function: tau*v n(1+n) hbarw/kT^2
!        call initialize_phi(temp)
         call initialize_phi(temp,nibz,ndyn,eivalibz,velocibz,qval,ne_dist)
!        call random_number(ne_dist)  ! alternatively can initialize with random numbers
        
         deallocate(rhs,qvalue_N,qvalue_U,mval)  ! qval is needed for initialization
 
         write(udist,*)'###### initialized distribution ######'
         call write_dist(ncol,ne_dist,udist)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' phi initialized before CG,        TIME IS ',cputim

         call precondition_all(ncol,coll,rhsi,ne_dist,qval)

         do al=1,3
            call cg_quadratic(ncol,coll,rhsi(:,al),ne_dist(:,al),funcmin,max_iter,n_dig_acc)
         enddo

         call switch_back(ncol,ne_dist,qval)

       endif

       write(udist,*)'####### final distribution ######'
       call write_dist(ncol,ne_dist,udist)

       call calculate_kappa_iter(nibz,ndyn,wibz,temp,ne_dist    ,eivalibz,velocibz,kappa) ! kappa(ndyn,3,3)
       open(ukap,file='kap_iter.dat')
       call write_kappa(ndyn,kappa,ukap)
       close(ukap)

call cpu_time(cputim)  !date_and_time(time=tim)
write(utimes,'(a,f12.4)')' After full BTE calculation,         TIME IS ',cputim

  endif

! find the hight-T limit (coefficient of 1/T in kappa): does not depend on T
  if(iter.eq.0) then  ! do it within the RTA
!        call kap_over_T(ksubset,nkc,kpc,ndyn,veloc,kap_ovr_t)

  endif

!           if (readv3.eq.2) then
!              call read_all_v3sq
!           endif

!           temp=tmin*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1
!           call put_collision_together(temp)  ! allocates and fills up Qvalue_N,U,qval,mval


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

 close(ulog)
 close(utimes)
 close(udist)


 end program thermal_split
!===============================================================

