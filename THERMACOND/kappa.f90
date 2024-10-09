!=========================================================
 program kappa_full
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
 use ios
 use kpoints
 use om_dos
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
 use mpi

 implicit none
 integer i,j,mx,la,ncm,ierr,i2,j2,k2,inside,indexg,itmp,igam,ig
 real(8), allocatable :: dkc(:),omega(:),sigp(:)
 real(8) , allocatable :: tau_klem_inv(:,:),kap(:,:,:),vgr(:,:)
 real(8) , allocatable :: c_matrix(:,:),inv_relax_t(:),evc(:),evc1(:), evc2(:), fw(:),taud(:),evlt(:)
 real(8), allocatable :: tau_inv_n(:,:), tau_inv_u(:,:), omega0(:),tse(:)!,v33s8(:)
 complex(8) , allocatable :: nself(:),uself(:),self(:,:,:,:),evct(:,:)
 real cputim
 character now*10,today*8,zone*5,ext*2,wxt*3

 integer ksubset(2), nv3_split, ksubset_bs(2)
 integer mergemode, num_files               
 integer convergence,iter_cont 
 integer safo,nss,nsl,nxyz,nn, narms,kvecop(48)
 integer i1,j1,k1,vs,vsdf
 integer ss,si,sa,ms,sd,kvecopsd(48)
 integer ds, kvecops(48)
 integer ni,i3,j3,k3,nz,ni0,sz,ssz,szz
 integer rs,klrs,col1,col2
 integer colkil,coll1,coll2
 real(8) temp,tempk,dq,qx,qy,qz,kappa_klem,kappa_q(3,3),q0(3),q1(3),om0,om1,vel0
 real(8) fs(3),qs(3), qsi(3),ve2(3), we2(3)
 real(8) kvecstars(3,48),qss(3),difvs
 real(8) nqsi(3), fsRTA(3),kvecstar(3,48),kvecstarsd(3,48)
 real(8) ommax,mass_cell,vgmax,mfp,sigma,qq(3),cross_section,boundary_rate
 complex(8) selfw
 integer nv3o
 real(8) kbss(3),ise 
 integer nkbss
 complex(8) nselfs,uselfs
 real(8) , allocatable ::evlbs(:,:),vlcbs(:,:,:)
 complex(8) , allocatable ::evcbs(:,:,:)
 integer kapRTA,tauf,kapd,tauef,tauefi,kappai,lwse,ll,lwse2
 real(8) normk,fm
 integer is,js,ks,s

!=========MPI=======================
real(8), allocatable :: P1sm_temp_mpi(:),P1sm_tempmpi(:)
real(8), allocatable :: P2sm_temp_mpi(:),P2sm_tempmpi(:)
!***real(8), allocatable :: p1sm_mpi(:,:,:,:,:), p2sm_mpi(:,:,:,:,:)
real(8), allocatable :: ise_temp_mpi(:),ise_tempmpi(:)
real(8), allocatable :: tse_temp_mpi(:),tse_tempmpi(:)
real(8), allocatable :: CM_temp_mpi(:),CM_tempmpi(:)

real(8), allocatable :: Qvalue_tempi(:),Qvaluetemp(:)
real(8), allocatable :: Qvalue_Ntempi(:),QvalueNtemp(:)
real(8), allocatable :: Qvalue_Utempi(:),QvalueUtemp(:)


integer, allocatable :: aksub_mpi(:), displs_mpi(:)
integer, allocatable :: caksub_mpi(:), cdispls_mpi(:)
integer            :: rank_mpi, nprocs_mpi
integer            :: errcode_mpi, ierr_mpi
integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
integer, parameter :: root_mpi=0, tag_mpi=0
integer            :: local_num_k_mpi
integer            :: values1_mpi, values2_mpi, ksub_mpi, smpi
integer            :: num_k_mpi, num_k_mpi2
integer            :: n1_mpi, n2_mpi, l1_mpi, l2_mpi, l3_mpi, n_mpi, k_mpi
integer            :: ncount1_mpi, ncount2_mpi,cncount1_mpi
logical            :: debugmpi
integer            :: tt,tto, sts_mpi, i_rankmpi,tetttst,ksubmpiv,sbs
character(99) filename_mpi, filename_mpii
integer            ::cmkn1,cmkn1t,cmkn2,cksubmpiv
character(len=100) :: filenamempi
double precision :: start_timempi, end_timempi, elapsed_timempi
integer          :: unit_numpi
integer          :: CM_1D, DS_1D
!===================================


call date_and_time(date=today,time=now,zone=zone)
call cpu_time(cputim)

  open(ulog,file='phonlog.dat',status='unknown')
  open(utimes,file='times.dat' ,status='unknown')
  open(ueiv,file='bs_freqs.dat',status='unknown')
  open(ucors,file='FBZ_freqs.dat',status='unknown')
  open(ugrun,file='bs_grun.dat',status='unknown')
  open(ugrun+1,file='FBZ_grun.dat',status='unknown')


 write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' Program kap8 was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(ulog,'(a,f12.4)')' STARTING TIME OF THE PROGRAM IS ',cputim
 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' Program kap8 was launched at ',today(1:4)//'/' &
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

20 format(2i6,2x,99(1x,f9.3))

! this is used to convert v33 to cm^_1
  const33 = 1d30*(hbar/uma)**1.5d0/(sqrt(ee*1d20/uma))**1.5d0*ee/h_plank/100/c_light


  call read_params
  print*,' file latdyn.params read'


  call read_input_fit
  print*,' file structure.params read'


  call read_lattice  ! reads lat_fc.dat


  call read_born
  print*,' file dielectric.params read'
 

  nkc = nc(1)*nc(2)*nc(3)
  dq = (volume_g /nkc)**(1d0/3)
  write(ulog,*)' nkc, dq(ang^-1)=', nkc,dq


  call read_fc234


  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After readings,                     TIME IS ',cputim


  ndyn = 3*natoms0

!---------------------------------------------------------------
  if (job .ge. 0)  then
!---------------------------------------------------------------
! phonon band structure & gruneisen params along the symmetry directions

  call make_kp_bs
  write(ulog,*)' Kpoints for band structure generated from kpbs.params'
  

  call allocate_eig_bs(ndyn,nkp_bs,ndyn) !3rd index is for eivecs needed for gruneisen
  write(*,*)' Eigenvalues allocated'


  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,ndyn,eigenvec_bs,ueiv,veloc)
  write(*,*)' Frequencies calculated'


  open(1008,file='nkpbs.dat')
       write(1008,*)nkp_bs
  close(1008)
  

  call gruneisen(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,eigenvec_bs,ugrun,grun_bs)
  write(*,*)' Gruneisen calculated'


! first generate the COARSE mesh, kpc and set the weights
  nkc = nc(1)*nc(2)*nc(3)
  allocate(kpc(3,nkc),mappos(nkc),wk(nkc))


! Need kpc for phase space integration 
  call make_kp_reg(nc,g1,g2,g3,shft,kpc,wk)        ! original one.
  write(*,*)' make_kp_reg called '


! call phase_space_lifetime_bs
  write(*,*)' Phase_space calculated; now deallocating'
 

  call deallocate_eig_bs
  call deallocate_kp_bs


  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After band structure,               TIME IS ',cputim

!---------------------------------------------------------------
! calculation of eigenvalues, Gruneisen in the FBZ, and DOS using the coarse mesh

  call get_negatives(nkc,kpc,mappos,npos)
  write(*,*)' get_negatives called '


! define the Irreducible BZ based on symmetry; output is nibz,kibz
  call get_weights(nc(1)*nc(2)*nc(3),kpc)
  write(*,*)' get_weights called '

!---------------------------------------------------------------
 
! call allocate_eig(ndyn,nibz,nkc)   ! last one (nkc) is for veloc
  call allocate_eig(ndyn,nkc,nkc)   ! last one (nkc) is for veloc
 

  allocate(dkc(nkc))
  dkc=0
  write(ulog,*)' eigenval allocated, entering the loop for DOS calculation'


  call get_frequencies(nkc,kpc,dkc,ndyn,eigenval,ndyn,eigenvec,ucors,veloc)


! gruneisen in the volume is needed for the thermal expansion coefficient
! call gruneisen(nkc,kpc,dkc,ndyn,eigenval,eigenvec,ugrun+1,grun)
! call gruneisen(nibz,kibz,dkc,ndyn,eigenval,eigenvec,ugrun,grun)


  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After om(k) calculation within FBZ, TIME IS ',cputim


! check to see if D(-k)=conjg(D(k))
  call check_mdyn(ndyn,nkc,kpc,eigenvec,eigenval) ! would not work with nibz
  write(*,*)' check_mdyn called'


  if(verbose) call write_eivecs(ndyn,nkc,kpc,eigenvec)


  call subst_eivecs(ndyn,nkc,eigenvec,kpc,npos,mappos)
  write(*,*)' subst_eivec called'


! calculate dos and joint-dos (added over all bands)
 
  call set_omdos(ndyn,wmesh)  ! allocates om and dos arrays of size mesh
  write(*,*)' set_omdos called'
 
  allocate(evc(nibz))
 
  dos = 0
  do la=1,ndyn
     j=0
     do i=1,nibz
        j=j+1
        evc(j)=eigenval(la,mapinv(i))
     enddo
     call calculate_dos_g(nibz,evc,wibz,wmesh,om,dos(la,:))  ! with gaussian broadening 
     dos(ndyn+1,:) = dos(ndyn+1,:) + dos(la,:)
  enddo


  call write_dos
  deallocate(evc)
  write(*,*) 'gaussian dos done'


! calculate dos using tetrahedron method
  call allocate_tetra(nc(1),nc(2),nc(3),wmesh,ndyn)
  call make_kp_reg_tet(nc(1),nc(2),nc(3),shft(1),shft(2),shft(3),kpc,wk)
  call calc_tet(nc(1),nc(2),nc(3),ndyn,kpc)
  write(*,*) 'tetrahedron dos done'
  deallocate(dos,dkc)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After calculate tetrahedron dos    TIME IS ',cputim

! call phase_space_lifetime_FBZ
  write(*,*) 'phase space dos done'


  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After calculation of phase space in FBZ,     TIME IS ',cputim

  deallocate(grun)
! call deallocate_tetra
    
  if (dos_bs_only ) stop  
  write(ulog,*) ' DOS_BS finished here '


!----------The anharmonic V3 matrix calculation in IBZ and kappa RTA and Direct using MPI--------------------------------
  nv3 = nkc*nibz*ndyn**3
!-------------------------------------------------------
  if ( job .eq. 1 ) then   ! MPI 
!-------------------------------------------------------
  write(ulog,*) 'entering cal v33 in IBZ split...'
  write(*,*) 'entering cal v33 in IBZ split...'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')'V33 calculation start,                TIME IS', cputim
       
 !Initialize MPI, get the local number of k_point
  call MPI_INIT(ierr_mpi)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs_mpi,ierr_mpi)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank_mpi,ierr_mpi)

! Create a unique filename for each rank
! write(filenamempi, '("timing_rank_", I0.2, ".txt")') rank_mpi

! Open the file for writing
! unit_numpi = 10 + rank_mpi
! open(unit=unit_numpi, file=filenamempi, status='replace', action='write')

  if (rank_mpi==0) then

      num_k_mpi=nibz

  endif


! At this point, only root knows these values.  Send them out.
  call MPI_Bcast(num_k_mpi,1,MPI_INTEGER,root_mpi,MPI_COMM_WORLD,ierr_mpi)

      local_num_k_mpi = num_k_mpi/nprocs_mpi
      write(*,*) 'num_k_mpi, local_num_k_mpi, nprocs_mpi',num_k_mpi, local_num_k_mpi, nprocs_mpi

! ksubset(1) and ksubset(2) are needed for the loops in the follwoing subroutines (calculate_w3_ibz_split_sq,calculate_FGR_mpi)
  ksubset(1) = rank_mpi * local_num_k_mpi + 1
  ksubset(2) = rank_mpi * local_num_k_mpi + local_num_k_mpi
       
  if (mod(num_k_mpi,nprocs_mpi).ne.0) then
     if  (rank_mpi.eq.(nprocs_mpi-1)) then
         ksubset(2) = rank_mpi * local_num_k_mpi + local_num_k_mpi + mod(num_k_mpi,nprocs_mpi)
     endif
  endif  
! write(*,*)'rank_mpi,ksubset1,ksubset2',rank_mpi,ksubset(1),ksubset(2)

  nv3_split=nkc*(ksubset(2)-ksubset(1)+1)*ndyn**3        
  call allocate_v33sq(nv3_split)

  allocate(v33s8((ksubset(2)-ksubset(1)+1),nkc,ndyn,ndyn,ndyn))
! call calculate_w3_fbz_split_sq(ksubset,nv3_split,ndyn,nkc,kpc,eigenval,eigenvec)
! start_timempi = MPI_Wtime()
  call calculate_w3_ibz_split_sq (ksubset,nv3_split,ndyn,nkc,kpc,eigenval,eigenvec,nv3o)
  write(*,*)'after v3'

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-v33: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately
 

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After v33 calculation,                TIME IS ',cputim

  ksub_mpi = (ksubset(2) - ksubset(1)) + 1

! write(*,*)'rank_mpi,ksubset1,ksubset2,ksub_mpi',rank_mpi,ksubset(1),ksubset(2), ksub_mpi
! writing v33 for each rank
! sts_mpi=1
! tt = 1333

! write(filename_mpi,fmt="(a,i3.3,a)") "v33.",rank_mpi,".dat"
! open(tt,file=filename_mpi,status='replace')

! do n1_mpi=1,ksub_mpi
!    do n2_mpi=1,nkc
!       do l1_mpi =1,ndyn
!          do l2_mpi=1,ndyn
!             do l3_mpi=1,ndyn

!                 write(tt,*) rank_mpi, n1_mpi, n2_mpi, l1_mpi, l2_mpi, l3_mpi, v33s8(n1_mpi,n2_mpi,l1_mpi,l2_mpi,l3_mpi)
!              enddo
!           enddo
!        enddo
!     enddo
! enddo

! close(tt)

  allocate(P1smpi((ksubset(2)-ksubset(1)+1),nkc,ndyn,ndyn,ndyn))
  allocate(P2smpi((ksubset(2)-ksubset(1)+1),nkc,ndyn,ndyn,ndyn))
        
  write(ulog,*) 'entering cal v33*delta...'
  write(*,*) 'entering cal v33*delta...'
 
! start_timempi = MPI_Wtime()
! changed to IBZ  call calculate_v3sq_delta (ksubset,ndyn,nkc,kpc)  
  call calculate_FGR_mpi2(ksubset,ndyn,nv3o)!,v33s8)
  write(*,*)'after FGR'
! call deallocate_iter

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-v33delta: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After v33_delta calculation,            TIME IS ',cputim

  deallocate(v33s8)
! call deallocate_v33sq

! writng P1sm and P2sm for each rank
! sts_mpi=1
! tt = 1333
! write(filename_mpi,fmt="(a,i3.3,a)") "P1smpi.",rank_mpi,".dat"
! write(filename_mpii,fmt="(a,i3.3,a)") "P2smpi.",rank_mpi,".dat"
! open(tt,file=filename_mpi,status='replace')
! open(tt+1,file=filename_mpii,status='replace')

! do n1_mpi=1,ksub_mpi
!    do n2_mpi=1,nkc
!       do l1_mpi =1,ndyn
!          do l2_mpi=1,ndyn
!             do l3_mpi=1,ndyn
!                 write(tt,*) rank_mpi, n1_mpi, n2_mpi, l1_mpi, l2_mpi, l3_mpi,P1smpi(n1_mpi,n2_mpi,l1_mpi,l2_mpi,l3_mpi)
!                 write(tt+1,*) rank_mpi, n1_mpi, n2_mpi, l1_mpi, l2_mpi, l3_mpi,P2smpi(n1_mpi,n2_mpi,l1_mpi,l2_mpi,l3_mpi)
!              enddo
!           enddo
!        enddo
!     enddo
! enddo

! close(tt)

!----------------------------------------------------------
  if (iso.eq.1) then     !Piso calculation in IBZ split 
!-----------------------------------------------------------
  write(ulog,*) 'entering Piso calculation in IBZ split'
  write(*,*) 'entering Piso calculation in IBZ split'
  write(ulog,*) 'job', job

  allocate (dist(nkc,ndyn),frequency(nkc,ndyn))
  allocate(Pisos((ksubset(2)-ksubset(1)+1),ndyn))
  tempk = tmin     ! use tmin in params.phon as temperature. tmax is not being used.
  call distribution(nkc,ndyn,tempk)  ! calculate BE distribution function for each state

  write(*,*)'natoms0', natoms0

  allocate(giso(natoms0))
  call cal_giso(giso)
  call P_iso_split(ksubset,giso,ndyn,eigenvec)
  deallocate(dist,frequency)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After Piso calculation,            TIME IS ',cputim
  write(ulog,*) 'Piso calculation done'
  write(ulog,*) 'Piso calculation done'
  endif
!------------------------------------------------------------------------------------------

!------------- Exact solution of BTE--------------------------------------------------------
  write(ulog,*) 'entering to the direct and iterative BTE solver'
  write(*,*) 'entering to the direct and iterative BTE solver'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' entering to the direct and iterative BTE solver,      TIME IS',cputim

! call allocate_iter (nibz,nkc,ndyn)    ! allocate memory
  call allocate_FGR (ksub_mpi,nkc,ndyn,ndyn,ndyn)
  write(ulog,*) 'allocate_iter done'

  tauf=3450
  kapRTA=6450
  open(tauf,file='tau_RTA.dat',status='unknown')
  open(kapRTA,file='kappa_RTA.dat',status='unknown')

  tempk = tmin     ! use tmin in params.phon as temperature. tmax is not being used.

  call distribution(nkc,ndyn,tempk)  ! calculate BE distribution function for each state
  write(ulog,*) 'distribution done'

  allocate(Qvaluempi((ksubset(2)-ksubset(1)+1),ndyn))
  allocate(Qvalue_Nmpi((ksubset(2)-ksubset(1)+1),ndyn))
  allocate(Qvalue_Umpi((ksubset(2)-ksubset(1)+1),ndyn))

!--------------------------------------------
  if (iso.eq.1) then   !reaind Piso if iso=true
!---------------------------------------------  
  call read_Piso(ndyn,Piso) 
  write(ulog,*) 'reading Piso done'
  call calculate_tauRTA_iso(ndyn,tempk)
  else
! start_timempi = MPI_Wtime()
  call calculate_tauRTA_mpi2(ksubset,ndyn,tempk)!,rank_mpi)        
! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-Q: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately
  endif

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Q and tau_RTA in IBZ done,    TIME IS',cputim
  write(ulog,*) 'Q and tau_RTA in IBZ done'


  deallocate(P1smpi,P2smpi)!,Qvalue_N,Qvalue_U,Qvalue_Nmpi,Qvalue_Umpi)
  deallocate(Piso)
  call deallocate_v33sq

  allocate(Qvalue_tempi((ksubset(2)-ksubset(1)+1)*ndyn))
  allocate(Qvalue_Ntempi((ksubset(2)-ksubset(1)+1)*ndyn))
  allocate(Qvalue_Utempi((ksubset(2)-ksubset(1)+1)*ndyn))

  Qvalue_tempi=0
  Qvalue_Ntempi=0
  Qvalue_Ntempi=0
  
  smpi=0 
  do n1_mpi=1,ksub_mpi
     do l1_mpi=1,ndyn
      
        smpi=smpi+1
        Qvalue_tempi(smpi) = Qvaluempi(n1_mpi,l1_mpi)
        Qvalue_Ntempi(smpi) = Qvalue_Nmpi(n1_mpi,l1_mpi)
        Qvalue_Utempi(smpi) = Qvalue_Umpi(n1_mpi,l1_mpi)

    enddo
  enddo

  deallocate(Qvaluempi,Qvalue_Nmpi,Qvalue_Umpi)
! Allocate array to gather ksubs from all processes
  allocate(aksub_mpi(nprocs_mpi))
  allocate(displs_mpi(nprocs_mpi))

  ksubmpiv=ksub_mpi*ndyn

! Gather ksubs from all processes
  call MPI_Allgather(ksubmpiv, 1, MPI_INTEGER, aksub_mpi, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  displs_mpi(1) = 0
  do i = 2, nprocs_mpi
     displs_mpi(i) = displs_mpi(i-1) + aksub_mpi(i-1)
!    write(*,*)'displs_mpi(i-1),displs_mpi(i),aksub_mpi(i-1),aksub_mpi(i),',i, displs_mpi(i-1),displs_mpi(i),aksub_mpi(i-1),aksub_mpi(i)
  end do

  allocate(Qvaluetemp(num_k_mpi*ndyn))
  allocate(QvalueNtemp(num_k_mpi*ndyn))
  allocate(QvalueUtemp(num_k_mpi*ndyn))
  ncount1_mpi=size(Qvaluempi)

! start_timempi = MPI_Wtime()

  call MPI_Allgatherv(Qvalue_tempi,ncount1_mpi,MPI_DOUBLE_PRECISION,Qvaluetemp,aksub_mpi,displs_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)
  call MPI_Allgatherv(Qvalue_Ntempi,ncount1_mpi,MPI_DOUBLE_PRECISION,QvalueNtemp,aksub_mpi,displs_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)
  call MPI_Allgatherv(Qvalue_Utempi,ncount1_mpi,MPI_DOUBLE_PRECISION,QvalueUtemp,aksub_mpi,displs_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-gatherv_RTA: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately


  smpi=0
  do n1_mpi=1,num_k_mpi
     do l1_mpi=1,ndyn

        smpi=smpi+1
        Qvalue(n1_mpi,l1_mpi) = Qvaluetemp(smpi)
        Qvalue_N(n1_mpi,l1_mpi) = QvalueNtemp(smpi)
        Qvalue_U(n1_mpi,l1_mpi) = QvalueUtemp(smpi)
!       write(*,*)n1_mpi,l1_mpi,Qvalue(n1_mpi,l1_mpi),smpi,Qvaluetemp(smpi)
     enddo
  enddo


  deallocate(Qvalue_tempi,Qvaluetemp,aksub_mpi,displs_mpi)
  deallocate(Qvalue_Ntempi,QvalueNtemp,Qvalue_Utempi,QvalueUtemp)

  call Qandtau(ndyn,tempk)

  allocate(FFF(nibz,ndyn,3))
! start_timempi = MPI_Wtime()
! calculates F_RTA(kibz) 
  call calculate_RTA_distributionfunction(nibz,nkc,kibz,ndyn,tempk,veloc,eigenval)
  write(ulog,*) 'FRTA done'

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' FRTA in IBZ done,    TIME IS',cputim

  call mode_kappa(ndyn,tempk,veloc,eigenval)
  call write_kappa(ndyn,tempk,veloc,eigenval,kapRTA,tauf)

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-FRTA and kappa_RTA: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' kappa_RTA calculation done,   TIME IS',cputim

  close(tauf)
  close(kapRTA)

  deallocate(Qvalue)

  write(ulog,*) 'Start direct approach using cg method'
  write(*,*) 'Start direct approach using cg method'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Start direct approach using cg method,                TIME IS',cputim

  tauef=7450
  kapd=8450
  open(tauef,file='tau_RTAandeff-cg.dat',status='unknown')
  open(kapd,file='kappa_direct-cg.dat',status='unknown')

  ksub_mpi = (ksubset(2) - ksubset(1)) + 1

  allocate(Collision_Matrixmpi(ksub_mpi*ndyn*3,nibz*ndyn*3))

!---------------------------------
  if (iso.eq.1) then !for iso scattering
!---------------------------------
      Piso=0
      call read_Piso(ndyn,Piso)
      write(ulog,*) 'reading Piso done'
      call CollisionMatrix_iso(ndyn)
  else
!     start_timempi = MPI_Wtime()
      call CollisionMatrix_mpi(ksubset,ndyn)
!     end_timempi = MPI_Wtime()
!     elapsed_timempi = end_timempi - start_timempi
!     write(unit_numpi, '(A, F6.2)') "Time after subroutine-CM: ", elapsed_timempi
!     flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately

  endif

! Allocate array to gather ksubs from all processes
  allocate(caksub_mpi(nprocs_mpi))
  allocate(cdispls_mpi(nprocs_mpi))    

  cksubmpiv=ksub_mpi*nibz*ndyn*3*ndyn*3

! Gather ksubs from all processes
  call MPI_Allgather(cksubmpiv, 1, MPI_INTEGER, caksub_mpi, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)


  cdispls_mpi(1) = 0
  do i = 2, nprocs_mpi
     cdispls_mpi(i) = cdispls_mpi(i-1) + caksub_mpi(i-1)
!    write(*,*)'cdispls_mpi(i-1),cdispls_mpi(i),caksub_mpi(i-1),caksub_mpi(i),',i, cdispls_mpi(i-1),cdispls_mpi(i),caksub_mpi(i-1),caksub_mpi(i)
  end do

! write(*,*) displs_mpi,aksub_mpi

  allocate(CM_temp_mpi(ksub_mpi*nibz*ndyn*3*ndyn*3))

  if (mod(num_k_mpi,nprocs_mpi).ne.0) then
      num_k_mpi2 = local_num_k_mpi * nprocs_mpi + mod(num_k_mpi,nprocs_mpi)
  else 
      num_k_mpi2= num_k_mpi
  endif

! write(*,*)num_k_mpi2,num_k_mpi,local_num_k_mpi, nprocs_mpi, mod(num_k_mpi,nprocs_mpi)

! allocate(CM_tempmpi(num_k_mpi2*nibz*ndyn*3*ndyn*3))
  

! allocate(CM_mpi(num_k_mpi2,nibz,ndyn,ndyn,ndyn))
! allocate(Collision_Matrix(num_k_mpi2*ndyn*3,nibz*ndyn*3))


!  values1_mpi=product(shape(p1sm_mpi(1,:,:,:,:)))
!  values2_mpi=product(shape(p2sm_mpi(1,:,:,:,:)))

  CM_temp_mpi =0
! P2sm_temp_mpi =0
  cmkn1=ksub_mpi*ndyn*3
  cmkn2=nibz*ndyn*3
  smpi=0
  do n1_mpi=1,cmkn1       !ksub_size2=ksubset(2)-ksubset(1)+1
     do n2_mpi=1,cmkn2

        smpi = smpi +1
        CM_temp_mpi(smpi) = Collision_Matrixmpi(n1_mpi,n2_mpi)
!       write(*,*)rank_mpi,cmkn1,cmkn2,n1_mpi, n2_mpi, Collision_Matrixmpi(n1_mpi,n2_mpi)         
     enddo
  enddo

  deallocate(Collision_Matrixmpi,P1,P2)
  cncount1_mpi=size(Collision_Matrixmpi)
  
  allocate(CM_tempmpi(num_k_mpi2*nibz*ndyn*3*ndyn*3))
  
 
! write(*,*)'rank_mpi,ncount1_mpi,ncount2_mpi,tetttst',rank_mpi,ncount1_mpi,ncount2_mpi,tetttst 

 
! start_timempi = MPI_Wtime()

  call MPI_Allgatherv(CM_temp_mpi,cncount1_mpi,MPI_DOUBLE_PRECISION,CM_tempmpi,caksub_mpi,cdispls_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-gatherv-CM: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately


  deallocate(CM_temp_mpi,cdispls_mpi,caksub_mpi)
  allocate(Collision_Matrix(num_k_mpi2*ndyn*3,nibz*ndyn*3))

  deallocate(CM_temp_mpi,cdispls_mpi,caksub_mpi)
  allocate(Collision_Matrix(num_k_mpi2*ndyn*3,nibz*ndyn*3))

! instead of using reshape function.
  cmkn1t=num_k_mpi2*ndyn*3
  smpi=0
  do n1_mpi=1,cmkn1t       !ksub_size2=ksubset(2)-ksubset(1)+1
     do n2_mpi=1,cmkn2

        smpi=smpi + !
        Collision_Matrix(n1_mpi,n2_mpi) = CM_tempmpi(smpi)

    enddo
  enddo

  deallocate(CM_tempmpi)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' collision-symmetry matrix(IBZ*IBZ) done TIME IS',cputim

  allocate(RHS(nibz*ndyn*3))

! start_timempi = MPI_Wtime()
  call cal_RHS_IBZ(ndyn,nkc,tempk,veloc,eigenval)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' RHS in IBZ done,                TIME IS',cputim


  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')'satrt cg,                TIME IS',cputim
! FFF=RESHAPE(FFF,(/nibz,ndyn,3/),order=[1,2,3])
! DirectSolution=RESHAPE(FFF,(/(nibz*ndyn*3)/))


  allocate(DirectSolution(nibz*ndyn*3))

  s=0
  do is=1,nibz
     do js=1,ndyn
        do ks=1,3
           s=s+1
           DirectSolution(s)=FFF(is,js,ks)
        enddo
     enddo
  enddo


  call cg_quadratic(nibz*ndyn*3,Collision_Matrix,RHS,DirectSolution,fm,100,2)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' cg done,                TIME IS',cputim

! FFF=RESHAPE(DirectSolution,(/nibz,ndyn,3/),order=[3,2,1])
  call mode_kappa_d(ndyn,tempk,veloc,eigenval,DirectSolution)  !Thermal conductivity and tau_eff
! call mode_kappa(ndyn,tempk,veloc,eigenval)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' thermal conductivty calculation, TIME IS',cputim

  call write_kappa(ndyn,tempk,veloc,eigenval,kapd,tauef)

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-RHS and kappa: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately


  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' direct approach using cg method done,                 TIME IS',cputim
  write(ulog,*) 'cg done'

  close(7450)
  close(8450)

call MPI_Finalize(ierr_mpi)
!---------------------------------------end of MPI part------------------------------------- 


!-------------------Self-energy calculation-MPI-----------------------------------
  elseif ( job .eq. 2 ) then
!---------------------------------------
  write(*,*)'self-energy calculation start'
  write(ulog,*) 'self-energy calculation start'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' self-energy calculation,               TIME IS',cputim

  lwse=531
  lwse2=553
  ll = 532

  call make_kp_bs
  allocate(evlbs(ndyn,nkp_bs),evcbs(ndyn,ndyn,nkp_bs),vlcbs(3,ndyn,nkp_bs),omega0(ndyn),tse(ndyn))
  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,evlbs,ndyn,evcbs,ll,vlcbs)
  
  open(lwse,file='lw-selfenergy.dat')
  open(lwse2,file='lw-selfenergy2.dat')
! do itmp=1,ntemp

  tempk=tmin+(tmax-tmin)*(itmp-1)**2/(ntemp-1.0000001d0)**2   ! this temp is in Kelvin
  temp = tempk*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1

! Initialize MPI, get the local number of k_point
  call MPI_INIT(ierr_mpi)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs_mpi,ierr_mpi)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank_mpi,ierr_mpi)


  if (rank_mpi==0) then

      num_k_mpi = nkp_bs

  endif


! At this point, only root knows these values.  Send them out.
  call MPI_Bcast(num_k_mpi,1,MPI_INTEGER,root_mpi,MPI_COMM_WORLD,ierr_mpi)

  local_num_k_mpi = num_k_mpi/nprocs_mpi
  write(*,*) 'num_k_mpi, local_num_k_mpi, nprocs_mpi',num_k_mpi, local_num_k_mpi, nprocs_mpi

! ksubset(1) and ksubset(2) are needed for the loops in the follwoing subroutines (calculate_w3_ibz_split_sq,calculate_FGR_mpi)
  ksubset(1) = rank_mpi * local_num_k_mpi + 1
  ksubset(2) = rank_mpi * local_num_k_mpi + local_num_k_mpi

  if (mod(num_k_mpi,nprocs_mpi).ne.0) then
     if  (rank_mpi.eq.(nprocs_mpi-1)) then
         ksubset(2) = rank_mpi * local_num_k_mpi + local_num_k_mpi + mod(num_k_mpi,nprocs_mpi)
     endif
  endif

  write(*,*)'rank_mpi,ksubset(1),ksubset(2)',rank_mpi,ksubset(1),ksubset(2)

  allocate(isempi((ksubset(2)-ksubset(1)+1),ndyn))
  allocate(tsempi((ksubset(2)-ksubset(1)+1),ndyn))

  sbs=0

! call selfenergyw3(ksubset,ndyn,temp,tempk,nkp_bs,kp_bs,dk_bs)

  do nkbss=ksubset(1),ksubset(2)
     kbss=kp_bs(:,nkbss)             
     sbs = sbs + 1
     do la=1,ndyn
        kbss=kp_bs(:,nkbss)
        omega0(la)=sqrt(abs(evlbs(la,nkbss))) * cnst
        call function_self_w3(kbss,la,omega0(la),temp,nselfs,uselfs)
        isempi(sbs,la)=2*aimag(nselfs+uselfs)
        tsempi(sbs,la)=1d10/c_light/isempi(sbs,la)
!       write(lwse,*)tempk,la,nkbss,kbss,dk_bs(nkbss),omega0(la),aimag(nselfs+uselfs),ise,tse(la),1/tse(la)
!       write(*,*) itmp,ntemp,temp, la, nkbss, nselfs, uselfs
     enddo
!    write(lwse2,*)tempk,dk_bs(nkbss),(omega0(la),la=1,ndyn),((1/tse(la)),la=1,ndyn)
!    write(*,*)tempk,dk_bs(nkbss),(omega0(la),la=1,ndyn),((1/tse(la)),la=1,ndyn)
  enddo

! enddo

  ksub_mpi = (ksubset(2) - ksubset(1)) + 1
! Allocate array to gather ksubsets from all processes
  allocate(aksub_mpi(nprocs_mpi))
  allocate(displs_mpi(nprocs_mpi))

  ksubmpiv=ksub_mpi*ndyn

! Gather ksubs from all processes
  call MPI_Allgather(ksubmpiv, 1, MPI_INTEGER, aksub_mpi, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  displs_mpi(1) = 0
  do i = 2, nprocs_mpi
     displs_mpi(i) = displs_mpi(i-1) + aksub_mpi(i-1)
!    write(*,*)'displs_mpi(i-1),displs_mpi(i),aksub_mpi(i-1),aksub_mpi(i),',i, displs_mpi(i-1),displs_mpi(i),aksub_mpi(i-1),aksub_mpi(i)
  end do

! write(*,*) displs_mpi,aksub_mpi

  allocate(ise_temp_mpi(ksub_mpi*ndyn))
  allocate(tse_temp_mpi(ksub_mpi*ndyn))

  if (mod(num_k_mpi,nprocs_mpi).ne.0) then
      num_k_mpi2 = local_num_k_mpi * nprocs_mpi + mod(num_k_mpi,nprocs_mpi)
  else
      num_k_mpi2= num_k_mpi
  endif

! write(*,*)num_k_mpi2,num_k_mpi,local_num_k_mpi, nprocs_mpi, mod(num_k_mpi,nprocs_mpi)

  allocate(ise_tempmpi(num_k_mpi2*ndyn))
  allocate(tse_tempmpi(num_k_mpi2*ndyn))

  allocate(ise_mpi(num_k_mpi2,ndyn))
  allocate(tse_mpi(num_k_mpi2,ndyn))


  ise_temp_mpi =0
  tse_temp_mpi =0

  smpi=0
  do n1_mpi=1,ksub_mpi       !ksub_size2=ksubset(2)-ksubset(1)+1
     do l3_mpi=1,ndyn

        smpi = smpi +1
        ise_temp_mpi(smpi) = isempi(n1_mpi,l3_mpi)
        tse_temp_mpi(smpi) = tsempi(n1_mpi,l3_mpi)

    enddo
  enddo


  ncount1_mpi=size(isempi)
  ncount2_mpi=size(tsempi)

  call MPI_Allgatherv(ise_temp_mpi,ncount1_mpi,MPI_DOUBLE_PRECISION,ise_tempmpi,aksub_mpi,displs_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)
  call MPI_Allgatherv(tse_temp_mpi,ncount2_mpi,MPI_DOUBLE_PRECISION,tse_tempmpi,aksub_mpi,displs_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)


! instead of using reshape function.
  smpi=0
  do n1_mpi=1,num_k_mpi2       !ksub_size2=ksubset(2)-ksubset(1)+1
     kbss=kp_bs(:,n1_mpi)
     do l1_mpi=1,ndyn
        omega0(l1_mpi)=sqrt(abs(evlbs(l1_mpi,n1_mpi))) * cnst
        smpi=smpi + 1
        ise_mpi(n1_mpi,l1_mpi) = ise_tempmpi(smpi)
        tse_mpi(n1_mpi,l1_mpi) = tse_tempmpi(smpi)
!       write(lwse,*)tempk,la,nkbss,kbss,dk_bs(nkbss),omega0(la),aimag(nselfs+uselfs),ise,tse(la),1/tse(la)
write(lwse,*)tempk,l1_mpi,n1_mpi,kbss,dk_bs(n1_mpi),omega0(l1_mpi),ise_mpi(n1_mpi,l1_mpi)/2,ise_mpi(n1_mpi,l1_mpi),tse_mpi(n1_mpi,l1_mpi),1/tse_mpi(n1_mpi,l1_mpi)

     enddo
! write(lwse2,*)tempk,dk_bs(nkbss),(omega0(la),la=1,ndyn),((1/tse(la)),la=1,ndyn) 
  write(lwse2,*)tempk,dk_bs(n1_mpi),(omega0(l1_mpi),l1_mpi=1,ndyn),((1/tse_mpi(n1_mpi,l1_mpi)),l1_mpi=1,ndyn)
  
  enddo



  close(lwse)
  close(lwse2)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Self-energy calculation done,     TIME IS',cputim
  write(ulog,*) ' Self-energy calculation finished here '
  write(*,*)' Self-energy calculation done'
call MPI_Finalize(ierr_mpi)
!----------------------------------------------------------------------------------------------



!-------------------------------------------------------
if ( job .eq. 3 ) then   ! MPI for very-large kpoints (with job=4)
!-------------------------------------------------------
  write(ulog,*) 'entering cal v33 in IBZ split...'
  write(*,*) 'entering cal v33 in IBZ split...'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')'V33 calculation start,                TIME IS', cputim
       
 !Initialize MPI, get the local number of k_point
  call MPI_INIT(ierr_mpi)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs_mpi,ierr_mpi)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank_mpi,ierr_mpi)

! Create a unique filename for each rank
! write(filenamempi, '("timing_rank_", I0.2, ".txt")') rank_mpi

! Open the file for writing
! unit_numpi = 10 + rank_mpi
! open(unit=unit_numpi, file=filenamempi, status='replace', action='write')

  if (rank_mpi==0) then

      num_k_mpi=nibz

  endif


 !At this point, only root knows these values.  Send them out.
  call MPI_Bcast(num_k_mpi,1,MPI_INTEGER,root_mpi,MPI_COMM_WORLD,ierr_mpi)

      local_num_k_mpi = num_k_mpi/nprocs_mpi
      write(*,*) 'num_k_mpi, local_num_k_mpi, nprocs_mpi',num_k_mpi, local_num_k_mpi, nprocs_mpi

 !ksubset(1) and ksubset(2) are needed for the loops in the follwoing subroutines (calculate_w3_ibz_split_sq,calculate_FGR_mpi)
  ksubset(1) = rank_mpi * local_num_k_mpi + 1
  ksubset(2) = rank_mpi * local_num_k_mpi + local_num_k_mpi
       
  if (mod(num_k_mpi,nprocs_mpi).ne.0) then
     if  (rank_mpi.eq.(nprocs_mpi-1)) then
         ksubset(2) = rank_mpi * local_num_k_mpi + local_num_k_mpi + mod(num_k_mpi,nprocs_mpi)
     endif
  endif  
! write(*,*)'rank_mpi,ksubset1,ksubset2',rank_mpi,ksubset(1),ksubset(2)

  nv3_split=nkc*(ksubset(2)-ksubset(1)+1)*ndyn**3        
  call allocate_v33sq(nv3_split)

  allocate(v33s8((ksubset(2)-ksubset(1)+1),nkc,ndyn,ndyn,ndyn))
! call calculate_w3_fbz_split_sq(ksubset,nv3_split,ndyn,nkc,kpc,eigenval,eigenvec)
! start_timempi = MPI_Wtime()
  call calculate_w3_ibz_split_sq (ksubset,nv3_split,ndyn,nkc,kpc,eigenval,eigenvec,nv3o)
  write(*,*)'after v3'

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-v33: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately
 

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After v33 calculation,                TIME IS ',cputim

  ksub_mpi = (ksubset(2) - ksubset(1)) + 1

! write(*,*)'rank_mpi,ksubset1,ksubset2,ksub_mpi',rank_mpi,ksubset(1),ksubset(2), ksub_mpi
! writing v33 for each rank
! sts_mpi=1
! tt = 1333

! write(filename_mpi,fmt="(a,i3.3,a)") "v33.",rank_mpi,".dat"
! open(tt,file=filename_mpi,status='replace')

! do n1_mpi=1,ksub_mpi
!    do n2_mpi=1,nkc
!       do l1_mpi =1,ndyn
!          do l2_mpi=1,ndyn
!             do l3_mpi=1,ndyn

!                 write(tt,*) rank_mpi, n1_mpi, n2_mpi, l1_mpi, l2_mpi, l3_mpi, v33s8(n1_mpi,n2_mpi,l1_mpi,l2_mpi,l3_mpi)
!              enddo
!           enddo
!        enddo
!     enddo
! enddo

! close(tt)

  allocate(P1smpi((ksubset(2)-ksubset(1)+1),nkc,ndyn,ndyn,ndyn))
  allocate(P2smpi((ksubset(2)-ksubset(1)+1),nkc,ndyn,ndyn,ndyn))
        
  write(ulog,*) 'entering cal v33*delta...'
  write(*,*) 'entering cal v33*delta...'
 
! start_timempi = MPI_Wtime()
! changed to IBZ  call calculate_v3sq_delta (ksubset,ndyn,nkc,kpc)  
  call calculate_FGR_mpi2(ksubset,ndyn,nv3o)!,v33s8)
  write(*,*)'after FGR'
! call deallocate_iter

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-v33delta: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After v33_delta calculation,            TIME IS ',cputim

  deallocate(v33s8)
! call deallocate_v33sq

! writng P1sm and P2sm for each rank
! sts_mpi=1
! tt = 1333
! write(filename_mpi,fmt="(a,i3.3,a)") "P1smpi.",rank_mpi,".dat"
! write(filename_mpii,fmt="(a,i3.3,a)") "P2smpi.",rank_mpi,".dat"
! open(tt,file=filename_mpi,status='replace')
! open(tt+1,file=filename_mpii,status='replace')

! do n1_mpi=1,ksub_mpi
!    do n2_mpi=1,nkc
!       do l1_mpi =1,ndyn
!          do l2_mpi=1,ndyn
!             do l3_mpi=1,ndyn
!                 write(tt,*) rank_mpi, n1_mpi, n2_mpi, l1_mpi, l2_mpi, l3_mpi,P1smpi(n1_mpi,n2_mpi,l1_mpi,l2_mpi,l3_mpi)
!                 write(tt+1,*) rank_mpi, n1_mpi, n2_mpi, l1_mpi, l2_mpi, l3_mpi,P2smpi(n1_mpi,n2_mpi,l1_mpi,l2_mpi,l3_mpi)
!              enddo
!           enddo
!        enddo
!     enddo
! enddo

! close(tt)

!----------------------------------------------------------
  if (iso.eq.1) then     !Piso calculation in IBZ split 
!-----------------------------------------------------------
  write(ulog,*) 'entering Piso calculation in IBZ split'
  write(*,*) 'entering Piso calculation in IBZ split'
  write(ulog,*) 'job', job

  allocate (dist(nkc,ndyn),frequency(nkc,ndyn))
  allocate(Pisos((ksubset(2)-ksubset(1)+1),ndyn))
  tempk = tmin     ! use tmin in params.phon as temperature. tmax is not being used.
  call distribution(nkc,ndyn,tempk)  ! calculate BE distribution function for each state

  write(*,*)'natoms0', natoms0

  allocate(giso(natoms0))
  call cal_giso(giso)
  call P_iso_split(ksubset,giso,ndyn,eigenvec)
  deallocate(dist,frequency)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After Piso calculation,            TIME IS ',cputim
  write(ulog,*) 'Piso calculation done'
  write(ulog,*) 'Piso calculation done'
  endif
!------------------------------------------------------------------------------------------

!------------- Exact solution of BTE--------------------------------------------------------
  write(ulog,*) 'entering to the direct and iterative BTE solver'
  write(*,*) 'entering to the direct and iterative BTE solver'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' entering to the direct and iterative BTE solver,      TIME IS',cputim

! call allocate_iter (nibz,nkc,ndyn)    ! allocate memory
  call allocate_FGR (ksub_mpi,nkc,ndyn,ndyn,ndyn)
  write(ulog,*) 'allocate_iter done'

  tauf=3450
  kapRTA=6450
  open(tauf,file='tau_RTA.dat',status='unknown')
  open(kapRTA,file='kappa_RTA.dat',status='unknown')

  tempk = tmin     ! use tmin in params.phon as temperature. tmax is not being used.

  call distribution(nkc,ndyn,tempk)  ! calculate BE distribution function for each state
  write(ulog,*) 'distribution done'

  allocate(Qvaluempi((ksubset(2)-ksubset(1)+1),ndyn))
  allocate(Qvalue_Nmpi((ksubset(2)-ksubset(1)+1),ndyn))
  allocate(Qvalue_Umpi((ksubset(2)-ksubset(1)+1),ndyn))

!--------------------------------------------
  if (iso.eq.1) then   !reaind Piso if iso=true
!---------------------------------------------  
  call read_Piso(ndyn,Piso) 
  write(ulog,*) 'reading Piso done'
  call calculate_tauRTA_iso(ndyn,tempk)
  else
! start_timempi = MPI_Wtime()
  call calculate_tauRTA_mpi2(ksubset,ndyn,tempk)!,rank_mpi)        
! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-Q: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately
  endif

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Q and tau_RTA in IBZ done,    TIME IS',cputim
  write(ulog,*) 'Q and tau_RTA in IBZ done'


  deallocate(P1smpi,P2smpi)!,Qvalue_N,Qvalue_U,Qvalue_Nmpi,Qvalue_Umpi)
  deallocate(Piso)
  call deallocate_v33sq

  allocate(Qvalue_tempi((ksubset(2)-ksubset(1)+1)*ndyn))
  allocate(Qvalue_Ntempi((ksubset(2)-ksubset(1)+1)*ndyn))
  allocate(Qvalue_Utempi((ksubset(2)-ksubset(1)+1)*ndyn))

  Qvalue_tempi=0
  Qvalue_Ntempi=0
  Qvalue_Ntempi=0
  
  smpi=0 
  do n1_mpi=1,ksub_mpi
     do l1_mpi=1,ndyn
      
        smpi=smpi+1
        Qvalue_tempi(smpi) = Qvaluempi(n1_mpi,l1_mpi)
        Qvalue_Ntempi(smpi) = Qvalue_Nmpi(n1_mpi,l1_mpi)
        Qvalue_Utempi(smpi) = Qvalue_Umpi(n1_mpi,l1_mpi)

    enddo
  enddo

  deallocate(Qvaluempi,Qvalue_Nmpi,Qvalue_Umpi)
! Allocate array to gather ksubs from all processes
  allocate(aksub_mpi(nprocs_mpi))
  allocate(displs_mpi(nprocs_mpi))

  ksubmpiv=ksub_mpi*ndyn

! Gather ksubs from all processes
  call MPI_Allgather(ksubmpiv, 1, MPI_INTEGER, aksub_mpi, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  displs_mpi(1) = 0
  do i = 2, nprocs_mpi
     displs_mpi(i) = displs_mpi(i-1) + aksub_mpi(i-1)
!    write(*,*)'displs_mpi(i-1),displs_mpi(i),aksub_mpi(i-1),aksub_mpi(i),',i, displs_mpi(i-1),displs_mpi(i),aksub_mpi(i-1),aksub_mpi(i)
  end do

  allocate(Qvaluetemp(num_k_mpi*ndyn))
  allocate(QvalueNtemp(num_k_mpi*ndyn))
  allocate(QvalueUtemp(num_k_mpi*ndyn))
  ncount1_mpi=size(Qvaluempi)

! start_timempi = MPI_Wtime()

  call MPI_Allgatherv(Qvalue_tempi,ncount1_mpi,MPI_DOUBLE_PRECISION,Qvaluetemp,aksub_mpi,displs_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)
  call MPI_Allgatherv(Qvalue_Ntempi,ncount1_mpi,MPI_DOUBLE_PRECISION,QvalueNtemp,aksub_mpi,displs_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)
  call MPI_Allgatherv(Qvalue_Utempi,ncount1_mpi,MPI_DOUBLE_PRECISION,QvalueUtemp,aksub_mpi,displs_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-gatherv_RTA: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately


  smpi=0
  do n1_mpi=1,num_k_mpi
     do l1_mpi=1,ndyn

        smpi=smpi+1
        Qvalue(n1_mpi,l1_mpi) = Qvaluetemp(smpi)
        Qvalue_N(n1_mpi,l1_mpi) = QvalueNtemp(smpi)
        Qvalue_U(n1_mpi,l1_mpi) = QvalueUtemp(smpi)
!       write(*,*)n1_mpi,l1_mpi,Qvalue(n1_mpi,l1_mpi),smpi,Qvaluetemp(smpi)
     enddo
  enddo


  deallocate(Qvalue_tempi,Qvaluetemp,aksub_mpi,displs_mpi)
  deallocate(Qvalue_Ntempi,QvalueNtemp,Qvalue_Utempi,QvalueUtemp)

  call Qandtau(ndyn,tempk)

  allocate(FFF(nibz,ndyn,3))
! start_timempi = MPI_Wtime()
! calculates F_RTA(kibz) 
  call calculate_RTA_distributionfunction(nibz,nkc,kibz,ndyn,tempk,veloc,eigenval)
  write(ulog,*) 'FRTA done'

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' FRTA in IBZ done,    TIME IS',cputim

  call mode_kappa(ndyn,tempk,veloc,eigenval)
  call write_kappa(ndyn,tempk,veloc,eigenval,kapRTA,tauf)

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-FRTA and kappa_RTA: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' kappa_RTA calculation done,   TIME IS',cputim

  close(tauf)
  close(kapRTA)

  deallocate(Qvalue)

  write(ulog,*) 'Start direct approach using cg method'
  write(*,*) 'Start direct approach using cg method'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Start direct approach using cg method,                TIME IS',cputim

! tauef=7450
! kapd=8450
! open(tauef,file='tau_RTAandeff-cg.dat',status='unknown')
! open(kapd,file='kappa_direct-cg.dat',status='unknown')


  ksub_mpi = (ksubset(2) - ksubset(1)) + 1

  allocate(Collision_Matrixmpi(ksub_mpi*ndyn*3,nibz*ndyn*3))

!---------------------------------
  if (iso.eq.1) then !for iso scattering
!---------------------------------
     Piso=0
     call read_Piso(ndyn,Piso)
     write(ulog,*) 'reading Piso done'
     call CollisionMatrix_iso(ndyn)
  else
! start_timempi = MPI_Wtime()
  call CollisionMatrix_mpi(ksubset,ndyn)
! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-CM: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately
  endif

! Allocate array to gather ksubs from all processes
  allocate(caksub_mpi(nprocs_mpi))
  allocate(cdispls_mpi(nprocs_mpi))    

  cksubmpiv=ksub_mpi*nibz*ndyn*3*ndyn*3

! Gather ksubs from all processes
  call MPI_Allgather(cksubmpiv, 1, MPI_INTEGER, caksub_mpi, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  cdispls_mpi(1) = 0
  do i = 2, nprocs_mpi
     cdispls_mpi(i) = cdispls_mpi(i-1) + caksub_mpi(i-1)
!    write(*,*)'cdispls_mpi(i-1),cdispls_mpi(i),caksub_mpi(i-1),caksub_mpi(i),',i, cdispls_mpi(i-1),cdispls_mpi(i),caksub_mpi(i-1),caksub_mpi(i)
  end do

! write(*,*) displs_mpi,aksub_mpi

  allocate(CM_temp_mpi(ksub_mpi*nibz*ndyn*3*ndyn*3))
  

  if (mod(num_k_mpi,nprocs_mpi).ne.0) then
      num_k_mpi2 = local_num_k_mpi * nprocs_mpi + mod(num_k_mpi,nprocs_mpi)
  else 
      num_k_mpi2= num_k_mpi
  endif

! write(*,*)num_k_mpi2,num_k_mpi,local_num_k_mpi, nprocs_mpi, mod(num_k_mpi,nprocs_mpi)

! allocate(CM_tempmpi(num_k_mpi2*nibz*ndyn*3*ndyn*3))
  

! allocate(CM_mpi(num_k_mpi2,nibz,ndyn,ndyn,ndyn))
! allocate(Collision_Matrix(num_k_mpi2*ndyn*3,nibz*ndyn*3))


! values1_mpi=product(shape(p1sm_mpi(1,:,:,:,:)))
! values2_mpi=product(shape(p2sm_mpi(1,:,:,:,:)))

  CM_temp_mpi =0
! P2sm_temp_mpi =0
  cmkn1=ksub_mpi*ndyn*3
  cmkn2=nibz*ndyn*3
  smpi=0
  do n1_mpi=1,cmkn1       !ksub_size2=ksubset(2)-ksubset(1)+1
     do n2_mpi=1,cmkn2

        smpi = smpi +1
        CM_temp_mpi(smpi) = Collision_Matrixmpi(n1_mpi,n2_mpi)
!       write(*,*)rank_mpi,cmkn1,cmkn2,n1_mpi, n2_mpi, Collision_Matrixmpi(n1_mpi,n2_mpi)         
     enddo
  enddo

  deallocate(Collision_Matrixmpi,P1,P2)
  cncount1_mpi=size(Collision_Matrixmpi)
  
  allocate(CM_tempmpi(num_k_mpi2*nibz*ndyn*3*ndyn*3))
  
 
! write(*,*)'rank_mpi,ncount1_mpi,ncount2_mpi,tetttst',rank_mpi,ncount1_mpi,ncount2_mpi,tetttst 

! start_timempi = MPI_Wtime()
  ! call MPI_Allgatherv(P1sm_temp_mpi,ncount1_mpi,MPI_DOUBLE_PRECISION,P1sm_tempmpi,aksub_mpi,displs_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)
   call MPI_Allgatherv(CM_temp_mpi,cncount1_mpi,MPI_DOUBLE_PRECISION,CM_tempmpi,caksub_mpi,cdispls_mpi,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr_mpi)

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-gatherv-CM: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately


  deallocate(CM_temp_mpi,cdispls_mpi,caksub_mpi)
! allocate(Collision_Matrix(num_k_mpi2*ndyn*3,nibz*ndyn*3))
  CM_1D=2109
  open(CM_1D,file='CM-1D.dat')
! instead of using reshape function.
  cmkn1t=num_k_mpi2*ndyn*3
  smpi=0
  do n1_mpi=1,cmkn1t       !ksub_size2=ksubset(2)-ksubset(1)+1
     do n2_mpi=1,cmkn2
    
        smpi=smpi + 1
        write(CM_1D,*) CM_tempmpi(smpi) 
              
    enddo
  enddo

  close(CM_1D)


  deallocate(CM_tempmpi)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' collision-symmetry matrix(IBZ*IBZ) done TIME IS',cputim

  allocate(RHS(nibz*ndyn*3))

! start_timempi = MPI_Wtime()
  call cal_RHS_IBZ(ndyn,nkc,tempk,veloc,eigenval)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' RHS in IBZ done,                TIME IS',cputim


  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')'satrt cg,                TIME IS',cputim
        !FFF=RESHAPE(FFF,(/nibz,ndyn,3/),order=[1,2,3])
       ! DirectSolution=RESHAPE(FFF,(/(nibz*ndyn*3)/))


!allocate(DirectSolution(nibz*ndyn*3))
  DS_1D=2106
  open(DS_1D,file='DS.dat')

  s=0
  do is=1,nibz
     do js=1,ndyn
        do ks=1,3
           s=s+1
      !    DirectSolution(s)=FFF(is,js,ks)
          write(DS_1D,*) FFF(is,js,ks)
        enddo
     enddo
  enddo

  close(DS_1D)

call MPI_Finalize(ierr_mpi)
!---------------------------------------end of MPI part------------------------------------- 


!----if CM is too big, this can be used as the second part to read CM and calculate kappa_direct-------------------
  elseif (job .eq. 4) then !(after running job=3)
!---------------------------------------------------
  write(ulog,*) 'Start direct approach using cg method'
  write(*,*) 'Start direct approach using cg method'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Start direct approach using cg method,                TIME IS',cputim

  call allocate_iter(nibz,nkc,ndyn)
  tauef=7450
  kapd=8450
  open(tauef,file='tau_RTAandeff-cg.dat',status='unknown')
  open(kapd,file='kappa_direct-cg.dat',status='unknown')
  tempk = tmin     ! use tmin in params.phon as temperature. tmax is not being used.

  allocate(Collision_Matrix(nibz*ndyn*3,nibz*ndyn*3),RHS(nibz*ndyn*3))

  call read_CM_1D_mpi(nibz,ndyn,Collision_Matrix)

  write(*,*)'after reading CM'

! start_timempi = MPI_Wtime()
  
  call read_RHS(nibz,ndyn,RHS)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' RHS in IBZ done,                TIME IS',cputim       
  write(*,*)'after reading RHS'

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')'satrt cg,                TIME IS',cputim
  
  allocate(DirectSolution(nibz*ndyn*3))

  call read_DS(nibz,ndyn,DirectSolution)

  write(*,*)'after reading DS'

  call cg_quadratic(nibz*ndyn*3,Collision_Matrix,RHS,DirectSolution,fm,100,2)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' cg done,                TIME IS',cputim

  allocate(FFF(nibz,ndyn,3))

! FFF=RESHAPE(DirectSolution,(/nibz,ndyn,3/),order=[3,2,1])
  call mode_kappa_d(ndyn,tempk,veloc,eigenval,DirectSolution)  !Thermal conductivity and tau_eff
! call mode_kappa(ndyn,tempk,veloc,eigenval)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' thermal conductivty calculation, TIME IS',cputim

  call write_kappa(ndyn,tempk,veloc,eigenval,kapd,tauef)

! end_timempi = MPI_Wtime()
! elapsed_timempi = end_timempi - start_timempi
! write(unit_numpi, '(A, F6.2)') "Time after subroutine-RHS and kappa: ", elapsed_timempi
! flush(unit=unit_numpi)  ! Ensure the data is written to the file immediately


  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' direct approach using cg method done,                 TIME IS',cputim
  write(ulog,*) 'cg done'
  close(7450)
  close(8450)
stop
!------------------------------------------------------------------------------------------



!-------------------Self-energy calculation wihout MPI and split----------------------------
  elseif ( job .eq. 5 ) then
!---------------------------------------
  write(*,*)'self-energy calculation start'
  write(ulog,*) 'self-energy calculation start'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' self-energy calculation,               TIME IS',cputim

  lwse=531
  lwse2=553
  ll=532

  call make_kp_bs
! write(*,*)' dk_bs =', dk_bs(: )
  allocate(evlbs(ndyn,nkp_bs),evcbs(ndyn,ndyn,nkp_bs),vlcbs(3,ndyn,nkp_bs),omega0(ndyn),tse(la))
  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,evlbs,ndyn,evcbs,ll,vlcbs)
  open(lwse,file='lw-selfenergy.dat')
  open(lwse2,file='lw-selfenergy2.dat')
  do itmp=1,ntemp

      tempk=tmin+(tmax-tmin)*(itmp-1)**2/(ntemp-1.0000001d0)**2   ! this temp is in Kelvin
      temp = tempk*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1

      do nkbss=1,nkp_bs
        do la=1,ndyn
!           do nkbss=1,nkp_bs

               kbss=kp_bs(:,nkbss)
               omega0(la)=sqrt(abs(evlbs(la,nkbss))) * cnst
               call function_self_w3(kbss,la,omega0(la),temp,nselfs,uselfs)
               ise=2*aimag(nselfs+uselfs)
               tse(la)=1d10/c_light/ise
               write(lwse,*)tempk,la,nkbss,kbss,dk_bs(nkbss),omega0(la),aimag(nselfs+uselfs),ise,tse(la),1/tse(la)
               write(*,*) itmp,ntemp,temp, la, nkbss, nselfs, uselfs
         enddo
         write(lwse2,*)tempk,dk_bs(nkbss),(omega0(la),la=1,ndyn),((1/tse(la)),la=1,ndyn)
      enddo
  enddo

  close(lwse)
  close(lwse2)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Self-energy calculation done,     TIME IS',cputim
  write(ulog,*) ' Self-energy calculation finished here '
  write(*,*)' Self-energy calculation done'
!------------------------------------------------------------------------------------------


!-------------------Self-energy calculation-split (part1)-----------------------------------
  elseif ( job .eq. 6 ) then
!---------------------------------------
  write(*,*)'self-energy calculation start'
  write(ulog,*) 'self-energy calculation start'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' self-energy calculation,               TIME IS',cputim

  lwse=531
  lwse2=553
  ll = 532

  call make_kp_bs
  allocate(evlbs(ndyn,nkp_bs),evcbs(ndyn,ndyn,nkp_bs),vlcbs(3,ndyn,nkp_bs),omega0(ndyn),tse(ndyn))
  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,evlbs,ndyn,evcbs,ll,vlcbs)
  call read_ksubset(ksubset)
  open(lwse,file='lw-selfenergy.dat')
  open(lwse2,file='lw-selfenergy2.dat')
  do itmp=1,ntemp

      tempk=tmin+(tmax-tmin)*(itmp-1)**2/(ntemp-1.0000001d0)**2   ! this temp is in Kelvin
      temp = tempk*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1

      do nkbss=ksubset(1),ksubset(2)
         do la=1,ndyn
             kbss=kp_bs(:,nkbss)
             omega0(la)=sqrt(abs(evlbs(la,nkbss))) * cnst
             call function_self_w3(kbss,la,omega0(la),temp,nselfs,uselfs)
             ise=2*aimag(nselfs+uselfs)
             tse(la)=1d10/c_light/ise
             write(lwse,*)tempk,la,nkbss,kbss,dk_bs(nkbss),omega0(la),aimag(nselfs+uselfs),ise,tse(la),1/tse(la)
!            write(*,*) itmp,ntemp,temp, la, nkbss, nselfs, uselfs
          enddo
          write(lwse2,*)tempk,dk_bs(nkbss),(omega0(la),la=1,ndyn),((1/tse(la)),la=1,ndyn)
          write(*,*)tempk,dk_bs(nkbss),(omega0(la),la=1,ndyn),((1/tse(la)),la=1,ndyn)
     enddo
  enddo

  close(lwse)
  close(lwse2)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Self-energy calculation done,     TIME IS',cputim
  write(ulog,*) ' Self-energy calculation finished here '
  write(*,*)' Self-energy calculation done'
!----------------------------------------------------------------------------------------------
!-------------------Self-energy calculation-split (part2)--------------------------------------
  elseif ( job .eq. 7 ) then
!---------------------------------------
  write(ulog,*) 'self-energy reaing part start split...'
  write(*,*) 'self-energy reaing part start split...'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')'self-energy reaing part start split...,                TIME IS', cputim


  call make_kp_bs
  call read_lw(ndyn,nkp_bs)

  call read_lw2(ndyn,nkp_bs)
 
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')'self-energy reaing part done, TIME IS', cputim
  write(ulog,*) 'Self-energy calculation reaidng part done here '
  write(*,*)'Self-energy reading part done'
!------------------------------------------------------------------------------------------------



!----------The anharmonic V3 matrix calculation in IBZ-------------------------------------
  nv3 = nkc*nibz*ndyn**3
!-------------------------------------------------------
  elseif ( job .eq. 8 ) then   ! v33 in IBZ split...
!-------------------------------------------------------
  write(ulog,*) 'entering cal v33 in IBZ split...'
  write(*,*) 'entering cal v33 in IBZ split...'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')'V33 calculation start,                TIME IS', cputim
       
  call read_ksubset(ksubset)  ! read index of initial and final kpoints set in the shellscript
  write(*,*)'after read_ksub'
  nv3_split=nkc*(ksubset(2)-ksubset(1)+1)*ndyn**3        
  call allocate_v33sq(nv3_split)

  allocate(v33s8((ksubset(2)-ksubset(1)+1),nkc,ndyn,ndyn,ndyn))
!  call calculate_w3_fbz_split_sq(ksubset,nv3_split,ndyn,nkc,kpc,eigenval,eigenvec)
   call calculate_w3_ibz_split_sq (ksubset,nv3_split,ndyn,nkc,kpc,eigenval,eigenvec,nv3o)
   write(*,*)'after v3'
   call cpu_time(cputim)  !date_and_time(time=tim)
   write(utimes,'(a,f12.4)')' After v33 calculation,                TIME IS ',cputim

      
   write(ulog,*) 'entering cal v33*delta...'
   write(*,*) 'entering cal v33*delta...'
! changed to IBZ        call calculate_v3sq_delta (ksubset,ndyn,nkc,kpc)  
  call calculate_FGR(ksubset,ndyn,nv3o)!,v33s8)
  write(*,*)'after FGR'
  call deallocate_iter

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After v33_delta calculation,            TIME IS ',cputim
! write(*,*)'iso',iso
!----------------------------------------------------------
  if (iso.eq.1) then     !Piso calculation in IBZ split 
!-----------------------------------------------------------
   write(ulog,*) 'entering Piso calculation in IBZ split'
   write(*,*) 'entering Piso calculation in IBZ split'
   write(ulog,*) 'job', job

   allocate (dist(nkc,ndyn),frequency(nkc,ndyn))
   allocate(Pisos((ksubset(2)-ksubset(1)+1),ndyn))
   tempk = tmin     ! use tmin in params.phon as temperature. tmax is not being used.
   call distribution(nkc,ndyn,tempk)  ! calculate BE distribution function for each state
   write(ulog,*) 'distribution done'
   write(*,*)'after distribution'

   allocate(giso(natoms0))

   call cal_giso(giso)

   call P_iso_split(ksubset,giso,ndyn,eigenvec)
   write(*,*)'after P_iso_split'
   deallocate(dist,frequency)
   call cpu_time(cputim)  !date_and_time(time=tim)
   write(utimes,'(a,f12.4)')' After Piso calculation,            TIME IS ',cputim
   write(ulog,*) 'Piso calculation done'
   write(ulog,*) 'Piso calculation done'
  endif
 
!------------------------------------------------------------------------------------------

 
!------------- Exact solution of BTE--------------------------------------------------------
  elseif (job .ge. 9) then
!------------------------------------------
  write(ulog,*) 'entering to the direct and iterative BTE solver'
  write(*,*) 'entering to the direct and iterative BTE solver'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' entering to the direct and iterative BTE solver,      TIME IS',cputim

  call allocate_iter (nibz,nkc,ndyn)    ! allocate memory
  call allocate_FGR (nibz,nkc,ndyn,ndyn,ndyn)
  write(ulog,*) 'allocate_iter done'

  tauf=3450
  kapRTA=6450
  open(tauf,file='tau_RTA.dat',status='unknown')
  open(kapRTA,file='kappa_RTA.dat',status='unknown')

  tempk = tmin     ! use tmin in params.phon as temperature. tmax is not being used.

  call distribution(nkc,ndyn,tempk)  ! calculate BE distribution function for each state
  write(ulog,*) 'distribution done'
!--------------------------------------------
  if (iso.eq.1) then   !reaind Piso if iso=true
!---------------------------------------------  
      call read_Piso(ndyn,Piso) 
      write(ulog,*) 'reading Piso done'
      call calculate_tauRTA_iso(ndyn,tempk)
  else
      call calculate_tauRTA(ndyn,tempk)        
  endif
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Q and tau_RTA in IBZ done,    TIME IS',cputim
  write(ulog,*) 'Q and tau_RTA in IBZ done'

  allocate(FFF(nibz,ndyn,3))
! call calculate_RTA_distributionfunction(nibz,kibz,ndyn,tempk,veloc,eigenval)  ! calculates F_RTA(kibz) 
  call calculate_RTA_distributionfunction(nibz,nkc,kibz,ndyn,tempk,veloc,eigenval)
  write(ulog,*) 'FRTA done'

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' FRTA in IBZ done,    TIME IS',cputim

  call mode_kappa(ndyn,tempk,veloc,eigenval)
  call write_kappa(ndyn,tempk,veloc,eigenval,kapRTA,tauf)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' kappa_RTA calculation done,   TIME IS',cputim

  close(tauf)
  close(kapRTA)
 
!-------------------------------------------------------
  if(job .eq. 10) then  ! Direct noniterative method-svd 
!-------------------------------------------------------
  write(ulog,*) 'Start direct approach'
  write(*,*) 'Start direct approach'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Start direct approach,                TIME IS ',cputim

  tauef=7450
  kapd=8450
  open(tauef,file='tau_RTAandeff.dat',status='unknown')
  open(kapd,file='kappa_direct.dat',status='unknown')

  allocate(DirectSolution(nibz*ndyn*3),sig(nibz*ndyn*3))
  allocate(Collision_Matrix(nibz*ndyn*3,nibz*ndyn*3),RHS(nibz*ndyn*3))
!-------------------------------------
  if (iso.eq.1) then !for iso scattering 
!-------------------------------------
     call read_Piso(ndyn,Piso)
     write(ulog,*) 'reading Piso done'

     call CollisionMatrix_iso(ndyn)
  else
     call CollisionMatrix(ndyn)
  endif
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' collision-symmetry matrix(IBZ*IBZ) done                TIME IS',cputim

  call cal_RHS_IBZ(ndyn,nkc,tempk,veloc,eigenval)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' RHS in IBZ done,                TIME IS',cputim

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')'satrt SVD,                TIME IS',cputim
  svdcut = 1d-9
  call svd_solver(nibz*ndyn*3,nibz*ndyn*3,Collision_Matrix,RHS,DirectSolution,sig,svdcut,errorr,ermax,sigma,'svd_output_IBZ2.dat')
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' SVD done,                TIME IS',cputim

  FFF=RESHAPE(DirectSolution,(/nibz,ndyn,3/),order=[3,2,1])
! call mode_kappa_d(ndyn,tempk,veloc,eigenval,DirectSolution)  !Thermal conductivity and tau_eff
  call mode_kappa(ndyn,tempk,veloc,eigenval)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' thermal conductivty calculation, TIME IS',cputim
  call write_kappa(ndyn,tempk,veloc,eigenval,kapd,tauef)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' direct approach done,                 TIME IS ',cputim
  write(ulog,*) 'SVD done'

  close(7450)
  close(8450)



!---------------------------------------------------------------
  elseif(job .eq. 11) then   ! Direct noniterative method-split-svd
!---------------------------------------------------------------
  write(ulog,*) 'direct approach-split start'
  write(*,*) 'direct approach-split start'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Start direct approach-split,                TIME IS',cputim

  call read_ksubset(ksubset)
  allocate(Collision_Matrix((ksubset(2)-ksubset(1)+1)*ndyn*3,nibz*ndyn*3))
  call CollisionMatrix_split(ksubset,ndyn)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' collision-symmetry matrix(IBZ*IBZ)-split done, TIME IS',cputim

!----------------------------------------------------------------------------
  elseif (job .eq. 12) then  !reading collision Matrix_split calculated by job=11  
!----------------------------------------------------------------------------
  write(ulog,*) 'direct approach-spit reading colliaion Matrix'
  write(*,*) 'direct approach-spiit reading collision matrix'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Start reading CollisionMatrix for direct approach-split,          TIME IS',cputim

  tauef=7450
  kapd=8450
  open(tauef,file='tau_RTAandeff.dat',status='unknown')
  open(kapd,file='kappa_direct.dat',status='unknown')

  allocate(DirectSolution(nibz*ndyn*3),sig(nibz*ndyn*3))
  allocate(Collision_Matrix(nibz*ndyn*3,nibz*ndyn*3),RHS(nibz*ndyn*3))

  call read_Collision_Matrix(ndyn,Collision_Matrix)
! call read_RHS(nibz,ndyn,RHS2)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' reading collision matrix, TIME IS',cputim

  call cal_RHS_IBZ(ndyn,nkc,tempk,veloc,eigenval)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' RHS in IBZ done,                TIME IS',cputim

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' start svd approach,                TIME IS ',cputim

  svdcut = 1d-9
  call svd_solver(nibz*ndyn*3,nibz*ndyn*3,Collision_Matrix,RHS,DirectSolution,sig,svdcut,errorr,ermax,sigma,'svd_output_IBZ2.dat')
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' svd done,                TIME IS ',cputim

  FFF=RESHAPE(DirectSolution,(/nibz,ndyn,3/),order=[3,2,1])
  call mode_kappa(ndyn,tempk,veloc,eigenval)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After kappa calculation,                TIME IS ',cputim
  call write_kappa(ndyn,tempk,veloc,eigenval,kapd,tauef)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' direct approach-split done,                 TIME IS',cputim
  write(ulog,*) 'direct approach done'


  close(7450)
  close(8450)



!-------------------------------------------------------
  elseif(job .eq. 13) then  ! Direct noniterative method-cg
!-------------------------------------------------------
  write(ulog,*) 'Start direct approach using cg method'
  write(*,*) 'Start direct approach using cg method'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' Start direct approach using cg method,                TIME IS',cputim

  tauef=7450
  kapd=8450
  open(tauef,file='tau_RTAandeff-cg.dat',status='unknown')
  open(kapd,file='kappa_direct-cg.dat',status='unknown')

  allocate(DirectSolution(nibz*ndyn*3))
  allocate(Collision_Matrix(nibz*ndyn*3,nibz*ndyn*3),RHS(nibz*ndyn*3))
! call CollisionMatrix(ndyn)
!---------------------------------
  if (iso.eq.1) then !for iso scattering        
!---------------------------------
      Piso=0
      call read_Piso(ndyn,Piso)
      write(ulog,*) 'reading Piso done'
      call CollisionMatrix_iso(ndyn)
  else
      call CollisionMatrix(ndyn)
  endif
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' collision-symmetry matrix(IBZ*IBZ) done TIME IS',cputim


  call cal_RHS_IBZ(ndyn,nkc,tempk,veloc,eigenval)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' RHS in IBZ done,                TIME IS',cputim


  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')'satrt cg,                TIME IS',cputim
! FFF=RESHAPE(FFF,(/nibz,ndyn,3/),order=[1,2,3])
! DirectSolution=RESHAPE(FFF,(/(nibz*ndyn*3)/))

  s=0
  do is=1,nibz
     do js=1,ndyn
        do ks=1,3
           s=s+1
           DirectSolution(s)=FFF(is,js,ks)
        enddo
     enddo
  enddo


  call cg_quadratic(nibz*ndyn*3,Collision_Matrix,RHS,DirectSolution,fm,100,2)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' cg done,                TIME IS',cputim


! FFF=RESHAPE(DirectSolution,(/nibz,ndyn,3/),order=[3,2,1])
  call mode_kappa_d(ndyn,tempk,veloc,eigenval,DirectSolution)  !Thermal conductivity and tau_eff
! call mode_kappa(ndyn,tempk,veloc,eigenval)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' thermal conductivty calculation, TIME IS',cputim

  call write_kappa(ndyn,tempk,veloc,eigenval,kapd,tauef)

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' direct approach using cg method done,                 TIME IS',cputim
  write(ulog,*) 'cg done'

  close(7450)
  close(8450)




!------------------------------------------------------
  elseif (job.eq.14) then  ! iterative solution here
!------------------------------------------------------
  write(ulog,*) 'iterative approach start'
  write(*,*) 'iterative approach start'
  write(ulog,*) 'job', job

  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' start iterative approach,                 TIME IS',cputim

  tauefi=9450
  kappai=4451
  open(tauefi,file='tau_RTAandeff_iter.dat',status='unknown')
  open(kappai,file='kappa_iter.dat',status='unknown')
 
  allocate(kappa_old(ndyn,3,3),kappa_new(ndyn,3,3),e(ndyn,3,3))
  F_RTA=FFF

  do i=1, max_iter     ! max_iter : maximum iteration allowed by params.phon

     write(ulog,*) '================= iteration',i
     call mode_kappa(ndyn,tempk,veloc,eigenval)
     kappa_old=kappa
     call cpu_time(cputim)  !date_and_time(time=tim)
     write(utimes,'(a,f12.4)')' kappa_old calculation done,              TIME IS',cputim

! calculate F value by using symmetry operations and using P, Qvalue in IBZ. 
    call update(ndyn)
    !------------------------------
    if (iso.eq.1) then !for iso scattering
    !-----------------------------        
       call update_iso(ndyn)
    else
       call update(ndyn)
    endif
    write(ulog,*) 'update F done'
    FFF(:,:,:) = FFF(:,:,:) + update_mix*(F2(:,:,:)-FFF(:,:,:))  ! update F with mixing
    write(ulog,*) 'update F done'
    call cpu_time(cputim)  !date_and_time(time=tim)
    write(utimes,'(a,f12.4)')' update F done,             TIME IS',cputim

    call mode_kappa(ndyn,tempk,veloc,eigenval)
    kappa_new=kappa
    write(ulog,*) 'mode_kappa done for kappa_new iteration',i
    write(kappai,*)'================= iteration',i
    call write_kappa(ndyn,tempk,veloc,eigenval,kappai,tauefi)
    call cpu_time(cputim)  !date_and_time(time=tim)
    write(utimes,'(a,f12.4)')' kappa_new calculation and writing it done ,    TIME IS',cputim

    write(*,*)'================= iteration',i
    normk=0
    e=0
    e=(kappa_new-kappa_old)**2/(ndyn*3*3)
    normk=sqrt( sum(e(:,:,:)))
    write(*,*)'normk2', normk
    if (normk .lt. 0.0001 ) then
        exit
    else
        cycle
    endif
    write(ulog,*) 'check_conv done...'
    call cpu_time(cputim)  !date_and_time(time=tim)
    write(utimes,'(a,f12.4)')' check_conv done ,TIME IS',cputim

  enddo

  call deallocate_FGR
  close(tauefi)
  close(kappai)
  call cpu_time(cputim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' iterative approach done,     TIME IS ',cputim


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
 close(uibz)
 close(bzf)
 close(fbz)
 close(ibs)
 close(ucors)

endif
endif
endif

 end program kappa_full



!===============================================================
 subroutine get_kindex(q,nkc,kpc,ik)
  implicit none
  integer ik,nkc,i,j
  real(8) kpc(3,nkc),q(3)

  ik=0
  mainlp: do i=1,nkc
     do j=1,3
        if (abs(q(j)-kpc(j,i)) .ge. 1d-4) exit
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

