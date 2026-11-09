!=========================================================
 program kappa_full
!! THERMACOND: phonon properties and lattice thermal conductivity from FOCEX force constants.
!! Inputs (same files as FOCEX): structure.params, latdyn.params, dielectric.params,
!! kpbs.params, lat_fc.dat, fc2.dat, fc3.dat.  latdyn.params line 8 selects the job:
!!   0  bands, Gruneisen along kpbs.params, DOS on the mesh (always done)
!!   1  0 + three-phonon linewidths along the band path, for every temperature
!!   2  0 + V3*delta in the irreducible zone, then for every temperature the RTA and the
!!        direct solution of the linearized BTE (line 9: 0 conjugate gradient, 1 SVD)
!!   3  0 + V3*delta, then for every temperature the RTA and the iterative solution
!! V3*delta does not depend on temperature and is computed once.  latdyn.params line 10:
!! 0 = one file per irreducible k-point in V3DELTA/ (reused by a later run on the same mesh),
!! 1 = kept in memory.
!! Serial and parallel runs are the same program:  tcond-2.2.x   or   mpirun -np N tcond-2.2.x
!! 2.2 (2026-11-09): mode-resolved and cumulative kappa (kappa_modes_T*, kappa_cumul_mfp_T*, kappa_cumul_freq_T*,
!! kappa_summary, kappa_wigner_modes .dat), tau_eff fixed, 5-figure output, times in h:m:s.
!! 2.1 (2026-11-09): |V3|^2 only where the tetrahedron delta weight is nonzero (v3_delta), and V3DELTA
!! files with the nonzero elements only (the dense files of 2.0 are still read).
!! Rank 0 works in the run folder; the other ranks work in .mpi_rank_NNN/ (removed at the end),
!! so only rank 0's output files are kept.
 use ios, only : ulog,uposcar,uparams,utraj,ufco,umap,umatrx,utimes,ufc1,ufc2,ufc3,ufc4,ufc,ufit1,ufit2,ufit3,ufit4,uborn,write_out
 use tc_ios
 use kpoints
 use om_dos
 use atoms_force_constants
 use phi3
 use phi3_sy
 use eigen
 use params, only : tolerance,margin,nconfigs,classical,fdfiles,threemtrx,lamin,lamax,nshells,include_fc,tmin,tmax,verbose
 use tc_params
 use constants
 use born, only : epsil,born_flag
 use exactBTE2
 use mod_ksubset2
 use tetrahedron
 use mpi
 use iso_c_binding
!$ use omp_lib
 implicit none
 integer nthreads
 integer(8) wclock0,wclock1,wrate
 character(64) thrsrc

 interface
    integer(c_int) function c_chdir(path) bind(C,name='chdir')
      import c_int, c_char
      character(kind=c_char) :: path(*)
    end function c_chdir
    integer(c_int) function c_mkdir(path,mode) bind(C,name='mkdir')
      import c_int, c_char
      character(kind=c_char) :: path(*)
      integer(c_int), value :: mode
    end function c_mkdir
    integer(c_int) function c_symlink(target,linkpath) bind(C,name='symlink')
      import c_int, c_char
      character(kind=c_char) :: target(*), linkpath(*)
    end function c_symlink
    type(c_ptr) function c_getcwd(buf,size) bind(C,name='getcwd')
      import c_ptr, c_char, c_size_t
      character(kind=c_char) :: buf(*)
      integer(c_size_t), value :: size
    end function c_getcwd
 end interface

 integer i,la,itmp,k,kl,ierr,nloc,s,is,js,ks,nit
 integer kglob(2),ksub(2)
 integer rank,nprocs
 logical root
 real(8) tempk,temp,normk,fm
 real cputim
 character now*10,today*8,zone*5,rankdir*20,runpath*400,v3dir*420
 real(8), allocatable :: dkc(:),evc(:)
 real(8), allocatable :: Qloc(:,:),QNloc(:,:),QUloc(:,:),CMloc(:,:),buf(:),gbuf(:)
 real(8), allocatable :: P1all(:,:,:,:,:),P2all(:,:,:,:,:),F2loc(:,:,:)
 integer, allocatable :: counts(:),displs(:)
 integer ukvt,kapRTA,tauf,kapd,tauef,kappai,tauefi,udist,ukw
 real(8) kRTA(3,3),kdir(3,3),resid,kWP(3,3),kWC(3,3),kblk(3,3)
 ! per IBZ point and mode, star summed (2.1, the mode-resolved and cumulative kappa): RTA (mode_kappa
 ! with F_RTA), Wigner kappa_P and kappa_C (half of each pair s,s' to each mode), RTA rates 1/tau
 real(8), allocatable :: cR_k(:,:,:,:),wP_k(:,:,:,:),wC_k(:,:,:,:),tauinvR(:,:)
 integer usum,uwm
! modes closer than degtol (1/cm) are treated as degenerate (|V3|^2 averaging, Wigner kappa_P).
! Symmetry-degenerate modes of Ge come out split by up to 1.5e-5 1/cm (4e-7 relative), while
! the closest non-degenerate ones are 1e-2 1/cm apart.
 real(8), parameter :: degtol=1d-3

 call MPI_INIT(ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
 root = (rank.eq.0)
 call get_cwd(runpath)
 v3dir = trim(runpath)//'/V3DELTA/'
 if (.not.root) call enter_rank_dir()

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(cputim)
 call system_clock(wclock0,wrate)

  open(ulog,file='phonlog.dat',status='unknown')
  open(utimes,file='times.dat' ,status='unknown')
  open(ueiv,file='bs_freqs.dat',status='unknown')
  open(ucors,file='FBZ_freqs.dat',status='unknown')
  open(ugrun,file='bs_grun.dat',status='unknown')
  open(ugrun+1,file='FBZ_grun.dat',status='unknown')

 write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' THERMACOND 2.2 launched at  ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(ulog,*)' MPI processes: ',nprocs
! OpenMP threads (common/omp_threads.f90): with more than one MPI process, 1 unless TCOND_NUM_THREADS is set
 call omp_threads_choose('TCOND_NUM_THREADS',nthreads,thrsrc)
 if (nprocs.gt.1 .and. trim(thrsrc).ne.'TCOND_NUM_THREADS') then
    nthreads=1 ; thrsrc='more than one MPI process'
!$  call omp_set_num_threads(1)
 endif
 write(ulog,'(a,i4,a)')' OpenMP threads per process: ',nthreads,'  ('//trim(thrsrc)//'; set TCOND_NUM_THREADS to choose)'
 if (root) write(*,'(a,i4,a)')' OpenMP threads per process: ',nthreads,'  ('//trim(thrsrc)//'; set TCOND_NUM_THREADS to choose)'
 write(utimes,'(a,i4,a)')' OpenMP threads per process: ',nthreads,'  ('//trim(thrsrc)//')'
 write(utimes,'(a)')' (CPU: time summed over the OpenMP threads; WALL: elapsed time; hours:minutes:seconds)'
 write(utimes,'(a,a)')' At the start of the Program,  CPU ',hms(dble(cputim))

! this is used to convert v33 to cm^_1
  const33 = 1d30*(hbar/uma)**1.5d0/(sqrt(ee*1d20/uma))**1.5d0*ee/h_plank/100/c_light

  call read_params
  call read_input_fit
  call read_lattice  ! reads lat_fc.dat
  call read_born
  call setup_nonanalytic
  nkc = nc(1)*nc(2)*nc(3)
  call read_fc234
  call timing('After readings')
  ndyn = 3*natoms0
  if (iso.eq.1 .and. job.ge.2) then
     if (root) write(*,*)'latdyn.params line 11: isotope scattering is not yet available in this version; set it to 0'
     call MPI_FINALIZE(ierr)
     stop
  endif
  if (root) write(*,'(a,i2,a,i4,a)')' THERMACOND: job ',job,' on ',nprocs,' process(es)'

!---------------------------------------------------------------
! job 0 part, always done: band structure and Gruneisen along kpbs.params, mesh, IBZ, DOS
  call make_kp_bs
  call allocate_eig_bs(ndyn,nkp_bs,ndyn) !3rd index is for eivecs needed for gruneisen
  call get_frequencies(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,ndyn,eigenvec_bs,ueiv,veloc,1)
  open(1008,file='nkpbs.dat')
  write(1008,*)nkp_bs
  close(1008)
  call gruneisen(nkp_bs,kp_bs,dk_bs,ndyn,eigenval_bs,eigenvec_bs,ugrun,grun_bs)
  allocate(kpc(3,nkc),mappos(nkc),wk(nkc))
  call make_kp_reg(nc,g01,g02,g03,shft,kpc,wk)   ! integer fractions of the primitive reciprocal vectors
  call deallocate_eig_bs
  call deallocate_kp_bs
  call timing('After band structure')

  call get_negatives(nkc,kpc,mappos,npos)
  call get_weights(nkc,kpc)            ! irreducible zone: nibz, kibz, mapinv, mapibz
  call allocate_eig(ndyn,nkc,nkc)      ! last one (nkc) is for veloc
  allocate(dkc(nkc))
  dkc=0
  call get_frequencies(nkc,kpc,dkc,ndyn,eigenval,ndyn,eigenvec,ucors,veloc,0)
  call timing('After om(k) calculation within FBZ')
  call check_mdyn(ndyn,nkc,kpc,eigenvec,eigenval)
  if(verbose) call write_eivecs(ndyn,nkc,kpc,eigenvec)
  call subst_eivecs(ndyn,nkc,eigenvec,kpc,npos,mappos)
  call mark_gamma_acoustic

  call set_omdos(ndyn,wmesh)  ! allocates om and dos arrays of size mesh
  allocate(evc(nibz))
  dos = 0
  do la=1,ndyn
     do i=1,nibz
        evc(i)=eigenval(la,mapinv(i))
     enddo
     call calculate_dos_g(nibz,evc,wibz,wmesh,om,dos(la,:))  ! with gaussian broadening
     dos(ndyn+1,:) = dos(ndyn+1,:) + dos(la,:)
  enddo
  call write_dos
  deallocate(evc)
  call allocate_tetra(nc(1),nc(2),nc(3),wmesh,ndyn)
  call make_kp_reg_tet(nc(1),nc(2),nc(3),shft(1),shft(2),shft(3),kpc,wk)
  call check_mesh_labels
  call calc_tet(nc(1),nc(2),nc(3),ndyn,kpc)
  deallocate(dos,dkc)
  deallocate(grun)
  call timing('After DOS (gaussian and tetrahedron)')
  write(ulog,*) ' job 0 part (bands, Gruneisen, DOS) done'

 select case (job)
 case (1)
    call linewidths_on_band_path()
 case (2,3)
    call v3_delta_ibz()
    call kappa_all_temperatures()
 end select

 call timing('End of the program')
 close(ulog)
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 if (root .and. nprocs.gt.1) call execute_command_line('rm -rf .mpi_rank_*')
 call MPI_FINALIZE(ierr)

 contains

!---------------------------------------------------------------
 subroutine check_mesh_labels()
!! every routine finds mesh points with get_k_info (index from the reduced coordinates on
!! g1,g2,g3, modulo G; Gamma-centred labelling for even meshes, as made by make_kp_reg and
!! make_kp_reg_tet). Stop if kpc(:,i), kpc(:,i)+G or the IBZ representatives are not found
!! at their own index, or if -k is not the inverse of k.
 integer ii,nn,a1,a2,a3,ins,nbad
 real(8) qq(3)
 nbad=0
 do ii=1,nkc
    call get_k_info(kpc(:,ii),NC,nn,a1,a2,a3,g1,g2,g3,ins)
    if (nn.ne.ii) nbad=nbad+1
    qq=kpc(:,ii)+g1-2*g2+g3
    call get_k_info(qq,NC,nn,a1,a2,a3,g1,g2,g3,ins)
    if (nn.ne.ii) nbad=nbad+1
    call get_k_info(-kpc(:,ii),NC,nn,a1,a2,a3,g1,g2,g3,ins)
    qq=kpc(:,nn)+kpc(:,ii)                         ! must be a reciprocal lattice vector
    qq=[dot_product(qq,v2a(r1)),dot_product(qq,v2a(r2)),dot_product(qq,v2a(r3))]/(2*pi)
    if (maxval(abs(qq-nint(qq))).gt.1d-6) nbad=nbad+1
 enddo
 do ii=1,nibz
    call get_k_info(kibz(:,ii),NC,nn,a1,a2,a3,g1,g2,g3,ins)
    if (nn.ne.mapinv(ii)) nbad=nbad+1
 enddo
 write(ulog,*)'CHECK_MESH_LABELS: k, k+G, -k and IBZ points found at their indices; errors=',nbad
 if (nbad.ne.0) error stop 'CHECK_MESH_LABELS: mesh and get_k_info disagree'
 end subroutine check_mesh_labels

 subroutine timing(what)
 character(*), intent(in) :: what
 call cpu_time(cputim)
 call system_clock(wclock1)
 write(utimes,'(a,a,a,a,a)') what,',  CPU ',hms(dble(cputim)),'   WALL ',hms(dble(wclock1-wclock0)/wrate)
 end subroutine timing

!---------------------------------------------------------------
 character(len=13) function hms(t)
!! seconds -> hhh:mm:ss.s
 real(8), intent(in) :: t
 integer h,m
 real(8) sec
 h=int(t/3600) ; m=int((t-3600d0*h)/60) ; sec=t-3600d0*h-60d0*m
 write(hms,'(i3,a,i2.2,a,f4.1)') h,':',m,':',sec
 if (sec.lt.9.95d0) hms(8:8)='0'
 end function hms

!---------------------------------------------------------------
 subroutine get_cwd(path)
 character(*), intent(out) :: path
 character(kind=c_char) :: cbuf(400)
 type(c_ptr) :: p
 integer n
 p = c_getcwd(cbuf,int(size(cbuf),c_size_t))
 path=''
 do n=1,size(cbuf)
    if (cbuf(n).eq.c_null_char) exit
    path(n:n)=cbuf(n)
 enddo
 end subroutine get_cwd

!---------------------------------------------------------------
 subroutine enter_rank_dir()
!! ranks > 0 run in .mpi_rank_NNN, with links to the input files of the run folder
 character(len=*), parameter :: inputs(7) = [character(len=17) :: 'structure.params', &
     'latdyn.params','dielectric.params','kpbs.params','lat_fc.dat','fc2.dat','fc3.dat']
 integer n,st
 write(rankdir,'(a,i3.3)') '.mpi_rank_',rank
 st = c_mkdir(trim(rankdir)//c_null_char, int(o'755',c_int))
 do n=1,size(inputs)
    st = c_symlink(trim(runpath)//'/'//trim(inputs(n))//c_null_char,  &
                   trim(rankdir)//'/'//trim(inputs(n))//c_null_char)
 enddo
 st = c_symlink(trim(runpath)//'/fc4.dat'//c_null_char, trim(rankdir)//'/fc4.dat'//c_null_char)
 if (c_chdir(trim(rankdir)//c_null_char).ne.0) then
    write(*,*)'rank ',rank,': cannot enter ',trim(rankdir)
    call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
 endif
 end subroutine enter_rank_dir

!---------------------------------------------------------------
 subroutine split_range(n,r1,r2)
!! the block of 1..n done by this rank (the first mod(n,nprocs) ranks get one more)
 integer, intent(in) :: n
 integer, intent(out) :: r1,r2
 integer q,m
 q=n/nprocs ; m=mod(n,nprocs)
 r1 = rank*q + min(rank,m) + 1
 r2 = r1 + q - 1 + merge(1,0,rank.lt.m)
 end subroutine split_range

!---------------------------------------------------------------
 subroutine gather_all(nper,loc,glob)
!! glob(1:nper*ntotal) on every rank = the ranks' loc(1:nper*nlocal) in rank order
 integer, intent(in) :: nper
 real(8), intent(in) :: loc(:)
 real(8), intent(out) :: glob(:)
 integer mine,r
 if (.not.allocated(counts)) allocate(counts(nprocs),displs(nprocs))
 mine=size(loc)
 call MPI_ALLGATHER(mine,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
 displs(1)=0
 do r=2,nprocs
    displs(r)=displs(r-1)+counts(r-1)
 enddo
 call MPI_ALLGATHERV(loc,mine,MPI_DOUBLE_PRECISION,glob,counts,displs,MPI_DOUBLE_PRECISION,  &
                     MPI_COMM_WORLD,ierr)
 end subroutine gather_all

!---------------------------------------------------------------
 subroutine gather_2d(a_loc,a_glob)
!! a_glob(nibz,n2) on every rank from the ranks' a_loc(nloc,n2) (rows = their IBZ block)
 real(8), intent(in) :: a_loc(:,:)
 real(8), intent(out) :: a_glob(:,:)
 real(8), allocatable :: b1(:),b2(:)
 integer n2
 n2=size(a_loc,2)
 allocate(b1(size(a_loc)),b2(size(a_glob)))
 b1=reshape(transpose(a_loc),[size(a_loc)])
 call gather_all(n2,b1,b2)
 a_glob=transpose(reshape(b2,[n2,size(a_glob,1)]))
 deallocate(b1,b2)
 end subroutine gather_2d

!---------------------------------------------------------------
 subroutine linewidths_on_band_path()
!! job 1: Im(self-energy) of every mode along the band path, for every temperature
 integer nb,b1,b2,lwse,lwse2,ll
 real(8), allocatable :: evlbs(:,:),vlcbs(:,:,:),om0(:),ise_loc(:,:),ise(:,:)
 complex(8), allocatable :: evcbs(:,:,:)
 complex(8) nselfs,uselfs
 real(8) kbss(3)
 lwse=531 ; lwse2=553 ; ll=532
 call make_kp_bs
 nb=nkp_bs
 allocate(evlbs(ndyn,nb),evcbs(ndyn,ndyn,nb),vlcbs(3,ndyn,nb),om0(ndyn))
 call get_frequencies(nb,kp_bs,dk_bs,ndyn,evlbs,ndyn,evcbs,ll,vlcbs,0)   ! k itself: V3 sums
 call split_range(nb,b1,b2)
 allocate(ise_loc(max(b2-b1+1,0),ndyn),ise(nb,ndyn))
 if (root) then
    open(lwse,file='lw-selfenergy.dat')
    open(lwse2,file='lw-selfenergy2.dat')
    write(lwse,*)'# T(K), band, k index, k(3), dk, omega(1/cm), Im Sigma(1/cm), 2 Im Sigma, tau(ps), 1/tau'
    write(lwse2,*)'# T(K), dk, omega(1/cm) of each band, then 2 Im Sigma (1/cm) of each band'
 endif
 do itmp=1,ntemp
    tempk=temperature(itmp)
    temp = tempk*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1
    do k=b1,b2
       kbss=kp_bs(:,k)
       do la=1,ndyn
          om0(la)=sqrt(abs(evlbs(la,k))) * cnst
          call function_self_w3(kbss,la,om0(la),temp,nselfs,uselfs)
          ise_loc(k-b1+1,la)=2*aimag(nselfs+uselfs)
       enddo
    enddo
    call gather_2d(ise_loc,ise)
    if (root) then
       do k=1,nb
          do la=1,ndyn
             om0(la)=sqrt(abs(evlbs(la,k))) * cnst
             write(lwse,'(f9.2,2i6,5(1x,f9.4),9(1x,g12.5))')tempk,la,k,kp_bs(:,k),dk_bs(k),om0(la),  &
                 ise(k,la)/2,ise(k,la),1d10/c_light/ise(k,la),c_light*ise(k,la)/1d10
          enddo
          write(lwse2,'(f9.2,1x,f9.4,999(1x,g12.5))')tempk,dk_bs(k),om0,ise(k,:)
       enddo
    endif
    call timing('linewidths along the band path')
 enddo
 if (root) then
    close(lwse)
    close(lwse2)
 endif
 end subroutine linewidths_on_band_path

!---------------------------------------------------------------
 real(8) function temperature(it)
 integer, intent(in) :: it
 if (ntemp.le.1) then
    temperature=tmin
 else
    temperature=tmin+(tmax-tmin)*(it-1)/dble(ntemp-1)
 endif
 end function temperature

!---------------------------------------------------------------
 character(len=440) function v3file(kk)
 integer, intent(in) :: kk
 write(v3file,'(a,a,i6.6,a)') trim(v3dir),'k',kk,'.bin'
 end function v3file

!---------------------------------------------------------------
 logical function read_v3(kk)
!! reads P1smpi(1,...),P2smpi(1,...) of IBZ point kk from V3DELTA if it was made on this mesh: the sparse
!! files of 2.1 (7-integer header ending in -2) or the dense ones of 2.0 (6 integers, then both arrays)
 integer, intent(in) :: kk
 integer u,hdr(6),st,fmt
 integer(8) fsize,dense
 logical ex
 read_v3=.false.
 inquire(file=trim(v3file(kk)),exist=ex,size=fsize)
 if (.not.ex) return
 open(newunit=u,file=trim(v3file(kk)),form='unformatted',access='stream',status='old',iostat=st)
 if (st.ne.0) return
 read(u,iostat=st) hdr
 if (st.eq.0 .and. all(hdr.eq.[nc(1),nc(2),nc(3),ndyn,nibz,kk])) then
    dense=24+16*int(nkc,8)*int(ndyn,8)**3
    if (fsize.eq.dense) then
       read(u,iostat=st) P1smpi(1,:,:,:,:),P2smpi(1,:,:,:,:)
       read_v3 = (st.eq.0)
    else
       read(u,iostat=st) fmt
       if (st.eq.0 .and. fmt.eq.-2) read_v3 = read_sparse(u,P1smpi) .and. read_sparse(u,P2smpi)
    endif
 endif
 close(u)
 end function read_v3

!---------------------------------------------------------------
 subroutine write_v3(kk)
!! 2.1: only the nonzero elements (about 2 %): header [nc, ndyn, nibz, kk, -2], then for P1smpi and P2smpi
!! the count (integer 8), the positions in P(1,:,:,:,:) (integer 8, column order) and the values (real 8)
 integer, intent(in) :: kk
 integer u
 open(newunit=u,file=trim(v3file(kk)),form='unformatted',access='stream',status='replace')
 write(u) [nc(1),nc(2),nc(3),ndyn,nibz,kk,-2]
 call write_sparse(u,P1smpi)
 call write_sparse(u,P2smpi)
 close(u)
 end subroutine write_v3

!---------------------------------------------------------------
 subroutine write_sparse(u,p)
 integer, intent(in) :: u
 real(8), intent(in) :: p(:,:,:,:,:)
 integer(8) m,lin
 integer w,a,b,c
 integer(8), allocatable :: idx(:)
 real(8), allocatable :: val(:)
 m=count(p(1,:,:,:,:).ne.0d0)
 allocate(idx(m),val(m))
 m=0 ; lin=0
 do c=1,ndyn
 do b=1,ndyn
 do a=1,ndyn
 do w=1,nkc
    lin=lin+1
    if (p(1,w,a,b,c).ne.0d0) then
       m=m+1 ; idx(m)=lin ; val(m)=p(1,w,a,b,c)
    endif
 enddo
 enddo
 enddo
 enddo
 write(u) m
 if (m.gt.0) write(u) idx,val
 end subroutine write_sparse

!---------------------------------------------------------------
 logical function read_sparse(u,p)
 integer, intent(in) :: u
 real(8), intent(inout) :: p(:,:,:,:,:)
 integer(8) m,k,t
 integer st,w,a,b,c
 integer(8), allocatable :: idx(:)
 real(8), allocatable :: val(:)
 read_sparse=.false.
 read(u,iostat=st) m
 if (st.ne.0 .or. m.lt.0) return
 allocate(idx(m),val(m))
 if (m.gt.0) then
    read(u,iostat=st) idx,val
    if (st.ne.0) return
 endif
 p(1,:,:,:,:)=0
 do k=1,m
    t=idx(k)-1
    w=int(mod(t,int(nkc,8)))+1 ; t=t/nkc
    a=int(mod(t,int(ndyn,8)))+1 ; t=t/ndyn
    b=int(mod(t,int(ndyn,8)))+1 ; c=int(t/ndyn)+1
    p(1,w,a,b,c)=val(k)
 enddo
 read_sparse=.true.
 end function read_sparse

!---------------------------------------------------------------
 subroutine read_whole(f,txt,ex)
!! the bytes of file f (ex = .false. if it does not exist)
 character(len=*), intent(in) :: f
 character(len=:), allocatable, intent(out) :: txt
 logical, intent(out) :: ex
 integer u,n
 inquire(file=f,exist=ex,size=n)
 if (.not.ex) then
    allocate(character(len=0) :: txt)
    return
 endif
 allocate(character(len=max(n,0)) :: txt)
 open(newunit=u,file=f,form='unformatted',access='stream',status='old')
 if (n.gt.0) read(u) txt
 close(u)
 end subroutine read_whole

!---------------------------------------------------------------
 subroutine mark_gamma_acoustic()
!! gamma_acoustic(la,nk): the 3 acoustic modes at Gamma of the mesh, the rigid translations
!! (translation_weight > translation_min, module phonon_rules; found from the eigenvector, not the
!! band index). omega -> 0 makes their V3 (1/sqrt(omega1 omega2 omega3)) and occupation diverge;
!! they carry no heat and have zero weight in the k-sum, so they are excluded from every sum:
!! V3 = 0, v = 0, F = 0. Their omega^2 is kept as computed.
 use phonon_rules, only : translation_weight, translation_min
 integer kk,l,n
 n=0
 do kk=1,nkc
    if (sqrt(sum(kpc(:,kk)**2)).gt.1d-10) cycle
    do l=1,ndyn
       if (translation_weight(natoms0,atom0(1:natoms0)%mass,eigenvec(:,l,kk)).gt.translation_min) then
          gamma_acoustic(l,kk)=.true.
          veloc(:,l,kk)=0
          n=n+1
          write(ulog,'(a,i4,a,g12.4)')' MARK_GAMMA_ACOUSTIC: Gamma, mode',l,' excluded from the sums, omega^2=', &
 &              eigenval(l,kk)
       endif
    enddo
 enddo
 if (n.ne.3) then
    write(ulog,*)' MARK_GAMMA_ACOUSTIC: found',n,' rigid translations at Gamma, expected 3'
    error stop 'MARK_GAMMA_ACOUSTIC: not 3 acoustic modes at Gamma (is Gamma on the mesh? force constants?)'
 endif
 end subroutine mark_gamma_acoustic

!---------------------------------------------------------------
 subroutine check_v3_inputs()
!! V3DELTA keeps copies (input_<name>) of the files V3*delta depends on. If any present file
!! differs from its copy (or one exists and the other not), the k######.bin files were made
!! with other force constants or another structure: they are removed and the copies renewed.
 character(len=*), parameter :: inputs(5) = [character(len=17) :: 'lat_fc.dat','fc2.dat', &
     'fc3.dat','structure.params','dielectric.params']
 character(len=:), allocatable :: a,b
 logical ea,eb,same,old
 integer n,u,kk
 same=.true.
 do n=1,size(inputs)
    call read_whole(trim(inputs(n)),a,ea)
    call read_whole(trim(v3dir)//'input_'//trim(inputs(n)),b,eb)
    if ((ea.neqv.eb) .or. (a.ne.b) .or. (len(a).ne.len(b))) same=.false.
 enddo
 if (same) then
    write(ulog,*)' V3DELTA: made from the same lat_fc/fc2/fc3/structure/dielectric files; its points are reused if on this mesh'
    return
 endif
 old=.false.
 do kk=1,nibz
    inquire(file=trim(v3file(kk)),exist=ea)
    if (ea) then
       old=.true.
       open(newunit=u,file=trim(v3file(kk)),status='old')
       close(u,status='delete')
    endif
 enddo
 if (old) then
    write(ulog,*)' V3DELTA: made from other lat_fc/fc2/fc3/structure/dielectric files: removed, recomputed'
    write(*,*)' V3DELTA was made from other lat_fc/fc2/fc3/structure/dielectric files: recomputed'
 endif
 do n=1,size(inputs)
    call read_whole(trim(inputs(n)),a,ea)
    if (ea) then
       open(newunit=u,file=trim(v3dir)//'input_'//trim(inputs(n)),form='unformatted', &
            access='stream',status='replace')
       write(u) a
       close(u)
    else
       inquire(file=trim(v3dir)//'input_'//trim(inputs(n)),exist=eb)
       if (eb) then
          open(newunit=u,file=trim(v3dir)//'input_'//trim(inputs(n)),status='old')
          close(u,status='delete')
       endif
    endif
 enddo
 end subroutine check_v3_inputs

!---------------------------------------------------------------
 subroutine v3_delta_ibz()
!! |V3|^2 * delta (tetrahedra), temperature independent: P1smpi (coalescence), P2smpi (decay)
 integer nv3o,nnew,nold
 call split_range(nibz,kglob(1),kglob(2))
 nloc=max(kglob(2)-kglob(1)+1,0)
 write(ulog,*)' this rank does IBZ points ',kglob,' of ',nibz
 call allocate_iter(nibz,nkc,ndyn)
 call cal_freq(nkc,ndyn)
 if (v3_store.eq.1) then                               ! all of this rank's points in memory
    allocate(v33s8(nloc,nkc,ndyn,ndyn,ndyn))
    allocate(P1smpi(nloc,nkc,ndyn,ndyn,ndyn),P2smpi(nloc,nkc,ndyn,ndyn,ndyn))
    P1smpi=0 ; P2smpi=0
    if (nloc.gt.0) call v3_delta(kglob)
    deallocate(v33s8)
    call timing('V3*delta in memory')
 else                                                  ! one file per IBZ point, reused if present
    if (root) ierr = c_mkdir(trim(v3dir)//c_null_char, int(o'755',c_int))
    if (root) call check_v3_inputs()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    allocate(v33s8(1,nkc,ndyn,ndyn,ndyn))
    allocate(P1smpi(1,nkc,ndyn,ndyn,ndyn),P2smpi(1,nkc,ndyn,ndyn,ndyn))
    nnew=0 ; nold=0
    do k=kglob(1),kglob(2)
       if (read_v3(k)) then
          nold=nold+1
          cycle
       endif
       ksub=[k,k]
       P1smpi=0 ; P2smpi=0
       call v3_delta(ksub)
       call write_v3(k)
       nnew=nnew+1
    enddo
    deallocate(v33s8)
    write(ulog,*)' V3*delta: ',nnew,' IBZ points computed, ',nold,' read from ',trim(v3dir)
    call timing('V3*delta on disk')
 endif
 end subroutine v3_delta_ibz

!---------------------------------------------------------------
 subroutine v3_delta(kr)
!! P1smpi, P2smpi of the IBZ points kr (2.1, 2026-11-09). Before 2.1 every |V3|^2 of the point was
!! computed (nkc*ndyn^3 matrix elements) and then weighted by the tetrahedron delta, which is zero for
!! about 98 % of them (Cu2O 18^3: 0.9-1.5 % of P1, P2 nonzero). Now: (1) the tetrahedron weights alone
!! (calculate_FGR_mpi2 with |V3|^2 = 1; they depend on the frequencies only); (2) v3_need = the elements
!! these weights need, closed over the degenerate sets; (3) only those |V3|^2; (4) their degenerate
!! average as before; (5) P = weight x |V3|^2 as calculate_FGR_mpi2 forms it; (6) the degenerate average
!! of P as before. Same P up to rounding: (sum_t c_t/6) F instead of sum_t (c_t/6 F).
 integer, intent(in) :: kr(2)
 integer nv3w,nv3c
 nv3w=(kr(2)-kr(1)+1)*nkc*ndyn**3
 v3_weights_only=.true.
 call calculate_FGR_mpi2(kr,ndyn,nv3w)                  ! (1) weights in P1smpi, P2smpi
 v3_weights_only=.false.
 call v3_needed(kr)                                     ! (2)
 write(ulog,'(a,2i7,a,f8.4,a)')' V3*delta: IBZ points ',kr,': ',100*count(v3_need)/dble(size(v3_need)), &
      ' % of the |V3|^2 elements needed'
 call calculate_w3_ibz_split_sq(kr,nv3w,ndyn,nkc,kpc,eigenval,eigenvec,nv3c)   ! (3)
 deallocate(v3_need)
 call average_degenerate_v3(kr)                         ! (4)
 call apply_v3_weights(kr)                              ! (5)
 call average_degenerate_p(kr)                          ! (6)
 end subroutine v3_delta

!---------------------------------------------------------------
 subroutine v3_needed(kr)
!! v3_need(s,n2,l1,l2,l3): the |V3|^2 elements the weights in P1smpi, P2smpi multiply. Process 1 of
!! calculate_FGR_mpi2 puts the element of q2 = w at P1(w), process 2 that of q2 = -w at P2(w)
!! (F_arbi(nk2n) with nk2n the index of -q2). average_degenerate_v3 mixes the partners of each
!! degenerate set at q1, q2 and q3 = -q1-q2, so the partners of a needed element are needed too.
 use lattice, only : g1,g2,g3
 integer, intent(in) :: kr(2)
 integer s,n1,n2,i1,i3,ia,ib,ic,ins,ng,w,l1,l2,l3
 integer g(ndyn),mw(nkc)
 real(8) q(3)
 if (allocated(v3_need)) deallocate(v3_need)
 allocate(v3_need(kr(2)-kr(1)+1,nkc,ndyn,ndyn,ndyn))
 v3_need=.false.
 do w=1,nkc
    q=-kpc(:,w)
    call get_k_info(q,NC,mw(w),ia,ib,ic,g1,g2,g3,ins)
 enddo
 s=0
 do n1=kr(1),kr(2)
    s=s+1
    do l3=1,ndyn
    do l2=1,ndyn
    do l1=1,ndyn
    do w=1,nkc
       if (P1smpi(s,w,l1,l2,l3).ne.0d0) v3_need(s,w,l1,l2,l3)=.true.
       if (P2smpi(s,w,l1,l2,l3).ne.0d0) v3_need(s,mw(w),l1,l2,l3)=.true.
    enddo
    enddo
    enddo
    enddo
    i1=mapinv(n1)
    call degenerate_sets(i1,g,ng)
    call close_axis(s,0,1,g,ng)
    do n2=1,nkc
       call degenerate_sets(n2,g,ng)
       call close_axis(s,n2,2,g,ng)
       q=-kpc(:,i1)-kpc(:,n2)
       call get_k_info(q,NC,i3,ia,ib,ic,g1,g2,g3,ins)
       call degenerate_sets(i3,g,ng)
       call close_axis(s,n2,3,g,ng)
    enddo
 enddo
 end subroutine v3_needed

!---------------------------------------------------------------
 subroutine close_axis(s,n2,ax,grp,n)
!! v3_need(s,n2,...) true for every member of a degenerate set along axis ax (1: l1, 2: l2, 3: l3) when it
!! is true for one of them; n2 = 0: every q2
 integer, intent(in) :: s,n2,ax,grp(ndyn),n
 integer k,a,b,m,qa,qb,q2
 logical(1) any1
 if (n.eq.ndyn) return                                  ! no degenerate set
 if (n2.eq.0) then
    qa=1 ; qb=nkc
 else
    qa=n2 ; qb=n2
 endif
 do k=1,n
    if (count(grp.eq.k).lt.2) cycle
    do q2=qa,qb
    do b=1,ndyn
    do a=1,ndyn
       any1=.false.
       do m=1,ndyn
          if (grp(m).ne.k) cycle
          select case (ax)
          case (1) ; any1 = any1 .or. v3_need(s,q2,m,a,b)
          case (2) ; any1 = any1 .or. v3_need(s,q2,a,m,b)
          case (3) ; any1 = any1 .or. v3_need(s,q2,a,b,m)
          end select
       enddo
       if (.not.any1) cycle
       do m=1,ndyn
          if (grp(m).ne.k) cycle
          select case (ax)
          case (1) ; v3_need(s,q2,m,a,b)=.true.
          case (2) ; v3_need(s,q2,a,m,b)=.true.
          case (3) ; v3_need(s,q2,a,b,m)=.true.
          end select
       enddo
    enddo
    enddo
    enddo
 enddo
 end subroutine close_axis

!---------------------------------------------------------------
 subroutine apply_v3_weights(kr)
!! P1smpi, P2smpi (the tetrahedron weights) times |V3|^2 as calculate_FGR_mpi2 forms it from v33s8:
!! x (2 pi)^2/nkc, and 0 for an acoustic mode (index <= 3) of q2 or q3 at Gamma (index NC/2+1 in
!! every direction). Process 1: q2 = kpc(w), q3 = q1+q2, into P1(w); process 2: q2 = kpc(j), Gamma
!! tested on -q2, q3 = q1+q2, into P2(index of -q2). Same code as there.
 use lattice, only : g1,g2,g3
 integer, intent(in) :: kr(2)
 integer s,n1,i,j,ii,jj,kk,ins,i2,j2,k2,i3,j3,k3,nk2,nk3
 integer nk2n(nkc)
 logical gam2(nkc),gam2n(nkc),gam3(nkc)
 real(8) q1(3),q2(3),f
 s=0
 do n1=kr(1),kr(2)
    s=s+1
    i=mapinv(n1)
    q1=kpc(:,i)
    do j=1,nkc
       q2=kpc(:,j)
       call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,ins)
       gam2(j) = i2.eq.NC(1)/2+1 .and. j2.eq.NC(2)/2+1 .and. k2.eq.NC(3)/2+1
       call get_k_info(-q2,NC,nk2n(j),i2,j2,k2,g1,g2,g3,ins)
       gam2n(j) = i2.eq.NC(1)/2+1 .and. j2.eq.NC(2)/2+1 .and. k2.eq.NC(3)/2+1
       call get_k_info(q1+q2,NC,nk3,i3,j3,k3,g1,g2,g3,ins)
       gam3(j) = i3.eq.NC(1)/2+1 .and. j3.eq.NC(2)/2+1 .and. k3.eq.NC(3)/2+1
    enddo
    do kk=1,ndyn
    do jj=1,ndyn
    do ii=1,ndyn
    do j=1,nkc
       f=v33s8(s,j,ii,jj,kk) * (2*pi)**2 / nkc
       if ((jj.le.3 .and. gam2(j)) .or. (kk.le.3 .and. gam3(j))) then
          P1smpi(s,j,ii,jj,kk)=0
       else
          P1smpi(s,j,ii,jj,kk)=P1smpi(s,j,ii,jj,kk)*f
       endif
       if ((jj.le.3 .and. gam2n(j)) .or. (kk.le.3 .and. gam3(j))) then
          P2smpi(s,nk2n(j),ii,jj,kk)=0
       else
          P2smpi(s,nk2n(j),ii,jj,kk)=P2smpi(s,nk2n(j),ii,jj,kk)*f
       endif
    enddo
    enddo
    enddo
    enddo
 enddo
 end subroutine apply_v3_weights

!---------------------------------------------------------------
 subroutine thermal_rates(it_first)
!! for the present temperature: P1,P2 (= P*smpi times occupations), Q of this rank's points
!! and (job 2) its rows of the collision matrix
 logical, intent(in) :: it_first
 integer kk
 if (it_first) then
    allocate(Qloc(nloc,ndyn),QNloc(nloc,ndyn),QUloc(nloc,ndyn))
    if (job.eq.2) allocate(CMloc(nloc*ndyn*3,nibz*ndyn*3))
 endif
 if (allocated(P1)) deallocate(P1,P2,Piso)
 if (allocated(Qvaluempi)) deallocate(Qvaluempi,Qvalue_Nmpi,Qvalue_Umpi)
 if (allocated(Collision_Matrixmpi)) deallocate(Collision_Matrixmpi)
 if (v3_store.eq.1) then
    call allocate_FGR(nloc,nkc,ndyn,ndyn,ndyn)
    allocate(Qvaluempi(nloc,ndyn),Qvalue_Nmpi(nloc,ndyn),Qvalue_Umpi(nloc,ndyn))
    if (nloc.gt.0) then
       call calculate_tauRTA_mpi2(kglob,ndyn,tempk)
       if (job.eq.2) then
          allocate(Collision_Matrixmpi(nloc*ndyn*3,nibz*ndyn*3))
          call CollisionMatrix_mpi(kglob,ndyn)
          CMloc=Collision_Matrixmpi
       endif
    endif
    Qloc=Qvaluempi ; QNloc=Qvalue_Nmpi ; QUloc=Qvalue_Umpi
 else
    call allocate_FGR(1,nkc,ndyn,ndyn,ndyn)
    allocate(Qvaluempi(1,ndyn),Qvalue_Nmpi(1,ndyn),Qvalue_Umpi(1,ndyn))
    if (job.eq.2) allocate(Collision_Matrixmpi(ndyn*3,nibz*ndyn*3))
    do kk=kglob(1),kglob(2)
       kl=kk-kglob(1)+1
       if (.not.read_v3(kk)) then
          write(*,*)'rank ',rank,': cannot read ',trim(v3file(kk))
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
       endif
       ksub=[kk,kk]
       call calculate_tauRTA_mpi2(ksub,ndyn,tempk)
       Qloc(kl,:)=Qvaluempi(1,:) ; QNloc(kl,:)=Qvalue_Nmpi(1,:) ; QUloc(kl,:)=Qvalue_Umpi(1,:)
       if (job.eq.2) then
          call CollisionMatrix_mpi(ksub,ndyn)
          CMloc((kl-1)*ndyn*3+1:kl*ndyn*3,:)=Collision_Matrixmpi
       endif
    enddo
 endif
 end subroutine thermal_rates

!---------------------------------------------------------------
 subroutine total_kappa(kt)
 real(8), intent(out) :: kt(3,3)
 kt=sum(kappa,dim=1)
 end subroutine total_kappa

!---------------------------------------------------------------
 subroutine kappa_all_temperatures()
!! for each temperature: RTA, then the direct (job 2) or iterative (job 3) solution
 integer n
 real(8), allocatable :: CMt(:,:),CMg(:,:),Ax(:)
 n=nibz*ndyn*3
 allocate(FFF(nibz,ndyn,3))
 allocate(kappa_k(nibz,ndyn,3,3),cR_k(nibz,ndyn,3,3),wP_k(nibz,ndyn,3,3),wC_k(nibz,ndyn,3,3),tauinvR(nibz,ndyn))
 if (root) then
    open(newunit=usum,file='kappa_summary.dat',status='replace')
    write(usum,'(a)')'# kappa (W/mK) versus T: average (xx+yy+zz)/3, then xx yy zz, for RTA = kappa_P (the RTA'// &
 &   ' with degenerate modes as blocks), the full solution (job 2 direct, job 3 iterative; same block form)'
    write(usum,'(a)')'# and Wigner = kappa_P + kappa_C (Simoncelli-Marzari-Mauri). Per mode and cumulative: kappa_modes_T*.dat,'// &
 &   ' kappa_cumul_mfp_T*.dat, kappa_cumul_freq_T*.dat'
    write(usum,'(a)')'# T(K)  RTA: avg xx yy zz   full: avg xx yy zz   Wigner: avg xx yy zz'
    open(newunit=uwm,file='kappa_wigner_modes.dat',status='replace')
    write(uwm,'(a)')'# Wigner kappa (W/mK) summed over the lowest bands at every q (TA1 = band 1, TA2 = band 2,'// &
 &   ' LA = band 3; by frequency): band 1, bands 1-2, bands 1-3, all bands (= the total)'
    write(uwm,'(a)')'# kappa_C of a pair of modes counted half to each. Columns: T(K), then for each of the 4: avg xx yy zz'
    write(uwm,'(a)')'# T(K)  TA1: avg xx yy zz   TA1+TA2: avg xx yy zz   TA1+TA2+LA: avg xx yy zz   all: avg xx yy zz'
    tauf=3450 ; kapRTA=6450 ; ukvt=6460
    open(tauf,file='tau_RTA.dat',status='unknown')
    open(kapRTA,file='kappa_RTA.dat',status='unknown')
    open(ukvt,file='kappa_vs_T.dat',status='unknown')
    ukw=6480
    open(ukw,file='kappa_wigner.dat',status='unknown')
    write(ukw,'(a)')'# Wigner (Simoncelli-Marzari-Mauri), RTA linewidths. kappa_P: pairs of degenerate modes'// &
 &                  ' (incl. s=s'', trace form); kappa_C: coherences; W/mK, xx yy zz xy xz yz'
    write(ukw,'(a)')'# T(K)  kappa_RTA(6)  kappa_P(6)  kappa_C(6)  kappa_P+kappa_C(6)'
    udist=6470
    open(udist,file='distribution_F.dat',status='unknown')
    write(udist,'(a)')'# FFF(ibz,mode,xyz) of sub.f90; RTA: -v hbar omega n0(n0+1)/(k_B T^2 Q). One block per T and solution'
    write(udist,'(a)')'# ibz  mode   F_x  F_y  F_z'
    if (job.eq.2) then
       tauef=7450 ; kapd=8450
       open(tauef,file='tau_RTAandeff.dat',status='unknown')
       open(kapd,file='kappa_direct.dat',status='unknown')
       write(ukvt,'(a)')'# degenerate modes as blocks: kappa_RTA = kappa_P (Wigner, trace form), kappa_direct = per-mode'
       write(ukvt,'(a)')'# value + kappa_P - per-mode RTA; then the per-mode values (basis-dependent at degeneracies)'
       write(ukvt,'(a)')'# T(K)  kappa_RTA(6)  kappa_direct(6)  |Ax-b|/|b|  per-mode: kappa_RTA(6) kappa_direct(6); xx yy zz xy xz yz, W/mK'
    else
       tauefi=9450 ; kappai=4451
       open(tauefi,file='tau_RTAandeff_iter.dat',status='unknown')
       open(kappai,file='kappa_iter.dat',status='unknown')
       write(ukvt,'(a)')'# degenerate modes as blocks: kappa_RTA = kappa_P (Wigner, trace form), kappa_iterative = per-mode'
       write(ukvt,'(a)')'# value + kappa_P - per-mode RTA; then the per-mode values (basis-dependent at degeneracies)'
       write(ukvt,'(a)')'# T(K)  kappa_RTA(6)  kappa_iterative(6)  iterations  per-mode: kappa_RTA(6) kappa_iterative(6); xx yy zz xy xz yz, W/mK'
    endif
 endif
 if (job.eq.2) allocate(RHS(n),DirectSolution(n))

 do itmp=1,ntemp
    tempk=temperature(itmp)
    write(ulog,*)'===================== T(K)=',tempk
    call distribution(nkc,ndyn,tempk)        ! Bose-Einstein occupations at this T
    call thermal_rates(itmp.eq.1)
    call gather_2d(Qloc,Qvalue)
    call gather_2d(QNloc,Qvalue_N)
    call gather_2d(QUloc,Qvalue_U)
    call Qandtau(ndyn,tempk)                 ! tau_N, tau_U, tau_total
    call calculate_RTA_distributionfunction(nibz,nkc,kibz,ndyn,tempk,veloc,eigenval)  ! FFF = F_RTA
    if (.not.allocated(F_RTA)) allocate(F_RTA(nibz,ndyn,3))
    F_RTA=FFF                                ! for tau_eff in mode_kappa (2.2)
    call mode_kappa(ndyn,tempk,veloc,eigenval)
    call total_kappa(kRTA)
    cR_k=kappa_k ; tauinvR=tauinv_tot
    if (root) call write_kappa(ndyn,tempk,veloc,eigenval,kapRTA,tauf)
    if (root) call write_F('RTA')
    if (root) then
       call wigner_kappa(kWP,kWC)
       write(ukw,'(f8.2,24(1x,es11.4))')tempk,(kRTA(i,i),i=1,3),kRTA(1,2),kRTA(1,3),kRTA(2,3),        &
 &        (kWP(i,i),i=1,3),kWP(1,2),kWP(1,3),kWP(2,3),(kWC(i,i),i=1,3),kWC(1,2),kWC(1,3),kWC(2,3), &
 &        (kWP(i,i)+kWC(i,i),i=1,3),kWP(1,2)+kWC(1,2),kWP(1,3)+kWC(1,3),kWP(2,3)+kWC(2,3)
       flush(ukw)
       call timing('Wigner')
    endif
    call timing('RTA')

    if (job.eq.2) then
       ! collision matrix rows from all ranks to rank 0, which solves M F = RHS
       allocate(CMt(n,nloc*ndyn*3))
       CMt=transpose(CMloc)
       if (root) then
          allocate(CMg(n,n))
       else
          allocate(CMg(1,1))
       endif
       call gather_rows(CMt,CMg)
       deallocate(CMt)
       if (.not.root) deallocate(CMg)
       if (root) then
          allocate(Collision_Matrix(n,n),Ax(n))
          Collision_Matrix=transpose(CMg)
          deallocate(CMg)
          call cal_RHS_IBZ(ndyn,nkc,tempk,veloc,eigenval)
          s=0
          do is=1,nibz ; do js=1,ndyn ; do ks=1,3
             s=s+1
             DirectSolution(s)=FFF(is,js,ks)       ! RTA solution as the starting point
          enddo ; enddo ; enddo
          if (solver_direct.eq.1) then
             if (.not.allocated(sig)) allocate(sig(n))
             svdcut=1d-9
             call svd_solver(n,n,Collision_Matrix,RHS,DirectSolution,sig,svdcut,errorr,ermax,fm,'svd_output_IBZ.dat')
          else
             call cg_quadratic(n,Collision_Matrix,RHS,DirectSolution,fm,100,2)
          endif
          Ax=matmul(Collision_Matrix,DirectSolution)
          resid=sqrt(sum((Ax-RHS)**2)/max(sum(RHS**2),1d-300))
          write(ulog,*)' direct solution: solver (0 CG, 1 SVD)=',solver_direct,'  |Ax-b|/|b|=',resid
          FFF=reshape(DirectSolution,[nibz,ndyn,3],order=[3,2,1])
          call mode_kappa(ndyn,tempk,veloc,eigenval)
          call total_kappa(kdir)
          call write_kappa(ndyn,tempk,veloc,eigenval,kapd,tauef)
          call write_F('direct')
          kblk=kdir+kWP-kRTA
          write(ukvt,'(f8.2,12(1x,es11.4),1x,es9.2,12(1x,es11.4))')tempk,(kWP(i,i),i=1,3),kWP(1,2),kWP(1,3),kWP(2,3), &
               (kblk(i,i),i=1,3),kblk(1,2),kblk(1,3),kblk(2,3),resid, &
               (kRTA(i,i),i=1,3),kRTA(1,2),kRTA(1,3),kRTA(2,3),(kdir(i,i),i=1,3),kdir(1,2),kdir(1,3),kdir(2,3)
          flush(ukvt)
          call write_modes()
          deallocate(Collision_Matrix,Ax)
       endif
       call timing('direct solution')
    else
       call iterate_bte(nit)
       call total_kappa(kdir)
       if (root) then
          call write_F('iterative')
          kblk=kdir+kWP-kRTA
          write(ukvt,'(f8.2,12(1x,es11.4),1x,i5,12(1x,es11.4))')tempk,(kWP(i,i),i=1,3),kWP(1,2),kWP(1,3),kWP(2,3), &
               (kblk(i,i),i=1,3),kblk(1,2),kblk(1,3),kblk(2,3),nit, &
               (kRTA(i,i),i=1,3),kRTA(1,2),kRTA(1,3),kRTA(2,3),(kdir(i,i),i=1,3),kdir(1,2),kdir(1,3),kdir(2,3)
          flush(ukvt)
          call write_modes()
       endif
       call timing('iterative solution')
    endif
 enddo
 if (root) then
    close(tauf) ; close(kapRTA) ; close(ukvt) ; close(udist) ; close(ukw) ; close(usum) ; close(uwm)
    if (job.eq.2) then
       close(tauef) ; close(kapd)
    else
       close(tauefi) ; close(kappai)
    endif
 endif
 end subroutine kappa_all_temperatures

!---------------------------------------------------------------
 subroutine average_degenerate_v3(kr)
!! |V3(q1 l1, q2 l2, q3 l3)|^2, q3=-q1-q2, averaged over each set of degenerate modes at q1, q2
!! and q3 (modes within degtol).  The eigenvectors of a degenerate set are fixed only up to a unitary rotation (they
!! differ from one compiler or LAPACK build to another); the sum of |V3|^2 over the set does not,
!! and by Schur's lemma the exact rate of every partner is that sum / its size.  The tetrahedron
!! weights are per band index, so without this the result depends on the rotation.
 use lattice, only : g1,g2,g3
 integer, intent(in) :: kr(2)
 integer s,n1,n2,i1,i3,ia,ib,ic,ins,ng
 integer g(ndyn)
 real(8) q3(3)
 s=0
 do n1=kr(1),kr(2)
    s=s+1
    i1=mapinv(n1)
    call degenerate_sets(i1,g,ng)                    ! l1 at q1
    call average_axis(v33s8,s,0,1,g,ng)
    do n2=1,nkc
       call degenerate_sets(n2,g,ng)                 ! l2 at q2
       call average_axis(v33s8,s,n2,2,g,ng)
       q3=-kpc(:,i1)-kpc(:,n2)
       call get_k_info(q3,NC,i3,ia,ib,ic,g1,g2,g3,ins)
       call degenerate_sets(i3,g,ng)                 ! l3 at q3
       call average_axis(v33s8,s,n2,3,g,ng)
    enddo
 enddo
 end subroutine average_degenerate_v3

!---------------------------------------------------------------
 subroutine average_degenerate_p(kr)
!! P1smpi, P2smpi = |V3|^2 x tetrahedron delta weight, averaged over the degenerate sets at q2 and
!! at q3 (q1+q2 for P1, q1-q2 for P2). The delta weights are per band index, so the |V3|^2
!! averaging alone leaves P different among the partners of a set, and the collision term
!! sum_l2 P F_l2 then depends on the basis of the set (Ge 4^3 iterative kappa: 0.3 % between
!! random rotations). Sums over a set, hence every rate Q and the RTA, are unchanged.
 use lattice, only : g1,g2,g3
 integer, intent(in) :: kr(2)
 integer s,n1,n2,i1,i3,ia,ib,ic,ins,ng
 integer g(ndyn)
 real(8) q3(3)
 s=0
 do n1=kr(1),kr(2)
    s=s+1
    i1=mapinv(n1)
    do n2=1,nkc
       call degenerate_sets(n2,g,ng)                 ! l2 at q2
       call average_axis(P1smpi,s,n2,2,g,ng)
       call average_axis(P2smpi,s,n2,2,g,ng)
       q3=kpc(:,i1)+kpc(:,n2)
       call get_k_info(q3,NC,i3,ia,ib,ic,g1,g2,g3,ins)
       call degenerate_sets(i3,g,ng)                 ! l3 at q1+q2
       call average_axis(P1smpi,s,n2,3,g,ng)
       q3=kpc(:,i1)-kpc(:,n2)
       call get_k_info(q3,NC,i3,ia,ib,ic,g1,g2,g3,ins)
       call degenerate_sets(i3,g,ng)                 ! l3 at q1-q2
       call average_axis(P2smpi,s,n2,3,g,ng)
    enddo
 enddo
 end subroutine average_degenerate_p

!---------------------------------------------------------------
 subroutine degenerate_sets(nk,grp,n)
!! grp(l): label of the degenerate set of mode l at FBZ point nk (frequencies ascending)
 integer, intent(in) :: nk
 integer, intent(out) :: grp(ndyn),n
 integer l
 n=1 ; grp(1)=1
 do l=2,ndyn
    if (abs(frequency(nk,l)-frequency(nk,l-1)).gt.degtol) n=n+1
    grp(l)=n
 enddo
 end subroutine degenerate_sets
 subroutine average_axis(arr,s,n2,ax,grp,n)
!! average arr(s,n2,l1,l2,l3) (v33s8, P1smpi or P2smpi) over the modes of each set along axis ax
!! (1: l1, 2: l2, 3: l3); n2=0 means every q2
 real(8), intent(inout) :: arr(:,:,:,:,:)
 integer, intent(in) :: s,n2,ax,grp(ndyn),n
 integer m,cnt,j1,lo,hi,jc,kc
 real(8), allocatable :: a(:,:,:)
 if (n.eq.ndyn) return                            ! no degeneracy
 lo=n2 ; hi=n2
 if (n2.eq.0) then
    lo=1 ; hi=nkc
 endif
 do j1=lo,hi
    do m=1,n
       cnt=count(grp.eq.m)
       if (cnt.eq.1) cycle
       select case (ax)
       case (1)
          allocate(a(1,ndyn,ndyn)) ; a(1,:,:)=0
          do jc=1,ndyn ; if (grp(jc).eq.m) a(1,:,:)=a(1,:,:)+arr(s,j1,jc,:,:) ; enddo
          do jc=1,ndyn ; if (grp(jc).eq.m) arr(s,j1,jc,:,:)=a(1,:,:)/cnt ; enddo
       case (2)
          allocate(a(ndyn,1,ndyn)) ; a(:,1,:)=0
          do jc=1,ndyn ; if (grp(jc).eq.m) a(:,1,:)=a(:,1,:)+arr(s,j1,:,jc,:) ; enddo
          do jc=1,ndyn ; if (grp(jc).eq.m) arr(s,j1,:,jc,:)=a(:,1,:)/cnt ; enddo
       case (3)
          allocate(a(ndyn,ndyn,1)) ; a(:,:,1)=0
          do kc=1,ndyn ; if (grp(kc).eq.m) a(:,:,1)=a(:,:,1)+arr(s,j1,:,:,kc) ; enddo
          do kc=1,ndyn ; if (grp(kc).eq.m) arr(s,j1,:,:,kc)=a(:,:,1)/cnt ; enddo
       end select
       deallocate(a)
    enddo
 enddo
 end subroutine average_axis

!---------------------------------------------------------------
 subroutine wigner_kappa(kP,kC)
!! Wigner thermal conductivity in the relaxation-time approximation, Simoncelli, Marzari and
!! Mauri, Nat. Phys. 15, 809 (2019):
!!   kappa_ab = 1/(N V) sum_q sum_{s,s'} (W_s+W_s')/4 (C_s/W_s + C_s'/W_s') Re(V^a_ss' V^b_s's)
!!              * G/((W_s-W_s')^2 + G^2),   G = (Gamma_s+Gamma_s')/2,
!! W angular frequencies, C_s the modal heat capacity, Gamma_s = 1/tau_s, V the velocity
!! operator (dD/dq)_ss'/(w_s+w_s') in the atom-position phase convention.  The s=s' terms are
!! C v v tau (the RTA); kP collects the pairs of degenerate modes (basis-independent trace
!! form), kC the others (coherences).  Sum over the whole mesh, no symmetrisation needed.
 use lattice, only : volume_r
 real(8), intent(out) :: kP(3,3),kC(3,3)
 complex(8) dyn(ndyn,ndyn),dd(ndyn,ndyn,3),vm(ndyn,ndyn,3),tmp(ndyn,ndyn),tau(ndyn,ndyn)
 real(8) w(ndyn),om(ndyn),cv(ndyn),gam(ndyn),nb,tcm,khat(3),vfac,pre,g,x
 real(8) nbe
 integer nk,ib,s,sp,al,be,a,c,l0
 tcm=tempk*k_b/(100*h_plank*c_light)
 vfac=cnst*1d-10*100*2*pi*c_light                   ! (dD/dk)/w -> m/s, as in get_freq_vg
 kP=0 ; kC=0
 wP_k=0 ; wC_k=0
 do nk=1,nkc
    ib=mapibz(nk)
    call set_dynamical_matrix(kpc(:,nk),dyn,ndyn,dd)
    call add_nonanalytic(kpc(:,nk),dyn,dd,ndyn)
    l0=0
    if (sqrt(sum(kpc(:,nk)**2)).lt.1d-10) l0=3     ! acoustic modes at Gamma
    do s=1,ndyn
       w(s)=sqrt(abs(eigenval(s,nk)))
       om(s)=w(s)*cnst*100*c_light*2*pi
       nb=nbe(w(s)*cnst,tcm,classical)
       cv(s)=(h_plank*w(s)*cnst*100*c_light)**2*nb*(nb+1)/(k_b*tempk**2)
       gam(s)=Qvalue(ib,s)*c_light*100/(nb*(nb+1)) ! 1/tau, as in calculate_RTA_distributionfunction
    enddo
    do al=1,3
       tmp=matmul(transpose(conjg(eigenvec(:,:,nk))),matmul(dd(:,:,al),eigenvec(:,:,nk)))
       tau=0                                       ! <e_s| tau_al |e_s'>, atom positions in the cell
       do a=1,natoms0 ; do c=1,3
          do sp=1,ndyn
             tau(:,sp)=tau(:,sp)+conjg(eigenvec(3*(a-1)+c,:,nk))*atompos(al,a)*eigenvec(3*(a-1)+c,sp,nk)
          enddo
       enddo ; enddo
       do s=1,ndyn ; do sp=1,ndyn
          vm(s,sp,al)=(tmp(s,sp)+ci*(eigenval(s,nk)-eigenval(sp,nk))*tau(s,sp))/(w(s)+w(sp)+1d-30)*vfac
       enddo ; enddo
    enddo
    do s=l0+1,ndyn ; do sp=l0+1,ndyn
       if (gam(s).le.0 .or. gam(sp).le.0) cycle
       g=(gam(s)+gam(sp))/2
       pre=(om(s)+om(sp))/4*(cv(s)/om(s)+cv(sp)/om(sp))*g/((om(s)-om(sp))**2+g*g)
       do al=1,3 ; do be=1,3
          x=pre*real(vm(s,sp,al)*vm(sp,s,be))
          if (abs(w(s)-w(sp))*cnst.lt.degtol) then
             kP(al,be)=kP(al,be)+x
             wP_k(ib,s,al,be)=wP_k(ib,s,al,be)+x/2 ; wP_k(ib,sp,al,be)=wP_k(ib,sp,al,be)+x/2
          else
             kC(al,be)=kC(al,be)+x
             wC_k(ib,s,al,be)=wC_k(ib,s,al,be)+x/2 ; wC_k(ib,sp,al,be)=wC_k(ib,sp,al,be)+x/2
          endif
       enddo ; enddo
    enddo ; enddo
 enddo
 kP=kP/(nkc*volume_r*1d-30) ; kC=kC/(nkc*volume_r*1d-30)
 wP_k=wP_k/(nkc*volume_r*1d-30) ; wC_k=wC_k/(nkc*volume_r*1d-30)
 end subroutine wigner_kappa

!---------------------------------------------------------------
 subroutine write_modes()
!! 2.1 (2026-11-09): mode-resolved and cumulative kappa at this temperature, from the per-IBZ-point tensors
!! (each star summed, so the components are symmetrized), with the totals of kappa_vs_T.dat and
!! kappa_wigner.dat:  RTA = kappa_P per mode (wP_k),  full = full solution + kappa_P - RTA per mode
!! (kappa_k + wP_k - cR_k, the block form of kappa_vs_T.dat),  Wigner = kappa_P + kappa_C (wP_k + wC_k).
!! Mean free path |v| tau in nm: tau = 1/tau_RTA for RTA and Wigner, 1/tau_eff of the full solution for full
!! (mode_kappa: tau_eff = tau_RTA (F_RTA.v)/(F.v)); 0 where the rate is not positive. Cumulative: modes sorted by their mean free path (each kind
!! separately) or by frequency (without the 3 acoustic modes at Gamma, gamma_acoustic), added from the smallest.
 integer u,i,ii,m,n,k,o,lst
 character(len=40) tl
 real(8), allocatable :: om(:),vn(:),mfp(:,:),c(:,:,:),cum(:)
 integer, allocatable :: ord(:,:)
 logical, allocatable :: ac(:)
 real(8) t(3),cs(4,3),wm(4,4)
 n=nibz*ndyn
 allocate(om(n),vn(n),mfp(n,3),c(n,4,3),ord(n,3),cum(4),ac(n))
 m=0
 do i=1,nibz ; do ii=1,ndyn
    m=m+1
    ac(m)=gamma_acoustic(ii,mapinv(i))
    om(m)=frequency(mapinv(i),ii)
    vn(m)=sqrt(sum(veloc(:,ii,mapinv(i))**2))*c_light       ! veloc is in units of c (veltest.dat: x c_light = m/s)
    t(1)=0 ; if (tauinvR(i,ii).gt.0) t(1)=1/tauinvR(i,ii)
    t(2)=0 ; if (tauinv_eff(i,ii).gt.0) t(2)=1/tauinv_eff(i,ii)
    mfp(m,1)=vn(m)*t(1)*1d9 ; mfp(m,2)=vn(m)*t(2)*1d9 ; mfp(m,3)=mfp(m,1)      ! nm
    do k=1,3
       c(m,k+1,1)=wP_k(i,ii,k,k)
       c(m,k+1,2)=kappa_k(i,ii,k,k)+wP_k(i,ii,k,k)-cR_k(i,ii,k,k)
       c(m,k+1,3)=wP_k(i,ii,k,k)+wC_k(i,ii,k,k)
    enddo
    c(m,1,:)=(c(m,2,:)+c(m,3,:)+c(m,4,:))/3
 enddo ; enddo
 if (abs(tempk-nint(tempk)).lt.1d-6) then
    write(tl,'(a,i0)') 'T',nint(tempk)
 else
    write(tl,'(a,f0.1)') 'T',tempk
 endif
 ! per mode
 open(newunit=u,file='kappa_modes_'//trim(tl)//'.dat',status='replace')
 write(u,'(a,f8.2,a)')'# T = ',tempk,' K: kappa (W/mK) of each IBZ point and mode, its star summed (sums = the totals of'// &
 &   ' kappa_summary.dat); avg = (xx+yy+zz)/3'
 write(u,'(a)')'# mean free paths mfp = |v| tau in nanometres (nm): tau_RTA for RTA and Wigner, tau_eff = tau_RTA (F_RTA.v)/(F.v) of the'// &
 &   ' full solution for full (mfp 0 where tau_eff <= 0 or v = 0); kappa_C half to each mode of a pair'
 write(u,'(a)')'# ibz mode omega(1/cm) omega(THz) |v|(m/s) tau_RTA(ps) tau_eff(ps) mfp_RTA(nm) mfp_eff(nm)'// &
 &   '  RTA: avg xx yy zz  full: avg xx yy zz  Wigner: avg xx yy zz'
 m=0
 do i=1,nibz ; do ii=1,ndyn
    m=m+1
    t=0
    if (tauinvR(i,ii).gt.0) t(1)=1d12/tauinvR(i,ii)
    if (tauinv_eff(i,ii).gt.0) t(2)=1d12/tauinv_eff(i,ii)
    write(u,'(i6,i4,7(1x,es11.4),12(1x,es11.4))') i,ii,om(m),om(m)*c_light*100*1d-12,vn(m),t(1),t(2),mfp(m,1),mfp(m,2), &
         c(m,:,1),c(m,:,2),c(m,:,3)
 enddo ; enddo
 close(u)
 ! cumulative versus mean free path, each kind sorted by its own mfp
 do k=1,3
    call sort_real(n,mfp(:,k),ord(:,k))
 enddo
 open(newunit=u,file='kappa_cumul_mfp_'//trim(tl)//'.dat',status='replace')
 write(u,'(a,f8.2,a)')'# T = ',tempk,' K: cumulative kappa (W/mK) versus mean free path (in nanometres, nm), modes sorted by'// &
 &   ' their mean free path, each kind by its own; the last line holds the totals'
 write(u,'(a)')'# mfp_RTA(nm)  RTA: avg xx yy zz   mfp_eff(nm)  full: avg xx yy zz   mfp_RTA(nm)  Wigner: avg xx yy zz'
 cs=0
 do o=1,n
    do k=1,3
       cs(:,k)=cs(:,k)+c(ord(o,k),:,k)
    enddo
    write(u,'(3(es11.4,4(1x,es11.4),2x))') (mfp(ord(o,k),k),cs(:,k),k=1,3)
 enddo
 close(u)
 ! cumulative versus frequency
 call sort_real(n,om,ord(:,1))
 open(newunit=u,file='kappa_cumul_freq_'//trim(tl)//'.dat',status='replace')
 write(u,'(a,f8.2,a)')'# T = ',tempk,' K: cumulative kappa (W/mK), modes sorted by frequency (the acoustic modes at Gamma left out);'// &
 &   ' the last line holds the totals'
 write(u,'(a)')'# omega(1/cm) omega(THz)  RTA: avg xx yy zz   full: avg xx yy zz   Wigner: avg xx yy zz (no acoustic modes at Gamma)'
 cs=0
 do o=1,n
    m=ord(o,1)
    if (ac(m)) cycle                                     ! acoustic modes at Gamma (omega ~ 0, no heat)
    cs=cs+c(m,:,:)
    write(u,'(2(es11.4,1x),3(4(1x,es11.4),1x))') om(m),om(m)*c_light*100*1d-12,cs(:,1),cs(:,2),cs(:,3)
 enddo
 close(u)
 ! Wigner over the lowest bands, and the summary line
 wm=0
 m=0
 do i=1,nibz ; do ii=1,ndyn
    m=m+1
    do lst=1,4
       if (lst.le.3 .and. ii.gt.lst) cycle
       wm(:,lst)=wm(:,lst)+c(m,:,3)
    enddo
 enddo ; enddo
 write(uwm,'(f8.2,16(1x,es11.4))') tempk,wm(:,1),wm(:,2),wm(:,3),wm(:,4)
 flush(uwm)
 do k=1,3
    cum=0
    do m=1,n
       cum=cum+c(m,:,k)
    enddo
    cs(:,k)=cum
 enddo
 write(usum,'(f8.2,12(1x,es11.4))') tempk,cs(:,1),cs(:,2),cs(:,3)
 flush(usum)
 deallocate(om,vn,mfp,c,ord,cum,ac)
 end subroutine write_modes

!---------------------------------------------------------------
 subroutine sort_real(n,x,ord)
!! ord: the permutation sorting x in increasing order (merge sort, stable)
 integer, intent(in) :: n
 real(8), intent(in) :: x(n)
 integer, intent(out) :: ord(n)
 integer, allocatable :: tmp(:)
 integer w,lo,mid,hi,i,j,m
 allocate(tmp(n))
 do i=1,n
    ord(i)=i
 enddo
 w=1
 do while (w.lt.n)
    lo=1
    do while (lo.le.n)
       mid=min(lo+w-1,n) ; hi=min(lo+2*w-1,n)
       i=lo ; j=mid+1 ; m=lo
       do while (i.le.mid .and. j.le.hi)
          if (x(ord(i)).le.x(ord(j))) then
             tmp(m)=ord(i) ; i=i+1
          else
             tmp(m)=ord(j) ; j=j+1
          endif
          m=m+1
       enddo
       do while (i.le.mid)
          tmp(m)=ord(i) ; i=i+1 ; m=m+1
       enddo
       do while (j.le.hi)
          tmp(m)=ord(j) ; j=j+1 ; m=m+1
       enddo
       ord(lo:hi)=tmp(lo:hi)
       lo=lo+2*w
    enddo
    w=2*w
 enddo
 deallocate(tmp)
 end subroutine sort_real

!---------------------------------------------------------------
 subroutine write_F(label)
!! the distribution function FFF in the IBZ at the present temperature, to distribution_F.dat
 character(*), intent(in) :: label
 integer kk
 write(udist,'(a,f9.2,2a)')'# T(K)= ',tempk,'  solution: ',label
 do kk=1,nibz ; do la=1,ndyn
    write(udist,'(i6,i4,3(1x,es11.4))')kk,la,FFF(kk,la,:)
 enddo ; enddo
 write(udist,*)
 end subroutine write_F

!---------------------------------------------------------------
 subroutine gather_rows(at_loc,at_glob)
!! at_glob(n, all rows) on rank 0 from the ranks' at_loc(n, their rows) (transposed blocks)
 real(8), intent(in) :: at_loc(:,:)
 real(8), intent(inout) :: at_glob(:,:)
 integer mine,r
 real(8) dummy(1)
 if (.not.allocated(counts)) allocate(counts(nprocs),displs(nprocs))
 mine=size(at_loc)
 call MPI_ALLGATHER(mine,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
 displs(1)=0
 do r=2,nprocs
    displs(r)=displs(r-1)+counts(r-1)
 enddo
 if (root) then
    call MPI_GATHERV(at_loc,mine,MPI_DOUBLE_PRECISION,at_glob,counts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 else
    call MPI_GATHERV(at_loc,mine,MPI_DOUBLE_PRECISION,dummy,counts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 endif
 end subroutine gather_rows

!---------------------------------------------------------------
 subroutine iterate_bte(niter)
!! job 3: F <- F_RTA + (1/Q) sum P*(F' -+ F''), mixed with update_mix, until kappa changes by
!! less than 1d-4 W/mK (at most max_iter iterations).  Needs P1,P2 of this rank's points:
!! kept from thermal_rates in memory mode, recomputed point by point from V3DELTA otherwise.
 integer, intent(out) :: niter
 integer kk
 real(8), allocatable :: kold(:,:,:)
 if (.not.allocated(F_RTA)) allocate(F_RTA(nibz,ndyn,3))
 if (.not.allocated(F2)) allocate(F2(nibz,ndyn,3),F2loc(nloc,ndyn,3))
 allocate(kold(ndyn,3,3))
 F_RTA=FFF
 do niter=1,max_iter
    kold=kappa
    if (v3_store.eq.1) then
       if (allocated(F2mpi)) deallocate(F2mpi)
       allocate(F2mpi(max(nloc,1),ndyn,3))
       if (nloc.gt.0) call update_mpi(ndyn,kglob)
       F2loc=F2mpi(1:nloc,:,:)
    else
       if (allocated(F2mpi)) deallocate(F2mpi)
       allocate(F2mpi(1,ndyn,3))
       do kk=kglob(1),kglob(2)
          kl=kk-kglob(1)+1
          if (.not.read_v3(kk)) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          ksub=[kk,kk]
          call calculate_tauRTA_mpi2(ksub,ndyn,tempk)   ! P1,P2 of this point at this T
          call update_mpi(ndyn,ksub)
          F2loc(kl,:,:)=F2mpi(1,:,:)
       enddo
    endif
    allocate(buf(nloc*ndyn*3),gbuf(nibz*ndyn*3))
    s=0
    do kl=1,nloc ; do la=1,ndyn ; do i=1,3
       s=s+1
       buf(s)=F2loc(kl,la,i)
    enddo ; enddo ; enddo
    call gather_all(ndyn*3,buf,gbuf)
    s=0
    do kk=1,nibz ; do la=1,ndyn ; do i=1,3
       s=s+1
       F2(kk,la,i)=gbuf(s)
    enddo ; enddo ; enddo
    deallocate(buf,gbuf)
    FFF = FFF + update_mix*(F2-FFF)
    call mode_kappa(ndyn,tempk,veloc,eigenval)
    if (root) then
       write(kappai,*)'================= T, iteration',tempk,niter
       call write_kappa(ndyn,tempk,veloc,eigenval,kappai,tauefi)
    endif
    normk=sqrt(sum((kappa-kold)**2)/(ndyn*9))
    write(ulog,*)' iteration ',niter,' change in kappa ',normk
    if (normk.lt.1d-4) exit
 enddo
 niter=min(niter,max_iter)
 end subroutine iterate_bte

 end program kappa_full
