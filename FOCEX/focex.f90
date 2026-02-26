!=====================================================

program FOCEX
!
! todo:
!
! add option of reading existing fc files and fitting the remaining ones.
! also writeout the fitted FCs in the format of other codes: ALAMODE, PHONOPY & SHENG
!
! Add the Coulomb contribution to FC3 and FC4
! Add reading of the stress tensor from DFT in case some adjustments of the equilibrium positions and
! lattice parameters are needed (requires high cutoff in the DFT calculations to be reliable)
!------------------------------------------------
!! Program to extract force constants from ab initio force-displacement data
!! eventually with imposed (linear) constraints of symmetry
!! INPUT FILES : structure.params, dielectric.params, default.params, POSCARi, OUTCARi (i=1,2,...)
!!               and for lattice dynamical calcs: kpbs.params and latdyn.params
!! in structure.params file user specifies the primitive cell and its atoms
!! the range of the force constants of desired ranks, 
!! whether or not translational, rotational and Huang invariance constraints are imposed, 
!! specifications(type,mass,name) of atoms in the primitive cell and their reduced coordinates 
!! in units of translations of the CONVENTIONAL cell.
!! In POSCARi file contains the super cell structure and atomic positions in VASP POSCAR format
!! Finally, FORCEDISPi contains the atomic displacement and corresponding forces on each atom in the
!! supercelli. It is obtained using a postprocessing utility (readoutcar) of VASP or QE output files.
!! OUTPUT FILES: log...dat contains log of the run and intermediate outputs
!! fcr_irr.dat fcr.dat : contain the obtained irreducible and full set of force constants of rank r
!! lat_fc.dat : coordinates of atoms surrounding the primitive cell within some range (15 Ang by default)
!! and the number of force constants of each rank that resulted from the SVD
!! svd-results.dat: contains the output of SVD algorithm, to check for errors, and condition number
!! and a bunch of other files for inspection of the run but not needed for postprocessing
!! such as amatrx.dat, maps.dat, etc...
!! by K. Esfarjani, March 2023
!------------------------------------------------
!* need to treat cases where equilibrium positions are not known but are
! calculated as outputs of this program. For this purpose, we need to be
! able to adjust the FCs to equilibrium points as they are calculated at the
! points defined in POSCAR : First need to solve F(u)=0 to get u0 and then use:
! phi3 = phi3 + u0 phi4; phi2 = phi2+ u0 phi3 + u0 u0 phi4/2;
! phi1 = phi1 + u0 phi2 + u0 u0 phi3/2 + u0 u0 u0 phi4/6; phi4 unchanged.
! and also apply the iterative Newton's algorithm to find the eq. positions
!* must also give as output the equilibrium volume by calculating the stress and bulkmod.
!
use ios
use lattice
use params
use atoms_force_constants
use ewald
use svd_stuff
use geometry
use fourier
use kpoints 
use born
use linalgb
use tetrahedron
use eigen !, only : ndyn,eigenval,eigenvec,allocate_eig,eigenval_bs,eigenvec_bs,allocate_eig
use constants, only : r15,ee,pi,eps0,cnst

implicit none
integer i,j,g,ti,rank,iunit,uio,nkernel,imax,n_hom,nt(maxrank),ntind(maxrank),tau,taup,nat,ncs(3),al,be,ia,jb
character xt*1,fn*7,now*10,today*8,zone*5,fni*11,bfl,bnp
logical ex,isperiodic
real(r15), allocatable :: xout(:),foldedk(:,:),rand(:),rand2(:),phi_naq0(:,:,:,:),phi_na(:,:,:,:,:)
complex(r15), allocatable :: auxr(:),auxg(:),aravg(:),axavg(:),ax1(:),ax2(:),dynmat(:,:),ddyn(:,:,:)
real(r15), allocatable :: mat(:,:),qmat(:)
integer, allocatable :: frc_constr(:),mp(:) 
real(r15) error,ermax,sd(4),sig, volmax,sf,dsf(3) ,q2(3),qr(3),largestcutoff,matr(3,3)
type(vector):: y1,y2,y3,k1,k2,k3
real tim
real(r15) gmax,ene,grad,pos(3,64) ,ewald_force(3,64),deteps3, etaewa,xx,asr
complex(r15) d2ew(3,3), d3ew(3,3,3), d0(3,3),sum
real(r15) q(3)

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(tim)

 open(utimes,file='times.dat' ,status='unknown')
 open(umap  ,file='maps.dat'  ,status='unknown')
 open(umatrx,file='amatrx.dat',status='unknown')
 open(ucor  ,file='corresp.dat',status='unknown')


 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' Program FOCEX was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(utimes,'(a,f12.4)')' At the start of the Program FOCEX the  TIME  IS ',tim
!
!@@@@@@@@@@@@@@@@@@@@@@  read inputs and set up cell information  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! read from structure.params the atoms in prim cell in reduced units, assign their mass, type, tau;
! also read specs of the prim cell, and flags for range of FC3 and FC4, enforcement of invariance relations
! outputs: lattparams,primlat,nshells,include_fc,itrans,irot,natom_prim_cell,mas,atname,
! atompos0 (reduced coords within the primitive cell) and object atom0(natom_prim_cell)
   call read_structure("log")

! reads dielectric constant and Born charges (=0 if unknown)
   call read_dielectric

   call read_latdyn

   call make_unitcell
! the number of shells is based on nshells in the structure.params file 
! inputs: natom_prim_cell,latticeparameters,primitivelattice,atom_type,atompos0
! outputs: r0i, atompos (cart coord of atoms and natoms neighbor shells),
! iatomcell, iatomcell0, iatomneighbor, iatomop, ... and symmetry operations
! sets nsmax, according to the default range for fc2s, rcut(2)

   coef=ee*1d10/eps0/volume_r0 ! to convert 1/ang^2 to eV/ang , includes cancellation of 4p (NA correction)
   ndyn=3*natom_prim_cell
   deteps3=det(epsil)**(1/3d0)
   a0_scale = (volume_r0**(1/3d0))
   np = 1 ! take up to np^th shell of G's in the NA-Parlinski term 

   write(ulog,7)'np,a0scale,deteps3,coef=',real(np),a0_scale,deteps3 , coef


   call find_WS_largest_SC(imax,volmax)
! of available supercells, takes the one with largest volume to eventually set
! the range of FC2= largest center to WS boundary distance; record rs1,rs2,rs3
! output is rs1,rs2,rs3 of this largest supercell

   call make_grids

   call write_atompos  !  atompos in .xyz format for visualization


 call cpu_time(tim)
 write(utimes,'(a,f12.4)')' make_unitcell, TIME                          IS ',tim

   call set_neighbor_list
! outputs: atom0%equilb_pos,shells%no_of_neighbors , rij, neighbors(:)%tau and n and maxshell
! maxshell is the actual # of shells within rcut(2). It is not nsmax; maxshell < maxneighbors

 write(ulog,*)' After set_neighbors: maxshells=',maxshells,' while maxneighbors=',maxneighbors
 write(ulog,*)' ************************* Writing the neighborlist ************************** '

   call write_neighbors
   etaewa=eta

! if input range is larger than lgridmax=cutoff(SC) or if fc2range=0, set the largest shell corresponding to that cutoff 
   largestcutoff=0
   do i=1,natom_prim_cell
      largestcutoff=max(largestcutoff,atom0(i)%shells(nshells(2,i))%radius) 
   enddo
   write(ulog,*)' Largest cutoff imposed instructure.params file is ',largestcutoff
   if(fc2range.eq.0 .or. largestcutoff.gt.lgridmax) then 
      call update_nshells(nrgrid,rgrid,lgridmax) ! neighborshells should stay within the WS cell of largest supercell < lgridmax
      write(ulog,*)' nshells2 updated according to cutoff = ',lgridmax
      write(ulog,*)' use largest WS; new updated nshells2 = ',nshells(2,1:natom_prim_cell)
   endif


! create a second subgrid of translation vectors in the supercell WS that has 
! full symmetry of the lattice for NA subtraction in case fc2range.ne.0
! output is nsubgrid, subgrid,subgrid_weights
   if(fc2range.ne.0) call make_subrgrid_symmetric(nrgrid,rgrid,rws26,nsubgrid)  


! define the dimensions of the FC arrays, for each rank if include_fc.ne.0
! collects identical force_constants in groups of shells, identifies the irreducible ones and
! finds the mapping between the irreducible ones and all FCs in that group
! inputs: outputs of force_constants_init: iatomop etc...
! outputs: map structure including for rank, the groups of indep FCs and the full FCs
! maxterms, nshells, ngroups, nterms, estimate of inv_constraints, nindepfc

   call setup_maps

   close(umap)


   write(ulog,*)"rank, indep FC terms , total FC terms , # of groups "
   do rank=1,maxrank
      write(ulog,*)rank,map(rank)%ntotind,map(rank)%ntot,map(rank)%ngr
   enddo

  call cpu_time(tim)
  write(utimes,'(a,f12.4)')' TIME after setup_maps and before write_correspondance IS ',tim

   do rank=1,maxrank
      ntind(rank)= map(rank)%ntotind
      nt   (rank)= map(rank)%ntot
   enddo
   call write_lat_fc(ntind,nt,'lat0_fc.dat')

!--------------------------------------------------------------------
! read FCs from a file if already present
! outputs: ngrnk(rank), irreducible fcs:fcrnk(rank,ngrnk(rank))
   do rank=1,maxrank
      if (include_fc(rank) .eq. 2 ) then
         write(xt,'(i1)')rank
         fni= 'fc'//xt//'_irr.dat'
         iunit = ufco+rank
         inquire ( file=fni, exist=ex)
         if (ex) then
            open(iunit,file=fni,status='old')
! this rank is to be read and is not included in svd extraction
            ngrnk(rank)=map(rank)%ntotind; allocate(fcrnk(rank,ngrnk(rank)))
!!!        call read_fcs(iunit,fn,rank,fcrnk(rank,:),ngrnk(rank))
         endif
      endif
   enddo

! if read, may need to check compatibility with given supercells?!?
!--------------------------------------------------------------------

! before imposing the constraints, first define the number of columns in matrix A
   if (.not.readfc) then
      do rank=2,2 !1,maxrank
         allocate(map(rank)%keep(map(rank)%ntotind))  ! which indep fc2s to keep based on rws26
      enddo
   endif

   do i=1,26
      call write_out(ulog,' rws26 ',rws26(:,i))
   enddo

! nkeptind is defined here
   if(fc2range.eq.0) then ! determine the range from the supercell WS 
      call setup_FC2_in_supercell(rws26) 
   else ! use default values found by setup_maps and the required shell numbers; exclude if range is beyond supercell
      map(2)%keep(:)=1  ! keep everything
      call exclude_beyond_sc(map(2)%ntotind,map(2)%keep(:))
      map(2)%nkeptind= sum(map(2)%keep(:)) ! not sum of the groups but sum of all rnk2 indep terms
   endif
! initialize nkeptind for all other ranks
   do rank=1,maxrank
      if(rank.ne.2) map(rank)%nkeptind = map(rank)%ntotind
   enddo

! write the connected pairs of atoms by FC2s ; used for display
   call write_springs

   write(ulog,*)' Old, and new number of fc2:map(2)%nkeptind=',map(2)%ntotind,map(2)%nkeptind

! get the number of independent force constants from setup_maps
   nindepfc=0
   do rank=1,maxrank
      if ( include_fc(rank) .eq. 1 ) then ! exclude them if they already exist and will be read
! for now, we put in amat all fc2 of number ntotind
         nindepfc = nindepfc + map(rank)%nkeptind 
         write(ulog,*)'rank,nindep_rnk,cumulative nindep_kept=',rank,map(rank)%nkeptind,nindepfc
      endif
   enddo
   write(ulog,*)' MAIN: total # of independent FCs of all rank, nindepfc=',nindepfc

!  map(2)%mat set to zero if (g,ti) not kept 
   do g=1,map(2)%ngr
   do ti=1,map(2)%ntind(g)
      write(ulog,'(a,9i4)')'g,ti,counter,keep,current=',g,ti,counteri(2,g,ti),map(2)%keep(counteri(2,g,ti)),current(2,g,ti)
      write(*   ,'(a,9i4)')'g,ti,counter,keep,current=',g,ti,counteri(2,g,ti),map(2)%keep(counteri(2,g,ti)),current(2,g,ti)
      map(2)%gr(g)%mat(:,ti) = map(2)%gr(g)%mat(:,ti)*map(2)%keep(counteri(2,g,ti))
   enddo
   enddo
!
!@@@@@@@@@@@@@@  read supercell force-displacements and setup the matrices @@@@@@@@@@@@@@@@
!

   allocate(frc_constr(fdfiles))  ! for each of the supercell files

 ! Read force-displacements from FORCEDISPi (OUTCARi) and construct aforce,bforce matrices
   call read_snaps_set_aforce(frc_constr,nconfigs)

  call cpu_time(tim)
  write(utimes,'(a,f12.4)')' TIME after read_snaps_set_aforce             IS ',tim

   close(ucor)

!--------------------------------------------------------------------
! depending on the constraints, setup the homogeneous part of the A and B matrices:
! if(born_flag.le.2) then ! use a real space treatment of ASR and rotational invariances

    call estimate_inv_constraints

    call set_translational_inv_constraints
    write(ulog,*)'Number of translational invariance constraints is=',transl_constraints

    call cpu_time(tim)
    write(utimes,'(a,f12.4)')' TIME after set_translational_inv_constraints IS ',tim

    call set_rotational_inv_constraints
    write(ulog,*)'Number of rotational invariance constraints is=',rot_constraints

    call cpu_time(tim)
    write(utimes,'(a,f12.4)')' TIME after set_rotational_inv_constraints    IS ',tim

    call set_huang_inv_constraints
    write(ulog,*)'Number of Huang    invariance constraints is=',huang_constraints

    call cpu_time(tim)
    write(utimes,'(a,f12.4)')' TIME after set_Huang_inv_constraints         IS ',tim

! now that aforce & bforce matrices are setup, we can do SVD and project on
! onto the kernel of homogeneous part of amatrix

! inputs: inv_constraints,force_constraints, aforce, bforce and invariance parts of amat
! allocates amat,bmat and fills them with atrans,arot,ahuang and aforce (same for b)
! call include_constraints_remove_zeros
! outputs are amat and bmat and their dimensions:dim_al,dim_ac, which go to SVD

! output is the homogeneous part of Amatr: ahom(inv_constraints,dim_ac)
    call homogeneous_constraints_overlap(inv_constraints) 

    write(ulog,*)'dim_al,dim_ac set to ',dim_al,dim_ac

! check if any include_fc=2, then read atompos and call setup_maps and read fcs
    if(readfc) then  
       call read_fcs 
! modify ahom,,aforce, bhom,bforce to move known stuff on the RHS
!      call modify_matrices 
    endif 

! give warnings if any column in amat is totally zero
! allocate(amat(inv_constraints+force_constraints,dim_ac),bmat(inv_constraints+force_constraints))
!   call check_zero_column(force_constraints,dim_ac,aforce)


    call cpu_time(tim)
    write(utimes,'(a,f12.4)')' TIME after include_constraints and 0 column  IS ',tim

! write the a and b matrices
   if (verbose) then
      write(umatrx,*)'#========= before call to svd, ahom is:',inv_constraints,dim_ac,'========='
      do i=1,inv_constraints
         write(umatrx,77)(ahom(i,j),j=1,dim_ac) !,bmat(i)
      enddo
      write(umatrx,*)'#========= before call to svd, aforce and bforce are:',dim_al,dim_ac,'========='
      do i=1,force_constraints
         write(umatrx,88)(aforce(i,j),j=1,dim_ac),bforce(i)
      enddo
   endif

77 format(g11.4,299(1x,f9.5))
88 format(299(f10.6))

!
!@@@@@@@@@@@@@@@@@@@@@@  Solve the linear equations using SVD and write the solution  @@@@@@@@@@@@@@@@@@@@@@@@@@
!
! do the SVD decomposition; beware amat is overwritten
   if(allocated(fcs)) deallocate(fcs)
   if (allocated(sigma)) deallocate(sigma)
   allocate(fcs(dim_ac),sigma(dim_ac))

   uio = 349
   if(enforce_inv.eq.1 .or. itemp.eq.1) then ! elimination will be used on ahom

      open(uio,file='elimination.dat',status='unknown')
! separate homogeneous part of amat from inhomogeneous part by finding the first non-zero elt of bmat
!      call find_first_nonzero(dim_al,bmat,n_hom)
      n_hom=inv_constraints
      write(ulog,*)'MAIN: Using kernel projection: enforce_inv and itemp are=',enforce_inv,itemp
      write(ulog,*)'MAIN: position of the first non-zero elt of bmat=',n_hom
      write(   *,*)'MAIN: position of the first non-zero elt of bmat=',n_hom
      allocate(kernelbasis(dim_ac,dim_ac))
!     if(allocated(a_hom)) deallocate(a_hom)
!     allocate(a_hom(n_hom,dim_ac),kernelbasis(dim_ac,dim_ac))
!     a_hom=amat(1:n_hom,1:dim_ac)
!     call write_out(uio,'#####     HOMOGENEOUS PART OF A ',amat(1:n_hom,1:dim_ac))
!     call get_kernel(n_hom,dim_ac,amat(1:n_hom,1:dim_ac),svdcut,nkernel,kernelbasis,uio)
      call get_kernel(n_hom,dim_ac,ahom,svdcut,nkernel,kernelbasis,uio)
      write(ulog,*)' Kernel of A, of dim ',n_hom,'x',dim_ac,' is of dimension ',nkernel
      write(ulog,*)' The number of eliminated components is ',dim_ac-nkernel
!     deallocate(a_hom)

   endif

   allocate(xout(dim_ac))  ! dummy for xout
   if (itemp.eq.1) then

      write(*,*) ' Temperature of ',tempk,' will be implemented with weighting by Boltzmann factor'
      if(allocated(qmat)) deallocate(qmat)
      if(allocated( mat)) deallocate( mat)
      allocate(mat(dim_ac,dim_ac),qmat(dim_ac))
      write(*,*) "Energy values from the OUTCAR1 are: ",energies
      write(*,*) "Dimension of energies are: ",shape(energies)
      write(*,*) "The value of nlines are: ",nlines
      write(*,*) "The value of nconfigs is: ",nconfigs
!     call implement_temperature(dim_al-n_hom,dim_ac,amat(n_hom+1:dim_al,1:dim_ac), &
      call implement_temperature(force_constraints,dim_ac,aforce, &
&            bforce,nconfigs,energies,nlines,tempk,mat,qmat)
      call solve_svd(dim_ac,dim_ac,mat,qmat,svdcut,xout,sigma,uio)
      call write_out(ulog,'TEMP: SVD solution before kernel projection',xout)
      call project_on(dim_ac,nkernel,xout,kernelbasis(1:dim_ac,1:nkernel),fcs)
      deallocate(kernelbasis,mat,qmat)

   else

      if(enforce_inv.eq.1) then 
! svd for force-displacement, followed by projection on the kernel of a_hom to satisfy invariances exactly
         write(ulog,*)'MAIN: Using svd of the force amatrix since enforce_inv and itemp are=',enforce_inv,itemp
         write(ulog,*)'MAIN: size of the homogeneous part is ',n_hom
!         call solve_svd(dim_al-n_hom,dim_ac,amat(n_hom+1:dim_al,1:dim_ac), &
!&             bmat(n_hom+1:dim_al),svdcut,bfrc,sigma,uio)
!        call svd_set(dim_al-n_hom,dim_ac,amat(n_hom+1:dim_al,1:dim_ac),  &
!&            bmat(n_hom+1:dim_al),xout,sigma,svdcut,error,ermax,sig,'svd-elim.dat')
         call svd_set(force_constraints,dim_ac,aforce,  &
 &            bforce,xout,sigma,svdcut,error,ermax,sig,'svd-elim.dat')
         write(ulog,'(a,f9.3,a)')'Percent error, || F_dft-F_fit || / || F_dft || =',sig,' %'
         call write_out(ulog,'Enforce=1: SVD solution before kernel projection',xout)
         write(ulog,*)'MAIN: Invariance violations '
         do i=1,n_hom
            write(ulog,*)i,dot_product(ahom(i,:),xout(:))
         enddo

         write(ulog,*)'MAIN: performing projection on the kernel basis'
         call project_on(dim_ac,nkernel,xout,kernelbasis(1:dim_ac,1:nkernel),fcs)
         call write_out(ulog,'Enforce=1: Solution After kernel projection',fcs)
         write(ulog,7)'Max and Average Change After kernel projection',maxval(abs(fcs-xout)), &
&                      length(fcs-xout)/dim_ac
         deallocate(kernelbasis,aforce,bforce,xout)

      else

         if (allocated(amat)) deallocate (amat)
         if (allocated(bmat)) deallocate (bmat)
         dim_al=inv_constraints+force_constraints
         allocate(amat(dim_al,dim_ac),bmat(dim_al))
         amat(1:inv_constraints,:) = ahom ; bmat(1:inv_constraints) = 0
         amat(1+inv_constraints:dim_al,:) = aforce ; bmat(1+inv_constraints:dim_al) = bforce
         deallocate (aforce,bforce)

         call svd_set(dim_al,dim_ac,amat,bmat,fcs,sigma,svdcut,error,ermax,sig,'svd-all.dat')
         write(ulog,'(a,f9.3,a)')'Percent error, || F_dft-F_fit || / || F_dft || =',sig,' %'
         call write_invariance_violations(ulog,dim_ac,fcs)
         deallocate(amat,bmat) 

      endif

   endif
   close(uio)
   close(umatrx)


   do rank=1,4
      write(xt,'(i1)')rank
      if(include_fc(rank).eq.1) then
         open(ufit1-1+rank ,file='fc'//xt//'_irr.dat',status='unknown')
         open( ufc1-1+rank ,file='fc'//xt//'.dat',status='unknown')
      endif
   enddo

   write(ulog,*)'map(ntotind)=',map(:)%ntotind

   call write_independent_fcs(dim_ac,sigma,sd,ulog)

   call write_output_fcs

   call create_extended_rgrid  ! to include all grid points covered by FC2s beyond supercell WS

   if (born_flag.ge.1) then

! for non-analytical terms: reset the ewald grids 
     call cuts(deteps3,a0_scale,born_flag,eta,rcutoff_ewa,gcutoff_ewa)
     call make_grid_shell_ewald(r01,r02,r03,rcutoff_ewa,gcutoff_ewa,epsil)
     write(6   ,4)'MAIN before LD: eta,rcut_ewa,gcut_ewa=',eta,rcutoff_ewa,gcutoff_ewa
     write(6   ,*)'nr_ewa,ng_ewa=',nr_ewald,ng_ewald
     write(ulog,4)'MAIN before LD: eta,rcut_ewa,gcut_ewa=',eta,rcutoff_ewa,gcutoff_ewa
     write(ulog,*)'nr_ewa,ng_ewa=',nr_ewald,ng_ewald

     if (.not. allocated(dyn_naq0)) then
        allocate(dyn_naq0(natom_prim_cell,3,3))
        call non_analq0(dyn_naq0,born_flag,1)
     endif
     do tau=1,natom_prim_cell
        call write_out(ulog,' dyn_naq0=',dyn_naq0(tau,:,:))
     enddo

!    allocate( phi_na(natom_prim_cell,natom_prim_cell,3,3,nrgrid_xtnd))
!    call phi_na_g2r(phi_na,nrgrid_xtnd,rgrid_xtnd,nggrid_xtnd,ggrid_xtnd,gws_weights_xtnd)
!    allocate( phi_na(natom_prim_cell,natom_prim_cell,3,3,nrgrid))
!    call phi_na_g2r(phi_na,nrgrid,rgrid,nggrid,ggrid,gws_weights)  ! this is wrong should not be used!
!    open(876,file='phi_na.dat')
!!   do i=1,nrgrid_xtnd
!    do i=1,nrgrid
!    do tau =1,natom_prim_cell
!    do taup=1,natom_prim_cell
!       write(876,33)'Rgrid, tau, taup =',i,tau,taup,rgrid(:,i),cart2red(rgrid(:,i),'r')
!       call write_out(876,'phi_na ',phi_na(tau,taup,:,:,i))
!    enddo 
!    enddo 
!    enddo 
!    close(876)

   endif

! represent the extracted (short-ranged) harmonic FCs on a 5D array on an extended grid
   call get_phiSR_5D

  call cpu_time(tim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After SVD and writing of FCs,          TIME  IS ',tim

   close(ufit1)
   close(ufit2)
   close(ufit3)
   close(ufit4)
   close(ufc1)
   close(ufc2)
   close(ufc3)
   close(ufc4)

   deallocate(sigma) 

!  call energy_strain(0.01,ene,grad) 
!  write(ulog,*)'for strain =0.01, ene,grad=',ene,grad
!
!@@@@@@@@@@@@@@@@@  Calculate phonons and thermodynamic properties  @@@@@@@@@@@@@@@
!
   write(*,*)'Entering set_band_structure '

   call set_band_structure 

  call cpu_time(tim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After phonon band structure,            TIME IS ',tim


  call date_and_time(date=today,time=now,zone=zone)
  call cpu_time(tim)
  write(utimes,'(a,f12.4)')' ENDING TIME OF THE PROGRAM                   IS ',tim
  write(ulog  ,'(a,f10.4)')' ENDING TIME OF THE PROGRAM                   IS ',tim


2 format(a,f7.3,a)
3 format(i5,2x,99(1x,g12.5))
4 format(a,3(1x,f9.4),1x,99(g12.4))
5 format(9(1x,f11.6))
6 format(i5,99(1x,g11.4))
7 format(a,2(2x,g12.5),2x,9(1x,g11.4))
8 format(i4,3(2x,f10.5),3x,i2,2x,'(',3(1x,i2),1x,')',2x,i4)
9 format(i5,9(1x,f9.4))
33 format(a,3(1x,i4),1x,99(g12.4))
99 format(1000(1x,g10.3))

 write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' Program  FOCEX ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' Program FOCEX ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

  write(ulog,7)' MAE and largest errors in force=',error,ermax
  write(ulog,7)' sd(rank) errors in FOCEX=',sd

  close(ulog)
  close(ufco)
  close(utimes)

  contains

!-------------------------------------------------------

 subroutine read_snaps_set_aforce(frc_cnstr,nconfigs)
!! read all POSCARi/OUTCARi files and the force displacements as inputs to matrices aforce,bforce
!! outputs: energies, aforce, bforce (and full displ & force arrays)
 use params, only : fdfiles
 use ios   , only : ulog,utimes
 use lattice
 use geometry
 use fourier
 use svd_stuff, only : energies,aforce,bforce,nindepfc,force_constraints
 use born, only : born_flag
 use atoms_force_constants, only : natom_super_cell,displ,force,energy,atom_sc
 use linalgb
 implicit none
 integer, intent(out) :: frc_cnstr(fdfiles),nconfigs
! real(r15), intent(out) :: energies(nconfigs),disp(3,natom_super_cell,nconfigs),force(3,natom_super_cell,nconfigs)
! real(r15), allocatable, intent(out) :: energies(:) !,disp(:,:,:),force(:,:,:)
 real(r15), allocatable :: afrc(:,:),bfrc(:) !,aux3(:,:,:),aux4(:,:,:)
 integer i,ncfg,j,nfctot 
 integer, allocatable :: nlin2(:) 
 character xt*1,poscar*7, outcar*10
 real tim

     allocate(nlin2(50000))  ! this is max # of snapshots
! this loop is to calculate # of force_constraints in order to allocate aforce
     nconfigs=0 ; nfctot=0
     structures: do i=1,fdfiles

        write(xt,'(i1)')i
        write(ulog,*)'=================== Reading supercell #',i,'========================'
        poscar='POSCAR'//xt
        call read_supercell(poscar) ! equilibrium atomic positions in supercell from POSCARi

        call make_r0g

! make sure of consistency POSCAR with input structure
        call identify_atoms_in_supercell

        call write_supercell(i)

! writing the coordinates of primitive cells in the supercell
! reference of each prim cell (n) is atom #1
        open(111,file='Bravaissupercell-'//xt//'.xyz')
        write(111,*) natom_super_cell/natom_prim_cell+4  ! also draw the 3 translation vectors
        write(111,*) 'super_cell ' ,i
        write(111,6)'Cs 0 0 0  '
        write(111,6)'Cs ',rs1
        write(111,6)'Cs ',rs2
        write(111,6)'Cs ',rs3
        do j=1, natom_super_cell
           if(atom_sc(j)%cell%tau .ne. 1) cycle
           write(111,6)'Si ',atom_sc(j)%equilibrium_pos  ,atom_sc(j)%cell%n ,atom_sc(j)%cell%tau
        enddo
        close(111)

! was OUTCAR before and changed to FORCEDISP
        outcar='FORCEDISP'//xt

        call count_configs(outcar,ncfg)
        nconfigs = nconfigs + ncfg
        write(ulog,*)'fd file#, Cumulative # of configs=',i,nconfigs

        frc_cnstr(i)=3*natom_super_cell*ncfg
        call allocate_edf(natom_super_cell,ncfg)
        write(*,*)' Allocation of energy, dsp and frc arrays DONE!'
        call read_force_position_data(outcar,ncfg,energy,displ,force,nlin2)  ! displ does not include equilibrium positions

        if(frc_cnstr(i) .ne. ncfg*natom_super_cell*3) then
           write(ulog,*)'READ_SNAPS: inconsistency in frc_cnstr,3n_sc*ncfg=',frc_cnstr(i),ncfg*natom_super_cell*3
           stop
        endif

  call cpu_time(tim)
  write(utimes,'(a,i3,a,f10.4)')' STRUCTURE ',i,' after read_force_position_data time IS ',tim

! 3 ways: if=0 ignore, if=1 just all if=2 use G-space ewald, if=3 or 4  use Fourier transforms
 !     write(ulog,*) 'last Forces BEFORE subtraction'
 !     write(*,*) 'last Forces BEFORE subtraction'
 !     call write_out(ulog,'Last Force',transpose(force(:,:,ncfg)))
 !     if(verbose) call write_forces_display(natom_super_cell,ncfg,displ,force,'beforesub')
 !

  if(born_flag .lt. 10)  call subtract_coulomb_force(born_flag,ncfg,displ,force)
! now frc contains the short-range (Non-Coulomb) part of the forces which we can fit
 !
 !      write(ulog,*) 'last Forces AFTER subtraction'
 !      write(*,*) 'last Forces AFTER subtraction'
 !      call write_out(ulog,'Last Force',transpose(force(:,:,ncfg)))
   ! if(verbose) call write_forces_display(natom_super_cell,ncfg,displ,force,'aftersub')

  call cpu_time(tim)
  write(utimes,'(a,i3,a,f10.4)')' STRUCTURE ',i,' after subtract_coulomb time IS ',tim

! now append new snapshot data to existing one
        if(i.eq.1) then
!       if(fdfiles.eq.1) then

            allocate(energies(ncfg))
            energies=energy
            allocate( aforce(frc_cnstr(i),nindepfc),bforce(frc_cnstr(i)) )
! put displ and force arrays (defined in atoms_force_constants) into aforce and bforce matrices
            call set_force_displacement_matrix(ncfg,frc_cnstr(i)  ,aforce,bforce)
            if(nconfigs.ne.ncfg) print*,'APPENDING ERROR, ifile=',i

        else ! append new arrays for i.ge.2 to the end of existing ones: energies,aforce,bforce

            allocate(afrc(frc_cnstr(i),nindepfc),bfrc(frc_cnstr(i)))
            write(ulog,*)' calling set_force_displacement_matrix, for FORCEDISP#',i
            write(ulog,*)' number of lines in afrc and bfrc=',frc_cnstr(i)
            call set_force_displacement_matrix(ncfg,frc_cnstr(i),afrc,bfrc)

!            energies=reshape(energies,shape=(/size(energies)+size(energy)/),pad=energy)
!            bforce  =reshape(bforce,shape=(/size(bforce)+3*natom_super_cell/),pad=bfrc)
!            aforce=reshape(aforce,  &
! &                 shape=(/size(aforce(:,1))+frc_cnstr(i),nindepfc/),  &
! &                 pad=afrc,order=(/2,1/))  ! in here, there is a pb with order
! or alternatively can use:

       call write_out(ulog,'energies before ',energies)
       call write_out(ulog,'energy   before ',energy)
!            call append_array(energies,energy,energies)
energies=reshape(energies,shape=(/size(energies)+size(energy)/),pad=energy)
       call write_out(ulog,'energies after ',energies)
       call write_out(ulog,'nindepfc  ',nindepfc)

     if(verbose)  call write_out(ulog,'bforce before ',bforce(1:5))
     if(verbose)  call write_out(ulog,'bfrc   before ',bfrc(1:5))
!          call append_array(bforce,bfrc,bforce)
 bforce=reshape(bforce,shape=(/size(bforce)+size(bfrc)/),pad=bfrc)
     if(verbose)  call write_out(ulog,'bforce after ',bforce(1:5))
     if(verbose)  call write_out(ulog,'bforce after ',bforce(1+frc_cnstr(i):5+frc_cnstr(i)))

     if(verbose)  call write_out(ulog,'aforce before ',aforce(1:5,1:5))
     if(verbose)  call write_out(ulog,'afrc   before ',afrc(1:5,1:5))
!           call append_array(aforce,afrc,aforce)
            aforce=reshape(transpose(aforce),shape=(/size(aforce,2), &
&           size(aforce,1)+size(afrc,1)/),pad=transpose(afrc),order=(/1,2/))
            aforce=transpose(aforce)
     if(verbose)  call write_out(ulog,'aforce after ',aforce(1:5,1:5))
     if(verbose)  call write_out(ulog,'aforce after ',aforce(1+frc_cnstr(i):5+frc_cnstr(i),1:5))

            deallocate(afrc,bfrc)

        endif


        deallocate( energy,displ,force)

        write(ulog,*)' set_force_displacement_matrix called, aforce being filled now'
        write(ulog,*)' nfctot,nfctot+frc_constr(i)=',nfctot,nfctot+frc_cnstr(i)
        nfctot = nfctot + frc_cnstr(i)

     enddo structures

! if(allocated(nlines)) deallocate(nlines)
 allocate(nlines(nconfigs))
 nlines=nlin2(1:nconfigs)
 deallocate(nlin2)

     write(ulog,*)'READ_SNAPS: total number of FD constraints(=dim_al)',nfctot

     force_constraints=sum(frc_cnstr)

     if (nfctot.ne.force_constraints) then
        write(ulog,*)'nfctot, force_constraints are not equal ',nfctot,force_constraints
        write(ulog,*)'the program will stop'
        write(*   ,*)'nfctot, force_constraints are not equal ',nfctot,force_constraints
        write(*   ,*)'the program will stop'
        stop
     endif

! extract dynamical matrix from Fourier transforms (uk,fk) using SVD
! the following extract_Fourier requires the knowledge of the IBZ points

!     call extract_fourier(nconfigs,natom_super_cell,displ,force,  &
!  &           nrgrid,rgrid,nggrid,ggrid,map_rtau_sc)

5 format(9(1x,f12.6))
6 format(a,3(1x,f12.6),2x,'(',3i3,')',i4)
9 format(i5,9(1x,f12.6))

 end subroutine read_snaps_set_aforce

!-------------------------------------------------------

 subroutine get_phiSR_5D
!! takes the dyn(G*) calculated at supercell reciprocal vectors G*; works for born_flag > 10
!! subtracts dyn^NA(G*), fourier back to update phi_periodic used later for fourier interpolation of the dynamical matrix
!! works correctly if fc2range=0 full SC is sampled; otherwise DYN_NA should be calculated 
!! on a subgrid conistent with nshell that has full symmetry of the crystal
 use ios
 use atoms_force_constants
 use params
 use lattice
 use constants
 use fourier
 use born
 use geometry, only : length
 implicit none
 integer tau,taup,al,be,g,t,ti,ired,j,k,l,igrid,nr(3),nboundary,ia,ib,ndyn,uioi,ik,ir
 complex(r15), allocatable :: auxr(:) ,aravg(:),auxg(:) ,agavg(:)
 real(r15) rr(3),rfold(3),dta(3),mi,mj,fc2,asr,q0(3),imag
 logical isperiodic
 complex(r15) dynmat(3*natom_prim_cell,3*natom_prim_cell) ,ddyn(3*natom_prim_cell,3*natom_prim_cell,3)
 complex(r15) dddn(natom_prim_cell,natom_prim_cell,3,3,3) 
 complex(r15) dn0(3,3),ddn0(3,3,3), dyn2(natom_prim_cell,natom_prim_cell,3,3,nggrid),test_dyn,er

 ndyn=3*natom_prim_cell

! initialize phi_bare and phi_sr with nrgrid_xtnd, and dyn_g with nggrid
  call allocate_fc_dyn(natom_prim_cell,nrgrid_xtnd,nggrid_xtnd)
! call allocate_fc_dyn(natom_prim_cell,nrgrid,nggrid)

 write(ulog,*)'GET_PHISR_5D: just allocated arrays of sizes natom0,nr,ng', &
&    natom_prim_cell,nrgrid_xtnd,nggrid_xtnd
 
! use bare fc2s from fitting to represent them on a 5D array called phi_bare
! it needs an extended rgrid, and the 5D structure is needed for Fourier transformations if needed.
! the supercell ggrid; it should be G-periodic with the old phase (step) convention; 
  call get_phi_bare_5D(natom_prim_cell,nrgrid_xtnd,rgrid_xtnd,phi_bare) 
! call get_phi_bare_5D(natom_prim_cell,nrgrid,rgrid,phi_bare) 

! BEWARE phi_bare is complex!!
  call write_fc2matrix_5d(natom_prim_cell,nrgrid_xtnd,rgrid_xtnd,real(phi_bare), &
&                        'trace_phi_lr.dat')

 fc2=maxval(abs(aimag(phi_bare)))
 write(*,*)" Largest imaginary part of phi_bare(should be ZERO!)=", fc2 
 if(fc2.gt.1d-6) stop


 if(born_flag.lt.10) then 
! no subtraction here; NA has already been subtracted from forces in subtract_coulomb_force 

    phi_sr = phi_bare  ! need to be initialized here as it's needed in set_dynamical_matrix_hard
    return

 elseif(born_flag.ge.10) then ! subtract phi_NA(R) from phi_bare(R)

!******** Now subtract dyn_na from dyn_g and fourier invert to get phi_sr ********

! need this call to calculate dyn_na(G*) in case born_flag>10
! on the extended ggrid it can be combined with LR fitted dynamical matrix to obtain the SR part
! but on the original ggrid, it would be commensurate with the supercell
  call dyn_na5D_on_grid(nggrid_xtnd,ggrid_xtnd)  ! this gives the 5D dyn_na, which satisfies ASR 
! call dyn_na5D_on_grid(nggrid,ggrid)  ! this gives the 5D dyn_na, which satisfies ASR 

  open(689,file='phina.dat')

! nrgrid_xtnd=nrgrid
! nggrid_xtnd=nggrid
  allocate(auxr(nrgrid_xtnd),aravg(nrgrid_xtnd))  
  allocate(auxg(nggrid_xtnd),agavg(nggrid_xtnd))  

! Inverse fourier dyn_na and subtract from phi_bare to get phi_sr
  do tau =1,natom_prim_cell
  do taup=1,natom_prim_cell
     do al=1,3
     do be=1,3
        auxg = dyn_na(tau,taup,al,be,:)
        call check_periodic(nggrid,ggrid,g0ws26,auxg,agavg,isperiodic)
 !      call check_periodic(nggrid_xtnd,ggrid_xtnd,g0ws26,auxg,agavg,isperiodic)

        if(.not.isperiodic) then
           call warn4(ulog,'dyn_na was not periodic ',real(dyn_na(tau,taup,al,be,:)))
           write(*,*)" Below are real and imaginary parts of the periodic version of dyn_na "
           write(*,*)"Real(dyn_na)"
           write(*,2)real(dyn_na(tau,taup,al,be,:))
           write(*,*)"Imag(dyn_na)"
           write(*,2)aimag(dyn_na(tau,taup,al,be,:))
           fc2=maxval(abs(dyn_na(tau,taup,al,be,:)-agavg(:)))
           write(ulog,6)"tau,taup,al,be,maxval(dyn_na-dyng_na_avg)=",tau,taup,al,be,fc2
           stop
        endif

! compare this fourier transform to phi_na
        call fourier_k2rxtnd(auxg,auxr,nggrid_xtnd,ggrid_xtnd,gws_weights_xtnd,nrgrid_xtnd,rgrid_xtnd)
   !    call fourier_k2r(auxg,auxr)
      
        do j=1,nrgrid_xtnd
           phi_sr(tau,taup,al,be,j) = phi_bare(tau,taup,al,be,j) - auxr(j)  ! do I need weights?
!          phi_sr(tau,taup,al,be,j) = phi_bare(tau,taup,al,be,j) - phi_na(tau,taup,al,be,j) 
            write(689,8)' ',tau,taup,al,be,j,phi_na(tau,taup,al,be,j) ,auxr(j)
        enddo

     enddo
     enddo
  enddo
  enddo

  close(689)

 endif

8 format(a,4i2,i4,3x,g10.3,2(1x,g10.3))

! get dyn_g from fcs and subtract dyn_na before FFT to compare to phi_sr
 do g=1,nggrid_xtnd
    call fourier_fc2fit_4D(ggrid_xtnd(:,g),dyn_g(:,:,:,:,g),dddn)  ! uses hard phase, needed for FFT
 !  call fourier_fc2fit_4D(ggrid(:,g),dyn_g(:,:,:,:,g),dddn)  ! uses hard phase, needed for FFT
 enddo
 dyn_g=dyn_g - dyn_na ! this is the dynamical matrix of the short-ranged forces on the extended grid

g = 5  ! as a test for the 5th reciprocal lattice vector
do tau=1,natom_prim_cell
do taup=1,natom_prim_cell
   do al=1,3
   do be=1,3
      ! Forward: already have dyn_g(tau,taup,al,be,g) from fourier_fc2fit_4D
      ! Backward: already have phi_sr from fourier_k2r
      ! Forward again at same grid point

      test_dyn=cmplx(0.0_r15,0.0_r15)
      do ir=1,nrgrid_xtnd
         test_dyn=test_dyn+cdexp( ci*(ggrid_xtnd(:,g).dot.rgrid_xtnd(:,ir))) * &
   &                   phi_sr(tau,taup,al,be,ir) * rws_weights_xtnd(ir) ! these weights are=1
!        test_dyn=test_dyn+cdexp( ci*(ggrid(:,g).dot.rgrid(:,ir))) * &
!  &                   phi_sr(tau,taup,al,be,ir) * rws_weights(ir) ! these weights are=1
      enddo
      
      er= dyn_g(tau,taup,al,be,g) - test_dyn

      if(abs(er).gt.1d-6) then
        write(*,*) 'Roundtrip test of phi_sr FAILED: tau,taup,al,be:', tau,taup,al,be,er
        write(ulog,*) 'Roundtrip test FAILED tau,taup,al,be:', tau,taup,al,be
        write(ulog,7) '  Original dyn_g, reconstructed, diff  :', dyn_g(tau,taup,al,be,g), test_dyn, er
      endif 

   enddo
   enddo
enddo
enddo

  deallocate(auxr,aravg)  
  deallocate(auxg,agavg)  

 imag=maxval(abs(aimag(phi_sr))) 
 write(ulog,*)'Largest imaginary part in phi_sr is=',imag
! BEWARE phi_sr is complex!!
  call write_fc2matrix_5d(natom_prim_cell,nrgrid_xtnd,rgrid_xtnd,real(phi_sr),'trace_phi_sr.dat')
! call write_fc2matrix_5d(natom_prim_cell,nrgrid,rgrid,real(phi_sr),'trace_phi_sr.dat')


2 format(199(f9.4))
3 format(3i5,199(f9.4))
4 format(4i5,199(f9.4))
5 format(a,3i5,99(f8.4))
6 format(a,4i5,99(f8.4))
7 format(a,99(2x,2(1x,g10.3)))

 end subroutine get_phiSR_5D

!-------------------------------------------------------
! Loop over fitted force constants
! do g = 1, map(2)%ngr
! do ti = 1, map(2)%ntind(g)
!    if (map(2)%keep(counteri(2,g,ti)) == 0) cycle
!    
!    do t = 1, map(2)%nt(g)
!       tau = map(2)%gr(g)%iat(1,t)
!       al = map(2)%gr(g)%ixyz(1,t)
!       j = map(2)%gr(g)%iat(2,t)
!       taup = iatomcell0(j)
!       be = map(2)%gr(g)%ixyz(2,t)
!       
!       ! This R is where the fitted FC is defined
!       rr = atompos(:,j) - atompos(:,taup)
!       
!       ! Compute phi_na at this specific R
!       phi_na_here = 0
!       do ig = 1, nggrid
!          call naq(ggrid(:,ig), tau, taup, 1e6, dna, ddna)
!          phi_na_here = phi_na_here + dna(al,be) * exp(-ci*dot_product(ggrid(:,ig), rr)) * gws_weights(ig)
!       enddo
!       
!       ! Subtract from fitted FC
!       fcs_SR(ired) = fcs(ired) - real(phi_na_here)
!    enddo
! enddo

!-------------------------------------------------------

 subroutine make_grids
!! makes grid of translation vectors within the Wigner-Seitz cell of the supercell, with smaller weights 
!! if they end up on the WS cell boundary, also makes a grid of supercell reciprocal vectors in the WS
!! of the unitcell, i.e. the First Brilloiun zone. 
!! rws26(output) are the 26 shortest vectors used to define the WS of the largest supercell
!! results are stored in rgrid(3,nrgrid) or ggrid(3,nggrid) along with their weights rws_weights & gws_weights
! rgrid(3,nrgrid) is the grid of translations from r0i inside the WS cell of rsi +boundaries

 use lattice
 use fourier
 use constants
 implicit none
 integer ngrd,cnt
 real(r15), allocatable :: grd(:,:),wei(:)
 integer ngrid
 type(vector) x01,x02,x03,x1,x2,x3  ! local; shortest primitive vectors 

 write(ulog,*)' ENTERING MAKE_GRIDS'
 nrgrid=10000
 allocate(grd(3,nrgrid),wei(nrgrid))
 matr=cart_to_prim ! to get reduced coordinates
 call get_26shortest_shell(r01,r02,r03,r0ws26,x01,x02,x03)
 call get_26shortest_shell(rs1,rs2,rs3,rws26,x1,x2,x3)
! call make_grid_weights_WS(r01,r02,r03,rs1,rs2,rs3,matr,nrgrid,grd,wei,rws26) 
 call make_grid_weights_WS(x01,x02,x03,x1,x2,x3,matr,nrgrid,grd,wei,rws26) 
 allocate(rgrid(3,nrgrid),rws_weights(nrgrid))
 rgrid=grd(:,1:nrgrid) ; rws_weights=wei(1:nrgrid)  ! reassigns in the newer versions of fortran
 deallocate(grd,wei)
 write(ulog,77)rws_weights
 write(ulog,*)' SUM of the RWS_WEIGHTS = ',sum(rws_weights)
 call write_lattice(nrgrid,rgrid,'r_supercell.xyz') ! grid of primitive translations in the WS supercell

     open(98,file='rgrid_raw.xyz')
     open(99,file='rgridWS.xyz')
     call show_ws_boundary(v2a(r01),v2a(r02),v2a(r03),r0ws26,19,'WSR0_boundary.xyz',lgridmax) 
     call show_ws_boundary(v2a(rs1),v2a(rs2),v2a(rs3),rws26 ,19,'WSR_boundary.xyz' ,lgridmax) 
 
     write(98,*)nrgrid
     write(98,28)"# name ,grid(cnt),weig(cnt),grid_red(cnt),cnt" 
     write(99,*)nrgrid
     write(99,*)"# name, cartesian grid, reduced grid "
     do cnt=1,nrgrid
        write(98,28)"Si ",rgrid(:,cnt),rws_weights(cnt),matmul(matr,rgrid(:,cnt)),length(rgrid(:,cnt)),cnt !,save_boundary(cnt)
        write(99,27)"Si ",rgrid(:,cnt),matmul(matr,rgrid(:,cnt)),length(rgrid(:,cnt))
     enddo
     close(98)
     close(99)
     write(ulog,*)'lgridmax=',lgridmax

27 format(a,99(1x,f10.4))
28 format(a,8(1x,f10.4),3i5)

 nggrid=10000
 allocate(grd(3,nggrid),wei(nggrid))
 matr=transpose(prim_to_cart)/(2*pi)
 call get_26shortest_shell(g01,g02,g03,g0ws26,x01,x02,x03)
 call get_26shortest_shell(gs1,gs2,gs3,gws26 ,x1,x2,x3)
! call make_grid_weights_WS(gs1,gs2,gs3,g01,g02,g03,matr,nggrid,grd,wei,'g',gws26,g0ws26) 
 call make_grid_weights_WS(x1,x2,x3,x01,x02,x03,matr,nggrid,grd,wei,g0ws26) 
 allocate(ggrid(3,nggrid),gws_weights(nggrid))
 ggrid=grd(:,1:nggrid) ; gws_weights=wei(1:nggrid)  ! reassigns in the newer versions of fortran
 deallocate(grd,wei)
 gws_weights = gws_weights * (volume_r0/volume_r) ! introduce 1/N since used for Fourier transforms
 write(ulog,77)gws_weights
 write(ulog,*)' SUM of the GWS_WEIGHTS = ',sum(gws_weights)
 call write_lattice(nggrid,ggrid,'g_supercell.xyz') ! grid of primitive translations in the WS supercell

     open(98,file='ggrid_raw.xyz')
     open(99,file='ggridWS.xyz')
     call show_ws_boundary(v2a(gs1),v2a(gs2),v2a(gs3), gws26,18, 'WSG_boundary.xyz',gmax) 
     call show_ws_boundary(v2a(g01),v2a(g02),v2a(g03),g0ws26,30,'WSG0_boundary.xyz',gmax) 
 
     write(98,28)"# name ,grid(cnt),weig(cnt),grid_red(cnt),cnt" 
     write(99,*)nggrid+26+8
     write(99,*)"# name, cartesian grid, reduced grid "
     do cnt=1,nggrid
        write(98,28)"Si ",ggrid(:,cnt),gws_weights(cnt),matmul(matr,ggrid(:,cnt)),length(ggrid(:,cnt)),cnt !,save_boundary(cnt)
        write(99,27)"Si ",ggrid(:,cnt),matmul(matr,ggrid(:,cnt)),length(ggrid(:,cnt))
     enddo
     do cnt=1,26
        write(99,27)"Bi ",g0ws26(:,cnt), matmul(matr,g0ws26(:,cnt))
     enddo
     write(99,*)'Ge   0 0 0 '
     write(99,27)'Ge ',g01
     write(99,27)'Ge ',g02
     write(99,27)'Ge ',g03
     write(99,27)'Ge ',g01+g02
     write(99,27)'Ge ',g03+g01
     write(99,27)'Ge ',g03+g02
     write(99,27)'Ge ',g03+g01+g02
     close(98)
     close(99)

77 format(299(1x,f9.4))
88 format(299(f10.6))

 end subroutine make_grids

end program FOCEX


