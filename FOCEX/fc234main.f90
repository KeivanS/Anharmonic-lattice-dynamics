program FC234
!
! todo:
! add option of reading existing fc files and fitting the remaining ones.
! also correct correspondingly the non-analytical part in THERMSPLIT
!
! add CS for fc3 fc4 so that the 2-body terms remain; or select the same 2-body terms as
! fc2's, for example as analytical derivatives of the morse potential fitted to the VASP data when
! only one atoms is moved at least 4 times +/- 0.1, +/- 0.05 (if not 6: +/- 0.15) plus a short-range term
! eventually may need to replace include_fc=2 by reading the fitted FC234's from an existing file
! and completing by the missing ones.
! ... CS can be implemented later. Switch to CS is FCs are too few and constraints too many.
! Reading born charges or fit without Born charges after subtraction of the Coulomb-ewald terms.
! in both cases, the file dielectric.params needs to be present and initialized (if needed, by Q_Born=0)
! Add the Coulomb contribution to FC3 and FC4
! Add reading of the stress tensor from DFT in case some adjustments of the equilibrium posititons and
! lattice parameters are needed (requires high cutoff in the DFT calculations to be reliable)
! Add the ability to eliminate some FCs using the invariance relations (needs testing) using svd
! we need to first subtract from the forces the long-range Coulombic part in order
! to insure the remaining part is short-ranged.
! NA part of FC2s=(Zi.q)(Zj.q)/eps q^2 corresponding to force: Fi=- sum_q (Zi.q)(Zj.q)/eps q^2 e^iq.Rij u_j
! need to calculate this force, subtract it from the DFT input force, and fit the remainder.
!------------------------------------------------
!! Program to extract force constants from ab initio force-displacement data
!! eventually with imposed (linear) constraints of symmetry
!! INPUT FILES : structure.params, dielectric.params, POSCARi, OUTCARi (i=1,2,...)
!! in structure.params file user specifies the specifications of the primitive cell and its atoms
!! the range of the force constants of rank 3 and 4 to be kept, whether or not translational,
!! rotational and Huang invariance constraints are imposed, specifications(type mass name) of
!! the atoms in the primitive cell and their reduced coordinates in units of translations of i
!! the conventional cell.
!! In POSCARi file the super cell structure and the equilibrium position of its atoms are specified
!! similar to VASP POSCAR format
!! Finally, OUTCARi constains the atomic displacement and corresponding forces on each atom in the
!! supercell are specified. It is obtained using a postprocessing utility of VASP or QE output files.
!! OUTPUT FILES: log.dat contains log of the run and intermediate outputs
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
! points defined in POSCAR : First need to solve F(u)=0 to get u and then use:
! phi3 = phi3 + u0 phi4; phi2 = phi2+ u0 phi3 + u0 u0 phi4/2;
! phi1 = phi1 + u0 phi2 + u0 u0 phi3/2 + u0 u0 u0 phi4/6; phi4 unchanged.
! and also apply the iterative Newton's algorithm to find the eq. positions
!* must also give as output the equilibrium volume by calculating the stress.
! check consistency of POSCAR and OUTCAR: differences must be less than 10%
!
use ios
use lattice
use params
use atoms_force_constants
use ewald
use svd_stuff
use geometry
use fourier
use kpoints !, only : kpc, wk, shift,nc,nkc,kibz,wibz
use born
implicit none
integer i,j,tau,ir,nfctot,rank,iunit,tra,rot,hua,ngr2, &
&      n_constr,n_kept,icfg,sc,nrgrd,ncfg,snap,ngrd,uio,nkernel,imax,n_hom
character xt*1,fn*7,poscar*7,outcar*7,now*10,today*8,zone*5
logical ex
real(8), allocatable :: afrc(:,:),bfrc(:),ewald_force(:,:),  &
&        us(:,:),vs(:,:),ws(:) !,project(:,:) !,grid(:,:),weights(:)newamat(:,:),newbmat(:),
real(8), allocatable :: newk(:,:) !,kibz(:,:),wibz(:)
real(8), allocatable :: mat(:,:),qmat(:),auxr(:),auxg(:),a_hom(:,:)
integer, allocatable :: frc_constr(:)
real(8) error,ermax,sd(4),sig, volmax
real(8) alfa,energy0,energy1,xx(3),ewaldsum,emad,deteps3,norm
real tim
type(vector) b01,b02,b03,as1,as2,as3,sx(26)

! for testing append....
  allocate(auxr(9)); auxr=10
  allocate(auxg(3)); auxg=3
!  allocate(qmat(7));
  print*,'auxr=',auxr
  print*,'auxg=',auxg

 auxr=reshape(auxr,shape=(/13/),pad=auxg)
! ! auxr=reshape(auxr,shape=(/size(auxr)+size(auxg)/),pad=auxg)
!  call pad_array(3,auxr(1:3),3,auxg,qmat)
!  print*,'qmat=',qmat
!  call pad_array(3,auxr(1:3),3,auxg,auxr)
! auxr=reshape(auxr,shape=(/size(auxr)+size(auxg)/),pad=auxg)
  write (*,*)"pad: auxr=",auxr

!  call append_array(2*auxr,auxg,qmat)
  print*,'auxr=',auxr
  print*,'qmat=',qmat

! write (ulog,*)"auxr"
! write (ulog,99)auxr
! write (ulog,*)'auxg'
! write (ulog,99)auxg
! write (ulog,*)'qmat'
! write (ulog,99)qmat
!  stop

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(tim)

 open(utimes,file='times.dat' ,status='unknown')
 open(umap  ,file='maps.dat'  ,status='unknown')
 open(umatrx,file='amatrx.dat',status='unknown')
 open(ucor  ,file='corresp.dat',status='unknown')
 open(ufco  ,file='lat_fc.dat',status='unknown')


 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' This program was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(utimes,'(a,f10.4)')' At the start of the Program FC234 the  TIME IS ',tim

  call set_maxterms   ! no need to read them anymore, the defaults are incremented if necessary


! read from structure.params the atoms in prim cell in reduced units, assign their mass, type, tau;
! also read specs of the prim cell, and flags for range of FC3 and FC4, enforcement of invariance relations
  call read_structure
! outputs: lattparams,primlat,nshells,include_fc,itrans,irot,natom_prim_cell,mas,atname,
! atompos0 (reduced coords within the primitive cell) and object atom0(natom_prim_cell)


! write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' Program FC234 was launched at ',today(1:4)//'/'  &
!&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
! write(ulog,'(a,f10.4)')' STARTING TIME OF THE PROGRAM                   IS ',tim

! reads dielectric constant and Born charges (=0 if unknown)
  call read_dielectric
  print*,' file dielectric.params read'

  call write_out(ulog,'latticeparameters',latticeparameters)
  call write_out(ulog,'primitive lattice',primitivelattice)
! write(ulog,*)'atompos0_d'
! do i=1,natom_prim_cell
!    call write_out(ulog,'atompos0 ',atompos0(:,i))
! enddo
  call write_out(ulog,'# of different atom types ',natom_type)
  call write_out(ulog,'# of each type in primcel ',natom)
  call write_out(ulog,'# mass of different types ',mas)

 call cpu_time(tim)
 write(utimes,'(a,f10.4)')' READ_INPUT CALLED, TIME                       IS ',tim


! inputs: natom_prim_cell,latticeparameters,primitivelattice,atom_type,atompos0
  call make_unitcell
! outputs: r0i, atompos (cart coord of atoms and natoms neighbor shells),
! iatomcell, iatomcell0, iatomneighbor, iatomop, ... and symmetry operations
! sets nsmax, according to the default range for fc2s, rcut(2) which is 15Ang

 call cpu_time(tim)
 write(utimes,'(a,f10.4)')' make_unitcell, TIME                          IS ',tim
 write(utimes,*)' make_unitcell called, atmpos, neighbs, iatmcell calculated'

! inputs: rcut(2) (and atompos through module atoms_force_constants)
  call set_neighbor_list(rcut(2),maxshell)
! outputs: atom0%equilb_pos,shells%no_of_neighbors , rij, neighbors(:)%tau and n and maxshell
! maxshell is the actual # of shells within rcut(2). It is not nsmax; maxshell < maxneighbors
 write(ulog,*)' After set_neighbors: maxshell=',maxshell,' while maxneighbors=',maxneighbors
 write(ulog,*)' ***************************************************************************** '
 write(ulog,*)' rcut(2) for FC2s set to be=',rcut(2),' corresponding to new maxshell=',maxshell
 write(ulog,*)' ************************* Writing the neighborlist ************************** '
 write(ulog,*)' Now the actual number of shells within the largest WS cell is set...'


  if (maxval(nshells(2,:)) .gt. maxshell) then
     nshells(2,:)=maxshell
     write(*   ,*)'nshells2 reduced from ',maxval(nshells(2,:)),' to ',maxshell
     write(ulog,*)'nshells2 reduced from ',maxval(nshells(2,:)),' to ',maxshell
  endif

  call write_neighbors

  call write_atompos


  if(fc2flag.eq.0) then  ! calculate the default # of shells of FC2 based on largest supercell

! of available superells, takes the one with largest volume to set the range of FC2= largest
! center to WS boundary distance; also collect all forces and displacements

     call find_WS_largest_SC(imax,volmax)
! output is the translation vectors of the largest supercell rsi

! create a grid of primitive translations inside the WS cell of supercell, rgrid(3,nrgrid), along with their weights
     call make_grid_weights_WS(r01,r02,r03,rs1,rs2,rs3,nrgrid,'r',rws26)
! rws26(output) are the 26 shortest vectors used to define the WS of the largest supercell
! in addition depending on 'r' or 'g' the grid of primitive translations inside the WS of supercell
! is calculated and stored in rgrid(3,nrgrid) or ggrid(3,nggrid) along with their weights rws_weights

     call update_nshells2(nrgrid,rgrid)
     write(ulog,*)'fc2flag=0; use largest WS; nshells2 updated to ',nshells(2,1:6)

  endif

!!!!!

!For each atompos in the loop for dynmat etc, we need to know to which supercell atom or its image it corresponds!!
!So need to define somemap(natoms) =super_cell index of that atom (used for the displacement)

!!!!!




! define the dimensions of the FC arrays, for each rank if include_fc.ne.0
! collects identical force_constants in groups of shells, identifies the irreducible ones and
! finds the mapping between the irreducible ones and all FCs in that group
  call setup_maps
! inputs: outputs of force_constants_init: iatomop etc...
! outputs: map structure including for rank, the groups of indep FCs and the full FCs
! maxterms, nshells, ngroups, nterms, estimate of inv_constraints, nindepfc

  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after setup_maps and before write_correspondance IS ',tim

  call write_lat_fc(map(:)%ntotind,map(:)%ntot)

!--------------------------------------------------------------------
! read FCs from a file if already present

! outputs: ngrnk(rank), irreducible fcs:fcrnk(rank,ngrnk(rank))
  do rank=1,4
  if (include_fc(rank) .eq. 2 ) then
     write(xt,'(i1)')rank
     fn = 'fc'//xt//'_irr.dat'
     iunit = ufco+rank
     inquire ( file=fn, exist=ex)
     if (ex) then
        open(iunit,file=fn,status='old')
! this rank is to be read and is not included in svd extraction
        ngrnk(rank)=map(rank)%ntotind; allocate(fcrnk(rank,ngrnk(rank)))
        call read_fcs_2(iunit,fn,rank,fcrnk(rank,:),ngrnk(rank))
     endif
  endif
  enddo

! if read, may need to check compatibility with given supercells?!?

!--------------------------------------------------------------------

! before imposing the constraints, first define the number of columns in matrix A
 allocate(keep_grp2(map(2)%ngr)) ! which groups to keep based on the vectors rws26
! keep_grp2(i)=1 is the list of FC2s within WS of supercell defined by 2ws26
! find the FC2s within the SC and their number, which is size_kept_fc2
 if(fc2flag.eq.0) then
     call setup_FC2_in_supercell !(keep_grp2,size_kept_fc2)
 else  ! use the default values found by setup_maps
     keep_grp2=1  ! keep everything
     size_kept_fc2=map(2)%ntotind !sum(keep_grp2)  not sum of the groups but sum of all rnk2 indep terms
 endif

! calculate the new nindepfc based on keep_grp2
! write(ulog,*)' new ngr(fc2) based size_kept_fc2=',map(2)%ngr-sum(map(2)%ntind)+size_kept_fc2
 write(ulog,*)' Old(*2) and new number of fc2:size_kept_fc2=',sum(map(2)%ntind),map(2)%ntotind,size_kept_fc2
 write(ulog,*)' Old and New # of kept groups_fc2=sum(keep_grp2)=',map(2)%ngr,sum(keep_grp2)

! change ngr from setup_maps
 nindepfc=0
 do rank=1,4
  if ( include_fc(rank) .eq. 1 ) then ! exclude them if they already exist and will be read

       if(rank.eq.2) then
          nindepfc = nindepfc + size_kept_fc2 ! this is the size used for FC2s; size in setup maps is ignored
       else
          nindepfc = nindepfc + map(rank)%ntotind
       endif
       write(ulog,*)'nindep_rnk,nindep_tot=',map(rank)%ntotind,nindepfc
   endif
 enddo

 write(ulog,*)' MAIN: total # of independent FCs of all rank, nindepfc=',nindepfc

 allocate(frc_constr(fdfiles))  ! for each of the supercell files

 ! Read force-displacements from FORCEDISPi (OUTCARi) and construct aforce,bforce matrices
 call read_snaps_set_aforce(frc_constr,nconfigs)

  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after read_snaps_set_aforce               IS ',tim

!
!--------------------------------------------------------------------
! depending on the constraints, setup the homogeneous part of the A and B matrices:
 if(born_flag.lt.2) then ! use a real space rteatment of ASR and rotational invariances

    call set_translational_inv_constraints
    write(ulog,*)'Number of translational invariance constraints is=',transl_constraints

    call cpu_time(tim)
    write(utimes,'(a,f10.4)')' TIME after set_translational_inv_constraints IS ',tim

    call set_rotational_inv_constraints
    write(ulog,*)'Number of rotational invariance constraints is=',rot_constraints

    call cpu_time(tim)
    write(utimes,'(a,f10.4)')' TIME after set_rotational_inv_constraints    IS ',tim

    call set_huang_inv_constraints
    write(ulog,*)'Number of Huang    invariance constraints is=',huang_constraints

    call cpu_time(tim)
    write(utimes,'(a,f10.4)')' TIME after set_Huang_inv_constraints         IS ',tim

  elseif(born_flag .eq. 2) then  ! set up ASR in k-space

  endif

! now that aforce & bforce matrices are setup, we can do SVD and project on
! onto the kernel of homogeneous part of amatrix

! inputs: inv_constraints,force_constraints, aforce, bforce and invariance parts of amat
  call include_constraints_remove_zeros
! outputs are amat and bmat and their dimensions:dim_al,dim_ac, which go to SVD

  write(ulog,*)'dim_al,dim_ac set to ',dim_al,dim_ac

! give warnings if any column in amat is totally zero
  call check_zero_column(dim_al,dim_ac,amat)


  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after include_constraints and zero column   IS   ',tim

! write the a and b matrices
  if (verbose) then
     write(umatrx,*)' before call to svd, amat and bmat are:'
     do i=1,dim_al
        write(umatrx,17)(amat(i,j),j=1,dim_ac),bmat(i)
     enddo
  endif
17 format(99(1x,g9.2))

! do the SVD decomposition; beware amat is overwritten
  if(allocated(fcs)) deallocate(fcs)
  if (allocated(sigma)) deallocate(sigma)
  allocate(fcs(dim_ac),sigma(dim_ac))

  uio = 349
  if(enforce_inv.eq.1 .or. itemp.eq.1) then ! elimination wil be used; need to save a_hom

      open(uio,file='elimination.dat',status='unknown')
! separate homogeneous part of amat from inhomogeneous part by finding the first non-zero elt of bmat
!      call find_first_nonzero(dim_al,bmat,n_hom)
      n_hom=inv_constraints
      write(ulog,*)'MAIN: Using kernel projection: enforce_inv and itemp are=',enforce_inv,itemp
      write(ulog,*)'MAIN: position of the first non-zero elt of bmat=',n_hom
      allocate(a_hom(n_hom,dim_ac),kernelbasis(dim_ac,dim_ac))
      a_hom=amat(1:n_hom,1:dim_ac)
      call write_out(uio,'#####     HOMOGENEOUS PART OF A ',a_hom)
      call get_kernel(n_hom,dim_ac,a_hom,svdcut,nkernel,kernelbasis,uio)
      write(ulog,*)' Kernel of A, of dim ',n_hom,'x',dim_ac,' is of dimension ',nkernel
      write(ulog,*)' The number of eliminated components is ',dim_ac-nkernel
      deallocate(a_hom)

  endif

  allocate(bfrc(dim_ac))  ! dummy for xout
  if (itemp.eq.1) then

      write(*,*) ' Temperature of ',tempk,' will be implemented'
      allocate(mat(dim_ac,dim_ac),qmat(dim_ac))
!      2023-07-03 16:24:12 ---> Checking the energy value as FOCEX code is not giving proper result (logged by Bikash. T.)
      write(*,*) "Energy values from the OUTCAR1 are: ",energies
      write(*,*) "Dimension of energies are: ",shape(energies)
      write(*,*) "The value of nlines are: ",nlines
      write(*,*) "The value of nconfigs is: ",nconfigs
      call implement_temperature(dim_al-n_hom,dim_ac,amat(n_hom+1:dim_al,1:dim_ac), &
&         bmat(n_hom+1:dim_al),nconfigs,energies,nlines,tempk,mat,qmat)
      call solve_svd(dim_ac,dim_ac,mat,qmat,svdcut,bfrc,sigma,uio)
      call write_out(ulog,'TEMP: SVD solution before kernel projection',bfrc)
      call project_on(dim_ac,nkernel,bfrc,kernelbasis(1:dim_ac,1:nkernel),fcs)
      deallocate(kernelbasis,mat,qmat)
  else

      if(enforce_inv.eq.1) then ! elimination wil be used; need to save a_hom
         write(ulog,*)'MAIN: Using svd of the force amatrix since enforce_inv and itemp are=',enforce_inv,itemp
!         call solve_svd(dim_al-n_hom,dim_ac,amat(n_hom+1:dim_al,1:dim_ac), &
!&             bmat(n_hom+1:dim_al),svdcut,bfrc,sigma,uio)
         call svd_set(dim_al-n_hom,dim_ac,amat(n_hom+1:dim_al,1:dim_ac),  &
 &            bmat(n_hom+1:dim_al),bfrc,sigma,svdcut,error,ermax,sig,'svd-elim.dat')
         call write_out(ulog,'Enforce=1: SVD solution before kernel projection',bfrc)
         write(ulog,*)'MAIN: performing projection on the kernel basis'

         call project_on(dim_ac,nkernel,bfrc,kernelbasis(1:dim_ac,1:nkernel),fcs)
         call write_out(ulog,'Enforce=1: Solution After kernel projection',bfrc)
         deallocate(kernelbasis)
      else
         call svd_set(dim_al,dim_ac,amat,bmat,fcs,sigma,svdcut,error,ermax,sig,'svd-all.dat')
         write(ulog,*)'After svdall, || F_dft-F_fit || / || F_dft || =',sig
      endif

  endif
  close(uio)

  call write_invariance_violations(ulog,dim_ac,fcs)

  if(include_fc(1).eq.1) open(ufit1 ,file='fc1_irr.dat',status='unknown')
  if(include_fc(2).eq.1) open(ufit2 ,file='fc2_irr.dat',status='unknown')
  if(include_fc(3).eq.1) open(ufit3 ,file='fc3_irr.dat',status='unknown')
  if(include_fc(4).eq.1) open(ufit4 ,file='fc4_irr.dat',status='unknown')
  if(include_fc(1).eq.1) open(ufc1 ,file='fc1.dat',status='unknown')
  if(include_fc(2).eq.1) open(ufc2 ,file='fc2.dat',status='unknown')
  if(include_fc(3).eq.1) open(ufc3 ,file='fc3.dat',status='unknown')
  if(include_fc(4).eq.1) open(ufc4 ,file='fc4.dat',status='unknown')

   call write_independent_fcs(dim_ac,sigma,sd,ulog)

   call write_output_fcs

  close(ufit1)
  close(ufit2)
  close(ufit3)
  close(ufit4)
  close(ufc1)
  close(ufc2)
  close(ufc3)
  close(ufc4)

   deallocate(sigma,amat,bmat,bfrc) !,fcs)

! calculate the phonon dispersion as output test: requires kpbs.in as input
  open(ugrun,file='mode_grun.dat',status='unknown')
  call set_dynmat !(map(2))

  call cpu_time(tim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After band structure,               TIME IS ',tim



  call date_and_time(date=today,time=now,zone=zone)
  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' ENDING TIME OF THE PROGRAM                   IS ',tim
  write(ulog,'(a,f10.4)')' ENDING TIME OF THE PROGRAM                   IS ',tim


3 format(i5,2x,9(1x,g12.5))
5 format(9(1x,f11.6))
7 format(a,2(2x,g12.5),2x,9(1x,g11.4))
8 format(i4,3(2x,f10.5),3x,i2,2x,'(',3(1x,i2),1x,')',2x,i4)
9 format(i5,9(1x,f9.4))
99 format(1000(1x,g10.3))

 write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' Program  FC234 ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' Program FC234 ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

  write(ulog,7)' Average and largest errors in force=',error,ermax
  write(ulog,7)' Average errors in FC234=',sd

  close(ulog)
  close(ufco)
  close(ucor)
  close(umap)
  close(umatrx)
  close(utimes)

  contains

 subroutine test(new,n,m)
 integer n,m
 real(8), allocatable :: new(:,:)
 allocate(new(n,m))
 new =2
 end subroutine test

!===================================================
 subroutine append_array2(a,b,c)
!! appends array b to the end of a and stores the result in c
 implicit none
! integer, intent(in) :: na,nb
! real(8), intent(in):: a(na),b(nb)
 real(8), intent(in):: a(:),b(:)
! real(8), allocatable, intent(inout):: c(:)
 real(8), allocatable :: c(:)
! real(8) :: c(size(a)+size(b))
 real(8), allocatable :: aux(:)
 integer ncc,naa,nbb

 naa=size(a);
 nbb=size(b);
 ncc=naa+nbb
! if (allocated(c)) deallocate (c)
 allocate(aux(ncc)) ! this is to allow calls like append_array(a,b,a)
! allocate(c(ncc))
 aux=reshape(a,shape=(/naa+nbb/),pad=b)
! aux=reshape(a,(/naa+nbb/),pad=b)
 if (allocated(c)) deallocate (c)
 allocate(c(ncc))
 c=aux
 deallocate(aux)
! c=reshape(a,(/na+nb/),pad=b,order=(/1,1/))

 end subroutine append_array2
!========================================
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
 use atoms_force_constants, only : natom_super_cell,displ,force,atom_sc
 implicit none
 integer, intent(out) :: frc_cnstr(fdfiles),nconfigs
! real(8), intent(out) :: energies(nconfigs),disp(3,natom_super_cell,nconfigs),force(3,natom_super_cell,nconfigs)
! real(8), allocatable, intent(out) :: energies(:) !,disp(:,:,:),force(:,:,:)
 real(8), allocatable :: dsp(:,:,:),frc(:,:,:),energy(:),afrc(:,:),bfrc(:),aux3(:,:,:),aux4(:,:,:),aux(:)
 integer i,imax,ncfg,at,j,nfctot,n1,n2
 type(vector) as1,as2,as3
 character xt*1,poscar*7, outcar*10
 real tim

! this loop is to calculate # of force_constraints in order to allocate aforce
     nconfigs=0 ; nfctot=0
     structures: do i=1,fdfiles

        write(xt,'(i1)')i
        write(ulog,*)'=================== Reading supercell #',i,'========================'
        poscar='POSCAR'//xt
        call read_crystal(poscar) ! equilibrium atomic positions in supercell from POSCARi

        call make_r0g

! make sure of consistency POSCAR with input structure
        call identify_atoms_in_supercell

! writing the coordinates of primitive cells in the supercell
! reference of each prim cell (n) is atom #1
        open(111,file='supercell-'//xt)
        write(111,5)rs1
        write(111,5)rs2
        write(111,5)rs3
        do j=1, natom_super_cell
           if(atom_sc(j)%cell%tau .ne. 1) cycle
           write(111,9)j,atom_sc(j)%cell%n(1)*r01+atom_sc(j)%cell%n(2)*r02  &
&                    +atom_sc(j)%cell%n(3)*r03,atom_sc(j)%equilibrium_pos
        enddo
        close(111)

! generating rgrid and ggrid for Ewald force calculations in the supercell
        if(born_flag.ne.0) then
           call make_grid_weights_WS(r01,r02,r03,rs1,rs2,rs3,nrgrid,'r',rws26)
           write(ulog,*)' original rgrid(:,j),rws_weights(j) '
           do j=1,nrgrid
              write(ulog,9)j,rgrid(:,j),rws_weights(j)
           enddo
           write(ulog,*)' sum of rws_weights(j) ',sum(rws_weights)
           call make_grid_weights_WS(gs1,gs2,gs3,g01,g02,g03,nggrid,'g',gws26)
           write(ulog,*)' original ggrid(:,j),gws_weights(j) '
           do j=1,nggrid
              write(ulog,9)j,ggrid(:,j),gws_weights(j)
           enddo
           write(ulog,*)' sum of gws_weights(j) ',sum(gws_weights)
        endif

! was OUTCAR before and changed to FORCEDISP
        outcar='FORCEDISP'//xt

        call count_configs(outcar,ncfg)
        nconfigs = nconfigs + ncfg

        frc_cnstr(i)=3*natom_super_cell*ncfg
        allocate( energy(ncfg),displ(3,natom_super_cell,ncfg),force(3,natom_super_cell,ncfg) )
        write(*,*)' Allocation of energy, dsp and frc arrays DONE!'
        call read_force_position_data(outcar,ncfg,energy,displ,force)

!       if(frc_cnstr(i) .ne. ncfg*natom_super_cell*3) then
!          write(ulog,*)'READ_SNAPS: inconsistency in frc_cnstr ',frc_cnstr(i),ncfg*natom_super_cell*3
!          stop
!       endif

  call cpu_time(tim)
  write(utimes,'(a,i3,a,f10.4)')' STRUCTURE ',i,' after read_force_position_data time IS ',tim

! 3 ways: if=0 ignore, if=1 use real space ewald, if=2 use Fourier transforms
        call subtract_coulomb_force(born_flag,ncfg,displ,force)
! now frc contains the short-range (Non-Coulomb) part of the force constants which we can fit

  call cpu_time(tim)
  write(utimes,'(a,i3,a,f10.4)')' STRUCTURE ',i,' after subtract_coulomb time IS ',tim

! now append new snapshot data to existing one
        if(i.eq.1) then
            allocate(energies(ncfg))
!            call allocate_pos(natom_super_cell,ncfg)   ! for displ and force
            energies=energy
!            force=frc
!            displ=dsp
            allocate( aforce(frc_cnstr(i),nindepfc),bforce(frc_cnstr(i)) )
            call set_force_displacement_matrix(frc_cnstr(i),aforce,bforce)
            if(nconfigs.ne.ncfg) print*,'APPENDING ERROR, ifile=',i
        else ! append new arrays for i.ge.2 to the end of existing ones: energies,force,disp
!            call append_array(energies,energy(1:ncfg),energies)
             energies=reshape(energies,shape=(/size(energies)+size(energy)/),pad=energy)
!             n1=size(energies); n2=size(energy)
!             allocate(aux(n1+n2))
!             aux(1:n1)=energies;
!             aux(1+n1:n2+n1)=energy
!             deallocate(energies); allocate(energies(n1+n2))
!             energies=aux
!             deallocate(aux)
            if(allocated (aux3)) deallocate(aux3)
            if(allocated (aux4)) deallocate(aux4)
            allocate(aux3(3,natom_super_cell,nconfigs+ncfg))  ! force
            allocate(aux4(3,natom_super_cell,nconfigs+ncfg))  ! displ
            aux3(:,:,nconfigs+1:nconfigs+ncfg)=force
            aux4(:,:,nconfigs+1:nconfigs+ncfg)=displ
!            do at=1,natom_super_cell
!            do j=1,3
!               call append_array(force(j,at,:),frc(j,at,1:ncfg),force(j,at,:))
!             n1=size(force(j,at,:)); n2=size(frc(j,at,:))
!             allocate(aux(n1+n2))
!             aux(1:n1)=force(1:n1);
!             aux(1+n1:n2+n1)=frc
!             deallocate(force); allocate(force(n1+n2))
!             force=aux
!             deallocate(aux)

!               call append_array(displ(j,at,:),dsp(j,at,1:ncfg),displ(j,at,:))
!             n1=size(displ); n2=size(dsp)
!             allocate(aux(n1+n2))
!             aux(1:n1)=displ;
!             aux(1+n1:n2+n1)=dsp
!             deallocate(force); allocate(displ(n1+n2))
!             displ=aux
!             deallocate(aux)

!               force(j,at,:)=reshape(force(j,at,:),(/size(force(j,at,:))+ncfg/),pad=frc(j,at,1:ncfg))
!               displ(j,at,:)=reshape(displ(j,at,:),(/size(displ(j,at,:))+ncfg/),pad=dsp(j,at,1:ncfg))

!            enddo
!            enddo

            allocate(afrc(frc_cnstr(i),nindepfc),bfrc(frc_cnstr(i)))
            write(ulog,*)' calling set_force_displacement_matrix, for str=',i
            call set_force_displacement_matrix(frc_cnstr(i),afrc,bfrc)
            do j=1,nindepfc
!               call append_array(aforce(:,j),afrc(:,j),aforce(:,j))
               aforce(:,j)=reshape(aforce(:,j),(/size(aforce(:,j))+3*natom_super_cell/),afrc(:,j))
            enddo
!            call append_array(bforce,bfrc,bforce)
            bforce=reshape(bforce,(/size(bforce)+3*natom_super_cell/),bfrc)
            deallocate(afrc,bfrc)
        endif


        write(ulog,*)' set_force_displacement_matrix called, aforce being filled now'
        write(ulog,*)' nfctot,nfctot+frc_constr(i)=',nfctot,nfctot+frc_cnstr(i)
        nfctot = nfctot + frc_cnstr(i)

     enddo structures

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
!  &           nrgrid,rgrid,rws_weights,nggrid,ggrid,gws_weights,map_rtau_sc)

5 format(9(1x,f12.6))
9 format(i5,9(1x,f12.6))

 end subroutine read_snaps_set_aforce


end program FC234

!===========================================================
!!! do the SVD of the overlap
!! allocate( us(n_constr,ngr),vs(ngr,ngr),ws(ngr),mw(ngr) )
!! us = a_hom
!! call svdcmp(us,n_constr,ngr,n_constr,ngr,ws,vs)

! allocate( us(n_constr,n_constr),vs(n_constr,n_constr),ws(n_constr) )
! us = overl
! call svdcmp(us,n_constr,n_constr,n_constr,n_constr,ws,vs)

! truncate and keep the highest eigenvalues of ws ! is that correct??
! call truncate(ngr,ws,n_kept,mw)
! write(ulog,*)' the largest  ws is=',ws(mw(1))
! write(ulog,*)' the smallest ws kept and its index is=',n_kept,ws(mw(n_kept))
! write(ulog,*)' their ratio is=',ws(mw(n_kept))/ws(mw(1))

! for the kept eigenvalues calculate the correction or the residual FCS
! by projecting the initial solution onto the kernel of the constraint matrix
! here we are constructing the projection matrix which should later multiply X0
! allocate( project(n_constr,n_constr) )
! project=0
! do i=1,n_constr
!    project(i,i)=1
! do j=1,n_constr
!    do n3=1,n_kept !+1,n_constr
!       project(i,j)=project(i,j)-vs(i,mw(n3))*vs(j,mw(n3))
!    enddo
! enddo
! enddo
! this was the kernel remains, but then if kernel truncation should be done at small thresholds!
! we need to multiply the final svd solution vector by project, to remain in the kernel of C

! stop

! nc(1)=2;nc(2)=3;nc(3)=4;
! nkc=nc(1)*nc(2)*nc(3)
! allocate(kpc(3,nkc),wk(nkc),newk(3,nkc))
! shift=0d0

! this is needed for the FBZ (Wigner-Seitz cell) construction; output is gshells(3,26)
!  call get_26shortest_shell(g01,g02,g03,gws26,b01,b02,b03)


!  call allocate_tetra(nc1,nc2,nc3,100)

!  call make_kp_reg(nc,g01,g02,g03,shift,kpc,wk)

!  call fit_in_BZ(nkc,kpc,gshells,newk)

!  call structure_factor(v2a(r01),nkc,newk,sig)
!  call structure_factor(v2a(r02),nkc,newk,sig)
!  call structure_factor(v2a(r03),nkc,newk,sig)

! stop

! rx1= v2a(r01)
! rx2= v2a(r02)
! rx3= v2a(r03)

!  do i=3,3
!     alfa = 0.1*2**i
!!    call ewald(alfa,natom_prim_cell,pos,natom_prim_cell,zborn,epsil,energy)
!     xx=(/0.000001d0,0d0,0d0/)
!    energy0=ewaldsum(xx,rx1,rx2,rx3,alfa,epsil,emad)
!     xx=(/0.001d0,0d0,0d0/)
!    energy1=ewaldsum(xx,rx1,rx2,rx3,alfa,epsil,emad)
!    write(ulog,*)'alfa, dE/dx==',alfa,(energy1-energy0)/0.001
!     write(ulog,*)'================================================='
!  enddo
!  stop

!do i=3,3
!    eta=sqrt(pi*deteps3)/(volume_r**0.3333)  !*(2**i)/10d0
!    rcutoff_ewa=8.60*sqrt(deteps3)/eta !*volume_r**0.3333
!    gcutoff_ewa=8.60*eta*2/sqrt(deteps3) !/volume_r**0.3333
!    rcutoff_ewa=9.60/eta !*volume_r**0.3333
!    gcutoff_ewa=9.60*eta*2 !4*pi/volume_r**0.3333
!
!    write(ulog,7)'eta,rcut,gcut,etaR,G/eta=',eta,rcutoff_ewa,   &
!    &                gcutoff_ewa,eta*rcutoff_ewa,gcutoff_ewa/eta
!
!
!  call test_ewald
!
!enddo
! this is needed for the FBZ (Wigner-Seitz cell) construction; output is gshells(3,26)
! call get_26shortest_shell(g01,g02,g03,gws26,b01,b02,b03)
!stop


! tra=0; rot=0; hua=0
! if (itrans .ne. 0) tra= transl_constraints
! if (irot   .ne. 0) rot= rot_constraints
! if (ihuang .ne. 0) hua= huang_constraints
! write(ulog,*)' Number of invariance constraints   =',inv_constraints
! write(ulog,*)' Size of the homogeneous part of A  =',tra+rot+hua
! if (inv_constraints.ne.tra+rot+hua) then
!     inv_constraints = tra+rot+hua
!     write(ulog,*)' Setting # of invariance constraints=',inv_constraints
! endif

! if (dim_ac .le. inv_constraints) then
!    write(ulog,*)' the # of indept FCs is too small for doing any elimination',dim_ac
!    call warn(ulog)
!    write(ulog,*)' the # of constraints is larger than the number of FCs ',inv_constraints
!    write(ulog,*)' you need to release some of the constraints or increase the range'
!    call warn(6)
!    write(6   ,*)' the # of constraints is larger than the number of FCs ',inv_constraints
!    write(6   ,*)' you need to release some of the constraints or increase the range'
!    stop
! endif
! we remove the following because if one starts with reference positons in POSCAR
! which are not equilibrium position, then FC1 are not zero and lines with u=0
! for now do it without reduction, put the full A matrix
! then do the reduction and check validity
! then impose other invariance relations
! then put higher order terms
!----------------------- END FOR TESTING -----------------------
!--------------------------- FOR TESTING -----------------------
! this part only does the svd on forces without including invariance relations
!  dim_al = force_constraints
!  dim_ac = ngr
!  call remove_zeros(dim_al,dim_ac, aforce, bforce, newdim_al)
!  allocate(newamat(newdim_al,dim_ac),newbmat(newdim_al))
!  newamat(:,:) = aforce(1:newdim_al,1:dim_ac)
!  newbmat(:  ) = bforce(1:newdim_al)
!  deallocate(aforce,bforce)
!  dim_al = newdim_al
!
!! give warnings if any column in aforce is totally zero
!  call check_zero_column(dim_al,dim_ac,newamat)
!
!  allocate(fcs(dim_ac),sigma(dim_ac))
!
!  call cpu_time(tim)
!  write(utimes,'(a,f10.4)')' TIME before SVD  of aforce                   IS   ',tim
!
!! do the SVD decomposition
!
!  call svd_set(dim_al,dim_ac,newamat,newbmat,fcs,sigma,svdcut,error,ermax,sig,'svd-force.dat')
!
!  write(ulog,*)'After svd, || F_dft-F_fit || / || F_dft || =',sig
!
!  call cpu_time(tim)
!  write(utimes,'(a,f10.4)')' TIME after SVD calculation of aforce         IS ',tim
!
!! Add the correction term to the old solution
!
!
! allocate(amat(dim_al,dim_ac),bmat(dim_al),btes(dim_al))
! amat(1:transl_constraints,1:ngr) = atransl(1:transl_constraints,1:ngr)
! bmat(1:transl_constraints)     = 0d0
! amat(transl_constraints + 1 : transl_constraints + rot_constraints,1:ngr) =  &
! &                        arot(1:rot_constraints,1:ngr)
! bmat(transl_constraints + 1 : transl_constraints + rot_constraints) = &
! &                        brot(1:rot_constraints)
!
! write all fcs into one array:
! allocate(fcs(dim_ac))
! i = 0
! fcs(i+1:i+ngroups(2)) = fcs_2(1:ngroups(2))
! i = i + ngroups(2)
! fcs(i+1:i+ngroups(3)) = fcs_3(1:ngroups(3))
! i = i + ngroups(3)
! fcs(i+1:i+ngroups(4)) = fcs_4(1:ngroups(4))
! write(ulog,*)' ###############  TESTING OF FCS  ############# '
! do i=1,ngr
!    write(ulog,*) i,fcs(i)
! enddo
! btes = matmul(amat,fcs)
! write(ulog,*)' ###############  TESTING OF INVARIANCES  ############# '
! do i=1,dim_al
!    write(ulog,*) i,bmat(i),btes(i)
! enddo
!
!--------------------------- FOR TESTING -----------------------
! this part is for dynamical matrix/ewald calculations -------------------------
!    nkmesh=natom_super_cell/natom_prim_cell
! number of primitive vectors within the WS of the supercell also =nrgrid
!    nrmesh=nkmesh
!     allocate(kmesh(3,nkmesh),rmesh(3,nrmesh))
!     call K_mesh(b01,b02,b03,gs1,gs2,gs3,gws26,nkmesh,kmesh)
!     call K_mesh(as1,as2,as3,r01,r02,r03,rws26,nrmesh,rmesh)
! this is for testing s(gi)=0 except for G0i and 0
!     call structure_factor(v2a(rs1),nkmesh,kmesh,sig)
!     call structure_factor(v2a(r01),nkmesh,kmesh,sig)
!     call structure_factor(v2a(r03),nkmesh,kmesh,sig)
! remove zero lines in amat and bmat and update dim_al and amat and bmat accordingly
!  call remove_zeros(dim_al,dim_ac)

!  allocate(newamat(newdim_al,dim_ac),newbmat(newdim_al))
!  newamat(:,:) = amat(1:newdim_al,1:dim_ac)
!  newbmat(:  ) = bmat(1:newdim_al)
!  deallocate(amat,bmat)
!  allocate(amat(newdim_al,dim_ac),bmat(newdim_al))
!  amat = newamat
!  bmat = newbmat
!  dim_al = newdim_al
!  deallocate(newamat,newbmat)
!do n3=1,5
!   eta=0.1*2**n3
!
!     allocate(force_born(3,natom_super_cell,nconfigs))
!     do icfg=1,1 !nconfigs
!        write(ulog,*)'iconfig=',icfg
!        call ewaldforce(natom_super_cell,displ(:,:,icfg),fewald(:,:))
!     enddo
!     deallocate(force_born)
!     write(ulog,*)'MAIN: force_born called and exited successfully'
!
!     do j=1,natom_super_cell
!        write(ulog,5)j,displ(:,j,icfg),force(:,j,icfg)-fewald(:,j)
!     enddo
!
!enddo
!
!stop

