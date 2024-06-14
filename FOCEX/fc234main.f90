!=====================================================

program FC234
!
! todo:
!
! add option of reading existing fc files and fitting the remaining ones.
! also writeout the fitted FCs in the format of other codes: ALAMODE, PHONOpy SHENG
!
! add CS for fc3 fc4 so that the 2-body terms remain; or select the same 2-body terms as
! fc2's, for example as analytical derivatives of the morse potential fitted to the VASP data when
! only one atoms is moved at least 4 times +/- 0.1, +/- 0.05 (if not 6: +/- 0.15) plus a short-range term
! eventually may need to replace include_fc=2 by reading the fitted FC234's from an existing file
! and completing by the missing ones.
! Add the Coulomb contribution to FC3 and FC4
! Add reading of the stress tensor from DFT in case some adjustments of the equilibrium positions and
! lattice parameters are needed (requires high cutoff in the DFT calculations to be reliable)
!------------------------------------------------
!! Program to extract force constants from ab initio force-displacement data
!! eventually with imposed (linear) constraints of symmetry
!! INPUT FILES : structure.params, dielectric.params, POSCARi, OUTCARi (i=1,2,...)
!! in structure.params file user specifies the primitive cell and its atoms
!! the range of the force constants of rank 3 and 4 to be kept, whether or not translational,
!! rotational and Huang invariance constraints are imposed, specifications(type mass name) of
!! the atoms in the primitive cell and their reduced coordinates in units of translations of
!! the CONVENTIONAL cell.
!! In POSCARi file the super cell structure and the equilibrium position of its atoms are specified
!! similar to VASP POSCAR format
!! Finally, OUTCARi contains the atomic displacement and corresponding forces on each atom in the
!! supercell. It is obtained using a postprocessing utility (readoutcar) of VASP or QE output files.
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
use linalgb
implicit none
integer i,j,g,ti,rank,iunit, uio,nkernel,imax,n_hom,nt(maxrank),ntind(maxrank)  ,tau,taup
character xt*1,fn*7,now*10,today*8,zone*5,fni*11
logical ex
real(r15), allocatable :: xout(:),foldedk(:,:)
real(r15), allocatable :: mat(:,:),qmat(:)
integer, allocatable :: frc_constr(:)
real(r15) error,ermax,sd(4),sig, volmax,sf,dsf(3) ,q2(3),qr(3)!,x(3)
! real(r15) dyn_coul(3,3) ,ddn(3,3,3),q(3),deteps3 
real tim

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(tim)

! a=reshape((/-5,3,2,2,7,0,1,-1,-2/),(/3,3/))
! b=(/2,14,-4/)
! call lin_solve(a,b,x)
! write(*,*)'x=',x
! stop
!
! allocate(mat(2,3),a_hom(3,3))
!
! do i=1,3
! do j=1,2
!    mat(j,i)=i+3*j
! enddo
! do j=1,3
!    a_hom(j,i)=10*(i+3*j-1)
! enddo
! enddo
! call write_out(6,' mat  ',mat)
! call write_out(6,'a_hom ',a_hom)
! call append_array(mat,a_hom,mat)
! call write_out(6,' mat  ',mat)
! stop

 open(utimes,file='times.dat' ,status='unknown')
 open(umap  ,file='maps.dat'  ,status='unknown')
 open(umatrx,file='amatrx.dat',status='unknown')
 open(ucor  ,file='corresp.dat',status='unknown')


 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' Program FOCEX was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(utimes,'(a,f10.4)')' At the start of the Program FOCEX the  TIME  IS ',tim


! read from structure.params the atoms in prim cell in reduced units, assign their mass, type, tau;
! also read specs of the prim cell, and flags for range of FC3 and FC4, enforcement of invariance relations
   call read_structure
! outputs: lattparams,primlat,nshells,include_fc,itrans,irot,natom_prim_cell,mas,atname,
! atompos0 (reduced coords within the primitive cell) and object atom0(natom_prim_cell)


! reads dielectric constant and Born charges (=0 if unknown)
   call read_dielectric

! the number of shells is based on nshells in the structure.params file and should be large enough
! to include all atoms in the supercell
   call make_unitcell
! inputs: natom_prim_cell,latticeparameters,primitivelattice,atom_type,atompos0
! outputs: r0i, atompos (cart coord of atoms and natoms neighbor shells),
! iatomcell, iatomcell0, iatomneighbor, iatomop, ... and symmetry operations
! sets nsmax, according to the default range for fc2s, rcut(2)

  call write_atompos  !  atompos in .xyz format for visualization

!  call find_WS_largest_SC(imax,volmax)
!  call write_out(6,' rsc1 ',rs1)
!  call write_out(6,' gsc1 ',gs1)
!  call make_grid_weights_WS(r01,r02,r03,rs1,rs2,rs3,nrgrid,'r',r0ws26,rws26)
!  call make_grid_weights_WS(gs1,gs2,gs3,g01,g02,g03,nggrid,'g',gws26,g0ws26)
!4 format(a,3(1x,f9.4),1x,g13.6)
!  do i=1,5
!  do j=1,5
!  do g=1,5
!     q2=i*gs1+j*gs2+g*gs3 ! + a2v(0.2d0*(/1,2,3/))
!     qr=reduce_g(q2)  ! matmul(transpose(prim_to_cart),q2)/(2*pi)
!     call structure_factor_recip(q2,nrgrid,rgrid,rws_weights,sf,dsf)
!     if (is_integer(qr)) then
!        if(abs(sf-1).gt.1d-5) write(6,4)'dyncoul ERROR: qred,sf(q)=',qr,sf
!     else 
!        if(abs(sf).gt.1d-5) write(6,4)'dyncoul ERROR: qred,sf(q)=',qr,sf
!     endif
!     write(6,4)'dyncoul test: qred,sf(q)=',qr,sf
!  enddo
!  enddo
!  enddo
!
! stop

! k1 added for testing
! allocate(foldedk(3,nkc))

! rs1=nc(1)*r01
! rs2=nc(2)*r02
! rs3=nc(3)*r03

! call read_latdyn
! allocate(kpc(3,nkc),wk(nkc))
! gs1=(1d0/nc(1))*g01; gs2=(1d0/nc(2))*g02; gs3=(1d0/nc(3))*g03; 
! call make_grid_weights_WS(gs1,gs2,gs3,g01,g02,g03,nggrid,'g',gws26,g0ws26)
! call make_kp_reg(nc,g01,g02,g03,shift,kpc,wk)
! kpc=ggrid; if(nkc.ne. nggrid) write(*,*)'ERROR in the test nkc,nggrid=',nkc,nggrid

!  call fold_in_WS_BZ(nkc,kpc,g0ws26,foldedk)   ! also calls get_weights
!  call get_kpfbz(nkc,kpc,g0ws26,foldedk) 

! Now calculate the weigths
! call get_weights3(nkc,ggrid) !foldedk) 

! stop

! of available supercells, takes the one with largest volume to eventually set
! the range of FC2= largest center to WS boundary distance; record rs1,rs2,rs3
! output is rs1,rs2,rs3 of this largest supercell
   call find_WS_largest_SC(imax,volmax)


 call cpu_time(tim)
 write(utimes,'(a,f10.4)')' make_unitcell, TIME                          IS ',tim

    call set_neighbor_list !(rcut(2),maxshell)

! outputs: atom0%equilb_pos,shells%no_of_neighbors , rij, neighbors(:)%tau and n and maxshell
! maxshell is the actual # of shells within rcut(2). It is not nsmax; maxshell < maxneighbors
 write(ulog,*)' After set_neighbors: maxshells=',maxshells,' while maxneighbors=',maxneighbors
 write(ulog,*)' ***************************************************************************** '
! write(ulog,*)' rcut(2) for FC2s set to be=',rcut(2),' corresponding to new maxshell=',maxshell
! write(ulog,*)' ************************* Writing the neighborlist ************************** '
 write(ulog,*)' Now the actual number of shells within the largest WS cell is set...'


!   if (fc2flag.ne.0 .and. maxval(nshells(2,:)) .gt. maxshells) then
!     nshells(2,:)=maxshells
!     write(*   ,*)'nshells2 reduced from ',maxval(nshells(2,:)),' to ',maxshells
!     write(ulog,*)'nshells2 reduced from ',maxval(nshells(2,:)),' to ',maxshells
!   endif

  call write_neighbors


! create a grid of primitive translations inside the WS cell of supercell+boundaries, rgrid(3,nrgrid),
! along with their weights  in case they are on the WS boundary
! rws26(output) are the 26 shortest vectors used to define the WS of the largest supercell
! in addition depending on 'r' or 'g' the grid of primitive translations inside the WS of supercell
! is calculated and stored in rgrid(3,nrgrid) or ggrid(3,nggrid) along with their weights rws_weights
     call make_grid_weights_WS(r01,r02,r03,rs1,rs2,rs3,nrgrid,'r',r0ws26,rws26)
     call write_lattice(size(rws26,2),rws26,'r_supercell.xyz') ! grid of primitive translations in the WS supercell
     write(ulog,*)' original rgrid(:,j),rws_weights(j) '
     do j=1,nrgrid
        write(ulog,9)j,rgrid(:,j),rws_weights(j),matmul(cart_to_prim,rgrid(:,j))
     enddo
     write(ulog,*)' sum of rws_weights(j) ',sum(rws_weights)

     call make_grid_weights_WS(gs1,gs2,gs3,g01,g02,g03,nggrid,'g',gws26,g0ws26)
     call write_lattice(size(gws26,2),gws26,'g_supercell.xyz') ! grid of gvectors of WS supercell in the WS of primitive
     write(ulog,*)' original ggrid(:,j),gws_weights(j) '
     do j=1,nggrid
        write(ulog,9)j,ggrid(:,j),gws_weights(j),matmul(transpose(prim_to_cart),ggrid(:,j))/(2*pi)
     enddo
     write(ulog,*)' sum of gws_weights(j) ',sum(gws_weights)

! check: structure factor should be 1 on the superlattice reciprocal vectors
   do i=1,nggrid
      call structure_factor_recip(ggrid(:,i),nrgrid,rgrid,rws_weights,sf,dsf)
      if( length(ggrid(:,i)) .lt. tolerance) then
         if(abs(sf-1).gt.tolerance ) then
            write(*,'(a,i3,99(1x,f9.4))')'SF=1 check for i,sf,ggrid_red(i)=',i,sf, &
&                  matmul(transpose(prim_to_cart),ggrid(:,i))/(2*pi)
            stop
         else
            cycle
         endif
      endif
      if(abs(sf).gt.1d-5 ) then
         write(*,'(a,i3,99(1x,f9.4))')'SF\=0 check for i,sf,ggrid_red(i)=',i,sf, &
&                  matmul(transpose(prim_to_cart),ggrid(:,i))/(2*pi)
!         stop
      endif
   enddo


! this tests fr2k etc
! allocate(auxr(nrgrid),auxg(nggrid))
! auxr=1
! call fr2k_5(nrgrid,auxr,rws_weights,nggrid,auxg)
! write(*,*) auxg
! call fk2r_5(nggrid,auxg,gws_weights,nrgrid,auxr)
! write(*,*) auxr
! deallocate(auxr,auxg)
! stop

  if(fc2flag.eq.0) then  ! calculate the default # of shells of FC2 based on largest supercell

     call update_nshells2(nrgrid,rgrid)
     write(ulog,*)'fc2flag=0; use largest WS; nshells2 updated to ',nshells(2,1:6)

  endif



! define the dimensions of the FC arrays, for each rank if include_fc.ne.0
! collects identical force_constants in groups of shells, identifies the irreducible ones and
! finds the mapping between the irreducible ones and all FCs in that group
  call setup_maps
! inputs: outputs of force_constants_init: iatomop etc...
! outputs: map structure including for rank, the groups of indep FCs and the full FCs
! maxterms, nshells, ngroups, nterms, estimate of inv_constraints, nindepfc

  write(ulog,*)"rank, indep FC terms , total FC terms , # of groups "
  do rank=1,maxrank
     write(ulog,*)rank,map(rank)%ntotind,map(rank)%ntot,map(rank)%ngr
  enddo

  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after setup_maps and before write_correspondance IS ',tim

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
        call read_fcs_2(iunit,fn,rank,fcrnk(rank,:),ngrnk(rank))
     endif
  endif
  enddo

! if read, may need to check compatibility with given supercells?!?
!--------------------------------------------------------------------

! before imposing the constraints, first define the number of columns in matrix A
 allocate(keep_fc2i(map(2)%ntotind))  ! which indep fc2s to keep based on the vectors rws26
! keep_fc2i(i)=1 is the list of FC2s within WS of supercell defined by 2ws26
! find the FC2s within the SC and their number, which is size_kept_fc2

  if(fc2flag.eq.0) then
     call setup_FC2_in_supercell !(keep_fc2i,size_kept_fc2)
  else  ! use the default values found by setup_maps
      keep_fc2i=1  ! keep everything
      size_kept_fc2=map(2)%ntotind !sum(keep_fc2i)  not sum of the groups but sum of all rnk2 indep terms
  endif

! calculate the new nindepfc based on keep_fc2i
! write(ulog,*)' new ngr(fc2) based size_kept_fc2=',map(2)%ngr-sum(map(2)%ntind)+size_kept_fc2
 write(ulog,*)' Old,old and new number of fc2:size_kept_fc2=',map(2)%ntotind,size_kept_fc2

! get the number of independent force constants from setup_maps
 nindepfc=0
 do rank=1,maxrank
  if ( include_fc(rank) .eq. 1 ) then ! exclude them if they already exist and will be read
! for now, we put in amat all fc2 of number ntotind
        if(rank.eq.2) then
           nindepfc = nindepfc + size_kept_fc2 ! this is the #of groups used for FC2s(not ti terms!); size in setup maps is ignored
        else
           nindepfc = nindepfc + map(rank)%ntotind
        endif
       write(ulog,*)'rank,nindep_rnk,cumulative nindep_kept=',rank,map(rank)%ntotind,nindepfc
   endif
 enddo
 write(ulog,*)' MAIN: total # of independent FCs of all rank, nindepfc=',nindepfc

! update map(2)%mat by setting the weight of non-kept fcs to zero
 do g=1,map(2)%ngr
 do ti=1,map(2)%ntind(g)
     write(ulog,*)'g,ti,counter,keep,current=',g,ti,counter2(g,ti),keep_fc2i(counter2(g,ti)),current2(g,ti)
     map(2)%gr(g)%mat(:,ti) = map(2)%gr(g)%mat(:,ti)*keep_fc2i(counter2(g,ti))
 enddo
 enddo

 allocate(frc_constr(fdfiles))  ! for each of the supercell files

 ! Read force-displacements from FORCEDISPi (OUTCARi) and construct aforce,bforce matrices
 call read_snaps_set_aforce(frc_constr,nconfigs)

  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after read_snaps_set_aforce             IS ',tim
!
!--------------------------------------------------------------------
! depending on the constraints, setup the homogeneous part of the A and B matrices:
! if(born_flag.le.2) then ! use a real space treatment of ASR and rotational invariances

    call estimate_inv_constraints

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

! now that aforce & bforce matrices are setup, we can do SVD and project on
! onto the kernel of homogeneous part of amatrix

! inputs: inv_constraints,force_constraints, aforce, bforce and invariance parts of amat
! allocates amat,bmat and fills them with atrans,arot,ahuang and aforce (same for b)
! call include_constraints_remove_zeros
! outputs are amat and bmat and their dimensions:dim_al,dim_ac, which go to SVD

! output is the homogeneous part of Amatr: ahom(inv_constraints,dim_ac)
    call homogeneous_constraints_overlap(inv_constraints) 

    write(ulog,*)'dim_al,dim_ac set to ',dim_al,dim_ac

! give warnings if any column in amat is totally zero
! allocate(amat(inv_constraints+force_constraints,dim_ac),bmat(inv_constraints+force_constraints))
!   call check_zero_column(force_constraints,dim_ac,aforce)


    call cpu_time(tim)
    write(utimes,'(a,f10.4)')' TIME after include_constraints and 0 column  IS   ',tim

! write the a and b matrices
  if (verbose) then
     write(umatrx,*)'#========= before call to svd, ahom is:',inv_constraints,dim_ac,'========='
     do i=1,inv_constraints
        write(umatrx,17)(ahom(i,j),j=1,dim_ac) !,bmat(i)
     enddo
     write(umatrx,*)'#========= before call to svd, aforce and bforce are:',dim_al,dim_ac,'========='
     do i=1,force_constraints
        write(umatrx,18)(aforce(i,j),j=1,dim_ac),bforce(i)
     enddo
  endif
17 format(299(1x,f9.4))
18 format(299(f10.6))

! do the SVD decomposition; beware amat is overwritten
  if(allocated(fcs)) deallocate(fcs)
  if (allocated(sigma)) deallocate(sigma)
  allocate(fcs(dim_ac),sigma(dim_ac))

  uio = 349
  if(enforce_inv.eq.1 .or. itemp.eq.1) then ! elimination wil be used on ahom

      open(uio,file='elimination.dat',status='unknown')
! separate homogeneous part of amat from inhomogeneous part by finding the first non-zero elt of bmat
!      call find_first_nonzero(dim_al,bmat,n_hom)
      n_hom=inv_constraints
      write(ulog,*)'MAIN: Using kernel projection: enforce_inv and itemp are=',enforce_inv,itemp
      write(ulog,*)'MAIN: position of the first non-zero elt of bmat=',n_hom
      write(   *,*)'MAIN: position of the first non-zero elt of bmat=',n_hom
!     if(allocated(a_hom)) deallocate(a_hom)
!     allocate(a_hom(n_hom,dim_ac),kernelbasis(dim_ac,dim_ac))
      allocate(kernelbasis(dim_ac,dim_ac))
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

      write(*,*) ' Temperature of ',tempk,' will be implemented'
      allocate(mat(dim_ac,dim_ac),qmat(dim_ac))
!      2023-07-03 16:24:12 ---> Checking the energy value as FOCEX code is not giving proper result (logged by Bikash. T.)
      write(*,*) "Energy values from the OUTCAR1 are: ",energies
      write(*,*) "Dimension of energies are: ",shape(energies)
      write(*,*) "The value of nlines are: ",nlines
      write(*,*) "The value of nconfigs is: ",nconfigs
!     call implement_temperature(dim_al-n_hom,dim_ac,amat(n_hom+1:dim_al,1:dim_ac), &
      call implement_temperature(force_constraints,dim_ac,aforce, &
&         bforce,nconfigs,energies,nlines,tempk,mat,qmat)
      call solve_svd(dim_ac,dim_ac,mat,qmat,svdcut,xout,sigma,uio)
      call write_out(ulog,'TEMP: SVD solution before kernel projection',xout)
      call project_on(dim_ac,nkernel,xout,kernelbasis(1:dim_ac,1:nkernel),fcs)
      deallocate(kernelbasis,mat,qmat)
  else

      if(enforce_inv.eq.1) then 
! elimination will be used; first svd for force-displacement, followed by projection on the kernal of a_hom
         write(ulog,*)'MAIN: Using svd of the force amatrix since enforce_inv and itemp are=',enforce_inv,itemp
         write(ulog,*)'MAIN: size of the homogeneous part is ',n_hom
!         call solve_svd(dim_al-n_hom,dim_ac,amat(n_hom+1:dim_al,1:dim_ac), &
!&             bmat(n_hom+1:dim_al),svdcut,bfrc,sigma,uio)
!        call svd_set(dim_al-n_hom,dim_ac,amat(n_hom+1:dim_al,1:dim_ac),  &
!&            bmat(n_hom+1:dim_al),xout,sigma,svdcut,error,ermax,sig,'svd-elim.dat')
         call svd_set(force_constraints,dim_ac,aforce,  &
 &            bforce,xout,sigma,svdcut,error,ermax,sig,'svd-elim.dat')
         call write_out(ulog,'Enforce=1: SVD solution before kernel projection',xout)
         write(ulog,*)'MAIN: performing projection on the kernel basis'

         call project_on(dim_ac,nkernel,xout,kernelbasis(1:dim_ac,1:nkernel),fcs)
         call write_out(ulog,'Enforce=1: Solution After kernel projection',fcs)
         deallocate(kernelbasis)
      else
         if (allocated(amat)) deallocate (amat)
         if (allocated(bmat)) deallocate (bmat)
         dim_al=inv_constraints+force_constraints
         allocate(amat(dim_al,dim_ac),bmat(dim_al))
         amat(1:inv_constraints,:) = ahom ; bmat(1:inv_constraints) = 0
         amat(1+inv_constraints:dim_al,:) = aforce ; bmat(1+inv_constraints:dim_al) = bforce
         deallocate (aforce,bforce)

         call svd_set(dim_al,dim_ac,amat,bmat,fcs,sigma,svdcut,error,ermax,sig,'svd-all.dat')
         write(ulog,'(a,f7.3,a)')'Percent error, || F_dft-F_fit || / || F_dft || =',sig,' %'
      endif

  endif
  close(uio)

  call write_invariance_violations(ulog,dim_ac,fcs)

  do rank=1,4
     write(xt,'(i1)')rank
     if(include_fc(rank).eq.1) then
        open(ufit1-1+rank ,file='fc'//xt//'_irr.dat',status='unknown')
        open( ufc1-1+rank ,file='fc'//xt//'.dat',status='unknown')
     endif
  enddo

  call write_independent_fcs(dim_ac,sigma,sd,ulog)

  call write_output_fcs

  call cpu_time(tim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After SVD and writing of FCs,           TIME IS ',tim

  close(ufit1)
  close(ufit2)
  close(ufit3)
  close(ufit4)
  close(ufc1)
  close(ufc2)
  close(ufc3)
  close(ufc4)

  deallocate(sigma,amat,bmat,xout) !,fcs)

! calculate the phonon dispersion as output test: requires kpbs.in as input
  call set_band_structure !(uband,ugrun)

  call cpu_time(tim)  !date_and_time(time=tim)
  write(utimes,'(a,f12.4)')' After phonon band structure,            TIME IS ',tim


  call date_and_time(date=today,time=now,zone=zone)
  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' ENDING TIME OF THE PROGRAM                   IS ',tim
  write(ulog,'(a,f10.4)')' ENDING TIME OF THE PROGRAM                   IS ',tim


2 format(a,f7.3,a)
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
        call write_supercell

        call make_r0g

! make sure of consistency POSCAR with input structure
        call identify_atoms_in_supercell

! writing the coordinates of primitive cells in the supercell
! reference of each prim cell (n) is atom #1
        open(111,file='Bravaissupercell-'//xt//'.xyz')
        write(111,*) natom_super_cell+4  ! also draw the 3 translation vectors
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

! generating rgrid and ggrid for Ewald force calculations in the supercell
! for each different FORCEDISPi file the supercell(rs1,rs2,rs3) is different, so
! the rgrid and ggrid need to be generated every time for subtracting Coulomb forces
        if(born_flag.ne.0) then
            call make_grid_weights_WS(r01,r02,r03,rs1,rs2,rs3,nrgrid,'r',r0ws26,rws26) !,rws_weights)
            write(ulog,*)' sum of rws_weights in SC# ',i,sum(rws_weights)
            call make_grid_weights_WS(gs1,gs2,gs3,g01,g02,g03,nggrid,'g',gws26,g0ws26) !,gws_weights)
            write(ulog,*)' sum of gws_weights in SC# ',i,sum(gws_weights)
        endif

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

! 3 ways: if=0 ignore, if=1 just all if=2 use real space ewald, if=3 or 4  use Fourier transforms
       write(ulog,*) 'last Forces BEFORE subtraction'
       write(*,*) 'last Forces BEFORE subtraction'
       call write_out(ulog,'Last Force',force(:,:,ncfg))

       call subtract_coulomb_force(born_flag,ncfg,displ,force)
! now frc contains the short-range (Non-Coulomb) part of the forces which we can fit

       write(ulog,*) 'last Forces AFTER subtraction'
       write(*,*) 'last Forces AFTER subtraction'
       call write_out(ulog,'Last Force',force(:,:,ncfg))

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


end program FC234


