!=====================================================

program FOCEX
!
! todo:
! add option of reading existing fc files and fitting the remaining ones.
! also writeout the fitted FCs in the format of other codes: ALAMODE, PHONOpy SHENG
!
! test ewald dynmat versus ewald force by finite difference? 
!
! add CS for fc3 fc4 so that the 2-body terms remain; or select the same 2-body terms as
! fc2's, for example as analytical derivatives of the morse potential fitted to the VASP data when
! only one atom is moved at least 4 times +/- 0.1, +/- 0.05 (if not 6: +/- 0.15) plus a short-range term
! eventually may need to replace include_fc=2 by reading the fitted FC234's from an existing file
! and completing by the missing ones.
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
use tetrahedron

implicit none
integer i,j,g,ti,rank,iunit, uio,nkernel,imax,n_hom,nt(maxrank),ntind(maxrank),tau,taup,nat,ncs(3)
character xt*1,fn*7,now*10,today*8,zone*5,fni*11
logical ex
real(r15), allocatable :: xout(:),foldedk(:,:),rand(:),rand2(:)
complex(r15), allocatable :: auxr(:),auxg(:),axavg(:),ax1(:),ax2(:)
real(r15), allocatable :: mat(:,:),qmat(:)
integer, allocatable :: frc_constr(:),mp(:) !,save_boundary(:)
real(r15) error,ermax,sd(4),sig, volmax,sf,dsf(3) ,q2(3),qr(3),largestcutoff,matr(3,3)
 type(vector):: y1,y2,y3,k1,k2,k3
real tim
real(r15) gmax,ene,grad
! real(r15) dyn_coul(3,3) ,ddn(3,3,3),q(3),deteps3 
!real(r15) r1(3),r2(3),r3(3),z1(3),z2(3),z3(3),rsh(3,26),r26(3,26),r026(3,26)

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(tim)

! tolerance = 0.001
! r1=(/1d0,1.732d0,0d0/)
! r2=(/2d0,0d0,0d0/)
! r3=(/0d0,0d0,1.5d0/)
!prim_to_cart(:,1)=r1
!prim_to_cart(:,2)=r2
!prim_to_cart(:,3)=r3
!call xmatinv(3,prim_to_cart,cart_to_prim,i)
! r01=a2v(r1)
! r02=a2v(r2)
! r03=a2v(r3)
! z1=(/4d0 ,0d0,0d0/)
! z2=(/0d0,5.196d0,0d0/)
! z3=(/0d0,0d0,3d0/)
! k1=a2v(z1)
! k2=a2v(z2)
! k3=a2v(z3)
! ngrd=500
! allocate(weig(ngrd),grd(3,ngrd))
! call get_26shortest_shell(r01,r02,r03,r026,x01,x02,x03)
! call get_26shortest_shell(k1,k2,k3,r26,x1,x2,x3)
! call make_grid_weights_WS(r01,r02,r03,k1,k2,k3,cart_to_prim,ngrd,grd,weig,r026,r26) 
! grd=grd(:,1:ngrd); weig = weig(1:ngrd)
! call write_out(ulog,'generated grid inside supercell ',transpose(grd))
! call show_ws_boundary(v2a(r01),v2a(r02),v2a(r03),r026,15,'WSR0_boundary.xyz',lgridmax) 
!    
! nat=6 ;ti=0
!   do i=-nat,nat
!   do j=-nat,nat
!   do g=-nat,nat
!      ti=ti+1
!      q2=i*r01+j*r02+g*r03 
!!     call fold_in_ws(q2,r26,qr)
!      qr = fold_fbz_ws(q2,r26)
!      write(123,3)ti,q2,qr
!   enddo
!   enddo
!   enddo
!
!   stop
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
!
!@@@@@@@@@@@@@@@@@@@@@@  read inputs and set up cell information  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! read from structure.params the atoms in prim cell in reduced units, assign their mass, type, tau;
! also read specs of the prim cell, and flags for range of FC3 and FC4, enforcement of invariance relations
! outputs: lattparams,primlat,nshells,include_fc,itrans,irot,natom_prim_cell,mas,atname,
! atompos0 (reduced coords within the primitive cell) and object atom0(natom_prim_cell)
   call read_structure


! reads dielectric constant and Born charges (=0 if unknown)
   call read_dielectric

! the number of shells is based on nshells in the structure.params file and should be large enough
! to include all atoms in the supercell
! inputs: natom_prim_cell,latticeparameters,primitivelattice,atom_type,atompos0
! outputs: r0i, atompos (cart coord of atoms and natoms neighbor shells),
! iatomcell, iatomcell0, iatomneighbor, iatomop, ... and symmetry operations
! sets nsmax, according to the default range for fc2s, rcut(2)
   call make_unitcell

 coef=ee*1d10/eps0/volume_r0 ! used for non-analytical part of the dynamical matrix

   call find_WS_largest_SC(imax,volmax)

   call make_grids

! check if any include_fc=2, then read atompos and call setup_maps and read fcs
   if(readfc)  call read_fcs  

   call write_atompos  !  atompos in .xyz format for visualization

!  call find_WS_largest_SC(imax,volmax)
!  call write_out(6,' rsc1 ',rs1)
!  call write_out(6,' gsc1 ',gs1)
!  call make_grid_weights_WS(r01,r02,r03,rs1,rs2,rs3,nrgrid,'r',r0ws26,rws26)
!  call make_grid_weights_WS(gs1,gs2,gs3,g01,g02,g03,nggrid,'g',gws26,g0ws26)
!4 format(a,3(1x,f9.4),1x,g13.6)
     
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

!  call read_latdyn
!  allocate(kpc(3,nkc),wk(nkc),mapinv(nkc),mapibz(nkc))
!  call allocate_tetra(nc(1),nc(2),nc(3),100)
!  call make_kp_reg_tet
!  call get_weights(nkc,kpc,mapinv,mapibz) 

! of available supercells, takes the one with largest volume to eventually set
! the range of FC2= largest center to WS boundary distance; record rs1,rs2,rs3
! output is rs1,rs2,rs3 of this largest supercell

 call cpu_time(tim)
 write(utimes,'(a,f10.4)')' make_unitcell, TIME                          IS ',tim

! outputs: atom0%equilb_pos,shells%no_of_neighbors , rij, neighbors(:)%tau and n and maxshell
! maxshell is the actual # of shells within rcut(2). It is not nsmax; maxshell < maxneighbors
 if(.not.readfc)  call set_neighbor_list

 write(ulog,*)' After set_neighbors: maxshells=',maxshells,' while maxneighbors=',maxneighbors
 write(ulog,*)' ************************* Writing the neighborlist ************************** '

  call write_neighbors


! along with their weights  in case they are on the WS boundary
! rws26(output) are the 26 shortest vectors used to define the WS of the largest supercell
! in addition depending on 'r' or 'g' the grid of primitive translations inside the WS of supercell
! is calculated and stored in rgrid(3,nrgrid) or ggrid(3,nggrid) along with their weights rws_weights
! create a grid of translations from r0i inside the WS cell of rsi +boundaries, rgrid(3,nrgrid),

! find the shortest basis bi among xi and the corresponding 26 shortest vectors
! call get_26shortest_shell(x01,x02,x03,s0x,b01,b02,b03)! s0x for WS of primcell
! call get_26shortest_shell(x1 ,x2 ,x3 ,sx ,b1 ,b2 ,b3) ! sx for WS of supercell
! x01=b01; x02=b02; x03=b03; x1=b1;  x2=b2;  x3=b3; 

! if(space.eq.'r' .or. space.eq.'R') then
!    matr=cart_to_prim
! elseif(space.eq.'g' .or. space.eq.'G') then
!    matr=transpose(prim_to_cart)/(2*pi)
! endif
!
! call make_grid_weights_WS(x01,x02,x03,x1,x2,x3,matr,ngrd,grd,weig,s0x,sx) 


!! generate primitive grid in SC
! n3reg(1)=4; n3reg(2)=4; n3reg(3)=2 ;
! nreg=nint(volume_r/volume_r0)
! if(nreg .ne. n3reg(1)*n3reg(2)*n3reg(3)) then
!    nreg = n3reg(1)*n3reg(2)*n3reg(3)
! endif
! write(ulog,*)'MAIN: there are ',nreg,' primitive vectors inside supercell'
! allocate(rgridreg(3,nreg),ggridreg(3,nreg))

! this is used to compare to the output of generate_grid_no_image
! call generate_grid_sc(nreg,rgridreg) !,ggridreg)
! call write_out(ulog,'# Primitive rgrid vectors ',transpose(rgridreg)) 

! allocate(mp(nrgrid))
! call generate_grid_noimage(nrgrid,rgrid,rws26,nat,mp)
! write(ulog,*)'NEWN(R)=',nat
!  do i=1,nat !nrgrid
!     call bring_to_ws_gx(rgridreg(:,i),rws26,rgridreg(:,i))
!     call write_out(ulog,'#rgridX vectors ',rgridreg(:,i))
!  enddo
!  deallocate(mp)

!  allocate(mp(nggrid))
! allocate(xout(nreg),rand(nreg),rand2(nreg),regweights(nreg),rgridreg(3,nreg),ggridreg(3,nreg)) !, rws_weights(nreg),gws_weights(nreg)
! q2=0 ! shift=0
! call generate_grid_noimage(nggrid,ggrid,g0ws26,nat,mp)
! call make_grid(n3reg,r01,r02,r03,q2,rgridreg,regweights)
! call make_grid(n3reg,gs1,gs2,gs3,q2,ggridreg,regweights)
! write(ulog,*)'NEWN(G)=',nat
! call write_out(ulog,'# Primitive ggrid vectors ',transpose(ggridreg))
!  do i=1,nat !ncs(1)*ncs(2)*ncs(3)  !nat !nggrid
!     call bring_to_ws_gx(ggridreg(:,i),g0ws26,ggridreg(:,i))
!     call write_out(ulog,'#ggridX vectors ',ggridreg(:,i))
!  enddo
! deallocate(rand)


!! TEST OF FOURIER!!  first make sure the grid is periodic
!     allocate(ax1(nrgrid),rand(2*nggrid),rand2(2*nggrid),auxr(nggrid),axavg(nggrid),ax2(nggrid))
!     call random_number(rand)
!     auxr=cmplx(rand(1:nggrid),rand(nggrid+1:2*nggrid))
!     call check_periodic(nggrid,ggrid,g0ws26,auxr,axavg,ex) ! rand 2 is periodic version of rand
!     call fourier_k2r(axavg,ax1) 
!     call fourier_r2k(ax1,ax2)  
!     write(*,*)' axavg before is '
!     write(*,77) axavg
!     write(*,*)' axavg after is '
!     write(*,77) ax2
!     write(*,*)'|difference|=', (maxval(abs(ax2-axavg))) 
!     stop
!! fourier on regular grid
!!call fr2k_5(nreg,rgridreg,rand,regweights,nreg,ggridreg,ax1)
!!call fk2r_5(nreg,ggridreg,ax1,regweights,nreg,rgridreg,rand2)
!    write(*,*)' rws_weights '
!    write(*,77) rws_weights
!    write(*,*)' gws_weights '
!    write(*,77) gws_weights
!    do i=1,0 !nrgrid
!       ax2=0
!       ax2(i)=1
!       if(length(rgrid(:,i)) .gt. tolerance ) call is_on_boundary(rgrid(:,i),rws26,j,tolerance)
!! fourier on WS grid
!     call fourier_r2k(ax2,ax1) 
!        write(*,75)' ax1 after FFT r2k is ',i,rgrid(:,i),j
!        write(*,77) ax1
!     call fourier_k2r(ax1,ax2)  
!        write(*,*)' ax2 after FFT k2r is '
!        write(*,77) ax2
!     enddo
!     deallocate(ax1,rand,rand2,axavg,auxr,ax2)
!
!!   stop
! 75 format(a,i4,3(f10.6),i9.2)

! TEST of EWALD dynamical matrix
!  call allocate_fc_dyn(natom_prim_cell,nrgrid,nggrid)
!  call calculate_dyn_na(nggrid,ggrid,dyn_na)   ! hard phase with ASR enforced
!  stop




 if(.not.readfc) then

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
! full symmetry of the primitive cell for NA subtraction in case fc2range.ne.0
! output is nsubgrid, subgrid,subgrid_weights
   if(fc2range.ne.0) call make_subrgrid_symmetric(nrgrid,rgrid,rws26,nsubgrid)  

!q2=0
!do i=0,25
!   q2=(i/100d0)*v2a(g01+g02)
!       call structure_factor_recip(q2,nsubgrid,subgrid,subgrid_weights,sf,dsf)
!       write(*,'(a,i3,99(1x,f9.4))')'SF check for i,sf,dsf,ggrid_red(i)=',i,q2,sf,dsf &
! &              ,matmul(transpose(prim_to_cart),q2)/(2*pi)
!       call structure_factor_recip(q2,nrgrid,rgrid,rws_weights,sf,dsf)
!       write(*,'(a,i3,99(1x,f9.4))')'SF check for i,sf,dsf,ggrid_red(i)=',i,q2,sf,dsf &
! &              ,matmul(transpose(prim_to_cart),q2)/(2*pi)
!enddo
! stop
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
! find the FC2s within the SC and their number, which is map(2)%nkeptind, 
!! by K1 on July 21 2024
! complete by symmetry, but remove anything that goes beyond WS cell of the superell
!! by K1 on July 21 2024

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
  write(utimes,'(a,f10.4)')' TIME after read_snaps_set_aforce             IS ',tim

  close(ucor)
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
        write(umatrx,77)(ahom(i,j),j=1,dim_ac) !,bmat(i)
     enddo
     write(umatrx,*)'#========= before call to svd, aforce and bforce are:',dim_al,dim_ac,'========='
     do i=1,force_constraints
        write(umatrx,88)(aforce(i,j),j=1,dim_ac),bforce(i)
     enddo
  endif

77 format(299(1x,f9.4))
88 format(299(f10.6))

!
!@@@@@@@@@@@@@@@@@@@@@@  Solve the linear equations using SVD and write the solution  @@@@@@@@@@@@@@@@@@@@@@@@@@
!
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
      if(allocated(qmat)) deallocate(qmat)
      if(allocated( mat)) deallocate( mat)
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
         write(ulog,'(a,f7.3,a)')'Percent error, || F_dft-F_fit || / || F_dft || =',sig,' %'
         call write_out(ulog,'Enforce=1: SVD solution before kernel projection',xout)
         write(ulog,*)'MAIN: Invariance violations '
         do i=1,n_hom
            write(ulog,*)i,dot_product(ahom(i,:),xout(:))
         enddo

         write(ulog,*)'MAIN: performing projection on the kernel basis'
         call project_on(dim_ac,nkernel,xout,kernelbasis(1:dim_ac,1:nkernel),fcs)
         call write_out(ulog,'Enforce=1: Solution After kernel projection',fcs)
         write(ulog,7)'Max and Average Change After kernel projection',maxval(abs(fcs-xout)), &
&                      sqrt(dot_product(fcs-xout,fcs-xout))/dim_ac
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
         write(ulog,'(a,f7.3,a)')'Percent error, || F_dft-F_fit || / || F_dft || =',sig,' %'
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

  call write_independent_fcs(dim_ac,sigma,sd,ulog)

  call write_output_fcs

  endif 

  call write_fc2matrix
! call write_fc2matrix_new

! Fourier transform fc2 at G_sc vectors and subtract the DYN_NA(G) and Fourier back  
! to get the short-range part
  call correct_fcs

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

  deallocate(sigma) 

!  call energy_strain(0.01,ene,grad) 
!  write(ulog,*)'for strain =0.01, ene,grad=',ene,grad
!
!@@@@@@@@@@@@@@@@@@@@@@  Calculate phonons and thermodynamic properties  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
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

 write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' Program  FOCEX ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' Program FOCEX ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

  write(ulog,7)' Average and largest errors in force=',error,ermax
  write(ulog,7)' Average errors in FOCEX=',sd

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

! generating rgrid and ggrid for Ewald force calculations in the supercell
! for each different FORCEDISPi file the supercell(rs1,rs2,rs3) is different, so
! the rgrid and ggrid need to be generated every time for subtracting Coulomb forces
!       if(born_flag.ne.0) then
!           call make_grid_weights_WS_old(r01,r02,r03,rs1,rs2,rs3,nrgrid,'r',r0ws26,rws26) !,rws_weights)
!           write(ulog,*)' sum of rws_weights in SC# ',i,sum(rws_weights)
!           call make_grid_weights_WS_old(gs1,gs2,gs3,g01,g02,g03,nggrid,'g',gws26,g0ws26) !,gws_weights)
!           write(ulog,*)' sum of gws_weights in SC# ',i,sum(gws_weights)
!       endif

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
 !     write(ulog,*) 'last Forces BEFORE subtraction'
 !     write(*,*) 'last Forces BEFORE subtraction'
 !     call write_out(ulog,'Last Force',transpose(force(:,:,ncfg)))
 !     if(verbose) call write_forces_display(natom_super_cell,ncfg,displ,force,'beforesub')
 !
 !     call subtract_coulomb_force(born_flag,ncfg,displ,force)
! now frc contains the short-range (Non-Coulomb) part of the forces which we can fit
 !
 !      write(ulog,*) 'last Forces AFTER subtraction'
 !      write(*,*) 'last Forces AFTER subtraction'
 !      call write_out(ulog,'Last Force',transpose(force(:,:,ncfg)))
       if(verbose) call write_forces_display(natom_super_cell,ncfg,displ,force,'aftersub')

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

 subroutine write_fc2matrix
!! takes the fitted FCs, writes into periodic 5D matrix form by fourier transforming back the dynamical matrix
 use svd_stuff
 use ios
 use atoms_force_constants
 use params
 use lattice
 use constants
 use fourier
 use born
 use geometry, only : length
 implicit none
 integer tau,taup,al,be,g,t,ti,ired,j,igrid,nr(3),nboundary,ndn,ixy(2),cnt,uio,l
 complex(r15) auxr(nrgrid),auxg(nggrid),aravg(nrgrid),agavg(nggrid),phase(nggrid)
 real(r15) rr(3),rfold(3),dta(3),fc2,asr,rij(3),error
 logical isperiodic,insid

 ndn = 3*natom_prim_cell
 call allocate_fc_dyn(natom_prim_cell,nrgrid,nggrid)
 uio=345
 open(uio,file='trace_fc2_lr.dat')
 write(uio,*)'# tau,taup,j,|r+taup-tau|,trace(fc2),igrid,r+taup-tau,red_fold(r+taup-tau),fc2(al,be)'

! if( born_flag.ne.0 ) call calculate_dyn_na(nggrid,ggrid,dyn_na)
! outputs are dyn_na(n0,n0,3,3,ng) and dyn_naq0(n0,3,3)

! use the bare fc2's from fitting to calculate the dynamical matrix (without masses)
 do g=1,nggrid
    call set_dynamical_matrix4d(ggrid(:,g),dyn_g(:,:,:,:,g),ndn)
 enddo

 phib=0
 write(ulog,*)'WRITE_FC2MATRIX: setting up bare fc2 in 5D form '
 do tau =1,natom_prim_cell
 do taup=1,natom_prim_cell
    dta = atompos(:,taup) - atompos(:,tau) 
 !  jloop: do j=1,natoms
 !     if(iatomcell0(j).ne.taup) cycle
 !     rij = atompos(:,j)-atompos(:,tau)
 !     call check_inside_ws(rij,rws26,insid,nboundary) 
 !     if (.not. insid) then
 !        cycle
 !     else
         cnt=0
         gridloop: do igrid=1,nrgrid
            rr=rgrid(:,igrid)+dta 
! added on Jan 14th
                 call check_inside_ws(rr ,rws26,insid,nboundary) 
                 if (.not. insid) then
                    rfold = fold_ws(rr,rws26,'r')  
                 else
                    rfold=rr
                 endif
 !          if ( rr .myeq. rij) then ! this rij is now inside
 !             rfold=rr
 !          else
 !             rfold = fold_ws(rr,rws26,'r')  
 !             if (.not. ( rfold .myeq. rij)) cycle gridloop
 !             write(ulog,9)'R+taup-tau not inside WS for tau,taup,igrid,j,rij,rfold= ',tau,taup,igrid,j,rij,rfold
 !          endif
            cnt=cnt+1  ! to make sure only one grid corresponds to a j 
            do g=1,map(2)%ngr
            do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
               ired = map(1)%nkeptind + current(2,g,ti)
               do t=1,map(2)%nt(g)
 !                if ( tau .ne. map(2)%gr(g)%iat(1,t) .or. j .ne. map(2)%gr(g)%iat(2,t) ) cycle 
                      if ( tau .ne. map(2)%gr(g)%iat(1,t) ) cycle 
                      j= map(2)%gr(g)%iat(2,t) 
                      if(iatomcell0(j).ne.taup) cycle
                      rij = atompos(:,j)-atompos(:,tau)
                      if(.not. (rij.myeq.rfold)) cycle

                  ixy(:) = map(2)%gr(g)%ixyz(:,t)
                  phib(tau,taup,ixy(1),ixy(2),igrid) = phib(tau,taup,ixy(1),ixy(2),igrid) + &
 &                          fcs(ired) * map(2)%gr(g)%mat(t,ti) 
               enddo
            enddo
            enddo

            fc2 = trace(real(phib(tau,taup,:,:,igrid))) ! this should satisfy ASR!
!         write(uio,3)tau,taup,j,length(rij),fc2,real(igrid),rij,cart2red(rij,'r'),real(phib(tau,taup,:,:,igrid)) 
         enddo gridloop
         if(cnt.ge.2 .or. cnt.eq.0) write(*   ,1)'ERROR: cnt FOR t,tp,j,igrid =',cnt,tau,taup,j,igrid
         if(cnt.ge.2 .or. cnt.eq.0) write(ulog,1)'ERROR: cnt FOR t,tp,j,igrid =',cnt,tau,taup,j,igrid
 !     endif
 !  enddo jloop

! now use Fourier transformation to get the periodic one from dyn_g
    do al=1,3
    do be=1,3
    ! the bare phib is not periodic; no need to check!
    !  call check_periodic(nrgrid,rgrid,rws26,phib(tau,taup,al,be,:),aravg,isperiodic)
!      call check_periodic(nggrid,ggrid,g0ws26,dyn_g(tau,taup,al,be,:),agavg,isperiodic)
!      if(.not.isperiodic) call warn4(ulog,'dyn_g was not periodic ',real(dyn_g(tau,taup,al,be,:)))
    !  write(*,*)'calling fourier_r2k ',tau,taup,al,be
    !  call fourier_r2k(aravg,agavg) ! use the periodic one anyways; it is the same as dynmat4d? 

! is it the same as dyn_g from dynamt4D? 
    !  error=maxval(abs(dyn_g(tau,taup,al,be,:)-agavg(:)))  
    !  if(error.gt.1d-5) then
    !     write(ulog,9)'tau,taup,al,be,error in dynmat0:dyn_g-agavg=',tau,taup,al,be,error
    !  endif

       write(*,*)'calling fourier_k2r to get periodic phi not the bare phi'
       call fourier_k2r(dyn_g(tau,taup,al,be,:),phip(tau,taup,al,be,:))
       phip(tau,taup,al,be,:)= phip(tau,taup,al,be,:) *rws_weights(:)  ! this weights is the correct fc2 set

    !  error = maxval(abs(phip(tau,taup,al,be,:)-aravg(:)))
    !  if(error.gt.1d-5) then
    !     write(ulog,9)'t,tp,a,b:FC error |phip-aravg| :',tau,taup,al,be,error
    !     write(ulog,77)real(phip(tau,taup,al,be,:)-aravg(:))
    !  endif

    enddo
    enddo

 !  do j=1,natoms
 !     if(iatomcell0(j).ne.taup) cycle
 !     rij = atompos(:,j)-atompos(:,tau)
       do igrid=1,nrgrid
          rr=rgrid(:,igrid)+dta 
 !        if (.not. (rr .myeq. rij)) cycle  ! this effectively folds pairs outside WS to inside
          fc2 = trace(real(phip(tau,taup,:,:,igrid))) ! this should satisfy ASR!
          rfold = fold_ws(rr,rws26,'r')  ! this should not happen!!
          write(uio,3)tau,taup,igrid,length(rr),fc2,length(rfold),cart2red(rr,'r'),real(phip(tau,taup,:,:,igrid)) 
       enddo
 enddo
 enddo
 close(uio)

 write(*,*)'ASR CHECK for long-range fc2_p matrix: tau,al,be,asr'
 write(ulog,*)'ASR CHECK for long-range fc2_p matrix: tau,al,be,asr'

 do tau =1,natom_prim_cell
 do al=1,3
 do be=1,3
    write(ulog,3)tau,al,be,sum(phip(tau,:,al,be,:))
 enddo
 enddo
 enddo

1 format(a,99(1x,i4))
2 format(a,99(f9.4))
3 format(3i5,99(f9.4))
5 format(a,5i5,9(f8.3))
9 format(a,4i4,9(f9.4))
77 format(299(1x,f9.4))

 write(*,*)'EXITING write_fc2matrix'
 write(ulog,*)'EXITING write_fc2matrix'

 end subroutine write_fc2matrix

!-------------------------------------------------------

 subroutine write_fc2matrix_new
!! takes the fitted FCs, writes output phib into 5D matrix form
 use svd_stuff
 use ios
 use atoms_force_constants
 use params
 use lattice
 use constants
 use fourier
 use born
 use geometry, only : length !,myeq
 implicit none
 integer tau,taup,al,be,g,t,ti,ired,j,k,l,igrid,nboundary,uio,cnt !,cnt(natom_prim_cell)
 complex(r15) auxr(nrgrid),auxg(nggrid),aravg(nrgrid),agavg(nggrid),phase(nggrid)
 real(r15) rfold(3),dta(3),rg(3),fc2,asr,zero,rij(3)
 real(r15),allocatable :: phi2(:,:,:,:),rr(:,:)
 integer,allocatable :: t2(:),igrd(:)
 logical isperiodic,insid
! complex(r15) dyn2(natom_prim_cell,natom_prim_cell,3,3,nggrid)
! complex(r15) ddyn2(natom_prim_cell,natom_prim_cell,3,3,nrgrid,3)

 zero=0
 call allocate_fc_dyn(natom_prim_cell,nrgrid,nggrid)

 write(ulog,*)'WRITE_FC2MATRIX_NEW: just allocated arrays of sizes natom0,nr,ng',natom_prim_cell,nrgrid,nggrid
 
! use bare fc2s from fitting to calculate the dynamical matrix WITHOUT MASSES on 
! the supercell ggrid; it should be G-periodic with the old phase (step) convention; 
! we then fourier interpolate with the fourier transform of this matrix, called phib
! phib is not directly equal to the fitted fc2! and must include the weights
 do g=1,nggrid
    call set_dynamical_matrix4D(ggrid(:,g),dyn_g(:,:,:,:,g),3*natom_prim_cell) ! ,ddyn2(:,:,:,:,g,:))
 enddo

! as a check, Fourier back to get fc2 in 5D-matrix form; print the trace vs distance
! this is extended over full WS cell of the supercell!
 write(ulog,*)'WRITE_FC2MATRIX: Checking periodicity of dyn_g before FT back'
 do tau =1,natom_prim_cell
 do taup=1,natom_prim_cell
 do al=1,3
 do be=1,3
    call check_periodic(nggrid,ggrid,g0ws26,dyn_g(tau,taup,al,be,:),agavg,isperiodic)
    if(.not.isperiodic) then  ! dyn_g should be periodic!!
       call warn4(ulog,'DYN_G not periodic, Re(D)=',real(dyn_g(tau,taup,al,be,:)))
       call write_out(ulog,'IM(DYN_G)            ',aimag(dyn_g(tau,taup,al,be,:)))
       stop
    endif
!   call fourier_k2r(dyn_g(tau,taup,al,be,:),phib(tau,taup,al,be,:)) 
    call fourier_k2r(agavg                  ,phip(tau,taup,al,be,:))  ! use the periodic one anyways
    call fourier_r2k(phip(tau,taup,al,be,:),dyn_g(tau,taup,al,be,:)) 
! the above is SC-periodic and includes all images, while the below is the same as the bare fc2
    if(maxval(abs(dyn_g(tau,taup,al,be,:)-agavg)) .gt. 1d-5) then
        call warn4(ulog,'DYN_G ne agavg =',real(dyn_g(tau,taup,al,be,:)-agavg(:)))
    endif
    phip(tau,taup,al,be,:)= phip(tau,taup,al,be,:) *rws_weights(:)  ! with this weights is the correct fc2 set
 enddo
 enddo
 enddo
 enddo

! should recover old dyn_g
! fc2=maxval(abs(dyn_g-dyn2))
! if(fc2.gt.1d-4) call warn3(ulog,'1:Dyn_g and dyn2 not the same(not periodic)!!!',fc2)

! phi2=0;
! dyn_g=0;phip=0
! do tau =1,natom_prim_cell
! do taup=1,natom_prim_cell
! do igrid=1,nrgrid 
!    do g=1,map(2)%ngr
!       tloop: do t=1,map(2)%nt(g)
!          if ( tau .ne. map(2)%gr(g)%iat(1,t) ) cycle tloop
!          j  =          map(2)%gr(g)%iat(2,t)
!          if (taup .ne. iatomcell0(j)) cycle tloop 
!          if(length(atompos(:,j)-atompos(:,taup)-rgrid(:,igrid)).gt.tolerance) cycle tloop
! !        if(length(fold_ws(atompos(:,j)-atompos(:,taup),rws26,'r')-rgrid(:,igrid)) &
!! &            .gt.tolerance) cycle tloop
!!         cnt(tau)=cnt(tau)+1
!!         rr(:,cnt(tau)) = atompos(:,j) - atompos(:,tau) ! this is just Rgrid+taup-tau 
!!         rg= atompos(:,j) - atompos(:,taup) ! this is just Rgrid
!
!! find igrid corresponding to term t
!      !   call findgrid(rg,rgrid,nrgrid,igrid)
!      !   if(igrid.eq.0) then
!      !      rfold=fold_ws(rg,rws26,'r')
!      !      call findgrid(rfold,rgrid,nrgrid,igrid)
!      !   endif
!
!!     !   t2  (cnt(tau))=taup
!!     !   igrd(cnt(tau))=igrid 
!          al = map(2)%gr(g)%ixyz(1,t)
!          be = map(2)%gr(g)%ixyz(2,t)
!          do ti=1,map(2)%ntind(g)  ! index of independent terms in that group g
!             ired = map(1)%ntotind + current(2,g,ti)
!!            phi2(tau,cnt(tau),al,be) = phi2(tau,cnt(tau),al,be) + fcs(ired) * map(2)%gr(g)%mat(t,ti) 
!             phip(tau,taup,al,be,igrid) = phip(tau,taup,al,be,igrid) + fcs(ired) * map(2)%gr(g)%mat(t,ti) 
!          enddo
!       enddo tloop
!    enddo
! enddo
! enddo
! enddo

 uio=345
 open(uio,file='trace_fc2_lr.dat')
 write(uio,*)'# tau,taup,j,|r+taup-tau|,trace(fc2),igrid,red_fold(r+taup-tau),fc2(al,be)'
 do tau =1,natom_prim_cell
 do taup=1,natom_prim_cell
     dta = atompos(:,taup) - atompos(:,tau) 
     do j=1,natoms
     if(iatomcell0(j).ne.taup) cycle
     rij = atompos(:,j)-atompos(:,tau)
     do igrid=1,nrgrid
        rg=rgrid(:,igrid)+dta 
        if (.not. (rg .myeq. rij)) cycle 
        call check_inside_ws(rg,rws26,insid,nboundary) 
        if(.not. insid) then
           rfold = fold_ws(rg,rws26,'r') 
        else
           rfold=rg
        endif
        fc2 = trace(real(phip(tau,taup,:,:,igrid))) ! this should satisfy ASR!
        if(abs(fc2).lt.0.0001) cycle
!   if(length(rg-rfold) .myeq. zero ) then
        write(uio,3)tau,taup,j,length(rfold),fc2,real(igrid)  ,cart2red(rfold,'r'),real(phip(tau,taup,:,:,igrid)) 
!   else
!      write(uio,3)tau,taup,igrid,length(rg),fc2,length(rfold),cart2red(rg,'r'),real(phip(tau,taup,:,:,igrid)) 
!   endif
 enddo
 enddo
 enddo
 enddo
 close(uio)

 write(ulog,*)'ASR CHECK for long-range fc2 matrix: tau,al,be,asr'
 do tau =1,natom_prim_cell
 do al=1,3
 do be=1,3
 asr=0
 do igrid=1,nrgrid
    asr=asr+sum(phip(tau,:,al,be,igrid))  ! this should satisfy ASR!
 enddo
    write(ulog,3)tau,al,be,asr
 enddo
 enddo
 enddo


2 format(a,99(f9.4))
3 format(3i5,99(f9.4))
4 format(a,4i5,9(f8.3))
5 format(a,5i5,9(f8.3))

 end subroutine write_fc2matrix_new

!-------------------------------------------------------

 subroutine correct_fcs
!! takes the dyn(G) calculated at supercell reciprocal vectors
!! subtracts dyn^NA, fourier back to update phib used later for fourier interpolation of the dynamical matrix
!! works correctly if fc2range=0 full SC is sampled; otherwise DYN_NA should be calculated 
!! on a subrid conistent with nshell that has full symmetry of the crystal
 use ios
 use atoms_force_constants
 use params
 use lattice
 use constants
 use fourier
 use born
 use geometry, only : length
 implicit none
 integer tau,taup,al,be,g,t,ti,ired,j,k,l,igrid,nr(3),nboundary,ia,ib,ndyn
 complex(r15) auxr(nrgrid),auxg(nggrid),aravg(nrgrid),agavg(nggrid),phase(nggrid)
 real(r15) rr(3),rfold(3),dta(3),mi,mj,fc2,asr,rij(3)
 logical isperiodic,insid
 complex(r15) dynmat(3*natom_prim_cell,3*natom_prim_cell) ,ddyn(3*natom_prim_cell,3*natom_prim_cell,3)
 complex(r15) dyn2(natom_prim_cell,natom_prim_cell,3,3,nggrid)
! complex(r15) phi(natom_prim_cell,natom_prim_cell,3,3,nrgrid)

 ndyn=3*natom_prim_cell

 if( born_flag.eq.0 ) return 

 call calculate_dyn_na(nggrid,ggrid,dyn_na)   ! hard phase with ASR enforced ; no mass factor
! outputs are dyn_na(n0,n0,3,3,ng) on ggrid and dyn_naq0(n0,3,3)
! can only subtract the g_ewald terms since R_ewald terms are short-ranged and do not really affect the NA part 

 write(ulog,*)'LARGEST DYN_NA=',maxval(abs(dyn_na))

! subtract the NA term
  dyn_g = dyn_g - dyn_na

! fix ASR by subtracting dyn_naq0 from diagonal <= already done in calculate_dyn_na
! do g=1,nggrid
! do tau=1,natom_prim_cell
!    dyn_g(tau,tau,:,:,g)=dyn_g(tau,tau,:,:,g)+dyn_naq0(tau,:,:)
! enddo
! enddo

! fourier transform back to make it short-range
 do tau =1,natom_prim_cell
 do taup=1,natom_prim_cell
 do al=1,3
 do be=1,3
    call fourier_k2r(dyn_g(tau,taup,al,be,:),phip(tau,taup,al,be,:)) 
    phip(tau,taup,al,be,:)= phip(tau,taup,al,be,:) *rws_weights(:)  
! with this weights it's the correct fc2 set (periodic not bare)
 enddo
 enddo
 enddo
 enddo

  write(ulog,*)'ASR CHECK for short-ranged FC2: tau,al,be,asr'
  do tau =1,natom_prim_cell
  do al=1,3
  do be=1,3
    asr=0
    do igrid=1,nrgrid
       asr=asr+sum(phip(tau,:,al,be,igrid)) 
    enddo
    write(ulog,3)tau,al,be,asr
  enddo
  enddo
  enddo
 
!  write this set of fcs into a file to compare with original fcs
  open(932,file='trace_fc2_sr.dat')
  write(932,*)'# tau,taup,j,|r+taup-tau|,trace(fc2), rij,red_fold(r+taup-tau),fc2_sr(al,be)'
!  open(933,file='trace_phi.dat')
!  write(933,*)'# tau,taup,j,|r+taup-tau|,trace(phi),igrid,red_fold(r+taup-tau),phi(al,be)'
  do tau =1,natom_prim_cell
  do taup=1,natom_prim_cell
     dta = atompos(:,taup)-atompos(:,tau)
     do j=1,natoms
     if(iatomcell0(j).ne.taup) cycle
     rij = atompos(:,j)-atompos(:,tau)
     do igrid=1,nrgrid
        rr=rgrid(:,igrid)+dta 
        if (.not. (rr .myeq. rij)) cycle 
        call check_inside_ws(rr,rws26,insid,nboundary) 
        if(.not. insid) then
           rfold = fold_ws(rr,rws26,'r') 
        else
           rfold=rr
        endif
        fc2= trace(real(phip(tau,taup,:,:,igrid)))
        if(abs(fc2).gt.0.0001) then
           write(932,3)tau,taup,j,length(rfold),fc2,real(igrid),rij,cart2red(rfold,'r'),real(phip(tau,taup,:,:,igrid)) 
        endif
 !       fc2= trace(real(phi(tau,taup,:,:,igrid)))
 !       if(abs(fc2).gt.0.0001) then
 !          write(933,3)tau,taup,j,length(rfold),fc2,real(igrid),cart2red(rfold,'r'),real(phi(tau,taup,:,:,igrid)) 
 !       endif
     enddo
     enddo
  enddo
  enddo
  close(932)
  close(933)


3 format(3i5,99(f9.4))
4 format(4i5,99(f9.4))
5 format(a,3i5,9(f8.3))

 end subroutine correct_fcs

!-------------------------------------------------------

 subroutine make_grids
 use lattice
 use fourier
 use constants
 implicit none
 integer ngrd,cnt
 real(r15), allocatable :: grd(:,:),wei(:)
 integer ngrid
 type(vector) x01,x02,x03,x1,x2,x3

 write(ulog,*)' ENTERING MAKE_GRIDS'
 nrgrid=10000
 allocate(grd(3,nrgrid),wei(nrgrid))
 matr=cart_to_prim
 call get_26shortest_shell(r01,r02,r03,r0ws26,x01,x02,x03)
 call get_26shortest_shell(rs1,rs2,rs3,rws26,x1,x2,x3)
! call make_grid_weights_WS(r01,r02,r03,rs1,rs2,rs3,matr,nrgrid,grd,wei,'r',r0ws26,rws26) 
 call make_grid_weights_WS(x01,x02,x03,x1,x2,x3,matr,nrgrid,grd,wei,'r',r0ws26,rws26) 
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
 call make_grid_weights_WS(x1,x2,x3,x01,x02,x03,matr,nggrid,grd,wei,'g',gws26,g0ws26) 
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


