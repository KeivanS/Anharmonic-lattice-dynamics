 module atoms_force_constants
! the type atom_id0 concerns the properties of the primitive unit cell
! it s atoms, masses, neighbors, shells and force constants, whereas the
! type atomic refers to atoms in the supercell, their displacements and forces.
 use constants, only : r15
 use geometry
 implicit none
 integer natom_super_cell,nmax, fc2range
 integer :: fcflag(8) = 0  ! fcflag(3:8), structure.params line 12 after fc2range; 0 (default) = rule of fc34_range, 1 = nshells

!-------------------------
 type cell_id         ! id of atoms: position within cell, and cell coordinates
    integer n(3)
    integer tau
    integer atomposindx  ! corresponding index in the atompos list
 end type
!-------------------------
 type shell
    integer no_of_neighbors  ! within that shell
    real(r15) radius              ! radius of that shell
    type(cell_id) :: neighbors(296)  ! id of atoms in that shell ! keep up to 296 neighbors
 end type
!-------------------------
 type atomic
!   type(tensor1) u       ! displacement from equilibrium
!   type(tensor1) force   ! force applied on it
    type(vector) equilibrium_pos
    type(cell_id) cell
    integer at_type
    character name*2
    real(r15) mass,charge(3,3)
 end type
!-------------------------
 type atom_id0   ! all characteristics of an atom in the primitive central cell
    integer at_type,tau !  (n1,n2,n3)=0 for atoms in the primitive cell
    character name*2
    real(r15) mass,charge(3,3)
    integer nshells       ! how many shells are included for that atom = actual dim(shell)
    type(vector) equilibrium_pos
    type(shell), allocatable ::  shells(:)  ! of dim atom0(i)%nshells
!   type(shell), pointer ::  shells(:)  ! of dim atom0(i)%nshells
!    type(tensor2), pointer :: phi2(:)  ! refers to atom numbers neighbor of i
!    type(tensor3), pointer :: phi3(:,:)
!    type(tensor4), pointer :: phi4(:,:,:)
 end type
!-------------------------
! everything with 0 refers to the primitive cell

 integer  natom_type
 integer, allocatable:: natom(:),atom_type(:), map_rtau_sc(:,:)!label_of_primitive_atoms(:),
 real(r15), allocatable:: atompos0(:,:),mas(:),force(:,:,:),displ(:,:,:),vel(:,:),energy(:)
 real(r15), allocatable:: cur(:,:)
 character(2), allocatable:: atname(:)
 type (atom_id0), allocatable :: atom0(:)
 type (atomic), allocatable :: atom_sc(:),atom_shl(:)

!---------------------------
! Module for saved data
!      module force_constants_module
! maximum shells of nearest neighbors (along radial direction) , and actual # of neighborshells
      integer maxneighbors,maxshells
! maximum number of atoms out to maxneighbors
      integer maxatoms,imaxatm,imaxnei
! op_matrix(k,j,i), matrix for the ith point operator acting on cartesian coordinates
! iop_matrix(k,j,i), matrix for the ith point operator acting on reduced coordinates
      integer iop_matrix(3,3,48)
      double precision op_matrix(3,3,48)
! op_kmatrix(k,j,i), matrix for the ith point operator acting on cartesian kvector
! iop_kmatrix(k,j,i), matrix for the ith point operator acting on reduced kvector
      integer iop_kmatrix(3,3,48)
      double precision op_kmatrix(3,3,48)
      integer lattpgcount
! isgopcount, number of operators in space group
      integer isgopcount
! isgop(i), point operation in ith space group operator
      integer isgop(48)
! sgfract(j,i), jth cartesian coordinate of fractional in ith space group operator
      double precision sgfract(3,48)
! iatomop(j,i), point operator that takes jth atom into ith atom
      integer iatomop(:,:)
! atomopfract(k,j,i), kth cartesian coordinate of fractional to accompany iatomop
      double precision atomopfract(:,:,:)
! natom_prim_cell, number of atoms in the primitive unit cell
      integer natom_prim_cell
! natoms, number of atoms out to maxneighbors
      integer natoms
! atompos(j,i), jth cartesian coordinate of ith atom
      real(r15), allocatable :: atompos(:,:)
! iatomcell(j,i), linear combination of basis vectors of the primitive lattice
! that takes us to the unit cell containing the ith atom
      integer, allocatable ::  iatomcell(:,:),iatomcell0(:)
! iatomcell0(i), identity of atom in unit cell at origin equivalent to ith atom
! iatomneighbor(j,i), nearest neighbor shell of jth atom in primitive unit
! cell at the origin that contains the ith atom
      integer iatomneighbor(:,:)
      allocatable iatomneighbor,iatomop,atomopfract

      logical has_inversion

!      end module force_constants_module
! -------------------------------------

contains

! subroutine allocate_disp_forc(n)  ! needed for the md code
!    integer n
!   allocate( disp(3,n),forc(3,n),vel(3,n),cur(3,n),ncel1(n),ncel2(n),ncel3(n) )
!    allocate( disp(3,n),forc(3,n),vel(3,n),cur(3,n) )
! end subroutine allocate_disp_forc

 subroutine allocate_edf(nat,ncfg)
    integer nat,ncfg
    allocate( displ(3,nat,ncfg),force(3,nat,ncfg),energy(ncfg) )
 end subroutine allocate_edf

 subroutine allocate_mass(n)
    integer n
    allocate( mas(n) ,natom(n), atname(n) )
 end subroutine allocate_mass

 subroutine allocate_primcell(n0)
 integer n0
!print*,'n0=',n0
 allocate( atom0(n0),atompos0(3,n0), atom_type(n0) )  !label_of_primitive_atoms(n0),
 end subroutine allocate_primcell

 subroutine allocate_supercell(n)
 integer n
 allocate( atom_sc(n) )
 end subroutine allocate_supercell

!------------------------------------
!subroutine set_neighbor_list (maxshell)
!  Generates a set of neighbors to the primitive cell (up to maxshell, which should be large enough)
!  and finds which shell corresponds to the cutoff rcut set by the WS cell
!  this is needed in order to know what force constants to associate
!  to a pair ij, or ijk or ijkl
!  if a neighbor j is known then one can get the vect(i,j) between
!  their equilibrium positions
! use atoms_force_constants , only : natom_prim_cell
!use params
!use lattice
!use ios
!use geometry
!implicit none
!! integer, parameter :: max=500
!integer, intent(in) :: maxshell
!real(r15), intent(in) :: rcut
!integer i0,j,shel_count,counter,msort(maxatoms),l
!real(r15) dist(maxatoms),d_old,d_min,rmax

!  allocate( atom0(1:natom_prim_cell)%shells(maxshells) )
! rmax = 0 ; dist = 1d10
! i0loop: do i0=1,natom_prim_cell
!   allocate( atom0(i0)%shells(0:maxshell) )
!   do j=1,natoms
!      call calculate_distance(i0,j,atompos,maxatoms,dist(j) )
!      if ( iatomneighbor(i0,j) .eq. nshells(2) .and. dist(j) .gt. rmax ) then
!      if ( iatomneighbor(i0,j) .eq. nshells(2,i0) .and. dist(j) .gt. rmax ) then
!           rmax = dist(j)
!      endif
!   enddo
!   call sort(natoms,dist,msort,maxatoms)

!   shel_count = -1
!   d_old = 0
!   write(ulog,*)' ========================================================='
!   write(ulog,*)' neighborlist of atom i0 =',i0

!   jloop: do j=1,min(500,natoms)
!      d_min = dist(msort(j))  ! start with the shortest distance
!      if ( (d_min .myeq. d_old) .and. (shel_count.ne.-1) ) then   ! same shell
!         counter = counter + 1
!      else    ! new shell
!         shel_count = shel_count+1
!         counter = 1
!         if ( shel_count .gt. maxshell ) then
!            write(ulog,*)maxshell,' shells completed'
!            write(ulog,*)'shel_count=',shel_count,' exceeded it for j=',j
!            exit jloop
!         endif
!         d_old = d_min
!      endif
!      if ( counter .gt. 296) then
!         write(ulog,*) ' counter in neighbors exceeded 296 ', counter
!         stop
!      endif
!      atom0(i0)%shells(shel_count)%no_of_neighbors = counter
!      atom0(i0)%shells(shel_count)%radius = d_min
!      atom0(i0)%shells(shel_count)%neighbors(counter)%tau = iatomcell0(msort(j))
!      atom0(i0)%shells(shel_count)%neighbors(counter)%n   = iatomcell(:,msort(j))
!      if(verbose) then
!        write(ulog,*)' ------------------------------------'
!         write(ulog,5)' count, j, of tau,n1,n2,n3,shel#,nb#,dij,rj =',j,msort(j),iatomcell0(msort(j)),   &
!                iatomcell(:,msort(j)),iatomneighbor(i0,msort(j)),counter,d_min,(atompos(l,msort(j)),l=1,3)
!      endif
!   enddo jloop

!   do shel_count = 0 , maxshell
!      write(ulog,*)'shell#, neigh#=',shel_count,atom0(i0)%shells(shel_count)%no_of_neighbors
!   enddo
!  also initialize atom0%equilibrium_pos

!  shift the position of the first atom to the  origin
! do i=1,natom_prim_cell
!    atompos(:,i) = atompos(:,i)-atompos(:,1)
!     atom0(i0)%equilibrium_pos = atompos0(1,i0)*r01 + atompos0(2,i0)*r02  &
! &                             + atompos0(3,i0)*r03
!   atom0(i0)%equilibrium_pos%x = atompos(1,i0)
!   atom0(i0)%equilibrium_pos%y = atompos(2,i0)
!   atom0(i0)%equilibrium_pos%z = atompos(3,i0)
! enddo
!  they might not be inside the primitive cell

!enddo i0loop
! rcut(2) = rmax

! format(a,8(2x,i4))
! format(a,' (',i4,2(1x,i4),1x,')')
! format(a,1x,i4,' (',f8.4,2(1x,f8.4),1x,')')
! format(a,3(2x,i4),' (',i4,2(1x,i4),1x,')',2(1x,i3),1x,f10.7,2x,9(1x,f9.6))
! format(a,1x,f10.7,2x,9(1x,f9.6))

!end subroutine set_neighbor_list

!--------------------------------
     subroutine test_ntau
     use ios, only : ulog
     implicit none
     integer i,tau1,n1(3),tau2,n2(3),iatom1,iatom2,err

     err=0
     do i=1,natoms
    !     call findatom_cart(v,iatom(irank))
        call get_n_tau_r    (atompos(:,i),n1,tau1)
        call get_n_tau_r_mod(atompos(:,i),n2,tau2)
        call find_atompos(n1,tau1,iatom1)
        call find_atompos(n2,tau2,iatom2)
        if(.not.(n1.myeq.n2) .or. tau1.ne.tau2) then
           write(*   ,4)'TEST_NTAU: i,i1,i2,tau1,tau2, n1,n2=', &
&                   i,iatom1,iatom2,tau1,tau2,n1,n2
           write(ulog,4)'TEST_NTAU: i,i1,i2,tau1,tau2, n1,n2=', &
&                   i,iatom1,iatom2,tau1,tau2,n1,n2
           err=1
        endif
     enddo
     if(err.eq.0) then
        write(ulog,*) 'Both get_n_tau_r and get_n_tau_r_mod yield the same results '
     else
        write(ulog,*) ' CAREFUL: get_n_tau_r and get_n_tau_r_mod yield different results! '
     endif
 4   format(a,3(1x,i5),3x,2i3,2(4x,3(1x,i2)))

     end subroutine test_ntau

 end module atoms_force_constants
!===========================================================
