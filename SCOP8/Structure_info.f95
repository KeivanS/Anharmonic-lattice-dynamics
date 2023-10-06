MODULE Structure_info
!!This module is used for reading info from lat_fc.dat and store those info
    USE lattice
    USE atoms_force_constants
    USE svd_stuff !for nterms(4), ngroups(4), etc

    IMPLICIT NONE
    PUBLIC
    INTEGER, ALLOCATABLE :: at_label(:)
    INTEGER :: atom_number,cell_number,tot_atom_number
    INTEGER :: variational_parameters_size(3)
    INTEGER,ALLOCATABLE,DIMENSION(:) :: atoms_fc2,atoms_fc3,atoms_fc4
    LOGICAL :: at_name_got = .False.
    !INTEGER,PARAMETER :: d=3
    REAL(8),ALLOCATABLE,DIMENSION(:,:) :: cell_vec,trans_vec,cell_save,trans_save
!---------------------------------------------------------------------------------
   TYPE PAtom !read in atom list in the center primitive cell, mainly for store mass and atom_type and tau
      REAL(8) :: pos_tau(d)
      INTEGER :: atom_type !atom_type is really the atom label for primitive cell i.e. Si#1 and Si#2
      INTEGER :: ttyp  !ttyp is really the element type of atom i.e. Si
      REAL(8) :: mass
      character(LEN=2) :: name
   END TYPE

   TYPE EAtom  !read every atoms, get and store the corresponding R vector and tau vector for later use
      INTEGER :: label_number
      REAL(8) :: x,y,z !this is only for d=3
      INTEGER :: type_tau,type_R,n1,n2,n3 !type_tau is the atom_type, type_R is the cell_number

      REAL(8),DIMENSION(d) :: R   !this is the 'cell_vec' for every atom, it's also needed during the calculation
      REAL(8),DIMENSION(d) :: tau !this is the 'pos_tau' for every atom,
   END TYPE
!---------------------------------------------------------------------------------
    TYPE(PAtom),DIMENSION(:),ALLOCATABLE :: iatom
    TYPE(EAtom),DIMENSION(:),ALLOCATABLE :: every_atom, prev_atom

    !d=3

CONTAINS
!************************************************************************************
!=============================original subs==========================================
subroutine read_params_org
!! legacy subroutine, not used
 use io2
 use om_dos
 use params
 use lattice
 use phi3
 use kpoints
 use atoms_force_constants
 implicit none
! integer i
 integer ulog2, scalelengths
 ulog2 = 444
OPEN(ulog2, FILE='read_params_log.dat',status='unknown')

  open(uparams,file='params.phon',status='old')
  write(6,*) 'READ_PARAMS: opening params.phon'
!  read(uparams,*) natom_type   ! # of different elements present in prim cell
!  write(6,*) 'READ_PARAMS: just read natom_type ',natom_type
! call allocate_mass(natom_type)
! read(uparams,*) (mas(i),i=1,natom_type)  ! in the same order as in POTCAR
! read(uparams,*) (atname(i),i=1,natom_type)  ! in the same order as in POTCAR
! read(uparams,*) natoms        ! # of atoms in primitive cell

  read(uparams,*)nc1,nc2,nc3,dos_bs_only ! coarse kmesh for big calculations
  NC(1)=nc1 ;NC(2)=nc2 ;NC(3)=nc3
  write(ulog2,*) 'READ_PARAMS: just read nib ',nc1
  read(uparams,*)lamin,lamax
  read(uparams,*)shftx,shfty,shftz ! kmesh for dos calculation
  read(uparams,*)wmesh,wmax  ! no of w mesh points and max frequency for DOS
  write(ulog2,*) 'READ_PARAMS: just read wmax ',wmax
  read(uparams,*)width,etaz  ! width of gaussian broadening for DOS,imag part of om
  write(ulog2,*) 'READ_PARAMS: just read width,etaz ',width,etaz
  read(uparams,*)tau0        ! large relaxation time: 1/tau0 to be added to gamma
  write(ulog2,*) 'READ_PARAMS: just read tau0 ',tau0
  read(uparams,*)verbose
  read(uparams,*)wshift      ! small shift added to make negative modes go away
  write(ulog2,*) 'READ_PARAMS: just read wshift ',wshift
  read(uparams,*)tmin,tmax,ntemp      ! temps in Kelvin to calculate thermal expansion and other thermal ppties
  read(uparams,*)iter,readv3,usetetra,isvd      ! if=1 then read from a file, else generate
  read(uparams,*)ncpu
  read(uparams,'(a)')v3path
  write(ulog2,*) 'READ_PARAMS: path is:',v3path
  read(uparams,*)max_iter, n_dig_acc
!, conv_error, conv_max_error, conv_diff, conv_max_diff, &
!  & conv_diff_kap, conv_max_diff_kap,conv_iter,update
   ! maximum iteration number for iterative solution, convergence criteria, update ratio
  read(uparams,*)v3_threshold        ! in 1/cm keep all v33 of norm above v3_threshold
  read(uparams,*)classical     ! if =1 then use classical distribution (kT/hw)
  read(uparams,*)cal_cross,qcros     ! if =1 then calculate the cross section at qcros (reduced U)
  read(uparams,*)threemtrx          ! if =1 then calculate the 3-phonon matrix elements
  read(uparams,*)scalelengths       ! multiply all lattice and atomic coordinates by scalelength

  close(uparams)

  write(ulog2,*) 'READ_PARAMS: read and kpoints allocated'
  write(ulog2,3)'nc1,nc2,nc3=',nc1,nc2,nc3
  write(ulog2,*)'wmesh, wmax=',wmesh,wmax
  write(ulog2,*)'readv3=',readv3  !,writev3
! etaz = max(etaz,wmax/wmesh)
! write(ulog2,*)' in case etaz is too small we use the following'
! write(ulog2,*)'width,etaz =',width,etaz
  write(ulog2,*)'Tmin,Tmax(K=',tmin,tmax
  if (classical .eq. 1) then
     write(ulog2,*)'Calculation is done in the CLASSICAL limit'
  else
     write(ulog2,*)'Calculation is done in the QUANTUM limit'
  endif
  if (dos_bs_only) then
     write(ulog2,*)' DOS will only be calculated (on the coarse mesh) ',nc1,nc2,nc3
!    write(ulog2,*)' and the program will stop after that!'
  endif
CLOSE(ulog2)
print*,' file params.phon read'

3 format(a,6(1x,i6))
 end subroutine read_params_org
 !--------------------------------------------------------------------------------
 subroutine read_input_fit_org
  !! legacy subroutine not used
 use io2
 use params
 use lattice
 use force_constants_module
 use atoms_force_constants
 implicit none
 integer i,counter,label
 real(8) junk
 character jjj*1

 integer ulog2
 ulog2 = 444
OPEN(ulog2, FILE='read_input_fit_log.dat',status='unknown')

 open(uparams,file='params.inp',status='old')

 read(uparams,*) latticeparameters   ! (a,b,c,alpha,beta,gamma)
 read(uparams,*) primitivelattice     ! prim latt vectors in terms of conventional above
 read(uparams,*) lattice_parameter
 read(uparams,*) maxneighbors
 read(uparams,*) include_fc    ! if=1 include this rank
! read(uparams,*) nshells       ! # of neigr-shells to use for each rank
 read(uparams,*)  i !junk
! read(uparams,*) junk
! read(uparams,*) junk
 read(uparams,*) tolerance  ! this really the flags...
 read(uparams,*) svdcut
 read(uparams,*) junk
 read(uparams,*) natom_type   ! # of different elements present in prim cell
!read(uparams,*) jjj
! allocate(natom(natom_type))
 write(ulog2,*) 'READ_INPUT_FIT: natom_type ',natom_type
 call allocate_mass(natom_type)
 read(uparams,*) (mas(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) (atname(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) natoms0        ! # of atoms in primitive cell
 read(uparams,*) nshells(1,1:natoms0)       ! # of neigr-shells to use for each rank
 read(uparams,*) nshells(2,1:natoms0)       ! # of neigr-shells to use for each rank
 read(uparams,*) nshells(3,1:natoms0)       ! # of neigr-shells to use for each rank
 read(uparams,*) nshells(4,1:natoms0)       ! # of neigr-shells to use for each rank
 write(ulog2,*)'reading ',natoms0,' atoms in the primitive cell'
 call allocate_primcell(natoms0)
! indices of atoms in primitive cell must be same order as in POSCAR1
 counter = 0
 natom = 0
 do i=1,natoms0
! positions here must be D format in conventional cell for use by fcs_init
    read(uparams,*) label,atom_type(label),atom0(label)%equilibrium_pos
! this is useless because natom(:) is read from POSCAR1
    if (atom_type(i).eq.counter) then
        natom(atom_type(i)) = natom(atom_type(i))+1
    else
        natom(atom_type(i)) = 1
    Endif
    counter = atom_type(i)
 write(ulog2,*)'reading  atom #',i, counter
 enddo

 write(ulog2,*)'reading done, closing the params.inp file'
 close(uparams)

 do i=1,natoms0
    atom0(i)%name = atname(atom_type(i))
    atom0(i)%at_type  = atom_type(i)
    atom0(i)%mass = mas(atom_type(i))
    atompos0(:,i) = atom0(i)%equilibrium_pos%component(:)
 enddo

 maxshells = maxneighbors

 do i=1,4
!   if (nshells(i) .gt. maxshells) then
    if (maxval(nshells(i,:)) .gt. maxshells) then
!      write(ulog2,*)' READ_INPUT: nshell> maxshells ',i,nshells(i),maxshells
       write(ulog2,*)' READ_INPUT: nshell> maxshells ',i,maxval(nshells(i,:)),maxshells
       write(ulog2,*)' Either increase maxshells or lower nshell for that rank'
       stop
    endif
 enddo
CLOSE(ulog2)
print*,' file params.inp read'

end subroutine read_input_fit_org
!------------------------------------------------------------------------------------
subroutine read_lattice_org
!! legacy subroutine not used
 use io2
 use params
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use constants
 implicit none
 character line*90,name*2
 integer i,j,tt,ttyp,natoms2,n1  !,at_count,n2,n3
INTEGER :: scalelengths
 real(8) mss !a1,a2,a3,om,a,b,c,
 type(vector) tau1,vv

 integer ulog2
 ulog2 = 444
OPEN(ulog2, FILE='read_lattice_log.dat',status='unknown')

open(ufc0,file='lat_fc.dat' ,status='old')

do while (line(1:22) .ne. '  Crystal data: transl')
 read(ufc0,'(a)') line
enddo
 read(ufc0,*) r1%component   ! coordinates of the primitive cell vectors
 read(ufc0,*) r2%component
 read(ufc0,*) r3%component


! now generate the lattice ppties in real and recip spaces
  call make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
  call make_reciprocal_lattice(g1,g2,g3,rr1,rr2,rr3)
  call calculate_volume(r1,r2,r3,volume_r)
  call calculate_volume(g1,g2,g3,volume_g)
  box(1) = length(r1)
  box(2) = length(r2)
  box(3) = length(r3)
  boxg(1) = length(g1)
  boxg(2) = length(g2)
  boxg(3) = length(g3)
  write(ulog2,3)' r1= ',r1
  write(ulog2,3)' r2= ',r2
  write(ulog2,3)' r3= ',r3
  write(ulog2,3)' box  = ',box
  write(ulog2,3)' g1= ',g1
  write(ulog2,3)' g2= ',g2
  write(ulog2,3)' g3= ',g3
  write(ulog2,3)' boxg = ',boxg
  write(ulog2,3)' volume_r,volume_g = ',volume_r,volume_g

 read(ufc0,*) line
! read(ufc,*) natoms0
! ############ check compatibiliy with input file "params.inp" ###########
 read(ufc0,*) n1
 if (n1 .ne. natoms0) then
    write(ulog2,*)' error in reading natoms0 ',n1,natoms0
    stop
 endif
! call allocate_primcell(natoms0)
 do i=1,natoms0
!  read(ufc1,*)atom0(i)%tau, atom_type(i),atom0(i)%equilibrium_pos,atom0(i)%mass
    read(ufc0,*)tt,name,ttyp,tau1,mss  ! this has the coordinates and atomic types and masses

!   atompos0(:,i) = atom0(i)%equilibrium_pos
!   vv%x= tau1.dot.g1 ; vv%y= tau1.dot.g2 ; vv%z= tau1.dot.g3
    vv%component(1) = tau1.dot.g1; vv%component(2) = tau1.dot.g2; vv%component(3) = tau1.dot.g3
    vv = (1/2d0/pi)*vv
    if(tt .eq. i .and. ttyp .eq. atom0(i)%at_type .and. atom0(i)%mass .eq. mss) then
!&     length(atom0(i)%equilibrium_pos-vv).lt.1d-4 .and.
      write(ulog2,*)'####### compatibility with params.inp checked  ######### ',i
    else
      write(ulog2,*)'READ_LATTICE: data in params.inp and lat_fc.dat are not compatible'
      write(ulog2,*)' i, atom type, mass '
      write(ulog2,*) i, tt, atom0(i)%at_type,ttyp, atom0(i)%mass,mss
! fix it so that if eq pos are the same module translarion vectors, it's OK
      write(ulog2,3)' reduced coordinates ',atom0(i)%equilibrium_pos,vv
      stop
    endif
 enddo

 read(ufc0,*) line
 read(ufc0,*) line
 read(ufc0,*) (include_fc(j),j=1,4)
 write(ulog2,*) (include_fc(j),j=1,4)
 read(ufc0,*) line
 write(ulog2,*) line
 read(ufc0,*) (nterms(j),j=1,4)
 write(ulog2,*) (nterms(j),j=1,4)
 read(ufc0,*) line
 write(ulog2,*) line
 read(ufc0,*) (ngroups(j),j=1,4)
 write(ulog2,*) (ngroups(j),j=1,4)
 read(ufc0,*) line
 write(ulog2,*) line
 read(ufc0,*) natoms2
 write(ulog2,*) natoms2

  maxatoms=2500; imaxat=1
  do while (imaxat.ne.0)
     maxatoms=maxatoms+300
     write(6,*)' maxatoms=',maxatoms
     write(ulog2,*)' maxatoms=',maxatoms
     if (allocated(iatomcell0)) deallocate(iatomcell0)
     if (allocated(iatomcell))  deallocate(iatomcell)
     if (allocated(atompos))    deallocate(atompos)
     allocate(atompos(3,maxatoms), iatomcell(3,maxatoms),iatomcell0(maxatoms))
! find symmetry operations of the crystal
     call force_constants_init(latticeparameters,primitivelattice,natoms0,  &
      &     atom_type,atompos0)  !(n,natom_type,atompos0)
  enddo

 if (natoms.ne.natoms2) then
    write(ulog2,*) "natoms read in lat_fc=",natoms2
    write(ulog2,*) "while natoms generated by fc_init=",natoms
    write(ulog2,*) "check the number of shells in params.inp, "
    write(ulog2,*) "it is probably too large if natoms in fc_init is larger than in lat_fc"
    write(ulog2,*) "in which case the program will stop"
 endif
 STOP
 deallocate(atompos, iatomcell,iatomcell0)
 allocate(atompos(3,natoms), iatomcell(3,natoms),iatomcell0(natoms))
 do j=1,natoms
    read(ufc0,*,end=99)i,atompos(:,i), iatomcell0(i),iatomcell(:,i)

    atompos(:,i)=atompos(:,i)  * scalelengths  !! added by K1 5/6/2015

    write(ulog2,4)i,atompos(:,i), iatomcell0(i),iatomcell(:,i)
 enddo
3 format(a,9(2x,f11.6))
4 format(i6,3(2x,g12.5),4(3x,i4))

CLOSE(ulog2)
 close(ufc0)
 write(ulog2,*)' COORDINATES read successfully '
 return

99 write(*,*) "End of lat_fc.dat file reached!"
   write(*,*) "check the number of shells in params.inp, it is probably too large!"
   stop

 end subroutine read_lattice_org
!===================================================================================
subroutine write_neighbors
!! legacy subroutine not used
! redundant module usage
! use atoms_force_constants
! use force_constants_module
! use params
! use ios
 implicit none
 integer i0,shel_count,j,nm(3),ta,nbmx
 real(8) dij

 write(ulog,*)' rcut(2) corresponding to nshell(2)=',rcut(2)
 write(ulog,*)' ************ Writing the neighborlist ************** '
 do i0=1,natoms0
    write(ulog,*)' Neighbors of atom number ',i0
    do shel_count=1,maxshells
       nbmx = atom0(i0)%shells(shel_count)%no_of_neighbors
       dij  = atom0(i0)%shells(shel_count)%rij
       do j=1,min(nbmx,500)  ! write the first 500 neighbors
          ta =  atom0(i0)%shells(shel_count)%neighbors(j)%tau
          nm =  atom0(i0)%shells(shel_count)%neighbors(j)%n
          write(ulog,3)' shell,radius,neighbor,atomid=',shel_count,dij,j,ta,nm
       enddo
       write(333,2)shel_count,dij,j-1
    enddo
 enddo

 write(ulog,*)' ************ End of the neighborlist ************** '
2 format(i4,2x,f8.4,2x,i4)
3 format(a,2x,i3,2x,f8.4,2x,i3,4x,i3,' (',3(1x,i4),')')
 end subroutine write_neighbors
!========================================================================
subroutine write_input_fit
!! legacy subroutine not used
    use io2
    implicit none
    integer i,counter,label,unit_params
    real(8) scal,junk
    character jjj*1

    unit_params=49
    open(unit_params,file='params.inp')

    write(unit_params,*) latticeparameters   ! (a,b,c,alpha,beta,gamma)
    write(unit_params,*) primitivelattice     ! prim latt vectors(3:3) in terms of conventional above
    write(unit_params,*) lattice_parameter !scale factor for lattparams (must be consistent with POSCAR data)
    scal = lattice_parameter
    write(unit_params,*) maxneighbors  !used in lat_fc.dat
    write(unit_params,*) include_fc    ! if=1 include this rank
! read(uparams,*) nshells       ! # of neigr-shells to use for each rank
    write(unit_params,*)  itrans,irot,ihuang,enforce_inv !not used? flags for including translational and rotational invce constraints.
! read(uparams,*) junk
! read(uparams,*) junk
    write(unit_params,*) tolerance,margin  ! this really the flags...&
                               !&tolerance for equating (two coordinates in general), margin for eliminating a FC
    write(unit_params,*) svdcut   !svd cutoff for the smallest eigenvalue to be included
    write(unit_params,*) fdfiles     !number of force-displacement files
    write(unit_params,*) natom_type   ! # of different elements present in prim cell, I call it 'atom_number'
!read(uparams,*) jjj
! allocate(natom(natom_type))

    write(unit_params,*) (mas(i),i=1,natom_type)  ! in the same order as in POTCAR, declared in module [atoms_force_constants]

    !how am I supposed to write this?? lmao fortran
    !write(unit_params,*) (atname(i),i=1,natom_type)  ! in the same order as in POTCAR, declared in module [atoms_force_constants]

    write(unit_params,*) natoms0        ! # of atoms in primitive cell,coordinates are in reduced units (of cubic cell)
 !nshells(4,20) are declared in module[params]
    write(unit_params,*) nshells(1,1:natoms0)       ! # of neigr-shells to use for each rank
    write(unit_params,*) nshells(2,1:natoms0)       ! # of neigr-shells to use for each rank
    write(unit_params,*) nshells(3,1:natoms0)       ! # of neigr-shells to use for each rank
    write(unit_params,*) nshells(4,1:natoms0)       ! # of neigr-shells to use for each rank

    do i=1,natoms0
    ! positions here must be D format in conventional cell for use by fcs_init
        label = at_label(i)
        write(unit_params,*) label,atom_type(label),atom0(label)%equilibrium_pos%component
        ! this is useless because natom(:) is read from POSCAR1
        if (atom_type(i).eq.counter) then
            natom(atom_type(i)) = natom(atom_type(i))+1
        else
            natom(atom_type(i)) = 1
        Endif
        counter = atom_type(i)
        WRITE(*,*) 'writing  atom #',i, counter
    enddo

    close(unit_params)
end subroutine write_input_fit
!-----------------------------------------------------------------------
 subroutine read_input_fit
 !! modified legacy subroutine, read params.inp
!! major
    use io2 !the uparams is causing trouble because it's declared both in [ios] & [io2] with different value
    !use params !redundant
    !use lattice
    !use force_constants_module !redundant
    !use atoms_force_constants
    implicit none
    integer i,counter,label,unit_params
    integer nouse1, nouse2, nouse3
    real(8) scal,junk
    character jjj*1

    unit_params=49
    open(unit_params,file='structure.params',status='old')

    !---hard-coded params that used to be set in 'params.inp'
    maxneighbors = 90 !I just default set it to 10 for now !Modify0726
    tolerance = 0.001; margin = 0.00001
    maxshells = 27 ! this is added from the fc234 main program to here !0726

    read(unit_params,*) latticeparameters   ! (a,b,c,alpha,beta,gamma)
    read(unit_params,*) primitivelattice     ! prim latt vectors(3:3) in terms of conventional above
    read(unit_params,*) lattice_parameter !scale factor for lattparams (must be consistent with POSCAR data)
    scal = lattice_parameter
!    read(unit_params,*) maxneighbors  !omitted in the new file but what value should it be?
    read(unit_params,*) include_fc    ! if=1 include this rank
    read(unit_params,*)  itrans,irot,ihuang,enforce_inv !not used? flags for including translational and rotational invce constraints.

    read(unit_params,*) nouse1, nouse2!skip the temperature line

!    read(unit_params,*) tolerance,margin  ! this really the flags...&
                               !&tolerance for equating (two coordinates in general), margin for eliminating a FC
!    read(unit_params,*) svdcut   !svd cutoff for the smallest eigenvalue to be included
    read(unit_params,*) fdfiles, verbose     !number of force-displacement files
    read(unit_params,*) natom_type   ! # of different elements present in prim cell, I call it 'atom_number'
    write(ulog,*) 'READ_INPUT_FIT: natom_type ',natom_type

    call allocate_mass(natom_type) !subroutine in module [atom_force_constants]
    read(unit_params,*) (mas(i),i=1,natom_type)  ! in the same order as in POTCAR, declared in module [atoms_force_constants]

    IF(.not.at_name_got) THEN
        read(unit_params,*) (atname(i),i=1,natom_type)  ! in the same order as in POTCAR, declared in module [atoms_force_constants]
        at_name_got = .TRUE.
    ENDIF

    read(unit_params,*) natoms0        ! # of atoms in primitive cell,coordinates are in reduced units (of cubic cell)

    read(unit_params,*) nouse3!skip the flag line

!    read(unit_params,*) nshells(1,1:natoms0)       ! # of neigr-shells to use for rank 1 is omitted in the new file, default set to 1
    do i=1, natoms0
        nshells(1,i) = 1
    end do

    read(unit_params,*) nshells(2,1:natoms0)       ! # of neigr-shells to use for each rank
    read(unit_params,*) nshells(3,1:natoms0)       ! # of neigr-shells to use for each rank
    read(unit_params,*) nshells(4,1:natoms0)       ! # of neigr-shells to use for each rank
    write(*,*)'reading ',natoms0,' atoms in the primitive cell'
!-------------missed-------------------
 write(*,*) svdcut,'    cutoff for smallest eigenvalue w to be included'
 write(*,*) tolerance,'   tolerance for equating two coordinates '
 write(*,*) include_fc,'  which ranks of FCs to include '
 write(*,*) nshells(1,:),'  how many shells to include for each rank '
 write(*,*) nshells(2,:),'  how many shells to include for each rank '
 write(*,*) nshells(3,:),'  how many shells to include for each rank '
 write(*,*) nshells(4,:),'  how many shells to include for each rank '
 write(*,*) itrans,irot,ihuang,enforce_inv,'  transl, rot and Huang invce, enforcing inv'

 latticeparameters(1:3) = latticeparameters(1:3)*scal
!-------------missed end---------------
!---------small modification----------------------
if(.not.allocated(at_label)) allocate(at_label(natoms0))
!-------------------------------------------------
    call allocate_primcell(natoms0)
    ! indices of atoms in primitive cell must be same order as in POSCAR1
    counter = 0
    natom = 0
    do i=1,natoms0
    ! positions here must be D format in conventional cell for use by fcs_init
        read(unit_params,*) label,atom_type(label),atom0(label)%equilibrium_pos%component
        at_label(i) = label
        ! this is useless because natom(:) is read from POSCAR1
        if (atom_type(i).eq.counter) then
            natom(atom_type(i)) = natom(atom_type(i))+1
        else
            natom(atom_type(i)) = 1
        Endif
        counter = atom_type(i)
        write(ulog,*)'reading  atom #',i, counter
        WRITE(*,*) 'reading  atom #',i, counter
    enddo

    write(ulog,*)'reading done, closing the params.inp file'
    WRITE(*,*) 'reading done, closing the params.inp file'
!    close(uparams)

    do i=1,natoms0
        atom0(i)%name = atname(atom_type(i))
        atom0(i)%at_type  = atom_type(i)
        atom0(i)%mass = mas(atom_type(i))
        atompos0(:,i) = atom0(i)%equilibrium_pos%component(:)
    enddo

    do i=1,4
        !if (nshells(i) .gt. maxshells) then
        if (maxval(nshells(i,:)) .gt. maxshells) then
        !write(ulog,*)' READ_INPUT: nshell> maxshells ',i,nshells(i),maxshells
            write(ulog,*)' READ_INPUT: nshell> maxshells ',i,maxval(nshells(i,:)),maxshells
            write(ulog,*)' Either increase maxshells or lower nshell for that rank'
        stop
        endif
    enddo
CLOSE(unit_params)
end subroutine read_input_fit
!==========================================================
    SUBROUTINE read_structure
    !! my subroutine to read lat_fc.dat and store info into corresponding TYPE
!! major
        IMPLICIT NONE

        CHARACTER(LEN=90) :: line
        !CHARACTER(LEN=2) :: name
        INTEGER :: unit_number,unit_log,i,j,direction, maxshell!MARK: 09/04/2023
        REAL(8) :: check(d)
        TYPE(vector) :: vv
        LOGICAL :: condition

        check=0
        unit_log=30
        OPEN(unit_log,file='logfile.txt',status='unknown')

        unit_number=50
        OPEN(unit_number, file='lat_fc.dat',status='old',action='read')

        DO WHILE(line(1:22).ne.'  Crystal data: transl') !notice the 2 space bar
            READ(unit_number,'(a)') line
        END DO

        ALLOCATE(trans_vec(d,d))
        ALLOCATE(trans_save(d,d))
        DO i=1,d
            READ(unit_number,*) trans_vec(:,i)
        END DO
        trans_save=trans_vec
!**********************************************************************
r1%component=trans_vec(:,1)
r2%component=trans_vec(:,2)
r3%component=trans_vec(:,3)
! now generate the lattice ppties in real and recip spaces
  call make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
  call make_reciprocal_lattice(g1,g2,g3,rr1,rr2,rr3)
  call calculate_volume(r1,r2,r3,volume_r)
  call calculate_volume(g1,g2,g3,volume_g)
  write(unit_log,3)' r1= ',r1
  write(unit_log,3)' r2= ',r2
  write(unit_log,3)' r3= ',r3
  write(unit_log,3)' box  = ',box
  write(unit_log,3)' g1= ',g1
  write(unit_log,3)' g2= ',g2
  write(unit_log,3)' g3= ',g3
  write(unit_log,3)' rr1= ',rr1
  write(unit_log,3)' rr2= ',rr2
  write(unit_log,3)' rr3= ',rr3
  write(unit_log,3)' boxg = ',boxg
  write(unit_log,3)' volume_r,volume_g = ',volume_r,volume_g
!**********************************************************************
        READ(unit_number,*) line
        ! read(ufc,*) natoms0
! ############ check compatibility with input file "params.inp" ###########
        READ(unit_number,*) atom_number
        IF(atom_number.ne.natoms0) THEN
            WRITE(unit_log,*) ' error in reading natoms0 ',atom_number,natoms0
            WRITE(*,*) ' error in reading natoms0 ',atom_number,natoms0
            STOP
        END IF
! allocate for primitive cell
        ALLOCATE(iatom(atom_number))
        DO i=1,atom_number
            READ(unit_number,*) iatom(i)%atom_type,iatom(i)%name,iatom(i)%ttyp,iatom(i)%pos_tau,iatom(i)%mass
            DO direction=1,d
                vv%component(direction)=iatom(i)%pos_tau.dot.g1
            END DO !direction/components loop
            vv%component=(1/2d0/pi)*vv%component
!#################check further compatibility with input file "params.inp" ##############
            condition=iatom(i)%atom_type.eq.i&
            & .and. iatom(i)%ttyp.eq.atom0(i)%at_type .and. atom0(i)%mass.eq.iatom(i)%mass
            IF(condition) THEN
                WRITE(unit_log,*) '####### compatibility with params.inp checked  ######### ',i
            ELSE
                WRITE(unit_log,*)'READ_LATTICE: data in params.inp and lat_fc.dat are not compatible'
                WRITE(unit_log,*)' i, atom type, mass '
                WRITE(unit_log,*) i,iatom(i)%atom_type,atom0(i)%at_type,iatom(i)%ttyp,atom0(i)%mass,iatom(i)%mass
            ! fix it so that if eq pos are the same module translarion vectors, it's OK
                WRITE(ulog,3)' reduced coordinates ',atom0(i)%equilibrium_pos%component,vv%component
                STOP
            END IF !end compatibility check
        END DO !iatom loop



        READ(unit_number,*) line
        READ(unit_number,*) line
        READ(unit_number,*) (include_fc(j),j=1,4)
        WRITE(unit_log,*) (include_fc(j),j=1,4)
        READ(unit_number,*) line
        WRITE(unit_log,*) line
        READ(unit_number,*) (nterms(j),j=1,4)
        WRITE(unit_log,*) (nterms(j),j=1,4)
        READ(unit_number,*) line
        WRITE(unit_log,*) line
        READ(unit_number,*) (ngroups(j),j=1,4)
        WRITE(unit_log,*) (ngroups(j),j=1,4)
        READ(unit_number,*) line
        WRITE(unit_log,*) line

        READ(unit_number,*) maxshell, tot_atom_number !MARK: 09/04/2023
        WRITE(unit_log,*)'Total Atom Number is', tot_atom_number
        ALLOCATE(every_atom(tot_atom_number))
!**********************************************************
        !this part has nothing to do with 'reading' and is directly copied from fc234 main program
        !notice that original comment reads:
        !"This does not work: fcinit must be called once or latparam(4:6) will be overwritten"
        maxatoms=2500; imaxat=1
        DO WHILE (imaxat.ne.0) !imaxat is a flag number, it equals 0 when some conditions are satisfied
            maxatoms=maxatoms+300 !???
            WRITE(6,*)' maxatoms=',maxatoms
            WRITE(unit_log,*)' maxatoms=',maxatoms
            IF (allocated(iatomcell0)) DEALLOCATE(iatomcell0)
            IF (allocated(iatomcell))  DEALLOCATE(iatomcell)
            IF (allocated(atompos))    DEALLOCATE(atompos)
            ALLOCATE(atompos(3,maxatoms), iatomcell(3,maxatoms),iatomcell0(maxatoms))
            CALL force_constants_init(latticeparameters,primitivelattice,natoms0,  &
            &     atom_type,atompos0)  !(n,natom_type,atompos0)
        END DO

        IF (natoms.ne.tot_atom_number) THEN !natoms is initialized in <force_constants_init>
            WRITE(unit_log,*) "natoms read in lat_fc=",tot_atom_number
            WRITE(unit_log,*) "while natoms generated by fc_init=",natoms
            WRITE(unit_log,*) "check the number of shells in params.inp, "
            WRITE(unit_log,*) "it is probably too large if natoms in fc_init is larger than in lat_fc"
            WRITE(unit_log,*) "in which case the program will stop"
        END IF
!********************************************************** !continue my original 'reading' now
        ! legacy variable atompos(j, i) is matched here
        Do i=1,tot_atom_number

            READ(unit_number,*) every_atom(i)%label_number,every_atom(i)%x,every_atom(i)%y,every_atom(i)%z,&
            &                   every_atom(i)%type_tau,every_atom(i)%n1,every_atom(i)%n2,every_atom(i)%n3
            every_atom(i)%R = every_atom(i)%n1*trans_vec(:,1) + &
            &                 every_atom(i)%n2*trans_vec(:,2) + every_atom(i)%n3*trans_vec(:,3)
            every_atom(i)%tau = iatom(every_atom(i)%type_tau)%pos_tau

            atompos(1,i)=every_atom(i)%x;atompos(2,i)=every_atom(i)%y;atompos(3,i)=every_atom(i)%z
            iatomcell0(i)=every_atom(i)%type_tau
            iatomcell(1,i)=every_atom(i)%n1;iatomcell(2,i)=every_atom(i)%n2;iatomcell(3,i)=every_atom(i)%n3
            WRITE(unit_log,*) i,atompos(:,i),iatomcell0(i),iatomcell(:,i)
        End Do
!############################this part is to save all the 'different' {R} in my formalism############################
        !store all the 'different' R in another array: cell_vec(:,i)
        ALLOCATE(cell_vec(d,tot_atom_number)) !the most number of cell_vec(3,:) is when every atom is in different unit cell
        ALLOCATE(cell_save(d,tot_atom_number))

        cell_vec(:,1)=every_atom(1)%R
        every_atom(1)%type_R=1 !center cell starts at R=0
        j=2
        Do i=2,tot_atom_number
            check = every_atom(i)%R-every_atom(i-1)%R
inner:      DO direction=1,d
                IF(check(direction).ne.0) THEN
                    cell_vec(:,j)=every_atom(i)%R
                    every_atom(i)%type_R=j
                    j=j+1
                    exit inner
                ELSE
                    every_atom(i)%type_R=j-1
                END IF
            END DO inner

        END DO

        cell_save=cell_vec
        cell_number=j-2   !j=j+1, starts from j=1, so -2

        !****************************************************************************************************
        !cell_vec(xyz,index) is labeled by index, the max index corresponds to the center unit-cell which is not
        !taken into account in cell number N;
        !cell_number=SIZE(cell_vec,DIM=2)-1

        !****************************************************************************************************

        variational_parameters_size(1)=d*atom_number    ! internal coordinates shift
        variational_parameters_size(2)=d*d       ! strain

        CLOSE(unit_number)
        WRITE(unit_log,*) 'Coordinates read successfully'

!DO i=1,tot_atom_number
!    WRITE(*,*) 'atom i',i,'type_R',every_atom(i)%type_R
!END DO
!STOP
3 format(a,9(2x,f11.6))
!4 format(i6,3(2x,g11.5),4(3x,i4))
    END SUBROUTINE read_structure
!==========================================================
    SUBROUTINE decide_ksize(cell_vec,nk_c1,nk_c2,nk_c3)
    !! old subroutine by me, not used
        IMPLICIT NONE
        REAL(8),ALLOCATABLE,DIMENSION(:,:),INTENT(IN) :: cell_vec
        INTEGER,ALLOCATABLE,DIMENSION(:) :: n1,n2,n3
        INTEGER :: i,nk_c1,nk_c2,nk_c3

        ALLOCATE(n1(tot_atom_number),n2(tot_atom_number),n3(tot_atom_number))

        DO i=1,tot_atom_number
            n1(i)=every_atom(i)%n1
            n2(i)=every_atom(i)%n2
            n3(i)=every_atom(i)%n3
        END DO
        nk_c1=MAXLOC(n1,DIM=1)
        nk_c1=n1(nk_c1)
        nk_c2=MAXLOC(n2,DIM=1)
        nk_c2=n2(nk_c2)
        nk_c3=MAXLOC(n3,DIM=1)
        nk_c3=n3(nk_c3)

        DEALLOCATE(n1,n2,n3)
        WRITE(*,*) nk_c1,nk_c2,nk_c3
    END SUBROUTINE decide_ksize
!==========================================================
    SUBROUTINE check_latfc
     !!check if <latfc.dat> is read correctly
        IMPLICIT NONE
        INTEGER :: i,direction
        REAL(8) :: x,y,z
        OPEN(23,FILE='check_read_latfc.dat')
        WRITE(23,*) 'i,x,y,z,type_tau,n1,n2,n3'
        WRITE(23,*) SIZE(every_atom)
        DO i=1,SIZE(every_atom)
            x = every_atom(i)%tau(1) + every_atom(i)%R(1)
            y = every_atom(i)%tau(2) + every_atom(i)%R(2)
            z = every_atom(i)%tau(3) + every_atom(i)%R(3)
            WRITE(23,6) i,x,y,z,every_atom(i)%type_tau,every_atom(i)%n1,&
            &every_atom(i)%n2,every_atom(i)%n3
        END DO
6 format(i7,3(7x,f14.10),4(i7))
        CLOSE(23)
    END SUBROUTINE check_latfc

END MODULE Structure_info
