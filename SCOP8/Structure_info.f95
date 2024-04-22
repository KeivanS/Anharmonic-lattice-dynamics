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
 !===================================================================================
subroutine write_neighbors
!!legacy subroutine 
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
    !! legacy subroutine 
    !! for generating 'params.inp', currently not used
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
    !! major I/O subroutine, read 'structure.params'
    use io2 
    !NOTE:the uparams is causing trouble because it's declared both in [ios] & [io2] with different value
    implicit none
    integer i,counter,label,unit_params
    integer nouse1, nouse2, nouse3
    real(8) scal,junk
    character jjj*1

    unit_params=49
    open(unit_params,file='structure.params',status='old')

    !---NOTE: hard-coded params that used to be set in 'params.inp'
    maxneighbors = 90 !I just default set it to 90 for now 
    tolerance = 0.001; margin = 0.00001
    maxshells = 90 !MODIFY:27 this is added from the fc234 main program to here 

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
    read(unit_params,*) nshells(5,1:natoms0)       ! # of neigr-shells to use for each rank
    read(unit_params,*) nshells(6,1:natoms0)       ! # of neigr-shells to use for each rank
    read(unit_params,*) nshells(7,1:natoms0)       ! # of neigr-shells to use for each rank
    read(unit_params,*) nshells(8,1:natoms0)       ! # of neigr-shells to use for each rank
    
    !-------------screen output for check------------------
    write(*,*)'reading ',natoms0,' atoms in the primitive cell'
    write(*,*) svdcut,'    cutoff for smallest eigenvalue w to be included'
    write(*,*) tolerance,'   tolerance for equating two coordinates '
    write(*,*) include_fc,'  which ranks of FCs to include '
    write(*,*) nshells(1,:),'  how many shells to include for each rank '
    write(*,*) nshells(2,:),'  how many shells to include for each rank '
    write(*,*) nshells(3,:),'  how many shells to include for each rank '
    write(*,*) nshells(4,:),'  how many shells to include for each rank '
    write(*,*) itrans,irot,ihuang,enforce_inv,'  transl, rot and Huang invce, enforcing inv'

    latticeparameters(1:3) = latticeparameters(1:3)*scal
    if(.not.allocated(at_label)) allocate(at_label(natoms0))
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
!-----------------------------------------------------------------------------------
SUBROUTINE read_structure
    !!major I/O subroutine to read and record 'lat_fc.dat'
    !!and integrate a part of fcinit code(legacy) for atom-force info
    IMPLICIT NONE

    CHARACTER(LEN=90) :: line
    !CHARACTER(LEN=2) :: name
    INTEGER :: unit_number,unit_log,i,j,direction, maxshell!MARK: 09/04/2023
    REAL(8) :: check(d)
    TYPE(vector) :: vv
    LOGICAL :: condition
    CHARACTER path_output*29

    check=0
    unit_log=30
    !MODIFY: add into output path
    path_output = 'output/'
    OPEN(unit_log,file=trim(path_output)//'logfile.txt',status='unknown')

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

    r1%component=trans_vec(:,1)
    r2%component=trans_vec(:,2)
    r3%component=trans_vec(:,3)
    ! now generate the lattice ppties in real and recip spaces
    call make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
    call make_reciprocal_lattice(g1,g2,g3,rr1,rr2,rr3)
    call calculate_volume(r1,r2,r3,volume_r)
    call calculate_volume(g1,g2,g3,volume_g)
    ! screen output for check
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
    !--------start to read file---------
    READ(unit_number,*) line
    ! read(ufc,*) natoms0
    !----check compatibility with input file "params.inp" ----
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
    !----check further compatibility with input file "params.inp" ----
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
    READ(unit_number,*) line
    READ(unit_number,*) (nterms(j),j=1,4)
    READ(unit_number,*) line
    READ(unit_number,*) (ngroups(j),j=1,4)
    READ(unit_number,*) line
    READ(unit_number,*) maxshell, tot_atom_number 

    WRITE(unit_log,*) (include_fc(j),j=1,4)
    WRITE(unit_log,*) line
    WRITE(unit_log,*) (nterms(j),j=1,4)
    WRITE(unit_log,*) line
    WRITE(unit_log,*) (ngroups(j),j=1,4)
    WRITE(unit_log,*) line
    WRITE(unit_log,*)'Total Atom Number is', tot_atom_number

    ALLOCATE(every_atom(tot_atom_number))
    !this part has nothing to do with 'reading' and is directly copied from fc234 main program
    !NOTE: increase maxatoms if there are more atoms listed in lat_fc.dat
    ! the original comment reads:
    !"This does not work: fcinit must be called once or latparam(4:6) will be overwritten"
    maxatoms=6500; imaxat=1
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
    !NOTE: this part is to save all the 'unique' {R} in my formalism
    !store all the 'unique' R in another array: cell_vec(:,i)
    !the largest possible # of cell_vec(3,:) is when every atom is in different unit cell
    ALLOCATE(cell_vec(d,tot_atom_number)) 
    ALLOCATE(cell_save(d,tot_atom_number))

    cell_vec(:,1)=every_atom(1)%R
    every_atom(1)%type_R=1 !center cell starts at R=0
    j=2
    Do i=2,tot_atom_number
        check = every_atom(i)%R-every_atom(i-1)%R
inner:  DO direction=1,d
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
    
    !NOTE:
    !cell_vec(xyz,index) is labeled by index, the max index corresponds to the center unit-cell which is not
    !taken into account in cell number N;
    !cell_number=SIZE(cell_vec,DIM=2)-1

    variational_parameters_size(1)=d*atom_number    ! internal coordinates shift
    variational_parameters_size(2)=d*d       ! strain

    CLOSE(unit_number)
    WRITE(unit_log,*) 'Coordinates read successfully'

3 format(a,9(2x,f11.6))
END SUBROUTINE read_structure
!-------------------------------------------------
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

!==========================================================================
END MODULE Structure_info
