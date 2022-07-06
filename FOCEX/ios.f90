!===========================================================
 subroutine read_input
! reads the param.inp file containing info on the atom types masses
! and coordinates within the primitive cell
 use ios
 use params
 use lattice
 use force_constants_module
 use atoms_force_constants
 use svd_stuff
 implicit none
 integer i,counter,label
 real(8) scal

 open(uparams,file='params.inp',status='old')

 read(uparams,*) latticeparameters   ! (a,b,c,alpha,beta,gamma)
 read(uparams,*) primitivelattice    ! prim latt vectors in terms of conventional above
 read(uparams,*) scal  ! scale factor for lattparams is read from params.inp
 read(uparams,*) maxneighbors       ! # of neigr-shells to use for each rank
 read(uparams,*) include_fc    ! if=0 do not include this rank
! read(uparams,*) nshells , radius      ! # of neigr-shells to use for each rank
! read(uparams,*) maxterms      ! max # of terms for each rank
! read(uparams,*) maxtermsindep ! max # of independent terms for each rank
 read(uparams,*) itrans,irot,ihuang,enforce_inv   ! translational and rotational invce flags (include if=1)
 read(uparams,*) tolerance, margin  ! tolerance for equating coords, margin for eliminating FCs
 read(uparams,*) svdcut    ! cutoff for smallest "eigenvalue" w to be included
 read(uparams,*) fdfiles , verbose ! number of force-displacement data files

 read(uparams,*) natom_type   ! # of different elements present in prim cell
! allocate(natom(natom_type))
 call allocate_mass(natom_type)
 read(uparams,*) (mas(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) (atname(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) natoms0        ! # of atoms in primitive cell
 read(uparams,*) nshells(1,1:natoms0)
 read(uparams,*) nshells(2,1:natoms0)
 read(uparams,*) nshells(3,1:natoms0)
 read(uparams,*) nshells(4,1:natoms0)

 write(ulog,*) svdcut,'    cutoff for smallest eigenvalue w to be included'
 write(ulog,*) tolerance,'   tolerance for equating two coordinates '
 write(ulog,*) include_fc,'  which ranks of FCs to include '
 write(ulog,*) nshells(1,:),'  how many shells to include for each rank '
 write(ulog,*) nshells(2,:),'  how many shells to include for each rank '
 write(ulog,*) nshells(3,:),'  how many shells to include for each rank '
 write(ulog,*) nshells(4,:),'  how many shells to include for each rank '
 write(ulog,*) itrans,irot,ihuang,enforce_inv,'  transl, rot and Huang invce, enforcing inv'

 latticeparameters(1:3) = latticeparameters(1:3)*scal
 call allocate_primcell(natoms0)
! indices of atoms in primitive cell must be same order as in POSCAR
 counter = 0
 natom = 0
 do i=1,natoms0
! positions here must be D format in conventional cell for use by fcs_init
    read(uparams,*) label,atom_type(label),atompos0(:,label)
    if (i.ne.label) then
       print*,' positions must be sorted according to the labels ',i,label
       write(ulog,*)' READ_INPUT: positions must be sorted like labels ',i,label
       write(ulog,*)' Check atomic coordinates in primitive cell'
       stop
    endif
    atom0(i)%name = atname(atom_type(label))
    atom0(i)%at_type  = atom_type(label)
    atom0(i)%mass = mas(atom_type(label))
!   atom0(i)%nshells  ! actual no of neighbors shells < maxneighbors
!   atom0(i)%shells   ! what is inside each shell
! Make sure the label tau is the same as the atom index in the array atompos
! this is useless because natom(:) is read from POSCAR
    if (atom_type(label).eq.counter) then
        natom(atom_type(label)) = natom(atom_type(label))+1
    else
        natom(atom_type(label)) = 1
    endif
    counter = atom_type(label)
 enddo
 ! Trying to create temperature insertion in the params.inp
read(uparams,*) temperature
close(uparams)

 do i=1,4
!   if (nshells(i) .gt. maxneighbors) then
    if (maxval(nshells(i,:)) .gt. maxneighbors) then
       write(ulog,*)' READ_INPUT: nshell> maxneighbors for rank=', &
&                    i,maxval(nshells(i,:)),maxneighbors
       write(ulog,*)' Increase maxneighbors or choose a smaller # of shells'
       stop
    endif
 enddo

 end subroutine read_input
!===========================================================
 subroutine read_crystal(poscar)
! the following reads the POSCAR file that is used by VASP
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 implicit none
 character line*90, poscar*(*)
 integer i
 real(8) om,latt_const,a,b,c,om0
 type(vector) pos
 logical exst

 n1min = 100000; n2min = 100000; n3min = 100000
 n1max =-100000; n2max =-100000; n3max =-100000
! open (uposcar,file='POSCAR',status='old')
 inquire(file=poscar,exist=exst)
 if(exst) then
      open (uposcar,file=poscar,status='old')
 else
      write(ulog,*)' poscar file ',poscar,' does not exist; check your files location and run again'
      stop
 endif
 read(uposcar,'(a)') line
 write(ulog,*)' VASP job title'
 write(ulog,'(a)') line
! this is for the super-cell and has nothing to do with the one in params.inp
! which is for the primitive cell
 read(uposcar,*)lattice_parameter
 read(uposcar,*)r1
 read(uposcar,*)r2
 read(uposcar,*)r3
 if (lattice_parameter .lt.0) then  ! it is the -volume of the cell
    volume_r = -lattice_parameter  ! need to find the volume(r1,r2,r3) for scaling
    call calculate_volume(r1,r2,r3,om)
    latt_const = (volume_r/om)**(1./3.)
 else
    latt_const = lattice_parameter
 endif
 r1=latt_const*r1; r2=latt_const*r2; r3=latt_const*r3
 box(1)=length(r1); box(2)=length(r2); box(3)=length(r3)
 call calculate_volume(r1,r2,r3,om)
 call calculate_volume(r01,r02,r03,om0)
 call write_out(ulog,' Volume of supercell ',om)
 call write_out(ulog,' Volume of primicell ',om0)

 write(ulog,*)' box size in 1st direction=',box(1)
 write(ulog,*)' box size in 2nd direction=',box(2)
 write(ulog,*)' box size in 3rd direction=',box(3)

! number of atoms participating in MD,
 read(uposcar,*)(natom(i),i=1,natom_type)
 natom_super_cell = sum(natom)
 if (.not. (natom_super_cell*om0/(natoms0*om) .myeq. 1d0 ) ) then
    write(ulog,*)' supercell inconsistency; check input coordinates again'
    write(ulog,*)' natoms0, om0=',natoms0,om0
    write(ulog,*)' natom_sc, om=',natom_super_cell,om
    stop
 endif
 write(ulog,*)" number of different elements and number of atoms of each type"
 write(ulog,*) natom_type,(natom(i),i=1,natom_type)
 write(ulog,*)" total number of atoms in the supercell= ", natom_super_cell
 write(ulog,*) "translation vectors of the supercell are"
 write(ulog,4) r1
 write(ulog,4) r2
 write(ulog,4) r3

 if ( allocated(atom_sc)) then
    deallocate(atom_sc)
 endif
 allocate(atom_sc(natom_super_cell))

! read equilibrium positions from POSCAR file and write in logfile ------------
! read(uposcar,'(a)') line
 read(uposcar,*) line
 if (line(1:1).eq.'s' .or. line(1:1).eq.'S' ) then
!   read(uposcar,'(a)') line
    read(uposcar,*) line
 endif
 if (line(1:1).eq.'d' .or. line(1:1).eq.'D' ) then
    do i=1,natom_super_cell
       read(uposcar,*) a,b,c
       atom_sc(i)%equilibrium_pos = a*r1+b*r2+c*r3
    enddo
 elseif (line(1:1).eq.'c' .or. line(1:1).eq.'C' ) then
    do i=1,natom_super_cell
       read(uposcar,*) pos
       atom_sc(i)%equilibrium_pos = latt_const*pos  ! pos
    enddo
!   if(scale.ne.1)write(ulog,*)'WARNING: read cartesian positions are multiplied by scal!'
 else
    write(ulog,*)'POSCAR: positions are not in direct or cartesian coordinates'
    stop
 endif
! now in atom_sc%equilibrium_pos we have the cartesian coordinates.

 close(uposcar)

 write(ulog,*)' POSCAR read successfully and closed'


 write(ulog,*)'RECIPROCAL_LATTICE: '
 call make_reciprocal_lattice(r1,r2,r3,g1,g2,g3)
 call write_out(ulog,'om ',om)
 call write_out(ulog,'g1 ',g1)
 call write_out(ulog,'g2 ',g2)
 call write_out(ulog,'g3 ',g3)

4 format(9(2x,f19.9))
!6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))

 end subroutine read_crystal
!===========================================================
 subroutine check_input_poscar_consistency
! see if all atoms in the supercell can be obtained from the
! atoms in the input file using translations of the primitive lattice
 use atoms_force_constants
 use ios
 use geometry
 use lattice
 use params
 implicit none
 integer i,j,k, counter,ier,n1,n2,n3
 real(8) a1,a2,a3
 type(vector) shift0,vec
 logical matched

! first describe r_i in terms of r0_i, and make sure the lin comb is integer
 write(ulog,*)' ----------------------------------------------------------'
 write(ulog,*)' CHECKING the commensurability between primcell and supercell'
 call check(r1,a1,a2,a3,ier,g01,g02,g03)
 write(ulog,4)'checking r1 in units of r0s,ier=',a1,a2,a3,ier
 if (ier .eq. 1) stop
 call check(r2,a1,a2,a3,ier,g01,g02,g03)
 write(ulog,4)'checking r2 in units of r0s,ier=',a1,a2,a3,ier
 if (ier .eq. 1) stop
 call check(r3,a1,a2,a3,ier,g01,g02,g03)
 write(ulog,4)'checking r3 in units of r0s,ier=',a1,a2,a3,ier
 if (ier .eq. 1) stop
 write(ulog,*)' COMMENSURABILITY CHECKED ---------------------------------'
4 format(a,3(2x,f8.3),2x,i1)
!7 format(a,2(1x,i4),9(2x,f9.4))
!8 format(a,3(1x,i4),9(2x,f9.4))

! do i=1,natoms0
!    tau1 = atom0(i)%at_type
!    do j=1,natom_super_cell
!       if (tau1 .eq. atom_sc(j)%cell%tau) then
! see if they can be mapped onto one another by a primitive translation vector
!          vec = atom_sc(j)%equilibrium_pos - atom0(i)%equilibrium_pos
!          call check(vec,a1,a2,a3,ier,g01,g02,g03)
!          if (ier .eq. 1) then  ! There might be a constant shift
!              shift0 = (a1-floor(a1)) * r01 + (a2-floor(a2)) * r02 + (a3-floor(a3)) * r03
!              write(ulog,8)'i,j,tau,shift=',i,j,tau1,shift0
!          endif
!       endif
!    enddo
! enddo

! return

 write(ulog,*)'FINDING possible translation vectors, trying them on other atoms'

 checkloop: do i=1,natoms0  !,natom_super_cell
    shift0 = atom0(i)%equilibrium_pos - atom_sc(1)%equilibrium_pos

! try this shift see if all atoms in SC can be mapped to the PC by it.
    SC: do k = 1,natom_super_cell
       matched = .false.
       PRIM: do j = 1,natoms0  ! one of the prim-cell atoms has to match
          vec =  shift0 + atom_sc(k)%equilibrium_pos - atom0(j)%equilibrium_pos
          call check(vec,a1,a2,a3,ier,g01,g02,g03)
!         if (ier .eq. 0) matched = .true.
          if (ier .eq. 0) then
            matched = .true.
            cycle SC
          else
            cycle PRIM
          endif
       enddo PRIM

       if (.not. matched) then
          cycle checkloop  !exit SC
       endif

    enddo SC
    if (matched) exit checkloop
 enddo checkloop

 if (i.ge.1 .and. i.le.natoms0) then
    call write_out(ulog,'THE shift vector ', shift0)
    write(ulog,*)' it maps atom 1 of the SC to atom ',i,' of the prim cell'
 else
    write(ulog,*)' NO SHIFT VECTOR COULD BE DEFINED!, check your coordinates again'
    stop
 endif

! shift all the atoms in the primitive cell by -shift0 so that they fall
! on the atoms in the supercell
! write(ulog,*)'The shifted (cartesian) positions in the prim cell are now:'
! do i=1,natoms0
!    atom0(i)%equilibrium_pos = atom0(i)%equilibrium_pos-shift0
!    call write_out(ulog,'shifted position ', atom0(i)%equilibrium_pos)
! enddo
! write(ulog,*)'The shifted (cartesian) positions for ATOMPOS are now:'
! do i=1,natoms
!    atompos(:,i) = atompos(:,i)- shift0
!    call write_out(ulog,'shifted position ', atompos(:,i))
! enddo

! shifting vector has been identified. We can now identify all atoms of the SC
      SC2: do k = 1,natom_super_cell
         counter = 0
         matched = .false.
         PRIM2: do j = 1,natoms0  ! one of the prim-cell atoms has to match
            vec =  atom_sc(k)%equilibrium_pos - atom0(j)%equilibrium_pos+shift0
            call check(vec,a1,a2,a3,ier,g01,g02,g03)
            if(verbose) write(ulog,5)'ier,a1,a2,a3=',ier,a1,a2,a3
            if (ier.eq.0) then ! found a matching vector
                matched = .true.
                counter = counter + 1
                atom_sc(k)%cell%tau = iatomcell0(j)
 if(verbose)   write(ulog,*)' k in SC has type tau,j=',k,atom_sc(k)%cell%tau,j
                n1 = nint(a1)
                n2 = nint(a2)
                n3 = nint(a3)
                if (n1.gt.n1max) n1max = n1
                if (n2.gt.n2max) n2max = n2
                if (n3.gt.n3max) n3max = n3
                if (n1.lt.n1min) n1min = n1
                if (n2.lt.n2min) n2min = n2
                if (n3.lt.n3min) n3min = n3
                atom_sc(k)%cell%n(1) = n1
                atom_sc(k)%cell%n(2) = n2
                atom_sc(k)%cell%n(3) = n3
                atom_sc(k)%mass = atom0(j)%mass
            endif
         enddo PRIM2

         if (counter.eq. 0) then
            call write_out(ulog,' ATOM which can not be mapped to the prim cell ',i)
            write(ulog,*)' check your coordinates '
            write(ulog,3)atom_sc(i)%equilibrium_pos
            stop
         elseif ( counter .ge. 2) then
            write(ulog,*)' CHECK: more than one atom matched! check your coordinates'
            stop
         endif
       enddo SC2

 write(ulog,*)'*******************************************************'
 write(ulog,*)'SINCE THE PROGRAM WAS NOT STOPPED, MAPPING-CHECK PASSED!'
 write(ulog,*)'*******************************************************'

3 format(9(1x,g13.6))
5 format(a,i2,9(2x,g11.4))

 end subroutine check_input_poscar_consistency
!===========================================================
 subroutine check_input_poscar_consistency_new
! see if all atoms in the supercell can be obtained from the
! atoms in the input file using translations of the primitive lattice
 use atoms_force_constants
 use ios
 use geometry
 use lattice
 use params
 implicit none
 integer i,j,k, counter,ier,n1,n2,n3
 real(8) a1,a2,a3
 type(vector) shift0,vec
 logical matched

! first describe r_i in terms of r0_i, and make sure the lin comb is integer
 write(ulog,*)' ----------------------------------------------------------'
 write(ulog,*)' CHECKING the commensurability between primcell and supercell'
 call check(r1,a1,a2,a3,ier,g01,g02,g03)
 write(ulog,4)'checking r1 in units of r0s,ier=',a1,a2,a3,ier
 if (ier .eq. 1) stop
 nr1(1)=nint(a1); nr1(2)=nint(a2); nr1(3)=nint(a3)
 call check(r2,a1,a2,a3,ier,g01,g02,g03)
 write(ulog,4)'checking r2 in units of r0s,ier=',a1,a2,a3,ier
 if (ier .eq. 1) stop
 nr2(1)=nint(a1); nr2(2)=nint(a2); nr2(3)=nint(a3)
 call check(r3,a1,a2,a3,ier,g01,g02,g03)
 write(ulog,4)'checking r3 in units of r0s,ier=',a1,a2,a3,ier
 if (ier .eq. 1) stop
 nr3(1)=nint(a1); nr3(2)=nint(a2); nr3(3)=nint(a3)
 write(ulog,*)' COMMENSURABILITY CHECKED ---------------------------------'
 write(ulog,*)' nr1=',nr1
 write(ulog,*)' nr2=',nr2
 write(ulog,*)' nr3=',nr3

4 format(a,3(2x,f8.3),2x,i1)
7 format(a,2(1x,i4),9(2x,f9.4))
!8 format(a,3(1x,i4),9(2x,f9.4))

! do i=1,natoms0
!    tau1 = atom0(i)%at_type
!    do j=1,natom_super_cell
!       if (tau1 .eq. atom_sc(j)%cell%tau) then
! see if they can be mapped onto one another by a primitive translation vector
!          vec = atom_sc(j)%equilibrium_pos - atom0(i)%equilibrium_pos
!          call check(vec,a1,a2,a3,ier,g01,g02,g03)
!          if (ier .eq. 1) then  ! There might be a constant shift
!              shift0 = (a1-floor(a1)) * r01 + (a2-floor(a2)) * r02 + (a3-floor(a3)) * r03
!              write(ulog,8)'i,j,tau,shift=',i,j,tau1,shift0
!          endif
!       endif
!    enddo
! enddo

! return

 write(ulog,*)'FINDING possible translation vectors, trying them on other atoms'
 write(ulog,*)'g01,g02,g03='
 write(ulog,*)g01
 write(ulog,*)g02
 write(ulog,*)g03

 checkloop: do i=1,natoms0  !,natom_super_cell
    shift0 = atom0(i)%equilibrium_pos - atom_sc(1)%equilibrium_pos
    write(ulog,7)'atom in PRIMCELL:',i,atom0(i)%at_type,atom0(i)%equilibrium_pos
    write(ulog,3)'trying shift vector=',shift0

! try this shift see if all atoms in SC can be mapped to the PC by it.
    SC: do k = 1,natom_super_cell
       matched = .false.
       PRIM: do j = 1,natoms0  ! one of the prim-cell atoms has to match
          vec =  shift0 + atom_sc(k)%equilibrium_pos - atom0(j)%equilibrium_pos
          call check(vec,a1,a2,a3,ier,g01,g02,g03)
          if(verbose) write(ulog,6)'* err,k_sc,r_sc,a1,a2,a3=',ier,k,atom_sc(k)%equilibrium_pos,a1,a2,a3
          if (ier .eq. 0) then
             matched = .true.
!            write(ulog,*)' atom ',i,' in primitice cell matched with atom ',k,' in supercell'
             cycle SC !exit PRIM
          else
             cycle PRIM
          endif
       enddo PRIM

       if (.not. matched) then  ! wrong shift , try another shift
           cycle checkloop     !      exit SC
       endif

    enddo SC
    if (k.ge.natom_super_cell .and. matched) then  ! all atoms were matched
       exit checkloop
!   if (matched) exit checkloop
    endif
 enddo checkloop

 if (i.ge.1 .and. i.le.natoms0) then
    call write_out(ulog,'THE shift vector ', shift0)
    write(ulog,*)' it maps atom 1 of the SC to atom ',i,' of the prim cell'
 else
    write(ulog,*)' NO SHIFT VECTOR COULD BE DEFINED!, check your coordinates again;i=',i
    stop
 endif

! shift all the atoms in the primitive cell by -shift0 so that they fall
! on the atoms in the supercell
! write(ulog,*)'The shifted (cartesian) positions in the prim cell are now:'
! do i=1,natoms0
!    atom0(i)%equilibrium_pos = atom0(i)%equilibrium_pos-shift0
!    call write_out(ulog,'shifted position ', atom0(i)%equilibrium_pos)
! enddo
! write(ulog,*)'The shifted (cartesian) positions for ATOMPOS are now:'
! do i=1,natoms
!    atompos(:,i) = atompos(:,i)- shift0
!    call write_out(ulog,'shifted position ', atompos(:,i))
! enddo

! shifting vector has been identified. We can now identify all atoms of the SC
      SC2: do k = 1,natom_super_cell
         counter = 0
         matched = .false.
         PRIM2: do j = 1,natoms0  ! one of the prim-cell atoms has to match
            vec =  atom_sc(k)%equilibrium_pos - atom0(j)%equilibrium_pos+shift0
            call check(vec,a1,a2,a3,ier,g01,g02,g03)
            if(verbose) write(ulog,5)'ier,a1,a2,a3=',ier,a1,a2,a3
            if (ier.eq.0) then ! found a matching vector
                matched = .true.
                counter = counter + 1
                atom_sc(k)%cell%tau = iatomcell0(j)
  if(verbose)  write(ulog,*)' k in SC has type tau,j=',k,atom_sc(k)%cell%tau,j
                n1 = nint(a1)
                n2 = nint(a2)
                n3 = nint(a3)
                if (n1.gt.n1max) n1max = n1
                if (n2.gt.n2max) n2max = n2
                if (n3.gt.n3max) n3max = n3
                if (n1.lt.n1min) n1min = n1
                if (n2.lt.n2min) n2min = n2
                if (n3.lt.n3min) n3min = n3
                atom_sc(k)%cell%n(1) = n1
                atom_sc(k)%cell%n(2) = n2
                atom_sc(k)%cell%n(3) = n3
                atom_sc(k)%mass = atom0(j)%mass
                atom_sc(k)%at_type = atom0(j)%at_type
            endif
         enddo PRIM2

         if (counter.eq. 0) then
            call write_out(ulog,' ATOM which can not be mapped to the prim cell ',i)
            write(ulog,*)' check your coordinates '
            write(ulog,3)atom_sc(i)%equilibrium_pos
            stop
         elseif ( counter .ge. 2) then
            write(ulog,*)' CHECK: more than one atom matched! check your coordinates'
            stop
         endif
       enddo SC2

 write(ulog,*)'*******************************************************'
 write(ulog,*)'SINCE THE PROGRAM WAS NOT STOPPED, MAPPING-CHECK PASSED!'
 write(ulog,*)'*******************************************************'

3 format(9(1x,g13.6))
5 format(a,i2,9(2x,g11.4))
6 format(a,i2,1x,i5,9(2x,g11.4))

 open(173,file='poscar.xyz')
 write(173,*) natom_super_cell
 write(173,*) 'poscar.xyz to visualize'
 i=0
 do n1=1,natom_type
      do i=1,natom_super_cell
           if(atom_sc(i)%at_type .eq. n1) then
    if(n1.eq.1)  then
         write(173,9)'N ',atom_sc(i)%equilibrium_pos
    elseif(n1.eq.2) then
         write(173,9)'O ',atom_sc(i)%equilibrium_pos
    elseif(n1.eq.3) then
         write(173,9)'C ',atom_sc(i)%equilibrium_pos
    elseif(n1.eq.4) then
         write(173,9)'Ge ',atom_sc(i)%equilibrium_pos
    elseif(n1.eq.5) then
         write(173,9)'Si ',atom_sc(i)%equilibrium_pos
    endif

         endif
    enddo
 enddo
 close(173)
9 format(a,9(2x,f19.9))

 end subroutine check_input_poscar_consistency_new
!===========================================================
 subroutine read_force_position_data(outcar,frc_constr)
 use ios
 use atoms_force_constants
 use geometry
 use params
 use lattice
 use svd_stuff
 implicit none
 integer i,t,j,frc_constr,junk, k, xyzconfgs,en_count
! type(vector) v
! real(8) x1,x2,x3,x4
 character line*99,outcar*(*)
 logical found,exst
real(8) en_nconfigs
! open(utraj,file='OUTCAR',status='old')
 inquire(file=outcar,exist=exst)
 if(exst) then
      open(utraj,file=outcar,status='old')
 else
      write(ulog,*)' outcar file ',outcar,' does not exist; check your files location and run again'
      stop
 endif

 t=0
 do j=1,11000000
    read(utraj,'(a)',end=99)line
    call findword('POSITION',line,found)
!    if (line(1:8) .eq. "POSITION" ) then
    if (found) then
       t = t+1
    endif
 enddo
99 write(ulog,*)' reached the end of OUTCAR file; number of configurations= ',t
 nconfigs = t
 if(t.eq.0) then
    write(ulog,*)' the word POSITION was not found in OUTCAR file, check it!'
    stop
 endif
 i=(j-1)/t - 2  ! this is the number of atoms read from OUTCAR
 if (i .ne. natom_super_cell ) then
    write(ulog,*)' number of atoms read .ne. no of atoms in POSCAR file',i,natom_super_cell
    write(ulog,*)' # of read lines in OUTCAR is=',j
    write(ulog,*)' check your POSCAR and OUTCAR again '
    write(ulog,*)' make sure the # of atoms is the same in both files'
    write(ulog,*)' there should be no blank lines at the end of OUTCAR'
    stop
 endif

! allocates displ and force arrays
 if ( allocated(displ)) deallocate(displ)
 if ( allocated(force)) deallocate(force)
 if ( allocated(energy)) deallocate(energy)
 call allocate_pos(natom_super_cell,nconfigs)
! allocate(energy(nconfigs)) ! was this so I changed it to 3*nconfigs*natom_super_cell because the number of line dim_al=3*nconfigs*natom
en_count=3*nconfigs*natom_super_cell
allocate(energy(en_count)) 
! now get the FORCES from OUTCAR file
 rewind(utraj)
 t=0
 k=0
 do j=1,11000000  ! sum over snapshots
    read(utraj,'(a)',end=88)line
!    if (line(1:8) .eq. "POSITION" ) then
    call findword('POSITION',line,found)
    if (found) then
       read(utraj,*) junk,en_nconfigs
       t = t+1
       do i=1,natom_super_cell
           read(utraj,*) displ(1:3,i,t),force(1:3,i,t)
           do xyzconfgs=k+1,3*natom_super_cell+k
               energy(xyzconfgs)=en_nconfigs
           enddo
       enddo
       k=k+3*natom_super_cell
    endif
 enddo
88 write(ulog,*)' reached the end of OUTCAR file after steps= ',t
 write(ulog,*)' last line read was '
 write(ulog,'(a)') line
 close(utraj)
 if (t .ne. nconfigs) then
    write(ulog,*)'ERROR in reading the force file OUTCAR'
    write(ulog,*)'nconfigs, # of read steps=',nconfigs,t
    stop
 endif

 call write_out(ulog,' last force ',force(:,natom_super_cell,t))
 call write_out(ulog,' last coord ',displ(:,natom_super_cell,t))

! force_constraints = nconfigs *natom_super_cell*3
 frc_constr = nconfigs *natom_super_cell*3

! call pos_out_consistency

 end subroutine read_force_position_data
!===========================================================
 subroutine write_independent_fcs(sd)
 use svd_stuff
 use atoms_force_constants
 use ios
 use params
 implicit none
 integer i,g,cnt2,cnt
 real(8), intent(out):: sd(4)

   write(ulog,*)'******* Violation of translational invariance relations: ',itrans
   do i=1,transl_constraints
      write(ulog,*)i,dot_product(atransl(i,:),fcs(:))
   enddo

   write(ulog,*)'******* Violation of rotational invariance relations: ',irot
   do i=1,rot_constraints
      write(ulog,*)i,dot_product(arot(i,:),fcs(:)),brot(i)
   enddo

!  if(.not.allocated(ahuang)) allocate( ahuang(dim_al,ngr))
   if (ihuang.ne.0) then
      write(ulog,*)'******* Violation of Huang invariance relations: '
      do i=1,huang_constraints
         write(ulog,*)i,dot_product(ahuang(i,:),fcs(:))
      enddo
   endif

!   call check_huang

!return

   write(ulog,*)' group, fcs(group), sigma(group) for all SVD'
   do i=1,dim_ac
      write(ulog,4) i,fcs(i),sigma(i),sigma(i)/sqrt(1.*dim_al)
   enddo
   write(ulog,*) "rank, average error"

   cnt2=0;sd=0
   do i=1,4
   !   if(include_fc(i).eq.0) cycle
   !   cnt=0
   !   do g=1,map(i)%ngr ! Upto line 670 to 678 is a comment that I want to have here for the case include_fc(2) .eq. 2
      !  do k=1,map(i)%ntind(g) ! This was commented previously
   !         cnt=cnt+1
   !         cnt2=cnt2+1 !cnt
   !         sd(i)=sd(i)+sigma(cnt2)/sqrt(1.*dim_al)
      !  enddo  ! This was commented previously
   !   enddo
   if (include_fc(i) .eq. 2) cycle !This is where I tried to do the change to see if it works
   cnt=0
      do g=1,map(i)%ngr
         !do k=1,map(i)%ntind(g)
            cnt=cnt+1
            cnt2=cnt2+1 !cnt
            sd(i)=sd(i)+sigma(cnt2)/sqrt(1.*dim_al)
         !enddo
      enddo

      if(cnt .ne. 0 ) sd(i)=sd(i)/cnt
      write(ulog,4) i,sd(i)
   enddo
   write(ulog,5)' !==== Error summary for each rank',(sd(i),i=1,4)

! if  (enforce_inv .eq. 0) then
! else
!   write(ulog,*)' group, fcs(group), sigma(group) for non-eliminated SVD'
!
!!  fc1 = fc1-matmul(a11ia12,fcs)
!   do i=1,n_indep
!      write(ulog,4) i,fc1(i)
!   enddo
!   do i=1,dim_ac
!      write(ulog,4) i,fcs(i),sigma(i),sigma(i)/sqrt(1.*dim_al)
!   enddo
!   sigma = fcs
!   deallocate(fcs)
!   allocate(fcs(ngr))
!   fcs(1:n_indep)=fc1; fcs(1+n_indep:ngr) = sigma
! endif
4 format(i6,3(2x,g15.8))
5 format(a,9(2x,g11.4))

 end subroutine write_independent_fcs
!===========================================================
 subroutine write_independent_fcs_2(sd)
 use svd_stuff
 use atoms_force_constants
 use ios
 use params
 implicit none
 integer i,k,g,cnt2,cnt
 real(8) sd(4)

 if  (enforce_inv .eq. 0) then
   write(ulog,*)'******* Violation of translational invariance relations: ',itrans
   do i=1,transl_constraints
      write(ulog,*)i,dot_product(atransl(i,:),fcs(:))
   enddo
   write(ulog,*)'******* Violation of rotational invariance relations: ',irot
   do i=1,rot_constraints
      write(ulog,*)i,dot_product(arot(i,:),fcs(:)),brot(i)
   enddo
   write(ulog,*)'******* Violation of Huang invariance relations: '
   do i=1,huang_constraints
      write(ulog,*)i,dot_product(ahuang(i,:),fcs(:))
   enddo
   call check_huang

   write(ulog,*)' group, fcs(group), sigma(group) for all SVD'
   do i=1,dim_ac
      write(ulog,4) i,fcs(i),sigma(i),sigma(i)/sqrt(1.*dim_al)
   enddo

   cnt2=0;sd=0
   do i=1,4
      if(include_fc(i).eq.0) then
         cnt=0
         cycle
      endif
      cnt=0
      do g=1,map(i)%ngr
         do k=1,map(i)%ntind(g)
            cnt=cnt+1
            cnt2=cnt2+1 !cnt
            sd(i)=sd(i)+sigma(cnt2)/sqrt(1.*dim_al)
         enddo
      enddo
      if(cnt .ne. 0 ) sd(i)=sd(i)/cnt
      write(ulog,4) i,sd(i)
   enddo
   write(ulog,5)' !==== Error summary for each rank',(sd(i),i=1,4)

 else
   write(ulog,*)' group, fcs(group), sigma(group) for non-eliminated SVD'

!  fc1 = fc1-matmul(a11ia12,fcs)
   do i=1,n_indep
      write(ulog,4) i,fc1(i)
   enddo
   do i=1,dim_ac
      write(ulog,4) i,fcs(i),sigma(i),sigma(i)/sqrt(1.*dim_al)
   enddo
   sigma = fcs
   deallocate(fcs)
   allocate(fcs(ngr))
   fcs(1:n_indep)=fc1; fcs(1+n_indep:ngr) = sigma
 endif
4 format(i6,3(2x,g15.8))
5 format(a,9(2x,g11.4))

 end subroutine write_independent_fcs_2
!===========================================================
 subroutine write_neighbors
 use atoms_force_constants
 use force_constants_module
 use params
 use ios
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
!============================================================
 subroutine write_output_fc2
 use svd_stuff
 use ios
 use force_constants_module
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
 integer rnk,t,ti,i,res,j,rs !ni(4),nt(4)
 integer iat(4),ixyz(4),g,ng,term,term2,cnt2,frm,k
 real(8) rij,bunit,one,fcd,trace,dij
 character frmt*2,goh*48,ln*1,geh*47

 bunit = ryd/ab/ab
 one =1d0

!6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
!7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(2(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
!9 format(9(2x,f19.10))

! first write the crystal data
! nt=0; ni=0
! do i=1,4
!    if(map(i)%ngr.gt.0) then
!      nt(i)= sum(map(i)%nt(:))
!      ni(i)= sum(map(i)%ntind(:))
!      if ( (ni(i).ne.map(i)%ntotind) .or. (nt(i).ne.map(i)%ntot) ) then
!         write(ulog,*)'WRITE_OUTPUT_FC2: ERROR!!'
!         write(ulog,*)'ni(i).ne.map(i)%ntotind ',ni(i),map(i)%ntotind
!         write(ulog,*)'nt(i).ne.map(i)%ntot    ',nt(i),map(i)%ntot
!    endif
! enddo
! call write_lat_fc(ni,nt)  ! same as call write_lat_fc(map(:)%ntotind,map(:)%ntot)
! write(*,*) "FCS VALUE is: ", fcs
!----------------------------------------
!write(*,*) "map(3)%ngrB and map(4)%ngrB: ", map(3)%ngr, map(4)%ngr
 res = 0
 if (include_fc(2) .ne. 2) then
 do rnk=1,4
  if ( include_fc(rnk) .ne. 0 ) then
    frm=30+rnk
    write(ln,'(i1)')rnk
    goh='(a1,i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,2(g14.8,2x),f5.2)'
    geh='(i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g14.8,f8.4,2x,f9.5)'
    write(frmt,'(i2)')30+rnk
    write(ulog,*)' FOR RANK=',rnk,' format=',frmt
    write(*,*)' FOR RANK=',rnk,' format=',frmt
    write(ufc1-1+rnk,*)'# RANK ',rnk,' tensors :term,group,(iatom,ixyz)_2 d^nU/dx_{i,alpha}^n'

  ng=map(rnk)%ngr ! number of groups
  cnt2=0
  term = 0
  term2= 0
  do g=1,map(rnk)%ngr  ! index of a given group
       if(g.gt.1) cnt2=cnt2+map(rnk)%ntind(g-1)
       !if(g .gt. 1 .and. rnk .eq. 2 .and. include_fc(rnk) .eq. 2) cnt2=cnt2-map(rnk)%ngr
      !if (include_fc(2) .eq. 2) then
      !   cnt2=cnt2-map(rnk-1)%ntind(g-1)
      !endif
    ! write in the log and fcn_fit.dat: cnt2+ti is the position of indep_fc of that rank
    do ti=1,map(rnk)%ntind(g)  ! index of independent terms in that group g
       iat(1:rnk)  = map(rnk)%gr(g)%iatind (:,ti)
       ixyz(1:rnk) = map(rnk)%gr(g)%ixyzind(:,ti)
       term = term+1
       if (rnk.eq.2) then
         rij = length(atompos(:,iat(1))-atompos(:,iat(2)))
       else
         rij = 0
       endif
   !    write(*,*) "The value of MAPERR is: ", cnt2+ti, g, ti
   !    write(*,*) "The value of MAPERR is: ", map(rnk)%err(cnt2+ti)
       write(ulog,goh) map(rnk)%err(cnt2+ti),g,ti,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+cnt2+ti),fcs(res+cnt2+ti)/ryd*ab**rnk,rij
       write(ufit1-1+rnk,geh) ti,g,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+cnt2+ti),one,rij

    enddo

    ! write in the fcn.dat file
    do t=1,map(rnk)%nt(g)  ! index of dependent terms in that group g
       iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
       ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
       term2= term2+1
       fcd = 0
       ! must find the corresponding index of the indep term t <-> ti
       do ti=1,map(rnk)%ntind(g)
          ! this is the index of the indep FC coming in the A*FC=b matrix product
          fcd = fcd + fcs(res+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
       enddo
       if( abs(fcd) .gt. margin) then
          write(ufc1-1+rnk,geh)t,g, (iat(j),ixyz(j),j=1,rnk),fcd,one
       endif
    enddo
  enddo
!  res = res+ndindp(rnk)
  res = res+map(rnk)%ntotind
  endif
 enddo
endif
! Chaning for the case IF_INCLUDE_FC(2) IS EQUAL TO 2--VERY CRUDE NEED TO ALLOCATE TO PLAY WITH RANK HERE
if (include_fc(2) .eq. 2) then
   res=0;
   do rnk=3,4
      if ( include_fc(rnk) .ne. 0 ) then
        frm=30+rnk
        write(ln,'(i1)')rnk
        goh='(a1,i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,2(g14.8,2x),f5.2)'
        geh='(i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g14.8,f8.4,2x,f9.5)'
        write(frmt,'(i2)')30+rnk
        write(ulog,*)' FOR RANK=',rnk,' format=',frmt
        write(*,*)' FOR RANK=',rnk,' format=',frmt
        write(ufc1-1+rnk,*)'# RANK ',rnk,' tensors :term,group,(iatom,ixyz)_2 d^nU/dx_{i,alpha}^n'
    
      ng=map(rnk)%ngr ! number of groups
      cnt2=0
      term = 0
      term2= 0
      do g=1,map(rnk)%ngr  ! index of a given group
           if(g.gt.1) cnt2=cnt2+map(rnk)%ntind(g-1)
           !if(g .gt. 1 .and. rnk .eq. 2 .and. include_fc(rnk) .eq. 2) cnt2=cnt2-map(rnk)%ngr
          !if (include_fc(2) .eq. 2) then
          !   cnt2=cnt2-map(rnk-1)%ntind(g-1)
          !endif
        ! write in the log and fcn_fit.dat: cnt2+ti is the position of indep_fc of that rank
        do ti=1,map(rnk)%ntind(g)  ! index of independent terms in that group g
           iat(1:rnk)  = map(rnk)%gr(g)%iatind (:,ti)
           ixyz(1:rnk) = map(rnk)%gr(g)%ixyzind(:,ti)
           term = term+1
           if (rnk.eq.2) then
             rij = length(atompos(:,iat(1))-atompos(:,iat(2)))
           else
             rij = 0
           endif
   !        write(*,*) "The value of MAPERR is: ", cnt2+ti, g, ti
       !    write(*,*) "The value of MAPERR is: ", map(rnk)%err(cnt2+ti)
           write(ulog,goh) map(rnk)%err(cnt2+ti),g,ti,(iat(j),ixyz(j),j=1,rnk),  &
           &     fcs(res+cnt2+ti),fcs(res+cnt2+ti)/ryd*ab**rnk,rij
           write(ufit1-1+rnk,geh) ti,g,(iat(j),ixyz(j),j=1,rnk),  &
           &     fcs(res+cnt2+ti),one,rij
    
        enddo
    
        ! write in the fcn.dat file
        do t=1,map(rnk)%nt(g)  ! index of dependent terms in that group g
           iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
           ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
           term2= term2+1
           fcd = 0
           ! must find the corresponding index of the indep term t <-> ti
           do ti=1,map(rnk)%ntind(g)
              ! this is the index of the indep FC coming in the A*FC=b matrix product
              fcd = fcd + fcs(res+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
           enddo
           if( abs(fcd) .gt. margin) then
              write(ufc1-1+rnk,geh)t,g, (iat(j),ixyz(j),j=1,rnk),fcd,one
           endif
        enddo
      enddo
    !  res = res+ndindp(rnk)
      endif
      res = res+map(rnk)%ntotind
     enddo
endif

write(ulog,*)'******* Trace for the harmonic FCs ********'
 open(456,file='trace_fc.dat')

! write the trace of FC2
 rnk=2
  iloop: do i=1,natoms0
  jloop: do j=1,natoms
     rij = length(atompos(:,i)-atompos(:,j))
     trace=0
!     rs=ndindp(1)
     rs=map(1)%ntotind
!     rs=0   ! I HAVE JUST SWITCHED THIS TO ZERO TO CHECK IF THE VALUE OF RS IS THE ISSUE
   !  if ( include_fc(rnk) .eq. 2 ) then
   if ( include_fc(rnk) .ne. 2 ) then ! was if (include_fc(rnk) .ne. 0) then. I changed it to 2 because thats what we need to check for
        ng=map(rnk)%ngr ! number of groups
        cnt2=0
        term2= 0
  do g=1,map(rnk)%ngr  ! index of a given group
       if(g.gt.1) cnt2=cnt2+map(rnk)%ntind(g-1)
   !    write(*,*) "VALUE OF CNT2 AT START IS: ", cnt2
    do t=1,map(rnk)%nt(g)  ! index of independent terms in that group g
       iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
       ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
  !    i=iat(1) ; j=iat(2)
       if (iat(1).ne.i) cycle !iloop
       if (iat(2).ne.j) cycle ! jloop
!      write(*,*)'i,j,term2,al,be=',i,j,term2,ixyz(1),ixyz(2)
       fcd = 0
       ! must find the corresponding index of the indep term t <-> ti
       do ti=1,map(rnk)%ntind(g)
   !       write(*,*) 'res,cnt2,ti,res+cnt2+ti=',res,cnt2,ti,res+cnt2+ti
          ! this is the index of the indep FC coming in the A*FC=b matrix product
      !    fcd = fcd + fcs(rs+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
          fcd = fcd + fcs(cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
       enddo
       dij = length(atompos(:,iat(1))-atompos(:,iat(2)))
       if (ixyz(1).eq.ixyz(2)) then
          term2= term2+1
          trace = trace+ fcd
  !       write(*,*)'al,term2,trace=',ixyz(1),term2,trace
       endif
    enddo
  enddo
       if(trace.ne.0) then
           write(ulog,8) i,j,dij,trace
           write(456,8) i,j,dij,trace
       endif
 !     if (term2.ne.3) write(456,*)'#ERROR: there are ',term2,' terms for rij=',dij

   endif

!if ( include_fc(2) .eq. 2 ) then
!      ng=map(rnk)%ngr ! number of groups
!      cnt2=0
!      term2= 0
!      res=0   ! res=map(1)%ntotind
!      write(*,*) "THE SHAPE OF FCS IS: ", size(fcs)
!!do rnk=3, 4
!do g=1,map(rnk)%ngr  ! index of a given group
!     if(g.gt.1) cnt2=cnt2+map(rnk)%ntind(g-1)
!!     write(*,*) "VALUE OF CNT2 AT START IS: ", cnt2
!  do t=1,map(rnk)%nt(g)  ! index of independent terms in that group g
!     iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
!     ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
!!    i=iat(1) ; j=iat(2)
!     if (iat(1).ne.i) cycle !iloop
!     if (iat(2).ne.j) cycle ! jloop
!      write(*,*)'i,j,term2,al,be=',i,j,term2,ixyz(1),ixyz(2)
!     fcd = 0
!     ! must find the corresponding index of the indep term t <-> ti
!     do ti=1,map(rnk)%ntind(g)
!   !     write(*,*) 'res,cnt2,ti,res+cnt2+ti=',res,cnt2,ti,res+cnt2+ti
!        ! this is the index of the indep FC coming in the A*FC=b matrix product
!    !    fcd = fcd + fcs(rs+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
!   !     fcd = fcd + fcs(res+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
!   !     write(*,*) "Length of fcs: ", size(fcs)
!        if (cnt2 .eq. size(fcs)) then
!   !         write(*,*) "Entered this RES EQUAL SIZE CONDITION"
!            cnt2=cnt2-1
!            fcd = fcd + fcs(res+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
!        endif
!        write(*,*) "THE VALUE OF RES, CNT2, TI is: ", res, cnt2, ti
!        if (cnt2 .ne. size(fcs)) then
!            fcd = fcd + fcs(res+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
!        endif
!     enddo
!     dij = length(atompos(:,iat(1))-atompos(:,iat(2)))
!     if (ixyz(1).eq.ixyz(2)) then
!        term2= term2+1
!        trace = trace+ fcd
!!       write(*,*)'al,term2,trace=',ixyz(1),term2,trace
!     endif
!  enddo
!enddo
!     if(trace.ne.0) then
!         write(ulog,8) i,j,dij,trace
!         write(456,8) i,j,dij,trace
!     endif
!!     if (term2.ne.3) write(456,*)'#ERROR: there are ',term2,' terms for rij=',dij
!   res=res+map(rnk)%ntotind  
!!   enddo
!endif

enddo jloop
enddo iloop

 close(456)

write(ulog,*)'***************** END OF FC Trace ******************'


!----------------------------------------
! 125  format(a)

 if (res.ne.ngr) then
    write(ulog,*)'WRITE_OUTPUT: sum(nterms),ngr=',res,ngr
    write(ulog,*)'WRITE_OUTPUT: they should be equal!'
 endif
 end subroutine write_output_fc2
!============================================================
!*************************************************BIKASH************************************************************************
!*************************************************BIKASH************************************************************************
!*************************************************BIKASH************************************************************************
 


 
!*************************************************BIKASH************************************************************************
!*************************************************BIKASH************************************************************************
!*************************************************BIKASH************************************************************************
 subroutine read_output_fc2
   use svd_stuff
   use ios
   use force_constants_module
   use atoms_force_constants
   use params
   use lattice
   use constants
   implicit none
   integer rnk,t,ti,i,res,j,rs !ni(4),nt(4)
   integer iat(4),ixyz(4),g,ng,term,term2,cnt2,frm, gi
   real(8) rij,bunit,one,fcd,trace,dij
   character frmt*2,goh*48,ln*1,geh*47, line*99
   logical ex
   bunit = ryd/ab/ab
   one =1d0
  !write(*,*) "This call is taking..."
  !6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
  !7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
  8 format(2(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
  !9 format(9(2x,f19.10))
  
  ! first write the crystal data
  ! nt=0; ni=0
  ! do i=1,4
  !    if(map(i)%ngr.gt.0) then
  !      nt(i)= sum(map(i)%nt(:))
  !      ni(i)= sum(map(i)%ntind(:))
  !      if ( (ni(i).ne.map(i)%ntotind) .or. (nt(i).ne.map(i)%ntot) ) then
  !         write(ulog,*)'WRITE_OUTPUT_FC2: ERROR!!'
  !         write(ulog,*)'ni(i).ne.map(i)%ntotind ',ni(i),map(i)%ntotind
  !         write(ulog,*)'nt(i).ne.map(i)%ntot    ',nt(i),map(i)%ntot
  !    endif
  ! enddo
  ! call write_lat_fc(ni,nt)  ! same as call write_lat_fc(map(:)%ntotind,map(:)%ntot)
  
  !----------------------------------------
   res = 0
   do rnk=1,4
    if ( include_fc(rnk) .ne. 0 ) then
      frm=30+rnk
      write(ln,'(i1)')rnk
      goh='(a1,i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,2(g14.8,2x),f5.2)'
      geh='(i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g14.8,f8.4,2x,f9.5)'
      write(frmt,'(i2)')30+rnk
      write(ulog,*)' FOR RANK=',rnk,' format=',frmt
      write(*,*)' FOR RANK=',rnk,' format=',frmt
      inquire(file="fc2.dat",exist=ex)
      if (ex) then
         if ( rnk .eq. 2) then
            if ( include_fc(rnk) .eq. 2) then
               write(*,*) "This call is taking FC"
               open(431 ,file="fc2.dat")
               read(431,'(A)') line
               write(*,*) line
            endif
         endif 
      endif
     ! write(ufc1-1+rnk,*)'# RANK ',rnk,' tensors :term,group,(iatom,ixyz)_2 d^nU/dx_{i,alpha}^n'
  !write(*,*) "Outside if is valid..."
    ng=map(rnk)%ngr ! number of groups
   ! write(*,*) "The value of ng is: ", map(rnk)%ngr
    cnt2=0
    term = 0
    term2= 0
    do g=1,map(rnk)%ngr  ! index of a given group
         if(g.gt.1) cnt2=cnt2+map(rnk)%ntind(g-1)
 ! write(*,*) "Looping in g..."
      ! write in the log and fcn_fit.dat: cnt2+ti is the position of indep_fc of that rank
      do ti=1,map(rnk)%ntind(g)  ! index of independent terms in that group g
         iat(1:rnk)  = map(rnk)%gr(g)%iatind (:,ti)
         ixyz(1:rnk) = map(rnk)%gr(g)%ixyzind(:,ti)
         term = term+1
         if (rnk.eq.2) then
           rij = length(atompos(:,iat(1))-atompos(:,iat(2)))
         else
           rij = 0
         endif
      !   write(ulog,goh) map(rnk)%err(cnt2+ti),g,ti,(iat(j),ixyz(j),j=1,rnk),  &
      !   &     fcs(res+cnt2+ti),fcs(res+cnt2+ti)/ryd*ab**rnk,rij
      !   write(ufit1-1+rnk,geh) ti,g,(iat(j),ixyz(j),j=1,rnk),  &
      !   &     fcs(res+cnt2+ti),one,rij
      enddo
! write(*,*) "Value for include_fc(rnk) is: ", include_fc(rnk) 
      ! read in from the fcn.dat file
if ( rnk .eq. 2) then
   if ( include_fc(rnk) .eq. 2) then
      do t=1,map(rnk)%nt(g)  ! index of dependent terms in that group g
!         iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
!         ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
         term2= term2+1
!         fcd = 0
!         write(*,*) "This call is also taking..."
         ! must find the corresponding index of the indep term t <-> ti, g <-> gi
!         do ti=1,map(rnk)%ntind(g)
            ! this is the index of the indep FC coming in the A*FC=b matrix product
!            fcd = fcd + fcs(res+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
!         enddo
!         if( abs(fcd) .gt. margin) then
!            write(ufc1-1+rnk,geh)t,g, (iat(j),ixyz(j),j=1,rnk),fcd,one
            read(431,*) ti,gi, (iat(j),ixyz(j),j=1,rnk),fcd,one
!         endif
 !           write(*,*) "The reading is done..."
      enddo
   endif
endif
close(431)
    enddo
  !  res = res+ndindp(rnk)
    res = res+map(rnk)%ntotind
    endif
   enddo
  
!  write(ulog,*)'******* Trace for the harmonic FCs ********'
   open(456,file='trace_fc.dat')
  
  ! write the trace of FC2
   rnk=2
    iloop: do i=1,natoms0
    jloop: do j=1,natoms
       rij = length(atompos(:,i)-atompos(:,j))
       trace=0
  !     rs=ndindp(1)
       rs=map(1)%ntotind
       if ( include_fc(rnk) .ne. 0 ) then
          ng=map(rnk)%ngr ! number of groups
          cnt2=0
          term2= 0
    do g=1,map(rnk)%ngr  ! index of a given group
         if(g.gt.1) cnt2=cnt2+map(rnk)%ntind(g-1)
      do t=1,map(rnk)%nt(g)  ! index of independent terms in that group g
         iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
         ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
    !    i=iat(1) ; j=iat(2)
         if (iat(1).ne.i) cycle !iloop
         if (iat(2).ne.j) cycle ! jloop
  !      write(*,*)'i,j,term2,al,be=',i,j,term2,ixyz(1),ixyz(2)
         fcd = 0
         ! must find the corresponding index of the indep term t <-> ti
         do ti=1,map(rnk)%ntind(g)
            ! this is the index of the indep FC coming in the A*FC=b matrix product
            fcd = fcd + fcs(rs+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
         enddo
         dij = length(atompos(:,iat(1))-atompos(:,iat(2)))
         if (ixyz(1).eq.ixyz(2)) then
            term2= term2+1
            trace = trace+ fcd
    !       write(*,*)'al,term2,trace=',ixyz(1),term2,trace
         endif
      enddo
    enddo
         if(trace.ne.0) then
             write(ulog,8) i,j,dij,trace
             write(456,8) i,j,dij,trace
         endif
   !     if (term2.ne.3) write(456,*)'#ERROR: there are ',term2,' terms for rij=',dij
  
       endif
    enddo jloop
    enddo iloop
  
   close(456)
  
  !write(ulog,*)'***************** END OF FC Trace ******************'
  !----------------------------------------
  ! 125  format(a)
  
   !if (res.ne.ngr) then
   !   write(ulog,*)'WRITE_OUTPUT: sum(nterms),ngr=',res,ngr
   !   write(ulog,*)'WRITE_OUTPUT: they should be equal!'
   !endif
   end subroutine read_output_fc2
 !*************************************************BIKASH************************************************************************
 !*************************************************BIKASH************************************************************************
 !*************************************************BIKASH************************************************************************
 subroutine write_output_fc
 use svd_stuff
 use ios
 use force_constants_module
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
 integer rank,t,i,res,j,term2
 real(8) rij,bunit,one,dij,trace,fcd
 !logical ex

 one = 1d0
 bunit = ryd/ab/ab

1 format(i6,1x,i5,1(3x,(i4,1x,i1)),3x,g15.8,2x,f5.2)
2 format(i6,1x,i5,2(3x,(i4,1x,i1)),3x,g15.8,2x,f5.2)
3 format(i6,1x,i5,3(3x,(i4,1x,i1)),3x,g15.8,2x,f5.2)
4 format(i6,1x,i5,4(3x,(i4,1x,i1)),3x,g15.8,2x,f5.2)
!6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
!7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(2(2x,i3),9(2x,g13.6))
18 format(a,2(2x,i3),9(2x,g13.6))
!9 format(9(2x,f19.10))
!12 format(i6,1x,i4,2(3x,(i3,1x,i1)),' [',i2,',(',3(1x,i2),') ]',3x,g13.6,2(2x,f7.3))
31 format(a1,i6,1x,i5,1(3x,(i3,1x,i1)),3x,2(g15.8,2x),f5.2)
32 format(a1,i6,1x,i4,2(3x,(i3,1x,i1)),' [',i3,',(',3(1x,i2),') ]',3x,2(g13.6,1x),2(2x,f7.3))
33 format(a1,i6,1x,i5,3(3x,(i3,1x,i1)),3x,2(g15.8,2x),f5.2)
34 format(a1,i6,1x,i5,4(3x,(i3,1x,i1)),3x,2(g15.8,2x),f5.2)

! first write the crystal data
 call write_lat_fc(ngroups,nterms)

 write(ulog,*)'********** INDEPENDENT FORCE CONSTANTS ARE *********'
 write(ulog,*)' term, group(term) , (atom,xyz)_rank , FC(eV/A^n), FC(Hart/aB^n), amplitude '

!----------------------------------------
 res = 0
 rank=1
 if ( include_fc(rank) .eq. 1 ) then
 write(ufc1,*)' RANK 1 tensors :term,group,iatom,ixyz dU/dx_{i,alpha}'
 do t=1,nterms(rank)
       write(ufc1,1)t,igroup_1(t), &
& iatomterm_1(1,t),ixyzterm_1(1,t),  &
& fcs(res+igroup_1(t))*ampterm_1(t),one
!& fcs(res+igroup_1(t)),ampterm_1(t)
 enddo

 do i=1,ngroups(rank)
   t = map_1(i)
!   write(ulog,1) t,igroup_1(t), iatomterm_1(1,t),ixyzterm_1(1,t),  &
! & fcs(res+igroup_1(t)),ampterm_1(t)
   write(ulog,31) err_1(igroup_1(t)),t,igroup_1(t), iatomterm_1(1,t),ixyzterm_1(1,t),  &
& fcs(res+igroup_1(t)),fcs(res+igroup_1(t))/ryd*ab,ampterm_1(t)
   write(ufit1,1) t,igroup_1(t), iatomterm_1(1,t),ixyzterm_1(1,t),  &
& fcs(res+igroup_1(t))*ampterm_1(t),one
!& fcs(res+igroup_1(t)),ampterm_1(t)
 enddo
 res = res + igroup_1(nterms(rank))
 endif
!----------------------------------------

 rank=2
 if ( include_fc(rank) .eq. 1 ) then
 write(ufc2,*)'# RANK 2 tensors :term,group,(iatom,ixyz)_2 d2U/dx_{i,alpha}^2'
 do t=1,nterms(rank)
    if( abs(fcs(res+igroup_2(t))*ampterm_2(t)) .gt. margin) then
       write(ufc2,2)t,igroup_2(t), &
& iatomterm_2(1,t),ixyzterm_2(1,t),iatomterm_2(2,t),ixyzterm_2(2,t),  &
& fcs(res+igroup_2(t))*ampterm_2(t),one
!& fcs(res+igroup_2(t)),ampterm_2(t)
    endif
 enddo

 do i=1,ngroups(rank)
   t = map_2(i)
   rij = length(atompos(:,iatomterm_2(1,t))-atompos(:,iatomterm_2(2,t)))
!   write(ulog,12) t,igroup_2(t),  &
!& iatomterm_2(1,t),ixyzterm_2(1,t),iatomterm_2(2,t),ixyzterm_2(2,t),  &
!& iatomcell0(iatomterm_2(2,t)),iatomcell(:,iatomterm_2(2,t)), &
!& fcs(res+igroup_2(t)),rij,ampterm_2(t)
   write(ulog,32) err_2(igroup_2(t)),t,igroup_2(t),  &
& iatomterm_2(1,t),ixyzterm_2(1,t),iatomterm_2(2,t),ixyzterm_2(2,t),  &
& iatomcell0(iatomterm_2(2,t)),iatomcell(:,iatomterm_2(2,t)), &
& fcs(res+igroup_2(t)),fcs(res+igroup_2(t))/ryd*ab*ab,rij,ampterm_2(t)
!   write(ufit2,2) t,igroup_2(t),  &
!& iatomterm_2(1,t),ixyzterm_2(1,t),iatomterm_2(2,t),ixyzterm_2(2,t),  &
!& fcs(res+igroup_2(t)),ampterm_2(t)
! fc2, ampterm has been changed to fc2*ampterm and 1 (basically the last column
! does not count) so that the pure fc is read in the general case where ampterm
! is really a matrix
   write(ufit2,2) t,igroup_2(t),  &
& iatomterm_2(1,t),ixyzterm_2(1,t),iatomterm_2(2,t),ixyzterm_2(2,t),  &
& fcs(res+igroup_2(t))*ampterm_2(t),one
 enddo

!----------------------------------------
write(ulog,*)'******* Trace for the harmonic FCs ********'
 open(456,file='trace_fc.dat')

! write the trace of FC2
 rank=2
 if ( include_fc(rank) .ne. 0 ) then
  iloop: do i=1,natoms0
  jloop: do j=1,natoms
        rij = length(atompos(:,i)-atompos(:,j))
        trace=0
        term2= 0
!       do g=1,ngroups(rank)
!          t = map_2(g)
        do t=1,nterms(rank)
           dij = length(atompos(:,iatomterm_2(1,t))-atompos(:,iatomterm_2(2,t)))
           if (iatomterm_2(1,t).ne.i) cycle !iloop
           if (iatomterm_2(2,t).ne.j) cycle ! jloop
           if (rij.ne.dij ) then
              write(*,18)'rij ne dij ',i,j,rij,dij
              stop
           endif
           write(*,*)'i,j,term2,al,be=',i,j,term2,ixyzterm_2(1,t),ixyzterm_2(2,t)
           fcd = fcs(res+igroup_2(t))*ampterm_2(t)
          if (ixyzterm_2(1,t).eq.ixyzterm_2(2,t)) then
             term2= term2+1
             trace = trace+ fcd
  !          write(*,*)'al,term2,trace=',ixyz(1),term2,trace
          endif
        enddo
        if(trace.ne.0) then
           write(ulog,8) i,j,rij,trace
           write(456,8) i,j,rij,trace
        endif
 !      if (term2.ne.3) write(456,*)'#ERROR: there are ',term2,' terms for rij=',dij

  enddo jloop
  enddo iloop
 endif

 close(456)

write(ulog,*)'***************** END OF FC Trace ******************'
!----------------------------------------

 res = res + igroup_2(nterms(rank))
 endif

 rank=3
 if ( include_fc(rank) .eq. 1 ) then
 write(ufc3,*)' RANK 3 tensors :term,group, (iatom,ixyz)_3 d3U/dx_{i,alpha}^3'
 do t=1,nterms(rank)
    if( abs(fcs(res+igroup_3(t))*ampterm_3(t)) .gt. margin) then
       write(ufc3,3)t,igroup_3(t), &
& iatomterm_3(1,t),ixyzterm_3(1,t),iatomterm_3(2,t),ixyzterm_3(2,t), &
& iatomterm_3(3,t),ixyzterm_3(3,t),  &
& fcs(res+igroup_3(t))*ampterm_3(t),one
!& fcs(res+igroup_3(t)),ampterm_3(t)
    endif
 enddo

 do i=1,ngroups(rank)
   t = map_3(i)
!   write(ulog,3) t,igroup_3(t),  &
!& iatomterm_3(1,t),ixyzterm_3(1,t),iatomterm_3(2,t),ixyzterm_3(2,t), &
!& iatomterm_3(3,t),ixyzterm_3(3,t),  &
!& fcs(res+igroup_3(t)),ampterm_3(t)
   write(ulog,33) err_3(igroup_3(t)),t,igroup_3(t),  &
& iatomterm_3(1,t),ixyzterm_3(1,t),iatomterm_3(2,t),ixyzterm_3(2,t), &
& iatomterm_3(3,t),ixyzterm_3(3,t),  &
& fcs(res+igroup_3(t)),fcs(res+igroup_3(t))/ryd*ab*ab*ab,ampterm_3(t)
   write(ufit3,3) t,igroup_3(t),  &
& iatomterm_3(1,t),ixyzterm_3(1,t),iatomterm_3(2,t),ixyzterm_3(2,t), &
& iatomterm_3(3,t),ixyzterm_3(3,t),  &
& fcs(res+igroup_3(t))*ampterm_3(t),one
!& fcs(res+igroup_3(t)),ampterm_3(t)
 enddo
 res = res + igroup_3(nterms(rank))
 endif
!----------------------------------------

 rank=4
 if ( include_fc(rank) .eq. 1 ) then
 write(ufc4,*)' RANK 4 tensors :term, group, (iatom,ixyz)_4 d4U/dx_{i,alpha}^4'
 do t=1,nterms(rank)
    if( abs(fcs(res+igroup_4(t))*ampterm_4(t)) .gt. margin) then
       write(ufc4,4)t,igroup_4(t), &
& iatomterm_4(1,t),ixyzterm_4(1,t),iatomterm_4(2,t),ixyzterm_4(2,t), &
& iatomterm_4(3,t),ixyzterm_4(3,t),iatomterm_4(4,t),ixyzterm_4(4,t), &
& fcs(res+igroup_4(t))*ampterm_4(t),one
!& fcs(res+igroup_4(t)),ampterm_4(t)
    endif
 enddo

 do i=1,ngroups(rank)
   t = map_4(i)
!   write(ulog,4) t,igroup_4(t),  &
!& iatomterm_4(1,t),ixyzterm_4(1,t),iatomterm_4(2,t),ixyzterm_4(2,t), &
!& iatomterm_4(3,t),ixyzterm_4(3,t),iatomterm_4(4,t),ixyzterm_4(4,t), &
!& fcs(res+igroup_4(t)),ampterm_4(t)
   write(ulog,34) err_4(igroup_4(t)),t,igroup_4(t),  &
& iatomterm_4(1,t),ixyzterm_4(1,t),iatomterm_4(2,t),ixyzterm_4(2,t), &
& iatomterm_4(3,t),ixyzterm_4(3,t),iatomterm_4(4,t),ixyzterm_4(4,t), &
& fcs(res+igroup_4(t)),fcs(res+igroup_4(t))/ryd*ab*ab*ab*ab,ampterm_4(t)
   write(ufit4,4) t,igroup_4(t),  &
& iatomterm_4(1,t),ixyzterm_4(1,t),iatomterm_4(2,t),ixyzterm_4(2,t), &
& iatomterm_4(3,t),ixyzterm_4(3,t),iatomterm_4(4,t),ixyzterm_4(4,t), &
& fcs(res+igroup_4(t))*ampterm_4(t),one
!& fcs(res+igroup_4(t)),ampterm_4(t)
 enddo
 res = res + igroup_4(nterms(rank))
 endif
!----------------------------------------
 if (res.ne.ngr) then
    write(ulog,*)'WRITE_OUTPUT: sum(nterms),ngr=',res,ngr
    write(ulog,*)'WRITE_OUTPUT: they should be equal!'
 endif
 end subroutine write_output_fc
!============================================================
 subroutine read_fcs_2(rank)   ! was this before subroutine read_fcs_2(iunit,rank)
 use svd_stuff
 use ios
 use force_constants_module
 use atoms_force_constants
 use params
 implicit none
 integer rank,iunit,t,res,i, a, b, c, d, e, f,reason !cnt, k, ti, j, g, l, newfile, rnk
 real amp, rijs !, gh
 !real(8), allocatable:: fc_ind(:)
 !character line*99
 logical ex

!newfile=431
!open(newfile,file="fc2.dat")
! read(newfile,'(a)') line
 !write(*,*) line
 t = 0
if ( rank .eq. 1) then

 res = 0
 do i=1,nterms(rank)
       read(iunit,*,err=91)t,igroup_1(t), &
& iatomterm_1(1,t),ixyzterm_1(1,t),  &
& fcs_1(igroup_1(t)),ampterm_1(t)
 enddo
 return
91 write(ulog,*)'READ_FCS: rank=',rank,' reached end of file after t=',t
 if (t.eq.0) then
    nterms(1) = 1
    iatomterm_1(1:1,1:1) = 1
    ixyzterm_1 (1:1,1:1) = 1
    igroup_1 (1:1) = 0
    ampterm_1(1:1) = 0
    include_fc(rank) = 0       ! if nothing, then exclude from fitting
 endif
!----------------------------------------
!write(*,*) "The value of map(2)%ngr is: ", map(2)%ngr
 elseif ( rank .eq. 2) then
!write(*,*) "Entering into INCLUDE_FC CONDITION 2: "
 if ( include_fc(rank) .eq. 2 ) then
!   write(*,*) "The value of map(2)%ngr and rank is: ", map(2)%ngr, rank
   allocate(fc_ind(map(2)%ngr))
   inquire ( file="fc2_fit.dat", exist=ex)
   if (ex) then
     open(473,file="fc2_fit.dat")
        do i=1, map(2)%ngr
           read(473,*,iostat=reason) a, b, c, d, e, f, fc_ind(i), amp, rijs
           if (reason > 0) then
              write(*,*) "something is wrong"
           elseif (reason < 0) then
              write(*,*) "END OF FILE REACHED"
           else
              write(*,*) a,b,c,d,e,f,fc_ind(i),amp,rijs
           endif
        enddo
     close(473)
  endif
 endif

 return

!elseif ( rank .eq. 2) then

! res =  igroup_1(nterms(1))
! do i=1,nterms(rank)
!       read(iunit,*,err=92)t,igroup_2(t), &
!& iatomterm_2(1,t),ixyzterm_2(1,t),iatomterm_2(2,t),ixyzterm_2(2,t),  &
!& fcs_2(igroup_2(t)),ampterm_2(t)
! enddo
! return
!92 write(ulog,*)'READ_FCS: rank=',rank,' reached end of file after t=',t
!----------------------------------------

elseif ( rank .eq. 3) then

 res =  igroup_1(nterms(1)) + igroup_2(nterms(2))
 do i=1,nterms(rank)
       read(iunit,*,err=93)t,igroup_3(t), &
& iatomterm_3(1,t),ixyzterm_3(1,t),iatomterm_3(2,t),ixyzterm_3(2,t), &
& iatomterm_3(3,t),ixyzterm_3(3,t),  &
& fcs_3(igroup_3(t)),ampterm_3(t)
 enddo
 return
93 write(ulog,*)'READ_FCS: rank=',rank,' reached end of file after t=',t
!----------------------------------------

elseif ( rank .eq. 4) then

 res =  igroup_1(nterms(1)) + igroup_2(nterms(2)) + igroup_3(nterms(3))
 do i=1,nterms(rank)
       read(iunit,*,err=94)t,igroup_4(t), &
& iatomterm_4(1,t),ixyzterm_4(1,t),iatomterm_4(2,t),ixyzterm_4(2,t), &
& iatomterm_4(3,t),ixyzterm_4(3,t),iatomterm_4(4,t),ixyzterm_4(4,t), &
& fcs_4(igroup_4(t)),ampterm_4(t)
 enddo
 return
94 write(ulog,*)'READ_FCS: rank=',rank,' reached end of file after t=',t
else

 write(ulog,*)' READ_FCS: rank must be from 1 to 4, not ',rank

endif

write(ulog,*)'READ_FCS: error!! should not have gotten here!'
if ( rank .ne. 1 .or. rank .ne. 2 .or. rank .ne. 3 .or. rank .ne. 4) then
stop  ! There was stop here
endif
 end subroutine read_fcs_2
!=================================================================
 subroutine calculate_and_write_displacements
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 implicit none
! character line*90
 integer i,j,t,step
 real(8) dc(3),dr(3)

 write(ulog,*)' t , particle#, cartesian disp & forces u1,u2,u3,f1,f2,f3'
! get displacements from positions
  do t=1,nconfigs
     do i=1,natom_super_cell
! first get direct coordinates, then take the distance using pbc, then
! transfrom back to cartesian coordinates
        dc(1) = displ(1,i,t) - atom_sc(i)%equilibrium_pos%x
        dc(2) = displ(2,i,t) - atom_sc(i)%equilibrium_pos%y
        dc(3) = displ(3,i,t) - atom_sc(i)%equilibrium_pos%z
        call cart_to_direct_aa(dc,dr)
        dr(1) = dr(1) - anint(dr(1))
        dr(2) = dr(2) - anint(dr(2))
        dr(3) = dr(3) - anint(dr(3))
        call direct_to_cart_aa(dr,displ(:,i,t))
     enddo
  enddo
  step=nconfigs-1
  if(verbose) step=1
  do t=1,nconfigs ,step           ! write first and last displ-force data
     do i=1,natom_super_cell
        write(ulog,6)t,i,(displ(j,i,t),j=1,3),(force(j,i,t),j=1,3)
     enddo
  enddo

 6 format(2(i5),3(1x,f10.6),3x,3(1x,g12.5))

 end subroutine calculate_and_write_displacements
!=================================================================
 subroutine identify_atoms_in_supercell
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 implicit none
! character line*90
 integer i,n1,n2,n3
 real(8) a(3),a1,a2,a3

! check mapping and assign label and address to each atom in the supercell
! this seems to be working for now....
 call check_input_poscar_consistency_new
 write(ulog,*)' n1min, n1max =',n1min,n1max
 write(ulog,*)' n2min, n2max =',n2min,n2max
 write(ulog,*)' n3min, n3max =',n3min,n3max
 write(ulog,*)' ASSIGNMENTS WERE DONE FOR ATOMS IN THE SUPER CELL '
 do i=1,natom_super_cell
!   write(ulog,*) '-----------------------  ATOM NUMBER, tau =',i,atom_sc(i)%cell%tau
!   call write_out(ulog,'      type ', atom_sc(i)%atom_type)
!   call write_out(ulog,'       tau ', atom_sc(i)%cell%tau)
!   call write_out(ulog,'      mass ', atom_sc(i)%mass)
!   call write_out(ulog,'     label ', atom_sc(i)%label_in_unit_cell)
    a1 = atom_sc(i)%equilibrium_pos .dot. g01
    a2 = atom_sc(i)%equilibrium_pos .dot. g02
    a3 = atom_sc(i)%equilibrium_pos .dot. g03
    n1 = nint(a1)
    n2 = nint(a2)
    n3 = nint(a3)
    write(ulog,8)'i,tau,type,n,r=',i,atom_sc(i)%cell%tau,atom_sc(i)%at_type,atom_sc(i)%cell%n,atom_sc(i)%equilibrium_pos
!   write(ulog,5)'   address is= ', atom_sc(i)%cell%n
!   call write_out(ulog,' eqlb posn ', atom_sc(i)%equilibrium_pos)
! if( atom_sc(i)%cell%n(1) .ne. n1) write(ulog,*)' address no match n1=',n1,a1
! if( atom_sc(i)%cell%n(2) .ne. n2) write(ulog,*)' address no match n2=',n2,a2
! if( atom_sc(i)%cell%n(3) .ne. n3) write(ulog,*)' address no match n3=',n3,a3
! this not needed as ai need not be integers in case the prim cell is not bravais
 enddo

 close(uposcar)
!3 format(i5,9(2x,g13.6))
!4 format(9(2x,f10.4))
!5 format(a,3(i6,2x))
!7 format(a,9(1x,f10.4))
8 format(a,i4,i3,i3,'(',3(i2,','),')',3(1x,f10.5))
 write(ulog,*)' POSCAR read successfully and closed'

 end subroutine identify_atoms_in_supercell
!=================================================================
      subroutine ustring(m,lineout,nrank,iatom,ixyz)
      implicit none
      integer i,m,n,nrank,iatom(nrank),ixyz(nrank)
      character lineout*80,xyz(3)
      data xyz/'x','y','z'/
!     write(*,*)' Entering ustring with rank=',nrank
      lineout(m+1:m+1)='d'
      m=m+1
      write(lineout(m+1:m+1),'(i1)')nrank
      m=m+1
      lineout(m+1:m+2)='U/'
      m=m+2
      do i=1,nrank
        lineout(m+1:m+2)='d'//xyz(ixyz(i))
        m=m+2
        n=iatom(i)
        if(n.lt.10)then
          write(lineout(m+1:m+1),'(i1)')n
          m=m+1
        else if(n.lt.100)then
          write(lineout(m+1:m+2),'(i2)')n
          m=m+2
        else if(n.lt.1000)then
          write(lineout(m+1:m+3),'(i3)')n
          m=m+3
        else
          write(lineout(m+1:m+4),'(i4)')n
          m=m+4
        endif
      enddo
!     write(*,*)' Exiting ustring'

      end subroutine ustring
!============================================================
 subroutine pos_out_consistency
 use geometry
 use ios
 use params
 use lattice
 use atoms_force_constants
 implicit none
 integer i
 real(8) dr(3),dc(3)
! checking consistency between POSCAR:atom_sc%equilibrium_pos and OUTCAR:displ(:,:,1)

 do i=1,natom_super_cell

    dc = atom_sc(i)%equilibrium_pos - displ(:,i,1)
    call cart_to_direct_aa(dc,dr)
    dr(1) = dr(1) - anint(dr(1))
    dr(2) = dr(2) - anint(dr(2))
    dr(3) = dr(3) - anint(dr(3))
    call direct_to_cart_aa(dr,dc)

!   if (.not.( dc .myeq. 0d0 ) ) then
    if ( abs(dc(1)).gt.tolerance .or. abs(dc(2)).gt.tolerance  &
    &                            .or. abs(dc(3)).gt.tolerance ) then
!    if ( dc(1).gt.0.1 .or. dc(2).gt.0.1 .or. dc(3).gt.0.1 ) then
       write(ulog,*)' atom # ',i,' in POSCAR and OUTCAR are different'
       write(ulog,*)' POSCAR:atom_SC=',atom_sc(i)%equilibrium_pos
       write(ulog,*)' OUTCAR: displ =',displ(:,i,1)
       write(ulog,*)' check your input files '
       write(*,*)' atom # ',i,' in POSCAR and OUTCAR are different'
       write(*,*)' check your input files '
       stop
    endif

 enddo

 end subroutine pos_out_consistency
!============================================================
 subroutine write_correspondance
 use params
 use atoms_force_constants
 use svd_stuff
 use ios
 use geometry
 use force_constants_module
 use lattice
 implicit none
 integer t,ired,nat,iatom,jatom,j
 integer al,be,taui,tauj,ni(3),nj(3),nsi(3),nsj(3)
 real(8) rij

 write(ulog,*)' HARMONIC FCs: correspondance between a given FC and all pairs '
 write(ulog,*)' Basically we describe the harmonic terms in the Taylor expansion'
 write(ulog,*)' Which atom pairs in the supercell are bridged by a FC element '
 write(ulog,*)' Below, ni refers to the translations of the supercell'
 write(ulog,*)'  i,   taui (nsi),  alpha ;     j,   tauj (nsj),   beta  :   term group  rij '
 write(ucor,*)'# i,   taui  (ni),  alpha ;     j,   tauj  (nj),   beta  :   term group  rij '

! do nat=1,natom_super_cell
  do nat=1,natoms0
  do al=1,3             ! this is how we order each line

!    taui = atom_sc(nat)%cell%tau
!    ni   = atom_sc(nat)%cell%n
     taui = iatomcell0(nat)
     ni(:)= iatomcell(:,nat)
     call findatom_sc(ni,nat,iatom)  ! now iatom is the supercell index
     nsi  = nint(matmul(r0g,dfloat(ni)))
     tloop2: do t=1,nterms(2)  ! ineq. terms and identify neighbors
           if ( (taui .eq. iatomterm_2(1,t)) .and.  &  !tau belongs to primcl
&               ( al  .eq. ixyzterm_2(1,t)) ) then
! this is the (tau,al) corresponding to that t, now what is its neighbor j?
!  (n1,tau1)=(iatomcell(:,iatomterm(2,t)),iatomcell0(iatomterm(2,t)))
               j    = iatomterm_2(2,t)
               tauj = iatomcell0(iatomterm_2(2,t))
               nj(:)= iatomcell(:,iatomterm_2(2,t)) + ni(:) ! translate by ni
               nsj  = nint(matmul(r0g,dfloat(nj)))
               be   = ixyzterm_2(2,t)
! Identify neighbor j within the SCell, find its displacement and add to ared
               call findatom_sc(nj,tauj,jatom)
               if (jatom.eq.0) then
        write(ulog,4)'WRITE_CORRESPONDANCE: jatom not found: tau,n ',tauj,nj
                  write(ulog,4)'for term ',t
                  write(ulog,4)'atom,xyz ',nat,al
                  stop
               endif
! if j and j+R(supercell) are in the same group but with opposite ampterms
! then their sum is cancelled in the force as they both have the same
! displacement. This will lead to errors in evaluation of the FCs if all
! terms in the group cancel in this way, and the
! corresponding group in FC2 will be evaluated to be zero. That will also
! affect the ASR and produce a violation of the sum rules.
               ired = igroup_2(t)  ! this is the corresponding index of ared
               rij = length(atompos(:,nat)-atompos(:,j))
             write(ucor,5)iatom,taui,ni,al,jatom,tauj,nj,be,t,ired,rij,ampterm_2(t)
             write(ulog,5)iatom,taui,nsi,al,jatom,tauj,nsj,be,t,ired,rij,ampterm_2(t)

           endif
     enddo tloop2

  enddo
  enddo

  write(ulog,*)' Correspondance of harmonic FCs and pairs of atoms established'

4 format(a,4(1x,i6))
5 format(2(i4,1x,' [ ',i4,' (',i2,',',i2,',',i2,') ] ',i1,2x),' : ',i6,1x,i4,2x,f7.3,2x,f7.3)

 end subroutine write_correspondance
!============================================================
 subroutine write_correspondance2
 use params
 use atoms_force_constants
 use svd_stuff
 use ios
 use geometry
 use force_constants_module
 use lattice
 implicit none
 integer t,nat,iatom,jatom,l,rnk
 integer al,be,taui,tauj,ni(3),nj(3),nsi(3),nsj(3)
 real(8) rij

 rnk = 2
 write(ulog,*)' HARMONIC FCs: correspondance between a given FC and all pairs '
 write(ulog,*)' Basically we describe the harmonic terms in the Taylor expansion'
 write(ulog,*)' Which atom pairs in the supercell are bridged by a FC element '
 write(ulog,*)' Below, ni refers to the translations of the supercell'
 write(ulog,*)'  i,    taui (nsi)        ;   j,    tauj (nsj)          :   term  group   rij '
 write(ucor,*)'# i,    taui (ni),  alpha ;   j,    tauj  (nj),   beta  :   term  group   rij '

! do nat=1,natom_super_cell
  do nat=1,natoms0
  do al=1,3             ! this is how we order each line

!    taui = atom_sc(nat)%cell%tau
!    ni   = atom_sc(nat)%cell%n
     taui = iatomcell0(nat)
     ni(:)= iatomcell(:,nat)
     call findatom_sc(ni,nat,iatom)  ! now iatom is the supercell index
     nsi  = nint(matmul(r0g,dfloat(ni)))

     do l=1,map(rnk)%ngr  ! sum over groups
!           if(l.gt.1) cnt2 = cnt2 + map(rnk)%ntind(l-1)
     do t=1,map(rnk)%nt(l) ! sum over all terms in that group
           if ( taui.eq. map(rnk)%gr(l)%iat(1,t) .and.  &
           &    al  .eq. map(rnk)%gr(l)%ixyz(1,t) ) then
              tauj = iatomcell0(map(rnk)%gr(l)%iat(2,t))
              nj(:)= iatomcell(:,map(rnk)%gr(l)%iat(2,t))+ni(:) ! trnslte by ni
              nsj  = nint(matmul(r0g,dfloat(nj)))
              be   = map(rnk)%gr(l)%ixyz(2,t)
! Identify neighbor j within the SCell,
              call findatom_sc(nj,tauj,jatom)
              if (jatom.eq.0) then
                  write(ulog,4)'WRITE_CORRESPONDANCE2:jatom not found: tau,n ',tauj,nj
                  write(ulog,4)'for rank,term ',rnk,t
                  write(ulog,4)'atom,xyz,cnfg ',nat,al
                  stop
              endif
! if j and j+R(supercell) are in the same group but with opposite ampterms
! then their sum is cancelled in the force as they both have the same
! displacement. This will lead to errors in evaluation of the FCs if all
! terms in the group cancel in this way, and the
! corresponding group in FC2 will be evaluated to be zero. That will also
! affect the ASR and produce a violation of the sum rules.
!              ired = igroup_2(t)  ! this is the corresponding index of ared
              rij = length(atompos(:,nat)-atompos(:,jatom))
!      write(ucor,5)iatom,taui,ni,al,jatom,tauj,nj,be,t,ired,rij,ampterm_2(t)
              write(ucor,6)iatom,taui,ni,al,jatom,tauj,nj,be,t,l,rij
 if(al.eq.1 .and. be.eq.1)  write(ulog,5)iatom,taui,nsi,jatom,tauj,nsj,t,l,rij
           endif
     enddo
     enddo
  enddo
  enddo

  write(ulog,*)' Correspondance of harmonic FCs and pairs of atoms established'

4 format(a,4(1x,i6))
5 format(2(i4,1x,' [ ',i4,' (',i2,',',i2,',',i2,') ] ',2x),':',i6,1x,i4,2x,f7.3,2x,f7.3)
6 format(2(i4,1x,' [ ',i4,' (',i2,',',i2,',',i2,') ] ',i1,2x),':',i6,1x,i4,2x,f7.3,2x,f7.3)

 end subroutine write_correspondance2
!============================================================
 subroutine warn(unt)
 implicit none
 integer unt

 write(unt,*)'********************************************************************'
 write(unt,*)'|                                                                  |'
 write(unt,*)'|                                                                  |'
 write(unt,*)'|      W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!     |'
 write(unt,*)'|      W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!     |'
 write(unt,*)'|      W    W  A    A  R    R  N N  N  II  N N  N  G       !!!     |'
 write(unt,*)'|      W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   !      |'
 write(unt,*)'|      WW  WW  A    A  R   R   N   NN  II  N   NN  G    G          |'
 write(unt,*)'|      W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!     |'
 write(unt,*)'|                                                                  |'
 write(unt,*)'|   There MIGHT be FC cancellation in the sum and perhaps errors   |'
 write(unt,*)'|                                                                  |'
 write(unt,*)'********************************************************************'
 end subroutine warn
!============================================================
 subroutine write_lat_fc(ngrps,ntrms)
 use svd_stuff
 use ios
 use force_constants_module
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
 real(8) rij
 integer i,j,ngrps(4),ntrms(4),j_sc

 write(ufco,*)' Crystal data: translation vectors of the primitive cell '
 write(ufco,9)r01
 write(ufco,9)r02
 write(ufco,9)r03
 write(ufco,*)' Crystal data: atoms in primitive cell: label,type,x,y,z,mass '
 write(ufco,*)natoms0
 do i=1,natoms0
   write(ufco,6)i,atom0(i)%name, atom0(i)%at_type,atom0(i)%equilibrium_pos,atom0(i)%mass
 enddo
 write(ufco,*)' Crystal data written ************************************'
 write(ufco,*)' Included ranks of FCs '
 write(ufco,*)include_fc(1),include_fc(2),include_fc(3),include_fc(4)
 write(ufco,*)' Number of FCs for each rank '
 write(ufco,*)ntrms (1),ntrms(2),ntrms(3),ntrms(4)
 write(ufco,*)' Number of independent FCs for each rank '
 write(ufco,*)ngrps(1),ngrps(2),ngrps(3),ngrps(4)
 write(ufco,*)' Neighborshell atoms: i,x,y,z,type_tau,n1,n2,n3 '
 write(ufco,*)natoms
 do i=1,natoms
   rij = length(atompos(:,i)-atompos(:,1))
   write(ufco,7)i,(atompos(j,i),j=1,3), iatomcell0(i),(iatomcell(j,i),j=1,3),rij
 enddo

 open(173,file='primlatt.xyz')
 write(173,*) natoms0
 write(173,*) ngrps(1),ngrps(2),ngrps(3),ngrps(4)
 do i=1,natoms0
   write(173,8)atom0(i)%name,atom0(i)%equilibrium_pos
!  write(173,7)atom0(i)%name,(atompos(j,i),j=1,3)
 enddo
 close(173)

 open(173,file='latfc.xyz')
 write(173,*) natoms
 write(173,*) ngrps(1),ngrps(2),ngrps(3),ngrps(4)
 do i=1,natoms
      call findatom_sc(iatomcell(:,i),iatomcell0(i),j_sc)  ! (=n(3),tau)
   write(173,8)atom0(iatomcell0(i))%name,(atompos(j,i),j=1,3),j_sc
 enddo
 close(173)

6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(a,3(2x,f19.10),i6)
9 format(9(2x,f19.10))
 end subroutine write_lat_fc
!============================================================

 !subroutine weighted_inversion(newamat,newbmat,Temp,energy_new,fc)
 !  integer i, j, l, dim_al_new, dim_ac_new
 !  real(8) V1, kB, Temp, kBT, norm_wts
 !  real(8), allocatable:: wts(:)
 !  real(8), allocatable:: fc(:), qmat(:)
 !  real(8), allocatable :: newamat(:,:),newbmat(:)
 !  real(8), allocatable:: energy_new(:),mat(:,:),mat_inverse(:,:) !an(:,:), bn(:)
 !  allocate(newamat(dim_al_new,dim_ac_new))!,newbmat(dim_ac_new))
 !  allocate(wts(dim_al_new),mat(dim_ac_new,dim_ac_new),mat_inverse(dim_ac_new,dim_ac_new),qmat(dim_ac_new))
   !allocate(an(dim_al_new,dim_ac_new),bn(dim_al_new),wts(dim_al_new),fc(dim_ac_new),energy_new(dim_al_new), &
   !qmat(dim_ac_new),mat(dim_ac_new,dim_ac_new),mat_inverse(dim_ac_new,dim_ac_new))
 !  dim_al_new=size(newamat,1)
 !  dim_ac_new=size(newbmat)
 !  write(*,*) "THE VALUE OF DIM_AL_NEW: ",dim_al_new
 !  write(*,*) "THE VALUE OF DIM_AC_NEW: ",dim_ac_new
 !  V1=1.0d0
 !  kB=0.00008617333262145d0
 !  kBT=kB*Temp
 !  norm_wts=0.0d0
 !  do i=1,dim_al_new
 !     wts(i)=exp(-1.0d0*(energy_new(i)-V1)/kBT)
 !     norm_wts=norm_wts+wts(i)
 !  enddo
 !  do i=1,dim_al_new
 !     wts(i)=wts(i)/norm_wts
 !  enddo
 !  qmat=0.0d0
 !  do i=1,dim_ac_new
 !     do l=1,dim_al_new
 !        qmat(i)=qmat(i)+wts(l)*newbmat(l)*newamat(l,i)
 !     enddo
 !     do j=1,dim_ac_new
 !        mat(i,j)=0.0d0
 !        do l=1,dim_al_new
 !           mat(i,j)=mat(i,j)+wts(l)*newamat(l,i)*newamat(l,j)
 !        enddo
 !     enddo
 !  enddo
 !  mat_inverse=0.0d0
 !  do i=1,dim_ac_new
 !     mat_inverse(i,i)=1.0d0
 !  enddo
 !  call inverse_real(mat,mat_inverse,dim_ac_new)
 !  fc=matmul(mat_inverse,qmat)
 !  write(*,*) "THE VALUE OF FC_WEIGHTED IS: ", fc
 !  deallocate(wts,qmat,mat,mat_inverse)
!end subroutine weighted_inversion
