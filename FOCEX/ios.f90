!===========================================================
 subroutine read_structure
! reads the param.inp file containing info on the atom types masses
! and coordinates within the primitive cell
 use ios
 use params
 use lattice
 use atoms_force_constants
 use svd_stuff
 implicit none
 integer i,counter,label
 real(8) scal
 character fdf*1,r234*3,invr*4,incl*4,it*1,zone*5,now*10,today*8
 real tim

 open(uparams,file='structure.params',status='old')

 read(uparams,*) latticeparameters   ! (a,b,c,alpha,beta,gamma)
 read(uparams,*) primitivelattice    ! prim latt vectors in terms of conventional above
 read(uparams,*) scal  ! scale factor for lattparams is read from params.inp
! read(uparams,*) maxneighbors       ! # of neigr-shells to use for each rank
 read(uparams,*) include_fc    ! if=0 do not include this rank
! read(uparams,*) nshells , radius      ! # of neigr-shells to use for each rank
! read(uparams,*) maxterms      ! max # of terms for each rank
! read(uparams,*) maxtermsindep ! max # of independent terms for each rank
 read(uparams,*) itrans,irot,ihuang,enforce_inv   ! translational and rotational invce flags (include if=1)
! read(uparams,*) tolerance, margin  ! tolerance for equating coords, margin for eliminating FCs
! read(uparams,*) svdcut    ! cutoff for smallest "eigenvalue" w to be included
 tolerance = 1d-4
! margin = 1d-5
 svdcut = 1d-9   ! default values
 read(uparams,*) itemp,tempk
 read(uparams,*) fdfiles , verbose ! number of force-displacement data files

 read(uparams,*) natom_type   ! # of different elements present in prim cell
! allocate(natom(natom_type))
 call allocate_mass(natom_type)
 read(uparams,*) (mas(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) (atname(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) natom_prim_cell        ! # of atoms in primitive cell
 read(uparams,*) fc2flag  ! if zero, take the default range from largest supercell
 nshells(1,1:natom_prim_cell)=0
 read(uparams,*) nshells(2,1:natom_prim_cell)
 read(uparams,*) nshells(3,1:natom_prim_cell)
 read(uparams,*) nshells(4,1:natom_prim_cell)

 maxneighbors=20 ! default for largest number of neighbor shells !maxval(nshells(2,:))
! if it's too small for rcut, nothing will happen and this value imposes the cutoff
! it it's large enough, we will go by rcut(2) in subroutine force_constants_init
! so we basically go by the lowest of the two. We choose a large maxneighbors for
! low-symmetry lattices
 if (itemp.eq.0) then
    if(fc2flag.eq.0) then
       it='0'
    else
       it='1'
    endif
 elseif(itemp.eq.1) then ! it='T'
    if(fc2flag.eq.0) then
       it='3'
    else
       it='4'
    endif
 endif
 write(fdf,'(i1)')fdfiles
 write(r234,'(3i1)')nshells(2,1),nshells(3,1),nshells(4,1)
 write(incl,'(4i1)')include_fc(:)
 invr='0000'
 if (itrans.ne.0) invr(1:1)='t'
 if (irot.ne.0) invr(2:2)='r'
 if (ihuang.ne.0) invr(3:3)='h'
 if (enforce_inv.ne.0) invr(4:4)='E'
 open(ulog  ,file='log'//fdf//it//'_'//r234//'_'//incl//invr//'.dat'   ,status='unknown')

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(tim)
 write(ulog,'(a,3x,a,3x,a,3x,a)')' Program FC234 was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(ulog,'(a,f10.4)')' STARTING TIME OF THE PROGRAM                   IS ',tim
 write(ulog,*)'===================================================================='

 write(ulog,*) svdcut,'    cutoff for smallest eigenvalue w to be included'
 write(ulog,*) tolerance,'   tolerance for equating two coordinates '
 write(ulog,*) include_fc,'  which ranks of FCs to include '
 write(ulog,*) fc2flag,'  if 0 default range of 15 ang is chosen for FC2 '
 write(ulog,*)' itemp,tempk = ',itemp,tempk
 write(ulog,*)' How many (not default) shells to include for ranks 2,3,4  '
 write(ulog,3)' 2 :', nshells(2,1:6)
 write(ulog,3)' 3 :', nshells(3,1:6)
 write(ulog,3)' 4 :', nshells(4,1:6)
 write(ulog,*) itrans,irot,ihuang,enforce_inv,'  transl, rot and Huang invce, enforcing inv'
 write(ulog,*) maxneighbors,'  default max neighbors  '
 write(ulog,*)' Reading ',fdfiles,' FORCEDISP & POSCAR files'
 write(ulog,*)' Included fcs and Imposed invariances: trans-rot-Huang=',invr
 write(ulog,*)'--------------------------------------------------'
 write(ulog,*)' Reading ',natom_type,' atom types with ',natom_prim_cell,'atoms in the primitive cell'

 latticeparameters(1:3) = latticeparameters(1:3)*scal
 call write_out(ulog,'Conventional Lattice parameters (1:3)',latticeparameters(1:3))
 call write_out(ulog,'Three angles in degrees              ',latticeparameters(4:6))
 call write_out(ulog,'Primitive lattice per conventional   ',primitivelattice)
 call allocate_primcell(natom_prim_cell)
! indices of atoms in primitive cell must be same order as in POSCAR
 counter = 0
 natom = 0
 write(ulog,*)'List of atoms: tau (or label), name,type,mass,reduced coordinates in conventional cell '
 do i=1,natom_prim_cell
! positions here must be D format in conventional cell for use by fcs_init
    read(uparams,*) label,atom_type(label),atompos0(:,label)
    if (i.ne.label) then
       print*,' positions must be sorted according to the labels ',i,label
       write(ulog,*)' READ_structure: positions must be sorted like labels ',i,label
       write(ulog,*)' Check atomic coordinates in primitive cell'
       stop
    endif
    atom0(i)%name = atname(atom_type(label))
    atom0(i)%at_type  = atom_type(label)
    atom0(i)%mass = mas(atom_type(label))
    atom0(i)%charge = 0 ! charge is=0 by default unless defined in dielectric.params
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
    write(ulog,4)i,atom0(i)%name,atom0(i)%at_type,atom0(i)%mass,atompos0(:,i)
 enddo


 close(uparams)

3 format(a,20(1x,i3))
4 format(i3,1x,a,i3,f8.3,3(1x,f10.5))

 end subroutine read_structure
!=========================================================
 subroutine read_dielectric
! reads ! born ! effective ! charge ! and ! dielectric ! constants ! from ! para.born
 use born
 use lattice
 use atoms_force_constants
 use ios
 integer i,j,k
 real(8) asr,aux(3,3)

 open(uborn,file='dielectric.params',status='old')
! allocate(zeu(3,3,natom_prim_cell))

 read(uborn,*) born_flag   ! if 0 use default
 write(ulog,*)' Born_flag  =',born_flag
 do i=1,3
    read(uborn,*)(epsil(i,j),j=1,3)
 end do

 aux=epsil
 call inverse_real(aux,epsinv,3)

 do k=1,natom_prim_cell
    do i=1,3
       read(uborn,*)(atom0(k)%charge(i,j),j=1,3)  !(zeu(i,j,k),j=1,3)
    end do
 end do
 close(uborn)

! impose ASR on Born charges: extra-charge is divided equally between atoms
 do i=1,3
 do j=1,3
      asr=0
      do k=1,natom_prim_cell
         asr=asr+atom0(k)%charge(i,j)  !zeu(i,j,k)
      end do
      asr=asr/natom_prim_cell
      atom0(:)%charge(i,j)=atom0(:)%charge(i,j) - asr
 end do
 end do

 write(ulog,*)' READ dielectric  constant='
 do i=1,3
    write(ulog,3)epsil(i,:)
 enddo
 write(ulog,*)' INVERSE dielectric  constant='
 do i=1,3
    write(ulog,3)epsinv(i,:)
 enddo

 do k=1,natom_prim_cell
 write(ulog,*)' ASR-enforced Born charges for atom ',k,atom0(k)%name
 do i=1,3
    write(ulog,3)atom0(k)%charge(i,:)
 enddo
 enddo

3 format(9(1x,g11.4))

 end subroutine read_dielectric
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
 integer i,j,n1,n2,n3,t
 real(8) a1,a2,a3,latt_const,om,a,b,c,dc(3),dr(3)
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
 write(ulog,*)' job title on line 1 of POSCAR is:'
 write(ulog,'(a)') line
! this is for the super-cell and has nothing to do with the one in params.inp
! which is for the primitive cell
 read(uposcar,*)lattice_parameter
 read(uposcar,*)rs1
 read(uposcar,*)rs2
 read(uposcar,*)rs3
 if (lattice_parameter .lt.0) then  ! it is the -volume of the cell
    volume_r = -lattice_parameter  ! need to find the volume(r1,r2,r3) for scaling
    call calculate_volume(rs1,rs2,rs3,om)
    latt_const = (volume_r/om)**(1./3.)
 else
    latt_const = lattice_parameter
 endif
 rs1=latt_const*rs1; rs2=latt_const*rs2; rs3=latt_const*rs3
 box(1)=length(rs1); box(2)=length(rs2); box(3)=length(rs3)
 call calculate_volume(rs1,rs2,rs3,volume_r)
!call calculate_volume(r01,r02,r03,volume_r0)
 call write_out(ulog,' Volume of supercell ',volume_r)
!call write_out(ulog,' Volume of primicell ',volume_r0)

 write(ulog,*)' box size in 1st direction=',box(1)
 write(ulog,*)' box size in 2nd direction=',box(2)
 write(ulog,*)' box size in 3rd direction=',box(3)

! number of atoms participating in MD,
 read(uposcar,*)(natom(i),i=1,natom_type)
 natom_super_cell = sum(natom)
! if (.not. (natom_super_cell*volume_r0/(natom_prim_cell*volume_r) .myeq. 1d0 ) ) then
!    write(ulog,*)' supercell inconsistency; check input coordinates again'
!    write(ulog,*)' natom_prim_cell, volume_r0=',natom_prim_cell,volume_r0
!    write(ulog,*)' natom_sc, volume_r=',natom_super_cell,volume_r
!    stop
! endif
 write(ulog,*)" number of different elements and number of atoms of each type"
 write(ulog,*) natom_type,(natom(i),i=1,natom_type)
 write(ulog,*)" total number of atoms in the supercell= ", natom_super_cell
 write(ulog,*) "translation vectors of the supercell are"
 write(ulog,4) rs1
 write(ulog,4) rs2
 write(ulog,4) rs3

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
       atom_sc(i)%equilibrium_pos = a*rs1+b*rs2+c*rs3
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


  if (.not. (natom_super_cell*volume_r0/(natom_prim_cell*volume_r) .myeq. 1d0 ) ) then
     write(ulog,*)' supercell inconsistency; check input coordinates again'
     write(ulog,*)' natom_prim_cell, volume_r0=',natom_prim_cell,volume_r0
     write(ulog,*)' natom_sc, volume_r=',natom_super_cell,volume_r
     stop
  endif

 write(ulog,*)'RECIPROCAL_LATTICE: '
 call make_reciprocal_lattice_2pi(rs1,rs2,rs3,gs1,gs2,gs3)
 call write_out(ulog,'volume_r ',volume_r)
 call write_out(ulog,'gsc1 ',gs1)
 call write_out(ulog,'gsc2 ',gs2)
 call write_out(ulog,'gsc3 ',gs3)

4 format(9(2x,f19.9))
6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))

 end subroutine read_crystal
!=================================================================
 subroutine identify_atoms_in_supercell
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 implicit none
 character line*90
 integer i,j,t,n1,n2,n3
 real(8) a(3),dr(3),a1,a2,a3

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
3 format(i5,9(2x,g13.6))
4 format(9(2x,f10.4))
5 format(a,3(i6,2x))
7 format(a,9(1x,f10.4))
8 format(a,i4,i3,i3,'(',3(i2,','),')',3(1x,f10.5))
 write(ulog,*)' POSCAR read successfully and closed'

 end subroutine identify_atoms_in_supercell
!===========================================================
 subroutine check_input_poscar_consistency_new
!! see if all atoms in the supercell can be obtained from the
!! atoms in the input file using translations of the primitive lattice
!! and assign atom_sc their type and other identities.
 use atoms_force_constants
 use ios
 use geometry
 use lattice
 use params
 implicit none
 integer i,j,k,ipr,isc, counter,ier,n1,n2,n3,tau1,nr1(3),nr2(3),nr3(3)
 real(8) a1,a2,a3
 type(vector) shift0,vec
 logical matched

! first describe r_i in terms of r0_i, and make sure the lin comb is integer
 write(ulog,*)' ----------------------------------------------------------'
 write(ulog,*)' CHECKING the commensurability between primcell and supercell'
 call check(rs1,a1,a2,a3,ier,g01,g02,g03)
 write(ulog,4)'checking r1 in units of r0s,ier=',a1,a2,a3,ier
 if (ier .eq. 1) stop
 nr1(1)=nint(a1); nr1(2)=nint(a2); nr1(3)=nint(a3)
 call check(rs2,a1,a2,a3,ier,g01,g02,g03)
 write(ulog,4)'checking r2 in units of r0s,ier=',a1,a2,a3,ier
 if (ier .eq. 1) stop
 nr2(1)=nint(a1); nr2(2)=nint(a2); nr2(3)=nint(a3)
 call check(rs3,a1,a2,a3,ier,g01,g02,g03)
 write(ulog,4)'checking r3 in units of r0s,ier=',a1,a2,a3,ier
 if (ier .eq. 1) stop
 nr3(1)=nint(a1); nr3(2)=nint(a2); nr3(3)=nint(a3)
 write(ulog,*)' COMMENSURABILITY CHECKED ---------------------------------'
 write(ulog,*)' nr1=',nr1
 write(ulog,*)' nr2=',nr2
 write(ulog,*)' nr3=',nr3

4 format(a,3(2x,f8.3),2x,i1)
7 format(a,2(1x,i4),9(2x,f9.4))
8 format(a,3(1x,i4),9(2x,f9.4))

! do i=1,natom_prim_cell
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

 checkloop: do i=1,natom_prim_cell  !,natom_super_cell
    shift0 = atom0(i)%equilibrium_pos - atom_sc(1)%equilibrium_pos
    write(ulog,7)'atom in PRIMCELL:',i,atom0(i)%at_type,atom0(i)%equilibrium_pos
    write(ulog,3)'trying shift vector=',shift0

! try this shift see if all atoms in SC can be mapped to the PC by it.
    SC: do k = 1,natom_super_cell
       matched = .false.
       PRIM: do j = 1,natom_prim_cell  ! one of the prim-cell atoms has to match
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

 if (i.ge.1 .and. i.le.natom_prim_cell) then
    call write_out(ulog,'THE shift vector ', shift0)
    write(ulog,*)' it maps atom 1 of the SC to atom ',i,' of the prim cell'
 else
    write(ulog,*)' NO SHIFT VECTOR COULD BE DEFINED!, check your coordinates again;i=',i
    stop
 endif

! shift all the atoms in the primitive cell by -shift0 so that they fall
! on the atoms in the supercell
! write(ulog,*)'The shifted (cartesian) positions in the prim cell are now:'
! do i=1,natom_prim_cell
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
         PRIM2: do j = 1,natom_prim_cell  ! one of the prim-cell atoms has to match
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
                atom_sc(k)%charge = atom0(j)%charge

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
 subroutine count_configs(outcar,ncfg)
 use ios , only : ulog,utraj
 use atoms_force_constants, only : natom_super_cell
 implicit none
 integer, intent(out) :: ncfg
 character, intent(in):: outcar*(*)
 integer t,j,i
 character line*99
 logical found,exst

 write(*,*)' opening OUTCAR file'

 inquire(file=outcar,exist=exst)
 if(exst) then
      open(utraj,file=outcar,status='old')
 else
      write(ulog,*)' outcar file ',outcar,' does not exist; check your files location and run again'
      stop
 endif

! first find the number of configurations: ncfg
 t=0
 do j=1,11000000
    read(utraj,'(a)',end=99)line
    call findword('POSITION',line,found)
    if (found) then
       t = t+1
    endif
 enddo
99 write(ulog,*)' reached the end of OUTCAR file; number of configurations= ',t
 close(utraj)
 ncfg = t
 write(*,*)' reached the end of OUTCAR file which had ',ncfg,' configurations '
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

 end subroutine count_configs
!===========================================================
 subroutine read_force_position_data(outcar,ncfg,energy,dsp,frc)
!! reads contents of OUTCARi; second lines are energies, frc_constr=3*natom_super_cell*ncfg
!! outputs are dsp(3,NSC,ncfg),frc(3,NSC,ncfg),energy(ncfg) for the OUTCARi file read
 use ios
 use atoms_force_constants
 use geometry
 use params
 use lattice
 use svd_stuff
 implicit none
 integer, intent(in) :: ncfg
 character, intent(in):: outcar*(*)
 integer n1,n2,i,t,j,k
 type(vector) v
 real(8) x1,x2,x3,x4,const,rr(3)
 real(8), save :: emin
 logical, save :: first_call=.true.
! real(8), allocatable, intent(inout) :: energy(:),dsp(:,:,:),frc(:,:,:)
 real(8), intent(out) :: energy(ncfg),dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg)

 character line*99
 logical found,exst

 write(*,*)' REopening OUTCAR file'

 inquire(file=outcar,exist=exst)
 if(exst) then
    open(utraj,file=outcar,status='old')
 else
    write(ulog,*)' outcar file ',outcar,' does not exist; check your files location and run again'
    stop
 endif


! now get the FORCES from OUTCAR file
 t=0 ; energy=0
 do j=1,ncfg
    read(utraj,'(a)',end=88)line
!    if (line(1:8) .eq. "POSITION" ) then
    call findword('POSITION',line,found)
    if (found) then
       t = t+1
       read(utraj,*) k,energy(t)  ! start with 1 since t=0
       if(t.ne.j ) then !.or. t-1.ne.k) then
          write(*,*)'Error in reading snapshot#s in OUTCAR?',j,t
       endif
       do i=1,natom_super_cell
           read(utraj,*) dsp(1:3,i,t),frc(1:3,i,t)
       enddo
       write(*,*)'j,t,k=',j,t,k
    endif
 enddo
88 write(ulog,*)' reached the end of OUTCAR file after ',t,' steps'
 close(utraj)
 write(ulog,*)' last line read was '
 write(ulog,'(a)') line
 close(utraj)
 if (t .ne. ncfg) then
    write(*,*)'ERROR in reading the force file OUTCAR'
    write(*,*)'ncfg, # of read steps=',ncfg,t
    stop
 endif

 write(*,*)'writing frc and dsp in the log file'
 call write_out(ulog,' last force ',frc(:,:,t))
 call write_out(ulog,' last coord ',dsp(:,:,t))


! get energy per primitive unit cell so that it does not depend on supercell size
 energy=energy/natom_super_cell*natom_prim_cell

! subtract lowest energy value ! get it from OUTCAR1; it is arbitrary anyways
 if( first_call ) then
    emin=minval(energy)
    first_call=.false.
 endif
 energy=energy-emin

 write(*,9)'Energy/primcell assigned ',energy
 write(*,*)'Calling calculate_and_write_displacements'
 call calculate_and_write_displacements(ncfg,dsp,frc)
! inputs: displ, including atom_sc%equilibrium_pos 
! outputs: displ - atom_sc%equilibrium_pos  

 write(*,*)'calling write_correspondance'
 call write_correspondance

 write(*,*)'exiting read_force_position_data '
9 format(a,200(1x,f7.3))

 end subroutine read_force_position_data
!===========================================================
 subroutine read_force_position_data2(outcar,frc_constr,ncfg,energy,dsp,frc)
!! reads contents of OUTCARi; second lines are energies, frc_constr=3*natom_super_cell*ncfg
!! outputs are nfcg,dsp(3,NSC,ncfg),frc(3,NSC,ncfg),energy(ncfg) for every outcar file read
 use ios
 use atoms_force_constants
 use geometry
 use params
 use lattice
 use svd_stuff
 implicit none
 integer, intent(out) :: frc_constr,ncfg
 character, intent(in):: outcar*(*)
 integer n1,n2,i,t,j,k
 type(vector) v
 real(8) x1,x2,x3,x4,const,rr(3),emin
! real(8), intent(inout) :: energy(ncfg),dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg)
 real(8), allocatable, intent(inout) :: energy(:),dsp(:,:,:),frc(:,:,:)
! integer, allocatable :: nlines(:)

 character line*99
 logical found,exst

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
    if (found) then
       t = t+1
    endif
 enddo
99 write(ulog,*)' reached the end of OUTCAR file; number of configurations= ',t
 ncfg = t
 if(t.eq.0) then
    write(ulog,*)' the word POSITION was not found in OUTCAR file, check it!'
    stop
 endif

! allocate( energy(ncfg),dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg) )
! if (ncfg.gt.5000) then
!    write(*,*)'READ_FORCE_POSITION_DATA: reading more than 5000 configurations ',ncfg
!    write(*,*)'either put in less than 5000 configurations or increase the size of energy in this subroutine'
!    stop
! endif

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
 if ( allocated(dsp)) deallocate(dsp)
 if ( allocated(frc)) deallocate(frc)
 if ( allocated(energy)) deallocate(energy)
! if ( allocated(nlines)) deallocate(nlines)

! call allocate_pos(natom_super_cell,ncfg) ! for displ and force
 allocate( energy(ncfg),dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg) )
! allocate( nlines(ncfg) )

! now get the FORCES from OUTCAR file
 rewind(utraj)
 t=0 ; energy=0
 do j=1,11000000
    read(utraj,'(a)',end=88)line
!    if (line(1:8) .eq. "POSITION" ) then
    call findword('POSITION',line,found)
    if (found) then
       t = t+1
       read(utraj,*) k,energy(t)  ! start with 1 since t=0
       if(t.ne.j .or. t.ne.k) then
          write(ulog,*)'Error in reading snapshot#s in OUTCAR?',k,j,t
       endif
       do i=1,natom_super_cell
           read(utraj,*) dsp(1:3,i,t),frc(1:3,i,t)
       enddo
!       nlines(t)=natom_super_cell
    endif
 enddo
88 write(ulog,*)' reached the end of OUTCAR file after steps= ',t
 write(ulog,*)' last line read was '
 write(ulog,'(a)') line
 close(utraj)
 if (t .ne. ncfg) then
    write(ulog,*)'ERROR in reading the force file OUTCAR'
    write(ulog,*)'ncfg, # of read steps=',ncfg,t
    stop
 endif

 call write_out(ulog,' last force ',frc(:,natom_super_cell,t))
 call write_out(ulog,' last coord ',dsp(:,natom_super_cell,t))

 frc_constr = ncfg *natom_super_cell*3

! get energy per unit cell so that it does not depend on supercell size
 energy=energy/natom_super_cell*natom_prim_cell

! subtract lowest energy value
 emin=minval(energy)
 energy=energy-emin

 call calculate_and_write_displacements(ncfg,dsp,frc)
! inputs: displ, including atom_sc%equilibrium_pos 
! outputs: displ - atom_sc%equilibrium_pos  

 call write_correspondance

 write(*,*)'exiting read_force_position_data2'
9 format(a,200(1x,f7.3))

 end subroutine read_force_position_data2
!===========================================================
 subroutine write_independent_fcs(n,sig,sd,ulog)
 use svd_stuff
 use atoms_force_constants
 use params
 implicit none
 integer, intent(in) :: ulog,n
 real(8), intent(in) :: sig(n)
 real(8), intent(out):: sd(4)
 integer i,k,g,cnt2,cnt


!   if(nindepfc.eq.dim_ac) then
!      call check_huang(fcs,nindepfc)
!   else
!      write(ulog,*)' **************ERROR in WRITE_INDEPENDENT_FCS ****************'
!      write(ulog,*)' cannot check hunag as nindepfc,dim_ac=',nindepfc,dim_ac
!   endif

   write(ulog,*)' group, fcs(group), sigma(group) for all SVD'
   do i=1,n
      write(ulog,4) i,fcs(i),sig(i),sig(i)/sqrt(1.*dim_al)
   enddo
   write(ulog,*) "rank, average error"

!! k1 -----------------
 map(2)%ngr = size_kept_fc2
!! k1 -----------------

   write(ulog,*)' map(2)%ngr=',map(2)%ngr
   write(ulog,*)' dim_ac    =',dim_ac,n

   cnt2=0;sd=0
   do i=1,4
      if(include_fc(i).ne.1) cycle
      cnt=0
      do g=1,map(i)%ngr
        do k=1,map(i)%ntind(g)
            cnt=cnt+1
            cnt2=cnt2+1 !cnt
            write(ulog,*)'rank,group,nind,cnt=',i,g,k,cnt
            sd(i)=sd(i)+sig(cnt2)/sqrt(1.*dim_al)
        enddo
      enddo
      if(cnt .ne. 0 ) sd(i)=sd(i)/cnt
      write(ulog,3)'Rank, sd(rank)=', i,sd(i)
   enddo
   write(ulog,5)' !==== Error summary for each rank',(sd(i),i=1,4)


3 format(a,i6,3(2x,g14.7))
4 format(i6,3(2x,g14.7))
5 format(a,9(2x,g11.4))

 end subroutine write_independent_fcs
!===========================================================
 subroutine write_neighbors
 use atoms_force_constants
! use force_constants_module
 use params
 use ios
 implicit none
 integer i0,shel_count,j,nm(3),ta,nbmx
 real(8) dij

 do i0=1,natom_prim_cell
    write(ulog,*)' ******************************'
    write(ulog,*)' Neighbors of atom number ',i0
    write(ulog,*)' ******************************'
    do shel_count=1,maxshell
       nbmx = atom0(i0)%shells(shel_count)%no_of_neighbors
       dij  = atom0(i0)%shells(shel_count)%rij
       do j=1,min(nbmx,500)
          ta =  atom0(i0)%shells(shel_count)%neighbors(j)%tau
          nm =  atom0(i0)%shells(shel_count)%neighbors(j)%n
          write(ulog,3)' shell,dij,nb#,tau,n=',shel_count,dij,j,ta,nm
       enddo
       write(ulog,2)'WRITE_NEIGHBORS: shel#,rij,nb#',shel_count,dij,j-1
    enddo
 enddo

 write(ulog,*)' ************ End of the neighborlist ************** '
2 format(a,i4,2x,f8.4,2x,i4)
3 format(a,2x,i3,2x,f8.4,2x,i3,4x,i3,' (',3(1x,i4),')')
 end subroutine write_neighbors
!============================================================
 subroutine write_output_fcs
 use svd_stuff
 use ios
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
 integer rnk,t,ti,i,res,j,rs
 integer iat(4),ixyz(4),g,ng,term,term2,cnt2,frm
 real(8) rij,bunit,one,fcd,trace,dij
 character frmt*2,goh*48,ln*1,geh*47

 bunit = ryd/ab/ab
 one =1d0

6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(2(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
9 format(9(2x,f19.10))

! first write the crystal data
 call write_lat_fc(map(:)%ntotind,map(:)%ntot)

!----------------------------------------
 res = 0
 ranks: do rnk=1,4
  if ( include_fc(rnk) .eq. 1 ) then
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
    groups: do g=1,map(rnk)%ntotind  ! index of a given group

!! K1 new change made on 2/13/23 ----------------
       if(rnk.eq.2) then
          if(keep_grp2(g).ne.1) cycle
       endif
!! K1 new change made on 2/13/23 ----------------

       if(g.gt.1) cnt2=cnt2+map(rnk)%ntind(g-1)

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
!          if( abs(fcd) .gt. margin) then
             write(ufc1-1+rnk,geh)t,g, (iat(j),ixyz(j),j=1,rnk),fcd,one
!          endif
       enddo
  enddo groups
  res = res+map(rnk)%ngr
  endif
 enddo ranks

write(ulog,*)'******* Trace for the harmonic FCs ********'
 open(456,file='trace_fc.dat')

! write the trace of FC2
 rnk=2
 iloop: do i=1,natom_prim_cell
 jloop: do j=1,natoms
     rij = length(atompos(:,i)-atompos(:,j))
     trace=0
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
              if (iat(1).ne.i .or. iat(2).ne.j) cycle !iloop
              write(*,*)'i,j,term2,al,be=',i,j,term2,ixyz(1),ixyz(2)
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
     endif
  enddo jloop
  enddo iloop

  close(456)

  write(ulog,*)'***************** END OF FC Trace ******************'

  if (res.ne.nindepfc) then
     write(ulog,*)'WRITE_OUTPUT_FCS: sum(nterms),ngr=',res,nindepfc
     write(ulog,*)'WRITE_OUTPUT_FCS: they should be equal!'
  endif

 end subroutine write_output_fcs
!============================================================
 subroutine read_fcs(iunit,fn,rank,fc,gr)
!! reads force constants from a file if the file exists
 use svd_stuff
 use ios
 use atoms_force_constants
 use params
 implicit none
 integer, intent(in) :: rank,iunit,gr
 integer t,ti,igr,res,i,a, b, c, d, e, f,reason,cnt2,term,g,j
 real amp, rijs,fc(gr)
 character line*99,fn*11
 logical ex

     inquire(file=fn,exist=ex)
     if (ex) then
        open(iunit ,file=fn ,status='old')
     else
        write(*   ,*)'For rank=',rank,' there is no FC data file present '
        write(*   ,*)'check your inputs and run the program again'
        write(ulog,*)'For rank=',rank,' there is no FC data file present '
        write(ulog,*)'check your inputs and run the program again'
        stop
     endif

     term = 0;   
     groups: do g=1,map(rank)%ntotind  ! index of a given group
       map(rank)%gr(g)%iatind(:,:)=0
       map(rank)%gr(g)%ixyzind(:,:)=0
       do ti=1,map(rank)%ntind(g)  ! index of independent terms in that group g
          term = term+1
          read(iunit,*,err=91) t,igr,(map(rank)%gr(igr)%iatind (j,t),  &
           &                          map(rank)%gr(igr)%ixyzind(j,t), j=1,rank), &
           &                          fc(term) 
       enddo
     enddo groups

91   write(ulog,*)'READ_FCS: rank=',rank,' reached end of file after ',term-1,' lines'
     if (term.eq.1) then
        write(ulog,*)'READ_FCS: the file exists but is empty! rank=',rank
        write(*   ,*)'READ_FCS: the file exists but is empty! rank=',rank
        stop
!       include_fc(rank) = 0       ! if nothing, then exclude from fitting
     endif
     close(iunit)

 end subroutine read_fcs
!============================================================
 subroutine read_fcs_2(iunit,rank,fc2,gr)
 use svd_stuff
 use ios
! use force_constants_module
 use atoms_force_constants
 use params
 implicit none
 integer rank,iunit,t,res,i,a, b, c, d, e, f,reason,cnt2,term,g,ti,j,gr
 real amp, rijs,fc2(gr)
 character line*99
 logical ex

! read(iunit,'(a)') line
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

elseif ( rank .eq. 2) then

 if ( include_fc(rank) .eq. 2 ) then
!   write(*,*) "The value of map(2)%ngr and rank is: ", map(2)%ngr, rank
   inquire ( file='fc2_irr.dat', exist=ex)
   if (ex) then
!     open(473,file='fc2_irr.dat',status='old')

     cnt2=0;   term = 0;   
     groups: do g=1,map(rank)%ngr  ! index of a given group
       map(2)%gr(g)%iatind(:,:)=0
       map(2)%gr(g)%ixyzind(:,:)=0

!! K1 new change made on 2/13/23 ----------------
!       if(rank.eq.2) then
!          if(keep_grp2(g).ne.1) cycle
!       endif
!! K1 new change made on 2/13/23 ----------------

       if(g.gt.1) cnt2=cnt2+map(2)%ntind(g-1)

       do t=1,map(2)%ntind(g)  ! index of independent terms in that group g
          term = term+1
           read(ufit2,*) ti,gr,map(2)%gr(gr)%iatind (1,ti),  &
           &                   map(2)%gr(gr)%ixyzind(1,ti),  &
           &                   map(2)%gr(gr)%iatind (2,ti),  &
           &                   map(2)%gr(gr)%ixyzind(2,ti),  & 
           &                   fc2(res+cnt2+ti),a,b 
       enddo

!        do i=1, map(2)%ngr
!           read(473,*,iostat=reason) a, b, c, d, e, f, fc_ind(i), amp, rijs
!           if (reason > 0) then
!              write(*,*) "something is wrong"
!           elseif (reason < 0) then
!              write(*,*) "END OF FILE REACHED"
!           else
!              write(*,*) a,b,c,d,e,f,fc_ind(i),amp,rijs
!           endif
!        enddo

        enddo groups
     close(ufit2)
  endif
 endif
 return
! this is how fc_irr were written
!    cnt2=0;   term = 0;    term2= 0
!    groups: do g=1,map(rnk)%ntotind  ! index of a given group

!! K1 new change made on 2/13/23 ----------------
!       if(rnk.eq.2) then
!          if(keep_grp2(g).ne.1) cycle
!       endif
!! K1 new change made on 2/13/23 ----------------

!       if(g.gt.1) cnt2=cnt2+map(rnk)%ntind(g-1)

    ! write in the log and fcn_fit.dat: cnt2+ti is the position of indep_fc of that rank
!       do ti=1,map(rnk)%ntind(g)  ! index of independent terms in that group g
!          iat(1:rnk)  = map(rnk)%gr(g)%iatind (:,ti)
!          ixyz(1:rnk) = map(rnk)%gr(g)%ixyzind(:,ti)
!          term = term+1
!          if (rnk.eq.2) then
!            rij = length(atompos(:,iat(1))-atompos(:,iat(2)))
!          else
!            rij = 0
!          endif
!          write(ufit1-1+rnk,geh) ti,g,(iat(j),ixyz(j),j=1,rnk),  &  ! for fcn_irr.dat
!       &     fcs(res+cnt2+ti),one,rij
!       enddo

! res =  igroup_1(nterms(1))
! do i=1,nterms(rank)
!       read(iunit,*,err=92)t,igroup_2(t), &
! & iatomterm_2(1,t),ixyzterm_2(1,t),iatomterm_2(2,t),ixyzterm_2(2,t),  &
! & fcs_2(igroup_2(t)),ampterm_2(t)
! enddo
! return
! 92 write(ulog,*)'READ_FCS: rank=',rank,' reached end of file after t=',t
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

stop

 end subroutine read_fcs_2
!=================================================================
 subroutine calculate_and_write_displacements(ncfg,dsp,frc)
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 implicit none
 integer, intent(in) :: ncfg
 real(8), intent(inout) :: dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg)
 character line*90
 integer i,j,t,step
 real(8) dc(3),dr(3)

 write(ulog,*)' t , particle#, cartesian disp & forces u1,u2,u3,f1,f2,f3'
 write(ulog,*)' Number of onfigurations NCFG=',ncfg
! get displacements from positions
 do t=1,ncfg
    do i=1,natom_super_cell
! first get direct coordinates, then take the distance using pbc, then
! transfrom back to cartesian coordinates
       dc(1) = dsp(1,i,t) - atom_sc(i)%equilibrium_pos%x
       dc(2) = dsp(2,i,t) - atom_sc(i)%equilibrium_pos%y
       dc(3) = dsp(3,i,t) - atom_sc(i)%equilibrium_pos%z
       call cart_to_direct_aa(dc,dr)
       dr(1) = dr(1) - anint(dr(1))
       dr(2) = dr(2) - anint(dr(2))
       dr(3) = dr(3) - anint(dr(3))
       call direct_to_cart_aa(dr,dsp(:,i,t))
    enddo
 enddo
 step=ncfg-1
 if(verbose) step=1
 do t=1,ncfg ,step           ! write first and last displ-force data
    do i=1,natom_super_cell
       write(ulog,6)t,i,(dsp(j,i,t),j=1,3),(frc(j,i,t),j=1,3)
    enddo
 enddo

 6 format(2(i5),3(1x,f10.6),3x,3(1x,g12.5))

 end subroutine calculate_and_write_displacements
!=================================================================
      subroutine ustring(m,lineout,nrank,iatom,ixyz)
      implicit none
      integer i,j,k,m,n,nrank,iatom(nrank),ixyz(nrank)
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
! checking consistency between POSCAR:atom_sc%equilibrium_pos and FORCEDISP:displ(:,:,1)

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
       write(ulog,*)' atom # ',i,' in POSCAR and FORCEDISP are different'
       write(ulog,*)' POSCAR:atom_SC=',atom_sc(i)%equilibrium_pos
       write(ulog,*)' FORCEDISP: displ =',displ(:,i,1)
       write(ulog,*)' check your input files '
       write(*,*)' atom # ',i,' in POSCAR and FORCEDISP are different'
       write(*,*)' check your input files '
       stop
    endif

 enddo

 end subroutine pos_out_consistency
!============================================================
 subroutine write_correspondance
!! write correspondance between supercell atoms and primitive structure
 use params
 use atoms_force_constants
 use svd_stuff
 use ios
 use geometry
! use force_constants_module
 use lattice
 implicit none
 integer t,ired,nat,iatom,jatom,g,rnk
 integer al,be,taui,tauj,ni(3),nj(3),nsi(3),nsj(3)
 real(8) rij

 rnk = 2
 write(ulog,*)' HARMONIC FCs: correspondance between a given FC and all pairs '
 write(ulog,*)' Basically we describe the harmonic terms in the Taylor expansion'
 write(ulog,*)' Which atom pairs in the supercell are bridged by a FC element '
 write(ulog,*)' Below, ni refers to the translations of the supercell'
 write(ulog,*)'# iatom, [  taui(nsi) ],  alpha_i=1 ;  j, [  tauj(nsj) ], alpha_j=1  :   term group  rij'
 write(ucor,*)'# iatom, [  taui(ni)  ],  alpha_i ;  jatom, [  tauj(nj), alpha_j  :   term group  rij'

! do nat=1,natom_super_cell
  do nat=1,natom_prim_cell
  do al=1,3             ! this is how we order each line

!    taui = atom_sc(nat)%cell%tau
!    ni   = atom_sc(nat)%cell%n
     taui = iatomcell0(nat)   ! taui should be = nat
     ni(:)= iatomcell(:,nat)
     call findatom_sc(ni,nat,iatom)  ! now iatom is the supercell index
     nsi  = nint(matmul(r0g,dfloat(ni)))

     do g=1,map(rnk)%ngr  ! sum over groups
!           if(l.gt.1) cnt2 = cnt2 + map(rnk)%ntind(l-1)
     do t=1,map(rnk)%nt(g) ! sum over all terms in that group
           if ( taui.eq. map(rnk)%gr(g)%iat(1,t) .and.  &
           &    al  .eq. map(rnk)%gr(g)%ixyz(1,t) ) then
              tauj = iatomcell0(map(rnk)%gr(g)%iat(2,t))
              nj(:)= iatomcell(:,map(rnk)%gr(g)%iat(2,t))+ni(:) ! trnslte by ni
              nsj  = nint(matmul(r0g,dfloat(nj)))
              be   = map(rnk)%gr(g)%ixyz(2,t)
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
              write(ucor,6)iatom,taui,ni,al,jatom,tauj,nj,be,t,g,rij
 if(al.eq.1 .and. be.eq.1)  write(ulog,5)iatom,taui,nsi,jatom,tauj,nsj,t,g,rij
           endif
     enddo
     enddo
  enddo
  enddo

  write(ulog,*)' Correspondance of harmonic FCs and pairs of atoms established'

4 format(a,4(1x,i6))
5 format(2(i4,1x,' [ ',i2,' (',i2,',',i2,',',i2,') ] ; ',1x),i6,1x,i4,2x,f7.3,2x,f7.3)
6 format(2(i4,1x,' [ ',i2,' (',i2,',',i2,',',i2,') ] ',i1,1x,' ; '),i6,1x,i4,2x,f7.3,2x,f7.3)

 end subroutine write_correspondance
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
! use force_constants_module
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
 real(8) rij
 integer i,j,ngrps(4),ntrms(4)

 write(ufco,*)' Crystal data: translation vectors of the primitive cell '
 write(ufco,9)r01
 write(ufco,9)r02
 write(ufco,9)r03
 write(ufco,*)' Crystal data: atoms in primitive cell: label,type,x,y,z,mass '
 write(ufco,*)natom_prim_cell
 do i=1,natom_prim_cell
   write(ufco,6)i,atom0(i)%name, atom0(i)%at_type,atom0(i)%equilibrium_pos,atom0(i)%mass
 enddo
 write(ufco,*)' Crystal data written ************************************'
 write(ufco,*)' Included ranks of FCs '
 write(ufco,*)  include_fc(:)
 write(ufco,*)' Number of FCs for each rank '
 write(ufco,*)  ntrms(:)
 write(ufco,*)' Number of independent FCs for each rank '
 write(ufco,*)  ngrps(:)
 write(ufco,*)' Neighborshell atoms: i,x,y,z,type_tau,n1,n2,n3 '
 write(ufco,*)  natoms
 do i=1,natoms
   rij = length(atompos(:,i)) !-atompos(:,1))
   write(ufco,7)i,(atompos(j,i),j=1,3), iatomcell0(i),(iatomcell(j,i),j=1,3),rij
 enddo

 open(173,file='primlatt.xyz')
 write(173,*) natom_prim_cell
 write(173,*) ngrps(:)
 do i=1,natom_prim_cell
   write(173,8)atom0(i)%name,atom0(i)%equilibrium_pos
 enddo
 close(173)

 open(173,file='latfc.xyz')
 write(173,*) natoms
 write(173,*) ngrps(:)
 do i=1,natoms
   write(173,8)atom0(iatomcell0(i))%name,(atompos(j,i),j=1,3)
 enddo
 close(173)

6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(a,6(2x,f19.10))
9 format(9(2x,f19.10))
 end subroutine write_lat_fc
!============================================================
 subroutine write_invariance_violations(ulog,n,fcs)
 use svd_stuff, only : atransl,arot,brot,ahuang, amat,bmat,transl_constraints, &
 &  rot_constraints,huang_constraints,inv_constraints,dim_al,itrans,irot,ihuang
 use ios, only : umatrx
 implicit none
 integer, intent(in) :: ulog,n
 real(8), intent(in) :: fcs(n)
 integer i
 real(8) prod,errmax,junk,err,num,denom


   write(ulog,3)'Dimension of extracted FCS is=',n
   write(*,3)'Dimension of extracted FCS is=',n
   write(ulog,*)' dimensions of ahuang=',size(ahuang(:,1)),size(ahuang(1,:))
   write(*,*)' dimensions of ahuang=',size(ahuang(:,1)),size(ahuang(1,:))
   write(*,*)' # of huang constraints is=',huang_constraints
   write(ulog,*)' # of huang constraints is=',huang_constraints

   write(ulog,3)'******* Violation of translational invariance relations: enforced? ',itrans
   do i=1,transl_constraints
      write(ulog,*)i,dot_product(atransl(i,:),fcs(:))
   enddo
   write(ulog,3)'******* Violation of rotational invariance relations: enforced? ',irot
   do i=1,rot_constraints
      write(ulog,*)i,dot_product(arot(i,:),fcs(:)),brot(i)
   enddo
   write(ulog,3)'******* Violation of Huang invariance relations: enforced? ',ihuang
   do i=1,huang_constraints
      write(ulog,*)i,dot_product(ahuang(i,:),fcs(:))
   enddo

   write(umatrx,*)'******* Violation of force-displacement relations: '
   errmax=-1
   err=0;num=0;denom=0
   do i=inv_constraints+1,dim_al
      prod = dot_product(amat(i,:),fcs(:))
      junk = abs(prod-bmat(i))
      write(umatrx,3)' ',i,prod,bmat(i),junk
      if (errmax .lt. junk) errmax=junk
      err=err+junk
      num=num+junk*junk
      denom=denom+prod*prod
   enddo
   err=err/(dim_al-inv_constraints)
   write(ulog,*)'Max and average errors in force-displacements,percent deviation=', &
   &    errmax,err,sqrt(num/denom)

3 format(a,i5,3(2x,g12.5))
 end subroutine write_invariance_violations
