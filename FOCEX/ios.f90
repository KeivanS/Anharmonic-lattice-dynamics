!===========================================================
 subroutine read_structure
! reads the param.inp file containing info on the atom types masses
! and coordinates within the primitive cell
 use ios
 use params
 use lattice
 use atoms_force_constants
 use svd_stuff
 use constants, only : r15
 implicit none
 integer i,counter,tau , bornflag
 real(r15) scal
 character fdf*1,r4*7,invr*4,incl*8,it*2,zone*5,now*10,today*8,born*2
 real tim

 open(uparams,file='default.params',status='old')
   read(uparams,*) tolerance    ! cutoff for smallest "eigenvalue" w to be included
   if(tolerance.eq.0) tolerance = 2d-3
   read(uparams,*) svdcut    ! cutoff for smallest "eigenvalue" w to be included
   if(svdcut.eq.0) svdcut = 1d-10   ! default values
   read(uparams,*) icutoff    ! cutoff for neighborshells
   if(icutoff.eq.0) icutoff = 7     ! default multiple of bond length; defines maxatoms
   read(uparams,*) maxterms(1)
   if(maxterms(1).eq.0) maxterms(1)=40 !100
   read(uparams,*) maxterms(2)
   if(maxterms(2).eq.0) maxterms(2)=4000 !500
   read(uparams,*) maxterms(3)
   if(maxterms(3).eq.0) maxterms(3)=5000 !1800
   read(uparams,*) maxterms(4)
   if(maxterms(4).eq.0) maxterms(4)=3000 !2000
   read(uparams,*) maxterms(5)
   if(maxterms(5).eq.0) maxterms(5)=3  
   read(uparams,*) maxterms(6)
   if(maxterms(6).eq.0) maxterms(6)=3 
   read(uparams,*) maxterms(7)
   if(maxterms(7).eq.0) maxterms(7)=3 
   read(uparams,*) maxterms(8)
   if(maxterms(8).eq.0) maxterms(8)=3 
   read(uparams,*) maxtermzero(1)
   if(maxtermzero(1).eq.0) maxtermzero(1)=20 !500
   read(uparams,*) maxtermzero(2)
   if(maxtermzero(2).eq.0) maxtermzero(2)=2000 !2000
   read(uparams,*) maxtermzero(3)
   if(maxtermzero(3).eq.0) maxtermzero(3)=5000!5000
   read(uparams,*) maxtermzero(4)
   if(maxtermzero(4).eq.0) maxtermzero(4)=3000 !8000
   read(uparams,*) maxtermzero(5)
   if(maxtermzero(5).eq.0) maxtermzero(5)=1
   read(uparams,*) maxtermzero(6)
   if(maxtermzero(6).eq.0) maxtermzero(6)=1
   read(uparams,*) maxtermzero(7)
   if(maxtermzero(7).eq.0) maxtermzero(7)=1
   read(uparams,*) maxtermzero(8)
   if(maxtermzero(8).eq.0) maxtermzero(8)=1
   read(uparams,*) maxtermsindep(1)
   if(maxtermsindep(1).eq.0) maxtermsindep(1)=10
   read(uparams,*) maxtermsindep(2)
   if(maxtermsindep(2).eq.0) maxtermsindep(2)=250
   read(uparams,*) maxtermsindep(3)
   if(maxtermsindep(3).eq.0) maxtermsindep(3)=150
   read(uparams,*) maxtermsindep(4)
   if(maxtermsindep(4).eq.0) maxtermsindep(4)=300
   read(uparams,*) maxtermsindep(5)
   if(maxtermsindep(5).eq.0) maxtermsindep(5)=1
   read(uparams,*) maxtermsindep(6)
   if(maxtermsindep(6).eq.0) maxtermsindep(6)=1
   read(uparams,*) maxtermsindep(7)
   if(maxtermsindep(7).eq.0) maxtermsindep(7)=1
   read(uparams,*) maxtermsindep(8)
   if(maxtermsindep(8).eq.0) maxtermsindep(8)=1
   read(uparams,*) maxgroups(1)
   if(maxgroups(1).eq.0) maxgroups(1)=10
   read(uparams,*) maxgroups(2)
   if(maxgroups(2).eq.0) maxgroups(2)=200
   read(uparams,*) maxgroups(3)
   if(maxgroups(3).eq.0) maxgroups(3)=150
   read(uparams,*) maxgroups(4)
   if(maxgroups(4).eq.0) maxgroups(4)=300
   read(uparams,*) maxgroups(5)
   if(maxgroups(5).eq.0) maxgroups(5)=1
   read(uparams,*) maxgroups(6)
   if(maxgroups(6).eq.0) maxgroups(6)=1
   read(uparams,*) maxgroups(7)
   if(maxgroups(7).eq.0) maxgroups(7)=1
   read(uparams,*) maxgroups(8)
   if(maxgroups(8).eq.0) maxgroups(8)=1
 close(uparams)

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
 read(uparams,*) itemp,tempk
 read(uparams,*) fdfiles , verbose ! number of force-displacement data files

 read(uparams,*) natom_type   ! # of different elements present in prim cell
! allocate(natom(natom_type))
 call allocate_mass(natom_type)
 read(uparams,*) (mas(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) (atname(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) natom_prim_cell        ! # of atoms in primitive cell
 read(uparams,*) fc2flag  ! if zero, take the default range from largest supercell
 allocate(nshells(8,natom_prim_cell))
 nshells(1,1:natom_prim_cell)=1
 read(uparams,*) nshells(2,1:natom_prim_cell)
 read(uparams,*) nshells(3,1:natom_prim_cell)
 read(uparams,*) nshells(4,1:natom_prim_cell)
 read(uparams,*) nshells(5,1:natom_prim_cell)
 read(uparams,*) nshells(6,1:natom_prim_cell)
 read(uparams,*) nshells(7,1:natom_prim_cell)
 read(uparams,*) nshells(8,1:natom_prim_cell)

! just for naming of log file
 open(321,file='dielectric.params',status='old')
 read(321,*)bornflag
 close(321)
 if(bornflag.le.0) then ! no Born correction
    born='B0'
 elseif(bornflag.eq.1) then
    born='BA'  ! A for add
 elseif(bornflag.eq.2) then
    born='BE'  ! E for EwaldE
 else
    born='BS'  ! B for Born (3 or 4)a ; S for subtract
 endif

 if (itemp.eq.0) then
    if(fc2flag.eq.0) then
       it='df'  ! for default
    else
       it='rd'  ! for read
    endif
 elseif(itemp.eq.1) then ! it='T'
    if(fc2flag.eq.0) then
       it='Td'   ! temperature+default
    else
       it='Tr'  ! temperature + read 
    endif
 endif
 write(fdf,'(i1)')fdfiles
 if (fc2flag.eq.1) then
    write(r4,'(i2.2,a1,i2.2,a1,i1)')nshells(2,1),'_',nshells(3,1),'_',nshells(4,1)
 else
    r4='0000000'
    write(r4(3:7),'(a1,i2.2,a1,i1)')'_',nshells(3,1),'_',nshells(4,1)
 endif
 write(incl,'(8i1)')include_fc(:)
 invr='0000'
 if (itrans.ne.0) invr(1:1)='t'
 if (irot.ne.0) invr(2:2)='r'
 if (ihuang.ne.0) invr(3:3)='h'
 if (enforce_inv.ne.0) invr(4:4)='E'
 open(ulog  ,file='log'//fdf//it//born//'_'//r4//'_'//incl//invr//'.dat'   ,status='unknown')

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(tim)
 write(ulog,'(a,3x,a,3x,a,3x,a)')' Program FC234 was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(ulog,'(a,f10.4)')' STARTING TIME OF THE PROGRAM                   IS ',tim
 write(ulog,*)'===================================================================='

 write(ulog,*) ' DEFAULT VALUES:'
 write(ulog,*) svdcut   ,'   cutoff for smallest eigenvalue w to be included for inversion'
 write(ulog,*) tolerance,'   tolerance(Ang) for equating two coordinates '
 write(ulog,*) rcutoff  ,'   cutoff multiple of dist/atom for neighbors of primitive cell '
 write(ulog,*) 'maxterms     =',maxterms
 write(ulog,*) 'maxtermzero  =',maxtermzero
 write(ulog,*) 'maxtermsindep=',maxtermsindep
 write(ulog,*) 'maxgroups    =',maxgroups
 write(ulog,*)' ---------------------------------------------------------'
 write(ulog,*) fc2flag,'  if 0 default range consistent with the largest supercell is chosen for FC2 '
 write(ulog,*) 'itemp,tempk  =',itemp,tempk
 write(ulog,*) 'Whether and How many (not default) shells to include for ranks 2,3,4,5,6,7,8  '
 write(ulog,*) 'rank,incl    = 1',include_fc(1)
       do i=2, maxrank
 write(ulog,*) 'rank,incl,shl=',i,include_fc(i), nshells(i,1:natom_prim_cell)
       enddo
 write(ulog,*) itrans,irot,ihuang,enforce_inv,'  transl, rot and Huang invce, enforcing inv'
 write(ulog,*) fdfiles,' FORCEDISP & POSCAR files to be read'
 write(ulog,*) 'Included fcs and Imposed invariances: trans-rot-Huang=',invr
 write(ulog,*) '--------------------------------------------------'
 write(ulog,*) 'Reading ',natom_type,' atom types with ',natom_prim_cell,'atoms in the primitive cell'

 latticeparameters(1:3) = latticeparameters(1:3)*scal
 call write_out(ulog,'Conventional Lattice parameters a,b,c  ',latticeparameters(1:3))
 call write_out(ulog,'Three angles alpha,beta,gamma (degrees) ',latticeparameters(4:6))
 call write_out(ulog,'Primitive lattice vectors in conventional units ',transpose(primitivelattice))
 call allocate_primcell(natom_prim_cell)
! indices of atoms in primitive cell must be same order as in POSCAR
 counter = 0
 natom = 0
 write(ulog,*)'Atoms: tau, name,type,mass,reduced coordinates in conventional cell '
 do i=1,natom_prim_cell
! positions here must be D format in conventional cell for use by fcs_init
! order does not matter; we use whatever is in this file
    read(uparams,*) tau,atom_type(tau),atompos0(:,tau)
    if (i.ne.tau) then
       print*,' positions must be sorted according to the labels 1,2,3... ',i,tau
       write(ulog,*)' READ_structure: positions must be sorted like labels tau=1,2,3.. ',i,tau
       write(ulog,*)' Check atomic coordinates in primitive cell'
       stop
    endif
    atom0(i)%name = atname(atom_type(tau))
    atom0(i)%at_type = atom_type(tau)
    atom0(i)%tau  = tau
    atom0(i)%mass = mas(atom_type(tau))
    atom0(i)%charge = 0 ! charge is=0 by default unless defined in dielectric.params
!   atom0(i)%nshells  ! actual no of neighbors shells < maxneighbors
!   atom0(i)%shells   ! what is inside each shell
! Make sure the label tau is the same as the atom index in the array atompos
! this is useless because natom(:) is read from POSCAR
    if (atom_type(tau).eq.counter) then
        natom(atom_type(tau)) = natom(atom_type(tau))+1
    else
        natom(atom_type(tau)) = 1
    endif
    counter = atom_type(tau)
    write(ulog,4)atom0(i)%tau,atom0(i)%name,atom0(i)%at_type,atom0(i)%mass,atompos0(:,i)
 enddo

 close(uparams)

   call write_out(6,'latticeparameters',latticeparameters)
   call write_out(6,'primitive lattice',primitivelattice)
   call write_out(6,'# of different atom types ',natom_type)
   call write_out(6,'# of each type in primcel ',natom)
   call write_out(6,'# mass of different types ',mas)
   call write_out(ulog,'# of different atom types ',natom_type)
   call write_out(ulog,'# of each type in primcel ',natom)
   call write_out(ulog,'# mass of different types ',mas)

3 format(a,90(1x,i3))
4 format(i3,1x,a,i3,f8.3,3(1x,f10.5))

 end subroutine read_structure
!=========================================================
 subroutine read_dielectric
!! reads born effective charge and dielectric constant from dielectric.params
 use born
 use lattice
 use atoms_force_constants
 use ios
 use constants, only : r15
 integer i,j,k
 real(r15) asr,aux(3,3)

 open(uborn,file='dielectric.params',status='old')
! allocate(zeu(3,3,natom_prim_cell))

 read(uborn,*) born_flag   ! if 0 use default
 if (born_flag.le.0) then
    write(ulog,*)' Born_flag  =',born_flag,' Born charges are not subtracted'
    write(6   ,*)' Born_flag  =',born_flag,' Born charges are not subtracted'
 elseif(born_flag.eq.1) then
    write(ulog,*)' Born_flag  =',born_flag,' Ewald force will be subtracted from input forces'
 else
    write(ulog,*)' Born_flag  =',born_flag,' non-analytical force will be subtracted from input forces'
    write(6   ,*)' Born_flag  =',born_flag,' non-analytical force will be subtracted from input forces'
 endif

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
!==========================================================
 subroutine read_latdyn
!! reads in ofrmation on kpoint mesh within the full BZ, their shift, temp range, whether classical or quantum
 use ios
 use om_dos
 use params
 use lattice
 use kpoints
 use constants, only : r15
 use atoms_force_constants
 implicit none
 real(r15) junk
 character line*99

  open(uparams,file='latdyn.params',status='old')
  write(6,*) 'READ_PARAMS: opening latdyn.params'
  read(uparams,*)nc
  write(*,*) 'READ_PARAMS: kmesh= ',nc
  nkc=nc(1)*nc(2)*nc(3)
  read(uparams,*)wmesh,wmax ! number of freq meshes for DOS and its upperbound 
  write(*,*) 'READ_PARAMS: wmax = ',wmax
  read(uparams,*)width  ! width of gaussian broadening for DOS
  write(*,*) 'READ_PARAMS: width= ',width
  read(uparams,*) line
  write(*,*) 'READ_PARAMS: just read verbose ',verbose
  read(uparams,*)tmin,tmax,ntemp      ! temps in Kelvin to calculate thermal expansion and other thermal ppties
  write(*,*) 'READ_PARAMS: tmin,tmax,ntemp=',tmin,tmax,ntemp
  read(uparams,*)classical     ! if =1 then use classical distribution (kT/hw)
! read(uparams,*)junk  !iter      ! if=1 then read from a file, else generate  ! sy added iter,split,calk
!  write(*,*)'iter,split,readv3,writev3,calk=',iter,split,readv3,writev3,calc_kappa
! read(uparams,*)junk  !ksub_size  ! each v33sq.xxx.dat file contains V33sq(ksub_size,nkc,ndn,ndn,ndn)
! read(uparams,'(a)')line  ! v3path              ! path to v33sq.dat files
! initialize some default parameters
!  read(uparams,*)calc_cross,q_cross  ! if =1 then calculate the cross section at qcros (reduced U)
!  read(uparams,*)lmicron     ! sample length in microns

  close(uparams)

  write(ulog,*) 'READ_PARAMS: read and kpoints allocated'
  write(ulog,3)'nc1,nc2,nc3=',nc
  write(ulog,3)'nkc        =',nkc
  write(ulog,*)'wmesh, wmax=',wmesh,wmax
  write(ulog,*)'Tmin,Tmax(K=',tmin,tmax
  if (classical .eq. 1) then
     write(ulog,*)'Calculation is done in the CLASSICAL limit'
  else
     write(ulog,*)'Calculation is done in the QUANTUM limit'
  endif

3 format(a,6(1x,i6))
 end subroutine read_latdyn
 !===========================================================
 subroutine read_supercell_vectors(poscar,volum,rm1,rm2,rm3)
!! the following reads the POSCAR file that is used by VASP
!! reads supercell translation vectors and atomic coordinates therein
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 use constants, only : r15
 implicit none
 character, intent(in) :: poscar*(*)
 real(r15), intent(out) :: volum,rm1(3),rm2(3),rm3(3)
 character line*90
 integer i
 real(r15) latt_const,om,a,b,c
 type(vector) pos
 logical exst

 inquire(file=poscar,exist=exst)
 if(exst) then
      open (uposcar,file=poscar,status='old')
 else
      write(ulog,*)' poscar file ',poscar,' does not exist; check your files location and run again'
      write(*   ,*)' poscar file ',poscar,' does not exist; check your files location and run again'
      stop
 endif
 read(uposcar,'(a)') line
 write(ulog,'(a)')' job title on line 1 of POSCAR is: '//line
 read(uposcar,*)lattice_parameter
 read(uposcar,*)rm1
 read(uposcar,*)rm2
 read(uposcar,*)rm3
 if (lattice_parameter .lt.0) then  ! it is the -volume of the cell
    volume_r = -lattice_parameter  ! need to find the volume(r1,r2,r3) for scaling
    call calculate_volume(rm1,rm2,rm3,om)
    latt_const = (volume_r/om)**(1./3.)
 else
    latt_const = lattice_parameter
 endif
 rm1=latt_const*rm1; rm2=latt_const*rm2; rm3=latt_const*rm3
 call calculate_volume(rm1,rm2,rm3,volum)
 close(uposcar)

 end subroutine read_supercell_vectors
 !===========================================================
 subroutine read_supercell(poscar)
!! the following reads the POSCAR file that is used by VASP
!! reads supercell translation vectors and atomic coordinates therein
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 use constants, only : r15
 implicit none
 character line*90, poscar*(*)
 integer i
 real(r15) latt_const,om,a,b,c
 type(vector) pos
 logical exst

! open (uposcar,file='POSCAR',status='old')
 inquire(file=poscar,exist=exst)
 if(exst) then
      open (uposcar,file=poscar,status='old')
 else
      write(ulog,*)' poscar file ',poscar,' does not exist; check your files location and run again'
      write(*   ,*)' poscar file ',poscar,' does not exist; check your files location and run again'
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
 if (line(1:1).eq.'s' .or. line(1:1).eq.'S' ) then  ! this for selective dynamics
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
       atom_sc(i)%equilibrium_pos = latt_const*pos  
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
 call write_out(ulog,'supercell volume ',volume_r)
 call write_out(ulog,'gsc1 ',gs1)
 call write_out(ulog,'gsc2 ',gs2)
 call write_out(ulog,'gsc3 ',gs3)

4 format(9(2x,f19.9))
6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))

 end subroutine read_supercell
!=================================================================
 subroutine identify_atoms_in_supercell
! useless extra check!
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 use constants, only : pi
 use constants, only : r15
 implicit none
 integer i !,n(3)
 real(r15) a(3)

! check mapping and assign label and address to each atom in the supercell
! this seems to be working for now....
 call check_input_poscar_consistency_new

! call write_out(6,' cart_to_prim ',cart_to_prim)
 write(ulog,*)' ASSIGNMENTS WERE DONE FOR ATOMS IN THE SUPER CELL '
! write(ulog,*)'i_sc,type(i),tau(i)      n(i)    ,                r_cart(i)              r_red(i) '
! open(777,file="supercell-1.xyz")
  do i=1,natom_super_cell
!    write(777,3)atom_sc(i)%cell%tau, atom_sc(i)%equilibrium_pos,  &
!&     matmul(cart_to_prim,v2a(atom_sc(i)%equilibrium_pos))
    write(ulog,8)' ',i,atom_sc(i)%at_type, &
&     atom_sc(i)%cell%tau,atom_sc(i)%cell%n,atom_sc(i)%equilibrium_pos,  &
&     matmul(cart_to_prim,v2a(atom_sc(i)%equilibrium_pos))
 enddo
! close(777)

3 format(i5,9(2x,f13.6))
4 format(9(2x,f10.4))
5 format(a,3(i6,2x))
7 format(a,9(1x,f10.4))
8 format(a,i4,i3,i3,'(',3(i2,','),')',9(1x,f10.5))
 write(ulog,*)' POSCAR read successfully and closed'

 end subroutine identify_atoms_in_supercell
!===========================================================
 subroutine check_input_poscar_consistency_new
!! see if all atoms in the supercell can be obtained from the
!! atoms in the input file structure.params using translation vectors
!! of the primitive lattice and assign atom_sc their type and
!! other features/attributes defining their identities.
!! assumes there is no rotation and supercell and primcell are "parallel"
 use atoms_force_constants
 use ios
 use geometry
 use lattice
 use params
 use constants, only : pi
 use constants, only : r15
 implicit none
 integer i,j,k,counter,ier,nt(3),n1,isave
 real(r15) a(3)
 type(vector) shift0,vec
 logical matched,all_matched

! first describe r_i in terms of r0_i, and make sure the lin comb is integer
 write(ulog,*)' ----------------------------------------------------------'
 write(ulog,*)' CHECKING the commensurability between primcell and supercell'
 call check_int(rs1,a,ier,g01,g02,g03)
 write(ulog,4)'checking rs1 in units of r0i,ier=',a,ier
 if (ier .eq. 1) stop
 n_sc(:,1)=nint(a)
 write(ulog,*)' direct coords of rs1=',n_sc(:,1)
 call check_int(rs2,a,ier,g01,g02,g03)
 write(ulog,4)'checking rs2 in units of r0i,ier=',a,ier
 if (ier .eq. 1) stop
 n_sc(:,2)=nint(a)
 write(ulog,*)' direct coords of rs2=',n_sc(:,2)
 call check_int(rs3,a,ier,g01,g02,g03)
 write(ulog,4)'checking rs3 in units of r0i,ier=',a,ier
 if (ier .eq. 1) stop
 n_sc(:,3)=nint(a)
 write(ulog,*)' direct coords of rs3=',n_sc(:,3)

 write(ulog,*)' COMMENSURABILITY CHECKED '
 write(ulog,*)' natom_super_cell=',natom_super_cell
 call write_out(ulog,' N vectors of the supercell ',dble(n_sc))
 call xmatinv(3,dble(n_sc),invn_sc,ier)
 if(ier.ne.0) then
    write(*,*)'n_sc matrix inversion returned error ',ier
    stop
 endif
 call write_out(ulog,' Inverse of n_sc matrix ',invn_sc)

4 format(a,3(2x,f8.3),2x,i1)
7 format(a,2(1x,i4),9(2x,f9.4))
8 format(a,3(1x,i4),9(2x,f9.4))

 write(ulog,*)'FINDING possible translation vectors, trying them on other atoms'
! write(ulog,*)'g01,g02,g03='
! write(ulog,*)g01
! write(ulog,*)g02
! write(ulog,*)g03
! write(ulog,*)' natom_prim_cell=', natom_prim_cell

 checkloop: do i=1,natom_prim_cell !,1,-1
    shift0 = atom0(i)%equilibrium_pos - atom_sc(1)%equilibrium_pos
    write(ulog,7)'atom in PRIMCELL & SC:',i,atom0(i)%at_type,atom0(i)%equilibrium_pos,atom_sc(1)%equilibrium_pos
    write(ulog,3)'trying shift vector to match atom_sc(1)=',shift0

! does it fall on any of the primitive cell atoms?
! try this shift see if all atoms in SC can be mapped to the PC by it.
    all_matched=.True.
    SC: do k = 2,natom_super_cell
       matched = .false.
       PRIM: do j = 1,natom_prim_cell  ! one of the prim-cell atoms has to match
          vec =  shift0 + atom_sc(k)%equilibrium_pos - atom0(j)%equilibrium_pos

! find its direct coordinates on r01,r02,r03

!          call check_int(vec,a,ier,g01,g02,g03)

          a(1)=vec.dot.g01/(2*pi); if(abs(a(1)-nint(a(1))).gt.1d-2)  cycle PRIM
          a(2)=vec.dot.g02/(2*pi); if(abs(a(2)-nint(a(2))).gt.1d-2)  cycle PRIM
          a(3)=vec.dot.g03/(2*pi); if(abs(a(3)-nint(a(3))).gt.1d-2)  cycle PRIM
          write(*,'(4i4,a,3(1x,f9.4),4x,a,3(1x,f9.4))')ier,i,j,k,' er,i0,j0,k_sc; vec=',vec,'vred=',a
          matched=.True.
          exit PRIM

!         if (ier .eq. 0) then
!            matched = .true.
!            write(ulog,'(a,i2,a,i4,a,9(1x,f6.3))')' atom ',j,' in primcell matched atom ',k,' in supercell a=',a,vec
!            write(   *,'(a,i2,a,i4,a,9(1x,f6.3))')' atom ',j,' in primcell matched atom ',k,' in supercell a=',a,vec
         !   cycle SC
! try on other supercell atoms
!            exit PRIM
!         else
!            cycle PRIM
!         endif
       enddo PRIM

       if (.not. matched) then  ! wrong shift , try another shift
           all_matched=.False.
           cycle checkloop     !      exit SC
       endif
    enddo SC

    if(all_matched) then
       write(*,*)'all_matched was true for i0=isave=',i
       write(ulog,*)'all_matched was true for i0=isave=',i
       isave=i
       exit checkloop
    endif
    if (k.ge.natom_super_cell .and. matched) then  ! all atoms were matched
       write(*,*)'all atoms were matched'
       isave=i
       exit checkloop
!   if (matched) exit checkloop
    endif
 enddo checkloop

! shift = r_primcell(i)-r_supercell(1) of direct coordinates a

 if (all_matched) then
    call write_out(ulog,'THE shift vector ', shift0)
    call write_out(ulog,'It maps atom 1 in SC to primitive atom whose number ', isave)
    call write_out(6   ,'THE shift vector ', shift0)
    call write_out(6   ,'It maps atom 1 in SC to primitive atom whose number ', isave)
 else
    write(ulog,*)' NO SHIFT VECTOR COULD BE DEFINED!, check your coordinates again'
    write(6   ,*)' NO SHIFT VECTOR COULD BE DEFINED!, check your coordinates again'
    stop
 endif

! shifting vector has been identified. We can now identify all SC atoms
! atom_sc(1) and atom_prim(i) have the same type and tau
! vector n of supercell atoms is measured wrt i0 (being on atom_sc(1))
                atom_sc(1)%name      = atom0(isave)%name
                atom_sc(1)%cell%n    = 0 !nint(a)  ! a is already integer!
                atom_sc(1)%at_type   = atom0(isave)%at_type
                atom_sc(1)%cell%tau  = atom0(isave)%tau ! should also equal isave
                atom_sc(1)%mass      = atom0(isave)%mass
                atom_sc(1)%charge    = atom0(isave)%charge

      SC2: do k = 2,natom_super_cell
         counter = 0
         matched = .false.
         PRIM2: do j = 1,natom_prim_cell  ! which prim-cell atom does k match?

! vec = atom_sc(k)-atom_sc(1)+atom_prim(isave)-atom_prim(j)  !! isave being fixed
         vec =  atom_sc(k)%equilibrium_pos - atom0(j)%equilibrium_pos+shift0
! get direct coordinates of vec
            call check_int(vec,a,ier,g01,g02,g03)
!           if(verbose) write(ulog,5)'j_prim,k_sc,ier,a=',j,k,ier,a
            if (ier.eq.0) then ! a is integer: j and k are connected by a primitive translation
                matched = .true.
                counter = counter + 1
                if(atom0(j)%tau .ne. iatomcell0(j)) call bomb("tau.ne.iatomcell0!")
                atom_sc(k)%name      = atom0(j)%name
                atom_sc(k)%cell%n    = nint(a) ! a must be integer
                atom_sc(k)%at_type   = atom0(j)%at_type
                atom_sc(k)%cell%tau  = iatomcell0(j) !j  !atom0(j)%tau =
                atom_sc(k)%mass      = atom0(j)%mass
                atom_sc(k)%charge    = atom0(j)%charge
                if(verbose) write(ulog,2)k,' in SC has attributes type,tau,n,r_1k-r5,reduced ', &
&                   atom_sc(k)%at_type ,atom_sc(k)%cell%tau ,atom_sc(k)%cell%n ,vec , &
&                   matmul(cart_to_prim,v2a(vec))
            !    exit PRIM2
                cycle SC2  ! go to the next k
            endif
! check k with the next j in primitive cell
         enddo PRIM2

         if (counter.eq. 0) then
            call write_out(ulog,'SC ATOM which can not be mapped to the prim cell ',k)
            write(ulog,*)' check your coordinates '
            write(ulog,3)atom_sc(k)%equilibrium_pos
            stop
         elseif ( counter .ge. 2) then
            write(ulog,*)' CHECK: more than one atom matched! check your coordinates'
            write(*   ,*)' CHECK: more than one atom matched! check your coordinates'
            stop
         endif
       enddo SC2

 write(ulog,*)'*******************************************************'
 write(ulog,*)'SINCE THE PROGRAM WAS NOT STOPPED, MAPPING-CHECK PASSED!'
 write(ulog,*)'NOW CHECKING CORRESPONDANCE WITH ATOMPOS ARRAY'
 write(ulog,*)'*******************************************************'

 write(*,*)'calling write_correspondance'
 call write_correspondance (isave,shift0)

2 format(i5,a,2(1x,i2),'(',3(i2),')',2(2x,3(1x,f8.3)))
3 format(9(1x,g13.6))
5 format(a,1x,3(4x,i9),9(1x,f9.4))
6 format(a,i2,1x,2i5,9(1x,f9.4))

 open(173,file='poscar.xyz')
 write(173,*) natom_super_cell
 write(173,*)'# poscar.xyz to visualize: name coords, type, n, tau'
 do i=1,natom_super_cell

    n1 =atom_sc(i)%at_type
    nt =atom_sc(i)%cell%n
    k  =atom_sc(i)%cell%tau
    if(n1.eq.1)  then
       write(173,9)'Ga ',atom_sc(i)%equilibrium_pos,n1,nt,k
    elseif(n1.eq.2) then
       write(173,9)'S  ',atom_sc(i)%equilibrium_pos,n1,nt,k
    elseif(n1.eq.3) then
       write(173,9)'P  ',atom_sc(i)%equilibrium_pos,n1,nt,k
    elseif(n1.eq.4) then
       write(173,9)'Ge ',atom_sc(i)%equilibrium_pos,n1,nt,k
    elseif(n1.eq.5) then
       write(173,9)'Si ',atom_sc(i)%equilibrium_pos,n1,nt,k
    endif

 enddo

 close(173)
9 format(a,3(1x,f12.5),3x,i3,'(',3i2,')',i5)

 end subroutine check_input_poscar_consistency_new
!===========================================================
 subroutine count_configs(outcar,ncfg)
 use ios , only : ulog,utraj
! use atoms_force_constants, only : natom_super_cell
 implicit none
 integer, intent(out) :: ncfg
 character, intent(in):: outcar*(*)
 integer t,j,i
 character line*99
 logical found,exst

 write(*,*)' opening FORCEDISP file'

 inquire(file=outcar,exist=exst)
 if(exst) then
      open(utraj,file=outcar,status='old')
 else
      write(ulog,*)' FORCEDISP file ',outcar,' does not exist; check your files location and run again'
      stop
 endif

! first find the number of configurations: ncfg
 t=0
 do j=1,11000000
    read(utraj,'(a)',end=99)line
    call findword('POSITION',line,found,i)
    if (found) then
       t = t+1
    else
    call findword('vasprun',line,found,i)
       if (found) then
          t = t+1
       endif
    endif
 enddo
99 write(ulog,*)' reached the end of FORCEDISP file; number of configurations= ',t
 close(utraj)
 ncfg = t
 write(*,*)' reached the end of FORCEDISP file which had ',ncfg,' configurations '
 if(t.eq.0) then
    write(ulog,*)' the word POSITION was not found in FORCEDISP file, check it!'
    stop
 endif

! i=(j-1)/t - 2  ! this is the number of atoms read from FORCEDISP
! if (i .ne. natom_super_cell ) then
!    write(ulog,*)' number of atoms read .ne. no of atoms in POSCAR file',i,natom_super_cell
!    write(ulog,*)' # of read lines in FORCEDISP is=',j-1
!    write(ulog,*)' check your POSCAR and FORCEDISP  again '
!    write(ulog,*)' make sure the # of atoms is the same in both files'
!    write(ulog,*)' there should be no blank lines at the end of FORCEDISP'
!    stop
! endif

 end subroutine count_configs
!===========================================================
 subroutine read_force_position_data(outcar,ncfg,engy,dsp,frc,nlin2)
!! reads contents of FORCEDISP; second lines are energies, frc_constr=3*natom_super_cell*ncfg
!! outputs are dsp(3,NSC,ncfg),frc(3,NSC,ncfg),engy(ncfg) for the FORCEDISP file read
 use ios
 use atoms_force_constants
 use geometry
 use params
 use lattice
 use svd_stuff
 use constants, only : r15
 implicit none
 integer, intent(in) :: ncfg
 integer, intent(out) :: nlin2(50000)
 character, intent(in):: outcar*(*)
 integer i,t,j,k,frm ! format of the outcar file
 real(r15), save :: emin
 logical, save :: first_call=.true.
! real(r15), allocatable, intent(inout) :: engy(:),dsp(:,:,:),frc(:,:,:)
 real(r15), intent(out) :: engy(ncfg),dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg)
 character line*99
 logical found,exst
 real(r15) dr(3),dc(3)

 write(*,*)' Re Opening FORCEDISP file'

 inquire(file=outcar,exist=exst)
 if(exst) then
    open(utraj,file=outcar,status='old')
 else
    write(ulog,*)' FORCEDISP file ',outcar,' does not exist; check your files location and run again'
    stop
 endif


! now get the FORCES from FORCEDISP file
 t=0 ; engy=0
 do j=1,ncfg
    read(utraj,'(a)',end=88)line
!    if (line(1:8) .eq. "POSITION" ) then
    call findword('POSITION',line,found,i)
    if (found) then
       frm=1  ! my format
       t = t+1
       read(utraj,*) k,engy(t)  ! start with 1 since t=0
       do i=1,natom_super_cell
           read(utraj,*) dsp(1:3,i,t),frc(1:3,i,t)
!          dsp(1:3,i,t)=dsp(1:3,i,t)-atom_sc(i)%equilibrium_pos
           dc=dsp(:,i,t)-v2a(atom_sc(i)%equilibrium_pos)
! eventually subtract any supercell translation vectors
        call cart_to_direct(dc,dr,gs1,gs2,gs3)
        dr(1) = dr(1) - anint(dr(1))
        dr(2) = dr(2) - anint(dr(2))
        dr(3) = dr(3) - anint(dr(3))
        call direct_to_cart(dr,dsp(:,i,t),rs1,rs2,rs3)
       enddo
       nlin2(t)=3*natom_super_cell
    else
       call findword('vasprun',line,found,i)
       if (found) then
          frm=2 ! Onishi format
          t=t+1
          call findword('(eV):',line,found,i)
          read(line((i+5):),*)engy(t)
          do i=1,natom_super_cell
             read(utraj,*) dsp(1:3,i,t),frc(1:3,i,t)
             dsp(1:3,i,t)=dsp(1:3,i,t)*ab !+atom_sc(i)%equilibrium_pos
             frc(1:3,i,t)=frc(1:3,i,t)*ryd
          enddo
          nlin2(t)=3*natom_super_cell
       endif
    endif
    write(*,*)'j=t,nlines,engy=',j,t,nlin2(t),engy(t)
    if(t.ne.j ) then !.or. t-1.ne.k) then
       write(*,*)'Error in reading snapshot#s in FORCEDISP?',j,t
    endif

 enddo
88 write(ulog,*)' reached the end of FORCEDISP file after ',t,' steps'
 close(utraj)
 write(ulog,*)' last line read was '
 write(ulog,'(a)') line
 close(utraj)
 if (t .ne. ncfg) then
    write(*,*)'ERROR in reading the force file FORCEDISP'
    write(*,*)'ncfg, # of read steps=',ncfg,t
    stop
 endif

 if(verbose) then
   write(ulog,*)'writing last snapshot frc and positions in the log file: 3xN_sc'
   call write_out(ulog,' last force ',frc(:,:,t))
   call write_out(ulog,' last coord ',dsp(:,:,t))
 endif

! get energy per primitive unit cell so that it does not depend on supercell size
! engy=engy/natom_super_cell*natom_prim_cell
!! energy in the exponent of Boltzmann should be proportional to supercell size

! subtract lowest energy value ! get it from FORCEDISP1; it is arbitrary anyways
 if( first_call ) then
    emin=minval(engy)
    first_call=.false.
 endif
 engy=engy-emin  ! this is for that supercell of fixed size; how to compare among different SCs?
 engy=engy/natom_super_cell*natom_prim_cell  ! to normalize among supercells of different size

 write(*,9)'Energy assigned ',engy
 write(ulog,9)'Energy assigned ',engy
 write(*,*)'Calling calculate_and_write_displacements'
 write(ulog,*)'Now writing all snapshots force and displacement in the log file'
 call calculate_and_write_displacements(ncfg,dsp,frc)
! input and output: displ, including atom_sc%equilibrium_pos

 write(*,*)'exiting read_force_position_data '
9 format(a,200(1x,g11.4))

 end subroutine read_force_position_data
!===========================================================
 subroutine write_independent_fcs(n,sig,sd,ulog)
!! write independent or irreducible force constants in file fci_irr.dat (i=rank)
 use svd_stuff
 use atoms_force_constants
 use params
 use constants, only : r15
 implicit none
 integer, intent(in) :: ulog,n
 real(r15), intent(in) :: sig(n)
 real(r15), intent(out):: sd(8)
 integer rnk,k,g,cnt2,cnt,ti,cnt4


   write(ulog,*)' indepterm, fcs, sigma, normalized sigma for all SVD'
   do k=1,n
      write(ulog,4) k,fcs(k),sig(k),sig(k)/sqrt(1.*dim_al)
   enddo
   write(ulog,*) "rank, average error"

!  k1 ----------------- for rank=2 we always work with size_kept_fc2
 if(fc2flag.eq.0) then
   write(ulog,*)'# of groups of rank 2, ngr, versus # of kept indep terms =',map(2)%ngr,sum(keep_fc2i)
   write(ulog,*)'# of indep terms versus # of kept terms=',map(2)%ntotind,size_kept_fc2
 endif
!  k1 -----------------

   write(ulog,*)' dim_ac, size(fcs) =',dim_ac,n

   cnt2=0;sd=0
   do rnk=1,maxrank

      write(ulog,*)'rank,# of groups=',rnk,map(rnk)%ngr
      if(include_fc(rnk).ne.1) cycle
      cnt=0
      do g=1,map(rnk)%ngr
   !    if(rnk.eq.2) then
   !       if( keep_fc2i(g).ne.1) cycle
   !    endif
        do ti=1,map(rnk)%ntind(g)
            cnt=cnt+1   ! counter of terms of given rank
            if(rnk.eq.2) then
               cnt=counter2(g,ti)
         !     cnt4=ti
         !     if(g.gt.1) cnt4=sum(map(2)%ntind(1:g-1)) +ti ! cumulative index up to ti in group g 
               cnt2 = map(1)%ntotind+ current2(g,ti) !sum(keep_fc2i(1:cnt4)) 
            else 
               cnt2=cnt2+1 ! cumulative counter over all ranks
            endif
            sd(rnk)=sd(rnk)+sig(cnt2)/sqrt(1.*dim_al)
            if(verbose) write(ulog,6)'rank,group,nind,cnt(rnk),cntot,sd(rnk)=' &
&                       ,rnk,g,ti,cnt,cnt2,sig(cnt2)/sqrt(1.*dim_al),sd(rnk)
        enddo
      enddo
      if(cnt .ne. 0 ) sd(rnk)=sd(rnk)/cnt
      write(ulog,3)'Rank, sd(rank)=', rnk,sd(rnk)

   enddo
   write(ulog,5)' !==== Error summary for each rank',(sd(rnk),rnk=1,4)


3 format(a,i6,3(2x,g14.7))
4 format(i6,3(2x,g14.7))
5 format(a,9(2x,g11.4))
6 format(a,5i4,9(2x,g11.4))

 end subroutine write_independent_fcs
!===========================================================
 subroutine write_neighbors
 use atoms_force_constants
 use params
 use ios
 use lattice, only : g01,g02,g03,check_int  !r01,r02,r03,
 use geometry
 use constants, only : r15
 implicit none
 integer, parameter :: mesh=600,mx=60
 integer i0,j0,shel_count,j,nm(3),n5(3),ta,nbmx,jj,ier
 real(r15) dij(natom_prim_cell,mx),rr(3),eps(3)
 real(r15) ds(natom_prim_cell,mesh),dmesh(mesh),w0(mx),dmax

 write(*,*)'mx=',mx
 write(ulog,*)'mx=',mx

 dmax=15d0
 do j=1,mesh
    dmesh(j)=dmax*j/dble(mesh)
 enddo
 open(345,file='pairs.dat')
 eps=1d-5
 dij=100;w0=1d0
 do i0=1,natom_prim_cell
    write(ulog,*)' ******************************'
    write(ulog,*)' Neighbors of atom number ',i0,' up to nshells=',atom0(i0)%nshells
    write(ulog,*)' ******************************'
!   nsh=size(atom0(1)%shells(:)%radius)
    write(*   ,*)' size  of shells=',atom0(i0)%nshells
    write(ulog,*)' i0, # of shells=',i0,atom0(i0)%nshells

    do shel_count=1,min(atom0(i0)%nshells , mx)
       write(*,*)' i0,mx,shell_count= ',i0,mx,shel_count
       if(shel_count.gt. atom0(i0)%nshells) cycle
       nbmx = atom0(i0)%shells(shel_count)%no_of_neighbors
       dij(i0,shel_count)  = atom0(i0)%shells(shel_count)%radius
       if(shel_count.eq.1) then
          write(ulog,7)'# i0,shel#,ngbrs ,radius    =',i0,shel_count,nbmx,dij(i0,shel_count)
          write(345,7)'# i0,shel#,ngbrs ,radius    =',i0,shel_count,nbmx,dij(i0,shel_count)
       else
          write(ulog,7)'# i0,shel#,ngbrs ,radius,gap=',i0,shel_count,nbmx,dij(i0,shel_count),dij(i0,shel_count)-dij(i0,shel_count-1)
          write(345,7)'# i0,shel#,ngbrs ,radius,gap=',i0,shel_count,nbmx,dij(i0,shel_count),dij(i0,shel_count)-dij(i0,shel_count-1)
       endif
       do j=1,min(nbmx,500)
          ta =  atom0(i0)%shells(shel_count)%neighbors(j)%tau
          nm =  atom0(i0)%shells(shel_count)%neighbors(j)%n
          jj =  atom0(i0)%shells(shel_count)%neighbors(j)%atomposindx
 !        write(ulog,6)' j,jj,tauj,nij=',j,jj,ta,nm
! consistency test:i must have tau,nm correspond to that atompos coordinates
          call check_int(a2v(atompos(:,jj)+eps),rr,ier,g01,g02,g03)
          n5=floor(rr)
          if (length(n5-nm).gt.1d-5) then
             write(*,*)' Error: n,tau does not correspond to floor(red(atompos)) ',nm,n5
             write(*,*)' iatomcell:tau,n=',iatomcell0(jj),iatomcell(:,jj)
             write(*,5)' pos(i)     = ',atom0(i0)%equilibrium_pos
             write(*,5)' ni(j)*r0i  = ',rr
             write(*,5)' atompos(j) = ',atompos(:,jj)
             write(*,6)' i0,n,tauj,j=',i0,nm,ta,jj
             write(*,5)' reduced coords of atom j=',rr
          !  stop
          endif
       enddo
    enddo
    call calculate_dos(mx,dij(i0,:),w0,mesh,dmesh,ds(i0,:))
 enddo

 do j=1,mesh
    write(345,9)dmesh(j),(ds(i0,j),i0=1,natom_prim_cell)
 enddo
 close(345)

 write(ulog,*)' ************ End of the neighbors list ************** '
2 format(a,i4,2x,f8.4,2x,i4)
3 format(a,2x,i3,2x,i3,2x,f10.4,2x,3i5,4x,' (',3(1x,i4),')')
4 format(a,2x,i3,2x,i3,2x,9f10.4)
5 format(a,9(1x,f9.4))
6 format(a,2x,3i5,4x,' (',3(i2),')')
7 format(a,2x,3i5,2x,99(1x,f8.3))
9 format(99(1x,f8.3))
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
 integer rnk,t,ti,i,res,j,rs,k
 integer iat(maxrank),ixyz(maxrank),g,ng,term,term2,cnt2,frm,cnt3,ntind(maxrank),ngroup(maxrank)
 real(r15) rij,bunit,one,fcd,trace_fc,dij 
! character frmt*2,goh*48,ln*1,geh*47
 character frmt*2,goh*60,ln*1,geh*60,lm*1

 bunit = ryd/ab/ab
 one =1d0

6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(2(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
9 format(9(2x,f19.10))

! first write the crystal data
 ngroup= map(:)%ntot
 ntind=map(:)%ntotind
 call write_lat_fc(ntind,ngroup,'lat_fc.dat') !map(:)%ntotind,map(:)%ntot)

!----------------------------------------
 res = 0
 ranks: do rnk=1,maxrank
  if ( include_fc(rnk) .eq. 1 ) then
    frm=30+rnk
    write(ln,'(i1)')rnk
    write(lm,'(i1)')rnk-1
    goh='(a1,i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g14.7,1x,f7.4,1x,15i3)'
    geh='(   i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g14.7,1x,f7.4,1x,15i3)'
!   write(*,*)'for rank ',rnk,' formats geh and goh are:'
!   write(*,*)geh
!   write(*,*)goh

!   geh='(i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g14.8,f8.4,2x,f9.5)'
    write(frmt,'(i2)')30+rnk
    write(ulog,*)' FOR RANK=',rnk,' format=',frmt
    write(*,*)' FOR RANK=',rnk,' format=',frmt
    write(ufc1-1+rnk,*)'# RANK ',rnk,' tensors :term,group,(iatom,ixyz)_2 d^nU/dx_{i,alpha}^n'

    ng=map(rnk)%ngr ! number of groups
    cnt2=0  ! cumulative # of indep terms of given rnk up to g and excluding group g
    term = 0
    term2= 0
    groups: do g=1,map(rnk)%ngr  ! index of a given group of given rank
         
       cnt2=0
       if(rnk.eq.2) then
          if(g.gt.1) cnt2=sum(keep_fc2i(1:counter2(g-1,map(2)%ntind(g-1))))
       else
          if(g.gt.1) cnt2 = sum(map(rnk)%ntind(1:g-1)) ! cumulative indep term index up to, but excluding  group g (used for ranks other than 2)
       endif
!         if(g.gt.1) s1=sum(ind(1:g-1)) 
!         cnt2 = sum(keep_fc2i(1:s1)) 
 !        cnt3=0
 ! write in the log and fcn_fit.dat: cnt2+ti is the position of indep_fc of that rank
       do ti=1,map(rnk)%ntind(g)  ! index of independent terms in that group g
 !        cnt3=cnt3+1
          iat(1:rnk)  = map(rnk)%gr(g)%iatind (:,ti)
          ixyz(1:rnk) = map(rnk)%gr(g)%ixyzind(:,ti)
          term = term+1
          if (rnk.eq.2) then
            rij = length(atompos(:,iat(1))-atompos(:,iat(2)))
          else
            rij = 0
          endif
!         write(ulog,goh) map(rnk)%err(cnt2+ti),g,ti,(iat(j),ixyz(j),j=1,rnk),  &
!      &     fcs(res+cnt2+ti),fcs(res+cnt2+ti)/ryd*ab**rnk,rij
!         write(ufit1-1+rnk,geh) ti,g,(iat(j),ixyz(j),j=1,rnk),  &
!      &     fcs(res+cnt2+ti),one,rij
          if(rnk.eq.2) then 
             if(keep_fc2i(counter2(g,ti)).eq.0) cycle
             write(ulog,goh) map(rnk)%err(current2(g,ti)),g,ti,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+current2(g,ti)),rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk)
             write(ufit1-1+rnk,geh) ti,g,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+current2(g,ti)),rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk)
          else
             write(ulog,goh) map(rnk)%err(cnt2+ti),g,ti,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+cnt2+ti),rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk)
!      &     fcs(res+cnt2+ti),fcs(res+cnt2+ti)/ryd*ab**rnk,rij
             write(ufit1-1+rnk,geh) ti,g,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+cnt2+ti),rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk)
!      &     fcs(res+cnt2+ti),one,rij
          endif
       enddo

    ! write in the fcn.dat file
       do t=1,map(rnk)%nt(g)  ! index of dependent terms in that group g
          iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
          ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
          term2= term2+1
          fcd = 0
  !       cnt3=0
       ! must find the corresponding index of the indep term t <-> ti
          do ti=1,map(rnk)%ntind(g)
  !          cnt3=cnt3+1
         ! this is the index of the indep FC coming in the A*FC=b matrix product
             if(rnk.eq.2) then
                if(keep_fc2i(counter2(g,ti)).eq.0) cycle
                fcd = fcd + fcs(res+current2(g,ti))*map(rnk)%gr(g)%mat(t,ti) !*keep_fc2i(counter2(g,ti))
             else
                fcd = fcd + fcs(res+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
             endif
          enddo
          if(rnk.eq.2) then
             write(ufc1-1+rnk,geh)t,g, (iat(j),ixyz(j),j=1,rnk),fcd,  &  !one
&              rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk)
          else
!          if( abs(fcd) .gt. margin) then
             write(ufc1-1+rnk,geh)t,g, (iat(j),ixyz(j),j=1,rnk),fcd,  &  !one
&              rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk)
!          endif
          endif
       enddo
  enddo groups

     if (fc2flag.eq.0 .and. rnk.eq.2) then
        res=res+ size_kept_fc2
     else
        res = res+map(rnk)%ntotind
     endif

  endif
 enddo ranks

write(ulog,*)'******* Trace for the harmonic FCs ********'
 open(456,file='trace_fc.dat')
 write(456,*)'# dij,trace_fc,tau_i,tau_j,j,atompos(:,i),atompos(:,j)'

! write the trace of FC2
 rnk=2; res=map(1)%ntotind
 iloop: do i=1,natom_prim_cell
 jloop: do j=1,natoms
     rij = length(atompos(:,i)-atompos(:,j))
     trace_fc=0
     rs=map(1)%ntotind
     if ( include_fc(rnk) .ne. 0 ) then
        ng=map(rnk)%ngr ! number of groups
  !     cnt2=0
        do g=1,map(rnk)%ngr  ! index of a given group
  !        if(g.gt.1) cnt2 = sum(keep_fc2i(1:sum(map(2)%ntind(1:g-1))))
           do t=1,map(rnk)%nt(g)  ! index of independent terms in that group g
              iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
              ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
              if (iat(1).ne.i .or. iat(2).ne.j) cycle !iloop
              fcd = 0
  !           cnt3=0
       ! must find the corresponding index of the indep term t <-> ti
              do ti=1,map(rnk)%ntind(g)
       !         cnt3=ti
       !         if(g.gt.1) cnt3=sum(map(2)%ntind(1:g-1))+ti
       ! this is the index of the indep FC coming in the A*FC=b matrix product
       !         fcd = fcd + fcs(res+cnt2+cnt3)*map(rnk)%gr(g)%mat(t,ti)*keep_fc2i(cnt2+ti)
                 fcd = fcd + fcs(res+current2(g,ti))*map(rnk)%gr(g)%mat(t,ti) !*keep_fc2i(counter2(g,ti))
                 dij = length(atompos(:,iat(1))-atompos(:,iat(2)))
              enddo
              if (ixyz(1).eq.ixyz(2)) then
                 trace_fc = trace_fc + fcd
              endif
           enddo
        enddo
! write(ulog,4)'i0,j,g,t,fcd,tr=',i,j,g,t,fcd,trace_fc
        if(abs(trace_fc).gt.1d-5) then  ! .and. ixyz(1).eq.1 .and. ixyz(2).eq.1) then
           write(ulog,8) i,j,dij,trace_fc
           write(456,3) dij,trace_fc,i,iatomcell0(j),j,atompos(:,i),atompos(:,j)
        endif
     endif
  enddo jloop
  enddo iloop

  close(456)

  write(ulog,*)'***************** END OF FC Trace ******************'

2 format(9(1x,f10.4))
3 format(2(1x,f10.4),3i5,9(1x,f10.4))
4 format(a,4i4,9(1x,f10.4))

! if (res.ne.nindepfc) then
!    write(ulog,*)'WRITE_OUTPUT_FCS: sum(nterms),ngr=',res,nindepfc
!    write(ulog,*)'WRITE_OUTPUT_FCS: the difference should be ',map(2)%ntotind-sum(keep_fc2i)
! endif

 end subroutine write_output_fcs
!============================================================
 subroutine write_springs
 use svd_stuff
 use ios
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
 integer rnk,t,ti,i,res,j,rs,k
 integer iat(maxrank),ixyz(maxrank),g,ng,term,term2,cnt2,frm,ntind(maxrank),ngroup(maxrank)
 real(r15) rij,bunit,one,fcd,trace_fc,dij , zero(3)
 character lm*1

 one =1d0
 zero = 0

 do i=1,natom_prim_cell
 enddo

 rnk=2; res=map(1)%ntotind
 iloop: do i=1,natom_prim_cell
    write(lm,'(i1)')i
    open(875+i,file='springs'//lm//'.dat')
    write(875+i,*)"# 0 0 0 ,atompos(:,iat(1))-atompos(:,iat(2)) , d12 "
 jloop: do j=1,natoms
     rij = length(atompos(:,i)-atompos(:,j))
     rs=map(1)%ntotind
     if ( include_fc(rnk) .ne. 0 ) then
        ng=map(rnk)%ngr ! number of groups
        term2= 0
        do g=1,map(rnk)%ngr  ! index of a given group
           do t=1,map(rnk)%nt(g)  ! index of independent terms in that group g
              iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
              ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
              if (iat(1).ne.i .or. iat(2).ne.j) cycle !iloop
              dij = length(atompos(:,iat(1))-atompos(:,iat(2)))
              do ti=1,map(rnk)%ntind(g)
     !           if(map(2)%gr(g)%mat(t,ti).ne.0.and.keep_fc2i(counter2(g,ti)).eq.1 .and.ixyz(1).eq.1.and.ixyz(2).eq.1) then
                 if(keep_fc2i(counter2(g,ti)).eq.1 .and.ixyz(1).eq.1.and.ixyz(2).eq.1) then
                     write(875+i,2)zero ,atompos(:,iat(1))-atompos(:,iat(2)),dij
                 endif
              enddo
           enddo
        enddo
     endif
  enddo jloop
     close(875+i)
  enddo iloop

2 format(9(1x,f10.4))
 end subroutine write_springs
!============================================================
 subroutine read_fcs_3(iunit,fn,rank,fc,nfc)
!! reads irreducible force constants of given rank from file fn if the latter exists
 use svd_stuff
 use ios
 use atoms_force_constants
 use params
 use constants, only : r15
 implicit none
 integer, intent(in) :: rank,iunit,nfc
 integer t,ti,igr,term,g,j
 real(r15) fc(nfc)
 character fn*11
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

     if (nfc.ne.map(rank)%ntotind) then
        write(ulog,*)' READ_FCS_3: rank=',rank
        write(ulog,*)' number of groups ',nfc,' is not the same as one by setup_maps ',map(rank)%ntotind
        stop
     endif

     term = 0;
     groups: do g=1,map(rank)%ngr ! %ntotind  ! index of a given group
       map(rank)%gr(g)%iatind(:,:)=0
       map(rank)%gr(g)%ixyzind(:,:)=0
       do ti=1,map(rank)%ntind(g)  ! index of independent terms in that group g
            term=term+1
            read(iunit,*,err=91) t,igr,(map(rank)%gr(igr)%iatind (j,t),  &
            &                           map(rank)%gr(igr)%ixyzind(j,t), j=1,rank), &
            &                           fc(term)
       enddo
     enddo groups

91   write(ulog,*)'READ_FCS_3: rank=',rank,' reached end of file after ',term-1,' lines'
     write(ulog,3)fc
     if (term.eq.1) then
        write(ulog,*)'READ_FCS_3: the file exists but is empty! rank=',rank
        write(*   ,*)'READ_FCS_3: the file exists but is empty! rank=',rank
        stop
!       include_fc(rank) = 0       ! if nothing, then exclude from fitting
     endif
     close(iunit)

3 format(99(1x,f9.5))
 end subroutine read_fcs_3
!============================================================
 subroutine read_fcs(iunit,fn,rank,fc2,nfc)
 use svd_stuff
 use ios
 use atoms_force_constants
 use params
 use constants, only : r15
 implicit none
 integer rank,iunit,t,res,i,a, b, cnt2,term,g,ti,nfc
 real(r15) fc2(nfc)
 character(len=*), intent(in) :: fn
 logical ex

 open(iunit,file=fn,status='old')

! read(iunit,'(a)') line
 t = 0
!----------------------------------------
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

 res =  igroup_1(nterms(1))

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
!          if(keep_fc2i(g).ne.1) cycle
!       endif
!! K1 new change made on 2/13/23 ----------------

       if(g.gt.1 .and. (keep_fc2i(g-1).eq.1)) cnt2=cnt2+map(2)%ntind(g-1)

       do t=1,map(2)%ntind(g)  ! index of independent terms in that group g
          term = term+1
           read(ufit2,*) ti,nfc,map(2)%gr(nfc)%iatind (1,ti),  &
           &                   map(2)%gr(nfc)%ixyzind(1,ti),  &
           &                   map(2)%gr(nfc)%iatind (2,ti),  &
           &                   map(2)%gr(nfc)%ixyzind(2,ti),  &
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
!          if(keep_fc2i(g).ne.1) cycle
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
! 
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

 end subroutine read_fcs
!=================================================================
 subroutine calculate_and_write_displacements(ncfg,dsp,frc)
!! equilibrium position already subtracted when reading the FORCEDISP files
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 use constants, only : r15
 implicit none
 integer, intent(in) :: ncfg
 real(r15), intent(inout) :: dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg)
 integer i,t,step
 real(r15) dc(3),dr(3)

 write(ulog,*)' t , particle#, cartesian disp & forces u1,u2,u3,f1,f2,f3'
 write(ulog,*)' Number of configurations NCFG=',ncfg
! get displacements from positions ! this is not needed!
!  do t=1,ncfg
!     do i=1,natom_super_cell
! ! first get direct coordinates, then take the distance using pbc, then
! ! transfrom back to cartesian coordinates
!        dc(1) = dsp(1,i,t) !- atom_sc(i)%equilibrium_pos%x
!        dc(2) = dsp(2,i,t) !- atom_sc(i)%equilibrium_pos%y
!        dc(3) = dsp(3,i,t) !- atom_sc(i)%equilibrium_pos%z
!        call cart_to_direct(dc,dr)
!        dr(1) = dr(1) - anint(dr(1))
!        dr(2) = dr(2) - anint(dr(2))
!        dr(3) = dr(3) - anint(dr(3))
! !      call direct_to_cart_aa(dr,dsp(:,i,t))
!        call direct_to_cart(dr,dsp(:,i,t))
! ! add back the positions
! !       dsp(1,i,t) = dsp(1,i,t) + atom_sc(i)%equilibrium_pos%x
! !       dsp(2,i,t) = dsp(2,i,t) + atom_sc(i)%equilibrium_pos%y
! !       dsp(3,i,t) = dsp(3,i,t) + atom_sc(i)%equilibrium_pos%z
!     enddo
!  enddo
 step=ncfg-1

 write(*,*)'ncfg,step=',ncfg,step
 if(verbose) step=1 ! every step, otherwise only the fist and last snapshots will be written
 if (step.eq.0) step=1
 write(*,*)'ncfg,step=',ncfg,step

 do t=1,ncfg ! write all steps ,step           ! write first and last displ-force data
    do i=1,natom_super_cell
       write(ulog,6)t,i,dsp(:,i,t),frc(:,i,t)
    enddo
 enddo

 6 format(2(i5),3(1x,f10.5),3x,3(1x,g12.5))

 end subroutine calculate_and_write_displacements
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
        else if(n.lt.10000)then
          write(lineout(m+1:m+4),'(i4)')n
          m=m+4
        else if(n.lt.100000)then
          write(lineout(m+1:m+5),'(i5)')n
          m=m+5
        else if(n.lt.1000000)then
          write(lineout(m+1:m+6),'(i6)')n
          m=m+6
        else if(n.lt.10000000)then
          write(lineout(m+1:m+7),'(i7)')n
          m=m+7
        else if(n.lt.100000000)then
          write(lineout(m+1:m+8),'(i8)')n
          m=m+8
        endif
      enddo
!     write(*,*)' Exiting ustring'

      end subroutine ustring
!============================================================
 subroutine pos_out_consistency
!! checks consistency between POSCAR:atom_sc%equilibrium_pos and FORCEDISP:displ(:,:,1)
!! allows a max difference of 0.2 ang per direction
 use geometry
 use ios, only : ulog
! use params
 use lattice
 use constants, only : r15
 use atoms_force_constants
 implicit none
 integer i
 real(r15) dr(3),dc(3)
 do i=1,natom_super_cell

    dc = atom_sc(i)%equilibrium_pos - displ(:,i,1)
    call cart_to_direct(dc,dr,gs1,gs2,gs3)
    dr(1) = dr(1) - anint(dr(1))
    dr(2) = dr(2) - anint(dr(2))
    dr(3) = dr(3) - anint(dr(3))
    call direct_to_cart(dr,dc,rs1,rs2,rs3)

!   if (.not.( dc .myeq. 0d0 ) ) then
!   if ( abs(dc(1)).gt.tolerance .or. abs(dc(2)).gt.tolerance  &
!   &                            .or. abs(dc(3)).gt.tolerance ) then
    if ( abs(dc(1)).gt.0.2 .or. abs(dc(2)).gt.0.2 .or. abs(dc(3)).gt.0.2 ) then
       write(ulog,*)' atom # ',i,' in POSCAR and 1st snapshot of FORCEDISP are different'
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
 subroutine write_correspondance(isave,shift0)
!! writes correspondance between atompos pair assigned to a spring in the chosen model and supercell pair
!! in which the first atom is in the primitive cell defined in structure.params, and the second one in poscar
!! isave is the atom0 identified with atom_sc(1), shift0=atom0(isave)-atom_sc(1) 
 use params
 use atoms_force_constants
 use svd_stuff
 use ios
 use geometry
 use lattice
 use constants, only : r15
 implicit none
 integer, intent(in) :: isave
 type(vector), intent(in) :: shift0   ! atom0(isave)=atom_sc(1)
 integer t,nat,j_sc,g,rnk,i0
 integer al,be,taui,tauj,ni(3),nj(3),j ,iat
 real(r15) rij,rn1(3),rsc(3)

 rnk = 2
 write(ulog,*)' HARMONIC FCs: correspondance between a given FC and all pairs '
 write(ulog,*)' Which atom pairs in the supercell are bridged by a FC element '
! write(ulog,*)'# iatom, [  taui(nsi) ],  alpha_i=1 ;  j, [  tauj(nsj) ], alpha_j=1  :   term group  rij'
 write(ulog,*)'# i taui ; j tauj  ; j_sc   term  group  rij'
 write(ucor,*)'# iatom0,[ taui(ni) ],alpha_i ; jatom,[ tauj(nj) ],alpha_j ; j_sc  term group  rij'

  do nat=1,natom_prim_cell

!    rn1=v2a(atom_sc(nat)%equilibrium_pos+shift0)  ! =rsc(k)-rsc(1)+atom0(isave)
!    rn1=atompos(:,nat)+v2a(shift0)  ! =rsc(k)-rsc(1)+atom0(isave)
!    call get_n_tau_r_mod(rn1,ni,taui)
     taui=iatomcell0(nat)
     ni  =iatomcell(:,nat)

!    call findatom_cart(rn1,iat)  ! iat is the corresponding atompos index
!    if (iat.le.0) then
!       write(ulog,*)'CORRESP: Supercell atom#',nat,' could not be indentified with any of the atompos'
!       write(ulog,*)'Increase the number of shells and rerun, otherwise the range of FCs is'
!       write(ulog,*)'determined by the values of nshells(2)=',nshells(2,:)
!       write(*   ,*)'CORRESP: Supercell atom#',nat,' could not be indentified with any of the atompos'
!       write(*   ,*)'Increase the number of shells and rerun, otherwise the range of FCs is'
!       write(*   ,*)'determined by the values of nshells(2)=',nshells(2,:)
!       write(*   ,*)'atom_sc(',nat,',)i+shift0=',matmul(cart_to_prim,rn1)
!    endif
!    if(nat.eq.1 .and. iat.ne.isave) then
!       write(*,*)'WRITE_CORRESPONDANCE: atom0 for iat=',iat,' and isave=',isave,'should be the same!'
!       stop
!    endif
!    taui=iatomcell0(iat)
!    ni  =iatomcell(:,iat)  ! this is the identity of atom_sc(nat)
!    if(taui.eq.0) then
!        write(*,*)'WRITE_CORRESPONDANCE: cannot find n_tau for atom#',nat,' in the primitive cell!'
!        stop
!    endif
!    write(*,*)'WRITE_CORR: before entering findatom_sc, taui,ni=',taui,ni

!    do i0=1,natom_prim_cell
!       rsc=v2a(atom_sc(1)%equilibrium_pos+atom0(i0)%equilibrium_pos-shift0) ! this is r_sc(k)=r_sc(1)+ri0-risave

!       call get_n_tau_r_mod(v2a(atom_sc(nat)%equilibrium_pos+shift0),ni,taui)
! this won't work if taui,ni are not defined
!    call findatom_sc(ni,taui,iat_sc)  ! now iat_sc is the atom index in the supercell
!    if (iatom.ne.nat) then  ! Atom i OF PRIMCELL DOES NOT NEED TO MATCH ATOM i OF SUPERCELL
!       write(ulog,*)'WRITE_CORRESPONDANCE: iatom_sc .ne. nat or taui:',iatom,nat
!       write(ulog,*)'this is not an error but it is unusual!'
!    endif
!    nsi  = nint(matmul(r0g,dfloat(ni)))  ! should be zero!

        do al=1,3             ! this is how we order each line

        do g=1,map(rnk)%ngr  ! sum over groups
!if (keep_fc2i(g).ne.1) then
!   write(ulog,*)' CORRESP: group#',g,' was not included because it went outside of the WS cell'
!   cycle  ! go to the next group
!endif
           do t=1,map(rnk)%nt(g) ! sum over all terms in that group
              if ( taui.eq. map(rnk)%gr(g)%iat(1,t) .and.  &
           &       al  .eq. map(rnk)%gr(g)%ixyz(1,t) ) then
! map()%gr% iat etc is the atom index from atompos, so we use iatomcell stuff to find (tau,n)
                 be= map(rnk)%gr(g)%ixyz(2,t)
                 j = map(rnk)%gr(g)%iat(2,t)
                 tauj = iatomcell0(j)
!  Below, ni refers to the translations of the supercell'
                 nj(:)= iatomcell(:,j) + ni(:) ! translate by ni
!                nsj  = nint(matmul(r0g,dfloat(nj)))
! Identify neighbor j within the SCell,
                 call findatom_sc(nj,tauj,j_sc)
                 if (j_sc.eq.0) then
                     write(ulog,4)'WRITE_CORRESPONDANCE:jatom_sc not found for j,tauj,beta,nj=',j,tauj,be,nj
                     write(ulog,4)'for atom0(i),ixyz,group,term ',nat,al,g,t
                     stop
                 endif
! if j and j+R(supercell) are in the same group but with opposite ampterms
! then their sum is cancelled in the force as they both have the same
! displacement. This will lead to undetermination in evaluation of the FCs if all
! terms in the group cancel in this way, and the
! corresponding group in FC2 will be evaluated to be zero. That will also
! affect the ASR and produce a violation of the sum rules, unless rotational invariances are imposed
                 rij = length(atom0(taui)%equilibrium_pos - atompos(:,j)) !  pos(:,jatom))
             
 if(al.eq.1 .and. be.eq.1)  write(ucor,6)nat,taui,ni,al,j,tauj,nj,be,j_sc,t,g,rij
 if(al.eq.1 .and. be.eq.1)  write(ulog,4)' ',nat,taui,j,tauj,j_sc,t,g,rij
              endif
           enddo
        enddo
        enddo
        write(ucor,*)'------------------------------'
        write(ulog,*)'------------------------------'
     enddo

  write(ulog,*)' Correspondance of harmonic FCs and pairs of supercell atoms established'
  write(ulog,*)' if chosen range in the model is too long-ranged, it is possible that more'
  write(ulog,*)' than one term in the fc list (group,term) are associated with the same supercell pair'
  write(ulog,*)' And if too small, some supercell pairs will not have any FCs associated with them'
! this can be taken care of by setting to zero the extra terms of longer distance in keep_fc2i

4 format(a,7(1x,i6),3x,f10.4)
5 format(2(i4,1x,' [ ',i2,' (',i2,',',i2,',',i2,') ] ; ',1x),i6,1x,i4,2x,f7.3,2x,f7.3)
6 format(2(i4,1x,' [ ',i2,' (',i2,',',i2,',',i2,') ] ',i1,1x,' ; '),i6,1x,i6,1x,i4,2x,f7.3,2x,f7.3)

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
 subroutine warn2(unt,mesg)
 implicit none
 integer, intent(in) :: unt
 character(len=*), intent(in) :: mesg

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
 write(unt,*)'|                                                                  |'
 write(unt,*)'********************************************************************'
 write(unt,*)'|   '//mesg
 end subroutine warn2
!============================================================
 subroutine write_lat_fc(ntindep,ntrms,fn)
 use svd_stuff
 use ios
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
 integer, intent(in) :: ntindep(maxrank),ntrms(maxrank)
 character(len=*), intent(in):: fn
 integer i,j,tau,nat,nm(3),largest
 real(r15) ri,rr(3)

! this is to reduce the size of lat_fc.dat
 if(fn .eq. 'lat_fc.dat') then
    call find_largest_atom_number(largest)
    if(largest.lt.1000) largest=min(6500,natoms)
 else 
    largest = natoms
 endif

 open(ufco,file=fn)
 write(ufco,*)' Crystal data: translation vectors of the primitive cell '
 write(ufco,9)r01
 write(ufco,9)r02
 write(ufco,9)r03
 write(ufco,*)' Crystal data: atoms in primitive cell: label,type,x,y,z,mass '
 write(ufco,*)natom_prim_cell
 do i=1,natom_prim_cell
   write(ufco,6)i,atom0(i)%name, atom0(i)%at_type,atom0(i)%equilibrium_pos,atom0(i)%mass
 enddo
 write(ufco,*)' Crystal data: groups ',map(1:8)%ngr,' indepterms=',map(1:8)%ntotind,' all terms=',map(1:8)%ntot
 write(ufco,*)' Included ranks of FCs '
 write(ufco,*)  include_fc(:)
 write(ufco,*)' Number of FCs for each rank '
 write(ufco,*)  ntrms(:)
 write(ufco,*)' Number of independent FCs for each rank '
 write(ufco,*)  ntindep(:)
 write(ufco,*)' maxshells,Neighborshell atoms: i,x,y,z,type_tau,n1,n2,n3 '
 write(ufco,*)  maxshells,largest,natoms
! write up to the used shells
 nat=0
 do i=1,largest !natoms
   ri = length(atompos(:,i)) !-atompos(:,1))
      nat=nat+1
      write(ufco,7)i,(atompos(j,i),j=1,3), iatomcell0(i),(iatomcell(j,i),j=1,3),ri
      nm = iatomcell(:,i)
      tau= iatomcell0(i)
      rr=nm(1)*r01+nm(2)*r02+nm(3)*r03
      if (length(rr+atom0(tau)%equilibrium_pos-atompos(:,i)).gt.1d-4) then
             write(*,*)' WRITE_LAT_FC error: n,tau does not correspond to atompos '
             write(*,8)' taui   = ',atom0(tau)%equilibrium_pos
             write(*,8)' ni*r0i = ',rr
             write(*,8)' atompos= ',atompos(:,i)
             write(*,5)' n,tau,j= ',nm,tau,i
          stop
      endif
 enddo
 close(ufco)

 open(173,file='primlatt.xyz')
 write(173,*) natom_prim_cell
 write(173,*)'# ', ntindep(:)
 do i=1,natom_prim_cell
   write(173,8)atom0(i)%name,atom0(i)%equilibrium_pos
 enddo
 close(173)

 open(173,file=fn//'.xyz')
 write(173,*) nat !oms
 write(173,*)'# ', ntindep(:)
 do i=1,natoms
!  if (length(atompos(:,i)).le.lmax) then  ! write up to lmax
       write(173,8)atom0(iatomcell0(i))%name,(atompos(j,i),j=1,3)
!  endif
 enddo
 close(173)

5 format(a,2x,3(i3),3x,i2,3x,i5,9(1x,f11.5))
6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(a,6(2x,f19.10))
9 format(9(2x,f19.10))
 end subroutine write_lat_fc
!============================================================
 subroutine read_latfc_allocate
 use svd_stuff
 use ios
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
! integer, intent(in) :: ntindep(maxrank),ntrms(maxrank)
! character(len=*), intent(in):: fn
 integer i,j,l,tau,nat,nm(3),largest,uio
 real(r15) ri,rr(3)
 character(LEN=99) line
 character(LEN=22) wrd1
 character(LEN=12) wrd2
 character(LEN=11) wrd3


 uio=987
 open(uio,file='lat_fc.dat',status='old')

 write(uio,*)line
 write(uio,*)line
 write(uio,*)line
 write(uio,*)line
 write(uio,*)line
 write(uio,*)line
 do i=1,natom_prim_cell
   write(uio,*)line
 enddo
 read(uio,4)wrd1,map(1:8)%ngr,wrd2,map(1:8)%ntotind,wrd3,map(1:8)%ntot

4 format(a22,4i4,a12,4i4,a11,4i5)
 read(uio,*)line
 read(uio,*) include_fc !(1:8)
 read(uio,*)line
 read(uio,*)line !  ntrms(:)
 read(uio,*)line
 read(uio,*)line !  ntindep(:)
 read(uio,*)line
 read(uio,*)  maxshells,largest,natoms
! write up to the used shells
 if(largest.ne.natoms)write(ulog,*)"WARNING ERROR in reading lat_fc.dat: largest,natoms='largest, natoms"
 nat=0
 do l=1,largest !natoms
   ri = length(atompos(:,i)) !-atompos(:,1))
      nat=nat+1
      read(uio,*)i,(atompos(j,i),j=1,3), iatomcell0(i),(iatomcell(j,i),j=1,3),ri
      nm = iatomcell(:,i)
      tau= iatomcell0(i)
      rr=nm(1)*r01+nm(2)*r02+nm(3)*r03
      if (length(rr+atom0(tau)%equilibrium_pos-atompos(:,i)).gt.1d-4) then
             write(*,*)' WRITE_LAT_FC error: n,tau does not correspond to atompos '
             write(*,8)' taui   = ',atom0(tau)%equilibrium_pos
             write(*,8)' ni*r0i = ',rr
             write(*,8)' atompos= ',atompos(:,i)
             write(*,5)' n,tau,j= ',nm,tau,i
          stop
      endif
 enddo
 close(uio)


5 format(a,2x,3(i3),3x,i2,3x,i5,9(1x,f11.5))
6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(a,6(2x,f19.10))
9 format(9(2x,f19.10))
 end subroutine read_latfc_allocate
!============================================================
 subroutine write_atompos
 use svd_stuff
 use ios
 use atoms_force_constants
 use params
 use lattice
 use constants
 implicit none
 integer i,iunit,tau,n(3)

 iunit=123

 open(iunit,file='atompos.xyz')
 write(iunit,*)natoms
 write(iunit,*)'# name, x , y , z , i,tau,n(3), reduced coordinates '
 do i=1,natoms
   tau = iatomcell0(i)
   n   = iatomcell(:,i)
   write(iunit,8)atom0(tau)%name,atompos(:,i) , i,tau,n,matmul(cart_to_prim,atompos(:,i))
 ! write(*    ,8)atom0(tau)%name,atompos(:,i) , i,tau,n,matmul(cart_to_prim,atompos(:,i))
 enddo
 close(iunit)

8 format(a,3(1x,f12.5),i5,i3,'(',3i2,')',3(1x,f6.3))
 end subroutine write_atompos
!============================================================
 subroutine write_supercell(j)
! use svd_stuff
! use ios
 use atoms_force_constants, only: atom_sc, natom_super_cell
! use params
 use geometry, only : v2a
 use lattice, only: cart_to_prim
 use constants, only : pi
 implicit none
 integer, intent(in) :: j
 integer i,iunit
 character(len=1) fdfil
 iunit=123
 write(fdfil,'(i1)')j
 open(iunit,file='supercell_'//fdfil//'.xyz')
 write(iunit,*)'# ', natom_super_cell
 write(iunit,*)'# name,  x(i)  ,  y(i)  ,  z(i)  , i ,tau,n(3), atompos(i), reduced coordinates '
 do i=1,natom_super_cell
    write(iunit,7)atom_sc(i)%name,atom_sc(i)%equilibrium_pos,i,atom_sc(i)%cell%tau, &
    atom_sc(i)%cell%n,atom_sc(i)%cell%atomposindx,matmul(cart_to_prim,v2a(atom_sc(i)%equilibrium_pos))
 enddo
 close(iunit)

7 format(a,3(1x,f12.5),i5,i3,'(',i2,2(',',i2),')',i5,3(1x,f9.5))

 end subroutine write_supercell
!============================================================
 subroutine write_invariance_violations(ulog,n,fcs)
 use svd_stuff, only : atransl,arot,brot,ahuang, amat,bmat,transl_constraints, &
 &  rot_constraints,huang_constraints,inv_constraints,dim_al,itrans,irot,ihuang
 use ios, only : umatrx
 use constants, only : r15
 implicit none
 integer, intent(in) :: ulog,n
 real(r15), intent(in) :: fcs(n)
 integer i
 real(r15) prod,errmax,junk,err,num,denom


   write(ulog,3)'Dimension of extracted FCS is=',n
   write(*   ,3)'Dimension of extracted FCS is=',n
   write(ulog,*)' dimensions of ahuang=',size(ahuang,1),'x',size(ahuang,2)
   write(*   ,*)' dimensions of ahuang=',size(ahuang,1),'x',size(ahuang,2)
   write(*,   *)' # of huang constraints is=',huang_constraints
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
   write(ulog,*)' if not=15, it means there were some that were exactly zero and removed'
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
   write(umatrx,2)'Max and average errors in force-displacements,percent deviation=', &
   &    errmax,err,sqrt(num/denom)*100
   write(ulog,2)'Max and average errors in force-displacements,percent deviation=', &
   &    errmax,err,sqrt(num/denom)*100

2 format(a,8(2x,g12.5))
3 format(a,i5,3(2x,g12.5))
 end subroutine write_invariance_violations
!============================================================
 subroutine find_largest_atom_number(largest)
!! find the largest atom number in the list map(2)%gr(g)%iat(2,t)
 use svd_stuff
 use ios, only : ulog
 implicit none
 integer, intent(out) :: largest
 integer i,g,t,j

    largest=0
    gloop: do g=1,map(2)%ngr
!      if(keep_fc2i(g).ne.1) cycle gloop
       tloop: do t=1,map(2)%nt(g)
          j  = map(2)%gr(g)%iat(2,t)
          if(j.gt. largest) largest=j
       enddo tloop
    enddo gloop
    write(*   ,*)' The largest atom number found in atompos was ',largest
    write(ulog,*)' The largest atom number found in atompos was ',largest

 end subroutine find_largest_atom_number
