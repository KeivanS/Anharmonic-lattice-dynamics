program FC234
!
! todo: supress redundant parameters lines maxshell. The latter
! should be based on the number of requested shells by fc2
! Same for error tolerances
! add option of reading existing fc files and fitting the remaining ones.
!
! if a line of amat is zero, then check whether the corresponding line
! in bmat is zero or not. if it is remove the line from amat, else
! print in log file the inconsistency and remove it anyways.
! CHECK THE iatomxyz PART OF FC4, THERE ARE SOME large NUMBERS IN IT!!
! also compare the si data with other papers (quteish...)
! allow some of the fcs of a given rank to be fixed and read from the
! file if include_fc=2, and the rest to be fitted and appended to the
! file for further use.
!
!------------------------------------------------
! program to extract force constants from ab initio force-displacement data
! eventually with (linear) constraints of symmetry
! INPUT FILES : param.in, POSCAR, OUTCAR
! in params.inp file user specifies the specifications of the primitive cell
! what potential parameters will be kept how many shells are included,
! whether or not transl and rot inv constraints are imposed
! the tolerance for equating two coordinates
! specifications(type mass name) of the atoms in prim cell and their red coord
! In POSCAR file the same specifications but for the super cell are read.
! Finally OUTCAR constains the force-displacements data
! OUTPUT FILES: log.dat contains log of the runs and intermediate outputs
! fcs.dat : contains the obtained force constants
! amatrx.dat : is the matrix A and column b to be solved by SVD
! svd-results.dat: contains the output of SVD algorithm, errors, array w
! maps.dat: contains the map taking the oned array of terms to atom i,xyz
! by K. Esfarjani December 2006
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
! show which FCs are included in a given group
!
! this new version separates the invariance constraints from the force-displacement
! data and enforces the former "exactly" using a separate svd algorithm.
! we need to first subtract from the forces the long-range Coulombic part in order
! to insure the remaining part is short-ranged.
! needs to read params.born to get dielectric constant and Born charges; the long-range
! part of the FCs is (Zi.q)(Zj.q)/eps q^2 corresponding to force: Fi=- sum_q (Zi.q)(Zj.q)/eps q^2 e^iq.Rij u_j
! need to calculate this force, subtract it from the DFT input force, and fit the remainder.
!
use lattice
use params
use ios
use atoms_force_constants
use force_constants_module
use svd_stuff
use geometry
implicit none
! type(vector) v
integer i,j,n3,rank,iunit,tra,rot,hua,eliminate_fc,fcalloc, nz_index,l ! g, t
character xt*1,fn*11,poscar*7,outcar*7
logical ex
real(8), allocatable:: mat(:,:),mat_inverse(:,:)
real(8), allocatable :: newamat(:,:),newbmat(:),afrc(:,:),bfrc(:),aux1(:,:),rm_zero_energy(:),fc(:),wts(:),qmat(:)
integer, allocatable :: frc_constr(:),auxi(:,:)
real(8) error,ermax,sd(4),sig,small,suml,V1, kB, Temp, kBT, norm_wts !,a1,a2,a3
real tim
 character now*10,today*8,zone*5

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(tim)

 open(ulog  ,file='log.dat'   ,status='unknown')
 open(utimes,file='times.dat' ,status='unknown')
 open(umap  ,file='maps.dat'  ,status='unknown')
 open(umatrx,file='amatrx.dat',status='unknown')
 open(ucor  ,file='corresp.dat',status='unknown')
 open(ufco  ,file='lat_fc.dat',status='unknown')
 open(ufit1 ,file='fc1_fit.dat',status='unknown')
 open(ufit2 ,file='fc2_fit.dat',status='unknown')
 open(ufit3 ,file='fc3_fit.dat',status='unknown')
 open(ufit4 ,file='fc4_fit.dat',status='unknown')

 write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' Program FC234 was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(ulog,'(a,f10.4)')' STARTING TIME OF THE PROGRAM                   IS ',tim
 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' This program was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(utimes,'(a,f10.4)')' At the start of the Program FC234 the  TIME IS ',tim

  call set_maxterms   ! no need to read them anymore, they are incremented if necessary

! read from params.inp the atoms in prim cell in reduced units, assign their
! mass, type, tau; also read specs of the prim cell, and max no of neighborshells
  call read_input
! outputs: lattparams,primlat,nshells,include_fc,itrans,irot,natoms0,mas,atname,
! atompos0 (reduced coords within the primitive cell)

  maxshells = maxneighbors

  call write_out(ulog,'latticeparameters',latticeparameters)
  call write_out(ulog,'primitive lattice',3,3,primitivelattice)
  write(ulog,*)'atompos0_d'
  do i=1,natoms0
!    call write_out(ulog,'atom0 ',3,atompos0(:,i))
     write(ulog,8)i,atompos0(:,i)
  enddo
  call write_out(ulog,'# of different atom types ',natom_type)
  call write_out(ulog,'# of each type in primcel ',natom)
  call write_out(ulog,'# mass of different types ',mas)

! get from lattparam and atompos0 previously read the cartesian coordinates
! of the atoms in the primitive cell(+neighbors) stored in atompos;
! also find symmetry of the crystal and its matrices and generates neighbors

 call cpu_time(tim)
 write(utimes,'(a,f10.4)')' READ_INPUT CALLED, TIME                       IS ',tim

! This does not work: fcinit must be called once or latparam(4:6) will be overwritten
  maxatoms=20; imaxat=1
  do while (imaxat.ne.0)
     maxatoms=maxatoms+300
     write(6,*)' maxatoms=',maxatoms
     write(ulog,*)' maxatoms=',maxatoms
     if (allocated(iatomcell0)) deallocate(iatomcell0)
     if (allocated(iatomcell))  deallocate(iatomcell)
     if (allocated(atompos))    deallocate(atompos)
     allocate(atompos(3,maxatoms), iatomcell(3,maxatoms),iatomcell0(maxatoms))
! inputs: natoms0,latticeparameters,primitivelattice,atom_type,atompos0
! outputs: r0i, atompos (cart coord of atoms and natoms neighbor shells),
! iatomcell, iatomcell0, iatomneighbor, iatomop, ...
     call force_constants_init(latticeparameters,primitivelattice,natoms0,  &
      &     atom_type,atompos0)  !(n,natom_type,atompos0)
  enddo

! reallocate atompos,iatomcell and iatomcell0
  allocate(auxi(3,natoms),aux1(3,natoms))

      auxi(1:3,1:natoms)=iatomcell(1:3,1:natoms) 
      aux1(1:3,1:natoms)=atompos(1:3,1:natoms)
      deallocate(atompos,iatomcell)
      allocate(atompos(3,natoms),iatomcell(3,natoms))
      atompos=aux1; iatomcell=auxi

      auxi(1,1:natoms)=iatomcell0(1:natoms) 
      deallocate(iatomcell0)
      allocate(iatomcell0(natoms))
      iatomcell0=auxi(1,:)             
      maxatoms=natoms

  deallocate(aux1,auxi)
       
 call cpu_time(tim)
 write(utimes,'(a,f10.4)')' FC_INIT CALLED, TIME                          IS ',tim
 write(utimes,*)' FC_init called, atmpos, neighbs, iatmcell calculated'

! r01,r02,r03 translation vectors of the primitive lattice were just defined
  call write_out(ulog,'r01=',r01)
  call write_out(ulog,'r02=',r02)
  call write_out(ulog,'r03=',r03)

  call write_out(ulog,' number of neighboring atoms, natoms',natoms)
  write(ulog,*)' i, atompos(i), tau(i), n(i) '
  do i=1,natoms
     write(ulog,8) i,atompos(:,i),iatomcell0(i),iatomcell(:,i)
  enddo

5 format(9(1x,f11.6))
8 format(i4,3(2x,f10.5),3x,i2,2x,'(',3(1x,i2),1x,')',2x,i4)
9 format(i5,9(1x,f9.4))

  call make_reciprocal_lattice(r01,r02,r03,g01,g02,g03) ! beware of the factor of 2pi in g0

  call set_neighbor_list
! inputs: atompos,
! outpus: atom0%equilb_pos,shells%no_of_neighbors , rij, neighbors(:)%tau and n

  call write_neighbors
!write(*,*) "The value of map(2)%ngr is: ", map(2)%ngr
! define the dimensions of the FC arrays, for each rank if include_fc.ne.0
! collects identical force_constants

  call setup_maps
! inputs: outputs of force_constants_init: iatomop etc...
! outputs: map structure including for rank, the groups of indep FCs and the full FCs
! maxterms, nshells, ngroups, nterms, inv_constraints, ngr
  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after setup_maps and before write_correspondance IS ',tim

!--------------------------------------------------------------------
! eventually read FCs from a file if already present
  do rank=1,4
     write(xt,'(i1)')rank
     fn = 'fc'//xt//'.dat'
     iunit = ufco+rank
     if (include_fc(rank) .eq. 2 ) then   ! this rank is not to be included
      fcalloc=0
!      write(*,*) "The value of RFCS2 allocation is: ", fcalloc
    ! in the fit, read the fcs from a file
      fn='fc'//xt//'_fit.dat'
        inquire(file="fc2.dat",exist=ex)
        if (ex) then
           !open(iunit ,file=fn ,status='old')
      !     call read_harm_fcs(2)
        else
           write(*,*)'For this rank there is no FC data file present ',rank
           write(*,*)'check your inputs and run the program again'
           write(ulog,*)'For this rank there is no FC data file present ',rank
           write(ulog,*)'check your inputs and run the program again'
           stop
        endif
     elseif (include_fc(rank) .eq. 1 ) then
        open(iunit ,file=fn ,status='unknown')
     endif
  enddo

!--------------------------------------------------------------------
! depending on the constraints, setup the lines of the A and B matrices:
  call set_translational_inv_constraints
  write(ulog,*)'Number of translational invce constraints is=',transl_constraints

  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after set_translational_inv_constraints IS ',tim

  call set_rotational_inv_constraints
  write(ulog,*)'Number of rotational invarnce constraints is=',rot_constraints

  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after set_rotational_inv_constraints    IS ',tim

  call set_huang_inv_constraints
  write(ulog,*)'Number of Huang    invariance constraints is=',huang_constraints

  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after set_Huang_inv_constraints         IS ',tim

 call write_lat_fc(map(:)%ntotind,map(:)%ntot) !ngroups,nterms)

!!  call homogeneous_constraints_overlap(n_constr)
!!
!!! do the SVD of the overlap
!! allocate( us(n_constr,ngr),vs(ngr,ngr),ws(ngr) )
!! us = ahom
!! call svdcmp(us,n_constr,ngr,n_constr,ngr,ws,vs)

! allocate( us(n_constr,n_constr),vs(n_constr,n_constr),ws(n_constr) )
! us = overl
! call svdcmp(us,n_constr,n_constr,n_constr,n_constr,ws,vs)

! truncate and keep the highest eigenvalues of ws
! call truncate(ngr,ws,n_kept)
! write(ulog,*)' the largest  ws is=',ws(1)
! write(ulog,*)' the smallest ws kept and its index is=',n_kept,ws(n_kept)
! write(ulog,*)' their ratio is=',ws(n_kept)/ws(1)

! for the kept eigenvalues calculate the correction or the residual FCS
! by projecting the initial solution onto the kernel of the constraint matrix
! here we are constructing the projection matrix which should later multiply X0
! allocate( project(n_constr,n_constr) )
! project=0
! do i=1,n_constr
!    project(i,i)=1
! do j=1,n_constr
!    do n3=1,n_kept !+1,n_constr
!       project(i,j)=project(i,j)-vs(i,n3)*vs(j,n3)
!    enddo
! enddo
! enddo

! we need to multiply the final svd solution vector by project, to remain in the kernel of C

!--------------------------------------------------------------------
! loop over different FP runs with different POSCARi-OUTCARi files
! supercell information enters here.
! this 1st loop is to calculate force_constraints in order to allocate aforce
 n3 = 0
 allocate(frc_constr(fdfiles))
 force_constraints=0
 do i=1,fdfiles

     write(xt,'(i1)')i
     poscar='POSCAR'//xt
     outcar='OUTCAR'//xt

! read the ideal crystal data (atoms in supercell at equilbrm) from POSCAR
     call read_crystal(poscar)
! inputs:
! outputs: supercell's lattice_parameter,r1i,g1i,boxxyz,om,natom_super_cell,
! atom_sc%equilibrium_pos
    write(ulog,*) 'supercell crystal read'

! make the matrix a_ij=r0i.dot.gj
     call make_r0g
    write(ulog,*) 'make r0g called'

! read the position-force data for a set of configurations
     call read_force_position_data(outcar,frc_constr(i))
    write(ulog,*) 'force - position read: i, frc_constr(i)= ',i,frc_constr(i)
! outputs: nconfigs,displ(=actual positions),force,frc_constr(i)

     force_constraints=force_constraints+frc_constr(i)
 enddo

 write(ulog,*)'Now allocating aforce of dims=',force_constraints,ngr
 allocate( aforce(force_constraints,ngr),bforce(force_constraints) )

!--------------------------------------------------------------------
! now that force_constraints is known, read and setup the aforce matrix
 structures: do i=1,fdfiles

     write(xt,'(i1)')i
     poscar='POSCAR'//xt
     outcar='OUTCAR'//xt

! read the ideal crystal data (atoms in supercell at equilbrm) from POSCAR
     call read_crystal(poscar)
! inputs:
! outputs: supercell's lattice_parameter,r1i,g1i,boxxyz,om,natom_super_cell,
! atom_sc%equilibrium_pos
    write(ulog,*) 'supercell crystal read'

! make the matrix a_ij=r0i.dot.gj
     call make_r0g
    write(ulog,*) 'make r0g called'

! read the position-force data for a set of configurations
     call read_force_position_data(outcar,frc_constr(i))
    write(ulog,*) 'force - position read for the second time '
! outputs: nconfigs,displ(=actual positions),force,frc_constr(i)

     call calculate_and_write_displacements
! inputs: displ, atom_sc%equilibrium_pos (get displ from actual coordinates)
! outputs: displ

     call identify_atoms_in_supercell
! inputs: atom0%eq_pos and atom_sc%eq_pos (finds the shift overlap atom0 with atom_sc)
! outputs: shift0, atom_sc%cell%n and tau, njmin,njmax

     call write_correspondance2

  call cpu_time(tim)
  write(utimes,'(a,i3,a,f10.4)')' STRUCTURE ',i,' after write_correspondance time IS ',tim

     write(ulog,*)' MAIN: total # of independent FCs of all rank, ngr=',ngr

     allocate(afrc(frc_constr(i),ngr),bfrc(frc_constr(i)))
     write(ulog,*)' calling set_force_displacement_matrix, for str=',i

     call set_force_displacement_matrix(frc_constr(i),afrc,bfrc)

     write(ulog,*)' set_force_displacement_matrix called, aforce being filled now'
     write(ulog,*)' n3,n3+frc_constr(i)=',n3,n3+frc_constr(i)
     aforce(n3+1:n3+frc_constr(i),1:ngr) = afrc(1:frc_constr(i),1:ngr)
     bforce(n3+1:n3+frc_constr(i)      ) = bfrc(1:frc_constr(i)      )
     n3 = n3 + frc_constr(i)
     deallocate(afrc,bfrc)

! writing the coordinates of primitive cells in the supercell
! reference of each prim cell is atom #1
     open(111,file='struct-supercell-'//xt)
     write(111,5)r1
     write(111,5)r2
     write(111,5)r3
     do j=1, natom_super_cell
        if(atom_sc(j)%cell%tau .ne. 1) cycle
        write(111,9)j,atom_sc(j)%cell%n(1)*r01+atom_sc(j)%cell%n(2)*r02  &
&                    +atom_sc(j)%cell%n(3)*r03,atom_sc(j)%equilibrium_pos
     enddo
     close(111)

 enddo structures

  if (n3.ne.force_constraints) then
     write(ulog,*)'n3, force_constraints are not equal ',n3,force_constraints
     write(ulog,*)'the program will stop'
     write(*   ,*)'n3, force_constraints are not equal ',n3,force_constraints
     write(*   ,*)'the program will stop'
     stop
  endif

  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after structure loop                    IS ',tim


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
!write(*,*) "CALL BEFORE INCLUDE_CONSTRAINTS"
  call include_constraints
! inputs: inv_constraints,force_constraints,
! outputs:  dim_al,dim_ac, amat, bmat
  !write(*,*) "CALL AFTER INCLUDE_CONSTRAINTS"
! we remove the following because if one starts with reference positons in POSCAR
! which are not equilibrium position, then FC1 are not zero and lines with u=0
! should still be kept, and not thrown away.
! the following removes lines=0 from amat (where all displacements are zero)

! SINCE I HAVE TO CREATE EXPLICIT INTERFACE FOR SUBROUTINE REMOVE_ZEROS IN ORDER FOR ME TO ACCESS
! AND REMOVE ENERGY WHERE THE LINES BMAT IS ZERO BUT I DONT WANT TO DO THAT BECAUSE I AM GOING TO ALTER
! THE FORMULATION OF REMOVE_ZERO SUBROUTINE. ONE OF THE WAY I SEE IS CHECK FOR ZERO COLUMN DIRECTLY HERE.

small=1d-7
nz_index=0
do i=1,dim_al    !shape(amat,1)
   suml = maxval(abs(amat(i,:)))
   if (suml .gt. small) then
      nz_index=nz_index+1
   endif
enddo
allocate(rm_zero_energy(nz_index))
nz_index=1
do i=1,dim_al
   suml=maxval(abs(amat(i,:)))
      if (suml .gt. small) then
         rm_zero_energy(nz_index)=energy(i)
         nz_index=nz_index+1      
      endif
enddo
allocate(mat(dim_ac,dim_ac),mat_inverse(dim_ac,dim_ac))
allocate(wts(dim_al),qmat(dim_ac))
V1=1.0d0
kB=0.00008617333262145d0
Temp=116040000.518120d0
kBT=kB*Temp
norm_wts=0.0d0
do i=1,dim_al
   wts(i)=exp(-1.0d0*(abs(energy(i))-V1)/kBT)
   write(*,*) "NORMALIZED LOOP WTS: ", wts(i)
   norm_wts=norm_wts+wts(i)
enddo
write(*,*) "THE NORMALIZED VALUE FOR WTS IS: ", norm_wts
do i=1,dim_al
   wts(i)=wts(i)/norm_wts
enddo
!wts=1.0d0
write(*,*) "THE VALUE OF WTS is: ", wts
qmat=0.0d0
do i=1,dim_ac
   do l=1,dim_al
      qmat(i)=qmat(i)+wts(l)*bmat(l)*amat(l,i)
   enddo
   do j=1,dim_ac
      mat(i,j)=0.0d0
      do l=1,dim_al
         mat(i,j)=mat(i,j)+wts(l)*amat(l,i)*amat(l,j)
      enddo
   enddo
enddo
mat_inverse=0.0d0
do i=1,dim_ac
   mat_inverse(i,i)=1.0d0
enddo
!call inverse_real(mat,mat_inverse,dim_ac)
!fc=matmul(mat_inverse,qmat)
!do i=1,dim_ac
!write(*,*) "THE VALUE OF FC_WEIGHTED IS: ", fc(i)
!enddo
!  deallocate(wts,qmat,mat,mat_inverse)

call remove_zeros(dim_al,dim_ac, amat, bmat,newdim_al)
!if ( include_fc(2) .ne. 2) then
 ! write(*,*) "THE NEWDIM_AL AND DIM_AC is: ",newdim_al, dim_ac
allocate(newamat(newdim_al,dim_ac),newbmat(newdim_al))
!endif
!if (include_fc(2) .eq. 2) then
!  allocate(newamat(newdim_al,map(2)%ngr:dim_ac),newbmat(map(2)%ngr:dim_al))
!   write(*,*) "The value for newdiv_al, map(3).ntotind+map(4).ntotind is: ", map(3)%ntotind+map(4)%ntotind
!   allocate(newamat(newdim_al,map(3)%ntotind+map(4)%ntotind),newbmat(map(2)%ngr:dim_al)) ! what is the dimension here map(2)%ngr:dim_al
   ! or 1:dim_al-map(2)%ngr`  
!endif
  newamat(:,:) = amat(1:newdim_al,1:dim_ac)
  newbmat(:  ) = bmat(1:newdim_al)
  deallocate(amat,bmat)
  allocate(amat(newdim_al,dim_ac),bmat(newdim_al))
  amat = newamat
  bmat = newbmat
  dim_al = newdim_al

  deallocate(newamat,newbmat)

! give warnings if any column in aforce is totally zero
  call check_zero_column(dim_al,dim_ac,amat)
!write(*,*) "CHECK ZERO COL, DIM_AC", dim_ac
 tra=0; rot=0; hua=0
 if (itrans .ne. 0) tra= transl_constraints
 if (irot   .ne. 0) rot= rot_constraints
 if (ihuang .ne. 0) hua= huang_constraints
 write(ulog,*)' Number of invariance constraints   =',inv_constraints
 write(ulog,*)' Size of the homogeneous part of A  =',tra+rot+hua
 if (inv_constraints.ne.tra+rot+hua) then
     inv_constraints = tra+rot+hua
     write(ulog,*)' Setting # of invariance constraints=',inv_constraints
 endif
 write(ulog,*)' Size of the unknown FCs of all rank=',dim_ac

 if (dim_ac .lt. (tra+rot+hua)) then
    write(ulog,*)' the # of indept FCs is too small for doing any elimination'
    eliminate_fc=0
 endif
 if (dim_ac .lt. (tra+rot+hua)+1) then
    call warn(ulog)
    write(ulog,*)' the # of constraints is too large and they are not consistent'
    write(ulog,*)' you need to release some of the constraints or increase the range'
    call warn(6)
    write(6,*)' the # of constraints is too large and they are not consistent'
    write(6,*)' you need to release some of the constraints or increase the range'
 !  stop
 endif
! we remove the following because if one starts with reference positons in POSCAR
! which are not equilibrium position, then FC1 are not zero and lines with u=0
! for now do it without reduction, put the full A matrix
! then do the reduction and check validity
! then impose other invariance relations
! then put higher order terms
!----------------------- END FOR TESTING -----------------------

  allocate(fcs(dim_ac),sigma(dim_ac))
  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after include_constraints  before SVD   IS   ',tim

! write the a and b matrices
  if (verbose) then
     write(umatrx,*)' before call to svd, amat and bmat are:'
 !    write(*,*) "DIM_AL, DIM_AC is: ", dim_al, dim_ac
     do i=1,dim_al
        write(umatrx,17)(amat(i,j),j=1,dim_ac),bmat(i) ! what is this and where do you get amat here?
     enddo
  endif
17 format(999(1x,g9.2))

! do the SVD decomposition
!  HERE I TRY TO ADD THE WEIGHTED VALUE BASED ON TEMPERATURE IN THE BELOW LINES BEFORE CALLING TO SVD SET
! SO THE AMAT WILL BE AT*A AND BMAT WILL BE AT*B
!   allocate(mat(dim_ac,dim_ac),mat_inverse(dim_ac,dim_ac))
!   allocate(wts(dim_al),qmat(dim_ac))
!   write(*,*) "THE VALUE OF RM_ZERO_ENERGY is: ",rm_zero_energy
!   V1=1.0d0
!   kB=0.00008617333262145d0
!   Temp=116040000.518120d0
!   kBT=kB*Temp
!   norm_wts=0.0d0
!   do i=1,dim_al
!      wts(i)=exp(-1.0d0*(abs(rm_zero_energy(i))-V1)/kBT)
!      norm_wts=norm_wts+wts(i)
!   enddo
!   do i=1,dim_al
!      wts(i)=wts(i)/norm_wts
!   enddo
!   wts=1.0d0
!   write(*,*) "THE VALUE OF WTS is: ", wts
!   qmat=0.0d0
!   do i=1,dim_ac
!      do l=1,dim_al
!         qmat(i)=qmat(i)+wts(l)*bmat(l)*amat(l,i)
!      enddo
!      do j=1,dim_ac
!         mat(i,j)=0.0d0
!         do l=1,dim_al
!            mat(i,j)=mat(i,j)+wts(l)*amat(l,i)*amat(l,j)
!         enddo
!      enddo
!   enddo
!   mat_inverse=0.0d0
!   do i=1,dim_ac
!      mat_inverse(i,i)=1.0d0
!   enddo
   !call inverse_real(mat,mat_inverse,dim_ac)
   !fc=matmul(mat_inverse,qmat)
   !do i=1,dim_ac
   !write(*,*) "THE VALUE OF FC_WEIGHTED IS: ", fc(i)
   !enddo
 !  deallocate(wts,qmat,mat,mat_inverse)
   do i=1,dim_al
      write(*,*) "AMAT NO SQ VAL: ",amat(i,:)
   enddo

   do i=1,dim_ac
      write(*,*) "AMAT SQ VAL: ", mat(i,:)
   enddo

   do i=1,dim_ac
      write(*,*) "BMAT VAL: ", bmat(i)
   enddo

   do i=1,dim_ac
      write(*,*) "QMAT VAL: ", qmat(i)
   enddo

   allocate(fc(dim_ac))
   call svd_set(dim_ac,dim_ac,mat,qmat,fc,sigma,svdcut,error,ermax,sig,'svd-all.dat') ! JUST SWITCH IT ON AND SEE THE RESULT

   do i=1,dim_ac
   write(*,*) "THE VALUE OF FC_WEIGHTED IS: ", fc(i)
   enddo

  call svd_set(dim_al,dim_ac,amat,bmat,fcs,sigma,svdcut,error,ermax,sig,'svd-all.dat')

  write(ulog,*)'After svd, || F_dft-F_fit || / || F_dft || =',sig

  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' TIME after SVD calculation                   IS ',tim

  call write_independent_fcs(sd)

  call write_output_fc2

  call date_and_time(date=today,time=now,zone=zone)
  call cpu_time(tim)
  write(utimes,'(a,f10.4)')' ENDING TIME OF THE PROGRAM                   IS ',tim
  write(ulog,'(a,f10.4)')' ENDING TIME OF THE PROGRAM                   IS ',tim

! if(enforce_inv.ne.0) then
!    call eliminate_fcs(sig)
!    write(ulog,*)'After ELIMINATE_FCS, sigma=',sig
! endif
deallocate(wts,qmat,mat,mat_inverse,rm_zero_energy,fc)
7 format(a,2(2x,g12.5),2x,4(1x,g11.4))

 write(ulog,'(a30,3x,a10,3x,a12,3x,a)')' Program  FC234 ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

 write(utimes,'(a30,3x,a10,3x,a12,3x,a)')' Program FC234 ended at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone

  write(ulog,7)' Average and largest errors in force=',error,ermax,sd

  close(ulog)
  close(ufit1)
  close(ufit2)
  close(ufit3)
  close(ufit4)
  close(ufc1)
  close(ufc2)
  close(ufc3)
  close(ufc4)
  close(ufco)
  close(ucor)
  close(umap)
  close(umatrx)
  close(utimes)

end program FC234
!===========================================================
 subroutine remove_zeros(nl,nc, amatr, bmatr,j)
! removes lines which are all zero from the amatrix
 use ios
 use params
 use svd_stuff  !for itrans,irot,ihuang
 implicit none
 integer nl,nc,i,j,nz,tra,rot,hua!,nz_index
 real(8) amatr(nl,nc),bmatr(nl),suml,small
 real(8), allocatable:: aa(:,:),bb(:)
 !real(8), allocatable:: rm_zero_energy(:) ! This is what I would like to create an array that removes energy for the lines whose bmat is zero
 integer, allocatable:: zz(:)
!write(*,*) "THE SIZE OF AMATR and NL is: ", size(amatr), NL
 small=1d-7 !10  displacements or cos(t) are never that small
 allocate(zz(nl))
 zz = 1
 nz = 0
 !nz_index=0
 do i=1,nl
! get the norm1 of each line
!   suml = 0
!   do j=1,nc
!      if ( abs(amatr(i,j)) .ge. suml ) suml = abs(amatr(i,j))
!   enddo
    suml = maxval(abs(amatr(i,:)))     ! sqrt(sum(amatr(i,:)*amatr(i,:))/nc)
    if (suml.lt.small) then
       if (abs(bmatr(i)) .lt. small) then
          if(verbose) write(ulog,*)' line #',i,' of Amatrix and Bmatrix is zero; it will be removed!'
       else
          write(ulog,*)' Inconsistency in line i of Amatrix which is zero!',i,suml
          write(ulog,*)' It will be removed Corresponding line in bmatr is=',bmatr(i)
          write(ulog,*)' Need to increase your range of FCs'
       endif
! record the line index for later removal
       zz(i) = 0
       nz = nz + 1
    endif
 enddo

!allocate(rm_zero_energy(nl-nz))
!nz_index=1
!do i=1,nl
!   suml = maxval(abs(amatr(i,:)))
!   if (suml .gt. small) then
!      rm_zero_energy(nz_index)=energy(i)
!      nz_index=nz_index+1
!   endif
!enddo

 write(ulog,*)' number of lines=0 in amat are=',nz
 allocate(aa(nl-nz,nc),bb(nl-nz))

 tra=0; rot=0; hua=0
 if (itrans .ne. 0) tra= transl_constraints
 if (irot   .ne. 0) rot= rot_constraints
 if (ihuang .ne. 0) hua= huang_constraints

 j = 0
 do i=1,nl
    if (zz(i).ne.0) then
       j = j+1
       aa(j,:) = amatr(i,:)
       bb(j) = bmatr(i)
    elseif(i.le.tra) then
       transl_constraints=transl_constraints-1
    elseif(i.le.tra+rot) then
       rot_constraints=rot_constraints-1
    elseif(i.le.tra+rot+hua) then
       huang_constraints=huang_constraints-1
    else
       force_constraints=force_constraints-1
    endif
 enddo
 write(ulog,4)'NEW tra,rot,hua,force_constaints=',transl_constraints,  &
 &       rot_constraints,huang_constraints,force_constraints

 amatr(1:j,:) = aa(1:j,:)
 bmatr(1:j)   = bb(1:j)
 amatr(j+1:nl,:) = 0
 bmatr(j+1:nl)   = 0
 if (j.ne.nl-nz) then
    write(ulog,*)'REMOVE ZEROS: Inconsistency! j.ne. nl.nz ',j,nl-nz
    stop
 endif
 write(ulog,*)'REMOVE ZEROS: new number of lines in amat=',j
 deallocate(aa,bb)
 4 format(a,4(i9))
 end subroutine remove_zeros
!===========================================================
 subroutine eliminate_fcs(sig)
! try to recalculate using elimination
 use svd_stuff
 use params
 use ios
 implicit none
 integer i,n_elim
 integer, allocatable:: fcmap(:)
 real(8), allocatable:: newa(:,:),newb(:,:),newc(:,:),newd(:,:),ainv(:,:),sig_new(:),fc_new(:),bmat_new(:),fc_elim(:)
 real(8) error,ermax,sd(4),sig,num,denom
 real tim

 if(inv_constraints.gt.dim_ac/2) then
    if(transl_constraints.gt.dim_ac/2) then
       return
    else
       n_elim=transl_constraints
    endif
 else
    n_elim=inv_constraints
 endif

 write(*   ,*)'number of constraints to be eliminated=',n_elim
 write(ulog,*)'number of constraints to be eliminated=',n_elim


 allocate (fcmap(dim_ac))
 i=n_elim
 allocate (newa(i,i),newb(i,dim_ac-i),newc(dim_al-i,i),newd(dim_al-i,dim_ac-i),ainv(i,i))
 allocate (fc_new(dim_ac-i),sig_new(dim_ac-i),bmat_new(dim_al-i),fc_elim(i))

 bmat_new=bmat(1+n_elim:dim_al)
 write(ulog,*)'Eliminate: dim_al,dim_ac,inv_constraints=',dim_al,dim_ac,inv_constraints

! decide which fcs to eliminate
! call sort(dim_ac,fcs,fcmap,dim_ac)
! keep the smallest ones and adjust the largest fcs to keep invariance relations satisfied.
 call invsort(dim_ac,fcs,fcmap)


! eliminate n_elim largest fcs defining A=[newa newb]  where newa(n_elim,n_elim) is inverted and substituted
!                                         [newc newd]
 do i=1,n_elim
    newa(1:n_elim       ,i)=amat(1:n_elim       ,fcmap(i))
    newc(1:dim_al-n_elim,i)=amat(1+n_elim:dim_al,fcmap(i))
 enddo
 do i=1+n_elim,dim_ac
    newb(1:n_elim       ,i-n_elim)=amat(1:n_elim       ,fcmap(i))
    newd(1:dim_al-n_elim,i-n_elim)=amat(1+n_elim:dim_al,fcmap(i))
 enddo

! now do the elimination by SVD-solving [D-C(1/A)B]y=b
 call inverse_real(newa,ainv,n_elim)

 newd=newd-matmul(newc,matmul(ainv,newb))

 call svd_set(dim_al-n_elim,dim_ac-n_elim,newd,bmat_new,fc_new,  &
&             sig_new,svdcut,error,ermax,sig,'svd-aux.out')

 deallocate(newa,newc,newd,bmat_new)

 call cpu_time(tim)
 write(utimes,'(a,f10.4)')' TIME after second SVD calculation             IS ',tim

! get the remaining FCs from x=-A_inv B y
 fc_elim=matmul(matmul(ainv,newb),fc_new)

! put back in the old order
 do i=1,n_elim
    fcs(fcmap(i))=fc_elim(i)          ! I DO NOT UNDERSTAND THIS PART - NO CHANGE HERE JUST WANT TO KNOW WHAT IS BEING DONE
 enddo
 do i=1+n_elim,dim_ac
    fcs(fcmap(i))=fc_new(i-n_elim)   ! AND WHAT IS THIS PART WHAT ARE WE DOING IT HERE - NO CHANGE HERE JUST WANT TO KNOW WHAT IS BEING DONE
 enddo

 call write_independent_fcs(sd)

 call write_output_fc2

 num=0; denom=0
 do i=1,dim_al
    num=num+(dot_product(amat(i,:),fcs(:))-bmat(i))**2
    denom=denom+(bmat(i))**2
 enddo
 sig=sqrt(num/denom)

 deallocate(newb,ainv,fc_elim,sig_new,fc_new)

 end subroutine eliminate_fcs
!============================================================
 subroutine invsort(n,r,mcor)
! sorts the first n elements of array r in descending order r(mcor(i)) is the ordered array
  implicit none
  integer n,i,j,temp
  real(8) r(n)
  integer mcor(n)

  do i=1,n
     mcor(i)=i
  enddo
  do i=1,n  ! was 1,n-1
     do j=i+1,n  ! was i+1,n
        if(r(mcor(i)).lt.r(mcor(j))) then
           temp=mcor(i)
           mcor(i)=mcor(j)
           mcor(j)=temp
        endif
     enddo
  enddo
!  r(mcor) is now the ordered array r(mcor(1)) > r(mcor(2)) > ... > r(mcor(n))
 end subroutine invsort