!===========================================================
 subroutine read_structure(fn1)
! reads the param.inp file containing info on the atom types masses
! and coordinates within the primitive cell
 use ios
 use params
 use lattice
 use atoms_force_constants
 use svd_stuff
 use born, only : born_flag
 use constants, only : r15
 implicit none
 integer i,counter,tau ,rnk,ios_cc,n_cc
 real(r15) scal,force_error
 character fdf*1,r4*7,invr*4,incl*8,it*2,zone*5,now*10,today*8,bornflg*3
 character, intent(in):: fn1*(*)
 character line3*200
 real tim

 open(uparams,file='default.params',status='old')
   read(uparams,*) tolerance    ! cutoff for smallest "eigenvalue" w to be included
   if(tolerance.eq.0) tolerance = 1d-4
!   read(uparams,*) svdcut    ! cutoff for smallest "eigenvalue" w to be included
!   if(svdcut.eq.0) svdcut = 1d-10   ! default values
   read(uparams,*) force_error
   if(force_error.eq.0) force_error=1d-5   ! default values
   svdcut=force_error*force_error  ! also = inverse of condition number
   read(uparams,*) icutoff    ! cutoff for neighborshells
   if(icutoff.eq.0) icutoff = 6     ! default multiple of bond length; defines maxatoms
   read(uparams,*) margin    ! cutoff for neighborshells
   if(margin.eq.0) margin = 1d-6     ! default margin to eliminate FCs
! line 5 (optional, absent = 0): 1 = relax the reference structure with the fitted FCs and re-expand the
! FCs about it (relax_reference.f90); everything after the fit then uses the relaxed structure
! (only the value 1 switches it on: older long default.params files have other numbers on line 5)
   read(uparams,*,iostat=ios_cc) relax_ref
   if(ios_cc.ne.0 .or. relax_ref.ne.1) relax_ref = 0
 close(uparams)

 open(uparams,file='structure.params',status='old')

 read(uparams,*) latticeparameters   ! (a,b,c,alpha,beta,gamma)
 read(uparams,*) primitivelattice    ! prim latt vectors in terms of conventional above
 read(uparams,'(a)') line3   ! scale factor for lattparams, then optionally convcoord
! Read line 3 as text so that the old one-value form ("1.0 scale") still works:
! convcoord=0 (coordinates in conventional units) unless an integer follows scal.
 read(line3,*) scal
 read(line3,*,iostat=ios_cc) scal, convcoord
 if(ios_cc.ne.0) convcoord=0
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
 read(uparams,'(a)') line3  ! fc2range [fc3flag ... fc8flag]; absent flags = 1 (nshells, as before)
 do n_cc = 7, 1, -1
    fcflag = 1
    read(line3,*,iostat=ios_cc) fc2range, fcflag(3:n_cc+1)
    if(ios_cc.eq.0) exit
 enddo
 if(ios_cc.ne.0) then
    write(*,*)'READ_STRUCTURE: line 12 must start with fc2range: ',trim(line3)
    stop
 endif
 allocate(nshells(8,natom_prim_cell))
 nshells(1,1:natom_prim_cell)=1
 read(uparams,*) nshells(2,1:natom_prim_cell)
 read(uparams,*) nshells(3,1:natom_prim_cell)
 read(uparams,*) nshells(4,1:natom_prim_cell)
 read(uparams,*) nshells(5,1:natom_prim_cell)
 read(uparams,*) nshells(6,1:natom_prim_cell)
 read(uparams,*) nshells(7,1:natom_prim_cell)
 read(uparams,*) nshells(8,1:natom_prim_cell)

! include_fc(rank)=2: those FCs are read by read_fcs (after the FC2 range is set) and go to the right-hand sides

! just for naming of log file
 open(321,file='dielectric.params',status='old')
 read(321,*)born_flag
 close(321)
! if(bornflag.le.0) then ! no Born correction
!    born='B0'
! elseif(bornflag.eq.1) then
!    born='BA'  ! A for add
! elseif(bornflag.eq.2) then
!    born='BE'  ! E for EwaldE
! else
!    born='BS'  ! B for Born (3 or 4)a ; S for subtract
! endif
 bornflg(1:1)='B'
 write(bornflg(2:3),'(i2.2)')born_flag

 if (itemp.eq.0) then
    if(fc2range.eq.0) then
       it='df'  ! for default
    else
       it='rd'  ! for read
    endif
 elseif(itemp.eq.1) then ! it='T'
    if(fc2range.eq.0) then
       it='Td'   ! temperature+default
    else
       it='Tr'  ! temperature + read 
    endif
 endif
 write(fdf,'(i1)')fdfiles
 if(include_fc(2).eq.1 .and. fc2range.eq.1) then
    write(r4(1:2),'(i2.2)')nshells(2,1)
 else
    r4(1:2)='00'
 endif
 if(include_fc(3).eq.1) then
    write(r4(3:5),'(a1,i2.2)')'_',nshells(3,1)
 else
    r4(3:5)='_00'
 endif
 if(include_fc(4).eq.1) then
    write(r4(6:7),'(a1,i1.1)')'_',nshells(4,1)
 else
    r4(6:7)='_0'
 endif
! write(incl,'(8i1)')include_fc(:)
 invr='0000'
 if (itrans.ne.0) invr(1:1)='t'
 if (irot.ne.0) invr(2:2)='r'
 if (ihuang.ne.0) invr(3:3)='h'
 if (enforce_inv.eq.1) invr(4:4)='E'
 if (enforce_inv.eq.2) invr(4:4)='L'   ! L1 / LASSO
 if (enforce_inv.eq.3) invr(4:4)='R'   ! LASSO + least-squares refit
 if (enforce_inv.lt.0 .or. enforce_inv.gt.3) then
    write(*,*)'READ_STRUCTURE: enforce_inv must be 0 (SVD on the stacked system),'
    write(*,*)' 1 (constrained fit in the kernel), 2 (LASSO) or 3 (LASSO + refit); read ',enforce_inv
    stop
 endif
 open(ulog  ,file=fn1//fdf//it//bornflg//'_'//r4//'_'//invr//'.dat'   ,status='unknown')

 call date_and_time(date=today,time=now,zone=zone)
 call cpu_time(tim)
 write(ulog,'(a,3x,a,3x,a,3x,a)')' Program FOCEX was launched at ',today(1:4)//'/'  &
&  //today(5:6)//'/'//today(7:8),now(1:2)//':'//now(3:4)//':'//now(5:10),zone
 write(ulog,'(a,f10.4)')' STARTING TIME OF THE PROGRAM                   IS ',tim
 write(ulog,*)'===================================================================='

 write(ulog,*) 'DEFAULT VALUES:'
 write(ulog,*) svdcut   ,'   cutoff for smallest eigenvalue w to be included for inversion'
 write(ulog,*) tolerance,'   tolerance(Ang) for equating two variables '
 write(ulog,*) rcutoff  ,'   cutoff multiple of dist/atom for neighbors of primitive cell '
 write(ulog,*)' ---------------------------------------------------------'
 write(ulog,*) fc2range,'  if 0 default range consistent with the largest supercell is chosen for FC2 '
 write(ulog,*) fcflag(3:8),'  fc3flag..fc8flag: 0 = FC2-range rule for 2-body and aligned clusters, else nshells'
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
 write(ulog,*)'Atoms: tau, name,type,mass,reduced coordinates (convcoord=',convcoord,')'
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
!==========================================================
 subroutine check_born_symmetry
!! Born charges and dielectric tensor against the space group (2026-10-23, user's request):
!! for every operation {R|t} and atom k with R x_k + t = x_k' + lattice vector, the tensors must
!! obey Z(k') = R Z(k) R^T (this includes the site symmetry, k'=k) and eps = R eps R^T. Both
!! indices of charge(:,:) are Cartesian, so the test does not depend on their order.
!! Deviation relative to max|Z| (resp. max|eps|): > 5e-2 stops, > 1e-4 gives a warning. Tensors on
!! the wrong atoms give O(1) (BZO: 0.50); DFT tensors that are not symmetrized give ~1e-2 (GeSe:
!! equivalent Ge with Z_yy 2.598 and 2.630, 1.0e-2), which only warns.
!! On a stop, every reassignment of the Z* blocks among atoms of the same species is tried; if one
!! obeys the symmetry, the reordered file is written to dielectric.params.symmetric and the
!! equivalent atom order of structure.params is printed (nothing is changed in the run).
!! Called after symmetry_init, only if born_flag /= 0.
 use ios, only : ulog
 use born, only : epsil, born_flag
 use lattice, only : cart_to_prim, prim_to_cart
 use atoms_force_constants, only : atom0, natom_prim_cell, op_matrix, isgop, isgopcount, sgfract
 use geometry, only : v2a
 use constants, only : r15
 use params, only : tol_sym
 implicit none
 integer isg, k, kp, kbad, kpbad, opbad, n, i, j, nperm, nfit
 integer, allocatable :: img(:,:), perm(:), best(:), idx(:)
 real(r15) r(3,3), d(3), zmax, emax, devmax, edevmax, dev, bestdev
 real(r15), parameter :: tol_stop = 5d-2, tol_warn = 1d-4
 integer, parameter :: maxperm = 1000000
 logical more

 if (born_flag .eq. 0) return
 n = natom_prim_cell
 allocate(img(n,isgopcount), perm(n), best(n), idx(n))
! img(k,isg): the atom k is taken to by operation isg
 do isg = 1, isgopcount
    r = op_matrix(:,:,isgop(isg))
    do k = 1, n
       img(k,isg) = 0
       do kp = 1, n
          if (atom0(kp)%at_type .ne. atom0(k)%at_type) cycle
          d = matmul(cart_to_prim, matmul(r, v2a(atom0(k)%equilibrium_pos)) + sgfract(:,isg) &
&                                  - v2a(atom0(kp)%equilibrium_pos))
          if (sqrt(sum(matmul(prim_to_cart, d - nint(d))**2)) .le. tol_sym) then
             img(k,isg) = kp ; exit
          endif
       enddo
       if (img(k,isg) .eq. 0) then
          write(ulog,*)'CHECK_BORN_SYMMETRY: no image of atom ',k,' under op ',isg
          error stop 'CHECK_BORN_SYMMETRY: atom not mapped by a space-group operation'
       endif
    enddo
 enddo
 zmax = 0
 do k = 1, n
    zmax = max(zmax, maxval(abs(atom0(k)%charge)))
 enddo
 if (zmax .eq. 0) zmax = 1
 emax = maxval(abs(epsil))
 edevmax = 0
 do isg = 1, isgopcount
    r = op_matrix(:,:,isgop(isg))
    edevmax = max(edevmax, maxval(abs(matmul(matmul(r, epsil), transpose(r)) - epsil)) / emax)
 enddo
 do k = 1, n
    perm(k) = k
 enddo
 devmax = deviation(perm, .true.)
 write(ulog,'(a,2es10.2)')' CHECK_BORN_SYMMETRY: max deviation of Z*, eps from the space group =', &
&     devmax, edevmax
 if (devmax .gt. tol_warn) then
    r = op_matrix(:,:,isgop(opbad))
    write(ulog,'(a,i3,a,i3,a,i3)')' worst: atom ',kpbad,' as image of atom ',kbad,' under op ',opbad
    write(ulog,*)' given Z* of atom ',kpbad,' (rows):'
    write(ulog,'(3f10.4)') transpose(atom0(kpbad)%charge)
    write(ulog,*)' R Z*(atom ',kbad,') R^T (rows):'
    write(ulog,'(3f10.4)') transpose(matmul(matmul(r, atom0(kbad)%charge), transpose(r)))
 endif
 if (devmax .gt. tol_stop .or. edevmax .gt. tol_stop) then
    write(*,'(a,2es10.2)')' CHECK_BORN_SYMMETRY: Z*, eps deviate from the crystal symmetry by', &
&        devmax, edevmax
    write(*,*)' dielectric.params: one block per atom in the order of structure.params, in the'
    write(*,*)' frame of structure.params (a along x, b in the xy plane); the log has the worst atom'
    if (edevmax .le. tol_stop) call search_order
    error stop 'CHECK_BORN_SYMMETRY: Born charges or dielectric tensor break the crystal symmetry'
 elseif (devmax .gt. tol_warn .or. edevmax .gt. tol_warn) then
    write(*,'(a,2es10.2)')' WARNING: Z*, eps deviate from the crystal symmetry by ',devmax,edevmax
    write(ulog,'(a,2es10.2)')' WARNING: Z*, eps deviate from the crystal symmetry by ',devmax,edevmax
 endif
 deallocate(img, perm, best, idx)

 contains
!----------------------------------------
 function deviation(p, track) result(dm)
!! max over ops and atoms of |Z(p(k')) - R Z(p(k)) R^T| / zmax, atom k carrying block p(k)
 integer, intent(in) :: p(:)
 logical, intent(in) :: track
 real(r15) dm, dv, rr(3,3)
 integer is, kk
 dm = 0
 do is = 1, isgopcount
    rr = op_matrix(:,:,isgop(is))
    do kk = 1, n
       dv = maxval(abs(atom0(p(img(kk,is)))%charge &
&                    - matmul(matmul(rr, atom0(p(kk))%charge), transpose(rr)))) / zmax
       if (dv .gt. dm) then
          dm = dv
          if (track) then
             kbad = kk ; kpbad = img(kk,is) ; opbad = is
          endif
       endif
    enddo
 enddo
 end function deviation
!----------------------------------------
 subroutine search_order
!! all reassignments of the blocks within each species (odometer of lexicographic permutations)
 real(r15) x
 nperm = 1
 do k = 1, n
    j = count(atom0(1:k)%at_type .eq. atom0(k)%at_type)     ! factorial of the species counts
    nperm = nperm * j
    if (nperm .gt. maxperm) then
       write(*,*)' too many reorderings to try (> ',maxperm,'): fix dielectric.params by hand'
       return
    endif
 enddo
 do k = 1, n
    perm(k) = k
 enddo
 bestdev = huge(1d0) ; nfit = 0
 more = .true.
 do while (more)
    x = deviation(perm, .false.)
    if (x .lt. tol_warn) nfit = nfit + 1
    if (x .lt. bestdev) then
       bestdev = x ; best = perm
    endif
    more = next_assignment(perm)
 enddo
 write(ulog,'(a,i8,a,es10.2,a,i6)')' CHECK_BORN_SYMMETRY: tried ',nperm,' block orders; best ', &
&     bestdev,'; orders obeying the symmetry: ',nfit
 if (bestdev .ge. tol_warn) then
    write(*,*)' no reordering of the Z* blocks within a species obeys the symmetry: the tensors'
    write(*,*)' themselves (values, frame or transposition) are wrong'
    return
 endif
 write(*,'(a,30i4)')' the symmetry is obeyed with the blocks of dielectric.params in the order', best
 write(*,*)' (atom k gets the block now in position best(k)); written to dielectric.params.symmetric'
 do k = 1, n               ! idx: structure.params order that fits the present dielectric.params
    idx(best(k)) = k
 enddo
 write(*,'(a,30i4)')' equivalently, keep dielectric.params and list the atoms of structure.params as', idx
 write(ulog,*)' CHECK_BORN_SYMMETRY: block order ',best,'; or atom order ',idx
 if (nfit .gt. 1) write(*,*)' (',nfit,' orders obey the symmetry; they give the same physics)'
 call write_reordered(best)
 end subroutine search_order
!----------------------------------------
 logical function next_assignment(p)
!! next permutation of the blocks among atoms of the same species: lexicographic next
!! permutation of each species' positions, odometer over species; false after the last one
 integer, intent(inout) :: p(:)
 integer t, m, sp(n), q(n), a, b, tmp
 next_assignment = .false.
 do t = maxval(atom0(1:n)%at_type), 1, -1
    m = 0
    do k = 1, n
       if (atom0(k)%at_type .eq. t) then
          m = m + 1 ; sp(m) = k ; q(m) = p(k)
       endif
    enddo
    if (m .lt. 2) cycle
    a = m - 1                                   ! lexicographic next permutation of q(1:m)
    do while (a .ge. 1)
       if (q(a) .lt. q(a+1)) exit
       a = a - 1
    enddo
    if (a .ge. 1) then
       b = m
       do while (q(b) .le. q(a))
          b = b - 1
       enddo
       tmp = q(a) ; q(a) = q(b) ; q(b) = tmp
       q(a+1:m) = q(m:a+1:-1)
       p(sp(1:m)) = q(1:m)
       next_assignment = .true.
       return
    endif
    q(1:m) = q(m:1:-1)                          ! wrap this species to its first order, carry on
    p(sp(1:m)) = q(1:m)
 enddo
 end function next_assignment
!----------------------------------------
 subroutine write_reordered(p)
!! dielectric.params with the Z* blocks (lines 5.., 3 per atom) in the order p; other lines kept
 integer, intent(in) :: p(:)
 character(len=400) line
 character(len=400), allocatable :: lines(:)
 integer u, nl, ios_, ll
 open(newunit=u, file='dielectric.params', status='old')
 nl = 0
 do
    read(u,'(a)',iostat=ios_) line
    if (ios_ .ne. 0) exit
    nl = nl + 1
 enddo
 allocate(lines(nl))
 rewind(u)
 do ll = 1, nl
    read(u,'(a)') lines(ll)
 enddo
 close(u)
 if (nl .lt. 4 + 3*n) then
    write(*,*)' dielectric.params has fewer lines than expected; reordered file not written'
    return
 endif
 open(newunit=u, file='dielectric.params.symmetric', status='replace')
 do ll = 1, 4
    write(u,'(a)') trim(lines(ll))
 enddo
 do k = 1, n
    do ll = 1, 3
       write(u,'(a)') trim(lines(4 + 3*(p(k)-1) + ll))
    enddo
 enddo
 do ll = 5 + 3*n, nl
    write(u,'(a)') trim(lines(ll))
 enddo
 close(u)
 end subroutine write_reordered

 end subroutine check_born_symmetry
!==========================================================
 subroutine read_latdyn
!! reads information on kpoint mesh within the full BZ, their shift, temp range, whether classical or quantum
 use ios
 use om_dos
 use params
 use lattice
 use kpoints
 use constants, only : r15
 use atoms_force_constants
 use geometry, only : length
 implicit none
 real(r15) junk
 character line*99

  open(uparams,file='latdyn.params',status='old')
  write(6,*) 'READ_PARAMS: opening latdyn.params'
  read(uparams,*)nc
  write(*,*) 'READ_PARAMS: kmesh= ',nc
  nkc=nc(1)*nc(2)*nc(3)
  read(uparams,*)shift, normal
  normal=normal/length(normal) 
  write(*,4) 'READ_PARAMS: shift, normal = ',shift,normal 
  read(uparams,*)wmesh,wmax ! number of freq meshes for DOS and its upperbound 
  write(*,*) 'READ_PARAMS: wmesh, wmax = ',wmesh,wmax
  read(uparams,*)width  ! width of gaussian broadening for DOS
  write(*,*) 'READ_PARAMS: width= ',width
  read(uparams,*) line
  write(*,*) 'READ_PARAMS: just read verbose ',verbose
  read(uparams,*)tmin,tmax !,ntemp      ! temps in Kelvin to calculate thermal expansion and other thermal ppties
  write(*,*) 'READ_PARAMS: tmin,tmax=',tmin,tmax !,ntemp
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
4 format(a,6(1x,f9.4))
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
 integer i,ierr
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

! Number of atoms of each type.
! VASP-5 and newer write a line of ELEMENT NAMES immediately before the counts;
! VASP-4 wrote the counts straight after the third lattice vector. Accept both:
! read the line as text and try to parse the counts out of it. Element names
! cannot be read as integers, so a failed parse means this was the name line and
! the counts are on the next one. Until 2026-08-14 only the VASP-4 layout was
! accepted, which is why the docs used to tell the user to delete line 6 of
! sc_snaps.x's poscar_000 -- no longer necessary.
 read(uposcar,'(a)') line
 read(line,*,iostat=ierr)(natom(i),i=1,natom_type)
 if (ierr.ne.0) then       ! that was the element-name line; counts follow
    read(uposcar,'(a)') line
    read(line,*,iostat=ierr)(natom(i),i=1,natom_type)
 endif
 if (ierr.ne.0) then
    write(ulog,*)'READ_SUPERCELL: could not read ',natom_type,' atom counts from ',poscar
    write(ulog,'(a)')' offending line: '//trim(line)
    write(ulog,*)' expected one integer per atom type, in the order of structure.params'
    write(*   ,*)'READ_SUPERCELL: could not read ',natom_type,' atom counts from ',poscar
    write(*   ,'(a)')' offending line: '//trim(line)
    stop
 endif
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


  if (abs(volume_r/natom_super_cell-volume_r0/natom_prim_cell) .gt. 1d-5*volume_r0) then
     write(ulog,*)' supercell inconsistency; check input coordinates again'
     write(ulog,5)' natom_prim_cell, volume_r0,density0=',natom_prim_cell,volume_r0,natom_prim_cell/volume_r0
     write(ulog,5)' natom_sc, volume_r,density=',natom_super_cell,volume_r,natom_super_cell/volume_r
     stop
  endif

 write(ulog,*)'RECIPROCAL_LATTICE: '
 call make_reciprocal_lattice_2pi(rs1,rs2,rs3,gs1,gs2,gs3)
 call write_out(ulog,'supercell volume ',volume_r)
 call write_out(ulog,'gsc1 ',gs1)
 call write_out(ulog,'gsc2 ',gs2)
 call write_out(ulog,'gsc3 ',gs3)

4 format(9(2x,f19.9))
5 format(a,i6,9(2x,f19.9))
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
 integer i,j !,n(3)
 real(r15) a(3)

! check mapping and assign label and address to each atom in the supercell
! this seems to be working for now....
 call check_input_poscar_consistency_new

! call write_out(6,' cart_to_prim ',cart_to_prim)
! index of each supercell atom (n,tau) in the atompos list; 0 if beyond the shells of atompos
 do i=1,natom_super_cell
    atom_sc(i)%cell%atomposindx = 0
    do j=1,natoms
       if (iatomcell0(j).eq.atom_sc(i)%cell%tau .and. all(iatomcell(:,j).eq.atom_sc(i)%cell%n)) then
          atom_sc(i)%cell%atomposindx = j
          exit
       endif
    enddo
 enddo

 write(ulog,*)' ASSIGNMENTS WERE DONE FOR ATOMS IN THE SUPER CELL '
! write(ulog,*)'i_sc,type(i),tau(i)      n(i)    ,                r_cart(i)              r_red(i) '
  do i=1,natom_super_cell
    write(ulog,8)' ',i,atom_sc(i)%at_type, &
&     atom_sc(i)%cell%tau,atom_sc(i)%cell%n,atom_sc(i)%equilibrium_pos,  &
&     matmul(cart_to_prim,v2a(atom_sc(i)%equilibrium_pos))
 enddo

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
 integer i,j,k,counter,ier,nt(3),n1,isave,itry
 real(r15) a(3),invn(3,3)
 type(vector) shift0,vec
 logical matched,all_matched,typed,typed_failed
 integer t,k0

! first describe r_i in terms of r0_i, and make sure the lin comb is integer
 write(ulog,*)' ----------------------------------------------------------'
 write(ulog,*)' CHECKING the commensurability between primcell and supercell'
 call check_int(rs1,a,ier,g01,g02,g03)
 write(ulog,4)'checking rs1 in units of r0i,ier=',a,ier
 if (ier .eq. 1) call stop_not_commensurate(1,a)
 n_sc(:,1)=nint(a)
 write(ulog,*)' direct coords of rs1=',n_sc(:,1)
 call check_int(rs2,a,ier,g01,g02,g03)
 write(ulog,4)'checking rs2 in units of r0i,ier=',a,ier
 if (ier .eq. 1) call stop_not_commensurate(2,a)
 n_sc(:,2)=nint(a)
 write(ulog,*)' direct coords of rs2=',n_sc(:,2)
 call check_int(rs3,a,ier,g01,g02,g03)
 write(ulog,4)'checking rs3 in units of r0i,ier=',a,ier
 if (ier .eq. 1) call stop_not_commensurate(3,a)
 n_sc(:,3)=nint(a)
 write(ulog,*)' direct coords of rs3=',n_sc(:,3)

 write(ulog,*)' COMMENSURABILITY CHECKED '
 write(ulog,*)' natom_super_cell=',natom_super_cell
 call write_out(ulog,' N vectors of the supercell ',dble(n_sc))
 call xmatinv(3,dble(n_sc),invn,ier)
 call write_out(ulog,' Inverse of n_sc matrix ',invn)
 if(ier.ne.0) then
    write(*,*)'n_sc matrix inversion returned error ',ier
    stop
 else
    invn_sc=nint(invn*det(n_sc))
 endif
 call write_out(ulog,' NINTeger part of inverse of n_sc matrix ',invn_sc)

 write(ulog,*)'FINDING possible translation vectors, trying them on other atoms'

! First try (typed): every POSCAR atom must also sit on a site of the species its counts line
! gives it -- this picks the right origin when the positions alone cannot (MgO: both rocksalt
! sublattices fit). If that fails, the second try (species from the positions only) is made just to
! report which atoms are which species; the run then stops (user's rule, 2026-10-24: a failed
! consistency check stops). MoTe2 test POSCAR: counts 64 128 (Mo Te, the order of structure.params)
! but the atoms were 128 Te then 64 Mo.
 typed_failed = .false.
 typeloop: do itry=1,2
 typed = itry.eq.1
! write(*,*)'# ier,i0,j0,k_sc; vec,vred='
 checkloop: do i=1,natom_prim_cell !,1,-1
    all_matched=.False.
    shift0 = atom0(i)%equilibrium_pos - atom_sc(1)%equilibrium_pos
    write(ulog,7)'atom in PRIMCELL & SC:',i,atom0(i)%at_type,atom0(i)%equilibrium_pos,atom_sc(1)%equilibrium_pos
    write(ulog,11)'trying shift vector to match atom_sc(1)=',i,v2a(shift0),cart2red(v2a(shift0),'r')
    if (typed .and. atom0(i)%at_type .ne. poscar_type(1)) cycle checkloop

! does it fall on any of the primitive cell atoms?
! try this shift see if all atoms in SC can be mapped to the PC by it.
    all_matched=.True.
    SC: do k = 2,natom_super_cell
       matched = .false.
       PRIM: do j = 1,natom_prim_cell  ! one of the prim-cell atoms has to match
          if (typed .and. atom0(j)%at_type .ne. poscar_type(k)) cycle PRIM
          vec =  shift0 + atom_sc(k)%equilibrium_pos - atom0(j)%equilibrium_pos

! In first loop (line 67-72), replace with:
          call check_int(vec, a, ier, g01, g02, g03)
          if(ier .ne. 0) cycle PRIM
          matched = .True.
          exit PRIM

! find its direct coordinates on r01,r02,r03
   !      a(1)=vec.dot.g01/(2*pi); if(abs(a(1)-nint(a(1))).gt.1d-2)  cycle PRIM
   !      a(2)=vec.dot.g02/(2*pi); if(abs(a(2)-nint(a(2))).gt.1d-2)  cycle PRIM
   !      a(3)=vec.dot.g03/(2*pi); if(abs(a(3)-nint(a(3))).gt.1d-2)  cycle PRIM
!     write(*,'(4i4,a,3(1x,f9.4),4x,a,3(1x,f9.4))')ier,i,j,k,' er,i0,j0,k_sc; vec=',vec,'vred=',a
   !      matched=.True.
   !      exit PRIM

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
 if (all_matched) exit typeloop
 if (typed) typed_failed = .true.
 enddo typeloop

! shift = r_primcell(i)-r_supercell(1) of direct coordinates a

 if (all_matched) then
    call write_out(ulog,'THE shift vector ', shift0)
    call write_out(ulog,'It maps atom 1 in SC to primitive atom whose number ', isave)
    call write_out(6   ,'THE shift vector ', shift0)
    call write_out(6   ,'It maps atom 1 in SC to primitive atom whose number ', isave)
 else
    write(ulog,*)' NO SHIFT VECTOR COULD BE DEFINED!, check your coordinates again'
    write(6   ,*)' NO SHIFT VECTOR COULD BE DEFINED!, check your coordinates again'
    error stop 'NO SHIFT VECTOR COULD BE DEFINED'
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
            if (typed .and. atom0(j)%at_type .ne. poscar_type(k)) cycle PRIM2

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
            write(*,*)'SC ATOM ',k,' can not be mapped to the primitive cell'
            error stop 'SC ATOM can not be mapped to the primitive cell'
         elseif ( counter .ge. 2) then
            write(ulog,*)' CHECK: more than one atom matched! check your coordinates'
            write(*   ,*)' CHECK: more than one atom matched! check your coordinates'
            stop
         endif
       enddo SC2

! the atoms fit the primitive cell only with species other than the counts line gives: report, stop
 if (typed_failed) then
    do i = 1, 2
       if (i.eq.1) k0 = ulog
       if (i.eq.2) k0 = 6
       write(k0,*)'ERROR: the POSCAR atoms do not sit on sites of the species given by its counts line'
       write(k0,*)'       (counts in the species order of structure.params):'
       write(k0,'(a,20(1x,a,i6))')'        counts line  :',(trim(atom0(minloc(abs(atom0(:)%at_type-t),1))%name), &
&            natom(t), t=1,natom_type)
       write(k0,*)'       species of the atoms from their positions (atom ranges):'
       j = 1
       do k = 2, natom_super_cell + 1
          if (k .gt. natom_super_cell) then
             write(k0,'(9x,i6,a,i6,2x,a)') j,' -',k-1,trim(atom_sc(j)%name)
          elseif (atom_sc(k)%at_type .ne. atom_sc(j)%at_type) then
             write(k0,'(9x,i6,a,i6,2x,a)') j,' -',k-1,trim(atom_sc(j)%name)
             j = k
          endif
       enddo
       write(k0,*)'       Reorder the atoms of POSCARi and of every FORCEDISPi block so that the species'
       write(k0,*)'       follow the order of structure.params, with the counts line to match.'
    enddo
    error stop 'POSCAR species order does not match its counts line (see above)'
 endif

 write(ulog,*)'*******************************************************'
 write(ulog,*)'SINCE THE PROGRAM WAS NOT STOPPED, MAPPING-CHECK PASSED!'
 write(ulog,*)'NOW CHECKING CORRESPONDANCE WITH ATOMPOS ARRAY'
 write(ulog,*)'*******************************************************'

 write(*,*)'calling write_correspondance'
 call write_correspondance (isave,shift0)

2 format(i5,a,2(1x,i2),'(',3(i2),')',2(2x,3(1x,f8.3)))
3 format(9(1x,g12.5))
4 format(a,3(2x,f8.3),2x,i1)
5 format(a,1x,3(4x,i9),9(1x,f9.4))
6 format(a,i2,1x,2i5,9(1x,f9.4))
7 format(a,2(1x,i4),9(2x,f9.4))
8 format(a,3(1x,i4),9(2x,f9.4))
9 format(a,3(1x,f12.5),3x,i3,'(',3i2,')',i5)
11 format(a,i6,99(1x,f9.4))

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

 contains

 subroutine stop_not_commensurate(k,a)
 integer, intent(in) :: k
 real(r15), intent(in) :: a(3)
 write(*,'(a,i1,a,3f10.4)')' supercell vector ',k,' of POSCAR in units of the primitive vectors =',a
 write(*,*)'is not integer: supercell and primitive cell are not commensurate in this frame.'
 write(*,*)'If these are nearly integers, the cell built from structure.params line 1'
 write(*,*)'(a along x, b in the xy plane) is rotated relative to POSCAR. To keep the POSCAR'
 write(*,*)'frame use line 1 = 1 1 1 90 90 90, line 2 = the cartesian primitive vectors,'
 write(*,*)'line 3 = 1.0 1 and atomic coordinates in primitive units.'
 stop
 end subroutine stop_not_commensurate

 integer function poscar_type(k)
!! species of atom k of POSCARi, from its counts line (species in the order of structure.params)
 integer, intent(in) :: k
 integer t
 poscar_type=0
 do t=1,natom_type
    if (k .le. sum(natom(1:t))) then
       poscar_type=t
       return
    endif
 enddo
 end function poscar_type

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

! energy per primitive cell first, then subtract the lowest per-cell energy of FORCEDISP1, so that the
! snapshots of supercells of different sizes are on the same scale (user, 2026-11-02; until then emin
! was the lowest TOTAL energy of supercell 1, subtracted from the totals of every supercell before the
! division, which shifted a supercell of another size by about E0(1)/N_cells)
 engy=engy/natom_super_cell*natom_prim_cell  ! eV per primitive cell
 if( first_call ) then
    emin=minval(engy)
    first_call=.false.
 endif
 engy=engy-emin

 write(*,9)'Energy per primitive cell relative to the lowest of FORCEDISP1 (eV): ',engy
 write(ulog,'(a,a,i6,a,i4)')' Energy assigned (eV per primitive cell, relative to the lowest of FORCEDISP1): ', &
 &    outcar,natom_super_cell/natom_prim_cell,' cells'
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
!! sd(rank) = mean standard error of the irreducible FCs of that rank (sig = standard error of each,
!! from the least-squares fit: enforce_inv=0 svd_set, 1 fit_kernel_basis; 0 = not estimated: LASSO,
!! itemp=1)
 use svd_stuff
 use atoms_force_constants
 use params
 use constants, only : r15
 implicit none
 integer, intent(in) :: ulog,n
 real(r15), intent(in) :: sig(n)
 real(r15), intent(out):: sd(8)
 integer rnk,k,g,cnt2,cnt,ti,res


   write(ulog,*)' indepterm, fcs, sigma, normalized sigma for all SVD'
   do k=1,n
      write(ulog,4) k,fcs(k),sig(k)
   enddo
   write(ulog,*) "rank, average error"

!  k1 ----------------- for rank=2 we always work with map(2)%nkeptind
! if(fc2range.eq.0) then
   write(ulog,*)'# of groups of rank 2, ngr, versus # of kept indep terms =',map(2)%ngr,sum(map(2)%keep(:))
   write(ulog,*)'# of indep terms versus # of kept terms=',map(2)%ntotind,map(2)%nkeptind
! endif
!  k1 -----------------

   write(ulog,*)' dim_ac, size(fcs) =',dim_ac,n

   cnt2=0;sd=0; res=0
   do rnk=1,maxrank

      write(ulog,*)'rank,# of groups=',rnk,map(rnk)%ngr
      if(include_fc(rnk).ne.1) cycle
      cnt=0
      do g=1,map(rnk)%ngr
        do ti=1,map(rnk)%ntind(g)
            cnt=cnt+1   ! counter of terms of given rank
            if(rnk.eq.2) then
               cnt=counteri(2,g,ti)
               cnt2 = map(1)%nkeptind+ current(2,g,ti) 
            else 
               cnt2=cnt2+1 ! cumulative counter over all ranks
               if(cnt2.ne.res+current(rnk,g,ti)) then
                  call warn2(ulog,'cnt2 ne res+current')
                  write(ulog,*)'cnt2,res+current=',cnt2,res+current(rnk,g,ti)
                  stop
               endif
            endif
            sd(rnk)=sd(rnk)+sig(cnt2)
            if(verbose) write(ulog,6)'rank,group,nind,cnt(rnk),cntot,sd(rnk)=' &
&                       ,rnk,g,ti,cnt,cnt2,sig(cnt2),sd(rnk)
        enddo
      enddo
      if(cnt .ne. 0 ) sd(rnk)=sd(rnk)/cnt
      write(ulog,3)'Rank, sd(rank)=', rnk,sd(rnk)

      res=res+map(rnk)%nkeptind
   enddo
   write(ulog,5)' !==== sd(rank) = mean standard error of the FCs of each rank=1,8 ',(sd(rnk),rnk=1,maxrank)
   if(all(sig.eq.0)) write(ulog,*)' sd(rank): no standard errors for this fit (LASSO or itemp=1)'


3 format(a,i6,3(2x,g14.7))
4 format(i6,3(2x,g14.7))
5 format(a,99(2x,g10.3))
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
 real(r15) dij(natom_prim_cell,mx),rr(3),eps(3),aux(mx),aux2(mesh)
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
 !     write(*,*)' i0,mx,shell_count= ',i0,mx,shel_count
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
    aux=dij(i0,:) 
    call calculate_dos(mx,aux,w0,mesh,dmesh,aux2)
    ds(i0,:)= aux2 
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
 integer rnk,t,ti,i,res,j,rs,k,taup
 integer iat(maxrank),ixyz(maxrank),g,ng,term,term2,cnt2,frm,cnt3,ntind(maxrank),ngroup(maxrank)
 real(r15) rij,bunit,one,fcd,trace_fc,dij(3) 
! character frmt*2,goh*48,ln*1,geh*47
 character frmt*2,lm*1,goh*60,ln*1,geh*60,goh2*71,geh2*71

 bunit = ryd/ab/ab
 one =1d0

6 format(2x,i5,1x,a2,2x,i5,9(2x,f19.10))
7 format(1(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
8 format(2(2x,i5),3(2x,f19.10),4(2x,i5),2x,f9.5)
9 format(9(2x,f19.10))

! first write the crystal data
 ngroup= map(:)%ntot
 ntind=map(:)%nkeptind
 write(ulog,*)'WRITE_OUTPUT_FCS: ntind=',ntind
 call write_lat_fc(ntind,ngroup,'lat_fc.dat') !map(:)%ntotind,map(:)%ntot,'lat_fc.dat') 

!----------------------------------------
 res = 0
 ranks: do rnk=1,maxrank
  if ( include_fc(rnk) .eq. 1 ) then
    frm=30+rnk
    write(ln,'(i1)')rnk
    write(lm,'(i1)')rnk-1
    goh='(a1,i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g15.8,1x,f7.4,1x,15i3)'
    geh='(   i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g15.8,1x,f7.4,1x,15i3)'
    goh2='(a1,i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g15.8,1x,f7.4,1x,4i3,6(1x,f7.4))'
    geh2='(   i6,1x,i5,'//ln//'(3x,i4,1x,i1),3x,g15.8,1x,f7.4,1x,4i3,6(1x,f7.4))'
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
          if(g.gt.1) cnt2=sum(map(2)%keep(1:counteri(2,g-1,map(2)%ntind(g-1))))
       else
          if(g.gt.1) cnt2 = counteri(rnk,g,1)-1 !sum(map(rnk)%ntind(1:g-1)) ! cumulative indep term index up to, but excluding  group g (used for ranks other than 2)
       endif
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
!         write(ulog,goh) map(rnk)%err(cnt2+ti),g,ti,(iat(j),ixyz(j),j=1,rnk),  &
!      &     fcs(res+cnt2+ti),fcs(res+cnt2+ti)/ryd*ab**rnk,rij
!         write(ufit1-1+rnk,geh) ti,g,(iat(j),ixyz(j),j=1,rnk),  &
!      &     fcs(res+cnt2+ti),one,rij
          if(rnk.eq.2) then 
             if(map(2)%keep(counteri(2,g,ti)).eq.0) cycle
             write(ulog,goh2) map(rnk)%err(current(2,g,ti)),g,ti,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+current(2,g,ti)),rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk), &
       &     matmul(cart_to_conv,atompos(:,iat(1))-atompos(:,iat(2)))
             write(ufit1-1+rnk,geh2) ti,g,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+current(2,g,ti)),rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk), &
       &     matmul(cart_to_conv,atompos(:,iat(1))-atompos(:,iat(2)))
          else
             write(ulog,goh) map(rnk)%err(cnt2+ti),g,ti,(iat(j),ixyz(j),j=1,rnk),  &
       &     fcs(res+current(rnk,g,ti)),rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk)
!      &     fcs(res+cnt2+ti),rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk)
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
       ! must find the corresponding index of the indep term t <-> ti
          do ti=1,map(rnk)%ntind(g)
         ! this is the index of the indep FC coming in the A*FC=b matrix product
             if(rnk.eq.2) then
                if(map(2)%keep(counteri(2,g,ti)).eq.0) cycle
                fcd = fcd + fcs(res+current(2,g,ti))*map(rnk)%gr(g)%mat(t,ti) 
             else
                fcd = fcd + fcs(res+cnt2+ti)*map(rnk)%gr(g)%mat(t,ti)
             endif
          enddo
          if(abs(fcd).lt.margin) cycle
          if(rnk.eq.2 ) then
             write(ufc1-1+rnk,geh2)t,g, (iat(j),ixyz(j),j=1,rnk),fcd,  &  !one
&              rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk),  &
       &     matmul(cart_to_conv,atompos(:,iat(1))-atompos(:,iat(2)))
          else
             write(ufc1-1+rnk,geh)t,g, (iat(j),ixyz(j),j=1,rnk),fcd,  &  !one
&              rij,(iatomcell0(iat(k)),iatomcell(:,iat(k)),k=2,rnk)
          endif
       enddo
  enddo groups

     res = res + map(rnk)%nkeptind
     write(ulog,*)'Rank=',rnk,' res= ',res

  endif
 enddo ranks

write(ulog,*)'******* Trace for the harmonic FCs ********'
 open(456,file='trace_fc.dat')
 write(456,*)'# tau_i,tau_j,j,dij,trace_fc,atompos(:,i),atompos(:,j),rj-ri'
 open(457,file='bond_fci.dat')
 write(457,*)'# g,iat(1),iat(2), dij, ri, rj-ri '

! write the trace of FC2
 rnk=2; res=map(1)%nkeptind
 iloop: do i=1,natom_prim_cell
     do taup=1,natom_prim_cell
 jloop: do j=1,natoms
     if(iatomcell0(j).ne.taup) cycle
     rij = length(atompos(:,j)-atompos(:,i))
     trace_fc=0
     rs=map(1)%nkeptind
     if ( include_fc(rnk) .ne. 0 ) then
        ng=map(rnk)%ngr ! number of groups
        do g=1,map(rnk)%ngr  ! index of a given group
           do t=1,map(rnk)%nt(g)  ! index of independent terms in that group g
              iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
              ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
              if (iat(1).ne.i .or. iat(2).ne.j) cycle !iloop
              fcd = 0
       ! must find the corresponding index of the indep term t <-> ti
              do ti=1,map(rnk)%ntind(g)
                 fcd = fcd + fcs(res+current(2,g,ti)) * map(rnk)%gr(g)%mat(t,ti) 
              enddo
              if (ixyz(1).eq.ixyz(2)) then
                 trace_fc = trace_fc + fcd
              endif
           enddo

           do ti=1,map(rnk)%ntind(g)
              iat(1:rnk)  = map(rnk)%gr(g)%iat (:,ti)
              ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,ti)
              if (iat(1).ne.i .or. iat(2).ne.j) cycle 
              if(map(2)%keep(counteri(2,g,ti)).eq.1 .and.ixyz(1).eq.1.and.ixyz(2).eq.1) then
                dij = atompos(:,iat(2))-atompos(:,iat(1))
                write(457,3) g,iat(1),iat(2),length(dij),atompos(:,iat(1)),dij
              endif
           enddo

        enddo
! write(ulog,4)'i0,j,g,t,fcd,tr=',i,j,g,t,fcd,trace_fc
        if(abs(trace_fc).gt.margin) then  ! .and. ixyz(1).eq.1 .and. ixyz(2).eq.1) then
           write(ulog,8) i,j,rij,trace_fc
!          write(456,3) i,iatomcell0(j),j,rij,trace_fc,atompos(:,i),atompos(:,j), &
           write(456,3) i,taup,j,rij,trace_fc,atompos(:,i),atompos(:,j), &
&                       cart2red(atompos(:,j)-atompos(:,i),'r')
        endif
     endif
  enddo jloop
    enddo
  enddo iloop

  close(456)
  close(457)

  write(ulog,*)'***************** END OF FC Trace ******************'

2 format(9(1x,f10.4))
3 format(3i5,99(1x,f10.4))
4 format(a,4i4,9(1x,f10.4))


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
 integer rnk,t,ti,i,j,k
 integer iat(maxrank),ixyz(maxrank),g,term,cnt2,ntind(maxrank),ngroup(maxrank)
 real(r15) bunit,one,fcd,trace_fc,dij , zero(3)
 character lm*1

 one =1d0
 zero = 0

 rnk=2
 if ( include_fc(rnk) .ne. 0 ) then
 iloop: do i=1,natom_prim_cell
    write(lm,'(i1)')i
    open(875+i,file='springs'//lm//'.dat')
    write(875+i,*)"# 0 0 0 ,atompos(:,iat(1))-atompos(:,iat(2)) , d12 "
! jloop: do j=1,natoms
        do g=1,map(rnk)%ngr  ! index of a given group
           do t=1,map(rnk)%nt(g)  ! index of independent terms in that group g
              iat(1:rnk)  = map(rnk)%gr(g)%iat (:,t)
              ixyz(1:rnk) = map(rnk)%gr(g)%ixyz(:,t)
              if (iat(1).ne.i ) cycle  !.or. iat(2).ne.j) cycle 
              dij = length(atompos(:,iat(1))-atompos(:,iat(2)))
              do ti=1,map(rnk)%ntind(g)
                 if( map(2)%keep(counteri(2,g,ti) ) .eq. 0 ) cycle 
                 if(ixyz(1).eq.1 .and. ixyz(2).eq.1) then
                     write(875+i,2)zero ,atompos(:,iat(1))-atompos(:,iat(2)),dij
                 endif
              enddo
           enddo
        enddo
!  enddo jloop
     close(875+i)
  enddo iloop
  endif

2 format(9(1x,f10.4))
 end subroutine write_springs
!============================================================
 subroutine read_fcs
!! include_fc(rank)=2: reads fcN.dat of those ranks and lat_fc.dat of the run that wrote them, matches every term to
!! this run's atoms through (sublattice, cell) (the two runs' atom numbering need not agree) and sets the irreducible
!! values xfull (module given_fcs) by least squares per group. Replaces the earlier read_fcs/read_fcs_3/read_fcs_old
!! (2026-11-03, user).
 use given_fcs, only : xfull, nfull, full_col, solve_small, sort_i8, find_i8
 use svd_stuff
 use params
 use ios, only : ulog
 use atoms_force_constants
 use lattice, only : r01, r02, r03
 use geometry, only : v2a
 implicit none
 integer, parameter :: ulat=191, ufc=192
 integer nat0, i, j, k, r, g, t, ti, nold, ia(8), ix(8), ios_, nlo(3), nhi(3), tau, nn(3), ncur, nread, nmiss, &
 &       nout, nt, nk, ib, jb, iold, cnt
 integer, allocatable :: otau(:), on(:,:), lookup(:,:,:,:), kidx(:), korder(:), tg(:), tt(:)
 integer(8), allocatable :: keys(:,:)
 integer(8) kk(2)
 real(r15) v(3), rold(3,3), val, resid, resmax, fmax, posold(3)
 real(r15), allocatable :: f(:,:), a(:,:), x(:), ata(:,:), atb(:)
 character(len=300) line
 character(len=1) rn
 logical ex
 integer, allocatable :: kept(:)

 nfull = 0
 do r = 1, maxrank
    if (include_fc(r) .ne. 0) nfull = nfull + map(r)%nkeptind
 enddo
 if (allocated(xfull)) deallocate(xfull)
 allocate(xfull(nfull)); xfull = 0
 do r = 5, maxrank
    if (include_fc(r) .eq. 2) then
       write(*,*)'READ_FCS: include_fc=2 is implemented for ranks 1-4 only; rank ',r
       write(ulog,*)'READ_FCS: include_fc=2 is implemented for ranks 1-4 only; rank ',r
       stop
    endif
 enddo

! ---- lat_fc.dat of the run that produced the FCs: cell check and atom number -> (tau, n)
 inquire(file='lat_fc.dat', exist=ex)
 if (.not. ex) then
    write(*,*)'READ_FCS: include_fc=2 needs lat_fc.dat (and fcN.dat) of the run that produced the FCs'
    write(ulog,*)'READ_FCS: include_fc=2 needs lat_fc.dat (and fcN.dat) of the run that produced the FCs'
    stop
 endif
 open(ulat, file='lat_fc.dat', status='old')
 read(ulat,'(a)') line
 do i = 1, 3
    read(ulat,*) rold(:,i)
 enddo
 if (maxval(abs(rold(:,1)-v2a(r01))) + maxval(abs(rold(:,2)-v2a(r02))) + maxval(abs(rold(:,3)-v2a(r03))) .gt. 1d-4) then
    write(*,*)'READ_FCS: the primitive vectors of lat_fc.dat differ from structure.params'
    write(ulog,*)'READ_FCS: the primitive vectors of lat_fc.dat differ from structure.params'
    stop
 endif
 read(ulat,'(a)') line
 read(ulat,*) nat0
 if (nat0 .ne. natom_prim_cell) then
    write(*,*)'READ_FCS: lat_fc.dat has ',nat0,' atoms in the primitive cell, structure.params ',natom_prim_cell
    stop
 endif
 do i = 1, nat0
    read(ulat,'(a)') line
 enddo
 do
    read(ulat,'(a)',iostat=ios_) line
    if (ios_ .ne. 0) then
       write(*,*)'READ_FCS: no atom list (maxshells,Neighborshell atoms) in lat_fc.dat'
       stop
    endif
    if (index(line,'Neighborshell atoms') .gt. 0) exit
 enddo
 read(ulat,*) i, nold, j
 allocate(otau(nold), on(3,nold))
 do i = 1, nold
    read(ulat,*) iold, posold, otau(i), on(:,i)
    if (iold .ne. i) then
       write(*,*)'READ_FCS: atom list of lat_fc.dat out of order at ',i
       stop
    endif
 enddo
 close(ulat)
 write(ulog,*)'READ_FCS: lat_fc.dat read, cell matches, atoms listed: ',nold

! ---- this run's atoms by (tau, n)
 nlo = minval(iatomcell(:,1:natoms),2); nhi = maxval(iatomcell(:,1:natoms),2)
 allocate(lookup(natom_prim_cell, nlo(1):nhi(1), nlo(2):nhi(2), nlo(3):nhi(3))); lookup = 0
 do i = 1, natoms
    lookup(iatomcell0(i), iatomcell(1,i), iatomcell(2,i), iatomcell(3,i)) = i
 enddo
 if (natoms .ge. 2**20) then
    write(*,*)'READ_FCS: more than 2**20-1 atoms in the neighbour list; term keys overflow'
    stop
 endif

 do r = 1, 4
    if (include_fc(r) .ne. 2) cycle
    write(rn,'(i1)') r
! keys of this run's terms: (iat(1:r), ixyz(1:r)) -> (g, t)
    ncur = sum(map(r)%nt(1:map(r)%ngr))
    allocate(keys(2,ncur), tg(ncur), tt(ncur), korder(ncur))
    k = 0
    do g = 1, map(r)%ngr
       do t = 1, map(r)%nt(g)
          k = k + 1
          keys(:,k) = term_key(r, map(r)%gr(g)%iat(1:r,t), map(r)%gr(g)%ixyz(1:r,t))
          tg(k) = g; tt(k) = t
       enddo
    enddo
    call sort_i8(ncur, keys, korder)
    allocate(f(maxval(map(r)%nt(1:map(r)%ngr)), map(r)%ngr)); f = 0
    inquire(file='fc'//rn//'.dat', exist=ex)
    if (.not. ex) then
       write(*,*)'READ_FCS: include_fc(',r,')=2 but fc'//rn//'.dat is missing'
       write(ulog,*)'READ_FCS: include_fc(',r,')=2 but fc'//rn//'.dat is missing'
       stop
    endif
    open(ufc, file='fc'//rn//'.dat', status='old')
    nread = 0; nmiss = 0; nout = 0
    do
       read(ufc,'(a)',iostat=ios_) line
       if (ios_ .ne. 0) exit
       if (len_trim(line) .eq. 0) cycle
       if (index(adjustl(line),'#') .eq. 1) cycle
       read(line,*,iostat=ios_) t, g, (ia(j), ix(j), j=1,r), val
       if (ios_ .ne. 0) cycle
       nread = nread + 1
       do j = 1, r                                   ! old atom number -> this run's atom
          if (ia(j) .lt. 1 .or. ia(j) .gt. nold) then
             write(*,*)'READ_FCS: fc'//rn//'.dat names atom ',ia(j),' not in lat_fc.dat (1..',nold,')'
             stop
          endif
          tau = otau(ia(j)); nn = on(:,ia(j))
          if (any(nn .lt. nlo) .or. any(nn .gt. nhi)) then
             ia(j) = 0
          else
             ia(j) = lookup(tau, nn(1), nn(2), nn(3))
          endif
       enddo
       if (any(ia(1:r) .eq. 0)) then
          nout = nout + 1; cycle
       endif
       kk = term_key(r, ia(1:r), ix(1:r))
       i = find_i8(ncur, keys, korder, kk)
       if (i .eq. 0) then
          nmiss = nmiss + 1; cycle
       endif
       f(tt(i), tg(i)) = val
    enddo
    close(ufc)
! irreducible values group by group (least squares; kept terms only for rank 2)
    resmax = 0; fmax = 0
    do g = 1, map(r)%ngr
       nt = map(r)%nt(g)
       allocate(kept(map(r)%ntind(g))); nk = 0
       do ti = 1, map(r)%ntind(g)
          if (full_col(r,g,ti) .gt. 0) then
             nk = nk + 1; kept(nk) = ti
          endif
       enddo
       if (nk .gt. 0) then
          allocate(a(nt,nk), x(nk), ata(nk,nk), atb(nk))
          do jb = 1, nk
             a(:,jb) = map(r)%gr(g)%mat(1:nt, kept(jb))
          enddo
          ata = matmul(transpose(a), a); atb = matmul(transpose(a), f(1:nt,g))
          call solve_small(nk, ata, atb, x)
          do jb = 1, nk
             xfull(full_col(r,g,kept(jb))) = x(jb)
          enddo
          resid = maxval(abs(matmul(a,x) - f(1:nt,g)))
          resmax = max(resmax, resid); fmax = max(fmax, maxval(abs(f(1:nt,g))))
          deallocate(a, x, ata, atb)
       endif
       deallocate(kept)
    enddo
    write(ulog,'(a,i2,a,i8,a,i8,a,i8,a,i8)')' READ_FCS: rank',r,': terms read ',nread,', matched ', &
 &        nread-nmiss-nout,', not in this run''s groups ',nmiss,', atoms outside its neighbour list ',nout
    write(ulog,'(a,i2,a,es10.2,a,es10.2)')' READ_FCS: rank',r,': max |mat x - f| over the terms =',resmax, &
 &        ' (max |f| =',fmax
    write(*,'(a,i2,a,i8,a,es10.2)')' READ_FCS: rank',r,' read from fc'//rn//'.dat: terms ',nread, &
 &        ', max residual of the irreducible fit',resmax
    if (nmiss+nout .gt. 0) then
       write(*,*)'READ_FCS: WARNING: ',nmiss+nout,' terms of fc'//rn//'.dat are not FCs of this run (ignored)'
       write(ulog,*)'READ_FCS: WARNING: ',nmiss+nout,' terms of fc'//rn//'.dat are not FCs of this run (ignored)'
    endif
    if (resmax .gt. 1d-4*max(fmax,1d-30)) then
       write(*,*)'READ_FCS: WARNING: the read FCs of rank ',r,' do not follow this run''s symmetry relations'
       write(ulog,*)'READ_FCS: WARNING: the read FCs of rank ',r,' do not follow this run''s symmetry relations'
    endif
    deallocate(keys, tg, tt, korder, f)
 enddo
 deallocate(otau, on, lookup)
 cnt = 0
 do r = 1, maxrank
    if (include_fc(r) .eq. 2) cnt = cnt + map(r)%nkeptind
 enddo
 write(ulog,*)'READ_FCS: irreducible FCs read (all read ranks): ',cnt,' of ',nfull

 contains
   function term_key(rr, at, dx) result(kk2)
!! (first atom and the directions ; the other atoms, 20 bits each): two integers, compared in that order
   integer, intent(in) :: rr, at(:), dx(:)
   integer(8) kk2(2)
   integer m
   kk2(1) = int(at(1),8)
   do m = 1, rr
      kk2(1) = kk2(1)*4_8 + int(dx(m)-1,8)
   enddo
   kk2(2) = 0
   do m = 2, rr
      kk2(2) = kk2(2)*1048576_8 + int(at(m),8)
   enddo
   end function term_key
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
! bring between -0.5 and 0.5
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
 subroutine warn3(unt,mesg,output)
 use constants, only: r15
 implicit none
 integer, intent(in) :: unt
 real(r15), intent(in) :: output
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
 write(unt,*)'|   '//mesg,output
 end subroutine warn3
!============================================================
 subroutine warn4(unt,mesg,output)
 use constants, only: r15
 implicit none
 integer, intent(in) :: unt
 real(r15), intent(in) :: output(:)
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
 write(unt,3)'|   '//mesg,output
3 format(a,999(1x,g11.4))
 end subroutine warn4
!============================================================
 subroutine warn5(unt,mesg,output)
 use constants, only: r15
 implicit none
 integer, intent(in) :: unt
 real(r15), intent(in) :: output(:,:)
 character(len=*), intent(in) :: mesg
 integer i

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
 do i=1,size(output,1)
    write(unt,3)'| '//mesg,output(i,:)
 enddo
3 format(a,999(1x,g11.4))
 end subroutine warn5
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
 write(ufco,*)' Crystal data: translation vectors of the primitive cell ', &
 &   '  (columns 4-6: the conventional cell vectors, Cartesian; older readers take columns 1-3)'
 write(ufco,9)r01,r1conv
 write(ufco,9)r02,r2conv
 write(ufco,9)r03,r3conv
 write(ufco,*)' Crystal data: atoms in primitive cell: label,type,x,y,z,mass '
 write(ufco,*)natom_prim_cell
 do i=1,natom_prim_cell
   write(ufco,6)i,atom0(i)%name, atom0(i)%at_type,atom0(i)%equilibrium_pos,atom0(i)%mass
 enddo

 write(ufco,*)' Included ranks of FCs '
 write(ufco,*)  include_fc(:)
 write(ufco,*)' Number of shells for each atom in each rank '
 do i=2,maxrank
    write(ufco,*)nshells(i,:)
 enddo
 write(ufco,*)' Crystal data: Number of groups for each rank ' 
 write(ufco,*)map(1:8)%ngr

 write(ufco,*)' Number of FCs for each rank '
 write(ufco,*)  ntrms(:)  ! same as  write(ufco,*)map(1:8)%ntot

 write(ufco,*)' Number of independent FCs for each rank '
 write(ufco,3)  map(1:8)%ntotind,map(1:8)%nkeptind  ! same as  write(ufco,*) ntindep(:) 
3 format(8i6,9x,8i6)

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


 uio=931
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
! read(uio,4)wrd1,map(1:8)%ngr,wrd2,map(1:8)%ntotind,wrd3,map(1:8)%ntot
 read(uio,4)wrd1,map(1:8)%ngr,wrd2,map(1:8)%nkeptind,wrd3,map(1:8)%ntot

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
 use svd_stuff, only : atransl,arot,brot,ahuang,bhuang, amat,bmat,transl_constraints, &
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
      write(ulog,*)i,dot_product(ahuang(i,:),fcs(:))-bhuang(i),bhuang(i)
   enddo

! not used since amat is overwritten after the SVD
!  write(umatrx,*)'******* Violation of force-displacement relations: '
!  errmax=-1
!  err=0;num=0;denom=0
!  do i=inv_constraints+1,dim_al
!     prod = dot_product(amat(i,:),fcs(:))
!     junk = abs(prod-bmat(i))
!     write(umatrx,3)' ',i,prod,bmat(i),junk
!     if (errmax .lt. junk) errmax=junk
!     err=err+junk
!     num=num+junk*junk
!     denom=denom+prod*prod
!  enddo
!  err=err/(dim_al-inv_constraints)
!  write(umatrx,2)'Max and average errors in force-displacements,percent deviation=', &
!  &    errmax,err,sqrt(num/denom)*100
!  write(ulog,2)'Max and average errors in force-displacements,percent deviation=', &
!  &    errmax,err,sqrt(num/denom)*100

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
       tloop: do t=1,map(2)%nt(g)
          j  = map(2)%gr(g)%iat(2,t)
          if(j.gt. largest) largest=j
       enddo tloop
    enddo gloop
    write(*   ,*)' The largest atom number found in atompos was ',largest
    write(ulog,*)' The largest atom number found in atompos was ',largest

 end subroutine find_largest_atom_number
!============================================================
 subroutine write_forces_display(nat,ncfg,dsp,frc,fn)
 use constants, only : r15
 use atoms_force_constants, only : atom_sc
 use geometry, only : v2a
 implicit none
 integer, intent(in) :: ncfg, nat
 real(r15), intent(in) :: dsp(3,nat,ncfg),frc(3,nat,ncfg)
 integer t,j
 character,intent(in) :: fn*(*)
 character num*3
    
 do t=1,ncfg
    write(num,'(i3.3)')t
    open(765,file=fn//num//'.dat')
    do j=1,nat
       write(765,4)atom_sc(j)%equilibrium_pos, dsp(:,j,t), frc(:,j,t)
    enddo
    close(765)
 enddo
4 format(9(1x,f10.5))

 end subroutine write_forces_display
!===================================================
 subroutine write_lattice(n,lat,fn)
!! writes the lattice vectors lat(3,n) in file fn in the xyz format
 use lattice
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(in) :: lat(3,n)
 character(*), intent(in) :: fn
 integer i
 real(r15) r2pi

 r2pi=(length(rs1)*length(rs2)*length(rs3))**0.67 /3.14 ! to rescale gs so that they are visible on the scale of rs 

 open(125,file=fn)
 write(125,*)n+15
 write(125,*)'#supercell lattice; rescaled G vectors followed by R-supercell vectors'
 write(125,*)'Ge   0 0 0 '
 write(125,4)'Ge ',r2pi*gs1 
 write(125,4)'Ge ',r2pi*gs2
 write(125,4)'Ge ',r2pi*gs3
 write(125,4)'Ge ',r2pi*(gs1+gs2)
 write(125,4)'Ge ',r2pi*(gs3+gs1)
 write(125,4)'Ge ',r2pi*(gs3+gs2)
 write(125,4)'Ge ',r2pi*(gs3+gs1+gs2)
 write(125,4)'Bi ',rs1
 write(125,4)'Bi ',rs2
 write(125,4)'Bi ',rs3
 write(125,4)'Bi ',rs1+rs2
 write(125,4)'Bi ',rs2+rs3
 write(125,4)'Bi ',rs3+rs1
 write(125,4)'Bi ',rs1+rs2+rs3
 do i=1,n
    write(125,4)'C  ',lat(:,i)
 enddo
 close(125)

4 format(a,9(1x,f11.5))

 end subroutine write_lattice
!===================================================

 subroutine write_fc2matrix_5d(nat,nr,rgrd,phi,fn)
!! takes the fitted FCs, writes output phi_bare into 5D matrix form
 use constants
 use atoms_force_constants
 use lattice, only : rws26 , fold_ws ,cart2red
 use geometry, only : length,trace
 use ios, only : ulog
 implicit none
 integer, intent(in) :: nat,nr
 real(r15), intent(in) :: rgrd(3,nr), phi(nat,nat,3,3,nr)
 character(*), intent(in) :: fn
 integer tau,taup,al,be,g,uio,igrid
 real(r15) dta(3),rg(3),fc2,asr,rfold(3) !,fold_ws


 uio=345
 open(uio,file=fn)
 write(uio,*)'# tau,taup,j,|r+taup-tau|,trace(fc2),igrid,red_fold(r+taup-tau),fc2(al,be)'
 do tau =1,nat
 do taup=1,nat
    dta = atompos(:,taup) - atompos(:,tau) 
    do g=1,nr
       rg=rgrd(:,g)+dta 
       rfold = fold_ws(rg,rws26) 
!      if(length(rg-rfold).gt.1d-3) then
!         write(ulog,3)'write_fc2: t,tp,ir,r+tp-t, folded =',tau,taup,g,rg,rfold
!      endif
       fc2 = trace(phi(tau,taup,:,:,g)) 
       if(abs(fc2).lt.0.00001) cycle
!      write(uio,3)tau,taup,g,length(rg   ),fc2,cart2red(rfold,'r'),real(phi(tau,taup,:,:,g)) 
       write(uio,3)' ',tau,taup,g,length(rfold),fc2,rfold,real(phi(tau,taup,:,:,g)) 
    enddo
 enddo
 enddo
 close(uio)

 write(ulog,*)'ASR CHECK for fc2 matrix: tau,al,be,asr going to '//fn
 do tau =1,nat
 do al=1,3
 do be=1,3
    asr=sum(phi(tau,:,al,be,:)) ! this should satisfy ASR!
    write(ulog,3)' ',tau,al,be,asr
 enddo
 enddo
 enddo

! check whether phi(tau,taup,al,be,r)=phi(taup,tau,be,al,-r) , assuming same weight
 call check_herm_sym_R(nat,nr,rgrd,phi)

3 format(a,3i5,99(f9.4))

 end subroutine write_fc2matrix_5D
!==========================================================
 subroutine write_all_eigenvalues(nk,kp,dk,n,eival,vg,uio)
 use constants
 use eigen, only : mysqrt
 use geometry
 implicit none
 integer, intent(in) :: n,uio,nk
 real(r15), intent(in) :: dk(nk),kp(3,nk),eival(n,nk),vg(3,n,nk)
 integer j,k

 write(uio,*)'# nk,nb,dk(i),kp(:,i), freq(nb,nk) ; velocity(nb,nk,:) ; |velocity| '
 do k=1,nk
 do j=1,n
    write(uio,3)k,j,dk(k),kp(:,k),cnst*mysqrt(eival(j,k)),vg(:,j,k)/1000,length(vg(:,j,k))/1000
 enddo
 enddo

2 format(i5,1x,4(1x,f8.3),999(1x,g11.4))
3 format(2i5,1x,4(1x,f8.3),999(1x,g11.4))

 end subroutine write_all_eigenvalues
