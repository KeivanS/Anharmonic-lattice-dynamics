!==========================================================
 subroutine read_params
 use ios
 use om_dos
 use params
 use lattice
 use phi3
 use kpoints
 use atoms_force_constants
 implicit none
! integer i

  open(uparams,file='latdyn.params',status='old')
  write(6,*) 'READ_PARAMS: opening latdyn.params'
  read(uparams,*)nc,dos_bs_only ! kmesh 
  write(*,*) 'READ_PARAMS: kmesh= ',nc
  read(uparams,*)shft !x,shfty,shftz ! shift in units of mesh
  read(uparams,*)wmesh,wmax  ! no of w mesh points and max frequency for DOS
  write(*,*) 'READ_PARAMS: just read wmesh,wmax ',wmesh,wmax
  read(uparams,*)width,etaz  ! width of gaussian broadening for DOS,imag part of om
  write(*,*) 'READ_PARAMS: just read width,etaz ',width,etaz
  tau0=1e-11
  write(*,*) 'READ_PARAMS: default tau0 set to ',tau0
  read(uparams,*)verbose
  write(*,*) 'READ_PARAMS: just read verbose ',verbose
  read(uparams,*)tmin,tmax,ntemp      ! temps in Kelvin to calculate thermal expansion and other thermal ppties
  write(*,*) 'READ_PARAMS: tmin,tmax,ntemp=',tmin,tmax,ntemp
  read(uparams,*)job      ! if=1 then read from a file, else generate  ! sy added iter,split,calk
  write(*,*)'iter,split,readv3,writev3,calk=',iter,split,readv3,writev3,calc_kappa
  read(uparams,*)ksub_size  ! each v33sq.xxx.dat file contains V33sq(ksub_size,nkc,ndn,ndn,ndn)
  read(uparams,'(a)') v3path              ! path to v33sq.dat files
  read(uparams,*)iso
! initialize some default parameters
 tolerance =0.0001
  max_iter=150
  conv_error=1E-5
  conv_max_error=1E-4
  conv_diff=1E-8
  conv_max_diff=1E-6
  conv_diff_kap=1E-5
  conv_max_diff_kap=1E-4
  Conv_iter=1E-6
  update_mix=0.4  ! mixing of new iteration
  write(*,*)'max_iter, conv_error, conv_max_error, conv_diff, conv_max_diff, conv_diff_kap, conv_max_diff_kap,conv_iter,update_mix'   
  write(*,*)max_iter, conv_error, conv_max_error, conv_diff, conv_max_diff, conv_diff_kap, conv_max_diff_kap,conv_iter,update_mix   
  v3_threshold=1E-7        ! in 1/cm keep all v33 of norm above v3_threshold
  write(*,*)'v3_threshold=',  v3_threshold  
  read(uparams,*)classical     ! if =1 then use classical distribution (kT/hw)
  read(uparams,*)calc_cross,q_cross  ! if =1 then calculate the cross section at qcros (reduced U)
  read(uparams,*)lmicron     ! sample length in microns

!!mine

  close(uparams)

  write(ulog,*) 'READ_PARAMS: read and kpoints allocated'
  write(ulog,3)'nc1,nc2,nc3=',nc
  write(ulog,*)'wmesh, wmax=',wmesh,wmax
  write(ulog,*)'read/writev3=',readv3,writev3
! etaz = max(etaz,wmax/wmesh)
! write(ulog,*)' in case etaz is too small we use the following'
! write(ulog,*)'width,etaz =',width,etaz
  write(ulog,*)'Tmin,Tmax(K=',tmin,tmax
  if (classical .eq. 1) then
     write(ulog,*)'Calculation is done in the CLASSICAL limit'
  else
     write(ulog,*)'Calculation is done in the QUANTUM limit'
  endif
  if (dos_bs_only) then
     write(ulog,*)' DOS will only be calculated (on the coarse mesh) ',nc
!    write(ulog,*)' and the program will stop after that!'
  endif

3 format(a,6(1x,i6))
 end subroutine read_params
!==========================================================
 subroutine read_input_fit
 use ios
 use params
 use lattice
 use atoms_force_constants
 implicit none
 integer i,counter,label
 real(8) scal,junk
 character jjj*1
 open(uparams,file='structure.params',status='old')

 read(uparams,*) latticeparameters   ! (a,b,c,alpha,beta,gamma)
 read(uparams,*) primitivelattice     ! prim latt vectors in terms of conventional above
 read(uparams,*) lattice_parameter
 latticeparameters(1:3) = latticeparameters(1:3)*lattice_parameter
 read(uparams,*) include_fc    ! if=1 include this rank
 read(uparams,*) junk ! itrans,irot,ihunag,all
 read(uparams,*) junk ! use temperature weighting, Tempk
 read(uparams,*) junk ! # of FORCEDISP files, verbosity
 read(uparams,*) natom_type   ! # of different elements present in prim cell
 write(ulog,*) 'READ_INPUT_FIT: natom_type ',natom_type
 call allocate_mass(natom_type) ! and atname and natom
 read(uparams,*) (mas(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) (atname(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) natoms0        ! # of atoms in primitive cell
 write(ulog,*)'reading ',natoms0,' atoms in the primitive cell'
 call allocate_primcell(natoms0)
 read(uparams,*) junk ! fc2flag
 read(uparams,*) junk ! nshells(2)
 read(uparams,*) junk ! nshells(3)
 read(uparams,*) junk ! nshells(4)
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
 write(ulog,*)'reading  atom #',i, counter
 enddo

 write(ulog,*)'reading done, closing the structure.params file'
 close(uparams)

 maxneighbors=50

 do i=1,natoms0
    atom0(i)%name = atname(atom_type(i))
    atom0(i)%at_type  = atom_type(i)
    atom0(i)%mass = mas(atom_type(i))
    atom0(i)%tau  = i  ! assumes labels are ordered in increasing order
    atompos0(:,i) = atom0(i)%equilibrium_pos
 enddo

! do i=1,4
!    if (nshells(i) .gt. maxshells) then
!       write(ulog,*)' READ_INPUT: nshell> maxshells ',i,nshells(i),maxshells
!       write(ulog,*)' Either increase maxshells or lower nshell for that rank'
!       stop
!    endif
! enddo

end subroutine read_input_fit
!==========================================================
 subroutine read_lattice
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 implicit none
 character line*90,name*2
 integer i,j,tt,ttyp,natoms2,n1  !,at_count,n2,n3
 real(8) mss !a1,a2,a3,om,a,b,c,
 type(vector) tau1,vv

open(ufc0,file='lat_fc.dat' ,status='old')

do while (line(1:22) .ne. '  Crystal data: transl')
 read(ufc0,'(a)') line
enddo
 read(ufc0,*) r1   ! coordinates of the primitive cell vectors
 read(ufc0,*) r2
 read(ufc0,*) r3
 r01=r1; r02=r2; r03=r3

! now generate the lattice ppties in real and recip spaces
  call make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
  call calculate_volume(r1,r2,r3,volume_r)
  call calculate_volume(g1,g2,g3,volume_g)
 g01=g1; g02=g2; g03=g3
  box(1) = length(r1)
  box(2) = length(r2)
  box(3) = length(r3)
  boxg(1) = length(g1)
  boxg(2) = length(g2)
  boxg(3) = length(g3)
  write(ulog,3)' r1= ',r1
  write(ulog,3)' r2= ',r2
  write(ulog,3)' r3= ',r3
  write(ulog,3)' box  = ',box
  write(ulog,3)' g1= ',g1
  write(ulog,3)' g2= ',g2
  write(ulog,3)' g3= ',g3
  write(ulog,3)' boxg = ',boxg
  write(ulog,3)' volume_r,volume_g = ',volume_r,volume_g

 read(ufc0,*) line
! read(ufc,*) natoms0
! ############ check compatibiliy with input file "params.inp" ###########
 read(ufc0,*) n1
 if (n1 .ne. natoms0) then
    write(ulog,*)' natoms0 in lat_fc.dat incompatible with structure.params ',n1,natoms0
    stop
 endif
! call allocate_primcell(natoms0)
 do i=1,natoms0
!  read(ufc1,*)atom0(i)%tau, atom_type(i),atom0(i)%equilibrium_pos,atom0(i)%mass
    read(ufc0,*)tt,name,ttyp,tau1,mss  ! atomic types names coordinates and and masses 
    vv%x= tau1.dot.g1 ; vv%y= tau1.dot.g2 ; vv%z= tau1.dot.g3
    vv = (1/2d0/pi)*vv  ! reduced coordinates of atoms in primcell
    if(tt .eq. i .and. ttyp .eq. atom0(i)%at_type .and. atom0(i)%mass .eq. mss) then
      write(ulog,*)'####### compatibility with structure.params checked  ######### ',i
    else
      write(ulog,*)'READ_LATTICE: structure.params and lat_fc.dat are not compatible'
      write(ulog,*)' i, atom type, mass '
      write(ulog,*) i,tt,  atom0(i)%at_type,ttyp, atom0(i)%mass,mss
! fix it so that if eq pos are the same module translarion vectors, it's OK
      write(ulog,3)' reduced coordinates ',atom0(i)%equilibrium_pos,vv
      stop
    endif
 enddo


 read(ufc0,*) line
 read(ufc0,*) line
 read(ufc0,*) (include_fc(j),j=1,4)
 write(ulog,*) (include_fc(j),j=1,4)
 read(ufc0,*) line
 write(ulog,*) line
 read(ufc0,*) (nterms(j),j=1,4)
 write(ulog,*) (nterms(j),j=1,4)
 read(ufc0,*) line
 write(ulog,*) line
 read(ufc0,*) (ngroups(j),j=1,4)
 write(ulog,*) (ngroups(j),j=1,4)
 read(ufc0,*) line
 write(ulog,*) line
 read(ufc0,*) nsmax,natoms2
 write(ulog,*) nsmax,natoms2

  maxatoms=natoms2-300 
  imaxat=1
  do while (imaxat.ne.0)
     maxatoms=maxatoms+300
     write(6,*)' maxatoms=',maxatoms
     write(ulog,*)' maxatoms=',maxatoms
     if (allocated(iatomcell0)) deallocate(iatomcell0)
     if (allocated(iatomcell))  deallocate(iatomcell)
     if (allocated(atompos))    deallocate(atompos)
     allocate(atompos(3,maxatoms), iatomcell(3,maxatoms),iatomcell0(maxatoms))
! maxneighbors is an output of force_constants_init
!      call force_constants_init(lattparams,primlatt,natoms_in, &
!     &     iatomtype,atompos_in)
!     call force_constants_init(natoms0,atom_type,atompos0)
     call force_constants_init(latticeparameters,primitivelattice,natoms0,  &
      &     atom_type,atompos0)  !(n,natom_type,atompos0)
  enddo

  write(ulog,*)'After force_constants_init, maxneighbors=',maxneighbors

 if (natoms.ne.natoms2) then
    write(ulog,*) "natoms read in lat_fc=",natoms2
    write(ulog,*) "while natoms generated by fc_init=",natoms
    write(ulog,*) "check the number of shells in structure.params, "
    write(ulog,*) "it is probably too large if natoms in fc_init is larger than in lat_fc"
    write(ulog,*) "in which case the program may stop"
 endif

 
 do j=1,natoms2
    read(ufc0,*,end=99)i,atompos(:,i), iatomcell0(i),iatomcell(:,i)
    write(ulog,4)i,atompos(:,i), iatomcell0(i),iatomcell(:,i)
 enddo
3 format(a,9(2x,f11.6))
4 format(i6,3(2x,g11.5),4(3x,i4))

 close(ufc0)
 write(ulog,*)' COORDINATES read successfully '
 return

99 write(*,*) "End of lat_fc.dat file reached!"
   write(*,*) "check the number of shells in params.inp, it is probably too large!"
   stop

 end subroutine read_lattice
!===========================================================
 subroutine read_fc234
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 implicit none
 character line*90
 integer t,rank,res,j,mx
 logical ex

! read(ufc2,'(a)') line
! read(ufc,*)include_fc(1),include_fc(2),include_fc(3),include_fc(4)
! read(ufc,'(a)') line
! read(ufc,*)nterms(1),nterms(2),nterms(3),nterms(4)
! read(ufc,'(a)') line
! read(ufc,*)ngroups(1),ngroups(2),ngroups(3),ngroups(4)
 res = 0
!----------------------------------------
! rank=1
! if ( include_fc(rank) .eq. 1 ) then
! mx = nterms(rank)
! allocate(iatomterm_1(rank,mx),ixyzterm_1(rank,mx),igroup_1(mx), &
! & ampterm_1(mx),fcs_1(ngroups(rank)))
! read(ufc,*)line
! do j=1,nterms(rank)
!       read(ufc,*)t,igroup_1(t), &
!& iatomterm_1(1,t),ixyzterm_1(1,t),  &
!& fcs_1(igroup_1(t)),ampterm_1(t)
! enddo
! res = res + igroup_1(nterms(rank))
! endif
!----------------------------------------
 rank=2
 inquire ( file='fc2.dat', exist=ex)
 if (ex) then
    open(ufc2,file='fc2.dat',status='old')
 else
    write(*,*)'READ_FC234: file fc2.dat does not exist'
    stop
 endif
 read(ufc2,'(a)') line
 do j=1,10000000
       read(ufc2,*,end=92)t
 enddo
92 mx=j-1
 if (mx.ne.nterms(rank)) then
    write(*,*)'READ_FC2: number of lines read=',mx,' inconsistent with ',nterms(2),' from lat_fc.dat'
    stop
 endif
 rewind(ufc2)
 write(ulog,*)'Rank=2, reading nterms(2)=',mx,' terms'
 allocate(iatomterm_2(rank,mx),ixyzterm_2(rank,mx),igroup_2(mx)) 
 allocate(fcs_2(mx),grun_fc(mx))
 read(ufc2,'(a)') line
 write(*  ,'(a)') line
 write(*,*) '********** FCs for rank=2',mx

 do j=1,mx
       read(ufc2,*)t,igroup_2(j), &
& iatomterm_2(1,j),ixyzterm_2(1,j),  &
& iatomterm_2(2,j),ixyzterm_2(2,j),  &
& fcs_2(j)
 enddo
 res = res + igroup_2(nterms(rank))
! endif
!----------------------------------------
 rank=3
 if ( include_fc(rank) .eq. 1 ) then
   inquire ( file='fc3.dat', exist=ex)
   if (ex) then
      open(ufc3,file='fc3.dat',status='old')
   else
      write(*,*)'READ_FC234: file fc3.dat does not exist'
      stop
   endif
   read(ufc3,'(a)') line
   do j=1,10000000
       read(ufc3,*,end=93)t
   enddo
93 mx=j-1
   if (mx.ne.nterms(rank)) then
      write(*,*)'READ_FC3: number of lines read=',mx,' inconsistent with ',nterms(3),' from lat_fc.dat'
      stop
   endif
   write(*,*) '********** FCs for rank=3',mx
   rewind(ufc3)
   write(ulog,*)'Rank=3, reading nterms(3)=',mx,' terms'
   allocate(iatomterm_3(rank,mx),ixyzterm_3(rank,mx),igroup_3(mx), fcs_3(mx))
   read(ufc3,'(a)') line
   write(*  ,'(a)') line
   do j=1,mx
       read(ufc3,*)t,igroup_3(j), &
&   iatomterm_3(1,j),ixyzterm_3(1,j),  &
&   iatomterm_3(2,j),ixyzterm_3(2,j),  &
&   iatomterm_3(3,j),ixyzterm_3(3,j),  &
&   fcs_3(j)
!   write(*,*) j,fcs_3(j)
   enddo
   res = res + igroup_3(nterms(rank))
 endif
!--------------------------------
 rank=4
 if ( include_fc(rank) .eq. 1 ) then
 inquire ( file='fc4.dat', exist=ex)
 if (ex) then
    open(ufc4,file='fc4.dat',status='old')
 else
    write(*,*)'READ_FC4: file fc4.dat does not exist'
    stop
 endif
 read(ufc4,'(a)') line
 do j=1,10000000
       read(ufc4,*,end=94)t
 enddo
94 mx=j-1
 if (mx.ne.nterms(rank)) then
    write(*,*)'READ_FC4: number of lines read=',mx,' inconsistent with ',nterms(4),' from lat_fc.dat'
    stop
 endif
 rewind(ufc4)
 write(ulog,*)'Rank=4, reading nterms(4)=',mx,' terms'
 allocate(iatomterm_4(rank,mx),ixyzterm_4(rank,mx),igroup_4(mx), fcs_4(mx) )
 read(ufc4,'(a)') line
 write(*  ,'(a)') line
 write(*,*) '********** FCs for rank=4',mx

 do j=1,mx
       read(ufc4,*)t,igroup_4(j), &
& iatomterm_4(1,j),ixyzterm_4(1,j),  &
& iatomterm_4(2,j),ixyzterm_4(2,j),  &
& iatomterm_4(3,j),ixyzterm_4(3,j),  &
& iatomterm_4(4,j),ixyzterm_4(4,j),  &
& fcs_4(j)
 enddo
 res = res + igroup_4(nterms(rank))
 endif
 write(ulog,*)'READ_FC234: done!, number of groups read is=',res

 end subroutine read_fc234
!==========================================================
 subroutine set_dynamical_matrix(kpt,dynmat,ndim,ddyn)
 use params
 use lattice
 use kpoints
 use atoms_force_constants
 use geometry
 use svd_stuff
 implicit none
 integer, intent(in) :: ndim
 complex(8), intent(out) :: dynmat(ndim,ndim),ddyn(ndim,ndim,3)
 real(8), intent(in) :: kpt(3)
 complex(8) junk
 real(8) mi,mj,all,rr(3),delt(3)
 integer i0,j,j0,al,be,i3,j3,t !,ired

! write(*,*)'set dynamical matrix: fc2=',size(fcs_2)
! write(*,3)fcs_2
!
! write(*,*)'set dynamical matrix: iatomterm=',size(iatomterm_2) !(1,:))
! write(*,2)iatomterm_2(1,:)
! write(*,2)iatomterm_2(2,:)
!
! write(*,*)'set dynamical matrix: ixyzterm=',size(ixyzterm_2) !(1,:))
! write(*,2)ixyzterm_2(1,:)
! write(*,2)ixyzterm_2(2,:)

2 format(9999(1x,i4))
3 format(9999(1x,g9.3))
4 format(a,4(1x,i5),9(2x,f9.4))
 delt = wshift !*wshift
 ddyn   = cmplx(0d0,0d0)
 dynmat = cmplx(0d0,0d0)
 do i0=1,natoms0
 do al=1,3
    i3 = al+3*(i0-1)
    mi = atom0(i0)%mass
!  write(*,*) 'i,al,mass=',i0,al,i3,mi
    tloop: do t=1,nterms(2)
       if ( i0 .eq. iatomterm_2(1,t) .and. al .eq. ixyzterm_2(1,t) ) then
          be =  ixyzterm_2(2,t)
          j  = iatomterm_2(2,t)
          j0 = iatomcell0(j)
          mj = atom0(j0)%mass
          j3 = be+3*(j0-1)
          rr = atompos(:,j)-atompos(:,j0)
!          ired = igroup_2(t)
!          junk = fcs_2(ired)*ampterm_2(t)* exp(ci*(kpt .dot. rr))/sqrt(mi*mj)
          junk = fcs_2(t)* exp(ci*(kpt .dot. rr))/sqrt(mi*mj)
!  write(ulog,4) 't,j,be,mass,dyn=',t,j,be,j3,mj,junk
          dynmat(i3,j3) = dynmat(i3,j3) + junk
          ddyn(i3,j3,:) = ddyn(i3,j3,:) + junk*ci*rr(:)
!         if (be.eq.al .and. j0.eq.i0 ) then
!         if (be.eq.al .and. (length(rr-r1) .myeq. 0d0)) then
! same type and same cartesian coord but not the same atom, just take the images
!     write(ulog,5)'i0,j0,j,i3,j3,rij=',i0,j0,j,i3,j3,rr
!            dynmat(i3,j3) = dynmat(i3,j3) - delt/2
!         endif
       endif
    enddo tloop
! this leaves  gamma eivals unchanged
!   dynmat(i3,i3) = dynmat(i3,i3) + delt(al)*(1-cos(kpt .dot. r1))
! but this one shifts everything up by delt=wshift
!   dynmat(i3,i3) = dynmat(i3,i3) + delt(al)
 enddo
 enddo
! write(*,*) ' Terms loop done, dynmat is defined '

5 format(a,5i5,9(f8.3))
 all = sum(cdabs(dynmat(:,:)))/(ndim*ndim)
! make sure it is hermitian
 do t=1,ndim
!    call write_out(6,'dynmat=',dynmat(t,:))
    if (abs(aimag(dynmat(t,t))) .gt. 9d-4*abs(real(dynmat(t,t))) ) then
       write(ulog,*)' dynmat is not hermitian on its diagonal'
       write(ulog,*)' diagonal element i=',t,dynmat(t,t)
!      stop
!   else
       dynmat(t,t) = cmplx(real(dynmat(t,t)),0d0)
    endif
  do j=t+1,ndim-1
    if (abs(aimag(dynmat(t,j))+aimag(dynmat(j,t))) .gt. 9d-4*all ) then
       write(ulog,*)' dynmat is not hermitian in AIMAG of its off-diagonal elts'
       write(ulog,*)' off-diagonal element i,j=',t,j,dynmat(t,j),dynmat(j,t)
       write(ulog,*)' comparing to avg(abs(dynmat))=',all
!      stop
    elseif(abs(real(dynmat(t,j))-real(dynmat(j,t))) .gt. 9d-4*all ) then
       write(ulog,*)' dynmat is not hermitian in REAL of its off-diagonal elts'
       write(ulog,*)' off-diagonal element i,j=',t,j,dynmat(t,j),dynmat(j,t)
       write(ulog,*)' comparing to avg(abs(dynmat))=',all
!      stop
    else
    endif
! enforcing it to be hermitian!!
    mi=(aimag(dynmat(t,j))-aimag(dynmat(j,t)))/2
    mj=(real (dynmat(t,j))+real (dynmat(j,t)))/2
    dynmat(t,j) = cmplx(mj, mi)
    dynmat(j,t) = cmplx(mj,-mi)
  enddo
 enddo

! enforce it to be real if all imaginary componenents are very small
 all = sum(abs(aimag(dynmat(:,:))))/(ndim*ndim)
 if (all .lt. sum(abs(real(dynmat(:,:))))/(ndim*ndim)*1d-10) then
  do t=1,ndim
  do j=1,ndim
     dynmat(t,j)=(dynmat(t,j)+conjg(dynmat(t,j)))/2d0
  enddo
  enddo
 endif

! write(*,*)'Exiting set_dynamical_matrix'
 end subroutine set_dynamical_matrix
!==========================================================
 subroutine diagonalize(n,mat,eival,nv,eivec,ier)
! n=size of mat; nv is the number of needed eigenvectors
 implicit none
 integer n,ier,zero,nv  !,i,j
 complex(8) mat(n,n),eivec(n,n)
 real(8)  eival(n),tol
! This is used by eigch
! real(8), allocatable :: w(:,:)
! integer, allocatable :: lw(:)
! This is used by ZHEGV
  real(8), allocatable :: rwork(:)
  complex(8), allocatable :: work(:)
  integer lwork


! n = size(mat(:,1))
 if (n .ne. size(eival) ) then
    write(   *,*)' EIGCH, size inconsistency:mat,eival=',n, size(eival)
    stop
 endif

! tol = -1.d0
! zero= 0  ! for now do not compute any eigenvectors
! ier=0
! allocate(w(n,7),lw(n))
! call eigch(mat,n,n,n,nv,tol,w,lw,eival,eivec,ier)
! deallocate(w,lw)
! return

  lwork = max(1,2*n-1)
  allocate(work(lwork),rwork(max(1,3*n-2)))
  if(nv.eq.0) then
    call zheev('N','U',n,mat,n,eival,work,lwork,rwork,ier)
  else
    call zheev('V','U',n,mat,n,eival,work,lwork,rwork,ier)
    eivec = mat
  endif
  deallocate(work,rwork)
4 format(i5,1x,3(f6.3,1x,f6.3,4x))

 end subroutine diagonalize
!==========================================================
 subroutine write_eigenvalues(i,dk,kp,eival,n,uio)
 use constants
 implicit none
 integer i,j,n,uio
 integer sig(n)
 real(8), dimension(n) :: eival
 real(8) kp(3),dk

! n = size(eival)
! write(ueiv,3)i-1,kp,(eival(j),j=1,n)
 do j=1,n
    sig(j) = 1
    if (eival(j).lt.0) then
       eival(j) = -eival(j)
       sig(j) = -1d0
    endif
 enddo
 !write(ueiv,3)i-1,kp,(cnst*sig(j)*sqrt(eival(j)),j=1,n)
 write(uio,3)i,dk,kp,(cnst*sig(j)*sqrt(eival(j)),j=1,n)

3 format(i5,1x,4(1x,f8.3),999(1x,g11.5))
 end subroutine write_eigenvalues
!==========================================================
 subroutine store_2d_to_1d(eigenval,eival1d,n,m)
 use ios
 implicit none
 integer n1d,n,m,i,j,k
 real(8), dimension(n*m) :: eival1d
 real(8), dimension(n,m) :: eigenval
! n   = size(eigenval(:,1))
! m   = size(eigenval(1,:))
 n1d = n*m !size(eival1d)
 if (n1d .ne. n*m ) then
    write(   *,*)' STORE, size inconsistency:1d,2d=',n1d,n,m
    write(ulog,*)' STORE, size inconsistency:1d,2d=',n1d,n,m
    stop
 endif

 write(ulog,*) 'eival1d********************',n,m
 k = 0
 do j=1,m
 do i=1,n
    k = k+1
    eival1d(k) = eigenval(i,j)
    if (eival1d(k) .lt. 0) then
       write(ulog,*)'NEGATIVE eival ',k,eival1d(k)
    endif
 enddo
 enddo

 end subroutine store_2d_to_1d
!==========================================================
 subroutine calculate_dos_g(mx,eival,wkp,mesh,omega,ds)
! use gaussian broadening
 use ios
 use params
 use om_dos
 use constants
 implicit none
 integer i,j,mx,mesh
 real(8) x,wkp(mx),delta,ds(mesh),omega(mesh),eival(mx)  !,cnst

! wkp=1d0/(nkx*nky*nkz)
! write(udos,*)'# wkp,width =',wkp,width

    do i=1,mesh
       ds(i) = 0
       do j=1,mx
          x = (cnst*sqrt(abs(eival(j))) - omega(i))/width   ! these are all in cm^-1
!         x = (eival(j) - ene*ene)/width/width
          if ( abs(x) .gt. 5 ) cycle    !sy- if abs(x) > 5, don't execute the below line.
          ds(i) = ds(i) + delta(x)/width*wkp(j) !/mx  !sy - width from params.phon, wkp is weighting factor for IBZ
       enddo
    enddo

3 format(i5,9(3x,g11.5))

 end subroutine calculate_dos_g
!==========================================================
 subroutine calculate_jdos(mx,eival,wkp,mesh,omega,udosj)
! use gaussian broadening to calculate the joint dos
 use ios
 use params
 use om_dos
 use constants
 implicit none
 integer i,j,k,mx,mesh,udosj
 real(8) xp,xm,wkp(mx),delta_g,dsp,dsm,omega(mesh),eival(mx),iself,tmp,one
 complex(8) oc2,omz

! wkp=1d0/(nkx*nky*nkz)
 one = 1d0
 write(udosj,*)'# width =',width
 tmp =0.07      ! the temp width =0.07/cm is equiv to 0.1 K
! then (1+n1+n2)=1 should produce dsp results

! open(udosj+1,file='jdos.discrete',status='unknown')
! do j=1,mx
! do k=j,mx
!    xp= cnst*(sqrt(abs(eival(j)))+sqrt(abs(eival(k))))
!    xm= cnst*(sqrt(abs(eival(j)))-sqrt(abs(eival(k))))
!    write(udosj+1,2) xp,one
!!    if (xm.gt.0)
!    write(udosj+1,2) abs(xm),-one
! enddo
! enddo
! close(udosj+1)

    write(udosj,*)'# i,omega(i),jdsm(i),jdsp(i), occupation-weighted jdos'
    do i=1,mesh
       dsp=0 ; dsm=0
       omz=2*cmplx(omega(i),-etaz)
       do j=1,mx
       do k=1,mx
          xp= (cnst*(sqrt(abs(eival(j)))+sqrt(abs(eival(k)))) - 2*omega(i))   ! these are all in cm^-1
          xm= (cnst*(sqrt(abs(eival(j)))-sqrt(abs(eival(k)))) - 2*omega(i))   ! these are all in cm^-1
          if ( abs(xp) .lt. 6 ) then
             dsp = dsp + delta_g(xp,width)*wkp(j)*wkp(k) !* width*sqrt(2*pi)
          endif
          if ( abs(xm) .lt. 6 ) then
             dsm = dsm + delta_g(xm,width)*wkp(j)*wkp(k) !* width*sqrt(2*pi)
          endif

          iself = aimag(oc2(tmp,omz,j,k))/pi*wkp(j)*wkp(k)*tmp
       enddo
       enddo
       write(udosj,3) i,2*omega(i),dsm,dsp,iself
    enddo
    close(udosj)

2 format(9(3x,g11.5))
3 format(i5,9(3x,g11.5))

 end subroutine calculate_jdos
!==========================================================
 subroutine lifetime_dos_l(q,ndn,omq,dosp,dosm)
! use lorentzian broadening to calculate the joint dos
! for given q, calculates sum_1,2  delta(q-q1-/+q2)delta(w_q-w1-/+w2)
 use ios
 use params
 use om_dos
 use constants
 use eigen
 use kpoints
 implicit none
 integer, intent(in) :: ndn
 real(8), intent(in) ::  q(3)
 real(8), intent(out) ::  omq(ndn),dosp(ndn),dosm(ndn)
 integer l,la,la1,la2
 real(8) q1(3),qp(3),qm(3),delta_l,   &  
 &       om1(ndn),omp(ndn),omm(ndn),evl1(ndn),evlp(ndn),evlm(ndn)
 complex(8) evc1(ndn,ndn),evcp(ndn,ndn),evcm(ndn,ndn)

 call get_freq(q,ndn,evl1,evc1)
 omq = sqrt(abs(evl1)) * cnst   ! convert to 1/cm
 write(*,8) 'omq=',omq

 dosp=0; dosm=0
 do l=1,nkc
    q1=kpc(:,l)
    qp=-q-q1
    qm=-q+q1
    call get_freq(q1,ndn,evl1,evc1)
    call get_freq(qp,ndn,evlp,evcp)
    call get_freq(qm,ndn,evlm,evcm)
    om1 = sqrt(abs(evl1)) * cnst 
    omp = sqrt(abs(evlp)) * cnst 
    omm = sqrt(abs(evlm)) * cnst 

    do la1=1,ndn
    do la2=1,ndn
    do la =1,ndn

       dosp(la)=dosp(la)+delta_l(omq(la)-om1(la1)-omp(la2),etaz) 
       dosm(la)=dosm(la)+delta_l(omq(la)-om1(la1)+omm(la2),etaz)

    enddo
    enddo
    enddo
 enddo

! normalize and convert to ps to get in in per unit cell per THz (or ps)
 dosp=dosp/nkc*100*c_light*1e-12
 dosm=dosm/nkc*100*c_light*1e-12

! write(*,8)'LIFETIME_DOS: q,omqdosp,dosm=',q,omq,dosp,dosm
8 format(a,99(1x,g10.3))
 end subroutine lifetime_dos_l
!===========================================================
 subroutine phase_space_lifetime_bs
 use ios
 use om_dos
 use kpoints
 use eigen
 implicit none
 integer i,la
 real(8) q(3),omq(ndyn),dosp(ndyn),dosm(ndyn)

 write(*,*)'Enetering phase_space_bs'
 open(udos,file='bs_phase_space.dat ',status='unknown')
 write(udos,*)'# dk,q,(omega(q,la),jdsm(la)[ps],jdsp(la)[ps], la=1,ndyn'

 do i=1,nkp_bs
    q=kp_bs(:,i)
    call lifetime_dos_l(q,ndyn,omq,dosp,dosm)
    write(udos,3) dk_bs(i),q,(omq(la),dosm(la),dosp(la),la=1,ndyn)
 enddo
 close(udos)

3 format(99(1x,g11.4))

 end subroutine phase_space_lifetime_bs
!===========================================================
 subroutine phase_space_lifetime_FBZ
!! calculates for  given q the integral sum_k delta(w(q) pm w2(k) - w3(-k-q)) 
 use ios
 use om_dos
 use kpoints
 use eigen
 implicit none
 integer i,la
 real(8) q(3),omq(ndyn),dosp(ndyn),dosm(ndyn)

 open(udos,file='FBZ_phase_space.dat ',status='unknown')
 write(udos,*)'# q,(omega(q,la),jdsm(la),jdsp(la), la=1,ndyn'

 do i=1,nkc
    q=kpc(:,i)
    call lifetime_dos_l(q,ndyn,omq,dosp,dosm)
    write(udos,3) q,(omq(la),dosm(la),dosp(la),la=1,ndyn)
 enddo
 close(udos)

3 format(99(1x,g11.4))

 end subroutine phase_space_lifetime_FBZ
!===========================================================
 subroutine calculate_v3  (ndn,nk,kp,eival,eivec)
! check correctness of v3 by numerical differentiation of dynmat
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units are in cm^-1
! this works for two k's (k2,k3) in the prim cell on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! this subroutine calculates |V3(k1,la1,k2,la2,k3,la3)|^2  with k1=-k2-k3+G  and all ki in the primitive G-cell
! to avoid ambiguity in the case of degenerate eigenstates, and to satisfy symmetry exactly,
! not to mention saving time, we calculate only for 0<q_i<G/2 and complete by symmetry
!
 use kpoints
 use ios
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ndn,l1,l2,l3,n2,n3,al,be,ga,i0,j0,k0,j,k,t,ta1,ta2,ta3
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside  !,mp1,mp3
 integer nk1,nk2,nk3,np,nq,ns  !,i5,j5,k5,mk1,mk2,mk3
 real(8) q1(3),q2(3),q3(3),mi,mj,mk,const,rr2(3),rr3(3),den,eival(ndn,nk),kp(3,nk) !,eivl
! real(8), allocatable :: v33sq_5(:,:,:,:,:)
 complex(8) eivec(ndn,ndn,nk),xx,eiqr

 const = ee*1d30   ! convert ev/A^3 to SI units (kg s^-2 m^-1)
 const = const *(hbar/uma)**1.5d0  ! that takes care of hbar^1.5 and the 3 masses in the denominator
 const = const /(sqrt(ee*1d20/uma))**1.5d0 ! that takes care of the 3 frequencies in the denominator
 const = const/hbar ! to convert energy(J) to frequency
 const = const /200/pi/c_light  ! to convert frequency to cm^-1
 write(ulog,*)'CALCULATE_V3: to convert v3 to cm^-1, const=',const
! const = ee*1d30*(hbar/uma/(200*pi*c_light))**(1.5d0) ! cnvert v3 to SI if freqs in cm^-1
! const = const/(100*h_plank*c_light) ! to convert from SI to cm^-1
! write(ulog,*)'const2=',const

 write(ulog,*)'writing V3(k2,k3,l1,l2,l3) ##########################'

 allocate(v33sq_5(nk,nk,ndn,ndn,ndn))
 v33sq_5 = 0d0 !cmplx(0d0,0d0)
 write(ulog,*)'V3: nc(123)=',nc

 loop2: do n2=1,nk   ! second argument k2 needs to be only in the IFBZ
! loop2: do np=1,npos
!   n2=mappos(np)
    q2 = kp(:,n2)
    call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(n2.ne.nk2) then
      write(ulog,*)'n2,nk2,inside=',n2,nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif
!   call get_k_info(-q2,NC,mk2,i2,j2,k2,g1,g2,g3,inside)
!   write(*,*)' nk2,mk2=',nk2,mk2

 loop3: do n3=1 ,nk
! loop3: do nq=1,npos
!   n3=mappos(nq)
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif
!   call get_k_info(-q3,NC,mk3,i2,j2,k2,g1,g2,g3,inside)

    q1 = -q2-q3
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
!   aux1=conjg(eigenvec(:,:,nk1))
!   rloop:do nt=1,npos
!      if(nk1.eq.mappos(nt)) then ! nk1 corresponds to a pos-eig
!        aux1=eigenvec(:,:,nk1)
!        exit rloop
!      endif
!   enddo rloop
!   call get_k_info(-q1,NC,mk1,i1,j1,k1,g1,g2,g3,inside)

  write(*,2)'nk123=',nk1,nk2,nk3

 do l2=1 ,ndn
 do l3=1 ,ndn

    np = nk2+(l2-1)*nk
    nq = nk3+(l3-1)*nk
    if(np.gt.nq) cycle

 do l1=1 ,ndn
    ns = nk1+(l1-1)*nk
    if(ns.gt.np) cycle

    xx = cmplx(0d0,0d0)
    tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop

       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rr2(:) = atompos(:,j) - atompos(:,j0)
       rr3(:) = atompos(:,k) - atompos(:,k0)
!  be careful: it has to be the translations R not atompos!
!&      exp(ci*(kp(:,n2) .dot. atompos(:,j) + kp(:,n3) .dot. atompos(:,k))) /  &
       eiqr = exp( ci* ((q2 .dot. rr2) + (q3 .dot. rr3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
&      eivec(ta1,l1,nk1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
! make sure for kp at the FBZ boundary it gives (-1) or 1

! n2 is already in the IBZ, n3 is not but for now we are using the full mesh
! and nibz=nkc;  therefore n1 is also in the full zone and is OK, otherwise
! for n1 and n3 it is necessary to rotate the eigenvector from the IBZ to its
! true position in the recirpocal space. (the op_kmatrix will do the job)

    enddo tloop

    den = sqrt(sqrt(abs(eival(l1,nk1)*eival(l2,nk2)*eival(l3,nk3))))
    if (den.ne.0) then
       xx = xx / den/2d0/sqrt(2d0) * const
!      write(*,6)'nk_i,la_i,xx=',nk1,nk2,nk3,l1,l2,l3,xx
    else   ! usually acoustic eigenvalues at k=0 are not strictly zero, but a small number
       write(ulog,*)'calculate_v3: denominator=0 !! '
       write(ulog,*)'n1,n2,n3,l123=', nk1,n2,n3,l1,l2,l3
       write(ulog,*)'q1=', q1
       write(ulog,*)'q2=', q2
       write(ulog,*)'q3=', q3
       write(ulog,*)'eivs123=',eival(l2,n2),eival(l3,nk3),eival(l1,nk1)
       stop
    endif
  
    v33sq_5(nk2,nk3,l1,l2,l3)=xx  *conjg(xx)
    v33sq_5(nk3,nk2,l1,l3,l2)=xx  *conjg(xx)
    v33sq_5(nk1,nk3,l2,l1,l3)=xx  *conjg(xx)
    v33sq_5(nk3,nk1,l2,l3,l1)=xx  *conjg(xx)
    v33sq_5(nk2,nk1,l3,l2,l1)=xx  *conjg(xx)
    v33sq_5(nk1,nk2,l3,l1,l2)=xx  *conjg(xx)
!   xx = conjg(xx)
!   v3(mk2,mk3,l1,l2,l3)=xx
!   v3(mk3,mk2,l1,l3,l2)=xx
!   v3(mk1,mk3,l2,l1,l3)=xx
!   v3(mk3,mk1,l2,l3,l1)=xx
!   v3(mk2,mk1,l3,l2,l1)=xx
!   v3(mk1,mk2,l3,l1,l2)=xx

 enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop3
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop2

 if( writev3.eq.1) then
   open(uv3,file='v3.dat',status='unknown',FORM='UNFORMATTED')
!  open(uv3,file='v3.dat',status='unknown')
   write(uv3)nc,ndn
!  write(uv3,*)nc,ndn
   do n2=1 ,nk
   do n3=1 ,nk
   do l1=1 ,ndn
   do l2=1 ,ndn
   do l3=1 ,ndn
      if (v33sq_5(n2,n3,l1,l2,l3).gt.v3_threshold**2) then
         write(uv3  )n2,n3,l1,l2,l3,v33sq_5(n2,n3,l1,l2,l3)
      endif
   enddo
   enddo
   enddo
   enddo
   enddo
   close(uv3)
 endif
! v3 should be of the order of sub eV (v3=1 eV if phi_3=1000 eV/A^3)

2 format(a,2(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g14.8))
6 format(a,6(i5),9(2x,g14.8))

 deallocate(v33sq_5)

 end subroutine calculate_v3
!===========================================================
 subroutine calculate_v3_new(ndn,nk,kp,eival,eivec)
! check correctness of v3 by numerical differentiation of dynmat
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units are in cm^-1
! this works for two k's (k2,k3) in the prim cell on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! this subroutine calculates V3(k1,la1,k2,la2,k3,la3)  with k1=-k2-k3+G  and all ki in the primitive G-cell
! to avoid ambiguity in the case of degenerate eigenstates, and to satisfy symmetry exactly,
! not to mention saving time, we calculate only for 0<q_i<G/2 and complete by symmetry
!
 use kpoints
 use ios
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ndn,l1,l2,l3,n2,n3,al,be,ga,i0,j0,k0,j,k,t,ta1,ta2,ta3
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,np,nq,ns,indx
 real(8) q1(3),q2(3),q3(3),mi,mj,mk,const,rr2(3),rr3(3),den,eival(ndn,nk),kp(3,nk) !,eivl
 complex(8) eivec(ndn,ndn,nk),xx,eiqr

 const = ee*1d30   ! convert ev/A^3 to SI units (kg s^-2 m^-1)
 const = const *(hbar/uma)**1.5d0  ! that takes care of hbar^1.5 and the 3 masses in the denominator
 const = const /(sqrt(ee*1d20/uma))**1.5d0 ! that takes care of the 3 frequencies in the denominator
 const = const/hbar ! to convert energy(J) to frequency
 const = const /200/pi/c_light  ! to convert frequency to cm^-1
 write(ulog,*)'CALCULATE_V3: to convert v3 to cm^-1, const=',const
! const = ee*1d30*(hbar/uma/(200*pi*c_light))**(1.5d0) ! cnvert v3 to SI if freqs in cm^-1
! const = const/(100*h_plank*c_light) ! to convert from SI to cm^-1
! write(ulog,*)'const2=',const

 write(ulog,*)'writing V3sq(k2,k3,l1,l2,l3) ##########################'

 v33sq = 0d0 ! cmplx(0d0,0d0)
 indx=0
 write(ulog,*)'V3: nc(123)=',nc

 loop2: do n2=1,nk   ! second argument k2 needs to be only in the IFBZ
! loop2: do np=1,npos
!   n2=mappos(np)
    q2 = kp(:,n2)
    call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(n2.ne.nk2) then
      write(ulog,*)'n2,nk2,inside=',n2,nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif
!   call get_k_info(-q2,NC,mk2,i2,j2,k2,g1,g2,g3,inside)
!   write(*,*)' nk2,mk2=',nk2,mk2

 loop3: do n3=1 ,nk
! loop3: do nq=1,npos
!   n3=mappos(nq)
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif
!   call get_k_info(-q3,NC,mk3,i2,j2,k2,g1,g2,g3,inside)

    q1 = -q2-q3
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
!   aux1=conjg(eigenvec(:,:,nk1))
!   rloop:do nt=1,npos
!      if(nk1.eq.mappos(nt)) then ! nk1 corresponds to a pos-eig
!        aux1=eigenvec(:,:,nk1)
!        exit rloop
!      endif
!   enddo rloop
!   call get_k_info(-q1,NC,mk1,i1,j1,k1,g1,g2,g3,inside)

  write(*,2)'nk123=',nk1,nk2,nk3

 do l2=1 ,ndn
 do l3=1 ,ndn

    np = nk2+(l2-1)*nk
    nq = nk3+(l3-1)*nk
    if(np.gt.nq) cycle

 do l1=1 ,ndn
    ns = nk1+(l1-1)*nk
    if(ns.gt.np) cycle

    xx = cmplx(0d0,0d0)
    tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop

       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rr2(:) = atompos(:,j) - atompos(:,j0)
       rr3(:) = atompos(:,k) - atompos(:,k0)
!  be careful: it has to be the translations R not atompos!
!&      exp(ci*(kp(:,n2) .dot. atompos(:,j) + kp(:,n3) .dot. atompos(:,k))) /  &
       eiqr = exp( ci* ((q2 .dot. rr2) + (q3 .dot. rr3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
&      eivec(ta1,l1,nk1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
! make sure for kp at the FBZ boundary it gives (-1) or 1

! n2 is already in the IBZ, n3 is not but for now we are using the full mesh
! and nibz=nkc;  therefore n1 is also in the full zone and is OK, otherwise
! for n1 and n3 it is necessary to rotate the eigenvector from the IBZ to its
! true position in the recirpocal space. (the op_kmatrix will do the job)

    enddo tloop

    den = sqrt(sqrt(abs(eival(l1,nk1)*eival(l2,nk2)*eival(l3,nk3))))
    if (den.ne.0) then
       xx = xx / den/2d0/sqrt(2d0) * const
!      write(*,6)'nk_i,la_i,xx=',nk1,nk2,nk3,l1,l2,l3,xx
    else   ! usually acoustic eigenvalues at k=0 are not strictly zero, but a small number
       write(ulog,*)'calculate_v3: denominator=0 !! '
       write(ulog,*)'n1,n2,n3,l123=', nk1,n2,n3,l1,l2,l3
       write(ulog,*)'q1=', q1
       write(ulog,*)'q2=', q2
       write(ulog,*)'q3=', q3
       write(ulog,*)'eivs123=',eival(l2,n2),eival(l3,nk3),eival(l1,nk1)
       stop
    endif

    if (abs(xx).gt.v3_threshold) then
       indx=indx+1
       v33sq(indx)=xx * conjg(xx)
       nq1(indx)= nk1
       nq2(indx)= nk2
       nq3(indx)= nk3
       la1(indx)= l1
       la2(indx)= l2
       la3(indx)= l3
    endif
 enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop3
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop2

 nv3 = indx
 write(ulog,*)' V33: total size of this array is =',nv3

 if( writev3.eq.1) then
   open(uv3,file='v33.dat',status='unknown',FORM='UNFORMATTED')
!  open(uv3,file='v33.dat',status='unknown')
   write(uv3)nv3
!  write(uv3,*)nv3
   do j=1 ,nv3
      write(uv3) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33sq(j)
!     write(uv3,6)' ',nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
   enddo
   close(uv3)
 endif
! v3 should be of the order of sub eV (v3=1 eV if phi_3=1000 eV/A^3)

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))

 end subroutine calculate_v3_new
!===========================================================
 subroutine calculate_v3_ibz(ndn,ni,ki,nk,kp,eival,eivec)
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units for V3 are in cm^-1
! this works for two k's (k2,k3) in the prim cell on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! this subroutine calculates V3(k1,la1,k2,la2,k3,la3)  with k1=-k2-k3+G  and all ki in the primitive G-cell
! with the restriction of k2 in the IRREDUCIBLE BZ defined by ki(3,ni) included in kp(3,nk)
!
 use kpoints
 use ios
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ni,ndn,l1,l2,l3,n2,n3,al,be,ga,i0,j0,k0,j,k,t,ta1,ta2,ta3
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,np,nq,ns,indx
 real(8) q1(3),q2(3),q3(3),mi,mj,mk,const,rr2(3),rr3(3),den,eival(ndn,nk)
 real(8) ki(3,ni),kp(3,nk)
 complex(8) eivec(ndn,ndn,nk),xx,eiqr

!*********************************************************************
!  NEED TO DEFINE EIVAL AND EIVEC OVER THE WHOLE ZONE
!*********************************************************************

 const = 1d30   ! convert ev/A^3 to eV ######## SI units (kg s^-2 m^-1)
 const = const *(hbar/uma)**1.5d0  ! that takes care of hbar^1.5 and the 3 masses in the denominator
 const = const /(sqrt(ee*1d20/uma))**1.5d0 ! that takes care of the 3 frequencies in the denominator
 const = const * ee  ! to convert energy from eV to energy Joules
 const = const/h_plank ! to convert energy(J) to frequency(Hz)
 const = const /100/c_light  ! to convert frequency(Hz) to cm^-1
 write(ulog,*)'CALCULATE_V3: to convert v3 to cm^-1, const=',const
 write(ulog,*)'writing V3(k2,k3,l1,l2,l3) ##########################'

 v33sq = 0d0 !cmplx(0d0,0d0)
 indx=0
 write(ulog,*)'V3: nc(123)=',nc

!************************************************
! second argument q2 needs to be only in the IFBZ
!************************************************
 loop2: do n2=1,ni
    q2 = ki(:,n2)  ! should be = kp(:,mapinv(n2))
    call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(mapinv(n2).ne.nk2) then
      write(ulog,*)'n2,mapinv(n2),nk2,inside=',n2,mapinv(n2),nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif
!   call get_k_info(-q2,NC,mk2,i2,j2,k2,g1,g2,g3,inside)
!   write(*,*)' nk2,mk2=',nk2,mk2

 loop3: do n3=1 ,nk
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif
!   call get_k_info(-q3,NC,mk3,i2,j2,k2,g1,g2,g3,inside)

    q1 = -q2-q3

    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)

  write(*,2)'nk123=',nk1,nk2,nk3

 do l2=1 ,ndn
 do l3=1 ,ndn

    np = nk2+(l2-1)*nk ! ni
    nq = nk3+(l3-1)*nk
!   if(np.gt.nq) cycle

 do l1=1 ,ndn
    ns = nk1+(l1-1)*nk
!   if(ns.gt.np) cycle

    xx = cmplx(0d0,0d0)
    tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop

       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rr2(:) = atompos(:,j) - atompos(:,j0)
       rr3(:) = atompos(:,k) - atompos(:,k0)
!  be careful: it has to be the translations R not atompos!
!&      exp(ci*(kp(:,n2) .dot. atompos(:,j) + kp(:,n3) .dot. atompos(:,k))) /  &
       eiqr = exp( ci* ((q2 .dot. rr2) + (q3 .dot. rr3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
&      eivec(ta1,l1,nk1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
! make sure for kp at the FBZ boundary it gives (-1) or 1

! n2 is already in the IBZ, n3 is not
! therefore n1 is also in the full zone and is OK, otherwise
! for n1 and n3 it is necessary to rotate the eigenvector from the IBZ to its
! true position in the recirpocal space. (the op_kmatrix will do the job)

    enddo tloop

    den = sqrt(sqrt(abs(eival(l1,nk1)*eival(l2,nk2)*eival(l3,nk3))))
    if (den.ne.0) then
       xx = xx / den/2d0/sqrt(2d0) * const
!      write(*,6)'nk_i,la_i,xx=',nk1,nk2,nk3,l1,l2,l3,xx
    else   ! usually acoustic eigenvalues at k=0 are not strictly zero, but a small number
       write(ulog,*)'calculate_v3: denominator=0 !! '
       write(ulog,*)'n1,n2,n3,l123=', nk1,nk2,nk3,l1,l2,l3
       write(ulog,*)'q1=', q1
       write(ulog,*)'q2=', q2
       write(ulog,*)'q3=', q3
       write(ulog,*)'eivs123=',eival(l1,nk1),eival(l2,nk2),eival(l3,nk3)
       stop
    endif

!    if (abs(xx).gt.v3_threshold) then
       indx=indx+1
       v33sq(indx)=xx*conjg(xx)
       nq1(indx)= nk1
       nq2(indx)= nk2
       nq3(indx)= nk3
       la1(indx)= l1
       la2(indx)= l2
       la3(indx)= l3
!      if (mod(indx,1000).eq.1) write(*,*) 'indx,v33=',indx,xx
!    endif
 enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop3
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop2

 nv3 = indx
 write(ulog,*)' V33: total size of this array is =',nv3

 if( writev3.eq.1) then
   open(uv3,file='v33.dat',status='unknown',FORM='UNFORMATTED')
   write(uv3)nv3
   do j=1 ,nv3
      write(uv3) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33sq(j)
   enddo
!  open(uv3,file='v33.dat',status='unknown')
!  write(uv3,*)nv3
!  do j=1 ,nv3
!     write(uv3,7) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
!  enddo
   close(uv3)
 endif
! v3 should be of the order of sub eV (v3=1 eV if phi_3=1000 eV/A^3)

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))

 end subroutine calculate_v3_ibz
!============================================================
 subroutine read_v3_all_formatted(nk2,nk3,ndn)
 use phi3
 use ios
 use lattice
 implicit none
 integer i,k2,k3,l1,l2,l3,nk2,nk3,ndn,unt ,ncc,n1,n2,n3,nd
 real(8) xx,yy

 unt = 111
 open(unt,file='v3.dat',status='old')

 read(unt,*,end=99) n1,n2,n3,nd
 write(ulog,*)' OPENING v3.dat, reading ',n1,n2,n3,nd
 ncc =  n1*n2*n3
 if (ncc .ne. nk2 .or. nd .ne. ndn ) then
    write(*,*) 'nk in file v3 is different from nc1*nc2*nc3',ncc,nk2
    write(ulog,*) 'nk in file v3 is different from nc1*nc2*nc3',ncc,nk2
    write(*,*) 'or ndyn in the file is different from that of the code ',nd,ndn
    write(ulog,*) 'or ndyn in the file is different from that of the code ',nd,ndn
 endif
 v33sq_5 = 0d0 !cmplx(0d0,0d0)
 do i=1,nk2*nk3*ndn*ndn*ndn
    read(unt,*,end=99) k2,k3,l1,l2,l3,xx,yy
    v33sq_5(k2,k3,l1,l2,l3)=xx*xx+yy*yy !cmplx(xx,yy)
 enddo
 close(unt)
 return

99 write(*,*)'v3 file end reached at line i=',i
   write(*,*)'k2,k3,l1,l2,l3=',k2,k3,l1,l2,l3
   stop

 end subroutine read_v3_all_formatted

!============================================================
 subroutine read_v3_all_unformatted(nk2,nk3,ndn)
 use phi3
 use ios
 use lattice
 implicit none
 integer i,k2,k3,l1,l2,l3,nk2,nk3,ndn,unt ,ncc,n1,n2,n3,nd

 unt = 111
!open(unt,file='v3.dat',status='old')
 open(unt,file='v3.dat',status='old',form='UNFORMATTED')

 read(unt,end=99) n1,n2,n3,nd
 write(ulog,*)' OPENING v3.dat, reading ',n1,n2,n3,nd
 ncc =  n1*n2*n3
 allocate(v33sq_5(nk2,nk3,ndn,ndn,ndn))
 if (ncc .ne. nk2 .or. nd .ne. ndn ) then
    write(*,*) 'nk in file v3 is different from nc1*nc2*nc3',ncc,nk2
    write(ulog,*) 'nk in file v3 is different from nc1*nc2*nc3',ncc,nk2
    write(*,*) 'or ndyn in the file is different from that of the code ',nd,ndn
    write(ulog,*) 'or ndyn in the file is different from that of the code ',nd,ndn
 endif
 v33sq_5 = 0d0 !cmplx(0d0,0d0)
 do i=1,nk2*nk3*ndn*ndn*ndn
    read(unt,end=99) k2,k3,l1,l2,l3,v33sq_5(k2,k3,l1,l2,l3)  !,yy 
 enddo
 close(unt)
 return

99 write(*,*)'v3 file end reached at line i=',i
   write(*,*)'k2,k3,l1,l2,l3=',k2,k3,l1,l2,l3
   stop

 end subroutine read_v3_all_unformatted
!============================================================
 subroutine read_v3sq_new_unformatted
 use phi3
 use ios
 use lattice
 implicit none
 integer j,unt
 real(8) xx,yy

 unt = 111
 open(unt,file='v33.dat',status='old',form='UNFORMATTED')
!open(unt,file='v33.dat',status='old')

 read(unt,end=99) nv3
! read(unt,*,end=99) nv3
 call allocate_v33sq(nv3)
 write(ulog,*)' OPENING v33.dat, reading ',nv3
 v33sq = 0d0 
 do j=1,nv3
    read(unt,end=99) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33sq(j)
 enddo
 close(unt)
 return

99 write(*,*)'v33 file end reached at line j=',j
   stop

 end subroutine read_v3sq_new_unformatted
!============================================================
 subroutine read_v3_new_formatted
 use phi3
 use ios
 use lattice
 implicit none
 integer j,unt

 unt = 111
 open(unt,file='v33.dat',status='old')

 read(unt,*,end=99) nv3
 call allocate_v33sq(nv3)
 write(ulog,*)' OPENING v33.dat, reading ',nv3
 v33sq = 0d0 ! cmplx(0d0,0d0)
 do j=1,nv3
    read(unt,*,end=99) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33sq(j) 
 enddo
 close(unt)
 return

99 write(*,*)'READ_V3_NEW: v33 file end reached at line j=',j
   stop

 end subroutine read_v3_new_formatted
!============================================================
 subroutine check_sym_of_v3(ndyn)
! check its invariance under permutations of (k_i,la_i)
! v3(k2,k3,l1,l2,l3)=v3(k3,k2,l1,l3,l2) and  v3(k2,k3,l1,l2,l3)=v3(-k2-k3,k3,l2,l1,l3)
! and  v3(k2,k3,l1,l2,l3)=v3(k2,-k2-k3,l3,l2,l1) ; and v3(-k_i,la_i)=conjg(v3(k_i,la_i))
! but the first k index only runs within the IFBZ, and the second one through the FBZ
 use kpoints
 use lattice
 use phi3      ! v3 is declared as : v3(nibz,nkc,ndyn,ndyn,ndyn)
 implicit none
 integer ndyn,l1,l2,l3,n1,n2,n3,mk2,mk3,i1,j1,k1,i2,j2,k2,i3,j3,k3,inside
 real(8) q1(3),tol,dif,xx,yy
! complex(8) 

 tol = 0.0001

 do l1=1,ndyn
 do l2=1,ndyn
 do l3=1,ndyn
 do n2=1,nkc
 do n3=1,nibz   !,nkc

  call get_k_info(-kpc(:,n2),NC,mk2,i2,j2,k2,g1,g2,g3,inside)
  call get_k_info(-kpc(:,n3),NC,mk3,i3,j3,k3,g1,g2,g3,inside)

  if (mk3.le.nibz) then
    xx = v33sq_5(mk3,mk2,l1,l3,l2)
    yy = v33sq_5(n3,n2,l1,l3,l2)
    dif = xx-yy !xx-conjg(yy)
    if ( abs( dif ) .gt. tol ) call err4(n3,n2,mk3,mk2,l1,l3,l2,xx,yy)
  endif

  if (n2.le.nibz) then
    dif= v33sq_5(n2,n3,l1,l2,l3)-v33sq_5(n3,n2,l1,l3,l2)
    if ( abs( dif ) .gt. tol ) call err3(1,n3,n2,l1,l3,l2,dif)
  endif

  q1(:) = -kpc(:,n2)-kpc(:,n3)   ! there is inversion symmetry in k-space
  call get_k_info(q1,NC,n1,i1,j1,k1,g1,g2,g3,inside)
  if (n1.le.nibz) then
    dif = v33sq_5(n3,n2,l1,l3,l2)-v33sq_5(n1,n2,l3,l1,l2)
    if ( abs( dif ) .gt. tol ) call err3(2,n3,n2,l1,l3,l2,dif)
    dif = v33sq_5(n3,n2,l1,l3,l2)-v33sq_5(n3,n1,l2,l3,l1)
    if ( abs( dif ) .gt. tol ) call err3(3,n3,n2,l1,l3,l2,dif)
  endif

 enddo
 enddo
 enddo
 enddo
 enddo

 end subroutine check_sym_of_v3
!===========================================================
 subroutine err3(i,n3,n2,l1,l3,l2,dif)
 integer i,n3,n2,l1,l3,l2
 real(8) dif
 write(100,3)' CHECK_SUM_OF_V3: ERROR index=',i,n3,n2,l1,l3,l2,dif
3 format(a,i2,2(1x,i6),2x,3(1x,i3),2x,9(1x,g11.5))
 end subroutine err3
!===========================================================
 subroutine err4(n3,n2,m3,m2,l1,l3,l2,x,y)
 integer n3,n2,l1,l3,l2,m3,m2
 real(8) x,y
 write(100,3)' CHECK_SUM_OF_V3 ERR4:n3,-n3,n2,-n2,l_i =',n3,m3,n2,m2,l1,l3,l2,x,y
3 format(a,2(1x,i5),2x,2(1x,i5),3x,3(1x,i3),2x,9(1x,f10.4))
 end subroutine err4
!===========================================================
 subroutine calculate_kappa(temp,kappau,kappat) !,kp,nkc)
 use kpoints
 use ios
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 implicit none
 integer w,ik,la !,nbmx,nbq(nbkmax),npar  ! ,nkp
! integer, allocatable:: nb(:,:)
 real(8) temp,kappau,kappat,omg,nbw,ckl,nbe,v23,gammau,gammat  !,vk(3)grad_intp(3),
 real(8) kcut,q(3)
 complex(8) uself,nself

! here one can afford to take a much finer mesh and interpolate or use the MP method

! need to define kcut for the fine mesh
 kcut = 2.5*length(kibz(:,1)-kibz(:,2))
 write(ulog,*)' KAPPA: kcut for fine mesh=',kcut

 kappat = 0; kappau = 0
 kloop: do ik=1,nbz_fine
    write(*,*)' kloop in kappa =',ik,nbz_fine

    q = kibz(:,ik)

!********  do something about this

!    call find_nearest_kpoints3(nibz_coarse,kibz,q,nbmx,nbq,kcut,npar)

!********  do something about this

    lloop: do la=1,ndyn
       omg = sqrt(abs(eigenval_f(la,ik))) * cnst     ! in same units as wmax and temp
! use the on-shell frequencies : find the corresponding w
       w = 1+int((omg/wmax)*wmesh)
       nbw = nbe(omg,temp,classical)
       ckl = (omg/temp)**2 * nbw*(nbw+1) * k_b      ! heat capacity per mode
       v23 = (veloc_f(:,la,ik) .dot. veloc_f(:,la,ik))/3
!      call function_self(q,la,omg,temp,nbmx,nbq,nself,uself,npar)
       call function_self_new(q,la,omg,temp,nself,uself)
!      uself = 0 ; nself = 0
       gammau = aimag(uself)
       gammat = aimag(uself+nself)
!      kappa = kappa + v23 /(2*gammak * c_light*100 + 1/tau0) *  &
!&                   ckl * (2*pi*c_light)**2
       kappat = kappat + v23 * ckl /(2*gammat*100 + 1/tau0/c_light) *  c_light
       kappau = kappau + v23 * ckl /(2*gammau*100 + 1/tau0/c_light) *  c_light
    enddo lloop
 enddo kloop
 kappau = kappau/nbz_fine*volume_g*1d30
 kappat = kappat/nbz_fine*volume_g*1d30

! deallocate(nb)
 end subroutine calculate_kappa
!===========================================================
 subroutine find_nearest_kpoints(nx,ny,nz,q,nbq,nbs,nbmx)
! this is very fast but works only for a cubic mesh
 use constants
! use kpoints
 use lattice
 use params
 implicit none
 real(8) q(3),a1,a2,a3,ep
 integer n1,n2,n3,nbmx,na1,na2,na3,e1,e2,e3,indexn,nx,ny,nz
 integer nbq(nbmx),nbs(nbmx)

 ep = 5*tolerance
 q = q - ep
! write q= ai*gi
 a1 = (q .dot. r1)/(2*pi)
 a2 = (q .dot. r2)/(2*pi)
 a3 = (q .dot. r3)/(2*pi)
 q = q + ep
! now find the 8 nearest kpoints to it
 n1 = floor(a1*nx)
 n2 = floor(a2*ny)
 n3 = floor(a3*nz)

! these are the indices of the 8 neighbors : kp(:,nbq(j)), j=1,8
 nbq(1) = indexn(n1  ,n2  ,n3  ,nx,ny,nz)
 nbq(2) = indexn(n1  ,n2  ,n3+1,nx,ny,nz)
 nbq(3) = indexn(n1  ,n2+1,n3  ,nx,ny,nz)
 nbq(4) = indexn(n1  ,n2+1,n3+1,nx,ny,nz)
 nbq(5) = indexn(n1+1,n2  ,n3  ,nx,ny,nz)
 nbq(6) = indexn(n1+1,n2  ,n3+1,nx,ny,nz)
 nbq(7) = indexn(n1+1,n2+1,n3  ,nx,ny,nz)
 nbq(8) = indexn(n1+1,n2+1,n3+1,nx,ny,nz)

! along each of the 3 directions find whether ai*Ni is closer to ni or ni+1
! (ei=+1;nai=ni+1) if ni+1 is chosen and (-1,ni) if ni is chosen
 if(a1*nx-n1 .gt. n1+1-a1*nx) then  ! include ni+2
   e1 = 1 ; na1 = n1+1
 else
   e1 =-1 ; na1 = n1
 endif
 if(a2*ny-n2 .gt. n2+1-a2*ny) then  ! include ni+2
   e2 = 1 ; na2 = n2+1
 else
   e2 =-1 ; na2 = n2
 endif
 if(a3*nz-n3 .gt. n3+1-a3*nz) then  ! include ni+2
   e3 = 1 ; na3 = n3+1
 else
   e3 =-1 ; na3 = n3
 endif

 nbq(9) = indexn(na1+e1,na2   ,na3   ,nx,ny,nz)
 nbq(10)= indexn(na1   ,na2+e2,na3   ,nx,ny,nz)
 nbq(11)= indexn(na1   ,na2   ,na3+e3,nx,ny,nz)
 nbq(12)= indexn(na1   ,na2+e2,na3+e3,nx,ny,nz)
 nbq(13)= indexn(na1+e1,na2   ,na3+e3,nx,ny,nz)
 nbq(14)= indexn(na1+e1,na2+e2,na3   ,nx,ny,nz)
 nbq(15)= indexn(na1+e1,na2+e2,na3+e3,nx,ny,nz)

 nbs(:)=1

 end subroutine find_nearest_kpoints
!===========================================================
 subroutine cubic_rates(q,la,temp,normal,umklapp)
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k )
! p(q,lambda)=sum_kla1la2 |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh
 use kpoints
 use ios
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp
 real(8), intent(out) :: normal,umklapp
 integer nq,nk1,l1,l2,iq,ik,jq,jk,kq,kk,i4,j4,k4,nk2,nk,inside  !,index_reg
 real(8) om1,om2,omq,k1(3),k2(3),nb1,nb2,nbe,nbq,delta_g,rate,w3  !,kin(3),absorption,emission

 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
 omq = sqrt(abs(eigenval(la,nq))) * cnst  ! should be the same as without mapibz
 nbq = nbe(omq,temp,classical)

 normal=0; umklapp=0
 k1loop: do nk1=1,nkc
    k1(:)= kpc(:,nk1)
    call get_k_info(k1,nc,nk,ik,jk,kk,g1,g2,g3,inside)
    if(nk.ne.nk1) then
       write(ulog,*)'GET_RATES: nk.ne.nk1 ',nk,nk1
       write(ulog,*)'GET_RATES: i,j,k,inside=',ik,jk,kk,inside
       stop
    endif

    bnd1loop: do l1=1,ndyn
     om1 = sqrt(abs(eigenval(l1,nk1))) * cnst ; nb1 = nbe(om1,temp,classical)

    bnd2loop: do l2=1,ndyn

! below are the q=>k1+k2 terms of creation or annihilation of q
      k2(:)= q(:)-k1(:)
      call get_k_info(k2,nc,nk2,i4,j4,k4,g1,g2,g3,inside)
      om2 = sqrt(abs(eigenval(l2,nk2))) * cnst ; nb2=nbe(om2,temp,classical)

      w3 =v33sq_5(nk1,nk2,la,l1,l2)
      rate = 0.5*w3 * ((nbq+1)*nb1*nb2-nbq*(nb1+1)*(nb2+1))  &
&            *delta_g(omq-om1-om2,etaz)
      if(inside.eq.1) then                   ! normal process
          normal = normal + rate
      else                                ! umklapp process
          umklapp = umklapp + rate
      endif

! below are the q+k1=>k2 terms of absorption or emission
      k2(:)= q(:)+k1(:)
      call get_k_info(k2,nc,nk2,i4,j4,k4,g1,g2,g3,inside) ! use this inside
      om2 = sqrt(abs(eigenval(l2,mapibz(nk2)))) * cnst ; nb2=nbe(om2,temp,classical)
      w3 =v33sq_5(nk1,nk2,la,l1,l2)
      rate = w3 * ((nbq+1)*(nb1+1)*nb2-nbq*nb1*(nb2+1))  &
&            *delta_g(omq+om1-om2,etaz)
      if(inside.eq.1) then                   ! normal process
          normal = normal + rate
      else                                ! umklappa process
          umklapp = umklapp + rate
      endif

    enddo bnd2loop
    enddo bnd1loop

 enddo k1loop

  normal = normal   /nkc*pi/4 ! (4d4*pi*pi*c_light*c_light * nkc)*pi/4
  umklapp = umklapp /nkc*pi/4 ! (4d4*pi*pi*c_light*c_light * nkc)*pi/4
  write(urate,3)nq,la,omq,normal,umklapp

3 format(i6,i6,3x,g11.5,1x,4(2x,2(1x,g11.5)))

 end subroutine cubic_rates
!============================================================
 subroutine function_self_old(q,la,omega,temp,nself,uself)
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh
 use kpoints
 use ios
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer nq,n1,l1,l2,iq,i3,jq,j3,kq,m3,inside,nk1,ik,jk,kk,nk2,mk2
 real(8) om1,om2,k1(3),k2(3),nb1,nb2,nbe,eta,tk,etas,tr,ti,w12,sr,si,vr,vi,w33,term
 complex(8) omz,sumz,self

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
 etas = eta(omega,tk)
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! om1 = sqrt(abs(eigenval(la,nq))) * cnst
! etaz= (deltak*vgmax)**2/omega !min(etaz,100/omega)
  write(ulog,4)'SELF: omega,etaz,(etas)=',omega,etaz,etas
  omz = cmplx(omega, etaz)   ! this is in cm^-1
! write(ulog,*)'FUNCTION_SELF: OMZ=',omz
 nself=0 ; uself=0

 ti=2*etaz*omega
!     write(117,4)' ======================== NEW KPOINT ', q

 k1loop: do n1=1,nkc
    k1(:)= kpc(:,n1)
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,inside)
    if (nk1.ne.n1) then
       write(ulog,*) 'Function self: nk1 .ne. n1 ',nk1,n1
       stop
    endif
!   mp1 = mapibz(n1)

    k2(:)=-kpc(:,nq)-k1(:)   !these are the normal processes if k2 is within the FBZ
    call get_k_info(k2,nc,nk2,i3,j3,m3,g1,g2,g3,inside)
!   mp2 = mapibz(nk2)
    call get_k_info(-k2,nc,mk2,i3,j3,m3,g1,g2,g3,inside) ! take inside from this line

    sumz = 0
    bnd1loop: do l1=1,ndyn
    bnd2loop: do l2=1,ndyn

!     om1 = sqrt(abs(eigenval(l1,mp1))) * cnst ; nb1=nbe(om1,temp)
!     om2 = sqrt(abs(eigenval(l2,mp2))) * cnst ; nb2=nbe(om2,temp)
      om1 = sqrt(abs(eigenval(l1,nk1))) * cnst ; nb1=nbe(om1,temp,classical)
      om2 = sqrt(abs(eigenval(l2,nk2))) * cnst ; nb2=nbe(om2,temp,classical)

      w33 = v33sq_5(nk1,nk2,la,l1,l2)
      if (abs(w33) .lt. v3_threshold**2)  cycle bnd2loop
!     omz = cmplx(omega, eta(om1,tk)+eta(om2,tk)+eta(omega,tk))   ! this is in cm^-1

      w12=(om1+om2)
      tr=(w12*w12-omega*omega+etaz*etaz)   ! also try it without etaz**2
      term= 2*w12 /(tr*tr+ti*ti)
      sr=term*tr
      si=term*ti

      w12=(om1-om2)
      tr=(w12*w12-omega*omega+etaz*etaz)
      term= 2*w12 /(tr*tr+ti*ti)
      vr=term*tr
      vi=term*ti
      sumz = sumz - w33 * cmplx((nb2+nb1+1)*sr+(nb2-nb1)*vr,(nb2+nb1+1)*si+(nb2-nb1)*vi)
! &              ((nb2+nb1+1)*(1d0/(om1+om2-omz)+ 1d0/(om1+om2+omz)) +  &
! &               (nb2-nb1)  *(1d0/(om1-om2-omz)+ 1d0/(om1-om2+omz)))
!&             ((nb2+nb1+1)*cmplx(sr,si) +  (nb2-nb1)  *cmplx(vr,vi))

!     write(117,4)' ',term,(nb2+nb1+1)*sr+(nb2-nb1)*vr,(nb2+nb1+1)*si+(nb2-nb1)*vi,sumz
    enddo bnd2loop
    enddo bnd1loop

    if(inside.eq.1) then                   ! normal process
        nself = nself + sumz
    else                                ! umklappa process
        uself = uself + sumz
    endif

!    call check_inside_fbz(k2,g1,g2,g3,inside)
!    call send_to_primcell(-q-k1,kk)
!    if (length(kk+q+k1).lt.1d-4) inside=.True.

 enddo k1loop

  nself = nself /nkc /2 *2*pi
  uself = uself /nkc /2 *2*pi
  self = nself + uself
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))

 end subroutine function_self_old
!============================================================
 subroutine function_self_new(q,la,omega,temp,nself,uself)
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use ios
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1
 integer, save :: cntr
 real(8) om1,eta,tk,v32,term,k1(3),k2(3),q1(3)
 complex(8) omz,self,ocfunc

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! om1 = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 write(ulog,5)'SELF_NEW: last cntr, omega,etaz=',cntr,omega,etaz
! omz = cmplx(omega,-etaz)   ! this is in cm^-1
 v32 = v3_threshold*v3_threshold
 nself=0 ; uself=0;
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 jq = nq+(la-1)*nkc
 cntr = 0
 nv3loop: do j=1,nv3
!  q must be the 2nd argument in v33
    if (nq.ne.nq2(j)) cycle nv3loop
    if (la.ne.la2(j)) cycle nv3loop
    term = v33sq(j) !v33(j)*conjg(v33(j))
    if (term.lt.v32) cycle nv3loop

    jk1 = nq1(j)+(la1(j)-1)*nkc
    jk2 = nq2(j)+(la2(j)-1)*nkc
    jk3 = nq3(j)+(la3(j)-1)*nkc
    if(jk2 .ne. jq) print*,'ERRRRRRRRRRORRRRRRRRRRE : jk2 ne jq ',jk2,jq

    k2=kpc(:,nq3(j))
!   q1=kpc(:,nq1(j))
    k1 = -q-k2
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,inside)

!   if(mod(j,36).eq.12) write(ulog,7)'k1,k2,k3=',nk1,nq,nq3(j),inside,k1,q,k2

    if(nk1.ne.nq1(j)) then
      write(ulog,3)'SELF_NEW: nk1.ne.nq1 ',nk1,nq1(j)
      write(ulog,3)'Check consistency between nki and those of v33.dat'
      write(ulog,3)'SELF_NEW: j,la,q,om=',j,la,q,omega
      write(ulog,4)'SELF_NEW: k1,q,k2=',k1,q,k2
      write(ulog,*)'SELF_NEW: nqi    =',nq1(j),nq2(j),nq3(j)
      write(ulog,4)'SELF_NEW:kpc(nqi)=',kpc(:,nq1(j)),kpc(:,nq2(j)),kpc(:,nq3(j))
      write(ulog,4)'SELF_NEW:la,la2(j)=',la,la2(j)
      write(*,3)'SELF_NEW: nk1.ne.nq1 ',nk1,nq1(j)
      write(*,3)'SELF_NEW: j,la,q,om=',j,la,q,omega
      stop
    endif

!    if (jq .eq. jk1) then   ! q is the 1st argument in v33
!       self = ocfunc(temp,omz,jk2,jk3)
!       sumz = sumz + term * self
!!      call get_k_info(kpc(:,nq1(j)),nc,nq,ik,jk,kk,g1,g2,g3,inside)
!    elseif (jq .eq. jk2) then  ! q is the 2nd argument in v33
!       self = ocfunc(temp,omz,jk1,jk3)
!       sumz = sumz + term * self
!!      call get_k_info(kpc(:,nq2(j)),nc,nq,ik,jk,kk,g1,g2,g3,inside)
!    elseif (jq .eq. jk3) then  ! q is the 3rd argument in v33
!       self = ocfunc(temp,omz,jk2,jk1)
!       sumz = sumz + term * self
!      call get_k_info(kpc(:,nq3(j)),nc,nq,ik,jk,kk,g1,g2,g3,inside)
!    else
!      cycle nv3loop
!    endif
!   write(118,4)' ',term,real(self),aimag(self),sumz
       cntr = cntr+1
       self = ocfunc(temp,omz,la1(j),nq1(j),la3(j),nq3(j))  !jk1,jk3)
       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif
!   endif
 enddo nv3loop
! if (cntr.ne.nkc) then
!    write(ulog,*)' K_WEIGHT INCONSISTENCY IN SELF_NEW: cntr.ne.nkc ',cntr,nkc
! endif
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_new
!============================================================
 subroutine function_self_sc(q,la,omega,temp,nself,uself)
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use ios
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1,cnt2
 integer, save :: cntr
 real(8) om1,eta,tk,v32,term,k1(3),k2(3),q1(3),error
 complex(8) omz,self,oldself,ocfunc

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! om1 = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 write(ulog,5)'SELF_SC: last cntr, omega,etaz=',cntr,omega,etaz
! omz = cmplx(omega,-etaz)   ! this is in cm^-1
 v32 = v3_threshold*v3_threshold
 error=1d9
 self=cmplx(0,etaz)

cnt2=0
do while(error.gt.etaz/50 .and. cnt2.lt.50)

 cnt2=cnt2+1
 nself=0 ; uself=0; oldself=self
 eta = max(etaz, aimag(self))   ! this is in cm^-1
 omz = cmplx(omega,-eta)    ! this is in cm^-1
!omz = cmplx(omega+real(self),-eta)    ! this is in cm^-1
!omz = omega+conjg(self)
 jq = nq+(la-1)*nkc
 cntr = 0
 nv3loop: do j=1,nv3
!  q must be the 2nd argument in v33
    if (nq.ne.nq2(j)) cycle nv3loop
    if (la.ne.la2(j)) cycle nv3loop
    term = v33sq(j) ! v33(j)*conjg(v33(j))
    if (term.lt.v32) cycle nv3loop

    jk1 = nq1(j)+(la1(j)-1)*nkc
    jk2 = nq2(j)+(la2(j)-1)*nkc
    jk3 = nq3(j)+(la3(j)-1)*nkc
    if(jk2 .ne. jq) print*,'ERRRRRRRRRRORRRRRRRRRRE : jk2 ne jq ',jk2,jq

    k2=kpc(:,nq3(j))
    k1 = -q-k2
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,inside)

    if(nk1.ne.nq1(j)) then
      write(ulog,3)'SELF_SC: nk1.ne.nq1 ',nk1,nq1(j)
      write(ulog,3)'Check consistency between nki and those of v33.dat'
      write(ulog,3)'SELF_SC: j,la,q,om=',j,la,q,omega
      write(ulog,4)'SELF_SC: k1,q,k2=',k1,q,k2
      write(ulog,*)'SELF_SC: nqi    =',nq1(j),nq2(j),nq3(j)
      write(ulog,4)'SELF_SC:kpc(nqi)=',kpc(:,nq1(j)),kpc(:,nq2(j)),kpc(:,nq3(j))
      write(ulog,4)'SELF_SC:la,la2(j)=',la,la2(j)
      write(*,3)'SELF_SC: nk1.ne.nq1 ',nk1,nq1(j)
      write(*,3)'SELF_SC: j,la,q,om=',j,la,q,omega
      stop
    endif

       cntr = cntr+1
       self = ocfunc(temp,omz,la1(j),nq1(j),la3(j),nq3(j))  !jk1,jk3)
       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif
 enddo nv3loop
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
 error = abs(self-oldself)
 write(ulog,5)'cnt,error,omz,self=',cnt2,error,omz,self
 if (aimag(self).lt.etaz) exit
enddo
 if (cnt2.ge.50) write(ulog,*)'SC CALCULATION OF THE SELF-ENERGY DID NOT CONVERGE ******************'

3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_sc
!============================================================
 subroutine three_ph_matrix(omega,wdth,res)
! frequencies, width are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! res=sum_q,la sum_12 0.5*|v3(q,la;1;2)|^2 delta(omega-w(q,la))  wk(q)
! it assumes the input momentum q is on the already-generated kmesh
 use kpoints
 use ios
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 real(8), intent(in) :: wdth,omega
 real(8), intent(out):: res
 integer j,nq,jq,la,jk1,jk2,jk3, nqq,nk1,inside
 real(8) omq,term,delta_g,sumz

 res = 0d0
 do nq=1,nibz  !nkc
 do la=1,ndyn
    call get_k_info(kibz(:,nq),nc,nqq,jk1,jk2,jk3,g1,g2,g3,inside)
    omq = sqrt(abs(eigenval(la,nqq))) * cnst
    if (abs(omega-omq).gt.wdth*6) cycle
!   jq = nq+(la-1)*nibz
    jq = nqq+(la-1)*nkc
    term = delta_g(omega-omq,wdth) * wibz(nq) !*wdth*sqrt(2*pi)
    sumz=0
    v3loop: do j=1,nv3
       if (la2(j) .ne. la) cycle v3loop
       if (nq2(j) .ne. nqq) cycle v3loop
       jk1 = nq1(j)+(la1(j)-1)*nkc
       jk2 = nq2(j)+(la2(j)-1)*nkc
       jk3 = nq3(j)+(la3(j)-1)*nkc
       sumz = sumz + term * v33sq(j)  !(v33(j)*conjg(v33(j)))
    enddo v3loop

    res = res+ sumz /nkc /2  *2*pi  ! should be comparable to iself in magnitude
  enddo
  enddo

3 format(i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))

 end subroutine three_ph_matrix
!===========================================================
! function ocfunc(temp,omz,j2,j3)  result(res)
 function ocfunc(temp,omz,l2,nk2,l3,nk3)  result(res)
! for given temperature and complex omz (in 1/cm) calculates the
! occupation-weighted joint DOS stored in res
! j2 and j3 are kpoint indices in the full FBZ
 use eigen
 use constants ! for cnst
 use kpoints   ! for nkc
 use params    ! for classical
 implicit none
 integer , intent(in) :: l2,l3,nk2,nk3
 integer j2,j3
 real(8), intent (in) :: temp
 real(8) om3,om2,nb3,nb2,nbe,resr,resi,et,om,delta_g,delta_l
 complex(8), intent (in) :: omz
 complex(8) res

!! nk2 = 1+mod(j2-1,nkc) ; l2 = 1+ (j2-nk2)/nkc
!! nk3 = 1+mod(j3-1,nkc) ; l3 = 1+ (j3-nk3)/nkc
! call nkla(j2,nkc,nk2,l2)
! call nkla(j3,nkc,nk3,l3)
 j2=(l2-1)*nkc+nk2
 j3=(l3-1)*nkc+nk3

 om3 = sqrt(abs(eigenval(l3,nk3))) * cnst ; nb3=nbe(om3,temp,classical)
 om2 = sqrt(abs(eigenval(l2,nk2))) * cnst ; nb2=nbe(om2,temp,classical)
! om=real(omz) ; et=-aimag(omz)

 if(j2.ne.j3) then   
 !88888888888888888888888888888888888888888888888888888
 ! need a factor of 2 since v3(q,1,2)^2 * (f(1,2)+f(2,1))
 ! the factor of 2 is needed if np<nq (upper half) is excluded in the case where
 ! both sums are done in the full FBZ
 !88888888888888888888888888888888888888888888888888888
!    resi=-pi*((nb2+nb3+1)*(-delta_l(om3+om2-om,et)+ delta_l(om3+om2+om,et)) +  &
! &            (nb2-nb3  )*(-delta_l(om3-om2-om,et)+ delta_l(om3-om2+om,et)) )
!    resr=real(-(nb2+nb3+1)*(1d0/(om3+om2-omz)+ 1d0/(om3+om2+omz))    &
! &            -(nb2-nb3  )*(1d0/(om3-om2-omz)+ 1d0/(om3-om2+omz)) )
     res =     -(nb2+nb3+1)*(1d0/(om3+om2-omz)+ 1d0/(om3+om2+omz))    &
  &            -(nb2-nb3  )*(1d0/(om3-om2-omz)+ 1d0/(om3-om2+omz))
 else
!    resi=-(2*nb2+1) * pi* (-delta_l(2*om2-om,et)+ delta_l(2*om2+om,et))
!    resr=-(2*nb2+1) * real(1d0/(2*om2-omz)+ 1d0/(2*om2+omz))
     res =-(2*nb2+1) *     (1d0/(2*om2-omz)+ 1d0/(2*om2+omz))
 endif
! res=cmplx(resr,resi)

 end function ocfunc
!===========================================================
 function oc2(temp,omz,j2,j3)  result(res)
! for given temperature and complex omz (in 1/cm) calculates the
! occupation-weighted joint DOS stored in res
! j2 and j3 are kpoint indices in the full FBZ
 use eigen
 use constants ! for cnst
 use kpoints   ! for nkc
 use params    ! for classical
 implicit none
 integer , intent(in) :: j2,j3
 integer l2,l3,nk2,nk3
 real(8), intent (in) :: temp
 real(8) om3,om2,nb3,nb2,nbe
 complex(8), intent (in) :: omz
 complex(8) res

! nk2 = 1+mod(j2-1,nkc) ; l2 = 1+ (j2-nk2)/nkc
! nk3 = 1+mod(j3-1,nkc) ; l3 = 1+ (j3-nk3)/nkc
 call nkla(j2,nkc,nk2,l2)
 call nkla(j3,nkc,nk3,l3)

 om3 = sqrt(abs(eigenval(l3,nk3))) * cnst ; nb3=nbe(om3,temp,classical)
 om2 = sqrt(abs(eigenval(l2,nk2))) * cnst ; nb2=nbe(om2,temp,classical)

 if(j2.ne.j3) then   ! need a factor of 2 since v3(q,1,2)^2 * (f(1,2)+f(2,1))
    res = (1d0/(om3-om2-omz)+ 1d0/(om3-om2+omz))
 else
    res = (1d0/(2*om2-omz)+ 1d0/(2*om2+omz))
 endif

 end function oc2
!===========================================================
 subroutine nkla(j,nk,nkj,lj)
! for given one-dimensional index j, and number of kpoints nk, gives the index
! of the actual kpoint, nkj and branch lj
 implicit none
 integer j,nk,nkj,lj

 nkj = 1 + mod(j-1,nk)
 lj  = 1 + (j-nkj)/nk

 end subroutine nkla
!===========================================================
 subroutine get_frequencies(nkp,kp,dk,ndn,eival,nv,eivec,uio,vg)
! also outputs the group velocities from HF theorem (exact) in units of c_light
 use params
 use ios
 use constants
 use born
 use geometry
 implicit none
 integer, intent(in) :: nkp,ndn,uio,nv ! no of wanted eivecs
 real(8), intent(in) :: kp(3,nkp),dk(nkp)
 real(8), intent(out):: eival(ndn,nkp),vg(3,ndn,nkp)
 complex(8), intent(out) :: eivec(ndn,ndn,nkp)
 integer i,j,k,l,ier,nd2,al,ll
 integer, allocatable :: mp(:)
 real(8), allocatable :: eivl(:)
 real(8) absvec,om
 complex(8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:)
 real(8) khat(3)
 character ext*3

! open files and write the group velocities
! do l=1,ndn
!    write(ext,'(i3.3)')l
!    open(uvel+l,file='veloc-'//ext//'.dat')
!    write(uvel+l,*)'# la,i1,i2,i3,nk,kp(nk),eival(l,mapibz(nk)),vg(l,nk),length(vg(l,nk))'
! enddo
 open(uvel,file='veloc.dat')
 write(uvel,*)'# la,nk,kp(nk),eival(l,nk),vg(l,nk),length(vg(l,nk))'
 nd2 = min(ndn,12)
!allocate(dynmat(ndn,ndn),mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3))
 allocate(mp(ndn),eivl(ndn),eivc(ndn,ndn))
! write(uio,*)'# dynmat allocated, printing: Kpoint, eival(j),j=1,ndyn)'

open(800,file='finitvel.dat')

 kloop: do i=1,nkp

!mine
write(800,*)i
 write(*,*)'GET_FREQUENCIES: before setting up dynmat, kp number ',i

! write(uio,*)'############ before setting up dynmat, kp number ',i
    call finitedif_vel(kp(:,i),ndn,vg(:,:,i),eivl,eivc) 

!   call set_dynamical_matrix(kp(:,i),dynmat,ndn,ddyn)
!
!!JS: call nonanalytical term for born effective change
!    if (kp(1,i)==0.d0 .AND. kp(2,i)==0.d0 .AND. kp(3,i)==0.d0) then
!        khat=kp(:,i)+1.0D-10
!    else
!        khat=kp(:,i)
!    endif
!    call nonanal(khat,dynmat,ndn,ddyn)
!
!    if (verbose) then
!       write(ulog,3)' ======================================================================'
!       write(ulog,3)' THE DYNAMICAL MATRIX IS:',i,kp(:,i)
!       do l=1,ndn
!          write(ulog,4)(dynmat(l,j),j=1,nd2)
!       enddo
!    endif
!
!    call diagonalize(ndn,dynmat,eivl,nv,eivc,ier)
!
! sort eivals in ascending order
    call sort(ndn,eivl(:),mp,ndn)

!   if (ier.ne.0 .or. verbose) then
!     write(   *,*)' ier=',ier
!     write(ulog,*)' ier=',ier, ' EIGENVALUES ARE...'
!     write(ulog,*)
!   endif

    do j=1,ndn
       eival(j,i) = eivl(mp(j))
       do l=1,nv
          eivec(j,l,i) = eivc(j,mp(l))
       enddo
    enddo

! if pure imaginary make eigenvectors real as they can be multiplied by any arbitrary number
!    do l=1,nv
!       absvec = sum(abs(real(eivec(:,l,i)))**2)
!       if (absvec .lt. 1d-1) then
!          eivec(:,l,i)=cmplx(0,1)*eivec(:,l,i)   ! this does not hurt anything
!       endif
!    enddo
!
!! put the acoustic bands in the first 3 branches 1:Ta 2:Ta 3:La
!! eivl(mp(i)) is the properly sorted array
!
!!   call sort_polarizations(ndn,eival(:,i),eivec(:,:,i),mp,kp(:,i))
!    do j=1,ndn
!       eivl(j)=eival(mp(j),i)
!       do l=1,nv
!          eivc(j,l) = eivec(j,mp(l),i)
!       enddo
!    enddo
!    call write_eigenvalues(i,dk(i),kp(:,i),eivl,ndn,uio)
!
!    if (verbose) then
!      do l=1,nv
!        write(ulog,3)' -----------  BAND ',l,eivl(l)
!        ll=l !mp(l)
!        do j=1,ndn/3
!           write(ulog,5)'atom ',j,eivec(3*(j-1)+1,ll,i),eivec(3*(j-1)+2,ll,i),eivec(3*(j-1)+3,ll,i)
!!          write(ulog,5)'atom ',j,eivc(3*(j-1)+1,l),eivc(3*(j-1)+2,l),eivc(3*(j-1)+3,l)
!         enddo
!      enddo
!    endif

! now calculate the group velocities from Hellman-Feynmann formula(dynmat=dummy)
!    do al=1,3
!    do l=1,ndn
!      vg(al,l,i)=0
!      do k=1,ndn
!         dynmat(k,l)=0
!         do j=1,ndn
!            dynmat(k,l)=dynmat(k,l)+ddyn(k,j,al)*eivec(j,l,i)
!         enddo
!      enddo
!      do k=1,ndn
!         vg(al,l,i)=vg(al,l,i)+dynmat(k,l)*conjg(eivec(k,l,i))
!      enddo
!      vg(al,l,i)=vg(al,l,i)/2/sqrt(abs(eival(l,i)))*cnst*1d-10*100*2*pi
!    enddo
!    enddo
!
!    do l=1,ndn
!       om=sqrt(abs(eival(l,i)))*cnst
!       write(uvel+l,6)l,i,kp(:,i),om,vg(:,l,i)*c_light,length(vg(:,l,i))*c_light
!    enddo
!
!!   do j=1,ndn
!!      eival(j,i)=eivl(j)
!!      do l=1,nv
!!         eivec(j,l,i)=eivc(j,l)
!!      enddo
!!   enddo
!
!     if(ier .ne. 0) stop
! write(uio,*)'==============='
! write(uio,4)(eivl(j),j=1,ndn)
! write(uio,*)'==============='
! do l=1,ndn
!    write(uio,4)(eivc(l,j),j=1,ndn)
! enddo
    call write_eigenvalues(i,dk(i),kp(:,i),eival(:,i),ndn,uio)
!   write(uio,7)i,kp(:,i),(eivl(j),j=1,ndn)

 enddo kloop
 

 do l=1,ndn
 do i=1,nkp
    om=sqrt(abs(eival(l,i)))*cnst
    write(uvel,6)l,i,kp(:,i),om,vg(:,l,i)*c_light,length(vg(:,l,i))*c_light
 enddo
 enddo
close(800)
! do l=1,ndn
!    close(uvel+l)
! enddo
close(uvel)
 deallocate(eivl,eivc,mp)  !,dynmat,ddyn)
 3 format(a,i5,9(1x,f19.9))
 4 format(99(1x,2(1x,g9.3),1x))
 5 format(a,i5,99(1x,2(1x,g9.3),1x))
 6 format(2i6,2x,99(1x,f11.3))
 7 format(i6,2x,3(1x,f7.3),1x,99(1x,f7.3))

 end subroutine get_frequencies
!==========================================================



!===========================================================
 subroutine get_frequencies2(nkp,kp,dk,ndn,eival,nv,eivec,uio,vg)
! also outputs the group velocities from HF theorem (exact) in units of c_light
 use params
 use ios
 use constants
 use born
 use geometry
 implicit none
 integer, intent(in) :: nkp,ndn,uio,nv ! no of wanted eivecs
 real(8), intent(in) :: kp(3,nkp),dk(nkp)
 real(8), intent(out):: eival(ndn,nkp),vg(3,ndn,nkp)
 complex(8), intent(out) :: eivec(ndn,ndn,nkp)
 integer i,j,k,l,ier,nd2,al,ll
 integer, allocatable :: mp(:)
 real(8), allocatable :: eivl(:)
 real(8) absvec,om
 complex(8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:)
 real(8) khat(3)
 character ext*3

! open files and write the group velocities
! do l=1,ndn
!    write(ext,'(i3.3)')l
!    open(uvel+l,file='veloc-'//ext//'.dat')
!    write(uvel+l,*)'#
!    la,i1,i2,i3,nk,kp(nk),eival(l,mapibz(nk)),vg(l,nk),length(vg(l,nk))'
! enddo
!*** open(uvel,file='veloc.dat')
!*** write(uvel,*)'# la,nk,kp(nk),eival(l,nk),vg(l,nk),length(vg(l,nk))'
 nd2 = min(ndn,12)
!allocate(dynmat(ndn,ndn),mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3))
 allocate(mp(ndn),eivl(ndn),eivc(ndn,ndn))
! write(uio,*)'# dynmat allocated, printing: Kpoint, eival(j),j=1,ndyn)'

!***open(800,file='finitvel.dat')

 kloop: do i=1,nkp

!mine
!***write(800,*)i
!*** write(*,*)'GET_FREQUENCIES: before setting up dynmat, kp number ',i

! write(uio,*)'############ before setting up dynmat, kp number ',i
    call finitedif_vel(kp(:,i),ndn,vg(:,:,i),eivl,eivc)

!   call set_dynamical_matrix(kp(:,i),dynmat,ndn,ddyn)
!
!!JS: call nonanalytical term for born effective change
!    if (kp(1,i)==0.d0 .AND. kp(2,i)==0.d0 .AND. kp(3,i)==0.d0) then
!        khat=kp(:,i)+1.0D-10
!    else
!        khat=kp(:,i)
!    endif
!    call nonanal(khat,dynmat,ndn,ddyn)
!
!    if (verbose) then
!       write(ulog,3)'
!       ======================================================================'
!       write(ulog,3)' THE DYNAMICAL MATRIX IS:',i,kp(:,i)
!       do l=1,ndn
!          write(ulog,4)(dynmat(l,j),j=1,nd2)
!       enddo
!    endif
!
!    call diagonalize(ndn,dynmat,eivl,nv,eivc,ier)
!
! sort eivals in ascending order
    call sort(ndn,eivl(:),mp,ndn)

!   if (ier.ne.0 .or. verbose) then
!     write(   *,*)' ier=',ier
!     write(ulog,*)' ier=',ier, ' EIGENVALUES ARE...'
!     write(ulog,*)
!   endif

    do j=1,ndn
       eival(j,i) = eivl(mp(j))
       do l=1,nv
          eivec(j,l,i) = eivc(j,mp(l))
       enddo
    enddo

! if pure imaginary make eigenvectors real as they can be multiplied by any
! arbitrary number
!    do l=1,nv
!       absvec = sum(abs(real(eivec(:,l,i)))**2)
!       if (absvec .lt. 1d-1) then
!          eivec(:,l,i)=cmplx(0,1)*eivec(:,l,i)   ! this does not hurt anything
!       endif
!    enddo
!
!! put the acoustic bands in the first 3 branches 1:Ta 2:Ta 3:La
!! eivl(mp(i)) is the properly sorted array
!
!!   call sort_polarizations(ndn,eival(:,i),eivec(:,:,i),mp,kp(:,i))
!    do j=1,ndn
!       eivl(j)=eival(mp(j),i)
!       do l=1,nv
!          eivc(j,l) = eivec(j,mp(l),i)
!       enddo
!    enddo
!    call write_eigenvalues(i,dk(i),kp(:,i),eivl,ndn,uio)
!
!    if (verbose) then
!      do l=1,nv
!        write(ulog,3)' -----------  BAND ',l,eivl(l)
!        ll=l !mp(l)
!        do j=1,ndn/3
!           write(ulog,5)'atom
!           ',j,eivec(3*(j-1)+1,ll,i),eivec(3*(j-1)+2,ll,i),eivec(3*(j-1)+3,ll,i)
!!          write(ulog,5)'atom
!',j,eivc(3*(j-1)+1,l),eivc(3*(j-1)+2,l),eivc(3*(j-1)+3,l)
!         enddo
!      enddo
!    endif

! now calculate the group velocities from Hellman-Feynmann formula(dynmat=dummy)
!    do al=1,3
!    do l=1,ndn
!      vg(al,l,i)=0
!      do k=1,ndn
!         dynmat(k,l)=0
!         do j=1,ndn
!            dynmat(k,l)=dynmat(k,l)+ddyn(k,j,al)*eivec(j,l,i)
!         enddo
!      enddo
!      do k=1,ndn
!         vg(al,l,i)=vg(al,l,i)+dynmat(k,l)*conjg(eivec(k,l,i))
!      enddo
!      vg(al,l,i)=vg(al,l,i)/2/sqrt(abs(eival(l,i)))*cnst*1d-10*100*2*pi
!    enddo
!    enddo
!
!    do l=1,ndn
!       om=sqrt(abs(eival(l,i)))*cnst
!       write(uvel+l,6)l,i,kp(:,i),om,vg(:,l,i)*c_light,length(vg(:,l,i))*c_light
!    enddo
!
!!   do j=1,ndn
!!      eival(j,i)=eivl(j)
!!      do l=1,nv
!!         eivec(j,l,i)=eivc(j,l)
!!      enddo
!!   enddo
!
!     if(ier .ne. 0) stop
! write(uio,*)'==============='
! write(uio,4)(eivl(j),j=1,ndn)
! write(uio,*)'==============='
! do l=1,ndn
!    write(uio,4)(eivc(l,j),j=1,ndn)
! enddo
    call write_eigenvalues(i,dk(i),kp(:,i),eival(:,i),ndn,uio)
!   write(uio,7)i,kp(:,i),(eivl(j),j=1,ndn)

 enddo kloop


!*** do l=1,ndn
!*** do i=1,nkp
!***    om=sqrt(abs(eival(l,i)))*cnst
!***    write(uvel,6)l,i,kp(:,i),om,vg(:,l,i)*c_light,length(vg(:,l,i))*c_light
!*** enddo
!*** enddo
!***close(800)
! do l=1,ndn
!    close(uvel+l)
! enddo
!****close(uvel)
 deallocate(eivl,eivc,mp)  !,dynmat,ddyn)
 3 format(a,i5,9(1x,f19.9))
 4 format(99(1x,2(1x,g9.3),1x))
 5 format(a,i5,99(1x,2(1x,g9.3),1x))
 6 format(2i6,2x,99(1x,f11.3))
 7 format(i6,2x,3(1x,f7.3),1x,99(1x,f7.3))

 end subroutine get_frequencies2
!=====================================================

!=====================================================
  subroutine finitedif_vel(q0,ndn,vgr,evl0,evc0)
! calculates the group velocities in units of c_light from finite difference
  use constants
  use params, only : verbose
  implicit none
  integer, intent(in) ::  ndn
  real(8), intent(in) ::  q0(3)
  real(8), intent(out) ::  vgr(3,ndn),evl0(ndn)
  complex(8), intent(out) ::  evc0(ndn,ndn)
  real(8) q1(3),dq,om0,om1
  real(8) evlp(ndn),evlm(ndn)
  complex(8) evct(ndn,ndn)
!mine
integer n,i

  dq=0.001

  call get_freq(q0,ndn,evl0,evc0)
 !write(*,*) evc0
  do i=1,3
     q1=q0; q1(i)=q0(i)+dq
     call get_freq(q1,ndn,evlp,evct)
     q1=q0; q1(i)=q0(i)-dq
     call get_freq(q1,ndn,evlm,evct)
     vgr(i,:)=(evlp-evlm)/2/dq / 2/sqrt(evl0) *cnst*1d-8 *2*pi !*c_light
  enddo
!  write(30,*)'********** q=',q0,'************************'
!  do i=1,3
!     write(30,5)'i,vgr(i,l)=',i,vgr(i,:)
!  enddo
!  write(30,*)'*******************************************'

!mine
if (verbose) then
   do n=1,ndn
      write(800,6)n,q0,vgr(:,n)*c_light
   enddo
endif

 4 format(a,3(1x,f9.3),a)
 5 format(a,i5,99(1x,f9.3))
 6 format(i5,99(1x,g11.4))

  end subroutine finitedif_vel
!============================================================
 subroutine gruneisen(nkp,kp,dk,ndn,eivl,eivc,ugr,grn)
! takes the eigenvalues (w^2) and eigenvectors calculated along some
! crystalline directions and calculates the corresponding mode
! gruneisen parameters
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use eigen
 use constants
 implicit none
 integer ik,i0,nkp,la,al,be,ga,j,k,j0,k0,ta1,ta2,t,ugr,ndn    ! ,i3,j3,k3
 real(8) mi,mj,rr3(3),rr2(3),qq(3),denom,qdotr,omk
 real(8) kp(3,nkp),dk(nkp),eivl(ndn,nkp)
 complex(8) zz,one,term
 complex(8) grn(ndn,nkp), eivc(ndn,ndn,nkp)

 one = cmplx(1d0,0d0)
! write(ugr,*)'# la,nk,dk(nk),kp(:,nk),om(la,nk),gruneisen(la,nk))'
 write(ugr,*)'# nk,dk(nk),kp(:,nk),gruneisen(la,nk).la=1,ndyn)'
 do ik=1,nkp
    qq(:) = kp(:,ik)
 do la=1,ndn
!   write(ulog,6)'Starting gruneisen parameters for kpoint & band# ',ik,la,kp(:,ik)
!   write(ulog,*)' i,la,t,fcs_3(t),term(2),6*eival(la,i),rr3(3),grun(la,i)'
    grn(la,ik) = 0
    denom = 6 * eivl(la,ik)
    omk = sqrt(abs(eivl(la,ik)))*cnst
    do i0=1,natoms0
         mi = atom0(i0)%mass
       tloop: do t=1,nterms(3)

         if ( i0 .ne. iatomterm_3(1,t) ) cycle tloop
         al = ixyzterm_3(1,t)
         be = ixyzterm_3(2,t)
         ga = ixyzterm_3(3,t)
!        i0 = iatomcell0(iatomterm_3(1,t))
         j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ; mj = atom0(j0)%mass
         k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k)
         ta1= al + 3*(i0-1)
         ta2= be + 3*(j0-1)

! rr2(:)=iatomcell(1,j)*r1+iatomcell(2,j)*r2+iatomcell(3,j)*r3
!  be careful: it has to be the translations R not atompos!
         rr2(:) = atompos(:,j) - atompos(:,j0)  ! R
         rr3(:) = atompos(:,k)                  ! R+tau
         qdotr =  ( qq .dot. rr2)
         zz = cdexp( ci * qdotr )
!! term = - ampterm_3(t)*fcs_3(igroup_3(t))*zz*eivc(ta1,la,i)*conjg(eivc(ta2,la,i))*rr3(ga)/sqrt(mi*mj)
         term = - fcs_3(t) * zz * eivc(ta2,la,ik)*conjg(eivc(ta1,la,ik))*rr3(ga)/sqrt(mi*mj)
         grn(la,ik) = grn(la,ik) + term
!        write(ulog,7)i,la,t,fcs_3(t),rr2,qdotr,zz,rr3,grn(la,ik)
       enddo tloop
    enddo
    grn(la,ik) = grn(la,ik)/denom
    if (aimag(grn(la,ik)) .gt. 1d-4) then
       write(ulog,*)' GRUNEISEN: ... has large imaginary part!! for nk,la=',ik,la,grn(la,ik)
!      stop
    endif
!   write(ugr,6)' ',la,ik,dk(ik),kp(:,ik),omk,real(grn(la,ik))
 enddo
    write(ugr,8)ik,dk(ik),kp(:,ik),(real(grn(la,ik)),la=1,ndn)
 enddo
5 format(4i7,9(1x,g10.4))
6 format(a,2i5,99(1x,g10.4))
7 format(i5,i5,i6,99(1x,g10.4))
8 format(i8,99(1x,g10.4))
! deallocate(eivl,eivc,kg,grn)
 end subroutine gruneisen
!============================================================
 subroutine gruneisen_fc
 use phi3
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use eigen
 use constants
 implicit none
 integer i0,al,be,ga,j,k,k0,t,t2,igr2,t2old,grold
 real(8) rr3(3),phi  !mi,mj,

 if( .not. allocated(grun_fc)) allocate(grun_fc(nterms(2)))
 write(ulog,*)'GRUNEISEN_FC: term , phi2_ij , -psi3_ijk , rr(ga) , grun'
 do igr2=1,ngroups(2)
 t2old = 0 ; grold = 0
 do t2=1,nterms(2)
    if( igroup_2(t2) .eq. grold ) cycle
    if( igroup_2(t2) .ne. igr2) cycle
    i0 = iatomcell0(iatomterm_2(1,t2))
    j  = iatomterm_2(2,t2)
    al = ixyzterm_2(1,t2)
    be = ixyzterm_2(2,t2)
!    phi= ampterm_2(t2)*fcs_2(igroup_2(t2))
    phi= fcs_2(t2)
    grun_fc(t2) = 0
!   grn = 0
    tloop: do t=1,nterms(3)
         if ( i0 .ne. iatomterm_3(1,t) ) cycle tloop
         if ( j  .ne. iatomterm_3(2,t) ) cycle tloop
         if ( al .ne. ixyzterm_3(1,t) ) cycle tloop
         if ( be .ne. ixyzterm_3(2,t) ) cycle tloop
         ga = ixyzterm_3(3,t)
         k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k)
         rr3(:) = atompos(:,k)
         grun_fc(t2) = grun_fc(t2) - rr3(ga)  / 6 /phi*fcs_3(t) 
!        grn = grn - rr3(ga)  / 6 /phi*fcs_3(t) !ampterm_3(t) *fcs_3(igroup_3(t))
!         write(ulog,6)'phi3: ', t ,phi,-ampterm_3(t)*fcs_3(igroup_3(t)) , rr3(ga),grn
!         write(ulog,6)'phi3: ', t ,phi,-fcs_3(t) , rr3(ga),grn
    enddo tloop
    grold = igr2
    write(ulog,7)igr2,i0,al,j,be,phi,grun_fc(t2)
 enddo
 enddo

6 format(a,i7,9(1x,g10.4))
7 format('****',i7,2(4x,'(',i3,',',i1,')'),9(1x,g10.4))
 end subroutine gruneisen_fc
!============================================================
 subroutine mechanical(bulk,c11,c44,dlogv)
 use ios
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use constants
 implicit none
 integer i0,al,be,j,t
 real(8) bulk,rija,rijb,c11,c44,dlogv

! write(ulog,*)'BULK_MOD: i0, al, j, be, rija,rijb,fcs_2(t),bulk'
 bulk=0; c11=0 ; c44=0
 do i0=1,natoms0
 do t=1,nterms(2)
    if ( i0 .ne. iatomterm_2(1,t) ) cycle
    al =  ixyzterm_2(1,t)
    be =  ixyzterm_2(2,t)
    j  = iatomterm_2(2,t)
    rija = atompos(al,j)-atompos(al,i0)
    rijb = atompos(be,j)-atompos(be,i0)
!   write(ulog,3)i0, al, j, be, rija,rijb,fcs_2(t),bulk
    bulk = bulk - fcs_2(t)*(1-2*grun_fc(t)*dlogv)*rija*rijb
    if (al.eq.1 .and. be.eq.1) then
       c11 = c11 - fcs_2(t)*(1-2*grun_fc(t)*dlogv)*rija*rijb
    endif
 enddo
 enddo
 bulk = bulk/volume_r/18
 c11 = c11/volume_r/2

 write(ulog,*)'BULK_MODULUS: in eV/A^3 ',bulk
 bulk = bulk*1d30*ee
 write(ulog,*)'BULK_MODULUS: in SI, Mbar ',bulk,bulk*1d-11
 c11 = c11*1d30*ee
 write(ulog,*)'C11 in SI, Mbar ',c11,c11*1d-11
3 format(4(i4),9(2x,f10.5))

 end subroutine mechanical
!============================================================
 subroutine calculate_thermal(nk,wk,ndyn,eival,grn,tmn,tmx,ual)
 use ios
 use lattice
 use constants
 use params
 implicit none
 integer nk,ndyn,b,k,ual,itemp,nat,ntmp,iter
 real(8) wk(nk),eival(ndyn,nk),nbx,cv,cv_nk,cv2,nbe,dlogv,dalpha
 real(8) temp,x,alpha,bulk_modulus,b0,tmn,tmx,gama,a1,c11,c44,etot,free,pres,pres0
 complex(8) grn(ndyn,nk)

! bulk modulus = a_0^2/V d^2E/da^2 evaluated at a_0 (equilibrium lattice parameter)
! call mechanical(b0,c11,c44)  ! this is the bulk_mod at T=0

! B(T,Veq)=B(T=0,Veq(T=0)) - pres0

 nat = ndyn/3

 write(ual,'(a120)')'# temperature(K) ,alpha (1/K) , Cv (J/K/mol), gama , E_tot(J/mol) , &
 &   E_free , Pressure(GPa) , P0 , Bulk_mod(GPa) '

 ntmp=ntemp !60
 do itemp=1,ntmp

    temp=tmn+(tmx-tmn)*(itemp-1)**3/(ntmp-1d0+1d-8)**3   ! this temp is in Kelvin
    if (temp.le.0) then
       write(ulog,*)'temperature not in the proper range!!',temp
       stop
    endif
    alpha=1d2; dalpha=1d9; iter=0
    do while (abs(dalpha) .gt. abs(alpha)/1000 .and. iter .lt. 50)
       iter=iter+1
       dlogv=temp*alpha
       call mechanical(b0,c11,c44,dlogv)  ! this is the bulk_mod at T=0
       call energies(nk,wk,ndyn,eival,grn,temp,etot,free,pres,pres0,cv2)
       bulk_modulus = b0 - pres0

       a1=0 ; cv=0 ; gama = 0
! in the sum over k, all quantities are per unitcell
       do k=1,nk
       do b=1,ndyn
          if(eival(b,k) .lt.0) then
             x=0
             write(ulog,*) 'ALPHA: negative eival for band,kp# ',b,k,eival(b,k)
             write(ulog,*) ' will use its absolute value instead!'
!      else
          endif
! to convert eV/A^2/mass to Hz we need to take the sqrt, multiply by cnst*100*c_light
          x=(h_plank*sqrt(abs(eival(b,k)))*cnst*100*c_light)/k_b/temp
          if (x.gt.60) then
              cv_nk = 0
          else
              nbx=nbe(x,1d0,classical)
              cv_nk = x*x*nbx*(1+nbx) !/4/sinh(x/2)/sinh(x/2)
          endif
          cv = cv + cv_nk*wk(k)   ! this is in J/K units : per unit cell
          a1 = a1 + grn(b,k)*cv_nk*wk(k)
       enddo
       enddo
       gama = a1 / cv
       cv = cv/nat*n_avog*k_b  ! this is per mole    )/(volume_r*1d-30)  !
       if (abs(cv2-cv).gt.1d-5 ) then
          write(ulog,*)'temp, cv from energies ne cv ',temp,cv2,cv
!      stop
       endif
       dalpha = a1/bulk_modulus /n_avog*nat/(volume_r*1d-30) - alpha
       alpha = a1/bulk_modulus /n_avog*nat/(volume_r*1d-30)
       write(*,4)'CALCULATE_THERMAL:i,T,a,da=',iter,temp,alpha,dalpha
    enddo
    write(ual,3)temp,alpha,cv,gama,etot,free,pres*1d-9,pres0*1d-9,bulk_modulus*1d-9
 enddo

3 format(9(2x,g11.5))
4 format(a,i3,9(2x,g11.5))

 end subroutine calculate_thermal
!============================================================
 subroutine read_grun(nk,ndyn,kp,grn)
 use geometry
 implicit none
 integer nk,ndyn,j,i,l
 real(8) grn(ndyn,nk),kp(3,nk),x,y,z  !eival(ndyn,nk),
! logical myeqz

 open(173,file='bulk-grundif-80.dat',status='old')
 do i=1,nk
    read(173,*)j,x,y,z,(grn(l,i),l=1,ndyn)
    if (abs(x-kp(1,i)).gt.0.05 *abs(x)) write(*,*)'kx not OK for i=',i,x,kp(1,i)
    if (abs(y-kp(2,i)).gt.0.05 *abs(y)) write(*,*)'ky not OK for i=',i,y,kp(2,i)
    if (abs(z-kp(3,i)).gt.0.05 *abs(z)) write(*,*)'kz not OK for i=',i,z,kp(3,i)
!    if (.not.(x.myeqz.kp(1,i))) print*,'kx not OK for i=',i,x,kp(1,i)
!    if (.not.(y.myeqz.kp(2,i))) print*,'ky not OK for i=',i,y,kp(2,i)
!    if (.not.(z.myeqz.kp(3,i))) print*,'kz not OK for i=',i,z,kp(3,i)
 enddo
 close(173)

 end subroutine read_grun
!============================================================
 subroutine energies(nk,wk,ndyn,eival,grn,temp,etot,free,pres,pres0,cv)
! calculate total and free energies within QHA, at a given temperature (temp in Kelvin)
 use ios
 use params
 use lattice
 use constants
 implicit none
 integer nk,ndyn,b,k,nat
 real(8) wk(nk),eival(ndyn,nk)
 real(8) temp,x,cv_nk,cv,hw,free,etot,pres,nbe,mdedv,pres0,nbx
 complex(8) grn(ndyn,nk)

    nat = ndyn/3
    if (temp.le.0) then
       write(ulog,*)'temperature not in the proper range!!',temp
       stop
    endif
    etot=0 ; cv=0 ; free=0 ; pres=0
    mdedv= 0.35388/(20.8**3-20.0**3)/ab**3 ! this is -dE/dV in eV/A^3
    mdedv = mdedv*1d+30*ee ! (in J/m^3 )
    do k=1,nk
    do b=1,ndyn
       if(eival(b,k) .lt.0) then
          x=0
          write(ulog,*) 'ALPHA: negative eival for band,kp# ',b,k,eival(b,k)
          write(ulog,*) ' will use its absolute value instead!'
!      else
       endif
! to convert eV/A^2/mass to Hz we need to take the sqrt, multiply by cnst*100*c_light
       x=(h_plank*sqrt(abs(eival(b,k)))*cnst*100*c_light)/k_b/temp
       if (x.gt.60) then
           cv_nk = 0
       else
           nbx=nbe(x,1d0,classical)
           cv_nk = x*x*nbx*(1+nbx) !/4/sinh(x/2)/sinh(x/2)
       endif
       hw = x*k_b*temp  ! hbar*omega in Joules
       cv  = cv + cv_nk*wk(k)   ! this is in J/K units : per unit cell
       etot= etot + hw * (nbe(x,1d0,classical) + 0.5d0) * wk(k)
       free= free + (0.5d0*hw + k_b*temp*log(1-exp(-x))) * wk(k)
       pres= pres + hw * (nbe(x,1d0,classical) + 0.5d0) * wk(k) * grn(b,k)
    enddo
    enddo
    cv = cv/nat*n_avog*k_b  ! this is per mole    )/(volume_r*1d-30)  !
    etot = etot/nat*n_avog   ! convert from Joule/cell to Joule per mole
    free = free/nat*n_avog   ! convert from Joule/cell to Joule per mole
    pres0= pres/(volume_r*1d-30) ! * 1d-8  ! 1d-8 is to convert to kbar
    pres = pres0+ mdedv
3 format(9(2x,g11.5))

 end subroutine energies
!============================================================
 subroutine get_group_velocity(n1,n2,n3,kp,ndyn,nibz,mapibz,eival,vg)
! reads the eigenvalues(frequencies) on a regular mesh (periodic in the FBZ)
! labeled here by kp and calculates the group velocities on this mesh
! output is vg(3,ndyn,nkpoint) in units of c_light
 use constants
 use lattice
 use ios
 implicit none
 integer nk,n1,n2,n3,la,ndyn,i1,i2,i3,ip1,ip2,ip3,im1,im2,im3,nibz
 integer ipp1,ipp2,ipp3,imm1,imm2,imm3
 integer mapibz(n1*n2*n3)
 real(8) kp(3,n1*n2*n3),eival(ndyn,nibz), vg(3,ndyn,n1*n2*n3),v(3),dw(3) ,om
 real(8), allocatable :: omg(:,:,:)
 character ext*2

 do la=1,ndyn
    write(ext,'(i2.2)')la
    open(uvel+la,file='veloc-'//ext//'.dat')
    write(uvel+la,*)'# la,i1,i2,i3,nk,kp(nk),eival(la,mapibz(nk)),vg(la,nk),length(vg(la,nk))'
 enddo
! first sore the eigenvalues on a 3d mesh for each band (PB: bands are sorted
! in increasing order and if they cross, there is a slight pb in v_group
 allocate ( omg(n1,n2,n3) )
 do la=1,ndyn
! first store eivals into a 3d array consistent with kp
    write(ulog,*)'# i1,i2,i3,nk,kp(:,nk),w(la,mapibz(nk)),vg(:,la,nk),|vg| for mode ',la
    nk=0
    do i1=1,n1
    do i2=1,n2
    do i3=1,n3
       nk = nk+1
! below om = 1/lambda in (ang^-1) and v_group=c_light 2pi*d (1/lambda)/dk and dk=2pi/N/R0
       omg(i1,i2,i3) = sqrt(abs(eival(la,mapibz(nk))))*cnst *100*1d-10
    enddo
    enddo
    enddo
! now calculate the derivatives by finite difference
    nk = 0
    do i1=1,n1
    do i2=1,n2
    do i3=1,n3
       nk = nk+1
       ip1 = 1+mod(i1,n1)
       ip2 = 1+mod(i2,n2)
       ip3 = 1+mod(i3,n3)
       im1 = mod(i1+n1-2,n1)+1
       im2 = mod(i2+n2-2,n2)+1
       im3 = mod(i3+n3-2,n3)+1
!      ipp1= 1+mod(ip1,n1)
!      ipp2= 1+mod(ip2,n2)
!      ipp3= 1+mod(ip3,n3)
!      imm1= mod(im1+n1-2,n1)+1
!      imm2= mod(im2+n2-2,n2)+1
!      imm3= mod(im3+n3-2,n3)+1
!      dw(1) = (omg(ip1,i2,i3)-omg(im1,i2,i3) - (omg(ipp1,i2,i3)-omg(imm1,i2,i3))/8)*4/3d0
!      dw(2) = (omg(i1,ip2,i3)-omg(i1,im2,i3) - (omg(i1,ipp2,i3)-omg(i1,imm2,i3))/8)*4/3d0
!      dw(3) = (omg(i1,i2,ip3)-omg(i1,i2,im3) - (omg(i1,i2,ipp3)-omg(i1,i2,imm3))/8)*4/3d0
       dw(1) = omg(ip1,i2,i3)-omg(im1,i2,i3)
       dw(2) = omg(i1,ip2,i3)-omg(i1,im2,i3)
       dw(3) = omg(i1,i2,ip3)-omg(i1,i2,im3)
! 2pi is not needed as 1/cm is for frequency not omega, but v=dw/dk
       vg(:,la,nk)= (n1*dw(1)*r1+n2*dw(2)*r2+n3*dw(3)*r3) / 2d0
!      if ( i1.eq.1) then
          om = sqrt(abs(eival(la,mapibz(nk))))*cnst
          write(uvel+la,3)la,i1,i2,i3,nk,kp(:,nk),om,vg(:,la,nk)*c_light,length(vg(:,la,nk))*c_light
!      endif
    enddo
    enddo
    write(uvel+la,*)' '
    enddo
 enddo

3 format(5(i6),99(1x,g11.5))
4 format(99(1x,g11.5))

 deallocate (omg)
 do la=1,ndyn
    close(uvel+la)
 enddo
 end subroutine get_group_velocity
!============================================================
 subroutine get_k_info(q,N,nk,i,j,k,g1,g2,g3,inside)
! for a vector q(3) in the primitive cell of the reciprocal space, defined on
! a mesh N(1),N(2),N(3), this subroutine finds the three indices of q
! and its number nk based on the triple loop
! nk=0 do i=0,N(1)-1 ; do j=0,N(2)-1; do k=0,N(3)-1; nk=nk+1 ;q=i*G1/N1+j*G2/N2+k*G3/N3

 use geometry
 implicit none
 integer, intent(in):: N(3)
 real(8), intent(in):: q(3)
 type(vector), intent(in):: g1,g2,g3
 integer, intent (out) :: nk,i,j,k,inside
 integer indexg

! this was for when the kpoints were between 0 and N-1
 !  call get_components_g(q,N,i,j,k,g1,g2,g3,inside)
 !  nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
 !  nk= indexg(i,j,k,-N(1)/2,N(1)/2-1,-N(2)/2,N(2)/2-1,-N(3)/2,N(3)/2-1)

! this is for when the k vectors are between -N/2 and N/2-1
  if (mod(N(1),2).eq.0 .and. mod(N(2),2).eq.0 .and. mod(N(3),2).eq.0 ) then
    call get_components_g2(q,N,i,j,k,g1,g2,g3,inside)
!   write(*,*)'for even ni, ijk=',i,j,k
  else   ! do the regulat thing
    call get_components_g(q,N,i,j,k,g1,g2,g3,inside)
  endif
  nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
! write(*,*)'the corresponding nk=',nk

  if (inside.eq.2) then
     write(*,*)' ERROR: vector ',q,' is not a kpoint'
     stop
  endif

 end subroutine get_k_info
!===========================================================
 subroutine collision_matrix(ndyn,nkp,kp,eival,temp,w0,c_matrix)
! this subroutine uses the Fermi-Golden-Rule to get the rate equations
! and construct the linearized collision operator for 3-phonon processes
! to which one can eventually add other scattering meachnisms, and
! diagonalize the matrix to get the rates (inverse of relaxation times)
! as eigenvalues. The latter will be used in Boltzmann equation to
! get the non-equilibruim part of the distribution function
! input is nq= # of q points in the FBZ times # of phonon branches
! temp is in 1/cm
! output is the collision matrix (C_matrix) determined from the rates by using FGR
 use constants
 use params
 implicit none
 integer nq,np,nk1,nk2,nk3,l1,l2,l3,ndyn,nkp,nk,inside
 real(8) c_matrix(nkp*ndyn,nkp*ndyn),ratepk,ratepq,rateqk
 real(8) w0,eival(ndyn,nkp),kp(3,nkp),nbe,temp,nbq,x1,nbp,x2,nbk,x3

 open(321,file='collision_mat.dat')
! write(321,*)'# n , nk , la , c_matrix(nq,nq) , rate(cm^-1) '
 write(321,*)'#  nk , DIAG[c_matrix(nk,la)] la=1,ndyn , in(cm^-1) '

 c_matrix=0

 nq = 0
 do nk1=1,nkp
 do l1=1,ndyn
    nq=nq+1    ! nq is the generic line index
    x1=sqrt(abs(eival(l1,nk1)))*cnst
    nbq = nbe(x1,temp,classical)

    np=0
    do nk2=1,nkp
    do l2=1,ndyn
       np=np+1  ! np is the generic column index
       x2=sqrt(abs(eival(l2,nk2)))*cnst
       nbp = nbe(x2,temp,classical)

       nk=0
       do nk3=1,nkp
       do l3=1,ndyn
          nk=nk+1  ! nk is the generic index over which the sums are made
          x3=sqrt(abs(eival(l3,nk3)))*cnst
          nbk = nbe(x3,temp,classical)

!         call FGR(la1,la2,la3,nk2,nk3,ndyn,nkp,kp,eival,temp,w0,rate1,inside)
!         call FGR(la2,la3,la1,nk3,nk1,ndyn,nkp,kp,eival,temp,w0,rate2,inside)
!         call FGR(la3,la1,la2,nk1,nk2,ndyn,nkp,kp,eival,temp,w0,rate3,inside)
          call FGR(l1,l2,l3,nk2,nk3,ndyn,nkp,kp,eival,temp,w0,ratepk,inside)
          call FGR(l2,l3,l1,nk3,nk1,ndyn,nkp,kp,eival,temp,w0,rateqk,inside)
          call FGR(l3,l1,l2,nk1,nk2,ndyn,nkp,kp,eival,temp,w0,ratepq,inside)

          if (np.eq.nq) then
!            c_matrix(nq,np)= c_matrix(nq,np) +0.5*rate1 +rate2
!            c_matrix(nq,np)= c_matrix(nq,np) +0.5*(rate1+rate2+rate3) ! should produce the same as above line
             c_matrix(nq,nq)= c_matrix(nq,nq) +0.5*ratepk* (1+nbp+nbk) + ratepq*(nbp-nbk)
          else
!            c_matrix(nq,np)= c_matrix(nq,np) - rate1 - rate2 + rate3
             c_matrix(nq,np)= c_matrix(nq,np) + (ratepq+ratepk)*(nbq-nbk) - rateqk*(1+nbq+nbk)
          endif

       enddo
       enddo

    enddo
    enddo

!   rateq=c_matrix(nq,nq)/nbq/(1+nbq)
!   write(321,3) nq,nk1,l1,c_matrix(nq,nq),rateq

 enddo
    nk=(nk1-1)*ndyn
    write(321,3) nk1,(c_matrix(nk+l1,nk+l1),l1=1,ndyn)
 enddo

 close(321)
3 format(1(1x,i9),99(1x,g11.5))

 end subroutine collision_matrix
!===========================================================
 subroutine FGR(l1,l2,l3,nk2,nk3,ndyn,nkp,kp,eival,temp,w0,rate,inside)
! this subroutines uses the FGR to calculate the rates needed in the collision matrix
! input is nq= # of q points in the FBZ times # of phonon branches
! output is rate(1,2,f)= # collision rate to go from state f to states (1,2)
! the delta function width is determined by the average frequency separation as
! delta(x,w0) = w0/pi/(x^2+w0^2) with w0=<w(q+1)-w(q)> =<vgroup> delta_q
! input (besides w0) is s1,s2,sf, integers representing the initial 2 phonon and the
! final one-phonon state : si=1,nk*ndyn
! output is the rate =  2pi/hbar |<i|V_3|j,k>|^2 delta(E_i-E_j-E_k) (all converted to cm^-1)
! to convert Joules to 1/cm need to divide by h_plank*100*c_light
 use constants
 use phi3
 use lattice
 use params
 use kpoints , only : nc
 implicit none
 integer l1,l2,l3,nk1,nk2,nk3,nkp,ndyn,i,j,k,inside
 real(8) w0,eival(ndyn,nkp),kp(3,nkp),nbe,temp,nb1,q(3),x1,rate,delta_l,x2,x3 ,w3!,nb2,nb3


! print*,'FGR: used broadening w0 for this calculation is=',w0,' cm^-1'

 q(:)=-kp(:,nk2)-kp(:,nk3)
 call get_k_info(q,NC,nk1,i,j,k,g1,g2,g3,inside)

 x1=sqrt(abs(eival(l1,nk1)))*cnst  ! convert frequencies to 1/cm as temp is also in 1/cm
 nb1 = nbe(x1,temp,classical)
 x2=sqrt(abs(eival(l2,nk2)))*cnst
! nb2 = nbe(x2,temp,classical)
 x3=sqrt(abs(eival(l3,nk3)))*cnst
! nb3 = nbe(x3,temp,classical)

 x1=(x1-x2-x3)
 w3 = v33sq_5(nk2,nk3,l1,l2,l3)
! rate = 2*pi*delta_l(x1,w0) * w3*conjg(w3) *2*pi * (nb1+1)*nb2*nb3
 rate = 2*pi*delta_l(x1,w0) * w3 *2*pi * (nb1+1)*nb1
 if ( rate.lt.0 ) then
    write(*,*)'FGR: rate<0 ',rate
    stop
 endif
! extra 2pi is needed to convert final results to cm^-1

 end subroutine FGR
!===========================================================
 subroutine write_collision_times(ndyn,nkp,eival,temp,inv_relax_t)
 use constants
 use params
 implicit none
 integer la,nk,nkp,ndyn,nq
 real(8) eival(ndyn,nkp),nbe,temp,inv_relax_t(ndyn*nkp),nbq,x1

 open(322,file='colmat_eivals.dat')
 write(322,*)'# n ,nk,la,eival_coll_matrx(nq) in cm^-1, inverse RT (cm^-1) and Thz'
 nq = 0
 do nk=1,nkp
 do la=1,ndyn
    nq=nq+1    ! nq is the generic line index

    x1=sqrt(abs(eival(la,nk)))*cnst
    nbq = nbe(x1,temp,classical)

    write(322,3) nq,nk,la,inv_relax_t(nq),inv_relax_t(nq)*100*c_light*1d-12
 enddo
 enddo
 close(322)
3 format(3(1x,i9),9(2x,g11.5))

 end subroutine write_collision_times
!===========================================================
 subroutine get_negatives(nk,kp,map,nps)
! this subroutine finds and eliminates the negatives in the kp mesh
! on output: if map(i), (i=1,nps) is the index of the needed kpoints
! the negative k's are kp(:,map(i)), i=nps+1,nk
 use kpoints
 use lattice
 use ios
 use params
 implicit none
 integer nk,i,map(nk),nps,i1,j1,k1,inside,mk,i2,j2,k2,k,match,cnt
 real(8) kp(3,nk) !,q(3)

 map=0; nps=0
 do i=1,nk
    call get_k_info( kp(:,i),NC, k,i1,j1,k1,g1,g2,g3,inside)
    if (k.ne.i) then
       write(ulog,*)'GET_NEGS: error k.ne.i ',k,i
       write(ulog,*)'GET_NEGS: i,j,k =      ',i1,j1,k1
       stop
    endif

    call get_k_info(-kp(:,i),NC,mk,i2,j2,k2,g1,g2,g3,inside)
! if -k is not in the negs list, then add it to the map, else disregard it
    if (mk.lt.i) then ! it is already included
       cycle
    else
       nps=nps+1
       map(nps)=i
    endif
    if(verbose) write(ulog,3)'GET_NEGS: k,-k,map=',i,mk,nps,i1,j1,k1,i2,j2,k2
 enddo
 write(ulog,*)'GET_NEGS: # of included kpoints=',nps
! now end the array map with the negative ones
 cnt=nps
 do i=1,nk
    match=0
    iloop: do i1=1,nps
       if(map(i1).eq.i) then
          match=1
          exit iloop
       endif
    enddo iloop
    if(match.ne.1) then
       cnt=cnt+1
       map(cnt)=i
       if(verbose) write(ulog,3)'negatives: cnt,map=',cnt,i
    endif
 enddo
3 format(a,3i5,2(3x,3i3))
! do i=1,nk
!    write(*,*)'i,map(i)=',i,map(i)
! enddo
 end subroutine get_negatives
!===========================================================
 subroutine subst_eivecs(ndyn,nk,eivec,kp,npos,map)
! this subroutine keeps the eivecs(q) for q in the npos list and replaces the
! other eivec(q') (for q'=-q) by the conjugate of eivec(q) in order to get
! rid of the ambiguity in the case of the degenerate case, and assure e(-q)=conjg(e(q)).
 use lattice  ! needed to access NC,g1,g2,g3
 use ios
 use kpoints , only : nc
 implicit none
 integer nk,ndyn,npos,map(nk),i1,j1,k1,n,inside,mk3,n3,la,mu
 real(8) kp(3,nk) !,q(3)
 complex(8) eivec(ndyn,ndyn,nk),aux(ndyn)

 write(ulog,*)'SUBST_EIVECS: nk,npos=',nk,npos,'=============='
 do n=npos+1,nk
    n3=map(n)      ! actual index of the kpoint
    call get_k_info(-kp(:,n3),NC,mk3,i1,j1,k1,g1,g2,g3,inside)
    do la=1,ndyn
       aux=eivec(:,la,n3)
!      write(ulog,*)'for kpoint,branch=',n3,la
       do mu=1,ndyn
          if(abs(eivec(mu,la,mk3) - conjg(aux(mu))) .gt. 1d-4) then
          if(abs(eivec(mu,la,mk3) + conjg(aux(mu))) .gt. 1d-4) then
!            write(ulog,*)mu,eivec(mu,la,mk3) , conjg(aux(mu))
          endif
          endif
       enddo
       eivec(:,la,mk3) = conjg(aux)
    enddo
 enddo
 write(ulog,*)'EXITING SUBST_EIVECS =========================='

 end subroutine subst_eivecs
!===========================================================
 subroutine check_mdyn(ndyn,nk,kp,eivec,eival)
! form the dynamical matrix from existing eigenvectors and eigenvalues
! and diagonalize and verify consistency :eivec(-q)=-conjg(eivec(q))
 use lattice  ! needed to access NC,g1,g2,g3
 use kpoints , only : nc
 implicit none
 integer ndyn,nk,j,k,l,ik,ier,i1,j1,k1,mi,inside
 real(8) kp(3,nk),eival(ndyn,nk),eivl(ndyn)
 complex(8) eivec(ndyn,ndyn,nk),dynm(ndyn,ndyn),eivc(ndyn,ndyn),d2(ndyn,ndyn)

 do ik=1,nk
    dynm=0
    do j=1,ndyn
    do k=1,ndyn
      do l=1,ndyn
        dynm(j,k)=dynm(j,k)+eival(l,ik)*eivec(j,l,ik)*conjg(eivec(k,l,ik))
      enddo
    enddo
    enddo
    call get_k_info(-kp(:,ik),NC,mi,i1,j1,k1,g1,g2,g3,inside)
    d2=0
    do j=1,ndyn
    do k=1,ndyn
      do l=1,ndyn
        d2(j,k)=d2(j,k)+eival(l,mi)*eivec(j,l,mi)*conjg(eivec(k,l,mi))
      enddo
    enddo
    enddo

    do j=1,ndyn
    do k=1,ndyn
       if(abs( d2(j,k)-conjg(dynm(j,k))) .gt.1d-5) then
         write(*,4)'CHECK_MDYN: j,k,d(-q),d(q)=',j,k,d2(j,k),dynm(j,k)
       endif
    enddo
    enddo
!   return

    call diagonalize(ndyn,dynm,eivl,ndyn,eivc,ier)
    do j=1,ndyn
       if(abs(eivl(j)-eival(j,ik)).gt.1d-4) then
          write(*,3)'CHECK_MDYN:j,eivl,eival=',j,eivl(j),eival(j,ik)
       endif
       do k=1,ndyn
         if(abs(eivc(k,j)-eivec(k,j,ik)).gt.1d-4) then
         if(abs(eivc(k,j)+eivec(k,j,ik)).gt.1d-4) then
            write(*,4)'CHECK_MDYN:j,k,eiv(j),eivecs(k,j)=',j,k,eivl(j),eivc(k,j),eivec(k,j,ik)
         endif
         endif
       enddo
    enddo
 enddo
3 format(a,i6,9(1x,f15.6))
4 format(a,2i6,f15.7,9(1x,f11.4))

 end subroutine check_mdyn
!===========================================================
 subroutine write_eivecs(ndyn,nk,kp,eigenvec)
 use lattice
 use ios
 use geometry
 use constants
 use kpoints , only : nc
 implicit none
 integer nk,ndyn,i,la,j,i2,j2,k2,inside,l,i3,j3,k3
 real(8) kp(3,nk),w(3)
 complex(8) eigenvec(ndyn,ndyn,nk)

 write(ulog,*)'WRITE_EIGVECS: eivc(q)-conjg(eivc(-q)),eivc(q),eivc(-q)'
do i=1,nk
! bring kp(j) into the FBZ

   call get_k_info(kp(:,i),NC,l,i3,j3,k3,g1,g2,g3,inside)
   if (i.ne.l) then
      write(ulog,*) 'i .ne. l ',i,l
      stop
   endif
!  call bring_to_cell_c(-kp(:,i),g1,g2,g3,r1,r2,r3,w)
!  w =w /2/pi
   w = -kp(:,i)
   call get_k_info(w,NC,j,i2,j2,k2,g1,g2,g3,inside)
!  write(ulog,5)'k,-k=',i3,j3,k3,' ** ',i2,j2,k2,kp(:,i),w

  if (j.lt.nk) then
   do la=1,ndyn
   do l=1,ndyn
      if(abs(eigenvec(l,la,i)-conjg(eigenvec(l,la,j))).gt.1d-4) then
      if(abs(eigenvec(l,la,i)+conjg(eigenvec(l,la,j))).gt.1d-4) then
!     write(ulog,8)'la,l,k,-k=', la,l,i,j,eigenvec(l,la,i)-conjg(eigenvec(l,la,j)) &
!     &     ,eigenvec(l,la,i),eigenvec(l,la,j)
      endif
      endif
   enddo
   enddo
  else
    write(ulog,4)'kpoint i in FBZ=',i,kp(:,i)
    write(ulog,4)'was taken to -k=',j,w
  endif
enddo
4 format(a,i5,2x,99(1x,f8.3))
5 format(a,3i5,a,3i5,2x,99(1x,f8.3))
8 format(a,4i6,99(1x,f8.3))

 end subroutine write_eivecs
!===========================================================
 subroutine mode_thermal_conductivity(nk,wk,tau_inv,veloc,omg,temp,kappa_q)
! temp is in cm^-1, veloc in c_light, omg in ev/ang^2/uma, tau_inv in cm_1
! input are the freqs, velocs , relaxation times and weights for each band
! output is the thermal conductivity for that band
 use constants
 use params
 use lattice
 implicit none
 integer, intent(in) ::  nk
 integer i,al,be
 real(8),intent(in):: tau_inv(nk),veloc(3,nk),omg(nk),wk(nk),temp
 real(8), intent(out) :: kappa_q(3,3)
 real(8) cv,x,tau,nbe,nbx

  kappa_q=0 ; cv=0
  do i=1,nk
     x=sqrt(omg(i))*cnst/temp   ! temp is in 1/cm
     if (x.gt.40) then
       cv = 0
     else
       nbx=nbe(x,1d0,classical)
       cv=nbx*(nbx+1)*x*x  !  x*x/4/sinh(x/2)/sinh(x/2)
     endif
     if (tau_inv(i) .eq. 0) then
       tau=1d9
     else
       tau=1/tau_inv(i)
     endif
     if(tau.lt.1d7) then  ! exclude gamma=0
!      kappa_q=kappa_q+cv*wk(i)*tau*sum(veloc(:,i)*veloc(:,i))/3
        do al=1,3
        do be=1,3
       kappa_q(al,be)=kappa_q(al,be)+cv*wk(i)*tau*veloc(al,i)*veloc(be,i)
        enddo
        enddo
     endif
  enddo
  kappa_q = kappa_q *k_b/100*c_light/volume_r*1d30  ! converted to SI (W/mK)

 end subroutine mode_thermal_conductivity
!===========================================================
 subroutine thermal_conductivity(nk,wk,ndn,tau_inv,veloc,omg,temp,kappa_q)
! temp is in cm^-1, veloc in c_light, omg in ev/ang^2/uma, tau_inv in cm_1
 use constants
 use params
 use lattice
 implicit none
 integer, intent(in) ::  nk,ndn
 integer i,la,al,be
 real(8),intent(in):: tau_inv(ndn,nk),veloc(3,ndn,nk),omg(ndn,nk),wk(nk),temp
 real(8), intent(out) :: kappa_q(3,3)
 real(8) cv,x,tau,nbe,nbx

  kappa_q=0 ; cv=0
  do i=1,nk
  do la=1,ndn
     x=sqrt(omg(la,i))*cnst/temp   ! temp is in 1/cm
     if (x.gt.60) then
       cv = 0
     else
       nbx=nbe(x,1d0,classical)
       cv=nbx*(nbx+1)*x*x  ! x*x/4/sinh(x/2)/sinh(x/2)
     endif
     tau=1/tau_inv(la,i)
     if(tau.lt.1d9) then  ! exclude gamma=0
!*********
! should use tau_inv=2*gamma
!*********
        do al=1,3
        do be=1,3
       kappa_q(al,be)=kappa_q(al,be)+cv*wk(i)*tau*veloc(al,la,i)*veloc(be,la,i)
        enddo
        enddo
!      kappa_c=kappa_c+   wk(i)*tau*sum(veloc(:,la,i)*veloc(:,la,i))
     endif
  enddo
  enddo
  kappa_q = kappa_q *k_b/100*c_light/volume_r*1d30  ! converted to SI (W/mK)
! kappa_c = kappa_c/3 *k_b/100*c_light/volume_r*1d30  ! converted to SI (W/mK)

 end subroutine thermal_conductivity
!===========================================================
 subroutine tau_klemens(temp,eival,grn,veloc,wmax,mass,tau_inv)
 use constants
 implicit none
 real(8) temp,eival,grn,veloc(3),wmax,tau_inv,mass,v2,kbt,w2

 v2 = sum(veloc*veloc)*c_light*c_light
 kbt = temp*100*c_light*h_plank  ! kT in Joules
 w2 = eival*cnst*cnst  ! in cm^-2
 tau_inv = kbt/(mass*uma*v2/2) /wmax*w2*grn*grn  ! this is in cm-1
! if (tau_inv .lt. 1d-5) tau_inv=1d30

 end subroutine tau_klemens
!===========================================================
 function eta(om,tempk)
 implicit none
 real(8) om,eta,tempk
 eta = 1d-4*om*om*0.01*(tempk/20)*(1+ 1d-6*om*om*om*0.2)
 end function eta
!===========================================================
 function cross_section(q,la,omega,eival,self) result(sigma)
! calculates the neutron cross section within BORN approx based on the eivals
! of the dynamical matrix and the 3-phonon self energy
 use constants
 implicit none
 real(8) q(3),omega,sigma,xnum,denom,omq,eival
 complex(8) self
 integer la

 omq = cnst*sqrt(eival)
 xnum = 2 *omq*aimag(self)
 denom = xnum*xnum + (omega*omega-omq*omq-2*omq*real(self))**2
 sigma = xnum/denom

 end function cross_section
!===========================================================
 subroutine nonanal(q,dynmat,ndim,ddyn)
 use constants
 use lattice
 use atoms_force_constants
 use born

 integer, intent(in) ::ndim
 complex(8), intent(inout) :: dynmat(ndim,ndim),ddyn(ndim,ndim,3)
 real(8), intent(in) :: q(3)
 real(8) zag(3),zbg(3),qeq,om0,ma,mb,rr(3)
 integer na,nb,i,j
 real(8) eps_scale,q2,gg

 eps_scale=8.8541878176D-12/1D10/ee
 gg=(length(g1)*length(g2)*length(g3))**0.33
! if ( gg .lt. 4*rho) then
! if (born_flag.eq.0)  rho = gg/4d0
  rho = gg/4d0
! write(30,*)'ADOPTED VALUE of RHO = ',rho
! endif
 rho2=rho*rho

 call calculate_volume(r1,r2,r3,om0)

 qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+ &
&       q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+ &
&       q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
 q2=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)

 do na = 1,natoms0
   ma = atom0(na)%mass
 do nb = 1,natoms0
   mb = atom0(nb)%mass
   rr = atompos(:,na)-atompos(:,nb)
   do i=1,3
!     zag(i) = q(1)*zeu(1,i,na)+q(2)*zeu(2,i,na)+q(3)*zeu(3,i,na)
!     zbg(i) = q(1)*zeu(1,i,nb)+q(2)*zeu(2,i,nb)+q(3)*zeu(3,i,nb)
      zag(i) = q(1)*atom0(na)%charge(1,i)+ &
      &        q(2)*atom0(na)%charge(2,i)+ &
      &        q(3)*atom0(na)%charge(3,i)
      zbg(i) = q(1)*atom0(nb)%charge(1,i)+ &
      &        q(2)*atom0(nb)%charge(2,i)+ &
      &        q(3)*atom0(nb)%charge(3,i)

   end do
   do i = 1,3
   do j = 1,3
!write(*,*)
!dynmat(i+3*(na-1),j+3*(nb-1)),
!4*pi*zag(i)*zbg(j)/(qeq*eps_scale)/om0/sqrt(ma*mb)*exp(-q2/rho2),ma,mb,zag(i),zbg(j)
!4*pi*zag(i)*zbg(j)/(qeq*eps_scale)/om0/sqrt(ma*mb)*exp(-q2/rho2)
      dynmat(i+3*(na-1),j+3*(nb-1)) = dynmat(i+3*(na-1),j+3*(nb-1))+ &
      & zag(i)*zbg(j)/(qeq*eps_scale)/om0/sqrt(ma*mb)*exp(-q2/rho2)
  !! ALSO NEED TO CORRECT THE GROUP VELOCITIES
   !  ddyn(i+3*(na-1),j+3*(nb-1)) = ddyn(i+3*(na-1),j+3*(nb-1))+ &
   !  & zag(i)*zbg(j)/(qeq*eps_scale)/om0/sqrt(ma*mb)*exp(-q2/rho2)
   end do
   end do
 end do
 end do

 end subroutine nonanal
!=========================================================
 subroutine read_born
! reads ! born ! effective ! charge ! and ! dielectric ! constants ! from ! para.born
 use born
 use lattice
 use atoms_force_constants
 use ios
 integer i,j,k

 open(uborn,file='dielectric.params',status='old')
 read(uborn,*) born_flag   ! if 0 use default
 do i=1,3
    read(uborn,*)(epsil(i,j),j=1,3)
 end do
 do k=1,natoms0
    do i=1,3
       read(uborn,*)atom0(k)%charge(i,:)
!       read(uborn,*)(zeu(i,j,k),j=1,3)
    end do
 end do
 close(uborn)

 end subroutine read_born
!============================================================
 subroutine matrix_elt(q1,q2,l1,l2,l3,w33,inside)
! for input modes (qi,li) i=1,2,3 calculates the 3-phonon matrix element w33(1,2,3) on the fly
 use kpoints
 use ios
 use phi3
 use svd_stuff
 use eigen
 use params
 use constants
 use atoms_force_constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: l1,l2,l3
 real(8), intent(in) :: q1(3),q2(3)
 integer, intent(out) :: inside
 complex(8), intent (out) :: w33
 integer i3,nk1,nk2,nk3,i0,al,be,ga,j,j0,k,k0,t,ta1,ta2,ta3,i1,j1,k1
 real(8) q3(3),mi,mj,mk,rr2(3),rr3(3),den
 complex(8) xx,eiqr

! write(*,*)'matrix_elt',q1,q2,l1,l2,l3
! write(*,*)'matrix_elt, const33=',const33
 call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
 call get_k_info(q2,NC,nk2,i1,j1,k1,g1,g2,g3,inside)
 w33 = cmplx(0d0,0d0)
 q3 = -q2-q1
 call get_k_info(q3,NC,nk3,i1,j1,k1,g1,g2,g3,inside)
! write(*,*)'matrix_elt: q3=',q3,nk3,inside
 xx = cmplx(0d0,0d0)

 tloop: do t=1,nterms(3)

! write(*,*)'matrix_elt: tloop=',t,nterms(3)
       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop ! first index must be in primitive cell
       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rr2(:) = atompos(:,j) - atompos(:,j0)
       rr3(:) = atompos(:,k) - atompos(:,k0)
       eiqr = exp( ci* ((q2 .dot. rr2) + (q3 .dot. rr3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
&      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)

 enddo tloop
 den = sqrt(8*sqrt(abs(eigenval(l1,nk1)*eigenval(l2,nk2)*eigenval(l3,nk3))))
 if (den.ne.0) then
    w33 = xx / den * const33
 else
    write(*,*)'MATRIX_ELT: den=0 ',den
    stop
 endif

 end subroutine matrix_elt
!===========================================================
 subroutine matrix_elt_full(q,q2,q3,omq,om2,om3,evq,ev2,ev3,w33)
! for input modes (qi,li) i=1,2,3 calculates the 3-phonon matrix element w33(1,2,3) on the fly
 use kpoints
 use ios
 use phi3
 use svd_stuff
 use eigen
 use params
 use constants
 use atoms_force_constants
 use om_dos
 use lattice
 implicit none
! integer, intent(in) :: la
 real(8), intent(in) :: q(3),q2(3),q3(3),omq,om2,om3
 complex(8), intent(in) :: evq(ndyn),ev2(ndyn),ev3(ndyn)
! integer, intent(out) :: inside
 complex(8), intent (out) :: w33
 integer i3,nk1,nk2,nk3,al,be,ga,j1,k1,i0,j0,k0,t,ta1,ta2,ta3
 real(8) mi,mj,mk,rr2(3),rr3(3),den
 complex(8) xx,eiqr

! write(*,*)'matrix_elt',q1,q2,l1,l2,l3
! write(*,*)'matrix_elt, const33=',const33
 w33 = cmplx(0d0,0d0)
 xx = cmplx(0d0,0d0)

 tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop ! first index must be in primitive cell
       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j1 = iatomterm_3(2,t)  ;  j0 = iatomcell0(j1) ;  mj = atom0(j0)%mass
       k1 = iatomterm_3(3,t)  ;  k0 = iatomcell0(k1) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rr2(:) = atompos(:,j1) - atompos(:,j0)
       rr3(:) = atompos(:,k1) - atompos(:,k0)
       eiqr = exp( ci* ((q2.dot. rr2) + (q3 .dot. rr3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
&      evq(ta1)*ev2(ta2)*ev3(ta3)
! &      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)

 enddo tloop
 den = sqrt(8*omq*om2*om3/(cnst*cnst*cnst))
 if (den.ne.0) then
    w33 = xx / den * const33
 else
    write(*,*)'MATRIX_ELT: den=0 ',den
    stop
 endif

 end subroutine matrix_elt_full
!===========================================================
 subroutine calculate_w3sq_ibz(ndn,ni,ki,nk,kp,eival,eivec)
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units for V3 are in cm^-1
! this works for two k's (k2,k3) in the prim cell on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! this subroutine calculates V3(k1,la1,k2,la2,k3,la3)  with k1=-k2-k3+G  and all ki in the primitive G-cell
! with the restriction of k2 in the IRREDUCIBLE BZ defined by ki(3,ni) included in kp(3,nk)
!
 use kpoints
 use ios
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ni,ndn,l1,l2,l3,n2,n3,j,k
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,indx
 real(8) q1(3),q2(3),q3(3),eival(ndn,nk)
 real(8) ki(3,ni),kp(3,nk)
 complex(8) eivec(ndn,ndn,nk),xx

 v33sq = 0d0 !cmplx(0d0,0d0)
 indx=0
 write(ulog,*)'V3: nc(123)=',nc

!************************************************
! second argument q2 needs to be only in the IFBZ
!************************************************
 loop2: do n2=1,ni
    q2 = ki(:,n2)  ! should be = kp(:,mapinv(n2))
    call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(mapinv(n2).ne.nk2) then
      write(ulog,*)'n2,mapinv(n2),nk2,inside=',n2,mapinv(n2),nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif

 loop3: do n3=1 ,nk
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif

    q1 = -q2-q3

    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)

  write(*,2)'nk123=',nk1,nk2,nk3

 do l2=1 ,ndn
 do l3=1 ,ndn
 do l1=1 ,ndn
    call matrix_elt(q2,q3,l2,l3,l1,xx,inside)
    indx=indx+1
    v33sq(indx)=xx*conjg(xx)
    nq1(indx)= nk1
    nq2(indx)= nk2
    nq3(indx)= nk3
    la1(indx)= l1
    la2(indx)= l2
    la3(indx)= l3
!   if (mod(indx,1000).eq.1) write(*,*) 'indx,v33=',indx,xx
 enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop3
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop2

 nv3 = indx
 write(ulog,*)' V33: total size of this array is =',nv3

 if( writev3.eq.1) then
   open(uv3,file='v33.dat',status='unknown',FORM='UNFORMATTED')
   write(uv3)nv3
!  open(uv3,file='v33.dat',status='unknown')
!  write(uv3,*)nv3
   do j=1 ,nv3
      write(uv3) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33sq(j)
   enddo
   close(uv3)
 endif

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))

 end subroutine calculate_w3sq_ibz
!============================================================
 subroutine function_self_w(q,la,omega,temp,nself,uself)
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use ios
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1,nk3,ins2,ifbz,l1,l3
 integer, save :: cntr
 real(8) om1,eta,tk,v32,term,k1(3),k2(3),q1(3),q3(3)
 complex(8) omz,self,ocfunc,oc2,s2,xx

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
! write(*,*) 'self_w: tempk=',tk
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! write(*,*) 'self_w: k_info: nq=',nq
! om1 = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz 
! write(ulog,5)'SELF_W: last cntr, omega,etaz=',cntr,omega,etaz
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0
 v32 = v3_threshold*v3_threshold

 cntr = 0
 do ifbz=1,nkc
    q3=kpc(:,ifbz)
    call get_k_info(q3,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
    if (nk3.ne.ifbz) then
       write(ulog,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       write(*,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       stop
    endif
    k1 = -q-q3
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,ins2)
!   if (ins2.ne.inside) then
!      write(ulog,*)'SELF_W: ins2 ne inside ',ins2,inside
!      stop
!   endif

    do l3=1,ndyn
    do l1=1,ndyn
!      write(*,*) 'self_w: enetering matrix_elt',q,q3,la,l3,l1
       call matrix_elt(q,q3,la,l3,l1,xx,inside)
       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)

       jk1 = nk1+(l1-1)*nkc
       jk3 = ifbz+(l3-1)*nkc

!      cntr = cntr+1
       self = ocfunc(temp,omz,l1,nk1,l3,nk3) !jk1,jk3)
 !     s2 = oc2(temp,omz,jk1,jk3)
 !     if (aimag(s2).gt.1/6/eta) write(222,6)ifbz,q3,q,la,l3,l1
       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif
    enddo
    enddo

 enddo
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
6 format(i6,6(1x,f9.4),2x,3(1x,i2))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w
!============================================================
 subroutine function_self_w2(q,la,omega,temp,nself,uself)
! calculates the self-energy on the fly for q-point in the generated kmesh and Lorentzian delta
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use ios
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1,nk3,ins2,ifbz,l1,l3
 integer, save :: cntr
 real(8) omq,om3,omk,eta,tk,v32,term,k1(3),k2(3),q1(3),q3(3)
 complex(8) omz,self,ocfunc,oc2,s2,xx

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
! write(*,*) 'self_w: tempk=',tk
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! write(*,*) 'self_w: k_info: nq=',nq
 omq = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz
! write(ulog,5)'SELF_W: last cntr, omega,etaz=',cntr,omega,etaz
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0
 v32 = v3_threshold*v3_threshold

 cntr = 0
 do ifbz=1,nkc
    q3=kpc(:,ifbz)
    call get_k_info(q3,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
    if (nk3.ne.ifbz) then
       write(ulog,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       write(*,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       stop
    endif
    k1 = -q-q3
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,ins2)

    do l3=1,ndyn
       om3 = sqrt(abs(eigenval(l3,nk3))) * cnst
    do l1=1,ndyn
       omk = sqrt(abs(eigenval(l1,nk1))) * cnst
!      call matrix_elt     (q,q3,la,l3,l1,xx,inside)
       call matrix_elt_full(q,q3,k1,omq,om3,omk,eigenvec(:,la,nq),eigenvec(:,l3,nk3),eigenvec(:,l1,nk1),xx)
       call check_inside_bz(k1,g1,g2,g3,inside)
       if(inside.ne.ins2) then
          write(*,*)'SELF_W:ERROR:inside .ne. ins2 ',inside,ins2,k1
          stop
       endif
       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)

       jk1 = nk1+(l1-1)*nkc
       jk3 = ifbz+(l3-1)*nkc

       self = ocfunc(temp,omz,l1,nk1,l3,nk3)
       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif
    enddo
    enddo

 enddo
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
6 format(i6,6(1x,f9.4),2x,3(1x,i2))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w2
!============================================================
 subroutine function_self_w3(q,la,omega,temp,nself,uself)
! this one is for arbitrary q
! uses lorentzian delta with 4 terms, and calculates both real and imaginary parts
! In this subroutine, which calculates like self_w on the fly, the second momentum
! is restricted to be on the kpc mesh, but the first and third momenta can be anywhere
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use ios
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer inside,ik,l3,l2
! integer, save :: cntr
 real(8) omq,eta,term,k3(3),k2(3),eivq(ndyn),eiv2(ndyn),eiv3(ndyn)
 real(8) etacut,arg,delta_l,om3,om2,nb3,nb2,nbe,ires,rres
 complex(8) omz,self,xx,evq(ndyn,ndyn),ev2(ndyn,ndyn),ev3(ndyn,ndyn)

! be careful band sorting should be avoided in this case; otherwise need to
! find the nearest kpoint in kpc, and sort eivq according to it
 call get_freq(q,ndyn,eivq,evq)
 omq=sqrt(abs(eivq(la)))*cnst
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0
 etacut = eta*3000 ! error due to non-zero eta is 1/9E6

 do ik=1,nkc
    k2=kpc(:,ik)
    ev2=eigenvec(:,:,ik)
    k3 = -q-k2
    call check_inside_bz(k3,g1,g2,g3,inside)

    call get_freq(k3,ndyn,eiv3,ev3)

    do l2=1,ndyn
       om2=sqrt(abs(eigenval(l2,ik)))*cnst
       nb2=nbe(om2,temp,classical)
    do l3=1,ndyn

! first check see if the energy is conserved
       om3=sqrt(abs(eiv3(l3)))*cnst
       nb3=nbe(om3,temp,classical)
       ires=0; rres=0

       arg=(omega-om3-om2)
       if (abs(arg) .lt. etacut) then
          ires=(nb3+nb2+1)*delta_l(arg,eta)
          rres=(nb3+nb2+1)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega+om3+om2)
       if (abs(arg) .lt. etacut) then
          ires=ires-(nb3+nb2+1)*delta_l(arg,eta)
          rres=rres-(nb3+nb2+1)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega-om3+om2)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb2-nb3)*delta_l(arg,eta)
          rres=rres+(nb2-nb3)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega-om2+om3)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb3-nb2)*delta_l(arg,eta)
          rres=rres+(nb3-nb2)*arg/(arg*arg+eta*eta)
       endif

       self=cmplx(rres,ires*pi)
 !     if (abs(self).lt.1d-6) cycle loop1   ! self is in cm

       call matrix_elt_full(q,k2,k3,omq,om2,om3,evq(:,la),ev2(:,l2),ev3(:,l3),xx)

       if (abs(xx).lt.v3_threshold) cycle 
       term =xx*conjg(xx)

       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif

 enddo
 enddo
 enddo
 nself = nself /nkc *pi
 uself = uself /nkc *pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w3
!============================================================
 subroutine function_self_w4(q,la,omega,temp,nself,uself)
! this one is for arbitrary q
! uses gaussian delta with 3 terms, and just calculates the imaginary part.
! In this subroutine, which calculates like self_w on the fly, the second momentum
! is restricted to be on the kpc mesh, but the first and third momenta can be anywhere
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use ios
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer inside,ik,l3,l2
! integer, save :: cntr
 real(8) omq,eta,term,k3(3),k2(3),eivq(ndyn),eiv2(ndyn),eiv3(ndyn)
 real(8) etacut,arg,delta_g,om3,om2,nb3,nb2,nbe,ires,rres
 complex(8) omz,self,xx,evq(ndyn,ndyn),ev2(ndyn,ndyn),ev3(ndyn,ndyn)

! be careful band sorting should be avoided in this case; otherwise need to
! find the nearest kpoint in kpc, and sort eivq according to it
 call get_freq(q,ndyn,eivq,evq)
 omq=sqrt(abs(eivq(la)))*cnst
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0
 etacut = eta*6 ! error due to non-zero eta is 1/9E6

 do ik=1,nkc
    k2=kpc(:,ik)
    ev2=eigenvec(:,:,ik)
    k3 = -q-k2
    call check_inside_bz(k3,g1,g2,g3,inside)

    call get_freq(k3,ndyn,eiv3,ev3)

    do l2=1,ndyn
       om2=sqrt(abs(eigenval(l2,ik)))*cnst
       nb2=nbe(om2,temp,classical)
    do l3=1,ndyn

! first check see if the energy is conserved
       om3=sqrt(abs(eiv3(l3)))*cnst
       nb3=nbe(om3,temp,classical)
       ires=0; rres=0

       arg=(omega-om3-om2)
       if (abs(arg) .lt. etacut) then
          ires=(nb3+nb2+1)*delta_g(arg,eta)
!         rres=(nb3+nb2+1)*arg/(arg*arg+eta*eta)
       endif

! this affects the lifetime of low-frequency phonons, can be removed...
! contributes only to noremal processes, and might yield a higher power of w!
!       arg=(omega+om3+om2)   
!       if (abs(arg) .lt. etacut) then
!          ires=ires-(nb3+nb2+1)*delta_g(arg,eta)
!!         rres=rres-(nb3+nb2+1)*arg/(arg*arg+eta*eta)
!       endif

       arg=(omega-om3+om2)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb2-nb3)*delta_g(arg,eta)
!         rres=rres+(nb2-nb3)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega-om2+om3)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb3-nb2)*delta_g(arg,eta)
!         rres=rres+(nb3-nb2)*arg/(arg*arg+eta*eta)
       endif

       self=cmplx(rres,ires*pi)
 !     if (abs(self).lt.1d-6) cycle loop1   ! self is in cm

       call matrix_elt_full(q,k2,k3,omq,om2,om3,evq(:,la),ev2(:,l2),ev3(:,l3),xx)

       if (abs(xx).lt.v3_threshold) cycle 
       term =xx*conjg(xx)

       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif

 enddo
 enddo
 enddo
 nself = nself /nkc *pi
 uself = uself /nkc *pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w4
!============================================================
 subroutine phase(ndn,ni,ki,nk,kp,eival,eivec)
! this subroutine computes for each q the sum_q2 delta(w-w_q2 \pm w_q3)
! and finds the kpoints available
 use kpoints
 use ios
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ni,ndn,l1,l2,l3,n2,n3,j,k
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk3,indx,igam
 real(8) q1(3),q2(3),q3(3),eival(ndn,nk)
 real(8) ki(3,ni),kp(3,nk)
 complex(8) eivec(ndn,ndn,nk),xx

 v33sq = 0d0 !cmplx(0d0,0d0)
 indx=0
 write(ulog,*)'V3: nc(123)=',nc

    q1=0
    call get_kindex(q1,nk,kp,igam)
    q1 = kp(:,igam+3)

 loop3: do n3=1 ,nk
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif

!! properly initilize q2 !!!!!!!!
    q1 = -q2-q3

    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)

  write(*,2)'nk13=',nk1,nk3

! do l2=1 ,ndn
! do l3=1 ,ndn
! do l1=1 ,ndn
!    call matrix_elt(q2,q3,l2,l3,l1,xx,inside)
! enddo
! enddo
! enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop3
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)

 nv3 = indx
 write(ulog,*)' V33: total size of this array is =',nv3

! if( writev3.eq.1) then
!   open(uv3,file='v33.dat',status='unknown',FORM='UNFORMATTED')
!   write(uv3)nv3
!   do j=1 ,nv3
!      write(uv3) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33sq(j)
!   enddo
!   close(uv3)
! endif

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))

 end subroutine phase
!===========================================================
 subroutine get_freq(kp,ndn,eival,eivec)
 use params
 use ios
 use constants
 use born
 implicit none
 integer, intent(in) :: ndn
 real(8), intent(in) :: kp(3)
 real(8), intent(out):: eival(ndn)
 complex(8), intent(out) :: eivec(ndn,ndn)
 integer i,j,k,l,ier,nd2,al
 integer, allocatable :: mp(:)
 real(8), allocatable :: eivl(:)
 real(8) absvec,vg(3,ndn)
 complex(8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:),temp(:,:)
 real(8) khat(3)

 nd2 = min(ndn,12)
 allocate(dynmat(ndn,ndn),mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3),temp(ndn,ndn))
! write(uio,*)'# dynmat allocated, printing: Kpoint, eival(j),j=1,ndyn)'

    call set_dynamical_matrix(kp,dynmat,ndn,ddyn)

!JS: call nonanalytical term for born effective change
    if (kp(1)==0.d0 .AND. kp(2)==0.d0 .AND. kp(3)==0.d0) then
        khat=kp(:)+1.0D-10
    else
        khat=kp(:)
    endif
    call nonanal(khat,dynmat,ndn,ddyn)

    if (verbose) then
       write(ulog,3)' ======================================================================'
       write(ulog,3)' THE DYNAMICAL MATRIX for nd , KP=',nd2,kp(:)
       do l=1,ndn
          write(ulog,8)(dynmat(l,j),j=1,nd2)
       enddo
    endif

! store dynmat in temp
    temp=dynmat
    call diagonalize(ndn,dynmat,eivl,ndn,eivc,ier)

! sort eivals in ascending order
    call sort(ndn,eivl,mp,ndn)
!   write(ulog,6)'map=',mp

    if (ier.ne.0 ) then
      write(   *,*)' ier=',ier
      write(ulog,*)' ier=',ier, ' EIGENVALUES ARE...'
      write(ulog,4) eivl
    endif

    do j=1,ndn
       eival(j) = abs(eivl(mp(j)))  ! so that all frequencies are positive
       do l=1,ndn
          eivec(l,j) = eivc(l,mp(j))
       enddo
    enddo

! if pure imaginary make eigenvectors real as they can be multiplied by any arbitrary number
    do l=1,ndn
       absvec = sum(abs(real(eivec(:,l))))
       if (absvec .lt. 1d-3) then
          eivec(:,l)=cmplx(0,1)*eivec(:,l)
       endif
    enddo

  if (verbose) then
    do l=1,ndn
        write(ulog,3)' GET_FREQ:-----------  BAND ',l,eival(l)
        do j=1,ndn/3
           write(ulog,5)'atom ',j,eivec(3*(j-1)+1,l),eivec(3*(j-1)+2,l),eivec(3*(j-1)+3,l)
        enddo
    enddo
  endif

! now calculate the square of frequencies based on dynmat
!   dynmat=temp
!   temp = matmul(dynmat,eivc)
!   temp = matmul(transpose(conjg(eivc)),temp)

!    do j=1,ndn
!      write(ulog,9)'e.D.e=',j, temp(j,:)

!     do k=1,ndn
!        eivl(k)=0
!        do l=1,ndn
!           eivl(k)=eivl(k)+dynmat(k,l)*eivec(l,j)
!        enddo
!     enddo
!     absvec=0  ! dummy for w^2
!     do k=1,ndn
!        absvec=absvec+ eivl(k)*conjg(eivec(k,j))
!     enddo
!     write(ulog,3)'j,om^2(j)=',j,absvec,eival(j)

!    enddo

! now calculate the group velocities from Hellman-Feynmann formula(dynmat=dummy)
!   do al=1,3
!   do l=1,ndn
!     vg(al,l)=0
!     do k=1,ndn
!        dynmat(k,l)=0
!        do j=1,ndn
!           dynmat(k,l)=dynmat(k,l)+ddyn(k,j,al)*eivec(j,l)
!        enddo
!     enddo
!     do k=1,ndn
!        vg(al,l)=vg(al,l)+dynmat(k,l)*conjg(eivec(k,l))
!     enddo
!     vg(al,l)=vg(al,l)/2/sqrt(abs(eival(l)))*cnst*1d-10*100*2*pi
!   enddo
!   enddo

 deallocate(eivl,eivc,dynmat,mp,ddyn,temp)

 3 format(a,i5,9(1x,f13.4))
 4 format(99(1x,2(1x,g9.3),1x))
 5 format(a,i5,99(1x,2(1x,g9.3),1x))
 6 format(a,99(1x,i5))
 7 format(i6,2x,3(1x,f7.3),1x,99(1x,f7.3))
 8 format(99(1x,2(1x,f9.3),1x))
 9 format(a,i5,99(1x,2(1x,f9.3),1x))

 end subroutine get_freq
!===========================================================
!subroutine sort_polarizations(n,nk,eivl,eivc,mp,q)
! put all longitudinal modes first, then transverse ones in the
! mapping array mp so that mp(1,...,n/3) contains the indices of
! the longitudinal branches and mp(n/3+1,...,n) has those of transverse modes

!use ios
!implicit none
!integer n,nk,i,j,l,mp(n),iless,igreat,ndeg
!real(8) eivl(n,nk),q(3,nk),x,y,z,rr,ri,qq,dpi,dpr,mean(n),dp,sq12,tmp
!complex(8) eivc(n,n,nk),mat(n,n)
!real(8), allocatable:: mm(:)
!integer, allocatable:: m2(:),aux(:)

!sq12=1/sqrt(2d0)
!allocate(aux(n))

!do l=2,nk
! compute eivec(l)*eivec(l+1)
!   mat=0
!   do i=1,n
!      do j=1,n
!         mat(i,j)=eivc(j,i,l-1)*conjg(eivc(j,i,l))
!         write(ulog,6)'kpoint#',l,i,j,mat(i,j)
!      enddo
!   enddo
!   aux=eivl(:,l) 
!   do i=1,n
!      if(abs(mat(i,j)*conjg(mat(i,j))).gt.0.5/n) then
!        ! switch i and j
!        tmp=eivl(i,l)
!        eivl(i,l)=eivl(j,l)
!        eivl(j,l)=tmp
!      endif
!   enddo
!enddo
! qq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
!4 format(a,99(1x,g10.4))
!5 format(a,i6,99(1x,g10.4))
!6 format(a,3i6,3(1x,2(1x,f5.2)),3x,99(1x,g10.4))
!7 format(a,i6,3(1x,1(1x,f5.2)),3x,99(1x,g10.4))
!write(ulog,4)'Q=',q,qq
!do i=1,n
!   mp(i)=i
!enddo
!if (qq.eq.0) return

! cntl=0; cntt=0
!do i=1,n    ! for each eigenstate
!   mean(i)=0
!   do j=1,n/3
!      x=real( eivc((j-1)*3+1,i))
!      y=real (eivc((j-1)*3+2,i))
!      z=real (eivc((j-1)*3+3,i))
!      rr=sqrt(x*x+y*y+z*z+1d-10)
!      dpr=(q(1)*x+q(2)*y+q(3)*z)/qq /rr
!      x=aimag(eivc((j-1)*3+1,i))
!      y=aimag(eivc((j-1)*3+2,i))
!      z=aimag(eivc((j-1)*3+3,i))
!      ri=sqrt(x*x+y*y+z*z+1d-10)
!      dpi=(q(1)*x+q(2)*y+q(3)*z)/qq /ri
!  !   if (rr.gt.ri ) then
!         dp= dpr
!  !   else
!  !      dp= dpi
!  !   endif
!      mean(i)= mean(i)+dp*dp  !abs(dp)
!  !   write(ulog,6)'jatom,eivec,rr,ri,dp,mean=',j,(eivc((j-1)*3+l,i),l=1,3),rr,ri,dp,mean(i)
!      write(ulog,7)'jatom,eivec,rr,dp,mean=',j,(real(eivc((j-1)*3+l,i)),l=1,3),rr,dp,mean(i)
!   enddo
!   mean(i)=mean(i)*3d0/n
!   write(ulog,*)'mode#, mean dotproduct=',i,mean(i)
!enddo

!call sort(n,-mean,mp,n)

!allocate(mm(n))
!do i=1,n
!   mm(i)=-mean(mp(i))
!enddo

!dp=abs( mm(n/3)-mm(n/3+1) )
! if (dp.gt.0.01) then
!   do i=1,n
!      write(ulog,*)'index, mode assignment=',i,mp(i)
!   enddo
!   deallocate(mm)
!   return
! else
! find the nearest non-zero gaps around N/3, and in this range, restore the old order
!   do i=0,n/3-2
!      dp=abs( mm(n/3-i)-mm(n/3-i-1) )
!      if (dp.gt.0.01) exit
!   enddo
!   iless=i
!   do i=0,2*n/3-2
!      dp=abs( mm(n/3+1+i)-mm(n/3+2+i) )
!      if (dp.gt.0.01) exit
!   enddo
!   igreat=i

!   deallocate(mm)

! there is zero gap between mm(n/3-iless) and mm(n/3+1+igreat)
! their number is 1+igreat-iless

!   ndeg=igreat+iless+2
!   write(ulog,*)'GAP:i<,i>,ndeg=',iless,igreat,ndeg
!   allocate(m2(ndeg),mm(ndeg),aux(ndeg))
!   do i=1,ndeg
!      mm(i)=dfloat(mp(n/3-iless+i-1))
!      write(ulog,*)'i,mp(n/3-iless+i-1)=',i,mp(n/3-iless+i-1)
!   enddo
!   call sort(ndeg,mm,m2,ndeg)
!   do i=1,ndeg
!      aux(i)=mp(n/3-iless+m2(i)-1)
!   enddo
!   do i=1,ndeg
!      mp(n/3-iless+i-1)=aux(i)
!      write(ulog,*)'i,mp(n/3-iless+i-1)=',i,mp(n/3-iless+i-1)
!   enddo
!   do i=1,n
!      write(ulog,*)'index, mode assignment=',i,mp(i)
!   enddo
!   deallocate(mm,aux,m2)

! endif

!end subroutine sort_polarizations
!===========================================================
!subroutine sort_polarizations2(n,eivl,eivc,mp,q)
! sort according to E(k).E(k+dk)

!use ios
!implicit none
!integer n,i,j,l,mp(n),iless,igreat,ndeg
!real(8) eivl(n),q(3),x,y,z,rr,ri,qq,dpi,dpr,mean(n),dp,sq12
!complex(8) eivc(n,n)
!real(8), allocatable:: mm(:)
!integer, allocatable:: m2(:),aux(:)

!sq12=1/sqrt(2d0)
!qq=sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
!4 format(a,99(1x,g10.4))
!5 format(a,i6,99(1x,g10.4))
!6 format(a,i6,3(1x,2(1x,f5.2)),3x,99(1x,g10.4))
!7 format(a,i6,3(1x,1(1x,f5.2)),3x,99(1x,g10.4))
!write(ulog,4)'Q=',q,qq
!do i=1,n
!   mp(i)=i
!enddo
!if (qq.eq.0) return

! cntl=0; cntt=0
!do i=1,n    ! for each eigenstate
!   mean(i)=0
!   do j=1,n/3
!      x=real( eivc((j-1)*3+1,i))
!      y=real (eivc((j-1)*3+2,i))
!      z=real (eivc((j-1)*3+3,i))
!      rr=sqrt(x*x+y*y+z*z+1d-10)
!      dpr=(q(1)*x+q(2)*y+q(3)*z)/qq /rr
!      x=aimag(eivc((j-1)*3+1,i))
!      y=aimag(eivc((j-1)*3+2,i))
!      z=aimag(eivc((j-1)*3+3,i))
!      ri=sqrt(x*x+y*y+z*z+1d-10)
!      dpi=(q(1)*x+q(2)*y+q(3)*z)/qq /ri
!  !   if (rr.gt.ri ) then
!         dp= dpr
!  !   else
!  !      dp= dpi
!  !   endif
!      mean(i)= mean(i)+dp*dp  !abs(dp)
!  !   write(ulog,6)'jatom,eivec,rr,ri,dp,mean=',j,(eivc((j-1)*3+l,i),l=1,3),rr,ri,dp,mean(i)
!      write(ulog,7)'jatom,eivec,rr,dp,mean=',j,(real(eivc((j-1)*3+l,i)),l=1,3),rr,dp,mean(i)
!   enddo
!   mean(i)=mean(i)*3d0/n
!   write(ulog,*)'mode#, mean dotproduct=',i,mean(i)
!enddo

!call sort(n,-mean,mp,n)

!allocate(mm(n))
!do i=1,n
!   mm(i)=-mean(mp(i))
!enddo

!dp=abs( mm(n/3)-mm(n/3+1) )
! if (dp.gt.0.01) then
!   do i=1,n
!      write(ulog,*)'index, mode assignment=',i,mp(i)
!   enddo
!   deallocate(mm)
!   return
! else
! find the nearest non-zero gaps around N/3, and in this range, restore the old order
!   do i=0,n/3-2
!      dp=abs( mm(n/3-i)-mm(n/3-i-1) )
!      if (dp.gt.0.01) exit
!   enddo
!   iless=i
!   do i=0,2*n/3-2
!      dp=abs( mm(n/3+1+i)-mm(n/3+2+i) )
!      if (dp.gt.0.01) exit
!   enddo
!   igreat=i

!   deallocate(mm)

! there is zero gap between mm(n/3-iless) and mm(n/3+1+igreat)
! their number is 1+igreat-iless

!   ndeg=igreat+iless+2
!   write(ulog,*)'GAP:i<,i>,ndeg=',iless,igreat,ndeg
!   allocate(m2(ndeg),mm(ndeg),aux(ndeg))
!   do i=1,ndeg
!      mm(i)=dfloat(mp(n/3-iless+i-1))
!      write(ulog,*)'i,mp(n/3-iless+i-1)=',i,mp(n/3-iless+i-1)
!   enddo
!   call sort(ndeg,mm,m2,ndeg)
!   do i=1,ndeg
!      aux(i)=mp(n/3-iless+m2(i)-1)
!   enddo
!   do i=1,ndeg
!      mp(n/3-iless+i-1)=aux(i)
!      write(ulog,*)'i,mp(n/3-iless+i-1)=',i,mp(n/3-iless+i-1)
!   enddo
!   do i=1,n
!      write(ulog,*)'index, mode assignment=',i,mp(i)
!   enddo
!   deallocate(mm,aux,m2)

! endif

!end subroutine sort_polarizations2
!===========================================================
 subroutine enforce_asr_simple(n,dyn)
 implicit none
 integer n,nat,i,j,io,jo,lo,iat,jat
 complex(8) dyn(n,n),sum
 
 nat=n/3
 do i=1,3
 do iat=1,nat
    io=i+3*(iat-1)
    do j=1,3
       sum=0
       jo=j+3*(iat-1)
       do jat=1,nat
          if (iat.eq.jat) cycle
          lo=j+3*(jat-1)
          sum=sum-dyn(io,lo)
       enddo
       dyn(io,jo)=sum
    enddo
 enddo
 enddo
 end subroutine enforce_asr_simple
