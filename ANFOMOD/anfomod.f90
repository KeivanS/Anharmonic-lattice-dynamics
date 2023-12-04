!============================================================
  module md_params
  real(8) pi,temperature_k,jj_avg(3),te,ke,current(3),potential, &
  &       dt,damp_rate,amplitude,t_unit,w_unit,timer,temperature,  &
  &       potential0,msda,relax_time,damping_time,omega,sigma,eps,rcut2,phi
  real(8), allocatable:: forc(:,:),disp(:,:),velo(:,:),cur(:,:)

  integer timestep,nstep,ene_period,xv_period,current_period,seed,  &
  &       nensemble,nfft_t,countr,potflag,relaxation_time,bc
  logical relax
  integer, parameter :: uene=60, ucurrent=50
  end module md_params
!============================================================
 program anharmonic_md
! performs MD simulation of an up to 4th order anharmonic potential
! use fc init and the modules and subroutines for initialization as before
! with the least number of changes.
! size of the super cell must be more than twice the range of interactions
! energies are in eV, lengths in Ang, mass in 1.66d-27 kg, time in 10.1786 fs
 use md_params
 use constants
 use params
 use atoms_force_constants
 use svd_stuff
 use lattice
 use ios
 implicit none
 character xt*1,fn*55,time*8
 integer iens,rank,iunit,mx,mxshl
 logical :: ex
 real(8) dt2,zero
 real tim,tim1,tim2

 zero=0
 call read_parameters ! needed for running the MD: dt,nsteps,nensembles...

 if (potflag .eq. 0) then
    open(ulog,file='log_tl.dat',status='unknown')
    open(uene,file='ene_tl.dat',status='unknown')
    open(ucurrent,file='cur_tl.dat',status='unknown')
    open(utraj,file='traj_tl.xyz',status='unknown')
  else
    open(ulog,file='log_lj.dat',status='unknown')
    open(uene,file='ene_lj.dat',status='unknown')
    open(ucurrent,file='cur_lj.dat',status='unknown')
    open(utraj,file='traj_lj.xyz',status='unknown')
 endif

 call cpu_time(tim1)
 write(ulog,'(a,f9.4)')' STARTING TIME OF THE ANFOMOD PROGRAM IS ',tim1

 write(ulog,*)'# of atom types=', natom_type
 write(ulog,*)'# of ensembles =', nensemble
 write(ulog,*)'relaxation time=', relaxation_time
 write(ulog,*)'delta t for MD =', dt
 write(ulog,*)'avg-temperature=', temperature_k
 write(ulog,*)'disp amplitude =', amplitude
 write(ulog,*)'# of FFT steps =', nfft_t
 write(ulog,*)'potential flag =', potflag
 write(ulog,*)'tolerance for= =', tolerance  ! for equating two coordinates
 write(ulog,*)'en,xv,jQ period=', ene_period,xv_period,current_period  ! output periods
 write(ulog,*)'seed, damp_rate=', seed, damp_rate

 write(uene,'(a108)')'#  step     time     Kinetic/N      Potential/N  &
&    TotalE/N     cumul_Temperature     cumul_current(3)    cumul_MSD (Ang)'
 write(ucurrent,'(a)')'# nstep,timer,jj0,sum(velo(1,:)),sum(velo(2,:)),sum(velo(3,:))'
! primitive lattice and its basis atoms, 
 call read_lattice

 call read_allfcs(potflag)

 call cpu_time(tim)
 write(ulog,'(a,f9.4)')' TIME after reading the FCs IS ',tim

 call read_supercell_coordinates  ! cartesian in Ang
 call make_r0g
 call check_super_cell_coordinates  ! consistency with primitive cell read from lat_fc.dat

 allocate( disp(3,natom_super_cell),forc(3,natom_super_cell),velo(3,natom_super_cell),cur(3,natom_super_cell) )
! if (potflag.ne.0) allocate( forc0(3,natom_super_cell,natom_super_cell) ) 

! calculate reference potential energy of atoms at equilibrium, which will be subtracted from actual potential
 if (potflag .ne. 0) call calculate_potential_lj0
 write(ulog,*)'After potential_lj0, phi=',phi  ! used to define amplitude

 nstep = 0
 dt2 = dt/2d0
ensemble_avge: do iens=1,nensemble

     msda=0 ; countr=0 ; jj_avg=0 ; temperature=0
 call cpu_time(tim)
 write(ulog,'(a,f11.4,i9)')'************* TIME , ensemble number =',tim,iens

! generate initial velocity and displacement of gaussian distribution
! for displacements, amplitude depends on on-site FC, phi, which we arbitrarily set phi to 10 eV/Ang^2
     amplitude=sqrt(k_b*temperature_k/ee/abs(phi))  ! in Ang 
     call move_atoms_gaussian(seed+20*iens,natom_super_cell,amplitude,disp,velo)
     call calculate_potential
     call updatev(zero) ! dt=0 just calculate ke/atom; no velocity update
     call write_output_and_average(current,ke,potential)

! the flag relax outs in contact with thermostat; if false, microcanonical run
     relax = .True.
     EQUILIBRATION_LOOP: do timestep=1,relaxation_time
        call updatexv
        call apply_bc
        call calculate_potential
        call updatev(dt2)
        te = potential + ke
     enddo EQUILIBRATION_LOOP

     call reset_momentum(natom_super_cell,velo) ! due to adding random forces  Ptot\=0
     call rescale_velocity(natom_super_cell,velo,temperature_k)

     relax = .False.
     temperature = 0 ; countr = 0 ; jj_avg = 0 ; msda = 0
     MD_LOOP: do timestep=1,nfft_t*current_period
        nstep = nstep + 1   ! cumulative number of steps (includes ens avg)
        timer = nstep*dt
        call updatexv
        call apply_bc
        call calculate_potential
        call updatev(dt2)
        call write_output_and_average(current,ke,potential)
     enddo MD_LOOP

!    write(ulog,*) nstep,' time steps written into jk and vk, loop=',iens

 enddo ensemble_avge

 call cpu_time(tim2)
 write(ulog,'(a,f11.4,i9)')' END TIME , ensemble number =',tim2,iens-1
 write(ulog,'(a,f14.4)')' TOTAL RUN TIME =',tim2-tim1
 write(*   ,'(a,f14.4)')' TOTAL RUN TIME =',tim2-tim1

 close(utraj)
 close(ulog)
 close(uene)
 close(ucurrent)

 3 format(i6,9(2x,g11.5))
 4 format(a,9(2x,g11.5))

 end program anharmonic_md
!============================================================
 subroutine calculate_potential
! for given displacements and force constants, this subroutine
! constructs the potential and the force function iand the current based
! on the equilibrium positions of atoms and superlattice specifications
! size of the super cell must be more than twice the range of interactions
 use OMP_LIB
 use lattice
 use md_params
 use geometry
 use atoms_force_constants
 use svd_stuff
 implicit none
 real(8) const,transl(3),random(3,natom_super_cell),fij,zero,phi2
 integer t,i,j,k,l,i0,j0,k0,l0,ni(3),nj(3),nk(3),nl(3),taui,tauj,tauk,taul
 integer al,be,ga,de
 real(8), allocatable, save:: f_ij(:,:),r_ij(:,:,:) 
! integer , allocatable :: jatom(:)
! if potential=phi u u /2 + psi u u u /6 + chi u u u u /24 then
! current = Rij vi uj (phi_ij /2 + psi u /3 + chi u u /8)
! kappa  = Volume/3 k_BT^2 \int dt <j(t)j(0)>
! n_atom = size(disp(3,:))
! write(ulog,*)' size of disp=number of atoms=',n_atom
! allocate(f_ij(natom_super_cell,natoms),r_ij(3,natom_super_cell,natoms) ) !,jatom(natoms))

 if (potflag.ne.0) then
    call calculate_potential_lj
    return
 endif

! otherwise potflag=0 and we use the polynomial potential read from fci.dat
! allocate only once!
 if(.not.allocated(f_ij)) allocate(f_ij(natom_super_cell,natoms),r_ij(3,natom_super_cell,natoms)) 

 f_ij = 0 ; r_ij = 0
 potential = 0  !potential0  !-1.435911d0*natom_super_cell
 forc = 0; current = 0; cur = 0; zero = 0
!$OMP PARALLEL DO !I is private by default
 I_LOOP: do i=1,natom_super_cell
   taui = atom_sc(i)%cell%tau
   ni   = atom_sc(i)%cell%n
!
! RANK 1 ------------------
    do t=1,nterms(1)
! this is the relative (wrt central prim cell) index of j
       i0   = iatomterm_1(1,t)
       if ( i0 .eq. taui ) then
          al = ixyzterm_1(1,t)
!         forc(al,i) = forc(al,i) - fcs_1(igroup_1(t))*ampterm_1(t)
!         potential = potential +  fcs_1(igroup_1(t))*ampterm_1(t)*disp(al,i)
          forc(al,i) = forc(al,i) - fcs_1(t)
          potential = potential +  fcs_1(t)*disp(al,i)
       endif
    enddo
! RANK 2 ------------------
    do t=1,nterms(2)
       i0   = iatomterm_2(1,t)
       if (i0 .eq. taui ) then   ! select the term such that i is in the primitive cell
          al   = ixyzterm_2(1,t)
          be   = ixyzterm_2(2,t)
          j0   = iatomterm_2(2,t)   ! this is the relative (wrt central prim cell) index of j
          tauj = atom_shl(j0)%cell%tau
          nj(:)= atom_shl(j0)%cell%n + ni(:) ! translate by ni
          call findatom_sc(nj,tauj,j)
!         phi2 =  fcs_2(igroup_2(t))*ampterm_2(t)
!         phi2 =  fcs_2(t)
          fij =  fcs_2(t) * disp(be,j)  ! if -disp(be,i)) !then take i.ne. j
          forc(al,i) = forc(al,i) - fij
          potential = potential + fij * disp(al,i)/2
          f_ij(i,j0) = f_ij(i,j0) - fij * velo(al,i)
!          f_ij(i,j0) = f_ij(i,j0) - phi2 * &
!&         (disp(al,i)*vel(be,j)+vel(al,i)*disp(be,j)) /4
  r_ij(:,i,j0) = atom_shl(j0)%equilibrium_pos-atom_shl(i0)%equilibrium_pos
!          write(ulog,8)i,j0,j,t,al,be, fij, r_ij(:,i,j0)
!         write(ulog,9)i,j0, r_ij(:,i,j0),f_ij(:,i,j0)
       endif
    enddo
! RANK 3 ------------------
    do t=1,nterms(3)
       i0   = iatomterm_3(1,t)
       if (i0 .eq. taui ) then
          al   = ixyzterm_3(1,t)
          be   = ixyzterm_3(2,t)
          ga   = ixyzterm_3(3,t)
          j0   = iatomterm_3(2,t)
          tauj = atom_shl(j0)%cell%tau
          nj(:)= atom_shl(j0)%cell%n + ni(:) ! translate by ni
          call findatom_sc(nj,tauj,j)
          k0   = iatomterm_3(3,t)
          tauk = atom_shl(k0)%cell%tau
          nk(:)= atom_shl(k0)%cell%n + ni(:) ! translate by ni
          call findatom_sc(nk,tauk,k)
!         fij = fcs_3(igroup_3(t))*ampterm_3(t) * disp(be,j)*disp(ga,k)/2
          fij = fcs_3(t) * disp(be,j)*disp(ga,k)/2
          forc(al,i) = forc(al,i) - fij
          potential = potential + fij * disp(al,i)/3
          f_ij(i,j0) = f_ij(i,j0) - fij * velo(al,i) *1/3d0
       endif
    enddo
! RANK 4 ------------------
    do t=1,nterms(4)
       i0   = iatomterm_4(1,t)
       if (i0 .eq. taui ) then
          al   = ixyzterm_4(1,t)
          be   = ixyzterm_4(2,t)
          ga   = ixyzterm_4(3,t)
          de   = ixyzterm_4(4,t)
          j0   = iatomterm_4(2,t)
          tauj = atom_shl(j0)%cell%tau
          nj(:)= atom_shl(j0)%cell%n + ni(:) ! translate by ni
          call findatom_sc(nj,tauj,j)
          k0   = iatomterm_4(3,t)
          tauk = atom_shl(k0)%cell%tau
          nk(:)= atom_shl(k0)%cell%n + ni(:) ! translate by ni
          call findatom_sc(nk,tauk,k)
          l0   = iatomterm_4(4,t)
          taul = atom_shl(l0)%cell%tau
          nl(:)= atom_shl(l0)%cell%n + ni(:) ! translate by ni
          call findatom_sc(nl,taul,l)
!         fij=fcs_4(igroup_4(t))*ampterm_4(t)*disp(be,j)*disp(ga,k)*disp(de,l)/6
          fij=fcs_4(t)*disp(be,j)*disp(ga,k)*disp(de,l)/6
          forc(al,i) = forc(al,i) - fij
          potential = potential + fij * disp(al,i)/4
          f_ij(i,j0) = f_ij(i,j0) - fij * velo(al,i) *3/2d0
       endif
    enddo

! current is calculated here: current(i) = 1/volume sum'_j (Fij.Vi) Rij /2
    do j0=1,natoms
          transl = r_ij(:,i,j0)
          if (length(transl) .myeq. zero) cycle
!          write(ulog,5)'i,j,rij,f.v=',i,j0, length(transl), transl,f_ij(i,j0)
          cur(:,i) = cur(:,i) - f_ij(i,j0) * transl/2
    enddo
    current(:) = current(:) + cur(:,i)
!    write(ucurrent,7)i,cur(:,i)
 enddo I_LOOP
!$OMP END PARALLEL DO

 if (relax) then ! apply heatbath contact
    const = sqrt(6*damp_rate*temperature/dt)
    call random_number(random); random = (2*random-1)*const
    do i=1,natom_super_cell
    do al=1,3
       forc(al,i)=(forc(al,i)+random(al,i)*sqrt(atom_sc(i)%mass)-atom_sc(i)%mass * &
&                      damp_rate*velo(al,i) ) / (1+dt*damp_rate/2)
    enddo
    enddo
 endif

 potential = potential/natom_super_cell
 current   = current  /omega
5 format(a,i6,i6,4(1x,f6.3),9(1x,g9.3))
8 format(6(1x,i4),1(1x,g9.3),1x,9(1x,f5.2))
9 format(2(1x,i4),3(1x,f5.2),3(1x,g9.3))
7 format((1x,i4),3(1x,g10.4))

! deallocate(f_ij,r_ij)

 end subroutine calculate_potential
!====================================================================
 subroutine read_parameters
 use md_params
 use params
 use atoms_force_constants
 use ios
 use constants, only : k_b,ee
 implicit none

 open(uparams,file='anh_md.params',status='old')

 read(uparams,*) natom_type
 read(uparams,*) nensemble
 read(uparams,*) relaxation_time ! in fs
 read(uparams,*) dt          ! in fs
 read(uparams,*) temperature_k
 read(uparams,*) amplitude   ! initial atomic amplitudes in angstroem
 read(uparams,*) nfft_t      ! number of times data is recorded every current_period steps
 read(uparams,*) potflag     ! 0 for taylor expansion
 read(uparams,*) tolerance   ! for equating two coordinates
 read(uparams,*) bc          ! 1 for boundary conditions 
 read(uparams,*) ene_period,xv_period,current_period  ! output periods
 read(uparams,*) seed, damping_time  ! relaxation time in same units as dt
 damp_rate = 1d0/damping_time     !damp_rate=inverse of relaxation time

 write(*,*)' READ_PARAMETERS: reading done!'
 allocate(natom(natom_type))
 close(uparams)
 temperature=k_b*temperature_k/ee ! convert Kelvin to eV

 end subroutine read_parameters
!====================================================================
 subroutine read_lattice
! reads the translation vectors of the primitive cell, its basis (or atoms)
! the atoms in its neighborhood consistent with Harold's code...
 use lattice
 use params
 use md_params
 use svd_stuff
 use atoms_force_constants
 use geometry
 use ios
 implicit none
 character line*80
 integer i,j

 open(ufc,file='lat_fc.dat',status='old')
 read(ufc,*)line
 read(ufc,*)r01
 read(ufc,*)r02
 read(ufc,*)r03
 read(ufc,*)line
 read(ufc,*)natom_prim_cell
 allocate(atom0(natom_prim_cell))
 do j=1,natom_prim_cell
    read(ufc,*)i,atom0(i)%name,atom0(i)%at_type,atom0(i)%equilibrium_pos,atom0(i)%mass
    atom0(i)%tau=i
 enddo
 read(ufc,*)line
 write(ulog,'(a)')line
 read(ufc,*)line
 write(ulog,'(a)')line
 read(ufc,*)include_fc(1),include_fc(2),include_fc(3),include_fc(4)
 write(ulog,*)include_fc(1),include_fc(2),include_fc(3),include_fc(4)
 read(ufc,*)line
 write(ulog,'(a)')line
 read(ufc,*)nterms(1),nterms(2),nterms(3),nterms(4)
 write(ulog,*)nterms(1),nterms(2),nterms(3),nterms(4)
 read(ufc,*)line
 write(ulog,'(a)')line
 read(ufc,*)ngroups(1),ngroups(2),ngroups(3),ngroups(4)  ! not relevant in this code
 write(ulog,*)ngroups(1),ngroups(2),ngroups(3),ngroups(4)
 read(ufc,*)line
 write(ulog,'(a)')line
 read(ufc,*) maxshell,natoms
 write(ulog,*) maxshell,natoms
 allocate(atom_shl(natoms))
 write(ulog,*)'******** READ NEIGHBOR-SHELL COORDINATES ARE *********'
 write(ulog,*)' i ,    equilibrium positions        ,         tau , n1 , n2 , n3 '
 do i=1,natoms
    read(ufc,*)j,atom_shl(j)%equilibrium_pos,atom_shl(j)%cell%tau,atom_shl(j)%cell%n
    write(ulog,4)j,atom_shl(j)%equilibrium_pos,atom_shl(j)%cell%tau,atom_shl(j)%cell%n
 enddo

 call make_reciprocal_lattice(r01,r02,r03,g01,g02,g03)

 close(ufc)
4 format(i6,4(1x,g10.4),4(2x,i3))

 end subroutine read_lattice
!====================================================================
 subroutine updatexv
 use atoms_force_constants
 use md_params
 implicit none
 integer i

 do i=1,natom_super_cell
    velo(:,i) = velo(:,i) + dt*forc(:,i)/atom_sc(i)%mass/2d0
    disp(:,i) = disp(:,i) + dt*velo(:,i)
 enddo

 end subroutine updatexv
!====================================================================
 subroutine updatev(delta)
 use atoms_force_constants
 use md_params
 implicit none
 real(8) :: delta
 integer i

 ke = 0
 do i=1,natom_super_cell
    velo(:,i) = velo(:,i) + delta*forc(:,i)/atom_sc(i)%mass
    ke = ke + atom_sc(i)%mass * sum(velo(:,i)*velo(:,i))
 enddo
 ke = ke /2d0/natom_super_cell

 end subroutine updatev
!====================================================================
 subroutine write_output_and_average(jj0,kine,pe)
! outputs the energies and coordinates and the heat current every n timesteps
 use atoms_force_constants
 use md_params
 use ios
 implicit none
 integer :: j
 real(8), intent(in) :: jj0(3),pe,kine
 real(8) tote
 character nam*2

 countr = countr + 1  ! cumulative steps for that ensemble
 temperature = temperature + kine*2d0/3  ! *2/space dimension
 tote=kine+pe
 jj_avg = jj_avg + jj0
 msda = msda + sqrt( sum(disp*disp)/(3d0*natom_super_cell) )


!  =============================
 if(mod(timestep-1,ene_period).eq.0) then


    write(*,4)    timestep,timer,kine,pe,tote,temperature/countr,jj_avg/countr
    write(uene,3) timestep,timer,kine,pe,tote,temperature/countr,jj_avg/countr,msda/countr

 endif
!  =============================
 if(mod(timestep-1,xv_period).eq.0) then

    write(utraj,4) natom_super_cell
    write(utraj,4) timestep,timer,kine,pe,tote,temperature/countr,jj_avg/countr
    do j=1,natom_super_cell
  !    write(utraj,6)j,atom_sc(j)%equilibrium_pos+disp(:,j),forc(:,j)!,cur(:,j)
       nam=atom0(atom_sc(j)%cell%tau)%name
       write(utraj,7)nam ,atom_sc(j)%equilibrium_pos+disp(:,j),forc(:,j),cur(:,j)
    enddo

 endif
!  =============================
 if(mod(timestep-1,current_period).eq.0) then

    nstep = nstep+1
    write(ucurrent,3)nstep,timer,jj0,sum(velo(1,:)),sum(velo(2,:)),sum(velo(3,:))

 endif
!  =============================

3 format(i9,9(1x,g13.7))
4 format(i9,9(1x,g11.5))
5 format(999(1x,f11.6))
6 format(i8,3(1x,f9.5),2x,999(1x,g9.3))
7 format(a,3(1x,f9.5),2x,999(1x,g9.3))

 end subroutine write_output_and_average
!====================================================================
 subroutine apply_bc
 use atoms_force_constants
 use lattice
 use md_params
 implicit none
 integer i
 real(8) dc(3),dr(3)

 if( bc.eq.1) then
     do i=1,natom_super_cell
! first get cartesian coordinates
        dc(1) = disp(1,i) + atom_sc(i)%equilibrium_pos%x
        dc(2) = disp(2,i) + atom_sc(i)%equilibrium_pos%y
        dc(3) = disp(3,i) + atom_sc(i)%equilibrium_pos%z
! then find the direct coordinates  take the distance using pbc, then
! transfrom back to cartesian coordinates
        call cart_to_direct(dc,gs1,gs2,gs3,dr)
! bring in [0,rs]
        dr = dr - floor(dr)
! bring between -box/2,box/2
!       dr = dr - box*nint(dr/box)
        call direct_to_cart(dr,rs1,rs2,rs3,dc)
        disp(1,i) = dc(1) - atom_sc(i)%equilibrium_pos%x
        disp(2,i) = dc(2) - atom_sc(i)%equilibrium_pos%y
        disp(3,i) = dc(3) - atom_sc(i)%equilibrium_pos%z
!       write(ulog,6)i,disp(:,i) !,forc(:,i)
     enddo
 endif

 end subroutine apply_bc
!============================================================
 subroutine read_fcs(iunit,rank,mx,iat,ixyz,fc)
 use atoms_force_constants
 use md_params
 use params
 use svd_stuff
 use ios , only : ulog
implicit none
 integer, intent(in) ::  rank,iunit
 integer, intent(inout) ::  mx
! integer, allocatable :: iat(:,:),ixyz(:,:)
! real(8), allocatable ::  fc(:)
 integer :: iat(rank,mx),ixyz(rank,mx)
 real(8) ::  fc(mx)
 integer t,i,junk,j
 character line*99

 write(*,*)'Entering read_fcs with rank,mx ',rank,mx
 t = 0
! mx = nterms(rank)
! allocate(iat(rank,mx),ixyz(rank,mx),fc(mx))
 if (ngroups(rank) .eq. 0) then
    allocate(fcs(1))
    nterms(rank) = 0
    fcs(:)=0
    include_fc(rank) = 0       ! if nothing, then exclude from formula
 else
    read(iunit,'(a)') line
    do i=1,mx               ! i is the index, t is just the number of terms in a given group
       read(iunit,*,end=91)t,junk, (iat(j,i),ixyz(j,i),j=1,rank),fc(i)
    enddo
 endif

91 mx=i-1 
write(ulog,*)'READ_FCS: rank=',rank,' reached end of file, new mx=',mx


 if ((i-1).gt.nterms(rank)) then
    write(*,*)'W: reading loop for rank ',rank,' exited after ',i-1,' terms!'
    write(*,*)'W: this is not consistent with nterms(rank)-',nterms(rank),' in lat_fc.dat file'
!   stop
 endif

 end subroutine read_fcs
!============================================================
 subroutine read_allfcs(potflag)
! which fcs are included, how many groups and terms for each rank, neighbor atoms to central primitive cell
 use ios, only: ufc,ulog 
 use svd_stuff
 use params, only: include_fc
 integer, intent(in) :: potflag
 integer rank,mx,iunit
 logical ex
 character xt*1,fn*20
 if (potflag.eq.0) then 

! get actual force constants from 4 files fc1,fc2,fc3,fc4
   do rank=1,4
      write(xt,'(i1)')rank
      fn = 'fc'//xt//'.dat'
      iunit = ufc+rank
      if (include_fc(rank) .ne. 0 ) then
         inquire(file=fn,exist=ex)
         if (ex) then
            open(iunit ,file=fn ,status='old')
            mx = nterms(rank)
            if (rank.eq.1 .and. mx.ne.0) then
               allocate(iatomterm_1(rank,mx),ixyzterm_1(rank,mx),fcs_1(mx))
               call read_fcs(iunit,rank,mx,iatomterm_1,ixyzterm_1,fcs_1)
               nterms(1)=mx
            elseif (rank.eq.2 .and. mx.ne.0) then
               allocate(iatomterm_2(rank,mx),ixyzterm_2(rank,mx),fcs_2(mx))
               call read_fcs(iunit,rank,mx,iatomterm_2,ixyzterm_2,fcs_2)
               nterms(2)=mx
            elseif (rank.eq.3 .and. mx.ne.0) then
               allocate(iatomterm_3(rank,mx),ixyzterm_3(rank,mx),fcs_3(mx))
               call read_fcs(iunit,rank,mx,iatomterm_3,ixyzterm_3,fcs_3)
               nterms(3)=mx
            elseif (rank.eq.4 .and. mx.ne.0) then
               allocate(iatomterm_4(rank,mx),ixyzterm_4(rank,mx),fcs_4(mx))
               call read_fcs(iunit,rank,mx,iatomterm_4,ixyzterm_4,fcs_4)
               nterms(4)=mx
            endif
         else
            write(*,*)'For this rank there is no FC data file present ',rank
            write(*,*)'check your inputs and run the program again'
            write(ulog,*)'For this rank there is no FC data file present ',rank
            write(ulog,*)'check your inputs and run the program again'
            stop
         endif
      endif
   enddo
 
 else
   write(*,*)'READ_ALLFCS: potflag is ',potflag,' and FCs are not read'
   return
 endif

 end subroutine read_allfcs
!============================================================
 subroutine read_supercell_coordinates
! the following reads the coordinates of the atoms in the supercell in VASP POSCAR format
 use md_params
 use lattice
 use geometry
 use atoms_force_constants
 use ios
 implicit none
 character line*90,coord_format*1
 integer i
 real(8) om,scale,a,b,c,pos(3)
! type(vector) pos

 open (uposcar,file='coordinates.sc',status='old')
 write(ulog,*)'OPENING the supercell coordinates file in POSCAR format ... reading the data'
 read(uposcar,'(a)') line
 write(ulog,'(a)') line

! read a scaling factor, and the translation vectors of the supercell
 read(uposcar,*)lattice_parameter
 read(uposcar,*)rs1
 read(uposcar,*)rs2
 read(uposcar,*)rs3
 if (lattice_parameter .lt.0) then  ! it is the -volume of the cell
    omega = -lattice_parameter  ! need to find the volume(r1,r2,r3) for scaling
    call calculate_volume(rs1,rs2,rs3,om)
    scale = (omega/om)**(1./3.)
 else
    scale = lattice_parameter
 endif
 rs1=scale*rs1; rs2=scale*rs2; rs3=scale*rs3
 call write_out(ulog,'rs1 ',rs1)
 call write_out(ulog,'rs2 ',rs2)
 call write_out(ulog,'rs3 ',rs3)
 box(1)=length(rs1); box(2)=length(rs2); box(3)=length(rs3)
 call calculate_volume(rs1,rs2,rs3,omega)
 write(ulog,*)' Volume of supercell is   =',omega
 write(ulog,*)' box length in 1st direction=',box(1)
 write(ulog,*)' box length in 2nd direction=',box(2)
 write(ulog,*)' box length in 3rd direction=',box(3)

! number of atoms of each type
 read(uposcar,*)(natom(i),i=1,natom_type)
! number of atoms participating in MD,
 natom_super_cell = sum(natom)
 write(ulog,*)" number of different elements and number of atoms of each type"
 write(ulog,*) natom_type,(natom(i),i=1,natom_type)
 write(ulog,*)" total number of atoms in the supercell= ", natom_super_cell
 write(ulog,*) "translation vectors of the supercell are"
 write(ulog,4) rs1
 write(ulog,4) rs2
 write(ulog,4) rs3
 allocate(atom_sc(natom_super_cell))


 read(uposcar,*) coord_format
! read equilibrium positions,types and masses from POSCAR file, write in logfile
 if (coord_format.eq.'D' .or. coord_format.eq.'d' .or.  &
 &   coord_format.eq.'R' .or. coord_format.eq.'r' ) then
! use below if direct coordinates (in units of r1,r2,r3)
     do i=1,natom_super_cell
        read(uposcar,*) a,b,c,atom_sc(i)%cell%tau,atom_sc(i)%mass
        atom_sc(i)%equilibrium_pos = a*rs1+b*rs2+c*rs3
     enddo
 elseif (coord_format.eq.'C' .or. coord_format.eq.'c' ) then
! use below if cartesian coordinates
     do i=1,natom_super_cell
        read(uposcar,*) pos,atom_sc(i)%cell%tau,atom_sc(i)%mass
        atom_sc(i)%equilibrium_pos = pos
!       atom_sc(j)%equilibrium_pos = scale*pos
!       do l=1,natoms0
!          if ( l .eq. atom_sc(j)%cell%tau ) then
!             atom_sc(j)%mass = atom0(l)%mass
!             exit
!          endif
!       enddo
     enddo
 endif

! now in atom_sc%equilibrium_pos we have the cartesian coordinates.
 close(uposcar)

 write(ulog,*)' supercell coordinates read successfully and closed'

 call make_reciprocal_lattice(rs1,rs2,rs3,gs1,gs2,gs3)

4 format(9(2x,f10.4))

 end subroutine read_supercell_coordinates
!====================================================================
 subroutine check_super_cell_coordinates
! check consistency of the atoms in the supercell and assign their n vector
 use md_params
 use params
 use lattice
 use geometry
 use atoms_force_constants
 use ios
 implicit none
 integer i,j,k,ipr,isc, counter,ier,nn(3)
 real(8) aa(3)
 type(vector) shift0,vec
 logical matched

! first describe r_i in terms of r0_i, and make sure the lin comb is integer
 write(ulog,*)' ----------------------------------------------------------'
 write(ulog,*)' CHECKING the commensurability between primcell and supercell'
 call check_int(rs1,aa,ier,g01,g02,g03)
 write(ulog,4)'checking rs1 in units of r0s,ier=',aa,ier
 if (ier .eq. 1) stop
 call check_int(rs2,aa,ier,g01,g02,g03)
 write(ulog,4)'checking rs2 in units of r0s,ier=',aa,ier
 if (ier .eq. 1) stop
 call check_int(rs3,aa,ier,g01,g02,g03)
 write(ulog,4)'checking rs3 in units of r0s,ier=',aa,ier
 if (ier .eq. 1) stop
 write(ulog,*)' COMMENSURABILITY CHECKED ---------------------------------'
4 format(a,3(2x,f8.3),2x,i1)

 write(ulog,*)'FINDING possible translation vectors, trying them on other atoms'
 checkloop: do i=1,natom_prim_cell  !,natom_super_cell
   shift0 = atom0(i)%equilibrium_pos   &
&       - atom_sc(1)%equilibrium_pos
   counter = 0
   PRIM:do j = 1,natom_prim_cell
      matched = .false.
      SC: do k = 1,natom_super_cell
         if(length(atom0(j)%equilibrium_pos - shift0 -    &
&           atom_sc(k)%equilibrium_pos) .lt. tolerance ) then
! found a matching vector... now identify the two atoms
            matched = .true.
            ipr=j   ! j and k have the same tau
            isc=k
!                       atom_sc(k)%cell%tau = atom0(j)%cell%tau
            exit SC
         endif
      enddo SC
      if (.not. matched) then  ! this tau was no good, try another tau
         exit PRIM
      else
         counter = counter + 1
         write(ulog,*)' atoms ',j,' in prim, and ',k,' in sc were matched,i=',i
         cycle PRIM
      endif
   enddo PRIM
   if (counter .eq. natom_prim_cell) then ! all the prim atoms were matched
       write(ulog,*)' tau matched all atoms of the primitive cell '
       call check_int(shift0,aa,ier,g01,g02,g03)
       aa = aa-floor(aa)
   endif
 enddo checkloop

 shift0 = aa(1)*r01 + aa(2)*r02 + aa(3)*r03
 write(ulog,*)'THE subtraction vector shift0 in red units ', shift0
! shift all the atoms in the primitive cell by -shift0 so that they fall
! on the atoms in the supercell
 write(ulog,*)'The shifted (cartesian) positions in the prim cell are now:'
 do i=1,natom_prim_cell
    atom0(i)%equilibrium_pos = atom0(i)%equilibrium_pos-shift0
    write(ulog,*)'i,shifted position ', i,atom0(i)%equilibrium_pos
 enddo

 write(ulog,*)'SC: i,    equilibrium position     ,  mass  , tau ,     n'
! now take every atom in the supercell and see if it can be mapped to any atom0
 ILOOP: do i=1,natom_super_cell
    matched = .false.
    counter = 0
    JLOOP: do j=1,natom_prim_cell
       vec = atom_sc(i)%equilibrium_pos - atom0(j)%equilibrium_pos
       call check_int(vec,aa,ier,g01,g02,g03)
       if( ier .eq. 0)  then         ! found a matching vector
           counter = counter + 1
           matched = .true.
! then assign its parameters such as type, mass, address...
           atom_sc(i)%cell%tau = atom0(j)%tau  ! atom_shl(j)%cell%tau 
           nn = nint(aa)
           atom_sc(i)%cell%n = nn
       endif
    enddo JLOOP
!   if (.not. matched) then
    if (counter.eq. 0) then
       write(ulog,*)' ATOM which can not be mapped to the prim cell ',i
       write(ulog,3)atom_sc(i)%equilibrium_pos
       stop
    elseif ( counter .ge. 2) then
       write(ulog,*)' CHECK: more than one atom matched! check your coordinates'
       stop
    endif
    write(ulog,5)i,atom_sc(i)%equilibrium_pos,atom_sc(i)%mass,  &
&                atom_sc(i)%cell%tau,atom_sc(i)%cell%n
 enddo ILOOP
 write(ulog,*)'*******************************************************'
 write(ulog,*)'SINCE THE PROGRAM WAS NOT STOPPED, MAPPING-CHECK PASSED!'
 write(ulog,*)'*******************************************************'
 if (.not.matched) then
    write(ulog,*)' SUPER_CELL_CHECK failed'
    write(*   ,*)' SUPER_CELL_CHECK failed'
    stop
 endif

3 format(9(1x,g12.6))
5 format(i6,4(1x,g12.6),4(2x,i3))

 end subroutine check_super_cell_coordinates

!============================================================
 subroutine move_atoms_gaussian(myseed,n,ampl,dsp,vl)
 use atoms_force_constants
 use md_params
 use constants, only : k_b,ee 
 use ios, only : ulog
 implicit none
 integer, intent(in) :: n,myseed
 real(8), intent(in) :: ampl
 real(8), intent(out) :: dsp(3,n),vl(3,n)
 integer i
 real(8) ampv2,avgmass

 do i=1,3
    call gaussian(myseed    ,n,dsp(i,:))
    call gaussian(myseed+800,n, vl(i,:))
 enddo
 dsp = dsp*ampl  ! in Ang 
 vl=vl*sqrt(k_b*temperature_k/ee) ! should be typical velocity amplitude but not centered at 0


! now that piavg=0, we set average vi^2 to 3kT/m_i 
 ampv2=0
 do i=1,n
    ampv2=ampv2+ atom_sc(i)%mass*(vl(:,i).dot.vl(:,i))/n  ! this will be rescaled to 3kT
 enddo
 write(ulog,*)'Target m*Velocity^2 amplitude vs 3kT(eV)=',ampv2,3*k_b*temperature_k/ee

 do i=1,n
    vl(:,i)=vl(:,i)*sqrt(3*atom_sc(i)%mass*k_b*temperature_k/ee/ampv2) 
 enddo

! if velocity unit = sqrt(eV/uma)=Ang/t_unit, then t_unit=ang/sqrt(eV/uma)=10.21641 fs
 write(ulog,*)' INITIAL DISPACEMENTS-VELOCITIES ARE:'
 do i=1,n
   write(ulog,3)i, dsp(:,i),vl(:,i)
 enddo

     open(123,file='dispvel0.dat')
     do i=1,n
        write(123,3)i,dsp(:,i),vl(:,i)
     enddo
     close(123)

3 format(i6,9(1x,g12.5))
 end subroutine move_atoms_gaussian
!============================================================
! function transl(n) result(y)
! use lattice
! integer, intent(in) :: n(3)
! real(8) y(3)
!
! y = n(1)*r01+n(2)*r02 + n(3)*r03
!
! end function transl
!============================================================
 subroutine calculate_potential_lj
! for given displacements and force constants, this subroutine
! constructs the potential and the force function iand the current based
! on the equilibrium positions of atoms and superlattice specifications
 use OMP_LIB
 use md_params
 use geometry
 use lattice
 use atoms_force_constants
 use ios, only : ulog
 implicit none
 real(8) random(3,natom_super_cell),fij(3)
 real(8) const,dr(3),dc(3),rij,r_ij(3),fdotv,arg,vij,dv,d2v
 integer al,i,j,thread_id,n_atom
! integer lx,ly,lz

 n_atom = size(disp,2)
! write(*,*)' size of disp=number of atoms=',n_atom
 potential = -potential0 ; forc = 0; current = 0; cur = 0

! !$OMP PARALLEL DO !I is private by default


! !$OMP PARALLEL PRIVATE(thread_id)
! 
!     thread_id = OMP_GET_THREAD_NUM()
! 
!     DO i=0,OMP_GET_MAX_THREADS()
!         IF (i == thread_id) THEN
!             PRINT *, "Hello from process: ", thread_id
!         END IF
!     END DO
! !$OMP END PARALLEL


! !  !#, OMP_GET_THREAD_NUM()

!$OMP PARALLEL PRIVATE(vij) SHARED(potential,forc,cur,current)
 ILOOP: do i=1,natom_super_cell
    do j=1,natom_super_cell  ! stupid, really need a neighbor_list here...
       if(i.eq.j ) cycle ! no self-interaction
       r_ij = v2a(atom_sc(i)%equilibrium_pos-atom_sc(j)%equilibrium_pos)
       dc   = r_ij + disp(:,i)-disp(:,j) 

! bring dr between -1/2 and 1/2 or r_ij between -box/2 and box/2
       call cart_to_direct(dc,gs1,gs2,gs3,dr)
       dr   = dr - anint(dr)
       call direct_to_cart(dr,rs1,rs2,rs3,dc)
       rij  = length(dc)   ! sqrt( dx*dx + dy*dy + dz*dz )

       if (rij .lt. 0.001*rcut2) then
         write(ulog,*)'LJ: atoms too close:i,j,rij ',i,j,rij,rcut2
         stop
       elseif(rij.gt.rcut2) then
         cycle
       endif

! vij is the pair interaction potential, and dv = v'(rij)/rij ; force(i) is sum_j -dv*r_ij
       call pot_force(rij,vij,dv,d2v)

       potential = potential + vij
       fij= -dv*dc   ! this is -d v_ij/d r_i = d v_ij/dr_j = -d v_ij/d rij (ri-rj)/rij
!       fij(:) = fij(:)+forc0(:,i,j)
       forc(:,i) = forc(:,i) + fij
       fdotv = fij(:) .dot. velo(:,i) /2
!       fdotv = fij(:) .dot. (velo(:,i) + velo(:,j)) /4   ! symmetrized expression
       cur (:,i) = cur(:,i) + fdotv * r_ij

    enddo
    current(:) = current(:) + cur(:,i)
!    write(ucurrent,7)i,cur(:,i)
 enddo ILOOP
!$OMP END PARALLEL 

 if (relax) then ! apply heatbath contact
    const = sqrt(6*damp_rate*temperature/dt)
    call random_number(random); random = (2*random-1)*const
!   write(ulog,*)'random numbers',random
    do i=1,natom_super_cell
    do al=1,3
       forc(al,i)=(forc(al,i)+random(al,i)*sqrt(atom_sc(i)%mass)-atom_sc(i)%mass * &
&                      damp_rate*velo(al,i) )/(1+dt*damp_rate/2)
    enddo
    enddo
 endif

 potential = potential/natom_super_cell/2
 current   = current  /omega
5 format(a,i6,i6,1(1x,f6.3),9(1x,g9.3))
7 format((1x,i4),3(1x,g10.4))

 end subroutine calculate_potential_lj
!============================================================
 subroutine calculate_potential_lj0
!! this subroutine constructs the potential and the force function and the current based
!! on the equilibrium positions of atoms and superlattice specifications
 use md_params
 use geometry
 use lattice
 use atoms_force_constants
 use ios, only : ulog
 implicit none
 real(8) sumfi(3),r_ij(3),red(3),arg,vij,dv,d2v,rij,fdotv
 integer i,j,n_atom

 n_atom = size(disp,2)
! set the potential parameters: energy and lengthscale and cutoff radius
 sigma=2.43d0; eps =5.1d0; rcut2=3.75
! sigma=4.; eps =1; rcut=1.75     ! for potflag=2
! sigma=1.414; eps =0; rcut=1.75  ! for potflag=3
 write(ulog,5)'LJ0: natom,potential parameters: sigma, eps, rcut=',n_atom,sigma,eps,rcut2

 potential0 = 0; velo=0; phi=0 !; current = 0; cur = 0
 do i=1,natom_super_cell
    sumfi = 0d0
    do j=1,natom_super_cell
       if(i.eq.j ) cycle ! no self-interaction
       r_ij = v2a(atom_sc(i)%equilibrium_pos - atom_sc(j)%equilibrium_pos)

! bring red between -1/2 and 1/2 or r_ij between -box/2 and box/2
       call cart_to_direct(r_ij,gs1,gs2,gs3,red)
       red   = (red - anint(red))
       call direct_to_cart(red,rs1,rs2,rs3,r_ij)

       rij = length(r_ij) 
       if (rij .lt. 0.001*rcut2) then
         write(ulog,*)'LJ0: atoms too close:i,j,rij ',i,j,rij
         stop
       elseif(rij.gt.rcut2) then
         cycle
       endif

! vij is the pair interaction potential, and dv = v'(rij)/rij ; force(i) is sum_j -dv*r_ij
       call pot_force(rij,vij,dv,d2v)

! phi_ii=sum_j (d2v-dv/rij) rij_a rij_b /rij^2 + d_ab dv/rij
       if (i.eq.1) then
          phi=phi+(d2v-dv/rij)/3+dv/rij
       endif
       potential0 = potential0 + vij
!       forc0(:,i,j) = -dv*r_ij ! we keep forc0 in case initial configuration is not equilibrium
       sumfi = sumfi -dv*r_ij ! total force on i 
!       fdotv = -dv* (r_ij .dot. velo(:,i))/2  ! velocity has  nnot been intialized yet
!       cur (:,i) = cur(:,i) + fdotv*r_ij

     enddo

!     current(:) = current(:) + cur(:,i)
     write(ulog,5)'LJ0: i,force(i)=',i,sumfi 

 enddo

 write(ulog,*)'reference PE=',potential0/ natom_super_cell

5 format(a,i6,9(1x,g9.3))

 end subroutine calculate_potential_lj0
!===========================================================
 subroutine pot_force(rij,vij,dv,d2v)
! vij is the pair interaction potential, and dv = v'(rij)/rij ; force(i) is sum_j -dv*r_ij
 use md_params
 use ios, only : ulog
 implicit none
 real(8), intent(in) :: rij
 real(8), intent(out) :: vij,dv,d2v
 real(8) arg

         if (potflag.eq.1) then
! this is for the LJ potential: eps*( (sigma/rij)**12 - (sigma/rij)**6 )
            arg=(sigma/rij)**6
            vij= eps*arg*(arg-1d0)
            dv = eps*arg*  (6-12d0*arg)/(rij*rij)
            d2v= eps*arg*6*(7-  26*arg)/(rij*rij)
         elseif(potflag.eq.2) then
! this is for the morse potential: eps*(2 exp(-2r/sigma)-exp(-r/sigma) )
            arg=exp(-rij/sigma)
            vij= eps*arg*(2*arg-1)
            dv = eps/sigma/rij*arg*(1-4*arg)
            d2v= eps/(sigma*sigma)*arg*(8*arg-1)
         elseif(potflag.eq.3) then
! this is for the quartic potential: y**2/2 + eps*(y**3/6+y**4/24) with y=rij-sigma
            arg = rij-sigma   ! eps=strength of anharmncty
            vij = arg*arg/2d0 *(1+eps*arg/3d0*(1+arg/4))
            dv  = (arg/rij * (1 + eps*arg/2d0*(1 + arg/3)))
            d2v = 1+eps*arg*(1+arg/2)
         else
            write(ulog,*)' potflag must either be 1,2 or 3, not',potflag
            stop
         endif

 end subroutine pot_force
!============================================================
      subroutine findatom_sc2(n3,tau,iatom)
! find atom in supercell
! arguments:
!     n(i) (input,output), linear combination of basis vectors of the primitive
!          lattice that takes us to the unit cell containing the ith atom
!     tau (input), identity of equivalent atom in unit cell at origin
!     iatom (output), location of atom in supercell. Returns zero if not found
      use geometry
      use lattice
      use md_params
      use params
      use atoms_force_constants
      use ios, only : ulog
      implicit none
      integer n3(3),tau,iatom,j
      type(vector) v,w

!      jloop: do j=1,natoms0
!         if ( tau .eq. atom0(j)%tau ) then
!            v=atom0(j)%equilibrium_pos
!            exit jloop
!         endif
!      enddo jloop

      v = atom0(tau)%equilibrium_pos
      v = v + n3(1)*r01 + n3(2)*r02 + n3(3)*r03
! v is in cartesian coordinates, now take it into the supercell (periodic BCs)
      call bring_to_cell_c(v,rs1,rs2,rs3,gs1,gs2,gs3,w)

! find its atom number
      iatom = 0
      do j=1,natom_super_cell
         if (atom_sc(j)%cell%tau  .ne. tau   ) cycle
         if (length(atom_sc(j)%equilibrium_pos - w) .gt. tolerance) cycle
!        if (atom_sc(j)%cell%n(1) .ne. n3(1) ) cycle
!        if (atom_sc(j)%cell%n(2) .ne. n3(2) ) cycle
!        if (atom_sc(j)%cell%n(3) .ne. n3(3) ) cycle
            iatom = j
            return
      enddo
      if (iatom.eq.0) then
         write(ulog,*)'FINDATOM_SC2:no atom found in super cell '
         write(ulog,*)'tau,n,position ',tau,n3,w
         stop
      endif

5 format(a,3(2x,f9.4))
      end subroutine findatom_sc2
!=====================================================================
      subroutine findatom_sc(n3,tau,iatom)
! find atom in data base
! arguments:
!     n(i) (input,output), linear combination of basis vectors of the primitive
!          lattice that takes us to the unit cell containing the ith atom
!     tau (input), identity of equivalent atom in unit cell at origin
!     iatom (output), location of atom in supercell. Returns zero if not found
      use geometry
      use lattice
      use atoms_force_constants
      use constants, only : pi
      use ios, only : ulog
      implicit none
      integer n3(3),tau,iatom,j,m(3)
      real(8) a(3),b(3),zero(3)
      type (vector) w
      zero = 0d0
      iatom = 0
      do j=1,natom_super_cell
         if (atom_sc(j)%cell%tau .eq. tau ) then
             m = atom_sc(j)%cell%n - n3
             w = m(1)*r01+m(2)*r02+m(3)*r03
             a(1)=gs1.dot.w /(2*pi)
             a(2)=gs2.dot.w /(2*pi)
             a(3)=gs3.dot.w /(2*pi)
             m = floor(a+0.00001)
             b = a - m
             if (b .myeq. zero ) then
                iatom = j
                return
             endif
         endif
      enddo
      if (iatom.eq.0) then
         write(ulog,*)'FINDATOM_SC:no atom found in super cell '
         write(ulog,*)'tau,n,position ',tau,n3
         stop
      endif
      end subroutine findatom_sc
!=====================================================================
 subroutine make_reciprocal_lattice(r1,r2,r3,g1,g2,g3)
! be careful it does not have the factor of 2*pi in it!
 use geometry
 use constants
 implicit none
!real(8), intent(in) :: r1(3),r2(3),r3(3)
!real(8), intent(out):: g1(3),g2(3),g3(3)
 type(vector), intent(in) :: r1,r2,r3
 type(vector), intent(out):: g1,g2,g3
 real(8) om

 om = abs(r1 .dot. (r2 .cross. r3))
 if (om.lt.1d-8) then
    write(*,*)'MAKE_RECIPROCAL_LATTICE_A: volume is zero; check your translation vectors'
    stop
 endif
 g1=(r2 .cross. r3)/(om /(2*pi))
 g2=(r3 .cross. r1)/(om /(2*pi))
 g3=(r1 .cross. r2)/(om /(2*pi))

 end subroutine make_reciprocal_lattice
!===============================================================
 subroutine gaussian(myseed,n,x)
!! generages n random numbers x of mean 0 and variance 1
 use md_params
 implicit none
 integer, intent(in) :: n,myseed
 real(8), intent(out) :: x(n)
 real(8) ksi(n),psi(n) 
 integer i,n2
 integer, allocatable:: iseed(:)

  call random_seed(size = n2)
  write(*,*)'n,n2=',n,n2
  allocate(iseed(n2))
  iseed=myseed+(/ (37*i,i=0,2*n2-2,2 ) /)
  call random_seed(put=iseed)
  write (*,*)'Generated seeds are:', iseed

! call random_init(.true.,.true.)
! call random_seed ( )   ! initialize based on time and date
 call random_number(ksi)
! call random_init(.false.,.true.)
 call random_number(psi)

 open(61,file='random.dat')
 write(61,*)'# random numbers ',n
 do i=1,n
    write(61,*) psi(i) , ksi(i)  ! uniform distribution in [0,1]
 enddo
 close(61)

 psi=psi*6.28318531
 ksi=sqrt(-log(1-ksi))
 write(*,*) ' generating ',n,' random number of gaussian distribution'

 open(51,file='gauss_0_1.dat')
 do i=1,n
    x(i)= cos(psi(i))*ksi(i)
    write(51,4) ' ',x(i)  ! x has a gaussian distribution
 enddo
 close(51)

 deallocate(iseed)

4 format(a,9(1x,f9.4))
 end subroutine gaussian
!===========================================================
 function is_integer(i)
!! checks if the variable is integer
 implicit none 
 real(8) i
 logical is_integer

 if (abs(i-nint(i)).lt.1d-5 ) then 
    is_integer=.true.
 else 
    is_integer=.false.
 endif
 end function is_integer
!===========================================================
 subroutine reset_momentum(n,vl)
! sets total momentum to zero. That reduces the kinetic energy as well.
 use atoms_force_constants
 use ios, only : ulog
 implicit none
 integer, intent(in) :: n
 real(8), intent(inout) :: vl(3,n)
 integer i
 real(8) pavg(3)

 pavg=0
 do i=1,n
    pavg=pavg+atom_sc(i)%mass*vl(:,i)/n
 enddo
 write(ulog,*)'Initial pavg=',pavg

 do i=1,n
    vl(:,i)=vl(:,i)-pavg/atom_sc(i)%mass
 enddo
 write(ulog,*)'sum of velocity components after subtraction=',sum(vl) 

 end subroutine reset_momentum
!===========================================================
 subroutine rescale_velocity(n,vl,tempk)
! rescale velocities such that mv^2 = 3 k_b Tempk
 use constants, only : k_b,ee
 use atoms_force_constants
 implicit none
 integer, intent(in) :: n
 real(8), intent(in)    :: tempk
 real(8), intent(inout) :: vl(3,n)
 real(8) ke_avg
 integer i

 ke_avg=0
 do i=1,n
    ke_avg=ke_avg + atom_sc(i)%mass*dot_product(vl(:,i),vl(:,i))/n
 enddo

 vl=vl*sqrt(3*k_b*tempk/ee/ke_avg)

 end subroutine rescale_velocity
