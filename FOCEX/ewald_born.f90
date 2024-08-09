!=====================================

 module ewald
 use constants, only : r15,pi,ci
 integer nr_ewald,ng_ewald
 real(r15), allocatable :: r_ewald(:,:),g_ewald(:,:)
 real(r15) rcutoff_ewa,gcutoff_ewa, eta

 contains

!============================================================
 subroutine make_grid_shell_ewald(x1,x2,x3,rcut,gcut,epsilo)
!! generates a grid of vectors "grid" linear combinations of xi within a
!! cutoff length. Sorted result is stored in r_ewald(3,nr_ewald)
!! also generates the same grid for the reciprocal space g_ewald(3,ng_ewald)
!! metric is epsilon in reciprocal space and epsinverse in real space
 use geometry
 use ios
 use lattice
 implicit none
 real(r15), intent(in ):: rcut,gcut,epsilo(3,3)
 type(vector), intent(in ):: x1,x2,x3
! integer, intent(inout):: ngrid
! real(r15), intent(inout):: grid(:,:)
! real(r15), allocatable, intent(inout):: grid(:,:)
 integer i1,i2,i3,m,cnt,maxx
 real(r15) v(3),rr,aux2(3,3),epsinv(3,3)
 real(r15), allocatable :: aux(:,:),lengths(:)
 integer, allocatable :: msort(:)
 type(vector) y1,y2,y3

! call apply_metric(x1,x2,x3,epsilo,t1,t2,t3)
 aux2=epsilo
 call inverse_real(aux2,epsinv,3)

 call get_upper_bounds(x1,x2,x3,rcut,maxx,m)
! call calculate_volume(x1,x2,x3,om0)
! max=nint(15d0/3d0*rcut**3/om0)+10
! m=nint((max/4.2)**0.333)+4

 allocate(aux(3,maxx),lengths(maxx),msort(maxx))
 lengths =1d20; aux=1d20
 write(ulog,4)'MAKE_GRID_SHELL, maxx(Rgrid), for rcut=',maxx,rcut
 write(*,4)'MAKE_GRID_SHELL, maxx(Rgrid), m for rcut=',maxx,float(m),rcut
 cnt=1
 do i1=-m,m
 do i2=-m,m
 do i3=-m,m

! generate vectors in a grid
    v= v2a(i1*x1 + i2*x2 + i3*x3)
    rr=sqrt(dot_product(v,matmul(epsinv,v)))
    if(rr.gt.rcut) cycle

    aux(:,cnt)=v
    lengths(cnt)=rr
    cnt=cnt+1
    if (cnt.gt.maxx) then
       write(ulog,*)'RLOOP: maxx size exceeded, need to increase variable maxx from ',maxx
       stop
    endif

 enddo
 enddo
 enddo

 nr_ewald=cnt-1
 write(ulog,5) 'MAKE_GRID_SHELL: within rcut=',rcut,' r_ewald generated ',nr_ewald,' vectors'

 if (allocated(r_ewald)) deallocate(r_ewald)  ! overwrite if previously called
 allocate(r_ewald(3,nr_ewald))

 call sort(nr_ewald,lengths,msort,maxx)

 open(ujunk,file='ewald_vecs.dat')
 write(ujunk,*)'MAKE RGRID:----first 200 vectors--------------',nr_ewald
 do i1=1,nr_ewald
    r_ewald(:,i1)=aux(:,msort(i1))
    if (i1.lt.200) write(ujunk,3)r_ewald(:,i1)
 enddo

 deallocate (aux, lengths, msort)

! now the G-vectors
! call calculate_volume(y1,y2,y3,om0)
! max=2*nint(12.6/3d0*gcut**3/om0)+10
! m=nint(max**0.333)*2
 call make_reciprocal_lattice_2pi(x1,x2,x3,y1,y2,y3)
 call get_upper_bounds(y1,y2,y3,gcut,maxx,m)
 write(ulog,4)'MAKE_GRID_SHELL, maxx(g_ewald), for gcut=',maxx,gcut
 write(*,4)'MAKE_GRID_SHELL, maxx(g_ewald), m for gcut=',maxx,float(m),gcut
 allocate(aux(3,maxx),lengths(maxx),msort(maxx))
 lengths =1d20; aux=1d20

 cnt=1
 do i1=-m,m
 do i2=-m,m
 do i3=-m,m

! generate vectors in a grid
    v= v2a(i1*y1 + i2*y2 + i3*y3)
    rr=sqrt(dot_product(v,matmul(epsilo,v)))
    if(rr.gt.gcut) cycle

    aux(:,cnt)=v
    lengths(cnt)=rr
    cnt=cnt+1
    if (cnt.gt.maxx .and. maxx.ne.1) then
       write(ulog,*)'GLOOP: maxx size exceeded, need to increase variable maxx from ',maxx
       stop
    endif

 enddo
 enddo
 enddo

 ng_ewald=cnt-1
 write(ulog,5) 'MAKE_GRID_SHELL: within gcut=',gcut,' g_ewald generated ',ng_ewald,' vectors'

 if (allocated(g_ewald)) deallocate(g_ewald)
 allocate(g_ewald(3,ng_ewald))

 call sort(ng_ewald,lengths,msort,maxx)

 write(ujunk,*)'MAKE GGRID:------------------------',ng_ewald
 do i1=1,ng_ewald
    g_ewald(:,i1)=aux(:,msort(i1))
    if (i1.lt.200) write(ujunk,3)g_ewald(:,i1)
 enddo

 deallocate (aux, lengths, msort)
 close(ujunk)

3 format(9(1x,f10.5))
4 format(a,i8,9(1x,g14.7))
5 format(a,g11.4,a,i8,9(1x,g14.7))

 end subroutine make_grid_shell_ewald
!============================================================
 subroutine dewapot(x,dpot)
!! calculates the -d/dx(sum_R 1/|R+x| -background) = force  with the corrected metric
!! Assumes translation vectors are r_ewald(:,igrid),g_ewald(:,igrid)
! use constants
 use lattice
 use params
 use born , only : epsil,epsinv
 implicit none
 real(r15), intent(out) :: dpot(3)
 real(r15), intent(in) :: x(3)
 integer igrid
 real(r15) termg,qpg(3),dd,my_erfc,geg,ep,vd(3)

 if (length(x).lt.1d-12) then
    write(*,*)'DEWAPOT: X is too small! ',x
    ep=1d-12
 else
    ep=0
 endif

    dpot=0
    do igrid=1,nr_ewald
       vd=eta*(x+r_ewald(:,igrid))
       dd=sqrt(dot_product(vd,matmul(epsinv,vd)))
       if(dd.lt.1d-12) dd=1d-12
       if (dd.gt.6) cycle
       dpot=dpot + eta*eta*(matmul(epsinv,vd) +  matmul(transpose(epsinv),vd))/2d0/dd/dd *  &
&           (my_erfc(dd)/dd+2/sqrt(pi)*exp(-dd*dd))/sqrt(det(epsil))
    enddo

! write(*,3)'dpotr=',dpot

    do igrid=2,ng_ewald
       qpg=g_ewald(:,igrid)
       geg=dot_product(qpg,matmul(epsil,qpg))
       if (geg.gt.100*eta*eta) exit
       termg=exp(-geg/4/eta/eta)/geg
       dpot=dpot + termg * qpg*sin(qpg .dot. x) *4*pi/volume_r  ! this is supercell volume
    enddo

! write(*,3)'dpotg=',dpot

    dpot=dpot + 4*pi/volume_r*matmul(epsinv,x)/3d0

! write(*,3)'dpot=',dpot
3 format(9(1x,g11.4))

 end subroutine dewapot
!===========================================================
 subroutine ewaldforce(n,dsp,fewald)
!! calculates the Coulomb force on atoms in supercell (in 1/Ang)
!! Assumes translation vectors are r_ewald(:,igrid),g_ewald(:,igrid)
 use constants, only : ee, eps0
 use lattice
 use params
 use atoms_force_constants
 use born
 use ios
 implicit none
 integer, intent(in) :: n  ! number of atoms in the supercell will be used
 real(r15), intent(in) :: dsp(3,n) ! displaced positions including equilibrium positions
 real(r15), intent(out):: fewald(3,n)
 integer tau1,tau2
 real(r15) dpot(3),r12(3),z1(3,3),z2(3,3),f1(3),f2(3) ,coef

 coef=ee/(4*pi*eps0)*1d10
 fewald=0
 do tau1=1,n  ! caculate the total force on tau1
    z1=atom_sc(tau1)%charge
 do tau2=1,n
    z2=atom_sc(tau2)%charge

    if (tau1.eq.tau2) cycle ! images of tau1 have zero force on tau1
    r12=dsp(:,tau1)-dsp(:,tau2) ! dsp includes equilibrium positions
! force on tau1 from tau2 and all its images
    call dewapot(r12,dpot)  ! has epsil-modified metric included
    f1=matmul(atom_sc(tau2)%charge,dpot)
    f2=matmul(atom_sc(tau1)%charge,f1)
    fewald(:,tau1)=fewald(:,tau1) + f2 * coef  ! to convert to eV/Ang

 enddo
 enddo

! now compare to the energy finite difference
! write(ulog,*)'******** TEST OF COULOMB FORCE ***********'
! allocate(dsp1(3,n))
! dsp1=dsp
  do tau1=1,n
!    call coulomb_energy(tau1,n,atom_sc,dsp,ecoul0) !assumes all other atoms are @eequilibrium
!    do al=1,3
!       dsp1(al,tau1)=dsp(al,tau1)+1d-5
!       call coulomb_energy(tau1,n,atom_sc,dsp1,ecoul1)
!       force1=(ecoul0-ecoul1)/1d-5
!       write(ulog,3)'al,tau,force,-de/dx=',al,tau1,frc(al,tau1),force1
        write(*,3)'force=',tau1,tau1,fewald(:,tau1)
!    enddo
  enddo
! deallocate(dsp1)

3 format(a,2i5,9(1x,g11.4))

 end subroutine ewaldforce
!===========================================================
 subroutine ewapot(x,pot)
! calculates the sum_R 1/|R+x| -background = pot (in 1/Ang no epsil)
! use constants , only: pi
 use lattice , only: volume_r
! use params, only : eta
 use geometry
 use born
 implicit none
 real(r15), intent(out) :: pot
 real(r15), intent(in) :: x(3)
 integer igrid
 real(r15) termg,qpg(3),dd,my_erfc,geg,ep,vd(3)

 if (length(x).lt.1d-12) then
    write(*,*)'EWAPOT: X is too small! ',x
    ep=1d-12
 else
    ep=0
 endif

    pot=0
    do igrid=1,nr_ewald
       vd=eta*(x+r_ewald(:,igrid))
       dd=sqrt(dot_product(vd,matmul(epsinv,vd)))+ep
       if (dd.gt.5) cycle  ! r_ewald is sorted but not x+r_ewald
       pot=pot + eta*my_erfc(dd)/dd/sqrt(det(epsil))
    enddo

    do igrid=2,ng_ewald
       qpg=g_ewald(:,igrid)
       geg=dot_product(qpg,matmul(epsil,qpg))
       if (geg.gt.100*eta*eta) exit
       termg=exp(-geg/(4*eta*eta))/geg
       pot=pot+ termg * cos(qpg .dot. x) *4*pi/volume_r
    enddo

    pot=pot - pi/volume_r/eta/eta   &
&   * (1+ 2/3d0*eta*eta* dot_product(x,matmul(epsinv,x)) )
! the second term above is the angular avge of G=0 extra-term


 end subroutine ewapot
!===========================================================
 subroutine coulomb_energy(i,n,atoms,dsp,e_coul)
! calculates the Coulomb interaction between atom i in primitive cell and the rest (in eV/Ang)
 use constants, only : ee,eps0
! use lattice
! use params
 use atoms_force_constants
! use born
 implicit none
 integer, intent(in) :: i,n
 type(atomic), intent(in) :: atoms(n)
 real(r15), intent(in) :: dsp(3,n)
 real(r15), intent(out) :: e_coul
 integer j
 real(r15) x(3),pot,zi(3,3),zj(3,3),emad,prod(3,3),trp

 call madelung2(emad)   ! sum'_R 1/R
! eps=trace(3,epsil)/3d0

    e_coul=0
    zi=atoms(i)%charge
    do j=1,n
       zj=atoms(j)%charge
       if(i.eq.j) then
          prod=matmul(zi,transpose(zi))
          trp=sum( (/ (prod(i,i),i=1,3)/) )   ! this is the trace function
          e_coul=e_coul+emad*trp/3d0
       else
!          x=atoms(i)%equilibrium_pos+dsp(:,i)-atoms(j)%equilibrium_pos-dsp(:,j)
          x=dsp(:,i)-dsp(:,j)
          prod=matmul(zi,transpose(zj))
          trp=sum( (/ (prod(i,i),i=1,3)/) )   ! this is the trace function
          call ewapot(x,pot)   ! sum_R 1/|R+x|
          e_coul=e_coul+pot*trp/3d0
       endif
    enddo

! e_coul=e_coul/eps ! for isotropic materials

 e_coul=e_coul*ee/(4*pi*eps0)*1d10  ! to convert to eV

 end subroutine coulomb_energy
!===========================================================
 subroutine madelung2(emadelung)
! Madelung energy from ewapot (in 1/Ang no epsil)
 use born
 use geometry
 implicit none
 real(r15), intent(out) :: emadelung
 real(r15) x(3),ene !det

 x=0; x(1)=1d-6
 call ewapot(x,ene)
 emadelung=ene-1/sqrt( dot_product(x,matmul(epsinv,x))*det(epsil) )

 end subroutine madelung2
!===========================================================
 subroutine test_ewald
 use ios , only : ulog
 implicit none
 integer al
 real(r15) x(3),dpot(3),r12(3),potp,potm,force1,pot0

 write(ulog,*)'**************** TEST OF COULOMB FORCE *********************'

 r12=(/0.1,0.2,0.5/)  !0.3*v2a(r01)
 write(ulog,3)' r12=',al,r12
 call dewapot(r12,dpot)

! now compare to the energy finite difference
 call ewapot(r12,pot0)
 do al=1,3
    x=r12
    x(al)=r12(al)+1d-4
    call ewapot(x,potp)
    x(al)=r12(al)-1d-4
    call ewapot(x,potm)
    force1=-(potp-potm)/2d-4
    write(ulog,3)'al,pot0,force,-de/dx=',al,pot0,dpot(al),force1
 enddo

3 format(a,i5,9(1x,g14.7))

 end subroutine test_ewald
!===========================================================
 subroutine subtract_coulomb_force(born_flag,ncfg,dsp,frc)
! dsp is the cartesian position not including the equilibrium positions
 use atoms_force_constants
! use constants, only : pi
 use fourier
 use lattice, only : rs1,rs2,rs3,g01,g02,g03,rws26, volume_r !, volume_r0
 use ios , only : ulog,write_out
 use params, only : verbose, tolerance
 use born , only : epsil
 implicit none
 integer, intent(in):: born_flag,ncfg
 real(r15), intent(inout) :: dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg)
! automatic arrays here (lost after the call, memory released)
 real(r15) frcr(3,natom_prim_cell,nrgrid),frcg(3,natom_prim_cell,nggrid)
 real(r15) disr(3,natom_prim_cell,nrgrid),disg(3,natom_prim_cell,nggrid)
 real(r15) fcoulg(3,natom_prim_cell,nggrid),fcoulr(3,natom_prim_cell,nrgrid)
 real(r15) ewald_force(3,natom_super_cell),auxr(nrgrid),auxg(nggrid)
! real(r15), allocatable :: frcr(:,:,:),frcg(:,:,:),disr(:,:,:),disg(:,:,:)
 real(r15) dyn_na(natom_prim_cell,natom_prim_cell,3,3,nggrid),ddn(3,3,3,nggrid)
 real(r15) phi_na(natom_prim_cell,natom_prim_cell,3,3,nrgrid)
 real(r15) rr(3),deteps3,asr(3,3), foldedr(3),pos(3,natom_super_cell) !rcutoff_ewa,eta,gcutoff_ewa,
 integer nsc,n3(3),tau,taup,al,be,icfg,i,j,l,jj,one

 one=1
 write(*,*)'SUBTRACT_COULOMB_FORCE: size of force is 3,',size(frc(1,:,1)),size(frc(1,1,:))

! Beginning of ewald part for deduction of ewald force in real space ------------------
  if(born_flag.le.0) then      ! skip if born_flag=<0 , otherwise, the NA term is always added
     return
  elseif(born_flag.eq.2) then  ! Ewald sum if=1, -FT(Dyn(K)*disp(k)) if=2

     write(ulog,*)'MAIN: going to call Ewaldforce '
! set cutoffs for Ewald sums ! this should come after the supercell is read!!!!
!     allocate(ewald_force(3,natom_super_cell))
     deteps3=det(epsil)**0.3333
     eta=sqrt(pi*deteps3)/(volume_r**0.3333) ! so that both R-sums and G-sums converge at the same rate
     rcutoff_ewa=10*sqrt(deteps3)/eta
     gcutoff_ewa=10*eta*2/sqrt(deteps3)
! generate real space and reciprocal space translation vectors of the supercell for Ewald sums
     call make_grid_shell_ewald(rs1,rs2,rs3,rcutoff_ewa,gcutoff_ewa,epsil)
! output is r_ewald(3,nr_ewald) and g_ewald(3,ng_ewald)

     do icfg=1,ncfg

        pos=v2a(atom_sc(:)%equilibrium_pos) + dsp(:,:,icfg)

! here ewaldforce requires the absolute cartesian coordinates of displaced atoms, so dsp=eq+displacement
!       call ewaldforce(natom_super_cell,dsp(:,:,icfg),ewald_force)
        call ewaldforce(natom_super_cell,pos          ,ewald_force)


        if ( verbose ) then
           write(ulog,*)' Ewald force for configuration #',icfg
           do j=1,natom_super_cell
              write(ulog,4)j, ewald_force(:,j)
           enddo
        endif

        frc(:,:,icfg)=frc(:,:,icfg) - ewald_force(:,:)

     enddo

!     deallocate(ewald_force)

  elseif(born_flag.eq.3) then  ! use non-analytical subtraction (in reciprocal space) F_SR(q)=F_DFT(q) - F_NA(q)

! fourier transform forces and displacements in the supercell, subtract the non-analytical term
! and fourier transform back before svd

!    allocate(map_rtau_sc(natom_prim_cell,nrgrid ))  ! needed for extract_fourier

! finds the mapping supercell atom<-> (tau,primcell translations defined in rgrid)
!    call find_map_rtau2sc(nrgrid,rgrid,map_rtau_sc)

!  erfc(4)=1.5d-8  so if eps=4/vol**.333 the real space terms become negligible
! exp(-18)=1.5d-8  so G larger than 5-6 shells will not contribute even if eps=4/vol^.3

!    allocate(frcr(3,natom_prim_cell,nrgrid),frcg(3,natom_prim_cell,nggrid))
!    allocate(disr(3,natom_prim_cell,nrgrid),disg(3,natom_prim_cell,nggrid))

     do tau=1,natom_prim_cell
     do taup=1,natom_prim_cell
        do j=1,nggrid
           call dyn_coulomb(tau,taup,ggrid(:,j),dyn_na(tau,taup,:,:,j),ddn(:,:,:,j))
        enddo
     enddo
     enddo

     config: do icfg=1,ncfg

! for each config convert force&disp to a 3 x N0 x Nrgrid format needed for FT
        do tau=1,natom_prim_cell
           do j=1,nrgrid
! find reduced coordinates of rgrid
              rr(1)=(rgrid(:,j) .dot. g01)/2/pi
              rr(2)=(rgrid(:,j) .dot. g02)/2/pi
              rr(3)=(rgrid(:,j) .dot. g03)/2/pi
              n3=nint(rr)
              if(length(rr-n3).gt.1d-4) then
                  write(*,5)' SUBTRACT_COULOMB:n3.ne.rr ',n3,rr
                  stop
              endif
              call findatom_sc(n3,tau,nsc)
              frcr(:,tau,j)=frc(:,nsc,icfg)
              disr(:,tau,j)=dsp(:,nsc,icfg)  !- v2a(atom_sc(nsc)%equilibrium_pos) ! subtract before FTT
           enddo
        enddo
! Fourer transform force and displacement to convert to 3 x N0 x Nggrid format 
        do tau=1,natom_prim_cell
           do al=1,3
               auxr=frcr(al,tau,:)
               call fourier_r2k(auxr,auxg)
               frcg(al,tau,:)=auxg
               auxr=disr(al,tau,:) ! equilibrium has been already subtracted 
               call fourier_r2k(auxr,auxg)
               disg(al,tau,:)=auxg
!              call fourier_r2k(disr(al,tau,:),disg(al,tau,:))
!              call fourier_r2k(frcr(al,tau,:),frcg(al,tau,:))
        !      write(*,*)'frcg&disg(29)=',al,tau,frcg(al,tau,29),disg(al,tau,29)
! call write_out(6,' frcg ',frcg)
! call write_out(6,' disg ',disg)
           enddo
        enddo

!write(*,*)' nggrid,nrgrid=',nggrid,nrgrid
!call write_out(6,' ggrid',ggrid)

! now calculate non-analytical force, fcoulg, in the form -D(k)u(k); requires u(k)
        fcoulg=0
        do tau=1,natom_prim_cell
           asr=0
           do taup=1,natom_prim_cell
              do j=1,nggrid
! Long-range Coulomb dynamical matrix for G-vector labeled by j
! call write_out(ulog,' SUBTRACT: dyn_na ',dyn_na)
                 if (length(ggrid(:,j)).lt.tolerance) then
                    asr=asr+dyn_na(tau,taup,:,:,j)
                 endif
! subtract coulomb
                 fcoulg(:,tau,j)=fcoulg(:,tau,j)-matmul(dyn_na(tau,taup,:,:,j),disg(:,taup,j)) 
              enddo
           enddo
           if(maxval(abs(asr)).gt. tolerance) call write_out(6,' asr should be zero ',asr)
! check asr on added term
!          do i=1,3
!          do j=1,3
!             if(abs(asr(i,j)).gt.tolerance)  write(ulog,*)'ASR in dyn_na broken ',tau,i,j,asr(i,j)
!          enddo
!          enddo
        enddo

! fourier transform back to get fcoulr
        do tau=1,natom_prim_cell
        do al=1,3
            auxg=fcoulg(al,tau,:)
            call fourier_k2r(auxg,auxr)
            fcoulr(al,tau,:)=auxr
    !       auxg=disg(al,tau,:)
    !       call fourier_k2r(auxg,auxr)
    !       disr(al,tau,:)=auxr
!           call fourier_k2r(frcg(al,tau,:),frcr(al,tau,:))
!           call fourier_k2r(disg(al,tau,:),disr(al,tau,:))
        enddo
        enddo
!       call write_out(ulog,' Coulomb force',fcoulr)

! from grid indices go back to supercell indices
        do j=1,nrgrid
           rr(1)=(rgrid(:,j) .dot. g01)/2/pi
           rr(2)=(rgrid(:,j) .dot. g02)/2/pi
           rr(3)=(rgrid(:,j) .dot. g03)/2/pi
           n3=nint(rr)
           do tau=1,natom_prim_cell
              call findatom_sc(n3,tau,nsc)
! write into the original array (no need for dsp since it was not changed)
              frc(:,nsc,icfg)=frc(:,nsc,icfg)-fcoulr(:,tau,j)
           enddo
        enddo

     enddo config

  elseif (born_flag.eq.4) then

! calculate directly the non-analytical part in real space , and subtract from the forces
     write(ulog,*)'NA term in real space on the Rgrid: tau,al;taup,be,phi_NA(Rgrid)='
     do tau=1,natom_prim_cell
     do taup=1,natom_prim_cell
        do j=1,nggrid
           call dyn_coulomb(tau,taup,ggrid(:,j),dyn_na(tau,taup,:,:,j),ddn(:,:,:,j)) ! should be zero for all ggrids but zero
        enddo
        do al=1,3
        do be=1,3
           call fourier_k2r(dyn_na(tau,taup,al,be,:),phi_na(tau,taup,al,be,:))
        enddo
        enddo
        write(ulog,8)'tau,al;taup,be,phi_NA=',tau,al,taup,be,phi_na(tau,taup,al,be,:) ! this constant is added phi_0,R for all R
     enddo
     enddo

     config3: do icfg=1,ncfg  

! convert dsp in supercell to disr on the WS inner grid
        do j=1,nrgrid
           rr(1)=(rgrid(:,j) .dot. g01)/2/pi
           rr(2)=(rgrid(:,j) .dot. g02)/2/pi
           rr(3)=(rgrid(:,j) .dot. g03)/2/pi
           n3=nint(rr)
           if(length(rr-n3).gt.1d-4) then
               write(*,5)' SUBTRACT_COULOMB:n3.ne.rr ',n3,rr
               stop
           endif
           do tau=1,natom_prim_cell
              call findatom_sc(n3,tau,nsc)
              disr(:,tau,j)=dsp(:,nsc,icfg)  !- v2a(atom_sc(nsc)%equilibrium_pos) ! subtract equilibrium
           enddo
        enddo

! Non-analytical Coulomb force in real space on Rgrid
        fcoulr=0
        do j=1,nrgrid
        do tau=1,natom_prim_cell
        do al=1,3
           do l=1,nrgrid
              rr=rgrid(:,j)-rgrid(:,l)
! find the corresponding rgrid number
              call fold_in_bz_new(rr,rws26,foldedr)
! jj corresponds to rgrid(j)-rgrid(l) folded back into the supercell
              call find_in_array(foldedr,nrgrid,rgrid,jj,tolerance)

              do taup=1,natom_prim_cell
              do be=1,3
                 fcoulr(al,tau,j)=fcoulr(al,tau,j)- disr(be,taup,l)*phi_na(tau,taup,al,be,jj)
              enddo
              enddo
           enddo
        enddo
        enddo
        enddo

! find correspondance between rgrid point and (tau,R)
        do j=1,nrgrid
           rr(1)=(rgrid(:,j) .dot. g01)/2/pi  ! find first the reduced coordinates
           rr(2)=(rgrid(:,j) .dot. g02)/2/pi
           rr(3)=(rgrid(:,j) .dot. g03)/2/pi
           n3=nint(rr)
           do tau=1,natom_prim_cell
              call findatom_sc(n3,tau,nsc)

              frc(:,nsc,icfg)=frc(:,nsc,icfg)-fcoulr(:,tau,j)
           enddo
        enddo

     enddo config3

  endif

4 format(i4,9(1x,g11.4))
5 format(a,3(i4),3(1x,f10.5))
8 format(a,2(i3,i2),99(1x,f10.4))


 end subroutine subtract_coulomb_force
!--------------------------------------------------
 function hfunc(v,y) result(h)
 use born , only : epsinv
 implicit none
 real(r15), intent(in) :: v(3),y
 integer al,be
 real(r15) h(3,3),t1,t2,my_erfc
 
 t1=my_erfc(y)/(y*y*y)
 t2=2/sqrt(pi)*exp(-y*y)/(y*y)
 do al=1,3
 do be=1,3
    h(al,be)=v(al)*v(be)/(y*y)*(3*t1+t2*(3+2*y*y))
 enddo
 enddo
 h=h-epsinv*(t1+t2)

 end function hfunc
!--------------------------------------------------
 subroutine ewald_2nd_deriv_hat(q,tau,taup,nr,ng,rgrid,ggrid,etaew,d2ew,d3ew)
!! calculates the second derivative of the Ewald potential to be used in the Non-analytical correction
!! input is the two atoms tau and taup, output is the 3x3 block D^EW_tau,taup(q) which goes to D^NA for q \to 0
! use constants, only : ci
 use params, only : tolerance
 use lattice
 use geometry, only : length
 use born , only : epsil,epsinv
 use atoms_force_constants, only : atom0,atompos
 use ios , only : ulog,write_out
 implicit none
 integer, intent(in) :: tau,taup,nr,ng
 real(r15), intent(in) :: q(3),rgrid(3,nr),ggrid(3,ng),etaew
 real(r15), intent(out) :: d2ew(3,3),d3ew(3,3,3)
 integer igrid,al,be,ga
 real(r15) termg,sd(3),del(3),dd,my_erfc,geg,ep,vd(3),dta(3),dhat(3,3),qpg(3),gscale, hout(3,3)

 dta=atompos(:,tau)-atompos(:,taup)
 gscale=6.28/volume_r0**0.3333

! real space term  
    dhat=0 ; d3ew=0
    do igrid=1,nr
       sd=etaew*(dta+rgrid(:,igrid))
       del=matmul(epsinv,sd)
       dd=sqrt(dot_product(sd,del))
       if(dd.lt.1d-6) cycle ! exclude self-interactions
       if (dd.gt.6) cycle
       hout= hfunc(del,dd)
!      write(6,4)'rgrid,dd,hout=',igrid,dd,hout(1,1)
       dhat=dhat - hout*cos(dot_product(q,rgrid(:,igrid))) 
       do ga=1,3
          d3ew(:,:,ga)=d3ew(:,:,ga) - hout*sin(dot_product(q,rgrid(:,igrid))) *rgrid(ga,igrid) 
       enddo
!      write(*,4)'R_sum: igrid,d,dhat_11=',igrid,dd, dhat(1,1)
    enddo
     if(tau.eq.taup) dhat=dhat - 4/3d0/sqrt(pi)*epsinv 
!    if(tau.eq.taup) d3ew=d3ew - 4/3d0/sqrt(pi)*epsinv   !!! TO BE CHECKED !!!
    dhat=dhat* etaew*etaew*etaew/sqrt(det(epsil)) 
    d3ew=d3ew* etaew*etaew*etaew/sqrt(det(epsil)) 
!   write(6,*)'Last R-term(1,1)=',hout(1,1)*etaew*etaew*etaew/sqrt(det(epsil)) 
!   call write_out(6,'total of Rsum terms ',dhat)

! reciprocal space term  
    do igrid=1,ng
       qpg=ggrid(:,igrid)+q
       if(length(qpg).lt.1d-10*gscale) cycle  ! exclude G+q=0
       geg=dot_product(qpg,matmul(epsil,qpg))
       if (geg.gt.100*etaew*etaew) exit  ! assumes ggrid is sorted 
       termg=exp(-geg/4/etaew/etaew)/geg
!      write(6,4)'ggrid,gg,term=',igrid,sqrt(geg),termg
       do al=1,3
       do be=1,3
          hout(al,be)= termg * qpg(al)*qpg(be)*cos(dot_product(qpg,dta)) 
!! NEED TO ADD THIRD DERIVATIVE !!
       enddo
       enddo
       dhat=dhat + hout*4*pi/volume_r0 
!      if(igrid.eq.1) call write_out(6,'G=0 term of ewald ',hout*180.9557368)
    enddo
!   write(*,3)'final EW2DERIV:dhat=',dhat
!   write(6,*)'Last G-term(1,1)=',hout(1,1)*4*pi/volume_r0 

    dhat=dhat /(4*pi*eps0)* ee*1d10
!   write(*,3)'final EW2DERIV:dhat_11 scaled by 1d10*ee/eps0=',dhat(1,1)

    d2ew=matmul(matmul(atom0(tau)%charge,dhat),transpose(atom0(taup)%charge)) * &
&          exp(ci*(q.dot.(atom0(taup)%equilibrium_pos-atom0(tau)%equilibrium_pos)))  ! new phase convention
 !  do al=1,3
 !  do be=1,3
 !     do a2=1,3
 !     do b2=1,3
 !        d2ew(al,be)=d2ew(al,be)+atom0(tau)%charge(al,a2)*dhat(a2,b2)*atom0(tau)%charge(be,b2)
 !     enddo
 !     enddo
 !  enddo
 !  enddo
     

3 format(a,99(1x,g14.7))
4 format(a,i4,99(1x,g11.4))

 end subroutine ewald_2nd_deriv_hat
!--------------------------------------------------

 end module ewald

