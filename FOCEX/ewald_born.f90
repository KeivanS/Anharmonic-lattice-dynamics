!===========================================================
! module born: ../common/born.f90

!=====================================

 module ewald_forces
!! FOCEX: Ewald forces and energies (subtracted from the DFT forces for born_flag 2,3,...);
!! the grids and the dynamical-matrix sums are in module ewald (../common/ewald.f90)
 use constants, only : r15,pi,ci
 use ios, only : ulog
 use ewald
 real(r15), allocatable :: dfc_hyb(:,:,:,:)   ! Phi_C - Phi_m in the current supercell (hybrid_harmonic); set per call
 contains
 subroutine dewapot(x,dpot)
!! calculates the -d/dx(sum_R 1/|R+x| -background) = force  with the corrected metric
!! Assumes translation vectors are r_ewald(:,igrid),g_ewald(:,igrid)
! use constants
 use lattice
 use geometry, only : det
 use params
 use born , only : epsil,epsinv
 implicit none
 real(r15), intent(out) :: dpot(3)
 real(r15), intent(in) :: x(3)
 integer igrid
 real(r15) termg,qpg(3),dd,my_erfc,geg,vd(3)

 if (length(x).lt.1d-12) then
    write(*,*)'DEWAPOT: X is too small! ',x
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

! write(*,3)'dpotr=',eta,dd,dpot

    do igrid=2,ng_ewald
       qpg=g_ewald(:,igrid)
       geg=dot_product(qpg,matmul(epsil,qpg))
       if (geg.gt.100*eta*eta) exit  ! exiting the loop because g_ewald is sorted
       termg=exp(-geg/4/eta/eta)/geg
       dpot=dpot + termg * qpg*sin(qpg .dot. x) *4*pi/volume_r  ! this is supercell volume
    enddo

! write(*,3)'termg,dpotg=',termg,dpot

!   dpot=dpot + 4*pi/volume_r*matmul(epsinv,x)/3d0

! write(*,3)'dpot=',dpot
3 format(a,9(1x,g11.4))

 end subroutine dewapot
!===========================================================
 subroutine d2ewapot(x,hpot)
!! hpot = grad grad phi(x) of the potential whose -grad is dewapot (same images, same cutoffs):
!! real space erfc(eta D)/(D sqrt det eps), D = sqrt(x.eps^-1.x); reciprocal (4pi/V) sum_G exp(-GeG/4eta^2)/GeG cos(G.x)
 use lattice
 use geometry, only : det
 use params
 use born , only : epsil,epsinv
 implicit none
 real(r15), intent(in) :: x(3)
 real(r15), intent(out) :: hpot(3,3)
 integer igrid,a
 real(r15) r(3),u(3),dd,d,e,g,g1,g2,qpg(3),geg,termg,ei(3,3),sd,my_erfc
 ei=(epsinv+transpose(epsinv))/2; sd=sqrt(det(epsil))
 hpot=0
 do igrid=1,nr_ewald
    r=x+r_ewald(:,igrid)
    u=matmul(ei,r)
    d=sqrt(dot_product(r,u)); dd=eta*d
    if(dd.lt.1d-12) cycle
    if (dd.gt.6) cycle
    e=my_erfc(dd); g=2*eta/sqrt(pi)*exp(-dd*dd)
    g1=-e/d**2-g/d
    g2=2*e/d**3+g*(2/d**2+2*eta*eta)
    do a=1,3
       hpot(:,a)=hpot(:,a)+(g2-g1/d)*u*u(a)/d**2/sd
       hpot(:,a)=hpot(:,a)+g1/d*ei(:,a)/sd
    enddo
 enddo
 do igrid=2,ng_ewald
    qpg=g_ewald(:,igrid)
    geg=dot_product(qpg,matmul(epsil,qpg))
    if (geg.gt.100*eta*eta) exit
    termg=exp(-geg/4/eta/eta)/geg
    do a=1,3
       hpot(:,a)=hpot(:,a) - termg*qpg*qpg(a)*cos(dot_product(qpg,x))*4*pi/volume_r
    enddo
 enddo
 end subroutine d2ewapot
!===========================================================
 subroutine hybrid_harmonic(n,dfc)
!! dfc(:,:,i,j) = Phi_C(i,j) - Phi_m(i,j) (eV/A^2) in the supercell at the reference: Phi_C = -e^2 Z_i H_ij Z_j^T the
!! Coulomb dipole FCs (the D_NA added back to D), Phi_m = -e^2 c_ij H_ij the harmonic part of ewaldforce's pair model,
!! H_ij = grad grad phi(r_ij) (d2ewapot, images summed); on-site blocks by the sum rule. The born_flag 2 subtraction
!! uses F_m(x) - F_m(x0) - dfc u: harmonic part = Phi_C exactly, anharmonic part from the pair model (user, 2026-11-03).
 use constants, only : ee, eps0,pi
 use atoms_force_constants
 use geometry, only : v2a
 use ios
 implicit none
 integer, intent(in) :: n
 real(r15), intent(out) :: dfc(3,3,n,n)
 integer i,j
 real(r15) h(3,3),zi(3,3),zj(3,3),e2,r(3),hf(3,3),dp(3),dm(3),err,hs
 e2=ee*1d10/(4*pi*eps0)
 dfc=0
 do i=1,n
    zi=atom_sc(i)%charge
    do j=1,n
       if (j.eq.i) cycle
       zj=atom_sc(j)%charge
       r=v2a(atom_sc(i)%equilibrium_pos)-v2a(atom_sc(j)%equilibrium_pos)
       call d2ewapot(r,h)
       dfc(:,:,i,j)=-e2*(matmul(zi,matmul(h,transpose(zj))) - sum(zi*zj)/3*h)
    enddo
 enddo
 do i=1,n
    do j=1,n
       if (j.ne.i) dfc(:,:,i,i)=dfc(:,:,i,i)-dfc(:,:,i,j)
    enddo
 enddo
! check of d2ewapot against differences of dewapot (one pair)
 r=v2a(atom_sc(1)%equilibrium_pos)-v2a(atom_sc(min(2,n))%equilibrium_pos)
 if (n.gt.1) then
    call d2ewapot(r,h)
    hs=1d-4
    do i=1,3
       r(i)=r(i)+hs; call dewapot(r,dp); r(i)=r(i)-2*hs; call dewapot(r,dm); r(i)=r(i)+hs
       hf(:,i)=-(dp-dm)/(2*hs)
    enddo
    err=maxval(abs(hf-h))/max(maxval(abs(h)),1d-30)
    write(ulog,'(a,es10.2)')' HYBRID_HARMONIC: d2ewapot vs differences of dewapot (relative) =',err
 endif
 write(ulog,'(a,es12.4,a)')' HYBRID_HARMONIC: max|Phi_C - Phi_m| =',maxval(abs(dfc)),' eV/A^2'
 end subroutine hybrid_harmonic
!===========================================================
 subroutine ewaldforce(n,pos,fewald)
!! calculates the Coulomb force on atoms in supercell (in 1/Ang)
!! Pair potential c_12 phi_eps(r_12), c_12 = tr(Z_1 Z_2^T)/3 (2026-11-03, user): a scalar coupling, so the force
!! derives from an energy and its direction is that of -grad phi; for Z = z 1 it is z_1 z_2 (point charges).
!! Before: Z_1 Z_2 (-grad phi), whose FCs were not symmetric for anisotropic Z.
!! Assumes translation vectors are r_ewald(:,igrid),g_ewald(:,igrid)
 use constants, only : ee, eps0,pi
 use lattice
 use params
 use atoms_force_constants
 use born
 use ios
 implicit none
 integer, intent(in) :: n  ! number of atoms in the supercell will be used
 real(r15), intent(in) :: pos(3,n) ! displaced positions including equilibrium positions
 real(r15), intent(out):: fewald(3,n)
 integer tau1,tau2
 real(r15) dpot(3),r12(3),z1(3,3),z2(3,3),efield(3) 


 fewald=0
 do tau1=1,n  ! caculate the total force on tau1
    z1=atom_sc(tau1)%charge
    efield=0
 do tau2=1,n
    z2=atom_sc(tau2)%charge

    if (tau1.eq.tau2) cycle ! images of tau1 have zero force on tau1 because of inversion symmetry
    r12=pos(:,tau1)-pos(:,tau2) ! pos includes equilibrium positions
! force on tau1 from tau2 and all its images inluding epsilon; beware of the sign!
    call dewapot(r12,dpot)  ! has epsil-modified metric included
    efield = efield + sum(z1*z2)/3 * dpot      ! tr(Z1 Z2^T)/3
 enddo
    fewald(:,tau1)= efield * ee*1d10/(4*pi*eps0)  ! to convert to eV/Ang
 enddo

! now compare to the energy finite difference
! write(ulog,*)'******** TEST OF COULOMB FORCE ***********'
! allocate(dsp1(3,n))
! dsp1=dsp
! do tau1=1,n
!    call coulomb_energy(tau1,n,atom_sc,dsp,ecoul0) !assumes all other atoms are @eequilibrium
!    do al=1,3
!       dsp1(al,tau1)=dsp(al,tau1)+1d-5
!       call coulomb_energy(tau1,n,atom_sc,dsp1,ecoul1)
!       force1=(ecoul0-ecoul1)/1d-5
!       write(ulog,3)'al,tau,force,-de/dx=',al,tau1,frc(al,tau1),force1
!       write(*,3)'force=',tau1,tau1,fewald(:,tau1)
!    enddo
! enddo
! deallocate(dsp1)

3 format(a,2i5,9(1x,g11.4))

 end subroutine ewaldforce
!===========================================================
 subroutine coulombforce  (n,pos,frc,dsp) 
!! calculates the Coulomb force on atoms in supercell (in 1/Ang)
!! Assumes translation vectors are rgrid(:,nrgrid),ggrid(:,nggrid)
 use constants, only : eps0scale,pi
 use lattice
 use params
 use atoms_force_constants
 use born
 use ios
 use fourier, only : nggrid,ggrid,nrgrid,rgrid,rws_weights,gws_weights,fourier_k2r,fourier_r2k
 implicit none
 integer, intent(in) :: n  ! number of atoms in the supercell will be used
 real(r15), intent(in) :: pos(3,n),dsp(3,n) ! displaced positions including equilibrium positions
 real(r15), intent(out):: frc(3,n)
 integer tau1,tau2,igrid
 real(r15) dpot(3),r12(3),z1(3),z2(3),efield(3) ,geg


 frc=0
 do tau1=1,n  ! caculate the total force on tau1
    efield=0
    do tau2=1,n
       if (tau1.eq.tau2) cycle ! images of tau1 have zero force on tau1 because of inversion symmetry
       r12=v2a(atom_sc(tau1)%equilibrium_pos)-v2a(atom_sc(tau2)%equilibrium_pos)
       g0loop: do igrid=2,nggrid ! exclude G=0
          z2=matmul(atom_sc(tau2)%charge,ggrid(:,igrid))
          geg=dot_product(ggrid(:,igrid),matmul(epsil,ggrid(:,igrid)))
          dpot=dpot + z2/geg *sin(ggrid(:,igrid) .dot. r12) * gws_weights(igrid) 
          efield = efield + dpot
       enddo g0loop
    enddo
    frc(:,tau1)= matmul(atom_sc(tau1)%charge,efield) *coef  ! this is force for u=0

    do tau2=1,n
       r12=pos(:,tau2)-v2a(atom_sc(tau2)%equilibrium_pos)-(pos(:,tau1)-v2a(atom_sc(tau1)%equilibrium_pos))
       g1loop: do igrid=2,nggrid ! exclude G=0
          geg=dot_product(ggrid(:,igrid),matmul(epsil,ggrid(:,igrid)))
          z1=matmul(atom_sc(tau1)%charge,ggrid(:,igrid))
          z2=matmul(atom_sc(tau2)%charge,ggrid(:,igrid))
          dpot=dpot + z2/geg *sin(ggrid(:,igrid) .dot. r12) * gws_weights(igrid) 
          efield = efield + dpot
       enddo g1loop
    enddo

 enddo

3 format(a,99(1x,g14.7))

 end subroutine coulombforce
!===========================================================
 subroutine ewapot(x,pot)
! calculates the sum_R 1/|R+x| -background = pot (in 1/Ang no epsil)
! use constants , only: pi
 use lattice , only: volume_r
 use geometry, only : det,length
 use born
 implicit none
 real(r15), intent(out) :: pot
 real(r15), intent(in) :: x(3)
 integer igrid
 real(r15) termg,qpg(3),dd,my_erfc,geg,vd(3)

 if (length(x).lt.1d-12) then
    write(*,*)'EWAPOT: X is too small! ',x
 endif

    pot=0
    do igrid=1,nr_ewald
       vd=eta*(x+r_ewald(:,igrid))
       dd=sqrt(dot_product(vd,matmul(epsinv,vd)))
       if(dd.lt.1d-12) dd=1d-12
       if (dd.gt.5) cycle  ! r_ewald is sorted but not x+r_ewald
       pot=pot + eta*my_erfc(dd)/dd/sqrt(det(epsil))
    enddo

    do igrid=2,ng_ewald
       qpg=g_ewald(:,igrid)
       geg=dot_product(qpg,matmul(epsil,qpg))
       if (geg.gt.100*eta*eta) exit
       termg=exp(-geg/(4*eta*eta))/geg
       pot=pot+ termg * cos(dot_product(qpg,x)) *4*pi/volume_r
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
 use geometry, only : det
 implicit none
 real(r15), intent(out) :: emadelung
 real(r15) x(3),ene 

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

 r12=(/0.1,0.2,0.5/)  
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
 stop
3 format(a,i5,9(1x,g14.7))

 end subroutine test_ewald
!===========================================================
 subroutine subtract_coulomb_force(born_flag,ncfg,dsp,frc)
!! subtracts the non-analytical long-range contribution of Coulomb forces from total
!! force to get and fit the short-range contribution. Then this NA part is added back
!! to the fitted short-range dynamical matrix in subroutine dyn_coulomb
! dsp is the cartesian position not including the equilibrium positions
 use atoms_force_constants
! use constants, only : pi
 use fourier, only : nggrid,ggrid,nrgrid,rgrid,rws_weights,fourier_k2r,fourier_r2k,gws_weights
 use geometry, only : length,v2a
 use lattice, only : rs1,rs2,rs3,r01,r02,r03,rws26,gws26,volume_r,volume_r0,a0_scale,asc_scale,g0ws26,fold_ws,cart2red
 use ios, only : ulog ,write_out
 use params, only : verbose, tolerance,coef
 use born , only : epsil,epsinv,dyn_na
 implicit none
 integer, intent(in):: born_flag,ncfg
 real(r15), intent(inout) :: dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg)
! automatic arrays here (lost after the call, memory released)
 integer nsc(natom_prim_cell,nrgrid)
!  complex(r15) frcr(3,natom_prim_cell,nrgrid),frcg(3,natom_prim_cell,nggrid)
 complex(r15) disr(3,natom_prim_cell,nrgrid),disg(3,natom_prim_cell,nggrid)
 complex(r15) fcoulg(3,natom_prim_cell,nggrid),fcoulr(3,natom_prim_cell,nrgrid)
 complex(r15) auxr(nrgrid),auxg(nggrid),aravg(nrgrid) ,agavg(nggrid) 
 real(r15) coulmb_force(3,natom_super_cell),coulmb_force0(3,natom_super_cell)
! real(r15), allocatable :: frcr(:,:,:),frcg(:,:,:),disr(:,:,:),disg(:,:,:)
 complex(r15) ddn(3,3,3,nggrid)
 real(r15) rr(3),deteps3,asr(3,3),pos(3,natom_super_cell) ,etaew,ftot,dfewa,dta(3),geg,rfold(3)
 real(r15) phi_na(natom_prim_cell,natom_prim_cell,3,3,nrgrid)
 integer n3(3),tau,taup,al,be,icfg,i,g,l,j,nsc1,ir,jr
 logical isperiodic

 if(born_flag.le.1 ) return  ! no subtraction for born_flag=1

 write(ulog,*)'born_flag=',born_flag
 if(born_flag.eq.1) write(ulog,*)'no subtraction but D^NA(Parlinski) is added to dyn_SR'
 if(born_flag.eq.2) write(ulog,*)'subtraction of deltaF(ewald) done from force, D^NA(_hat) added to dyn_SR'
 if(born_flag.eq.3) write(ulog,*)'subtraction of - D^NA_K(_hat) u_K from force, D^NA(_hat) added to dyn_SR'
 if(born_flag.eq.4) write(ulog,*)'subtraction of - D^NA_K(Gnew) u_K from force, D^NA(Gnew) added to dyn_SR'
 if(born_flag.eq.6) write(ulog,*)'no subtraction here but D^NA_struct_factor added '
 if(born_flag.eq.7) write(ulog,*)'subtraction of fourier summation of coulomb force, D^NA(_hat) added to dyn_SR'
 if(born_flag.eq.8) write(ulog,*)'subtraction of -phi_NA_t,Rtp *u_Rtp , D^NA(_hat) added to dyn_SR'
 if(born_flag.eq.9) write(ulog,*)'subtraction of force from Fourier of Coulomb (no Ewald) D^NA(_hat) added to dyn_SR'


    write(*,*)'SUBTRACT_COULOMB_FORCE: sizes of force are:3,',  &
&            size(frc(1,:,1)),size(frc(1,1,:))

! ------------------------------  SOME PRELIMINARY CALCULATIONS  ----------------------------------
! erfc(4)=1.5d-8  so if eps=4/vol**.333 the real space terms become negligible
! exp(-18)=1.5d-8  so G larger than 5-6 shells will not contribute even if eps=4/vol^.3
! set cutoffs for Ewald sums ! this should come after the supercell is read!!!!
! to subtract coulomb forces from ewald sums, the superccell translation vectors are needed
! output is r_ewald(3,nr_ewald) and g_ewald(3,ng_ewald)  for ewald summations
  deteps3=det(epsil)**(1/3d0)
  asc_scale = (volume_r **(1/3d0))
  a0_scale  = (volume_r0**(1/3d0))
  if ( born_flag.eq.2) then
     call cuts(deteps3,asc_scale,born_flag,eta,rcutoff_ewa,gcutoff_ewa)
     call make_grid_shell_ewald(rs1,rs2,rs3,rcutoff_ewa,gcutoff_ewa,epsil)
  else
     call cuts(deteps3, a0_scale,born_flag,eta,rcutoff_ewa,gcutoff_ewa)
     call make_grid_shell_ewald(r01,r02,r03,rcutoff_ewa,gcutoff_ewa,epsil)
  endif

  write(ulog,6)'eta,rcut_ewa,gcut_ewa=',eta,rcutoff_ewa,gcutoff_ewa
  write(6   ,6)'eta,rcut_ewa,gcut_ewa=',eta,rcutoff_ewa,gcutoff_ewa
  write(ulog,*)'nr_ewa,ng_ewa=',nr_ewald,ng_ewald
  write(6   ,*)'nr_ewa,ng_ewa=',nr_ewald,ng_ewald

 if(born_flag.gt.9) return  ! only subtract from forces for bflag=2,3,4; 

! mapping from (rgrid,tau) to supercell: nsc(tau,jgrid)
  call nsc_from_rgrid(nrgrid,rgrid,nsc) 

! calculate non-analytical part DYN_NA(q) and subtract from dny_na(G*) ; no mass factor
!
! hybrid (born_flag 2): Phi_C - Phi_m of THIS supercell, computed on every call (one call per POSCARi/FORCEDISPi
! pair). Kept from the first call until 2026-11-03, a second supercell read past the end of the first one's
! array (GaN, 24- and 72-atom supercells: garbage forces, fit errors changing from run to run).
  if (born_flag.eq.2) then
     if (allocated(dfc_hyb)) deallocate(dfc_hyb)
     allocate(dfc_hyb(3,3,natom_super_cell,natom_super_cell))
     call hybrid_harmonic(natom_super_cell,dfc_hyb)
  endif

! ----------------   NOW DO SUBTRACTION FOR EVERY CONFIGURATION    ----------------
  config_sum: do icfg=1,ncfg

      if( born_flag.eq.2 .or. born_flag.eq.9) then ! subtract ewald force from forces

         write(ulog,*)'MAIN: going to call Ewaldforce '

         do i=1,natom_super_cell
            pos(:,i)=v2a(atom_sc(i)%equilibrium_pos)
         enddo 
!        if( born_flag.eq.2 ) then ! subtract ewald force from forces
            call ewaldforce  (natom_super_cell,pos,coulmb_force0) ! is it necessary to subtract?
!        elseif( born_flag.eq.9 ) then ! subtract coulomb force from forces
!           call coulombforce  (natom_super_cell,pos,coulmb_force0) ! is it necessary to subtract?
!        endif
! yes because the residual force pi0 is also subtrated from forces when fitting. so this should be linear in u

         do i=1,natom_super_cell
            pos(:,i)=v2a(atom_sc(i)%equilibrium_pos) + dsp(:,i,icfg)
         enddo 
!        if( born_flag.eq.2 ) then ! subtract ewald force from forces
            call ewaldforce  (natom_super_cell,pos,coulmb_force) ! is it necessary to subtract?
!        elseif( born_flag.eq.9 ) then ! subtract coulomb force from forces
!           call coulombforce  (natom_super_cell,pos,coulmb_force) ! is it necessary to subtract?
!        endif

!        if ( verbose ) then
            write(ulog,*)' Ewald force F,F0,F-F0,dF_ewa/Ftot for configuration #',icfg
            dfewa=0; ftot=0
            do i=1,natom_super_cell
               write(ulog,3)i, coulmb_force(:,i),coulmb_force0(:,i), &
 &                             coulmb_force(:,i)-coulmb_force0(:,i), &
 &            length(coulmb_force(:,i)-coulmb_force0(:,i))/(1d-12+length(frc(:,i,icfg)))
               dfewa=dfewa+length(coulmb_force(:,i)-coulmb_force0(:,i))
               ftot=ftot+length(frc(:,i,icfg))
            enddo
!        endif
         write(ulog,5)'######## sum(Ftot),sum(dfewa),Sdfewa/SFtot=',ftot,dfewa,dfewa/ftot
         coulmb_force = coulmb_force- coulmb_force0 
! hybrid (born_flag 2): replace the harmonic part of the pair model by the Coulomb dipole one, - (Phi_C - Phi_m) u
         if (born_flag.eq.2) then
            do i=1,natom_super_cell
               do j=1,natom_super_cell
                  coulmb_force(:,i)=coulmb_force(:,i)-matmul(dfc_hyb(:,:,i,j),dsp(:,j,icfg))
               enddo
            enddo
         endif
         write(982,*)' nsc,rgrid,tau, coulmb_force on Rgrid (redundant) ====== BF, nrgrid=',born_flag,nrgrid
         do j=1,nrgrid
         do tau=1,natom_prim_cell
            nsc1=nsc(tau,j)
            write(982,4) nsc1,j,tau, coulmb_force(:,nsc1)
         enddo
         enddo

         frc(:,:,icfg)=frc(:,:,icfg) -  coulmb_force(:,:) 

      elseif(born_flag.eq.3 .or. born_flag.eq.4 ) then   ! subtract FT of -D_NA(t,t'G) * u(t',G)

         write(ulog,*)'SUBTRACT_COULOMB: calculating dyn_na on ggrid for bflag=',born_flag
         call dyn_na5D_on_grid(nggrid,ggrid)  ! this gives the 5D dyn_na(G*), which satisfies ASR 

! for each config convert disp and force to a 3 x N0 x Nrgrid format needed for FT
         do j=1,nrgrid
         do tau=1,natom_prim_cell
            nsc1=nsc(tau,j)
            disr  (:,tau,j)=dsp(:,nsc1,icfg)  ! should be periodic if j on boundary
 !          fcoulr(:,tau,j)=frc(:,nsc1,icfg) 
         enddo
         enddo

! Fourier transform displacement and force to convert to 3 x N0 x Nggrid format 
         do tau=1,natom_prim_cell
         do al=1,3
            auxr=disr(al,tau,:)           
            call check_periodic(nrgrid,rgrid,rws26,auxr,aravg,isperiodic)
            if(isperiodic) then
              call fourier_r2k(auxr,auxg)   
            else
              write(ulog,*)'icfg=',icfg,' disr was not periodic!!!'
              stop
            ! call fourier_r2k(aravg,auxg)  
            endif 
            disg(al,tau,:)=auxg           

         enddo
         enddo

! subtract  -D_NA(t,t',G)*u(t',G) from fcoulg(t,G)
         fcoulg = 0
         do tau=1,natom_prim_cell
         do al=1,3
            do taup=1,natom_prim_cell
            do be=1,3
               fcoulg(al,tau,:)=fcoulg(al,tau,:) - dyn_na(tau,taup,al,be,:)*disg(be,taup,:) 
            enddo
            enddo
         enddo
         enddo

! fourier transform back to real space and get fcoulr
         do tau=1,natom_prim_cell
         do al=1,3
            auxg=fcoulg(al,tau,:)    
            call check_periodic(nggrid,ggrid,g0ws26,auxg,agavg,isperiodic)
            if(isperiodic) then
               call fourier_k2r(auxg,auxr)  
            else
               write(ulog,*)'icfg=',icfg,' fcoulg was not periodic!!!',tau,al
               call write_out(ulog,'fcoulg was not periodic!!!',auxg)
               call write_out(ulog,'fcoulg was not periodic!!!',agavg)
               call write_out(ulog,'fcoulg-agavg ',agavg-auxg)
               stop
     !         call fourier_k2r(agavg,auxr)   
            endif 
            fcoulr(al,tau,:)=auxr / rws_weights  ! if on boundary, still should count with full weight   
         enddo
         enddo

! go back from grid indices to supercell indices
         write(983,*)' nsc,rgrid,tau, NA_force(:,nsc)========= BF, nrgrid=',born_flag,nrgrid
         do j=1,nrgrid
         do tau=1,natom_prim_cell
            nsc1=nsc(tau,j)
            write(983,4) nsc1,j,tau,real(fcoulr(:,tau,j)), aimag(fcoulr(:,tau,j))
            frc(:,nsc1,icfg)=frc(:,nsc1,icfg)-real(fcoulr(:,tau,j))
            if(maxval(abs(aimag(fcoulr(:,tau,j)))).gt.1d-6) then
               write(*   ,5)'ERROR:tau,rgrid,nsc,fcoulr has IMAGINARY part!='  &
 &                                ,tau,j,nsc1,aimag(fcoulr(:,tau,j))
               write(ulog,5)'ERROR:tau,rgrid,nsc,fcoulr has IMAGINARY part!='  &
 &                                ,tau,j,nsc1,aimag(fcoulr(:,tau,j))
               stop
            endif
         enddo
         enddo

      elseif (born_flag.eq.7) then

! calculate the Force using fourier summation and subtract from the force
         fcoulr =0
         do tau=1,natom_prim_cell
         do ir=1,nrgrid
            do taup=1,natom_prim_cell
               dta=rgrid(:,ir)+v2a(atom0(taup)%equilibrium_pos) - v2a(atom0(tau)%equilibrium_pos) 
            do jr=2,nggrid
               geg=dot_product(ggrid(:,jr),matmul(epsil,ggrid(:,jr)))
               rr=matmul(atom0(taup)%charge,ggrid(:,jr)) / geg *   &
   &                 sin(ggrid(:,jr).dot.dta) * gws_weights(jr) * coef
               fcoulr(:,tau,ir)=fcoulr(:,tau,ir)+matmul( atom0(tau)%charge,rr)
            enddo
            enddo
         enddo
         enddo

         write(987,*)' nsc,rgrid,tau, Coulmb_force (BF=8) on Rgrid (redundant) ==== nrgrid=',nrgrid
         do ir=1,nrgrid
         do tau=1,natom_prim_cell
            nsc1=nsc(tau,ir)
            write(987,4) nsc1,ir,tau,real(fcoulr(:,tau,ir)), aimag(fcoulr(:,tau,ir))
            frc(:,nsc1,icfg)=frc(:,nsc1,icfg)-real(fcoulr(:,tau,ir))
            if(maxval(abs(aimag(fcoulr(:,tau,ir)))).gt.1d-6) then
               write(*   ,5)'ERROR:tau,rgrid,nsc,fcoulr has IMAGINARY part!='  &
 &                                ,tau,ir,nsc1,aimag(fcoulr(:,tau,ir))
               write(ulog,5)'ERROR:tau,rgrid,nsc,fcoulr has IMAGINARY part!='  &
 &                                ,tau,ir,nsc1,aimag(fcoulr(:,tau,ir))
               stop
            endif
         enddo
         enddo

      elseif (born_flag.eq.8) then

! calculate the phi_na using fourier summation and subtract from the force
         call phi_na_g2r(phi_na,nrgrid,rgrid,nggrid,ggrid,gws_weights)
         fcoulr =0
         do tau=1,natom_prim_cell
         do al=1,3
            do ir=1,nrgrid
            do jr=1,nrgrid
               rr=rgrid(:,ir)-rgrid(:,jr)
! fold rr back in the rgrid
               rfold = fold_ws(rr,rws26) 
! find the index of rfold
               i=0
               lloop: do l=1,nrgrid
                  if(length(rgrid(:,l)-rfold).lt.tolerance) then
                     i=l
                     exit lloop
                  endif
               enddo lloop
               if(i.eq.0) then
                  write(*,*)'SUBTRACT: BF=8, i=0, rfold not found ',cart2red(rfold,'r')
                  stop
               endif
               do taup=1,natom_prim_cell
               do be=1,3
          fcoulr(al,tau,ir)=fcoulr(al,tau,ir)-phi_na(tau,taup,al,be,i)*disr(be,taup,jr) 
               enddo
               enddo
            enddo
            enddo
         enddo
         enddo

         write(988,*)' nsc,rgrid,tau, Coulmb_force (BF=8) on Rgrid (redundant) ==== nrgrid=',nrgrid
         do ir=1,nrgrid
         do tau=1,natom_prim_cell
            nsc1=nsc(tau,ir)
            write(988,4) nsc1,ir,tau,real(fcoulr(:,tau,ir)), aimag(fcoulr(:,tau,ir))
            frc(:,nsc1,icfg)=frc(:,nsc1,icfg)-real(fcoulr(:,tau,ir))
            if(maxval(abs(aimag(fcoulr(:,tau,ir)))).gt.1d-6) then
               write(*   ,5)'ERROR:tau,rgrid,nsc,fcoulr has IMAGINARY part!='  &
 &                                ,tau,ir,nsc1,aimag(fcoulr(:,tau,ir))
               write(ulog,5)'ERROR:tau,rgrid,nsc,fcoulr has IMAGINARY part!='  &
 &                                ,tau,ir,nsc1,aimag(fcoulr(:,tau,ir))
               stop
            endif
         enddo
         enddo

      endif

  enddo config_sum
  if (allocated(dfc_hyb)) deallocate(dfc_hyb)

3 format(i4,99(1x,f9.4))
4 format(3i4,99(2x,3(1x,g10.3)))
5 format(a,33(1x,g11.4))
6 format(a,9(1x,f10.5))
8 format(a,2(i3,i2),999(1x,f9.4))

 end subroutine subtract_coulomb_force
!--------------------------------------------------
 end module ewald_forces

!===========================================================
