 module ewald
 integer nr_ewald,ng_ewald
 real(8), allocatable :: r_ewald(:,:),g_ewald(:,:)

 contains

 subroutine dyn_coulomb(tau,taup,j,dyn)
!! calculates the q component of the Coulomb force to subtract from fourier transform of total force
! sdyn(tau,taup) = z(tau)^al,ga.q^ga  z(taup)^be.de.q^de  &
!    &          1/eps0/(q.epsr.q)/volume_r q * u(qtaup,de)
 use atoms_force_constants
 use lattice, only : volume_r
 use fourier, only : ggrid
 use born, only : epsil
 use constants , only : eps0,ee
 implicit none
 integer, intent(in) :: tau,taup,j
 real(8), intent(out) :: dyn(3,3)  ! element (tau,taup)(ggrid(:,j)) of the Coulomb dynamical mat
 integer al,be
 real(8) zqa(3),zqb(3),denom,q(3),coef

 coef=ee/(eps0*1d10) ! to convert 1/ang^2 to eV/ang , includes cancellation of 4pi
 q=ggrid(:,j) 
 zqa=matmul(atom0(tau)%charge,q)
 zqb=matmul(atom0(taup)%charge,q)
 denom=dot_product(q,matmul(epsil,q))
 
 do al=1,3
 do be=1,3
    dyn(al,be)=zqa(al)*zqb(be)/denom/volume_r*coef
 enddo
 enddo

 end subroutine dyn_coulomb
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
 real(8), intent(in ):: rcut,gcut,epsilo(3,3)
 type(vector), intent(in ):: x1,x2,x3
! integer, intent(inout):: ngrid
! real(8), intent(inout):: grid(:,:)
! real(8), allocatable, intent(inout):: grid(:,:)
 integer i1,i2,i3,m,cnt,max
 real(8) v(3),rr,om0,aux2(3,3),epsinv(3,3)
 real(8), allocatable :: aux(:,:),lengths(:)
 integer, allocatable :: msort(:)
 type(vector) y1,y2,y3,t1,t2,t3

! call apply_metric(x1,x2,x3,epsilo,t1,t2,t3)
 aux2=epsilo
 call inverse_real(aux2,epsinv,3)

 call get_upper_bounds(x1,x2,x3,rcut,max,m)
! call calculate_volume(x1,x2,x3,om0)
! max=nint(15d0/3d0*rcut**3/om0)+10
! m=nint((max/4.2)**0.333)+4

 allocate(aux(3,max),lengths(max),msort(max))
 lengths =1d20; aux=1d20
 write(ulog,4)'MAKE_GRID_SHELL, max(Rgrid), for rcut=',max,rcut
 cnt=1
 do i1=-m,m
 do i2=-m,m
 do i3=-m,m

! generate vectors in a grid
    v= v2a(i1*x1 + i2*x2 + i3*x3)
!    rl=length(v)
    rr=dot_product(v,matmul(epsinv,v))
!    rr=sqrt(dot_product(v,matmul(epsinv,v)))
    if(rr.gt.rcut) cycle

    aux(:,cnt)=v
    lengths(cnt)=rr
    cnt=cnt+1
    if (cnt.gt.max) then
       write(ulog,*)'RLOOP: max size exceeded, need to increase variable max from ',max
       stop
    endif

 enddo
 enddo
 enddo

 nr_ewald=cnt-1
 write(ulog,5) 'MAKE_GRID_SHELL: within rcut=',rcut,' r_ewald generated ',nr_ewald,' vectors'

 if (allocated(r_ewald)) deallocate(r_ewald)  ! overwrite if previously called
 if (allocated(msort)) deallocate(msort)
 allocate(r_ewald(3,nr_ewald),msort(max))

 call sort(nr_ewald,lengths,msort,max)

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
 call get_upper_bounds(y1,y2,y3,gcut,max,m)
 write(ulog,4)'MAKE_GRID_SHELL, max(g_ewald), for gcut=',max,gcut
 allocate(aux(3,max),lengths(max),msort(max))
 lengths =1d20; aux=1d20

 cnt=1
 do i1=-m,m
 do i2=-m,m
 do i3=-m,m

! generate vectors in a grid
    v= v2a(i1*y1 + i2*y2 + i3*y3)
!    rr=length(v)
    rr=sqrt(dot_product(v,matmul(epsilo,v)))
    if(rr.gt.gcut) cycle

    aux(:,cnt)=v
    lengths(cnt)=rr
    cnt=cnt+1
    if (cnt.gt.max .and. max.ne.1) then
       write(ulog,*)'GLOOP: max size exceeded, need to increase variable max from ',max
       stop
    endif

 enddo
 enddo
 enddo

 ng_ewald=cnt-1
 write(ulog,5) 'MAKE_GRID_SHELL: within gcut=',gcut,' g_ewald generated ',ng_ewald,' vectors'

 if (allocated(g_ewald)) deallocate(g_ewald)
 if (allocated(msort)) deallocate(msort)
 allocate(g_ewald(3,ng_ewald),msort(max))

 call sort(ng_ewald,lengths,msort,max)

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
 use constants
 use lattice
 use params
 use born , only : epsil,epsinv
 implicit none
 real(8), intent(out) :: dpot(3)
 real(8), intent(in) :: x(3)
 integer igrid
 real(8) termg,qpg(3),dd,erfc,geg,ep,vd(3)

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
       dpot=dpot + eta*eta*matmul(epsinv,vd)/dd/dd * (erfc(dd)/dd+2/sqrt(pi)*exp(-dd*dd))/sqrt(det(epsil))
    enddo

! write(*,3)'dpotr=',dpot

    do igrid=2,ng_ewald
       qpg=g_ewald(:,igrid)
       geg=dot_product(qpg,matmul(epsil,qpg))
       if (geg.gt.100*eta*eta) exit
       termg=exp(-geg/4/eta/eta)/geg
       dpot=dpot + termg * qpg*sin(qpg .dot. x) *4*pi/volume_r
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
 use constants
 use lattice
 use params
 use atoms_force_constants
 use born
 use ios
 implicit none
 integer, intent(in) :: n  ! number of atoms in the supercell will be used
 real(8), intent(in) :: dsp(3,n) ! displaced positions
 real(8), intent(out):: fewald(3,n)
 integer tau1,tau2 
 real(8) dpot(3),r12(3),z1(3,3),z2(3,3),f1(3),f2(3) ,coef

 coef=ee/(4*pi*eps0)*1d10
 fewald=0
 do tau1=1,n  ! caculate the total force on tau1
    z1=atom_sc(tau1)%charge
 do tau2=1,n
    z2=atom_sc(tau2)%charge

    if (tau1.eq.tau2) cycle ! images of tau1 have zero force on tau1
!    r12=atom_sc(tau1)%equilibrium_pos+dsp(:,tau1)-atom_sc(tau2)%equilibrium_pos-dsp(:,tau2)
    r12=dsp(:,tau1)-dsp(:,tau2)
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
 use constants
 use lattice
 use params
 use geometry
 use born 
 implicit none
 real(8), intent(out) :: pot
 real(8), intent(in) :: x(3)
 integer igrid
 real(8) termg,qpg(3),dd,erfc,geg,ep,vd(3)

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
       pot=pot + eta*erfc(dd)/dd/sqrt(det(epsil))
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
 use constants
! use lattice
! use params
 use atoms_force_constants
! use born
 implicit none
 integer, intent(in) :: i,n
 type(atomic), intent(in) :: atoms(n)
 real(8), intent(in) :: dsp(3,n)
 real(8), intent(out) :: e_coul
 integer j
 real(8) x(3),pot,zi(3,3),zj(3,3),trace,emad

 call madelung2(emad)   ! sum'_R 1/R
! eps=trace(3,epsil)/3d0

    e_coul=0
    zi=atoms(i)%charge
    do j=1,n
       zj=atoms(j)%charge
       if(i.eq.j) then
          e_coul=e_coul+emad*trace(matmul(zi,transpose(zi)))/3d0
       else
!          x=atoms(i)%equilibrium_pos+dsp(:,i)-atoms(j)%equilibrium_pos-dsp(:,j)
          x=dsp(:,i)-dsp(:,j)
          call ewapot(x,pot)   ! sum_R 1/|R+x|
          e_coul=e_coul+pot*trace(matmul(zi,transpose(zj)))/3d0
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
 real(8), intent(out) :: emadelung
 real(8) x(3),ene !det

 x=0; x(1)=1d-6
 call ewapot(x,ene)
 emadelung=ene-1/sqrt( dot_product(x,matmul(epsinv,x))*det(epsil) )

 end subroutine madelung2
!===========================================================
 subroutine test_ewald
! use constants
! use lattice
! use geometry
! use params
! use atoms_force_constants
! use born
 use ios , only : ulog
 implicit none
 integer tau1,tau2,al
 real(8) x(3),dpot(3),r12(3),potp,potm,force1,pot0

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
 use atoms_force_constants
 use constants, only : pi
 use fourier
 use lattice, only : rs1,rs2,rs3,g01,g02,g03 , volume_r
 use ios , only : ulog
 use params, only : verbose
 use born , only : epsil
 implicit none
 integer, intent(in):: born_flag,ncfg
 real(8), intent(inout) :: dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg) ! dsp(:,:,:),frc(:,:,:)
 real(8), allocatable :: ewald_force(:,:),frcr(:,:),frcg(:,:),disr(:,:),disg(:,:)
 real(8) rr(3),dyn_coul(3,3),rcutoff_ewa,eta,gcutoff_ewa,deteps3
 integer nsc,n3(3),tau,taup,al,icfg,j

 write(*,*)'SUBTRACT_COULOMB_FORCE: size of force is 3,',size(frc(1,:,1)),size(frc(1,1,:))
 
! Beginning of ewald part for deduction of ewald force in real space ------------------
  if(born_flag.eq.0) then      ! skip if born_flag=0 
     return
  elseif(born_flag.eq.1) then  !  use ewald sums if=1, -FT(Dyn(K)*disp(k)) if=2 

     write(ulog,*)'MAIN: going to call Ewaldforce '
! set cutoffs for Ewald sums ! this should come after the supercell is read!!!!
     allocate(ewald_force(3,natom_super_cell))
     deteps3=det(epsil)**0.3333
     eta=sqrt(pi*deteps3)/(volume_r**0.3333)
     rcutoff_ewa=10*sqrt(deteps3)/eta
     gcutoff_ewa=10*eta*2/sqrt(deteps3)
! generate real space and reciprocal space translation vectors of the supercell for Ewald sums
     call make_grid_shell_ewald(rs1,rs2,rs3,rcutoff_ewa,gcutoff_ewa,epsil)
! output is r_ewald(3,nr_ewald) and g_ewald(3,ng_ewald) 

     do icfg=1,ncfg

! here displ is the absolute cartesian coordinates of displaced atoms
        call ewaldforce(natom_super_cell,dsp(:,:,icfg),ewald_force)

        if ( verbose ) then
           write(ulog,*)' Ewald force for configuration #',icfg
           do j=1,natom_super_cell
              write(ulog,4)j, ewald_force(:,j)
           enddo
        endif

        frc(:,:,icfg)=frc(:,:,icfg) - ewald_force(:,:)

     enddo

     deallocate(ewald_force)

  elseif(born_flag.eq.2) then  !  use dynmat_ewald*u_k (in reciprocal space)  

! fourier transform forces and displacements in the supercell, subtract the q=0 ewald term
! and fourier transform back before svd

!    allocate(map_rtau_sc(natom_prim_cell,nrgrid ))  ! needed for extract_fourier

! finds the mapping supercell atom<-> (tau,primcell translations defined in rgrid)
!    call find_map_rtau2sc(nrgrid,rgrid,map_rtau_sc)

!  erfc(4)=1.5d-8  so if eps=4/vol**.333 the real space terms become negligible
! exp(-18)=1.5d-8  so G larger than 5-6 shells will not contribute even if eps=4/vol^.3 

     allocate(frcr(3,nrgrid),frcg(3,nggrid))
     allocate(disr(3,nrgrid),disg(3,nggrid))
     do icfg=1,ncfg
        do tau=1,natom_prim_cell
       
           do j=1,nrgrid
              rr(1)=(rgrid(:,j) .dot. g01)/2/pi
              rr(2)=(rgrid(:,j) .dot. g02)/2/pi
              rr(3)=(rgrid(:,j) .dot. g03)/2/pi
              n3=nint(rr) 
              call findatom_sc(n3,tau,nsc) 
              frcr(:,j)=frc(:,nsc,icfg)
              disr(:,j)=dsp(:,nsc,icfg)
           enddo
           do al=1,3
              call fourier_r2k(frcr(al,:),frcg(al,:))
              call fourier_r2k(disr(al,:),disg(al,:))
           enddo

! now subtract ewald force in the form -D(k)u(k); requires u(k)
           do taup=1,natom_prim_cell
              if(taup.eq.tau) cycle  ! tau and tau' must be different
        
              do j=1,nggrid
! Long-range Coulomb dynamical matrix for G-vector labeled by j
                 call dyn_coulomb(tau,taup,j,dyn_coul)
! subtract coulomb
                 frcg(:,j)=frcg(:,j)+matmul(dyn_coul,disg(:,j))  ! + because -grad(force)=phi
              enddo
           enddo
! fourier transform back
           do al=1,3
              call fourier_k2r(frcg(al,:),frcr(al,:))
              call fourier_k2r(disg(al,:),disr(al,:))
           enddo
           do j=1,nrgrid
              rr(1)=(rgrid(:,j) .dot. g01)/2/pi
              rr(2)=(rgrid(:,j) .dot. g02)/2/pi
              rr(3)=(rgrid(:,j) .dot. g03)/2/pi
              n3=nint(rr) 
              call findatom_sc(n3,tau,nsc) 
! write into the original array 
              frc(:,nsc,icfg)=frc(:,nsc,icfg)-frcr(:,j)
           enddo

        enddo
     enddo
     deallocate(frcr,frcg,disr,disg)

  endif

4 format(i4,9(1x,g11.4))

 end subroutine subtract_coulomb_force
!--------------------------------------------------
 
 end module ewald 

!==============================================
