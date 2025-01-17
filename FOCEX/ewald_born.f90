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
 integer i1,i2,i3,mxshl,cnt,maxx,m
 real(r15) v(3),rr,aux2(3,3),epsinv(3,3)
 real(r15), allocatable :: aux(:,:),lengths(:)
 integer, allocatable :: msort(:)
 type(vector) y1,y2,y3

! call apply_metric(x1,x2,x3,epsilo,t1,t2,t3)
 aux2=epsilo
 call inverse_real(aux2,epsinv,3)

 call get_upper_bounds(x1,x2,x3,rcut,maxx,mxshl)
! call calculate_volume(x1,x2,x3,om0)
! max=nint(15d0/3d0*rcut**3/om0)+10
! m=nint((max/4.2)**0.333)+4

 allocate(aux(3,maxx),lengths(maxx),msort(maxx))
 lengths =1d20; aux=1d20
 write(ulog,4)'MAKE_GRID_SHELL, maxx(Rgrid), m,for rcut=',maxx,float(mxshl),rcut
 write(*   ,4)'MAKE_GRID_SHELL, maxx(Rgrid), m for rcut=',maxx,float(mxshl),rcut
 cnt=1
 shelloop: do m=0,mxshl
 do i1=-m,m
 do i2=-m,m
 do i3=-m,m
    if(iabs(i1).ne.m.and.iabs(i2).ne.m.and.iabs(i3).ne.m)cycle

! generate vectors in a grid
    v= v2a(i1*x1 + i2*x2 + i3*x3)
    rr=sqrt(dot_product(v,matmul(epsinv,v)))
    if(rr.gt.rcut) cycle

    aux(:,cnt)=v
    lengths(cnt)=rr
    cnt=cnt+1
    if (cnt.gt.maxx) then
       write(ulog,*)'RLOOP: maxx size exceeded, need to increase variable maxx from ',maxx
       exit shelloop
    endif

 enddo
 enddo
 enddo
 enddo shelloop

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
 call get_upper_bounds(y1,y2,y3,gcut,maxx,mxshl)
 write(ulog,4)'MAKE_GRID_SHELL, maxx(g_ewald), for gcut=',maxx,gcut
 write(*,4)'MAKE_GRID_SHELL, maxx(g_ewald), mxshl for gcut=',maxx,float(mxshl),gcut
 allocate(aux(3,maxx),lengths(maxx),msort(maxx))
 lengths =1d20; aux=1d20

 cnt=1
 gshelloop: do m=0,mxshl
 do i1=-m,m
 do i2=-m,m
 do i3=-m,m
    if(iabs(i1).ne.m.and.iabs(i2).ne.m.and.iabs(i3).ne.m)cycle

! generate vectors in a grid
    v= v2a(i1*y1 + i2*y2 + i3*y3)
    rr=sqrt(dot_product(v,matmul(epsilo,v)))
    if(rr.gt.gcut) cycle

    aux(:,cnt)=v
    lengths(cnt)=rr
    cnt=cnt+1
    if (cnt.gt.maxx .and. maxx.ne.1) then
       write(ulog,*)'GLOOP: maxx size exceeded, need to increase variable maxx from ',maxx
       exit gshelloop
!      stop
    endif

 enddo
 enddo
 enddo
 enddo gshelloop

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
 subroutine dewapot_g(x,dpot)
!! calculates the -d/dx(sum_R 1/|R+x| -background) = force  with the corrected metric
!! Assumes translation vectors are r_ewald(:,igrid),g_ewald(:,igrid)
! use constants
 use lattice
 use params
 use ios, only: ulog,write_out
 use born , only : epsil,epsinv
 implicit none
 real(r15), intent(out) :: dpot(3)
 real(r15), intent(in) :: x(3)
 integer igrid
 real(r15) termg,qpg(3),dd,my_erfc,geg,ep,vd(3)

    dpot=0
    do igrid=2,ng_ewald
       qpg=g_ewald(:,igrid)
       geg=dot_product(qpg,matmul(epsil,qpg))
       if (geg.gt.100*eta*eta) exit
       termg=exp(-geg/4/eta/eta)/geg
       dpot=dpot + termg * qpg*sin(qpg .dot. x) *4*pi/volume_r  ! this is supercell volume
    enddo
    if( geg.lt. 36*eta*eta ) then
       call warn3(ulog,'DEWAPOT_G sum not converged ',termg)
    endif

! now add the G=0 term
    dpot=dpot + 4*pi/volume_r*matmul(epsinv,x)/3d0

! write(*,3)'dpot=',dpot
3 format(9(1x,g11.4))

 end subroutine dewapot_g
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
 subroutine ewaldforce(n,pos,fewald)
!! calculates the Coulomb force on atoms in supercell (in 1/Ang)
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
 real(r15) dpot(3),r12(3),z1(3,3),z2(3,3) 


 fewald=0
 do tau1=1,n  ! caculate the total force on tau1
    z1=atom_sc(tau1)%charge
 do tau2=1,n
    z2=atom_sc(tau2)%charge

    if (tau1.eq.tau2) cycle ! images of tau1 have zero force on tau1 because of inversion symmetry
    r12=pos(:,tau1)-pos(:,tau2) ! pos includes equilibrium positions
! force on tau1 from tau2 and all its images inluding epsilon; beware of the sign!
    call dewapot(r12,dpot)  ! has epsil-modified metric included
    fewald(:,tau1)=fewald(:,tau1) + matmul(z1,matmul(z2,dpot)) * ee/(4*pi*eps0)*1d10  ! to convert to eV/Ang

 enddo
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
 subroutine ewaldforce_G(n,pos,fewald)
!! calculates the Coulomb force on atoms in supercell (in 1/Ang)
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
 real(r15) dpot(3),r12(3),z1(3,3),z2(3,3) 


 fewald=0
 do tau1=1,n  ! caculate the total force on tau1
    z1=atom_sc(tau1)%charge
 do tau2=1,n
    z2=atom_sc(tau2)%charge

    if (tau1.eq.tau2) cycle ! images of tau1 have zero force on tau1 because of inversion symmetry
    r12=pos(:,tau1)-pos(:,tau2) ! pos includes equilibrium positions
! force on tau1 from tau2 and all its images inluding epsilon; beware of the sign!
    call dewapot_g(r12,dpot)  ! has epsil-modified metric included
    fewald(:,tau1)=fewald(:,tau1) + matmul(z1,matmul(z2,dpot)) * ee/(4*pi*eps0)*1d10  ! to convert to eV/Ang

 enddo
 enddo

3 format(a,2i5,9(1x,g11.4))

 end subroutine ewaldforce_G
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
!! subtracts the non-analytical long-range contribution of Coulomb forces from total
!! force to get and fit the short-range contribution. Then this NA part is added back
!! to the fitted short-range dynamical matrix in subroutine dyn_coulomb
! dsp is the cartesian position not including the equilibrium positions
 use atoms_force_constants
! use constants, only : pi
 use fourier
 use lattice, only : rs1,rs2,rs3,g01,g02,g03,rws26,gws26, volume_r,volume_r0 !,cart2red_r
 use ios, only : ulog ,write_out
 use params, only : verbose, tolerance
 use born , only : epsil,dyn_naq0,dyn_na
 implicit none
 integer, intent(in):: born_flag,ncfg
 real(r15), intent(inout) :: dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg)
! automatic arrays here (lost after the call, memory released)
 integer nsc(natom_prim_cell,nrgrid)
 complex(r15) frcr(3,natom_prim_cell,nrgrid),frcg(3,natom_prim_cell,nggrid)
 complex(r15) disr(3,natom_prim_cell,nrgrid),disg(3,natom_prim_cell,nggrid)
 complex(r15) fcoulg(3,natom_prim_cell,nggrid),fcoulr(3,natom_prim_cell,nrgrid)
 complex(r15) auxr(nrgrid),auxg(nggrid),aravg(nrgrid) ,agavg(nggrid) 
 real(r15) ewald_force(3,natom_super_cell),ewald_force0(3,natom_super_cell)
! real(r15), allocatable :: frcr(:,:,:),frcg(:,:,:),disr(:,:,:),disg(:,:,:)
 complex(r15) ddn(3,3,3,nggrid)
 real(r15) rr(3),deteps3,asr(3,3),pos(3,natom_super_cell) 
 integer n3(3),tau,taup,al,be,icfg,i,j,l,jj,nsc1
 logical isperiodic

 if(born_flag.lt.1) return      ! skip if born_flag=<1 , otherwise, the NA term is always added

! allocate( dyn_na(natom_prim_cell,natom_prim_cell,3,3,nggrid) )

 write(*,*)'SUBTRACT_COULOMB_FORCE: size of force is 3,',size(frc(1,:,1)),size(frc(1,1,:))

! ------------------------------  SOME PRELIMINARY CALCULATIONS  ----------------------------------
! erfc(4)=1.5d-8  so if eps=4/vol**.333 the real space terms become negligible
! exp(-18)=1.5d-8  so G larger than 5-6 shells will not contribute even if eps=4/vol^.3
! set cutoffs for Ewald sums ! this should come after the supercell is read!!!!
     deteps3=det(epsil)**0.3333
     eta=sqrt(pi*deteps3)/(volume_r**0.3333) ! so both R-sums and G-sums converge at the same rate
     rcutoff_ewa=8*sqrt(deteps3)/eta
     gcutoff_ewa=8/sqrt(deteps3)*eta*2
     write(ulog,5)'SUBTRACT: deteps^1/3, eta=',n3,deteps3,eta
     write(ulog,5)'SUBTRACT: rcutoff,gcutoff=',n3,rcutoff_ewa,gcutoff_ewa

! generate real space and reciprocal space translation vectors of the supercell for Ewald sums
! output is r_ewald(3,nr_ewald) and g_ewald(3,ng_ewald)
     call make_grid_shell_ewald(rs1,rs2,rs3,rcutoff_ewa,gcutoff_ewa,epsil)

     call nsc_from_rgrid(nrgrid,rgrid,nsc) ! mapping from (rgrid,tau) to supercell: nsc(tau,jgrid)

! calculation of DYN_NA(tau,taup,G) (no mass denominator) for multiplication by disp(taup,G) to compare with -force_ewald(tau,G)
     if (born_flag.eq.4) then ! calculate non-analytical part phi_NA in real space and subtract from forces

       write(ulog,*)'NA term in real space on the Rgrid: tau,al;taup,be,phi_NA(Rgrid)='
       if(allocated(dyn_naq0)) then
           deallocate(dyn_naq0)
       endif
       allocate( dyn_naq0(natom_prim_cell,3,3) )
!      allocate( dyn_na0(natom_prim_cell,natom_prim_cell,3,3,nggrid) )
!      call phi_na_0
!      dyn_na=0
!      do tau =1,natom_prim_cell
!      do taup=1,natom_prim_cell
!      do j=1,nggrid
!         dyn_na(tau,taup,:,:,j)=dyn_na(tau,taup,:,:,j)+dyn_naq0(tau,:,:)* &
! & gws_weights(j)*exp(ci*dot_product(ggrid(:,j),(atompos(:,taup)-atompos(:,tau))))
!      enddo
!      enddo
!      enddo
       
     else  ! calculate non-analytical part DYN_NA(q); no mass factor

       do tau=1,natom_prim_cell
       do taup=1,natom_prim_cell
       do j=1,nggrid
          call dyn_coulomb(tau,taup,ggrid(:,j),dyn_na(tau,taup,:,:,j),ddn(:,:,:,j))
       enddo
       enddo
       enddo

     endif
       
     do j=1,nggrid
     do tau =1,natom_prim_cell
     do taup=1,natom_prim_cell
        write(ulog,*)'BFLAG,j, tau,taup= ',born_flag,j,tau,taup
        call write_out(ulog,'SUBTRACT:  dyn_NA(tau,taup) ',(dyn_na(tau,taup,:,:,j)))
     enddo
     enddo
     enddo

! ----------------   NOW DO SUBTRACTION FOR EVERY CONFIGURATION    --------------------
     config_sum: do icfg=1,ncfg

! for each config convert disp to a 3 x N0 x Nrgrid format needed for FT
        do j=1,nrgrid
        do tau=1,natom_prim_cell
           nsc1=nsc(tau,j)
           disr(:,tau,j)=dsp(:,nsc1,icfg)  ! this should be periodic if j on boundary
        enddo
        enddo
! Fourier transform displacement to convert to 3 x N0 x Nggrid format 
        do tau=1,natom_prim_cell
        do al=1,3
           auxr=disr(al,tau,:)           !           auxr=frcr(al,tau,:)
! is auxr, defined on rgrid, periodic of period rws26?
           call check_periodic(nrgrid,rgrid,rws26,auxr,aravg,isperiodic)
           if(isperiodic) then
              call fourier_r2k(auxr,auxg)   !           call fourier_r2k(auxr,auxg)
           else
              write(ulog,*)'icfg=',icfg,' disr was not periodic!!!'
              call fourier_r2k(aravg,auxg)   !           call fourier_r2k(auxr,auxg)
           endif 
           disg(al,tau,:)=auxg           !           frcg(al,tau,:)=auxg
        enddo
        call write_out(ulog,' disg(tau) ',transpose(disg(:,tau,:)))
        enddo



        if    (born_flag.eq.1 .or. born_flag.eq.5) then    ! - - - - - - - - - - - - - - - 

            cycle    ! no subtraction
 
        elseif(born_flag.eq.2) then    ! - - - - - - - - - - - - - - -  
! Ewald force_g is the long-range part which is subtracted in real space 

            write(ulog,*)'MAIN: going to call Ewaldforce_G '

            do i=1,natom_super_cell
               pos(:,i)=v2a(atom_sc(i)%equilibrium_pos)
            enddo 
            call ewaldforce_g(natom_super_cell,pos,ewald_force0)

! ewaldforce_g requires the absolute cartesian coordinates of atoms, so pos=eq+displacement
            do i=1,natom_super_cell
               pos(:,i)=v2a(atom_sc(i)%equilibrium_pos) + dsp(:,i,icfg)
            enddo 
            call ewaldforce_g(natom_super_cell,pos,ewald_force)

!            if ( verbose ) then
               write(ulog,*)' Ewald force F,F0,F-F0 for configuration #',icfg
               do j=1,natom_super_cell
                  write(ulog,4)j, ewald_force(:,j),ewald_force0(:,j),ewald_force(:,j)-ewald_force0(:,j)
               enddo
!            endif

            frc(:,:,icfg)=frc(:,:,icfg) - ( ewald_force(:,:)- ewald_force0(:,:) )

! compare to -D_NA(t,t',G)*u(t',g)
            fcoulg=0
            do tau=1,natom_prim_cell
            do j=1,nggrid
               do taup=1,natom_prim_cell
                 fcoulg(:,tau,j)=fcoulg(:,tau,j)-real(matmul(dyn_na(tau,taup,:,:,j),disg(:,taup,j))) 
               enddo
            enddo
            call write_out(6,' fcoulg(1) ',transpose(fcoulg(:,tau,:)))
            enddo
! fourier transform back to real space and get fcoulr
            do tau=1,natom_prim_cell
            do al=1,3
               auxg=fcoulg(al,tau,:)          !       auxg=disg(al,tau,:)
               call check_periodic(nggrid,ggrid,gws26,auxg,agavg,isperiodic)
               if(isperiodic) then
                  call fourier_k2r(auxg,auxr)   !           call fourier_r2k(auxr,auxg)
               else
                  write(ulog,*)'icfg=',icfg,' fcoulg was not periodic!!!'
                  call fourier_k2r(agavg,auxr)   !           call fourier_r2k(auxr,auxg)
               endif 
               fcoulr(al,tau,:)=auxr      
            enddo
            enddo

! subtract fcoulr from forces; first from grid indices go back to supercell indices
            do j=1,nrgrid
            do tau=1,natom_prim_cell
               nsc1=nsc(tau,j)
!              frc(:,nsc1,icfg)=frc(:,nsc1,icfg)-fcoulr(:,tau,j)
               write(ulog,5)'tau,j,nsc1,fcoulr,df_ewald=',tau,j,nsc1, fcoulr(:,tau,j),ewald_force(:,nsc1)-ewald_force0(:,nsc1)
            enddo
            enddo



        elseif(born_flag.eq.3) then    ! - - - - - - - - - - - - - - -  
! NA subtraction in reciprocal space F=F - (-FFT(Dyn(K)*disp(k)) )

! bflag=3 takes only Zq.Zq/qEq while bflag=4 takes the sum_G (K=q+G) Zk.ZK exp(-KEK/4e^2) /KEK
! i.e. dyn_NA_ewald either with only G sum (bflag=4) or full R and G sums (bflag=5)

! now calculate non-analytical force, fcoulg, in the form -D(k)*disg(k)
           fcoulg=0
           do tau=1,natom_prim_cell
              asr=0
              do j=1,nggrid
              do taup=1,natom_prim_cell
! Long-range Coulomb dynamical matrix for G-vector labeled by j
                 if (length(ggrid(:,j)).lt.tolerance) then  ! ASR is only for G=0
                    asr=asr+dyn_na(tau,taup,:,:,j)
                 endif
                 fcoulg(:,tau,j)=fcoulg(:,tau,j)-real(matmul(dyn_na(tau,taup,:,:,j),disg(:,taup,j))) 
              enddo
              enddo
              if(maxval(abs(asr)).gt. tolerance) then
                 call warn(6)
                 call write_out(6,'bflag=3: dyn_na asr should be zero ',asr)
              endif
           enddo

! fourier transform back to real space and get fcoulr
           do tau=1,natom_prim_cell
           do al=1,3
              auxg=fcoulg(al,tau,:)          !       auxg=disg(al,tau,:)
              call check_periodic(nggrid,ggrid,gws26,auxg,agavg,isperiodic)
              if(isperiodic) then
                 call fourier_k2r(auxg,auxr)   !           call fourier_r2k(auxr,auxg)
              else
                 call fourier_k2r(agavg,auxr)   !           call fourier_r2k(auxr,auxg)
              endif 
              fcoulr(al,tau,:)=auxr          !       disr(al,tau,:)=auxr
           enddo
           enddo

! subtract fcoulr from forces; first from grid indices go back to supercell indices
           do j=1,nrgrid
           do tau=1,natom_prim_cell
              nsc1=nsc(tau,j)
              frc(:,nsc1,icfg)=frc(:,nsc1,icfg)-fcoulr(:,tau,j)
           enddo
           enddo


        elseif(born_flag.eq.4) then    ! - - - - - - - - - - - - - - - ! NA subtraction in real space F(tau)=F(tau) - (-phi_NA(tau,R+taup)*disp(taup)) )

! convert dsp in supercell to disr on the WS inner grid
           do j=1,nrgrid
           do tau=1,natom_prim_cell
              nsc1=nsc(tau,j)
              disr(:,tau,j)=dsp(:,nsc1,icfg) 
           enddo
           enddo

! Non-analytical Coulomb force in real space on Rgrid
           fcoulr=0
           do j=1,nrgrid
           do tau=1,natom_prim_cell
              do l=1,nrgrid
       !         rr=rgrid(:,l)-rgrid(:,j)
! find the corresponding rgrid number
       !         call bring_to_ws_r(rr,rws26,foldedr)
! jj corresponds to rgrid(j)-rgrid(l) folded back into the supercell WS cell
       !         call find_in_array(3,foldedr,nrgrid,rgrid,jj,tolerance)  not needed since phi_na is constant over cells
       !         if (jj.eq.0 .or. jj .gt. nrgrid) then
       !            write(ulog,6)'SUBRTRACT: could not find array ',foldedr,reduce_r(foldedr)
       !            write(   *,6)'SUBRTRACT: could not find array ',foldedr,reduce_r(foldedr)
       !            stop
       !         endif
                 do taup=1,natom_prim_cell
   !                fcoulr(:,tau,j)=fcoulr(:,tau,j)- matmul(real(phi_na(tau,taup,:,:,jj)),disr(:,taup,l))
!                   fcoulr(:,tau,j)=fcoulr(:,tau,j)- matmul(phi_naq0(tau,taup,:,:),disr(:,taup,l)) *volume_r0/volume_r 
                 enddo
              enddo
           enddo
           enddo

! find correspondance between rgrid point and (tau,R)
           do j=1,nrgrid
           do tau=1,natom_prim_cell
              nsc1=nsc(tau,j)
              frc(:,nsc1,icfg)=frc(:,nsc1,icfg)-fcoulr(:,tau,j)
           enddo
           enddo

        endif

     enddo config_sum


4 format(i4,99(1x,f9.4))
5 format(a,3(i4),33(1x,f10.5))
6 format(a,9(1x,f10.5))
8 format(a,2(i3,i2),999(1x,f9.4))

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
 function coulg(tau,taup,g) result(h)
 use born , only : epsinv,epsil
 use constants, only : r15
 use lattice, only : volume_r0
 use atoms_force_constants, only : natom_prim_cell,atom0
 use geometry, only : length
 use params, only : coef
 implicit none
 integer, intent(in) :: tau,taup 
 real(r15), intent(in) :: g(3)
 integer al,be
 real(r15) h(3,3),geg,zg1(3),zg2(3)
 
 if(length(g).lt.1d-6) then
    h = matmul(matmul(atom0(tau )%charge,epsinv),transpose(atom0(taup)%charge))*coef/3 
 else
    zg1 = matmul(atom0(tau )%charge,g)
    zg2 = matmul(atom0(taup)%charge,g)
    geg = dot_product(g,matmul(epsil,g))
    do al=1,3
    do be=1,3
       h(al,be)= zg1(al)*zg2(be)/geg *coef
    enddo
    enddo
 endif

 end function coulg
!--------------------------------------------------
 subroutine ewald_2nd_deriv_hat(q,tau,taup,nr,ng,rgrid,ggrid,etaew,d2ew,d3ew)
!! calculates the second derivative of the Ewald potential to be used in the Non-analytical correction; step(non-smooth) phase
!! input is the two atoms tau and taup, output is the 3x3 block D^EW_tau,taup(q) which goes to D^NA for q \to 0
 use params, only : tolerance
 use lattice
 use geometry, only : length
 use born , only : epsil,epsinv
 use atoms_force_constants, only : atom0,atompos
 use ios , only : ulog,write_out
 implicit none
 integer, intent(in) :: tau,taup,nr,ng
 real(r15), intent(in) :: q(3),rgrid(3,nr),ggrid(3,ng),etaew
 complex(r15), intent(out) :: d2ew(3,3),d3ew(3,3,3)
 integer igrid,al,be,ga
 real(r15) termg,sd(3),del(3),dd,my_erfc,geg,ep,vd(3),dta(3),qpg(3),gscale,qal(3),qbe(3)
 complex(r15) dhat(3,3), hout(3,3)

 dta=atompos(:,taup)-atompos(:,tau)  ! choice for the smooth phase
 gscale=6.28/volume_r0**0.3333

! real space term  
    dhat=0 ; d3ew=0
    do igrid=1,nr
       sd=etaew*(dta+rgrid(:,igrid))
       del=matmul(epsinv,sd)
       dd=sqrt(dot_product(sd,del))
       if(dd.lt.1d-6) cycle ! exclude self-interactions
       if (dd.gt.6) cycle
       hout= hfunc(del,dd) * exp(ci*dot_product(q,rgrid(:,igrid))) 
!      write(6,4)'rgrid,dd,hout=',igrid,dd,hout(1,1)
!      dhat=dhat - hout * cos(dot_product(q,rgrid(:,igrid))) 
       dhat=dhat - hout
       do ga=1,3
!         d3ew(:,:,ga)=d3ew(:,:,ga) - matmul(matmul(atom0(tau )%charge,hout),transpose(atom0(taup)%charge))  &
!  &                                * sin(dot_product(q,rgrid(:,igrid))) *rgrid(ga,igrid) 
!         d3ew(:,:,ga)=d3ew(:,:,ga) - hout * sin(dot_product(q,rgrid(:,igrid))) *rgrid(ga,igrid) 
          d3ew(:,:,ga)=d3ew(:,:,ga) - ci*rgrid(ga,igrid) *hout  
       enddo
!      write(*,4)'R_sum: igrid,d,dhat_11=',igrid,dd, dhat(1,1)
    enddo
     if(tau.eq.taup) dhat=dhat - 4/3d0/sqrt(pi)*epsinv 
!    if(tau.eq.taup) d3ew=d3ew - 4/3d0/sqrt(pi)*epsinv   !!! TO BE CHECKED !!! not needed since does not depend on q
    dhat=dhat* etaew*etaew*etaew/sqrt(det(epsil)) 
    d3ew=d3ew* etaew*etaew*etaew/sqrt(det(epsil)) 
!   write(6,*)'Last R-term(1,1)=',hout(1,1)*etaew*etaew*etaew/sqrt(det(epsil)) 
!   call write_out(6,'total of Rsum terms ',dhat)

! reciprocal space term  
    do igrid=1,ng
       qpg=ggrid(:,igrid)+q
       if(length(qpg).lt.1d-5*gscale) cycle  ! exclude G+q=0
       geg=dot_product(qpg,matmul(epsil,qpg))
       if (geg.gt.100*etaew*etaew) exit  ! assumes ggrid is sorted 
       termg=exp(-geg/4/etaew/etaew)/geg
!      write(6,4)'ggrid,gg,term=',igrid,sqrt(geg),termg
!      qal(:)=matmul(atom0(tau )%charge,qpg)
!      qbe(:)=matmul(atom0(taup)%charge,qpg)
       do al=1,3
       do be=1,3
!         hout(al,be)= termg * qal(al)*qbe(be)*cos(dot_product(qpg,dta)) 
!         hout(al,be)= termg * qpg(al)*qpg(be)*cos(dot_product(qpg,dta)) 
          hout(al,be)= termg * qpg(al)*qpg(be)*exp(-ci*dot_product(qpg,dta)) 
!! NEED TO ADD THIRD DERIVATIVE !!
       enddo
       enddo
       dhat=dhat + hout*4*pi/volume_r0 
!      if(igrid.eq.1) call write_out(6,'G=0 term of ewald ',hout*180.9557368)
    enddo
!   write(*,3)'final EW2DERIV:dhat=',dhat
!   write(6,*)'Last G-term(1,1)=',hout(1,1)*4*pi/volume_r0 

    dhat=dhat /(4*pi*eps0)* ee*1d10
!   call write_out(6,'final EW2DERIV:dhat scaled by 1d10*ee/eps0 ',dhat)

    d2ew= matmul(matmul(atom0(tau)%charge,dhat),transpose(atom0(taup)%charge))
 !  d2ew= d2ew*exp(-ci*(q.dot.dta))   ! smooth phase convention
     
! Fix d3ew here ; this is copied from old NA PArlinski term
 !      rtp  = v2a(atom0(taup)%equilibrium_pos-atom0(tau)%equilibrium_pos)
 !      phase= exp(ci*(q.dot.rtp))
 !      mysf = sf * phase
 !      mydsf=(dsf + ci*rtp*sf) * phase
 !         do ga=1,3

 !            ddyn(al,be,ga) = dyn(al,be) * ( &
 !    &         atom0(tau )%charge(al,ga)/dot_product(atom0(tau )%charge(al,:),q) + &
 !    &         atom0(taup)%charge(be,ga)/dot_product(atom0(taup)%charge(be,:),q) - &
 !    &         dqeq(ga)/qeq  + ci*rtp(ga) + mydsf(ga)/mysf ) 

 !         enddo








3 format(a,99(1x,g14.7))
4 format(a,i4,99(1x,g11.4))

 end subroutine ewald_2nd_deriv_hat
!--------------------------------------------------
 subroutine ewald_2nd_deriv_Gnew(q,tau,taup,ng,ggrid,etaew,d2ew,d3ew)
!! calculates the dynamical matrix associated with the NA term (only G sum) 
!! input is the two atoms tau and taup, output is the 3x3 block D^EW_tau,taup(q) 
!! which goes to D^NA for q \to 0 ; G are reciprocal lattice vectors (not supercell)
!! uses the new phase convention
 use params, only : tolerance,coef
 use lattice
 use fourier, only : nrgrid,rgrid,rws_weights
 use geometry, only : length
 use born , only : epsil,epsinv
 use atoms_force_constants, only : atom0,atompos
 use ios , only : ulog,write_out
 implicit none
 integer, intent(in) :: tau,taup,ng
 real(r15), intent(in) :: q(3),ggrid(3,ng),etaew
 complex(r15), intent(out) :: d2ew(3,3),d3ew(3,3,3)
 complex(r15) phase
 integer igrid,al,be,ga
 real(r15) termg,geg,ep,dta(3),qpg(3),gscale,hout(3,3),qal(3),qbe(3)

 dta=atompos(:,taup)-atompos(:,tau)
 gscale=6.28/volume_r0**0.3333

    d2ew=0 ; d3ew=0
! reciprocal space term  
    gloop: do igrid=1,ng
       qpg=ggrid(:,igrid)+q
       phase = 1 ! exp(-ci*dot_product(ggrid(:,igrid),dta))
       if(length(qpg).lt.1d-6*gscale) then 
          cycle
!          d2ew = d2ew+ 1/3d0*matmul(matmul(atom0(tau )%charge,epsinv),  &
! &               transpose(atom0(taup)%charge)) * phase 
       else
          geg=dot_product(qpg,matmul(epsil,qpg))
          if (geg.gt.90 *etaew*etaew) exit gloop  ! assumes ggrid is sorted 
          termg=exp(-geg/4/etaew/etaew)/geg                
          qal(:)=matmul(atom0(tau )%charge,qpg)
          qbe(:)=matmul(atom0(taup)%charge,qpg)
 !        write(6,4)'ggrid,gg,term=',igrid,sqrt(geg),termg
          do al=1,3
          do be=1,3
             hout(al,be)= termg * qal(al)*qbe(be)*exp(-ci*dot_product(qpg,dta)) 
!! NEED TO ADD THE q-DERIVATIVE for group velocities !!
          enddo
          enddo
          d2ew=d2ew + hout * phase
       endif
!      write(ulog,4)'G,d_Ew(g)=',igrid,real(d2ew)
    enddo gloop

    d2ew=d2ew * coef 
    d3ew=d3ew * coef 
3 format(a,99(1x,g14.7))
4 format(a,i4,99(1x,g11.4))

 end subroutine ewald_2nd_deriv_Gnew
!--------------------------------------------------
 subroutine ewald_2nd_deriv_G(q,tau,taup,ng,ggrid,etaew,d2ew,d3ew)
!! calculates the dynamical matrix associated with the NA term (only G sum) 
!! input is the two atoms tau and taup, output is the 3x3 block D^EW_tau,taup(q) which goes to D^NA for q \to 0
! use constants, only : ci
 use params, only : tolerance,coef
 use lattice
 use fourier, only : nrgrid,rgrid,rws_weights
 use geometry, only : length
 use born , only : epsil,epsinv
 use atoms_force_constants, only : atom0,atompos
 use ios , only : ulog,write_out
 implicit none
 integer, intent(in) :: tau,taup,ng
 real(r15), intent(in) :: q(3),ggrid(3,ng),etaew
 complex(r15), intent(out) :: d2ew(3,3),d3ew(3,3,3)
 complex(r15)mysf,mydsf(3)
 integer igrid,al,be,ga
 real(r15) termg,geg,ep,dta(3),qpg(3),gscale,hout(3,3),qal(3),qbe(3)

 dta=atompos(:,taup)-atompos(:,tau)
 gscale=6.28/volume_r0**0.3333

    d2ew=0 ; d3ew=0
! reciprocal space term  
    gloop: do igrid=1,ng
!      qpg=ggrid(:,igrid)+q
       qpg=ggrid(:,igrid)
       if(length(qpg).lt.1d-10*gscale) then ! cycle  ! exclude G+q=0
         d2ew = d2ew+ 1/3d0*matmul(matmul(atom0(tau )%charge,epsinv),  &
&               transpose(atom0(taup)%charge)) 
       else
         qal(:)=matmul(atom0(tau )%charge,qpg)
         qbe(:)=matmul(atom0(taup)%charge,qpg)
         geg=dot_product(qpg,matmul(epsil,qpg))
         if (geg.gt.100*etaew*etaew) exit gloop  ! assumes ggrid is sorted 
         termg=exp(-geg/4/etaew/etaew)/geg
!        write(6,4)'ggrid,gg,term=',igrid,sqrt(geg),termg

         call structure_factor_complx(q+ggrid(:,igrid),nrgrid,rgrid,rws_weights,mysf,mydsf)

         do al=1,3
         do be=1,3
!           hout(al,be)= termg * qal(al)*qbe(be)*cos(dot_product(ggrid(:,igrid),dta)) 
            hout(al,be)= termg * qal(al)*qbe(be)
!! NEED TO ADD THIRD DERIVATIVE !!
         enddo
         enddo
!        dhat=dhat + hout
         d2ew=d2ew + hout*mysf * exp(ci*(ggrid(:,igrid).dot.dta))
       endif
    enddo gloop

!   d2ew=dhat * coef * exp(-ci*(q.dot.dta))
    d2ew=d2ew * coef * exp(ci*(q.dot.dta))

3 format(a,99(1x,g14.7))
4 format(a,i4,99(1x,g11.4))

 end subroutine ewald_2nd_deriv_G

!--------------------------------------------------
 subroutine phi_na_g(phi_na,nrg,rgrid,ngg,ggrid,weighg)
!! calculates the non-analytical Coulomb FCs in real space grid within the supercell from its Fourier transform
 use lattice, only : volume_r0
 use geometry, only : length
 use born , only : epsil
 use atoms_force_constants, only : atom0,atompos,natom_prim_cell,natom_super_cell
! use ios , only : ulog,write_out
 use params, only : coef 
 implicit none
 integer, intent(in) :: ngg,nrg
 integer  igrid,tau,taup,al,be,i
 real(r15), intent(in) :: ggrid(3,ngg),rgrid(3,nrg),weighg(ngg)
 real(r15), intent(out) :: phi_na(natom_prim_cell,natom_prim_cell,3,3,nrg)
 real(r15) geg,dta(3),gscale, gtau(3),gtaup(3)

 gscale=6.28/volume_r0**0.3333

 phi_na=0
! reciprocal space term  
 do igrid=1,ngg
    if(length(ggrid(:,igrid)).lt.1d-10*gscale) cycle  ! exclude G=0
    geg=dot_product(ggrid(:,igrid),matmul(epsil,ggrid(:,igrid)))
!   write(6,4)'ggrid,gg,term=',igrid,sqrt(geg),termg
    do tau =1,natom_prim_cell
    do taup=1,natom_prim_cell
       gtau (:)=matmul(atom0(tau )%charge,ggrid(:,igrid))
       gtaup(:)=matmul(atom0(taup)%charge,ggrid(:,igrid))
    do i=1,nrg
       dta=rgrid(:,i)+atompos(:,taup)-atompos(:,tau)
    do al=1,3
    do be=1,3
       phi_na(tau,taup,al,be,i)= phi_na(tau,taup,al,be,i)+ gtau(al)*gtaup(be)/geg *  &
    &                 weighg(igrid) ! * cos(dot_product(ggrid(:,igrid),dta))
    enddo
    enddo
    enddo
    enddo
    enddo
 enddo

 phi_na = phi_na * coef / dble(natom_super_cell)

3 format(a,99(1x,g14.7))
4 format(a,i4,99(1x,g11.4))

 end subroutine phi_na_g
!--------------------------------------------------
 subroutine phi_na_0(phi_naq0)
!! calculates the non-analytical Coulomb FCs in real space grid within the supercell from its Fourier transform
 use lattice, only : volume_r0
 use born , only : epsinv
 use atoms_force_constants, only : atom0,natom_prim_cell
 use params, only : coef
 implicit none
 real(r15), intent(out) :: phi_naq0(natom_prim_cell,natom_prim_cell,3,3)
 integer tau,taup

 do tau =1,natom_prim_cell
 do taup=1,natom_prim_cell
    phi_naq0(tau,taup,:,:) = matmul(matmul(atom0(tau )%charge,epsinv),transpose(atom0(taup)%charge))/3d0 
 enddo
 enddo

 phi_naq0 = phi_naq0 * coef   ! be careful should be supercell volume there!

 end subroutine phi_na_0
!--------------------------------------------------

 end module ewald

