!! module ewald: Ewald grids and sums for the non-analytical dynamical matrix (shared by FOCEX
!! and THERMACOND).  FOCEX's Ewald forces are in FOCEX/ewald_born.f90, module ewald_forces.
 module ewald
 use constants, only : r15,pi,ci
 use ios, only : ulog
 integer nr_ewald,ng_ewald
 real(r15), allocatable :: r_ewald(:,:),g_ewald(:,:)
 real(r15) rcutoff_ewa,gcutoff_ewa, eta
! E = 0 (no macroscopic field; 2026-11-02, for the elastic constants from the long-wave limit): when
! set, ewald_2nd_deriv_hat keeps of its G = 0 term only the analytic part (exp(-x)-1)(Z.q)(Z.q)/(q.eps.q),
! x = q.eps.q/(4 eta^2), i.e. the non-analytic (macroscopic-field) term of D(q -> 0) is left out
 logical :: ewald_no_macro = .false.

 contains

 subroutine cuts(eps,a0,bflg,et,rcut,gcut)
   real(r15), intent(in):: eps,a0
   integer, intent(in):: bflg
   real(r15), intent(out):: et,rcut,gcut
   et  =1
   rcut=0
   gcut=0
   if(bflg.lt.1) then
     return
   elseif( mod(bflg,10).eq.1 ) then ! equal terms in real and reciprocal for 1,2,3
     et=sqrt(pi*eps)/a0
     rcut=4*a0
     gcut=8/a0
   elseif(mod(bflg,10).le.3 .or. mod(bflg,10).eq.6) then ! equal terms in real and reciprocal for 2,3,6
     et=sqrt(2*pi*eps)/a0  ! this is to make the R&G sums equally convergent
     rcut=8*sqrt(eps)/et    *2.0  ! *2.0 is to generate more than needed!
     gcut=8/sqrt(eps)*et*2  *2.3
   else ! only sum in the reciprocal for 4 and above for _Gnew
     et=sqrt(pi*eps)/a0 *3.1 ! this is to use a larger cutoff for G vectors in Gnew and still converge
     rcut=8*sqrt(eps)/et    ! this mesh is not used for bflag>3 
     gcut=8/sqrt(eps)*et*2  ! need more G vectors to converge ewald sums
   endif
   write(ulog,7)'CUTS:bf,det,a0,et,rc,gc=',bflg,eps,a0,et,rcut,gcut

7 format(a,i5,99(g11.4))
 end subroutine cuts
!============================================================
 subroutine make_grid_shell_ewald(x1,x2,x3,rcut,gcut,epsilo)
!! generates a grid of vectors "grid" linear combinations of xi within a
!! cutoff length rcut. Sorted result is stored in r_ewald(3,nr_ewald)
!! also generates the same grid for the reciprocal space g_ewald(3,ng_ewald)
!! metric is epsilon in reciprocal space and epsinverse in real space
 use geometry
 use ios
 implicit none
 real(r15), intent(in ):: rcut,gcut,epsilo(3,3)
 type(vector), intent(in ):: x1,x2,x3
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
! write(ulog,4)'MAKE_GRID_SHELL, maxx(Rgrid), m,for rcut=',maxx,float(mxshl),rcut
 write(ulog,4)'MAKE_GRID_SHELL, maxx(Rgrid), m for rcut=',maxx,float(mxshl),rcut

 call run_grid(mxshl,maxx,x1,x2,x3,epsinv,aux,lengths,rcut,nr_ewald) 

 write(ulog,5) 'MAKE_GRID_SHELL: within rcut=',rcut,' r_ewald generated ',nr_ewald,' vectors'

 if (allocated(r_ewald)) deallocate(r_ewald)  ! overwrite if previously called
 allocate(r_ewald(3,nr_ewald))

 call sort(nr_ewald,lengths,msort,maxx)

  open(ujunk,file='ewald_vecs.dat')
  write(ujunk,*)'MAKE RGRID:----first 200 vectors--------------',nr_ewald
  do i1=1,nr_ewald
     r_ewald(:,i1)=aux(:,msort(i1))
     if (i1.lt.200) write(ujunk,3)r_ewald(:,i1),length(r_ewald(:,i1))
  enddo

 deallocate (aux, lengths, msort)

! now the G-vectors
 call make_reciprocal_lattice_2pi(x1,x2,x3,y1,y2,y3)
 call get_upper_bounds(y1,y2,y3,gcut,maxx,mxshl)
! write(ulog,4)'MAKE_GRID_SHELL, maxx(g_ewald), for gcut=',maxx,gcut
 write(*,4)'MAKE_GRID_SHELL, maxx(g_ewald), mxshl for gcut=',maxx,float(mxshl),gcut
 allocate(aux(3,maxx),lengths(maxx),msort(maxx))
 lengths =1d20; aux=1d20

 call run_grid(mxshl,maxx,y1,y2,y3,epsilo,aux,lengths,gcut,ng_ewald) 

 write(ulog,5) 'MAKE_GRID_SHELL: within gcut=',gcut,' g_ewald generated ',ng_ewald,' vectors'

 if (allocated(g_ewald)) deallocate(g_ewald)
 allocate(g_ewald(3,ng_ewald))

 call sort(ng_ewald,lengths,msort,maxx)

  write(ujunk,*)'MAKE GGRID:------------------------',ng_ewald
  do i1=1,ng_ewald
     g_ewald(:,i1)=aux(:,msort(i1))
     if (i1.lt.200) write(ujunk,3)g_ewald(:,i1),length(g_ewald(:,i1))
  enddo
!
  deallocate (aux, lengths, msort)
  close(ujunk)

3 format(9(1x,f10.5))
4 format(a,i8,9(1x,g14.7))
5 format(a,g11.4,a,i8,9(1x,g14.7))

 end subroutine make_grid_shell_ewald
!============================================================
 subroutine run_grid(mxshl,maxx,y1,y2,y3,epsilo,aux,lengths,gcut,ng_ewald) 
!! maxx and mxshl are initial guesses (upper bounds) for the size of the arrays and # of shells
!! given 3 translation vectors y1,y2,y3 and metric defined by epsilo, generates a 
!! grid, aux(3,maxx), of size ng_ewald within the sphere of radius gcut; 
!! lengths(maxx) is later used for sorting aux according to their lengths
 use geometry, only : v2a,vector
 use ios, only : ulog
 implicit none
 integer, intent(in):: mxshl,maxx
 integer, intent(out):: ng_ewald
 real(r15), intent(in ):: gcut,epsilo(3,3)
 real(r15), intent(out):: aux(3,maxx),lengths(maxx)
 type(vector), intent(in ):: y1,y2,y3
 integer i1,i2,i3,cnt,m
 real(r15) v(3),rr

 cnt=0; lengths=1d9; aux=1d8
 gshelloop: do m=0,mxshl
 do i1=-m,m
 do i2=-m,m
 do i3=-m,m
    if(iabs(i1).ne.m.and.iabs(i2).ne.m.and.iabs(i3).ne.m)cycle

! generate vectors in a grid
    v= i1*v2a(y1) + i2*v2a(y2) + i3*v2a(y3)
    rr=sqrt(dot_product(v,matmul(epsilo,v)))
    if(rr.gt.gcut) cycle
    cnt=cnt+1
    if (cnt.gt.maxx .and. maxx.ne.1) then
       write(ulog,*)'GLOOP: maxx size exceeded, need to increase variable maxx from ',maxx
       cnt=cnt-1
       exit gshelloop
    else
       aux(:,cnt)=v
       lengths(cnt)=rr
    endif

 enddo
 enddo
 enddo
 enddo gshelloop

 ng_ewald=cnt

 end subroutine run_grid
!============================================================
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
!! calculates the second derivative of the Ewald potential to be used in the Non-analytical correction; 
!! step(non-smooth) phase; no Born charge and no mass denominator
!! input: two atoms tau and taup, output is the 3x3 block D^EW_{tau,taup}(q) which goes to D^NA for q \to 0
 use constants, only : eps0scale,pi
 use geometry, only : v2a,det,trace
 use born , only : epsil,epsinv
 use atoms_force_constants, only : atom0
 implicit none
 integer, intent(in) :: tau,taup,nr,ng
 real(r15), intent(in) :: q(3),rgrid(3,nr),ggrid(3,ng),etaew
 complex(r15), intent(out) :: d2ew(3,3),d3ew(3,3,3)
 integer igrid,al,be,ga
 real(r15) sd(3),del(3),dd,my_erfc,geg,ep,vd(3),dta(3),qpg(3), &
&          qal(3),qbe(3),dum(3),ep3,tol,sumcut,dd2
 complex(r15) hout(3,3),gout(3,3,3),zz,gsum,rpgsum,consta,hmac(3,3),gmac(3,3,3)
 integer, external :: delta_k

 dta = v2a(atom0(taup)%equilibrium_pos)-v2a(atom0(tau)%equilibrium_pos)
 tol=1d-14; sumcut=sqrt(log(1d0/tol))
! write(ulog,*)'_hat : nr,Rcut=',nr,sumcut/etaew*sqrt(trace(epsil)/3)
! write(ulog,*)'_hat : ng,Gcut=',ng,2*sumcut*etaew/sqrt(trace(epsil)/3)
    d2ew=cmplx(0d0,0d0) !;  epsinv*coef/3 !; d2ew=
    d3ew=cmplx(0d0,0d0) 

! reciprocal space term  
    do igrid=1,ng  ! G=0 is needed to assure periodicity
       qpg=ggrid(:,igrid)+q
       geg=dot_product(qpg,matmul(epsil,qpg))
!      if (geg.gt.100*etaew*etaew) cycle 
       if (geg.gt.4*sumcut*sumcut*etaew*etaew) cycle 
       call naq(qpg,tau,taup,etaew,hout,gout)
       if (ewald_no_macro .and. dot_product(ggrid(:,igrid),ggrid(:,igrid)).lt.1d-20) then
          call naq(qpg,tau,taup,1d6,hmac,gmac)       ! the bare G = 0 term (exp -> 1), removed
          hout = hout - hmac
          gout = gout - gmac
       endif
       d2ew = d2ew + hout
       do ga=1,3
          d3ew(:,:,ga)=d3ew(:,:,ga) + gout(:,:,ga)  
       enddo
    enddo
    gsum=trace(d2ew)

    ep3 = etaew*etaew*etaew/sqrt(det(epsil))/eps0scale

! real space term  
    hout=0
    do igrid=1,nr
       sd=(dta+rgrid(:,igrid))
       del=matmul(epsinv,sd)
       dd=sqrt(dot_product(sd,del))
       dd2=sqrt(dot_product(rgrid(:,igrid),matmul(epsinv,rgrid(:,igrid))))
       if(dd.lt.1d-8) cycle ! exclude self-interactions
 !     if (etaew*dd.gt.6) cycle
       if (etaew*dd2.gt.sumcut) cycle ! to keep sums even in Rgrid
       zz = exp( ci*dot_product(q,rgrid(:,igrid)))  !-dta)) 
       hout= hfunc(etaew*del,etaew*dd) * zz *ep3
       d2ew=d2ew - hout
       do ga=1,3
          d3ew(:,:,ga)=d3ew(:,:,ga) - ci*rgrid(ga,igrid) *hout  
       enddo
    enddo
    rpgsum=trace(d2ew)
!   write(ulog,5)'G_hat:tau,taup,q, gsum,rpgsum,const=',tau,taup,q, gsum,rpgsum,consta 

    if(tau.eq.taup) d2ew=d2ew - 4/3d0/sqrt(pi)*epsinv *ep3  
    consta=trace(d2ew)

!  write(*,3)'q_red,G_sum, G+R_sum, G+R+Const trace=',cart2red(q,'g'), gsum,rpgsum,consta

3 format(a,99(1x,g11.4))
4 format(a,i4,99(1x,g11.4))
5 format(a,2i4,99(1x,g11.4))

 end subroutine ewald_2nd_deriv_hat
!--------------------------------------------------
 subroutine ewald_2nd_deriv_Gnew(q,tau,taup,ng,ggrid,etaew,d2ew,d3ew)
!! calculates the dynamical matrix associated with the NA term (only G sum) 
!! input is q, two atoms tau and taup, output is the 3x3 block D^EW_{tau,taup}(q) 
!! which goes to D^NA for q \to 0 ; G are reciprocal lattice vectors (not supercell)
!! uses the step phase convention
!! D_{t,t'}(q)=sum_G exp[i(q-G).dta] f(q-G) with f(q)=q^al q^be/qeq exp[-qeq/4/eta^2]
! T(q,x)=sum_R e^iq.R /|X+R| =4pi/vol sum_G exp[i(G-q).X]/eps|q-G|^2 exp[-q-G.eps.q-G/4/eta^2] 
! D(q)=\nabla_X^2 T(q,t'-t) with X=t'-t
 use geometry, only : trace
 use born , only : epsil
 use atoms_force_constants, only : natom_prim_cell
 use ios , only : ulog
 implicit none
 integer, intent(in) :: tau,taup,ng
 real(r15), intent(in) :: q(3),ggrid(3,ng),etaew
 complex(r15), intent(out) :: d2ew(3,3),d3ew(3,3,3)
 integer, external :: delta_k
 complex(r15)  hout(3,3),gout(3,3,3),gsum
 integer igrid,al,be,ga,tau2
 real(r15) geg,qpg(3),tol,sumcut

 tol=1d-12; sumcut=sqrt(log(1d0/tol))
 write(ulog,*)'_Gnew: ng,Gcut=',ng,2*sumcut*etaew/sqrt(trace(epsil)/3)
    d2ew= cmplx(0d0,0d0) ! ; epsinv*coef/3 ; 
    d3ew=0
    gloop: do igrid=1,ng
       qpg=ggrid(:,igrid)+q
       geg=dot_product(qpg,matmul(epsil,qpg))
       if (geg.lt.1d-20) cycle gloop
!      if (geg.gt.100 *etaew*etaew) cycle gloop  
       if (geg.gt.4*sumcut*sumcut*etaew*etaew) cycle gloop

       call naq(qpg,tau,taup,etaew,hout,gout)

       d2ew=d2ew + hout 
       do ga=1,3
          d3ew(:,:,ga)=d3ew(:,:,ga) + gout(:,:,ga)
       enddo
!      write(6,4)'rgrid,geg,hout=',igrid,geg,hout !(1,1)
    enddo gloop
    gsum=trace(d2ew)
    write(ulog,4)'last term in G-sum in ewa_2nd_der_Gnew=',igrid,d2ew(1,:),hout(1,:)
    write(ulog,5)'G_new:tau,taup,q, trace=',tau,taup,q, gsum

!  return

! subtract the self-interaction before the final charge multiplication 
    if(tau.eq.taup) then
       do tau2=1,natom_prim_cell
       g0loop: do igrid=2,ng ! exclude G=0
          qpg=ggrid(:,igrid)
          geg=dot_product(qpg,matmul(epsil,qpg))
!         if (geg.gt.100 *etaew*etaew) cycle g0loop  
          if (geg.gt.4*sumcut*sumcut*etaew*etaew) cycle g0loop

          call naq(qpg,tau,tau2,etaew,hout,gout)

          d2ew=d2ew - hout 
          do ga=1,3
             d3ew(:,:,ga)=d3ew(:,:,ga) - gout(:,:,ga)
          enddo
       enddo g0loop
!      write(ulog,4)'last term in G-sum in ewa_2nd_der_q=0 =',igrid,d2ew(1,:),hout(1,:)
       enddo
    endif

    write(ulog,5)'G_new:tau,taup,q, trace=',tau,taup,q, trace(d2ew)

3 format(a,99(1x,g14.7))
4 format(a,i4,99(1x,g11.4))
5 format(a,2i4,99(1x,g11.4))

 end subroutine ewald_2nd_deriv_Gnew

 end module ewald
