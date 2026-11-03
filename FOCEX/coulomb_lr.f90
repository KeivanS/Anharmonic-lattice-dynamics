!============================================================
 module coulomb_lr
!! Long-range (Born-charge) corrections of the invariance rows for born_flag 2 (2026-11-03, user's design).
!! born_flag 2 subtracts from the DFT forces the force of the Born-dipole FCs Phi_B (+ the anharmonic part of a
!! Coulomb pair model, see subtract_coulomb_force), so the fit sees x_fit = x_total - x_sub and the rows
!! A x_total = b become A x_fit = b - A x_sub. x_sub has no FC1 (the reference force is subtracted) and its FC2 is
!! Phi_B(t a, u b; R) = e^2 Z_t T(R + x_u - x_t) Z_u^T, T = -grad grad phi_eps (the D_NA added back, E = 0).
!!   Huang:            A x_sub = S_B(ab,gd) - S_B(gd,ab),  S(ab,gd) = sum_{t in cell, j} Phi(t a, j b) r_g r_d
!!   rank-2 rotation:  A x_sub = sum_{a,m} M(t; a1, a, m) eps(a, m, nu),  M = sum_j Phi(t a1, j a) r_m
!! Both lattice sums are conditionally convergent; they are taken from the small-q expansion of
!! C_tu(q) = sum_R Phi_B exp(i q.(R + x_u - x_t)) (Ewald, E = 0: the bare G = 0 term left out):
!!   S(ab,gd) = - d^2/dq_g dq_d sum_tu Re C_ab(q) at 0,   M(t;a1,a,m) = d/dq_m sum_u Im C_t,u(q)_(a1,a) at 0,
!! by central differences with Richardson extrapolation. Same kernel as REFERENCE/exact_model/longwave/
!! coul_brackets.py (checked there: independent of the Ewald parameter to 4e-11).
 use constants, only : r15, pi, ci, eps0scale
 implicit none
 private
 public :: lr_c, lr_brackets, lr_moments

 logical, save :: ready=.false.
 integer, save :: nrs=0, ngs=0
 real(r15), allocatable, save :: rlist(:,:), glist(:,:)
 real(r15), save :: lam, sdet, eih(3,3), ei(3,3), eps(3,3), vol
 integer, save :: nat
 real(r15), allocatable, save :: xpos(:,:), zc(:,:,:)

 contains
!------------------------------------------------------------
 subroutine setup
 use born, only : epsil
 use lattice, only : r01, r02, r03, volume_r0
 use atoms_force_constants, only : atom0, natom_prim_cell
 use geometry, only : v2a, det
 real(r15) a(3,3), b(3,3), w(3), v(3,3), rwork(7), rc, gc, dmax, r(3)
 complex(r15) cv(3,3), work(12)
 integer i, j, k, n(3), info, m
 if (ready) return
 nat = natom_prim_cell
 allocate(xpos(3,nat), zc(3,3,nat))
 do i=1,nat
    xpos(:,i) = v2a(atom0(i)%equilibrium_pos)
    zc(:,:,i) = atom0(i)%charge
 enddo
 a(:,1)=v2a(r01); a(:,2)=v2a(r02); a(:,3)=v2a(r03)
 vol = volume_r0
 eps = (epsil+transpose(epsil))/2
 cv = eps
 call zheev('V','U',3,cv,3,w,work,12,rwork,info)       ! (zhegv.f; no LAPACK linked)
 v = real(cv)
 eih = matmul(v, matmul(reshape((/1/sqrt(w(1)),0d0,0d0, 0d0,1/sqrt(w(2)),0d0, 0d0,0d0,1/sqrt(w(3))/),(/3,3/)), &
 &     transpose(v)))
 ei = matmul(eih,eih)
 sdet = sqrt(det(eps))
 lam = 2.5d0/vol**(1d0/3)
 dmax = 0
 do i=1,nat
 do j=1,nat
    dmax = max(dmax, sqrt(sum((xpos(:,j)-xpos(:,i))**2)))
 enddo
 enddo
 rc = 6.5d0*sqrt(maxval(w))/lam + dmax        ! erfc(6.5) = 4e-20
 gc = 13d0*lam/sqrt(minval(w)) + 1             ! exp(-42) ; +1 for |q|
! reciprocal vectors as rows: b(i,:).a(:,j) = 2 pi delta_ij
 call inverse3(a, b); b = 2*pi*b
 do i=1,3
    n(i) = int(rc*sqrt(sum(b(i,:)**2))/(2*pi)) + 1
 enddo
 allocate(rlist(3,(2*n(1)+1)*(2*n(2)+1)*(2*n(3)+1)))
 m = 0
 do i=-n(1),n(1); do j=-n(2),n(2); do k=-n(3),n(3)
    r = i*a(:,1)+j*a(:,2)+k*a(:,3)
    if (sqrt(sum(r**2)).gt.rc) cycle
    m = m+1; rlist(:,m) = r
 enddo; enddo; enddo
 nrs = m
 do i=1,3
    n(i) = int(gc*sqrt(sum(a(:,i)**2))/(2*pi)) + 1
 enddo
 allocate(glist(3,(2*n(1)+1)*(2*n(2)+1)*(2*n(3)+1)))
 m = 0
 do i=-n(1),n(1); do j=-n(2),n(2); do k=-n(3),n(3)
    r = i*b(1,:)+j*b(2,:)+k*b(3,:)
    if (sqrt(sum(r**2)).gt.gc) cycle
    m = m+1; glist(:,m) = r
 enddo; enddo; enddo
 ngs = m
 ready = .true.
 end subroutine setup
!------------------------------------------------------------
 subroutine inverse3(a, ai)
 real(r15), intent(in) :: a(3,3)
 real(r15), intent(out) :: ai(3,3)
 real(r15) d
 ai(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2); ai(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3); ai(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
 ai(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3); ai(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1); ai(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
 ai(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1); ai(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2); ai(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
 d = a(1,1)*ai(1,1)+a(1,2)*ai(2,1)+a(1,3)*ai(3,1)
 ai = ai/d
 end subroutine inverse3
!------------------------------------------------------------
 subroutine t_real(r, t)
!! -grad grad [erfc(lam D)/D]/sqrt(det eps), D = |eps^-1/2 r|
 real(r15), intent(in) :: r(3)
 real(r15), intent(out) :: t(3,3)
 real(r15) x(3), rho, e, g, f1, f2, h(3,3)
 integer i
 x = matmul(eih, r); rho = sqrt(sum(x**2))
 e = erfc(lam*rho); g = 2*lam/sqrt(pi)*exp(-(lam*rho)**2)
 f1 = -e/rho**2 - g/rho
 f2 = 2*e/rho**3 + g*(2/rho**2 + 2*lam**2)
 do i=1,3
    h(:,i) = x*x(i)/rho**2*(f2-f1/rho)
    h(i,i) = h(i,i) + f1/rho
 enddo
 t = -matmul(eih, matmul(h, eih))/sdet
 end subroutine t_real
!------------------------------------------------------------
 subroutine lr_c(q, c)
!! C_tu(q) (eV/A^2), full-distance phase, E = 0, blocks c(3,3,t,u); Born charges Z(displacement, field)
 real(r15), intent(in) :: q(3)
 complex(r15), intent(out) :: c(3,3,nat,nat)
 complex(r15) k(3,3)
 real(r15) d(3), r(3), t(3,3), kv(3), kek, fac
 integer it, iu, m, i
 call setup
 do it=1,nat
 do iu=1,nat
    d = xpos(:,iu)-xpos(:,it)
    k = 0
    do m=1,nrs
       r = rlist(:,m)+d
       if (sqrt(sum(r**2)).lt.1d-8) cycle
       call t_real(r, t)
       k = k + t*exp(ci*dot_product(q,r))
    enddo
    do m=1,ngs
       kv = q+glist(:,m); kek = dot_product(kv, matmul(eps,kv))
       if (kek.lt.1d-14) cycle
       fac = 4*pi/vol/kek*exp(-kek/(4*lam**2))
       if (sum(glist(:,m)**2).lt.1d-20) fac = fac - 4*pi/vol/kek     ! E = 0: bare G = 0 term removed
       do i=1,3
          k(:,i) = k(:,i) + fac*kv*kv(i)*exp(-ci*dot_product(glist(:,m),d))
       enddo
    enddo
    if (it.eq.iu) k = k - 4*lam**3/(3*sqrt(pi))*ei/sdet               ! Gaussian self term
    c(:,:,it,iu) = matmul(zc(:,:,it), matmul(k, transpose(zc(:,:,iu))))/eps0scale
 enddo
 enddo
 end subroutine lr_c
!------------------------------------------------------------
 subroutine csum(q, s)
!! sum over t,u of C_tu(q)
 real(r15), intent(in) :: q(3)
 complex(r15), intent(out) :: s(3,3)
 complex(r15) c(3,3,nat,nat)
 integer it, iu
 call lr_c(q, c)
 s = 0
 do it=1,nat
 do iu=1,nat
    s = s + c(:,:,it,iu)
 enddo
 enddo
 end subroutine csum
!------------------------------------------------------------
 subroutine lr_brackets(s)
!! s(a,b,g,d) = sum_{t in cell, j} Phi_B(t a, j b) r_g r_d (eV), Phi_B the Born-dipole FCs (E = 0)
 real(r15), intent(out) :: s(3,3,3,3)
 real(r15) h, eg(3), ed(3), d1(3,3), d2(3,3)
 complex(r15) f0(3,3)
 integer g, dd
 call setup
 h = 0.04d0/vol**(1d0/3)
 call csum((/0d0,0d0,0d0/), f0)
 do g=1,3
 do dd=g,3
    eg = 0; eg(g) = 1; ed = 0; ed(dd) = 1
    call second(h, d1); call second(2*h, d2)
    s(:,:,g,dd) = -(4*d1-d2)/3
    s(:,:,dd,g) = s(:,:,g,dd)
 enddo
 enddo
 contains
   subroutine second(hh, out)
   real(r15), intent(in) :: hh
   real(r15), intent(out) :: out(3,3)
   complex(r15) p(3,3), m(3,3), pm(3,3), mp(3,3)
   if (g.eq.dd) then
      call csum(hh*eg, p); call csum(-hh*eg, m)
      out = real(p+m-2*f0)/hh**2
   else
      call csum(hh*(eg+ed), p); call csum(-hh*(eg+ed), m)
      call csum(hh*(eg-ed), pm); call csum(-hh*(eg-ed), mp)
      out = real(p+m-pm-mp)/(4*hh**2)
   endif
   end subroutine second
 end subroutine lr_brackets
!------------------------------------------------------------
 subroutine lr_moments(mo)
!! mo(t, a1, a, m) = sum_j Phi_B(t a1, j a) (x_j - x_t)_m (eV/A), t in the primitive cell
 real(r15), intent(out) :: mo(:,:,:,:)
 complex(r15), allocatable :: cp(:,:,:,:), cm(:,:,:,:), c2p(:,:,:,:), c2m(:,:,:,:)
 real(r15) h, e(3), d1(3,3), d2(3,3)
 integer m, it, iu
 call setup                                  ! sets nat
 allocate(cp(3,3,nat,nat), cm(3,3,nat,nat), c2p(3,3,nat,nat), c2m(3,3,nat,nat))
 h = 0.04d0/vol**(1d0/3)
 do m=1,3
    e = 0; e(m) = 1
    call lr_c(h*e, cp); call lr_c(-h*e, cm); call lr_c(2*h*e, c2p); call lr_c(-2*h*e, c2m)
    do it=1,nat
       d1 = 0; d2 = 0
       do iu=1,nat
          d1 = d1 + aimag(cp(:,:,it,iu)-cm(:,:,it,iu))/(2*h)
          d2 = d2 + aimag(c2p(:,:,it,iu)-c2m(:,:,it,iu))/(4*h)
       enddo
       mo(it,:,:,m) = (4*d1-d2)/3
    enddo
 enddo
 end subroutine lr_moments

 end module coulomb_lr
