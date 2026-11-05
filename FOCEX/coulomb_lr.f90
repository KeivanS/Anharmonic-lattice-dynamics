!============================================================
 module coulomb_lr
!! Long-range (Born-charge) corrections of the invariance rows for born_flag 2 (2026-11-03, user's design).
!! born_flag 2 subtracts from the DFT forces the force of the Coulomb dipole FCs Phi_C (+ the anharmonic part of a
!! Coulomb pair model, see subtract_coulomb_force), so the fit sees x_fit = x_total - x_sub and the rows
!! A x_total = b become A x_fit = b - A x_sub. x_sub has no FC1 (the reference force is subtracted) and its FC2 is
!! Phi_C(t a, u b; R) = e^2 Z_t T(R + x_u - x_t) Z_u^T, T = -grad grad phi_eps (the D_NA added back, E = 0).
!!   Huang:            A x_sub = S_C(ab,gd) - S_C(gd,ab),  S(ab,gd) = sum_{t in cell, j} Phi(t a, j b) r_g r_d
!!   rank-2 rotation:  A x_sub = sum_{a,m} M(t; a1, a, m) eps(a, m, nu),  M = sum_j Phi(t a1, j a) r_m
!! Both lattice sums are conditionally convergent; they are taken from the small-q expansion of
!! C_tu(q) = sum_R Phi_C exp(i q.(R + x_u - x_t)) (Ewald, E = 0: the bare G = 0 term left out):
!!   S(ab,gd) = - d^2/dq_g dq_d sum_tu Re C_ab(q) at 0,   M(t;a1,a,m) = d/dq_m sum_u Im C_t,u(q)_(a1,a) at 0,
!! by central differences with Richardson extrapolation. Same kernel as REFERENCE/exact_model/longwave/
!! coul_brackets.py (checked there: independent of the Ewald parameter to 4e-11).
!! Rank-3 rotational rows (lr_rot3_pair, lr_rot3_onsite): b = -A_full x_sub, x_sub = Phi_C (rank 2) + the FC3 of the
!! pair model c_ij phi_eps, c_ij = tr(Z_i Z_j^T)/3 (rank 3). Row (i1 a1, i2 a2; nu), r = X_i1 - X_i2, e2 = 1/eps0scale:
!!   pair (i2 /= i1):  b = e2 [ c phi_(a1 a2 a)(r) r_mu eps(a,mu,nu) + P_(g a2) eps(g,a1,nu) + P_(a1 g) eps(g,a2,nu) ],
!!                     P = Z_i1 H Z_i2^T - c H, H = grad grad phi_eps (bare: one pair of the crystal);
!!   on-site (i2 = i1): b = - sum_(l /= i1) b_pair(i1,l) (ASR of x_sub), a conditionally convergent sum, taken as the
!!                     E = 0 Ewald lattice sum. Its first part uses c phi_(a1 a2 a) r_mu eps(a,mu,nu) + c (H eps terms)
!!                     = c D_eps H [eps J_nu - J_nu eps], J_nu(a,mu) = eps(a,mu,nu) (rotating r at fixed eps =
!!                     rotating eps the other way), so it is the eps-derivative of the E = 0 sum of c H.
!! For scalar Z and eps all of it vanishes (the pair model is then rotationally invariant and Phi_C = Phi_m).
 use constants, only : r15, pi, ci, eps0scale
 implicit none
 private
 public :: lr_c, lr_brackets, lr_moments, lr_moments_pairs, lr_rot3_pair, lr_rot3_onsite

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
 real(r15) a(3,3), b(3,3), w(3), rc, gc, dmax, r(3)
 integer i, j, k, n(3), m
 if (ready) return
 nat = natom_prim_cell
 allocate(xpos(3,nat), zc(3,3,nat))
 do i=1,nat
    xpos(:,i) = v2a(atom0(i)%equilibrium_pos)
    zc(:,:,i) = atom0(i)%charge
 enddo
 a(:,1)=v2a(r01); a(:,2)=v2a(r02); a(:,3)=v2a(r03)
 vol = volume_r0
 call set_eps((epsil+transpose(epsil))/2, w)
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
 subroutine set_eps(e, w)
!! eps and the derived eih = eps^-1/2, ei = eps^-1, sdet = sqrt(det eps); w its eigenvalues
 use geometry, only : det
 real(r15), intent(in) :: e(3,3)
 real(r15), intent(out) :: w(3)
 real(r15) v(3,3), rwork(7)
 complex(r15) cv(3,3), work(12)
 integer info
 eps = e
 cv = eps
 call zheev('V','U',3,cv,3,w,work,12,rwork,info)       ! (zhegv.f; no LAPACK linked)
 v = real(cv)
 eih = matmul(v, matmul(reshape((/1/sqrt(w(1)),0d0,0d0, 0d0,1/sqrt(w(2)),0d0, 0d0,0d0,1/sqrt(w(3))/),(/3,3/)), &
 &     transpose(v)))
 ei = matmul(eih,eih)
 sdet = sqrt(det(eps))
 end subroutine set_eps
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
 complex(r15) k(3,3,nat,nat)
 integer it, iu
 call lr_k(q, k)
 do it=1,nat
 do iu=1,nat
    c(:,:,it,iu) = matmul(zc(:,:,it), matmul(k(:,:,it,iu), transpose(zc(:,:,iu))))
 enddo
 enddo
 end subroutine lr_c
!------------------------------------------------------------
 subroutine lr_k(q, kk)
!! e2 sum_R T(R + x_u - x_t) exp(i q.(R + x_u - x_t)), T = -grad grad phi_eps (E = 0; r = 0 left out), eV/A^2
 real(r15), intent(in) :: q(3)
 complex(r15), intent(out) :: kk(3,3,nat,nat)
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
    kk(:,:,it,iu) = k/eps0scale
 enddo
 enddo
 end subroutine lr_k
!------------------------------------------------------------
 subroutine bare_h(r, h, h3)
!! h = grad grad phi_eps(r), h3(a,b,c) = its gradient; phi_eps = 1/(sqrt(det eps) |eps^-1/2 r|), unscreened, one pair
 real(r15), intent(in) :: r(3)
 real(r15), intent(out) :: h(3,3), h3(3,3,3)
 real(r15) x(3), rho, hx(3,3), tx(3,3,3)
 integer i, j, k, a, b, c
 x = matmul(eih, r); rho = sqrt(sum(x**2))
 do i=1,3
 do j=1,3
    hx(i,j) = 3*x(i)*x(j)/rho**5
    if (i.eq.j) hx(i,j) = hx(i,j) - 1/rho**3
    do k=1,3
       tx(i,j,k) = -15*x(i)*x(j)*x(k)/rho**7
       if (i.eq.j) tx(i,j,k) = tx(i,j,k) + 3*x(k)/rho**5
       if (i.eq.k) tx(i,j,k) = tx(i,j,k) + 3*x(j)/rho**5
       if (j.eq.k) tx(i,j,k) = tx(i,j,k) + 3*x(i)/rho**5
    enddo
 enddo
 enddo
 h = matmul(eih, matmul(hx, eih))/sdet           ! eih symmetric: d/dr_a = sum_i eih(i,a) d/dx_i
 h3 = 0
 do a=1,3; do b=1,3; do c=1,3
    do i=1,3; do j=1,3; do k=1,3
       h3(a,b,c) = h3(a,b,c) + eih(i,a)*eih(j,b)*eih(k,c)*tx(i,j,k)
    enddo; enddo; enddo
 enddo; enddo; enddo
 h3 = h3/sdet
 end subroutine bare_h
!------------------------------------------------------------
 subroutine levi(lc)
 integer, intent(out) :: lc(3,3,3)
 lc = 0
 lc(1,2,3) =  1 ; lc(2,3,1) =  1 ; lc(3,1,2) =  1
 lc(1,3,2) = -1 ; lc(2,1,3) = -1 ; lc(3,2,1) = -1
 end subroutine levi
!------------------------------------------------------------
 subroutine lr_rot3_pair(r, zi, zj, b)
!! b(a1,a2,nu): right-hand side of the rank-3 rotational row (i1 a1, i2 a2; nu), i2 /= i1, r = X_i1 - X_i2 (eV/A)
 real(r15), intent(in) :: r(3), zi(3,3), zj(3,3)
 real(r15), intent(out) :: b(3,3,3)
 real(r15) h(3,3), h3(3,3,3), p(3,3), c
 integer lc(3,3,3), a1, a2, nu, a, mu, g
 call setup
 call levi(lc)
 call bare_h(r, h, h3)
 c = sum(zi*zj)/3
 p = matmul(zi, matmul(h, transpose(zj)))        ! c W = c (h3.r eps + H eps terms), plus P = Z H Z^T - c H
 b = 0
 do nu=1,3; do a1=1,3; do a2=1,3
    do a=1,3; do mu=1,3
       b(a1,a2,nu) = b(a1,a2,nu) + c*h3(a1,a2,a)*r(mu)*lc(a,mu,nu)
    enddo; enddo
    do g=1,3
       b(a1,a2,nu) = b(a1,a2,nu) + p(g,a2)*lc(g,a1,nu) + p(a1,g)*lc(g,a2,nu)
    enddo
 enddo; enddo; enddo
 b = b/eps0scale
 end subroutine lr_rot3_pair
!------------------------------------------------------------
 subroutine lr_rot3_onsite(bo, ulog)
!! bo(t,a1,a2,nu): right-hand side of the rank-3 rotational row (t a1, t a2; nu), t in the primitive cell (eV/A):
!! - sum_(l /= t) b_pair(t,l) as E = 0 lattice sums (see the module header). Logs a check of the eps-derivative
!! identity on one pair (bare): c D_eps H [eps J - J eps] = c (h3.r eps + H eps terms).
 integer, intent(in) :: ulog
 real(r15), intent(out) :: bo(:,:,:,:)
 complex(r15), allocatable :: k0(:,:,:,:), kp(:,:,:,:), km(:,:,:,:)
 real(r15) e0(3,3), de(3,3), w(3), he, p(3,3), dk(3,3), cc, jn(3,3), h(3,3), h3(3,3,3), hp(3,3), hm(3,3), r(3), &
 &         wa(3,3), wb(3,3), err
 integer lc(3,3,3), it, iu, nu, a1, a2, g, a, mu
 call setup
 call levi(lc)
 allocate(k0(3,3,nat,nat), kp(3,3,nat,nat), km(3,3,nat,nat))
 e0 = eps
 he = 1d-4
 call lr_k((/0d0,0d0,0d0/), k0)
 bo = 0
 err = 0
 do nu=1,3
    jn = lc(:,:,nu)
    de = matmul(e0,jn) - matmul(jn,e0)
    call set_eps(e0+he*de, w); call lr_k((/0d0,0d0,0d0/), kp)
    call set_eps(e0-he*de, w); call lr_k((/0d0,0d0,0d0/), km)
    call set_eps(e0, w)
    do it=1,nat
       dk = 0; p = 0
       do iu=1,nat
          cc = sum(zc(:,:,it)*zc(:,:,iu))/3
          dk = dk + cc*real(kp(:,:,it,iu)-km(:,:,it,iu))/(2*he)
          p = p + matmul(zc(:,:,it), matmul(real(k0(:,:,it,iu)), transpose(zc(:,:,iu)))) - cc*real(k0(:,:,it,iu))
       enddo
! sum_l e2 c H = -sum_u c K, sum_l e2 P = -sum_u (Z K Z^T - c K): bo = -sum_l b_pair = D_eps(sum_u c K) + P eps terms
       do a1=1,3; do a2=1,3
          bo(it,a1,a2,nu) = dk(a1,a2)
          do g=1,3
             bo(it,a1,a2,nu) = bo(it,a1,a2,nu) + p(g,a2)*lc(g,a1,nu) + p(a1,g)*lc(g,a2,nu)
          enddo
       enddo; enddo
    enddo
! check of the identity on one bare pair (nearest atom 2, else an image of atom 1)
    r = xpos(:,1) - xpos(:,min(2,nat)) + rlist(:,2)
    call bare_h(r, h, h3)
    call set_eps(e0+he*de, w); call bare_h(r, hp, h3)
    call set_eps(e0-he*de, w); call bare_h(r, hm, h3)
    call set_eps(e0, w); call bare_h(r, h, h3)
    wa = (hp-hm)/(2*he)
    wb = 0
    do a1=1,3; do a2=1,3
       do a=1,3; do mu=1,3
          wb(a1,a2) = wb(a1,a2) + h3(a1,a2,a)*r(mu)*lc(a,mu,nu)
       enddo; enddo
       do g=1,3
          wb(a1,a2) = wb(a1,a2) + h(g,a2)*lc(g,a1,nu) + h(a1,g)*lc(g,a2,nu)
       enddo
    enddo; enddo
    err = max(err, maxval(abs(wa-wb))/max(maxval(abs(h)),1d-30))
 enddo
 write(ulog,'(a,es10.2)')' LR_ROT3_ONSITE: eps-derivative identity on a bare pair, max|D_eps H - (h3.r + H eps)|/max|H| =',err
 deallocate(k0, kp, km)
 end subroutine lr_rot3_onsite
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
!! s(a,b,g,d) = sum_{t in cell, j} Phi_C(t a, j b) r_g r_d (eV), Phi_C the Coulomb dipole FCs (E = 0)
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
!! mo(t, a1, a, m) = sum_j Phi_C(t a1, j a) (x_j - x_t)_m (eV/A), t in the primitive cell
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
!------------------------------------------------------------
 subroutine lr_moments_pairs(mo)
!! mo(t, u, a1, a, m) = sum_{j in sublattice u} Phi_C(t a1, j a) (x_j - x_t)_m (eV/A), t, u in the
!! primitive cell: lr_moments before the sum over u (2026-11-04; relax_reference needs them for the
!! internal displacements, r = F r0 + u0_u - u0_t)
 real(r15), intent(out) :: mo(:,:,:,:,:)
 complex(r15), allocatable :: cp(:,:,:,:), cm(:,:,:,:), c2p(:,:,:,:), c2m(:,:,:,:)
 real(r15) h, e(3)
 integer m, it, iu
 call setup
 allocate(cp(3,3,nat,nat), cm(3,3,nat,nat), c2p(3,3,nat,nat), c2m(3,3,nat,nat))
 h = 0.04d0/vol**(1d0/3)
 do m=1,3
    e = 0; e(m) = 1
    call lr_c(h*e, cp); call lr_c(-h*e, cm); call lr_c(2*h*e, c2p); call lr_c(-2*h*e, c2m)
    do it=1,nat
    do iu=1,nat
       mo(it,iu,:,:,m) = (4*aimag(cp(:,:,it,iu)-cm(:,:,it,iu))/(2*h) - aimag(c2p(:,:,it,iu)-c2m(:,:,it,iu))/(4*h))/3
    enddo
    enddo
 enddo
 end subroutine lr_moments_pairs

 end module coulomb_lr
