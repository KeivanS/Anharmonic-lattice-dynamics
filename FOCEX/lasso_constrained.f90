!===============================================================================
! Constrained LASSO for the force constants, followed by a constrained least-squares
! refit on the selected terms, with the L1 weight chosen by cross-validation.
!
!   minimise  1/(2N) ||A x - b||^2 + lam sum_j w_j |x_j|   subject to  C x = d
!
! d = 0 except for rows with a right-hand side (Huang rows of a stressed reference): then
! x = x0 + K y with x0 a particular solution, and z = w*x0 + Q t below (z0 = w*x0).
!
! x are the irreducible force constants themselves (not the kernel coefficients), so
! the L1 term makes the FORCE CONSTANTS sparse while the invariances hold exactly.
! w_j = rms of column j of A (standardisation: the penalty does not depend on the
! units of each rank). The constraint C x = 0 enters as x in range(K), K the kernel
! basis FOCEX builds (get_kernel).
!
! Variables: z = w*x (standardised FCs). range(W K) has the orthonormal basis Q
! (left singular vectors of W K), z = Q t, and the problem becomes
!   minimise  1/2 t^T G t - c^T t + lam ||Q t||_1,   G = At^T At, c = At^T q,
!   At = (A/w) Q / sqrt(N), q = b / sqrt(N).
! Solved by a primal-dual interior-point method (Mehrotra predictor-corrector) on
! the QP form, see ipm_classo, to a duality gap of tol relative. The support S is
! |z_i| > sel_tol max|z|. Interior-point zeros are small but not exact: at small lam
! a few terms that are zero in the exact solution pass this test (against LARS on
! unconstrained problems: 57 kept for 46 at lam = 1e-4 lam_max, no true term lost; the
! dual multipliers separate them only by a factor 2 there, and tol below ~1e-12 is not
! reachable in double precision). The rule errs on the side of keeping terms: a lost
! true term makes the refit collapse, an extra one is fitted to ~0 by it. ADMM (Boyd et
! al., Found. Trends Mach. Learn. 3, 1 (2011)) was tried first: on FOCEX's exact cubic
! model it had not converged after 200000 iterations at 3 of 5 lam; the interior point
! needs 11-16 and agrees with cvxpy/Clarabel to 1e-11 in the objective.
!
! Refit (refit=.true., enforce_inv=3): least squares over x with x_j = 0 off S and
! C x = 0, i.e. t in the null space
! of Q(S^c,:). This removes the shrinkage bias of the LASSO ("relaxed lasso",
! Meinshausen, Comput. Stat. Data Anal. 52, 374 (2007)). On FOCEX's exact cubic
! model the unconstrained selection followed by this refit lost true FCs and the
! refit collapsed (kernel 3); the constrained selection plus refit gave the exact FCs.
!
! lam grid: lam_max = ||Q c||_inf (at and above it t = 0 solves the problem: nu = Q c
! is a feasible dual point), down to epsratio*lam_max, log-spaced (extended by two decades
! at a time, to 1e-8, while the CV minimum is at its end: classo_cv). k-fold cross-
! validation, interleaved rows; the lam of smallest mean held-out error of the
! REFITTED model is chosen, and the final answer is the refit on all rows there.
!
! Uses the external svdcmp (svd_ridge.f90). Validated standalone (scratch lassodev,
! drv_classo) on the exact cubic model dumped from FOCEX: objective = cvxpy/Clarabel
! to 1e-11, refit = numpy least squares to 1e-9, CV + refit returns the exact FCs.
!
! K. Esfarjani, 2026
!===============================================================================
 module lasso_constrained
 implicit none
 private
 public classo_cv, classo_path, ipm_classo, refit_support, classo_setup
 integer, parameter :: dp = kind(1d0)
 real(dp), parameter :: sel_tol = 1d-8, null_tol = 1d-9, refit_cut = 1d-13, w_floor = 1d-8, &
&     eps_min = 1d-8

 contains
!-------------------------------------------------------------------------------
 subroutine classo_setup(nrow, ncol, nker, A, b, Kb, nfold, fold, w, Q, Gf, cf, bf, nf, info)
!! column weights w, orthonormal basis Q of range(W Kb), and per-fold Gram data in t:
!! Gf(:,:,f) = sum over rows of fold f of a a^T, a = Q^T (A(i,:)/w)/sqrt(N), etc.
!! info=1 if W Kb is rank deficient (should not happen: Kb has full column rank).
 integer, intent(in) :: nrow, ncol, nker, nfold, fold(nrow)
 real(dp), intent(in) :: A(nrow,ncol), b(nrow), Kb(ncol,nker)
 real(dp), intent(out) :: w(ncol), Q(ncol,nker), Gf(nker,nker,nfold), cf(nker,nfold), bf(nfold)
 integer, intent(out) :: nf(nfold), info
 real(dp), allocatable :: V(:,:), sv(:), row(:)
 real(dp) sn
 integer i, j, f
 info = 0
 sn = sqrt(dble(nrow))
 do j = 1, ncol
    w(j) = sqrt(dot_product(A(:,j), A(:,j)))/sn
 enddo
! floor: columns the snapshots do not see (Ge: 0 and 2e-20 next to genuine ones down to
! 3e-6 of the largest) would make x = (Q t)/w divide round-off by ~1e-20 (Ge: FC2 of 153
! where the SVD has 1e-6)
 if (maxval(w) <= 0) then
    w = 1
 else
    w = max(w, w_floor*maxval(w))
 endif
 allocate(V(nker,nker), sv(nker), row(nker))
 do j = 1, nker
    Q(:,j) = w*Kb(:,j)
 enddo
 call svdcmp(Q, ncol, nker, sv, V)
 if (minval(sv) <= 1d-12*maxval(sv)) info = 1
 Gf = 0 ; cf = 0 ; bf = 0 ; nf = 0
 do i = 1, nrow
    f = fold(i) ; nf(f) = nf(f) + 1
    row = matmul(A(i,:)/w, Q)/sn
    do j = 1, nker
       if (row(j) == 0) cycle
       Gf(:,j,f) = Gf(:,j,f) + row*row(j)
    enddo
    cf(:,f) = cf(:,f) + row*b(i)/sn
    bf(f) = bf(f) + b(i)**2/nrow
 enddo
 deallocate(V, sv, row)
 end subroutine classo_setup
!-------------------------------------------------------------------------------
 subroutine ipm_classo(n, nx, G, c, Q, z0, lam, tol, maxit, t, z, gap, nit, conv)
!! primal-dual interior point (Mehrotra predictor-corrector) for
!!   min 1/2 t^T G t - c^T t + lam ||z0 + Q t||_1   (Q: nx x n, orthonormal columns; z0 = 0
!!   unless some invariance rows have a right-hand side, x = x0 + K y)
!! as the QP  min 1/2 t^T G t - c^T t + lam 1^T s,  a = s - Q t >= 0,  b = s + Q t >= 0,
!! multipliers al, be >= 0. KKT: G t - c + Q^T (al - be) = 0, al + be = lam,
!! al a = be b = 0. Newton step reduced to (G + Q^T D Q) dt = rhs, D = 4 A B/(A+B),
!! A = al/a, B = be/b (positive definite: Q has full column rank).
!! Stops when the dual residuals are below tol (relative) and the duality gap
!! al^T a + be^T b <= tol * (|terms of the objective| + ||c||_inf^2/<G_jj>; the second
!! term keeps the scale finite where the solution tends to 0). gap: final relative gap.
 integer, intent(in) :: n, nx, maxit
 real(dp), intent(in) :: G(n,n), c(n), Q(nx,n), z0(nx), lam, tol
 real(dp), intent(out) :: t(n), z(nx), gap
 integer, intent(out) :: nit
 logical, intent(out) :: conv
 real(dp) s(nx), al(nx), be(nx), a(nx), bv(nx), rt(n), rs(nx), ra(nx), rb(nx)
 real(dp) dt(n), ds(nx), da(nx), db(nx), dal(nx), dbe(nx), A1(nx), B1(nx), Dpl(nx), Dm(nx), Dd(nx)
 real(dp) M(n,n), L(n,n), daa(nx), dba(nx), dala(nx), dbea(nx)
 real(dp) mu, mua, sig, st, sca, zs, tr, obj
 integer it, j
 conv = .false.
! lam >= ||Q c||_inf: t = 0 is the solution (nu = Q c is dual feasible); an interior
! point cannot reach that boundary, so return it directly
 if (all(z0 == 0) .and. maxval(abs(matmul(Q, c))) <= lam) then
    t = 0 ; z = 0 ; gap = 0 ; nit = 0 ; conv = .true. ; return
 endif
 tr = 0
 do j = 1, n ; tr = tr + G(j,j) ; enddo
 zs = max(maxval(abs(c))/max(tr/n, tiny(1d0)), tiny(1d0))
 t = 0 ; s = zs + abs(z0) ; al = lam/2 ; be = lam/2
 conv = .false.
 do it = 1, maxit
    nit = it
    z = z0 + matmul(Q, t) ; a = s - z ; bv = s + z
    rt = matmul(G, t) - c + matmul(al - be, Q)
    rs = lam - al - be
    mu = (dot_product(al, a) + dot_product(be, bv))/(2*nx)
    obj = 0.5d0*dot_product(t, matmul(G, t)) - dot_product(c, t) + lam*sum(abs(z))
    sca = 0.5d0*abs(dot_product(t, matmul(G, t))) + abs(dot_product(c, t)) + lam*sum(abs(z)) + maxval(abs(c))*zs
    gap = 2*nx*mu/sca
    if (maxval(abs(rt)) <= tol*maxval(abs(c)) .and. maxval(abs(rs)) <= tol*lam .and. gap <= tol) then
       conv = .true. ; exit
    endif
    A1 = al/a ; B1 = be/bv ; Dpl = A1 + B1 ; Dm = A1 - B1 ; Dd = 4*A1*B1/Dpl
    do j = 1, n
       M(:,j) = G(:,j) + matmul(Dd*Q(:,j), Q)
    enddo
    call chol_fac(n, M, 0d0, L)
! predictor (affine scaling)
    ra = -al*a ; rb = -be*bv
    call newton
    daa = da ; dba = db ; dala = dal ; dbea = dbe
    st = steplen(1d0)
    mua = (dot_product(al + st*dal, a + st*da) + dot_product(be + st*dbe, bv + st*db))/(2*nx)
    sig = (mua/mu)**3
! corrector
    ra = sig*mu - al*a - dala*daa ; rb = sig*mu - be*bv - dbea*dba
    call newton
    st = steplen(0.99d0)
! beyond what double precision resolves the step turns non-finite: keep the last iterate
    if (.not. (all(abs(dt) <= huge(1d0)) .and. all(abs(ds) <= huge(1d0)) .and. &
&       all(abs(dal) <= huge(1d0)) .and. all(abs(dbe) <= huge(1d0)) .and. st == st)) exit
    t = t + st*dt ; s = s + st*ds ; al = al + st*dal ; be = be + st*dbe
 enddo
 z = z0 + matmul(Q, t)
 contains
    subroutine newton
    real(dp) e(nx), f(nx), qdt(nx), rhs(n)
    e = ra/a + rb/bv - rs
    f = ra/a - rb/bv - Dm*e/Dpl
    rhs = -rt - matmul(f, Q)
    call chol_sol(n, L, rhs, dt)
    qdt = matmul(Q, dt)
    ds = (e + Dm*qdt)/Dpl
    da = ds - qdt ; db = ds + qdt
    dal = (ra - al*da)/a ; dbe = (rb - be*db)/bv
    end subroutine newton
    real(dp) function steplen(frac)
! largest step <= 1 keeping a, b, al, be positive, times frac
    real(dp), intent(in) :: frac
    real(dp) smax
    smax = 1/frac
    smax = min(smax, mstep(a, da), mstep(bv, db), mstep(al, dal), mstep(be, dbe))
    steplen = min(1d0, frac*smax)
    end function steplen
    real(dp) function mstep(v, dv)
    real(dp), intent(in) :: v(nx), dv(nx)
    integer i
    mstep = huge(1d0)
    do i = 1, nx
       if (dv(i) < 0) mstep = min(mstep, -v(i)/dv(i))
    enddo
    end function mstep
 end subroutine ipm_classo
!-------------------------------------------------------------------------------
 subroutine refit_support(n, nx, G, c, Q, z0, sel, t)
!! least squares 1/2 t^T G t - c^T t over t with (z0 + Q t)_i = 0 for .not. sel(i):
!! t = tp + N s, N the null space of B = Q(not sel,:), tp = pinv(B) (-z0(not sel)) (0 if z0 = 0);
!! s = pinv(N^T G N) N^T (c - G tp) (cutoff refit_cut relative on the eigenvalues of N^T G N;
!! minimum norm in the rest).
 integer, intent(in) :: n, nx
 real(dp), intent(in) :: G(n,n), c(n), Q(nx,n), z0(nx)
 logical, intent(in) :: sel(nx)
 real(dp), intent(out) :: t(n)
 real(dp), allocatable :: B(:,:), Nb(:,:), Hm(:,:), Vh(:,:), wh(:), hv(:), s(:), tp(:)
 integer m, i, j, p
 m = count(.not. sel)
 t = 0
 allocate(tp(n)) ; tp = 0
 if (m > 0 .and. any(z0 /= 0)) call shift_part(tp)
 if (m == 0) then
    p = n ; allocate(Nb(n,n)) ; Nb = 0
    do j = 1, n ; Nb(j,j) = 1 ; enddo
 else
! null space of B = Q(not sel,:) by a pivoted Householder QR of B^T (svdcmp did not
! converge on some of these B, Ge example: their singular values cluster at 0 and 1)
    allocate(B(n,m)) ; j = 0
    do i = 1, nx
       if (sel(i)) cycle
       j = j + 1 ; B(:,j) = Q(i,:)
    enddo
    allocate(Nb(n,n))
    call null_space_qr(n, m, B, null_tol, p, Nb)
    deallocate(B)
    if (p == 0) then ; t = tp ; return ; endif
 endif
 allocate(Hm(p,p), Vh(p,p), wh(p), hv(p), s(p))
 Hm = matmul(transpose(Nb(:,1:p)), matmul(G, Nb(:,1:p)))
 hv = matmul(c - matmul(G, tp), Nb(:,1:p))
 call svdcmp(Hm, p, p, wh, Vh)                ! Hm symmetric PSD: H = U W V^T
 s = 0
 do j = 1, p
    if (wh(j) <= refit_cut*maxval(wh)) cycle
    s = s + Vh(:,j)*dot_product(Hm(:,j), hv)/wh(j)
 enddo
 t = tp + matmul(Nb(:,1:p), s)
 deallocate(Nb, Hm, Vh, wh, hv, s, tp)
 contains
    subroutine shift_part(tq)
! a least-squares tq of Q(i,:) tq = -z0(i), i off the support (any one: two differ by a
! vector of the null space N, which the refit optimises over). Basic solution of a
! Householder QR with column pivoting of these rows, rank cut null_tol as in null_space_qr.
! (An SVD of B B^T, B = Q(off,:)^T, did not converge in 18 of 242 calls on the exact wurtzite
! model: its eigenvalues cluster at 0 and 1.)
    real(dp), intent(out) :: tq(n)
    real(dp), allocatable :: Cm(:,:), g(:), cn(:), v(:), y(:), col(:)
    integer, allocatable :: piv(:)
    real(dp) cmax0, al, vn
    integer ii, k, j, jp, r, it
    allocate(Cm(m,n), g(m), cn(n), v(m), y(n), col(m), piv(n))
    j = 0
    do ii = 1, nx
       if (sel(ii)) cycle
       j = j + 1 ; Cm(j,:) = Q(ii,:) ; g(j) = -z0(ii)
    enddo
    do k = 1, n ; piv(k) = k ; cn(k) = sqrt(sum(Cm(:,k)**2)) ; enddo
    cmax0 = maxval(cn) ; r = 0
    do k = 1, min(m, n)
       do j = k, n ; cn(j) = sqrt(sum(Cm(k:m,j)**2)) ; enddo
       jp = k - 1 + maxloc(cn(k:n), 1)
       if (cn(jp) <= null_tol*cmax0) exit
       if (jp /= k) then
          col = Cm(:,k) ; Cm(:,k) = Cm(:,jp) ; Cm(:,jp) = col
          it = piv(k) ; piv(k) = piv(jp) ; piv(jp) = it
       endif
       v = 0 ; v(k:m) = Cm(k:m,k)
       al = -sign(sqrt(sum(v(k:m)**2)), v(k))
       v(k) = v(k) - al
       vn = sqrt(sum(v(k:m)**2)) ; v(k:m) = v(k:m)/vn
       do j = k, n
          Cm(k:m,j) = Cm(k:m,j) - 2*v(k:m)*dot_product(v(k:m), Cm(k:m,j))
       enddo
       g(k:m) = g(k:m) - 2*v(k:m)*dot_product(v(k:m), g(k:m))
       r = k
    enddo
    y = 0                                    ! R11 y(1:r) = (H^T g)(1:r), rest 0
    do k = r, 1, -1
       y(k) = (g(k) - dot_product(Cm(k,k+1:r), y(k+1:r)))/Cm(k,k)
    enddo
    tq = 0
    do k = 1, n ; tq(piv(k)) = y(k) ; enddo
    deallocate(Cm, g, cn, v, y, col, piv)
    end subroutine shift_part
 end subroutine refit_support
!-------------------------------------------------------------------------------
 subroutine classo_path(n, nx, G, c, Q, z0, nlam, lam, tol, maxit, Tl, Tr, nsel, nits, nconv)
!! constrained LASSO at every lam of the decreasing grid (warm starts) and the refit
!! on each support: Tl(:,k) LASSO, Tr(:,k) refit, nsel(k) = |support|, nits(k)
!! interior-point iterations, nconv = number of lam where it did not converge.
 integer, intent(in) :: n, nx, nlam, maxit
 real(dp), intent(in) :: G(n,n), c(n), Q(nx,n), z0(nx), lam(nlam), tol
 real(dp), intent(out) :: Tl(n,nlam), Tr(n,nlam)
 integer, intent(out) :: nsel(nlam), nits(nlam), nconv
 real(dp) t(n), z(nx), gap
 logical sel(nx), selold(nx), conv
 integer k
 nconv = 0 ; selold = .false.
 do k = 1, nlam
    call ipm_classo(n, nx, G, c, Q, z0, lam(k), tol, maxit, t, z, gap, nits(k), conv)
    if (.not. conv) then
       nconv = nconv + 1
       write(*,'(a,i4,a,es10.3,a,es10.3,a,es10.3)') ' WARNING: constrained LASSO: no convergence at lam', k, &
&         ' lam/||Qc||', lam(k)/maxval(abs(matmul(Q, c))), ' gap', gap, ' |t|', maxval(abs(t))
    endif
    Tl(:,k) = t
    sel = abs(z) > sel_tol*maxval(abs(z)) ; nsel(k) = count(sel)
    if (k > 1 .and. all(sel .eqv. selold)) then
       Tr(:,k) = Tr(:,k-1)
    elseif (nsel(k) == 0 .and. all(z0 == 0)) then
       Tr(:,k) = 0
    else
       call refit_support(n, nx, G, c, Q, z0, sel, Tr(:,k))
    endif
    selold = sel
 enddo
 end subroutine classo_path
!-------------------------------------------------------------------------------
 subroutine classo_cv_res(nrow, ncol, nker, A, b, Kb, x0, nfold, nlam, epsratio, tol, maxit, refit, &
&                     x, lam_best, nsel_best, info, uio)
!! the whole procedure. refit: x = refit on all rows at the lam of smallest refit CV
!! error; otherwise x = the constrained LASSO solution at the lam of smallest LASSO CV error.
!! info: 0 ok, 1 W Kb rank deficient, 2 the interior point did not converge somewhere (see log),
!! 3 CV minimum at the smallest lam of the grid. Table written to unit uio (if > 0).
 integer, intent(in) :: nrow, ncol, nker, nfold, nlam, maxit, uio
 real(dp), intent(in) :: A(nrow,ncol), b(nrow), Kb(ncol,nker), x0(ncol), epsratio, tol
 logical, intent(in) :: refit
 real(dp), intent(out) :: x(ncol), lam_best
 integer, intent(out) :: nsel_best, info
 real(dp), allocatable :: z0(:), w(:), Q(:,:), Gf(:,:,:), cf(:,:), G(:,:), c(:), Tl(:,:), Tr(:,:), &
&     lam(:), el(:,:), er(:,:), mr(:), sr(:), ml(:)
 real(dp) bf(nfold), bt, lmax
 integer fold(nrow), nf(nfold), nsel(nlam), nits(nlam), nconv, i, k, f, ibest, inf0
 allocate(w(ncol), Q(ncol,nker), Gf(nker,nker,nfold), cf(nker,nfold), G(nker,nker), c(nker), &
&     Tl(nker,nlam), Tr(nker,nlam), lam(nlam), el(nlam,nfold), er(nlam,nfold), mr(nlam), &
&     sr(nlam), ml(nlam))
 do i = 1, nrow
    fold(i) = mod(i-1, nfold) + 1
 enddo
 call classo_setup(nrow, ncol, nker, A, b - matmul(A, x0), Kb, nfold, fold, w, Q, Gf, cf, bf, nf, inf0)
 info = inf0
 allocate(z0(ncol)) ; z0 = w*x0
 G = sum(Gf, 3) ; c = sum(cf, 2) ; bt = sum(bf)
 lmax = maxval(abs(matmul(Q, c)))
 do k = 1, nlam
    lam(k) = lmax*exp(log(epsratio)*dble(k-1)/dble(max(nlam-1,1)))
 enddo
 if (uio > 0) write(uio,'(a,i7,a,i6,a,i6,a,es12.4)') '# constrained LASSO: rows', nrow, &
&    '  FCs', ncol, '  kernel', nker, '  lam_max', lmax
 do f = 1, nfold
! training Gram = total - fold; lam scaled to the training rows' 1/(2 N_train) form
    call classo_path(nker, ncol, (G - Gf(:,:,f))*nrow/(nrow-nf(f)), (c - cf(:,f))*nrow/(nrow-nf(f)), &
&        Q, z0, nlam, lam, tol, maxit, Tl, Tr, nsel, nits, nconv)
    if (nconv > 0) info = max(info, 2)
    if (uio > 0) write(uio,'(a,i3,a,i4,a,i8)') '# fold', f, '  lam not converged', nconv, &
&       '  IPM iterations', sum(nits)
    do k = 1, nlam
       el(k,f) = vgerr(Tl(:,k)) ; er(k,f) = vgerr(Tr(:,k))
    enddo
 enddo
 do k = 1, nlam
    ml(k) = sum(el(k,:))/nfold ; mr(k) = sum(er(k,:))/nfold
    sr(k) = sqrt(sum((er(k,:) - mr(k))**2)/max(nfold-1,1)/nfold)
 enddo
 if (.not. refit) mr = ml                    ! select on the CV error of the model returned
 ibest = minloc(mr, 1)
! round-off ties (exact data: refit errors ~1e-15 over a range of lam): the largest lam
! within 1e-13 <b^2> of the minimum (1e-10 picked a lam where the exact wurtzite model's refit
! error was still 1e-10 <b^2>, with 3 true FCs off the support)
 do k = 1, ibest
    if (mr(k) <= mr(ibest) + 1d-13*bt) then ; ibest = k ; exit ; endif
 enddo
 if (ibest == nlam) info = max(info, 3)
 lam_best = lam(ibest)
! all rows, down to lam_best
 call classo_path(nker, ncol, G, c, Q, z0, ibest, lam(1:ibest), tol, maxit, Tl(:,1:ibest), &
&     Tr(:,1:ibest), nsel(1:ibest), nits(1:ibest), nconv)
 if (nconv > 0) info = max(info, 2)
 nsel_best = nsel(ibest)
 if (refit) then
    x = x0 + matmul(Q, Tr(:,ibest))/w
 else
    x = x0 + matmul(Q, Tl(:,ibest))/w
 endif
 if (uio > 0) then
    write(uio,'(a)') '#   k        lam/lam_max  CV err LASSO   CV err refit   std err   (all rows:) selected  IPM it'
    do k = 1, nlam
       if (k <= ibest) then
          write(uio,'(i5,es16.6,3es15.6,i10,i9)') k, lam(k)/lmax, ml(k), mr(k), sr(k), nsel(k), nits(k)
       else
          write(uio,'(i5,es16.6,3es15.6)') k, lam(k)/lmax, ml(k), mr(k), sr(k)
       endif
    enddo
    write(uio,'(a,i4,a,es12.4,a,i6,a,i2,a,l1)') '# selected k', ibest, '  lam/lam_max', lam_best/lmax, &
&       '  FCs kept', nsel_best, '  info', info, '  refit ', refit
 endif
 deallocate(z0, w, Q, Gf, cf, G, c, Tl, Tr, lam, el, er, mr, sr, ml)
 contains
    real(dp) function vgerr(tt)
! mean squared held-out residual per row of fold f: (t^T G_f t - 2 c_f^T t + b_f^T b_f) N / n_f
    real(dp), intent(in) :: tt(:)
    vgerr = max(dot_product(tt, matmul(Gf(:,:,f), tt)) - 2*dot_product(cf(:,f), tt) + bf(f), 0d0) &
&           *nrow/max(nf(f),1)
    end function vgerr
 end subroutine classo_cv_res
!-------------------------------------------------------------------------------
 subroutine classo_cv(nrow, ncol, nker, A, b, Kb, x0, svdcut, nfold, nlam, epsratio, tol, maxit, refit, &
&                     x, lam_best, nsel_best, info, uio)
!! The whole procedure (see classo_cv_res), after restricting x to the part of the
!! invariance kernel the data resolve: right singular vectors of A Kb with singular
!! value > svdcut * the largest (the cut svd_set makes, so the unresolved directions
!! are zero as in the SVD fit). Without it, FCs the snapshots do not see (zero or
!! 1e-20 columns of A, Ge example) carry ~zero L1 weight w_j and grow without bound
!! through the invariance couplings (Ge: phonons ~1e8 1/cm).
 integer, intent(in) :: nrow, ncol, nker, nfold, nlam, maxit, uio
 real(dp), intent(in) :: A(nrow,ncol), b(nrow), Kb(ncol,nker), x0(ncol), svdcut, epsratio, tol
 logical, intent(in) :: refit
 real(dp), intent(out) :: x(ncol), lam_best
 integer, intent(out) :: nsel_best, info
 real(dp), allocatable :: Kr(:,:)
 real(dp) eps
 integer nres, nl
 allocate(Kr(ncol,nker))
 call resolved_kernel(nrow, ncol, nker, A, Kb, svdcut, nres, Kr)
 if (uio > 0) write(uio,'(a,i6,a,i6,a,es10.2)') '# kernel', nker, '  resolved by the data', nres, &
&    '  (svdcut', svdcut, ')'
 write(*,'(a,i6,a,i6)') ' LASSO: invariance kernel', nker, ', resolved by the data', nres
 if (nres == 0) then
    x = x0 ; lam_best = 0 ; nsel_best = 0 ; info = 0 ; return
 endif
! the CV minimum at the smallest lam of the grid (info 3; exact data: always): the grid is
! extended by two decades at the same density and the CV redone, down to eps_min lam_max
 eps = epsratio ; nl = nlam
 do
    call classo_cv_res(nrow, ncol, nres, A, b, Kr(:,1:nres), x0, nfold, nl, eps, tol, maxit, refit, &
&        x, lam_best, nsel_best, info, uio)
    if (info /= 3 .or. eps <= eps_min*(1 + 1d-12)) exit
    nl = nl + nint(2*(nlam - 1)/log10(1/epsratio))
    eps = max(eps*1d-2, eps_min)
    write(*,'(a,es9.2,a,i4,a)') ' LASSO: CV minimum at the smallest lam; grid extended to', eps, &
&       ' lam_max (', nl, ' points)'
    if (uio > 0) write(uio,'(a,es9.2,a,i4,a)') '# CV minimum at the smallest lam; grid extended to', &
&       eps, ' lam_max (', nl, ' points)'
 enddo
 deallocate(Kr)
 end subroutine classo_cv
!-------------------------------------------------------------------------------
 subroutine resolved_kernel(nrow, ncol, nker, A, Kb, svdcut, nres, Kr)
!! Kr(:,1:nres) = Kb V_r, V_r the right singular vectors of A Kb with sigma > svdcut
!! sigma_max. A Kb is reduced row by row to its triangular factor R by Givens
!! rotations (A Kb is never stored; R has the singular values of A Kb, so nothing
!! is squared), then R = U S V^T by svdcmp.
 integer, intent(in) :: nrow, ncol, nker
 real(dp), intent(in) :: A(nrow,ncol), Kb(ncol,nker), svdcut
 integer, intent(out) :: nres
 real(dp), intent(out) :: Kr(ncol,nker)
 real(dp), allocatable :: R(:,:), row(:), V(:,:), sv(:), tmp(:)
 real(dp) cg, sg, h
 integer i, k, j
 allocate(R(nker,nker), row(nker), V(nker,nker), sv(nker), tmp(nker))
 R = 0
 do i = 1, nrow
    row = matmul(A(i,:), Kb)
    do k = 1, nker
       if (row(k) == 0) cycle
       h = sqrt(R(k,k)**2 + row(k)**2)
       cg = R(k,k)/h ; sg = row(k)/h
       tmp(k:nker) = cg*R(k,k:nker) + sg*row(k:nker)
       row(k:nker) = -sg*R(k,k:nker) + cg*row(k:nker)
       R(k,k:nker) = tmp(k:nker)
    enddo
 enddo
 call svdcmp(R, nker, nker, sv, V)
 nres = 0
 do j = 1, nker
    if (sv(j) > svdcut*maxval(sv)) then
       nres = nres + 1 ; Kr(:,nres) = matmul(Kb, V(:,j))
    endif
 enddo
 deallocate(R, row, V, sv, tmp)
 end subroutine resolved_kernel
!-------------------------------------------------------------------------------
 subroutine null_space_qr(n, m, Bt, tol, p, Nb)
!! orthonormal basis Nb(:,1:p) of the null space of B (m x n), given Bt = B^T (n x m,
!! overwritten): Householder QR of Bt with column pivoting; the rank r is where the
!! largest remaining column norm drops below tol (relative to the largest initial
!! one); the null space is spanned by columns r+1..n of the orthogonal factor.
 integer, intent(in) :: n, m
 real(dp), intent(inout) :: Bt(n,m)
 real(dp), intent(in) :: tol
 integer, intent(out) :: p
 real(dp), intent(out) :: Nb(n,n)
 real(dp) V(n,min(n,m)), cn(m), cmax0, al, vn, e(n), col(n)
 integer k, j, jp, r, i
 r = 0
 do j = 1, m ; cn(j) = sqrt(sum(Bt(:,j)**2)) ; enddo
 cmax0 = maxval(cn)
 if (m > 0 .and. cmax0 > 0) then
    do k = 1, min(n, m)
       do j = k, m ; cn(j) = sqrt(sum(Bt(k:n,j)**2)) ; enddo
       jp = k - 1 + maxloc(cn(k:m), 1)
       if (cn(jp) <= tol*cmax0) exit
       if (jp /= k) then
          col = Bt(:,k) ; Bt(:,k) = Bt(:,jp) ; Bt(:,jp) = col
       endif
! reflector I - 2 v v^T taking Bt(k:n,k) to -sign*|.| e_k
       V(:,k) = 0
       V(k:n,k) = Bt(k:n,k)
       al = -sign(sqrt(sum(V(k:n,k)**2)), V(k,k))
       V(k,k) = V(k,k) - al
       vn = sqrt(sum(V(k:n,k)**2))
       V(k:n,k) = V(k:n,k)/vn
       do j = k, m
          Bt(k:n,j) = Bt(k:n,j) - 2*V(k:n,k)*dot_product(V(k:n,k), Bt(k:n,j))
       enddo
       r = k
    enddo
 endif
 p = n - r
 do j = 1, p
    e = 0 ; e(r+j) = 1
    do i = r, 1, -1                      ! H_1 ... H_r e_(r+j)
       e = e - 2*V(:,i)*dot_product(V(:,i), e)
    enddo
    Nb(:,j) = e
 enddo
 end subroutine null_space_qr
!-------------------------------------------------------------------------------
 subroutine chol_fac(n, G, shift, L)
 integer, intent(in) :: n
 real(dp), intent(in) :: G(n,n), shift
 real(dp), intent(out) :: L(n,n)
 integer i, j
 real(dp) s
 L = 0
 do j = 1, n
    s = G(j,j) + shift - dot_product(L(j,1:j-1), L(j,1:j-1))
    L(j,j) = sqrt(max(s, tiny(1d0)))
    do i = j+1, n
       L(i,j) = (G(i,j) - dot_product(L(i,1:j-1), L(j,1:j-1)))/L(j,j)
    enddo
 enddo
 end subroutine chol_fac
!-------------------------------------------------------------------------------
 subroutine chol_sol(n, L, r, x)
 integer, intent(in) :: n
 real(dp), intent(in) :: L(n,n), r(n)
 real(dp), intent(out) :: x(n)
 integer i
 do i = 1, n
    x(i) = (r(i) - dot_product(L(i,1:i-1), x(1:i-1)))/L(i,i)
 enddo
 do i = n, 1, -1
    x(i) = (x(i) - dot_product(L(i+1:n,i), x(i+1:n)))/L(i,i)
 enddo
 end subroutine chol_sol

 end module lasso_constrained
