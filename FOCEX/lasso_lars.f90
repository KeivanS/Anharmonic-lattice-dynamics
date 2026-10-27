!===============================================================================
! LASSO solution path by the LARS-lasso homotopy (Efron, Hastie, Johnstone,
! Tibshirani, Ann. Statist. 32, 407 (2004)), in Gram form:
!
!   minimise  1/2 w^T G w - c^T w + alpha ||w||_1      (G = A^T A, c = A^T y)
!
! for alpha down to the smallest value of a given grid. The path w(alpha) is
! piecewise linear; the solutions at the grid values are interpolated exactly while
! the path is followed, so only n x nlam numbers are stored. No ridge: LARS needs
! none (a column in the span of the active set never enters).
!
! Same steps as scikit-learn 1.9.1 sklearn/linear_model/_least_angle.py,
! _lars_path_solver with method='lasso' and a Gram matrix (BSD-3-Clause,
! (c) 2007-2026 The scikit-learn developers), written afresh with an index list
! instead of row/column swaps:
!  - alpha = C = max_j |c_j - (G w)_j|, all active correlations equal C;
!  - direction on the active set A: d = AA G_AA^{-1} s, AA = (s^T G_AA^{-1} s)^(-1/2);
!  - step gamma = min( C/AA, min_pos (C-cov_j)/(AA-a_j), min_pos (C+cov_j)/(AA+a_j) ),
!    a_j = G_jA d for inactive j; lasso modification: if an active coefficient
!    reaches zero first (gamma = min_pos -w_A/d), it is dropped (sign kept for the
!    rest of the step) and the Cholesky factor of G_AA is rebuilt;
!  - a column whose Cholesky pivot is below piv_tol*sqrt(G_jj) (scikit-learn uses
!    1e-7 absolute, which assumes unit-norm columns; relative here: the column's
!    distance from the span of the active set over its norm) lies in that span: it is flagged and takes no part in entry
!    or step length until a drop changes the span, so no direction with G_AA
!    singular is ever taken (scikit-learn zeroes its correlation instead);
!  - a column enters only when its |correlation| has reached the active level
!    (within lev_tol), and a step-length ratio with a zero denominator (|AA -+ a_j|
!    <= den_tol AA: the column moves with the active set) is no event. Without these
!    a not-yet-flagged duplicate of an active column ended a step on a 0/0 ratio and
!    the next column was added below the level (found on duplicated test columns);
!  - stop at the last grid value (interpolating the last segment), when alpha would
!    increase by more than inc_tol (1e-6) relative (scikit-learn's "early stopping"
!    stops at any increase; smaller rises are round-off and are clamped), or after
!    max_steps.
!===============================================================================
 module lasso_lars
 implicit none
 private
 public lars_path_gram_f
 integer, parameter :: dp = kind(1d0)
 real(dp), parameter :: piv_tol = 1d-7, tiny32 = 1.17549435d-38
 real(dp), parameter :: lev_tol = 1d-9, den_tol = 1d-10, inc_tol = 1d-6, gam_tol = 1d-12

 contains
!-------------------------------------------------------------------------------
 subroutine lars_path_gram_f(n, G, c, nlam, lam, max_steps, Wgrid, nbp, info)
!! Solutions at the decreasing grid lam(1:nlam): Wgrid(:,k) = w(lam(k)), by linear
!! interpolation between the two breakpoints of the path that bracket lam(k)
!! (exact: the path is linear in alpha between breakpoints). nbp: breakpoints passed.
!! info: 0 reached lam(nlam); 1 alpha increased (stopped early); 2 max_steps;
!! 3 no column can enter. For info /= 0 the grid points below the last breakpoint get
!! the last breakpoint's w, which is NOT their solution.
 integer, intent(in) :: n, nlam, max_steps
 real(dp), intent(in) :: G(n,n), c(n), lam(nlam)
 real(dp), intent(out) :: Wgrid(n, nlam)
 integer, intent(out) :: nbp, info
 real(dp) wprev(n), aprev
 integer kgrid
 integer act(n), sgn(n), na, k, j, jmax, it, idrop
 logical isact(n), skip(n), drop, added
 real(dp) w(n), cov(n), Lc(n,n), d(n), u(n), aj, Clev, Cj, AA, gam, g1, z, zmin, diag, s, prev_alpha

 w = 0 ; na = 0 ; isact = .false. ; drop = .false. ; nbp = 0 ; info = 2 ; skip = .false.
 prev_alpha = huge(1d0) ; kgrid = 1 ; wprev = 0 ; aprev = huge(1d0)
 do it = 0, max_steps
    cov = c - matmul(G, w)
! alpha = the common correlation of the active set (= max over all columns)
    if (na > 0) then
       Clev = maxval(abs(cov(act(1:na))))
    else
       Clev = maxval(abs(cov))
    endif
! alpha must not increase along the path; a rise up to inc_tol (relative) is round-off
! in the active correlations (ill-conditioned G_AA, many steps): clamp it and go on
    if (it > 0 .and. Clev > prev_alpha .and. Clev <= prev_alpha*(1 + inc_tol)) Clev = prev_alpha
    nbp = nbp + 1
! grid points in [Clev, aprev): interpolate between the previous breakpoint and this one
    do while (kgrid <= nlam)
       if (lam(kgrid) < Clev) exit
       if (nbp == 1) then
          Wgrid(:,kgrid) = 0                         ! above max|c| the solution is zero
       else
          s = (aprev - lam(kgrid))/(aprev - Clev)
          Wgrid(:,kgrid) = wprev + s*(w - wprev)
       endif
       kgrid = kgrid + 1
    enddo
    if (kgrid > nlam) then ; info = 0 ; return ; endif
    wprev = w ; aprev = Clev
    if (it == max_steps) info = 2
    if (it > 0 .and. Clev > prev_alpha*(1 + inc_tol)) info = 1
    if (info == 1 .or. it == max_steps) then
       do while (kgrid <= nlam) ; Wgrid(:,kgrid) = w ; kgrid = kgrid + 1 ; enddo
       return
    endif
    prev_alpha = Clev
! add the inactive column of largest correlation (not after a drop), unless it is
! degenerate with the active set; then try the next one
    added = .false.
    if (.not. drop) then
       do
          jmax = 0
          do j = 1, n
             if (isact(j) .or. skip(j)) cycle
             if (jmax == 0) then
                jmax = j
             elseif (abs(cov(j)) > abs(cov(jmax))) then
                jmax = j
             endif
          enddo
          if (jmax == 0) exit
! only a column that has reached the active level may enter (after a step that ended
! on an event of the active set itself, e.g. the end of a degenerate column's ratio)
          if (na > 0 .and. abs(cov(jmax)) < Clev*(1 - lev_tol)) then
             jmax = 0 ; exit
          endif
          call chol_row(n, na, G, act, jmax, Lc, diag)
          if (diag >= piv_tol*sqrt(G(jmax,jmax))) then
             na = na + 1 ; act(na) = jmax ; isact(jmax) = .true.
             sgn(na) = int(sign(1d0, cov(jmax))) ; Lc(na,na) = diag ; added = .true.
             exit
          endif
          skip(jmax) = .true.
       enddo
    endif
    if (na == 0) then
       info = 3 ; do while (kgrid <= nlam) ; Wgrid(:,kgrid) = w ; kgrid = kgrid + 1 ; enddo ; return
    endif
! direction: u = G_AA^{-1} s, AA = (s.u)^(-1/2), d = AA u
    do
       u(1:na) = dble(sgn(1:na))
       call chol_solve_f(n, na, Lc, u)
       AA = 1/sqrt(max(dot_product(dble(sgn(1:na)), u(1:na)), tiny(1d0)))
       d(1:na) = AA*u(1:na)
! Lasso sign condition (Efron et al. 2004, eq. 30): a column that has just entered
! must move with the sign of its correlation; otherwise its entry is a numerical tie
! (the "one at a time" condition of their Theorem 1 fails): take it out for this step
       if (.not. added .or. d(na)*sgn(na) > 0) exit
       isact(act(na)) = .false. ; skip(act(na)) = .true. ; na = na - 1 ; added = .false.
       if (na == 0) then
          info = 3 ; do while (kgrid <= nlam) ; Wgrid(:,kgrid) = w ; kgrid = kgrid + 1 ; enddo ; return
       endif
    enddo
! step length: an inactive correlation catching up, or the end of the path
    gam = Clev/AA
    do j = 1, n
       if (isact(j) .or. skip(j)) cycle           ! degenerate columns cannot enter
       aj = 0
       do k = 1, na
          aj = aj + G(j, act(k))*d(k)
       enddo
       Cj = cov(j)
! a zero denominator (to round-off) means column j moves with the active set: no event
! a zero-length event (g1 <= gam_tol C/AA) is a column already at the level, e.g. one
! just dropped (eq. 35 of Efron et al.), not one arriving there: no event
       if (abs(AA - aj) > den_tol*AA) then
          g1 = (Clev - Cj)/(AA - aj) ; if (g1 > gam_tol*Clev/AA .and. g1 < gam) gam = g1
       endif
       if (abs(AA + aj) > den_tol*AA) then
          g1 = (Clev + Cj)/(AA + aj) ; if (g1 > gam_tol*Clev/AA .and. g1 < gam) gam = g1
       endif
    enddo
! lasso modification: an active coefficient reaching zero first is dropped
    drop = .false. ; zmin = huge(1d0) ; idrop = 0
    do k = 1, na
       z = -w(act(k))/(d(k) + tiny32)
       if (z > 0 .and. z < zmin) then ; zmin = z ; idrop = k ; endif
    enddo
    if (zmin < gam) then ; gam = zmin ; drop = .true. ; endif
    do k = 1, na
       w(act(k)) = w(act(k)) + gam*d(k)
    enddo
    if (drop) then                               ! remove act(idrop), rebuild the factor
       w(act(idrop)) = 0
       isact(act(idrop)) = .false.
       do k = idrop, na-1
          act(k) = act(k+1) ; sgn(k) = sgn(k+1)
       enddo
       na = na - 1
       do k = 1, na
          call chol_row(n, k-1, G, act, act(k), Lc, diag)
          Lc(k,k) = diag
       enddo
       skip = .false.                            ! the span changed: re-test degeneracy
    endif
 enddo
 end subroutine lars_path_gram_f
!-------------------------------------------------------------------------------
 subroutine chol_row(n, na, G, act, j, Lc, diag)
!! row na+1 of the Cholesky factor of G_AA when column j joins act(1:na)
 integer, intent(in) :: n, na, act(n), j
 real(dp), intent(in) :: G(n,n)
 real(dp), intent(inout) :: Lc(n,n)
 real(dp), intent(out) :: diag
 integer k, i
 real(dp) s
 do k = 1, na                                   ! forward substitution L l = G_Aj
    s = G(act(k), j)
    do i = 1, k-1
       s = s - Lc(k,i)*Lc(na+1,i)
    enddo
    Lc(na+1,k) = s/Lc(k,k)
 enddo
 s = G(j,j) - dot_product(Lc(na+1,1:na), Lc(na+1,1:na))
 diag = sqrt(max(abs(s), epsilon(1d0)))
 end subroutine chol_row
!-------------------------------------------------------------------------------
 subroutine chol_solve_f(n, na, Lc, u)
!! u <- (L L^T)^{-1} u for the leading na x na block
 integer, intent(in) :: n, na
 real(dp), intent(in) :: Lc(n,n)
 real(dp), intent(inout) :: u(n)
 integer i, k
 do i = 1, na
    do k = 1, i-1
       u(i) = u(i) - Lc(i,k)*u(k)
    enddo
    u(i) = u(i)/Lc(i,i)
 enddo
 do i = na, 1, -1
    do k = i+1, na
       u(i) = u(i) - Lc(k,i)*u(k)
    enddo
    u(i) = u(i)/Lc(i,i)
 enddo
 end subroutine chol_solve_f

 end module lasso_lars
