!******************************************************************************
! Fit diagnostics (user, 2026-11-11): statistics of the force fit and what they say about over-fitting
! and about the completeness of the model, written to the log, to the screen (conclusions only) and at
! the end of the svd file of the fit (svd-kernel.dat, svd-lasso.dat, svd-all.dat).
!
! svd_set calls fit_stats when diag_on is set by the caller (the force fits only). With U, the filter
! factors phi_k = w_k^3/sqrt(w_k^6+wcut^6) of svbksb and the residual r = A x - b:
!   hat matrix H = U diag(phi) U^T, leverages h_i = sum_k phi_k U_ik^2, p_eff = trace H = sum phi_k;
!   leave-one-row-out residual  r_i/(1-h_i)                          (PRESS of the rows);
!   leave-one-snapshot-out residual (I - H_bb)^-1 r_b for the rows b of one snapshot (the snapshot is
!   left out of the fit and predicted by the others: the held-out force error).
! For a correct model and independent noise, E[PRESS]/E[RSS] = 1/(1-p/m)^2 (uniform leverage), so the
! held-out error exceeds the training error by the factor 1/(1-p_eff/m) at most; a larger factor means
! over-fitting (or correlated rows), a small factor with a large residual means the residual is not
! variance but bias: the model is incomplete (ranges, ranks) or the forces are noisy.
! Thresholds: p/m > pn_warn (0.2, fewer than 5 force rows per free parameter: +25 % expected held-out
! error) and > pn_strong (0.5: +100 %) warn that LASSO (enforce_inv=3) may be / is needed.
!******************************************************************************
 module fit_diag
 use constants, only : r15
 implicit none
 private
 public :: diag_on, diag_label, d_refit, add_snapshot_rows, fit_stats, fit_report, fit_diag_reset, set_basic

 integer, allocatable, save :: snap_rows(:)          ! rows of each snapshot, in the order of aforce
 integer, save :: nsnap = 0
 logical, save :: diag_on = .false.                  ! set by the caller of a force fit
 character(len=60), save :: diag_label = ''
 integer, save :: d_refit = -1                       ! svd_set noise refit: -1 not tried, 0 rejected, 1 kept
 logical, save :: d_done = .false., d_snap_done = .false.
 integer, save :: d_m = 0, d_p = 0, d_nlev1 = 0
 logical, save :: d_pbound = .false.               ! p counts all FCs (invariances not eliminated)
 integer, save :: d_snap_why = 0                     ! 1 no snapshot list, 2 rows differ, 3 too costly, 4 singular
 real(r15), save :: d_peff = 0, d_rss = 0, d_bb = 0, d_press_row = 0, d_press_snap = 0, d_rss_snap = 0
 integer, save :: d_nsing = 0                        ! snapshots that alone fix some parameters (left out)
 real(r15), parameter :: pn_warn = 0.2_r15, pn_strong = 0.5_r15, cv_warn = 1.25_r15, eps_model = 0.05_r15
 real(r15), parameter :: snap_cost_max = 5e11_r15    ! flops above which the snapshot LOO is skipped

 contains
!------------------------------------------------------------------------------
 subroutine add_snapshot_rows(ncfg, nrows)
!! ncfg more snapshots of nrows force rows each, appended to aforce in that order
 integer, intent(in) :: ncfg, nrows
 integer, allocatable :: tmp(:)
 integer i
 if (.not. allocated(snap_rows)) allocate(snap_rows(0))
 allocate(tmp(nsnap + ncfg))
 tmp(1:nsnap) = snap_rows(1:nsnap)
 do i = 1, ncfg
    tmp(nsnap+i) = nrows
 enddo
 call move_alloc(tmp, snap_rows)
 nsnap = nsnap + ncfg
 end subroutine add_snapshot_rows
!------------------------------------------------------------------------------
 subroutine fit_diag_reset
 d_done = .false.; d_snap_done = .false.; d_refit = -1; d_pbound = .false.
 d_m = 0; d_p = 0; d_peff = 0; d_rss = 0; d_bb = 0; d_press_row = 0; d_press_snap = 0; d_nlev1 = 0
 d_rss_snap = 0; d_nsing = 0; d_snap_why = 0
 end subroutine fit_diag_reset
!------------------------------------------------------------------------------
 subroutine set_basic(m, p)
!! only the sizes (fits without the hat matrix: enforce_inv=0 with the invariance rows in the SVD); p is then
!! the number of FCs, an upper bound of the free parameters (the invariance rows are not eliminated)
 integer, intent(in) :: m, p
 d_m = m; d_p = p; d_pbound = .true.
 end subroutine set_basic
!------------------------------------------------------------------------------
 subroutine fit_stats(m, n, u, phi, r, b)
!! u(m,n): the U of the SVD of the fitted matrix; phi(n): filter factors; r = A x - b; b: right-hand side
 integer, intent(in) :: m, n
 real(r15), intent(in) :: u(m,n), phi(n), r(m), b(m)
 real(r15), allocatable :: h(:), pb(:), us(:,:), hb(:,:), e(:)
 integer i, k, ib, i1, i2, nb
 real(r15) cost
 logical sing
 real(r15), allocatable :: flag(:), rb(:)

 d_done = .true.
 d_m = m; d_p = n; d_peff = sum(phi); d_rss = sum(r**2); d_bb = sum(b**2)
 allocate(h(m))
 do i = 1, m
    h(i) = sum(phi*u(i,:)**2)
 enddo
 d_press_row = 0; d_nlev1 = 0
 do i = 1, m
    if (h(i) .gt. 1 - 1d-8) then
       d_nlev1 = d_nlev1 + 1                  ! row fitted exactly by its own parameters: not predicted
       cycle
    endif
    d_press_row = d_press_row + (r(i)/(1 - h(i)))**2
 enddo
 deallocate(h)

 d_snap_done = .false.; d_snap_why = 1
 if (.not. allocated(snap_rows)) return
 d_snap_why = 2
 if (sum(snap_rows) .ne. m) return              ! not the force rows alone (enforce_inv=0: invariance rows too)
 d_snap_why = 3
 cost = real(n,r15)*sum(real(snap_rows,r15)**2)
 if (cost .gt. snap_cost_max) return
 d_snap_why = 0
 allocate(pb(nsnap), flag(nsnap), rb(nsnap))
!$omp parallel do schedule(dynamic) private(ib,i1,i2,nb,us,hb,e,k,sing)
 do ib = 1, nsnap
    i1 = sum(snap_rows(1:ib-1)) + 1; nb = snap_rows(ib); i2 = i1 + nb - 1
    allocate(us(nb,n), hb(nb,nb), e(nb))
    do k = 1, n
       us(:,k) = u(i1:i2,k)*phi(k)
    enddo
    hb = -matmul(us, transpose(u(i1:i2,:)))
    do k = 1, nb
       hb(k,k) = hb(k,k) + 1
    enddo
    e = r(i1:i2)
    rb(ib) = sum(e**2)
    call solve_gauss(nb, hb, e, sing)
    if (sing) then
       pb(ib) = 0; flag(ib) = 1
    else
       pb(ib) = sum(e**2); flag(ib) = 0
    endif
    deallocate(us, hb, e)
 enddo
!$omp end parallel do
 d_press_snap = 0; d_rss_snap = 0
 do ib = 1, nsnap                               ! in snapshot order: independent of the threads
    if (flag(ib) .ne. 0) cycle
    d_press_snap = d_press_snap + pb(ib)
    d_rss_snap = d_rss_snap + rb(ib)
 enddo
 d_nsing = nint(sum(flag))
! a snapshot whose rows alone fix some parameters cannot be predicted by the others: left out; the test
! needs most snapshots
 d_snap_done = (2*d_nsing .le. nsnap)
 if (.not. d_snap_done) d_snap_why = 4
 deallocate(pb, flag, rb)
 end subroutine fit_stats
!------------------------------------------------------------------------------
 subroutine solve_gauss(n, a, x, sing)
!! a x = b (b in x on entry), Gaussian elimination with partial pivoting; a destroyed
 integer, intent(in) :: n
 real(r15), intent(inout) :: a(n,n), x(n)
 logical, intent(out) :: sing
 integer i, k, p
 real(r15) amax, f, t, row(n)
 sing = .false.
 amax = maxval(abs(a))
 do k = 1, n
    p = k - 1 + maxloc(abs(a(k:n,k)),1)
    if (abs(a(p,k)) .le. 1d-10*amax) then
       sing = .true.; return
    endif
    if (p .ne. k) then
       row = a(k,:); a(k,:) = a(p,:); a(p,:) = row
       t = x(k); x(k) = x(p); x(p) = t
    endif
    do i = k+1, n
       f = a(i,k)/a(k,k)
       if (f .ne. 0) then
          a(i,k:n) = a(i,k:n) - f*a(k,k:n)
          x(i) = x(i) - f*x(k)
       endif
    enddo
 enddo
 do k = n, 1, -1
    x(k) = (x(k) - dot_product(a(k,k+1:n), x(k+1:n)))/a(k,k)
 enddo
 end subroutine solve_gauss
!------------------------------------------------------------------------------
 subroutine fit_report(svdfile, nx, x, sig, nrank_cols, eps_fit)
!! the block of statistics and conclusions, to the log, the screen (conclusions) and the end of svdfile.
!! x, sig: fitted irreducible FCs and their standard errors; nrank_cols(rank): their number per rank, in
!! column order; eps_fit: relative residual ||F_dft-F_fit||/||F_dft|| of the final FCs (percent)
 use ios, only : ulog
 character(len=*), intent(in) :: svdfile
 integer, intent(in) :: nx, nrank_cols(8)
 real(r15), intent(in) :: x(nx), sig(nx), eps_fit
 integer, parameter :: maxl = 60
 character(len=200) :: lines(maxl), concl(30)
 integer nl, nc, k, r, c0, nsmall(8), nsig(8), us, ios_, lastrank
 real(r15) pn, peffn, eps, eps_row, eps_snap, ratio, ratio0, frac(8)
 logical over

 nl = 0; nc = 0
 call addl(' ==================== FIT DIAGNOSTICS ('//trim(diag_label)//') ====================')
 if (d_m .le. 0) then
    call addl(' no force rows: nothing to report')
    call flush_out
    return
 endif
 pn = real(d_p,r15)/d_m
 write(lines(nl+1),'(a,i10,a,i8,a)') ' force rows m                                   =',d_m,'   (',nsnap,' snapshots)'
 nl = nl + 1
 if (d_pbound) then
    write(lines(nl+1),'(a,i10)') ' FCs p (invariances not eliminated: upper bound)=',d_p; nl = nl + 1
 else
    write(lines(nl+1),'(a,i10)') ' free parameters p (after the invariances)      =',d_p; nl = nl + 1
 endif
 write(lines(nl+1),'(a,f10.4,a,f6.2,a)') ' p/m                                            =',pn, &
 &     '   (',1/max(pn,1d-30),' rows per parameter)'
 nl = nl + 1
 eps = -1; eps_row = -1; eps_snap = -1; ratio = -1; ratio0 = -1
 if (d_done) then
    peffn = d_peff/d_m
    write(lines(nl+1),'(a,f10.2,a,f8.4)') ' effective parameters p_eff = trace(hat matrix) =',d_peff,'   p_eff/m =',peffn
    nl = nl + 1
! ratios to the training residual of this SVD (its right-hand side may be shifted: b - A xpart, or the
! residual of a LASSO fit), applied to the relative residual eps_fit of the final FCs
    eps = eps_fit
    ratio0 = 1/max(1 - peffn, 1d-6)
    eps_row = eps*sqrt(d_press_row/max(d_rss,tiny(1d0)))
    write(lines(nl+1),'(a,f10.3,a)') ' training residual ||F_dft-F_fit||/||F_dft||    =',eps,' %'; nl = nl + 1
    write(lines(nl+1),'(a,f10.3,a)') ' leave-one-row-out residual                     =',eps_row,' %'; nl = nl + 1
    if (d_nlev1 .gt. 0) then
       write(lines(nl+1),'(a,i8,a)') '   (',d_nlev1,' rows with leverage 1 left out: fitted by their own parameters)'
       nl = nl + 1
    endif
    if (d_snap_done) then
       ratio = sqrt(d_press_snap/max(d_rss_snap,tiny(1d0)))
       eps_snap = eps*ratio
       write(lines(nl+1),'(a,f10.3,a)') ' leave-one-snapshot-out (held-out) residual     =',eps_snap,' %'; nl = nl + 1
       if (d_nsing .gt. 0) then
          write(lines(nl+1),'(a,i6,a)') '   (',d_nsing,' snapshots fix some parameters alone: left out of this test)'
          nl = nl + 1
       endif
       write(lines(nl+1),'(a,f10.3,a,f6.3,a)') ' held-out / training                            =',ratio, &
 &          '   (correct model, independent noise: about ',ratio0,')'
       nl = nl + 1
    else
       ratio = sqrt(d_press_row/max(d_rss,tiny(1d0)))
       write(lines(nl+1),'(a,f10.3,a,f6.3,a)') ' leave-one-row-out / training                   =',ratio, &
 &          '   (correct model, independent noise: about ',ratio0,')'
       nl = nl + 1
       select case (d_snap_why)
       case (1) ; call addl(' (leave-one-snapshot-out not computed: no snapshot list)')
       case (2) ; call addl(' (leave-one-snapshot-out not computed: the SVD rows are not the force rows alone)')
       case (3) ; call addl(' (leave-one-snapshot-out not computed: too costly, n*sum(rows^2) > 5e11)')
       case (4) ; call addl(' (leave-one-snapshot-out not computed: most snapshots fix some parameters alone)')
       end select
    endif
 endif
 if (.not. d_done) then
    write(lines(nl+1),'(a,f10.3,a)') ' relative residual ||F_dft-F_fit||/||F_dft||    =',eps_fit,' %'; nl = nl + 1
 endif

! per rank: FCs smaller than twice their standard error
 call addl(' rank   fitted FCs   |x| < 2 sigma   fraction')
 c0 = 0; nsmall = 0; nsig = 0; frac = -1; lastrank = 0
 do r = 1, 8
    if (nrank_cols(r) .le. 0) cycle
    lastrank = r
    do k = c0+1, min(c0+nrank_cols(r), nx)
       if (sig(k) .gt. 0) then
          nsig(r) = nsig(r) + 1
          if (abs(x(k)) .lt. 2*sig(k)) nsmall(r) = nsmall(r) + 1
       endif
    enddo
    if (nsig(r) .gt. 0) frac(r) = real(nsmall(r),r15)/nsig(r)
    if (nsig(r) .gt. 0) then
       write(lines(nl+1),'(i5,i13,i16,f11.3)') r, nrank_cols(r), nsmall(r), frac(r)
    else
       write(lines(nl+1),'(i5,i13,a)') r, nrank_cols(r), '      (no standard errors)'
    endif
    nl = nl + 1
    c0 = c0 + nrank_cols(r)
 enddo

! ---- conclusions
 if (pn .gt. pn_strong) then
    call addc('WARNING: p/m > 0.5 (fewer than 2 force rows per free parameter): the least-squares fit is')
    call addc('  over-parametrized; LASSO (enforce_inv=3) is strongly recommended, or shorter ranges / more snapshots.')
 elseif (pn .gt. pn_warn) then
    call addc('WARNING: p/m > 0.2 (fewer than 5 force rows per free parameter): risk of over-fitting;')
    call addc('  LASSO (enforce_inv=3) may be necessary (or shorter ranges / more snapshots).')
 endif
 over = .false.
 if (ratio .gt. 0) then
    if (ratio .gt. cv_warn) then
       over = .true.
       write(concl(nc+1),'(a,f6.1,a)') 'OVER-FITTING: forces of held-out data are fitted ',100*(ratio-1), &
 &          ' % worse than the training data.'
       nc = nc + 1
       if (ratio .gt. 1.1*ratio0) then
          call addc('  More than independent noise explains: parameters are fitted to individual snapshots')
          call addc('  (snapshots too similar, or too few). Use LASSO (enforce_inv=3), shorter ranges or more snapshots.')
       else
          call addc('  As expected from p_eff/m for this many parameters: use LASSO (enforce_inv=3) or more snapshots.')
       endif
    else
       call addc('Not over-fitted: the held-out error is close to the training error.')
    endif
 endif
 if (.not. over .and. eps_fit .gt. 100*eps_model) then
    write(concl(nc+1),'(a,f7.2,a)') 'MODEL INCOMPLETE (or noisy forces): the residual of ',eps_fit, &
 &       ' % is systematic, not over-fitting.'
    nc = nc + 1
    call addc('  Missing terms (FC ranges too short, a missing rank such as FC4) or DFT force noise of this size.')
    call addc('  Longer ranges / more ranks lower training AND held-out residuals if it is the model; a check of')
    call addc('  the DFT force precision (EDIFF, ENCUT, k-mesh, FFT grid) tells the noise apart.')
    if (d_refit .eq. 0) call addc('  The SVD noise test rejected a lower-rank refit: the misfit behaves as model error.')
 endif
 do r = 2, 8          ! not rank 1: FC1 of a relaxed reference is zero within its error by design
    if (frac(r) .lt. 0) cycle
    if (frac(r) .gt. 0.5) then
       write(concl(nc+1),'(a,i1,a,f5.1,a)') 'Rank ',r,': ',100*frac(r), &
 &          ' % of its FCs are below 2 sigma: the data do not resolve most of them;'
       nc = nc + 1
       call addc('  its range is longer than the data support (LASSO, enforce_inv=3, would remove them).')
    elseif (r .eq. lastrank .and. r .ge. 3 .and. frac(r) .lt. 0.25 .and. .not.over .and. eps_fit .gt. 100*eps_model) then
       write(concl(nc+1),'(a,i1,a,f5.1,a)') 'Rank ',r,': only ',100*frac(r), &
 &          ' % of its FCs are below 2 sigma: the data support all of them; its range (or the next rank) can be extended.'
       nc = nc + 1
    endif
 enddo
 if (nc .eq. 0) call addc('No warnings.')
 call addl(' Conclusions:')
 do k = 1, nc
    call addl('  '//trim(concl(k)))
 enddo
 call addl(' ==========================================================================')
 call flush_out

 contains
   subroutine addl(s)
   character(len=*), intent(in) :: s
   if (nl .lt. maxl) then
      nl = nl + 1; lines(nl) = s
   endif
   end subroutine addl
   subroutine addc(s)
   character(len=*), intent(in) :: s
   if (nc .lt. size(concl)) then
      nc = nc + 1; concl(nc) = s
   endif
   end subroutine addc
   subroutine flush_out
   integer i
   do i = 1, nl
      write(ulog,'(a)') trim(lines(i))
   enddo
   write(*,'(a)') trim(lines(1))
   do i = 1, nc
      write(*,'(a)') '  '//trim(concl(i))
   enddo
   if (len_trim(svdfile) .gt. 0) then
      open(newunit=us, file=svdfile, status='unknown', position='append', iostat=ios_)
      if (ios_ .eq. 0) then
         do i = 1, nl
            write(us,'(a)') trim(lines(i))
         enddo
         close(us)
      endif
   endif
   end subroutine flush_out
 end subroutine fit_report
 end module fit_diag
