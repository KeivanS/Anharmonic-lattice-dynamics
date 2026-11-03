!============================================================
 module given_fcs
!! include_fc(rank)=2 (structure.params line 4; 2026-11-03, user): the FCs of that rank are read, not fitted.
!! Their irreducible values xfull go to the right-hand sides: force rows b - A_read x_read (set_force_displacement_matrix),
!! translational and rotational rows -C_read x_read (invariance_ranks); Huang rows contain FC2 only and are dropped
!! when FC2 is read. After the fit merge_given_fcs (invariance_ranks) puts fitted and read values into one vector, so
!! everything downstream runs as for a full fit.
!! Reading (read_fcs in ios.f90): fcN.dat (all terms) and lat_fc.dat of the earlier run. lat_fc.dat gives each atom number of
!! that run as (sublattice tau, cell n), so the terms are matched to this run's atoms by (tau, n), whatever the numbering.
!! The irreducible values of a group g are the least-squares solution of mat(t,ti) x(ti) = f(t) over its terms t
!! (f = 0 for a term absent from fcN.dat: written only above the margin); the residual is logged (consistency of the
!! read FCs with this run's symmetry analysis). Ranks 1-4; same primitive cell required (checked).
!! Layout: full_col(rnk,g,ti) = index among the kept irreducible FCs of all ranks with include_fc /= 0, in rank order.
 use constants, only : r15
 implicit none
 private
 public :: xfull, nfull, full_col, solve_small, sort_i8, find_i8

 real(r15), allocatable, save :: xfull(:)
 integer, save :: nfull = 0

 contains
!------------------------------------------------------------
 integer function full_col(rnk, g, ti)
!! index of the kept irreducible term ti of group g of rank rnk in the full layout, 0 if not a term
 use svd_stuff
 use params
 integer, intent(in) :: rnk, g, ti
 integer r
 full_col = 0
 if (include_fc(rnk) .eq. 0) return
 if (rnk .eq. 2) then
    if (map(2)%keep(counteri(2,g,ti)) .ne. 1) return
 endif
 full_col = current(rnk, g, ti)
 do r = 1, rnk-1
    if (include_fc(r) .ne. 0) full_col = full_col + map(r)%nkeptind
 enddo
 end function full_col
!------------------------------------------------------------
 subroutine solve_small(n, a, b, x)
!! a x = b, Gaussian elimination with partial pivoting (n small; a symmetric positive semi-definite normal matrix)
 integer, intent(in) :: n
 real(r15), intent(inout) :: a(n,n), b(n)
 real(r15), intent(out) :: x(n)
 integer i, k, p
 real(r15) m, tmp(n), tb, amax
 amax = maxval(abs(a))
 do k = 1, n
    p = k - 1 + maxloc(abs(a(k:n,k)),1)
    if (p .ne. k) then
       tmp = a(k,:); a(k,:) = a(p,:); a(p,:) = tmp
       tb = b(k); b(k) = b(p); b(p) = tb
    endif
    if (abs(a(k,k)) .le. 1d-12*amax) then
       a(k,k) = 1; b(k) = 0; a(k+1:n,k) = 0          ! direction not fixed by the terms: 0
       cycle
    endif
    do i = k+1, n
       m = a(i,k)/a(k,k)
       a(i,k:n) = a(i,k:n) - m*a(k,k:n); b(i) = b(i) - m*b(k)
    enddo
 enddo
 do k = n, 1, -1
    x(k) = (b(k) - dot_product(a(k,k+1:n), x(k+1:n)))/a(k,k)
 enddo
 end subroutine solve_small
!------------------------------------------------------------
 logical function le2(a, b)
 integer(8), intent(in) :: a(2), b(2)
 le2 = a(1) .lt. b(1) .or. (a(1) .eq. b(1) .and. a(2) .le. b(2))
 end function le2
!------------------------------------------------------------
 subroutine sort_i8(n, k, ord)
!! ord: permutation sorting the pairs k(:,i) ascending (merge sort)
 integer, intent(in) :: n
 integer(8), intent(in) :: k(2,n)
 integer, intent(out) :: ord(n)
 integer, allocatable :: tmp(:)
 integer w, lo, mid, hi, i, j, m
 allocate(tmp(n))
 do i = 1, n
    ord(i) = i
 enddo
 w = 1
 do while (w .lt. n)
    lo = 1
    do while (lo .le. n)
       mid = min(lo+w-1, n); hi = min(lo+2*w-1, n)
       i = lo; j = mid+1; m = lo
       do while (i .le. mid .and. j .le. hi)
          if (le2(k(:,ord(i)), k(:,ord(j)))) then
             tmp(m) = ord(i); i = i+1
          else
             tmp(m) = ord(j); j = j+1
          endif
          m = m+1
       enddo
       do while (i .le. mid)
          tmp(m) = ord(i); i = i+1; m = m+1
       enddo
       do while (j .le. hi)
          tmp(m) = ord(j); j = j+1; m = m+1
       enddo
       ord(lo:hi) = tmp(lo:hi)
       lo = lo + 2*w
    enddo
    w = 2*w
 enddo
 deallocate(tmp)
 end subroutine sort_i8
!------------------------------------------------------------
 integer function find_i8(n, k, ord, key)
!! index i with k(i) = key (binary search through ord), 0 if absent
 integer, intent(in) :: n, ord(n)
 integer(8), intent(in) :: k(2,n), key(2)
 integer lo, hi, mid
 find_i8 = 0
 lo = 1; hi = n
 do while (lo .le. hi)
    mid = (lo+hi)/2
    if (all(k(:,ord(mid)) .eq. key)) then
       find_i8 = ord(mid); return
    elseif (le2(k(:,ord(mid)), key)) then
       lo = mid+1
    else
       hi = mid-1
    endif
 enddo
 end function find_i8

 end module given_fcs
