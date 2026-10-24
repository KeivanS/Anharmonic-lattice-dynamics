!===========================================================================
 module invariance_ranks
!! Translational and rotational invariance rows of the fit for every rank in one loop.
!! set_translational_inv_constraints and set_rotational_inv_constraints (the default since
!! 2026-10-23) replace the rank-by-rank routines, kept as set_*_inv_constraints_old in
!! setup_amatrix.f90 (FOCEX_INVARIANCE=old selects them, see focex.f90).
!!
!! The first index of a force constant is an atom i1 of the primitive cell, the others are its
!! neighbours. For rank n, with all indices but the last fixed (the "prefix"):
!!  translational  sum_{i_n} Phi^n(i1 a1, ..., i_n a_n) = 0
!!                 row key (i1..i_{n-1} ; a1..a_n)
!!  rotational     sum_{i_n a_n mu} Phi^n(i1 a1, ..., i_n a_n) X^mu_{i_n} eps(a_n,mu,nu)
!!               + sum_{k=1..n-1} sum_g Phi^{n-1}(i1 a1, .., i_k g, .., i_{n-1} a_{n-1}) eps(g,a_k,nu) = 0
!!                 row key (i1..i_{n-1} ; a1..a_{n-1} ; nu)
!! As in the old routines, only prefixes whose atoms i2..i_{n-1} lie within nshells(n,i1) are
!! used: beyond the range of rank n the rotational rows would force combinations of rank n-1 to
!! vanish, an artefact of the truncation. Every direction combination is used (the old
!! translational rows had a2>=a1, a3>=a2); repeated and proportional rows are removed by
!! compare2previous_lines, as before. Huang invariance (rank 2 only) is unchanged.
!!
!! One pass over the terms of every group adds (key, column, value) contributions; sorting the
!! keys gathers each row. Column of independent term ti of group g of rank n: column_of.
 use constants, only : r15
 implicit none
 private
 public set_translational_inv_constraints, set_rotational_inv_constraints

 integer, allocatable :: ckey(:,:), ccol(:)
 real(r15), allocatable :: cval(:)
 integer ncon, keylen

 contains
!---------------------------------------------------------------------------
 subroutine set_translational_inv_constraints
!! outputs atransl, btransl, transl_constraints
 call build_inv_rows(.true., .false.)
 end subroutine set_translational_inv_constraints
!---------------------------------------------------------------------------
 subroutine set_rotational_inv_constraints
!! outputs arot, brot, rot_constraints
 call build_inv_rows(.false., .true.)
 end subroutine set_rotational_inv_constraints
!---------------------------------------------------------------------------
 subroutine build_inv_rows(do_tr, do_rot)
 use svd_stuff
 use params
 use atoms_force_constants
 use ios
 implicit none
 logical, intent(in) :: do_tr, do_rot
 real(r15), allocatable :: arow(:,:), rrow(:,:)
 integer ntr, nro, rnk, lc(3,3,3), g, t, ti, i1, c, nu, mu, k, ak, l, a
 integer, allocatable :: key(:)
 real(r15) cf

 lc = 0
 lc(1,2,3) =  1 ; lc(2,3,1) =  1 ; lc(3,1,2) =  1
 lc(1,3,2) = -1 ; lc(2,1,3) = -1 ; lc(3,2,1) = -1
 write(ulog,*)' BUILD_INV_ROWS (all ranks): translational, rotational =', do_tr, do_rot
 call check_columns
 ntr = 0 ; nro = 0
 allocate(arow(16,nindepfc), rrow(16,nindepfc))

 do rnk = 1, maxrank
    if (include_fc(rnk) .eq. 0) cycle
    keylen = 2*rnk - 1
    allocate(key(keylen))
!------------------------------------ translational, rank rnk
    if (do_tr) then
    call start_contributions
    do g = 1, map(rnk)%ngr
    do t = 1, map(rnk)%nt(g)
       if (.not. prefix_ok(rnk, map(rnk)%gr(g)%iat(:,t), rnk)) cycle
       key(1:rnk-1) = map(rnk)%gr(g)%iat(1:rnk-1,t)
       key(rnk:2*rnk-1) = map(rnk)%gr(g)%ixyz(1:rnk,t)
       do ti = 1, map(rnk)%ntind(g)
          c = column_of(rnk, g, ti)
          if (c .gt. 0) call add(key, c, map(rnk)%gr(g)%mat(t,ti))
       enddo
    enddo
    enddo
    call gather_rows(arow, ntr, 'translational', rnk)
    endif
!------------------------------------ rotational, rank rnk: the rank-rnk part ...
    if (do_rot) then
    call start_contributions
    do g = 1, map(rnk)%ngr
    do t = 1, map(rnk)%nt(g)
       if (.not. prefix_ok(rnk, map(rnk)%gr(g)%iat(:,t), rnk)) cycle
       l = map(rnk)%gr(g)%iat(rnk,t) ; a = map(rnk)%gr(g)%ixyz(rnk,t)
       key(1:rnk-1) = map(rnk)%gr(g)%iat(1:rnk-1,t)
       key(rnk:2*rnk-2) = map(rnk)%gr(g)%ixyz(1:rnk-1,t)
       do nu = 1, 3
          cf = 0
          do mu = 1, 3
             cf = cf + atompos(mu,l)*lc(a,mu,nu)
          enddo
          if (cf .eq. 0) cycle
          key(2*rnk-1) = nu
          do ti = 1, map(rnk)%ntind(g)
             c = column_of(rnk, g, ti)
             if (c .gt. 0) call add(key, c, map(rnk)%gr(g)%mat(t,ti)*cf)
          enddo
       enddo
    enddo
    enddo
!------------------------------------ ... and the rank rnk-1 part: slot k of the prefix takes
!                                        direction ak, the term has g = ixyz(k); eps(g,ak,nu)
    if (rnk .ge. 2) then
    if (include_fc(rnk-1) .ne. 0) then
       do g = 1, map(rnk-1)%ngr
       do t = 1, map(rnk-1)%nt(g)
          if (.not. prefix_ok(rnk-1, map(rnk-1)%gr(g)%iat(:,t), rnk)) cycle
          do k = 1, rnk-1
          do ak = 1, 3
          do nu = 1, 3
             if (lc(map(rnk-1)%gr(g)%ixyz(k,t), ak, nu) .eq. 0) cycle
             key(1:rnk-1) = map(rnk-1)%gr(g)%iat(1:rnk-1,t)
             key(rnk:2*rnk-2) = map(rnk-1)%gr(g)%ixyz(1:rnk-1,t)
             key(rnk-1+k) = ak
             key(2*rnk-1) = nu
             do ti = 1, map(rnk-1)%ntind(g)
                c = column_of(rnk-1, g, ti)
                if (c .gt. 0) call add(key, c, map(rnk-1)%gr(g)%mat(t,ti)*lc(map(rnk-1)%gr(g)%ixyz(k,t),ak,nu))
             enddo
          enddo
          enddo
          enddo
       enddo
       enddo
    endif
    endif
    call gather_rows(rrow, nro, 'rotational', rnk)
    endif
    deallocate(key)
 enddo

 ! include_fc=1 ranks only for now: the right-hand sides are zero. For a rank read from file
! (include_fc=2) its terms are to go to these right-hand sides, b = -sum(coef * known FC)
 if (do_tr) then
    transl_constraints = ntr
    if (allocated(atransl)) deallocate(atransl)
    if (allocated(btransl)) deallocate(btransl)
    allocate(atransl(ntr,nindepfc), btransl(ntr))
    atransl = arow(1:ntr,:) ; btransl = 0
    write(ulog,*)' SET_TRANSLATIONAL_INV_CONSTRAINTS: rows (all ranks) =', ntr
 endif
 if (do_rot) then
    rot_constraints = nro
    if (allocated(arot)) deallocate(arot)
    if (allocated(brot)) deallocate(brot)
    allocate(arot(nro,nindepfc), brot(nro))
    arot = rrow(1:nro,:) ; brot = 0
    write(ulog,*)' SET_ROTATIONAL_INV_CONSTRAINTS: rows (all ranks) =', nro
 endif
 deallocate(arow, rrow)
 if (allocated(ckey)) deallocate(ckey, ccol, cval)

 contains

   logical function prefix_ok(r, iat, rshell)
!! the first atom is in the primitive cell; atoms 2..r of a rank-r term that form a prefix of a
!! rank-rshell row lie within nshells(rshell, first atom) (rshell = r for the rank-r rows, r+1
!! for the rank-r part of the rank-(r+1) rotational rows)
   integer, intent(in) :: r, iat(:), rshell
   integer m2, last
   prefix_ok = .false.
   i1 = iat(1)
   if (i1 .lt. 1 .or. i1 .gt. natom_prim_cell) return
   if (i1 .ne. iatomcell0(i1)) return
   last = min(r, rshell-1)          ! atoms of the prefix: 1..rshell-1
   do m2 = 2, last
      if (iatomneighbor(i1, iat(m2)) .gt. nshells(rshell, i1)) return
   enddo
   prefix_ok = .true.
   end function prefix_ok

   subroutine gather_rows(rows, nrow, label, r)
!! sorts the contributions by key, adds each key's contributions into one row and keeps it
!! if it is not zero and not a multiple of a row kept already (any rank)
   real(r15), allocatable, intent(inout) :: rows(:,:)
   integer, intent(inout) :: nrow
   character(*), intent(in) :: label
   integer, intent(in) :: r
   integer, allocatable :: ord(:)
   real(r15), allocatable :: a1d(:), zero(:), tmp(:,:)
   integer i, j, n0
   logical new
   allocate(ord(ncon), a1d(nindepfc), zero(nindepfc)) ; zero = 0
   call sort_keys(ord)
   n0 = nrow
   i = 1
   do while (i .le. ncon)
      a1d = 0
      j = i
      do while (j .le. ncon)
         if (any(ckey(:,ord(j)) .ne. ckey(:,ord(i)))) exit
         a1d(ccol(ord(j))) = a1d(ccol(ord(j))) + cval(ord(j))
         j = j + 1
      enddo
      i = j
      if (a1d .myeq. zero) cycle
      call compare2previous_lines(size(rows,1), nindepfc, a1d, rows, nrow, new)
      if (.not. new) cycle
      if (nrow .eq. size(rows,1)) then              ! grow the row store
         allocate(tmp(2*size(rows,1), nindepfc)) ; tmp = 0
         tmp(1:nrow,:) = rows(1:nrow,:)
         call move_alloc(tmp, rows)
      endif
      nrow = nrow + 1
      rows(nrow,:) = a1d
   enddo
   write(ulog,'(a,a,a,i2,a,i8,a,i8)')'  ', label, ' rank', r, ': contributions', ncon, ', new rows', nrow-n0
   deallocate(ord, a1d, zero)
   end subroutine gather_rows

 end subroutine build_inv_rows
!---------------------------------------------------------------------------
 integer function column_of(rnk, g, ti)
!! column of independent term ti of group g of rank rnk in the fit, 0 if it is not a column:
!! the rank is not in the model, or an FC2 group not kept (beyond the supercell).
!! include_fc=2 (FCs read from file) is not treated yet: their terms will go to the right-hand
!! side of the rows instead of a column.
 use svd_stuff
 use params
 implicit none
 integer, intent(in) :: rnk, g, ti
 integer r
 column_of = 0
 if (include_fc(rnk) .ne. 1) return
 if (rnk .eq. 2) then
    if (map(2)%keep(counteri(2,g,ti)) .ne. 1) return
 endif
 column_of = current(rnk, g, ti)
 do r = 1, rnk-1
    if (include_fc(r) .eq. 1) column_of = column_of + map(r)%nkeptind
 enddo
 end function column_of
!---------------------------------------------------------------------------
 subroutine check_columns
!! the columns of column_of must fill 1..nindepfc
 use svd_stuff
 use params
 use ios
 implicit none
 integer r, n
 n = 0
 do r = 1, maxrank
    if (include_fc(r) .eq. 1) n = n + map(r)%nkeptind
 enddo
 if (n .ne. nindepfc) then
    write(ulog,*)'SET_INVARIANCE_CONSTRAINTS_RANKS: columns of included ranks', n, ' .ne. nindepfc', nindepfc
    write(*,*)'SET_INVARIANCE_CONSTRAINTS_RANKS: columns of included ranks', n, ' .ne. nindepfc', nindepfc
    stop
 endif
 end subroutine check_columns
!---------------------------------------------------------------------------
 subroutine start_contributions
 implicit none
 if (allocated(ckey)) deallocate(ckey, ccol, cval)
 allocate(ckey(keylen,1024), ccol(1024), cval(1024))
 ncon = 0
 end subroutine start_contributions
!---------------------------------------------------------------------------
 subroutine add(key, c, v)
 implicit none
 integer, intent(in) :: key(:), c
 real(r15), intent(in) :: v
 integer, allocatable :: k2(:,:), c2(:)
 real(r15), allocatable :: v2(:)
 integer n
 if (v .eq. 0) return
 n = size(ccol)
 if (ncon .eq. n) then
    allocate(k2(keylen,2*n), c2(2*n), v2(2*n))
    k2(:,1:n) = ckey ; c2(1:n) = ccol ; v2(1:n) = cval
    call move_alloc(k2, ckey) ; call move_alloc(c2, ccol) ; call move_alloc(v2, cval)
 endif
 ncon = ncon + 1
 ckey(:,ncon) = key ; ccol(ncon) = c ; cval(ncon) = v
 end subroutine add
!---------------------------------------------------------------------------
 subroutine sort_keys(ord)
!! ord: permutation putting the keys in lexicographic order (merge sort, stable)
 implicit none
 integer, intent(out) :: ord(:)
 integer, allocatable :: w(:)
 integer i, width, lo, mid, hi
 do i = 1, ncon
    ord(i) = i
 enddo
 allocate(w(ncon))
 width = 1
 do while (width .lt. ncon)
    lo = 1
    do while (lo .le. ncon)
       mid = min(lo+width-1, ncon) ; hi = min(lo+2*width-1, ncon)
       if (mid .lt. hi) call merge_runs(lo, mid, hi)
       lo = lo + 2*width
    enddo
    width = 2*width
 enddo
 deallocate(w)
 contains
   subroutine merge_runs(lo, mid, hi)
   integer, intent(in) :: lo, mid, hi
   integer i, j, k
   i = lo ; j = mid+1 ; k = lo
   do while (i .le. mid .and. j .le. hi)
      if (key_le(ord(i), ord(j))) then
         w(k) = ord(i) ; i = i+1
      else
         w(k) = ord(j) ; j = j+1
      endif
      k = k+1
   enddo
   do while (i .le. mid) ; w(k) = ord(i) ; i = i+1 ; k = k+1 ; enddo
   do while (j .le. hi)  ; w(k) = ord(j) ; j = j+1 ; k = k+1 ; enddo
   ord(lo:hi) = w(lo:hi)
   end subroutine merge_runs
 end subroutine sort_keys
!---------------------------------------------------------------------------
 logical function key_le(i, j)
 implicit none
 integer, intent(in) :: i, j
 integer m
 do m = 1, keylen
    if (ckey(m,i) .lt. ckey(m,j)) then
       key_le = .true. ; return
    elseif (ckey(m,i) .gt. ckey(m,j)) then
       key_le = .false. ; return
    endif
 enddo
 key_le = .true.
 end function key_le

 end module invariance_ranks
