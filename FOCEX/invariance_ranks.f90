!===========================================================================
 module invariance_ranks
!! Translational and rotational invariance rows of the fit for every rank in one loop.
!! set_translational_inv_constraints and set_rotational_inv_constraints (the default since
!! 2026-10-23) replaced the rank-by-rank routines set_*_inv_constraints_old (removed 2026-11-01 after the
!! exact-model and DFT checks recorded in CLAUDE.md; in git history before that date).
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
!! born_flag 2 (2026-11-03, user): the fit sees FC2 - Phi_B (Born-dipole FCs subtracted with the forces,
!! no FC1 subtracted), so the rank-2 rotational rows get the right-hand side -sum_{a,mu} M(i1;a1,a,mu) eps(a,mu,nu),
!! M = sum_j Phi_B(i1 a1, j a) (X_j - X_i1)_mu from coulomb_lr. Rank-3 rows (i1 a1, i2 a2; nu): b = -A_full x_sub with
!! x_sub's FC3 = that of the pair model tr(Z_i Z_j^T)/3 phi_eps (lr_rot3_pair, lr_rot3_onsite).
 use constants, only : r15
 implicit none
 private
 public set_translational_inv_constraints, set_rotational_inv_constraints, write_map_dump, column_of, merge_given_fcs

 integer, allocatable :: ckey(:,:), ccol(:)
 real(r15), allocatable :: cval(:)
 integer ncon, keylen
 integer, allocatable :: rkey(:,:)     ! rotational rows kept: rank and key (i1, a1, nu for rank 2)
 real(r15) :: rdrop                    ! max |C_read x_read| of rows without fitted columns (read FCs' own violation)

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
 use born, only : born_flag
 use coulomb_lr, only : lr_moments, lr_rot3_pair, lr_rot3_onsite
 use given_fcs, only : xfull, full_col
 implicit none
 logical, intent(in) :: do_tr, do_rot
 real(r15), allocatable :: arow(:,:), rrow(:,:), mo(:,:,:,:), trhs(:), rrhs(:)
 real(r15) b3(3,3,3), b2max, b3max
 integer i2, a2, count_log, count_log2
 integer ntr, nro, rnk, lc(3,3,3), g, t, ti, i1, c, nu, mu, k, ak, l, a, ir, a1
 integer, allocatable :: key(:)
 real(r15) cf

 lc = 0
 lc(1,2,3) =  1 ; lc(2,3,1) =  1 ; lc(3,1,2) =  1
 lc(1,3,2) = -1 ; lc(2,1,3) = -1 ; lc(3,2,1) = -1
 write(ulog,*)' BUILD_INV_ROWS (all ranks): translational, rotational =', do_tr, do_rot
 call check_columns
 ntr = 0 ; nro = 0
 allocate(arow(16,nindepfc), rrow(16,nindepfc), trhs(16), rrhs(16)); trhs = 0; rrhs = 0
 rdrop = 0
 if (allocated(rkey)) deallocate(rkey)
 allocate(rkey(0:15,16)); rkey = 0

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
          if (c .gt. 0) then
             call add(key, c, map(rnk)%gr(g)%mat(t,ti))
          elseif (include_fc(rnk) .eq. 2 .and. full_col(rnk,g,ti) .gt. 0) then     ! read FC -> right-hand side
             call add(key, -1, map(rnk)%gr(g)%mat(t,ti)*xfull(full_col(rnk,g,ti)))
          endif
       enddo
    enddo
    enddo
    call gather_rows(arow, trhs, ntr, 'translational', rnk)
    endif
!------------------------------------ rotational, rank rnk: the rank-rnk part ...
! No rank-1 rows: the torque of FC1 over the primitive cell, sum_i Phi_i x X_i, is not an
! invariance of a crystal (moving atom i by a lattice vector L changes it by Phi_i x L). On exact
! model FCs with FC1 /= 0 these were the only rotational rows violated; dropped 2026-10-26 (user).
! FC1 still enters the rank-2 rows below.
    if (do_rot .and. rnk .ge. 2) then
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
             if (c .gt. 0) then
                call add(key, c, map(rnk)%gr(g)%mat(t,ti)*cf)
             elseif (include_fc(rnk) .eq. 2 .and. full_col(rnk,g,ti) .gt. 0) then
                call add(key, -1, map(rnk)%gr(g)%mat(t,ti)*cf*xfull(full_col(rnk,g,ti)))
             endif
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
                if (c .gt. 0) then
                   call add(key, c, map(rnk-1)%gr(g)%mat(t,ti)*lc(map(rnk-1)%gr(g)%ixyz(k,t),ak,nu))
                elseif (include_fc(rnk-1) .eq. 2 .and. full_col(rnk-1,g,ti) .gt. 0) then
                   call add(key, -1, map(rnk-1)%gr(g)%mat(t,ti)*lc(map(rnk-1)%gr(g)%ixyz(k,t),ak,nu) &
 &                          *xfull(full_col(rnk-1,g,ti)))
                endif
             enddo
          enddo
          enddo
          enddo
       enddo
       enddo
    endif
    endif
    call gather_rows(rrow, rrhs, nro, 'rotational', rnk)
    endif
    deallocate(key)
 enddo

! right-hand sides: -C_read x_read of the read ranks (include_fc=2; zero otherwise), plus the Coulomb part below
 if (any(include_fc .eq. 2)) write(ulog,'(a,es12.4)') &
 &   ' BUILD_INV_ROWS: read FCs: max |C x_read| of rows with no fitted FC (their own violation; born_flag 2: '// &
 &   'Coulomb part not included) =', rdrop
 if (do_tr) then
    transl_constraints = ntr
    if (allocated(atransl)) deallocate(atransl)
    if (allocated(btransl)) deallocate(btransl)
    allocate(atransl(ntr,nindepfc), btransl(ntr))
    atransl = arow(1:ntr,:) ; btransl = trhs(1:ntr)
    write(ulog,*)' SET_TRANSLATIONAL_INV_CONSTRAINTS: rows (all ranks) =', ntr
 endif
 if (do_rot) then
    rot_constraints = nro
    if (allocated(arot)) deallocate(arot)
    if (allocated(brot)) deallocate(brot)
    allocate(arot(nro,nindepfc), brot(nro))
    arot = rrow(1:nro,:) ; brot = rrhs(1:nro)
    write(ulog,*)' SET_ROTATIONAL_INV_CONSTRAINTS: rows (all ranks) =', nro
    if (born_flag.eq.2 .and. include_fc(2).ne.0) then           ! fitted or read (include_fc=2)
       allocate(mo(natom_prim_cell,3,3,3))
       call lr_moments(mo)
       do ir = 1, nro
          if (rkey(0,ir) .ne. 2) cycle
          i1 = rkey(1,ir); a1 = rkey(2,ir); nu = rkey(3,ir)
          cf = 0
          do a = 1, 3
          do mu = 1, 3
             cf = cf + mo(i1,a1,a,mu)*lc(a,mu,nu)
          enddo
          enddo
          brot(ir) = brot(ir) - cf
       enddo
       b2max = maxval(abs(brot))
       write(ulog,'(a,es12.4,a)')'  born_flag 2: rank-2 rotational right-hand sides from the Born-dipole FCs, max|b| =', &
 &        b2max,' eV/A'
       deallocate(mo)
    endif
    if (born_flag.eq.2 .and. include_fc(2).ne.0 .and. include_fc(3).ne.0) then
       allocate(mo(natom_prim_cell,3,3,3))
       call lr_rot3_onsite(mo, ulog)
       b3max = 0; count_log = 0; count_log2 = 0
       do ir = 1, nro
          if (rkey(0,ir) .ne. 3) cycle
          i1 = rkey(1,ir); i2 = rkey(2,ir); a1 = rkey(3,ir); a2 = rkey(4,ir); nu = rkey(5,ir)
          if (i2 .eq. i1) then
             cf = mo(i1,a1,a2,nu)
          else
             call lr_rot3_pair(atompos(:,i1)-atompos(:,i2), atom0(iatomcell0(i1))%charge, &
 &                             atom0(iatomcell0(i2))%charge, b3)
             cf = b3(a1,a2,nu)
          endif
          brot(ir) = brot(ir) + cf
          b3max = max(b3max, abs(cf))
          if (abs(cf).gt.1d-8 .and. ((i2.eq.i1 .and. count_log.lt.30) .or. (i2.ne.i1 .and. count_log2.lt.30))) then
             if (i2.eq.i1) count_log = count_log + 1          ! a sample of rows, for checking
             if (i2.ne.i1) count_log2 = count_log2 + 1
             write(ulog,'(a,5i5,3f12.6,es16.8)')'  ROT3_RHS i1,i2,a1,a2,nu, X_i2-X_i1, b =',i1,i2,a1,a2,nu, &
 &              atompos(:,i2)-atompos(:,i1), cf
          endif
       enddo
       write(ulog,'(a,es12.4,a)')'  born_flag 2: rank-3 rotational right-hand sides (Phi_B and the pair-model FC3), max|b| =', &
 &        b3max,' eV/A^2'
       deallocate(mo)
    endif
 endif
 deallocate(arow, rrow, trhs, rrhs)
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

   subroutine gather_rows(rows, rhs, nrow, label, r)
!! sorts the contributions by key, adds each key's contributions into one row and keeps it
!! if it is not zero and not a multiple of a row kept already (any rank). Contributions with column -1 are
!! read FCs (include_fc=2) times their value: their sum s gives the row's right-hand side -s.
   real(r15), allocatable, intent(inout) :: rows(:,:), rhs(:)
   integer, intent(inout) :: nrow
   character(*), intent(in) :: label
   integer, intent(in) :: r
   integer, allocatable :: ord(:), ktmp(:,:)
   real(r15), allocatable :: a1d(:), zero(:), tmp(:,:), tmp1(:)
   real(r15) sr
   integer i, j, n0, i0
   logical new
   allocate(ord(ncon), a1d(nindepfc), zero(nindepfc)) ; zero = 0
   call sort_keys(ord)
   n0 = nrow
   i = 1
   do while (i .le. ncon)
      a1d = 0; sr = 0
      j = i
      i0 = i
      do while (j .le. ncon)
         if (any(ckey(:,ord(j)) .ne. ckey(:,ord(i)))) exit
         if (ccol(ord(j)) .lt. 0) then
            sr = sr + cval(ord(j))
         else
            a1d(ccol(ord(j))) = a1d(ccol(ord(j))) + cval(ord(j))
         endif
         j = j + 1
      enddo
      i = j
      if (a1d .myeq. zero) then
         rdrop = max(rdrop, abs(sr))
         cycle
      endif
      call compare2previous_lines(size(rows,1), nindepfc, a1d, rows, nrow, new)
      if (.not. new) cycle
      if (nrow .eq. size(rows,1)) then              ! grow the row store
         allocate(tmp(2*size(rows,1), nindepfc)) ; tmp = 0
         tmp(1:nrow,:) = rows(1:nrow,:)
         call move_alloc(tmp, rows)
         allocate(tmp1(2*size(rhs))) ; tmp1 = 0
         tmp1(1:nrow) = rhs(1:nrow)
         call move_alloc(tmp1, rhs)
      endif
      nrow = nrow + 1
      rows(nrow,:) = a1d
      rhs(nrow) = -sr
      if (label .eq. 'rotational') then                ! key of the row, for its right-hand side
         if (nrow .gt. size(rkey,2)) then
            allocate(ktmp(0:15,2*size(rkey,2))) ; ktmp = 0
            ktmp(:,1:size(rkey,2)) = rkey
            call move_alloc(ktmp, rkey)
         endif
         rkey(0,nrow) = r
         rkey(1:keylen,nrow) = ckey(:,ord(i0))
      endif
   enddo
   write(ulog,'(a,a,a,i2,a,i8,a,i8)')'  ', label, ' rank', r, ': contributions', ncon, ', new rows', nrow-n0
   deallocate(ord, a1d, zero)
   end subroutine gather_rows

 end subroutine build_inv_rows
!---------------------------------------------------------------------------
 integer function column_of(rnk, g, ti)
!! column of independent term ti of group g of rank rnk in the fit, 0 if it is not a column:
!! the rank is not in the model, or an FC2 group not kept (beyond the supercell).
!! include_fc=2 (FCs read from file): not a column (0); build_inv_rows puts their terms on the right-hand side.
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
 subroutine write_map_dump(fn)
!! the map arrays at full precision, for tests against known force constants (FOCEX_DUMP_MAP=1).
!! Stream file: per included rank: int32 rank, ngr; per group: int32 nt, ntind, columns(ntind)
!! (column_of, 0 = not a column), iat(rank,nt), ixyz(rank,nt); float64 mat(nt,ntind)
!! (column-major). A term t of group g is sum_ti mat(t,ti) * x(column(ti)).
 use svd_stuff
 use params
 implicit none
 character(*), intent(in) :: fn
 integer rnk, g, ti
 open(352, file=fn, access='stream', form='unformatted', status='replace')
 do rnk = 1, maxrank
    if (include_fc(rnk) .ne. 1) cycle
    write(352) int(rnk,4), int(map(rnk)%ngr,4)
    do g = 1, map(rnk)%ngr
       write(352) int(map(rnk)%nt(g),4), int(map(rnk)%ntind(g),4), &
&           (int(column_of(rnk,g,ti),4), ti=1,map(rnk)%ntind(g))
       write(352) int(map(rnk)%gr(g)%iat(1:rnk,1:map(rnk)%nt(g)),4), &
&           int(map(rnk)%gr(g)%ixyz(1:rnk,1:map(rnk)%nt(g)),4)
       write(352) real(map(rnk)%gr(g)%mat(1:map(rnk)%nt(g),1:map(rnk)%ntind(g)),8)
    enddo
 enddo
 close(352)
 end subroutine write_map_dump
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

!---------------------------------------------------------------------------
 subroutine merge_given_fcs(uio)
!! include_fc=2: fitted FCs (columns column_of) and read FCs (xfull) into one vector in the full layout (full_col),
!! fcs and sigma (0 for the read ones) reallocated; the read ranks are then marked include_fc=1 and nindepfc,
!! dim_ac set to the full size, so that the output and everything after it work as for a full fit.
 use svd_stuff
 use params
 use given_fcs, only : xfull, nfull, full_col
 implicit none
 integer, intent(in) :: uio
 real(r15), allocatable :: f(:), s(:)
 integer r, g, ti, c, k
 allocate(f(nfull), s(nfull)); f = 0; s = 0
 do r = 1, maxrank
    if (include_fc(r) .eq. 0) cycle
    do g = 1, map(r)%ngr
       do ti = 1, map(r)%ntind(g)
          k = full_col(r, g, ti)
          if (k .eq. 0) cycle
          if (include_fc(r) .eq. 1) then
             c = column_of(r, g, ti)
             f(k) = fcs(c); s(k) = sigma(c)
          else
             f(k) = xfull(k)
          endif
       enddo
    enddo
 enddo
 call move_alloc(f, fcs); call move_alloc(s, sigma)
 do r = 1, maxrank
    if (include_fc(r) .eq. 2) then
       include_fc(r) = 1
       write(uio,*)'MERGE_GIVEN_FCS: rank ',r,' (read) merged with the fitted ranks; marked include_fc=1 from here on'
    endif
 enddo
 nindepfc = nfull; dim_ac = nfull
 end subroutine merge_given_fcs

 end module invariance_ranks
