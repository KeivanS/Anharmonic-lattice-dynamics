!******************************************************************************
! Range of the clusters of rank >= 3 (user's rule, 2026-10-24), used by
! size_force_constants and collect_force_constants. structure.params line 12 is
! "fc2range [fc3flag ... fc8flag]" (fcflag(3:8), default 1 when absent).
! fcflag(rank) /= 0: the cluster is kept if all its pairs lie within nshells(rank),
! as before. fcflag(rank) = 0: the cluster is kept if
!  - every pair of it is a kept FC2 pair (its FC2 group has a member inside the
!    WS cell of the supercell: the test of setup_FC2_in_supercell for fc2range=0,
!    of exclude_beyond_sc otherwise), and
!  - all its pairs lie within nshells(rank), or it has 2 distinct atoms
!    (Phi_iij, Phi_iiij, Phi_iijj), or its 3 or more distinct atoms lie on one line
!    (every triangle of zero area).
! Reason: interatomic potentials are mostly pair terms plus a smaller correction, so
! the two-body parts of FC3/FC4 can reach as far as FC2; LASSO removes the small ones.
! "Within nshells" uses the shell counts exactly as the original loop does: the
! shell of atom a seen from the primitive atom equivalent to b (via find_atompos).
!******************************************************************************
 module fc34_range
 use constants, only : r15
 implicit none
 logical, allocatable :: pair_keep(:,:)   ! (primitive atom, atom of the list): kept FC2 pair
 real(r15) :: fc34_r2max = 0     ! largest kept FC2 pair distance (Ang), for the pruning
 real(r15), parameter :: col_tol = 1d-5
 contains
!------------------------------------------------------------------------------
 subroutine fc34_init
!! builds pair_keep from the FC2 groups; call after map(2) is filled, before rank 3
 use atoms_force_constants
 use svd_stuff, only : map
 use params, only : include_fc
 use lattice, only : rws26
 use ios, only : ulog
 integer g, ti, t, i0, j, m, nb, rnk, npair
 integer icell(3)
 real(r15) rij(3)
 logical inside, gkeep
 if (allocated(pair_keep)) deallocate(pair_keep)
 if (all(fcflag(3:8).ne.0 .or. include_fc(3:8).eq.0)) return
 if (include_fc(2).eq.0 .or. map(2)%ngr.eq.0) then
    write(ulog,*)'FC34_INIT: fcflag=0 needs the FC2 groups (include_fc(2)/=0); using nshells for all ranks'
    write(*   ,*)'FC34_INIT: fcflag=0 needs the FC2 groups (include_fc(2)/=0); using nshells for all ranks'
    fcflag(3:8) = 1
    return
 endif
 allocate(pair_keep(natom_prim_cell,natoms)); pair_keep = .false.
 do g = 1, map(2)%ngr
    gkeep = .false.
    tiloop: do ti = 1, map(2)%ntind(g)
       do t = 1, map(2)%nt(g)
          if (fc2range.eq.0 .and. map(2)%gr(g)%mat(t,ti).eq.0) cycle
          i0 = map(2)%gr(g)%iat(1,t)
          j  = map(2)%gr(g)%iat(2,t)
          rij = atompos(:,j) - atompos(:,i0)
          call check_inside_ws(rij,rws26,inside,nb)
          if (inside) then ; gkeep = .true. ; exit tiloop ; endif
       enddo
    enddo tiloop
    if (.not.gkeep) cycle
    do t = 1, map(2)%nt(g)
       i0 = map(2)%gr(g)%iat(1,t)
       j  = map(2)%gr(g)%iat(2,t)
       icell = iatomcell(:,j) - iatomcell(:,i0)
       call find_atompos(icell, iatomcell0(j), m)
       if (m.ne.0) pair_keep(iatomcell0(i0),m) = .true.
       icell = -icell
       call find_atompos(icell, iatomcell0(i0), m)
       if (m.ne.0) pair_keep(iatomcell0(j),m) = .true.
    enddo
 enddo
 fc34_r2max = 0 ; npair = 0
 do t = 1, natom_prim_cell
    do j = 1, natoms
       if (.not.pair_keep(t,j) .or. iatomcell0(j).eq.t .and. all(iatomcell(:,j).eq.0)) cycle
       npair = npair + 1
       fc34_r2max = max(fc34_r2max, sqrt(sum((atompos(:,j)-atompos(:,t))**2)))
    enddo
 enddo
 do rnk = 3, 8
    if (include_fc(rnk).ne.0) write(ulog,'(a,i2,a,i2)')' FC34_INIT: rank',rnk,' fcflag=',fcflag(rnk)
 enddo
 write(ulog,*)'FC34_INIT: kept FC2 pairs (from the primitive atoms)=',npair,' largest distance=',fc34_r2max
 write(*   ,*)'FC34_INIT: kept FC2 pairs (from the primitive atoms)=',npair,' largest distance=',fc34_r2max
 end subroutine fc34_init
!------------------------------------------------------------------------------
 logical function fc34_active(nrank)
!! the rule of this module applies to rank nrank
 use atoms_force_constants, only : fcflag
 integer, intent(in) :: nrank
 fc34_active = .false.
 if (nrank.lt.3 .or. nrank.gt.8 .or. .not.allocated(pair_keep)) return
 fc34_active = fcflag(nrank).eq.0
 end function fc34_active
!------------------------------------------------------------------------------
 integer function shell_rel(iref, ioth)
!! shell index of atom ioth seen from the primitive atom equivalent to iref, as in
!! collect_force_constants; huge if the translated atom is not in the list
 use atoms_force_constants
 integer, intent(in) :: iref, ioth
 integer icell(3), m
 icell = iatomcell(:,ioth) - iatomcell(:,iref)
 call find_atompos(icell, iatomcell0(ioth), m)
 if (m .eq. 0) then
    shell_rel = huge(1)
 else
    shell_rel = iatomneighbor(iatomcell0(iref), m)
 endif
 end function shell_rel
!------------------------------------------------------------------------------
 logical function pair_kept(ia, ib)
!! the pair (ia, ib) of atoms of the list is a kept FC2 pair
 use atoms_force_constants
 integer, intent(in) :: ia, ib
 integer icell(3), m
 pair_kept = .false.
 icell = iatomcell(:,ib) - iatomcell(:,ia)
 call find_atompos(icell, iatomcell0(ib), m)
 if (m.ne.0) pair_kept = pair_keep(iatomcell0(ia), m)
 if (pair_kept) return
 icell = -icell
 call find_atompos(icell, iatomcell0(ia), m)
 if (m.ne.0) pair_kept = pair_keep(iatomcell0(ib), m)
 end function pair_kept
!------------------------------------------------------------------------------
 logical function old_ok(n, at, nsh)
!! the original test on a cluster at(1:n), at(1) the atom in the cell at origin:
!! pair (1,q) seen from at(1), pair (p,q), p >= 2, seen from at(q); all <= nsh
 integer, intent(in) :: n, at(n), nsh
 integer p, q
 old_ok = .false.
 do q = 2, n
    if (shell_rel(at(1), at(q)) .gt. nsh) return
    do p = 2, q-1
       if (shell_rel(at(q), at(p)) .gt. nsh) return
    enddo
 enddo
 old_ok = .true.
 end function old_ok
!------------------------------------------------------------------------------
 subroutine distinct(n, at, nd, d)
 integer, intent(in) :: n, at(n)
 integer, intent(out) :: nd, d(n)
 integer i
 nd = 0
 do i = 1, n
    if (any(d(1:nd) .eq. at(i))) cycle
    nd = nd + 1 ; d(nd) = at(i)
 enddo
 end subroutine distinct
!------------------------------------------------------------------------------
 logical function collinear(nd, d)
!! every atom of d(1:nd) lies on the line through its farthest pair
 use atoms_force_constants, only : atompos
 integer, intent(in) :: nd, d(nd)
 integer p, q, ia, ib
 real(r15) dmax, r, u(3), v(3), c(3)
 dmax = -1 ; ia = d(1) ; ib = d(1)
 do p = 1, nd
    do q = p+1, nd
       r = sqrt(sum((atompos(:,d(p))-atompos(:,d(q)))**2))
       if (r .gt. dmax) then ; dmax = r ; ia = d(p) ; ib = d(q) ; endif
    enddo
 enddo
 collinear = .true.
 u = atompos(:,ib) - atompos(:,ia)
 do p = 1, nd
    v = atompos(:,d(p)) - atompos(:,ia)
    c = [u(2)*v(3)-u(3)*v(2), u(3)*v(1)-u(1)*v(3), u(1)*v(2)-u(2)*v(1)]
    if (sqrt(sum(c**2)) .gt. col_tol*sum(u**2)) then ; collinear = .false. ; return ; endif
 enddo
 end function collinear
!------------------------------------------------------------------------------
 logical function all_pairs_kept(nd, d)
 integer, intent(in) :: nd, d(nd)
 integer p, q
 all_pairs_kept = .false.
 do p = 1, nd
    do q = p+1, nd
       if (.not.pair_kept(d(p), d(q))) return
    enddo
 enddo
 all_pairs_kept = .true.
 end function all_pairs_kept
!------------------------------------------------------------------------------
 logical function fc34_partial(n, at, nsh)
!! necessary condition on a partial cluster at(1:n) (the loop fills the remaining
!! positions with at(n)): some completion can pass fc34_final
 integer, intent(in) :: n, at(n), nsh
 integer nd, d(n)
 call distinct(n, at, nd, d)
 fc34_partial = all_pairs_kept(nd, d)
 if (.not.fc34_partial) return
 fc34_partial = nd.le.2
 if (fc34_partial) return
 fc34_partial = old_ok(n, at, nsh)
 if (fc34_partial) return
 fc34_partial = collinear(nd, d)
 end function fc34_partial
!------------------------------------------------------------------------------
 logical function fc34_final(n, at, nsh)
!! the rule of this module on a complete cluster at(1:n)
 integer, intent(in) :: n, at(n), nsh
 fc34_final = fc34_partial(n, at, nsh)
 end function fc34_final
!------------------------------------------------------------------------------
 subroutine write_clusters
!! clusters<rank>.dat for rank >= 2: one line per edge of each distinct atom set
!! (first atom in the primitive cell, positions relative to it), for display:
!! x y z dx dy dz  nd  type  group  tau   (type 0 within nshells, 1 two-body, 2 on a line,
!! 3 the other ones -- not possible under the rule; tau the primitive atom)
 use atoms_force_constants
 use svd_stuff, only : map
 use params, only : include_fc, nshells
 integer rnk, g, t, nd, d(8), p, q, ityp, u, nset, k
 integer, allocatable :: seen(:,:)
 real(r15) r0(3)
 character xt*1
 do rnk = 2, 8
    if (include_fc(rnk).eq.0) cycle
    if (map(rnk)%ngr.eq.0) cycle
    write(xt,'(i1)') rnk
    u = 9301
    open(u, file='clusters'//xt//'.dat')
    write(u,'(a)')'# x y z dx dy dz (Ang, from the primitive atom) nd type group tau;'// &
 &     ' type 0 within nshells, 1 two-body, 2 on a line; tau the primitive atom'
    do g = 1, map(rnk)%ngr     ! an atom set belongs to one group: no duplicates across groups
       allocate(seen(8, map(rnk)%nt(g))) ; nset = 0
       tloop: do t = 1, map(rnk)%nt(g)
          call distinct(rnk, map(rnk)%gr(g)%iat(:,t), nd, d)
          if (nd.lt.2) cycle
          call sort_int(nd, d)
          d(nd+1:) = 0
          do k = 1, nset
             if (all(seen(:,k).eq.d)) cycle tloop
          enddo
          nset = nset + 1 ; seen(:,nset) = d
          if (rnk.eq.2 .or. old_ok(rnk, map(rnk)%gr(g)%iat(:,t), &
 &            nshells(rnk, iatomcell0(map(rnk)%gr(g)%iat(1,t))))) then
             ityp = 0
          elseif (nd.eq.2) then
             ityp = 1
          elseif (collinear(nd, d)) then
             ityp = 2
          else
             ityp = 3
          endif
          r0 = atompos(:, map(rnk)%gr(g)%iat(1,t))
          do p = 1, nd
             do q = p+1, nd
                write(u,'(6(1x,f10.4),4(1x,i5))') atompos(:,d(p))-r0, &
 &                  atompos(:,d(q))-atompos(:,d(p)), nd, ityp, g, iatomcell0(map(rnk)%gr(g)%iat(1,t))
             enddo
          enddo
       enddo tloop
       deallocate(seen)
    enddo
    close(u)
 enddo
 contains
 subroutine sort_int(n, a)
 integer, intent(in) :: n
 integer, intent(inout) :: a(n)
 integer i, j, x
 do i = 2, n
    x = a(i) ; j = i - 1
    do while (j.ge.1)
       if (a(j).le.x) exit
       a(j+1) = a(j) ; j = j - 1
    enddo
    a(j+1) = x
 enddo
 end subroutine sort_int
 end subroutine write_clusters
 end module fc34_range
