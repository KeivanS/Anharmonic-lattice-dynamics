!******************************************************************************
! Clusters of the ranks that are read (include_fc(rank)=2; user, 2026-11-11).
! The groups of a read rank are built from the atom clusters listed in its fcN.dat, whatever line 12
! (fcNflag) and the shell numbers of structure.params say: setup_maps neither adds a cluster to the read
! list nor drops one from it. A cluster is the multiset of the rank's atoms of a term of fcN.dat, as
! (sublattice tau, cell n) through the lat_fc.dat of the run that wrote it, up to a lattice translation.
! During the enumeration of setup_maps_fc34 (size_ and collect_force_constants_fc34) a candidate atom j
! is followed only if it belongs to some read cluster that contains the primitive atom iatom0 at the
! origin (rdc_partial), and a complete candidate is kept only if it is a read cluster (rdc_final).
! rdc_report then checks, from the finished maps, that every read cluster became a group (none dropped)
! and counts the clusters this run's symmetry put into the same groups although the file does not list
! them (the file was written with another symmetry, or their terms were below the write-out margin).
! fcN.dat written since 2026-11-11 lists also the groups none of whose terms reaches the margin (by their
! largest term), so every group of the writing run is read back; older files lack those (FCs < margin).
! A read cluster that reaches beyond this run's neighbour list stops the run (raise icutoff, line 3 of
! default.params): it cannot be represented, and dropping it is not allowed.
!******************************************************************************
 module fc_read_clusters
 use constants, only : r15
 implicit none
 private
 public :: rdc_init, rdc_active, rdc_partial, rdc_final, rdc_report

 type rdset
    integer :: n = 0
    integer(8), allocatable :: key(:,:)    ! (4, n) canonical cluster codes, sorted lexicographically
    logical, allocatable :: found(:)       ! the enumeration reached it
    logical, allocatable :: allow(:,:)     ! (primitive atom, atom of the list): in a read cluster with it
 end type rdset
 type(rdset), save :: rd(4)
 logical, save :: on(8) = .false.
 integer, parameter :: ncmax = 63          ! largest |cell offset| inside a cluster (7 bits per component)

 contains
!------------------------------------------------------------------------------
 logical function rdc_active(nrank)
!! rank nrank is read and its clusters come from fcN.dat
 integer, intent(in) :: nrank
 rdc_active = .false.
 if (nrank.ge.1 .and. nrank.le.8) rdc_active = on(nrank)
 end function rdc_active
!------------------------------------------------------------------------------
 logical function rdc_partial(nrank, iatom0, j)
!! atom j of the list can join a cluster whose first atom is the primitive atom iatom0 at the origin
 integer, intent(in) :: nrank, iatom0, j
 rdc_partial = rd(nrank)%allow(iatom0, j)
 end function rdc_partial
!------------------------------------------------------------------------------
 logical function rdc_final(nrank, at)
!! the atoms at(1:nrank) of the list form a read cluster
 use atoms_force_constants, only : iatomcell0, iatomcell
 integer, intent(in) :: nrank, at(nrank)
 integer tau(nrank), nc(3,nrank), k, i
 integer(8) key(4)
 do k = 1, nrank
    tau(k) = iatomcell0(at(k)); nc(:,k) = iatomcell(:,at(k))
 enddo
 call canonical(nrank, tau, nc, key, i)
 rdc_final = .false.
 if (i .ne. 0) return                      ! wider than any read cluster
 i = locate(rd(nrank), key)
 if (i .eq. 0) return
 rd(nrank)%found(i) = .true.
 rdc_final = .true.
 end function rdc_final
!------------------------------------------------------------------------------
 subroutine rdc_init
!! reads the clusters of every rank with include_fc=2 from fcN.dat and the lat_fc.dat of the run that wrote them
 use params, only : include_fc
 use ios, only : ulog
 use atoms_force_constants, only : natoms, natom_prim_cell, iatomcell0, iatomcell
 use lattice, only : r01, r02, r03
 use geometry, only : v2a
 integer, parameter :: ulat=193, ufc=194
 integer r, i, j, k, a, nold, nat0, iold, ios_, t, g, ia(8), ix(8), nkey, mxkey, nlo(3), nhi(3), nn(3), idx, ierr
 integer, allocatable :: otau(:), on3(:,:), lookup(:,:,:,:), ord(:)
 integer(8), allocatable :: keys(:,:), tmp(:,:)
 integer tau(8), nc(3,8)
 integer(8) key(4)
 real(r15) rold(3,3), posold(3), val
 character(len=300) line
 character(len=1) rn
 logical ex

 on = .false.
 if (all(include_fc .ne. 2)) return
 do r = 5, 8
    if (include_fc(r) .eq. 2) then
       write(*,*)'RDC_INIT: include_fc=2 is implemented for ranks 1-4 only; rank ',r
       write(ulog,*)'RDC_INIT: include_fc=2 is implemented for ranks 1-4 only; rank ',r
       stop
    endif
 enddo

! ---- lat_fc.dat of the run that wrote the FCs: atom number -> (tau, n)
 inquire(file='lat_fc.dat', exist=ex)
 if (.not. ex) then
    write(*,*)'RDC_INIT: include_fc=2 needs lat_fc.dat (and fcN.dat) of the run that produced the FCs'
    write(ulog,*)'RDC_INIT: include_fc=2 needs lat_fc.dat (and fcN.dat) of the run that produced the FCs'
    stop
 endif
 open(ulat, file='lat_fc.dat', status='old')
 read(ulat,'(a)') line
 do i = 1, 3
    read(ulat,*) rold(:,i)
 enddo
 if (maxval(abs(rold(:,1)-v2a(r01))) + maxval(abs(rold(:,2)-v2a(r02))) + maxval(abs(rold(:,3)-v2a(r03))) .gt. 1d-4) then
    write(*,*)'RDC_INIT: the primitive vectors of lat_fc.dat differ from structure.params'
    write(ulog,*)'RDC_INIT: the primitive vectors of lat_fc.dat differ from structure.params'
    stop
 endif
 read(ulat,'(a)') line
 read(ulat,*) nat0
 if (nat0 .ne. natom_prim_cell) then
    write(*,*)'RDC_INIT: lat_fc.dat has ',nat0,' atoms in the primitive cell, structure.params ',natom_prim_cell
    stop
 endif
 do i = 1, nat0
    read(ulat,'(a)') line
 enddo
 do
    read(ulat,'(a)',iostat=ios_) line
    if (ios_ .ne. 0) then
       write(*,*)'RDC_INIT: no atom list (maxshells,Neighborshell atoms) in lat_fc.dat'
       stop
    endif
    if (index(line,'Neighborshell atoms') .gt. 0) exit
 enddo
 read(ulat,*) i, nold, j
 allocate(otau(nold), on3(3,nold))
 do i = 1, nold
    read(ulat,*) iold, posold, otau(i), on3(:,i)
 enddo
 close(ulat)

! ---- this run's atoms by (tau, n)
 nlo = minval(iatomcell(:,1:natoms),2); nhi = maxval(iatomcell(:,1:natoms),2)
 allocate(lookup(natom_prim_cell, nlo(1):nhi(1), nlo(2):nhi(2), nlo(3):nhi(3))); lookup = 0
 do i = 1, natoms
    lookup(iatomcell0(i), iatomcell(1,i), iatomcell(2,i), iatomcell(3,i)) = i
 enddo

 do r = 1, 4
    if (include_fc(r) .ne. 2) cycle
    write(rn,'(i1)') r
    inquire(file='fc'//rn//'.dat', exist=ex)
    if (.not. ex) then
       write(*,*)'RDC_INIT: include_fc(',r,')=2 but fc'//rn//'.dat is missing'
       write(ulog,*)'RDC_INIT: include_fc(',r,')=2 but fc'//rn//'.dat is missing'
       stop
    endif
    mxkey = 4096; nkey = 0
    allocate(keys(4,mxkey))
    open(ufc, file='fc'//rn//'.dat', status='old')
    do
       read(ufc,'(a)',iostat=ios_) line
       if (ios_ .ne. 0) exit
       if (len_trim(line) .eq. 0) cycle
       if (index(adjustl(line),'#') .eq. 1) then
          if (index(line,'margin=') .gt. 0) write(ulog,'(a)') ' RDC_INIT: fc'//rn//'.dat header:'//trim(line)
          cycle
       endif
       read(line,*,iostat=ios_) t, g, (ia(j), ix(j), j=1,r), val
       if (ios_ .ne. 0) cycle
       do j = 1, r
          if (ia(j) .lt. 1 .or. ia(j) .gt. nold) then
             write(*,*)'RDC_INIT: fc'//rn//'.dat names atom ',ia(j),' not in lat_fc.dat (1..',nold,')'
             stop
          endif
          tau(j) = otau(ia(j)); nc(:,j) = on3(:,ia(j))
       enddo
       call canonical(r, tau(1:r), nc(:,1:r), key, ierr)
       if (ierr .ne. 0) then
          write(*,*)'RDC_INIT: a cluster of fc'//rn//'.dat spans more than ',ncmax,' cells along a primitive vector'
          write(ulog,*)'RDC_INIT: a cluster of fc'//rn//'.dat spans more than ',ncmax,' cells along a primitive vector'
          stop
       endif
       if (nkey .eq. mxkey) then
          allocate(tmp(4,2*mxkey)); tmp(:,1:nkey) = keys(:,1:nkey)
          call move_alloc(tmp, keys); mxkey = 2*mxkey
       endif
       nkey = nkey + 1; keys(:,nkey) = key
    enddo
    close(ufc)
! sort, remove duplicates
    allocate(ord(nkey))
    call sort_keys(nkey, keys, ord)
    allocate(tmp(4,nkey)); k = 0
    do i = 1, nkey
       if (k .gt. 0) then
          if (all(keys(:,ord(i)) .eq. tmp(:,k))) cycle
       endif
       k = k + 1; tmp(:,k) = keys(:,ord(i))
    enddo
    rd(r)%n = k
    if (allocated(rd(r)%key)) deallocate(rd(r)%key, rd(r)%found, rd(r)%allow)
    allocate(rd(r)%key(4,k), rd(r)%found(k), rd(r)%allow(natom_prim_cell,natoms))
    rd(r)%key = tmp(:,1:k); rd(r)%found = .false.; rd(r)%allow = .false.
    deallocate(tmp, ord, keys)
! allow(tau_a, j): j is a member of a read cluster translated so that its member a sits in the cell at the origin
    do i = 1, rd(r)%n
       call decode(r, rd(r)%key(:,i), tau, nc)
       do a = 1, r
          do j = 1, r
             nn = nc(:,j) - nc(:,a)
             idx = 0
             if (all(nn .ge. nlo) .and. all(nn .le. nhi)) idx = lookup(tau(j), nn(1), nn(2), nn(3))
             if (idx .eq. 0) then
                write(*,*)'RDC_INIT: a cluster of fc'//rn//'.dat reaches beyond this run''s neighbour list:'
                write(*,*)'          raise icutoff (default.params line 3); a read FC is never dropped'
                write(ulog,*)'RDC_INIT: a cluster of fc'//rn//'.dat reaches beyond this run''s neighbour list:'
                write(ulog,*)'          raise icutoff (default.params line 3); a read FC is never dropped'
                write(ulog,*)'          cluster (tau, cell):',(tau(k), nc(:,k), k=1,r)
                stop
             endif
             rd(r)%allow(tau(a), idx) = .true.
          enddo
       enddo
    enddo
    on(r) = .true.
    write(ulog,'(a,i2,a,i8,a)')' RDC_INIT: rank',r,': ',rd(r)%n, &
 &       ' clusters read from fc'//rn//'.dat; the groups of this rank are built from them only'
    write(*,'(a,i2,a,i8,a)')' RDC_INIT: rank',r,': ',rd(r)%n,' clusters read from fc'//rn//'.dat'
 enddo
 deallocate(otau, on3, lookup)
 end subroutine rdc_init
!------------------------------------------------------------------------------
 subroutine rdc_report
!! after setup_maps: every read cluster must be a group (none dropped); clusters of the groups that the
!! file does not list are counted (added by this run's symmetry)
 use svd_stuff, only : map
 use ios, only : ulog
 use atoms_force_constants, only : iatomcell0, iatomcell
 integer r, g, t, k, nmissing, nadded, ierr, nseen, i
 integer tau(8), nc(3,8)
 integer(8) key(4)
 integer(8), allocatable :: added(:,:)
 do r = 1, 4
    if (.not. on(r)) cycle
    nmissing = count(.not. rd(r)%found)
    if (nmissing .gt. 0) then
       write(*,*)'RDC_REPORT: ERROR: ',nmissing,' clusters of the read fc file of rank ',r,' were not built into groups'
       write(ulog,*)'RDC_REPORT: ERROR: ',nmissing,' clusters of the read fc file of rank ',r,' were not built into groups'
       do i = 1, rd(r)%n
          if (rd(r)%found(i)) cycle
          call decode(r, rd(r)%key(:,i), tau, nc)
          write(ulog,*)'   missing cluster (tau, cell):',(tau(k), nc(:,k), k=1,r)
       enddo
       stop
    endif
    nadded = 0; nseen = 0
    allocate(added(4, max(1,sum(map(r)%nt(1:map(r)%ngr)))))
    do g = 1, map(r)%ngr
       do t = 1, map(r)%nt(g)
          do k = 1, r
             tau(k) = iatomcell0(map(r)%gr(g)%iat(k,t)); nc(:,k) = iatomcell(:,map(r)%gr(g)%iat(k,t))
          enddo
          call canonical(r, tau(1:r), nc(:,1:r), key, ierr)
          if (ierr .eq. 0) then
             if (locate(rd(r), key) .ne. 0) cycle
          endif
          do i = 1, nseen
             if (all(added(:,i) .eq. key)) exit
          enddo
          if (i .le. nseen) cycle
          nseen = nseen + 1; added(:,nseen) = key
       enddo
    enddo
    nadded = nseen
    deallocate(added)
    write(ulog,'(a,i2,a,i8,a,i8,a)')' RDC_REPORT: rank',r,': all ',rd(r)%n,' read clusters are groups; ',nadded, &
 &       ' clusters not in the file share their groups (symmetry images; their FCs follow from the read ones)'
    if (nadded .gt. 0) write(*,'(a,i2,a,i8,a)')' RDC_REPORT: rank',r,': ',nadded, &
 &       ' clusters not listed in the read file are symmetry images of read ones (see the log)'
 enddo
 end subroutine rdc_report
!------------------------------------------------------------------------------
 subroutine canonical(r, tau, nc, key, ierr)
!! translation-invariant code of the multiset of atoms (tau(k), nc(:,k)), k=1..r: for each member as the origin,
!! the sorted atom codes; the lexicographically smallest of them. ierr=1 if an offset exceeds ncmax.
 integer, intent(in) :: r, tau(r), nc(3,r)
 integer(8), intent(out) :: key(4)
 integer, intent(out) :: ierr
 integer a, k, l
 integer(8) c(4), s
 ierr = 0
 key = huge(1_8)
 do a = 1, r
    c = -1_8
    do k = 1, r
       if (any(abs(nc(:,k)-nc(:,a)) .gt. ncmax)) then
          ierr = 1; return
       endif
       c(k) = code(tau(k), nc(:,k)-nc(:,a))
    enddo
    do k = 2, r                                   ! insertion sort, r <= 4
       s = c(k); l = k - 1
       do while (l .ge. 1)
          if (c(l) .le. s) exit
          c(l+1) = c(l); l = l - 1
       enddo
       c(l+1) = s
    enddo
    if (lexless(c, key)) key = c
 enddo
 end subroutine canonical
!------------------------------------------------------------------------------
 integer(8) function code(tau, n)
 integer, intent(in) :: tau, n(3)
 code = ((int(tau,8)*128_8 + int(n(1)+64,8))*128_8 + int(n(2)+64,8))*128_8 + int(n(3)+64,8)
 end function code
!------------------------------------------------------------------------------
 subroutine decode(r, key, tau, nc)
 integer, intent(in) :: r
 integer(8), intent(in) :: key(4)
 integer, intent(out) :: tau(8), nc(3,8)
 integer k
 integer(8) c
 do k = 1, r
    c = key(k)
    nc(3,k) = int(mod(c,128_8)) - 64; c = c/128_8
    nc(2,k) = int(mod(c,128_8)) - 64; c = c/128_8
    nc(1,k) = int(mod(c,128_8)) - 64; c = c/128_8
    tau(k) = int(c)
 enddo
 end subroutine decode
!------------------------------------------------------------------------------
 logical function lexless(a, b)
 integer(8), intent(in) :: a(4), b(4)
 integer k
 lexless = .false.
 do k = 1, 4
    if (a(k) .lt. b(k)) then
       lexless = .true.; return
    elseif (a(k) .gt. b(k)) then
       return
    endif
 enddo
 end function lexless
!------------------------------------------------------------------------------
 integer function locate(s, key)
!! index of key in the sorted list s%key, 0 if absent (binary search)
 type(rdset), intent(in) :: s
 integer(8), intent(in) :: key(4)
 integer lo, hi, mid
 locate = 0
 lo = 1; hi = s%n
 do while (lo .le. hi)
    mid = (lo + hi)/2
    if (all(s%key(:,mid) .eq. key)) then
       locate = mid; return
    elseif (lexless(s%key(:,mid), key)) then
       lo = mid + 1
    else
       hi = mid - 1
    endif
 enddo
 end function locate
!------------------------------------------------------------------------------
 subroutine sort_keys(n, keys, ord)
!! ord: the order of the columns of keys(4,n), lexicographically ascending (bottom-up merge sort, stable)
 integer, intent(in) :: n
 integer(8), intent(in) :: keys(4,n)
 integer, intent(out) :: ord(n)
 integer, allocatable :: w(:)
 integer width, lo, mid, hi, i, j, k
 do i = 1, n
    ord(i) = i
 enddo
 allocate(w(n))
 width = 1
 do while (width .lt. n)
    lo = 1
    do while (lo .le. n)
       mid = min(lo + width - 1, n); hi = min(lo + 2*width - 1, n)
       i = lo; j = mid + 1; k = lo
       do while (i .le. mid .and. j .le. hi)
          if (lexless(keys(:,ord(j)), keys(:,ord(i)))) then
             w(k) = ord(j); j = j + 1
          else
             w(k) = ord(i); i = i + 1
          endif
          k = k + 1
       enddo
       do while (i .le. mid)
          w(k) = ord(i); i = i + 1; k = k + 1
       enddo
       do while (j .le. hi)
          w(k) = ord(j); j = j + 1; k = k + 1
       enddo
       lo = lo + 2*width
    enddo
    ord = w
    width = 2*width
 enddo
 deallocate(w)
 end subroutine sort_keys
 end module fc_read_clusters
