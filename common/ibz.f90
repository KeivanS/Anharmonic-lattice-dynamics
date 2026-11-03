!============================================================
 module ibz_rules
!! Irreducible wedge of a k-mesh, shared by FOCEX and THERMACOND (2026-11-03, user: one routine, THERMACOND's).
!! The algorithm is THERMACOND's get_weights of 2026-10-22:
!!  - operations: the Laue group of the crystal (laue_operations): the rotations of its space-group operations
!!    (screw axes and glide planes included) and their products with inversion (time reversal, omega(-k) = omega(k)).
!!    Operations of the lattice point group outside it are left out: for a crystal of lower symmetry than its lattice
!!    they would make inequivalent k-points equivalent (the old FOCEX get_weights used them all).
!!  - each mesh point not yet assigned starts a star; every star member is mapped onto its mesh point modulo G, so a
!!    star is one IBZ point with weight = its number of mesh points (weights normalized to 1).
!!  - representative: the member, folded into the Wigner-Seitz cell, with the largest u.k for a fixed generic
!!    direction u (0.9, sqrt2-1, sqrt5-2): {k: u.k >= u.Rk for all R} is a fundamental domain of the point group, so
!!    the IBZ points fill one connected wedge around u. Ties (zone boundary) by the largest u.k, then (x, y, z).
!! Mesh lookup: reduced coordinates k.r_i/2pi modulo 1 on a 1e-6 grid, sorted; works for any regular mesh.
 implicit none
 private
 public :: laue_operations, ibz_reduce

 real(8), parameter :: uwedge(3) = (/ 0.9d0, 0.41421356237309505d0, 0.23606797749978970d0 /)

 contains
!------------------------------------------------------------
 subroutine laue_operations(nlat, op, nsg, isg, nl, rot)
!! rot(:,:,1:nl): the operations op(:,:,i), i = 1..nlat (lattice point group, Cartesian), that are the rotation
!! part of a space-group operation (isg(1:nsg)) or minus one
 integer, intent(in) :: nlat, nsg, isg(:)
 real(8), intent(in) :: op(3,3,*)
 integer, intent(out) :: nl
 real(8), intent(out) :: rot(3,3,48)
 integer i, j
 logical laue
 nl = 0
 do i = 1, nlat
    laue = .false.
    do j = 1, nsg
       if (isg(j) .eq. i) laue = .true.
       if (all(abs(op(:,:,i)+op(:,:,isg(j))) .lt. 1d-6)) laue = .true.
    enddo
    if (.not. laue) cycle
    nl = nl + 1
    rot(:,:,nl) = op(:,:,i)
 enddo
 end subroutine laue_operations
!------------------------------------------------------------
 subroutine ibz_reduce(nk, kp, r1, r2, r3, nrot, rot, ulog, nibz, kibz, wibz, mapibz, mapinv)
!! kp(3,nk): the mesh (Cartesian, 1/A, with 2pi); r1, r2, r3: primitive vectors (Cartesian, A); rot: Laue operations.
!! Outputs: nibz, kibz(3,nibz) (representatives, Cartesian, in the WS cell), wibz(nibz) (normalized), mapibz(nk)
!! (IBZ point of each mesh point), mapinv(nibz) (mesh index of each representative, modulo G)
 integer, intent(in) :: nk, nrot, ulog
 real(8), intent(in) :: kp(3,nk), r1(3), r2(3), r3(3), rot(3,3,nrot)
 integer, intent(out) :: nibz, mapibz(nk), mapinv(nk)
 real(8), allocatable, intent(out) :: kibz(:,:), wibz(:)
 real(8) pi, rr(3,3), gg(3,3), q(3), qws(3), krep(3), best, f, sw
 real(8), allocatable :: k2(:,:), w2(:)
 integer(8), allocatable :: keys(:,:)
 integer, allocatable :: ord(:)
 integer i, j, l, jrep
 pi = 4*atan(1d0)
 rr(:,1) = r1; rr(:,2) = r2; rr(:,3) = r3
 call reciprocal(rr, gg)                         ! gg(:,i).rr(:,j) = 2 pi delta_ij
 allocate(keys(3,nk), ord(nk), k2(3,nk), w2(nk))
 do i = 1, nk
    keys(:,i) = mesh_key(kp(:,i))
 enddo
 call sort_keys(nk, keys, ord)
 nibz = 0; mapibz = 0; mapinv = 0; w2 = 0
 kploop: do i = 1, nk
    if (mapibz(i) .ne. 0) cycle kploop
    nibz = nibz + 1
    best = -huge(1d0); jrep = i; krep = kp(:,i)
    lloop: do l = 1, nrot
       q = matmul(rot(:,:,l), kp(:,i))
       j = find_key(nk, keys, ord, mesh_key(q))
       if (j .eq. 0) cycle lloop                     ! star member not on the mesh (shifted mesh)
       if (mapibz(j) .eq. 0) then
          mapibz(j) = nibz
          w2(nibz) = w2(nibz) + 1
       elseif (mapibz(j) .ne. nibz) then
          write(ulog,*)'IBZ_REDUCE: mesh point ',j,' is in the stars of IBZ points ',mapibz(j),nibz
          stop 'IBZ_REDUCE: inconsistent stars, see the log'
       endif
       call fold_ws(q, qws)
       f = dot_product(uwedge, qws)
       if (f .gt. best+1d-8 .or. (abs(f-best) .le. 1d-8 .and. lexi_greater(qws, krep))) then
          best = f; jrep = j; krep = qws
       endif
    enddo lloop
    mapinv(nibz) = jrep
    k2(:,nibz) = krep
 enddo kploop
 if (allocated(kibz)) deallocate(kibz)
 if (allocated(wibz)) deallocate(wibz)
 allocate(kibz(3,nibz), wibz(nibz))
 kibz = k2(:,1:nibz); wibz = w2(1:nibz)
 sw = sum(wibz); wibz = wibz/sw
 write(ulog,*)'IBZ_REDUCE: ',nk,' mesh points, ',nrot,' Laue operations, ',nibz,' points in the irreducible wedge'
 deallocate(keys, ord, k2, w2)

 contains
   function mesh_key(qq) result(kk)
   real(8), intent(in) :: qq(3)
   integer(8) kk(3)
   integer m
   do m = 1, 3
      kk(m) = modulo(nint(dot_product(qq, rr(:,m))/(2*pi)*1d6, 8), 1000000_8)
   enddo
   end function mesh_key

   subroutine fold_ws(qq, w)
!! w = qq + G with the smallest |w|; on the zone boundary the image with the largest u.w (then (x,y,z))
   real(8), intent(in) :: qq(3)
   real(8), intent(out) :: w(3)
   real(8) r(3), v(3), b0(3), tol
   integer n1, n2, n3
   do n1 = 1, 3
      r(n1) = dot_product(qq, rr(:,n1))/(2*pi)
   enddo
   r = r - nint(r)
   b0 = matmul(gg, r); w = b0
   do n1 = -2, 2; do n2 = -2, 2; do n3 = -2, 2
      v = b0 + n1*gg(:,1) + n2*gg(:,2) + n3*gg(:,3)
      tol = 1d-6*max(sum(w*w), 1d-6)      ! lat_fc.dat gives the cell to ~1e-10: equal lengths differ by more
      if (sum(v*v) .lt. sum(w*w)-tol) then
         w = v
      elseif (abs(sum(v*v)-sum(w*w)) .le. tol) then
         if (dot_product(uwedge,v) .gt. dot_product(uwedge,w)+1d-8 .or.                            &
 &          (abs(dot_product(uwedge,v)-dot_product(uwedge,w)) .le. 1d-8 .and. lexi_greater(v,w))) w = v
      endif
   enddo; enddo; enddo
   end subroutine fold_ws
 end subroutine ibz_reduce
!------------------------------------------------------------
 subroutine reciprocal(rr, gg)
 real(8), intent(in) :: rr(3,3)
 real(8), intent(out) :: gg(3,3)
 real(8) v, pi
 pi = 4*atan(1d0)
 gg(:,1) = cross(rr(:,2), rr(:,3)); gg(:,2) = cross(rr(:,3), rr(:,1)); gg(:,3) = cross(rr(:,1), rr(:,2))
 v = dot_product(rr(:,1), gg(:,1))
 gg = 2*pi*gg/v
 contains
   function cross(a, b) result(c)
   real(8), intent(in) :: a(3), b(3)
   real(8) c(3)
   c = (/ a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1) /)
   end function cross
 end subroutine reciprocal
!------------------------------------------------------------
 logical function lexi_greater(a, b)
!! a > b in the order (x, then y, then z), with a tolerance
 real(8), intent(in) :: a(3), b(3)
 integer m
 lexi_greater = .false.
 do m = 1, 3
    if (a(m) .gt. b(m)+1d-8) then
       lexi_greater = .true.; return
    elseif (a(m) .lt. b(m)-1d-8) then
       return
    endif
 enddo
 end function lexi_greater
!------------------------------------------------------------
 logical function key_le(a, b)
 integer(8), intent(in) :: a(3), b(3)
 integer m
 key_le = .true.
 do m = 1, 3
    if (a(m) .lt. b(m)) return
    if (a(m) .gt. b(m)) then
       key_le = .false.; return
    endif
 enddo
 end function key_le
!------------------------------------------------------------
 subroutine sort_keys(n, k, ord)
!! ord: permutation sorting the keys k(:,i) (merge sort)
 integer, intent(in) :: n
 integer(8), intent(in) :: k(3,n)
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
          if (key_le(k(:,ord(i)), k(:,ord(j)))) then
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
 end subroutine sort_keys
!------------------------------------------------------------
 integer function find_key(n, k, ord, key)
!! index i with k(:,i) = key, 0 if none
 integer, intent(in) :: n, ord(n)
 integer(8), intent(in) :: k(3,n), key(3)
 integer lo, hi, mid
 find_key = 0
 lo = 1; hi = n
 do while (lo .le. hi)
    mid = (lo+hi)/2
    if (all(k(:,ord(mid)) .eq. key)) then
       find_key = ord(mid); return
    elseif (key_le(k(:,ord(mid)), key)) then
       lo = mid+1
    else
       hi = mid-1
    endif
 enddo
 end function find_key

 end module ibz_rules
