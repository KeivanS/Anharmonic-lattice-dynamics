!============================================================
 module point_group_rules
!! Point groups from their 3x3 Cartesian operations (2026-11-02, user): multiplication table, classes,
!! character table (computed, not looked up), reduction of any representation into irreducible ones.
!! Independent of phonons: the caller passes the operations and, for a representation, its character
!! for every operation (phonons: displacement patterns; electrons: Bloch states at a k-point with the
!! little group's operations; IR/Raman: the vector and symmetric-tensor characters below).
!!
!! Characters by the class-multiplication method: with classes K_1..K_r (sizes h_k) and
!! K_i K_j = sum_k c_ijk K_k, the vectors w_k = h_k chi(K_k)/d of an irrep (dimension d) satisfy
!! M_i w = w_i w, M_i(j,k) = c_ijk. In the metric diag(1/h) these M_i are normal, so the Hermitian and
!! anti-Hermitian parts of D^-1/2 M_i D^1/2 (D = diag(h)) share their eigenvectors; one random real
!! combination of them is diagonalized with zheev. chi_k = d w_k/h_k with w_1 = 1 (identity class) and
!! d^2 = |G| / sum_k |w_k|^2/h_k. Checked on return: orthogonality of rows and columns, sum d^2 = |G|.
!!
!! Labels (not Mulliken's numbering; user 2026-11-02: type, degeneracy and parity are what matter):
!! A/B for d = 1 (A: chi = +1 on every proper rotation of the highest order; groups without rotations,
!! and cubic groups (several principal axes): A), E for d = 2, T for d = 3; g/u from the inversion; '
!! and '' from sigma_h (mirror normal to the principal axis) when there is no inversion; then a serial
!! number among irreps with the same symbol, left out when the symbol occurs once.
!! Complex 1-d irreps come in conjugate pairs (C3, C4, C6, S4, S6, C3h, C4h, C6h, T, Th); for real
!! (time-reversal-invariant) problems such a pair is one 2-fold degenerate level: pg_partner(j) gives
!! the conjugate, and the label is 1E/2E-like ("Ea", "Eb").
!! Magnetic groups (antiunitary operations) and double groups (spin) are not treated: pass the unitary
!! operations; corepresentations would be built on top of this table.
 use constants, only : r15
 implicit none
 private
 public :: pg_setup, pg_reduce, pg_write_table, pg_op_name, pg_chi_vector, pg_chi_sym2
 integer, public :: pg_nop = 0, pg_ncls = 0, pg_nirr = 0
 integer, public, allocatable :: pg_mult(:,:), pg_cls(:), pg_hcls(:), pg_dim(:), pg_partner(:)
 real(r15), public, allocatable :: pg_ops(:,:,:)
 complex(r15), public, allocatable :: pg_chi(:,:)          ! (irrep, class)
 character(len=12), public, allocatable :: pg_clsname(:), pg_label(:)
 real(r15), parameter :: tol = 1d-6

 contains
!------------------------------------------------------------
 function det3(s) result(d)
 real(r15), intent(in) :: s(3,3)
 real(r15) d
 d = s(1,1)*(s(2,2)*s(3,3)-s(2,3)*s(3,2)) - s(1,2)*(s(2,1)*s(3,3)-s(2,3)*s(3,1)) &
 &  + s(1,3)*(s(2,1)*s(3,2)-s(2,2)*s(3,1))
 end function det3
!------------------------------------------------------------
 subroutine pg_op_name(s, name, order, axis)
!! type of a point operation: E, C2, C3, C4, C6, i, m (mirror), S3, S4, S6; order n of its proper
!! part R = det(s) s (rotation by 2 pi/n); axis = rotation axis of R (mirror: its normal)
 real(r15), intent(in) :: s(3,3)
 character(len=*), intent(out) :: name
 integer, intent(out) :: order
 real(r15), intent(out) :: axis(3)
 real(r15) r(3,3), d, cs, c1(3), c2(3), c3(3), pi
 integer j
 pi = 4*atan(1d0)
 d = sign(1d0, det3(s))
 r = d*s
 cs = max(-1d0, min(1d0, (r(1,1)+r(2,2)+r(3,3)-1)/2))
 if (cs .gt. 1-tol) then
    order = 1
 else
    order = nint(2*pi/acos(cs))
 endif
 axis = 0
 if (order .gt. 1) then                 ! null vector of R - 1
    do j = 1, 3 ; r(j,j) = r(j,j) - 1 ; enddo
    c1 = cross(r(1,:), r(2,:)) ; c2 = cross(r(2,:), r(3,:)) ; c3 = cross(r(3,:), r(1,:))
    axis = c1
    if (sum(c2**2) .gt. sum(axis**2)) axis = c2
    if (sum(c3**2) .gt. sum(axis**2)) axis = c3
    axis = axis/sqrt(sum(axis**2))
 endif
 if (d .gt. 0) then
    if (order .eq. 1) then
       name = 'E'
    else
       write(name,'(a,i0)') 'C', order
    endif
 else                                   ! s = -R = i R
    select case (order)
    case (1) ; name = 'i'
    case (2) ; name = 'm'
    case (3) ; name = 'S6'
    case (4) ; name = 'S4'
    case (6) ; name = 'S3'
    case default ; write(name,'(a,i0)') 'iC', order
    end select
 endif
 contains
   function cross(a, b) result(c)
   real(r15), intent(in) :: a(3), b(3)
   real(r15) c(3)
   c = (/ a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1) /)
   end function cross
 end subroutine pg_op_name
!------------------------------------------------------------
 function pg_chi_vector(s) result(c)
!! character of the vector representation (x, y, z: IR activity, electric-dipole transitions)
 real(r15), intent(in) :: s(3,3)
 real(r15) c
 c = s(1,1) + s(2,2) + s(3,3)
 end function pg_chi_vector
!------------------------------------------------------------
 function pg_chi_sym2(s) result(c)
!! character of the symmetric second-rank tensors (xx, xy, ...: Raman activity)
 real(r15), intent(in) :: s(3,3)
 real(r15) c
 c = ((s(1,1)+s(2,2)+s(3,3))**2 + sum(s*transpose(s)))/2
 end function pg_chi_sym2
!------------------------------------------------------------
 subroutine pg_setup(nop, ops, ulog, ier)
!! the group of the nop operations ops(3,3,nop) (Cartesian, orthogonal): multiplication table, classes
!! and character table, in the module variables. ier = 0 on success; 1: not closed; 2: table checks failed.
 integer, intent(in) :: nop, ulog
 real(r15), intent(in) :: ops(3,3,nop)
 integer, intent(out) :: ier
 integer a, b, c, k, i, j, l, ie, nr, info, ord, ordmax, ninv, hmax
 integer, allocatable :: inv(:), rep(:), cnt(:,:,:), key(:), perm(:), newc(:)
 integer kc2, i1, i2
 real(r15), allocatable :: hv(:), wr(:), rwork(:)
 complex(r15), allocatable :: m(:,:,:), h(:,:), work(:), w(:,:), chi(:,:)
 real(r15) s(3,3), ax(3), d2, dev, rnd, x, axp(3), cpr
 character(len=12) nm
 character(len=12), allocatable :: base(:)
 integer, allocatable :: ser(:)
 logical, allocatable :: done(:)

 ier = 0
 if (allocated(pg_ops)) deallocate(pg_ops, pg_mult, pg_cls, pg_hcls, pg_clsname, pg_chi, pg_dim, &
 &                                  pg_partner, pg_label)
 pg_nop = nop
 allocate(pg_ops(3,3,nop), pg_mult(nop,nop), pg_cls(nop), inv(nop))
 pg_ops = ops
! multiplication table: ops(a) ops(b) = ops(pg_mult(a,b)); identity ie
 ie = 0
 do a = 1, nop
    if (maxval(abs(ops(:,:,a) - reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/)))) .lt. tol) ie = a
    do b = 1, nop
       s = matmul(ops(:,:,a), ops(:,:,b))
       pg_mult(a,b) = 0
       do c = 1, nop
          if (maxval(abs(s - ops(:,:,c))) .lt. tol) then
             pg_mult(a,b) = c ; exit
          endif
       enddo
       if (pg_mult(a,b) .eq. 0) then
          write(ulog,*)' POINT_GROUP: the operations are not closed (product of ',a,b,')'
          ier = 1 ; return
       endif
    enddo
 enddo
 if (ie .eq. 0) then
    write(ulog,*)' POINT_GROUP: no identity among the operations'
    ier = 1 ; return
 endif
 do a = 1, nop
    do b = 1, nop
       if (pg_mult(a,b) .eq. ie) inv(a) = b
    enddo
 enddo
! classes: conjugates b a b^-1; class 1 = identity
 pg_cls = 0 ; nr = 0
 allocate(rep(nop))
 pg_cls(ie) = 1 ; nr = 1 ; rep(1) = ie
 do a = 1, nop
    if (pg_cls(a) .ne. 0) cycle
    nr = nr + 1 ; rep(nr) = a
    do b = 1, nop
       pg_cls(pg_mult(pg_mult(b,a),inv(b))) = nr
    enddo
 enddo
! fixed class order (independent of the order of the operations): E, C2, C3, C4, C6, i, m, S6, S4, S3;
! the C2 about and the mirror normal to the principal axis (highest-order rotation) first in their type
 ordmax = 1 ; axp = 0
 do a = 1, nop
    call pg_op_name(ops(:,:,a), nm, ord, ax)
    if (nm(1:1) .eq. 'C' .and. ord .gt. ordmax) then
       ordmax = ord ; axp = ax
    endif
 enddo
 allocate(key(nr), perm(nr), newc(nr))
 do k = 1, nr
    call pg_op_name(ops(:,:,rep(k)), nm, ord, ax)
    select case (trim(nm))
    case ('E') ; key(k) = 0
    case ('C2') ; key(k) = 10
    case ('C3') ; key(k) = 20
    case ('C4') ; key(k) = 30
    case ('C6') ; key(k) = 40
    case ('i') ; key(k) = 50
    case ('m') ; key(k) = 60
    case ('S6') ; key(k) = 70
    case ('S4') ; key(k) = 80
    case ('S3') ; key(k) = 90
    case default ; key(k) = 95
    end select
    if ((nm .eq. 'C2' .or. nm .eq. 'm') .and. abs(dot_product(ax, axp)) .lt. 1-1d-6) key(k) = key(k) + 1
    perm(k) = k
 enddo
 do i = 1, nr-1                                  ! stable sort of the classes by key
    do j = 1, nr-i
       if (key(perm(j+1)) .lt. key(perm(j))) then
          l = perm(j) ; perm(j) = perm(j+1) ; perm(j+1) = l
       endif
    enddo
 enddo
 do k = 1, nr
    newc(perm(k)) = k
 enddo
 pg_cls = newc(pg_cls)
 rep(1:nr) = rep(perm)
 deallocate(key, perm, newc)
 pg_ncls = nr
 allocate(pg_hcls(nr), pg_clsname(nr))
 do k = 1, nr
    pg_hcls(k) = count(pg_cls .eq. k)
    call pg_op_name(ops(:,:,rep(k)), nm, ord, ax)
    write(pg_clsname(k),'(i0,a)') pg_hcls(k), trim(nm)
    if (pg_hcls(k) .eq. 1) pg_clsname(k) = nm
 enddo
! class multiplication constants c_ijk: for a fixed z in K_k, the pairs x in K_i, y in K_j with xy = z
 allocate(cnt(nr,nr,nr), m(nr,nr,nr), h(nr,nr), hv(nr), wr(nr), work(4*nr), rwork(3*nr), w(nr,nr))
 cnt = 0
 do k = 1, nr
    do a = 1, nop
    do b = 1, nop
       if (pg_mult(a,b) .eq. rep(k)) cnt(pg_cls(a), pg_cls(b), k) = cnt(pg_cls(a), pg_cls(b), k) + 1
    enddo
    enddo
 enddo
! A_i = D^-1/2 M_i D^1/2 (normal); H = sum_i (r_i Herm(A_i) + s_i AntiHerm(A_i)), r, s fixed "random"
 hv = pg_hcls
 h = 0
 do i = 1, nr
    do j = 1, nr
    do k = 1, nr
       m(j,k,i) = cnt(i,j,k) * sqrt(hv(k)/hv(j))
    enddo
    enddo
    rnd = sin(1.7d0*i + 0.3d0)           ! fixed pseudo-random weights (reproducible)
    x = cos(2.9d0*i + 1.1d0)
    h = h + rnd*(m(:,:,i) + conjg(transpose(m(:,:,i))))/2 &
 &        + x*(m(:,:,i) - conjg(transpose(m(:,:,i))))/(2*(0d0,1d0))
 enddo
 call zheev('V', 'U', nr, h, nr, wr, work, 4*nr, rwork, info)
 if (info .ne. 0) then
    write(ulog,*)' POINT_GROUP: zheev info=',info
    ier = 2 ; return
 endif
! eigenvectors u -> w = D^1/2 u, w_1 = 1, chi_k = d w_k/h_k
 pg_nirr = nr
 allocate(pg_chi(nr,nr), pg_dim(nr), chi(nr,nr))
 do j = 1, nr
    w(:,j) = h(:,j) * sqrt(hv)
    w(:,j) = w(:,j)/w(1,j)
    d2 = nop / sum(abs(w(:,j))**2/hv)
    pg_dim(j) = nint(sqrt(d2))
    chi(j,:) = pg_dim(j) * w(:,j)/hv
 enddo
! order: by dimension, then real before complex, then by chi of the classes (stable)
 call sort_irreps(nr, chi, pg_dim)
 pg_chi = chi
! checks: sum d^2 = |G|; row and column orthogonality
 dev = abs(sum(pg_dim**2) - nop)
 do i = 1, nr
    do j = 1, nr
       dev = max(dev, abs(sum(hv*conjg(pg_chi(i,:))*pg_chi(j,:))/nop - merge(1d0, 0d0, i.eq.j)))
       dev = max(dev, abs(sum(conjg(pg_chi(:,i))*pg_chi(:,j))*hv(i)/nop - merge(1d0, 0d0, i.eq.j)))
    enddo
 enddo
 if (dev .gt. 1d-6) then
    write(ulog,*)' POINT_GROUP: character table checks failed, max deviation ',dev
    ier = 2
 endif
! conjugate partners (complex irreps)
 allocate(pg_partner(nr))
 do i = 1, nr
    pg_partner(i) = i
    if (maxval(abs(aimag(pg_chi(i,:)))) .lt. 1d-6) cycle
    do j = 1, nr
       if (maxval(abs(pg_chi(j,:) - conjg(pg_chi(i,:)))) .lt. 1d-6) pg_partner(i) = j
    enddo
 enddo
! labels
 ordmax = 1 ; ninv = 0 ; hmax = 0
 do k = 2, nr
    call pg_op_name(ops(:,:,rep(k)), nm, ord, ax)
    if (nm .eq. 'i') ninv = k
    if (nm(1:1) .eq. 'C') ordmax = max(ordmax, ord)
 enddo
 do k = 2, nr                                   ! size of the largest class of highest-order rotations
    call pg_op_name(ops(:,:,rep(k)), nm, ord, ax)
    if (nm(1:1) .eq. 'C' .and. ord .eq. ordmax) hmax = max(hmax, pg_hcls(k))
 enddo
 axp = 0
 do k = 2, nr                                   ! principal axis: a highest-order proper rotation
    call pg_op_name(ops(:,:,rep(k)), nm, ord, ax)
    if (nm(1:1) .eq. 'C' .and. ord .eq. ordmax) then
       axp = ax ; exit
    endif
 enddo
 allocate(pg_label(nr), base(nr), ser(nr), done(nr))
 do i = 1, nr
    select case (pg_dim(i))
    case (1)
       base(i) = 'A'
       do k = 2, nr
          call pg_op_name(ops(:,:,rep(k)), nm, ord, ax)
          if (nm(1:1) .eq. 'C' .and. ord .eq. ordmax .and. abs(pg_chi(i,k) - 1) .gt. 1d-6 &
 &            .and. hmax .le. 2) base(i) = 'B'          ! cubic (several principal axes): all A
       enddo
       if (pg_partner(i) .ne. i) base(i) = 'E'
    case (2) ; base(i) = 'E'
    case (3) ; base(i) = 'T'
    case default ; write(base(i),'(a,i0)') 'd', pg_dim(i)
    end select
    if (ninv .gt. 0) then
       if (real(pg_chi(i,ninv)) .gt. 0) then
          base(i) = trim(base(i))//'g'
       else
          base(i) = trim(base(i))//'u'
       endif
    else                                         ! sigma_h: mirror with normal along the principal axis
       do k = 2, nr
          call pg_op_name(ops(:,:,rep(k)), nm, ord, ax)
          cpr = abs(dot_product(ax, axp))
          if (nm .eq. 'm' .and. (cpr .gt. 1-1d-6 .or. ordmax .eq. 1)) then
             if (real(pg_chi(i,k))/pg_dim(i) .gt. 0) then
                base(i) = trim(base(i))//"'"
             else
                base(i) = trim(base(i))//"''"
             endif
             exit
          endif
       enddo
    endif
 enddo
 done = .false.
 do i = 1, nr
    if (done(i)) cycle
    l = 0
    do j = i, nr
       if (base(j) .eq. base(i) .and. .not. done(j)) then
          if (pg_partner(j) .ne. j .and. pg_partner(j) .lt. j) cycle   ! numbered with its partner
          l = l + 1 ; ser(j) = l ; done(j) = .true.
          if (pg_partner(j) .ne. j) then
             ser(pg_partner(j)) = l ; done(pg_partner(j)) = .true.
          endif
       endif
    enddo
 enddo
 if (ordmax .eq. 6) then          ! hexagonal: E1 has chi(C2 about the principal axis) = -2 (as the tables)
    kc2 = 0
    do k = 2, nr
       call pg_op_name(ops(:,:,rep(k)), nm, ord, ax)
       if (nm .eq. 'C2' .and. abs(dot_product(ax, axp)) .gt. 1-1d-6) kc2 = k
    enddo
    if (kc2 .gt. 0) then
       do i1 = 1, nr                 ! real 2-d irreps, and complex pairs through their first member
       if (pg_partner(i1) .lt. i1) cycle
       do i2 = i1+1, nr
          if (pg_partner(i2) .lt. i2) cycle
          if (base(i1) .eq. base(i2) .and. base(i1)(1:1) .eq. 'E' .and. pg_dim(i1) .eq. pg_dim(i2) .and. &
 &            real(pg_chi(i1,kc2)) .gt. 0 .and. real(pg_chi(i2,kc2)) .lt. 0) then
             l = ser(i1) ; ser(i1) = ser(i2) ; ser(i2) = l
             ser(pg_partner(i1)) = ser(i1) ; ser(pg_partner(i2)) = ser(i2)
          endif
       enddo
       enddo
    endif
 endif
 kc2 = 0                          ! cubic: T1 has chi = +1 on the C4 class (O, Oh) or the S4 class (Td)
 do k = 2, nr
    call pg_op_name(ops(:,:,rep(k)), nm, ord, ax)
    if (nm .eq. 'C4') kc2 = k
 enddo
 if (kc2 .eq. 0) then
    do k = 2, nr
       call pg_op_name(ops(:,:,rep(k)), nm, ord, ax)
       if (nm .eq. 'S4') kc2 = k
    enddo
 endif
 if (kc2 .gt. 0) then
    do i1 = 1, nr
    do i2 = i1+1, nr
       if (base(i1) .eq. base(i2) .and. pg_dim(i1) .eq. 3 .and. pg_dim(i2) .eq. 3 .and. &
 &         real(pg_chi(i1,kc2)) .lt. 0 .and. real(pg_chi(i2,kc2)) .gt. 0) then
          l = ser(i1) ; ser(i1) = ser(i2) ; ser(i2) = l
       endif
    enddo
    enddo
 endif
 do i = 1, nr
    if (count(base .eq. base(i) .and. ser .gt. 1) .eq. 0) then   ! the only one with this symbol: no number
       pg_label(i) = base(i)(1:1)//trim(base(i)(2:))
    else
       write(pg_label(i),'(a,i0,a)') base(i)(1:1), ser(i), trim(base(i)(2:))
    endif
    if (pg_partner(i) .ne. i) pg_label(i) = trim(pg_label(i))//merge('a','b',i.lt.pg_partner(i))
 enddo
 deallocate(inv, rep, cnt, m, h, hv, wr, work, rwork, w, chi, base, ser, done)
 end subroutine pg_setup
!------------------------------------------------------------
 subroutine sort_irreps(nr, chi, dim)
!! irreps ordered by dimension, real before complex, then by -Re chi over the classes (stable)
 integer, intent(in) :: nr
 complex(r15), intent(inout) :: chi(nr,nr)
 integer, intent(inout) :: dim(nr)
 integer i, j, k, it
 complex(r15) ct(nr)
 logical sw
 do i = 1, nr-1
    do j = 1, nr-i
       sw = .false.
       if (dim(j+1) .lt. dim(j)) then
          sw = .true.
       elseif (dim(j+1) .eq. dim(j)) then
          if (maxval(abs(aimag(chi(j+1,:)))) .lt. 1d-6 .and. maxval(abs(aimag(chi(j,:)))) .ge. 1d-6) then
             sw = .true.
          elseif ((maxval(abs(aimag(chi(j+1,:)))) .lt. 1d-6) .eqv. (maxval(abs(aimag(chi(j,:)))) .lt. 1d-6)) then
             do k = 2, nr
                if (abs(real(chi(j+1,k)) - real(chi(j,k))) .gt. 1d-6) then
                   sw = real(chi(j+1,k)) .gt. real(chi(j,k)) ; exit
                endif
                if (abs(aimag(chi(j+1,k)) - aimag(chi(j,k))) .gt. 1d-6) then
                   sw = aimag(chi(j+1,k)) .gt. aimag(chi(j,k)) ; exit
                endif
             enddo
          endif
       endif
       if (sw) then
          ct = chi(j,:) ; chi(j,:) = chi(j+1,:) ; chi(j+1,:) = ct
          it = dim(j) ; dim(j) = dim(j+1) ; dim(j+1) = it
       endif
    enddo
 enddo
 end subroutine sort_irreps
!------------------------------------------------------------
 subroutine pg_reduce(chiop, mult)
!! multiplicities of the irreps in a representation with character chiop(g) for every operation g
!! (in the order passed to pg_setup): n_j = (1/|G|) sum_g conjg(chi_j(g)) chiop(g)
 complex(r15), intent(in) :: chiop(pg_nop)
 real(r15), intent(out) :: mult(pg_nirr)
 integer j, g
 complex(r15) s
 do j = 1, pg_nirr
    s = 0
    do g = 1, pg_nop
       s = s + conjg(pg_chi(j, pg_cls(g))) * chiop(g)
    enddo
    mult(j) = real(s)/pg_nop
 enddo
 end subroutine pg_reduce
!------------------------------------------------------------
 subroutine pg_write_table(u, prefix)
!! the character table: classes (size + operation type), one row per irrep (complex values as re+im i)
 integer, intent(in) :: u
 character(len=*), intent(in) :: prefix
 integer i, k
 character(len=24) cell
 character(len=600) line
 write(line,'(a,a8)') prefix, 'irrep'
 do k = 1, pg_ncls
    write(cell,'(a12)') adjustr(pg_clsname(k))
    line = trim(line)//cell(1:12)
 enddo
 write(u,'(a)') trim(line)
 do i = 1, pg_nirr
    write(line,'(a,a8)') prefix, adjustr(pg_label(i)(1:8))
    do k = 1, pg_ncls
       if (abs(aimag(pg_chi(i,k))) .lt. 1d-6) then
          write(cell,'(f12.3)') real(pg_chi(i,k))
       else
          write(cell,'(f6.2,sp,f5.2,a)') real(pg_chi(i,k)), aimag(pg_chi(i,k)), 'i'
          cell = adjustr(cell(1:12))
       endif
       line = trim(line)//cell(1:12)
    enddo
    write(u,'(a)') trim(line)
 enddo
 end subroutine pg_write_table
 end module point_group_rules
