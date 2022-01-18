 !============================================================================
 ! functional programming: substitute for map(f, x)
 PURE ELEMENTAL REAL FUNCTION square(x)
    REAL, INTENT(in) :: x
    square = x**2
 END FUNCTION square
 !=============================================================================
 subroutine findword(word,line,found)
! says if the string "word" exists in string "line"
 implicit none
 logical found
 character(*) line, word
 integer i,l,k
 l = len_trim(line); k=len_trim(word)
 found = .False.
 do i=1,l-k+1
    if(line(i:i+k-1) .eq. word(1:k)) then
       found = .True.
       exit
    endif
 enddo
 end subroutine findword

!======================================================================
 subroutine make_r0g
! matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in the supercell
! matmul(r0g,v) gives the 3 reduced coordinates of an atom in the supercell if its
! reduced coordinates in the primitive cell are given by v
 use lattice
 implicit none

 r0g(1,1) = r01 .dot. g1
 r0g(2,1) = r01 .dot. g2
 r0g(3,1) = r01 .dot. g3
 r0g(1,2) = r02 .dot. g1
 r0g(2,2) = r02 .dot. g2
 r0g(3,2) = r02 .dot. g3
 r0g(1,3) = r03 .dot. g1
 r0g(2,3) = r03 .dot. g2
 r0g(3,3) = r03 .dot. g3

 end subroutine make_r0g
!============================================================
 subroutine check(r,a1,a2,a3,ier,g1,g2,g3)
! subroutine to check whether r is an integer multiple of (r01,r02,r03)
! output ai are the coefficients of its linear combination on this basis
! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use ios
 implicit none

 type(vector) :: r,g1,g2,g3
 real(8) a1,a2,a3,ep
 integer ier

 ep = tolerance ! 0.001 ! displacements larger than 0.001 A are not tolerated
!call write_out(ulog,'CHECK: R ',r)
 a1 = r .dot. g1
 a2 = r .dot. g2
 a3 = r .dot. g3
 if(abs(a1-nint(a1)).ge.ep .or.abs(a2-nint(a2)).ge.ep   &
&   .or.abs(a3-nint(a3)).ge.ep ) then
!  write(ulog,*) ' R is not a multiple of r0s , check your inputs '
!  write(ulog,*) ' R.r0i=',a1,a2,a3
   ier = 1
!  stop
 else
!  write(ulog,*) ' R is a multiple of r0s'
!  write(ulog,3) ' n1,n2,n3  =',a1,a2,a3
   ier = 0
 endif
3 format(a,3(1x,g12.6))

 end subroutine check

!======================================================================
  function fcut(x,x1,x2) result(y)
  implicit none
  real(8), intent(in) :: x,x1,x2
  real(8) y

  if (x.lt.x1) then
    y = 1
  elseif(x.gt.x2) then
    y = 0
  else
    y = 0.5*(1+cos(3.1415926*(x-x1)/(x2-x1)))
  endif

  end function fcut
!=====================================================================
 subroutine calculate_distance(i,j,atompos,maxatoms,rij)
 implicit none
 integer i,j,maxatoms
 real(8) atompos(3,maxatoms),rij

  rij = sqrt( (atompos(1,i)-atompos(1,j))*(atompos(1,i)-atompos(1,j)) +  &
&             (atompos(2,i)-atompos(2,j))*(atompos(2,i)-atompos(2,j)) +  &
&             (atompos(3,i)-atompos(3,j))*(atompos(3,i)-atompos(3,j)) )

 end subroutine calculate_distance
!============================================================
 subroutine sort(n,r,mcor,maxat)
! sorts the first n elements of array r in ascending order r(mcor(i)) is the ordered array
! make sure to initialize the array r to a huge number if n\=maxat
  implicit none
  integer maxat,n,i,j,temp
  real(8) r(maxat)
  integer mcor(maxat)

  do i=1,maxat ! n
     mcor(i)=i
  enddo
  do i=1,n  ! was 1,n-1
     do j=i+1,maxat  ! was i+1,n
        if(r(mcor(i)).gt.r(mcor(j))) then
           temp=mcor(i)
           mcor(i)=mcor(j)
           mcor(j)=temp
        endif
     enddo
  enddo
!  r(mcor) is now the ordered array!
 end subroutine sort
!============================================================
subroutine calculate_volume(r1,r2,r3,om)
use geometry
implicit none
real(8) om
type(vector) cross,r1,r2,r3

cross = r1 .cross. r2
om = abs(r3 .dot. cross)

end subroutine calculate_volume
!===============================================================
 subroutine make_reciprocal_lattice(r1,r2,r3,g1,g2,g3)
! be careful it does not have the factor of 2*pi in it!
 use geometry
 use constants
 use ios
 implicit none
 type(vector) :: r1,r2,r3,g1,g2,g3
 real(8) om

 om = r1 .dot. (r2 .cross. r3)
! write(ulog,*)'RECIPROCAL_LATTICE: '
 g1=(r2 .cross. r3)/om
 g2=(r3 .cross. r1)/om
 g3=(r1 .cross. r2)/om
! call write_out(ulog,'om ',om)
! call write_out(ulog,'g1 ',g1)
! call write_out(ulog,'g2 ',g2)
! call write_out(ulog,'g3 ',g3)
 end subroutine make_reciprocal_lattice
!===============================================================
 subroutine make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
! be careful it does not have the factor of 2*pi in it!
 use geometry
 use constants
 use ios
 implicit none
 type(vector) :: r1,r2,r3,g1,g2,g3
 real(8) om

 om = r1 .dot. (r2 .cross. r3)
 write(ulog,*)'RECIPROCAL_LATTICE: '
 g1=2*pi*(r2 .cross. r3)/om
 g2=2*pi*(r3 .cross. r1)/om
 g3=2*pi*(r1 .cross. r2)/om
 call write_out(ulog,'om ',om)
 call write_out(ulog,'g1 ',g1)
 call write_out(ulog,'g2 ',g2)
 call write_out(ulog,'g3 ',g3)

 end subroutine make_reciprocal_lattice_2pi
!==========================================================
 function delta(x) result(y)
 use constants
 implicit none
 real(8) x,y

 y = exp(-x*x/2d0)/sqrt(2*pi)

 end function delta
!==========================================================
 function delta_g(x,w) result(y)
 use constants
 implicit none
 real(8) x,y,w

 y = exp(-x*x/2d0/w/w)/sqrt(2*pi)/w

 end function delta_g
!==========================================================
 function delta_l(x,w) result(y)
 use constants
 implicit none
 real(8) x,y,w

 y = w/(x*x+w*w)/pi

 end function delta_l
!===========================================================
 function nbe(w,t,m) result(y)
 implicit none
 real(8) , intent(in) :: w,t
 real(8) y,z
 integer m

 if(t.le. 0) then
   write(*,*)'ERROR in N_BE: t=<0 ',t
   stop
 endif
 if(w.le. 0) then
   write(*,*)'ERROR in N_BE: w=<0 ',w
   stop
 endif

 z = w/t
!WRITE(*,*)'z=w/t=',z
 if(m.eq.0) then ! quantum calculation
   if (z.gt.60) then
     y = 0
   elseif(z.lt.0.0010) then
     y = 1/(z*(1+z/2*(1+z/3)))
     !y = 1/(z*(1+z/2*(1+z/3*(1+z/4))))
   else
     y = 1/(exp(z) - 1)
   endif
 elseif(m.eq.1) then ! classical calculation
   y=1/z
 else
   print*,'NBE: m must be either 0 or 1, not ',m
   stop
 endif

 end function nbe
!===========================================================
    subroutine find_root(fcn,x1,x2,root)
    implicit none
    real(8) :: x1,x2,root,fcn,error,fcnroot,fx1
    integer counter
    external fcn

    counter = 0
    error = 0.0000001
    fx1=fcn(x1)
    if(fx1*fcn(x2).gt.0) then
       print*,'bad guess for initial root brackets ',x1,x2,fx1
       stop
    endif
    fcnroot= 100000000

    do while(dabs(fcnroot).gt.error )
       counter = counter + 1
       if (counter.gt.50) then
          print*,'find_root did not converge after 50 iterations'
          print*,'guessed root,f(root)=',root,fcnroot
          stop
       endif
       root=(x1+x2)/2.
       fcnroot= fcn(root)
       if (fcnroot*fx1 .lt. 0 ) then
          x2= root
       else
          x1= root
          fx1 = fcnroot
       endif

    enddo

    end subroutine find_root
!===========================================================
 subroutine check_inside_bz(q,g1,g2,g3,inside)
 use geometry
! this one is in accordane with make_kp_reg with even nki, and zero shift
! between 0,gi before shifting by (g1+g2+g3)/2
 implicit none
 type(vector),intent(in) :: g1,g2,g3
 real(8),intent(in) :: q(3)
 real(8) :: qx,qy,qz,q2(3)
 integer n
 integer inside

! first shift q by  (g1+g2+g3)/2
 q2=q+ 0.5d0*(g1+g2+g3)

! then find its reduced coordinates
 call get_direct_components(q2,qx,qy,qz,g1,g2,g3,inside)

 end subroutine check_inside_bz
!===========================================================
 subroutine check_inside_fbz(q,g1,g2,g3,inside)
 use geometry
 use constants
 implicit none
 logical inside
 type(vector) g1,g2,g3
 real(8) qdotg,gg,gt(3),sm ,q(3)
 sm = -.000001  ! on the boundary is also counted as inside, but once only

 inside = .false.
!------------------- along the diagonal 100
       qdotg = ( q .dot. g1) + sm
       gg = g1 .dot. g1
       if(abs(qdotg) .gt. gg/2) return  ! this is the Bragg condition |q.G| < G.G/2
       qdotg = ( q .dot. g2) + sm
       gg = g2 .dot. g2
       if(abs(qdotg) .gt. gg/2) return
       qdotg = ( q .dot. g3) + sm
       gg = g3 .dot. g3
       if(abs(qdotg) .gt. gg/2) return
!------------------- along the diagonal 110
       gt = g1+g2
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1-g2
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1+g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1-g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return

       gt = g3+g2
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g3-g2
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
!------------------- along the diagonal 111
       gt = g1+g2+g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1+g2-g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1-g2+g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
       gt = g1-g2-g3
       qdotg = ( q .dot. gt) + sm
       gg = gt .dot. gt
       if(abs(qdotg) .gt. gg/2) return
!-------------------
 inside = .true.

 end subroutine check_inside_fbz
!===========================================================
 subroutine check_inside_irred_fbz2(k,inside)
! assuming k is inside FBZ, this checks whether it is inside the irreducible FBZ
! the case of FCC and hexagonal are treated here.
 use lattice
 use constants
 implicit none
 real(8) k(3),kdg,q2
 logical inside,inx

! this is for FCC
! if( (abs(k(3)).lt.1d-9 .or. (k(3) >= 0)) .and. (k(3) <= k(2)) .and. (k(2) <= k(1)) &
!&    .and. (k(1) <= boxg(1))  .and.  (k(1)+k(2)+k(3) <= 3.5d0*boxg(1)/2d0) ) then
!    inside = .true.  ! this goes 15% beyond the FBZ boundary
! else
!    inside = .false.
!    write(444,*)k
! endif

! this is for a 2d hexagonal lattice (ex: graphene)
 inx = .false. ; inside=.false.
 if ( abs(k(3)).lt.1d-3 ) then
    kdg = k(1)*g1%component(1) + k(2)*g1%component(2)
    q2 =  g1%component(1)*g1%component(1) + g1%component(2)*g1%component(2)
    if( kdg .lt. q2/2 .and. kdg .ge.0) inx=.true.
    kdg = k(1)*g2%component(1) + k(2)*g2%component(2)
    q2 =  g2%component(1)*g2%component(1) + g2%component(2)*g2%component(2)
    if ( kdg .lt. q2/2 .and. kdg .ge.0) inside=inx
!   write(6,9)'k,g01,g02=',k,g1,g2,kdg,q2
 endif
!if (inside) write(444,*)"k HEX=",k
9 format(a,99(1x,f6.3))
 end subroutine check_inside_irred_fbz2
!===========================================================
 subroutine check_inside_irred_fbz(k,inside)
! assuming kp is inside FBZ, this checks whether it is inside the irreducible FBZ of FCC
 use lattice
 use constants
 implicit none
 real(8) k(3)
 logical inside

 if( (abs(k(3)).lt.1d-9 .or. (k(3) >= 0)) .and. (k(3) <= k(2)) .and. (k(2) <= k(1)) &
&    .and. (k(1) <= boxg(1))  .and.  (k(1)+k(2)+k(3) <= 3d0*boxg(1)/2d0) ) then
    inside = .true.
 else
    inside = .false.
    write(444,*)k
 endif

 end subroutine check_inside_irred_fbz
!===========================================================
 subroutine select_IBZ(nbz,kbz,nib,kib)
! for a given kpoint mesh kbz within the FBZ, this routine stores
! those which are in the irreducible FBZ into kib
 use lattice
 implicit none
 integer nbz,i,nib
 real(8) kbz(3,nbz),q(3),kib(3,nbz)
 logical inside
! real(8), allocatable kib(:,:),ki(:,:)
!  allocate(ki(3,nkfine))

 nib = 0
 kib = 0
 do i=1,nbz
    q(:) = kbz(:,i)
    call check_inside_irred_fbz2(q,inside)
    if(inside) then
       nib = nib+1
       kib(:,nib) = q(:)
    endif
 enddo

! allocate(kib(3,ni))
! do i=1,ni
!    kib(:,i)=ki(:,i)
! enddo
! deallocate(ki)

 end subroutine select_IBZ
!============================================================
!============================================================
 subroutine get_k_info(q,N,nk,i,j,k,g1,g2,g3,inside)
! for a vector q(3) in the primitive cell of the reciprocal space, defined on
! a mesh N(1),N(2),N(3), this subroutine finds the three indices of q
! and its number nk based on the triple loop
! nk=0 do i=0,N(1)-1 ; do j=0,N(2)-1; do k=0,N(3)-1; nk=nk+1 ;q=i*G1/N1+j*G2/N2+k*G3/N3

 use geometry
 implicit none
 integer, intent(in):: N(3)
 real(8), intent(in):: q(3)
 type(vector), intent(in):: g1,g2,g3
 integer, intent (out) :: nk,i,j,k,inside
 integer indexg

! this was for when the kpoints were between 0 and N-1

    call get_components_g(q,N,i,j,k,g1,g2,g3,inside)
    nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
 !  nk= indexg(i,j,k,-N(1)/2,N(1)/2-1,-N(2)/2,N(2)/2-1,-N(3)/2,N(3)/2-1)
 !WRITE(*,*)'The value of nk=',nk
    return

! this is for when the k vectors are between -N/2 and N/2-1
  if (mod(N(1),2).eq.0 .and. mod(N(2),2).eq.0 .and. mod(N(3),2).eq.0 ) then
! use this if shifted by -0.5(g1+g2+g3)
WRITE(*,*)'<get_components_g_centered> is needed'
    call get_components_g_centered (q,N,i,j,k,g1,g2,g3,inside)
!   write(*,*)'for even ni, ijk=',i,j,k
  else   ! do the regulat thing
WRITE(*,*)'<get_components_g> is needed'
    call get_components_g(q,N,i,j,k,g1,g2,g3,inside)
  endif
  nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
! write(*,*)'the corresponding nk=',nk

 end subroutine get_k_info
!============================================================
subroutine get_k_info_cent(q,N,nk,i,j,k,g1,g2,g3,inside)
! for a vector q(3) in the primitive cell of the CENTERED-reciprocal space, defined on
! a mesh N(1),N(2),N(3), this subroutine finds the three indices of q
! and its number nk based on the triple loop
! nk=0 do i=0,N(1)-1 ; do j=0,N(2)-1; do k=0,N(3)-1; nk=nk+1 ;q=i*G1/N1+j*G2/N2+k*G3/N3

 use geometry
 implicit none
 integer, intent(in):: N(3)
 real(8), intent(in):: q(3)
 type(vector), intent(in):: g1,g2,g3
 integer, intent (out) :: nk,i,j,k,inside
 integer indexg

! this was for when the kpoints were between 0 and N-1
    call get_components_g_centered (q,N,i,j,k,g1,g2,g3,inside)
    nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
 !  nk= indexg(i,j,k,-N(1)/2,N(1)/2-1,-N(2)/2,N(2)/2-1,-N(3)/2,N(3)/2-1)
    return

! this is for when the k vectors are between -N/2 and N/2-1
  if (mod(N(1),2).eq.0 .and. mod(N(2),2).eq.0 .and. mod(N(3),2).eq.0 ) then
! use this if shifted by -0.5(g1+g2+g3)
    call get_components_g_centered (q,N,i,j,k,g1,g2,g3,inside)
!   write(*,*)'for even ni, ijk=',i,j,k
  else   ! do the regulat thing
    call get_components_g(q,N,i,j,k,g1,g2,g3,inside)
  endif
  nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
! write(*,*)'the corresponding nk=',nk

 end subroutine get_k_info_cent
!============================================================
 function indexg(i,j,k,nil,nih,njl,njh,nkl,nkh) result(n)
! finds the index n of the kpoint defined with 3 loop indices
! ijk going in general from nil to nih, njl to njh nkl to nkh
! n=0; do i1=nil,nih; do j=njl,njh; do k=nkl,nkh; n=n+1
 implicit none
 integer i,j,k,nil,nih,njl,njh,nkl,nkh,n

 n = (k-nkl+1) + (j-njl)*(nkh-nkl+1) + (i-nil)*(njh-njl+1)*(nkh-nkl+1)

 end function indexg
!============================================================
 function indexn(i,j,k,n1,n2,n3) result(n)
! finds the index n of the coarse kpoint defined with 3 loop indices
! ijk going from -ni to ni
 implicit none
 integer i,j,k,n1,n2,n3,n

 n = (k+n3+1) + (j+n2)*(2*n3+1) + (i+n1)*(2*n2+1)*(2*n3+1)

 end function indexn
!============================================================
 function index_reg(i,j,k,n1,n2,n3) result(n)
! finds the index n of the coarse kpoint defined with 3 loop indices
! ijk going from 1 to ni , i being the outer loop and k the inner one
 implicit none
 integer i,j,k,n1,n2,n3,n

 n = k + (j-1)*n3 + (i-1)*n2*n3

 end function index_reg
!===========================================================
 subroutine make_sorted_gs(g1,g2,g3,nshell,gg)
! from the basis (g1,g2,g3) generate the nshell shortest linear combinations including 0
 use geometry
! use io2
 implicit none
 integer nshell,i,j,k,n5,ik,ns
 type(vector) :: g1,g2,g3
 real(8) gg(3,nshell)
 real(8), allocatable :: g5(:,:),gs(:)
 integer, allocatable :: xmap(:)

 ns = (2*nshell+1)**(0.3333)
 n5 = (2*ns+1)**3
 allocate ( g5(3,n5),gs(n5),xmap(n5) )

 ik = 0
 do i=-ns,ns
 do j=-ns,ns
 do k=-ns,ns
    ik = ik+1
    g5(:,ik) = i*g1 + j*g2 + k*g3
    gs(ik) = length(g5(:,ik))
 enddo
 enddo
 enddo

! now sort g2 in ascending order
 call sort(n5,gs,xmap,n5)
 do i=1,nshell
    gg(:,i) = g5(:,xmap(i))
    write(30,3)i,gg(:,i),length(gg(:,i))
!WRITE(*,*)'For CHECK purpose',length(gg(:,i))
 enddo

3 format(i6,9(2x,g10.4))

 deallocate (g5,xmap,gs)
 end subroutine make_sorted_gs
!============================================================
! subroutine get_direct_components(q,n,qx,qy,qz,g1,g2,g3,inside)
 subroutine get_direct_components(q,qx,qy,qz,g1,g2,g3,inside)
! for a given q-vector, it finds its direct components assuming it was
! created as: q=qx*g1+qy*g2+qz*g3
! if the integer variable inside=1 then q is inside the primcell, defined
! by [0:g1[,[0:g2[,[0:g3[ ; if inside=0 it's outside
 use geometry
 implicit none
 real(8), intent(in):: q(3)
 type(vector),intent(in):: g1,g2,g3
 type(vector)rr1,rr2,rr3
! integer, intent(in):: n(3)
 integer, intent(out):: inside
 real(8), intent(out):: qx,qy,qz
 real(8) epsl

 epsl=5d-10

 call make_reciprocal_lattice(g1,g2,g3,rr1,rr2,rr3)  ! no factor of 2pi
 inside = 1

 qx = (q.dot.rr1) + epsl
 if (qx.lt.0 .or. qx.ge.1) inside=0
! write(*,7)'i,qdotr,aux=',i,qdr,aux

 qy = (q.dot.rr2) + epsl
 if (qy.lt.0 .or. qy.ge.1) inside=0

 qz = (q.dot.rr3) + epsl
 if (qz.lt.0 .or. qz.ge.1) inside=0

7 format(a,i5,9(1x,g11.5))
 end subroutine get_direct_components
!============================================================
 subroutine comp1(q,rr,n,i,insid)
! inside would be =1 if 0<q<1
 use geometry
 implicit none
 real(8) q(3),aux,qdr,epsl
 type(vector) rr
 integer n,i,insid

 epsl=5d-10
!WRITE(*,*)'Check the rr value', rr%component !should they be 0?
!WRITE(*,*)'Check the q value', q
 qdr = (q.dot.rr) + epsl
! include from -0.00005 to 0.99995 included (include 0.0 and exclude 1.0)
 aux = qdr-floor(qdr)
 if (qdr.lt.0 .or. qdr.ge.1) insid = 0
 !if (aux .gt. 1e-4)     insid = 2
  !floor() always round down, thus make a difference with int() when it's negative

 i = nint(1+ n*aux) !nint() rounds its argument to the nearest integer
! write(*,7)'i,qdotr,aux=',i,qdr,aux
 if (i.eq.n+1) i=1

7 format(a,i5,9(1x,g11.5))
 end subroutine comp1
!============================================================
 subroutine comp_c(q,rr,n,i,insid)
! inside would be =1 if -0.5<q<0.5
 use geometry
 implicit none
 real(8) q(3),aux,qdr,epsl
 type(vector) rr
 integer n,i,insid

 qdr = (q.dot.rr)
! include from -0.00005 to 0.99995 included (include 0.0 and exclude 1.0)
 if (qdr.lt.-0.5d0 .or. qdr.ge.0.5d0) insid=0

! in order to get i, bring qdr between 0 and 1 by adding 0.5
 aux = 0.5d0+qdr-floor(qdr+0.5d0)
 i = nint(1+ n*aux)
! write(*,7)'i,qdotr,aux=',i,qdr,aux
 if (i.eq.n+1) i=1

7 format(a,i5,9(1x,g11.5))
 end subroutine comp_c
!============================================================
 subroutine get_components_g(q,n,i,j,k,gg1,gg2,gg3,inside)
! for a given q-vector, it finds its integer components assuming it was
! created as: q=(i-1)/n1*g1 + (j-1)/n2*g2 + (k-1)/n3*g3; i=1,n1 j=1,n2 k=1,n3
! as the ones created in make_kp_reg with zero shift
! if the integer variable inside=1 then q is inside the primcell
! if q is outside the prim_cell, then inside=0 but produces the i,j,k of its
! image inside
 use lattice
 use geometry
 implicit none
 real(8), intent(in):: q(3)
 type(vector),intent(in):: gg1,gg2,gg3!NO PROBLEM MARK
! type(vector)rx1,rx2,rx3
 integer, intent(in):: n(3)
 integer, intent(out):: i,j,k,inside

!These below are originally commented out, but in that case if <make_reciprocal_lattice>
!hasn't been called once before, the rr1,rr2,rr3 are all 0

! call make_reciprocal_lattice(g1,g2,g3,rr1,rr2,rr3)
! i = nint(1+ n(1)* q.dot.r1)
! j = nint(1+ n(2)* q.dot.r2)
! k = nint(1+ n(3)* q.dot.r3)
! if q.dot.r is larger than 1, we want to have its fractional part
 inside = 1

 call comp1(q,rr1,n(1),i,inside)
 call comp1(q,rr2,n(2),j,inside)
 call comp1(q,rr3,n(3),k,inside)

 end subroutine get_components_g
!============================================================
 subroutine get_components_g_centered(q,n,i,j,k,gg1,gg2,gg3,inside)
! for a given q-vector, it finds its integer components assuming it was
! created as: q=(i-1-n1/2)/n1*g1 + (j-1-n2/2)/n2*g2 + (k-1-n3/2)/n3*g3;
! i=1,n1 j=1,n2 k=1,n3
! as the ones created in make_kp_reg with zero shift
! if the integer variable inside=1 then q is inside the primcell, if inside=0 it's outside
! it works even if there's a shift less than 0.5, and if q is oustide the prim_cell
 use geometry
 use lattice
 implicit none
 real(8), intent(in):: q(3)
 type(vector),intent(in):: gg1,gg2,gg3
 integer, intent(in):: n(3)
 integer, intent(out):: i,j,k,inside
! real(8) w(3)

! w = q+0.5d0*(gg1+gg2+gg3)
! write(*,8)'q=',q
! write(*,8)'w=',w

 inside = 1

 call comp_c(q,rr1,n(1),i,inside)
 call comp_c(q,rr2,n(2),j,inside)
 call comp_c(q,rr3,n(3),k,inside)

! write(*,*)'w-components,inside are=',i,j,k ,inside

 end subroutine get_components_g_centered
!===========================================================
 subroutine set_neighbork(nkp,kp,nb,nbk,nbmax,kcutoff,nshl,gg)
! finds the list of neighbors to the nkp vectors kp within kcutoff but
! a maximum of nbmax per vector. Output is nb(i,l) which is the
! label of the lth neighbor of the ith vector, and nbk(i)( nbmax),
! which is the total no of neighbors of ith vector within kcutoff
! periodic structure for the kmesh is assumed : supercell is (g1,g2,g3)
! if no pbc should be applied, set nshl = 1
 use geometry
 use lattice
 use ios
 implicit none
 integer nkp,i,j,nbmax,l,sh,nshl,sh0
 integer nb(nkp,nbmax),nbk(nkp)
 real(8) kp(3,nkp),q2(3),qw(3),gg(3,nshl),dismin,kcutoff,kdist

 write(ulog,*)'SET_NEIGHBORK: kp#, neighbr#, kpnghbr#, kp,kdist,deltak,sh0'
 nb = 0
 do i=1,nkp
    l = 0
    jloop: do j=1,nkp
       if (i.eq.j) cycle jloop
       dismin = 1d9
       do sh=1,nshl
          q2 = kp(:,i)-kp(:,j)+gg(:,sh)
          kdist = length(q2)
          if(kdist.lt.dismin) then
             dismin = kdist
             qw = q2
             sh0 = sh
          endif
       enddo
! this has selected among kp(j) and its images the nearest one to kp(i)

       if (dismin .lt. kcutoff ) then
        if ( l.lt. nbmax) then
          l = l+1
          nb(i,l)=j
          write(ulog,6)i,l,nb(i,l),kp(:,i),dismin,qw,sh0
        else
!          write(ulog,*)'neighbor # too small, increase nbmax from ',nbmax
        endif
       else
!          write(ulog,*)'neighbor # too small, or decrease kcutoff from ',kcutoff
!          stop
       endif
    enddo jloop
    nbk(i) = l
 enddo

! nbmax=max(l,1)  ! in case l = 0 (one kpoint)
! if (nbmax .eq.1) then
!    nb(1,1)=1
! endif

 write(ulog,*)'SET_NEIGHBORK: largest no of neighbors=nbkmax,nshl=',nbmax,nshl
6 format(3(i8),7(2x,f9.4),1x,i4)

 end subroutine set_neighbork
!===========================================================
 subroutine send_to_fbz(kp,nshl,gg,q,ns)
! folds the kpoint kp into the FBZ and stores the result in q
! also returns ns : ns=1 means it was inside, ns=16 means too far , it failed
 use lattice
!! use io2
 implicit none
 integer, intent(in) :: nshl
 real(8), intent(in) :: kp(3),gg(3,nshl)
 integer, intent(out) :: ns
 real(8), intent(out) :: q(3)
 logical inside

 foldloop: do ns=1,nshl
     q = kp(:)+gg(:,ns)
     call check_inside_fbz(q,g1,g2,g3,inside)
     if(inside) then
       return
     endif
 enddo foldloop

3 format(3(i6),9(2x,f8.4))

 end subroutine send_to_fbz
!===========================================================
 subroutine send_to_primcell(kp,q)
! folds the kpoint kp into the primitive cell (between 0 and G) stores the result in q
 use lattice
!! use io2
 implicit none
 real(8), intent(in) :: kp(3)
 real(8), intent(out) :: q(3)
 real(8) a1,a2,a3

 a1 = (kp .dot. r1)/2/pi
 a2 = (kp .dot. r2)/2/pi
 a3 = (kp .dot. r3)/2/pi
 a1 = a1-floor(a1)
 a2 = a2-floor(a2)
 a3 = a3-floor(a3)
 q=a1*g1+a2*g2+a3*g3

 end subroutine send_to_primcell
!============================================================
 subroutine inverse_real(a,b,n)
      implicit none
      integer n,imax,k,j,i,ii
      real(8) ::  amax
      real(8) a(n,n),b(n,n),atmp,btmp,amult,div

      b=0d0
      do i=1,n
         b(i,i)=1d0
      enddo

  MAIN_LOOP: do k=1,n

         amax=-1.d-31
         imax=1
         if (k .ne. n) then

           do i=k,n
             if (abs(a(i,k)) .ge. amax) then
               amax=abs(a(i,k))
               imax=i
             endif
           enddo
           if (imax.ne.k) then
             do j=1,n
               atmp=a(imax,j)
               a(imax,j)=a(k,j)
               a(k,j)=atmp
               btmp=b(imax,j)
               b(imax,j)=b(k,j)
               b(k,j)=btmp
             enddo
           endif

         endif

         div=a(k,k)

         if(abs(div).gt.1.d-14) then

           do ii=1,n
             a(k,ii)=a(k,ii)/div
             b(k,ii)=b(k,ii)/div
           enddo
!          a(k,:)=a(k,:)/div
!          b(k,:)=b(k,:)/div

         endif

         do i=1,n
            amult=a(i,k)
            if (i.eq.k) cycle
            do j=1,n
               a(i,j)=a(i,j)-amult*a(k,j)
               b(i,j)=b(i,j)-amult*b(k,j)
            enddo
         enddo

   enddo MAIN_LOOP

 end subroutine inverse_real
!===============================================================
  subroutine mean_sd(x,mean,sd)
  implicit none
  integer n,i
  real(8) x(:),mean,sd,sd2

  n=size(x)
  mean = sum(x)/n
  sd2  = sum(x*x)/n
  sd   = sqrt(sd2 - mean*mean)

  end subroutine mean_sd
!=====================================================
  subroutine histogram(m,x,mesh)
! calculates the histogram of the data set x (needs not be sorted)
! and writes the distribution function in the file 'histo.dat'
  implicit none
  integer, intent(in):: m,mesh
  real(8), intent(in):: x(m)
  integer i,j,cnt,unit
  real(8) xmin,xmax,dx,sume,cume,e(mesh)

  unit=123
  open(unit,file='histo.dat')
  write(unit,*)'# j,  x(j) , e(j) , e(j)/sume/dx,accumulated e '

  cnt = size(x)
  xmax= maxval(x)
  xmin= minval(x)
  dx  = (xmax-xmin)/mesh
  write(*,5)'HISTO: size,xmin,xmax,dx=',cnt,xmin,xmax,dx

  e=0
  do i=1,cnt
! j=position of the bin to which x(i) belongs
     j= int((x(i)-xmin)/dx)+1
     if(j.eq.mesh+1) j=mesh
     if (j.lt.1 .or. j.gt.mesh) then
        write(*,*)'j is out of range ',j
        stop
     endif
!     write(unit,4)j,xmin+(j-0.5)*dx,e(j),e(j)/sume/dx,cume
  enddo
 4 format(i8,9(2x,g12.6))
 5 format(a,i8,9(2x,g12.6))

  close(unit)

  end subroutine histogram
!============================================================
 function trace(n,mat)
 implicit none
 integer n,i
 real(8) trace,mat(n,n)

 trace=0
 do i=1,n
    trace=trace+mat(i,i)
 enddo

 end function trace
!===========================================================
 subroutine dos_l(ene,nk,eig,afunc,dosl,adosl)
 implicit none
 integer nk,i,j
 real(8) eig(nk),afunc(nk),dosl,adosl,ene,width,delta_l

 width=0.03

 dosl=0; adosl=0
 do i=1,nk
    dosl = dosl+delta_l(ene-eig(i),width)/nk
    adosl=adosl+delta_l(ene-eig(i),width)*afunc(i)/nk
 enddo

 end subroutine dos_l
!===================================================================
 subroutine check_mdyn(ndn,nk,kp,eivec,eival)
! form the dynamical matrix from existing eigenvectors and eigenvalues
! and diagonalize and verify consistency :eivec(-q)=-conjg(eivec(q))
 use lattice  ! needed to access NC,g1,g2,g3
 use MatrixDiagonalize
 implicit none
 integer ndn,nk,j,k,l,i,ier,i1,j1,k1,mi,inside
 real(8) kp(3,nk),eival(ndn,nk),eivl(ndn)
 complex(8) eivec(ndn,ndn,nk),dynm(ndn,ndn),eivc(ndn,ndn),d2(ndn,ndn)

 do i=1,nk
    dynm=0
    do j=1,ndn
    do k=1,ndn
      do l=1,ndn
        dynm(j,k)=dynm(j,k)+eival(l,i)*eival(l,i)*cnst*cnst*  &
&                           eivec(j,l,i)*conjg(eivec(k,l,i))
      enddo
    enddo
    enddo
    call get_k_info(-kp(:,i),NC,mi,i1,j1,k1,g1,g2,g3,inside)
    d2=0
    do j=1,ndn
    do k=1,ndn
      do l=1,ndn
        d2(j,k)=d2(j,k)+eival(l,mi)*eival(l,mi)*cnst*cnst*     &
&                            eivec(j,l,mi)*conjg(eivec(k,l,mi))
      enddo
    enddo
    enddo

    do j=1,ndn
    do k=1,ndn
       if(abs( d2(j,k)-conjg(dynm(j,k))) .gt.1d-5) then
         write(*,4)'CHECK_MDYN: j,k,d(-q),d(q)=',j,k,d2(j,k),dynm(j,k)
       endif
    enddo
    enddo
!   return

    call diagonalize(ndn,dynm,eivl,ndn,eivc,ier)
    do j=1,ndn
       if( abs(cnst*sqrt(abs(eivl(j)))-eival(j,i)).gt.1d-4) then
          write(*,3)'CHECK_MDYN:j,eivl,eival=',j,cnst*sqrt(abs(eivl(j))),eival(j,i)
       endif
       do k=1,ndn
         if(abs(eivc(k,j)-eivec(k,j,i)).gt.1d-4) then
         if(abs(eivc(k,j)+eivec(k,j,i)).gt.1d-4) then
            write(*,4)'CHECK_MDYN:j,k,eiv(j),eivecs(k,j)=',j,k,eival(j,i),eivc(k,j),eivec(k,j,i)
         endif
         endif
       enddo
    enddo
 enddo
3 format(a,i6,9(1x,f15.6))
4 format(a,2i6,f15.7,9(1x,f11.4))

 end subroutine check_mdyn
!===========================================================
 subroutine set_dynamical_matrix(kpt,dynmat,ndim,ddyn) !I wrote my own code for this
 !use params
 !use lattice
 use kpoints
 use atoms_force_constants
 !use geometry
 use svd_stuff
 implicit none
 integer, intent(in) :: ndim
 complex(8), intent(out) :: dynmat(ndim,ndim),ddyn(ndim,ndim,3)
 real(8), intent(in) :: kpt(3)
 complex(8) junk
 real(8) mi,mj,all,rr(3),delt(3)
 integer i0,j,j0,al,be,i3,j3,t,u_log !,ired

4 format(a,4(1x,i5),9(2x,f9.4))
9 format(i6,99(1x,f9.3))

!---------------replace 'ulog' file------------------------
u_log=44
open(u_log,file='dynamic matrix log.dat',status='unknown',action='write')
!----------------------------------------------------------

 delt = wshift !*wshift
 ddyn   = cmplx(0d0,0d0)
 dynmat = cmplx(0d0,0d0)
 do i0=1,natoms0 !natoms0 should be the atom numbers within a unit cell
 do al=1,3
    i3 = al+3*(i0-1) !absolute order number of the term  with al direction, i0 atom type
    mi = atom0(i0)%mass !the mass of i0th atom
! write(ulog,*) 'i,al,mass=',i0,al,i3,mi
    tloop: do t=1,nterms(2) !nterms should be the total number(size) of force constant fc2
       if ( i0 .eq. iatomterm_2(1,t) .and. al .eq. ixyzterm_2(1,t) ) then !find the corresponding phi with same index: atom i0,direction al
          be = ixyzterm_2(2,t)   !al for the second atom
          j  = iatomterm_2(2,t)  !label number of the second atom
          j0 = iatomcell0(j)     !corresponding atom number in the first unit-cell
          mj = atom0(j0)%mass
          j3 = be+3*(j0-1)       !
          rr = atompos(:,j)-atompos(:,j0)
!          ired = igroup_2(t)
!          junk = fcs_2(ired)*ampterm_2(t)* exp(ci*(kpt .dot. rr))/sqrt(mi*mj)
          junk = fcs_2(t)* exp(ci*(kpt .dot. rr))/sqrt(mi*mj)
!  write(ulog,4) 't,j,be,mass,dyn=',t,j,be,j3,mj,junk
          dynmat(i3,j3) = dynmat(i3,j3) + junk
          ddyn(i3,j3,:) = ddyn(i3,j3,:) + junk*ci*rr(:)
       endif
    enddo tloop
! this leaves  gamma eivals unchanged
!   dynmat(i3,i3) = dynmat(i3,i3) + delt(al)*(1-cos(kpt .dot. r1))
! but this one shifts everything up by delt=wshift
    dynmat(i3,i3) = dynmat(i3,i3) + delt(al)
 enddo
 enddo

 if (verbose) then
  write(u_log,*)'SET_DYNAMICAL_MATRIX: d(dyn)/dk is:'
  do t=1,3
     write(u_log,*)'=======component of v ',t
     do j=1,ndim
        write(u_log,9) j, ddyn(j,:,t)
     enddo
  enddo
 endif

5 format(a,5i5,9(f8.3))
 all = sum(cdabs(dynmat(:,:)))/(ndim*ndim)
! make sure it is hermitian
 do t=1,ndim
    if (abs(aimag(dynmat(t,t))) .gt. 9d-4*abs(real(dynmat(t,t))) ) then
       write(u_log,*)' dynmat is not hermitian on its diagonal'
       write(u_log,*)' diagonal element i=',t,dynmat(t,t)
!      stop
!   else
       dynmat(t,t) = cmplx(real(dynmat(t,t)),0d0)
    endif
  do j=t+1,ndim-1
    if (abs(aimag(dynmat(t,j))+aimag(dynmat(j,t))) .gt. 9d-4*all ) then
       write(u_log,*)' dynmat is not hermitian in AIMAG of its off-diagonal elts'
       write(u_log,*)' off-diagonal element i,j=',t,j,dynmat(t,j),dynmat(j,t)
       write(u_log,*)' comparing to avg(abs(dynmat))=',all
!      stop
    elseif(abs(real(dynmat(t,j))-real(dynmat(j,t))) .gt. 9d-4*all ) then
       write(u_log,*)' dynmat is not hermitian in REAL of its off-diagonal elts'
       write(u_log,*)' off-diagonal element i,j=',t,j,dynmat(t,j),dynmat(j,t)
       write(u_log,*)' comparing to avg(abs(dynmat))=',all
!      stop
    else
    endif
! enforcing it to be hermitian!!
    mi=(aimag(dynmat(t,j))-aimag(dynmat(j,t)))/2
    mj=(real (dynmat(t,j))+real (dynmat(j,t)))/2
    dynmat(t,j) = cmplx(mj, mi)
    dynmat(j,t) = cmplx(mj,-mi)
  enddo
 enddo

! enforce it to be real if all imaginary componenents are very small
! all = sum(abs(aimag(dynmat(:,:))))/(ndim*ndim)
! if (all .lt. sum(abs(real(dynmat(:,:))))/(ndim*ndim)*1d-10) then
!  do t=1,ndim
!  do j=1,ndim
!     dynmat(t,j)=(dynmat(t,j)+conjg(dynmat(t,j)))/2d0
!  enddo
!  enddo
! endif

 end subroutine set_dynamical_matrix
!===========================================================
 subroutine nonanal(q,dynmat,ndim,ddyn)
 use constants
 use lattice
 use atoms_force_constants
 use born

 integer, intent(in) ::ndim
 complex(8), intent(inout) :: dynmat(ndim,ndim),ddyn(ndim,ndim,3)
 real(8), intent(in) :: q(3)
 real(8) zag(3),zbg(3),qeq,om0,ma,mb,rr(3),dqeq(3)
 integer na,nb,i,j,al
 real(8) eps_scale,q2,gg,term

 eps_scale=8.8541878176D-12/1D10/ee
 gg=(length(g1)*length(g2)*length(g3))**0.33
! if ( gg .lt. 4*rho) then
 if (born_flag.eq.0)  rho = gg/4d0
! write(30,*)'ADOPTED VALUE of RHO = ',rho
! endif
 rho2=rho*rho

 call calculate_volume(r1,r2,r3,om0)

 q2=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
 qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+ &
&       q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+ &
&       q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
 dqeq(1) = 2*epsil(1,1)*q(1)+(epsil(1,2)+epsil(2,1))*q(2)+(epsil(1,3)+epsil(3,1))*q(3)
 dqeq(2) = 2*epsil(2,2)*q(2)+(epsil(2,3)+epsil(3,2))*q(3)+(epsil(2,1)+epsil(1,2))*q(1)
 dqeq(3) = 2*epsil(3,3)*q(3)+(epsil(3,1)+epsil(1,3))*q(1)+(epsil(3,2)+epsil(2,3))*q(2)

 term=1d0/(qeq*eps_scale)/om0*exp(-q2/rho2)

 do na = 1,natoms0
   ma = atom0(na)%mass
 do nb = 1,natoms0
   mb = atom0(nb)%mass
   rr = atompos(:,na)-atompos(:,nb)
   do i=1,3
      zag(i) = q(1)*zeu(1,i,na)+q(2)*zeu(2,i,na)+q(3)*zeu(3,i,na)
      zbg(i) = q(1)*zeu(1,i,nb)+q(2)*zeu(2,i,nb)+q(3)*zeu(3,i,nb)
   end do
   do i = 1,3
   do j = 1,3
      dynmat(i+3*(na-1),j+3*(nb-1)) = dynmat(i+3*(na-1),j+3*(nb-1))+ &
      & zag(i)*zbg(j)*term/sqrt(ma*mb)

     do al=1,3
        ddyn(i+3*(na-1),j+3*(nb-1),al) = ddyn(i+3*(na-1),j+3*(nb-1),al)+term/sqrt(ma*mb) *  &
 &          (zeu(al,i,na)*zbg(j)+zeu(al,j,nb)*zag(i)-dqeq(al)*zag(i)*zbg(j)/qeq)
     enddo

   end do
   end do
 end do
 end do

 end subroutine nonanal
!=========================================================
!--------------------------------------------------------
 subroutine get_freq(kp,ndn,vg,eival,eivec)
! for given kp vector it calculates the dynamical matrix, eigenvalues and eigenvectors and
! corresponding group velocity using the Hellman-Feynman theorem.
 use params
 use io2
 use constants
 use born
 use MatrixDiagonalize
 implicit none
 integer, intent(in) :: ndn
 real(8), intent(in) :: kp(3)
 real(8), intent(out):: eival(ndn)
 complex(8), intent(out) :: eivec(ndn,ndn)
 integer i,j,k,l,ier,nd2,al
 integer, allocatable :: mp(:)
 real(8), allocatable :: eivl(:)
 real(8) absvec,vg(3,ndn)
 complex(8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:),temp(:,:)
 real(8) khat(3)

 nd2 = min(ndn,12)
 allocate(dynmat(ndn,ndn),mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3),temp(ndn,ndn))


    call set_dynamical_matrix(kp,dynmat,ndn,ddyn)

!JS: call nonanalytical term for born effective change
    if (kp(1)==0.d0 .AND. kp(2)==0.d0 .AND. kp(3)==0.d0) then
        khat=kp(:)+1.0D-10
    else
        khat=kp(:)
    endif
    call nonanal(khat,dynmat,ndn,ddyn)

    if (verbose) then
       write(ulog,3)' ======================================================================'
       write(ulog,3)' THE DYNAMICAL MATRIX for KP=',ndn,kp(:)
       do l=1,ndn
          write(ulog,8)(dynmat(l,j),j=1,nd2)
       enddo
       write(ulog,3)' ======================================================================'
       do al=1,3
          write(ulog,3)' THE k-DERIVATIVE OF THE DYNAMICAL MATRIX =',al
          do l=1,ndn
             write(ulog,8)(ddyn(l,j,al),j=1,nd2)
          enddo
       enddo
    endif

! store dynmat in temp
    temp=dynmat
    call diagonalize(ndn,dynmat,eivl,ndn,eivc,ier)

! sort eivals in ascending order
    call sort(ndn,eivl,mp,ndn)
!   write(ulog,6)'map=',mp

    if (ier.ne.0 .or. verbose) then
      write(   *,*)' ier=',ier
      write(ulog,*)' ier=',ier, ' EIGENVALUES ARE...'
      write(ulog,4) eivl
    endif

    do j=1,ndn
       eival(j) = eivl(mp(j))  ! so that all frequencies are positive
       do l=1,ndn
          eivec(l,j) = eivc(l,mp(j))
       enddo
    enddo

! if pure imaginary make eigenvectors real as they can be multiplied by any arbitrary number
    do l=1,ndn
       absvec = sum(abs(real(eivec(:,l))))
       if (absvec .lt. 1d-3) then
          eivec(:,l)=cmplx(0,1)*eivec(:,l)
       endif
    enddo

  if (verbose) then

    do l=1,ndn
        write(ulog,3)' GET_FREQ:-----------  BAND ',l,eival(l)
        do j=1,ndn/3
           write(ulog,5)'atom ',j,eivec(3*(j-1)+1,l),eivec(3*(j-1)+2,l),eivec(3*(j-1)+3,l)
        enddo
    enddo

! now calculate the square of frequencies based on dynmat
     dynmat=temp
     temp = matmul(dynmat,eivc)
     temp = matmul(transpose(conjg(eivc)),temp)
     write(ulog,*)' TEST: below square of eivals should appear in the diagonal if e.D.e is right'
     do j=1,ndn
       write(ulog,9)'e.D.e=',j, temp(j,:)
     enddo
!    write(ulog,*)' three components of ddyn are '
!    do al=1,3
!    do j=1,ndn
!      write(ulog,9)'al=',al, ddyn(j,:,al)
!    enddo
!    enddo

  endif

! now calculate the group velocities from Hellman-Feynmann formula(for each component al=1,3)
    vg=0
    do al=1,3

! first sandwich ddyn between eigenvectors
!      temp=matmul(ddyn(:,:,al),eivc)
!      dynmat=matmul(transpose(conjg(eivc)),temp)
!      do l=1,ndn
!         write(ulog,9)'al,d(Dyn)/dk=',al, dynmat(l,:)
!         vg(al,l)=dynmat(mp(l),mp(l))/(2*sqrt(abs(eival(l))))*cnst*1d-10*100*2*pi !*c_light
!      enddo


    do l=1,ndn
      vg(al,l)=0
      do k=1,ndn
         dynmat(k,l)=0
         do j=1,ndn
            dynmat(k,l)=dynmat(k,l)+ddyn(k,j,al)*eivec(j,l)
         enddo
      enddo
      do k=1,ndn
         vg(al,l)=vg(al,l)+dynmat(k,l)*conjg(eivec(k,l))
      enddo
      vg(al,l)=vg(al,l)/2/sqrt(abs(eival(l)))*cnst*1d-10*100*2*pi
    enddo
    enddo

 deallocate(eivl,eivc,dynmat,mp,ddyn,temp)

 3 format(a,i5,9(1x,f19.9))
 4 format(99(1x,2(1x,g9.3),1x))
 5 format(a,i5,99(1x,2(1x,g9.3),1x))
 6 format(a,99(1x,i5))
 7 format(i6,2x,3(1x,f7.3),1x,99(1x,f7.3))
 8 format(99(1x,2(1x,f9.3),1x))
 9 format(a,i5,99(1x,2(1x,f9.3),1x))

 end subroutine get_freq
 !---------------------------------------------------------
 subroutine dir2cart_g(q,k)
! takes q in direct coordinates and outputs k in cartesian
 use geometry
 use lattice
 real(8) q(3),k(3)

 k(:) = q(1)*g1 + q(2)*g2 + q(3)*g3

 end subroutine dir2cart_g
!===========================================================
!=====================================================
  subroutine finitedif_vel(q0,ndn,vgr,evl0,evc0)
! calculates the group velocities in units of c_light from finite difference
! it would give zero near band crossings, thus HF is a better way to do it
  use constants
  implicit none
  integer ndn,i,j
  real(8) q0(3),vgr(3,ndn),q1(3),dq,om0,om1,aux(3,ndn),aux0(3,ndn)
  real(8) evl0(ndn),evlp(ndn),evlm(ndn)
  complex(8) evc0(ndn,ndn),evct(ndn,ndn)

  dq=3d-4  ! this is the ideal dq (at least for si)

   call get_freq(q0,ndn,aux0,evl0,evc0)

!  return

  do i=1,3
     q1=q0; q1(i)=q0(i)+dq
     call get_freq(q1,ndn,aux,evlp,evct)
     q1=q0; q1(i)=q0(i)-dq
     call get_freq(q1,ndn,aux,evlm,evct)
!    vgr(i,:)=(evlp-evlm)/2/dq / 2/sqrt(abs(evl0)) *cnst*1d-8 *2*pi !*c_light
     vgr(i,:)=(evlp-evlm)/2/dq   *1d-8 *2*pi !*c_light
  enddo
!  write(30,6)'q,dq,om=',q0,dq,evl0
!  do i=1,3
!     write(30,5)'i,fd_vg(i),hf_vg(i)=',i,vgr(i,:)*c_light,aux0(i,:)*c_light
!  enddo
!  write(30,*)'*******************************************'


 4 format(a,3(1x,f9.3),a)
 5 format(a,i5,99(1x,f9.3))
 6 format(a,3(1x,f9.3),1x,g12.6,9(1x,g10.4))

  end subroutine finitedif_vel
!============================================================
      subroutine findatom_sc(n3,tau,iatom)
! finds the index of the atom in supercell with identity (n3[modulo supercell],tau)
! arguments:
!     n3(3) (input), linear combination of basis vectors of the primitive
!          lattice that takes us to the unit cell containing the ith atom
!     tau (input), identity of equivalent atom in unit cell at origin
!     iatom (output), index of equivalent atom in supercell. Returns zero if not found
      use geometry
      use lattice
      use ios
      use atoms_force_constants
      implicit none
      integer n3(3),tau,iatom,j,m(3),l,isin
      real(8) a(3),b(3),zero(3)

! first check to see if it is one of the sc atoms
      zero = 0d0
      iatom = 0
      jloop: do j=1,natom_super_cell
         if (atom_sc(j)%cell%tau .eq. tau ) then
             m = atom_sc(j)%cell%n - n3
! find direct coordinates of n3 - n(j) and see if it is integer
             a = matmul(r0g,dfloat(m))
! if the 3 components of a are integer, then iatom=j
!             m = floor(a+0.00001)
!             b = a - m
!             if (b .myeq. zero ) then
!                iatom = j
!                return
!             endif
              if(abs(a(1)-nint(a(1))).lt.0.0001 .and.  &
    &            abs(a(2)-nint(a(2))).lt.0.0001 .and.  &
    &            abs(a(3)-nint(a(3))).lt.0.0001) then
                  iatom=j
                  exit jloop
              endif

         endif
      enddo jloop

      end subroutine findatom_sc
!======================================================================
subroutine ocfunc_sy(temp,omz,l2,nk2,l3,nk3,delta1,delta2,delta3,delta4,delta_tot)
! for given temperature and complex omz (in 1/cm) calculates the
! occupation-weighted joint DOS stored in res
! j2 and j3 are kpoint indices in the full FBZ
 use eigen
 use kpoints   ! for nkc
!use constants ! for cnst
!use params    ! for classical

 implicit none
 integer , intent(in) :: l2,l3,nk2,nk3
 integer j2,j3
 real(8), intent (in) :: temp
 real(8) om3,om2,nb3,nb2,nbe,resr,resi,et,om,delta_g,delta_l
 complex(8), intent (in) :: omz
 complex(8) delta1, delta2, delta3, delta4, delta_tot

!! nk2 = 1+mod(j2-1,nkc) ; l2 = 1+ (j2-nk2)/nkc
!! nk3 = 1+mod(j3-1,nkc) ; l3 = 1+ (j3-nk3)/nkc
! call nkla(j2,nkc,nk2,l2)
! call nkla(j3,nkc,nk3,l3)
 j2=(l2-1)*nkc+nk2
 j3=(l3-1)*nkc+nk3

 om3 = eigenval(l3,nk3) ; nb3=nbe(om3,temp,classical)
 om2 = eigenval(l2,nk2) ; nb2=nbe(om2,temp,classical)
! om=real(omz) ; et=-aimag(omz)

 delta1=(nb2+nb3+1)*(1d0/(om2+om3+omz))
 delta2=(nb2+nb3+1)*(1d0/(om2+om3-omz))
 delta3=(nb2-nb3  )*(1d0/(om3-om2+omz))
 delta4=(nb2-nb3  )*(1d0/(om3-om2-omz))
 delta_tot = delta1+delta2+delta3+delta4


 end subroutine ocfunc_sy
!===========================================================
function ocfunc(temp,omz,l2,nk2,l3,nk3)  result(res)
! for given temperature and complex omz (in 1/cm) calculates the
! occupation-weighted joint DOS stored in res
! j2 and j3 are kpoint indices in the full FBZ
 use eigen
 use kpoints   ! for nkc
! use params    ! for classical
! use constants ! for cnst
 implicit none
 integer , intent(in) :: l2,l3,nk2,nk3
 integer j2,j3
 real(8), intent (in) :: temp
 real(8) om3,om2,nb3,nb2,nbe,resr,resi,et,om,delta_g,delta_l
 complex(8), intent (in) :: omz
 complex(8) res

!! nk2 = 1+mod(j2-1,nkc) ; l2 = 1+ (j2-nk2)/nkc
!! nk3 = 1+mod(j3-1,nkc) ; l3 = 1+ (j3-nk3)/nkc
! call nkla(j2,nkc,nk2,l2)
! call nkla(j3,nkc,nk3,l3)
 j2=(l2-1)*nkc+nk2
 j3=(l3-1)*nkc+nk3

 om3 = eigenval(l3,nk3)  ; nb3=nbe(om3,temp,classical)
 om2 = eigenval(l2,nk2)  ; nb2=nbe(om2,temp,classical)
!--------my revision-------
! om3 = ABS(eigenval(l3,nk3)) ; nb3=nbe(om3,temp,classical)
! om2 = ABS(eigenval(l2,nk2)) ; nb2=nbe(om2,temp,classical)
!--------------------------

! om3 = sqrt(abs(eigenval(l3,nk3))) * cnst ; nb3=nbe(om3,temp,classical)
! om2 = sqrt(abs(eigenval(l2,nk2))) * cnst ; nb2=nbe(om2,temp,classical)
! om=real(omz) ; et=-aimag(omz)

 if(j2.ne.j3) then
 !88888888888888888888888888888888888888888888888888888
 ! need a factor of 2 since v3(q,1,2)^2 * (f(1,2)+f(2,1))
 ! the factor of 2 is needed if np<nq (upper half) is excluded in the case where
 ! both sums are done in the full FBZ
 !88888888888888888888888888888888888888888888888888888
!    resi=-pi*((nb2+nb3+1)*(-delta_l(om3+om2-om,et)+ delta_l(om3+om2+om,et)) +  &
! &            (nb2-nb3  )*(-delta_l(om3-om2-om,et)+ delta_l(om3-om2+om,et)) )
!    resr=real(-(nb2+nb3+1)*(1d0/(om3+om2-omz)+ 1d0/(om3+om2+omz))    &
! &            -(nb2-nb3  )*(1d0/(om3-om2-omz)+ 1d0/(om3-om2+omz)) )
     res =     -(nb2+nb3+1)*(1d0/(om3+om2-omz)+ 1d0/(om3+om2+omz))    &
  &            -(nb2-nb3  )*(1d0/(om3-om2-omz)+ 1d0/(om3-om2+omz))
 else
!    resi=-(2*nb2+1) * pi* (-delta_l(2*om2-om,et)+ delta_l(2*om2+om,et))
!    resr=-(2*nb2+1) * real(1d0/(2*om2-omz)+ 1d0/(2*om2+omz))
     res =-(2*nb2+1) *     (1d0/(2*om2-omz)+ 1d0/(2*om2+omz))
 endif
! res=cmplx(resr,resi)

 end function ocfunc
!-------------------------------------------------------------------------------------------
function ocfunc_36(temp,omz,om2, om3)  result(res)
! for given temperature and complex omz (in 1/cm) calculates the
! occupation-weighted joint DOS stored in res
 use eigen
 use kpoints   ! for nkc
 implicit none
 integer j2,j3
 real(8), intent (in) :: temp, om2, om3
 real(8) nb3,nb2,nbe,resr,resi,et,om,delta_g,delta_l
 complex(8), intent (in) :: omz
 complex(8) res


 nb3=nbe(om3,temp,classical)
 nb2=nbe(om2,temp,classical)


 if(om2.ne.om3) then
     res =     -(nb2+nb3+1)*(1d0/(om3+om2-omz)+ 1d0/(om3+om2+omz))    &
  &            -(nb2-nb3  )*(1d0/(om3-om2-omz)+ 1d0/(om3-om2+omz))
 else
     res =-(2*nb2+1) * (1d0/(2*om2-omz)+ 1d0/(2*om2+omz))
 endif

 end function ocfunc_36
