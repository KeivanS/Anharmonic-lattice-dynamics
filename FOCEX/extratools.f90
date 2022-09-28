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
!======================================================================================
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
      integer, intent(in):: n3(3),tau
      integer j,m(3)
      integer, intent(out):: iatom
      real(8) a(3),zero(3)

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
 subroutine make_r0g
! matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
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
!3 format(a,3(1x,g13.6))

 end subroutine check
!============================================================
 subroutine check_d(r,a1,a2,a3,ier)
! subroutine to check whether r is an integer vector
! output ai are the coefficients of its linear combination on this basis
! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use ios
 implicit none
 type(vector) :: r
 real(8) a1,a2,a3,ep
 integer ier

 ep = tolerance ! 1d-3
!call write_out(ulog,'CHECK: R ',r)
 a1 = r%x
 a2 = r%y
 a3 = r%z
 if(abs(a1-nint(a1)).ge.ep .or.abs(a2-nint(a2)).ge.ep   &
&   .or.abs(a3-nint(a3)).ge.ep ) then
   ier = 1
 else
   ier = 0
 endif

 end subroutine check_d
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
 if(m.eq.0) then ! quantum calculation
   if (z.gt.60) then
     y = 0
   elseif(z.lt.0.0010) then
     y = 1/(z*(1+z/2*(1+z/3)))
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
 type(vector) g1,g2,g3
 real(8) q(3),qx,qy,qz,q2(3)
! integer n
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
    kdg = k(1)*g1%x + k(2)*g1%y
    q2 =  g1%x*g1%x + g1%y*g1%y
    if( kdg .lt. q2/2 .and. kdg .ge.0) inx=.true.
    kdg = k(1)*g2%x + k(2)*g2%y
    q2 =  g2%x*g2%x + g2%y*g2%y
    if ( kdg .lt. q2/2 .and. kdg .ge.0) inside=inx
!   write(6,9)'k,g01,g02=',k,g1,g2,kdg,q2
 endif
!if (inside) write(444,*)"k HEX=",k
!9 format(a,99(1x,f6.3))
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
 function indexg(i,j,k,nil,njl,njh,nkl,nkh) result(n)
! function indexg(i,j,k,nil,nih,njl,njh,nkl,nkh) result(n)
! finds the index n of the kpoint defined with 3 loop indices
! ijk going in general from nil to nih, njl to njh nkl to nkh
! n=0; do i1=nil,nih; do j=njl,njh; do k=nkl,nkh; n=n+1
 implicit none
 integer i,j,k,nil,njl,njh,nkl,nkh,n

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
! function index_reg(i,j,k,n1,n2,n3) result(n)
!! finds the index n of the coarse kpoint defined with 3 loop indices
!! ijk going from 1 to ni , i being the outer loop and k the inner one
! implicit none
! integer i,j,k,n1,n2,n3,n
!
! n = k + (j-1)*n3 + (i-1)*n2*n3
!
! end function index_reg
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

 ns = nint((2*nshell+1)**(0.3333))
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
 enddo

3 format(i6,9(2x,g11.4))

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

!7 format(a,i5,9(1x,g12.5))
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

 qdr = (q.dot.rr) + epsl
! include from -0.00005 to 0.99995 included (include 0.0 and exclude 1.0)
 if (qdr.lt.0 .or. qdr.ge.1) insid=0
 aux = qdr-floor(qdr)
 i = nint(1+ n*aux)
! write(*,7)'i,qdotr,aux=',i,qdr,aux
 if (i.eq.n+1) i=1

!7 format(a,i5,9(1x,g12.5))
 end subroutine comp1
!============================================================
 subroutine comp_c(q,rr,n,i,insid)
! inside would be =1 if -0.5<q<0.5
 use geometry
 implicit none
 real(8) q(3),aux,qdr
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

!7 format(a,i5,9(1x,g12.5))
 end subroutine comp_c
!============================================================
 subroutine get_components_g(q,n,i,j,k,inside)
! subroutine get_components_g(q,n,i,j,k,gg1,gg2,gg3,inside)
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
! type(vector),intent(in):: gg1,gg2,gg3
! type(vector)rx1,rx2,rx3
 integer, intent(in):: n(3)
 integer, intent(out):: i,j,k,inside


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
 subroutine get_components_g_centered(q,n,i,j,k,inside)
! subroutine get_components_g_centered(q,n,i,j,k,gg1,gg2,gg3,inside)
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
! type(vector),intent(in):: gg1,gg2,gg3
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

!3 format(3(i6),9(2x,f8.4))

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
  integer n
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
  real(8) xmin,xmax,dx,e(mesh)

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
! 4 format(i8,9(2x,g13.6))
 5 format(a,i8,9(2x,g13.6))

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
