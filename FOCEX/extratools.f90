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
!! finds the index of the atom in supercell with identity (n3[modulo supercell],tau)
!! arguments:
!!     n3(3) (input), linear combination of basis vectors of the primitive
!!          lattice that takes us to the unit cell containing the ith atom
!!     tau (input), identity of equivalent atom in unit cell at origin
!!     iatom (output), index of equivalent atom in supercell. Returns zero if not found
      use geometry
      use lattice
      use ios
      use atoms_force_constants
      implicit none
      integer, intent(in) :: n3(3),tau
      integer, intent(out) :: iatom
      integer j,m(3),l,isin
      real(8) a(3),b(3),zero(3)

! first check to see if it is one of the sc atoms
      zero = 0d0
      iatom = 0
      jloop: do j=1,natom_super_cell
         if (atom_sc(j)%cell%tau .eq. tau ) then
             m = atom_sc(j)%cell%n - n3
! find direct coordinates of n3 - n(j) in the supercell and see if it is integer
! to do this form (sum_i m(i)*r0(i)).gsc(j) this should be integer
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
! coordinates are v(1)*r01+v(2)*r02+v(3)*r03
 use lattice
 use constants , only : pi
 implicit none

 r0g(1,1) = r01 .dot. gs1
 r0g(2,1) = r01 .dot. gs2
 r0g(3,1) = r01 .dot. gs3
 r0g(1,2) = r02 .dot. gs1
 r0g(2,2) = r02 .dot. gs2
 r0g(3,2) = r02 .dot. gs3
 r0g(1,3) = r03 .dot. gs1
 r0g(2,3) = r03 .dot. gs2
 r0g(3,3) = r03 .dot. gs3
 r0g=r0g/(2*pi)

 end subroutine make_r0g
!======================================================================
 subroutine make_rg(x1,x2,x3,q1,q2,q3,n)
! matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
! matmul(r0g,v) gives the 3 reduced coordinates of an atom in the supercell if its
! reduced coordinates in the primitive cell are given by v
 use geometry
 implicit none
 real(8) x1(3),x2(3),x3(3),q1(3),q2(3),q3(3),n(3,3)

 n(1,1) = x1 .dot. q1
 n(1,2) = x1 .dot. q2
 n(1,3) = x1 .dot. q3
 n(2,1) = x2 .dot. q1
 n(2,2) = x2 .dot. q2
 n(2,3) = x2 .dot. q3
 n(3,1) = x3 .dot. q1
 n(3,2) = x3 .dot. q2
 n(3,3) = x3 .dot. q3

 end subroutine make_rg
!============================================================
 subroutine check(r,a1,a2,a3,ier,g1,g2,g3)
!! subroutine to check whether r is an integer multiple of (r01,r02,r03)
!! output ai are the coefficients of its linear combination on this basis
!! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use constants
 use ios
 implicit none

 type(vector) :: r,g1,g2,g3
 real(8) a1,a2,a3,ep
 integer ier

 ep = tolerance ! 0.001 ! displacements larger than 0.001 A are not tolerated
!call write_out(ulog,'CHECK: R ',r)
 a1 = (r .dot. g1)/(2*pi)  ! g's have a 2 pi
 a2 = (r .dot. g2)/(2*pi)
 a3 = (r .dot. g3)/(2*pi)
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
3 format(a,3(1x,g13.6))

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
! r(mcor(1)) < r(mcor(2)) < r(mcor(3)) < ... < r(mcor(n))  
  implicit none
  integer maxat,n,i,j,temp
  real(8) r(maxat)
  integer mcor(maxat)

  do i=1,maxat ! n
     mcor(i)=i
  enddo
  do i=1,n-1  ! was 1,n-1
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
!===============================================================
 subroutine make_reciprocal_lattice_a(r1,r2,r3,g1,g2,g3)
! be careful it does not have the factor of 2*pi in it!
 use geometry
 use constants
 use ios
 implicit none
 real(8), intent(in) :: r1(3),r2(3),r3(3)
 real(8), intent(out):: g1(3),g2(3),g3(3)
 real(8) om

 om = abs(r1 .dot. (r2 .cross. r3))
 if (om.lt.1d-8) then
    write(*,*)'MAKE_RECIPROCAL_LATTICE_A: volume is zero; check your translation vectors'
    stop
 endif
 g1=(r2 .cross. r3)/om
 g2=(r3 .cross. r1)/om
 g3=(r1 .cross. r2)/om

 end subroutine make_reciprocal_lattice_a
!===============================================================
 subroutine make_reciprocal_lattice_v(r1,r2,r3,g1,g2,g3)
! be careful it does not have the factor of 2*pi in it!
 use geometry
 use constants
 use ios
 implicit none
 type(vector) :: r1,r2,r3,g1,g2,g3
 real(8) om

 om = abs(r1 .dot. (r2 .cross. r3))
 if (om.lt.1d-8) then
    write(*,*)'MAKE_RECIPROCAL_LATTICE_V: volume is zero; check your translation vectors'
    stop
 endif
! write(ulog,*)'RECIPROCAL_LATTICE: '
 g1=(r2 .cross. r3)/om
 g2=(r3 .cross. r1)/om
 g3=(r1 .cross. r2)/om
! call write_out(ulog,'om ',om)
! call write_out(ulog,'g1 ',g1)
! call write_out(ulog,'g2 ',g2)
! call write_out(ulog,'g3 ',g3)

 end subroutine make_reciprocal_lattice_v
!===============================================================
 subroutine make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
! be careful it does not have the factor of 2*pi in it!
 use geometry
 use constants
 use ios
 implicit none
 type(vector) :: r1,r2,r3,g1,g2,g3
 real(8) om

 om = abs(r1 .dot. (r2 .cross. r3))
 if (om.lt.1d-8) then
    write(*,*)'MAKE_RECIPROCAL_LATTICE_2PI: volume is zero; check your translation vectors'
    stop
 endif
! write(ulog,*)'RECIPROCAL_LATTICE: '
 g1=2*pi*(r2 .cross. r3)/om
 g2=2*pi*(r3 .cross. r1)/om
 g3=2*pi*(r1 .cross. r2)/om
! call write_out(ulog,'om ',om)
! call write_out(ulog,'g1 ',g1)
! call write_out(ulog,'g2 ',g2)
! call write_out(ulog,'g3 ',g3)

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
    kdg = k(1)*g01%x + k(2)*g01%y
    q2 =  g01%x*g01%x + g01%y*g01%y
    if( kdg .lt. q2/2 .and. kdg .ge.0) inx=.true.
    kdg = k(1)*g02%x + k(2)*g02%y
    q2 =  g02%x*g02%x + g02%y*g02%y
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
 use ios
 use params
 implicit none
 integer, intent(in) ::  nshell
 type(vector), intent(in) :: g1,g2,g3
 real(8), intent(out) :: gg(3,nshell)
 integer i,j,k,n5,ik,ns
 real(8), allocatable :: g5(:,:),gs(:)
 integer, allocatable :: xmap(:)

 ns = (2*nshell+1)**(0.3333) ! this is larger than the radius corresponding to nshell points
 n5 = (2*ns+1)**3 ! this therefore an upper bound to the needed number
 allocate ( g5(3,n5) )

 ik = 0
 do i=-ns,ns
 do j=-ns,ns
 do k=-ns,ns
    ik = ik+1
    g5(:,ik) = i*g1 + j*g2 + k*g3
!   gs(ik) = length(g5(:,ik))
 enddo
 enddo
 enddo

 call sort_grid(n5,g5)

 gg(:,1:nshell)=g5(:,1:nshell)
 if(verbose) then
    write(ulog,*)' Sorted G vectors and their length' 
    do i=1,nshell
       write(ulog,3)i,g5(:,i),length(g5(:,i))
    enddo
 endif

3 format(i6,9(2x,g11.4))

 deallocate (g5)
 end subroutine make_sorted_gs
!===========================================================
 subroutine sort_grid(nshell,gg)
! from the basis (g1,g2,g3) generate the nshell shortest linear combinations including 0
 use geometry
 implicit none
 integer nshell,i,j,k,n5,ik,ns
 real(8) gg(3,nshell)
 real(8) g5(3,nshell),gs(nshell)
 integer xmap(nshell)

 do i=1,nshell
    gs(i) = length(gg(:,i))
 enddo

! now sort g2 in ascending order
 call sort(nshell,gs,xmap,nshell)

 g5=gg
 do i=1,nshell
    gg(:,i) = g5(:,xmap(i))
 enddo

 end subroutine sort_grid
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

 call make_reciprocal_lattice_v(g1,g2,g3,rr1,rr2,rr3)  ! no factor of 2pi
 inside = 1

 qx = (q.dot.rr1) + epsl
 if (qx.lt.0 .or. qx.ge.1) inside=0
! write(*,7)'i,qdotr,aux=',i,qdr,aux

 qy = (q.dot.rr2) + epsl
 if (qy.lt.0 .or. qy.ge.1) inside=0

 qz = (q.dot.rr3) + epsl
 if (qz.lt.0 .or. qz.ge.1) inside=0

7 format(a,i5,9(1x,g12.5))
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

7 format(a,i5,9(1x,g12.5))
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

7 format(a,i5,9(1x,g12.5))
 end subroutine comp_c
!============================================================
! subroutine get_components_g(q,n,i,j,k,gg1,gg2,gg3,inside)
 subroutine get_components_g(q,n,i,j,k,inside)  ! g's  are primitive g0i
! for a given q-vector, it finds its integer components assuming it was
! created as: q=(i-1)/n1*g1 + (j-1)/n2*g2 + (k-1)/n3*g3; i=1,n1 j=1,n2 k=1,n3
! as the ones created in make_kp_reg with zero shift
! if the integer variable inside=1 then q is inside the primcell
! if q is outside the prim_cell, then inside=0 but produces the i,j,k of its
! image inside
 use lattice
 use geometry
 use constants, only : pi
 implicit none
 real(8), intent(in):: q(3)
! type(vector),intent(in):: gg1,gg2,gg3
! type(vector)rx1,rx2,rx3
 integer, intent(in):: n(3)
 integer, intent(out):: i,j,k,inside
 real(8) i2pi

 i2pi=1d0/(2*pi)
! call make_reciprocal_lattice(g1,g2,g3,rr1,rr2,rr3)
! i = nint(1+ n(1)* q.dot.r1)
! j = nint(1+ n(2)* q.dot.r2)
! k = nint(1+ n(3)* q.dot.r3)
! if q.dot.r is larger than 1, we want to have its fractional part
 inside = 1

 call comp1(q,i2pi*r01,n(1),i,inside)
 call comp1(q,i2pi*r02,n(2),j,inside)
 call comp1(q,i2pi*r03,n(3),k,inside)

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
 use constants, only : pi
 implicit none
 real(8), intent(in):: q(3)
 type(vector),intent(in):: gg1,gg2,gg3
 integer, intent(in):: n(3)
 integer, intent(out):: i,j,k,inside
! real(8) w(3)
 real(8) i2pi

 i2pi=1d0/(2*pi)

! w = q+0.5d0*(gg1+gg2+gg3)
! write(*,8)'q=',q
! write(*,8)'w=',w

 inside = 1

 call comp_c(q,i2pi*r01,n(1),i,inside)
 call comp_c(q,i2pi*r02,n(2),j,inside)
 call comp_c(q,i2pi*r03,n(3),k,inside)

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
 subroutine send_to_primcell(kp,q)
! folds the kpoint kp into the primitive cell (between 0 and G) stores the result in q
 use lattice
!! use io2
 implicit none
 real(8), intent(in) :: kp(3)
 real(8), intent(out) :: q(3)
 real(8) a1,a2,a3

 a1 = (kp .dot. r01)/2/pi
 a2 = (kp .dot. r02)/2/pi
 a3 = (kp .dot. r03)/2/pi
 a1 = a1-floor(a1)
 a2 = a2-floor(a2)
 a3 = a3-floor(a3)
 q=a1*g01+a2*g02+a3*g03

 end subroutine send_to_primcell
!============================================================
 subroutine find_in_grid(v,n,grid,y1,y2,y3,nj,tolerance)
! checks whether the vector v belongs to the array "grid" modulo its period within the tolerance;
! nj is its index in the grid; if nj=0, the vector did not belong to the grid
 use geometry
 implicit none
 integer nj,j,n
 real(8) v(3),grid(3,n),tolerance,y1(3),y2(3),y3(3),i1,i2,i3

 nj=0
 do j=1,n
    i1=(v-grid(:,j)).dot. y1
    if(abs(i1-nint(i1)).gt.tolerance) cycle
    i2=(v-grid(:,j)).dot. y2
    if(abs(i2-nint(i2)).gt.tolerance) cycle
    i3=(v-grid(:,j)).dot. y3
    if(abs(i3-nint(i3)).gt.tolerance) cycle
    nj=j
    exit
 enddo

 end subroutine find_in_grid
!============================================================
 subroutine find_in_array(v,n,grid,nj,tolerance)
! checks whether the vector v belongs to the array "grid" within the tolerance;
! nj is its index in the grid; if nj=0, the vector did not belong to the grid
 use geometry
 implicit none
 integer nj,j,n
 real(8) v(3),grid(3,n),tolerance

 nj=0
 do j=1,n
    if(length(v-grid(:,j)) .lt. tolerance) then
       nj=j
       exit
    endif
 enddo

 end subroutine find_in_array
!============================================================
 subroutine inverse_real(a,b,n)
!! inverts a real matrix a(nxn)
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
 4 format(i8,9(2x,g13.6))
 5 format(a,i8,9(2x,g13.6))

  close(unit)

  end subroutine histogram
!============================================================
 function trace(mat)
 implicit none
 integer i
 real(8) trace,mat(:,:)

 trace=0
 do i=1,size(mat(1,:))
    trace=trace+mat(i,i)
 enddo

 end function trace
!===========================================================
 recursive function determinant(mat) result(det)
 implicit none
 integer n,i,j,k
! real(8) det,mat(n,n),cofact(n-1,n-1)
 real(8), dimension(:,:), intent(in)  :: mat
 real(8), allocatable :: cofact(:,:)
 real(8) det

 n=size(mat(1,:))
 allocate(cofact(n-1,n-1))

 if (n.eq.1) then
    det=mat(1,1)
 else
    det=0
    do i=1,n
       do j=1,n
       do k=2,n
          if (j.lt.i) then
             cofact(k-1,j)=mat(k,j)
          elseif(j.gt.i) then
             cofact(k-1,j-1)=mat(k,j)
          endif
       enddo
       enddo
       det=det+mat(1,i)*determinant(cofact)*(-1)**(i-1)
    enddo
 endif
 end function determinant
!===========================================================
     DOUBLE PRECISION FUNCTION erfc(x) result(S15ADF)
!     MARK 5A REVISED - NAG COPYRIGHT 1976
!     MARK 5C REVISED
!     MARK 11.5(F77) REVISED. (SEPT 1985.)
!     COMPLEMENT OF ERROR FUNCTION ERFC(X)
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
!     .. Local Scalars ..
      DOUBLE PRECISION                 T, XHI, XLO, Y
!     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP
!     .. Data statements ..
!     PRECISION DEPENDENT CONSTANTS
! 08   DATA XLO/-4.5D0/
! 12   DATA XLO/-5.25D0/
! 14   DATA XLO/-5.75D0/
      DATA XLO/-6.25D0/
! 18   DATA XLO/-6.5D0/
!
!     RANGE DEPENDENT CONSTANTS
      DATA XHI/ 2.66D+1 /
!     XHI = LARGEST X SUCH THAT EXP(-X*X) .GT. MINREAL (ROUNDED DOWN)
! R1   DATA XHI/13.0D0/
! R2   DATA XHI/9.5D0/
! R3   DATA XHI/13.0D0/
! R4   DATA XHI/25.0D0/
! R5   DATA XHI/26.0D0/
!     .. Executable Statements ..
!
!     TEST EXTREME EXITS
      IF (X.GE.XHI) then
           s15adf=0
           return
      endif
      IF (X.LE.XLO) then
           s15adf=2
           return
      endif
!
!     EXPANSION ARGUMENT
      T = 1.0D0 - 7.5D0/(ABS(X)+3.75D0)
!
!      * EXPANSION (0021) *
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 08E
! 08   Y = (((((((((((+3.1475326D-5)*T-1.3874589D-4)*T-6.4127909D-6)
! 08  *    *T+1.7866301D-3)*T-8.2316935D-3)*T+2.4151896D-2)
! 08  *    *T-5.4799165D-2)*T+1.0260225D-1)*T-1.6357229D-1)
! 08  *    *T+2.2600824D-1)*T-2.7342192D-1)*T + 1.4558972D-1
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 12E
! 12   Y = ((((((((((((((((-4.21661579602D-8*T-8.63384346353D-8)
! 12  *    *T+6.06038693567D-7)*T+5.90655413508D-7)
! 12  *    *T-6.12872971594D-6)*T+3.73223486059D-6)
! 12  *    *T+4.78645837248D-5)*T-1.52546487034D-4)
! 12  *    *T-2.55222360474D-5)*T+1.80299061562D-3)
! 12  *    *T-8.22062412199D-3)*T+2.41432185990D-2)
! 12  *    *T-5.48023263289D-2)*T+1.02604312548D-1)
! 12  *    *T-1.63571895545D-1)*T+2.26008066898D-1)
! 12  *    *T-2.73421931495D-1)*T + 1.45589721275D-1
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 14E
! 14   Y = (((((((((((((((-2.2356173494379D-9
! 14  *    *T+4.5302502889845D-9)*T+2.5918103316137D-8)
! 14  *    *T-6.3684846832869D-8)*T-1.7642194353331D-7)
! 14  *    *T+6.4907607131235D-7)*T+7.4296952017617D-7)
! 14  *    *T-6.1758018478516D-6)*T+3.5866167916231D-6)
! 14  *    *T+4.7895180610590D-5)*T-1.5246364229106D-4)
! 14  *    *T-2.5534256252531D-5)*T+1.8029626230333D-3)
! 14  *    *T-8.2206213481002D-3)*T+2.4143223946968D-2)
! 14  *    *T-5.4802326675661D-2)*T+1.0260431203382D-1
! 14   Y = (((Y*T-1.6357189552481D-1)*T+2.2600806691658D-1)
! 14  *    *T-2.7342193149541D-1)*T + 1.4558972127504D-1
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 16E
      Y = (((((((((((((((+3.328130055126039D-10        &
     &    *T-5.718639670776992D-10)*T-4.066088879757269D-9)       &
     &    *T+7.532536116142436D-9)*T+3.026547320064576D-8)       &
     &    *T-7.043998994397452D-8)*T-1.822565715362025D-7)       &
     &    *T+6.575825478226343D-7)*T+7.478317101785790D-7)       &
     &    *T-6.182369348098529D-6)*T+3.584014089915968D-6)       &
     &    *T+4.789838226695987D-5)*T-1.524627476123466D-4)       &
     &    *T-2.553523453642242D-5)*T+1.802962431316418D-3)       &
     &    *T-8.220621168415435D-3)*T+2.414322397093253D-2
      Y = (((((Y*T-5.480232669380236D-2)*T+1.026043120322792D-1)       &
     &    *T-1.635718955239687D-1)*T+2.260080669166197D-1)       &
     &    *T-2.734219314954260D-1)*T + 1.455897212750385D-1
!
!     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 18E
! 18   Y = (((((((((((((((-1.58023488119651697D-11
! 18  *    *T-4.94972069009392927D-11)*T+1.86424953544623784D-10)
! 18  *    *T+6.29796246918239617D-10)*T-1.34751340973493898D-9)
! 18  *    *T-4.84566988844706300D-9)*T+9.22474802259858004D-9)
! 18  *    *T+3.14410318645430670D-8)*T-7.26754673242913196D-8)
! 18  *    *T-1.83380699508554268D-7)*T+6.59488268069175234D-7)
! 18  *    *T+7.48541685740064308D-7)*T-6.18344429012694168D-6)
! 18  *    *T+3.58371497984145357D-6)*T+4.78987832434182054D-5)
! 18  *    *T-1.52462664665855354D-4)*T-2.55353311432760448D-5
! 18   Y = ((((((((Y*T+1.80296241673597993D-3)
! 18  *    *T-8.22062115413991215D-3)
! 18  *    *T+2.41432239724445769D-2)*T-5.48023266949776152D-2)
! 18  *    *T+1.02604312032198239D-1)*T-1.63571895523923969D-1)
! 18  *    *T+2.26008066916621431D-1)*T-2.73421931495426482D-1)*T +
! 18  *     1.45589721275038539D-1
!
      S15ADF = EXP(-X*X)*Y
      IF (X.LT.0.0D0) then
            S15ADF = 2.0D0 - S15ADF
      endif

      END function erfc !s15adf
!==================================================
      SUBROUTINE choldc(a,n,p)
! Choleski decomposition of a=LLT; p contains the diagonals of L
! upper triangular part of a is used; Lower part of L is stored in lower part of a
      implicit none
      INTEGER n
      REAL(8) a(n,n),p(n)
      INTEGER i,j,k
      REAL(8) sum
      do i=1,n
        do j=i,n
          sum=a(i,j)
          do  k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
          enddo
          if(i.eq.j)then
            if(sum.le.0.) then
                write(*,*) 'choldc failed; at i=',i,sum
                stop !read(*,*)
            endif
            p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
      enddo
      enddo

      END subroutine choldc
!  (C) Copr. 1986-92 Numerical Recipes Software !+!).
!========================================
 subroutine make_grid(x01,x02,x03,x1,x2,x3,n,grid)
!! generates a grid of vectors "grid" linear combinations of x0i but
!! contained in the Wigner-Seitz cell of x1,x2,x3, including the boundaries
 use geometry
 implicit none
 integer, intent(in) :: n
! real(8), intent(in ):: x01(3),x02(3),x03(3),x1(3),x2(3),x3(3)
 type(vector), intent(in ):: x01,x02,x03,x1,x2,x3
 real(8), intent(out):: grid(3,n)
 real(8) v(3),a1,a2,a3
 type(vector) y1,y2,y3
 integer i1,i2,i3,m,cnt


 call make_reciprocal_lattice_v(x1,x2,x3,y1,y2,y3)
 m=50; cnt=1
 do i1=-m,m
 do i2=-m,m
 do i3=-m,m

! generate vectors in a grid
    v= v2a(i1*x01 + i2*x02 + i3*x03)

! find reduced coordinates of a vector on the grid v=sum_i ai*xi
    a1 = v .dot. y1
    a2 = v .dot. y2
    a3 = v .dot. y3

! see if -0.5<ai=<0.5 to see if it's within the Wigner-Seitz cell of the supercell
    if (a1.gt.0.50001d0 .or. a1.le.-0.49999d0) cycle
    if (a2.gt.0.50001d0 .or. a2.le.-0.49999d0) cycle
    if (a3.gt.0.50001d0 .or. a3.le.-0.49999d0) cycle
    grid(:,cnt)=v
    cnt=cnt+1
    write(98,3)v, a1,a2,a3
    if(cnt.gt.n+1) then
       write(*,*)'too many vectors generated',cnt,' compared to ',n
       stop
    endif

 enddo
 enddo
 enddo

3 format(9(1x,f10.5))
 end subroutine make_grid
!============================================================
 subroutine get_upper_bounds(x1,x2,x3,rcut,maxp,m)
!! finds maxp, the number of grid points in the sphere of radius rcut and m, the needed mesh
!! outputs are m and maxp: upper bounds in 1D and in 3D (used for 3 loops and allocation)
 use geometry
 implicit none
 type(vector),intent(in) :: x1,x2,x3
 integer,intent(out) :: m,maxp ! upper boundary in 1D and in 3D (used for 3 loops and allocation)
 real(8),intent(in) :: rcut
 integer m1,m2
 real(8) om0,ls

 call calculate_volume(x1,x2,x3,om0)
 ls=om0**0.333333 ! typical unitcell size
 maxp=nint(12.56/3d0*(1.3*(rcut+ls))**3/om0)+1  ! upperbound to the number of unitcells in the sphere

! find shortest translation vector length
 ls=min(length(x1),length(x2))
 ls=min(ls,length(x3))
 write(*,3)'GET_UPPER_BOUNDS: shortest translation length=',ls

 m1=nint(rcut/ls)+2
 m2=(nint((3*maxp/12.56)**0.33333)+1)
 m=max(m1,m2)
! maxp=(2*m+1)**3

 write(*,4)'GET_UPPER_BOUNDS:m1,m,max,rcut=',m1,m,maxp,rcut

3 format(a,9(1x,g14.7))
4 format(a,3i8,9(1x,g14.7))

 end subroutine get_upper_bounds
!============================================================
 subroutine apply_metric(x1,x2,x3,eps,t1,t2,t3)
!! generates vectors ti such that ti=sqrt(eps)xi
 use geometry
 use ios
 implicit none
 type(vector), intent(in ):: x1,x2,x3
! type(vector), intent(out):: t1,t2,t3
 real(8), intent(out) :: t1(3),t2(3),t3(3)
 real(8) eps(3,3),eps2(3,3),sqeps(3,3),a1(3)
 integer i,j

 eps2=eps  ! to keep eps
 call write_out(ulog,'Before Choleski',eps2)
 call choldc(eps2,3,a1)
 do i=1,3
    sqeps(i,i)=a1(i)
    do j=i+1,3
       sqeps(j,i)= eps2(j,i)
       sqeps(i,j)=sqeps(j,i)
    enddo
 enddo

 call write_out(ulog,'After Choleski; upper diag',sqeps)

  eps2=matmul(sqeps,transpose(sqeps))
  call write_out(ulog,'sq*transp(sq)',eps2)



  t1=matmul(transpose(sqeps),v2a(x1))
  t2=matmul(transpose(sqeps),v2a(x2))
  t3=matmul(transpose(sqeps),v2a(x3))

 end subroutine apply_metric
!============================================================
 subroutine make_grid_weights_WS(x01,x02,x03,x1,x2,x3,ngrd,space,sx)
!! generates a grid of vectors "rgrid" or "ggrid" (depending on space) linear combinations
!! of x0i but contained in the Wigner-Seitz cell of x1,x2,x3, with corresponding weights
!! works for fourier interpolation in the reciprocal space because get_stars works for q-points
!! space is r (for real space) or g (for reciprocal space)
!! the 26 shortest translation vectors sx for defining the WS cell of xi is also output
 use geometry
 use ios, only : ulog
 use lattice
 use params
 use fourier
 implicit none
 type(vector), intent(in):: x01,x02,x03,x1,x2,x3
 real(8), intent(out):: sx(3,26)
 integer, intent(out):: ngrd
 character(len=1), intent(in) :: space
 integer i1,i2,i3,m,cnt,j,in,n,k,nb,jmax,narms,j2,nboundary,insid,nbtot
 type(vector) b01,b02,b03,b1,b2,b3
 real(8) v(3),vred(3),w(3),wred(3),s0x(3,26),a1,a2,a3,t1,t2,t3,omx,omx0,prev,minleng,tol
 integer, allocatable :: save_boundary(:)
 real(8), allocatable :: aux(:,:),grd(:,:),weig(:)
! type(vector) y1,y2,y3

! find the shortest basis bi among xi and the corresponding 26 shortest vectors
 call get_26shortest_shell(x01,x02,x03,s0x,b01,b02,b03)! s0x for WS of primcell
 call get_26shortest_shell(x1 ,x2 ,x3 ,sx ,b1 ,b2 ,b3) ! sx for WS of supercell

 tol = tolerance*length(b1)

! y's are used to find the reduced coordinates on the basis xi
! call make_reciprocal_lattice_v(x1,x2,x3,y1,y2,y3)

! find the number of grid points inside the supercell
 call calculate_volume(x1,x2,x3,omx)
 call calculate_volume(x01,x02,x03,omx0)
 n=nint(omx/omx0)
 write(ulog,4)' grid size,om0,om=',n,omx0,omx

! this is an upperbound to the points inside and on the WS boundaries
 n=n+nint(8*n**0.67+12*n**0.3333+24)

! aux contains the grid vectors inside or on the boundary
 allocate(aux(3,n))
 write(ulog,*)' allocated grid size, including boundary points=',n

! create the grid within the WS supercell xi
 m=nint(n**0.3333) !1+floor(length(sx(:,26))/length(sx(:,1))/2) ! ratio of the longest to smallest lengths
 cnt=0
 write(ulog,*)' upper index to loop over to create grid =',m
 do i1=-m,m
 do i2=-m,m
 i3loop: do i3=-m,m

! generate a grid of vectors from primitive translations x0i
    v= v2a(dble(i1)*b01 + dble(i2)*b02 + dble(i3)*b03)

! see if v is in the WS cell; throw away v outside of WS cell
    call check_inside_ws(v,sx,insid)

    if (insid.eq.0) cycle i3loop

! take only those points inside or on the boundary of WS (corresponding to in=1)
    cnt=cnt+1
    write(*,4)'new vector ',cnt,v
    aux(:,cnt)=v

 enddo i3loop
 enddo
 enddo

 ngrd=cnt
 write(ulog,4)' actual grid size=',ngrd
! write(ulog,4)' Number of vectors defining the boundary=',nboundary

 allocate(grd(3,ngrd),save_boundary(ngrd),weig(ngrd))
 save_boundary=0
 grd(:,1:ngrd)=aux(:,1:ngrd)
 deallocate(aux)

! find the weights of each point on the grid (and its boundary)
 weig=1
 nbtot=0
 write(ulog,*)' grid vectors, number on boundary, weight =========='
 do cnt=1,ngrd

! identify vectors on the boundary
    call is_on_boundary(grd(:,cnt), sx,nboundary,tol)
    if(nboundary.ne.0) then
       save_boundary(cnt)=nboundary
       nbtot=nbtot+1
!       weig(cnt)=1d0/nboundary ! vectors on many boundaries have a smaller weight
    endif
    write(ulog,6)'i,grid(i), on how many boundaries ',cnt,nboundary,grd(:,cnt)
 enddo
 write(ulog,*)'Of ',ngrd,' grid points ',nbtot,' are on the boundary'
 write(ulog,6)'CHECK: sum of saveboundary=',nbtot,sum(save_boundary)

 call find_ws_weights(ngrd,grd,save_boundary,sx,weig)

! normalize weig (should normalize to 1)
 weig=weig*omx0/omx

 write(ulog,*)'sum of the weights=',sum(weig(1:ngrd))
 if(abs(sum(weig)-1).gt.tolerance) then
    write(*,*)'WARNING: weights not normalized to 1 ',sum(weig(1:ngrd))
    write(ulog,*)'WARNING: weights not normalized to 1 ',sum(weig(1:ngrd))
 endif


 if(space.eq.'r' .or. space.eq.'R') then
    if (allocated(rgrid)) deallocate(rgrid)
    if (allocated(rws_weights)) deallocate(rws_weights) 
    allocate(rgrid(3,ngrd), rws_weights(ngrd))
    nrgrid=ngrd
    rgrid = grd
    rws_weights = weig 
    open(98,file='rgrid_raw.dat')
    open(99,file='rgridWS.xyz')
 elseif(space.eq.'g' .or. space.eq.'G') then
    if (allocated(ggrid)) deallocate(ggrid)
    if (allocated(gws_weights)) deallocate(gws_weights) 
    allocate(ggrid(3,ngrd), gws_weights(ngrd))
    nggrid=ngrd
    ggrid = grd
    gws_weights = weig / (volume_r/volume_r0) ! introduce 1/N since used for Fourier transforms
    open(98,file='ggrid_raw.dat')
    open(99,file='ggridWS.xyz')
 endif

 write(99,*)ngrd
 write(99,*)" "
 do cnt=1,ngrd
    write(98,6)" ",cnt,save_boundary(cnt),grd(:,cnt),weig(cnt)
    write(99,7)"Si ",grd(:,cnt)
 enddo
 close(98)
 close(99)
 deallocate(weig,grd,save_boundary)

 write(ulog,*) ' EXITING make_grid_weights_WS'
 
!weights=0d0; weights(1:ngrd)=1 ; grid(:,1:ngrd)=grd(:,ngrd)
!do cnt=1,ngrd
!   v=grd(:,cnt)
!   if (save_boundary(cnt).eq.1) then ! take vectors on the WS boundary

! find reduced coordinates of vector v on the grid v=sum_i ai*xi
!         vred(1) = v .dot. y1
!         vred(2) = v .dot. y2
!         vred(3) = v .dot. y3
!         call get_allstars_k(vred,primitivelattice,narms,kvecstar,kvecop)
!         write(97,4)'***** narms for vred,v=',narms,vred,v
!         do i1=1,narms
!            write(97,3)i1,kvecstar(:,i1)
!         enddo

! count how many of the stars are a grid vector on the boundary differing by a xi
!         i3=0  ! counts the boundary stars which are also on the grid
!         armloop: do i1=1,narms

!            w= v2a(kvecstar(1,i1)*x1+kvecstar(2,i1)*x2+kvecstar(3,i1)*x3)
!            call find_in_array(w,ngrd,grd,i2,1d-5)
!            if(i2.eq.0) cycle ! the star #i1 must also be grid vector
!              if(save_boundary(i2).eq.1) then  ! it is also on the boundary

! check the difference and count if the difference is a xi vector (work with reduced)
!                  wred=kvecstar(:,i1)-vred  ! this is reduced coords

!         do k=1,3
!           n=nint(wred(k))
!           if(ncmp(wred(k)-n).ne.0)cycle armloop ! if wred non-integer skip to next arm
!           if(k.eq.3) then
!              i3=i3+1
!              cycle armloop  ! goto next sym_op iff v3=integer : i.e. G-vector
!           endif
!         enddo

!              endif

!         enddo armloop

! the weight should be the inverse of the number of arms in the grid which differ by a G
!         weights(cnt)=(1d0)/i3
!         write(ulog,6)'cnt,# of stars, v, boundary, weight:',cnt,i3,v,vred,weights(cnt)
!   endif

!   write(99,3)cnt,weights(cnt),v, vred
!enddo

2 format(9(1x,f10.5))
3 format(i5,1x,f8.5,3x,3(1x,g10.3),3x,3(1x,f8.4))
4 format(a,i5,9(1x,f13.5))
5 format(3(1x,f10.5),i3,1x,g11.4)
6 format(a,2i5,99(1x,f10.5))
7 format(a,9(1x,g12.5))

 end subroutine make_grid_weights_WS
!============================================================
 subroutine find_WS_largest_SC(imax,volmax)
!! scans over all available supercell POSCARs and picks the one with largest volume
!! also finds the number of force constraints for each supercell
 use params, only : fdfiles
 use ios   , only : ulog
 use lattice, only : volume_r
 implicit none
 integer, intent(out) :: imax  ! index of largest supercell
 real(8), intent(out) :: volmax
 integer i
 character xt*1,poscar*7

     volmax=0             !  pick the SC with largest volume
     do i=1,fdfiles

        write(xt,'(i1)')i
        poscar='POSCAR'//xt
! read the atomic coordinates of atoms in supercell (at equilibrium) from POSCAR
        call read_crystal(poscar)
        write(ulog,*) 'supercell ',i,' read, volume=',volume_r

! find the largest volume
        if(volume_r .gt. volmax) then
           volmax=volume_r
           imax=i
        endif

     enddo
     write(ulog,*) 'largest volume is found for supercell ',i,' volmax=',volmax

! work with this POSCAR
     write(xt,'(i1)')imax
     poscar='POSCAR'//xt
     call read_crystal(poscar)

 end subroutine find_WS_largest_SC
!===========================================================
 subroutine count_configs(outcar,ncfg)
 use ios , only : ulog,utraj
 use atoms_force_constants, only : natom_super_cell
 implicit none
 integer, intent(out) :: ncfg
 character, intent(in):: outcar*(*)
 integer t,j,i
 character line*99
 logical found,exst

 write(*,*)' opening OUTCAR file'

 inquire(file=outcar,exist=exst)
 if(exst) then
      open(utraj,file=outcar,status='old')
 else
      write(ulog,*)' outcar file ',outcar,' does not exist; check your files location and run again'
      stop
 endif

! first find the number of configurations: ncfg
 t=0
 do j=1,11000000
    read(utraj,'(a)',end=99)line
    call findword('POSITION',line,found)
    if (found) then
       t = t+1
    endif
 enddo
99 write(ulog,*)' reached the end of OUTCAR file; number of configurations= ',t
 close(utraj)
 ncfg = t
 write(*,*)' reached the end of OUTCAR file which had ',ncfg,' configurations '
 if(t.eq.0) then
    write(ulog,*)' the word POSITION was not found in OUTCAR file, check it!'
    stop
 endif

 i=(j-1)/t - 2  ! this is the number of atoms read from OUTCAR
 if (i .ne. natom_super_cell ) then
    write(ulog,*)' number of atoms read .ne. no of atoms in POSCAR file',i,natom_super_cell
    write(ulog,*)' # of read lines in OUTCAR is=',j
    write(ulog,*)' check your POSCAR and OUTCAR again '
    write(ulog,*)' make sure the # of atoms is the same in both files'
    write(ulog,*)' there should be no blank lines at the end of OUTCAR'
    stop
 endif

 end subroutine count_configs
!===========================================================
 subroutine read_force_position_data(outcar,ncfg,energy,dsp,frc)
!! reads contents of OUTCARi; second lines are energies, frc_constr=3*natom_super_cell*ncfg
!! outputs are dsp(3,NSC,ncfg),frc(3,NSC,ncfg),energy(ncfg) for the OUTCARi file read
 use ios
 use atoms_force_constants
 use geometry
 use params
 use lattice
 use svd_stuff
 implicit none
 integer, intent(in) :: ncfg
 character, intent(in):: outcar*(*)
 integer n1,n2,i,t,j,k
 type(vector) v
 real(8) x1,x2,x3,x4,const,rr(3)
 real(8), save :: emin
 logical, save :: first_call=.true.
! real(8), allocatable, intent(inout) :: energy(:),dsp(:,:,:),frc(:,:,:)
 real(8), intent(out) :: energy(ncfg),dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg)

 character line*99
 logical found,exst

 write(*,*)' REopening OUTCAR file'

 inquire(file=outcar,exist=exst)
 if(exst) then
    open(utraj,file=outcar,status='old')
 else
    write(ulog,*)' outcar file ',outcar,' does not exist; check your files location and run again'
    stop
 endif


! now get the FORCES from OUTCAR file
 t=0 ; energy=0
 do j=1,ncfg
    read(utraj,'(a)',end=88)line
!    if (line(1:8) .eq. "POSITION" ) then
    call findword('POSITION',line,found)
    if (found) then
       t = t+1
       read(utraj,*) k,energy(t)  ! start with 1 since t=0
       if(t.ne.j .or. t-1.ne.k) then
          write(*,*)'Error in reading snapshot#s in OUTCAR?',k,j,t
       endif
       do i=1,natom_super_cell
           read(utraj,*) dsp(1:3,i,t),frc(1:3,i,t)
       enddo
       write(*,*)'j,t,k=',j,t,k
    endif
 enddo
88 write(ulog,*)' reached the end of OUTCAR file after ',t,' steps'
 close(utraj)
 write(ulog,*)' last line read was '
 write(ulog,'(a)') line
 close(utraj)
 if (t .ne. ncfg) then
    write(*,*)'ERROR in reading the force file OUTCAR'
    write(*,*)'ncfg, # of read steps=',ncfg,t
    stop
 endif

 write(*,*)'writing frc and dsp in the log file'
 call write_out(ulog,' last force ',frc(:,:,t))
 call write_out(ulog,' last coord ',dsp(:,:,t))


! get energy per primitive unit cell so that it does not depend on supercell size
 energy=energy/natom_super_cell*natom_prim_cell

! subtract lowest energy value ! get it from OUTCAR1; it is arbitrary anyways
 if( first_call ) then
    emin=minval(energy)
    first_call=.false.
 endif
 energy=energy-emin

 write(*,9)'Energy/primcell assigned ',energy
 write(*,*)'Calling calculate_and_write_displacements'
 call calculate_and_write_displacements(ncfg,dsp,frc)
! inputs: displ, including atom_sc%equilibrium_pos 
! outputs: displ - atom_sc%equilibrium_pos  

 write(*,*)'calling write_correspondance'
 call write_correspondance

 write(*,*)'exiting read_force_position_data '
9 format(a,200(1x,f7.3))

 end subroutine read_force_position_data
!===========================================================
 subroutine read_force_position_data2(outcar,frc_constr,ncfg,energy,dsp,frc)
!! reads contents of OUTCARi; second lines are energies, frc_constr=3*natom_super_cell*ncfg
!! outputs are nfcg,dsp(3,NSC,ncfg),frc(3,NSC,ncfg),energy(ncfg) for every outcar file read
 use ios
 use atoms_force_constants
 use geometry
 use params
 use lattice
 use svd_stuff
 implicit none
 integer, intent(out) :: frc_constr,ncfg
 character, intent(in):: outcar*(*)
 integer n1,n2,i,t,j,k
 type(vector) v
 real(8) x1,x2,x3,x4,const,rr(3),emin
! real(8), intent(inout) :: energy(ncfg),dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg)
 real(8), allocatable, intent(inout) :: energy(:),dsp(:,:,:),frc(:,:,:)
! integer, allocatable :: nlines(:)

 character line*99
 logical found,exst

 inquire(file=outcar,exist=exst)
 if(exst) then
      open(utraj,file=outcar,status='old')
 else
      write(ulog,*)' outcar file ',outcar,' does not exist; check your files location and run again'
      stop
 endif

 t=0
 do j=1,11000000
    read(utraj,'(a)',end=99)line
    call findword('POSITION',line,found)
    if (found) then
       t = t+1
    endif
 enddo
99 write(ulog,*)' reached the end of OUTCAR file; number of configurations= ',t
 ncfg = t
 if(t.eq.0) then
    write(ulog,*)' the word POSITION was not found in OUTCAR file, check it!'
    stop
 endif

! allocate( energy(ncfg),dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg) )
! if (ncfg.gt.5000) then
!    write(*,*)'READ_FORCE_POSITION_DATA: reading more than 5000 configurations ',ncfg
!    write(*,*)'either put in less than 5000 configurations or increase the size of energy in this subroutine'
!    stop
! endif

 i=(j-1)/t - 2  ! this is the number of atoms read from OUTCAR
 if (i .ne. natom_super_cell ) then
    write(ulog,*)' number of atoms read .ne. no of atoms in POSCAR file',i,natom_super_cell
    write(ulog,*)' # of read lines in OUTCAR is=',j
    write(ulog,*)' check your POSCAR and OUTCAR again '
    write(ulog,*)' make sure the # of atoms is the same in both files'
    write(ulog,*)' there should be no blank lines at the end of OUTCAR'
    stop
 endif

! allocates displ and force arrays
 if ( allocated(dsp)) deallocate(dsp)
 if ( allocated(frc)) deallocate(frc)
 if ( allocated(energy)) deallocate(energy)
! if ( allocated(nlines)) deallocate(nlines)

! call allocate_pos(natom_super_cell,ncfg) ! for displ and force
 allocate( energy(ncfg),dsp(3,natom_super_cell,ncfg),frc(3,natom_super_cell,ncfg) )
! allocate( nlines(ncfg) )

! now get the FORCES from OUTCAR file
 rewind(utraj)
 t=0 ; energy=0
 do j=1,11000000
    read(utraj,'(a)',end=88)line
!    if (line(1:8) .eq. "POSITION" ) then
    call findword('POSITION',line,found)
    if (found) then
       t = t+1
       read(utraj,*) k,energy(t)  ! start with 1 since t=0
       if(t.ne.j .or. t.ne.k) then
          write(ulog,*)'Error in reading snapshot#s in OUTCAR?',k,j,t
       endif
       do i=1,natom_super_cell
           read(utraj,*) dsp(1:3,i,t),frc(1:3,i,t)
       enddo
!       nlines(t)=natom_super_cell
    endif
 enddo
88 write(ulog,*)' reached the end of OUTCAR file after steps= ',t
 write(ulog,*)' last line read was '
 write(ulog,'(a)') line
 close(utraj)
 if (t .ne. ncfg) then
    write(ulog,*)'ERROR in reading the force file OUTCAR'
    write(ulog,*)'ncfg, # of read steps=',ncfg,t
    stop
 endif

 call write_out(ulog,' last force ',frc(:,natom_super_cell,t))
 call write_out(ulog,' last coord ',dsp(:,natom_super_cell,t))

 frc_constr = ncfg *natom_super_cell*3

! get energy per unit cell so that it does not depend on supercell size
 energy=energy/natom_super_cell*natom_prim_cell

! subtract lowest energy value
 emin=minval(energy)
 energy=energy-emin

 call calculate_and_write_displacements(ncfg,dsp,frc)
! inputs: displ, including atom_sc%equilibrium_pos 
! outputs: displ - atom_sc%equilibrium_pos  

 call write_correspondance

 end subroutine read_force_position_data2
!========================================
 subroutine update_nshells2(ng,grd) !,,x01,x02,x03 )
!! finds the longest vector in grd(3,ng) and resets nshells(2,:) accordingly
 use geometry
 use ios
! use lattice, only : volume_r,volume_r0
 use atoms_force_constants
 use params, only : nshells
 implicit none
! type(vector), intent(in):: x01,x02,x03
 integer, intent(in):: ng
 real(8), intent(in):: grd(3,ng)
 integer i,shel_count,i0,m,i1,i2,i3,insid,cnt
! type(vector) b01,b02,b03,b1,b2,b3
 real(8) lmax
! real(8), allocatable :: aux(:,:)

! allocate(sx(3,26),s0x(3,26))

! find the shortest basis bi among xi and the corresponding 26 shortest vectors
! call get_26shortest_shell(x01,x02,x03,s0x,b01,b02,b03)! s0x for WS of primcell
! call get_26shortest_shell(x1 ,x2 ,x3 ,sx ,b1 ,b2 ,b3) ! sx for WS of supercell

! y's are used to find the reduced coordinates on the basis xi
! call make_reciprocal_lattice_v(x1,x2,x3,y1,y2,y3)

      ! find the number of grid points inside the supercell
!      n=nint(volume_r/volume_r0)
!      write(ulog,4)' grid size,vol0,vol_sc=',n,volume_r,volume_r0

      ! this is an upperbound to the points inside and on the WS boundaries
!      nmax=n+nint(8*n**0.67+12*n**0.3333+24)

      ! aux contains the grid vectors inside or on the boundary
!      allocate(aux(3,nmax))
!      write(ulog,*)' allocated grid size, including boundary points=',nmax

      ! create the grid within the WS supercell xi
!      m=nint(nmax**0.3333) !1+floor(length(sx(:,26))/length(sx(:,1))/2) ! ratio of the longest to smallest lengths
!      cnt=0
!      lmax=0
!      write(ulog,*)' upper index to loop over to create grid =',m
!      do i1=-m,m
!      do i2=-m,m
!      i3loop: do i3=-m,m

      ! generate a grid of vectors from primitive translations x0i
!       v = v2a(dble(i1)*x01 + dble(i2)*x02 + dble(i3)*x03)

      ! see if v is in the WS cell; throw away v outside of WS cell
!       call check_inside_ws(v,v2a(sx),insid)

!       if (insid.eq.0) cycle i3loop

      ! take only those points inside or on the boundary of WS (corresponding to in=1)
!       cnt=cnt+1
!       write(*,4)'new vector ',cnt,v
!       aux(:,cnt)=v
!       if(length(v).gt.lmax) lmax=length(v)
!      enddo i3loop
!      enddo
!      enddo

!      write(ulog,4)' actual grid size=',cnt
! write(ulog,4)' Number of vectors defining the boundary=',nboundary
!      deallocate(aux)



! find lamx, length of the longest vector in grid defined by grd
 lmax=0
 do i=1,ng
    if( length(grd(:,i)) .gt. lmax  ) lmax = length(grd(:,i))
 enddo
 write(ulog,4)' longest vector length in the grid=',lmax

 do i0=1,natom_prim_cell
    shell_loop: do shel_count = 0 , maxneighbors
       if(lmax .myeq. atom0(i0)%shells(shel_count)%rij) then
          nshells(2,i0) = shel_count
          exit shell_loop
       endif
    enddo shell_loop
 enddo

 write(ulog,5)' ****************************************'
 write(ulog,5)' UPDATED SHELLS of rank 2 ',nshells(2,1:natom_prim_cell)
 write(ulog,5)' ****************************************'

4 format(a,9(1x,f13.5))
5 format(a,20(1x,i2))

 end subroutine update_nshells2
!==========================================================
 subroutine find_ws_weights(ngrd,grd,save_boundary,transl,weight)
 use geometry
 implicit none
 integer, intent(in) :: ngrd
 real(8), dimension(3,ngrd), intent(in) :: grd
 real(8), dimension(3,26)  , intent(in) :: transl
 integer, dimension(ngrd)  , intent(in) :: save_boundary
 real(8), dimension(ngrd) , intent(out) :: weight
 real(8) x(3)
 integer i,j,t !ngrd,ntransl,

! ngrd=size(grd(1,:))
 weight=1
! if grid is on the boundary, find to how many boundary points it is mapped by a transl vector
 do i=1,ngrd
    if (save_boundary(i).ne.0) then
       do t=1,26
          x=grd(:,i)+transl(:,t)
          gloop: do j=1,ngrd  ! is x on the boundary? if so add +1 to the weight of i
             if ( length(x-grd(:,j)) .lt. 1d-3 ) then ! x is also a grid point
                if ( save_boundary(j).ne.0 ) then  ! x is also on the boundary
                   weight(i)=weight(i)+1
                endif
             endif
          enddo gloop
       enddo
    endif
 enddo
 weight=1d0/weight
 do i=1,ngrd
    write(*,*)i,' th grid point has weight=',weight(i)
 enddo

 end subroutine find_ws_weights
!===========================================================
 subroutine is_on_boundary(v,grid,indx,tolerance)
!! checks whether the vector v is on the FBZ boundary within the tolerance
!! where FBZ is defined by the grid(3,26)
!! index is the number of boundaries to which v belongs
!! index=0 means not on boundary
 use geometry
 implicit none
 real(8), intent(in) :: v(3),grid(3,26),tolerance
 integer, intent(out) :: indx
 integer j
 real(8) vgoverg2

 indx=0
 do j=1,26
    vgoverg2 = 2 * dot_product(v,grid(:,j))/dot_product(grid(:,j),grid(:,j))
    if ( abs(vgoverg2 - 1) .lt. tolerance) then
       indx=indx+1
    endif
 enddo

 end subroutine is_on_boundary
!============================================================
 subroutine belongs_to_grid(v,n,grid,nj,tolerance)
! checks whether the vector v belongs to the array "grid" within the 1d-4;
! nj is its index in the grid; if nj=0, the vector did not belong to the grid
 use geometry
 implicit none
 integer, intent(in) :: n
 integer, intent(out) :: nj
 real(8), intent(in) :: v(3),grid(3,n),tolerance
 integer j

 nj=0
 do j=1,n
    if(length(v-grid(:,j)) .lt. tolerance) then
       nj=j
       exit
    endif
 enddo

 end subroutine belongs_to_grid
!==========================================================
 subroutine myreduce(q,nq)
! reduces a 3x3 matrix to tridiagonal form by adding/subtracting lines and columns so that the determinant does not change
 implicit none
 integer, intent(out) :: nq(3)
 real(8), intent(in)  :: q(3,3)
 real(8) c(3),d(3)
 logical sing

 call qrdcmp(q,3,3,c,d,sing)

 nq=nint(d)
! nq(1)=nint(q(1,1))
! nq(2)=nint(q(2,2))
! nq(3)=nint(q(3,3))
 end subroutine myreduce
!========================================
      SUBROUTINE qrdcmp(a,n,np,c,d,sing)
      implicit none
      INTEGER n,np
      REAL(8) a(np,np),c(n),d(n)
      LOGICAL sing
      INTEGER i,j,k
      REAL(8) scale,sigma,sum,tau

      sing=.false.
      scale=0.
      do k=1,n-1
        do i=k,n
          scale=max(scale,abs(a(i,k)))
        enddo
        if(scale.eq.0.)then
          sing=.true.
          c(k)=0.
          d(k)=0.
        else
          do i=k,n
            a(i,k)=a(i,k)/scale
          enddo
          sum=0.
          do i=k,n
            sum=sum+a(i,k)**2
          enddo
          sigma=sign(sqrt(sum),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          do j=k+1,n
            sum=0.
            do i=k,n
              sum=sum+a(i,k)*a(i,j)
            enddo
            tau=sum/c(k)
            do i=k,n
              a(i,j)=a(i,j)-tau*a(i,k)
            enddo
          enddo
        endif
      enddo
      d(n)=a(n,n)
      if(d(n).eq.0.)sing=.true.

      END subroutine qrdcmp
 !===========================================================
      subroutine get_allstars_k(kvec,primlatt,narms,kvecstar,kvecop)
!     kvec(i) (input), ith dimensionless component of k vector
!     narms (output), number of star vectors associated with kvec
!     kvecstar(3,1:narms), all the stars of kvec
!     kvecop(i), the symmetry operation number for the star vecor i
      use atoms_force_constants
!      use force_constants_module
      implicit none
      integer i,j,k,m,n,narms,kvecop(48),ier,ncmp
      double precision kvec(3),kvecstar(3,48),primlatt(3,3),v(3),  &
     &     v2(3),v3(3),kvecstarp(3,48),sum
      narms=0
!      print*,'lattpgcount=',lattpgcount

      iloop: do i=1,lattpgcount
! apply symmetry operation to k to get v=kstar=opk(i)*k
!        call xvmlt(op_kmatrix(1,1,i),kvec,v,3,3,3)
        v=matmul(op_kmatrix(:,:,i),kvec)
! find the reduced coordinates of v and store in v2
!        call xmatmlt(v,primlatt,v2,1,3,3,1,3,1)
        v2=matmul(v,primlatt)
! now check if v2 - any_previous_v2 is zero
! if so, skip; else store this v as a new star vector
        do j=1,narms
        ! subtract previous_v2(=kvecstarp) from v2; result is v3
!          call xvsub(v2,kvecstarp(1,j),v3,3)
          v3=v2-kvecstarp(:,j)
! if all v3(i)=0 cycle and go to the next arm, else, it's a new star
          sum=0
          do k=1,3
             sum=sum+abs(v3(k))
!            n=nint(v3(k))
!            if(ncmp(v3(k)-n).ne.0)exit  ! if v3 non-integer skip to next arm
!            if(abs(ncmp(v3(k)).eq.0)exit  ! if it's not new v3 non-integer skip to next arm
!            if(k.eq.3)cycle iloop  ! goto next sym_op iff v3=integer : i.e. G-vector
          enddo
          if(sum.lt.1d-5) cycle iloop
        enddo

! store this v as a new star vector
        narms=narms+1
        kvecstar(1:3,narms)=v(1:3)
        kvecstarp(1:3,narms)=v2(1:3)
        kvecop(narms)=i
      enddo iloop

      end subroutine get_allstars_k
!========================================
 subroutine fold_in_WS_BZ(nk,kp,gshel,foldedk)
!! for the BZ defined by (gshel) this subroutine generates foldedk obtained
!! from folding of kp into the FBZ,
 use lattice , only : r01, r02, r03
 use geometry
 use kpoints  !, only : nkc,kpc
 implicit none
 integer, intent(in) :: nk
 real(8), intent(in) :: kp(3,nk),gshel(3,26)
 real(8), intent(out) :: foldedk(3,nk)
 real(8) prim2cart(3,3)  !,gmax  !cart2prim(3,3) ,
 integer i,inside

 prim2cart(:,1)=r01
 prim2cart(:,2)=r02
 prim2cart(:,3)=r03

 open(134,file='kfolding.dat')
 write(134,*)'# i       kp(:,i)          folded_k(:,i)'
 do i=1,nk
! for each kpoint, if it is not within the FBZ, fold it into FBZ
    call check_inside_ws(kp(:,i),gshel,inside)
    if(inside .eq. 1) then
       foldedk(:,i)=kp(:,i)
    else
! need to fold kp in the FBZ
       call fold_in_bz_new(kp(:,i),gshel,foldedk(:,i))
    endif
    write(134,4)i,kp(:,i),foldedk(:,i)
 enddo
 close(134)
4 format(i6,2(3x,3(1x,f9.4)))

! Now calculate the weigths
 call get_weights3(nk,foldedk,prim2cart,nibz)

 end subroutine fold_in_WS_BZ
!========================================
 subroutine get_kpfbz(gg1,gg2,gg3,gshel,kpfull)
!! this routine translates kp grid into the formal FBZ. Output to KPFBZ.DAT
 use kpoints !, only : nkc,kpc
 use geometry
 implicit none
 type(vector), intent(in) :: gg1(3),gg2(3),gg3(3)
 real(8), intent(in) :: gshel(3,26)
 real(8), intent(out) :: kpfull(3,nkc)
 integer i , inside
 real(8) gp(3,3),q(3)
 integer :: unit = 100

 do i=1,nkc
   call fold_in_bz_new(kpc(:,i),gshel,kpfull(:,i))
 enddo

 open(unit,file='KPFBZ.DAT',status='unknown')
 write(unit,'(A)') &
 & "#   n              kpfbz(:,n)          |kpfbz(:,n)|           kp(:,n)          |kp(:,n)|"
 do i=1,nkc
   write(unit,'(I10,4X,4(1x,g11.5),4X,4(1x,g11.5))') i, &
 &      kpfull(:,i),sqrt(sum(kpfull(:,i)**2)), kpc(:,i),sqrt(sum(kpc(:,i)**2))
 enddo
 close(unit)

 end subroutine get_kpfbz
!========================================
 subroutine get_short_rlatvecs(g1,g2,g3,q)
!! given an arbitrary set of reciprocal lattice vectors (g1,g2,g3), find
!! alternate k-space basis vectors (q(:,1),q(:,2),q(:,3)) as short as possible
!
 use constants        !! pi
 use geometry
 use ios
 implicit none
 real(8), intent(in) :: g1(3),g2(3),g3(3)
 real(8), intent(OUT) :: q(3,3)
 real(8)              :: g(3,3),r(3,3),k(3),v(3),q2(3),q2max,k2,r1(3),r2(3),r3(3)
 integer              :: i,j,i1,i2,i3,index(3)
 real(8) :: delta = 1d-6
 real(8) :: zero2 = 1d-12

 call make_reciprocal_lattice_a(g1,g2,g3,r1,r2,r3)

 g = reshape((/g1(1),g1(2),g1(3),g2(1),g2(2),g2(3),g3(1),g3(2),g3(3)/),(/3,3/))
 r = reshape((/r1(1),r1(2),r1(3),r2(1),r2(2),r2(3),r3(1),r3(2),r3(3)/),(/3,3/))

 q = g
 q2 = sum(q**2,1)
 q2max = maxval(q2)
 write(ulog,*)'GET_SHORTY:q2max=',q2max
 do i=1,3
   index(i) = floor(sum(r(:,i)**2)*q2max/(2d0*pi)+delta)
 enddo
 write(ulog,*)'GET_SHORTY:index=',index

 do i3=0,index(3)
   do i2=-index(2),index(2)
     kpoint: do i1=-index(1),index(1)
       if (i1.eq.0.and.i2.eq.0.and.i3.eq.0) cycle
       k = matmul(g,(/i1,i2,i3/))
       k2 = k.dot.k
       if (k2.ge.q2max) cycle
       ! check if colinear
       do i=1,3
         v = k-q(:,i)*((k.dot.q(:,i))/q2(i))
         if ((v.dot.v).gt.zero2) cycle
         if (k2.ge.q2(i)) cycle kpoint
         q(:,i) = k
         q2(i) = k2
         q2max = maxval(q2)
         cycle kpoint
       enddo
       ! check if coplanar
       do i=1,3
         j = modulo(i,3)+1
         v = k-q(:,i)*((k.dot.q(:,i))/(q2(i)))
         v = v-q(:,j)*((v.dot.q(:,j))/(q2(j)))
         if ((v.dot.v).gt.zero2) cycle
         if (k2.ge.q2(i) .and. k2.ge.q2(j)) cycle kpoint
         if (q2(i).gt.q2(j)) j = i
         q(:,j) = k
         q2(j) = k2
         q2max = maxval(q2)
         cycle kpoint
       enddo
       ! neither colinar or coplanar
       i = maxloc(q2,1)
       q(:,i) = k
       q2(i) = k2
       q2max = maxval(q2)
     enddo kpoint
   enddo
 enddo

 end subroutine get_short_rlatvecs
!==========================================================
 subroutine get_weights3(nk,kp,prim2cart,nibz) !,kibz,wibz)
!! takes as input kp(3,nk) in cartesian and finds based on crystal symmetries, the
!! weights associated with each kpoint, and then normalizes them
!! nibz is the final number of kpoints stored in kibz in the irreducible BZ
!! the other important output is the mapping of the full kmesh onto the ones
!! folded in the irreducible zone : mapibz(j=1:nk) is the index of the k in the IBZ
!! corresponding to the argument j
!! for i=1,nibz mapinv(i) gives the index k of the kpoints generated from
! n1,n2,n3 loops
 use lattice, only : primitivelattice,r01,r02,r03 
 use kpoints, only : kibz, wibz
 use constants
 use geometry
 use params
 use ios
 implicit none
 integer, intent(in):: nk
 real(8), intent(in):: kp(3,nk) ,prim2cart(3,3) !r1(3),r2(3),r3(3)
 integer, intent(out):: nibz
! real(8), intent(out), allocatable :: kibz(:,:) ,wibz(:)
 integer, allocatable :: mcor(:),mapibz(:),mapinv(:)
 real(8), allocatable :: k2(:,:),lg(:),w2(:)
 real(8) zro,q(3),kvecstar(3,48),sw,skc(3),rr1(3),rr2(3),rr3(3)
 integer nkibz,i,j,l,narms,kvecop(48),aux
 logical exists,foundit

 rr1=v2a(r01)/2/pi
 rr2=v2a(r02)/2/pi
 rr3=v2a(r03)/2/pi

 open(uibz,file='KPOINT.IBZ',status='unknown')

 allocate(k2(3,nk),w2(nk),mapibz(nk),mapinv(nk))
 zro=0d0  ! shift
 write(ulog,*)'GET_WEIGHTS3: generating kpoints in the irreducible FBZ '
 write(*,*)'MAKE_KP_IBZ: nk,wk(nk),k2(:,nk)'

 nibz=0 ; mapibz=0; w2=1

! initialize mapvinv with identity so that later elements can be switched
 do i=1,nk
    mapinv(i)=i
 enddo

! main loop to identify points in the FBZ
 kploop: do i=1,nk
!    q = kp(:,i)
! kp is in cartesian coordinates, we need to convert it to reduced units:
!    q(1)=(kp(:,i)  .dot. r1) /2/pi
!    q(2)=(kp(:,i)  .dot. r2) /2/pi
!    q(3)=(kp(:,i)  .dot. r3) /2/pi
! below the cartesian components of kp are needed
    call getstar(kp(:,i),primitivelattice,narms,kvecstar,kvecop)

! see if already exists among the previously assigned ones
    exists=.False.
    if(verbose) write(ulog,4)'list of kvecstar(l),l=1,narms for kp_red=',i,q

    lloop: do l=1,narms

        if(verbose)   write(ulog,4)'stars are:',l,kvecstar(:,l)

! set weight for the first kpoint where nibz=0
        jloop: do j=1,nibz
          if (k2(:,j) .myeq. kvecstar(:,l)) then
! first bring the star in the FBZ, then compare to the existing points
             exists=.True.
             w2(j)=w2(j)+1d0
             mapibz(i)=j
!            write(ulog,4)' this kpoint turned out to exist ',j,k2(:,j)
             exit lloop
          endif
        enddo jloop
    enddo lloop

!   write(ulog,4)' kpoint, folded one  ',i,kp(:,i),exists,j,k2(:,j)
    if(exists ) then
       cycle kploop
    else
       nibz=nibz+1
!  choose the kpoint star in the first quadrant: at least  works for cubic systems
    !  call getstar(kp(:,i),primitivelattice,narms,kvecstar,kvecop)
    !  foundit=.False.
    !  strloop: do l=1,narms
    !     q= kvecstar(:,l)
    !     if ( q(1).ge.q(2) .and. q(2).ge.q(3) .and. q(3).ge.-1d-5) then
    !        k2(:,nibz)=q
    !        foundit=.True.
    !        exit strloop
    !     endif
    !  enddo strloop

       if (nibz.eq.1) w2(nibz) = 1
       mapibz(i) = nibz

! here need to switch two elements of mapinv; so we add this line
       aux=mapinv(nibz)
       mapinv(nibz)=i
       mapinv(i)=aux

    !  if (.not. foundit) then
          k2(:,nibz)=kp(:,i)
          write(ulog,4)'new vector*:',nibz,w2(nibz),k2(:,nibz),length(k2(:,nibz))
          write(*,4)'new vector*:',nibz,w2(nibz),k2(:,nibz),length(k2(:,nibz))
    !  else
    !     write(ulog,4)'new vector :',nibz,w2(nibz),k2(:,nibz),length(k2(:,nibz))
    !     write(*,4)'new vector :',nibz,w2(nibz),k2(:,nibz),length(k2(:,nibz))
    !  endif
    endif
 ! test for FCC
 !  q(:) = kp(:,i)
 !  if ( q(1).ge.q(2) .and. q(2).ge.q(3) .and. q(3).ge.-1d-5) then
 !     write(bzf,3)i,q
 !  endif

 enddo kploop
 write(ulog,*)'GET_WEIGHTS3: generated ',nibz,' points in the irreducible FBZ'

! define kibz and wibz ---------------------------------------------
 if ( allocated(wibz)) deallocate(wibz)
 if ( allocated(kibz)) deallocate(kibz)
 allocate(wibz(nibz),kibz(3,nibz))
 wibz(1:nibz) = w2(1:nibz)
 kibz(:,1:nibz)=k2(:,1:nibz)
! k2 is in reduced units
! do j=1,nibz
!    kibz(:,j)=k2(1,j)*g1+k2(2,j)*g2+k2(3,j)*g3
! enddo

 deallocate(k2,w2)
 sw = sum(wibz(1:nibz))
 wibz = wibz/sw

! sort and write out the kpoints and their weights-------------------------
 allocate(lg(nibz),mcor(nibz))
 do i=1,nibz
    lg(i)=length(kibz(:,i))
 enddo
 call sort(nibz,lg,mcor,nibz)

!
! this sorting and writing is donw in the main routine
! if ( allocated(lg)) deallocate(lg)
! if ( allocated(mcor)) deallocate(mcor)
! allocate(lg(nk),mcor(nk))
! do i=1,nk
!    lg(i)=length(kp(:,i))
! enddo
! call sort(nk,lg,mcor,nk)
! write(fbz,*)'#i,l,mapibz(l),kp(l),length(kp(l)),kibz(mapibz(l)),wibz(mapibz(l)) l=mcor(i)'
! do i=1,nk
!     j=mcor(i)
!!    write(fbz,3)i,kp(:,mcor(i)),length(kp(:,mcor(i))),kibz(:,mapibz(mcor(i))),wibz(mapibz(mcor(i)))
!    write(fbz,2)i,j,mapibz(j),kp(:,j),length(kp(:,j)),kibz(:,mapibz(j)),wibz(mapibz(j))
! enddo
! deallocate(mcor,lg)
! close(fbz)

 write(uibz,*)'#i,kibz(:,i),wibz(i),kibz_reduced,length(kibz(i))'
 open(345,file='mapinv.dat')
 open(346,file='NEWKP.dat')
 write(345,*)'# i,mapinv(i)(i=1,nibz),mapibz(i)'
 write(346,*)'# i,        newkp(i)         newkp_reduced(i)'
 do i=1,nk !,nibz
    if(i.le.nibz) then
       j=mcor(i)
!       q = matmul(transpose(prim2cart),kibz(:,j))
       call reduce(kibz(:,j),rr1,rr2,rr3,q)
       write(uibz,3)i,kibz(:,j),wibz(j),q,length(kibz(:,j))
    endif
    write(345,*)i,mapinv(i),mapibz(i)
!    skc = matmul(transpose(prim2cart),kp(:,i))
!    call reduce(kp(:,i),gshells(:,1),gshells(:,2),gshells(:,3),skc)
    call reduce(kp(:,i),rr1,rr2,rr3,skc)
    write(346,3)i,kp(:,i),skc
 enddo
 close(345)
 close(346)
 close(uibz)

2 format(3i7,2x,3(1x,f9.4),2x,f9.4,2x,3(1x,f9.5),3x,f14.8)
3 format(i7,2x,3(1x,f9.4),2x,f9.5,2x,3(1x,f9.5),3x,f9.5)
4 format(a,i7,2x,f9.5,2x,99(1x,f9.5))
5 format(a,2x,99(1x,f9.5))

 end subroutine get_weights3
!==========================================================
 subroutine get_weights4(nk,kp,ngibz,mapibz,gibz,wgibz)
!! takes as input kp(3,nk) in cartesian and finds based on crystal symmetries, the
!! weights associated with each kpoint, and then normalizes them
!! nibz is the final number of kpoints stored in kibz in the irreducible BZ
!! the other important output is the mapping of the full kmesh onto the ones
!! folded in the irreducible zone : mapibz(j=1:nk) is the index of the k in the IBZ
!! corresponding to the argument j
!! for i=1,nibz mapinv(i) gives the index k of the kpoints generated from
! n1,n2,n3 loops
 use lattice, only : primitivelattice,r01,r02,r03 
 use kpoints, only : kibz, wibz
 use constants
 use geometry
 use params
 use ios
 implicit none
 integer, intent(in):: nk
 real(8), intent(in):: kp(3,nk) !,prim2cart(3,3) !r1(3),r2(3),r3(3)
 integer, intent(out):: ngibz,mapibz(nk)
 real(8), intent(out), allocatable :: gibz(:,:) ,wgibz(:)
 integer, allocatable :: mapinv(:)  !mcor(:),mapibz(:),
 real(8), allocatable :: k2(:,:),lg(:),w2(:)
 real(8) zro,q(3),kvecstar(3,48),sw,skc(3),rr1(3),rr2(3),rr3(3)
 integer nkibz,i,j,l,narms,kvecop(48),aux
 logical exists,foundit


 open(uibz,file='GMESH.IBZ',status='unknown')

 allocate(k2(3,nk),w2(nk),mapinv(nk))
 zro=0d0  ! shift
 write(ulog,*)'GET_WEIGHTS3: generating kpoints in the irreducible FBZ '
 write(*,*)'MAKE_KP_IBZ: nk,wk(nk),k2(:,nk)'

 ngibz=0 ; mapibz=0; w2=1

! initialize mapvinv with identity so that later elements can be switched
 do i=1,nk
    mapinv(i)=i
 enddo

! main loop to identify points in the FBZ
 kploop: do i=1,nk

! below the cartesian components of kp are needed
    call getstar(kp(:,i),primitivelattice,narms,kvecstar,kvecop)

! see if already exists among the previously assigned ones
    exists=.False.
    if(verbose) write(ulog,4)'list of kvecstar(l),l=1,narms for kp_red=',i,q

    lloop: do l=1,narms

        if(verbose)   write(ulog,4)'stars are:',l,kvecstar(:,l)

! set weight for the first kpoint where nibz=0
        jloop: do j=1,ngibz
          if (k2(:,j) .myeq. kvecstar(:,l)) then
! first bring the star in the FBZ, then compare to the existing points
             exists=.True.
             w2(j)=w2(j)+1d0
             mapibz(i)=j
!            write(ulog,4)' this kpoint turned out to exist ',j,k2(:,j)
             exit lloop
          endif
        enddo jloop
    enddo lloop

!   write(ulog,4)' kpoint, folded one  ',i,kp(:,i),exists,j,k2(:,j)
    if(exists ) then
       cycle kploop
    else
       ngibz=ngibz+1
       if (ngibz.eq.1) w2(ngibz) = 1
       mapibz(i) = ngibz
! here need to switch two elements of mapinv; so we add this line
       aux=mapinv(ngibz)
       mapinv(ngibz)=i
       mapinv(i)=aux
       k2(:,ngibz)=kp(:,i)
       write(ulog,4)'new vector*:',ngibz,w2(ngibz),k2(:,ngibz),length(k2(:,ngibz))
       write(*,4)'new vector*:',ngibz,w2(ngibz),k2(:,ngibz),length(k2(:,ngibz))
    endif

 enddo kploop
 write(ulog,*)'GET_WEIGHTS3: generated ',ngibz,' points in the irreducible FBZ'

! define kibz and wibz ---------------------------------------------
 if ( allocated(wgibz)) deallocate(wgibz)
 if ( allocated(gibz)) deallocate(gibz)
 allocate(wgibz(ngibz),gibz(3,ngibz))
 wgibz(1:ngibz) = w2(1:ngibz)
 gibz(:,1:ngibz)=k2(:,1:ngibz)
 deallocate(k2,w2)
 sw = sum(wgibz(1:ngibz))
 wgibz = wgibz/sw

 write(uibz,*)'#i,kibz(:,i),wibz(i),kibz_reduced,length(kibz(i))'
 open(345,file='mapinv.dat')
 open(346,file='NEWKP.dat')
 write(345,*)'# i,mapinv(i)(i=1,nibz),mapibz(i)'
 write(346,*)'# i,        newkp(i)         newkp_reduced(i)'
 do i=1,ngibz
    write(345,*)i,mapinv(i),mapibz(i)
!    skc = matmul(transpose(prim2cart),kp(:,i))
!    call reduce(kp(:,i),gshells(:,1),gshells(:,2),gshells(:,3),skc)
    call reduce(kp(:,i),rr1,rr2,rr3,skc)
    write(346,3)i,kp(:,i),skc
    write(uibz,3)i,gibz(:,i),wgibz(i)
 enddo
 close(345)
 close(346)
 close(uibz)

2 format(3i7,2x,3(1x,f9.4),2x,f9.4,2x,3(1x,f9.5),3x,f14.8)
3 format(i7,2x,3(1x,f9.4),3x,f9.5,2x,3(1x,f9.5),3x,f9.5)
4 format(a,i7,2x,f9.5,2x,99(1x,f9.5))
5 format(a,2x,99(1x,f9.5))

 end subroutine get_weights4
!====================================
 subroutine structure_factor(r,nk,kp,st)
!! takes kpoints within the FBZ as input, and calculates sum_k cos(k.R) for given R
 implicit none
 integer, intent(in) :: nk
 real(8), intent(in) :: r(3),kp(3,nk)
 real(8), intent(out) :: st
 integer i

 st=0
 do i=1,nk
    st=st+cos(dot_product(r,kp(:,i)))
 enddo
 st=st/nk

 write(*,5)'Structure factor for ', r,st
5 format(a,3(1x,f10.4),3x,g11.5)

 end subroutine structure_factor
!====================================
 subroutine K_mesh(g01,g02,g03,g1,g2,g3,gshells,nk,kmesh)
!! given the translation vectors of the primitive and supercell reciprocal space,
!! this subroutine calculates the mesh of translation vectors within the Wigner-Seitz
!! cell of the supercell, including the border lines.
!! kmesh are the linear combination of gi's folded into the WS cell of g0i
 use geometry
 use ios , only : ulog, write_out
! use kpoints , only : kmesh
 use lattice , only : volume_r0,volume_r, make_rg
 use constants , only : pi
 implicit none
 type(vector), intent(in) :: g01,g02,g03,g1,g2,g3
 real(8), intent(in) :: gshells(3,26)
 integer, intent(in) :: nk
 real(8), intent(out) :: kmesh(3,nk)
 real(8), allocatable:: aux(:,:)
 real(8) q(3),nij(3,3),qred(3)
 type(vector) r01,r02,r03,kred,rr1,rr2,rr3
 integer i,j,k,l,mxi,inside,cnt,mxa
 logical exists,is_integer

 mxi=7
 mxa=(2*mxi+1)**3
 allocate(aux(3,mxa))
 aux=0

 open(347,file='KMESH.SC')
 write(347,*)'# large vectors, vol=',volume_r0
 write(347,3)'# ',g01
 write(347,3)'# ',g02
 write(347,3)'# ',g03
 write(347,*)'# small vectors, vol,ratio=',volume_r,volume_r/volume_r0
 write(347,3)'# ',g1
 write(347,3)'# ',g2
 write(347,3)'# ',g3
! close(347)
 write(347,*)'# i , Kmesh(i) , reduced coordinates '

 call make_reciprocal_lattice_v(g1,g2,g3,rr1,rr2,rr3)  ! rr=rsc/2pi  or g0/2pi
! get reduced coordinates of g0i in terms of gi
 call make_rg(g01,g02,g03,rr1,rr2,rr3,nij)  ! nij(i,j)=g0i.rj

 write(ulog,*) 'reduced coordinates of g0i on gi'
 call write_out(ulog,' g01 on gi ',nij(1,:))
 call write_out(ulog,' g02 on gi ',nij(2,:))
 call write_out(ulog,' g03 on gi ',nij(3,:))
 call write_out(ulog,' determinant ',det(nij))
 call write_out(ulog,' volume ratio ',volume_r/volume_r0)



 cnt=0
 do i=-mxi,mxi
 do j=-mxi,mxi
 do k=-mxi,mxi
    q=i*v2a(g1)+j*v2a(g2)+k*v2a(g3)
    call check_inside_ws(q,gshells,inside)
    if(inside .eq. 1) then
       write(ulog,5)"Q is inside; cnt ",cnt,q
! make sure it does not exist
       exists=.false.
       lloop: do l=1,cnt
          if (length(q-aux(:,l)).lt. 1d-8) then
             exists=.true.
             exit lloop
          endif
          call reduce(q-aux(:,l),rr1,rr2,rr3,qred) ! qred is reduced components on gi
!          if ((qred.dot.qred) .lt. 1d-8 ) then
!             exists=.true.
!             exit lloop
!          elseif(is_integer(qred(1)) .and. is_integer(qred(2)) .and. is_integer(qred(3)) ) then
!             exists=.true.
!             exit lloop
!          endif
       enddo lloop
       if(exists) then
          cycle
       else
         if (cnt.ge.mxa) then
          write(*,*)"KMESH: size of array aux exceeded, please increase and try again"
          stop
         else
          cnt=cnt+1
          aux(:,cnt)=q
          write(ulog,5)"KMESH: new vector q,qred=",cnt,q,qred
         endif
       endif
    endif
 enddo
 enddo
 enddo

 if( nk.ne.cnt ) write(ulog,*)'K_MESH: count is NOT equal to nk!!',cnt,nk
! allocate(kmesh(3,nk))
 kmesh(:,1:nk)=aux(:,1:nk)
 deallocate (aux)

! call make_reciprocal_lattice_v(g01,g02,g03,rr1,rr2,rr3)

 do i=1,nk
    call reduce(kmesh(:,i),rr1,rr2,rr3,kred)
    write(347,4) i,kmesh(:,i),kred%x,kred%y,kred%z
 enddo
 close(347)
3 format(a,2(3x,3(1x,g11.4)))
4 format(i8,2(3x,3(1x,g11.4)))
5 format(a,i8,2(3x,3(1x,g10.3)))


 end subroutine K_mesh
!===========================================================
 subroutine check_inside_ws(q,gshel,inside)
!! checks if q is inside the the Wigner-Seitz cell defined by neighborshells
!! vectors gshel, It is considered inside if on the WS facet/edge/vertex
 use geometry
 use constants
 implicit none
 integer, intent(out) :: inside
! type(vector), intent(in):: g1,g2,g3
 real(8), intent(in):: q(3),gshel(3,26)
 real(8) qdotg,gg,sm,gmax
 integer i,j,k,n

 sm = 1d-4  ! a point on the boundary is also counted as inside

 inside = 0

! construct a lattice out of Gs and take the shortest vectors
! this has to be modified for 1D or 2D systems or if one G is much larger than others
 do i=1,26

    qdotg = ( q .dot. gshel(:,i))
    gg = gshel(:,i) .dot. gshel(:,i)
! check the Bragg condition |q.G| < G.G/2
    if(abs(qdotg) .gt. gg/2 + sm) return  ! keep if on the border

 enddo

! qdotg was smaller than all gg/2: q is therefore inside
 inside = 1

 end subroutine check_inside_ws
!===========================================================
 function is_integer(i)
!! checks if the variable is integer
 implicit none
 real(8) i
 logical is_integer

 if (abs(i-nint(i)).lt.1d-5 ) then
    is_integer=.true.
 else
    is_integer=.false.
 endif
 end function is_integer
!===========================================================
 subroutine get_26shortest_shell(gg1,gg2,gg3,gshells,x1,x2,x3)
!! get the 26 nonzero shortest vectors to define the FBZ, and also to be used for check_inside_fbz
!! outputs gshells and the shortest basis x1,x2,x3 which may replace gg1,gg2,gg3
 use geometry
 use ios , only : ulog
 implicit none
 type(vector), intent(in) :: gg1,gg2,gg3
 type(vector), intent(out):: x1,x2,x3  ! the 3 shortest vectors
 real(8)     , intent(out):: gshells(3,26)
 real(8) , allocatable :: aux(:,:),leng(:)
 integer , allocatable :: mcor(:)
 real(8) junk,vol,gt(3),gmax,shortest(3)
 integer size,i,j,k,n,cnt,radius,diameter

! find the longest G
 gmax=max(length(gg1),length(gg2))
 gmax=max(gmax,length(gg3))
 call calculate_volume(gg1,gg2,gg3,vol)
 junk=vol**0.33333   ! typical short length
 radius = (nint(gmax/junk)+1)
 diameter=2*radius+1
 size= diameter**3 - 1   ! to exclude 0
 write(*,*)'26SHORTEST_SHELL: diameter**3-1=size=',size
 allocate(leng(size),aux(3,size),mcor(size))

! construct a lattice out of Gs and take the 27 shortest vectors
! this has to be modified for 1D or 2D systems or if one G is much larger than others
 cnt=0
 do i=-radius,radius
 do j=-radius,radius
 do k=-radius,radius
    if(i*i+j*j+k*k .eq. 0) cycle
    cnt=cnt+1
    gt=i*v2a(gg1)+j*v2a(gg2)+k*v2a(gg3)
    leng(cnt)=length(gt)
    aux(:,cnt)=gt
 enddo
 enddo
 enddo
 write(*,*)'26SHORTEST_SHELL: counted=',cnt

 call sort(cnt,leng,mcor,cnt)

! find the shortest 3 vectors
 x1= a2v(aux(:,mcor(1)))

! make sure x1 and "x2" are independent
 j=2
 vol=0
 do while  (abs(vol).lt.1d-7*length(x1)**2)
    x2 = a2v(aux(:,mcor(j))) ! 2nd shortest linearly indept from x1
    x3 = x1 .cross. x2
    vol=length(x3)   ! thisis really the area
    write(*,*)'26SHORTEST_SHELL2: j, area=',j,vol
    j=j+1
 enddo

 vol=0
 do while  (abs(vol).lt.1d-10*length(x1)**3)
    x3 = a2v(aux(:,mcor(j)))
    call calculate_volume(x1,x2,x3,vol)
    write(*,*)'26SHORTEST_SHELL3: j, volume=',j,vol
    j=j+1
 enddo
 call calculate_volume(x1,x2,x3,vol)
 write(*,3)'26SHORTEST_SHELL: final smallest vector1=',x1
 write(*,3)'26SHORTEST_SHELL: final smallest vector2=',x2
 write(*,3)'26SHORTEST_SHELL: final smallest vector3=',x3
 write(*,*)'26SHORTEST_SHELL: final smallest volume=',vol

3 format(a,9(1x,g11.4))
4 format(i6,9(1x,g11.4))

 deallocate(leng,aux,mcor)
 allocate(leng(26),aux(3,26),mcor(26))
 cnt=0
 do i=-1,1
 do j=-1,1
 do k=-1,1
    if(i*i+j*j+k*k .eq. 0) cycle
    cnt=cnt+1
    gt=i*(x1)+j*(x2)+k*(x3)
    leng(cnt)=length(gt)
    aux(:,cnt)=gt
 enddo
 enddo
 enddo
 call sort(cnt,leng,mcor,26)
 write(ulog,*)' 26 shortest shell vectors are :'
 do i=1,cnt
    gshells(:,i)=aux(:,mcor(i))
    write(ulog,4)i,gshells(:,i), length(gshells(:,i))
 enddo

 deallocate(leng,aux,mcor)

 end subroutine get_26shortest_shell
 !============================================================
 function bring_to_prim_cell_cav(a) result(w)
! takes cartesian coordinates and returns cart coordinates within the primcell
 use lattice
 use geometry
 implicit none
 type(vector) v,w, bring_to_prim_cell_cv
 real(8) a(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w =  bring_to_prim_cell_cv(v)

 end function bring_to_prim_cell_cav
!-----------------------------------------
 function bring_to_super_cell_cav(a) result(w)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
 use geometry
 implicit none
 type(vector) v,w, bring_to_super_cell_cv
 real(8) a(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w = bring_to_super_cell_cv(v)

 end function bring_to_super_cell_cav
!-----------------------------------------
 function bring_to_prim_cell_caa(a) result(b)
! takes cartesian coordinates and returns cart coordinates within the primcell
 use lattice
 use geometry
 implicit none
 type(vector) v,w, bring_to_prim_cell_cv
 real(8) a(3),b(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w =  bring_to_prim_cell_cv(v)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end function bring_to_prim_cell_caa
!-----------------------------------------
 function bring_to_super_cell_caa(a) result(b)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
 use geometry
 implicit none
 type(vector) v,w, bring_to_super_cell_cv
 real(8) a(3),b(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w = bring_to_super_cell_cv(v)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end function bring_to_super_cell_caa
!-----------------------------------------
 function bring_to_prim_cell_cv(v) result(w)
! takes cartesian coordinates and returns cart coordinates within the primcell
 use lattice
 use geometry
 implicit none
 type(vector) v,w
 real(8) a1,a2,a3

! get direct coordinates first
 a1 = v .dot. g01
 a2 = v .dot. g02
 a3 = v .dot. g03
! bring into cell: between 0 and cell
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
! convert to cartesian coordinates
 w = a1*r01+a2*r02+a3*r03

 end function bring_to_prim_cell_cv
!-----------------------------------------
 function bring_to_super_cell_cv(v) result(w)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
 use geometry
 use constants, only : pi
 implicit none
 type(vector) v,w
 real(8) a1,a2,a3

! get direct coordinates first
 a1 = (v .dot. g01)/(2*pi)
 a2 = (v .dot. g02)/(2*pi)
 a3 = (v .dot. g03)/(2*pi)
! bring into cell: between 0 and cell
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
! convert to cartesian coordinates
 w = a1*r01+a2*r02+a3*r03

 end function bring_to_super_cell_cv
!-----------------------------------
 subroutine cart_to_direct_v(v,w)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
 use geometry
 use constants, only : pi
 implicit none
 type(vector) v,w

 w%x = (v .dot. g01)/(2*pi)
 w%y = (v .dot. g02)/(2*pi)
 w%z = (v .dot. g03)/(2*pi)

 end subroutine cart_to_direct_v
!-----------------------------------
 subroutine cart_to_direct_av(a,w)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
 use geometry
 implicit none
 real(8) a(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call cart_to_direct_v(v,w)

 end subroutine cart_to_direct_av
!-----------------------------------
 subroutine cart_to_direct_aa(a,b)
! takes cart coordinates and returns cart coordinates within the supercell
 use geometry
 use lattice
 implicit none
 real(8) a(3),b(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call cart_to_direct_v(v,w)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end subroutine cart_to_direct_aa
!-----------------------------------
 subroutine direct_to_cart_v(v,w)
! takes direct coordinates and returns cart coordinates
 use geometry
 use lattice
 implicit none
 type(vector) v,w

 w%x = v%x*r01%x + v%y*r02%x + v%z*r03%x
 w%y = v%x*r01%y + v%y*r02%y + v%z*r03%y
 w%z = v%x*r01%z + v%y*r02%z + v%z*r03%z

 end subroutine direct_to_cart_v
!-----------------------------------
 subroutine direct_to_cart_av(a,w)
! takes direct coordinates and returns cart coordinates
 use geometry
 use lattice
 implicit none
 real(8) a(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call direct_to_cart_v(v,w)

 end subroutine direct_to_cart_av
!-----------------------------------
 subroutine direct_to_cart_aa(a,b)
! takes direct coordinates and returns cart coordinates
 use geometry
 use lattice
 implicit none
 real(8) a(3),b(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call direct_to_cart_v(v,w)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end subroutine direct_to_cart_aa
!------------------------------------
 subroutine dir2cart_g(q,k)
! takes q in direct coordinates and outputs k in cartesian
 use geometry
 use lattice
 real(8) q(3),k(3)

 k(:) = q(1)*g01 + q(2)*g02 + q(3)*g03

 end subroutine dir2cart_g
!-------------------------------------
 subroutine find_map_rtau2sc(nr,rmesh,map_rtau_sc)
!! finds the index of the supercell atom corresponding to primitive translation vector
!! defined in rmesh (modulo supercell translations) and primitive atom tau
!! map_rtau_sc(tau,igrid)=k=atom index in supercell 
 use lattice , only : gs1,gs2,gs3 !rs1,rs2,rs3,
 use ios , only : ulog
 use geometry
 use constants, only : pi
 use atoms_force_constants, only : atom_sc, atom0, natom_super_cell,natom_prim_cell
 implicit none
 integer, intent(in) :: nr
 real(8), intent(in) :: rmesh(3,nr)
 integer, intent (out) :: map_rtau_sc(natom_prim_cell,natom_super_cell/natom_prim_cell)
 integer i,tau,k,nsc,l
 type(vector) diff
 real(8) dist(3),dred(3)
 logical is_integer

 kloop:do k=1,natom_super_cell
    tau= atom_sc(k)%cell%tau
    diff= atom_sc(k)%equilibrium_pos - atom0(tau)%equilibrium_pos
! see if diff is a translation vector , and which rmesh vector is it modulo rsc
! first fold it into the supercell WS cell, then compare to rmesh
    nsc=0
    rloop: do i=1,nr
       dist=v2a(diff) - rmesh(:,i)
!       call reduce_v2(dist,gs1,gs2,gs3,dred)
       call reduce(dist,gs1,gs2,gs3,dred)  ! find dist in units of supercell
       dred=dred/(2*pi)
       if(is_integer(dred(1)) .and. is_integer(dred(2)) .and. is_integer(dred(3)) ) then
! make sure not previously assigned
           do l=1,nsc
              if(k.eq.map_rtau_sc(tau,nsc)) cycle rloop
           enddo
           if (nsc.eq.natom_super_cell/natom_prim_cell) then
              write(*,*)'find_map_rtau2sc: atoms all assigned, nsc=',nsc,' exceeding its limit'
              stop
           endif
           nsc=nsc+1
           map_rtau_sc(tau,nsc)=k  ! what's the use of this mapping array? meaning of nsc?
           write(ulog,4)'FIND_MAP_RTAU2SC: tau,isc,d(isc;i0,tau)=',nsc,tau,k,diff
       endif
    enddo rloop
!    write(*,3)'FIND_MAP_RTAU2SC: dred=',dred
!    dred=mod(dred,1d0)
!    write(*,4)'FIND_MAP_RTAU2SC: tau,igridr,isc,rsc_tau-r0_tau[SC]=',tau,i,k,dred
!    if(dot_product(dred,dred) .lt. 1d-8) then
!       map_rtau_sc(tau,i)=k
 enddo kloop

3 format(a,9(1x,f9.4))
4 format(a,3(i7),6(1x,f9.4)) !1x,f9.4))

 end subroutine find_map_rtau2sc
!============================================================
 subroutine invsort(n,r,mcor)
! sorts the first n elements of array r in descending order r(mcor(i)) is the ordered array
! r(mcor) is the ordered array r(mcor(1)) > r(mcor(2)) > ... > r(mcor(n))
  implicit none
  integer n,i,j,temp
  real(8) r(n)
  integer mcor(n)

  do i=1,n
     mcor(i)=i
  enddo
  do i=1,n-1
     do j=i+1,n
        if(r(mcor(i)).lt.r(mcor(j))) then
           temp=mcor(i)
           mcor(i)=mcor(j)
           mcor(j)=temp
        endif
     enddo
  enddo
 end subroutine invsort
!===================================================
 subroutine project_on(ndim,nbasis,v,basis,vp)

!! projects the vector v on of the space spanned by the set of n orthonormal
!! vectors basis (if not the basis has to be orthonormalized first)
!! formula is: vpo=sum_{c=1,n} (v.dot.basis) basis
 implicit none
 integer, intent(in) :: ndim,nbasis
 real(8), intent(in) :: basis(ndim,nbasis),v(ndim)
 real(8), intent(out) :: vp(ndim)
 integer c,i,j
 real(8) res
 
 ! make sure basis is orthonormal
  do i=1,nbasis
  do j=i,nbasis
     res=dot_product(basis(:,i),basis(:,j))
     if (i.eq.j .and. abs(res-1d0).gt.1d-5 ) then
        write(*,*)'Basis not normalized,i,vi^2=',i,res
     elseif ( abs(res).gt.1d-5 ) then
        write(*,*)'Basis not normalized,i,j,vi.vj=',i,j,res
     endif
  enddo
 enddo

 vp=0
 do c=1,nbasis
    vp = vp + dot_product(v,basis(:,c)) * basis(:,c)
 enddo

 end subroutine project_on
!===================================================
 subroutine project_out(ndim,n,v,basis,vpo)
!! projects out the vector v out of the space spanned by the set of n orthonormal
!! vectors basis (if not the basis has to be orthonormalized first)
!! formula is: vpo=v - sum_{c=1,n} (v.dot.basis) basis
 implicit none
 integer, intent(in) :: ndim,n
 real(8), intent(in) :: basis(ndim,n),v(ndim)
 real(8), intent(out) :: vpo(ndim)
 integer c

 vpo=v
 do c=1,n
    vpo = vpo - dot_product(v,basis(:,c)) * basis(:,c)
 enddo

 end subroutine project_out
!===================================================
 subroutine find_first_nonzero(dim,bmat,n_hom)
!! finds the position of the first non-zero element in array bmat
 implicit none
 integer, intent(in):: dim
 real(8), intent(in) :: bmat(dim)
 integer, intent(out) :: n_hom
 integer i

 n_hom=0
 do i=1,dim
    if(bmat(i).eq.0) then
       n_hom=n_hom+1
    else
       exit
    endif
 enddo
 
 end subroutine find_first_nonzero
 !===================================================
 subroutine pad_array(na,a,nb,b,c)
!! extends array a by appending array b to the end of it 
 implicit none
 integer, intent(in) :: na,nb
 real(8), intent(in) :: a(na)
 real(8), intent(in) :: b(nb)
 real(8), allocatable, intent(inout) :: c(:)
 real(8), allocatable :: aux(:)
 integer ncc

 ncc=na+nb
! allocate(aux(ncc))
! aux=reshape(a,shape=(/naa+nbb/),pad=b)
 if(allocated(c)) deallocate(c)
 allocate(c(ncc))
 c(1:na)=a(1:na); c(na+1:na+nb)=b(1:nb)
! return

 end subroutine pad_array
!===================================================
 subroutine pad_array2(a,b)
!! appends array b to the end of a and stores the result in a
 implicit none
 real(8), allocatable, intent(inout) :: a(:)
!real(8), allocatable, intent(in) :: b(:)
!real(8), intent(inout) :: a(:)
 real(8), intent(in) :: b(:)
 real(8), allocatable :: aux(:)
 integer ncc,naa,nbb

 naa=size(a);
 nbb=size(b);
 ncc=naa+nbb

 allocate(aux(ncc)) ! this is to allow calls like append_array(a,b,a)
 aux=reshape(a,shape=(/naa+nbb/),pad=b)
! aux=reshape(a,(/naa+nbb/),pad=b)
! if (allocated(c)) 
 deallocate (a)
 allocate(a(ncc))
 a=aux
 deallocate(aux)

 end subroutine pad_array2
!===================================================
 subroutine append_array(a,b,c)
!! appends array b to the end of a and stores the result in c
 implicit none
! integer, intent(in) :: na,nb
! real(8), intent(in):: a(na),b(nb)
 real(8), intent(in):: a(:),b(:)
! real(8), allocatable, intent(inout):: c(:)
 real(8), allocatable :: c(:)
! real(8) :: c(size(a)+size(b))
 real(8), allocatable :: aux(:)
 integer ncc,naa,nbb

 naa=size(a);
 nbb=size(b);
 ncc=naa+nbb
! if (allocated(c)) deallocate (c)
 allocate(aux(ncc)) ! this is to allow calls like append_array(a,b,a)
 aux=reshape(a,shape=(/naa+nbb/),pad=b)
! aux=reshape(a,(/naa+nbb/),pad=b)
 if (allocated(c)) deallocate (c)
 allocate(c(ncc))
 c=aux
 deallocate(aux)
! c=reshape(a,(/na+nb/),pad=b,order=(/1,1/))


 end subroutine append_array
 

