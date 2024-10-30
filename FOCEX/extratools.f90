!======================================================================r
 subroutine findword(word,line,found,i)
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
!=======================================================================
 subroutine nsc_from_rgrid(ngrid,rgrid,nsc)
!! for rgrid in the WS of supercell and tau in primcell, finds the corresponding 
!! supercell atom index nsc: maps (jgrid,tau) to nsc , may repeat as ngrid*tau > nsc
 use geometry
 use lattice, only : g01,g02,g03,cart2red_r
 use constants, only : r15, pi
 use atoms_force_constants, only : natom_prim_cell
 implicit none
 integer, intent(in) :: ngrid
 real(r15), intent(in) :: rgrid(3,ngrid)
 integer, intent(out) :: nsc(natom_prim_cell,ngrid)
 integer j,n1,n3(3),tau
 real(r15) rr(3)

   do j=1,ngrid
      rr=cart2red_r(rgrid(:,j)) 
  !   rr(1)=(rgrid(:,j) .dot. g01)/2/pi  ! find first the reduced coordinates
  !   rr(2)=(rgrid(:,j) .dot. g02)/2/pi
  !   rr(3)=(rgrid(:,j) .dot. g03)/2/pi
      n3=nint(rr)
      if(length(rr-n3).gt.1d-4) then
         write(*,5)' NSC_FROM_RGRID: n3.ne.rr ',n3,rr
         stop
      endif
      do tau=1,natom_prim_cell
         call findatom_sc(n3,tau,n1)
         nsc(tau,j)=n1
      enddo
   enddo

5 format(a,3(i4),3(1x,f10.5))
 end subroutine nsc_from_rgrid
!=======================================================================
 subroutine findatom_sc(n3,tau,iatomsc)
!! finds the index of the atom in supercell with identity (n3[modulo supercell],tau)
!! arguments:
!!     n3(3) (input), linear combination of basis vectors of the primitive
!!          lattice that takes us to the unit cell containing the ith atom
!!     tau (input), identity of equivalent atom in unit cell at origin
!!     iatomsc (output), index of equivalent atom in supercell. Returns zero if not found
      use geometry
      use lattice
      use ios
      use atoms_force_constants
      use constants , only:pi,r15
      implicit none
      integer, intent(in) :: n3(3),tau
      integer, intent(out) :: iatomsc
      integer j,m(3)  !,al
      real(r15) mi(3)  !,b(3)

      iatomsc = 0
! write(*,*) '************** entering findatom_sc with tau,n=',tau,n3
      jloop: do j=1,natom_super_cell
! write(*,*)'FINDATOM_SC:j_sc,tauj_sc,nj_sc=',j,atom_sc(j)%cell%tau,atom_sc(j)%cell%n
         if (atom_sc(j)%cell%tau .eq. tau ) then
             m = atom_sc(j)%cell%n - n3
             mi=matmul(r0g,m)  ! m is either 0 or an integer (modulo supercell)
             if (is_integer(mi)) then
                 iatomsc=j
                 return
             endif
         endif

!            m=floor(mi+0.001)
!            b=mi-m
!            if(length(b).lt.1d-3) then
!                iatom=j
!                return
!            endif
!        endif
      enddo jloop

      if(iatomsc.eq.0) then
         write(*,*)'FINDATOM_SC: no SC atom could be found for tau,n=',tau,n3
         stop
      endif

4 format(a,9(1x,i3))
5 format(a,9(1x,f9.4))

 end subroutine findatom_sc
!=======================================================================
 subroutine findatom_sc2(n3,tau,iatom)
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
      use constants , only:pi,r15
      implicit none
      integer, intent(in) :: n3(3),tau
      integer, intent(out) :: iatom
      integer j,m(3),al
      real(r15) mi(3)
!     logical is_integer

      iatom = 0
 write(*,*) '************** entering findatom_sc2 with tau,n=',tau,n3
      jloop: do j=1,natom_super_cell
 write(*,*)'FINDATOM_SC2:j_sc,tauj_sc,nj_sc=',j,atom_sc(j)%cell%tau,atom_sc(j)%cell%n
     !   cycle jloop
! find n,tau of atom j and compare to the input n3,tau, modulo nr
         if (atom_sc(j)%cell%tau .eq. tau ) then
             m = atom_sc(j)%cell%n - n3
!write(*,4) ' j,tau,m=',j,tau,m
! bring it to the supercell i.e. make it modulo n_sc : use invn_sc matrix
             mi=matmul(invn_sc,m)
 write(*,*)'   FINDATOM_SC2:m,mi=',m,mi
!write(*,5) ' invrsem=',mi
! mi are the coeffs of m vector in units or n_sc,
! if it has integer components then n3 is good
            do al=1,3
               if(is_integer(mi(al))) then
 write(*,*)'   FINDATOM_SC2:OK for al,mi(al) =',al,mi(al)
                  cycle
               else
                  cycle jloop  ! this mi is no good, try the next atom
               endif
            enddo   ! if we exit the loop here, then mi is good
            iatom=j
 write(*,*)'   FINDATOM_SC2:OK for atom j, exiting ',j
            exit jloop
         endif
      enddo jloop

      if(iatom.eq.0) then
         write(*,*)'FINDATOM_SC2: no SC atom could be found for tau,n=',tau,n3
         stop
      endif

4 format(a,9(1x,i3))
5 format(a,9(1x,f9.4))

 end subroutine findatom_sc2
!===================================================
 subroutine get_n_tau_r(r,n,tau)
!! for a vector position r, finds the (tau,n) FROM ATOMPOS list
!! if r(3) is not in the atompos list, then tau=0
 use constants, only : r15
 use atoms_force_constants
 use params, only : tolerance
 integer, intent(out) :: n(3),tau
 real(r15), intent(in) :: r(3)
 integer i

 tau=0
 do  i=1,natoms
     if(length(r-atompos(:,i)) .lt. tolerance) then
        tau=iatomcell0(i)
        n  =iatomcell(:,i)
        exit
     endif
 enddo

 end subroutine get_n_tau_r
!===================================================
 subroutine get_n_tau_r_mod(r,n,tau)
!! for a vector position r, finds the (tau,n) regardless of ATOMPOS list
!! assumes atoms in atompos are of reduced units between 0 and 1 times (r01,r02,r03)
! just brings it to the primcell and finds the tau.
! It is not the (tau,n) defining atompos!! or is it??
 use constants, only : r15
 use atoms_force_constants
 use params, only : tolerance
 use lattice , only : g01,g02,g03,check_int,prim_to_cart,cart_to_prim !,r01,r02,r03
 use geometry
 integer, intent(out) :: n(3),tau
 real(r15), intent(in) :: r(3)
 integer i0,ier
 real(r15) a(3),v(3),p(3)
 external unitcell

! bring it to unitcell, i.e. between 0 and r0i 
 call unitcell(cart_to_prim,prim_to_cart,r,p)
 v=r-p      ! by construction, this should be an integer multiple of r0i
 n = nint(matmul(cart_to_prim,v))   ! these are its reduced coordinates
 tau=0
 do  i0=1,natom_prim_cell
     a=(p-atompos(:,i0))
     if( length(a) .lt. tolerance) then
        tau=i0  ! i=iatomcell0(i) for i =< natom_prim_cell
        return
!     else
!        write(*,5)'i0,r_cell-pos(i),red=',dble(i0),a,matmul(cart_to_prim,a)
     endif
 enddo
 write(*,4)'GET_N_TAU_R_MOD: atom not found; tau,n,r,r_red=',tau,n,r,matmul(cart_to_prim,r)
 return

 call check_int(a2v(r),a,ier,g01,g02,g03)  !  a is r in reduced units
 n=floor(a)
 v=a-n  ! v components are in [0,1[

 tau=0
 iloop:do  i0=1,natom_prim_cell
     p=matmul(cart_to_prim,v2a(atom0(i0)%equilibrium_pos))  ! reduced coordinates of i0
     if( v .myeq. p) then  
        tau=i0  ! i=iatomcell0(i) for i =< natom_prim_cell
        exit iloop
     else
        write(*,5)'i0,r_red,pos(i0)_red,diff,n=',dble(i0),a,p,a-p,dble(n)
     endif
 enddo iloop

4 format(a,4i4,99(1x,f9.3))
5 format(a,99(1x,f9.3))
 end subroutine get_n_tau_r_mod
!===================================================
 subroutine get_n_tau(j_sc,n,tau)
!! for a supercell atom j_sc, finds the (tau,n)
 use atoms_force_constants
 integer, intent(out) :: n(3),tau
 integer, intent(in) :: j_sc

 n=atom_sc(j_sc)%cell%n
 tau=atom_sc(j_sc)%cell%tau

 end subroutine get_n_tau
!======================================================================
  function fcut(x,x1,x2) result(y)
 use constants, only : r15
  implicit none
  real(r15), intent(in) :: x,x1,x2
  real(r15) y

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
 use constants, only : r15
 implicit none
 integer, intent(in):: i,j,maxatoms
 real(r15), intent(in) :: atompos(3,maxatoms)
 real(r15), intent(out) :: rij

  rij = sqrt( (atompos(1,i)-atompos(1,j))*(atompos(1,i)-atompos(1,j)) +  &
&             (atompos(2,i)-atompos(2,j))*(atompos(2,i)-atompos(2,j)) +  &
&             (atompos(3,i)-atompos(3,j))*(atompos(3,i)-atompos(3,j)) )

 end subroutine calculate_distance
!============================================================
 subroutine sort(n,r,mcor,maxat)
! sorts the first n elements of array r in ascending order r(mcor(i)) is the ordered array
! make sure to initialize the array r to a huge number if n\=maxat
! r(mcor(1)) < r(mcor(2)) < r(mcor(3)) < ... < r(mcor(n))
 use constants, only : r15
  implicit none
  integer,intent(in):: maxat,n
  integer i,j,temp
  real(r15), intent(in) :: r(maxat)
  integer, intent(out) :: mcor(maxat)

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
 use constants, only : r15
 use geometry
 use ios
 implicit none
 real(r15), intent(in) :: r1(3),r2(3),r3(3)
 real(r15), intent(out):: g1(3),g2(3),g3(3)
 real(r15) om

 om = (r1 .dot. (r2 .cross. r3))
 if (abs(om).lt.1d-8) then
    write(*,*)'MAKE_RECIPROCAL_LATTICE_A: volume is zero; check your translation vectors'
    stop
 endif
 g1=(r2 .cross. r3)/om
 g2=(r3 .cross. r1)/om
 g3=(r1 .cross. r2)/om
 om=abs(om)

 end subroutine make_reciprocal_lattice_a
!===============================================================
 subroutine make_reciprocal_lattice_v(r1,r2,r3,g1,g2,g3)
! be careful it does not have the factor of 2*pi in it!
 use geometry
 use constants
 use ios
 implicit none
 type(vector),intent(in) :: r1,r2,r3
 type(vector),intent(out) :: g1,g2,g3
 real(r15) om

 om = (r1 .dot. (r2 .cross. r3))
 if (abs(om).lt.1d-8) then
    write(*,*)'MAKE_RECIPROCAL_LATTICE_V: volume is zero; check your translation vectors'
    write(*,*)'r1=',r1
    write(*,*)'r2=',r2
    write(*,*)'r3=',r3
    stop
 endif
! write(ulog,*)'RECIPROCAL_LATTICE: '
 g1=(r2 .cross. r3)/om
 g2=(r3 .cross. r1)/om
 g3=(r1 .cross. r2)/om
 om=abs(om)

 end subroutine make_reciprocal_lattice_v
!===============================================================
 subroutine make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
! be careful it does not have the factor of 2*pi in it!
 use geometry
 use constants
 use ios
 implicit none
 type(vector) :: r1,r2,r3,g1,g2,g3
 real(r15) om

 om = (r1 .dot. (r2 .cross. r3))
 if (abs(om)/length(r1).lt.1d-8) then
    write(*,*)'MAKE_RECIPROCAL_LATTICE_2PI: volume is zero; check your translation vectors'
    stop
 endif
! write(ulog,*)'RECIPROCAL_LATTICE: '
 g1=2*pi*(r2 .cross. r3)/om
 g2=2*pi*(r3 .cross. r1)/om
 g3=2*pi*(r1 .cross. r2)/om
 om=abs(om)

 end subroutine make_reciprocal_lattice_2pi
!==========================================================
 function delta_k(i,j) result(k)  ! kronecker delta for integers
 use constants
 implicit none
 integer i,j,k

 if (i.eq.j) then 
   k=1
 else
   k=0
 endif

 end function delta_k
!==========================================================
 function delta(x) result(y)
 use constants
 implicit none
 real(r15) x,y

 y = exp(-x*x/2d0)/sqrt(2*pi)

 end function delta
!==========================================================
 function delta_g(x,w) result(y)
 use constants
 implicit none
 real(r15) x,y,w

 y = exp(-x*x/2d0/w/w)/sqrt(2*pi)/w

 end function delta_g
!==========================================================
 function delta_l(x,w) result(y)
 use constants
 implicit none
 real(r15) x,y,w

 y = w/(x*x+w*w)/pi

 end function delta_l
!===========================================================
 function nbe(w,t,m) result(y)
 use constants, only : r15
 implicit none
 real(r15) , intent(in) :: w,t
 integer, intent(in):: m
 real(r15) y,z,x

 x=w
 if(t.le. 0) then
   write(*,*)'ERROR in N_BE: t=<0 ',t
   stop
 endif
 if(w.lt. 0) then
   write(*,*)'ERROR in N_BE: w=<0 ',w
   stop
 elseif(w.eq.0) then
   x=0.001*t
 endif

 z = x/t
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
 use constants, only : r15
    implicit none
    real(r15) :: x1,x2,root,fcn,error,fcnroot,fx1
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
 subroutine check_inside_irred_fbz2d(k,inside)
! assuming k is inside FBZ, this checks whether it is inside the irreducible FBZ
! the case of FCC and hexagonal are treated here.
 use lattice
 use constants
 implicit none
 real(r15) , intent(in) :: k(3)
 real(r15) kdg,q2
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
 end subroutine check_inside_irred_fbz2d
!===========================================================
 subroutine check_inside_irred_fbz(k,inside)
! assuming kp is inside FBZ, this checks whether it is inside the irreducible FBZ of FCC
 use lattice
 use constants
 implicit none
 real(r15) k(3)
 logical, intent(out) :: inside

 if( (abs(k(3)).lt.1d-9 .or. (k(3) >= 0)) .and. (k(3) <= k(2)) .and. (k(2) <= k(1)) &
&    .and. (k(1) <= boxg(1))  .and.  (k(1)+k(2)+k(3) <= 3d0*boxg(1)/2d0) ) then
    inside = .true.
 else
    inside = .false.
    write(444,*)k
 endif

 end subroutine check_inside_irred_fbz
!===========================================================
 function indexg(i,j,k,nil,nih,njl,njh,nkl,nkh) result(n)
! finds the index n of the kpoint defined with 3 loop indices
! ijk going in general from nil to nih, njl to njh nkl to nkh
! n=0; do i1=nil,nih; do j=njl,njh; do k=nkl,nkh; n=n+1
 implicit none
 integer, intent(in) :: i,j,k,nil,nih,njl,njh,nkl,nkh
 integer n

 n = (k-nkl+1) + (j-njl)*(nkh-nkl+1) + (i-nil)*(njh-njl+1)*(nkh-nkl+1)

 end function indexg
!============================================================
 function indexn(i,j,k,n1,n2,n3) result(n)
! finds the index n of the coarse kpoint defined with 3 loop indices
! ijk going from -ni to ni
 implicit none
 integer, intent(in) :: i,j,k,n1,n2,n3
 integer n

 n = (k+n3+1) + (j+n2)*(2*n3+1) + (i+n1)*(2*n2+1)*(2*n3+1)

 end function indexn
!============================================================
 function index_reg(i,j,k,n1,n2,n3) result(n)
! finds the index n of the kpoints defined with 3 loop indices
! ijk going from 1 to ni , i being the outer loop and k the inner one
 implicit none
 integer, intent(in) :: i,j,k,n1,n2,n3
 integer n

 n = k + (j-1)*n3 + (i-1)*n2*n3

 end function index_reg
!===========================================================
 subroutine make_sorted_gs(g1,g2,g3,nshell,gg)
! from the basis (g1,g2,g3) generate the nshell shortest linear combinations including 0
 use geometry
 use ios
 use params
 use constants, only : r15
 implicit none
 integer, intent(in) ::  nshell
 type(vector), intent(in) :: g1,g2,g3
 real(r15), intent(out) :: gg(3,nshell)
 integer i,j,k,n5,ik,ns
 real(r15), allocatable :: g5(:,:)

 ns = nint((2*nshell+1)**(0.3333)) ! this is larger than the radius corresponding to nshell points
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
 use constants, only : r15
 implicit none
 integer nshell,i
 real(r15) gg(3,nshell)
 real(r15) g5(3,nshell),gs(nshell)
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
 use constants, only : r15
 use geometry
 implicit none
 real(r15), intent(in):: q(3)
 type(vector),intent(in):: g1,g2,g3
 type(vector)rr1,rr2,rr3
! integer, intent(in):: n(3)
 integer, intent(out):: inside
 real(r15), intent(out):: qx,qy,qz
 real(r15) epsl

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
! inside would be =1 if 0=<q<1
 use constants, only : r15
 use geometry
 implicit none
 integer, intent(in):: n
 integer, intent(out):: i,insid
 real(r15),intent(in) :: q(3)
 type(vector), intent(in) :: rr
 real(r15) aux,qdr,epsl

 epsl=5d-8

 qdr = (q.dot.rr) + epsl
! include from -0.00005 to 0.99995 included (include 0.0 and exclude 1.0)
 if (qdr.lt.0 .or. qdr.ge.1) then
   insid=0
 else
   insid=1
 endif
 aux = qdr-floor(qdr)
 i = nint(1+ n*aux)
! write(*,7)'i,qdotr,aux=',i,qdr,aux
 if (i.eq.n+1) i=1

7 format(a,i5,9(1x,g12.5))
 end subroutine comp1
!============================================================
 subroutine comp_c(q,rr,n,i,insid)
! inside would be =1 if -0.5<q<0.5
 use constants, only : r15
 use geometry
 implicit none
 integer, intent(in):: n
 integer, intent(out):: i,insid
 real(r15),intent(in) :: q(3)
 type(vector), intent(in) :: rr
 real(r15) aux,qdr

 qdr = (q.dot.rr)
! include from -0.00005 to 0.99995 included (include 0.0 and exclude 1.0)
 if (qdr.lt.-0.5d0 .or. qdr.ge.0.5d0) then
   insid=0
 else
   insid=1
 endif

! in order to get i, bring qdr between 0 and 1 by adding 0.5
 aux =  0.5d0+qdr-floor(qdr+0.5d0)
 i = nint(1+ n*aux)
! write(*,7)'i,qdotr,aux=',i,qdr,aux
 if (i.eq.n+1) i=1

7 format(a,i5,9(1x,g12.5))
 end subroutine comp_c
!============================================================
 subroutine get_components_g(q,n,i,j,k,inside)  ! g's  are primitive g0i
! for a given q-vector, it finds its integer components assuming it was
! created as: q=(i-1)/n1*g1 + (j-1)/n2*g2 + (k-1)/n3*g3; i=1,n1 j=1,n2 k=1,n3
! as the ones created in make_kp_reg with zero shift
! if the integer variable inside=1 then q is inside the primcell
! if q is outside the prim_cell, then inside=0 but produces the i,j,k of its
! image inside
 use lattice
 use geometry
 use constants, only : pi,r15
 implicit none
 real(r15), intent(in):: q(3)
 integer, intent(in):: n(3)
 integer, intent(out):: i,j,k,inside
 real(r15) i2pi

 i2pi=1d0/(2*pi)
! if q.dot.r is larger than 1, we want to have its fractional part
 inside = 1

 call comp1(q,i2pi*r01,n(1),i,inside)
 call comp1(q,i2pi*r02,n(2),j,inside)
 call comp1(q,i2pi*r03,n(3),k,inside)

 end subroutine get_components_g
!============================================================
 subroutine get_components_g_centered(q,n,i,j,k,inside) !gg1,gg2,gg3,inside)
! for a given q-vector, it finds its integer components assuming it was
! created as: q=(i-1-n1/2)/n1*g1 + (j-1-n2/2)/n2*g2 + (k-1-n3/2)/n3*g3;
! i=1,n1 j=1,n2 k=1,n3
! as the ones created in make_kp_reg with zero shift
! if the integer variable inside=1 then q is inside the primcell, if inside=0 it's outside
! it works even if there's a shift less than 0.5, and if q is oustide the prim_cell
 use geometry
 use lattice
 use constants, only : pi,r15
 implicit none
 real(r15), intent(in):: q(3)
! type(vector),intent(in):: gg1,gg2,gg3
 integer, intent(in):: n(3)
 integer, intent(out):: i,j,k,inside
! real(r15) w(3)
 real(r15) i2pi

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
 use constants, only : r15
 implicit none
 integer nkp,i,j,nbmax,l,sh,nshl,sh0
 integer nb(nkp,nbmax),nbk(nkp)
 real(r15) kp(3,nkp),q2(3),qw(3),gg(3,nshl),dismin,kcutoff,kdist

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
 use constants, only : r15
!! use io2
 implicit none
 real(r15), intent(in) :: kp(3)
 real(r15), intent(out) :: q(3)
 real(r15) a1,a2,a3

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
 use constants, only : r15
 implicit none
 integer nj,j,n
 real(r15) v(3),grid(3,n),tolerance,y1(3),y2(3),y3(3),i1,i2,i3

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
 subroutine find_in_array(n1,v,n,grid,nj,tolerance)
! checks whether the vector v belongs to the array "grid" within the tolerance;
! nj is its index in the grid; if nj=0, the vector did not belong to the grid
 use geometry
 use constants, only : r15
 implicit none
 integer, intent(in) :: n1,n
 integer, intent(out) :: nj
 real(r15), intent(in) :: v(n1),grid(n1,n),tolerance
 integer j

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
 use constants, only : r15
      implicit none
      integer, intent(in) :: n
      integer imax,k,j,i,ii
      real(r15), intent(inout) :: a(n,n)
      real(r15), intent(out) :: b(n,n)
      real(r15) amax,atmp,btmp,amult,div


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
 use constants, only : r15
  implicit none
  integer n
  real(r15) x(:),mean,sd,sd2

  n=size(x)
  mean = sum(x)/n
  sd2  = sum(x*x)/n
  sd   = sqrt(sd2 - mean*mean)

  end subroutine mean_sd
!=====================================================
  subroutine histogram(m,x,mesh)
! calculates the histogram of the data set x (needs not be sorted)
! and writes the distribution function in the file 'histo.dat'
 use constants, only : r15
  implicit none
  integer, intent(in):: m,mesh
  real(r15), intent(in):: x(m)
  integer i,j,cnt,unit
  real(r15) xmin,xmax,dx,e(mesh)

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
!===========================================================
 recursive function determinant(mat) result(det)
 use constants, only : r15
 implicit none
 integer n,i,j,k
! real(r15) det,mat(n,n),cofact(n-1,n-1)
 real(r15), dimension(:,:), intent(in)  :: mat
 real(r15), allocatable :: cofact(:,:)
 real(r15) det

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
     FUNCTION my_erfc(x) result(S15ADF)
     use constants, only : r15
!     MARK 5A REVISED - NAG COPYRIGHT 1976
!     MARK 5C REVISED
!     MARK 11.5(F77) REVISED. (SEPT 1985.)
!     COMPLEMENT OF ERROR FUNCTION ERFC(X)
!
!     .. Scalar Arguments ..
      real(r15), intent(in) :: X
!     .. Local Scalars ..
      real(r15)   T, XHI, XLO, Y, s15adf
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

      END function my_erfc !s15adf
!===========================================================
     DOUBLE PRECISION FUNCTION my_old_erfc(x) result(S15ADF)
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

      END function my_old_erfc !s15adf
!==================================================
      SUBROUTINE choldc(a,n,p)
! Choleski decomposition of a=LLT; p contains the diagonals of L
! upper triangular part of a is used; Lower part of L is stored in lower part of a
 use constants, only : r15
      implicit none
      INTEGER n
      real(r15) a(n,n),p(n)
      INTEGER i,j,k
      real(r15) sum
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
 subroutine make_grid2(x01,x02,x03,x1,x2,x3,n,grid)
!! generates a grid of vectors "grid" linear combinations of x0i but
!! contained in the Wigner-Seitz cell of x1,x2,x3, including the boundaries
 use geometry
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
! real(r15), intent(in ):: x01(3),x02(3),x03(3),x1(3),x2(3),x3(3)
 type(vector), intent(in ):: x01,x02,x03,x1,x2,x3
 real(r15), intent(out):: grid(3,n)
 real(r15) v(3),a1,a2,a3
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
 end subroutine make_grid2
!============================================================
 subroutine get_upper_bounds(x1,x2,x3,rcut,maxatm,nsmx)
!! finds maxp, the number of grid points in the sphere of radius rcut and m, the needed mesh
!! outputs are nsmx and maxatm: upper bounds in 1D and in 3D (used for 3 loops and allocation)
 use geometry
 use atoms_force_constants, only : natom_prim_cell
 use ios, only : ulog
 use constants, only : r15
 implicit none
 type(vector),intent(in) :: x1,x2,x3
 integer,intent(out) :: nsmx,maxatm ! upper boundary in 1D and in 3D (used for 3 loops and allocation)
 real(r15),intent(in) :: rcut
 integer m1,m2
 real(r15) om0,ls,lmax

 call calculate_volume(x1,x2,x3,om0)
 if(om0.myeq.0d0) then
    write(*,*)'CALVOL: Volume is zero! check your basis vectors!'
    write(*,'(a,9(1x,g12.5))')'r1,r2,r3=',x1,x2,x3
    stop
 endif
 ls=om0**0.333333 ! typical unitcell size
 maxatm=(nint(12.56/3d0*(1.3*(rcut+ls))**3/om0)+1)*natom_prim_cell  ! upperbound to the number of atoms in the sphere

 lmax= max(length(x3),length(x1),length(x2))
! find shortest translation vector length
 ls=om0/lmax/lmax
 ls=min(ls,length(x3),length(x1),length(x2))
 write(ulog,3)'GET_UPPER_BOUNDS: lower bound to translation vectors length=',ls

 m1=nint(rcut/ls)+2
 m2=(nint((3/12.56*(maxatm/natom_prim_cell))**0.33333)+1)
 nsmx=max(m1,m2)
! maxp=(2*m+1)**3

 write(ulog,4)'GET_UPPER_BOUNDS:rcut,(rcut+ls)/2,nsmx,maxatom=',rcut,m1,nsmx,maxatm

3 format(a,9(1x,g14.7))
4 format(a,f9.4,1x,3i8)

 end subroutine get_upper_bounds
!============================================================
 subroutine apply_metric(x1,x2,x3,eps,t1,t2,t3)
!! generates vectors ti such that ti=sqrt(eps)xi
 use geometry
 use ios
 use constants, only : r15
 implicit none
 type(vector), intent(in ):: x1,x2,x3
! type(vector), intent(out):: t1,t2,t3
 real(r15), intent(out) :: t1(3),t2(3),t3(3)
 real(r15) eps(3,3),eps2(3,3),sqeps(3,3),a1(3)
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
 subroutine make_grid_weights_WS(x01,x02,x03,x1,x2,x3,matr,ngrd,grd,weig,space,s0x,sx) 
!! generates a grid "grd" from primitive translation vectors x0i
!! but contained in the Wigner-Seitz cell of x1,x2,x3, with corresponding weights
!! used for fourier interpolation 
!! also outputs the 26 shortest translation vectors s0x and sx used for 
!! defining the WS cell of x0i and xi 
 use geometry
 use ios, only : ulog
 use constants, only : pi,r15
 use lattice
 use params, only : tolerance
 implicit none
 type(vector), intent(in):: x01,x02,x03,x1,x2,x3
 real(r15), intent(in):: sx(3,26),s0x(3,26)
 integer, intent(inout):: ngrd
 real(r15), intent(inout):: weig(ngrd),grd(3,ngrd)
 real(r15), intent(in):: matr(3,3)
 character(len=1), intent(in) :: space
 integer i1,i2,i3,m,cnt,ns,nboundary,nbtot,ns2
 logical insid
 real(r15) v(3),omx,omx0,tol,y1,y2,y3 ,gmax ,vfold(3)
 integer, allocatable :: save_boundary(:)
 real(r15), allocatable :: aux(:,:)

! It is assumed x0i and xi are the shortest basis vectors

! y's are used to find the reduced coordinates on the basis xi
! call make_reciprocal_lattice_v(x1,x2,x3,y1,y2,y3)

! find the number of grid points inside the supercell
 call calculate_volume(x1,x2,x3,omx)
 call calculate_volume(x01,x02,x03,omx0)
 if((omx .myeq. 0d0) .or. (omx0 .myeq. 0d0)) then
    write(*,*)'CALVOL: Volume is zero! check your basis vectors!'
    write(*,'(a,9(1x,g12.5))')'r01,r02,r03=',x01,x02,x03
    write(*,'(a,9(1x,g12.5))')'r1,r2,r3=',x1,x2,x3
    stop
 endif
 tol = tolerance*omx0**0.33333

 ns=nint(omx/omx0)
 write(ulog,6)' grid size with Xs ,om0,om=',ns,ngrd,omx0,omx


! this is an upperbound to the points inside and on the WS boundaries
 ns=ns+nint(14*ns**0.67+36*ns**0.3333+24 )

 m=maxval(abs(n_sc))+1
 write(ulog,*)' upper index to loop over to create grid =',m
 ns2 = (2*m+1)**3

 write(ulog,*)'First and second ns estimates =',ns,ns2
 allocate(aux(3,ns))
 write(ulog,*)' allocated grid size, including boundary points=',ns !(2*m+1)**3
! aux contains the grid vectors inside or on the boundary

! create the grid within the WS supercell xi
 cnt=0
 shelloop: do m=0,maxval(abs(n_sc))+3
 do i1=-m,m
 do i2=-m,m
 do i3=-m,m
    if(iabs(i1).ne.m.and.iabs(i2).ne.m.and.iabs(i3).ne.m)cycle

! generate a grid of vectors from primitive translations x0i
    v= v2a(dble(i1)*x01 + dble(i2)*x02 + dble(i3)*x03)

! see if v is in the WS cell; throw away v outside of WS cell
!   write(*,4)'* trying vector ',cnt,i1,i2,i3,v,matmul(matr,v)
    call check_inside_ws(v,sx,insid,nboundary)
!   if(insid) then
!      write(*,4)' inside vector ',cnt,i1,i2,i3
!   else
!      cycle
!   endif
! skip if it is outside
    if (.not. insid) cycle 

    vfold= fold_ws(v,sx,space)
    write(234,2)v,vfold

! take only those points inside or on the boundary of WS (corresponding to in=1)
    cnt=cnt+1
    write(*,4)'* New vector ',cnt,i1,i2,i3,v,matmul(matr,v)
    aux(:,cnt)=v
    if(cnt.eq.ns) then
       write(*,*)'Max size of aux reached ',cnt,ns,' stopping the loop!'
       write(ulog,*)'Max size of aux reached ',cnt,ns,' stopping the loop!'
       stop 
    endif

 enddo
 enddo
 enddo
 enddo shelloop

 ngrd=cnt

! if(allocated(grd)) deallocate(grd)
! if(allocated(weig)) deallocate(weig)
! if(allocated(save_boundary)) deallocate(save_boundary)
! allocate(grd(3,ngrd),save_boundary(ngrd),weig(ngrd))
 allocate(save_boundary(ngrd))
 save_boundary=0
 grd(:,1:ngrd)=aux(:,1:ngrd) !; weig = weig(1:ngrd)
 deallocate(aux)

! find the weights of each point on the grid (and its boundary)
 weig=1
 nbtot=0
 write(ulog,*)' grid vectors, number on boundary, coordinates  =========='
 do cnt=1,ngrd
! which grd on the boundary? nboundary=on how many facets;=0 means not on boundary
    call is_on_boundary(grd(:,cnt), sx,nboundary,tolerance) 
    if(nboundary.ne.0) then
       save_boundary(cnt)=nboundary
       nbtot=nbtot+1
    endif
    write(ulog,6)'i, nboundaries, grid(i), grid_reduced ',cnt,nboundary, &
&                 grd(:,cnt),matmul(matr,grd(:,cnt)),length(grd(:,cnt))
 enddo
 write(ulog,*)'Of ',ngrd,' grid points ',nbtot,' are on the boundary'

 call find_ws_weights(ngrd,grd,save_boundary,sx,weig)

 if(abs(sum(weig)-1).gt.tolerance) then
    write(*,*)'WARNING: weights not normalized to 1 ',sum(weig(1:ngrd))
    write(ulog,*)'WARNING: weights not normalized to 1 ',sum(weig(1:ngrd)) 
 endif

 deallocate(save_boundary)

 write(ulog,*) ' EXITING make_grid_weights_WS'

2 format(9(1x,f10.4))
3 format(i5,1x,f8.5,3x,3(1x,g10.3),3x,3(1x,f8.4))
4 format(a,i5,2x,3i2,2x,9(1x,f13.5))
5 format(3(1x,f10.4),i3,1x,g11.4)
6 format(a,2i5,99(1x,f10.4))
7 format(a,9(1x,f10.4))
8 format(a,7(1x,f10.4),3i5)

 end subroutine make_grid_weights_WS
!============================================================
 subroutine make_grid_weights_WS_old(x01,x02,x03,x1,x2,x3,ngrd,space,s0x,sx) !,weig)
!! generates a grid "rgrid" or "ggrid" (depending on space) from primitive translation vectors
!! x0i but contained in the Wigner-Seitz cell of x1,x2,x3, with corresponding weights
!! works for fourier interpolation in the reciprocal space because get_stars works for q-points
!! space is r (for real space) or g (for reciprocal space)
!! the 26 shortest translation vectors sx for defining the WS cell of xi is also output
 use geometry
 use ios, only : ulog
 use constants, only : pi,r15
 use lattice
 use params, only : tolerance
 use fourier   ! grid points and their weights defined in Fourier
 implicit none
 type(vector), intent(inout):: x01,x02,x03,x1,x2,x3
 character(len=1), intent(in) :: space
 real(r15), intent(out):: sx(3,26),s0x(3,26)
 integer, intent(out):: ngrd
 integer i1,i2,i3,m,cnt,ns,nboundary,nbtot
 logical insid
 type(vector) b01,b02,b03,b1,b2,b3
 real(r15) v(3),omx,omx0,tol,matr(3,3),y1,y2,y3 ,gmax !,vfold(3)
 integer, allocatable :: save_boundary(:)
 real(r15), allocatable :: aux(:,:),grd(:,:)
 real(r15), allocatable :: weig(:)

! find the shortest basis bi among xi and the corresponding 26 shortest vectors
 call get_26shortest_shell(x01,x02,x03,s0x,b01,b02,b03)! s0x for WS of primcell
 call get_26shortest_shell(x1 ,x2 ,x3 ,sx ,b1 ,b2 ,b3) ! sx for WS of supercell
 x01=b01; x02=b02; x03=b03; x1=b1;  x2=b2;  x3=b3; 

 if(space.eq.'r' .or. space.eq.'R') then
    matr=cart_to_prim
 elseif(space.eq.'g' .or. space.eq.'G') then
    matr=transpose(prim_to_cart)/(2*pi)
 endif

 tol = tolerance*length(b1)

! y's are used to find the reduced coordinates on the basis xi
! call make_reciprocal_lattice_v(x1,x2,x3,y1,y2,y3)

! find the number of grid points inside the supercell
 call calculate_volume(b1,b2,b3,omx)
 call calculate_volume(b01,b02,b03,omx0)
 ns=nint(omx/omx0)
 write(ulog,6)' grid size with Bs ,om0,om=',ns,ngrd,omx0,omx

 call calculate_volume(x1,x2,x3,omx)
 call calculate_volume(x01,x02,x03,omx0)
 ns=nint(omx/omx0)
 write(ulog,6)' grid size with Xs ,om0,om=',ns,ngrd,omx0,omx

! this is an upperbound to the points inside and on the WS boundaries
 ns=ns+nint(6*ns**0.67+12*ns**0.3333+8 )

 m=maxval(abs(n_sc))+1
 write(ulog,*)' upper index to loop over to create grid =',m
 ns = (2*m+1)**3
 allocate(aux(3,ns))
 write(ulog,*)' allocated grid size, including boundary points=',(2*m+1)**3
! aux contains the grid vectors inside or on the boundary

! create the grid within the WS supercell xi
 cnt=0
 shelloop: do m=0,maxval(abs(n_sc))+3
 do i1=-m,m
 do i2=-m,m
 do i3=-m,m
    if(iabs(i1).ne.m.and.iabs(i2).ne.m.and.iabs(i3).ne.m)cycle
! generate a grid of vectors from primitive translations x0i
    v= v2a(dble(i1)*x01 + dble(i2)*x02 + dble(i3)*x03)

!   vfold= fold_ws(v,sx)
!   write(234,2)v,vfold

! see if v is in the WS cell; throw away v outside of WS cell
!   write(*,4)'* trying vector ',cnt,i1,i2,i3,v,matmul(matr,v)
    call check_inside_ws(v,sx,insid,nboundary)
!   if(insid) then
!      write(*,4)' inside vector ',cnt,i1,i2,i3
!   else
!      cycle
!   endif
! skip if it is outside
    if (.not. insid) cycle 

! take only those points inside or on the boundary of WS (corresponding to in=1)
    cnt=cnt+1
    write(*,4)'* New vector ',cnt,i1,i2,i3,v,matmul(matr,v)
    aux(:,cnt)=v
    if(cnt.eq.ns) then
       write(*,*)'Max size of aux reached ',cnt,ns,' stopping the loop!'
       write(ulog,*)'Max size of aux reached ',cnt,ns,' stopping the loop!'
       stop !exit i1loop
    endif

 enddo
 enddo
 enddo
 enddo shelloop

 ngrd=cnt

 if(allocated(grd)) deallocate(grd)
 if(allocated(weig)) deallocate(weig)
 if(allocated(save_boundary)) deallocate(save_boundary)
 allocate(grd(3,ngrd),save_boundary(ngrd),weig(ngrd))
 save_boundary=0
 grd(:,1:ngrd)=aux(:,1:ngrd)
 deallocate(aux)
! find the weights of each point on the grid (and its boundary)
 weig=1
 nbtot=0
 write(ulog,*)' grid vectors, number on boundary, coordinates  =========='
 do cnt=1,ngrd
! identify vectors on the boundary
    call is_on_boundary(grd(:,cnt), sx,nboundary,tolerance)  ! nboundary=on how many facets;=0 means not on boundary
    if(nboundary.ne.0) then
       save_boundary(cnt)=nboundary
       nbtot=nbtot+1
    endif
    write(ulog,6)'i,on how many boundaries, grid(i), grid_reduced ',cnt,nboundary,grd(:,cnt),matmul(matr,grd(:,cnt))
 enddo
 write(ulog,*)'Of ',ngrd,' grid points ',nbtot,' are on the boundary'

 call find_ws_weights(ngrd,grd,save_boundary,sx,weig)

 if(space.eq.'r' .or. space.eq.'R') then
    if (allocated(rgrid)) deallocate(rgrid)
    if (allocated(rws_weights)) deallocate(rws_weights)
    allocate(rgrid(3,ngrd), rws_weights(ngrd) ) !,i1r(ngrd),i2r(ngrd),i3r(ngrd))
    nrgrid=ngrd
    rgrid = grd
!    i1r=i1p(1:ngrd); i2r=i2p(1:ngrd); i3r=i3p(1:ngrd)
! normalize weights by volume ratio (should normalize to 1 for real space vectors)
    weig = weig*omx0/omx
    rws_weights = weig
    write(ulog,*)'sum of the weights=',sum(rws_weights) 
    if(abs(sum(rws_weights)-1).gt.tolerance) then
       write(*,*)'WARNING: Rweights not normalized to 1 ',sum(rws_weights)
       write(ulog,*)'WARNING: Rweights not normalized to 1 ',sum(rws_weights) 
    endif
    open(98,file='rgrid_raw.xyz')
    open(99,file='rgridWS.xyz')
!    call show_ws_boundary(v2a(x1),v2a(x2),v2a(x3),sx,15,'WSR_boundary.xyz',lgridmax) 

    do cnt=1,ngrd
       write(98,8)"Si ",grd(:,cnt),weig(cnt),matmul(matr,grd(:,cnt)),cnt,save_boundary(cnt)
       write(99,7)"Si ",grd(:,cnt),matmul(matr,grd(:,cnt))
    enddo
    close(98)
    close(99)
    write(ulog,*)'lgridmax=',lgridmax

 elseif(space.eq.'g' .or. space.eq.'G') then
    if (allocated(ggrid)) deallocate(ggrid)
    if (allocated(gws_weights)) deallocate(gws_weights)
    allocate(ggrid(3,ngrd), gws_weights(ngrd))
    nggrid=ngrd
    ggrid = grd
! normalize weights (should normalize to vol/vol0)
    gws_weights = weig ! / (volume_r/volume_r0) ! introduce 1/N since used for Fourier transforms
    write(ulog,*)'sum of the weights=',sum(weig)
    write(ulog,*)'Gweights normalization is ',sum(weig) 
    open(98,file='ggrid_raw.xyz')
    open(99,file='ggridWS.xyz')
!    call show_ws_boundary(v2a(x1),v2a(x2),v2a(x3),sx,15,'WSG_boundary.xyz',gmax) 

    write(98,8)"# name ,grid(cnt),weig(cnt),grid_red(cnt),cnt,save_boundary(cnt)"
    write(99,*)ngrd+26+8
    write(99,*)"# name, cartesian grid, reduced grid "
    do cnt=1,ngrd
       write(98,8)"Si ",grd(:,cnt),weig(cnt),matmul(matr,grd(:,cnt)),cnt,save_boundary(cnt)
       write(99,7)"Si ",grd(:,cnt),matmul(matr,grd(:,cnt))
    enddo
    do cnt=1,26
       write(99,7)"Bi ",sx(:,cnt), matmul(matr,sx(:,cnt))
    enddo
    write(99,*)'Ge   0 0 0 '
    write(99,7)'Ge ',x01
    write(99,7)'Ge ',x02
    write(99,7)'Ge ',x03
    write(99,7)'Ge ',x01+x02
    write(99,7)'Ge ',x03+x01
    write(99,7)'Ge ',x03+x02
    write(99,7)'Ge ',x03+x01+x02
    close(98)
    close(99)

 endif

 deallocate(grd,save_boundary,weig)

 write(ulog,*) ' EXITING make_grid_weights_WS'

2 format(9(1x,f10.4))
3 format(i5,1x,f8.5,3x,3(1x,g10.3),3x,3(1x,f8.4))
4 format(a,i5,2x,3i2,2x,9(1x,f13.5))
5 format(3(1x,f10.4),i3,1x,g11.4)
6 format(a,2i5,99(1x,f10.4))
7 format(a,9(1x,f10.4))
8 format(a,7(1x,f10.4),3i5)

 end subroutine make_grid_weights_WS_old
!============================================================
! subroutine make_periodic_WS(array,n,grid)
!!! an array defined on the grid inside the WS cell, should be periodic. IF not, this subroutine makes it!
! integer, intent(in) :: n
! real(r15), intent(inout) :: array(n),grid(3,n)
!
!! find boundary points that are equivalent, and take their average
!
! end subroutine make_periodic_WS
!============================================================
 subroutine make_subrgrid_symmetric(n,grd,sx,nsub)
!! makes a subgrid of the R-translations grd that has the full symmetry of the crystal
 use fourier
 use geometry
 use params, only : tolerance
 use lattice, only : cart_to_prim
 use ios, only : ulog,write_out
 use atoms_force_constants , only : lattpgcount, op_matrix
 implicit none
 integer, intent(in) :: n
 integer, intent(out) :: nsub
 real(r15), intent(in) :: grd(3,n),sx(3,26)
 real(r15), allocatable :: subg(:,:)
 integer, allocatable :: save_boundary(:)
 real(r15) v(3)
 integer i,j,k,g,nboundary,nbtot
 logical belongs,is_sub

 allocate(subg(3,n))
 nsub=0
 gridloop: do g=1,n
   !  write(ulog,6)'rgrid vectors ',n,g, grd(:,g)
! make sure all arms of grd are in the grd set
      is_sub=.false.
      grouploop: do i=1,lattpgcount
! apply symmetry operation to grd
!call write_out(ulog,'op_matrix ',op_matrix(:,:,i))
         v=matmul(op_matrix(:,:,i),grd(:,g))
         belongs=.false.
         kloop: do k=1,n
            if(length(v-grd(:,k)).lt.tolerance) then  ! arm belongs to the set, OK!
               belongs=.true.
 !    write(ulog,6)'g,k,v=',g,k,v
               exit kloop
            endif
         enddo kloop
         if (belongs) then ! v=S(grd(g)) belongs to the grid; try other point group operations 
            cycle grouploop
         else
            cycle gridloop
         endif
      enddo grouploop
      is_sub=.true. ! all point group symmetries mapped grd to another grd
      nsub=nsub+1
      subg(:,nsub)=grd(:,g)
 enddo gridloop
  
! now find the weights
 allocate(subgrid(3,nsub),save_boundary(nsub),subgrid_weights(nsub))
 subgrid = subg(:,1:nsub)
 deallocate(subg)
 save_boundary=0
 subgrid_weights=1
 nbtot=0
 do g=1,nsub
! identify vectors on the boundary
    call is_on_boundary(subgrid(:,g),sx,nboundary,tolerance)  ! nboundary=on how many facets;=0 means not on boundary
    if(nboundary.ne.0) then
       save_boundary(g)=nboundary
       nbtot=nbtot+1
    endif
    write(ulog,6)'i,on how many boundaries, subgrid(i), subgrid_reduced ',g,nboundary,subgrid(:,g),matmul(cart_to_prim,subgrid(:,g))
 enddo
 write(ulog,*)'Of ',nsub,' subgrid points ',nbtot,' are on the boundary'
 write(ulog,6)'CHECK: sum of sub_ saveboundary=',nbtot,sum(save_boundary)
! now find the weight of subg set
 call find_ws_weights(nsub,subgrid,save_boundary,sx,subgrid_weights)
 subgrid_weights=subgrid_weights/sum(subgrid_weights(1:nsub))

    open(98,file='subgrid_raw.xyz'); 
    open(99,file='subgridWS.xyz')
    do g=1,nsub
       write(98,8)"Si ",subgrid(:,g),subgrid_weights(g),matmul(cart_to_prim,subgrid(:,g)),g,save_boundary(g)
       write(99,7)"Si ",subgrid(:,g),matmul(cart_to_prim,subgrid(:,g))
    enddo
    close(98); close(99)

 deallocate(save_boundary)

6 format(a,2i5,99(1x,f10.5))
7 format(a,9(1x,f12.5))
8 format(a,7(1x,f12.5),3i5)
 
 end subroutine make_subrgrid_symmetric
!============================================================
 subroutine find_cutoff_from_largest_SC (volmax,rmax,r1,r2,r3,rshells)
!! scans over all available supercell POSCARs and picks the one with
!! largest volume, and check its consistency with the primitive cell
 use params, only : fdfiles
 use ios   , only : ulog
 use lattice, only : volume_r,cart_to_prim
 use constants, only : r15
 use atoms_force_constants , only : natom_super_cell,atom_sc
 use geometry, only : vector ,v2a,a2v
 implicit none
! integer, intent(out) :: imax  ! index of largest supercell
 real(r15), intent(out) :: r1(3),r2(3),r3(3),volmax,rmax,rshells(3,26)
 real(r15) volum,rm1(3),rm2(3),rm3(3),lmax,sx(3,26)
 type(vector) v1,v2,v3
 integer i,imx
 character xt*1,poscar*7

     volmax=0;rmax=0             !  pick the SC with largest length rmax
     do i=1,fdfiles
        write(xt,'(i1)')i
        poscar='POSCAR'//xt
! read the atomic coordinates of atoms in supercell (at equilibrium) from POSCAR
        call read_supercell_vectors(poscar,volum,r1,r2,r3) 
        call cube(r1,r2,r3,sx)
        call show_ws_boundary(r1,r2,r3,sx,13,'transl.xyz',lmax) 
        write(ulog,*) 'supercell ',i,' read, rmax,volume=',lmax,volum
        if(lmax.gt.rmax) then
          imx=i
          rmax=lmax; rm1=r1; rm2=r2; rm3=r3
        endif
     end do
!    call get_26shortest_shell(rm1,rm2,rm3,rshells,r1,r2,r3)
     call get_26shortest_shell(a2v(rm1),a2v(rm2),a2v(rm3),rshells,v1,v2,v3)
     r1=v2a(v1);r2=v2a(v2);r3=v2a(v3)
     write(ulog,*)'FIND_CUTOFF_FROM LARGEST_SC: rmax=',rmax

 end subroutine find_cutoff_from_largest_SC
!============================================================
 subroutine cube(r1,r2,r3,sx)
 use constants, only : r15
 implicit none
 real(r15), intent(in) :: r1(3),r2(3),r3(3)
 real(r15), intent(out) :: sx(3,26)
 integer i,j,k,cnt

 cnt=0
 do i=-1,1
 do j=-1,1
 do k=-1,1
    if (i*i+j*j+k*k .eq.0) cycle
    cnt=cnt+1
    sx(:,cnt)=i*r1+j*r2+k*r3
 enddo
 enddo
 enddo

 end subroutine cube
!============================================================
 subroutine find_WS_largest_SC(imax,volmax)
!! scans over all available supercell POSCARs and picks the one with
!! largest volume, and check its consistency with the primitive cell
 use params, only : fdfiles
 use ios   , only : ulog
 use lattice, only : volume_r,cart_to_prim
 use constants, only : r15
 use atoms_force_constants , only : natom_super_cell,atom_sc
 use geometry, only : v2a
 implicit none
 integer, intent(out) :: imax  ! index of largest supercell
 real(r15), intent(out) :: volmax
 integer i
 character xt*1,poscar*7

     volmax=0             !  pick the SC with largest volume
     do i=1,fdfiles

        write(xt,'(i1)')i
        poscar='POSCAR'//xt
! read the atomic coordinates of atoms in supercell (at equilibrium) from POSCAR
        call read_supercell(poscar)
        write(ulog,*) 'supercell ',i,' read, volume=',volume_r

! find the largest volume
        if(volume_r .gt. volmax) then
           volmax=volume_r
           imax=i
        endif

     enddo

     write(ulog,*) 'largest volume is found for ',i,' th supercell; volmax=',volmax

! work with this POSCAR
     write(xt,'(i1)')imax
     poscar='POSCAR'//xt
     call read_supercell(poscar)

! open(745,file='supercell.xyz')
! write(745,*) natom_super_cell
! write(745,*)' name, x , y , z , i,tau,n(3), jatompos, reduced coordinates '
! do i=1,natom_super_cell
!! name has not been assigned yet
!   write(745,7)'Si ',atom_sc(i)%equilibrium_pos,i,atom_sc(i)%cell%tau,atom_sc(i)%cell%n  &
!&      ,atom_sc(i)%cell%atomposindx,matmul(cart_to_prim,v2a(atom_sc(i)%equilibrium_pos))
! enddo
! close(745)

7 format(a,3(1x,f12.5),i5,i3,'(',3i2,')',i5,3(1x,f6.3))
     call check_input_poscar_consistency_new

 end subroutine find_WS_largest_SC
!========================================
 subroutine update_nshells(ng,grd,lcut) 
!! resets nshells(2,:) according to the rgrid vectors at the boundary of the WS of supercell
!! if input range is larger than cutoff(SC) or if fc2range=0, sets nshells(2,:) corresponding to that cutoff 
 use geometry
 use ios
 use lattice, only : volume_r,volume_r0 
 use atoms_force_constants
 use params, only : nshells !,lmax
 use constants, only : r15
 implicit none
 integer, intent(in):: ng
 real(r15), intent(in):: grd(3,ng),lcut
 integer i,shel_count,i0,shelmax(natom_prim_cell),shlmx

!! find lmax, length of the shortest vector on the boundary of rgrid
! lmax=0
! do i=1,ng
!    if( length(grd(:,i)) .gt. lmax  ) lmax = length(grd(:,i))
! enddo
 write(ulog,4)' UPDATING NSHELLS up to a cutoff length of ',lcut

! find shelmax, the largest shell
 shelmax=0
 do i0=1,natom_prim_cell
    shlmx=atom0(i0)%nshells ! is supposedly already too large (> lcut)
    write(*,*)'i0,shlmx,lcut=',i0,shlmx,lcut
    write(ulog,*)'i0,shlmx,lcut=',i0,shlmx,lcut
    shloop: do shel_count=0,shlmx
!      write(ulog,*)'i0, shell_count,radius=',i0,shel_count,atom0(i0)%shells(shel_count)%radius
       if( atom0(i0)%shells(shel_count)%radius .gt. lcut*1.001 ) then
          shelmax(i0)=shel_count
          exit shloop
       endif
! did not reach lcut
       shelmax(i0)=shlmx
    enddo shloop
    write(ulog,*)'i0,shelmax(i0)=',i0,shelmax(i0)
 enddo

! include all the shells up to distance lcut
 do i0=1,natom_prim_cell
    shell_loop: do shel_count = 0 , shelmax(i0) 
! we update nshells if lcut is reached, otherwise if lcut > rcut no update
       if(lcut .myeq. atom0(i0)%shells(shel_count)%radius) then
          nshells(2,i0) = shel_count
          write(ulog,*)'UPDATE_NSHELLS: for atom ',i0,' the nshell is updated to ',shel_count 
          exit shell_loop
       else ! take the farthest shell available within rcut
          nshells(2,i0)=shelmax(i0)
       endif
    enddo shell_loop
    atom0(i0)%nshells=shel_count
 enddo

 write(ulog,5)' ****************************************'
 write(ulog,5)' UPDATED SHELLS of rank 2 ',nshells(2,1:natom_prim_cell)
 write(ulog,5)' ****************************************'

4 format(a,9(1x,f13.5))
5 format(a,20(1x,i3))

 end subroutine update_nshells
!==========================================================
 subroutine find_ws_weights(ngrd,grd,save_boundary,transl,weight)
!! finds the # of times a point within FBZ is on its boundaries
 use geometry
 use constants, only : r15
 implicit none
 integer, intent(in) :: ngrd
 real(r15), dimension(3,ngrd), intent(in) :: grd
 real(r15), dimension(3,26)  , intent(in) :: transl
 integer, dimension(ngrd)  , intent(in) :: save_boundary
 real(r15), dimension(ngrd) , intent(out) :: weight
 real(r15) x(3)
 integer i,j,t !ngrd,ntransl,

! ngrd=size(grd(1,:))
 weight=1
! if grid point is on the boundary, find to how many boundary points it is mapped by a transl vector
 do i=1,ngrd
    if (save_boundary(i).ne.0) then  ! i is on the boundary
       do t=1,26   ! count on how many facets or boundaries
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
    write(*,*)i,' th grid point has weight;1/weight=',weight(i),1d0/weight(i)
 enddo
 write(*,*)' sum of the weights=om/om0=',sum(weight)

 end subroutine find_ws_weights
!============================================================
 subroutine belongs_to_grid(v,n,grid,nj,tolerance)
! checks whether the vector v belongs to the array "grid" within the 1d-4;
! nj is its index in the grid; if nj=0, the vector did not belong to the grid
 use geometry
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 integer, intent(out) :: nj
 real(r15), intent(in) :: v(3),grid(3,n),tolerance
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
 use constants, only : r15
 implicit none
 integer, intent(out) :: nq(3)
 real(r15), intent(in)  :: q(3,3)
 real(r15) c(3),d(3)
 logical sing

 call qrdcmp(q,3,3,c,d,sing)

 nq=nint(d)
! nq(1)=nint(q(1,1))
! nq(2)=nint(q(2,2))
! nq(3)=nint(q(3,3))
 end subroutine myreduce
!========================================
      SUBROUTINE qrdcmp(a,n,np,c,d,sing)
 use constants, only : r15
      implicit none
      INTEGER n,np
      real(r15) a(np,np),c(n),d(n)
      LOGICAL sing
      INTEGER i,j,k
      real(r15) scale,sigma,sum,tau

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
 use constants, only : r15
      use atoms_force_constants
      implicit none
      integer i,j,k,narms,kvecop(48)
      real(r15) kvec(3),kvecstar(3,48),primlatt(3,3),v(3),  &
     &     v2(3),v3(3),kvecstarp(3,48),sum
      narms=0
!      print*,'lattpgcount=',lattpgcount

      iloop: do i=1,lattpgcount
! apply symmetry operation to k to get v=kstar=opk(i)*k
!        call xvmlt(op_kmatrix(1,1,i),kvec,v,3,3,3)
        v=matmul(iop_kmatrix(:,:,i),kvec)
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
 subroutine fold_in_WS(kp,gshel,foldedk)
!! for the BZ defined by (gshel) this subroutine generates foldedk obtained
!! from folding of kp into the WS Wigner-Seitz cell of FBZ,
 use lattice , only : cart2red_g ! prim_to_cart  !r01, r02, r03,
 use geometry
 use kpoints  !, only : nkc,kpc
 use constants, only : r15,pi
 implicit none
 real(r15), intent(in) :: kp(3),gshel(3,26)
 real(r15), intent(out) :: foldedk(3)
 real(r15) kr(3),fkr(3) !,gmax  prim2cart(3,3) !cart2prim(3,3) ,
 integer i

    call bring_to_ws (kp,gshel,foldedk)

 end subroutine fold_in_WS
!========================================
 subroutine fold_in_WS_all(nk,kp,gshel,foldedk)
!! for the BZ defined by (gshel) this subroutine generates foldedk obtained
!! from folding of kp into the WS Wigner-Seitz cell of FBZ,
 use lattice , only : cart2red_g ! prim_to_cart  !r01, r02, r03,
 use geometry
 use kpoints  !, only : nkc,kpc
 use constants, only : r15,pi
 implicit none
 integer, intent(in) :: nk
 real(r15), intent(in) :: kp(3,nk),gshel(3,26)
 real(r15), intent(out) :: foldedk(3,nk)
 real(r15) kr(3),fkr(3) !,gmax  prim2cart(3,3) !cart2prim(3,3) ,
 integer i

 open(134,file='kfolding.dat')
 write(134,*)'# i       kp(:,i)    k_red      folded_k(:,i)   folded_k_red'
 do i=1,nk
! need to fold kp in the FBZ
!   call fold_in_bz_new(kp(:,i),gshel,foldedk(:,i))
    call bring_to_ws_g  (kp(:,i),gshel,foldedk(:,i))
    kr = cart2red_g(kp(:,i))  !matmul(transpose(prim_to_cart),kp(:,i))/(2*pi)
    fkr = cart2red_g(foldedk(:,i))  !matmul(transpose(prim_to_cart),foldedk(:,i))/(2*pi)
    write(134,4)i,kp(:,i),kr,foldedk(:,i),fkr
 enddo
 close(134)
4 format(i6,4(3x,3(1x,f9.4)))

 end subroutine fold_in_WS_all
!========================================
 subroutine fold_in_bz_old(q,gshel,qin)
!! takes the reduced coordinates and puts it in the WS of a cube of length=1
!! then converts the reduced coordinates back to real ones
!! folds the vector q into the WS cell of FBZ and stores the result in qin
! use lattice
 use constants, only : r15,pi
 use ios , only : ulog
 use geometry
 use lattice, only : prim_to_cart, cart_to_prim
 implicit none
 real(r15), intent(in):: gshel(3,26)
 real(r15), intent(in) :: q(3)
 real(r15), intent(out) :: qin(3)
 integer j,nboundary
 real(r15) leng,b1(3),b2(3)
 logical inside

! bring q in [-0.5,0.5] times G1,G2,G3
 b1=matmul(transpose(prim_to_cart),q)/(2*pi) ! reduced coordinates of q
 b2=b1-nint(b1)  ! bring it in [-0.5:0.5]
 qin=matmul(cart_to_prim,b2)*2*pi  ! back to cartesian coords

 do j=1,14 
    if( dot_product(qin,gshel(:,j)) .gt. dot_product(gshel(:,j),gshel(:,j))/2 ) qin=qin-gshel(:,j)
 enddo

! verify it is inside and the rest of 12 vectors do not change the answer
 inside= .false.
 call check_inside_ws(qin,gshel,inside,nboundary)
 write(ulog,5)' inside,q,q_red=',merge(1,0,inside),qin, matmul(transpose(prim_to_cart),qin)/(2*pi) 
 if(.not. inside) then
    write(ulog,*)' now trying the remaining 12 vectors'  
    do j=15,26 
       if( dot_product(qin,gshel(:,j)) .gt. dot_product(gshel(:,j),gshel(:,j))/2 ) then
          qin=qin-gshel(:,j)
          write(ulog,4)'j, new q,q_red=',j, qin,matmul(transpose(prim_to_cart),qin)/(2*pi) 
       endif
    enddo
 endif


      return


      if ( dot_product(b2,(/ 0, 1, 1/)).gt.2/2 ) b2=b2-(/ 0, 1, 1/)
      if ( dot_product(b2,(/ 1, 0, 1/)).gt.2/2 ) b2=b2-(/ 1, 0, 1/)
      if ( dot_product(b2,(/ 1, 1, 0/)).gt.2/2 ) b2=b2-(/ 1, 1, 0/)
      if ( dot_product(b2,(/ 0,-1,-1/)).gt.2/2 ) b2=b2-(/ 0,-1,-1/)
      if ( dot_product(b2,(/-1, 0,-1/)).gt.2/2 ) b2=b2-(/-1, 0,-1/)
      if ( dot_product(b2,(/-1,-1, 0/)).gt.2/2 ) b2=b2-(/-1,-1, 0/)

      if ( dot_product(b2,(/ 0,-1, 1/)).gt.2/2 ) b2=b2-(/ 0,-1, 1/)
      if ( dot_product(b2,(/-1, 0, 1/)).gt.2/2 ) b2=b2-(/-1, 0, 1/)
      if ( dot_product(b2,(/ 1,-1, 0/)).gt.2/2 ) b2=b2-(/ 1,-1, 0/)
      if ( dot_product(b2,(/ 0, 1,-1/)).gt.2/2 ) b2=b2-(/ 0, 1,-1/)
      if ( dot_product(b2,(/ 1, 0,-1/)).gt.2/2 ) b2=b2-(/ 1, 0,-1/)
      if ( dot_product(b2,(/-1, 1, 0/)).gt.2/2 ) b2=b2-(/-1, 1, 0/)

      if ( dot_product(b2,(/ 1, 1, 1/)).gt.3/2 ) b2=b2-(/ 1, 1, 1/)
      if ( dot_product(b2,(/-1, 1, 1/)).gt.3/2 ) b2=b2-(/-1, 1, 1/)
      if ( dot_product(b2,(/ 1,-1, 1/)).gt.3/2 ) b2=b2-(/ 1,-1, 1/)
      if ( dot_product(b2,(/ 1, 1,-1/)).gt.3/2 ) b2=b2-(/ 1, 1,-1/)
      if ( dot_product(b2,(/ 1,-1,-1/)).gt.3/2 ) b2=b2-(/ 1,-1,-1/)
      if ( dot_product(b2,(/-1, 1,-1/)).gt.3/2 ) b2=b2-(/-1, 1,-1/)
      if ( dot_product(b2,(/-1,-1, 1/)).gt.3/2 ) b2=b2-(/-1,-1, 1/)
      if ( dot_product(b2,(/-1,-1,-1/)).gt.3/2 ) b2=b2-(/-1,-1,-1/)

      write(456,4) 'q_red,qin_red=',b1,b2
      qin=matmul(cart_to_prim,b2)*2*pi  ! back to cartesian coords

 inside= .false.
 call check_inside_ws(qin,gshel,inside,nboundary)
 if(inside) then
    return 
 else
    write(ulog,4)'fold_in_bz_new: starting q,qred=',q,b1
    write(ulog,4)'fold_in_bz_new: final  qin,qinred=',qin,b2
    write(*,*)'fold_in_bz_new: still not inside FBZ!'
    write(ulog,*)'fold_in_bz_new: still not inside FBZ!'
 !   stop
 endif

 return
! qin=q
! leng=length(q)
 do j=1,26 !while (leng.gt.gmin)
    if( dot_product(qin,gshel(:,j)) .gt. dot_product(gshel(:,j),gshel(:,j))/2 ) qin=qin-gshel(:,j)
!   if(length(q-gshel(:,j)).lt.leng) then
!      qin=q-gshel(:,j)
!      leng=length(qin)
!   endif
 enddo

 call check_inside_ws(qin,gshel,inside,nboundary)
 if(inside) then
    return
 else
    write(*,4)'fold_in_bz_new: starting q=',q
    write(*,4)'fold_in_bz_new: final  qin=',qin
    write(*,*)'fold_in_bz_new: still not inside FBZ!'
    write(ulog,*)'fold_in_bz_new: still not inside FBZ!'
    stop
 endif

3 format(3(i6),9(2x,f8.4))
4 format(a,9(1x,g11.4))
5 format(a,i5,9(1x,g11.4))

 end subroutine fold_in_bz_old
!-------------------------------------------
  subroutine map_ibz(nk,kp,mapi,ni)
!! checks whether kp is inside IBZ; mapi(jibz)=index of k in FBZ
!! the first ni are in IBZ and the rest is in FBZ
  use constants, only : r15
  implicit none
  integer nk,mapi(nk),i,j,l,ni
  real(r15) kp(3,nk)
  logical inside

  l=0 ; j=0 ; mapi=0
  do i=1,nk
     call check_inside_irred_fbz(kp(:,i),inside)
     if (inside) then
        j=j+1
        mapi(j) = i
     else
        l=l+1
        mapi(nk-l+1) = i
     endif
  enddo
  ni = j
  if(j+l .ne. nk)then
     write(*,*)' pb in map_ibz j,l,nk=',j,l,nk
  endif
  end subroutine map_ibz
!----------------------------------------------
 subroutine foldall_in_fbz_new(nk,kp,gg,nbz,kz,mapz)
!! folds the k-mesh defined by kp into the FBZ and stores the result in kz
 use lattice
 use ios
 use constants, only : r15
 implicit none
 integer, intent(in) :: nk
 integer, intent(out) :: nbz,mapz(nk)
 real(r15), intent(in) :: kp(3,nk),gg(3,26)
 real(r15), intent(out) :: kz(3,nk)
 integer i,ns,j,nboundary
 logical inside
 real(r15) qaux(3),dmin,dist


! write(ulog,*)'FOLD IN FBZ: ik,shell,i,q ====================='
 nbz=0
 do i=1,nk
 foldloop: do ns=1,26
     qaux = kp(:,i)+gg(:,ns)
     call check_inside_ws(qaux,gg,inside,nboundary)
     if(inside) then
! make sure it is different from previous points: find smallest neighbor distance
        dmin = 1d9
        jloop: do j=1,nbz
           dist = length(qaux-kz(:,j))
           if( dist .lt. dmin) then
              dmin = dist
           endif
        enddo jloop
! if smallest distance is not 0 then accept the point
        if ( dmin .gt. 1.d-6 ) then
           nbz=nbz+1
           mapz(nbz) = i
           kz(:,nbz) = qaux(:)
!          write(ulog,3)i,ns,nbz,qaux
        endif
     endif
 enddo foldloop
 enddo
 write(ulog,*)'FOLD IN FBZ:  DONE! ==============================='

3 format(3(i6),9(2x,f8.4))

 end subroutine foldall_in_fbz_new
!----------------------------------------------
  subroutine map_ibz2(nk,kp,mapi,ni)
 use constants, only : r15
  implicit none
  integer nk,mapi(nk),i,j,l,ni
  real(r15) kp(3,nk)
  logical inside

  l=0 ; j=0 ; mapi=0
  do i=1,nk
     call check_inside_irred_fbz(kp(:,i),inside)
     if (inside) then
        j=j+1
        mapi(j) = i
     else
        l=l+1
        mapi(nk-l+1) = i
     endif
  enddo
  ni = j
  if(j+l .ne. nk)then
     write(*,*)' pb in map_ibz j,l,nk=',j,l,nk
  endif
  end subroutine map_ibz2
!----------------------------------------------
  subroutine get_kindex(q,nkc,kpc,ik)
  use constants, only : r15
  use  geometry
  implicit none
  integer ik,nkc,i,j
  real(r15) kpc(3,nkc),q(3)

  ik=0
  mainlp: do i=1,nkc
     do j=1,3
        if (.not. (abs(q(j)-kpc(j,i)) .myeq. 0d0) ) exit
        if (j.eq.3) then
           ik=i
           return !exit mainlp
        endif
     enddo
  enddo mainlp
  if (ik.eq.0) then
     write(*,*)'GET_INDEX: could not find the index of',q
     stop
  endif

  end subroutine get_kindex
!========================================
 subroutine get_kpfbz(nk,kp,gshel,kpfold,fn)
!! this routine translates kp grid into the formal FBZ. Output to KPFBZ.DAT
! use kpoints !, only : nkc,kpc
 use geometry
 use constants, only : r15
 use lattice, only : fold_ws,cart2red_g
 implicit none
 integer, intent(in) :: nk
 real(r15), intent(in) :: gshel(3,26),kp(3,nk)
 real(r15), intent(out) :: kpfold(3,nk)
 integer i
 character(*) fn
 integer :: unit = 793

 do i=1,nk
!   call fold_in_bz_old(kp(:,i),gshel,kpfold(:,i))
   kpfold(:,i)=fold_ws(kp(:,i),gshel,'g')
 enddo

! open(unit,file='KPFBZ.DAT',status='unknown')
 open(unit,file=fn)
 write(unit,'(A)') &
 & "#   n    folded kp(:,n)  reduced_folded kp(:,n)    kp(:,n)  reduced_kp(:,n)"
 do i=1,nk
      write(unit,'(I10,4(3x,3(1x,f9.4)))') i, &
 &    kpfold(:,i),cart2red_g(kpfold(:,i)), kp(:,i), cart2red_g(kp(:,i))
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
 real(r15), intent(in) :: g1(3),g2(3),g3(3)
 real(r15), intent(OUT) :: q(3,3)
 real(r15)              :: g(3,3),r(3,3),k(3),v(3),q2(3),q2max,k2,r1(3),r2(3),r3(3)
 integer              :: i,j,i1,i2,i3,index(3)
 real(r15) :: delta = 1d-6
 real(r15) :: zero2 = 1d-12

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
 subroutine get_weights(nk,kp ,mapinv,mapibz)
!! takes as input kp(3,nk) in cartesian and finds based on crystal symmetries, the
!! weights associated with each kpoint, and then normalizes them
!! nibz is the final number of kpoints stored in kibz in the irreducible BZ
!! the other important output is the mapping of the full kmesh onto the ones
!! folded in the irreducible zone : mapibz(j=1:nk) is the index of the k in the IBZ list
!! corresponding to the argument j
!! for i=1,nibz mapinv(i) gives the index of the corresponding kpoint in the original kp list
 use lattice, only : primitivelattice,r01,r02,r03 ,cart2red_g,red2cart_g,prim_to_conv,fold_ws ,g0ws26
 use kpoints, only : kibz, wibz, nibz
 use constants, only : pi,r15
 use geometry
 use params
 use ios
 implicit none
 integer, intent(in):: nk
 real(r15), intent(in):: kp(3,nk)
 integer, intent(out):: mapibz(nk),mapinv(nk)
 real(r15), allocatable :: k2(:,:),w2(:)
 real(r15) kvecstar(3,48),kf(3),kred(3),weig,sw,q(3)
 integer i,j,s,k,narms,kvecop(48),aux,ibz

 write(*,*) ' entering get_weights'
 open(uibz,file='KPOINT.IBZ',status='unknown')
 open(ufbz,file='KPOINT.FBZ',status='unknown')
 do i=1,nk
    mapinv(i)=i
 enddo
 allocate(w2(nk),k2(3,nk))
 write(*,*) ' w2,k2,mapinv,mapibz allocated'

 nibz=0; w2=0
 iloop: do i=1,nk
    kred=cart2red_g(kp(:,i))
! first fold
 !  kf  = fold_ws(kp(:,i),g0ws26) 
 !  kred=cart2red_g(kf) 
! find stars of kp in reduced units 
    call getstarp(kred,narms,kvecstar,kvecop)

!   write(*,*) ' called getstar i,narms,kred=',i,narms,kred
!   do s=1,narms
!      write(*,3)s,kvecstar(:,s),red2cart_g(kvecstar(:,s))
!   enddo

! compare the stars of i to existing IBZ vectors and skip i if found
! k2 contains the list of IBZs in reduced units
    armloop: do s=1,narms
       do k=1,nibz  ! kvecstar and k2 list are in reduced coordinates
          q=(kvecstar(:,s)-k2(:,k))  
 !        write(*,5)'i,arm,k,nibz,s(k(i))-kibz=',i,s,k,nibz,q
          if (is_integer(q)) then
               mapibz(i)=k 
               w2(k)=w2(k)+1        ! weight of this ibz should be increased
               cycle iloop
           endif
       enddo
! star(s) was not in the ibz list, try other stars
    enddo armloop

! none of the stars was in the ibz list; add i to the ibz list
    nibz=nibz+1
    k2(:,nibz)=kred ! save as reduced coordinates
    mapibz(i) = nibz
    w2(nibz)  = 1
! here need to switch two elements of mapinv; so we add this line
    aux=mapinv(nibz)
    mapinv(nibz)=i
    mapinv(i)=aux
    write(*   ,5)'* GET_WEIGHTS: New vector; nibz,i,mapinv(i), mapibz(i), k in reduced units :', nibz,i,mapinv(i),mapibz(i),kred 
    write(ulog,5)'* GET_WEIGHTS: New vector; nibz,i,mapinv(i), mapibz(i), k in reduced units :', nibz,i,mapinv(i),mapibz(i),kred 
       
 enddo iloop

 write(*,*)'exited iloop with NIBZ=',nibz,' and allocating now WIBZ and KIBZ '
 if ( allocated(wibz)) deallocate(wibz)
 if ( allocated(kibz)) deallocate(kibz)
 allocate(wibz(nibz),kibz(3,nibz))
 wibz=w2(1:nibz)
 sw = sum(wibz(1:nibz))
 wibz = wibz/sw

! of all the stars of k2 select the one whose reduced coords satisfies 0<q1<q2<q3
 do i=1,nibz
    call getstarp(k2(:,i),narms,kvecstar,kvecop)  ! k2 is already reduced
    sloop: do s=1,narms
       q=kvecstar(:,s)
       if ( 0.le.q(1) .and. q(1).le.q(2) .and. q(2).le.q(3) ) then
          kibz(:,i)=red2cart_g(q)
          exit sloop
       endif
! take the first if the above not satisfied
       kibz(:,i)=red2cart_g(kvecstar(:,1)) 
    enddo sloop     
 enddo 

! k(mapinv(1:nibz) is kibz; all in reduced coordinates
! indexing is ik = i1 + n1*i2 + n1*n2*i3 where i1=(0:n1-1) etc

! write the kpoints and their weights-------------------------
 write(uibz,*)'#i,kibz(:,i),kibz_reduced,wibz(i)'
 write(ufbz,*)'#i,kp(i),kp_reduced,folded_WS(kp(i)),reduced_folded_WS(kp(i))'
 open(345,file='mapinv.dat')
 write(345,*)'# i,mapinv(i)(i=1,nfbz),mapibz(i)'
 do i=1,nk !,nibz
    if(i.le.nibz) then
 !     kibz(:,i)=red2cart_g(k2(:,i)) 
 !     kred= fold_ws(kibz(:,i),g0ws26)  kibz is already folded
       write(uibz,3)i,kibz(:,i),cart2red_g(kibz(:,i)),  wibz(i) ! kred id dummy!
    endif
! also fold kp and kibz into in FZ and write 
    kf  = fold_ws(kp(:,i),g0ws26,'g') 
    write(ufbz,3)i,kp(:,i),cart2red_g(kp(:,i)), kf,cart2red_g(kf)
    write(345,*)i,mapinv(i),mapibz(i)
 enddo
 close(345)
 close(uibz)
 close(ufbz)
 deallocate(w2,k2)

3 format(i7,9(3x,3(1x,f9.4)))
4 format(a,2i7,2x,f9.5,2x,99(1x,f9.5))
5 format(a,4i7,2x,99(1x,f9.5))

 end subroutine get_weights
!==========================================================
 subroutine get_weights2(nk,kp ,mapibz)
!! takes as input kp(3,nk) in cartesian and finds based on crystal symmetries, the
!! weights associated with each kpoint, and then normalizes them
!! nibz is the final number of kpoints stored in kibz in the irreducible BZ
!! the other important output is the mapping of the full kmesh onto the ones
!! folded in the irreducible zone : mapibz(j=1:nk) is the index of the k in the IBZ list
!! corresponding to the argument j
!! for i=1,nibz mapinv(i) gives the index of the corresponding kpoint in the original kp list
 use lattice, only : primitivelattice,r01,r02,r03 ,cart2red_g,red2cart_g,prim_to_conv,fold_ws ,g0ws26
 use kpoints, only : kibz, wibz , nibz,wk
 use constants, only : pi,r15
 use geometry
 use params
 use ios
 implicit none
 integer, intent(in):: nk
 real(r15), intent(in):: kp(3,nk)
 integer, intent(out):: mapibz(nk)
! real(r15), intent(out), allocatable :: kibz(:,:) ,wibz(:)
 integer, allocatable :: mcor(:),mapinv(:)
 real(r15), allocatable :: k2(:,:),lg(:),w2(:)
 real(r15) q(3),kvecstar(3,48),sw,skc(3),qout(3),kfold(3,nk),weig
 integer i,j,l,k,narms,kvecop(48),aux
 logical processed(nk)

 open(uibz,file='KPOINT.IBZ',status='unknown')
 open(ufbz,file='KPOINT.FBZ',status='unknown')

 allocate(kibz(3,nk),wibz(nk),mapinv(nk))
 wibz=0d0  
 write(ulog,*)'GET_WEIGHTS2: folding the kpoints in the WS cell '
 do i=1,nk
    kfold(:,i) = fold_ws(kp(:,i),g0ws26,'g')
    wk(i)=1d0/nk
    write(ufbz,3)i,kfold(:,i),wk(i),cart2red_g(kfold(:,i))
    mapinv(i)=i
 enddo
 close(ufbz)
 processed=.False.

 do i=1,nk
 if(.not.processed(i)) then
! find stars of folded vectors
    call getstarp(cart2red_g(kfold(:,i)),narms,kvecstar,kvecop)

    weig=0
! mark all kpoints in the stars as processed
    armloop: do l=1,narms
!  take the star vector to cartesian coordinates and fold into WS cell
       qout = fold_ws(red2cart_g(kvecstar(:,l)),g0ws26,'g')
! compare the stars to existing IBZ vectors and mark processed if found
       do k=1,nk 
          if(.not.processed(k)) then
            if (kfold(:,k) .myeq. qout) then
               processed(k)=.true.
               mapibz(k)=nibz+1  ! use nibz+1 as it will be incremented later
               weig=weig+1
! here need to switch two elements of mapinv; so we add this line
               aux=mapinv(nibz+1)
               mapinv(nibz+1)=k
               mapinv(k)=aux
            endif
          endif
       enddo
    enddo armloop

    nibz=nibz+1
    kibz(:,nibz)=kfold(:,i)
    wibz(nibz)=dble(weig)
    q = cart2red_g(kibz(:,nibz))    ! to get it in reduced coordinates of g0i
    mapibz(i)=nibz
    write(ulog,4)'new vector*:',nibz,narms,wibz(nibz),kibz(:,nibz),length(kibz(:,nibz)),q
    write(*   ,4)'new vector*:',nibz,narms,wibz(nibz),kibz(:,nibz),length(kibz(:,nibz)),q
 endif
 enddo

 kibz=kibz(:,1:nibz)
 wibz=wibz(  1:nibz)
 sw = sum(wibz(1:nibz))
 wibz = wibz/sw

! indexing is ik = i1 + n1*i2 + n1*n2*i3 where i1=(0:n1-1) etc

! sort kibz(1:nibz) in order of increasing length 
! allocate(k2(3,nibz),lg(nibz),mcor(nibz),w2(nibz))
! do i=1,nibz
!    lg(i)=length(k2(:,i))
! enddo
! call sort(nibz,lg,mcor,nibz)

! define new kibz and wibz ---------------------------------------------
! do i=1,nibz
!    wibz(i) = w2(mcor(i))
!    kibz(:,i)=k2(:,mcor(i))
! enddo
! deallocate(k2,w2,lg,mcor)

! sort and write out the kpoints and their weights-------------------------
 write(uibz,*)'#i,kibz(:,i),wibz(i),kibz_reduced,length(kibz(i))'
 open(345,file='mapinv.dat')
 write(345,*)'# i,mapinv(i)(i=1,nibz),mapibz(i)'
 do i=1,nk !,nibz
    if(i.le.nibz) then
       write(uibz,3)i,kibz(:,i),wibz(i),cart2red_g(kibz(:,i)),length(kibz(:,i))
    endif
    write(345,*)i,mapinv(i),mapibz(i)
 enddo
 close(345)
 close(uibz)

3 format(i7,2x,3(1x,f9.4),2x,f9.5,2x,3(1x,f9.5),3x,f9.5)
4 format(a,2i7,2x,f9.5,2x,99(1x,f9.5))

 end subroutine get_weights2
!==========================================================
 subroutine get_weights3(nk,kp ,mapibz)
!! takes as input kp(3,nk) in cartesian and finds based on crystal symmetries, the
!! weights associated with each kpoint, and then normalizes them
!! nibz is the final number of kpoints stored in kibz in the irreducible BZ
!! the other important output is the mapping of the full kmesh onto the ones
!! folded in the irreducible zone : mapibz(j=1:nk) is the index of the k in the IBZ list
!! corresponding to the argument j
!! for i=1,nibz mapinv(i) gives the index of the corresponding kpoint in the original kp list
 use lattice, only : primitivelattice,r01,r02,r03 ,cart2red_g,red2cart_g,prim_to_conv,fold_ws ,g0ws26
 use kpoints, only : kibz, wibz , nibz
 use constants, only : pi,r15
 use geometry
 use params
 use ios
 implicit none
 integer, intent(in):: nk
 real(r15), intent(in):: kp(3,nk)
 integer, intent(out):: mapibz(nk)
! real(r15), intent(out), allocatable :: kibz(:,:) ,wibz(:)
 integer, allocatable :: mcor(:),mapinv(:)
 real(r15), allocatable :: k2(:,:),lg(:),w2(:)
 real(r15) zro,q(3),kvecstar(3,48),sw,skc(3),qout(3)
 integer i,j,l,narms,kvecop(48),aux
 logical exists

!  interface v2a
!     module procedure fold_ws
!  end interface
 open(uibz,file='KPOINT.IBZ',status='unknown')

 allocate(k2(3,nk),w2(nk),mapinv(nk))
 zro=0d0  ! shift
 write(ulog,*)'GET_WEIGHTS3: generating kpoints in the irreducible FBZ '
 write(*,*)'MAKE_KP_IBZ: nk,wk(nk),k2(:,nk)'

 nibz=0 ; mapibz=0; w2=1; k2=0

! initialize mapvinv with identity so that later elements can be switched
 do i=1,nk
    mapinv(i)=i
 enddo

! main loop to identify points in the FBZ
 kploop: do i=1,nk

! below the reduced components of kp are needed, output is also in reduced
    call getstarp(cart2red_g(kp(:,i)),narms,kvecstar,kvecop)
  ! call getstar (kp(:,i),prim_to_conv,narms,kvecstar,kvecop)

! see if already exists among the previously assigned ones
    exists=.False.

    armloop: do l=1,narms

        if(verbose) write(ulog,4)'stars in red&cart are:',i,l,kvecstar(:,l),red2cart_g(kvecstar(:,l))

!  take the star vector to cartesian coordinates and fold into WS cell
        qout = fold_ws(red2cart_g(kvecstar(:,l)),g0ws26,'g')

! set weight for the first kpoint where nibz=0
        jloop: do j=1,nibz
          if (k2(:,j) .myeq. qout) then 
             exists=.True.
             w2(j)=w2(j)+1d0
             mapibz(i)=j
             exit armloop
          endif
        enddo jloop

    enddo armloop

!   write(ulog,4)' kpoint, folded one  ',i,kp(:,i),exists,j,k2(:,j)
    if(exists ) then

       cycle kploop

    else

       nibz=nibz+1  ! didn't exist: new kvector in the ibz list

       if (nibz.eq.1) w2(nibz) = 1
       mapibz(i) = nibz

! here need to switch two elements of mapinv; so we add this line
       aux=mapinv(nibz)
       mapinv(nibz)=i
       mapinv(i)=aux

!      call bring_to_ws_gx(kp(:,i),g0ws26,k2(:,nibz))
       k2(:,nibz) = fold_ws(kp(:,i),g0ws26,'g') 

 !     k2(:,nibz)=kp(:,i)
       q = cart2red_g(k2(:,nibz))    ! to get it in reduced coordinates of g0i
       write(ulog,4)'new vector*:',nibz,narms,w2(nibz),k2(:,nibz),length(k2(:,nibz)),q
       write(*   ,4)'new vector*:',nibz,narms,w2(nibz),k2(:,nibz),length(k2(:,nibz)),q

    endif

 enddo kploop
 write(ulog,*)'GET_WEIGHTS3: generated ',nibz,' points in the irreducible FBZ'

! sort k2(1:nibz) in order of increasing length 
 allocate(lg(nibz),mcor(nibz))
 do i=1,nibz
    lg(i)=length(k2(:,i))
 enddo
 call sort(nibz,lg,mcor,nibz)

! define kibz and wibz ---------------------------------------------
 if ( allocated(wibz)) deallocate(wibz)
 if ( allocated(kibz)) deallocate(kibz)
 allocate(wibz(nibz),kibz(3,nibz))
 do i=1,nibz
    wibz(i) = w2(mcor(i))
    kibz(:,i)=k2(:,mcor(i))
 enddo
! k2 is in reduced units
! do j=1,nibz
!    kibz(:,j)=k2(1,j)*g1+k2(2,j)*g2+k2(3,j)*g3
! enddo

 deallocate(k2,w2)
 sw = sum(wibz(1:nibz))
 wibz = wibz/sw

! sort and write out the kpoints and their weights-------------------------
 write(uibz,*)'#i,kibz(:,i),wibz(i),kibz_reduced,length(kibz(i))'
 open(345,file='mapinv.dat')
! open(346,file='NEWKP.dat')
 write(345,*)'# i,mapinv(i)(i=1,nibz),mapibz(i)'
! write(346,*)'# i,        newkp(i)         newkp_reduced(i)'
 do i=1,nk !,nibz
    if(i.le.nibz) then
       j=mcor(i)
       write(uibz,3)i,kibz(:,j),wibz(j),cart2red_g(kibz(:,j)),length(kibz(:,j))
    endif
    write(345,*)i,mapinv(i),mapibz(i)
!    write(346,3)i,kp(:,i),reduce_g(kp(:,i))
 enddo
 close(345)
! close(346)
 close(uibz)

2 format(3i7,2x,3(1x,f9.4),2x,f9.4,2x,3(1x,f9.5),3x,f14.8)
3 format(i7,2x,3(1x,f9.4),2x,f9.5,2x,3(1x,f9.5),3x,f9.5)
4 format(a,2i7,2x,f9.5,2x,99(1x,f9.5))
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
 use lattice, only : primitivelattice  !,r01,r02,r03
! use kpoints, only : kibz, wibz
 use constants
 use geometry
 use params
 use ios
 implicit none
 integer, intent(in):: nk
 real(r15), intent(in):: kp(3,nk) !,prim2cart(3,3) !r1(3),r2(3),r3(3)
 integer, intent(out):: ngibz,mapibz(nk)
 real(r15), intent(out), allocatable :: gibz(:,:) ,wgibz(:)
 integer, allocatable :: mapinv(:)  !mcor(:),mapibz(:),
 real(r15), allocatable :: k2(:,:),w2(:)
 real(r15) zro,q(3),kvecstar(3,48),sw,skc(3),rr1(3),rr2(3),rr3(3)
 integer i,j,l,narms,kvecop(48),aux
 logical exists


 open(uibz,file='GMESH.IBZ',status='unknown')

 allocate(k2(3,nk),w2(nk),mapinv(nk))
 zro=0d0  ! shift
 write(ulog,*)'GET_WEIGHTS4: generating kpoints in the irreducible FBZ '
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
 write(ulog,*)'GET_WEIGHTS4: generated ',ngibz,' points in the irreducible FBZ'

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
 subroutine structure_factor_reg(v,n,grid,st,dst)
!! takes kpoints within the FBZ as input, and calculates sum_k cos(k.R) for given R
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(in) :: v(3),grid(3,n)
 real(r15), intent(out) :: st,dst(3)
 integer i

 st=0;dst=0
 do i=1,n
     st= st+cos(dot_product(v,grid(:,i))) / n
    dst=dst-sin(dot_product(v,grid(:,i))) * grid(:,i) / n
 enddo

!  write(*,5)'Structure factor for r=', r,st,dst
5 format(a,9(1x,f10.4),3x,g12.5)

 end subroutine structure_factor_reg
!====================================
 subroutine structure_factor(r,nk,kp,wei,st,dst)
!! takes kpoints within the FBZ as input, and calculates sum_k cos(k.R) for given R
 use constants, only : r15
 implicit none
 integer, intent(in) :: nk
 real(r15), intent(in) :: r(3),kp(3,nk),wei(nk)
 real(r15), intent(out) :: st,dst(3)
 integer i

 st=0;dst=0
 do i=1,nk
     st= st+cos(dot_product(r,kp(:,i))) * wei(i)
    dst=dst-sin(dot_product(r,kp(:,i))) * wei(i) * kp(:,i)
 enddo

!  write(*,5)'Structure factor for r=', r,st,dst
5 format(a,9(1x,f10.4),3x,g12.5)

 end subroutine structure_factor
!====================================
  subroutine structure_factor_recip(q,ngrid,grid,weig,sf,dsf)
 !! takes primitive vectors within the WS of the supercell as input, and calculates
 !! 1/N_R sum_R cos(q.R) for given q
  use lattice, only : cart_to_prim ,prim_to_cart !,volume_r,volume_r0
  use constants, only : pi,r15
  use params, only : tolerance
  implicit none
  integer, intent(in) :: ngrid
  real(r15), intent(in) :: q(3),grid(3,ngrid),weig(ngrid)
  real(r15), intent(out) :: sf,dsf(3)
  integer i
  real(r15) phase

!  if(abs(sum(weig)-1).gt.tolerance) then
!     write(*,*)'STRUCTURE_FACTOR_RECIP: weights should add up to 1, not ',sum(weig)
!!! quick fix: for now we rescale weights to add to 1
!     scal=1/sum(weig)
!     scalw=weig*scal
!!! quick fix: for now we rescale weights to add to 1
!  endif


  sf=0; dsf=0
!   write(*,4)'Structure factor_recip q,qred=',q,matmul(transpose(prim_to_cart),q)/(2*pi)
  do i=1,ngrid
     phase = dot_product(q,grid(:,i)) 
     sf = sf + cos(phase) * weig(i)
!     write(*,5)'Structure factor_recip for i,phase/2pi,cos(phase),sf,grd_red ',i,phase/(2*pi),cos(phase),sf, &
!&               matmul(cart_to_prim,grid(:,i))
    dsf = dsf - sin(phase) * weig(i) * grid(:,i) 
  enddo
  sf = sf/sum(weig)
  dsf=dsf/sum(weig)

 4 format(a,3(1x,f7.3),3x,3(1x,f7.3))
 5 format(a,i3,3(1x,f9.4),3x,3(1x,f7.3))

  end subroutine structure_factor_recip
 !====================================
  subroutine structure_factor_complx(q,ngrid,grid,weig,sf,dsf)
 !! takes primitive vectors within the WS of the supercell as input, and calculates
 !! 1/N_R sum_R cos(q.R) for given q
  use lattice, only : cart_to_prim ,prim_to_cart
  use constants, only : pi,r15,ci
  use params, only : tolerance
  implicit none
  integer, intent(in) :: ngrid
  real(r15), intent(in) :: q(3),grid(3,ngrid),weig(ngrid)
  complex(r15), intent(out) :: sf,dsf(3)
  integer i
  real(r15) phase

  sf=0; dsf=0
  do i=1,ngrid
     phase = dot_product(q,grid(:,i)) 
     sf = sf + exp(ci*phase) * weig(i)
    dsf = dsf +exp(ci*phase) * weig(i) * grid(:,i) * ci
  enddo

  end subroutine structure_factor_complx
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
 use constants, only : r15
! use constants , only : pi
 implicit none
 type(vector), intent(in) :: g01,g02,g03,g1,g2,g3
 real(r15), intent(in) :: gshells(3,26)
 integer, intent(in) :: nk
 real(r15), intent(out) :: kmesh(3,nk)
 real(r15), allocatable:: aux(:,:)
 real(r15) q(3),nij(3,3),qred(3)
 type(vector) kred,rr1,rr2,rr3
 integer i,j,k,l,mxi,cnt,mxa,nboundary
 logical exists,inside

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
    call check_inside_ws(q,gshells,inside,nboundary)
    if(inside ) then
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
 subroutine is_on_boundary(v,grid,indx,tol)
!! checks whether the vector v is on the FBZ WS cell boundary within the tolerance
!! where FBZ is defined by the grid(3,26)
!! indx is the number of boundaries to which v belongs
!! indx=0 means not on boundary
 use geometry
 use constants, only : r15
 implicit none
 real(r15), intent(in) :: v(3),grid(3,26),tol
 integer, intent(out) :: indx
 integer j
 real(r15) vgoverg2

 indx=0
 do j=1,13
    if(length(grid(:,j)).lt.tol) then
        write(*,*)'grid vector j is zero!! ',j,grid(:,j)
        stop
    endif
    vgoverg2 = 2 * dot_product(v,grid(:,j))/dot_product(grid(:,j),grid(:,j))
    if ( abs(abs(vgoverg2) - 1) .lt. tol ) then
       indx=indx+1
    endif
 enddo

 end subroutine is_on_boundary
!===========================================================
 subroutine check_periodic(n,grid,r26,array,aravg,isperiodic)
!! check whether the array(n) defined on the grid(3,n) is periodic of period one of r26
!! if not, aravg will contain the symmetrized array, which is periodic 
 use constants, only : r15
 use params, only : tolerance
 use geometry, only : length
 use ios, only: ulog,write_out
 implicit none
 integer, intent(in):: n
 real(r15), intent(in) :: grid(3,n),r26(3,26)
 complex(r15), intent(in) :: array(n)
 complex(r15), intent(out) :: aravg(n)
 complex(r15) suma
 logical, intent(out):: isperiodic
! nbound(i) the the number of boundaries the grid(i) is on; cnt(i,nb) contains the label of r26 boundary vectors
! cnt(i,10) is the  counter of boundary indices for grid(i)
 integer i,j,l,nbound(n),cnt(n,20),m,nboundr
 real(r15) vgoverg2,v(3),arscale,tol

 arscale=sum(abs(array))/dble(n)
 tol=tolerance*arscale
 aravg=array
 isperiodic=.true.
 do i=1,n
    call is_on_boundary(grid(:,i),r26,nboundr,tol)
! check whether grid(i) is on the boundary
    nbound(i)=0 ; cnt(i,:)=0 ; 
    do j=1,26  ! all 26 are needed ; direction is important!
       vgoverg2 = 2 * dot_product(r26(:,j),grid(:,i))/dot_product(r26(:,j),r26(:,j))
       if ( abs(vgoverg2 - 1) .lt. tol ) then ! it's on the booundary
          nbound(i)=nbound(i)+1   ! count on how many boundaries
          cnt(i,nbound(i))=j    
       endif
    enddo
!    write(ulog,*)'site ',i,' is on ',nbound(i),nboundr,' boundaries corresponding to vectors ',cnt(i,1:nbound(i))
 enddo   

!   if(nbound.eq.0) cycle ! was not on the boundary; try next grid point

! it is on the intersection of nbound surfaces 
! find array indices(&values) on those boundaries to symmetrize
!   cnt(i,11)=i
!   write(ulog,*)' CNT(i,:)=',i,cnt(i,:)
 do i=1,n
    if(nbound(i).eq.0) then ! it's inside but not on boundary
        cycle
    else
        suma=array(i) ! initialize the averaging    
        cnt(i,10)= 1  ! counts the weight; should be = nbound(i)+1 at the end  of the loop; used for averaging
        cnt(i,11)= i  ! cnt(i,10+m) is the index l of mth grid point image of i on another boundary r26(j); j=cnt(i,m)  
    endif
    do m=1,nbound(i)  ! these exclude the grid(i) itself
       v=grid(:,i)-r26(:,cnt(i,m))   
! which grid index corresponds to v?
       lloop: do l=1,n
          if(length(v-grid(:,l)).lt. tolerance) then
! accumulate for averaging
             suma=suma+array(l)
             cnt(i,10)=cnt(i,10)+1 ! counter of images of grid(i) on other boundaries
             cnt(i,10+cnt(i,10))=l ! index of the image grids
             exit lloop
          endif
       enddo lloop
    enddo
    if(nbound(i)+1.ne.cnt(i,10)) then
        write(*,*)'grid i=',i,' nbound+1=',nbound(i)+1,'cnt(i,10)=',cnt(i,10),' diff=',nbound(i)+1-cnt(i,10)
    endif
    aravg(i)=suma/cnt(i,10)  ! this is the average of array(i)
    write(175,3)'grid(i), image, array, final avg,dif=',i,l,array(i),aravg(i),array(i)-aravg(i)
 !  write(ulog,*)'i, images=',i,cnt(i,10), cnt(i,11:10+cnt(i,10)) ,nbound(i)+1, cnt(i,1:1+nbound(i)) ! cnt(i,ind) is the rws index
!   write(ulog,2)' array(l)=',array(cnt(i,11:10+cnt(i,10))) , aravg(i)
    if (abs(array(i)-aravg(i)).gt.1d-6) then
       isperiodic=.false.
    endif
 enddo
! if(.not.isperiodic) then
!     call write_out(ulog,'CHECK_PERIODIC: array is not periodic ',array)
!     call write_out(ulog,'CHECK_PERIODIC:         average array ',aravg)
! endif

2 format(a,9(1x,g11.4))
3 format(a,2i4,9(1x,f10.4))

 end subroutine check_periodic
!===========================================================
 subroutine check_inside_ws(q,gshel,inside,nboundary)
!! checks if q is inside the the Wigner-Seitz cell defined by neighborshells
!! vectors gshel, It is also considered inside if on the WS facet/edge/vertex
 use ios, only : ulog
 use constants, only : r15
 implicit none
 logical, intent(out) :: inside
 integer, intent(out) :: nboundary
 real(r15), intent(in):: q(3),gshel(3,26)
 real(r15) qdotg,gg2,tol
 integer i

 tol=0.002  ! a point on the boundary is also counted as inside

 inside = .false.
 nboundary = 0
! write(ulog,*)'CHECK_INSIDE: q=',q

! this has to be modified for 1D or 2D systems 
 do i=1,13
    qdotg = dot_product(q,gshel(:,i))
    gg2   = dot_product(gshel(:,i),gshel(:,i))/2d0
! write(ulog,4)'i,qdotg,gg/2=',i,qdotg,gg/2,gshel(:,i)
! check the Bragg condition |q.G| < G.G/2
    if(abs(qdotg) .gt. gg2 * (1+tol)) then  ! it falls outside
       nboundary=i
       return  ! exit if it is outside
    elseif(abs(qdotg) .gt. gg2 * (1-tol)) then  ! it's on the boundary
       nboundary= nboundary+1
    endif 
 enddo

! qdotg was =< than ALL gg/2: q is therefore inside or on boundary
 inside = .true.

! call is_on_boundary(q,gshel,nboundary,tol)

4 format(a,i6,9(1x,f11.4))

 end subroutine check_inside_ws
!===========================================================
 subroutine check_inside_ws_noimage(q,gshel,inside,nboundary)
!! checks if q is inside the the Wigner-Seitz cell defined by neighborshells
!! vectors gshel, It is considered inside if on the WS facet/edge/vertex
! use geometry, only : dot
 use ios, only : ulog
 use constants, only : r15
 implicit none
 logical, intent(out) :: inside
 integer, intent(out) :: nboundary
 real(r15), intent(in):: q(3),gshel(3,26)
 real(r15) qdotg,gg2,tol
 integer i

 tol=0.002  ! a point on the boundary is also counted as inside

 inside = .false.
 nboundary = 0
! write(ulog,*)'CHECK_INSIDE: q=',q

 do i=1,13
    qdotg = dot_product(q,gshel(:,i))
    gg2   = dot_product(gshel(:,i),gshel(:,i))/2d0
! write(ulog,4)'i,qdotg,gg/2=',i,qdotg,gg/2,gshel(:,i)
! check the Bragg condition |q.G| < G.G/2
!   if(abs(qdotg) .gt. gg2 * (1+tol)) return  ! exit if it is outside 
    if(qdotg .ge. gg2 .or. qdotg .lt. -gg2 ) return  ! exit if it is outside 
 enddo

! qdotg was smaller than ALL gg/2: q is therefore inside but not on boundary
 inside = .true.

 call is_on_boundary(q,gshel,nboundary,tol)

4 format(a,i6,9(1x,f11.4))

 end subroutine check_inside_ws_noimage
!===========================================================
 subroutine get_26shortest_shell(gg1,gg2,gg3,gshells,x1,x2,x3)
!! get the 26 nonzero shortest vectors to define the FBZ, and also to be used for check_inside_fbz
!! outputs gshells and the shortest basis x1,x2,x3 which may replace gg1,gg2,gg3
 use geometry
 use ios , only : ulog
 use constants, only : r15
 implicit none
 type(vector), intent(in) :: gg1,gg2,gg3
 type(vector), intent(out):: x1,x2,x3  ! the 3 shortest vectors
 type(vector) xg1,xg2,xg3 
 real(r15) , intent(out):: gshells(3,26)
 real(r15) , allocatable :: aux(:,:),leng(:)
 integer , allocatable :: mcor(:)
 real(r15) junk,vol,gt(3),gmax,y1,y2,y3
 integer d3,i,j,k,cnt,radius,diameter

! find the longest G
 gmax=max(length(gg1),length(gg2))
 gmax=max(gmax,length(gg3))
 call calculate_volume(gg1,gg2,gg3,vol)
 if(vol.myeq.0d0) then
    write(*,*)'CALVOL: Volume is zero! check your basis vectors!'
    write(*,'(a,9(1x,g12.5))')'r1,r2,r3=',gg1,gg2,gg3
    stop
 endif
 junk=vol**0.33333   ! typical short length
 radius = (nint(gmax/junk)+1)
 diameter=2*radius+1
 write(*,*)'gmax,vol,junk,radius,diameter=',gmax,vol,junk,radius,diameter
 d3= diameter**3 - 1   ! to exclude 0
 write(*,*)'26SHORTEST_SHELL: diameter**3-1=d3=',d3
 allocate(leng(d3),aux(3,d3),mcor(d3))

! construct a lattice out of Gs and take the 27 shortest vectors
! this has to be modified for 1D or 2D systems or if one G is much larger than others
 cnt=0
 leng=1d5
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

 call sort(cnt,leng,mcor,d3)

! find the shortest 3 vectors
 x1= a2v(aux(:,mcor(1)))

! make sure x1 and "x2" are independent
 j=2
 vol=0
 do while  (abs(vol).lt.1d-7*length(x1)**2)
    x2 = a2v(aux(:,mcor(j))) ! 2nd shortest linearly indept from x1
    x3 = x1 .cross. x2
    vol=length(x3)   ! this is really the area
    j=j+1
 enddo
 write(*,*)'26SHORTEST_SHELL2: j, area=',j,vol

 vol=0
 do while  (abs(vol).lt.1d-10*length(x1)**3)
    x3 = a2v(aux(:,mcor(j)))
    call calculate_volume(x1,x2,x3,vol)
    j=j+1
 enddo
 write(*,*)'26SHORTEST_SHELL3: j, volume=',j,vol
 call calculate_volume(x1,x2,x3,vol)
 write(*,3)'26SHORTEST_SHELL: final smallest vector1=',x1
 write(*,3)'26SHORTEST_SHELL: final smallest vector2=',x2
 write(*,3)'26SHORTEST_SHELL: final smallest vector3=',x3
 write(*,*)'26SHORTEST_SHELL: final smallest volume=',vol
 call make_reciprocal_lattice_v(x1,x2,x3,xg1,xg2,xg3)

3 format(a,3('(',3(1x,f11.4),')'))
4 format(i6,9(1x,f11.4))

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
!call sort(cnt,leng,mcor,26)
 write(ulog,3)' 26 (unsorted) shortest shell vectors for :',x1,x2,x3
 do i=1,cnt
    gshells(:,i)=aux(:,i) !mcor(i))
    y1=dot_product(gshells(:,i),v2a(xg1))
    y2=dot_product(gshells(:,i),v2a(xg2))
    y3=dot_product(gshells(:,i),v2a(xg3))
    write(ulog,4)i,gshells(:,i), length(gshells(:,i)),y1,y2,y3
 enddo

 deallocate(leng,aux,mcor)

 end subroutine get_26shortest_shell
 !============================================================
 function bring_to_prim_cell_cav(a) result(w)
! takes cartesian coordinates and returns cart coordinates within the primcell
 use lattice
 use geometry
 use constants, only : r15
 implicit none
 type(vector) v,w, bring_to_prim_cell_cv
 real(r15) a(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w =  bring_to_prim_cell_cv(v)

 end function bring_to_prim_cell_cav
!-----------------------------------------
 function bring_to_super_cell_cav(a) result(w)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
 use geometry
 use constants, only : r15
 implicit none
 type(vector) v,w, bring_to_super_cell_cv
 real(r15) a(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w = bring_to_super_cell_cv(v)

 end function bring_to_super_cell_cav
!-----------------------------------------
 function bring_to_prim_cell_caa(a) result(b)
! takes cartesian coordinates and returns cart coordinates within the primcell
 use lattice
 use geometry
 use constants, only : r15
 implicit none
 type(vector) v,w, bring_to_prim_cell_cv
 real(r15) a(3),b(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w =  bring_to_prim_cell_cv(v)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end function bring_to_prim_cell_caa
!-----------------------------------------
 function bring_to_super_cell_caa(a) result(b)
! takes cart coordinates and returns cart coordinates within the supercell
 use lattice
 use geometry
 use constants, only : r15
 implicit none
 type(vector) v,w, bring_to_super_cell_cv
 real(r15) a(3),b(3)

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 w = bring_to_super_cell_cv(v)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end function bring_to_super_cell_caa
!-----------------------------------------
 function bring_to_prim_cell_cv(v) result(w)
! takes cartesian coordinates and returns cart coordinates within the primcell
 use lattice
 use geometry
 use constants, only : r15
 implicit none
 type(vector) v,w
 real(r15) a1,a2,a3

! get direct coordinates first
 a1 = v .dot. g01/(2*pi)
 a2 = v .dot. g02/(2*pi)
 a3 = v .dot. g03/(2*pi)
! bring into cell: between 0 and cell
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
! convert to cartesian coordinates
 w = a1*r01+a2*r02+a3*r03

 end function bring_to_prim_cell_cv
!-----------------------------------------
 function fold_to_WS_super_cell_r(v) result(w)
! takes cart coordinates and returns cart coordinates within the supercell 0=< v < R_sc
 use lattice
 use geometry
 use constants, only : pi,r15
 implicit none
 real(r15), intent(in) :: v(3)
 real(r15) a1,a2,a3,w(3),d(3)

! get direct coordinates first
! d=cart2red_rs(v)
 a1 = (v .dot. gs1)/(2*pi)
 a2 = (v .dot. gs2)/(2*pi)
 a3 = (v .dot. gs3)/(2*pi)
! fold to WS cell
 a1 = a1 - nint(a1)
 a2 = a2 - nint(a2)
 a3 = a3 - nint(a3)
! convert back to cartesian coordinates
 w = a1*rs1+a2*rs2+a3*rs3

 end function fold_to_WS_super_cell_r
!-----------------------------------------
 function bring_to_super_cell_cv(v) result(w)
! takes cart coordinates and returns cart coordinates within the supercell 0=< v < R_sc
 use lattice
 use geometry
 use constants, only : pi,r15
 implicit none
 type(vector) v,w
 real(r15) a1,a2,a3

! get direct coordinates first
 a1 = (v .dot. gs1)/(2*pi)
 a2 = (v .dot. gs2)/(2*pi)
 a3 = (v .dot. gs3)/(2*pi)
! bring into cell: between 0 and cell
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
! convert to cartesian coordinates
 w = a1*rs1+a2*rs2+a3*rs3

 end function bring_to_super_cell_cv
!-------------------------------------
 subroutine find_map_rtau2sc(nr,rmesh,map_rtau_sc)
!! finds the index of the supercell atom corresponding to primitive translation vector
!! defined in rmesh (modulo supercell translations) and primitive atom tau
!! map_rtau_sc(tau,igrid)=k=atom index in supercell
 use lattice , only : gs1,gs2,gs3 !rs1,rs2,rs3,
 use ios , only : ulog
 use geometry
 use constants, only : pi,r15
 use atoms_force_constants, only : atom_sc, atom0, natom_super_cell,natom_prim_cell
 implicit none
 integer, intent(in) :: nr
 real(r15), intent(in) :: rmesh(3,nr)
 integer, intent (out) :: map_rtau_sc(natom_prim_cell,natom_super_cell/natom_prim_cell)
 integer i,tau,k,nsc,l
 type(vector) diff
 real(r15) dist(3),dred(3)
!logical is_integer

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
 use constants, only : r15
  implicit none
  integer n,i,j,temp
  real(r15) r(n)
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
 use constants, only : r15
 implicit none
 integer, intent(in) :: ndim,nbasis
 real(r15), intent(in) :: basis(ndim,nbasis),v(ndim)
 real(r15), intent(out) :: vp(ndim)
 integer c,i,j
 real(r15) res

 ! make sure basis is orthonormal
  do i=1,nbasis
  do j=i,nbasis
     res=dot_product(basis(:,i),basis(:,j))
     if (i.eq.j .and. abs(res-1d0).gt.1d-5 ) then
        write(*,*)'Basis not normalized,i,vi^2=',i,res
     elseif (i.ne.j .and. abs(res).gt.1d-5 ) then
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
 use constants, only : r15
 implicit none
 integer, intent(in) :: ndim,n
 real(r15), intent(in) :: basis(ndim,n),v(ndim)
 real(r15), intent(out) :: vpo(ndim)
 integer c

 vpo=v
 do c=1,n
    vpo = vpo - dot_product(v,basis(:,c)) * basis(:,c)
 enddo

 end subroutine project_out
!===================================================
 subroutine find_first_nonzero(dim,bmat,n_hom)
!! finds the position of the first non-zero element in array bmat
 use constants, only : r15
 implicit none
 integer, intent(in):: dim
 real(r15), intent(in) :: bmat(dim)
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
 use constants, only : r15
 implicit none
 integer, intent(in) :: na,nb
 real(r15), intent(in) :: a(na)
 real(r15), intent(in) :: b(nb)
 real(r15), allocatable, intent(inout) :: c(:)
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
 use constants, only : r15
 implicit none
 real(r15), allocatable, intent(inout) :: a(:)
!real(r15), allocatable, intent(in) :: b(:)
!real(r15), intent(inout) :: a(:)
 real(r15), intent(in) :: b(:)
 real(r15), allocatable :: aux(:)
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
! subroutine append_array(a,b,c)
!!! appends array b to the end of a and stores the resu     lt in c
! implicit none
!! integer, intent(in) :: na,nb
!! real(r15), intent(in):: a(na),b(nb)
!! real(r15), allocatable, intent(inout):: c(:)
!! real(r15) :: c(size(a)+size(b))
! real(r15), intent(in):: a(:),b(:)
! real(r15), allocatable :: c(:)
! real(r15), allocatable :: aux(:)
! integer nc,na,nb

! na=size(a);
! nb=size(b);
! nc=na+nb
! allocate(aux(nc)) ! this is to allow calls like append_array(a,b,a)
! aux=reshape(a,shape=(/na+nb/),pad=b)
! if (allocated(c)) deallocate (c)
! allocate(c(nc))
! c=aux
! deallocate(aux)

!! cal also use in the 1D case:
!! c=reshape(a,(/na+nb/),pad=b,order=(/1,1/))
!
! end subroutine append_array
!===================================================
  subroutine swap(a,b)
 use constants, only : r15
  real(r15), intent(inout) :: a,b
  real(r15) c
  c=a
  a=b
  b=c
  end subroutine swap
!===================================================
      subroutine xmatinv(n,xmatin,xmatout,ier)
 use constants, only : r15
 use ios, only : write_out
      implicit none

!! invert a n by n matrix

!     double precision xmatin(3,3),xmatout(3,3),buffer(3,3),x
      integer, intent(in) :: n
      real(r15), dimension(n,n),intent(in) :: xmatin
      real(r15), dimension(n,n),intent(out) :: xmatout
      integer, intent(out) :: ier
      real(r15) x
      real(r15) buffer(n,n) !size(xmatin,1),size(xmatin,2))
! , allocatable:: buffer(:,:) ! dimension(size(xmatin,1),size(xmatin,2)) :: buffer
      integer i,j,indx(n) !size(xmatin,1)),n,i,j
      external lubksb,ludcmp

! dimension of matrix
!     n=size(xmatin,1)   ! it's a square matrix!
!     n=3
!     write(6,*)'Entering XMATINV with matrix size=',n
!     allocate(buffer(n,n))
! clear error flag
      ier=0
      buffer=xmatin  ! keep matin at exit

      xmatout=0
      do i=1,n
        xmatout(i,i)=1
      enddo
! decomposition
!     write(6,*)'XMATINV: calling LUDCMP for size=',n
      call ludcmp(buffer,n,n,indx,x)
! singular matrix
      if(x.eq.0.0d0)then
        ier=1
        write(*,*)' XMATINV: INVERSION ERROR; ier.ne.0 ',ier
        call write_out(6,'input matrix for inversion ',xmatin)
        stop !return
      endif
! inverse matrix
      do j=1,n
        call lubksb(buffer,n,n,indx,xmatout(:,j))
      enddo
   !  deallocate(buffer)

      end subroutine xmatinv
!------------------------------------------------------------------
! The following routines are from Numerical Recipes
      subroutine ludcmp(a,n,np,indx,d)
 use constants, only : r15
      implicit none
      integer , intent(in) :: np,n
      integer , intent(out) :: indx(n)
      real(r15), intent(inout) :: a(np,np),d
      real(r15) tiny
      parameter (tiny=1.0d-20)
      real(r15) vv(np),aamax,dum,sum
      integer i,j,k,imax,ncmp

      d=1
      do i=1,n
        aamax=0
        do j=1,n
          if(dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
        enddo
        if(ncmp(aamax).eq.0)then
! singular matrix
          d=0
          return
        endif
        vv(i)=1/aamax
      enddo
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if(dum.ge.aamax)then
            imax=i
            aamax=dum
          endif
        enddo
        if(j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.0d0)a(j,j)=tiny
        if(j.ne.n)then
          dum=1/a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
      end subroutine ludcmp
!-----------------------------------------------
      subroutine lubksb(a,n,np,indx,b)
 use constants, only : r15
      implicit none
      integer , intent(in) :: np,n
      integer , intent(in) :: indx(n)
      real(r15), intent(in) :: a(np,np)
      real(r15), intent(inout) :: b(n) 
      real(r15) sumx
      integer ii,i,j,ll

      ii=0
      do i=1,n
        ll=indx(i)
        sumx=b(ll)
        b(ll)=b(i)
        if(ii.ne.0)then
          do j=ii,i-1
            sumx=sumx-a(i,j)*b(j)
          enddo
        else if(sumx.ne.0.0d0)then
          ii=i
        endif
        b(i)=sumx
      enddo
      do i=n,1,-1
        sumx=b(i)
        if(i.lt.n)then
          do j=i+1,n
            sumx=sumx-a(i,j)*b(j)
          enddo
        endif
        b(i)=sumx/a(i,i)
      enddo

      end subroutine lubksb

!===========================================================
subroutine sort3(n1,n2,n3,mp)
!! sorts the 3 integers n1,n2,n3 in increasing order!! that's stupid!
implicit none
integer n1,n2,n3,mp(3)

if( n2.ge.n1) then
  if (n3.ge.n2) then
     mp(1)=n1
     mp(2)=n2
     mp(3)=n3
  elseif( n3.ge.n1) then
     mp(1)=n1
     mp(2)=n3
     mp(3)=n2
  else
     mp(1)=n3
     mp(2)=n1
     mp(3)=n2
  endif
elseif( n3.ge.n1) then
     mp(1)=n2
     mp(2)=n1
     mp(3)=n3
 elseif(n3.ge.n2) then
     mp(1)=n2
     mp(2)=n3
     mp(3)=n1
 else
     mp(1)=n3
     mp(2)=n2
     mp(3)=n1
 endif
end subroutine sort3

!===========================================================
 subroutine check_inside_fbz_old(q,g1,g2,g3,inside)
 use geometry
 use constants
 implicit none
 real(r15), intent(in):: q(3)
 type(vector), intent(in):: g1,g2,g3
 logical, intent(out) :: inside
 real(r15) qdotg,gg,gt(3),sm
 sm = -.0001 * (g1.dot.g1) !dot_product(g1,g1) ! on the boundary is also counted as inside, but once only

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

 end subroutine check_inside_fbz_old
!================================================
 subroutine symmetrize4(n,mat4) 
! symmetrizes wrt to swap of indices 1<->2 and 3<->4 
 use constants, only : r15
 use linalgb
 implicit none
 integer, intent(in) :: n       
 real(r15), intent(inout) :: mat4(n,n,n,n) !:,:,:,:) 
 real(r15), allocatable:: mean(:,:),mea2(:,:) 
 integer i,j

! n=size(mat4,1)
! if(n.ne.3) then
!   write(*,*)' Symmetrize4: n is ',n
!   stop
! endif
 allocate(mean(n,n),mea2(n,n))
 do i=1,n
    mean=mat4(i,i,:,:)
    call symmetrize2(n,mean)
    mat4(i,i,:,:)=mean
 do j=i+1,n
    mean=mat4(i,j,:,:)
    mea2=mat4(j,i,:,:)
    call symmetrize2(n,mean)
    call symmetrize2(n,mea2)
    mat4(i,j,:,:)=(mean+mea2)/2
    mat4(j,i,:,:)=(mean+mea2)/2
 enddo
 enddo
! now wrt ab<->gd
! do i=1,n
! do j=1,n
!    mean=mat4(:,:,i,j)
!    call symmetrize2(n,mean)
!    mat4(:,:,i,j)=mean
! enddo
! enddo

 deallocate(mean)
 end subroutine symmetrize4
!================================================
 subroutine convert_to_voigt(a,c)
 use ios, only: ulog,write_out
 use constants, only : r15
 implicit none
 real(r15), intent(in):: a(3,3,3,3)
 real(r15), intent(out)::c(6,6)
 integer voigt,al,be,ga,de,i,j

 c=0
 do al=1,3
 do be=al,3
 do ga=1,3
 do de=ga,3
    i=voigt(al,be); j=voigt(ga,de)
    c(i,j)=a(al,be,ga,de)
 enddo
 enddo
 enddo
 enddo
 end subroutine convert_to_voigt
!================================================
 subroutine grad(n,q0,f,gf)
!! computes gradient of the scalar function f(q) at point q0 in the n-dimensional space
 use constants, only : r15
 implicit none

 integer, intent(in):: n
 real(r15), intent(in):: q0(n)
 real(r15), intent(out)::gf(n)  
 real(r15), external :: f
 real(r15) dq,fp,fm,qp(n),qm(n)
 integer i

 dq=1d-4
 qp=q0; qm=q0
 do i=1,n
    qp(i)=q0(i)+dq
    qm(i)=q0(i)-dq
    fp=f(qp) ! or call function(qp,fp)
    fm=f(qm)
    gf(i)=(fp-fm)/(2*dq)
    qp(i)=q0(i)-dq
    qm(i)=q0(i)+dq
 enddo
   
 end subroutine grad
!================================================
 subroutine grad_sub(sub,n,q0,f,gf,h)
!! computes gradient of the scalar function f(q) at point q0 in the n-dimensional space
 use constants, only : r15
 implicit none

 interface 
   subroutine sub(n,qin,fout)
   use constants, only : r15
   implicit none
   integer, intent(in) :: n
   real(r15), intent(in) :: qin(n)
   real(r15), intent(out) :: fout
   end subroutine sub
 end interface

 integer, intent(in):: n
 real(r15), intent(in):: q0(n)
 real(r15), intent(in), optional :: h
 real(r15), intent(out)::gf(n)  
 real(r15), external :: f
 real(r15) dq,fp,fm,qp(n),qm(n)
 integer i

 If(present(h)) then
   dq=h
 else
   dq=1d-4  ! default step size
 endif

 qp=q0; qm=q0
 do i=1,n
    qp(i)=q0(i)+dq
    qm(i)=q0(i)-dq
    call sub(n,qp,fp)
    call sub(n,qm,fm)
    gf(i)=(fp-fm)/(2*dq)
    qp(i)=q0(i)-dq
    qm(i)=q0(i)+dq
 enddo
   
 end subroutine grad_sub

!=========================================
 subroutine select_corner(n,ks,k)
!! among the star ks, selects the one with largest kx,ky,kz such that kx<ky<kz
 use constants, only : r15
 use geometry, only : length
  use lattice, only : cart_to_prim,cart2red_g
 implicit none
 integer, intent(in) :: n
 real(r15), intent(in) :: ks(3,n)
 real(r15), intent(out):: k(3)
 real(r15) kred(3,n),leng(n),kc(3)
 integer i,mcor(n),ier

 do i=1,n
    kred(:,i)=cart2red_g(ks(:,i))   ! get its reduced vector
    leng(i)=dot_product(kred(:,i),(/1,1,1/)) ! find longest vector along (1,1,1)
 enddo
 call sort(n,leng,mcor,n)

 ier=1
 iloop: do i=0,n-1
    k=kred(:,mcor(n-i))  ! among the longest parallel to, select the one with kz>ky>kx
    if(k(3).ge.k(2) .and. k(2).ge.k(1)) then
       kc=k
       ier=0
       exit iloop
    endif
 enddo iloop
 if(ier.eq.0) then
   k=matmul(cart_to_prim,kc)   ! back to cartesian
 else
   write(*,*)'SELECT_CORNER MAY HAVE ERROR'
   k=matmul(cart_to_prim,kred(:,mcor(n)))   ! back to cartesian
 endif

 end subroutine select_corner
!=========================================
 subroutine bring_to_ws(qin,gs,qout)
!! takes qin (in cartesian) and applies translations defined by gs to bring it in the WS cell of gs
!! excludes one of the boundaries 0=< g <1
 use constants, only : r15
 use lattice, only : cart2red_g
 implicit none
 real(r15), intent(in) :: gs(3,26),qin(3)
 real(r15), intent(out):: qout(3)
 integer i,nboundary
 logical insid

 qout=qin
 do i=1,26
    do while( 2*dot_product(qout,gs(:,i)).gt.dot_product(gs(:,i),gs(:,i)) )
       qout=qout-gs(:,i)   
       write(*,4)'<- i, q,qred,gred=',i,qout,cart2red_g(qout),cart2red_g(gs(:,i))   
    enddo
 enddo
 write(*,4)'final, q,qred=',i,qout,cart2red_g(qout)

 call check_inside_ws(qout,gs,insid,nboundary) 
 if(.not.insid) then
   write(*,4)'BRING_TO_WS_GX ERROR: final q=',i,qout,cart2red_g(qout)
   stop
 endif

4 format(a,i4,99(1x,f8.4)) 
 end subroutine bring_to_ws
!=========================================
 subroutine bring_to_ws_g(qin,gs,qout)
!! takes qin and applies translations defined by gs to bring it in the WS cell of gs
 use constants, only : r15
 real(r15), intent(in) :: gs(3,26),qin(3)
 real(r15), intent(out):: qout(3)
 integer i,nboundary
 logical insid

 qout=qin
 do i=1,26
    do while( 2*dot_product(qout,gs(:,i)).gt.dot_product(gs(:,i),gs(:,i)) )
       qout=qout-gs(:,i)   
       write(*,4)'<- i, q,qred=',i,qout !,reduce_g(qout)
    enddo
  ! do while( 2*dot_product(qout,gs(:,i)).lt.-dot_product(gs(:,i),gs(:,i)) )
  !    qout=qout+gs(:,i)   
! !    write(*,4)'-> i, q,qred=',i,qout,reduce_g(qout)
  ! enddo
 enddo
 write(*,4)'final, q,qred=',i,qout !,reduce_g(qout)

 call check_inside_ws(qout,gs,insid,nboundary) 
 if(.not.insid) then
   write(*,*)'BRING_TO_WS_G ERROR: final q=',qout
   stop
 endif

4 format(a,i4,9(1x,f10.5)) 
 end subroutine bring_to_ws_g
!=========================================
 subroutine show_ws_boundary(r1,r2,r3,sx,m,fn,rmax) 
!! for the 26 nearest vectors , generates a fine grid and outputs the points nearest to the WS cell boundary 
!! m is the gridsize, fn is the output filename
 use constants, only : r15
 use geometry, only : length
 use ios, only : ulog
 implicit none
 integer, intent(in) :: m 
 real(r15), intent(in) :: sx(3,26),r1(3),r2(3),r3(3)
 real(r15), intent(out) :: rmax
 character(*), intent(in) :: fn
 real(r15) x(3)
 integer i,j,k,nb,iout,nboundary
 logical insid

 iout=567
 open(iout,file=fn)
 rmax=0; nboundary=0
 do i=-2*m,2*m
 do j=-2*m,2*m
 do k=-2*m,2*m
    x= (i*r1+j*r2+k*r3)/(2d0*m)
    call check_inside_ws(x,sx,insid,nb) 
    if(insid .and. nb.gt.0) then
          rmax=max(rmax,length(x))
          nboundary=nboundary+1
   !      write(*,*)'rmax,l(x)=',rmax,length(x)
          write(iout,3)" Si ", x 
    else
       cycle
    endif
 enddo
 enddo
 enddo
 close(iout)
 write(ulog,*)'SHOW_WS_BOUNDARY generated ', nboundary,' points on the boundary with m=',m
 write(ulog,*)'SHOW_WS_BOUNDARY rmax=',rmax
3 format(a,9(1x,f9.4))
4 format(5i4,9(1x,f9.4))

 end subroutine show_ws_boundary
!=====================================
 subroutine find_sym_op(k1,k2,mat)
!! for two vectors k1,k2 find the symmetry matrix that took k1 to k2: k2=mat*k1
 use constants, only : r15,pi
 use lattice, only : primitivelattice,prim_to_cart,cart2red_g
 use atoms_force_constants, only : op_kmatrix,lattpgcount
 use geometry, only : length
 use ios, only: write_out
 implicit none
 real(r15), intent(in) :: k1(3),k2(3)
 real(r15), intent(out):: mat(3,3)
 integer i
 real(r15) v(3),v2(3)

 write(*,3)'FIND_SYMOP: finding rotation matrix taking ',k1,' to ',k2
 write(*,3)'FIND_SYMOP: taking reduced vector k1 ',cart2red_g(k1),' to k2 ',cart2red_g(k2)
 if(abs(length(k1)-length(k2)).gt.1d-4) then
    write(*,3)' the 2 vectors have different lengths ',length(k1),length(k2)
    stop
 endif
 call write_out(6,' primitive lattice ',primitivelattice)
 call write_out(6,' prim_to_cart/2pi  ',prim_to_cart/(2*pi))
 do i=1,lattpgcount
    mat=matmul(matmul(transpose(primitivelattice),op_kmatrix(:,:,i)),transpose(prim_to_cart))/(2*pi)
!  matmul(transpose(cart_to_prim),matmul(op_kmatrix(:,:,g),transpose(prim_to_cart))))
    call write_out(6,' trying rotation matrix ',op_kmatrix(:,:,i))
    call write_out(6,' nonreduce rotation matrix ',mat)
    write(6,5)i, matmul(mat,k1),k2
    write(6,5)i,cart2red_g(matmul(mat,k1)),cart2red_g(k2)
    if(length(k2-matmul(mat,k1)) .lt. (1+length(k1)+length(k2))*1d-5 ) return
 enddo
 write(*,3)'FIND_SYMOP: no sym op could be found to map ',k1,' to ',k2
! stop
 
3 format(a,3(1x,f10.4),a,3(1x,f10.4))
5 format(i4,9(1x,f10.4))
 end subroutine find_sym_op

!=====================================
 subroutine find_sym_op2(k1,k2,mat)
!! for two reduced vectors k1,k2 find the symmetry matrix that took k1 to k2: k2=mat*k1
 use constants, only : r15,pi
 use lattice, only : primitivelattice,prim_to_cart,cart2red_g
 use atoms_force_constants, only : op_kmatrix,lattpgcount
 use geometry, only : length
 use ios, only: write_out
 implicit none
 real(r15), intent(in) :: k1(3),k2(3)
 real(r15), intent(out):: mat(3,3)
 integer i
 real(r15) v(3),v2(3)

 write(*,3)'FIND_SYMOP2: finding rotation matrix taking ',k1,' to ',k2
 write(*,3)'FIND_SYMOP2: taking reduced vector k1 ',cart2red_g(k1),' to k2 ',cart2red_g(k2)
! if(abs(length(k1)-length(k2)).gt.1d-4) then
!    write(*,3)' the 2 vectors have different lengths ',length(k1),length(k2)
!    stop
! endif
! call write_out(6,' primitive lattice ',primitivelattice)
! call write_out(6,' prim_to_cart/2pi  ',prim_to_cart/(2*pi))
 do i=1,lattpgcount
!   mat=matmul(matmul(transpose(primitivelattice),op_kmatrix(:,:,i)),transpose(prim_to_cart))/(2*pi)
    mat=op_kmatrix(:,:,i)
    call write_out(6,' trying rotation matrix ',op_kmatrix(:,:,i))
!   call write_out(6,' nonreduce rotation matrix ',mat)
    write(6,5)i, matmul(mat,k1),k2
    if(length(k2-matmul(mat,k1)) .lt. (1+length(k1)+length(k2))*1d-5 ) return
 enddo
 write(*,3)'FIND_SYMOP2: no sym op could be found to map ',k1,' to ',k2
! stop
 
3 format(a,3(1x,f10.4),a,3(1x,f10.4))
5 format(i4,9(1x,f10.4))
 end subroutine find_sym_op2
!=====================================
 subroutine generate_grid_noimage(n,grid,sx,newn,map)
!! of the grid vectors inside WS cell defined by sx, it only keeps one on the boundary so
!! that all weights are equal to vol(1-grid)/vol(sx)

 use constants, only : r15,pi
 use lattice, only : primitivelattice,prim_to_cart !,reduce_g
 use geometry, only : length
 use ios, only: write_out,ulog
 implicit none
 integer, intent(in) :: n
 integer, intent(out) :: newn,map(n)  ! map(1:newn) are the indices of the grid to keep
 real(r15), intent(in) :: grid(3,n),sx(3,26)
 integer i,j
 real(r15) gdots, s2  

 newn=0; map=0
 write(ulog,*)'entering generategrid_noimage with n=',n
! call write_out(ulog,'grid points are ',transpose(grid))
! call write_out(ulog,'WS is defined by these 26 vectors',transpose(sx))
 iloop: do i=1,n 
    write(ulog,2)'i,r=',i,grid(:,i)
    jloop: do j=1,26
       s2    = dot_product(sx(:,j),sx(:,j))
       gdots =2d0*dot_product(grid(:,i),sx(:,j))/s2
       if (gdots .le. -1d0) then ! .or. gdots .ge. 1d0) then
          cycle iloop
       endif 
    enddo jloop
! didn't cycle iloop so keep it, it's inside
    newn=newn+1
    map(newn)=i
    write(ulog,3)'GRID_NO_IMAGE: NEW VECTOR ',newn,map(newn),grid(:,i) 
 enddo iloop

2 format(a,i5,9(1x,f10.4))
3 format(a,2i5,9(1x,f10.4))

 end subroutine generate_grid_noimage
!=======================================================================
 subroutine generate_grid_sc(n,rgrid) !,ggrid)
!! generates grid of primitive translation vectors inside the supercell excluding their images
      use geometry
      use lattice
      use ios
      use atoms_force_constants
      use constants , only:r15
      implicit none
      integer, intent(in) :: n
      real(r15), intent(out) :: rgrid(3,n) !,ggrid(3,n)
      real(r15) r0(3)
      integer j,cnt

      rgrid = 0;cnt=0
      do j=1,natom_super_cell
         if (atom_sc(j)%cell%tau .eq. 1 ) then
             cnt=cnt+1 
             if(cnt.eq.1) r0=v2a(atom_sc(j)%equilibrium_pos)
             rgrid(:,cnt)= v2a(atom_sc(j)%equilibrium_pos)-r0
         endif
      enddo 

      if(n.ne.cnt) then
         write(*,*)'GENERATE_GRID_SC: did not count the same number :n,cnt',n,cnt
         stop
      endif

 end subroutine generate_grid_sc
!=======================================================================
! Function to check if a k-point in reduced coordinates is within the wedge
function is_in_wedge(kpoint) result(in_wedge)
  use constants , only:r15
  implicit none
  real(r15), intent(in) :: kpoint(3)
  logical :: in_wedge

  ! Example conditions for a cubic system
  if (kpoint(1) >= 0.0d0 .and. kpoint(1) <= kpoint(2) .and. &
  &   kpoint(2) <= kpoint(3) .and. kpoint(3) <= 0.5d0) then
    in_wedge = .true.
  else
    in_wedge = .false.
  end if
end function is_in_wedge
!============================================================
 subroutine get_components_g0(q,n,i,j,k,g1,g2,g3,inside)
! for a given q-vector, it finds its integer components assuming it was
! created as: q=(i-1)/n1*g1 + (j-1)/n2*g2 + (k-1)/n3*g3; i=1,n1 j=1,n2 k=1,n3
! as the ones created in make_kp_reg with zero shift
! if the integer variable inside=1 then q is inside the primcell
! if q is outside the prim_cell, then inside=0 but produces the i,j,k of its
! image inside
 use geometry
 implicit none
 real(8), intent(in):: q(3)
 type(vector),intent(in):: g1,g2,g3
 type(vector)rr1,rr2,rr3
 integer, intent(in):: n(3)
 integer, intent(out):: i,j,k,inside


 call make_reciprocal_lattice_v(g1,g2,g3,rr1,rr2,rr3)
! i = nint(1+ n(1)* q.dot.r1)
! j = nint(1+ n(2)* q.dot.r2)
! k = nint(1+ n(3)* q.dot.r3)
! if q.dot.r is larger than 1, we want to have its fractional part
 inside = 1

 call comp1(q,rr1,n(1),i,inside)
 call comp1(q,rr2,n(2),j,inside)
 call comp1(q,rr3,n(3),k,inside)

 end subroutine get_components_g0
!============================================================
 subroutine get_components_g2(q,n,i,j,k,g1,g2,g3,inside)
! for a given q-vector, it finds its integer components assuming it was
! created as: q=(i-1-n1/2)/n1*g1 + (j-1-n2/2)/n2*g2 + (k-1-n3/2)/n3*g3; 
! i=1,n1 j=1,n2 k=1,n3
! as the ones created in make_kp_reg with zero shift
! if the integer variable inside=1 then q is inside the primcell, if inside=0 it's outside
! it works even if there's a shift less than 0.5, and if q is oustide the prim_cell
 use geometry
 implicit none
 real(8), intent(in):: q(3)
 type(vector),intent(in):: g1,g2,g3
 integer, intent(in):: n(3)
 integer, intent(out):: i,j,k,inside
 real(8) w(3)

 w = q+0.5d0*(g1+g2+g3)
! write(*,8)'q=',q 
! write(*,8)'w=',w 
 call get_components_g0(w,n,i,j,k,g1,g2,g3,inside)
! write(*,*)'w-components,inside are=',i,j,k ,inside

 end subroutine get_components_g2
!===========================================================
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
 !  call get_components_g0(q,N,i,j,k,g1,g2,g3,inside)
 !  nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
 !  nk= indexg(i,j,k,-N(1)/2,N(1)/2-1,-N(2)/2,N(2)/2-1,-N(3)/2,N(3)/2-1)

! this is for when the k vectors are between -N/2 and N/2-1
  if (mod(N(1),2).eq.0 .and. mod(N(2),2).eq.0 .and. mod(N(3),2).eq.0 ) then
    call get_components_g2(q,N,i,j,k,g1,g2,g3,inside)
!   write(*,*)'for even ni, ijk=',i,j,k
  else   ! do the regulat thing
    call get_components_g0(q,N,i,j,k,g1,g2,g3,inside)
  endif
  nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
! write(*,*)'the corresponding nk=',nk

  if (inside.eq.2) then
     write(*,*)' ERROR: vector ',q,' is not a kpoint'
     stop
  endif

 end subroutine get_k_info
!===========================================================
 subroutine findgrid(r,grid,ngrid,igrid)
!! find  to which grid index does the vector r correspond, igrid=0 means not found
 use geometry !, only : length,myeq
 use constants, only : r15
 implicit none
 integer, intent(in):: ngrid
 integer, intent(out):: igrid
 real(r15), intent(in):: r(3),grid(3,ngrid)
 integer l

 igrid=0
 do l=1,ngrid
    if(length(grid(:,l)-r) .lt. 2d-3) then
       igrid=l
       exit
    endif
 enddo 

 end subroutine findgrid
