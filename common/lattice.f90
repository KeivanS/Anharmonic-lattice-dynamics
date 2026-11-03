 module lattice
 use geometry
 use constants
 implicit none
 type(vector) gs1,gs2,gs3,rs1,rs2,rs3  ! rri=translation vectors of the supercell
 type(vector) r1conv,r2conv,r3conv,g1conv,g2conv,g3conv
 type(vector) r01,r02,r03,g01,g02,g03  ! tr vect of prim cell and its recip spce
 real(r15) volume_r,volume_g,volume_r0,volume_g0,a0_scale,asc_scale
 real(r15) lattice_parameter,latticeparameters(6),primitivelattice(3,3)
 real(r15),dimension(3,3):: conv_to_cart,cart_to_conv,prim_to_conv,conv_to_prim,  &
  &                       prim_to_cart,cart_to_prim ,cart_to_redpr,cart_to_redsc
 real(r15) box(3),boxg(3),density
 real(r15) r0g(3,3)
 real(r15) gws26(3,26),rws26(3,26) ! superlattice shells defining the WS of SL
 real(r15) g0ws26(3,26),r0ws26(3,26) , lgridmax
 integer n_sc(3,3),invn_sc(3,3)  
 integer, allocatable :: meshr3d(:,:),meshg3d(:,:)

  interface make_rg
     module procedure make_rga,make_rgv
  end interface

  interface check_int
     module procedure check3,check_a
  end interface

  interface cart_to_direct
     module procedure cart_to_direct_aa,cart_to_direct_av,cart_to_direct_v
  end interface

  interface direct_to_cart
     module procedure direct_to_cart_aa,direct_to_cart_av,direct_to_cart_v
  end interface

  interface bring_to_cell_c
     module procedure bring_to_cell_cv,bring_to_cell_ca
  end interface

  interface bring_to_cell_d
     module procedure bring_to_cell_dv,bring_to_cell_da
  end interface

 contains

!-------------------------------------------
 function cart2red(x,space) result(red)
!! converts cartesian coordinates of a vector to its reduced coordinates
 use constants, only : r15,pi
 implicit none
 real(r15) , intent(in) :: x(3)
 character(len=1), intent(in) :: space
 real(r15)  red(3)
 if (space.eq.'r' .or. space .eq. 'R') then
     red=matmul(cart_to_prim,x)
 elseif (space.eq.'g' .or. space .eq. 'G') then
     red=matmul(transpose(prim_to_cart),x)/(2*pi)
 else
     write(*,*)'CART2RED: space should either be r or g not ',space
     stop
 endif
   
 end function cart2red
!-------------------------------------------
 function red2cart(x,space) result(cartx)
!! converts reduced coordinates of a real space vector to its cartesian coordinates
 use constants, only : r15,pi
 implicit none
 real(r15) , intent(in) :: x(3)
 character(len=1), intent(in) :: space
 real(r15)  cartx(3)
 if (space.eq.'r' .or. space .eq. 'R') then
     cartx=matmul(prim_to_cart,x)
 elseif (space.eq.'g' .or. space .eq. 'G') then
     cartx=matmul(transpose(cart_to_prim),x)*(2*pi)
 else
     write(*,*)'RED2CART: space should either be r or g not ',space
     stop
 endif
   
 end function red2cart
!-------------------------------------------
 function red2cart_r(x) result(cartx)
!! converts reduced coordinates of a real space vector to its cartesian coordinates
 use constants, only : r15
 implicit none
 real(r15) , intent(in) :: x(3)
 real(r15)  cartx(3)
 cartx=matmul(prim_to_cart,x)
 end function red2cart_r
!-------------------------------------------
 function red2cart_g(q) result(cartq)
!! converts reduced coordinates of a kpoint to its cartesian coordinates
 use constants, only : pi,r15
 implicit none
 real(r15) , intent(in) :: q(3)
 real(r15)  cartq(3)
 cartq=matmul(transpose(cart_to_prim),q)*(2*pi)
 end function red2cart_g
!-------------------------------------------
 function cart2red_g(q) result(redq)
!! converts cartesian coordinates of a kpoint to its reduced coordinates
 use constants, only : pi,r15
 implicit none
 real(r15) , intent(in) :: q(3)
 real(r15)  redq(3)
 redq=matmul(transpose(prim_to_cart),q)/(2*pi)
 end function cart2red_g
!-------------------------------------------
 function cart2red_r(x) result(redx)
 use constants, only : pi,r15
 implicit none
 real(r15) , intent(in) :: x(3)
 real(r15) redx(3)
 redx=matmul(cart_to_prim,x)
 end function cart2red_r
!-------------------------------------------
function fold_ws(q, gs) result(qout)
! folds vector q in the WS cell defined by gs vectors
   use constants, only : r15
   use geometry, only : length
   use params, only : tolerance
   implicit none
   real(r15), intent(in) :: q(3), gs(3,26)
   real(r15) :: qout(3), qtest(3), dist_min, dist
   integer :: i, j, niter, nboundary
   logical :: inside, changed
   
   qout = q
   
   do niter = 1, 100  ! Increase max iterations
      call check_inside_ws(qout, gs, inside, nboundary)
      if (inside) return
      
      ! Find the G that brings q closest to origin
      dist_min = length(qout)
      changed = .false.
      
      do i = 1, 26
         qtest = qout - gs(:,i)
         dist = length(qtest)
         
         if (dist < dist_min - tolerance) then
            qout = qtest
            dist_min = dist
            changed = .true.
         endif
      enddo
      
      ! If no improvement, we're stuck
      if (.not. changed) then
         write(*,*) 'FOLD_WS: Cannot reduce distance further'
         write(*,*) 'q =', qout, ' dist =', dist_min
         exit
      endif
   enddo
   
   if (.not. inside) then
      write(*,*) 'FOLD_WS: Failed after', niter, 'iterations'
      write(*,*) 'q =', qout
      stop
   endif
   
end function fold_ws
!-------------------------------------------
 function fold_ws1(q,gs) result(qout)
!! this function folds the vector q in the Wigner-Seitz cell defined by gs(3,1:26)
 use constants, only : r15
 use geometry, only : v2a,length
 use params, only : tolerance
 implicit none
 real(r15), intent(in):: q(3),gs(3,26)
 real(r15) x,c(3),qin(3),test_vec(3),dist,dmin,qout(3)
 integer i,j,k,m,nboundary,niter
 logical inside,folded

! first bring it in -G/2;G2 by using the reduced coordinates
  qin=q
  do niter=1,20
     call check_inside_ws(qin,gs,inside,nboundary)
     if (inside) then
        qout=qin    
        return
     endif

! try folding by subtracting tnearest face-normal component
     do i=1,26
        x=dot_product(qin,gs(:,i))/dot_product(gs(:,i),gs(:,i))
        if(abs(x) > 0.5_r15) then 
           qin=qin-gs(:,i)*nint(x)
        endif
     enddo

  enddo

  write(*,*)' FOLD_WS : failed to fold vector q=',qin
  stop

 end function fold_ws1
!-------------------------------------------
 function fold_ws_old(q,gs,space) result(qout)
!! for WS cell defined by gs(3,1:26), this folds q in the ws cell of gsi
 use constants, only : r15
 use geometry, only : v2a,length
 use params, only : tolerance
 implicit none
 real(r15), intent(in):: q(3),gs(3,26)
 real(r15) x,c(3),qin(3),test_vec(3),dist,dmin,qout(3)
 character(len=1), intent(in) :: space
 integer i,j,k,m,nboundary
 logical inside,folded

! first bring it in -G/2;G2 by using the reduced coordinates
  qin=q
  do i=1,26
     x=dot_product(qin,gs(:,i))/dot_product(gs(:,i),gs(:,i))
     qin=qin-gs(:,i)*nint(x)
     call check_inside_ws(qin,gs,inside,nboundary)
     if (inside) then
        qout=qin    
        return
     endif
  enddo

  dmin=length(qin)
  do i = 1, 26
     dist = length(qin+gs(:,i))
     if(dist.lt.dmin) then
        dmin=dist
        qout = qin+gs(:,i)
     endif
  enddo

! make sure it is inside
  call check_inside_ws(qout,gs,inside,nboundary)
  if(inside) then
     return
  else
     write(*,*)' final vector is not inside!!',qout
     stop
  endif

   folded = .false.
   k=0 ; dmin=1d8
   do while (.not. folded .and. k.lt.4)
      folded = .true.
      k=k+1
      write(*,*)'FOLD_WS: k=',k
      call check_inside_ws(q,gs,inside,nboundary)
      if (inside) then
         qin=q    
         return
      endif

      do i = 1, 26
         qin = q+gs(:,i)*k
         dist = length(qin) !   )) sqrt(sum(folded_vec**2))

      !  call check_inside_ws(qin,gs,inside,nboundary)
      !  if (inside) then
      !    return
      !  endif
!        folded = .false.
          
          if (dist < dmin) then
            dmin = dist
            qout=qin
            folded = .false.
          end if
      end do
   end do

 return
 
! first find its reduced units
! c=cart2red_g(q)
  !    if(x.gt.0.5001) then
  !      m=nint(x) 
  !      qin=qin-m*gs(:,i)
  
 do i=1,26
    x=dot_product(qin,gs(:,i))/dot_product(gs(:,i),gs(:,i))
    m=nint(x) !=floor(x) !mod(x,1) !m=x-y
    qin=qin-m*gs(:,i)  
 enddo
 
 end function fold_ws_old

!-------------------------------------------
 subroutine make_r0g
! matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
! matmul(r0g,v) gives the 3 reduced coordinates of v in units of the supercell translations 
! coordinates are v(1)*r01+v(2)*r02+v(3)*r03; r0g is the inver of N matrix  R_sc=N.R_0
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
!-------------------------------------------
! subroutine make_rg(x1,x2,x3,q1,q2,q3,n)
!  matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
!  matmul(r0g,v) gives the 3 reduced coordinates of an atom in the supercell if its
!  reduced coordinates in the primitive cell are given by v
! use geometry
! implicit none
! real(r15) x1(3),x2(3),x3(3),q1(3),q2(3),q3(3),n(3,3)

! n(1,1) = x1 .dot. q1
! n(1,2) = x1 .dot. q2
! n(1,3) = x1 .dot. q3
! n(2,1) = x2 .dot. q1
! n(2,2) = x2 .dot. q2
! n(2,3) = x2 .dot. q3
! n(3,1) = x3 .dot. q1
! n(3,2) = x3 .dot. q2
! n(3,3) = x3 .dot. q3

! end subroutine make_rg
! -------------------------------------------
 subroutine check_a(r,a,ier,g1,g2,g3)
!! subroutine to check whether r is an integer multiple of (r01,r02,r03)
!! output ai are the coefficients of its linear combination on this basis
!! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use constants, only : pi
 use ios
 implicit none

 type(vector) :: r,g1,g2,g3
 real(r15) a(3)
 integer ier,n(3)
! logical is_integer

! lengths smaller than tolerance are considered to be zero
!call write_out(ulog,'CHECK: R ',r)
 a(1) = (r .dot. g1)/(2*pi)  ! g's have a 2 pi
 a(2) = (r .dot. g2)/(2*pi)
 a(3) = (r .dot. g3)/(2*pi)
 n=floor(a)
! if(is_integer_a(a)) then
 if(is_integer(a)) then
!  write(ulog,*) ' R is a multiple of r0s'
!  write(ulog,3) ' n1,n2,n3  =',a1,a2,a3
   ier = 0
 else
   ier = 1
!  write(ulog,*) ' R is not a multiple of r0s , check your inputs '
!  write(ulog,3) 'ier=1;r and r.gi/2pi=',r,a
!  stop
 endif
3 format(a,9(1x,f11.6))

 end subroutine check_a
!-------------------------------------------
 subroutine check3(r,a1,a2,a3,ier,g1,g2,g3)
!! subroutine to check whether r is an integer multiple of (r01,r02,r03)
!! output ai are the coefficients of its linear combination on this basis
!! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use constants, only : pi
 use ios
 implicit none

 type(vector) :: r,g1,g2,g3
 real(r15) a1,a2,a3
 integer ier

! lengths smaller than tolerance are considered to be zero
! call write_out(ulog,'CHECK: R ',r)
 a1 = (r .dot. g1)/(2*pi)  ! g's have a 2 pi
 a2 = (r .dot. g2)/(2*pi)
 a3 = (r .dot. g3)/(2*pi)
 if(abs(a1-nint(a1)).ge.tolerance .or.abs(a2-nint(a2)).ge.tolerance   &
&   .or.abs(a3-nint(a3)).ge.tolerance ) then
   ier = 1
!  write(ulog,*) ' R is not a multiple of r0s , check your inputs '
!  write(ulog,3) 'ier=1;r and r.gi/2pi=',r,a1,a2,a3
!  stop
 else
!  write(ulog,*) ' R is a multiple of r0s'
!  write(ulog,3) ' n1,n2,n3  =',a1,a2,a3
   ier = 0
 endif
3 format(a,9(1x,f11.6))

 end subroutine check3
!-------------------------------------------
 subroutine check_d(r,a1,a2,a3,ier)
! subroutine to check whether r is an integer vector
! output ai are the coefficients of its linear combination on this basis
! ier=0 means r is an integer multiple of the basis
 use geometry
 use params
 use ios
 implicit none
 type(vector) :: r
 real(r15) a1,a2,a3
 integer ier

!call write_out(ulog,'CHECK: R ',r)
 a1 = r%x
 a2 = r%y
 a3 = r%z
 if(abs(a1-nint(a1)).ge.tolerance .or.abs(a2-nint(a2)).ge. tolerance   &
&   .or.abs(a3-nint(a3)).ge.tolerance ) then
   ier = 1
   write(ulog,3) 'ier=1; r is not integer ',a1,a2,a3
 else
   ier = 0
 endif
3 format(a,9(1x,g11.4))

 end subroutine check_d
!-----------------------------------------
 subroutine transform_input_structure
!! takes input variables latticeparameters(6),primitivelattice(3,3) to construct the primitive cell
 use ios
! use constants, only : pi
! use atoms_force_constants, only : xmatinv
 implicit none
 external xmatinv
 real(r15) latp(6),d2r,ca,cb,cg,sg
 integer ier

!k1
!      d2r = 4d0*datan(1d0)/180d0 ! 3.1415926535897932384626/180d0
!k1
      conv_to_cart=0
      prim_to_conv=primitivelattice
!     call write_out(6,'Prim_to_conv=prim vectors in reduced units of conv in lines ',prim_to_conv)
!     call write_out(ulog,'Prim_to_conv=prim vectors in reduced units of conv in lines ',prim_to_conv)
      call xmatinv(3,prim_to_conv,conv_to_prim,ier)
!     call write_out(6,'conv_to_prim has primitive reciprocal vectors in units of recip conv ',conv_to_prim)
!     call write_out(ulog,'conv_to_prim has primitive reciprocal vectors in units of recip conv ',conv_to_prim)
      if(ier.ne.0) then
         write(*,*)'prim_to_conv does not have an inverse!'
         stop
      endif
      latp=latticeparameters ! for short
!      latp(4:6) = latp(4:6)*d2r   ! convert angles from degree to radian
! cos and sin of the angles, which are in degrees. dcosd/dsind are a GNU extension
! that some compilers do not have (Rivanna's gfortran, 2026-10-30); this is standard.
      d2r = pi/180d0
      ca = cos(latp(4)*d2r)
      cb = cos(latp(5)*d2r)
      cg = cos(latp(6)*d2r)
      sg = sin(latp(6)*d2r)
! exact at the right angles, as dcosd/dsind were (cos(90*pi/180) is 6d-17, not 0)
      if(abs(ca).lt.1d-15) ca=0
      if(abs(cb).lt.1d-15) cb=0
      if(abs(cg).lt.1d-15) cg=0
      if(abs(sg-1).lt.1d-15) sg=1
      conv_to_cart(1,1)=latp(1)
      conv_to_cart(1,2)=latp(2)*cg
      conv_to_cart(2,2)=latp(2)*sg
      conv_to_cart(1,3)=latp(3)*cb
      conv_to_cart(2,3)=latp(3)*(ca   &
     &     -cg*cb)  /sg
      conv_to_cart(3,3)=sqrt(latp(3)**2-conv_to_cart(1,3)**2   &
     &     -conv_to_cart(2,3)**2)
      write(*,*)'conv_to_cart has conventional vectors on its columns='
      write(*,*)'prim_to_cart has primitive vectors on its columns'
      write(*,*)'cart_to_conv has conventional reciprocal vectors on its lines (no 2pi)'
      write(*,*)'cart_to_prim has primitive reciprocal vectors on its lines (no 2pi) '
      call write_out(ulog,' Conventional cell(conv_to_cart) on columns ',conv_to_cart)
! go to reciprocal space by matrix inversion
      call xmatinv(3,conv_to_cart,cart_to_conv,ier)
      call write_out(ulog,' Conv. reciprocal lattice on lines  (no 2pi)',cart_to_conv)

! the usual way
      r1conv=conv_to_cart(:,1)
      r2conv=conv_to_cart(:,2)
      r3conv=conv_to_cart(:,3)
! this has an extra factor of 2pi wrt cart_to_conv
      call make_reciprocal_lattice_2pi(r1conv,r2conv,r3conv,g1conv,g2conv,g3conv)
      call write_out(ulog,' Conventional reciprocal(1) with 2pi',g1conv)
      call write_out(ulog,' Conventional reciprocal(2) with 2pi',g2conv)
      call write_out(ulog,' Conventional reciprocal(3) with 2pi',g3conv)
!     call write_out(ulog,' Conv recip on lines/2pi ',cart_to_conv)
! prim_to_cart(i,j) is the ith cartesian coordinate of the jth primitive translation vector
      prim_to_cart=matmul(conv_to_cart,prim_to_conv)
      call write_out(ulog,' Primitive translations: prim_to_cart ( on columns ) ',prim_to_cart)
! i.e. multiplied by primitive reduced coordinates will give you cartesian coordinates
      r01 = prim_to_cart(:,1)
      r02 = prim_to_cart(:,2)
      r03 = prim_to_cart(:,3)

      call make_reciprocal_lattice_2pi(r01,r02,r03,g01,g02,g03)

      call xmatinv(3,prim_to_cart,cart_to_prim,ier)
      if(ier.ne.0)then
        write(6,*)'Error in transform_input_structure: primitive_lattice '   &
     &       //'is singular'
        stop
      endif
      call write_out(ulog,' Primitive reciprocal: cart_to_prim (no 2pi, on lines)  ',cart_to_prim)

 end subroutine transform_input_structure
!-----------------------------------------
 subroutine make_unitcell
!! lattice parameters and atompos0 read from structure.params, it
!! constructs the primitive vectors r0i and reciprocal vectors g0i
!! initializes the object atom0%equilibrium_pos, mass, name, type and charge
!! finds the symmetry operations of the unit cell
!! calculates number and position of the atoms within rcut(2)=15 of the primitive cell
 use ios
 use params
 use geometry
 use atoms_force_constants
 use constants
 implicit none
 integer i,mx1
 real(r15) particle_density,volcut,r1(3),r2(3),r3(3),volmax,rshells(3,26),rc1
 external unitcell,symmetry_init

! generate r0i and g0i from latticeparameters, primitivelattice
  call transform_input_structure
  call calculate_volume(r01,r02,r03,volume_r0)
  call calculate_volume(g01,g02,g03,volume_g0)

  write(ulog,*)' PRIMITIVE CELL *************************'
  call write_out(ulog,' r01',r01)
  call write_out(ulog,' r02',r02)
  call write_out(ulog,' r03',r03)
  call write_out(ulog,' Volume of primitive cell ',volume_r0)
  call write_out(ulog,' g01',g01)
  call write_out(ulog,' g02',g02)
  call write_out(ulog,' g03',g03)
  call write_out(ulog,' Volume of reciprocal primitive cell ',volume_g0)
  call write_out(6,' primitive lattice (lines)',transpose(prim_to_cart))
  call write_out(6,' reciprocal lattice (lines, no 2pi) ',cart_to_prim)


  box(1) = length(r01)
  box(2) = length(r02)
  box(3) = length(r03)
  boxg(1) = length(g01)
  boxg(2) = length(g02)
  boxg(3) = length(g03)
  write(ulog,3)' Lengths of r0i = ',box
  write(ulog,3)' Lengths of g0i = ',boxg
! write(ulog,3)' volume_r0,volume_g0 = ',volume_r0,volume_g0


  write(ulog,*)' Atomic positions as read; convcoord=',convcoord,' (0: conventional, else primitive units)'
  do i=1,natom_prim_cell
     write(ulog,8)atom0(i)%name,atom0(i)%at_type,atom0(i)%tau,atompos0(:,i),atom0(i)%mass
  enddo

! convert to cartesian
  write(ulog,*)' reduced Atomic positions in conventional, then primitive units, mass ======== '
  do i=1,natom_prim_cell
     if(convcoord.eq.0) then ! atompos0 is in units of the conventional vectors (line 1)
        atompos0(:,i)=matmul(conv_to_cart,atompos0(:,i))
     else                    ! atompos0 is in units of the primitive vectors (line 2)
        atompos0(:,i)=matmul(conv_to_cart,matmul(prim_to_conv,atompos0(:,i)))
     endif
     write(ulog,8)' '//atom0(i)%name,atom0(i)%at_type,atom0(i)%tau,  &
&    matmul(cart_to_conv,atompos0(:,i)),matmul(cart_to_prim,atompos0(:,i)),atom0(i)%mass

! bring to [0 r0i[
     call unitcell(cart_to_prim,prim_to_cart,atompos0(:,i), atompos0(:,i))
     atom0(i)%equilibrium_pos = atompos0(:,i)
  enddo

  write(ulog,*)' Cartesian Atomic positions, mass ====================== '
  do i=1,natom_prim_cell
     write(ulog,8)'After unitcell '//atom0(i)%name,atom0(i)%at_type,atom0(i)%tau,atom0(i)%equilibrium_pos,atom0(i)%mass
  enddo

  density=sum(atom0(:)%mass)/volume_r0*uma*1d27
  particle_density=natom_prim_cell/volume_r0  
  write(ulog,*)'       density (uma/Ang^3)=',sum(atom0(:)%mass)/volume_r0
  write(ulog,*)'       density (  g/ cm^3)=',density
  write(ulog,*)' particle density (1/Ang^3)=',particle_density
  write(ulog,*)' COORDINATES read successfully '
  write(   *,*)' COORDINATES read successfully '


  call find_cutoff_from_largest_SC (volmax,rcutoff,r1,r2,r3,rshells)

! nshells should be within this cutoff
  write(ulog,*)' the suggested value of rcutoff from SC is ',rcutoff
! icutoff is an option to increase it if need be
  rcutoff=max(rcutoff, icutoff*volume_r0**0.333)
  write(ulog,*)' rcutoff updated from icutoff is ',rcutoff
! use this volume for atompos; the largest supercell will be fully included in this volume
  volcut=4*pi/3d0*(rcutoff*1.3)**3   
  maxatoms = nint(volcut * particle_density)  
! this seems reasonably large, but can be increased if more is needed
  maxneighbors=100
  maxshells = icutoff   ! let the user decide it

! FCs will exist at the most up to maxshells; it's ok if maxshell is too large
! if (fc2range .eq.0 ) then
!    maxshells=nint(rcutoff / min(length(r01),length(r02),length(r03)))+1
! else
!    maxshells=max(maxval(nshells)+2,icutoff )  ! go with what the user prescribed
! endif


!  mx1=nint(4*pi/3*icutoff**3)
!  write(ulog,*)'reduced scutoff=',icutoff,' suggests maxatoms=',mx1
!! Convert to actual distances
!  rc1 = 1d0/particle_density**0.333333 * icutoff  ! based on average bond length
! rcutoff=max(rcutoff,rc1)  ! to be able to change using the default if too small 
! write(ulog,*)' if you need a larger value, increase the value of icutoff in default.params from 10'  



!  everything is based on maxneighbors, maxatoms and maxshells will be updated accordingly
! maxneighbors is the largest number of neighborshells; largest in a glass or 
! huge unit cell
! maxneighbors=min(3*maxshells,nint(((2*maxshells+1)**3)/48d0))

  write(ulog,*)'Based on this rcutoff, volcut=4pi/3*(1.3rcut)^3, maxatoms=volcut*density'
  write(ulog,*)'the default  maxatom=', maxatoms
  write(ulog,*)'maxshells set to    =', maxshells
  write(ulog,*)'Start   maxneighbors=', maxneighbors
  write(ulog,*)'*****************************************************************'

  imaxatm=1
     if(allocated(atompos   )) deallocate(atompos   )
     if(allocated(iatomcell )) deallocate(iatomcell )
     if(allocated(iatomcell0)) deallocate(iatomcell0)
     allocate(atompos(3,maxatoms), iatomcell(3,maxatoms),iatomcell0(maxatoms))

  imaxnei=1

  do while(imaxnei.eq.1)

! find symmetry operations; the above allocation is for this subroutine
! requires maxneighbors to generate atoms within these shells
     call symmetry_init(natom_prim_cell,atom_type,atompos0)
  enddo
!    ----------------        DO not want to increase maxshells
! do while(imaxshl.eq.1)
!    maxshells=maxshells+5
!    write(*,*)' updating maxshells to ',maxshells
!    call force_constants_init(natom_prim_cell,atom_type,atompos0)
! enddo
  write(*   ,*)' maxatoms, natoms=',maxatoms,natoms
  write(ulog,*)' maxatoms, natoms=',maxatoms,natoms

  call test_ntau
! cat_to_prim had reciprocal lattice vectors on its lines
! call write_out(ulog,'reciprocal lattice vectors ',cart_to_prim)


3 format(a,9(2x,g13.6))
4 format(i6,3(2x,g12.5),4(3x,i4))
6 format(a,(2x,g12.5),a,4(3x,i6))
8 format(a,2i3,9(1x,g13.6))

 end subroutine make_unitcell
!-----------------------------------------
! subroutine get_components(q,n,i,j,k )  !,g1,g2,g3)
!! for a given q-vector, it finds its integer components assuming it was
!! created as: nk=0 do i=0,N(1)-1 ; do j=0,N(2)-1; do k=0,N(3)-1; nk=nk+1 ;q=i*G1/N1+j*G2/N2+k*G3/N3
!! as the ones created in make_kp_reg with zero shift
!! it works even if there's a shift less than 0.5
! implicit none
! real(r15) q(3) !,g1(3),g2(3),g3(3),r1(3),r2(3),r3(3)
!! type(vector) g1,g2,g3,r1,r2,r3
! integer i,j,k,n(3)
!
!! call make_reciprocal_lattice(g1,g2,g3,r1,r2,r3)
! i = nint(1+ n(1)* (q.dot.r1)/2/pi )
! j = nint(1+ n(2)* (q.dot.r2)/2/pi )
! k = nint(1+ n(3)* (q.dot.r3)/2/pi )
!
! end subroutine get_components
!-----------------------------------------
 subroutine make_rga(x1,x2,x3,q1,q2,q3,n)   ! n(i,j)=xi.qj
! matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
! matmul(r0g,v) gives the 3 reduced coordinates of an atom in the supercell if its
! reduced coordinates in the primitive cell are given by v
 use geometry
 implicit none
 real(r15), intent(in ) :: x1(3),x2(3),x3(3),q1(3),q2(3),q3(3)
 real(r15), intent(out) :: n(3,3)

 n(1,1) = x1 .dot. q1
 n(1,2) = x1 .dot. q2
 n(1,3) = x1 .dot. q3
 n(2,1) = x2 .dot. q1
 n(2,2) = x2 .dot. q2
 n(2,3) = x2 .dot. q3
 n(3,1) = x3 .dot. q1
 n(3,2) = x3 .dot. q2
 n(3,3) = x3 .dot. q3

 end subroutine make_rga
!-----------------------------------------
 subroutine make_rgv(x1,x2,x3,q1,q2,q3,n)   ! n(i,j)=xi.qj
! matmul(r0g,n) gives the 3 reduced coordinates of the primcell of index n in units of the supercell ri's
! matmul(r0g,v) gives the 3 reduced coordinates of an atom in the supercell if its
! reduced coordinates in the primitive cell are given by v
 use geometry
 implicit none
 type(vector), intent(in) :: x1,x2,x3,q1,q2,q3
 real(r15), intent(out) :: n(3,3)

 n(1,1) = x1 .dot. q1
 n(1,2) = x1 .dot. q2
 n(1,3) = x1 .dot. q3
 n(2,1) = x2 .dot. q1
 n(2,2) = x2 .dot. q2
 n(2,3) = x2 .dot. q3
 n(3,1) = x3 .dot. q1
 n(3,2) = x3 .dot. q2
 n(3,3) = x3 .dot. q3

 end subroutine make_rgv
!-----------------------------------------
 subroutine set_neighbor_list
!! Generates a set of neighbors to the primitive cell up to
!! and calculates the corresponding number of neighborshells
! this is needed in order to know what force constants to associate
! to a pair ij, or ijk or ijkl
! if a neighbor j is known then one can get the vect(i,j) between
! their equilibrium positions
 use atoms_force_constants !, only : natom_prim_cell, atom0,maxatoms
 use params
 use ios
 use geometry
 implicit none
 integer i0,j,jj,s,l,k,shel_count,minat,cnt,nbs,mxshell
 integer nb(0:maxneighbors),msort(maxatoms)  ! because dist is of size natoms
! real(r15) dist(maxatoms),d_old,d_min(0:maxneighbors),rmax,dij,d_new
 real(r15) dist(maxatoms),d_old,d_min(0:maxatoms),rmax,dij,d_new

 call write_out(ulog,'SET_NEIGHBOR_LIST: entering with maxneighbors ',maxneighbors)
 call write_out(ulog,'SET_NEIGHBOR_LIST: entering with maxshells    ',maxshells )
 call write_out(ulog,'SET_NEIGHBOR_LIST: entering with natoms       ',natoms)
 write(ulog,*)' the neighbor list will stop at maxshells '
! allocate( atom0(natom_prim_cell)%shells(0:mxshell)) !%neighbors(maxneighbors) )
! allocate( atom0(natom_prim_cell)%shells(0:mxshell)%neighbors(maxneighbors) )
 rmax = 0  ! is the largest radius among the natoms neighbors
 mxshell = 0

 i0loop: do i0=1,natom_prim_cell
    dist = 1d10
! calculate pair i0-j distances, and sort them
    cnt=0
    do j=1,natoms
       dist(j)=length(atompos(:,i0)-atompos(:,j))
    enddo
    call sort(natoms,dist,msort,maxatoms)

! count the corresponding # of shells and nb, # of atoms within that shell, for each i0
    shel_count = -1
    d_old = 1e4
    d_min=0
    nb=1
    minat=min(natoms,350)  ! about 7 shells at the most: 7^3=343
    write(ulog,*)' ========================================================='
    write(ulog,*)' first ',minat,' neighborlist of atom i0 =',i0

    jloop: do j=1,natoms !minat
       d_new = dist(msort(j))  ! start in order of increasing distances from d=0

       if ( abs(d_new-d_old).gt.tolerance ) then   ! new shell
          if(shel_count.eq.maxneighbors) exit jloop
          shel_count = shel_count+1
          nb(shel_count)=1
          d_old = d_new
          d_min(shel_count)=d_new
 !        write(*,4)' NEW shell , dij= ',shel_count,d_new
 !        if ( shel_count .gt. maxneighbors ) then
          if ( d_new .gt. rcutoff ) then
             write(   *,*) shel_count,' shells completed up to distance ',d_new
             write(   *,*)'It  corresponds to ',msort(j),' th nearest atom'
       !     write(   *,*)'The largest # of shells is set to ',shel_count-1 ,' which should be ',maxneighbors
             write(ulog,*) shel_count,' shells completed up to distance ',d_new
             write(ulog,*)'It  corresponds to ',msort(j),' th nearest atom'
!             write(ulog,*)'The largest # of shells is set to ',shel_count-1 ,' which should be ',maxneighbors
             exit jloop
          endif

          iatomneighbor(i0,msort(j))=shel_count

! after maxshell iatomneighbor values are not reliable ! FIXED !

       else  ! same shell
          nb(shel_count)=nb(shel_count)+1  ! counts how many atoms within shell=shel_count
          if(j.lt.90) write(*,9)'shell#, nb#, j, dij= ',shel_count,nb(shel_count),j,d_new
          iatomneighbor(i0,msort(j))=shel_count
       endif
    enddo jloop

!   shel_count=shel_count-1  ! why is that necessary?
    atom0(i0)%nshells=shel_count
    write(*,*)'allocating atom0(',i0,')%shells(:) to size ',shel_count
!    if(allocated(atom0(i0)%shells)) deallocate(atom0(i0)%shells)
    allocate(atom0(i0)%shells(0:shel_count))  ! %neighbors dimension is fixed to 296
    cnt=0  ! counts atoms from nearest to farthest, shell by shell
    do s=0,shel_count
       atom0(i0)%shells(s)%no_of_neighbors = nb(s)
       atom0(i0)%shells(s)%radius = d_min(s)
       do l=1,nb(s)  ! l is the index of the neighbor of i0 in that shell

          cnt=cnt+1  ! counts the total # of atoms up to this shell
          jl:do j=1,natoms
       !     if(cnt.eq.msort(j)) then
             jj=msort(j)
             if(d_min(s) .myeq. length(atompos(:,jj)-atompos(:,i0))) then
! make sure it is not already assigned
               do k=1,l-1
                  if( jj.eq.atom0(i0)%shells(s)%neighbors(k)%atomposindx) cycle jl
               enddo

               atom0(i0)%shells(s)%neighbors(l)%tau = iatomcell0(jj)
               atom0(i0)%shells(s)%neighbors(l)%n   = iatomcell(:,jj)
               atom0(i0)%shells(s)%neighbors(l)%atomposindx = jj
      !        iatomneighbor(i0,msort(j))= s ! atom0(i0)%shells(s)%neighbors(l)
!    if(verbose) then
!        write(ulog,*)' ------------------------------------'
    write(ulog,5)'shl,nb#,atomj,tauj,nj,shel#,dij,Rj =',s,l,jj,iatomcell0(jj),   &
    &             iatomcell(:,jj),iatomneighbor(i0,jj),d_min(s),atompos(:,jj)
!    endif
               exit jl
             endif
          enddo jl
       enddo
!      nbs= atom0(i0)%shells(s)%no_of_neighbors
!      write(ulog,*)atom0(i0)%shells(s)%neighbors(1:nbs)
    enddo

! mxshell is the largest nearest neighbor shell number, and will replace initially guessed maxneighbors
    if(mxshell.le.shel_count) then
        mxshell=shel_count
        write(*   ,*)'i0,shel_count,mxshell=',i0,shel_count,mxshell
        write(ulog,*)'i0,shel_count,mxshell=',i0,shel_count,mxshell
    endif

  if (verbose) then
    do s = 0 , shel_count
       nbs= atom0(i0)%shells(s)%no_of_neighbors
       write(ulog,7)'shell#, neigh#, rij/rij(1),rij=',s,atom0(i0)%shells(s)%no_of_neighbors   &
            ,atom0(i0)%shells(s)%radius/atom0(i0)%shells(1)%radius,atom0(i0)%shells(s)%radius
       write(ulog,*)'    nb#,jatom,tau,     n    '
       do l=1,nbs
          write(ulog,8)' --  ',l,atom0(i0)%shells(s)%neighbors(l)%atomposindx  &
&                ,atom0(i0)%shells(s)%neighbors(l)%tau  &
&                ,(atom0(i0)%shells(s)%neighbors(l)%n(j),j=1,3)
       enddo
    enddo
  endif

! also initialize atom0%equilibrium_pos
    atom0(i0)%equilibrium_pos%x = atompos(1,i0)
    atom0(i0)%equilibrium_pos%y = atompos(2,i0)
    atom0(i0)%equilibrium_pos%z = atompos(3,i0)

 enddo i0loop

 write(ulog,*)'SET_NEIGHBOR_LIST: largest#of neighbor shells=', mxshell
 maxneighbors=mxshell
 write(ulog,*)'# Maxneighbors updated to ',maxneighbors


2 format(a,8(2x,i4))
3 format(a,' (',i4,2(1x,i4),1x,')')
4 format(a,1x,i4,' (',f8.4,2(1x,f8.4),1x,')')
5 format(a,4(1x,i4),' (',i2,2(1x,i2),')',1x,i4,9(1x,f9.4))
6 format(a,1x,f10.7,2x,9(1x,f9.6))
7 format(a,1x,i3,1x,i4,2(1x,f11.5),2x,999('[',i5,',',i2,'(',3(i2),')]'))
8 format(a,i3,i5,1x,i2,' (',i2,2(1x,i2),')')
9 format(a,3i4,3x,9(1x,f9.3))

 end subroutine set_neighbor_list
!-------------------------------------
 subroutine cart_to_direct_v(v,w,g01,g02,g03)
! takes cart coordinates and returns cart coordinates within the supercell
 use geometry
 use constants, only : pi,r15
 implicit none
 type(vector), intent(in) :: v
 type(vector), intent(in) :: g01,g02,g03
 type(vector), intent(out) :: w

 w%x = (v .dot. g01)/(2*pi)
 w%y = (v .dot. g02)/(2*pi)
 w%z = (v .dot. g03)/(2*pi)

 end subroutine cart_to_direct_v
!-----------------------------------
 subroutine cart_to_direct_av(a,w,g01,g02,g03)
! takes cart coordinates and returns cart coordinates within the supercell
 use geometry
 use constants, only : r15
 implicit none
 type(vector), intent(out) :: w
 type(vector), intent(in) :: g01,g02,g03
 real(r15), intent(in) ::  a(3)
 type(vector) v

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call cart_to_direct_v(v,w,g01,g02,g03)

 end subroutine cart_to_direct_av
!-----------------------------------
 subroutine cart_to_direct_aa(a,b,g01,g02,g03)
! takes cart coordinates and returns cart coordinates within the supercell
 use geometry
 use constants, only : r15
 implicit none
 real(r15), intent(in) ::  a(3)
 type(vector), intent(in) :: g01,g02,g03
 real(r15), intent(out) ::  b(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call cart_to_direct_v(v,w,g01,g02,g03)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end subroutine cart_to_direct_aa
!-----------------------------------
 subroutine direct_to_cart_v(v,w,r01,r02,r03)
! takes direct coordinates and returns cart coordinates
 use geometry
 implicit none
 type(vector), intent(in) :: v
 type(vector), intent(in) :: r01,r02,r03
 type(vector), intent(out) :: w

 w%x = v%x*r01%x + v%y*r02%x + v%z*r03%x
 w%y = v%x*r01%y + v%y*r02%y + v%z*r03%y
 w%z = v%x*r01%z + v%y*r02%z + v%z*r03%z

 end subroutine direct_to_cart_v
!-----------------------------------
 subroutine direct_to_cart_av(a,w,r01,r02,r03)
! takes direct coordinates and returns cart coordinates
 use geometry
 use constants, only : r15
 implicit none
 real(r15), intent(in) ::  a(3)
 type(vector), intent(in) :: r01,r02,r03
 type(vector), intent(out) :: w
 type(vector) v

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call direct_to_cart_v(v,w,r01,r02,r03)

 end subroutine direct_to_cart_av
!-----------------------------------
 subroutine direct_to_cart_aa(a,b,r01,r02,r03)
! takes direct coordinates and returns cart coordinates
 use geometry
 use constants, only : r15
 implicit none
 real(r15), intent(in) ::  a(3)
 type(vector), intent(in) :: r01,r02,r03
 real(r15), intent(out) ::  b(3)
 type(vector) v,w

 v%x = a(1) ; v%y = a(2) ; v%z = a(3)
 call direct_to_cart_v(v,w,r01,r02,r03)
 b(1) = w%x ; b(2) = w%y ; b(3) = w%z

 end subroutine direct_to_cart_aa
!------------------------------------
 subroutine dir2cart_g(q,k)
! takes q in direct coordinates and outputs k in cartesian
 use geometry
! use lattice
 use constants, only : r15
 real(r15), intent(in) ::  q(3)
 real(r15), intent(out) ::  k(3)

 k(:) = q(1)*g01 + q(2)*g02 + q(3)*g03

 end subroutine dir2cart_g
!------------------------------------
 subroutine dirconv2cart(q,k)
! takes q in direct coordinates of the conventional cell and outputs k in cartesian
 use geometry
! use lattice
 use ios
 use constants, only : r15
 real(r15), intent(in) ::  q(3)
 real(r15), intent(out) ::  k(3)

 k(:) = q(1)*g1conv + q(2)*g2conv + q(3)*g3conv

! k2=matmul(cart_to_prim,q)
! call write_out(ulog,'k1=',k1)
! call write_out(ulog,'k2=',k2)

 end subroutine dirconv2cart
!-----------------------------------------
 function bring_to_cell_dv(v) result(w)
! takes direct coordinates and returns direct coordinates
 implicit none
 type(vector) v,w
 real(r15) a1,a2,a3

 a1 = v%x ; a2 = v%y ; a3 = v%z
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
 w%x = a1 ; w%y = a2 ; w%z = a3

 end function bring_to_cell_dv
!-----------------------------------------
 function bring_to_cell_da(a) result(w)
! takes direct coordinates and returns direct coordinates
 implicit none
 type(vector) w
 real(r15) a(3)

 a(1) = a(1) - floor(a(1))
 a(2) = a(2) - floor(a(2))
 a(3) = a(3) - floor(a(3))
 w%x = a(1) ; w%y = a(2) ; w%z = a(3)

 end function bring_to_cell_da
!-----------------------------------------
 subroutine bring_to_cell_ca(a,r1,r2,r3,g1,g2,g3,b)
! takes cart coordinates and returns cart coordinates within the supercell
! assumes r_i.g_j=delta_ij (no factor of 2pi included)
 implicit none
 real(r15), intent(in) :: a(3)
 real(r15), intent(out) :: b(3)
 type(vector), intent(in) :: g1,g2,g3,r1,r2,r3
 type(vector) v,w
 real(r15) a1,a2,a3

 v%x = a(1);v%y = a(2);v%z = a(3)
! get direct coordinates first
 a1 = v .dot. g1
 a2 = v .dot. g2
 a3 = v .dot. g3
! bring into cell: between 0 and cell
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
! convert to cartesian coordinates
 w = a1*r1+a2*r2+a3*r3
 b(1)=w%x ; b(2)=w%y ; b(3)=w%z

 end subroutine bring_to_cell_ca
!-----------------------------------
 subroutine bring_to_cell_cv(v,r1,r2,r3,g1,g2,g3,w)
! takes cart coordinates and returns cart coordinates within the supercell
! assumes r_i.g_j=delta_ij (no factor of 2pi included)
 implicit none
 type(vector), intent(in) :: v,g1,g2,g3,r1,r2,r3
 type(vector), intent(out) :: w
 real(r15) a1,a2,a3

! get direct coordinates first
 a1 = v .dot. g1
 a2 = v .dot. g2
 a3 = v .dot. g3
!print*,'v=',v
!print*,'g1,g2,g3=',g1,g2,g3
!print*,'a1,a2,a3=',a1,a2,a3
! bring into cell: between 0 and cell
 a1 = a1 - floor(a1)
 a2 = a2 - floor(a2)
 a3 = a3 - floor(a3)
!print*,'NEW a1,a2,a3=',a1,a2,a3
! convert to cartesian coordinates
 w = a1*r1+a2*r2+a3*r3

 end subroutine bring_to_cell_cv
!-----------------------------------------

 end module lattice
!============================================================
