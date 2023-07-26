!******************************************************************************
! initialize data
      subroutine force_constants_init(natoms_in,iatomtype,atompos_in)
!! given atomic positions in the primitive cell and the translation vectors,
!! it calculates the symmetry operations matrices of the space group

! arguments:
!     lattparams(i) (input), lattice parameters a,b,c,alpha,beta,gamma
!     primlatt(i,j) (input), ith dimensionless coordinate of jth basis vector
!          of the primitive lattice
!     natoms_in (input), number of atoms in the primitive unit cell
!     iatomtype(i) (input), type of ith atom, numbered 1,2,3,etc.
!     atompos_in(j,i) (input), jth dimensionless coordinate of the ith atom
! also uses nsmax for shell looping and rcut and maxneighbors and maxatoms

      use lattice
      use atoms_force_constants !force_constants_module
      use ios !, only : ulog
      use params , only : rcut,nsmax
      implicit none

      integer, intent(in) :: natoms_in
      real(8), intent(in) :: atompos_in(3,natoms_in)
      integer, intent(in) :: iatomtype(natoms_in)
      integer i,j,k,m,n,i1,i2,i3,j1,j2,j3,nd2save,ncmp,nshell,   &
     &     ier,ifractcount,iatom,   &
     &     iatom2,itype2,iatom3,ipntop,   &
     &     itype,isg,iatom0,iop,icell(3),   &
     &     iop_matrix(3,3,48) , natneighbor
      logical foundit,foundone
! maxneighbors= largest # of neighbor shells <  50 usually unless very low symmetry
      real(8) d2save(maxneighbors),r(3),d2, tempi(3,3),latp(6),   & !d2r,   &
     &     atomposconv(3,natoms_in),fract(3),v(3),v2(3),temp(3,3),junk
!    &     lattparams(6),primlatt(3,3),
!    &     conv_to_cart(3,3),prim_to_cart(3,3),cart_to_prim(3,3),   &
!    &     prim_to_conv(3,3),conv_to_prim(3,3), &


!-------------------------------------------------------------------------------
! some preliminary stuff
! copy input arguments into global variables
      natom_prim_cell=natoms_in
      atomposconv(1:3,1:natom_prim_cell)=atompos_in(1:3,1:natom_prim_cell)
      call write_out(6,'atompos_in(primcell)= ',atompos_in)

! allocate memory
      if(allocated(iatomneighbor)) deallocate(iatomneighbor)
      if(allocated(atomopfract))   deallocate(atomopfract)
      if(allocated(iatomop))       deallocate(iatomop)
      allocate(iatomneighbor(natom_prim_cell,maxatoms))
      allocate(iatomop(natom_prim_cell,natom_prim_cell), &
 &       atomopfract(3,natom_prim_cell,natom_prim_cell))
      iatomneighbor=maxneighbors+1
      iatomop=0


 3 format(3(1x,f10.5))

! get primitive lattice parameters in cartesian coordinates
!     call xmatmlt(conv_to_cart,prim_to_conv,prim_to_cart,3,3,3,3,3,3)

      call xmatinv(prim_to_cart,cart_to_prim,ier)
      if(ier.ne.0)then
        write(6,*)'Error in force_constants_init: primitive_lattice '   &
     &       //'is singular'
        stop
      endif
      write(*,*)'cart_to_prim has reciprocal vectors to primitive in its lines'
      write(*,3) cart_to_prim(1,:)
      write(*,3) cart_to_prim(2,:)
      write(*,3) cart_to_prim(3,:)

! added by k1 --------
! prim_to_cart(i,j) is the ith cartesian coordinate of the jth translation
! vector of the primitive lattice.
!     call write_out(ulog,'translations of the prim_cell ',3,3,prim_to_cart)
!     r01 = prim_to_cart(:,1)
!     r02 = prim_to_cart(:,2)
!     r03 = prim_to_cart(:,3)
!
!     call write_out(ulog,'reciprocal lattice vectors ',3,3,transpose(cart_to_prim))

! get atomic positions in cartesian coordinates and the neighbor list
      natoms=natom_prim_cell
      iatomcell(1:3,1:natoms)=0
      do i=1,natom_prim_cell
        write(*,*)'before xvmlt ',i
        call xvmlt(conv_to_cart,atomposconv(1,i),atompos(1,i),3,3,3)
!        atompos(:,i)=matmul(conv_to_cart,atomposconv(:,i))
        write(*,*)'before unitcell ',i
        call unitcell(cart_to_prim,prim_to_cart,atompos(1,i), atompos(1,i))
        iatomcell0(i)=i
        iatomneighbor(i,i)=0
      enddo
      write(*,*)'before dlatmat2'
!-------------------------------------------------------------------------------
! find symmetry of crystal.
! find point group of lattice
!k1   call dlatmat2(prim_to_cart,1d-6,lattpgcount,iop_matrix)
      call dlatmat2(prim_to_cart,1d-5,lattpgcount,iop_matrix)
      do i=1,lattpgcount ! number of symmetries of the bravais lattice
        do j=1,3
        do k=1,3
          temp(k,j)=iop_matrix(k,j,i)
        enddo
        enddo
        call xmatmlt(temp,cart_to_prim,temp,3,3,3,3,3,3)
        call xmatmlt(prim_to_cart,temp,op_matrix(1,1,i),3,3,3,3,3,3)
        write(ulog,*) 'Symmetry operation (op_matrix) #',i
        write(ulog,'(3(1x,f9.5))') op_matrix(1,:,i)
        write(ulog,'(3(1x,f9.5))') op_matrix(2,:,i)
        write(ulog,'(3(1x,f9.5))') op_matrix(3,:,i)
      enddo
      write(*,*)'after lattpgcount loop'

! find transformation matrices for k vectors
      do i=1,lattpgcount
         do j=1,3
         do k=1,3
           temp(k,j)=iop_matrix(k,j,i)
         enddo
         enddo
         call xmatmlt(temp,conv_to_prim,temp,3,3,3,3,3,3)
         call xmatmlt(prim_to_conv,temp,temp,3,3,3,3,3,3)
         call xmatinv(temp,tempi,ier)
         temp = tempi
         if(ier.ne.0) then
           write(*,*)'After xmatinv in line 127: ier=',ier
           call bomb
         endif

         write(*,*)' NOW GENERATING THE STARS OF K; lattpgcount=',i
         do j=1,3
         do k=1,3

            op_kmatrix(k,j,i)=nint(temp(j,k))
            if(ncmp(temp(j,k)-op_kmatrix(k,j,i)).ne.0) then
               write(ulog,*)'Op_kmatrix',temp(j,k),op_kmatrix(k,j,i)
      write(ulog,*)'****************************************************'
      write(ulog,*)'WARNING!! SOMETHNG WRONG IN EVALUATION OF Op_kmatrix'
      write(ulog,*)'****************************************************'
               call bomb
            endif

         enddo
         enddo
         write(ulog,*) 'Symmetry operation (op_kmatrix) #',i
         write(ulog,'(3(1x,f9.5))') op_kmatrix(1,:,i)
         write(ulog,'(3(1x,f9.5))') op_kmatrix(2,:,i)
         write(ulog,'(3(1x,f9.5))') op_kmatrix(3,:,i)
      enddo

      call write_out(ulog,' total # of symmetry ops ',lattpgcount)

! find elements of space group
! count them
      isgopcount=0
! try each element of point group of lattice
      ipntoploop: do ipntop=1,lattpgcount
! operate on atom 1.  v contains the coordinates of the atom after
! the operation
        call xvmlt(op_matrix(1,1,ipntop),atompos(1,1),v,3,3,3)
! try to map the rotated atom 1 onto every other atom of the same type
        iatomloop2: do iatom=1,natom_prim_cell
          if(iatomtype(iatom).eq.iatomtype(1))then
! find the fractional translation required in the space group element
            call xvsub(atompos(1,iatom),v,fract,3)
            call unitcell(cart_to_prim,prim_to_cart, fract,fract)
! now try this space group element on every other atom
            iatom2loop2: do iatom2=1,natom_prim_cell
              itype2=iatomtype(iatom2)
! operate on the atom: rotation followed by a fractional translation
! v contains the coordinates of the atom after the operation by the space
! group element.
              call xvmlt(op_matrix(1,1,ipntop),atompos(1,iatom2),v2, 3,3,3)
              call xvadd(fract,v2,v2,3)
              call unitcell(cart_to_prim,prim_to_cart,v2,v2)
! try to find another atom of the same type with these coordinates
              iatom3loop2: do iatom3=1,natom_prim_cell
                if(iatomtype(iatom3).eq.itype2)then
                  do i=1,3
! this one isn't it, try another one
                    if(ncmp(v2(i)-atompos(i,iatom3)).ne.0)  cycle iatom3loop2
                  enddo
! we found it, try the space group element on another atom
                  cycle iatom2loop2
                endif
! try to find another atom (iatom3)
              enddo iatom3loop2
! did not find it, try mapping atom 1 onto another atom
              cycle iatomloop2
! try the space group element on the next atom (iatom2)
            enddo iatom2loop2
! space group element was successful for every atom, save it
! count them
            isgopcount=isgopcount+1
! save rotation part of element
            isgop(isgopcount)=ipntop
! save translational part of element
            do i=1,3
              sgfract(i,isgopcount)=fract(i)
            enddo
! do next element of point group
            cycle ipntoploop
! try mapping atom 1 onto another atom
          endif
        enddo iatomloop2
! this element is not in space group
! do next element of point group (ipg)
      enddo ipntoploop

      write(ulog,*)' # of space group translation vectors=',isgopcount
      do isg=1,isgopcount
         write(ulog,4)' isg, isgop, sgfract=',isg,isgop(isg),sgfract(:,isg)
      enddo
 4    format(a,2i4,9(1x,g11.4))

!----------------------------------------------------------------------------
! find operators that take atoms into atoms in the primitive unit cell
! do each atom
      do iatom=1,natom_prim_cell
        itype=iatomtype(iatom)
! do each element of space group
        do isg=1,isgopcount
          ipntop=isgop(isg)
! operate on position of atom iatom
          call xvmlt(op_matrix(1,1,ipntop),atompos(1,iatom),v,3,3,3)
          call xvadd(sgfract(1,isg),v,v2,3)
          call unitcell(cart_to_prim,prim_to_cart,v2,v2)
! look for atom
          iatom2loop3: do iatom2=1,natom_prim_cell
            if(iatomtype(iatom2).eq.itype)then
              do i=1,3
                if(ncmp(v2(i)-atompos(i,iatom2)).ne.0)  cycle iatom2loop3
              enddo
! found it.  save it if not already previously found
              if(iatomop(iatom,iatom2).eq.0)then
                iatomop(iatom,iatom2)=ipntop
                call xvsub(atompos(1,iatom),v, atomopfract(1,iatom,iatom2),3)
                exit iatom2loop3
              endif
! try to match another atom (iatom2)
            endif
          enddo iatom2loop3
! next space group element
        enddo
! next atom (iatom)
      enddo

      write(ulog,*)' Atoms that are transformed to each other under space group ops'
      do iatom=1,natom_prim_cell
      do iatom2=1,natom_prim_cell
         i=iatomop(iatom,iatom2)
         if (i.ne.0) then
            write(ulog,5)'atoms # i,j, symmetry op#, sgfract=',iatom,iatom2,i,sgfract(:,i)
         endif
      enddo
      enddo
5 format(a,3i4,9(1x,g11.4))
      write(ulog,*)'--------------------------------------------------------- '
      write(ulog,*)'nsmax=',nsmax

!-------------------------------------------------------------------------------
! collect nearest neighbor atoms (up to a cutoff distance rcut(2))
! do each atom in primitive unit cell
      primloop: do iatom0=1,natom_prim_cell
        nd2save=0
        m=0
! collect distances to nearest neighbors
! do one shell of unit cells at a time
        shelloop: do n=0,nsmax
          foundone=.false.
          do i1=-n,n
          do i2=-n,n
          do i3=-n,n
! walk only on the facets of the cube of length 2n+1
            if(iabs(i1).ne.n.and.iabs(i2).ne.n.and.iabs(i3).ne.n)cycle
! do each atom in unit cell
            iloop: do i=1,natom_prim_cell
! skip the atom at the origin
              if(n.eq.0.and.i.eq.iatom0)cycle
! position of atom and distance squared
              r(1:3)=atompos(1:3,i)
              d2=0
              do j=1,3
                r(j)=atompos(j,i)+i1*prim_to_cart(j,1)  &
     &               +i2*prim_to_cart(j,2)+i3*prim_to_cart(j,3) -atompos(j,iatom0)
                d2=d2+r(j)**2
              enddo

! K1 reject if larger than rcut(2) .or. n > nshell
     !        if (d2 .gt. rcut(2)*rcut(2)) cycle

! did we find any new ones in this shell?
              if(m.eq.0)then
                foundone=.true.
              else if(ncmp(d2-d2save(nd2save)).eq.0 .or. d2.lt.d2save(nd2save))then
                foundone=.true.
              endif
! compare distance to atom with previously found distances
              do j=1,nd2save
! same distance: do another atom
                if(ncmp(d2-d2save(j)).eq.0)cycle iloop
! new distance: insert into list
                if(d2.lt.d2save(j))then
                  do k=maxneighbors,j+1,-1
                    d2save(k)=d2save(k-1)
                  enddo
                  if(nd2save.lt.maxneighbors) then
                     nd2save=nd2save+1
                  else
                   !  write(ulog,*)'nd2save >= maxneighbors; array d2save is going to overflow!!',nd2save
                   !  write(*,*)'nd2save >= maxneighbors; array d2save is going to overflow!!',nd2save
                     !exit !stop
                  endif

                  d2save(j)=d2
                  cycle iloop
                endif
              enddo
! new distance: insert at end of list if list is not already full
              if(nd2save.lt.maxneighbors)then
       !   the above is the cutoff, equivalent to  if(d2.lt.rcut(2)*rcut(2)) then
                nd2save=nd2save+1
                d2save(nd2save)=d2
       write(*,*)'nd2save,d2save=',nd2save,d2save
              else
                nshell=n-1
                exit shelloop  ! to avoid exceeding the size of d2save array
       !        write(ulog,*)'nd2save >= maxneighbors; or d2 too large!!',nd2save,d2
       !        write(*,*)'nd2save >= maxneighbors; or d2 too large!!',nd2save,d2
          !     stop
              endif
            enddo iloop
          enddo
          enddo
          enddo
! no new atom found in shell:  save previous shell and exit loop
          if(.not.foundone.and.(nd2save.eq.maxneighbors.or. d2.ge.rcut(2)*rcut(2))) then
            nshell=n-1
            exit shelloop
          endif
! list is filled
          if(nd2save.eq.maxneighbors.and.m.eq.0)m=1  
! we reached the last shell before we finished the list
! K1 below commented as maxshell is not used but rcut
       !  if(n.eq.maxneighbors)then
          if(n.eq.nsmax)then
            write(*,*)'Error:  maxneighbors or rcut(2) not large enough ',maxneighbors,rcut(2)
            write(*,*)'cutoff is imposed then by nsmax= ',nsmax
            nshell=nsmax
            write(*,*)'FORCE_CONSTANTS_INIT: nshell=',nshell
            exit shelloop
          endif
        enddo shelloop
        write(ulog,*)'FORCE_CONSTANTS_INIT: nshell=',nshell
!-----------------------------------------------------------------------------
! generate positions of neighbors
! do each shell of unit cells
        do n=0,nshell
          do i1=-n,n
          icell(1)=i1
          do i2=-n,n
          icell(2)=i2
          do i3=-n,n
          icell(3)=i3
            if(iabs(i1).ne.n.and.iabs(i2).ne.n.and.iabs(i3).ne.n)cycle
! do each atom in unit cell
            do i=1,natom_prim_cell
! skip atom at origin
              if(n.eq.0.and.i.eq.iatom0)cycle
! position of atom
              r(1:3)=atompos(1:3,i)
              d2=0
              do j=1,3
                r(j)=atompos(j,i)+i1*prim_to_cart(j,1)  &
     &               +i2*prim_to_cart(j,2)+i3*prim_to_cart(j,3)
                d2=d2+(r(j)-atompos(j,iatom0))**2
              enddo
! find in list
              do m=1,nd2save
                if(ncmp(d2-d2save(m)).eq.0)then
                  call findatom(icell,i,iatom)
                  if(iatom.eq.0)then
                    natoms=natoms+1
                    if(natoms.gt.maxatoms)then
                      write(6,*)'Warning: in force_constants_init:  '   &
     &                     //'the value of maxatoms needs to be larger.'
                      imaxat=1
                      return
                    endif
                    atompos(1:3,natoms)=r(1:3)
                    iatomcell(1:3,natoms)=icell(1:3)
                    iatomcell0(natoms)=i
                    iatomneighbor(iatom0,natoms)=m
                  else
                    iatomneighbor(iatom0,iatom)=m
                  endif
                endif
              enddo
            enddo
          enddo
          enddo
          enddo
        enddo
! next inequivalent atom
      enddo primloop
      maxneighbors=nd2save
      write(ulog,*) 'FC_INIT: maxneighbors now changed to ',maxneighbors
      imaxat=0
      end subroutine force_constants_init
      !******************************************************************************
!! get force constants
      subroutine force_constants(nrank,amat,iatomd,ixyzd,   &
     &     ntermsindep,iatomtermindep,ixyztermindep,   &
     &     ntermall,iatomtermall,ixyztermall,   &
     &     ntermszero,iatomtermzero,ixyztermzero,   &
     &     maxrank,maxterms,maxtermsindep,maxtermszero, &
     &     ierz,iert,ieri) !,ierg)
! find relationships between terms of the form d^nU/dxidxjdxk...
! where U is the total energy of the crystal and xi is a coordinate (x=x,y,z)
! of the ith atom.

!! arguments:
!!     nrank (input), order of derivative
!!     iatomd(i) (input), atom at ith location in denominator
!!     ixyzd(i) (input), coordinate (1,2,3=x,y,z) at ith location in denominator
!!     ntermsindep (output), number of independent terms
!!     iatomtermindep(i,j) (output), atom at ith location in denominator of
!!          jth independent term
!!     ixyztermindep(i,j) (output), coordinate at ith location in denomonitor
!!          of jth independent term
!!     ntermall (output), total number of nonzero terms
!!     iatomtermall(i,j) (output), atom at ith location in denominator of
!!          jth nonzero term
!!     ixyztermall(i,j) (output), coordinate at ith location in denomonitor of
!!          jth nonzero term
!!     amat(i,j) (output), coefficient of jth independent term in equation for
!!          kth nonzero term
!!     ntermszero (output), number of zero terms generated
!!     iatomtermzero(i,j) (output), atom at ith location in denominator of
!!          jth zero term
!!     ixyztermzero(i,j) (output), coordinate at ith location in denomonitor
!!          of jth zero term
!!     maxrank (input), number of rows in arrays iatomtermindep, ixyztermindep,
!!          iatomtermall, ixyztermall, iatomtermzero and ixyztermzero
!!     maxterms (input), number of columns in arrays iatomtermall and
!!          ixyztermall
!!     maxtermsindep (input), number of columns in arrays iatomtermindep and
!!          ixyztermindep
!!     maxtermszero (input), number columns in arrays iatomtermzero and
!!          ixyztermzero

!      use force_constants_module
      use atoms_force_constants
      implicit none
      integer, intent(in):: nrank,maxrank,maxterms,maxtermsindep,maxtermszero
      integer, intent(out):: ierz,iert,ieri,ntermsindep,ntermall,ntermszero
      integer, intent(in):: iatomd(nrank),ixyzd(nrank)
     integer, intent(out)::  &
     &   iatomtermindep(maxrank,maxterms)    ,ixyztermindep(maxrank,maxterms), &
     &   iatomtermall  (maxrank,maxterms)    ,ixyztermall  (maxrank,maxterms), &
     &   iatomtermzero (maxrank,maxtermszero),ixyztermzero (maxrank,maxtermszero)
      integer i,j,k,m,n,iatom(nrank),ncmp,ixyz(nrank),   &
     &   iatomterm(maxrank,maxterms),ixyzterm(maxrank,maxterms),jterm,   &
     &     icell(3),iv(3),isg,msave,nn,   &
     &     irank,k2,neqs,ifactorial,msave2,j2,   &
     &     npermute,ixyzfirst(nrank),ixyz4(nrank),   &
     &     mapdep(maxterms),mapterms(maxterms),   &
     &     maptermsindep(maxterms),mapindepterms(maxterms)
      logical foundit(maxterms),firstone,zero(maxterms)
      double precision v(3),amp
      real(8), intent(out):: amat(maxterms,maxtermsindep)
      integer, allocatable :: ipermute(:,:)
      real(8), allocatable :: eqs(:,:),eqs2(:,:)

      allocate(eqs(maxterms,maxterms),eqs2(maxterms,maxterms))
      iert=0 ; ieri=0 ; ierz=0

! check if input values are valid
      if(nrank.gt.maxrank)then
        write(6,*)'Error in force_constants: nrank > maxrank'
        stop
      endif
      do i=1,nrank
        if(ixyzd(i).lt.1.or.ixyzd(i).gt.3)then
          write(6,*)'Error in force_constants:  invalid value in ixyzd'
          stop
        endif
        if(iatomd(i).lt.1.or.iatomd(i).gt.natoms)then
          write(6,*)'Error in force_constants:  invalid value in iatomd'
          stop
        endif
      enddo
! get permutations of nrank items
      npermute=ifactorial(nrank)
      n=1
      do i=2,nrank
        n=n*i
      enddo
      allocate(ipermute(nrank,n))
      call permutations(nrank,ipermute,nrank,n)

      eqs=0
      ntermall=1
      call unique_force_constant(nrank,iatomd,ixyzd,iatomterm,ixyzterm)
      jterm=1
      neqs=0
      ntermszero=0
      jtermloop: do while(jterm.le.ntermall)
!        write(6,*)jterm,ntermall
! do each space-group operator
      sgloop: do isg=1,isgopcount
        neqs=neqs+1
        if(neqs.gt.maxterms)then
          write(6,*)'Warning: in force_constants: maxterms too small'
          iert=iert+1
          deallocate(ipermute,eqs,eqs2)
          return
        endif
        eqs(neqs:maxterms,1:ntermall)=0
        eqs(neqs,jterm)=1
! operate on each item in denominator
        do irank=1,nrank
          call xvmlt(op_matrix(1,1,isgop(isg)),  &
     &         atompos(1,iatomterm(irank,jterm)),v,3,3,3)
          call xvadd(sgfract(1,isg),v,v,3)   ! K1: could be:  ! v=v+sgfract(:,isg)
! find atom
          call findatom2(v,iatom(irank))
          if(iatom(irank).eq.0)then
            write(6,*)'Error0 in force_constants: atom not found'
            write(6,*)'Increase the # of NN shells or check POSCAR-OUTCAR'
            write(6,*)'noncorresponding vector=',v
            stop
          endif
        enddo
! rotate coordinates
        do k=1,nrank
          do i=1,3
            if(ncmp(op_matrix(i,ixyzterm(k,jterm),isgop(isg))).ne.0)then
              ixyz4(k)=i
              ixyzfirst(k)=i
              exit
            endif
            if(i.eq.3)stop 'This cannot happen!'
          enddo
        enddo
        ixyz4(nrank)=ixyz4(nrank)-1
        whileloop: do while(.true.)
          kloop: do k=nrank,1,-1
            do while(.true.)
              ixyz4(k)=ixyz4(k)+1
              if(ixyz4(k).gt.3)then
                if(k.eq.1)then
                  exit whileloop
                else
                  cycle kloop
                endif
              endif
              if(ncmp(op_matrix(ixyz4(k),ixyzterm(k,jterm),isgop(isg))) .ne.0)exit
            enddo
            do k2=k+1,nrank
              ixyz4(k2)=ixyzfirst(k2)
            enddo
            exit
          enddo kloop
! amplitude
          amp=1
          do k=1,nrank
            amp=amp*op_matrix(ixyz4(k),ixyzterm(k,jterm),isgop(isg))
          enddo
          if(ncmp(amp).eq.0)stop 'This cannot happen!'
! get unique denominator
          call unique_force_constant(nrank,iatom,ixyz4,iatom,ixyz)
! look for variable
          iloop: do i=1,ntermall
            do k=1,nrank
              if(iatomterm(k,i).ne.iatom(k).or.ixyzterm(k,i).ne.ixyz(k))  exit
! foundit: add in amp
              if(k.eq.nrank)then
                eqs(neqs,i)=eqs(neqs,i)-amp
                exit iloop
              endif
            enddo
! did not find it: add it and add in amp: this calculates ntermall
            if(i.eq.ntermall)then
              ntermall=ntermall+1
              if(ntermall.gt.maxterms)then
           write(6,*) 'Warning2: in force_constants: maxterms too small',maxterms
                iert=iert+1
                deallocate(ipermute,eqs,eqs2)
                return
              endif
              iatomterm(1:nrank,ntermall)=iatom(1:nrank)
              ixyzterm(1:nrank,ntermall)=ixyz(1:nrank)
              eqs(neqs,ntermall)=eqs(neqs,ntermall)-amp
            endif
          enddo iloop
        enddo whileloop
! next space-group operator
      enddo sgloop
! solve simultaneous equations
      call xrowop2(eqs,neqs,ntermall,maxterms,maxterms)
! find number of equations
      iloop2: do i=1,neqs
        do j=1,ntermall
          if(ncmp(eqs(i,j)).ne.0)exit
          if(j.eq.ntermall)then
            neqs=i-1
            exit iloop2
          endif
        enddo
      enddo iloop2
! next seed
      jterm=jterm+1
      enddo jtermloop
! put terms in order
      foundit=.false.
      eqs2=0
      do i=1,ntermall
        firstone=.true.
        do j=1,ntermall
          if(foundit(j))cycle
          if(firstone)then
            n=j
            firstone=.false.
          else
            do k=1,nrank
              if(iatomterm(k,j).gt.iatomterm(k,n))then
                n=j
                exit
              else if(iatomterm(k,j).lt.iatomterm(k,n))then
                exit
              endif
              if(k.eq.nrank)then
                do m=1,nrank
                  if(ixyzterm(m,j).gt.ixyzterm(m,n))then
                    n=j
                    exit
                  else if(ixyzterm(m,j).lt.ixyzterm(m,n))then
                    exit
                  endif
                enddo
              endif
            enddo
          endif
        enddo
        foundit(n)=.true.
        mapterms(i)=n
        eqs2(1:neqs,i)=eqs(1:neqs,n)
      enddo
      call xrowop2(eqs2,neqs,ntermall,maxterms,maxterms)
! find independent and dependent terms
      mapdep=0
      do i=1,neqs
        do j=1,ntermall
          if(ncmp(eqs2(i,j)).ne.0)then
            mapdep(j)=i
            exit
          endif
          if(j.eq.ntermall)then
            write(6,*)'Error in force_constants: zero equation'
            stop
          endif
        enddo
      enddo
! save independent terms: this loop calculates ntermsindep
      ntermsindep=0
      ntiloop: do i=ntermall,1,-1
        if(mapdep(i).eq.0)then
          ntermsindep=ntermsindep+1
          if(ntermsindep.gt.maxtermsindep)then
            write(6,*) 'Warning: in force_constants: maxtermsindep too small'
            ieri=ieri+1
            deallocate(ipermute,eqs,eqs2)
            return  !stop
          endif
          mapindepterms(ntermsindep)=i
          maptermsindep(i)=ntermsindep
          iatomtermindep(1:nrank,ntermsindep)= iatomterm(1:nrank,mapterms(i))
          ixyztermindep(1:nrank,ntermsindep)=  ixyzterm(1:nrank,mapterms(i))
        endif
      enddo ntiloop
! no independent terms: all terms are zero
      if(ntermsindep.eq.0)then
        if(ntermall.gt.maxtermszero)then
         write(6,*)'Warning: in force_constants: maxtermszero too small'
          ierz=ierz+1
          deallocate(ipermute,eqs,eqs2)
          return !stop
        endif
        ntermszero=ntermall
        iatomtermzero(1:nrank,1:ntermall)=iatomterm(1:nrank,1:ntermall)
        ixyztermzero(1:nrank,1:ntermall)=ixyzterm(1:nrank,1:ntermall)
        ntermall=0
        deallocate(eqs,eqs2)
        return
      endif
! look for zero terms
      zero=.false.
      do i=1,ntermall
        if(mapdep(i).ne.0)then
          do j=1,ntermsindep
            if(ncmp(eqs2(mapdep(i),mapindepterms(j))).ne.0)exit
            if(j.eq.ntermsindep)then
              zero(i)=.true.
            endif
          enddo
        endif
      enddo
! do each term
      m=0
      do n=ntermall,1,-1
        nn=mapterms(n)
! find all permutations
        msave=m
        permuteloop: do i=1,npermute
          do j=1,nrank
            iatomtermall(j,m+1)=iatomterm(ipermute(j,i),nn)
            ixyztermall(j,m+1)=ixyzterm(ipermute(j,i),nn)
          enddo
! look for it
          do k=msave+1,m
            do j=1,nrank
              if(iatomtermall(j,k).ne.iatomtermall(j,m+1))exit
              if(ixyztermall(j,k).ne.ixyztermall(j,m+1))exit
              if(j.eq.nrank)cycle permuteloop
            enddo
          enddo
! new term: add it
          m=m+1
          if(m.eq.maxterms)then
            write(6,*)'Warning: in fcs: maxterms too small',maxterms
            iert=iert+1
            deallocate(ipermute,eqs,eqs2)
            return
!            call bomb
          endif
        enddo permuteloop
! find all translations
        msave2=m
        translateloop: do i=msave+1,msave2
          do j2=1,nrank
            if(iatomcell0(iatomtermall(j2,i)).eq.iatomtermall(j2,i))  cycle
            iv(1:3)=iatomcell(1:3,iatomtermall(j2,i))
! move atoms and identify them
            do k=1,nrank
              do k2=1,3
                icell(k2)=iatomcell(k2,iatomtermall(k,i))-iv(k2)
              enddo
              call findatom(icell,iatomcell0(iatomtermall(k,i)), iatomtermall(k,m+1))
              if(iatomtermall(k,m+1).eq.0)then
                write(6,*)'Error in force_constants: atom not found'
                stop
              endif
              ixyztermall(k,m+1)=ixyztermall(k,i)
            enddo
! look for it
            do k=msave+1,m
              do j=1,nrank
                if(iatomtermall(j,k).ne.iatomtermall(j,m+1))exit
                if( ixyztermall(j,k).ne. ixyztermall(j,m+1))exit
                if(j.eq.nrank)cycle translateloop
              enddo
            enddo
! new term: add it
            m=m+1
            if(m.eq.maxterms)then
              write(6,*)'Warning3: in fcs: maxterms too small',maxterms
              iert=iert+1
              deallocate(ipermute,eqs,eqs2)
              return
!              call bomb
            endif
          enddo
        enddo translateloop
! zero terms
        if(zero(n))then
          if(ntermszero+m-msave.gt.maxtermszero)then
         write(6,*)'Warning: in force_constants: maxtermszero too small'
            ierz=ierz+1
            deallocate(ipermute,eqs,eqs2)
            return     !       stop
          endif
          iatomtermzero(1:nrank,ntermszero+1:ntermszero+m-msave)=   &
     &         iatomtermall(1:nrank,msave+1:m)
          ixyztermzero(1:nrank,ntermszero+1:ntermszero+m-msave)=   &
     &         ixyztermall(1:nrank,msave+1:m)
          ntermszero=ntermszero+m-msave
          m=msave
! independent terms
        else if(mapdep(n).eq.0)then
          amat(msave+1:m,maptermsindep(n))=1
! dependent terms
        else
          do i=1,ntermsindep
            amat(msave+1:m,i)=-eqs2(mapdep(n),mapindepterms(i))
          enddo
        endif
      enddo
      ntermall=m
      deallocate(ipermute,eqs,eqs2)
      ierz=0
      iert=0
      ieri=0

      end subroutine force_constants
!****************************************************************************
!! bring a force constant to a unique form: atoms
      subroutine unique_force_constant(nrank,iatomin,ixyzin,iatomout,ixyzout)

!! arguments:
!!     nrank (input), order of derivative
!!     iatomin(i) (input), atom at ith location in denominator
!!     ixyzin(i) (input), x,y,z coordinate of atom
!!     iatomout(i) (output), atom at ith location in denominator
!!     ixyzout(i) (output), x,y,z coordinate of atom

      use atoms_force_constants
!      use force_constants_module
      implicit none
      integer, intent(in):: nrank,iatomin(nrank),ixyzin(nrank)
      integer, intent(out):: iatomout(nrank),ixyzout(nrank)
      integer k,m,iatomtemp(nrank),irank,iv(3),icell(3),iatom,k2,k3,  &
     &     iatomtemp2(nrank),ixyztemp(nrank),ixyztemp2(nrank),ixyz
      logical firsttime

      firsttime=.true.
! bring atom in each item to unit cell at origin
      irankloop: do irank=1,nrank
        iv(1:3)=iatomcell(1:3,iatomin(irank))
! move atoms and identify them
        do k=1,nrank
          do m=1,3
            icell(m)=iatomcell(m,iatomin(k))-iv(m)
          enddo
          call findatom(icell,iatomcell0(iatomin(k)),iatomtemp(k))
          if(iatomtemp(k).eq.0)then
            write(6,*)'Error in force_constants: atom not found'
            stop
          endif
          ixyztemp(k)=ixyzin(k)
! put atoms in order
          do k2=1,k-1
            if(iatomtemp(k).lt.iatomtemp(k2).or.  &
     &           (iatomtemp(k).eq.iatomtemp(k2).and.  &
     &           ixyztemp(k).lt.ixyztemp(k2))) then
              iatom=iatomtemp(k)
              ixyz=ixyztemp(k)
              do k3=k,k2+1,-1
                iatomtemp(k3)=iatomtemp(k3-1)
                ixyztemp(k3)=ixyztemp(k3-1)
              enddo
              iatomtemp(k2)=iatom
              ixyztemp(k2)=ixyz
              exit
            endif
          enddo
        enddo
!compared with saved term
        if(firsttime)then
          firsttime=.false.
          iatomtemp2(1:nrank)=iatomtemp(1:nrank)
          ixyztemp2(1:nrank)=ixyztemp(1:nrank)
        else
          do k=1,nrank
            if(iatomtemp(k).lt.iatomtemp2(k))then
              iatomtemp2(1:nrank)=iatomtemp(1:nrank)
              ixyztemp2(1:nrank)=ixyztemp(1:nrank)
              cycle irankloop
            else if(iatomtemp(k).gt.iatomtemp2(k))then
              cycle irankloop
            endif
            if(k.eq.nrank)then
              do m=1,nrank
                if(ixyztemp(m).lt.ixyztemp2(m))then
                  iatomtemp2(1:nrank)=iatomtemp(1:nrank)
                  ixyztemp2(1:nrank)=ixyztemp(1:nrank)
                  cycle irankloop
                else if(ixyztemp(m).gt.ixyztemp2(m))then
                  cycle irankloop
                endif
              enddo
            endif
          enddo
        endif
      enddo irankloop
      iatomout(1:nrank)=iatomtemp2(1:nrank)
      ixyzout(1:nrank)=ixyztemp2(1:nrank)

      end subroutine unique_force_constant

!****************************************************************************
! routines
      subroutine findatom(icell,icell0,iatom)
!! find atom in data base
!! arguments:
!!     icell(i) (input), linear combination of basis vectors of the primitive
!!          lattice that takes us to the unit cell containing the ith atom
!!     icell0 (input), identity of equivalent atom in unit cell at origin
!!     iatom (output), location of atom in data base.  Returns zero if not found
!!          in data base
      use atoms_force_constants
!      use force_constants_module
      implicit none
      integer, intent(in):: icell(3),icell0
      integer, intent(out):: iatom
      integer i,j

      do i=1,natoms
        if(icell0.ne.iatomcell0(i))cycle
        do j=1,3
          if(icell(j).ne.iatomcell(j,i))exit
          if(j.eq.3)then
            iatom=i
            return
          endif
        enddo
      enddo
      iatom=0
      end subroutine findatom
!--------------------------------------------------------------------------------
      subroutine findatom2(pos,iatom)
!! find atom in data base
!! arguments:
!!     pos(i) (input), ith cartesian coordinate of atomic position
!!     iatom (output), location of atom in data base.  Returns zero if not found
!!          in data base
      use atoms_force_constants
!      use force_constants_module
      implicit none
      integer iatom,i,j,ncmp
      double precision pos(3)
      do i=1,natoms
        do j=1,3
          if(ncmp(atompos(j,i)-pos(j)).ne.0)exit
          if(j.eq.3)then
            iatom=i
            return
          endif
        enddo
      enddo
      iatom=0
      end subroutine findatom2
!--------------------------------------------------------------------------------
      function ncmp(x)
      implicit none
!
!!	COMPARE X WITH ZERO
!!	NCMP=0 IF X IS CLOSE ENOUGH TO ZERO
!!	NCMP=1 OTHERWISE
!!	X IS REAL
!
      integer ncmp
      real(8) x !,delta
!     data delta/1.e-3/
      ncmp=0
!     if(abs(x).gt.delta)ncmp=1
      if(abs(x).gt.1d-5)ncmp=1  ! was originally 1d-6 ! K1
      return
      end function ncmp
!--------------------------------------------------------------------------------
      subroutine xmatinv(xmatin,xmatout,ier)
      implicit none

!! invert a 3 by 3 matrix

      double precision xmatin(3,3),xmatout(3,3),buffer(3,3),x
      integer indx(3),ier,n,i,j

! dimension of matrix
      n=3
! clear error flag
      ier=0
      do i=1,n
        do j=1,n
          xmatout(i,j)=0
          buffer(i,j)=xmatin(i,j)
        enddo
        xmatout(i,i)=1
      enddo
! decomposition
      call ludcmp(buffer,n,n,indx,x)
! singular matrix
      if(x.eq.0.0d0)then
        ier=1
        return
      endif
! inverse matrix
      do j=1,n
        call lubksb(buffer,n,n,indx,xmatout(1,j))
      enddo

      end subroutine xmatinv
!------------------------------------------------------------------
! The following routines are from Numerical Recipes
      subroutine ludcmp(a,n,np,indx,d)
      implicit none
      integer nmax,np,n
      double precision tiny
      parameter (nmax=3,tiny=1.0d-20)
      double precision a(np,np),vv(nmax),d,aamax,dum,sum
      integer indx(n),i,j,k,imax,ncmp
      d=1
      do i=1,n
        aamax=0
        do j=1,n
          if(dabs(a(i,j)).gt.aamax)aamax=dabs(a(i,j))
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
      implicit none
      integer n,np
      double precision a(np,np),b(n),sum
      integer indx(n),ii,i,j,ll
      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if(ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if(sum.ne.0.0d0)then
          ii=i
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do j=i+1,n
            sum=sum-a(i,j)*b(j)
          enddo
        endif
        b(i)=sum/a(i,i)
      enddo

      end subroutine lubksb
!------------------------------------------------------------------------------
      subroutine xmatmlt(x1,x2,x3,nrow1,ncol1,ncol2,nr1,nr2,nr3)

!! multiply two real matrices, x3=x1*x2
! double precision version
! arguments:
!     x1,x2 (input), first and second matrix
!     x3 (output), product x1*x2
!     nrow1 (input), number of rows in x1, also the number of rows in x3
!     ncol1 (input), number of columns in x1, also the number of
!          rows in x2
!     ncol2 (input), number of columns in x2, also the number of
!          columns in x3
!     nr1 (input), number of rows in the physical array x1
!     nr2 (input), number of rows in the physical array x2
!     nr3 (input), number of rows in the physical array x3

      implicit none
      integer i,j,k,nrow1,ncol1,ncol2,nr1,nr2,nr3
      double precision x1(nr1,ncol1),x2(nr2,ncol2),x3(nr3,ncol2)
      real(8), allocatable :: x(:,:)

      allocate(x(nrow1,ncol2))
      do i=1,ncol2
      do j=1,nrow1
        x(j,i)=0
        do k=1,ncol1
          x(j,i)=x(j,i)+x1(j,k)*x2(k,i)
        enddo
      enddo
      enddo
      do i=1,ncol2
      do j=1,nrow1
        x3(j,i)=x(j,i)
      enddo
      enddo
      deallocate(x)

      end subroutine xmatmlt
!-------------------------------------------------------------------------------
      subroutine xvmlt(x,v1,v2,nrow,ncol,nr)

!! multiply a double precision vector by a double precision matrix, v2=x*v1
! arguments:
!     x (input), matrix
!     v1 (input), vector
!     v2 (output), product x*v1
!     nrow (input), number of rows in x, also the number of rows in v2
!     ncol (input), number of columns in x, also the number of rows in v1
!     nr (input), number of rows in the physical array x

      implicit none
      integer nrow,ncol,nr,i,j
      double precision x(nr,ncol),v1(ncol),v2(nrow)
      real(8), allocatable:: v(:)

      allocate(v(nrow))
      do i=1,nrow
        v(i)=0
        do j=1,ncol
          v(i)=v(i)+x(i,j)*v1(j)
        enddo
      enddo
      v2(1:nrow)=v(1:nrow)
      deallocate(v)

      end subroutine xvmlt
!------------------------------------------------------------------------------
      subroutine unitcell(cart_to_prim,prim_to_cart,v1,v2)

!! bring a point v1 into v2 in the unit cell at the origin, between 0 and R
!! cart_to_prim has the reciprocal vectors in lines 1,2,3
!! prim_to_cart has the translation vectors in real space in columns 1,2,3

      implicit none
      double precision cart_to_prim(3,3),prim_to_cart(3,3),v1(3),  &
     &     v2(3),buff(3)
      integer i,ncmp
! change coordinates of point to linear combination of basis vectors of the
! primitive lattice
      call xmatmlt(cart_to_prim,v1,buff,3,3,1,3,3,3)
! in the unit cell at the origin, the coefficient must be greater than or
! equal to zero and less than one.
      do i=1,3
        do while(buff(i).gt.1.0d0.or.ncmp(buff(i)-1).eq.0)
          buff(i)=buff(i)-1
        enddo
        do while(buff(i).lt.0.0d0.and.ncmp(buff(i)).ne.0)
          buff(i)=buff(i)+1
        enddo
      enddo
! return to cartesian coordinates
      call xmatmlt(prim_to_cart,buff,v2,3,3,1,3,3,3)

      end subroutine unitcell
!----------------------------------------------------------------------------
      subroutine dlatmat2(cart,eps,nmatrices,matrices)

!! find symmetry matrices for a given BRAVAIS lattice

!! arguments:
!!     cart(i,j) (input), ith cartesian component of jth basis vector
!!     eps (input), tolerance for length
!!     nmatrices (output), number of matrices
!!     matrices(i,j,k) (output), kth matrix

      implicit none

      integer nmax
      parameter(nmax=400)

      real(8), intent(in)  :: cart(3,3),eps
      integer, intent(out) :: nmatrices,matrices(3,3,48)
      integer n,i,j,j1,j2,j3,k,m,i1,i2,i3,nshort(3),ndet,itrans(3,3),  &
     &      ichoose(3),ishort(3,nmax,3)
      real(8) vshort(3,nmax,3),v(3),xmax,vlength,x,d,abc(3,3),dshort(nmax,3)
      logical foundone,tried(48,48)

! some initialization
      dshort=0
      ishort=0
      vshort=0
      nshort=0
      do i=1,3
        abc(i,i)=vlength(3,cart(1,i))
      enddo
      do i=2,3
        do j=1,i-1
          call xvsub(cart(1,i),cart(1,j),v,3)
          abc(i,j)=vlength(3,v)
          abc(j,i)=abc(i,j)
        enddo
      enddo
      i=0
      foundone=.true.
! longest lattice parameter
      xmax=dmax1(abc(1,1),abc(2,2),abc(3,3))+eps
! try each shell until every vector in a shell is longer than the longest
! lattice parameter
      do while(foundone)
        i=i+1
        foundone=.false.
! find all lattice vectors in shell
        do j1=-i,i
        do j2=-i,i
        do j3=-i,i
        if(iabs(j1).eq.i.or.iabs(j2).eq.i.or.iabs(j3).eq.i)then
! length of lattice vector
          v=0
          do k=1,3
            v(k)=v(k)+j1*cart(k,1)
            v(k)=v(k)+j2*cart(k,2)
            v(k)=v(k)+j3*cart(k,3)
          enddo
          d=vlength(3,v)
! if shorter than longest lattice parameter, then do next shell too
          if(d.lt.xmax)foundone=.true.
! check each lattice parameter a,b,c
          do k=1,3
! equal to length of lattice parameter to within tolerance
            if(dabs(d-abc(k,k)).lt.eps)then
! count them
              nshort(k)=nshort(k)+1
              if(nshort(k).gt.nmax)then
                write(6,*)'Error in dlatmat:  nmax too small'
                call bomb
              endif
! length
              dshort(nshort(k),k)=d
! dimensionless coordinates
              ishort(1,nshort(k),k)=j1
              ishort(2,nshort(k),k)=j2
              ishort(3,nshort(k),k)=j3
! cartesian coordinates
              vshort(1,nshort(k),k)=v(1)
              vshort(2,nshort(k),k)=v(2)
              vshort(3,nshort(k),k)=v(3)
            endif
          enddo
! next vector in shell
        endif
        enddo
        enddo
        enddo
! next shell
      enddo

! try mappings of basis vectors onto vectors the "same" length
      nmatrices=1
      do i1=1,nshort(1)
      ichoose(1)=i1
      itrans(1:3,1)=ishort(1:3,i1,1)
      do i2=1,nshort(2)
      ichoose(2)=i2
      itrans(1:3,2)=ishort(1:3,i2,2)
      i3loop: do i3=1,nshort(3)
      ichoose(3)=i3
      itrans(1:3,3)=ishort(1:3,i3,3)
! determinant of the transformation matrix must be equal to 1
      if(iabs(ndet(itrans)).ne.1)cycle
! lengths of differences of lattice vectors must match to within tolerance
      do i=2,3
        do j=1,i-1
          call xvsub(vshort(1,ichoose(i),i),vshort(1,ichoose(j),j),v,3)
          x=vlength(3,v)
          if(dabs(x-abc(i,j)).gt.eps)cycle i3loop
        enddo
      enddo
! found a transformation:  count them and save it
! if this is the identity op, just put it into the first matrix where we
! have reserved a place for it
      if(itrans(1,1)+itrans(2,2)+itrans(3,3).eq.3)then
        matrices(1:3,1:3,1)=itrans(1:3,1:3)
      else
        nmatrices=nmatrices+1
        if(nmatrices.gt.48)then
          write(6,*)'Error in dlatmat2: more than 48 point operators'
          call bomb
        endif
        matrices(1:3,1:3,nmatrices)=itrans(1:3,1:3)
      endif
! next mapping
      enddo i3loop
      enddo
      enddo
! find any additional matrices by multiplication
      foundone=.true.
      tried=.false.
      do while(foundone)
        foundone=.false.
        do i=1,nmatrices
          do j=1,nmatrices
            if(.not.tried(i,j))then
              tried(i,j)=.true.
!             call matmlt(matrices(1,1,i),matrices(1,1,j),itrans)
              itrans=matmul(matrices(:,:,i),matrices(:,:,j))
              kloop: do k=1,nmatrices
                mloop: do m=1,3
                do n=1,3
                  if(matrices(n,m,k).ne.itrans(n,m))exit mloop
                  if(m.eq.3.and.n.eq.3)exit kloop
                enddo
                enddo mloop
                if(k.eq.nmatrices)then
                  foundone=.true.
                  nmatrices=nmatrices+1
                  if(nmatrices.gt.48)then
                    write(6,*)'Error in dlatmat2: '//  &
     &                   'more than 48 point operators'
                    call bomb
                  endif
                  matrices(1:3,1:3,nmatrices)=itrans(1:3,1:3)
                endif
              enddo kloop
            endif
          enddo
        enddo
      enddo

      end subroutine dlatmat2
!--------------------------------------------------------------------------------
      subroutine xvsub(v1,v2,v3,nrow)

!! subtract two real vectors: v3=v1-v2
! double precision version
! arguments:
!     v1,v2 (input), vectors
!     v3 (output), vector v1-v2
!     nrow (input), number of rows in each vector

      implicit none
      integer i,nrow
      double precision v1(nrow),v2(nrow),v3(nrow)
      do i=1,nrow
        v3(i)=v1(i)-v2(i)
      enddo

      end subroutine xvsub
!------------------------------------------------------------------------------
      subroutine xvadd(v1,v2,v3,nrow)

!! add two real vectors: v3=v1+v2
! double precision version
! arguments:
!     v1,v2 (input), vectors
!     v3 (output), vector v1+v2
!     nrow (input), number of rows in each vector

      implicit none
      integer i,nrow
      double precision v1(nrow),v2(nrow),v3(nrow)
      do i=1,nrow
        v3(i)=v1(i)+v2(i)
      enddo

      end subroutine xvadd
!-------------------------------------------------------------------------------
      function vlength(n,v)
      implicit none
      integer n
      double precision v(n),vlength,x
      integer i
      x=0
      do i=1,n
        x=x+v(i)**2
      enddo
      vlength=dsqrt(x)
      end function vlength
!-------------------------------------------------------------------------------
      subroutine bomb
      implicit none
      write(6,'(a)')'This program has bombed.'
      write(6,'(a)')'exit'
      stop
      end subroutine bomb
!------------------------------------------------------------------------------
      function ndet(mat)
      implicit none
!
!!	FIND THE DETERMINANT OF A 3X3 MATRIX MAT
!
      integer ndet,mat(3,3)
      ndet=mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))  &
     & -mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1))  &
     & +mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))

      end function ndet
!------------------------------------------------------------------------------
      subroutine permutations(n,ipermutations,nr,nc)
!! find all permutations of n objects
!! arguments:
!!     n (input), number of objects
!!     ipermutations(i,j) (output), object in ith location of jth permutation
!!     nr (input), number of rows in ipermutations
!!     nc (input), number of columns in ipermutations

      implicit none
      integer nr,nc
      integer n,ipermutations(nr,nc),ip(10),np,ifactorial, i,j,k,m !,mm,nexc !,ip2(10)
      logical used(10)

! check for valid input
      if(n.gt.10)then
        write(6,*)'error in permutations: n is too large'
        call bomb
      endif
      if(n.gt.nr)then
        write(6,*)'error in permutations: nr is too small'
        call bomb
      endif
      np=ifactorial(n)
      if(np.gt.nc)then
        write(6,*)'error in permutations: nc is too small'
        call bomb
      endif
! first permutation
      do i=1,n
        ip(i)=i
      enddo
      ipermutations(1:n,1)=ip(1:n)
! done if only one object
      if(n.eq.1)return
! mark which objects are used in the first n-1 locations
      used(1:n-1)=.true.
      used(n)=.false.
! count permutations
      np=1
! try changing the object in a location, starting at location n-1
    loop1: do while (np.lt.nc)
      do i=n-1,1,-1
! moving this object: not at this location anymore
        used(ip(i))=.false.
! look for a vacant location
        do j=ip(i)+1,n
! found one
          if(.not.used(j))then
! put object in that location
            ip(i)=j
            used(j)=.true.
! put the remaining objects in the remaining locations
            do k=i+1,n
              do m=1,n
                if(.not.used(m))then
                  ip(k)=m
                  if(k.ne.n)used(m)=.true.
                  exit
                endif
                if(m.eq.n)then
                  write(*,*)'ERROR in rearranging: m=n ',m
                  call bomb
                endif
              enddo
            enddo
! count permutations and save it
            np=np+1
            ipermutations(1:n,np)=ip(1:n)
            cycle loop1 !goto 1
          endif
        enddo
      enddo
    enddo loop1

      end subroutine permutations
!-------------------------------------------------------------------------------
!! factorial:  ifactorial = n!
      function ifactorial(n)
      implicit none
      integer ifactorial,n,i,m
      if(n.lt.0)then
        call bomb
      else if(n.eq.0)then
        ifactorial=1
      else
        m=1
        do i=2,n
          m=m*i
        enddo
        ifactorial=m
      endif

      end function ifactorial
!-------------------------------------------------------------------------------
      subroutine xrowop2(zna,nrow,ncol,nnrow,nncol)
      implicit none
!
!!	DO ROW OPERATIONS ON MATRIX ZNA TO BRING IT TO UPPER TRIANGULAR FORM
!!	DOUBLE PRECISION NUMBERS
!
      integer nnrow,nncol
      double precision zna(nnrow,nncol),zntemp
      integer nrow,ncol,nout,j,k,l,ncmp
!
!	TRANSFORM MATRIX BY ROWS
      nout=0
   do j=1,nrow
5     if(j+nout.gt.ncol)return
!	PUT NON-ZERO ELEMENTS ON DIAGONAL
      if(ncmp(zna(j,j+nout)).ne.0)goto 4
      do  k=j+1,nrow
      if(ncmp(zna(k,j+nout)).ne.0)goto 3
      enddo
!	IF COLUMN HAS ALL ZEROES, GO TO NEXT COLUMN
      nout=nout+1
      goto 5
!	INTERCHANGE ROWS
3     do l=1,ncol
      zntemp=zna(j,l)
      zna(j,l)=zna(k,l)
      zna(k,l)=zntemp
      enddo
!	NORMALIZE ELEMENTS IN ROW
4     zntemp=zna(j,j+nout)
      do k=j+nout,ncol
      zna(j,k)=zna(j,k)/zntemp
      enddo
!	REDUCE ALL OTHER ROWS
      do  k=1,nrow
        if(k.eq.j) cycle
        if(ncmp(zna(k,j+nout)).eq.0) cycle
        do  l=1,ncol
          if(l.eq.j+nout) cycle
          zna(k,l)=zna(k,l)-zna(j,l)*zna(k,j+nout)
        enddo
        zna(k,j+nout)=0.
      enddo
  enddo

      end subroutine xrowop2
!-------------------------------------------------------------------------------
      subroutine getstar(kvec,primlatt,narms,kvecstar,kvecop)
!! find atom in data base
!! arguments:
!!     kvec(i) (input), ith dimensionless component of k vector
!!     narms (output), number of star vectors associated with kvec
!!     kvecstar(3,1:narms), all the stars of kvec
!!     kvecop(i), the symmetry operation number for the star vector i
      use atoms_force_constants
!      use force_constants_module
      implicit none
      integer, intent(out):: narms,kvecop(48)
      real(8), intent(in) :: kvec(3),primlatt(3,3)
      real(8), intent(out):: kvecstar(3,48)
      integer i,j,k,n,ncmp
      real(8) v(3),v2(3),v3(3),kvecstarp(3,48)

      narms=0
!      print*,'lattpgcount=',lattpgcount
      iloop: do i=1,lattpgcount
! apply symmetry operation to k to get v=kstar
        call xvmlt(op_kmatrix(1,1,i),kvec,v,3,3,3)
! find the reduced coordinates of v and store in v2
        call xmatmlt(v,primlatt,v2,1,3,3,1,3,1)

! now check if v2 - any_previous_v2 is integer (differ by a G vector)
! if so, skip; else store this v as a new star vector
        do j=1,narms
        ! subtract previous_v2(=kvecstarp) from v2; result is v3
          call xvsub(v2,kvecstarp(1,j),v3,3)
          do k=1,3
            n=nint(v3(k))
            if(ncmp(v3(k)-n).ne.0)exit
            if(k.eq.3)cycle iloop  ! goto next sym_op iff v3=integer
          enddo
        enddo
        narms=narms+1
        kvecstar(1:3,narms)=v(1:3)
        kvecstarp(1:3,narms)=v2(1:3)
        kvecop(narms)=i
      enddo iloop

      end subroutine getstar
!****************************************************************************
       subroutine collect_force_constants(nrank,nshell,   &
     &     maxrank,maxterms,maxtermsindep,maxtermszero,maxgroups,  &
     &     ngroups,amat,  &
     &     ntermsindep,iatomtermindep,ixyztermindep,    &
     &     nterm6,     iatomterm,     ixyzterm,    &
     &     ntermszero, iatomtermzero, ixyztermzero, &
     &     ierz,iert,ieri,ierg)

!     subroutine collect_force_constants(nrank,nshell,ngroups,  &
!     &     ntermsindep,iatomtermindep,ixyztermindep,nterm6,  &
!     &     iatomterm,ixyzterm,amat,ntermszero,iatomtermzero,  &
!     &     ixyztermzero,maxrank,maxterms,maxtermsindep,maxtermszero,  &
!     &     maxgroups,ierz,iert,ieri,ierg)
!! get force constants for a given rank out to a given nearest neighbor shell
! arguments:
!     nrank (input), order or derivative
!     nshell (input), shell of nearest-neighbor atoms to be included
!     ngroups (output), number of groups of terms.  Within each group, the
!          terms are related by symmetry.
!     For the ith group:
!     ntermsindep(i) (output), number of independent terms
!     iatomtermindep(k,j,i) (output), atom at kth location in denominator of
!          jth independent term (k=1:rank)
!     ixyztermindep(k,j,i) (output), coordinate at kth location in denomonitor
!          of jth independent term
!     nterm6(i) (output), total number of terms
!     iatomterm(k,j,i) (output), atom at kth location in denominators of
!          jth term
!     ixyzterm(k,j,i) (output), coordinate at kth location in denomonitor of
!          jth term
!     amat(k,j,i) (output), coefficient of jth independent term in equation
!          for kth term
!     ntermszero (output), number of zero terms
!     iatomtermzero(k,j) (output), atom at kth location in denominator of
!          jth zero term
!     ixyztermzero(k,j) (output), coordinate at kth location in denomonitor
!          of jth zero term
!     maxrank (input), number of rows in arrays iatomtermindep, ixyztermindep,
!          iatomterm, and ixyzterm
!     maxterms (input), number of rows in array amat and number of columns
!          in arrays iatomterm and ixyzterm
!     maxtermsindep (input), number of columns
!          in arrays amat, iatomtermindep, and ixyztermindep
!     maxtermszero (input), number of columns
!          in arrays iatomtermzero and ixyztermzero
!     maxgroups (input), maximum number of groups

      use atoms_force_constants
!      use force_constants_module
      use ios , only : ulog
      implicit none
      integer, intent(in) :: nrank,maxrank,maxterms,maxtermsindep,maxgroups,  &
     &     maxtermszero,nshell(natom_prim_cell)
!  make sure natom_prim_cell < 20
      integer, intent(out):: ierz,iert,ieri,ierg,ngroups,ntermszero,nterm6(maxgroups), &
     &     ntermsindep(maxgroups),iatomtermindep(maxrank,maxterms,maxgroups),  &
     &     ixyztermindep(maxrank,maxterms,maxgroups), &
     &     iatomterm(maxrank,maxterms,maxgroups),ixyzterm(maxrank,maxterms,maxgroups),  &
     &     iatomtermzero(maxrank,maxtermszero),ixyztermzero(maxrank,maxtermszero)
      real(8), intent(out):: amat(maxterms,maxtermsindep,maxgroups)
      integer i,j,k,m,n0,ncount,iatom,iatom0,iatomd(nrank),ixyz,ngroupsave,  &
     &     ixyzd(nrank),icell(3),iatomd2(nrank),ixyzd2(nrank),ntermszerosave
      logical firsttime

      ierz=0; iert=0; ierg=0; ieri=0

      write(*,*)' INPUTS OF COLLECT_FORCE_CONSTANTS '
      write(*,*)' nrank=',nrank
      do i=1,natom_prim_cell
         write(*,*)' iatom,nshell(iatom)=',i,nshell(i)
      enddo
      write(*,*)' mxtrm,mxtrmindp,maxtermszero,maxgroups=',maxterms,maxtermsindep, &
      &           maxtermszero,maxgroups 

! check if input values are valid
      if(nrank.gt.maxrank)then
        write(6   ,*)'Error in collect_force_constants: nrank \= maxrank',nrank,maxrank
        write(ulog,*)'Error in collect_force_constants: nrank \= maxrank',nrank,maxrank
        stop
      endif

      ncount=0
      ngroups=0
      amat=0
      ntermszero=0
! do each atom in unit cell
      iatom0loop: do iatom0=1,natom_prim_cell
! find first nearest neighbor in list
        do iatom=1,natoms
          if(iatomneighbor(iatom0,iatom).le.nshell(iatom0))exit
          if(iatom.eq.natoms)then
            write(6,*)'Error in collect_force_constants: first '  &
     &           //'nearest neighbor not found'
            stop
          endif
        enddo
! try each atom for each position in denominator
! begin with first nearest neighbor
        iatomd(1:nrank)=iatom
        iatomd(nrank)=iatom-1
        iatomd(1)=iatom0
        firsttime=.true.
        nextatomloop: do while(.true.)
! next set of atoms
          if(firsttime)then
            firsttime=.false.
          else
            if(nrank.eq.1)exit nextatomloop
          endif
! try each position in denominator
          iloop: do i=nrank,2,-1
! try each atom in that position
            jloop: do j=iatomd(i)+1,natoms
! nearest neighbor to atom in cell at origin?
              if(iatomneighbor(iatom0,j).le.nshell(iatom0))then
! nearest neighbor to all other atoms in denominator?
                do k=2,i-1
                  do m=1,3
                    icell(m)=iatomcell(m,iatomd(k))-iatomcell(m,j)
                  enddo
                  call findatom(icell,iatomcell0(iatomd(k)),m)
                  if(m.eq.0)cycle jloop
                  if(iatomneighbor(iatomcell0(j),m).gt.nshell(iatom0)) cycle jloop
                enddo
! yes: put atom in denominator and get terms
                iatomd(i:nrank)=j
                exit iloop
              endif
            enddo jloop
! done: try next atom in cell at origin
            if(i.eq.2)cycle iatom0loop
          enddo iloop
! check if found yet
          ixyzd=1
          call unique_force_constant(nrank,iatomd,ixyzd,iatomd2,ixyzd2)
          do k=1,ngroups
          iloop2: do i=1,nterm6(k)
            do j=1,nrank
              if(iatomterm(j,i,k).ne.iatomd2(j))cycle iloop2
            enddo
! found it: try next set of atoms
            cycle nextatomloop
          enddo iloop2
          enddo
          iloop4: do i=1,ntermszero
            do j=1,nrank
              if(iatomtermzero(j,i).ne.iatomd2(j))cycle iloop4
            enddo
! found it: try next set of atoms
            cycle nextatomloop
          enddo iloop4
! did not find it: generate terms
! all possible sets of coordinates
          ngroupsave=ngroups
          ntermszerosave=ntermszero
          ixyzloop: do ixyz=1,3**nrank
            m=ixyz-1
            do i=1,nrank
              ixyzd(i)=mod(m,3)+1
              m=m/3
            enddo
! check if found yet
            call unique_force_constant(nrank,iatomd,ixyzd, iatomd2,ixyzd2)
            do k=ngroupsave+1,ngroups
            iloop3: do i=1,nterm6(k)
              do j=1,nrank
                if(iatomterm(j,i,k).ne.iatomd2(j))cycle iloop3
                if(ixyzterm(j,i,k).ne.ixyzd2(j))cycle iloop3
              enddo
! found it: try next set of coordinates
              cycle ixyzloop
            enddo iloop3
            enddo
           iloop5: do i=ntermszerosave+1,ntermszero
             do j=1,nrank
               if(iatomtermzero(j,i).ne.iatomd2(j))cycle iloop5
               if(ixyztermzero(j,i).ne.ixyzd2(j))cycle iloop5
             enddo
! found it: try next set of coordinates
             cycle ixyzloop
           enddo iloop5
! did not find it: generate terms
            ngroups=ngroups+1

            if(ngroups.ge.maxgroups)then
              write(6,*)' maxgroups is too small ',maxgroups
              ierg=ierg+1
              return
            endif

            if (ntermszero.ge.maxtermszero) then
              write(*,3)' mxt0 =< nt0=',maxtermszero,ntermszero
              write(ulog,3)' mxt0 =< nt0=',maxtermszero,ntermszero
              ierz=5
              return
            endif


            call force_constants(nrank,amat(1,1,ngroups),iatomd2,ixyzd2,  &
     &           ntermsindep(ngroups),iatomtermindep(1,1,ngroups),ixyztermindep(1,1,ngroups), &
     &           nterm6(ngroups),iatomterm(1,1,ngroups),ixyzterm(1,1,ngroups),  &
     &           n0,iatomtermzero(1,ntermszero+1),ixyztermzero(1,ntermszero+1),  &
     &           maxrank,maxterms,maxtermsindep,maxtermszero-ntermszero,  &
     &           ierz,iert,ieri) !,ierg)

            write(*,3)'OUTPUT OF FORCE_CONSTANTS FOR RANK/ngroups=',nrank,ngroups
            write(*,3)'OUTPUT OF FORCE_CONSTANTS: ntermsindep=',ntermsindep(1:ngroups)
            write(*,3)'OUTPUT OF FORCE_CONSTANTS: nterms     =',nterm6(1:ngroups)
            write(*,3)'OUTPUT OF FORCE_CONSTANTS: ntermszero =',ntermszero
!            write(ulog,3)'OUTPUT OF FORCE_CONSTANTS FOR RANK/ngroups=',nrank,ngroups
!            write(ulog,3)'OUTPUT OF FORCE_CONSTANTS: ntermsindep=',ntermsindep(1:ngroups)
!            write(ulog,3)'OUTPUT OF FORCE_CONSTANTS: nterms     =',nterm6(1:ngroups)
!            write(ulog,3)'OUTPUT OF FORCE_CONSTANTS: ntermszero =',ntermszero

            if (ierz.ne.0) then
              write(*,3)' maxtermszero is too small,nt0=',maxtermszero,ierz,ntermszero
              return
            endif
            if (iert.ne.0) then
              write(*,3)' maxterms is too small,nterms',maxterms,iert,sum(nterm6)
              return
            endif
            if (ieri.ne.0) then
              write(*,3)' maxtermsindep is too small,ntindep',maxtermsindep,ieri,sum(ntermsindep)
              return
            endif
            ntermszero=ntermszero+n0
! no nonzero terms in group
            if(nterm6(ngroups).eq.0)then
              ngroups=ngroups-1
              ncount=ncount+1
            endif
! next set of coordinates
          enddo ixyzloop
! next set of atoms
        enddo nextatomloop
      enddo iatom0loop

      write(ulog,*)' ngroups=',ngroups
      write(ulog,3)' ntermsindep        (1:ngroups)=', ntermsindep(1:ngroups)
      write(ulog,3)' nterm6             (1:ngroups)=', nterm6(1:ngroups)
      write(ulog,3)' iatomtrmindep(rnk,1,1:ngroups)=', iatomtermindep(nrank,1,1:ngroups)
      write(ulog,3)' ixyztermindep(rnk,1,1:ngroups)=', ixyztermindep(nrank,1,1:ngroups)
      write(ulog,3)' iatomterm    (rnk,1,1:ngroups)=', iatomterm(nrank,1,1:ngroups)
      write(ulog,3)' ixyzterm     (rnk,1,1:ngroups)=', ixyzterm(nrank,1,1:ngroups)
      write(ulog,3)' ntermszero =', ntermszero

      ierz=0 ; iert=0 ; ieri=0 ; ierg=0

 3    format(a,99(i4))
      end subroutine collect_force_constants
!****************************************************************************
