!===============================================================================
 subroutine setup_maps
!! maps the output of collect_fcs to the new arrays for each rank,
!! with the proper (allocatable) dimensions
 use svd_stuff
 use params
 use atoms_force_constants
 use ios
 use geometry
 implicit none
 character lineout*80,frmt*90
 integer i0,j,i,n,m,mxzero,mx,mxi,l,rnk,ichkloop   &
 &      ,ntermszero,nd,ntindp,ierz,iert,ieri,ierg,mzero(maxrank),mxgrps
 integer, allocatable:: iatmtermzero(:,:),ixyztermzero(:,:)
! type arraysize
!    integer rnk
! end type arraysize

 write(*,*) 'Entering setup_maps routine with read nshells as:'
 call write_out(ulog,' #include_FCs ',include_fc)
 call write_out(ulog,' #   maxterms ',maxterms)
 call write_out(ulog,' #r=1 nshells ',nshells(1,:))
 call write_out(ulog,' #r=2 nshells ',nshells(2,:))
 call write_out(ulog,' #r=3 nshells ',nshells(3,:))
 call write_out(ulog,' #r=4 nshells ',nshells(4,:))
 call write_out(ulog,' #r=5 nshells ',nshells(5,:))
 call write_out(ulog,' #r=6 nshells ',nshells(6,:))
 call write_out(ulog,' #r=7 nshells ',nshells(7,:))
 call write_out(ulog,' #r=8 nshells ',nshells(8,:))
 write(ulog,*)'# trans,rot,huang,enforce_inv=',itrans,irot,ihuang,enforce_inv

!! K1 2/12/23
!   mzero(1:4) = maxtermzero(1:4)
    mzero = maxtermzero  ! 3/20/24
!! K1 2/12/23

 nd = sum(include_fc)
 if (nd.gt.4) then
    write(ulog,*)' SETUP_MAPS: nd>4, check your include_fc array ',nd
    write(*   ,*)' SETUP_MAPS: nd>4, check your include_fc array ',nd
!    stop
 endif


 write(ulog,*)' SETUP_MAPS: ******************************************'

! loop over ranks that are supposed to be included
! -----------------------------------
 rankloop: do rnk=1,maxrank !,1,-1
  if ( include_fc(rnk) .ne. 0 ) then

    ierz=1; iert=1;  ieri=1; ierg=0
    mxgrps = maxgroups(rnk)
    mx     = maxterms(rnk)
    mxi    = maxtermsindep(rnk)
    mxzero = maxtermzero(rnk)

    write(ulog,4)' ******************** FOR rank *****************= ',rnk
    write(*   ,4)' ***** FOR rank = ',rnk

    ichkloop=0
! start from a small maxterms and maxtermszero and increase by 100 or more if too small
  checkloop: do while (ierg+ierz+iert+ieri.ne.0)
     ichkloop=ichkloop+1
     write(*   ,4) 'Checkloop: icheckloop=',ichkloop
     write(*   ,4) 'maxtrmzero,maxgrp=',mxzero,mxgrps
!     write(*   ,4) 'ntermsindep=',ntermsindep
!     write(*   ,4) 'nterms     =',nterms
     write(ulog,4) 'iteration,maxtrmzero,maxgrp,maxternindep,maxterms=',ichkloop,mxzero,mxgrps,mxi,mx
!     write(ulog,4) 'ntermsindp=',ntia
!     write(ulog,4) 'nterms    =',nta
!     mx=maxval(nta) ; mxi=maxval(ntia)
!    allocate(mapmat(mx,mxi,maxgroups) &
!    &     ,iatmtrm(rnk,mx,maxgroups),ixyztrm(rnk,mx,maxgroups)  &
!    &     ,iatmtermindp(rnk,mx,maxgroups),ixyztermindp(rnk,mx,maxgroups)  )
!  if(ichkloop .eq.1) then
     if(.not.allocated(nterm))       allocate(nterm(mxgrps))
     if(.not.allocated(ntermsindep)) allocate(ntermsindep(mxgrps))
     if(.not.allocated(mapmat))      allocate( mapmat(mx,mxi,mxgrps)); mapmat=0
     if(.not.allocated(iatmtrm))     allocate( iatmtrm(rnk,mx,mxgrps));iatmtrm=0
     if(.not.allocated(ixyztrm))     allocate( ixyztrm(rnk,mx,mxgrps));ixyztrm=0
     if(.not.allocated(iatmtermindp)) allocate( iatmtermindp(rnk,mx,mxgrps));iatmtermindp=0
     if(.not.allocated(ixyztermindp)) allocate( ixyztermindp(rnk,mx,mxgrps));ixyztermindp=0
     if(.not.allocated(iatmtermzero)) allocate( iatmtermzero(rnk,mxzero));iatmtermzero=0
     if(.not.allocated(ixyztermzero)) allocate( ixyztermzero(rnk,mxzero));ixyztermzero=0
! k1 2/12/23
! k1 3/5/23
!     if(.not.allocated(iatmtermzero)) allocate( iatmtermzero(rnk,ntermszero))
!     if(.not.allocated( ixyztermzero)) allocate(  ixyztermzero(rnk,ntermszero))
!     if(.not.allocated(iatmtermzero)) allocate(iatmtermzero(rnk,mzero(rnk)))
!     if(.not.allocated(ixyztermzero )) allocate( ixyztermzero(rnk,mzero(rnk)))
! k2 2/12/23
!     write(ulog,4)' calling collect_force_constants ier:ztig=',ierz,iert,ieri,ierg

     call collect_force_constants(rnk,nshells(rnk,:),  &
     &       rnk,mx,mxi,mxzero,mxgrps,   &
     &       ngroups(rnk),mapmat,  &
     &       ntermsindep,iatmtermindp,ixyztermindp,  &
     &       nterm,iatmtrm,ixyztrm,   &
     &       ntermszero,iatmtermzero,ixyztermzero,   &
     &       ierz,iert,ieri,ierg)

!     write(ulog,4)' collect_force_constants called with ier:ztig=',ierz,iert,ieri,ierg
!     write(ulog,*)' ntermszero,ngroups,ntindep,nt=',ntermszero,ngroups,sum(ntermsindep(1:ngroups(rnk))) &
!     &   ,sum(nterm(1:ngroups(rnk)))

 4 format(a,9(i5))

     if (ierz.ne.0) then
   !      write(ulog,*)' before mxtermszero, nterms0=',mxzero ,ntermszero
!       if(ichkloop.eq.1) then
          mxzero = mxzero*2 !+ 50*rnk !,ntermszero)
!       else
!          mxzero=mxzero+ntermszero
!       endif
         if (allocated(iatmtermzero)) deallocate(iatmtermzero)
         if (allocated(ixyztermzero)) deallocate(ixyztermzero)
     endif
     write(ulog,*)' mxtermszero  is now ',mxzero
     if (iert.ne.0) then
    !     write(ulog,*)' before mxterms,sum(nterms)=',mx,sum(nterm)
!       if(ichkloop.eq.1) then
         mx = mx*2 !+100*natom_prim_cell*nshells(rnk,1)**(rnk-1) !,sum(nterm(1:ngroups(rnk))))
!       else
!          mx=mx+sum(nterm(1:ngroups(rnk)))
!       endif
!         write(ulog,*)' ntermall=',sum(nterm)
         if (allocated(mapmat)) deallocate(mapmat)
         if (allocated(iatmtrm)) deallocate(iatmtrm)
         if (allocated(ixyztrm)) deallocate(ixyztrm)
         if (allocated(iatmtermindp)) deallocate(iatmtermindp)
         if (allocated(ixyztermindp)) deallocate(ixyztermindp)
     endif
     write(ulog,*)' maxterms      is now ',mx
   !  write(ulog,*)' before mxtermsindep,sum(ntermsindep)=',mxi,sum(ntermsindep)
     if (ieri.ne.0) then
!       if(ichkloop.eq.1) then
         mxi= mxi*2 !+5*natom_prim_cell*nshells(rnk,1)**(rnk-1) !,sum(ntermsindep(1:ngroups(rnk))))
!       else
!          mxi=mxi+sum(ntermsindep(1:ngroups(rnk)))
!       endif
         if (allocated(mapmat)) deallocate(mapmat)
!         if (allocated(iatmtrm)) deallocate(iatmtrm)
!         if (allocated(ixyztrm)) deallocate(ixyztrm)
!         if (allocated(iatmtermindp)) deallocate(iatmtermindp)
!         if (allocated(ixyztermindp)) deallocate(ixyztermindp)
     endif
     write(ulog,*)' maxtermsindep is now ',mxi
     if (ierg.ne.0) then
    !     write(ulog,*)' before maxgroups, ngroups=',mxgrps,ngroups(rnk)
!       if(ichkloop.eq.1) then
         mxgrps= mxgrps*2 !+2**(rnk-1)*natom_prim_cell*nshells(rnk,1) !,ngroups(rnk))
!       else
!         mxgrps=ngroups(rnk)
!       endif
         if (allocated(nterm)) deallocate(nterm)
         if (allocated(ntermsindep)) deallocate(ntermsindep)
         if (allocated(mapmat)) deallocate(mapmat)
         if (allocated(iatmtrm)) deallocate(iatmtrm)
         if (allocated(ixyztrm)) deallocate(ixyztrm)
         if (allocated(iatmtermindp)) deallocate(iatmtermindp)
         if (allocated(ixyztermindp)) deallocate(ixyztermindp)
      !   if (allocated(nterm)) deallocate(ixyztermindp)
     endif
     write(ulog,*)' maxgroups     is now ',mxgrps

  enddo checkloop

     write(frmt,'(a,i3,a)') '(a,i1,a,',natom_prim_cell,'(1x,i2),a,i8,a,i5,a)'   

!    write(umap,'(a,i1,a,20(1x,i2),a,i8,a,i5,a)')   &
     write(umap,frmt)    &
&     '************ rank ',rnk,', shell ',nshells(rnk,1:natom_prim_cell),', groups ',ngroups(rnk)
     if(ngroups(rnk).gt.0) then

        allocate(map(rnk)%gr   (ngroups(rnk)),    &
        &        map(rnk)%nt   (ngroups(rnk)),    &
        &        map(rnk)%ntind(ngroups(rnk)) )
        map(rnk)%ngr      = ngroups(rnk)
        map(rnk)%nt   (:) = nterm      (1:ngroups(rnk))
        map(rnk)%ntind(:) = ntermsindep(1:ngroups(rnk))
        map(rnk)%ntotind = sum(map(rnk)%ntind(:))
        map(rnk)%ntot    = sum(map(rnk)%nt   (:))
        allocate(map(rnk)%err(map(rnk)%ntotind))
        write(ulog,*) 'ngroups(rnk)=',ngroups(rnk)

        do i=1,ngroups(rnk)
!            if (rnk .ge. 2) then
!            ! take the first indep term in that group
!               i0=iatmtermindp(1,1,i)  ! map(2)%gr(i)%iat(1,1)
!               j =iatmtermindp(2,1,i)  ! map(2)%gr(i)%iat(2,1)  ! this j is for atompos
!               r12=length(atompos(:,j)-atompos(:,i0))
!            endif
!           write(ulog,9) 'i, ntindep(i)=',i,map(rnk)%ntind(:)
!           write(ulog,9) 'i,      nt(i)=',i,map(rnk)%nt(:)

           allocate(map(rnk)%gr(i)%mat (nterm(i),ntermsindep(i)),   &
       &         map(rnk)%gr(i)%iat (rnk,nterm(i)),                 &
       &         map(rnk)%gr(i)%ixyz(rnk,nterm(i)),                 &
       &         map(rnk)%gr(i)%iatind (rnk,ntermsindep(i))  ,      &
       &         map(rnk)%gr(i)%ixyzind(rnk,ntermsindep(i))  )
        enddo


        write(umap,*)'rank=',rnk,' groups ',map(rnk)%ngr,' indepterms=',map(rnk)%ntotind,' all terms=',map(rnk)%ntot
        write(umap,'(a,i2,a,99(i5))')'rank=',rnk,' indepterms=',map(rnk)%ntind(:)
        write(umap,'(a,i2,a,99(i5))')'rank=',rnk,' all  terms=',map(rnk)%nt(:)

        write(ulog,*)' allocation of MAP done! '
        map(rnk)%err = 'C'

! copy from temp array into real arrays with actual size
        do i=1,map(rnk)%ngr    !  = ngroups(rnk)
          map(rnk)%gr(i)%mat (:,:)   = mapmat(1:nterm(i),1:ntermsindep(i),i)
          map(rnk)%gr(i)%iat (:,:)   = iatmtrm(1:rnk,1:nterm(i),i)
          map(rnk)%gr(i)%ixyz(:,:)   = ixyztrm (1:rnk,1:nterm(i),i)
          map(rnk)%gr(i)%iatind (:,:)= iatmtermindp(1:rnk,1:ntermsindep(i),i)
          map(rnk)%gr(i)%ixyzind(:,:)= ixyztermindp (1:rnk,1:ntermsindep(i),i)
        enddo

        m=0
        do i=0,nshells(2,1)
           m = m + atom0(1)%shells(i)%no_of_neighbors
        enddo
        write(ulog,*)'For rank=2 atm 1 has ',m,'inclusive nghbrs within ',nshells(rnk,1),' shell'
        inv_constraints = inv_constraints + natom_prim_cell*12*m
        write(ulog,*)' rank=',rnk,' # of groups=',map(rnk)%ngr
        write(ulog,9)' ntind=',map(rnk)%ntind(:)
        write(ulog,9)' nterm=',map(rnk)%nt(:)
        write(ulog,*)' Cumulative invce_cnstrnts for this shell estimated to be ', inv_constraints

     endif

     deallocate(nterm,iatmtrm,ixyztrm,mapmat,  &
  &     iatmtermzero,ixyztermzero,ntermsindep,iatmtermindp,ixyztermindp)

     write(ulog,*)'iattrm,ixyztrm,mapmat, iattrm0,ixyztrm0,iattrmindp,ixyztrmindp deallocated'

  else
     map(rnk)%ngr=0
     write(ulog,*) 'Allocating map(rank=',rnk,')'
     allocate(map(rnk)%nt(1),map(rnk)%ntind(1) )
     map(rnk)%nt(1)=0
     map(rnk)%ntind(1)=0
  endif
 enddo rankloop

 write(*,*) 'SETUP_MAPS: Exited the main rank loop , and going to write by calling ustring'

 do rnk=1,maxrank
  if ( include_fc(rnk) .ne. 0 ) then
    if(ngroups(rnk).gt.0) then
      n=7+4*rnk
!   do i=1,min(5000,nterm(maxgroups))   ! for now write out at most 5000 terms
      do j=1,map(rnk)%ngr
         write(umap,'(a,2i4)')' ============== RANK, GROUP ===============',rnk,j
         write(umap,'(a,i4)')' ---- # OF INDEPENDENT TERMS ------ ',map(rnk)%ntind(j)
         do i=1,map(rnk)%ntind(j)
            m=0
            lineout=' '
!           write(ulog,*) 'iat=',map(rnk)%gr(j)%iatind(1:rnk,i)
!           write(ulog,*) 'ixy=',map(rnk)%gr(j)%ixyzind(1:rnk,i)
            call ustring(m,lineout,rnk,map(rnk)%gr(j)%iatind(:,i),map(rnk)%gr(j)%ixyzind(:,i))
            i0=map(rnk)%gr(j)%iatind(1,i)
            if(rnk.ge.2) then
               l=map(rnk)%gr(j)%iatind(2,i)
               write(umap,'(4x,a,f11.5)') lineout(1:m)//'       r12=',length(atompos(:,l)-atompos(:,i0))
            endif
         enddo
         write(umap,'(a,i4)')' ---- # OF ALL TERMS, MAT   ------- ',map(rnk)%nt(j)
         ntindp = map(rnk)%ntind(j)
         do i=1,map(rnk)%nt(j)
            m=0
            lineout=' '
!           write(ulog,*) 'iat=',map(rnk)%gr(j)%iat(1:rnk,i)
!           write(ulog,*) 'ixy=',map(rnk)%gr(j)%ixyz(1:rnk,i)
            call ustring(m,lineout,rnk,map(rnk)%gr(j)%iat(:,i),map(rnk)%gr(j)%ixyz(:,i))
            write(umap,'(4x,a,90f7.3)') lineout(1:25),map(rnk)%gr(j)%mat(i,1:ntindp)
         enddo
      enddo
    endif
  endif
 enddo !rankloop

 nindepfc = 0
 do rnk=1,maxrank
    if ( include_fc(rnk) .eq. 1 ) then ! exclude them if they already exist and will be read
! get the total number of independent and full terms for each rank
       call get_dim(map(rnk),map(rnk)%ntotind,map(rnk)%ntot)
       write(ulog,*)'SETUP_MAPS: RANK, NDINDEP, NDFULL=',rnk,map(rnk)%ntotind,map(rnk)%ntot

       nindepfc = nindepfc + map(rnk)%ntotind

    endif
 enddo
 write(ulog,*)'END OF SETUP_MAPS: nindepfc before imposing new cutoff=',nindepfc
 write(*,*) 'exiting setup_maps routine'

9 format(a,999(i4))
44 format(a,3i5,2x,f11.5)
 end subroutine setup_maps
!===============================================================================
 subroutine include_constraints_remove_zeros
! now that we have the 3 matrices atransl,arot,ahuang and aforce, here we decide
! how to put them together to do the SVD
! outputs are amat and bmat and their dimensions:dim_al,dim_ac, which go to SVD
 use svd_stuff
 use params
 use ios
 use linalgb
 use constants, only : r15
 implicit none
 integer nlout,ierr
 real(r15), allocatable:: bout(:),aout(:,:)
! real(r15), allocatable:: a11(:,:),b1(:),a12(:,:),a11i(:,:),b2(:)

!  dim_al = force_constraints
!  allocate(amat(dim_al,dim_ac),bmat(dim_al))
!  amat=aforce ; bmat=bforce
!  deallocate(aforce,bforce,stat=ierr)
!  if(ierr.ne.0) then
!     write(*,*)'translation deallocation error in remove_zeros'
!     stop
!  endif

! start with translational invariance matrix
   dim_al=transl_constraints+rot_constraints+huang_constraints+force_constraints
   dim_ac = nindepfc
   allocate(amat(dim_al,dim_ac),bmat(dim_al))


!-----------------------------------------------------------------------
! in this case, invariances are obtained in a rms sense and are not exact.
! in principle zero lines have already been removed from the invariance part
! of amat, and the call to remove zeros in not necessary!


   dim_al=0
   if (itrans .ne. 0) then
      allocate(aout(transl_constraints,dim_ac),bout(transl_constraints))
      call remove_zeros(transl_constraints,dim_ac,atransl,btransl,  &
      &        nlout,aout,bout)

      amat(dim_al+1:dim_al+nlout,1:dim_ac)=aout(1:nlout,1:dim_ac)
      bmat(dim_al+1:dim_al+nlout)=bout(1:nlout)
      dim_al= dim_al+ nlout
      write(ulog,*)'INCLUDE_CONSTRAINTS :c(transl),dim_hom(A)=',nlout,dim_al
      deallocate(aout,bout,stat=ierr)
      if(ierr.ne.0) then
         write(*,*)'rotation deallocation error in remove_zeros'
         stop
      endif
      transl_constraints=nlout

      write(ulog,*)'After remove_zeros from transl: NEW # of lines=',nlout
      write(*   ,*)'After remove_zeros from transl: NEW # of lines=',nlout

   endif
!     write(*,*)' going to call append_array 0'
!     call append_array(amat,aout,amat)
!     call append_array(bmat,bout,bmat)

!      deallocate(atransl,btransl)
!     allocate(amat(dim_al,dim_ac),bmat(dim_al))
!     amat(1:dim_al,1:dim_ac)=aout(1:nlout,1:dim_ac)
!     bmat(1:dim_al)=bout(1:nlout)
!     transl_constraints=nlout
!     write(ulog,*)'INCLUDE_CONSTRAINTS :c(transl),dim_hom(A)=',nlout,dim_al
!     write(*   ,*)'INCLUDE_CONSTRAINTS :c(transl),dim_hom(A)=',nlout,dim_al
!     deallocate(aout,bout,stat=ierr)
!     if(ierr.ne.0) then
!        write(*,*)'translation deallocation error in remove_zeros'
!        stop
!     endif
!  endif
   if (irot   .ne. 0) then
      allocate(aout(rot_constraints,dim_ac),bout(rot_constraints))
      call remove_zeros(rot_constraints,dim_ac,arot,brot,  &
      &        nlout,aout,bout)


!write(ulog,*)'After remzeros from rot: NEW # of lines=',nlout
!      deallocate(arot,brot)

! amat=reshape(transpose(amat),shape=(/size(amat,2),size(amat,1)+size(aout,1)/),pad=    transpose(aout),order=(/1,2/))
! amat=transpose(amat)
! bmat=reshape(bmat,shape=(/size(bmat)+size(bout)/),pad=bout)

!     write(*,*)' going to call append_array 2'
!     call append_array(amat,aout,amat)
!     call append_array(bmat,bout,bmat)

      amat(dim_al+1:dim_al+nlout,1:dim_ac)=aout(1:nlout,1:dim_ac)
      bmat(dim_al+1:dim_al+nlout)=bout(1:nlout)
      dim_al= dim_al + nlout
      write(ulog,*)'INCLUDE_CONSTRAINTS :c(rotatn),dim_hom(A)=',nlout,dim_al
      deallocate(aout,bout,stat=ierr)
      if(ierr.ne.0) then
         write(*,*)'rotation deallocation error in remove_zeros'
         stop
      endif
      rot_constraints=nlout
   endif

   if (ihuang .ne. 0) then
      allocate(aout(huang_constraints,dim_ac),bout(huang_constraints))
      call remove_zeros(huang_constraints,dim_ac,ahuang,bhuang,   &
      &        nlout,aout,bout)
      amat(dim_al+1:dim_al+nlout,1:dim_ac)=aout(1:nlout,1:dim_ac)
      bmat(dim_al+1:dim_al+nlout)=bout(1:nlout)
      dim_al= dim_al+  nlout
      write(ulog,*)'INCLUDE_CONSTRAINTS :c(huang ),dim_hom(A)=',nlout,dim_al
      deallocate(aout,bout,stat=ierr)
      if(ierr.ne.0) then
         write(*,*)'Huang deallocation error in remove_zeros'
         stop
      endif

 write(ulog,*)'After remzeros from Huang: NEW # of lines=',nlout
!     deallocate(ahuang,bhuang)

! amat=reshape(transpose(amat),shape=(/size(amat,2),size(amat,1)+size(aout,1)/),pad=    transpose(aout),order=(/1,2/))
! amat=transpose(amat)
! bmat=reshape(bmat,shape=(/size(bmat)+size(bout)/),pad=bout)

!     write(*,*)' going to call append_array 1'
!     call append_array(amat,aout,amat)
!     call append_array(bmat,bout,bmat)

      huang_constraints=nlout
   endif
!

   inv_constraints=dim_al
   write(ulog,3)' REMOVE_ZEROS: # of invce cnstrnts=',transl_constraints,rot_constraints,huang_constraints
   write(ulog,*)' Total # of invariance constraints=',inv_constraints

      allocate(aout(force_constraints,dim_ac),bout(force_constraints))
   ! no need to remove zeros for forces
      call remove_zeros(force_constraints,dim_ac,aforce,bforce,   &
      &        nlout,aout,bout)
!     deallocate(aforce,bforce)


!     write(*,*)' force_constraints, before appending ',force_constraints
!     write(*,*)' force_constraints, after  appending ',force_constraints

!     write(*,*)' number of lines before appending forces=',size(amat,1),size(aout,1)
!     write(*,*)' going to call append_array 3'
!     call append_array(amat,aout,amat)
!     call append_array(bmat,bout,bmat)
!     call append_array(amat,aforce,amat)
!     call append_array(bmat,bforce,bmat)
      amat(dim_al+1:dim_al+nlout,1:dim_ac)=aout(1:nlout,1:dim_ac)
      bmat(dim_al+1:dim_al+nlout)=bout(1:nlout)
      dim_al= dim_al+  nlout
      write(ulog,*)'INCLUDE_CONSTRAINTS :c(force ),dim_tot(A)=',nlout,dim_al
      force_constraints=nlout
      deallocate(aout,bout,stat=ierr)
      if(ierr.ne.0) then
         write(*,*)'force deallocation error in remove_zeros'
         stop
      endif


   write(ulog,*)' size(a) nlines,ncolumn=(dim(a1d))=',dim_al,dim_ac
   write(*,*)' Exiting include_constraints: size(a) nlines,ncolumn=(dim(a1d))=',dim_al,dim_ac
   deallocate (aforce,bforce,stat=ierr)
   if(ierr.ne.0) then
      write(*,*)'aforce,bforce deallocation error in remove_zeros'
      stop
   endif
! aforce bforce used here as dummy vars to remove deleted lines from amat and bmat
   allocate (aforce(dim_al,dim_ac),bforce(dim_al))
   aforce=amat(1:dim_al,1:dim_ac)
   bforce=bmat(1:dim_al)

! reallocate amat and bmat and rewrite in them
   if (allocated(amat)) deallocate (amat)
   if (allocated(bmat)) deallocate (bmat)
   allocate(amat(dim_al,dim_ac),bmat(dim_al))
   amat = aforce ; bmat = bforce
   deallocate (aforce,bforce)

!  call write_out(123,' amatrix ',amat)
!  call write_out(124,' bmatrix ',bmat)

3 format(a,9(i4))

 end subroutine include_constraints_remove_zeros
!===============================================================================
 subroutine homogeneous_constraints_overlap(n_constr)
!! now that we have the 3 matrices atransl,arot,ahuang,
!! we put them together to do the SVD
!! output s ahom and its dimensions:n_constr,dim_ac, which go to SVD
 use svd_stuff
 use params
 use ios, only : ulog
 implicit none
 integer , intent(out) :: n_constr

!-----------------------------------------------------------------------
! in this case, invariances are obtained in a rms sense and are not exact.
   dim_hom= 0
   dim_ac = nindepfc
! we put the invariance constraints into A
   if (itrans .ne. 0) dim_hom= dim_hom+ transl_constraints
   if (irot   .ne. 0) dim_hom= dim_hom+    rot_constraints
   if (ihuang .ne. 0) dim_hom= dim_hom+  huang_constraints
   allocate(ahom(dim_hom,dim_ac))
   ahom = 0d0
   write(ulog,*)' size(ahom) nlines,ncolumn=(dim(a1d))=',dim_hom,dim_ac
   n_constr = 1
   if (itrans .ne. 0) then
      ahom(n_constr:n_constr-1+transl_constraints,1:nindepfc) =  &
&                    atransl(1:transl_constraints,1:nindepfc)
      n_constr = n_constr + transl_constraints
   endif
   if (irot .ne. 0) then
      ahom(n_constr:n_constr-1+   rot_constraints,1:nindepfc) =  &
&                          arot(1:rot_constraints,1:nindepfc)
      n_constr = n_constr + rot_constraints
   endif
   if (ihuang .ne. 0) then
      ahom(n_constr:n_constr-1+ huang_constraints,1:nindepfc) =  &
&                      ahuang(1:huang_constraints,1:nindepfc)
      n_constr = n_constr + huang_constraints
   endif
   n_constr=n_constr-1
   inv_constraints=n_constr
   write(ulog,*)' ahom is setup with ',inv_constraints,' invariance constraints lines'
   if (dim_hom.ne.inv_constraints) then
      write(ulog,*)' dim_hom, inv_constraints are not equal! ', dim_hom,inv_constraints
      stop
   endif

! now form the overlap matrix from the constraints
!   allocate(overl(n_constr,n_constr))
!   do i=1,n_constr
!   do j=1,n_constr
!      overl(i,j)=dot_product(ahom(i,:),ahom(j,:))
!   enddo
!   enddo

 end subroutine homogeneous_constraints_overlap
! ============================================
 function voigt(i,j)
 integer i,j,voigt
 if (i.eq.j) then
    voigt=i
 else
    voigt=9-i-j
 endif
 end function voigt
! ============================================
 subroutine estimate_inv_constraints
!! there are 3+9n0+27n0*N3+81n0*N4^2 translational and also rotational invce constraints+15 Huang
 use atoms_force_constants
 use ios
 use geometry
 use params
 use svd_stuff
 implicit none
 integer i,m,rnk

 transl_constraints=3; rot_constraints=3; huang_constraints=0  ! 3 is for safety and for rank=1
 do rnk=2,4
! count the neighbors within a shell to get an estimate of invariance constraints
   m=0
   do i=0,nshells(rnk,1)
      m = m + atom0(1)%shells(i)%no_of_neighbors
   enddo
   write(ulog,*)' ESTIMATE_CONSTR:atm 1 has ',m,'inclusve nghbrs within ',nshells(rnk,1),' th shell'
!  if(itrans.eq.1)then
     if(include_fc(rnk).eq.1) transl_constraints= transl_constraints + (3**rnk)*natom_prim_cell*(m**(rnk-2))
!  endif
!  if(irot.eq.1)then
     if(include_fc(rnk).eq.1) rot_constraints = rot_constraints +  3**rnk*natom_prim_cell*m**(rnk-2)
!  endif
 enddo
!if(ihuang.eq.1)then
    huang_constraints = 15
!endif
 write(ulog,*)' Cumulative transl_cnstrnts estimated to be ', transl_constraints
 write(ulog,*)' Cumulative    rot_cnstrnts estimated to be ', rot_constraints

 inv_constraints=transl_constraints+rot_constraints+huang_constraints
 write(ulog,*)' Cumulative invce_cnstrnts  estimated to be ', inv_constraints


 end subroutine estimate_inv_constraints
! ============================================
 subroutine get_dim(mapd,ndimindep,ndimfull)
!! calculate the number of in/dependent fcs for all groups
 use svd_stuff
 implicit none
 integer g,ndimindep,ndimfull
 type(fulldmatrix) mapd

 ndimindep=0 ; ndimfull=0
 do g=1,mapd%ngr
    ndimindep = ndimindep + mapd%ntind(g)
    ndimfull  = ndimfull  + mapd%nt(g)
 enddo
 end subroutine get_dim
! ============================================
 subroutine check_zero_column(line,col,afrc)
! check if there is any column that is zero (the corresponding FC can not be extracted)
 use ios
 use params
 use svd_stuff
 use constants, only : r15
 implicit none
 integer line,col,i,j,noforc,nosym,rnk,dum
 real(r15) afrc(line,col)

 write(ulog,*)'==== CHECK_ZERO_COLUMN: lines and columns are:',line,col
 main: do i=1,col
! separate the contribution of symmetries from the force-displacements
! nosym counts the number of non-zero symmetry constraints
    nosym = 0
    col_loop: do j=1,line-force_constraints  ! homogeneouds part
!      write(*,*)'j=',j,line-force_constraints
       if (abs(afrc(j,i)) .gt. 1d-8) then
          nosym = 1+nosym
          exit col_loop
       endif
    enddo col_loop
    if (nosym .eq. 0) then  ! column i was all = zero in the symmetry constraints
       write(6,*)' The FC #',i,' is not involved in the symmetry constraints'
       write(6,*)' Thus it might only be fixed from the force-displacements'
       write(ulog,*)' The FC #',i,' is not involved in the symmetry constraints'
    endif

! now the contribution of force-displacements
! noforc counts the number of non-zero force-displacement lines
    noforc = 0
    col_loop2: do j=1+line-force_constraints,line
       if (abs(afrc(j,i)) .gt. 1d-8) then
          noforc = 1+noforc
          exit col_loop2
       endif
    enddo col_loop2
    if (noforc .eq. 0) then  ! column i was all = zero in the symmetry constraints
       write(ulog,*)' The FC #',i,' is not involved in the force-displacements data'
       write(6,*)' The FC #',i,' is not involved in the force-displacements data'
       write(6,*)' Thus it might only be fixed from the enforced symmetry constraints'
    endif

    if ((noforc .eq. 0) .and. (nosym .eq. 0)) then  ! column i was all = zero
       call warn(6)
       write(6,*)' The column no. ',i,' of Amatrix is 0, the corresponding '// &
&                ' FC will not be evaluated correctly because it does not appear in' // &
&                ' any of the constraints or force-displacement data '// &
&                ' perhaps a larger supercell or one of different symmetry is needed'
! for a given col number,i ,get corresponding rank and  ntindp number=dum
       call find_t(i,rnk,dum)
       if(include_fc(rnk).ne.1) then
           write(ulog,*)'ARE YOU SURE THIS RANK WAS INCLUDED?! col,rnk,gr=',i,rnk,dum
       endif
       map(rnk)%err(dum) = '*'
       write(ulog,*)' WARNING: column no. ',i,' of Amatrix is 0, the corresponding '// &
&                ' FC will not be evaluated correctly because it does not appear in' // &
&                ' any of the constraints or force-displacement data '// &
&                ' perhaps a larger supercell or one of different symmetry is needed', &
&                ' rank and group no. are ',rnk,dum

    endif
 enddo main
 write(ulog,*)' ============  END OF CHECK_ZERO_COLUMN =============='

 end subroutine check_zero_column
!=============================================================================
 subroutine find_t(t,rnk,nti)
!! for a given column index t in {1,...,ngr} find its corresponding rank and ntindp number nti
 use svd_stuff
 use ios  ! for the ulog file
 implicit none
 integer, intent(in):: t
 integer, intent(out):: rnk,nti
 integer i,m(4),gr

 m=0
 do i=1,4
    do gr=1,map(i)%ngr
       m(i)=m(i)+map(i)%ntind(gr)
    enddo
 enddo
! write(ulog,*)'FIND_T: ntindep(rnk)=',m

 if (t.ge.1 .and. t.le.m(1) ) then   ! because by default nterms(1)=1
    rnk=1
    nti = t
 elseif (t.gt.m(1) .and. t.le.m(1)+m(2) ) then
    rnk=2
    nti = t-m(1)
 elseif (t.gt.m(1)+m(2) .and. t.le.m(1)+m(2)+m(3) ) then
    rnk=3
    nti = t-m(1)-m(2)
 elseif (t.gt.m(1)+m(2)+m(3) .and. t.le. m(1)+m(2)+m(3)+m(4) ) then
    rnk=4
    nti = t-m(1)-m(2)-m(3)
 else
    print*,'FIND_T: error: t=',t,' must be less than ',  &
    &    m(1)+m(2)+m(3)+m(4)
    write(ulog,*)'FIND_T: error: t=',t,' must be less than ',  &
    &    m(1)+m(2)+m(3)+m(4)
    stop
 endif


! if (t.ge.1 .and. t.le.map(1)%ngr ) then   ! because by default nterms(1)=1
!    rnk=1
!    gr = t
! elseif (t.gt.map(1)%ngr .and. t.le.map(1)%ngr+map(2)%ngr ) then
!    rnk=2
!    gr = t-map(1)%ngr
! elseif (t.gt.map(1)%ngr+map(2)%ngr .and. t.le.map(1)%ngr+map(2)%ngr+map(3)%ngr ) then
!    rnk=3
!    gr = t-map(1)%ngr-map(2)%ngr
! elseif (t.gt.map(1)%ngr+map(2)%ngr+map(3)%ngr .and. t.le.  &
! &            map(1)%ngr+map(2)%ngr+map(3)%ngr+map(4)%ngr ) then
!    rnk=4
!    gr = t-map(1)%ngr-map(2)%ngr-map(3)%ngr
! else
!    print*,'FIND_T: error: t=',t,' must be less than ',  &
!    &    map(1)%ngr+map(2)%ngr+map(3)%ngr+map(4)%ngr
!    write(ulog,*)'FIND_T: error: t=',t,' must be less than ',  &
!    &    map(1)%ngr+map(2)%ngr+map(3)%ngr+map(4)%ngr
!    stop
! endif

 end subroutine find_t
!============================================================
 subroutine compare2previous_lines(m,n,a1d,a,counter,new)
!! compares line a1d to all lines of a matrix up to counter, and confirms whether a1d is a new line or not
 use geometry
 use params
 use ios
 use constants, only : r15
 implicit none
 integer j,t,counter,m,n,j0 !,i
 logical new
 real(r15) junk,zero
 real(r15), dimension(n) :: a1d,b
 real(r15), dimension(m,n) :: a

 zero = 0d0
 new=.true.
 j0=0
    do j=1,n  ! = no of columns in ared; it's really ngr
       if ( a1d(j) .myeq. zero ) cycle
       j0 = j   ! it is the first nonzero index (column)
       exit
    enddo
    if(j0.eq.0) then
       new=.false.  ! all a1d was zero
       return
    endif

    checkloop: do t=1,counter
! check to see if they are multiples of each other
!   write(ulog,*)' comparing with line ',t,' in ared matrix'
  !       if (j.gt.n .and. verbose ) then
  !          write(ulog,*)' this line of ared1d is all=0'
  !          write(ulog,4)(a1d(i),i=1,n)
! !          stop
  !       endif
          junk = a(t,j0)/a1d(j0)
          b = a1d*junk   ! this is to avoid repeating amat and -amat
          if ( a(t,:) .myeq. b(:) ) then
             new=.false.
             return !exit checkloop
          endif
    enddo checkloop
    if(verbose .and. new) then
         write(ulog,5)'NEW LINE:',counter+1,a1d
    endif

4 format(99(1x,g11.4))
5 format(a,i5,99(1x,f7.3))
 end subroutine compare2previous_lines

!===============================================================================
 subroutine check_huang(phi2,ngr)
! this subroutine checks Huang invariance relations:
! sum_r ( phi_0,r^al,be r^ga r^de -  phi_0,r^ga,de r^al r^be ) = 0
 use params
 use atoms_force_constants
 use svd_stuff
 use ios
 use constants, only : r15
! use geometry
! use force_constants_module
 implicit none
 integer, intent(in) :: ngr
 real(r15), intent(in) :: phi2(ngr)
 integer i,j,al,be,ga,de,t,ti,cnt2,mterm,rnk,voigt,g,vab,vcd,res
 real(r15) huang1,huang2,rr(3),stress(6)

 write(ulog,*)'CHECKING HUANG INVARIANCE RELATIONS       ###############################'
 res = sum(map(1)%ntind(:))
 write(ulog,*)' # of indep terms of rank 1 = res=',res
 rnk = 2   !-----------------------------------------------------
 write(ulog,*)' i0 , al , be , ga , de ,mterm, I_huang '

 do al=1,3
 do be=al,3
  do ga=1,3
  do de=ga,3
    vab = voigt(al,be)
    vcd = voigt(ga,de)
    if (vab .le. vcd) cycle
    huang1 = 0
    huang2 = 0
    mterm = 0
    do i=1,natom_prim_cell
      cnt2=0
      do g=1,map(2)%ngr
         if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
      do t=1,map(2)%nt(g)
         if ( map(2)%gr(g)%iat(1,t).ne.i .or. map(2)%gr(g)%ixyz(1,t).ne.al  &
       & .or. map(2)%gr(g)%ixyz(2,t).ne.be ) cycle
         mterm = mterm+1
         j  =  map(2)%gr(g)%iat(2,t)
         rr = atompos(:,j)-atompos(:,i)
         do ti =1,map(2)%ntind(g)
            huang1 = huang1 + rr(ga)*rr(de)*map(2)%gr(g)%mat(t,ti)*phi2(res+cnt2+ti)
         enddo
      enddo
      enddo
      cnt2=0
      do g=1,map(2)%ngr
         if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
      do t=1,map(2)%nt(g)
         if ( map(2)%gr(g)%iat(1,t).ne.i .or. map(2)%gr(g)%ixyz(1,t).ne.ga   &
    &    .or. map(2)%gr(g)%ixyz(2,t).ne.de ) cycle
         mterm = mterm+1
         j  =  map(2)%gr(g)%iat(2,t)
         rr = atompos(:,j)-atompos(:,i)
         do ti =1,map(2)%ntind(g)
            huang2 = huang2 - rr(al)*rr(be)*map(2)%gr(g)%mat(t,ti)*phi2(res+cnt2+ti)
         enddo
      enddo
      enddo
    enddo
    write(ulog,5)i,al,be,ga,de,vab,vcd,mterm,huang1+huang2

    stress(vab)=stress(vab)+huang1+huang2

  enddo
  enddo
 enddo
 enddo
 write(ulog,*)'HUANG INVARIANCE RELATIONS CHECKED        ###############################'
 write(ulog,4)'STRESS Tensor=',stress

4 format(a,9(1x,g11.4))
5 format(8(1x,i4),9(1x,g13.6))

 end subroutine check_huang
! ============================================
 subroutine truncate(n,ws,n_kept,mw)
! find a gap in the ws and drop the small ws values
 use constants, only : r15
 implicit none
 integer n,n_kept,i,i1
 real(r15) ws(n), logws(n),gap(n-1),maxgap1,wmax
 integer mw(n)

 wmax=maxval(ws)

! sort the ws in descending order (largest first)
 call invsort(n,ws,mw)
 ws=ws/wmax
 logws=log(ws)/log(10.)  ! decimal log of ws
 if (logws(1) .ne.0) then
    write(*,*)'TRUNCATE: ws(1),logws(1)=',ws(1),logws(1)
    stop
 endif
 write(*,*)'TRUNCATE: dimension of ws is n=',n

! calculate the level separations
 do i=1,n-1
!   gap(i)=logws(i)-logws(i+1)
    gap(i)=logws(mw(i))-logws(mw(i+1))
    write(*,5)'TRUNCATE: i, gap(i) logws, ws=',i,gap(i),logws(mw(i)),ws(mw(i))
 enddo
5 format(a,i4,9(1x,g11.4))

 do i=1,n
    if (logws(mw(i)).gt. -5) cycle
    exit
 enddo
 n_kept=i  ! any logws above -5, i.e. ws>10^-5 is kept

 maxgap1=maxval(gap(i+1:n-1))
 i1=maxloc(gap(i+1:n-1),n-i-1)
 n_kept=i1

 ws=ws*wmax

! if(i1.lt.n_kept) then
!! the largest gap is before n_kept , ie larger than 10^-6
!! use the second largest gap
!   i2=maxloc(gap(i1+1:n-1))
! else
!! the largest gap comes after, and so that's where we will truncate
!   n_kept=i1
!   return
! endif


 end subroutine truncate
!===========================================================
 subroutine setup_FC2_in_supercell(rws26)
!! selects only the independent FC2s which are in the WS cell of the supercell, including its boundary, 
!! defined by rws26(3,nrgrid). If map(rnk)%keep(i)=1, then it is in the WS and is kept
!! for fc2range=0 everything in WS and all its images by symmetry are kept
!! for fc2range\=0 all the requested shells are kept provided one of them is within the WS of the supercell
! use lattice
 use ios
 use atoms_force_constants
 use geometry
 use svd_stuff
 use params
 use constants, only : r15
 implicit none
 integer g,ti,t,i0,j,cnt,nboundary
 real(r15), intent(in) :: rws26(3,26)
 real(r15) rij(3),dij(999),transl(3),rij0(3)
 logical inside

 write(ulog,*)'SETUP_FC2: dimension of map(2)%keep(:)=',size(map(2)%keep(:))
 map(2)%keep(:)=0
 map(2)%nkeptind=0  ! this is the size in amatrix of the FC2s not eliminated<map(2)%ntind.
 do g=1,map(2)%ngr
    do ti=1,map(2)%ntind(g)
    cnt= counteri(2,g,ti)
    tloop: do t =1,map(2)%nt(g) ! loop over t until one of them is in WS; id the (g,ti) pair
       i0=map(2)%gr(g)%iat(1,t)
       j =map(2)%gr(g)%iat(2,t)  ! this j is for atompos
       rij=atompos(:,j)-atompos(:,i0) ! is this within the WS cell?
!      write(ulog,*)'g,t,cnt,rij=',g,ti,cnt,rij
       call check_inside_ws(rij,rws26,inside,nboundary)
! Do we have elements of the same group, i.e. same distance, one being inside and one outside?
       if (inside .and. map(2)%gr(g)%mat(t,ti) .ne.0) then
         map(2)%keep(cnt)=1 
         write(ulog,5)'Group,ti ',g,ti,' term of distance rij=',length(rij),' is kept'     
         exit tloop
       endif
    enddo tloop
    if(map(2)%keep(cnt).eq.0) write(ulog,5)'Group,ti ',g,ti,' term of distance rij=',length(rij),' NOT kept'     
    enddo
 enddo
 map(2)%nkeptind=sum(map(2)%keep(:))  ! number of independent FC2 terms kept

 write(ulog,*)'size of groups of FC2s:ngr, size(kept indep terms)=',map(2)%ngr,size(map(2)%keep)
 write(ulog,*)'number of kept FC2 independent terms=',map(2)%nkeptind, 'out of ',map(2)%ntotind
 write(ulog,*)'Total # of independent terms kept   =',sum(map(2)%keep)
 write(ulog,*)'# i, map(2)%keep(i)'

5 format(a,2i4,a,f9.4,a,3i5)

 end subroutine setup_FC2_in_supercell
!==============================================
 subroutine exclude_beyond_sc(n,keep)
!! selects only the independent FC2s which are in the WS cell of the supercell, including its boundary, 
!! defined by rws26(3,nrgrid).
 use lattice
 use ios
 use atoms_force_constants
 use geometry
 use svd_stuff
 use params
 use constants, only : r15
 implicit none
 integer, intent(inout) :: n,keep(n)
 integer g,ti,t,i0,j,cnt,nboundary
 real(r15) rij(3),dij(999),transl(3),rij0(3)
 logical inside

 write(ulog,*)'EXCLUDE_BEYOND_SC: dimension of keep=',n
 keep=0
 do g=1,map(2)%ngr
    do ti=1,map(2)%ntind(g)
    tloop: do t =1,map(2)%nt(g)
       i0=map(2)%gr(g)%iat(1,t)
       j =map(2)%gr(g)%iat(2,t)  ! this j is for atompos
       rij=atompos(:,j)-atompos(:,i0) ! is this within the WS cell?
       cnt= counteri(2,g,ti)
!      write(ulog,4)'g,t,cnt,rij=',g,ti,cnt,rij
       call check_inside_ws(rij,rws26,inside,nboundary)
!      write(ulog,*)'inside is ',inside
! Do we have elements of the same group, i.e. same distance, one being inside and one outside?
       if(inside) then  ! keep if at least one of the terms is inside
          keep(cnt)=1 
          write(ulog,5)'Group,ti ',g,ti,' term of distance rij=',length(rij),' is kept'     
          exit tloop
       endif
    enddo tloop
    enddo
 enddo

 write(ulog,*)'Total # of independent terms kept, map2%ntotind =',sum(keep),map(2)%ntotind

4 format(a,3i4,9f9.4)
5 format(a,2i4,a,f9.4,a,3i5)

 end subroutine exclude_beyond_sc
!==============================================
 subroutine implement_temperature(m3,n,amat,bmat,nconfg,ene,nlines,tempk,mat,qmat)
! input is the inhomogeneous part of amat; output is a square matrix mat(n,n) and qmat(n)
 use constants
 implicit none
 integer, intent(in) :: m3,n,nconfg,nlines(nconfg) ! # of lines for each configuration
 real(r15), intent(in) :: amat(m3,n),bmat(m3), ene(nconfg),tempk
 real(r15), intent(out):: mat(n,n),qmat(n)
 real(r15) wei,kt
 integer s,nat3,nat,nl,res,i1,in

 kt=tempk*k_b/ee  ! kT in eV
 if( nint((m3*1d0)/(3*nconfg)).ne.(m3*1d0)/(3*nconfg) ) then
    write(*,*)'m3 not divisible by 3*number of snaps ',m3,nconfg
 else
    nat3=(m3)/nconfg
    nat=nat3/3
    write(*,*)'TEMPERATURE: We have ',nat,' atoms!'
 endif
 nl=(m3)/nconfg
 if(sum(nlines).ne.m3) then
   write(*,*)' IMPLEMENT T: sum(nlines) \= force_constraints:',sum(nlines),m3  
   stop
 endif

 res=0
 mat=0
 qmat=0
 do s=1,nconfg
    wei=exp(-ene(s)/kt)
    if(s.ne.1) res=res+nlines(s-1)
    i1=res+1
    in=res+nlines(s)
    mat = mat+wei*matmul(transpose(amat(i1:in,1:n)),amat(i1:in,1:n))
    qmat=qmat+wei*matmul(transpose(amat(i1:in,1:n)),bmat(i1:in    ))
 enddo

 end subroutine implement_temperature
