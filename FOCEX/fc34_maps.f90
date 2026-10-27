!******************************************************************************
! setup_maps_fc34, size_force_constants_fc34, collect_force_constants_fc34:
! copies of setup_maps (others.f90), size_force_constants and
! collect_force_constants (force_constant_symmetries.f90) with the fcflag=0 range
! rule of module fc34_range added (user, 2026-10-26: keep the original routines
! unchanged, give the new version its own name). focex.f90 calls setup_maps_fc34
! iff some fcflag(3:8) is 0; otherwise the original setup_maps runs.
!******************************************************************************
!===============================================================================
 subroutine setup_maps_fc34
!! maps the output of collect_fcs to the new arrays for each rank,
!! with the proper (allocatable) dimensions
 use svd_stuff
 use params
 use atoms_force_constants
 use ios
 use geometry
 use fc34_range
 implicit none
 character lineout*80,frmt*90
 integer i0,j,i,n,m,mxzero,mx,mxi,l,rnk,ichkloop,mxnt,mxnti   &
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
  if (rnk.eq.3) call fc34_init   ! kept FC2 pairs for the fcflag=0 rule (module fc34_range)
  if ( include_fc(rnk) .ne. 0 ) then

    ierz=0; iert=0;  ieri=0; ierg=0

    write(ulog,4)' ******************** FOR rank *****************= ',rnk
    write(*   ,4)' ***** FOR rank = ',rnk

 4 format(a,9(i5))

! ---------------------------------------------------------------------------
! First pass: find the exact sizes. This allocates nterm(1:ngroups) and
! ntermsindep(1:ngroups) and returns the largest of each, so that every array
! below is allocated once, at the size actually needed. No guessing, and no
! restart of the symmetry analysis.
! ---------------------------------------------------------------------------
     call size_force_constants_fc34(rnk,nshells(rnk,:),ngroups(rnk),  &
     &       ntermszero,mxnt,mxnti,mx,mxi,mxzero)

     write(ulog,4)'EXACT SIZES: rank,ngroups,max nterm,max ntindep,ntermszero=', &
     &             rnk,ngroups(rnk),mxnt,mxnti,ntermszero
     write(*   ,4)'EXACT SIZES: rank,ngroups,max nterm,max ntindep,ntermszero=', &
     &             rnk,ngroups(rnk),mxnt,mxnti,ntermszero

! One slot of headroom each: collect_force_constants counts a candidate group
! before it finds that all its terms vanish (and removes it again), so ngroups
! reaches the final number + 1 if that happens after the last group (MgO, rank 3);
! it rejects on ngroups.gt.maxgroups and on ntermszero.ge.maxtermszero. mx and
! mxi are the work bounds that size_force_constants found to be sufficient.
     mxgrps = max(ngroups(rnk),1) + 1
     mxzero = max(ntermszero ,1) + 1

! ---------------------------------------------------------------------------
! Second pass: allocate once, at the size actually needed, then fill.
! ---------------------------------------------------------------------------
     allocate( mapmat(mx,mxi,mxgrps)); mapmat=0
     allocate( iatmtrm(rnk,mx,mxgrps));iatmtrm=0
     allocate( ixyztrm(rnk,mx,mxgrps));ixyztrm=0
     allocate( iatmtermindp(rnk,mx,mxgrps));iatmtermindp=0
     allocate( ixyztermindp(rnk,mx,mxgrps));ixyztermindp=0
     allocate( iatmtermzero(rnk,mxzero));iatmtermzero=0
     allocate( ixyztermzero(rnk,mxzero));ixyztermzero=0
     if(allocated(nterm))       deallocate(nterm)
     if(allocated(ntermsindep)) deallocate(ntermsindep)
     allocate(nterm(mxgrps),ntermsindep(mxgrps)); nterm=0; ntermsindep=0

     call collect_force_constants_fc34(rnk,nshells(rnk,:),  &
     &       rnk,mx,mxi,mxzero,mxgrps,   &
     &       ngroups(rnk),mapmat,  &
     &       ntermsindep,iatmtermindp,ixyztermindp,  &
     &       nterm,iatmtrm,ixyztrm,   &
     &       ntermszero,iatmtermzero,ixyztermzero,   &
     &       ierz,iert,ieri,ierg)

     if (ierz+iert+ieri+ierg .ne. 0) then
        write(ulog,4)'SETUP_MAPS: sizes from size_force_constants were rejected'// &
     &               ' by collect_force_constants, ier z,t,i,g=',ierz,iert,ieri,ierg
        write(*   ,4)'SETUP_MAPS: sizing disagreement, ier z,t,i,g=',ierz,iert,ieri,ierg
        stop
     endif

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
 call write_clusters   ! clusters<rank>.dat, for display

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
 end subroutine setup_maps_fc34
!****************************************************************************
      subroutine size_force_constants_fc34(nrank,nshell,ngroups, &
     &     ntermszero,maxnterm,maxntermsindep,mxwork,mxiwork,mxzwork)
!! Determines the EXACT sizes of the force-constant arrays of rank nrank, before
!! any of them is allocated. It walks the same group-discovery loop as
!! collect_force_constants, but keeps only what is needed to recognise a term
!! that has already been seen, and throws away the coefficient matrix. Its
!! outputs are meant to be used to allocate the arrays that
!! collect_force_constants then fills, so that no array is ever over-allocated
!! and the symmetry analysis never has to be restarted with a larger guess.
!!
!! arguments:
!!     nrank (input), order of the derivative
!!     nshell (input), neighbor shells to include, per atom of the primitive cell
!!     ngroups (output), number of groups of symmetry-related terms
!!     ntermszero (output), total number of terms that vanish by symmetry
!!     maxnterm, maxntermsindep (output), the largest per-group count of all
!!          terms and of independent terms
!!     mxwork, mxiwork, mxzwork (output), work-array bounds that were sufficient
!!          for every call to force_constants. These are what the caller must
!!          allocate, NOT maxnterm: inside force_constants the same bound
!!          maxterms also limits the number of symmetry equations neqs, which is
!!          generally larger than the number of terms in the group.
!!
!! The per-group counts themselves are returned in the module arrays
!! svd_stuff::nterm and svd_stuff::ntermsindep, which this routine allocates to
!! exactly ngroups elements.
!!
!! The per-group work arrays handed to force_constants are grown and the call
!! retried when they turn out to be too small. That retry costs one group, not
!! the whole analysis.
      use atoms_force_constants
      use fc34_range
      use svd_stuff, only : nterm,ntermsindep
      use ios , only : ulog
      use constants, only : r15
      implicit none
      integer, intent(in) :: nrank,nshell(natom_prim_cell)
      integer, intent(out):: ngroups,ntermszero,maxnterm,maxntermsindep
      integer, intent(out):: mxwork,mxiwork,mxzwork

! growable store of the terms of every group found so far, kept flat:
! group g occupies keyat(:,gofs(g)+1 : gofs(g)+gcnt(g))
      integer, allocatable :: keyat(:,:),keyxyz(:,:),gofs(:),gcnt(:)
      integer, allocatable :: zat(:,:),zxyz(:,:)
      integer, allocatable :: ntg(:),ntig(:)
! per-group scratch handed to force_constants
      integer, allocatable :: iatind(:,:),ixyzind(:,:),iatall(:,:),ixyzall(:,:)
      integer, allocatable :: iat0(:,:),ixyz0(:,:)
      real(r15), allocatable :: amloc(:,:)

      integer i,j,k,m,n0,iatom,iatom0,iatomd(nrank),ixyz,ngroupsave,nzsave
      integer ixyzd(nrank),icell(3),iatomd2(nrank),ixyzd2(nrank)
      integer nkey,mxkey,mxgrp,mxzero,mxt,mxti,mxt0,ierz,iert,ieri,itry,g,i1,i2
      integer ntloc,ntiloc
      logical firsttime,seen,fc34_on

      if(nrank.gt.8)then
        write(ulog,*)'SIZE_FORCE_CONSTANTS: nrank too large ',nrank
        stop
      endif

! initial capacities; these only affect how often we reallocate, never the result
      mxkey = 4096 ; mxgrp = 256 ; mxzero = 1024
      allocate(keyat(nrank,mxkey),keyxyz(nrank,mxkey))
      allocate(gofs(mxgrp),gcnt(mxgrp),ntg(mxgrp),ntig(mxgrp))
      allocate(zat(nrank,mxzero),zxyz(nrank,mxzero))
! per-group scratch: start modest, grown on demand below
      mxt = 64 ; mxti = 16 ; mxt0 = 64
      call alloc_scratch

      ngroups=0 ; ntermszero=0 ; nkey=0
      fc34_on = fc34_active(nrank)

      iatom0loop: do iatom0=1,natom_prim_cell
        do iatom=1,natoms
          if(iatomneighbor(iatom0,iatom).le.nshell(iatom0))exit
          if(iatom.eq.natoms)then
            write(6,*)'SIZE_FORCE_CONSTANTS: first nearest neighbor not found'
            stop
          endif
        enddo
        iatomd(1:nrank)=iatom
        iatomd(nrank)=iatom-1
        iatomd(1)=iatom0
        firsttime=.true.
        nextatomloop: do while(.true.)
          if(firsttime)then
            firsttime=.false.
          else
            if(nrank.eq.1)exit nextatomloop
          endif
          iloop: do i=nrank,2,-1
            jloop: do j=iatomd(i)+1,natoms
              if(fc34_on) then        ! fcflag=0 range rule, module fc34_range
                iatomd(i)=j
                if(.not.fc34_partial(i,iatomd(1:i),nshell(iatom0))) cycle jloop
                iatomd(i:nrank)=j
                exit iloop
              endif
              if(iatomneighbor(iatom0,j).le.nshell(iatom0))then
                do k=2,i-1
                  do m=1,3
                    icell(m)=iatomcell(m,iatomd(k))-iatomcell(m,j)
                  enddo
                  call find_atompos(icell,iatomcell0(iatomd(k)),m)
                  if(m.eq.0)cycle jloop
                  if(iatomneighbor(iatomcell0(j),m).gt.nshell(iatom0)) cycle jloop
                enddo
                iatomd(i:nrank)=j
                exit iloop
              endif
            enddo jloop
            if(i.eq.2)cycle iatom0loop
          enddo iloop

! fcflag=0 range rule on the complete cluster
          if(fc34_on) then
            if(.not.fc34_final(nrank,iatomd,nshell(iatom0))) cycle nextatomloop
          endif
! already known, as an atom set?
          ixyzd=1
          call unique_force_constant(nrank,iatomd,ixyzd,iatomd2,ixyzd2)
          do i=1,nkey
            seen=.true.
            do j=1,nrank
              if(keyat(j,i).ne.iatomd2(j))then
                 seen=.false. ; exit
              endif
            enddo
            if(seen) cycle nextatomloop
          enddo
          do i=1,ntermszero
            seen=.true.
            do j=1,nrank
              if(zat(j,i).ne.iatomd2(j))then
                 seen=.false. ; exit
              endif
            enddo
            if(seen) cycle nextatomloop
          enddo

          ngroupsave=ngroups
          nzsave=ntermszero
          ixyzloop: do ixyz=1,3**nrank
            m=ixyz-1
            do i=1,nrank
              ixyzd(i)=mod(m,3)+1
              m=m/3
            enddo
            call unique_force_constant(nrank,iatomd,ixyzd, iatomd2,ixyzd2)
! already generated in this pass over the coordinate sets?
            do g=ngroupsave+1,ngroups
              i1=gofs(g)+1 ; i2=gofs(g)+gcnt(g)
              do i=i1,i2
                seen=.true.
                do j=1,nrank
                  if(keyat(j,i).ne.iatomd2(j).or.keyxyz(j,i).ne.ixyzd2(j))then
                     seen=.false. ; exit
                  endif
                enddo
                if(seen) cycle ixyzloop
              enddo
            enddo
            do i=nzsave+1,ntermszero
              seen=.true.
              do j=1,nrank
                if(zat(j,i).ne.iatomd2(j).or.zxyz(j,i).ne.ixyzd2(j))then
                   seen=.false. ; exit
                endif
              enddo
              if(seen) cycle ixyzloop
            enddo

! new group: generate its terms, growing the scratch until it fits
            itry=0
            retry: do
              itry=itry+1
              if(itry.gt.40)then
                write(ulog,*)'SIZE_FORCE_CONSTANTS: scratch growth did not converge'
                stop
              endif
              call force_constants(nrank,amloc,iatomd2,ixyzd2,  &
     &             ntiloc,iatind,ixyzind, ntloc,iatall,ixyzall,  &
     &             n0,iat0,ixyz0, nrank,mxt,mxti,mxt0, ierz,iert,ieri)
              if(ierz.eq.0 .and. iert.eq.0 .and. ieri.eq.0) exit retry
              if(iert.ne.0) mxt  = 2*mxt
              if(ieri.ne.0) mxti = 2*mxti
              if(ierz.ne.0) mxt0 = 2*mxt0
              if(mxti.gt.mxt) mxti=mxt
              call alloc_scratch
            enddo retry

            ngroups=ngroups+1
            call grow_groups(ngroups)
            ntg (ngroups)=ntloc
            ntig(ngroups)=ntiloc

! record this group's terms for the duplicate tests, and the new zero terms
            call grow_keys(nkey+ntg(ngroups))
            gofs(ngroups)=nkey
            gcnt(ngroups)=ntg(ngroups)
            do i=1,ntg(ngroups)
              keyat (1:nrank,nkey+i)=iatall (1:nrank,i)
              keyxyz(1:nrank,nkey+i)=ixyzall(1:nrank,i)
            enddo
            nkey=nkey+ntg(ngroups)

            call grow_zeros(ntermszero+n0)
            do i=1,n0
              zat (1:nrank,ntermszero+i)=iat0 (1:nrank,i)
              zxyz(1:nrank,ntermszero+i)=ixyz0(1:nrank,i)
            enddo
            ntermszero=ntermszero+n0

! a group with no surviving term is not a group
            if(ntg(ngroups).eq.0)then
              nkey=nkey-gcnt(ngroups)
              ngroups=ngroups-1
            endif
          enddo ixyzloop
        enddo nextatomloop
      enddo iatom0loop

      if(allocated(nterm))       deallocate(nterm)
      if(allocated(ntermsindep)) deallocate(ntermsindep)
      allocate(nterm(max(ngroups,1)),ntermsindep(max(ngroups,1)))
      nterm=0 ; ntermsindep=0
      if(ngroups.gt.0)then
        nterm      (1:ngroups)=ntg (1:ngroups)
        ntermsindep(1:ngroups)=ntig(1:ngroups)
      endif
      maxnterm       = 0
      maxntermsindep = 0
      if(ngroups.gt.0)then
        maxnterm       = maxval(nterm(1:ngroups))
        maxntermsindep = maxval(ntermsindep(1:ngroups))
      endif

      mxwork  = mxt
      mxiwork = mxti
      mxzwork = mxt0

      write(ulog,3)'SIZE_FORCE_CONSTANTS: rank,ngroups,ntermszero=',nrank,ngroups,ntermszero
      write(ulog,3)'SIZE_FORCE_CONSTANTS: sum(nt),sum(ntind),max nt,max ntind=', &
     &      sum(nterm),sum(ntermsindep),maxnterm,maxntermsindep
      write(ulog,3)'SIZE_FORCE_CONSTANTS: work bounds mxt,mxti,mxt0=',mxwork,mxiwork,mxzwork
      write(*   ,3)'SIZE_FORCE_CONSTANTS: rank,ngroups,ntermszero=',nrank,ngroups,ntermszero

      deallocate(keyat,keyxyz,gofs,gcnt,zat,zxyz,ntg,ntig)
      deallocate(iatind,ixyzind,iatall,ixyzall,iat0,ixyz0,amloc)

 3    format(a,99(i8))

      contains

      subroutine alloc_scratch
      if(allocated(iatind )) deallocate(iatind,ixyzind,iatall,ixyzall,iat0,ixyz0,amloc)
      allocate(iatind(nrank,mxt),ixyzind(nrank,mxt))
      allocate(iatall(nrank,mxt),ixyzall(nrank,mxt))
      allocate(iat0(nrank,mxt0),ixyz0(nrank,mxt0))
      allocate(amloc(mxt,mxti))
      end subroutine alloc_scratch

      subroutine grow_keys(need)
      integer, intent(in) :: need
      integer, allocatable :: t1(:,:),t2(:,:)
      if(need.le.mxkey) return
      do while(mxkey.lt.need)
         mxkey=2*mxkey
      enddo
      allocate(t1(nrank,mxkey),t2(nrank,mxkey)); t1=0; t2=0
      t1(:,1:nkey)=keyat(:,1:nkey) ; t2(:,1:nkey)=keyxyz(:,1:nkey)
      call move_alloc(t1,keyat) ; call move_alloc(t2,keyxyz)
      end subroutine grow_keys

      subroutine grow_zeros(need)
      integer, intent(in) :: need
      integer, allocatable :: t1(:,:),t2(:,:)
      if(need.le.mxzero) return
      do while(mxzero.lt.need)
         mxzero=2*mxzero
      enddo
      allocate(t1(nrank,mxzero),t2(nrank,mxzero)); t1=0; t2=0
      t1(:,1:ntermszero)=zat(:,1:ntermszero) ; t2(:,1:ntermszero)=zxyz(:,1:ntermszero)
      call move_alloc(t1,zat) ; call move_alloc(t2,zxyz)
      end subroutine grow_zeros

      subroutine grow_groups(need)
      integer, intent(in) :: need
      integer, allocatable :: t1(:),t2(:),t3(:),t4(:)
      integer nold
      if(need.le.mxgrp) return
      nold=mxgrp
      do while(mxgrp.lt.need)
         mxgrp=2*mxgrp
      enddo
      allocate(t1(mxgrp),t2(mxgrp),t3(mxgrp),t4(mxgrp)); t1=0;t2=0;t3=0;t4=0
      t1(1:nold)=gofs(1:nold) ; t2(1:nold)=gcnt(1:nold)
      t3(1:nold)=ntg (1:nold) ; t4(1:nold)=ntig(1:nold)
      call move_alloc(t1,gofs) ; call move_alloc(t2,gcnt)
      call move_alloc(t3,ntg ) ; call move_alloc(t4,ntig)
      end subroutine grow_groups

      end subroutine size_force_constants_fc34
!****************************************************************************
       subroutine collect_force_constants_fc34(nrank,nshell,   &
     &     maxrank,maxterms,maxtermsindep,maxtermszero,maxgroups,  &
     &     ngroups,amat,  &
     &     ntermsindep,iatomtermindep,ixyztermindep,    &
     &     ntermsall,     iatomterm,     ixyzterm,    &
     &     ntermszero, iatomtermzero, ixyztermzero, &
     &     ierz,iert,ieri,ierg)

!     subroutine collect_force_constants(nrank,nshell,ngroups,  &
!     &     ntermsindep,iatomtermindep,ixyztermindep,ntermsall,  &
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
!     ntermsall(i) (output), total number of terms
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
      use fc34_range
!      use force_constants_module
      use ios , only : ulog
 use constants, only : r15
      implicit none
      integer, intent(in) :: nrank,maxrank,maxterms,maxtermsindep,maxgroups,  &
     &     maxtermszero,nshell(natom_prim_cell)
      integer, intent(out):: ierz,iert,ieri,ierg,ngroups,ntermszero,ntermsall(maxgroups), &
     &     ntermsindep(maxgroups),iatomtermindep(maxrank,maxterms,maxgroups),  &
     &     ixyztermindep(maxrank,maxterms,maxgroups), &
     &     iatomterm(maxrank,maxterms,maxgroups),ixyzterm(maxrank,maxterms,maxgroups),  &
     &     iatomtermzero(maxrank,maxtermszero),ixyztermzero(maxrank,maxtermszero)
      real(r15), intent(out):: amat(maxterms,maxtermsindep,maxgroups)
      integer i,j,k,m,n0,ncount,iatom,iatom0,iatomd(nrank),ixyz,ngroupsave,  &
     &     ixyzd(nrank),icell(3),iatomd2(nrank),ixyzd2(nrank),ntermszerosave
      logical firsttime,fc34_on

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
      fc34_on = fc34_active(nrank)
      amat=0
      ntermszero=0
! do each atom in unit cell
      iatom0loop: do iatom0=1,natom_prim_cell
! find first nearest neighbor in list
        do iatom=1,natoms
          if(iatomneighbor(iatom0,iatom).le.nshell(iatom0))exit
          if(iatom.eq.natoms)then
            write(6,*)'COLLECT_FORCE_CONSTANTS: first nearest neighbor not found'
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
              if(fc34_on) then        ! fcflag=0 range rule, module fc34_range
                iatomd(i)=j
                if(.not.fc34_partial(i,iatomd(1:i),nshell(iatom0))) cycle jloop
                iatomd(i:nrank)=j
                exit iloop
              endif
! nearest neighbor to atom in cell at origin?
              if(iatomneighbor(iatom0,j).le.nshell(iatom0))then
! nearest neighbor to all other atoms in denominator?
                do k=2,i-1
                  do m=1,3
                    icell(m)=iatomcell(m,iatomd(k))-iatomcell(m,j)
                  enddo
                  call find_atompos(icell,iatomcell0(iatomd(k)),m)
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
! fcflag=0 range rule on the complete cluster
          if(fc34_on) then
            if(.not.fc34_final(nrank,iatomd,nshell(iatom0))) cycle nextatomloop
          endif
! check if found yet
          ixyzd=1
          call unique_force_constant(nrank,iatomd,ixyzd,iatomd2,ixyzd2)
          do k=1,ngroups
          iloop2: do i=1,ntermsall(k)
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
            iloop3: do i=1,ntermsall(k)
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

            if(ngroups.gt.maxgroups)then     ! slots 1..maxgroups are valid
!! can add to the size of arrays ixyzterm*, iatomterm* 
              write(6,*)' maxgroups is too small ',maxgroups
              ierg=ierg+1
              return
            endif

            if (ntermszero.ge.maxtermszero) then
!! can add to the size of arrays iatomtermzero,ixyztermzero
              write(*,3)' mxt0 =< nt0=',maxtermszero,ntermszero
              write(ulog,3)' mxt0 =< nt0=',maxtermszero,ntermszero
              ierz=5
              return
            endif


            call force_constants(nrank,amat(1,1,ngroups),iatomd2,ixyzd2,  &
     &           ntermsindep(ngroups),iatomtermindep(1,1,ngroups),ixyztermindep(1,1,ngroups), &
     &           ntermsall(ngroups),iatomterm(1,1,ngroups),ixyzterm(1,1,ngroups),  &
     &           n0,iatomtermzero(1,ntermszero+1),ixyztermzero(1,ntermszero+1),  &
     &           maxrank,maxterms,maxtermsindep,maxtermszero-ntermszero,  &
     &           ierz,iert,ieri) !,ierg)

            write(*,3)'OUTPUT OF FORCE_CONSTANTS FOR atom,xyz RANK/ngroups=',iatom0,ixyz,nrank,ngroups
!           write(*,3)'OUTPUT OF FORCE_CONSTANTS: ntermszero =',ntermszero
!            write(ulog,3)'OUTPUT OF FORCE_CONSTANTS FOR RANK/ngroups=',nrank,ngroups
!            write(ulog,3)'OUTPUT OF FORCE_CONSTANTS: ntermsindep=',ntermsindep(1:ngroups)
!            write(ulog,3)'OUTPUT OF FORCE_CONSTANTS: nterms     =',ntermsall(1:ngroups)
!            write(ulog,3)'OUTPUT OF FORCE_CONSTANTS: ntermszero =',ntermszero

            if (ierz.ne.0) then
              write(*,3)' maxtermszero is too small,nt0=',maxtermszero,ierz,ntermszero
              return
            endif
            if (iert.ne.0) then
              write(*,3)' maxterms is too small,nterms',maxterms,iert,sum(ntermsall)
              return
            endif
            if (ieri.ne.0) then
              write(*,3)' maxtermsindep is too small,ntindep',maxtermsindep,ieri,sum(ntermsindep)
              return
            endif
            ntermszero=ntermszero+n0
! no nonzero terms in group
            if(ntermsall(ngroups).eq.0)then
              ngroups=ngroups-1
              ncount=ncount+1
            endif
! next set of coordinates
          enddo ixyzloop
! next set of atoms
        enddo nextatomloop
        write(*,3)'OUTPUT OF FORCE_CONSTANTS FOR atom,xyz RANK/ngroups=',iatom0,ixyz,nrank,ngroups
        write(*,3)'OUTPUT OF FORCE_CONSTANTS: ntermsindep=',ntermsindep(1:ngroups)
        write(*,3)'OUTPUT OF FORCE_CONSTANTS: nterms     =',ntermsall(1:ngroups)
        write(*,3)'OUTPUT OF FORCE_CONSTANTS: ntermszero =',ntermszero
      enddo iatom0loop

!      write(ulog,*)' ngroups=',ngroups
!      write(ulog,3)' ntermsindep        (1:ngroups)=', ntermsindep(1:ngroups)
!      write(ulog,3)' ntermsall             (1:ngroups)=', ntermsall(1:ngroups)
!      write(ulog,3)' iatomtrmindep(rnk,1,1:ngroups)=', iatomtermindep(nrank,1,1:ngroups)
!      write(ulog,3)' ixyztermindep(rnk,1,1:ngroups)=', ixyztermindep(nrank,1,1:ngroups)
!      write(ulog,3)' iatomterm    (rnk,1,1:ngroups)=', iatomterm(nrank,1,1:ngroups)
!      write(ulog,3)' ixyzterm     (rnk,1,1:ngroups)=', ixyzterm(nrank,1,1:ngroups)
!      write(ulog,3)' ntermszero =', ntermszero

      ierz=0 ; iert=0 ; ieri=0 ; ierg=0

 3    format(a,99(i6))
      end subroutine collect_force_constants_fc34
