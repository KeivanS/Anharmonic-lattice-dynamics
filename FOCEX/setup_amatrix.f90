 ! In the following 3 subroutines, only the treatment of FC2s has changed: instead
 ! of using all map(2)%ngr, we only use the groups for which keep_grp2(g)=1
 ! and map(2)%ngr should be replaced by map(2)%nindepfc !! REALLY ??
 subroutine set_translational_inv_constraints
!! outputs atransl, transl_constraints part of amatrix
!! used for svd and also check the violation of translational invariance constraints
!! once the FC2s are calculated
! no need to change if fc2 is already present, bmat is also always zero
 use svd_stuff
 use params
 use atoms_force_constants
 use ios
 use geometry
 use constants, only : r15
! use force_constants_module
 implicit none
 integer i0,j,k,ired,counter,rnk,res
 integer al,be,ga,de,cnt2,g,t,ti
 real(r15), allocatable:: atemp(:,:),btemp(:),ared1d(:),zero(:)
 logical new

 write(*,*) 'entering translational invce constraints routine with estimate ',transl_constraints
 write(ulog,*)' SETTING TRANSLATIONAL INVARIANCE CONSTRAINTS part of AMATRIX ************'

! set dimensions of ared matrix according to the number of constraints
! exact constraints from translational invariance

! only if include_fc=1 then add the ngroups in the number of columns
! nindepfc = (sum(ngroups(:)*include_fc(:)*(2-include_fc(:))))
! nindepfc = (sum(ngroups(:)))
 dim_al = transl_constraints
! call write_out(ulog,'Estimated # of invariance constraints',inv_constraints)
 call write_out(ulog,'Total # of irreducible FCs=# of columns in Ared',nindepfc)
 allocate( atemp(dim_al,nindepfc),btemp(dim_al),ared1d(nindepfc),zero(nindepfc) )

! As for the mappings: if for a given term t=1,map(2)%nt(g) , we are looking
! for its corresponding atoms, we can find them according to the following
! say rank=2  ;
!   map(2)%gr(g)%iat(1,t) =(0,taui)  map(2)%gr(g)%ixyz(1,t)=al
!   map(2)%gr(g)%iat(2,t) =(nj,tauj)  map(2)%gr(g)%ixyz(2,t)=be
! Really the (n1,tau1)=(iatomcell(:,iatomterm(2,t)),iatomcell0(iatomterm(2,t)))
! We need to sum over the independent terms and find the corresponding FCs
! by using do ti=1,map(2)%ntind(g)  and multiply the indep fcs by map(2)%gr(g)%mat(t,ti)
! In this part supercell information is not needed,
!------------------------------------------------------------------
!       TRANSLATIONAL INVARIANCES
!------------------------------------------------------------------

 atemp = 0d0 ; btemp=0 ; counter = 0 ; zero = 0d0 ; res = 0; !rcut = 0

! include translatnal invce constraints, otherwise calculate deviations from it
 write(ulog,*)' TREATING TRANSLATIONAL INVARIANCES ======================='

 rnk = 1  !********************************************
 if ( include_fc(rnk) .ne. 0 ) then  ! 2 is also OK; needed for violation calculations
!if ( include_fc(rnk) .eq. 1 ) then
! mx=nterms(rnk)
! write(ulog,*)' Rank, alloc size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
! constraint is : sum_i0 phi1_{i0,al} = 0 for all al=1,2,3
!if (mx.ne.0) then
 do al=1,3
    ared1d = zero
    cnt2=0
    gloop1: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors
      if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
      do t=1,map(rnk)%nt(g)
         if ( map(rnk)%gr(g)%ixyz(1,t) .ne. al ) cycle
         do ti=1,map(rnk)%ntind(g)
            ired = cnt2+ti    ! this is the corresponding index of i in ared
            if (ired.gt.nindepfc .or. ired.lt.1) then
               write(ulog,*)'rnk=1:ired=',ired,'> nindepfc=',nindepfc
               stop
            endif
            ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)
         enddo
!     write(ulog,6)'i0,term,red_indx,ared=',iatomterm_1(1,t),t,ired,ared1d(ired)
      enddo
    enddo gloop1

!   if(verbose) write(ulog,11)'Rank1: counter,ared=',counter,ared1d(res+1:size(ared1d))  !res+map(rnk)%ntotind)

! compare with previous lines of ared; write if ared1d not repeated
    if ( ared1d .myeqz. zero ) cycle
    call compare2previous_lines(dim_al,nindepfc,ared1d,atemp,counter,new)
    if (new) then
       counter = counter+1
       atemp(counter,:) = ared1d(:)
       btemp(counter) = 0d0
       if(verbose) write(umatrx,7) ared1d(res+1:size(ared1d)) !res+map(rnk)%ntotind)
       if (counter.eq.dim_al-1) write(ulog,*)' DIM_AL TOO SMALL tr_inv rnk=',rnk
    endif
 enddo
 if(map(rnk)%ngr.gt.0) res = res + map(rnk)%ntotind

 write(ulog,*)' RANK=1; residual index is=',res
 write(ulog,*)' NUMBER OF TRANSLATIONAL CONSTRAINTS of rank=1 is ', counter
 endif

 rnk = 2  !********************************************
 if ( include_fc(rnk) .ne. 0 ) then
!if ( include_fc(rnk) .eq. 1 ) then
! mx=nterms(rnk)
! write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
! constraint is : sum_{n1 tau1} phi2_{tau al,n1 tau1 be} = 0 for all al,be,tau

 do i0=1,natom_prim_cell
 do al=1,3
 do be=al,3  ! avoid redundancies in the constraints

    ared1d = zero
    cnt2=0  ! cnt2 counts the previous independent terms

    gloop2: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors

!! K1 new change made on 2/13/23 ----------------
      if(keep_grp2(g).ne.1) cycle gloop2
!! K1 new change made on 2/13/23 ----------------

      if(g.gt.1) then
         if (keep_grp2(g-1).eq.1) cnt2=cnt2+map(rnk)%ntind(g-1)
      endif
!! K1 new change made on 9/21/23 ----------------
      do t=1,map(rnk)%nt(g)
        if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or.  &
 & map(rnk)%gr(g)%ixyz(1,t) .ne. al .or. map(rnk)%gr(g)%ixyz(2,t) .ne. be ) cycle
        do ti=1,map(rnk)%ntind(g)
          ired = res+cnt2+ti    ! this is the corresponding index of i in ared
!         if(verbose)   write(ulog,*)'TRANS_INVCE: rnk=2 g,ti,ired=',g,ti,ired
          if (ired.gt.nindepfc .or. ired.lt.1) then
             write(ulog,*)'TRANS_INVCE: rnk=2 g,ti,ired=',g,ti,ired,'> nindepfc =',nindepfc
             stop
          endif
          ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)
        enddo
!     write(ulog,6)'term,ti,red_indx,ared=',t,ti,ired,ared1d(ired)
      enddo
    enddo gloop2

!   write(ulog,11)'Rank2: counter,ared=',counter,ared1d(res+1:size(ared1d)) !res+map(rnk)%ntotind)

! compare with previous lines of ared; write if ared1d not repeated
    if ( ared1d .myeqz. zero ) cycle
    call compare2previous_lines(dim_al,nindepfc,ared1d,atemp,counter,new)
    if (new) then
       counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter) = 0d0
       if(verbose) write(umatrx,7) ared1d(res+1:res+ size_kept_fc2) !map(rnk)%ntotind)
       if (counter.eq.dim_al-1) write(ulog,*)' DIM_AL TOO SMALL tr_inv rnk=',rnk
    endif
 enddo
 enddo
 enddo
! res = res + sum(map(rnk)%ntind(:))

!! K1 new change made on 2/13/23 ----------------

! if(map(rnk)%ngr.gt.0) res  = res  + sum(map(rnk  )%ntind(:))
 if(map(rnk)%ngr.gt.0) res  = res  +  size_kept_fc2  ! res counts the cu,mulative # of previous groups of previous rank

!! K1 new change made on 2/13/23 ----------------

   write(ulog,*)' RANK=2; residual index is=',res
   write(ulog,*)' NUMBER OF CUMULATIVE TRANSLATIONAL CONSTRAINTS up to rank 2 is ', counter
 endif

 rnk = 3  !********************************************
 if ( include_fc(rnk) .ne. 0 ) then
!if ( include_fc(rnk) .eq. 1 ) then
! mx=nterms(rnk)
! write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
! constraint is : sum_{k} phi3_{i0 al,j be,k ga} = 0
! for all al,be,ga and i0,j

 do i0=1,natom_prim_cell
 do j=1,natoms ! take only j atoms within the nshells'th neighbor shell of i0
    if ( iatomneighbor(i0,j) .gt. nshells(rnk,i0) ) cycle
 !  if ( iatomneighbor(i0,j) .eq. nshells(rnk,i0) ) then
 !       rcut(rnk) = length(atompos(:,j)-atompos(:,i0))
 !  endif
 do al=1,3
 do be=al,3  ! avoid redundancies in the constraints
 do ga=be,3  ! avoid redundancies in the constraints

    ared1d = zero
    cnt2=0
    gloop3: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors
      if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
      do t=1,map(rnk)%nt(g)
        if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or. map(rnk)%gr(g)%ixyz(1,t) .ne. al &
 &      .or. map(rnk)%gr(g)%iat(2,t) .ne. j  .or. map(rnk)%gr(g)%ixyz(2,t) .ne. be &
 &                                           .or. map(rnk)%gr(g)%ixyz(3,t) .ne. ga ) cycle
        do ti=1,map(rnk)%ntind(g)
          ired = res+cnt2+ti    ! this is the corresponding index of i in ared
          if (ired.gt.nindepfc .or. ired.lt.1) then
             write(ulog,*)'TRANS_INVCE: rnk=3 ired=',ired,'> nindepfc=',nindepfc
             stop
          endif
          ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)
        enddo
!     write(ulog,6)'i0,term,red_indx,ared=',iatomterm_1(1,t),t,ired,ared1d(ired)
      enddo
    enddo gloop3

!   write(ulog,11)'Rank3: counter,ared=',counter,ared1d(res+1:size(ared1d))  !res+map(rnk)%ntotind)

! compare with previous lines of ared; write if ared1d not repeated
    if ( ared1d .myeqz. zero ) cycle
    call compare2previous_lines(dim_al,nindepfc,ared1d,atemp,counter,new)
    if (new) then
       counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter) = 0d0
       if(verbose) write(umatrx,7) ared1d(res+1:res+map(rnk)%ntotind)
       if (counter.eq.dim_al-1) write(ulog,*)' DIM_AL TOO SMALL tr_inv rnk=',rnk
    endif
 enddo
 enddo
 enddo
 enddo
 enddo
 if(map(rnk  )%ngr.gt.0) res  = res  + map(rnk)%ntotind
 write(ulog,*)' RANK=3; residual index is=',res
 write(ulog,*)' NUMBER OF CUMULATIVE TRANSLATIONAL CONSTRAINTS up to rank 3 is ', counter
 endif

 rnk = 4  !********************************************
 if ( include_fc(rnk) .ne. 0 ) then
!if ( include_fc(rnk) .eq. 1 ) then
! mx=nterms(rnk)
! write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
! constraint is : sum_{n3 tau3} phi4_{tau al,n1 tau1 be,n2 tau2 ga,n3 tau3 de}=0
! for all al,be,ga,de and tau n1 tau1 n2 tau2

 do i0=1,natom_prim_cell
 do j=1,natoms ! take only j atoms within the nshells'th neighbor shell of i0
    if ( iatomneighbor(i0,j) .gt. nshells(rnk,i0) ) cycle
 !  if ( iatomneighbor(i0,j) .eq. nshells(rnk,i0) ) then
 !       rcut(rnk) = length(atompos(:,j)-atompos(:,i0))
 !  endif
 do k=1,natoms ! take only j atoms within the nshells'th neighbor shell of i0
    if ( iatomneighbor(i0,k) .gt. nshells(rnk,i0) ) cycle
 do al=1,3
 do be=al,3  ! avoid redundancies in the constraints
 do ga=be,3  ! avoid redundancies in the constraints
 do de=ga,3  ! avoid redundancies in the constraints

    ared1d = zero
    cnt2=0
    gloop4: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors
      if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
      do t=1,map(rnk)%nt(g)
        if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or. map(rnk)%gr(g)%ixyz(1,t) .ne. al &
 &      .or. map(rnk)%gr(g)%iat(2,t) .ne. j  .or. map(rnk)%gr(g)%ixyz(2,t) .ne. be &
 &      .or. map(rnk)%gr(g)%iat(3,t) .ne. k  .or. map(rnk)%gr(g)%ixyz(3,t) .ne. ga &
 &                                           .or. map(rnk)%gr(g)%ixyz(4,t) .ne. de ) cycle
        do ti=1,map(rnk)%ntind(g)
          ired = res+cnt2+ti    ! this is the corresponding index of i in ared
          if (ired.gt.nindepfc .or. ired.lt.1) then
             write(ulog,*)'rnk=4:ired=',ired,'> nindepfc=',nindepfc
             stop
          endif
          ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)
        enddo
!     write(ulog,6)'i0,term,red_indx,ared=',iatomterm_1(1,t),t,ired,ared1d(ired)
      enddo
    enddo gloop4

!   write(ulog,11)'Rank4: counter,ared=',counter,ared1d(res+1:size(ared1d))  !res+map(rnk)%ntotind)

! compare with previous lines of ared; write if ared1d not repeated
    if ( ared1d .myeqz. zero ) cycle
    call compare2previous_lines(dim_al,nindepfc,ared1d,atemp,counter,new)
    if (new) then
       counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter) = 0d0
       if(verbose) write(umatrx,7) ared1d(res+1:size(ared1d))  !res+map(rnk)%ntotind)
       if (counter.eq.dim_al-1) write(ulog,*)' DIM_AL TOO SMALL tr_inv rnk=',rnk
    endif
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 if(map(rnk)%ngr.gt.0) res = res + map(rnk)%ntotind
 write(ulog,*)' RANK=4 res=tot# of independent terms=',res
 write(ulog,*)' NUMBER OF CUMULATIVE TRANSLATIONAL CONSTRAINTS up to rank 4 is ', counter
 endif

 transl_constraints = counter
 write(umatrx,*)'**********************************************'
 write(umatrx,*)' TRANSLATIONAL constr part of amatrix size has ',counter,' lines'
 write(umatrx,*)'**********************************************'
 write(ulog,*)' UPDATED TOTAL NUMBER OF TRANSLATIONAL CONSTRAINTS =', transl_constraints
 allocate(atransl(transl_constraints,nindepfc),btransl(transl_constraints))
 atransl(1:transl_constraints,1:nindepfc)=atemp(1:transl_constraints,1:nindepfc)
! to change if include_fc .ne.0 *****************************
 btransl=btemp
! to change if include_fc .ne.0 *****************************
 deallocate( atemp,ared1d,zero,btemp )
 write(*,*) 'exiting translational invce constraints routine'

6 format(a,3(1x,i3),66(2x,f7.3))
7 format(266(1x,f7.3))
11 format(a,1(1x,i3),66(2x,f7.3))

 end subroutine set_translational_inv_constraints
!===========================================================================
 subroutine set_rotational_inv_constraints
!! outputs arot, brot, rot_constraints
!! used for svd and also check the violation of rotational invariance constraints
!! once the FC2s are calculated
 use svd_stuff
 use params
 use atoms_force_constants
 use ios
 use geometry
 use constants, only : r15
! use force_constants_module
 implicit none
 integer i0,j,k,l,t,ired,counter,rnk,res,res1,g,ti,cnt2,cnt3
 integer al,be,ga,de,mu,nu
 integer lc(3,3,3)
 real(r15), allocatable:: ared1d(:),zero(:),aux(:)
 real(r15), allocatable:: atemp(:,:),btemp(:)
 real(r15) junk
 logical new

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!* BEWARE: for include_fc=2; you may need to add res to cnt3 when defining ired
!* res being from the previous rank
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 write(*,*) 'entering rotational invce constraints routine with estimate ',rot_constraints
 write(ulog,*)' SETTING ROTATIONAL INVARIANCE CONSTRAINTS  part of AMATRIX  ******************'

! set dimensions of ared matrix according to the number of constraints
!------------------------------------------------------------------
!       ROTATIONAL INVARIANCES
!------------------------------------------------------------------

 dim_al=rot_constraints  ! # of lines from estimate
 allocate( atemp(dim_al,nindepfc),btemp(dim_al),ared1d(nindepfc),zero(nindepfc) )
! initialize the Levi-Civita symbol here (called lc)
 lc = 0; counter = 0; zero = 0d0; atemp = 0d0; btemp = 0d0; junk = 0d0
 lc(1,2,3) =  1 ; lc(2,3,1) =  1 ; lc(3,1,2) =  1
 lc(1,3,2) = -1 ; lc(2,1,3) = -1 ; lc(3,2,1) = -1

 res = 0 ; res1 = 0
 rnk = 1  !********************************************
 if ( include_fc(rnk) .ne. 0 ) then  ! sum of the torques is zero
!   mx=nterms(rnk)
!   write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
! constraint is : sum_i0 phi1_{i0,al} .cross. X_{i0,be} = 0
   do ga=1,3
      ared1d = zero
      cnt2=0
      gloop1: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors
        if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
!         write(*,*)'g,cnt2=',g,cnt2
        do t=1,map(rnk)%nt(g)
          i0 = map(rnk)%gr(g)%iat(1,t)   ! should be in the primitive cell
!         write(*,*)'ga,t,i0=',ga,t,i0

!! ADDED LINE BY K1 on JULY 10th &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

          if(i0 .ne. iatomcell0(i0) ) cycle  ! restrict to atoms in the primitive cell

!! ADDED LINE BY K1 on JULY 10th &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

          al = map(rnk)%gr(g)%ixyz(1,t)
          do ti=1,map(rnk)%ntind(g)
            ired = cnt2+ti    ! this is the corresponding index of i in ared
            if (ired.gt.nindepfc .or. ired.lt.1) then
               write(ulog,*)'rnk=1:ired=',ired,'> nindepfc=',nindepfc
               stop
            endif
            do be=1,3
              ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)*atompos(be,i0)*lc(al,be,ga)
            enddo
          enddo
!          if(verbose) write(ulog,5)'i0,g,term,ired,ared=',i0,g,t,ired,ared1d(ired)
        enddo
      enddo gloop1
!      write(ulog,11)'al,ared=',al,ared1d

! compare with previous lines of ared; write if ared1d not repeated
      if ( ared1d .myeqz. zero ) cycle
      call compare2previous_lines(dim_al,nindepfc,ared1d,atemp,counter,new)
      if (new) then
         counter = counter+1; atemp(counter,:) = ared1d(:);
         if(verbose) write(umatrx,7) ared1d(res+1:size(ared1d))  !res+map(rnk)%ntotind),btemp(counter)
         if (counter.eq.dim_al-1) write(ulog,*)'DIM_AL TOO SMALL rot_inv rnk=',rnk
         if(include_fc(2).ne.2) then
            btemp(counter) = 0d0
         else
            btemp(counter)=0 !!! .... to be completed
! constraint is for each i0,al,nu: phi1_{i0,al} eps(al,be,nu)= - sum_j phi_{i0,j) phi_{i0,j}^{be,al} eps(al,ga,nu) Rj^ga















         endif
      endif
   enddo
   if(map(rnk)%ngr.gt.0) res = res + map(rnk)%ntotind
   write(ulog,*)' RANK=1; residual index is=',res
   write(ulog,*)' NUMBER OF CUMULATIVE ROTATIONAL CONSTRAINTS up to rank 1 is ', counter
 endif

 rnk = 2  !********************************************
 if ( include_fc(rnk) .ne. 0 ) then ! .and. include_fc(1) .ne. 0 ) then
!   mx=nterms(rnk)
   allocate(aux(ngroups(1)))
!   write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
! constraint is : sum_{n1 tau1} phi2_{tau al,n1 tau1 be} X_{n1 tau1 be'} +
! phi1_{tau be} delta_{al be'} symm wrt be<-> be'  for all tau,al,be,be'

   do i0=1,natom_prim_cell
   do al=1,3
   do de=1,3
      ared1d = zero
      cnt2=0
      junk=0

      if(include_fc(2).eq.1) then
      gloop12: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors

!! K1 new change made on 2/13/23 ----------------
      if(keep_grp2(g).ne.1) cycle gloop12
        if(g.gt.1) then
           if (keep_grp2(g-1).eq.1) cnt2=cnt2+map(rnk)%ntind(g-1)
        endif
!! K1 new change made on 9/21/23 ----------------
        do t=1,map(rnk)%nt(g)
          if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or. map(rnk)%gr(g)%ixyz(1,t) .ne. al ) cycle
          j  = map(rnk)%gr(g)%iat(2,t) ;       be = map(rnk)%gr(g)%ixyz(2,t)
          do ti=1,map(rnk)%ntind(g)
            ired = res+cnt2+ti    ! this is the corresponding index of i in ared
            if (ired.gt.nindepfc .or. ired.lt.1) then
               write(ulog,*)'rnk=2:ired=',ired,'> nindepfc=',nindepfc
               stop
            endif
            do ga=1,3
              ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)*atompos(ga,j)*lc(be,ga,de)
            enddo
          enddo
!          if(verbose) write(ulog,5)'i0,g,t,ired,ared=',i0,g,t,ired,ared1d(ired)
        enddo
      enddo gloop12

   elseif (include_fc(2) .eq. 2) then
!if (include_fc(2) .eq. 2) then
! bring this piece to the RHS( bmat)
!btemp(counter)=-sum(ared1d)
!endif
      gloop13: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors
         if(keep_grp2(g).ne.1) cycle gloop13
         if(g.gt.1) then
            if (keep_grp2(g-1).eq.1) cnt2=cnt2+map(rnk)%ntind(g-1)
         endif
!! K1 new change made on 9/21/23 ----------------
        do t=1,map(rnk)%nt(g)
          if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or. map(rnk)%gr(g)%ixyz(1,t) .ne. al ) cycle
          j  = map(rnk)%gr(g)%iat(2,t) ;       be = map(rnk)%gr(g)%ixyz(2,t)
          do ti=1,map(rnk)%ntind(g)
            ired=cnt2+ti
          !  ired = res+cnt2+ti    ! this is the corresponding index of i in ared
            if (ired.gt.nindepfc .or. ired.lt.1) then
               write(ulog,*)'rnk=2:ired=',ired,'> nindepfc =',nindepfc
               stop
            endif
            do ga=1,3
               junk=junk-fc_ind(ired)*atompos(ga,j)*lc(be,ga,de) !UNCOMMENT !map(rnk)%gr(g)%mat(t,ti)*atompos(ga,j)*lc(be,ga,de)
              !ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)*atompos(ga,j)*lc(be,ga,de)
            enddo
          enddo
          if(verbose) write(ulog,5)'i0,g,t,ired,ared=',i0,g,t,ired,ared1d(ired)
        enddo
      enddo gloop13


   endif


! TO CHECK------------
!   ared1d = -ared1d
! TO CHECK------------
      aux=0; !junk = 0   ! junk goes to bmat (it already includes the - sign)
      cnt3 = 0
      if(include_fc(1) .eq. 1 .or. include_fc(1) .eq. 2) then


!!
!! why is bikash adding the condition .or. include_fc(1) .eq. 2 to the above?
!!


      do g=1,map(rnk-1)%ngr  ! ineq. terms and identify neighbors
           if(g.gt.1) then
              if (keep_grp2(g-1).eq.1) cnt3=cnt3+map(rnk)%ntind(g-1)
           endif
!! K1 new change made on 9/21/23 ----------------
!       if(g.gt.1) cnt3 = cnt3 + map(rnk-1)%ntind(g-1)  ! this is wrong it should be rnk not rnk-1 !
        do t=1,map(rnk-1)%nt(g)
           if ( map(rnk-1)%gr(g)%iat(1,t) .ne.i0 ) cycle
           be = map(rnk-1)%gr(g)%ixyz(1,t)
           do ti=1,map(rnk-1)%ntind(g)
              ired = cnt3+ti    ! this is the corresponding index of i in ared
            if (ired.gt.nindepfc .or. ired.lt.1) then
               write(ulog,*)'rnk=2 ired=',ired,'> nindepfc=',nindepfc
               stop
            endif
              ared1d(ired) = ared1d(ired) + map(rnk-1)%gr(g)%mat(t,ti)*lc(be,al,de)
           enddo
        enddo
      enddo
      endif
!     if(include_fc(1) .eq. 2) then  ! put it on the right hand side in b
!       do g=1,map(rnk-1)%ngr       ! ngroups(1)
!       do t=1,map(rnk-1)%nt(g)        !   t = map_1(j)
!          junk = junk - aux(igroup_1(t))*fcs_1(igroup_1(t))
!       enddo
!     endif
! now set ared1d properly
!     if(include_fc(1) .eq. 1) then
   !      ared1d( 1:ngroups(1) ) = aux
!     endif
      if ( ared1d .myeqz. zero ) cycle
      call compare2previous_lines(dim_al,nindepfc,ared1d,atemp,counter,new)
      if (new) then
         counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter)=junk
         if(verbose) write(umatrx,7) ared1d(res+1:res+size_kept_fc2) !size(ared1d))!res+map(rnk)%ntotind),btemp(counter)
       if (counter.eq.dim_al-1) write(ulog,*)'DIM_AL TOO SMALL rot_inv rnk=',rnk
      endif
   enddo
   enddo
   enddo
   deallocate(aux)

!! K1 new change made on 2/13/23 ----------------

   if(include_fc(2) .eq. 1) then
! if(map(rnk  )%ngr.gt.0) res  = res  + sum(map(rnk  )%ntind(:))
     if(map(rnk  )%ngr.gt.0) res  = res  +  size_kept_fc2
     if(map(rnk-1)%ngr.gt.0) res1 = res1 + sum(map(rnk-1)%ntind(:))
   elseif(include_fc(2) .eq. 2) then
      if(map(rnk-1)%ngr.gt.0) then
         write(*,*) "RNK and RES1 are: ", rnk, res1
         !res1 = res1 + map(1)%ntotind(:)
      endif
   endif

!! K1 new change made on 2/13/23 ----------------

   write(ulog,*)' RANK=2; residual index is=',res
   write(ulog,*)' NUMBER OF CUMULATIVE ROTATIONAL CONSTRAINTS up to rank 2 is ', counter
 endif

 rnk = 3  !********************************************
 if ( include_fc(rnk) .ne. 0 ) then
!   mx=nterms(rnk)
   allocate(aux(map(2)%ntotind))
!!    allocate(aux(ngroups(2)))
!!   write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
   write(ulog,*)' allocated size of aux(:) is=',map(2)%ntotind

   do i0=1,natom_prim_cell
   do j =1,natoms
      if ( iatomneighbor(i0,j) .gt. nshells(rnk,i0) ) cycle
   do al=1,3
   do be=1,3
   do nu=1,3
      ared1d = zero
      cnt2 = 0
      gloop33: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors
        if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
        do t=1,map(rnk)%nt(g)
          if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or. map(rnk)%gr(g)%ixyz(1,t) .ne. al ) cycle
          if ( map(rnk)%gr(g)%iat(2,t) .ne. j  .or. map(rnk)%gr(g)%ixyz(2,t) .ne. be ) cycle
          k  = map(rnk)%gr(g)%iat(3,t) ;       ga = map(rnk)%gr(g)%ixyz(3,t)
          do ti=1,map(rnk)%ntind(g)
            ired = res+cnt2+ti    ! this is the corresponding index of i in ared
            if (ired.gt.nindepfc .or. ired.lt.1) then
               write(ulog,*)'ROT_INVCE: rnk=3 ired=',ired,'> nindepfc=',nindepfc
               stop
            endif
            do de=1,3
              ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)*atompos(de,k)*lc(ga,de,nu)
            enddo
          enddo
!          if(verbose) write(ulog,5)'i0,g,t,ired,ared=',i0,g,t,ired,ared1d(ired)
        enddo
      enddo gloop33

! TO CHECK------------
!   ared1d = -ared1d
! TO CHECK------------
      aux = 0; junk = 0   ! junk goes to bmat (it already includes the - sign)
      cnt3 = 0
      if(include_fc(2) .eq. 1 ) then !.or. include_fc(2) .eq. 2) then
      do g=1,map(rnk-1)%ngr  ! ineq. terms and identify neighbors
        if(g.gt.1) cnt3 = cnt3 + map(rnk-1)%ntind(g-1)
        do t=1,map(rnk-1)%nt(g)
           if ( map(rnk-1)%gr(g)%iat(1,t) .ne.i0 .or. map(rnk-1)%gr(g)%iat(2,t) .ne.j ) cycle
           do ti=1,map(rnk-1)%ntind(g)
              ired = res1+cnt3+ti    ! this is the corresponding index of i in ared
            if (ired.gt.nindepfc .or. ired.lt.1) then
               write(ulog,*)'rnk=3 ired=',ired,'> nindepfc=',nindepfc
               stop
            endif
              if ( map(rnk-1)%gr(g)%ixyz(2,t) .eq. be ) then
                 ga = map(rnk-1)%gr(g)%ixyz(1,t)
                 ared1d(ired) = ared1d(ired) + map(rnk-1)%gr(g)%mat(t,ti)*lc(ga,al,nu)
              elseif ( map(rnk-1)%gr(g)%ixyz(1,t) .eq. al ) then
                 ga = map(rnk-1)%gr(g)%ixyz(2,t)
                 ared1d(ired) = ared1d(ired) + map(rnk-1)%gr(g)%mat(t,ti)*lc(ga,be,nu)
              endif
           enddo
        enddo
      enddo
      endif

      if(include_fc(2) .eq. 2 ) then ! So in this case include_fc(2) equals 1 or 2 are separated see exactly above
    !     call read_fcs_2(2)
!         write(*,*) "value of map(rnk-1)%ngr is: ", map(rnk-1)%ngr
           do g=1,map(rnk-1)%ngr  ! ineq. terms and identify neighbors
             if(g.gt.1) cnt3 = cnt3 + map(rnk-1)%ntind(g-1)
             do t=1,map(rnk-1)%nt(g)
                if ( map(rnk-1)%gr(g)%iat(1,t) .ne.i0 .or. map(rnk-1)%gr(g)%iat(2,t) .ne.j ) cycle
                do ti=1,map(rnk-1)%ntind(g)
                   ired = cnt3+ti    ! this is the corresponding index of i in ared
         !          ired = res1+ cnt3+ti    ! this is the corresponding index of i in ared ! was this before
!                 write(*,*) "value of ired is: ", ired
                 if (ired.gt.nindepfc .or. ired.lt.1) then
                    write(ulog,*)'rnk=3 ired=',ired,'> nindepfc=',nindepfc
                    stop
                 endif
                   if ( map(rnk-1)%gr(g)%ixyz(2,t) .eq. be ) then
                      ga = map(rnk-1)%gr(g)%ixyz(1,t)
                      aux(ired) = aux(ired) + map(rnk-1)%gr(g)%mat(t,ti)*lc(ga,al,nu)
                     ! ared1d(ired) = ared1d(ired) + map(rnk-1)%gr(g)%mat(t,ti)*lc(ga,al,nu) ! was this before
                   elseif ( map(rnk-1)%gr(g)%ixyz(1,t) .eq. al ) then
                      ga = map(rnk-1)%gr(g)%ixyz(2,t)
                       aux(ired) = aux(ired) + map(rnk-1)%gr(g)%mat(t,ti)*lc(ga,be,nu)
                     ! ared1d(ired) = ared1d(ired) + map(rnk-1)%gr(g)%mat(t,ti)*lc(ga,be,nu) ! was like this before
                   endif
                enddo
             enddo
           enddo
!           write(*,*) "Junk before is given by this: ", junk
           junk=junk-dot_product(aux,fc_ind) ! This is added here!make sure fc_ind is global
!           write(*,*) "Junk after this is given by: ", junk
           endif

!     if(include_fc(2) .eq. 2) then  ! put it on the right hand side in b
!       do k=1,ngroups(2)
!          t = map_2(k)
!          ired = igroup_2(t)
!          junk = junk - aux(ired)*fcs_2(igroup_2(t))
!       enddo
!     endif
!     if(include_fc(2) .eq. 1) then
   !      ared1d( res-ngroups(2)+1:res ) = aux
!     endif
      if ( ared1d .myeqz. zero ) cycle
      call compare2previous_lines(dim_al,nindepfc,ared1d,atemp,counter,new)
      if (new) then
         counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter)=junk
         if(verbose) write(umatrx,7) ared1d(res+1:res+map(rnk)%ntotind),btemp(counter)
       if(verbose) write(ulog,5)'ROT:new term for:i0,j,al,be,nu=',i0,j,al,be,real(nu)
       if (counter.eq.dim_al-1) write(ulog,*)'DIM_AL TOO SMALL rot_inv rnk=',rnk
      endif
   enddo
   enddo
   enddo
   enddo
   enddo
   deallocate(aux)

!! K1 new change made on 2/13/23 ----------------

   if(map(rnk)%ngr.gt.0) res = res + map(rnk)%ntotind
   if(map(rnk-1)%ngr.gt.0) res1  = res1  +  size_kept_fc2

!! K1 new change made on 2/13/23 ----------------

   write(ulog,*)' RANK=3; residual index is=',res
   write(ulog,*)' NUMBER OF CUMULATIVE ROTATIONAL CONSTRAINTS up to rank 3 is ', counter
 endif

 rnk = 4  !********************************************
 if ( include_fc(rnk) .ne. 0 ) then
!   mx=nterms(rnk)
   allocate(aux(ngroups(3)))
!   write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx

   do i0=1,natom_prim_cell
   do j =1,natoms
      if ( iatomneighbor(i0,j) .gt. nshells(rnk,i0) ) cycle
   do k =1,natoms
      if ( iatomneighbor(i0,k) .gt. nshells(rnk,i0) ) cycle
   do al=1,3
   do be=1,3
   do ga=1,3
   do nu=1,3
      ared1d = zero
      cnt2 = 0
      gloop14: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors
        if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
        do t=1,map(rnk)%nt(g)
          if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or. map(rnk)%gr(g)%ixyz(1,t) .ne. al ) cycle
          if ( map(rnk)%gr(g)%iat(2,t) .ne. j  .or. map(rnk)%gr(g)%ixyz(2,t) .ne. be ) cycle
          if ( map(rnk)%gr(g)%iat(3,t) .ne. k  .or. map(rnk)%gr(g)%ixyz(3,t) .ne. ga ) cycle
          l  = map(rnk)%gr(g)%iat(4,t) ;       de = map(rnk)%gr(g)%ixyz(4,t)
          do ti=1,map(rnk)%ntind(g)
            ired = res+cnt2+ti    ! this is the corresponding index of i in ared
            if (ired.gt.nindepfc .or. ired.lt.1) then
               write(ulog,*)'rnk=4:ired=',ired,'> nindepfc=',nindepfc
               stop
            endif
            do mu=1,3
              ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)*atompos(mu,l)*lc(de,mu,nu)
            enddo
          enddo
        enddo
      enddo gloop14

! TO CHECK------------
!   ared1d = -ared1d
! TO CHECK------------
      aux = 0; junk = 0   ! junk goes to bmat (it already includes the - sign)
      cnt3 = 0
      if(include_fc(3) .eq. 1 .or. include_fc(3) .eq. 2) then
      do g=1,map(rnk-1)%ngr  ! ineq. terms and identify neighbors
        if(g.gt.1) cnt3 = cnt3 + map(rnk-1)%ntind(g-1)
        do t=1,map(rnk-1)%nt(g)
           if ( map(rnk-1)%gr(g)%iat(1,t) .ne.i0 .or. map(rnk-1)%gr(g)%iat(2,t) .ne.j &
       &   .or. map(rnk-1)%gr(g)%iat(3,t) .ne.k ) cycle
           do ti=1,map(rnk-1)%ntind(g)
              ired = res1+cnt3+ti    ! this is the corresponding index of i in ared
            if (ired.gt.nindepfc .or. ired.lt.1) then
               write(ulog,*)'rnk=4 ired=',ired,'> nindepfc=',nindepfc
               stop
            endif
              if ( map(rnk-1)%gr(g)%ixyz(2,t) .eq. be .and. map(rnk-1)%gr(g)%ixyz(3,t) .eq. ga ) then
                 de = map(rnk-1)%gr(g)%ixyz(1,t)
                 ared1d(ired) = ared1d(ired) + map(rnk-1)%gr(g)%mat(t,ti)*lc(de,al,nu)
              elseif ( map(rnk-1)%gr(g)%ixyz(1,t) .eq. al .and. map(rnk-1)%gr(g)%ixyz(3,t) .eq. ga ) then
                 de = map(rnk-1)%gr(g)%ixyz(2,t)
                 ared1d(ired) = ared1d(ired) + map(rnk-1)%gr(g)%mat(t,ti)*lc(de,be,nu)
              elseif ( map(rnk-1)%gr(g)%ixyz(1,t) .eq. al .and. map(rnk-1)%gr(g)%ixyz(2,t) .eq. be ) then
                 de = map(rnk-1)%gr(g)%ixyz(3,t)
                 ared1d(ired) = ared1d(ired) + map(rnk-1)%gr(g)%mat(t,ti)*lc(de,ga,nu)
              endif
           enddo
        enddo
      enddo
      endif

!     if(include_fc(3) .eq. 2) then  ! put it on the right hand side in b
!       do l=1,ngroups(3)
!          t = map_3(l)
!          ired = igroup_3(t)
!          junk = junk - aux(ired)*fcs_3(igroup_3(t))
!       enddo
!     endif
!     if(include_fc(3) .eq. 1) then
   !      ared1d( res-ngroups(3)+1:res ) = aux
!     endif
      if ( ared1d .myeqz. zero ) cycle
      call compare2previous_lines(dim_al,nindepfc,ared1d,atemp,counter,new)
      if (new) then
         counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter)=junk
         if(verbose) write(umatrx,7) ared1d(res+1:res+map(rnk)%ntotind),btemp(counter)
       if (counter.eq.dim_al-1) write(ulog,*)'DIM_AL TOO SMALL rot_inv rnk=',rnk
      endif
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo
   enddo
   deallocate(aux)
   if(map(rnk)%ngr.gt.0) res = res + map(rnk)%ntotind
   write(ulog,*)' RANK=4; residual index is=',res
   write(ulog,*)' NUMBER OF CUMULATIVE ROTATIONAL CONSTRAINTS up to rank 4 is ', counter
 endif

 rot_constraints = counter
 write(ulog,*)' UPDATED TOTAL NUMBER OF ROTATIONAL INVARIANCE CONSTRAINTS =', rot_constraints
 allocate(arot(rot_constraints,nindepfc),brot(rot_constraints))
 arot(1:rot_constraints,1:nindepfc)=atemp(1:counter,1:nindepfc)
 brot(1:rot_constraints)=btemp(1:counter)

 if(res.ne.nindepfc) then
    write(ulog,*)' res=',res,' must be equal to # of indepdt FCs ',nindepfc
    write(*   ,*)' res=',res,' must be equal to # of indepdt FCs ',nindepfc
    stop
 endif
  write(umatrx,*)'**********************************************'
  write(umatrx,*)' ROTATIONAL constr part of amatrix size has ',counter,' lines'
  write(umatrx,*)'**********************************************'

 deallocate ( atemp,btemp,ared1d,zero )
 write(*,*) 'exiting rotational invce constraints routine'

3 format(a,2(1x,i6),66(1x,i3))
4 format(a,66(1x,i3))
5 format(a,4(1x,i3),66(2x,f7.3))
6 format(a,3(1x,i3),66(2x,f7.3))
7 format(266(1x,f7.3))
8 format(i8,1x,i8,3x,f25.15,3x,f25.15)
9 format(i8,1x,99(1x,f8.4))
11 format(a,1(1x,i3),66(2x,f7.3))
13 format(a,3(1x,i3),66(2x,f7.3))
15 format(a,5(1x,i3),66(2x,f7.3))
17 format(a,7(1x,i3),66(2x,f7.3))

 end subroutine set_rotational_inv_constraints
!======================================================================================
 subroutine set_huang_inv_constraints
!! The constraint is: sum_Rtautau' phi(tau,al,r+tau'be)(r+tau'-tau)_ga (r+tau'-tau)_de
!! is invariant under(al,be) <-> (ga,de)
 use params
 use atoms_force_constants
 use svd_stuff
 use ios
 use constants, only : r15
! use geometry
! use force_constants_module
 implicit none
 integer i,j,al,be,ga,de,t,ti,mterm,cnt2,g,rnk,voigt,vab,vcd,counter,ired,res,dim_h
 real(r15) huang,rr(3)
 real(r15), allocatable:: ared1d(:),zero(:)

 huang_constraints=15  ! this is only for rank 2

   rnk = 2;
!  res = sum(map(1)%ntind(:))
   if(map(1)%ngr.gt.0) then
     res = map(1)%ntotind
   else
     res=0
   endif
   write(ulog,*)' SETTING HUANG INVARIANCE CONSTRAINTS  part of AMATRIX *****************'
   dim_h = huang_constraints
   write(ulog,*)' allocating ',dim_h,' lines and ',nindepfc,' columns'
   write(*   ,*)' allocating ',dim_h,' lines and ',nindepfc,' columns'
   allocate( ahuang(dim_h,nindepfc),ared1d(nindepfc),zero(nindepfc),bhuang(dim_h) )
   counter = 0; zero = 0d0 ; ahuang = 0d0

   do al=1,3
   do be=al,3
   do ga=1,3
   do de=ga,3
      vab = voigt(al,be)
      vcd = voigt(ga,de)
      if ( vab.le.vcd ) cycle
      huang = 0
      mterm = 0
      ared1d = zero  !! K1 :: this comes after the al,be,ga,de loop !!
      atomloop: do i=1,natom_prim_cell
      cnt2=0     ! cnt2 counts the previous independent terms
      gloop: do g=1,map(2)%ngr

!! K1 new change made on 2/13/23 ----------------
      if(keep_grp2(g).ne.1) cycle gloop
!! K1 new change made on 2/13/23 ----------------

      if(g.gt.1) then
         if (keep_grp2(g-1).eq.1) cnt2=cnt2+map(rnk)%ntind(g-1)
      endif
!! K1 new change made on 9/21/23 ----------------
      do t=1,map(2)%nt(g)
         if ( map(2)%gr(g)%iat(1,t).ne.i ) cycle
         if ( map(2)%gr(g)%ixyz(1,t).eq.al .and. map(2)%gr(g)%ixyz(2,t).eq.be ) then
            mterm = mterm+1
            j  =  map(2)%gr(g)%iat(2,t)
!           tauj=iatomcell0(j)
!           nj=iatomcell(:,j)
            rr = atompos(:,j)-atompos(:,i)
            do ti =1,map(2)%ntind(g)
               ired = res+cnt2+ti  ! the corresponding index of i in ared
               ared1d(ired) = ared1d(ired) + rr(ga)*rr(de)*map(2)%gr(g)%mat(t,ti)
!               huang = huang + rr(ga)*rr(de)*map(2)%gr(g)%mat(t,ti)*fcs(res+cnt2+ti)
            enddo
         elseif ( map(2)%gr(g)%ixyz(1,t).eq.ga .and. map(2)%gr(g)%ixyz(2,t).eq.de ) then
            mterm = mterm+1
            j  =  map(2)%gr(g)%iat(2,t)
            rr = atompos(:,j)-atompos(:,i)
            do ti =1,map(2)%ntind(g)
               ired = res+cnt2+ti  ! the corresponding index of i in ared
               ared1d(ired) = ared1d(ired) - rr(al)*rr(be)*map(2)%gr(g)%mat(t,ti)
            enddo
         endif
      enddo
      enddo gloop
      enddo atomloop

      counter = counter+1
      if(verbose) write(ulog,6)'counter,mterm,vab,vcd=',counter,mterm,al,be,ga,de,vab,vcd
      if (counter.eq.dim_h+1) then
         write(ulog,*)' DIM_H TOO SMALL, huang_inv=15, counter= ',counter
         write(*   ,*)' DIM_H TOO SMALL, huang_inv=15, counter= ',counter
         stop
      endif
      ahuang(counter,:) = ared1d(:)
      if(verbose) write(umatrx,7) ared1d(res+1:res+size_kept_fc2) !map(2)%ntotind)
      bhuang(counter)=0d0
   enddo
   enddo
   enddo
   enddo

   deallocate( ared1d,zero )
  write(umatrx,*)'**********************************************'
  write(umatrx,*)' HUANG constr part of amatrix size has ',counter,' lines'
  write(umatrx,*)'**********************************************'

6 format(a,8(1x,i4),9(2x,f7.3))
7 format(266(1x,f9.3))

 end subroutine set_huang_inv_constraints
!===========================================================================
 subroutine set_force_displacement_matrix(ncfgs,frc_constr,afrc,bfrc)
!! outputs afrc and bfrc matrices afrc,bfrc from the knowledge of arrays force and displ for a given supercell
!! equilibrium positions already subtracted from displ
 use svd_stuff
 use params
 use atoms_force_constants
 use ios
 use geometry
 use constants, only : r15
 implicit none
 integer, intent(in) :: frc_constr,ncfgs
 real(r15), intent(out) :: afrc(frc_constr,nindepfc),bfrc(frc_constr)
 integer l,t,ired,counter,rnk,nat,confg,jatom,katom,latom,res,cnt2
 integer al,be,ga,de,taui,tauj,tauk,taul,ni(3),nj(3),nk(3),nl(3),ti
 real(r15), allocatable:: ared1d(:),aux(:)

 write(*,*) 'entering set_force_displacements routine,nindepfc=',nindepfc
 write(ulog,*)' SETUP_AMATRICES : ************************************'
! l is group index, t: general term index, and ti: independent term index
! this convention is followed in all setting constraint or amat subroutines
! set dimensions of ared matrix according to the number of constraints
!------------------------------------------------------------------
!       FORCE-DISPLACEMENT CONSTRAINTS
!------------------------------------------------------------------

! allocate(aforce(force_constraints,ngr),bforce(force_constraints))
 allocate( ared1d(nindepfc) )
 afrc = 0d0; bfrc = 0d0

5 format(a,9(i5))
 write(ulog,6)'********************************************************'
 write(ulog,6)' Now making the force-displacement part of the A matrix'
 write(ulog,*)' Force_constraints      =', frc_constr
! write(ulog,*)' number of configs,3N_SC=',size(nlines),nlines
!if(frc_constr.ne.sum(nlines)) then
!   write(*,*)'check parameters frc_constr natom_super_cell ncfgs=',frc_constr,natom_super_cell,ncfgs
!   stop
!endif

! if(frc_constr.ne.3*natom_super_cell*ncfgs) then
!    write(*,*)'check parameters frc_constr natom_super_cell ncfgs=',frc_constr,natom_super_cell,ncfgs
!    write(*,*)'are your supercels of different size??'
! endif
! now compute the force constraints on the FCs. The equation is:
! -F_ia = phi1_ia + phi2_{ia,jb} u_jb + 0.5 phi3_{ia,jb,kc} u_jb u_kc + ...
!   iatomterm(1,t)=(0,tau)   iatomterm(2,t)=(n1,tau1)
!   ixyzterm (1,t)=alpha     ixyzterm (2,t)=beta
! Really the (n1,tau1)=(iatomcell(:,iatomterm(2,t)),iatomcell0(iatomterm(2,t)))
  counter = 0
  do confg=1,ncfgs
  do nat=1,natom_super_cell
  do al=1,3             ! this is how we order each line

        counter = counter+1      ! this is basically the line number
        bfrc(counter) = -force(al,nat,confg)
! find the corresponding atom in the primitive cell
        taui = atom_sc(nat)%cell%tau
        ni   = atom_sc(nat)%cell%n
!       write(ulog,3)'line#,confg,i,alfa,tau,n =',counter,confg,nat,al,taui,ni

! now set up ared1d: rank by rank
        ared1d = 0d0 ; res=0
        rnk=1  !********************************************
        if ( include_fc(rnk) .ne. 0 ) then

!           allocate(aux(ndindp(rnk))) ! size of aux =  # of independent FCs for that rank
           allocate(aux(map(rnk)%ntotind))  ! size of aux =  # of independent FCs for that rank
           cnt2= 0        ! cnt2 is the size of previous groups
           aux = 0
           do l=1,map(rnk)%ngr  ! sum over groups
              if(l.gt.1) cnt2 = cnt2 + map(rnk)%ntind(l-1)
              do t=1,map(rnk)%nt(l) ! sum over all terms in that group
                 if ( taui.eq. map(rnk)%gr(l)%iat(1,t) .and.  &
             &      al  .eq. map(rnk)%gr(l)%ixyz(1,t) ) then
                    do ti=1,map(rnk)%ntind(l)
! this is the index of the indep FC coming in the A*FC=b matrix product
                       ired = cnt2+ti
!                       write(ulog,5)'l,cnt2,ti,ired=',l,cnt2,ti,ired
                       aux(ired) = map(rnk)%gr(l)%mat(t,ti)
                    enddo
                 endif
              enddo
           enddo

           if ( include_fc(rnk) .eq. 1 ) then
!              ared1d(res+1:res+ndindp(rnk)) = aux(1:ndindp(rnk))
!              res = res + ndindp(rnk)
              ared1d(res+1:res+map(rnk)%ntotind) = aux(1:map(rnk)%ntotind)
              res = res + map(rnk)%ntotind
!          elseif ( include_fc(rnk) .eq. 2 ) then
!             do j=1,map(rnk)%ngr
!                t = map_1(j)
!                bfrc(counter)=bfrc(counter)-fcs_1(igroup_1(t))*aux(igroup_1(t))
!             enddo
           endif
           deallocate(aux)
        endif

        rnk=2  !********************************************
        if ( include_fc(rnk) .ne. 0 ) then

! the force on atom i is from j in the supercell and all of its images for which the corresponding FC needs to be identified
!           allocate(aux(ndindp(rnk)))
           allocate(aux(size_kept_fc2)) !map(rnk)%ntotind))
! for given (ni,taui,al) find its neighbors within nshells(2)
           cnt2= 0        ! cnt2 is the size of previous groups
           aux = 0
           groups: do l=1,map(rnk)%ngr  ! sum over groups

!! K1 new change made on 2/13/23 ----------------
           if(keep_grp2(l).ne.1) cycle groups
!! K1 new change made on 2/13/23 ----------------
           if(l.gt.1) then
              if (keep_grp2(l-1).eq.1) cnt2=cnt2+map(rnk)%ntind(l-1)
           endif
!! K1 new change made on 9/21/23 ----------------
           l8: do t=1,map(rnk)%nt(l) ! sum over all terms in that group
              if ( taui.eq. map(rnk)%gr(l)%iat(1,t) .and.  &
          &        al  .eq. map(rnk)%gr(l)%ixyz(1,t) ) then
                 tauj = iatomcell0(map(rnk)%gr(l)%iat(2,t))
                 nj(:)= iatomcell(:,map(rnk)%gr(l)%iat(2,t)) + ni(:) ! translate by ni
                 be   = map(rnk)%gr(l)%ixyz(2,t)
! Identify neighbor j within the SCell, and its images, find its displacement and add to ared
                 call findatom_sc(nj,tauj,jatom)
                 if (jatom.eq.0) then
                    write(ulog,4)'SET_ARED:jatom not found: tau,n ',tauj,nj
                    write(ulog,4)'for rank,term ',rnk,t
                    write(ulog,4)'atom,xyz,cnfg ',nat,al,confg
                    stop
                 endif
    !! K1 NEW --------------------------
  !!!  !            rij=length(atompos(:,nat)-atompos(:,jatom))
  !!!  !            if (rij.gt.radius(rnk)) cycle l8
! need to find how many images of j are connected to i with the same force constant.
    !! K1 NEW --------------------------
! --------------    Aug 6th 2018 --------------------
! j is now identified within the supercell, but we need to add the FCs of
! all of its images if the range is larger than the supercell size
! So need to translate nj by the supercell translation vectors, subtract ni, and identify the
! corresponding FCs and add them to aux
                 do ti=1,map(rnk)%ntind(l)
! this is the index of the indep FC coming in the A*FC=b matrix product
                    ired = cnt2+ti ! this is the corresponding index of aux
!                   write(ulog,5)'l,cnt2,ti,ired=',l,cnt2,ti,ired
! this ads all equivalent FCs of that group or shell (ie atoms on edges or vertices of the WS of supercell
                    aux(ired) = aux(ired) + displ(be,jatom,confg)*map(rnk)%gr(l)%mat(t,ti)
                 enddo
  !             write(ulog,22)'term,red_indx,ared=',t,iatomterm_2(1,t),  &
  ! &    ixyzterm_2(1,t),iatomterm_2(2,t),ixyzterm_2(2,t),ired,ared1d(ired)
              endif
           enddo l8
           enddo groups

           if ( include_fc(rnk) .eq. 1 ) then

!! K1 new change made on 2/13/23 ----------------

!              ared1d(res+1:res+map(rnk)%ntotind) = aux(1:map(rnk)%ntotind)
!              res = res + map(rnk)%ntotind
               ared1d(res+1:res+size_kept_fc2) = aux(1:size_kept_fc2)
               res  = res  +  size_kept_fc2

!! K1 new change made on 2/13/23 ----------------

           elseif ( include_fc(rnk) .eq. 2 ) then

               cnt2=0
               do l=1,map(rnk)%ngr
                  if(keep_grp2(l).ne.1) cycle
                  if(l.gt.1) then
                     if (keep_grp2(l-1).eq.1) cnt2=cnt2+map(rnk)%ntind(l-1)
                  endif
!! K1 new change made on 9/21/23 ----------------
               do t=1,map(rnk)%nt(l) ! sum over all terms in that group
                  if ( taui.eq. map(rnk)%gr(l)%iat (1,t) .and.  &
            &          al  .eq. map(rnk)%gr(l)%ixyz(1,t) )         then
                     tauj = iatomcell0(map(rnk)%gr(l)%iat(2,t))
                     nj(:)= iatomcell(:,map(rnk)%gr(l)%iat(2,t)) + ni(:) ! translate by ni
                     be   = map(rnk)%gr(l)%ixyz(2,t)
  ! Identify neighbor j within the SCell, and its images, find its displacement and add to ared
                     call findatom_sc(nj,tauj,jatom)
                     if (jatom.eq.0) then
                        write(ulog,4)'SET_ARED:jatom not found: tau,n ',tauj,nj
                        write(ulog,4)'for rank,term ',rnk,t
                        write(ulog,4)'atom,xyz,cnfg ',nat,al,confg
                        stop
                     endif
                     do ti=1,map(rnk)%ntind(l)
                        ired = cnt2+ti ! this is the corresponding index of aux
                        bfrc(counter)=bfrc(counter)-fcrnk(rnk,ired)*displ(be,jatom,confg)
                     enddo
                  endif
               enddo
               enddo

           endif

!              res = res + ndindp(rnk)
!          elseif ( include_fc(rnk) .eq. 2 ) then
!             do l=1,map(rnk)%ngr
!               t = map_2(l)
!               bfrc(counter)=bfrc(counter)-fcs_2(igroup_2(t))*aux(igroup_2(t))
! shouldn't this line simply be: bf = bf - fcs_2(j) * aux(j) ??
!             enddo
           deallocate(aux)
        endif

        rnk=3  !********************************************
        if ( include_fc(rnk) .ne. 0 ) then

!           allocate(aux(ndindp(rnk)))
           allocate(aux(map(rnk)%ntotind))
! for given (ni,taui,al) find its two neighbors within nshells(3)
           cnt2= 0        ! cnt2 is the size of previous groups
           aux = 0
           do l=1,map(rnk)%ngr  ! sum over groups
              if(l.gt.1) cnt2 = cnt2 + map(rnk)%ntind(l-1)
           do t=1,map(rnk)%nt(l) ! sum over all terms in that group
              if ( taui.eq. map(rnk)%gr(l)%iat(1,t) .and.  &
              &    al  .eq. map(rnk)%gr(l)%ixyz(1,t) ) then
                 tauj = iatomcell0(map(rnk)%gr(l)%iat(2,t))
                 nj(:)= iatomcell(:,map(rnk)%gr(l)%iat(2,t)) + ni(:)
                 be   = map(rnk)%gr(l)%ixyz(2,t)
                 tauk = iatomcell0(map(rnk)%gr(l)%iat(3,t))
                 nk(:)= iatomcell(:,map(rnk)%gr(l)%iat(3,t)) + ni(:)
                 ga   = map(rnk)%gr(l)%ixyz(3,t)
! Identify neighbor j within the SCell, find its displacement and add to ared
                 call findatom_sc(nj,tauj,jatom)
                 if (jatom.eq.0) then
                    write(ulog,4)'SET_ARED:jatom not found: tau,n ',tauj,nj
                    write(ulog,4)'for rank,term ',rnk,t
                    write(ulog,4)'atom,xyz,cnfg ',nat,al,confg
                    stop
                 endif
                 call findatom_sc(nk,tauk,katom)
                 if (katom.eq.0) then
                    write(ulog,4)'SET_ARED:katom not found: tau,n ',tauk,nk
                    write(ulog,4)'for rank,term ',rnk,t
                    write(ulog,4)'atom,xyz,cnfg ',nat,al,confg
                    stop
                 endif
    !! K1 NEW --------------------------
    !!           rij=length(atompos(:,nat)-atompos(:,jatom))
    !!           rik=length(atompos(:,nat)-atompos(:,katom))
    !!           rjk=length(atompos(:,jatom)-atompos(:,katom))
    !!   ! cycle if two of the pair distances is larger than radius
    !!           if (rij.gt.radius(rnk)) then
    !!              if (rik.gt.radius(rnk) .or. rjk.gt.radius(rnk) ) cycle
    !!           else
    !!              if (rik.gt.radius(rnk) .and. rjk.gt.radius(rnk) ) cycle
    !!           endif
    !! K1 NEW --------------------------
                 do ti=1,map(rnk)%ntind(l)
! this is the index of the indep FC coming in the A*FC=b matrix product
                    ired = cnt2+ti ! this is the corresponding index of aux
!                   write(ulog,5)'l,cnt2,ti,ired=',l,cnt2,ti,ired
                    aux(ired) = aux(ired) + map(rnk)%gr(l)%mat(t,ti) *  &
  &                 displ(be,jatom,confg)*displ(ga,katom,confg)/2d0
                 enddo
              endif
           enddo
           enddo

           if ( include_fc(rnk) .eq. 1 ) then
!              ared1d(res+1:res+ndindp(rnk)) = aux(1:ndindp(rnk))
!              res = res + ndindp(rnk)
              ared1d(res+1:res+map(rnk)%ntotind) = aux(1:map(rnk)%ntotind)
              res = res + map(rnk)%ntotind
           endif
           deallocate(aux)
        endif

        rnk=4  !********************************************
        if ( include_fc(rnk) .ne. 0 ) then

!           allocate(aux(ndindp(rnk)))
           allocate(aux(map(rnk)%ntotind))
           cnt2= 0        ! cnt2 is the size of previous groups
           aux = 0
           do l=1,map(rnk)%ngr  ! sum over groups
              if(l.gt.1) cnt2 = cnt2 + map(rnk)%ntind(l-1)
           do t=1,map(rnk)%nt(l) ! sum over all terms in that group
              if ( taui.eq. map(rnk)%gr(l)%iat(1,t) .and.  &
              &    al  .eq. map(rnk)%gr(l)%ixyz(1,t) ) then
                 tauj = iatomcell0(map(rnk)%gr(l)%iat(2,t))
                 nj(:)= iatomcell(:,map(rnk)%gr(l)%iat(2,t)) + ni(:)
                 be   = map(rnk)%gr(l)%ixyz(2,t)
                 tauk = iatomcell0(map(rnk)%gr(l)%iat(3,t))
                 nk(:)= iatomcell(:,map(rnk)%gr(l)%iat(3,t)) + ni(:)
                 ga   = map(rnk)%gr(l)%ixyz(3,t)
                 taul = iatomcell0(map(rnk)%gr(l)%iat(4,t))
                 nl(:)= iatomcell(:,map(rnk)%gr(l)%iat(4,t)) + ni(:)
                 de   = map(rnk)%gr(l)%ixyz(4,t)
! Identify neighbor j within the SCell, find its displacement and add to ared
                 call findatom_sc(nj,tauj,jatom)
                 if (jatom.eq.0) then
                    write(ulog,4)'SET_ARED:jatom not found: tau,n ',tauj,nj
                    write(ulog,4)'for rank,term ',rnk,t
                    write(ulog,4)'atom,xyz,cnfg ',nat,al,confg
                    stop
                 endif
                 call findatom_sc(nk,tauk,katom)
                 if (katom.eq.0) then
                    write(ulog,4)'SET_ARED:katom not found: tau,n ',tauk,nk
                    write(ulog,4)'for rank,term ',rnk,t
                    write(ulog,4)'atom,xyz,cnfg ',nat,al,confg
                    stop
                 endif
                 call findatom_sc(nl,taul,latom)
                 if (latom.eq.0) then
                    write(ulog,4)'SET_ARED:katom not found: tau,n ',taul,nl
                    write(ulog,4)'for rank,term ',rnk,t
                    write(ulog,4)'atom,xyz,cnfg ',nat,al,confg
                    stop
                 endif
    !! K1 NEW --------------------------
    !!           rij=length(atompos(:,nat)-atompos(:,jatom))
    !!           rik=length(atompos(:,nat)-atompos(:,katom))
    !!           rjk=length(atompos(:,jatom)-atompos(:,katom))
    !!   ! cycle if two of the pair distances is larger than radius
    !!           if (rij.gt.radius(rnk)) then
    !!              if (rik.gt.radius(rnk) .or. rjk.gt.radius(rnk) ) cycle
    !!           else
    !!              if (rik.gt.radius(rnk) .and. rjk.gt.radius(rnk) ) cycle
    !!           endif
    !! K1 NEW --------------------------
                 do ti=1,map(rnk)%ntind(l)
! this is the index of the indep FC coming in the A*FC=b matrix product
                    ired = cnt2+ti ! this is the corresponding index of aux
!                   write(ulog,5)'l,cnt2,ti,ired=',l,cnt2,ti,ired
                    aux(ired) = aux(ired) + map(rnk)%gr(l)%mat(t,ti) *  &
  &  displ(be,jatom,confg)*displ(ga,katom,confg)*displ(de,latom,confg)/6d0
                 enddo
              endif
           enddo
           enddo

           if ( include_fc(rnk) .eq. 1 ) then
!              ared1d(res+1:res+ndindp(rnk)) = aux(1:ndindp(rnk))
!              res = res + ndindp(rnk)
              ared1d(res+1:res+map(rnk)%ntotind) = aux(1:map(rnk)%ntotind)
              res = res + map(rnk)%ntotind
           endif
           deallocate(aux)
        endif

! now write into the ared matrix ( we assume force terms are not redundant
! otherwise after the SVD, w matrix will have very small or zero elements)
        afrc(counter,:) = ared1d(:)

!        if(verbose) write(umatrx,7)(afrc(counter,l),l=1,ngr),bfrc(counter)
!make sure the index ired is not repeated for each rank and is increased
!from a given rank to a higher rank

 enddo
 enddo
 enddo

 write(umatrx,*)'**********************************************'
 write(umatrx,*)' Force-displ constraints part of amatrix size is ',counter
 write(umatrx,*)'**********************************************'
 deallocate(ared1d)

 write(*,*) 'exiting set_force_displacements routine'

4 format(a,4(1x,i6))
6 format(a,3(1x,i3),66(2x,f7.3))
7 format(266(1x,g9.2))
12 format(a,2(1x,i3),66(2x,g11.4))
22 format(a,6(1x,i3),66(2x,f9.5))

 end subroutine set_force_displacement_matrix
