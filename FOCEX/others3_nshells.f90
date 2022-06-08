!=============================================================================
 subroutine set_translational_inv_constraints
!! sets the translational invariance part of the amatrx
 use svd_stuff
 use params
 use atoms_force_constants
 use ios
 use geometry
 use force_constants_module
 implicit none
 integer i0,j,k,ired,counter,rnk,res
 integer al,be,ga,de,cnt2,g,t,ti
 real(8), allocatable:: atemp(:,:),btemp(:),ared1d(:),zero(:)
 logical new
 if(itrans .ne. 1) return
 write(*,*) 'entering translational invce constraints routine'
 write(ulog,*)' SETTING TRANSLATIONAL INVARIANCE CONSTRAINTS *****************'

! set dimensions of ared matrix according to the number of constraints
! exact constraints from translational invariance

! only if include_fc=1 then add the ngroups in the number of columns
! ngr = (sum(ngroups(:)*include_fc(:)*(2-include_fc(:))))
! ngr = (sum(ngroups(:)))
 dim_al = inv_constraints
 call write_out(ulog,'Estimated # of invariance constraints',inv_constraints)
 call write_out(ulog,'Total # of irreducible FCs=# of columns in Ared',ngr)
 allocate( atemp(dim_al,ngr),btemp(dim_al),ared1d(ngr),zero(ngr) )

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

 atemp = 0d0 ; counter = 0; zero = 0d0 ; res = 0; rcut = 0

! include translatnal invce constraints, otherwise calculate deviations from it
 write(ulog,*)' TREATING TRANSLATIONAL INVARIANCES ======================='

 rnk = 1  !------------------------------------------
 if ( include_fc(rnk) .ne. 0 ) then
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
          if (ired.gt.ngr .or. ired.lt.1) then
             write(ulog,*)'TI:rnk=1:ired=',ired,'>ngr=',ngr
             stop
          endif
          ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)
       enddo
!     write(ulog,6)'i0,term,red_indx,ared=',iatomterm_1(1,t),t,ired,ared1d(ired)
    enddo
    enddo gloop1

!   write(ulog,11)'al,ared=',al,ared1d

! compare with previous lines of ared; write if ared1d not repeated
    if ( ared1d .myeqz. zero ) cycle
    call compare2previous_lines(dim_al,ngr,ared1d,atemp,counter,new)
    if (new) then
       counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter) = 0d0
!      write(umatrx,7)(atemp(counter,dum),dum=1,ngr)
       if (counter.eq.dim_al-1) write(ulog,*)' DIM_AL TOO SMALL tr_inv rnk=',rnk
    endif
 enddo
 if(map(rnk)%ngr.gt.0) res = res + sum(map(rnk)%ntind(:)) !
!endif
 write(ulog,*)' RANK=1; residual index is=',res
 endif

 rnk = 2  !------------------------------------------
!if ( include_fc(rnk) .ne. 0 ) then
if ( include_fc(rnk) .eq. 1 ) then
! mx=nterms(rnk)
! write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
! constraint is : sum_{n1 tau1} phi2_{tau al,n1 tau1 be} = 0 for all al,be,tau

 do i0=1,natoms0
 do al=1,3
 do be=al,3  ! avoid redundancies in the constraints

    ared1d = zero
    cnt2=0  ! cnt2 counts the previous independent terms
    gloop2: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors
      if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
      do t=1,map(rnk)%nt(g)
        if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or.  &
 & map(rnk)%gr(g)%ixyz(1,t) .ne. al .or. map(rnk)%gr(g)%ixyz(2,t) .ne. be ) cycle
        do ti=1,map(rnk)%ntind(g)
          ired = res+cnt2+ti    ! this is the corresponding index of i in ared
          if (ired.gt.ngr .or. ired.lt.1) then
             write(ulog,*)'TI:rnk=2:ired=',ired,'>ngr=',ngr
             stop
          endif
          ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)
        enddo
!     write(ulog,6)'term,ti,red_indx,ared=',t,ti,ired,ared1d(ired)
      enddo
    enddo gloop2

! compare with previous lines of ared; write if ared1d not repeated
    if ( ared1d .myeqz. zero ) cycle
    call compare2previous_lines(dim_al,ngr,ared1d,atemp,counter,new)
    if (new) then
       counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter) = 0d0
!      write(umatrx,7)(atemp(counter,dum),dum=1,ngr)
       if (counter.eq.dim_al-1) write(ulog,*)' DIM_AL TOO SMALL tr_inv rnk=',rnk
    endif
 enddo
 enddo
 enddo
! res = res + sum(map(rnk)%ntind(:))
 if(map(rnk  )%ngr.gt.0 .and. include_fc(2).eq.1) res  = res  + sum(map(rnk  )%ntind(:)) !
 write(ulog,*)' RANK=2; residual index is=',res
 endif

 rnk = 3  !------------------------------------------
 if ( include_fc(rnk) .ne. 0 ) then
!if ( include_fc(rnk) .eq. 1 ) then
! mx=nterms(rnk)
! write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
! constraint is : sum_{k} phi3_{i0 al,j be,k ga} = 0
! for all al,be,ga and i0,j

 do i0=1,natoms0
 do j=1,natoms ! take only j atoms within the nshells'th neighbor shell of i0
!   if ( iatomneighbor(i0,j) .gt. nshells(rnk) ) cycle
!   if ( iatomneighbor(i0,j) .eq. nshells(rnk) ) then
    if ( iatomneighbor(i0,j) .gt. nshells(rnk,i0) ) cycle
    if ( iatomneighbor(i0,j) .eq. nshells(rnk,i0) ) then
         rcut(rnk) = length(atompos(:,j)-atompos(:,i0))
    endif
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
          if (ired.gt.ngr .or. ired.lt.1) then
             write(ulog,*)'TI:rnk=3:ired=',ired,'>ngr=',ngr
             stop
          endif
          ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)
        enddo
!     write(ulog,6)'i0,term,red_indx,ared=',iatomterm_1(1,t),t,ired,ared1d(ired)
      enddo
    enddo gloop3
!   write(ulog,15)'i0,j,al,be,ga,ared=',i0,j,al,be,ga,ared1d

! compare with previous lines of ared; write if ared1d not repeated
    if ( ared1d .myeqz. zero ) cycle
    call compare2previous_lines(dim_al,ngr,ared1d,atemp,counter,new)
    if (new) then
       counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter) = 0d0
!      write(umatrx,7)(atemp(counter,dum),dum=1,ngr)
       if (counter.eq.dim_al-1) write(ulog,*)' DIM_AL TOO SMALL tr_inv rnk=',rnk
    endif
 enddo
 enddo
 enddo
 enddo
 enddo
! res = res + sum(map(rnk)%ntind(:))
 if(map(rnk  )%ngr.gt.0) res  = res  + sum(map(rnk  )%ntind(:))
 write(ulog,*)' RANK=3; residual index is=',res
 write(ulog,*)' RANK=3 cnter for tr cnstr=',counter
 endif

 rnk = 4  !------------------------------------------
 if ( include_fc(rnk) .ne. 0 ) then
!if ( include_fc(rnk) .eq. 1 ) then
! mx=nterms(rnk)
! write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
! constraint is : sum_{n3 tau3} phi4_{tau al,n1 tau1 be,n2 tau2 ga,n3 tau3 de}=0
! for all al,be,ga,de and tau n1 tau1 n2 tau2

 do i0=1,natoms0
 do j=1,natoms ! take only j atoms within the nshells'th neighbor shell of i0
!   if ( iatomneighbor(i0,j) .gt. nshells(rnk) ) cycle
!   if ( iatomneighbor(i0,j) .eq. nshells(rnk) ) then
    if ( iatomneighbor(i0,j) .gt. nshells(rnk,i0) ) cycle
    if ( iatomneighbor(i0,j) .eq. nshells(rnk,i0) ) then
         rcut(rnk) = length(atompos(:,j)-atompos(:,i0))
    endif
 do k=1,natoms ! take only j atoms within the nshells'th neighbor shell of i0
!   if ( iatomneighbor(i0,k) .gt. nshells(rnk) ) cycle
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
          if (ired.gt.ngr .or. ired.lt.1) then
             write(ulog,*)'TI:rnk=4:ired=',ired,'>ngr=',ngr
             stop
          endif
          ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)
        enddo
!     write(ulog,6)'i0,term,red_indx,ared=',iatomterm_1(1,t),t,ired,ared1d(ired)
      enddo
    enddo gloop4
!   write(ulog,17)'i0,j,k,al,be,ga,de,ared=',i0,j,k,al,be,ga,de,ared1d

! compare with previous lines of ared; write if ared1d not repeated
    if ( ared1d .myeqz. zero ) cycle
    call compare2previous_lines(dim_al,ngr,ared1d,atemp,counter,new)
    if (new) then
       counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter) = 0d0
!      write(umatrx,7)(atemp(counter,dum),dum=1,ngr)
       if (counter.eq.dim_al-1) write(ulog,*)' DIM_AL TOO SMALL tr_inv rnk=',rnk
    endif
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
!res = res + sum(map(rnk)%ntind(:))
 if(map(rnk  )%ngr.gt.0) res  = res  + sum(map(rnk  )%ntind(:))
 write(ulog,*)' RANK=4 res=tot# of groups=',res
 write(ulog,*)' RANK=4 cnter for tr cnstr=',counter
 endif

 transl_constraints = counter
 write(ulog,*)' NUMBER OF TRANSLATIONAL CONSTRAINTS =', transl_constraints
 allocate(atransl(transl_constraints,ngr))
 atransl(1:transl_constraints,1:ngr)=atemp(1:transl_constraints,1:ngr)
 deallocate( atemp,ared1d,zero )
 write(*,*) 'exiting translational invce constraints routine'

!6 format(a,3(1x,i3),66(2x,f7.3))
!7 format(66(1x,f7.3))
!11 format(a,1(1x,i3),66(2x,f7.3))

 end subroutine set_translational_inv_constraints
!===========================================================================
 subroutine set_rotational_inv_constraints
!! sets the rotational invariance part of the amatrx
 use svd_stuff
 use params
 use atoms_force_constants
 use ios
 use geometry
 use force_constants_module
 implicit none
 integer i0,j,k,l,t,ired,counter,rnk,res,res1,g,ti,cnt2,cnt3
 integer al,be,ga,de,mu,nu
 integer lc(3,3,3)
 real(8), allocatable:: ared1d(:),zero(:),aux(:)
 real(8), allocatable:: atemp(:,:),btemp(:)
 real(8) junk
 logical new
if(irot .ne. 1) return
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! BEWARE: for include_fc=2; you may need to add res to cnt3 when defining ired
! res being from the previous rank
!-----------------------------------------------------------------
!-----------------------------------------------------------------

 write(*,*) 'entering rotational invce constraints routine'
 write(ulog,*)' SETTING ROTATIONAL INVARIANCE CONSTRAINTS  ******************'

! set dimensions of ared matrix according to the number of constraints
!------------------------------------------------------------------
!       ROTATIONAL INVARIANCES
!------------------------------------------------------------------

 allocate( atemp(dim_al,ngr),btemp(dim_al),ared1d(ngr),zero(ngr) )
 write(ulog,*)' TREATING ROTATIONAL INVARIANCES ======================='
! initialize the Levi-Civita symbol here (called lc)
 lc = 0; counter = 0; zero = 0d0; atemp = 0d0; btemp = 0d0; junk = 0d0
 lc(1,2,3) =  1 ; lc(2,3,1) =  1 ; lc(3,1,2) =  1
 lc(1,3,2) = -1 ; lc(2,1,3) = -1 ; lc(3,2,1) = -1

 res = 0 ; res1 = 0
 rnk = 1  !------------------------------------------
 !if ( include_fc(rnk) .ne. 0 ) then
   if ( include_fc(rnk) .eq. 1 ) then  ! sum of the torques is zero
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

! ADDED LINE BY K1 on JULY 10th &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

          if(i0 .ne. iatomcell0(i0) ) cycle  ! restrict to atoms in the primitive cell

! ADDED LINE BY K1 on JULY 10th &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          al = map(rnk)%gr(g)%ixyz(1,t)
          do ti=1,map(rnk)%ntind(g)
            ired = cnt2+ti    ! this is the corresponding index of i in ared
            if (ired.gt.ngr .or. ired.lt.1) then
               write(ulog,*)'rnk=1:ired=',ired,'>ngr=',ngr
               stop
            endif
            do be=1,3
              ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)*atompos(be,i0)*lc(al,be,ga)
            enddo
          enddo
          if(verbose) write(ulog,5)'i0,g,term,ired,ared=',i0,g,t,ired,ared1d(ired)
        enddo
      enddo gloop1
!      write(ulog,11)'al,ared=',al,ared1d

! compare with previous lines of ared; write if ared1d not repeated
      if ( ared1d .myeqz. zero ) cycle
      call compare2previous_lines(dim_al,ngr,ared1d,atemp,counter,new)
      if (new) then
         counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter) = 0d0
!        write(umatrx,7)(atemp(counter,dum),dum=1,ngr),btemp(counter)
         if (counter.eq.dim_al-1) write(ulog,*)'DIM_AL TOO SMALL rot_inv rnk=',rnk
      endif
   enddo
   if(map(rnk)%ngr.gt.0) res = res + sum(map(rnk)%ntind(:))
   write(ulog,*)' RANK=1; residual index is=',res
   write(ulog,*)' RANK=1 cnter for rot cntr=',counter
 endif

 rnk = 2  !------------------------------------------
 if ( include_fc(rnk) .ne. 0 ) then ! .and. include_fc(1) .ne. 0 ) then
!   mx=nterms(rnk)
   allocate(aux(ngroups(1)))
!   write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx
! constraint is : sum_{n1 tau1} phi2_{tau al,n1 tau1 be} X_{n1 tau1 be'} +
! phi1_{tau be} delta_{al be'} symm wrt be<-> be'  for all tau,al,be,be'

   do i0=1,natoms0
   do al=1,3
   do de=1,3
      ared1d = zero
      cnt2=0
      junk=0;

      if (include_fc(2) .eq. 1) then
      gloop12: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors
        if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
        do t=1,map(rnk)%nt(g)
          if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or. map(rnk)%gr(g)%ixyz(1,t) .ne. al ) cycle
          j  = map(rnk)%gr(g)%iat(2,t) ;       be = map(rnk)%gr(g)%ixyz(2,t)
          do ti=1,map(rnk)%ntind(g)
            ired = res+cnt2+ti    ! this is the corresponding index of i in ared
            if (ired.gt.ngr .or. ired.lt.1) then
               write(ulog,*)'rnk=2:ired=',ired,'>ngr=',ngr
               stop
            endif
            do ga=1,3
              ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)*atompos(ga,j)*lc(be,ga,de) 
            enddo
          enddo
          if(verbose) write(ulog,5)'i0,g,t,ired,ared=',i0,g,t,ired,ared1d(ired)
        enddo
      enddo gloop12
   elseif (include_fc(2) .eq. 2) then
!if (include_fc(2) .eq. 2) then
! bring this piece to the RHS( bmat)
!btemp(counter)=-sum(ared1d)
!endif
!      call read_fcs_2(2)
      gloop13: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors
        if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
        do t=1,map(rnk)%nt(g)
          if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or. map(rnk)%gr(g)%ixyz(1,t) .ne. al ) cycle
          j  = map(rnk)%gr(g)%iat(2,t) ;       be = map(rnk)%gr(g)%ixyz(2,t)
          do ti=1,map(rnk)%ntind(g)
            ired=cnt2+ti
          !  ired = res+cnt2+ti    ! this is the corresponding index of i in ared
            if (ired.gt.ngr .or. ired.lt.1) then
               write(ulog,*)'rnk=2:ired=',ired,'>ngr=',ngr
               stop
            endif
            do ga=1,3
               junk=junk-fc_ind(ired)*atompos(ga,j)*lc(be,ga,de) !UNCOMMENT      !map(rnk)%gr(g)%mat(t,ti)*atompos(ga,j)*lc(be,ga,de) 
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
      aux=0;! junk = 0   ! junk goes to bmat (it already includes the - sign)
      cnt3 = 0
      if(include_fc(1) .eq. 1 .or. include_fc(1) .eq. 2 .or. include_fc(2) .eq. 2) then    ! what if I add include_fc(2) .eq. 2 here? ...
      do g=1,map(rnk-1)%ngr  ! ineq. terms and identify neighbors
        if(g.gt.1) cnt3 = cnt3 + map(rnk-1)%ntind(g-1)
        do t=1,map(rnk-1)%nt(g)
           if ( map(rnk-1)%gr(g)%iat(1,t) .ne.i0 ) cycle
           be = map(rnk-1)%gr(g)%ixyz(1,t)
           do ti=1,map(rnk-1)%ntind(g)
              ired = cnt3+ti    ! this is the corresponding index of i in ared
            if (ired.gt.ngr .or. ired.lt.1) then
               write(ulog,*)'rnk=2 ired=',ired,'>ngr=',ngr
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
      call compare2previous_lines(dim_al,ngr,ared1d,atemp,counter,new)
      if (new) then
         counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter)=junk
!        write(umatrx,7)(atemp(counter,dum),dum=1,ngr),btemp(counter)
       if (counter.eq.dim_al-1) write(ulog,*)'DIM_AL TOO SMALL rot_inv rnk=',rnk
      endif
   enddo
   enddo
   enddo
   deallocate(aux)
   if(include_fc(2) .eq. 1) then
      if (map(rnk)%ngr .gt. 0) then
!         write(*,*) "RNK and RES1 are: ", rnk, res1
         !res1=res1+map(1)%ntotind(:)
      endif
      if (map(rnk-1)%ngr .gt. 0) then
!         write(*,*) "RNK and RES1 are: ", rnk, res1
         !res1 = res1 + map(1)%ntotind(:)
      endif
   elseif(include_fc(2) .eq. 2) then
      if(map(rnk-1)%ngr.gt.0) then
!         write(*,*) "RNK and RES1 are: ", rnk, res1
         !res1 = res1 + map(1)%ntotind(:)
      endif
   endif
   write(ulog,*)' RI: RANK=2; residual index is=',res
   write(ulog,*)' RANK=2 cnter for rot cntr=',counter
 endif

 rnk = 3  !------------------------------------------
 if ( include_fc(rnk) .ne. 0 ) then
!   mx=nterms(rnk)
   allocate(aux(map(2)%ntotind))
 !  allocate(aux(ngroups(2))) ! It was this before
  write(ulog,*)' allocated size of aux is=',map(2)%ntotind

   do i0=1,natoms0
   do j =1,natoms
!      if ( iatomneighbor(i0,j) .gt. nshells(rnk) ) cycle
      if ( iatomneighbor(i0,j) .gt. nshells(rnk,i0) ) cycle
   do al=1,3
   do be=1,3
   do nu=1,3
      ared1d = zero
      cnt2 = 0
      gloop33: do g=1,map(rnk)%ngr  ! ineq. terms and identify neighbors !two gloop13 here defined previously need to change previous one
        if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
        do t=1,map(rnk)%nt(g)
          if ( map(rnk)%gr(g)%iat(1,t) .ne. i0 .or. map(rnk)%gr(g)%ixyz(1,t) .ne. al ) cycle
          if ( map(rnk)%gr(g)%iat(2,t) .ne. j  .or. map(rnk)%gr(g)%ixyz(2,t) .ne. be ) cycle
          k  = map(rnk)%gr(g)%iat(3,t) ;       ga = map(rnk)%gr(g)%ixyz(3,t)
          do ti=1,map(rnk)%ntind(g)
            ired = res+cnt2+ti    ! this is the corresponding index of i in ared
            if (ired.gt.ngr .or. ired.lt.1) then
               write(ulog,*)'rnk=3:ired=',ired,'>ngr=',ngr
               stop
            endif
            do de=1,3
              ared1d(ired) = ared1d(ired) + map(rnk)%gr(g)%mat(t,ti)*atompos(de,k)*lc(ga,de,nu)
            enddo
          enddo
          if(verbose) write(ulog,5)'i0,g,t,ired,ared=',i0,g,t,ired,ared1d(ired)
        enddo
      enddo gloop33

! TO CHECK------------
!   ared1d = -ared1d
! TO CHECK------------
      aux = 0; junk = 0   ! junk goes to bmat (it already includes the - sign)
      cnt3 = 0
      if(include_fc(2) .eq. 1 ) then
    !  if(include_fc(2) .eq. 1 .or. include_fc(2) .eq. 2) then ! It was like this before
      do g=1,map(rnk-1)%ngr  ! ineq. terms and identify neighbors
        if(g.gt.1) cnt3 = cnt3 + map(rnk-1)%ntind(g-1)
        do t=1,map(rnk-1)%nt(g)
           if ( map(rnk-1)%gr(g)%iat(1,t) .ne.i0 .or. map(rnk-1)%gr(g)%iat(2,t) .ne.j ) cycle
           do ti=1,map(rnk-1)%ntind(g)
              ired = res1+cnt3+ti    ! this is the corresponding index of i in ared
            if (ired.gt.ngr .or. ired.lt.1) then
               write(ulog,*)'rnk=3 ired=',ired,'>ngr=',ngr
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
!write(*,*) "Value of include_fc(2) is: ", include_fc(2)
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
                 if (ired.gt.ngr .or. ired.lt.1) then
                    write(ulog,*)'rnk=3 ired=',ired,'>ngr=',ngr
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
           junk=junk-dot_product(aux,fc_ind) ! This is added here make sure fc_ind is declared global
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
      call compare2previous_lines(dim_al,ngr,ared1d,atemp,counter,new)
      if (new) then
         counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter)=junk
!        write(umatrx,7)(atemp(counter,dum),dum=1,ngr),btemp(counter)
       if(verbose) write(ulog,5)'ROT:new term for:i0,j,al,be,nu=',i0,j,al,be,real(nu)
       if (counter.eq.dim_al-1) write(ulog,*)'DIM_AL TOO SMALL rot_inv rnk=',rnk
      endif
   enddo
   enddo
   enddo
   enddo
   enddo
   deallocate(aux)
   if(map(rnk  )%ngr.gt.0) res  = res  + sum(map(rnk  )%ntind(:))
   if(map(rnk-1)%ngr.gt.0 .and. include_fc(2).eq.1) res1 = res1 + sum(map(rnk-1)%ntind(:))
   write(ulog,*)' RANK=3; residual index is=',res
   write(ulog,*)' RANK=3 cnter for rot cntr=',counter
 endif

 rnk = 4  !------------------------------------------
 if ( include_fc(rnk) .ne. 0 ) then
!   mx=nterms(rnk)
   allocate(aux(ngroups(3)))
!   write(ulog,*)' allocated size of igroup arrays, mx=nterms(rnk) is=',rnk,mx

   do i0=1,natoms0
   do j =1,natoms
!     if ( iatomneighbor(i0,j) .gt. nshells(rnk) ) cycle
      if ( iatomneighbor(i0,j) .gt. nshells(rnk,i0) ) cycle
   do k =1,natoms
!     if ( iatomneighbor(i0,k) .gt. nshells(rnk) ) cycle
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
            if (ired.gt.ngr .or. ired.lt.1) then
               write(ulog,*)'rnk=4:ired=',ired,'>ngr=',ngr
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
            if (ired.gt.ngr .or. ired.lt.1) then
               write(ulog,*)'rnk=4 ired=',ired,'>ngr=',ngr
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
      call compare2previous_lines(dim_al,ngr,ared1d,atemp,counter,new)
      if (new) then
         counter = counter+1; atemp(counter,:) = ared1d(:); btemp(counter)=junk
!        write(umatrx,7)(atemp(counter,dum),dum=1,ngr),btemp(counter)
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
 ! res = res + sum(map(rnk)%ntind(:))
   if(map(rnk  )%ngr.gt.0) res  = res  + sum(map(rnk  )%ntind(:))
   write(ulog,*)' RANK=4; residual index is=',res
   write(ulog,*)' RANK=4 cnter for rot cntr=',counter
 endif

 rot_constraints = counter
 write(ulog,*)' NUMBER OF ROTATIONAL CONSTRAINTS =', rot_constraints
 allocate(arot(rot_constraints,ngr),brot(rot_constraints))
 arot(1:rot_constraints,1:ngr)=atemp(1:counter,1:ngr)
 brot(1:rot_constraints)=btemp(1:counter)

 write(ulog,*)' res=',res,' must be equal to # of indepdt FCs ',ngr
 inv_constraints = rot_constraints + transl_constraints
 write(ulog,*)' TOTAL NUMBER OF INVARIANCE CONSTRAINTS =', inv_constraints
 write(umatrx,*)'**********************************************'
 write(umatrx,*)' sym constr part of amatrix size has ',inv_constraints,' lines'
 write(umatrx,*)'**********************************************'

 deallocate ( atemp,btemp,ared1d,zero )
 write(*,*) 'exiting rotational invce constraints routine'

!3 format(a,2(1x,i6),66(1x,i3))
!4 format(a,66(1x,i3))
5 format(a,4(1x,i3),66(2x,f7.3))
!6 format(a,3(1x,i3),66(2x,f7.3))
!7 format(66(1x,f7.3))
!8 format(i8,1x,i8,3x,f25.15,3x,f25.15)
!9 format(i8,1x,99(1x,f8.4))
!11 format(a,1(1x,i3),66(2x,f7.3))
!13 format(a,3(1x,i3),66(2x,f7.3))
!15 format(a,5(1x,i3),66(2x,f7.3))
!17 format(a,7(1x,i3),66(2x,f7.3))

 end subroutine set_rotational_inv_constraints
!===========================================================================
 subroutine set_force_displacement_matrix(frc_constr,afrc,bfrc)
 use svd_stuff
 use params
 use atoms_force_constants
 use ios
 use geometry
 use force_constants_module
 implicit none
 integer l,t,ired,counter,rnk,nat,confg,jatom,katom,latom,res,cnt2, a !, b, c, d, e, f, i, reason
 integer al,be,ga,de,taui,tauj,tauk,taul,ni(3),nj(3),nk(3),nl(3),ti,frc_constr
 real(8), allocatable:: ared1d(:),aux(:)!, fc_ind(:)
 real(8) afrc(frc_constr,ngr),bfrc(frc_constr),rij !, amp, rijs
! real(8) junk
logical ex !new

 write(*,*) 'entering set_force_displacements routine'
 write(ulog,*)' SETUP_AMATRICES : ************************************'
! l is group index, t: general term index, and ti: independent term index
! this convention is followed in all setting constraint or amat subroutines
! set dimensions of ared matrix according to the number of constraints
!------------------------------------------------------------------
!       FORCE-DISPLACEMENT CONSTRAINTS
!------------------------------------------------------------------

! allocate(aforce(force_constraints,ngr),bforce(force_constraints))
 allocate( ared1d(ngr) )
 !allocate(rfcs_2(map(2)%ntotind))
 afrc = 0d0; bfrc = 0d0

!5 format(a,9(i5))
 write(ulog,6)'********************************************************'
 write(ulog,6)' Now making the force-displacement part of the A matrix'
 write(ulog,*)' Force_constraints     =', frc_constr
! now compute the force constraints on the FCs. The equation is:
! -F_ia = phi1_ia + phi2_{ia,jb} u_jb + 0.5 phi3_{ia,jb,kc} u_jb u_kc + ...
!   iatomterm(1,t)=(0,tau)   iatomterm(2,t)=(n1,tau1)
!   ixyzterm (1,t)=alpha     ixyzterm (2,t)=beta
! Really the (n1,tau1)=(iatomcell(:,iatomterm(2,t)),iatomcell0(iatomterm(2,t)))
  counter = 0
!  if ( include_fc(2) .eq. 2) then
!   write(*,*) "ALLOCATION CONDITION FC_IND"
 !  allocate(fc_ind(map(2)%ngr))
!  endif
  do confg=1,nconfigs
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
        rnk=1  !------------------------------------------
        if ( include_fc(rnk) .ne. 0 ) then

!          allocate(aux(ndindp(rnk))) ! size of aux =  # of independent FCs for that rank
           allocate(aux(map(rnk)%ntotind))  ! size of aux =  # of independent FCs for that rank
           cnt2= 0        ! cnt2 is the size of previous groups
           aux = 0
           do l=1,map(rnk)%ngr  ! sum over groups
              if(l.gt.1) cnt2 = cnt2 + map(rnk)%ntind(l-1)
              do t=1,map(rnk)%nt(l) ! sum over all therms in that group
                 if ( taui.eq. map(rnk)%gr(l)%iat(1,t) .and.  &
             &      al  .eq. map(rnk)%gr(l)%ixyz(1,t) ) then
                    do ti=1,map(rnk)%ntind(l)
! this is the index of the indep FC coming in the A*FC=b matrix product
                       ired = cnt2+ti
                       write(ulog,5)'l,cnt2,ti,ired=',l,cnt2,ti,ired
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

!write(*,*) 'Rank, res=',rnk,res ! write was on here so I commented this

        rnk=2  !------------------------------------------
      if ( include_fc(rnk) .ne. 0 ) then

! the force on atom i is from j in the supercell and all of its images for which the corresponding FC needs to be identified
!           allocate(aux(ndindp(rnk)))
           allocate(aux(map(rnk)%ntotind))  ! ntotind is the total number of independent terms
! for given (ni,taui,al) find its neighbors within nshells(2)
           cnt2= 0        ! cnt2 is the size of previous groups
           aux = 0
           groups: do l=1,map(rnk)%ngr  ! sum over groups
              if(l.gt.1) cnt2 = cnt2 + map(rnk)%ntind(l-1)
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
    ! K1 NEW --------------------------
  !              rij=length(atompos(:,nat)-atompos(:,jatom))
  !              if (rij.gt.radius(rnk)) cycle l8
! need to find how many images of j are connected to i with the same force constant.
    ! K1 NEW --------------------------
                 do ti=1,map(rnk)%ntind(l)
! this is the index of the indep FC coming in the A*FC=b matrix product
                    ired = cnt2+ti ! this is the corresponding index of aux
                    write(ulog,5)'l,cnt2,ti,ired=',l,cnt2,ti,ired
                    aux(ired) = aux(ired) + displ(be,jatom,confg)*map(rnk)%gr(l)%mat(t,ti)
                 enddo
 rij=atom_sc(nat)%equilibrium_pos .dot. atom_sc(nat)%equilibrium_pos
! if(l.ge.11 .and. l.le.14 .and. displ(be,jatom,confg).ne.0) then
!      write(*,9) nat,jatom,l,t, aux(11:14),rij
! if(l.ge.4 .and. l.le.7 .and. displ(be,jatom,confg).ne.0) then
!      write(*,9) nat,jatom,l,t, aux(4:7),rij
!endif
  !             write(ulog,22)'term,red_indx,ared=',t,iatomterm_2(1,t),  &
  ! &    ixyzterm_2(1,t),iatomterm_2(2,t),ixyzterm_2(2,t),ired,ared1d(ired)
              endif
           enddo l8
           enddo groups
 !          write(*,*) "DISPLACEMENT: ",displ(be,jatom,confg)
           if ( include_fc(rnk) .eq. 1 ) then
              ared1d(res+1:res+map(rnk)%ntotind) = aux(1:map(rnk)%ntotind)
!              ared1d(res+1:res+ndindp(rnk)) = aux(1:ndindp(rnk))
              res = res + map(rnk)%ntotind
!              res = res + ndindp(rnk)
              !if ( rnk .eq. 2) then
              !write(*,*) "AUX FOR INCLUDE_FC 1: ", aux
              !endif
          elseif ( include_fc(rnk) .eq. 2 ) then
!            write(*,*) "Check for include_fc(rnk) 2: "
!             inquire ( file="fc2_fit.dat", exist=ex)
!             if (ex) then
!               open(473,file="fc2_fit.dat")
!                  do i=1, map(2)%ngr
!                     read(473,*,iostat=reason) a, b, c, d, e, f, fc_ind(i), amp, rijs
!                     if (reason > 0) then
!                        write(*,*) "something is wrong"
!                     elseif (reason < 0) then
!                        write(*,*) "END OF FILE REACHED"
!                     else
!                        write(*,*) a,b,c,d,e,f,fc_ind(i),amp,rijs
!                     endif
!                  enddo
!               close(473)
!             endif
         !   write(*,*) "The value for rfcsind_2 call is: ", fc_ind
         !   allocate(map_2(map(rnk)%ngr))
         !    write(*,*) "The value for aux is: ", aux
  !          call read_fcs_2(2)
         !   write(*,*) "DISPLACEMENT: ",displ(be,jatom,confg)
!            write(*,*) "VALUE OF BFRC BEFORE: ", bfrc(counter)
            ! do l=1,map(rnk)%ngr !! HAVE TO UNCOMMENT THIS LINE ! IT WAS SOMETHING LIKE THIS WHY WE ARE DOING IT HERE? DOES IT NOT DO map(rnk)%ngr times the dot prod?
            !   t = map_2(l)
            !   write(*,*) "Check shape of aux and rfcsind_2: ", shape(aux), shape(rfcsind_2)
!               write(*,*) "Check for the value of BFRCS counter before: ", bfrc(counter)
!               if (counter .eq. 1) then
!               endif
   !            write(*,*) "BEFORE BFRC COUNTER: ", bfrc(counter)
            !   write(*,*) "VALUE OF BFRC BEFORE: ", bfrc(counter) 
               bfrc(counter)=bfrc(counter)-dot_product(aux,fc_ind) ! Have to uncomment this line
       !        write(*,*) "COUNTER AND BFRCS VALUE: ", counter, bfrc(counter)
            !  bfrc(counter)=bfrc(counter)-fcs_2(igroup_2(t))*aux(igroup_2(t))
! shouldn't this line simply be: bf = bf - fcs_2(j) * aux(j) ??
            ! enddo   !! HAVE TO UNCOMMENT THIS LINE
!             write(*,*) "VALUE OF FC_IND: ", fc_ind
!             write(*,*) "VALUE OF BFRC AFTER: ", bfrc(counter)
!             write(*,*) "VALUE OF AUX: ", aux
!             write(*,*) "VALUE OF DOT_PRODUCT: ", dot_product(aux,fc_ind)
           endif
      !     write(*,*) "The shape of bfrc(counter) is: ", size(bfrc)
           deallocate(aux)
        endif

!write(*,*) 'Rank, res=',rnk,res ! write was on here so I commented this


        rnk=3  !------------------------------------------
        if ( include_fc(rnk) .ne. 0 ) then

!           allocate(aux(ndindp(rnk)))
           allocate(aux(map(rnk)%ntotind))
! for given (ni,taui,al) find its two neighbors within nshells(3)
           cnt2= 0        ! cnt2 is the size of previous groups
           aux = 0
           do l=1,map(rnk)%ngr  ! sum over groups
              if(l.gt.1) cnt2 = cnt2 + map(rnk)%ntind(l-1)
           do t=1,map(rnk)%nt(l) ! sum over all therms in that group
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
    ! K1 NEW --------------------------
    !           rij=length(atompos(:,nat)-atompos(:,jatom))
    !           rik=length(atompos(:,nat)-atompos(:,katom))
    !           rjk=length(atompos(:,jatom)-atompos(:,katom))
    !   ! cycle if two of the pair distances is larger than radius
    !           if (rij.gt.radius(rnk)) then
    !              if (rik.gt.radius(rnk) .or. rjk.gt.radius(rnk) ) cycle
    !           else
    !              if (rik.gt.radius(rnk) .and. rjk.gt.radius(rnk) ) cycle
    !           endif
    ! K1 NEW --------------------------
                 do ti=1,map(rnk)%ntind(l)
! this is the index of the indep FC coming in the A*FC=b matrix product
                    ired = cnt2+ti ! this is the corresponding index of aux
                    write(ulog,5)'l,cnt2,ti,ired=',l,cnt2,ti,ired
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


!write(*,*) 'Rank, res=',rnk,res  ! There was Rank, res print out so I commented this section

        rnk=4  !------------------------------------------
        if ( include_fc(rnk) .ne. 0 ) then

!           allocate(aux(ndindp(rnk)))
           allocate(aux(map(rnk)%ntotind))
           cnt2= 0        ! cnt2 is the size of previous groups
           aux = 0
           do l=1,map(rnk)%ngr  ! sum over groups
              if(l.gt.1) cnt2 = cnt2 + map(rnk)%ntind(l-1)
           do t=1,map(rnk)%nt(l) ! sum over all therms in that group
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
    ! K1 NEW --------------------------
    !           rij=length(atompos(:,nat)-atompos(:,jatom))
    !           rik=length(atompos(:,nat)-atompos(:,katom))
    !           rjk=length(atompos(:,jatom)-atompos(:,katom))
    !   ! cycle if two of the pair distances is larger than radius
    !           if (rij.gt.radius(rnk)) then
    !              if (rik.gt.radius(rnk) .or. rjk.gt.radius(rnk) ) cycle
    !           else
    !              if (rik.gt.radius(rnk) .and. rjk.gt.radius(rnk) ) cycle
    !           endif
    ! K1 NEW --------------------------
                 do ti=1,map(rnk)%ntind(l)
! this is the index of the indep FC coming in the A*FC=b matrix product
                    ired = cnt2+ti ! this is the corresponding index of aux
                    write(ulog,5)'l,cnt2,ti,ired=',l,cnt2,ti,ired
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
5 format(a,9(1x,i6))
6 format(a,3(1x,i3),66(2x,f7.3))
!7 format(66(1x,g9.2))
!9 format(4(i4),66(1x,g9.2))
!12 format(a,2(1x,i3),66(2x,g11.4))
!22 format(a,6(1x,i3),66(2x,f9.5))

 end subroutine set_force_displacement_matrix
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
 character lineout*80
 integer j,i,n,m,mx,mxi,mxzero,rnk   &
 &      ,ntermszero,nd,ntindp,ierz,iert,ieri,ierg,jnk
 integer, allocatable:: iatomtermzero(:,:),ixyztermzero(:,:)
! logical new

 write(*,*) 'entering setup_maps routine'
!write(*,*) "The value of include_fc(2) isbikash: ", include_fc(2)
 call write_out(ulog,' #include_FCs ',include_fc)
 call write_out(ulog,' #   maxterms ',maxterms)
 call write_out(ulog,' #r=1 nshells ',nshells(1,:))
 call write_out(ulog,' #r=2 nshells ',nshells(2,:))
 call write_out(ulog,' #r=3 nshells ',nshells(3,:))
 call write_out(ulog,' #r=4 nshells ',nshells(4,:))
 write(ulog,*)'# trans,rot,huang,enforce_inv=',itrans,irot,ihuang,enforce_inv
 nd = sum(include_fc)
 if (nd.gt.4) then
    write(ulog,*)'SETUP_MAPS: nd>4, check your include_fc array ',nd
    write(*   ,*)'SETUP_MAPS: nd>4, check your include_fc array ',nd,include_fc(2)
    if (include_fc(2) .eq. 2) then
!      write(*,*) "SUM OF INCLUDE FORCE_CONSTANT IS GREATER THAN 4"
!   call read_harm_fcs(2)
    endif
!    stop
 endif

 write(ulog,*)' SETUP_MAPS: ********************************'
! loop over ranks that are supposed to be included
!if(include_fc(2) .eq. 2) then  ! This is what you have to do in order for to get force constant and declare fc_ind as global
!   write(*,*) "ENTERING SETUP MAPS READ_FCS_2-ROUTINE: ", map(2)%ntotind
!   call read_fcs_2(2)
!endif
!write(*,*) "THE FC_IND-BIKASH: ", fc_ind
! -----------------------------------
 rankloop: do rnk=1,4 !,1,-1
  if ( include_fc(rnk) .ne. 0 ) then
!write(*,*) "This condition is satisfied:::"
    mx=maxterms(rnk)
    mxi = maxtermsindep(rnk)
    write(ulog,4)'rank,mx,mxi,maxgrps=',rnk,mx,mxi,maxgroups
    mxzero = 90
    iert=1 ; ieri=1 ; ierz=1
!write(*,*) "The value for ierg is: ", ierg, iert, ieri, ierz
! start from a small maxterms and maxtermszero and increase by 100 if too small
  checkloop: do while (ierg+ierz+iert+ieri.ne.0)
!  write(*,*) "IF IERGIERZ is satisfied::"
!write(*,*) "VALUES: nterm, ntermsindep, mapmat, iatmtrm, ixyzterm, iatomtermindep, ixyztermindep, iatomtermzero, ixyztermzero", &
!& nterm, ntermsindep, mapmat, iatmtrm, iatomtermindep, ixyztermindep, iatomtermzero, ixyztermzero
     write(*,4) 'maxtrmzero,maxtrm,maxtrmindp,maxgrp=',mxzero,mx,mxi,maxgroups
!    allocate(mapmat(mx,mxi,maxgroups) &
!    &     ,iatmtrm(rnk,mx,maxgroups),ixyztrm(rnk,mx,maxgroups)  &
!    &     ,iatomtermindep(rnk,mx,maxgroups),ixyztermindep(rnk,mx,maxgroups)  )
     if(.not.allocated(nterm))       allocate(nterm(maxgroups))
     if(.not.allocated(ntermsindep)) allocate(ntermsindep(maxgroups))
     if(.not.allocated(mapmat))      allocate( mapmat(mx,mxi,maxgroups))
     if(.not.allocated(iatmtrm))     allocate( iatmtrm(rnk,mx,maxgroups))
     if(.not.allocated(ixyztrm))     allocate( ixyztrm(rnk,mx,maxgroups))
     if(.not.allocated(iatomtermindep)) allocate( iatomtermindep(rnk,mx,maxgroups))
     if(.not.allocated(ixyztermindep)) allocate( ixyztermindep(rnk,mx,maxgroups))
     if(.not.allocated(iatomtermzero)) allocate( iatomtermzero(rnk,mxzero))
     if(.not.allocated(ixyztermzero )) allocate( ixyztermzero(rnk,mxzero))

     write(ulog,4)' calling collect_force_constants for rank,ier:ztig=',rnk,ierz,iert,ieri,ierg
     call collect_force_constants(rnk,nshells(rnk,:),ngroups(rnk),ntermsindep,  &
&            iatomtermindep,ixyztermindep,nterm,iatmtrm,ixyztrm,   &
&            mapmat,ntermszero,iatomtermzero,ixyztermzero,rnk,mx,  &
&            mxi,mxzero,maxgroups,ierz,iert,ieri,ierg)
     write(ulog,4)' collect_force_constants called with ier:ztig=',ierz,iert,ieri,ierg
 4 format(a,9(i5))

     if (ierz.ne.0) then
   !     mxzero = mxzero+30*(4**rnk)
         mxzero = mxzero+333*rnk
         write(ulog,*)' mxtermszero   increased to ',mxzero
         if (allocated(iatomtermzero)) deallocate(iatomtermzero)
         if (allocated(ixyztermzero )) deallocate(ixyztermzero)
     endif
     if (iert.ne.0) then
!        mx = mx*2 !+(5**rnk)*natoms0
         mx = mx+30*rnk*natoms0
         write(ulog,*)' maxterms      increased to ',mx
         if (allocated(mapmat)) deallocate(mapmat)
         if (allocated(iatmtrm)) deallocate(iatmtrm)
         if (allocated(ixyztrm)) deallocate(ixyztrm)
         if (allocated(iatomtermindep)) deallocate(iatomtermindep)
         if (allocated(ixyztermindep)) deallocate(ixyztermindep)
     endif
     if (ieri.ne.0) then
!        mxi= mxi*2 ! +50*(2**rnk)*natoms0
         mxi= mxi+10*rnk*natoms0
         write(ulog,*)' maxtermsindep increased to ',mxi
         if (allocated(mapmat)) deallocate(mapmat)
     endif
     if (ierg.ne.0) then
!        maxgroups= maxgroups*2 !+100*rnk*natoms0
         maxgroups= maxgroups+10*rnk*natoms0
         write(ulog,*)' maxgroups     increased to ',maxgroups
         if (allocated(mapmat)) deallocate(mapmat)
         if (allocated(iatmtrm)) deallocate(iatmtrm)
         if (allocated(ixyztrm)) deallocate(ixyztrm)
         if (allocated(iatomtermindep)) deallocate(iatomtermindep)
         if (allocated(ixyztermindep)) deallocate(ixyztermindep)
         if (allocated(nterm)) deallocate(nterm)
         if (allocated(ntermsindep)) deallocate(ntermsindep)
     endif
      !if ( rnk .eq. 2) then
      !   if(include_fc(rnk) .eq. 2) then  ! This is what you have to do in order for to get force constant and declare fc_ind as global
      !      call read_fcs_2(2)
      !   endif
      !endif

  enddo checkloop

     write(umap,'(a,i1,a,10(1x,i2),a,i8,a,i5,a)')   &
&     '************ rank ',rnk,', shell ',nshells(rnk,:),', groups ',ngroups(rnk)
     if(ngroups(rnk).gt.0) then
       write(umap,'(a,i2,a,99(i5))')'rank=',rnk,' indepterms=',ntermsindep(1:ngroups(rnk))
       write(umap,'(a,i2,a,99(i5))')'rank=',rnk,' all  terms=',nterm(1:ngroups(rnk))
       write(ulog,*)'Rank, size of allocated array, real size=',rnk,mx,nterm(ngroups(rnk))

        allocate(map(rnk)%gr   (ngroups(rnk)),    &
        &        map(rnk)%nt   (ngroups(rnk)),    &
        &        map(rnk)%ntind(ngroups(rnk)) )
        map(rnk)%ngr      = ngroups(rnk)
        map(rnk)%nt   (:) = nterm      (1:ngroups(rnk))
        map(rnk)%ntind(:) = ntermsindep(1:ngroups(rnk))
        jnk=sum(map(rnk)%ntind(:))
        allocate(map(rnk)%err(jnk))
        write(ulog,*) 'ngroups(rnk)=',ngroups(rnk)
        do i=1,ngroups(rnk)
           write(ulog,*) 'i, nterm(i),ntindp(i)=',i,nterm(i),ntermsindep(i)
           allocate(map(rnk)%gr(i)%mat (nterm(i),ntermsindep(i)),   &
       &         map(rnk)%gr(i)%iat (rnk,nterm(i)),                 &
       &         map(rnk)%gr(i)%ixyz(rnk,nterm(i)),                 &
       &         map(rnk)%gr(i)%iatind (rnk,ntermsindep(i))  ,      &
       &         map(rnk)%gr(i)%ixyzind(rnk,ntermsindep(i))  )
        enddo
        write(ulog,*)' allocation of MAP done! '
        map(rnk)%err = 'C'

! copy from temp array into real arrays with actual size
        do i=1,map(rnk)%ngr    !  = ngroups(rnk)
          map(rnk)%gr(i)%mat (:,:) = mapmat(1:nterm(i),1:ntermsindep(i),i)
          map(rnk)%gr(i)%iat (:,:)   = iatmtrm(1:rnk,1:nterm(i),i)
          map(rnk)%gr(i)%ixyz(:,:)   = ixyztrm (1:rnk,1:nterm(i),i)
          map(rnk)%gr(i)%iatind (:,:)= iatomtermindep(1:rnk,1:ntermsindep(i),i)
          map(rnk)%gr(i)%ixyzind(:,:)= ixyztermindep (1:rnk,1:ntermsindep(i),i)
        enddo

        m=0
        do i=0,nshells(rnk,1)
           m = m + atom0(1)%shells(i)%no_of_neighbors
        enddo
        write(*,*)'value of m, natoms0 and inv_constraints is: ',m,natoms0, inv_constraints
        write(ulog,*)'atm 1 has ',m,'inclusve nghbrs within ',nshells(rnk,1),' shell'
        inv_constraints = inv_constraints + natoms0*12*m
        write(ulog,*)' rank=',rnk,' groups=',ngroups(rnk)
        write(ulog,'(a,99(i4))')' nterm=',nterm(1:ngroups(rnk))
        write(ulog,*)' Cumulative invce_cnstrnts for this shell ', inv_constraints

     endif

     if ( rnk .eq. 2) then
         if( include_fc(rnk) .eq. 2) then
   !         write(*,*) "INFO-NGR: ", map(rnk)%ngr
            call read_fcs_2(rnk)
         endif
     endif

     deallocate(iatmtrm,ixyztrm,mapmat,  &
  &     iatomtermzero,ixyztermzero,iatomtermindep,ixyztermindep)

     write(ulog,*)'iattrm,ixyztrm,mapmat, iattrm0,ixyztrm0,iattrmindp,ixyztrmindp deallocated'

  else
     map(rnk)%ngr=0 ! So with this it goes to map(rnk)%ngr in this line and assign it to zero
     write(ulog,*) 'Allocating map(rank=',rnk,')'
     allocate(map(rnk)%nt(1),map(rnk)%ntind(1) )
     map(rnk)%nt(1)=0
     map(rnk)%ntind(1)=0
  endif

 enddo rankloop

!if(include_fc(2) .eq. 2) then  ! This is what you have to do in order for to get force constant and declare fc_ind as global
!   call read_fcs_2(2) ! You need to put this here
!   write(*,*) "Check FC_IND-Bikash: ",fc_ind
!endif

 write(*,*) 'SETUP_MAPS: Exited the main rank loop , and going to write by calling ustring'

 do rnk=1,4
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
            write(umap,'(4x,a)') lineout(1:m)
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

 ngr = 0
 do rnk=1,4
!    if ( include_fc(rnk) .ne. 0 ) then
! get the total number of independent and full terms for each rank
!       call get_dim(map(rnk),map(rnk)%ntotind,map(rnk)%ntot)
!       write(ulog,*)'SETUP_MAPS: RANK, NDINDEP, NDFULL=',rnk,map(rnk)%ntotind,map(rnk)%ntot
!       ngr = ngr + map(rnk)%ntotind
!    endif

    if ( include_fc(rnk) .ne. 0 ) then
      ! get the total number of independent and full terms for each rank
             call get_dim(map(rnk),map(rnk)%ntotind,map(rnk)%ntot)
             write(ulog,*)'SETUP_MAPS: RANK, NDINDEP, NDFULL=',rnk,map(rnk)%ntotind,map(rnk)%ntot
             if ( include_fc(rnk) .eq. 1 ) ngr = ngr + map(rnk)%ntotind
    endif

 enddo
 write(ulog,*)'END OF SETUP_MAPS: NGR=',ngr
 write(*,*) 'exiting setup_maps routine'

 end subroutine setup_maps
!===============================================================================
 subroutine include_constraints
!! now that we have the 3 matrices atransl,arot,ahuang and aforce, here we decide
!! how to put them together to do the SVD
!! outputs are amat and bmat and their dimensions:dim_al,dim_ac, which go to SVD
 use svd_stuff
 use params
 use ios
 implicit none
 integer n_constr
! real(8), allocatable:: w(:),v(:,:)
! real(8), allocatable:: a11(:,:),b1(:),a12(:,:),a11i(:,:),b2(:)
!write(*,*) "THE FORCE_CONSTRAINTS VALUE IS: ", force_constraints
!-----------------------------------------------------------------------
! in this case, invariances are obtained in a rms sense and are not exact.
   dim_al= force_constraints
   dim_ac = ngr
! we put the invariance constraints into A
   write(ulog,*)'INCLUDE_CONSTRAINTS :c(forces),dim(A)=',force_constraints ,dim_al
   write(ulog,*)'INCLUDE_CONSTRAINTS :c(transl),dim(A)=',transl_constraints,dim_al
   write(ulog,*)'INCLUDE_CONSTRAINTS :c(rotatn),dim(A)=',rot_constraints   ,dim_al
   write(ulog,*)'INCLUDE_CONSTRAINTS :c(huang ),dim(A)=',huang_constraints ,dim_al
!   write(*,*) "The dimension of the AMATRIX, dima_al, dim_ac is: ", dim_al, ngr
   if (itrans .ne. 0) dim_al= dim_al+ transl_constraints
   if (irot   .ne. 0) dim_al= dim_al+    rot_constraints
   if (ihuang .ne. 0) dim_al= dim_al+  huang_constraints
   !if ( include_fc(2) .ne. 2) then
   allocate(amat(dim_al,dim_ac),bmat(dim_al))
   !endif
   !if (include_fc(2) .eq. 2) then
   !   allocate(amat(dim_al,map(2)%ngr:dim_ac),bmat(map(2)%ngr:dim_al))
   !endif
   amat = 0d0 ; bmat = 0d0
   write(ulog,*)' size(a) nlines,ncolumn=(dim(a1d))=',dim_al,dim_ac
   n_constr = 1
! BEWARE these change if some of the FCS are read in and not fitted, bmat has to change
   if (itrans .ne. 0) then
    !  if (include_fc(2) .ne. 2) then
      !write(*,*) "Size for AMAT is: n_constr, n_constr-1+transl_constraints, ngr ", n_constr, n_constr-1+transl_constraints, ngr
      amat(n_constr:n_constr-1+transl_constraints,1:ngr) =  &
&                    atransl(1:transl_constraints,1:ngr)
      bmat(n_constr:n_constr-1+transl_constraints)       = 0d0
      n_constr = n_constr + transl_constraints
     ! endif
      !if (include_fc(2) .eq. 2) then
         !write(*,*) "Size for AMAT is: n_constr, n_constr-1+transl_constraints, ngr ", n_constr, n_constr-1+transl_constraints, ngr
       !  amat(n_constr:n_constr-1+transl_constraints,map(2)%ngr:ngr) =  &
   !&                    atransl(1:transl_constraints,map(2)%ngr:ngr)
    !     bmat(map(2)%ngr:map(2)%ngr-1+transl_constraints)       = 0d0
     !    n_constr = n_constr + transl_constraints
      !endif
   endif
   if (irot .ne. 0) then
    !  if ( include_fc(2) .ne. 2) then
      amat(n_constr:n_constr-1+   rot_constraints,1:ngr) =  &
&                          arot(1:rot_constraints,1:ngr)
      bmat(n_constr:n_constr-1+rot_constraints) = brot(1:rot_constraints)
      n_constr = n_constr + rot_constraints
     ! endif
      !if ( include_fc(2) .eq. 2) then
     !    amat(n_constr:n_constr-1+   rot_constraints,map(2)%ngr:ngr) =  &
! &                          arot(1:rot_constraints,map(2)%ngr:ngr)
  !    bmat(map(2)%ngr:map(2)%ngr-1+rot_constraints) = brot(map(2)%ngr:map(2)%ngr+rot_constraints)
   !   n_constr = n_constr + rot_constraints
    !  endif
   endif
   if (ihuang .ne. 0) then
   !   if ( include_fc(2) .ne. 2) then
      amat(n_constr:n_constr-1+ huang_constraints,1:ngr) =  &
&                          ahuang(1:huang_constraints,1:ngr)
      bmat(n_constr:n_constr-1+huang_constraints)        = 0d0
      n_constr = n_constr + huang_constraints
    !  endif

     ! if ( include_fc(2) .eq. 2) then
   !      amat(n_constr:n_constr-1+ huang_constraints,map(2)%ngr:ngr) =  &
   ! &                          ahuang(1:huang_constraints,map(2)%ngr:ngr)
   !      bmat(map(2)%ngr:map(2)%ngr-1+huang_constraints)        = 0d0
   !      n_constr = n_constr + huang_constraints
   !      endif
   endif
!   if (include_fc(2) .ne. 2) then
   amat(n_constr:n_constr-1+force_constraints,1:ngr) =  &
 &                  aforce(1:force_constraints,1:ngr)
   bmat(n_constr:n_constr-1+force_constraints) = bforce(1:force_constraints)
  ! endif
 !  write(*,*) "Value for force_constraints is: ", force_constraints
  ! if (include_fc(2) .eq. 2) then
  !    amat(n_constr:n_constr-1+force_constraints,map(2)%ngr:ngr) =  &
!&                  aforce(1:force_constraints,map(2)%ngr:ngr)
 !  bmat(map(2)%ngr:force_constraints) = bforce(map(2)%ngr:force_constraints)
  ! endif

 deallocate (aforce,bforce)

 write(ulog,*)' aforce,bforce written into amat,bmat, and deallocated'

 end subroutine include_constraints
!===============================================================================
 subroutine homogeneous_constraints_overlap(n_constr)
! now that we have the 3 matrices atransl,arot,ahuang,
! we put them together to do the SVD
! outputs are amat and bmat and their dimensions:dim_al,dim_ac, which go to SVD
 use svd_stuff
 use params
 use ios
 implicit none
 integer i,j,n_constr

!-----------------------------------------------------------------------
! in this case, invariances are obtained in a rms sense and are not exact.
   dim_hom= 0
   dim_ac = ngr
! we put the invariance constraints into A
   if (itrans .ne. 0) dim_hom= dim_hom+ transl_constraints
   if (irot   .ne. 0) dim_hom= dim_hom+    rot_constraints
   if (ihuang .ne. 0) dim_hom= dim_hom+  huang_constraints
   allocate(ahom(dim_hom,dim_ac))
   ahom = 0d0
   write(ulog,*)' size(ahom) nlines,ncolumn=(dim(a1d))=',dim_hom,dim_ac
   n_constr = 1
   if (itrans .ne. 0) then
      ahom(n_constr:n_constr-1+transl_constraints,1:ngr) =  &
&                    atransl(1:transl_constraints,1:ngr)
      n_constr = n_constr + transl_constraints
   endif
   if (irot .ne. 0) then
      ahom(n_constr:n_constr-1+   rot_constraints,1:ngr) =  &
&                          arot(1:rot_constraints,1:ngr)
      n_constr = n_constr + rot_constraints
   endif
   if (ihuang .ne. 0) then
      ahom(n_constr:n_constr-1+ huang_constraints,1:ngr) =  &
&                          ahuang(1:huang_constraints,1:ngr)
      n_constr = n_constr + huang_constraints
   endif
   n_constr=n_constr-1
   write(ulog,*)' ahom is setup with ',n_constr,' lines'
   if (dim_hom.ne.n_constr) then
      write(ulog,*)' dim_hom, n_constr are not equal! ', dim_hom,n_constr
      stop ! The program goes upto here and stops. This line is 
   endif

! now form the overlap matrix from the constraints
   allocate(overl(n_constr,n_constr))
   do i=1,n_constr
   do j=1,n_constr
      overl(i,j)=dot_product(ahom(i,:),ahom(j,:))
   enddo
   enddo

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
 use atoms_force_constants
 use ios
 use geometry
 use params
 use svd_stuff
 implicit none
 integer i,m,rnk

 do rnk=1,4
! count the neighbors within a shell to get an estimate of invariance constraints
   m=0
   do i=0,nshells(rnk,1)
      m = m + atom0(1)%shells(i)%no_of_neighbors
   enddo
   write(ulog,*)'atm 1 has ',m,'inclusve nghbrs within ',nshells(rnk,1),' shell'
   inv_constraints = inv_constraints + natoms0*6*m
   write(ulog,*)' Cumulative invce_cnstrnts for this shell ', inv_constraints
 enddo

 end subroutine estimate_inv_constraints
! ============================================
 subroutine get_dim(mapd,ndimindep,ndimfull)
!! calculate the number of in/dependent fcs for all groups to fillup map%ntot and %ntotind
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
!! check if there is any column that is zero (the corresponding FC can not be extracted)
 use ios
 use params
 use svd_stuff
 implicit none
 integer line,col,i,j,noforc,nosym,rnk,dum
 real(8) afrc(line,col)

 write(ulog,*)'==== CHECK_ZERO_COLUMN: lines and columns are:',line,col
 main: do i=1,col
! separate the contribution of symmetries from the force-displacements
    nosym = 0
    col_loop: do j=1,line-force_constraints
write(*,*)'j=',j,line-force_constraints
       if (abs(afrc(j,i)) .gt. 1d-8) then
          nosym = 1+nosym
          exit col_loop
       endif
    enddo col_loop
! now the contribution of force-displacements
    noforc = 0
    col_loop2: do j=1+line-force_constraints,line
       if (abs(afrc(j,i)) .gt. 1d-8) then
          noforc = 1+noforc
          exit col_loop2
       endif
    enddo col_loop2

!   if (nosym .eq. 0) then  ! column i was all = zero in the symmetry constraints
!      write(6,*)' The FC #',i,' is not involved in the symmetry constraints'
!      write(ulog,*)' The FC #',i,' is not involved in the symmetry constraints'
!   endif
    if (noforc .eq. 0) then  ! column i was all = zero in the symmetry constraints
       write(6,*)' The FC #',i,' is not involved in the force-displacements data'
       write(6,*)' Thus it might only be fixed from the enforced symmetry constraints'
       write(ulog,*)' The FC #',i,' is not involved in the force-displacements data'
       write(ulog,*)' Thus it might only be fixed from the enforced symmetry constraints'
    endif
    if ((noforc .eq. 0) .and. (nosym .eq. 0)) then  ! column i was all = zero in the symmetry constraints
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
       write(ulog,*)' WARNINIG: column no. ',i,' of Amatrix is 0, the corresponding '// &
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
!! for a given t in {1,...,ngr} find its corresponding rank and ntindp number nti
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
 write(ulog,*)'FIND_T: ntindep(rnk)=',m

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
!use atoms_force_constants
 use geometry
 use params
 use ios
 implicit none
 integer j,t,counter,m,n,j0,i
 logical new
 real(8) junk,zero
 real(8), dimension(n) :: a1d,b
 real(8), dimension(m,n) :: a

 zero = 0d0
 new=.true.
 j0=0

    checkloop: do t=1,counter
! check to see if they are multiples of each other
!   write(ulog,*)' comparing with line ',t,' in ared matrix'
          do j=1,n  ! = no of columns in ared; it's really ngr
             if ( a1d(j) .myeqz. zero ) cycle
             j0 = j   ! it is the first nonzero index
!    write(ulog,*)' found j0=',j0
             exit
          enddo
          if (j.gt.n .and. verbose ) then
             write(ulog,*)' this line of ared1d is all=0'
             write(ulog,4)(a1d(i),i=1,n)
!            stop
          endif
          junk = a(t,j0)/a1d(j0)
          b = a1d*junk   ! this is to avoid repeating amat and -amat
          if ( a(t,:) .myeqz. b(:) ) then
             new=.false.
             exit checkloop
          endif
    enddo checkloop
    if(verbose) then
    if (new) then
       write(ulog,5)'NEW ',counter+1,a1d
    else
       write(ulog,5)'COMPARE2PREVIOUS: LINE EXISTED ',counter+1,a1d
    endif
    endif

4 format(99(1x,g11.4))
5 format(a,i5,99(1x,f7.3))
 end subroutine compare2previous_lines
!======================================================================================
 subroutine set_huang_inv_constraints
 use params
 use atoms_force_constants
 use svd_stuff
 use ios
! use geometry
 use force_constants_module
 implicit none
 integer i,j,al,be,ga,de,t,ti,mterm,cnt2,g,rnk,voigt,vab,vcd,counter,ired,res  !,nterm
 real(8) huang,rr(3)
 real(8), allocatable:: ared1d(:),zero(:)

 if (ihuang .ne. 1) return

   rnk = 2; huang_constraints=15
!  res = sum(map(1)%ntind(:))
   if(map(1)%ngr.gt.0) then
     res = sum(map(1)%ntind(:))
   else
     res=0
   endif
   write(ulog,*)' SETTING HUANG INVARIANCE CONSTRAINTS *****************'
   dim_al = huang_constraints
!   ngr = (sum(ngroups(:)))
   allocate( ahuang(dim_al,ngr),ared1d(ngr),zero(ngr) )
   counter = 0; zero = 0d0 ; ahuang = 0d0
   ared1d = zero

   do al=1,3
   do be=al,3
   do ga=1,3
   do de=ga,3
      vab = voigt(al,be)
      vcd = voigt(ga,de)
      if ( vab.le.vcd ) cycle
      huang = 0
      mterm = 0
      atomloop: do i=1,natoms0
      cnt2=0     ! cnt2 counts the previous independent terms
      do g=1,map(2)%ngr
         if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
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
      enddo
    enddo atomloop

      counter = counter+1
      if(verbose) write(ulog,6)'counter,mterm,vab,vcd=',counter,mterm,al,be,ga,de,vab,vcd
      if (counter.eq.dim_al+1) then
         write(ulog,*)' DIM_AL TOO SMALL, huang_inv=15, counter= ',counter
         stop
      endif
      ahuang(counter,:) = ared1d(:)
!     write(umatrx,7)(ahuang(counter,dum),dum=1,ngr)

   enddo
   enddo
   enddo
   enddo

   deallocate( ared1d,zero )

6 format(a,8(1x,i4),9(2x,f7.3))
!7 format(66(1x,g9.2))

 end subroutine set_huang_inv_constraints
!===============================================================================
 subroutine check_huang
! this subroutine checks Huang invariance relations:
! sum_r ( phi_0,r^al,be r^ga r^de -  phi_0,r^ga,de r^al r^be ) = 0
 use params
 use atoms_force_constants
 use svd_stuff
 use ios
! use geometry
 use force_constants_module
 implicit none
 integer i,j,al,be,ga,de,t,ti,cnt2,mterm,rnk,voigt,g,vab,vcd,res
 real(8) huang,rr(3)

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
    huang = 0
    mterm = 0
    do i=1,natoms0
      cnt2=0
      do g=1,map(2)%ngr
         if(g.gt.1) cnt2 = cnt2 + map(rnk)%ntind(g-1)
      do t=1,map(2)%nt(g)
         if ( map(2)%gr(g)%iat(1,t).ne.i .or. map(2)%gr(g)%ixyz(1,t).ne.al  &
       & .or. map(2)%gr(g)%ixyz(2,t).ne.be ) cycle
         mterm = mterm+1
         j  =  map(2)%gr(g)%iat(2,t)
!         tauj=iatomcell0(j)
!         nj=iatomcell(:,j)
         rr = atompos(:,j)-atompos(:,i)
         do ti =1,map(2)%ntind(g)
            huang = huang + rr(ga)*rr(de)*map(2)%gr(g)%mat(t,ti)*fcs(res+cnt2+ti)
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
            huang = huang - rr(al)*rr(be)*map(2)%gr(g)%mat(t,ti)*fcs(res+cnt2+ti)
         enddo
      enddo
      enddo
    enddo
    write(ulog,5)i,al,be,ga,de,vab,vcd,mterm,huang

   enddo
   enddo
   enddo
   enddo
 write(ulog,*)'HUANG INVARIANCE RELATIONS CHECKED        ###############################'
5 format(8(1x,i4),9(1x,g13.6))

 end subroutine check_huang
! ============================================
 subroutine truncate(n,ws,n_kept)
!! find a gap in the ws and drop the small ws values
 implicit none
 integer n,n_kept,i,i1
 real(8) ws(n), logws(n),gap(n-1),maxgap1,wmax

 wmax=maxval(ws)
 ws=ws/wmax
 logws=log(ws)/log(10.)  ! decimal log of ws
 if (logws(1) .ne.0) then
    write(*,*)'TRUNCATE: ws(1),logws(1)=',ws(1),logws(1)
    stop
 endif
 write(*,*)'TRUNCATE: dimension of ws is n=',n

! calculate the level separations
 do i=1,n-1
    gap(i)=logws(i)-logws(i+1)
    write(*,5)'TRUNCATE: i, ws,logws, gap(i)',i,ws(i),logws(i),gap(i)
 enddo
5 format(a,i4,9(1x,g11.4))

 do i=1,n
    if (logws(i).gt. -5) cycle
    exit
 enddo
 n_kept=i  ! any logws above -5, i.e. ws>10^-5 is kept

 maxgap1=maxval(gap(i+1:n-1))
 i1=maxloc(gap(i+1:n-1),n-i-1)
 n_kept=i1

 ws=ws*wmax

! if(i1.lt.n_kept) then
! the largest gap is before n_kept , ie larger than 10^-6
! use the second largest gap
!   i2=maxloc(gap(i1+1:n-1))
!else
! the largest gap comes after, and so that's where we will truncate
!   n_kept=i1
!   return
! endif

 end subroutine truncate