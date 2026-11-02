  module svd_stuff
 use constants, only : r15
! stuff related to Harold's subroutines, mappings and the svd matrices

 type groupmatrix
    real(r15), allocatable:: mat(:,:)
    integer, allocatable:: iat(:,:),ixyz(:,:),iatind(:,:),ixyzind(:,:)
 end type groupmatrix
 type fulldmatrix
    real(r15), allocatable:: amatrx(:,:),transinv(:,:),rotinv(:,:)
    integer alines,tranlines,rotlines
    type(groupmatrix), allocatable :: gr(:)  ! groups of fcs
    integer, allocatable:: nt(:),ntind(:),keep(:)  ! full and independent fc terms for each group
    character(1), allocatable:: err(:)
    integer ngr,ntotind,ntot,nkeptind  ! number of groups, total number of independent terms and full terms
 end type fulldmatrix

 integer maxrank,itrans,irot,ihuang,enforce_inv
 integer, allocatable :: nlines(:),tind2(:)
 parameter(maxrank=8)
 integer nterms(maxrank),ngroups(maxrank),ngrnk(maxrank)
 integer, allocatable :: nterm(:),ntermsindep(:)
 integer, allocatable :: iatmtermindp(:,:,:),ixyztermindp(:,:,:)
 integer, allocatable :: iatmtrm(:,:,:),ixyztrm(:,:,:)
 real(r15), allocatable :: mapmat(:,:,:)
 real(r15) svdcut 
 integer, allocatable:: iatomterm_1(:,:),ixyzterm_1(:,:),igroup_1(:),map_1(:)
 integer, allocatable:: iatomterm_2(:,:),ixyzterm_2(:,:),igroup_2(:),map_2(:)
 integer, allocatable:: iatomterm_3(:,:),ixyzterm_3(:,:),igroup_3(:),map_3(:)
 integer, allocatable:: iatomterm_4(:,:),ixyzterm_4(:,:),igroup_4(:),map_4(:)
 character(1), allocatable:: err_1(:),err_2(:),err_3(:),err_4(:)
 type(fulldmatrix) map(maxrank)
 real(r15), allocatable:: fcrnk(:,:)
 real(r15), allocatable:: ampterm_1(:),fcs_1(:)
 real(r15), allocatable:: ampterm_2(:),fcs_2(:),grun_fc(:)
 real(r15), allocatable:: ampterm_3(:),fcs_3(:)
 real(r15), allocatable:: ampterm_4(:),fcs_4(:)
 real(r15), allocatable:: amat(:,:),bmat(:),sigma(:),fcs(:),ahom(:,:),bhom(:),kernelbasis(:,:)
 real(r15), allocatable:: a11ia12(:,:),fc1(:)
 real(r15), allocatable:: atransl(:,:),btransl(:),arot(:,:),brot(:),ahuang(:,:),bhuang(:), &
 &                      aforce(:,:),bforce(:),energies(:) !,afrc(:,:),bfrc(:),engy(:)
 real(r15), allocatable:: fc_ind(:) ! This is for getting force constant by reading FC2.dat file
 integer inv_constraints,force_constraints,dim_al,dim_ac,n_indep,newdim_al,dim_hom   &
&        ,transl_constraints, rot_constraints, huang_constraints,nindepfc
contains

 function counteri(rnk,g,ti) result(cnt)
!! counteri(rnk,g,ti) counts the number of independent terms of rank=rnk from all groups  
!! up to group g-1 and the extra ti terms in group g
 implicit none
 integer, intent(in) :: rnk,g,ti
 integer cnt

 if(g.eq.1) then
    cnt=ti
 else
    cnt=sum(map(rnk)%ntind(1:g-1))+ti
 endif

 end function counteri
!------------------------------------------------------
 function current(rnk,g,ti) result(cnt)
!! current(rnk,g,ti) is the location of the the KEPT independent term ti in group g, of rank
!! rnk in the fc list; it is smaller or equal to counteri
 implicit none
 integer, intent(in) :: rnk,g,ti
 integer cnt,c2

 c2=counteri(rnk,g,ti)
 if(rnk.eq.2) then
   if( c2 .gt. size(map(rnk)%keep(:)) ) then
     write(*,*)'g,ti,counter2,size=',g,ti,c2,size(map(rnk)%keep(:))
     stop
   else
     cnt = sum(map(rnk)%keep(1:c2))
   endif
 else
   cnt = c2  ! everything is kept
 endif

 end function current
!------------------------------------------------------
! subroutine set_maxterms
!   maxterms(1)=40 !100
!   maxterms(2)=4000 !500
!   maxterms(3)=5000 !1800
!   maxterms(4)=3000 !2000
!   maxtermzero(1)=20 !500
!   maxtermzero(2)=1000 !2000
!   maxtermzero(3)=5000!5000
!   maxtermzero(4)=3000 !8000
!   maxtermsindep(1)=10
!   maxtermsindep(2)=150
!   maxtermsindep(3)=150
!   maxtermsindep(4)=300
!   maxgroups(1)=10
!   maxgroups(2)=100
!   maxgroups(3)=150
!   maxgroups(4)=300
! end subroutine set_maxterms

!--------------------------------------------
 subroutine remove_zeros_all(nl,nc)
! removes lines which are all zero from the amatrix
 use ios
 use params
! use svd_stuff  !for itrans,irot,ihuang
 implicit none
 integer, intent(in) :: nc
 integer, intent(inout) :: nl
! real(r15), dimension(:,:), intent(inout) :: amatr
! real(r15), dimension(:)  , intent(inout) :: bmatr
! real(r15), intent(inout) :: amatr(nl,nc)
! real(r15), intent(inout) :: bmatr(nl)
 integer i,j,nz,tra,rot,hua,tra_out,rot_out,hua_out
 real(r15), allocatable:: aa(:,:),bb(:) !at(:,:),bt(:),ar(:,:),br(:),ah(:,:),bh(:),
 real(r15) small,suml
 integer, allocatable:: zz(:)

 small=1d-6 !10  displacements or cos(t) are never that small

 tra=0; rot=0; hua=0
 if (itrans .ne. 0) tra_out= transl_constraints
 if (irot   .ne. 0) rot_out= rot_constraints
 if (ihuang .ne. 0) hua_out= huang_constraints

! zero lines have already been removed in the corresponding subroutines!
! call remove_zeros(tra,nc,amat(1:tra,1:nc),bmat(1:tra),  &
! &  tra_out,at(1:tra_out,1:nc),bt(1:tra_out))
! call remove_zeros(rot,nc,amat(tra+1:tra+rot,1:nc),bmat(tra+1:tra+rot),  &
! &  rot_out,ar(1:rot_out,1:nc),br(1:rot_out))
! call remove_zeros(hua,nc,amat(tra+rot+1:tra+rot+hua,1:nc),bmat(tra+rot+1:tra+rot+hua), &
! &  hua_out,ah(1:hua_out,1:nc),bh(1:hua_out))

 nl=tra_out+rot_out+hua_out+force_constraints
! allocate(aa(nl,nc),bb(nl))
! aa(1:tra_out,1:nc)=at
! bb(1:tra_out   )=bt
! aa(tra_out+1:tra_out+rot_out,1:nc)=ar
! bb(tra_out+1:tra_out+rot_out   )=br
! aa(tra_out+rot_out+1:tra_out+rot_out+hua_out,1:nc)=ah
! bb(tra_out+rot_out+1:tra_out+rot_out+hua_out   )=bh
! aa(tra_out+rot_out+hua_out+1:nl,1:nc)=aforce
! bb(tra_out+1:tra_out   )=bt

! call remove_zeros(rot,nc,amat(tra+1:tra+rot,1:nc),bmat(tra+1:tra+rot), &
! & rot_out,aa(tra_out+1:tra_out+rot_out,1:nc),bb(tra_out+1:tra_out+rot_out))
! call remove_zeros(hua,nc,amat(tra+rot+1:tra+rot+hua,1:nc),bmat(tra+tro+1:tra+rot+hua), &
! &  hua_out,aa(tra_out+rot_out+1:tra_out+rot_out+hua_out,1:nc), &
! & bb(tra_out+rot_out+1:tra_out+rot_out+hua_out))

 allocate(zz(nl))
 zz = 1
 nz = 0
 do i=1,nl
! get the norm1 of each line
!   suml = 0
!   do j=1,nc
!      if ( abs(amatr(i,j)) .ge. suml ) suml = abs(amatr(i,j))
!   enddo
    suml = maxval(abs(amat(i,:)))     ! sqrt(sum(amatr(i,:)*amatr(i,:))/nc)
    if (suml.lt.small) then
       if (abs(bmat(i)) .lt. small) then
          if(verbose) write(ulog,*)' line #',i,' of Amatrix and Bmatrix is zero; it will be removed!'
       else
          write(ulog,2)' Inconsistency in line i of Amatrix which is zero!',i,suml
          write(ulog,3)' It will be removed; Corresponding line in bmatr is=',bmat(i)
          write(ulog,*)' Need to increase your range of FCs'
       endif
! record the line index for later removal
       zz(i) = 0
       nz = nz + 1  ! nz counts the number of zero lines
    endif
 enddo
 write(ulog,*)' number of lines=0 in amat are=',nz
 allocate(aa(nl-nz,nc),bb(nl-nz))

 deallocate(amat,bmat)
 allocate(amat(nl-nz,nc),bmat(nl-nz))
 amat = aa
 bmat = bb
 write(ulog,*)'REMOVE ZEROS: new number of lines in amat=',nl
 deallocate(aa,bb,zz)

 return

 j = 0
 do i=1,nl
    if (zz(i).ne.0) then
       j = j+1
       aa(j,:) = amat(i,:)
       bb(j) = bmat(i)
    elseif(i.le.transl_constraints) then
       tra_out=tra_out - 1 !nsl_constraints=transl_constraints-1
    elseif(i.le.transl_constraints+rot_constraints) then
       rot_out=rot_out - 1 !_constraints=rot_constraints-1
    elseif(i.le.transl_constraints+rot_constraints+huang_constraints) then
       hua_out = hua_out-1 !ng_constraints=huang_constraints-1
    else
       force_constraints=force_constraints-1
    endif
 enddo
 if (sum(zz).ne.nl-nz) then
    write(ulog,*)'REMOVE ZEROS: Inconsistency! sum(zz).ne. nl.nz ',sum(zz),nl-nz
    stop
 endif
 if (j.ne.nl-nz) then
    write(ulog,*)'REMOVE ZEROS: Inconsistency! j.ne. nl.nz ',j,nl-nz
    stop
 endif

 inv_constraints = tra_out+rot_out+hua_out !transl_constraints+rot_constraints+huang_constraints

 write(ulog,4)'After remove_zeros: NEW tra,rot,hua,force_constaints=',tra_out,  &
 &       rot_out,hua_out,force_constraints

 nl=j
 deallocate(amat,bmat)
 allocate(amat(nl,nc),bmat(nl))
 amat = aa
 bmat = bb
 write(ulog,*)'REMOVE ZEROS: new number of lines in amat=',nl
 deallocate(aa,bb,zz)

 2 format(a,i9,g11.4)
 3 format(a,g11.4)
 4 format(a,4(i9))

 end subroutine remove_zeros_all

!--------------------------------------------
 subroutine remove_zeros(nl,nc,matin,vectin,nlout,matout,vectout)
!! removes lines which are all zero from the matrix mat and the vector v
 use ios , only : ulog
 use params, only : verbose
 implicit none
 integer, intent(out) :: nlout
 integer, intent(in) :: nl,nc
 real(r15), dimension(nl,nc), intent(in) :: matin
 real(r15), dimension(nl)   , intent(in) :: vectin
 real(r15), allocatable, dimension(:,:), intent(out) :: matout
 real(r15), allocatable, dimension(:)  , intent(out) :: vectout
 !real(r15), dimension(nl,nc), intent(out) :: matout
 !real(r15), dimension(nl   ), intent(out) :: vectout
 integer i,j,nz
 real(r15) small,suml
 integer zz(nl)

 small=1d-6 !10  displacements or cos(t) are never that small

 zz = 1
 nz = 0
 do i=1,nl
    suml = maxval(abs(matin(i,:)))     ! sqrt(sum(amatr(i,:)*amatr(i,:))/nc)
!   if(verbose) write(*,*)' line#, suml ',i,suml
    if (suml.lt.small) then
       if (abs(vectin(i)) .lt. small) then
!         if(verbose) write(ulog,*)' line #',i,' of Amatrix and Bmatrix is zero; it will be removed!'
       else
          write(ulog,3)' Inconsistency in line i of Amatrix and Bmatrix!',i,suml,vectin(i)
          write(ulog,*)' It will be removed; increase range of FCs, or include higher orders of FCs'
       endif
! record the line index for later removal
       zz(i) = 0
       nz = nz + 1  ! nz counts the number of zero lines
    endif
 enddo
 write(ulog,*)' number of lines=0 in matin are=',nz
 write(*   ,*)' number of lines=0 in matin are=',nz
 nlout=nl-nz
 write(ulog,4)'After removing zeros: NEW # of lines=',nlout
 write(*,4)'After removing zeros: NEW # of lines=',nlout
 if (sum(zz).ne.nl-nz) then
    write(ulog,*)'REMOVE ZEROS: Inconsistency! sum(zz).ne. nl.nz ',sum(zz),nl-nz
    write(*   ,*)'REMOVE ZEROS: Inconsistency! sum(zz).ne. nl.nz ',sum(zz),nl-nz
    stop
 endif

 if(allocated(matout)) deallocate (matout)
 if(allocated(vectout)) deallocate (vectout)
 allocate(matout(nlout,nc),vectout(nlout))

 write(*   ,*)'REMOVE ZEROS: matout,vectout allocated; nlout,nc=',nlout,nc

 j = 0
 do i=1,nl
    if (zz(i).ne.0) then
       j = j+1
       matout(j,:) = matin(i,:)
       vectout(j ) = vectin(i)
    endif
 enddo

! deallocate(matin,vectin)


 3 format(a,i9,2(1x,g11.4))
 4 format(a,4(i9))

 end subroutine remove_zeros

  end module svd_stuff

!===========================================================
