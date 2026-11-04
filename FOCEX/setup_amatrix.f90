! The rank-by-rank invariance routines set_translational_inv_constraints_old and
! set_rotational_inv_constraints_old were removed on 2026-11-01 (validated replacements in invariance_ranks.f90;
! they remain in git history). Here: the Huang rows and the force-displacement matrix.
!======================================================================================
 subroutine set_huang_inv_constraints
!! The constraint is: S(al,be,ga,de) - S(ga,de,al,be) = -2 Omega (d_al,be sig_ga,de - d_ga,de sig_al,be),
!! S(al,be,ga,de) = sum_Rtautau' phi(tau,al,r+tau'be)(r+tau'-tau)_ga (r+tau'-tau)_de, Omega the
!! volume of the primitive cell, sig = (1/V) dE/d(strain) the stress of the reference structure
!! (tensile > 0), read from stress.params (3 lines, GPa); without that file sig = 0 and S is
!! invariant under (al,be) <-> (ga,de) as before. Holds with FC1 /= 0 too: checked on exact model
!! FC2 to 2e-13 (wurtzite, monoclinic HfO2, 2026-10-26).
! born_flag 2: the fit sees FC2 - Phi_B (the Born-dipole FCs are subtracted with the forces), so the right-hand
! side becomes b - (S_B(ab,gd) - S_B(gd,ab)), S_B from coulomb_lr (2026-11-03, user).
 use params
 use atoms_force_constants
 use svd_stuff
 use ios
 use lattice, only : volume_r0
 use constants, only : r15
 use born, only : born_flag
 use coulomb_lr, only : lr_brackets
 use invariance_ranks, only : column_of
! use geometry
! use force_constants_module
 implicit none
 integer i,j,al,be,ga,de,t,ti,mterm,cnt,cnt2,g,rnk,voigt,vab,vcd,counter,ired,res,dim_h
 real(r15) huang,rr(3),term,sig(3,3),dl(3,3),schange
 real(r15), allocatable:: ared1d(:),zero(:)
 logical ex, lrcor
 real(r15) sb(3,3,3,3)

 huang_constraints=15  ! this is only for rank 2
 if (include_fc(2) .ne. 1) then        ! FC2 not fitted (absent or read, include_fc=2): the rows have no unknowns
    huang_constraints = 0
    if (allocated(ahuang)) deallocate(ahuang)
    if (allocated(bhuang)) deallocate(bhuang)
    allocate(ahuang(0,nindepfc), bhuang(0))
    write(ulog,*)' HUANG: FC2 is not fitted (include_fc(2)=',include_fc(2),'): no Huang rows'
    return
 endif
 lrcor = (born_flag.eq.2)
 sb = 0
 if (lrcor) then
    call lr_brackets(sb)
    write(ulog,*)' HUANG: born_flag 2, right-hand side corrected by the Born-dipole brackets S_B(ab,gd)-S_B(gd,ab)'
    write(ulog,'(a,es12.4,a)')'  max|S_B| =',maxval(abs(sb)),' eV'
 endif

! stress of the reference structure, GPa -> eV/A^3
 sig = 0
 inquire(file='stress.params', exist=ex)
 if (ex) then
    open(9,file='stress.params',status='old')
    do i=1,3
       read(9,*) sig(i,:)
    enddo
    close(9)
    write(ulog,*)' HUANG: stress of the reference (GPa) read from stress.params:'
 else
    write(ulog,*)' HUANG: no stress.params, the reference is taken as stress free (GPa):'
 endif
 do i=1,3
    write(ulog,'(3f14.6)') sig(i,:)
 enddo
 call pointgroup_average2(sig,schange)
 write(ulog,'(a,es10.2,a)')'  stress averaged over the point group (change',schange,' GPa), used:'
 do i=1,3
    write(ulog,'(3f14.6)') sig(i,:)
 enddo
 sig = sig/160.21766208d0
 dl = 0 ; dl(1,1) = 1 ; dl(2,2) = 1 ; dl(3,3) = 1

   rnk = 2;
!  res = sum(map(1)%ntind(:))
   if(map(1)%ngr.gt.0) then
     res = map(1)%nkeptind 
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
      ared1d = zero  !! K1 :: this comes after the al,be,ga,de loop !!
      atomloop: do i=1,natom_prim_cell
      gloop: do g=1,map(2)%ngr
      do ti =1,map(2)%ntind(g)
         ired = column_of(2,g,ti)                      ! 0 if the term is not kept
         if (ired .eq. 0) cycle
      do t=1,map(2)%nt(g)
         if ( map(2)%gr(g)%iat(1,t).ne.i ) cycle
            j  =  map(2)%gr(g)%iat(2,t)
            rr = atompos(:,j)-atompos(:,i)
            term = map(2)%gr(g)%mat(t,ti) 
         if     ( map(2)%gr(g)%ixyz(1,t).eq.al .and. map(2)%gr(g)%ixyz(2,t).eq.be ) then
               ared1d(ired) = ared1d(ired) + rr(ga)*rr(de)*term
         elseif ( map(2)%gr(g)%ixyz(1,t).eq.ga .and. map(2)%gr(g)%ixyz(2,t).eq.de ) then
               ared1d(ired) = ared1d(ired) - rr(al)*rr(be)*term
         endif
      enddo
      enddo
      enddo gloop
      enddo atomloop

      counter = counter+1
      if(verbose) write(ulog,6)'counter,vab,vcd=',counter,al,be,ga,de,vab,vcd
      if (counter.eq.dim_h+1) then
         write(ulog,*)' DIM_H TOO SMALL, huang_inv=15, counter= ',counter
         write(*   ,*)' DIM_H TOO SMALL, huang_inv=15, counter= ',counter
         stop
      endif
      ahuang(counter,:) = ared1d(:)
      if(verbose) write(umatrx,7) ared1d(res+1:res+map(2)%nkeptind) 
      bhuang(counter) = -2*volume_r0*(dl(al,be)*sig(ga,de) - dl(ga,de)*sig(al,be)) &
 &                      - (sb(al,be,ga,de) - sb(ga,de,al,be))
      if (lrcor) write(ulog,'(a,4i2,a,2es14.6)')'  HUANG row',al,be,ga,de,': -2V(d s - d s), Born-dipole part =', &
 &        -2*volume_r0*(dl(al,be)*sig(ga,de) - dl(ga,de)*sig(al,be)), sb(al,be,ga,de) - sb(ga,de,al,be)
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
! which groups to keep is settled in setup_maps
 use svd_stuff
 use params
 use atoms_force_constants
 use ios
 use geometry
 use constants, only : r15
 use invariance_ranks, only : column_of
 use given_fcs, only : xfull, full_col
 implicit none
 integer, intent(in) :: frc_constr,ncfgs
 real(r15), intent(out) :: afrc(frc_constr,nindepfc),bfrc(frc_constr)
 integer al,g,t,ti,ired,line,rnk,nat,confg,rk,ifull
 integer ga(8),tau(8),j_sc(8),nj(3,8)
 real(r15) ux

 write(*,*) 'entering set_force_displacements routine,nindepfc=',nindepfc
 write(ulog,*)' SETUP_AMATRICES : ************************************'
! l is group index, t: general term index, and ti: independent term index
! this convention is followed in all setting constraint or amat subroutines
! set dimensions of ared matrix according to the number of constraints
!------------------------------------------------------------------
!       FORCE-DISPLACEMENT CONSTRAINTS
!------------------------------------------------------------------

! allocate(aforce(force_constraints,ngr),bforce(force_constraints))
 afrc = 0d0; bfrc = 0d0

 write(ulog,*)'********************************************************'
 write(ulog,*)' Now making the force-displacement part of the A matrix'
 write(ulog,*)' Force_constraints, # of independent FCs      =', frc_constr,nindepfc

! now compute the force constraints on the FCs. The equation is:
! -F_ia = phi1_ia + phi2_{ia,jb} u_jb + 0.5 phi3_{ia,jb,kc} u_jb u_kc + ...
!   iatomterm(1,t)=(0,tau)   iatomterm(2,t)=(n1,tau1)
!   ixyzterm (1,t)=alpha     ixyzterm (2,t)=beta
! Really the (n1,tau1)=(iatomcell(:,iatomterm(2,t)),iatomcell0(iatomterm(2,t)))
  line = 0
  do confg=1,ncfgs
  do nat=1,natom_super_cell
  do al=1,3             ! this is how we order each line

     line = line+1      ! this is basically the line number
     bfrc(line) = -force(al,nat,confg)
! find the corresponding atom in the primitive cell
     tau(1) = atom_sc(nat)%cell%tau
     ga (1) = al
     nj(:,1)= atom_sc(nat)%cell%n
!    write(*,4)'line#,confg,i,alfa,tau,n =',line,confg,nat,al,taui,ni

! now set up ared1d: rank by rank. Column of a fitted term: column_of; a read term (include_fc=2) goes to the
! right-hand side with its value xfull(full_col) (given_fcs; 2026-11-03)
     do rnk=1,maxrank   !********************************************
        if (include_fc(rnk) .eq. 0) cycle

        do g=1,map(rnk)%ngr  ! sum over groups

           do t=1,map(rnk)%nt(g) ! sum over all terms in that group

              if ( tau(1).ne. map(rnk)%gr(g)%iat (1,t)   .or.  &
          &        ga (1).ne. map(rnk)%gr(g)%ixyz(1,t) ) cycle

              ux=1
              do rk=2,rnk
                 tau(rk) = iatomcell0( map(rnk)%gr(g)%iat (rk,t))
                 ga(rk)   =            map(rnk)%gr(g)%ixyz(rk,t)
                 nj(:,rk)= iatomcell(:,map(rnk)%gr(g)%iat (rk,t)) + nj(:,1)
! Identify neighbor j within the SCell, find its displacement
                 call findatom_sc(nj(:,rk),tau(rk),j_sc(rk))
                 ux=ux*displ(ga(rk),j_sc(rk),confg)/(rk-1)
              enddo

              do ti=1,map(rnk)%ntind(g)
                 if ( include_fc(rnk) .eq. 1 ) then
                    ired = column_of(rnk,g,ti)        ! 0: FC2 term not kept
                    if (ired .gt. 0) afrc(line,ired)=afrc(line,ired) + ux*map(rnk)%gr(g)%mat(t,ti)
                 else
                    ifull = full_col(rnk,g,ti)
                    if (ifull .gt. 0) bfrc(line)=bfrc(line) - ux*map(rnk)%gr(g)%mat(t,ti)*xfull(ifull)
                 endif
              enddo
           enddo
        enddo

     enddo

 enddo
 enddo
 enddo


 write(umatrx,*)'**********************************************'
 write(umatrx,*)' Force-displ constraints part of amatrix size is ',line
 write(umatrx,*)'**********************************************'

 write(*,*) 'exiting set_force_displacements routine'

4 format(a,99(1x,i6))
7 format(266(1x,g9.2))

 end subroutine set_force_displacement_matrix
