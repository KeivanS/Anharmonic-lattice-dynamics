!===========================================================

 subroutine non_anal(q,tau,taup,dn,ddn,bflag,charge)
!! this is used for the NA term which will be added to the dynamical matrix (3n0 x 3n0)
!! calculates for pair tau,taup the 3x3 matrix of the NA term excluding mass denominator
!! to be added to the short-range part ; uses hard phase , ASR correction NOT included
 use lattice, only: volume_r0,volume_r
 use params, only : coef
 use atoms_force_constants
 use fourier, only : nrgrid,rgrid,rws_weights,ggrid,nggrid,nreg,rgridreg
 use ewald
 use born 
 use geometry, only : length,trace
 implicit none
 integer, intent(in) :: tau,taup,bflag,charge
 real(r15), intent(in) :: q(3)
 complex(r15), intent(out) :: dn(3,3) ,ddn(3,3,3)
 integer al,be,ga
 real(r15) sf,dsf(3),qeq,dqeq(3),etaew,dta(3)

! if (length(q) < 1e-12_r15) then
!    dn  = cmplx(0d0,0d0)
!    ddn = cmplx(0d0,0d0)
!    return
! endif
!  write(*,'(a,999(f6.2))') 'kp_nonanal=', q

 etaew = eta

 if(bflag.le.0) then 

     dn=0; ddn=0

 elseif(mod(bflag,10).eq.1) then  ! naq with damp=etaew and sum over -np,np shells

     etaew = sqrt(trace(epsil)/3d0) / volume_r0**0.3333
     call dyn_coulomb_pure_step(q,tau,taup,etaew,dn,ddn,np)

 elseif(mod(bflag,10).eq.4 .or. mod(bflag,10).eq.5 ) then ! use full G-space ewald summation)  

     call ewald_2nd_deriv_Gnew(q,tau,taup,ng_ewald,g_ewald,etaew,dn,ddn)

 elseif(mod(bflag,10).eq.6 ) then ! use structure factor instead of Parlinski

  !  call dyn_coulomb_structfactor(q,tau,taup,dn,ddn)
     call ewald_2nd_deriv_hat(q,tau,taup,nr_ewald,ng_ewald,r_ewald,g_ewald,etaew,dn,ddn)

 elseif(mod(bflag,10).eq.9 ) then ! this is the pure coulomb FT term

     call naq(q,tau,taup,1000000d0,dn,ddn) 

 else ! mod(bflag,10).eq.2 3 7 8  ) 

     call ewald_2nd_deriv_hat(q,tau,taup,nr_ewald,ng_ewald,r_ewald,g_ewald,etaew,dn,ddn)

 endif

 if (charge.ne.0) then
    dn = matmul(matmul(atom0(tau)%charge,dn),transpose(atom0(taup)%charge)) 
    do ga=1,3
       ddn(:,:,ga)=matmul(matmul(atom0(tau)%charge,ddn(:,:,ga)),transpose(atom0(taup)%charge)) 
    enddo
 endif

! if(length(q).lt.0.1 .and. tau.eq.1 .and. taup.eq.1) then
! write(*,3) 'BEFORE ASR: dn(1,1),dn(2,2),dn(3,3)=',real(dn(1,1)),real(dn(2,2)),real(dn(3,3))
! endif

3 format(a,9(1x,g11.4))

 end subroutine non_anal
!----------------------------------------------
 subroutine non_analq0(dn0,bflag,charge)
!! iHandles the q=0 diagonal block corrections to the NA terms (3n0,3,3)
!! to be subtracted from any q.ne.0 component; uses hard phase 
 use lattice, only: volume_r0,volume_r
 use params, only : coef
 use atoms_force_constants
 use fourier, only : nrgrid,rgrid,rws_weights,ggrid,nggrid,nreg,rgridreg
 use ewald
 use born
 use geometry, only : length,trace
 use ios, only : ulog, write_out
 implicit none
 integer, intent(in) :: bflag,charge
 real(r15), intent(out) :: dn0(natom_prim_cell,3,3)
 integer al,be,ga,tau,taup
 real(r15) q0(3),etaew
 complex(r15) dn(3,3),ddn(3,3,3) 

 q0=0.0_r15
 dn =cmplx(0.0_r15,0.0_r15)
 dn0=cmplx(0.0_r15,0.0_r15)
 do tau=1,natom_prim_cell
    do taup=1,natom_prim_cell
     
       etaew = eta  

       if(mod(bflag,10).eq.1) then  

          etaew = sqrt(trace(epsil)/3d0) / volume_r0**0.3333
          call dyn_coulomb_pure_step(q0,tau,taup,etaew,dn,ddn,np)
  write(*,4) 'dn0=',dn

       elseif(mod(bflag,10).eq.4 .or. mod(bflag,10).eq.5 ) then ! use full G-space ewald summation)  

          call ewald_2nd_deriv_Gnew(q0,tau,taup,ng_ewald,g_ewald,etaew,dn,ddn)

       elseif(mod(bflag,10).eq.6 ) then ! use structure factor instead of Parlinski

          call dyn_coulomb_structfactor(q0,tau,taup,dn,ddn)

       elseif(mod(bflag,10).eq.9 ) then ! this is the pure coulomb FT term

          call naq(q0,tau,taup,1000000d0,dn,ddn) 

       else !(mod(bflag,10).eq.2 .or. mod(bflag,10).eq.3 .or. mod(bflag,10).eq.6 ) then ! 

          call ewald_2nd_deriv_hat(q0,tau,taup,nr_ewald,ng_ewald,r_ewald,g_ewald,etaew,dn,ddn)

       endif

       if (charge.ne.0) then
          dn = matmul(matmul(atom0(tau)%charge,dn),transpose(atom0(taup)%charge)) 
       endif
       dn0(tau,:,:)=dn0(tau,:,:)+real(dn)

    enddo
    write(ulog,*)'NON_ANALQ0: ========== tau, ',tau
    call write_out(ulog,'dn0 ',dn0(tau,:,:))
 enddo

4 format(a,9(1x,2f9.4))
 end subroutine non_analq0
!----------------------------------------------
 subroutine dyn_coulomb_pure_step(q,tau,taup,damp,dyn,ddyn,np)
!! calculates the q component of the Coulomb dynamical matrix to add if born_flag\=0 
! dyn(tau,taup) = [z(tau)^al,ga.q^ga][z(taup)^be.de.q^de] * sf(q)/eps0/(q.epsr.q)/volume_r0
! DOES NOT INCLUDE THE MASS DENOMINATOR, but has the soft phase, so we remove it and adopt the step phase
! use atoms_force_constants
 use geometry
 use lattice, only :  volume_r0,g01,g02,g03 ,g0ws26 !, volume_r, cart2red_g ,r0ws26
 use ios, only : ulog,write_out
 use ewald
! use born, only : born_flag
 implicit none
 integer, intent(in) :: tau,taup,np
 real(r15), intent(in) :: q(3),damp
 complex(r15), intent(out) :: dyn(3,3),ddyn(3,3,3) 
 integer i1,i2,i3 
 real(r15) gpq(3)
 complex(r15) d3(3,3) ,dd3(3,3,3) 

! call write_out(6,'purestep   q  =',q)
! damp= sqrt(2 * trace(epsil))/volume_r0**0.3333  ! so that q.eps.q/(4*damp^2) = q^2 a^2/4  
 call naq(q,tau,taup,damp,dyn,ddyn) ! this is in step phase
! call write_out(6,'after naq: purestep   q  =',q)
! call write_out(6,'first naq   dn  =',dyn)
 do i1=-np,np
 do i2=-np,np
 do i3=-np,np
    if ( i1*i1 + i2*i2 + i3*i3 .eq.0 ) cycle
    gpq=v2a( i1*g01 +i2*g02 + i3*g03 ) + q
    call naq(gpq,tau,taup,damp,d3,dd3) ! this is in step phase
!  write(6,3)'next naq 1,i2,i3,  dn  =',i1,i2,i3,d3
    dyn=dyn+d3
    ddyn=ddyn+dd3
 enddo
 enddo
 enddo

! call write_out(6,'END of purestep   q  =',q)
! dyn= matmul(matmul(atom0(tau)%charge,dyn),transpose(atom0(taup)%charge))
! do ga=1,3
!    ddyn(:,:,ga)= matmul(matmul(atom0(tau)%charge,ddyn(:,:,ga)),transpose(atom0(taup)%charge))
! enddo
! call write_out(6,'purestep Z*dyn*Z=',dyn)

2 format(99(1x,g11.4))
3 format(a,3i3,99(1x,g11.4))

 end subroutine dyn_coulomb_pure_step
!----------------------------------------------
 subroutine dyn_na5D_on_grid(ng,grid) 
!! calculates the NA part of the dynamical matrix (on G* grid )
 use atoms_force_constants, only : natom_prim_cell
 use born, only : dyn_na, dyn_naq0,born_flag
 use ios, only : ulog,write_out
 use params, only: verbose
 use lattice, only : cart2red
 use constants, only : r15,ci
 implicit none
 integer, intent(in):: ng
 real(r15), intent(in):: grid(3,ng)
 integer tau,taup,g
 real(r15) q0(3)
 complex(r15) dn0(3,3),ddn0(3,3,3),asr

   write(ulog,*)'DYN_NA5D_ON_GRID: BORN_FLAG, ng ',born_flag,ng

! get q=0 NA term; has to be subtracted from the NA term at arbitrary q to enforce ASR
!  if(allocated(dyn_naq0)) deallocate(dyn_naq0)
!  allocate( dyn_naq0(natom_prim_cell,3,3) )
  if(allocated(dyn_na  )) deallocate( dyn_na )
  allocate( dyn_na(natom_prim_cell,natom_prim_cell,3,3,ng) )
!  call non_analq0(dyn_naq0,born_flag,1)  ! calculates the q=0 component; born charges multiplied

! calculate the Long-range NA term on the G* grid
   dyn_na=cmplx(0d0,0d0)
   do tau =1,natom_prim_cell
     write(*,*)'DYN_na5D_on_grid: subtracting dyn_naq0 from dyn_na tau=',tau
     do g=1,ng
        do taup=1,natom_prim_cell
           call non_anal(grid(:,g),tau,taup,dn0,ddn0,born_flag,1)
           dyn_na(tau,taup,:,:,g)=dn0 
        enddo
        dyn_na(tau,tau,:,:,g)=dyn_na(tau,tau,:,:,g) - dyn_naq0(tau,:,:) 
     enddo
   enddo

   do g=1,ng
   do tau =1,natom_prim_cell
   do taup=1,natom_prim_cell
      write(ulog,3)'DYN_na5D_on_grid: tau,taup,grid,Gvector=',tau,taup,g,grid(:,g),cart2red(grid(:,g),'g')
      call write_out(ulog,'dyn_na(tau,taup,:,:,g) ',dyn_na(tau,taup,:,:,g))
   enddo
   enddo
   enddo

   call check_herm_sym_G(natom_prim_cell,ng,grid,dyn_na)

   write(*,*)'EXITING DYN_NA5D_ON_GRID'

3 format(a,3i4,9(1x,f9.4))

 end subroutine dyn_na5D_on_grid
!----------------------------------------------
 subroutine set_dynamical_NA_3N3N(kp,dynmat,ndim,ddyn)
!! uses the hard phase convention, with mass term
 use params
 use atoms_force_constants
 use geometry
 use born
 use fourier
 use constants, only : r15,ci
 use ios, only : ulog
 use lattice, only : cart2red
 implicit none
 integer, intent(in) :: ndim
 real(r15), intent(in) :: kp(3)
 complex(r15), intent(out) :: dynmat(ndim,ndim) ,ddyn(ndim,ndim,3) 
 complex(r15) junk,dd(3),dn(3,3),ddn(3,3,3),djunk(3),dn0(natom_prim_cell,3,3)
 integer tau,taup,al,be,ga,ia,jb,igrid
 real(r15) mi,mj,rr(3),rfold(3),q0(3)

 if(born_flag.le.0) return

! fourier transform the real-space fcs stored in phi_sr and add NA term
! use the convoluted dyn_na (g=0 removed) and add NA term
! dn0 = 0  ! used for ASR
! q0 = 0
! do tau=1,natom_prim_cell
! do taup=1,natom_prim_cell
!    call non_anal(q0,tau,taup,dn,ddn,born_flag)
!    dn0 (tau,:,:)= dn0 (tau,:,:) + dn(:,:) 
! enddo
! enddo
! write(*,*)'entering set_dynmat_NA_3N3N ',cart2red(kp,'g')
 do tau=1,natom_prim_cell
 do taup=1,natom_prim_cell
    
    dn=0;ddn=0
    call non_anal(kp,tau,taup,dn,ddn,born_flag,1)
    if (tau.eq.taup) then
      dn=dn-dyn_naq0(tau,:,:) 
    endif
!   if(length(kp).lt.0.1 .and. tau.eq.1 .and. taup.eq.1) then
!     write(*,3) 'AFTER  ASR: dn(1,1), dn(2,2), dn(3,3)=', real(dn(1,1)), real(dn(2,2)), real(dn(3,3))
!   endif

    mi = atom0(tau )%mass
    mj = atom0(taup)%mass

    do al=1,3
    do be=1,3

       ia=al+3*(tau -1)
       jb=be+3*(taup-1)
       dynmat(ia,jb)= (dn(al,be))  /sqrt(mi*mj)
       ddyn(ia,jb,:)= (ddn(al,be,:))   /sqrt(mi*mj)

    enddo
    enddo
 
 enddo
 enddo

3 format(a,9(1x,g11.4))

 end subroutine set_dynamical_NA_3N3N

!===========================================================

subroutine dyn_na5D_on_fine_grid(nggrid_fine, ggrid_fine, dyn_na_fine)
!! Computes dyn_na on a fine reciprocal grid
!! This is similar to your existing dyn_na5D_on_grid but for fine grid
 use atoms_force_constants, only: natom_prim_cell
 use born
 use constants, only: r15
 use geometry, only: length
 use ios, only : ulog,write_out
 implicit none
 
 integer, intent(in) :: nggrid_fine
 real(r15), intent(in) :: ggrid_fine(3, nggrid_fine)
 complex(r15), intent(out) :: dyn_na_fine(natom_prim_cell, natom_prim_cell, 3, 3, nggrid_fine)
 
 integer :: ig, tau, taup
 real(r15) :: q(3),q0(3)
 complex(r15) :: dn(3,3), ddn(3,3,3),dnq0(natom_prim_cell,3,3)
 
 dyn_na_fine = cmplx(0.0_r15, 0.0_r15)
 q0=0.0_r15
 do tau = 1, natom_prim_cell
    dnq0 = cmplx(0.0_r15, 0.0_r15)
    do taup= 1, natom_prim_cell
       call non_anal(q0, tau, taup, dn, ddn, born_flag, 1)
       dnq0(tau, :, :)= dnq0(tau, :, :)+dn
    enddo
    write(ulog,*)'FINE_GRID: dyn_na(q=0) for tau= ',tau
    call write_out(ulog,' d_na(q=0,tau) ',dnq0(tau, :, :))
 enddo

 do ig = 1, nggrid_fine
    q = ggrid_fine(:, ig)
    
!   ! Skip q=0 (or very small q)
!   if (length(q) < 1d-6) cycle
    
    do tau = 1, natom_prim_cell
    do taup = 1, natom_prim_cell
       
       ! Compute NA term at this q
       call non_anal(q, tau, taup, dn, ddn, born_flag, 1)
       
       ! Apply ASR correction for diagonal blocks
       if (tau == taup) then
          dn = dn - dnq0(tau, :, :)
       endif
       
       dyn_na_fine(tau, taup, :, :, ig) = dn
       
    enddo
    enddo
 enddo
 
end subroutine dyn_na5D_on_fine_grid

!===========================================================

subroutine compute_phi_na_realspace(phi_na, nrgrid_xtnd, rgrid_xtnd)
!! Computes phi_na(tau,taup,al,be,R) on extended real-space grid
!! by Fourier transforming dyn_na(G*) computed on a fine reciprocal grid
 use atoms_force_constants, only : natom_prim_cell
 use lattice, only : g01,g02,g03,gs1,gs2,gs3,g0ws26,prim_to_cart,volume_r,volume_r0
 use born
 use geometry
 use fourier, only: nggrid, ggrid, gws_weights
 use constants, only: r15, ci, pi
 use ios
 implicit none
 
 integer, intent(in) :: nrgrid_xtnd
 real(r15), intent(in) :: rgrid_xtnd(3, nrgrid_xtnd)
 complex(r15), intent(out) :: phi_na(natom_prim_cell, natom_prim_cell, 3, 3, nrgrid_xtnd)
 
 ! Fine grid variables
 integer :: nggrid_fine, m_factor
 real(r15), allocatable :: ggrid_fine(:,:), gws_weights_fine(:)
 complex(r15), allocatable :: dyn_na_fine(:,:,:,:,:)
 
 ! Grid generation variables
 type(vector) :: gsp1, gsp2, gsp3
 real(r15) :: matg(3,3)
 real(r15), allocatable :: grd_temp(:,:), weig_temp(:)
 
 ! Loop variables
 integer :: ig, ir, tau, taup, al, be
 integer :: ir_zero
 real(r15) :: R(3),asr
 
 write(ulog,*) 'COMPUTE_PHI_NA_REALSPACE: Starting...'

 ! ========================================
 ! Generate fine grid 
 ! ========================================
 
 m_factor = 2
 gsp1 = gs1 / real(m_factor, r15)
 gsp2 = gs2 / real(m_factor, r15)
 gsp3 = gs3 / real(m_factor, r15)

 write(ulog,*) 'Original gs1:', gs1
 write(ulog,*) 'Fine gsp1:', gsp1
 write(ulog,*) 'Refinement factor m_factor:', m_factor
 
 ! Matrix for converting to reduced coordinates
 matg = transpose(prim_to_cart) / (2*pi)
 
 ! Allocate with generous upper bound
 nggrid_fine = nggrid * (m_factor**3) + 1000
 allocate(grd_temp(3, nggrid_fine))
 allocate(weig_temp(nggrid_fine))
 
 ! Generate the fine grid within the FBZ
 ! gsp1,gsp2,gsp3 are the fine primitive basis
 ! g01,g02,g03 define the FBZ boundary
 call make_grid_weights_WS(gsp1, gsp2, gsp3, g01, g02, g03, matg, &
       &                  nggrid_fine, grd_temp, weig_temp, g0ws26)
 
 write(ulog,*) 'Original ggrid size:', nggrid
 write(ulog,*) 'Fine ggrid size:', nggrid_fine
 write(ulog,*) 'Ratio:', real(nggrid_fine)/real(nggrid)
 
 ! Allocate final fine grid arrays
 allocate(ggrid_fine(3, nggrid_fine))
 allocate(gws_weights_fine(nggrid_fine))
 allocate(dyn_na_fine(natom_prim_cell, natom_prim_cell, 3, 3, nggrid_fine))
 
 ! Copy from temporary arrays
 ggrid_fine = grd_temp(:, 1:nggrid_fine)
 gws_weights_fine = weig_temp(1:nggrid_fine)
 
 ! Apply volume normalization (like in make_grids)
 gws_weights_fine = gws_weights_fine * m_factor**3 * (volume_r0 / volume_r)
 
 write(ulog,*) 'Sum of fine gws_weights:', sum(gws_weights_fine)
 write(ulog,*) 'Normalizing to 1 anyways '

 gws_weights_fine = gws_weights_fine / sum(gws_weights_fine)
 
 deallocate(grd_temp, weig_temp)
 
 ! ========================================
 ! Step 3: Compute dyn_na on fine grid
 ! ========================================
 write(ulog,*) 'Computing dyn_na on fine grid...'

 call dyn_na5D_on_fine_grid(nggrid_fine, ggrid_fine, dyn_na_fine)
 
 ! ========================================
 ! Step 4: Fourier transform to extended real-space grid
 ! ========================================
 write(ulog,*) 'Fourier transforming to extended real space...'
 
 phi_na = cmplx(0.0_r15, 0.0_r15)
 
 do ir = 1, nrgrid_xtnd
    R = rgrid_xtnd(:, ir)
    
    do tau = 1, natom_prim_cell
    do taup = 1, natom_prim_cell
    do al = 1, 3
    do be = 1, 3
       
       ! Sum over fine reciprocal grid
       do ig = 1, nggrid_fine
          phi_na(tau, taup, al, be, ir) = phi_na(tau, taup, al, be, ir) + &
             exp(-ci * dot_product(ggrid_fine(:,ig), R)) * &
             dyn_na_fine(tau, taup, al, be, ig) * &
             gws_weights_fine(ig)
       enddo
       
    enddo
    enddo
    enddo
    enddo
    
    write(*,4)'i,R, phi_na(11xx,R),phi(12xx,R)=',ir,length(R),phi_na(1,1,1,1,ir),phi_na(1,2,1,1,ir)
    ! Progress indicator
    if (mod(ir, max(nrgrid_xtnd/10, 1)) == 0) then
       write(ulog,*) '  Progress:', ir, '/', nrgrid_xtnd
    endif
 enddo
 
 write(ulog,*) 'COMPUTE_PHI_NA_REALSPACE: Done!'
 write(ulog,*) 'Max |phi_na|:', maxval(abs(phi_na))
 write(ulog,*) 'Min |phi_na|:', minval(abs(phi_na))
 
 ! Check if phi_na is real (should be, approximately)
 write(ulog,*) 'Max imag(phi_na) (should be ~0):', maxval(abs(aimag(phi_na)))

4 format(a,i5,9(1x,g11.4))
 
 deallocate(ggrid_fine, gws_weights_fine, dyn_na_fine)


! Find R=0 in extended grid
 ir_zero = 0
 do ir = 1, nrgrid_xtnd
    if (length(rgrid_xtnd(:,ir)) < 1d-6) then
       ir_zero = ir
       exit
    endif
 enddo

 if (ir_zero == 0) then
    write(ulog,*) 'ERROR: R=0 not found in extended grid!'
    stop
 endif

! Finally, check and impose asr if needed
    do tau = 1, natom_prim_cell
    do al = 1, 3
    do be = 1, 3
 
       asr = sum(phi_na(tau,:,al,be,:)) 
       write(ulog,*) 'asr1 check for phi_na: ',tau,al,be,asr 
       asr = sum(phi_na(:,tau,al,be,:)) 
       write(ulog,*) 'asr2 check for phi_na: ',tau,al,be,asr 
       write(ulog,*) 'imposing asr for tau,al,be=',tau,al,be
       phi_na(tau,tau,al,be,ir_zero) = phi_na(tau,tau,al,be,ir_zero) - asr 

    enddo
    enddo
    enddo

end subroutine compute_phi_na_realspace

!--------------------------------------------------
 subroutine phi_na_g2r(phi_na,nr,rgrid,ng,ggrid,weighg)
!! calculates the non-analytical Coulomb FCs in real space grid within the supercell from its Fourier transform
 use lattice, only : volume_r0,cart2red
 use constants, only : r15,ci
 use geometry, only : length,v2a
 use born , only : epsil,epsinv
 use atoms_force_constants, only : atom0,natom_prim_cell,natom_super_cell
 use ios , only : ulog,write_out
 use params, only : coef,tolerance 
 implicit none
 integer, intent(in) :: ng,nr
 integer  ig,tau,taup,al,be,ir,ir0
 real(r15), intent(in) :: ggrid(3,ng),rgrid(3,nr),weighg(ng)
 real(r15), intent(out) :: phi_na(natom_prim_cell,natom_prim_cell,3,3,nr)
 real(r15) geg,dta(3),gscale, gtau(3),gtaup(3),asr

 gscale=6.28/volume_r0**0.3333

 phi_na=0
 do ir=1,nr
 do tau =1,natom_prim_cell
 do taup=1,natom_prim_cell

    dta=rgrid(:,ir)+v2a(atom0(taup)%equilibrium_pos) - v2a(atom0(tau)%equilibrium_pos) 
    do ig=1,ng
       if(length(ggrid(:,ig)).gt.1d-8*gscale) then
          geg=dot_product(ggrid(:,ig),matmul(epsil,ggrid(:,ig)))
          gtau (:)=matmul(atom0(tau )%charge,ggrid(:,ig))
          gtaup(:)=matmul(atom0(taup)%charge,ggrid(:,ig))
          do al=1,3
          do be=1,3
          phi_na(tau,taup,al,be,ir)= phi_na(tau,taup,al,be,ir)+ gtau(al)*gtaup(be) &
    &      /geg  * coef * weighg(ig) * cos(dot_product(ggrid(:,ig),dta)) 
          enddo
          enddo
!       else ! G=0 case
!          phi_na(tau,taup,:,:,ir)= phi_na(tau,taup,:,:,ir)+ matmul(atom0(tau)%charge, &
!  &              matmul(epsinv,transpose(atom0(taup)%charge))) * coef * weighg(ig)  
       endif
    enddo

 enddo
 enddo
 enddo

! apply ASR
 do tau =1,natom_prim_cell
 do al=1,3
 do be=1,3
    asr= 0
    do taup=1,natom_prim_cell
    do ir=1,nr
       if(tau == taup .and. length(rgrid(:,ir)) < tolerance) cycle 
       asr=asr+phi_na(tau,taup,al,be,ir)
    enddo
    enddo

 ! Find ir0 where rgrid = 0
    ir0 = 0
    do ir = 1, nr
       if(length(rgrid(:,ir)) < tolerance) then
          ir0 = ir
          exit
       endif
    enddo

    write(*,*)'PHI_NA_G2R: ir0=',ir0 
    phi_na(tau,tau,al,be,ir0) = -asr
 enddo
 enddo
 enddo

! do ir=1,nr
! do tau =1,natom_prim_cell
! do taup=1,natom_prim_cell
!    write(ulog,4)'ir,tau,taup,R,R_red=',ir,tau,taup,rgrid(:,ir),cart2red(rgrid(:,ir),'r')
!    call write_out(ulog,'phi_na ',phi_na(tau,taup,:,:,ir))
! enddo
! enddo
! enddo
3 format(a,99(1x,g14.7))
4 format(a,3i4,99(1x,f8.4))

 end subroutine phi_na_g2r
!--------------------------------------------------
 subroutine phi_na_0(phi_naq0)
!! calculates the non-analytical Coulomb FCs in real space grid within the supercell from its Fourier transform
 use born , only : epsinv
 use atoms_force_constants, only : atom0,natom_prim_cell
 use params, only : coef
 use constants, only : r15,ci
 implicit none
 real(r15), intent(out) :: phi_naq0(natom_prim_cell,natom_prim_cell,3,3)
 integer tau,taup

 do tau =1,natom_prim_cell
 do taup=1,natom_prim_cell
    phi_naq0(tau,taup,:,:) = matmul(matmul(atom0(tau )%charge,epsinv),transpose(atom0(taup)%charge))/3d0 
 enddo
 enddo

 phi_naq0 = phi_naq0 * coef   ! be careful should there be supercell volume there?

 end subroutine phi_na_0
!----------------------------------------------
 subroutine naq(q,tau,taup,eta,dna,ddna)
!! computes the NA term : (4pi/om0) (qxq)/(q.eps.q) exp(-iq.(taup-tau)) exp(-(q.eps.q)/4eta^2) and its q-derivative
!! this is in the hard phase convention
 use constants, only : r15,ci
 use params, only : coef
 use atoms_force_constants
 use geometry, only : v2a,length
 use born, only: epsil, epsinv
 implicit none
 integer, intent(in) :: tau,taup
 real(r15), intent(in) :: q(3),eta 
 complex(r15), intent(out) :: dna(3,3),ddna(3,3,3)
 integer al,be,ga
 integer , external :: delta_k 
 real(r15) dta(3),term,qeq,dum(3),q0(3)
 complex(r15) zz

    if(length(q).lt.1d-12 ) then
        dna=epsinv/3*coef
!        q0=(/1.0e-12_r15,0.0_r15,0.0_r15/) ! to treat q=0 case
!        dna=cmplx(0d0,0d0) ;
        ddna=cmplx(0d0,0d0) ;
        return
    else 
        q0=q
    endif

    dta = v2a(atom0(taup)%equilibrium_pos)-v2a(atom0(tau)%equilibrium_pos)
    zz = exp(-ci*dot_product(q0,dta) )    ! step phase convention
    dum= matmul((epsil+transpose(epsil)),q0)
    qeq = dot_product(q0,matmul(epsil,q0))
    term = coef * exp(-qeq/(4*eta*eta))/qeq  
    do al=1,3
    do be=1,3

       dna(al,be) = q0(al)*q0(be) * term * zz

       do ga=1,3
          ddna(al,be,ga)= -dna(al,be)*(ci*dta(ga) + dum(ga)*(1/qeq + 1/(4*eta*eta))) &
  &     + term * zz * (q0(al)*delta_k(be,ga)+q0(be)*delta_k(al,ga))
       enddo

    enddo
    enddo

 end subroutine naq
!--------------------------------------------------
