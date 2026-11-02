!===========================================================

!----------------------------------------------
!----------------------------------------------
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
!--------------------------------------------------
