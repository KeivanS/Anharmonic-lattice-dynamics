!! Non-analytical (Born charge) part of the dynamical matrix, shared by FOCEX and THERMACOND
!! (moved from FOCEX/aux_na.f90).  born_flag selects the Ewald variant (see FOCEX manual);
!! born_flag 6 uses dyn_coulomb_structfactor, which needs FOCEX's supercell grids.
 subroutine non_anal(q,tau,taup,dn,ddn,bflag,charge)
!! this is used for the NA term which will be added to the dynamical matrix (3n0 x 3n0)
!! calculates for pair tau,taup the 3x3 matrix of the NA term excluding mass denominator
!! to be added to the short-range part ; uses hard phase , ASR correction NOT included
 use lattice, only: volume_r0
 use atoms_force_constants
 use ewald
 use born 
 use geometry, only : trace
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

 subroutine non_analq0(dn0,bflag,charge)
!! iHandles the q=0 diagonal block corrections to the NA terms (3n0,3,3)
!! to be subtracted from any q.ne.0 component; uses hard phase 
 use lattice, only: volume_r0
 use atoms_force_constants
 use ewald
 use born
 use geometry, only : trace
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

 subroutine dyn_coulomb_pure_step(q,tau,taup,damp,dyn,ddyn,np)
!! calculates the q component of the Coulomb dynamical matrix to add if born_flag\=0 
! dyn(tau,taup) = [z(tau)^al,ga.q^ga][z(taup)^be.de.q^de] * sf(q)/eps0/(q.epsr.q)/volume_r0
! DOES NOT INCLUDE THE MASS DENOMINATOR, but has the soft phase, so we remove it and adopt the step phase
! use atoms_force_constants
 use geometry
 use lattice, only :  g01,g02,g03
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

 subroutine set_dynamical_NA_3N3N(kp,dynmat,ndim,ddyn)
!! uses the hard phase convention, with mass term
 use params
 use atoms_force_constants
 use geometry
 use born
 use constants, only : r15,ci
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
