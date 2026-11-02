!============================================================
 module phonon_rules
!! Rules for frequencies, eigenvectors and velocities shared by FOCEX and THERMACOND (FOCEX 9.1).
!! Each code builds and diagonalizes its own dynamical matrix (freq_routine: omega^2 sorted
!! ascending, eigenvectors, Hellmann-Feynman velocities, in that code's units); this module only
!! decides where it is evaluated.
!! - omega at the k-point itself; eigenvectors and velocities at k + 3e-4 |G| u (freq_vel), which
!!   separates degenerate bands so that both are unique (FOCEX, and THERMACOND's band path).
!! - Gamma as the limit q = 1e-6 along a direction (gamma_limit): the non-analytic term has no value
!!   at q = 0 itself.
!! - the 3 acoustic modes at Gamma are the rigid translations: mass-weighted translation weight of
!!   the eigenvector > 0.999 (translation_weight); found this way, not by band index, because an
!!   unstable mode sorts below them.
 use constants, only : r15
 use geometry, only : length
 implicit none
! a fixed direction whose components are mutually incommensurate: sqrt(2)-1, sqrt(3)-1, sqrt(5)-2.
! Shifting along it leaves the point on no symmetry plane or axis of any lattice.
 real(r15), parameter :: uinc(3)= (/ 0.41421356237309505d0,  &
 &                                   0.73205080756887729d0,  &
 &                                   0.23606797749978970d0 /)
 real(r15), parameter :: translation_min=0.999d0
! |k| below which a point is Gamma, for gamma_limit and freq_vel alike (2026-11-02: a band path ended at
! |k| = 2.7e-10, which freq_vel took for Gamma and gamma_limit did not -> long-wave limit along uinc)
 real(r15), parameter :: gamma_tol=1d-5

 abstract interface
    subroutine freq_routine(kp,ndn,vg,eival,eivec)
    import r15
    integer, intent(in) :: ndn
    real(r15), intent(in) :: kp(3)
    real(r15), intent(out) :: eival(ndn),vg(3,ndn)
    complex(r15), intent(out) :: eivec(ndn,ndn)
    end subroutine freq_routine
 end interface

 contains
!------------------------------------------------------------
 subroutine gamma_limit(i,nkp,kp,ipath,q,dir)
!! the point where k-point i of kp is evaluated, and the direction of the eigenvector shift.
!! Away from Gamma: q = kp(:,i), dir = uinc. At Gamma (|k| < gamma_tol): q = 1e-6 dir, with dir the band-path segment
!! (ipath=1; incoming if the previous point differs, else outgoing) or uinc.
 integer, intent(in) :: i,nkp,ipath
 real(r15), intent(in) :: kp(3,nkp)
 real(r15), intent(out) :: q(3),dir(3)

 q=kp(:,i); dir=uinc
 if (length(q).lt.gamma_tol) then
    if (ipath.eq.1 .and. nkp.gt.1) then
       dir=0
       if (i.gt.1) dir=kp(:,i)-kp(:,i-1)
       if (length(dir).lt.1d-12 .and. i.lt.nkp) dir=kp(:,i+1)-kp(:,i)
       if (length(dir).lt.1d-12) dir=uinc
    endif
    dir=dir/length(dir)
    q=1d-6*dir
 endif
 end subroutine gamma_limit
!------------------------------------------------------------
 function translation_weight(nat,mass,ev) result(proj)
!! weight of the eigenvector ev(3*nat) in the subspace of rigid translations (mass weighted):
!! 1 for a pure translation, 0 for a mode orthogonal to all three
 integer, intent(in) :: nat
 real(r15), intent(in) :: mass(nat)
 complex(r15), intent(in) :: ev(3*nat)
 real(r15) proj,mtot
 complex(r15) s
 integer al,tau

 mtot=sum(mass)
 proj=0
 do al=1,3
    s=0
    do tau=1,nat
       s=s+sqrt(mass(tau)/mtot)*ev(al+3*(tau-1))
    enddo
    proj=proj+abs(s)**2
 enddo
 end function translation_weight
!------------------------------------------------------------
 subroutine freq_vel(getf,qq,shdir,gsize,nat,mass,ndn,vgr,evl0,evc0,ulog)
!! frequencies (omega^2, evl0) at qq itself; eigenvectors and Hellmann-Feynman group velocities at
!! qq + 3e-4 gsize shdir/|shdir| (gsize = volume of the reciprocal primitive cell^(1/3)).
!! At the Gamma limit (|qq| < 1e-5, see gamma_limit) the rigid translations get omega^2 = 0 exactly:
!! they are not vibrations (the value before is written to ulog).
 procedure(freq_routine) :: getf
 integer, intent(in) :: nat,ndn,ulog
 real(r15), intent(in) :: qq(3),shdir(3),gsize,mass(nat)
 real(r15), intent(out) :: vgr(3,ndn),evl0(ndn)
 complex(r15), intent(out) :: evc0(ndn,ndn)
 real(r15) qr(3),v2(3,ndn)
 complex(r15) evct(ndn,ndn)
 integer l

 qr = qq + 3d-4*gsize*shdir/length(shdir)
 call getf(qr,ndn,vgr,evl0,evc0)
 call getf(qq,ndn,v2,evl0,evct)
 if (length(qq).lt.gamma_tol) then
    do l=1,ndn
       if (translation_weight(nat,mass,evct(:,l)).gt.translation_min) then
          write(ulog,'(a,i4,a,g12.4,a)')' FREQ_VEL: Gamma, mode',l,' is a rigid translation, omega^2=', &
 &               evl0(l),' set to 0'
          evl0(l)=0
       endif
    enddo
 endif
 end subroutine freq_vel

 end module phonon_rules
