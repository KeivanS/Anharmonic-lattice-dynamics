!============================================================
 module gruneisen_rules
!! Mode Gruneisen parameters, one routine for FOCEX and THERMACOND (2026-11-01; FOCEX's formula):
!!   gamma_s(q) = -<e_s|W(q)|e_s> / (6 omega_s^2)
!!   W(q)_{ta,t'b} = sum_{j in t', k} Psi(ta,jb,kg) y_k^g exp(i q.r_j) / sqrt(m_t m_t')
!! y_k = x_k + xi(k0): position plus the internal relaxation of its sublattice per unit isotropic strain,
!! xi = -Gam Q (internal_strain_trace); r_j = the phase vector of the calling code's dynamical matrix
!! (FOCEX x_j - x_t, THERMACOND the lattice vector R_j). Each code passes its cubic terms once
!! (w_terms_set) and supplies its dynamical matrix (dynmat_routine, non-analytic term included).
!! At Gamma the 3 rigid translations (omega^2 = 0 exactly, set by freq_vel) are 0/0: their gamma is the
!! long-wave limit along the direction of approach (gruneisen_longwave; derivation and the exact-model
!! test in FOCEX/notes/gruneisen_q0.tex, FOCEX/REFERENCE/exact_model/gruneisen_q0.py).
 use constants, only : r15, ci
 implicit none
 private
 public :: dynmat_routine, w_terms_set, internal_strain_trace, mode_gruneisen, gruneisen_longwave

 integer, save :: nw=0
 integer, allocatable, save :: wrow(:),wcol(:)
 real(r15), allocatable, save :: wval(:),wph(:,:)

 abstract interface
    subroutine dynmat_routine(kp,ndn,dyn)
    import r15
    integer, intent(in) :: ndn
    real(r15), intent(in) :: kp(3)
    complex(r15), intent(out) :: dyn(ndn,ndn)
    end subroutine dynmat_routine
 end interface

 contains
!------------------------------------------------------------
 subroutine w_terms_set(n,row,col,val,ph)
!! the cubic terms of W: W(row,col) += val exp(i q.ph), val = Psi y / sqrt(m m')
 integer, intent(in) :: n,row(n),col(n)
 real(r15), intent(in) :: val(n),ph(3,n)
 if (allocated(wrow)) deallocate(wrow,wcol,wval,wph)
 allocate(wrow(n),wcol(n),wval(n),wph(3,n))
 nw=n; wrow=row; wcol=col; wval=val; wph=ph
 end subroutine w_terms_set
!------------------------------------------------------------
 subroutine internal_strain_trace(nat,n2,row,col,fc2,rij,xi)
!! relaxation of the sublattices per unit isotropic strain (eta_ab = eta delta_ab), as FOCEX's
!! get_phi_zeta_Xi: xi(s) = -sum_s' Gam(s,s') Q(s'), Q(ta) = sum_{j g} Phi(ta,jg) (x_j - x_t)_g,
!! Gam = inverse of Phi(q=0) = sum_R Phi with atom 1 held fixed (xi of atom 1 = 0).
!! row = a+3(t-1), col = g+3(t'-1), rij = x_j - x_t for each FC2 term
 integer, intent(in) :: nat,n2,row(n2),col(n2)
 real(r15), intent(in) :: fc2(n2),rij(3,n2)
 real(r15), intent(out) :: xi(3*nat)
 real(r15) phi0(3*nat,3*nat),qtr(3*nat)
 real(r15), allocatable :: g(:,:),gi(:,:)
 integer t,ga,n

 n=3*nat
 phi0=0; qtr=0
 do t=1,n2
    ga=mod(col(t)-1,3)+1
    phi0(row(t),col(t))=phi0(row(t),col(t))+fc2(t)
    qtr(row(t))=qtr(row(t))+fc2(t)*rij(ga,t)
 enddo
 xi=0
 if (nat.gt.1) then
    allocate(g(n-3,n-3),gi(n-3,n-3))
    g=phi0(4:n,4:n)
    call inverse_real(g,gi,n-3)
    xi(4:n)=-matmul(gi,qtr(4:n))
    deallocate(g,gi)
 endif
 end subroutine internal_strain_trace
!------------------------------------------------------------
 subroutine w_matrix(q,ndn,w)
 integer, intent(in) :: ndn
 real(r15), intent(in) :: q(3)
 complex(r15), intent(out) :: w(ndn,ndn)
 integer t
 w=0
 do t=1,nw
    w(wrow(t),wcol(t))=w(wrow(t),wcol(t))+wval(t)*exp(ci*dot_product(q,wph(:,t)))
 enddo
 end subroutine w_matrix
!------------------------------------------------------------
 subroutine w_taylor(qh,n,ndn,w)
!! coefficient of q^n in W(q qh): sum val (i qh.ph)^n / n!
 integer, intent(in) :: ndn,n
 real(r15), intent(in) :: qh(3)
 complex(r15), intent(out) :: w(ndn,ndn)
 integer t
 real(r15) fact
 fact=1; if (n.eq.2) fact=2
 w=0
 do t=1,nw
    w(wrow(t),wcol(t))=w(wrow(t),wcol(t))+wval(t)*(ci*dot_product(qh,wph(:,t)))**n/fact
 enddo
 end subroutine w_taylor
!------------------------------------------------------------
 subroutine herm_eig(n,a,e,v)
!! eigenvalues (ascending) and eigenvectors of the Hermitian part of a
 integer, intent(in) :: n
 complex(r15), intent(in) :: a(n,n)
 real(r15), intent(out) :: e(n)
 complex(r15), intent(out) :: v(n,n)
 complex(r15) work(max(1,4*n))
 real(r15) rwork(max(1,3*n-2))
 integer info
 v=(a+conjg(transpose(a)))/2
 call zheev('V','U',n,v,n,e,work,4*n,rwork,info)
 if (info.ne.0) write(*,*)'GRUNEISEN_RULES: zheev info=',info
 end subroutine herm_eig
!------------------------------------------------------------
 subroutine gruneisen_longwave(getdyn,qdir,ndn,nat,mass,gsize,a,b,gam,ulog)
!! long-wave limit of the Gruneisen parameters of the 3 acoustic branches along the unit vector qh:
!!   D(q qh) = D0 + q D1 + q^2 D2 + ...  (D0 includes the non-analytic term of direction qh)
!!   W(q qh) = W0 + q W1 + q^2 W2 + ...
!!   T = rigid translations, sqrt(m_t/M) delta; P = 1 - T T^+; G = (P D0 P)^+
!!   A = T+ D2 T - T+ D1 G D1 T                      omega_s^2   = q^2 a_s + O(q^3)
!!   B = T+ W2 T - T+ W1 G D1 T - T+ D1 G W1 T + T+ D1 G W0 G D1 T      <e_s|W|e_s> = q^2 b_s + ...
!!   gamma_s = -b_s/(6 a_s), b_s from B inside each degenerate subspace of A; a ascending (TA, TA, LA).
!! D0, D1, D2 from the code's own D at t = +-h, +-2h, +-3h along qh (h = 1e-3 gsize), even and odd
!! parts fitted exactly to degree 4 and 5; W0, W1, W2 exactly from the cubic terms.
 procedure(dynmat_routine) :: getdyn
 integer, intent(in) :: ndn,nat,ulog
 real(r15), intent(in) :: qdir(3),mass(nat),gsize
 real(r15), intent(out) :: a(3),b(3),gam(3)
 real(r15) qh(3)
 complex(r15) dp(ndn,ndn,3),dm(ndn,ndn,3),d0(ndn,ndn),d1(ndn,ndn),d2(ndn,ndn),ev(ndn,ndn,3),od(ndn,ndn,3)
 complex(r15) w0(ndn,ndn),w1(ndn,ndn),w2(ndn,ndn),g(ndn,ndn),pdp(ndn,ndn),vv(ndn,ndn)
 complex(r15) tt(ndn,3),am(3,3),bm(3,3),u(3,3),sub(3,3),vs(3,3),gd1t(ndn,3),gw1t(ndn,3)
 real(r15) h,m(3,3),mi(3,3),lam(ndn),lmax,mtot,es(3),chk(4),scd,scw
 integer k,i,j,t,al,n1,n2,ns

 qh=qdir/sqrt(dot_product(qdir,qdir))   ! unit vector (a and b scale with |qdir|^2, gamma does not)
 h=1d-3*gsize
 do k=1,3
    call getdyn( k*h*qh,ndn,dp(:,:,k))
    call getdyn(-k*h*qh,ndn,dm(:,:,k))
    ev(:,:,k)=(dp(:,:,k)+dm(:,:,k))/2            ! c0 + c2 t^2 + c4 t^4
    od(:,:,k)=(dp(:,:,k)-dm(:,:,k))/(2*k*h)      ! c1 + c3 t^2 + c5 t^4
 enddo
! [1 u^2 u^4] for u = 1,2,3 (t = u h); its inverse gives (c0, c2 h^2, c4 h^4) and (c1, c3 h^2, c5 h^4)
 do k=1,3
    m(k,:)=(/ 1d0, dble(k*k), dble(k**4) /)
 enddo
 call inverse_real(m,mi,3)
 d0=0; d1=0; d2=0
 do k=1,3
    d0=d0+mi(1,k)*ev(:,:,k)
    d2=d2+mi(2,k)*ev(:,:,k)/h**2
    d1=d1+mi(1,k)*od(:,:,k)
 enddo
 d0=(d0+conjg(transpose(d0)))/2; d1=(d1+conjg(transpose(d1)))/2; d2=(d2+conjg(transpose(d2)))/2
 call w_taylor(qh,0,ndn,w0); call w_taylor(qh,1,ndn,w1); call w_taylor(qh,2,ndn,w2)
 w0=(w0+conjg(transpose(w0)))/2; w1=(w1+conjg(transpose(w1)))/2; w2=(w2+conjg(transpose(w2)))/2

 mtot=sum(mass)
 tt=0
 do t=1,nat
 do al=1,3
    tt(al+3*(t-1),al)=sqrt(mass(t)/mtot)
 enddo
 enddo
! G = pseudo-inverse of P D0 P (the 3 translations are its zero modes)
 pdp=d0-matmul(tt,matmul(conjg(transpose(tt)),d0))
 pdp=pdp-matmul(matmul(pdp,tt),conjg(transpose(tt)))
 call herm_eig(ndn,pdp,lam,vv)
 lmax=maxval(abs(lam))
 g=0; ns=0
 do i=1,ndn
    if (abs(lam(i)).lt.1d-9*lmax) then
       ns=ns+1; cycle
    endif
    do j=1,ndn
       g(:,j)=g(:,j)+vv(:,i)*conjg(vv(j,i))/lam(i)
    enddo
 enddo
 if (ns.ne.3) write(ulog,'(a,i3,a)')' GRUNEISEN_LONGWAVE: WARNING: ',ns, &
 &   ' near-zero eigenvalues of P D0 P (expected 3: a soft optical mode at Gamma?)'
 gd1t=matmul(g,matmul(d1,tt)); gw1t=matmul(g,matmul(w1,tt))
 am=matmul(conjg(transpose(tt)),matmul(d2,tt)) - matmul(conjg(transpose(matmul(d1,tt))),gd1t)
 bm=matmul(conjg(transpose(tt)),matmul(w2,tt)) - matmul(conjg(transpose(matmul(w1,tt))),gd1t) &
 &  - matmul(conjg(transpose(matmul(d1,tt))),gw1t) + matmul(conjg(transpose(gd1t)),matmul(w0,gd1t))
 scd=maxval(abs(d0)); scw=max(maxval(abs(w0)),tiny(1d0))
 chk(1)=maxval(abs(matmul(d0,tt)))/scd; chk(2)=maxval(abs(matmul(w0,tt)))/scw
 chk(3)=maxval(abs(matmul(conjg(transpose(tt)),matmul(d1,tt))))/scd
 chk(4)=maxval(abs(matmul(conjg(transpose(tt)),matmul(w1,tt))))/scw
 write(ulog,'(a,3f9.5,a,4es10.2)')' GRUNEISEN_LONGWAVE: qh=',qh, &
 &   '  |D0 T|, |W0 T|, |T+D1T|, |T+W1T| (relative, ~0)=',chk
! a and the degenerate subspaces of A; B diagonalized inside each
 call herm_eig(3,am,a,u)
 n1=1
 do while (n1.le.3)
    n2=n1
    do while (n2.lt.3)
       if (abs(a(n2+1)-a(n1)).gt.1d-6*maxval(abs(a))) exit
       n2=n2+1
    enddo
    sub(1:n2-n1+1,1:n2-n1+1)=matmul(conjg(transpose(u(:,n1:n2))),matmul(bm,u(:,n1:n2)))
    call herm_eig(n2-n1+1,sub(1:n2-n1+1,1:n2-n1+1),es(1:n2-n1+1),vs(1:n2-n1+1,1:n2-n1+1))
    b(n1:n2)=es(1:n2-n1+1)
    n1=n2+1
 enddo
 do i=1,3
    if (a(i).gt.0) then
       gam(i)=-b(i)/(6*a(i))
    else
       gam(i)=0
       write(ulog,*)' GRUNEISEN_LONGWAVE: WARNING: acoustic omega^2/q^2 <= 0 along qh, gamma set to 0',a(i)
    endif
 enddo
 write(ulog,'(a,3es14.6,a,3f12.6)')' GRUNEISEN_LONGWAVE: omega^2/q^2=',a,'  gamma=',gam
 end subroutine gruneisen_longwave
!------------------------------------------------------------
 subroutine mode_gruneisen(getdyn,nkp,kp,ipath,ndn,nat,mass,gsize,eivl,eivc,grn,ulog)
!! gamma of every mode at every k-point, written in the band order of eivl, eivc (get_frequencies).
!! gamma uses the eigenvectors of D at the point itself (user, 2026-11-01: the 3e-4 shift of freq_vel is
!! only a device to separate degenerate bands), with W diagonalized inside each degenerate set.
!! Those translations get the long-wave limit along the direction gamma_limit uses at that point
!! (band path segment for ipath=1, else uinc).
!! Near Gamma, 0 < |k| V^(1/3) < qrad (user, 2026-11-01: 0.05; V = primitive-cell volume, and
!! V^(1/3) = 2 pi/gsize), the 3 acoustic modes (eigenvectors closest to rigid translations) also take
!! the long-wave value along k/|k|: the standard formula is a q^2 cancellation there.
 use phonon_rules, only : gamma_limit, translation_weight
 procedure(dynmat_routine) :: getdyn
 integer, intent(in) :: nkp,ipath,ndn,nat,ulog
 real(r15), intent(in) :: kp(3,nkp),eivl(ndn,nkp),mass(nat),gsize
 complex(r15), intent(in) :: eivc(ndn,ndn,nkp)
 complex(r15), intent(out) :: grn(ndn,nkp)
 real(r15), parameter :: qrad=0.05d0
 complex(r15) w(ndn,ndn),dk(ndn,ndn),ev(ndn,ndn),sub(ndn,ndn),u(ndn,ndn)
 real(r15) q(3),dir(3),a(3),b(3),g3(3),qk,tw(ndn),pi2,lam(ndn),gm(ndn),ws(ndn),tol,ov,best
 integer ik,la,nt,itr(ndn),iac(3),i,j,nnear,n1,n2,ng,ordc(ndn),ordm(ndn),mm(ndn),ic,im,k,jb
 logical used(ndn),istr(ndn),acou(ndn)

 pi2=8*atan(1d0)
 nnear=0
 do ik=1,nkp
! the mode's own eigenvectors at the point where omega^2 is evaluated (k; at Gamma 1e-6 along dir):
! inside a set of degenerate omega^2, W is diagonalized, so gamma does not depend on the basis
    call gamma_limit(ik,nkp,kp,ipath,q,dir)
    call getdyn(q,ndn,dk)
    call herm_eig(ndn,dk,lam,ev)
    call w_matrix(q,ndn,w)
    w=(w+conjg(transpose(w)))/2
    tol=1d-7*maxval(abs(lam))
    gm=0
    n1=1
    do while (n1.le.ndn)
       n2=n1
       do while (n2.lt.ndn)
          if (abs(lam(n2+1)-lam(n1)).gt.tol) exit
          n2=n2+1
       enddo
       ng=n2-n1+1
       sub(1:ng,1:ng)=matmul(conjg(transpose(ev(:,n1:n2))),matmul(w,ev(:,n1:n2)))
       call herm_eig(ng,sub(1:ng,1:ng),ws(1:ng),u(1:ng,1:ng))
       ev(:,n1:n2)=matmul(ev(:,n1:n2),u(1:ng,1:ng))      ! basis that diagonalizes W in the set
       do i=1,ng
          if (lam(n1+i-1).ne.0) gm(n1+i-1)=-ws(i)/(6*lam(n1+i-1))
       enddo
       n1=n2+1
    enddo
! rigid translations: at Gamma (eivl=0, freq_vel) and, near it, the acoustic modes (user's radius
! 0 < |k| V^(1/3) < qrad, V^(1/3) = 2 pi/gsize): the long-wave limit along the direction of approach
    nt=count(eivl(:,ik).eq.0)
    qk=sqrt(dot_product(kp(:,ik),kp(:,ik)))
    acou=.false.
    if (nt.eq.3 .or. (nt.eq.0 .and. qk.gt.1d-12 .and. qk*pi2/gsize.lt.qrad)) then
       do la=1,ndn
          tw(la)=translation_weight(nat,mass,ev(:,la))
       enddo
       do i=1,3                         ! the 3 modes closest to rigid translations
          iac(i)=maxloc(tw,1); tw(iac(i))=-1
       enddo
       do i=1,2                         ! in order of increasing omega^2 (as a_s)
       do j=i+1,3
          if (lam(iac(j)).lt.lam(iac(i))) then
             la=iac(i); iac(i)=iac(j); iac(j)=la
          endif
       enddo
       enddo
       if (nt.eq.3) then
          call gruneisen_longwave(getdyn,dir,ndn,nat,mass,gsize,a,b,g3,ulog)
       else
          call gruneisen_longwave(getdyn,kp(:,ik)/qk,ndn,nat,mass,gsize,a,b,g3,ulog)
          nnear=nnear+1
       endif
       do i=1,3
          gm(iac(i))=g3(i); acou(iac(i))=.true.
       enddo
    elseif (nt.gt.0) then
       write(ulog,*)' MODE_GRUNEISEN: WARNING: ',nt,' modes with omega=0 at k=',kp(:,ik),'; gamma 0'
    endif
! map onto the code's bands: the code's translations (eivl=0) to the acoustic modes, the rest by omega^2
! (both sorted); inside a degenerate set by the largest overlap with the code's eigenvector
    istr=(eivl(:,ik).eq.0)
    used=.false.; mm=0
    if (nt.eq.3) then
       k=0
       do la=1,ndn
          if (istr(la)) then
             k=k+1; mm(la)=iac(k); used(iac(k))=.true.
          endif
       enddo
    endif
    ic=0
    do la=1,ndn                       ! code bands that are not translations, by increasing eivl
       if (istr(la)) cycle
       ic=ic+1; ordc(ic)=la
    enddo
    do i=1,ic-1
    do j=i+1,ic
       if (eivl(ordc(j),ik).lt.eivl(ordc(i),ik)) then
          la=ordc(i); ordc(i)=ordc(j); ordc(j)=la
       endif
    enddo
    enddo
    im=0
    do k=1,ndn                        ! own modes not yet used, by increasing lam (already sorted)
       if (used(k)) cycle
       im=im+1; ordm(im)=k
    enddo
    i=1
    do while (i.le.ic)
! the run of own modes degenerate with ordm(i)
       j=i
       do while (j.lt.ic)
          if (abs(lam(ordm(j+1))-lam(ordm(i))).gt.tol) exit
          j=j+1
       enddo
       do k=i,j                       ! code band ordc(k): own mode of the run with the largest overlap
          best=-1; jb=ordm(i)
          do n1=i,j
             if (used(ordm(n1))) cycle
             ov=abs(dot_product(ev(:,ordm(n1)),eivc(:,ordc(k),ik)))**2
             if (ov.gt.best) then
                best=ov; jb=ordm(n1)
             endif
          enddo
          mm(ordc(k))=jb; used(jb)=.true.
       enddo
       i=j+1
    enddo
    do la=1,ndn
       grn(la,ik)=gm(mm(la))
       if (abs(grn(la,ik)).gt.1d3) grn(la,ik)=0   ! as before: very large gamma replaced by 0
    enddo
 enddo
 write(ulog,'(a,i6,a,f6.3,a)')' MODE_GRUNEISEN: ',nnear,' k-points with 0 < |k| V^(1/3) <',qrad, &
 &   ': acoustic gamma from the long-wave limit along k'
 end subroutine mode_gruneisen

 end module gruneisen_rules
