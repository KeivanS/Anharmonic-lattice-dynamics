!============================================================
 subroutine relax_reference
!! Relaxes the reference structure with the fitted force constants and re-expands the FCs about the
!! relaxed one (user, 2026-11-02; default.params line 5 = 1). ALATDYN paper, section on elastic constants:
!!   S_{R tau al} = eta_{al be} (R+tau)_be + u0_{tau al}                      (static displacement)
!!   K = Phi + Psi S + Chi S S/2   (effective FC2; SCOP8: trialfc2_value),  Psi' = Psi + Chi S,
!!   Pi' = Pi + Phi S + Psi S S/2 + Chi S S S/6                               (re-expanded FC1)
!! with S_j - S_tau in place of S_j (same by ASR, independent of the image chosen for an atom).
!! E(eta,u0) is minimized by Newton steps: gradient dE/du0 = Pi'(eta,u0),
!!   dE/deta = V0 sigma0 + int_0^1 [A(K(t q)) eta + Q(K(t q))^T u0] dt   (Simpson, exact up to FC4),
!! Hessian [[A, Q^T],[Q, phi]] of K at the current point (paper: phi, Q, A; reference positions x);
!! sigma0 from stress.params (GPa, tensile > 0), else the FC1 estimate sum_tau Pi_tau tau/V (image-
!! dependent when Pi /= 0). eta symmetric (6 Voigt components), u0(atom 1) = 0.
!! On exit: atompos, lattice vectors, grids strained by M = 1+eta and shifted by u0; fcs = Pi', K, Psi'
!! projected on the irreducible FCs of the same groups; FC4 unchanged. stress.params is then read as 0.
!! Coulomb (born_flag 2, 2026-11-04, user): the fitted FC2/FC3/FC4 are then the short-range part; the Coulomb
!! FCs Phi_C, Psi_C, Chi_C are added. Phi_C (Z* tensors and eps, E = 0, Ewald; coulomb_lr) enters through its
!! lattice sums: phi_c (q = 0 matrix, internal relaxation), q_c (moments, strain-displacement coupling) and a_c
!! (brackets, strain-strain); a deformation F transforms the position factors (q_c F^T, F a_c F^T). Psi_C and
!! Chi_C are term lists of the FC3/FC4 form (psi_c, chi_c): empty (zero) for now, so K_C = Phi_C.
!! Exact Coulomb (born_flag 2 with Z* = z 1, 2026-11-04, user): the route above keeps the Coulomb energy to second order
!! about the reference, which is not rotationally invariant and misses its higher orders (exact HfO2r model with point
!! charges: strain off by 1.4e-4, u0 by 6e-3 A). Instead E = E_SR(X) + E_C(X) - E_C(0) with E_SR the Taylor series of
!! the fitted FCs with FC1 - F_C(0) and stress.params - sigma_C(0) (short-range alone, rotationally invariant) and
!! E_C the Ewald energy of the pair model that the
!! born_flag-2 subtraction uses (ewaldforce: c_ij = tr(Z_i Z_j^T)/3, phi_eps) at the deformed cell and positions
!! (ewald_cell). X = (6 Voigt strains, u0 of atoms 2..), G(X) = dE/dX = 0 solved by modified Broyden (bro90 of SCOP8,
!! module broyden_mix) on f = -H0^-1 G, H0 the Hessian at the reference with Phi_C. dE_SR/du0 = Pi' (pi_prime, short-range
!! terms). dE_SR/deta: the stress formula (elastic_born_huang) holds only at force-free states, so the short-range stress
!! is integrated (dsdt, Coulomb off) along the path where the short-range forces vanish, u0*(eta), and then carried from
!! u0* to u0 at fixed eta by int Q^T du (Q with reference vectors; exact by Simpson, Q is quadratic in t). At that
!! force-free state dE/deta = V sigma F^-T. For anisotropic Z* the second-order route is kept (warning).
 use ios
 use params
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use fourier
 use constants
 use born, only : born_flag
 use coulomb_lr, only : lr_c, lr_brackets, lr_moments, lr_moments_pairs
 use born, only : epsil
 use ewald_cell, only : ewald_cell_efw
 use broyden_mix
 implicit none
 integer, parameter :: vg(2,6)=reshape([1,1, 2,2, 3,3, 2,3, 1,3, 1,2],[2,6])
 integer nprim,r,g,t,ti,i,j,k,l,n,it,np,nu,nvar,al,be,ga,de,a,b,c,d,tau,taup,kk,ll,ier,ip,itr
 integer ofs(8),nfct(4)
 integer, allocatable :: iat1(:),ix1(:),iat2(:,:),ix2(:,:),iat3(:,:),ix3(:,:),iat4(:,:),ix4(:,:)
 integer, allocatable :: g2(:),t2(:),g3(:),t3(:),g1(:),t1(:)
 integer, allocatable :: pidx(:,:),pt(:),pj(:),p3of4(:)
 real(r15), allocatable :: v1(:),v2(:),v3(:),v4(:),x(:,:),kmat(:,:,:),kmat0(:,:,:),pip(:,:),psip(:)
 real(r15) pnet(3)   ! net force removed by pi_prime (from the FC2-FC4 terms when their ASR is not exact)
 real(r15) pnet1(3)  ! sum of the fitted FC1, removed before relaxing
! Coulomb FCs (born_flag 2): harmonic Phi_C through its lattice sums; Psi_C, Chi_C as term lists (empty for now)
 logical coul,zon
 real(r15), allocatable :: phi_c(:,:),q_c(:,:,:),mo_c(:,:,:,:),psi_c(:),chi_c(:),mop_c(:,:,:,:,:)
 complex(r15), allocatable :: c0_c(:,:,:,:)
 integer, allocatable :: iat3c(:,:),ix3c(:,:),iat4c(:,:),ix4c(:,:)
 real(r15) a_c(3,3,3,3),sb_c(3,3,3,3),f1(3,3)
 integer nfct3c,nfct4c
! exact Coulomb route (coulx): short-range force-free path state, Coulomb couplings, gradients at the reference
 logical coulx
 integer nx
 real(r15) fpth(3,3),spth(3,3),e2c,alph,zan,zmx,ztr,zz(3,3)
 real(r15), allocatable :: upth(:,:),cmat(:,:),gc0(:),gu0h(:,:)
 integer, allocatable :: imgsg(:,:)
 integer isg,ntn
 real(r15) tnew(3),tref(3),sgmax
 real(r15), allocatable :: phi_(:,:),q_(:,:,:),u0(:,:),xc(:,:),um(:,:),dum3(:,:),ustart(:,:)
 real(r15) eta(3,3),e6(6),s0(3,3),xi0(3,3),a_(3,3,3,3),gsim(3,3),w(3),dd(3),dl(3),em(3,3),minv(3,3)
 real(r15) gconv,uconv,sgpa,dmax,prmax,rel,v0,fac,pmax0,sig_in(3,3),kdrop,resmax(3),fmax(3),detm
 logical ex
 integer(8), allocatable :: hkey(:)
 integer, allocatable :: hval(:)
 integer nh
 real(r15) wsimp(3),tsimp(3)
 real(r15) ft(3,3),gm(3,3),sig(3,3),crel(3,3,3,3),vol,g6(6),h6(6,6),ek(3,3),el(3,3),dstr(3,3),ss(3,3), &
 &         k1(3,3),k2(3,3),k3(3,3),k4(3,3),t0,fm(3,3),jj,rot(3,3)
 integer ind6(6)
 integer, parameter :: nsub=8
 data wsimp/0.1666666666666667d0,0.6666666666666667d0,0.1666666666666667d0/, tsimp/0d0,0.5d0,1d0/

 nprim=natom_prim_cell
 write(ulog,*)'RELAX_REFERENCE: relaxing the reference structure with the fitted FCs (ALATDYN paper, S=eta(R+tau)+u0)'
 write(*,*)'RELAX_REFERENCE: relaxing the reference structure with the fitted FCs'
 sgpa=1d30*ee*1d-9          ! eV/A^3 -> GPa
 v0=volume_r0

! offsets of each rank in fcs, as in write_output_fcs
 ofs=0; n=0
 do r=1,maxrank
    ofs(r)=n
    if(include_fc(r).eq.1) n=n+map(r)%nkeptind
 enddo

! full term lists of ranks 1-4 (first atom in the primitive cell) with their values
 nfct=0
 do r=1,4
    if(include_fc(r).ne.1) cycle
    do g=1,map(r)%ngr
    do t=1,map(r)%nt(g)
       if(map(r)%gr(g)%iat(1,t).le.nprim) nfct(r)=nfct(r)+1
    enddo
    enddo
 enddo
 allocate(iat1(nfct(1)),ix1(nfct(1)),v1(nfct(1)),g1(nfct(1)),t1(nfct(1)))
 allocate(iat2(2,nfct(2)),ix2(2,nfct(2)),v2(nfct(2)),g2(nfct(2)),t2(nfct(2)))
 allocate(iat3(3,nfct(3)),ix3(3,nfct(3)),v3(nfct(3)),g3(nfct(3)),t3(nfct(3)))
 allocate(iat4(4,nfct(4)),ix4(4,nfct(4)),v4(nfct(4)))
 do r=1,4
    if(include_fc(r).ne.1) cycle
    n=0
    do g=1,map(r)%ngr
    do t=1,map(r)%nt(g)
! terms whose first atom is outside the primitive cell are translations of in-cell ones (skipped, as in
! get_phi_zeta_Xi)
       if(map(r)%gr(g)%iat(1,t).gt.nprim) cycle
       n=n+1
       fac=0
       do ti=1,map(r)%ntind(g)
          if(r.eq.2) then
             if(map(2)%keep(counteri(2,g,ti)).eq.0) cycle
          endif
          fac=fac+map(r)%gr(g)%mat(t,ti)*fcs(ofs(r)+current(r,g,ti))
       enddo
       select case(r)
       case(1); iat1(n)=map(1)%gr(g)%iat(1,t); ix1(n)=map(1)%gr(g)%ixyz(1,t); v1(n)=fac; g1(n)=g; t1(n)=t
       case(2); iat2(:,n)=map(2)%gr(g)%iat(1:2,t); ix2(:,n)=map(2)%gr(g)%ixyz(1:2,t); v2(n)=fac; g2(n)=g; t2(n)=t
       case(3); iat3(:,n)=map(3)%gr(g)%iat(1:3,t); ix3(:,n)=map(3)%gr(g)%ixyz(1:3,t); v3(n)=fac; g3(n)=g; t3(n)=t
       case(4); iat4(:,n)=map(4)%gr(g)%iat(1:4,t); ix4(:,n)=map(4)%gr(g)%ixyz(1:4,t); v4(n)=fac
       end select
    enddo
    enddo
 enddo
 write(ulog,*)'RELAX_REFERENCE: number of FC terms of ranks 1-4 with the first atom in the primitive cell =',nfct

! sum_tau Pi_tau = 0 (rank-1 translational invariance) is exact only with enforce_inv >= 1. A net force
! cannot be relaxed (it would translate the crystal), so it is removed from FC1 before relaxing: the
! average force per atom is subtracted from every atom (user, 2026-11-04). The fitted FC1 is kept in
! fc1_unrelaxed.dat; fc1.dat is then the relaxed FC1, with this correction.
 pnet1=0
 do n=1,nfct(1)
    pnet1(ix1(n))=pnet1(ix1(n))+v1(n)
 enddo
 pnet1=pnet1/nprim
 write(ulog,'(a,3es12.4)')' RELAX_REFERENCE: fitted sum_tau FC1(tau) (eV/A; 0 if the ASR is exact) =',pnet1*nprim
 do n=1,nfct(1)
    v1(n)=v1(n)-pnet1(ix1(n))
 enddo
 pnet1=0
 do n=1,nfct(1)
    pnet1(ix1(n))=pnet1(ix1(n))+v1(n)
 enddo
 write(ulog,'(a,3es12.4)')' RELAX_REFERENCE: sum_tau FC1(tau) after removing the average force per atom   =',pnet1
 if(maxval(abs(pnet1)).gt.1d-12) call warn3(ulog, &
 &  'RELAX_REFERENCE: sum of FC1 not zero after the correction (an atom lacks that FC1 component)',maxval(abs(pnet1)))

! reference positions
 allocate(x(3,natoms))
 x=atompos(:,1:natoms)

! pairs (tau,j) carrying K: those of FC2, then those appearing in FC3/FC4 only
 allocate(pidx(nprim,natoms)); pidx=0
 np=0
 allocate(pt(nfct(2)+nfct(3)+nfct(4)),pj(nfct(2)+nfct(3)+nfct(4)))
 do n=1,nfct(2)
    call addpair(iat2(1,n),iat2(2,n))
 enddo
 nu=np   ! pairs of the FC2 map
 do n=1,nfct(3)
    call addpair(iat3(1,n),iat3(2,n))
 enddo
 do n=1,nfct(4)
    call addpair(iat4(1,n),iat4(2,n))
 enddo
 write(ulog,*)'RELAX_REFERENCE: atom pairs of FC2, and of FC3/FC4 beyond them =',nu,np-nu
 allocate(kmat(3,3,np),kmat0(3,3,np),pip(3,nprim),psip(nfct(3)))

! Coulomb FCs (born_flag 2: the fit saw the forces minus those of Phi_C, so its FC2 is the short-range part)
 zon=.false.
 do tau=1,nprim
    if(maxval(abs(atom0(tau)%charge)).gt.0) zon=.true.
 enddo
 coul = born_flag.eq.2 .and. zon
 f1=0; do i=1,3; f1(i,i)=1; enddo
 allocate(phi_c(3*nprim,3*nprim),q_c(3*nprim,3,3),mop_c(nprim,nprim,3,3,3)); phi_c=0; q_c=0; a_c=0; mop_c=0
 nfct3c=0; nfct4c=0
 allocate(iat3c(3,nfct3c),ix3c(3,nfct3c),psi_c(nfct3c),iat4c(4,nfct4c),ix4c(4,nfct4c),chi_c(nfct4c))
 if(coul) then
    allocate(c0_c(3,3,nprim,nprim),mo_c(nprim,3,3,3))
    call lr_c((/0d0,0d0,0d0/),c0_c)                  ! sum_R Phi_C(t,u;R) without the self term (E = 0)
    do tau=1,nprim
    do taup=1,nprim
       phi_c(3*tau-2:3*tau,3*taup-2:3*taup)=real(c0_c(:,:,tau,taup))
    enddo
    enddo
    do tau=1,nprim            ! self term Phi_C(t,t;0) = -sum of the others (ASR), after all blocks are set
    do taup=1,nprim
       phi_c(3*tau-2:3*tau,3*tau-2:3*tau)=phi_c(3*tau-2:3*tau,3*tau-2:3*tau)-real(c0_c(:,:,tau,taup))
    enddo
    enddo
    call lr_moments(mo_c)                              ! sum_j Phi_C(t a1, j a) (x_j - x_t)_m
    do tau=1,nprim
       do al=1,3
          q_c(al+3*(tau-1),:,:)=mo_c(tau,al,:,:)
       enddo
    enddo
    call lr_moments_pairs(mop_c)                       ! the same per sublattice pair (t,u)
    call lr_brackets(sb_c)                             ! sum_{t,j} Phi_C(t a, j b) r_g r_d (eV)
    do a=1,3; do b=1,3; do c=1,3; do d=1,3
       a_c(a,c,b,d)=-sb_c(a,b,c,d)/2                   ! same index order as a_ of build_aq
    enddo; enddo; enddo; enddo
    write(ulog,'(a,3es12.4)')' RELAX_REFERENCE: Coulomb Phi_C added (born_flag 2): max|phi_c|, max|q_c|, max|a_c| =', &
 &        maxval(abs(phi_c)),maxval(abs(q_c)),maxval(abs(a_c))
    write(ulog,'(a,3es12.4)')' RELAX_REFERENCE: checks of phi_c: max|phi-phi^T|, max|sum_tau phi_c|, max|sum_tau q_c| =', &
 &        maxval(abs(phi_c-transpose(phi_c))),maxval(abs(sum(reshape(phi_c,[3,nprim,3*nprim]),dim=2))), &
 &        maxval(abs(sum(reshape(q_c,[3,nprim,3,3]),dim=2)))
    write(ulog,*)'RELAX_REFERENCE: Coulomb Psi_C, Chi_C: not included yet (zero); number of terms =',nfct3c,nfct4c
    deallocate(c0_c,mo_c)
 endif
 coulx=.false.
 if(coul) then
    zan=0; zmx=0
    do tau=1,nprim
       zz=atom0(tau)%charge; zmx=max(zmx,maxval(abs(zz)))
       ztr=(zz(1,1)+zz(2,2)+zz(3,3))/3
       do i=1,3
          zz(i,i)=zz(i,i)-ztr
       enddo
       zan=max(zan,maxval(abs(zz)))
    enddo
    coulx = zan.le.1d-6*max(1d0,zmx)
    if(coulx) then
       write(ulog,'(a,es10.2,a)')' RELAX_REFERENCE: Z* = z 1 (max|Z - tr Z/3| = ',zan, &
 &         '): Coulomb part exact (Ewald at the deformed cell and positions), solved by Broyden'
    else
       call warn3(ulog,'RELAX_REFERENCE: anisotropic Z*: Coulomb only to second order (Phi_C); max|Z - tr Z/3| =',zan)
    endif
 endif
 if(zon .and. born_flag.ne.0 .and. .not.coul) then
    call warn3(ulog,'RELAX_REFERENCE: Born charges with this born_flag: the relaxation has no Coulomb term; born_flag =', &
 &       dble(born_flag))
 endif

! FC3 term of the first three atoms and directions of each FC4 term (hash table)
 nh=1
 do while(nh.lt.2*max(1,nfct(3)))
    nh=2*nh
 enddo
 allocate(hkey(0:nh-1),hval(0:nh-1)); hkey=-1; hval=0
 do n=1,nfct(3)
    call hput(key3(iat3(:,n),ix3(:,n)),n)
 enddo
 allocate(p3of4(nfct(4))); kk=0
 do n=1,nfct(4)
    p3of4(n)=hget(key3(iat4(1:3,n),ix4(1:3,n)))
    if(p3of4(n).eq.0) kk=kk+1
 enddo
 if(kk.gt.0) write(ulog,*)'RELAX_REFERENCE: FC4 terms without an FC3 term of their first 3 atoms (Psi'' not updated there)=',kk

! reference stress
 call read_stress_params(ulog,sig_in,ex)
 xi0=0
 do n=1,nfct(1)
    do be=1,3
       xi0(ix1(n),be)=xi0(ix1(n),be)+v1(n)*x(be,iat1(n))
    enddo
 enddo
 xi0=(xi0+transpose(xi0))/2
 if(ex) then
    s0=sig_in*v0/sgpa                 ! V0 sigma0 (eV)
    call write_out(ulog,' RELAX_REFERENCE: reference stress from stress.params (GPa) ',sig_in)
 else
    s0=xi0
    call write_out(ulog,' RELAX_REFERENCE: no stress.params; reference stress = sum_tau Pi tau/V (GPa) ',xi0/v0*sgpa)
    call warn3(ulog,'RELAX_REFERENCE: stress from FC1 only (image-dependent if FC1 /= 0); give stress.params',0d0)
 endif

! The strain Hessian is not -1/2 sum K d d (paper's A; only its part symmetric in the two position
! indices is fixed by the FCs): it is taken from FOCEX's elastic constants (Born-Huang brackets + stress
! terms, elastic_born_huang, exact on the exact model at force-free states). So every state used is
! force-free: stage 1 relaxes u0 in the reference cell (the stress changes by int Q^T du, exact: Q is
! image-independent); stage 2 takes Newton steps in a symmetric strain about the current state, with u0
! relaxed at each point, and integrates the stress (PK2 about the step's start) along the step (RK4):
!   dS/dt = J F^-1 F^-1 c_rel(t) F^-T F^-T : (de F + F de)/2,  F = 1 + t de,  sigma = F S F^T / J,
!   c_rel = c(brackets of K(t), stress sigma(t)) - Sym(Q Gamma Q)/V   (relaxed internal coordinates).
 allocate(u0(3,nprim),phi_(3*nprim,3*nprim),q_(3*nprim,3,3),xc(3,natoms),um(3,nprim),dum3(3,nprim),ustart(3,nprim))
 dstr=0; pmax0=0
 ft=0; do i=1,3; ft(i,i)=1; enddo
 u0=0; eta=0
 gm=s0
 call hess6_reference
 if(coulx) then
! Short-range part = the fitted FCs with FC1 and stress of the short-range energy alone: FC1 - Coulomb force and
! stress.params - Coulomb stress at the reference. Only then is it rotationally invariant, which the stress formula
! on its force-free path needs (with the total FC1 and stress the short-range dE/deta came out 0.17 eV per cell off,
! HfO2r_q, 2026-11-05). The Coulomb part is then E_C(X) - E_C(0).
    call coulx_setup
    coul=.false.   ! stage 1 then relaxes u0 with the short-range forces only: start of the path u0*(eta)
    write(ulog,'(a)')' RELAX_REFERENCE: stage 1 with the short-range terms only (start of the short-range force-free path)'
 endif
 write(ulog,'(a)')' RELAX_REFERENCE: stage 1, internal relaxation in the reference cell: iter, max|FC1| (eV/A)'
 do itr=0,50
    call kfill(1d0); call pi_prime
    uconv=0
    if(nprim.gt.1) uconv=maxval(abs(pip(:,2:nprim)))   ! atom 1 fixed: its force is -sum of the others (ASR)
    if(itr.eq.0) pmax0=uconv
    write(ulog,'(a,i4,1x,es12.4)')' RELAX_REFERENCE: ',itr,uconv
    if(uconv.lt.1d-12 .or. nprim.eq.1) exit
    call build_aq(x,.true.,f1)
    call solve_u(dum3)
! stress change over the step, Simpson in t (K quadratic in t up to FC4)
    um=u0
    do k=1,3
       u0=um+tsimp(k)*dum3
       call kfill(1d0); call build_aq(x,.false.,f1)
       do a=1,3
       do b=1,3
          do tau=2,nprim
          do al=1,3
             gm(a,b)=gm(a,b)+wsimp(k)*q_(al+3*(tau-1),a,b)*dum3(al,tau)
          enddo
          enddo
       enddo
       enddo
    enddo
    u0=um+dum3
 enddo
 call write_out(ulog,' RELAX_REFERENCE: antisymmetric part of dE/deta after stage 1 (eV, 0 by rotational invariance)', &
 &              (gm-transpose(gm))/2)
 sig=(gm+transpose(gm))/2/v0                 ! eV/A^3
 call write_out(ulog,' RELAX_REFERENCE: stress after the internal relaxation (GPa) ',sig*sgpa)

 if(coulx) then
    allocate(upth(3,nprim))
    fpth=f1; upth=u0; spth=sig
    call relax_broyden
 else
 write(ulog,'(a)')' RELAX_REFERENCE: stage 2, strain: iter, max|stress| (GPa), max|strain step|, max|FC1| (eV/A)'
 do itr=0,60
    call elastic_at(ft,u0,sig,crel,vol)
    gconv=maxval(abs(sig))*sgpa
    write(ulog,'(a,i4,3(1x,es12.4))')' RELAX_REFERENCE: ',itr,gconv,maxval(abs(dstr)),uconv
    if(gconv.lt.1d-8) exit
    if(itr.eq.60) then
       call warn3(ulog,'RELAX_REFERENCE: no convergence in 60 strain steps; max|stress| (GPa) =',gconv)
       call warn3(   6,'RELAX_REFERENCE: no convergence in 60 strain steps; max|stress| (GPa) =',gconv)
       exit
    endif
! Newton step: H_kl = V [c_rel + d^2(sigma:e^2/2)] (E^k, E^l), g_k = V sigma:E^k
    do k=1,6
       call eta_of(unit6(k),ek)
       g6(k)=vol*sum(sig*ek)
       do l=1,6
          call eta_of(unit6(l),el)
          h6(k,l)=vol*(contract(crel,ek,el)+sum(sig*(matmul(ek,el)+matmul(el,ek)))/2)
       enddo
    enddo
    call ludcmp(h6,6,6,ind6,fac)
    call lubksb(h6,6,6,ind6,g6)
    call eta_of(-g6,dstr)
! stress along the step
    ss=sig; ustart=u0
    do kk=1,nsub
       t0=dble(kk-1)/nsub
       call dsdt(t0,ss,ustart,k1)
       call dsdt(t0+0.5d0/nsub,ss+k1/(2*nsub),ustart,k2)
       call dsdt(t0+0.5d0/nsub,ss+k2/(2*nsub),ustart,k3)
       call dsdt(t0+1d0/nsub,ss+k3/nsub,ustart,k4)
       ss=ss+(k1+2*k2+2*k3+k4)/(6*nsub)
    enddo
    fm=dstr; do i=1,3; fm(i,i)=fm(i,i)+1; enddo
    jj=det3(fm)
    sig=matmul(fm,matmul(ss,transpose(fm)))/jj
    ft=matmul(fm,ft)
    u0=ustart
    call relax_u(ft,u0)
! remove the rotation (polar decomposition ft = R U): ft -> U, u0 -> R^T u0, sigma -> R^T sigma R
    call polar(ft,rot)
    ft=matmul(transpose(rot),ft)
    u0=matmul(transpose(rot),u0)
    sig=matmul(transpose(rot),matmul(sig,rot))
    sig=(sig+transpose(sig))/2
 enddo
 endif
 eta=ft; do i=1,3; eta(i,i)=eta(i,i)-1; enddo
 call kfill(1d0); call pi_prime
 if(coulx) call add_coul_pip   ! FC1 at the relaxed state = total gradient (short-range + exact Coulomb)
 uconv=0
 if(nprim.gt.1) uconv=maxval(abs(pip(:,2:nprim)))
 e6=[eta(1,1),eta(2,2),eta(3,3),eta(2,3),eta(1,3),eta(1,2)]

 write(ulog,'(a,6(1x,es20.12))')' RELAX_REFERENCE: relaxed strain eta (Voigt xx yy zz yz xz xy) =',e6
 do tau=1,nprim
    write(ulog,'(a,i4,3(1x,es20.12))')' RELAX_REFERENCE: tau, u0(tau) (A) =',tau,u0(:,tau)
 enddo
 write(ulog,'(a,es12.4,a,es12.4)')' RELAX_REFERENCE: max|FC1| (atoms 2..) before ',pmax0,' after ',uconv
 write(ulog,'(a,3(1x,es12.4))')' RELAX_REFERENCE: net force per atom from the FC2-FC4 terms at the relaxed state, '// &
 &     'removed (eV/A; 0 if their ASR is exact) =',pnet
 write(ulog,'(a,es12.4)')' RELAX_REFERENCE: max|stress| after (GPa) =',maxval(abs(sig))*sgpa

! re-expanded FCs at the relaxed state (kmat, pip already there): Psi' = Psi + Chi (S_l - S_tau)
 do n=1,nfct(3)
    psip(n)=v3(n)
 enddo
 do n=1,nfct(4)
    if(p3of4(n).eq.0) cycle
    dl=dvec(iat4(1,n),iat4(4,n))
    psip(p3of4(n))=psip(p3of4(n))+v4(n)*dl(ix4(4,n))
 enddo
 kdrop=0
 do ip=nu+1,np
    kdrop=max(kdrop,maxval(abs(kmat(:,:,ip))))
 enddo
 write(ulog,'(a,es12.4)')' RELAX_REFERENCE: largest effective FC2 on pairs beyond the FC2 range (dropped, eV/A^2) =',kdrop

! back to the irreducible FCs of the same groups (least squares in each group)
 resmax=0; fmax=0
 if(include_fc(1).eq.1) call project(1)
 if(include_fc(2).eq.1) call project(2)
 if(include_fc(3).eq.1) call project(3)
 write(ulog,'(a,3(1x,es12.4))')' RELAX_REFERENCE: max |full FC - symmetrized| ranks 1-3 =',resmax
 write(ulog,'(a,3(1x,es12.4))')' RELAX_REFERENCE: max |full FC| ranks 1-3                =',fmax

! geometry: x -> M x + u0(tau), reciprocal vectors -> M^-T g
 em=eta
 do i=1,3
    em(i,i)=em(i,i)+1
 enddo
 call xmatinv(3,em,minv,ier)
 detm=em(1,1)*(em(2,2)*em(3,3)-em(2,3)*em(3,2))-em(1,2)*(em(2,1)*em(3,3)-em(2,3)*em(3,1)) &
 &   +em(1,3)*(em(2,1)*em(3,2)-em(2,2)*em(3,1))
 do i=1,natoms
    atompos(:,i)=matmul(em,x(:,i))+u0(:,iatomcell0(i))
 enddo
! space-group images of the reference atoms (old positions, old fractional translations), for the
! fractional translations of the relaxed structure below
 allocate(imgsg(nprim,isgopcount))
 do isg=1,isgopcount
    do tau=1,nprim
       imgsg(tau,isg)=0
       do taup=1,nprim
          if(atom0(taup)%at_type.ne.atom0(tau)%at_type) cycle
          dl=matmul(cart_to_prim,matmul(op_matrix(:,:,isgop(isg)),atompos0(:,tau))+sgfract(:,isg)-atompos0(:,taup))
          if(sqrt(sum(matmul(prim_to_cart,dl-nint(dl))**2)).le.tol_sym) then
             imgsg(tau,isg)=taup; exit
          endif
       enddo
    enddo
 enddo
 do tau=1,nprim
    atompos0(:,tau)=matmul(em,atompos0(:,tau))+u0(:,tau)
    atom0(tau)%equilibrium_pos=atompos0(:,tau)
 enddo
 conv_to_cart=matmul(em,conv_to_cart); prim_to_cart=matmul(em,prim_to_cart)
 cart_to_conv=matmul(cart_to_conv,minv); cart_to_prim=matmul(cart_to_prim,minv)
 cart_to_redpr=matmul(cart_to_redpr,minv); cart_to_redsc=matmul(cart_to_redsc,minv)
 r1conv=matmul(em,v2a(r1conv)); r2conv=matmul(em,v2a(r2conv)); r3conv=matmul(em,v2a(r3conv))
 r01=matmul(em,v2a(r01)); r02=matmul(em,v2a(r02)); r03=matmul(em,v2a(r03))
 rs1=matmul(em,v2a(rs1)); rs2=matmul(em,v2a(rs2)); rs3=matmul(em,v2a(rs3))
 g1conv=matmul(transpose(minv),v2a(g1conv)); g2conv=matmul(transpose(minv),v2a(g2conv))
 g3conv=matmul(transpose(minv),v2a(g3conv))
 g01=matmul(transpose(minv),v2a(g01)); g02=matmul(transpose(minv),v2a(g02)); g03=matmul(transpose(minv),v2a(g03))
 gs1=matmul(transpose(minv),v2a(gs1)); gs2=matmul(transpose(minv),v2a(gs2)); gs3=matmul(transpose(minv),v2a(gs3))
 rws26=matmul(em,rws26); r0ws26=matmul(em,r0ws26)
 gws26=matmul(transpose(minv),gws26); g0ws26=matmul(transpose(minv),g0ws26)
 if(allocated(rgrid)) rgrid=matmul(em,rgrid)
 if(allocated(subgrid)) subgrid=matmul(em,subgrid)
 if(allocated(ggrid)) ggrid=matmul(transpose(minv),ggrid)
 volume_r0=volume_r0*detm; volume_r=volume_r*detm
! fractional translations (Cartesian) of the relaxed structure: the strain changes them (GaN, strain zz
! -2.3e-3: the 6_3 translation c/2 off by 6e-3 A > tol_sym, so gamma_mode_symmetry found no images; user's
! report, 2026-11-04). t = mean over the atoms of x'(image) - S x', each reduced to the nearest lattice vector.
 sgmax=0
 do isg=1,isgopcount
    tnew=0; ntn=0
    do tau=1,nprim
       if(imgsg(tau,isg).eq.0) cycle
       dl=atompos0(:,imgsg(tau,isg))-matmul(op_matrix(:,:,isgop(isg)),atompos0(:,tau))
       if(ntn.eq.0) tref=dl
       dd=matmul(cart_to_prim,dl-tref)
       tnew=tnew+dl-matmul(prim_to_cart,dble(nint(dd)))
       ntn=ntn+1
    enddo
    if(ntn.eq.0) cycle
    tnew=tnew/ntn
    dd=matmul(cart_to_prim,tnew-matmul(em,sgfract(:,isg)))
    tnew=tnew-matmul(prim_to_cart,dble(nint(dd)))      ! the lattice image nearest the strained old one
    sgmax=max(sgmax,maxval(abs(tnew-sgfract(:,isg))))
    sgfract(:,isg)=tnew
 enddo
 deallocate(imgsg)
 write(ulog,'(a,es12.4)')' RELAX_REFERENCE: fractional translations updated to the relaxed cell; max change (A) =',sgmax
 volume_g0=volume_g0/detm; volume_g=volume_g/detm
 box=[length(r01),length(r02),length(r03)]; boxg=[length(g01),length(g02),length(g03)]
 density=density/detm
 ref_relaxed=.true.

 call write_out(ulog,' RELAX_REFERENCE: relaxed primitive vectors (columns, A) ',prim_to_cart)
 write(ulog,*)'RELAX_REFERENCE: volume of the primitive cell before, after (A^3) =',v0,volume_r0
 do tau=1,nprim
    write(ulog,'(a,i4,3(1x,f14.8))')' RELAX_REFERENCE: tau, relaxed position (A) =',tau,atompos0(:,tau)
 enddo
 call write_poscar_relaxed
 write(*,'(a,6(1x,es10.3))')' RELAX_REFERENCE: done; strain (Voigt) =',e6

 deallocate(phi_c,q_c,mop_c,iat3c,ix3c,psi_c,iat4c,ix4c,chi_c)
 if(allocated(upth)) deallocate(upth)
 if(allocated(cmat)) deallocate(cmat,gc0,gu0h)
 deallocate(iat1,ix1,v1,g1,t1,iat2,ix2,v2,g2,t2,iat3,ix3,v3,g3,t3,iat4,ix4,v4,x,pidx,pt,pj,kmat,kmat0, &
 &          pip,psip,hkey,hval,p3of4,u0,phi_,q_,xc,um,dum3,ustart)

 contains

! Exact-Coulomb route (coulx): relaxation by Broyden on G(X) = dE/dX, X = (6 Voigt strains, u0 of atoms 2..)
 subroutine relax_broyden
 integer it,m,k,l,tau,al
 real(r15), allocatable :: xv(:),gv(:),fv(:),hlu(:,:)
 integer, allocatable :: ihl(:)
 real(r15) ek(3,3),el(3,3),ar(3,3,3,3),cc(3,3,3,3),dd8,fmaxv,smaxv,fnew(3,3),volx,uz(3,nprim),gtmp(3,nprim)
 logical conv
 allocate(xv(nx),gv(nx),fv(nx),hlu(nx,nx),ihl(nx))
! preconditioner: Hessian at the reference (atoms fixed strain block, Q coupling, phi), Coulomb to second order
 coul=.true.; eta=0; u0=0
 call kfill(1d0); call build_aq(x,.true.,f1); ar=a_/v0; call elastic_born_huang(ar,s0/v0,cc)
 hlu=0
 do k=1,6
    call eta_of(unit6(k),ek)
    do l=1,6
       call eta_of(unit6(l),el)
       hlu(k,l)=v0*(contract(cc,ek,el)+sum(s0/v0*(matmul(ek,el)+matmul(el,ek)))/2)
    enddo
    do tau=2,nprim
    do al=1,3
       m=6+3*(tau-2)+al
       hlu(k,m)=sum(q_(al+3*(tau-1),:,:)*ek); hlu(m,k)=hlu(k,m)
    enddo
    enddo
 enddo
 if(nx.gt.6) hlu(7:nx,7:nx)=phi_(4:3*nprim,4:3*nprim)
 coul=.false.
 call ludcmp(hlu,nx,nx,ihl,dd8)
 call broyden_init(nx,200,1d0)
 write(ulog,'(a)')' RELAX_REFERENCE: Broyden, exact Coulomb: iter, max|dE/de_k|/V (GPa), max|dE/du0| (eV/A)'
 xv=0; conv=.false.
 do it=1,199
    call gradfun(xv,gv,fnew,volx)
    fmaxv=0
    if(nx.gt.6) fmaxv=maxval(abs(gv(7:nx)))
    smaxv=maxval(abs(gv(1:6)))/volx*sgpa
    write(ulog,'(a,i4,2(1x,es12.4))')' RELAX_REFERENCE: ',it,smaxv,fmaxv
    if(fmaxv.lt.1d-10 .and. smaxv.lt.1d-8) then
       conv=.true.; exit
    endif
    fv=-gv
    call lubksb(hlu,nx,nx,ihl,fv)
    call bro90(nx,xv,fv,it)
 enddo
 if(.not.conv) then
    call warn3(ulog,'RELAX_REFERENCE: Broyden not converged in 199 steps; max|dE/de|/V (GPa) =',smaxv)
    call warn3(   6,'RELAX_REFERENCE: Broyden not converged in 199 steps; max|dE/de|/V (GPa) =',smaxv)
 endif
 call broyden_free
 uz=0
 if(nx.gt.6) uz(:,2:nprim)=reshape(xv(7:nx),[3,nprim-1])
 call coul_grad(fnew,uz,fv,gtmp)      ! fv reused: Coulomb gradient at the relaxed state
 write(ulog,'(a,6es16.8)')' RELAX_REFERENCE: Coulomb dE/de_k at the reference (eV)          ',gc0(1:6)
 write(ulog,'(a,6es16.8)')' RELAX_REFERENCE: Coulomb dE/de_k at the relaxed state (eV)      ',fv(1:6)
 write(ulog,'(a,6es16.8)')' RELAX_REFERENCE: short-range dE/de_k at the relaxed state (eV)  ',gv(1:6)-fv(1:6)
 ft=fnew
 u0=0
 if(nx.gt.6) u0(:,2:nprim)=reshape(xv(7:nx),[3,nprim-1])
 eta=ft; do m=1,3; eta(m,m)=eta(m,m)-1; enddo
 sig=0                                   ! stress of the relaxed state from dE/de_k (= V sigma at a force-free state)
 do m=1,6
    if(m.le.3) then
       sig(vg(1,m),vg(2,m))=gv(m)/volx
    else
       sig(vg(1,m),vg(2,m))=gv(m)/volx/2; sig(vg(2,m),vg(1,m))=sig(vg(1,m),vg(2,m))
    endif
 enddo
 deallocate(xv,gv,fv,hlu,ihl)
 end subroutine relax_broyden

! Coulomb couplings and gradient at the reference; FC1 and s0 made short-range (v1 - g_C,u(0), s0 - g_C,eta(0))
 subroutine coulx_setup
 integer k,l,n1
 real(r15) uz(3,nprim),gq0(3,3),rest(3,nprim)
 nx=6+3*(nprim-1)
 allocate(cmat(nprim,nprim),gc0(nx),gu0h(3,nprim))
 e2c=ee*1d10/(4*pi*eps0)
 do k=1,nprim
 do l=1,nprim
    cmat(k,l)=sum(atom0(k)%charge*atom0(l)%charge)/3      ! tr(Z_k Z_l^T)/3, as ewaldforce
 enddo
 enddo
 alph=sqrt(pi)*(dble(nprim)/v0**2)**(1d0/6)
 uz=0
 call coul_grad(f1,uz,gc0,gu0h,gq0)
 call ewald_check
 rest=gu0h
 do n1=1,nfct(1)
    v1(n1)=v1(n1)-gu0h(ix1(n1),iat1(n1))
    rest(ix1(n1),iat1(n1))=0
 enddo
 write(ulog,'(a,es12.4)')' RELAX_REFERENCE: Coulomb force at the reference (eV/A), max =',maxval(abs(gu0h))
 if(maxval(abs(rest)).gt.1d-8) call warn3(ulog, &
 &  'RELAX_REFERENCE: Coulomb force on an FC1 component not in the FC1 groups (not subtracted), max =',maxval(abs(rest)))
 call write_out(ulog,' RELAX_REFERENCE: Coulomb stress at the reference (GPa) ',(gq0+transpose(gq0))/2/v0*sgpa)
 s0=s0-(gq0+transpose(gq0))/2
 gm=s0
 end subroutine coulx_setup

! total gradient at X: short-range (Taylor of the fitted FCs) + exact Coulomb change
 subroutine gradfun(xv,gv,fnew,volx)
 real(r15), intent(in) :: xv(:)
 real(r15), intent(out) :: gv(:),fnew(3,3),volx
 real(r15) em1(3,3),fi(3,3),gq(3,3),uu(3,nprim),du(3,nprim),ek(3,3),gcv(nx),gu(3,nprim)
 integer k,tau,al,ie
 call eta_of(xv(1:6),em1); fnew=f1+em1
 call sr_path_to(fnew)                         ! (fpth,upth,spth): short-range force-free state at fnew
 volx=v0*det3(fnew); call xmatinv(3,fnew,fi,ie)
 gq=volx*matmul(spth,transpose(fi))            ! dE_SR/deta there (forces 0)
 uu=0
 if(nx.gt.6) uu(:,2:nprim)=reshape(xv(7:nx),[3,nprim-1])
 du=uu-upth
 eta=em1
 if(maxval(abs(du)).gt.0) then                 ! carry dE/deta from u0* to u0: int Q^T du (Simpson)
    do k=1,3
       u0=upth+tsimp(k)*du
       call kfill(1d0); call build_aq(x,.false.,f1)
       do tau=1,nprim
       do al=1,3
          gq=gq+wsimp(k)*q_(al+3*(tau-1),:,:)*du(al,tau)
       enddo
       enddo
    enddo
 endif
 u0=uu; call kfill(1d0); call pi_prime
 do k=1,6
    call eta_of(unit6(k),ek); gv(k)=sum(gq*ek)
 enddo
 if(nx.gt.6) gv(7:nx)=reshape(pip(:,2:nprim),[3*(nprim-1)])
 call coul_grad(fnew,uu,gcv,gu)
 gv=gv+gcv
 end subroutine gradfun

! moves the short-range force-free state (fpth,upth,spth) to F = fnew along F(t) = (1 + t dstr) fpth (Coulomb off)
 subroutine sr_path_to(fnew)
 real(r15), intent(in) :: fnew(3,3)
 real(r15) fi(3,3),ssl(3,3),fmm(3,3),q1(3,3),q2(3,3),q3(3,3),q4(3,3),tt,jd
 integer kk,ns,ie,m
 call xmatinv(3,fpth,fi,ie)
 dstr=matmul(fnew,fi); do m=1,3; dstr(m,m)=dstr(m,m)-1; enddo
 if(maxval(abs(dstr)).lt.1d-15) return
 ft=fpth; ssl=spth; ustart=upth
 ns=max(nsub,ceiling(maxval(abs(dstr))/2d-3))
 do kk=1,ns
    tt=dble(kk-1)/ns
    call dsdt(tt,ssl,ustart,q1)
    call dsdt(tt+0.5d0/ns,ssl+q1/(2*ns),ustart,q2)
    call dsdt(tt+0.5d0/ns,ssl+q2/(2*ns),ustart,q3)
    call dsdt(tt+1d0/ns,ssl+q3/ns,ustart,q4)
    ssl=ssl+(q1+2*q2+2*q3+q4)/(6*ns)
 enddo
 fmm=dstr; do m=1,3; fmm(m,m)=fmm(m,m)+1; enddo
 jd=det3(fmm)
 spth=matmul(fmm,matmul(ssl,transpose(fmm)))/jd
 spth=(spth+transpose(spth))/2
 fpth=fnew
 u0=ustart; call relax_u(fnew,u0); upth=u0
 end subroutine sr_path_to

! exact Coulomb gradient at positions fmat x + uu, cell fmat a: gcv = (dE/de_k, dE/du0 atoms 2..), gu = dE/du0 all atoms
 subroutine coul_grad(fmat,uu,gcv,gu,gqo)
 real(r15), intent(in) :: fmat(3,3),uu(3,nprim)
 real(r15), intent(out) :: gcv(:),gu(3,nprim)
 real(r15), intent(out), optional :: gqo(3,3)
 real(r15) ec,fc(3,nprim),wc(3,3),fi(3,3),gq(3,3),ek(3,3),v(3)
 integer m,a,ie
 call coul_energy(fmat,uu,alph,ec,fc,wc)
 gu=-fc
 call xmatinv(3,fmat,fi,ie)
 gq=matmul(wc,transpose(fi))            ! dE/deta at fixed u0: W F^-T + sum_i f_i (F^-1 u_i)
 do m=1,nprim
    v=matmul(fi,uu(:,m))
    do a=1,3
       gq(a,:)=gq(a,:)+fc(a,m)*v
    enddo
 enddo
 do m=1,6
    call eta_of(unit6(m),ek); gcv(m)=sum(gq*ek)
 enddo
 if(nx.gt.6) gcv(7:nx)=reshape(gu(:,2:nprim),[3*(nprim-1)])
 if(present(gqo)) gqo=gq
 end subroutine coul_grad

 subroutine coul_energy(fmat,uu,al0,ec,fc,wc)
 real(r15), intent(in) :: fmat(3,3),uu(3,nprim),al0
 real(r15), intent(out) :: ec,fc(3,nprim),wc(3,3)
 real(r15) pos(3,nprim),cell(3,3)
 integer m
 do m=1,nprim
    pos(:,m)=matmul(fmat,x(:,m))+uu(:,m)
 enddo
 cell(:,1)=matmul(fmat,v2a(r01)); cell(:,2)=matmul(fmat,v2a(r02)); cell(:,3)=matmul(fmat,v2a(r03))
 call ewald_cell_efw(nprim,pos,cell,cmat,epsil,al0,ec,fc,wc)
 ec=ec*e2c; fc=fc*e2c; wc=wc*e2c
 end subroutine coul_energy

! checks of ewald_cell at a deformed, displaced state: energy independent of the Ewald parameter; analytic gradient
! = 4-point differences of the energy (h = 1e-4)
 subroutine ewald_check
 real(r15) f2(3,3),ut(3,nprim),g1(nx),g2(nx),gu(3,nprim),ea,eb,ep(4),fc(3,nprim),wc(3,3),ek(3,3),h,fp(3,3),up(3,nprim)
 integer k,j,m
 real(r15), parameter :: st(4)=[-2d0,-1d0,1d0,2d0]
 f2=reshape([1d-2,2d-3,1d-3, 2d-3,-5d-3,3d-3, 1d-3,3d-3,7d-3],[3,3])
 do m=1,3; f2(m,m)=f2(m,m)+1; enddo
 ut=0
 do m=2,nprim
    ut(:,m)=2d-2*[sin(1d0*m),cos(2d0*m),sin(3d0*m+1)]
 enddo
 call coul_energy(f2,ut,alph,ea,fc,wc)
 call coul_energy(f2,ut,1.3d0*alph,eb,fc,wc)
 call coul_grad(f2,ut,g1,gu)
 h=1d-4
 do k=1,nx
    do j=1,4
       fp=f2; up=ut
       if(k.le.6) then
          call eta_of(unit6(k),ek); fp=f2+st(j)*h*ek
       else
          up(mod(k-7,3)+1,(k-7)/3+2)=up(mod(k-7,3)+1,(k-7)/3+2)+st(j)*h
       endif
       call coul_energy(fp,up,alph,ep(j),fc,wc)
    enddo
    g2(k)=(ep(1)-8*ep(2)+8*ep(3)-ep(4))/(12*h)
 enddo
 write(ulog,'(a,3es12.4)')' RELAX_REFERENCE: Ewald check: E (eV), E(1.3 alpha) - E, max|grad|          =', &
 &     ea,eb-ea,maxval(abs(g1))
 write(ulog,'(a,2es12.4)')' RELAX_REFERENCE: Ewald check: max|analytic - numerical| strain, u0 (eV, eV/A) =', &
 &     maxval(abs(g1(1:6)-g2(1:6))),maxval(abs(g1(7:nx)-g2(7:nx)))
 end subroutine ewald_check

! FC1 at the relaxed state for the exact-Coulomb route: Pi' (short-range) + Coulomb change of the forces
 subroutine add_coul_pip
 real(r15) gcv(nx),gu(3,nprim)
 call coul_grad(ft,u0,gcv,gu)
 pip=pip+gu                     ! v1 holds the short-range FC1 (coulx_setup)
 end subroutine add_coul_pip

 subroutine addpair(i1,j1)
 integer, intent(in) :: i1,j1
 if(pidx(i1,j1).ne.0) return
 np=np+1; pidx(i1,j1)=np; pt(np)=i1; pj(np)=j1
 end subroutine addpair

 function key3(ia,ix) result(kk8)
 integer, intent(in) :: ia(3),ix(3)
 integer(8) kk8
 kk8=((int(ia(1)-1,8)*natoms+(ia(2)-1))*natoms+(ia(3)-1))*27+(ix(1)-1)*9+(ix(2)-1)*3+(ix(3)-1)
 end function key3

 subroutine hput(kk8,v)
 integer(8), intent(in) :: kk8
 integer, intent(in) :: v
 integer h
 h=int(modulo(kk8*2654435761_8,int(nh,8)))
 do while(hkey(h).ne.-1 .and. hkey(h).ne.kk8)
    h=mod(h+1,nh)
 enddo
 hkey(h)=kk8; hval(h)=v
 end subroutine hput

 function hget(kk8) result(v)
 integer(8), intent(in) :: kk8
 integer v,h
 v=0
 h=int(modulo(kk8*2654435761_8,int(nh,8)))
 do while(hkey(h).ne.-1)
    if(hkey(h).eq.kk8) then
       v=hval(h); return
    endif
    h=mod(h+1,nh)
 enddo
 end function hget

! d^2E/de_k de_l (eV) at the reference (F = 1, u0 = 0, atoms held fixed), Voigt e (eta_of), as the Newton
! step builds it: V [c(brackets, s0) (E^k,E^l) + s0:(E^k E^l + E^l E^k)/2], and its Coulomb part (the a_c
! brackets of Phi_C). Printed for comparison with the exact model (REFERENCE/exact_model/relax/
! coulomb_hess6.py).
 subroutine hess6_reference
 real(r15) ar(3,3,3,3),cc(3,3,3,3),cc0(3,3,3,3),ek(3,3),el(3,3),hh(6,6),hc(6,6)
 logical csave
 integer k,l
 call kfill(1d0)
 call build_aq(x,.false.,f1); ar=a_/v0; call elastic_born_huang(ar,s0/v0,cc)
 csave=coul; coul=.false.
 call build_aq(x,.false.,f1); ar=a_/v0; call elastic_born_huang(ar,s0/v0,cc0)
 coul=csave
 do k=1,6
    call eta_of(unit6(k),ek)
    do l=1,6
       call eta_of(unit6(l),el)
       hh(k,l)=v0*(contract(cc,ek,el)+sum(s0/v0*(matmul(ek,el)+matmul(el,ek)))/2)
       hc(k,l)=v0*contract(cc-cc0,ek,el)
    enddo
 enddo
 write(ulog,'(a)')' RELAX_REFERENCE: d2E/de_k de_l at the reference, atoms fixed (eV; Voigt xx yy zz yz xz xy), total:'
 do k=1,6
    write(ulog,'(a,6es16.8)')' RELAX_REFERENCE: H6 ',hh(k,:)
 enddo
 if(coul) then
    write(ulog,'(a)')' RELAX_REFERENCE: its Coulomb part (a_c brackets of Phi_C):'
    do k=1,6
       write(ulog,'(a,6es16.8)')' RELAX_REFERENCE: H6C ',hc(k,:)
    enddo
 endif
 end subroutine hess6_reference

 subroutine eta_of(e,et)
 real(r15), intent(in) :: e(6)
 real(r15), intent(out) :: et(3,3)
 integer m
 et=0
 do m=1,6
    et(vg(1,m),vg(2,m))=e(m); et(vg(2,m),vg(1,m))=e(m)
 enddo
 end subroutine eta_of

 integer function ncomb(m)
 integer, intent(in) :: m
 ncomb=1
 if(m.gt.3) ncomb=2
 end function ncomb

 subroutine comb(m,i0,a0,b0)
 integer, intent(in) :: m,i0
 integer, intent(out) :: a0,b0
 if(i0.eq.1) then
    a0=vg(1,m); b0=vg(2,m)
 else
    a0=vg(2,m); b0=vg(1,m)
 endif
 end subroutine comb

! S_j - S_i at the current eta, u0 (scaled by the module-level factor through eta,u0 arguments)
 function dvec(i0,j0) result(dv)
 integer, intent(in) :: i0,j0
 real(r15) dv(3)
 dv=matmul(eta,x(:,j0)-x(:,i0))+u0(:,iatomcell0(j0))-u0(:,iatomcell0(i0))
 end function dvec

! effective FC2 at the state t*(eta,u0): K = Phi + t Psi d + t^2 Chi d d/2
 subroutine kfill(ts)
 real(r15), intent(in) :: ts
 integer m,p
 real(r15) d1(3),d2(3)
 kmat=0
 do m=1,nfct(2)
    p=pidx(iat2(1,m),iat2(2,m))
    kmat(ix2(1,m),ix2(2,m),p)=kmat(ix2(1,m),ix2(2,m),p)+v2(m)
 enddo
 do m=1,nfct(3)
    p=pidx(iat3(1,m),iat3(2,m))
    d1=dvec(iat3(1,m),iat3(3,m))
    kmat(ix3(1,m),ix3(2,m),p)=kmat(ix3(1,m),ix3(2,m),p)+ts*v3(m)*d1(ix3(3,m))
 enddo
 do m=1,nfct(4)
    p=pidx(iat4(1,m),iat4(2,m))
    d1=dvec(iat4(1,m),iat4(3,m)); d2=dvec(iat4(1,m),iat4(4,m))
    kmat(ix4(1,m),ix4(2,m),p)=kmat(ix4(1,m),ix4(2,m),p)+ts*ts*v4(m)*d1(ix4(3,m))*d2(ix4(4,m))/2
 enddo
! Coulomb Psi_C, Chi_C (empty for now; their pairs must then be in pidx, added with addpair)
 do m=1,nfct3c
    p=pidx(iat3c(1,m),iat3c(2,m))
    d1=dvec(iat3c(1,m),iat3c(3,m))
    kmat(ix3c(1,m),ix3c(2,m),p)=kmat(ix3c(1,m),ix3c(2,m),p)+ts*psi_c(m)*d1(ix3c(3,m))
 enddo
 do m=1,nfct4c
    p=pidx(iat4c(1,m),iat4c(2,m))
    d1=dvec(iat4c(1,m),iat4c(3,m)); d2=dvec(iat4c(1,m),iat4c(4,m))
    kmat(ix4c(1,m),ix4c(2,m),p)=kmat(ix4c(1,m),ix4c(2,m),p)+ts*ts*chi_c(m)*d1(ix4c(3,m))*d2(ix4c(4,m))/2
 enddo
 end subroutine kfill

! A(ab,cd) = -1/2 sum K(tau a, j c) r_b r_d ; Q(tau al; a,b) = sum_j K(tau al, j a) r_b ; phi(tau al, tau' be)
! fd: deformation of the positions xx relative to the reference x (for the Coulomb lattice sums)
! Coulomb (Phi_C fixed): r = fd r0 + du, du = u_t' - u_t the internal displacements in xx (u = xx - fd x on
! the primitive atoms; 0 when xx = x). 2026-11-04: the du terms were missing, so the strain moments and
! brackets of Phi_C were those of the unrelaxed atoms (HfO2r with charges: strain off by 1.4e-4).
 subroutine build_aq(xx,wphi,fd)
 real(r15), intent(in) :: xx(:,:),fd(3,3)
 logical, intent(in) :: wphi
 integer p,m1,m2,m3,m4,i1,j1,n1,n2,t1,t2
 real(r15) rr(3),uu(3,nprim),du(3),fm(3)
 a_=0; q_=0
 if(wphi) phi_=0
 do p=1,np
    i1=pt(p); j1=pj(p)
    rr=xx(:,j1)-xx(:,i1)
    do m1=1,3
    do m3=1,3
       if(kmat(m1,m3,p).eq.0) cycle
       do m2=1,3
          q_(m1+3*(i1-1),m3,m2)=q_(m1+3*(i1-1),m3,m2)+kmat(m1,m3,p)*rr(m2)
          do m4=1,3
             a_(m1,m2,m3,m4)=a_(m1,m2,m3,m4)-kmat(m1,m3,p)*rr(m2)*rr(m4)/2
          enddo
       enddo
       if(wphi) phi_(m1+3*(i1-1),m3+3*(iatomcell0(j1)-1))=phi_(m1+3*(i1-1),m3+3*(iatomcell0(j1)-1))+kmat(m1,m3,p)
    enddo
    enddo
 enddo
 if(coul) then   ! Phi_C: r -> fd r in the moments and brackets
    do m1=1,3*nprim
       q_(m1,:,:)=q_(m1,:,:)+matmul(q_c(m1,:,:),transpose(fd))
    enddo
    do m1=1,3; do m2=1,3; do m3=1,3; do m4=1,3
       do n1=1,3; do n2=1,3
          a_(m1,m2,m3,m4)=a_(m1,m2,m3,m4)+fd(m2,n1)*fd(m4,n2)*a_c(m1,n1,m3,n2)
       enddo; enddo
    enddo; enddo; enddo; enddo
    do t1=1,nprim
       uu(:,t1)=xx(:,t1)-matmul(fd,x(:,t1))
    enddo
    do t1=1,nprim
    do t2=1,nprim
       du=uu(:,t2)-uu(:,t1)
       if(maxval(abs(du)).eq.0) cycle
       do m1=1,3
       do m3=1,3
          do m2=1,3
             q_(m1+3*(t1-1),m3,m2)=q_(m1+3*(t1-1),m3,m2)+phi_c(m1+3*(t1-1),m3+3*(t2-1))*du(m2)
          enddo
          fm=matmul(fd,mop_c(t1,t2,m1,m3,:))           ! sum_{j in t2} Phi_C(t1 m1, j m3) (fd r0)
          do m2=1,3; do m4=1,3
             a_(m1,m2,m3,m4)=a_(m1,m2,m3,m4)-(fm(m2)*du(m4)+du(m2)*fm(m4) &
 &                           +phi_c(m1+3*(t1-1),m3+3*(t2-1))*du(m2)*du(m4))/2
          enddo; enddo
       enddo
       enddo
    enddo
    enddo
    if(wphi) phi_=phi_+phi_c
 endif
 end subroutine build_aq

! Pi' = Pi + Phi d + Psi d d/2 + Chi d d d/6
 subroutine pi_prime
 integer m
 real(r15) d1(3),d2(3),d3(3)
 pip=0
 do m=1,nfct(1)
    pip(ix1(m),iat1(m))=pip(ix1(m),iat1(m))+v1(m)
 enddo
 do m=1,nfct(2)
    d1=dvec(iat2(1,m),iat2(2,m))
    pip(ix2(1,m),iat2(1,m))=pip(ix2(1,m),iat2(1,m))+v2(m)*d1(ix2(2,m))
 enddo
 do m=1,nfct(3)
    d1=dvec(iat3(1,m),iat3(2,m)); d2=dvec(iat3(1,m),iat3(3,m))
    pip(ix3(1,m),iat3(1,m))=pip(ix3(1,m),iat3(1,m))+v3(m)*d1(ix3(2,m))*d2(ix3(3,m))/2
 enddo
 do m=1,nfct(4)
    d1=dvec(iat4(1,m),iat4(2,m)); d2=dvec(iat4(1,m),iat4(3,m)); d3=dvec(iat4(1,m),iat4(4,m))
    pip(ix4(1,m),iat4(1,m))=pip(ix4(1,m),iat4(1,m))+v4(m)*d1(ix4(2,m))*d2(ix4(3,m))*d3(ix4(4,m))/6
 enddo
! Coulomb: Phi_C d = phi_c u0 + q_c : eta (lattice sums), + Psi_C d d/2 + Chi_C d d d/6 (empty for now)
 if(coul) then
    do m=1,3*nprim
       pip(mod(m-1,3)+1,(m-1)/3+1)=pip(mod(m-1,3)+1,(m-1)/3+1)+dot_product(phi_c(m,:),reshape(u0,[3*nprim])) &
 &                                  +sum(q_c(m,:,:)*eta)
    enddo
 endif
 do m=1,nfct3c
    d1=dvec(iat3c(1,m),iat3c(2,m)); d2=dvec(iat3c(1,m),iat3c(3,m))
    pip(ix3c(1,m),iat3c(1,m))=pip(ix3c(1,m),iat3c(1,m))+psi_c(m)*d1(ix3c(2,m))*d2(ix3c(3,m))/2
 enddo
 do m=1,nfct4c
    d1=dvec(iat4c(1,m),iat4c(2,m)); d2=dvec(iat4c(1,m),iat4c(3,m)); d3=dvec(iat4c(1,m),iat4c(4,m))
    pip(ix4c(1,m),iat4c(1,m))=pip(ix4c(1,m),iat4c(1,m))+chi_c(m)*d1(ix4c(2,m))*d2(ix4c(3,m))*d3(ix4c(4,m))/6
 enddo
! A net force on the cell cannot be relaxed (it would translate the crystal); with atom 1 held fixed it would
! push equivalent atoms apart (GaN, enforce_inv 0: FC1 summing to -2e-5 eV/A left u0 = 1e-6 A between the two
! Ga, 2026-11-04). FC1 is corrected before relaxing (above); what is left here comes from Phi S + Psi S S/2 + ...
! when the FC2-FC4 translational invariance is not exact. It is not relaxed: sum_tau Pi'_tau = 0.
 pnet=sum(pip,dim=2)/nprim
 do m=1,nprim
    pip(:,m)=pip(:,m)-pnet
 enddo
 end subroutine pi_prime

 function unit6(m) result(e)
 integer, intent(in) :: m
 real(r15) e(6)
 e=0; e(m)=1
 end function unit6

 function contract(cc,e1,e2) result(v)
 real(r15), intent(in) :: cc(3,3,3,3),e1(3,3),e2(3,3)
 real(r15) v
 integer m1,m2,m3,m4
 v=0
 do m1=1,3; do m2=1,3; do m3=1,3; do m4=1,3
    v=v+e1(m1,m2)*cc(m1,m2,m3,m4)*e2(m3,m4)
 enddo; enddo; enddo; enddo
 end function contract

 function det3(m) result(dt)
 real(r15), intent(in) :: m(3,3)
 real(r15) dt
 dt=m(1,1)*(m(2,2)*m(3,3)-m(2,3)*m(3,2))-m(1,2)*(m(2,1)*m(3,3)-m(2,3)*m(3,1))+m(1,3)*(m(2,1)*m(3,2)-m(2,2)*m(3,1))
 end function det3

! Newton step for u0 (atoms 2..): phi_red du = -Pi'  (phi_ and pip of the current state)
 subroutine solve_u(dun)
 real(r15), intent(out) :: dun(3,nprim)
 integer m,nn,ii,jx
 real(r15), allocatable :: hh(:,:),bb(:)
 integer, allocatable :: ix(:)
 real(r15) dd8
 dun=0
 nn=3*(nprim-1)
 if(nn.eq.0) return
 allocate(hh(nn,nn),bb(nn),ix(nn))
 hh=phi_(4:3*nprim,4:3*nprim)
 do m=2,nprim
    bb(3*(m-2)+1:3*(m-1))=-pip(:,m)
 enddo
 call ludcmp(hh,nn,nn,ix,dd8)
 call lubksb(hh,nn,nn,ix,bb)
 do m=2,nprim
    dun(:,m)=bb(3*(m-2)+1:3*(m-1))
 enddo
 deallocate(hh,bb,ix)
 end subroutine solve_u

! u0 relaxed (Pi' = 0 on atoms 2..) at the deformation fmat (positions fmat x + u0)
 subroutine relax_u(fmat,uu)
 real(r15), intent(in) :: fmat(3,3)
 real(r15), intent(inout) :: uu(3,nprim)
 integer m,iu
 real(r15) dun(3,nprim)
 eta=fmat; do m=1,3; eta(m,m)=eta(m,m)-1; enddo
 u0=uu
 do iu=1,50
    call kfill(1d0); call pi_prime
    uconv=0
    if(nprim.gt.1) uconv=maxval(abs(pip(:,2:nprim)))
    if(uconv.lt.1d-12 .or. nprim.eq.1) exit
    call build_aq(x,.true.,f1)
    call solve_u(dun)
    u0=u0+dun
 enddo
 if(uconv.gt.1d-9) call warn3(ulog,'RELAX_REFERENCE: internal relaxation not converged, max|FC1| =',uconv)
 uu=u0
 end subroutine relax_u

! relaxed elastic constants c_rel (eV/A^3) at the force-free state (fmat, uu) with Cauchy stress sg (eV/A^3)
 subroutine elastic_at(fmat,uu,sg,cr,vv)
 real(r15), intent(in) :: fmat(3,3),uu(3,nprim),sg(3,3)
 real(r15), intent(out) :: cr(3,3,3,3),vv
 integer m,m1,m2,m3,m4,nn
 real(r15) araw(3,3,3,3),cc(3,3,3,3),qs(3*nprim,3,3),dd8
 real(r15), allocatable :: gam(:,:)
 eta=fmat; do m=1,3; eta(m,m)=eta(m,m)-1; enddo
 u0=uu
 call kfill(1d0)
 do m=1,natoms
    xc(:,m)=matmul(fmat,x(:,m))+u0(:,iatomcell0(m))
 enddo
 vv=v0*det3(fmat)
 call build_aq(xc,.true.,fmat)
 araw=a_/vv
 call elastic_born_huang(araw,sg,cc)
 do m=1,3*nprim
    qs(m,:,:)=(q_(m,:,:)+transpose(q_(m,:,:)))/2
 enddo
 cr=cc
 nn=3*(nprim-1)
 if(nn.gt.0) then
    allocate(gam(nn,nn))
    call xmatinv(nn,phi_(4:3*nprim,4:3*nprim),gam,m)
    do m1=1,3; do m2=1,3; do m3=1,3; do m4=1,3
       cr(m1,m2,m3,m4)=cr(m1,m2,m3,m4)-dot_product(qs(4:3*nprim,m1,m2),matmul(gam,qs(4:3*nprim,m3,m4)))/vv
    enddo; enddo; enddo; enddo
    deallocate(gam)
 endif
 end subroutine elastic_at

! dS/dt along F = 1 + t dstr from the state (ft, ustart); S: PK2 about that state (eV/A^3, per V of the state)
 subroutine dsdt(tt,sk,ust,dsk)
 real(r15), intent(in) :: tt,sk(3,3)
 real(r15), intent(inout) :: ust(3,nprim)
 real(r15), intent(out) :: dsk(3,3)
 real(r15) ff(3,3),fi(3,3),jt,sg(3,3),cr(3,3,3,3),vv,del(3,3),cl(3,3,3,3),fs(3,3),uu(3,nprim)
 integer m,ie,i1,i2,i3,i4,j1,j2,j3,j4
 ff=tt*dstr; do m=1,3; ff(m,m)=ff(m,m)+1; enddo
 call xmatinv(3,ff,fi,ie)
 jt=det3(ff)
 sg=matmul(ff,matmul(sk,transpose(ff)))/jt
 fs=matmul(ff,ft)
 uu=ust
 call relax_u(fs,uu)
 call elastic_at(fs,uu,sg,cr,vv)
! Lagrangian strain rate of F(t) = 1 + t dstr: (dstr^T F + F^T dstr)/2. Until 2026-11-05 (dstr F + F dstr)/2, right only
! for a symmetric dstr (the Newton steps); the Broyden path step F_new F_old^-1 - 1 is not symmetric (HfO2r_q:
! short-range dE/deta off by 0.17 eV per cell).
 del=(matmul(transpose(dstr),ff)+matmul(transpose(ff),dstr))/2
 cl=0
 do i1=1,3; do i2=1,3; do i3=1,3; do i4=1,3
    do j1=1,3; do j2=1,3; do j3=1,3; do j4=1,3
       cl(i1,i2,i3,i4)=cl(i1,i2,i3,i4)+fi(i1,j1)*fi(i2,j2)*fi(i3,j3)*fi(i4,j4)*cr(j1,j2,j3,j4)
    enddo; enddo; enddo; enddo
 enddo; enddo; enddo; enddo
 cl=cl*jt
 do i1=1,3; do i2=1,3
    dsk(i1,i2)=0
    do i3=1,3; do i4=1,3
       dsk(i1,i2)=dsk(i1,i2)+cl(i1,i2,i3,i4)*del(i3,i4)
    enddo; enddo
 enddo; enddo
 end subroutine dsdt

! rotation r of the polar decomposition fmat = r u (iteration r <- (r + r^-T)/2)
 subroutine polar(fmat,r)
 real(r15), intent(in) :: fmat(3,3)
 real(r15), intent(out) :: r(3,3)
 real(r15) ri(3,3)
 integer m,ie
 r=fmat
 do m=1,50
    call xmatinv(3,r,ri,ie)
    ri=(r+transpose(ri))/2
    if(maxval(abs(ri-r)).lt.1d-15) then
       r=ri; exit
    endif
    r=ri
 enddo
 end subroutine polar

! least-squares irreducible FCs of each group from the full re-expanded values; written into fcs
 subroutine project(rk)
 integer, intent(in) :: rk
 integer gg,m,ni,i1,j1,p,nr,tt
 integer, allocatable :: kept(:),ind(:),rows(:),gl(:),tl(:)
 real(r15), allocatable :: mm(:,:),rhs(:),full(:),mt(:,:)
 real(r15) dd8,fit
 select case(rk)
 case(1); allocate(gl(nfct(1)),tl(nfct(1))); gl=g1; tl=t1
 case(2); allocate(gl(nfct(2)),tl(nfct(2))); gl=g2; tl=t2
 case(3); allocate(gl(nfct(3)),tl(nfct(3))); gl=g3; tl=t3
 end select
 allocate(rows(nfct(rk)))
 m=1
 do gg=1,map(rk)%ngr
! stored (in-cell) terms of group gg: they were stored group by group
    nr=0
    do while(m.le.nfct(rk))
       if(gl(m).ne.gg) exit
       nr=nr+1; rows(nr)=m; m=m+1
    enddo
    if(nr.eq.0) cycle
    allocate(kept(map(rk)%ntind(gg)))
    ni=0
    do i1=1,map(rk)%ntind(gg)
       if(rk.eq.2) then
          if(map(2)%keep(counteri(2,gg,i1)).eq.0) cycle
       endif
       ni=ni+1; kept(ni)=i1
    enddo
    if(ni.eq.0) then
       deallocate(kept); cycle
    endif
    allocate(mm(ni,ni),rhs(ni),ind(ni),full(nr),mt(nr,ni))
    do i1=1,nr
       n=rows(i1); tt=tl(n)
       select case(rk)
       case(1); full(i1)=pip(ix1(n),iat1(n))
       case(2); p=pidx(iat2(1,n),iat2(2,n)); full(i1)=kmat(ix2(1,n),ix2(2,n),p)
       case(3); full(i1)=psip(n)
       end select
       do j1=1,ni
          mt(i1,j1)=map(rk)%gr(gg)%mat(tt,kept(j1))
       enddo
    enddo
    mm=matmul(transpose(mt),mt); rhs=matmul(transpose(mt),full)
    call ludcmp(mm,ni,ni,ind,dd8)
    call lubksb(mm,ni,ni,ind,rhs)
    do i1=1,nr
       fit=dot_product(mt(i1,:),rhs)
       resmax(rk)=max(resmax(rk),abs(fit-full(i1)))
       fmax(rk)=max(fmax(rk),abs(full(i1)))
    enddo
    do i1=1,ni
       fcs(ofs(rk)+current(rk,gg,kept(i1)))=rhs(i1)
    enddo
    deallocate(kept,mm,rhs,ind,full,mt)
 enddo
 deallocate(rows,gl,tl)
 end subroutine project

 subroutine write_poscar_relaxed
 integer ut,ty,m,cnt(natom_type)
 ut=371
 open(ut,file='POSCAR_relaxed',status='unknown')
 write(ut,'(a)')'reference relaxed by FOCEX with the fitted FCs (RELAX_REFERENCE)'
 write(ut,'(a)')'1.0'
 do m=1,3
    write(ut,'(3(1x,f20.12))')prim_to_cart(:,m)
 enddo
 cnt=0
 do m=1,nprim
    cnt(atom0(m)%at_type)=cnt(atom0(m)%at_type)+1
 enddo
 write(ut,'(99(1x,a))')(trim(atname(ty)),ty=1,natom_type)
 write(ut,'(99(1x,i5))')cnt
 write(ut,'(a)')'Cartesian'
 do ty=1,natom_type
 do m=1,nprim
    if(atom0(m)%at_type.eq.ty) write(ut,'(3(1x,f20.12))')atompos0(:,m)
 enddo
 enddo
 close(ut)
 write(ulog,*)'RELAX_REFERENCE: relaxed structure written to POSCAR_relaxed'
 end subroutine write_poscar_relaxed

 end subroutine relax_reference
