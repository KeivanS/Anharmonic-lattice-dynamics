!==========================================================
 subroutine jdos_sy(ndn,nk,kp,eig,temperature)
! Calculate iself except for matrix element contribution
! jdos1p = sum(1,2) (1+n1+n2)*eta/((w1+w2+w)^2-eta^2)
! jdos1n = sum(1,2) (1+n1+n2)*eta/((w1+w2-w)^2-eta^2)
! jdos2p = sum(1,2) (n2-n1)*eta/((w1-w2+w)^2-eta^2)
! jdos2n = sum(1,2) (n2-n1)*eta/((w1-w2-w)^2-eta^2)
! jdos_tot = jdos1p+jdos1n+jdos2p+jdos2n

 use io2
 use om_dos
 use kpoints
 use eigen

 implicit none

 integer i,la,ndn,nk
 real(8) q(3),omq(ndn),eig(ndn,nk),kp(3,nk)
 real(8) temperature, jdos1p(ndn),jdos1n(ndn),jdos2p(ndn),jdos2n(ndn),jdos_tot(ndn)

 open(udos,file='lifetime_dos_occ.dat ',status='unknown')
 write(udos,*) 'dk, q, omq, jdos1p, jdos1n, jdos2p, jdos2n, jdos_tot(la)'
! temperature=300d0

 do i=1,nk
       q(:)=kp(:,i)
       do la=1,ndn
          omq(la) = eig(la,i)
       enddo
       call lifetime_dos_sy(temperature,q,omq,jdos1p,jdos1n,jdos2p,jdos2n,jdos_tot)
       write(udos,3) q,(omq(la),jdos1n(la),jdos2n(la),jdos_tot(la),la=1,ndyn)
       write(ulog,*) 'i, nk', i, nk
 enddo
 close(udos)

3 format(99(1x,g11.4))

 end subroutine jdos_sy

!==========================================================
 subroutine lifetime_dos_sy(temperature,q,omq,jdos1p,jdos1n,jdos2p,jdos2n,jdos_tot)
! use Lorentzian broadening to calculate the delta function inside the self energy term
! for given q and mode index la, calculates sum_1,2
 use io2
 use params
 use om_dos
 use constants
 use eigen
 use kpoints
 implicit none
 integer i,j,k,l,nq,la,i1,j1,k1,nq1,ip,jp,kp,im,jm,km,nqp,nqm,la1,la2,inside
 real(8) q(3),q1(3),q2(3),jdos1p(ndyn),jdos1n(ndyn),jdos2p(ndyn),jdos2n(ndyn),jdos_tot(ndyn)
 real(8) om1,om2,omq(ndyn),delta_l,evl2(ndyn),vg(3,ndyn)
 complex(8) oc2,omz,evc2(ndyn,ndyn)

 real(8) temperature, n_be1, n_be2


 jdos1p=0; jdos1n=0; jdos2p=0; jdos2n=0; jdos_tot=0
 do l=1,nkc
    q1(:)=kpc(:,l)
    call get_k_info(q1,NC,nq1,i1,j1,k1,g1,g2,g3,inside)
    if(l.ne.nq1) then
      write(ulog,*)'n1,nq1,inside=',l,nq1,inside
      write(ulog,*)'q1=',q1
      write(ulog,*)'ijk1=',i1,j1,k1
      stop
    endif
    q2=-q-q1          ! enforce momentum conservation
    call get_freq(q2,ndyn,vg,evl2,evc2)

    do la1=1,ndyn
    do la2=1,ndyn
    do la =1,ndyn

       om1 = eigenval(la1,l)
       om2 = evl2(la2)

    if(temperature.gt.0) then

       call BS_distribution_sy(temperature,om1,n_be1)
       call BS_distribution_sy(temperature,om2,n_be2)

       !write(debug,8) 'temperature,om1,n_be1,om2,n_be2', temperature, om1, n_be1, om2, n_be2
       !write(*,8) 'temperature,om1,n_be1,om2,n_be2', temperature, om1, n_be1, om2, n_be2

       jdos1p(la)=jdos1p(la)+(1.0d0+n_be1+n_be2)*delta_l(om1+om2+omq(la),etaz)
       jdos1n(la)=jdos1n(la)+(1.0d0+n_be1+n_be2)*delta_l(om1+om2-omq(la),etaz)
       jdos2p(la)=jdos2p(la)+(n_be2-n_be1)*delta_l(om1-om2+omq(la),etaz)
       jdos2n(la)=jdos2n(la)+(n_be2-n_be1)*delta_l(om1-om2-omq(la),etaz)

    else ! for negative temperature, do not include BE factors

       jdos1p(la)=jdos1p(la)+delta_l(om1+om2+omq(la),etaz)/nkc
       jdos1n(la)=jdos1n(la)+delta_l(om1+om2-omq(la),etaz)/nkc
       jdos2p(la)=jdos2p(la)+delta_l(om1-om2+omq(la),etaz)/nkc
       jdos2n(la)=jdos2n(la)+delta_l(om1-om2-omq(la),etaz)/nkc

    endif

       jdos_tot(la)=jdos1p(la)+jdos1n(la)+jdos2p(la)+jdos2n(la)

       !write(debug,8) 'jdos1p, 1n, 2p, 2n, tot', jdos1p(la), jdos1n(la), jdos2p(la), jdos2n(la), jdos_tot(la)
       !write(*,8) 'jdos1p, 1n, 2p, 2n, tot', jdos1p(la), jdos1n(la), jdos2p(la), jdos2n(la), jdos_tot(la)

    enddo
    enddo
    enddo
    !write(*,*) 'l, nkc', l, nkc

 enddo
10 format(99(1x,es10.4))
! write(*,8)'LIFETIME_DOS: q,omqdosp,dosm=',q,omq,dosp,dosm
8 format(a,99(1x,g10.3))
 end subroutine lifetime_dos_sy


!============================================================
 subroutine BS_distribution_sy(T,omega,n_be)
! calculate Bose-Einstein distribution for a certain temperature and phonon frequency
! T is temperature in K, omega is phonon frequency in cm^-1, n is distribution function
 use constants

 implicit none
 real(8) T, omega, n_be, omega_Hz

 omega_Hz=omega*c_light*1.0d2        ! convert from cm^-1 to Hz

 n_be=1d0/(EXP(h_plank*omega_Hz/(k_b*T))-1)

 end subroutine BS_distribution_sy

!============================================================
 subroutine matrix_kpt_sy
! print out k2 points, omega, sum(1,3) V33(k1,la1,k2,la2,k3,la3)^2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer j,nq,jq,la,jk1,jk2,jk3,nqq2,nk1,inside,q2,l2
 real(8) omq2(ndyn),sumz(ndyn)

 open(udos,file='v33square_k.dat', status='unknown')
 write(udos,*) 'i, q, omq, v33square'

 do q2=1,nibz
 sumz=0
 do l2=1,ndyn
    call get_k_info(kibz(:,q2),nc,nqq2,jk1,jk2,jk3,g1,g2,g3,inside)
    omq2(l2) = eigenval(l2,nqq2)
    !write(debug,*) 'start v3loop...q2,l2,nv3', q2, l2, nv3
    v3loop: do j=1,nv3
      !write(debug,*) 'j, la2(j), l2, nq2(j), nqq2', j, la2(j), l2, nq2(j), nqq2
      if (la2(j) .ne. l2) cycle v3loop
      if (nq2(j) .ne. nqq2) cycle v3loop
      sumz(l2) = sumz(l2) + (v33(j)*conjg(v33(j)))
      !write(debug,*) 'sumz, v33(j)', sumz, v33(j)
    enddo v3loop

    !write(debug,*) 'mapinv(q2),kibz(:,q2), wibz(q2)', mapinv(q2), kibz(:,q2), wibz(q2)
    sumz(l2)=sumz(l2)*wibz(q2)

    !write(ulog,*) 'matrix_kpt_sy, q2, nibz', q2, nibz
    !write(debug,*) 'kibz(:,q2), wibz(q2)', kibz(:,q2), wibz(q2)

  enddo
  write(udos,8) mapinv(q2), kibz(:,q2), (omq2(l2),sumz(l2),l2=1,ndyn)
  enddo

 close(udos)

8 format(1x,9(1x,g10.4))

 end subroutine matrix_kpt_sy


!===========================================================
 subroutine calculate_w3_ibz_split_sq(ibz_subset,nv3_split,ndn,nk,kp)
! This is a modified calculate_w3_ibz by sy.
! It reads ibz_subset (ibz_subset(1)-starting index, ibz_subset(2)-end index in nibz)
! This subroutine will not write header in output file. This should be done by another code for merging v33.dat

! Below is for the original calculate_w3_ibz
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units for V3 are in cm^-1
! this works for two k's (k2,k3) in the prim cell on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! this subroutine calculates V3(k1,la1,k2,la2,k3,la3)  with k1=-k2-k3+G  and all ki in the primitive G-cell
! with the restriction of k2 in the IRREDUCIBLE BZ defined by ki(3,ni) included in kp(3,nk)
!
 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ndn,l1,l2,l3,n2,n1,j,k,nibzsub
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,indx
 integer ibz_subset(2),nv3_split
 real(8) q1(3),q2(3),q3(3),eival(ndn,nk) !,sh(3)
 real(8) kp(3,nk)  !ki(3,ni),
 complex(8) eivec(ndn,ndn,nk),xx
 real cputim
 character(4) cibz1,cibz2
 character(6) cnk

! v33 = cmplx(0d0,0d0)
 v33sq=0d0
 indx=0
! sh=(-0.5d0)*(g1+g2+g3) ! this is the shift applied to get kibz
 write(ulog,*)'entered v3 subroutine...'
 write(ulog,*)'V3: nc(123)=',nc
 write(ulog,*)'starting ibz point=',ibz_subset(1)
 write(ulog,*)'ending   ibz point=', ibz_subset(2)
 nibzsub=ibz_subset(2)-ibz_subset(1)+1
 write(ulog,*)'number of ibz subset points=', nibzsub
 write(cibz1,'(i4.4)') ibz_subset(1)
 write(cibz2,'(i4.4)') ibz_subset(2)
 write(cnk,'(i6.6)') nk

   open(uv3,file='v33-'//cnk//'-'//cibz1//'-'//cibz2//'.dat',status='unknown') !,FORM='UNFORMATTED')
   write(uv3,*)nv3_split,ndn,ibz_subset(1),ibz_subset(2)

!************************************************
! First argument q1 needs to be only in the IFBZ
!************************************************
 loop2: do n1=ibz_subset(1),ibz_subset(2)       ! calculate v33 only for subset of k points
    q1 = kibz(:,n1) !kp(:,mapinv(n2))  ! is also ki(:,n2)  ! should be
    !write(debug,*) 'q2=' q2
    call get_k_info_cent(q1,NC,nk1,i2,j2,k2,g1,g2,g3,inside)
     if(mapinv(n1).ne.nk1) then
       write(ulog,*)'n1,mapinv(n1),nk1,inside=',n1,mapinv(n1),nk1,inside
       write(ulog,*)'q1=',q1
       write(ulog,*)'ijk=',i2,j2,k2
       stop
     endif

 call cpu_time(cputim)  !date_and_time(time=tim)
 write(*,1) 'W3_IBZ_SPLIT: nibz,nk1,q...',n1,nk1,q1
 write(ulog,1) 'entering nibz,nk1,q...',n1,nk1,q1
 write(ulog,'(a,f12.4)')'  TIME IS ',cputim


 loop3: do n2=1 ,nk
    q2 = kp(:,n2)   ! third argument is in the whole FBZ coarse mesh
    !write(debug,*) 'q3=', q3
    call get_k_info_cent(q2-shft,NC,nk2,i3,j3,k3,g1,g2,g3,inside)
    if(n2.ne.nk2) then
      write(ulog,*)'n2,nk2,inside=',n2,nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i3,j3,k3
      stop
    endif

    q3 = -q2-q1

    call get_k_info_cent(q3-shft,NC,nk3,i1,j1,k1,g1,g2,g3,inside)

 ! write(ulog,2)'nk123=',nk1,nk2,nk3
  !write(ulog,*)'q2,q3,q1',q2,q3,q1
 do l1=1 ,ndn
 do l2=1 ,ndn
 do l3=1 ,ndn
    call matrix_elt(q1,q2,l1,l2,l3,xx,inside)
    indx=indx+1
!    v33(indx)=xx
    v33sq(indx)=xx*conjg(xx)

!   write(uv3,5) n1,n2,l1,l2,l3,v33sq(indx)

     nq1(indx)= (nk1-1)*ndn+l1
     nq2(indx)= (nk2-1)*ndn+l2
     nq3(indx)= (nk3-1)*ndn+l3
     write(uv3,5) nq1(indx),nq2(indx),nq3(indx),v33sq(indx)
! reverse mapping is: la=mod(nq,ndn) ; nk=1+(nq-la)/ndn
!   nq2(indx)= nk2
!   nq3(indx)= nk3
!   la1(indx)= l1
!   la2(indx)= l2
!   la3(indx)= l3
    !write(debug,*),'l1,l3,l2,indx,v33',l1,l3,l2,indx,v33(indx)
!   if (mod(indx,1000).eq.1) write(*,*) 'indx,v33=',indx,xx
 enddo
 enddo
 enddo
 enddo loop3
 enddo loop2

 nv3 = indx
 write(ulog,*)' V33: total size of this array is =',nv3

! if( writev3.eq.1) then
! if the v33 are calculated, they are automatically written into a file
   !write(uv3,*) nv3_split
!   do j=1 ,nv3_split
!   enddo
   close(uv3)
! endif

1 format(a,2(1x,i5),9(1x,g11.5))
2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))

 end subroutine calculate_w3_ibz_split_sq
!-------------------------------------------------------------------
 subroutine read_ksubset (ksubset,nib)
! this subroutine reads ksubset.inp
 use io2
 implicit none
 integer ksubset_f,tmp,nib
 integer ksubset(2)

 ksubset_f=9998
 ksubset=0
 open(ksubset_f,file='ksubset.inp', status='old')
 read(ksubset_f,*) ksubset
 close(ksubset_f)
 if (ksubset(2).lt.ksubset(1)) then
     write(ulog,*) '!!!!!!!!!!!WARNING!!!!!!!!KSUBSET=', ksubset(1), ksubset(2)
     write(ulog,*) 'ckeck your ksubset input'
     tmp=ksubset(1)
     ksubset(1)=ksubset(2)
     ksubset(2)=tmp
 endif
 if (ksubset(2) .gt. nib) then
     write(ulog,*) '!!!!!!!!!!!WARNING!!!!!!!!KSUBSET=', ksubset(1), ksubset(2)
     write(ulog,*) 'ckeck your ksubset input; k(2)>nibz; setting k(2)=nibz'
     ksubset(2)=nib
 endif
 end subroutine read_ksubset

!===========================================================
 subroutine mode_thermal_conductivity_sy(temperature,lambda,nk,kpc,wk,tau_inv,tau_inv_n,tau_inv_u,veloc,omg,temp,kappa_q)
! temp is in cm^-1, veloc in c_light, omg in ev/ang^2/uma, tau_inv in cm_1
! input are the freqs, velocs , relaxation times and weights for each band
! output is the thermal conductivity for that band
 use io2
 use constants
 use params
 use lattice
 implicit none
 integer, intent(in) ::  nk
 integer i,al,be
 real(8),intent(in):: tau_inv(nk),tau_inv_n(nk),tau_inv_u(nk),veloc(3,nk),omg(nk),wk(nk),temp,kpc(3,nk)
 real(8), intent(out) :: kappa_q(3,3)
 real(8) cv,x,tau,nbe,nbx
 real(8) frequency, kappa_temp(3,3),temperature,mfp, kappa_tot
 integer ii,jj, lambda

  kappa_q=0 ; cv=0
  do i=1,nk
     frequency=sqrt(omg(i))*cnst
     x=sqrt(omg(i))*cnst/temp   ! temp is in 1/cm
     if (x.gt.40) then
       cv = 0
     else
       nbx=nbe(x,1d0,classical)
       cv=nbx*(nbx+1)*x*x  !  x*x/4/sinh(x/2)/sinh(x/2)
     endif
     if (tau_inv(i) .eq. 0) then
       tau=1d9
     else
       tau=1/tau_inv(i)
     endif
     if(tau.lt.1d7) then  ! exclude gamma=0
!      kappa_q=kappa_q+cv*wk(i)*tau*sum(veloc(:,i)*veloc(:,i))/3
        do al=1,3
        do be=1,3
            kappa_temp(al,be)=cv*wk(i)*tau*veloc(al,i)*veloc(be,i)*k_b/100*c_light/volume_r*1d30
            kappa_q(al,be)=kappa_q(al,be)+cv*wk(i)*tau*veloc(al,i)*veloc(be,i)

        enddo
        enddo

       mfp=length(veloc(:,i)) * tau * 1d7    ! nanometer
       kappa_tot=(kappa_temp(1,1)+kappa_temp(2,2)+kappa_temp(3,3))/3

       write(k_kappa,7) temperature, i, kpc(:,i), lambda, frequency, tau_inv(i), tau_inv_n(i), tau_inv_u(i), &
& veloc(:,i), mfp, cv, nbx, (kappa_temp(ii,ii),ii=1,3), ((kappa_temp(ii,jj),jj=ii+1,3),ii=1,2), kappa_tot

     endif
  enddo
  kappa_q = kappa_q *k_b/100*c_light/volume_r*1d30  ! converted to SI (W/mK)
7 format(99(g14.8,5x))

 end subroutine mode_thermal_conductivity_sy


!============================================================
 subroutine function_self_w2_sy(q,la,omega,temp,nself,uself)
! calculates the self-energy on the fly for q-point in the generated kmesh and Lorentzian delta
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1,nk3,ins2,ifbz,l1,l3
 integer, save :: cntr
 real(8) omq,om3,omk,eta,tk,term,k1(3),k2(3),q1(3),q3(3)
 complex(8) omz,self,ocfunc,oc2,s2,xx


 integer k_FBZ,q_ibz
 real(8) v33_square_sum,tempk
 complex(8) delta_tot_sum, delta1_sum, delta2_sum, delta3_sum, delta4_sum
 complex(8) delta_tot, delta1, delta2, delta3, delta4
 v33_square_sum=0; delta_tot_sum=0; delta1_sum=0; delta2_sum=0; delta3_sum=0; delta4_sum=0

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
! write(*,*) 'self_w: tempk=',tk
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! write(*,*) 'self_w: k_info: nq=',nq
 omq = eigenval(la,nq)
 eta= etaz
! write(ulog,5)'SELF_W: last cntr, omega,etaz=',cntr,omega,etaz
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0

 write(ulog,*) 'entered self energy calculation... at q and la', nq, la

 cntr = 0
 do ifbz=1,nkc
    q3=kpc(:,ifbz)
    call get_k_info(q3,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
    if (nk3.ne.ifbz) then
       write(ulog,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       write(*,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       stop
    endif
    k1 = -q-q3
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,ins2)

    do l3=1,ndyn
       om3 = eigenval(l3,nk3)
    do l1=1,ndyn
       omk = eigenval(l1,nk1)
!      call matrix_elt     (q,q3,la,l3,l1,xx,inside)
       call matrix_elt_full(q,q3,k1,omq,om3,omk,eigenvec(:,la,nq),eigenvec(:,l3,nk3),eigenvec(:,l1,nk1),xx)
       call check_inside_bz(k1,g1,g2,g3,inside)
       if(inside.ne.ins2) then
          write(*,*)'SELF_W:ERROR:inside .ne. ins2 ',inside,ins2,k1
          stop
       endif

       ! for raw data output
       call ocfunc_sy(temp,omz,l1,nk1,l3,nk3,delta1,delta2,delta3,delta4,delta_tot)

       v33_square_sum=v33_square_sum+xx*conjg(xx)
       delta1_sum=delta1_sum+delta1
       delta2_sum=delta2_sum+delta2
       delta3_sum=delta3_sum+delta3
       delta4_sum=delta4_sum+delta4
       delta_tot_sum=delta_tot_sum+delta_tot
       ! end for raw data output


       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)

       jk1 = nk1+(l1-1)*nkc
       jk3 = ifbz+(l3-1)*nkc

        self = ocfunc(temp,omz,l1,nk1,l3,nk3)
       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif
    enddo
    enddo

 enddo
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself


! Write Raw data.....
! this is to make sure k and all its stars have the same self energy
 tempk=temp*(100*h_plank*c_light)/k_b
! do k_FBZ=1, nkc
!   write(*,*) 'k_FBZ, mapibz(k_FBZ), nq', k_FBZ, mapibz(k_FBZ), nq
!   if(mapibz(k_FBZ) .eq. q_ibz) then   ! IBZ -> FBZ
!     write(*,*) 'matched--------------------------------------------------------'
!     write(self_detail,8) tempk, k_FBZ, kpc(:,k_FBZ), la, omega, &  ! ki, kpoint, lambda, omega
!& v33_square_sum, real(delta1_sum), aimag(delta1_sum), real(delta2_sum), aimag(delta2_sum) &
!& , real(delta3_sum), aimag(delta3_sum), real(delta4_sum), aimag(delta4_sum), real(delta_tot_sum), aimag(delta_tot_sum)
!   endif
! enddo

   write(self_detail,8) tempk, q, la, omega, &
& v33_square_sum, real(delta1_sum), aimag(delta1_sum), real(delta2_sum), aimag(delta2_sum) &
& , real(delta3_sum), aimag(delta3_sum), real(delta4_sum), aimag(delta4_sum), real(delta_tot_sum), aimag(delta_tot_sum)


! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
6 format(i6,6(1x,f9.4),2x,3(1x,i2))
7 format(a,4i3,99(1x,f7.3))
8 format(99(5x,g16.8))

 end subroutine function_self_w2_sy


!============================================================
 subroutine function_self_w3_sy(nq, q,la,omega,temp)
! this one is for arbitrary q
! uses lorentzian delta with 4 terms, and calculates both real and imaginary parts
! In this subroutine, which calculates like self_w on the fly, the second momentum
! is restricted to be on the kpc mesh, but the first and third momenta can be anywhere
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la, nq
 real(8), intent(in) :: q(3),temp,omega
 complex(8) nself,uself
 integer inside,ik,l3,l2
! integer, save :: cntr
 real(8) omq,eta,term,k3(3),k2(3),eivq(ndyn),eiv2(ndyn),eiv3(ndyn),vg(3,ndyn)
 real(8) etacut,arg,delta_l,om3,om2,nb3,nb2,nbe,ires,rres
 complex(8) omz,self,xx,evq(ndyn,ndyn),ev2(ndyn,ndyn),ev3(ndyn,ndyn)

 real(8) v33_square_sum, tempk
 real(8) delta1_sum(2), delta2_sum(2), delta3_sum(2), delta4_sum(2), delta_tot(2)

 tempk = temp/k_b*(100*h_plank*c_light)


! be careful band sorting should be avoided in this case; otherwise need to
! find the nearest kpoint in kpc, and sort eivq according to it
write(*,*) 'la, nq, q, omega, tempk',la,nq,q,omega,tempk

 call get_freq(q,ndyn,vg,eivq,evq)
 omq=eivq(la)

write(*,*) 'omega, omq', omega, omq

 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0
 v33_square_sum=0; delta1_sum=0; delta2_sum=0; delta3_sum=0; delta4_sum=0; delta_tot=0

 etacut = eta*3000 ! error due to non-zero eta is 1/9E6

 do ik=1,nkc
    k2=kpc(:,ik)
    ev2=eigenvec(:,:,ik)
    k3 = -q-k2
    call check_inside_bz(k3,g1,g2,g3,inside)

    call get_freq(k3,ndyn,vg,eiv3,ev3)

    do l2=1,ndyn
       om2=eigenval(l2,ik)
       nb2=nbe(om2,temp,classical)
    do l3=1,ndyn

! first check see if the energy is conserved
       om3=eiv3(l3)
       nb3=nbe(om3,temp,classical)
       ires=0; rres=0

       arg=(omega-om3-om2)        ! omega2+omega3-omega
       if (abs(arg) .lt. etacut) then
          ires=(nb3+nb2+1)*delta_l(arg,eta)
          rres=(nb3+nb2+1)*arg/(arg*arg+eta*eta)
          delta2_sum(1)=delta2_sum(1)+(nb3+nb2+1)*arg/(arg*arg+eta*eta)  ! real part
          delta2_sum(2)=delta2_sum(2)+(nb3+nb2+1)*delta_l(arg,eta)   ! imaginary part
       endif

       arg=(omega+om3+om2)
       if (abs(arg) .lt. etacut) then
          ires=ires-(nb3+nb2+1)*delta_l(arg,eta)   !! Why ires-?
          rres=rres-(nb3+nb2+1)*arg/(arg*arg+eta*eta)
          delta1_sum(1)=delta1_sum(1)+(nb3+nb2+1)*arg/(arg*arg+eta*eta)
          delta1_sum(2)=delta1_sum(2)+(nb3+nb2+1)*delta_l(arg,eta)
       endif

       arg=(omega-om3+om2)      ! omega2-omega3+omega
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb2-nb3)*delta_l(arg,eta)
          rres=rres+(nb2-nb3)*arg/(arg*arg+eta*eta)
          delta3_sum(1)=delta3_sum(1)+(nb2-nb3)*arg/(arg*arg+eta*eta)
          delta3_sum(2)=delta3_sum(2)+(nb2-nb3)*delta_l(arg,eta)
       endif

       arg=(omega-om2+om3)      ! omega2-omega3-omega
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb3-nb2)*delta_l(arg,eta)
          rres=rres+(nb3-nb2)*arg/(arg*arg+eta*eta)
          delta4_sum(1)=delta4_sum(1)+(nb3-nb2)*arg/(arg*arg+eta*eta)
          delta4_sum(2)=delta4_sum(2)+(nb2-nb3)*delta_l(arg,eta)
       endif

       self=cmplx(rres,ires*pi)
 !     if (abs(self).lt.1d-6) cycle loop1   ! self is in cm

       call matrix_elt_full(q,k2,k3,omq,om2,om3,evq(:,la),ev2(:,l2),ev3(:,l3),xx)

       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)
       v33_square_sum=v33_square_sum+term

       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif

 enddo
 enddo
 enddo
 nself = nself /nkc *pi
 uself = uself /nkc *pi
 self  = nself + uself


 delta1_sum(1)=delta1_sum(1)/nkc*pi
 delta1_sum(2)=delta1_sum(2)/nkc*pi**2
 delta2_sum(1)=delta2_sum(1)/nkc*pi
 delta2_sum(2)=delta2_sum(2)/nkc*pi**2
 delta3_sum(1)=delta3_sum(1)/nkc*pi
 delta3_sum(2)=delta3_sum(2)/nkc*pi**2
 delta4_sum(1)=delta4_sum(1)/nkc*pi
 delta4_sum(2)=delta4_sum(2)/nkc*pi**2

 delta_tot=delta1_sum+delta2_sum+delta3_sum+delta4_sum


 write(iself_bs,8) tempk, la, nq, q, v33_square_sum, delta1_sum, delta2_sum, delta3_sum, delta4_sum, delta_tot, &
&  1d0/(2d0*aimag(nself)), 1d0/(2d0*aimag(uself)), 1d0/(2d0*aimag(self))

8 format(99(5x,g16.8))

! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w3_sy






!===========================================================
 subroutine ocfunc_sy(temp,omz,l2,nk2,l3,nk3,delta1,delta2,delta3,delta4,delta_tot)
! for given temperature and complex omz (in 1/cm) calculates the
! occupation-weighted joint DOS stored in res
! j2 and j3 are kpoint indices in the full FBZ
 use eigen
 use constants ! for cnst
 use kpoints   ! for nkc
 use params    ! for classical
 implicit none
 integer , intent(in) :: l2,l3,nk2,nk3
 integer j2,j3
 real(8), intent (in) :: temp
 real(8) om3,om2,nb3,nb2,nbe,resr,resi,et,om,delta_g,delta_l
 complex(8), intent (in) :: omz
 complex(8) delta1, delta2, delta3, delta4, delta_tot

!! nk2 = 1+mod(j2-1,nkc) ; l2 = 1+ (j2-nk2)/nkc
!! nk3 = 1+mod(j3-1,nkc) ; l3 = 1+ (j3-nk3)/nkc
! call nkla(j2,nkc,nk2,l2)
! call nkla(j3,nkc,nk3,l3)
 j2=(l2-1)*nkc+nk2
 j3=(l3-1)*nkc+nk3

 om3 = eigenval(l3,nk3) ; nb3=nbe(om3,temp,classical)
 om2 = eigenval(l2,nk2) ; nb2=nbe(om2,temp,classical)
! om=real(omz) ; et=-aimag(omz)

 delta1=(nb2+nb3+1)*(1d0/(om2+om3+omz))
 delta2=(nb2+nb3+1)*(1d0/(om2+om3-omz))
 delta3=(nb2-nb3  )*(1d0/(om3-om2+omz))
 delta4=(nb2-nb3  )*(1d0/(om3-om2-omz))
 delta_tot = delta1+delta2+delta3+delta4


 end subroutine ocfunc_sy
!===========================================================

subroutine selfenergy2(ksubset,nk,kp,ndn,eigenval,temp)

!use exactBTE
use exactBTE2
use constants
use io2
use lattice
use params
use om_dos
use phi3_sy
use phi3
implicit none

 integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside, ncl
 integer nk1,nk2,nk3, nk1n, nk2n
 integer ksubset(2),nk
 real(8) q1(3), q2(3), q3(3), q1n(3), q2n(3)
 real(8) kp(3,nk)

 integer i,n1,n2,n3,l1,l2,l3,ndn,indx
 real(8) temp     ! temperature in cm^-1
 real(8) eigenval(ndn,nk)
 real(8) omega1, omega2, omega3, V3sq1, V3sq2, col1, col2, nbe, nb2, nb3
 real(8) delta_l
 integer iself2

 iself2=9011

! open(ucol,file='col_matrix.dat',STATUS='unknown',FORM='unformatted')
 open(iself2,file='iself2.dat',STATUS='unknown')

write(*,*) 'entering selfenergy2...'
write(*,*) 'v3_threshold=',v3_threshold

ncl=(ksubset(2)-ksubset(1)+1)*nk*ndn**3
write(iself2,*) ncl

iselfenergy=0
indx=0
write(iself2,*) 'indx, nkF1(indx), laF1(indx), q1, omega, iself(indx)'
loop1 : do n1=ksubset(1),ksubset(2)


    q1=kp(:,n1)
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)

write(*,*) 'q1,nk1=',q1,nk1

    if (n1 .ne. nk1) then
       write(ulog,*) 'n1,nk1,inside=',n1,nk1,inside
       write(ulog,*) 'q1=',q1
       write(ulog,*) 'ijk1=',i1,j1,k1
       stop
    endif

  do l1=1,ndn

indx=indx+1
write(*,*) 'index=',indx
!nkF1(indx)=n1
!laF1(indx)=l1

write(*,*) 'index=',indx

    loop2 : do n2=1,nk

          q2=kp(:,n2)
          call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
          if (n2 .ne. nk2) then
               write(ulog,*) 'n2,nk2,inside=',n2,nk2,inside
               write(ulog,*) 'q2=',q2
               write(ulog,*) 'ijk3=',i2,j2,k2
               stop
           endif

!write(*,*) 'q2,nk2=',q2,nk2

           q3=-q1-q2
           call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)

!write(*,*) 'q3,nk3=',q3,nk3

!           do l3=1, ndn
           do l2=1, ndn
           do l3=1, ndn

                 omega1=sqrt(eigenval(l1,nk1)) * cnst    !linear freq.
                 omega2=sqrt(eigenval(l2,nk2)) * cnst
                 omega3=sqrt(eigenval(l3,nk3)) * cnst

                 q1n=-q1
                 q2n=-q2
                 call get_k_info(q1n,NC,nk1n,i1,j1,k1,g1,g2,g3,inside)
                 call get_k_info(q2n,NC,nk2n,i2,j2,k2,g1,g2,g3,inside)
write(debug,*) '---nk3,q3/nk2,q2/nk1,q1/nk1n,q1n/nk2n,q2n'
write(debug,8) nk3,q3
write(debug,8) nk2,q2
write(debug,8) nk1,q1
write(debug,8) nk1n,q1n
write(debug,8) nk2n,q2n

!                 call find_V3(nk1n,nk2,nk3,l1,l2,l3,V3sq1)   ! return V3sq1=abs(V3)^2* (2pi)^2

!                  V3sq1=v33_md(nk1n,nk2,l1,l2,l3)*conjg(v33_md(nk1n,nk2,l1,l2,l3)) * (2*pi)**2 / nk

!if (abs(v33_md(nk1,nk2,l1,l2,l3)) >= v3_threshold) then
V3sq1=v33sq1(nk1,nk2,l1,l2,l3) * (2*pi)**2 / nk


nb2=nbe(omega2,temp,classical)
nb3=nbe(omega3,temp,classical)
iselfenergy(n1,l1)=iselfenergy(n1,l1) + (pi/2) * V3sq1 * (2*(nb2-nb3)* &
&  delta_l(omega1+omega2-omega3,etaz) + (1+nb2+nb3)*delta_l(omega1-omega2-omega3,etaz)) / &
& (2*pi)


!endif

            enddo
            enddo
!            enddo
     enddo loop2

write(iself2,10) indx, n1, l1, q1, omega1, iselfenergy(n1,l1)

enddo
enddo loop1

close(iself2)

7 format(6(i6),2(2x,g14.8))
8 format(i6,3(2x,g14.8))
10 format(99(1x,g11.4))

end subroutine selfenergy2

!============================================================
 subroutine read_v3_new_unformatted_sy(nk,ndn)

! This is modified version of read_v3_new_unformatted.
! It will use multi-dimensional array for v3

 use phi3
 use phi3_sy
 use io2
 use lattice
 implicit none
 integer j,unt
 real(8) xx,yy
 integer nk,ndn

 integer n1,n2,n3,l1,l2,l3

 unt = 111
 open(unt,file='v33.dat',status='old',form='UNFORMATTED')
!open(unt,file='v33.dat',status='old')

 read(unt,end=99) nv3
! read(unt,*,end=99) nv3
 call allocate_v33_sy(nk,nk,ndn,ndn,ndn)
 write(ulog,*)' OPENING v33.dat, reading ',nv3
 v33_md= cmplx(0d0,0d0)
 do j=1,nv3
!    read(unt,*,end=99) n1,n2,n3,l1,l2,l3,xx,yy
    read(unt,end=99) n1,n2,n3,l1,l2,l3,xx,yy
    if (V33_md(n1,n2,l1,l2,l3) .ne. 0) then
       write(ulog,*) 'V33_md duplicate'
       stop
    endif

    v33_md(n1,n2,l1,l2,l3)=cmplx(xx,yy)
 enddo
 close(unt)
 return

99 write(*,*)'v33 file end reached at line j=',j
   stop

 end subroutine read_v3_new_unformatted_sy



!============================================================
 subroutine read_v3_new_unformatted_sq(nk,ndn)

! This is modified version of read_v3_new_unformatted.
! It will use multi-dimensional array for v3

 use phi3
 use phi3_sy
 use io2
 use lattice
 implicit none
 integer j,unt
 real(8) xx,yy
 integer nk,ndn

 integer n1,n2,n3,l1,l2,l3

 unt = 111
! open(unt,file='v33sq.dat',status='old',form='UNFORMATTED')
open(unt,file='v33sq.dat',status='old')

! read(unt,end=99) nv3
 read(unt,*,end=99) nv3
 call allocate_v33_sq(nk,nk,ndn,ndn,ndn)
 write(ulog,*)' OPENING v33.dat, reading ',nv3
 v33sq1=0
 do j=1,nv3
   read(unt,*,end=99) n1,n2,l1,l2,l3,xx
!    read(unt,end=99) n1,n2,l1,l2,l3,xx
    if (V33sq1(n1,n2,l1,l2,l3) .ne. 0 ) then
       write(ulog,*) 'V33_sq duplicate'
       stop
    endif
    !write(*,*) 'n1,n2,l1,l2,l3,xx,yy',n1,n2,l1,l2,l3,xx,yy

    v33sq1(n1,n2,l1,l2,l3)=xx
    !write(*,*) 'v33sq1',v33sq1(n1,n2,l1,l2,l3)

 enddo
 close(unt)
 return

99 write(*,*)'v33 file end reached at line j=',j
   stop

 end subroutine read_v3_new_unformatted_sq


!===========================================================
subroutine cal_F2(nk,kp,ndn)

use constants
use exactBTE2
use params
use lattice
use io2

implicit none

integer nk,ndn
integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
integer nk1,nk2,nk3, nk3_proc1, nk3_proc2,xyz
real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3)
real(8) kp(3,nk)
real(8) Avalue(nk,ndn,3)

integer n1,n2,n3,l1,l2,l3


Avalue=0
loop1 : do n1=1,nk

    q1=kp(:,n1)
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
    if (n1 .ne. nk1) then
       write(ulog,*) 'n1,nk1,inside=',n1,nk1,inside
       write(ulog,*) 'q1=',q1
       write(ulog,*) 'ijk1=',i1,j1,k1
       stop
    endif


    do l1=1,ndn
    loopxyz : do xyz=1,3

       loop2: do n2=1,nk

          q2=kp(:,n2)
          call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
          if (n2 .ne. nk2) then
               write(ulog,*) 'n2,nk2,inside=',n2,nk2,inside
               write(ulog,*) 'q2=',q2
               write(ulog,*) 'ijk3=',i2,j2,k2
               stop
           endif

           q3_proc1=q1+q2       ! for class1 process
           q3_proc2=q1-q2       ! for class2 process

           call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside)
           call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside)

           do l2=1,ndn
           do l3=1,ndn

               Avalue(n1,l1,xyz)=Avalue(n1,l1,xyz) + P1(n1,n2,l1,l2,l3)*(F1(nk3_proc1,l3,xyz)-F1(n2,l2,xyz)) &
               &  + 0.5*P2(n1,n2,l1,l2,l3)*(F1(n2,l2,xyz)+F1(nk3_proc2,l3,xyz))

           enddo
           enddo
        enddo loop2

        F2(n1,l1,xyz)=F_RTA(n1,l1,xyz) + Avalue(n1,l1,xyz)/Qvalue(n1,l1)
        diff(n1,l1,xyz)=(F2(n1,l1,xyz)-F1(n1,l1,xyz))/F2(n1,l1,xyz)

    enddo loopxyz
    enddo
enddo loop1

write(*,*) F_RTA(1,1,1), F1(1,1,1), F2(1,1,1)


end subroutine cal_F2


!==============================================================
 subroutine matrix_elt_sy(q1,q2,l1,l2,l3,w33_proc1,w33_proc2,inside)
! for input modes (qi,li) i=1,2,3 calculates the 3-phonon matrix element w33(1,2,3) on the fly
 use kpoints
 use io2
 use phi3
 use svd_stuff
 use eigen
 use params
 use constants
 use atoms_force_constants
 use force_constants_module
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: l1,l2,l3
 real(8), intent(in) :: q1(3),q2(3)
 integer, intent(out) :: inside
 complex(8), intent (out) :: w33_proc1, w33_proc2
 integer i3,nk1,nk2,nk3,i0,al,be,ga,j,j0,k,k0,t,ta1,ta2,ta3,i1,j1,k1, nk3_proc1, nk3_proc2
 real(8) q3(3),mi,mj,mk,rx2(3),rx3(3),den, q3_proc1(3), q3_proc2(3), den_proc1, den_proc2
 complex(8) xx,eiqr, xx_proc1, xx_proc2, eiqr_proc1, eiqr_proc2

! write(*,*)'matrix_elt',q1,q2,l1,l2,l3
! write(*,*)'matrix_elt, const33=',const33
 call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
 call get_k_info(q2,NC,nk2,i1,j1,k1,g1,g2,g3,inside)
 w33_proc1 = cmplx(0d0,0d0)
 w33_proc2=cmplx(0d0,0d0)
q3_proc1 = q1+q2
q3_proc2 = q1-q2

 call get_k_info(q3_proc1,NC,nk3_proc1,i1,j1,k1,g1,g2,g3,inside)
 call get_k_info(q3_proc2,NC,nk3_proc2,i1,j1,k1,g1,g2,g3,inside)

! write(*,*)'matrix_elt: q3=',q3,nk3,inside
 xx_proc1 = cmplx(0d0,0d0)
 xx_proc2 = cmplx(0d0,0d0)

 tloop: do t=1,nterms(3)

! write(*,*)'matrix_elt: tloop=',t,nterms(3)
       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop ! first index must be in primitive cell
       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j) - atompos(:,j0)
       rx3(:) = atompos(:,k) - atompos(:,k0)
       eiqr_proc1 = exp( ci* ((q2 .dot. rx2) + (q3_proc1 .dot. rx3)) )
       xx_proc1 = xx_proc1 + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr_proc1 * &
&      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3_proc1)

       eiqr_proc2 = exp( ci* ((q2 .dot. rx2) + (q3_proc2 .dot. rx3)) )
       xx_proc2 = xx_proc2 + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr_proc2 * &
&      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3_proc2)


 enddo tloop
 den_proc1 = sqrt(8*sqrt(abs(eigenval(l1,nk1)*eigenval(l2,nk2)*eigenval(l3,nk3_proc1))))
 if (den_proc1.ne.0) then
    w33_proc1 = xx_proc1 / den_proc1 * const33
 else
    write(*,*)'MATRIX_ELT: den=0 ',den_proc1
    stop
 endif
 den_proc2 = sqrt(8*sqrt(abs(eigenval(l1,nk1)*eigenval(l2,nk2)*eigenval(l3,nk3_proc2))))
 if (den_proc2.ne.0) then
    w33_proc2 = xx_proc2 / den_proc2 * const33
 else
    write(*,*)'MATRIX_ELT: den=0 ',den_proc2
    stop
 endif

 end subroutine matrix_elt_sy



!=========================================================================================
subroutine distribution(nk,ndn,tempk)
! calculate frequency and distribution function and store.
! frequency is in cm^-1


use exactBTE2
use eigen
use params
use constants

implicit none

integer nk,ndn
integer i,ii
real(8) nbe,tempk,temp, omega

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

do i=1,nk
do ii=1,ndn
frequency(i,ii)=sqrt(eigenval(ii,i)) * cnst
dist(i,ii)=nbe(frequency(i,ii),temp,classical)
enddo
enddo

end subroutine distribution


!==========================================================================================
subroutine cal_freq(nk,ndn)
! calculate frequency and distribution function and store.
! frequency is in cm^-1

use exactBTE2
use eigen
use params
use constants

implicit none

integer nk,ndn
integer i,ii
real(8) nbe,tempk,temp, omega


do i=1,nk
do ii=1,ndn
frequency(i,ii)=sqrt(eigenval(ii,i)) * cnst
enddo
enddo

end subroutine cal_freq



!==================================================================================
subroutine cal_eiqr2(nk,ndn)
! calculate eiqr2, eivec2 for efficient calculation of matrix elements.
use kpoints
use eigen
use pre_matrix_elt
use lattice
use params
 use atoms_force_constants
 use force_constants_module
 use io2
 use phi3
 use svd_stuff
use constants

implicit none

integer nk, max_atomnum, ndn
integer i,j,xyz,ta,j0


max_atomnum = maxval(iatomterm_3(2,:))
write(*,*) 'max_atomnum=',max_atomnum
call allocate_eiqr2(nk,max_atomnum)
call allocate_eivec2(max_atomnum,ndn,nk)
! determine maximum of atom numbers in fc3.dat
! allocate eiqr2, eivec2

eiqr2(:,:)=0
eivec2(:,:,:,:)=0

do i=1,nk
  do j=1,max_atomnum
    j0 = iatomcell0(j)
    eiqr2(i,j)=cdexp(ci * dot_product(kpc(:,i),atompos(:,j)-atompos(:,j0)))
  enddo
enddo


do j=1,max_atomnum
  do xyz=1,3
     ta= xyz + 3*(iatomcell0(j)-1)
     eivec2(j,xyz,:,:)=eigenvec(ta,:,:)
  enddo
enddo

write(*,*) 'eiqr2(1,1)=',eiqr2(1,1)
write(*,*) 'eivec2(1,1,1,1)=',eivec2(1,1,1,1)

!eiqr2(nk2,j) where nk2 is q vector id, j is atom number
!eivec2(j,be,l2,nk2) where j is atom number, be is xyz



end subroutine cal_eiqr2





!============================================================================================
double precision function wtime()
! function for wall time clock. Reading large files takes more time than cpu time.
  integer(8) C,R,M
    CALL SYSTEM_CLOCK (COUNT=C, COUNT_RATE=R, COUNT_MAX=M)
    wtime = dble(C)/dble(R)
    return
  end



!============================================================================================
subroutine mode_thermal_RTA(tempk,nk,ndn,kp)
! calculate thermal conductivity using the direct formula with RTA.
! To test this gives the exactly same results from K1's original code



 use kpoints
 use io2
 use phi3
 use lattice
! use geometry
! use atoms_force_constants
! use force_constants_module
! use svd_stuff
 use params
 use constants
 use om_dos
! use coll_matrix
 use exactBTE2
 use phi3_sy
 use eigen



implicit none

integer i,ii,indx,ndn,j,jj,kk, n1,n2,n3,l1,l2,l3
integer i1, j1, k1
integer nk, nk_subset, k_shift, ksub_size2
integer nk1, nk2, nk3, nk2n
real(8) kp(3,nk), q1(3), q2(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
real(8) temp, tempk, nb1   ! in cm^-1
integer ksubset(2)
real(8) omega1, omega2, omega3, omega3_proc1, omega3_proc2, V3sq1, V3sq2, col1, col2

integer i2,ii2,j2,jj2,k2,kk2, inside1, inside2, inside, inside_proc1, inside_proc2
integer inside_proc1n, inside_proc2n
integer i3,j3,k3,nk3_proc1,nk3_proc2

 character(99) filename, filename_temp

!real cputim1, cputim2, cputim3, cputim0, cputim1_wall, cputim2_wall, cputim0_wall, wtime
real cputim1, cputim2, cputim3, cputim0
real(8) cputim1_wall, cputim2_wall, cputim0_wall, wtime

integer ucol, unt, nv3_2, n1_2, n2_2, l1_2, l2_2, l3_2

real(8), allocatable :: v33sq_sy(:,:,:,:,:),tauinv_RTA(:,:)
real(8) tauinv_RTA_N(nk,ndn), tauinv_RTA_U(nk,ndn)
real(8) kappa_k_RTA2(nk,ndn,3,3), kappa_RTA2(ndn,3,3), nbe, delta_l

real(8) tauinv_temp1, tauinv_temp2

unt=112
ucol=9010

allocate(v33sq_sy(nk,nk,ndn,ndn,ndn))
v33sq=0
 allocate(tauinv_RTA(nk,ndn))
 tauinv_RTA=0


indx=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

!call cpu_time(cputim0)
!cputim0_wall=wtime()
write(*,*) 'entering mode_thermal_RTA'


open(unt+11,file='iself.check.dat',status='unknown')
open(unt+12,file='dist_diff.dat',status='unknown')
open(unt+13,file='dist_diff2.dat',status='unknown')



loop1: do n1=1,nk

    q1=kp(:,n1)
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
    if (n1 .ne. nk1) then
       write(ulog,*) 'n1,nk1,inside=',n1,nk1,inside
       write(ulog,*) 'q1=',q1
       write(ulog,*) 'ijk1=',i1,j1,k1
       stop
    endif

  do l1=1,ndn

    loop2 : do n2=1,nk

          q2=kp(:,n2)
          call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
          if (n2 .ne. nk2) then
               write(ulog,*) 'n2,nk2,inside=',n2,nk2,inside
               write(ulog,*) 'q2=',q2
               write(ulog,*) 'ijk3=',i2,j2,k2
               stop
           endif

                 q3=-q1-q2
                 call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)

                 q3_proc1=q1+q2
                 q3_proc1n=-q1-q2
                 q3_proc2=q1-q2
                 q3_proc2n=-q1+q2
                 call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside_proc1)
                 call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside_proc2)
                 call get_k_info(q3_proc1n,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside_proc1n)
                 call get_k_info(q3_proc2n,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside_proc2n)


           do l2=1, ndn
           do l3=1, ndn

                 V3sq1=v33sq_sy(n1,n2,l1,l2,l3) * (2*pi)**2 / nk

                 omega1=frequency(nk1,l1)
                 omega2=frequency(nk2,l2)
                 omega3=frequency(nk3,l3)
                 omega3_proc1=frequency(nk3_proc1,l3)
                 omega3_proc2=frequency(nk3_proc2,l3)

! --- with below, results are exactly same as one from K1's original code.
!                 tauinv_RTA(n1,l1)=tauinv_RTA(n1,l1) - (1./2)*v3sq1 * &
!&                ((dist(nk3,l3)-dist(nk2,l2))*delta_l(omega1+omega2-omega3,etaz)  &
!&                - (dist(nk3,l3)-dist(nk2,l2))*delta_l(omega1-omega2+omega3,etaz)  &
!&                - (1+dist(nk2,l2)+dist(nk3,l3))*delta_l(omega1-omega2-omega3,etaz)  &
!&                + (1+dist(nk2,l2)+dist(nk3,l3))*delta_l(omega1+omega2+omega3,etaz))

!write(*,*) 'n1,n2,l1,l2,l3,nk1n,nk2n',n1,n2,l1,l2,l3,nk1n,nk2n
!write(*,*) 'v33sq(nk1n,l1,nk2n,l2,l3)',v33sq(nk1n,l1,nk2n,l2,l3)
!write(*,*) 'v33sq(nk1n,l1,nk2,l2,l3)',v33sq(nk1n,l1,nk2,l2,l3)
!write(*,*) 'v33sq=',v33sq(nk1n,l1,nk2n,l2,l3), v33sq(nk1n,l1,nk2,l2,l3)

! --- below is from collision matrix. direct derivation by neglecting coupling terms
                tauinv_temp1 = (1./2)* V3sq1 * (2*dist(n2,l2)*(dist(nk3,l3)+1)/(dist(n1,l1)+1))*delta_l(omega1+omega2-omega3,etaz)
                tauinv_temp2 = (1./2)* V3sq1 * ((dist(n2,l2)+1)*(dist(nk3,l3)+1)/(dist(n1,l1)+1))*delta_l(omega1-omega2-omega3,etaz)


!--- below is from eq. 3.68
!                 tauinv_temp1 = (1./2)* V3sq1 * (2*(dist(n2,l2)-dist(nk3,l3)))*delta_l(omega1+omega2-omega3,etaz)
!                 tauinv_temp2 = (1./2)* V3sq1 * (1+dist(n2,l2)+dist(nk3,l3))*delta_l(omega1-omega2-omega3,etaz)

                 if (inside .eq. 1) then   ! normal process
                     tauinv_RTA_N(n1,l1)=tauinv_RTA_N(n1,l1) + tauinv_temp1
                 else
                     tauinv_RTA_U(n1,l1)=tauinv_RTA_U(n1,l1) + tauinv_temp1
                 endif

                 if (inside .eq. 1) then   ! normal process
                     tauinv_RTA_N(n1,l1)=tauinv_RTA_N(n1,l1) + tauinv_temp2
                 else
                     tauinv_RTA_U(n1,l1)=tauinv_RTA_U(n1,l1) + tauinv_temp2
                 endif


                 tauinv_RTA(n1,l1)=tauinv_RTA(n1,l1) + tauinv_temp1 + tauinv_temp2


! below is for displaying viloence of detailed balance
if (abs((dist(n2,l2)-dist(nk3,l3))- dist(n2,l2)*(dist(nk3,l3)+1)/(dist(n1,l1)+1)) > 1d5 ) then
write(unt+12,11) nk1,nk2,nk3,l1,l2,l3,dist(n1,l1),dist(n2,l2),dist(nk3,l3), &
& abs((dist(n2,l2)-dist(nk3,l3))- dist(n2,l2)*(dist(nk3,l3)+1)/(dist(n1,l1)+1)), delta_l(omega1+omega2-omega3,etaz)
endif

if (abs(1+dist(n2,l2)+dist(nk3,l3)- (dist(n2,l2)+1)*(dist(nk3,l3)+1)/(dist(n1,l1)+1)) > 1d5 ) then
write(unt+13,11) nk1,nk2,nk3,l1,l2,l3,dist(n1,l1),dist(n2,l2),dist(nk3,l3), &
& abs(1+dist(n2,l2)+dist(nk3,l3)- (dist(n2,l2)+1)*(dist(nk3,l3)+1)/(dist(n1,l1)+1)), delta_l(omega1-omega2-omega3,etaz)
endif

            enddo
            enddo
!            enddo
     enddo loop2

!write(*,*) tauinv_RTA(n1,l1)

enddo
enddo loop1

! print out selfenergy
do l1=1,ndn
do n1=1,nk
write(unt+11,11) l1, n1, tempk, frequency(n1,l1), 0.5*tauinv_RTA(n1,l1), 0.5*tauinv_RTA_N(n1,l1), 0.5*tauinv_RTA_U(n1,l1)
enddo
enddo

close(unt+11)
close(unt+12)
close(unt+13)

!---------- calculate thermal conductivity
nkF1_loop: do i=1,nk

      laF1_loop: do ii=1,ndn

          omega1=frequency(i,ii)
          nb1=dist(i,ii)

          alpha_loop: do j=1,3
          beta_loop: do jj=1,3

          kappa_k_RTA2(i,ii,j,jj) = (omega1*c_light*100.)**2.0 * (veloc(j,ii,i)*c_light) * &
     &     (veloc(jj,ii,i)*c_light) * nb1 * (nb1+1) * (1/(tauinv_RTA(i,ii)*c_light*100.))
       kappa_k_RTA2(i,ii,j,jj) = kappa_k_RTA2(i,ii,j,jj) * h_plank**2.0 / (nk*volume_r/1d30) / (k_b*tempk**2.0)
               kappa_RTA2(ii,j,jj)=kappa_RTA2(ii,j,jj) + kappa_k_RTA2(i,ii,j,jj)

          enddo beta_loop
          enddo alpha_loop

       enddo laF1_loop

    !   write(*,*) kappa_RTA2(1,2,2), kappa_RTA2(5,2,2)
enddo nkF1_loop

!----------- print out kapp_RTA2

open(ucol,file='kappa_RTA.dat',status='unknown')
write(ucol,*) 'xx,yy,zz, xy,xz,yz'
do ii=1,ndn
   write(ucol,9) ii,kappa_RTA2(ii,1,1),kappa_RTA2(ii,2,2),kappa_RTA2(ii,3,3),  &
   &     kappa_RTA2(ii,1,2),kappa_RTA2(ii,1,3),kappa_RTA2(ii,2,3)
enddo




close(unt)
close(ucol)






7 format(5(i6),2(2x,g14.8))
8 format(i6,3(2x,g14.8))
9 format(i3,6(2x,g14.8))
10 format(6(i6),9(2x,g14.8))
11 format(2(i6),9(2x,g14.8))



end subroutine mode_thermal_RTA





!============================================================================================

!---------------------------------------------------------------------------------------
subroutine calculate_v3sq_delta (ksubset,ndn,nk,kp)
!- read v33sq from v33sq from v33sq.xxx.dat
!- calculate delta function for coalescence and decay processes
!- write v33sq * delta to files

use io2
use exactBTE2
use mod_ksubset2
use params
use constants
use phi3
use lattice
use tetrahedron
use kpoints

implicit none

integer i,ii,indx,ndn,j,jj,kk
integer nk, nk_subset, k_shift, ksub_size2
real(8) kp(3,nk), q1(3), q2(3), q2n(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
real(8) temp, tempk, nb1   ! in cm^-1
integer ksubset(2)
integer nx, ny, nz

integer i2,ii2,j2,jj2,k2,kk2, inside1, inside2, inside, inside1n, inside2n
integer i3,j3,k3,nk3_proc1,nk3_proc2, nk3, nk2
integer n1,n2,l1,l2,l3,nk2n,i2n,j2n,k2n
integer j_tet, k_tet,w

character(99) filename_temp, filename

integer uQ, nv3_2,unt

real(8) cputim0, cputim1, cputim2, cputim3, cputim0_wall, cputim1_wall, wtime
real(8) omega1, V3sq1
real(8) F_arbi(nk)
real(8), allocatable :: v33sq_sy(:,:,:,:,:), P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)
! P1sm and P2sm are v33sq*delta(proc1) and v33sq*delta(proc2)

nx=NC(1); ny=NC(2); nz=NC(3)
unt=9400

call allocate_iter1(nk,nkc,ndn)
call cal_freq(nk,ndn)
write(ulog,*) 'cal_freq done'

write(ulog,*) 'entering calculate_v3sq_delta ...'

! calculate indx using nk and ksubset
! calculate ksub_size2
if (ksubset(2) .eq. nk) then   ! if this is the last file in v33sq.xxx.dat
   ksub_size2=ksubset(2)-ksubset(1)+1
   indx=(ksubset(1)-1)/ksub_size + 1
   if (mod(ksubset(1)-1,ksub_size) .ne. 0 ) then
       write(ulog,*) 'wrong indx. stop'
       stop
   endif

else
    ksub_size2=ksubset(2)-ksubset(1)+1
    indx=ksubset(2)/ksub_size           ! ksub_size is from params.phon
    if (mod(ksubset(2),ksub_size) .ne. 0) then
       write(ulog,*) 'wrong indx. stop'
       stop
    endif

endif

k_shift=ksubset(1)-1

write(ulog,*) 'Check indx is same as folder name, indx=',indx



! Read v33sq from v33sq.xxx.dat in the main directory
write(filename_temp,fmt="(a,i3.3,a)") "v33sq.",indx,".dat"
filename=trim(v3path)//filename_temp
open(unt,file=filename,status='old',form='unformatted')
read(unt) nv3_2

if (allocated(v33sq_sy)) then
   deallocate(v33sq_sy)
endif

allocate(v33sq_sy(ksub_size2,nk,ndn,ndn,ndn))

do n1=1,ksub_size2
do l1=1,ndn
do n2=1,nk
do l2=1,ndn
do l3=1,ndn
       read(unt) v33sq_sy(n1,n2,l1,l2,l3)
enddo
enddo
enddo
enddo
enddo
close(unt)


! calculate delta function
if (allocated(P1sm)) then
    deallocate(P1sm)
endif
if (allocated(P2sm)) then
    deallocate(P2sm)
endif

allocate(P1sm(ksub_size2,nk,ndn,ndn,ndn))
allocate(P2sm(ksub_size2,nk,ndn,ndn,ndn))

call cpu_time(cputim1)

do n1=1,ksub_size2

    i=n1+k_shift

    laF1_loop: do ii=1,ndn

    write(ulog,*) 'i,ii=',i,ii

    omega1=frequency(i,ii)        ! set fixed omega1 in delta(omega1-omega2+-omega3)

           laF2_loop: do jj=1,ndn
           laF3_loop: do kk=1,ndn


           ! process 1
           ! set arbitrary function here. (distribution function...)
           ! call eigen_tet using w2-w3
           ! call weight_tet

           ! set arbitrary function, F_arbi(j)
           do j=1,nk   ! for set arbitrary function

                q1=kp(:,i)
                q2=kp(:,j)

!                q3=-q1-q2
                q3_proc1=q1+q2
!                q3_proc2=q1-q2

                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
!                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)


                V3sq1=v33sq_sy(i-k_shift,j,ii,jj,kk) * (2*pi)**2 / nk
                if (jj<=3 .and. i2.eq.NC(1)/2+1 .and. j2.eq.NC(2)/2+1 .and. k2.eq.NC(3)/2+1 ) then  ! Acoustic modes at Gamma point
                    V3sq1=0.0                                                                        ! to enforce matrix_elt = 0
                endif
                if (kk<=3 .and. i3.eq.NC(1)/2+1 .and. j3.eq.NC(2)/2+1 .and. k3.eq.NC(3)/2+1 ) then
                    V3sq1=0.0
                endif

                F_arbi(j)=V3sq1

            enddo

            ! eigen_tet: set omega2+omega3
            do j_tet=1,nk*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                    w = tet(j_tet)%p(k_tet)%i                          !label for kpt(3,i)

                    q1=kp(:,i)
                    q2=kp(:,w)

                    q3=-q1-q2
                    q3_proc1=q1+q2

                    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                    call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)

!assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
                    tet(j_tet)%p(k_tet)%w=-frequency(w,jj)+frequency(nk3_proc1,kk)
               enddo
            enddo

           call weight_tet(omega1)   ! calculate weighting factors in tetrahedron method

           do j_tet=1,nk*6                       !iterate over all tetrahedrons
               do k_tet=1,4                      !iterate over four corners of tetrahedron
                   w = tet(j_tet)%p(k_tet)%i                ! k point
                   P1sm(i-k_shift,w,ii,jj,kk) = P1sm(i-k_shift,w,ii,jj,kk) + &
&                  tet(j_tet)%p(k_tet)%c/6d0 * F_arbi(w)
               enddo
           enddo


           ! process 2
           ! set arbitrary function here. (distribution function...)
           ! call eigen_tet using w2-w3
           ! call weight_tet

           do j=1,nk   ! for set arbitrary function

                q1=kp(:,i)
                q2=kp(:,j)
                q2n=-kp(:,j)

!                q3=-q1-q2n
!                q3_proc1=q1+q2
                q3_proc2=q1+q2

                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                call get_k_info(q2n,NC,nk2n,i2n,j2n,k2n,g1,g2,g3,inside)
!                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)


                V3sq1=v33sq_sy(i-k_shift,j,ii,jj,kk) * (2*pi)**2 / nk
! Acoustic modes at Gamma point
                if (jj<=3 .and. i2n.eq.NC(1)/2+1 .and. j2n.eq.NC(2)/2+1 .and. k2n.eq.NC(3)/2+1 ) then
                    V3sq1=0.0                  ! to enforce matrix_elt = 0
                endif
                if (kk<=3 .and. i3.eq.NC(1)/2+1 .and. j3.eq.NC(2)/2+1 .and. k3.eq.NC(3)/2+1 ) then
                    V3sq1=0.0
                endif

                F_arbi(nk2n)=V3sq1  ! no 2pi at the front since delta is in linear frequency.

            enddo

            ! eigen_tet: set omega2+omega3
            do j_tet=1,nk*6             !iterate over all tetrahedrons
               do k_tet=1,4             !iterate over four corners of tetrahedron
                    w = tet(j_tet)%p(k_tet)%i                          !label for kpt(3,i)

                    q1=kp(:,i)
                    q2=kp(:,w)

!                    q3=-q1-q2
                    q3_proc2=q1-q2

!                    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                    call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside1)

!assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
                    tet(j_tet)%p(k_tet)%w=frequency(w,jj)+frequency(nk3_proc2,kk)
               enddo
            enddo

! k1           call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors
           call weight_tet(omega1)   ! calculate weighting factors

           do j_tet=1,nk*6                       !iterate over all tetrahedrons
               do k_tet=1,4                      !iterate over four corners of tetrahedron
                   w = tet(j_tet)%p(k_tet)%i                ! k point
                   P2sm(i-k_shift,w,ii,jj,kk) = P2sm(i-k_shift,w,ii,jj,kk) + tet(j_tet)%p(k_tet)%c/6d0 * F_arbi(w)
               enddo
           enddo

           enddo laF3_loop
           enddo laF2_loop


    enddo laF1_loop


    call cpu_time(cputim2)
    write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1

enddo



! write P1sm and P2sm
open(unt+1,file='v33sq_delta.dat',status='unknown',form='unformatted')
write(unt+1) nv3_2
!write(*,*) 'nv3_2=',nv3_2
do n1=1,ksub_size2
do l1=1,ndn
do n2=1,nk
do l2=1,ndn
do l3=1,ndn
   write(unt+1) P1sm(n1,n2,l1,l2,l3), P2sm(n1,n2,l1,l2,l3)
!write(*,*) 'n1,n2,l1,l2,l3',n1,n2,l1,l2,l3
!write(*,*) 'P1sm,P2sm',P1sm(n1,n2,l1,l2,l3),P2sm(n1,n2,l1,l2,l3)

enddo
enddo
enddo
enddo
!write(*,*) 'n1,n2,l1,l2,l3',n1,n2,l1,l2,l3
!write(*,*) 'P1sm,P2sm',P1sm(n1,n2,l1,l2,l3),P2sm(n1,n2,l1,l2,l3)
enddo

close(unt+1)




7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during calculate_v3sq_delta...'
!   stop
!endif



end subroutine calculate_v3sq_delta


!==========================================================
 subroutine read_params
 use io2
 use om_dos
 use params
 use lattice
 use phi3
 use kpoints
 use atoms_force_constants
 implicit none
! integer i

  open(uparams,file='params.phon',status='old')
  write(6,*) 'READ_PARAMS: opening params.phon'
!  read(uparams,*) natom_type   ! # of different elements present in prim cell
!  write(6,*) 'READ_PARAMS: just read natom_type ',natom_type
! call allocate_mass(natom_type)
! read(uparams,*) (mas(i),i=1,natom_type)  ! in the same order as in POTCAR
! read(uparams,*) (atname(i),i=1,natom_type)  ! in the same order as in POTCAR
! read(uparams,*) natoms        ! # of atoms in primitive cell

  read(uparams,*)nc1,nc2,nc3,dos_bs_only ! coarse kmesh for big calculations
  NC(1)=nc1 ;NC(2)=nc2 ;NC(3)=nc3
  write(ulog,*) 'READ_PARAMS: just read nib ',nc1
  read(uparams,*)lamin,lamax
  read(uparams,*)shftx,shfty,shftz ! kmesh for dos calculation
  read(uparams,*)wmesh,wmax  ! no of w mesh points and max frequency for DOS
  write(ulog,*) 'READ_PARAMS: just read wmax ',wmax
  read(uparams,*)width,etaz  ! width of gaussian broadening for DOS,imag part of om
  write(ulog,*) 'READ_PARAMS: just read width,etaz ',width,etaz
  read(uparams,*)tau0        ! large relaxation time: 1/tau0 to be added to gamma
  write(ulog,*) 'READ_PARAMS: just read tau0 ',tau0
  read(uparams,*)verbose
  read(uparams,*)wshift      ! small shift added to make negative modes go away
  write(ulog,*) 'READ_PARAMS: just read wshift ',wshift
  read(uparams,*)tmin,tmax,ntemp      ! temps in Kelvin to calculate thermal expansion and other thermal ppties
  read(uparams,*)iter,readv3,usetetra,isvd      ! if=1 then read from a file, else generate
  read(uparams,*)ncpu
  read(uparams,'(a)')v3path
  write(ulog,*) 'READ_PARAMS: path is:',v3path
  read(uparams,*)max_iter, n_dig_acc 
!, conv_error, conv_max_error, conv_diff, conv_max_diff, &
!  & conv_diff_kap, conv_max_diff_kap,conv_iter,update
   ! maximum iteration number for iterative solution, convergence criteria, update ratio
  read(uparams,*)v3_threshold        ! in 1/cm keep all v33 of norm above v3_threshold
  read(uparams,*)classical     ! if =1 then use classical distribution (kT/hw)
  read(uparams,*)cal_cross,qcros     ! if =1 then calculate the cross section at qcros (reduced U)
  read(uparams,*)threemtrx          ! if =1 then calculate the 3-phonon matrix elements
  read(uparams,*)scalelengths       ! multiply all lattice and atomic coordinates by scalelength

  close(uparams)

  write(ulog,*) 'READ_PARAMS: read and kpoints allocated'
  write(ulog,3)'nc1,nc2,nc3=',nc1,nc2,nc3
  write(ulog,*)'wmesh, wmax=',wmesh,wmax
  write(ulog,*)'readv3=',readv3  !,writev3
! etaz = max(etaz,wmax/wmesh)
! write(ulog,*)' in case etaz is too small we use the following'
! write(ulog,*)'width,etaz =',width,etaz
  write(ulog,*)'Tmin,Tmax(K=',tmin,tmax
  if (classical .eq. 1) then
     write(ulog,*)'Calculation is done in the CLASSICAL limit'
  else
     write(ulog,*)'Calculation is done in the QUANTUM limit'
  endif
  if (dos_bs_only) then
     write(ulog,*)' DOS will only be calculated (on the coarse mesh) ',nc1,nc2,nc3
!    write(ulog,*)' and the program will stop after that!'
  endif

3 format(a,6(1x,i6))
 end subroutine read_params
!==========================================================
 subroutine read_input_fit
 use io2
 use params
 use lattice
 use force_constants_module
 use atoms_force_constants
 implicit none
 integer i,counter,label
 real(8) junk
 character jjj*1
 open(uparams,file='params.inp',status='old')

 read(uparams,*) latticeparameters   ! (a,b,c,alpha,beta,gamma)
 read(uparams,*) primitivelattice     ! prim latt vectors in terms of conventional above
 read(uparams,*) lattice_parameter
 read(uparams,*) maxneighbors
 read(uparams,*) include_fc    ! if=1 include this rank
! read(uparams,*) nshells       ! # of neigr-shells to use for each rank
 read(uparams,*)  i !junk
! read(uparams,*) junk
! read(uparams,*) junk
 read(uparams,*) tolerance  ! this really the flags...
 read(uparams,*) svdcut
 read(uparams,*) junk
 read(uparams,*) natom_type   ! # of different elements present in prim cell
!read(uparams,*) jjj
! allocate(natom(natom_type))
 write(ulog,*) 'READ_INPUT_FIT: natom_type ',natom_type
 call allocate_mass(natom_type)
 read(uparams,*) (mas(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) (atname(i),i=1,natom_type)  ! in the same order as in POTCAR
 read(uparams,*) natoms0        ! # of atoms in primitive cell
 read(uparams,*) nshells(1,1:natoms0)       ! # of neigr-shells to use for each rank
 read(uparams,*) nshells(2,1:natoms0)       ! # of neigr-shells to use for each rank
 read(uparams,*) nshells(3,1:natoms0)       ! # of neigr-shells to use for each rank
 read(uparams,*) nshells(4,1:natoms0)       ! # of neigr-shells to use for each rank
 write(ulog,*)'reading ',natoms0,' atoms in the primitive cell'
 call allocate_primcell(natoms0)
! indices of atoms in primitive cell must be same order as in POSCAR1
 counter = 0
 natom = 0
 do i=1,natoms0
! positions here must be D format in conventional cell for use by fcs_init
    read(uparams,*) label,atom_type(label),atom0(label)%equilibrium_pos
! this is useless because natom(:) is read from POSCAR1
    if (atom_type(i).eq.counter) then
        natom(atom_type(i)) = natom(atom_type(i))+1
    else
        natom(atom_type(i)) = 1
    Endif
    counter = atom_type(i)
 write(ulog,*)'reading  atom #',i, counter
 enddo

 write(ulog,*)'reading done, closing the params.inp file'
 close(uparams)

 do i=1,natoms0
    atom0(i)%name = atname(atom_type(i))
    atom0(i)%at_type  = atom_type(i)
    atom0(i)%mass = mas(atom_type(i))
    atompos0(:,i) = atom0(i)%equilibrium_pos  
 enddo

 maxshells = maxneighbors

 do i=1,4
!   if (nshells(i) .gt. maxshells) then
    if (maxval(nshells(i,:)) .gt. maxshells) then
!      write(ulog,*)' READ_INPUT: nshell> maxshells ',i,nshells(i),maxshells
       write(ulog,*)' READ_INPUT: nshell> maxshells ',i,maxval(nshells(i,:)),maxshells
       write(ulog,*)' Either increase maxshells or lower nshell for that rank'
       stop
    endif
 enddo

end subroutine read_input_fit
!==========================================================
 subroutine read_lattice
 use io2
 use params
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use constants
 implicit none
 character line*90,name*2
 integer i,j,tt,ttyp,natoms2,n1  !,at_count,n2,n3
 real(8) mss !a1,a2,a3,om,a,b,c,
 type(vector) tau1,vv

open(ufc0,file='lat_fc.dat' ,status='old')

do while (line(1:22) .ne. '  Crystal data: transl')
 read(ufc0,'(a)') line
enddo
 read(ufc0,*) r1   ! coordinates of the primitive cell vectors
 read(ufc0,*) r2
 read(ufc0,*) r3

 r1= scalelengths * r1   !! added by K1 5/6/2015
 r2= scalelengths * r2   !! added by K1 5/6/2015
 r3= scalelengths * r3   !! added by K1 5/6/2015

! now generate the lattice ppties in real and recip spaces
  call make_reciprocal_lattice_2pi(r1,r2,r3,g1,g2,g3)
  call make_reciprocal_lattice(g1,g2,g3,rr1,rr2,rr3)
  call calculate_volume(r1,r2,r3,volume_r)
  call calculate_volume(g1,g2,g3,volume_g)
  box(1) = length(r1)
  box(2) = length(r2)
  box(3) = length(r3)
  boxg(1) = length(g1)
  boxg(2) = length(g2)
  boxg(3) = length(g3)
  write(ulog,3)' r1= ',r1
  write(ulog,3)' r2= ',r2
  write(ulog,3)' r3= ',r3
  write(ulog,3)' box  = ',box
  write(ulog,3)' g1= ',g1
  write(ulog,3)' g2= ',g2
  write(ulog,3)' g3= ',g3
  write(ulog,3)' boxg = ',boxg
  write(ulog,3)' volume_r,volume_g = ',volume_r,volume_g

 read(ufc0,*) line
! read(ufc,*) natoms0
! ############ check compatibiliy with input file "params.inp" ###########
 read(ufc0,*) n1
 if (n1 .ne. natoms0) then
    write(ulog,*)' error in reading natoms0 ',n1,natoms0
    stop
 endif
! call allocate_primcell(natoms0)
 do i=1,natoms0
!  read(ufc1,*)atom0(i)%tau, atom_type(i),atom0(i)%equilibrium_pos,atom0(i)%mass
    read(ufc0,*)tt,name,ttyp,tau1,mss  ! this has the coordinates and atomic types and masses

    tau1 = scalelengths * tau1  !! added by K1 5/6/2015

!   atompos0(:,i) = atom0(i)%equilibrium_pos
    vv%x= tau1.dot.g1 ; vv%y= tau1.dot.g2 ; vv%z= tau1.dot.g3
    vv = (1/2d0/pi)*vv
    if(tt .eq. i .and. ttyp .eq. atom0(i)%at_type .and. atom0(i)%mass .eq. mss) then
!&     length(atom0(i)%equilibrium_pos-vv).lt.1d-4 .and.
      write(ulog,*)'####### compatibility with params.inp checked  ######### ',i
    else
      write(ulog,*)'READ_LATTICE: data in params.inp and lat_fc.dat are not compatible'
      write(ulog,*)' i, atom type, mass '
      write(ulog,*) i, tt, atom0(i)%at_type,ttyp, atom0(i)%mass,mss
! fix it so that if eq pos are the same module translarion vectors, it's OK
      write(ulog,3)' reduced coordinates ',atom0(i)%equilibrium_pos,vv
      stop
    endif
 enddo

 read(ufc0,*) line
 read(ufc0,*) line
 read(ufc0,*) (include_fc(j),j=1,4)
 write(ulog,*) (include_fc(j),j=1,4)
 read(ufc0,*) line
 write(ulog,*) line
 read(ufc0,*) (nterms(j),j=1,4)
 write(ulog,*) (nterms(j),j=1,4)
 read(ufc0,*) line
 write(ulog,*) line
 read(ufc0,*) (ngroups(j),j=1,4)
 write(ulog,*) (ngroups(j),j=1,4)
 read(ufc0,*) line
 write(ulog,*) line
 read(ufc0,*) natoms2
 write(ulog,*) natoms2

  maxatoms=2500; imaxat=1
  do while (imaxat.ne.0)
     maxatoms=maxatoms+300
     write(6,*)' maxatoms=',maxatoms
     write(ulog,*)' maxatoms=',maxatoms
     if (allocated(iatomcell0)) deallocate(iatomcell0)
     if (allocated(iatomcell))  deallocate(iatomcell)
     if (allocated(atompos))    deallocate(atompos)
     allocate(atompos(3,maxatoms), iatomcell(3,maxatoms),iatomcell0(maxatoms))
! find symmetry operations of the crystal
     call force_constants_init(latticeparameters,primitivelattice,natoms0,  &
      &     atom_type,atompos0)  !(n,natom_type,atompos0)
  enddo

 if (natoms.ne.natoms2) then
    write(ulog,*) "natoms read in lat_fc=",natoms2
    write(ulog,*) "while natoms generated by fc_init=",natoms
    write(ulog,*) "check the number of shells in params.inp, "
    write(ulog,*) "it is probably too large if natoms in fc_init is larger than in lat_fc"
    write(ulog,*) "in which case the program will stop"
 endif
 deallocate(atompos, iatomcell,iatomcell0)
 allocate(atompos(3,natoms), iatomcell(3,natoms),iatomcell0(natoms))
 do j=1,natoms
    read(ufc0,*,end=99)i,atompos(:,i), iatomcell0(i),iatomcell(:,i)

    atompos(:,i)=atompos(:,i)  * scalelengths  !! added by K1 5/6/2015

    write(ulog,4)i,atompos(:,i), iatomcell0(i),iatomcell(:,i)
 enddo
3 format(a,9(2x,f11.6))
4 format(i6,3(2x,g11.5),4(3x,i4))

 close(ufc0)
 write(ulog,*)' COORDINATES read successfully '
 return

99 write(*,*) "End of lat_fc.dat file reached!"
   write(*,*) "check the number of shells in params.inp, it is probably too large!"
   stop

 end subroutine read_lattice
!===========================================================
 subroutine read_fc23
 use io2
 use params
 use lattice
 use geometry
 use force_constants_module
 use atoms_force_constants
 use svd_stuff
 implicit none
 character line*90
 integer t,rank,res,j,mx

  open(ufc2,file='fc2.dat'    ,status='old')
  open(ufc3,file='fc3.dat'    ,status='old')

! read(ufc2,'(a)') line
! read(ufc,*)include_fc(1),include_fc(2),include_fc(3),include_fc(4)
! read(ufc,'(a)') line
! read(ufc,*)nterms(1),nterms(2),nterms(3),nterms(4)
! read(ufc,'(a)') line
! read(ufc,*)ngroups(1),ngroups(2),ngroups(3),ngroups(4)
 res = 0
!----------------------------------------
! rank=1
! if ( include_fc(rank) .eq. 1 ) then
! mx = nterms(rank)
! allocate(iatomterm_1(rank,mx),ixyzterm_1(rank,mx),igroup_1(mx), &
! & ampterm_1(mx),fcs_1(ngroups(rank)))
! read(ufc,*)line
! do j=1,nterms(rank)
!       read(ufc,*)t,igroup_1(t), &
!& iatomterm_1(1,t),ixyzterm_1(1,t),  &
!& fcs_1(igroup_1(t)),ampterm_1(t)
! enddo
! res = res + igroup_1(nterms(rank))
! endif
!----------------------------------------
 rank=2
! if ( include_fc(rank) .eq. 1 ) then
 mx = nterms(rank)
 read(ufc2,'(a)') line
 do j=1,mx
    read(ufc2,*,end=92)t
 enddo
92 mx=j-1
 rewind(ufc2)
 nterms(rank)=mx
 write(ulog,*)'Rank=2, reading nterms(2)=',mx,' terms'
 allocate(iatomterm_2(rank,mx),ixyzterm_2(rank,mx),igroup_2(mx)) !,ampterm_2(mx))
! allocate(fcs_2(ngroups(rank)))
 allocate(fcs_2(mx),grun_fc(mx))
 read(ufc2,'(a)') line
 write(*  ,'(a)') line
 write(*,*) '********** FCs for rank=2',mx

 do j=1,mx
       read(ufc2,*)t,igroup_2(j), &
& iatomterm_2(1,j),ixyzterm_2(1,j),  &
& iatomterm_2(2,j),ixyzterm_2(2,j),  &
!& fcs_2(igroup_2(j)),ampterm_2(j)
& fcs_2(j)
 enddo
 res = res + igroup_2(nterms(rank))
! endif
!----------------------------------------
 rank=3
 if ( include_fc(rank) .eq. 1 ) then
 mx = nterms(rank)
 write(*,*) '********** FCs for rank=3',mx
 read(ufc3,'(a)') line
 do j=1,nterms(rank)
       read(ufc3,*,end=93)t
 enddo
93 mx=j-1
 rewind(ufc3)
 nterms(rank)=mx
 write(ulog,*)'Rank=3, reading nterms(3)=',mx,' terms'
 allocate(iatomterm_3(rank,mx),ixyzterm_3(rank,mx),igroup_3(mx) ) !,ampterm_3(mx),fcs_3(ngroups(rank)))
 allocate(fcs_3(mx))
 read(ufc3,'(a)') line
 write(*  ,'(a)') line
 do j=1,mx
!       read(ufc3,*)t,igroup_3(t), &
!& iatomterm_3(1,t),ixyzterm_3(1,t),  &
!& iatomterm_3(2,t),ixyzterm_3(2,t),  &
!& iatomterm_3(3,t),ixyzterm_3(3,t),  &
!& fcs_3(igroup_3(t)),ampterm_3(t)
       read(ufc3,*)t,igroup_3(j), &
& iatomterm_3(1,j),ixyzterm_3(1,j),  &
& iatomterm_3(2,j),ixyzterm_3(2,j),  &
& iatomterm_3(3,j),ixyzterm_3(3,j),  &
& fcs_3(j)
!& fcs_3(igroup_3(j)),ampterm_3(j)
 write(*,*) j,fcs_3(j)
 enddo
 res = res + igroup_3(nterms(rank))
 endif
 write(ulog,*)'READ_FC23: done!, number of groups read is=',res
 write(ulog,*)'READ_FC23:',fcs_2

 close(ufc2)
 close(ufc3)
 end subroutine read_fc23
!==========================================================
 subroutine set_dynamical_matrix(kpt,dynmat,ndim,ddyn)
 use params
 use lattice
 use kpoints
 use atoms_force_constants
 use geometry
 use svd_stuff
 implicit none
 integer, intent(in) :: ndim
 complex(8), intent(out) :: dynmat(ndim,ndim),ddyn(ndim,ndim,3)
 real(8), intent(in) :: kpt(3)
 complex(8) junk
 real(8) mi,mj,all,rr(3),delt(3)
 integer i0,j,j0,al,be,i3,j3,t !,ired

4 format(a,4(1x,i5),9(2x,f9.4))
9 format(i6,99(1x,f9.3))

 delt = wshift !*wshift
 ddyn   = cmplx(0d0,0d0)
 dynmat = cmplx(0d0,0d0)
 do i0=1,natoms0
 do al=1,3
    i3 = al+3*(i0-1)
    mi = atom0(i0)%mass
! write(ulog,*) 'i,al,mass=',i0,al,i3,mi
    tloop: do t=1,nterms(2)
       if ( i0 .eq. iatomterm_2(1,t) .and. al .eq. ixyzterm_2(1,t) ) then
          be =  ixyzterm_2(2,t)
          j  = iatomterm_2(2,t)
          j0 = iatomcell0(j)
          mj = atom0(j0)%mass
          j3 = be+3*(j0-1)
          rr = atompos(:,j)-atompos(:,j0)
!          ired = igroup_2(t)
!          junk = fcs_2(ired)*ampterm_2(t)* exp(ci*(kpt .dot. rr))/sqrt(mi*mj)
          junk = fcs_2(t)* exp(ci*(kpt .dot. rr))/sqrt(mi*mj)
!  write(ulog,4) 't,j,be,mass,dyn=',t,j,be,j3,mj,junk
          dynmat(i3,j3) = dynmat(i3,j3) + junk
          ddyn(i3,j3,:) = ddyn(i3,j3,:) + junk*ci*rr(:)
!         if (be.eq.al .and. j0.eq.i0 ) then
!         if (be.eq.al .and. (length(rr-r1) .myeq. 0d0)) then
! same type and same cartesian coord but not the same atom, just take the images
!     write(ulog,5)'i0,j0,j,i3,j3,rij=',i0,j0,j,i3,j3,rr
!            dynmat(i3,j3) = dynmat(i3,j3) - delt/2
!         endif
       endif
    enddo tloop
! this leaves  gamma eivals unchanged
!   dynmat(i3,i3) = dynmat(i3,i3) + delt(al)*(1-cos(kpt .dot. r1))
! but this one shifts everything up by delt=wshift
    dynmat(i3,i3) = dynmat(i3,i3) + delt(al)
 enddo
 enddo

! if (verbose) then
!  write(ulog,*)'SET_DYNAMICAL_MATRIX: d(dyn)/dk is:'
!  do t=1,3
!     write(ulog,*)'=======component of v ',t
!     do j=1,ndim
!        write(ulog,9) j, ddyn(j,:,t)
!     enddo
!  enddo
! endif

5 format(a,5i5,9(f8.3))
 all = sum(cdabs(dynmat(:,:)))/(ndim*ndim)
! make sure it is hermitian
 do t=1,ndim
    if (abs(aimag(dynmat(t,t))) .gt. 9d-4*abs(real(dynmat(t,t))) ) then
       write(ulog,*)' dynmat is not hermitian on its diagonal'
       write(ulog,*)' diagonal element i=',t,dynmat(t,t)
!      stop
!   else
       dynmat(t,t) = cmplx(real(dynmat(t,t)),0d0)
    endif
  do j=t+1,ndim-1
    if (abs(aimag(dynmat(t,j))+aimag(dynmat(j,t))) .gt. 9d-4*all ) then
       write(ulog,*)' dynmat is not hermitian in AIMAG of its off-diagonal elts'
       write(ulog,*)' off-diagonal element i,j=',t,j,dynmat(t,j),dynmat(j,t)
       write(ulog,*)' comparing to avg(abs(dynmat))=',all
!      stop
    elseif(abs(real(dynmat(t,j))-real(dynmat(j,t))) .gt. 9d-4*all ) then
       write(ulog,*)' dynmat is not hermitian in REAL of its off-diagonal elts'
       write(ulog,*)' off-diagonal element i,j=',t,j,dynmat(t,j),dynmat(j,t)
       write(ulog,*)' comparing to avg(abs(dynmat))=',all
!      stop
    else
    endif
! enforcing it to be hermitian!!
    mi=(aimag(dynmat(t,j))-aimag(dynmat(j,t)))/2
    mj=(real (dynmat(t,j))+real (dynmat(j,t)))/2
    dynmat(t,j) = cmplx(mj, mi)
    dynmat(j,t) = cmplx(mj,-mi)
  enddo
 enddo

! enforce it to be real if all imaginary componenents are very small
! all = sum(abs(aimag(dynmat(:,:))))/(ndim*ndim)
! if (all .lt. sum(abs(real(dynmat(:,:))))/(ndim*ndim)*1d-10) then
!  do t=1,ndim
!  do j=1,ndim
!     dynmat(t,j)=(dynmat(t,j)+conjg(dynmat(t,j)))/2d0
!  enddo
!  enddo
! endif

 end subroutine set_dynamical_matrix
!==========================================================
 subroutine diagonalize(n,mat,eival,nv,eivec,ier)
! n=size of mat; nv is the number of needed eigenvectors
 implicit none
 integer n,ier,zero,nv  !,i,j
 complex(8) mat(n,n),eivec(n,n)
 real(8)  eival(n),tol
! This is used by eigch
! real(8), allocatable :: w(:,:)
! integer, allocatable :: lw(:)
! This is used by ZHEGV
  real(8), allocatable :: rwork(:)
  complex(8), allocatable :: work(:)
  integer lwork


! n = size(mat(:,1))
 if (n .ne. size(eival) ) then
    write(   *,*)' EIGCH, size inconsistency:mat,eival=',n, size(eival)
    stop
 endif

! tol = -1.d0
! zero= 0  ! for now do not compute any eigenvectors
! ier=0
! allocate(w(n,7),lw(n))
! call eigch(mat,n,n,n,nv,tol,w,lw,eival,eivec,ier)
! deallocate(w,lw)
! return

  lwork = max(1,2*n-1)
  allocate(work(lwork),rwork(max(1,3*n-2)))
  if(nv.eq.0) then
    call zheev('N','U',n,mat,n,eival,work,lwork,rwork,ier)
  else
    call zheev('V','U',n,mat,n,eival,work,lwork,rwork,ier)
    eivec = mat
  endif
  deallocate(work,rwork)
4 format(i5,1x,3(f6.3,1x,f6.3,4x))

 end subroutine diagonalize
!==========================================================
 subroutine write_eigenvalues(i,dk,kp,eival,n,uio)
 use constants
 implicit none
 integer i,j,n,uio
 integer sig(n)
 real(8), dimension(n) :: eival
 real(8) kp(3),dk,mysqrt

! n = size(eival)
! write(ueiv,3)i-1,kp,(eival(j),j=1,n)
! do j=1,n
!    sig(j) = 1
!    if (eival(j).lt.0) then
!       eival(j) = -eival(j)
!       sig(j) = -1d0
!    endif
! enddo
! write(ueiv,3)i-1,kp,(cnst*sig(j)*sqrt(eival(j)),j=1,n)
 write(uio,3)i-1,dk,kp,(cnst*mysqrt(eival(j)),j=1,n)

3 format(i5,1x,4(1x,f8.3),999(1x,g11.5))
 end subroutine write_eigenvalues
!==========================================================
 subroutine calculate_dos(mx,eival,wkp,mesh,omega,ds)
! use gaussian broadening
 use io2
 use params
 use om_dos
 use constants
 implicit none
 integer i,j,mx,mesh
 real(8) x,wkp(mx),delta,ds(mesh),omega(mesh),eival(mx)  !,cnst

! wkp=1d0/(nkx*nky*nkz)
! write(udos,*)'# wkp,width =',wkp,width

    do i=1,mesh
       ds(i) = 0
       do j=1,mx
          x = (eival(j) - omega(i))/width   ! these are all in cm^-1
!         x = (eival(j) - ene*ene)/width/width
          if ( abs(x) .gt. 5 ) cycle
          ds(i) = ds(i) + delta(x)/width*wkp(j) !/mx
       enddo
    enddo

3 format(i5,9(3x,g11.5))

 end subroutine calculate_dos
!==========================================================
 subroutine calculate_jdos(mx,eival,wkp,mesh,omega,udosj)
! use gaussian broadening to calculate the joint dos
 use io2
 use params
 use om_dos
 use constants
 implicit none
 integer i,j,k,mx,mesh,udosj
 real(8) xp,xm,wkp(mx),delta_g,dsp,dsm,omega(mesh),eival(mx),iself,tmp,one
 complex(8) oc2,omz

! wkp=1d0/(nkx*nky*nkz)
 one = 1d0
 write(udosj,*)'# width =',width
 tmp =0.07      ! the temp width =0.07/cm is equiv to 0.1 K
! then (1+n1+n2)=1 should produce dsp results

! open(udosj+1,file='jdos.discrete',status='unknown')
! do j=1,mx
! do k=j,mx
!    xp= cnst*(sqrt(abs(eival(j)))+sqrt(abs(eival(k))))
!    xm= cnst*(sqrt(abs(eival(j)))-sqrt(abs(eival(k))))
!    write(udosj+1,2) xp,one
!!    if (xm.gt.0)
!    write(udosj+1,2) abs(xm),-one
! enddo
! enddo
! close(udosj+1)

    write(udosj,*)'# i,omega(i),jdsm(i),jdsp(i), occupation-weighted jdos'
    do i=1,mesh
       dsp=0 ; dsm=0
       omz=2*cmplx(omega(i),-etaz)
       do j=1,mx
       do k=1,mx
          xp= (cnst*(sqrt(abs(eival(j)))+sqrt(abs(eival(k)))) - 2*omega(i))   ! these are all in cm^-1
          xm= (cnst*(sqrt(abs(eival(j)))-sqrt(abs(eival(k)))) - 2*omega(i))   ! these are all in cm^-1
          if ( abs(xp) .lt. 6 ) then
             dsp = dsp + delta_g(xp,width)*wkp(j)*wkp(k) !* width*sqrt(2*pi)
          endif
          if ( abs(xm) .lt. 6 ) then
             dsm = dsm + delta_g(xm,width)*wkp(j)*wkp(k) !* width*sqrt(2*pi)
          endif

          iself = aimag(oc2(tmp,omz,j,k))/pi*wkp(j)*wkp(k)*tmp
       enddo
       enddo
       write(udosj,3) i,2*omega(i),dsm,dsp,iself
    enddo
    close(udosj)

2 format(9(3x,g11.5))
3 format(i5,9(3x,g11.5))

 end subroutine calculate_jdos
!==========================================================
 subroutine lifetime_dos(q,la,omq,dosp,dosm)
! use gaussian broadening to calculate the joint dos
! for given q and mode index la, calculates sum_1,2
! delta(1-q1-/+q2)delta(w_q-w1-/+w2)
 use io2
 use params
 use om_dos
 use constants
 use eigen
 use kpoints
 implicit none
 integer i,j,k,l,la,nq,i1,j1,k1,nq1,ip,jp,kp,im,jm,km,nqp,nqm,la1,la2,inside
 real(8) q(3),q1(3),qp(3),qm(3),dosp,dosm,om1,omp,omm,omq,delta_l
 complex(8) oc2,omz

 dosp=0; dosm=0
 do l=1,nkc
    q1=kpc(:,l)
    call get_k_info_cent(q1-shft,NC,nq1,i1,j1,k1,g1,g2,g3,inside)
    if(l.ne.nq1) then
      write(ulog,*)'n1,nq1,inside=',l,nq1,inside
      write(ulog,*)'q1=',q1
      write(ulog,*)'ijk1=',i1,j1,k1
      stop
    endif
    qp=-q-q1
    qm=-q+q1
    call get_k_info_cent(qp-shft,NC,nqp,ip,jp,kp,g1,g2,g3,inside)
    call get_k_info_cent(qm-shft,NC,nqm,im,jm,km,g1,g2,g3,inside)
    do la1=1,ndyn
       om1 = eigenval(la1,nq1)
       do la2=1,ndyn
          omp = eigenval(la2,nqp)
          omm = eigenval(la2,nqm)
          dosp=dosp+delta_l(omq-om1-omp,etaz)
          dosm=dosm+delta_l(omq-om1+omm,etaz)
       enddo
    enddo
 enddo
 end subroutine lifetime_dos
!===========================================================
 subroutine phase_space_lifetime(i1,i2,nk,kp,eival)
 use io2
 use om_dos
 use kpoints
 use eigen
 use constants
 implicit none
 integer i,la,nk,i1,i2
 real(8) q(3),omq,dosp,dosm,kp(3,nk),eival(ndyn,nk)  !,vq(3,ndyn)
 complex(8) eivc(ndyn,ndyn)
 character cik1*4,cik2*4


 write(cik1,'(i4.4)') i1
 write(cik2,'(i4.4)') i2

 open(udos,file='lifetime_dos-'//cik1//'-'//cik2//'.dat ',status='unknown')
 write(udos,*)'# la,q,omega(q,la),jdsm(i),jdsp(i)'

    do i=i1,i2
       q=kibz(:,i)
!       call get_freq(q,ndyn,vq,eivl,eivc)
       do la=1,ndyn
          omq = eival(la,i)
          call lifetime_dos(q,la,omq,dosp,dosm)
          write(udos,3) la,q,omq,dosm/nkc,dosp/nkc
       enddo
    enddo
 close(udos)

3 format(i5,9(3x,g10.4))

 end subroutine phase_space_lifetime
!===========================================================
 subroutine read_v35(ksubibz,ndn)
! this version is used when direct access is needed to v33 for given k_indices and bands
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units are in cm^-1
! this works for k1 in the IBZ and k2 in the FBZ on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! to avoid ambiguity in the case of degenerate eigenstates, and to satisfy symmetry exactly,
! not to mention saving time, we calculate only for 0<q_i<G/2 and complete by symmetry
!
 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use params
 use constants
 implicit none
 integer, intent(in) :: ksubibz(2),ndn
 integer l1,l2,l3,n1,n2,nk1,nk2,nibzsub,lm1,lm2,lm3,np1
! integer np1,np2,np3
 character(4) cibz1,cibz2
 character(6) cnk
 character(132) filname

 write(ulog,*)'calculating V3(k2,k3,l1,l2,l3) ##########################'
 write(cibz1,'(i4.4)') ksubibz(1)
 write(cibz2,'(i4.4)') ksubibz(2)
 write(cnk,'(i6.6)') nkc
 nibzsub=ksubibz(2)-ksubibz(1)+1

 filname=trim(v3path)//'v35-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
 write(ulog,*)'opening the file ',filname
 open(uv3,file=trim(filname),status='OLD' ,FORM='UNFORMATTED')
 read(uv3)nkc,n2,lm1,lm2 !ksubibz(1),ksubibz(2)
 if (n2.ne.ndn) then
    write(ulog,*)'READ_V3: inconsistency in ndn:reading ',n2,' instead of ',ndn
    stop
 endif
 if (lm1.ne.ksubibz(1) .or. lm2.ne.ksubibz(2) ) then
    write(ulog,*)'READ_V3: inconsistency in ksub:reading ',lm1,lm2,' instead of ',ksubibz
    stop
 endif

! v3 = cmplx(0d0,0d0)
! not initializing might be better if some array elements are not properly initialized

 loop1: do n1=ksubibz(1),ksubibz(2)
    np1=n1-ksubibz(1)+1
 loop2: do n2=1,nkc   ! second argument k2 needs to be only in the IFBZ
 do lm1=1 ,ndn
 do lm2=1 ,ndn
!   np1= l1+(nk1-1)*ndn
!   np2= l2+(nk2-1)*ndn
!   if(np1.gt.np2) cycle
 do lm3=1 ,ndn
!   np3= l3+(nk3-1)*ndn
!   if(np2.gt.np3) cycle

    read(uv3)nk1,nk2,l1,l2,l3,v3sq(np1,nk2,l1,l2,l3)

 enddo
 enddo
 enddo
 enddo loop2
 enddo loop1

 close(uv3)

 end subroutine read_v35
!===========================================================
 subroutine calculate_v35(ksubibz,ndn,nib,eivlibz,eivcibz,nk,eival,eivec)
! this version is used when direct access is needed to v33 for given k_indices and bands
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units are in cm^-1
! this works for k1 in the IBZ and k2 in the FBZ on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! to avoid ambiguity in the case of degenerate eigenstates, and to satisfy symmetry exactly,
! not to mention saving time, we calculate only for 0<q_i<G/2 and complete by symmetry
!
 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use params
 use constants
 implicit none
 integer, intent(in) :: ksubibz(2),ndn,nk,nib
 integer l1,l2,l3,n2,n3,al,be,ga,i0,j0,k0,j,k,t,ta1,ta2,ta3
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,np1,np2,np3,n1 ,nibzsub
 real(8) q1(3),q2(3),q3(3),mi,mj,mk,rx2(3),rx3(3),den,eivlibz(ndn,nib),eival(ndn,nk)
 complex(8) eivec(ndn,ndn,nk),xx,eiqr,eivcibz(ndn,ndn,nib)
 character(4) cibz1,cibz2
 character(6) cnk
 character(132) filname

 write(ulog,*)'calculating V3(k2,k3,l1,l2,l3) ##########################'
 write(cibz1,'(i4.4)') ksubibz(1)
 write(cibz2,'(i4.4)') ksubibz(2)
 write(cnk,'(i6.6)') nkc
 nibzsub=ksubibz(2)-ksubibz(1)+1

 filname=trim(v3path)//'v35-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
 write(ulog,*)'opening the file ',filname
 open(uv3,file=trim(filname),status='unknown' ,FORM='UNFORMATTED')
 write(uv3)nkc,ndn,ksubibz(1),ksubibz(2)

! v3 = cmplx(0d0,0d0)
! not initializing might be better if some array elements are not properly initialized

 loop1: do n1=ksubibz(1),ksubibz(2)
    np1=n1-ksubibz(1)+1
    q1 = kibz(:,n1) ! =kpc(:,mapinv(n1))
    call get_k_info_cent(q1,NC,nk1,i2,j2,k2,g1,g2,g3,inside)
    if(mapinv(n1).ne.nk1) then
      write(ulog,*)'mapinv(n1),nk1,inside=',mapinv(n1),nk1,inside
      write(ulog,*)'q1=',q1
      write(ulog,*)'ijk1=',i2,j2,k2
      stop
    endif
    write(*,2)'nibz,nk,kibz=',n1,nk1,q1
    write(ulog,2)'nibz,nk,kibz=',n1,nk1,q1

 loop2: do n2=1,nkc   ! second argument k2 needs to be only in the IFBZ
    q2 = kpc(:,n2)
    call get_k_info_cent(q2-shft,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(n2.ne.nk2) then
      write(ulog,*)'n2,nk2,inside=',n2,nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif

    q3=-q1-q2
    call get_k_info_cent(q3-shft,NC,nk3,i2,j2,k2,g1,g2,g3,inside)

!   write(*,2)'nk123=',nk1,nk2,nk3

 do l1=1 ,ndn
 do l2=1 ,ndn

!   np1= l1+(nk1-1)*ndn
!   np2= l2+(nk2-1)*ndn
!   if(np1.gt.np2) cycle    ! permutation symmetry can not be used that easily
! as the thwo k1 and k2 meshes differ. In case they overlap, in that overlap region one can
! use permiation symmetry which can be headache.

 do l3=1 ,ndn
!   np3= l3+(nk3-1)*ndn
!   if(np2.gt.np3) cycle

    xx = cmplx(0d0,0d0)
    tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop

       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j) - atompos(:,j0)
       rx3(:) = atompos(:,k) - atompos(:,k0)
!  be careful: it has to be the translations R not atompos!
!&      exp(ci*(kp(:,n2) .dot. atompos(:,j) + kp(:,n3) .dot. atompos(:,k))) /  &
       eiqr = exp( ci* ((q2 .dot. rx2) + (q3 .dot. rx3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
!&      eivec(ta1,l1,nk1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
&      eivcibz(ta1,l1,n1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
! make sure for kp at the FBZ boundary it gives (-1) or 1

! n2 is already in the IBZ, n3 is not but for now we are using the full mesh
! and nibz=nkc;  therefore n1 is also in the full zone and is OK, otherwise
! for n1 and n3 it is necessary to rotate the eigenvector from the IBZ to its
! true position in the recirpocal space. (the op_kmatrix will do the job)

    enddo tloop

!   den = sqrt(sqrt(abs(eival(l1,nk1)*eival(l2,nk2)*eival(l3,nk3))))
    den = sqrt(eivlibz(l1,n1)*eival(l2,nk2)*eival(l3,nk3))/cnst/sqrt(cnst)
    if (den.ne.0) then
       xx = xx / den/2d0/sqrt(2d0) * const33
!      write(*,6)'nk_i,la_i,xx=',nk1,nk2,nk3,l1,l2,l3,xx
    else   ! usually acoustic eigenvalues at k=0 are not strictly zero, but a small number
       write(ulog,*)'calculate_v3: denominator=0 !! '
       write(ulog,*)'n1,n2,n3,l123=', n1,n2,nk3,l1,l2,l3
       write(ulog,*)'q1=', q1
       write(ulog,*)'q2=', q2
       write(ulog,*)'q3=', q3
       write(ulog,*)'eivs123=',eivlibz(l1,n1),eival(l2,n2),eival(l3,nk3)
       stop
    endif

    v3sq(np1,nk2,l1,l2,l3)=xx*conjg(xx)
    write(uv3)n1,nk2,l1,l2,l3,v3sq(np1,nk2,l1,l2,l3)

 enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop2
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop1

 close(uv3)

2 format(a,2(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g14.8))
6 format(a,6(i5),9(2x,g14.8))

 end subroutine calculate_v35
!===========================================================
 subroutine calculate_number_of_v3_terms
 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use params
 use constants
 use eigen
 implicit none
 integer i1,i2,i3,l1,l2,l3,k1,k2,k3,inside,n,ndm,isym,narms,kvecop(48),ier,ncmp,mapp(3)
 integer nk1,nk2,nk3,i6,j6,k6,nkc3,j
 real(8) q1(3),q2(3),q3(3),w(3,3),v1(3),v2(3),w1(3),w2(3),w3(3),wdr
 real(8) kvecstar(3,48),primlatt(3,3),v(3),kvecstarp(3,48)
 real(8), allocatable:: triplet(:,:,:),triplets(:,:,:)
 integer, allocatable:: ik1(:),ik2(:),ik3(:)

 nkc3=nkc*nkc*nkc/6
 allocate(triplet(3,3,nkc3),ik1(nkc3),ik2(nkc3),ik3(nkc3))
      
 !mapping is: i=lambda+(ik-1)*ndyn ; inverse mapp is: lambda=1+mod(i-1,ndyn) ; ik=1+(i-lambda)/ndyn
 ndm=ndyn*nkc

! for now we only search for irreducible triplets of kpoints (regardless of the bands)
 n=0
! do i1=1,ndm
!    l1=1+mod(i1-1,ndyn) ; k1=1+(i1-l1)/ndyn ; q1=kpc(:,k1)
! do i2=i1,ndm
!    l2=1+mod(i2-1,ndyn) ; k2=1+(i2-l2)/ndyn ; q2=kpc(:,k2)
! main: do i3=i2,ndm
!    l3=1+mod(i3-1,ndyn) ; k3=1+(i3-l3)/ndyn ; q3=kpc(:,k3)
 do k1=1,ndm
 do k2=k1,ndm
 main: do k3=k2,ndm
    q1=kpc(:,k1); q2=kpc(:,k2); q3=kpc(:,k3)

! they have to add up to zero or G
!   w1=0 ;  call get_k_info_cent(w1,NC,nk1,i6,j6,k6,g1,g2,g3,inside)
!   w2=q1+q2+q3;  call get_k_info_cent(w2,NC,nk2,i6,j6,k6,g1,g2,g3,inside)
!   if (nk1.ne.nk2) cycle main 
    w2=q1+q2+q3; 
    wdr = (w2.dot.rr1)
    if (.not. (wdr-floor(wdr) .myeq. 0d0)) cycle main 
    wdr = (w2.dot.rr2)
    if (.not. (wdr-floor(wdr) .myeq. 0d0)) cycle main 
    wdr = (w2.dot.rr3)
    if (.not. (wdr-floor(wdr) .myeq. 0d0)) cycle main 
    
      narms=0
      symloop: do isym=1,lattpgcount
! apply symmetry operation to k to get v=kstar & find the reduced coordinates of v and store in v2            
        w(:,1)=matmul(op_kmatrix(:,:,isym),q1)
        call get_k_info_cent(w(:,1),NC,nk1,i6,j6,k6,g1,g2,g3,inside)

        w(:,2)=matmul(op_kmatrix(:,:,isym),q2)
        call get_k_info_cent(w(:,2),NC,nk2,i6,j6,k6,g1,g2,g3,inside)

        w(:,3)=matmul(op_kmatrix(:,:,isym),q3)
        call get_k_info_cent(w(:,3),NC,nk3,i6,j6,k6,g1,g2,g3,inside)

! now sort nk1,nk2,nk3 
        call sort3(nk1,nk2,nk3,mapp)  
! mapp(i) gives the index if the ith smallest integer: mapp(1) =nk1 if nk1 is the smallest 
        w1=kpc(:,mapp(1)) ;w2=kpc(:,mapp(2)) ;w3=kpc(:,mapp(3))   
! these are the sorted rotated vectors in case their rotated image goes outside the kpc list

! now compare the new triplet (w1,w2,w3) to the old set generated up to now
        do j=1,n   ! loop over existing triplets
           if (length(w1-triplet(:,1,j)) .myeq. 0d0) then
           if (length(w2-triplet(:,2,j)) .myeq. 0d0) then
           if (length(w3-triplet(:,3,j)) .myeq. 0d0) then ! the triplet has already been created
              cycle symloop
           endif
           endif
           endif
        enddo
     enddo symloop
! if this point reached it means the new triplet was not encountered: so it is new and we add it to the list
        narms=narms+1 ! this is the number of stars of the initial triplet
        n=n+1
        triplet(:,1,n)=w1
        triplet(:,2,n)=w2
        triplet(:,3,n)=w3
        ik1(n)=k1 ; ik2(n)=k2 ; ik3(n)=k3 
        write(567,*)k1,k2,k3,n

 enddo main
 enddo
 enddo
 
 stop








!       call xmatmlt(w1,primlatt,v1,1,3,3,1,3,1)
!       call xmatmlt(w2,primlatt,v2,1,3,3,1,3,1)
!       call xmatmlt(w3,primlatt,v3,1,3,3,1,3,1)
! now check if vi - any_previous_vi is integer (differ by a G vector)
! if so, skip; else store this v as a new star vector
!       do j=1,narms
!       ! subtract previous_v2(=kvecstarp) from v2; result is v3
!         call xvsub(v2,kvecstarp(1,j),v3,3)
!         do k=1,3
!           n=nint(v3(k))
!           if(ncmp(v3(k)-n).ne.0)exit
!           if(k.eq.3)cycle iloop  ! goto next sym_op iff v3=integer
!         enddo
!       enddo
!       narms=narms+1
!       kvecstar(1:3,narms)=v(1:3)
!       kvecstarp(1:3,narms)=v2(1:3)
!       kvecop(narms)=i

 end subroutine calculate_number_of_v3_terms
!===========================================================
 subroutine calculate_v35r(ksubibz,ndn,nib,eivlibz,eivcibz,nk,eival,eivec)
! this version is used when direct access is needed to v33 for given k_indices and bands
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units are in cm^-1
! this works for k1 in the IBZ and k2 in the FBZ on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! to avoid ambiguity in the case of degenerate eigenstates, and to satisfy symmetry exactly,
! not to mention saving time, we calculate only for 0<q_i<G/2 and complete by symmetry
!
 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use params
 use constants
 implicit none
 integer, intent(in) :: ksubibz(2),ndn,nk,nib
 integer l1,l2,l3,n2,n3,al,be,ga,i0,j0,k0,j,k,t,ta1,ta2,ta3
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,np1,np2,np3,n1 ,nibzsub
 real(8) q1(3),q2(3),q3(3),mi,mj,mk,rx2(3),rx3(3),den,eivlibz(ndn,nib),eival(ndn,nk)
 complex(8) eivec(ndn,ndn,nk),xx,eiqr,eivcibz(ndn,ndn,nib)
 character(4) cibz1,cibz2
 character(6) cnk
 character(132) filname


! any triplet ki,li is sorted first 

 write(ulog,*)'calculating V3(k2,k3,l1,l2,l3) ##########################'
 write(cibz1,'(i4.4)') ksubibz(1)
 write(cibz2,'(i4.4)') ksubibz(2)
 write(cnk,'(i6.6)') nkc
 nibzsub=ksubibz(2)-ksubibz(1)+1

 filname=trim(v3path)//'v35-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
 write(ulog,*)'opening the file ',filname
 open(uv3,file=trim(filname),status='unknown' ,FORM='UNFORMATTED')
 write(uv3)nkc,ndn,ksubibz(1),ksubibz(2)

! v3 = cmplx(0d0,0d0)
! not initializing might be better if some array elements are not properly initialized

 loop1: do n1=ksubibz(1),ksubibz(2)
    np1=n1-ksubibz(1)+1
    q1 = kibz(:,n1) ! =kpc(:,mapinv(n1))
    call get_k_info_cent(q1,NC,nk1,i2,j2,k2,g1,g2,g3,inside)
    if(mapinv(n1).ne.nk1) then
      write(ulog,*)'mapinv(n1),nk1,inside=',mapinv(n1),nk1,inside
      write(ulog,*)'q1=',q1
      write(ulog,*)'ijk1=',i2,j2,k2
      stop
    endif
    write(*,2)'nibz,nk,kibz=',n1,nk1,q1
    write(ulog,2)'nibz,nk,kibz=',n1,nk1,q1

 loop2: do n2=1,nkc   ! second argument k2 needs to be only in the IFBZ
    q2 = kpc(:,n2)
    call get_k_info_cent(q2-shft,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(n2.ne.nk2) then
      write(ulog,*)'n2,nk2,inside=',n2,nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif

    q3=-q1-q2
    call get_k_info_cent(q3-shft,NC,nk3,i2,j2,k2,g1,g2,g3,inside)

!   write(*,2)'nk123=',nk1,nk2,nk3

 do l1=1 ,ndn
 do l2=1 ,ndn

!   np1= l1+(nk1-1)*ndn
!   np2= l2+(nk2-1)*ndn
!   if(np1.gt.np2) cycle    ! permutation symmetry can not be used that easily
! as the thwo k1 and k2 meshes differ. In case they overlap, in that overlap region one can
! use permiation symmetry which can be headache.

 do l3=1 ,ndn
!   np3= l3+(nk3-1)*ndn
!   if(np2.gt.np3) cycle

    xx = cmplx(0d0,0d0)
    tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop

       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j) - atompos(:,j0)
       rx3(:) = atompos(:,k) - atompos(:,k0)
!  be careful: it has to be the translations R not atompos!
!&      exp(ci*(kp(:,n2) .dot. atompos(:,j) + kp(:,n3) .dot. atompos(:,k))) /  &
       eiqr = exp( ci* ((q2 .dot. rx2) + (q3 .dot. rx3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
!&      eivec(ta1,l1,nk1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
&      eivcibz(ta1,l1,n1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
! make sure for kp at the FBZ boundary it gives (-1) or 1

! n2 is already in the IBZ, n3 is not but for now we are using the full mesh
! and nibz=nkc;  therefore n1 is also in the full zone and is OK, otherwise
! for n1 and n3 it is necessary to rotate the eigenvector from the IBZ to its
! true position in the recirpocal space. (the op_kmatrix will do the job)

    enddo tloop

!   den = sqrt(sqrt(abs(eival(l1,nk1)*eival(l2,nk2)*eival(l3,nk3))))
    den = sqrt(eivlibz(l1,n1)*eival(l2,nk2)*eival(l3,nk3))/cnst/sqrt(cnst)
    if (den.ne.0) then
       xx = xx / den/2d0/sqrt(2d0) * const33
!      write(*,6)'nk_i,la_i,xx=',nk1,nk2,nk3,l1,l2,l3,xx
    else   ! usually acoustic eigenvalues at k=0 are not strictly zero, but a small number
       write(ulog,*)'calculate_v3: denominator=0 !! '
       write(ulog,*)'n1,n2,n3,l123=', n1,n2,nk3,l1,l2,l3
       write(ulog,*)'q1=', q1
       write(ulog,*)'q2=', q2
       write(ulog,*)'q3=', q3
       write(ulog,*)'eivs123=',eivlibz(l1,n1),eival(l2,n2),eival(l3,nk3)
       stop
    endif

    v3sq(np1,nk2,l1,l2,l3)=xx*conjg(xx)
    write(uv3)n1,nk2,l1,l2,l3,v3sq(np1,nk2,l1,l2,l3)

 enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop2
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop1

 close(uv3)

2 format(a,2(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g14.8))
6 format(a,6(i5),9(2x,g14.8))

 end subroutine calculate_v35r
!===========================================================
 subroutine calculate_v3_fly  (ksubibz,ndn)
! does not include calls to get_k_info, especially for the 3rd q-vector
! this version is used when direct access is needed to v33 for given k_indices and bands
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units are in cm^-1
! this works for k1 in the IBZ and k2 in the FBZ on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! to avoid ambiguity in the case of degenerate eigenstates, and to satisfy symmetry exactly,
! not to mention saving time, we calculate only for 0<q_i<G/2 and complete by symmetry
!
 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use params
 use constants
 use eigen
 implicit none
 integer, intent(in) :: ksubibz(2),ndn
 integer l1,l2,l3,n1,n2,n3,al,be,ga,i0,j0,k0,j,k,t,ta1,ta2,ta3,nibzsub
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside,np1
 real(8) q1(3),q2(3),q3(3),mi,mj,mk,rx2(3),rx3(3),den,eivl3(ndn),vq3(3,ndn)
 complex(8) eivc3(ndn,ndn),xx,eiqr
 character(4) cibz1,cibz2
 character(6) cnk
 character(132) filname

 write(ulog,*)'calculating V3(k2,k3,l1,l2,l3) ##########################'
 write(cibz1,'(i4.4)') ksubibz(1)
 write(cibz2,'(i4.4)') ksubibz(2)
 write(cnk,'(i6.6)') nkc
 nibzsub=ksubibz(2)-ksubibz(1)+1

 filname=trim(v3path)//'v35-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
 write(ulog,*)'opening the file ',filname
 open(uv3,file=trim(filname),status='unknown',FORM='UNFORMATTED')
 write(uv3)nkc,ndn,ksubibz(1),ksubibz(2)

! v3 = cmplx(0d0,0d0)
! not initializing might be better if some array elements are not properly initialized

 loop1: do n1=ksubibz(1),ksubibz(2)
    np1=n1-ksubibz(1)+1
    q1 = kibz(:,n1)   !kpc(:,mapinv(n1))

 loop2: do n2=1,nkc   ! second argument k2 needs to be only in the IFBZ
    q2 = kpc(:,n2)
    q3=-q1-q2

    call get_freq(q3,ndn,vq3,eivl3,eivc3)

 do l1=1 ,ndn
 do l2=1 ,ndn
 do l3=1 ,ndn

    xx = cmplx(0d0,0d0)
    tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop

       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j) - atompos(:,j0)
       rx3(:) = atompos(:,k) - atompos(:,k0)
!  be careful: it has to be the translations R not atompos!
!&      exp(ci*(kp(:,n2) .dot. atompos(:,j) + kp(:,n3) .dot. atompos(:,k))) /  &
       eiqr = exp( ci* ((q2 .dot. rx2) + (q3 .dot. rx3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
&      eivecibz(ta1,l1,n1)*eigenvec(ta2,l2,n2)*eivc3(ta3,l3)
! make sure for kp at the FBZ boundary it gives (-1) or 1

! n2 is already in the IBZ, n3 is not but for now we are using the full mesh
! and nibz=nkc;  therefore n1 is also in the full zone and is OK, otherwise
! for n1 and n3 it is necessary to rotate the eigenvector from the IBZ to its
! true position in the recirpocal space. (the op_kmatrix will do the job)

    enddo tloop

    den = sqrt(sqrt(abs(eivalibz(l1,n1)*eigenval(l2,n2)*eivl3(l3))))
    if (den.ne.0) then
       xx = xx / den/2d0/sqrt(2d0) * const33
!      write(*,6)'nk_i,la_i,xx=',nk1,nk2,nk3,l1,l2,l3,xx
    else   ! usually acoustic eigenvalues at k=0 are not strictly zero, but a small number
       write(ulog,*)'calculate_v3: denominator=0 !! '
       write(ulog,*)'n1,n2,l123=', n1,n2,l1,l2,l3
       write(ulog,*)'q1=', q1
       write(ulog,*)'q2=', q2
       write(ulog,*)'q3=', q3
       write(ulog,*)'eivs123=',eivalibz(l1,n1),eigenval(l2,n2),eivl3(l3)
       stop
    endif

    v3sq(np1,n2,l1,l2,l3)=xx*conjg(xx)
    write(uv3)n1,n2,l1,l2,l3,v3sq(np1,n2,l1,l2,l3)

 enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop2
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop1

 close(uv3)

2 format(a,2(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g14.8))
6 format(a,6(i5),9(2x,g14.8))

 end subroutine calculate_v3_fly
!===========================================================
 subroutine calculate_v3_new(ndn,nk,kp,eival,eivec)
! check correctness of v3 by numerical differentiation of dynmat
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units are in cm^-1
! this works for two k's (k2,k3) in the prim cell on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! this subroutine calculates V3(k1,la1,k2,la2,k3,la3)  with k1=-k2-k3+G  and all ki in the primitive G-cell
! to avoid ambiguity in the case of degenerate eigenstates, and to satisfy symmetry exactly,
! not to mention saving time, we calculate only for 0<q_i<G/2 and complete by symmetry
!
 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ndn,l1,l2,l3,n2,n3,al,be,ga,i0,j0,k0,j,k,t,ta1,ta2,ta3
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,np,nq,ns,indx
 real(8) q1(3),q2(3),q3(3),mi,mj,mk,const,rx2(3),rx3(3),den,eival(ndn,nk),kp(3,nk) !,eivl
 complex(8) eivec(ndn,ndn,nk),xx,eiqr

 const = ee*1d30   ! convert ev/A^3 to SI units (kg s^-2 m^-1)
 const = const *(hbar/uma)**1.5d0  ! that takes care of hbar^1.5 and the 3 masses in the denominator
 const = const /(sqrt(ee*1d20/uma))**1.5d0 ! that takes care of the 3 frequencies in the denominator
 const = const/hbar ! to convert energy(J) to frequency
 const = const /200/pi/c_light  ! to convert frequency to cm^-1
 write(ulog,*)'CALCULATE_V3: to convert v3 to cm^-1, const=',const
! const = ee*1d30*(hbar/uma/(200*pi*c_light))**(1.5d0) ! cnvert v3 to SI if freqs in cm^-1
! const = const/(100*h_plank*c_light) ! to convert from SI to cm^-1
! write(ulog,*)'const2=',const

 write(ulog,*)'writing V3(k2,k3,l1,l2,l3) ##########################'

 v33 = cmplx(0d0,0d0)
 indx=0
 write(ulog,*)'V3: nc(123)=',nc

 loop2: do n2=1,nk   ! second argument k2 needs to be only in the IFBZ
! loop2: do np=1,npos
!   n2=mappos(np)
    q2 = kp(:,n2)
    call get_k_info_cent(q2-shft,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(n2.ne.nk2) then
      write(ulog,*)'n2,nk2,inside=',n2,nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif
!   call get_k_info(-q2,NC,mk2,i2,j2,k2,g1,g2,g3,inside)
!   write(*,*)' nk2,mk2=',nk2,mk2

 loop3: do n3=1 ,nk
! loop3: do nq=1,npos
!   n3=mappos(nq)
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    call get_k_info_cent(q3-shft,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif
!   call get_k_info(-q3,NC,mk3,i2,j2,k2,g1,g2,g3,inside)

    q1 = -q2-q3
    call get_k_info_cent(q1-shft,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
!   aux1=conjg(eigenvec(:,:,nk1))
!   rloop:do nt=1,npos
!      if(nk1.eq.mappos(nt)) then ! nk1 corresponds to a pos-eig
!        aux1=eigenvec(:,:,nk1)
!        exit rloop
!      endif
!   enddo rloop
!   call get_k_info(-q1,NC,mk1,i1,j1,k1,g1,g2,g3,inside)

  write(*,2)'nk123=',nk1,nk2,nk3

 do l2=1 ,ndn
 do l3=1 ,ndn

    np = nk2+(l2-1)*nk
    nq = nk3+(l3-1)*nk
    if(np.gt.nq) cycle

 do l1=1 ,ndn
    ns = nk1+(l1-1)*nk
    if(ns.gt.np) cycle

    xx = cmplx(0d0,0d0)
    tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop

       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j) - atompos(:,j0)
       rx3(:) = atompos(:,k) - atompos(:,k0)
!  be careful: it has to be the translations R not atompos!
!&      exp(ci*(kp(:,n2) .dot. atompos(:,j) + kp(:,n3) .dot. atompos(:,k))) /  &
       eiqr = exp( ci* ((q2 .dot. rx2) + (q3 .dot. rx3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
&      eivec(ta1,l1,nk1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
! make sure for kp at the FBZ boundary it gives (-1) or 1

! n2 is already in the IBZ, n3 is not but for now we are using the full mesh
! and nibz=nkc;  therefore n1 is also in the full zone and is OK, otherwise
! for n1 and n3 it is necessary to rotate the eigenvector from the IBZ to its
! true position in the recirpocal space. (the op_kmatrix will do the job)

    enddo tloop

    den = sqrt(sqrt(abs(eival(l1,nk1)*eival(l2,nk2)*eival(l3,nk3))))
    if (den.ne.0) then
       xx = xx / den/2d0/sqrt(2d0) * const
!      write(*,6)'nk_i,la_i,xx=',nk1,nk2,nk3,l1,l2,l3,xx
    else   ! usually acoustic eigenvalues at k=0 are not strictly zero, but a small number
       write(ulog,*)'calculate_v3: denominator=0 !! '
       write(ulog,*)'n1,n2,n3,l123=', nk1,n2,n3,l1,l2,l3
       write(ulog,*)'q1=', q1
       write(ulog,*)'q2=', q2
       write(ulog,*)'q3=', q3
       write(ulog,*)'eivs123=',eival(l2,n2),eival(l3,nk3),eival(l1,nk1)
       stop
    endif

    if (abs(xx).gt.v3_threshold) then
       indx=indx+1
       v33(indx)=xx
       nq1(indx)= nk1
       nq2(indx)= nk2
       nq3(indx)= nk3
       la1(indx)= l1
       la2(indx)= l2
       la3(indx)= l3
    endif
 enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop3
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop2

 nv3 = indx
 write(ulog,*)' V33: total size of this array is =',nv3

 if( readv3.eq.0) then
! k1 if( writev3.eq.1) then
   open(uv3,file='v33.dat',status='unknown',FORM='UNFORMATTED')
!  open(uv3,file='v33.dat',status='unknown')
   write(uv3)nv3
!  write(uv3,*)nv3
   do j=1 ,nv3
      write(uv3) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
!     write(uv3,6)' ',nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
   enddo
   close(uv3)
 endif
! v3 should be of the order of sub eV (v3=1 eV if phi_3=1000 eV/A^3)

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))

 end subroutine calculate_v3_new
!===========================================================
 subroutine calculate_v3_ibz(ndn,ni,ki,nk,kp,eival,eivec)
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units for V3 are in cm^-1
! this works for two k's (k2,k3) in the prim cell on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! this subroutine calculates V3(k1,la1,k2,la2,k3,la3)  with k1=-k2-k3+G  and all ki in the primitive G-cell
! with the restriction of k2 in the IRREDUCIBLE BZ defined by ki(3,ni) included in kp(3,nk)
!
 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ni,ndn,l1,l2,l3,n2,n3,al,be,ga,i0,j0,k0,j,k,t,ta1,ta2,ta3
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,np,nq,ns,indx
 real(8) q1(3),q2(3),q3(3),mi,mj,mk,const,rx2(3),rx3(3),den,eival(ndn,nk)
 real(8) ki(3,ni),kp(3,nk)
 complex(8) eivec(ndn,ndn,nk),xx,eiqr

!*********************************************************************
!  NEED TO DEFINE EIVAL AND EIVEC OVER THE WHOLE ZONE
!*********************************************************************

 const = 1d30   ! convert ev/A^3 to eV ######## SI units (kg s^-2 m^-1)
 const = const *(hbar/uma)**1.5d0  ! that takes care of hbar^1.5 and the 3 masses in the denominator
 const = const /(sqrt(ee*1d20/uma))**1.5d0 ! that takes care of the 3 frequencies in the denominator
 const = const * ee  ! to convert energy from eV to energy Joules
 const = const/h_plank ! to convert energy(J) to frequency(Hz)
 const = const /100/c_light  ! to convert frequency(Hz) to cm^-1
 write(ulog,*)'CALCULATE_V3: to convert v3 to cm^-1, const=',const
 write(ulog,*)'writing V3(k2,k3,l1,l2,l3) ##########################'

 v33 = cmplx(0d0,0d0)
 indx=0
 write(ulog,*)'V3: nc(123)=',nc

!************************************************
! second argument q2 needs to be only in the IFBZ
!************************************************
 loop2: do n2=1,ni
    q2 = ki(:,n2)  ! should be = kp(:,mapinv(n2))
    call get_k_info_cent(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(mapinv(n2).ne.nk2) then
      write(ulog,*)'n2,mapinv(n2),nk2,inside=',n2,mapinv(n2),nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif
!   call get_k_info(-q2,NC,mk2,i2,j2,k2,g1,g2,g3,inside)
!   write(*,*)' nk2,mk2=',nk2,mk2

 loop3: do n3=1 ,nk
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    call get_k_info_cent(q3-shft,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif
!   call get_k_info(-q3,NC,mk3,i2,j2,k2,g1,g2,g3,inside)

    q1 = -q2-q3

    call get_k_info_cent(q1-shft,NC,nk1,i1,j1,k1,g1,g2,g3,inside)

  write(*,2)'nk123=',nk1,nk2,nk3

 do l2=1 ,ndn
 do l3=1 ,ndn

    np = nk2+(l2-1)*nk ! ni
    nq = nk3+(l3-1)*nk
!   if(np.gt.nq) cycle

 do l1=1 ,ndn
    ns = nk1+(l1-1)*nk
!   if(ns.gt.np) cycle

    xx = cmplx(0d0,0d0)
    tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop

       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j) - atompos(:,j0)
       rx3(:) = atompos(:,k) - atompos(:,k0)
!  be careful: it has to be the translations R not atompos!
!&      exp(ci*(kp(:,n2) .dot. atompos(:,j) + kp(:,n3) .dot. atompos(:,k))) /  &
       eiqr = exp( ci* ((q2 .dot. rx2) + (q3 .dot. rx3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
&      eivec(ta1,l1,nk1)*eivec(ta2,l2,nk2)*eivec(ta3,l3,nk3)
! make sure for kp at the FBZ boundary it gives (-1) or 1

! n2 is already in the IBZ, n3 is not
! therefore n1 is also in the full zone and is OK, otherwise
! for n1 and n3 it is necessary to rotate the eigenvector from the IBZ to its
! true position in the recirpocal space. (the op_kmatrix will do the job)

    enddo tloop

    den = sqrt(sqrt(abs(eival(l1,nk1)*eival(l2,nk2)*eival(l3,nk3))))
    if (den.ne.0) then
       xx = xx / den/2d0/sqrt(2d0) * const
!      write(*,6)'nk_i,la_i,xx=',nk1,nk2,nk3,l1,l2,l3,xx
    else   ! usually acoustic eigenvalues at k=0 are not strictly zero, but a small number
       write(ulog,*)'calculate_v3: denominator=0 !! '
       write(ulog,*)'n1,n2,n3,l123=', nk1,nk2,nk3,l1,l2,l3
       write(ulog,*)'q1=', q1
       write(ulog,*)'q2=', q2
       write(ulog,*)'q3=', q3
       write(ulog,*)'eivs123=',eival(l1,nk1),eival(l2,nk2),eival(l3,nk3)
       stop
    endif

    if (abs(xx).gt.v3_threshold) then
       indx=indx+1
       v33(indx)=xx
       nq1(indx)= nk1
       nq2(indx)= nk2
       nq3(indx)= nk3
       la1(indx)= l1
       la2(indx)= l2
       la3(indx)= l3
       if (mod(indx,1000).eq.1) write(*,*) 'indx,v33=',indx,xx
    endif
 enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop3
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop2

 nv3 = indx
 write(ulog,*)' V33: total size of this array is =',nv3

 if( readv3.eq.0) then
! k1 if( writev3.eq.1) then
   open(uv3,file='v33.dat',status='unknown',FORM='UNFORMATTED')
   write(uv3)nv3
   do j=1 ,nv3
      write(uv3) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
!     write(uv3,7) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
   enddo
   close(uv3)
 endif
! v3 should be of the order of sub eV (v3=1 eV if phi_3=1000 eV/A^3)

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))

 end subroutine calculate_v3_ibz
!============================================================
 subroutine read_v3_new_unformatted_k1(ksubs,nv33split,nk)
! this is used for the new version of split mode, with the full index (la,k) of states
! it reads the square of the matrix elements
 use phi3
 use io2
 use lattice
 implicit none
 integer, intent(in)::ksubs(2),nv33split
 integer j,unt,k1,k2,nv33,nk,ndn
 character(4) cibz1,cibz2
 character(6) cnk

 unt = 111
 write(cibz1,'(i4.4)') ksubs(1)
 write(cibz2,'(i4.4)') ksubs(2)
 write(cnk,'(i6.6)') nk
 open(unt,file='v33-'//cnk//'-'//cibz1//'-'//cibz2//'.dat',status='old',form='UNFORMATTED')

   read(unt,end=99)nv33,ndn,k1,k2
   write(ulog,*)' OPENING v33.dat, reading ',nv33,k1,k2

! check consistency
   if( nv33.ne.nv33split .or. k1.ne.ksubs(1) .or. k2.ne.ksubs(2) ) then
      write(ulog,*)'READ_V3_NEW_UNFORMATTED_K1: inconsistency in header'
      write(ulog,*)'nv3,ksubs is supposed to be:',nv33split,ksubs
      stop
   endif

   v33sq= 0d0
   do j=1 ,nv33split
      read(unt,end=99) nq1(j),nq2(j),nq3(j),v33sq(j)
!     read(unt,*,end=99) nq1(j),nq2(j) !,l1(j),l2(j),l3(j),v33sq(j)
   enddo

 close(unt)
 return

99 write(*,*)'v33 file end reached at line j=',j
   write(ulog,*)'READ_V3_NEW_UNFORMATTED_k1:v33 file end reached at line j=',j
   stop

 end subroutine read_v3_new_unformatted_k1
!============================================================
 subroutine calculate_kappa(temp,kappau,kappat) !,kp,nkc)
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 implicit none
 integer w,ik,la !,nbmx,nbq(nbkmax),npar  ! ,nkp
! integer, allocatable:: nb(:,:)
 real(8) temp,kappau,kappat,omg,nbw,ckl,nbe,v23,gammau,gammat  !,vk(3)grad_intp(3),
 real(8) kcut,q(3)
 complex(8) uself,nself

! here one can afford to take a much finer mesh and interpolate or use the MP method

! need to define kcut for the fine mesh
 kcut = 2.5*length(kibz(:,1)-kibz(:,2))
 write(ulog,*)' KAPPA: kcut for fine mesh=',kcut

 kappat = 0; kappau = 0
 kloop: do ik=1,nibz
    write(*,*)' kloop in kappa =',ik,nbz_fine

    q = kibz(:,ik)

!********  do something about this

!    call find_nearest_kpoints3(nibz_coarse,kibz,q,nbmx,nbq,kcut,npar)

!********  do something about this

    lloop: do la=1,ndyn
       omg = sqrt(abs(eivalibz(la,ik))) * cnst     ! in same units as wmax and temp
! use the on-shell frequencies : find the corresponding w
       w = 1+int((omg/wmax)*wmesh)
       nbw = nbe(omg,temp,classical)
       ckl = (omg/temp)**2 * nbw*(nbw+1) * k_b      ! heat capacity per mode
       v23 = (velocibz(:,la,ik) .dot. velocibz(:,la,ik))/3
!      call function_self(q,la,omg,temp,nbmx,nbq,nself,uself,npar)
       call function_self_new(q,la,omg,temp,nself,uself)
!      uself = 0 ; nself = 0
       gammau = aimag(uself)
       gammat = aimag(uself+nself)
!      kappa = kappa + v23 /(2*gammak * c_light*100 + 1/tau0) *  &
!&                   ckl * (2*pi*c_light)**2
       kappat = kappat + v23 * ckl /(2*gammat*100 + 1/tau0/c_light) *  c_light
       kappau = kappau + v23 * ckl /(2*gammau*100 + 1/tau0/c_light) *  c_light
    enddo lloop
 enddo kloop
 kappau = kappau/nbz_fine*volume_g*1d30
 kappat = kappat/nbz_fine*volume_g*1d30

! deallocate(nb)
 end subroutine calculate_kappa
!===========================================================
 subroutine find_nearest_kpoints(nx,ny,nz,q,nbq,nbs,nbmx)
! this is very fast but works only for a cubic mesh
 use constants
! use kpoints
 use lattice
 use params
 implicit none
 real(8) q(3),a1,a2,a3,ep
 integer n1,n2,n3,nbmx,na1,na2,na3,e1,e2,e3,indexn,nx,ny,nz
 integer nbq(nbmx),nbs(nbmx)

 ep = 5*tolerance
 q = q - ep
! write q= ai*gi
 a1 = (q .dot. r1)/(2*pi)
 a2 = (q .dot. r2)/(2*pi)
 a3 = (q .dot. r3)/(2*pi)
 q = q + ep
! now find the 8 nearest kpoints to it
 n1 = floor(a1*nx)
 n2 = floor(a2*ny)
 n3 = floor(a3*nz)

! these are the indices of the 8 neighbors : kp(:,nbq(j)), j=1,8
 nbq(1) = indexn(n1  ,n2  ,n3  ,nx,ny,nz)
 nbq(2) = indexn(n1  ,n2  ,n3+1,nx,ny,nz)
 nbq(3) = indexn(n1  ,n2+1,n3  ,nx,ny,nz)
 nbq(4) = indexn(n1  ,n2+1,n3+1,nx,ny,nz)
 nbq(5) = indexn(n1+1,n2  ,n3  ,nx,ny,nz)
 nbq(6) = indexn(n1+1,n2  ,n3+1,nx,ny,nz)
 nbq(7) = indexn(n1+1,n2+1,n3  ,nx,ny,nz)
 nbq(8) = indexn(n1+1,n2+1,n3+1,nx,ny,nz)

! along each of the 3 directions find whether ai*Ni is closer to ni or ni+1
! (ei=+1;nai=ni+1) if ni+1 is chosen and (-1,ni) if ni is chosen
 if(a1*nx-n1 .gt. n1+1-a1*nx) then  ! include ni+2
   e1 = 1 ; na1 = n1+1
 else
   e1 =-1 ; na1 = n1
 endif
 if(a2*ny-n2 .gt. n2+1-a2*ny) then  ! include ni+2
   e2 = 1 ; na2 = n2+1
 else
   e2 =-1 ; na2 = n2
 endif
 if(a3*nz-n3 .gt. n3+1-a3*nz) then  ! include ni+2
   e3 = 1 ; na3 = n3+1
 else
   e3 =-1 ; na3 = n3
 endif

 nbq(9) = indexn(na1+e1,na2   ,na3   ,nx,ny,nz)
 nbq(10)= indexn(na1   ,na2+e2,na3   ,nx,ny,nz)
 nbq(11)= indexn(na1   ,na2   ,na3+e3,nx,ny,nz)
 nbq(12)= indexn(na1   ,na2+e2,na3+e3,nx,ny,nz)
 nbq(13)= indexn(na1+e1,na2   ,na3+e3,nx,ny,nz)
 nbq(14)= indexn(na1+e1,na2+e2,na3   ,nx,ny,nz)
 nbq(15)= indexn(na1+e1,na2+e2,na3+e3,nx,ny,nz)

 nbs(:)=1

 end subroutine find_nearest_kpoints
!===========================================================
 subroutine function_self_old(q,la,omega,temp,nself,uself)
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer nq,n1,l1,l2,iq,i3,jq,j3,kq,m3,inside,nk1,ik,jk,kk,nk2,mk2
 real(8) om1,om2,k1(3),k2(3),nb1,nb2,nbe,eta,tk,etas,tr,ti,w12,sr,si,vr,vi,term
 complex(8) omz,sumz,w33,self

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
 etas = eta(omega,tk)
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! om1 = sqrt(abs(eigenval(la,nq))) * cnst
! etaz= (deltak*vgmax)**2/omega !min(etaz,100/omega)
  write(ulog,4)'SELF: omega,etaz,(etas)=',omega,etaz,etas
  omz = cmplx(omega, etaz)   ! this is in cm^-1
! write(ulog,*)'FUNCTION_SELF: OMZ=',omz
 nself=0 ; uself=0

 ti=2*etaz*omega
!     write(117,4)' ======================== NEW KPOINT ', q

 k1loop: do n1=1,nkc
    k1(:)= kpc(:,n1)
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,inside)
    if (nk1.ne.n1) then
       write(ulog,*) 'Function self: nk1 .ne. n1 ',nk1,n1
       stop
    endif
!   mp1 = mapibz(n1)

    k2(:)=-kpc(:,nq)-k1(:)   !these are the normal processes if k2 is within the FBZ
    call get_k_info(k2,nc,nk2,i3,j3,m3,g1,g2,g3,inside)
!   mp2 = mapibz(nk2)
    call get_k_info(-k2,nc,mk2,i3,j3,m3,g1,g2,g3,inside) ! take inside from this line

    sumz = 0
    bnd1loop: do l1=1,ndyn
    bnd2loop: do l2=1,ndyn

!     om1 = sqrt(abs(eigenval(l1,mp1))) * cnst ; nb1=nbe(om1,temp)
!     om2 = sqrt(abs(eigenval(l2,mp2))) * cnst ; nb2=nbe(om2,temp)
      om1 = sqrt(abs(eigenval(l1,nk1))) * cnst ; nb1=nbe(om1,temp,classical)
      om2 = sqrt(abs(eigenval(l2,nk2))) * cnst ; nb2=nbe(om2,temp,classical)

      w33 = v3(nk1,nk2,la,l1,l2)
      if (abs(w33) .lt. v3_threshold)  cycle bnd2loop
!     omz = cmplx(omega, eta(om1,tk)+eta(om2,tk)+eta(omega,tk))   ! this is in cm^-1

      w12=(om1+om2)
      tr=(w12*w12-omega*omega+etaz*etaz)   ! also try it without etaz**2
      term= 2*w12 /(tr*tr+ti*ti)
      sr=term*tr
      si=term*ti

      w12=(om1-om2)
      tr=(w12*w12-omega*omega+etaz*etaz)
      term= 2*w12 /(tr*tr+ti*ti)
      vr=term*tr
      vi=term*ti
      term = w33*conjg(w33)
      sumz = sumz - term * cmplx((nb2+nb1+1)*sr+(nb2-nb1)*vr,(nb2+nb1+1)*si+(nb2-nb1)*vi)
! &              ((nb2+nb1+1)*(1d0/(om1+om2-omz)+ 1d0/(om1+om2+omz)) +  &
! &               (nb2-nb1)  *(1d0/(om1-om2-omz)+ 1d0/(om1-om2+omz)))
!&             ((nb2+nb1+1)*cmplx(sr,si) +  (nb2-nb1)  *cmplx(vr,vi))

!     write(117,4)' ',term,(nb2+nb1+1)*sr+(nb2-nb1)*vr,(nb2+nb1+1)*si+(nb2-nb1)*vi,sumz
    enddo bnd2loop
    enddo bnd1loop

    if(inside.eq.1) then                   ! normal process
        nself = nself + sumz
    else                                ! umklappa process
        uself = uself + sumz
    endif

!    call check_inside_fbz(k2,g1,g2,g3,inside)
!    call send_to_primcell(-q-k1,kk)
!    if (length(kk+q+k1).lt.1d-4) inside=.True.

 enddo k1loop

  nself = nself /nkc /2 *2*pi
  uself = uself /nkc /2 *2*pi
  self = nself + uself
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))

 end subroutine function_self_old
!============================================================
 subroutine function_self_new(q,la,omega,temp,nself,uself)
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1
 integer, save :: cntr
 real(8) om1,eta,tk,term,k1(3),k2(3),q1(3)
 complex(8) omz,self,ocfunc

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! om1 = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 write(ulog,5)'SELF_NEW: last cntr, omega,etaz=',cntr,omega,etaz
! omz = cmplx(omega,-etaz)   ! this is in cm^-1
 nself=0 ; uself=0;
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 jq = nq+(la-1)*nkc
 cntr = 0
 nv3loop: do j=1,nv3
!  q must be the 2nd argument in v33
    if (nq.ne.nq2(j)) cycle nv3loop
    if (la.ne.la2(j)) cycle nv3loop
    if (abs(v33(j)).lt.v3_threshold) cycle nv3loop
    term = v33(j)*conjg(v33(j))

    jk1 = nq1(j)+(la1(j)-1)*nkc
    jk2 = nq2(j)+(la2(j)-1)*nkc
    jk3 = nq3(j)+(la3(j)-1)*nkc
    if(jk2 .ne. jq) print*,'ERRRRRRRRRRORRRRRRRRRRE : jk2 ne jq ',jk2,jq

    k2=kpc(:,nq3(j))
!   q1=kpc(:,nq1(j))
    k1 = -q-k2
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,inside)

!   if(mod(j,36).eq.12) write(ulog,7)'k1,k2,k3=',nk1,nq,nq3(j),inside,k1,q,k2

    if(nk1.ne.nq1(j)) then
      write(ulog,3)'SELF_NEW: nk1.ne.nq1 ',nk1,nq1(j)
      write(ulog,3)'Check consistency between nki and those of v33.dat'
      write(ulog,3)'SELF_NEW: j,la,q,om=',j,la,q,omega
      write(ulog,4)'SELF_NEW: k1,q,k2=',k1,q,k2
      write(ulog,*)'SELF_NEW: nqi    =',nq1(j),nq2(j),nq3(j)
      write(ulog,4)'SELF_NEW:kpc(nqi)=',kpc(:,nq1(j)),kpc(:,nq2(j)),kpc(:,nq3(j))
      write(ulog,4)'SELF_NEW:la,la2(j)=',la,la2(j)
      write(*,3)'SELF_NEW: nk1.ne.nq1 ',nk1,nq1(j)
      write(*,3)'SELF_NEW: j,la,q,om=',j,la,q,omega
      stop
    endif

!    if (jq .eq. jk1) then   ! q is the 1st argument in v33
!       self = ocfunc(temp,omz,jk2,jk3)
!       sumz = sumz + term * self
!!      call get_k_info(kpc(:,nq1(j)),nc,nq,ik,jk,kk,g1,g2,g3,inside)
!    elseif (jq .eq. jk2) then  ! q is the 2nd argument in v33
!       self = ocfunc(temp,omz,jk1,jk3)
!       sumz = sumz + term * self
!!      call get_k_info(kpc(:,nq2(j)),nc,nq,ik,jk,kk,g1,g2,g3,inside)
!    elseif (jq .eq. jk3) then  ! q is the 3rd argument in v33
!       self = ocfunc(temp,omz,jk2,jk1)
!       sumz = sumz + term * self
!      call get_k_info(kpc(:,nq3(j)),nc,nq,ik,jk,kk,g1,g2,g3,inside)
!    else
!      cycle nv3loop
!    endif
!   write(118,4)' ',term,real(self),aimag(self),sumz
       cntr = cntr+1
       self = ocfunc(temp,omz,la1(j),nq1(j),la3(j),nq3(j))  !jk1,jk3)
       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif
!   endif
 enddo nv3loop
! if (cntr.ne.nkc) then
!    write(ulog,*)' K_WEIGHT INCONSISTENCY IN SELF_NEW: cntr.ne.nkc ',cntr,nkc
! endif
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_new
!============================================================
 subroutine function_self_tetra(nqs,q,la,omega,temp,nself,uself)
! the tetrahedron method requires direct access to v3 for given k,la indices.
! we therefore assume v3sq is calculated and used
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the 2nd momentum k2 is on the already-generated kmesh and
! is the second argument of V33(q,k2,k3) where k2 covers the FBZ and k3=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 use tetrahedron
 implicit none
 integer, intent(in) :: la,nqs
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,nk1,nk2,nk3,l1,l2,l3
 integer, save :: cntr
 real(8) om2,nb2,om3,nb3,nbe,resm,resp,eta,tk,k1(3),k3(3),q1(3),q3(3)
 complex(8) omz,self,ocfunc
 real(8), allocatable :: eivp(:),eivm(:),funcp(:),funcm(:),junkp(:),junkm(:)
 integer, allocatable:: ins(:)

 allocate(eivp(nkc),funcp(nkc),eivm(nkc),funcm(nkc),junkp(nkc),junkm(nkc),ins(nkc))

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
! call get_k_info_cent(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! jq = (nq-1)*ndyn+la  !nq+(la-1)*nkc

! omq = sqrt(abs(eigenval(la,nq))) * cnst
! write(ulog,5)'SELF_NEW: last cntr, omega,etaz=',cntr,omega,etaz
 nself=0 ; uself=0;
!cntr = 0

 do l2=1,ndyn
 do l3=1,ndyn

    do nk2=1,nkc
       k3=-q-kpc(:,nk2)
       call get_k_info_cent(k3-shft,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
       ins(nk2)=inside
       om3 = eigenval(l3,nk3) ; nb3=nbe(om3,temp,classical)
       om2 = eigenval(l2,nk2) ; nb2=nbe(om2,temp,classical)

! define the arguments of delta(omega-om1-om3)
       eivp(nk2)=om2+om3
       eivm(nk2)=om2-om3
! define the weighting function
       funcp(nk2)=v3sq(nqs,nk2,la,l2,l3)*(1+nb2+nb3)
       funcm(nk2)=v3sq(nqs,nk2,la,l2,l3)*(nb3-nb2)
    enddo

    call tet_sum(omega,nkc,eivp,funcp,resp,junkp)  ! res= sum_k func(k)*delta(omega-eiv(k))
    call tet_sum(omega,nkc,eivm,funcm,resm,junkm)  ! res= sum_k func(k)*delta(omega-eiv(k))

    do nk2=1,nkc
!   cntr = cntr+1
       if(ins(nk2).eq.1) then                ! normal process
!        nself = nself + cmplx(0,resp)/2
         nself = nself + cmplx(0,junkp(nk2))/2
       else                                ! umklapp process
!        uself = uself + cmplx(0,res)/2
         uself = uself + cmplx(0,junkp(nk2))/2
       endif

!   cntr = cntr+1
       if(ins(nk2).eq.1) then                ! normal process
!        nself = nself + cmplx(0,resm)
         nself = nself + cmplx(0,junkm(nk2))
       else                                ! umklapp process
!        uself = uself + cmplx(0,res)
         nself = nself + cmplx(0,junkm(nk2))
       endif
    enddo

 enddo
 enddo

 nself = nself  *2*pi*pi
 uself = uself  *2*pi*pi
! self  = nself + uself

 deallocate(eivp,eivm,funcp,funcm,junkp,junkm,ins)

! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_tetra
!============================================================
 subroutine function_self_35(nqs,q,la,omega,temp,nself,uself)
! this version uses the 5-component array v33_squared
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kibz mesh and
! the second argument of V33(q,k2,k3) k2 covers the FBZ and similarly k3=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la,nqs
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,n2,nk3,l2,l3
 integer, save :: cntr
 real(8) om1,eta,tk,v32,k2(3),k3(3),v33sqr
 complex(8) omz,self,ocfunc

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm

! input q must be a IBZ vector, nq is its eivalibz index
! call get_k_info_cent(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)

 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 write(ulog,5)'SELF_35: nqs, omega,etaz=',nqs,omega,etaz
 v32 = v3_threshold*v3_threshold
 nself=0 ; uself=0;
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 cntr = 0

 k2loop: do n2=1,nkc

    k2=kpc(:,n2)
    k3=-q-k2
    call get_k_info_cent(k3-shft,nc,nk3,iq,jq,kq,g1,g2,g3,inside)
!   write(ulog,7)'ins,nq,n2,nk3,q,k2,k3=',inside,nqs,n2,nk3,q,k2,k3

    do l2=1,ndyn
    do l3=1,ndyn

       v33sqr = v3sq(nqs,n2,la,l2,l3)
       if(v33sqr.lt.v32) cycle
       self = ocfunc(temp,omz,l2,n2,l3,nk3)
! write(ulog,*)'nqs,n2,nk3,l2,l3=',nqs,n2,nk3,l2,l3
       if(inside.eq.1) then                ! normal process
          nself = nself + v33sqr * self
       else                                ! umklapp process
          uself = uself + v33sqr * self
       endif
    enddo
    enddo
 enddo k2loop
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i4,9(2x,3(1x,f7.3)))

 end subroutine function_self_35
!============================================================
 subroutine function_self_new_sq(q,la,omega,temp,nself,uself)
! this version uses the general index in v33sq
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1,nk3,l1,l2,l3
 integer, save :: cntr
 real(8) om1,eta,tk,v32,k2(3),k3(3),q1(3) !,sh(3)
 complex(8) omz,self,ocfunc

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
! sh=(-0.5d0)*(g1+g2+g3)

! input q must be a IBZ vector
 call get_k_info_cent(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
 jq = (nq-1)*ndyn+la  !nq+(la-1)*nkc

! om1 = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 write(ulog,5)'SELF_NEW: last cntr, omega,etaz=',cntr,omega,etaz
! omz = cmplx(omega,-etaz)   ! this is in cm^-1
 v32 = v3_threshold*v3_threshold
 nself=0 ; uself=0;
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 cntr = 0
 nv3loop: do j=1,nv3split
!  q must be the 2nd argument in v33
    if ( jq .ne. nq1(j) ) cycle nv3loop
    if ( v33sq(j) < v32) cycle nv3loop

    l3=1+mod(nq3(j)-1,ndyn) ; jk3=1+(nq3(j)-l3)/ndyn
    l2=1+mod(nq2(j)-1,ndyn) ; jk2=1+(nq2(j)-l2)/ndyn
    l1=1+mod(nq1(j)-1,ndyn) ; jk1=1+(nq1(j)-l1)/ndyn  ! this is (q,la):
    k2=kpc(:,jk2)
    k3 = -q-k2
    call get_k_info_cent(k3-shft,nc,nk3,ik,jk,kk,g1,g2,g3,inside)

    if(jk1.ne.nq) then
      write(ulog,3)'SELF_NEW: nq.ne.jk1 ',nq,jk1
      stop
    endif
    if(nk3.ne.jk3) then
      write(ulog,3)'SELF_NEW: nk3.ne.jk3 ',nk3,jk3
      write(ulog,3)'Check consistency between nki and those of v33.dat'
      write(ulog,3)'SELF_NEW:j,l,q,om=',j,la,q,omega
      write(ulog,4)'SELF_NEW: q,k2,k3=',q,k2,k3
      write(ulog,*)'SELF_NEW: nqi    =',nq1(j),nq2(j),nq3(j)
      write(ulog,*)'SELF_NEW:l1,l2,l3=',l1,l2,l3
      write(ulog,4)'SELF_NEW:kpc(jki)=',kpc(:,jk1),kpc(:,jk2),kpc(:,jk3)
!      write(ulog,4)'SELF_NEW:la,la2(j)=',la,la2(j)
      write(*,3)'SELF_NEW: nk3.ne.jk3 ',nk3,jk3
      write(*,3)'SELF_NEW: j,la,q,om =',j,la,q,omega
      stop
    endif

    cntr = cntr+1
    self = ocfunc(temp,omz,l2,jk2,l3,jk3)
    if(inside.eq.1) then                ! normal process
        nself = nself + v33sq(j) * self
    else                                ! umklapp process
        uself = uself + v33sq(j) * self
    endif
 enddo nv3loop
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_new_sq
!============================================================
 subroutine function_self_sc(q,la,omega,temp,nself,uself)
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1,cnt2
 integer, save :: cntr
 real(8) om1,eta,tk,term,k1(3),k2(3),q1(3),error
 complex(8) omz,self,oldself,ocfunc

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! om1 = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 write(ulog,5)'SELF_SC: last cntr, omega,etaz=',cntr,omega,etaz
! omz = cmplx(omega,-etaz)   ! this is in cm^-1
 error=1d9
 self=cmplx(0,etaz)

cnt2=0
do while(error.gt.etaz/50 .and. cnt2.lt.50)

 cnt2=cnt2+1
 nself=0 ; uself=0; oldself=self
 eta = max(etaz, aimag(self))   ! this is in cm^-1
 omz = cmplx(omega,-eta)    ! this is in cm^-1
!omz = cmplx(omega+real(self),-eta)    ! this is in cm^-1
!omz = omega+conjg(self)
 jq = nq+(la-1)*nkc
 cntr = 0
 nv3loop: do j=1,nv3
!  q must be the 2nd argument in v33
    if (nq.ne.nq2(j)) cycle nv3loop
    if (la.ne.la2(j)) cycle nv3loop
    if (abs(v33(j)).lt.v3_threshold) cycle nv3loop

    term = v33(j)*conjg(v33(j))

    jk1 = nq1(j)+(la1(j)-1)*nkc
    jk2 = nq2(j)+(la2(j)-1)*nkc
    jk3 = nq3(j)+(la3(j)-1)*nkc
    if(jk2 .ne. jq) print*,'ERRRRRRRRRRORRRRRRRRRRE : jk2 ne jq ',jk2,jq

    k2=kpc(:,nq3(j))
    k1 = -q-k2
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,inside)

    if(nk1.ne.nq1(j)) then
      write(ulog,3)'SELF_SC: nk1.ne.nq1 ',nk1,nq1(j)
      write(ulog,3)'Check consistency between nki and those of v33.dat'
      write(ulog,3)'SELF_SC: j,la,q,om=',j,la,q,omega
      write(ulog,4)'SELF_SC: k1,q,k2=',k1,q,k2
      write(ulog,*)'SELF_SC: nqi    =',nq1(j),nq2(j),nq3(j)
      write(ulog,4)'SELF_SC:kpc(nqi)=',kpc(:,nq1(j)),kpc(:,nq2(j)),kpc(:,nq3(j))
      write(ulog,4)'SELF_SC:la,la2(j)=',la,la2(j)
      write(*,3)'SELF_SC: nk1.ne.nq1 ',nk1,nq1(j)
      write(*,3)'SELF_SC: j,la,q,om=',j,la,q,omega
      stop
    endif

       cntr = cntr+1
       self = ocfunc(temp,omz,la1(j),nq1(j),la3(j),nq3(j))  !jk1,jk3)
       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif
 enddo nv3loop
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
 error = abs(self-oldself)
 write(ulog,5)'cnt,error,omz,self=',cnt2,error,omz,self
 if (aimag(self).lt.etaz) exit
enddo
 if (cnt2.ge.50) write(ulog,*)'SC CALCULATION OF THE SELF-ENERGY DID NOT CONVERGE ******************'

3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_sc
!============================================================
 subroutine three_ph_matrix(omega,wdth,res)
! frequencies, width are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! res=sum_q,la sum_12 0.5*|v3(q,la;1;2)|^2 delta(omega-w(q,la))  wk(q)
! it assumes the input momentum q is on the already-generated kmesh
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 real(8), intent(in) :: wdth,omega
 real(8), intent(out):: res
 integer j,nq,jq,la,jk1,jk2,jk3, nqq,nk1,inside
 real(8) omq,term,delta_g,sumz

 res = 0d0
 do nq=1,nibz  !nkc
 do la=1,ndyn
    call get_k_info(kibz(:,nq),nc,nqq,jk1,jk2,jk3,g1,g2,g3,inside)
!   omq = sqrt(abs(eigenval(la,nqq))) * cnst
    omq = eivalibz(la,nq)
    if (abs(omega-omq).gt.wdth*6) cycle
!   jq = nq+(la-1)*nibz
    jq = nqq+(la-1)*nkc
    term = delta_g(omega-omq,wdth) * wibz(nq) !*wdth*sqrt(2*pi)
    sumz=0
    v3loop: do j=1,nv3
       if (la2(j) .ne. la) cycle v3loop
       if (nq2(j) .ne. nqq) cycle v3loop
       jk1 = nq1(j)+(la1(j)-1)*nkc
       jk2 = nq2(j)+(la2(j)-1)*nkc
       jk3 = nq3(j)+(la3(j)-1)*nkc
       sumz = sumz + term * (v33(j)*conjg(v33(j)))
    enddo v3loop

    res = res+ sumz /nkc /2  *2*pi  ! should be comparable to iself in magnitude
  enddo
  enddo

3 format(i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))

 end subroutine three_ph_matrix
!===========================================================
! function ocfunc(temp,omz,j2,j3)  result(res)
 function ocfunc(temp,omz,l2,nk2,l3,nk3)  result(res)
! for given temperature and complex omz (in 1/cm) calculates the
! occupation-weighted joint DOS stored in res
! j2 and j3 are kpoint indices in the full FBZ
 use eigen
 use constants ! for cnst
 use kpoints   ! for nkc
 use params    ! for classical
 implicit none
 integer , intent(in) :: l2,l3,nk2,nk3
 integer j2,j3
 real(8), intent (in) :: temp
 real(8) om3,om2,nb3,nb2,nbe,resr,resi,et,om,delta_g,delta_l
 complex(8), intent (in) :: omz
 complex(8) res

!! nk2 = 1+mod(j2-1,nkc) ; l2 = 1+ (j2-nk2)/nkc
!! nk3 = 1+mod(j3-1,nkc) ; l3 = 1+ (j3-nk3)/nkc
! call nkla(j2,nkc,nk2,l2)
! call nkla(j3,nkc,nk3,l3)
 j2=(l2-1)*nkc+nk2
 j3=(l3-1)*nkc+nk3

 om3 = eigenval(l3,nk3)  ; nb3=nbe(om3,temp,classical)
 om2 = eigenval(l2,nk2)  ; nb2=nbe(om2,temp,classical)
! om3 = sqrt(abs(eigenval(l3,nk3))) * cnst ; nb3=nbe(om3,temp,classical)
! om2 = sqrt(abs(eigenval(l2,nk2))) * cnst ; nb2=nbe(om2,temp,classical)
! om=real(omz) ; et=-aimag(omz)

 if(j2.ne.j3) then
 !88888888888888888888888888888888888888888888888888888
 ! need a factor of 2 since v3(q,1,2)^2 * (f(1,2)+f(2,1))
 ! the factor of 2 is needed if np<nq (upper half) is excluded in the case where
 ! both sums are done in the full FBZ
 !88888888888888888888888888888888888888888888888888888
!    resi=-pi*((nb2+nb3+1)*(-delta_l(om3+om2-om,et)+ delta_l(om3+om2+om,et)) +  &
! &            (nb2-nb3  )*(-delta_l(om3-om2-om,et)+ delta_l(om3-om2+om,et)) )
!    resr=real(-(nb2+nb3+1)*(1d0/(om3+om2-omz)+ 1d0/(om3+om2+omz))    &
! &            -(nb2-nb3  )*(1d0/(om3-om2-omz)+ 1d0/(om3-om2+omz)) )
     res =     -(nb2+nb3+1)*(1d0/(om3+om2-omz)+ 1d0/(om3+om2+omz))    &
  &            -(nb2-nb3  )*(1d0/(om3-om2-omz)+ 1d0/(om3-om2+omz))
 else
!    resi=-(2*nb2+1) * pi* (-delta_l(2*om2-om,et)+ delta_l(2*om2+om,et))
!    resr=-(2*nb2+1) * real(1d0/(2*om2-omz)+ 1d0/(2*om2+omz))
     res =-(2*nb2+1) *     (1d0/(2*om2-omz)+ 1d0/(2*om2+omz))
 endif
! res=cmplx(resr,resi)

 end function ocfunc
!===========================================================
 function oc2(temp,omz,j2,j3)  result(res)
! for given temperature and complex omz (in 1/cm) calculates the
! occupation-weighted joint DOS stored in res
! j2 and j3 are kpoint indices in the full FBZ
 use eigen
 use constants ! for cnst
 use kpoints   ! for nkc
 use params    ! for classical
 implicit none
 integer , intent(in) :: j2,j3
 integer l2,l3,nk2,nk3
 real(8), intent (in) :: temp
 real(8) om3,om2,nb3,nb2,nbe
 complex(8), intent (in) :: omz
 complex(8) res

! nk2 = 1+mod(j2-1,nkc) ; l2 = 1+ (j2-nk2)/nkc
! nk3 = 1+mod(j3-1,nkc) ; l3 = 1+ (j3-nk3)/nkc
 call nkla(j2,nkc,nk2,l2)
 call nkla(j3,nkc,nk3,l3)

 om3 = eigenval(l3,nk3) ; nb3=nbe(om3,temp,classical)
 om2 = eigenval(l2,nk2) ; nb2=nbe(om2,temp,classical)
! om3 = sqrt(abs(eigenval(l3,nk3))) * cnst ; nb3=nbe(om3,temp,classical)
! om2 = sqrt(abs(eigenval(l2,nk2))) * cnst ; nb2=nbe(om2,temp,classical)

 if(j2.ne.j3) then   ! need a factor of 2 since v3(q,1,2)^2 * (f(1,2)+f(2,1))
    res = (1d0/(om3-om2-omz)+ 1d0/(om3-om2+omz))
 else
    res = (1d0/(2*om2-omz)+ 1d0/(2*om2+omz))
 endif

 end function oc2
!===========================================================
 subroutine nkla(j,nk,nkj,lj)
! for given one-dimensional index j, and number of kpoints nk, gives the index
! of the actual kpoint, nkj and branch lj assuming the outer loop is over bands
! and the inner loop is over kpoints
 implicit none
 integer j,nk,nkj,lj

 nkj = 1 + mod(j-1,nk)
 lj  = 1 + (j-nkj)/nk

 end subroutine nkla
!===========================================================
 subroutine get_frequencies(nkp,kp,dk,ndn,eival,nv,eivec,vg)
! also outputs the group velocities from HF theorem (exact) in units of c_light
 use params
 use io2
 use constants
 use born
 use geometry
 implicit none
 integer, intent(in) :: nkp,ndn,nv ! no of wanted eivecs
 real(8), intent(in) :: kp(3,nkp),dk(nkp)
 real(8), intent(out):: eival(ndn,nkp),vg(3,ndn,nkp)
 complex(8), intent(out) :: eivec(ndn,ndn,nkp)
 integer i,j,k,l,ier,nd2,al,ll
 integer, allocatable :: mp(:) !Map for ascending band sorting
 real(8), allocatable :: eivl(:)
 real(8) absvec,om
 complex(8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:)
 real(8) mysqrt
 character ext*3

 integer, allocatable :: mp1(:,:) !Map for projection band sorting
 real(8), allocatable :: eival_tmp(:) !Temp matrices for reordering eigenvalues
 complex(8), allocatable :: eivec_tmp(:,:)
 real(8), allocatable :: vg_tmp(:,:)

! open files and write the group velocities
! do l=1,ndn
!    write(ext,'(i3.3)')l
!    open(uvel+l,file='veloc-'//ext//'.dat')
!    write(uvel+l,*)'# la,i1,i2,i3,nk,kp(nk),eival(l,mapibz(nk)),vg(l,nk),length(vg(l,nk))'
! enddo
 nd2 = min(ndn,12)
!allocate(dynmat(ndn,ndn),mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3))
 allocate(mp(ndn),eivl(ndn),eivc(ndn,ndn))
! write(uio,*)'# dynmat allocated, printing: Kpoint, eival(j),j=1,ndyn)'

 kloop: do i=1,nkp

! write(uio,*)'############ before setting up dynmat, kp number ',i
    call finitedif_vel(kp(:,i),ndn,vg(:,:,i),eivl,eivc)

    do j=1,ndn
!      eival(j,i) = sqrt(abs(eivl(j)))*cnst 
       eival(j,i) = eivl(j)
       do l=1,nv
          eivec(j,l,i) = eivc(j,l)  !mp(l))
       enddo
    enddo


 enddo kloop

! Band structure projection sorting, added by Bolin
!   allocate(mp1(ndn,nkp))
!    write(*,*) "Projection Band Sorting for band structures!"
!   call band_sort_bs(nkp,ndn,kp,dk,eival,eivec,mp1)
!   allocate(eival_tmp(ndn),eivec_tmp(ndn,ndn),vg_tmp(3,ndn))
!   do i=1,nkp
!      eival_tmp=eival(:,i)
!      eivec_tmp=eivec(:,:,i)
!      vg_tmp=vg(:,:,i)
!      do j=1,ndn
!         eival(j,i)=eival_tmp(mp1(j,i))
!         eivec(:,j,i)=eivec_tmp(:,mp1(j,i))
!         vg(:,j,i)=vg_tmp(:,mp1(j,i))
!      end do
!      call write_eigenvalues(i,dk(i),kp(:,i),eival(:,i),ndn,uio)
!   end do
!   deallocate(mp1,eival_tmp,eivec_tmp,vg_tmp)

 deallocate(eivl,eivc,mp)  !,dynmat,ddyn)

 3 format(a,i5,9(1x,f19.9))
 4 format(99(1x,2(1x,g9.3),1x))
 5 format(a,i5,99(1x,2(1x,g9.3),1x))
 6 format(2i6,2x,99(1x,f9.3))
 7 format(i6,2x,3(1x,f7.3),1x,99(1x,f7.3))

 end subroutine get_frequencies
!=====================================================
  subroutine finitedif_vel(q0,ndn,vgr,evl0,evc0)
! calculates the group velocities in units of c_light from finite difference
! it would give zero near band crossings, thus HF is a better way to do it
  use constants
  implicit none
  integer ndn,i,j
  real(8) q0(3),vgr(3,ndn),q1(3),dq,om0,om1,aux(3,ndn),aux0(3,ndn)
  real(8) evl0(ndn),evlp(ndn),evlm(ndn)
  complex(8) evc0(ndn,ndn),evct(ndn,ndn)

  dq=3d-4  ! this is the ideal dq (at least for si)

   call get_freq(q0,ndn,aux0,evl0,evc0)
 
!  return

  do i=1,3
     q1=q0; q1(i)=q0(i)+dq
     call get_freq(q1,ndn,aux,evlp,evct)
     q1=q0; q1(i)=q0(i)-dq
     call get_freq(q1,ndn,aux,evlm,evct)
!    vgr(i,:)=(evlp-evlm)/2/dq / 2/sqrt(abs(evl0)) *cnst*1d-8 *2*pi !*c_light
     vgr(i,:)=(evlp-evlm)/2/dq   *1d-8 *2*pi !*c_light
  enddo
!  write(30,6)'q,dq,om=',q0,dq,evl0
!  do i=1,3
!     write(30,5)'i,fd_vg(i),hf_vg(i)=',i,vgr(i,:)*c_light,aux0(i,:)*c_light
!  enddo
!  write(30,*)'*******************************************'


 4 format(a,3(1x,f9.3),a)
 5 format(a,i5,99(1x,f9.3))
 6 format(a,3(1x,f9.3),1x,g12.6,9(1x,g10.4))

  end subroutine finitedif_vel
!============================================================
 subroutine gruneisen(nkp,kp,dk,ndn,eivl,eivc,ugr,grn)
! takes the eigenvalues (w^2) and eigenvectors calculated along some
! crystalline directions and calculates the corresponding mode
! gruneisen parameters
 use io2
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use eigen
 use constants
 implicit none
 integer ik,i0,nkp,la,al,be,ga,j,k,j0,k0,ta1,ta2,t,ugr,ndn    ! ,i3,j3,k3
 real(8) mi,mj,rx3(3),rx2(3),qq(3),denom,qdotr,omk,mysqrt
 real(8) kp(3,nkp),dk(nkp),eivl(ndn,nkp)
 complex(8) zz,one,term
 complex(8) grn(ndn,nkp), eivc(ndn,ndn,nkp)

 one = cmplx(1d0,0d0)
 write(ugr,*)'# la,nk,dk(nk),kp(:,nk),om(la,nk),gruneisen(la,nk))'
 do ik=1,nkp
    qq(:) = kp(:,ik)
 do la=1,ndn
!   write(ulog,6)'Starting gruneisen parameters for kpoint & band# ',ik,la,kp(:,ik)
!   write(ulog,*)' i,la,t,fcs_3(t),term(2),6*eival(la,i),rr3(3),grun(la,i)'
    grn(la,ik) = 0
    denom = 6 * eivl(la,ik)**2/cnst/cnst
    omk = eivl(la,ik)
    do i0=1,natoms0
         mi = atom0(i0)%mass
       tloop: do t=1,nterms(3)

         if ( i0 .ne. iatomterm_3(1,t) ) cycle tloop
         al = ixyzterm_3(1,t)
         be = ixyzterm_3(2,t)
         ga = ixyzterm_3(3,t)
!        i0 = iatomcell0(iatomterm_3(1,t))
         j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ; mj = atom0(j0)%mass
         k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k)
         ta1= al + 3*(i0-1)
         ta2= be + 3*(j0-1)

! rr2(:)=iatomcell(1,j)*r1+iatomcell(2,j)*r2+iatomcell(3,j)*r3
!  be careful: it has to be the translations R not atompos!
         rx2(:) = atompos(:,j) - atompos(:,j0)  ! R
         rx3(:) = atompos(:,k)                  ! R+tau
         qdotr =  ( qq .dot. rx2)
         zz = cdexp( ci * qdotr )
!! term = - ampterm_3(t)*fcs_3(igroup_3(t))*zz*eivc(ta1,la,i)*conjg(eivc(ta2,la,i))*rr3(ga)/sqrt(mi*mj)
         term = - fcs_3(t) * zz * eivc(ta2,la,ik)*conjg(eivc(ta1,la,ik))*rx3(ga)/sqrt(mi*mj)
         grn(la,ik) = grn(la,ik) + term
!        write(ulog,7)i,la,t,fcs_3(t),rr2,qdotr,zz,rr3,grn(la,ik)
       enddo tloop
    enddo
    grn(la,ik) = grn(la,ik)/denom
    if (aimag(grn(la,ik)) .gt. 1d-4) then
       write(ulog,*)' GRUNEISEN: ... has large imaginary part!! for nk,la=',ik,la,grn(la,ik)
!      stop
    endif
    write(ugr,6)' ',la,ik,dk(ik),kp(:,ik),omk,real(grn(la,ik))
 enddo
!   write(ugr,8)ik,dk(ik),kp(:,ik),(real(grn(la,ik)),la=1,ndn)
 enddo

5 format(4i7,9(1x,g10.4))
6 format(a,2i5,99(1x,g10.4))
7 format(i5,i5,i6,99(1x,g10.4))
8 format(i8,99(1x,g10.4))
! deallocate(eivl,eivc,kg,grn)
 
 end subroutine gruneisen
!============================================================
 subroutine gruneisen_fc
 use phi3
 use io2
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use eigen
 use constants
 implicit none
 integer i0,al,be,ga,j,k,k0,t,t2,igr2,t2old,grold
 real(8) rx3(3),phi  !mi,mj,

 if( .not. allocated(grun_fc)) allocate(grun_fc(nterms(2)))
 write(ulog,*)'GRUNEISEN_FC: term , phi2_ij , grun'
 do igr2=1,ngroups(2)
 t2old = 0 ; grold = 0 ; grun_fc=0
 do t2=1,nterms(2)
    if( igroup_2(t2) .eq. grold ) cycle
    if( igroup_2(t2) .ne. igr2) cycle
    i0 = iatomcell0(iatomterm_2(1,t2))
    j  = iatomterm_2(2,t2)
    al = ixyzterm_2(1,t2)
    be = ixyzterm_2(2,t2)
!    phi= ampterm_2(t2)*fcs_2(igroup_2(t2))
    phi= fcs_2(t2)
! k1    grun_fc(t2) = 0
!   grn = 0
    tloop: do t=1,nterms(3)
         if ( i0 .ne. iatomterm_3(1,t) ) cycle tloop
         if ( j  .ne. iatomterm_3(2,t) ) cycle tloop
         if ( al .ne. ixyzterm_3(1,t) ) cycle tloop
         if ( be .ne. ixyzterm_3(2,t) ) cycle tloop
         ga = ixyzterm_3(3,t)
         k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k)
         rx3(:) = atompos(:,k)
         grun_fc(t2) = grun_fc(t2) - rx3(ga)  / 6 /phi*fcs_3(t)
!        grn = grn - rr3(ga)  / 6 /phi*fcs_3(t) !ampterm_3(t) *fcs_3(igroup_3(t))
!         write(ulog,6)'phi3: ', t ,phi,-ampterm_3(t)*fcs_3(igroup_3(t)) , rr3(ga),grn
!         write(ulog,6)'phi3: ', t ,phi,-fcs_3(t) , rr3(ga),grn
    enddo tloop
    grold = igr2
    write(ulog,7)igr2,i0,al,j,be,phi,grun_fc(t2)
 enddo
 enddo

6 format(a,i7,9(1x,g10.4))
7 format('****',i7,2(4x,'(',i3,',',i1,')'),9(1x,g10.4))
 end subroutine gruneisen_fc
!============================================================
 subroutine mechanical(bulk,c11,c44,dlogv)
 use io2
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use constants
 implicit none
 integer i0,al,be,j,t
 real(8) bulk,rija,rijb,c11,c44,dlogv

! write(ulog,*)'BULK_MOD: i0, al, j, be, rija,rijb,fcs_2(t),bulk'
 bulk=0; c11=0 ; c44=0
 do i0=1,natoms0
 do t=1,nterms(2)
    if ( i0 .ne. iatomterm_2(1,t) ) cycle
    al =  ixyzterm_2(1,t)
    be =  ixyzterm_2(2,t)
    j  = iatomterm_2(2,t)
    rija = atompos(al,j)-atompos(al,i0)
    rijb = atompos(be,j)-atompos(be,i0)
!   write(ulog,3)i0, al, j, be, rija,rijb,fcs_2(t),bulk
    bulk = bulk - fcs_2(t)*(1-2*grun_fc(t)*dlogv)*rija*rijb
    if (al.eq.1 .and. be.eq.1) then
       c11 = c11 - fcs_2(t)*(1-2*grun_fc(t)*dlogv)*rija*rijb
    endif
 enddo
 enddo
 bulk = bulk/volume_r/18
 c11 = c11/volume_r/2

 write(ulog,*)'BULK_MODULUS: in eV/A^3 ',bulk
 bulk = bulk*1d30*ee
 write(ulog,*)'BULK_MODULUS: in SI, Mbar ',bulk,bulk*1d-11
 c11 = c11*1d30*ee
 write(ulog,*)'C11 in SI, Mbar ',c11,c11*1d-11
3 format(4(i4),9(2x,f10.5))

 end subroutine mechanical
!============================================================
 subroutine calculate_thermal(nk,wk,ndyn,eival,grn,tmn,tmx,ual,veloc)
 use io2
 use lattice
 use constants
 use params
 implicit none
 integer nk,ndyn,b,k,ual,itemp,nat,ntmp,iter,al,be
 real(8) wk(nk),eival(ndyn,nk),veloc(3,ndyn,nk),nbx,cv,cv_nk,cv2,nbe,dlogv,dalpha,kapoverl(3,3)
 real(8) temp,x,alpha,bulk_modulus,b0,tmn,tmx,gama,a1,c11,c44,etot,free,pres,pres0,magv
 complex(8) grn(ndyn,nk)

! bulk modulus = a_0^2/V d^2E/da^2 evaluated at a_0 (equilibrium lattice parameter)
! call mechanical(b0,c11,c44)  ! this is the bulk_mod at T=0

! B(T,Veq)=B(T=0,Veq(T=0)) - pres0

 nat = ndyn/3

 write(ual,'(a132)')'# temperature(K) ,alpha (1/K) , Cv (J/K/mol), gama , E_tot(J/mol) , &
 &   E_free , Pressure(GPa) , P0 , Bulk_mod(GPa) , kappa/L(nm) '

 ntmp=ntemp !60
 do itemp=1,ntmp

    temp=tmn+(tmx-tmn)*(itemp-1)**3/(ntmp-1d0+1d-8)**3   ! this temp is in Kelvin
    if (temp.le.0) then
       write(ulog,*)'temperature not in the proper range!!',temp
       stop
    endif
    alpha=1d2; dalpha=1d9; iter=0
    do while (abs(dalpha) .gt. abs(alpha)/1000 .and. iter .lt. 50)
       iter=iter+1
       dlogv=temp*alpha
       call mechanical(b0,c11,c44,dlogv)  ! this is the bulk_mod at T=0
       call energies(nk,wk,ndyn,eival,grn,temp,etot,free,pres,pres0,cv2)
       bulk_modulus = b0 - pres0

       a1=0 ; cv=0 ; gama = 0; kapoverl=0
! in the sum over k, all quantities are per unitcell
       do k=1,nk
       do b=1,ndyn
!         if(eival(b,k) .lt.0) then
!            x=0
!            write(ulog,*) 'ALPHA: negative eival for band,kp# ',b,k,eival(b,k)
!            write(ulog,*) ' will use its absolute value instead!'
!      else
!         endif
! to convert eV/A^2/mass to Hz we need to take the sqrt, multiply by cnst*100*c_light
          x=(h_plank*eival(b,k)*100*c_light)/k_b/temp
          if (x.gt.60) then
              cv_nk = 0
          else
              nbx=nbe(x,1d0,classical)
              cv_nk = x*x*nbx*(1+nbx) !/4/sinh(x/2)/sinh(x/2)
          endif
          cv = cv + cv_nk*wk(k)   ! this is in J/K units : per unit cell
          a1 = a1 + grn(b,k)*cv_nk*wk(k)
          magv=sqrt(dot_product(veloc(:,b,k),veloc(:,b,k)))
          do al=1,3
          do be=1,3
             kapoverl(al,be)=kapoverl(al,be)+cv_nk*wk(k)*veloc(al,b,k)*veloc(be,b,k)/(magv+1d-20)
          enddo
          enddo
       enddo
       enddo
! multiplied by MFP(nm), it is the thermal conductivity
       kapoverl = kapoverl *k_b/volume_r*1d30*c_light*1d-9
       gama = a1 / cv
       cv = cv/nat*n_avog*k_b  ! this is per mole    )/(volume_r*1d-30)  !
       if (abs(cv2-cv).gt.1d-5 ) then
          write(ulog,4)'CALCULATE_THERMAL:temp, cv from energies ne cv ',iter,temp,cv2,cv
!      stop
       endif
       dalpha = a1/bulk_modulus /n_avog*nat/(volume_r*1d-30) - alpha
       alpha = a1/bulk_modulus /n_avog*nat/(volume_r*1d-30)
       write(*,4)'CALCULATE_THERMAL:i,T,a,da=',iter,temp,alpha,dalpha,magv,kapoverl
    enddo
    write(ual,3)temp,alpha,cv,gama,etot,free,pres*1d-9,pres0*1d-9,bulk_modulus*1d-9,kapoverl
 enddo

3 format(99(2x,g10.4))
4 format(a,i3,9(2x,g11.5))

 end subroutine calculate_thermal
!============================================================
 subroutine energies(nk,wk,ndyn,eival,grn,temp,etot,free,pres,pres0,cv)
! calculate total and free energies within QHA, at a given temperature (temp in Kelvin)
 use io2
 use params
 use lattice
 use constants
 implicit none
 integer nk,ndyn,b,k,nat
 real(8) wk(nk),eival(ndyn,nk)
 real(8) temp,x,cv_nk,cv,hw,free,etot,pres,nbe,mdedv,pres0,nbx
 complex(8) grn(ndyn,nk)

    nat = ndyn/3
    if (temp.le.0) then
       write(ulog,*)'temperature not in the proper range!!',temp
       stop
    endif
    etot=0 ; cv=0 ; free=0 ; pres=0
    mdedv= 0.35388/(20.8**3-20.0**3)/ab**3 ! this is -dE/dV in eV/A^3
    mdedv = mdedv*1d+30*ee ! (in J/m^3 )
    do k=1,nk
    do b=1,ndyn
       if(eival(b,k) .lt.0) then
          x=0
          write(ulog,*) 'ALPHA: negative eival for band,kp# ',b,k,eival(b,k)
          write(ulog,*) ' will use its absolute value instead!'
!      else
       endif
! to convert eV/A^2/mass to Hz we need to take the sqrt, multiply by cnst*100*c_light
       x=(h_plank*eival(b,k)*100*c_light)/k_b/temp
       if (x.gt.60) then
           cv_nk = 0
       else
           nbx=nbe(x,1d0,classical)
           cv_nk = x*x*nbx*(1+nbx) !/4/sinh(x/2)/sinh(x/2)
       endif
       hw = x*k_b*temp  ! hbar*omega in Joules
       cv  = cv + cv_nk*wk(k)   ! this is in J/K units : per unit cell
       etot= etot + hw * (nbe(x,1d0,classical) + 0.5d0) * wk(k)
       free= free + (0.5d0*hw + k_b*temp*log(1-exp(-x))) * wk(k)
       pres= pres + hw * (nbe(x,1d0,classical) + 0.5d0) * wk(k) * grn(b,k)
    enddo
    enddo
    cv = cv/nat*n_avog*k_b  ! this is per mole    )/(volume_r*1d-30)  !
    etot = etot/nat*n_avog   ! convert from Joule/cell to Joule per mole
    free = free/nat*n_avog   ! convert from Joule/cell to Joule per mole
    pres0= pres/(volume_r*1d-30) ! * 1d-8  ! 1d-8 is to convert to kbar
    pres = pres0+ mdedv
3 format(9(2x,g11.5))

 end subroutine energies
!============================================================
  subroutine get_kindex(q,nkc,kpc,ik)
  use  geometry
  implicit none
  integer ik,nkc,i,j
  real(8) kpc(3,nkc),q(3)

  ik=0
  mainlp: do i=1,nkc
     do j=1,3
        if (.not. (abs(q(j)-kpc(j,i)) .myeq. 0d0) ) exit
        if (j.eq.3) then
           ik=i
           exit mainlp
        endif
     enddo
  enddo mainlp
  if (ik.eq.0) then
     write(*,*)'GET_INDEX: could not find the index of',q
     stop
  endif

  end subroutine get_kindex
!============================================================
 subroutine get_k_info(q,N,nk,i,j,k,g1,g2,g3,inside)
! for a vector q(3) in the primitive cell of the reciprocal space, defined on
! a mesh N(1),N(2),N(3), this subroutine finds the three indices of q
! and its number nk based on the triple loop
! nk=0 do i=0,N(1)-1 ; do j=0,N(2)-1; do k=0,N(3)-1; nk=nk+1 ;q=i*G1/N1+j*G2/N2+k*G3/N3

 use geometry
 implicit none
 integer, intent(in):: N(3)
 real(8), intent(in):: q(3)
 type(vector), intent(in):: g1,g2,g3
 integer, intent (out) :: nk,i,j,k,inside
 integer indexg

! this was for when the kpoints were between 0 and N-1
    call get_components_g(q,N,i,j,k,g1,g2,g3,inside)
    nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
 !  nk= indexg(i,j,k,-N(1)/2,N(1)/2-1,-N(2)/2,N(2)/2-1,-N(3)/2,N(3)/2-1)
    return

! this is for when the k vectors are between -N/2 and N/2-1
  if (mod(N(1),2).eq.0 .and. mod(N(2),2).eq.0 .and. mod(N(3),2).eq.0 ) then
! use this if shifted by -0.5(g1+g2+g3)
    call get_components_g_centered (q,N,i,j,k,g1,g2,g3,inside)
!   write(*,*)'for even ni, ijk=',i,j,k
  else   ! do the regulat thing
    call get_components_g(q,N,i,j,k,g1,g2,g3,inside)
  endif
  nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
! write(*,*)'the corresponding nk=',nk

 end subroutine get_k_info
!============================================================
 subroutine get_k_info_cent(q,N,nk,i,j,k,g1,g2,g3,inside)
! for a vector q(3) in the primitive cell of the CENTERED-reciprocal space, defined on
! a mesh N(1),N(2),N(3), this subroutine finds the three indices of q
! and its number nk based on the triple loop
! nk=0 do i=0,N(1)-1 ; do j=0,N(2)-1; do k=0,N(3)-1; nk=nk+1 ;q=i*G1/N1+j*G2/N2+k*G3/N3

 use geometry
 implicit none
 integer, intent(in):: N(3)
 real(8), intent(in):: q(3)
 type(vector), intent(in):: g1,g2,g3
 integer, intent (out) :: nk,i,j,k,inside
 integer indexg

! this was for when the kpoints were between 0 and N-1
    call get_components_g_centered (q,N,i,j,k,g1,g2,g3,inside)
    nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
 !  nk= indexg(i,j,k,-N(1)/2,N(1)/2-1,-N(2)/2,N(2)/2-1,-N(3)/2,N(3)/2-1)
    return

! this is for when the k vectors are between -N/2 and N/2-1
  if (mod(N(1),2).eq.0 .and. mod(N(2),2).eq.0 .and. mod(N(3),2).eq.0 ) then
! use this if shifted by -0.5(g1+g2+g3)
    call get_components_g_centered (q,N,i,j,k,g1,g2,g3,inside)
!   write(*,*)'for even ni, ijk=',i,j,k
  else   ! do the regulat thing
    call get_components_g(q,N,i,j,k,g1,g2,g3,inside)
  endif
  nk= indexg(i,j,k,1,N(1),1,N(2),1,N(3))
! write(*,*)'the corresponding nk=',nk

 end subroutine get_k_info_cent
!===========================================================
 subroutine get_k_info3(q,nkt,ex) !,i1,j1,k1,gg1,gg2,gg3,inside)
! scans the kpoints to identify the index of the input q
 use geometry
 use kpoints
 use params
 use io2
 implicit none
 real(8) q(3),ll
 integer inside,i1,j1,k1,nkt,mc(3),i
 logical ex

3 format(a,i6,2x,3(1x,g9.3),2x,f10.4,L3)

 nkt=0; ex=.false.
 loop: do i=1,nkc
!    if( (q(1).myeq.kpc(1,i)) .and. (q(2).myeq.kpc(2,i)) .and. (q(3).myeq.kpc(3,i)) ) then
    ll=length(q-kpc(:,i))
!   if( length(q-kpc(:,i)) .lt. 1d-4) then
!       write(*,3)'i,k,q,k-q=',i,kpc(:,i),q,ll
    if( ll .lt. 1d-4) then
       nkt=i
       ex=.true.
       exit loop
    endif
 enddo loop
! write(*,*)'nkt=',nkt
 if (nkt.eq.0) then
   write(ulog,3)'GET_K_INFO3: qpoint not found! nq,q,l,exist?',nkt,q,length(q),ex
!  stop
 endif
 end subroutine get_k_info3
!-------------------------------------
 subroutine collision_matrix(ndyn,nkp,kp,eival,temp,w0,c_matrix)
! this subroutine uses the Fermi-Golden-Rule to get the rate equations
! and construct the linearized collision operator for 3-phonon processes
! to which one can eventually add other scattering meachnisms, and
! diagonalize the matrix to get the rates (inverse of relaxation times)
! as eigenvalues. The latter will be used in Boltzmann equation to
! get the non-equilibruim part of the distribution function
! input is nq= # of q points in the FBZ times # of phonon branches
! temp is in 1/cm
! output is the collision matrix (C_matrix) determined from the rates by using FGR
 use constants
 use params
 implicit none
 integer nq,np,nk1,nk2,nk3,l1,l2,l3,ndyn,nkp,nk,inside
 real(8) c_matrix(nkp*ndyn,nkp*ndyn),ratepk,ratepq,rateqk
 real(8) w0,eival(ndyn,nkp),kp(3,nkp),nbe,temp,nbq,x1,nbp,x2,nbk,x3

 open(321,file='collision_mat.dat')
! write(321,*)'# n , nk , la , c_matrix(nq,nq) , rate(cm^-1) '
 write(321,*)'#  nk , DIAG[c_matrix(nk,la)] la=1,ndyn , in(cm^-1) '

 c_matrix=0

 nq = 0
 do nk1=1,nkp
 do l1=1,ndyn
    nq=nq+1    ! nq is the generic line index
    x1=eival(l1,nk1)
    nbq = nbe(x1,temp,classical)

    np=0
    do nk2=1,nkp
    do l2=1,ndyn
       np=np+1  ! np is the generic column index
       x2=eival(l2,nk2) !sqrt(abs(eival(l2,nk2)))*cnst
       nbp = nbe(x2,temp,classical)

       nk=0
       do nk3=1,nkp
       do l3=1,ndyn
          nk=nk+1  ! nk is the generic index over which the sums are made
          x3=eival(l3,nk3)  !sqrt(abs(eival(l3,nk3)))*cnst
          nbk = nbe(x3,temp,classical)

!         call FGR(la1,la2,la3,nk2,nk3,ndyn,nkp,kp,eival,temp,w0,rate1,inside)
!         call FGR(la2,la3,la1,nk3,nk1,ndyn,nkp,kp,eival,temp,w0,rate2,inside)
!         call FGR(la3,la1,la2,nk1,nk2,ndyn,nkp,kp,eival,temp,w0,rate3,inside)
          call FGR(l1,l2,l3,nk2,nk3,ndyn,nkp,kp,eival,temp,w0,ratepk,inside)
          call FGR(l2,l3,l1,nk3,nk1,ndyn,nkp,kp,eival,temp,w0,rateqk,inside)
          call FGR(l3,l1,l2,nk1,nk2,ndyn,nkp,kp,eival,temp,w0,ratepq,inside)

          if (np.eq.nq) then
!            c_matrix(nq,np)= c_matrix(nq,np) +0.5*rate1 +rate2
!            c_matrix(nq,np)= c_matrix(nq,np) +0.5*(rate1+rate2+rate3) ! should produce the same as above line
             c_matrix(nq,nq)= c_matrix(nq,nq) +0.5*ratepk* (1+nbp+nbk) + ratepq*(nbp-nbk)
          else
!            c_matrix(nq,np)= c_matrix(nq,np) - rate1 - rate2 + rate3
             c_matrix(nq,np)= c_matrix(nq,np) + (ratepq+ratepk)*(nbq-nbk) - rateqk*(1+nbq+nbk)
          endif

       enddo
       enddo

    enddo
    enddo

!   rateq=c_matrix(nq,nq)/nbq/(1+nbq)
!   write(321,3) nq,nk1,l1,c_matrix(nq,nq),rateq

 enddo
    nk=(nk1-1)*ndyn
    write(321,3) nk1,(c_matrix(nk+l1,nk+l1),l1=1,ndyn)
 enddo

 close(321)
3 format(1(1x,i9),99(1x,g11.5))

 end subroutine collision_matrix
!===========================================================
 subroutine FGR(l1,l2,l3,nk2,nk3,ndyn,nkp,kp,eival,temp,w0,rate,inside)
! this subroutines uses the FGR to calculate the rates needed in the collision matrix
! input is nq= # of q points in the FBZ times # of phonon branches
! output is rate(1,2,f)= # collision rate to go from state f to states (1,2)
! the delta function width is determined by the average frequency separation as
! delta(x,w0) = w0/pi/(x^2+w0^2) with w0=<w(q+1)-w(q)> =<vgroup> delta_q
! input (besides w0) is s1,s2,sf, integers representing the initial 2 phonon and the
! final one-phonon state : si=1,nk*ndyn
! output is the rate =  2pi/hbar |<i|V_3|j,k>|^2 delta(E_i-E_j-E_k) (all converted to cm^-1)
! to convert Joules to 1/cm need to divide by h_plank*100*c_light
 use constants
 use phi3
 use lattice
 use params
 implicit none
 integer l1,l2,l3,nk1,nk2,nk3,nkp,ndyn,i,j,k,inside
 real(8) w0,eival(ndyn,nkp),kp(3,nkp),nbe,temp,nb1,q(3),x1,rate,delta_l,x2,x3 !,nb2,nb3
 complex(8) w3

! print*,'FGR: used broadening w0 for this calculation is=',w0,' cm^-1'

 q(:)=-kp(:,nk2)-kp(:,nk3)
 call get_k_info(q,NC,nk1,i,j,k,g1,g2,g3,inside)

 x1=sqrt(abs(eival(l1,nk1)))*cnst  ! convert frequencies to 1/cm as temp is also in 1/cm
 nb1 = nbe(x1,temp,classical)
 x2=sqrt(abs(eival(l2,nk2)))*cnst
! nb2 = nbe(x2,temp,classical)
 x3=sqrt(abs(eival(l3,nk3)))*cnst
! nb3 = nbe(x3,temp,classical)

 x1=(x1-x2-x3)
 w3 = v3(nk2,nk3,l1,l2,l3)
! rate = 2*pi*delta_l(x1,w0) * w3*conjg(w3) *2*pi * (nb1+1)*nb2*nb3
 rate = 2*pi*delta_l(x1,w0) * w3*conjg(w3) *2*pi * (nb1+1)*nb1
 if ( rate.lt.0 ) then
    write(*,*)'FGR: rate<0 ',rate
    stop
 endif
! extra 2pi is needed to convert final results to cm^-1

 end subroutine FGR
!===========================================================
 subroutine write_collision_times(ndyn,nkp,eival,temp,inv_relax_t)
 use constants
 use params
 implicit none
 integer la,nk,nkp,ndyn,nq
 real(8) eival(ndyn,nkp),nbe,temp,inv_relax_t(ndyn*nkp),nbq,x1

 open(322,file='colmat_eivals.dat')
 write(322,*)'# n ,nk,la,eival_coll_matrx(nq) in cm^-1, inverse RT (cm^-1) and Thz'
 nq = 0
 do nk=1,nkp
 do la=1,ndyn
    nq=nq+1    ! nq is the generic line index

    x1=sqrt(abs(eival(la,nk)))*cnst
    nbq = nbe(x1,temp,classical)

    write(322,3) nq,nk,la,inv_relax_t(nq),inv_relax_t(nq)*100*c_light*1d-12
 enddo
 enddo
 close(322)
3 format(3(1x,i9),9(2x,g11.5))

 end subroutine write_collision_times
!===========================================================
 subroutine get_negatives(nk,kp,map,nps)
! this subroutine finds and eliminates the negatives in the kp mesh
! on output: if map(i), (i=1,nps) is the index of the needed kpoints
! the negative k's are kp(:,map(i)), i=nps+1,nk
 use kpoints
 use lattice
 use io2
 use params
 implicit none
 integer nk,i,map(nk),nps,i1,j1,k1,inside,mk,i2,j2,k2,k,match,cnt
 real(8) kp(3,nk) !,q(3)
 logical exists

 map=0; nps=0
 do i=1,nk
    call get_k_info_cent( kp(:,i)-shft,NC, k,i1,j1,k1,g1,g2,g3,inside)
    if (k.ne.i) then
       write(ulog,*)'GET_NEGS: error k.ne.i ',k,i
       write(ulog,*)'GET_NEGS: i,j,k =      ',i1,j1,k1
       stop
    endif

    call get_k_info3(-kp(:,i),mk,exists)
! if -k is not in the negs list, then add it to the map, else disregard it
 if (exists) then  ! if it doesn't exist, we won't substitute the eigenvector
    if (mk.lt.i) then ! it is already included
       cycle
    else
       nps=nps+1
       map(nps)=i
    endif
    if(verbose) write(ulog,3)'GET_NEGS: i(k),i(-k),npos, k,-k=',i,mk,nps,i1,j1,k1
 endif
 enddo
 write(ulog,*)'GET_NEGS: # of included kpoints=',nps

! do i=1,nk
!    write(ulog,*)'ik,mapos(ik)=',i,map(i)
! enddo

return

! now end the array map with the negative ones
 cnt=nps
 do i=1,nk
    match=0
    iloop: do i1=1,nps
       if(map(i1).eq.i) then
          match=1
          exit iloop
       endif
    enddo iloop
    if(match.ne.1) then
       cnt=cnt+1
       map(cnt)=i
       if(verbose) write(ulog,3)'negatives: cnt,map=',cnt,i
    endif
 enddo
3 format(a,3i5,2(3x,3i3))
 end subroutine get_negatives
!===========================================================
 subroutine subst_eivecs(ndn,nk,eival,eivec,kp,nposi,map)
! this subroutine keeps the eivecs(q) for q in the nposi list and replaces the
! other eivec(q') (for q'=-q) by the conjugate of eivec(q) in order to get
! rid of the ambiguity in the case of the degenerate case, and assure e(-q)=conjg(e(q)).
 use lattice  ! needed to access NC,g1,g2,g3
 use io2
 use kpoints
 implicit none
 integer nk,ndn,nposi,map(nk),i1,j1,k1,n,inside,mk3,n3,la,mu
 real(8) kp(3,nk) ,eival(ndn)
 complex(8) eivec(ndn,ndn,nk),aux(ndn)
 logical exists

 write(ulog,*)'SUBST_EIVECS: nk,nposi=',nk,nposi,'=============='
! do n=1,nk
!    write(ulog,*)'ik,map(ik)=',n,map(n)
! enddo
! do n=nposi+1,nk
 do n=1,nposi  ! these are the kpoints that have a negative in the kpc list
    n3=map(n)      ! actual index of the kpoint
!   call get_k_info(-kp(:,n3)-shft,NC,mk3,i1,j1,k1,g1,g2,g3,inside)
!write(ulog,*)'npos,ik,k ',n,n3,kp(:,n3)
    call get_k_info3(-kp(:,n3),mk3,exists)
  if (exists) then
    do la=1,ndn
       aux=eivec(:,la,n3)
!      write(ulog,*)'for kpoint,branch=',n3,la
       do mu=1,ndn
          if(abs(eivec(mu,la,mk3) - conjg(aux(mu))) .gt. 1d-4) then
          if(abs(eivec(mu,la,mk3) + conjg(aux(mu))) .gt. 1d-4) then
            write(ulog,4)'WILL SWITCH:mu,k,eivl,ev(-k),ev*(k) =',mu,kp(:,n3),eival(mu),eivec(mu,la,mk3) , conjg(aux(mu))
          endif
          endif
       enddo
       eivec(:,la,mk3) = conjg(aux)
    enddo
  endif
 enddo
 write(ulog,*)'EXITING SUBST_EIVECS =========================='

4 format(a,i4,99(1x,f9.3))

 end subroutine subst_eivecs
!===========================================================
 subroutine check_mdyn(ndn,nk,kp,eivec,eival)
! form the dynamical matrix from existing eigenvectors and eigenvalues
! and diagonalize and verify consistency :eivec(-q)=-conjg(eivec(q))
 use lattice  ! needed to access NC,g1,g2,g3
 implicit none
 integer ndn,nk,j,k,l,i,ier,i1,j1,k1,mi,inside
 real(8) kp(3,nk),eival(ndn,nk),eivl(ndn)
 complex(8) eivec(ndn,ndn,nk),dynm(ndn,ndn),eivc(ndn,ndn),d2(ndn,ndn)

 do i=1,nk
    dynm=0
    do j=1,ndn
    do k=1,ndn
      do l=1,ndn
        dynm(j,k)=dynm(j,k)+eival(l,i)*eival(l,i)*cnst*cnst*  &
&                           eivec(j,l,i)*conjg(eivec(k,l,i))
      enddo
    enddo
    enddo
    call get_k_info(-kp(:,i),NC,mi,i1,j1,k1,g1,g2,g3,inside)
    d2=0
    do j=1,ndn
    do k=1,ndn
      do l=1,ndn
        d2(j,k)=d2(j,k)+eival(l,mi)*eival(l,mi)*cnst*cnst*     &
&                            eivec(j,l,mi)*conjg(eivec(k,l,mi))
      enddo
    enddo
    enddo

    do j=1,ndn
    do k=1,ndn
       if(abs( d2(j,k)-conjg(dynm(j,k))) .gt.1d-5) then
         write(*,4)'CHECK_MDYN: j,k,d(-q),d(q)=',j,k,d2(j,k),dynm(j,k)
       endif
    enddo
    enddo
!   return

    call diagonalize(ndn,dynm,eivl,ndn,eivc,ier)
    do j=1,ndn
       if( abs(cnst*sqrt(abs(eivl(j)))-eival(j,i)).gt.1d-4) then
          write(*,3)'CHECK_MDYN:j,eivl,eival=',j,cnst*sqrt(abs(eivl(j))),eival(j,i)
       endif
       do k=1,ndn
         if(abs(eivc(k,j)-eivec(k,j,i)).gt.1d-4) then
         if(abs(eivc(k,j)+eivec(k,j,i)).gt.1d-4) then
            write(*,4)'CHECK_MDYN:j,k,eiv(j),eivecs(k,j)=',j,k,eival(j,i),eivc(k,j),eivec(k,j,i)
         endif
         endif
       enddo
    enddo
 enddo
3 format(a,i6,9(1x,f15.6))
4 format(a,2i6,f15.7,9(1x,f11.4))

 end subroutine check_mdyn
!===========================================================
 subroutine write_eivecs(ndn,nk,kp,eigenvec)
 use lattice
 use io2
 use geometry
 use constants
 use kpoints ! for shft
 implicit none
 integer nk,ndn,i,la,j,i2,j2,k2,inside,l,i3,j3,k3
 real(8) kp(3,nk),w(3)
 complex(8) eigenvec(ndn,ndn,nk)

 write(ulog,*)'WRITE_EIGVECS: eivc(q)-conjg(eivc(-q)),eivc(q),eivc(-q)'
do i=1,nk
! bring kp(j) into the FBZ

   call get_k_info_cent(kp(:,i)-shft,NC,l,i3,j3,k3,g1,g2,g3,inside)
   if (i.ne.l) then
      write(ulog,*) 'i .ne. l ',i,l
      stop
   endif
!  call bring_to_cell_c(-kp(:,i),g1,g2,g3,r1,r2,r3,w)
!  w =w /2/pi
   w = -kp(:,i)
   call get_k_info_cent(w-shft,NC,j,i2,j2,k2,g1,g2,g3,inside)
!  write(ulog,5)'k,-k=',i3,j3,k3,' ** ',i2,j2,k2,kp(:,i),w

  if (j.lt.nk) then
   do la=1,ndn
   do l=1,ndn
      if(abs(eigenvec(l,la,i)-conjg(eigenvec(l,la,j))).gt.1d-4) then
      if(abs(eigenvec(l,la,i)+conjg(eigenvec(l,la,j))).gt.1d-4) then
          write(ulog,8)'la,l,k,-k=', la,l,i,j,eigenvec(l,la,i)-conjg(eigenvec(l,la,j)) &
      &     ,eigenvec(l,la,i),eigenvec(l,la,j)
      endif
      endif
   enddo
   enddo
  else
    write(ulog,4)'kpoint i in FBZ=',i,kp(:,i)
    write(ulog,4)'was taken to -k=',j,w
  endif
enddo
4 format(a,i5,2x,99(1x,f8.3))
5 format(a,3i5,a,3i5,2x,99(1x,f8.3))
8 format(a,4i6,99(1x,f8.3))

 end subroutine write_eivecs
!===========================================================
 subroutine mode_thermal_conductivity(nk,wk,tau_inv,veloc,omg,temp,kappa_q)
! temp is in cm^-1, veloc in c_light, omg in 1/cm, tau_inv in cm_1
! input are the freqs, velocs , relaxation times and weights for each band
! output is the thermal conductivity for that band
 use constants
 use params
 use lattice
 implicit none
 integer, intent(in) ::  nk
 integer i,al,be
 real(8),intent(in):: tau_inv(nk),veloc(3,nk),omg(nk),wk(nk),temp
 real(8), intent(out) :: kappa_q(3,3)
 real(8) cv,x,tau,nbe,nbx

  kappa_q=0 ; cv=0
  do i=1,nk
     x=omg(i)/temp   ! temp is in 1/cm
     if (x.gt.40) then
       cv = 0
     else
       nbx=nbe(x,1d0,classical)
       cv=nbx*(nbx+1)*x*x  !  x*x/4/sinh(x/2)/sinh(x/2)
     endif
     tau=1/tau_inv(i)
     if(tau.lt.1d7) then  ! exclude gamma=0
!      kappa_q=kappa_q+cv*wk(i)*tau*sum(veloc(:,i)*veloc(:,i))/3
        do al=1,3
        do be=1,3
       kappa_q(al,be)=kappa_q(al,be)+cv*wk(i)*tau*veloc(al,i)*veloc(be,i)
        enddo
        enddo
     endif
  enddo
  kappa_q = kappa_q *k_b/100*c_light/volume_r*1d30  ! converted to SI (W/mK)

 end subroutine mode_thermal_conductivity
!===========================================================
 subroutine mode_thermal_conductivity_single_k(wj,tau_invj,velocj,eivj,temp,kappa_qj)
! this is the contribution of a single mode j=(k,la) with given freq omj,velocj,tau_invj
! temp is in cm^-1, veloc in c_light, omg in ev/ang^2/uma, tau_inv in cm_1
! input are the freqs, velocs , relaxation times and weights for each band
! output is the thermal conductivity for that band
 use constants
 use params
 use lattice
 implicit none
 integer i,al,be
 real(8),intent(in):: wj,tau_invj,velocj(3),eivj,temp
 real(8), intent(out) :: kappa_qj(3,3)
 real(8) cv,x,tau,nbe,nbx

     kappa_qj=0d0 ; cv=0d0
!    x=sqrt(abs(eivj))*cnst/temp   ! temp is in 1/cm
     x=eivj/temp   ! temp is in 1/cm
     if (x.gt.40) then
       cv = 0d0
     else
       nbx=nbe(x,1d0,classical)
       cv=nbx*(nbx+1)*x*x  !  x*x/4/sinh(x/2)/sinh(x/2)
     endif
     tau=1/tau_invj
     if(abs(tau).lt.1d7) then  ! exclude gamma=0
!      kappa_q=kappa_q+cv*wk(i)*tau*sum(veloc(:,i)*veloc(:,i))/3
        do al=1,3
        do be=1,3
           kappa_qj(al,be)=kappa_qj(al,be)+cv*wj*tau*velocj(al)*velocj(be)
        enddo
        enddo
     endif
     kappa_qj = kappa_qj *k_b/100*c_light/volume_r*1d30  ! converted to SI (W/mK)

 end subroutine mode_thermal_conductivity_single_k
!===========================================================
 subroutine thermal_conductivity(nk,wk,ndn,tau_inv,veloc,omg,temp,kappa_q)
! temp is in cm^-1, veloc in c_light, omg in ev/ang^2/uma, tau_inv in cm_1
 use constants
 use params
 use lattice
 implicit none
 integer, intent(in) ::  nk,ndn
 integer i,la,al,be
 real(8),intent(in):: tau_inv(ndn,nk),veloc(3,ndn,nk),omg(ndn,nk),wk(nk),temp
 real(8), intent(out) :: kappa_q(3,3)
 real(8) cv,x,tau,nbe,nbx

  kappa_q=0 ; cv=0
  do i=1,nk
  do la=1,ndn
     x=omg(la,i)/temp   ! temp is in 1/cm
     if (x.gt.60) then
       cv = 0
     else
       nbx=nbe(x,1d0,classical)
       cv=nbx*(nbx+1)*x*x  ! x*x/4/sinh(x/2)/sinh(x/2)
     endif
     tau=1/tau_inv(la,i)
     if(tau.lt.1d9) then  ! exclude gamma=0
!*********
! should use tau_inv=2*gamma
!*********
        do al=1,3
        do be=1,3
       kappa_q(al,be)=kappa_q(al,be)+cv*wk(i)*tau*veloc(al,la,i)*veloc(be,la,i)
        enddo
        enddo
!      kappa_c=kappa_c+   wk(i)*tau*sum(veloc(:,la,i)*veloc(:,la,i))
     endif
  enddo
  enddo
  kappa_q = kappa_q *k_b/100*c_light/volume_r*1d30  ! converted to SI (W/mK)
! kappa_c = kappa_c/3 *k_b/100*c_light/volume_r*1d30  ! converted to SI (W/mK)

 end subroutine thermal_conductivity
!===========================================================
 subroutine tau_klemens(temp,eival,grn,veloc,wmax,mass,tau_inv)
 use constants
 implicit none
 real(8) temp,eival,grn,veloc(3),wmax,tau_inv,mass,v2,kbt,w2

 v2 = sum(veloc*veloc)*c_light*c_light
 kbt = temp*100*c_light*h_plank  ! kT in Joules
 w2 = eival*cnst*cnst  ! *** TO BE CHECKED **** ! in cm^-2
 tau_inv = kbt/(mass*uma*v2/2) /wmax*w2*grn*grn  ! this is in cm-1
! if (tau_inv .lt. 1d-5) tau_inv=1d30

 end subroutine tau_klemens
!===========================================================
 function eta(om,tempk)
 implicit none
 real(8) om,eta,tempk
 eta = 1d-4*om*om*0.01*(tempk/20)*(1+ 1d-6*om*om*om*0.2)
 end function eta
!===========================================================
 function cross_section(q,la,omega,eival,self) result(sigma)
! calculates the neutron cross section within BORN approx based on the eivals
! of the dynamical matrix and the 3-phonon self energy
 use constants
 implicit none
 real(8) q(3),omega,sigma,xnum,denom,omq,eival,mysqrt
 complex(8) self
 integer la

! omq = cnst*mysqrt(eival)
 !omq = cnst*sqrt(abs(eival))
 omq = eival
 xnum = 2 *omq*aimag(self)
 denom = xnum*xnum + (omega*omega-omq*omq-2*omq*real(self))**2
 sigma = xnum/denom

 end function cross_section
!===========================================================
 subroutine nonanal(q,dynmat,ndim,ddyn)
 use constants
 use lattice
 use atoms_force_constants
 use born

 integer, intent(in) ::ndim
 complex(8), intent(inout) :: dynmat(ndim,ndim),ddyn(ndim,ndim,3)
 real(8), intent(in) :: q(3)
 real(8) zag(3),zbg(3),qeq,om0,ma,mb,rr(3),dqeq(3)
 integer na,nb,i,j,al
 real(8) eps_scale,q2,gg,term

 eps_scale=8.8541878176D-12/1D10/ee
 gg=(length(g1)*length(g2)*length(g3))**0.33
! if ( gg .lt. 4*rho) then
 if (born_flag.eq.0)  rho = gg/4d0
! write(30,*)'ADOPTED VALUE of RHO = ',rho
! endif
 rho2=rho*rho

 call calculate_volume(r1,r2,r3,om0)

 q2=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
 qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+ &
&       q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+ &
&       q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
 dqeq(1) = 2*epsil(1,1)*q(1)+(epsil(1,2)+epsil(2,1))*q(2)+(epsil(1,3)+epsil(3,1))*q(3)
 dqeq(2) = 2*epsil(2,2)*q(2)+(epsil(2,3)+epsil(3,2))*q(3)+(epsil(2,1)+epsil(1,2))*q(1)
 dqeq(3) = 2*epsil(3,3)*q(3)+(epsil(3,1)+epsil(1,3))*q(1)+(epsil(3,2)+epsil(2,3))*q(2)

 term=1d0/(qeq*eps_scale)/om0*exp(-q2/rho2)

 do na = 1,natoms0
   ma = atom0(na)%mass
 do nb = 1,natoms0
   mb = atom0(nb)%mass
   rr = atompos(:,na)-atompos(:,nb)
   do i=1,3
      zag(i) = q(1)*zeu(1,i,na)+q(2)*zeu(2,i,na)+q(3)*zeu(3,i,na)
      zbg(i) = q(1)*zeu(1,i,nb)+q(2)*zeu(2,i,nb)+q(3)*zeu(3,i,nb)
   end do
   do i = 1,3
   do j = 1,3
      dynmat(i+3*(na-1),j+3*(nb-1)) = dynmat(i+3*(na-1),j+3*(nb-1))+ &
      & zag(i)*zbg(j)*term/sqrt(ma*mb)

     do al=1,3
        ddyn(i+3*(na-1),j+3*(nb-1),al) = ddyn(i+3*(na-1),j+3*(nb-1),al)+term/sqrt(ma*mb) *  &
 &          (zeu(al,i,na)*zbg(j)+zeu(al,j,nb)*zag(i)-dqeq(al)*zag(i)*zbg(j)/qeq)
     enddo

   end do
   end do
 end do
 end do

 end subroutine nonanal
!=========================================================
 subroutine read_born
! reads ! born ! effective ! charge ! and ! dielectric ! constants ! from ! para.born
 use born
 use lattice
 use atoms_force_constants
 use ios
 integer i,j,k

 open(uborn,file='params.born',status='old')
 read(uborn,*) rho,born_flag   ! if 0 use default
 do i=1,3
    read(uborn,*)(epsil(i,j),j=1,3)
 end do
 do k=1,natoms0
    do i=1,3
       read(uborn,*)(zeu(i,j,k),j=1,3)
    end do
 end do

 end subroutine read_born
!============================================================
 subroutine matrix_elt(q1,q2,l1,l2,l3,w33,inside)
! for input modes (qi,li) i=1,2,3 calculates the 3-phonon matrix element w33(1,2,3) on the fly
 use kpoints
 use io2
 use phi3
 use svd_stuff
 use eigen
 use params
 use constants
 use atoms_force_constants
 use force_constants_module
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: l1,l2,l3
 real(8), intent(in) :: q1(3),q2(3)
 integer, intent(out) :: inside
 complex(8), intent (out) :: w33
 integer i3,nk1,nk2,nk3,i0,al,be,ga,j,j0,k,k0,t,ta1,ta2,ta3,i1,j1,k1
 real(8) q3(3),mi,mj,mk,rx2(3),rx3(3),den
 complex(8) xx,eiqr

! write(*,*)'matrix_elt',q1,q2,l1,l2,l3
! write(*,*)'matrix_elt, const33=',const33

!***************************
! here we assume q1 is IBZ (non-shifted) but q2 is FBZ and shifted
!***************************

 call get_k_info_cent(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
 call get_k_info_cent(q2-shft,NC,nk2,i1,j1,k1,g1,g2,g3,inside)
 w33 = cmplx(0d0,0d0)
 q3 = -q2-q1
 call get_k_info_cent(q3-shft,NC,nk3,i1,j1,k1,g1,g2,g3,inside)
! write(*,*)'matrix_elt: q3=',q3,nk3,inside
 xx = cmplx(0d0,0d0)

 tloop: do t=1,nterms(3)

! write(*,*)'matrix_elt: tloop=',t,nterms(3)
       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop ! first index must be in primitive cell
       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j) - atompos(:,j0)
       rx3(:) = atompos(:,k) - atompos(:,k0)
       eiqr = exp( ci* ((q2 .dot. rx2) + (q3 .dot. rx3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
!&      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)
&      eivecibz(ta1,l1,mapibz(nk1))*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)

 enddo tloop
! den = sqrt(8*eigenval(l1,nk1)*eigenval(l2,nk2)*eigenval(l3,nk3))
 den = sqrt(8*eivalibz(l1,mapibz(nk1))*eigenval(l2,nk2)*eigenval(l3,nk3))
 if (den.ne.0) then
! the last term cnst^1.5 is new and corrects for previous cnst*sqrt(eiv)
    w33 = xx / den * const33 *sqrt(cnst*cnst*cnst)  
 else
    write(*,*)'MATRIX_ELT: den=0 ',den
    stop
 endif

 end subroutine matrix_elt
!===========================================================
 subroutine matrix_elt_full(q,q2,q3,omq,om2,om3,evq,ev2,ev3,w33)
! for input modes (qi,li) i=1,2,3 calculates the 3-phonon matrix element w33(1,2,3) on the fly
 use kpoints
 use io2
 use phi3
 use svd_stuff
 use eigen
 use params
 use constants
 use atoms_force_constants
 use force_constants_module
 use om_dos
 use lattice
 implicit none
! integer, intent(in) :: la
 real(8), intent(in) :: q(3),q2(3),q3(3),omq,om2,om3
 complex(8), intent(in) :: evq(ndyn),ev2(ndyn),ev3(ndyn)
! integer, intent(out) :: inside
 complex(8), intent (out) :: w33
 integer i3,nk1,nk2,nk3,al,be,ga,j1,k1,i0,j0,k0,t,ta1,ta2,ta3
 real(8) mi,mj,mk,rx2(3),rx3(3),den
 complex(8) xx,eiqr

! write(*,*)'matrix_elt',q1,q2,l1,l2,l3
! write(*,*)'matrix_elt, const33=',const33
 w33 = cmplx(0d0,0d0)
 xx = cmplx(0d0,0d0)

 tloop: do t=1,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (i0.ne.i3) cycle tloop ! first index must be in primitive cell
       al = ixyzterm_3(1,t)
       be = ixyzterm_3(2,t)
       ga = ixyzterm_3(3,t)
       mi = atom0(i0)%mass
       j1 = iatomterm_3(2,t)  ;  j0 = iatomcell0(j1) ;  mj = atom0(j0)%mass
       k1 = iatomterm_3(3,t)  ;  k0 = iatomcell0(k1) ;  mk = atom0(k0)%mass
       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
       rx2(:) = atompos(:,j1) - atompos(:,j0)
       rx3(:) = atompos(:,k1) - atompos(:,k0)
       eiqr = exp( ci* ((q2.dot. rx2) + (q3 .dot. rx3)) )
       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
&      evq(ta1)*ev2(ta2)*ev3(ta3)
! &      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)

 enddo tloop
 den = sqrt(8*omq*om2*om3/(cnst*cnst*cnst))
 if (den.ne.0) then
    w33 = xx / den * const33
 else
    write(*,*)'MATRIX_ELT: den=0 ',den
    stop
 endif

 end subroutine matrix_elt_full
!===========================================================
 subroutine calculate_w3_ibz(ndn,ni,ki,nk,kp,eival,eivec)
! this subroutine computes the sum_R2,R3,Tau1,Tau2,Tau3,al,be,ga of  (hbar/2)**1.5
! Psi_R2,R3,Tau1,Tau2,Tau3^al,be,ga e^{i k2.R2+i k3.R3} /sqrt(m_tau1*m_tau2*m_tau3) *
! e_la1^{tau1,al}(k1)*e_la2^{tau2,be}(k2)*e_la3^{tau3,ga}(k3)/sqrt(w_la1(k1)*w_la2(k2)*w_la3(k3))
! with k1=-k2-k3+G so that all k1,k2,k3 are within the primitive cell of the recip space i.e.
! between 0 and gi; i=1,2,3 ; units for V3 are in cm^-1
! this works for two k's (k2,k3) in the prim cell on the k-mesh, defined by sum_i=1,3 ni/Ni*gi
! this subroutine calculates V3(k1,la1,k2,la2,k3,la3)  with k1=-k2-k3+G  and all ki in the primitive G-cell
! with the restriction of k2 in the IRREDUCIBLE BZ defined by ki(3,ni) included in kp(3,nk)
!
 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ni,ndn,l1,l2,l3,n2,n3,j,k
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,indx
 real(8) q1(3),q2(3),q3(3),eival(ndn,nk)
 real(8) ki(3,ni),kp(3,nk)
 complex(8) eivec(ndn,ndn,nk),xx

 v33 = cmplx(0d0,0d0)
 indx=0
 write(ulog,*)'V3: nc(123)=',nc

!************************************************
! second argument q2 needs to be only in the IFBZ
!************************************************
 loop2: do n2=1,ni
    q2 = ki(:,n2)  ! should be = kp(:,mapinv(n2))
    call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(mapinv(n2).ne.nk2) then
      write(ulog,*)'n2,mapinv(n2),nk2,inside=',n2,mapinv(n2),nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif

 loop3: do n3=1 ,nk
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif

    q1 = -q2-q3

    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)

  write(*,2)'nk123=',nk1,nk2,nk3

 do l2=1 ,ndn
 do l3=1 ,ndn
 do l1=1 ,ndn
    call matrix_elt(q2,q3,l2,l3,l1,xx,inside)
    indx=indx+1
    v33(indx)=xx
    nq1(indx)= nk1
    nq2(indx)= nk2
    nq3(indx)= nk3
    la1(indx)= l1
    la2(indx)= l2
    la3(indx)= l3
!   if (mod(indx,1000).eq.1) write(*,*) 'indx,v33=',indx,xx
 enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop3
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)
 enddo loop2

 nv3 = indx
 write(ulog,*)' V33: total size of this array is =',nv3

 if( readv3.eq.0) then
! k1 if( writev3.eq.1) then
   open(uv3,file='v33.dat',status='unknown',FORM='UNFORMATTED')
   write(uv3)nv3
!  open(uv3,file='v33.dat',status='unknown')
!  write(uv3,*)nv3
   do j=1 ,nv3
      write(uv3) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
!     write(uv3,7) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
   enddo
   close(uv3)
 endif

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))

 end subroutine calculate_w3_ibz
!============================================================
 subroutine function_self_w(q,la,omega,temp,nself,uself)
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1,nk3,ins2,ifbz,l1,l3
 integer, save :: cntr
 real(8) om1,eta,tk,term,k1(3),k2(3),q1(3),q3(3)
 complex(8) omz,self,ocfunc,oc2,s2,xx

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
! write(*,*) 'self_w: tempk=',tk
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! write(*,*) 'self_w: k_info: nq=',nq
! om1 = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz
! write(ulog,5)'SELF_W: last cntr, omega,etaz=',cntr,omega,etaz
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0

 cntr = 0
 do ifbz=1,nkc
    q3=kpc(:,ifbz)
    call get_k_info(q3,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
    if (nk3.ne.ifbz) then
       write(ulog,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       write(*,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       stop
    endif
    k1 = -q-q3
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,ins2)
!   if (ins2.ne.inside) then
!      write(ulog,*)'SELF_W: ins2 ne inside ',ins2,inside
!      stop
!   endif

    do l3=1,ndyn
    do l1=1,ndyn
!      write(*,*) 'self_w: enetering matrix_elt',q,q3,la,l3,l1
       call matrix_elt(q,q3,la,l3,l1,xx,inside)
       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)

       jk1 = nk1+(l1-1)*nkc
       jk3 = ifbz+(l3-1)*nkc

!      cntr = cntr+1
       self = ocfunc(temp,omz,l1,nk1,l3,nk3) !jk1,jk3)
 !     s2 = oc2(temp,omz,jk1,jk3)
 !     if (aimag(s2).gt.1/6/eta) write(222,6)ifbz,q3,q,la,l3,l1
       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif
    enddo
    enddo

 enddo
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
6 format(i6,6(1x,f9.4),2x,3(1x,i2))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w
!============================================================
 subroutine function_self_w2(q,la,omega,temp,nself,uself)
! calculates the self-energy on the fly for q-point in the generated kmesh and Lorentzian delta
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer iq,jq,kq,inside,ik,jk,kk,j,nq,jk1,jk2,jk3,nk1,nk3,ins2,ifbz,l1,l3
 integer, save :: cntr
 real(8) omq,om3,omk,eta,tk,term,k1(3),k2(3),q1(3),q3(3)
 complex(8) omz,self,ocfunc,oc2,s2,xx

 tk = temp*100*c_light*h_plank/k_b   ! 1.44 K = 1 /cm
! write(*,*) 'self_w: tempk=',tk
 call get_k_info(q,nc,nq,iq,jq,kq,g1,g2,g3,inside)
! write(*,*) 'self_w: k_info: nq=',nq
 omq = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz
! write(ulog,5)'SELF_W: last cntr, omega,etaz=',cntr,omega,etaz
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0

 cntr = 0
 do ifbz=1,nkc
    q3=kpc(:,ifbz)
    call get_k_info(q3,nc,nk3,ik,jk,kk,g1,g2,g3,inside)
    if (nk3.ne.ifbz) then
       write(ulog,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       write(*,*)'SELF_W: nk3 ne ifbz ',nk3,ifbz
       stop
    endif
    k1 = -q-q3
    call get_k_info(k1,nc,nk1,ik,jk,kk,g1,g2,g3,ins2)

    do l3=1,ndyn
       om3 = sqrt(abs(eigenval(l3,nk3))) * cnst
    do l1=1,ndyn
       omk = sqrt(abs(eigenval(l1,nk1))) * cnst
!      call matrix_elt     (q,q3,la,l3,l1,xx,inside)
       call matrix_elt_full(q,q3,k1,omq,om3,omk,eigenvec(:,la,nq),eigenvec(:,l3,nk3),eigenvec(:,l1,nk1),xx)
       call check_inside_bz(k1,g1,g2,g3,inside)
       if(inside.ne.ins2) then
          write(*,*)'SELF_W:ERROR:inside .ne. ins2 ',inside,ins2,k1
          stop
       endif
       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)

       jk1 = nk1+(l1-1)*nkc
       jk3 = ifbz+(l3-1)*nkc

       self = ocfunc(temp,omz,l1,nk1,l3,nk3)
       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif
    enddo
    enddo

 enddo
 nself = nself /nkc /2 *2*pi
 uself = uself /nkc /2 *2*pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
6 format(i6,6(1x,f9.4),2x,3(1x,i2))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w2
!============================================================
 subroutine function_self_w3(q,la,omega,temp,nself,uself)
! this one is for arbitrary q
! uses lorentzian delta with 4 terms, and calculates both real and imaginary parts
! In this subroutine, which calculates like self_w on the fly, the second momentum
! is restricted to be on the kpc mesh, but the first and third momenta can be anywhere
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer inside,ik,l3,l2
! integer, save :: cntr
 real(8) omq,eta,term,k3(3),k2(3),eivq(ndyn),eiv2(ndyn),eiv3(ndyn),vg(3,ndyn)
 real(8) etacut,arg,delta_l,om3,om2,nb3,nb2,nbe,ires,rres
 complex(8) omz,self,xx,evq(ndyn,ndyn),ev2(ndyn,ndyn),ev3(ndyn,ndyn)

! be careful band sorting should be avoided in this case; otherwise need to
! find the nearest kpoint in kpc, and sort eivq according to it
 call get_freq(q,ndyn,vg,eivq,evq)
 omq=eivq(la)
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0
 etacut = eta*3000 ! error due to non-zero eta is 1/9E6

 do ik=1,nkc
    k2(:)=kpc(:,ik)
    ev2=eigenvec(:,:,ik)
    k3 = -q-k2
    call check_inside_bz(k3,g1,g2,g3,inside)

    call get_freq(k3,ndyn,vg,eiv3,ev3)

    do l2=1,ndyn
       om2=eigenval(l2,ik)
       nb2=nbe(om2,temp,classical)
    do l3=1,ndyn

! first check see if the energy is conserved
       om3=eiv3(l3)
       nb3=nbe(om3,temp,classical)
       ires=0; rres=0

       arg=(omega-om3-om2)
       if (abs(arg) .lt. etacut) then
          ires=(nb3+nb2+1)*delta_l(arg,eta)
          rres=(nb3+nb2+1)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega+om3+om2)
       if (abs(arg) .lt. etacut) then
          ires=ires-(nb3+nb2+1)*delta_l(arg,eta)
          rres=rres-(nb3+nb2+1)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega-om3+om2)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb2-nb3)*delta_l(arg,eta)
          rres=rres+(nb2-nb3)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega-om2+om3)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb3-nb2)*delta_l(arg,eta)
          rres=rres+(nb3-nb2)*arg/(arg*arg+eta*eta)
       endif

       self=cmplx(rres,ires*pi)
 !     if (abs(self).lt.1d-6) cycle loop1   ! self is in cm

       call matrix_elt_full(q,k2,k3,omq,om2,om3,evq(:,la),ev2(:,l2),ev3(:,l3),xx)

       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)

       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif

 enddo
 enddo
 enddo
 nself = nself /nkc *pi
 uself = uself /nkc *pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w3
!============================================================
 subroutine function_self_w4(q,la,omega,temp,nself,uself)
! this one is for arbitrary q
! uses gaussian delta with 3 terms, and just calculates the imaginary part.
! In this subroutine, which calculates like self_w on the fly, the second momentum
! is restricted to be on the kpc mesh, but the first and third momenta can be anywhere
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use io2
 use phi3
 use eigen
 use params
 use constants
 use om_dos
 use lattice
 implicit none
 integer, intent(in) :: la
 real(8), intent(in) :: q(3),temp,omega
 complex(8), intent(out) :: nself,uself
 integer inside,ik,l3,l2
! integer, save :: cntr
 real(8) omq,eta,term,k3(3),k2(3),eivq(ndyn),eiv2(ndyn),eiv3(ndyn),vg(3,ndyn)
 real(8) etacut,arg,delta_g,om3,om2,nb3,nb2,nbe,ires,rres
 complex(8) omz,self,xx,evq(ndyn,ndyn),ev2(ndyn,ndyn),ev3(ndyn,ndyn)

! be careful band sorting should be avoided in this case; otherwise need to
! find the nearest kpoint in kpc, and sort eivq according to it
 call get_freq(q,ndyn,vg,eivq,evq)
 omq=eivq(la)
 eta= etaz !min(etaz,1000/omega) ! (deltak*vgmax)**2/omega !
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0
 etacut = eta*6 ! error due to non-zero eta is 1/9E6

 do ik=1,nkc
    k2=kpc(:,ik)
    ev2=eigenvec(:,:,ik)
    k3 = -q-k2
    call check_inside_bz(k3,g1,g2,g3,inside)

    call get_freq(k3,ndyn,vg,eiv3,ev3)

    do l2=1,ndyn
       om2=eigenval(l2,ik)
       nb2=nbe(om2,temp,classical)
    do l3=1,ndyn

! first check see if the energy is conserved
       om3=eiv3(l3)
       nb3=nbe(om3,temp,classical)
       ires=0; rres=0

       arg=(omega-om3-om2)
       if (abs(arg) .lt. etacut) then
          ires=(nb3+nb2+1)*delta_g(arg,eta)
!         rres=(nb3+nb2+1)*arg/(arg*arg+eta*eta)
       endif

! this affects the lifetime of low-frequency phonons, can be removed...
! contributes only to noremal processes, and might yield a higher power of w!
!       arg=(omega+om3+om2)
!       if (abs(arg) .lt. etacut) then
!          ires=ires-(nb3+nb2+1)*delta_g(arg,eta)
!!         rres=rres-(nb3+nb2+1)*arg/(arg*arg+eta*eta)
!       endif

       arg=(omega-om3+om2)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb2-nb3)*delta_g(arg,eta)
!         rres=rres+(nb2-nb3)*arg/(arg*arg+eta*eta)
       endif

       arg=(omega-om2+om3)
       if (abs(arg) .lt. etacut) then
          ires=ires+(nb3-nb2)*delta_g(arg,eta)
!         rres=rres+(nb3-nb2)*arg/(arg*arg+eta*eta)
       endif

       self=cmplx(rres,ires*pi)
 !     if (abs(self).lt.1d-6) cycle loop1   ! self is in cm

       call matrix_elt_full(q,k2,k3,omq,om2,om3,evq(:,la),ev2(:,l2),ev3(:,l3),xx)

       if (abs(xx).lt.v3_threshold) cycle
       term =xx*conjg(xx)

       if(inside.eq.1) then                ! normal process
          nself = nself + term * self
       else                                ! umklapp process
          uself = uself + term * self
       endif

 enddo
 enddo
 enddo
 nself = nself /nkc *pi
 uself = uself /nkc *pi
 self  = nself + uself
! 2 pi is necessary because |V|^2 /(hbar*om) is in cm^-1 ; thus we
! need to multiply it by 100*c_light*h_plank to convert it to J, but
! there is also another extra factor of 1/hbar in the denominator, which
! multiplied by h_plank leaves an extra 2*pi and dividing by the remaining
! 100*c_light converts the inverse time (in 1/s) to inverse length (in 1/cm)
!
! 2pi for om to Hertz, c_light for Hertz to m^-1, 100 for m^-1 to cm^-1
! for meV, need to multiply Hz by h_plank/ee*1000
! write(uslf,3)nq,la,kpc(:,nq),omega,self, nself,uself
3 format(a,i6,i6,3x,3(1x,f8.3),1x,g11.5,1x,9(2x,2(1x,g10.4)))
4 format(a,9(1x,g10.4))
5 format(a,i8,99(1x,g10.4))
7 format(a,4i3,99(1x,f7.3))

 end subroutine function_self_w4
!============================================================
 subroutine phase(ndn,ni,ki,nk,kp,eival,eivec)
! this subroutine computes for each q the sum_q2 delta(w-w_q2 \pm w_q3)
! and finds the kpoints available
 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
 use svd_stuff
 use params
 use constants
 implicit none
 integer nk,ni,ndn,l1,l2,l3,n2,n3,j,k
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk3,indx,igam
 real(8) q1(3),q2(3),q3(3),eival(ndn,nk)
 real(8) ki(3,ni),kp(3,nk)
 complex(8) eivec(ndn,ndn,nk),xx

 v33 = cmplx(0d0,0d0)
 indx=0
 write(ulog,*)'V3: nc(123)=',nc

    q1=0
    call get_kindex(q1,nk,kp,igam)
    q1 = kp(:,igam+3)

 loop3: do n3=1 ,nk
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif

!! properly initilize q2 !!!!!!!!
    q1 = -q2-q3

    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)

  write(*,2)'nk13=',nk1,nk3

! do l2=1 ,ndn
! do l3=1 ,ndn
! do l1=1 ,ndn
!    call matrix_elt(q2,q3,l2,l3,l1,xx,inside)
! enddo
! enddo
! enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop3
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)

 nv3 = indx
 write(ulog,*)' V33: total size of this array is =',nv3

 if( readv3.eq.0) then
! k1 if( writev3.eq.1) then
   open(uv3,file='v33.dat',status='unknown',FORM='UNFORMATTED')
   write(uv3)nv3
   do j=1 ,nv3
      write(uv3) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
   enddo
   close(uv3)
 endif

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))

 end subroutine phase
!===========================================================
 subroutine get_freq(kp,ndn,vg,eival,eivec)
 use params
 use io2
 use constants
 use born
 implicit none
 integer, intent(in) :: ndn
 real(8), intent(in) :: kp(3)
 real(8), intent(out):: eival(ndn)
 complex(8), intent(out) :: eivec(ndn,ndn)
 integer i,j,k,l,ier,nd2,al
 integer, allocatable :: mp(:)
 real(8), allocatable :: eivl(:)
 real(8) absvec,vg(3,ndn)
 complex(8), allocatable:: dynmat(:,:),eivc(:,:),ddyn(:,:,:),temp(:,:)
 real(8) khat(3),margn
 integer ndeg,lower,upper,nv

 nd2 = min(ndn,12)
 margn=1.d-4 
 nv=0

 allocate(dynmat(ndn,ndn),mp(ndn),eivl(ndn),eivc(ndn,ndn),ddyn(ndn,ndn,3),temp(ndn,ndn))
! write(uio,*)'# dynmat allocated, printing: Kpoint, eival(j),j=1,ndyn)'
! write(uio,*)'############ before setting up dynmat, kp number ',i

    call set_dynamical_matrix(kp,dynmat,ndn,ddyn)

!JS: call nonanalytical term for born effective change
    if (kp(1)==0.d0 .AND. kp(2)==0.d0 .AND. kp(3)==0.d0) then
        khat=kp(:)+1.0D-10
    else
        khat=kp(:)
    endif
    call nonanal(khat,dynmat,ndn,ddyn)

    if (verbose) then
       write(ulog,*)' ======================================================================'
       write(ulog,3)' THE DYNAMICAL MATRIX for ndyn, KP=',ndn,kp(:)
       do l=1,ndn
          write(ulog,8)(dynmat(l,j),j=1,nd2)
       enddo
       write(ulog,*)' ======================================================================'
       do al=1,3
          write(ulog,*)' THE k-DERIVATIVE OF THE DYNAMICAL MATRIX =',al
          do l=1,ndn
             write(ulog,8)(ddyn(l,j,al),j=1,nd2)
          enddo
       enddo
    endif

! store dynmat in temp
    temp=dynmat
    call diagonalize(ndn,dynmat,eivl,ndn,eivc,ier)

! sort eivals in ascending order
    call sort(ndn,eivl,mp,ndn)
!   write(ulog,6)'eival sorting_map=',mp

    if (ier.ne.0 .or. verbose) then
      write(   *,*)' ier=',ier
      write(ulog,*)' ier=',ier, ' EIGENVALUES ARE...'
      write(ulog,4) eivl
    endif

    do j=1,ndn
       eival(j) = cnst*sqrt(abs(eivl(mp(j)))) ! so that all frequencies are positive
!      eival(j) = abs(eivl(mp(j))) ! so that all frequencies are positive
       do l=1,ndn
          eivec(l,j) = eivc(l,mp(j))
       enddo
    enddo

! if pure imaginary make eigenvectors real as they can be multiplied by any arbitrary number
!   do l=1,ndn
!      absvec = sum(abs(real(eivec(:,l))))
!      if (absvec .lt. 1d-3) then
!         eivec(:,l)=cmplx(0,1)*eivec(:,l)
!      endif
!   enddo

  if (verbose) then

     do l=1,ndn
        write(ulog,*)' GET_FREQ:-----------  BAND ',l,eival(l)
        do j=1,ndn/3
           write(ulog,9)'atom ',j,eivec(3*(j-1)+1,l),eivec(3*(j-1)+2,l),eivec(3*(j-1)+3,l)
        enddo
     enddo

! now calculate the square of frequencies based on dynmat
     dynmat=temp
     temp = matmul(dynmat,eivec)
     temp = matmul(transpose(conjg(eivec)),temp)
     write(ulog,*)' TEST: below square of eivals should appear in the diagonal if e.D.e is right'
     do j=1,ndn
       write(ulog,9)'e.D.e,eivl=',j, temp(j,j),eivl(mp(j))
     enddo
!    write(ulog,*)' three components of ddyn are '
!    do al=1,3
!    do j=1,ndn
!      write(ulog,9)'al=',al, ddyn(j,:,al)
!    enddo
!    enddo

  endif

 deallocate(eivl,eivc,dynmat,mp,ddyn,temp)
 return

     do al=1,3
        temp = matmul(transpose(conjg(eivec)),matmul(ddyn(:,:,al),eivec))
!       write(ulog,*)' ======================================================================'
!       write(ulog,3)' DDYN MATRIX in the eivec basis for alpha, KP=',al,kp(:)
!       do l=1,ndn
!          write(ulog,8)(temp(l,j),j=1,ndn)
!       enddo
        call diagonalize(ndn,temp,vg(al,:),nv,eivc,ier)
!       write(ulog,*)' ======================================================================'
!       write(ulog,*)' THE eigenvalues of the above matrix =',al
        vg(al,:)=vg(al,:)/2/eival(:)*cnst*cnst*1d-10*100*2*pi
!       write(ulog,8) (vg(al,l),l=1,ndn)
     enddo

!------------------------------------------------------------------------------
! now calculate the group velocities from Hellman-Feynmann formula(for each component al=1,3)

! identify the degenerate states and diagonalize ddyn in that space to get group velocities
    lower=1; upper=1
    do while(upper.le.ndn)

! form and diagonalize ddyn from lower to upper
       if (eival(upper)-eival(lower) .gt. margn*eival(upper) ) then !diagonalize the (lower:upper-1) degenerate block

          if(verbose) then
             write(ulog,2)'GET_FREQ:lower,upper-1,eival=',lower,upper-1,eival(lower:upper-1)
          endif
          deallocate (temp,eivc)
          ndeg=upper-lower
          allocate(temp(ndeg,ndeg),eivc(ndeg,ndeg))
          do al=1,3
             temp = matmul(transpose(conjg(eivec(:,lower:upper-1))),matmul(ddyn(:,:,al),eivec(:,lower:upper-1)))
             call diagonalize(ndeg,temp,vg(al,lower:upper-1),nv,eivc,ier)
             vg(al,lower:upper-1)=vg(al,lower:upper-1)/2/eival(lower:upper-1)*cnst*cnst*1d-10*100*2*pi
          enddo

          if (upper.ne.ndn) then ! define the indices of the next block
             lower=upper
             upper=upper+1
          else ! upper=ndn: this is the last eigenvalue, last block

             if(verbose) then
                write(ulog,2)'GET_FREQ:lower,upper  ,eival=',upper,upper,eival(upper:upper)
             endif
             deallocate (temp,eivc)
             ndeg=1
             allocate(temp(ndeg,ndeg),eivc(ndeg,ndeg))
             do al=1,3
                temp = matmul(transpose(conjg(eivec(:,upper:upper))),matmul(ddyn(:,:,al),eivec(:,upper:upper)))
                call diagonalize(ndeg,temp,vg(al,upper:upper),nv,eivc,ier)
                vg(al,upper:upper)=vg(al,upper:upper)/2/eival(upper:upper)*cnst*cnst*1d-10*100*2*pi
             enddo
             exit

          endif
       else  ! they are degenerate, go to the next eival
          if (upper.lt.ndn) then
             upper=upper+1
          else  ! reached the end of array

             if(verbose) then
                write(ulog,2)'GET_FREQ:lower,upper  ,eival=',lower,upper,eival(lower:upper)
             endif
             deallocate (temp,eivc)
             ndeg=upper-lower+1
             allocate(temp(ndeg,ndeg),eivc(ndeg,ndeg))
             do al=1,3
                temp = matmul(transpose(conjg(eivec(:,lower:upper))),matmul(ddyn(:,:,al),eivec(:,lower:upper)))
                call diagonalize(ndeg,temp,vg(al,lower:upper),nv,eivc,ier)
                vg(al,lower:upper)=vg(al,lower:upper)/2/eival(lower:upper)*cnst*cnst*1d-10*100*2*pi
             enddo
             exit

          endif
       endif

    enddo


    if(verbose) then
      write(ulog,3)'GET_FREQ: vgx_hf=',ndn,(vg(1,l)*c_light,l=1,ndn)
      write(ulog,3)'GET_FREQ: vgy_hf=',ndn,(vg(2,l)*c_light,l=1,ndn)
      write(ulog,3)'GET_FREQ: vgz_hf=',ndn,(vg(3,l)*c_light,l=1,ndn)
    endif

!   do al=1,3
!   do l=1,ndn

!         dynmat=0
!         do l=1,ndn
!         do k=1,ndn
!            do j=1,ndn
!               dynmat(k,l)=dynmat(k,l)+ddyn(k,j,al)*eivec(j,l)
!            enddo
!         enddo
!         enddo

!         do l=1,ndn
!           vg(al,l)=0
!           do k=1,ndn
!              vg(al,l)=vg(al,l)+dynmat(k,l)*conjg(eivec(k,l))
!           enddo
!           vg(al,l)=vg(al,l)/2/sqrt(abs(eivl(mp(l))))*cnst*1d-10*100*2*pi
!         enddo

!   enddo
!   enddo

 deallocate(eivl,eivc,dynmat,mp,ddyn,temp)

 2 format(a,2(1x,i5),9(1x,f9.3))
 3 format(a,i5,99(1x,f11.4))
 4 format(99(1x,2(1x,g9.3),1x))
 5 format(a,i5,99(1x,2(1x,g9.3),1x))
 6 format(a,99(1x,i5))
 7 format(i6,2x,3(1x,f7.3),1x,99(1x,f7.3))
 8 format(99(1x,2(1x,f7.3),1x))
 9 format(a,i5,99(1x,2(1x,f9.3),1x))

 end subroutine get_freq
!===========================================================
 subroutine enforce_asr_simple(n,dyn)
 implicit none
 integer n,nat,i,j,io,jo,lo,iat,jat
 complex(8) dyn(n,n),sum

 nat=n/3
 do i=1,3
 do iat=1,nat
    io=i+3*(iat-1)
    do j=1,3
       sum=0
       jo=j+3*(iat-1)
       do jat=1,nat
          if (iat.eq.jat) cycle
          lo=j+3*(jat-1)
          sum=sum-dyn(io,lo)
       enddo
       dyn(io,jo)=sum
    enddo
 enddo
 enddo
 end subroutine enforce_asr_simple
!===========================================================
subroutine band_sort_bs(nkp,ndyn,kp,dk,eival,eivec,map)
implicit none

integer, intent(in) :: nkp,ndyn
real(8), intent(in) :: kp(3,nkp)
real(8), intent(in) :: eival(ndyn,nkp)
complex(8), intent(in) :: eivec(ndyn,ndyn,nkp)
integer, intent(out) :: map(ndyn,nkp)

integer i,j,k,l
integer ntbd !Number of connections to be determined

real(8) bmax !Maximum band derivative max(dE/dk)
real(8) dk(nkp)
real(8) fsaf !Safety factor

complex(8), allocatable :: overlap(:,:)
real(8), allocatable :: overlap_q(:,:)

!First estimate the maximum band derivative bmax=max(dE/dk)
!bmax=0
!fsaf=10
!allocate(dk(2:nkp))
!do i=2,nkp
!   dk(i)=sqrt((kp(1,i)-kp(1,i-1))**2+(kp(2,i)-kp(2,i-1))**2+(kp(3,i)-kp(3,i-1))**2)
!   do j=1,ndyn
!      bmax = max(bmax,abs(eival(j,i)-eival(j,i-1))/dk(i))
!   end do
!end do

!Start from the first k-point, which will be used as the reference for other k points
do i=1,ndyn
    map(i,1)=i
end do

!from the second k-point to the last, calculate the overlap matrix and do band sorting
allocate (overlap(ndyn,ndyn),overlap_q(ndyn,ndyn))

do i=2,nkp
    ntbd=ndyn
    overlap=matmul(transpose(dconjg(eivec(:,:,i-1))),eivec(:,:,i))
    do j=1,ndyn
    do k=1,ndyn
        overlap_q(j,k)=overlap(j,k)*dconjg(overlap(j,k))
    end do
    end do

!Step 0: if two bands with energy difference larger than fsaf*bmax, the bands are not supposed to be connected
!    write(*,*) "energy window for band sorting:", fsaf*bmax*dk(i)
!    do j=1,ndyn-1
!    do k=j+1,ndyn
!       if (abs(eival(j,i-1)-eival(k,i)) .ge. fsaf*bmax*dk(i)) then
!          overlap_q(j,k)=0
!       end if
!    end do
!    end do

!Step 1: find matrix elements larger than 0.6, which indicates a strong connectivity
    do j=1,ndyn
    do k=1,ndyn
        if (overlap_q(j,k) .ge. 0.6) then
            do l=1,ndyn
                if ( map(l,i-1) .eq. j) then
                    map(l,i)=k
                end if
            end do
                ntbd=ntbd-1
                overlap_q(j,:)=0  !Connected ones will not be examined again
                overlap_q(:,k)=0
        end if
    end do
    end do

    if (ntbd .ne. 0) then !If there are unresolved connections remaining
!Step 2: find the largest remaining matrix elements in each row, and connect eigenvalues connected by that.
    do j=1,ndyn
        if (maxval(overlap_q(j,:)) .ge. 0.1) then
kkloop:     do k=1,ndyn
               if (overlap_q(j,k) .eq. maxval(overlap_q(j,:))) then
                  if (overlap_q(j,k) .lt. 0.3) then
                     write(*,*) "Warning: the overlap matrix element is &
   & smaller than 0,3, there may be ambiguity for band connectivity"

 write(*,*) "k-point:",kp(:,i),"index of kp:",i,"Matrix element:",overlap_q(j,k)
                  end if
                  do l=1,ndyn
                     if ( map(l,i-1) .eq. j) then
                        map(l,i)=k
                     end if
                  end do
                  ntbd=ntbd-1
                  overlap_q(j,:)=0
                  overlap_q(:,k)=0
                  exit kkloop
               end if
            end do kkloop
        end if
    end do
    end if

    if (ntbd .ne. 0) then !If there are still unresolved connections remaining
!Step 3: connect the remaining pairs in a ascending order
    do j=1,ndyn
    do k=1,ndyn
        if (overlap_q(j,k) .ne. 0) then
write(*,*) "Warning: the overlap matrix element is smaller than 0.3, &
& there may be ambiguity to determine the band connectivity"
write(*,*) "k-point:",kp(:,i),"index of kp:",i,"Matrix element:",overlap_q(j,k)
           do l=1,ndyn
              if ( map(l,i-1) .eq. j) then
                 map(l,i)=k
              end if
           end do
        ntbd=ntbd-1
        overlap_q(j,:)=0
        overlap_q(:,k)=0
        end if
    end do
    end do
    end if

        write(*,*) 'map for',i
        write(*,*) map(:,i)

    if (ntbd .ne. 0) then
        write(*,*) 'error: band sorting failed'
        write(*,*) 'k-point:',i,kp(:,i),ntbd
    end if
end do

deallocate(overlap,overlap_q)

end subroutine band_sort_bs
!=========================================
 function mysqrt(x) result(y)
 implicit none
 real(8), intent(in) :: x
 real(8) y  !, intent(out):: y

 if(x.ge.0) then
   y=sqrt(x)
 else
   y=-sqrt(-x)
 endif

 end function mysqrt
!=========================================
 subroutine get_v3_indices(i,ni,nk,nd1,nd2,nd3,n1,n2,l1,l2,l3)
! if index i is generated from n1=1,ni;n2=1,nk,l1=1,nd1;l2=1,nd2;l3=1,nd3
! this routine produces the values of n1,n2,l1,l2,l3
 implicit none
 integer, intent(in) :: i,ni,nk,nd1,nd2,nd3
 integer, intent(out):: n1,n2,l1,l2,l3
 integer i4,i3,i2,i1

 call remainder(i ,nd3,l3,i1)

 call remainder(i1+1,nd2,l2,i2)
 call remainder(i2+1,nd1,l1,i3)
 call remainder(i3+1,nk,n2,n1)
! call remainder(i4,nk,n2,i4)
 n1=n1+1
! l3=mod(i,nd3)

! i3=i-l3
! l2=mod((i3)/nd3,nd2)+1

! i2=i3-l2*nd3
! l1=mod(i2/nd2,nd1)+1

! i1=i2-l3*nd2*nd3
! n2=mod(i1/nd3,nk)+1

! i0=i1-n2*nd2*nd3*nk
! n1=mod(i0/nk,ni)+1

 end subroutine get_v3_indices
!=========================================
 subroutine remainder(i,nd,l,j)
! calculates for the general index i, in blocks of nd, the block index j and the remainder in the block l
 implicit none
 integer i,j,l,nd

 l=1+mod(i-1,nd)
 j=(i-l)/nd

 end subroutine remainder
!========================================= 
 subroutine right_hand_side(n,b,temp)
! this subroutine computes the right hand side of the BE: rhs=(dn/dt)collisions=-dn/tau
! it is given by rhs(k,alpha)=v_k^alpha dn_k/dT
! the first third of the components give the response to dT/dx, dT/dy, and then to dT/dz
! units are 1/temp, i.e. in cm
 use eigen
 use kpoints
 use params
 use constants
 implicit none
 integer n,i,j,nk,la,al !,onedim
 real(8) temp,nbe,nbx,om,dndt,b(n,3)

 if(n.ne.ndyn*nkc) then
    write(*,*)'RHS: input size n and actual dimension of b dont match ',n,ndyn*nkc
    stop    
 endif
 do al=1,3    ! fix this part
    i=0
    do nk=1,nkc
    j=mapinv(nk)
    do la=1,ndyn
       i=onedim(la,j) !i+1
       om=eigenval(la,j)
       nbx=nbe(om,temp,classical)
       dndt=nbx*(nbx+1)*om/temp/temp  ! careful of the units (in cm) 
       b(i,al)=veloc(al,la,j)*dndt*c_light      
    enddo
    enddo
 enddo

 end subroutine right_hand_side
!========================================= 
 subroutine read_collision(temp)  !ncol,nlin,coll,temp)
! reads the diagonal and non-diagonal part of all the collision matrices written into files
! in split mode, matrices are nkc lines by ksubsize columns 
 use params
 use exactBTE2
 use io2
 use eigen     ! for ndyn 
 use kpoints   ! for nkc 
 use phi3      ! for v3path
 implicit none
 integer uq,i,ii,li,lli,j,ksub(2),ks(2),cnti,dk,nibz_proc,nklast_cpu,icpu
 real(8) temp,tmp,om
 character fn*99,fn2*99,cibz1*4,cibz2*4,cnk*6,line*99

 if ( ncpu.eq.1) then
    nibz_proc=nibz  ! number of IBZ kpoints per cpu, except for the last one
 elseif (mod(nibz,ncpu) .eq. 0) then
    nibz_proc=nibz/ncpu  ! number of IBZ kpoints per cpu, except for the last one
 else
    nibz_proc=int(nibz/real(ncpu))+1  ! number of IBZ kpoints per cpu, except for the last one
 endif
 nklast_cpu=nibz-(ncpu-1)*nibz_proc ! number of IBZ kpoints for the last cpu < nibz_proc
 write(ulog,*)ncpu,' processors each having ',nibz_proc,' kpoints, except for'
 write(ulog,*)'the last one only getting ',nklast_cpu

 call allocate_iter0(nibz,nkc,ndyn)
 
do icpu=1,ncpu

    ksub(1)=(icpu-1)*nibz_proc+1
    if ( i.ne.ncpu ) then
       ksub(2)=nibz_proc*icpu 
    else
       ksub(2)=nibz
       if (nklast_cpu.ne.ksub(2)-ksub(1)+1) then
          write(ulog,*)'READ_COLL: inconsistency in last ksubs, nklast_cpu=',nklast_cpu 
          write(ulog,*)'READ_COLL: ksub(1),ksub(2)=',ksub
          stop
       endif
    endif

 write(cibz1,'(i4.4)') ksub(1)
 write(cibz2,'(i4.4)') ksub(2)
 write(cnk,'(i6.6)') nkc
 fn=trim(v3path)//'qval-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
 fn2=trim(v3path)//'coll-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
 uQ=9402
 open(uQ,file=fn,status='old')
 open(uQ+1,file=fn2,status='old')
 read(uQ  ,*)ks,tmp
 if(ks(1).ne.ksub(1) .or. ks(2).ne.ksub(2) ) then
    write(ulog,*)'READ_COLLISION: ksub in call and ksub read do not match ',ksub,ks
    stop
 elseif( abs(temp-tmp) .gt.1d-4) then
    write(ulog,*)'READ_COLLISION: temp in call and temp read do not match ',temp,tmp
    stop
 endif
 read(uQ+1,*)ks,tmp
 if(ks(1).ne.ksub(1) .or. ks(2).ne.ksub(2) ) then
    write(ulog,*)'READ_COLLISION: ksub in call and ksub read do not match ',ksub,ks
    stop
 elseif( abs(temp-tmp) .gt.1d-4) then
    write(ulog,*)'READ_COLLISION: temp in call and temp read do not match ',temp,tmp
    stop
 endif
 read(uQ,*)line


 dk=ksub(2)-ksub(1)+1

 cnti=0
 do i=ksub(1),ksub(2)
 do li=1,ndyn
    cnti=cnti+1
    if( onedim(li,i-ksub(1)+1).ne.cnti ) then
       write(ulog,*)'READ_COLL: 1d(li,i),cnti =',onedim(li,i),cnti
    endif
!   read(uq,*)ii,lli,om,qvalue_n(ii,lli),qvalue_u(ii,lli)
!   qval(cnti)=qvalue_n(ii,lli)+qvalue_u(ii,lli)
    read(uq,*)ii,lli,om,qval(cnti) 
 enddo
 enddo
 do ii=1,nkc*ndyn
    read(uq+1,*)j, (mval(j,i),i=1,ndyn*(ksub(2)-ksub(1)+1))
 enddo

 close(uq); close(uq+1)

enddo

 end subroutine read_collision
!========================================= 
 subroutine calculate_sub_collision(temp,ksub)
! for a given k subset, calculates the collision matrix into a
! single array (nkc*ndyn by nksub*ndyn) for the temperature of interest
! A_k;ij = 2pi/hbar^2 delta(wk-wi-wj) v3sq(i,j,k) (nk+1)*ni*nj
 use kpoints
 use eigen
 use params
 use exactBTE2
 use phi3
 use om_dos
 use constants
 use io2
 implicit none
 integer i,ii,mi,j,jj,li,lj,lk,nk,i1,j1,k1,inside,ksub(2),cnti,cntj
 real(8) qi(3),qj(3),qk(3),temp,q_12,m_12,v33_2
 real(8) nbe, delta_l,nbi,nbj,nbl,oi,oj,ok,qn,qu,n_rate,u_rate,tau_ps,rate
 integer uQ
 character fn*99,fn2*99,cibz1*4,cibz2*4,cnk*6

 write(cibz1,'(i4.4)') ksub(1)
 write(cibz2,'(i4.4)') ksub(2)
 write(cnk,'(i6.6)') nkc
 fn=trim(v3path)//'qval-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
 fn2=trim(v3path)//'coll-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
 uQ=9402
 open(uQ,file=fn,status='unknown')
 open(uQ+1,file=fn2,status='unknown')
 write(uQ  ,*)ksub,temp
 write(uQ+1,*)ksub,temp
 write(uQ,*)'#nk,la,omega,Qval,N-rate,U-rate,total(1/cm),tau(ps),la(nm)'

 call allocate_iter0(ksub(2)-ksub(1)+1,nkc,ndyn)

   if (usetetra.eq.1) then  ! use tetrahedron method



   else  ! use lorentzian or gaussian broadening
      
! first the diagonal elements Q_value are calculated
      Qvalue_N=0; Qvalue_U=0; Mval=0
      cnti=0  
      do ii=ksub(1),ksub(2)
         i=mapinv(ii)
         qi=kpc(:,i)  !kibz(:,i)
         mi=ii-ksub(1)+1
      do li=1,ndyn
!        oi=eivalibz(li,mapinv(i))   ;  nbi=nbe(oi,temp,classical)
         oi=eigenval(li,i)           ;  nbi=nbe(oi,temp,classical)
         cnti=cnti+1

         cntj=0
         do jj=1,nkc
! we introduce j=mapinv(jj) so that j=the first nibz kpc's are in the IBZ 
! because for the collision matrix, unlike v33, the order in kpc must be according to IBZ, and then the rest
            j=mapinv(jj)   ! kibz(jj)=kpc(mapinv(jj))
            qj=kpc(:,j)
            qk=-qi-qj        
            call get_k_info_cent(qk-shft,NC,nk,i1,j1,k1,g1,g2,g3,inside)
         do lj=1,ndyn
            oj=eigenval(lj,j )          ;  nbj=nbe(oj,temp,classical)
            cntj=cntj+1

 if(cnti.eq.cntj) then
   if(abs(oj-oi).gt.1d-7) then
      write(ulog,9)'SUB_COLL:cnti,oi,oj,qi,qj=',cnti,oi,oj,qi,qj
      stop
   endif
 endif
            do lk=1,ndyn
               ok=eigenval(lk,nk)       ;  nbl=nbe(ok,temp,classical)
! 
! Q_kq=0.5*\sum_2 (A_kq;2 + A_2k;q + A_1q;k) diagonal elements (inverse relaxation times*n(n+1))
! M_kq=    \sum_2 (A_kq;2 - A_2k;q - A_1q;k) off diagonal elements 
               call cal_QM(mi,li,j,lj,nk,lk,oi,oj,ok,nbi,nbj,nbl,q_12,m_12)
! in cal_QM, the first index is in IBZ, while the second is in kpc, just like in v33
               if (inside.eq.1) then
!               Qvalue_N(i,li)=Qvalue_N(i,li)+v3sq(i,j,li,lj,lk)* &
!&                     (0.5*delta_l(oi-oj-ok,etaz)*nbj*nbl*(nbi+1) + &
!&                          delta_l(oi-oj+ok,etaz)*nbi*nbl*(nbj+1)   ) 
                  Qvalue_N(mi,li)=Qvalue_N(mi,li)+q_12 *4*pi*pi/nkc  
               else
                  Qvalue_U(mi,li)=Qvalue_U(mi,li)+q_12 *4*pi*pi/nkc
               endif
               if(cnti.ne.cntj) Mval(cntj,cnti)=mval(cntj,cnti)+m_12 *4*pi*pi
            enddo

         enddo
         enddo
         qn=qvalue_n(mi,li) ; qu=qvalue_u(mi,li)
         qval(cnti)=qn+qu  !qvalue_n(i,li)+qvalue_u(i,li)
         n_rate=qn/nbi/(1+nbi); u_rate=qu/nbi/(1+nbi)
         tau_ps=1d+10/c_light/(u_rate+n_rate)
         write(uq,3)i,li,oi,qval(cnti),n_rate,u_rate,n_rate+u_rate,tau_ps,tau_ps*length(velocibz(:,li,i))*c_light*1d-3
!        write(uq+2,3)i,li,oi,qvalue_n(mi,li)+qvalue_u(mi,li),mval(cnti,cnti)
      enddo
      enddo

      do li=1,ndyn
      do i=ksub(1),ksub(2)
         mi=i-ksub(1)+1
         oi=eigenval(li,mapinv(i))              ;  nbi=nbe(oi,temp,classical)
         rate=(qvalue_n(mi,li)+qvalue_u(mi,li))/nbi/(nbi+1)
         tau_ps=1d+10/c_light/rate
         write(uq+3,3)li,i,oi,rate,qvalue_n(mi,li)/nbi/(nbi+1),qvalue_u(mi,li)/nbi/(nbi+1),tau_ps 
      enddo
      enddo

! now write the data into a file
      do j=1,nkc*ndyn
         write(uq+1,4)j, (mval(j,i),i=1,ndyn*(ksub(2)-ksub(1)+1))
      enddo
   endif

 3 format(2i6,9(1x,g10.4))
 4 format(i6,9999(1x,g11.5))
 9 format(a,i6,2(1x,g12.6),2(2x,3(1x,f5.2)))

 close(uq)
 close(uq+1)

 end subroutine calculate_sub_collision
!========================================= 
 subroutine put_collision_together(temp)
! the v33 files for each ksubset are already read. This reads the collision matrix into a
! single array mval(nkc*ndyn by nibz*ndyn) and qval(nibz*ndyn) (diagonal elements) at T=temp
! A_ij;k = 2pi/hbar^2 delta(wk-wi-wj) v3sq(i,j,k) (nk+1)*ni*nj
 use kpoints
 use eigen
 use params
 use exactBTE2
 use phi3
 use om_dos
 implicit none
 integer i,j,jj,li,lj,lk,nk,i1,j1,k1,inside,cnti,cntj
 real(8) qi(3),qj(3),qk(3),temp,q_12,m_12,v33_2
 real(8) nbe, delta_l,nbi,nbj,nbl,oi,oj,ok,qn,qu,n_rate,u_rate,tau_ps

integer uQ
uQ=9402
open(uQ  ,file='Qvalue.dat',status='unknown')
open(uQ+1,file='Mvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Qval,N-rate,U-rate,total(1/cm),tau(ps),la(nm)'


! read the different v3sq files from all ksubsets to generate v3sq(nibz,nkc,ndyn,ndyn,ndyn)
! obviously it assumes a run with all v33 calculations was already done
! call read_all_v3sq

 call allocate_iter0(nibz,nkc,ndyn)  ! allocate mval,qval,qval_U&N
! first the diagonal elements (inverse relaxation times)
! Q_kq=0.5* \sum_2 (A_kq;2 + A_2k;q + A_1q;k) 
   if (usetetra.eq.1) then  ! use tetrahedron method



   else  ! use lorentzian or gaussian broadening
      
! first the diagonal elements Q_value are calculated
      Qvalue_N=0; Qvalue_U=0 ; Mval=0
      cnti=0
      do i=1,nibz
         qi=kibz(:,i)
      do li=1,ndyn
!        oi=eivalibz(li,mapinv(i))   ;  nbi=nbe(oi,temp,classical)
         oi=eivalibz(li,i)   ;  nbi=nbe(oi,temp,classical)
         cnti=cnti+1 

         cntj=0
         do jj=1,nkc
! we introduce j=mapinv(jj) so that j=the first nibz kpc's are in the IBZ 
            j=mapinv(jj)
            qj=kpc(:,j)  
            qk=-qi-qj        
            call get_k_info(qk,NC,nk,i1,j1,k1,g1,g2,g3,inside)
         do lj=1,ndyn
            oj=eigenval(lj,j )          ;  nbj=nbe(oj,temp,classical)
            cntj=cntj+1

            do lk=1,ndyn
               ok=eigenval(lk,nk)          ;  nbl=nbe(ok,temp,classical)
               call cal_QM(i,li,j,lj,nk,lk,oi,oj,ok,nbi,nbj,nbl,q_12,m_12)
               if (inside.eq.1) then
!                      Qvalue_N(i,li)=Qvalue_N(i,li)+v3sq(i,j,li,lj,lk)* &
!&                     (0.5*delta_l(oi-oj-ok,etaz)*nbj*nbl*(nbi+1) + &
!&                          delta_l(oi-oj+ok,etaz)*nbi*nbl*(nbj+1)   ) 
                  Qvalue_N(i,li)=Qvalue_N(i,li)+q_12 *4*pi*pi/nkc
               else
                  Qvalue_U(i,li)=Qvalue_U(i,li)+q_12 *4*pi*pi/nkc
               endif
               if(onedim(li,i).ne.cnti .or.onedim(lj,j).ne.cntj ) then
                  write(ulog,*)'PUT_COLL: 1d(li,i),cnti 1d(lj,j),cntj=',onedim(li,i),cnti,onedim(lj,j),cntj 
               endif
!              Mval(cntj,cnti)=Mval(cntj,cnti)+m_12
               if(cnti.ne.cntj) Mval(cntj,cnti)=mval(cntj,cnti)+m_12 *4*pi*pi
            enddo

         enddo
         enddo
         qn=qvalue_n(i,li) ; qu=qvalue_u(i,li)
         qval(cnti)=qn+qu !qvalue_n(i,li)+qvalue_u(i,li)
         n_rate=qn/nbi/(1+nbi); u_rate=qu/nbi/(1+nbi)
         tau_ps=1d10/c_light/(u_rate+n_rate)
         write(uq,3)i,li,qval(cnti),n_rate,u_rate,n_rate+u_rate,tau_ps,tau_ps*length(velocibz(:,li,i))*c_light*1d-3
      enddo
      enddo

      do j=1,nkc*ndyn
         write(uq+1,4)j, (mval(j,i),i=1,cnti)
      enddo

   endif

 3 format(2i6,99(1x,g10.4))
 4 format(i6,9999(1x,g10.4))
 close(uq) ; close(uq+1)

 end subroutine put_collision_together
!=========================================
 subroutine setup_cg_matrices(ni,nk,mval,qval,rhs,coll,rhsi)
! this subroutine assumes the columns of the collision matrix I(i,j) for i in FBZ
! and j in IBZ are collected (over j subsets), and setsup the minimization matrices
! A and b such that A_ij*phi_j=b_i  with (i,j) in IBZ
! A_ij=sum_l coll_lj*coll*li  ; bi=sum_l rhs_l coll_li
! now A and b can be used in the CG minimization routine
 implicit none
 integer ni,nk,i,j,l
 real(8) coll(ni,ni),rhs(nk),A(ni,ni),b(ni),mval(nk,ni),qval(ni),rhsi(ni)

 A=0; b=0
 do i=1,ni
    b(i)=dot_product(rhs,coll(:,i))
    do j=1,ni
       a(i,j)=dot_product(coll(:,i),coll(:,j))
    enddo
 enddo

 end subroutine setup_cg_matrices
!========================================= 
 subroutine rearrange_into_ibz(ni,nk,mval,qval,b,coll,bi)
! takes mval qval defined on fbz*ibz and projects them on ibz stored into coll
! same for b mapped to bi
 implicit none
 integer ni,nk,i,j,al
 real(8) mval(nk,ni),coll(ni,ni),qval(ni),b(nk,3),bi(ni,3)

 coll=matmul(transpose(mval),mval)
 do i=1,ni
    do al=1,3
       bi(i,al)=dot_product(b(:,al),mval(:,i))
    enddo
!   coll(i,i)=coll(i,i)+qval(i)*qval(i)    ! qval already added to mval
!   do j=1,ni
!      coll(i,j)=coll(i,j)+2*mval(i,j)*qval(i)
!   enddo
 enddo

! do al=1,3
! do j=1,ni
!    bi(j,al)=bi(j,al)+b(j,al)*qval(j)
! enddo
! enddo

 end subroutine rearrange_into_ibz
!========================================= 
 subroutine initialize_phi(temp,nk,ndn,eiv,vel,qval,nedist)
! initializes the distribution function in the IBZ according to the RTA
! qval=n(n+1)/tau; ne_dist*(-grad_alpha T) = phi_RTA; ne_dist=v*tau*hbarom/T^2
! units are in cm^2 as they should
! use kpoints
 use params
! use eigen
! use phi3
 use constants
! use exactbte2
 implicit none
 integer i,l,n,al,nk,ndn
 real(8) ni,nbe,temp,eiv(ndn,nk),vel(3,ndn,nk),qval(ndn*nk),nedist(nk*ndn,3)

 nedist=0 
 do al=1,3
    n=0
    do i=1,nk
    do l=1,ndn
       n=n+1
       ni=nbe(eiv(l,i),temp,classical)
       nedist(n,al)=vel(al,l,i)*c_light/qval(n)*ni*(ni+1)*eiv(l,i)/temp/temp
    enddo
    enddo
 enddo
! ne_dist_rta=ne_dist

 end subroutine initialize_phi
!========================================= 
 subroutine read_all_v3sq
! used for readv3.eq.2 assuming all v3 matrix elements were calculated in split mode and stored
 use phi3
 use eigen
 use kpoints
 use params
 use io2
 implicit none
 integer i,j,l1,l2,l3,ksub(2),n1,nk,n2,lm1,lm2,lm3,nibz_proc,nklast_cpu,nv3_split,cnt
 integer nk1,nk2
 character fn*99,cibz1*4,cibz2*4,cnk*6

 nibz_proc=int(nibz/real(ncpu))+1  ! number of IBZ kpoints per cpu, except for the last one
 nklast_cpu=nibz-(ncpu-1)*nibz_proc ! number of IBZ kpoints for the last cpu < nibz_proc
 write(ulog,*)ncpu,' processors each having ',nibz_proc,' kpoints, except for'
 write(ulog,*)'the last one only getting ',nklast_cpu

 if(iter.ne.0) then  ! use this format for full iterative solution
    allocate(v3sq(nibz,nkc,ndyn,ndyn,ndyn))
 else ! this is used for RTA
    nv3_split=nibz*nkc*ndyn*ndyn*ndyn
    call allocate_v33(nv3_split)
 endif

 cnt=0
 do i=1,ncpu

    ksub(1)=nibz_proc*(i-1)+1 
    if ( i.ne.ncpu ) then
       ksub(2)=nibz_proc*i 
    else
       ksub(2)=nibz
       if (nklast_cpu.ne.ksub(2)-ksub(1)+1) then
          write(ulog,*)'READ_ALL: inconsistency in last ksubs, nklast_cpu=',nklast_cpu 
          write(ulog,*)'READ_ALL: ksub(1),ksub(2)=',ksub
          stop
       endif
    endif

    write(ulog,*)'READING ALL of V3SQ, ksubs=',ksub
    write(cibz1,'(i4.4)') ksub(1)
    write(cibz2,'(i4.4)') ksub(2)
    write(cnk,'(i6.6)') nkc

  if(iter.ne.0) then  ! use this format for full iterative solution

    fn=trim(v3path)//'v35-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
    write(ulog,*)'opening the file ',fn
    open(uv3,file=trim(fn),status='OLD') !,FORM='UNFORMATTED')
    read(uv3,*)nk,n2,lm1,lm2
    if (n2.ne.ndyn) then
       write(ulog,*)'READ_V3: inconsistency in ndn:reading ',n2,' instead of ',ndyn
       stop
    endif
    if (lm1.ne.ksub(1) .or. lm2.ne.ksub(2) ) then
       write(ulog,*)'READ_V3: inconsistency in ksub:reading ',lm1,lm2,' instead of ',ksub
       stop
    endif
    write(ulog,*)'the file ',fn,' was successfully opened'

    v3sq=1d50  ! large numbers would signal using v3sq which was not or improperly read
    loop1: do n1=ksub(1),ksub(2)
    loop2: do n2=1,nkc   ! second argument k2 needs to be only in the IFBZ
    do lm1=1 ,ndyn
    do lm2=1 ,ndyn
    do lm3=1 ,ndyn
       read(uv3,*)nk1,nk2,l1,l2,l3,v3sq(nk1,nk2,l1,l2,l3)
       if (n1.ne.nk1) then
          write(ulog,*)'READ_V3: inconsistency in n1:reading ',nk1,' instead of ',n1
          stop
       endif
    enddo
    enddo
    enddo
    enddo loop2
    enddo loop1
    close(uv3)
 
  else   ! use this format for RTA, where v3 is a single index array

    fn=trim(v3path)//'v33-'//cnk//'-'//cibz1//'-'//cibz2//'.dat'
    write(ulog,*)'opening the file ',fn
    open(uv3,file=trim(fn),status='OLD') !,FORM='UNFORMATTED')
    read(uv3,*)nv3_split,n2,lm1,lm2
    if (n2.ne.ndyn) then
       write(ulog,*)'READ_V3: inconsistency in ndn:reading ',n2,' instead of ',ndyn
       stop
    endif
    if (lm1.ne.ksub(1) .or. lm2.ne.ksub(2) ) then
       write(ulog,*)'READ_V3: inconsistency in ksub:reading ',lm1,lm2,' instead of ',ksub
       stop
    endif
    l3=nibz_proc*nkc*n2*n2*n2
    if (nv3_split.ne.l3 ) then
       write(ulog,*)'READ_V3: inconsistency nv3split:reading ',nv3_split,' instead of ',l3
       stop
    endif
    do j=1 ,nv3_split
      cnt=cnt+1
! this needs to be checked as everytime j starts from 1
      read(uv3,*,end=99) nq1(cnt),nq2(cnt),nq3(cnt),v33sq(cnt)
    enddo
    close(uv3)

  endif
    write(*   ,*)'READ_ALL done for subset # ',i
 enddo

 return
99 write(*,*)'READ_ALL_V3SQ: file end reached at cnt=',cnt
   write(ulog,*)'READ_ALL_V3SQ: file end reached at cnt=',cnt
   stop

 end subroutine read_all_v3sq
!===========================================================================
subroutine cal_Q(nk,ndn)
! calculates the diagonal elements of the collision matrix (1/tau)
use exactBTE2

implicit none

integer i,ii,indx,ndn,j,jj,kk
integer nk

integer uQ
uQ=9402
open(uQ,file='Qvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Q'

Qvalue(:,:)=0
indx=1

write(*,*) 'entering cal_Q'

nkF1_loop: do i=1,nk ! this is usually the IBZ points

    laF1_loop: do ii=1,ndn

           nkF2_loop: do j=1,nk
           laF2_loop: do jj=1,ndn
           laF3_loop: do kk=1,ndn

                Qvalue(i,ii)=Qvalue(i,ii) + P1(i,j,ii,jj,kk) + 0.5*P2(i,j,ii,jj,kk)

           enddo laF3_loop
           enddo laF2_loop
           enddo nkF2_loop

           write(uQ,*) i,ii,Qvalue(i,ii)

      enddo laF1_loop
enddo nkF1_loop



end subroutine cal_Q
!========================================================================================
subroutine cal_QM(ni1,l1,n2,l2,n3,l3,om1,om2,om3,nb1,nb2,nb3,q_12,m_12)
! calculates the collision matrix elements Q and M
! in the split_k calculation ni1 starts from 1 to ksub(2)-ksub(1)+1
 use phi3
 use params      
 use om_dos     ! for etaz
implicit none
real(8), intent(in):: om1,om2,om3,nb1,nb2,nb3
real(8) a_123,a_231,a_312,nbe, delta_l,V33_2
real(8), intent(out):: q_12,m_12 !temp
integer, intent(in):: ni1,l1,n2,l2,n3,l3

!n1=mapinv(ni1)   ! that is the FBZ index
!om1=eivalibz(l1,ni1)
!om2=eigenval(l2,n2)
!om3=eigenval(l3,n3)

!nb1=nbe(om1,temp,classical)
!nb2=nbe(om2,temp,classical)
!nb3=nbe(om3,temp,classical)
!
v33_2=v3sq(ni1,n2,l1,l2,l3)
if(v33_2.lt.v3_threshold*v3_threshold) then
  q_12=0d0; m_12=0d0
else
  a_123=delta_l(om1-om2-om3,etaz)*(nb1+1)*nb2*nb3
  a_231=delta_l(om2-om3-om1,etaz)*(nb2+1)*nb3*nb1
  a_312=delta_l(om3-om1-om2,etaz)*(nb3+1)*nb1*nb2
  q_12=0.5*v33_2*( a_123+a_231+a_312) 
  m_12=    v33_2*(-a_123-a_231+a_312) 
endif

end subroutine cal_QM

!=====================================================
  subroutine calculate_kappa_iter(nk,ndn,wk,temp,ne_dist,eiv,vel,kap) ! kappa(ndyn,3,3)
  use params
  use lattice
  use constants
! use kpoints
  use eigen
  implicit none
  integer nk,ndn,i,j,al,be,la,n !,onedim
  real(8) ni,nbe,om,x,temp
  real(8) ne_dist(ndn*nk,3),eiv(ndn,nk),vel(3,ndn,nk),kap(ndn,3,3),wk(nk)

  kap=0
  do al=1,3
  do be=1,3
     n=0
     do la=1,ndn
     do i=1,nk
         n=onedim(la,i) !n+1        
         x=eiv(la,i)/temp
         ni=nbe(x,1d0,classical)
         kap(la,al,be)=kap(la,al,be)+vel(al,la,i)*ni*(ni+1)*wk(i)*eiv(la,i)*ne_dist(n,be)
     enddo
     enddo
  enddo
  enddo
  kap=kap/volume_r*1d30*k_B/100  

  end subroutine calculate_kappa_iter
!=====================================================
  subroutine self_energy(nks,k,la,omega,temp,nself,uself,iread,usetetra)
  implicit none
  integer iread,la,usetetra,nks
  complex(8) nself,uself
  real(8) temp,omega,k(3)

if (usetetra .eq. 1) then 
  call function_self_tetra(nks,k,la,omega,temp,nself,uself)
else
  if(iread .lt. 2) then
!    write(*,*)' entering self_new '
!    call function_self_new_sq(k,la,omega,temp,nself,uself)
    call function_self_35(nks,k,la,omega,temp,nself,uself)
!    call function_self_new(k,la,omega,temp,nself,uself)
!    call function_self_sc(k,la,omega,temp,nself,uself)
  elseif(iread.eq.2)then  ! on the fly but for k belonging to kpoint-mesh kpc
!    write(*,*)' entering self_w '
!    write(*,*)' args :k,la,om,temp=',k,la,omega,temp
!     call function_self_w2(k,la,omega,temp,nself,uself)
     call function_self_w2_sy(k,la,omega,temp,nself,uself)   ! sy
  elseif(iread.eq.3)then  ! on the fly but for arbitrary kpoint
     call function_self_w3(k,la,omega,temp,nself,uself)
  elseif(iread.eq.4)then  ! only imaginary part using gaussian for arbitrary kpoint
     call function_self_w4(k,la,omega,temp,nself,uself)
  endif
endif

  end subroutine self_energy
!=====================================================
  subroutine write_dist(n,dist,udst)
use eigen
use kpoints
  implicit none
  integer  n,i,j,udst,la,ik,l !,onedim
  real(8) dist(n,3)

  i=0
! do i=1,n
  do ik=1,nibz  !nkc
  do la=1,ndyn
!     i=onedim(la,ik) !
     i=i+1
 !   la=nband(i,ndyn) ; ik=nkpt(i,ndyn)
     write(udst,4)i,(dist(i,j),j=1,3),(sign(dist(i,l),velocibz(l,la,ik))/(1d-13+dist(i,l)),l=1,3)
  enddo
  enddo
 
4 format(i6,2x,3(1x,f8.3),2x,3(1x,f4.1))

  end subroutine write_dist
!=====================================================
  subroutine write_kappa(ndn,kap,unt)
  implicit none
  integer i,j,l,ndn,unt
  real(8) kap(ndn,3,3)

  write(unt,*)'#la,kap_xx(la),kap_yy(la),kap_zz(la),kap_xy(la),kap_xz(la),kap_yz(la)'
  do l=1,ndn
     write(unt,4)l,(kap(l,i,i),i=1,3),((kap(l,i,j),j=i+1,3),i=1,3)
  enddo
  write(unt,4)l,(sum(kap(:,i,i)),i=1,3),((sum(kap(:,i,j)),j=i+1,3),i=1,3)
 
4 format(i6,9(1x,g10.4))
  end subroutine write_kappa
!================================================
 subroutine precondition_all(n,coll,rhs,dist,qval)
! takes as input the rhs,collision matrix starting distribution function and 
! preconditions it by dividing/multiplying by sqrt(n(n+1)/tau)=sqrt(qval)
 implicit none
 integer i,j,n,al
 real(8) coll(n,n),rhs(n,3),qval(n),dist(n,3)

 do al=1,3
    do i=1,n
       dist(i,al)=dist(i,al)*sqrt(qval(i))
       rhs (i,al)=rhs (i,al)/sqrt(qval(i))
    enddo
 enddo
 do i=1,n
    do j=1,n
       coll(i,j)=coll(i,j)/sqrt(qval(i)*qval(j))
    enddo
 enddo

 end subroutine precondition_all
!================================================
 subroutine switch_back(n,dist,qval)
! takes as input the solution to cg minimization and changes back to the old form
! by multiplying by sqrt(n(n+1)/tau)=sqrt(qval)
 implicit none
 integer i,j,n,al
 real(8) qval(n),dist(n,3)

 do al=1,3
    do i=1,n
       dist(i,al)=dist(i,al)/sqrt(qval(i))
    enddo
 enddo

 end subroutine switch_back
!================================================
 subroutine memory_usage
 integer mb,complex_mem,real_mem
 real_mem=8
 complex_mem=16
 mb=1024*1024

!  real_mem*n/mb ! is the memory usage in Mb for a real array of size n

 end subroutine memory_usage
!=====================================================
 subroutine rta_sub(ksub,tempk)
! for subset of kpoints IBZ (first index of v33) in the set ksub, computes the self-energy
! lifetimes, mfps and thermal conductivity of the corresponding modes, 
! writes into a file for later collection
 use kpoints
 use eigen
 use params
 use lattice
 use constants
 use io2
 use exactbte2
 use phi3  ! for readv3
 implicit none 
 integer ksub(2),i,j,l,la,jks,k
 real(8), allocatable :: evc(:),omega(:)
 complex(8), allocatable :: nself(:),uself(:)
 real(8) temp,tempk,mfp,kappa_q(3,3),trace
 character ctemp*3,cnk*6,cibz1*4,cibz2*4

 allocate(evc(ksub(1):ksub(2)),omega(ndyn),nself(ndyn),uself(ndyn),kappa_rta(ndyn,3,3))

 temp = tempk*k_b/(100*h_plank*c_light) ! to convert from SI to cm^-1

 write(cibz1,'(i4.4)') ksub(1)
 write(cibz2,'(i4.4)') ksub(2)
 write(cnk,'(i6.6)') nkc
 write(ctemp,'(i3.3)')floor(tempk)

 open(ukap,file= 'kappa-'//ctemp//'-'//cnk//'-'//cibz1//'-'//cibz2//'.dat',status='unknown')
 open(umod,file='modkap-'//ctemp//'-'//cnk//'-'//cibz1//'-'//cibz2//'.dat',status='unknown')
 open(uslf+1  ,file='rself-'//ctemp//'-'//cnk//'-'//cibz1//'-'//cibz2//'.temp')
 open(uslf+500,file='iself-'//ctemp//'-'//cnk//'-'//cibz1//'-'//cibz2//'.temp')
 write(ukap,*)'#tempk,(kappa(la),la=1,ndyn),sum(kappa)'
 write(umod,*)'#la,tempk,(kappa_rta(la,i,i),i=1,3),((kappa_rta(la,i,j),j=i+1,3),i=1,2)'

 kappa_rta=0
 modeloop: do la=1,ndyn
    write(uslf+1  ,2)'#la,nk,tempk,omega(k,la),real(self(la)),normal,umklapp '
    write(uslf+500,2)'#la,nk,tempk,omega(k,la),2*imag(self(la)),2*normal,2*umklapp,tau(ps),MFP(nm),Tr(kap)/3,kap123-654 '
! SELF-ENEGY calculation on the ibz sub-mesh
    do j=ksub(1),ksub(2)
       jks=j-ksub(1)+1
       omega(la) = eivalibz(la,j)
       call self_energy(jks,kibz(:,j),la,omega(la),temp,nself(la),uself(la),readv3,usetetra)
       evc(j)=2*aimag(nself(la)+uself(la))   ! used for inverse relaxation time
       mfp=length(velocibz(:,la,j)) / evc(j) * 1d7   ! in nanometers

       do k=1,nkc  ! this loop selects only the stars of j (therefore no shifts allowed)
          if(j.ne.mapibz(k)) cycle 

! contribution to the thermal conductivity including the weight associated with kibz
!      call mode_thermal_conductivity_single_k(wibz(j),evc(j),velocibz(:,la,j),eivalibz(la,j),temp,kappa_q)
          call mode_thermal_conductivity_single_k(wk(k),evc(j),veloc(:,la,k),eivalibz(la,j),temp,kappa_q)

          if(verbose) write(ulog,3)la,j,sum(kappa_q,MASK= i==l),kappa_q

! Thermal conductivity from the RTA to the solution of Boltzmann equation

          kappa_rta(la,:,:) = kappa_rta(la,:,:) + kappa_q(:,:)    ! here used as a dummy array

       enddo

       write(uslf+1  ,8)la,j,tempk,omega(la),real (nself(la)+uself(la)), &
&           real (nself(la)),real (uself(la))
       write(uslf+500,8)la,j,tempk,omega(la),2*aimag(nself(la)+uself(la)), &
&           2*aimag(nself(la)),2*aimag(uself(la)),1d10/c_light/evc(j),mfp,   &
&           trace(3,kappa_rta(la,:,:))/3,(kappa_rta(la,i,i),i=1,3),((kappa_rta(la,i,l),l=i+1,3),i=1,3)

       write(ulog,4)la,sum(kappa_rta(la,:,:),MASK= i==l),kappa_rta(la,:,:)
       write(*   ,*) 'self_energy done; la,jks=',la,jks

    enddo
    write(uslf+1  ,*)' '
    write(uslf+500,*)' '

    omega(la)=(kappa_rta(la,1,1)+kappa_rta(la,2,2)+kappa_rta(la,3,3))/3
    write(umod,4)la,tempk,(kappa_rta(la,i,i),i=1,3),((kappa_rta(la,i,j),j=i+1,3),i=1,2)

 enddo modeloop

 write(ukap,6)tempk,(omega(la),la=1,ndyn),sum(omega)

 deallocate(evc,omega,nself,uself,kappa_rta)
 close(uslf+1  )
 close(uslf+500)
 close(ukap)
 close(umod)

2 format(a,2x,f19.12,9(2x,g11.5))
3 format(i3,1x,i6,1x,g11.5,3x,99(1x,g9.3))
4 format(i3,1x,g11.5,3x,99(1x,g9.3))
6 format(99(2x,g10.4))
8 format(2i6,2x,f7.2,2x,f8.3,3(1x,g11.4),2(2x,g9.3),3x,f8.3,2x,9(1x,f7.2))

 end subroutine rta_sub
!=====================================================
!function onedim(nb,nk)
! calculates the running index onedim
! calculates the running index onedim
!use eigen
!implicit none
!integer nk,nb,onedim
!! onedim=(nb-1)*nkc2+nk !if outermost loop is j=1,ndyn and innermost loop is i=1,nkc2
!onedim=(nk-1)*ndyn+nb   !if innermost loop is j=1,ndyn and outermost loop is i=1,nkc2
!end function onedim
!=====================================================
!function nband(onedm)
!use eigen
!implicit none
!integer nband,onedm
!nband = int((onedm-1)/nkc2)+1
!end function nband
!=====================================================
!function nkpt(onedm)
!use eigen
!implicit none
!integer nkpt,nband,onedm
!nkpt = mod(onedm-1,nkc2)+1
! nband = int(onedm-nkpt)/nkc2+1
!end function nkpt

