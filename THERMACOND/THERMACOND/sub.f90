!==========================================================
 subroutine jdos_bs_sy
! Calculate iself except for matrix element contribution along specified sym. lines
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

 integer i,la
 real(8) q(3),omq(ndyn)
 real(8) temperature, jdos1p(ndyn),jdos1n(ndyn),jdos2p(ndyn),jdos2n(ndyn),jdos_tot(ndyn) 

 open(udos,file='lifetime_dos_bs_occ.dat ',status='unknown')
 write(udos,*) 'dk_bs, q, omq, jdos1p, jdos1n, jdos2p, jdos2n, jdos_tot(la)'
 temperature=200d0

 do i=1,nkp_bs
       q=kp_bs(:,i)
       do la=1,ndyn
          omq(la) = sqrt(abs(eigenval_bs(la,i))) * cnst
       enddo
       call lifetime_dos_sy(temperature,q,omq,jdos1p,jdos1n,jdos2p,jdos2n,jdos_tot)
       write(udos,3) dk_bs(i),q,(omq(la),jdos1p(la),jdos1n(la),jdos2p(la),jdos2n(la),jdos_tot(la),la=1,ndyn)
       write(ulog,*) 'i, nkp_bs', i, nkp_bs
 enddo
 close(udos)

3 format(99(1x,g11.4))

 end subroutine jdos_bs_sy
!============================================================
 subroutine jdos_sy
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

 integer i,la
 real(8) q(3),omq(ndyn)
 real(8) temperature, jdos1p(ndyn),jdos1n(ndyn),jdos2p(ndyn),jdos2n(ndyn),jdos_tot(ndyn)

 open(udos,file='lifetime_dos_occ.dat ',status='unknown')
 write(udos,*) 'i,q, omq, jdos1p, jdos1n, jdos2p, jdos2n, jdos_tot(la)'
 temperature=200d0

 do i=1,nkc
       q=kpc(:,i)
       do la=1,ndyn
          omq(la) = sqrt(abs(eigenval_bs(la,i))) * cnst
       enddo
       call lifetime_dos_sy(temperature,q,omq,jdos1p,jdos1n,jdos2p,jdos2n,jdos_tot)
       write(udos,3) i,q,(omq(la),jdos1p(la),jdos1n(la),jdos2p(la),jdos2n(la),jdos_tot(la),la=1,ndyn)
       write(ulog,*) 'i, nkc', i, nkc
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
 real(8) om1,om2,omq(ndyn),delta_l,evl2(ndyn)
 complex(8) oc2,omz,evc2(ndyn,ndyn)

 real(8) temperature, n_be1, n_be2


 jdos1p=0; jdos1n=0; jdos2p=0; jdos2n=0; jdos_tot=0 
 do l=1,nkc
    q1=kpc(:,l)
    call get_k_info(q1,NC,nq1,i1,j1,k1,g1,g2,g3,inside)
    if(l.ne.nq1) then
      write(ulog,*)'n1,nq1,inside=',l,nq1,inside
      write(ulog,*)'q1=',q1
      write(ulog,*)'ijk1=',i1,j1,k1
      stop
    endif
    q2=-q-q1          ! enforce momentum conservation
    call get_freq(q2,ndyn,evl2,evc2)

    do la1=1,ndyn
    do la2=1,ndyn
    do la =1,ndyn

       om1 = sqrt(abs(eigenval(la1,l))) * cnst
       om2 = sqrt(abs(evl2(la2))) * cnst
       
       call BS_distribution_sy(temperature,om1,n_be1)
       call BS_distribution_sy(temperature,om2,n_be2)
 
       !write(debug,8) 'temperature,om1,n_be1,om2,n_be2', temperature, om1, n_be1, om2, n_be2
       !write(*,8) 'temperature,om1,n_be1,om2,n_be2', temperature, om1, n_be1, om2, n_be2      
 
       jdos1p(la)=jdos1p(la)+(1.0d0+n_be1+n_be2)*delta_l(om1+om2+omq(la),etaz)
       jdos1n(la)=jdos1n(la)+(1.0d0+n_be1+n_be2)*delta_l(om1+om2-omq(la),etaz)
       jdos2p(la)=jdos2p(la)+(n_be2-n_be1)*delta_l(om1-om2+omq(la),etaz)
       jdos2n(la)=jdos2n(la)+(n_be2-n_be1)*delta_l(om1-om2-omq(la),etaz)
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

 n_be=1d0/(EXP(hbar*omega_Hz/(k_b*T))-1)
 
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
    omq2(l2) = sqrt(abs(eigenval(l2,nqq2))) * cnst
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
 subroutine calculate_w3_ibz_split(ibz_subset,nv3_split,ndn,ni,ki,nk,kp,eival,eivec)
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
 integer nk,ni,ndn,l1,l2,l3,n2,n3,j,k
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,indx
 integer ibz_subset(2),nv3_split
 real(8) q1(3),q2(3),q3(3),eival(ndn,nk)
 real(8) ki(3,ni),kp(3,nk)
 complex(8) eivec(ndn,ndn,nk),xx
 real cputim
 
 v33 = cmplx(0d0,0d0)
 indx=0
 write(ulog,*)'entered v3 subroutine...'
 write(ulog,*)'V3: nc(123)=',nc
 write(ulog,*)'starting ibz point=',ibz_subset(1)
 write(ulog,*)'ending ibz point=', ibz_subset(2)

!************************************************
! second argument q2 needs to be only in the IFBZ
!************************************************
 loop2: do n2=ibz_subset(1),ibz_subset(2)       ! calculate v33 only for subset of k points
    q2 = ki(:,n2)  ! should be = kp(:,mapinv(n2))
    !write(debug,*) 'q2=' q2
    call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(mapinv(n2).ne.nk2) then
      write(ulog,*)'n2,mapinv(n2),nk2,inside=',n2,mapinv(n2),nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif

 call cpu_time(cputim)  !date_and_time(time=tim)
 write(ulog,*) 'entering q2...',n2,q2
 write(ulog,'(a,f12.4)')'  TIME IS ',cputim
   

 loop3: do n3=1 ,nk
    q3 = kp(:,n3)   ! third argument is in the whole FBZ coarse mesh
    !write(debug,*) 'q3=', q3
    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
    if(n3.ne.nk3) then
      write(ulog,*)'n3,nk3,inside=',n3,nk3,inside
      write(ulog,*)'q3=',q3
      write(ulog,*)'ijk3=',i3,j3,k3
      stop
    endif

    q1 = -q2-q3

    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)

 ! write(ulog,2)'nk123=',nk1,nk2,nk3
  !write(ulog,*)'q2,q3,q1',q2,q3,q1
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
    !write(debug,*),'l1,l3,l2,indx,v33',l1,l3,l2,indx,v33(indx)
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

 if( writev3.eq.1) then
   open(uv3,file='v33.dat',status='unknown',FORM='UNFORMATTED')
   !open(uv3,file='v33.dat',status='unknown')
   write(uv3)nv3_split
   !write(uv3,*) nv3_split
   do j=1 ,nv3_split
      write(uv3) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
   !   write(uv3,7) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
   enddo
   close(uv3)
 endif

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))

 end subroutine calculate_w3_ibz_split


!-------------------------------------------------------------------
 subroutine read_ksubset (ksubset)
! this subroutine reads ksubset.inp
 use io2
 implicit none
 integer ksubset_f
 integer ksubset(2)
  
 ksubset_f=9998
 ksubset=0
 open(ksubset_f,file='ksubset.inp', status='old')
 read(ksubset_f,*) ksubset
 close(ksubset_f)
! write(debug,*) 'ksubset=', ksubset(1), ksubset(2)
 end subroutine read_ksubset


!-------------------------------------------------------------------
 subroutine read_merge (mergemode, num_files)
! this subroutine reads ksubset.inp

 use io2

 implicit none
 integer i,merge_f
 integer mergemode, num_files

 merge_f=9998
 mergemode=0
 open(merge_f,file='merge.inp', status='old')
 read(merge_f,*) mergemode
 read(merge_f,*) num_files
 close(merge_f)
! write(debug,*) 'mergemode, num_files', mergemode, num_files
 end subroutine read_merge


!-----------------------------------------------------------------
 subroutine merge_v3(num_files)
 
 use phi3
 use io2
 use lattice
 implicit none
 integer i,num_files
 integer j,unt
 real(8) xx,yy
 character(99) filename
 integer nv3_split, nv3_split_acc

 write(*,*) 'Now in merge_v3'
 write(*,*) 'nv3=',nv3

 unt = 111
 open(unt+1,file='v33.dat',status='unknown',FORM='UNFORMATTED')
 ! open(unt+1,file='v33.dat',status='unknown')
 !write(unt+1,*) nv3
 write(unt+1) nv3 

 nv3_split_acc=0
 do i=1,num_files 
    write(filename,fmt="(a,i3.3,a)") "v33_",i,".dat"
    open(unt,file=filename,status='old',FORM='UNFORMATTED')
    !open(unt,file=filename,status='old')
    !read(unt,*) nv3_split
    read(unt) nv3_split

    nv3_split_acc=nv3_split_acc+nv3_split    
    
    write(*,*) 'i,nv3_split',i,nv3_split

    call allocate_v33(nv3_split)
    write(ulog,*)' OPENING v33.dat, reading ',nv3_split
    v33= cmplx(0d0,0d0)
    do j=1,nv3_split
      read(unt) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),xx,yy
      !read(unt,*) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),xx,yy
      !write(unt+1,7) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),xx,yy       
      v33(j)=cmplx(xx,yy)
      write(unt+1) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33(j)
    enddo
    close(unt)
    call deallocate_v33
 enddo
 close(unt+1)
7 format(6(i6),9(2x,g14.8))

 !enddo
 !close(unt)
 !return

!99 write(*,*)'READ_V3_NEW: v33 file end reached at line j=',j
!   stop

 if (nv3_split_acc .ne. nv3) then
   write(ulog,*) 'WARNING: nv3_split_acc .ne. nv3'
 endif

 end subroutine merge_v3


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
 subroutine function_self_w2_sy(q_ibz,q,la,omega,temp,nself,uself)
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
 real(8) omq,om3,omk,eta,tk,v32,term,k1(3),k2(3),q1(3),q3(3)
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
 omq = sqrt(abs(eigenval(la,nq))) * cnst
 eta= etaz
! write(ulog,5)'SELF_W: last cntr, omega,etaz=',cntr,omega,etaz
 omz = cmplx(omega,-eta)   ! this is in cm^-1
 nself=0 ; uself=0
 v32 = v3_threshold*v3_threshold


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
 tempk=temp*(100*h_plank*c_light)/k_b
 do k_FBZ=1, nkc
   write(*,*) 'k_FBZ, mapibz(k_FBZ), nq', k_FBZ, mapibz(k_FBZ), nq
   if(mapibz(k_FBZ) .eq. q_ibz) then   ! IBZ -> FBZ
     write(*,*) 'matched--------------------------------------------------------'
     write(self_detail,8) tempk, k_FBZ, kpc(:,k_FBZ), la, omega, &  ! ki, kpoint, lambda, omega
& v33_square_sum, real(delta1_sum), aimag(delta1_sum), real(delta2_sum), aimag(delta2_sum) &
& , real(delta3_sum), aimag(delta3_sum), real(delta4_sum), aimag(delta4_sum), real(delta_tot_sum), aimag(delta_tot_sum)
   endif
 enddo

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
6 format(i6,6(1x,f9.4),2x,3(1x,i2))
7 format(a,4i3,99(1x,f7.3))

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
 real(8) omq,eta,term,k3(3),k2(3),eivq(ndyn),eiv2(ndyn),eiv3(ndyn)
 real(8) etacut,arg,delta_l,om3,om2,nb3,nb2,nbe,ires,rres
 complex(8) omz,self,xx,evq(ndyn,ndyn),ev2(ndyn,ndyn),ev3(ndyn,ndyn)

 real(8) v33_square_sum, tempk
 real(8) delta1_sum(2), delta2_sum(2), delta3_sum(2), delta4_sum(2), delta_tot(2)

 tempk = temp/k_b*(100*h_plank*c_light)


! be careful band sorting should be avoided in this case; otherwise need to
! find the nearest kpoint in kpc, and sort eivq according to it
write(*,*) 'la, nq, q, omega, tempk',la,nq,q,omega,tempk

 call get_freq(q,ndyn,eivq,evq)
 omq=sqrt(abs(eivq(la)))*cnst

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

    call get_freq(k3,ndyn,eiv3,ev3)

    do l2=1,ndyn
       om2=sqrt(abs(eigenval(l2,ik)))*cnst
       nb2=nbe(om2,temp,classical)
    do l3=1,ndyn

! first check see if the energy is conserved
       om3=sqrt(abs(eiv3(l3)))*cnst
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

 om3 = sqrt(abs(eigenval(l3,nk3))) * cnst ; nb3=nbe(om3,temp,classical)
 om2 = sqrt(abs(eigenval(l2,nk2))) * cnst ; nb2=nbe(om2,temp,classical)
! om=real(omz) ; et=-aimag(omz)

 delta1=(nb2+nb3+1)*(1d0/(om2+om3+omz))
 delta2=(nb2+nb3+1)*(1d0/(om2+om3-omz))
 delta3=(nb2-nb3  )*(1d0/(om3-om2+omz))
 delta4=(nb2-nb3  )*(1d0/(om3-om2-omz))
 delta_tot = delta1+delta2+delta3+delta4


 end subroutine ocfunc_sy
!===========================================================



!=====================================================================
! read splitted collision matrix and merge them into n files.
! n should be specified in mergeP.inp file
! each merged file shouldn't exceed limit in memory
subroutine merge_Pmatrix

implicit none

integer Pmerge, num_pfiles
integer, allocatable :: merge_file(:,:)
integer i,j,k
character(99) filename1,filename2
integer ncol
integer, allocatable :: n1(:),n2(:),l1(:),l2(:),l3(:)
real(8), allocatable :: col1(:),col2(:)

Pmerge=10000

open(Pmerge,file='Pmerge.inp',status='unknown')
read(Pmerge,*) num_pfiles
allocate(merge_file(num_pfiles,2))

do i=1,num_pfiles
   read(Pmerge,*) merge_file(i,1), merge_file(i,2)
enddo
close(Pmerge)

do i=1,num_pfiles

   write(*,*) 'merge P files',i,'/',num_pfiles

   write(filename1,fmt="(a,i3.3,a)") "Pmatrix.",i,".dat"
   open(Pmerge+1,file=filename1,status='unknown')
!   open(Pmerge+1,file=filename1,status='unknown',form='unformatted')


   do j=merge_file(i,1),merge_file(i,2)
 
        write(filename2,fmt="(a,i3.3,a)") "col_matrix.",j,".dat"
        open(Pmerge+2,file=filename2,status='old')
        read(Pmerge+2,*) ncol

!        open(Pmerge+2,file=filename2,status='old',form='unformatted')
!        read(Pmerge+2) ncol

        allocate(n1(ncol),n2(ncol),l1(ncol),l2(ncol),l3(ncol),col1(ncol),col2(ncol))

        do k=1,ncol
             read(Pmerge+2,*) n1(k),n2(k),l1(k),l2(k),l3(k),col1(k),col2(k)
!             read(Pmerge+2) n1(k),n2(k),l1(k),l2(k),l3(k),col1(k),col2(k)
!             read(Pmerge+2) col1(k),col2(k)
        enddo

        do k=1,ncol
             write(Pmerge+1,8) n1(k),n2(k),l1(k),l2(k),l3(k),col1(k),col2(k)
!             write(Pmerge+1) n1(k),n2(k),l1(k),l2(k),l3(k),col1(k),col2(k)
!             write(Pmerge+1) col1(k),col2(k)
        enddo
        
        close(Pmerge+2)
        deallocate(n1,n2,l1,l2,l3,col1,col2)

   enddo

   close(Pmerge+1)

enddo

deallocate(merge_file)

8 format (99(2x,g12.6))

! do i=1,num_pfiles
! make file name for merged P file using i
! open merged p file

! do j=merge_file(i,1), merge_file(i,2)
! make file name using merge_file
! open file
! read
! write
! close
! enddo

! close merged p file
! enddo




end subroutine merge_Pmatrix

!=====================================================================


!=====================================================================
! read v33_split
! calculate collision matrix
! save it to files
subroutine Pmatrix (ksubset,nk,kp,ndn,eigenval,temp)

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

 implicit none

 integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
 integer nk1,nk2,nk3, nk1n, nk2n, nk3n, nk3_proc1, nk3_proc2
 integer ksubset(2),nk
 real(8) q1(3), q2(3), q3(3), q1n(3), q2n(3), q3n(3), q3_proc1(3), q3_proc2(3)
 real(8) kp(3,nk)

 integer i,n1,n2,n3,l1,l2,l3,ndn
 real(8) temp     ! temperature in cm^-1
 real(8) eigenval(ndn,nk)
 real(8) omega1, omega2, omega3, V3sq1, V3sq2, col1, col2, omega3_proc1, omega3_proc2 
 real(8) proc1(nk,nk,ndn,ndn,ndn), proc2(nk,nk,ndn,ndn,ndn)


 integer ucol
 ucol=9010

 write(*,*) 'entering cal Pmatrix...'

! open(ucol,file='col_matrix.dat',STATUS='unknown',FORM='unformatted')
 open(ucol,file='col_matrix.dat',STATUS='unknown')

ncol=(ksubset(2)-ksubset(1)+1)*nk*ndn**3 
write(ucol,*) ncol 

loop1 : do n1=ksubset(1),ksubset(2)

 
    q1=kp(:,n1)
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
    if (n1 .ne. nk1) then
       write(ulog,*) 'n3,nk3,inside=',n3,nk3,inside
       write(ulog,*) 'q3=',q3
       write(ulog,*) 'ijk3=',i3,j3,k3
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



           do l2=1, ndn
           do l3=1, ndn

                 omega1=sqrt(eigenval(l1,nk1)) * cnst    !linear freq.
                 omega2=sqrt(eigenval(l2,nk2)) * cnst 
                 omega3=sqrt(eigenval(l3,nk3)) * cnst 

                 q3_proc1=q1+q2
                 q3_proc2=q1-q2
                 call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside)
                 call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside)
                 
                 omega3_proc1=sqrt(eigenval(l3,nk3_proc1))*cnst
                 omega3_proc2=sqrt(eigenval(l3,nk3_proc2))*cnst

!                 call find_V3(nk1n,nk2,nk3,l1,l2,l3,V3sq2)    ! return V3sq2

                 V3sq1=v33sq1(nk1,nk2,l1,l2,l3) * (2*pi)**2 / nk 
!V33_md(nk1,nk2,l1,l2,l3)=V33_md(nk1,nk2,nk3n,l1,l2,l3) since q3=-q1-q2 when calculating V33_md
                 V3sq2=v33sq1(nk1,nk2,l1,l2,l3) * (2*pi)**2 / nk

                 call cal_P(temp,omega1,omega2,omega3_proc1,omega3_proc2,V3sq1,col1,col2)   ! return P1, P2

                 proc1(nk1,nk2,l1,l2,l3)=col1
                 proc2(nk1,nk2,l1,l2,l3)=col2


            enddo
            enddo
!            enddo
     enddo loop2
write(ulog,*) 'nk1,l1', nk1, l1

enddo
enddo loop1

do n1=ksubset(1),ksubset(2)
do n2=1,nk
do l1=1,ndn
do l2=1,ndn
do l3=1,ndn
                 write(ucol,7) n1, n2, l1, l2, l3, proc1(n1,n2,l1,l2,l3), proc2(n1,n2,l1,l2,l3)
                 ! write(ucol) nk1, nk2, nk3, l1, l2, l3, col1, col2
enddo
enddo
enddo
enddo
enddo


close(ucol)

7 format(5(i6),2(2x,g14.8))
8 format(i6,3(2x,g14.8))

! loop for q3: ksubset(1) to ksubset(2)
! loop for q2 in entire FBZ
!  calculate q1
!  get omega1, omega2, omega3

!  get nk1 for -q1, nk2 for -q2, nk3 for q3
!  find V3(nk1,nk2,nk3,la1,la2,la3)
!  calculate P1

!  get nk1 for -q1, nk2 for q2, nk3 for q3
!  find V3(nk1,nk2,nk3,la1,la2,la3)
!  calculate P2

! end all loops
! Write P1 and P2 

end subroutine Pmatrix



!===========================================================
subroutine cal_P(temp,omega1,omega2,omega3_proc1,omega3_proc2,V3sq1,col1,col2)

! resulting V3sq1, V3sq2 is in cm^-1

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
 use om_dos
! use coll_matrix
! use exactBTE2

implicit none

real(8) temp,omega1,omega2,omega3,V3sq1,V3sq2,col1,col2, omega3_proc1, omega3_proc2
real(8) nb1,nb2,nb3_proc1, nb3_proc2
real(8) nbe, delta_l

nb1=nbe(omega1,temp,classical)
nb2=nbe(omega2,temp,classical)
nb3_proc1=nbe(omega3_proc1,temp,classical)
nb3_proc2=nbe(omega3_proc2,temp,classical)


col1=nb1*nb2*(nb3_proc1+1) * V3sq1 * delta_l(omega1+omega2-omega3_proc1,etaz)  ! no 2pi at the front since delta is in linear frequency.
col2=nb1*(nb2+1)*(nb3_proc1+1) * V3sq1 * delta_l(omega1-omega2-omega3_proc1,etaz)

!write(*,*) 'om in cal_p', omega1,omega2,omega3_proc1,omega3_proc2
!write(*,*) 'nb in cal_p', nb1,nb2,nb3_proc1,nb3_proc2

end subroutine cal_P



!===========================================================
subroutine cal_P2(n1,l1,n2,l2,n3,l3,V3sq1,V3sq2,col1,col2)

! resulting V3sq1, V3sq2 is in cm^-1

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
 use om_dos
 use exactBTE2
! use coll_matrix

implicit none

real(8) temp,omega1,omega2,omega3,V3sq1,V3sq2,col1,col2
real(8) nb1,nb2, nb3
real(8) nbe, delta_l
integer n1,l1,n2,l2,n3,l3


nb1=dist(n1,l1)
nb2=dist(n2,l2)
nb3=dist(n3,l3)

omega1=frequency(n1,l1)
omega2=frequency(n2,l2)
omega3=frequency(n3,l3)

col1=nb1*nb2*(nb3+1) * V3sq1 * delta_l(omega1+omega2-omega3,etaz)  ! no 2pi at the front since delta is in linear frequency.
col2=nb1*(nb2+1)*(nb3+1) * V3sq2 * delta_l(omega1-omega2-omega3,etaz)


end subroutine cal_P2


!===========================================================
subroutine cal_P3(n1,l1,n2,l2,n3_proc1,n3_proc2,l3,V3sq1,V3sq2,col1,col2)

! resulting V3sq1, V3sq2 is in cm^-1

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
 use om_dos
 use exactBTE2
! use coll_matrix

implicit none

real(8) temp,omega1,omega2,omega3,V3sq1,V3sq2,col1,col2
real(8) nb1,nb2, nb3, nb3_proc1, nb3_proc2, omega3_proc1, omega3_proc2
real(8) nbe, delta_l
integer n1,l1,n2,l2,n3,l3, n3_proc1, n3_proc2


nb1=dist(n1,l1)
nb2=dist(n2,l2)
nb3=dist(n3,l3)
nb3_proc1=dist(n3_proc1,l3)
nb3_proc2=dist(n3_proc2,l3)

omega1=frequency(n1,l1)
omega2=frequency(n2,l2)
omega3=frequency(n3,l3)
omega3_proc1=frequency(n3_proc1,l3)
omega3_proc2=frequency(n3_proc2,l3)

col1=nb1*nb2*(nb3_proc1+1) * V3sq1 * delta_l(omega1+omega2-omega3_proc1,etaz)  ! no 2pi at the front since delta is in linear frequency.
col2=nb1*(nb2+1)*(nb3_proc2+1) * V3sq2 * delta_l(omega1-omega2-omega3_proc2,etaz)


end subroutine cal_P3




!===========================================================
subroutine find_V3 (n1,n2,n3,l1,l2,l3,V3sq1)

 use kpoints
 use io2
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use force_constants_module
! use svd_stuff
 use params
 use constants
 use om_dos
! use coll_matrix

implicit none

integer n1,n2,n3,l1,l2,l3
real(8) V3sq1

integer i,match

match=0
i=0

find : do

    i=i+1
    if (n1.eq.nq1(i) .and. n2.eq.nq2(i) .and. n3.eq.nq3(i) .and. l1.eq.la1(i) .and. l2.eq.la2(i) .and. l3.eq.la3(i) ) then
        V3sq1 = V33(i)*conjg(V33(i)) *  (2*pi)**2   ! V3 in Jiv's thesis = V3 in Keivan's code * 2pi 
        exit find
    endif

    if (i .gt. nv3) then
        exit find
    endif

enddo find
! den_proc1 = sqrt(8*sqrt(abs(eigenval(l1,nk1)*eigenval(l2,nk2)*eigenval(l3,nk3_proc1))))
! if (den.ne.0) then
!    w33_proc1 = xx_proc1 / den_proc1 * const33
! else
!    write(*,*)'MATRIX_ELT: den=0 ',den
!    stop
! endif

end subroutine find_V3





!===========================================================
 subroutine calculate_w3_fbz_split_sq(fbz_subset,nv3_split,ndn,nk,kp,eival,eivec)
! v33 (q1,q2,q3,la1,la2,la3), ksubset is in q1, q2 is all kpoints in FBZ, q3 is from momentum conservation 
! Full BZ calculated for iterative solution

! This is a modified calculate_w3_ibz by sy.
! It reads fbz_subset (fbz_subset(1)-starting index, fbz_subset(2)-end index in nibz)
! This subroutine will write header in output file. 

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
 use phi3_sy
 implicit none
 integer nk,ndn,l1,l2,l3,n2,n3,j,k,n1
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,indx
 integer fbz_subset(2),nv3_split
 real(8) q1(3),q2(3),q3(3),eival(ndn,nk)
 real(8) kp(3,nk)
 complex(8) eivec(ndn,ndn,nk),xx
 real cputim, cputim2, cputim3, cputim_matrix_elt
 real(8) v33sq1_temp

 integer utemp
 utemp=100


 call cal_eiqr2(nk,ndn)    ! calculate eiqr2, eivec2

 write(*,*)'exited cal_eiqr2...'

 v33 = cmplx(0d0,0d0)
 indx=0
 write(ulog,*)'entered v3 subroutine...'
 write(ulog,*)'V3: nc(123)=',nc
 write(ulog,*)'starting fbz point=',fbz_subset(1)
 write(ulog,*)'ending fbz point=', fbz_subset(2)

   open(uv3,file='v33sq.dat',status='unknown',FORM='UNFORMATTED')
!   open(uv3,file='v33sq.dat',status='unknown')

   nv3=(fbz_subset(2)-fbz_subset(1)+1)*nk*ndn**3
   write(uv3) nv3
!   write(uv3,*) nv3


open(utemp,file='timeanal.dat',status='unknown')
call cpu_time(cputim2)
write(utemp,*) 'before loop1',cputim2

!************************************************
! second argument q2 needs to be only in the IFBZ
!************************************************
 loop1: do n1=fbz_subset(1),fbz_subset(2)       ! calculate v33 only for subset of k points
    q1 = kp(:,n1)  ! should be = kp(:,mapinv(n2))
    write(*,*) 'in V3 cal, n1=',n1
    write(ulog,*) 'n1=', n1
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
!    if(mapinv(n2).ne.nk2) then
      if(n1.ne.nk1) then
      write(ulog,*)'n1,nk1,inside=',n1,nk1,inside
      write(ulog,*)'q1=',q1
      write(ulog,*)'ijk1=',i1,j1,k1
      stop
    endif

 call cpu_time(cputim)  !date_and_time(time=tim)
 write(ulog,*) 'entering q1...',n1,q1
 write(ulog,'(a,f12.4)')'  TIME IS ',cputim
 write(utemp,*) 'entering q1...',n1,q1
 write(utemp,'(a,f12.4)')'  TIME IS ',cputim

 cputim_matrix_elt=0.0


do l1=1,ndn

 loop2: do n2=1 ,nk
    q2 = kp(:,n2)   ! third argument is in the whole FBZ coarse mesh
    !write(debug,*) 'q3=', q3
    call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(n2.ne.nk2) then
      write(ulog,*)'n2,nk2,inside=',n2,nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif

    q3 = -q2-q1

   call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)


 ! write(ulog,2)'nk123=',nk1,nk2,nk3
  !write(ulog,*)'q2,q3,q1',q2,q3,q1
 do l2=1 ,ndn
 do l3=1 ,ndn

call cpu_time(cputim2)

!!! CAUTION: use matrix_elt_simplified only when basis atoms have the same mass
!!! Otherwise, use matrix_elt_simplified2
call matrix_elt_simplified2(q1,q2,l1,l2,l3,xx,inside,ndn)    ! same as matrix_elt but efficiency is improved.
!call matrix_elt_simplified2((q1,q2,l1,l2,l3,xx,inside,ndn)

call cpu_time(cputim3)
cputim_matrix_elt=cputim_matrix_elt + (cputim3-cputim2)

v33sq1_temp=xx*conjg(xx)

write(uv3) v33sq1_temp
!write(uv3,9) nk1,nk2,nk3,l1,l2,l3,v33sq1_temp

    indx=indx+1
!    v33(indx)=xx
!    nq1(indx)= nk1
!    nq2(indx)= nk2
!    nq3(indx)= nk3
!    la1(indx)= l1
!    la2(indx)= l2
!    la3(indx)= l3
    !write(debug,*),'l1,l3,l2,indx,v33',l1,l3,l2,indx,v33(indx)
!   if (mod(indx,1000).eq.1) write(*,*) 'indx,v33=',indx,xx
! enddo
 enddo
 enddo
!    write(ulog,5)n2,n3,l1-1,l2-1,l3-1,v3(n2,n3,l1-1,l2-1,l3-1)
 enddo loop2
!    write(*,5)n2,n3-1,l1-1,l2-1,l3-1,v3(n2,n3-1,l1-1,l2-1,l3-1)

 enddo

 write(ulog,*) 'leaving q1...',n1,q1
 write(ulog,'(a,f12.4)')'  cputime for matrix_elt ',cputim_matrix_elt





 enddo loop1

 nv3 = indx

 if (nv3 .ne. (fbz_subset(2)-fbz_subset(1)+1)*nk*ndn**3) then
   write(ulog,*) 'nv3 does not match in v3 calculation'
   stop
 endif

write(ulog,*)' V33: total size of this array is =',nv3


close(uv3)

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))
8 format(5(i6),2(2x,g14.8))
9 format(99(2x,g14.3))
 end subroutine calculate_w3_fbz_split_sq




!===========================================================================
subroutine read_col_matrix(nk,ndn)

use exactBTE2

implicit none

integer ucol, i, n1,n2,n3,l1,l2,l3
integer nk,ndn
ucol=9200



open(ucol,file='col_matrix.dat',status='unknown')
!read(ucol,*) ncol
ncol=nk**2*ndn**3

do i=1,ncol
   read(ucol,*) n1,n2,l1,l2,l3,P1(n1,n2,l1,l2,l3),P2(n1,n2,l1,l2,l3)
enddo

close(ucol)


end subroutine read_col_matrix



!===========================================================================
subroutine cal_Q(nk,ndn)

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

nkF1_loop: do i=1,nk

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



!=============================================================================
subroutine RTA(nk,kp,ndn,tempk,veloc,eigenval)
! calculate F value using only diagonal terms in collision matrix (equivalent as RTA)

use constants
use exactBTE2
use params
use lattice

implicit none

integer i,ii,iii,indx,nk,ndn
real(8) nbe,temp,tempk,nb,omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk), kp(3,nk)

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

nkF1_loop: do i=1,nk
    
    laF1_loop: do ii=1,ndn

          !omega=sqrt(eigenval(ii,i)) * cnst
         
          !nb=nbe(omega,temp,classical)
          omega=frequency(i,ii)
          nb=dist(i,ii)        
         
          xyz_loop: do iii=1,3

                F_RTA(i,ii,iii)= -(c_light*veloc(iii,ii,i)) * h_plank * (omega*c_light*100) * nb * (nb+1) / (k_b*tempk**2) / (Qvalue(i,ii)*c_light*100)
                ! veloc in c0, h_plank in J*s, omega in cm^-1 (linear freq.), k_b in J/K, tempk in K, Qvalue in cm^-1
                ! resulting F1 is in m/K

!write(*,*) 'indx,F1',indx,F1(indx,iii)
!                tau(i,ii,iii)= nb*(nb+1)/(Qvalue(i,ii))      ! tau is in cm

!                write(uRTA,10) tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i)

          enddo xyz_loop

    enddo laF1_loop

enddo nkF1_loop

!F1(:,:,:)=F_RTA(:,:,:) ! for iteration

6 format(2i6,2x,99(1x,f9.3))

end subroutine RTA



!============================================================================
subroutine mode_thermal_iter(nki,nk,ndn,tempk,veloc,eigenval) 
! calculate mode thermal conductivity for RTA and full solution of BTE

use constants
use exactBTE2
use params
use lattice
use io2
use kpoints !!mine
implicit none

integer i,ii,j,jj,indx,nk,ndn
real(8) nbe, temp, tempk, nbb, omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk)

real(8) tot_kap
integer nki
!integer ukapiter
!ukapiter = 9400
tot_kap=0

kappa=0
kappa_k=0
kappa_RTA=0
kappa_k_RTA=0
diff_kap=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1


!write(debug,*) 'i,ii,j,jj,om,veloc_j,nb,F2(indx,jj)'

write(ulog,*) 'entering mode_thermal_iter...'

nkF1_loop: do i=1,nki

!             q=kibbz(:,i)
!            call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
!             call kstarsinfbz1(q,nk,kp,narms,kvecstar,kvecop)
!             do s=1,narms
!                call syop_f(kvecop(s),F1(,l2,:),Fsq2)


      laF1_loop: do ii=1,ndn

          !!omega=sqrt(eigenval(ii,i)) * cnst !!mine command
          omega=sqrt(eigenval(ii,mapinv(i))) * cnst

          nbb=nbe(omega,temp,classical)

          alpha_loop: do j=1,3
          beta_loop: do jj=1,3

               !!kappa_k(i,ii,j,jj) = -h_plank * (omega*c_light*100) * (veloc(j,ii,i)*c_light) * nb * (nb+1) * F1(i,ii,jj)
               kappa_k(mapinv(i),ii,j,jj) = -h_plank * (omega*c_light*100) * (veloc(j,ii,mapinv(i))*c_light) * nbb * (nbb+1) * F1(mapinv(i),ii,jj)
               kappa_k(mapinv(i),ii,j,jj) = kappa_k(mapinv(i),ii,j,jj) / ( nk*volume_r/1d30) 
               kappa(ii,j,jj)=kappa(ii,j,jj) + kappa_k(mapinv(i),ii,j,jj)  

               diff_kap(mapinv(i),ii,j,jj)= -h_plank * (omega*c_light*100) * (veloc(j,ii,mapinv(i))*c_light) * nbb * (nbb+1) * (F1(mapinv(i),ii,jj)-F1_old(mapinv(i),ii,jj))
               diff_kap(mapinv(i),ii,j,jj)= diff_kap(mapinv(i),ii,j,jj) / (volume_r/1d30)   ! didn't divide by nk so that order of mode_kapp ~ order of total_kapp            

               !!kappa_k_RTA(i,ii,j,jj) = -h_plank * (omega*c_light*100) * (veloc(j,ii,i)*c_light) * nb * (nb+1) * F_RTA(i,ii,jj)
               !!kappa_k_RTA(i,ii,j,jj) = -h_plank * (omega*c_light*100) * (veloc(j,ii,mapinv(i))*c_light) * nbb * (nbb+1) * wibz(i)  * F_RTA(mapinv(i),ii,jj)
               kappa_k_RTA(mapinv(i),ii,j,jj) = -h_plank * (omega*c_light*100) * (veloc(j,ii,mapinv(i))*c_light) * nbb * (nbb+1) *  F_RTA(mapinv(i),ii,jj)
               kappa_k_RTA(mapinv(i),ii,j,jj) = kappa_k_RTA(mapinv(i),ii,j,jj) / (nk*volume_r/1d30)
               kappa_RTA(ii,j,jj)=kappa_RTA(ii,j,jj) + kappa_k_RTA(mapinv(i),ii,j,jj)

!write(debug,*) '--i,ii,j,jj'
!write(debug,*) i,ii,j,jj
!write(debug,*) 'i,ii,j,jj,indx,om,veloc_j,nb,F2(indx,jj)'
!write(debug,3) i,ii,j,jj,indx,omega,veloc(j,ii,i),nb,F1(i,ii,jj)
               write(129,*)i,mapinv(i),ii,j,jj, omega, veloc(j,ii,mapinv(i)),F1(mapinv(i),ii,jj),F_RTA(mapinv(i),ii,jj),nbb,kappa_k_RTA(mapinv(i),ii,j,jj),kappa_k(mapinv(i),ii,j,jj)
               write(129,*)kappa_RTA(ii,j,jj)
          enddo beta_loop

          if (veloc(j,ii,mapinv(i)) .eq. 0.0) then       ! Gamma point has zero velocity 
                tauinv_eff(mapinv(i),ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F1(mapinv(i),ii,j)/(0.01*c_light)
          else
                tauinv_eff(mapinv(i),ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F1(mapinv(i),ii,j)/(veloc(j,ii,mapinv(i))*c_light)
          end if

          enddo alpha_loop

       enddo laF1_loop
!enddo
enddo nkF1_loop


!open(ukapiter,file='mode_kappa_iter.dat',status='unknown')
!open(ukapiter+1,file='kappa_iter.dat',status='unknown')

!write(ukapiter,*) 'nk, la, al, be, mode_kap, mode_kap_RTA'
!write(ukapiter+1,*) 'la, kxx,kyy,kzz, kxy,kxz,kyz, k_RTA'

!do j=1,3
!do jj=1,3
!do ii=1,ndn
!do i=1,nk

!write(ukapiter,8) i,ii,j,jj,kappa_k(i,ii,j,jj),kappa_k_RTA(i,ii,j,jj)

!enddo

!write(ukapiter+1,*) ii,j,jj,kappa(ii,j,jj)
!if (j.eq.jj) then
!tot_kap=tot_kap+kappa(ii,j,jj)
!endif

!enddo
!enddo
!enddo

!do ii=1,ndn
!write(ukapiter+1,8) ii,(kappa(ii,j,j),j=1,3),((kappa(ii,j,jj),jj=j+1,3),j=1,2), (kappa_RTA(ii,j,j),j=1,3),((kappa_RTA(ii,j,jj),jj=j+1,3),j=1,2)
!enddo

!write(*,*) 'tot_kap=', tot_kap/3

3 format(99(1x,g11.4))
8 format(99(1x,g11.4))
end subroutine mode_thermal_iter
!====================================================================================

!============================================================================
subroutine mode_thermal_iter2(nki,nk,kibbz,kp,ndn,tempk,veloc,eigenval) !!Saf!!
! calculate mode thermal conductivity for RTA and full solution of BTE

use constants
use exactBTE2
use params
use lattice
use io2
use kpoints !!mine
implicit none

integer i,ii,j,jj,indx,nk,ndn
real(8) nbe, temp, tempk, nbb, omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk)

real(8) tot_kap,FRTA(3),q(3),kp(3,nk),kibbz(3,nki),kvecstar(3,48),F12(3)
integer nki,s, narms,kvecop(48),n0,i3,j3,k3, inside
real(8) velocc(3)
!integer ukapiter
!ukapiter = 9400
tot_kap=0

kappa=0
kappa_k=0
kappa_RTA=0
kappa_k_RTA=0
diff_kap=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1


!write(debug,*) 'i,ii,j,jj,om,veloc_j,nb,F2(indx,jj)'

write(ulog,*) 'entering mode_thermal_iter...'

nkF1_loop: do i=1,nki

             q=kibbz(:,i)
!            call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
             call kstarsinfbz1(q,nk,kp,narms,kvecstar,kvecop)


      laF1_loop: do ii=1,ndn
           do s=1,narms
          call syop_f(kvecop(s),F_RTA(mapinv(i),ii,:),FRTA)
          call syop_f(kvecop(s),veloc(:,ii,mapinv(i)),velocc)
          call syop_f(kvecop(s),F1(mapinv(i),ii,:),F12)
          call get_k_info(kvecstar(:,s),NC,n0,i3,j3,k3,g1,g2,g3,inside)
          !!omega=sqrt(eigenval(ii,i)) * cnst !!mine command
          omega=sqrt(eigenval(ii,n0)) * cnst

          nbb=nbe(omega,temp,classical)

          alpha_loop: do j=1,3
          beta_loop: do jj=1,3

               !!kappa_k(i,ii,j,jj) = -h_plank * (omega*c_light*100) * (veloc(j,ii,i)*c_light) * nb * (nb+1) * F1(i,ii,jj)
               kappa_k(n0,ii,j,jj) = -h_plank * (omega*c_light*100) * (velocc(j)*c_light) * nbb * (nbb+1) * F12(jj)
               kappa_k(n0,ii,j,jj) = kappa_k(n0,ii,j,jj) / (nk*volume_r/1d30)
               kappa(ii,j,jj)=kappa(ii,j,jj) + kappa_k(n0,ii,j,jj)


              !!diff_kap(i,ii,j,jj)= -h_plank * (omega*c_light*100) * (veloc(j,ii,i)*c_light) * nb * (nb+1) * (F1(i,ii,jj)-F1_old(i,ii,jj))
               diff_kap(n0,ii,j,jj)= -h_plank * (omega*c_light*100) * (velocc(j)*c_light) * nbb * (nbb+1) *(F12(jj)-F1_old(n0,ii,jj))
               diff_kap(n0,ii,j,jj)= diff_kap(n0,ii,j,jj) /(volume_r/1d30)   ! didn't divide by nk so that order of mode_kapp ~ order of total_kapp

               !!kappa_k_RTA(i,ii,j,jj) = -h_plank * (omega*c_light*100) * (veloc(j,ii,i)*c_light) * nb * (nb+1) * F_RTA(i,ii,jj)
               !!kappa_k_RTA(i,ii,j,jj) = -h_plank * (omega*c_light*100) * (veloc(j,ii,mapinv(i))*c_light) * nbb * (nbb+1) * wibz(i)  * F_RTA(mapinv(i),ii,jj)
              !@ kappa_k_RTA(mapinv(i),ii,j,jj) = -h_plank * (omega*c_light*100) *(veloc(j,ii,mapinv(i))*c_light) * nbb * (nbb+1) *  F_RTA(mapinv(i),ii,jj)
!               kappa_k_RTA(n0,ii,j,jj) = -h_plank * (omega*c_light*100)*(veloc(j,ii,n0)*c_light) * nbb * (nbb+1) *  FRTA(jj)
               kappa_k_RTA(n0,ii,j,jj) = -h_plank * (omega*c_light*100)*(velocc(j)*c_light) * nbb * (nbb+1) *  FRTA(jj)
               kappa_k_RTA(n0,ii,j,jj) = kappa_k_RTA(n0,ii,j,jj) /(nk*volume_r/1d30)
               kappa_RTA(ii,j,jj)=kappa_RTA(ii,j,jj) +kappa_k_RTA(n0,ii,j,jj)

!write(debug,*) '--i,ii,j,jj'
!write(debug,*) i,ii,j,jj
!write(debug,*) 'i,ii,j,jj,indx,om,veloc_j,nb,F2(indx,jj)'
!write(debug,3) i,ii,j,jj,indx,omega,veloc(j,ii,i),nb,F1(i,ii,jj)

          enddo beta_loop

!          if (veloc(j,ii,n0) .eq. 0.0) then       ! Gamma point has zero velocity
!                tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F1(n0,ii,j)/(0.01*c_light)
           if (velocc(j) .eq. 0.0) then       ! Gamma point has zero velocity
                tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F12(j)/(0.01*c_light)
 
         else
!                tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F1(n0,ii,j)/(veloc(j,ii,n0)*c_light)
               tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F12(j)/(velocc(j)*c_light)
          end if

          enddo alpha_loop
         enddo
       enddo laF1_loop
!enddo
enddo nkF1_loop


!open(ukapiter,file='mode_kappa_iter.dat',status='unknown')
!open(ukapiter+1,file='kappa_iter.dat',status='unknown')

!write(ukapiter,*) 'nk, la, al, be, mode_kap, mode_kap_RTA'
!write(ukapiter+1,*) 'la, kxx,kyy,kzz, kxy,kxz,kyz, k_RTA'

!do j=1,3
!do jj=1,3
!do ii=1,ndn
!do i=1,nk

!write(ukapiter,8) i,ii,j,jj,kappa_k(i,ii,j,jj),kappa_k_RTA(i,ii,j,jj)

!enddo

!write(ukapiter+1,*) ii,j,jj,kappa(ii,j,jj)
!if (j.eq.jj) then
!tot_kap=tot_kap+kappa(ii,j,jj)
!endif

!enddo
!enddo
!enddo

!do ii=1,ndn
!write(ukapiter+1,8) ii,(kappa(ii,j,j),j=1,3),((kappa(ii,j,jj),jj=j+1,3),j=1,2),
!(kappa_RTA(ii,j,j),j=1,3),((kappa_RTA(ii,j,jj),jj=j+1,3),j=1,2)
!enddo

!write(*,*) 'tot_kap=', tot_kap/3

3 format(99(1x,g11.4))
8 format(99(1x,g11.4))
end subroutine mode_thermal_iter2
!=====================================================================================




!============================================================================
subroutine mode_thermal_noniter22(nki,nk,kibbz,kp,ndn,tempk,veloc,eigenval,Fi_x) !!Saf!!
! calculate mode thermal conductivity for RTA and full solution of BTE

use constants
use exactBTE2
use params
use lattice
use io2
use kpoints !!mine
implicit none

integer i,ii,j,jj,indx,nk,ndn
real(8) nbe, temp, tempk, nbb, omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk)

real(8) tot_kap,FRTA(3),q(3),kp(3,nk),kibbz(3,nki),kvecstar(3,48),F12(3)
integer nki,s, narms,kvecop(48),n0,i3,j3,k3, inside
real(8) FF1(nki,ndn,3),Fi_x(nki*ndn*3), velocc(3)
integer i1,j1,j2,r
!integer ukapiter
!ukapiter = 9400
tot_kap=0

kappa=0
kappa_k=0
kappa_RTA=0
kappa_k_RTA=0
diff_kap=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1


!write(debug,*) 'i,ii,j,jj,om,veloc_j,nb,F2(indx,jj)'

write(ulog,*) 'entering mode_thermal_iter...'


FF1=0
r=0
do i1=1,nki
   do j1=1,ndn
      do j2=1,3
!         if (s.lt.kl) then
         r=r+1
         FF1(i1,j1,j2)=Fi_x(r)
!         endif
       enddo
    enddo
enddo



nkF1_loop: do i=1,nki

             q=kibbz(:,i)
!            call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
             call kstarsinfbz1(q,nk,kp,narms,kvecstar,kvecop)
!             do s=1,narms
!                call syop_f(kvecop(s),F1(,l2,:),Fsq2)


      laF1_loop: do ii=1,ndn
           do s=1,narms
         call syop_f(kvecop(s),F_RTA(mapinv(i),ii,:),FRTA)
        ! laF1_loop: do ii=1,ndn 
          call syop_f(kvecop(s),FF1(i,ii,:),F12)

          call syop_f(kvecop(s),veloc(:,ii,mapinv(i)),velocc)

          call get_k_info(kvecstar(:,s),NC,n0,i3,j3,k3,g1,g2,g3,inside)
          !!omega=sqrt(eigenval(ii,i)) * cnst !!mine command
          omega=sqrt(eigenval(ii,n0)) * cnst

          nbb=nbe(omega,temp,classical)

          alpha_loop: do j=1,3
          beta_loop: do jj=1,3

               !!kappa_k(i,ii,j,jj) = -h_plank * (omega*c_light*100) * (veloc(j,ii,i)*c_light) * nb * (nb+1) * F1(i,ii,jj)
!               kappa_k(n0,ii,j,jj) = -h_plank * (omega*c_light*100)*(veloc(j,ii,n0)*c_light) * nbb * (nbb+1) * F12(jj)
                kappa_k(n0,ii,j,jj) = -h_plank * (omega*c_light*100)*(velocc(j)*c_light) * nbb * (nbb+1) * F12(jj)
               kappa_k(n0,ii,j,jj) = kappa_k(n0,ii,j,jj) / (nk*volume_r/1d30)
               kappa(ii,j,jj)=kappa(ii,j,jj) + kappa_k(n0,ii,j,jj)

           !    diff_kap(n0,ii,j,jj)= -h_plank * (omega*c_light*100)*(veloc(j,ii,n0)*c_light) * nbb * (nbb+1) *(F12(jj)-F1_old(n0,ii,jj))
           !    diff_kap(n0,ii,j,jj)= diff_kap(n0,ii,j,jj) /(volume_r/1d30)   !didn't divide by nk so that order of mode_kapp ~ order of total_kapp

               !!kappa_k_RTA(i,ii,j,jj) = -h_plank * (omega*c_light*100) *
               !(veloc(j,ii,i)*c_light) * nb * (nb+1) * F_RTA(i,ii,jj)
               !!kappa_k_RTA(i,ii,j,jj) = -h_plank * (omega*c_light*100) *
               !(veloc(j,ii,mapinv(i))*c_light) * nbb * (nbb+1) * wibz(i)  *
               !F_RTA(mapinv(i),ii,jj)
              !@ kappa_k_RTA(mapinv(i),ii,j,jj) = -h_plank * (omega*c_light*100)
              !*(veloc(j,ii,mapinv(i))*c_light) * nbb * (nbb+1) *
              !F_RTA(mapinv(i),ii,jj)
!               kappa_k_RTA(n0,ii,j,jj) = -h_plank *(omega*c_light*100)*(veloc(j,ii,n0)*c_light) * nbb * (nbb+1) *  FRTA(jj)
               kappa_k_RTA(n0,ii,j,jj) = -h_plank * (omega*c_light*100)*(velocc(j)*c_light) * nbb * (nbb+1) *  FRTA(jj)
               kappa_k_RTA(n0,ii,j,jj) = kappa_k_RTA(n0,ii,j,jj)/(nk*volume_r/1d30)
               kappa_RTA(ii,j,jj)=kappa_RTA(ii,j,jj) +kappa_k_RTA(n0,ii,j,jj)

!write(debug,*) '--i,ii,j,jj'
!write(debug,*) i,ii,j,jj
!write(debug,*) 'i,ii,j,jj,indx,om,veloc_j,nb,F2(indx,jj)'
!write(debug,3) i,ii,j,jj,indx,omega,veloc(j,ii,i),nb,F1(i,ii,jj)
               write(129,*)i,mapinv(i),ii,j,jj,omega,veloc(j,ii,mapinv(i)),F1(mapinv(i),ii,jj),FRTA(jj),nbb,kappa_k_RTA(mapinv(i),ii,j,jj),kappa_k(mapinv(i),ii,j,jj)
               write(129,*)kappa_RTA(ii,j,jj)
          enddo beta_loop

          !if (veloc(j,ii,n0) .eq. 0.0) then       ! Gamma point has zero velocity
           !     tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F1(n0,ii,j)/(0.01*c_light)
            if (velocc(j) .eq. 0.0) then       ! Gamma point has zero velocity
                 tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F12(j)/(0.01*c_light)

           else
            !    tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F1(n0,ii,j)/(veloc(j,ii,n0)*c_light)
                tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F12(j)/(velocc(j)*c_light)
           end if

          enddo alpha_loop
         enddo
       enddo laF1_loop
!enddo
enddo nkF1_loop


!open(ukapiter,file='mode_kappa_iter.dat',status='unknown')
!open(ukapiter+1,file='kappa_iter.dat',status='unknown')

!write(ukapiter,*) 'nk, la, al, be, mode_kap, mode_kap_RTA'
!write(ukapiter+1,*) 'la, kxx,kyy,kzz, kxy,kxz,kyz, k_RTA'

!do j=1,3
!do jj=1,3
!do ii=1,ndn
!do i=1,nk

!write(ukapiter,8) i,ii,j,jj,kappa_k(i,ii,j,jj),kappa_k_RTA(i,ii,j,jj)

!enddo

!write(ukapiter+1,*) ii,j,jj,kappa(ii,j,jj)
!if (j.eq.jj) then
!tot_kap=tot_kap+kappa(ii,j,jj)
!endif

!enddo
!enddo
!enddo

!do ii=1,ndn
!write(ukapiter+1,8) ii,(kappa(ii,j,j),j=1,3),((kappa(ii,j,jj),jj=j+1,3),j=1,2),
!(kappa_RTA(ii,j,j),j=1,3),((kappa_RTA(ii,j,jj),jj=j+1,3),j=1,2)
!enddo

!write(*,*) 'tot_kap=', tot_kap/3

3 format(99(1x,g11.4))
8 format(99(1x,g11.4))
end subroutine mode_thermal_noniter22
!=====================================================================================


!============================================================================
subroutine mode_thermal_noniter23(nki,nk,kibbz,kp,ndn,tempk,veloc,eigenval,Fi_x) !!Saf!!
! calculate mode thermal conductivity for RTA and full solution of BTE

use constants
use exactBTE2
use params
use lattice
use io2
use kpoints !!mine
implicit none

integer i,ii,j,jj,indx,nk,ndn
real(8) nbe, temp, tempk, nbb, omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk)

real(8) tot_kap,FRTA(3),q(3),kp(3,nk),kibbz(3,nki),kvecstar(3,48),F12(3)
integer nki,s, narms,kvecop(48),n0,i3,j3,k3, inside
real(8) FF1(nki,ndn,3),Fi_x(nki*ndn*3), velocc(3), veloci(3,ndn,nk)
integer i1,j1,j2,r
!integer ukapiter
!ukapiter = 9400
tot_kap=0

kappa=0
kappa_k=0
kappa_RTA=0
kappa_k_RTA=0
diff_kap=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1


!write(debug,*) 'i,ii,j,jj,om,veloc_j,nb,F2(indx,jj)'

write(ulog,*) 'entering mode_thermal_iter...'


FF1=0
r=0
do i1=1,nki
   do j1=1,ndn
      do j2=1,3
!         if (s.lt.kl) then
         r=r+1
         FF1(i1,j1,j2)=Fi_x(r)
!         endif
       enddo
    enddo
enddo

call veloc_avergeIBZ_invsyFBZ(nk,kp,nki,kibbz,ndn,veloc,veloci)

nkF1_loop: do i=1,nki

             q=kibbz(:,i)
!            call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
             call kstarsinfbz1(q,nk,kp,narms,kvecstar,kvecop)
!             do s=1,narms
!                call syop_f(kvecop(s),F1(,l2,:),Fsq2)


      laF1_loop: do ii=1,ndn
           do s=1,narms
         call syop_f(kvecop(s),F_RTA(mapinv(i),ii,:),FRTA)
        ! laF1_loop: do ii=1,ndn
          call syop_f(kvecop(s),FF1(i,ii,:),F12)

          call syop_f(kvecop(s),veloci(:,ii,i),velocc)

          call get_k_info(kvecstar(:,s),NC,n0,i3,j3,k3,g1,g2,g3,inside)
          !!omega=sqrt(eigenval(ii,i)) * cnst !!mine command
          omega=sqrt(eigenval(ii,n0)) * cnst

          nbb=nbe(omega,temp,classical)

          alpha_loop: do j=1,3
          beta_loop: do jj=1,3

               !!kappa_k(i,ii,j,jj) = -h_plank * (omega*c_light*100) *
               !(veloc(j,ii,i)*c_light) * nb * (nb+1) * F1(i,ii,jj)
!               kappa_k(n0,ii,j,jj) = -h_plank *
!               (omega*c_light*100)*(veloc(j,ii,n0)*c_light) * nbb * (nbb+1) *
!               F12(jj)
                kappa_k(n0,ii,j,jj) = -h_plank * (omega*c_light*100)*(velocc(j)*c_light) * nbb * (nbb+1) * F12(jj)
               kappa_k(n0,ii,j,jj) = kappa_k(n0,ii,j,jj) / (nk*volume_r/1d30)
               kappa(ii,j,jj)=kappa(ii,j,jj) + kappa_k(n0,ii,j,jj)

           !    diff_kap(n0,ii,j,jj)= -h_plank *
           !    (omega*c_light*100)*(veloc(j,ii,n0)*c_light) * nbb * (nbb+1)
           !    *(F12(jj)-F1_old(n0,ii,jj))
           !    diff_kap(n0,ii,j,jj)= diff_kap(n0,ii,j,jj) /(volume_r/1d30)
           !    !didn't divide by nk so that order of mode_kapp ~ order of
           !    total_kapp

               !!kappa_k_RTA(i,ii,j,jj) = -h_plank * (omega*c_light*100) *
               !(veloc(j,ii,i)*c_light) * nb * (nb+1) * F_RTA(i,ii,jj)
               !!kappa_k_RTA(i,ii,j,jj) = -h_plank * (omega*c_light*100) *
               !(veloc(j,ii,mapinv(i))*c_light) * nbb * (nbb+1) * wibz(i)  *
               !F_RTA(mapinv(i),ii,jj)
              !@ kappa_k_RTA(mapinv(i),ii,j,jj) = -h_plank * (omega*c_light*100)
              !*(veloc(j,ii,mapinv(i))*c_light) * nbb * (nbb+1) *
              !F_RTA(mapinv(i),ii,jj)
!               kappa_k_RTA(n0,ii,j,jj) = -h_plank
!               *(omega*c_light*100)*(veloc(j,ii,n0)*c_light) * nbb * (nbb+1) *
!               FRTA(jj)
               kappa_k_RTA(n0,ii,j,jj) = -h_plank * (omega*c_light*100)*(velocc(j)*c_light) * nbb * (nbb+1) *  FRTA(jj)
               kappa_k_RTA(n0,ii,j,jj) = kappa_k_RTA(n0,ii,j,jj)/(nk*volume_r/1d30)
               kappa_RTA(ii,j,jj)=kappa_RTA(ii,j,jj) +kappa_k_RTA(n0,ii,j,jj)

!write(debug,*) '--i,ii,j,jj'
!write(debug,*) i,ii,j,jj
!write(debug,*) 'i,ii,j,jj,indx,om,veloc_j,nb,F2(indx,jj)'
!write(debug,3) i,ii,j,jj,indx,omega,veloc(j,ii,i),nb,F1(i,ii,jj)
!***           write(129,*)i,mapinv(i),ii,j,jj,omega,veloc(j,ii,mapinv(i)),F1(mapinv(i),ii,jj),FRTA(jj),nbb,kappa_k_RTA(mapinv(i),ii,j,jj),kappa_k(mapinv(i),ii,j,jj)
!***           write(129,*)kappa_RTA(ii,j,jj)
          enddo beta_loop

          !if (veloc(j,ii,n0) .eq. 0.0) then       ! Gamma point has zero
          !velocity
           !     tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F1(n0,ii,j)/(0.01*c_light)
            if (velocc(j) .eq. 0.0) then       ! Gamma point has zero velocity
                 tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F12(j)/(0.01*c_light)

           else
            !    tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F1(n0,ii,j)/(veloc(j,ii,n0)*c_light)
                tauinv_eff(n0,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*F12(j)/(velocc(j)*c_light)
           end if

          enddo alpha_loop
         enddo
       enddo laF1_loop
!enddo
enddo nkF1_loop


!open(ukapiter,file='mode_kappa_iter.dat',status='unknown')
!open(ukapiter+1,file='kappa_iter.dat',status='unknown')

!write(ukapiter,*) 'nk, la, al, be, mode_kap, mode_kap_RTA'
!write(ukapiter+1,*) 'la, kxx,kyy,kzz, kxy,kxz,kyz, k_RTA'

!do j=1,3
!do jj=1,3
!do ii=1,ndn
!do i=1,nk

!write(ukapiter,8) i,ii,j,jj,kappa_k(i,ii,j,jj),kappa_k_RTA(i,ii,j,jj)

!enddo

!write(ukapiter+1,*) ii,j,jj,kappa(ii,j,jj)
!if (j.eq.jj) then
!tot_kap=tot_kap+kappa(ii,j,jj)
!endif

!enddo
!enddo
!enddo

!do ii=1,ndn
!write(ukapiter+1,8) ii,(kappa(ii,j,j),j=1,3),((kappa(ii,j,jj),jj=j+1,3),j=1,2),
!(kappa_RTA(ii,j,j),j=1,3),((kappa_RTA(ii,j,jj),jj=j+1,3),j=1,2)
!enddo

!write(*,*) 'tot_kap=', tot_kap/3

3 format(99(1x,g11.4))
8 format(99(1x,g11.4))
end subroutine mode_thermal_noniter23
!=====================================================================================





















!============================================================================
subroutine mode_thermal_noniter2(nk,ndn,tempk,veloc,eigenval,F_MRHS2)
! calculate mode thermal conductivity for RTA and full solution of BTE

use constants
use exactBTE2
use params
use lattice
use io2

implicit none

integer i,ii,j,jj,indx,nk,ndn
real(8) nbe, temp, tempk, nb, omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk)

real(8) tot_kap
real(8) F_MRHS2(nk*ndn*3),FF1(nk,ndn,3)
integer i1,j1,j2,s,kl
!integer ukapiter
!ukapiter = 9400
tot_kap=0

kappa=0
kappa_k=0
kappa_RTA=0
kappa_k_RTA=0
diff_kap=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1


!write(debug,*) 'i,ii,j,jj,om,veloc_j,nb,F2(indx,jj)'

write(ulog,*) 'entering mode_thermal_iter...'
kl=nk*ndn*3
FF1=0
s=0
do i1=1,nk
   do j1=1,ndn
      do j2=1,3
!         if (s.lt.kl) then
         s=s+1
         FF1(i1,j1,j2)=F_MRHS2(s)
!         endif
       enddo
    enddo
enddo




nkF1_loop: do i=1,nk

      laF1_loop: do ii=1,ndn

          omega=sqrt(eigenval(ii,i)) * cnst
          nb=nbe(omega,temp,classical)

          alpha_loop: do j=1,3
          beta_loop: do jj=1,3

               kappa_k(i,ii,j,jj) = -h_plank * (omega*c_light*100) * (veloc(j,ii,i)*c_light) * nb * (nb+1) * FF1(i,ii,jj)
               kappa_k(i,ii,j,jj) = kappa_k(i,ii,j,jj) / ( nk*volume_r/1d30)
               kappa(ii,j,jj)=kappa(ii,j,jj) + kappa_k(i,ii,j,jj)

               !!!diff_kap(i,ii,j,jj)= -h_plank * (omega*c_light*100) * (veloc(j,ii,i)*c_light) * nb * (nb+1) * (F1(i,ii,jj)-F1_old(i,ii,jj))
               !!!diff_kap(i,ii,j,jj)= diff_kap(i,ii,j,jj) / (volume_r/1d30)   !didn't divide by nk so that order of mode_kapp ~ order of total_kapp

               kappa_k_RTA(i,ii,j,jj) = -h_plank * (omega*c_light*100) *(veloc(j,ii,i)*c_light) * nb * (nb+1) * F_RTA(i,ii,jj)
               kappa_k_RTA(i,ii,j,jj) = kappa_k_RTA(i,ii,j,jj) /(nk*volume_r/1d30)
               kappa_RTA(ii,j,jj)=kappa_RTA(ii,j,jj) + kappa_k_RTA(i,ii,j,jj)

              write(129,*) i,ii,j,jj, omega, veloc(j,ii,i),FF1(i,ii,jj),F_RTA(i,ii,jj)

!write(debug,*) '--i,ii,j,jj'
!write(debug,*) i,ii,j,jj
!write(debug,*) 'i,ii,j,jj,indx,om,veloc_j,nb,F2(indx,jj)'
!write(debug,3) i,ii,j,jj,indx,omega,veloc(j,ii,i),nb,F1(i,ii,jj)


          enddo beta_loop

          if (veloc(j,ii,i) .eq. 0.0) then       ! Gamma point has zero velocity
                tauinv_eff(i,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*FF1(i,ii,j)/(0.01*c_light)
          else
                tauinv_eff(i,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*FF1(i,ii,j)/(veloc(j,ii,i)*c_light)
          end if

          enddo alpha_loop

       enddo laF1_loop
enddo nkF1_loop


!open(ukapiter,file='mode_kappa_iter.dat',status='unknown')
!open(ukapiter+1,file='kappa_iter.dat',status='unknown')

!write(ukapiter,*) 'nk, la, al, be, mode_kap, mode_kap_RTA'
!write(ukapiter+1,*) 'la, kxx,kyy,kzz, kxy,kxz,kyz, k_RTA'

!do j=1,3
!do jj=1,3
!do ii=1,ndn
!do i=1,nk

!write(ukapiter,8) i,ii,j,jj,kappa_k(i,ii,j,jj),kappa_k_RTA(i,ii,j,jj)

!enddo

!write(ukapiter+1,*) ii,j,jj,kappa(ii,j,jj)
!if (j.eq.jj) then
!tot_kap=tot_kap+kappa(ii,j,jj)
!endif

!enddo
!enddo
!enddo

!do ii=1,ndn
!write(ukapiter+1,8) ii,(kappa(ii,j,j),j=1,3),((kappa(ii,j,jj),jj=j+1,3),j=1,2),
!(kappa_RTA(ii,j,j),j=1,3),((kappa_RTA(ii,j,jj),jj=j+1,3),j=1,2)
!enddo

!write(*,*) 'tot_kap=', tot_kap/3

3 format(99(1x,g11.4))
8 format(99(1x,g11.4))
end subroutine mode_thermal_noniter2
!=================================================                                      


































!=====================================================================================
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

 integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
 integer nk1,nk2,nk3, nk1n, nk2n
 integer ksubset(2),nk
 real(8) q1(3), q2(3), q3(3), q1n(3), q2n(3)
 real(8) kp(3,nk)

 integer i,n1,n2,n3,l1,l2,l3,ndn,indx
 real(8) temp     ! temperature in cm^-1
 real(8) eigenval(ndn,nk)
 real(8) omega1, omega2, omega3, V3sq1, V3sq2, col1, col2, nbe, nb2, nb3

 real(8) delta_l
 !real(8) Qvalue

 integer iself2
 iself2=9011

! open(ucol,file='col_matrix.dat',STATUS='unknown',FORM='unformatted')
 open(iself2,file='iself2.dat',STATUS='unknown')

write(*,*) 'entering selfenergy2...'
write(*,*) 'v3_threshold=',v3_threshold

ncol=(ksubset(2)-ksubset(1)+1)*nk*ndn**3
write(iself2,*) ncol

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
iselfenergy(n1,l1)=iselfenergy(n1,l1) + (pi/2) * V3sq1 * (2*(nb2-nb3)*delta_l(omega1+omega2-omega3,etaz) + (1+nb2+nb3)*delta_l(omega1-omega2-omega3,etaz)) / (2*pi)
!iselfenergy(n1,l1)=iselfenergy(n1,l1) + (pi/2) * V3sq1 * ((nb3-nb2)*(delta_l(omega1+omega2-omega3,etaz)-delta_l(omega1-omega2+omega3,etaz)) + (1+nb2+nb3)*(-delta_l(omega1-omega2-omega3,etaz)) + (1+nb2+nb3)*delta_l(omega1+omega2+omega3,etaz)) / (2*pi)

!iselfenergy(n1,l1)=iselfenergy(n1,l1) + (pi/2) * V3sq1 * ((nb3-nb2)*(delta_l(omega1+omega2-omega3,etaz)-delta_l(omega1-omega2+omega3,etaz)) + (1+nb2+nb3)*(-delta_l(omega1-omega2-omega3,etaz)) ) / (2*pi)

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

               Avalue(n1,l1,xyz)=Avalue(n1,l1,xyz) + P1(n1,n2,l1,l2,l3)*(F1(nk3_proc1,l3,xyz)-F1(n2,l2,xyz)) + 0.5*P2(n1,n2,l1,l2,l3)*(F1(n2,l2,xyz)+F1(nk3_proc2,l3,xyz))

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
subroutine check_conv(i,nk,kp,ndn,tempk,convergence,iter_cont)

use constants
use exactBTE2
use params
use lattice
use io2
use phi3
use mod_ksubset2
use eigen

implicit none

integer nk,ndn,i,convergence,iter_cont
integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
integer nk1,nk2,nk3, nk3_proc1, nk3_proc2,xyz,xyz2
real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3)
real(8) kp(3,nk)
real(8) Avalue(nk,ndn,3)
integer n1,n2,n3,l1,l2,l3
real(8) leftBTE, rightBTE  ! lefthand term in BTE (3.78)
real(8) omega1, nb1, nbe, temp, tempk, norm_diff, max_diff, norm_error, max_error, norm_diff_kap, max_diff_kap

integer nk_subset, k_shift, indx, ksubset(2), ksub_size2
integer max_diff_kap_nk, max_diff_kap_la, max_diff_kap_xyz, max_diff_kap_xyz2


temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

write(*,*) 'entering check_conv...'

Avalue=0
indx=0
loop1 : do n1=1,nk

    q1=kp(:,n1)
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
    if (n1 .ne. nk1) then
       write(ulog,*) 'n1,nk1,inside=',n1,nk1,inside
       write(ulog,*) 'q1=',q1
       write(ulog,*) 'ijk1=',i1,j1,k1
       stop
    endif

 if (mod(n1,ksub_size) .eq. 1) then   ! start point of each v33.xxx.dat

    if(allocated(P1)) then
        call deallocate_iter2
    endif

    indx=indx+1
    if (indx*ksub_size > nk) then
        ksub_size2=nk - (indx-1)*ksub_size
        ksubset(1)=ksub_size*(indx-1)+1    ! starting point of k subset
        ksubset(2)=nk
    else
        ksub_size2=ksub_size
        ksubset(1)=ksub_size*(indx-1)+1
        ksubset(2)=ksub_size*indx
    endif
    k_shift=ksub_size*(indx-1)    ! P1 of n1,n2,l1,l2,l3 will be saved in P1(n1-k_shift,n2,l1,l2,l3).

    call allocate_iter2(ksub_size2,nk,ndn,ndn,ndn)
    call Pmatrix3(indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)

!    write(*,*) 'Pmatrix3 done,n1=',n1
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

               Avalue(n1,l1,xyz)=Avalue(n1,l1,xyz) + P1(n1-k_shift,n2,l1,l2,l3)*(F1(nk3_proc1,l3,xyz)-F1(n2,l2,xyz)) + 0.5*P2(n1-k_shift,n2,l1,l2,l3)*(F1(n2,l2,xyz)+F1(nk3_proc2,l3,xyz))

           enddo
           enddo
        enddo loop2

        omega1=frequency(n1,l1)
        nb1=dist(n1,l1)

        rightBTE = (F1(n1,l1,xyz)*Qvalue(n1,l1) - Avalue(n1,l1,xyz)) * c_light*100   ! righthand term in BTE
        leftBTE = -(veloc(xyz,l1,n1)*c_light) * h_plank * (omega1*c_light*100) * nb1*(nb1+1) / (k_b*tempk**2)               ! lefthand term in BTE
        error(n1,l1,xyz)=rightBTE-leftBTE

    enddo loopxyz
    enddo


enddo loop1

!write(*,*) 'F1,rightBTE',F1(1,1,1),rightBTE
!write(*,*) 'leftBTE',leftBTE


!calculate norm(diff) and max(diff)
!check convergence
norm_error=0
max_error=0
norm_diff=0
max_diff=0
norm_diff_kap=0
max_diff_kap=0

do n1=1,nk
do l1=1,ndn
do xyz=1,3

norm_error=norm_error+(error(n1,l1,xyz)**2/(nk*ndn*3))   ! normalize with nk*ndn*3
norm_diff=norm_diff+(diff(n1,l1,xyz)**2/(nk*ndn*3))

if (abs(error(n1,l1,xyz)) .gt. max_error) then
    max_error=abs(error(n1,l1,xyz))
endif

if (abs(diff(n1,l1,xyz)) .gt. max_diff) then
    max_diff=abs(diff(n1,l1,xyz))
endif

do xyz2=1,3
norm_diff_kap=norm_diff_kap+(diff_kap(n1,l1,xyz,xyz2)**2/(nk*ndn*3*3))
if (abs(diff_kap(n1,l1,xyz,xyz2)) .gt. abs(max_diff_kap)) then
    max_diff_kap=diff_kap(n1,l1,xyz,xyz2)
    max_diff_kap_nk=n1
    max_diff_kap_la=l1
    max_diff_kap_xyz=xyz
    max_diff_kap_xyz2=xyz2
endif
enddo
!write(*,*) 'n1,l1,xyz,F1',n1,l1,xyz,F1(n1,l1,xyz)

enddo
enddo
enddo

norm_error=sqrt(norm_error)
norm_diff=sqrt(norm_diff)

convergence=0
if (norm_error.lt.conv_error .and. max_error.lt.conv_max_error .and. norm_diff.lt.conv_diff .and. max_diff.lt.conv_max_diff .and. norm_diff_kap.lt.conv_diff_kap .and. max_diff_kap.lt.conv_max_diff_kap) then
convergence=1
iter_cont=iter_cont+1
else
convergence=0
iter_cont=0
endif

write(ulog,'(a2,99(2x,g12.6))') '++',i,norm_diff_kap,max_diff_kap,max_diff_kap_nk,max_diff_kap_la,max_diff_kap_xyz,max_diff_kap_xyz2



end subroutine check_conv





!============================================================
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
 real(8) q3(3),mi,mj,mk,rr2(3),rr3(3),den, q3_proc1(3), q3_proc2(3), den_proc1, den_proc2
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
       rr2(:) = atompos(:,j) - atompos(:,j0)
       rr3(:) = atompos(:,k) - atompos(:,k0)
       eiqr_proc1 = exp( ci* ((q2 .dot. rr2) + (q3_proc1 .dot. rr3)) )
       xx_proc1 = xx_proc1 + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr_proc1 * &
&      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3_proc1)

       eiqr_proc2 = exp( ci* ((q2 .dot. rr2) + (q3_proc2 .dot. rr3)) )
       xx_proc2 = xx_proc2 + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr_proc2 * &
&      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3_proc2)


 enddo tloop
 den_proc1 = sqrt(8*sqrt(abs(eigenval(l1,nk1)*eigenval(l2,nk2)*eigenval(l3,nk3_proc1))))
 if (den_proc1.ne.0) then
    w33_proc1 = xx_proc1 / den_proc1 * const33
 else
    write(*,*)'MATRIX_ELT: den=0 ',den
    stop
 endif
 den_proc2 = sqrt(8*sqrt(abs(eigenval(l1,nk1)*eigenval(l2,nk2)*eigenval(l3,nk3_proc2))))
 if (den_proc2.ne.0) then
    w33_proc2 = xx_proc2 / den_proc2 * const33
 else
    write(*,*)'MATRIX_ELT: den=0 ',den
    stop
 endif

 end subroutine matrix_elt_sy



!=========================================================================================
subroutine write_thermal(iteration,nk,ndn,tempk,veloc,eigenval,kp)

use constants
use exactBTE2
use params
use lattice
use io2

implicit none

integer i,ii,j,jj,indx,nk,ndn,iteration

integer ukapiter,uRTA
real(8) tempk, nb, omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk), kp(3,nk)

ukapiter = 9401
uRTA=9405

write(ulog,*) 'entering write_thermal...'


open(uRTA,file='iself_RTA.dat',status='replace')

open(ukapiter,file='mode_kappa_iter.dat',status='replace')
open(ukapiter+1,file='kappa_iter.dat',status='unknown')

write(ukapiter,*) 'nk, la, al, be, mode_kap, mode_kap_RTA'
write(ukapiter+1,*) 'iteration',iteration
write(ukapiter+1,*) 'la, kxx,kyy,kzz, kxy,kxz,kyz, k_RTA'

do j=1,3
do jj=1,3
do ii=1,ndn
do i=1,nk

write(ukapiter,3) i,ii,j,jj,kappa_k(i,ii,j,jj),kappa_k_RTA(i,ii,j,jj)

enddo

!if (j.eq.jj) then
!tot_kap=tot_kap+kappa(ii,j,jj)
!endif

enddo
enddo
enddo

do ii=1,ndn
write(ukapiter+1,8) ii,(kappa(ii,j,j),j=1,3),((kappa(ii,j,jj),jj=j+1,3),j=1,2), (kappa_RTA(ii,j,j),j=1,3),((kappa_RTA(ii,j,jj),jj=j+1,3),j=1,2)
enddo
!write(*,*) 'tot_kap=', tot_kap/3

close(ukapiter)


nkF1_loop: do i=1,nk
    laF1_loop: do ii=1,ndn
          omega=frequency(i,ii)
          nb=dist(i,ii)

           if(tauinv_U(i,ii) .eq. 0) then
           tauinv_U(i,ii)=1.0d10
           endif
           if(tauinv_N(i,ii) .eq. 0) then
           tauinv_N(i,ii)=1.0d10
           endif

           write(uRTA,10) tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i),(kappa_k_RTA(i,ii,j,j),j=1,3),(kappa_k(i,ii,j,j),j=1,3)
!          write(uRTA,10) tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i),(kappa_k_RTA(i,ii,j,j),j=1,3),(kappa_k(i,ii,j,j),j=1,3),&
!&         (tauinv_eff(i,ii,j),j=1,3)

!!-- temporary
if (tauinv_U(i,ii) .eq. 0.0) then
write(ulog,*) 'zero tauinv_U, i, ii', i, ii
endif

    enddo laF1_loop
enddo nkF1_loop

close(uRTA)


10 FORMAT(99(1x,g14.8))
3 format(99(1x,g11.4))
8 format(i2,99(1x,g11.4))


end subroutine write_thermal
!===========================================================


!=========================================================================================
subroutine write_thermal_noniter2(nk,ndn,tempk,veloc,eigenval,kp) !!Saf!!

use constants
use exactBTE2
use params
use lattice
use io2

implicit none

integer i,ii,j,jj,indx,nk,ndn,iteration

integer ukapiter,uRTA
real(8) tempk, nb, omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk), kp(3,nk)

ukapiter = 9401
uRTA=9405

write(ulog,*) 'entering write_thermal...'


open(uRTA,file='iself_RTA_noniter.dat',status='replace')

open(ukapiter,file='mode_kappa_noniter.dat',status='replace')
open(ukapiter+1,file='kappa_noniter.dat',status='unknown')

write(ukapiter,*) 'nk, la, al, be, mode_kap, mode_kap_RTA'
!write(ukapiter+1,*) 'iteration',iteration
write(ukapiter+1,*) 'la, kxx,kyy,kzz, kxy,kxz,kyz, k_RTA'

do j=1,3
do jj=1,3
do ii=1,ndn
do i=1,nk

write(ukapiter,3) i,ii,j,jj,kappa_k(i,ii,j,jj),kappa_k_RTA(i,ii,j,jj)

enddo

!if (j.eq.jj) then
!tot_kap=tot_kap+kappa(ii,j,jj)
!endif

enddo
enddo
enddo

do ii=1,ndn
write(ukapiter+1,8) ii,(kappa(ii,j,j),j=1,3),((kappa(ii,j,jj),jj=j+1,3),j=1,2),(kappa_RTA(ii,j,j),j=1,3),((kappa_RTA(ii,j,jj),jj=j+1,3),j=1,2)
enddo
!write(*,*) 'tot_kap=', tot_kap/3

close(ukapiter)


nkF1_loop: do i=1,nk
    laF1_loop: do ii=1,ndn
          omega=frequency(i,ii)
          nb=dist(i,ii)

           if(tauinv_U(i,ii) .eq. 0) then
           tauinv_U(i,ii)=1.0d10
           endif
           if(tauinv_N(i,ii) .eq. 0) then
           tauinv_N(i,ii)=1.0d10
           endif

           write(uRTA,10) tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i),(kappa_k_RTA(i,ii,j,j),j=1,3),(kappa_k(i,ii,j,j),j=1,3)
!          write(uRTA,10)
!          tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i),(kappa_k_RTA(i,ii,j,j),j=1,3),(kappa_k(i,ii,j,j),j=1,3),&
!&         (tauinv_eff(i,ii,j),j=1,3)

!!-- temporary
if (tauinv_U(i,ii) .eq. 0.0) then
write(ulog,*) 'zero tauinv_U, i, ii', i, ii
endif

    enddo laF1_loop
enddo nkF1_loop

close(uRTA)


10 FORMAT(99(1x,g14.8))
3 format(99(1x,g11.4))
8 format(i2,99(1x,g11.4))


end subroutine write_thermal_noniter2
!===================================================



































!========================================================================================
subroutine Pmatrix2 (ksubset,nk,kp,ndn,eigenval,temp)

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

 implicit none

 integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
 integer nk1,nk2,nk3, nk1n, nk2n, nk3n, nk3_proc1, nk3_proc2, nr1,nr2,lr1,lr2,lr3
 integer ksubset(2),nk
 real(8) q1(3), q2(3), q3(3), q1n(3), q2n(3), q3n(3), q3_proc1(3), q3_proc2(3)
 real(8) kp(3,nk)

 integer i,n1,n2,n3,l1,l2,l3,ndn
 real(8) temp     ! temperature in cm^-1
 real(8) eigenval(ndn,nk)
 real(8) omega1, omega2, omega3, V3sq1, V3sq2, col1, col2, omega3_proc1, omega3_proc2
 real(8) xx


 integer ucol, unt, nv3_2
 ucol=9010
 unt=112

! open(unt,FILE='v33sq.dat',STATUS='old',FORM='UNFORMATTED')
 open(unt,file='v33sq.dat',status='old')
write(*,*) 'open unt...'

! read(unt) nv3_2
 read(unt,*) nv3_2
 write(*,*) 'read nv3_2'

 write(*,*) 'entering cal Pmatrix...'

! open(ucol,file='col_matrix.dat',STATUS='unknown',FORM='unformatted')
 open(ucol,file='col_matrix.dat',STATUS='unknown')

ncol=(ksubset(2)-ksubset(1)+1)*nk*ndn**3
write(ucol,*) ncol
!write(ucol) ncol

loop1 : do n1=ksubset(1),ksubset(2)


    q1=kp(:,n1)
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
    if (n1 .ne. nk1) then
       write(ulog,*) 'n3,nk3,inside=',n3,nk3,inside
       write(ulog,*) 'q3=',q3
       write(ulog,*) 'ijk3=',i3,j3,k3
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



           do l2=1, ndn
           do l3=1, ndn
!write(*,*) 'before reading nr1,nr2,,...'
      !           read(unt,*) nr1,nr2,lr1,lr2,lr3,xx
!                  read(unt) nr1,nr2,lr1,lr2,lr3,xx
!                  read(unt) xx
                   read(unt,*) xx
!write(*,*) 'read nr1,nr2,lr1...'

!                 if (nr1.ne.n1 .or. nr2.ne.n2 .or. lr1.ne.l1 .or. lr2.ne.l2 .or. lr3.ne.l3) then
!                      write(ulog,*) 'n1,n2,l1,l2,l3 read from V3 does not match. exit'
!                      stop
!                 else

                 V3sq1=xx * (2*pi)**2 / nk
                 V3sq2=xx * (2*pi)**2 / nk

                 omega1=sqrt(eigenval(l1,nk1)) * cnst    !linear freq.
                 omega2=sqrt(eigenval(l2,nk2)) * cnst
                 omega3=sqrt(eigenval(l3,nk3)) * cnst

                 q3_proc1=q1+q2
                 q3_proc2=q1-q2
                 call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside)
                 call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside)

                 omega3_proc1=sqrt(eigenval(l3,nk3_proc1))*cnst
                 omega3_proc2=sqrt(eigenval(l3,nk3_proc2))*cnst

!                 call find_V3(nk1n,nk2,nk3,l1,l2,l3,V3sq2)    ! return V3sq2

!                 V3sq1=v33sq1(nk1,nk2,l1,l2,l3) * (2*pi)**2 / nk
!V33_md(nk1,nk2,l1,l2,l3)=V33_md(nk1,nk2,nk3n,l1,l2,l3) since q3=-q1-q2 when calculating V33_md
!                 V3sq2=v33sq1(nk1,nk2,l1,l2,l3) * (2*pi)**2 / nk

!                 call cal_P(temp,omega1,omega2,omega3_proc1,omega3_proc2,V3sq1,col1,col2)   ! return P1, P2

!!k1
        !         call cal_P2(n1,l1,n2,l2,nk3_proc1,l3,V3sq1,col1,col2)
!!k1
!                 proc1(nk1,nk2,l1,l2,l3)=col1
!                 proc2(nk1,nk2,l1,l2,l3)=col2

!                 write(ucol,7) n1, n2, l1, l2, l3, col1, col2
!                  write(ucol) n1,n2,l1,l2,l3,col1,col2
!                  write(ucol) col1,col2
                   write(ucol,*) col1,col2               
!                 endif

            enddo
            enddo
!            enddo
     enddo loop2
write(ulog,*) 'nk1,l1', nk1, l1

enddo
enddo loop1

!do n1=ksubset(1),ksubset(2)
!do n2=1,nk
!do l1=1,ndn
!do l2=1,ndn
!do l3=1,ndn
!                 write(ucol,7) n1, n2, l1, l2, l3, proc1(n1,n2,l1,l2,l3), proc2(n1,n2,l1,l2,l3)
                 ! write(ucol) nk1, nk2, nk3, l1, l2, l3, col1, col2
!enddo
!enddo
!enddo
!enddo
!enddo

close(ucol)
close(unt)

7 format(5(i6),2(2x,g14.8))
8 format(i6,3(2x,g14.8))


end subroutine Pmatrix2




!============================================================================================
subroutine cal_Qm(nk,ndn)

use io2
use exactBTE2
use mod_ksubset2

implicit none

integer i,ii,indx,ndn,j,jj,kk,indx2
integer nk, nk_subset, k_shift

integer uQ
uQ=9402
open(uQ,file='Qvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Q'

Qvalue(:,:)=0
indx=1
indx2=1

nkF1_loop: do i=1,nk

    if (i .eq. ksubset2(indx2,1)) then
         
         k_shift=ksubset2(indx2,1)-1
         nk_subset=ksubset2(indx2,2)-ksubset2(indx2,1)+1
         call allocate_iter2(nk_subset,nk,ndn,ndn,ndn)
         call read_Pmatrix(indx2,nk,ndn,k_shift)

    endif    


    laF1_loop: do ii=1,ndn

           nkF2_loop: do j=1,nk
           laF2_loop: do jj=1,ndn
           laF3_loop: do kk=1,ndn

                Qvalue(i,ii)=Qvalue(i,ii) + P1(i-k_shift,j,ii,jj,kk) + 0.5*P2(i-k_shift,j,ii,jj,kk)

           enddo laF3_loop
           enddo laF2_loop
           enddo nkF2_loop

           write(uQ,*) i,ii,Qvalue(i,ii)

      enddo laF1_loop

      if (i .eq. ksubset2(indx2,2)) then
          indx2=indx2+1
          call deallocate_iter2
      endif

enddo nkF1_loop


if(indx2-1 .ne. num_pfiles) then
   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
   stop
endif


end subroutine cal_Qm




!===========================================================
subroutine cal_F2m(nk,kp,ndn)

use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2

implicit none

integer nk,ndn
integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
integer nk1,nk2,nk3, nk3_proc1, nk3_proc2,xyz
real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3)
real(8) kp(3,nk)
real(8) Avalue(nk,ndn,3)

integer n1,n2,n3,l1,l2,l3
integer indx2,k_shift,nk_subset


indx2=1
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

    if (n1 .eq. ksubset2(indx2,1)) then

         k_shift=ksubset2(indx2,1)-1
         nk_subset=ksubset2(indx2,2)-ksubset2(indx2,1)+1
         call allocate_iter2(nk_subset,nk,ndn,ndn,ndn)
         call read_Pmatrix(indx2,nk,ndn,k_shift)

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

               Avalue(n1,l1,xyz)=Avalue(n1,l1,xyz) + P1(n1-k_shift,n2,l1,l2,l3)*(F1(nk3_proc1,l3,xyz)-F1(n2,l2,xyz)) + 0.5*P2(n1-k_shift,n2,l1,l2,l3)*(F1(n2,l2,xyz)+F1(nk3_proc2,l3,xyz))

           enddo
           enddo
        enddo loop2

        F2(n1,l1,xyz)=F_RTA(n1,l1,xyz) + Avalue(n1,l1,xyz)/Qvalue(n1,l1)
        diff(n1,l1,xyz)=(F2(n1,l1,xyz)-F1(n1,l1,xyz))/F2(n1,l1,xyz)

    enddo loopxyz
    enddo

    if (n1 .eq. ksubset2(indx2,2)) then
         indx2=indx2+1
         call deallocate_iter2
    endif

enddo loop1

if(indx2-1 .ne. num_pfiles) then
   write(ulog,*) 'error in reading col_matrix.dat during cal_F2m...'
   stop
endif


write(*,*) F_RTA(1,1,1), F1(1,1,1), F2(1,1,1)


end subroutine cal_F2m




!===================================================================================================
subroutine read_ksubset2

use mod_ksubset2

implicit none

integer uksubset2,i

uksubset2=11000

open(uksubset2,file='ksubset2.inp',status='unknown')
read(uksubset2,*) num_pfiles

call allocate_ksubset2(num_pfiles)

do i=1,num_pfiles
   read(uksubset2,*) ksubset2(i,1),ksubset2(i,2)
enddo

close(uksubset2)

end subroutine read_ksubset2




!====================================================================================
subroutine read_Pmatrix(indx2,nk,ndn,k_shift)

use mod_ksubset2
use exactBTE2

implicit none
integer indx2, uP,nk,ndn
character(99) filename1
integer indx,i, nP
integer n1,n2,l1,l2,l3,k_shift

uP=12000

write(filename1,fmt="(a,i3.3,a)") "Pmatrix.",indx2,".dat"
!open(uP,file=filename1,status='unknown')
open(uP,file=filename1,status='unknown',form='unformatted')
nP=(ksubset2(indx2,2)-ksubset2(indx2,1)+1)*nk*ndn**3

!do i=1,nP
!    read(uP,*) n1,n2,l1,l2,l3,P1(n1-k_shift,n2,l1,l2,l3),P2(n1-k_shift,n2,l1,l2,l3)
!    read(uP) n1,n2,l1,l2,l3,P1(n1-k_shift,n2,l1,l2,l3),P2(n1-k_shift,n2,l1,l2,l3)
!enddo

do n1=ksubset2(indx2,1),ksubset2(indx2,2)
do l1=1,ndn
do n2=1,nk
do l2=1,ndn
do l3=1,ndn
read(uP) P1(n1-k_shift,n2,l1,l2,l3),P2(n1-k_shift,n2,l1,l2,l3)
enddo
enddo
enddo
enddo
enddo


end subroutine read_Pmatrix



!==========================================================================================
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



!============================================================================================
subroutine cal_Qm2(nk,ndn,kp,tempk)

use io2
use exactBTE2
use mod_ksubset2
use params
use constants
use phi3
use lattice

implicit none

integer i,ii,indx,ndn,j,jj,kk
integer nk, nk_subset, k_shift, ksub_size2
real(8) kp(3,nk), q1(3), q2(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
real(8) temp, tempk, nb1   ! in cm^-1
integer ksubset(2)

integer i2,ii2,j2,jj2,kk2, inside1, inside2, inside, inside1n, inside2n
integer i3,j3,k3,nk3_proc1,nk3_proc2, nk3

integer uQ

real cputim1, cputim2, cputim3

uQ=9402
open(uQ,file='Qvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Q'

!open(uQ+1,file='Pmatrix_paral.dat',status='unknown')

Qvalue(:,:)=0    ! in cm^-1
Qvalue_N(:,:)=0
Qvalue_U(:,:)=0
tauinv_N(:,:)=0
tauinv_U(:,:)=0
indx=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

!write(*,*) 'entering cal_Qm2...'

call cpu_time(cputim1)
write(*,*) 'entering cal_Qm2...'

nkF1_loop: do i=1,nk

 ! call Pmatrix3(i)
     ! -> read V33.xxx.dat
     ! -> allocate P1,P2 subset
     ! -> calculate P1, P2
     ! -> close V33.xxx.dat
 if (mod(i,ksub_size) .eq. 1) then   ! start point of each v33.xxx.dat

    if(allocated(P1)) then
        call deallocate_iter2
    endif

    indx=indx+1
    if (indx*ksub_size > nk) then
        ksub_size2=nk - (indx-1)*ksub_size
        ksubset(1)=ksub_size*(indx-1)+1    ! starting point of k subset
        ksubset(2)=nk
    else
        ksub_size2=ksub_size
        ksubset(1)=ksub_size*(indx-1)+1
        ksubset(2)=ksub_size*indx
    endif
    k_shift=ksub_size*(indx-1)    ! P1 of n1,n2,l1,l2,l3 will be saved in P1(n1-k_shift,n2,l1,l2,l3).

    call allocate_iter2(ksub_size2,nk,ndn,ndn,ndn)
    call Pmatrix3(indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)

    call cpu_time(cputim2)
    write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1

!    do i2=ksubset(1),ksubset(2)
!    do ii2=1,ndn
!    do j2=1,nk
!    do jj2=1,ndn
!    do kk2=1,ndn
!       write(uQ+1,7) i2,j2,ii2,jj2,kk2,P1(i2-k_shift,j2,ii2,jj2,kk2),P2(i2-k_shift,j2,ii2,jj2,kk2)
!    enddo
!    enddo
!    enddo
!    enddo
!    enddo

!    write(*,*) 'In cal_Qm2,Pmatrix3 done,i=',i
!    write(*,*) 'In cal_Qm2,k_shift',k_shift

 endif

 !    if (i .eq. ksubset2(indx2,1)) then

 !        k_shift=ksubset2(indx2,1)-1
 !        nk_subset=ksubset2(indx2,2)-ksubset2(indx2,1)+1
 !        call allocate_iter2(nk_subset,nk,ndn,ndn,ndn)
 !        call read_Pmatrix(indx2,nk,ndn,k_shift)

 !   endif

    laF1_loop: do ii=1,ndn

           nkF2_loop: do j=1,nk
           laF2_loop: do jj=1,ndn
           laF3_loop: do kk=1,ndn

                Qvalue(i,ii)=Qvalue(i,ii) + P1(i-k_shift,j,ii,jj,kk) + 0.5*P2(i-k_shift,j,ii,jj,kk)

                q1=kp(:,i)
                q2=kp(:,j)

                q3=-q1-q2
                q3_proc1=q1+q2
                q3_proc2=q1-q2
                q3_proc1n=-q1-q2
                q3_proc2n=-q1+q2

                call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)
                call get_k_info(q3_proc1n,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1n)
                call get_k_info(q3_proc2n,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2n)

                if (inside1 .eq. 1) then    ! normal process
                    Qvalue_N(i,ii)=Qvalue_N(i,ii) + P1(i-k_shift,j,ii,jj,kk)
                else
                    Qvalue_U(i,ii)=Qvalue_U(i,ii) + P1(i-k_shift,j,ii,jj,kk)
                endif

                if (inside2 .eq. 1) then   
                    Qvalue_N(i,ii)=Qvalue_N(i,ii) + 0.5*P2(i-k_shift,j,ii,jj,kk)
                else
                    Qvalue_U(i,ii)=Qvalue_U(i,ii) + 0.5*P2(i-k_shift,j,ii,jj,kk)
                endif

! determine k vector
! determine i,j,k is normal or Umklapp
! If normal, Qvalue_n = Qvalue_n + P1 + 0.5*P2
! If Umklapp, Qvalue_u = Qvalue_u + P1 + 0.5*P2

           enddo laF3_loop
           enddo laF2_loop
           enddo nkF2_loop

           write(uQ,*) i,ii,Qvalue(i,ii)

           nb1=dist(i,ii)
           tauinv_N(i,ii) = Qvalue_N(i,ii)/(nb1*(nb1+1)) * c_light*100    ! tauinv_N in 1/sec
           tauinv_U(i,ii) = Qvalue_U(i,ii)/(nb1*(nb1+1)) * c_light*100
           tauinv_tot(i,ii) = tauinv_N(i,ii) + tauinv_U(i,ii)

      enddo laF1_loop

! give a sign to move to next v33.xxx.dat file


!      if (i .eq. ksubset2(indx2,2)) then
!          indx2=indx2+1
!          call deallocate_iter2
!      endif

enddo nkF1_loop

close(uQ)
!close(uQ+1)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif


end subroutine cal_Qm2
!!==========================================================
subroutine get_k_infos(q,nq,nk,kp,ins)

use lattice 
!use kpoints
use params
use io2

implicit none
real(8) q(3), lk
real(8) kp(3,nk)
integer nq,nk,i 
logical ins

nq=0
ins=.false.
loop:  do i=1,nk
   lk =length(q-kp(:,i))
   if (lk.lt.1d-3) then
   nq = i 
   ins=.true.
   exit loop
   endif 
  enddo loop

 if (nq.eq.0) then
 write(ulog,*)'q was not found', q
!!! write(*,*)'q was not found', q
 !stop
 endif

end subroutine get_k_infos
!!==========================================================
!subroutine 

!end subroutine
!!=========================================================
subroutine get_kvecopn(q,k,spf,nk,ms)

use constants
use exactBTE2
use lattice
use params
use kpoints
use io2
implicit none
real(8) q(3), kvecstar(3,48),xp, k(3)
integer n,s,narms,sp,spf,kvecop(48),nk,sd,sk
integer i3,j3,k3, inside, ms

          call getstar(q,primitivelattice,narms,kvecstar,kvecop)
!!!***         call get_weights(nk,q,narms,kvecstar,kvecop)
          !write(*,*)s, narmsp1
         call get_k_info(k,NC,sk,i3,j3,k3,g1,g2,g3,inside)
!!!         write(*,*)'q3-k', k , i3,j3,k3,inside
 loop: do s=1,narms

             sp=kvecop(s)
          call get_k_info(kvecstar(:,s),NC,sd,i3,j3,k3,g1,g2,g3,inside)
!!!          write(*,*)'s,kvecstar,sd',s,kvecstar(:,s),sd,sk,i3,j3,k3,inside
           
          !  write(*,*)s,ss1,narmsp1,sp1,q3_proc1_i(:),kvecstarp1(:,ss1)
            !xp1(:) =MATMUL(op_kmatrix(:,:,sp1),q3_proc1_i(:))
         !!!    write(*,*) 's,sp,narms',s ,sp, narms
            !if (xp1(1).eq.q3_proc1_i(1) .and. xp1(2).eq.q3_proc1_i(2) .and.
            !xp1(3).eq.q3_proc1_i(3)) then
            xp=length(kvecstar(:,s)-k)
      !     if (kvecstarp1(1,ss1).eq.q3_proc1_i(1) .and.
      !     kvecstarp1(2,ss1).eq.q3_proc1_i(2) .and.
      !     kvecstarp1(3,ss1).eq.q3_proc1_i(3)) then
          !  if (xp.lt.1d-4) then
             if (sd.eq.sk) then
              spf=sp
              ms =inside
    !!##          write(*,*)'safoura-loop1, sp1f=sp1', sp, spf
              exit loop
             endif
          !write(*,*)s,ss1,narmsp1
          enddo loop


!          call getstar(q3_proc2_i,primitivelattice,narmsp2,kvecstarp2,kvecopp2)
! loopss2: do ss2=1,narmsp2

!             sp2=kvecopp2(ss2)
             !xp2(:) =MATMUL(op_kmatrix(:,:,sp2),q3_proc2_i(:))
!             write(*,*) 's,ss2,sp2', s,ss2,sp2
             !if (xp2(1).eq.q3_proc2_i(1) .and. xp2(2).eq.q3_proc2_i(2) .and.
             !xp2(3).eq.q3_proc2_i(3)) then
!             xp2=sqrt((kvecstarp2(1,ss2)-q3_proc2(1))**2+(kvecstarp2(2,ss2)-q3_proc2(2))**2 +(kvecstarp2(3,ss2)-q3_proc2(3))**2)
             !if (kvecstarp2(1,ss2).eq.q3_proc2_i(1) .and.
             !kvecstarp2(2,ss2).eq.q3_proc2_i(2) .and.
             !kvecstarp2(3,ss2).eq.q3_proc2_i(3)) then
!             if (xp2.lt.ssigma) then
!              sp2f=sp2
 !            write(*,*)'safoura-loop2, sp2f=sp2',sp2,sp2f
!              exit loopss2
!             endif
!          enddo loopss2



end subroutine get_kvecopn
!======================================================
!==========================================================
subroutine get_kvecopn2(q,k,nk,kp,spf) !!Saf!!
!! find kvecop number for q3_proc1,2
use constants
use exactBTE2
use lattice
use params
use kpoints
use io2
implicit none
real(8) q(3), kvecstar(3,48),xp, k(3)
integer n,s,narms,sp,spf,kvecop(48),nk,sd,sk
integer i3,j3,k3, inside, ms
real(8) kp(3,nk)

      !@@    call getstar(q,primitivelattice,narms,kvecstar,kvecop)
         call kstarsinfbz1(q,nk,kp,narms,kvecstar,kvecop)
!!!***         call get_weights(nk,q,narms,kvecstar,kvecop)
          !write(*,*)s, narmsp1
         call get_k_info(k,NC,sk,i3,j3,k3,g1,g2,g3,inside)
!!!         write(*,*)'q3-k', k , i3,j3,k3,inside
 loop: do s=1,narms

             sp=kvecop(s)
          call get_k_info(kvecstar(:,s),NC,sd,i3,j3,k3,g1,g2,g3,inside)
!!!          write(*,*)'s,kvecstar,sd',s,kvecstar(:,s),sd,sk,i3,j3,k3,inside

          !  write(*,*)s,ss1,narmsp1,sp1,q3_proc1_i(:),kvecstarp1(:,ss1)
            !xp1(:) =MATMUL(op_kmatrix(:,:,sp1),q3_proc1_i(:))
         !!!    write(*,*) 's,sp,narms',s ,sp, narms
            !if (xp1(1).eq.q3_proc1_i(1) .and. xp1(2).eq.q3_proc1_i(2) .and.
            !xp1(3).eq.q3_proc1_i(3)) then
         !@@   xp=length(kvecstar(:,s)-k)
      !     if (kvecstarp1(1,ss1).eq.q3_proc1_i(1) .and.
      !     kvecstarp1(2,ss1).eq.q3_proc1_i(2) .and.
      !     kvecstarp1(3,ss1).eq.q3_proc1_i(3)) then
          !  if (xp.lt.1d-4) then
             if (sd.eq.sk) then
              spf=sp
              ms =inside
    !!##          write(*,*)'safoura-loop1, sp1f=sp1', sp, spf
              exit loop
             endif
          !write(*,*)s,ss1,narmsp1
          enddo loop


!          call getstar(q3_proc2_i,primitivelattice,narmsp2,kvecstarp2,kvecopp2)
! loopss2: do ss2=1,narmsp2

!             sp2=kvecopp2(ss2)
             !xp2(:) =MATMUL(op_kmatrix(:,:,sp2),q3_proc2_i(:))
!             write(*,*) 's,ss2,sp2', s,ss2,sp2
             !if (xp2(1).eq.q3_proc2_i(1) .and. xp2(2).eq.q3_proc2_i(2) .and.
             !xp2(3).eq.q3_proc2_i(3)) then
!             xp2=sqrt((kvecstarp2(1,ss2)-q3_proc2(1))**2+(kvecstarp2(2,ss2)-q3_proc2(2))**2
!             +(kvecstarp2(3,ss2)-q3_proc2(3))**2)
             !if (kvecstarp2(1,ss2).eq.q3_proc2_i(1) .and.
             !kvecstarp2(2,ss2).eq.q3_proc2_i(2) .and.
             !kvecstarp2(3,ss2).eq.q3_proc2_i(3)) then
!             if (xp2.lt.ssigma) then
!              sp2f=sp2
 !            write(*,*)'safoura-loop2, sp2f=sp2',sp2,sp2f
!              exit loopss2
!             endif
!          enddo loopss2
                                                                        
end subroutine get_kvecopn2
!===============================================



!==========================================================
subroutine syop_f(s,ff,fs) !!Saf!!
!unsing symmetry operations to transform a vector from IBZ to FBZ
use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none

integer s,ier
real(8) fs(3), ff(3), tempi(3,3)
 

  fs(:) = MATMUL(op_kmatrix(:,:,s),ff(:))


end subroutine syop_f
!==========================================================
!==========================================================
subroutine syop_f2(s,ff,fs2) !!Saf!!
!using inverse symetry operations to transform a vector
use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none

integer s, ier
real(8) ff(3), fs2(3), tempi(3,3)


  call xmatinv(op_kmatrix(:,:,s),tempi,ier)

  fs2(:) = MATMUL(tempi,ff(:))

end subroutine syop_f2
!=========================================================
!=====================================================
  subroutine finitedif_vvel2(q,ndn,nk,kp,s,ss)
! calculates the group velocities in units of c_light from finite difference
  use constants
  use kpoints
  implicit none
  integer ndn,i
  real(8) q0(3),vgr(3,ndn),q1(3),dq,om0,om1
  real(8) evl0(ndn),evlp(ndn),evlm(ndn)
  complex(8) evc0(ndn,ndn),evct(ndn,ndn)

  integer n,nn,ni, narms,kvecops(48),nk,s,ss
  integer i3,j3,k3, inside
  real(8) kvecstars(3,48),q(3),kp(3,nk)
  
  dq=0.001

  !call get_freq(q0,ndn,evl0,evc0)
!s=1
  call kstarsinfbz1(q,nk,kp,narms,kvecstars,kvecops)
  do n=1,narms
  s=s+1
  call syop_f(kvecops(n),q,q0)
  !q0=q00
 call get_k_info(q0,NC,ni,i3,j3,k3,g1,g2,g3,inside)
 call get_freq(q0,ndn,evl0,evc0) 
 do i=1,3
     q1=q0; q1(i)=q0(i)+dq
     call get_freq(q1,ndn,evlp,evct)
     q1=q0; q1(i)=q0(i)-dq
     call get_freq(q1,ndn,evlm,evct)
    
     vgr(i,:)=(evlp-evlm)/2/dq / 2/sqrt(evl0) *cnst*1d-8 *2*pi !*c_light
  !write(500) n,i,vgr
 enddo
!s=s+1

 do nn=1,ndn
ss=ss+1 
 write(500,*)ss,s, n,narms,ni,nn,kvecops(n),vgr(:,nn)*c_light
!ss=ss+1
 enddo
!s=s+1
enddo
!  write(30,*)'********** q=',q0,'************************'
!  do i=1,3
!     write(30,5)'i,vgr(i,l)=',i,vgr(i,:)
!  enddo
!  write(30,*)'*******************************************'

 4 format(a,3(1x,f9.3),a)
 5 format(a,i5,99(1x,f9.3))

  end subroutine finitedif_vvel2
!============================================================
!=====================================================
  subroutine RTA_vvel2(q,ndn,nk,kp,s,tempk)
! calculates the group velocities in units of c_light from finite difference
  use constants
  use kpoints
  use exactBTE2
  use params
  use lattice
  implicit none
  integer ndn,i
  real(8) q0(3),vgr(3,ndn),q1(3),dq,om0,om1
  real(8) evl0(ndn),evlp(ndn),evlm(ndn)
  complex(8) evc0(ndn,ndn),evct(ndn,ndn)

  integer n,nn,ni, narms,kvecops(48),nk,s
  integer i3,j3,k3, inside,iii
  real(8) kvecstars(3,48),q(3),kp(3,nk)
  real(8) temp,tempk,nbb,omega
 
  dq=0.001
  temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1
  !call get_freq(q0,ndn,evl0,evc0)
!s=1
  call kstarsinfbz1(q,nk,kp,narms,kvecstars,kvecops)
  do n=1,narms
  call syop_f(kvecops(n),q,q0)
  !q0=q00
 call get_k_info(q0,NC,ni,i3,j3,k3,g1,g2,g3,inside)
 call get_freq(q0,ndn,evl0,evc0)
 do i=1,3
     q1=q0; q1(i)=q0(i)+dq
     call get_freq(q1,ndn,evlp,evct)
     q1=q0; q1(i)=q0(i)-dq
     call get_freq(q1,ndn,evlm,evct)

     vgr(i,:)=(evlp-evlm)/2/dq / 2/sqrt(evl0) *cnst*1d-8 *2*pi !*c_light
  !write(500) n,i,vgr
  enddo
s=s+1
 do nn=1,ndn
!  write(500,*)s, n,narms,ni,nn,kvecops(n),vgr(:,nn)*c_light
          omega=frequency(ni,nn)
          nbb=dist(ni,nn)

          xyz_loop: do iii=1,3

                F_RTA(ni,nn,iii)= -(c_light*vgr(iii,nn)) * h_plank*(omega*c_light*100) * nbb * (nbb+1) / (k_b*tempk**2) / (Qvalue(ni,nn)*c_light*100)
          enddo xyz_loop          
 write(600,*)s, n,narms,ni,nn,kvecops(n),F_RTA(ni,nn,:) 
enddo
!s=s+1
enddo
!  write(30,*)'********** q=',q0,'************************'
!  do i=1,3
!     write(30,5)'i,vgr(i,l)=',i,vgr(i,:)
!  enddo
!  write(30,*)'*******************************************'

 4 format(a,3(1x,f9.3),a)
 5 format(a,i5,99(1x,f9.3))

  end subroutine RTA_vvel2
!============================================================



!=============================================================================
subroutine RTAs2(nk,kp,ndn,tempk,veloc,eigenval)
! calculate F value using only diagonal terms in collision matrix (equivalent as
! RTA)

use constants
use exactBTE2
use params
use lattice

implicit none

integer i,ii,iii,indx,nk,ndn
real(8) nbe,temp,tempk,nb,omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk), kp(3,nk)
!mine
integer s

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1
s=0
nkF1_loop: do i=1,nk

    laF1_loop: do ii=1,ndn

          !omega=sqrt(eigenval(ii,i)) * cnst


          !nb=nbe(omega,temp,classical)
          omega=frequency(i,ii)
          nb=dist(i,ii)

          xyz_loop: do iii=1,3

                F_RTA(i,ii,iii)= -(c_light*veloc(iii,ii,i)) * h_plank *(omega*c_light*100) * nb * (nb+1) / (k_b*tempk**2) / (Qvalue(i,ii)*c_light*100)
                ! veloc in c0, h_plank in J*s, omega in cm^-1 (linear freq.),
                ! k_b in J/K, tempk in K, Qvalue in cm^-1
                ! resulting F1 is in m/K

!write(*,*) 'indx,F1',indx,F1(indx,iii)
!                tau(i,ii,iii)= nb*(nb+1)/(Qvalue(i,ii))      ! tau is in cm

!                write(uRTA,10)
!                tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i)
          !!!     write(170,*)i,ii,iii,F_RTA(i,ii,iii)
          enddo xyz_loop
          s=s+1
          write(170,*)s,i,ii,F_RTA(i,ii,:)
    enddo laF1_loop

enddo nkF1_loop

!F1(:,:,:)=F_RTA(:,:,:) ! for iteration



end subroutine RTAs2
!============================================================

!===========================================================
subroutine cal_F2m2(nk,kp,nibbz,kibbz,ndn,tempk) !!Saf!!
! calculate F value based on F value from previous iteration.
! In this subroutine, k1, k2 are from IBZ, but by finding their stars, they are in FBZ.
! Finally, F2 are calculated in FBZ


use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none

integer nk,ndn, nibbz
integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
integer nk1,nk2,nk3, nk3_proc1, nk3_proc2,xyz
real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3), q2n(3), q3(3), qq3_proc1(3)
real(8) q3_proc1_i(3), q3_proc2_i(3)
real(8) kp(3,nk)

real(8) Avalue(nk,ndn,3)
!real(8) Avalue(nibbz,ndn,3)
real(8) kibbz(3,nibbz)   
real(8) kvecstar(3,48), kvecstarp1(3,48), kvecstarp2(3,48)
real(8) kvecstar1(3,48), q10(3)  
real(8) Fsp1(3), Fsp2(3), Fsq2(3)   
real xp1,xp2, ssigma  

integer n1,n2,n3,l1,l2,l3, nk2n
integer indx,nk_subset
integer ksub_size2, ksubset(2)
integer narms, s,sq2,ss1,ss2,sp1,sp2,sp1f,sp2f,narmsp1, narmsp2, ns, kvecop(48),kvecopp1(48),kvecopp2(48) 
integer kvecop1(48),narms1,s1
integer nk3_proc1_i, nk3_proc2_i, ns_star_i 
integer nq 
logical ins 
integer ms1, ms2


real cputim1,cputim2
real(8) cputim1_wall, cputim2_wall, wtime

real(8) temp, tempk
temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1


call cpu_time(cputim1)
cputim1_wall = wtime()
write(ulog,*) 'entering cal_F2m2...'


indx=0
Avalue=0


loop1 : do n1=1,nibbz

           q10=kibbz(:,n1)

           call kstarsinfbz1(q10,nk,kp,narms1,kvecstar1,kvecop1)
           do s1=1,narms1
              q1=kvecstar1(:,s1)
              call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
           do l1=1,ndn
             
           loopxyz : do xyz=1,3

              loop2: do n2=1,nibbz

                        q2=kibbz(:,n2)
                        call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                        call kstarsinfbz1(q2,nk,kp,narms,kvecstar,kvecop)

                     do s=1,narms

                        q3_proc1=q1+kvecstar(:,s)   ! for class1 process, coalescence process
                        q3_proc2=q1-kvecstar(:,s)   ! for class2 process, decay process
      
                        call get_k_info(kvecstar(:,s),NC,ns,i3,j3,k3,g1,g2,g3,inside)
                        call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside)
                        call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside)
            
     
                     do l2=1,ndn
                     do l3=1,ndn
           
                    Avalue(nk1,l1,xyz)=Avalue(nk1,l1,xyz) + P1(nk1,ns,l1,l2,l3)*(F1(nk3_proc1,l3,xyz)-F1(ns,l2,xyz)) + 0.5*P2(nk1,ns,l1,l2,l3)*(F1(ns,l2,xyz)+F1(nk3_proc2,l3,xyz))
                    
                     enddo
                     enddo
                  enddo  !loop s (narms)
             enddo loop2
       
                F2(nk1,l1,xyz)=F_RTA(nk1,l1,xyz) + Avalue(nk1,l1,xyz)/Qvalue(nk1,l1)

                diff(n1,l1,xyz)=(F2(n1,l1,xyz)-F1(n1,l1,xyz))/F2(n1,l1,xyz)

         enddo loopxyz
   
      enddo

!    if (n1 .eq. ksubset2(indx2,2)) then
!         indx2=indx2+1
!         call deallocate_iter2
!    endif
enddo
enddo loop1



3 format(g11.4,2(3x,3(g11.4,1x)))
!if(indx-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_F2m...'
!   stop
!endif
!   stop
!endif
4 format(2i3,99(1x,g11.4))
8 format(3i3,99(1x,g11.4))



end subroutine cal_F2m2
!==================================

!===========================================================
subroutine cal_F2m23(nk,kp,nibbz,kibbz,ndn,tempk) !!Saf!!
! calculate F value based on F value from previous iteration.
! This subroutine uses symmetry operations to calculate F in FBZ through F in IBZ.
! Summation over k1 and k2 in IBZ, but in the case of k2, F(k2 in FBZ) = S*F(k in IBZ)
! Finally, F2 are calculated in IBZ 


use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none

integer nk,ndn, nibbz
integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
integer nk1,nk2,nk3, nk3_proc1, nk3_proc2,xyz
real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3), q2n(3), q3(3), qq3_proc1(3)
real(8) q3_proc1_i(3), q3_proc2_i(3)
real(8) kp(3,nk)

real(8) Avalue(nk,ndn,3)
!real(8) Avalue(nibbz,ndn,3)
real(8) kibbz(3,nibbz)  
real(8) kvecstar(3,48), kvecstarp1(3,48), kvecstarp2(3,48)
real(8) kvecstar1(3,48), q10(3)  
real(8) Fsp1(3), Fsp2(3), Fsq2(3)   
real xp1,xp2, ssigma  

integer n1,n2,n3,l1,l2,l3, nk2n
integer indx,nk_subset
integer ksub_size2, ksubset(2)
integer narms, s,sq2,ss1,ss2,sp1,sp2,sp1f,sp2f,narmsp1, narmsp2, ns,kvecop(48),kvecopp1(48),kvecopp2(48) 
integer kvecop1(48),narms1,s1
integer nk3_proc1_i, nk3_proc2_i, ns_star_i 
integer nq,n2i0 
logical ins 
integer ms1, ms2


real cputim1,cputim2
real(8) cputim1_wall, cputim2_wall, wtime

real(8) temp, tempk
temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1


call cpu_time(cputim1)
cputim1_wall = wtime()
write(ulog,*) 'entering cal_F2m23...'


indx=0
Avalue=0

loop1 : do n1=1,nibbz

           q1=kibbz(:,n1)
 
           call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)


       do l1=1,ndn
    loopxyz : do xyz=1,3

       loop2: do n2=1,nibbz

                q2=kibbz(:,n2)
                n2i0=mapinv(n2)
                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                call kstarsinfbz1(q2,nk,kp,narms,kvecstar,kvecop)

                do s=1,narms

                   q3_proc1=q1+kvecstar(:,s)   ! for class1 process, coalescence process
                   q3_proc2=q1-kvecstar(:,s)   ! for class2 process, decay process

                   call get_k_info(kvecstar(:,s),NC,ns,i3,j3,k3,g1,g2,g3,inside)
                   call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside)
                   call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside)

                   !  find the index of q3_proc in the IBZ
                   nk3_proc1_i = mapibz(nk3_proc1)
                   nk3_proc2_i = mapibz(nk3_proc2)
                   ns_star_i = mapibz(ns)

                   ! q3_proc1,2 in IBZ
                   q3_proc1_i=kibbz(:,nk3_proc1_i)
                   q3_proc2_i=kibbz(:,nk3_proc2_i)

                   !kvecop for q2
                   sq2= kvecop(s)
                   ! find kvecop number for q3_proc1,2
                   call get_kvecopn2(q3_proc1_i,q3_proc1,nk,kp,sp1f)
                   call get_kvecopn2(q3_proc2_i,q3_proc2,nk,kp,sp2f)

 
                  do l2=1,ndn
                   do l3=1,ndn

                    !F(FBZ)=S*F(IBZ)
                    call syop_f(sq2,F1(n2i0,l2,:),Fsq2)
                    call syop_f(sp1f,F1(mapinv(nk3_proc1_i),l3,:),Fsp1)
                    call syop_f(sp2f,F1(mapinv(nk3_proc2_i),l3,:),Fsp2)

                    Avalue(mapinv(n1),l1,xyz)=Avalue(mapinv(n1),l1,xyz)+P1(mapinv(n1),ns,l1,l2,l3)*(Fsp1(xyz)-Fsq2(xyz))+0.5*P2(mapinv(n1),ns,l1,l2,l3)*(Fsq2(xyz)+Fsp2(xyz))

                  enddo
                 enddo
           enddo !loop s(narms)
        enddo loop2

              F2(mapinv(n1),l1,xyz)=F_RTA(mapinv(n1),l1,xyz) + Avalue(mapinv(n1),l1,xyz)/Qvalue(mapinv(n1),l1)


    enddo loopxyz

    enddo

!    if (n1 .eq. ksubset2(indx2,2)) then
!         indx2=indx2+1
!         call deallocate_iter2
!    endif
!@@enddo
enddo loop1



3 format(g11.4,2(3x,3(g11.4,1x)))
!if(indx-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_F2m...'
!   stop
!endif
!   stop
!endif
4 format(2i3,99(1x,g11.4))
8 format(3i3,99(1x,g11.4))

66 format(10i6,2x,99(1x,f9.3))

end subroutine cal_F2m23
!==================================      


!====================================
subroutine cal_collisionM(nk,kp,ndn,C)


use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none



integer n1,l1,n2,l2,l3,nk,ndn,s1,s2,nl3,nk3_p1,nk3_p2
integer i3,j3,k3, inside,i,j,s,kl,d1,kl1,kl2,ii,jj
integer nl2,n22
!real(8) pp,pp2
real(8) kp(3,nk),q1(3),q2(3),q3_p1(3),q3_p2(3)
real(8) C(nk*ndn*3,nk*ndn*3) 

open(301,file='c.dat')
open(302,file='c2.dat')

s1=0
s2=0
d1=0
!d2=0
!pp=0
C=0
!kl1=nk*ndn
kl=nk*ndn*3

!loopi: do i=1,kl

loopq1:  do n1=1,nk

!loopi: do i=1,kl
!          s=1
   loopl1:  do l1=1,ndn
if (s1.lt.kl) then
       s1=s1+1
!       if (s1.eq.1) then 
!         s2=1
!        else 
         s2=0
         d1=0
!        endif 
     loopq2:   do n2=1,nk
      loopl2:    do l2=1,ndn
   if (n1.eq.n2 .and. l1.eq.l2) then
!       s1=s1+1
       !s2=s2+1
         
!                 do n22=1,nk
                   do nl2=1,ndn
       loopl3:     do l3=1,ndn
                    do n22=1,nk
                             
                        C(s1,s2+1) = C(s1,s2+1) + (-P1(n1,n22,l1,nl2,l3) - 0.5*P2(n1,n22,l1,nl2,l3))
!                        s2=s1
                        write(301,*)'n1,l1,n2,l2,n22,nl2,l3,s1,s1,s2', n1,l1,n2,l2,n22,nl2,l3,s1,s1,s2,P1(n1,n22,l1,nl2,l3),P2(n1,n22,l1,nl2,l3),C(s1,s1)
!, C(s1,s2) 
                   !s2=s1
                    enddo  !loopl3
                   enddo  loopl3
                enddo
                write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s2+1,C(s1,s2+1)
                d1=s2+2
!                C(s1,d1)=C(s1,s2+1)
                 C(s1,d1)=0
                write(302,*)'n1,l1,n2,l2,s1,d1,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
                d1=s2+3
!                C(s1,d1)=C(s1,s2+1)
                 C(s1,d1)=0
                write(302,*)'n1,l1,n2,l2,s1,d1,C',n1,l1,n2,l2,s1,d1,C(s1,d1)               
!                s1=d1 
                s2=d1
                !  do i=1,3              
                 !    write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s1,C(s1,s1)
                 !  enddo
else 
    !  s1=s1+1    
  if (s2.lt.kl) then
      s2=s2+1
     q1 = kp(:,n1)
     q2 = kp(:,n2)
     q3_p1 = q1+q2
     q3_p2 = q1-q2

 !    call get_k_info(q1,NC,n01,i3,j3,k3,g1,g2,g3,inside)
 !    call get_k_info(q2,NC,n02,i3,j3,k3,g1,g2,g3,inside)

     call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
!     call get_k_info(q3_p2,NC,nk3_p2,i3,j3,k3,g1,g2,g3,inside)
     
loopnl3:  do nl3=1,ndn
   
!                C(s1,s2) = C(s1,s2) + (-P1(n1,n2,l1,l2,nl3) + P2(n1,n2,l1,l2,nl3))
   !P1(n1,nk3_p1,l1,l2,nl3) + P2(n1,n2,l1,l2,nl3))
                  C(s1,s2) = C(s1,s2) + (-P1(n1,n2,l1,l2,nl3) + P1(n1,nk3_p1,l1,l2,nl3) + P2(n1,n2,l1,l2,nl3))

               write(301,*)'n1,l1,n2,l2,nl3,s1,s2',n1,l1,n2,l2,nl3,s1,s2,P1(n1,n2,l1,l2,nl3),P1(n1,nk3_p1,l1,l2,nl3),P2(n1,n2,l1,l2,nl3),C(s1,s2)
          enddo loopnl3
            write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,s2,C(s1,s2)
            d1=s2+1
            C(s1,d1)=C(s1,s2)
            write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
            d1=s2+2
            C(s1,d1)=C(s1,s2)  
            write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
            s2=d1
!     C() = pp2
else
 cycle
endif
endif

 !  do ii=1,kl

 !      C(s1+1,ii)=C(s1,ii)
 !      write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1+1,ii,C(s1+1,ii)
 !   enddo
 !   do jj=1,kl
 !      C(s1+2,ii)=C(s1,ii)
 !      write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1+2,ii,C(s1+2,ii)
 !   enddo

!    s1=s1+2

   enddo  loopl2

   enddo  loopq2
    if ((s1+2).le.kl) then    
    do ii=1,kl
!      if ((s1+1).lt.kl) then 
       C(s1+1,ii)=C(s1,ii)
       write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,s1+1,ii,C(s1+1,ii)
 !   endif
    enddo
    do jj=1,kl
  !   if ((s1+2).lt.kl) then
       C(s1+2,jj)=C(s1,jj)
       write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,s1+2,jj,C(s1+2,jj)
   !   endif  
    enddo
   s1=s1+2
  endif
! s1=s1+2    

!enddo loopi
endif
enddo loopl1
!enddo loopi
!if (s1.lt.kl) then
! s1=s1+1
! s2=1
enddo loopq1

close(301)
close(302)

end subroutine cal_collisionM
!======================================

!=======================================
subroutine cal_collisionM2(nk,kp,ndn,C) !!Saf!!
! collision matrix (FBZ*FBZ) with considering 3Dimension

use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none



integer n1,l1,n2,l2,l3,nk,ndn,s1,s2,nl3,nk3_p1,nk3_p2
integer i3,j3,k3, inside,s,kl,kl1,kl2,ii,jj
integer nl2,n22
integer i,j,e1,e2,d1,d2
integer ij,ji
!real(8) pp,pp2
real(8) kp(3,nk),q1(3),q2(3),q3_p1(3),q3_p2(3)
real(8) C(nk*ndn*3,nk*ndn*3)

open(301,file='c_FBZ.dat')
open(302,file='c2_FBZ.dat')

s1=0
s2=0
d1=0
d2=0
e1=0
e2=0
!pp=0
C=0
!kl1=nk*ndn
kl=nk*ndn*3

!loopi: do i=1,kl

loopq1:  do n1=1,nk

!loopi: do i=1,kl
!          s=1
   loopl1:  do l1=1,ndn
if (s1.lt.kl) then
       s1=s1+1
!       if (s1.eq.1) then
!         s2=1
!        else
         s2=0
         !d1=0
!        endif
     loopq2:   do n2=1,nk
      loopl2:    do l2=1,ndn
   if (n1.eq.n2 .and. l1.eq.l2) then
!       s1=s1+1
       !s2=s2+1

!                 do n22=1,nk
                   do nl2=1,ndn
       loopl3:     do l3=1,ndn
                    do n22=1,nk
                        s2=s1
                        C(s1,s2) = C(s1,s2) + (-P1(n1,n22,l1,nl2,l3) - 0.5*P2(n1,n22,l1,nl2,l3))
!                        s2=s1
                        write(301,*)'n1,l1,n2,l2,n22,nl2,l3,s1,s1,s2', n1,l1,n2,l2,n22,nl2,l3,s1,s2,P1(n1,n22,l1,nl2,l3),P2(n1,n22,l1,nl2,l3),C(s1,s2)
!, C(s1,s2)
                   !s2=s1
                    enddo  !loopl3
                   enddo  loopl3
                enddo
               ! write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s2,C(s1,s2)
                
                e1=s1
                e2=s2
                do i=1,3
                 s1=e1+(i-1)
                 s2=e2
                  do j=1,3
                    s2=e2+(j-1)
                     if (s2.eq.s1) then
                        C(s1,s2) = C(e1,e2)
                         write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s2,C(s1,s2)
                      else 
                        C(s1,s2) = 0
                        write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s2,C(s1,s2)
                      endif
                   enddo
                 enddo
                s1=e1
!                s2=e2

              !write(*,*)'n*****', s1,s2,e1,e2
                 !d1=s2+2
!                C(s1,d1)=C(s1,s2+1)
                 !C(s1,d1)=0
                !write(302,*)'n1,l1,n2,l2,s1,d1,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
                !d1=s2+3
!                C(s1,d1)=C(s1,s2+1)
                ! C(s1,d1)=0
               ! write(302,*)'n1,l1,n2,l2,s1,d1,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
!                s1=d1
                !s2=d1
                !  do i=1,3
                 !    write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s1,C(s1,s1)
                 !  enddo
else
    !  s1=s1+1
  if (s2.lt.kl) then
      s2=s2+1
 !     s1=e1
     q1 = kp(:,n1)
     q2 = kp(:,n2)
     q3_p1 = q1+q2
     q3_p2 = q1-q2

 !    call get_k_info(q1,NC,n01,i3,j3,k3,g1,g2,g3,inside)
 !    call get_k_info(q2,NC,n02,i3,j3,k3,g1,g2,g3,inside)

     call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
!     call get_k_info(q3_p2,NC,nk3_p2,i3,j3,k3,g1,g2,g3,inside)

loopnl3:  do nl3=1,ndn

!                C(s1,s2) = C(s1,s2) + (-P1(n1,n2,l1,l2,nl3) +
!                P2(n1,n2,l1,l2,nl3))
   !P1(n1,nk3_p1,l1,l2,nl3) + P2(n1,n2,l1,l2,nl3))
                  C(s1,s2) = C(s1,s2) + (-P1(n1,n2,l1,l2,nl3) + P1(n1,nk3_p1,l1,l2,nl3) + P2(n1,n2,l1,l2,nl3))

               write(301,*)'n1,l1,n2,l2,nl3,s1,s2',n1,l1,n2,l2,nl3,s1,s2,P1(n1,n2,l1,l2,nl3),P1(n1,nk3_p1,l1,l2,nl3),P2(n1,n2,l1,l2,nl3),C(s1,s2)
          enddo loopnl3
           ! write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,s2,C(s1,s2)
         

                d1=s1
                d2=s2
              !write(*,*)d1,d2,s1,s2
                do ii=1,3
                 s1=d1+(ii-1)
                 s2=d2
                  do jj=1,3
                    s2=d2+(jj-1)
                     if ((s1.eq.d1 .and. s2.eq.d2) .or.(s1.eq.d1+1 .and. s2.eq.d2+1) .or. (s1.eq.d1+2 .and. s2.eq.d2+2)) then
                        C(s1,s2) = C(d1,d2)
                        write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s2,C(s1,s2)
                      else
                        C(s1,s2) = 0
                        write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s2,C(s1,s2)
                      endif
                   enddo
                 enddo
              s1=d1
           !  d1=s2+1
          !  C(s1,d1)=C(s1,s2)
          !  write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
           ! d1=s2+2
          !  C(s1,d1)=C(s1,s2)
          !  write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
           ! s2=d1
!     C() = pp2
else
 cycle
endif
endif

 !  do ii=1,kl

 !      C(s1+1,ii)=C(s1,ii)
 !      write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1+1,ii,C(s1+1,ii)
 !   enddo
 !   do jj=1,kl
 !      C(s1+2,ii)=C(s1,ii)
 !      write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1+2,ii,C(s1+2,ii)
 !   enddo

!    s1=s1+2

   enddo  loopl2

   enddo  loopq2
!    if ((s1+2).le.kl) then
!    do ii=1,kl
!      if ((s1+1).lt.kl) then
!       C(s1+1,ii)=C(s1,ii)
!       write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,s1+1,ii,C(s1+1,ii)
 !   endif
!    enddo
!    do jj=1,kl
  !   if ((s1+2).lt.kl) then
!       C(s1+2,jj)=C(s1,jj)
!       write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,s1+2,jj,C(s1+2,jj)
   !   endif
!    enddo
!   s1=s1+1
!  endif
! s1=s1+2

!enddo loopi
endif
s1=s1+2
enddo loopl1
!enddo loopi
!if (s1.lt.kl) then
! s1=s1+1
! s2=1
enddo loopq1

close(301)
close(302)



open(350,file='Cmatrix_FBZ.dat')

do ij=1,kl
   do ji=1,kl

  write(350,*) ij,ji, C(ij,ji)
   enddo
enddo

close(350)





end subroutine cal_collisionM2
!=====================================


!=============================================================================
subroutine cal_RHS(nk,kp,ndn,tempk,veloc,eigenval,RHSS) !!Saf!!
!(equivalent as RTA)
! Right hand side of the linearized BTE in FBZ

use constants
use exactBTE2
use params
use lattice

implicit none

integer i,ii,iii,indx,nk,ndn,s,kl,j
real(8) nbe,temp,tempk,nb,omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk), kp(3,nk), RHSS(nk*ndn*3)


temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

open(320,file='RHSS_FBZ.dat')

kl=nk*ndn*3
s=0
nkF1_loop: do i=1,nk

    laF1_loop: do ii=1,ndn

          !omega=sqrt(eigenval(ii,i)) * cnst


          !nb=nbe(omega,temp,classical)
          omega=frequency(i,ii)
          nb=dist(i,ii)

          xyz_loop: do iii=1,3
              !  do s=1,kl
                s=s+1
                RHSS(s)= -(c_light*veloc(iii,ii,i)) * h_plank * (omega*c_light*100) * nb * (nb+1) / (k_b*tempk**2) / (c_light*100)
                ! veloc in c0, h_plank in J*s, omega in cm^-1 (linear freq.),
                ! k_b in J/K, tempk in K, Qvalue in cm^-1
                ! resulting F1 is in m/K

!write(*,*) 'indx,F1',indx,F1(indx,iii)
!                tau(i,ii,iii)= nb*(nb+1)/(Qvalue(i,ii))      ! tau is in cm

!                write(uRTA,10)
!                tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i)
               !s=s+1
                write(320,*)s,i,ii,iii,RHSS(s)
              !enddo
          enddo xyz_loop

    enddo laF1_loop

enddo nkF1_loop

!F1(:,:,:)=F_RTA(:,:,:) ! for iteration

close(320)

open(402,file='RHSsy_FBZ.dat')
do j=1,kl

  write(402,*)RHSS(j)
enddo
close(402)


end subroutine cal_RHS
!===========================================================



!====================================
subroutine cal_collisionM_FBZ_x(nk,kp,ndn,C) !!Saf!!
!collision matrix (FBZ*FBZ) without considering 3dimension


use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none



integer n1,l1,n2,l2,l3,nk,ndn,s1,s2,nl3,nk3_p1,nk3_p2
integer i3,j3,k3, inside,i,s,kl
integer nl2,n22
!real(8) pp,pp2
real(8) kp(3,nk),q1(3),q2(3),q3_p1(3),q3_p2(3)
real(8) C(nk*ndn,nk*ndn), d

open(401,file='c_FBZx.dat')
open(403,file='c2_FBZx.dat')

s1=0
s2=0
!d=0
!pp=0
C=0
kl=nk*ndn

!loopi: do i=1,kl

loopq1:  do n1=1,nk

!loopi: do i=1,kl
!          s=1
   loopl1:  do l1=1,ndn
if (s1.lt.kl) then
       s1=s1+1
!       if (s1.eq.1) then
!         s2=1
!        else
         s2=0
!        endif
     loopq2:   do n2=1,nk
      loopl2:    do l2=1,ndn
   if (n1.eq.n2 .and. l1.eq.l2) then
!       s1=s1+1
       !s2=s2+1

!                 do n22=1,nk
                   do nl2=1,ndn
       loopl3:     do l3=1,ndn
                    do n22=1,nk

                        C(s1,s1) = C(s1,s1) + (-P1(n1,n22,l1,nl2,l3) - 0.5*P2(n1,n22,l1,nl2,l3))
                        s2=s1
                        write(401,*)'n1,l1,n2,l2,n22,nl2,l3,s1,s1,s2', n1,l1,n2,l2,n22,nl2,l3,s1,s1,s2,P1(n1,n22,l1,nl2,l3),P2(n1,n22,l1,nl2,l3),C(s1,s1)
!, C(s1,s2)
                   !s2=s1
                    enddo  !loopl3
                   enddo  loopl3
                enddo

                     write(403,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s1,C(s1,s1)
else
    !  s1=s1+1
  if (s2.lt.kl) then
      s2=s2+1
     q1 = kp(:,n1)
     q2 = kp(:,n2)
     q3_p1 = q1+q2
     q3_p2 = q1-q2

 !    call get_k_info(q1,NC,n01,i3,j3,k3,g1,g2,g3,inside)
 !    call get_k_info(q2,NC,n02,i3,j3,k3,g1,g2,g3,inside)

     call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
!     call get_k_info(q3_p2,NC,nk3_p2,i3,j3,k3,g1,g2,g3,inside)

loopnl3:  do nl3=1,ndn

!                C(s1,s2) = C(s1,s2) + (-P1(n1,n2,l1,l2,nl3) +
!                P2(n1,n2,l1,l2,nl3))
   !P1(n1,nk3_p1,l1,l2,nl3) + P2(n1,n2,l1,l2,nl3))
                  C(s1,s2) = C(s1,s2) + (-P1(n1,n2,l1,l2,nl3) + P1(n1,nk3_p1,l1,l2,nl3) + P2(n1,n2,l1,l2,nl3))

               write(401,*)'n1,l1,n2,l2,nl3,s1,s2',n1,l1,n2,l2,nl3,s1,s2,P1(n1,n2,l1,l2,nl3),P1(n1,nk3_p1,l1,l2,nl3),P2(n1,n2,l1,l2,nl3),C(s1,s2)
          enddo loopnl3
              write(403,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,s2,C(s1,s2)
!     C() = pp2
else
 cycle
endif
endif
   enddo  loopl2

   enddo  loopq2
!enddo loopi
endif
enddo loopl1
!enddo loopi
!if (s1.lt.kl) then
! s1=s1+1
! s2=1
enddo loopq1

close(401)
close(403)

end subroutine cal_collisionM_FBZ_x
!======================================


!=============================================================================
subroutine cal_RHS_x(nk,kp,ndn,tempk,veloc,eigenval,RHSSx) !!Saf!!
! (equivalent as RTA)
! Right hand side of the linearized BTE in just x-direction.

use constants
use exactBTE2
use params
use lattice

implicit none

integer i,ii,iii,indx,nk,ndn,s,kl,j
real(8) nbe,temp,tempk,nb,omega
real(8) veloc(3,ndn,nk), eigenval(ndn,nk), kp(3,nk), RHSSx(nk*ndn*1)


temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

open(321,file='RHSS_FBZx.dat')

kl=nk*ndn*1
s=0
nkF1_loop: do i=1,nk

    laF1_loop: do ii=1,ndn

          !omega=sqrt(eigenval(ii,i)) * cnst


          !nb=nbe(omega,temp,classical)
          omega=frequency(i,ii)
          nb=dist(i,ii)

          xyz_loop: do iii=1,1
              !  do s=1,kl
                s=s+1
                RHSSx(s)= -(c_light*veloc(iii,ii,i)) * h_plank * (omega*c_light*100) * nb * (nb+1) / (k_b*tempk**2)/ (c_light*100)

                ! veloc in c0, h_plank in J*s, omega in cm^-1 (linear freq.),
                ! k_b in J/K, tempk in K, Qvalue in cm^-1
                ! resulting F1 is in m/K

!write(*,*) 'indx,F1',indx,F1(indx,iii)
!                tau(i,ii,iii)= nb*(nb+1)/(Qvalue(i,ii))      ! tau is in cm

!                write(uRTA,10)
!                tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i)
               !s=s+1
                write(321,*)s,i,ii,iii,RHSSx(s)
              !enddo
          enddo xyz_loop

    enddo laF1_loop

enddo nkF1_loop

!F1(:,:,:)=F_RTA(:,:,:) ! for iteration

close(321)

open(404,file='RHSsy_FBZx.dat')
do j=1,kl

  write(404,*)RHSSx(j)
enddo
close(404)




end subroutine cal_RHS_x
!===========================================================

!===========================================================
subroutine cal_collisionM2_IBZ(nk,kp,nibbz,kibbz,ndn,C)


use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none



integer n1,l1,n2,l2,l3,nk,ndn,s1,s2,nl3,nk3_p1,nk3_p2,nibbz
integer i3,j3,k3, inside,s,kl,kl1,kl2,ii,jj,kil
integer nl2,n22
integer i,j,e1,e2,d1,d2
integer ij,ji
integer n1i
!real(8) pp,pp2
real(8) kp(3,nk),kibbz(3,nibbz),q1(3),q2(3),q3_p1(3),q3_p2(3)
real(8) C(nibbz*ndn*3,nk*ndn*3)

open(331,file='c_IBZ.dat')
open(332,file='c2_IBZ.dat')

s1=0
s2=0
d1=0
d2=0
e1=0
e2=0
!pp=0
C=0
!kl1=nk*ndn
kl=nk*ndn*3
kil=nibbz*ndn*3
!loopi: do i=1,kl

loopq1:  do n1=1,nibbz
            n1i=mapinv(n1)
!loopi: do i=1,kl
!          s=1
   loopl1:  do l1=1,ndn
if (s1.lt.kil) then
       s1=s1+1
!       if (s1.eq.1) then
!         s2=1
!        else
         s2=0
         !d1=0
!        endif
     loopq2:   do n2=1,nk
      loopl2:    do l2=1,ndn
   if (n1i.eq.n2 .and. l1.eq.l2) then
!       s1=s1+1
       !s2=s2+1

!                 do n22=1,nk
                   do nl2=1,ndn
       loopl3:     do l3=1,ndn
                    do n22=1,nk
                        s2=s1
                        C(s1,s2) = C(s1,s2) + (-P1(n1i,n22,l1,nl2,l3) - 0.5*P2(n1i,n22,l1,nl2,l3))
!                        s2=s1
                        write(331,*)'n1,l1,n2,l2,n22,nl2,l3,s1,s1,s2', n1,l1,n2,l2,n22,nl2,l3,s1,s2,P1(n1,n22,l1,nl2,l3),P2(n1,n22,l1,nl2,l3),C(s1,s2)
!, C(s1,s2)
                   !s2=s1
                    enddo  !loopl3
                   enddo  loopl3
                enddo
               ! write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s2,C(s1,s2)

                e1=s1
                e2=s2
                do i=1,3
                 s1=e1+(i-1)
                 s2=e2
                  do j=1,3
                    s2=e2+(j-1)
                     if (s2.eq.s1) then
                        C(s1,s2) = C(e1,e2)
                         write(332,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                      else
                        C(s1,s2) = 0
                        write(332,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                      endif
                   enddo
                 enddo
                s1=e1
!                s2=e2

              !write(*,*)'n*****', s1,s2,e1,e2
                 !d1=s2+2
!                C(s1,d1)=C(s1,s2+1)
                 !C(s1,d1)=0
                !write(302,*)'n1,l1,n2,l2,s1,d1,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
                !d1=s2+3
!                C(s1,d1)=C(s1,s2+1)
                ! C(s1,d1)=0
               ! write(302,*)'n1,l1,n2,l2,s1,d1,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
!                s1=d1
                !s2=d1
                !  do i=1,3
                 !    write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s1,C(s1,s1)
                 !  enddo
else
    !  s1=s1+1
  if (s2.lt.kl) then
      s2=s2+1
 !     s1=e1
     q1 = kibbz(:,n1)
     q2 = kp(:,n2)
     q3_p1 = q1+q2
     q3_p2 = q1-q2

 !    call get_k_info(q1,NC,n01,i3,j3,k3,g1,g2,g3,inside)
 !    call get_k_info(q2,NC,n02,i3,j3,k3,g1,g2,g3,inside)

     call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
!     call get_k_info(q3_p2,NC,nk3_p2,i3,j3,k3,g1,g2,g3,inside)

loopnl3:  do nl3=1,ndn

!                C(s1,s2) = C(s1,s2) + (-P1(n1,n2,l1,l2,nl3) +
!                P2(n1,n2,l1,l2,nl3))
   !P1(n1,nk3_p1,l1,l2,nl3) + P2(n1,n2,l1,l2,nl3))
                  C(s1,s2) = C(s1,s2) + (-P1(n1i,n2,l1,l2,nl3) + P1(n1i,nk3_p1,l1,l2,nl3) + P2(n1i,n2,l1,l2,nl3))

               write(331,*)'n1,l1,n2,l2,nl3,s1,s2',n1,l1,n2,l2,nl3,s1,s2,P1(n1,n2,l1,l2,nl3),P1(n1,nk3_p1,l1,l2,nl3),P2(n1,n2,l1,l2,nl3),C(s1,s2)
          enddo loopnl3
           ! write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,s2,C(s1,s2)


                d1=s1
                d2=s2
              !write(*,*)d1,d2,s1,s2
                do ii=1,3
                 s1=d1+(ii-1)
                 s2=d2
                  do jj=1,3
                    s2=d2+(jj-1)
                     if ((s1.eq.d1 .and. s2.eq.d2) .or.(s1.eq.d1+1 .and. s2.eq.d2+1) .or. (s1.eq.d1+2 .and. s2.eq.d2+2)) then
                        C(s1,s2) = C(d1,d2)
                        write(332,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                      else
                        C(s1,s2) = 0
                        write(332,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                      endif
                   enddo
                 enddo
              s1=d1
           !  d1=s2+1
          !  C(s1,d1)=C(s1,s2)
          !  write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
           ! d1=s2+2
          !  C(s1,d1)=C(s1,s2)
          !  write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
           ! s2=d1
!     C() = pp2
else
 cycle
endif
endif

 !  do ii=1,kl

 !      C(s1+1,ii)=C(s1,ii)
 !      write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1+1,ii,C(s1+1,ii)
 !   enddo
 !   do jj=1,kl
 !      C(s1+2,ii)=C(s1,ii)
 !      write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1+2,ii,C(s1+2,ii)
 !   enddo

!    s1=s1+2

   enddo  loopl2

   enddo  loopq2
!    if ((s1+2).le.kl) then
!    do ii=1,kl
!      if ((s1+1).lt.kl) then
!       C(s1+1,ii)=C(s1,ii)
!       write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,s1+1,ii,C(s1+1,ii)
 !   endif
!    enddo
!    do jj=1,kl
  !   if ((s1+2).lt.kl) then
!       C(s1+2,jj)=C(s1,jj)
!       write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,s1+2,jj,C(s1+2,jj)
   !   endif
!    enddo
!   s1=s1+1
!  endif
! s1=s1+2

!enddo loopi
endif
s1=s1+2
enddo loopl1
!enddo loopi
!if (s1.lt.kl) then
! s1=s1+1
! s2=1
enddo loopq1

close(331)
close(332)



open(351,file='Cmatrix_IBZ.dat')

do ij=1,kil
   do ji=1,kl

  write(351,*) ij,ji, C(ij,ji)
   enddo
enddo

close(351)





end subroutine cal_collisionM2_IBZ
!=====================================


!===========================================================
subroutine cal_collisionM2_IBZ2(nk,kp,nibbz,kibbz,ndn,C) !!Saf!!
!collision matrix (IBZ*FBZ)

use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none



integer n1,l1,n2,l2,l3,nk,ndn,s1,s2,nl3,nk3_p1,nk3_p2,nibbz
integer i3,j3,k3, inside,s,kl,kl1,kl2,ii,jj,kil
integer nl2,n22
integer i,j,e1,e2,d1,d2
integer ij,ji
integer n1i
!real(8) pp,pp2
real(8) kp(3,nk),kibbz(3,nibbz),q1(3),q2(3),q3_p1(3),q3_p2(3)
real(8) C(nibbz*ndn*3,nk*ndn*3)

open(331,file='c_IBZ.dat')
open(332,file='c2_IBZ.dat')

s1=0
s2=0
d1=0
d2=0
e1=0
e2=0
!pp=0
C=0
!kl1=nk*ndn
kl=nk*ndn*3
kil=nibbz*ndn*3
!loopi: do i=1,kl

loopq1:  do n1=1,nibbz
            n1i=mapinv(n1)
!loopi: do i=1,kl
!          s=1
   loopl1:  do l1=1,ndn
if (s1.lt.kil) then
       s1=s1+1
!       if (s1.eq.1) then
!         s2=1
!        else
         s2=0
         !d1=0
!        endif
     loopq2:   do n2=1,nk
      loopl2:    do l2=1,ndn
   if (n1i.eq.n2 .and. l1.eq.l2) then
!       s1=s1+1
       s2=s2+1
                
!                 do n22=1,nk
                   do nl2=1,ndn
       loopl3:     do l3=1,ndn
                    do n22=1,nk
                       ! s2=s2+1
                        C(s1,s2) = C(s1,s2) + (-P1(n1i,n22,l1,nl2,l3) - 0.5*P2(n1i,n22,l1,nl2,l3))
!                        s2=s1
                        write(331,*)'n1,l1,n2,l2,n22,nl2,l3,s1,s1,s2',n1,l1,n2,l2,n22,nl2,l3,s1,s2,P1(n1,n22,l1,nl2,l3),P2(n1,n22,l1,nl2,l3),C(s1,s2)
!, C(s1,s2)
                   !s2=s1
                    enddo  !loopl3
                   enddo  loopl3
                enddo
               ! write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s2,C(s1,s2)

                e1=s1
                e2=s2
                do i=1,3
                 s1=e1+(i-1)
                 s2=e2
                  do j=1,3
                    s2=e2+(j-1)
       !              if (s2.eq.s1) then
                      if ((s1.eq.e1 .and. s2.eq.e2) .or.(s1.eq.e1+1 .and. s2.eq.e2+1) .or. (s1.eq.e1+2 .and. s2.eq.e2+2)) then
                        C(s1,s2) = C(e1,e2)
                         write(332,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                      !   write(*,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                      else
                        C(s1,s2) = 0
                        write(332,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                      !  write(*,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                      endif
                   enddo
                 enddo
                s1=e1
!                s2=e2

              !write(*,*)'n*****', s1,s2,e1,e2
                 !d1=s2+2
!                C(s1,d1)=C(s1,s2+1)
                 !C(s1,d1)=0
                !write(302,*)'n1,l1,n2,l2,s1,d1,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
                !d1=s2+3
!                C(s1,d1)=C(s1,s2+1)
                ! C(s1,d1)=0
               ! write(302,*)'n1,l1,n2,l2,s1,d1,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
!                s1=d1
                !s2=d1
                !  do i=1,3
                 !    write(302,*)'n1,l1,n2,l2,s1,s1,C',n1,l1,n2,l2,s1,s1,C(s1,s1)
                 !  enddo
else
    !  s1=s1+1
  if (s2.lt.kl) then
      s2=s2+1
 !     s1=e1
     q1 = kibbz(:,n1)
     q2 = kp(:,n2)
     q3_p1 = q1+q2
     q3_p2 = q1-q2

 !    call get_k_info(q1,NC,n01,i3,j3,k3,g1,g2,g3,inside)
 !    call get_k_info(q2,NC,n02,i3,j3,k3,g1,g2,g3,inside)

     call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
!     call get_k_info(q3_p2,NC,nk3_p2,i3,j3,k3,g1,g2,g3,inside)

loopnl3:  do nl3=1,ndn

!                C(s1,s2) = C(s1,s2) + (-P1(n1,n2,l1,l2,nl3) +
!                P2(n1,n2,l1,l2,nl3))
   !P1(n1,nk3_p1,l1,l2,nl3) + P2(n1,n2,l1,l2,nl3))
                  C(s1,s2) = C(s1,s2) + (-P1(n1i,n2,l1,l2,nl3) + P1(n1i,nk3_p1,l1,l2,nl3) + P2(n1i,n2,l1,l2,nl3))

               write(331,*)'n1,l1,n2,l2,nl3,s1,s2',n1,l1,n2,l2,nl3,s1,s2,P1(n1,n2,l1,l2,nl3),P1(n1,nk3_p1,l1,l2,nl3),P2(n1,n2,l1,l2,nl3),C(s1,s2)
          enddo loopnl3
           ! write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,s2,C(s1,s2)


                d1=s1
                d2=s2
              !write(*,*)d1,d2,s1,s2
                do ii=1,3
                 s1=d1+(ii-1)
                 s2=d2
                  do jj=1,3
                    s2=d2+(jj-1)
                     if ((s1.eq.d1 .and. s2.eq.d2) .or.(s1.eq.d1+1 .and. s2.eq.d2+1) .or. (s1.eq.d1+2 .and. s2.eq.d2+2)) then
                        C(s1,s2) = C(d1,d2)
                        write(332,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                   !     write(*,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                      else
                        C(s1,s2) = 0
                        write(332,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                  !      write(*,*)'n1,n1i,l1,n2,l2,s1,s1,C',n1,n1i,l1,n2,l2,s1,s2,C(s1,s2)
                      endif
                   enddo
                 enddo
              s1=d1
           !  d1=s2+1
          !  C(s1,d1)=C(s1,s2)
          !  write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
           ! d1=s2+2
          !  C(s1,d1)=C(s1,s2)
          !  write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1,d1,C(s1,d1)
           ! s2=d1
!     C() = pp2
else
 cycle
endif
endif

 !  do ii=1,kl

 !      C(s1+1,ii)=C(s1,ii)
 !      write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1+1,ii,C(s1+1,ii)
 !   enddo
 !   do jj=1,kl
 !      C(s1+2,ii)=C(s1,ii)
 !      write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,n2,l2,s1+2,ii,C(s1+2,ii)
 !   enddo

!    s1=s1+2

   enddo  loopl2

   enddo  loopq2
!    if ((s1+2).le.kl) then
!    do ii=1,kl
!      if ((s1+1).lt.kl) then
!       C(s1+1,ii)=C(s1,ii)
!       write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,s1+1,ii,C(s1+1,ii)
 !   endif
!    enddo
!    do jj=1,kl
  !   if ((s1+2).lt.kl) then
!       C(s1+2,jj)=C(s1,jj)
!       write(302,*)'n1,l1,n2,l2,s1,s2,C',n1,l1,s1+2,jj,C(s1+2,jj)
   !   endif
!    enddo
!   s1=s1+1
!  endif
! s1=s1+2

!enddo loopi
endif
s1=s1+2
enddo loopl1
!enddo loopi
!if (s1.lt.kl) then
! s1=s1+1
! s2=1
enddo loopq1

close(331)
close(332)



open(351,file='Cmatrix_IBZ.dat')

do ij=1,kil
   do ji=1,kl

  write(351,*) ij,ji, C(ij,ji)
   enddo
enddo

close(351)





end subroutine cal_collisionM2_IBZ2
!=====================================


!=============================================================================
subroutine cal_RHS_IBZ(nibbz,kibbz,ndn,tempk,veloc,eigenval,RHSS) !!Saf!!
! (equivalent as RTA) 
! Right hand side of the linearized BTE in just IBZ.

use constants
use exactBTE2
use params
use lattice
use kpoints 

implicit none

integer i,ii,iii,indx,ndn,s,kl,j,nibbz,n0
real(8) nbe,temp,tempk,nbb,omega
real(8) veloc(3,ndn,nibbz), eigenval(ndn,nibbz), RHSS(nibbz*ndn*3), kibbz(3,nibbz)


temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

open(322,file='RHSS_IBZ.dat')

kl=nibbz*ndn*3
s=0
nkF1_loop: do i=1,nibbz

          n0 = mapinv(i)

    laF1_loop: do ii=1,ndn

          !omega=sqrt(eigenval(ii,i)) * cnst


          !nb=nbe(omega,temp,classical)
          omega=frequency(n0,ii)
          nbb=dist(n0,ii)

          xyz_loop: do iii=1,3
              !  do s=1,kl
                s=s+1
                !RHSS(s)= -(c_light*veloc(iii,ii,n0)) * h_plank * (omega*c_light*100) * nbb * (nbb+1) / (k_b*tempk**2) / (c_light*100)
                RHSS(s)= (c_light*veloc(iii,ii,n0)) * h_plank * (omega*c_light*100) * nbb * (nbb+1) / (k_b*tempk**2) / (c_light*100)
                ! veloc in c0, h_plank in J*s, omega in cm^-1 (linear freq.),
                ! k_b in J/K, tempk in K, Qvalue in cm^-1
                ! resulting F1 is in m/K

!write(*,*) 'indx,F1',indx,F1(indx,iii)
!                tau(i,ii,iii)= nb*(nb+1)/(Qvalue(i,ii))      ! tau is in cm

!                write(uRTA,10)
!                tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i)
               !s=s+1
                write(322,*)s,i,n0,ii,iii,RHSS(s)
              !enddo
          enddo xyz_loop

    enddo laF1_loop

enddo nkF1_loop

!F1(:,:,:)=F_RTA(:,:,:) ! for iteration

close(322)

open(412,file='RHS_IBZ.dat')
do j=1,kl

  write(412,*)RHSS(j)
enddo
close(412)




end subroutine cal_RHS_IBZ
!===========================================================


!==========================================================
subroutine sy_matrix_IBZ_FBZ(nibbz,kibbz,nk,kp,ndn,S)


use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none

integer nk,ndn, nibbz,narms, kvecops(48)
integer i,j,i3,j3,k3
integer n,ni0,na,n0,l,s1,s2,d,kil,kl,d1,d2
integer iii,r,jj
real(8) q(3),kvecstars(3,48)
real(8) kp(3,nk),kibbz(3,nibbz), fiRTA(nibbz*ndn*3), S(nk*ndn*3,nibbz*ndn*3)
real(8) f_fbz(nk*ndn*3),fsRTA(3)

open(440,file='SS.dat')
open(441,file='firta.dat')
open(443,file='f_fsRTA.dat')
s1=0
s2=1
S=0
d1=0
d2=0
kil=nibbz*ndn*3
kl=nk*ndn*3
r=0
do n=1,nibbz
   q=kibbz(:,n)
   ni0=mapinv(n)
  ! i=i+1
   call kstarsinfbz1(q,nk,kp,narms,kvecstars,kvecops)

  do l=1,ndn
    do iii=1,3
      r=r+1
      fiRTA(r)=F_RTA(ni0,l,iii)
      write(441,*) n,ni0,l,r,iii,fiRTA(r)
  !  r=r+1
    enddo
    do na=1,narms
   !    call get_k_info(kvecstars(:,na),NC,n0,i3,j3,k3,g1,g2,g3,inside)
     !  nz=mapibz(n0)
      call syop_f(kvecops(na),F_RTA(ni0,l,:),fsRTA)
    
      write(443,*) n,ni0,l,na,narms,kvecops(na),F_RTA(ni0,l,:),fsRTA
!      write(443,*)op_kmatrix(:,:,kvecops(na))
     
      if (s1.le.(kl-3)) then
       s1=s1+1
       d1=s1
       d2=s2
       do i=1,3
        s1=s1+(i-1)
       ! s2=d2
          do j=1,3
             s2=s2+(j-1)
             S(s1,s2) = op_kmatrix(i,j,kvecops(na))      
             write(440,*) n,l,na,narms,d,s1,s2,S(s1,s2)
            s2=d2
           enddo
         s1=d1
       enddo
        s1=s1+2

 !      S(s1,s2+1) = op_kmatrix(1,j,na) 
 !      S(s1,s2+2) = op_kmatrix(1,j,na)
       
!       S(s1+1,s2) = op_kmatrix(2,1,na)
!       S(s1+1,s2+1) = op_kmatrix(2,2,na)
!       S(s1+1,s2+2) = op_kmatrix(2,3,na)
       
!       S(s1+2,s2) = op_kmatrix(3,1,na)
!       S(s1+2,s2+1) = op_kmatrix(3,2,na)
!       S(s1+2,s2+2) = op_kmatrix(3,3,na)
       else
       exit
      endif

      enddo
  !  s2=s2+3
  if (s2.le.(kil-3)) then
   s2=s2+3
   else
   cycle
   endif
  ! r=r+1
   enddo
   !s2=s2+3
enddo


close(440)
close(441)
close(443)

open(442,file='f_fbz.dat')

f_fbz = MATMUL(S,fiRTA)
 r=0
do n=1,nk
  do l=1,ndn
    do iii=1,3
      r=r+1
!do jj=1,kl

write(442,4)r,f_fbz(r),F_RTA(n,l,iii) ,f_fbz(r)-F_RTA(n,l,iii) 

enddo
enddo
enddo

close(442)
4 format(i5,3(1x,g12.5))

end subroutine sy_matrix_IBZ_FBZ 
!==========================================================






!==========================================================
subroutine sy_matrix_IBZ_FBZ2(nibbz,kibbz,nk,kp,ndn,S) !!Saf!!
!Symmetry matrix(FBZ*IBZ)

use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none

integer nk,ndn, nibbz,narms, kvecops(48)
integer i,j,i3,j3,k3,inside
integer n,ni0,na,n0,l,s1,s2,d,kil,kl,d1,d2
integer iii,r,jj,ll,ro,col
real(8) q(3),kvecstars(3,48)
real(8) kp(3,nk),kibbz(3,nibbz), fiRTA(nibbz*ndn*3), S(nk*ndn*3,nibbz*ndn*3)
real(8) f_fbz(nk*ndn*3),fsRTA(3)

open(440,file='SS.dat')
open(441,file='firta.dat')
open(443,file='f_fsRTA.dat')
s1=0
s2=1
S=0
d1=0
d2=0
kil=nibbz*ndn*3
kl=nk*ndn*3
r=0
do n=1,nibbz
   q=kibbz(:,n)
   ni0=mapinv(n)
  ! i=i+1
   call kstarsinfbz1(q,nk,kp,narms,kvecstars,kvecops)

!!  do l=1,ndn
  !  do iii=1,3
  !    r=r+1
  !    fiRTA(r)=F_RTA(ni0,l,iii)
  !    write(441,*) n,ni0,l,r,iii,fiRTA(r)
  !  r=r+1
  !  enddo
  !enddo
   write(440,*)'narms',narms
    do na=1,narms
       call get_k_info(kvecstars(:,na),NC,n0,i3,j3,k3,g1,g2,g3,inside)
     !  nz=mapibz(n0)
!      call syop_f(kvecops(na),F_RTA(ni0,l,:),fsRTA)

!      write(443,*) n,ni0,l,na,narms,kvecops(na),F_RTA(ni0,l,:),fsRTA
!      write(443,*)op_kmatrix(:,:,kvecops(na))

   !   if (s1.le.(kl-3)) then
   !    s1=s1+1
   
       s1=((n0-1)*3*ndn)+1
       s2=((n-1)*3*ndn)+1
       write(*,*)'n, na, narms,ni0,n0,s1,s2', n, na, narms,ni0,n0,s1,s2
       !d1=s1
       !d2=s2
     do ll=1,ndn
       d2=s2
        do i=1,3
        !s1=s1+(i-1)
       ! s2=d2
          do j=1,3
             s2=d2+(j-1)
             S(s1,s2) = op_kmatrix(i,j,kvecops(na))
             write(440,*) n,ll,na,narms,ni0,n0,s1,s2,S(s1,s2)
           ! s2=d2
        !     s2=s2+1
           enddo
!         if (s1.lt.kl) then
         s1=s1+1
!         else
 !        exit
  !       endif
       enddo
       !s1=s1+3
!       if (s2.lt.kil) then
       s2=s2+1
!       else
!       cycle
!       endif
       enddo
 !      S(s1,s2+1) = op_kmatrix(1,j,na)
 !      S(s1,s2+2) = op_kmatrix(1,j,na)

!       S(s1+1,s2) = op_kmatrix(2,1,na)
!       S(s1+1,s2+1) = op_kmatrix(2,2,na)
!       S(s1+1,s2+2) = op_kmatrix(2,3,na)

!       S(s1+2,s2) = op_kmatrix(3,1,na)
!       S(s1+2,s2+1) = op_kmatrix(3,2,na)
!       S(s1+2,s2+2) = op_kmatrix(3,3,na)
!       else
!       exit
!      endif

      enddo
  !  s2=s2+3
 ! if (s2.le.(kil-3)) then
 !  s2=s2+3
 !  else
 !  cycle
 !  endif
  ! r=r+1
  ! enddo
   !s2=s2+3
enddo


close(440)
!close(441)
close(443)

open(445,file='SS2.dat')
  do ro=1,kl
    do col=1,kil
       write(445,*)ro,col, S(ro,col) 
    enddo
  enddo
close(445)


do n=1,nibbz
!   q=kibbz(:,n)
   ni0=mapinv(n)

   do l=1,ndn
     do iii=1,3
       r=r+1
       fiRTA(r)=F_RTA(ni0,l,iii)
       write(441,*) n,ni0,l,r,iii,fiRTA(r)
   !  r=r+1
     enddo
   enddo

enddo

close(441)

open(442,file='f_fbz2.dat')

f_fbz = MATMUL(S,fiRTA)
 r=0
do n=1,nk
  do l=1,ndn
    do iii=1,3
      r=r+1
!do jj=1,kl

!write(442,4)r,f_fbz(r),n,l,iii,F_RTA(n,l,iii) ,f_fbz(r)-F_RTA(n,l,iii)
write(442,*) r,f_fbz(r),n,l,iii,F_RTA(n,l,iii) ,f_fbz(r)-F_RTA(n,l,iii)
enddo
enddo
enddo

close(442)
4 format(i5,3(1x,g12.5))

end subroutine sy_matrix_IBZ_FBZ2
!==========================================================



!==========================================================
subroutine sycollisionM(nibbz,kibbz,ndn,nk,kp,Cs) !!Saf!!
!Collision matrix(IBZ*FBZ)*Symmerty matrix(FBZ*IBZ)=Cs(IBZ*IBZ) :output

use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none

integer n1,l1,n2,l2,ns, s1,s2
integer ndn,nibbz,narms,nk
integer i,j,kil,ii,jj,d
real(8) kp(3,nk),kibbz(3,nibbz)
real(8) q(3),q1(3),q2(3),q3_p1(3),q3_p2(3)
integer i3,j3,k3,inside,nk3_p1
real(8) Cs(nibbz*ndn*3,nibbz*ndn*3), C_d, C_of, C
real(8) kvecstars(3,48)
integer kvecops(48),n0,n22,nl2,l3,nl3,ni0,n1i,d1,d2
s1=0
s2=0
Cs=0
C=0
kil=nibbz*ndn*3

open(701,file='sycmatrix.dat')

    do n1=1,nibbz    
       n1i=mapinv(n1)
    do l1=1,ndn
      s1=s1+1
      d1=s1

      do n2=1,nibbz
           q=kibbz(:,n2)
           ni0=mapinv(n2)
           call kstarsinfbz1(q,nk,kp,narms,kvecstars,kvecops)
      write(701,*) '**************',n2,narms    
       do l2=1,ndn  
          s2=s2+1
          d2=s2
         !C_d=0
         !C_of=0
       write(701,*)'n1,  n1i,  l1,  n2,  ni0,  n0,  l2,  narms,  ns,  s1,  s2,  C_d,  C_of,   Cs(s1,s2)'
       do ns=1,narms
          call get_k_info(kvecstars(:,ns),NC,n0,i3,j3,k3,g1,g2,g3,inside)
      ! do l2=1,ndn
     !    C_d=0
     !    C_of=0
        ! C2=0
        ! C3=0
        ! C4=0
        ! C5=0
           if (n1i.eq.n0 .and. l1.eq.l2) then
               C=0
                 do nl2=1,ndn
                   do l3=1,ndn
                    do n22=1,nk
                        
                        C = C + (-P1(n1i,n22,l1,nl2,l3) - 0.5*P2(n1i,n22,l1,nl2,l3))

                    enddo  
                   enddo  
                enddo
            ! C2 = C_d
            ! C3 = C_d
    
           else
           
         q1 = kibbz(:,n1)
         q2 = kvecstars(:,ns)
         q3_p1 = q1+q2
         q3_p2 = q1-q2
   
        ! write(*,*)'here', n1,n1i,l1,n2,ns,narms,n0,l2    
  
         call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
               C=0
      loopnl3:  do nl3=1,ndn

                  C = C + (-P1(n1i,n0,l1,l2,nl3) + P1(n1i,nk3_p1,l1,l2,nl3) + P2(n1i,n0,l1,l2,nl3))
                enddo loopnl3
         !    C4 = C_of
         !    C5 = C_of
           endif
         s2=d2
         s1=d1
         do i=1,3
           s2=d2
           s1=d1+(i-1)
            do j=1,3
                s2=d2+(j-1)
                ! i=1
                ! j=1
                Cs(s1,s2) = Cs(s1,s2) + (C * op_kmatrix(i,j,kvecops(ns)))
                !write(*,*)'n1,n1i,l1,n2,l2,ni0,narms,ns,n0,s1,s2, C'
                write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,C_d,C_of,op_kmatrix(i,j,kvecops(ns)),kvecops(ns), Cs(s1,s2)
                !s2=s2+1
            enddo
           !s1=s1+1
         enddo 
        !s1=d1
        !if (s2.lt.kil) then
        !s2=s2+1
        !else
        !cycle
        !endif 
       enddo   !narms loop
!      write(*,*)'after summation over narms'
!       d2=s2
!       d1=s1
!         do i=1,3
!           do j=1,3
!             s2=d2+(j-1)
!             Cs(s1,s2) = C * op_kmatrix(i,j,kvecops(ns))
!             write(*,*) s1,s2,C,ns, op_kmatrix(1,1,kvecops(1)), Cs(s1,s2),i,j
!           enddo
!          s1=s1+1
!         enddo

!        s1=d1
!        if (s2.lt.kil) then
!        s2=s2+1
 !       else
 !       cycle
 !       endif

       enddo    !l2 loop
!      s2=d2
!      s1=s1+1
      enddo     !n2 loop

!   if (s1.le.(kil-3)) then
!    s1=s1+3  
!    else
!    exit
!    endif
   s2=0
   !s1=s1+1 
  enddo
   enddo

close(701)

open(601,file='Cs_sycollisionm.dat')
do ii=1,kil
   do jj=1,kil
      write(601,*) ii,jj, Cs(ii,jj)

   enddo
enddo
close(601)



end subroutine sycollisionM
!===========================================================


!==========================================================
subroutine veloc_avergeIBZ_invsyFBZ(nk,kp,nibbz,kibbz,ndn,veloc,veloci) !!Saf!!

use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none


integer i,j,ii,n,jj
integer nibbz,nk,ni0,n0,narms,ndn
integer i3,j3,k3,inside, kvecops(48)
real(8) kp(3,nk),kibbz(3,nibbz),qi(3)
real(8) veli(3),veli_a(3),veloci(3,ndn,nibbz)
real(8) kvecstars(3,48)
real(8) veloc(3,ndn,nk)
real(8) vel0(3)
veli_a=0
veli=0
veloci=0
open(790,file='velIBZ_aver.dat')
do i=1,nibbz
   qi=kibbz(:,i)
   ni0=mapinv(i)
!  call kstarsinfbz1(qi,nk,kp,narms,kvecstars,kvecops)
   call getstar(qi,primitivelattice,narms,kvecstars,kvecops)
  do j=1,ndn
  veli_a=0
  veli=0
  do n=1,narms
   call get_k_info(kvecstars(:,n),NC,n0,i3,j3,k3,g1,g2,g3,inside)
   call syop_f2(kvecops(n),veloc(:,j,n0),veli)
  ! write(*,*)'syop_f2'
   call syop_f(kvecops(n),veloc(:,j,ni0),vel0)
 write(790,*)i,j,n,narms,veloc(:,j,n0)*c_light,kvecops(n), veli*c_light
 write(790,*)i,j,n,narms,vel0*c_light,kvecops(n),veloc(:,j,ni0)*c_light
 write(790,*)'********************************************************************'
   do ii=1,3
      veli_a(ii) = veli_a(ii) + (veli(ii)/narms)
   enddo

enddo
!   do jj=1,3
    veloci(:,j,i) = veli_a(:)
!   enddo
  write(790,*)i,j, veloci(:,j,i)

enddo
enddo

close(790)


end subroutine veloc_avergeIBZ_invsyFBZ
!==========================================================













!===================================
subroutine kstarsinfbz(nk,kp,nibbz,kibbz,ds,kvecstars,kvecops) !!Saf!!


use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none

integer nk, nibbz,narms,kvecop(48), kvecops(48)
integer i,ds,s,j,sl
real(8) q(3),qk(3),kvecstar(3,48),kvecstars(3,48)
real(8) kibbz(3,nibbz), kp(3,nk)
sl=0
do i=1,nibbz
ds=0    
    q=kibbz(:,i)
    call getstar(q,primitivelattice,narms,kvecstar,kvecop)
    do s=1,narms
       do j=1,nk
          qk=kp(:,j)
          if (qk .myeq. kvecstar(:,s)) then    
             ds=ds+1
             kvecstars(:,ds)=kvecstar(:,s)
             kvecops(ds)=kvecop(s)
             sl=sl+1
!***            write(190,*)sl,i,ds,kvecops(ds),kvecstars(:,ds)
          endif
        enddo
     enddo
enddo



end subroutine kstarsinfbz
!=================================== subroutine above just for one q from ibz
subroutine kstarsinfbz1(q,nk,kp,ds,kvecstars,kvecops) !!Saf!!


use constants
use exactBTE2
use params
use lattice
use io2
use mod_ksubset2
use phi3
use kpoints
use force_constants_module
implicit none

integer nk, nibbz,narms,kvecop(48), kvecops(48)
integer i,ds,s,j,sl
real(8) q(3),qk(3),kvecstar(3,48),kvecstars(3,48)
real(8) kp(3,nk)
!sl=0
!do i=1,nibbz
ds=0
!    q=kibbz(:,i)
    call getstar(q,primitivelattice,narms,kvecstar,kvecop)
    do s=1,narms
       do j=1,nk
          qk=kp(:,j)
          if (qk .myeq. kvecstar(:,s)) then
             ds=ds+1
             kvecstars(:,ds)=kvecstar(:,s)
             kvecops(ds)=kvecop(s)
!             sl=sl+1
!            write(190,*)sl,i,ds,kvecops(ds),kvecstars(:,ds)
          endif
        enddo
     enddo
!enddo



end subroutine kstarsinfbz1


!================================test-RTA
!subroutine test_RTA()

!use constants
!use exactBTE2
!use params
!use lattice
!use io2
!use mod_ksubset2
!use phi3
!use kpoints
!use force_constants_module
!implicit none

!do i=nk


!end subroutine test_RTA
!========================================================================================
subroutine Pmatrix3 (indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)

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

 integer indx, ksub_size2, k_shift, ksubset(2)
 integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
 integer nk1,nk2,nk3, nk1n, nk2n, nk3n, nk3_proc1, nk3_proc2, nr1,nr2,lr1,lr2,lr3
 integer nk
 real(8) q1(3), q2(3), q3(3), q1n(3), q2n(3), q3n(3), q3_proc1(3), q3_proc2(3)
 real(8) kp(3,nk)

 integer i,n1,n2,n3,l1,l2,l3,ndn
 real(8) temp     ! temperature in cm^-1
 real(8) omega1, omega2, omega3, V3sq1, V3sq2, col1, col2, omega3_proc1, omega3_proc2
 real(8) xx

 character(99) filename, filename_temp

 integer ucol, unt, nv3_2
 integer nn1,nn2,ll1,ll2,ll3
 
 real(8), allocatable :: v33sq(:,:,:,:,:)
! real cputim1,cputim2,cputim0, cputim1_wall, cputim2_wall, cputim0_wall, wtime
 real cputim1, cputim2, cputim0
 real(8)  cputim1_wall, cputim2_wall, cputim0_wall, wtime

 ucol=9010
 unt=112

 P1=0.0; P2=0.0

!call cpu_time(cputim0)
!cputim0_wall = wtime()

 allocate(v33sq(ksub_size2,nk,ndn,ndn,ndn))
 v33sq=0

! write(filename_temp,fmt="(a,i3.3,a)") "v33sq.",indx,".dat"
! filename=trim(v3path)//filename_temp
! open(unt,file=filename,status='old',form='unformatted')
! read(unt) nv3_2
!
! do n1=1,ksub_size2
! do l1=1,ndn
! do n2=1,nk
! do l2=1,ndn
! do l3=1,ndn
!     read(unt) v33sq(n1,n2,l1,l2,l3)
! enddo
! enddo
! enddo
! enddo
! enddo


!call cpu_time(cputim1)
!cputim1_wall = wtime()
!write(ulog,*) indx,' read v3sq.dat done, TIME=',nint(cputim1-cputim0),'WALL TIME=',nint(cputim1_wall-cputim0_wall)

loop1 : do n1=ksubset(1),ksubset(2)


    q1=kp(:,n1)
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
    if (n1 .ne. nk1) then
       write(ulog,*) 'n3,nk3,inside=',n3,nk3,inside
       write(ulog,*) 'q3=',q3
       write(ulog,*) 'ijk3=',i3,j3,k3
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

           q2n=-q2
           call get_k_info(q2n,NC,nk2n,i2,j2,k2,g1,g2,g3,inside)

           q3=-q1-q2
           call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)

           q3_proc1=q1+q2
           call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside)
           q3_proc2=q1+q2n
           call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside)

!           write(*,*) 'nk2,nk2n',nk2,nk2n           

           do l2=1, ndn
           do l3=1, ndn

                 V3sq1=v33sq(n1,n2,l1,l2,l3) * (2*pi)**2 / nk
!                 V3sq2=v33sq(n1,nk2n,l1,l2,l3) * (2*pi)**2 / nk

!                 omega1=frequency(nk1,l1)
!                 omega2=frequency(nk2,l2)
!                 omega3=frequency(nk3,l3)

                  call cal_P2(n1,l1,n2,l2,nk3,l3,V3sq1,V3sq1,col1,col2)
!                 call cal_P3(n1,l1,n2,l2,nk3_proc1,nk3_proc2,l3,V3sq1,V3sq2,col1,col2)
                 P1(n1-k_shift,n2,l1,l2,l3)=col1
!                 P2(n1-k_shift,n2,l1,l2,l3)=col2
                 if (P2(n1-k_shift,nk2n,l1,l2,l3) .ne. 0.0) then
                     write(ulog,*) 'P2 already calculated'
                     stop
                 endif
                 P2(n1-k_shift,nk2n,l1,l2,l3)=col2

            enddo
           enddo
!            enddo
     enddo loop2

enddo
enddo loop1

call cpu_time(cputim2)
cputim2_wall = wtime()
write(ulog,*) '   cal_P matrix done,  TIME=',nint(cputim2-cputim1),'WALL TIME=',nint(cputim2_wall-cputim1_wall)

!deallocate(v33sq)
close(unt)

7 format(5(i6),2(2x,g14.8))
8 format(i6,3(2x,g14.8))


end subroutine Pmatrix3



!==============================================================
subroutine check_conv2(i,nk,kp,ndn,tempk,convergence,iter_cont)

use constants
use exactBTE2
use params
use lattice
use io2
use phi3
use mod_ksubset2
use eigen

implicit none

integer nk,ndn,i,convergence,iter_cont
integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
integer nk1,nk2,nk3,xyz,xyz2
real(8) q1(3), q2(3)
real(8) kp(3,nk)
real(8) Avalue(nk,ndn,3)
integer n1,n2,n3,l1,l2,l3
real(8) leftBTE, rightBTE  ! lefthand term in BTE (3.78)
real(8) omega1, nb1, nbe, temp, tempk, norm_diff, max_diff, norm_error, max_error, norm_diff_kap, max_diff_kap

integer nk_subset, k_shift, indx, ksubset(2), ksub_size2
integer max_diff_kap_nk, max_diff_kap_la, max_diff_kap_xyz, max_diff_kap_xyz2

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

write(ulog,*) 'entering check_conv...'

!calculate norm(diff) and max(diff)
!check convergence
norm_diff_kap=0
max_diff_kap=0

do n1=1,nk
do l1=1,ndn
do xyz=1,3


do xyz2=1,3
norm_diff_kap=norm_diff_kap+(diff_kap(n1,l1,xyz,xyz2)**2/(nk*ndn*3*3))
if (abs(diff_kap(n1,l1,xyz,xyz2)) .gt. abs(max_diff_kap)) then
    max_diff_kap=diff_kap(n1,l1,xyz,xyz2)
    max_diff_kap_nk = n1
    max_diff_kap_la = l1
    max_diff_kap_xyz = xyz
    max_diff_kap_xyz2 = xyz2
endif
enddo
!write(*,*) 'n1,l1,xyz,F1',n1,l1,xyz,F1(n1,l1,xyz)

enddo
enddo
enddo

norm_error=sqrt(norm_error)
norm_diff=sqrt(norm_diff)

convergence=0
!if (norm_error.lt.conv_error .and. max_error.lt.conv_max_error .and. norm_diff.lt.conv_diff .and. max_diff.lt.conv_max_diff .and. norm_diff_kap.lt.conv_diff_kap .and. max_diff_kap.lt.conv_max_diff_kap) then
!convergence=1
!iter_cont=iter_cont+1
!else
!convergence=0
!iter_cont=0
!endif

if (norm_diff_kap.lt.conv_diff_kap .and. max_diff_kap.lt.conv_max_diff_kap) then
convergence=1
iter_cont=iter_cont+1
else
convergence=0
iter_cont=0
endif

write(ulog,'(a2,99(2x,g12.6))') '++',i,norm_diff_kap,max_diff_kap,max_diff_kap_nk,max_diff_kap_la,max_diff_kap_xyz,max_diff_kap_xyz2, kappa_k(max_diff_kap_nk,max_diff_kap_la,max_diff_kap_xyz,max_diff_kap_xyz2)*nk



end subroutine check_conv2
!===================================================================


!==================================================================
!==============================================================
subroutine check_conv23(i,nk,kp,nibbz,kibbz,ndn,tempk,convergence,iter_cont) !!Saf!!

use constants
use exactBTE2
use params
use lattice
use io2
use phi3
use mod_ksubset2
use eigen
use kpoints !!mine


implicit none

integer nk,ndn,i,convergence,iter_cont
integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
integer nk1,nk2,nk3,xyz,xyz2
real(8) q1(3), q2(3)
real(8) kp(3,nk)
real(8) Avalue(nk,ndn,3)
integer n1,n2,n3,l1,l2,l3
real(8) leftBTE, rightBTE  ! lefthand term in BTE (3.78)
real(8) omega1, nb1, nbe, temp, tempk, norm_diff, max_diff, norm_error,max_error, norm_diff_kap, max_diff_kap

integer nk_subset, k_shift, indx, ksubset(2), ksub_size2
integer max_diff_kap_nk, max_diff_kap_la, max_diff_kap_xyz, max_diff_kap_xyz2

integer s,n0,narms,kvecop(48),nibbz
real(8) kvecstar(3,48),kibbz(3,nibbz),q(3)


temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

write(ulog,*) 'entering check_conv...'

!calculate norm(diff) and max(diff)
!check convergence
norm_diff_kap=0
max_diff_kap=0

!do n1=1,nk
do n1=1,nibbz

             q=kibbz(:,n1)
!            call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
             call kstarsinfbz1(q,nk,kp,narms,kvecstar,kvecop)


do l1=1,ndn
do s=1,narms
   call get_k_info(kvecstar(:,s),NC,n0,i3,j3,k3,g1,g2,g3,inside)

do xyz=1,3


do xyz2=1,3
norm_diff_kap=norm_diff_kap+(diff_kap(n0,l1,xyz,xyz2)**2/(nk*ndn*3*3))
if (abs(diff_kap(n0,l1,xyz,xyz2)) .gt. abs(max_diff_kap)) then
    max_diff_kap=diff_kap(n0,l1,xyz,xyz2)
    max_diff_kap_nk = n0
    max_diff_kap_la = l1
    max_diff_kap_xyz = xyz
    max_diff_kap_xyz2 = xyz2
endif
enddo
!write(*,*) 'n1,l1,xyz,F1',n1,l1,xyz,F1(n1,l1,xyz)

enddo
enddo
enddo
enddo

norm_error=sqrt(norm_error)
norm_diff=sqrt(norm_diff)

convergence=0
!if (norm_error.lt.conv_error .and. max_error.lt.conv_max_error .and.
!norm_diff.lt.conv_diff .and. max_diff.lt.conv_max_diff .and.
!norm_diff_kap.lt.conv_diff_kap .and. max_diff_kap.lt.conv_max_diff_kap) then
!convergence=1
!iter_cont=iter_cont+1
!else
!convergence=0
!iter_cont=0
!endif

if (norm_diff_kap.lt.conv_diff_kap .and. max_diff_kap.lt.conv_max_diff_kap) then
convergence=1
iter_cont=iter_cont+1
else
convergence=0
iter_cont=0
endif

write(ulog,'(a2,99(2x,g12.6))') '++',i,norm_diff_kap,max_diff_kap,max_diff_kap_nk,max_diff_kap_la,max_diff_kap_xyz,max_diff_kap_xyz2,kappa_k(max_diff_kap_nk,max_diff_kap_la,max_diff_kap_xyz,max_diff_kap_xyz2)*nk



end subroutine check_conv23



!============================================================











!============================================================
 subroutine matrix_elt_simplified(q1,q2,l1,l2,l3,w33,inside,ndn)
! Exactly same as matrix_elt. But algebraic part has been simplified for reducing cpu time.
!!!!!!!!! WARNING!!!!!!! This subroutine works only when atomic masses are same!!!!!!!!!!!
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
 use pre_matrix_elt

 implicit none
 integer, intent(in) :: l1,l2,l3,ndn
 real(8), intent(in) :: q1(3),q2(3)
 integer, intent(out) :: inside
 complex(8), intent (out) :: w33
 integer i3,nk1,nk2,nk3,i0,al,be,ga,j,j0,k,k0,t,ta1,ta2,ta3,i1,j1,k1,num_basis
 real(8) q3(3),mi,mj,mk,rr2(3),rr3(3),den
 complex(8) xx,eiqr
 complex(8), allocatable :: eivec3(:,:,:,:,:,:)
 real cputim,cputim2

 integer xyz1,xyz2,xyz3,basis1,basis2,basis3

 num_basis=ndn/3
 allocate(eivec3(3,3,3,num_basis,num_basis,num_basis))


! call cpu_time(cputim)
! write(*,*) '----start matrix_elt at 0'
 
! write(*,*)'matrix_elt',q1,q2,l1,l2,l3
! write(*,*)'matrix_elt, const33=',const33
 call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
 call get_k_info(q2,NC,nk2,i1,j1,k1,g1,g2,g3,inside)
 w33 = cmplx(0d0,0d0)
 q3 = -q2-q1
 call get_k_info(q3,NC,nk3,i1,j1,k1,g1,g2,g3,inside)


 do xyz1=1,3
 do xyz2=1,3
 do xyz3=1,3
 do basis1=1,num_basis
 do basis2=1,num_basis
 do basis3=1,num_basis

 eivec3(xyz1,xyz2,xyz3,basis1,basis2,basis3)=eigenvec(xyz1+3*(basis1-1),l1,nk1)*eigenvec(xyz2+3*(basis2-1),l2,nk2)*&
& eigenvec(xyz3+3*(basis3-1),l3,nk3)

 enddo
 enddo
 enddo
 enddo
 enddo
 enddo



! write(*,*)'matrix_elt: q3=',q3,nk3,inside
 xx = cmplx(0d0,0d0)

! call cpu_time(cputim2)
! write(*,*) 'before tloop at ',cputim2-cputim

 tloop: do t=1,nterms(3)

! write(*,*)'matrix_elt: tloop=',t,nterms(3)

       i3 = iatomterm_3(1,t)
       if (iatomcell0(i3) .ne. i3) cycle tloop ! first index must be in primitive cell
       j  = iatomterm_3(2,t)
       k  = iatomterm_3(3,t)

!       xx = xx + fcs_3(t) * eiqr2(nk2,j) * eiqr2(nk3,k) * &
!&      eivec2(i3,ixyzterm_3(1,t),l1,nk1)*eivec2(j,ixyzterm_3(2,t),l2,nk2)*eivec2(k,ixyzterm_3(3,t),l3,nk3)

       xx = xx + fcs_3(t) * eiqr2(nk2,j) * eiqr2(nk3,k) * &
&      eivec3(ixyzterm_3(1,t),ixyzterm_3(2,t),ixyzterm_3(3,t),iatomcell0(i3),iatomcell0(j),iatomcell0(k))


enddo tloop

!call cpu_time(cputim)
!write(*,*) 'after tloop at ',cputim-cputim2

xx=xx/sqrt(atom0(iatomcell0(i3))%mass**3)   ! mass is same.
!       i3 = iatomterm_3(1,t)
!       i0 = iatomcell0(i3)
!       if (i0.ne.i3) cycle tloop ! first index must be in primitive cell
!       al = ixyzterm_3(1,t)
!       be = ixyzterm_3(2,t)
!       ga = ixyzterm_3(3,t)
!       mi = atom0(i0)%mass
!       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
!       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
!       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
!       rr2(:) = atompos(:,j) - atompos(:,j0)
!       rr3(:) = atompos(:,k) - atompos(:,k0)
!       eiqr = exp( ci* ((q2 .dot. rr2) + (q3 .dot. rr3)) )
!       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
!&      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)

! enddo tloop
 den = sqrt(8*sqrt(abs(eigenval(l1,nk1)*eigenval(l2,nk2)*eigenval(l3,nk3))))
 if (den.ne.0) then
    w33 = xx / den * const33
 else
    write(*,*)'MATRIX_ELT: den=0 ',den
    stop
 endif


 deallocate(eivec3)


 end subroutine matrix_elt_simplified




!============================================================
 subroutine matrix_elt_simplified2(q1,q2,l1,l2,l3,w33,inside,ndn)
 ! same as matrix_elt_simplified. But, use when basis atoms have different mass.

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
 use pre_matrix_elt

 implicit none
 integer, intent(in) :: l1,l2,l3,ndn
 real(8), intent(in) :: q1(3),q2(3)
 integer, intent(out) :: inside
 complex(8), intent (out) :: w33
 integer i3,nk1,nk2,nk3,i0,al,be,ga,j,j0,k,k0,t,ta1,ta2,ta3,i1,j1,k1,num_basis
 real(8) q3(3),mi,mj,mk,rr2(3),rr3(3),den
 complex(8) xx,eiqr
 complex(8), allocatable :: eivec3(:,:,:,:,:,:)
 real cputim,cputim2

 integer xyz1,xyz2,xyz3,basis1,basis2,basis3

 num_basis=ndn/3
 allocate(eivec3(3,3,3,num_basis,num_basis,num_basis))


! call cpu_time(cputim)
! write(*,*) '----start matrix_elt at 0'

! write(*,*)'matrix_elt',q1,q2,l1,l2,l3
! write(*,*)'matrix_elt, const33=',const33
 call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
 call get_k_info(q2,NC,nk2,i1,j1,k1,g1,g2,g3,inside)
 w33 = cmplx(0d0,0d0)
 q3 = -q2-q1
 call get_k_info(q3,NC,nk3,i1,j1,k1,g1,g2,g3,inside)


 do xyz1=1,3
 do xyz2=1,3
 do xyz3=1,3
 do basis1=1,num_basis
 do basis2=1,num_basis
 do basis3=1,num_basis

 eivec3(xyz1,xyz2,xyz3,basis1,basis2,basis3)=eigenvec(xyz1+3*(basis1-1),l1,nk1)*eigenvec(xyz2+3*(basis2-1),l2,nk2)*&
& eigenvec(xyz3+3*(basis3-1),l3,nk3)

 enddo
 enddo
 enddo
 enddo
 enddo
 enddo



! write(*,*)'matrix_elt: q3=',q3,nk3,inside
 xx = cmplx(0d0,0d0)

! call cpu_time(cputim2)
! write(*,*) 'before tloop at ',cputim2-cputim

 tloop: do t=1,nterms(3)

! write(*,*)'matrix_elt: tloop=',t,nterms(3)

       i3 = iatomterm_3(1,t)
       i0 = iatomcell0(i3)
       if (iatomcell0(i3) .ne. i3) cycle tloop ! first index must be in primitive cell
       j  = iatomterm_3(2,t)
       k  = iatomterm_3(3,t)

       mi = atom0(i0)%mass
       j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
       k0 = iatomcell0(k) ;  mk = atom0(k0)%mass

!       xx = xx + fcs_3(t) * eiqr2(nk2,j) * eiqr2(nk3,k) * &
!&      eivec2(i3,ixyzterm_3(1,t),l1,nk1)*eivec2(j,ixyzterm_3(2,t),l2,nk2)*eivec2(k,ixyzterm_3(3,t),l3,nk3)

       xx = xx + fcs_3(t) /sqrt(mi*mj*mk) * eiqr2(nk2,j) * eiqr2(nk3,k) * &
&      eivec3(ixyzterm_3(1,t),ixyzterm_3(2,t),ixyzterm_3(3,t),iatomcell0(i3),iatomcell0(j),iatomcell0(k))


enddo tloop

!call cpu_time(cputim)
!write(*,*) 'after tloop at ',cputim-cputim2

!       i3 = iatomterm_3(1,t)
!       i0 = iatomcell0(i3)
!       if (i0.ne.i3) cycle tloop ! first index must be in primitive cell
!       al = ixyzterm_3(1,t)
!       be = ixyzterm_3(2,t)
!       ga = ixyzterm_3(3,t)
!       mi = atom0(i0)%mass
!       j  = iatomterm_3(2,t)  ;  j0 = iatomcell0(j) ;  mj = atom0(j0)%mass
!       k  = iatomterm_3(3,t)  ;  k0 = iatomcell0(k) ;  mk = atom0(k0)%mass
!       ta1= al + 3*(i0-1)     ;  ta2= be + 3*(j0-1) ;  ta3= ga + 3*(k0-1)
!       rr2(:) = atompos(:,j) - atompos(:,j0)
!       rr3(:) = atompos(:,k) - atompos(:,k0)
!       eiqr = exp( ci* ((q2 .dot. rr2) + (q3 .dot. rr3)) )
!       xx = xx + fcs_3(t)  /sqrt(mi*mj*mk) * eiqr * &
!&      eigenvec(ta1,l1,nk1)*eigenvec(ta2,l2,nk2)*eigenvec(ta3,l3,nk3)

! enddo tloop
 den = sqrt(8*sqrt(abs(eigenval(l1,nk1)*eigenval(l2,nk2)*eigenval(l3,nk3))))
 if (den.ne.0) then
    w33 = xx / den * const33
 else
    write(*,*)'MATRIX_ELT: den=0 ',den
    stop
 endif


 deallocate(eivec3)


 end subroutine matrix_elt_simplified2





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
write(*,*) 'memory size required=',max_atomnum*3*ndn*nk
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

real(8), allocatable :: v33sq(:,:,:,:,:),tauinv_RTA(:,:)
real(8) tauinv_RTA_N(nk,ndn), tauinv_RTA_U(nk,ndn)
real(8) kappa_k_RTA2(nk,ndn,3,3), kappa_RTA2(ndn,3,3), nbe, delta_l

real(8) tauinv_temp1, tauinv_temp2

unt=112
ucol=9010

allocate(v33sq(nk,nk,ndn,ndn,ndn))
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

! if (mod(n1,ksub_size) .eq. 1) then   ! start point of each v33.xxx.dat


!call cpu_time(cputim0)
!cputim0_wall=wtime()


!    indx=indx+1
!    if (indx*ksub_size > nk) then
!        ksub_size2=nk - (indx-1)*ksub_size
!        ksubset(1)=ksub_size*(indx-1)+1    ! starting point of k subset
!        ksubset(2)=nk
!    else
!        ksub_size2=ksub_size
!        ksubset(1)=ksub_size*(indx-1)+1
!        ksubset(2)=ksub_size*indx
!    endif
!    k_shift=ksub_size*(indx-1)    ! P1 of n1,n2,l1,l2,l3 will be saved in P1(n1-k_shift,n2,l1,l2,l3).

    ! READ V33sq.xxx.dat and store it to V33sq
!   if (allocated(v33sq)) then
!   deallocate(v33sq)
!   endif


!   allocate(v33sq(ksub_size2,nk,ndn,ndn,ndn)) 
!   V33sq=0.0
! write(filename_temp,fmt="(a,i3.3,a)") "v33sq.",indx,".dat"
! filename=trim(v3path)//filename_temp
! open(unt,file=filename,status='old',form='unformatted')
! read(unt) nv3_2
! do n1_2=1,ksub_size2
! do l1_2=1,ndn
! do n2_2=1,nk
! do l2_2=1,ndn
! do l3_2=1,ndn
!     read(unt) v33sq(n1_2,n2_2,l1_2,l2_2,l3_2)
!     read(unt,*) nn1,nn2,ll1,ll2,ll3,v33sq(nn1-k_shift,nn2,ll1,ll2,ll3)
!     if (n1 .ne. nn1-k_shift) then
!         write(*,*) 'ERRROR while reading v33'
!         stop
!    endif
! enddo
! enddo
! enddo
! enddo
! enddo
! call cpu_time(cputim1)
! cputim1_wall = wtime()
! write(ulog,*) indx, ' read v3sq.dat done, TIME=', nint(cputim1-cputim0),'WALL TIME=', nint(cputim1_wall-cputim0_wall)
! write(*,*) 'read v33sq.dat done'
!  ---- End of reading V33sq.dat

! endif

    q1=kp(:,n1)
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
    if (n1 .ne. nk1) then
       write(ulog,*) 'n1,nk1,inside=',n1,nk1,inside
       write(ulog,*) 'q1=',q1
       write(ulog,*) 'ijk1=',i1,j1,k1
       stop
    endif

!if (n1 .eq. 21) then
!write(*,*) 'q1',q1
!write(*,*) 'nk1,i1,j1,k1 of n1=21',nk1,i1,j1,k1
!call get_k_info(-q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
!write(*,*) '-q1',-q1
!write(*,*) 'nk1,i1,j1,k1 of -q1', nk1,i1,j1,k1
!stop
!endif


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

!           write(*,*) 'q2=',q2
!           call get_k_info(-q2,NC,nk2n,i2,j2,k2,g1,g2,g3,inside)
!           write(*,*) 'nk2,nk2n=',nk2,nk2n
           !write(*,*) 'q2=',q2


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


!if (n1 .eq. 445) then
!if (n1 .eq. 667) then
!write(ulog,*) 'q1=',q1
!write(ulog,*) 'inside_proc1=',inside_proc1
!endif

!100 format(5(2x,i3))


           do l2=1, ndn
           do l3=1, ndn

                 V3sq1=v33sq(n1,n2,l1,l2,l3) * (2*pi)**2 / nk

                 omega1=frequency(nk1,l1)
                 omega2=frequency(nk2,l2)
                 omega3=frequency(nk3,l3)
                 omega3_proc1=frequency(nk3_proc1,l3)
                 omega3_proc2=frequency(nk3_proc2,l3)

!if (n1 .eq. 667) then
!write(ulog,*) 'nk2,inside',n2,inside
!endif

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

!call cpu_time(cputim2)
!cputim2_wall = wtime()
!write(ulog,*) '   cal_P matrix done,  TIME=',cputim2-cputim1,'WALL TIME=',cputim2_wall-cputim1_wall
!write(*,*) '    cal_P matrix done, time=', cputim2-cputim1

do l1=1,ndn
do n1=1,nk
write(unt+11,11) l1, n1, tempk, frequency(n1,l1), 0.5*tauinv_RTA(n1,l1), 0.5*tauinv_RTA_N(n1,l1), 0.5*tauinv_RTA_U(n1,l1)  ! print out selfenergy
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

               kappa_k_RTA2(i,ii,j,jj) = (omega1*c_light*100.)**2.0 * (veloc(j,ii,i)*c_light) * (veloc(jj,ii,i)*c_light) * nb1 * (nb1+1) * (1/(tauinv_RTA(i,ii)*c_light*100.))
               kappa_k_RTA2(i,ii,j,jj) = kappa_k_RTA2(i,ii,j,jj) * h_plank**2.0 / (nk*volume_r/1d30) / (k_b*tempk**2.0)
               kappa_RTA2(ii,j,jj)=kappa_RTA2(ii,j,jj) + kappa_k_RTA2(i,ii,j,jj)

          enddo beta_loop
          enddo alpha_loop

       enddo laF1_loop

    !   write(*,*) kappa_RTA2(1,2,2), kappa_RTA2(5,2,2)
enddo nkF1_loop

!write(*,*) 'veloc at gamma'
!write(*,*) veloc(:,1,43)
!write(*,*) veloc(:,2,43)
!write(*,*) veloc(:,3,43)
!write(*,*) veloc(:,4,43)
!write(*,*) veloc(:,5,43)
!write(*,*) veloc(:,6,43)

!----------- print out kapp_RTA2

open(ucol,file='kappa_RTA.dat',status='unknown')
write(ucol,*) 'xx,yy,zz, xy,xz,yz'
do ii=1,ndn
   write(ucol,9) ii,kappa_RTA2(ii,1,1),kappa_RTA2(ii,2,2),kappa_RTA2(ii,3,3),kappa_RTA2(ii,1,2),kappa_RTA2(ii,1,3),kappa_RTA2(ii,2,3)
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
subroutine cal_Qm3(nk,ndn,kp,tempk)
! calculate sum of diagonal terms in collision matrix
! same as cal_Qm2 but calculate Pmatrix and save it to memory


use io2
use exactBTE2
use mod_ksubset2
use params
use constants
use phi3
use lattice

implicit none

integer i,ii,indx,ndn,j,jj,kk
integer nk, nk_subset, k_shift, ksub_size2
real(8) kp(3,nk), q1(3), q2(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
real(8) temp, tempk, nb1   ! in cm^-1
integer ksubset(2)

integer i2,ii2,j2,jj2,kk2, inside1, inside2, inside, inside1n, inside2n
integer i3,j3,k3,nk3_proc1,nk3_proc2, nk3

integer uQ

real cputim1, cputim2, cputim3

uQ=9402
open(uQ,file='Qvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Q'

!open(uQ+1,file='Pmatrix_paral.dat',status='unknown')

Qvalue(:,:)=0    ! in cm^-1
Qvalue_N(:,:)=0
Qvalue_U(:,:)=0
tauinv_N(:,:)=0
tauinv_U(:,:)=0
indx=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

!write(*,*) 'entering cal_Qm2...'

call cpu_time(cputim1)
write(*,*) 'entering cal_Qm2...'
nkF1_loop: do i=1,nk

 ! call Pmatrix3(i)
     ! -> read V33.xxx.dat
     ! -> allocate P1,P2 subset
     ! -> calculate P1, P2
     ! -> close V33.xxx.dat
 if (mod(i,ksub_size) .eq. 1) then   ! start point of each v33.xxx.dat

!    if(allocated(P1)) then
!        call deallocate_iter2
!    endif

    indx=indx+1
    if (indx*ksub_size > nk) then
        ksub_size2=nk - (indx-1)*ksub_size
        ksubset(1)=ksub_size*(indx-1)+1    ! starting point of k subset
        ksubset(2)=nk
    else
        ksub_size2=ksub_size
        ksubset(1)=ksub_size*(indx-1)+1
        ksubset(2)=ksub_size*indx
    endif
    k_shift=ksub_size*(indx-1)    ! P1 of n1,n2,l1,l2,l3 will be saved in P1(n1-k_shift,n2,l1,l2,l3).

!    call allocate_iter2(ksub_size2,nk,ndn,ndn,ndn)
    call Pmatrix4(indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)

    call cpu_time(cputim2)
    write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1

!    do i2=ksubset(1),ksubset(2)
!    do ii2=1,ndn
!    do j2=1,nk
!    do jj2=1,ndn
!    do kk2=1,ndn
!       write(uQ+1,7) i2,j2,ii2,jj2,kk2,P1(i2-k_shift,j2,ii2,jj2,kk2),P2(i2-k_shift,j2,ii2,jj2,kk2)
!    enddo
!    enddo
!    enddo
!    enddo
!    enddo

!    write(*,*) 'In cal_Qm2,Pmatrix3 done,i=',i
!    write(*,*) 'In cal_Qm2,k_shift',k_shift
 endif

 !    if (i .eq. ksubset2(indx2,1)) then

 !        k_shift=ksubset2(indx2,1)-1
 !        nk_subset=ksubset2(indx2,2)-ksubset2(indx2,1)+1
 !        call allocate_iter2(nk_subset,nk,ndn,ndn,ndn)
 !        call read_Pmatrix(indx2,nk,ndn,k_shift)

 !   endif

    laF1_loop: do ii=1,ndn

           nkF2_loop: do j=1,nk
           laF2_loop: do jj=1,ndn
           laF3_loop: do kk=1,ndn

                Qvalue(i,ii)=Qvalue(i,ii) + P1(i,j,ii,jj,kk) + 0.5*P2(i,j,ii,jj,kk)

                q1=kp(:,i)
                q2=kp(:,j)

                q3=-q1-q2
                q3_proc1=q1+q2
                q3_proc2=q1-q2
!                q3_proc1n=-q1-q2
!                q3_proc2n=-q1+q2

                call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)
!                call get_k_info(q3_proc1n,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1n)
!                call get_k_info(q3_proc2n,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2n)

                if (inside1 .eq. 1) then    ! normal process in coalescence process
                    Qvalue_N(i,ii)=Qvalue_N(i,ii) + P1(i,j,ii,jj,kk)
                else
                    Qvalue_U(i,ii)=Qvalue_U(i,ii) + P1(i,j,ii,jj,kk)
                endif

                if (inside2 .eq. 1) then    ! normal process in decay process
                    Qvalue_N(i,ii)=Qvalue_N(i,ii) + 0.5*P2(i,j,ii,jj,kk)
                else
                    Qvalue_U(i,ii)=Qvalue_U(i,ii) + 0.5*P2(i,j,ii,jj,kk)
                endif

! determine k vector
! determine i,j,k is normal or Umklapp
! If normal, Qvalue_n = Qvalue_n + P1 + 0.5*P2
! If Umklapp, Qvalue_u = Qvalue_u + P1 + 0.5*P2

           enddo laF3_loop
           enddo laF2_loop
           enddo nkF2_loop

           write(uQ,*) i,ii,Qvalue(i,ii)

           nb1=dist(i,ii)
           tauinv_N(i,ii) = Qvalue_N(i,ii)/(nb1*(nb1+1)) * c_light*100    ! tauinv_N in 1/sec
           tauinv_U(i,ii) = Qvalue_U(i,ii)/(nb1*(nb1+1)) * c_light*100
           tauinv_tot(i,ii) = tauinv_N(i,ii) + tauinv_U(i,ii)

      enddo laF1_loop

! give a sign to move to next v33.xxx.dat file


!      if (i .eq. ksubset2(indx2,2)) then
!          indx2=indx2+1
!          call deallocate_iter2
!      endif

enddo nkF1_loop

close(uQ)
!close(uQ+1)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif


end subroutine cal_Qm3



!========================================================================================
subroutine Pmatrix4 (indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)
! calculate collision matrix P1 and P2


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

 integer indx, ksub_size2, k_shift, ksubset(2)
 integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside, i2n, j2n, k2n
 integer nk1,nk2,nk3, nk1n, nk2n, nk3n, nk3_proc1, nk3_proc2, nr1,nr2,lr1,lr2,lr3
 integer nk
 real(8) q1(3), q2(3), q3(3), q1n(3), q2n(3), q3n(3), q3_proc1(3), q3_proc2(3)
 real(8) kp(3,nk)

 integer i,n1,n2,n3,l1,l2,l3,ndn
 real(8) temp     ! temperature in cm^-1
 real(8) omega1, omega2, omega3, V3sq1, V3sq2, col1, col2, omega3_proc1, omega3_proc2
 real(8) xx

 character(99) filename, filename_temp

 integer ucol, unt, nv3_2
 integer nn1,nn2,ll1,ll2,ll3

 real(8), allocatable :: v33sq(:,:,:,:,:)
! real cputim1,cputim2,cputim0, cputim1_wall, cputim2_wall, cputim0_wall, wtime
 real cputim1, cputim2, cputim0
 real(8)  cputim1_wall, cputim2_wall, cputim0_wall, wtime

 ucol=9010
 unt=112


call cpu_time(cputim0)
cputim0_wall = wtime()

 allocate(v33sq(ksub_size2,nk,ndn,ndn,ndn))
 v33sq=0

 write(filename_temp,fmt="(a,i3.3,a)") "v33sq.",indx,".dat"
 filename=trim(v3path)//filename_temp
 open(unt,file=filename,status='old',form='unformatted')
 read(unt) nv3_2

 do n1=1,ksub_size2
 do l1=1,ndn
 do n2=1,nk
 do l2=1,ndn
 do l3=1,ndn
     read(unt) v33sq(n1,n2,l1,l2,l3)
 enddo
 enddo
 enddo
 enddo
 enddo


call cpu_time(cputim1)
cputim1_wall = wtime()
write(ulog,*) indx,' read v3sq.dat done, TIME=',nint(cputim1-cputim0),'WALL TIME=',nint(cputim1_wall-cputim0_wall)

loop1 : do n1=ksubset(1),ksubset(2)


    q1=kp(:,n1)
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
    if (n1 .ne. nk1) then
       write(ulog,*) 'n3,nk3,inside=',n3,nk3,inside
       write(ulog,*) 'q3=',q3
       write(ulog,*) 'ijk3=',i3,j3,k3
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

           q2n=-q2
           call get_k_info(q2n,NC,nk2n,i2n,j2n,k2n,g1,g2,g3,inside)

           q3_proc1=q1+q2
           call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside)
           q3_proc2=q1+q2n
           call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside)
           q3=-q1-q2
           call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)


           do l2=1, ndn
           do l3=1, ndn

                     V3sq1=v33sq(n1-k_shift,n2,l1,l2,l3) * (2*pi)**2 / nk

                 ! below is for correction for process including gamma point.
                 ! even number should be used for k-mesh size to correctly determine gamma point.

                 if (l2<=3 .and. i2.eq.NC(1)/2+1 .and. j2.eq.NC(2)/2+1 .and. k2.eq.NC(3)/2+1 ) then  ! Acoustic modes at Gamma point
                    V3sq1=0.0                                                                        ! to enforce matrix_elt = 0
!                    write(ulog,9) 'q2,l2',q2,l2
                 endif
                 if (l3<=3 .and. i3.eq.NC(1)/2+1 .and. j3.eq.NC(3)/2+1 .and. k3.eq.NC(3)/2+1 ) then
                    V3sq1=0.0
!                    write(ulog,9) 'q3,l3',q3,l3
                 endif

!                 V3sq2=v33sq(n1-k_shift,nk2n,l1,l2,l3) * (2*pi)**2 / nk

!                 omega1=frequency(nk1,l1)
!                 omega2=frequency(nk2,l2)
!                 omega3=frequency(nk3,l3)

                  call cal_P2(n1,l1,n2,l2,nk3,l3,V3sq1,V3sq1,col1,col2)   ! calculate collision matrix elements
!                 call cal_P3(n1,l1,n2,l2,nk3_proc1,nk3_proc2,l3,V3sq1,V3sq2,col1,col2) 
                 P1(n1,n2,l1,l2,l3)=col1                                  ! collision matrix element for process 1
!                 P2(n1-k_shift,n2,l1,l2,l3)=col2
                 if (P2(n1,nk2n,l1,l2,l3) .ne. 0.0) then
                     write(ulog,*) 'P2 already calculated'
                     stop
                 endif
                 P2(n1,nk2n,l1,l2,l3)=col2    ! for P2(q1,q2), v3(q1,-q2) should be used.

            enddo
           enddo
!            enddo
     enddo loop2

enddo
enddo loop1

call cpu_time(cputim2)
cputim2_wall = wtime()
write(ulog,*) '   cal_P matrix done,  TIME=',nint(cputim2-cputim1),'WALL TIME=',nint(cputim2_wall-cputim1_wall)

deallocate(v33sq)
close(unt)

7 format(5(i6),2(2x,g14.8))
8 format(i6,3(2x,g14.8))
9 format(99(2x,g14.8))

end subroutine Pmatrix4





!============================================================================================
subroutine cal_Qm3_tet(nk,ndn,kp,tempk)
! calculate sum of diagonal terms in collision matrix
! same as cal_Qm2 but calculate Pmatrix and save it to memory


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
real(8), allocatable :: v33sq(:,:,:,:,:), P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)

nx=NC(1); ny=NC(2); nz=NC(3)

uQ=9402
unt=9400
open(uQ,file='Qvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Q'

!open(uQ+1,file='Pmatrix_paral.dat',status='unknown')

Qvalue(:,:)=0    ! in cm^-1
Qvalue_N(:,:)=0
Qvalue_U(:,:)=0
tauinv_N(:,:)=0
tauinv_U(:,:)=0
indx=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

!write(*,*) 'entering cal_Qm2...'

call cpu_time(cputim1)
write(*,*) 'entering cal_Qm2...'
nkF1_loop: do i=1,nk

 ! call Pmatrix3(i)
     ! -> read V33.xxx.dat
     ! -> allocate P1,P2 subset
     ! -> calculate P1, P2
     ! -> close V33.xxx.dat
 if (mod(i,ksub_size) .eq. 1) then   ! start point of each v33.xxx.dat

    if(allocated(P1sm) .or. allocated(P2sm)) then
        deallocate(P1sm)
        deallocate(P2sm)
    endif

    indx=indx+1
    if (indx*ksub_size > nk) then
        ksub_size2=nk - (indx-1)*ksub_size
        ksubset(1)=ksub_size*(indx-1)+1    ! starting point of k subset
        ksubset(2)=nk
    else
        ksub_size2=ksub_size
        ksubset(1)=ksub_size*(indx-1)+1
        ksubset(2)=ksub_size*indx
    endif
    k_shift=ksub_size*(indx-1)    ! P1 of n1,n2,l1,l2,l3 will be saved in P1(n1-k_shift,n2,l1,l2,l3).

!    call Pmatrix4(indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)
    call cpu_time(cputim0)
    cputim0_wall = wtime()

    if (allocated(v33sq)) then
          deallocate(v33sq)
    endif
    allocate(v33sq(ksub_size2,nk,ndn,ndn,ndn))
    allocate(P1sm(ksub_size2,nk,ndn,ndn,ndn))
    allocate(P2sm(ksub_size2,nk,ndn,ndn,ndn))

    v33sq=0d0
    P1sm=0d0
    P2sm=0d0

    write(filename_temp,fmt="(a,i3.3,a)") "v33sq.",indx,".dat"
    filename=trim(v3path)//filename_temp
    open(unt,file=filename,status='old',form='unformatted')
    read(unt) nv3_2

!    write(*,*) 'nv3_2=',nv3_2

    do n1=1,ksub_size2
    do l1=1,ndn
    do n2=1,nk
    do l2=1,ndn
    do l3=1,ndn
       read(unt) v33sq(n1,n2,l1,l2,l3)
!write(*,*) 'n1,l1,n2,l2,l3',n1,l1,n2,l2,l3
    enddo
    enddo
    enddo
    enddo
    enddo

    call cpu_time(cputim1)
    cputim1_wall = wtime()
    write(ulog,*) indx,' read v3sq.dat done, TIME=',nint(cputim1-cputim0),'WALL TIME=',nint(cputim1_wall-cputim0_wall)

  endif
 
! Here call calc_tet
! Calculate Q_N, Q_U, Q_tot

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

      
                V3sq1=v33sq(i-k_shift,j,ii,jj,kk) * (2*pi)**2 / nk
                if (jj<=3 .and. i2.eq.NC(1)/2+1 .and. j2.eq.NC(2)/2+1 .and. k2.eq.NC(3)/2+1 ) then  ! Acoustic modes at Gamma point
                    V3sq1=0.0                                                                        ! to enforce matrix_elt = 0
                endif
                if (kk<=3 .and. i3.eq.NC(1)/2+1 .and. j3.eq.NC(2)/2+1 .and. k3.eq.NC(3)/2+1 ) then
                    V3sq1=0.0
                endif

                F_arbi(j)=dist(i,ii)*dist(j,jj)*(dist(nk3_proc1,kk)+1)*V3sq1   ! no 2pi at the front since delta is in linear frequency.

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

                    tet(j_tet)%p(k_tet)%w=-frequency(w,jj)+frequency(nk3_proc1,kk)     !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
               enddo
            enddo

           call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors

           do j_tet=1,nk*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                   w = tet(j_tet)%p(k_tet)%i                ! k point
                   P1sm(i-k_shift,w,ii,jj,kk) = P1sm(i-k_shift,w,ii,jj,kk) + tet(j_tet)%p(k_tet)%c/6d0 * F_arbi(w)      ! calculate P1
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


                V3sq1=v33sq(i-k_shift,j,ii,jj,kk) * (2*pi)**2 / nk
                if (jj<=3 .and. i2n.eq.NC(1)/2+1 .and. j2n.eq.NC(2)/2+1 .and. k2n.eq.NC(3)/2+1 ) then  ! Acoustic modes at Gamma point
                    V3sq1=0.0                                                                        ! to enforce matrix_elt = 0
                endif
                if (kk<=3 .and. i3.eq.NC(1)/2+1 .and. j3.eq.NC(2)/2+1 .and. k3.eq.NC(3)/2+1 ) then
                    V3sq1=0.0
                endif

                F_arbi(nk2n)=dist(i,ii)*(dist(j,jj)+1)*(dist(nk3_proc2,kk)+1)*V3sq1        ! no 2pi at the front since delta is in linear frequency.

            enddo

            ! eigen_tet: set omega2+omega3
            do j_tet=1,nk*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                    w = tet(j_tet)%p(k_tet)%i                          !label for kpt(3,i)

                    q1=kp(:,i)
                    q2=kp(:,w)

!                    q3=-q1-q2
                    q3_proc2=q1-q2

!                    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                    call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside1)

                    tet(j_tet)%p(k_tet)%w=frequency(w,jj)+frequency(nk3_proc2,kk)     !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
               enddo
            enddo

           call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors

           do j_tet=1,nk*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                   w = tet(j_tet)%p(k_tet)%i                ! k point
                   P2sm(i-k_shift,w,ii,jj,kk) = P2sm(i-k_shift,w,ii,jj,kk) + tet(j_tet)%p(k_tet)%c/6d0 * F_arbi(w)      ! calculate P1
               enddo
           enddo

           do j=1,nk

                Qvalue(i,ii)=Qvalue(i,ii) + P1sm(i-k_shift,j,ii,jj,kk) + 0.5*P2sm(i-k_shift,j,ii,jj,kk)

                q1=kp(:,i)
                q2=kp(:,j)

                q3=-q1-q2
                q3_proc1=q1+q2
                q3_proc2=q1-q2

                call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)

                if (inside1 .eq. 1) then    ! normal process in coalescence process
                    Qvalue_N(i,ii)=Qvalue_N(i,ii) + P1sm(i-k_shift,j,ii,jj,kk)
                else
                    Qvalue_U(i,ii)=Qvalue_U(i,ii) + P1sm(i-k_shift,j,ii,jj,kk)
                endif

                if (inside2 .eq. 1) then    ! normal process in decay process
                    Qvalue_N(i,ii)=Qvalue_N(i,ii) + 0.5*P2sm(i-k_shift,j,ii,jj,kk)
                else
                    Qvalue_U(i,ii)=Qvalue_U(i,ii) + 0.5*P2sm(i-k_shift,j,ii,jj,kk)
                endif
 
           enddo

 
           enddo laF3_loop
           enddo laF2_loop

           write(uQ,*) i,ii,Qvalue(i,ii)

           nb1=dist(i,ii)
           tauinv_N(i,ii) = Qvalue_N(i,ii)/(nb1*(nb1+1)) * c_light*100    ! tauinv_N in 1/sec
           tauinv_U(i,ii) = Qvalue_U(i,ii)/(nb1*(nb1+1)) * c_light*100
           tauinv_tot(i,ii) = tauinv_N(i,ii) + tauinv_U(i,ii)



    enddo laF1_loop



    call cpu_time(cputim2)
    write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1


enddo nkF1_loop



close(uQ)
!close(uQ+1)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif


end subroutine cal_Qm3_tet




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
real(8), allocatable :: v33sq(:,:,:,:,:), P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)
! P1sm and P2sm are v33sq*delta(proc1) and v33sq*delta(proc2)

nx=NC(1); ny=NC(2); nz=NC(3)
unt=9400

call allocate_iter1(nk,ndn)
call cal_freq(nk,ndn)
write(ulog,*) 'cal_freq done'

write(ulog,*) 'entering calculate_v3sq_delta ...'

! calculate indx using nk and ksubset
! calculate ksub_size2
if (ksubset(2) .eq. nk) then   ! if this is the last file in v33sq.xxx.dat
   ksub_size2=ksubset(2)-ksubset(1)+1
   indx=(ksubset(1)-1)/ksub_size + 1
   if (mod(ksubset(1)-1,ksub_size) .ne. 0 ) then
       write(ulog,*) 'wrong indx. stop',ksubset(1),ksub_size
       stop
   endif

else
    ksub_size2=ksubset(2)-ksubset(1)+1
    indx=ksubset(2)/ksub_size           ! ksub_size is from params.phon
    if (mod(ksubset(2),ksub_size) .ne. 0) then
       write(ulog,*) 'wrong indx. stop',ksubset(2),ksub_size
       stop
    endif

endif

k_shift=ksubset(1)-1

write(ulog,*) 'Check indx is same as folder name, indx=',indx



! Read v33sq from v33sq.xxx.dat in the main directory
write(filename_temp,fmt="(a,i3.3,a)") "v33sq.",indx,".dat"
filename=trim(v3path)//filename_temp
write(ulog,*) 'opening file ',filename
open(unt,file=filename,status='old',form='unformatted')
read(unt) nv3_2

if (allocated(v33sq)) then
   deallocate(v33sq)
endif

allocate(v33sq(ksub_size2,nk,ndn,ndn,ndn))

do n1=1,ksub_size2
do l1=1,ndn
do n2=1,nk
do l2=1,ndn
do l3=1,ndn
       read(unt) v33sq(n1,n2,l1,l2,l3)
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


                V3sq1=v33sq(i-k_shift,j,ii,jj,kk) * (2*pi)**2 / nk
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

                    tet(j_tet)%p(k_tet)%w=-frequency(w,jj)+frequency(nk3_proc1,kk)     !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
               enddo
            enddo

           call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors in tetrahedron method

           do j_tet=1,nk*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                   w = tet(j_tet)%p(k_tet)%i                ! k point
                   P1sm(i-k_shift,w,ii,jj,kk) = P1sm(i-k_shift,w,ii,jj,kk) + tet(j_tet)%p(k_tet)%c/6d0 * F_arbi(w)      ! calculate P1
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


                V3sq1=v33sq(i-k_shift,j,ii,jj,kk) * (2*pi)**2 / nk
                if (jj<=3 .and. i2n.eq.NC(1)/2+1 .and. j2n.eq.NC(2)/2+1 .and. k2n.eq.NC(3)/2+1 ) then  ! Acoustic modes at Gamma point
                    V3sq1=0.0                                                                        ! to enforce matrix_elt = 0
                endif
                if (kk<=3 .and. i3.eq.NC(1)/2+1 .and. j3.eq.NC(2)/2+1 .and. k3.eq.NC(3)/2+1 ) then
                    V3sq1=0.0
                endif

                F_arbi(nk2n)=V3sq1        ! no 2pi at the front since delta is in linear frequency.

            enddo

            ! eigen_tet: set omega2+omega3
            do j_tet=1,nk*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                    w = tet(j_tet)%p(k_tet)%i                          !label for kpt(3,i)

                    q1=kp(:,i)
                    q2=kp(:,w)

!                    q3=-q1-q2
                    q3_proc2=q1-q2

!                    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                    call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside1)

                    tet(j_tet)%p(k_tet)%w=frequency(w,jj)+frequency(nk3_proc2,kk)     !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
               enddo
            enddo

           call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors

           do j_tet=1,nk*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                   w = tet(j_tet)%p(k_tet)%i                ! k point
                   P2sm(i-k_shift,w,ii,jj,kk) = P2sm(i-k_shift,w,ii,jj,kk) + tet(j_tet)%p(k_tet)%c/6d0 * F_arbi(w)      ! calculate P1
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
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif



end subroutine calculate_v3sq_delta





!============================================================================================
subroutine cal_Qm3_tet2(nk,ndn,kp,tempk)
! calculate sum of diagonal terms in collision matrix
! same as cal_Qm2 but calculate Pmatrix and save it to memory


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
real(8), allocatable :: v33sq(:,:,:,:,:), P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)

nx=NC(1); ny=NC(2); nz=NC(3)

uQ=9402
unt=9400
open(uQ,file='Qvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Q'

!open(uQ+1,file='Pmatrix_paral.dat',status='unknown')

Qvalue(:,:)=0    ! in cm^-1
Qvalue_N(:,:)=0
Qvalue_U(:,:)=0
tauinv_N(:,:)=0
tauinv_U(:,:)=0
indx=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

!write(*,*) 'entering cal_Qm2...'

call cpu_time(cputim1)
write(*,*) 'entering cal_Qm2...'
nkF1_loop: do i=1,nk

 ! call Pmatrix3(i)
     ! -> read V33.xxx.dat
     ! -> allocate P1,P2 subset
     ! -> calculate P1, P2
     ! -> close V33.xxx.dat
 if (mod(i,ksub_size) .eq. 1) then   ! start point of each v33.xxx.dat

    if(allocated(P1sm) .or. allocated(P2sm)) then
        deallocate(P1sm)
        deallocate(P2sm)
    endif

    indx=indx+1
    if (indx*ksub_size > nk) then
        ksub_size2=nk - (indx-1)*ksub_size
        ksubset(1)=ksub_size*(indx-1)+1    ! starting point of k subset
        ksubset(2)=nk
    else
        ksub_size2=ksub_size
        ksubset(1)=ksub_size*(indx-1)+1
        ksubset(2)=ksub_size*indx
    endif
    k_shift=ksub_size*(indx-1)    ! P1 of n1,n2,l1,l2,l3 will be saved in P1(n1-k_shift,n2,l1,l2,l3).

!    call Pmatrix4(indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)
    call cpu_time(cputim0)
    cputim0_wall = wtime()

    allocate(P1sm(ksub_size2,nk,ndn,ndn,ndn))
    allocate(P2sm(ksub_size2,nk,ndn,ndn,ndn))

    P1sm=0d0
    P2sm=0d0

! read v33sq_delta.xxx.dat

    write(filename_temp,fmt="(a,i3.3,a)") "v33sq_delta.",indx,".dat"
    filename=trim(v3path)//filename_temp
    open(unt,file=filename,status='old',form='unformatted')
    read(unt) nv3_2


    do n1=1,ksub_size2
    do l1=1,ndn
    do n2=1,nk
    do l2=1,ndn
    do l3=1,ndn
       read(unt) P1sm(n1,n2,l1,l2,l3), P2sm(n1,n2,l1,l2,l3)
    enddo
    enddo
    enddo
    enddo
    enddo

    call cpu_time(cputim1)
    cputim1_wall = wtime()
    write(ulog,*) indx,' read v3sq_delta.dat done, TIME=',nint(cputim1-cputim0),'WALL TIME=',nint(cputim1_wall-cputim0_wall)

  endif

! Here call calc_tet
! Calculate Q_N, Q_U, Q_tot

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
                q2n=-kp(:,j)

                q3=-q1-q2
                q3_proc1=q1+q2
                q3_proc2=q1-q2

                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)


                F_arbi(j)=dist(i,ii)*dist(j,jj)*(dist(nk3_proc1,kk)+1)   ! no 2pi at the front since delta is in linear frequency.
                P1(i,j,ii,jj,kk) = P1sm(i-k_shift,j,ii,jj,kk) * F_arbi(j)    

                F_arbi(j)=dist(i,ii)*(dist(j,jj)+1)*(dist(nk3_proc2,kk)+1)   
                P2(i,j,ii,jj,kk) = P2sm(i-k_shift,j,ii,jj,kk) * F_arbi(j)     

                Qvalue(i,ii)=Qvalue(i,ii) + P1(i,j,ii,jj,kk) + 0.5*P2(i,j,ii,jj,kk)

                if (inside1 .eq. 1) then    ! normal process in coalescence process
                    Qvalue_N(i,ii)=Qvalue_N(i,ii) + P1(i,j,ii,jj,kk)
                else
                    Qvalue_U(i,ii)=Qvalue_U(i,ii) + P1(i,j,ii,jj,kk)
                endif

                if (inside2 .eq. 1) then    ! normal process in decay process
                    Qvalue_N(i,ii)=Qvalue_N(i,ii) + 0.5*P2(i,j,ii,jj,kk)
                else
                    Qvalue_U(i,ii)=Qvalue_U(i,ii) + 0.5*P2(i,j,ii,jj,kk)
                endif

           enddo


           enddo laF3_loop
           enddo laF2_loop

           nb1=dist(i,ii)

           write(uQ,*) i,ii,Qvalue(i,ii),Qvalue(i,ii)/(nb1*(nb1+1))

           tauinv_N(i,ii) = Qvalue_N(i,ii)/(nb1*(nb1+1)) * c_light*100    ! tauinv_N in 1/sec
           tauinv_U(i,ii) = Qvalue_U(i,ii)/(nb1*(nb1+1)) * c_light*100
           tauinv_tot(i,ii) = tauinv_N(i,ii) + tauinv_U(i,ii)



    enddo laF1_loop



    call cpu_time(cputim2)
    write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1


enddo nkF1_loop



close(uQ)
!close(uQ+1)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif


end subroutine cal_Qm3_tet2


