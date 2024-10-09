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
 use ios
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
      sumz(l2) = sumz(l2) + v33sq(j)  !(v33(j)*conjg(v33(j)))
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
 use ios
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
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
 
 v33sq = 0d0 !cmplx(0d0,0d0)
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
    v33sq(indx)=xx*conjg(xx)
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
   open(uv3,file='v33sq.dat',status='unknown',FORM='UNFORMATTED')
   !open(uv3,file='v33.dat',status='unknown')
   write(uv3)nv3_split
   !write(uv3,*) nv3_split
   do j=1 ,nv3_split
      write(uv3) nq1(j),nq2(j),nq3(j),la1(j),la2(j),la3(j),v33sq(j)
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
 use ios
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



!============================================================
 subroutine function_self_w2_sy(q_ibz,q,la,omega,temp,nself,uself)
! calculates the self-energy on the fly for q-point in the generated kmesh and Lorentzian delta
! frequencies, temperatures are in cm^-1 (same unit as wmax)
! It uses the following formula for q in the IFBZ (no restriction on k in the sum):
! sigma(q,lambda,w)=sum_k |v3(q,la;k,la1;-k-q,la2)|2 * (n1+n2+1/... + n1-n2/...)
! it assumes the input momentum q is on the already-generated kmesh and
! is the second argument of V33(k1,q,k2) where k2 covers the FBZ and k1=-q-k2
 use kpoints
 use ios
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
 use ios
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
! read v33_split
! calculate collision matrix
! save it to files
subroutine Pmatrix (ksubset,nk,kp,ndn,eigenval,temp)

 use kpoints
 use ios
 use phi3
 use lattice
! use geometry
! use atoms_force_constants
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
 real(8) q1(3), q2(3), q3(3), q1n(3), q3n(3), q3_proc1(3), q3_proc2(3)
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
               write(ulog,*) 'ijk2=',i2,j2,k2
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
 use ios
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
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
 use ios
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
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
 use ios
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
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

!! v33 = cmplx(0d0,0d0)
! v33sq = 0d0
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
    q1 = kp(:,n1)  ! should be = kp(:,mapinv(n1))
    write(*,*) 'in V3 cal, n1=',n1
    write(ulog,*) 'n1=', n1
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
!    if(mapinv(n1).ne.nk1) then
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


call matrix_elt_simplified(q1,q2,l1,l2,l3,xx,inside,ndn)    ! same as matrix_elt but efficiency is improved.

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


!===========================================================
 subroutine calculate_w3_ibz_split_sq(ibz_subset,nv3_split,ndn,nk,kp,eival,eivec,nv3s)!,v33sq_8)
! v33 (q1,q2,q3,la1,la2,la3), ksubset is in q1, q2 is all kpoints in FBZ, q3 is from momentum conservation 
! Full BZ calculated for iterative solution

! This is a modified calculate_w3_ibz by sy.
! It reads fbz_subset (fbz_subset(1)-starting index, ibz_subset(2)-end index in nibz)
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
 use ios
 use phi3
 use lattice
 use geometry
 use atoms_force_constants
 use svd_stuff
 use params
 use constants
 use phi3_sy
 implicit none
 integer nk,ndn,l1,l2,l3,n2,n3,j,k,n1
 integer i1,i2,i3,j1,j2,j3,k1,k2,k3,inside
 integer nk1,nk2,nk3,indx
 integer ibz_subset(2),nv3_split
 real(8) q1(3),q2(3),q3(3),eival(ndn,nk)
 real(8) kp(3,nk)
 complex(8) eivec(ndn,ndn,nk),xx
 real cputim, cputim2, cputim3, cputim_matrix_elt
 real(8) v33sq1_temp!, v33sq_8((ibz_subset(2)-ibz_subset(1)+1)*nk*ndn*ndn*ndn)

 integer utemp
 integer nv3s,s
 utemp=100


 call cal_eiqr2(nk,ndn)    ! calculate eiqr2, eivec2

 !write(*,*)'exited cal_eiqr2...'

!! v33 = cmplx(0d0,0d0)
! v33sq = 0d0
 indx=0
 write(ulog,*)'entered v3 subroutine...'
 write(ulog,*)'V3: nc(123)=',nc
 write(ulog,*)'starting fbz point=',ibz_subset(1)
 write(ulog,*)'ending fbz point=', ibz_subset(2)

!***   open(uv3,file='v33sq.dat',status='unknown',FORM='UNFORMATTED')
!   open(uv3,file='v33sq.dat',status='unknown')

   nv3=(ibz_subset(2)-ibz_subset(1)+1)*nk*ndn**3
!***   write(uv3) nv3
!   write(uv3,*) nv3
nv3s=nv3
v33s8=0
open(utemp,file='timeanal.dat',status='unknown')
call cpu_time(cputim2)
write(utemp,*) 'before loop1',cputim2
s=0
!************************************************
! second argument q2 needs to be only in the IFBZ
!************************************************
 loop1: do n1=ibz_subset(1),ibz_subset(2)       ! calculate v33 only for subset of k points
s=s+1
    q1 = kp(:,mapinv(n1))  ! should be = kp(:,mapinv(n1))
    !write(*,*) 'in V3 cal, n1=',n1
    write(ulog,*) 'n1=', n1
    call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
      if(mapinv(n1).ne.nk1) then
!     if(n1.ne.nk1) then
      write(ulog,*)'n1,nk1,inside=',n1,nk1,inside
      write(ulog,*)'q1=',q1
      write(ulog,*)'ijk1=',i1,j1,k1
      stop
    endif
!write(*,*)'first loop after checking mapinv(n1).e.nk1'
 call cpu_time(cputim)  !date_and_time(time=tim)
 write(ulog,*) 'entering q1...',n1,q1
 write(ulog,'(a,f12.4)')'  TIME IS ',cputim
 write(utemp,*) 'entering q1...',n1,q1
 write(utemp,'(a,f12.4)')'  TIME IS ',cputim

 cputim_matrix_elt=0.0


do l1=1,ndn

 loop2: do n2=1,nk
    q2 = kp(:,n2)   ! third argument is in the whole FBZ coarse mesh
    !write(debug,*) 'q3=', q3
    call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
    if(n2.ne.nk2) then
      write(ulog,*)'n2,nk2,inside=',n2,nk2,inside
      write(ulog,*)'q2=',q2
      write(ulog,*)'ijk2=',i2,j2,k2
      stop
    endif
!write(*,*)'second k-loop after checking mapinv(n1).e.nk1'
    q3 = -q2-q1

   call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)


 ! write(ulog,2)'nk123=',nk1,nk2,nk3
  !write(ulog,*)'q2,q3,q1',q2,q3,q1
 do l2=1 ,ndn
 do l3=1 ,ndn

call cpu_time(cputim2)


call matrix_elt_simplified(q1,q2,l1,l2,l3,xx,inside,ndn)    ! same as matrix_elt but efficiency is improved.

call cpu_time(cputim3)
cputim_matrix_elt=cputim_matrix_elt + (cputim3-cputim2)

v33sq1_temp=xx*conjg(xx)
!write(*,*)xx,v33sq1_temp
!s=s+1
!v33sq_8(s) = v33sq1_temp
  v33s8(s,n2,l1,l2,l3)=v33sq1_temp !v33sq_8(s)
!***write(*,*)s,n2,l1,l2,l3,v33sq1_temp

!***write(uv3) v33sq1_temp
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

 if (nv3 .ne. (ibz_subset(2)-ibz_subset(1)+1)*nk*ndn**3) then
   write(ulog,*) 'nv3 does not match in v3 calculation'
   stop
 endif

write(ulog,*)' V33: total size of this array is =',nv3


!***close(uv3)

2 format(a,3(3x,i5),9(2x,g11.5))
3 format(a,2(3x,3i2,1x,i5),9(2x,g11.5))
4 format(a,2(i5),3(2x,3(1x,i2)))
5 format(5(i6),9(2x,g13.7))
6 format(a,6(i5),9(2x,g14.8))
7 format(6(i6),9(2x,g14.8))
8 format(5(i6),2(2x,g14.8))
9 format(99(2x,g14.3))
 end subroutine calculate_w3_ibz_split_sq



!=============================================================================
subroutine calculate_RTA_distributionfunction(nk,nkk,kp,ndn,tempk,veloc,eigenval)!,FRTA)
! calculate F value using only diagonal terms in collision matrix (equivalent as RTA)

use constants
use exactBTE2
use params
use lattice
use kpoints, only : mapinv

implicit none

integer i,ii,iii,indx,nk,nkk,ndn
real(8) nbe,temp,tempk,nb,omega!,FRTA(nk,ndn,3)
real(8) veloc(3,ndn,nkk), eigenval(ndn,nk), kp(3,nk)

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

nkF1_loop: do i=1,nk  ! it's done in the IBZ
    
    laF1_loop: do ii=1,ndn

          !omega=sqrt(eigenval(ii,i)) * cnst
         
          !nb=nbe(omega,temp,classical)
          omega=frequency(mapinv(i),ii)
          nb=dist(mapinv(i),ii)        
         
          xyz_loop: do iii=1,3

                FFF(i,ii,iii)= -(c_light*veloc(iii,ii,mapinv(i))) * h_plank * (omega*c_light*100) * &
&                                nb * (nb+1) / (k_b*tempk**2) / (Qvalue(i,ii)*c_light*100)

              !write(*,*) i,ii,iii,FFF
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

end subroutine calculate_RTA_distributionfunction

!============================================================================
subroutine mode_kappa(ndn,tempk,veloc,eigenval)!,FFF) !!Saf!!
!***subroutine mode_thermal_iter2(ndn,tempk,veloc,eigenval)!,FFF) !!Saf!!
! calculate mode thermal conductivity for RTA and full solution of BTE

use constants
use exactBTE2
use params
use lattice
use ios
use kpoints 
implicit none

integer i,ii,j,jj,indx,ndn
integer s, narms,kvecop(48),n0,i3,j3,k3, inside
real(8) nbe, temp, tempk, nb0, omega,const1,const2
real(8) veloc(3,ndn,nkc), eigenval(ndn,nkc)
real(8) tot_kap,FRTA(3),q(3)!,FFF(nibz,ndn,3)
real(8) kvecstar(3,48),F12(3)
real(8) velocc(3)
integer n0i
integer ukapiter
ukapiter = 9400
tot_kap=0

kappa=0
!kappa_k=0
!kappa_RTA=0
!kappa_k_RTA=0
!diff_kap=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1
const1 = h_plank *c_light *100* c_light / (nkc*volume_r/1d30)
const2 = k_b*tempk**2/(h_plank*c_light*100)/c_light



write(ulog,*) 'entering mode_thermal_iter...'

nkF1_loop: do i=1,nibz

             q=kibz(:,i)

             call kstarsinfbz1(q,nkc,kpc,narms,kvecstar,kvecop)

laF1_loop: do ii=1,ndn
           tauinv_eff=0
           do s=1,narms

          !*** call syop_f(kvecop(s),F_RTA(i,ii,:),FRTA)
            call syop_f(kvecop(s),FFF(i,ii,:),F12)
            call syop_f(kvecop(s),veloc(:,ii,mapinv(i)),velocc)
           !***call syop_f(kvecop(s),F1(i,ii,:),F12)

          call get_k_info(kvecstar(:,s),NC,n0,i3,j3,k3,g1,g2,g3,inside)
          
          omega=sqrt(eigenval(ii,n0)) * cnst

          nb0=nbe(omega,temp,classical)

          alpha_loop: do j=1,3
          beta_loop: do jj=1,3

               
               kappa(ii,j,jj)=kappa(ii,j,jj) - const1 * omega * velocc(j) * nb0 * (nb0+1) * F12(jj)

             !***  kappa_RTA(ii,j,jj)=kappa_RTA(ii,j,jj) - const1 * omega * velocc(j) * nbb * (nbb+1) * FRTA(jj)


          enddo beta_loop
            if (velocc(j) .eq. 0.0) then       ! Gamma point has zero velocity
                  tauinv_eff(i,ii)= tauinv_eff(i,ii)-const2*F12(j) / omega / 0.01
            else
                  tauinv_eff(i,ii)= tauinv_eff(i,ii)-const2*F12(j) / omega / velocc(j)  /(narms*3d0)
            end if

          enddo alpha_loop
         enddo
       enddo laF1_loop
!enddo
enddo nkF1_loop



! open(ukapiter+1,file='kappa_svd.dat',status='unknown')
!***open(ukapiter,file='tauinv_eff.dat',status='unknown')
!***write(ukapiter,*) 'nibz, la, tainv_eff'
!***, kap , kap_RTA

!***tot_kap=0
!***do ii=1,ndn
!**do i=1,nk
!***   do i=1,nibz
!***    write(ukapiter,8)i, ii,tauinv_eff(i,ii)
!***enddo
!***enddo
!*** write(ukapiter,3)(kappa(ii,j,j),j=1,3),((kappa(ii,j,jj),jj=j+1,3),j=1,2)
!***   do j=1,3
!***      tot_kap=tot_kap+kappa(ii,j,j)/3
!***   enddo
!*** enddo
!***   write(*,*) 'tot_kap=', tot_kap

!***close(ukapiter)




3 format(99(1x,g11.4))
8 format(99(1x,g11.4))
end subroutine mode_kappa
!=====================================================================================



!============================================================================
subroutine mode_kappa_d(ndn,tempk,veloc,eigenval,Fi_x) !!Saf!!
!subroutine mode_thermal_noniter23(nki,nk,kibbz,kp,ndn,tempk,veloc,eigenval,Fi_x) !!Saf!!
! calculate mode thermal conductivity for RTA and full solution of BTE

use constants
use exactBTE2
use params
use lattice
use ios
use kpoints 
implicit none

integer i,ii,j,jj,indx,ndn
integer s, narms,kvecop(48),n0,ni,i3,j3,k3, inside
real(8) nbe, temp, tempk, nb0, omega ,const1,const2
real(8) veloc(3,ndn,nkc), eigenval(ndn,nkc)
real(8) tot_kap,FRTA(3),q(3)
real(8) kvecstar(3,48),F12(3)
real(8) Fi_x(nibz*ndn*3), velocc(3), veloci(3,ndn,nkc)
integer i1,j1,j2,r
integer ukapiter
ukapiter = 9400
tot_kap=0

kappa=0
!kappa_k=0
!kappa_RTA=0
!kappa_k_RTA=0
!diff_kap=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1
const1 = h_plank *c_light *100* c_light / (nkc*volume_r/1d30)
const2 = k_b*tempk**2/(h_plank*c_light*100)/c_light

!write(debug,*) 'i,ii,j,jj,om,veloc_j,nb,F2(indx,jj)'

write(ulog,*) 'entering mode_thermal_iter...'


FFF=0
r=0
do i1=1,nibz
   do j1=1,ndn
      do j2=1,3
!         if (s.lt.kl) then
         r=r+1
         FFF(i1,j1,j2)=Fi_x(r)
!         endif
       enddo
    enddo
enddo

call veloc_avergeIBZ_invsyFBZ(nkc,kpc,nibz,kibz,ndn,veloc,veloci)

nkF1_loop: do i=1,nibz

             q=kibz(:,i)
!            call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
             call kstarsinfbz1(q,nkc,kpc,narms,kvecstar,kvecop)
!             do s=1,narms
!                call syop_f(kvecop(s),F1(,l2,:),Fsq2)


      laF1_loop: do ii=1,ndn
           tauinv_eff=0
           do s=1,narms
!        call syop_f(kvecop(s),F_RTA(mapinv(i),ii,:),FRTA)
      !***   call syop_f(kvecop(s),F_RTA(i,ii,:),FRTA)
        ! laF1_loop: do ii=1,ndn
          call syop_f(kvecop(s),FFF(i,ii,:),F12)

          call syop_f(kvecop(s),veloci(:,ii,i),velocc)

          call get_k_info(kvecstar(:,s),NC,n0,i3,j3,k3,g1,g2,g3,inside)
          !!omega=sqrt(eigenval(ii,i)) * cnst !!mine command
          omega=sqrt(eigenval(ii,n0)) * cnst

          nb0=nbe(omega,temp,classical)

          alpha_loop: do j=1,3
          beta_loop: do jj=1,3

           kappa    (ii,j,jj)=kappa    (ii,j,jj) - const1 * omega * velocc(j) * nb0 * (nb0+1) * F12(jj)
!***           kappa_RTA(ii,j,jj)=kappa_RTA(ii,j,jj) - const1 * omega * velocc(j) * nbb * (nbb+1) * FRTA(jj)

          enddo beta_loop

            if (velocc(j) .eq. 0.0) then       ! Gamma point has zero velocity
                 tauinv_eff(i,ii)= tauinv_eff(i,ii)-const2*F12(j) / omega / 0.01
           else
                 tauinv_eff(i,ii)= tauinv_eff(i,ii)-const2*F12(j) / omega / velocc(j)  /(narms*3d0)
           end if

          enddo alpha_loop
         enddo
       enddo laF1_loop
enddo nkF1_loop


! open(ukapiter+1,file='kappa_svd.dat',status='unknown')
!***@ open(ukapiter,file='tauinv_eff.dat',status='unknown')
!***@ write(ukapiter,*) 'nibz, la, tainv_eff' 
!***, kap , kap_RTA

!***@tot_kap=0
!***@do ii=1,ndn
!**do i=1,nk
!***@   do i=1,nibz
!***@   write(ukapiter,8)i, ii,tauinv_eff(i,ii)
!***@enddo
!***@enddo
!***   write(ukapiter,3)(kappa(ii,j,j),j=1,3),((kappa(ii,j,jj),jj=j+1,3),j=1,2)
!***   do j=1,3
!***      tot_kap=tot_kap+kappa(ii,j,j)/3
!***   enddo
!*** enddo 
!***   write(*,*) 'tot_kap=', tot_kap

!***@close(ukapiter)
!close(ukapiter+1)

3 format(99(1x,g11.4))
8 format(2i5,99(1x,g11.4))
end subroutine mode_kappa_d
!=====================================================================================



!============================================================================
subroutine mode_thermal_noniter2(nk,ndn,tempk,veloc,eigenval,F_MRHS2)
! calculate mode thermal conductivity for RTA and full solution of BTE

use constants
use exactBTE2
use params
use lattice
use ios

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
!               tauinv_eff(i,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*FF1(i,ii,j)/(0.01*c_light)
          else
!               tauinv_eff(i,ii,j)=-(k_b*tempk**2)/(h_plank*(omega*c_light*100))*FF1(i,ii,j)/(veloc(j,ii,i)*c_light)
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


 subroutine read_v3_new_unformatted_sq(nki,nk,ndn)

! This is modified version of read_v3_new_unformatted.
! It will use multi-dimensional array for v3

 use phi3
 use phi3_sy
 use ios
 use lattice
 implicit none
 integer j,unt
 real(8) xx,yy
 integer nki,nk,ndn

 integer n1,n2,n3,l1,l2,l3

 unt = 111
! open(unt,file='v33sq.dat',status='old',form='UNFORMATTED')
open(unt,file='v33sq.dat',status='old')

! read(unt,end=99) nv3
 read(unt,*,end=99) nv3
 call allocate_v33_sq1(nki,nk,ndn,ndn,ndn)
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



!=========================================================================================
subroutine write_thermal(iteration,nk,ndn,tempk,veloc,eigenval,kp)

use constants
use exactBTE2
use params
use lattice
use ios

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

write(ukapiter+1,*) ' '
write(ukapiter+1,*) 'tot_kap='
write(ukapiter+1,*) ' '

write(ukapiter+1,78) (sum(kappa(:,j,j)),j=1,3) !,((sum(kappa(:,j,ajj)),jj=j+1,3),j=1,2) &
! &                    (sum(kappa_RTA(:,j,j)),j=1,3),((sum(kappa_RTA(:,j,ajj)),jj=j+1,3),j=1,2)

78 format(2(3x,2(3(1x,f9.3),2x)))

close(ukapiter)
close(ukapiter+1)


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
 subroutine write_kappa(ndn,tempk,veloc,eigenval,ukapiter,uRTA) !!Saf!!
!subroutine write_thermal_noniter2(ndn,tempk,veloc,eigenval) !!Saf!!
! nk = nibz
use constants
use exactBTE2
use params
use lattice
use ios
use kpoints
implicit none

integer i,ii,j,jj,indx,ndn,iteration

integer ukapiter,uRTA
real(8) tempk, nb0, omega ,tot_kap
real(8) veloc(3,ndn,nkc), eigenval(ndn,nkc)

!ukapiter = 9401
!uRTA=9405
write(ulog,*) 'entering write_thermal...'


!***open(uRTA,file='iself_RTA_noniter.dat',status='replace')
!***open(uRTA,file='tau_RTAandeff.dat',status='replace')
!***open(ukapiter,file='mode_kappa_noniter.dat',status='replace')
!***open(ukapiter+1,file='kappa_noniter.dat',status='unknown')

!***write(ukapiter,*) 'nk, la, al, be, mode_kap, mode_kap_RTA'
!write(ukapiter+1,*) 'iteration',iteration
write(ukapiter,*) 'la, kxx,kyy,kzz, kxy,kxz,kyz, k_RTA'

tot_kap=0
do j=1,3
do jj=1,3
do ii=1,ndn
!do i=1,nk
!write(ukapiter,3) i,ii,j,jj,kappa_k(i,ii,j,jj),kappa_k_RTA(i,ii,j,jj)
!enddo

 if (j.eq.jj) then
    tot_kap=tot_kap+kappa(ii,j,jj)/3
 endif

enddo
enddo
enddo

do ii=1,ndn
write(ukapiter,8) ii,(kappa(ii,j,j),j=1,3),((kappa(ii,j,jj),jj=j+1,3),j=1,2)!,(kappa_RTA(ii,j,j),j=1,3),((kappa_RTA(ii,j,jj),jj=j+1,3),j=1,2)
enddo
write(ukapiter,*) ' '
write(ukapiter,*) 'tot_kap=' ,tot_kap
write(ukapiter,*) ' '

write(ukapiter,78) (sum(kappa(:,j,j)),j=1,3) !,((sum(kappa(:,j,ajj)),jj=j+1,3),j=1,2) &
! &                    (sum(kappa_RTA(:,j,j)),j=1,3),((sum(kappa_RTA(:,j,ajj)),jj=j+1,3),j=1,2)
write(*,*) (sum(kappa(:,j,j)),j=1,3)
78 format(2(3x,2(3(1x,f9.3),2x)))

!***close(ukapiter)
!***close(ukapiter+1)

write(uRTA,*)'tempk,i,kibz(:,i),ii,omega,nb0,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),1/tauinv_eff(i,ii),veloc(:,ii,mapinv(i))'

nkF1_loop: do i=1,nibz
    laF1_loop: do ii=1,ndn
          omega=frequency(mapinv(i),ii)
          nb0=dist(mapinv(i),ii)

           if(tauinv_U(i,ii) .eq. 0) then
           tauinv_U(i,ii)=1.0d10
           endif
           if(tauinv_N(i,ii) .eq. 0) then
           tauinv_N(i,ii)=1.0d10
           endif

           write(uRTA,10) tempk,i,kibz(:,i),ii,omega,nb0,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),1/tauinv_eff(i,ii),veloc(:,ii,mapinv(i))!,(kappa_k(i,ii,j,j),j=1,3)
!,(kappa_k_RTA(i,ii,j,j),j=1,3),(kappa_k(i,ii,j,j),j=1,3)
!          write(uRTA,10)
!          tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i),(kappa_k_RTA(i,ii,j,j),j=1,3),(kappa_k(i,ii,j,j),j=1,3),&
!&         (tauinv_eff(i,ii,j),j=1,3)

!!-- temporary
if (tauinv_U(i,ii) .eq. 0.0) then
write(ulog,*) 'zero tauinv_U, i, ii', i, ii
endif

    enddo laF1_loop
enddo nkF1_loop

!***close(uRTA)


10 FORMAT(99(1x,g14.8))
3 format(99(1x,g11.4))
8 format(i2,99(1x,g11.4))


end subroutine write_kappa
!===================================================








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

!!==========================================================
subroutine get_k_infos(q,nq,nk,kp,ins)

use lattice 
!use kpoints
use params
use ios

implicit none
integer, intent(in) :: nk 
integer, intent(out) :: nq 
integer i 
real(8)  lk
real(8), intent(in) :: q(3),kp(3,nk)
logical, intent(out):: ins

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
!==========================================================

!==========================================================
subroutine get_kvecopn2(q,k,nk,kp,spf) !!Saf!!
!! find kvecop number for q3_proc1,2
use constants
use exactBTE2
use lattice
use params
use kpoints
use ios
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
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
implicit none

integer s,ier
real(8) fs(3), ff(3), tempi(3,3)
 

  fs(:) = MATMUL(op_kmatrix(:,:,s),ff(:))


end subroutine syop_f
!==========================================================

!===========================================================
subroutine update(ndn)!!Saf!!
!subroutine cal_F2m23(nk,kp,nibbz,kibbz,ndn,tempk) !!Saf!!
! calculate F value based on F value from previous iteration.
! This subroutine uses symmetry operations to calculate F in FBZ through F in IBZ.
! Summation over k1 and k2 in IBZ, but in the case of k2, F(k2 in FBZ) = S*F(k in IBZ)
! Finally, F2 are calculated in IBZ 


use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
implicit none

integer ndn 
integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
integer nk1,nk2,nk3, nk3_proc1, nk3_proc2,xyz
real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3), q3(3), qq3_proc1(3)
real(8) q3_proc1_i(3), q3_proc2_i(3)


real(8) Avalue(nibz,ndn,3)
!real(8) Avalue(nibbz,ndn,3)  
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
!temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1


call cpu_time(cputim1)
cputim1_wall = wtime()
write(ulog,*) 'entering cal_F2m23...'


indx=0
Avalue=0

loop1 : do n1=1,nibz

           q1=kibz(:,n1)
 
           call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)


       do l1=1,ndn
    loopxyz : do xyz=1,3

       loop2: do n2=1,nibz

                q2=kibz(:,n2)
                n2i0=mapinv(n2)
                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                call kstarsinfbz1(q2,nkc,kpc,narms,kvecstar,kvecop)

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
                   q3_proc1_i=kibz(:,nk3_proc1_i)
                   q3_proc2_i=kibz(:,nk3_proc2_i)

                   !kvecop for q2
                   sq2= kvecop(s)
                   ! find kvecop number for q3_proc1,2
                   call get_kvecopn2(q3_proc1_i,q3_proc1,nkc,kpc,sp1f)
                   call get_kvecopn2(q3_proc2_i,q3_proc2,nkc,kpc,sp2f)

 
                  do l2=1,ndn
                   do l3=1,ndn

                    !F(FBZ)=S*F(IBZ)
!                    call syop_f(sq2,F1(n2i0,l2,:),Fsq2)
!                    call syop_f(sp1f,F1(mapinv(nk3_proc1_i),l3,:),Fsp1)
!                    call syop_f(sp2f,F1(mapinv(nk3_proc2_i),l3,:),Fsp2)
!changetoIBZ
                    call syop_f(sq2,FFF(n2,l2,:),Fsq2)
                    call syop_f(sp1f,FFF(nk3_proc1_i,l3,:),Fsp1)
                    call syop_f(sp2f,FFF(nk3_proc2_i,l3,:),Fsp2)



!changetoIBZ       Avalue(mapinv(n1),l1,xyz)=Avalue(mapinv(n1),l1,xyz)+P1(mapinv(n1),ns,l1,l2,l3)*(Fsp1(xyz)-Fsq2(xyz))+0.5*P2(mapinv(n1),ns,l1,l2,l3)*(Fsq2(xyz)+Fsp2(xyz))
                   Avalue(n1,l1,xyz)=Avalue(n1,l1,xyz)+P1(n1,ns,l1,l2,l3)*(Fsp1(xyz)-Fsq2(xyz))+0.5*P2(n1,ns,l1,l2,l3)*(Fsq2(xyz)+Fsp2(xyz))
                  enddo
                 enddo
           enddo !loop s(narms)
        enddo loop2

!changetoIBZ   F2(mapinv(n1),l1,xyz)=F_RTA(mapinv(n1),l1,xyz) + Avalue(mapinv(n1),l1,xyz)/Qvalue(mapinv(n1),l1)
               F2(n1,l1,xyz)=F_RTA(n1,l1,xyz) + Avalue(n1,l1,xyz)/Qvalue(n1,l1)

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

end subroutine update
!==================================      


!===========================================================
subroutine update_iso(ndn)!!Saf!!
!subroutine cal_F2m23(nk,kp,nibbz,kibbz,ndn,tempk) !!Saf!!
! calculate F value based on F value from previous iteration.
! This subroutine uses symmetry operations to calculate F in FBZ through F in
! IBZ.
! Summation over k1 and k2 in IBZ, but in the case of k2, F(k2 in FBZ) = S*F(k
! in IBZ)
! Finally, F2 are calculated in IBZ


use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
implicit none

integer ndn
integer i1,j1,k1, i2,j2,k2, i3,j3,k3, inside
integer nk1,nk2,nk3, nk3_proc1, nk3_proc2,xyz
real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3), q3(3), qq3_proc1(3)
real(8) q3_proc1_i(3), q3_proc2_i(3)


real(8) Avalue(nibz,ndn,3)
!real(8) Avalue(nibbz,ndn,3)
real(8) kvecstar(3,48), kvecstarp1(3,48), kvecstarp2(3,48)
real(8) kvecstar1(3,48), q10(3)
real(8) Fsp1(3), Fsp2(3), Fsq2(3)
real xp1,xp2, ssigma

integer n1,n2,n3,l1,l2,l3, nk2n
integer indx,nk_subset
integer ksub_size2, ksubset(2)
integer narms, s,sq2,ss1,ss2,sp1,sp2,sp1f,sp2f,narmsp1, narmsp2,ns,kvecop(48),kvecopp1(48),kvecopp2(48)
integer kvecop1(48),narms1,s1
integer nk3_proc1_i, nk3_proc2_i, ns_star_i
integer nq,n2i0
logical ins
integer ms1, ms2


real cputim1,cputim2
real(8) cputim1_wall, cputim2_wall, wtime

real(8) temp, tempk
!temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1


call cpu_time(cputim1)
cputim1_wall = wtime()
write(ulog,*) 'entering cal_F2m23...'


indx=0
Avalue=0

loop1 : do n1=1,nibz

           q1=kibz(:,n1)

           call get_k_info(q1,NC,nk1,i1,j1,k1,g1,g2,g3,inside)


       do l1=1,ndn
    loopxyz : do xyz=1,3

       loop2: do n2=1,nibz

                q2=kibz(:,n2)
                n2i0=mapinv(n2)
                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                call kstarsinfbz1(q2,nkc,kpc,narms,kvecstar,kvecop)

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
                   q3_proc1_i=kibz(:,nk3_proc1_i)
                   q3_proc2_i=kibz(:,nk3_proc2_i)

                   !kvecop for q2
                   sq2= kvecop(s)
                   ! find kvecop number for q3_proc1,2
                   call get_kvecopn2(q3_proc1_i,q3_proc1,nkc,kpc,sp1f)
                   call get_kvecopn2(q3_proc2_i,q3_proc2,nkc,kpc,sp2f)


                  do l2=1,ndn
                   do l3=1,ndn

                    !F(FBZ)=S*F(IBZ)
!                    call syop_f(sq2,F1(n2i0,l2,:),Fsq2)
!                    call syop_f(sp1f,F1(mapinv(nk3_proc1_i),l3,:),Fsp1)
!                    call syop_f(sp2f,F1(mapinv(nk3_proc2_i),l3,:),Fsp2)
!changetoIBZ
                    call syop_f(sq2,FFF(n2,l2,:),Fsq2)
                    call syop_f(sp1f,FFF(nk3_proc1_i,l3,:),Fsp1)
                    call syop_f(sp2f,FFF(nk3_proc2_i,l3,:),Fsp2)



!changetoIBZ
!Avalue(mapinv(n1),l1,xyz)=Avalue(mapinv(n1),l1,xyz)+P1(mapinv(n1),ns,l1,l2,l3)*(Fsp1(xyz)-Fsq2(xyz))+0.5*P2(mapinv(n1),ns,l1,l2,l3)*(Fsq2(xyz)+Fsp2(xyz))
                   Avalue(n1,l1,xyz)=Avalue(n1,l1,xyz)+P1(n1,ns,l1,l2,l3)*(Fsp1(xyz)-Fsq2(xyz))+0.5*P2(n1,ns,l1,l2,l3)*(Fsq2(xyz)+Fsp2(xyz))
                  enddo
                 ! Avalue(n1,l1,xyz)=Avalue(n1,l1,xyz)+Piso(n1,ns,l1,l2)*Fsq2(xyz)
                 enddo
           enddo !loop s(narms)
        enddo loop2
               Avalue(n1,l1,xyz)=Avalue(n1,l1,xyz)+Piso(n1,l1)*Fsq2(xyz)
!changetoIBZ   F2(mapinv(n1),l1,xyz)=F_RTA(mapinv(n1),l1,xyz) + Avalue(mapinv(n1),l1,xyz)/Qvalue(mapinv(n1),l1)
               F2(n1,l1,xyz)=F_RTA(n1,l1,xyz) + Avalue(n1,l1,xyz)/Qvalue(n1,l1)

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

end subroutine update_iso
!==================================





!=============================================================================
subroutine cal_RHS_IBZ(ndn,nkk,tempk,veloc,eigenval)!,RHSS) !!Saf!!
! (equivalent as RTA) 
! Right hand side of the linearized BTE in just IBZ.

use constants
use exactBTE2
use params
use lattice
use kpoints 

implicit none

integer i,ii,iii,indx,ndn,s,kl,j,n0,nkk
real(8) nbe,temp,tempk,nb0,omega
real(8) veloc(3,ndn,nkk), eigenval(ndn,nibz)


temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

open(322,file='RHSS_IBZ.dat')

kl=nibz*ndn*3
s=0
nkF1_loop: do i=1,nibz

          n0 = mapinv(i)

    laF1_loop: do ii=1,ndn

          !omega=sqrt(eigenval(ii,i)) * cnst


          !nb=nbe(omega,temp,classical)
          omega=frequency(n0,ii)
          nb0=dist(n0,ii)

          xyz_loop: do iii=1,3
              !  do s=1,kl
                s=s+1
                !RHSS(s)= -(c_light*veloc(iii,ii,n0)) * h_plank * (omega*c_light*100) * nbb * (nbb+1) / (k_b*tempk**2) / (c_light*100)
                RHS(s)= (c_light*veloc(iii,ii,n0)) * h_plank * (omega*c_light*100) * nb0 * (nb0+1) / (k_b*tempk**2) / (c_light*100)
                ! veloc in c0, h_plank in J*s, omega in cm^-1 (linear freq.),
                ! k_b in J/K, tempk in K, Qvalue in cm^-1
                ! resulting F1 is in m/K

!write(*,*) 'indx,F1',indx,F1(indx,iii)
!                tau(i,ii,iii)= nb*(nb+1)/(Qvalue(i,ii))      ! tau is in cm

!                write(uRTA,10)
!                tempk,i,kp(:,i),ii,omega,nb,1/tauinv_N(i,ii),1/tauinv_U(i,ii),1/tauinv_tot(i,ii),veloc(:,ii,i)
               !s=s+1
                write(322,*)s,i,n0,ii,iii,RHS(s)
              !enddo
          enddo xyz_loop

    enddo laF1_loop

enddo nkF1_loop

!F1(:,:,:)=F_RTA(:,:,:) ! for iteration

close(322)

open(412,file='RHS_IBZ.dat')
do j=1,kl

  write(412,*)RHS(j)
enddo
close(412)




end subroutine cal_RHS_IBZ
!===========================================================


!==========================================================
!subroutine sycollisionM(nibbz,kibbz,ndn,nk,kp,Cs) !!Saf!!
subroutine sycollisionM_split2(ksbst,nibbz,kibbz,ndn,nk,kp,Cs) !!Saf!!
!Collision matrix(IBZ*FBZ)*Symmerty matrix(FBZ*IBZ)=Cs(IBZ*IBZ) :output

use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
implicit none

integer n1,l1,n2,l2,ns, s1,s2,ksbst(2)
integer ndn,nibbz,narms,nk
integer i,j,kil,kil1,kil2,ii,jj,d
real(8), intent(in) :: kp(3,nk),kibbz(3,nibbz)
real(8) q(3),q1(3),q2(3),q3_p1(3),q3_p2(3)
integer i3,j3,k3,inside,nk3_p1
real(8) Cs((ksbst(2)-ksbst(1)+1)*ndn*3,nibbz*ndn*3)
real(8) C_d, C_of, C
real(8) kvecstars(3,48)
integer kvecops(48),n0,n22,nl2,l3,nl3,ni0,n1i,d1,d2
s1=0
s2=0
Cs=0
C=0
! kil=nibbz*ndn*3
kil1=(ksbst(2)-ksbst(1)+1)*ndn*3
kil2=nibbz*ndn*3

open(701,file='sycmatrix.dat')

!   do n1=1,nibbz    
    do n1=ksbst(1),ksbst(2)
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
                        
! K1                    C = C + (-P1(n1i,n22,l1,nl2,l3) - 0.5*P2(n1i,n22,l1,nl2,l3))
                        C = C + (-P1(n1 ,n22,l1,nl2,l3) - 0.5*P2(n1 ,n22,l1,nl2,l3))

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

! K1              C = C + (-P1(n1i,n0,l1,l2,nl3) + P1(n1i,nk3_p1,l1,l2,nl3) + P2(n1i,n0,l1,l2,nl3))
                  C = C + (-P1(n1 ,n0,l1,l2,nl3) + P1(n1 ,nk3_p1,l1,l2,nl3) + P2(n1 ,n0,l1,l2,nl3))
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
!! K1  write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,C_d,C_of,op_kmatrix(i,j,kvecops(ns)),kvecops(ns), Cs(s1,s2)
                write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,op_kmatrix(i,j,kvecops(ns)),kvecops(ns), Cs(s1,s2)
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
do ii=1,kil1
   do jj=1,kil2
      write(601,*) ii,jj, Cs(ii,jj)

   enddo
enddo
close(601)



end subroutine sycollisionM_split2

!==========================================================
subroutine CollisionMatrix(ndn)!(nibbz,kibbz,ndn,nk,kp)!,Cs) !!Saf!!
!subroutine sycollisionM(ksbst,nibbz,kibbz,ndn,nk,kp,Cs) !!Saf!!
!Collision matrix(IBZ*FBZ)*Symmerty matrix(FBZ*IBZ)=Cs(IBZ*IBZ) :output

use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
implicit none

integer n1,l1,n2,l2,ns, s1,s2
integer ndn,narms
integer i,j,kil,ii,jj,d
!***real(8), intent(in) :: 
real(8) q(3),q1(3),q2(3),q3_p1(3),q3_p2(3)
integer i3,j3,k3,inside,nk3_p1
!***real(8) Cs(nibbz*ndn*3,nibbz*ndn*3)
real(8) C_d, C_of, C
real(8) kvecstars(3,48)
integer kvecops(48),n0,n22,nl2,l3,nl3,ni0,n1i,d1,d2
s1=0
s2=0
!***Cs=0
C=0
Collision_Matrix=0

 kil=nibz*ndn*3

open(701,file='sycmatrix.dat')

    do n1=1,nibz    
       n1i=mapinv(n1)
    do l1=1,ndn
      s1=s1+1
      d1=s1

      do n2=1,nibz
           q=kibz(:,n2)
           ni0=mapinv(n2)
           call kstarsinfbz1(q,nkc,kpc,narms,kvecstars,kvecops)
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
                    do n22=1,nkc
                        
! K1                    C = C + (-P1(n1i,n22,l1,nl2,l3) - 0.5*P2(n1i,n22,l1,nl2,l3))
                        C = C + (-P1(n1 ,n22,l1,nl2,l3) - 0.5*P2(n1 ,n22,l1,nl2,l3))

                    enddo  
                   enddo  
                enddo
            ! C2 = C_d
            ! C3 = C_d
    
           else
           
         q1 = kibz(:,n1)
         q2 = kvecstars(:,ns)
         q3_p1 = q1+q2
         q3_p2 = q1-q2
   
        ! write(*,*)'here', n1,n1i,l1,n2,ns,narms,n0,l2    
  
         call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
               C=0
      loopnl3:  do nl3=1,ndn

! K1              C = C + (-P1(n1i,n0,l1,l2,nl3) + P1(n1i,nk3_p1,l1,l2,nl3) + P2(n1i,n0,l1,l2,nl3))
                  C = C + (-P1(n1 ,n0,l1,l2,nl3) + P1(n1 ,nk3_p1,l1,l2,nl3) + P2(n1 ,n0,l1,l2,nl3))
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
                Collision_Matrix(s1,s2) = Collision_Matrix(s1,s2) + (C * op_kmatrix(i,j,kvecops(ns)))
                !write(*,*)'n1,n1i,l1,n2,l2,ni0,narms,ns,n0,s1,s2, C'
!! K1  write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,C_d,C_of,op_kmatrix(i,j,kvecops(ns)),kvecops(ns), Cs(s1,s2)
                write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,op_kmatrix(i,j,kvecops(ns)),kvecops(ns), Collision_Matrix(s1,s2)
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
      write(601,*) ii,jj, Collision_Matrix(ii,jj)
   enddo
enddo
close(601)



end subroutine CollisionMatrix
!===========================================================
!==========================================================
subroutine CollisionMatrix_mpi(ksbst,ndn)!Cs) !!Saf!!
!subroutine sycollisionM(ksbst,nibbz,kibbz,ndn,nk,kp,Cs) !!Saf!!
!Collision matrix(IBZ*FBZ)*Symmerty matrix(FBZ*IBZ)=Cs(IBZ*IBZ) :output

use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
implicit none

integer n1,l1,n2,l2,ns, s1,s2
integer ndn,narms
integer i,j,kil,ii,jj,d,kil1,kil2
real(8) q(3),q1(3),q2(3),q3_p1(3),q3_p2(3)
integer i3,j3,k3,inside,nk3_p1
real(8) C_d, C_of, C
real(8) kvecstars(3,48)
integer kvecops(48),n0,n22,nl2,l3,nl3,ni0,n1i,d1,d2
integer ksbst(2),k_shift,ksubsize3

s1=0
s2=0
!Cs=0
C=0
! kil=nibbz*ndn*3
Collision_Matrixmpi=0
!!!!Collision_Matrixmpi
!kil1=((ksbst(2)-ksbst(1))+1)*ndn*3
!kil2=nibz*ndn*3


ksubsize3=(ksbst(2)-ksbst(1))+1
k_shift=ksbst(1)-1

!nkI1_loop: do ibz=1,ksubsize3

!open(701,file='sycmatrix.dat')
 !** do n1=ksbst(1),ksbst(2)
   do n1=1,ksubsize3
!    do n1=1,nibbz
       n1i=mapinv(n1+k_shift)
    do l1=1,ndn
      s1=s1+1
      d1=s1

      do n2=1,nibz
           q=kibz(:,n2)
           ni0=mapinv(n2)
           call kstarsinfbz1(q,nkc,kpc,narms,kvecstars,kvecops)
 !     write(701,*) '**************',n2,narms
       do l2=1,ndn
          s2=s2+1
          d2=s2
         !C_d=0
         !C_of=0
 !      write(701,*)'n1,  n1i,  l1,  n2,  ni0,  n0,  l2,  narms,  ns,  s1,  s2,C_d,  C_of,   Cs(s1,s2)'
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
                    do n22=1,nkc

! K1                    C = C + (-P1(n1i,n22,l1,nl2,l3) - 0.5*P2(n1i,n22,l1,nl2,l3))
                        C = C + (-P1(n1 ,n22,l1,nl2,l3) - 0.5*P2(n1,n22,l1,nl2,l3))

                    enddo
                   enddo
                enddo
            ! C2 = C_d
            ! C3 = C_d

           else

         q1 = kibz(:,n1+k_shift)
         q2 = kvecstars(:,ns)
         q3_p1 = q1+q2
         q3_p2 = q1-q2

        ! write(*,*)'here', n1,n1i,l1,n2,ns,narms,n0,l2

         call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
               C=0
      loopnl3:  do nl3=1,ndn

! K1              C = C + (-P1(n1i,n0,l1,l2,nl3) + P1(n1i,nk3_p1,l1,l2,nl3) +  P2(n1i,n0,l1,l2,nl3))
                  C = C + (-P1(n1 ,n0,l1,l2,nl3) + P1(n1 ,nk3_p1,l1,l2,nl3) + P2(n1 ,n0,l1,l2,nl3))
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
                Collision_Matrixmpi(s1,s2) = Collision_Matrixmpi(s1,s2) + (C * op_kmatrix(i,j,kvecops(ns)))
                !write(*,*)'n1,n1i,l1,n2,l2,ni0,narms,ns,n0,s1,s2, C'
!! K1 write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,C_d,C_of,op_kmatrix(i,j,kvecops(ns)),kvecops(ns),Cs(s1,s2)
        !        write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,op_kmatrix(i,j,kvecops(ns)),kvecops(ns), Collision_Matrix(s1,s2)
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

!close(701)

!open(601,file='Cs_sycollisionm.dat',status='unknown',form='unformatted')
!do ii=1,kil1
!   do jj=1,kil2
!      write(601) Collision_Matrix(ii,jj)
!   enddo
!enddo
!close(601)



end subroutine CollisionMatrix_mpi
!===========================================================


!==========================================================
subroutine CollisionMatrix_split(ksbst,ndn)!Cs) !!Saf!!
!subroutine sycollisionM(ksbst,nibbz,kibbz,ndn,nk,kp,Cs) !!Saf!!
!Collision matrix(IBZ*FBZ)*Symmerty matrix(FBZ*IBZ)=Cs(IBZ*IBZ) :output

use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
implicit none

integer n1,l1,n2,l2,ns, s1,s2
integer ndn,narms
integer i,j,kil,ii,jj,d,kil1,kil2
real(8) q(3),q1(3),q2(3),q3_p1(3),q3_p2(3)
integer i3,j3,k3,inside,nk3_p1
real(8) C_d, C_of, C
real(8) kvecstars(3,48)
integer kvecops(48),n0,n22,nl2,l3,nl3,ni0,n1i,d1,d2
integer ksbst(2)

s1=0
s2=0
!Cs=0
C=0
! kil=nibbz*ndn*3
Collision_Matrix=0
kil1=((ksbst(2)-ksbst(1))+1)*ndn*3
kil2=nibz*ndn*3

open(701,file='sycmatrix.dat')
  do n1=ksbst(1),ksbst(2)
!    do n1=1,nibbz
       n1i=mapinv(n1)
    do l1=1,ndn
      s1=s1+1
      d1=s1

      do n2=1,nibz
           q=kibz(:,n2)
           ni0=mapinv(n2)
           call kstarsinfbz1(q,nkc,kpc,narms,kvecstars,kvecops)
      write(701,*) '**************',n2,narms
       do l2=1,ndn
          s2=s2+1
          d2=s2
         !C_d=0
         !C_of=0
       write(701,*)'n1,  n1i,  l1,  n2,  ni0,  n0,  l2,  narms,  ns,  s1,  s2,C_d,  C_of,   Cs(s1,s2)'
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
                    do n22=1,nkc

! K1                    C = C + (-P1(n1i,n22,l1,nl2,l3) - 0.5*P2(n1i,n22,l1,nl2,l3))
                        C = C + (-P1(n1 ,n22,l1,nl2,l3) - 0.5*P2(n1,n22,l1,nl2,l3))

                    enddo
                   enddo
                enddo
            ! C2 = C_d
            ! C3 = C_d

           else

         q1 = kibz(:,n1)
         q2 = kvecstars(:,ns)
         q3_p1 = q1+q2
         q3_p2 = q1-q2

        ! write(*,*)'here', n1,n1i,l1,n2,ns,narms,n0,l2

         call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
               C=0
      loopnl3:  do nl3=1,ndn

! K1              C = C + (-P1(n1i,n0,l1,l2,nl3) + P1(n1i,nk3_p1,l1,l2,nl3) +  P2(n1i,n0,l1,l2,nl3))
                  C = C + (-P1(n1 ,n0,l1,l2,nl3) + P1(n1 ,nk3_p1,l1,l2,nl3) + P2(n1 ,n0,l1,l2,nl3))
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
                Collision_Matrix(s1,s2) = Collision_Matrix(s1,s2) + (C * op_kmatrix(i,j,kvecops(ns)))
                !write(*,*)'n1,n1i,l1,n2,l2,ni0,narms,ns,n0,s1,s2, C'
!! K1 write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,C_d,C_of,op_kmatrix(i,j,kvecops(ns)),kvecops(ns),Cs(s1,s2)
                write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,op_kmatrix(i,j,kvecops(ns)),kvecops(ns), Collision_Matrix(s1,s2)
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

open(601,file='Cs_sycollisionm.dat',status='unknown',form='unformatted')
do ii=1,kil1
   do jj=1,kil2
      write(601) Collision_Matrix(ii,jj)
   enddo
enddo
close(601)



end subroutine CollisionMatrix_split
!===========================================================



!==========================================================
subroutine CollisionMatrix_iso(ndn)!(nibbz,kibbz,ndn,nk,kp)!,Cs) !!Saf!!
!subroutine sycollisionM(ksbst,nibbz,kibbz,ndn,nk,kp,Cs) !!Saf!!
!Collision matrix(IBZ*FBZ)*Symmerty matrix(FBZ*IBZ)=Cs(IBZ*IBZ) :output

use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
implicit none

integer n1,l1,n2,l2,ns, s1,s2
integer ndn,narms
integer i,j,kil,ii,jj,d
!***real(8), intent(in) ::
real(8) q(3),q1(3),q2(3),q3_p1(3),q3_p2(3)
integer i3,j3,k3,inside,nk3_p1
!***real(8) Cs(nibbz*ndn*3,nibbz*ndn*3)
real(8) C_d, C_of, C
real(8) kvecstars(3,48)
integer kvecops(48),n0,n22,nl2,l3,nl3,ni0,n1i,d1,d2
s1=0
s2=0
!***Cs=0
C=0
Collision_Matrix=0

 kil=nibz*ndn*3

open(701,file='sycmatrix.dat')

    do n1=1,nibz
       n1i=mapinv(n1)
    do l1=1,ndn
      s1=s1+1
      d1=s1

      do n2=1,nibz
           q=kibz(:,n2)
           ni0=mapinv(n2)
           call kstarsinfbz1(q,nkc,kpc,narms,kvecstars,kvecops)
      write(701,*) '**************',n2,narms
       do l2=1,ndn
          s2=s2+1
          d2=s2
         !C_d=0
         !C_of=0
       write(701,*)'n1,  n1i,  l1,  n2,  ni0,  n0,  l2,  narms,  ns,  s1,  s2, C_d,  C_of,   Cs(s1,s2)'
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
                    do n22=1,nkc

! K1                    C = C + (-P1(n1i,n22,l1,nl2,l3) -  0.5*P2(n1i,n22,l1,nl2,l3))
                        C = C + (-P1(n1 ,n22,l1,nl2,l3) - 0.5*P2(n1,n22,l1,nl2,l3))!+Piso(n1,n22,l1,nl2))

                    enddo
                   enddo
                enddo
            ! C2 = C_d
            ! C3 = C_d

           !     do nl2=1,ndn
           !       do n22=1,nkc
                     
                     C = C - Piso(n1,l1)

            !      enddo
            !    enddo

           else

         q1 = kibz(:,n1)
         q2 = kvecstars(:,ns)
         q3_p1 = q1+q2
         q3_p2 = q1-q2

        ! write(*,*)'here', n1,n1i,l1,n2,ns,narms,n0,l2

         call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
               C=0
      loopnl3:  do nl3=1,ndn

! K1              C = C + (-P1(n1i,n0,l1,l2,nl3) + P1(n1i,nk3_p1,l1,l2,nl3) +  P2(n1i,n0,l1,l2,nl3))
                  C = C + (-P1(n1 ,n0,l1,l2,nl3) + P1(n1 ,nk3_p1,l1,l2,nl3) + P2(n1 ,n0,l1,l2,nl3))
                enddo loopnl3
         !    C4 = C_of
         !    C5 = C_of

           !     do nl2=1,ndn
            !      do n22=1,nkc

                     C = C + Piso(n1,l1)

             !     enddo
              !  enddo


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
                Collision_Matrix(s1,s2) = Collision_Matrix(s1,s2) + (C * op_kmatrix(i,j,kvecops(ns)))
                !write(*,*)'n1,n1i,l1,n2,l2,ni0,narms,ns,n0,s1,s2, C'
!! K1 write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,C_d,C_of,op_kmatrix(i,j,kvecops(ns)),kvecops(ns), Cs(s1,s2)^M
                write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,op_kmatrix(i,j,kvecops(ns)),kvecops(ns), Collision_Matrix(s1,s2)
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
      write(601,*) ii,jj, Collision_Matrix(ii,jj)
   enddo
enddo
close(601)



end subroutine CollisionMatrix_iso
!===========================================================



!==========================================================
subroutine CollisionMatrix_iso_split(ndn)!(nibbz,kibbz,ndn,nk,kp)!,Cs) !!Saf!!
!subroutine sycollisionM(ksbst,nibbz,kibbz,ndn,nk,kp,Cs) !!Saf!!
!Collision matrix(IBZ*FBZ)*Symmerty matrix(FBZ*IBZ)=Cs(IBZ*IBZ) :output

use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
implicit none

integer n1,l1,n2,l2,ns, s1,s2
integer ndn,narms
integer i,j,kil,ii,jj,d
!***real(8), intent(in) ::
real(8) q(3),q1(3),q2(3),q3_p1(3),q3_p2(3)
integer i3,j3,k3,inside,nk3_p1
!***real(8) Cs(nibbz*ndn*3,nibbz*ndn*3)
real(8) C_d, C_of, C
real(8) kvecstars(3,48)
integer kvecops(48),n0,n22,nl2,l3,nl3,ni0,n1i,d1,d2
s1=0
s2=0
!***Cs=0
C=0
Collision_Matrix=0

 kil=nibz*ndn*3

open(701,file='sycmatrix.dat')

    do n1=1,nibz
       n1i=mapinv(n1)
    do l1=1,ndn
      s1=s1+1
      d1=s1

      do n2=1,nibz
           q=kibz(:,n2)
           ni0=mapinv(n2)
           call kstarsinfbz1(q,nkc,kpc,narms,kvecstars,kvecops)
      write(701,*) '**************',n2,narms
       do l2=1,ndn
          s2=s2+1
          d2=s2
         !C_d=0
         !C_of=0
       write(701,*)'n1,  n1i,  l1,  n2,  ni0,  n0,  l2,  narms,  ns,  s1,  s2, C_d,  C_of,   Cs(s1,s2)'
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
                    do n22=1,nkc

! K1                    C = C + (-P1(n1i,n22,l1,nl2,l3) - 0.5*P2(n1i,n22,l1,nl2,l3))
                        C = C + (-P1(n1 ,n22,l1,nl2,l3) - 0.5*P2(n1,n22,l1,nl2,l3))!+Piso(n1,n22,l1,nl2))

                    enddo
                   enddo
                enddo
            ! C2 = C_d
            ! C3 = C_d

                !do nl2=1,ndn
                !  do n22=1,nkc

                     C = C - Piso(n1,l1)

                 ! enddo
                !enddo

           else

         q1 = kibz(:,n1)
         q2 = kvecstars(:,ns)
         q3_p1 = q1+q2
         q3_p2 = q1-q2

        ! write(*,*)'here', n1,n1i,l1,n2,ns,narms,n0,l2

         call get_k_info(q3_p1,NC,nk3_p1,i3,j3,k3,g1,g2,g3,inside)
               C=0
      loopnl3:  do nl3=1,ndn

! K1              C = C + (-P1(n1i,n0,l1,l2,nl3) + P1(n1i,nk3_p1,l1,l2,nl3) + P2(n1i,n0,l1,l2,nl3))
                  C = C + (-P1(n1 ,n0,l1,l2,nl3) + P1(n1 ,nk3_p1,l1,l2,nl3) + P2(n1 ,n0,l1,l2,nl3))
                enddo loopnl3
         !    C4 = C_of
         !    C5 = C_of

               ! do nl2=1,ndn
                !  do n22=1,nkc

                     C = C + Piso(n1,l1)

                 ! enddo
                !enddo


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
                Collision_Matrix(s1,s2) = Collision_Matrix(s1,s2) + (C * op_kmatrix(i,j,kvecops(ns)))
                !write(*,*)'n1,n1i,l1,n2,l2,ni0,narms,ns,n0,s1,s2, C'
!! K1
!write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,C_d,C_of,op_kmatrix(i,j,kvecops(ns)),kvecops(ns), Cs(s1,s2)^M
                write(701,*)n1,n1i,l1,n2,ni0,n0,l2,narms,ns,s1,s2,op_kmatrix(i,j,kvecops(ns)),kvecops(ns),Collision_Matrix(s1,s2)
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
      write(601,*) ii,jj, Collision_Matrix(ii,jj)
   enddo
enddo
close(601)



end subroutine CollisionMatrix_iso_split
!===========================================================



!===========================================================
subroutine cal_g(attname,g)

 use ios
 use params
 use lattice
 use atoms_force_constants


implicit none
! character(99) :: Ge
 character*2 :: H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,&
             Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe,Cs,Ba,La,&
                 Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,&
                 At,Rn,Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,&
                 Rg,Cn,Uuq,Uuh


 !character*2 isoname(287)
!character isoname(287)
character*2, allocatable:: isoname(:)
character*2 attname
real(8) isomass(287), isof(287), mave,g !geatom(natom_type),g
integer numiso,nn,ii,isosize
isosize=287
allocate(isoname(isosize) )

    isoname(1)= 'Ag'
    isomass(1)=106.905095d0
    isof(1)=51.84d0
    isoname(2)= 'Ag'
    isomass(2)=108.904754d0
    isof(2)=48.16d0
    isoname(3)= 'Al'
    isomass(3)=26.981541d0
    isof(3)=100.0d0
    isoname(4)= 'Ar'
    isomass(4)=35.967546d0
    isof(4)=0.34d0
    isoname(5)= 'Ar'
    isomass(5)=37.962732d0
    isof(5)=0.063d0
    isoname(6)= 'Ar'
    isomass(6)=39.962383d0
    isof(6)=99.6d0
    isoname(7)= 'As'
    isomass(7)=74.921596d0
    isof(7)=100.0d0
    isoname(8)= 'Au'
    isomass(8)=196.96656d0
    isof(8)=100.0d0
    isoname(9)= 'B'
    isomass(9)=10.012938d0
    isof(9)=19.8d0
    isoname(10)= 'B'
    isomass(10)=11.009305d0
    isof(10)=80.2d0
    isoname(11)= 'Ba'
    isomass(11)=129.906277d0
    isof(11)=0.11d0
    isoname(12)= 'Ba'
    isomass(12)=131.905042d0
    isof(12)=0.1d0
    isoname(13)= 'Ba'
    isomass(13)=133.90449d0
    isof(13)=2.42d0
    isoname(14)= 'Ba'
    isomass(14)=134.905668d0
    isof(14)=6.59d0
    isoname(15)= 'Ba'
    isomass(15)=135.904556d0
    isof(15)=7.85d0
    isoname(16)= 'Ba'
    isomass(16)=136.905816d0
    isof(16)=11.23d0
    isoname(17)= 'Ba'
    isomass(17)=137.905236d0
    isof(17)=71.7d0
    isoname(18)= 'Be'
    isomass(18)=9.012183d0
    isof(18)=100.0d0
    isoname(19)= 'Bi'
    isomass(19)=208.980388d0
    isof(19)=100.0d0
    isoname(20)= 'Br'
    isomass(20)=78.918336d0
    isof(20)=50.69d0
    isoname(21)= 'Br'
    isomass(21)=80.91629d0
    isof(21)=49.31d0
    isoname(22)= 'C'
    isomass(22)=12.0d0
    isof(22)=98.9d0
    isoname(23)= 'C'
    isomass(23)=13.003355d0
    isof(23)=1.1d0
    isoname(24)= 'Ca'
    isomass(24)=39.962591d0
    isof(24)=96.95d0
    isoname(25)= 'Ca'
    isomass(25)=41.958622d0
    isof(25)=0.65d0
    isoname(26)= 'Ca'
    isomass(26)=42.95877d0
    isof(26)=0.14d0
    isoname(27)= 'Ca'
    isomass(27)=43.955485d0
    isof(27)=2.086d0
    isoname(28)= 'Ca'
    isomass(28)=45.953689d0
    isof(28)=0.004d0
    isoname(29)= 'Ca'
    isomass(29)=47.952532d0
    isof(29)=0.19d0
    isoname(30)= 'Cd'
    isomass(30)=105.906461d0
    isof(30)=1.25d0
    isoname(31)= 'Cd'
    isomass(31)=107.904186d0
    isof(31)=0.89d0
    isoname(32)= 'Cd'
    isomass(32)=109.903007d0
    isof(32)=12.49d0
    isoname(33)= 'Cd'
    isomass(33)=110.904182d0
    isof(33)=12.8d0
    isoname(34)= 'Cd'
    isomass(34)=111.902761d0
    isof(34)=24.13d0
    isoname(35)= 'Cd'
    isomass(35)=112.904401d0
    isof(35)=12.22d0
    isoname(36)= 'Cd'
    isomass(36)=113.903361d0
    isof(36)=28.73d0
    isoname(37)= 'Cd'
    isomass(37)=115.904758d0
    isof(37)=7.49d0
    isoname(38)= 'Ce'
    isomass(38)=135.90714d0
    isof(38)=0.19d0
    isoname(39)= 'Ce'
    isomass(39)=137.905996d0
    isof(39)=0.25d0
    isoname(40)= 'Ce'
    isomass(40)=139.905442d0
    isof(40)=88.48d0
    isoname(41)= 'Ce'
    isomass(41)=141.909249d0
    isof(41)=11.08d0
    isoname(42)= 'Cl'
    isomass(42)=34.968853d0
    isof(42)=75.77d0
    isoname(43)= 'Cl'
    isomass(43)=36.965903d0
    isof(43)=24.23d0
    isoname(44)= 'Co'
    isomass(44)=58.933198d0
    isof(44)=100.0d0
    isoname(45)= 'Cr'
    isomass(45)=49.946046d0
    isof(45)=4.35d0
    isoname(46)= 'Cr'
    isomass(46)=51.94051d0
    isof(46)=83.79d0
    isoname(47)= 'Cr'
    isomass(47)=52.940651d0
    isof(47)=9.5d0
    isoname(48)= 'Cr'
    isomass(48)=53.938882d0
    isof(48)=2.36d0
    isoname(49)= 'Cs'
    isomass(49)=132.905433d0
    isof(49)=100.0d0
    isoname(50)= 'Cu'
    isomass(50)=62.929599d0
    isof(50)=69.17d0
    isoname(51)= 'Cu'
    isomass(51)=64.927792d0
    isof(51)=30.83d0
    isoname(52)= 'Dy'
    isomass(52)=155.924287d0
    isof(52)=0.06d0
    isoname(53)= 'Dy'
    isomass(53)=157.924412d0
    isof(53)=0.1d0
    isoname(54)= 'Dy'
    isomass(54)=159.925203d0
    isof(54)=2.34d0
    isoname(55)= 'Dy'
    isomass(55)=160.926939d0
    isof(55)=18.9d0
    isoname(56)= 'Dy'
    isomass(56)=161.926805d0
    isof(56)=25.5d0
    isoname(57)= 'Dy'
    isomass(57)=162.928737d0
    isof(57)=24.9d0
    isoname(58)= 'Dy'
    isomass(58)=163.929183d0
    isof(58)=28.2d0
    isoname(59)= 'Er'
    isomass(59)=161.928787d0
    isof(59)=0.14d0
    isoname(60)= 'Er'
    isomass(60)=163.929211d0
    isof(60)=1.61d0
    isoname(61)= 'Er'
    isomass(61)=165.930305d0
    isof(61)=33.6d0
    isoname(62)= 'Er'
    isomass(62)=166.932061d0
    isof(62)=22.95d0
    isoname(63)= 'Er'
    isomass(63)=167.932383d0
    isof(63)=26.8d0
    isoname(64)= 'Er'
    isomass(64)=169.935476d0
    isof(64)=14.9d0
    isoname(65)= 'Eu'
    isomass(65)=150.91986d0
    isof(65)=47.8d0
    isoname(66)= 'Eu'
    isomass(66)=152.921243d0
    isof(66)=52.2d0
    isoname(67)= 'F'
    isomass(67)=18.998403d0
    isof(67)=100.0d0
    isoname(68)= 'Fe'
    isomass(68)=53.939612d0
    isof(68)=5.8d0
    isoname(69)= 'Fe'
    isomass(69)=55.934939d0
    isof(69)=91.72d0
    isoname(70)= 'Fe'
    isomass(70)=56.935396d0
    isof(70)=2.2d0
    isoname(71)= 'Fe'
    isomass(71)=57.933278d0
    isof(71)=0.28d0
    isoname(72)= 'Ga'
    isomass(72)=68.925581d0
    isof(72)=60.1d0
    isoname(73)= 'Ga'
    isomass(73)=70.924701d0
    isof(73)=39.9d0
    isoname(74)= 'Gd'
    isomass(74)=151.919803d0
    isof(74)=0.2d0
    isoname(75)= 'Gd'
    isomass(75)=153.920876d0
    isof(75)=2.18d0
    isoname(76)= 'Gd'
    isomass(76)=154.822629d0
    isof(76)=14.8d0
    isoname(77)= 'Gd'
    isomass(77)=155.92213d0
    isof(77)=20.47d0
    isoname(78)= 'Gd'
    isomass(78)=156.923967d0
    isof(78)=15.65d0
    isoname(79)= 'Gd'
    isomass(79)=157.924111d0
    isof(79)=24.84d0
    isoname(80)= 'Gd'
    isomass(80)=159.927061d0
    isof(80)=21.86d0
    isoname(81)='Ge'
    isomass(81)=69.92425d0
    isof(81)=20.5d0
    isoname(82)='Ge'
    isomass(82)=71.92208d0
    isof(82)=27.4d0
    isoname(83)='Ge'
    isomass(83)=72.923464d0
    isof(83)=7.8d0
    isoname(84)='Ge'
    isomass(84)=73.921179d0
    isof(84)=36.5d0
    isoname(85)='Ge'
    isomass(85)=75.921403d0
    isof(85)=7.8d0
    isoname(86)= 'H'
    isomass(86)=1.007825d0
    isof(86)=99.99d0
    isoname(87)= 'H'
    isomass(87)=2.014102d0
    isof(87)=0.015d0
    isoname(88)= 'He'
    isomass(88)=3.016029d0
    isof(88)=0.0001d0
    isoname(89)= 'He'
    isomass(89)=4.002603d0
    isof(89)=100.0d0
    isoname(90)= 'Hf'
    isomass(90)=173.940065d0
    isof(90)=0.16d0
    isoname(91)= 'Hf'
    isomass(91)=175.94142d0
    isof(91)=5.2d0
    isoname(92)= 'Hf'
    isomass(92)=176.943233d0
    isof(92)=18.6d0
    isoname(93)= 'Hf'
    isomass(93)=177.94371d0
    isof(93)=27.1d0
    isoname(94)= 'Hf'
    isomass(94)=178.945827d0
    isof(94)=13.74d0
    isoname(95)= 'Hf'
    isomass(95)=179.946561d0
    isof(95)=35.2d0
    isoname(96)= 'Hg'
    isomass(96)=195.965812d0
    isof(96)=0.15d0
    isoname(97)= 'Hg'
    isomass(97)=197.96676d0
    isof(97)=10.1d0
    isoname(98)= 'Hg'
    isomass(98)=198.968269d0
    isof(98)=17.0d0
    isoname(99)= 'Hg'
    isomass(99)=199.968316d0
    isof(99)=23.1d0
    isoname(100)= 'Hg'
    isomass(100)=200.970293d0
    isof(100)=13.2d0
    isoname(101)= 'Hg'
    isomass(101)=201.970632d0
    isof(101)=29.65d0
    isoname(102)= 'Hg'
    isomass(102)=203.973481d0
    isof(102)=6.8d0
    isoname(103)= 'Ho'
    isomass(103)=164.930332d0
    isof(103)=100.0d0
    isoname(104)= 'I'
    isomass(104)=126.904477d0
    isof(104)=100.0d0
    isoname(105)= 'In'
    isomass(105)=112.904056d0
    isof(105)=4.3d0
    isoname(106)= 'In'
    isomass(106)=114.903875d0
    isof(106)=95.7d0
    isoname(107)= 'Ir'
    isomass(107)=190.960603d0
    isof(107)=37.3d0
    isoname(108)= 'Ir'
    isomass(108)=192.962942d0
    isof(108)=62.7d0
    isoname(109)= 'K'
    isomass(109)=38.963708d0
    isof(109)=93.2d0
    isoname(110)= 'K'
    isomass(110)=39.963999d0
    isof(110)=0.012d0
    isoname(111)= 'K'
    isomass(111)=40.961825d0
    isof(111)=6.73d0
    isoname(112)= 'Kr'
    isomass(112)=77.920397d0
    isof(112)=0.35d0
    isoname(113)= 'Kr'
    isomass(113)=79.916375d0
    isof(113)=2.25d0
    isoname(114)= 'Kr'
    isomass(114)=81.913483d0
    isof(114)=11.6d0
    isoname(115)= 'Kr'
    isomass(115)=82.914134d0
    isof(115)=11.5d0
    isoname(116)= 'Kr'
    isomass(116)=83.911506d0
    isof(116)=57.0d0
    isoname(117)= 'Kr'
    isomass(117)=85.910614d0
    isof(117)=17.3d0
    isoname(118)= 'La'
    isomass(118)=137.907114d0
    isof(118)=0.09d0
    isoname(119)= 'La'
    isomass(119)=138.906355d0
    isof(119)=99.91d0
    isoname(120)= 'Li'
    isomass(120)=6.015123d0
    isof(120)=7.42d0
    isoname(121)= 'Li'
    isomass(121)=7.016005d0
    isof(121)=92.58d0
    isoname(122)= 'Lu'
    isomass(122)=174.940785d0
    isof(122)=97.4d0
    isoname(123)= 'Lu'
    isomass(123)=175.942694d0
    isof(123)=2.6d0
    isoname(124)= 'Mg'
    isomass(124)=23.985045d0
    isof(124)=78.9d0
    isoname(125)= 'Mg'
    isomass(125)=24.985839d0
    isof(125)=10.0d0
    isoname(126)= 'Mg'
    isomass(126)=25.982595d0
    isof(126)=11.1d0
    isoname(127)= 'Mn'
    isomass(127)=54.938046d0
    isof(127)=100.0d0
    isoname(128)= 'Mo'
    isomass(128)=91.906809d0
    isof(128)=14.84d0
    isoname(129)= 'Mo'
    isomass(129)=93.905086d0
    isof(129)=9.25d0
    isoname(130)= 'Mo'
    isomass(130)=94.905838d0
    isof(130)=15.92d0
    isoname(131)= 'Mo'
    isomass(131)=95.904676d0
    isof(131)=16.68d0
    isoname(132)= 'Mo'
    isomass(132)=96.906018d0
    isof(132)=9.55d0
    isoname(133)= 'Mo'
    isomass(133)=97.905405d0
    isof(133)=24.13d0
    isoname(134)= 'Mo'
    isomass(134)=99.907473d0
    isof(134)=9.63d0
    isoname(135)= 'N'
    isomass(135)=14.003074d0
    isof(135)=99.63d0
    isoname(136)= 'N'
    isomass(136)=15.000109d0
    isof(136)=0.37d0
    isoname(137)= 'Na'
    isomass(137)=22.98977d0
    isof(137)=100.0d0
    isoname(138)= 'Nb'
    isomass(138)=92.906378d0
    isof(138)=100.0d0
    isoname(139)= 'Nd'
    isomass(139)=141.907731d0
    isof(139)=27.13d0
    isoname(140)= 'Nd'
    isomass(140)=142.909823d0
    isof(140)=12.18d0
    isoname(141)= 'Nd'
    isomass(141)=143.910096d0
    isof(141)=23.8d0
    isoname(142)= 'Nd'
    isomass(142)=144.912582d0
    isof(142)=8.3d0
    isoname(143)= 'Nd'
    isomass(143)=145.913126d0
    isof(143)=17.19d0
    isoname(144)= 'Nd'
    isomass(144)=147.916901d0
    isof(144)=5.76d0
    isoname(145)= 'Nd'
    isomass(145)=149.9209d0
    isof(145)=5.64d0
    isoname(146)= 'Ne'
    isomass(146)=19.992439d0
    isof(146)=90.6d0
    isoname(147)= 'Ne'
    isomass(147)=20.993845d0
    isof(147)=0.26d0
    isoname(148)= 'Ne'
    isomass(148)=21.991384d0
    isof(148)=9.2d0
    isoname(149)= 'Ni'
    isomass(149)=57.935347d0
    isof(149)=68.27d0
    isoname(150)= 'Ni'
    isomass(150)=59.930789d0
    isof(150)=26.1d0
    isoname(151)= 'Ni'
    isomass(151)=60.931059d0
    isof(151)=1.13d0
    isoname(152)= 'Ni'
    isomass(152)=61.928346d0
    isof(152)=3.59d0
    isoname(153)= 'Ni'
    isomass(153)=63.927968d0
    isof(153)=0.91d0
    isoname(154)= 'O'
    isomass(154)=15.994915d0
    isof(154)=99.76d0
    isoname(155)= 'O'
    isomass(155)=16.999131d0
    isof(155)=0.038d0
    isoname(156)= 'O'
    isomass(156)=17.999159d0
    isof(156)=0.2d0
    isoname(157)= 'Os'
    isomass(157)=183.952514d0
    isof(157)=0.02d0
    isoname(158)= 'Os'
    isomass(158)=185.953852d0
    isof(158)=1.58d0
    isoname(159)= 'Os'
    isomass(159)=186.955762d0
    isof(159)=1.6d0
    isoname(160)= 'Os'
    isomass(160)=187.95585d0
    isof(160)=13.3d0
    isoname(161)= 'Os'
    isomass(161)=188.958156d0
    isof(161)=16.1d0
    isoname(162)= 'Os'
    isomass(162)=189.958455d0
    isof(162)=26.4d0
    isoname(163)= 'Os'
    isomass(163)=191.961487d0
    isof(163)=41.0d0
    isoname(164)= 'P'
    isomass(164)=30.973763d0
    isof(164)=100.0d0
    isoname(165)= 'Pb'
    isomass(165)=203.973037d0
    isof(165)=1.4d0
    isoname(166)= 'Pb'
    isomass(166)=205.974455d0
    isof(166)=24.1d0
    isoname(167)= 'Pb'
    isomass(167)=206.975885d0
    isof(167)=22.1d0
    isoname(168)= 'Pb'
    isomass(168)=207.976641d0
    isof(168)=52.4d0
    isoname(169)= 'Pd'
    isomass(169)=101.905609d0
    isof(169)=1.02d0
    isoname(170)= 'Pd'
    isomass(170)=103.904026d0
    isof(170)=11.14d0
    isoname(171)= 'Pd'
    isomass(171)=104.905075d0
    isof(171)=22.33d0
    isoname(172)= 'Pd'
    isomass(172)=105.903475d0
    isof(172)=27.33d0
    isoname(173)= 'Pd'
    isomass(173)=107.903894d0
    isof(173)=26.46d0
    isoname(174)= 'Pd'
    isomass(174)=109.905169d0
    isof(174)=11.72d0
    isoname(175)= 'Pr'
    isomass(175)=140.907657d0
    isof(175)=100.0d0
    isoname(176)= 'Pt'
    isomass(176)=189.959937d0
    isof(176)=0.01d0
    isoname(177)= 'Pt'
    isomass(177)=191.961049d0
    isof(177)=0.79d0
    isoname(178)= 'Pt'
    isomass(178)=193.962679d0
    isof(178)=32.9d0
    isoname(179)= 'Pt'
    isomass(179)=194.964785d0
    isof(179)=33.8d0
    isoname(180)= 'Pt'
    isomass(180)=195.964947d0
    isof(180)=25.3d0
    isoname(181)= 'Pt'
    isomass(181)=197.967879d0
    isof(181)=7.2d0
    isoname(182)= 'Rb'
    isomass(182)=84.9118d0
    isof(182)=72.17d0
    isoname(183)= 'Rb'
    isomass(183)=86.909184d0
    isof(183)=27.84d0
    isoname(184)= 'Re'
    isomass(184)=184.952977d0
    isof(184)=37.4d0
    isoname(185)= 'Re'
    isomass(185)=186.955765d0
    isof(185)=62.6d0
    isoname(186)= 'Rh'
    isomass(186)=102.905503d0
    isof(186)=100.0d0
    isoname(187)= 'Ru'
    isomass(187)=95.907596d0
    isof(187)=5.52d0
    isoname(188)= 'Ru'
    isomass(188)=97.905287d0
    isof(188)=1.88d0
    isoname(189)= 'Ru'
    isomass(189)=98.905937d0
    isof(189)=12.7d0
    isoname(190)= 'Ru'
    isomass(190)=99.904218d0
    isof(190)=12.6d0
    isoname(191)= 'Ru'
    isomass(191)=100.905581d0
    isof(191)=17.0d0
    isoname(192)= 'Ru'
    isomass(192)=101.904348d0
    isof(192)=31.6d0
    isoname(193)= 'Ru'
    isomass(193)=103.905422d0
    isof(193)=18.7d0
    isoname(194)= 'S'
    isomass(194)=31.972072d0
    isof(194)=95.02d0
    isoname(195)= 'S'
    isomass(195)=32.971459d0
    isof(195)=0.75d0
    isoname(196)= 'S'
    isomass(196)=33.967868d0
    isof(196)=4.21d0
    isoname(197)= 'S'
    isomass(197)=35.967079d0
    isof(197)=0.02d0
    isoname(198)= 'Sb'
    isomass(198)=120.903824d0
    isof(198)=57.3d0
    isoname(199)= 'Sb'
    isomass(199)=122.904222d0
    isof(199)=42.7d0
    isoname(200)= 'Sc'
    isomass(200)=44.955914d0
    isof(200)=100.0d0
    isoname(201)= 'Se'
    isomass(201)=73.922477d0
    isof(201)=0.9d0
    isoname(202)= 'Se'
    isomass(202)=75.919207d0
    isof(202)=9.0d0
    isoname(203)= 'Se'
    isomass(203)=76.919908d0
    isof(203)=7.6d0
    isoname(204)= 'Se'
    isomass(204)=77.917304d0
    isof(204)=23.5d0
    isoname(205)= 'Se'
    isomass(205)=79.916521d0
    isof(205)=49.6d0
    isoname(206)= 'Se'
    isomass(206)=81.916709d0
    isof(206)=9.4d0
    isoname(207)= 'Si'
    isomass(207)=27.976928d0
    isof(207)=92.23d0
    isoname(208)= 'Si'
    isomass(208)=28.976496d0
    isof(208)=4.67d0
    isoname(209)= 'Si'
    isomass(209)=29.973772d0
    isof(209)=3.1d0
    isoname(210)= 'Sm'
    isomass(210)=143.912009d0
    isof(210)=3.1d0
    isoname(211)= 'Sm'
    isomass(211)=146.914907d0
    isof(211)=15.0d0
    isoname(212)= 'Sm'
    isomass(212)=147.914832d0
    isof(212)=11.3d0
    isoname(213)= 'Sm'
    isomass(213)=148.917193d0
    isof(213)=13.8d0
    isoname(214)= 'Sm'
    isomass(214)=149.917285d0
    isof(214)=7.4d0
    isoname(215)= 'Sm'
    isomass(215)=151.919741d0
    isof(215)=26.7d0
    isoname(216)= 'Sm'
    isomass(216)=153.922218d0
    isof(216)=22.7d0
    isoname(217)= 'Sn'
    isomass(217)=111.904826d0
    isof(217)=0.97d0
    isoname(218)= 'Sn'
    isomass(218)=113.902784d0
    isof(218)=0.65d0
    isoname(219)= 'Sn'
    isomass(219)=114.903348d0
    isof(219)=0.36d0
    isoname(220)= 'Sn'
    isomass(220)=115.901744d0
    isof(220)=14.7d0
    isoname(221)= 'Sn'
    isomass(221)=116.902954d0
    isof(221)=7.7d0
    isoname(222)= 'Sn'
    isomass(222)=117.901607d0
    isof(222)=24.3d0
    isoname(223)= 'Sn'
    isomass(223)=118.90331d0
    isof(223)=8.6d0
    isoname(224)= 'Sn'
    isomass(224)=119.902199d0
    isof(224)=32.4d0
    isoname(225)= 'Sn'
    isomass(225)=121.90344d0
    isof(225)=4.6d0
    isoname(226)= 'Sn'
    isomass(226)=123.905271d0
    isof(226)=5.6d0
    isoname(227)= 'Sr'
    isomass(227)=83.913428d0
    isof(227)=0.56d0
    isoname(228)= 'Sr'
    isomass(228)=85.909273d0
    isof(228)=9.86d0
    isoname(229)= 'Sr'
    isomass(229)=86.908902d0
    isof(229)=7.0d0
    isoname(230)= 'Sr'
    isomass(230)=87.905625d0
    isof(230)=82.58d0
    isoname(231)= 'Ta'
    isomass(231)=179.947489d0
    isof(231)=0.012d0
    isoname(232)= 'Ta'
    isomass(232)=180.948014d0
    isof(232)=99.99d0
    isoname(233)= 'Tb'
    isomass(233)=158.92535d0
    isof(233)=100.0d0
    isoname(234)= 'Te'
    isomass(234)=119.904021d0
    isof(234)=0.096d0
    isoname(235)= 'Te'
    isomass(235)=121.903055d0
    isof(235)=2.6d0
    isoname(236)= 'Te'
    isomass(236)=122.904278d0
    isof(236)=0.91d0
    isoname(237)= 'Te'
    isomass(237)=123.902825d0
    isof(237)=4.82d0
    isoname(238)= 'Te'
    isomass(238)=124.904435d0
    isof(238)=7.14d0
    isoname(239)= 'Te'
    isomass(239)=125.90331d0
    isof(239)=18.95d0
    isoname(240)= 'Te'
    isomass(240)=127.904464d0
    isof(240)=31.69d0
    isoname(241)= 'Te'
    isomass(241)=129.906229d0
    isof(241)=33.8d0
    isoname(242)= 'Th'
    isomass(242)=232.038054d0
    isof(242)=100.0d0
    isoname(243)= 'Ti'
    isomass(243)=45.952633d0
    isof(243)=8.0d0
    isoname(244)= 'Ti'
    isomass(244)=46.951765d0
    isof(244)=7.3d0
    isoname(245)= 'Ti'
    isomass(245)=47.947947d0
    isof(245)=73.8d0
    isoname(246)= 'Ti'
    isomass(246)=48.947871d0
    isof(246)=5.5d0
    isoname(247)= 'Ti'
    isomass(247)=49.944786d0
    isof(247)=5.4d0
    isoname(248)= 'Tl'
    isomass(248)=202.972336d0
    isof(248)=29.52d0
    isoname(249)= 'Tl'
    isomass(249)=204.97441d0
    isof(249)=70.48d0
    isoname(250)= 'Tm'
    isomass(250)=168.934225d0
    isof(250)=100.0d0
    isoname(251)= 'U'
    isomass(251)=234.040947d0
    isof(251)=0.006d0
    isoname(252)= 'U'
    isomass(252)=235.043925d0
    isof(252)=0.72d0
    isoname(253)= 'U'
    isomass(253)=238.050786d0
    isof(253)=99.27d0
    isoname(254)= 'V'
    isomass(254)=49.947161d0
    isof(254)=0.25d0
    isoname(255)= 'V'
    isomass(255)=50.943963d0
    isof(255)=99.75d0
    isoname(256)= 'W'
    isomass(256)=179.946727d0
    isof(256)=0.13d0
    isoname(257)= 'W'
    isomass(257)=181.948225d0
    isof(257)=26.3d0
    isoname(258)= 'W'
    isomass(258)=182.950245d0
    isof(258)=14.3d0
    isoname(259)= 'W'
    isomass(259)=183.950953d0
    isof(259)=30.67d0
    isoname(260)= 'W'
    isomass(260)=185.954377d0
    isof(260)=28.6d0
    isoname(261)= 'Xe'
    isomass(261)=123.905894d0
    isof(261)=0.1d0
    isoname(262)= 'Xe'
    isomass(262)=125.904281d0
    isof(262)=0.09d0
    isoname(263)= 'Xe'
    isomass(263)=127.903531d0
    isof(263)=1.91d0
    isoname(264)= 'Xe'
    isomass(264)=128.90478d0
    isof(264)=26.4d0
    isoname(265)= 'Xe'
    isomass(265)=129.90351d0
    isof(265)=4.1d0
    isoname(266)= 'Xe'
    isomass(266)=130.905076d0
    isof(266)=21.2d0
    isoname(267)= 'Xe'
    isomass(267)=131.904148d0
    isof(267)=26.9d0
    isoname(268)= 'Xe'
    isomass(268)=133.905395d0
    isof(268)=10.4d0
    isoname(269)= 'Xe'
    isomass(269)=135.907219d0
    isof(269)=8.9d0
    isoname(270)= 'Y'
    isomass(270)=88.905856d0
    isof(270)=100.0d0
    isoname(271)= 'Yb'
    isomass(271)=167.933908d0
    isof(271)=0.13d0
    isoname(272)= 'Yb'
    isomass(272)=169.934774d0
    isof(272)=3.05d0
    isoname(273)= 'Yb'
    isomass(273)=170.936338d0
    isof(273)=14.3d0
    isoname(274)= 'Yb'
    isomass(274)=171.936393d0
    isof(274)=21.9d0
    isoname(275)= 'Yb'
    isomass(275)=172.938222d0
    isof(275)=16.12d0
    isoname(276)= 'Yb'
    isomass(276)=173.938873d0
    isof(276)=31.8d0
    isoname(277)= 'Yb'
    isomass(277)=175.942576d0
    isof(277)=12.7d0
    isoname(278)= 'Zn'
    isomass(278)=63.929145d0
    isof(278)=48.6d0
    isoname(279)= 'Zn'
    isomass(279)=65.926035d0
    isof(279)=27.9d0
    isoname(280)= 'Zn'
    isomass(280)=66.927129d0
    isof(280)=4.1d0
    isoname(281)= 'Zn'
    isomass(281)=67.924846d0
    isof(281)=18.8d0
    isoname(282)= 'Zn'
    isomass(282)=69.925325d0
    isof(282)=0.6d0
    isoname(283)= 'Zr'
    isomass(283)=89.904708d0
    isof(283)=51.45d0
    isoname(284)= 'Zr'
    isomass(284)=90.905644d0
    isof(284)=11.27d0
    isoname(285)= 'Zr'
    isomass(285)=91.905039d0
    isof(285)=17.17d0
    isoname(286)= 'Zr'
    isomass(286)=93.906319d0
    isof(286)=17.33d0
    isoname(287)='Zr'
    isomass(287)=95.908272d0
    isof(287)=2.78d0


!write(*,*)'isoname(287)',isoname(287)



 numiso=size(isoname)
 mave=0
write(*,*) 'natom_type', natom_type
!do nn=1,natom_type
!write(*,*)'atname(nn)', atname(nn)
!write(*,*)'numiso', numiso
     do ii=1,numiso
! write(*,*)'isoname(ii)', isoname(ii),atname(nn)
          if (isoname(ii).eq.attname) then
           mave=mave+isomass(ii)*isof(ii)
write(*,*) ii, isoname(ii), attname, mave
        endif
     enddo
!enddo

g=0
mave=mave/100
!do nn=1,natom_type
   do ii=1,numiso
       if(isoname(ii).eq.attname) then
           g=g+isof(ii)*(1-(isomass(ii)/mave))**2
       endif
   enddo
!enddo
!g=0
!do nn=1,natom_type

!g=g+geatom(nn)

!enddo

g=g/100
write(*,*)'cal_gggg',g,mave

end subroutine cal_g
!==========================================================
subroutine cal_giso(gis)

  use ios
  use params
  use lattice
  use atoms_force_constants


implicit none

real(8) g(natom_type),gis(natoms)
integer i,j

write(*,*)'natom_type-natoms0',natom_type, natoms0

          do i=1,natom_type
           call cal_g(atname(i),g(i))
           write(*,*)'g',atname(i),g(i)
          enddo
          write(*,*)'safoura-after-calg'
          do i=1,natoms0
             do j=1,natom_type
               if (atom_type(i).eq.j) then
                 gis(i)=g(j)
                 write(*,*)'natom,natom_type',i,j,g(j),gis(i) 
               endif
            enddo
         enddo


end subroutine cal_giso
!===========================================================




!==========================================================
subroutine P_iso(g,ndn,evec)


use ios
use exactBTE2
use mod_ksubset2
use params
use constants
use phi3
use lattice
use tetrahedron
use kpoints
implicit none



integer n1,i,l1,n2,l2,ndn
real(8) omega1,omega2,g!,Piso(nkc,nkc,ndn,ndn)
real(8) q1(3),q2(3),nb1,nb2
complex(8)  evec(ndn,ndn,nkc)
integer j_tet, k_tet,w
integer nx, ny, nz
!integer s
nx=NC(1); ny=NC(2); nz=NC(3)


!g=?
Piso=0


do n1=1,nibz
   i=mapinv(n1)
  do l1=1,ndn

     nb1=dist(i,l1)
     omega1=frequency(i,l1)

     do n2=1,nkc
        do l2=1,ndn

           nb2=dist(n2,l2)
           omega2=frequency(n2,l2)


             ! eigen_tet: set omega2+omega3
             do j_tet=1,nkc*6                       !iterate over all tetrahedrons
                 do k_tet=1,4       !iterate over four corners of tetrahedron
                      w = tet(j_tet)%p(k_tet)%i !label for kpt(3,i)

                      q1=kpc(:,i)
                      q2=kpc(:,w)

                    ! q3=-q1-q2
                    ! q3_proc1=q1+q2

                    !call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                    !call
                    !get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)



                    tet(j_tet)%p(k_tet)%w=frequency(w,l2)! CHECK frequency(w,jj)!+frequency(nk3_proc1,kk) !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
                    ! s=frequency(w,l2)
                enddo
             enddo

             call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors in tetrahedron method

           do j_tet=1,nkc*6                       !iterate over all tetrahedrons
                 do k_tet=1,4   !iterate over four corners of tetrahedron
                     w = tet(j_tet)%p(k_tet)%i                ! k point
                     !nb2=dist(w,l2)
                     !omega2=frequency(w,l2)
           !          Piso(n1,w,l1,l2) =  nb1*(nb2+1)*omega1*omega2*g* tet(j_tet)%p(k_tet)%c/6d0 *(abs(dot_product(evec(:,l1,i),evec(:,l2,w))))**2!* F_arbi(w)
                     !Piso(n1,w,l1,l2) =  nb1*(nb2+1)*omega1*omega2*g* tet(j_tet)%p(k_tet)%c/6d0*(abs(sum(conjg(eigenvec(:,l1,n1))*eigenvec(:,l2,w))))**2

           !          write(*,*)n1,w,l1,l2,Piso(n1,w,l1,l2) 
                 enddo
           enddo


      enddo
    enddo
  enddo
enddo


endsubroutine P_iso
!=========================================================



!==========================================================
subroutine P_iso_split(ksubset,g,ndn,evec)


use ios
use exactBTE2
use mod_ksubset2
use params
use constants
use phi3
use lattice
use tetrahedron
use kpoints
use atoms_force_constants
implicit none



integer n1,i,l1,n2,l2,ndn
real(8) omega1,omega2,g(natoms0) !,g(natom_type)!,Piso(nkc,nkc,ndn,ndn)
real(8) q1(3),q2(3),nb1,nb2,gs
complex(8)  evec(ndn,ndn,nkc)
integer j_tet, k_tet,w
integer nx, ny, nz,ai
integer ksubset(2), indx, k_shift,ii,jj
integer ksub_size2
nx=NC(1); ny=NC(2); nz=NC(3)


!g=?
Pisos=0


if (ksubset(2) .eq. nibz) then   ! if this is the last file in v33sq.xxx.dat
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


!do n1=1,nibz
do n1=1,ksub_size2
   i=mapinv(n1+k_shift)
  do l1=1,ndn

     nb1=dist(i,l1)
     omega1=frequency(i,l1)

     do n2=1,nkc
        do l2=1,ndn

           nb2=dist(n2,l2)
           !omega2=frequency(n2,l2)

      !**      do ai=1,natoms0
             ! eigen_tet: set omega2+omega3
             do j_tet=1,nkc*6                       !iterate over all tetrahedrons
                 do k_tet=1,4       !iterate over four corners of tetrahedron
                      w = tet(j_tet)%p(k_tet)%i !label for kpt(3,i)

                      q1=kpc(:,i)
                      q2=kpc(:,w)

                    ! q3=-q1-q2
                    ! q3_proc1=q1+q2

                    !call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                    !call
                    !get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)



                    tet(j_tet)%p(k_tet)%w=frequency(w,l2)! CHECK!frequency(w,jj)!+frequency(nk3_proc1,kk) !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
                    ! s=frequency(w,l2)
                enddo
             enddo

             call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors in tetrahedron method

           do j_tet=1,nkc*6                       !iterate over all tetrahedrons
                 do k_tet=1,4   !iterate over four corners of tetrahedron
                     w = tet(j_tet)%p(k_tet)%i                ! k point
    !                 nb2=dist(w,l2)
    !                 omega2=frequency(w,l2)
                      gs=0
                     do ai=1,natoms0

                       gs = gs + g(ai) * (abs(dot_product(evec(((ai-1)*3+1):((ai-1)*3+3),l1,i),evec(((ai-1)*3+1):((ai-1)*3+3),l2,w))))**2
                     enddo


                    Pisos(n1,l1) = Pisos(n1,l1)  + gs * tet(j_tet)%p(k_tet)%c/6d0 !&
                          !*(abs(dot_product(evec(((ai-1)*3+1):((ai-1)*3+3),l1,i),evec(((ai-1)*3+1):((ai-1)*3+3),l2,w))))**2!* F_arbi(w)
                     !Piso(n1,w,l1,l2) =  nb1*(nb2+1)*omega1*omega2*g*
                     !tet(j_tet)%p(k_tet)%c/6d0*(abs(sum(conjg(eigenvec(:,l1,n1))*eigenvec(:,l2,w))))**2

           !          write(*,*)n1,w,l1,l2,Piso(n1,w,l1,l2)
           !      enddo
           enddo
         enddo
       !write(*,*)n1,i,n2,w,l1,l2, Pisos(n1,w,l1,l2)
!       Pisos(n1,l1) = Pisos(n1,w) * nb1*(nb2+1)*omega1**2* pi/(2*nkc)
      enddo
    enddo
       Pisos(n1,l1) = Pisos(n1,l1) *omega1**2* pi/(2*nkc)
  enddo
enddo


open(901,file='Piso.dat',status='unknown',form='unformatted')
do ii=1,ksub_size2
   do l1=1,ndn
!   do jj=1,nkc
!   do l2=1,ndn
      write(901) Pisos(ii,l1)
!   enddo
!   enddo
   enddo
enddo

!close(901)


endsubroutine P_iso_split
!=========================================================




!===========================================================

subroutine read_Piso(ndn,piss) !!Saf!!
! this subroutine reads the Cs_sycollisionm.dat
 use ios
 use phi3
 use kpoints
 implicit none
 integer ndn,n1,n2,l1,l2!kil,i,j,kil1,s,kil2
 character(99) filename_temp, filename
 !real(8), allocatable :: col_sy(:,:)
 integer indx, pis, ksubset(2), ksub_size2,k_shift !colsy
 real(8) xxs, piss(nibz,ndn)!, cs(nibz*ndn*3,nibz*ndn*3)

!kil=nibz*ndn*3
!kil1=ksub_size*ndn*3
pis=1004
indx=0
piss=0
!    if(allocated(col_sy) ) then
!        deallocate(col_sy)
!    endif

!    allocate(col_sy(kil,kil))
 do while (indx .lt. ((nibz/ksub_size)+mod(nibz,ksub_size)))

     indx=indx+1

!open(colsy,file='ksubset.inp', status='old') Cs_sycollisionm.k001.dat

    write(filename_temp,fmt="(a,i3.3,a)") "Piso.",indx,".dat"
    filename=trim(v3path)//filename_temp
    open(pis,file=filename,status='old',form='unformatted')
!    open(colsy,file=filename,status='old')

    if (indx .eq. ((nibz/ksub_size)+mod(nibz,ksub_size))) then
       ksub_size2=nibz-((indx-1)*ksub_size)
     else
      ksub_size2=ksub_size
    endif
!ksub_size2=ksubset(2)-ksubset(1)+1
k_shift = (indx-1)*ksub_size
!s=(kil1*(indx-1))+1
do n1=1,ksub_size2
   ! s=s+1
   do l1=1,ndn
    !do n2=1,nkc
    !do l2=1,ndn

     read(pis) piss((n1+k_shift),l1) !xxs !col_sy(s,j)

        !piss((n1+k_shift),n2,l1,l2) = xxs

    ! write(*,*) s, j, cs(s,j)
    !enddo
   !enddo
  enddo
enddo
! if (s.lt.kil1) then
!    s=s+1
! endif
! enddo
!write(*,*)'%%%%%%%%%%%%%%%%%%out'
! close(colsy)
enddo



!do n1=1,nibz
   ! s=s+1
!   do l1=1,ndn
!    do n2=1,nkc
!    do l2=1,ndn

!     Piso(n1,n2,l1,l2)= piss(n1,n2,l1,l2) !xxs !col_sy(s,j)

        !piss((n1+k_shift),n2,l1,l2) = xxs

    ! write(*,*) s, j, cs(s,j)
!    enddo
!   enddo
!  enddo
!enddo




end subroutine read_Piso
!==========================================================


!==========================================================
subroutine selfenergyw3(kbsubset,ndn,tmp,tmpk,nkpbss,kpbss,dkbss) !!Saf!!


use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants

implicit none

  integer kbsubset(2)
  integer nbss,ndn,la,lll,nkpbss
  complex(8) nselfss,uselfss
  real(8) , allocatable ::evlbbs(:,:),vlcbbs(:,:,:)!,giso(:)
  complex(8) , allocatable ::evcbbs(:,:,:)
  real(8) omegab0(ndn),kbbss(3),kpbss(3,nkpbss),dkbss(nkpbss)
  real(8) tmp,tmpk

lll = 5322
 !        call make_kp_bs


         allocate(evlbbs(ndn,nkpbss),evcbbs(ndn,ndn,nkpbss),vlcbbs(3,ndn,nkpbss))!,omega0(ndyn),tse(ndyn))
         call get_frequencies(nkpbss,kpbss,dkbss,ndn,evlbbs,ndn,evcbbs,lll,vlcbbs)


      do nbss=kbsubset(1),kbsubset(2)

              do la=1,ndn
                 kbbss=kp_bs(:,nbss)
                 omegab0(la)=sqrt(abs(evlbbs(la,nbss))) * cnst

                call function_self_w3(kbbss,la,omegab0(la),tmp,nselfss,uselfss)

                isempi(nbss,la)=2*aimag(nselfss+uselfss)
                 tsempi(nbss,la)=1d10/c_light/isempi(nbss,la)
!                 write(lwse,*)tempk,la,nkbss,kbss,dk_bs(nkbss),omega0(la),aimag(nselfs+uselfs),ise,tse(la),1/tse(la)
                 write(*,*)'nbss,la', nbss,la !nselfs, uselfs
               enddo
!               write(lwse2,*)tempk,dk_bs(nkbss),(omega0(la),la=1,ndyn),((1/tse(la)),la=1,ndyn)
              ! write(*,*)tempk,dk_bs(nkbss),(omega0(la),la=1,ndyn),((1/tse(la)),la=1,ndyn)
            enddo

 !       enddo

end subroutine selfenergyw3
!==========================================================



!===========================================================

subroutine read_lw(ndn,nks) !!Saf!!
! this subroutine reads the Cs_sycollisionm.dat
 use ios
 use phi3
 use kpoints
 implicit none
 integer ndn,n1,l1,nks!,l2!kil,i,j,kil1,s,kil2
 character(99) filename_temp, filename
 !real(8), allocatable :: col_sy(:,:)
 integer indx, lws, ksubset(2), ksub_size2!,k_shift !colsy
 !real(8) xxs, piss(nibz,nkc,ndn,ndn)!, cs(nibz*ndn*3,nibz*ndn*3)
real(8) ls1,ls2,ls3, ls4,ls5,ls6,ls7,ls8,ls9,ls10,ls11,ls12
!kil=nibz*ndn*3
!kil1=ksub_size*ndn*3
lws=1004
indx=0
!piss=0
!    if(allocated(col_sy) ) then
!        deallocate(col_sy)
!    endif
!open(1005,file=)
open(1005,file='lw-selfenergy.dat')
!    allocate(col_sy(kil,kil))
 do while (indx .lt. ((nks/ksub_size)+mod(nks,ksub_size)))

     indx=indx+1

!open(colsy,file='ksubset.inp', status='old') Cs_sycollisionm.k001.dat

    write(filename_temp,fmt="(a,i3.3,a)") "lws.",indx,".dat"
    filename=trim(v3path)//filename_temp
    open(lws,file=filename,status='old')!,form='unformatted')
!    open(colsy,file=filename,status='old')

    if (indx .eq. ((nks/ksub_size)+mod(nks,ksub_size))) then
       ksub_size2=nks-((indx-1)*ksub_size)
     else
      ksub_size2=ksub_size
    endif
!ksub_size2=ksubset(2)-ksubset(1)+1
!k_shift = (indx-1)*ksub_size
!s=(kil1*(indx-1))+1
do n1=1,ksub_size2
   ! s=s+1
   do l1=1,ndn
!    do n2=1,nkc
!    do l2=1,ndn

     read(lws,*) ls1,ls2,ls3, ls4,ls5,ls6,ls7,ls8,ls9,ls10,ls11,ls12
        !piss((n1+k_shift),n2,l1,l2) = xxs
     write(1005,*)ls1,ls2,ls3, ls4,ls5,ls6,ls7,ls8,ls9,ls10,ls11,ls12
    ! write(*,*) s, j, cs(s,j)
  !  enddo
 !  enddo
  enddo
enddo
! if (s.lt.kil1) then
!    s=s+1
! endif
! enddo
!write(*,*)'%%%%%%%%%%%%%%%%%%out'
! close(colsy)
enddo

close(1005)

!do n1=1,nibz
   ! s=s+1
!   do l1=1,ndn
!    do n2=1,nkc
!    do l2=1,ndn

!     Piso(n1,n2,l1,l2)= piss(n1,n2,l1,l2) !xxs !col_sy(s,j)

        !piss((n1+k_shift),n2,l1,l2) = xxs

    ! write(*,*) s, j, cs(s,j)
!    enddo
!   enddo
!  enddo
!enddo




end subroutine read_lw
!==========================================================



!===========================================================

subroutine read_lw2(ndn,nks) !!Saf!!
! this subroutine reads the Cs_sycollisionm.dat
 use ios
 use phi3
 use kpoints
 implicit none
 integer ndn,n1,l1,nks!,l2!kil,i,j,kil1,s,kil2
 character(99) filename_temp, filename
 !real(8), allocatable :: col_sy(:,:)
 integer indx, lws, ksubset(2), ksub_size2!,k_shift !colsy
 !real(8) xxs, piss(nibz,nkc,ndn,ndn)!, cs(nibz*ndn*3,nibz*ndn*3)
real(8) ls1,ls2,ls3, ls4,ls5,ls6,ls7,ls8,ls9,ls10,ls11,ls12
 real(8) lwssf2(2+2*ndn)
!kil=nibz*ndn*3
!kil1=ksub_size*ndn*3
lws=1004
indx=0
!piss=0
!    if(allocated(col_sy) ) then
!        deallocate(col_sy)
!    endif
!open(1005,file=)
open(1005,file='lw-selfenergy2.dat')
!    allocate(col_sy(kil,kil))
 do while (indx .lt. ((nks/ksub_size)+mod(nks,ksub_size)))

     indx=indx+1

!open(colsy,file='ksubset.inp', status='old') Cs_sycollisionm.k001.dat

    write(filename_temp,fmt="(a,i3.3,a)") "lws2.",indx,".dat"
    filename=trim(v3path)//filename_temp
    open(lws,file=filename,status='old')!,form='unformatted')
!    open(colsy,file=filename,status='old')

    if (indx .eq. ((nks/ksub_size)+mod(nks,ksub_size))) then
       ksub_size2=nks-((indx-1)*ksub_size)
     else
      ksub_size2=ksub_size
    endif
!ksub_size2=ksubset(2)-ksubset(1)+1
!k_shift = (indx-1)*ksub_size
!s=(kil1*(indx-1))+1
do n1=1,ksub_size2
   ! s=s+1
   !do l1=1,ndn
!    do n2=1,nkc
!    do l2=1,ndn

     read(lws,*) lwssf2 !ls1,ls2,ls3, ls4,ls5,ls6,ls7,ls8,ls9,ls10,ls11,ls12
        !piss((n1+k_shift),n2,l1,l2) = xxs
     write(1005,*) lwssf2 !ls1,ls2,ls3, ls4,ls5,ls6,ls7,ls8,ls9,ls10,ls11,ls12
    ! write(*,*) s, j, cs(s,j)
  !  enddo
 !  enddo
  !enddo
enddo
! if (s.lt.kil1) then
!    s=s+1
! endif
! enddo
!write(*,*)'%%%%%%%%%%%%%%%%%%out'
! close(colsy)
enddo

close(1005)

!do n1=1,nibz
   ! s=s+1
!   do l1=1,ndn
!    do n2=1,nkc
!    do l2=1,ndn

!     Piso(n1,n2,l1,l2)= piss(n1,n2,l1,l2) !xxs !col_sy(s,j)

        !piss((n1+k_shift),n2,l1,l2) = xxs

    ! write(*,*) s, j, cs(s,j)
!    enddo
!   enddo
!  enddo
!enddo




end subroutine read_lw2
!==========================================================





!===========================================================

subroutine read_Collision_Matrix(ndn,cs) !!Saf!!
! this subroutine reads the Cs_sycollisionm.dat
 use ios
 use phi3
 use kpoints
 implicit none
 integer ndn,kil,i,j,kil1,s,kil2
 character(99) filename_temp, filename
 !real(8), allocatable :: col_sy(:,:)
 integer indx, colsy
 real(8) xxs, cs(nibz*ndn*3,nibz*ndn*3)

kil=nibz*ndn*3
kil1=ksub_size*ndn*3
colsy=1004
indx=0

!    if(allocated(col_sy) ) then
!        deallocate(col_sy)
!    endif

!    allocate(col_sy(kil,kil))
 do while (indx .lt. ((nibz/ksub_size)+mod(nibz,ksub_size)))

     indx=indx+1

!open(colsy,file='ksubset.inp', status='old')
!Cs_sycollisionm.k001.dat

    write(filename_temp,fmt="(a,i3.3,a)") "Cs_sycollisionm.",indx,".dat"
    filename=trim(v3path)//filename_temp
    open(colsy,file=filename,status='old',form='unformatted')
!    open(colsy,file=filename,status='old')

    if (indx .eq. ((nibz/ksub_size)+mod(nibz,ksub_size))) then
       kil2=kil-((indx-1)*kil1)
     else
      kil2=kil1
    endif


s=(kil1*(indx-1))+1
 do i=1,kil2
   ! s=s+1
    do j=1,kil

      read(colsy) xxs !col_sy(s,j)
     
        cs(s,j) = xxs

     write(*,*) s, j, cs(s,j)
    enddo
! if (s.lt.kil1) then
    s=s+1
! endif
 enddo
!write(*,*)'%%%%%%%%%%%%%%%%%%out'
! close(colsy)
enddo

end subroutine read_Collision_Matrix
!==========================================================
!==========================================================
 subroutine read_RHS(nik,ndn,RHSi)!!Saf!!
! this subroutine reads RHS_IBZ.dat
 use ios
 implicit none
 integer rhs_f
 integer i,nik,ndn,kil
 !real(8), allocatable :: RHSi(:)
 real(8) xx, RHSi(nik*ndn*3)
!write(*,*)'here%%%%%%%%%%#########'

 kil=nik*ndn*3
! write(*,*)'here%%%%%%%%%%@@@@@@@@@'
 rhs_f=5999
! RHS2=0
!write(*,*)'here%%%%%%%%%%@@@@@@@@@'
  !  if(allocated(RHS2)) then
  !      deallocate(RHS2)
  !  endif

 !   allocate(RHSi(kil))
!write(*,*)'here%%%%%%%%%%'
RHSi=0
xx=0
 open(rhs_f,file='RHS_IBZ.dat', status='old')
 do i=1,kil
 read(rhs_f,*) xx !RHS2(i)
 RHSi(i)=xx
 !write(*,*) RHSi(i)
 enddo
 close(rhs_f)

end subroutine read_RHS
!==========================================================
!==========================================================
 subroutine read_CM_1D_mpi(nik,ndn,CM3D)!!Saf!!
! this subroutine reads RHS_IBZ.dat
 use ios
 implicit none
 integer CM_f
 integer i,j,s,nik,ndn,kil
 !real(8), allocatable :: RHSi(:)
 real(8) xx, CM3D(nik*ndn*3,nik*ndn*3)
!write(*,*)'here%%%%%%%%%%#########'

 kil=nik*ndn*3
! write(*,*)'here%%%%%%%%%%@@@@@@@@@'
 CM_f=5999
! RHS2=0
!write(*,*)'here%%%%%%%%%%@@@@@@@@@'
  !  if(allocated(RHS2)) then
  !      deallocate(RHS2)
  !  endif

 !   allocate(RHSi(kil))
!write(*,*)'here%%%%%%%%%%'
CM3D=0
xx=0
s=0
 open(CM_f,file='CM-1D.dat', status='old')
 do i=1,kil
    do j=1,kil
       s=s+1
       read(CM_f,*) xx !RHS2(i)
       CM3D(i,j)=xx
 !write(*,*) RHSi(i)
     enddo
 enddo
 close(CM_f)

end subroutine read_CM_1D_mpi
!==========================================================
!==========================================================
 subroutine read_DS(nik,ndn,DS1d)!!Saf!!
! this subroutine reads RHS_IBZ.dat
 use ios
 implicit none
 integer DS_f
 integer i,nik,ndn,kil
 !real(8), allocatable :: RHSi(:)
 real(8) xx, DS1d(nik*ndn*3)
!write(*,*)'here%%%%%%%%%%#########'

 kil=nik*ndn*3
! write(*,*)'here%%%%%%%%%%@@@@@@@@@'
 DS_f=5999
! RHS2=0
!write(*,*)'here%%%%%%%%%%@@@@@@@@@'
  !  if(allocated(RHS2)) then
  !      deallocate(RHS2)
  !  endif

 !   allocate(RHSi(kil))
!write(*,*)'here%%%%%%%%%%'
DS1d=0
xx=0
 open(DS_f,file='DS.dat', status='old')
 do i=1,kil
 read(DS_f,*) xx !RHS2(i)
 DS1d(i)=xx
 !write(*,*) RHSi(i)
 enddo
 close(DS_f)

end subroutine read_DS
!==========================================================







!==========================================================
subroutine veloc_avergeIBZ_invsyFBZ(nk,kp,nibbz,kibbz,ndn,veloc,veloci) !!Saf!!

use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
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











!=================================== subroutine above just for one q from ibz
subroutine kstarsinfbz1(q,nk,kp,ds,kvecstars,kvecops) !!Saf!!


use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
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
!use ios
!use mod_ksubset2
!use phi3
!use kpoints
!use atoms_force_constants
!implicit none

!do i=nk


!end subroutine test_RTA
!========================================================================================
subroutine Pmatrix3 (indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)

 use kpoints
 use ios
 use phi3
 use lattice
! use geometry
! use atoms_force_constants
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
 
! real(8), allocatable :: v33sq_5(:,:,:,:,:)
! real cputim1,cputim2,cputim0, cputim1_wall, cputim2_wall, cputim0_wall, wtime
 real cputim1, cputim2, cputim0
 real(8)  cputim1_wall, cputim2_wall, cputim0_wall, wtime

 ucol=9010
 unt=112

 P1=0.0; P2=0.0

call cpu_time(cputim1)
cputim1_wall = wtime()
call cpu_time(cputim0)
cputim0_wall = wtime()

 allocate(v33sq_5(ksub_size2,nk,ndn,ndn,ndn))
 v33sq_5=0

 write(filename_temp,fmt="(a,i3.3,a)") "v33sq.",indx,".dat"
 filename=trim(v3path)//filename_temp
 open(unt,file=filename,status='old',form='unformatted')
 read(unt) nv3_2

 do n1=1,ksub_size2
 do l1=1,ndn
 do n2=1,nk
 do l2=1,ndn
 do l3=1,ndn
     read(unt) v33sq_5(n1,n2,l1,l2,l3)
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

                 V3sq1=v33sq_5(n1,n2,l1,l2,l3) * (2*pi)**2 / nk
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

deallocate(v33sq_5)
close(unt)

7 format(5(i6),2(2x,g14.8))
8 format(i6,3(2x,g14.8))


end subroutine Pmatrix3



!==================================================================
!==============================================================
subroutine check_conv23(i,nk,kp,nibbz,kibbz,ndn,tempk,convergence,iter_cont) !!Saf!!

use constants
use exactBTE2
use params
use lattice
use ios
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

 use kpoints
 use ios
 use phi3
 use svd_stuff
 use eigen
 use params
 use constants
 use atoms_force_constants
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


 end subroutine matrix_elt_simplified


!==================================================================================
subroutine cal_eiqr2(nk,ndn)
! calculate eiqr2, eivec2 for efficient calculation of matrix elements.
use kpoints
use eigen
use pre_matrix_elt
use lattice
use params
 use atoms_force_constants
 use ios
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
 use ios
 use phi3
 use lattice
! use geometry
! use atoms_force_constants
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

real(8), allocatable :: tauinv_RTA(:,:)  !v33sq_5(:,:,:,:,:),
real(8) tauinv_RTA_N(nk,ndn), tauinv_RTA_U(nk,ndn)
real(8) kappa_k_RTA2(nk,ndn,3,3), kappa_RTA2(ndn,3,3), nbe, delta_l

real(8) tauinv_temp1, tauinv_temp2

unt=112
ucol=9010

allocate(v33sq_5(nk,nk,ndn,ndn,ndn))
v33sq_5=0
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
               write(ulog,*) 'ijk2=',i2,j2,k2
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

                 V3sq1=v33sq_5(n1,n2,l1,l2,l3) * (2*pi)**2 / nk

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

deallocate(v33sq_5)

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





!========================================================================================
subroutine Pmatrix4 (indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)
! calculate collision matrix P1 and P2


 use kpoints
 use ios
 use phi3
 use lattice
! use geometry
! use atoms_force_constants
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

! real(8), allocatable :: v33sq_5(:,:,:,:,:)
! real cputim1,cputim2,cputim0, cputim1_wall, cputim2_wall, cputim0_wall, wtime
 real cputim1, cputim2, cputim0
 real(8)  cputim1_wall, cputim2_wall, cputim0_wall, wtime

 ucol=9010
 unt=112


call cpu_time(cputim0)
cputim0_wall = wtime()

 allocate(v33sq_5(ksub_size2,nk,ndn,ndn,ndn))
 v33sq_5= 0d0

 write(filename_temp,fmt="(a,i3.3,a)") "v33sq.",indx,".dat"
 filename=trim(v3path)//filename_temp
 open(unt,file=filename,status='old',form='unformatted')
 read(unt) nv3_2

 do n1=1,ksub_size2
 do l1=1,ndn
 do n2=1,nk
 do l2=1,ndn
 do l3=1,ndn
     read(unt) v33sq_5(n1,n2,l1,l2,l3)
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
               write(ulog,*) 'ijk2=',i2,j2,k2
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

                     V3sq1=v33sq_5(n1-k_shift,n2,l1,l2,l3) * (2*pi)**2 / nk

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

deallocate(v33sq_5)
close(unt)

7 format(5(i6),2(2x,g14.8))
8 format(i6,3(2x,g14.8))
9 format(99(2x,g14.8))

end subroutine Pmatrix4


!---------------------------------------------------------------------------------------
subroutine calculate_v3sq_delta (ksubset,ndn,nk,kp)
!- read v33sq from v33sq from v33sq.xxx.dat
!- calculate delta function for coalescence and decay processes
!- write v33sq * delta to files

use ios
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
real(8), allocatable :: P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)  ! v33sq_5(:,:,:,:,:), 
! P1sm and P2sm are v33sq*delta(proc1) and v33sq*delta(proc2)

nx=NC(1); ny=NC(2); nz=NC(3)
unt=9400

call allocate_iter(nk,nk,ndn)
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

if (allocated(v33sq_5)) then
   deallocate(v33sq_5)
endif

allocate(v33sq_5(ksub_size2,nk,ndn,ndn,ndn))

do n1=1,ksub_size2
do l1=1,ndn
do n2=1,nk
do l2=1,ndn
do l3=1,ndn
       read(unt) v33sq_5(n1,n2,l1,l2,l3)
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


                V3sq1=v33sq_5(i-k_shift,j,ii,jj,kk) * (2*pi)**2 / nk
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


                V3sq1=v33sq_5(i-k_shift,j,ii,jj,kk) * (2*pi)**2 / nk
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


!! K1    call cpu_time(cputim2)
!!    write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1  ! cputim1 unknown!

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


deallocate(v33sq_5)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif



end subroutine calculate_v3sq_delta

!---------------------------------------------------------------------------------------
subroutine calculate_FGR(ksubset,ndn,nv3oi)!,v33sq_8) !  ,nk,kp)
!- read v33sq from v33sq from v33sq.xxx.dat
!- calculate delta function for coalescence and decay processes
!- write v33sq * delta to files

use ios
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
integer nk_subset, k_shift, ksub_size2
real(8) q1(3), q2(3), q2n(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
real(8) temp, tempk, nb1   ! in cm^-1
integer ksubset(2)
integer nx, ny, nz

integer i2,ii2,j2,jj2,k2,kk2, inside1, inside2, inside, inside1n, inside2n
integer i3,j3,k3,nk3_proc1,nk3_proc2, nk3, nk2
integer n1,n2,l1,l2,l3,nk2n,i2n,j2n,k2n
integer j_tet, k_tet,w

character(99) filename_temp, filename

integer uQ, nv3_2,unt,nv3oi,s

real(8) cputim0, cputim1, cputim2, cputim3, cputim0_wall, cputim1_wall, wtime
real(8) omega1, V3sq1!,v33sq_8((ksubset(2)-ksubset(1)+1)*nkc*ndn*ndn*ndn)
real(8) F_arbi(nkc)
real(8), allocatable :: P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)  ! v33sq_5(:,:,:,:,:), 
! P1sm and P2sm are v33sq*delta(proc1) and v33sq*delta(proc2)

nx=NC(1); ny=NC(2); nz=NC(3)
unt=9400

call allocate_iter(nibz,nkc,ndn)
call cal_freq(nkc,ndn)
write(ulog,*) 'cal_freq done'

write(ulog,*) 'entering calculate_FGR ...'

! calculate indx using nk and ksubset
! calculate ksub_size2
if (ksubset(2) .eq. nibz) then   ! if this is the last file in v33sq.xxx.dat
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
!***write(filename_temp,fmt="(a,i3.3,a)") "v33sq.",indx,".dat"
!***filename=trim(v3path)//filename_temp
!***write(ulog,*) 'opening file ',filename
!***open(unt,file=filename,status='old',form='unformatted')
!***read(unt) nv3_2
nv3_2=nv3oi
!****if (allocated(v33sq_5)) then
!****   deallocate(v33sq_5)
!****endif

!****allocate(v33sq_5(ksub_size2,nkc,ndn,ndn,ndn))
!s=0
!do n1=1,ksub_size2
!do l1=1,ndn
!do n2=1,nkc
!do l2=1,ndn
!do l3=1,ndn
!    s=s+1
!    v33sq_5(n1,n2,l1,l2,l3)=v33sq_8(s)
!enddo
!enddo
!enddo
!enddo
!enddo
!***close(unt)


! calculate delta function
if (allocated(P1sm)) then
    deallocate(P1sm)
endif
if (allocated(P2sm)) then
    deallocate(P2sm)
endif

allocate(P1sm(ksub_size2,nkc,ndn,ndn,ndn))
allocate(P2sm(ksub_size2,nkc,ndn,ndn,ndn))


do n1=1,ksub_size2

    i=mapinv(n1+k_shift)   ! n1 is the ibz index, i is the FBZ index

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
           do j=1,nkc  ! for set arbitrary function

                q1=kpc(:,i)
                q2=kpc(:,j)

!                q3=-q1-q2
                q3_proc1=q1+q2
!                q3_proc2=q1-q2

                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
!                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)


                V3sq1=v33s8(n1,j,ii,jj,kk) * (2*pi)**2 / nkc
                write(*,*)'v3sq1',v3sq1
                if (jj<=3 .and. i2.eq.NC(1)/2+1 .and. j2.eq.NC(2)/2+1 .and. k2.eq.NC(3)/2+1 ) then  ! Acoustic modes at Gamma point
                    V3sq1=0.0                                                                        ! to enforce matrix_elt = 0
                endif
                if (kk<=3 .and. i3.eq.NC(1)/2+1 .and. j3.eq.NC(2)/2+1 .and. k3.eq.NC(3)/2+1 ) then
                    V3sq1=0.0
                endif

                F_arbi(j)=V3sq1  

            enddo

            ! eigen_tet: set omega2+omega3
            do j_tet=1,nkc*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                    w = tet(j_tet)%p(k_tet)%i                          !label for kpt(3,i)

                    q1=kpc(:,i)
                    q2=kpc(:,w)

                    q3=-q1-q2
                    q3_proc1=q1+q2

                    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                    call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)

                    tet(j_tet)%p(k_tet)%w=-frequency(w,jj)+frequency(nk3_proc1,kk)     !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
               enddo
            enddo

           call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors in tetrahedron method

           do j_tet=1,nkc*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                   w = tet(j_tet)%p(k_tet)%i                ! k point
                   P1sm(n1,w,ii,jj,kk) = P1sm(n1,w,ii,jj,kk) + tet(j_tet)%p(k_tet)%c/6d0 * F_arbi(w) 
               enddo
           enddo


           ! process 2
           ! set arbitrary function here. (distribution function...)
           ! call eigen_tet using w2-w3
           ! call weight_tet

           do j=1,nkc  ! for set arbitrary function

                q1=kpc(:,i)
                q2=kpc(:,j)
                q2n=-kpc(:,j)

!                q3=-q1-q2n
!                q3_proc1=q1+q2
                q3_proc2=q1+q2

                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                call get_k_info(q2n,NC,nk2n,i2n,j2n,k2n,g1,g2,g3,inside)
!                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)


                V3sq1=v33s8(n1,j,ii,jj,kk) * (2*pi)**2 / nkc
                if (jj<=3 .and. i2n.eq.NC(1)/2+1 .and. j2n.eq.NC(2)/2+1 .and. k2n.eq.NC(3)/2+1 ) then  ! Acoustic modes at Gamma point
                    V3sq1=0.0                                                                        ! to enforce matrix_elt = 0
                endif
                if (kk<=3 .and. i3.eq.NC(1)/2+1 .and. j3.eq.NC(2)/2+1 .and. k3.eq.NC(3)/2+1 ) then
                    V3sq1=0.0
                endif

                F_arbi(nk2n)=V3sq1        ! no 2pi at the front since delta is in linear frequency.

            enddo

            ! eigen_tet: set omega2+omega3
            do j_tet=1,nkc*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                    w = tet(j_tet)%p(k_tet)%i                          !label for kpt(3,i)

                    q1=kpc(:,i)
                    q2=kpc(:,w)

!                    q3=-q1-q2
                    q3_proc2=q1-q2

!                    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                    call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside1)

                    tet(j_tet)%p(k_tet)%w=frequency(w,jj)+frequency(nk3_proc2,kk)     !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
               enddo
            enddo

           call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors

           do j_tet=1,nkc*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                   w = tet(j_tet)%p(k_tet)%i                ! k point
                   P2sm(n1,w,ii,jj,kk) = P2sm(n1,w,ii,jj,kk) + tet(j_tet)%p(k_tet)%c/6d0 * F_arbi(w)
               enddo
           enddo

           enddo laF3_loop
           enddo laF2_loop


    enddo laF1_loop


!! K1    call cpu_time(cputim2)
!!    write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1  ! cputim1 unknown!

enddo



! write P1sm and P2sm
open(unt+1,file='v33sq_delta.dat',status='unknown',form='unformatted')
write(unt+1) nv3_2
!write(*,*) 'nv3_2=',nv3_2
do n1=1,ksub_size2
do l1=1,ndn
do n2=1,nkc
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


!***deallocate(v33sq_5)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif



end subroutine calculate_FGR

!---------------------------------------------------------------------------------------
subroutine calculate_FGR_mpi2(ksubset,ndn,nv3oi)!,v33sq_8) !  ,nk,kp)
!- read v33sq from v33sq from v33sq.xxx.dat
!- calculate delta function for coalescence and decay processes
!- write v33sq * delta to files

use ios
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
integer nk_subset, k_shift, ksub_size2
real(8) q1(3), q2(3), q2n(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
real(8) temp, tempk, nb1   ! in cm^-1
integer ksubset(2)
integer nx, ny, nz

integer i2,ii2,j2,jj2,k2,kk2, inside1, inside2, inside, inside1n, inside2n
integer i3,j3,k3,nk3_proc1,nk3_proc2, nk3, nk2
integer n1,n2,l1,l2,l3,nk2n,i2n,j2n,k2n
integer j_tet, k_tet,w

character(99) filename_temp, filename

integer uQ, nv3_2,unt,nv3oi,s

real(8) cputim0, cputim1, cputim2, cputim3, cputim0_wall, cputim1_wall, wtime
real(8) omega1, V3sq1!,v33sq_8((ksubset(2)-ksubset(1)+1)*nkc*ndn*ndn*ndn)
real(8) F_arbi(nkc)
!*****real(8), allocatable :: P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)  ! v33sq_5(:,:,:,:,:), 
! P1sm and P2sm are v33sq*delta(proc1) and v33sq*delta(proc2)

nx=NC(1); ny=NC(2); nz=NC(3)
unt=9400

call allocate_iter(nibz,nkc,ndn)
call cal_freq(nkc,ndn)
write(ulog,*) 'cal_freq done'

write(ulog,*) 'entering calculate_FGR ...'

! calculate indx using nk and ksubset
! calculate ksub_size2
!***if (ksubset(2) .eq. nibz) then   ! if this is the last file in v33sq.xxx.dat
!***   ksub_size2=ksubset(2)-ksubset(1)+1
!***   indx=(ksubset(1)-1)/ksub_size + 1
!***   if (mod(ksubset(1)-1,ksub_size) .ne. 0 ) then
!***       write(ulog,*) 'wrong indx. stop',ksubset(1),ksub_size
!***       stop
!***   endif

!***else
!***    ksub_size2=ksubset(2)-ksubset(1)+1
!***    indx=ksubset(2)/ksub_size           ! ksub_size is from params.phon
!***    if (mod(ksubset(2),ksub_size) .ne. 0) then
!***       write(ulog,*) 'wrong indx. stop',ksubset(2),ksub_size
!***       stop
!***    endif

!***endif

 ksub_size2=ksubset(2)-ksubset(1)+1

k_shift=ksubset(1)-1

!***write(ulog,*) 'Check indx is same as folder name, indx=',indx



! Read v33sq from v33sq.xxx.dat in the main directory
!***write(filename_temp,fmt="(a,i3.3,a)") "v33sq.",indx,".dat"
!***filename=trim(v3path)//filename_temp
!***write(ulog,*) 'opening file ',filename
!***open(unt,file=filename,status='old',form='unformatted')
!***read(unt) nv3_2
nv3_2=nv3oi
!****if (allocated(v33sq_5)) then
!****   deallocate(v33sq_5)
!****endif

!****allocate(v33sq_5(ksub_size2,nkc,ndn,ndn,ndn))
!s=0
!do n1=1,ksub_size2
!do l1=1,ndn
!do n2=1,nkc
!do l2=1,ndn
!do l3=1,ndn
!    s=s+1
!    v33sq_5(n1,n2,l1,l2,l3)=v33sq_8(s)
!enddo
!enddo
!enddo
!enddo
!enddo
!***close(unt)


! calculate delta function
!*****if (allocated(P1sm)) then
!*****    deallocate(P1sm)
!*****endif
!*****if (allocated(P2sm)) then
!*****    deallocate(P2sm)
!*****endif

!*****allocate(P1sm(ksub_size2,nkc,ndn,ndn,ndn))
!*****allocate(P2sm(ksub_size2,nkc,ndn,ndn,ndn))


do n1=1,ksub_size2

    i=mapinv(n1+k_shift)   ! n1 is the ibz index, i is the FBZ index

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
           do j=1,nkc  ! for set arbitrary function

                q1=kpc(:,i)
                q2=kpc(:,j)

!                q3=-q1-q2
                q3_proc1=q1+q2
!                q3_proc2=q1-q2

                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
!                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)


                V3sq1=v33s8(n1,j,ii,jj,kk) * (2*pi)**2 / nkc
           !     write(*,*)'v3sq1',v3sq1
                if (jj<=3 .and. i2.eq.NC(1)/2+1 .and. j2.eq.NC(2)/2+1 .and. k2.eq.NC(3)/2+1 ) then  ! Acoustic modes at Gamma point
                    V3sq1=0.0                                                                        ! to enforce matrix_elt = 0
                endif
                if (kk<=3 .and. i3.eq.NC(1)/2+1 .and. j3.eq.NC(2)/2+1 .and. k3.eq.NC(3)/2+1 ) then
                    V3sq1=0.0
                endif

                F_arbi(j)=V3sq1  

            enddo

            ! eigen_tet: set omega2+omega3
            do j_tet=1,nkc*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                    w = tet(j_tet)%p(k_tet)%i                          !label for kpt(3,i)

                    q1=kpc(:,i)
                    q2=kpc(:,w)

                    q3=-q1-q2
                    q3_proc1=q1+q2

                    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                    call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)

                    tet(j_tet)%p(k_tet)%w=-frequency(w,jj)+frequency(nk3_proc1,kk)     !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
               enddo
            enddo

           call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors in tetrahedron method

           do j_tet=1,nkc*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                   w = tet(j_tet)%p(k_tet)%i                ! k point
                   P1smpi(n1,w,ii,jj,kk) = P1smpi(n1,w,ii,jj,kk) + tet(j_tet)%p(k_tet)%c/6d0 * F_arbi(w) 
               enddo
           enddo


           ! process 2
           ! set arbitrary function here. (distribution function...)
           ! call eigen_tet using w2-w3
           ! call weight_tet

           do j=1,nkc  ! for set arbitrary function

                q1=kpc(:,i)
                q2=kpc(:,j)
                q2n=-kpc(:,j)

!                q3=-q1-q2n
!                q3_proc1=q1+q2
                q3_proc2=q1+q2

                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                call get_k_info(q2n,NC,nk2n,i2n,j2n,k2n,g1,g2,g3,inside)
!                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)


                V3sq1=v33s8(n1,j,ii,jj,kk) * (2*pi)**2 / nkc
                if (jj<=3 .and. i2n.eq.NC(1)/2+1 .and. j2n.eq.NC(2)/2+1 .and. k2n.eq.NC(3)/2+1 ) then  ! Acoustic modes at Gamma point
                    V3sq1=0.0                                                                        ! to enforce matrix_elt = 0
                endif
                if (kk<=3 .and. i3.eq.NC(1)/2+1 .and. j3.eq.NC(2)/2+1 .and. k3.eq.NC(3)/2+1 ) then
                    V3sq1=0.0
                endif

                F_arbi(nk2n)=V3sq1        ! no 2pi at the front since delta is in linear frequency.

            enddo

            ! eigen_tet: set omega2+omega3
            do j_tet=1,nkc*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                    w = tet(j_tet)%p(k_tet)%i                          !label for kpt(3,i)

                    q1=kpc(:,i)
                    q2=kpc(:,w)

!                    q3=-q1-q2
                    q3_proc2=q1-q2

!                    call get_k_info(q3,NC,nk3,i3,j3,k3,g1,g2,g3,inside)
                    call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside1)

                    tet(j_tet)%p(k_tet)%w=frequency(w,jj)+frequency(nk3_proc2,kk)     !assigns eigenvalue for each band. eival(l,w) l=band#,w=kpt#
               enddo
            enddo

           call weight_tet(nx,ny,nz,omega1)   ! calculate weighting factors

           do j_tet=1,nkc*6                       !iterate over all tetrahedrons
               do k_tet=1,4                                         !iterate over four corners of tetrahedron
                   w = tet(j_tet)%p(k_tet)%i                ! k point
                   P2smpi(n1,w,ii,jj,kk) = P2smpi(n1,w,ii,jj,kk) + tet(j_tet)%p(k_tet)%c/6d0 * F_arbi(w)
               enddo
           enddo

           enddo laF3_loop
           enddo laF2_loop


    enddo laF1_loop


!! K1    call cpu_time(cputim2)
!!    write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1  ! cputim1 unknown!

enddo



! write P1sm and P2sm
!*****open(unt+1,file='v33sq_delta.dat',status='unknown',form='unformatted')
!*****write(unt+1) nv3_2
!write(*,*) 'nv3_2=',nv3_2
!*****do n1=1,ksub_size2
!*****do l1=1,ndn
!*****do n2=1,nkc
!*****do l2=1,ndn
!*****do l3=1,ndn
!*****   write(unt+1) P1sm(n1,n2,l1,l2,l3), P2sm(n1,n2,l1,l2,l3)
!write(*,*) 'n1,n2,l1,l2,l3',n1,n2,l1,l2,l3
!write(*,*) 'P1sm,P2sm',P1sm(n1,n2,l1,l2,l3),P2sm(n1,n2,l1,l2,l3)

!*****enddo
!*****enddo
!*****enddo
!*****enddo
!write(*,*) 'n1,n2,l1,l2,l3',n1,n2,l1,l2,l3
!write(*,*) 'P1sm,P2sm',P1sm(n1,n2,l1,l2,l3),P2sm(n1,n2,l1,l2,l3)
!*****enddo

!*****close(unt+1)


!***deallocate(v33sq_5)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif



end subroutine calculate_FGR_mpi2



!============================================================================================
subroutine cal_Qm3_tet2(nk,ndn,kp,tempk)
! calculate sum of diagonal terms in collision matrix
! same as cal_Qm2 but calculate Pmatrix and save it to memory


use ios
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
real(8) kp(3,nk), q1(3), q2(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
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
real(8), allocatable :: P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)  

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
write(*,*) 'entering cal_Qm3_tet2...'
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
    write(ulog,*) indx,' read v3sq_delta.dat done, TIME=',(cputim1-cputim0),'WALL TIME=',(cputim1_wall-cputim0_wall)

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



!============================================================================================
subroutine calculate_tauRTA(ndn,tempk)
!! calculates sum of diagonal terms in collision matrix
!! calculates Pmatrix and save it to memory
!! BEWARE ******
!! the first index ibz is in IBZ, if you want to get the corresponding FBZ index 
!! use i=mapinv(ibz)
!!
use ios
use exactBTE2
use mod_ksubset2
use params
use constants
use phi3
use lattice , only : g1,g2,g3
use tetrahedron
use kpoints

implicit none

integer i,ii,indx,ndn,j,jj,kk
integer nk_subset, k_shift, ksub_size2
real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
real(8) temp, tempk, nb1   ! in cm^-1
integer ksubset(2)
integer nx, ny, nz
integer i2,ii2,j2,jj2,k2,kk2, inside1, inside2, inside, inside1n, inside2n
integer i3,j3,k3,nk3_proc1,nk3_proc2, nk3, nk2,nk1,ibz,i1,j1,k1
integer n1,n2,l1,l2,l3,nk2n,i2n,j2n,k2n
integer j_tet, k_tet,w

character(99) filename_temp, filename

integer uQ, nv3_2,unt

real(8) cputim0, cputim1, cputim2, cputim3, cputim0_wall, cputim1_wall, wtime
real(8) omega1, V3sq1
real(8) F_arbi(nkc)
real(8), allocatable :: P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)  

nx=NC(1); ny=NC(2); nz=NC(3)

uQ=9402
unt=9400
open(uQ,file='Qvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Q'

Qvalue(:,:)=0    ! in cm^-1
Qvalue_N(:,:)=0
Qvalue_U(:,:)=0
tauinv_N(:,:)=0
tauinv_U(:,:)=0
indx=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

!write(*,*) 'entering cal_Qm2...'

call cpu_time(cputim1)
write(*,*) 'entering calculate_tauRTA...'

nkI1_loop: do ibz=1,nibz
                i=mapinv(ibz)
                call get_k_info(kibz(:,ibz) ,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
                if (nk1.ne.i) then
                   write(*,*)'calculate_tauRTA error: i.ne.nk1 ',i,nk1
                   stop
                endif

 ! call Pmatrix3(i)
     ! -> read V33.xxx.dat
     ! -> allocate P1,P2 subset
     ! -> calculate P1, P2
     ! -> close V33.xxx.dat
 if (mod(ibz,ksub_size) .eq. 1) then   ! start point of each v33.xxx.dat

    if(allocated(P1sm) .or. allocated(P2sm)) then
        deallocate(P1sm)
        deallocate(P2sm)
    endif

    indx=indx+1
    if (indx*ksub_size > nibz) then
        ksub_size2=nibz - (indx-1)*ksub_size
        ksubset(1)=ksub_size*(indx-1)+1    ! starting point of k subset
        ksubset(2)=nibz
    else
        ksub_size2=ksub_size
        ksubset(1)=ksub_size*(indx-1)+1
        ksubset(2)=ksub_size*indx
    endif
    k_shift=ksub_size*(indx-1)    ! P1 of n1,n2,l1,l2,l3 will be saved in P1(n1-k_shift,n2,l1,l2,l3).

!    call Pmatrix4(indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)
    call cpu_time(cputim0)
    cputim0_wall = wtime()

    allocate(P1sm(ksub_size2,nkc,ndn,ndn,ndn))
    allocate(P2sm(ksub_size2,nkc,ndn,ndn,ndn))

    P1sm=0d0
    P2sm=0d0

! read v33sq_delta.xxx.dat

    write(filename_temp,fmt="(a,i3.3,a)") "v33sq_delta.",indx,".dat"
    filename=trim(v3path)//filename_temp
    open(unt,file=filename,status='old',form='unformatted')
    read(unt) nv3_2


    do n1=1,ksub_size2
    do l1=1,ndn
    do n2=1,nkc
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
    write(ulog,*) indx,' read v3sq_delta.dat done, TIME=',(cputim1-cputim0),'WALL TIME=',(cputim1_wall-cputim0_wall)

  endif

! Here call calc_tet
! Calculate Q_N, Q_U, Q_tot

    laI1_loop: do ii=1,ndn

    write(ulog,*) 'i,ii=',i,ii

    omega1=frequency(i,ii)        ! set fixed omega1 in delta(omega1-omega2+-omega3)

           laF2_loop: do jj=1,ndn
           laF3_loop: do kk=1,ndn


           ! process 1
           ! set arbitrary function here. (distribution function...)
           ! call eigen_tet using w2-w3
           ! call weight_tet

           ! set arbitrary function, F_arbi(j)
           do j=1,nkc   ! for set arbitrary function
                q1=kpc(:,i)
                q2=kpc(:,j)

                q3=-q1-q2
                q3_proc1=q1+q2
                q3_proc2=q1-q2

                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                if (nk2.ne.j) then
                   write(*,*)'calculate_tauRTA error: j.ne.nk2 ',j,nk2
                   stop
                endif
                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)


                F_arbi(j)=dist(i,ii)*dist(j,jj)*(dist(nk3_proc1,kk)+1)   ! no 2pi at the front since delta is in linear frequency.
                P1(ibz,j,ii,jj,kk) = P1sm(ibz-k_shift,j,ii,jj,kk) * F_arbi(j)    

                F_arbi(j)=dist(i,ii)*(dist(j,jj)+1)*(dist(nk3_proc2,kk)+1)   
                P2(ibz,j,ii,jj,kk) = P2sm(ibz-k_shift,j,ii,jj,kk) * F_arbi(j)     

! K1   Qvalue(ibz,ii)=Qvalue(ibz,ii) + P1(ibz,j,ii,jj,kk) + 0.5*P2(ibz,j,ii,jj,kk)
                Qvalue(ibz,ii)=Qvalue(ibz,ii) + ( P1(ibz,j,ii,jj,kk)  &
                &              + 0.5*P2(ibz,j,ii,jj,kk) ) 

                if (inside1 .eq. 1) then    ! normal process in coalescence process
                    Qvalue_N(ibz,ii)=Qvalue_N(ibz,ii) + P1(ibz,j,ii,jj,kk)
                else
                    Qvalue_U(ibz,ii)=Qvalue_U(ibz,ii) + P1(ibz,j,ii,jj,kk)
                endif

                if (inside2 .eq. 1) then    ! normal process in decay process
                    Qvalue_N(ibz,ii)=Qvalue_N(ibz,ii) + 0.5*P2(ibz,j,ii,jj,kk)
                else
                    Qvalue_U(ibz,ii)=Qvalue_U(ibz,ii) + 0.5*P2(ibz,j,ii,jj,kk)
                endif

           enddo


           enddo laF3_loop
           enddo laF2_loop

           nb1=dist(i,ii)

           write(uQ,*) ibz,ii,Qvalue(ibz,ii),Qvalue(ibz,ii)/(nb1*(nb1+1))

! tauinv_N in 1/sec
           tauinv_N(ibz,ii) = Qvalue_N(ibz,ii)/(nb1*(nb1+1)) * c_light*100    
           tauinv_U(ibz,ii) = Qvalue_U(ibz,ii)/(nb1*(nb1+1)) * c_light*100
           tauinv_tot(ibz,ii) = tauinv_N(ibz,ii) + tauinv_U(ibz,ii)



    enddo laI1_loop



    call cpu_time(cputim2)
    write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1


enddo nkI1_loop



close(uQ)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif


end subroutine calculate_tauRTA
!=========================================================


!============================================================================================
subroutine calculate_tauRTA_mpi2(ibz_subset,ndn,tempk)!,rankmpi)!,P1sm,P2sm)
!! calculates sum of diagonal terms in collision matrix
!! calculates Pmatrix and save it to memory
!! BEWARE ******
!! the first index ibz is in IBZ, if you want to get the corresponding FBZ index 
!! use i=mapinv(ibz)
!!
use ios
use exactBTE2
use mod_ksubset2
use params
use constants
use phi3
use lattice , only : g1,g2,g3
use tetrahedron
use kpoints

implicit none

integer i,ii,indx,ndn,j,jj,kk
integer nk_subset, k_shift, ksub_size2
real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
real(8) temp, tempk, nb1   ! in cm^-1
integer ksubset(2)
integer nx, ny, nz
integer i2,ii2,j2,jj2,k2,kk2, inside1, inside2, inside, inside1n, inside2n
integer i3,j3,k3,nk3_proc1,nk3_proc2, nk3, nk2,nk1,ibz,i1,j1,k1
integer n1,n2,l1,l2,l3,nk2n,i2n,j2n,k2n
integer j_tet, k_tet,w

character(99) filename_temp, filename

integer uQ, nv3_2,unt

real(8) cputim0, cputim1, cputim2, cputim3, cputim0_wall, cputim1_wall, wtime
real(8) omega1, V3sq1
real(8) F_arbi(nkc)
!***real(8), allocatable :: P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)  
integer ibz_subset(2), rankmpi, ksubsize3


nx=NC(1); ny=NC(2); nz=NC(3)

uQ=9402
unt=9400
open(uQ,file='Qvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Q'

Qvaluempi(:,:)=0    ! in cm^-1
Qvalue_Nmpi(:,:)=0
Qvalue_Umpi(:,:)=0
tauinv_N(:,:)=0
tauinv_U(:,:)=0
indx=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

!write(*,*) 'entering cal_Qm2...'

call cpu_time(cputim1)
write(*,*) 'entering calculate_tauRTA...'

ksubsize3=(ibz_subset(2)-ibz_subset(1))+1
k_shift=ibz_subset(1)-1

nkI1_loop: do ibz=1,ksubsize3
!nkI1_loop: do ibz=1,nibz
                i=mapinv(ibz+k_shift)
                call get_k_info(kibz(:,ibz+k_shift) ,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
                if (nk1.ne.i) then
                   write(*,*)'calculate_tauRTA error: i.ne.nk1 ',i,nk1
                   stop
                endif

 ! call Pmatrix3(i)
     ! -> read V33.xxx.dat
     ! -> allocate P1,P2 subset
     ! -> calculate P1, P2
     ! -> close V33.xxx.dat
!*** if (mod(ibz,ksub_size) .eq. 1) then   ! start point of each v33.xxx.dat

!***    if(allocated(P1sm) .or. allocated(P2sm)) then
!***        deallocate(P1sm)
!***        deallocate(P2sm)
!***    endif

!***    indx=indx+1
!***    if (indx*ksub_size > nibz) then
!***        ksub_size2=nibz - (indx-1)*ksub_size
!***        ksubset(1)=ksub_size*(indx-1)+1    ! starting point of k subset
!***        ksubset(2)=nibz
!***    else
!***        ksub_size2=ksub_size
!***        ksubset(1)=ksub_size*(indx-1)+1
!***        ksubset(2)=ksub_size*indx
!***    endif
!***    k_shift=ksub_size*(indx-1)    ! P1 of n1,n2,l1,l2,l3 will be saved in P1(n1-k_shift,n2,l1,l2,l3).

!    call Pmatrix4(indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)
!***    call cpu_time(cputim0)
!***    cputim0_wall = wtime()

!***    allocate(P1sm(ksub_size2,nkc,ndn,ndn,ndn))
!***    allocate(P2sm(ksub_size2,nkc,ndn,ndn,ndn))

!***    P1sm=0d0
!***    P2sm=0d0

! read v33sq_delta.xxx.dat

!***    write(filename_temp,fmt="(a,i3.3,a)") "v33sq_delta.",indx,".dat"
!***    filename=trim(v3path)//filename_temp
!***    open(unt,file=filename,status='old',form='unformatted')
!***    read(unt) nv3_2


!***    do n1=1,ksub_size2
!***    do l1=1,ndn
!***    do n2=1,nkc
!***    do l2=1,ndn
!***    do l3=1,ndn
!***       read(unt) P1sm(n1,n2,l1,l2,l3), P2sm(n1,n2,l1,l2,l3)
!***    enddo
!***    enddo
!***    enddo
!***    enddo
!***    enddo

 !***   call cpu_time(cputim1)
 !***   cputim1_wall = wtime()
 !***   write(ulog,*) indx,' read v3sq_delta.dat done, TIME=',(cputim1-cputim0),'WALL TIME=',(cputim1_wall-cputim0_wall)

 !*** endif

! Here call calc_tet
! Calculate Q_N, Q_U, Q_tot
!ksubsize3=(ibz_subset(2)-ibz_subset(1))+1
!k_shift=ksubsize3*rankmpi

    laI1_loop: do ii=1,ndn

    write(ulog,*) 'i,ii=',i,ii

    omega1=frequency(i,ii)        ! set fixed omega1 in delta(omega1-omega2+-omega3)

           laF2_loop: do jj=1,ndn
           laF3_loop: do kk=1,ndn


           ! process 1
           ! set arbitrary function here. (distribution function...)
           ! call eigen_tet using w2-w3
           ! call weight_tet

           ! set arbitrary function, F_arbi(j)
           do j=1,nkc   ! for set arbitrary function
                q1=kpc(:,i)
                q2=kpc(:,j)

                q3=-q1-q2
                q3_proc1=q1+q2
                q3_proc2=q1-q2

                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                if (nk2.ne.j) then
                   write(*,*)'calculate_tauRTA error: j.ne.nk2 ',j,nk2
                   stop
                endif
                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)


                F_arbi(j)=dist(i,ii)*dist(j,jj)*(dist(nk3_proc1,kk)+1)   ! no 2pi at the front since delta is in linear frequency.
           !***     P1(ibz,j,ii,jj,kk) = P1sm(ibz-k_shift,j,ii,jj,kk) * F_arbi(j)    
                    P1(ibz,j,ii,jj,kk) = P1smpi(ibz,j,ii,jj,kk) * F_arbi(j)
                   ! P1(ibz,j,ii,jj,kk) = P1smpi(ibz-k_shift,j,ii,jj,kk) * F_arbi(j)
                   
                F_arbi(j)=dist(i,ii)*(dist(j,jj)+1)*(dist(nk3_proc2,kk)+1)   
           !***     P2(ibz,j,ii,jj,kk) = P2sm(ibz-k_shift,j,ii,jj,kk) * F_arbi(j)     
                    P2(ibz,j,ii,jj,kk) = P2smpi(ibz,j,ii,jj,kk) * F_arbi(j)
                 !   P2(ibz,j,ii,jj,kk) = P2smpi(ibz-k_shift,j,ii,jj,kk) * F_arbi(j)

! K1   Qvalue(ibz,ii)=Qvalue(ibz,ii) + P1(ibz,j,ii,jj,kk) + 0.5*P2(ibz,j,ii,jj,kk)
                Qvaluempi(ibz,ii)=Qvaluempi(ibz,ii) + ( P1(ibz,j,ii,jj,kk)  &
                &              + 0.5*P2(ibz,j,ii,jj,kk) ) 


!write(*,*)Qvaluempi(ibz,ii),ibz,ii

                if (inside1 .eq. 1) then    ! normal process in coalescence process
                    Qvalue_Nmpi(ibz,ii)=Qvalue_Nmpi(ibz,ii) + P1(ibz,j,ii,jj,kk)
                else
                    Qvalue_Umpi(ibz,ii)=Qvalue_Umpi(ibz,ii) + P1(ibz,j,ii,jj,kk)
                endif

                if (inside2 .eq. 1) then    ! normal process in decay process
                    Qvalue_Nmpi(ibz,ii)=Qvalue_Nmpi(ibz,ii) + 0.5*P2(ibz,j,ii,jj,kk)
                else
                    Qvalue_Umpi(ibz,ii)=Qvalue_Umpi(ibz,ii) + 0.5*P2(ibz,j,ii,jj,kk)
                endif

           enddo


           enddo laF3_loop
           enddo laF2_loop

           nb1=dist(i,ii)

          ! write(uQ,*) ibz,ii,Qvalue(ibz,ii),Qvalue(ibz,ii)/(nb1*(nb1+1))

! tauinv_N in 1/sec
          ! tauinv_N(ibz,ii) = Qvalue_N(ibz,ii)/(nb1*(nb1+1)) * c_light*100    
          ! tauinv_U(ibz,ii) = Qvalue_U(ibz,ii)/(nb1*(nb1+1)) * c_light*100
          ! tauinv_tot(ibz,ii) = tauinv_N(ibz,ii) + tauinv_U(ibz,ii)


write(uQ,*)Qvaluempi(ibz,ii),ibz,ii,ibz+k_shift,ii
!write(*,*)Qvaluempi(ibz,ii),ibz,ii,ibz+k_shift,ii
    enddo laI1_loop



   ! call cpu_time(cputim2)
   ! write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1


enddo nkI1_loop



close(uQ)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif


end subroutine calculate_tauRTA_mpi2
!=========================================================

!============================================================================================
subroutine Qandtau(ndn,tempk)!,P1sm,P2sm)
!! calculates sum of diagonal terms in collision matrix
!! calculates Pmatrix and save it to memory
!! BEWARE ******
!! the first index ibz is in IBZ, if you want to get the corresponding FBZ index
!! use i=mapinv(ibz)
!!
use ios
use exactBTE2
use mod_ksubset2
use params
use constants
use phi3
use lattice , only : g1,g2,g3
use tetrahedron
use kpoints

implicit none

integer i,ii,indx,ndn!,j,jj,kk
!integer nk_subset, k_shift, ksub_size2
!real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
real(8) temp, tempk, nb1   ! in cm^-1
!integer ksubset(2)
!integer nx, ny, nz
!integer i2,ii2,j2,jj2,k2,kk2, inside1, inside2, inside, inside1n, inside2n
!integer i3,j3,k3,nk3_proc1,nk3_proc2, nk3, nk2,nk1,ibz!,i1,j1,k1
 integer nk1,ibz,i1,j1,k1,inside
!integer n1,n2,l1,l2,l3,nk2n,i2n,j2n,k2n
!integer j_tet, k_tet,w

!character(99) filename_temp, filename

integer uQ, nv3_2,unt

real(8) cputim0, cputim1, cputim2, cputim3, cputim0_wall, cputim1_wall, wtime
real(8) omega1, V3sq1
!real(8) F_arbi(nkc)
!***real(8), allocatable :: P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)

!nx=NC(1); ny=NC(2); nz=NC(3)

uQ=9402
unt=9400
open(uQ,file='Qvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Q'


tauinv_N(:,:)=0
tauinv_U(:,:)=0


temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

!write(*,*) 'entering cal_Qm2...'

call cpu_time(cputim1)
write(*,*) 'entering calculate_tauRTA...'

nkI1_loop: do ibz=1,nibz
                i=mapinv(ibz)
                call get_k_info(kibz(:,ibz) ,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
                if (nk1.ne.i) then
                   write(*,*)'calculate_tauRTA error: i.ne.nk1 ',i,nk1
                   stop
                endif

    laI1_loop: do ii=1,ndn

    write(ulog,*) 'i,ii=',i,ii

    omega1=frequency(i,ii)        ! set fixed omega1 in delta(omega1-omega2+-omega3)
     
            
           nb1=dist(i,ii)

           write(uQ,*) ibz,ii,Qvalue(ibz,ii),Qvalue(ibz,ii)/(nb1*(nb1+1))

! tauinv_N in 1/sec
           tauinv_N(ibz,ii) = Qvalue_N(ibz,ii)/(nb1*(nb1+1)) * c_light*100
           tauinv_U(ibz,ii) = Qvalue_U(ibz,ii)/(nb1*(nb1+1)) * c_light*100
           tauinv_tot(ibz,ii) = tauinv_N(ibz,ii) + tauinv_U(ibz,ii)



    enddo laI1_loop



   ! call cpu_time(cputim2)
   ! write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1


enddo nkI1_loop



close(uQ)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif


end subroutine Qandtau
!=========================================================



!============================================================================================
subroutine calculate_tauRTA_iso(ndn,tempk)
!! calculates sum of diagonal terms in collision matrix
!! calculates Pmatrix and save it to memory
!! BEWARE ******
!! the first index ibz is in IBZ, if you want to get the corresponding FBZ index
!! use i=mapinv(ibz)
!!
use ios
use exactBTE2
use mod_ksubset2
use params
use constants
use phi3
use lattice , only : g1,g2,g3
use tetrahedron
use kpoints

implicit none

integer i,ii,indx,ndn,j,jj,kk
integer nk_subset, k_shift, ksub_size2
real(8) q1(3), q2(3), q3_proc1(3), q3_proc2(3), q3(3), q3_proc1n(3), q3_proc2n(3)
real(8) temp, tempk, nb1   ! in cm^-1
integer ksubset(2)
integer nx, ny, nz
integer i2,ii2,j2,jj2,k2,kk2, inside1, inside2, inside, inside1n, inside2n
integer i3,j3,k3,nk3_proc1,nk3_proc2, nk3, nk2,nk1,ibz,i1,j1,k1
integer n1,n2,l1,l2,l3,nk2n,i2n,j2n,k2n
integer j_tet, k_tet,w

character(99) filename_temp, filename

integer uQ, nv3_2,unt

real(8) cputim0, cputim1, cputim2, cputim3, cputim0_wall, cputim1_wall, wtime
real(8) omega1, V3sq1
real(8) F_arbi(nkc)
real(8), allocatable :: P1sm(:,:,:,:,:), P2sm(:,:,:,:,:)

nx=NC(1); ny=NC(2); nz=NC(3)

uQ=9402
unt=9400
open(uQ,file='Qvalue.dat',status='unknown')
write(uQ,*) 'nk,la,Q'

Qvalue(:,:)=0    ! in cm^-1
Qvalue_N(:,:)=0
Qvalue_U(:,:)=0
tauinv_N(:,:)=0
tauinv_U(:,:)=0
indx=0

temp = tempk*k_b/(100*h_plank*c_light)   ! convert to cm^-1

!write(*,*) 'entering cal_Qm2...'

call cpu_time(cputim1)
write(*,*) 'entering calculate_tauRTA...'

nkI1_loop: do ibz=1,nibz
                i=mapinv(ibz)
                call get_k_info(kibz(:,ibz) ,NC,nk1,i1,j1,k1,g1,g2,g3,inside)
                if (nk1.ne.i) then
                   write(*,*)'calculate_tauRTA error: i.ne.nk1 ',i,nk1
                   stop
                endif

 ! call Pmatrix3(i)
     ! -> read V33.xxx.dat
     ! -> allocate P1,P2 subset
     ! -> calculate P1, P2
     ! -> close V33.xxx.dat
 if (mod(ibz,ksub_size) .eq. 1) then   ! start point of each v33.xxx.dat

    if(allocated(P1sm) .or. allocated(P2sm)) then
        deallocate(P1sm)
        deallocate(P2sm)
    endif

    indx=indx+1
    if (indx*ksub_size > nibz) then
        ksub_size2=nibz - (indx-1)*ksub_size
        ksubset(1)=ksub_size*(indx-1)+1    ! starting point of k subset
        ksubset(2)=nibz
    else
        ksub_size2=ksub_size
        ksubset(1)=ksub_size*(indx-1)+1
        ksubset(2)=ksub_size*indx
    endif
    k_shift=ksub_size*(indx-1)    ! P1 of n1,n2,l1,l2,l3 will be saved in P1(n1-k_shift,n2,l1,l2,l3).

!    call Pmatrix4(indx,ksub_size2,ksubset,nk,kp,ndn,temp,k_shift)
    call cpu_time(cputim0)
    cputim0_wall = wtime()

    allocate(P1sm(ksub_size2,nkc,ndn,ndn,ndn))
    allocate(P2sm(ksub_size2,nkc,ndn,ndn,ndn))

    P1sm=0d0
    P2sm=0d0

! read v33sq_delta.xxx.dat

    write(filename_temp,fmt="(a,i3.3,a)") "v33sq_delta.",indx,".dat"
    filename=trim(v3path)//filename_temp
    open(unt,file=filename,status='old',form='unformatted')
    read(unt) nv3_2


    do n1=1,ksub_size2
    do l1=1,ndn
    do n2=1,nkc
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
    write(ulog,*) indx,' read v3sq_delta.dat done, TIME=',(cputim1-cputim0),'WALL TIME=',(cputim1_wall-cputim0_wall)

  endif

! Here call calc_tet
! Calculate Q_N, Q_U, Q_tot

    laI1_loop: do ii=1,ndn

    write(ulog,*) 'i,ii=',i,ii

    omega1=frequency(i,ii)        ! set fixed omega1 in delta(omega1-omega2+-omega3)

           laF2_loop: do jj=1,ndn
           laF3_loop: do kk=1,ndn


           ! process 1
           ! set arbitrary function here. (distribution function...)
           ! call eigen_tet using w2-w3
           ! call weight_tet

           ! set arbitrary function, F_arbi(j)
           do j=1,nkc   ! for set arbitrary function
                q1=kpc(:,i)
                q2=kpc(:,j)

                q3=-q1-q2
                q3_proc1=q1+q2
                q3_proc2=q1-q2

                call get_k_info(q2,NC,nk2,i2,j2,k2,g1,g2,g3,inside)
                if (nk2.ne.j) then
                   write(*,*)'calculate_tauRTA error: j.ne.nk2 ',j,nk2
                   stop
                endif
                call get_k_info(q3_proc1,NC,nk3_proc1,i3,j3,k3,g1,g2,g3,inside1)
                call get_k_info(q3_proc2,NC,nk3_proc2,i3,j3,k3,g1,g2,g3,inside2)


                F_arbi(j)=dist(i,ii)*dist(j,jj)*(dist(nk3_proc1,kk)+1)   ! no 2pi at the front since delta is in linear frequency.
                P1(ibz,j,ii,jj,kk) = P1sm(ibz-k_shift,j,ii,jj,kk) * F_arbi(j)

                F_arbi(j)=dist(i,ii)*(dist(j,jj)+1)*(dist(nk3_proc2,kk)+1)
                P2(ibz,j,ii,jj,kk) = P2sm(ibz-k_shift,j,ii,jj,kk) * F_arbi(j)

! K1   Qvalue(ibz,ii)=Qvalue(ibz,ii) + P1(ibz,j,ii,jj,kk) + 0.5*P2(ibz,j,ii,jj,kk)
                Qvalue(ibz,ii)=Qvalue(ibz,ii) + ( P1(ibz,j,ii,jj,kk)  &
                &              + 0.5*P2(ibz,j,ii,jj,kk) )

                if (inside1 .eq. 1) then    ! normal process in coalescence process
                    Qvalue_N(ibz,ii)=Qvalue_N(ibz,ii) + P1(ibz,j,ii,jj,kk)
                else
                    Qvalue_U(ibz,ii)=Qvalue_U(ibz,ii) + P1(ibz,j,ii,jj,kk)
                endif

                if (inside2 .eq. 1) then    ! normal process in decay process
                    Qvalue_N(ibz,ii)=Qvalue_N(ibz,ii) + 0.5*P2(ibz,j,ii,jj,kk)
                else
                    Qvalue_U(ibz,ii)=Qvalue_U(ibz,ii) + 0.5*P2(ibz,j,ii,jj,kk)
                endif

           enddo


           enddo laF3_loop
           enddo laF2_loop

         !  do jj=1,ndn
         !     do j=1,nkc

           !write(*,*)Qvalue_N(ibz,ii),Piso(ibz,j,ii,jj)
           Qvalue(ibz,ii)=Qvalue(ibz,ii) + Piso(ibz,ii)
           Qvalue_N(ibz,ii)=Qvalue_N(ibz,ii) + Piso(ibz,ii)
           Qvalue_U(ibz,ii)=Qvalue_U(ibz,ii) + Piso(ibz,ii)

           !write(*,*)ibz,j,ii,jj, Qvalue_N(ibz,ii),Piso(ibz,j,ii,jj)
          !    enddo
          !enddo


           nb1=dist(i,ii)

           write(uQ,*) ibz,ii,Qvalue(ibz,ii),Qvalue(ibz,ii)/(nb1*(nb1+1))

! tauinv_N in 1/sec
           tauinv_N(ibz,ii) = Qvalue_N(ibz,ii)/(nb1*(nb1+1)) * c_light*100
           tauinv_U(ibz,ii) = Qvalue_U(ibz,ii)/(nb1*(nb1+1)) * c_light*100
           tauinv_tot(ibz,ii) = tauinv_N(ibz,ii) + tauinv_U(ibz,ii)



    enddo laI1_loop



    call cpu_time(cputim2)
    write(*,*) indx,'th ksubset for Pmatrix3 done...', cputim2-cputim1


enddo nkI1_loop



close(uQ)

7 format(5(i6),2(2x,g14.8))

!if(indx2-1 .ne. num_pfiles) then
!   write(ulog,*) 'error in reading col_matrix.dat during cal_Qm...'
!   stop
!endif


end subroutine calculate_tauRTA_iso

!=========================================================






!==========================================================
subroutine syop_f2(s,ff,fs2) !!Saf!!
!using inverse symetry operations to transform a vector
use constants
use exactBTE2
use params
use lattice
use ios
use mod_ksubset2
use phi3
use kpoints
use atoms_force_constants
implicit none

integer s, ier
real(8) ff(3), fs2(3), tempi(3,3)


  call xmatinv(op_kmatrix(:,:,s),tempi,ier)

  fs2(:) = MATMUL(tempi,ff(:))

end subroutine syop_f2

!==========================================================
