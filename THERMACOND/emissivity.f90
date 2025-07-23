!==================================================
subroutine renormalizedphononfrequency

use constants
use exactBTE2
use params
use lattice
use kpoints

implicit none

integer i,j,nk,ndn
real(8) temp,tempk,omega,ise,omega_re
complex(8) nselfs,uselfs
real(8) q(3)

nk_loop: do i=1,nk  !  in the IBZ

           q=kibz(:,i)
   l_loop: do j=1,ndn

             omega=frequency(mapinv(i),j)
             !nb=dist(mapinv(i),ii)        
         
             call function_self_w3(q,j,omega,temp,nselfs,uselfs)
             ise = aimag(nselfs+uselfs)
             
             omega_re = sqrt(omega**2 + 2* omega * ise)              
 
             enddo l_loop
          enddo nk_loop

end subroutine renormalizedphononfrequency
!=====================================================

!============================================================
 subroutine dielectric(n,eival,eivec,mesh,om,epsilon0)
! takes the eigenvalues and eigenvectors at the gamma point and calculates the susceptibility
 use ios
 use atoms_force_constants, only : natoms0, atom0 !natom_prim_cell, atom0
 use lattice, only: volume_r !volume_r0
 use params
 use born, only : epsil
 use om_dos, only : width
 use constants, only : ee,pi,eps0,cnst,ci !r15,ee,pi,eps0,cnst,ci
 use kpoints, only : normald
 implicit none
 integer, intent(in):: n,mesh
 real(8),intent(in) :: om(mesh),eival(n)
 complex(8),intent(in) :: eivec(n,n)
 real(8),intent(out) :: epsilon0(3,3)
 complex(8) chi(3,3),chinormal !,intent(out) :: chi(3,3,mesh)
 real(8) d1,d2,z(n,3),reflectivity
 integer i,j,k,b,al,be,ga,de,tau,taup

 if(maxval(abs(aimag(eivec))).gt.1d-6*maxval(abs(eivec))) then
    write(ulog,*)'DIELECTRIC: eivec at gamma has imaginary parts '
    call write_out(ulog,' imag(eivec) ',aimag(eivec))
    call write_out(   6,' imag(eivec) ',aimag(eivec))
    stop
 endif

! calculate Born chargs in the eivec basis
 write(ulog,2)'#================== IR intensities versus  frequency(1/cm) ====================='
 z=0
 do b=4,n
 do de=1,3
    do ga=1,3
    do tau =1,natoms0  !natom_prim_cell
       j=ga+3*(tau-1)
       z(b,de)=z(b,de) + atom0(tau)%charge(de,ga)*eivec(j,b) / sqrt(atom0(tau)%mass)
    enddo
    enddo
    write(ulog,2) sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:))
 enddo
    write(*,*)'dot zz', sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:)) ,z(b,:)
 enddo
 write(ulog,2)'#=============================================================================='

! calculate and write susceptibility on a frequency mesh
 open(345,file='chi_real.dat')
 open(346,file='chi_imag.dat')
 write(345,*)'# omega(1/cm) , chi_xx , chi_xy , chi_xz  , chi_yx  ,  chi_yy  ,  chi_yz  , chi_zx , chi_zy , chi_zz , R'
 write(346,*)'# omega(1/cm) , chi_xx ,  chi_xy  ,  chi_xz  ,  chi_yx  ,  chi_yy  ,  chi_yz  ,  chi_zx  , chi_zy  , chi_zz'
 do i=1,mesh
    chi=0
    do al=1,3
    do be=1,3
    do b=4,n
       chi(al,be)=chi(al,be)+ z(b,al)*z(b,be)  / &
 &                            (-om(i)*om(i)+eival(b)*cnst*cnst-ci*width*width)

    write(*,*)'i,al,be,b,n',i,al,be,be,n,om(i),width,(z(b,al)*z(b,be)),chi(al,be)
    enddo
    enddo
    enddo
    chi=chi/volume_r*ee/eps0*1d10 * cnst*cnst
    if(i.eq.1) epsilon0=epsil+real(chi)
    chinormal=dot_product(normald,matmul((chi+epsil),normald))
    reflectivity= abs((sqrt(chinormal)-1)/(sqrt(chinormal)+1))**2
    write(345,2)om(i),real(chi),reflectivity
    write(346,2)om(i),aimag(chi)
 enddo


2 format(999(1x,g11.4))

 end subroutine dielectric
!===============
!============================================================
 subroutine dielectric_se(n,eival,eivec,mesh,om,tm,epsilon0)
! takes the eigenvalues and eigenvectors at the gamma point and calculates the susceptibility
 use ios
 use atoms_force_constants, only : natoms0, atom0  !natom_prim_cell, atom0
 use lattice, only: volume_r !volume_r0
 use params
 use born, only : epsil
 use om_dos, only : width
 use constants, only : ee,pi,eps0,cnst,ci !r15,ee,pi,eps0,cnst,ci
 use kpoints, only : normald
 implicit none
 integer, intent(in):: n,mesh
 real(8),intent(in) :: om(mesh),eival(n)
 complex(8),intent(in) :: eivec(n,n)
 real(8),intent(out) :: epsilon0(3,3)
 complex(8) chi(3,3),chinormal,eps(3,3) !,intent(out) :: chi(3,3,mesh)
 real(8) d1,d2,z(n,3),reflectivity, emissivity,ise,rse,tm,qe(3)
 complex(8) nselfs,uselfs
 integer i,j,k,b,al,be,ga,de,tau,taup!,qe
 

 if(maxval(abs(aimag(eivec))).gt.1d-6*maxval(abs(eivec))) then
    write(ulog,*)'DIELECTRIC: eivec at gamma has imaginary parts '
    call write_out(ulog,' imag(eivec) ',aimag(eivec))
    call write_out(   6,' imag(eivec) ',aimag(eivec))
    stop
 endif

 write(*,*) "volume=",volume_r
 write(*,*) "EIVAL, EIVEC are :"
 do tau=1,n
    write(6,4) eival(tau)
    do j=1,n
       write(6,4) eivec(j,tau)
    enddo
 enddo
 write(*,*) "BORN CHARGES are :"
 do tau=1,natoms0 
   do i=1,3
    write(6,4) atom0(tau)%charge(i,:)
   enddo
 enddo
4 format(99(1x,f9.4))










    !write(*,*)'b,de',b,de!,ga,tau!,eivec(j,b)!,atom0(tau)%charge(de,ga),sqrt(atom0(tau)%mass)
! calculate Born chargs in the eivec basis
 write(ulog,2)'#================== IR intensities versus  frequency(1/cm) ====================='
 z=0
 do b=4,n
 do de=1,3
    !write(*,*)'b,de',b,de
    do ga=1,3
!    write(*,*)'b,de,ga',b,de,ga, natoms0 !natom_prim_cell
    do tau =1,natoms0 !natom_prim_cell
       !write(*,*)'b,de,ga,tau',b,de,ga,tau
       j=ga+3*(tau-1)
       z(b,de)=z(b,de) + atom0(tau)%charge(de,ga)*eivec(j,b) / sqrt(atom0(tau)%mass)
    !write(*,*)'b,de',b,de!,ga,tau!,eivec(j,b)!,atom0(tau)%charge(de,ga),sqrt(atom0(tau)%mass)
    enddo
    enddo
    !write(*,*)'b,de',b,de,z(b,de)
   ! write(ulog,2) sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:))
 enddo
 !M(b,:) = z(b,:)/(2*sqrt(eival(b))*cnst) !?????
 write(*,*)'dot zz', sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:)) ,z(b,:)
enddo
 write(ulog,2)'#=============================================================================='

! calculate and write susceptibility on a frequency mesh
 open(345,file='chi_real.dat')
 open(346,file='chi_imag.dat')
 write(345,*)'# omega(1/cm) , chi_xx , chi_xy , chi_xz  , chi_yx  ,  chi_yy  ,  chi_yz  , chi_zx , chi_zy , chi_zz , R'
 write(346,*)'# omega(1/cm) , chi_xx ,  chi_xy  ,  chi_xz  ,  chi_yx  ,  chi_yy  ,  chi_yz  ,  chi_zx  , chi_zy  , chi_zz'
 do i=1,mesh
    chi=0
    do al=1,3
    do be=1,3
    do b=4,n
!    write(*,*)'i,al,be,b,n',i,al,be,be,n

       chi(al,be)=chi(al,be)+ z(b,al)*z(b,be)  / &
 &                            (-om(i)*om(i)+eival(b)*cnst*cnst-ci*width*width)

    write(*,*)'i,al,be,b,n',i,al,be,be,n,om(i),width,(z(b,al)*z(b,be)),chi(al,be)
    enddo
    enddo
    enddo
    chi=chi/volume_r*ee/eps0*1d10 * cnst*cnst
    if(i.eq.1) epsilon0=epsil+real(chi)
    chinormal=dot_product(normald,matmul((chi+epsil),normald))
    reflectivity= abs((sqrt(chinormal)-1)/(sqrt(chinormal)+1))**2
    
    emissivity = 1.0d0 - reflectivity

    write(345,2)om(i),real(chi),reflectivity, emissivity
    write(346,2)om(i),aimag(chi)
 enddo
stop
! calculate and write susceptibility on a frequency mesh
 open(3451,file='chis_realse.dat')
 open(3461,file='chis_imagse.dat')
 write(3451,*)'# omega(1/cm) , chi_xx , chi_xy , chi_xz  , chi_yx  ,  chi_yy  ,  chi_yz  , chi_zx , chi_zy , chi_zz , R'
 write(3461,*)'# omega(1/cm) , chi_xx ,  chi_xy  ,  chi_xz  ,  chi_yx  ,  chi_yy  ,  chi_yz  ,  chi_zx  , chi_zy  , chi_zz'
 open(3471,file='eps_realse.dat')
 open(3481,file='eps_imagse.dat')
 write(3471,*)'# omega(1/cm) , eps_xx , eps_xy , eps_xz  , eps_yx  ,  eps_yy  ,  eps_yz  , eps_zx , eps_zy , eps_zz'
 write(3481,*)'# omega(1/cm) , eps_xx ,  eps_xy  ,  eps_xz  ,  eps_yx  ,  eps_yy  ,  eps_yz  ,  eps_zx  , eps_zy  , eps_zz'
 
 chi=0
 do i=1,mesh
    chi=0
       qe=0d0
!       write(*,*)'i---------------',i
       call function_self_w3(qe,b,om(i),tm,nselfs,uselfs)
       ise = aimag(nselfs+uselfs)
       rse = real(nselfs+uselfs)
       write(*,*)'ise,rse', ise, rse
write(*,*)'i---------------',i
    do al=1,3
    do be=1,3
    do b=4,n
!       qe=0d0
!       write(*,*)'i---------------',i
!       call function_self_w3(qe,b,om(i),tm,nselfs,uselfs)
!       ise = aimag(nselfs+uselfs)
!       rse = real(nselfs+uselfs)
!       write(*,*)'ise,rse', ise, rse
       chi(al,be)=chi(al,be)+ z(b,al)*z(b,be)  / &      !M(b,al) * M(b,be) ?????
 &                            (om(i)*om(i)-eival(b)*cnst*cnst-2*sqrt(eival(b))*cnst*(rse-ci*ise))
    
    !   eps(al,be)=eps_inf* delta_l(al,be) + chi(al,be)
      write(*,*)'i-end-------------',i,al,be,b,z(b,al),z(b,be),chi(al,be)
    enddo
    enddo
    enddo
    chi=chi/volume_r*ee/eps0*1d10 * cnst*cnst
    eps=chi+epsil
    !write(*,*)'i-end-------------',i, chi, eps
    if(i.eq.1) epsilon0=epsil+real(chi)
    chinormal=dot_product(normald,matmul((chi+epsil),normald))
    !write(*,*)' chinormal',i,chinormal
    reflectivity= abs((sqrt(chinormal)-1)/(sqrt(chinormal)+1))**2
    
    !ieps=aimag(eps)
    !reps=real(eps)
     
    emissivity = 1.0d0 - reflectivity
   
    write(*,*)'i,om,emissivity',i,om(i),real(chi),reflectivity, emissivity 
    write(3451,2)om(i),real(chi),reflectivity, emissivity
    write(3461,2)om(i),aimag(chi)
    write(3471,2)om(i),real(eps)
    write(3481,2)om(i),aimag(eps)
 enddo


2 format(999(1x,g11.4))

 end subroutine dielectric_se
!===============




!============================================================
 subroutine dielectric_sempi(wmsubset,n,eival,eivec,mesh,om,tm,epsilon0)
! takes the eigenvalues and eigenvectors at the gamma point and calculates the susceptibility
 use ios
 use atoms_force_constants, only : natoms0, atom0  !natom_prim_cell, atom0
 use lattice, only: volume_r !volume_r0
 use params
 use born, only : epsil
 use om_dos, only : width
 use constants, only : ee,pi,eps0,cnst,ci !r15,ee,pi,eps0,cnst,ci
 use kpoints, only : normald
 use phi3, only :rsempi_em,isempi_em,reflectivitympi,chimpi,nmpi,kmpi,reflectivitympink
 implicit none
 integer, intent(in):: n,mesh
 real(8),intent(in) :: om(mesh),eival(n)
 complex(8),intent(in) :: eivec(n,n)
 real(8),intent(out) :: epsilon0(3,3)
 complex(8) chi(3,3),chinormal,eps(3,3) !,intent(out) :: chi(3,3,mesh)
 real(8) d1,d2,z(n,3),reflectivity, emissivity,ise,rse,tm,qe(3),eps1,eps2
 complex(8) nselfs,uselfs
 integer i,j,k,b,al,be,ga,de,tau,taup,ii,s!,qe
 integer wmsubset(2) ,wmsubsetsize3,wm_shift

 !if(maxval(abs(aimag(eivec))).gt.1d-6*maxval(abs(eivec))) then
 !   write(ulog,*)'DIELECTRIC: eivec at gamma has imaginary parts '
 !   call write_out(ulog,' imag(eivec) ',aimag(eivec))
 !   call write_out(   6,' imag(eivec) ',aimag(eivec))
 !   stop
 !endif

! write(*,*) "volume=",volume_r
! write(*,*) "EIVAL, EIVEC are :"
! do tau=1,n
!    write(6,4) eival(tau)
!    do j=1,n
!       write(6,4) eivec(j,tau)
!    enddo
! enddo
! write(*,*) "BORN CHARGES are :"
! do tau=1,natoms0 
!   do i=1,3
!    write(6,4) atom0(tau)%charge(i,:)
!   enddo
! enddo
!4 format(99(1x,f9.4))


    !write(*,*)'b,de',b,de!,ga,tau!,eivec(j,b)!,atom0(tau)%charge(de,ga),sqrt(atom0(tau)%mass)
! calculate Born chargs in the eivec basis
 write(ulog,2)'#================== IR intensities versus  frequency(1/cm) ====================='
 z=0
 do b=4,n
 do de=1,3
    !write(*,*)'b,de',b,de
    do ga=1,3
!    write(*,*)'b,de,ga',b,de,ga, natoms0 !natom_prim_cell
    do tau =1,natoms0 !natom_prim_cell
       !write(*,*)'b,de,ga,tau',b,de,ga,tau
       j=ga+3*(tau-1)
       z(b,de)=z(b,de) + atom0(tau)%charge(de,ga)*eivec(j,b) / sqrt(atom0(tau)%mass)
    !write(*,*)'b,de',b,de!,ga,tau!,eivec(j,b)!,atom0(tau)%charge(de,ga),sqrt(atom0(tau)%mass)
    enddo
    enddo
    !write(*,*)'b,de',b,de,z(b,de)
   ! write(ulog,2) sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:))
 enddo
 !M(b,:) = z(b,:)/(2*sqrt(eival(b))*cnst) !?????
 !write(*,*)'dot zz', sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:)) ,z(b,:)
enddo
 write(ulog,2)'#=============================================================================='

! calculate and write susceptibility on a frequency mesh
! open(345,file='chi_real.dat')
! open(346,file='chi_imag.dat')
! write(345,*)'# omega(1/cm) , chi_xx , chi_xy , chi_xz  , chi_yx  ,  chi_yy  ,  chi_yz  , chi_zx , chi_zy , chi_zz , R'
! write(346,*)'# omega(1/cm) , chi_xx ,  chi_xy  ,  chi_xz  ,  chi_yx  ,  chi_yy  ,  chi_yz  ,  chi_zx  , chi_zy  , chi_zz'
! do i=1,mesh
!    chi=0
!    do al=1,3
!    do be=1,3
!    do b=4,n
!    write(*,*)'i,al,be,b,n',i,al,be,be,n

!       chi(al,be)=chi(al,be)+ z(b,al)*z(b,be)  / &
! &                            (-om(i)*om(i)+eival(b)*cnst*cnst-ci*width*width)

!    write(*,*)'i,al,be,b,n',i,al,be,be,n,om(i),width,(z(b,al)*z(b,be)),chi(al,be)
!    enddo
!    enddo
!    enddo
!    chi=chi/volume_r*ee/eps0*1d10 * cnst*cnst
!    if(i.eq.1) epsilon0=epsil+real(chi)
!    chinormal=dot_product(normald,matmul((chi+epsil),normald))
!    reflectivity= abs((sqrt(chinormal)-1)/(sqrt(chinormal)+1))**2
    
!    emissivity = 1.0d0 - reflectivity

!    write(345,2)om(i),real(chi),reflectivity, emissivity
!    write(346,2)om(i),aimag(chi)
! enddo
!stop
! calculate and write susceptibility on a frequency mesh
! open(3451,file='chis_realse.dat')
! open(3461,file='chis_imagse.dat')
! write(3451,*)'# omega(1/cm) , chi_xx , chi_xy , chi_xz  , chi_yx  ,  chi_yy  ,  chi_yz  , chi_zx , chi_zy , chi_zz , R'
! write(3461,*)'# omega(1/cm) , chi_xx ,  chi_xy  ,  chi_xz  ,  chi_yx  ,  chi_yy  ,  chi_yz  ,  chi_zx  , chi_zy  , chi_zz'
! open(3471,file='eps_realse.dat')
! open(3481,file='eps_imagse.dat')
! write(3471,*)'# omega(1/cm) , eps_xx , eps_xy , eps_xz  , eps_yx  ,  eps_yy  ,  eps_yz  , eps_zx , eps_zy , eps_zz'
! write(3481,*)'# omega(1/cm) , eps_xx ,  eps_xy  ,  eps_xz  ,  eps_yx  ,  eps_yy  ,  eps_yz  ,  eps_zx  , eps_zy  , eps_zz'
 

wmsubsetsize3=(wmsubset(2)-wmsubset(1))+1
wm_shift=wmsubset(1)-1

qe=0d0
qe(3)=0.0001
  
 chimpi=0
 isempi_em =0
 rsempi_em =0
 reflectivitympi =0
 do i=1,wmsubsetsize3 
    !chi=0
       ii=i+wm_shift
    write(*,*)'i---------------',i 
   do b=4,n
      ! qe=0d0
       s=b-3
!       write(*,*)'i---------------',i
       call function_self_w3(qe,b,om(ii),tm,nselfs,uselfs)
       isempi_em(i,s) = aimag(nselfs+uselfs)
       rsempi_em(i,s) = real(nselfs+uselfs)
       write(*,*)'i,ise,rse',i,ii,om(ii),b, isempi_em(i,s), rsempi_em(i,s)
   
    enddo
 enddo


 do i=1,wmsubsetsize3 
    !chi=0
       ii=i+wm_shift
write(*,*)'i---------------',i
    do al=1,3
    do be=1,3
    do b=4,n
       qe=0d0
       s=b-3
!       write(*,*)'i---------------',i
!       call function_self_w3(qe,b,om(ii),tm,nselfs,uselfs)
!       isempi_em(i,s) = aimag(nselfs+uselfs)
!       rsempi_em(i,s) = real(nselfs+uselfs)
!       write(*,*)'i,ise,rse',i,ii,om(ii),b, isempi_em(i,s), rsempi_em(i,s)


!       call function_self_w3(qe,b,om(i),tm,nselfs,uselfs)
!       ise = aimag(nselfs+uselfs)
!       rse = real(nselfs+uselfs)
!       write(*,*)'ise,rse', ise, rse
       chimpi(i,al,be)=chimpi(i,al,be)- z(b,al)*z(b,be)  / &      !M(b,al) * M(b,be) ?????
 &                            (om(ii)*om(ii)-eival(b)*cnst*cnst-2*sqrt(eival(b))*cnst*(rsempi_em(i,s)-ci*isempi_em(i,s)))
    
    !   eps(al,be)=eps_inf* delta_l(al,be) + chi(al,be)
      !write(*,*)'i-end-------------',i,al,be,b,z(b,al),z(b,be),chi(al,be)
    enddo
    enddo
    enddo
    chimpi(i,:,:)=chimpi(i,:,:)/volume_r*ee/eps0*1d10 * cnst*cnst
!    epsmpi(i,:,:)=chimpi(i,:,:)+epsil
    !write(*,*)'i-end-------------',i, chi, eps
    if(i.eq.1) epsilon0=epsil+real(chimpi(i,:,:))
    chinormal=dot_product(normald,matmul((chimpi(i,:,:)+epsil),normald))
    eps1 = real(chinormal)
    eps2 = aimag(chinormal)
    nmpi(i) = sqrt((eps1 + sqrt(eps1**2 + eps2**2))/2)
    kmpi(i) = sqrt((- eps1 + sqrt(eps1**2 + eps2**2))/2)
    reflectivitympink(i) = ((nmpi(i) - 1)**2 + kmpi(i)**2)/((nmpi(i) + 1)**2 + kmpi(i)**2) 
    !write(*,*)' chinormal',i,chinormal
    reflectivitympi(i)= abs((sqrt(chinormal)-1)/(sqrt(chinormal)+1))**2
    


     
    !ieps=aimag(eps)
    !reps=real(eps)
     
    !emissivity = 1.0d0 - reflectivity
    
 !   write(*,*)'i,om,emissivityinsideeee',i,om(ii),real(chimpi(i,:,:)),reflectivitympi(i), 1-reflectivitympi(i)  
    write(*,*)'i,om,n,k,R',i,om(ii),nmpi(i), kmpi(i), reflectivitympink(i)
    write(*,*)'i,om,insideeeeee', i,om(ii), aimag(chimpi(i,:,:))
    write(*,*)'i,om,emissivity',i,om(ii),real(chimpi(i,:,:))!,reflectivity, emissivity 
    !write(3451,2)om(i),real(chi),reflectivity, emissivity
    !write(3461,2)om(i),aimag(chi)
    !write(3471,2)om(i),real(eps)
    !write(3481,2)om(i),aimag(eps)
 enddo


2 format(999(1x,g11.4))

 end subroutine dielectric_sempi
!===============

!============================================================
 subroutine dielectric_sempi2(wmsubset,n,eival,eivec,mesh,om,tm,epsilon0)
! takes the eigenvalues and eigenvectors at the gamma point and calculates the susceptibility
 use ios
 use atoms_force_constants, only : natoms0, atom0  !natom_prim_cell, atom0
 use lattice, only: volume_r !volume_r0
 use params
 use born, only : epsil
 use om_dos, only : width
 use constants, only : ee,pi,eps0,cnst,ci !r15,ee,pi,eps0,cnst,ci
 use kpoints, only : normald
 use phi3, only :rsempi_em,isempi_em,reflectivitympi,chimpi,nmpi,kmpi,reflectivitympink
 implicit none
 integer, intent(in):: n,mesh
 real(8),intent(in) :: om(mesh),eival(n)
 complex(8),intent(in) :: eivec(n,n)
 real(8),intent(out) :: epsilon0(3,3)
 complex(8) chi(3,3),chinormal,eps(3,3) !,intent(out) :: chi(3,3,mesh)
 real(8) d1,d2,z(n,3),reflectivity, emissivity,ise,rse,tm,qe(3),eps1,eps2,rekk,wes
 complex(8) nselfs,uselfs
 integer i,j,k,b,al,be,ga,de,tau,taup,ii,s!,qe
 integer wmsubset(2) ,wmsubsetsize3,wm_shift
 real(8), allocatable :: imgse(:), xkk(:), omm(:) 
 
 if(maxval(abs(aimag(eivec))).gt.1d-6*maxval(abs(eivec))) then
    write(ulog,*)'DIELECTRIC: eivec at gamma has imaginary parts '
    call write_out(ulog,' imag(eivec) ',aimag(eivec))
    call write_out(   6,' imag(eivec) ',aimag(eivec))
    stop
 endif

! write(*,*) "volume=",volume_r
! write(*,*) "EIVAL, EIVEC are :"
! do tau=1,n
!    write(6,4) eival(tau)
!    do j=1,n
!       write(6,4) eivec(j,tau)
!    enddo
! enddo
! write(*,*) "BORN CHARGES are :"
! do tau=1,natoms0 
!   do i=1,3
!    write(6,4) atom0(tau)%charge(i,:)
!   enddo
! enddo
!4 format(99(1x,f9.4))


    !write(*,*)'b,de',b,de!,ga,tau!,eivec(j,b)!,atom0(tau)%charge(de,ga),sqrt(atom0(tau)%mass)
! calculate Born chargs in the eivec basis
 write(ulog,2)'#================== IR intensities versus  frequency(1/cm) ====================='
 z=0
 do b=4,n
 do de=1,3
    !write(*,*)'b,de',b,de
    do ga=1,3
!    write(*,*)'b,de,ga',b,de,ga, natoms0 !natom_prim_cell
    do tau =1,natoms0 !natom_prim_cell
       !write(*,*)'b,de,ga,tau',b,de,ga,tau
       j=ga+3*(tau-1)
       z(b,de)=z(b,de) + atom0(tau)%charge(de,ga)*eivec(j,b) / sqrt(atom0(tau)%mass)
    !write(*,*)'b,de',b,de!,ga,tau!,eivec(j,b)!,atom0(tau)%charge(de,ga),sqrt(atom0(tau)%mass)
    enddo
    enddo
    !write(*,*)'b,de',b,de,z(b,de)
   ! write(ulog,2) sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:))
 enddo
 !M(b,:) = z(b,:)/(2*sqrt(eival(b))*cnst) !?????
 !write(*,*)'dot zz', sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:)) ,z(b,:)
enddo
 write(ulog,2)'#=============================================================================='

! calculate and write susceptibility on a frequency mesh
! open(345,file='chi_real.dat')
! open(346,file='chi_imag.dat')
! write(345,*)'# omega(1/cm) , chi_xx , chi_xy , chi_xz  , chi_yx  ,  chi_yy  ,  chi_yz  , chi_zx , chi_zy , chi_zz , R'
! write(346,*)'# omega(1/cm) , chi_xx ,  chi_xy  ,  chi_xz  ,  chi_yx  ,  chi_yy  ,  chi_yz  ,  chi_zx  , chi_zy  , chi_zz'
! do i=1,mesh
!    chi=0
!    do al=1,3
!    do be=1,3
!    do b=4,n
!    write(*,*)'i,al,be,b,n',i,al,be,be,n

!       chi(al,be)=chi(al,be)+ z(b,al)*z(b,be)  / &
! &                            (-om(i)*om(i)+eival(b)*cnst*cnst-ci*width*width)

!    write(*,*)'i,al,be,b,n',i,al,be,be,n,om(i),width,(z(b,al)*z(b,be)),chi(al,be)
!    enddo
!    enddo
!    enddo
!    chi=chi/volume_r*ee/eps0*1d10 * cnst*cnst
!    if(i.eq.1) epsilon0=epsil+real(chi)
!    chinormal=dot_product(normald,matmul((chi+epsil),normald))
!    reflectivity= abs((sqrt(chinormal)-1)/(sqrt(chinormal)+1))**2
    
!    emissivity = 1.0d0 - reflectivity

!    write(345,2)om(i),real(chi),reflectivity, emissivity
!    write(346,2)om(i),aimag(chi)
! enddo
!stop
! calculate and write susceptibility on a frequency mesh
! open(3451,file='chis_realse.dat')
! open(3461,file='chis_imagse.dat')
! write(3451,*)'# omega(1/cm) , chi_xx , chi_xy , chi_xz  , chi_yx  ,  chi_yy  ,  chi_yz  , chi_zx , chi_zy , chi_zz , R'
! write(3461,*)'# omega(1/cm) , chi_xx ,  chi_xy  ,  chi_xz  ,  chi_yx  ,  chi_yy  ,  chi_yz  ,  chi_zx  , chi_zy  , chi_zz'
! open(3471,file='eps_realse.dat')
! open(3481,file='eps_imagse.dat')
! write(3471,*)'# omega(1/cm) , eps_xx , eps_xy , eps_xz  , eps_yx  ,  eps_yy  ,  eps_yz  , eps_zx , eps_zy , eps_zz'
! write(3481,*)'# omega(1/cm) , eps_xx ,  eps_xy  ,  eps_xz  ,  eps_yx  ,  eps_yy  ,  eps_yz  ,  eps_zx  , eps_zy  , eps_zz'
 

wmsubsetsize3=(wmsubset(2)-wmsubset(1))+1
wm_shift=wmsubset(1)-1

!call calculate_v35(ksubibz,ndn,nib,eivlibz,eivcibz,nk,eival,eivec)


 chimpi=0
 isempi_em =0
 rsempi_em =0
 reflectivitympi =0
!ss=0
 do i=1,wmsubsetsize3 
    ii=i+wm_shift
write(*,*)'i---------------',i
       do b=4,n
       qe=0d0
       s=b-3
 !      ss=ss+1
!       write(*,*)'i---------------',i
       !!!!!!call function_self_w3(qe,b,om(ii),tm,nselfs,uselfs)
       wes=(om(ii)/cnst )**2
       write(*,*)'wes=(om(ii)/cnst )**2',wes       
       call function_self_tetra(1,qe,b,wes,tm,nselfs,uselfs)      

       isempi_em(i,s) = aimag(nselfs+uselfs)
  !     xkk(ss) = om(ii)*om(ii)
  !     imgse(ss)=isempi_em(i,s)
!       rsempi_em(i,s) = real(nselfs+uselfs)
       write(*,*)'i,ise,rse',i,ii,om(ii),b, isempi_em(i,s)!, rsempi_em(i,s)
      enddo
enddo
  
allocate(imgse(wmsubsetsize3),xkk(wmsubsetsize3),omm(wmsubsetsize3))

do s = 1, n - 3
    do i = 1, wmsubsetsize3
        imgse(i) = isempi_em(i, s)
        xkk(i) = om(i)**2
    enddo

    do i = 1, wmsubsetsize3
        ii = i + wm_shift
        omm(i) = om(ii)**2
        call kk_im2re(omm(i), wmsubsetsize3, imgse, xkk, rekk)
        rsempi_em(i, s) = rekk
        write(*,*)'rsempi_em(i, s) ',s,i,rsempi_em(i, s)

    enddo
enddo
deallocate(imgse,xkk)



 do i=1,wmsubsetsize3 
    !chi=0
       ii=i+wm_shift
write(*,*)'i---------------',i
    do al=1,3
    do be=1,3
    do b=4,n
       qe=0d0
       s=b-3
!       write(*,*)'i---------------',i
       !!!!!!call function_self_w3(qe,b,om(ii),tm,nselfs,uselfs)
        
 !      call function_self_tetra(1,qe,b,om(ii),tm,nselfs,uselfs)      

 !      isempi_em(i,s) = aimag(nselfs+uselfs)
 !      rsempi_em(i,s) = real(nselfs+uselfs)
 !      write(*,*)'i,ise,rse',i,ii,om(ii),b, isempi_em(i,s), rsempi_em(i,s)


!       call function_self_w3(qe,b,om(i),tm,nselfs,uselfs)
!       ise = aimag(nselfs+uselfs)
!       rse = real(nselfs+uselfs)
!       write(*,*)'ise,rse', ise, rse
       chimpi(i,al,be)=chimpi(i,al,be)- z(b,al)*z(b,be)  / &      !M(b,al) * M(b,be) ?????
 &                            (om(ii)*om(ii)-eival(b)*cnst*cnst-2*sqrt(eival(b))*cnst*(rsempi_em(i,s)-ci*isempi_em(i,s)))
    
    !   eps(al,be)=eps_inf* delta_l(al,be) + chi(al,be)
      !write(*,*)'i-end-------------',i,al,be,b,z(b,al),z(b,be),chi(al,be)
    enddo
    enddo
    enddo
    chimpi(i,:,:)=chimpi(i,:,:)/volume_r*ee/eps0*1d10 * cnst*cnst
!    epsmpi(i,:,:)=chimpi(i,:,:)+epsil
    !write(*,*)'i-end-------------',i, chi, eps
    if(i.eq.1) epsilon0=epsil+real(chimpi(i,:,:))
    chinormal=dot_product(normald,matmul((chimpi(i,:,:)+epsil),normald))
    eps1 = real(chinormal)
    eps2 = aimag(chinormal)
    nmpi(i) = sqrt((eps1 + sqrt(eps1**2 + eps2**2))/2)
    kmpi(i) = sqrt((- eps1 + sqrt(eps1**2 + eps2**2))/2)
    reflectivitympink(i) = ((nmpi(i) - 1)**2 + kmpi(i)**2)/((nmpi(i) + 1)**2 + kmpi(i)**2) 
    !write(*,*)' chinormal',i,chinormal
    reflectivitympi(i)= abs((sqrt(chinormal)-1)/(sqrt(chinormal)+1))**2
    


     
    !ieps=aimag(eps)
    !reps=real(eps)
     
    !emissivity = 1.0d0 - reflectivity
    
 !   write(*,*)'i,om,emissivityinsideeee',i,om(ii),real(chimpi(i,:,:)),reflectivitympi(i), 1-reflectivitympi(i)  
    write(*,*)'i,om,n,k,R',i,om(ii),nmpi(i), kmpi(i), reflectivitympink(i)
    write(*,*)'i,om,insideeeeee', i,om(ii), aimag(chimpi(i,:,:))
    write(*,*)'i,om,emissivity',i,om(ii),real(chimpi(i,:,:))!,reflectivity, emissivity 
    !write(3451,2)om(i),real(chi),reflectivity, emissivity
    !write(3461,2)om(i),aimag(chi)
    !write(3471,2)om(i),real(eps)
    !write(3481,2)om(i),aimag(eps)
 enddo


2 format(999(1x,g11.4))

 end subroutine dielectric_sempi2
!===============
!============================================================
 subroutine dielectric_setetra(n,eival,eivec,mesh,om,tm,epsilon0)
! takes the eigenvalues and eigenvectors at the gamma point and calculates the susceptibility
 use ios
 use atoms_force_constants, only : natoms0, atom0  !natom_prim_cell, atom0
 use lattice, only: volume_r !volume_r0
 use params
 use born, only : epsil
 use om_dos, only : width
 use constants, only : ee,pi,eps0,cnst,ci !r15,ee,pi,eps0,cnst,ci
 use kpoints, only : normald
 !use phi3, only :rsempi_em,isempi_em,reflectivitympi,chimpi,nmpi,kmpi,reflectivitympink
 implicit none
 integer, intent(in):: n,mesh
 real(8),intent(in) :: om(mesh),eival(n)
 complex(8),intent(in) :: eivec(n,n)
 real(8),intent(out) :: epsilon0(3,3)
 complex(8) chi(3,3),chinormal,eps(3,3),z(n,3) !,intent(out) :: chi(3,3,mesh)
 real(8) d1,d2,reflectivity, emissivity,ise,rse,tm,qe(3),eps1,eps2,rekk,wes
 complex(8) nselfs,uselfs
 real(8) ise_em(mesh,n-3),rse_em(mesh,n-3), n_es,k_es,reflectivity_nk
 integer i,j,k,b,al,be,ga,de,tau,taup,ii,s!,qe
! integer wmsubset(2) ,wmsubsetsize3,wm_shift
 real(8), allocatable :: imgse(:), xkk(:), omm(:) 
 
! if(maxval(abs(aimag(eivec))).gt.1d-6*maxval(abs(eivec))) then
!    write(ulog,*)'DIELECTRIC: eivec at gamma has imaginary parts '
!    call write_out(ulog,' imag(eivec) ',aimag(eivec))
!    call write_out(   6,' imag(eivec) ',aimag(eivec))
!    stop
! endif

open(3433,file='egnvec-zzz.dat')

! calculate Born chargs in the eivec basis
 write(ulog,2)'#================== IR intensities versus  frequency(1/cm) ====================='
 z=0
 do b=4,n
 do de=1,3
    !write(*,*)'b,de',b,de
    do ga=1,3
!    write(*,*)'b,de,ga',b,de,ga, natoms0 !natom_prim_cell
    do tau =1,natoms0 !natom_prim_cell
       !write(*,*)'b,de,ga,tau',b,de,ga,tau
       j=ga+3*(tau-1)
       z(b,de)=z(b,de) + atom0(tau)%charge(de,ga)*eivec(j,b) / sqrt(atom0(tau)%mass)
    write(3433,*)'b,de,ga,tau',b,de,ga,tau,eivec(j,b),atom0(tau)%charge(de,ga),sqrt(atom0(tau)%mass)
    enddo
    enddo
    write(3433,*)'#================================',b
    write(3433,*)'b,de,z(b,de)',b,de,z(b,de)
   ! write(ulog,2) sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:))
 enddo
 !M(b,:) = z(b,:)/(2*sqrt(eival(b))*cnst) !?????
 !write(*,*)'dot zz', sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:)) ,z(b,:)
enddo
 write(ulog,2)'#=============================================================================='

! calculate and write susceptibility on a frequency mesh
 open(3451,file='chis_realse.dat')
 open(3461,file='chis_imagse.dat')
 open(3471,file='optical-constantsee.dat')
 open(3481,file='selfenergy-imag.dat')
 open(3481,file='selfenergy-real.dat')
 open(3381,file='selfenergy.dat')

!call calculate_v35(ksubibz,ndn,nib,eivlibz,eivcibz,nk,eival,eivec)
qe=0d0
qe(3)=0.0001
 ise_em =0
 rse_em =0
 reflectivity =0
 do i=1,mesh !wmsubsetsize3 
    !ii=i+wm_shift
!write(*,*)'i---------------',i
       do b=4,n
      ! qe=0d0
       s=b-3
       
!       wes=(om(i)/cnst )**2
 !      write(*,*)'wes=(om(ii)/cnst )**2',wes       
       call function_self_tetra(1,qe,b,om(i),tm,nselfs,uselfs)      

       ise_em(i,s) = aimag(nselfs+uselfs)
       
       write(*,*)'i,ise,rse',i,s,om(i),b, ise_em(i,s)!, rsempi_em(i,s)
       write(3481,2)i,s,om(i),b, ise_em(i,s) 
     enddo
enddo
  
allocate(imgse(mesh),xkk(mesh),omm(mesh))

do s = 1, n - 3
    do i = 1, mesh 
        imgse(i) = ise_em(i, s)
        xkk(i) = om(i)**2
    enddo

    do i = 1,mesh  
        !ii = i + wm_shift
        omm(i) = om(i)**2
        call kk_im2re(omm(i), mesh, imgse, xkk, rekk)
        rse_em(i, s) = rekk
        write(3421,2) s,i,om(i),rse_em(i, s)

    enddo
enddo
deallocate(imgse,xkk,omm)


do i=1,mesh
   do b=4,n
    s=b-3
    write(3381,2)i,b,b-3,om(i),ise_em(i,s),rse_em(i, s)
    enddo
enddo



 do i=1,mesh !wmsubsetsize3 
    chi=0
    !   ii=i+wm_shift
write(*,*)'i---------------',i
    do al=1,3
    do be=1,3
    do b=4,n
       qe=0d0
       s=b-3


       chi(al,be)=chi(al,be)- z(b,al)*z(b,be)  / &      !M(b,al) * M(b,be) ?????
 &                            (om(i)*om(i)-eival(b)**2 -2*eival(b)*(rse_em(i,s)-ci*ise_em(i,s)))
 ! &                            (om(i)*om(i)-eival(b)*cnst*cnst-2*sqrt(eival(b))*cnst*(rse_em(i,s)-ci*ise_em(i,s)))
    
    enddo
    enddo
    enddo
    chi=chi/volume_r*ee/eps0*1d10 * cnst*cnst
    if(i.eq.1) epsilon0=epsil+real(chi)
    chinormal=dot_product(normald,matmul((chi+epsil),normald))
    eps1 = real(chinormal)
    eps2 = aimag(chinormal)
    n_es = sqrt((eps1 + sqrt(eps1**2 + eps2**2))/2)
    k_es = sqrt((- eps1 + sqrt(eps1**2 + eps2**2))/2)
    reflectivity_nk = ((n_es - 1)**2 + k_es**2)/((n_es + 1)**2 + k_es**2) 
    !write(*,*)' chinormal',i,chinormal
    reflectivity= abs((sqrt(chinormal)-1)/(sqrt(chinormal)+1))**2
    


     
    
 !   write(*,*)'i,om,emissivityinsideeee',i,om(ii),real(chimpi(i,:,:)),reflectivitympi(i), 1-reflectivitympi(i)  
    write(*,*)'i,om,n,k,R',i,om(i),n_es, k_es, reflectivity_nk,1-reflectivity_nk
    write(*,*)'i,om,insideeeeee', i,om(i), aimag(chi)
    write(*,*)'i,om,emissivity',i,om(i),real(chi)!,reflectivity, emissivity 
    write(3451,2)om(i),real(chi),reflectivity, 1-reflectivity
    write(3461,2)om(i),aimag(chi)
    write(3471,2)i,om(i),n_es, k_es, reflectivity_nk,1-reflectivity_nk
 enddo


2 format(999(1x,g11.4))

 end subroutine dielectric_setetra
!===============
!============================================================
 subroutine dielectric_setetrampirzero(n,eival,eivec,mesh,om,tm,ise_em,epsilon0)
! takes the eigenvalues and eigenvectors at the gamma point and calculates the susceptibility
 use ios
 use atoms_force_constants, only : natoms0, atom0  !natom_prim_cell, atom0
 use lattice, only: volume_r !volume_r0
 use params
 use born, only : epsil
 use om_dos, only : width
 use constants, only : ee,pi,eps0,cnst,ci !r15,ee,pi,eps0,cnst,ci
 use kpoints, only : normald
 !use phi3, only :rsempi_em,isempi_em,reflectivitympi,chimpi,nmpi,kmpi,reflectivitympink
 implicit none
 integer, intent(in):: n,mesh
 real(8),intent(in) :: om(mesh),eival(n)
 complex(8),intent(in) :: eivec(n,n)
 real(8),intent(out) :: epsilon0(3,3)
 complex(8) chi(3,3),chinormal,eps(3,3),z(n,3) !,intent(out) :: chi(3,3,mesh)
 real(8) d1,d2,reflectivity, emissivity,ise,rse,tm,qe(3),eps1,eps2,rekk,wes
 complex(8) nselfs,uselfs
 real(8) ise_em(mesh,n-3),rse_em(mesh,n-3), n_es,k_es,reflectivity_nk
 integer i,j,k,b,al,be,ga,de,tau,taup,ii,s!,qe
! integer wmsubset(2) ,wmsubsetsize3,wm_shift
 real(8), allocatable :: imgse(:), xkk(:), omm(:) 
 
! if(maxval(abs(aimag(eivec))).gt.1d-6*maxval(abs(eivec))) then
!    write(ulog,*)'DIELECTRIC: eivec at gamma has imaginary parts '
!    call write_out(ulog,' imag(eivec) ',aimag(eivec))
!    call write_out(   6,' imag(eivec) ',aimag(eivec))
!    stop
! endif

open(3433,file='egnvec-zzz.dat')

! calculate Born chargs in the eivec basis
 write(ulog,2)'#================== IR intensities versus  frequency(1/cm) ====================='
 z=0
 do b=4,n
 do de=1,3
    !write(*,*)'b,de',b,de
    do ga=1,3
!    write(*,*)'b,de,ga',b,de,ga, natoms0 !natom_prim_cell
    do tau =1,natoms0 !natom_prim_cell
       !write(*,*)'b,de,ga,tau',b,de,ga,tau
       j=ga+3*(tau-1)
       z(b,de)=z(b,de) + atom0(tau)%charge(de,ga)*eivec(j,b) / sqrt(atom0(tau)%mass)
    write(3433,*)'b,de,ga,tau',b,de,ga,tau,eivec(j,b),atom0(tau)%charge(de,ga),sqrt(atom0(tau)%mass)
    enddo
    enddo
    write(3433,*)'#================================',b
    write(3433,*)'b,de,z(b,de)',b,de,z(b,de)
   ! write(ulog,2) sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:))
 enddo
 !M(b,:) = z(b,:)/(2*sqrt(eival(b))*cnst) !?????
 !write(*,*)'dot zz', sqrt(eival(b))*cnst, dot_product(z(b,:),z(b,:)) ,z(b,:)
enddo
 write(ulog,2)'#=============================================================================='

! calculate and write susceptibility on a frequency mesh
 open(3451,file='chis_realse.dat')
 open(3461,file='chis_imagse.dat')
 open(3471,file='optical-constantsee.dat')
 open(3481,file='selfenergy-imag.dat')
 open(3481,file='selfenergy-real.dat')
 open(3381,file='selfenergy.dat')

!call calculate_v35(ksubibz,ndn,nib,eivlibz,eivcibz,nk,eival,eivec)
qe=0d0
qe(3)=0.0001
! ise_em =0
 rse_em =0
 reflectivity =0
! do i=1,mesh !wmsubsetsize3 
    !ii=i+wm_shift
!write(*,*)'i---------------',i
!       do b=4,n
      ! qe=0d0
!       s=b-3
       
!       wes=(om(i)/cnst )**2
 !      write(*,*)'wes=(om(ii)/cnst )**2',wes       
!       call function_self_tetra(1,qe,b,om(i),tm,nselfs,uselfs)      

!       ise_em(i,s) = aimag(nselfs+uselfs)
       
!       write(*,*)'i,ise,rse',i,s,om(i),b, ise_em(i,s)!, rsempi_em(i,s)
!       write(3481,2)i,s,om(i),b, ise_em(i,s) 
!     enddo
!enddo
  
allocate(imgse(mesh),xkk(mesh),omm(mesh))

do s = 1, n - 3
    do i = 1, mesh 
        imgse(i) = ise_em(i, s)/(2*pi)
        xkk(i) = om(i)**2
    enddo

    do i = 1,mesh  
        !ii = i + wm_shift
        omm(i) = om(i)**2
        call kk_im2re(omm(i), mesh, imgse, xkk, rekk)
        rse_em(i, s) = rekk
        write(3421,2) s,i,om(i),rse_em(i, s)

    enddo
enddo
deallocate(imgse,xkk,omm)


do i=1,mesh
   do b=4,n
    s=b-3
    write(3381,2)i,b,b-3,om(i),ise_em(i,s),rse_em(i, s)
    enddo
enddo



 do i=1,mesh !wmsubsetsize3 
    chi=0
    !   ii=i+wm_shift
write(*,*)'i---------------',i
    do al=1,3
    do be=1,3
    do b=4,5!n
       !qe=0d0
       s=b-3


       chi(al,be)=chi(al,be)- z(b,al)*z(b,be)  / &      !M(b,al) * M(b,be) ?????
 &                            (om(i)*om(i)-eival(b)**2 -2*eival(b)*(rse_em(i,s)-ci*(ise_em(i,s)/(2*pi))))
 ! &                            (om(i)*om(i)-eival(b)*cnst*cnst-2*sqrt(eival(b))*cnst*(rse_em(i,s)-ci*ise_em(i,s)))
    
    enddo
    enddo
    enddo
    chi=chi/volume_r*ee/eps0*1d10 * cnst*cnst
    if(i.eq.1) epsilon0=epsil+real(chi)
    chinormal=dot_product(normald,matmul((chi+epsil),normald))
    eps1 = real(chinormal)
    eps2 = aimag(chinormal)
    n_es = sqrt((eps1 + sqrt(eps1**2 + eps2**2))/2)
    k_es = sqrt((- eps1 + sqrt(eps1**2 + eps2**2))/2)
    reflectivity_nk = ((n_es - 1)**2 + k_es**2)/((n_es + 1)**2 + k_es**2) 
    !write(*,*)' chinormal',i,chinormal
    reflectivity= abs((sqrt(chinormal)-1)/(sqrt(chinormal)+1))**2
    


     
    
 !   write(*,*)'i,om,emissivityinsideeee',i,om(ii),real(chimpi(i,:,:)),reflectivitympi(i), 1-reflectivitympi(i)  
    write(*,*)'i,om,n,k,R',i,om(i),n_es, k_es, reflectivity_nk,1-reflectivity_nk
    write(*,*)'i,om,insideeeeee', i,om(i), aimag(chi)
    write(*,*)'i,om,emissivity',i,om(i),real(chi)!,reflectivity, emissivity 
    write(3451,2)om(i),real(chi),reflectivity, 1-reflectivity
    write(3461,2)om(i),aimag(chi)
    write(3471,2)i,om(i),n_es, k_es, reflectivity_nk,1-reflectivity_nk
 enddo


2 format(999(1x,g11.4))

 end subroutine dielectric_setetrampirzero
!===============

