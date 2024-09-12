module cells
!! contains primitice and supercell translation vectors, as well as
!! the equilibrium position of their atoms
real(8) r01(3),r02(3),r03(3),r1(3),r2(3),r3(3),rsc1(3),rsc2(3),rsc3(3)
real(8) g01(3),g02(3),g03(3),gsc1(3),gsc2(3),gsc3(3)
real(8) prim_to_cart(3,3),cart_to_prim(3,3),prim_to_cart0(3,3),cart_to_prim0(3,3)
real(8), allocatable:: pos_conv_red(:,:),pos_prim(:,:),pos_sc(:,:),mass(:),mass_sc(:) !,kp(:,:)
integer, allocatable:: ipos(:,:)
real(8) vol0,vol_sc
character(2), allocatable :: atname(:),atname_sc(:)
integer natom0(10),ntype,natoms0,natom_sc,cell,nsnap,nrand,n_sc,iaux

end module cells
!========================
module ios
integer, parameter:: ulog=9,usnap=8,umode=7,uval=11
end module ios
!========================
module pairpot

    real(8) a4,d4,req,alfa
    contains

     function pot(r) result(V)
     implicit none
     real(8), intent(in):: r
     real(8) V

! morse potential V(r)=D(1-exp(-a(r-re)))^2
     d4=10 ; a4=req/1.0; !req=2.27 !2.251666
     V=d4*(1-exp(-a4*(r-req)))**2
! 1/r^4 potential produces eesonable vibrational modes and frequencies than morse
     V=1/r**4

     end function pot
    !==================================
     function d1v(r) result(v1)
     implicit none
     real(8), intent(in):: r
     real(8) v1,findif !fd_deriv

!     v1=2*a*D*(1-exp(-a*(r-req)))*exp(-a*(r-req))
      v1=findif(r,pot)

     end function d1v
    !!==================================
     function d2v(r) result(v2)
     implicit none
     real(8), intent(in):: r
     real(8) v2 ,findif

!     v2=2*a*a*D*(2*exp(-a*(r-req))-1)*exp(-a*(r-req))
      v2=findif(r,d1v)

     end function d2v

end module pairpot
!========================
program generate_snapshots
!! this program generates a surpercell and snapshots of it where
!! atoms have been randomly moved along (close to) normal modes by a given amplitude
!! given the information on the primitive lattice, supercell size and
!! displacements amplitude, translated in terms of a temperature
!! INPUT FILES: cell.inpi for the primitive cell format is similar to params.inp
!!              supercell.inp contains the supercell transl vect in terms of the primitiv
!!              snaps_temp.inp contains the temperature corresponding to that temperature
use cells
use ios
use pairpot
implicit none
integer i,i0,j,s,s2,k, mode,nmodes
real(8), allocatable:: phi(:),ksi(:),normal_modes(:,:),freq_thz(:),disp(:,:),vel(:,:),gaussian(:),dgaussian(:)
real(8) tempk,wavgcm,zero(3)
integer, allocatable:: seed(:)
character(3) csn

   open(ulog,file='log.dat')
   call read_primitive
   call read_supercell
   call read_snaps(wavgcm,tempk,cell,nsnap)
   write(ulog,*)'Nearest neighbor distance=',req

   nmodes=3*natom_sc
   zero = 0
   write(*,*)'there are ',nmodes,' modes '
   write(ulog,*)'there are ',nmodes,' modes '

   allocate(freq_thz(nmodes),normal_modes(nmodes,nmodes))
   call get_normal_modes(natom_sc,pos_sc,mass_sc,normal_modes,freq_thz,wavgcm)
   nrand = nmodes*nmodes  ! don't exclude acoustic (translational) modes

! if snapshots are displacements along normal modes
   allocate(ksi(nrand),phi(nrand),seed(nrand),disp(3,natom_sc),vel(3,natom_sc)) !random numbers for each mode
   allocate(gaussian(nrand),dgaussian(nrand)) !modes))
   write(ulog,*)'there are ',nrand,' random numbers generated '

   call random_seed ( )   ! initialize based on time and date
   call random_number(ksi)
   call random_number(phi)
   open(60,file='random.dat')
   write(60,*)'# random numbers ',nrand
   do i=1,nrand
      write(60,*) phi(i) , ksi(i)
   enddo
   close(60)
   phi=phi*6.28318531
   ksi=sqrt(-2*log(1-ksi))
   write(*,*) ' generating ',nrand,' random number of gaussian distribution'
   open(50,file='gauss_0_1.dat')
   do i=1,nrand
      gaussian(i)= cos(phi(i))*ksi(i)
     dgaussian(i)= sin(phi(i))*ksi(i)
      write(50,4) ' ',gaussian(i),dgaussian(i)
   enddo
   close(50)

   open(usnap,file='snapshots.xyz')
! writing first (unmoved) snapshot
   write(usnap,*) natom_sc
   write(usnap,*)'# Snapshot 0'
   s=0
   write(csn,'(i3.3)')s
!  open(70+s,file='snap_'//csn//'.xyz')
   open(70+s,file='poscar_'//csn)
   write(70+s,*)'# undistorted supercell '
   write(70+s,*) ' 1 '
   write(70+s,3) rsc1 
   write(70+s,3) rsc2 
   write(70+s,3) rsc3 
   write(70+s,*) (n_sc*natom0(j),j=1,ntype) !natom_sc
   write(70+s,*)' C '
   do i=1,natom_sc
      write(usnap,4) atname_sc(i),pos_sc(:,i)
      write(70+s ,4)' ',pos_sc(:,i),zero
   enddo

   open(60,file='sqdisp.dat')
   open(61,file='vel.dat')

         snaploop: do s=1,nsnap 

           write(csn,'(i3.3)')s
           open(70+s,file='poscar_'//csn)
           write(70+s,*)'# distorted supercell number ',s
           write(70+s,*) ' 1 '
           write(70+s,3) rsc1 
           write(70+s,3) rsc2 
           write(70+s,3) rsc3 
           write(70+s,*) (n_sc*natom0(j),j=1,ntype) !natom_sc
           write(70+s,*)' C '

           write(usnap,*)'# Snapshot ',s

           disp=0; vel=0
           do i=1,natom_sc
           do j=1,3
           do mode=4,nmodes
              disp(j,i)=disp(j,i)+normal_modes(3*(i-1)+j,mode)  &
        &      /freq_thz(mode)*sqrt(tempk/mass_sc(i))* gaussian(nmodes*(s-1)+mode)
              vel(j,i) = vel(j,i)+normal_modes(3*(i-1)+j,mode)  &
        &                     *sqrt(tempk/mass_sc(i))*dgaussian(nmodes*(s-1)+mode)
           enddo
           enddo
           enddo
           disp=1e10*sqrt(1.38e-23/1.66e-27)/1e12/6.28318531*disp  !=0.14528*disp in Ang
           vel =1e10*sqrt(1.38e-23/1.66e-27)/1e12/6.28318531*vel*100   !=14.528*disp in m/s


!          write(60,*)'# snap ',s
           do i=1,natom_sc
              do j=1,3
                 write(60,3) disp(j,i) !sqrt(disp(i,1)**2+disp(i,2)**2+disp(i,3)**2)
                 write(61,3) vel(j,i) !sqrt(disp(i,1)**2+disp(i,2)**2+disp(i,3)**2)
              enddo
              write(usnap,4) atname_sc(i),pos_sc(:,i)+disp(:,i),vel(:,i)
              write(70+s ,4) ' ' ,pos_sc(:,i)+disp(:,i),vel(:,i)
           enddo
           close(70+s)

         enddo snaploop

   close(60)
   close(61)

   deallocate(ksi,phi,gaussian,seed,disp,vel)

   close(usnap)
   close(ulog)

3 format(9(1x,f15.8))
4 format(a,2(2x,3(1x,f14.7)))


end program generate_snapshots
!===============================
subroutine read_primitive
use cells
use ios
use pairpot
implicit none
real(8) n1(3),n2(3),n3(3),latp(6),ms(10),scal,pi !, dcosd,dsind
integer i,j,k,ier

pi=4d0*datan(1d0)
!! first read conventional cell; format is similar to params.inp
open(10,file='cell.inp',status='old')
read(10,*)latp !=a,b,c,al,be,ga
read(10,*) n1,n2,n3
read(10,*) scal

 latp(1:3) = latp(1:3)*scal   ! convert angles from degree to radian
 latp(4:6)=latp(4:6)*pi/180
 r1(1)=latp(1)
 r2(1)=latp(2)*dcos(latp(6))
 r2(2)=latp(2)*dsin(latp(6))
 r3(1)=latp(3)*dcos(latp(5))
 r3(2)=latp(3)*(dcos(latp(4))-dcos(latp(6))*dcos(latp(5)))/dsin(latp(6))
 r3(3)=sqrt(latp(3)**2-r3(1)**2-r3(2)**2)

write(ulog,*) 'Conventional cell'
write(ulog,4) 'r1=',r1
write(ulog,4) 'r2=',r2
write(ulog,4) 'r3=',r3


!! generate primitive cell
 r01=n1(1)*r1+n1(2)*r2+n1(3)*r3
 r02=n2(1)*r1+n2(2)*r2+n2(3)*r3
 r03=n3(1)*r1+n3(2)*r2+n3(3)*r3

write(ulog,*) 'Primitive cell'
write(ulog,4) 'r01=',r01
write(ulog,4) 'r02=',r02
write(ulog,4) 'r03=',r03

 prim_to_cart0(:,1) = r01
 prim_to_cart0(:,2) = r02
 prim_to_cart0(:,3) = r03

 call xmatinv(3,prim_to_cart0,cart_to_prim0,ier)

 g01= cart_to_prim0(1,:)
 g02= cart_to_prim0(2,:)
 g03= cart_to_prim0(3,:)

write(ulog,*) 'Primitive cell reciprocal lattice'
write(ulog,4) 'g01=',g01
write(ulog,4) 'g02=',g02
write(ulog,4) 'g03=',g03

! read primitive atoms
read(10,*) ntype
read(10,*) (natom0(j),j=1,ntype)
read(10,*) (ms(j),j=1,ntype)

natoms0=0
do j=1,ntype
   natoms0=natoms0+natom0(j)
enddo
write(ulog,*)' # of atom types read      =',ntype
write(ulog,*)' Total # of atoms in the primitive cell=',natoms0

allocate(pos_conv_red(3,natoms0),pos_prim(3,natoms0),mass(natoms0),atname(natoms0))
k=0
do j=1,ntype
   do i=1,natom0(j)
      k=k+1
      read(10,*)atname(k), pos_conv_red(:,k)
      mass(k)=ms(j)
   enddo
enddo

! this pos is in units of conventional vectors; convert to cartesian
write(ulog,*)' Cartesian coordinates of the primitive atoms '
do k=1,natoms0
   pos_prim(:,k)=pos_conv_red(1,k)*r1+pos_conv_red(2,k)*r2+pos_conv_red(3,k)*r3
   write(ulog,3)pos_prim(:,k)
enddo
!req=1d9
!do k=1,natoms0
!   do j=1,natoms0
!      scal=sqrt(dot_product((pos_prim(k,:)-pos_prim(j,:)),(pos_prim(k,:)-pos_prim(j,:))))
!      if (scal.gt.1d-6 .and. scal .lt. req) req=scal
!   enddo
!enddo

3 format(3(1x,g14.8))
4 format(a,2x,3(1x,f14.8))

close(10)

end subroutine read_primitive
!===============================
subroutine read_supercell
! reads supercell translation vectors in terms of primitive vectors
! and generates coordinateds of atoms in the supercell
use cells
use ios
use pairpot
implicit none
integer n1(3),n2(3),n3(3),i,j,k,i1,i2,i3,mxshl,l,cnt,i0
real(8) r23(3),r_sc,rr(3),x,y,z,eps,rij

open(10,file='supercell.inp',status='old')
read(10,*) n1
read(10,*) n2
read(10,*) n3
close(10)

! rsc1=n1(1)*r01+n1(2)*r02+n1(3)*r03
! rsc2=n2(1)*r01+n2(2)*r02+n2(3)*r03
! rsc3=n3(1)*r01+n3(2)*r02+n3(3)*r03
 rsc1=n1(1)*r1+n1(2)*r2+n1(3)*r3
 rsc2=n2(1)*r1+n2(2)*r2+n2(3)*r3
 rsc3=n3(1)*r1+n3(2)*r2+n3(3)*r3

 prim_to_cart(:,1) = rsc1
 prim_to_cart(:,2) = rsc2
 prim_to_cart(:,3) = rsc3

 call xmatinv(3,prim_to_cart,cart_to_prim,i1 )

write(ulog,*) 'Super cell'
write(ulog,*) n1
write(ulog,*) n2
write(ulog,*) n3
write(ulog,4) 'rsc1=',rsc1
write(ulog,4) 'rsc2=',rsc2
write(ulog,4) 'rsc3=',rsc3


! generate all the atoms with r0i, translations of the primitive cell
! and select those that fall inside the SC

! first calculate the cell volumes
 call cross_product(r02,r03 ,r23)
 vol0  =abs(dot_product(r01 ,r23))
 if(vol0.eq.0) then
    write(ulog,*)'READ_SUPERCELL ERROR: vol0=0; change your primitive cell!'
    stop
 endif
 call cross_product(rsc2,rsc3,r23)
 vol_sc=abs(dot_product(rsc1,r23))
 if(vol_sc.eq.0) then
    write(ulog,*)'READ_SUPERCELL ERROR: vol_sc=0; change your supercell!'
    stop
 endif

 gsc1= cart_to_prim(1,:)
 gsc2= cart_to_prim(2,:)
 gsc3= cart_to_prim(3,:)

write(ulog,*) 'Super cell reciprocal vectors (no 2pi)'
write(ulog,4) 'gsc1=',gsc1
write(ulog,4) 'gsc2=',gsc2
write(ulog,4) 'gsc3=',gsc3

! call cross_product(rsc2,rsc3,b1)
! call cross_product(rsc3,rsc1,b2)
! call cross_product(rsc1,rsc2,b3)
! b1=b1/vol_sc
! b2=b2/vol_sc
! b3=b3/vol_sc

 r_sc = vol_sc/vol0
 write(ulog,*) 'Super cell to primitive volume ratio=',r_sc
 n_sc=nint(r_sc)
 natom_sc=n_sc*natoms0
 write(ulog,*) 'n_sc, natom_sc=',n_sc,natom_sc
 allocate(pos_sc(3,natom_sc),mass_sc(natom_sc),atname_sc(natom_sc),ipos(natoms0,n_sc))
 open(22,file='SC.xyz',status='unknown')
 write(22,*) natom_sc
 write(22,*)'junk'
 mxshl=nint(sqrt(max(dot_product(n1,n1),dot_product(n2,n2),dot_product(n3,n3))+1d0))+3 ! to be safe!
 write(ulog,*)'maxshell=',mxshl

eps=1d-5
cnt=0  ! counts the atoms in the supercell
i =0   ! counts the atoms in the primitive cell
do j=1,ntype
   do i0=1,natom0(j)
      i=i+1

k=0 ! counts the primitive cells
do l=0,mxshl
do i1=-l,l
do i2=-l,l
do i3=-l,l
! walk only on the facets of the cube of length 2n+1
   if(iabs(i1).ne.l.and.iabs(i2).ne.l.and.iabs(i3).ne.l)cycle

   rr=i1*r01+i2*r02+i3*r03 + pos_prim(:,i)
!   call bring_to_cell(cart_to_prim,prim_to_cart,rr,r23)
! check to see if rr is inside the supercell
   x=dot_product(rr,gsc1)
   if (x.lt.-eps .or. x.ge.1-eps) cycle
   y=dot_product(rr,gsc2)
   if (y.lt.-eps .or. y.ge.1-eps) cycle
   z=dot_product(rr,gsc3)
   if (z.lt.-eps .or. z.ge.1-eps) cycle

   k=k+1  ! vector rr is in the supercell
   cnt=cnt+1
   ipos(i,k)=cnt
   pos_sc(:,cnt)=rr(:)
   mass_sc(cnt)=mass(i)
   atname_sc(cnt)=atname(i)
   write(22,4) atname_sc(cnt),pos_sc(:,cnt) !,x,y,z,rr

enddo
enddo
enddo
enddo

   enddo
enddo


if(k.ne.n_sc) then
  write(ulog,*)'ERROR in supercell generation, n_sc,k=',n_sc,k
  stop
endif

! find the shortest bond length and store in req
req=1000000
do i=1,natoms0
    do j=1,natom_sc
        rr=pos_sc(:,j)-pos_prim(:,i)
        rij=sqrt(dot_product(rr,rr))
        if (rij.lt.req .and. rij.gt.1d-3) req=rij
    end do
end do
write(ulog,*)'shortest bond length=',req
req=req*1.01
write(ulog,*)'req=',req

4 format(a,2x,3(1x,f14.8),9(1x,g10.4))

close(22)

end subroutine read_supercell
!===============================
subroutine read_snaps(wavgcm,tempk,cell,nsnap)
!use pairpot, only : alfa
implicit none
real(8) tempk,wavgcm
integer cell,nsnap

open(10,file='snaps.inp',status='old')
read(10,*) wavgcm  ! average freq in 1/cm
read(10,*) tempk   ! temperature in K for canonical sampling
read(10,*) nsnap   ! number of desired snapshots, typically 50 should be good
close(10)

end subroutine read_snaps
!===============================
subroutine cross_product(a,b,c)
implicit none
real(8), intent(in) :: a(3),b(3)
real(8) c(3)

c(1)=a(2)*b(3)-a(3)*b(2)
c(2)=a(3)*b(1)-a(1)*b(3)
c(3)=a(1)*b(2)-a(2)*b(1)

end subroutine cross_product
!===============================
subroutine get_normal_modes(natom,ps,mas,modes,freq,wavgcm)
!! calculates the nromal modes of a crystal defined by atoms ps(natom,3))
!! and rescales the frequencies to get an average-squared equal to wavgcm^2
use cells
use ios
use pairpot
implicit none
integer, intent(in ):: natom
real(8), intent(in ):: ps(3,natom),mas(natom),wavgcm
real(8), intent(out):: freq(3*natom),modes(3*natom,3*natom)
real(8) fc(3*natom,3*natom),eival(3*natom)
! complex(8) eivec(3*natom,3*natom),cfc(3*natom,3*natom)
real(8) x,y,z,rij,rr(3),rp(3),rcut,cnst,sumw2,wscale,v1,v2
integer i,j,al,be,del,it_num,rot_num,it_max

!! construct the force constant matrix defined by rij_a*rij_be/rij^6 for in eV/Ang^2
!! only a couple of neighbors defined by rcut
!! the strength of the FCs are fixed so that the average frequency squared is wavgcm^2 in 1/cm

 it_max=100
 cnst=521.11  ! to convert sqrt(eV/Ang^2) to 1/cm
 rcut=req*1.1 !2.7 !0.9 * vol0**0.333
 write(ulog,*)'RCUT=',rcut

fc=0
do i=1,natom
do j=i+1,natom

   rr=ps(:,i)-ps(:,j)

! impose periodic BCs
   x=dot_product(gsc1,rr)
   y=dot_product(gsc2,rr)
   z=dot_product(gsc3,rr)
   x=x-nint(x) ! makes -0.5<x<0.5
   y=y-nint(y) ! makes -0.5<x<0.5
   z=z-nint(z) ! makes -0.5<x<0.5
   rp=x*rsc1+y*rsc2+z*rsc3 ! this has the shortest distance between i&j
   rij=sqrt(dot_product(rp,rp))

   v2=(d2v(rij)-d1v(rij)/rij)
   v1=(d1v(rij)/rij)

   if(rij.lt.rcut) then
     write(*,19)i,j,rij,v1,v2
     do al=1,3
     do be=1,3
        fc(3*(i-1)+al,3*(j-1)+be)=-rp(al)*rp(be)/rij/rij*v2 + del(al,be)*v1
        fc(3*(j-1)+be,3*(i-1)+al)=fc(3*(i-1)+al,3*(j-1)+be)
     enddo
     write(*,9) fc(3*(i-1)+al,3*(j-1)+1),fc(3*(i-1)+al,3*(j-1)+2),fc(3*(i-1)+al,3*(j-1)+3)
     enddo
   endif

enddo
enddo

! apply ASR
    do i=1,natom
    do al=1,3
    do be=1,3
       do j=1,natom
          if(i.eq.j) cycle
          fc(3*(i-1)+al,3*(i-1)+be)=fc(3*(i-1)+al,3*(i-1)+be) - fc(3*(i-1)+al,3*(j-1)+be)
       enddo
    enddo
    enddo
    enddo

    open(90,file='fcs.dat')
    do i=1,3*natom
       write(90,9)fc(i,:)
    enddo
    close(90)

! apply mass normalization
    do i=1,natom
    do j=1,natom
       do al=1,3
       do be=1,3
          fc(3*(i-1)+al,3*(j-1)+be)=fc(3*(i-1)+al,3*(j-1)+be)  /sqrt(mas(i)*mas(j))
       enddo
       enddo
    enddo
    enddo

! normalize trace or sum of eivals to 1
!   sumw2=0
!   do i=1,3*natom
!      sumw2=sumw2 + fc(i,i) !3*(i-1)+be,3*(i-1)+be)
!   enddo
!   sumw2=sumw2/(3*natom)
!   fc=fc/sumw2

! now diagonalize to get normal modes
! cfc=cmplx(fc,0)
! call EIGCH(cfc,3*natom,3*natom,3*natom,3*natom,1d-8, eival, eivec, IER)
 call jacobi_eigenvalue ( 3*natom, fc, it_max, modes, eival, it_num, rot_num )

 open(umode,file='modes.dat')
 open(uval,file='freqs.dat')

 if (it_num+1.eq.it_max) then !ier.ne.0) then
    write(ulog,*)'DIAGONALIZATION: jacobi did not converge, increase it_max'
    write( *  ,*)'DIAGONALIZATION: jacobi did not converge, increase it_max'
    stop
 else
    write(umode,*)"# UNSCALED EIGENVALUES & EIGENVECTORS OF THE FC MATRIX after ",rot_num,it_num
    sumw2=0
    do i=1,3*natom
       write(umode,*)i,eival(i),sqrt(abs(eival(i)))*cnst/33
       sumw2 = sumw2+eival(i)
       do j=1,natom
          write(umode,9)(modes(3*(j-1)+al,i),al=1,3)
       enddo
    enddo
 endif
 write(ulog,*)'trace eivals/3N (eV/Ang^2) is ',sumw2/3/natom

! normalize so that average eival is Wavg^2
 eival=eival/sumw2*(3*natom)* (wavgcm/cnst)**2
 write(uval,*) '# i,eival(i),omega(i) [THz] after renormalizing by [THz,1/cm] ',wavgcm/33,wavgcm
 do i=1,3*natom
    write(uval,*)i,eival(i),sqrt(abs(eival(i)))*cnst/33
 enddo

! wscale=wavg/(sqrt(sumw2/(3*natom-3))*cnst)
! write(ulog,*)'average frequency in cm^-1 /wavg is ',sqrt(sumw2/(3*natom-3))*cnst / wavg

 write(ulog,*)'Three lowest eigenvalues are'
 write(ulog,*)"1 ",eival(1)
 write(ulog,*)"2 ",eival(2)
 write(ulog,*)"3 ",eival(3)
! write(ulog,*)"1 ",eival(3*natom-0)
! write(ulog,*)"2 ",eival(3*natom-1)
! write(ulog,*)"3 ",eival(3*natom-2)

 write(ulog,*)'Frequency of all other modes in THz'
 do i=4,3*natom
   if(eival(i).gt.0) then
      freq(i)=sqrt(eival(i))*cnst/33 ! to convert to 1/cm and then to THz
      write(ulog,*)i,freq(i)
!      modes(:,i)=eivec(:,i)
   else
      write(ulog,*)"NEGATIVE EIGENVALUE!!!!!!!!!!!"
      write(ulog,*)i,eival(i)
      stop
   endif
 enddo

 close(umode)
 close(uval)

9 format(99(1x,f9.3))
19 format(2i5,99(1x,g9.3))


end subroutine get_normal_modes
!==================================
 function del(i,j) result(de)
 implicit none
 integer, intent(in):: i,j
 integer de
 if (i.eq.j) then
    de=1
 else
    de=0
 endif
 end function del
!==================================
subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )

!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gene Golub, Charles VanLoan,
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996,
!    ISBN: 0-8018-4513-X,
!    LC: QA188.G65.
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    real ( kind = rk ) A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    integer IT_MAX, the maximum number of iterations.
!
!  Output:
!
!    real ( kind = rk ) V(N,N), the matrix of eigenvectors.
!
!    real ( kind = rk ) D(N), the eigenvalues, in descending order.
!
!    integer IT_NUM, the total number of iterations.
!
!    integer ROT_NUM, the total number of rotations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) bw(n)
  real ( kind = rk ) c
  real ( kind = rk ) d(n)
  real ( kind = rk ) g
  real ( kind = rk ) gapq
  real ( kind = rk ) h
  integer i
  integer it_max
  integer it_num
  integer j
  integer k
  integer l
  integer m
  integer p
  integer q
  integer rot_num
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) tau
  real ( kind = rk ) term
  real ( kind = rk ) termp
  real ( kind = rk ) termq
  real ( kind = rk ) theta
  real ( kind = rk ) thresh
  real ( kind = rk ) v(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) zw(n)

  do j = 1, n
    do i = 1, n
      v(i,j) = 0.0D+00
    end do
    v(j,j) = 1.0D+00
  end do

  do i = 1, n
    d(i) = a(i,i)
  end do

  bw(1:n) = d(1:n)
  zw(1:n) = 0.0D+00
  it_num = 0
  rot_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
    thresh = 0.0D+00
    do j = 1, n
      do i = 1, j - 1
        thresh = thresh + a(i,j) ** 2
      end do
    end do

    thresh = sqrt ( thresh ) / real ( 4 * n, kind = rk )

    if ( thresh == 0.0D+00 ) then
      exit
    end if

    do p = 1, n
      do q = p + 1, n

        gapq = 10.0D+00 * abs ( a(p,q) )
        termp = gapq + abs ( d(p) )
        termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
        if ( 4 < it_num .and. &
             termp == abs ( d(p) ) .and. &
             termq == abs ( d(q) ) ) then

          a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
        else if ( thresh <= abs ( a(p,q) ) ) then

          h = d(q) - d(p)
          term = abs ( h ) + gapq

          if ( term == abs ( h ) ) then
            t = a(p,q) / h
          else
            theta = 0.5D+00 * h / a(p,q)
            t = 1.0D+00 / ( abs ( theta ) + sqrt ( 1.0D+00 + theta * theta ) )
            if ( theta < 0.0D+00 ) then
              t = - t
            end if
          end if

          c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
          s = t * c
          tau = s / ( 1.0D+00 + c )
          h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
          zw(p) = zw(p) - h
          zw(q) = zw(q) + h
          d(p) = d(p) - h
          d(q) = d(q) + h

          a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
          do j = 1, p - 1
            g = a(j,p)
            h = a(j,q)
            a(j,p) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = p + 1, q - 1
            g = a(p,j)
            h = a(j,q)
            a(p,j) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = q + 1, n
            g = a(p,j)
            h = a(q,j)
            a(p,j) = g - s * ( h + g * tau )
            a(q,j) = h + s * ( g - h * tau )
          end do
!
!  Accumulate information in the eigenvector matrix.
!
          do j = 1, n
            g = v(j,p)
            h = v(j,q)
            v(j,p) = g - s * ( h + g * tau )
            v(j,q) = h + s * ( g - h * tau )
          end do

          rot_num = rot_num + 1

        end if

      end do
    end do

    bw(1:n) = bw(1:n) + zw(1:n)
    d(1:n) = bw(1:n)
    zw(1:n) = 0.0D+00

  end do
!
!  Restore upper triangle of input matrix.
!
  do j = 1, n
    do i = 1, j - 1
      a(i,j) = a(j,i)
    end do
  end do
!
!  Ascending sort the eigenvalues and eigenvectors.
!
  do k = 1, n - 1

    m = k

    do l = k + 1, n
      if ( d(l) < d(m) ) then
        m = l
      end if
    end do

    if ( m /= k ) then

      t    = d(m)
      d(m) = d(k)
      d(k) = t

      w(1:n)   = v(1:n,m)
      v(1:n,m) = v(1:n,k)
      v(1:n,k) = w(1:n)

    end if

  end do


end subroutine jacobi_eigenvalue
!------------------------------------------------------------------------------
      subroutine xmatinv(n,xmatin,xmatout,ier)
      implicit none

!! invert a n by n matrix
      integer, intent(in) :: n
      integer, intent(out) :: ier
      real(8), intent(in) :: xmatin(n,n)
      real(8), intent(out) :: xmatout(n,n)
      real(8) x,buffer(n,n)
      integer indx(3),i,j

! clear error flag
      ier=0
      do i=1,n
        do j=1,n
          xmatout(i,j)=0
          buffer(i,j)=xmatin(i,j)
        enddo
        xmatout(i,i)=1
      enddo
! decomposition
      call ludcmp(buffer,n,n,indx,x)
! singular matrix
      if(x.eq.0.0d0)then
        ier=1
        return
      endif
! inverse matrix
      do j=1,n
        call lubksb(buffer,n,n,indx,xmatout(1,j))
      enddo

      end subroutine xmatinv
!-----------------------------------------------------------
! The following routines are from Numerical Recipes
      subroutine ludcmp(a,n,np,indx,d)
      implicit none
      integer nmax,np,n
      real(8) tiny
      parameter (nmax=3,tiny=1.0d-20)
      real(8) a(np,np),vv(nmax),d,aamax,dum,sum
      integer indx(n),i,j,k,imax,ncmp
      d=1
      do i=1,n
        aamax=0
        do j=1,n
          if(dabs(a(i,j)).gt.aamax)aamax=dabs(a(i,j))
        enddo
        if(ncmp(aamax).eq.0)then
! singular matrix
          d=0
          return
        endif
        vv(i)=1/aamax
      enddo
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if(dum.ge.aamax)then
            imax=i
            aamax=dum
          endif
        enddo
        if(j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.0d0)a(j,j)=tiny
        if(j.ne.n)then
          dum=1/a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
      end subroutine ludcmp
!---------------------------------------------------
      subroutine lubksb(a,n,np,indx,b)
      implicit none
      integer n,np
      real(8) a(np,np),b(n),sum
      integer indx(n),ii,i,j,ll
      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if(ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if(sum.ne.0.0d0)then
          ii=i
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do j=i+1,n
            sum=sum-a(i,j)*b(j)
          enddo
        endif
        b(i)=sum/a(i,i)
      enddo
      end subroutine lubksb
!----------------------------------------------------------------------------
      subroutine bring_to_cell(cart_to_prim,prim_to_cart,v1,v2)
!! bring a point into the unit cell at the origin

      implicit none
      real(8) , intent(in):: cart_to_prim(3,3),prim_to_cart(3,3),v1(3)
      real(8) , intent(out):: v2(3)
      real(8) buff(3)
      integer i,ncmp
! change coordinates of point to linear combination of basis vectors of the
! primitive lattice
      buff=matmul(cart_to_prim,v1)
! in the unit cell at the origin, the coefficient must be greater than or
! equal to zero and less than one.
      do i=1,3
        do while(buff(i).gt.1.0d0.or.ncmp(buff(i)-1).eq.0)
          buff(i)=buff(i)-1
        enddo
        do while(buff(i).lt.0.0d0.and.ncmp(buff(i)).ne.0)
          buff(i)=buff(i)+1
        enddo
      enddo
! return to cartesian coordinates
      v2= matmul(prim_to_cart,buff)

      end subroutine bring_to_cell
!--------------------------------------------------------------------------
      function ncmp(x)
!! COMPARE X WITH ZERO;	NCMP=0 IF X IS CLOSE ENOUGH TO ZERO; NCMP=1 OTHERWISE
      implicit none
      integer ncmp
      real(8) x,delta
      data delta/1.e-6/
      ncmp=0
      if(abs(x).gt.delta)ncmp=1

      end function ncmp
!==================================
 function fd_deriv(x,func) result(deriv)
 implicit none
 real(8), intent(in):: x
 real(8) h,deriv ,func
 external func
 h=1d-6
 deriv=(func(x+h)-func(x-h))/(2*h)
 end function fd_deriv
!==================================
  function findif(x,func) result(deriv)
!! derivative by Finite differentiation
  implicit none
  real(8), intent(in):: x
!  external func
  interface
      function func(x) result(y)
         real(8), intent(in) :: x
         real(8) y
      end function func
  end interface

 real(8) h,deriv
 h=1d-6
 deriv=(func(x+h)-func(x-h))/(2*h)
 end function findif
!==================================
