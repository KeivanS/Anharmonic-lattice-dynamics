!=================================================================
 module workfft_t
    implicit none
    integer n_previous_t,nfft_t
    real(8), allocatable:: work_t(:)
    contains
      subroutine allocfft_t(msh)
        integer msh
        allocate( work_t(0:msh-1) )
      end subroutine allocfft_t

      subroutine fftini_t( msh )
      implicit none
      integer msh
      msh  = 2**(int(log(real(msh))/log(2d0) ))
      call allocfft_t(msh)
      n_previous_t = 0  !make sure sine/cosine table is recalculated
      return
      end subroutine fftini_t

 end module workfft_t
!=====================================================
 module data
 implicit none
 integer nfft,natom,nensemble
 real(8) pi,eta,dw,wmax,w,dt
 real(8), allocatable:: current(:,:,:)

   contains

   subroutine allocate_data(n,s)
   integer n,s
      allocate(current(3,n,s))
   end subroutine allocate_data

 end module data
!=====================================================
 program read
! reads the current time series from a file and computes the integral
! int_0^inf <j(t) j(0)> dt
! kappa = volume / (3 k_B T^2) int_0^inf <j(t).j(0)> dt
! so the output of this program needs to be multiplied by volume/k_B T^2 
 use workfft_t
 use data
 implicit none
 integer i,j,k

 pi = 4.d0*datan(1.d0)

 call read_traj

 eta = 3.

 call std_correlation
 print*,' std correlation done!'

! call pad_correlation
! print*,' pad correlation done!'

 call my_correlation
 print*,' fft correlation done!'

 end program read
!=====================================================
 subroutine read_traj
 use data
 implicit none
 integer i,j,k,s,al
 real(8) x,y,avg(3)
 character line*90

 open(20,file='anh_md.params',status='old')
   read(20,*) natom
   read(20,*) nensemble
   read(20,*) natom
   read(20,*) dt
   read(20,*) natom
   read(20,*) natom
   read(20,*) nfft
   read(20,*) natom
   read(20,*) natom,j,k
 close(20)
 dt = dt*k   ! this is the time difference between two successive data points

 print*,'nfft, nensemble read ', nfft, nensemble
 call allocate_data(nfft, nensemble)

 open(10,file='current.dat',status='old')
! open(10,file='cur2.dat',status='old')
 avg = 0
 do s=1,nensemble

    do i=1,nfft
       read(10,*)j,x,(current(al,i,s),al=1,3)
    enddo
! for each set of the ensemble subtract the average
    avg(1) = sum(current(1,:,s))/(1.*nfft)
    avg(2) = sum(current(2,:,s))/(1.*nfft)
    avg(3) = sum(current(3,:,s))/(1.*nfft)
    current(1,:,s) = current(1,:,s) - avg(1)
    current(2,:,s) = current(2,:,s) - avg(2)
    current(3,:,s) = current(3,:,s) - avg(3)

 enddo

! do i=1,nfft
!    write(15,*)i,current(i)
! enddo
 print*,'Trajectory read'
 close(10)

 end subroutine read_traj
!=====================================================
 subroutine std_correlation
 use data
 implicit none
 integer i,j,k,s
 complex(8) y,z,x
 real(8) , allocatable :: corr(:)

 allocate ( corr(0:nfft/2) )
 print*,'eta,pi=',eta,pi
 open(16,file='std_cor.dat')
 open(30,file='std_kappa.dat')

 corr = 0
 do i=0,nfft/2
    do s=1,nensemble
    do j=1,nfft/2
      corr(i) = corr(i) + current(1,i+j,s)*current(1,j,s) &
&                       + current(2,i+j,s)*current(2,j,s) &
&                       + current(3,i+j,s)*current(3,j,s) 
    enddo
    enddo
 enddo
 corr = corr/(nfft*nensemble/2.)/3*dt

 do i=0,nfft/2
    write(16,2)i*dt,corr(i),corr(i)*exp(-eta*i*2/float(nfft)),corr(i)*(1+cos(2*pi*i/float(nfft)))/2
 enddo

3 format(i9,1x,g10.4,99(3x,g10.4,1x,g10.4))
2 format(g10.4,99(3x,g10.4,1x,g10.4))

 dw = 2*pi/(nfft*dt)
 wmax = dw*nfft
 do w=0,wmax,dw
    x=0 ; y=0 ; z=0
    do j = 0,nfft/2
       x = x + corr(j) * exp(cmplx(0, j*w) )
       y = y + corr(j) * exp(cmplx(0, j*w) ) * exp(-eta*j*2/float(nfft))
       z = z + corr(j) * exp(cmplx(0, j*w) ) * (1+cos(2*pi*j/float(nfft)))/2.
    enddo
    write(30,2)w,x,y,z
 enddo

 deallocate ( corr )
 close(16)
 close(30)

 end subroutine std_correlation
!=====================================================
 subroutine pad_correlation
 use data
 implicit none
 integer i,j,k,ipj,count,s
 complex(8) y,z,x
 real(8) junk(3)
 real(8) , allocatable :: corr(:)
 complex(8) , allocatable :: kappa(:)

 allocate ( corr(0:nfft), kappa(0:nfft) )
 print*,'eta,pi=',eta,pi
 open(17,file='pad_cor.dat')
 open(31,file='pad_kappa.dat')

 corr = 0
 do i=0,nfft-1
    count = 0
    do s=1,nensemble
    do j=1,nfft
       ipj=i+j
       if (ipj .gt. nfft) then
          junk=0 ! =current(ipj)
       else
          count = count+1
          junk = current(:,ipj,s)
       endif
       corr(i) = corr(i) + junk(1)*current(1,j,s) + junk(2)*current(2,j,s) + &
&                          junk(3)*current(3,j,s)
    enddo
    enddo
    corr(i) = corr(i)/count/3.d0*dt
    write(17,3)i*dt,corr(i),corr(i)*exp(-eta*i/float(nfft)),corr(i)*(1+cos(pi*i/float(nfft)))/2.
 enddo

3 format(g10.4,99(3x,g10.4,1x,g10.4))

 dw = 2*pi/(nfft*dt)
 wmax = dw*nfft
 do w=0,wmax,dw
    x=0 ; y=0 ; z=0
    do j = 0,nfft
       x = x + corr(j) * exp(cmplx(0, j*w) )
       y = y + corr(j) * exp(cmplx(0, j*w) ) * exp(-eta*j/float(nfft))
       z = z + corr(j) * exp(cmplx(0, j*w) ) * (1+cos(pi*j/float(nfft)))/2.
    enddo
    write(31,3)w,x,y,z
 enddo

 deallocate ( corr, kappa )
 close(17)
 close(31)

 end subroutine pad_correlation
!=====================================================
 subroutine my_correlation
 use data
 use workfft_t
 implicit none
 integer i,j,k,ipj,count,s,al,be
 complex(8) y,z,x
 real(8) junk,wunit,tunit
 real(8) , allocatable :: corr(:,:,:)
 complex(8) , allocatable :: kappa(:,:,:)
 complex(8) ksum
! complex(8), allocatable:: jc(:),ftjc(:),corl(:),ftcorl(:)
 complex(8), allocatable:: zt(:),zw(:)

 nfft_t = nfft * 2   ! *2 needed for my_correlation2

 call fftini_t(nfft_t)

 allocate ( zt(0:2*nfft-1), zw(0:nfft-1) )
 allocate ( corr(3,3,0:2*nfft-1),kappa(3,3,0:nfft-1) )
 print*,'eta,pi=',eta,pi
 tunit = dt
 wunit = 2*pi/float(nfft) /tunit
 print*,'wunit,nfft=',wunit,nfft_t
 open(17,file='fft_cor.dat')
 open(31,file='fft_kappa.dat')

 corr = 0d0 ; kappa = cmplx(0d0,0d0)
 do al = 1,3
 do be = 1,3
 do s=1,nensemble

    call fft_correlation(nfft,current(al,:,s),current(be,:,s),zt,zw)
    corr(al,be,:)  = corr(al,be,:)  + real(zt)
    kappa(al,be,:) = kappa(al,be,:) + zw

 enddo
 enddo
 enddo
 corr  = corr /nensemble*dt
 kappa = kappa/nensemble*dt

 do j = 0,nfft-1
    ksum = (kappa(1,1,j)+kappa(2,2,j)+kappa(3,3,j))/3
    write(17,3) j,tunit*j,corr(:,:,j)
    write(31,4) j,wunit*j,ksum,kappa(:,:,j)
 enddo

    deallocate ( zt, zw, corr, kappa )
    deallocate ( work_t )
3 format(i9,1x,g10.4,1x,99(1x,g10.4))
4 format(i9,1x,g10.4,1x,88(1x,g10.4))
    close(17)
    close(31)

 end subroutine my_correlation
!========================================================================
 subroutine fft_correlation(n,cura,curb,zt,zw)
! takes the two time series cura and curb and computes, using padding and fft
! methods, the autocorrelation of them and its fourier transform (zt and zw)
 use data
 use workfft_t
 implicit none
 integer i,j,k,ipj,count,s,al,be,n
 complex(8) y,z,x
 real(8) junk,wunit
 real(8) cura(n),curb(n)
 complex(8) zt(0:2*n-1), zw(0:n-1) 
 complex(8), allocatable:: jc(:),ftjc(:) !,corl(:),ftcorl(:)

 allocate ( ftjc(0:2*n-1) )

 zt = cmplx(0,0) ; zw = cmplx(0,0)

! pad j=current with zeros
    allocate ( jc(0:2*n-1) )
    jc(0:n-1) = cmplx(cura(1:n),0)
    jc(n:2*n-1) = 0
    call cfftdt(jc,nfft_t,1,ftjc)

    jc(0:n-1) = cmplx(curb(1:n),0)
    jc(n:2*n-1) = 0
    call cfftdt(jc,nfft_t,1,zt)

    ftjc = ftjc*conjg(zt) !*sqrt(nfft_t*1.)
!   deallocate(jc)
!   allocate ( corl(0:2*n-1) )
! back from frequency to time
    call cfftdt(ftjc,nfft_t,2,jc)

! properly normalize (ftjc), and Tuckey-filter(corl)
    zt = cmplx(0,0) 
    do j = 0,n-1
       x= jc(j) /(n - j - 0)
       zt(j) = x*0.5*(1+cos(pi*(j)/(n-1.)))
!      ftjc(j) = x 
    enddo

!   allocate ( jc(0:n-1) )
    call cfftdt(zt(0:n-1),n,1,zw)
!   call cfftdt(ftjc(0:n-1),n,1,jc)
!   deallocate(jc)
    deallocate ( ftjc )
3 format(i9,1x,g10.4,9(3x,g10.4,1x,g10.4))

 end subroutine fft_correlation
!========================================================================
      subroutine cfftdt( cf, n, iopt, cc )
!***********************************************************************
!*                                                                     *
!*  purpose                                                            *
!*      complex fast fourier transformation. radix 2.                  *
!*      in double precision.                                           *
!*                                                                     *
!*  usage                                                              *
!*      call  cfftd1( f, n, iopt, c )                                  *
!*                                                                     *
!*  description of parameters                                          *
!*    (input)                                                          *
!*    f(1:n) - 1-dim complex array containing the original data        *
!*               to be fourier-transformed.                            *
!*               this array will not be conserved.                     *
!*      n      - integer showing the data length, .ge. 8.              *
!*               power of 2.                                           *
!*     iopt      1, for forward fourier transformation.                *
!*               2, for inverse fourier transformation.                *
!*                                                                     *
!*    (result)                                                         *
!*    c(1:n) - 1-dim complex array containing the fourier transform    *
!*               of f(0:n-1).                                          *
!*                                                                     *
!*    (work_t area)                                                    *
!*      work_t(n)  - real array of length n.                           *
!*               this area is used for a table of trigonometric        *
!*               functions.                                            *
!*                                                                     *
!*  remark                                                             *
!*      if you call cfftd1 repeatedly with the same data length n,     *
!*      be sure not to modify the array work_t(n), to avoid            *
!*      recalculation of the sin and cos functions.                    *
!*                                                                     *
!*  subroutines required                                               *
!*      the following child subroutines are called from cfftd1.        *
!*        cfftc1, cfftc2, cfftc3, cfftc4.                              *
!*                                                                     *
!*  copyright                                                          *
!*      y. onodera, august, 1989, rewritten m. sluiter, feb 17, 2000   *
!*                                                                     *
!***********************************************************************
      use workfft_t
      implicit none
      integer n,n4,n2,np,alpha,beta1,b2,i,j,k,iopt,key
      real(8) f(0:n/2-1,0:3),c(0:n/2-1,0:3),sn,cs,fn,wx,wy,step
      real(8) twopi
      complex(8) cf(n),cc(n)
      save sn,cs,n4,n2,fn,np,wx,wy,step

      twopi = 8d0 * atan(1.d0)
      if ( n .lt. 8 ) goto 999
      if ( n .eq. n_previous_t ) then
        if ( work_t(1) .eq. sn .and. work_t(n4-1) .eq. cs ) goto 10
      else
        fn = real(n)
        np = log( fn ) * 1.45
        if ( n .ne. 2**np ) goto 999
        n_previous_t = n
        n4 = n / 4
        n2 = n / 2
        wx = 1 / fn
        wy = - wx
        step = twopi * wx
      endif
! sin-cos table.
! sin( 2 pi j / n ) = work_t(j),
! cos( 2 pi j / n ) = work_t(n4-j).
      do j=0,n4-1
        work_t(j) = sin( j * step )
      enddo
      work_t(n4) = 1
      sn = work_t(1)
      cs = work_t(n4-1)

10    continue
! MS modification, map complex CF onto real F
      i=0
      do j=0,3
        do k=0,n2-2,2
          i        = i+1
          f(k,j)   = real( cf(i) )
          f(k+1,j) = aimag( cf(i) )
        enddo
      enddo

      do j=0, n4-1
        k = j + n4
! real part
        c(j,0) = f(2*j,0) + f(2*j,2)
        c(j,1) = f(2*j,0) - f(2*j,2)
        c(k,0) = f(2*j,1) + f(2*j,3)
        c(k,1) = f(2*j,1) - f(2*j,3)
! imaginary part
        c(j,2) = f(2*j+1,0) + f(2*j+1,2)
        c(j,3) = f(2*j+1,0) - f(2*j+1,2)
        c(k,2) = f(2*j+1,1) + f(2*j+1,3)
        c(k,3) = f(2*j+1,1) - f(2*j+1,3)
      enddo
      if( iopt .eq. 2 ) then
        c(:,0) = c(:,0) * wx
        c(:,1) = c(:,1) * wx
        c(:,2) = c(:,2) * wy
        c(:,3) = c(:,3) * wy
      endif
      key = 1
! key = -1, when the data is stored in the array f.
! key = +1, when the data is stored in the array c.
      alpha = 2
      beta1  = n4
! outermost loop begins.
! alpha * beta1 = n/2
20    continue
! - - beta1 >= alpha,  earlier stages.
      if( key .lt. 0 ) then
        call cfftc1( alpha, beta1, f, c, work_t )
      else
        call cfftc1( alpha, beta1, c, f, work_t )
      endif
        key = - key
        alpha = alpha + alpha
        beta1 = beta1 / 2
        if ( beta1 .ge. alpha ) goto 20
! when alpha exceeds beta1, the data is transposed.
        b2 = beta1 + beta1
      if( key .lt. 0 ) then
        call cfftc2( alpha, b2, f, c )
        call cfftc2( alpha, b2, f(0,2), c(0,2) )
      else
        call cfftc2( alpha, b2, c, f )
        call cfftc2( alpha, b2, c(0,2), f(0,2) )
      endif
      key = - key
      if ( beta1 .ne. 1 ) then
! later stages, alpha > beta1.
30      continue
        if ( key .lt. 0 ) then
          call cfftc3( alpha, beta1, f, c, work_t, work_t(n4), work_t(n2) )
        else
          call cfftc3( alpha, beta1, c, f, work_t, work_t(n4), work_t(n2) )
        endif
        key = - key
        alpha = alpha + alpha
        beta1 = beta1 / 2
        if( beta1 .ge. 2 ) goto 30
      endif
! final stage: alpha = n/2, beta1 = 1
      if ( key .ge. 0 ) f = c
      call cfftc4( n2, f, c, work_t )
      if( iopt .eq. 2 ) then
        do j=1, n2-1, 2
          c(j,:) = -c(j,:)
        enddo
      endif
! exit
! MS modification, map real C onto complex CC
      i=0
      do j=0,3
        do k=0,n2-2,2
          i        = i+1
          cc(i)    = cmplx( c(k,j),c(k+1,j) )
        enddo
      enddo
      return

999   write(*,'(a,i5)')'CFFTD1: error [n<8 or n><2**i], n=',n
      stop
      end subroutine cfftdt
!==============================================================
      subroutine cfftc1( n, m, f, g, s )
! child subroutine used in earlier ( n <= m ) stages of complex fft
! n           2**l
! m           N / ( 2 n )
! f(0:n-1)    input data
! g(0:n-1)    output data
! s(0:n/4)    table of sin( 2 pi j / n )
      implicit none
      integer n,m,j,k,n2
      real(8) f(0:2*m-1,0:n-1,2),g(0:m-1,0:n-1,0:3),s(0:m-1,0:*), &
     &   sn,cs,u,v
! for k=0:  cs=1, sn=0
      do j=0,m-1
        g(j,0,0) = f(j,0,1) + f(j+m,0,1)
        g(j,0,1) = f(j,0,1) - f(j+m,0,1)
        g(j,0,2) = f(j,0,2) + f(j+m,0,2)
        g(j,0,3) = f(j,0,2) - f(j+m,0,2)
      enddo
      if ( n .eq. 1 )  return
      n2 = n / 2
! for k=n/2: cs=0, sn=1
      do j=0,m-1
        g(j,n2,0) = f(j,n2,1) + f(j+m,n2,2)
        g(j,n2,1) = f(j,n2,1) - f(j+m,n2,2)
        g(j,n2,2) = f(j,n2,2) - f(j+m,n2,1)
        g(j,n2,3) = f(j,n2,2) + f(j+m,n2,1)
      enddo
      if ( n .eq. 2 )  return

      do k=1,n-1
        if ( k .ne. n2 ) then
          if ( k .lt. n2 ) then
            cs = s(0,n2-k)  ! = cos( 2 pi m k / n ) = cos( pi k / n )
            sn = s(0,k)     ! = sin( 2 pi m k / n ) = sin( pi k / n )
          else
            cs = - s(0,k-n2)
            sn = s(0,n-k)
          endif
          do j=0,m-1
            u = sn * f(j+m,k,2) + cs * f(j+m,k,1)
            v = cs * f(j+m,k,2) - sn * f(j+m,k,1)
            g(j,k,0) = f(j,k,1) + u
            g(j,k,1) = f(j,k,1) - u
            g(j,k,2) = f(j,k,2) + v
            g(j,k,3) = f(j,k,2) - v
          enddo
        endif
      enddo
      return
      end subroutine cfftc1
!==============================================================
      subroutine cfftc2( n, m, f, g )
! child subroutine for transposing matrices
      implicit none
      integer n,m,i,j
      real(8) f(m,n),g(n,m)

      if ( n .ge. m ) then
        do j=1,m
          do i=1,n
            g(i,j) = f(j,i)
          enddo
        enddo
      else
        do i=1,n
          do j=1,m
            g(i,j) = f(j,i)
          enddo
        enddo
      endif
      return
      end subroutine cfftc2
!==============================================================
      subroutine cfftc3( a, b, f, g, s, cs, sn )
! child subroutine used in later ( a > b ) stages of complex fft
! a           2**l
! b           n / ( 2 a )
! f(0:n-1)    input data
! g(0:n-1)    output data
! s(0:n/4)    table of sin( 2 pi j / n )
      implicit none
      integer a,b,j,k,a2
      real(8) f(0:a-1,0:b-1,0:3),g(0:a-1,0:1,0:b-1,2),s(0:b-1,0:*), &
     &   cs(0:a-1),sn(0:a-1),u,v

      a2 = a / 2
      do k=0,a2-1
        j = k + a2
        cs(k) = s(0,a2-k)
        sn(j) = cs(k)
        sn(k) = s(0,k)
        cs(j) = - sn(k)
      enddo
      do j=0,b-1
        do k=0,a-1
          u = sn(k) * f(k,j,3) + cs(k) * f(k,j,1)
          v = cs(k) * f(k,j,3) - sn(k) * f(k,j,1)
          g(k,0,j,1) = f(k,j,0) + u
          g(k,1,j,1) = f(k,j,0) - u
          g(k,0,j,2) = f(k,j,2) + v
          g(k,1,j,2) = f(k,j,2) - v
        enddo
      enddo
      return
      end subroutine cfftc3
!==============================================================
      subroutine cfftc4( n2, f, g, s )
! child subroutine used in the final stage ( a=n/2, b=1 )
! of complex fft.
      implicit none
      integer k, n2, n4
      real(8) f(0:n2-1,0:3),g(0:2*n2-1,2),s(0:*),u,v,w,x

      n4 = n2 / 2   != n / 4
      do k=0,n4-1
        u = s(n4-k) * f(k,1) + s(k) * f(k,3)
        v = s(n4-k) * f(k,3) - s(k) * f(k,1)
        g( 2*k, 1) = f(k,0) + u
        g( 2*k, 2) = f(k,0) - u
        g(2*k+1,1) = f(k,2) + v
        g(2*k+1,2) = f(k,2) - v
      enddo
      do k=n4,n2-1
        w = s(k-n4) * f(k,1) - s(n2-k) * f(k,3)
        x = s(k-n4) * f(k,3) + s(n2-k) * f(k,1)
        g( 2*k, 1) = f(k,0) - w
        g( 2*k, 2) = f(k,0) + w
        g(2*k+1,1) = f(k,2) - x
        g(2*k+1,2) = f(k,2) + x
      enddo

      return
      end subroutine cfftc4
