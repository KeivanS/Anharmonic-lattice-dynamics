 module workfft
    integer n_previous,nfft
!   complex(8), allocatable:: cwk(:),fwk(:)
    real(8), allocatable:: work(:)
    contains
      subroutine allocfft(msh)
        integer msh
!       if(allocated(cwk)) deallocate(cwk,fwk,work)
!       allocate( cwk(msh),fwk(msh),work(0:msh-1) )
        allocate( work(0:msh-1) )
      end subroutine allocfft

      subroutine fftini( msh )
 !    use workfft
      implicit none
      integer msh
!k1   msh  = 2**(int(log(real(msh))/log(2d0) + 0.999d0))
      msh  = 2**(int(log(real(msh))/log(2d0) ))
      call allocfft(msh)
      n_previous = 0  !make sure sine/cosine table is recalculated
      return
      end subroutine fftini

 end module workfft
!------------------------------------
 module units
 integer,parameter :: ulog=40,uparam=10,utraj=50,ucur0=60
 end module units
!===========================================================
 module params
 real(8) :: pi,t_unit,w_unit,time
 integer mode_no,relax_time,no_of_blocks
 end module params
!===========================================================
 module correlation
 complex(8), allocatable :: z(:),zhat(:),zt(:),zthat(:),kap(:),zkappa(:)
 integer :: nz
 
   contains
    
    subroutine allocate_correl(n)
    integer n
    allocate( z(0:n-1),zhat(0:n-1),zt(0:n-1),zthat(0:n-1),kap(0:n-1),zkappa(0:n-1) )
    end subroutine allocate_correl
    
 end module correlation
!===========================================================
module current
implicit none

 integer atoms,natom,nsteps      
 real(8),allocatable :: currents(:,:), javg(:) ,currentk(:,:)

      contains

!---------
   subroutine read_trajectory
   use units
   use params
   use correlation
   implicit none
   integer :: i,j,k,nfile,ucurr
   real(8) :: mean,sd
   real(8) :: time0
   character char*3

   open(99,file='traj.out',status='old')
   read(99,*)char,natom
   write(ulog,*) 'natom=',natom
   close(99)

   nfile = natom/10  ! each file will contain 10 currents

   open(utraj,FILE='current.out',status='old')

   do i=1,nfile+1
     ucurr = ucur0+i
     write(char,'(i3.3)')i
     open(ucurr,FILE='current'//char//'.out',status='old')
   enddo

! read trajectory once, and determine Number of timesteps
    nsteps = 0
    fst_read :Do i=1,1000000000
         read(utraj,*,end=99) time
    enddo fst_read
99  rewind(utraj)
    write(*   ,*)' nstep reads=',i-1

    nsteps = i-1-relax_time 
!k1   nz = 2*nsteps

    write(ulog,*)' nsteps=',nsteps
    write(*   ,*)' nsteps=',nsteps
    atoms = 5
    allocate ( currents(nsteps,atoms) ,javg(nsteps) , currentk(nsteps,natom) )

! skip the first relax_time steps
    do j=1,relax_time
       read(utraj,*) time0    
       do i=1,nfile+1
          ucurr = ucur0 + i
          read(ucurr,*) time0    
       enddo
    enddo  

! read and store all current information
    do j=1,nsteps
       read(utraj,*) time,javg(j),(currents(j,i),i=1,atoms)
       do i=1,nfile
          ucurr = ucur0 + i
          read(ucurr,*) time,(currentk(j,k),k=1+(i-1)*10 ,i*10)
       enddo
       ucurr = ucur0+nfile+1
       read(ucurr,*) time,(currentk(j,k),k=1+nfile*10 ,natom)
    enddo  

! subtract their averages
    call stats(nsteps,javg,mean,sd)
    javg = javg - mean
    write(ulog,5)'JAVG: dim,mean,sd,mean/sd = ',nsteps,mean,sd,mean/sd
    do i=1,atoms
       call stats(nsteps,currents(:,i),mean,sd)
       currents(:,i)=currents(:,i)-mean
       write(ulog,5)'Real space currents on atom N/2**i: i, mean/sd = ',i-1,mean,sd,mean/sd
!      currents(:,i)=currents(:,i)-sum(currents(1:nsteps,i))/nsteps
    enddo
    do i=1,natom
       call stats(nsteps,currentk(:,i),mean,sd)
       currentk(:,i)=currentk(:,i)-mean
       write(ulog,5)'Real space currents in mode i: i, mean/sd = ',i,mean,sd,mean/sd
!      currentk(:,i)=currentk(:,i)-sum(currentk(1:nsteps,i))/nsteps
    enddo

   close(utraj)
   do i=1,nfile+1
     ucurr = ucur0+i
     close(ucurr)
   enddo

! write out in log.dat file the normal mode frequencies

   do i=1,natom
      write(ulog,4)i,2*sin(i*pi/2d0/(natom+1))
   enddo
   4 format('omega',i6,f12.5)
   5 format(a,i6,9(1x,f12.5))

   end subroutine read_trajectory

end module current
!===========================================================
 program spectrum
 use params
 use current
 use correlation
 use workfft
 use units
! use numerical_libraries
 implicit none
 integer i,ip
 real(8) sd,mean,correlation_time
 complex(8), allocatable:: kappa_tot(:)
!include 'DXMLDEF.FOR'
! external dfftrf
! record /DXML_D_FFT_STRUCTURE/ FFT_STRUCT

 open(ulog,file='log.dat',status='unknown')

!call test
!stop

 call read_params

 call read_trajectory

 t_unit = time/(relax_time+nsteps)
 w_unit = pi*2/nsteps/t_unit
!k1  w_unit = pi/2/nsteps/t_unit

 nfft=nsteps / no_of_blocks
 call fftini(nfft)
 no_of_blocks = nsteps/nfft
 w_unit = pi*2/nfft/t_unit
 write(ulog,*)'nfft        =',nfft
 write(ulog,*)'mode_no     =',mode_no
 write(ulog,*)'no_of_blocks=',no_of_blocks
 write(ulog,*)'tunit,wunit =',t_unit,w_unit
 allocate(kappa_tot(0:nfft-1))
 nz = nfft
 call allocate_correl(nz)
 print*,' files read and allocation done!'
 print*,' mode_no=',mode_no
 write(*,*)'nfft     =',nfft
 write(*,*)'no_of_blocks=',no_of_blocks
 do i=1,nfft
    write(ulog,4) i,currentk(i,mode_no) 
 enddo
4 format(i9,9(1x,g11.5))

 call get_correl_block(nsteps,currentk(:,mode_no))
 call write_out_r(0,nfft-1,real(zkappa),'rkappa.dat',w_unit)
 call write_out_r(0,nfft-1,aimag(zkappa),'ikappa.dat',w_unit)
 call write_out_c(0,nfft-1,zt,'correl.dat',t_unit)

 call get_correl_block(nsteps,javg)
 call write_out_r(0,nfft-1,real(zkappa),'rkappavg.dat',w_unit)
 call write_out_r(0,nfft-1,aimag(zkappa),'ikappavg.dat',w_unit)
 call write_out_c(0,nfft-1,zt,'coravg.dat',t_unit)

 call stats(nsteps,javg,mean,sd)
 correlation_time = 1 + 2d0/sd*sum(zt)
 write(ulog,*)" Statistics for Javg "
 write(ulog,*)" mean =",mean
 write(ulog,*)"  sd  =",sd
 write(ulog,*)" tau  =",correlation_time 

!kappa_tot = 0
!do i=1,natom 
!   call get_correl_block(nsteps,currentk(:,i))
!   kappa_tot = kappa_tot + zkappa
!enddo
!call write_out_r(0,nfft-1,real(kappa_tot),'rkappa_tot.dat',w_unit)
!call write_out_r(0,nfft-1,aimag(kappa_tot),'ikappa_tot.dat',w_unit)

! call get_correl1(nz,javg)
! call write_out(1,nz,xt,'zavg.dat',t_unit)
! call write_out(1,nz,xthat,'zavghat.dat',w_unit)
!
! call get_correl2
! call write_out(1,nz,xt,'z.dat',t_unit)
! call write_out(1,nz,xthat,'zhat.dat',w_unit)

 close(ulog)

 end program spectrum
!===========================================================
 subroutine read_params
 use units
 use params
! use correlation
 implicit none

 pi   = 4d0*atan(1d0)

 open(uparam,file='cond-params.dat',status='old')

 read(uparam,*) mode_no
 read(uparam,*) relax_time
 read(uparam,*) no_of_blocks

 close(uparam)

1 format(a,9(1x,i5))
3 format(a,3(1x,g10.4))

 end subroutine read_params
!===========================================================
 subroutine write_out_c(m,n,dos,fn,scale)
 implicit none
 integer i,n,m
 real(8) scale
 complex(8) dos(m:n)
 character*(*) fn

 open(30,file=fn,status='unknown')

 do i=m+(n-m)/2+1,n
      write(30,3) i,(i-n+m-1)*scale,dos(i)      
 enddo
 do i=m,m+(n-m)/2
      write(30,3) i,(i-0)*scale,dos(i)      
 enddo

 close(30)
3 format(i9,8(2x,g11.5))

 end subroutine write_out_c
!===========================================================
 subroutine write_out_r(m,n,dos,fn,scale)
 implicit none
 integer i,n,m
 real(8) dos(m:n),scale
 character*(*) fn

 open(30,file=fn,status='unknown')

 do i=m+(n-m)/2+1,n
      write(30,3) i,(i-n+m-1)*scale,dos(i)      
 enddo
 do i=m,m+(n-m)/2
      write(30,3) i,(i-0)*scale,dos(i)      
 enddo

 close(30)
3 format(i9,8(2x,g11.5))

 end subroutine write_out_r
!===========================================================
 subroutine write_out(m,n,dos,fn,scale)
 implicit none
 integer i,n,m
 real(8) dos(m:n),scale
 character*(*) fn

 open(30,file=fn,status='unknown')

 do i=m,n
      write(30,3) i,i*scale,dos(i)      
 enddo

 close(30)
3 format(i9,8(2x,g11.5))

 end subroutine write_out
!===========================================================
 subroutine get_correl(n,current,zkappa)
! this subroutine takes the current array of dimension n, and calculates
! the fourier transform of its causal autocorrelation function:
! zkappa(w) = sum_1^n  e^{iw_c t}  <j(t+x) j(x)> assuming j is periodic
! with w_c = w + i eta (eta=0+)
 use params
 use units
! use correlation
 use workfft
 implicit none
 integer w,t,i,ip,n
 real(8) current(n),sss,y(nfft)
 complex(8) zkappa(0:nfft-1),zcur(0:nfft-1),zcurhat(0:nfft-1)
    zcur(0:nfft-1)=current(1:nfft)
    y = 0

!   time_loop: do t=1,n
!       avg_loop: do i=1,n
!           ip=mod(i+t-1,n)+1
!           xt(t) = xt(t) + current(i)*current(ip) 
!       enddo avg_loop
!       xt(t) = xt(t)/n
!   enddo time_loop
!   call write_out(1,n,xt,'z3.dat',t_unit)

! Now take the FFT and write out the FFT of the current in file

 
 call cfftd1(zcur,nfft,1 ,zcurhat) ! direct fft (2 for inverse fft)
 open(8,file='curhat.dat')
 do i=0,nfft-1
  y(i+1) = abs(zcurhat(i))*abs(zcurhat(i))    ! y(w)=|J(w)|^2
  write(8,4) i,i*t_unit,zcur(i),i*w_unit,zcurhat(i),y(i+1)
 enddo
!call write_out_r(1,nfft,y,'zhat.dat',w_unit) ! this contains x(w)=|J(w)|^2
 close(8)

! get the correlator in time domain by using inverse FFT
!zcur(0:nfft-1) = cmplx(y(1:nfft),0d0)
!call cfftd1(zcur,nfft,2 ,zkappa)   ! zkappa=<J(0)J(t)> should be the same as xt above
!call write_out(0,nfft-1,real(zkappa),'z.dat',t_unit)

 zkappa = 0
 W_loop: do w=0,nfft-1 
     sss = 0
     do t=0,w-1
        sss = sss + abs(zcurhat(t))**2/(t-w)
     enddo
     do t=w+1,nfft-1
        sss = sss + abs(zcurhat(t))**2/(t-w)
     enddo
     zkappa(w) = cmplx( pi * abs(zcurhat(w))**2 , sss)
 enddo W_loop

4 format(i9,9(1x,g11.5))

 end subroutine get_correl
!===========================================================
 subroutine get_correl_block(n,current) ! zkappa, and zt are the outputs
! this subroutine takes the current array of dimension n, and calculates
! the fourier transform of its causal autocorrelation function:
! zkappa(w) = sum_1^n  e^{iw_c t}  <j(t+x) j(x)> assuming j is periodic
! with w_c = w + i eta (eta=0+)
! only the time average is taken by dividing the total time into no_of_blocks
! blocks and average the current over them. 
 use params
 use units
 use correlation
 use workfft
 implicit none
 integer nb,n
 real(8) current(n)
 complex(8) zcur(0:nfft-1),zcurhat(0:nfft-1)

! As a test do the calculation only for one block and write the conductivity
! in "kap1.dat", and the correlation function in "correl.dat"
   
! Now do the block averaging 
 zkappa = 0
 zt = 0
 block_loop: do nb=1,no_of_blocks

    zcur(0:nfft-1)=current(1+(nb-1)*nfft:nb*nfft)
    call kappa(nfft,zcur)   ! output is kap(w) and z(t)
    zkappa = zkappa + kap
    zt = zt + z 

 enddo block_loop 
 zkappa = zkappa/no_of_blocks  ! ensemble averaging performed
 zt = zt/no_of_blocks

4 format(i9,9(1x,g11.5))

 end subroutine get_correl_block
!===========================================================
 subroutine get_correl3(n,current)
 use params
 use units
 use correlation
 implicit none
 integer w,t,i,ip,n
 real(8) current(n),sss
 complex(8) kappa(n)
! kappa(w) = sum_1^n  e&iwt  <j(t+x) j(x)> assuming j is periodic

    time_loop: do t=1,n
          zt(t)=0
            avg_loop: do i=1,n
                ip=mod(i+t,n)+1
            zt(t) = zt(t) + current(i)*current(ip) 
            enddo avg_loop
            zt(t) = zt(t)/n
    enddo time_loop

    call write_out(1,nz,zt,'z3.dat',t_unit)

! Now take the FFT
    kappa = 0
    W_loop: do w=1,n !steps
      do t=1,n
        kappa(w) = kappa(w) + exp(cmplx(0,-(w*t*2*pi)/n))*zt(t)
      enddo
        kappa(w) = kappa(w)/n
    enddo W_loop
 zthat(1:n)=real(kappa(1:n)) !sqrt( real(kappa(1:n))**2+Imag(kappa(1:n))**2 )
 zthat(n+1:n+n)=Imag(kappa(1:n))

 write(ulog,*)' xnorm=',sum(kappa)

4 format(a,i9,5(1x,g11.5))

 end subroutine get_correl3
!===========================================================
 subroutine get_correl4(nn,current)
 use params
 use units
 use correlation
 implicit none
 integer w,t,i,ip,n,nn
 real(8) current(nn),sss
 complex(8) kappa(nn)
! kappa(w) = sum_1^n  e&iwt  <j(t+x) j(x)> assuming j is periodic
 n=nn !/2

    time_loop: do t=1,n
          zt(t)=0
            avg_loop: do i=1,n
                ip=i+t  !mod(i+t,n)+1
                zt(t) = zt(t) + current(i)*current(ip) 
            enddo avg_loop
            zt(t) = zt(t)/n
zt(t)=cos(pi*2*t/n)
    enddo time_loop
 t_unit =1
 call write_out(1,nz,zt,'z3.dat',t_unit)

! Now take the FFT
    kappa = 0
    W_loop: do w=1,n !steps
      do t=1,n
        kappa(w) = kappa(w) + zt(t)*cexp(cmplx(0,-(w*t*2*pi)/n))
      enddo
!       kappa(w) = kappa(w)/n
 print*,w,kappa(w)
    enddo W_loop

 zthat(1:n)=real(kappa(1:n)) !sqrt( real(kappa(1:n))**2+Imag(kappa(1:n))**2 )
!zthat(n+1:n+n)=imag(kappa(1:n))

 write(ulog,*)' xnorm=',sum(kappa)

4 format(a,i9,5(1x,g11.5))

!w_unit =1
! call write_out(1,nz,zthat,'z3hat.dat',w_unit)

 end subroutine get_correl4
!===========================================================
 subroutine test
 use workfft
 use params
 implicit none
 integer, parameter:: mesh=32
 complex(8) fin(0:mesh-1),fout(0:mesh-1)
 integer i,iopt
 real(8) twopi

 twopi = 8d0 * datan(1.d0)
 t_unit = 1
 w_unit = twopi/mesh/t_unit
 iopt = 1
 fin = 0
 do i=0,mesh-1
! fin(i) =   cos(15d0*i*twopi/mesh)
 enddo
 fin(3)=32
 call write_out(0,mesh-1,real(fin),'cur.dat',t_unit)
 call fftini(mesh)
 call cfftd1(fin,mesh,2 ,fout)
 call write_out_c(0,mesh-1,fout,'curhat.dat',w_unit)
!open(8,file='curhat.dat')
!do i=0,(mesh-1)/2
! write(8,4) i,i*t_unit,fin(i),i*w_unit,fout(i)
!enddo
!do i=(mesh-1)/2+1,mesh-1
! write(8,4) i,i*t_unit,fin(i),(i-mesh)*w_unit,fout(i)
!enddo
!close(8)

4 format(i6,9(1x,g11.5))

 end subroutine test
!------------------------------------
      subroutine cfftd1( cf, n, iopt,                    cc )
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
!*    f(0:n-1) - 1-dim complex array containing the original data      *
!*               to be fourier-transformed.                            *
!*               this array will not be conserved.                     *
!*      n      - integer showing the data length, .ge. 8.              *
!*               power of 2.                                           *
!*     iopt      1, for forward fourier transformation.                *
!*               2, for inverse fourier transformation.                *
!*                                                                     *
!*    (result)                                                         *
!*    c(0:n-1) - 1-dim complex array containing the fourier transform  *
!*               of f(0:n-1).                                          *
!*                                                                     *
!*    (work area)                                                      *
!*      work(n)  - real array of length n.                                *
!*               this area is used for a table of trigonometric        *
!*               functions.                                            *
!*                                                                     *
!*  remark                                                             *
!*      if you call cfftd1 repeatedly with the same data length n,     *
!*      be sure not to modify the array work(n), to avoid recalculation   *
!*      of the sin and cos functions.                                  *
!*                                                                     *
!*  subroutines required                                               *
!*      the following child subroutines are called from cfftd1.        *
!*        cfftc1, cfftc2, cfftc3, cfftc4.                              *
!*                                                                     *
!*  copyright                                                          *
!*      y. onodera, august, 1989, rewritten m. sluiter, feb 17, 2000   *
!*                                                                     *
!***********************************************************************
      use workfft
      implicit none
      integer n,n4,n2,np,alpha,beta,b2,i,j,k,iopt,key
      real(8) f(0:n/2-1,0:3),c(0:n/2-1,0:3),sn,cs,fn,wx,wy,step
      real(8) twopi
      complex(8) cf(n),cc(n)
      save sn,cs,n4,n2,fn,np,wx,wy,step

      twopi = 8d0 * atan(1.d0)
      if ( n .lt. 8 ) goto 999
      if ( n .eq. n_previous ) then
        if ( work(1) .eq. sn .and. work(n4-1) .eq. cs ) goto 10
      else
        fn = real(n)
        np = log( fn ) * 1.45
        if ( n .ne. 2**np ) goto 999
        n_previous = n
        n4 = n / 4
        n2 = n / 2
        wx = 1 / fn
        wy = - wx
        step = twopi * wx
      endif
! sin-cos table.
! sin( 2 pi j / n ) = work(j),
! cos( 2 pi j / n ) = work(n4-j).
      do j=0,n4-1
        work(j) = sin( j * step )
      enddo
      work(n4) = 1
      sn = work(1)
      cs = work(n4-1)

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
      beta  = n4
! outermost loop begins.
! alpha * beta = n/2
20    continue
! - - beta >= alpha,  earlier stages.
      if( key .lt. 0 ) then
        call cfftc1( alpha, beta, f, c, work )
      else
        call cfftc1( alpha, beta, c, f, work )
      endif
        key = - key
        alpha = alpha + alpha
        beta = beta / 2
        if ( beta .ge. alpha ) goto 20
! when alpha exceeds beta, the data is transposed.
        b2 = beta + beta
      if( key .lt. 0 ) then
        call cfftc2( alpha, b2, f, c )
        call cfftc2( alpha, b2, f(0,2), c(0,2) )
      else
        call cfftc2( alpha, b2, c, f )
        call cfftc2( alpha, b2, c(0,2), f(0,2) )
      endif
      key = - key
      if ( beta .ne. 1 ) then
! later stages, alpha > beta.
30      continue
        if ( key .lt. 0 ) then
          call cfftc3( alpha, beta, f, c, work, work(n4), work(n2) )
        else
          call cfftc3( alpha, beta, c, f, work, work(n4), work(n2) )
        endif
        key = - key
        alpha = alpha + alpha
        beta = beta / 2
        if( beta .ge. 2 ) goto 30
      endif
! final stage: alpha = n/2, beta = 1
      if ( key .ge. 0 ) f = c
      call cfftc4( n2, f, c, work )
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
      end subroutine cfftd1

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
!========================================================================
 subroutine kappa(nft,zcur)  ! kap and z are the outputs
! this subroutine takes the current array of dimension n, and calculates
! the fourier transform of its causal autocorrelation function:
! zkappa(w) = sum_1^n  e^{iw_c t}  <j(t+x) j(x)> assuming j is periodic
! with w_c = w + i eta (eta=0+)
! only the time average is taken by dividing the total time into no_of_blocks
! blocks and average the current over them. 
 use params
 use units
 use correlation
 use workfft
 implicit none
 integer w,t,i,nft
 real(8) sss
 complex(8) zcur(0:nft-1),zcurhat(0:nft-1)
 
 call cfftd1(zcur,nft,1 ,zcurhat) ! direct fft (2 for inverse fft)

 kap = 0
 do w=0,nft-1 
    sss = 0
    do t=0,w-1
       sss = sss + abs(zcurhat(t))**2/(t-w)
    enddo
    do t=w+1,nft-1
       sss = sss + abs(zcurhat(t))**2/(t-w)
    enddo
    kap(w) = cmplx( pi * abs(zcurhat(w))**2 , sss)
 enddo 
! call write_out_c(0,nft-1,zkappa,'kap1.dat',w_unit)
! now get the time correlation function
 do i=0,nft-1
    zhat(i) = zcurhat(i)*zcurhat(nft-i-1)
 enddo
 z = 0
 call cfftd1(zhat,nft,2 ,z) ! inverse fft (2 for inverse fft)
! call write_out_c(0,nft-1,x,'correl.dat',t_unit)

 end subroutine kappa
!========================================================================
 subroutine stats(m,x,mean,sd)
 use units
 implicit none
 integer m
 real(8):: mean,sd,errorbar
 real(8) :: x(m)

! calculate the statistical properties
 mean = sum(x)/m
 sd = sqrt(sum((x-mean)*(x-mean))/(m-1))
 errorbar = sd/sqrt(1d0*m)
! write(ulog,*) ' mean = ',mean
! write(ulog,*) '  sd  = ', sd 
! write(ulog,*) ' error on mean  = ',errorbar

4 format(2(i9),9(1x,g11.5))

 end subroutine stats

