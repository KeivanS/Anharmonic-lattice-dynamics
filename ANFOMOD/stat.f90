   program stat
! read a stream of data a(i) from stdin; 
! perform statistical analysis, see statfor.pdf
      implicit none
      integer i,m,k,l
      real(8), allocatable :: a(:)
      real(8) x
      character(66) filename

      read(*,*)filename
      open(10,file=filename,status='old')

      do i=1,10000000
       read(10,*,end=1) k,x
      enddo
    1 m=i-1
      rewind(10)
      allocate(a(m))
      do i=1,m
       read(10,*,end=1) k,x, a(i)
      enddo
      close(10)

      call corr(a,i)
      call blocking(a,i)
      call histo(a,i)
  
      end program stat
!============================================
      subroutine corr(a,n)
! estimated mean error and autocorrelation
      implicit none
      integer mc,i,k,n
      parameter(mc=200)
      real*8 a(n),ave,v,rkappa,c
! mean
      ave=0.d0
      do i=1,n
       ave=ave+a(i)
      enddo
      ave=ave/n
      write(*,'(''average   '',f20.5)')ave
! variance
      v=0.d0
      do i=1,n
       v=v+(a(i)-ave)**2
      enddo
      v=v/(n-1)
      write(*,'(''variance '',e20.5)')v
! autocorrelation (up to a mc steps) and correlation time
      open(2,file='corr.out')
      rkappa=1.d0
      do i=1,mc
       c=0.d0
       do k=1,n-i
        c=c+(a(k)-ave)*(a(k+i)-ave)
       enddo
       c=c/(n-i)/v
       write(2,*)i,c
       rkappa=rkappa+2*c
      enddo  
      close(2)
      rkappa=max(1.d0,rkappa)
      write(*,'(''t corr   '',f20.5)')rkappa
! effective number of data
      write(*,'(''n eff    '',f20.5)')n/rkappa
! error of mean
      write(*,'(''sigma    '',f20.5)')sqrt(v*rkappa/n)
      return
      end
!============================================
      subroutine blocking(a,n)
! blocking analysis
      implicit none
      integer i,k,l,n,nblk,large,isize,isize_step,minleft,nsizes
      parameter(minleft=20,nsizes=100)
      real*8 a(n),ab,error,average,average2,average3,average4
      open(2,file='blocking.out')
      large=n/minleft ! want at least minleft blocks left
      isize_step=max(1,large/nsizes) ! want at most ~nsizes block sizes
! loop on block size
      do isize=1,large,isize_step
! # blocks
       nblk=n/isize
       k=0
       average=0.d0
       average2=0.d0
       average3=0.d0
       average4=0.d0
       do i=1,nblk
! block averages
        ab=0.d0
        do l=1,isize
         k=k+1
         ab=ab+a(k)
        enddo
        ab=ab/isize
        average=average+ab
        average2=average2+ab*ab
        average3=average3+ab*ab*ab
        average4=average4+ab*ab*ab*ab
       enddo
       average=average/nblk
       average2=average2/nblk
       average3=average3/nblk
       average4=average4/nblk
! estimated error of mean at this block size
       error=sqrt((average2-average**2)/(nblk-1))
       write(2,5)isize,error,average,average2,average3,average4
      enddo
5     format(i5,9(1x,g11.5))
      close(2)
      return
      end
!============================================
      subroutine histo(a,n)
      implicit none
      integer i,j,n,m,nbin
      parameter(nbin=21)
      real*8 a(n),h(0:nbin+1),a_min,a_max,delta
      do j=0,nbin+1
       h(j)=0.d0
      enddo
! min and max
      a_min=a(1)
      a_max=a(1)
      do i=2,n
       a_min=min(a_min,a(i))
       a_max=max(a_max,a(i))
      enddo
! bin size
      delta=(a_max-a_min)/nbin
! histogram
      do i=1,n
       j=nint(0.5d0+(a(i)-a_min)/delta)
       h(j)=h(j)+1.d0
      enddo
      h(1)=h(1)+h(0)
      h(nbin)=h(nbin)+h(nbin+1)
! write
      open(2,file='histo.out')
      do j=1,nbin
       write(2,*)a_min+(j-0.5d0)*delta,h(j)/(n*delta)
      enddo
      close(2)
      return
      end
