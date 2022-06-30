!===========================================================
 subroutine svd_set(m3,n,a,b,x,sig,svdcut,error,ermax,sigma,fnsvd)
! use mappings
 implicit none
 integer i,j,m3,n,k,uio
 real(8), intent(in) :: a(m3,n),b(m3)
! real(8) u(m3,n),v(n,n),w(n)
 real(8), allocatable :: u(:,:),v(:,:),w(:)
 real(8), intent(out) :: error,ermax,x(n),sig(n)
 real(8) prod,wmax,wmin,svdcut,junk,sigma,num,denom
 character(*) fnsvd
! m=dim_a ; n=ngr
!
! do i=1,m
! do j=1,n
!    call random_number( prod )
!    a(i,j) = 2*prod - 1
! enddo
!    call random_number( prod )
!    b(i) = 2*prod - 1
! enddo
 uio = 345

! open(uio,file='svd-results.dat',status='unknown')
 open(uio,file=fnsvd,status='unknown')

 write(uio,*)' Matrix A and vector b are '
! do i=1,m3
!   write(uio,9)i,(a(i,j),j=1,n),b(i)
! enddo
!9 format(i8,1x,99(1x,f7.3))
8 format(i8,1x,i8,3x,f25.15,3x,f25.15)
7 format(i8,1x,99(1x,g16.9))

 write(uio,*)' A and b written, now doing the SVD '

 allocate( u(m3,n),v(n,n),w(n) )

 u = a
 call svdcmp(u,m3,n,m3,n,w,v)

 write(uio,*)' product of U^T*U is: (no lines means it is identity) '
  do k=1,n
  do j=1,n
     prod = 0
     do i=1,m3
        prod = prod+u(i,j)*u(i,k)
     enddo
! if it is not 0 or 1 write...
     if( abs(abs(prod-0.5)-0.5) .gt. 1d-12) write(uio,8)k,j,prod
  enddo
  enddo

!write(uio,*)' product of V*V^T is: (no lines means it is identity) '
  do k=1,n
  do j=1,n
     prod = 0
     do i=1,n
        prod = prod+v(j,i)*v(k,i)
     enddo
! if it is not 0 or 1 write...
     if( abs(abs(prod-0.5)-0.5) .gt. 1d-12) write(uio,8)k,j,prod
  enddo
  enddo
 write(uio,*)' W (better not be too small) is:  '
   do k=1,n
      write(uio,*)'k,w(k)=',k,w(k)
   enddo
 write(uio,*)' Its condition number w_max/w_min should not be too large ',maxval(w)/minval(w)

!write(uio,*)' SVD product of U*W*V^T, a is: '
  do i=1,m3
  do j=1,n
     prod = 0
     do k=1,n
        prod = prod+v(j,k)*u(i,k)*w(k)
     enddo
!    write(uio,8)i,j,prod,a(i,j)
     if( abs(prod-a(i,j)) .gt. 1d-10) write(uio,8)i,j,prod,a(i,j)
  enddo
  enddo

  wmax = maxval(abs(w))
  wmin = svdcut*wmax
   write(uio,*)' Maxval of w is=',wmax
   write(uio,*)' w larger than ',wmin,' was used in the inversion'
! this is to eliminate the v vectors with small w
  do j=1,n
     if(abs(w(j)).lt.wmin) w(j)=0d0
  enddo
  call svbksb(u,w,v,m3,n,m3,n,b,x)

  write(uio,*)' results of SVD solution, variance, error are: x,sigma,Dx'
  do j=1,n
     sig(j) = 0
     do k=1,n
        if (w(k).gt.wmin) then
           sig(j) = sig(j) + v(j,k)*v(j,k) / w(k)/w(k)
        endif
     enddo
     sig(j) = sqrt(sig(j))
      write(uio,7) j,x(j),sig(j),sig(j)/sqrt(m3*1.)
  enddo

 deallocate(u,v,w)

!write(uio,*)' residual of SVD solution is: Ax,b,Ax-b,(Ax-b)/Ax'
  error = 0; ermax = 0; num=0; denom=0
  do i=1,m3
     prod = 0
     do j=1,n
        prod = prod + a(i,j)*x(j)
     enddo
     junk = abs(prod-b(i))
     num=num     + junk*junk
     denom=denom + b(i)*b(i)
     if ( ermax .lt. junk ) ermax = junk
      write(uio,6) i, prod,b(i), prod-b(i), (prod-b(i))/prod
!    if (prod.ne.0) then
!       error = error + abs((prod-b(i))/prod)
!    elseif(b(i).ne.0) then
!       error = error + abs((prod-b(i))/b(i))
!    endif
     error = error + junk
  enddo
  error = error /(1.*m3)
!write(*,*), "VAL OF ERROR IS: ",error
  sigma=sqrt(num/denom)
  write(uio,3)' Average, largest errors in force,percent deviation=',error,ermax,sigma
6 format(i6,3(1x,g13.6),3x,g11.4)
3 format(a,3(1x,g13.6))
  close(uio)

 end subroutine svd_set
!===================================================
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      implicit none
      INTEGER m,mp,n,np !,NMAX
      REAL(8) a(mp,np),v(np,np),w(np)
!     PARAMETER (NMAX=500)
      INTEGER i,its,j,jj,k,l,nm
!     REAL(8) anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag2
      REAL(8) anorm,c,f,g,h,s,scale,x,y,z,rv1(n),pythag2
      g=0.d0
      scale=0.d0
      anorm=0.d0
      L25: do i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.d0
        s=0.d0
        scale=0.d0
        if(i.le.m)then
          do k=i,m
            scale=scale+abs(a(k,i))
          enddo
          if(scale.ne.0.d0)then
            do k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
            enddo
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do j=l,n
              s=0.d0
              do k=i,m
                s=s+a(k,i)*a(k,j)
              enddo
              f=s/h
              do k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
              enddo
            enddo
            do k=i,m
              a(k,i)=scale*a(k,i)
            enddo
          endif
        endif
        w(i)=scale *g
        g=0.d0
        s=0.d0
        scale=0.d0
        if((i.le.m).and.(i.ne.n))then
          do k=l,n
            scale=scale+abs(a(i,k))
          enddo
          if(scale.ne.0.d0)then
            do k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
            enddo
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do k=l,n
              rv1(k)=a(i,k)/h
            enddo
            do j=l,m
              s=0.d0
              do k=l,n
                s=s+a(j,k)*a(i,k)
              enddo
              do k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
              enddo
            enddo
            do k=l,n
              a(i,k)=scale*a(i,k)
            enddo
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
      enddo L25
      L32: do i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.d0)then
            do j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
            enddo
            do j=l,n
              s=0.d0
              do k=l,n
                s=s+a(i,k)*v(k,j)
              enddo
              do k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
              enddo
            enddo
          endif
          do j=l,n
            v(i,j)=0.d0
            v(j,i)=0.d0
          enddo
        endif
        v(i,i)=1.d0
        g=rv1(i)
        l=i
      enddo L32
      L39: do i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do j=l,n
          a(i,j)=0.d0
        enddo
        if(g.ne.0.d0)then
          g=1.0/g
          do j=l,n
            s=0.d0
            do k=l,m
              s=s+a(k,i)*a(k,j)
            enddo
            f=(s/a(i,i))*g
            do k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
            enddo
          enddo
          do j=i,m
            a(j,i)=a(j,i)*g
          enddo
        else
          do j= i,m
            a(j,i)=0.d0
          enddo
        endif
        a(i,i)=a(i,i)+1.d0
      enddo L39
      L49:do k=n,1,-1
        do its=1,30
          L41: do l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  exit L41 !goto 1
          enddo L41
!1         c=0.0
          c=0.d0
          s=1.d0
          L43: do i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) exit L43 ! goto 2
            g=w(i)
            h=pythag2(f,g)
            w(i)=h
            h=1.d0/h
            c= (g*h)
            s=-(f*h)
            do j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
            enddo
          enddo L43
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.d0)then
              w(k)=-z
              do j=1,n
                v(j,k)=-v(j,k)
              enddo
            endif
            cycle L49  !goto 3
          endif
          if(its.eq.30) write(*,*) 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.d0*h*y)
          g=pythag2(f,1d0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag2(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
            enddo
            z=pythag2(f,h)
            w(j)=z
            if(z.ne.0.d0)then
              z=1.d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
            enddo
          enddo
          rv1(l)=0.d0
          rv1(k)=f
          w(k)=x
        enddo
! 3       continue
      enddo L49

      END SUBROUTINE svdcmp
!  (C) Copr. 1986-92 Numerical Recipes Software !+!).
!===================================================
      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      implicit none
      INTEGER m,mp,n,np
      REAL(8) b(mp),u(mp,np),v(np,np),w(np),x(np)
!     PARAMETER (NMAX=500)
      INTEGER i,j,jj
!     REAL(8) s,tmp(NMAX)
      REAL(8) s,tmp(n)

      do j=1,n
        s=0d0
        if(w(j).ne.0d0)then
          do i=1,m
            s=s+u(i,j)*b(i)
          enddo
          s=s/w(j)
        endif
        tmp(j)=s
      enddo
      do j=1,n
        s=0d0
        do jj=1,n
          s=s+v(j,jj)*tmp(jj)
        enddo
        x(j)=s
      enddo

      END SUBROUTINE svbksb
!  (C) Copr. 1986-92 Numerical Recipes Software !+!)
!===================================================
      FUNCTION pythag2(a,b)
      implicit none
      REAL(8) a,b,pythag2
      REAL(8) absa,absb,rat

      absa=abs(a)
      absb=abs(b)

      if(absa.gt.absb)then
      rat = absb/absa
      pythag2=absa*sqrt(1d0+rat*rat)
      else
        if(absb.eq.0d0)then
          pythag2=0d0
        else
          rat = absa/absb
          pythag2=absb*sqrt(1d0+rat*rat)
        endif
      endif

      END FUNCTION pythag2
!  (C) Copr. 1986-92 Numerical Recipes Software !+!).
