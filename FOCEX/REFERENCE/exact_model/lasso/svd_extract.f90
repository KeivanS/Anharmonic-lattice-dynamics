 module constants
 integer, parameter :: r15 = kind(1d0)
 end module constants
      SUBROUTINE svdcmp(a,m,n,w,v)
!! performs SVD decomposition of A in  the form A=U*W*transpose(V)
!! on output U is written into A
 use constants, only : r15
      implicit none
      INTEGER, intent(in) :: m,n
      real(r15), intent(inout) :: a(m,n)
      real(r15), intent(out) :: v(n,n),w(n)
      INTEGER i,its,j,jj,k,l,nm
      real(r15) anorm,c,f,g,h,s,scale,x,y,z,rv1(n),pythag2

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
      FUNCTION pythag2(a,b)
 use constants, only : r15
      implicit none
      real(r15) a,b,pythag2
      real(r15) absa,absb,rat

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
