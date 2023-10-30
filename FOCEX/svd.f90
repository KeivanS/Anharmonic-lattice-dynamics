!===========================================================
 subroutine svd_set(m3,n,a,b,x,sigma,svdcut,error,ermax,sig,fnsvd)
!! performs a Singular Value Decomposition of the matrix a(m3,n) in order to solve the linear
!! system ax=b. SVD results will be written in the file fnsvd
 implicit none
 integer , intent(in):: m3,n
 real(8), intent(inout) :: a(m3,n)
 real(8), intent(in) :: b(m3),svdcut
! real(8), allocatable :: u(:,:),v(:,:),w(:) 
 real(8) v(n,n),w(n),w2(n,n) , ax(m3)
 real(8), intent(out) :: error,ermax,x(n),sigma(n),sig
 character(LEN=*), intent(in) :: fnsvd
 integer i,j,k,uio,umat
 real(8) prod,wmax,wmin,wcut,junk,num,denom

 uio = 345
 umat= 346
 open(uio,file=fnsvd,status='unknown')
 open(umat,file="matrix_"//fnsvd,status='unknown')
! eye=0
! do i=1,n
!    eye(i,i)=1d0   ! make the identity matrix
! enddo

 write(uio,*)' A and b written, now doing the SVD '

! allocate( v(n,n),w(n),u(m3,n) )

!u = a
 call svdcmp(a,m3,n,w,v)

! write svd matrices in a file
 write(uio,*)' W (better not be too small) is:  '
 write(umat,*)' W (better not be too small) is:  '
 do k=1,n
    write(uio,*)'k,w(k)=',k,w(k)
    write(umat,*)'k,w(k)=',k,w(k)
 enddo
 wmin=minval(abs(w))
 wmax=maxval(abs(w))
 wcut = wmax*svdcut
! write(uio,*)' Its condition number w_max/w_min = ',wmax/wmin
 write(uio,*)' w larger than ',wcut,' will be used in the inversion'
! write(umat,*)' Its condition number w_max/w_min = ',wmax/wmin
 write(umat,*)' w larger than ',wcut,' will be used in the inversion'

! write(umat,*)' Matrix u is  '
! do i=1,m3
!   write(umat,9)i,(u(i,j),j=1,n)
! enddo
! write(umat,*)' Matrix v is  '
! do i=1,n
!   write(umat,9)i,(v(i,j),j=1,n)
! enddo

! test for orthogonality of u
! aux=matmul(transpose(u),u)
! prod = maxval(abs(aux-eye))
! if (prod .gt. 1d-8) then
!   write(umat,*) ' SVD: transpose(u)*u is not identity ',prod
!   stop
! endif

! test for orthogonality of v
! aux=matmul(transpose(v),v)
! prod = maxval(abs(aux-eye))
! if (prod .gt. 1d-8) then
!   write(umat,*) ' SVD: transpose(v)*v is not identity ',prod
!   stop
! endif

! test for correctness of SVD by comparing A to u*w*transpose(v)
! write(uio,*)' SVD product of U*W*V^T, a is: '
!  do i=1,m3
!  do j=1,n
!     prod = 0
!     do k=1,n
!        prod = prod+v(j,k)*u(i,k)*w(k)
!     enddo
!     if( abs(prod-a(i,j)) .gt. 1d-10) write(uio,8)'i,j,UW(V^T),a(i,j)=',i,j,prod,a(i,j)
!  enddo
!  enddo


! this is to eliminate the v vectors with small w
  do j=1,n
     if(abs(w(j)).lt.wcut) w(j)=0d0
  enddo

! Solve the system using x=V(1/W)(U^T) b; after eliminating w=0 terms
  call svbksb(a,w,v,m3,n,b,x)

  write(uio,*)' results of SVD solution, variance, error are: x,sigma,Dx'
  w2=0
  do j=1,n
     sigma(j) = 0
     w2(j,j)=w(j)
     do k=1,n
        if (w(k).ge.wcut) then
           sigma(j) = sigma(j) + v(j,k)*v(j,k) / w(k)/w(k)
        endif
     enddo
     sigma(j) = sqrt(sigma(j))
     write(uio,7) j,x(j),sigma(j),sigma(j)/sqrt(m3*1.)
  enddo

! deallocate(v,w,u) !,aux,eye)

! reconstruct a from u,v,w2: a=u*w2*v^T
 w2=matmul(w2,transpose(v))
 a=matmul(a,w2) !matmul(w2,transpose(v)))
 write(uio,*)' residual of SVD solution is: Ax,b,|Ax-b|,|Ax-b|/Ax'
!  error = 0; ermax = 0; num=0; denom=0
  ax=matmul(a,x)
  num=dot_product((ax-b),(ax-b))
  denom=dot_product(ax,ax)
  error=sum(abs(ax-b))/m3
  ermax=maxval(abs(ax-b))
  do i=1,m3
     write(uio,*)' Ax-b ',i,ax(i)-b(i)
  enddo
  sig=sqrt(num/denom)*100
  write(uio,3)' Average, largest errors in force,percent deviation=',error,ermax,sig

6 format(i6,3(1x,g13.6),3x,g11.4)
3 format(a,3(1x,g13.6))
9 format(i8,1x,999(1x,f7.3))
8 format(a,i8,1x,i8,3(1x,f17.10))
7 format(i8,1x,99(1x,g11.4))

  close(uio)
  close(umat)

 end subroutine svd_set
!===================================================
 subroutine solve_by_elimination(m3,n,m0,a,b,x,sigma,svdcut,nconfigs,nlines,energies,itemp,tempk)
!! for a system of linear homogeneous equations, A(m3,n)*X(n)=b , this subroutine finds the
!! range of A and eliminates dependent variables X using SVD, first by solving the
!! inhomogeneous part of Ax=b (i.e. b\=0), and then projecting the solution on the
!! kernel of the homogeneous part (i.e. b=0)
 implicit none
 integer, intent(in) ::  m3,n,m0 ,itemp,nconfigs,nlines(nconfigs)  ! the first m0 lines form the homogeneous part of a
 real(8), intent(in) ::  a(m3,n),b(m3),svdcut ,tempk,energies(nconfigs)
 real(8), intent(out) :: x(n),sigma(n)
 real(8), allocatable :: ahom(:,:),ainhom(:,:),kernelbasis(:,:)
 real(8), allocatable :: mat(:,:),qmat(:)
 integer i,uio,nkernel
 real(8) norm,error,ermax,num,denom,prod,junk,sig

 write(*,*)'SOLVE_BY_ELIMINATION: First solving the kernel of the homogeneous part'

 uio = 345
 open(uio,file='elimination.dat',status='unknown')
 if(m0.ne.0) then  ! there is a homogeneous part

    norm=maxval(abs(b(1:m0)))
    if(norm.gt.1d-8) then
       write(uio,*)'SOLVE_BY_ELIMINATION: the homogeneous part of B does not seem to be zero ',norm
       stop
    endif
    allocate (kernelbasis(n,n),ahom(1:m0,1:n))  ! dummy for the homogeneous part
    ahom=a(1:m0,1:n)
    call get_kernel(m0,n,ahom,svdcut,nkernel,kernelbasis,uio)
    deallocate (ahom)
 endif

 write(*,*)'SOLVE_BY_ELIMINATION: Now solving the inhomogeneous part'
 allocate (ainhom(m3-m0,n))
 ainhom=a(m0+1:m3,1:n)

 if(itemp.eq.1) then  ! implement temperature

    write(*,*) ' Temperature of ',tempk,' will be implemented'
    allocate(mat(n,n),qmat(n))
    write(uio,*)' Boltzmann weighting before solving for the inhomogeneous part'
    call implement_temperature(m3,n,ainhom,b(m0+1:m3),nconfigs,energies,nlines,tempk,mat,qmat)
    call solve_svd(n,n,mat,qmat,svdcut,x,sigma,uio)

 else

    call solve_svd(m3-m0,n,ainhom,b(m0+1:m3),svdcut,x,sigma,uio)

 endif
 write(uio,*)' Solution of the inhomogeneous part is=',x

 deallocate(ainhom)

 write(uio,*)' After solving for the inhomogeneous part, now we project on the kernel of a_hom'
 call project_on(n,nkernel,x,kernelbasis(1:n,1:nkernel),x)
 write(uio,*)' Solution after kernel projection   is=',x

 deallocate(kernelbasis)

! this part can be removed
 write(uio,*)' SOLVE_BY_ELIMINATION: residual of SVD2 solution is: Ax,b,|Ax-b|,|Ax-b|/Ax'
 error = 0; ermax = 0; num=0; denom=0
 do i=1,m3
    prod = dot_product(a(i,:),x)
    junk = abs(prod-b(i))
    num=num     + junk*junk
    denom=denom + prod*prod
    if ( ermax .lt. junk ) ermax = junk
    write(uio,6) i, prod,b(i), junk, junk/prod
    error = error + junk
 enddo
 error = error /(1.*m3)
 sig=sqrt(num/denom)
 write(uio,3)' ELIMINATION: Average, largest errors in force,percent deviation=',error,ermax,sig*100,' %'
 write(uio,*)' After elimination, || F_dft-F_fit || / || F_dft || =',sig

 close(uio)

7 format(i8,1x,99(1x,g13.6))
6 format(i6,3(1x,g13.6),3x,g11.4)
4 format(a,99(1x,g11.4))
3 format(a,3(1x,g11.4),a)

 end subroutine solve_by_elimination
!===================================================
      SUBROUTINE svdcmp(a,m,n,w,v)
!! performs SVD decomposition of A in  the form A=U*W*transpose(V)
!! on output U is written into A
      implicit none
      INTEGER, intent(in) :: m,n
      REAL(8), intent(inout) :: a(m,n)
      real(8), intent(out) :: v(n,n),w(n)
      INTEGER i,its,j,jj,k,l,nm
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
      SUBROUTINE isvbksb(u,w,v,m,n,b,x)
      implicit none
      INTEGER m,n
      REAL(8) b(m),u(m,n),v(n,n),w(n),x(n)
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

      END SUBROUTINE isvbksb
!=================
      SUBROUTINE svbksb(u,w,v,m,n,b,x)
!! back substitutes u,w,v to solve for x, given b after small W terms are set to zero
      implicit none
      INTEGER, intent(in) :: m,n
      REAL(8), intent(in) :: b(m),u(m,n),v(n,n),w(n)
      REAL(8), intent(out) :: x(n)
      INTEGER j
      REAL(8) tmp(n)

      tmp=matmul(transpose(u),b)
      do j=1,n
        if(w(j).ne.0d0)then
           tmp(j)=tmp(j)/w(j)
         endif
      enddo
      x= matmul(v,tmp)

      END SUBROUTINE svbksb
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
!===================================================
 subroutine solve_svd(ndim,n,a,b,svdcut,xout,sig,uio)
!! solves the linear problem ax=b  using svd
!! inputs are a,b,svdcut,a,b and their dimensions; output is xout
!! write resutls and sd=sig in file with unit=uio
 use ios
 use params, only : verbose
 implicit none
 integer, intent(in) :: ndim,n,uio
 real(8), intent(in) :: b(ndim),svdcut
 real(8), intent(in) :: a(ndim,n)
 real(8), intent(out):: xout(n) ,sig(n)
 integer i,k
 real(8), allocatable :: u(:,:),w(:),v(:,:)
 real(8) wmax,wmin,wcut

 allocate( u(ndim,n),w(n),v(n,n) )
 write(uio,*)'SOLVE_SVD: amatrix is of dimensions ',ndim,n
 if (verbose) call write_out(uio,'=========================== SOLVE_SVD: amatrix ',a)
 u=a ! to keep a intact
 call svdcmp(u,ndim,n,w,v)  ! BEWARE: a will be overwritten !

! make small w equal to zero
 wmax=maxval(abs(w))
 wmin=minval(abs(w))
 wcut=wmax*svdcut
 write(uio,4)'SOLVE_SVD: largest w, condition number is ',wmax,wmax/wmin
 write(uio,4)'SOLVE_SVD: will only keep svs larger than ',wcut
 call write_out(uio,'=========================== SOLVE_SVD: W ',w)

! this is to eliminate the v vectors with small w
 do i=1,n
    if(abs(w(i)).lt.wcut) w(i)=0d0
 enddo

 call svbksb(u,w,v,ndim,n,b,xout)

 write(uio,*)' SOLVE_SVD: solution, variance, error are: x,sigma,Dx'
 do i=1,n
    sig(i) = 0
    do k=1,n
       if (w(k).ge.wcut) then
          sig(i) = sig(i) + v(i,k)*v(i,k) / w(k)/w(k)
       endif
    enddo
    sig(i) = sqrt(sig(i))
    write(uio,7) i,xout(i),sig(i),sig(i)/sqrt(ndim*1.)
 enddo

 deallocate(u,v,w)

4 format(a,99(1x,g11.4))
7 format(i8,1x,g13.6,99(1x,g11.4))

 end subroutine solve_svd
!======================================================
 subroutine get_kernel(m0,n,amat,svdcut,nkernel,kernel,uio)
 use ios
 implicit none
 integer, intent(in) :: m0,n,uio
 integer, intent(out):: nkernel
 real(8), intent(in) :: amat(m0,n)
 real(8), intent(out):: kernel(n,n)
 integer j,nb
 integer,allocatable :: r(:)
 real(8),allocatable :: u(:,:),v(:,:),w(:)
 real(8) wmin,wmax,wcut,svdcut

 write(*,*)'GET_KERNEL: calculating the kernel space of the homogeneous part of amat'

 allocate(u(m0,n),v(n,n),w(n),r(n))

! SVD the homogeneous part of A to find the kernel basis and project x2 in it
    u=amat(1:m0,1:n)
    call svdcmp(u,m0,n,w,v)

    wmax=maxval(abs(w))
    wmin=minval(abs(w))
    wcut=wmax*svdcut
!    write(uio,4)'SOLVE_BY_ELIMINATION: largest w, condition number is ',wmax,wmax/wmin
    write(uio,4)'Maxval of w is=',wmax
    write(uio,4)'SOLVE_BY_ELIMINATION: will only keep svds larger than wmax*1d-9=',wcut

    r=0
    do j=1,n
       write(uio,*)'j,w(j)=',j,w(j)
       if(abs(w(j)).lt.wcut) then  ! kernel space
          w(j)=0d0
          r(j)=1  ! r is a map of zero elements of w , where x can be arbitrary (free)
       else
          r(j)=0  ! 1-r is a map of non-zero elements of w , x should be zero in this subspace
       endif
    enddo
    nkernel = sum(r)
    write(uio,*)' Kernel of A has dimension=',nkernel
    write(uio,*)' This is the number of independent variables or  equations '
    write(uio,*)' The number of eliminated components is ',n-nkernel

! identify the Kernel of w
! The space of nonzero w should have X=0 since we are solving a homogeneous set of equations
! but if w is zero, then x is free to move in that subspace
! also store vector components that belong to the kernel of w (where x is free)
! the solution of the inhomogeneous part belongs to the Kernel of the Homogeneous part.

    call write_out(uio,'#####  GET_KERNEL: Umatrix',u)
    call write_out(uio,'#####  GET_KERNEL: Vmatrix',v)
    kernel=0
    nb=0
    do j=1,n
! this is to eliminate the v vectors with small w, and count the range of a
       if(r(j).eq.1) then
         nb=nb+1
         kernel(:,nb)=v(:,j)
       endif
    enddo
    if(nb.ne.nkernel) then
      write(*  ,*)'ERROR in GET_KERNEL, nb,nkernel=',nb,nkernel
      write(uio,*)'ERROR in GET_KERNEL, nb,nkernel=',nb,nkernel
      stop
    endif
    call write_out(uio,'#####  GET_KERNEL: kernel_basis',kernel)

4 format(a,99(1x,g11.4))

 end subroutine get_kernel
