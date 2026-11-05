!============================================================
 module ewald_cell
!! Ewald energy, forces and strain derivative of the pair model E = 1/2 sum_ij c_ij phi_eps(x_j - x_i + R) in a periodic
!! cell of any shape (2026-11-04, for relax_reference with born_flag 2). phi_eps is the potential of dewapot
!! (ewald_born.f90): real space erfc(a D)/(D sqrt det eps), D = sqrt(r.eps^-1.r); reciprocal
!! (4 pi/V) sum_{G/=0} exp(-G.eps.G/4a^2)/(G.eps.G) cos(G.r); self term -a/(sqrt(pi) sqrt det eps) sum_i c_ii;
!! background -pi/(2 V a^2) sum_ij c_ij (0 for a neutral cell). eps is held fixed under the deformation.
!! Energy in units of c times e^2/(4 pi eps0) per A: the caller multiplies by that factor.
 use constants, only : r15, pi
 implicit none
 contains

 subroutine ewald_cell_efw(n,pos,cell,cc,eps,alpha,e,f,w)
!! pos(3,n) Cartesian (A), cell(3,3) lattice vectors in columns, cc(n,n) symmetric couplings, alpha (1/A).
!! e: energy; f(3,n) = -dE/dpos; w(a,b) = dE/d eps_ab for the homogeneous deformation pos -> (1+eps) pos,
!! cell -> (1+eps) cell of the current configuration (eps fixed dielectric tensor not deformed).
 integer, intent(in) :: n
 real(r15), intent(in) :: pos(3,n),cell(3,3),cc(n,n),eps(3,3),alpha
 real(r15), intent(out) :: e,f(3,n),w(3,3)
 real(r15) r0(3),ei(3,3),sd,vol,b(3,3),tre,lmin,rc,gc,r(3),u(3),d,ad,ph,dph,gv(3),eg(3),s,wg,dwg,s2,fac
 real(r15) cs(n),sn(n),t(n),csum
 integer nr(3),ng(3),i,j,i1,i2,i3,a
 call inv3(eps,ei); ei=(ei+transpose(ei))/2
 sd=sqrt(det3x(eps)); vol=abs(det3x(cell))
 call inv3(cell,b); b=2*pi*transpose(b)          ! reciprocal vectors in columns: b(:,i).cell(:,j) = 2 pi delta_ij
 tre=eps(1,1)+eps(2,2)+eps(3,3)
 lmin=det3x(eps)/tre**2                         ! lower bound of the smallest eigenvalue of eps
 rc=6.5d0*sqrt(tre)/alpha                       ! a D > 6.5 beyond rc (D >= |r|/sqrt(lambda_max))
 gc=13d0*alpha/sqrt(lmin)                       ! G.eps.G > 169 a^2 beyond gc
 do i=1,3
    nr(i)=int(rc*sqrt(sum(b(:,i)**2))/(2*pi))+2
    ng(i)=int(gc*sqrt(sum(cell(:,i)**2))/(2*pi))+1
 enddo
 e=0; f=0; w=0
! real space
 do i=1,n
 do j=1,n
    if(cc(i,j).eq.0) cycle
    r0=pos(:,j)-pos(:,i)
    r0=r0-matmul(cell,dble(nint(matmul(transpose(b),r0)/(2*pi))))   ! nearest image, so |frac| <= 1/2
    do i1=-nr(1),nr(1); do i2=-nr(2),nr(2); do i3=-nr(3),nr(3)
       r=r0+i1*cell(:,1)+i2*cell(:,2)+i3*cell(:,3)
       if(sum(r*r).lt.1d-20) cycle                                   ! the atom itself (self term below)
       u=matmul(ei,r)
       d=sqrt(dot_product(r,u)); ad=alpha*d
       if(ad.gt.6.5d0) cycle
       ph=erfc(ad)/(d*sd)
       dph=(-erfc(ad)/d**2-2*alpha/sqrt(pi)*exp(-ad*ad)/d)/sd      ! d phi / d D
       e=e+cc(i,j)*ph/2
       f(:,i)=f(:,i)+cc(i,j)*dph/d*u
       do a=1,3
          w(a,:)=w(a,:)+cc(i,j)*dph/d*u(a)*r/2
       enddo
    enddo; enddo; enddo
 enddo
 enddo
! reciprocal space
 do i1=-ng(1),ng(1); do i2=-ng(2),ng(2); do i3=-ng(3),ng(3)
    if(i1.eq.0 .and. i2.eq.0 .and. i3.eq.0) cycle
    gv=i1*b(:,1)+i2*b(:,2)+i3*b(:,3)
    eg=matmul(eps,gv); s=dot_product(gv,eg)
    if(s.gt.(gc**2)*lmin) cycle
    wg=4*pi/vol*exp(-s/(4*alpha**2))/s
    dwg=-wg*(1/(4*alpha**2)+1/s)                       ! d wg / d s at fixed V
    do i=1,n
       ph=dot_product(gv,pos(:,i)); cs(i)=cos(ph); sn(i)=sin(ph)
    enddo
    s2=0
    do i=1,n
       t(i)=dot_product(cc(i,:),cs)*cs(i)+dot_product(cc(i,:),sn)*sn(i)   ! sum_j c_ij cos(G.(x_j-x_i))
       s2=s2+t(i)
       fac=dot_product(cc(i,:),sn)*cs(i)-dot_product(cc(i,:),cs)*sn(i)   ! sum_j c_ij sin(G.(x_j-x_i))
       f(:,i)=f(:,i)-wg*fac*gv
    enddo
    e=e+wg*s2/2
    do a=1,3
       w(a,a)=w(a,a)-wg*s2/2
       w(a,:)=w(a,:)-dwg*s2*gv(a)*eg
    enddo
 enddo; enddo; enddo
! self and background terms
 csum=sum(cc)
 do i=1,n
    e=e-alpha*cc(i,i)/(sqrt(pi)*sd)
 enddo
 e=e-pi/(2*vol*alpha**2)*csum
 do a=1,3
    w(a,a)=w(a,a)+pi/(2*vol*alpha**2)*csum
 enddo
 end subroutine ewald_cell_efw

 function det3x(m) result(dt)
 real(r15), intent(in) :: m(3,3)
 real(r15) dt
 dt=m(1,1)*(m(2,2)*m(3,3)-m(2,3)*m(3,2))-m(1,2)*(m(2,1)*m(3,3)-m(2,3)*m(3,1))+m(1,3)*(m(2,1)*m(3,2)-m(2,2)*m(3,1))
 end function det3x

 subroutine inv3(m,mi)
 real(r15), intent(in) :: m(3,3)
 real(r15), intent(out) :: mi(3,3)
 real(r15) dt
 dt=det3x(m)
 mi(1,1)=m(2,2)*m(3,3)-m(2,3)*m(3,2); mi(1,2)=m(1,3)*m(3,2)-m(1,2)*m(3,3); mi(1,3)=m(1,2)*m(2,3)-m(1,3)*m(2,2)
 mi(2,1)=m(2,3)*m(3,1)-m(2,1)*m(3,3); mi(2,2)=m(1,1)*m(3,3)-m(1,3)*m(3,1); mi(2,3)=m(1,3)*m(2,1)-m(1,1)*m(2,3)
 mi(3,1)=m(2,1)*m(3,2)-m(2,2)*m(3,1); mi(3,2)=m(1,2)*m(3,1)-m(1,1)*m(3,2); mi(3,3)=m(1,1)*m(2,2)-m(1,2)*m(2,1)
 mi=mi/dt
 end subroutine inv3

 end module ewald_cell
