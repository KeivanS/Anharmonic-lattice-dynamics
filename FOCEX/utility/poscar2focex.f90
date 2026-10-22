!! poscar2focex -- POSCAR -> FOCEX structure.params, POSCAR_focex and kpbs.params, with the cell
!! and the band-structure path in the convention of Setyawan and Curtarolo,
!! Comput. Mater. Sci. 49, 299 (2010) (as in pymatgen's KPathSetyawanCurtarolo).
!!
!!    poscar2focex.x POSCAR [-o structure.params] [-p POSCAR_focex] [-k kpbs.params]
!!                          [-s symprec] [-f]
!!
!! 1. Niggli reduction of the input cell (Krivy and Gruber 1976, with the epsilon tests of
!!    Grosse-Kunstleve et al. 2004).
!! 2. Symmetry: the lattice operations (integer matrices that keep the metric of the reduced
!!    cell) and, among them, those of the crystal (that map every atom onto an atom of the
!!    same type after some translation; tolerance symprec, default 1e-3 A).
!! 3. Crystal system from the rotations of the crystal; conventional cell from its rotation
!!    axes (shortest lattice vectors along them); centring from the lattice points of that
!!    cell; then the Setyawan-Curtarolo choice of axes:
!!      cubic           a = b = c along the 4-fold (or, for T and Th, the 2-fold) axes
!!      tetragonal      c along the 4-fold axis, a the shortest lattice vector normal to it
!!      orthorhombic    along the three 2-fold axes, a < b < c; for a base-centred lattice c
!!                      is normal to the centred face and a < b
!!      hexagonal       c along the 6- or 3-fold axis, a the shortest lattice vector normal to
!!                      it, b = a turned by +120 degrees; rhombohedral lattices in the
!!                      obverse hexagonal setting
!!      monoclinic      a along the 2-fold axis (unique axis a, as in Setyawan-Curtarolo),
!!                      b, c normal to it with alpha < 90 degrees; a centred lattice is
!!                      written with the centring vector (a+b)/2
!!      triclinic       the reduced cell, with the signs that make the reciprocal angles all
!!                      acute or all obtuse
!!    The primitive cell is the Setyawan-Curtarolo one (e.g. fcc (0,1/2,1/2),(1/2,0,1/2),
!!    (1/2,1/2,0) in conventional units).
!! 4. Writes, in the frame FOCEX builds from a,b,c,alpha,beta,gamma (a along x, b in xy):
!!      structure.params  line 1 a b c alpha beta gamma (A, degrees) of the conventional
!!                        cell, line 2 the primitive vectors in conventional units, line 3
!!                        "1.0 0", atoms in conventional units; masses and names from the
!!                        POSCAR, the other values the defaults of FOCEX/example/Ge, marked
!!                        "default, check"
!!      POSCAR_focex      the primitive cell (build the DFT supercells from it)
!!      kpbs.params       the Setyawan-Curtarolo path in primitive reduced coordinates
!!    The atom order of the POSCAR is kept.  A VASP-5 POSCAR (with the element line) is needed.
 program poscar2focex
 implicit none
 integer, parameter :: dp=kind(1d0), maxat=2000, maxop=48, maxkp=40
 real(dp), parameter :: pi=3.14159265358979323846_dp, keytol=1d-5
 character(len=256) :: fin,fout,fpos,fkp,arg,title
 character(len=8) :: species(100)
 integer :: counts(100),ntyp,nat,typ(maxat)
 real(dp) :: A(3,3),F(3,maxat),B(3,3),FB(3,maxat),symprec
 logical :: force,ex
 integer :: i,k,iarg,nlat,ncry,ncen
 integer :: Mlat(3,3,maxop),Wcry(3,3,maxop),Nconv(3,3),Mp(3,3)
 real(dp) :: tcry(3,maxop),C(3,3),T(3,3),Pc(3,3),Q(3,3),latp(6),Fr(3,3),Pr(3,3)
 real(dp) :: fp(3,maxat),xc(3,maxat),x(3)
 character(len=16) :: system,bravais,scname
 character(len=4) :: kname(maxkp)
 real(dp) :: kred(3,maxkp)
 integer :: nk,npath
 character(len=4) :: path(maxkp)

 fout='structure.params' ; fpos='POSCAR_focex' ; fkp='kpbs.params'
 symprec=1d-3 ; force=.false. ; fin=''
 iarg=1
 do while (iarg.le.command_argument_count())
    call get_command_argument(iarg,arg)
    select case (trim(arg))
    case ('-o') ; iarg=iarg+1 ; call get_command_argument(iarg,fout)
    case ('-p') ; iarg=iarg+1 ; call get_command_argument(iarg,fpos)
    case ('-k') ; iarg=iarg+1 ; call get_command_argument(iarg,fkp)
    case ('-s') ; iarg=iarg+1 ; call get_command_argument(iarg,arg) ; read(arg,*) symprec
    case ('-f') ; force=.true.
    case ('-h','--help')
       write(*,*)'usage: poscar2focex.x POSCAR [-o structure.params] [-p POSCAR_focex] [-k kpbs.params]'
       write(*,*)'                         [-s symprec(A), default 1e-3] [-f (overwrite)]'
       stop
    case default ; fin=arg
    end select
    iarg=iarg+1
 enddo
 if (fin.eq.'') error stop 'usage: poscar2focex.x POSCAR [-o ..] [-p ..] [-k ..] [-s symprec] [-f]; -h for help'
 do i=1,3
    if (i.eq.1) arg=fout
    if (i.eq.2) arg=fpos
    if (i.eq.3) arg=fkp
    inquire(file=trim(arg),exist=ex)
    if (ex .and. .not.force) then
       write(*,*) trim(arg),' exists; use -f to overwrite it'
       error stop 1
    endif
 enddo

 call read_poscar(fin)
 call niggli(A,B,Mp)                       ! B = Mp A
 call to_frac_of(B,nat,F,A,FB)             ! fractional coordinates in the reduced cell
 call lattice_ops(B,nlat,Mlat)
 call crystal_ops(B,nat,FB,typ,nlat,Mlat,ncry,Wcry,tcry)
 call crystal_system(ncry,Wcry,system)
 call conventional(B,ncry,Wcry,system,Nconv,bravais,ncen)
 call sc_primitive(bravais,T)              ! primitive rows in conventional units
 C=matmul(dble(Nconv),B)
 Pc=matmul(T,C)
 ! the S-C primitive cell must generate the lattice of B
 Q=matmul(Pc,inv3(B))
 if (maxval(abs(Q-anint(Q))).gt.1d-6 .or. abs(abs(det3(anint(Q)))-1).gt.1d-6) &
      error stop 'internal error: the Setyawan-Curtarolo primitive cell is not a basis of the lattice'
 if (trim(bravais).eq.'RHL') call rhl_check(Pc)
 ! rotate into the FOCEX frame of the conventional cell
 call lattice_params(C,latp)
 call focex_frame(latp,Fr)
 Q=matmul(inv3(C),Fr)                      ! C Q = Fr
 if (maxval(abs(matmul(Q,transpose(Q))-eye3())).gt.1d-8 .or. det3(Q).lt.0) &
      error stop 'cannot rotate the cell into the FOCEX frame (left-handed conventional cell?)'
 Pr=matmul(Pc,Q)
 ! atoms: fractional in the primitive cell, then conventional units x_conv = f_prim T
 do i=1,nat
    x=matmul(FB(:,i),B)                    ! Cartesian (input frame)
    fp(:,i)=matmul(x,inv3(Pc))
    fp(:,i)=fp(:,i)-floor(fp(:,i)+1d-10)
    xc(:,i)=matmul(fp(:,i),T)
 enddo
 call canonical_setting(T,latp)
 call sc_kpath(bravais,C,Pc,scname)

 call write_structure(fout,latp,T)
 call write_poscar(fpos,Pr)
 call write_kpbs(fkp)
 call report()

 contains

!--------------------------------------------------------------- input
 subroutine read_poscar(fn)
 character(len=*), intent(in) :: fn
 character(len=512) :: line
 real(dp) :: s,X0(3,maxat)
 logical :: cart
 integer :: u,ios,n,i0
 open(newunit=u,file=trim(fn),status='old',iostat=ios)
 if (ios.ne.0) then
    write(*,*) 'cannot open ',trim(fn) ; error stop 1
 endif
 read(u,'(a)') title
 read(u,*) s
 do i=1,3
    read(u,*) A(i,:)
 enddo
 if (s.lt.0) s=(-s/abs(det3(A)))**(1d0/3d0)
 A=A*s
 read(u,'(a)') line
 if (len_trim(line).eq.0 .or. scan(adjustl(line),'0123456789').eq.1) &
      error stop 'no element line (line 6) in the POSCAR: a VASP-5 POSCAR is needed'
 call split_words(line,ntyp,species)
 read(u,'(a)') line
 counts=0
 read(line,*,iostat=ios) counts(1:ntyp)
 if (ios.ne.0) error stop 'line 7 of the POSCAR must hold one number of atoms per element'
 nat=sum(counts(1:ntyp))
 if (nat.gt.maxat) error stop 'too many atoms'
 read(u,'(a)') line
 line=adjustl(line)
 if (line(1:1).eq.'S' .or. line(1:1).eq.'s') then
    read(u,'(a)') line ; line=adjustl(line)
 endif
 cart = (index('CcKk',line(1:1)).gt.0)
 n=0
 do i=1,ntyp
    do i0=1,counts(i)
       n=n+1
       typ(n)=i
       read(u,*) X0(:,n)
    enddo
 enddo
 close(u)
 do i=1,nat
    if (cart) then
       F(:,i)=matmul(X0(:,i)*s,inv3(A))
    else
       F(:,i)=X0(:,i)
    endif
 enddo
 end subroutine read_poscar

 subroutine split_words(line,n,w)
 character(len=*), intent(in) :: line
 integer, intent(out) :: n
 character(len=*), intent(out) :: w(:)
 integer :: p,q,lt
 n=0 ; p=1 ; lt=len_trim(line)
 do
    do while (p.le.lt)                       ! skip blanks (no short-circuit .and. in Fortran)
       if (line(p:p).ne.' ' .and. line(p:p).ne.char(9)) exit
       p=p+1
    enddo
    if (p.gt.lt) exit
    if (line(p:p).eq.'!' .or. line(p:p).eq.'#') exit
    q=p
    do while (q.le.lt)
       if (line(q:q).eq.' ' .or. line(q:q).eq.char(9)) exit
       q=q+1
    enddo
    n=n+1
    w(n)=line(p:q-1)
    p=q
 enddo
 end subroutine split_words

!--------------------------------------------------------------- Niggli reduction
 subroutine niggli(A0,Bn,M)
!! Bn = M A0 (rows), Niggli-reduced and right-handed, M integer unimodular
 real(dp), intent(in) :: A0(3,3)
 real(dp), intent(out) :: Bn(3,3)
 integer, intent(out) :: M(3,3)
 real(dp) :: L(3,3),e,Ag,Bg,Cg,xi,eta,zeta,v(3),tt
 integer :: it,s(3),want,i1,j1,k1,sg(2)
 logical :: done
 sg=[1,-1]
 L=A0
 e=1d-5*abs(det3(L))**(2d0/3d0)
 do it=1,1000
    call gm(L,Ag,Bg,Cg,xi,eta,zeta)
    if (Ag.gt.Bg+e .or. (abs(Ag-Bg).le.e .and. abs(xi).gt.abs(eta)+e)) then        ! A1
       v=L(1,:) ; L(1,:)=-L(2,:) ; L(2,:)=-v ; L(3,:)=-L(3,:) ; cycle
    endif
    if (Bg.gt.Cg+e .or. (abs(Bg-Cg).le.e .and. abs(eta).gt.abs(zeta)+e)) then      ! A2
       v=L(2,:) ; L(2,:)=-L(3,:) ; L(3,:)=-v ; L(1,:)=-L(1,:) ; cycle
    endif
    s=[isg(xi,e),isg(eta,e),isg(zeta,e)]
    want=-1
    if (s(1)*s(2)*s(3).eq.1) want=1                                                 ! A3 / A4
    done=.false.
    do i1=1,2 ; do j1=1,2 ; do k1=1,2
       if (done) cycle
       v=[sg(j1)*sg(k1)*xi, sg(i1)*sg(k1)*eta, sg(i1)*sg(j1)*zeta]
       if ((want.eq.1 .and. all(v.gt.e)) .or. (want.eq.-1 .and. all(v.le.e))) then
          L(1,:)=sg(i1)*L(1,:) ; L(2,:)=sg(j1)*L(2,:) ; L(3,:)=sg(k1)*L(3,:)
          done=.true.
       endif
    enddo ; enddo ; enddo
    call gm(L,Ag,Bg,Cg,xi,eta,zeta)
    if (abs(xi).gt.Bg+e .or. (abs(xi-Bg).le.e .and. 2*eta.lt.zeta-e) .or. (abs(xi+Bg).le.e .and. zeta.lt.-e)) then
       L(3,:)=L(3,:)-sign(1d0,xi)*L(2,:) ; cycle                                    ! A5
    endif
    if (abs(eta).gt.Ag+e .or. (abs(eta-Ag).le.e .and. 2*xi.lt.zeta-e) .or. (abs(eta+Ag).le.e .and. zeta.lt.-e)) then
       L(3,:)=L(3,:)-sign(1d0,eta)*L(1,:) ; cycle                                   ! A6
    endif
    if (abs(zeta).gt.Ag+e .or. (abs(zeta-Ag).le.e .and. 2*xi.lt.eta-e) .or. (abs(zeta+Ag).le.e .and. eta.lt.-e)) then
       L(2,:)=L(2,:)-sign(1d0,zeta)*L(1,:) ; cycle                                  ! A7
    endif
    tt=xi+eta+zeta+Ag+Bg
    if (tt.lt.-e .or. (abs(tt).le.e .and. 2*(Ag+eta)+zeta.gt.e)) then
       L(3,:)=L(3,:)+L(1,:)+L(2,:) ; cycle                                          ! A8
    endif
    exit
 enddo
 if (it.gt.1000) error stop 'Niggli reduction did not converge'
 if (det3(L).lt.0) L=-L
 M=nint(matmul(L,inv3(A0)))
 if (abs(abs(det3(dble(M)))-1).gt.1d-6) error stop 'Niggli reduction gave a non-unimodular transformation'
 Bn=matmul(dble(M),A0)
 end subroutine niggli

 subroutine gm(L,Ag,Bg,Cg,xi,eta,zeta)
 real(dp), intent(in) :: L(3,3)
 real(dp), intent(out) :: Ag,Bg,Cg,xi,eta,zeta
 Ag=dot_product(L(1,:),L(1,:)) ; Bg=dot_product(L(2,:),L(2,:)) ; Cg=dot_product(L(3,:),L(3,:))
 xi=2*dot_product(L(2,:),L(3,:)) ; eta=2*dot_product(L(1,:),L(3,:)) ; zeta=2*dot_product(L(1,:),L(2,:))
 end subroutine gm

 integer function isg(x,e)
 real(dp), intent(in) :: x,e
 isg=0
 if (x.gt.e) isg=1
 if (x.lt.-e) isg=-1
 end function isg

!--------------------------------------------------------------- symmetry
 subroutine to_frac_of(Bn,n,F0,A0,Fn)
!! fractional coordinates Fn in cell Bn of the points with fractional coordinates F0 in A0
 real(dp), intent(in) :: Bn(3,3),F0(3,maxat),A0(3,3)
 integer, intent(in) :: n
 real(dp), intent(out) :: Fn(3,maxat)
 integer :: i0
 do i0=1,n
    Fn(:,i0)=matmul(matmul(F0(:,i0),A0),inv3(Bn))
    Fn(:,i0)=Fn(:,i0)-floor(Fn(:,i0)+1d-10)
 enddo
 end subroutine to_frac_of

 subroutine lattice_ops(Bn,nl,Ml)
!! the integer matrices M (entries -2..2) with M G M^T = G: point operations of the lattice.
!! A lattice vector with integer coordinates n (rows of Bn) goes to n M^T... the point with
!! fractional coordinates f goes to f' = M^T f; the Cartesian operation is Bn^T M^T Bn^-T
 real(dp), intent(in) :: Bn(3,3)
 integer, intent(out) :: nl,Ml(3,3,maxop)
 real(dp) :: G(3,3),G2(3,3),tolg
 integer :: m(9),M3(3,3),id,ii
 G=matmul(Bn,transpose(Bn))
 tolg=2*symprec*sqrt(maxval([(G(ii,ii),ii=1,3)]))
 nl=0
 do id=0,5**9-1
    ii=id
    do k=1,9
       m(k)=mod(ii,5)-2 ; ii=ii/5
    enddo
    M3=reshape(m,[3,3])
    if (abs(idet3(M3)).ne.1) cycle
    G2=matmul(matmul(dble(M3),G),transpose(dble(M3)))
    if (maxval(abs(G2-G)).gt.tolg) cycle
    nl=nl+1
    if (nl.gt.maxop) error stop 'more than 48 lattice operations: tolerance too large?'
    Ml(:,:,nl)=M3
 enddo
 end subroutine lattice_ops

 subroutine crystal_ops(Bn,n,Fn,ty,nl,Ml,nop,W,t)
!! lattice operations W = M^T (on fractional coordinates) that, with a translation t, map the
!! crystal onto itself: W f_i + t = f_j (mod 1) with type(j) = type(i), within symprec
 real(dp), intent(in) :: Bn(3,3),Fn(3,maxat)
 integer, intent(in) :: n,ty(maxat),nl,Ml(3,3,maxop)
 integer, intent(out) :: nop,W(3,3,maxop)
 real(dp), intent(out) :: t(3,maxop)
 integer :: il,j0,i0,k0
 real(dp) :: tt(3),y(3),d(3)
 logical :: okall,found
 nop=0
 do il=1,nl
    do j0=1,n
       if (ty(j0).ne.ty(1)) cycle
       tt=Fn(:,j0)-matmul(transpose(dble(Ml(:,:,il))),Fn(:,1))
       okall=.true.
       do i0=1,n
          y=matmul(transpose(dble(Ml(:,:,il))),Fn(:,i0))+tt
          found=.false.
          do k0=1,n
             if (ty(k0).ne.ty(i0)) cycle
             d=y-Fn(:,k0) ; d=d-anint(d)
             if (sqrt(sum(matmul(d,Bn)**2)).lt.symprec) then
                found=.true. ; exit
             endif
          enddo
          if (.not.found) then
             okall=.false. ; exit
          endif
       enddo
       if (okall) then
          nop=nop+1
          W(:,:,nop)=transpose(Ml(:,:,il))
          t(:,nop)=tt-floor(tt+1d-10)
          exit
       endif
    enddo
 enddo
 end subroutine crystal_ops

 integer function rot_order(Wm)
!! order of the proper part det(W) W: 1, 2, 3, 4 or 6, from its trace
 integer, intent(in) :: Wm(3,3)
 integer :: tr
 tr=idet3(Wm)*(Wm(1,1)+Wm(2,2)+Wm(3,3))
 select case (tr)
 case (3) ; rot_order=1
 case (-1) ; rot_order=2
 case (0) ; rot_order=3
 case (1) ; rot_order=4
 case (2) ; rot_order=6
 case default ; error stop 'internal error: not a crystallographic rotation'
 end select
 end function rot_order

 subroutine crystal_system(nop,W,sys)
 integer, intent(in) :: nop,W(3,3,maxop)
 character(len=*), intent(out) :: sys
 integer :: n2,n3,n4,n6,ic
 n2=0 ; n3=0 ; n4=0 ; n6=0
 do ic=1,nop
    select case (rot_order(W(:,:,ic)))
    case (2) ; n2=n2+1
    case (3) ; n3=n3+1
    case (4) ; n4=n4+1
    case (6) ; n6=n6+1
    end select
 enddo
 if (n3.ge.8) then
    sys='cubic'
 elseif (n6.ge.1) then
    sys='hexagonal'
 elseif (n3.ge.2) then
    sys='trigonal'
 elseif (n4.ge.1) then
    sys='tetragonal'
 elseif (n2.ge.3) then
    sys='orthorhombic'
 elseif (n2.ge.1) then
    sys='monoclinic'
 else
    sys='triclinic'
 endif
 end subroutine crystal_system

 subroutine axis_of(Wm,u)
!! u = the shortest lattice vector (integer coordinates in the reduced cell) along the axis of
!! the proper part of Wm: the solution of (Wp - 1) u = 0
 integer, intent(in) :: Wm(3,3)
 integer, intent(out) :: u(3)
 integer :: Wp(3,3),R(3,3),c(3),i1,j1,g
 Wp=idet3(Wm)*Wm
 R=Wp ; do i1=1,3 ; R(i1,i1)=R(i1,i1)-1 ; enddo
 u=0
 do i1=1,3 ; do j1=i1+1,3
    c=[R(i1,2)*R(j1,3)-R(i1,3)*R(j1,2), R(i1,3)*R(j1,1)-R(i1,1)*R(j1,3), R(i1,1)*R(j1,2)-R(i1,2)*R(j1,1)]
    if (any(c.ne.0)) then
       u=c ; exit
    endif
 enddo ; if (any(u.ne.0)) exit ; enddo
 if (all(u.eq.0)) error stop 'internal error: no rotation axis'
 g=max(1,igcd(igcd(abs(u(1)),abs(u(2))),abs(u(3))))
 u=u/g
 do i1=1,3                                  ! a fixed sign: first non-zero component positive
    if (u(i1).ne.0) then
       if (u(i1).lt.0) u=-u
       exit
    endif
 enddo
 end subroutine axis_of

 subroutine axes_of_order(nop,W,ord,na,ax,iop)
!! the distinct axes (na of them) of the crystal rotations of the given order; iop: an
!! operation of that order about each axis
 integer, intent(in) :: nop,W(3,3,maxop),ord
 integer, intent(out) :: na,ax(3,maxop),iop(maxop)
 integer :: ic,u(3),j0
 logical :: new
 na=0
 do ic=1,nop
    if (rot_order(W(:,:,ic)).ne.ord) cycle
    call axis_of(W(:,:,ic),u)
    new=.true.
    do j0=1,na
       if (all(ax(:,j0).eq.u)) new=.false.
    enddo
    if (new) then
       na=na+1 ; ax(:,na)=u ; iop(na)=ic
    endif
 enddo
 end subroutine axes_of_order

 subroutine shortest_normal(Bn,u,v)
!! v: the shortest lattice vector (integer coordinates) normal to the lattice vector u
 real(dp), intent(in) :: Bn(3,3)
 integer, intent(in) :: u(3)
 integer, intent(out) :: v(3)
 integer :: n1,n2,n3,w(3)
 real(dp) :: cu(3),cw(3),best,l
 cu=matmul(dble(u),Bn) ; best=huge(1d0) ; v=0
 do n1=-4,4 ; do n2=-4,4 ; do n3=-4,4
    w=[n1,n2,n3]
    if (all(w.eq.0)) cycle
    cw=matmul(dble(w),Bn)
    l=sqrt(sum(cw**2))
    if (abs(dot_product(cw,cu)).gt.1d-6*l*sqrt(sum(cu**2))) cycle
    if (l.lt.best-1d-8 .or. (abs(l-best).le.1d-8 .and. lexi_gt(cw,matmul(dble(v),Bn)))) then
       best=l ; v=w
    endif
 enddo ; enddo ; enddo
 if (all(v.eq.0)) error stop 'internal error: no lattice vector normal to the axis'
 end subroutine shortest_normal

 logical function lexi_gt(a1,b1)
 real(dp), intent(in) :: a1(:),b1(:)
 integer :: m
 lexi_gt=.false.
 do m=1,3
    if (a1(m).gt.b1(m)+1d-8) then
       lexi_gt=.true. ; return
    elseif (a1(m).lt.b1(m)-1d-8) then
       return
    endif
 enddo
 end function lexi_gt

 logical function is_lattice_point(Nc,f)
!! is the point with fractional coordinates f in the cell Nc (rows, integer coordinates in the
!! reduced cell) a lattice point?
 integer, intent(in) :: Nc(3,3)
 real(dp), intent(in) :: f(3)
 real(dp) :: y(3)
 y=matmul(f,dble(Nc))
 is_lattice_point = all(abs(y-anint(y)).lt.1d-6)
 end function is_lattice_point

!--------------------------------------------------------------- conventional cell
 subroutine conventional(Bn,nop,W,sys,Nc,brav,ncentr)
!! Nc: conventional cell (rows, integer coordinates in the reduced cell Bn), in the
!! Setyawan-Curtarolo choice; brav: CUB FCC BCC TET BCT ORC ORCF ORCI ORCC HEX RHL MCL MCLC TRI
 real(dp), intent(in) :: Bn(3,3)
 integer, intent(in) :: nop,W(3,3,maxop)
 character(len=*), intent(in) :: sys
 integer, intent(out) :: Nc(3,3),ncentr
 character(len=*), intent(out) :: brav
 integer :: na,ax(3,maxop),iop(maxop),u(3),v(3),i1,j1,k1,ord,p(3)
 real(dp) :: lens(3),cv(3,3)
 select case (trim(sys))
 case ('cubic')
    call axes_of_order(nop,W,4,na,ax,iop)
    if (na.eq.0) call axes_of_order(nop,W,2,na,ax,iop)        ! T, Th: 2-fold axes along <100>
    if (na.ne.3) error stop 'internal error: cubic crystal without three 4- or 2-fold axes'
    do i1=1,3 ; Nc(i1,:)=ax(:,i1) ; enddo
    if (idet3(Nc).lt.0) Nc(3,:)=-Nc(3,:)
    ncentr=abs(idet3(Nc))
    select case (ncentr)
    case (1) ; brav='CUB'
    case (2) ; brav='BCC'
    case (4) ; brav='FCC'
    case default ; error stop 'internal error: cubic cell with unexpected centring'
    end select
 case ('tetragonal')
    call axes_of_order(nop,W,4,na,ax,iop)
    u=ax(:,1)
    call shortest_normal(Bn,u,v)
    Nc(1,:)=v
    Nc(2,:)=matmul(W(:,:,iop(1))*idet3(W(:,:,iop(1))),v)      ! a turned by the 4-fold axis
    Nc(3,:)=u
    if (idet3(Nc).lt.0) Nc(3,:)=-Nc(3,:)
    ncentr=abs(idet3(Nc))
    select case (ncentr)
    case (1) ; brav='TET'
    case (2) ; brav='BCT'
    case default ; error stop 'internal error: tetragonal cell with unexpected centring'
    end select
 case ('orthorhombic')
    call axes_of_order(nop,W,2,na,ax,iop)
    if (na.ne.3) error stop 'internal error: orthorhombic crystal without three 2-fold axes'
    do i1=1,3 ; Nc(i1,:)=ax(:,i1) ; enddo
    ncentr=abs(idet3(Nc))
    cv=matmul(dble(Nc),Bn)
    lens=[(sqrt(sum(cv(i1,:)**2)),i1=1,3)]
    select case (ncentr)
    case (1,4)
       brav='ORC' ; if (ncentr.eq.4) brav='ORCF'
       call order3(lens,p)
    case (2)
       if (is_lattice_point(Nc,[0.5d0,0.5d0,0.5d0])) then
          brav='ORCI' ; call order3(lens,p)
       else                                  ! base-centred: c normal to the centred face
          brav='ORCC'
          if (is_lattice_point(Nc,[0.5d0,0.5d0,0d0])) then
             k1=3
          elseif (is_lattice_point(Nc,[0d0,0.5d0,0.5d0])) then
             k1=1
          else
             k1=2
          endif
          i1=mod(k1,3)+1 ; j1=mod(k1+1,3)+1
          if (lens(i1).le.lens(j1)) then
             p=[i1,j1,k1]
          else
             p=[j1,i1,k1]
          endif
       endif
    case default ; error stop 'internal error: orthorhombic cell with unexpected centring'
    end select
    Nc=Nc(p,:)
    if (idet3(Nc).lt.0) Nc(3,:)=-Nc(3,:)
 case ('hexagonal','trigonal')
    ord=6 ; if (trim(sys).eq.'trigonal') ord=3
    call axes_of_order(nop,W,ord,na,ax,iop)
    u=ax(:,1)
    call shortest_normal(Bn,u,v)
    call axes_of_order(nop,W,3,na,ax,iop)    ! a 3-fold rotation about the same axis
    Nc(1,:)=v
    Nc(2,:)=matmul(W(:,:,iop(1))*idet3(W(:,:,iop(1))),v)
    Nc(3,:)=u
    cv=matmul(dble(Nc),Bn)
    if (dot_product(cross(cv(1,:),cv(2,:)),cv(3,:)).lt.0) then  ! turn by +120 degrees
       Nc(2,:)=matmul(W(:,:,iop(1))*idet3(W(:,:,iop(1))),Nc(2,:))
    endif
    cv=matmul(dble(Nc),Bn)
    if (dot_product(cross(cv(1,:),cv(2,:)),cv(3,:)).lt.0) error stop 'internal error: hexagonal axes'
    ncentr=abs(idet3(Nc))
    select case (ncentr)
    case (1) ; brav='HEX'
    case (3)
       brav='RHL'
       if (.not.is_lattice_point(Nc,[2d0/3,1d0/3,1d0/3])) then   ! reverse -> obverse setting
          Nc(1,:)=-Nc(1,:) ; Nc(2,:)=-Nc(2,:)
       endif
       if (.not.is_lattice_point(Nc,[2d0/3,1d0/3,1d0/3])) error stop 'internal error: rhombohedral centring'
    case default ; error stop 'internal error: hexagonal cell with unexpected centring'
    end select
 case ('monoclinic')
    call axes_of_order(nop,W,2,na,ax,iop)
    call monoclinic_cell(Bn,ax(:,1),Nc,ncentr)
    brav='MCL' ; if (ncentr.eq.2) brav='MCLC'
 case default
    call triclinic_cell(Bn,Nc)
    ncentr=1 ; brav='TRI'
 end select
 end subroutine conventional

 subroutine order3(l,p)
!! p: indices of l in increasing order (stable)
 real(dp), intent(in) :: l(3)
 integer, intent(out) :: p(3)
 integer :: i1,j1,t1
 p=[1,2,3]
 do i1=1,2 ; do j1=i1+1,3
    if (l(p(j1)).lt.l(p(i1))-1d-10) then
       t1=p(i1) ; p(i1)=p(j1) ; p(j1)=t1
    endif
 enddo ; enddo
 end subroutine order3

 subroutine monoclinic_cell(Bn,u,Nc,ncentr)
!! a along the 2-fold axis u; b, c a basis of the lattice vectors normal to a, with the
!! centring (if any) (a+b)/2, alpha < 90 degrees, |b| <= |c| where the centring allows it,
!! b and c as short as possible; right-handed
 real(dp), intent(in) :: Bn(3,3)
 integer, intent(in) :: u(3)
 integer, intent(out) :: Nc(3,3),ncentr
 integer :: n1,n2,n3,nv,vv(3,2000),i1,j1,best(3,3),trial(3,3),d0
 real(dp) :: cu(3),cw(3),l,sc,bsc,cb(3),ccv(3),cosal
 cu=matmul(dble(u),Bn)
 nv=0
 do n1=-4,4 ; do n2=-4,4 ; do n3=-4,4
    if (n1.eq.0 .and. n2.eq.0 .and. n3.eq.0) cycle
    cw=matmul(dble([n1,n2,n3]),Bn) ; l=sqrt(sum(cw**2))
    if (abs(dot_product(cw,cu)).gt.1d-6*l*sqrt(sum(cu**2))) cycle
    nv=nv+1 ; vv(:,nv)=[n1,n2,n3]
 enddo ; enddo ; enddo
 ! the index of the lattice in the cell (a, plane basis) is 1 or 2; find the plane basis
 ! of smallest |det| first: that is the index
 d0=huge(1)
 do i1=1,nv ; do j1=1,nv
    trial(1,:)=u ; trial(2,:)=vv(:,i1) ; trial(3,:)=vv(:,j1)
    if (idet3(trial).ne.0) d0=min(d0,abs(idet3(trial)))
 enddo ; enddo
 if (d0.ne.1 .and. d0.ne.2) error stop 'internal error: monoclinic cell with unexpected centring'
 ncentr=d0
 bsc=huge(1d0) ; best=0
 do i1=1,nv ; do j1=1,nv
    trial(1,:)=u ; trial(2,:)=vv(:,i1) ; trial(3,:)=vv(:,j1)
    if (abs(idet3(trial)).ne.ncentr) cycle
    if (ncentr.eq.2) then
       if (.not.is_lattice_point(trial,[0.5d0,0.5d0,0d0])) cycle
    endif
    cb=matmul(dble(vv(:,i1)),Bn) ; ccv=matmul(dble(vv(:,j1)),Bn)
    cosal=dot_product(cb,ccv)/sqrt(sum(cb**2)*sum(ccv**2))
    if (cosal.le.1d-8) cycle                                     ! alpha < 90 degrees
    ! as short as possible, then |b| <= |c|
    sc=sum(cb**2)+sum(ccv**2)
    if (sqrt(sum(cb**2)).gt.sqrt(sum(ccv**2))+1d-8) sc=sc*(1+1d-6)
    if (sc.lt.bsc-1d-8 .or. (abs(sc-bsc).le.1d-8 .and. lexi_gt(cb,matmul(dble(best(2,:)),Bn)))) then
       bsc=sc ; best=trial
    endif
 enddo ; enddo
 if (bsc.eq.huge(1d0)) error stop 'internal error: no monoclinic cell with alpha < 90 degrees'
 Nc=best
 if (idet3(Nc).lt.0) Nc(1,:)=-Nc(1,:)                            ! right-handed: turn a
 end subroutine monoclinic_cell

 subroutine triclinic_cell(Bn,Nc)
!! the reduced cell with the sign pattern (as pymatgen: identity, (-a,b,-c), (-a,-b,c),
!! (a,-b,-c), the last that fits) whose reciprocal angles are all acute or all obtuse
 real(dp), intent(in) :: Bn(3,3)
 integer, intent(out) :: Nc(3,3)
 integer :: s(3,4),is,i1,trial(3,3)
 real(dp) :: rc(3,3),ka(3)
 s(:,1)=[1,1,1] ; s(:,2)=[-1,1,-1] ; s(:,3)=[-1,-1,1] ; s(:,4)=[1,-1,-1]
 Nc=0 ; do i1=1,3 ; Nc(i1,i1)=1 ; enddo
 do is=1,4
    trial=0
    do i1=1,3 ; trial(i1,i1)=s(i1,is) ; enddo
    rc=recip(matmul(dble(trial),Bn))
    ka=[angle(rc(2,:),rc(3,:)),angle(rc(1,:),rc(3,:)),angle(rc(1,:),rc(2,:))]
    if (all(ka.le.90+1d-8) .or. all(ka.gt.90-1d-8)) Nc=trial
 enddo
 end subroutine triclinic_cell

!--------------------------------------------------------------- unique atom setting
 subroutine canonical_setting(Tm,lpar)
!! A crystal whose point group is smaller than that of its lattice (pyrite, wurtzite, CdI2 ...)
!! can be put in the same conventional cell in several inequivalent ways (e.g. wurtzite with
!! c or -c).  Of the proper rotations R that keep the conventional cell and the lattice
!! (det R = +1, so a chiral crystal is not mirrored), take the one that gives the smallest list
!! of (type, primitive fractional coordinates) -- sorted, from the atom of type 1 that makes
!! it smallest -- and apply it to the atoms, with that atom at the origin.  Two choices with
!! the same list describe the crystal identically; between them the list of positions in the
!! order of the input atoms decides, so that the output, atom by atom, does not depend on the
!! orientation, basis or origin of the input cell.
 real(dp), intent(in) :: Tm(3,3),lpar(6)
 real(dp) :: Gc(3,3),Cm(3,3),Ti(3,3),best(4,maxat),key(4,maxat),R(3,3),o(3),obest(3)
 real(dp) :: best2(4,maxat),key2(4,maxat)
 integer :: m(9),id,ii,i0,ibest(3,3),Ri(3,3)
 logical :: have
 call focex_frame(lpar,Cm)
 Gc=matmul(Cm,transpose(Cm)) ; Ti=inv3(Tm)
 have=.false. ; ibest=0
 do id=0,3**9-1
    ii=id
    do k=1,9
       m(k)=mod(ii,3)-1 ; ii=ii/3
    enddo
    Ri=reshape(m,[3,3])                            ! acts on conventional coordinates: x' = Ri x
    if (idet3(Ri).ne.1) cycle
    R=dble(Ri)
    if (maxval(abs(matmul(matmul(transpose(R),Gc),R)-Gc)).gt.1d-6*maxval(abs(Gc))) cycle
    ! the primitive vectors (rows of Tm, in conventional units) must go to lattice vectors
    if (maxval(abs(matmul(matmul(Tm,transpose(R)),Ti)-anint(matmul(matmul(Tm,transpose(R)),Ti)))).gt.1d-6) cycle
    call setting_key(R,Ti,key,key2,o)
    if (.not.have .or. better(key,key2,best,best2)) then
       best=key ; best2=key2 ; ibest=Ri ; obest=o ; have=.true.
    endif
 enddo
 do i0=1,nat
    xc(:,i0)=matmul(dble(ibest),xc(:,i0))
    fp(:,i0)=matmul(xc(:,i0),Ti)-obest
    fp(:,i0)=fp(:,i0)-floor(fp(:,i0)+keytol)
    where (fp(:,i0).ge.1d0-keytol) fp(:,i0)=0d0
    xc(:,i0)=matmul(fp(:,i0),Tm)
 enddo
 end subroutine canonical_setting

 subroutine setting_key(R,Ti,key,key2,obest)
!! key: sorted (type, primitive fractional coordinates) of the atoms turned by R, from the atom
!! of type 1 that makes it smallest; key2: the same positions in the input order (decides
!! between origins with the same key); obest: the position of that atom
 real(dp), intent(in) :: R(3,3),Ti(3,3)
 real(dp), intent(out) :: key(4,maxat),key2(4,maxat),obest(3)
 real(dp) :: p(3,maxat),trial(4,maxat),trial2(4,maxat),o(3)
 integer :: i0,j0
 logical :: first
 do i0=1,nat
    p(:,i0)=matmul(matmul(R,xc(:,i0)),Ti)
 enddo
 first=.true.
 do j0=1,nat
    if (typ(j0).ne.typ(1)) cycle
    o=p(:,j0)
    do i0=1,nat
       trial(1,i0)=typ(i0)
       trial(2:4,i0)=p(:,i0)-o
       trial(2:4,i0)=trial(2:4,i0)-floor(trial(2:4,i0)+keytol)      ! in [0,1), -1e-7 -> 0
       where (trial(2:4,i0).ge.1d0-keytol) trial(2:4,i0)=0d0
    enddo
    trial2=trial
    call sort_rows(trial)
    if (first) then
       key=trial ; key2=trial2 ; obest=o ; first=.false.
    elseif (better(trial,trial2,key,key2)) then
       key=trial ; key2=trial2 ; obest=o
    endif
 enddo
 end subroutine setting_key

 subroutine sort_rows(a)
!! sort the columns 1..nat of a(4,:) in lexicographic order
 real(dp), intent(inout) :: a(4,maxat)
 real(dp) :: tmp(4)
 integer :: i0,j0
 do i0=2,nat
    tmp=a(:,i0) ; j0=i0-1
    do while (j0.ge.1)
       if (.not.col_less(tmp,a(:,j0))) exit
       a(:,j0+1)=a(:,j0) ; j0=j0-1
    enddo
    a(:,j0+1)=tmp
 enddo
 end subroutine sort_rows

 logical function col_less(x1,x2)
!! lexicographic x1 < x2, values within keytol (1e-5, fractional) counting as equal, so that
!! noise in the input positions cannot decide the setting
 real(dp), intent(in) :: x1(:),x2(:)
 integer :: m1
 col_less=.false.
 do m1=1,size(x1)
    if (x1(m1).lt.x2(m1)-keytol) then
       col_less=.true. ; return
    elseif (x1(m1).gt.x2(m1)+keytol) then
       return
    endif
 enddo
 end function col_less

 logical function better(k1,k1b,k2,k2b)
!! (k1, then k1b) < (k2, then k2b)
 real(dp), intent(in) :: k1(4,maxat),k1b(4,maxat),k2(4,maxat),k2b(4,maxat)
 better = key_less(k1,k2) .or. (.not.key_less(k2,k1) .and. key_less(k1b,k2b))
 end function better

 logical function key_less(k1,k2)
 real(dp), intent(in) :: k1(4,maxat),k2(4,maxat)
 integer :: i0
 key_less=.false.
 do i0=1,nat
    if (col_less(k1(:,i0),k2(:,i0))) then
       key_less=.true. ; return
    elseif (col_less(k2(:,i0),k1(:,i0))) then
       return
    endif
 enddo
 end function key_less

!--------------------------------------------------------------- Setyawan-Curtarolo
 subroutine sc_primitive(brav,Tm)
!! primitive vectors (rows) in units of the conventional ones
 character(len=*), intent(in) :: brav
 real(dp), intent(out) :: Tm(3,3)
 Tm=eye3()
 select case (trim(brav))
 case ('FCC','ORCF')
    Tm=reshape([0d0,.5d0,.5d0, .5d0,0d0,.5d0, .5d0,.5d0,0d0],[3,3],order=[2,1])
 case ('BCC','BCT','ORCI')
    Tm=reshape([-.5d0,.5d0,.5d0, .5d0,-.5d0,.5d0, .5d0,.5d0,-.5d0],[3,3],order=[2,1])
 case ('ORCC')
    Tm=reshape([.5d0,-.5d0,0d0, .5d0,.5d0,0d0, 0d0,0d0,1d0],[3,3],order=[2,1])
 case ('MCLC')
    Tm=reshape([.5d0,.5d0,0d0, -.5d0,.5d0,0d0, 0d0,0d0,1d0],[3,3],order=[2,1])
 case ('RHL')
    ! the lattice vectors (2,1,1)/3, (-1,1,1)/3, (-1,-2,1)/3 of the obverse hexagonal cell, in
    ! this (right-handed) order: they are then the Setyawan-Curtarolo vectors, turned.  pymatgen
    ! takes them as (-1,1,1),(2,1,1),(-1,-2,1)/3, a left-handed set, and keeps the fractional
    ! coordinates when it rebuilds the right-handed S-C cell, i.e. it mirrors the crystal
    Tm=reshape([2d0,1d0,1d0, -1d0,1d0,1d0, -1d0,-2d0,1d0],[3,3],order=[2,1])/3
 end select
 end subroutine sc_primitive

 subroutine rhl_check(Pp)
!! the rhombohedral primitive vectors must be the right-handed set of equal lengths and angles
!! that Setyawan-Curtarolo (and pymatgen) build from a and alpha
 real(dp), intent(in) :: Pp(3,3)
 real(dp) :: lp(6)
 call lattice_params(Pp,lp)
 if (abs(lp(1)-lp(2)).gt.1d-6*lp(1) .or. abs(lp(1)-lp(3)).gt.1d-6*lp(1) .or.  &
     abs(lp(4)-lp(5)).gt.1d-6 .or. abs(lp(4)-lp(6)).gt.1d-6 .or. det3(Pp).le.0) &
     error stop 'internal error: rhombohedral primitive cell'
 end subroutine rhl_check

 subroutine sc_kpath(brav,Cc,Pp,name)
!! sub-case, high-symmetry points (primitive reduced coordinates: nk, kname, kred) and path
!! (npath, path), as in pymatgen's KPathSetyawanCurtarolo
 character(len=*), intent(in) :: brav
 real(dp), intent(in) :: Cc(3,3),Pp(3,3)
 character(len=*), intent(out) :: name
 real(dp) :: lc(6),lp(6),a,b,c,al,eta,nu,zeta,phi,delta,mu,psi,rho,omega,kg,rc(3,3),e,s
 call lattice_params(Cc,lc) ; call lattice_params(Pp,lp)
 a=lc(1) ; b=lc(2) ; c=lc(3) ; al=lc(4)*pi/180
 e=1d-6
 nk=0 ; npath=0
 select case (trim(brav))
 case ('CUB')
    name='CUB'
    call kp('G',0,0,0) ; call kp('X',0d0,.5d0,0d0) ; call kp('R',.5d0,.5d0,.5d0) ; call kp('M',.5d0,.5d0,0d0)
    call pa('G X M G R X|M R')
 case ('FCC')
    name='FCC'
    call kp('G',0,0,0) ; call kp('K',3d0/8,3d0/8,3d0/4) ; call kp('L',.5d0,.5d0,.5d0)
    call kp('U',5d0/8,.25d0,5d0/8) ; call kp('W',.5d0,.25d0,.75d0) ; call kp('X',.5d0,0d0,.5d0)
    call pa('G X W K G L U W L K|U X')
 case ('BCC')
    name='BCC'
    call kp('G',0,0,0) ; call kp('H',.5d0,-.5d0,.5d0) ; call kp('P',.25d0,.25d0,.25d0) ; call kp('N',0d0,0d0,.5d0)
    call pa('G H N G P H|P N')
 case ('TET')
    name='TET'
    call kp('G',0,0,0) ; call kp('A',.5d0,.5d0,.5d0) ; call kp('M',.5d0,.5d0,0d0)
    call kp('R',0d0,.5d0,.5d0) ; call kp('X',0d0,.5d0,0d0) ; call kp('Z',0d0,0d0,.5d0)
    call pa('G X M G Z R A Z|X R|M A')
 case ('BCT')
    if (c.lt.a) then
       name='BCT1' ; eta=(1+c**2/a**2)/4
       call kp('G',0,0,0) ; call kp('M',-.5d0,.5d0,.5d0) ; call kp('N',0d0,.5d0,0d0)
       call kp('P',.25d0,.25d0,.25d0) ; call kp('X',0d0,0d0,.5d0) ; call kp('Z',eta,eta,-eta)
       call kp('Z1',-eta,1-eta,eta)
       call pa('G X M G Z P N Z1 M|X P')
    else
       name='BCT2' ; eta=(1+a**2/c**2)/4 ; zeta=a**2/(2*c**2)
       call kp('G',0,0,0) ; call kp('N',0d0,.5d0,0d0) ; call kp('P',.25d0,.25d0,.25d0)
       call kp('Sg',-eta,eta,eta) ; call kp('Sg1',eta,1-eta,-eta) ; call kp('X',0d0,0d0,.5d0)
       call kp('Y',-zeta,zeta,.5d0) ; call kp('Y1',.5d0,.5d0,-zeta) ; call kp('Z',.5d0,.5d0,-.5d0)
       call pa('G X Y Sg G Z Sg1 N P Y1 Z|X P')
    endif
 case ('ORC')
    name='ORC'
    call kp('G',0,0,0) ; call kp('R',.5d0,.5d0,.5d0) ; call kp('S',.5d0,.5d0,0d0) ; call kp('T',0d0,.5d0,.5d0)
    call kp('U',.5d0,0d0,.5d0) ; call kp('X',.5d0,0d0,0d0) ; call kp('Y',0d0,.5d0,0d0) ; call kp('Z',0d0,0d0,.5d0)
    call pa('G X S Y G Z U R T Z|Y T|U X|S R')
 case ('ORCF')
    s=1/a**2-1/b**2-1/c**2
    if (abs(s).le.e/a**2) then
       name='ORCF3'
    elseif (s.gt.0) then
       name='ORCF1'
    else
       name='ORCF2'
    endif
    if (name.eq.'ORCF2') then
       phi=(1+c**2/b**2-c**2/a**2)/4 ; eta=(1+a**2/b**2-a**2/c**2)/4 ; delta=(1+b**2/a**2-b**2/c**2)/4
       call kp('G',0,0,0) ; call kp('C',.5d0,.5d0-eta,1-eta) ; call kp('C1',.5d0,.5d0+eta,eta)
       call kp('D',.5d0-delta,.5d0,1-delta) ; call kp('D1',.5d0+delta,.5d0,delta) ; call kp('L',.5d0,.5d0,.5d0)
       call kp('H',1-phi,.5d0-phi,.5d0) ; call kp('H1',phi,.5d0+phi,.5d0) ; call kp('X',0d0,.5d0,.5d0)
       call kp('Y',.5d0,0d0,.5d0) ; call kp('Z',.5d0,.5d0,0d0)
       call pa('G Y C D X G Z D1 H C|C1 Z|X H1|H Y|L G')
    else
       zeta=(1+a**2/b**2-a**2/c**2)/4 ; eta=(1+a**2/b**2+a**2/c**2)/4
       call kp('G',0,0,0) ; call kp('A',.5d0,.5d0+zeta,zeta) ; call kp('A1',.5d0,.5d0-zeta,1-zeta)
       call kp('L',.5d0,.5d0,.5d0) ; call kp('T',1d0,.5d0,.5d0) ; call kp('X',0d0,eta,eta)
       call kp('X1',1d0,1-eta,1-eta) ; call kp('Y',.5d0,0d0,.5d0) ; call kp('Z',.5d0,.5d0,0d0)
       if (name.eq.'ORCF1') then
          call pa('G Y T Z G X A1 Y|T X1|X A Z|L G')
       else
          call pa('G Y T Z G X A1 Y|X A Z|L G')
       endif
    endif
 case ('ORCI')
    name='ORCI'
    zeta=(1+a**2/c**2)/4 ; eta=(1+b**2/c**2)/4 ; delta=(b**2-a**2)/(4*c**2) ; mu=(a**2+b**2)/(4*c**2)
    call kp('G',0,0,0) ; call kp('L',-mu,mu,.5d0-delta) ; call kp('L1',mu,-mu,.5d0+delta)
    call kp('L2',.5d0-delta,.5d0+delta,-mu) ; call kp('R',0d0,.5d0,0d0) ; call kp('S',.5d0,0d0,0d0)
    call kp('T',0d0,0d0,.5d0) ; call kp('W',.25d0,.25d0,.25d0) ; call kp('X',-zeta,zeta,zeta)
    call kp('X1',zeta,1-zeta,-zeta) ; call kp('Y',eta,-eta,eta) ; call kp('Y1',1-eta,eta,-eta)
    call kp('Z',.5d0,.5d0,-.5d0)
    call pa('G X L T W R X1 Z G Y S W|L1 Y|Y1 Z')
 case ('ORCC')
    name='ORCC' ; zeta=(1+a**2/b**2)/4
    call kp('G',0,0,0) ; call kp('A',zeta,zeta,.5d0) ; call kp('A1',-zeta,1-zeta,.5d0)
    call kp('R',0d0,.5d0,.5d0) ; call kp('S',0d0,.5d0,0d0) ; call kp('T',-.5d0,.5d0,.5d0)
    call kp('X',zeta,zeta,0d0) ; call kp('X1',-zeta,1-zeta,0d0) ; call kp('Y',-.5d0,.5d0,0d0)
    call kp('Z',0d0,0d0,.5d0)
    call pa('G X S R A Z G Y X1 A1 T Y|Z T')
 case ('HEX')
    name='HEX'
    call kp('G',0,0,0) ; call kp('A',0d0,0d0,.5d0) ; call kp('H',1d0/3,1d0/3,.5d0)
    call kp('K',1d0/3,1d0/3,0d0) ; call kp('L',.5d0,0d0,.5d0) ; call kp('M',.5d0,0d0,0d0)
    call pa('G M K G A L H A|L M|K H')
 case ('RHL')
    al=lp(4)*pi/180
    if (lp(4).lt.90) then
       name='RHL1' ; eta=(1+4*cos(al))/(2+4*cos(al)) ; nu=0.75d0-eta/2
       call kp('G',0,0,0) ; call kp('B',eta,.5d0,1-eta) ; call kp('B1',.5d0,1-eta,eta-1)
       call kp('F',.5d0,.5d0,0d0) ; call kp('L',.5d0,0d0,0d0) ; call kp('L1',0d0,0d0,-.5d0)
       call kp('P',eta,nu,nu) ; call kp('P1',1-nu,1-nu,1-eta) ; call kp('P2',nu,nu,eta-1)
       call kp('Q',1-nu,nu,0d0) ; call kp('X',nu,0d0,-nu) ; call kp('Z',.5d0,.5d0,.5d0)
       call pa('G L B1|B Z G X|Q F P1 Z|L P')
    else
       name='RHL2' ; eta=1/(2*tan(al/2)**2) ; nu=0.75d0-eta/2
       call kp('G',0,0,0) ; call kp('F',.5d0,-.5d0,0d0) ; call kp('L',.5d0,0d0,0d0)
       call kp('P',1-nu,-nu,1-nu) ; call kp('P1',nu,nu-1,nu-1) ; call kp('Q',eta,eta,eta)
       call kp('Q1',1-eta,-eta,-eta) ; call kp('Z',.5d0,-.5d0,.5d0)
       call pa('G P Z Q G F P1 Q1 L Z')
    endif
 case ('MCL')
    name='MCL'
    eta=(1-b*cos(al)/c)/(2*sin(al)**2) ; nu=0.5d0-eta*c*cos(al)/b
    call kp('G',0,0,0) ; call kp('A',.5d0,.5d0,0d0) ; call kp('C',0d0,.5d0,.5d0) ; call kp('D',.5d0,0d0,.5d0)
    call kp('D1',.5d0,0d0,-.5d0) ; call kp('E',.5d0,.5d0,.5d0) ; call kp('H',0d0,eta,1-nu)
    call kp('H1',0d0,1-eta,nu) ; call kp('H2',0d0,eta,-nu) ; call kp('M',.5d0,eta,1-nu)
    call kp('M1',.5d0,1-eta,nu) ; call kp('M2',.5d0,eta,-nu) ; call kp('X',0d0,.5d0,0d0)
    call kp('Y',0d0,0d0,.5d0) ; call kp('Y1',0d0,0d0,-.5d0) ; call kp('Z',.5d0,0d0,0d0)
    call pa('G Y H C E M1 A X H1|M D Z|Y D')
 case ('MCLC')
    rc=recip(Pp)
    kg=angle(rc(1,:),rc(2,:))
    s=b*cos(al)/c+b**2*sin(al)**2/a**2
    if (kg.gt.90+e) then
       name='MCLC1'
    elseif (abs(kg-90).le.e) then
       name='MCLC2'
    elseif (s.lt.1-e) then
       name='MCLC3'
    elseif (s.gt.1+e) then
       name='MCLC5'
    else
       name='MCLC4'
    endif
    select case (trim(name))
    case ('MCLC1','MCLC2')
       zeta=(2-b*cos(al)/c)/(4*sin(al)**2) ; eta=0.5d0+2*zeta*c*cos(al)/b
       psi=0.75d0-a**2/(4*b**2*sin(al)**2) ; phi=psi+(0.75d0-psi)*b*cos(al)/c
       call kp('G',0,0,0) ; call kp('N',.5d0,0d0,0d0) ; call kp('N1',0d0,-.5d0,0d0)
       call kp('F',1-zeta,1-zeta,1-eta) ; call kp('F1',zeta,zeta,eta) ; call kp('F2',-zeta,-zeta,1-eta)
       if (name.eq.'MCLC2') call kp('F3',1-zeta,-zeta,1-eta)
       call kp('I',phi,1-phi,.5d0) ; call kp('I1',1-phi,phi-1,.5d0) ; call kp('L',.5d0,.5d0,.5d0)
       call kp('M',.5d0,0d0,.5d0) ; call kp('X',1-psi,psi-1,0d0) ; call kp('X1',psi,1-psi,0d0)
       call kp('X2',psi-1,-psi,0d0) ; call kp('Y',.5d0,.5d0,0d0) ; call kp('Y1',-.5d0,-.5d0,0d0)
       call kp('Z',0d0,0d0,.5d0)
       if (name.eq.'MCLC1') then
          call pa('G Y F L I|I1 Z F1|Y X1|X G N|M G')
       else
          call pa('G Y F L I|I1 Z F1|N G M')
       endif
    case ('MCLC3','MCLC4')
       mu=(1+b**2/a**2)/4 ; delta=b*c*cos(al)/(2*a**2)
       zeta=mu-0.25d0+(1-b*cos(al)/c)/(4*sin(al)**2) ; eta=0.5d0+2*zeta*c*cos(al)/b
       phi=1+zeta-2*mu ; psi=eta-2*delta
       call kp('G',0,0,0) ; call kp('F',1-phi,1-phi,1-psi) ; call kp('F1',phi,phi-1,psi)
       call kp('F2',1-phi,-phi,1-psi) ; call kp('H',zeta,zeta,eta) ; call kp('H1',1-zeta,-zeta,1-eta)
       call kp('H2',-zeta,-zeta,1-eta) ; call kp('I',.5d0,-.5d0,.5d0) ; call kp('M',.5d0,0d0,.5d0)
       call kp('N',.5d0,0d0,0d0) ; call kp('N1',0d0,-.5d0,0d0) ; call kp('X',.5d0,-.5d0,0d0)
       call kp('Y',mu,mu,delta) ; call kp('Y1',1-mu,-mu,-delta) ; call kp('Y2',-mu,-mu,-delta)
       call kp('Y3',mu,mu-1,delta) ; call kp('Z',0d0,0d0,.5d0)
       if (name.eq.'MCLC3') then
          call pa('G Y F H Z I F1|H1 Y1 X G N|M G')
       else
          call pa('G Y F H Z I|H1 Y1 X G N|M G')
       endif
    case ('MCLC5')
       zeta=(b**2/a**2+(1-b*cos(al)/c)/sin(al)**2)/4 ; eta=0.5d0+2*zeta*c*cos(al)/b
       mu=eta/2+b**2/(4*a**2)-b*c*cos(al)/(2*a**2) ; nu=2*mu-zeta ; rho=1-zeta*a**2/b**2
       omega=(4*nu-1-b**2*sin(al)**2/a**2)*c/(2*b*cos(al)) ; delta=zeta*c*cos(al)/b+omega/2-0.25d0
       call kp('G',0,0,0) ; call kp('F',nu,nu,omega) ; call kp('F1',1-nu,1-nu,1-omega)
       call kp('F2',nu,nu-1,omega) ; call kp('H',zeta,zeta,eta) ; call kp('H1',1-zeta,-zeta,1-eta)
       call kp('H2',-zeta,-zeta,1-eta) ; call kp('I',rho,1-rho,.5d0) ; call kp('I1',1-rho,rho-1,.5d0)
       call kp('L',.5d0,.5d0,.5d0) ; call kp('M',.5d0,0d0,.5d0) ; call kp('N',.5d0,0d0,0d0)
       call kp('N1',0d0,-.5d0,0d0) ; call kp('X',.5d0,-.5d0,0d0) ; call kp('Y',mu,mu,delta)
       call kp('Y1',1-mu,-mu,-delta) ; call kp('Y2',-mu,-mu,-delta) ; call kp('Y3',mu,mu-1,delta)
       call kp('Z',0d0,0d0,.5d0)
       call pa('G Y F L I|I1 Z H F1|H1 Y1 X G N|M G')
    end select
 case ('TRI')
    rc=recip(Pp)
    lp=[0d0,0d0,0d0,angle(rc(2,:),rc(3,:)),angle(rc(1,:),rc(3,:)),angle(rc(1,:),rc(2,:))]
    if (lp(4).gt.90 .and. lp(5).gt.90 .and. lp(6).ge.90-e) then
       name='TRI1a'
       call kp('G',0,0,0) ; call kp('L',.5d0,.5d0,0d0) ; call kp('M',0d0,.5d0,.5d0) ; call kp('N',.5d0,0d0,.5d0)
       call kp('R',.5d0,.5d0,.5d0) ; call kp('X',.5d0,0d0,0d0) ; call kp('Y',0d0,.5d0,0d0) ; call kp('Z',0d0,0d0,.5d0)
    else
       name='TRI1b'
       call kp('G',0,0,0) ; call kp('L',.5d0,-.5d0,0d0) ; call kp('M',0d0,0d0,.5d0) ; call kp('N',-.5d0,-.5d0,.5d0)
       call kp('R',0d0,-.5d0,.5d0) ; call kp('X',0d0,-.5d0,0d0) ; call kp('Y',.5d0,0d0,0d0) ; call kp('Z',-.5d0,0d0,.5d0)
    endif
    call pa('X G Y|L G Z|N G M|R G')
 case default
    error stop 'internal error: unknown Bravais lattice'
 end select
 end subroutine sc_kpath

 subroutine kp(l,x1,x2,x3)
 character(len=*), intent(in) :: l
 class(*), intent(in) :: x1,x2,x3
 nk=nk+1 ; kname(nk)=l ; kred(:,nk)=[num(x1),num(x2),num(x3)]
 end subroutine kp

 real(dp) function num(x1)
 class(*), intent(in) :: x1
 select type (x1)
 type is (integer) ; num=dble(x1)
 type is (real(dp)) ; num=x1
 class default ; num=0
 end select
 end function num

 subroutine pa(s)
!! the path, segments separated by '|', laid out as one list (consecutive repeats removed), as
!! focex-agent.py does for kpbs.params
 character(len=*), intent(in) :: s
 character(len=4) :: w(60)
 integer :: nw,iw
 character(len=len(s)) :: s2
 s2=s
 do iw=1,len(s2)
    if (s2(iw:iw).eq.'|') s2(iw:iw)=' '
 enddo
 call split_words(s2,nw,w)
 npath=0
 do iw=1,nw
    if (npath.gt.0) then
       if (trim(path(npath)).eq.trim(w(iw))) cycle
    endif
    npath=npath+1 ; path(npath)=w(iw)
 enddo
 end subroutine pa

!--------------------------------------------------------------- output
 subroutine write_structure(fn,lpar,Tm)
 character(len=*), intent(in) :: fn
 real(dp), intent(in) :: lpar(6),Tm(3,3)
 integer :: u,i0,r
 character(len=400) :: s
 character(len=12) :: rr
 open(newunit=u,file=trim(fn),status='replace')
 write(u,'(6(f15.10,1x),a)') lpar,' # a b c alpha beta gamma of the conventional cell'
 s=''
 do i0=1,3 ; do r=1,3
    s=trim(s)//' '//trim(gfmt(Tm(i0,r)))
 enddo ; enddo
 write(u,'(a)') trim(adjustl(s))//'  # primitive vectors R01 R02 R03 in units of the conventional vectors'
 write(u,'(a)') '1.0 0              # scale for line 1; 0 = coordinates in units of the conventional vectors'
 write(u,'(a)') '1 1 1 0 0 0 0 0    # default, check: include rank r=1..8 in the fit (0/1)'
 write(u,'(a)') '1 1 0 0            # default, check: invariances: translational, rotational, Huang, enforcement'
 write(u,'(a)') '0 300              # default, check: itemp (0 = no Boltzmann weighting), tempk (K)'
 write(u,'(a)') '1 .false.          # default, check: number of POSCARi/FORCEDISPi pairs, verbosity'
 write(u,'(i0,a)') ntyp,'                  # number of atom types'
 s=''
 do i0=1,ntyp
    write(rr,'(f12.4)') mass_of(species(i0))
    s=trim(s)//' '//trim(adjustl(rr))
 enddo
 write(u,'(a)') trim(adjustl(s))//'  # mass of each atom type (amu)'
 s=''
 do i0=1,ntyp ; s=trim(s)//' '//trim(species(i0)) ; enddo
 write(u,'(a)') trim(adjustl(s))//'  # name of each atom type'
 write(u,'(i0,a)') nat,'                  # number of atoms in the primitive cell'
 write(u,'(a)') '0                  # default, check: fc2range (0 = FC2 range from the supercell, 1 = next line)'
 write(u,'(a)') trim(repeat_int(11))//'  # default, check: shells for rank 2 (one per atom)'
 write(u,'(a)') trim(repeat_int(3))//'  # default, check: shells for rank 3'
 do r=4,8
    write(u,'(a,i0)') trim(repeat_int(0))//'  # default, check: shells for rank ',r
 enddo
 do i0=1,nat
    if (i0.eq.1) then
       write(u,'(i0,1x,i0,3(1x,f15.10),a)') i0,typ(i0),clean(xc(:,i0)), &
            '  # atom number, atom type, coordinates in units of CONVENTIONAL vectors'
    else
       write(u,'(i0,1x,i0,3(1x,f15.10))') i0,typ(i0),clean(xc(:,i0))
    endif
 enddo
 close(u)
 end subroutine write_structure

 function repeat_int(v) result(s)
 integer, intent(in) :: v
 character(len=8000) :: s
 character(len=12) :: w
 integer :: i0
 s=''
 write(w,'(i0)') v
 do i0=1,nat ; s=trim(s)//' '//trim(w) ; enddo
 s=adjustl(s)
 end function repeat_int

 subroutine write_poscar(fn,Pm)
 character(len=*), intent(in) :: fn
 real(dp), intent(in) :: Pm(3,3)
 integer :: u,i0
 character(len=400) :: s
 open(newunit=u,file=trim(fn),status='replace')
 write(u,'(a)') trim(title)//' (Setyawan-Curtarolo primitive cell, FOCEX frame; poscar2focex)'
 write(u,'(a)') '1.0'
 do i0=1,3
    write(u,'(3(2x,f16.10))') clean(Pm(i0,:))
 enddo
 s='' ; do i0=1,ntyp ; s=trim(s)//'  '//trim(species(i0)) ; enddo
 write(u,'(a)') trim(s)
 s='' ; do i0=1,ntyp ; write(s(len_trim(s)+3:),'(i0)') counts(i0) ; enddo
 write(u,'(a)') trim(s)
 write(u,'(a)') 'Direct'
 do i0=1,nat
    write(u,'(3(2x,f16.10))') clean(fp(:,i0))
 enddo
 close(u)
 end subroutine write_poscar

 subroutine write_kpbs(fn)
 character(len=*), intent(in) :: fn
 integer :: u,i0,j0
 open(newunit=u,file=trim(fn),status='replace')
 write(u,'(a)') '1     # units of the k-points: 0 conventional, 1 primitive reduced coordinates'
 write(u,'(a)') '40    # number of k-points along each direction'
 write(u,'(i0,a)') npath-1,'     # number of directions; Setyawan-Curtarolo '//trim(scname)//' path'
 do i0=1,npath
    do j0=1,nk
       if (trim(kname(j0)).eq.trim(path(i0))) exit
    enddo
    write(u,'(a4,3(1x,f12.8))') kname(j0),clean(kred(:,j0))
 enddo
 close(u)
 end subroutine write_kpbs

 subroutine report()
 real(dp) :: lp(6)
 write(*,'(a)') trim(fin)//': '//trim(title)
 write(*,'(a,i0,a,i0,a,i0,a)') '  ',nat,' atoms; ',nlat,' lattice operations, ',ncry,' crystal operations'
 write(*,'(a)') '  crystal system '//trim(system)//', Setyawan-Curtarolo lattice '//trim(bravais)// &
      ' ('//trim(scname)//')'
 call lattice_params(A,lp) ; write(*,'(a,3f11.5,a,3f9.3)') '  input cell        :',lp(1:3),' A  ',lp(4:6)
 call lattice_params(B,lp) ; write(*,'(a,3f11.5,a,3f9.3)') '  Niggli cell       :',lp(1:3),' A  ',lp(4:6)
 call lattice_params(Pc,lp) ; write(*,'(a,3f11.5,a,3f9.3)') '  primitive (S-C)   :',lp(1:3),' A  ',lp(4:6)
 write(*,'(a,3f11.5,a,3f9.3)') '  conventional (S-C):',latp(1:3),' A  ',latp(4:6)
 write(*,'(a)') '  wrote '//trim(fout)//', '//trim(fpos)//' (build the DFT supercells from it) and '//trim(fkp)
 end subroutine report

 function gfmt(x) result(s)
!! a number with at most 12 decimals, trailing zeros removed (like %.12g for these values)
 real(dp), intent(in) :: x
 character(len=24) :: s
 integer :: l
 write(s,'(f24.12)') clean1(x)
 s=adjustl(s)
 l=len_trim(s)
 do while (s(l:l).eq.'0')
    s(l:l)=' ' ; l=l-1
 enddo
 if (s(l:l).eq.'.') s(l:l)=' '
 if (trim(s).eq.'-0') s='0'
 end function gfmt

 real(dp) function clean1(x)
 real(dp), intent(in) :: x
 clean1=anint(x*1d10)/1d10
 if (abs(clean1).lt.1d-12) clean1=0
 end function clean1

 function clean(v) result(w)
 real(dp), intent(in) :: v(:)
 real(dp) :: w(3)
 integer :: i0
 do i0=1,3 ; w(i0)=clean1(v(i0)) ; enddo
 end function clean

 real(dp) function mass_of(el)
 character(len=*), intent(in) :: el
 character(len=2), parameter :: sym(103)=[character(len=2) :: 'H','He','Li','Be','B','C','N','O','F', &
  'Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu',   &
  'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',     &
  'In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho',     &
  'Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',     &
  'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']
 real(dp), parameter :: ms(103)=[1.0079d0,4.0026d0,6.9410d0,9.0122d0,10.8110d0,12.0107d0,14.0067d0, &
  15.9994d0,18.9984d0,20.1797d0,22.9898d0,24.3050d0,26.9815d0,28.0855d0,30.9738d0,32.0650d0,         &
  35.4530d0,39.9480d0,39.0983d0,40.0780d0,44.9559d0,47.8670d0,50.9415d0,51.9961d0,54.9380d0,         &
  55.8450d0,58.9332d0,58.6934d0,63.5460d0,65.4090d0,69.7230d0,72.6400d0,74.9216d0,78.9600d0,         &
  79.9040d0,83.7980d0,85.4678d0,87.6200d0,88.9059d0,91.2240d0,92.9064d0,95.9400d0,98.0000d0,         &
  101.0700d0,102.9055d0,106.4200d0,107.8682d0,112.4110d0,114.8180d0,118.7100d0,121.7600d0,           &
  127.6000d0,126.9045d0,131.2930d0,132.9055d0,137.3270d0,138.9055d0,140.1160d0,140.9076d0,           &
  144.2420d0,145.0000d0,150.3600d0,151.9640d0,157.2500d0,158.9254d0,162.5000d0,164.9303d0,           &
  167.2590d0,168.9342d0,173.0400d0,174.9670d0,178.4900d0,180.9479d0,183.8400d0,186.2070d0,           &
  190.2300d0,192.2170d0,195.0840d0,196.9666d0,200.5900d0,204.3833d0,207.2000d0,208.9804d0,           &
  210.0000d0,210.0000d0,220.0000d0,223.0000d0,226.0000d0,227.0000d0,232.0381d0,231.0359d0,           &
  238.0289d0,237.0000d0,244.0000d0,243.0000d0,247.0000d0,247.0000d0,251.0000d0,252.0000d0,           &
  257.0000d0,258.0000d0,259.0000d0,262.0000d0]
 character(len=8) :: e0
 integer :: i0,p
 e0=el ; p=index(e0,'_') ; if (p.gt.0) e0=e0(1:p-1)
 do i0=1,103
    if (trim(sym(i0)).eq.trim(e0)) then
       mass_of=ms(i0) ; return
    endif
 enddo
 write(*,*) 'unknown element ',trim(el),' on line 6 of the POSCAR'
 error stop 1
 end function mass_of

!--------------------------------------------------------------- small linear algebra
 function eye3() result(E)
 real(dp) :: E(3,3)
 E=0 ; E(1,1)=1 ; E(2,2)=1 ; E(3,3)=1
 end function eye3

 real(dp) function det3(M)
 real(dp), intent(in) :: M(3,3)
 det3=M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))-M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1)) &
     +M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
 end function det3

 integer function idet3(M)
 integer, intent(in) :: M(3,3)
 idet3=M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))-M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1)) &
      +M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
 end function idet3

 function inv3(M) result(R)
 real(dp), intent(in) :: M(3,3)
 real(dp) :: R(3,3),d
 d=det3(M)
 R(1,1)=M(2,2)*M(3,3)-M(2,3)*M(3,2) ; R(1,2)=M(1,3)*M(3,2)-M(1,2)*M(3,3) ; R(1,3)=M(1,2)*M(2,3)-M(1,3)*M(2,2)
 R(2,1)=M(2,3)*M(3,1)-M(2,1)*M(3,3) ; R(2,2)=M(1,1)*M(3,3)-M(1,3)*M(3,1) ; R(2,3)=M(1,3)*M(2,1)-M(1,1)*M(2,3)
 R(3,1)=M(2,1)*M(3,2)-M(2,2)*M(3,1) ; R(3,2)=M(1,2)*M(3,1)-M(1,1)*M(3,2) ; R(3,3)=M(1,1)*M(2,2)-M(1,2)*M(2,1)
 R=R/d
 end function inv3

 function cross(a1,b1) result(c1)
 real(dp), intent(in) :: a1(:),b1(:)
 real(dp) :: c1(3)
 c1=[a1(2)*b1(3)-a1(3)*b1(2), a1(3)*b1(1)-a1(1)*b1(3), a1(1)*b1(2)-a1(2)*b1(1)]
 end function cross

 real(dp) function angle(a1,b1)
 real(dp), intent(in) :: a1(:),b1(:)
 angle=acos(max(-1d0,min(1d0,dot_product(a1,b1)/sqrt(sum(a1**2)*sum(b1**2)))))*180/pi
 end function angle

 function recip(M) result(R)
!! reciprocal vectors (rows, with 2 pi) of the cell M (rows)
 real(dp), intent(in) :: M(3,3)
 real(dp) :: R(3,3)
 R=2*pi*transpose(inv3(M))
 end function recip

 integer function igcd(a1,b1)
 integer, intent(in) :: a1,b1
 integer :: x1,y1,t1
 x1=a1 ; y1=b1
 do while (y1.ne.0)
    t1=mod(x1,y1) ; x1=y1 ; y1=t1
 enddo
 igcd=x1                                    ! gcd(0,0) = 0
 end function igcd

 subroutine lattice_params(M,lp)
 real(dp), intent(in) :: M(3,3)
 real(dp), intent(out) :: lp(6)
 lp(1)=sqrt(sum(M(1,:)**2)) ; lp(2)=sqrt(sum(M(2,:)**2)) ; lp(3)=sqrt(sum(M(3,:)**2))
 lp(4)=angle(M(2,:),M(3,:)) ; lp(5)=angle(M(1,:),M(3,:)) ; lp(6)=angle(M(1,:),M(2,:))
 end subroutine lattice_params

 subroutine focex_frame(lp,Fm)
!! the conventional vectors (rows) as FOCEX and SC-SNAPS build them: a along x, b in xy
 real(dp), intent(in) :: lp(6)
 real(dp), intent(out) :: Fm(3,3)
 real(dp) :: al,be,ga,cx,cy
 al=lp(4)*pi/180 ; be=lp(5)*pi/180 ; ga=lp(6)*pi/180
 cx=lp(3)*cos(be) ; cy=lp(3)*(cos(al)-cos(ga)*cos(be))/sin(ga)
 Fm(1,:)=[lp(1),0d0,0d0]
 Fm(2,:)=[lp(2)*cos(ga),lp(2)*sin(ga),0d0]
 Fm(3,:)=[cx,cy,sqrt(lp(3)**2-cx**2-cy**2)]
 end subroutine focex_frame

 end program poscar2focex
