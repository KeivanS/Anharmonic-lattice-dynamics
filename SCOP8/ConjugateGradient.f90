MODULE CG
    USE VA_math
    USE check
    IMPLICIT NONE

    INTEGER :: small_it_number=300

CONTAINS
!---------------------------------------------------------------------------
SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
    IMPLICIT NONE
    REAL(8) ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
    EXTERNAL func
    PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
    REAL(8) dum,fu,q,r,u,ulim
    fa=func(ax)
    fb=func(bx)
    if(fb.gt.fa )then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
    end if

    cx=bx+GOLD*(bx-ax)
    fc=func(cx)
1   if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
            fu=func(u)
            if(fu.lt.fc)then
                ax=bx
                fa=fb
                bx=u
                fb=fu
                return
            else if(fu.gt.fb)then
                cx=u
                fc=fu
                return
            endif
            u=cx+GOLD*(cx-bx)
            fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
            fu=func(u)
            if(fu.lt.fc)then
                bx=cx
                cx=u
                u=cx+GOLD*(cx-bx)
                fb=fc
                fc=fu
                fu=func(u)
            endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
            u=ulim
            fu=func(u)
        else
            u=cx+GOLD*(cx-bx)
            fu=func(u)
        end if
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
    end if
    return
END SUBROUTINE mnbrak
!-----------------------------------------------------------------------------------------
FUNCTION brent(ax,bx,cx,f,tol,xmin)
    INTEGER ITMAX_1
    REAL(8) brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
    EXTERNAL f
    PARAMETER (CGOLD=.3819660,ZEPS=1.0e-10)
    INTEGER iter
    REAL(8) a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    ITMAX_1=small_it_number
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.
    fx=f(x)
    fv=fx
    fw=fx
    do iter=1,ITMAX_1
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.*(q-r)
            if(q.gt.0.) p=-p
            q=abs(q)
            etemp=e
            e=d
            if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) goto 1
            d=p/q
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
        endif
1       if(x.ge.xm) then
            e=a-x
        else
            e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
            u=x+d
        else
            u=x+sign(tol1,d)
        end if
        fu=f(u)
        if(fu.le.fx) then
            if(u.ge.x) then
                a=x
            else
                b=x
            endif
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
        else
            if(u.lt.x) then
                a=u
            else
                b=u
            endif
            if(fu.le.fw .or. w.eq.x) then
                v=w
                fv=fw
                w=u
                fw=fu
            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
                v=u
                fv=fu
            endif
        endif
    enddo

    Write(*,*) 'WARNING: BRENT exceed maximum iterations'
3   xmin=x
    brent=fx
    return
END FUNCTION brent
!--------------------------------------------------------------------------------
SUBROUTINE linmin(p,xi,n,fret)!p->x,n->mx,xi->f
    IMPLICIT NONE
    INTEGER n,NMAX
    REAL(8) fret,p(n),xi(n),TOL
    PARAMETER (NMAX=2000,TOL=1.e-4)
    INTEGER j,ncom
    REAL(8) ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX)
    COMMON /f1com/ pcom,xicom,ncom
    ncom=n
    do j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
    enddo
    ax=0.
    xx=1.
    call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
    fret=brent(ax,xx,bx,f1dim,TOL,xmin)
    do j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
    enddo
    return
END SUBROUTINE linmin
!--------------------------------------------------------------------------------
FUNCTION f1dim(x)
    IMPLICIT NONE
    INTEGER NMAX
    REAL(8) f1dim,func,x
    PARAMETER(NMAX=2000)
    INTEGER j,ncom
    REAL(8) :: pcom(NMAX),xicom(NMAX),xt(NMAX)
    COMMON /f1com/ pcom,xicom,ncom
    REAL(8),ALLOCATABLE :: x_new(:),f_new(:),x_new2(:),f_new2(:)
    REAL(8) :: L1_norm

    do j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
    enddo

    ALLOCATE(x_new(ncom))
    DO j=1,ncom
        x_new(j)=xt(j)
    END DO

    ALLOCATE(f_new(SUM(variational_parameters_size)))

    !this is a comment: f1dim=func(xt) => calculate new free energy
    !NEED CHECK HERE
    !comment below block for test purpose
    !uncomment below block for actual running
    !---------update F:free energy--------------
    IF(unleashed) THEN
        CALL decompose_variants(x_new)
    ELSE
        CALL release_x(fixed_params,x_new,x_new2)
        CALL decompose_variants(x_new2)
    END IF

    CALL GetV_avg_And_GradientV_avg(kvector)!get new <V>
    CALL combine_gradients(f_new) !get gradients as a 1D array f(:)
    CALL UpdateTrialFC2 !get new eff fc2
    CALL GetEigen(kvector) !get new eivals
    CALL GetF0_and_V0 !get new F0 and <V0>
    F_trial = F0+V_avg-V0
    f1dim=REAL(F_trial)
    !-------------------------------------------
    IF(unleashed) THEN
        CALL printGradients(x_new)
        CALL calculate_threshold(x_new,f_new,L1_norm)
    ELSE
        CALL printGradients(x_new2)
        CALL select_xy(fixed_params,f_new,f_new2)
        CALL calculate_threshold(x_new,f_new2,L1_norm)
    END IF
    !---------------for test only---------------
!    f1dim=testF(xt,10)
    !-------------------------------------------

    !---------------check_small_convergence------------
    OPEN(55,FILE='convergence_small.dat',POSITION='append')
    WRITE(55,*) L1_norm,REAL(F_trial)
    CLOSE(55)
    !--------------------------------------------------

    DEALLOCATE(x_new,f_new)
    return
END FUNCTION f1dim
!----------------------------------------------------------------------------
SUBROUTINE update_fgh(mx,f,g,h)!update f,g,h
    IMPLICIT NONE
    INTEGER mx
    REAL(8) f(mx),g(mx),h(mx)
    INTEGER j
    REAL(8) dgg,gam,gg

    gg=0.
    dgg=0.
    DO j=1,mx
        gg=gg+g(j)**2
!        dgg=dgg+f(j)**2 !This statement for Fletcher-Reeves.
        dgg=dgg+(f(j)+g(j))*f(j) !This statement for Polak-Ribiere.
    END DO
    IF(gg.eq.0.) RETURN !unlikely

    gam=dgg/gg
    DO j=1,mx
        g(j)=-f(j) !g(:) is updated here, which is the -gradients of last loop
        h(j)=g(j)+gam*h(j)!h(:) is updated here
        f(j)=h(j) !gradients f(:) is updated here
    ENDDO

END SUBROUTINE update_fgh
!----------------------------------------------------------------------------
FUNCTION testF(x,mx)
    IMPLICIT NONE
    INTEGER mx
    INTEGER i
    REAL(8) x(mx), testF
    SELECTCASE(mx)
        CASE(1)
            testF = x(1)**4+x(1)**3+x(1)**2-9*x(1)
        CASE(2)
            testF = (x(1)+x(2)-3)**2+(x(1)-2*x(2))**2
        CASE(10)
            testF = (x(1)-1)**4
            DO i=2,10
                testF = testF + i*(x(i)-x(i-1)-1)**4
            END DO
    ENDSELECT

    return
END FUNCTION testF


FUNCTION testdF(x,mx)
    IMPLICIT NONE
    INTEGER mx
    INTEGER i
    REAL(8) x(mx), testdF(mx)

    SELECTCASE(mx)
        CASE(1)
            testdF(1) = 4*x(1)**3+3*x(1)**2+2*x(1)-9
        CASE(2)
            testdF(1) = 2*(x(1)+x(2)-3)+2*(x(1)-2*x(2))
            testdF(2) = 2*(x(1)+x(2)-3)-4*(x(1)-2*x(2))
        CASE(10)
            testdF(1) = -4*2*(x(2)-x(1)-1)**3 + 4*(x(1)-1)**3
            DO i=2,9
                testdF(i) = -4*i*(x(i+1)-x(i)-1)**3 + 4*i*(x(i)-x(i-1)-1)**3
            END DO
            testdF(10) = 4*10*(x(10)-x(9)-1)**3
    ENDSELECT
    return
END FUNCTION testdF


SUBROUTINE testCG
    IMPLICIT NONE
    INTEGER j,mx,i
    INTEGER counter
    REAL(8) :: fp,threshold,fret
    REAL(8),DIMENSION(:),ALLOCATABLE :: x,f,g,h

    mx=10 !test cases
    ALLOCATE(x(mx),f(mx),g(mx),h(mx))
!    CALL RANDOM_NUMBER(x)
!    x = 10*x
    x=[(i,i=1,10,1)]
    counter = 0
    fp = testF(x,mx)
    f = testdF(x,mx)
    g = -f
    h = g
    threshold = 100
    WRITE(*,*) 'f=',f
    WRITE(*,*) 'x=',x
    DO WHILE(threshold > 1e-8)
        f = h
        CALL linmin(x,f,mx,fret)
        fp = fret
        f = testdF(x,mx)
        WRITE(*,*) 'F=',fp
        WRITE(*,*) 'x=',x
        CALL update_fgh(mx,f,g,h)
        threshold = 0d0
        DO j=1,mx
            threshold = threshold + ABS(f(j)*x(j))
        END DO
        threshold = threshold/mx
        WRITE(*,*)'threshold=',threshold
!        PAUSE 'press to continue'
        counter = counter + 1
    END DO
    !----------------------------
    WRITE(*,*)'==== Finished! ===='
    WRITE(*,*) 'x=',x
    WRITE(*,*) '# of iterations: ',counter
END SUBROUTINE testCG
!----------------------------------------------------------------------------
!SUBROUTINE frprmn(p,n,ftol,iter,fret)!p->x,n->mx,xi->f
!    INTEGER iter,n,NMAX,ITMAX
!    REAL fret,ftol,p(n),EPS,func
!    EXTERNAL func
!    PARAMETER (NMAX=n,ITMAX=max_it_number,EPS=1.e-10)
!    INTEGER its,j
!    REAL dgg,fp,gam,gg,g(NMAX),h(NMAX),xi(NMAX)
!    fp=func(p)
!    call dfunc(p,xi)
!    do j=1,n
!        g(j)=-xi(j)
!        h(j)=g(j)
!        xi(j)=h(j)
!    enddo
!    do its=1,ITMAX
!        iter=its
!        call linmin(p,xi,n,fret)
!        if(2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS))return
!        fp=fret
!        call dfunc(p,xi) !replaced with <VA_math> scheme
!
!        gg=0.
!        dgg=0.
!        do  j=1,n
!            gg=gg+g(j)**2
!!           dgg=dgg+xi(j)**2 !This statement for Fletcher-Reeves.
!            dgg=dgg+(xi(j)+g(j))*xi(j) !This statement for Polak-Ribiere.
!        enddo
!        if(gg.eq.0.)return !unlikely
!        gam=dgg/gg
!        do j=1,n
!            g(j)=-xi(j)
!            h(j)=g(j)+gam*h(j)
!            xi(j)=h(j)
!        enddo
!    enddo
!    pause 'frprmn maximum iterations exceeded'
!    return
!END SUBROUTINE
!---------------------------------------------------------------------
!combine GradientV_### into 1D array
SUBROUTINE combine_gradients(f)
    IMPLICIT NONE
    INTEGER :: i,temp
    INTEGER :: tau1,atom2,direction1,direction2
    REAL(8),DIMENSION(:),INTENT(out) :: f

    !gradients of free energy w.r.t atomic deviation
    i=0
    DO while(i<variational_parameters_size(1))
        i=i+1
        IF(MOD(i,d).eq.0) THEN
            f(i)=GradientV_utau(d,INT(i/d))
        ELSE
            f(i)=GradientV_utau(MOD(i,d),INT(i/d+1))
        END IF
    END DO

    !gradients of free energy w.r.t strain
    DO while(i<variational_parameters_size(1)+variational_parameters_size(2))
        i=i+1
        temp=i-variational_parameters_size(1)
        IF(MOD(temp,d).eq.0) THEN
            f(i)=GradientV_eta(INT(temp/d),d)
        ELSE
            f(i)=GradientV_eta(INT(temp/d+1),MOD(temp,d))
        END IF
    END DO

    !gradients of free energy w.r.t <YY>, then subtract 1/2 effective fc2(which should give 0)
    !...as a substitute of gradients w.r.t effective fc2 itself
    DO while(i<SUM(variational_parameters_size))
        i=i+1
        temp=i-variational_parameters_size(1)-variational_parameters_size(2)
        tau1=eff_fc2_index(temp)%iatom_number
        atom2=eff_fc2_index(temp)%jatom_number
        direction1=eff_fc2_index(temp)%iatom_xyz
        direction2=eff_fc2_index(temp)%jatom_xyz

        f(i)=GradientV_cor(tau1,atom2)%phi(direction1,direction2)-&
        &0.5*trialfc2_value(tau1,atom2)%phi(direction1,direction2)
    END DO

END SUBROUTINE combine_gradients
!---------------------------------------------------------------------
SUBROUTINE calculate_threshold(x,f,L1_norm) !for every small iterations in f1dim
    IMPLICIT NONE
    INTEGER :: i,j
    REAL(8),DIMENSION(:),INTENT(in) :: x,f
    REAL(8),INTENT(out) :: L1_norm

    L1_norm = 0d0

    !threshold for all variables, just for reference
    IF(unleashed) THEN
        DO i=1,variational_parameters_size(1)+variational_parameters_size(2)
            L1_norm = L1_norm + ABS(f(i)*x(i))
        END DO
        DO i=variational_parameters_size(1)+variational_parameters_size(2)+1, SIZE(x)
            j = i - variational_parameters_size(1) - variational_parameters_size(2)
            L1_norm = L1_norm + ABS(f(i)*YY_record(j))
        END DO
        L1_norm = L1_norm/SIZE(GradientF_trial)/(temperature*100*h_plank*c_light*6.242d+18) !dimensionless
    ELSE
        L1_norm = normGradients(free_var,x,f)
        L1_norm = L1_norm/SIZE(free_var)/(temperature*100*h_plank*c_light*6.242d+18)
    END IF

END SUBROUTINE calculate_threshold
!---------------------------------------------------------------------
END MODULE CG
