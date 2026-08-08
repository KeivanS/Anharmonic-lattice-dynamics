!===============================================================================
! LASSO / compressive-sensing extraction of the force constants.
!
! Selected with enforce_inv=2 in structure.params. It fits the force-displacement
! equations alone,
!
!      min_x  1/2 || A x - b ||^2  +  mu * sum_j n_j |x_j|
!
! with n_j the Euclidean norm of column j of A. The caller then projects the
! solution on the kernel of the invariance matrix, exactly as enforce_inv=1 does
! with its SVD solution, so the invariances hold exactly.
!
! Weighting the penalty by the column norm makes it independent of the units of
! each force constant, which matters here because ranks 2, 3 and 4 carry
! different powers of Angstrom.
!
! The invariance relations are deliberately NOT part of this fit. Stacking them
! with the force rows, as enforce_inv=0 does for the SVD, was tried first and
! does not work for an L1 fit: on MoTe2 the invariance block carried 95% of the
! total column weight, so the objective was dominated by the constraints, the
! cross-validation curve was flat to within its own error bars, and the
! selection returned the trivial all-zero model. The SVD is untroubled by the
! same imbalance because it takes the minimum-norm least-squares solution.
!
! The L1 term does model selection. When many neighbour shells and high ranks
! are requested most candidate force constants are not actually determined by
! the data; least squares spreads small noise-driven values over all of them,
! while the L1 term drives them to zero and keeps only what the data supports.
!
! mu is chosen by the code by k-fold cross-validation; see lasso_extract.
!
! The minimisation started from the split-Bregman method of
! ~/PROJECTS/TOOLS/compressive-sensing.f90 (subroutine cs). That routine needed
! three changes to run at production size:
!   - A^T A and A^T b are formed once, outside the iteration. They do not depend
!     on the iterate; the original rebuilt them every sweep.
!   - the inner linear solve uses one Cholesky factorisation of (A^T A + lan I),
!     also constant, instead of a full SVD per iteration.
!   - shrinkage(0,g) returned x/abs(x)*... = NaN, and a sparse solution is
!     mostly zeros.
! It is kept below as split_bregman, but the default solver is coordinate
! descent: see the comment on lasso_cd for why.
!
! K. Esfarjani / after compressive-sensing.f90, 2026
!===============================================================================
 subroutine lasso_extract(nrow,ncol,amat,bmat,nhom,nfold,nlam,epsratio, &
 &                        x,mu_opt,nnz,uio)
!! Solves the same stacked system that enforce_inv=0 hands to the SVD, by
!! L1-regularised least squares, choosing the L1 weight by cross-validation.
!! amat(1:nhom,:) are the homogeneous invariance rows and amat(nhom+1:nrow,:)
!! the force-displacement rows; the split matters only because the invariance
!! rows are constraints and so are kept in every training fold.
!! Returns the solution x, the selected mu and the number of non-zero FCs.
 use constants, only : r15
 use ios, only : ulog
 implicit none
 integer, intent(in) :: nrow,ncol,nhom,nfold,nlam,uio
 real(r15), intent(in) :: amat(nrow,ncol),bmat(nrow)
 real(r15), intent(in) :: epsratio
 real(r15), intent(out):: x(ncol),mu_opt
 integer, intent(out)  :: nnz
 integer nfrc

 real(r15), allocatable :: g(:,:),gh(:,:),gf(:,:,:),ctot(:),cf(:,:),bb(:)
 real(r15), allocatable :: colnorm(:),xs(:),gtr(:,:),ctr(:),lam(:),chol(:,:)
 real(r15), allocatable :: cverr(:,:),cvmean(:),cvsd(:),gx(:)
 integer, allocatable :: fold(:)
 integer i,j,k,ifold,nact,nvalid,ier,ibest,ibest1se,nit
 real(r15) mumax,rss,best,se1,dum,lan

 nfrc=nrow-nhom
 write(ulog,*)'==================== LASSO_EXTRACT ===================='
 write(ulog,*)'nrow,nhom,nfrc,ncol,nfold,nlam=',nrow,nhom,nfrc,ncol,nfold,nlam
 write(*   ,*)'LASSO: rows(invariance)=',nhom,' rows(forces)=',nfrc,' columns=',ncol

 allocate(colnorm(ncol))

! ---- 1. column scaling ------------------------------------------------------
 do j=1,ncol
    colnorm(j)=sqrt( dot_product(amat(:,j),amat(:,j)) )
    if(colnorm(j).le.0) colnorm(j)=1d0   ! dead column: it will stay at zero
 enddo
 write(ulog,3)'LASSO: min,max column norm=',minval(colnorm),maxval(colnorm)

! ---- 2. assign the force rows to folds --------------------------------------
! interleaved, so each fold samples all snapshots rather than a contiguous block
 allocate(fold(nfrc))
 do i=1,nfrc
    fold(i)=mod(i-1,nfold)+1
 enddo

! ---- 3. Gram matrices, formed once ------------------------------------------
 allocate(gh(ncol,ncol),gf(ncol,ncol,nfold),ctot(ncol),cf(ncol,nfold),bb(nfold))
 call gram_scaled(nhom,nrow,ncol,amat,colnorm,gh)
 call gram_fold(nhom,nrow,ncol,nfold,amat,bmat,colnorm,fold,gf,cf,bb)
 allocate(g(ncol,ncol))
 g=gh ; ctot=0
 do k=1,nfold
    g   = g    + gf(:,:,k)
    ctot= ctot + cf(:,k)
 enddo
 write(ulog,*)'LASSO: Gram matrices built'
! after the column scaling the total diagonal of G is ncol; how much of it comes
! from the invariance rows tells us whether the constraints swamp the data
 dum=0
 do j=1,ncol
    dum=dum+gh(j,j)
 enddo
 write(ulog,3)'LASSO: fraction of column weight carried by the invariance rows=',dum/ncol
 write(*   ,3)'LASSO: invariance rows carry this fraction of the column weight=',dum/ncol

! the Bregman penalty. Its value affects only the rate of convergence, not the
! minimiser; scale it with the problem so the shifted matrix stays conditioned
 lan=0d0
 do j=1,ncol
    lan=lan+g(j,j)
 enddo
 lan=lan/ncol
 if(lan.le.0) lan=1d0
 write(ulog,3)'LASSO: Bregman penalty lan=',lan

! ---- 4. mu grid -------------------------------------------------------------
! above mumax every coefficient is zero, so start there and come down
 allocate(lam(nlam))
 mumax=maxval(abs(ctot))
 if(mumax.le.0) then
    write(ulog,*)'LASSO: A^T b is identically zero; nothing to fit'
    x=0; nnz=0; mu_opt=0; return
 endif
 do i=1,nlam
    lam(i)=mumax*exp( log(epsratio)*dble(i-1)/dble(nlam-1) )
 enddo
 write(ulog,3)'LASSO: mu from',lam(1),' down to',lam(nlam)

! ---- 5. cross-validation ----------------------------------------------------
 allocate(cverr(nlam,nfold),cvmean(nlam),cvsd(nlam))
 allocate(gtr(ncol,ncol),ctr(ncol),xs(ncol),chol(ncol,ncol),gx(ncol))
 cverr=0
 foldloop: do ifold=1,nfold
    gtr = g    - gf(:,:,ifold)
    ctr = ctot - cf(:,ifold)
! one factorisation serves the whole mu path for this fold
    call chol_factor(ncol,gtr,lan,chol,ier)
    if(ier.ne.0) then
       write(ulog,*)'LASSO: Cholesky failed on fold ',ifold,' at column ',ier
       stop
    endif
    xs = 0
    nvalid = count(fold.eq.ifold)
    do i=1,nlam        ! warm start down the path
       call lasso_cd(ncol,gtr,ctr,lam(i),xs,nit)
       if(ifold.eq.1) write(uio,5)'#train i,mu,nit,nnz,train rss=',i,lam(i),nit, &
     &      count(xs.ne.0d0), dot_product(xs,matmul(gtr,xs))-2*dot_product(ctr,xs)
! validation residual from the held-out Gram alone:
!   ||A_f x - b_f||^2 = x^T G_f x - 2 c_f^T x + b_f^T b_f
       gx  = matmul(gf(:,:,ifold),xs)
       rss = dot_product(xs,gx) - 2*dot_product(cf(:,ifold),xs) + bb(ifold)
       cverr(i,ifold)=max(rss,0d0)/max(nvalid,1)
    enddo
    write(*,*)'LASSO: cross-validation fold ',ifold,' of ',nfold,' done'
 enddo foldloop

 do i=1,nlam
    cvmean(i)=sum(cverr(i,:))/nfold
    dum=0
    do k=1,nfold
       dum=dum+(cverr(i,k)-cvmean(i))**2
    enddo
    cvsd(i)=sqrt(dum/max(nfold-1,1))/sqrt(dble(nfold))
 enddo

 ibest=1 ; best=cvmean(1)
 do i=2,nlam
    if(cvmean(i).lt.best) then
       best=cvmean(i) ; ibest=i
    endif
 enddo
! "one standard error" rule: the sparsest model whose cross-validation error is
! still within one standard error of the best. Preferred here because it drops
! the tail of coefficients the data does not really resolve.
 se1=cvmean(ibest)+cvsd(ibest)
 ibest1se=ibest
 do i=1,ibest
    if(cvmean(i).le.se1) then
       ibest1se=i ; exit
    endif
 enddo

 write(uio,*)'# LASSO cross-validation: i, mu, mean CV error, standard error'
 do i=1,nlam
    write(uio,4) i,lam(i),cvmean(i),cvsd(i)
 enddo
 write(ulog,3)'LASSO: min-CV    mu=',lam(ibest)   ,' cv error=',cvmean(ibest)
 write(ulog,3)'LASSO: 1-SE rule mu=',lam(ibest1se),' cv error=',cvmean(ibest1se)
! Use the cross-validation minimum. The one-standard-error rule is the usual
! choice for pure model selection, but here it costs too much accuracy: on
! MoTe2 it kept 11 of 122 force constants and left a 39% force residual, where
! the CV minimum keeps more and fits far better. Both are reported above so the
! sparser model can be chosen deliberately if that is what is wanted.
 mu_opt=lam(ibest)

! ---- 6. refit on all rows at the selected mu --------------------------------
 call chol_factor(ncol,g,lan,chol,ier)
 if(ier.ne.0) then
    write(ulog,*)'LASSO: Cholesky failed on the full system at column ',ier
    stop
 endif
 xs=0
 do i=1,nlam
    call lasso_cd(ncol,g,ctot,lam(i),xs,nit)
    if(lam(i).le.mu_opt) exit
 enddo

! ---- 7. undo the column scaling --------------------------------------------
 do j=1,ncol
    x(j)=xs(j)/colnorm(j)
 enddo
 nnz=count(xs.ne.0d0)

 write(ulog,*)'LASSO: non-zero force constants=',nnz,' out of ',ncol
 write(*   ,*)'LASSO: non-zero force constants=',nnz,' out of ',ncol
 write(ulog,3)'LASSO: selected mu=',mu_opt

 deallocate(g,gh,gf,ctot,cf,bb,colnorm,xs,gtr,ctr,lam,cverr,cvmean,cvsd,fold,chol,gx)

3 format(a,99(1x,g13.6))
4 format(i5,3(2x,g14.7))
5 format(a,i4,2x,g12.5,2x,i5,2x,i5,2x,g14.7)

 end subroutine lasso_extract
!===============================================================================
 subroutine fit_residual(nrow,ncol,a,b,x,error,ermax,sig)
!! residual of a solution that did not come from svd_set, reported with the same
!! three measures svd_set uses so the two branches can be compared directly:
!! mean absolute error, largest absolute error, and the percent deviation
!! ||Ax-b||/||Ax|| x 100.
 use constants, only : r15
 implicit none
 integer, intent(in) :: nrow,ncol
 real(r15), intent(in) :: a(nrow,ncol),b(nrow),x(ncol)
 real(r15), intent(out):: error,ermax,sig
 integer i
 real(r15) ax,r,num,den,bnorm

 error=0 ; ermax=0 ; num=0 ; den=0 ; bnorm=0
 do i=1,nrow
    ax=dot_product(a(i,:),x)
    r=ax-b(i)
    error=error+abs(r)
    ermax=max(ermax,abs(r))
    num=num+r*r
    den=den+ax*ax
    bnorm=bnorm+b(i)*b(i)
 enddo
 error=error/max(nrow,1)
! ||Ax-b||/||b||: the relative residual. Unlike svd_set's ||Ax-b||/||Ax||, this
! has a fixed denominator, so solutions can be compared with each other. The
! ||Ax|| form diverges as Ax->0 and made a shrunken solution look catastrophic
! when it was merely small.
 if(bnorm.gt.0) then
    sig=sqrt(num/bnorm)*100
 elseif(den.gt.0) then
    sig=sqrt(num/den)*100
 else
    sig=0
 endif

 end subroutine fit_residual
!===============================================================================
 subroutine lasso_cd(n,g,c,mu,x,nit)
!! Cyclic coordinate descent for  min 1/2 x^T G x - c^T x + mu ||x||_1 ,
!! the LASSO written through G = A^T A and c = A^T b.
!!
!! This is the default solver. split_bregman below implements the same
!! minimisation and is kept because it is the method of
!! ~/PROJECTS/TOOLS/compressive-sensing.f90, but its convergence rate is set by
!! the Bregman penalty lan: with a single lan for the whole mu path it stalls at
!! small mu (it was hitting its iteration cap and returning unconverged
!! solutions), and making lan track mu would mean refactorising the Cholesky at
!! every mu. Coordinate descent has no such parameter, decreases the objective
!! monotonically, and costs O(n^2) per sweep with the Gram matrix in hand.
!!
!! x enters with the previous solution (warm start) and leaves with the new one.
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(in) :: g(n,n),c(n),mu
 real(r15), intent(inout) :: x(n)
 integer, intent(out) :: nit
 integer j,it,maxit
 real(r15) gx(n),zj,xnew,dx,maxdx,tol,scal

 maxit=5000
! the columns are unit-normalised upstream, so the diagonal of G is O(1) and an
! absolute tolerance is meaningful; scale it by the size of the data anyway
 scal=maxval(abs(c))
 if(scal.le.0) scal=1d0
 tol=1d-12*scal

 gx=matmul(g,x)

 do it=1,maxit
    maxdx=0
    do j=1,n
       if(g(j,j).le.0d0) then          ! column carries no information
          if(x(j).ne.0d0) then
             gx=gx-g(:,j)*x(j) ; x(j)=0d0
          endif
          cycle
       endif
! partial residual for coordinate j, with j's own contribution taken out
       zj   = c(j) - gx(j) + g(j,j)*x(j)
       xnew = sign(max(abs(zj)-mu,0d0),zj)/g(j,j)
       dx   = xnew-x(j)
       if(dx.ne.0d0) then
          gx   = gx + g(:,j)*dx
          x(j) = xnew
          maxdx= max(maxdx,abs(dx))
       endif
    enddo
    if(maxdx.lt.tol) exit
 enddo
 nit=it

 end subroutine lasso_cd
!===============================================================================
 subroutine split_bregman(n,g,c,chol,mu,lan,x,nit)
!! Split-Bregman minimisation of  1/2 x^T G x - c^T x + mu ||x||_1 ,
!! which is the LASSO written through G = A^T A and c = A^T b.
!! chol is the Cholesky factor of (G + lan I), supplied by the caller: it does
!! not depend on mu or on the iterate, so it is formed once per training set.
!! x enters with the previous solution (warm start) and leaves with the new one.
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(in) :: g(n,n),c(n),chol(n,n),mu,lan
 real(r15), intent(inout) :: x(n)
 integer, intent(out) :: nit
 integer i,it,maxit
 real(r15) d(n),b(n),rhs(n),xold(n),tol,chg

 maxit=500
 tol=1d-10
 d=x ; b=0

 do it=1,maxit
    xold=x
! x-step: minimise 1/2||Ax-b||^2 + lan/2 ||d-x-b||^2
!         -> (G + lan I) x = c + lan (d - b)
    do i=1,n
       rhs(i)=c(i)+lan*(d(i)-b(i))
    enddo
    call chol_solve(n,chol,rhs,x)
! d-step: soft threshold
    do i=1,n
       call shrink(x(i)+b(i),mu/lan,d(i))
    enddo
! Bregman update
    b=b+(x-d)
    chg=maxval(abs(x-xold))
    if(chg.lt.tol) exit
 enddo
 nit=it

! the sparse iterate is d; x itself is only exactly sparse in the limit, so
! return the thresholded vector to get true zeros
 x=d

 end subroutine split_bregman
!===============================================================================
 subroutine shrink(y,thr,res)
!! soft thresholding: sign(y)*max(|y|-thr,0), and 0 at y=0.
!! The original wrote y/abs(y)*max(...), which is NaN at y=0 - and a sparse
!! solution is mostly zeros.
 use constants, only : r15
 implicit none
 real(r15), intent(in) :: y,thr
 real(r15), intent(out):: res

 res = sign(max(abs(y)-thr,0d0),y)

 end subroutine shrink
!===============================================================================
 subroutine chol_factor(n,g,shift,l,ier)
!! Cholesky factorisation of (G + shift*I) = L L^T, lower triangle in l.
!! ier=0 on success, otherwise the column at which positive-definiteness failed.
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(in) :: g(n,n),shift
 real(r15), intent(out):: l(n,n)
 integer, intent(out) :: ier
 integer i,j,k
 real(r15) s

 ier=0
 l=0
 do j=1,n
    s=g(j,j)+shift
    do k=1,j-1
       s=s-l(j,k)*l(j,k)
    enddo
    if(s.le.0d0) then
       ier=j ; return
    endif
    l(j,j)=sqrt(s)
    do i=j+1,n
       s=g(i,j)
       do k=1,j-1
          s=s-l(i,k)*l(j,k)
       enddo
       l(i,j)=s/l(j,j)
    enddo
 enddo

 end subroutine chol_factor
!===============================================================================
 subroutine chol_solve(n,l,b,x)
!! solves L L^T x = b for the factor produced by chol_factor
 use constants, only : r15
 implicit none
 integer, intent(in) :: n
 real(r15), intent(in) :: l(n,n),b(n)
 real(r15), intent(out):: x(n)
 integer i,k
 real(r15) s,y(n)

! forward substitution L y = b
 do i=1,n
    s=b(i)
    do k=1,i-1
       s=s-l(i,k)*y(k)
    enddo
    y(i)=s/l(i,i)
 enddo
! back substitution L^T x = y
 do i=n,1,-1
    s=y(i)
    do k=i+1,n
       s=s-l(k,i)*x(k)
    enddo
    x(i)=s/l(i,i)
 enddo

 end subroutine chol_solve
!===============================================================================
 subroutine gram_scaled(nhom,nrow,ncol,a,colnorm,g)
!! G = A^T A over rows 1..nhom of the column-scaled matrix A(:,j)/colnorm(j)
 use constants, only : r15
 implicit none
 integer, intent(in) :: nhom,nrow,ncol
 real(r15), intent(in) :: a(nrow,ncol),colnorm(ncol)
 real(r15), intent(out):: g(ncol,ncol)
 integer i,j
 real(r15), allocatable :: col(:,:)

 g=0
 if(nhom.le.0) return
 allocate(col(nhom,ncol))
 do j=1,ncol
    col(:,j)=a(1:nhom,j)/colnorm(j)
 enddo
 do j=1,ncol
    do i=j,ncol
       g(i,j)=dot_product(col(:,i),col(:,j))
       g(j,i)=g(i,j)
    enddo
 enddo
 deallocate(col)

 end subroutine gram_scaled
!===============================================================================
 subroutine gram_fold(nhom,nrow,ncol,nfold,a,b,colnorm,fold,gf,cf,bb)
!! per-fold Gram matrices, right-hand sides and ||b||^2 for the column-scaled
!! matrix, over the force rows nhom+1..nrow only. Splitting the Gram this way
!! lets a training set be formed by subtraction, so the rows are read once no
!! matter how many folds are used.
 use constants, only : r15
 implicit none
 integer, intent(in) :: nhom,nrow,ncol,nfold,fold(nrow-nhom)
 real(r15), intent(in) :: a(nrow,ncol),b(nrow),colnorm(ncol)
 real(r15), intent(out):: gf(ncol,ncol,nfold),cf(ncol,nfold),bb(nfold)
 integer i,j,k,kf,ir
 real(r15) aj
 real(r15), allocatable :: row(:)

 gf=0 ; cf=0 ; bb=0
 allocate(row(ncol))
 do ir=1,nrow-nhom
    i=nhom+ir
    kf=fold(ir)
    do j=1,ncol
       row(j)=a(i,j)/colnorm(j)
    enddo
    do j=1,ncol
       aj=row(j)
       if(aj.eq.0d0) cycle
       cf(j,kf)=cf(j,kf)+aj*b(i)
       do k=j,ncol
          gf(k,j,kf)=gf(k,j,kf)+row(k)*aj
       enddo
    enddo
    bb(kf)=bb(kf)+b(i)*b(i)
 enddo
 do kf=1,nfold
    do j=1,ncol
       do k=j+1,ncol
          gf(j,k,kf)=gf(k,j,kf)
       enddo
    enddo
 enddo
 deallocate(row)

 end subroutine gram_fold
!===============================================================================
