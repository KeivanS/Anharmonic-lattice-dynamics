!===============================================================================
! Constrained fitting of the force constants: least squares (or LASSO) performed
! INSIDE the kernel of the invariance matrix, so that the translational,
! rotational and Huang relations hold exactly.
!
! Why this exists. enforce_inv=1 used to solve the force equations over all
! dim_ac force constants and then orthogonally project that solution on the
! kernel with project_on. An orthogonal projection is not a constrained fit:
!
!    projection      x_p = K K^T x_unconstrained
!    constrained fit x_c = K argmin_y || A K y - b ||
!
! The two agree only if A^T A acts as a multiple of the identity. In general the
! projection discards the information about how the residual varies inside the
! directions it removes, and can be arbitrarily bad. Measured on MoTe2 (122 force
! constants, 87 independent invariance relations, so a 35-dimensional kernel),
! with the relative residual ||Ax-b||/||b|| on the force equations:
!
!    unconstrained fit                    ~13 %
!    orthogonal projection on the kernel  ~3400 %
!    refit inside the kernel basis        ~35 %
!
! K, the orthonormal kernel basis returned by get_kernel, IS the elimination of
! the dependent force constants: with nker free parameters y, x = K y produces
! all dim_ac force constants at once, the eliminated ones included, so there is
! no separate back-substitution step. It is the same elimination one would get
! by row-reducing the constraint matrix, but in an orthonormal and therefore
! far better conditioned basis.
!
! Nothing in svd_ridge.f90, extratools.f90 or lasso.f90 is modified by this
! file; get_kernel, svd_set, project_on and lasso_extract are used as they are.
! Note that svd_set overwrites its matrix argument (svdcmp returns U in it), so
! the reduced matrix built here is a local copy and the caller's aforce is left
! untouched.
!
! K. Esfarjani, 2026
!===============================================================================
 subroutine particular_solution(m,n,c,d,svdcut,xp,resid)
!! Minimum-norm solution xp of the invariance rows C xp = d (C: m x n), by SVD with the
!! cut svdcut used for the kernel; xp is orthogonal to the kernel, so the solutions of
!! C x = d are x = xp + K y. Needed when some rows have a right-hand side: the Huang rows
!! of a stressed reference (stress.params). resid = max|C xp - d| (a check that the rows
!! are consistent).
 use constants, only : r15
 implicit none
 integer, intent(in) :: m,n
 real(r15), intent(in) :: c(m,n),d(m),svdcut
 real(r15), intent(out):: xp(n),resid
 real(r15), allocatable :: u(:,:),v(:,:),w(:),dd(:)
 integer j,mm
 mm = max(m,n)                ! svdcmp on a square-or-tall copy: zero rows added if m < n
 allocate(u(mm,n),v(n,n),w(n),dd(mm))
 u = 0 ; u(1:m,:) = c ; dd = 0 ; dd(1:m) = d
 call svdcmp(u,mm,n,w,v)
 xp = 0
 do j=1,n
    if (abs(w(j)) .gt. svdcut*maxval(abs(w))) xp = xp + v(:,j)*dot_product(u(:,j),dd)/w(j)
 enddo
 resid = maxval(abs(matmul(c,xp)-d))
 deallocate(u,v,w,dd)
 end subroutine particular_solution
!===============================================================================
 subroutine fit_kernel_basis(nrow,ncol,nker,a,b,kernel,svdcut,x,sigx,error,ermax,sig,fnsvd)
!! Constrained least squares. Minimises ||A x - b|| over the x that satisfy the
!! invariance relations exactly, by fitting the nker free parameters y of
!! x = K y and returning x.
!!   nrow  : number of force equations
!!   ncol  : number of force constants (dimension of x)
!!   nker  : dimension of the kernel = number of free parameters
!!   a,b   : force-displacement matrix and right-hand side (not modified)
!!   kernel: kernel(ncol,>=nker), orthonormal basis from get_kernel
!!   sigx  : standard error of each x(j): sqrt([K Cov(y) K^T]_jj), Cov(y) from svd_set (residual
!!           variance s^2 = ||A K y - b||^2/(nrow-nker))
 use constants, only : r15
 use ios, only : ulog
 implicit none
 integer, intent(in) :: nrow,ncol,nker
 real(r15), intent(in) :: a(nrow,ncol),b(nrow),kernel(ncol,nker),svdcut
 real(r15), intent(out):: x(ncol),sigx(ncol),error,ermax,sig
 character(len=*), intent(in) :: fnsvd
 real(r15), allocatable :: ared(:,:),y(:),sigred(:),vs(:,:),kvs(:,:)
 integer j

 write(ulog,*)'FIT_KERNEL_BASIS: constrained least squares in the kernel basis'
 write(ulog,*)'FIT_KERNEL_BASIS: rows,columns,free parameters=',nrow,ncol,nker

 if(nker.le.0) then
    write(ulog,*)'FIT_KERNEL_BASIS: kernel is empty; the invariances leave no'
    write(ulog,*)' free parameter. Check itrans/irot/ihuang and the FC ranges.'
    x=0 ; sigx=0 ; error=0 ; ermax=0 ; sig=0
    return
 endif

 allocate(ared(nrow,nker),y(nker),sigred(nker),vs(nker,nker))
! A K : the force equations expressed in the free parameters
 ared = matmul(a,kernel(:,1:nker))
! svd_set overwrites ared, which is why it is a local copy
 call svd_set(nrow,nker,ared,b,y,sigred,vs,svdcut,error,ermax,sig,fnsvd)
! back to the full set of force constants, eliminated ones included
 x = matmul(kernel(:,1:nker),y)
! Cov(x) = K vs vs^T K^T: standard errors from the rows of K vs
 deallocate(ared)
 allocate(kvs(ncol,nker))
 kvs = matmul(kernel(:,1:nker),vs)
 do j=1,ncol
    sigx(j)=sqrt(dot_product(kvs(j,:),kvs(j,:)))
 enddo

 write(ulog,*)'FIT_KERNEL_BASIS: done'
 deallocate(y,sigred,vs,kvs)

 end subroutine fit_kernel_basis
!===============================================================================
 subroutine lasso_support_errors(nrow,ncol,nker,a,b,kernel,x,svdcut,sigx)
!! Standard errors of a LASSO solution x (enforce_inv=2,3; 2026-11-06, user: sd was 0 there).
!! Support S as in lasso_constrained: |w_j x_j| > 1d-8 max, w_j = column rms of a (floored).
!! Least squares on S with the invariances exact: x + K N z, N = null space of the rows of K
!! off S (null_space_qr, as refit_support), by fit_kernel_basis on the basis K N; s^2 from the
!! residual of that refit. enforce_inv=3: x is that refit, so this is its standard error.
!! enforce_inv=2: the error the refit would have; neither includes the selection itself or
!! the LASSO shrinkage. x is not changed.
 use constants, only : r15
 use ios, only : ulog
 use lasso_constrained, only : null_space_qr
 implicit none
 integer, intent(in) :: nrow,ncol,nker
 real(r15), intent(in) :: a(nrow,ncol),b(nrow),kernel(ncol,nker),x(ncol),svdcut
 real(r15), intent(out):: sigx(ncol)
 real(r15), allocatable :: w(:),bt(:,:),nb(:,:),kn(:,:),dx(:)
 logical, allocatable :: sel(:)
 real(r15) error,ermax,sig
 integer i,j,m,p

 sigx=0
 if(nker.le.0) return
 allocate(w(ncol),sel(ncol))
 do j=1,ncol
    w(j)=sqrt(dot_product(a(:,j),a(:,j))/nrow)
 enddo
 if(maxval(w).le.0) then
    w=1
 else
    w=max(w,1d-8*maxval(w))
 endif
 sel = abs(w*x) .gt. 1d-8*maxval(abs(w*x))
 m = count(.not.sel)
 allocate(nb(nker,nker))
 if(m.eq.0) then
    p=nker ; nb=0
    do j=1,nker ; nb(j,j)=1 ; enddo
 else
    allocate(bt(nker,m)) ; i=0
    do j=1,ncol
       if(sel(j)) cycle
       i=i+1 ; bt(:,i)=kernel(j,:)
    enddo
    call null_space_qr(nker,m,bt,1d-9,p,nb)
    deallocate(bt)
 endif
 write(ulog,'(a,i7,a,i7,a,i7)')'LASSO_SUPPORT_ERRORS: FCs kept',count(sel),' of',ncol, &
 &     ', free parameters of the refit on them',p
 if(p.gt.0) then
    allocate(kn(ncol,p),dx(ncol))
    kn = matmul(kernel,nb(:,1:p))
    call fit_kernel_basis(nrow,ncol,p,a,b-matmul(a,x),kn,svdcut,dx,sigx,error,ermax,sig,'svd-lasso.dat')
    write(ulog,'(a,g12.4)')'LASSO_SUPPORT_ERRORS: max |refit - LASSO| (not applied)=',maxval(abs(dx))
    deallocate(kn,dx)
 endif
 deallocate(w,sel,nb)

 end subroutine lasso_support_errors
!===============================================================================
 subroutine lasso_kernel_basis(nrow,ncol,nker,a,b,kernel,x0,nfold,nlam,epsratio,refit, &
 &                             x,mu_opt,nnz,uio)
!! LASSO on the force constants x themselves with the invariances exact
!! (x in x0 + range(kernel); x0 = particular solution, 0 unless some invariance rows have a
!! right-hand side, see particular_solution), then a constrained least-squares refit on the selected
!! terms if refit (enforce_inv=3); the L1 weight by nfold cross-validation of the
!! model returned (lasso_constrained.f90). nnz = number of force constants selected.
!! The refit removes the LASSO shrinkage; use it when the forces are accurate. With
!! noisy forces or too few snapshots it fits the selected terms to the noise:
!! enforce_inv=2 (no refit) is then the safer choice.
!! FOCEX_LASSO=old selects the previous method, LASSO on the kernel coefficients y
!! (lasso_kernel_basis_old): on the exact cubic model that one returned 98/80/376
!! non-zero FC2/FC3/FC4 for 11/20/43 true ones, and FC4 wrong.
 use constants, only : r15
 use ios, only : ulog
 use lasso_constrained, only : classo_cv
 use svd_stuff, only : svdcut
 implicit none
 integer, intent(in) :: nrow,ncol,nker,nfold,nlam,uio
 real(r15), intent(in) :: a(nrow,ncol),b(nrow),kernel(ncol,nker),x0(ncol),epsratio
 logical, intent(in)   :: refit
 real(r15), intent(out):: x(ncol),mu_opt
 integer, intent(out)  :: nnz
 character(len=16) envv
 integer info

 call get_environment_variable('FOCEX_LASSO',envv)
 if(trim(envv).eq.'old') then
    write(ulog,*)'LASSO_KERNEL_BASIS: FOCEX_LASSO=old, LASSO on the kernel coefficients'
    write(*   ,*)'LASSO_KERNEL_BASIS: FOCEX_LASSO=old, LASSO on the kernel coefficients'
    if(refit) write(*,*)'WARNING: FOCEX_LASSO=old has no refit; enforce_inv=3 acts as 2'
! L1 on the kernel coefficients y of x = x0 + K y
    call lasso_kernel_basis_old(nrow,ncol,nker,a,b-matmul(a,x0),kernel,nfold,nlam,epsratio, &
 &       x,mu_opt,nnz,uio)
    x = x + x0
    return
 endif
 write(ulog,*)'LASSO_KERNEL_BASIS: constrained LASSO on the force constants + refit'
 if(nker.le.0) then
    write(ulog,*)'LASSO_KERNEL_BASIS: kernel is empty; nothing to fit'
    x=x0 ; nnz=0 ; mu_opt=0
    return
 endif
 call classo_cv(nrow,ncol,nker,a,b,kernel,x0,svdcut,nfold,nlam,epsratio,1d-10,300,refit, &
 &    x,mu_opt,nnz,info,uio)
 write(ulog,'(a,i6,a,i6,a,es12.4,a,i2)')'LASSO: force constants selected=',nnz,' of ',ncol, &
 &    '  lam=',mu_opt,'  info=',info
 write(*   ,'(a,i6,a,i6,a,i2)')'LASSO: force constants selected=',nnz,' of ',ncol,'  info=',info
 if(info.eq.2) then
    write(ulog,*)'WARNING: LASSO: the interior point did not converge at some lam (see above)'
    write(*   ,*)'WARNING: LASSO: the interior point did not converge at some lam (see above)'
 elseif(info.eq.3) then
    write(ulog,*)'WARNING: LASSO: the CV minimum is at the smallest lam of the grid (epsratio too large?)'
    write(*   ,*)'WARNING: LASSO: the CV minimum is at the smallest lam of the grid (epsratio too large?)'
 elseif(info.eq.1) then
    write(ulog,*)'WARNING: LASSO: W K is rank deficient'
    write(*   ,*)'WARNING: LASSO: W K is rank deficient'
 endif

 end subroutine lasso_kernel_basis
!===============================================================================
 subroutine lasso_kernel_basis_old(nrow,ncol,nker,a,b,kernel,nfold,nlam,epsratio, &
 &                             x,mu_opt,nnz,uio)
!! Same as fit_kernel_basis but the reduced problem is solved by LASSO instead
!! of by SVD, so that the free parameters the data does not support are set to
!! zero. The invariances still hold exactly, since x = K y for whatever y.
!!
!! Caveat worth keeping in mind: the sparsity is in the free parameters y, not
!! in the force constants x. x = K y is a dense combination, so a sparse y does
!! not by itself mean few non-zero force constants. Making x itself sparse while
!! keeping the invariances exact is a generalized-lasso problem and is not what
!! this routine solves.
 use constants, only : r15
 use ios, only : ulog
 implicit none
 integer, intent(in) :: nrow,ncol,nker,nfold,nlam,uio
 real(r15), intent(in) :: a(nrow,ncol),b(nrow),kernel(ncol,nker),epsratio
 real(r15), intent(out):: x(ncol),mu_opt
 integer, intent(out)  :: nnz
 real(r15), allocatable :: ared(:,:),y(:)

 write(ulog,*)'LASSO_KERNEL_BASIS: LASSO in the kernel basis'

 if(nker.le.0) then
    write(ulog,*)'LASSO_KERNEL_BASIS: kernel is empty; nothing to fit'
    x=0 ; nnz=0 ; mu_opt=0
    return
 endif

 allocate(ared(nrow,nker),y(nker))
 ared = matmul(a,kernel(:,1:nker))
! nhom=0: there are no constraint rows left, they are built into the basis
 call lasso_extract(nrow,nker,ared,b,0,nfold,nlam,epsratio,y,mu_opt,nnz,uio)
 x = matmul(kernel(:,1:nker),y)

 write(ulog,*)'LASSO_KERNEL_BASIS: non-zero free parameters=',nnz,' of ',nker
 deallocate(ared,y)

 end subroutine lasso_kernel_basis_old
!===============================================================================
