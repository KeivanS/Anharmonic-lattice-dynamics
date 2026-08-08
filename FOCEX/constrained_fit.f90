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
 subroutine fit_kernel_basis(nrow,ncol,nker,a,b,kernel,svdcut,x,error,ermax,sig,fnsvd)
!! Constrained least squares. Minimises ||A x - b|| over the x that satisfy the
!! invariance relations exactly, by fitting the nker free parameters y of
!! x = K y and returning x.
!!   nrow  : number of force equations
!!   ncol  : number of force constants (dimension of x)
!!   nker  : dimension of the kernel = number of free parameters
!!   a,b   : force-displacement matrix and right-hand side (not modified)
!!   kernel: kernel(ncol,>=nker), orthonormal basis from get_kernel
 use constants, only : r15
 use ios, only : ulog
 implicit none
 integer, intent(in) :: nrow,ncol,nker
 real(r15), intent(in) :: a(nrow,ncol),b(nrow),kernel(ncol,nker),svdcut
 real(r15), intent(out):: x(ncol),error,ermax,sig
 character(len=*), intent(in) :: fnsvd
 real(r15), allocatable :: ared(:,:),y(:),sigred(:)

 write(ulog,*)'FIT_KERNEL_BASIS: constrained least squares in the kernel basis'
 write(ulog,*)'FIT_KERNEL_BASIS: rows,columns,free parameters=',nrow,ncol,nker

 if(nker.le.0) then
    write(ulog,*)'FIT_KERNEL_BASIS: kernel is empty; the invariances leave no'
    write(ulog,*)' free parameter. Check itrans/irot/ihuang and the FC ranges.'
    x=0 ; error=0 ; ermax=0 ; sig=0
    return
 endif

 allocate(ared(nrow,nker),y(nker),sigred(nker))
! A K : the force equations expressed in the free parameters
 ared = matmul(a,kernel(:,1:nker))
! svd_set overwrites ared, which is why it is a local copy
 call svd_set(nrow,nker,ared,b,y,sigred,svdcut,error,ermax,sig,fnsvd)
! back to the full set of force constants, eliminated ones included
 x = matmul(kernel(:,1:nker),y)

 write(ulog,*)'FIT_KERNEL_BASIS: done'
 deallocate(ared,y,sigred)

 end subroutine fit_kernel_basis
!===============================================================================
 subroutine lasso_kernel_basis(nrow,ncol,nker,a,b,kernel,nfold,nlam,epsratio, &
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

 end subroutine lasso_kernel_basis
!===============================================================================
