 program drv
! standalone: constrained LASSO (lasso_constrained.f90) on dumped A, b, K, x0; args: epsratio nlam refit(0/1)
 use lasso_constrained, only : classo_cv
 implicit none
 integer, parameter :: dp = kind(1d0)
 integer nr, nc, nk, nsel, info, nlam, ir
 real(dp), allocatable :: A(:,:), b(:), K(:,:), x0(:), x(:)
 real(dp) eps, lam
 character(32) arg
 call get_command_argument(1, arg) ; read(arg,*) eps
 call get_command_argument(2, arg) ; read(arg,*) nlam
 call get_command_argument(3, arg) ; read(arg,*) ir
 open(1, file='lasso_in.bin', access='stream', form='unformatted', status='old')
 read(1) nr, nc, nk
 allocate(A(nr,nc), b(nr), K(nc,nk), x0(nc), x(nc))
 read(1) A, b, K, x0
 close(1)
 open(2, file='lasso_table.dat')
 call classo_cv(nr, nc, nk, A, b, K, x0, 1d-10, 5, nlam, eps, 1d-10, 300, ir == 1, x, lam, nsel, info, 2)
 close(2)
 open(3, file='lasso_out.bin', access='stream', form='unformatted', status='replace')
 write(3) x
 close(3)
 print '(a,i6,a,i2,a,es10.3)', ' selected', nsel, ' info', info, ' lam', lam
 end program drv
