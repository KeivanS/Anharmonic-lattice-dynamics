module params
 use constants, only : r15
 real(r15) tolerance,margin,scalelengths,alfaborn,tempk,rcutoff,coef
 integer nconfigs,classical,fdfiles,cal_cross,threemtrx,lamin,lamax,ncpu,n_dig_acc,itemp,convcoord
 integer , allocatable :: nshells(:,:)  !nshells(8,natom_prim_cell) up to which shell to include for each rank of FC (read from input file)
 integer include_fc(8), icutoff !,nsmax  ! whether to include FCs of that rank ,max# of shells looping
 real(r15) tmin,tmax,qcros(3),svdc,lmax ! lmax is the cutoff length of FC2 limited by the supercell WS
 logical verbose,readfc
! LASSO (enforce_inv=2). The L1 weight lasso_mu is not an input: it is chosen by
! cross-validation inside lasso_extract and reported here. The three settings
! below control that search and rarely need changing.
 integer  :: lasso_nfold   = 5      ! folds of the cross-validation
 integer  :: lasso_nlambda = 40     ! points on the mu path
 real(r15):: lasso_epsratio= 1d-4   ! smallest mu, as a fraction of the largest
 real(r15):: lasso_mu               ! the selected value (output)

end module params
!============================================================
