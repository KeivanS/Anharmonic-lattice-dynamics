module params
 use constants, only : r15
 real(r15) tolerance,margin,scalelengths,alfaborn,tempk,rcutoff,coef
! same atom position (Angstrom, Cartesian distance): accepting a symmetry operation and finding the images
! of atoms under it use this one value (2026-11-02, user: a ZnO relaxed 8e-5 A off wurtzite passed the
! symmetry finder at 1e-3 but failed the image lookup at 1e-4)
 real(r15), parameter :: tol_sym = 1d-3
 integer nconfigs,classical,fdfiles,cal_cross,threemtrx,lamin,lamax,ncpu,n_dig_acc,itemp,convcoord
 integer , allocatable :: nshells(:,:)  !nshells(8,natom_prim_cell) up to which shell to include for each rank of FC (read from input file)
 integer include_fc(8), icutoff !,nsmax  ! whether to include FCs of that rank ,max# of shells looping
 real(r15) tmin,tmax,qcros(3),svdc,lmax ! lmax is the cutoff length of FC2 limited by the supercell WS
 logical verbose
! LASSO (enforce_inv=2). The L1 weight lasso_mu is not an input: it is chosen by
! cross-validation inside lasso_extract and reported here. The three settings
! below control that search and rarely need changing.
 integer  :: lasso_nfold   = 5      ! folds of the cross-validation
 integer  :: lasso_nlambda = 40     ! points on the mu path
 real(r15):: lasso_epsratio= 1d-4   ! smallest mu, as a fraction of the largest
 real(r15):: lasso_mu               ! the selected value (output)
! lasso_rule (default.params line 6, absent = 0; 2026-11-06, user): how lam is chosen. 0 = smallest CV
! error; 1 = noise test, the largest lam whose refit is consistent with the force noise (classo_cv_res)
 integer  :: lasso_rule = 0
! lasso_l2 (default.params line 7, absent or 0 = off): L2 term with prior width tau_r = lasso_l2 x rms of the
! rank-r FCs of the least-squares fit (classo_cv_res)
 real(r15):: lasso_l2 = 0
! svd_noise_check: svd_set compares the measured relative noise of its fit with line 2's guess and raises the
! cut (svd_ridge.f90). Only for fits of the forces themselves: focex.f90 sets it around them (off for fits of
! residuals, e.g. lasso_support_errors, whose relative residual is ~1 by construction)
 logical  :: svd_noise_check = .false.
! relax_ref (default.params line 5, absent = 0): 1 = relax the reference with the fitted FCs and use the
! re-expanded FCs downstream (FOCEX relax_reference); ref_relaxed is set once that is done
 integer  :: relax_ref = 0
 logical  :: ref_relaxed = .false.

end module params
!============================================================
