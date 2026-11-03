# Exact test of born_flag 2 (Coulomb subtraction), 2026-11-03

What FOCEX does with born_flag 2 (ewald_born.f90, coulomb_lr.f90):
- Subtracted from the forces: F_sub = F_m(x) - F_m(x0) - (Phi_B - Phi_m) u, where
  - F_m is the Ewald force of point charges with strength c_ij = tr(Z_i Z_j^T)/3, screened by eps;
  - Phi_m = -e^2 c_ij H_ij is its harmonic part (H = grad grad phi, image-summed);
  - Phi_B = -e^2 Z_i H_ij Z_j^T is the Born-dipole term (the same D_NA that is added back to D).
  So the harmonic part subtracted = the one added back, and the anharmonic part is the trace model's.
- Invariance rows of the short-range fit, x_fit = x_total - x_sub, become A x_fit = -A x_sub:
  - translational: unchanged (F_sub is exactly translation-invariant);
  - Huang: b = -2V(d_ab s_gd - d_gd s_ab) - (S_B(ab,gd) - S_B(gd,ab)), with S_B = sum Phi_B r_g r_d;
  - rank-2 rotational: b = -sum_{a,mu} M(i;a1,a,mu) eps(a,mu,nu), M = sum_j Phi_B(i a1, j a)(x_j - x_i)_mu.
  S_B and M are conditionally convergent: small-q expansion of the E = 0 Ewald C(q) (lr_brackets,
  lr_moments; central differences, Richardson).
  - FC1 is not subtracted: the fitted FC1 = total FC1, i.e. it includes -F_C(x0).
- The mech.dat checks (Huang, rank-2 rotational) add the same Phi_B brackets/moments.
- Rank-3 and higher rotational rows: homogeneous. Exact for scalar Z (the trace model is then rotationally
  invariant and Phi_B = Phi_m); for anisotropic Z not yet corrected (see "Open").

## Scripts

- `mkcoulomb.py BASE OUT [q1,q2,..] [eps]`: truth = the exact short-range model of BASE (model3, stressed
  reference) + point charges per species in scalar eps, Ewald-summed. Writes FORCEDISP1 (model + Coulomb),
  stress.params (model + Coulomb, (1/V) dE/deta by differences), dielectric.params (born_flag 2), coulomb.json
  (F0_coulomb, Coulomb stress). Self-checks: Ewald lambda-independence, F = -dE/dx.
- `compare.py BASE RUN...`: max |x_run - x_base| of fcN_irr.dat per rank, absolute and relative.
- For scalar charges FOCEX's subtraction equals the true Coulomb force minus its reference value, so the fit
  must give BASE's short-range FC2.., and FC1 shifted by -F_C(x0).

## Results (FOCEX built from the working tree of 2026-11-03; enforce_inv = 1)

Runs were made with the right-hand sides switched off one at a time, to show that each is needed
(test-only switches, since removed from the code).

ZnO (run_ZnO, ranks 1-4, charges +2/-2, eps 5; F_C(x0) 1.8e-4 eV/A):

| run | fit | FC2 rel | FC3 rel | FC4 rel |
|---|---|---|---|---|
| both RHS | 0.000 % | 5.7e-9 | 3.3e-7 | 8.2e-8 (printed digits) |
| Huang RHS off | 0.026 % | 6.0e-5 | 1.4e-3 | 1.3e-2 |
| rotational RHS off | 0.000 % | 1.6e-7 | 3.7e-6 | 6.0e-5 |
| both off | 0.026 % | 6.0e-5 | 1.4e-3 | 1.3e-2 |

rel = max |x - x_base| / max |x_base| of the rank (compare.py). FC1 shift (both RHS) 1.811e-4 = -F_C(x0).

HfO2 (cubic-truncated model, `mkcubic.py HfO2 cubic`, ranks 1-3, charges +2/-1, eps 5;
F_C(x0) up to 0.376 eV/A, Coulomb stress yz -0.53 GPa):

| run | fit | FC2 rel | FC3 rel |
|---|---|---|---|
| no charges (base) | 0.000 % | | |
| both RHS | 0.000 % | 1.3e-8 | 4.5e-8 (printed digits) |
| Huang RHS off | 0.095 % | 7.2e-4 | 3.4e-3 |
| rotational RHS off | 0.510 % | 4.9e-3 | 2.8e-2 |
| both off | 0.523 % | 5.2e-3 | 3.1e-2 |

- FC1 shift (both RHS) = -F_C(x0) of the three inequivalent atoms to all 7 printed digits.
- mech.dat checks (both RHS): Huang 1.4e-8 GPa, rotational 1.8e-13 eV/A.
- The full HfO2 model with ranks 1-4 (7977 unknowns) takes ~9 h, mostly the Numerical Recipes SVD in
  get_kernel; not needed: the rank-4 rows have no Coulomb term for scalar charges, and ZnO tests them.
- User's ZnO DFT data (born_flag 2, enforce_inv 0; before = FOCEX before the hybrid, same inputs): fit
  1.720 -> 1.269 %; long-wave max |B - B^T| (E = 0) 1.36 -> 0.017 GPa; mech.dat Huang check (fitted + Born
  brackets) 4.7e-5 GPa. The long-wave C is still unphysical (rho v^2 < 0 along several directions).

## Open

- Rank-3 rotational rows for anisotropic Z: right-hand side -A_full x_sub (not yet implemented, no exact
  test with anisotropic Z yet).
