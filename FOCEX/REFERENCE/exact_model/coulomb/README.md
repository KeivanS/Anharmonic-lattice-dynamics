# Exact test of born_flag 2 (Coulomb subtraction), 2026-11-03

What FOCEX does with born_flag 2 (ewald_born.f90, coulomb_lr.f90):
- Subtracted from the forces: F_sub = F_m(x) - F_m(x0) - (Phi_C - Phi_m) u, where
  - F_m is the Ewald force of point charges with strength c_ij = tr(Z_i Z_j^T)/3, screened by eps;
  - Phi_m = -e^2 c_ij H_ij is its harmonic part (H = grad grad phi, image-summed);
  - Phi_C = -e^2 Z_i H_ij Z_j^T is the Coulomb dipole term (the same D_NA that is added back to D).
  So the harmonic part subtracted = the one added back, and the anharmonic part is the trace model's.
- Invariance rows of the short-range fit, x_fit = x_total - x_sub, become A x_fit = -A x_sub:
  - translational: unchanged (F_sub is exactly translation-invariant);
  - Huang: b = -2V(d_ab s_gd - d_gd s_ab) - (S_C(ab,gd) - S_C(gd,ab)), with S_C = sum Phi_C r_g r_d;
  - rank-2 rotational: b = -sum_{a,mu} M(i;a1,a,mu) eps(a,mu,nu), M = sum_j Phi_C(i a1, j a)(x_j - x_i)_mu.
  S_C and M are conditionally convergent: small-q expansion of the E = 0 Ewald C(q) (lr_brackets,
  lr_moments; central differences, Richardson).
  - FC1 is not subtracted: the fitted FC1 = total FC1, i.e. it includes -F_C(x0).
- The mech.dat checks (Huang, rank-2 rotational) add the same Phi_C brackets/moments.
- Rank-3 rotational rows (i1 a1, i2 a2; nu): b = -A_full x_sub (lr_rot3_pair, lr_rot3_onsite), x_sub = Phi_C (rank 2)
  + the pair model's FC3 (rank 3), r = X_i1 - X_i2:
  - pair, i2 /= i1 (bare, one pair of the crystal): b = e2 [c phi_(a1 a2 a) r_mu eps(a,mu,nu) + (Z H Z^T)_(g a2) eps(g,a1,nu)
    + (Z H Z^T)_(a1 g) eps(g,a2,nu)];
  - on-site, i2 = i1: b = -sum_l b_pair(i1,l), as E = 0 Ewald lattice sums; the phi''' part through the identity
    (J r).grad H - J H - H J^T = D_eps H [eps J - J eps] (rotating r at fixed eps = rotating eps the other way).
  - Zero for scalar Z and eps. Rank 4 and higher: homogeneous (exact for scalar eps; not corrected for tensor eps).

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

## Rank-3 rotational right-hand side: checks (2026-11-03)

- Scalar charges (HfO2 cubic, `r3_on`): max |b| = 8e-16; FCs = the run without it, byte for byte.
- User's ZnO (tensor Z and eps; scratch z_b2_r3): max |b| 0.278 eV/A^2; fit 1.269 -> 1.232 %; long-wave
  |B - B^T| 0.017 GPa, Huang check 4.4e-5 GPa.
- Independent Python (scratch rot3chk/): 30 pair rows from mpmath derivatives of phi_eps: |b - b_FOCEX| <= 2.3e-7
  (positions printed to 1e-6 A); on-site rows by rotating the lattice positions at fixed eps (own Ewald, E = 0,
  Richardson in the angle) instead of FOCEX's eps-derivative: <= 5.2e-7 (3e-6 relative); FOCEX's ASR-enforced Z.
  A double-precision finite-difference third derivative was not enough: the row is the difference of two terms of
  ~7 eV/A^2 and ends at ~0.04.
- No exact end-to-end test with anisotropic Z: the hybrid's subtracted FC3 is not the rotational partner of Phi_C,
  so the true short-range FC3 has a 1/r^4 tail and no finite-range fit is exact.

## include_fc = 2 (read FCs), tests (2026-11-03; scratch coul_cubic/rf_*)

HfO2 cubic exact model; fcN.dat and lat_fc.dat from the full fits c_base / c_on (rel = vs the full fit):

| run | read | FC1 rel | FC2 rel | FC3 rel |
|---|---|---|---|---|
| rf_base, rf_on | none | byte-identical to the full fits (force rows rewritten) | | |
| rf_b_r3 | FC3 | 2.1e-7 | 2.6e-7 | read |
| rf_b_r2 | FC2 | 2.6e-5 | read | 1.1e-4 |
| rf_o_r2 (charges) | FC2 | 3.9e-6 | read | 1.3e-4 |
| rf_o_r12 (charges) | FC1, FC2 | read | read | 3.1e-6 |

- Every term of the read fcN.dat matched (atoms via lat_fc.dat's (tau, n)); least-squares residual of the
  irreducible values <= 2e-16.
- These runs used 7-digit files. The FC3 error with FC2 read is the rounding of the read FC2 (p8_*, p15_*: FC2 read from a
  full fit written with 7, 8 and, in a scratch build only, 15 digits):

  | digits of fc2.dat | FC1 rel | FC3 rel |
  |---|---|---|
  | 7 | 2.6e-5 | 1.1e-4 |
  | 8 (FOCEX since 2026-11-03) | 1.4e-5 | 5.3e-5 |
  | 15 | 5.8e-10 | 2.4e-9 |

  FC3 is that sensitive to FC2 (re-rounding fc2.dat at the 3e-7 level moved FC3 by 4.7e-4); the mechanism was not
  analysed. FC3 read (8 digits): FC1, FC2 off by 2e-8, 4e-8.
- Atom numbering (p8_r2perm, p8_r3perm; scratch rot3chk/permute.py): the earlier run's atoms 13..6500 renumbered by a
  random permutation in lat_fc.dat and the read fcN.dat: every output (fc1-3, fc*_irr, lat_fc.dat, mech.dat)
  byte-identical to the unpermuted read; the rewritten read rank is in this run's numbering (= the original file).
- The rewritten fcN.dat of a read rank is byte-identical to the one read (same numbering).

## Open

- Rank-4 rotational rows for tensor eps (the trace model is then not rotationally invariant at rank 4).
