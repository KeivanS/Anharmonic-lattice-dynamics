# Ballistic conductance test on exact chain models (2026-11-03)

Test of `conductance.dat` (FOCEX 9.5; dynamical.f90). Requested by K. Esfarjani; scripts and runs by
Claude. Every number below was printed by `compare.py`.

## What FOCEX computes

M(ω)/A = 2π Σ_λ ∫ d³k/(2π)³ δ(ω − ω_kλ) |v_kλ·n|/2: the number of ballistic channels per Å² of cross
section normal to n (n from `latdyn.params` line 2, normalized).

- |v·n|/2 summed over the full zone equals the sum of v·n over v·n > 0, because v(−k) = −v(k). The plain
  sum of v·n over the zone is zero, which is why that weight was not used.
- With ω in 1/cm the code evaluates M/A = [(1/N) Σ_k δ(ν − ν_k) |v·n|/2] / (c V), with v and c in m/s and
  V in Å²·cm. The k sum uses the tetrahedron method (`tet_sum`) on the mesh `kpc`.
- Frequencies and velocities on the mesh come from the IBZ by the Laue operation R of each mesh point:
  v(Rk) = R v(k). This also holds for R = −S (time reversal). `ibz_reduce` returns the operation in its
  optional argument `oparm`.
- In 1D, M·A_cell is the number of branches that cross ω with v > 0, an integer. With transverse
  dispersion it becomes a continuous function.

## Models (`setup_chain.py`, model of ../model3.py with the harmonic term only)

- Lattice: tetragonal a × a × 2.5 Å, one Ar atom, chains along z.
- Interaction: E = Σ L s over the pairs with r < rL (s = squared distance), with L = K0 (1 − s0/rL²)².
- These springs are isotropic, so the three acoustic branches are degenerate everywhere.
- The exact reference is `ref.py`: a fine-mesh histogram of the same formula, computed from the
  analytic dispersion, independent of FOCEX.

| run | a (Å) | rL (Å) | neighbours (Å) | expected |
|---|---|---|---|---|
| B | 8.4 | 8.6 | 2.5, 5.0, 7.5 along the chain; 8.4 between chains (L = 0.25% of the first) | M·a² = 3 up to ν(π) = 49.34 1/cm; 6 up to the band maximum 51.03 (at k = 0.654 π/c); 0 above |
| D | 3.5 | 3.6 | 2.5 along; 3.5 between (1.1%) | ≈ 3, nearly integer |
| C | 3.5 | 4.6 | 2.5, 3.5, 4.30 | continuous |
| C123 | as C | | normal (1,2,3)/√14, fixed by no operation | continuous |

## Inputs

- **Displacements:** Gaussian, 0.04 Å per component, not from SC_SNAPS. SC_SNAPS builds its sampling
  model from springs with r < 2 r_nn (`sc_snaps.f90:557`) and does not sample modes below fmin (lines
  98–100). In B every such spring lies along the chain, so the snapshots had no x or y displacements.
  FOCEX then fitted transverse FC2 = 0, the elastic tensor was singular, and `xmatinv` stopped.
- **Mesh:** B uses 8×8×96 (the 3→6 window is 3% of ω); C and C123 use 12×12×48; D uses 8×8×48.

## Results (`run_chains.sh WORK`, then `compare.py` in WORK)

- **B, plateau 8–48.5 1/cm:** M·a² mean 2.9990 (min 2.82, max 3.18).
- **B, window 49.8–50.6 1/cm:** mean 6.0059 (5.83–6.27). Above 51.03 FOCEX gives 0.
- **D, plateau 5–26 1/cm:** mean 2.9998 (2.87–3.13).
- **C:** rms |FOCEX − exact| / mean(exact) = 0.0080 over 3–46 1/cm.
- **C123:** 0.0086. The values differ from C (e.g. 2.405 against 3.000 at 38 1/cm), so the projection
  onto n and the rotation of v are both exercised.
- **Remaining error:** the scatter about the plateaus is that of the k mesh. The biggest deviations are at
  the band edges and below about 5 1/cm, where the mesh resolves the dispersion poorly.

## Regression

- **Chains B and C123, tagged 9.5 binary against the new code:** only conductance.dat (it held zeros
  before), times and the line numbers of a runtime warning differ.
- **Ge example (12³ mesh):** the same. The tetrahedron DOS plus conductance block takes 1.7 s, against 0.9 s
  for the DOS alone.
- **THERMACOND Ge 4³:** outputs byte-identical. Only times, the run path and the random-seed line differ.
