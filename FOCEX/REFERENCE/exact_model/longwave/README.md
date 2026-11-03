# Elastic constants and stress from the sound speeds (2026-11-02)

Tests of FOCEX's `elastic_longwave` (FOCEX/dynamical.f90) on the exact model of `../` (no Born charges).

- `longwave_check.py MAT` (Python, exact Hessian of the model): the stressed reference of run_MAT with the
  internal coordinates relaxed at fixed lattice (force-free, stressed); long-wave acoustic tensor along the 13
  directions with the optical modes folded out; fit of the 36 brackets; deviatoric stress from their
  asymmetry against the exact stress; C from B - d_ik s_jl against the exact relaxed c (strain Hessian with the
  internal coordinates relaxed). Logs `longwave_check_{ZnO,HfO2,CaF2}.log`: s_dev to 1e-13 GPa, C to 5e-10
  GPa; without the stress term C is 15 (ZnO), 5 (HfO2), 15 (CaF2) GPa off. (CaF2: s_dev = 0 by symmetry, so
  its relative residual is 0/0.)
- `mkstressed.py ZnO DIR`: FOCEX inputs at that state (same 100 displacements, exact forces, stress.params =
  exact stress); the x, y coordinates were then reset to the strings of run_ZnO/structure.params (the
  relaxation is along z; round-off 1e-13). FOCEX 9.4 + elastic_longwave: `focex_stressed.log`, s_dev
  (0.03609, 0.03609, -0.07219) GPa = exact, C = exact (printed digits).
- `focex_relax.log`: run_ZnO with relax_reference (default.params line 5 = 1, as ../relax): s_dev 0, C =
  the exact relaxed strain Hessian of ../relax/paper_A_check_ZnO_min.log (V = 45.2995 A^3) to 1e-3 GPa.

```bash
python longwave_check.py ZnO > longwave_check_ZnO.log     # ~3 s
python mkstressed.py ZnO /tmp/zno_stressed                 # then reset x, y as above
cd /tmp/zno_stressed && ~/ALADYN/bin/focex-9.5.x > focex.out   # mech.dat: ELASTIC_LONGWAVE lines
```

## coul_brackets.py (2026-11-03)

Brackets of the Born-charge (dipole) force constants alone, Ewald with E = 0, for the user's ZnO
(Z*, eps_inf and structure hard-coded from its dielectric.params / lat_fc.dat). `python3 coul_brackets.py`
(~1.5 min; `neutral` removes the small net charge). Prints: Ewald-parameter check (4e-11), the 36 brackets,
their deviatoric stress (-0.108, -0.108, +0.216 GPa) and the antisymmetric part of Q (0.134 eV/A; it has the
residual-force form, pi_z = -0.133 Zn / +0.133 O eV/A, see CLAUDE.md).
