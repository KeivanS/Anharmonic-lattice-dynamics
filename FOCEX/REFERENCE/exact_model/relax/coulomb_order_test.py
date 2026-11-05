"""coulomb_order_test.py COULOMBRUN MINJSON: is FOCEX's relaxed state (born_flag 2) the minimum of
   exact short-range model + Coulomb energy to 2nd order about the reference?           (2026-11-04)
FOCEX's relaxation uses the fitted FCs (= the exact short-range model, mkcoulomb.py subtracts the full
Ewald force) plus the Coulomb terms FC1/stress (reference), Phi_C, q_c, a_c: the Coulomb energy to 2nd
order in (u, eta); Psi_C, chi_C = 0. Evaluated at FOCEX's state q_F (strain from POSCAR_relaxed, u0 from
the log, as compare_relaxed.py), with the functions of model_min_c.py:
  g_exact = grad_model(q_F) + grad_coul(q_F)                       (true model)
  g_2     = grad_model(q_F) + grad_coul(0) + H_C q_F                (Coulomb to 2nd order)
H_C q_F = d/dt grad_coul(t q_F) at t = 0, 4-point central differences (t = +-0.01, +-0.02).
If FOCEX minimizes the 2nd-order energy, g_2 ~ 0 (difference accuracy) while g_exact is not; then
the difference to the exact minimum (MINJSON, model_min_c.py) is the Coulomb energy beyond 2nd order."""
import sys, os, re, json, numpy as np
EM = os.path.expanduser('~/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/REFERENCE/exact_model')
crun, mj = sys.argv[1], sys.argv[2]
src = open(f'{EM}/relax/model_min_c.py').read().split('q = np.zeros(n)\nHi = None')[0]
argv = sys.argv; sys.argv = ['model_min_c.py', crun, '/dev/null']
exec(src)                       # Xi, latx, ewald, grad_model, grad_coul, eta_of, N, n, base, ...
sys.argv = argv
# FOCEX's relaxed state
d = crun
log = open([os.path.join(d, f) for f in os.listdir(d) if f.startswith('log') and f.endswith('.dat')][0]).read()
conv, prim, atoms = read_structure(f'{d}/structure.params')
Pn = np.array([[float(v) for v in l.split()[:3]] for l in open(f'{d}/POSCAR_relaxed').read().splitlines()[2:5]])
Fm = np.linalg.solve(prim, Pn).T
e6 = np.array([Fm[0, 0] - 1, Fm[1, 1] - 1, Fm[2, 2] - 1, Fm[1, 2], Fm[0, 2], Fm[0, 1]])
u0 = np.array([[float(v) for v in l.split('=')[1].split()[1:]] for l in re.findall(r'RELAX_REFERENCE: tau, u0\(tau\) \(A\) =.*', log)])
tau = []
for x in Xi:
    t = [k for k, (_, p) in enumerate(atoms) if np.allclose(np.round(np.linalg.solve(prim.T, x - p)), np.linalg.solve(prim.T, x - p), atol=1e-6)]
    tau.append(t[0])
tau = np.array(tau)
def state(e, u0p):
    return np.concatenate([e, u0p[tau].ravel()])
qF = state(e6, u0)
mn = json.load(open(mj)); qX = state(np.array(mn['strain']), np.array(mn['u0']))
def report(name, g):
    print(f'  {name:46s} max|dE/deta| {np.abs(g[:6]).max():.3e} eV   max|force| {np.abs(g[6:]).max():.3e} eV/A', flush=True)
g0c = grad_coul(np.zeros(n))
def hc_dir(q, t=0.01):
    gp1, gm1 = grad_coul(t * q), grad_coul(-t * q)
    gp2, gm2 = grad_coul(2 * t * q), grad_coul(-2 * t * q)
    return (8 * (gp1 - gm1) - (gp2 - gm2)) / (12 * t)
for name, q in (('FOCEX relaxed state', qF), ('exact minimum (model_min_c)', qX)):
    gm, gc = grad_model(q), grad_coul(q)
    hq = hc_dir(q)
    print(name, ': strain', np.round(q[:6], 7).tolist())
    report('exact model, full Ewald', gm + gc)
    report('exact model, Ewald to 2nd order', gm + g0c + hq)
    report('Coulomb beyond 2nd order (full - 2nd order)', gc - g0c - hq)
# accuracy of the directional difference: same with t halved
hq2 = hc_dir(qF, 0.005)
print('difference step check at q_F: max|H_C q (t=0.01) - H_C q (t=0.005)| =', np.abs(hc_dir(qF) - hq2).max())
