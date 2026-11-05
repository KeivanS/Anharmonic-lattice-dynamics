"""coulomb_rotation_test.py COULOMBRUN: is the Coulomb energy to 2nd order about the reference rotationally
invariant at FOCEX's relaxed state?                                                      (2026-11-04)
State: deformation F (3x3, any) and displacements U; X = Xi F^T + U, lattice latx F^T. A rigid rotation of
the deformed crystal, X -> X (1+w)^T with w antisymmetric, is F -> (1+w) F, U -> U (1+w)^T. The exact Ewald
energy does not change (torque 0). Its 2nd-order Taylor E2(q) = E(0) + g.q + q.H.q/2 in q = (F - 1, U) has
torque  T_w = d/ds E2(q + s v_w) = v.grad E(0) + d/dt [v.grad E(t q)]_{t=0},  v_w = (w F, U w^T),
computed by central differences of the Ewald energy (inner s: +-1e-4, 4-point; outer t: +-0.01, +-0.02).
If |T_w| is comparable to the 0.15 eV strain mismatch, a 2nd-order
Coulomb energy is not a unique target: FOCEX (polar decomposition, Lagrangian elastic constants) and the
F = 1 + eta parametrization of coulomb_order_test.py differ at 3rd order, and Psi_C is needed."""
import sys, os, re, json, numpy as np
EM = os.path.expanduser('~/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/REFERENCE/exact_model')
crun = sys.argv[1]
src = open(f'{EM}/relax/model_min_c.py').read().split('q = np.zeros(n)\nHi = None')[0]
argv = sys.argv; sys.argv = ['model_min_c.py', crun, '/dev/null']
exec(src)
sys.argv = argv
d = crun
log = open([os.path.join(d, f) for f in os.listdir(d) if f.startswith('log') and f.endswith('.dat')][0]).read()
conv, prim, atoms = read_structure(f'{d}/structure.params')
Pn = np.array([[float(v) for v in l.split()[:3]] for l in open(f'{d}/POSCAR_relaxed').read().splitlines()[2:5]])
FF = np.linalg.solve(prim, Pn).T
u0 = np.array([[float(v) for v in l.split('=')[1].split()[1:]] for l in re.findall(r'RELAX_REFERENCE: tau, u0\(tau\) \(A\) =.*', log)])
tau = []
for x in Xi:
    t = [k for k, (_, p) in enumerate(atoms) if np.allclose(np.round(np.linalg.solve(prim.T, x - p)), np.linalg.solve(prim.T, x - p), atol=1e-6)]
    tau.append(t[0])
UF = u0[np.array(tau)]
def EC(Fm, U):
    return ewald(Xi @ Fm.T + U, latx @ Fm.T, lam)[0]
def dd(Fm, U, vF, vU, h=1e-4):          # directional derivative of E_C at (F, U) along (vF, vU)
    f = lambda s: EC(Fm + s * vF, U + s * vU)
    return (8 * (f(h) - f(-h)) - (f(2 * h) - f(-2 * h))) / (12 * h)
I = np.eye(3); DF = FF - I
def torque2(vF, vU, t=0.01):            # v.grad E2 at q_F = v.grad E(0) + d/dt [v.grad E(t q_F)]
    g = lambda s: dd(I + s * DF, s * UF, vF, vU)
    return dd(I, 0 * UF, vF, vU) + (8 * (g(t) - g(-t)) - (g(2 * t) - g(-2 * t))) / (12 * t)
def exact(vF, vU):
    return dd(FF, UF, vF, vU)
gens = {'w_yz': (1, 2), 'w_xz': (0, 2), 'w_xy': (0, 1)}
print('FOCEX relaxed F - 1 (should be symmetric):\n', np.round(DF, 7))
for name, (a, b) in gens.items():
    w = np.zeros((3, 3)); w[a, b] = 1; w[b, a] = -1
    vF, vU = w @ FF, UF @ w.T
    print(f'{name}: torque of exact Ewald {exact(vF, vU):.3e} eV,  of Ewald to 2nd order {torque2(vF, vU):.3e} eV', flush=True)
