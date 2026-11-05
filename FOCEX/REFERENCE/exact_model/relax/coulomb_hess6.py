"""coulomb_hess6.py COULOMBRUN [FOCEX_RUNDIR]: exact d^2E/de_k de_l (eV) at the reference, atoms fixed (u = 0),
F = 1 + eta(e), Voigt e (xx yy zz yz xz xy, off-diagonal eta_ab = eta_ba = e), for the exact short-range
model and for the Ewald point charges (model_min_c.py functions; 4-point central differences of their
strain gradients, h = 1e-3). With FOCEX_RUNDIR, compared with FOCEX's log lines "RELAX_REFERENCE: H6"
(total) and "H6C" (Coulomb part from the a_c brackets), printed by hess6_reference (2026-11-04)."""
import sys, os, re, numpy as np
EM = os.path.expanduser('~/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/REFERENCE/exact_model')
crun = sys.argv[1]; frun = sys.argv[2] if len(sys.argv) > 2 else None
src = open(f'{EM}/relax/model_min_c.py').read().split('q = np.zeros(n)\nHi = None')[0]
argv = sys.argv; sys.argv = ['model_min_c.py', crun, '/dev/null']
exec(src)
sys.argv = argv
def h6(gfun, h=1e-3):
    H = np.zeros((6, 6))
    for k in range(6):
        e = np.zeros(n); e[k] = h
        H[:, k] = ((8 * (gfun(e) - gfun(-e)) - (gfun(2 * e) - gfun(-2 * e))) / (12 * h))[:6]
    return 0.5 * (H + H.T)
np.set_printoptions(precision=6, suppress=True, linewidth=150)
Hm = h6(grad_model); print('exact model part (eV):\n', Hm, flush=True)
Hc = h6(grad_coul);  print('exact Coulomb (Ewald) part (eV):\n', Hc, flush=True)
print('exact total (eV):\n', Hm + Hc)
if frun:
    log = open([os.path.join(frun, f) for f in os.listdir(frun) if f.startswith('log') and f.endswith('.dat')][0]).read()
    Ft = np.array([[float(v) for v in l.split()] for l in re.findall(r'RELAX_REFERENCE: H6 (.*)', log)])
    Fc = np.array([[float(v) for v in l.split()] for l in re.findall(r'RELAX_REFERENCE: H6C (.*)', log)])
    print('FOCEX total (eV):\n', Ft); print('FOCEX Coulomb part (eV):\n', Fc)
    print('FOCEX - exact: total      max', np.abs(Ft - Hm - Hc).max(), '\n', Ft - Hm - Hc)
    print('               Coulomb    max', np.abs(Fc - Hc).max(), '\n', Fc - Hc)
    print('               short-range (total - Coulomb) max', np.abs(Ft - Fc - Hm).max())
