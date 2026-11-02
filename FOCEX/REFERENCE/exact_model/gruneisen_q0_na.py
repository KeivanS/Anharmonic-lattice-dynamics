"""gruneisen_q0_na.py: the long-wave Gruneisen limit with a non-analytic (Born-charge) term in D, on the
exact ZnO model with chosen charges: Z(Zn) = diag(2.0, 2.0, 2.1) = -Z(O) (sum rule exact), eps = diag(5, 5, 5.5).
  D(q) = D_SR(q) + NA(qh),  NA_{ta,t'b} = (4 pi e^2/V) (qh.Z_t)_a (qh.Z_t')_b / (qh.eps.qh) / sqrt(m m')
(rigid-ion form, constant along a ray; W has no NA part, as in the codes). The limit uses D0 = D_SR0 + NA(qh),
D1, D2 of D_SR. Checks direct -> limit as q -> 0, and how much the NA term changes the acoustic gamma."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gruneisen_q0 as G
from model import read_structure

E2 = 14.399645                                   # e^2/(4 pi eps0), eV A
s = G.setup('ZnO', relax=True)
conv, prim, atoms = read_structure(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ZnO', 'structure.params'))
V = abs(np.linalg.det(prim))
Zzn = np.diag([2.0, 2.0, 2.1])
Z = [Zzn if sp == 1 else -Zzn for sp, _ in atoms]
eps = np.diag([5.0, 5.0, 5.5])
m = s['m']; nat = s['nat']


def na(qh):
    M = np.zeros((3 * nat, 3 * nat))
    den = qh @ eps @ qh
    for t in range(nat):
        for tp in range(nat):
            M[3 * t:3 * t + 3, 3 * tp:3 * tp + 3] = 4 * np.pi * E2 / V * np.outer(qh @ Z[t], qh @ Z[tp]) / den / np.sqrt(m[t] * m[tp])
    return M


orig_blocks, orig_taylor = G.blocks, G.taylor
use_na = True


def blocks(s_, q, kind):
    M = orig_blocks(s_, q, kind)
    if kind == 'D' and use_na and np.linalg.norm(q) > 0:
        M = M + na(q / np.linalg.norm(q))
    return M


def taylor(s_, qh, kind, n):
    M = orig_taylor(s_, qh, kind, n)
    if kind == 'D' and n == 0 and use_na:
        M = M + na(qh)
    return M


G.blocks, G.taylor = blocks, taylor
print(f'ZnO model + NA term (Z_zz(Zn) 2.1, eps_zz 5.5), with relaxation; V = {V:.3f} A^3')
for name, v in (('a1', conv[0]), ('c', conv[2]), ('uinc', G.UINC)):
    qh = v / np.linalg.norm(v)
    use_na = True
    a, b, gl, ck = G.gamma_limit(s, qh)
    use_na = False
    a0, b0, gl0, _ = G.gamma_limit(s, qh)
    use_na = True
    print(f' direction {name}: limit with NA {np.round(gl, 6)} (omega^2/q^2 {np.round(a, 5)}); without NA {np.round(gl0, 6)}')
    for q in (1e-2, 1e-3, 3e-4):
        w2, num, gd = G.gamma_direct(s, q * qh)
        print(f'    q={q:6.0e}: direct {np.round(gd[:3], 6)}  max |direct - limit| {np.abs(gd[:3] - gl).max():.2e}')
