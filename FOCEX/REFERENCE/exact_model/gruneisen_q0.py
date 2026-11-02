"""gruneisen_q0.py MAT [relax]: test of the long-wave (q -> 0) limit of the mode Gruneisen parameters of
the acoustic branches on the exact model (model3.py), against direct evaluation at small q.

Exact inputs: FC2 = analytic Hessian; W = sum_k Psi(i,j,k) y_k = [H(X + l y) - H(X - l y)] / 2l, exact
because the model is a polynomial of degree 4 (FC5 = 0). y_k = x_k (no internal relaxation, THERMACOND)
or x_k + xi_k (FOCEX: xi = -Gam Q, the relaxation of the sublattices under isotropic strain).

Conventions (FOCEX): D(q)_{ta,t'b} = sum_{j in t'} Phi(ta,jb) e^{iq.(x_j - x_t)} / sqrt(m m'), same for W;
gamma = -<e|W|e> / (6 omega^2).
Long wave along qh: D(q) = D0 + q D1 + q^2 D2 + ..., Dn = (1/n!) sum (i qh.r)^n Phi/sqrt(mm'), same for W.
T = rigid translations (sqrt(m_t/M) delta), P = 1 - T T^+, G = (P D0 P)^+ (pseudo-inverse).
  A = T^+ D2 T - T^+ D1 G D1 T                      omega_s^2 = q^2 a_s + O(q^3)
  B = T^+ W2 T - T^+ W1 G D1 T - T^+ D1 G W1 T + T^+ D1 G W0 G D1 T      <e_s|W|e_s> = q^2 b_s + O(q^3)
  gamma_s(qh) = -b_s / (6 a_s), B diagonalized inside each degenerate subspace of A.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import read_structure, cluster
from model3 import Model

CFG = {'CaF2': (5.0, 4.2), 'ZnO': (4.72, 3.5), 'HfO2': (3.87, 3.22)}
MASS = {'Ca': 40.078, 'F': 18.9984, 'Zn': 65.409, 'O': 15.9994, 'Hf': 178.49}
UINC = np.array([np.sqrt(2) - 1, np.sqrt(3) - 1, np.sqrt(5) - 2])


def setup(mat, relax):
    here = os.path.dirname(os.path.abspath(__file__))
    conv, prim, atoms = read_structure(f'{here}/{mat}/structure.params')
    SL = open(f'{here}/{mat}/structure.params').read().split('\n')
    names = SL[9].split('#')[0].split()
    m = np.array([MASS[names[sp - 1]] for sp, _ in atoms])
    rL, Ra = CFG[mat]
    R = max(np.linalg.norm(p) for _, p in atoms) + rL + 0.5
    X, S, TAU, NN = cluster(prim, atoms, R)
    nat = len(atoms)
    cen = [int(np.where((TAU == t) & (np.abs(NN).sum(1) == 0))[0][0]) for t in range(nat)]
    mod = Model(rL, Ra).build(X, lat=None, centers=cen)
    N = len(X)
    H0 = mod.hessian(X).reshape(N, 3, N, 3)
    # check FC2 against an independent route: the force derivative
    rows = [(cen[t], a) for t in range(nat) for a in range(3)]
    Phi = np.array([H0[i, a] for i, a in rows])                      # (3nat, N, 3)
    # Phi(q=0) without mass, Gam = inverse on dofs 4.. (atom 1 fixed) as FOCEX get_phi_zeta_Xi
    phi0 = np.zeros((3 * nat, 3 * nat)); qiu = np.zeros((3 * nat, 3, 3))
    for r, (i, a) in enumerate(rows):
        for j in range(N):
            tp = TAU[j]
            phi0[r, 3 * tp:3 * tp + 3] += Phi[r, j]
            qiu[r] += np.outer(Phi[r, j], X[j] - X[i])                 # Q(ta; g,l) = sum Phi(ta,jg) r_l
    qiu = (qiu + qiu.transpose(0, 2, 1)) / 2
    gam = np.zeros_like(phi0)
    if nat > 1:
        gam[3:, 3:] = np.linalg.inv(phi0[3:, 3:])
    xi_tr = -gam @ np.einsum('sii->s', qiu)                          # du(t,g) per unit isotropic strain
    c = xi_tr.reshape(nat, 3) if relax else np.zeros((nat, 3))
    y = X + c[TAU]
    lam = 1e-3
    Wf = ((mod.hessian(X + lam * y) - mod.hessian(X - lam * y)) / (2 * lam)).reshape(N, 3, N, 3)
    Wr = np.array([Wf[i, a] for i, a in rows])
    return dict(X=X, TAU=TAU, cen=cen, nat=nat, m=m, Phi=Phi, W=Wr, phi0=phi0, xi=c, mod=mod)


def blocks(s, q, kind):
    """D(q) (or W(q)) in FOCEX's phase convention; kind 'D' or 'W'"""
    F = s['Phi'] if kind == 'D' else s['W']
    nat, X, TAU, cen, m = s['nat'], s['X'], s['TAU'], s['cen'], s['m']
    M = np.zeros((3 * nat, 3 * nat), complex)
    for t in range(nat):
        r = X - X[cen[t]]
        ph = np.exp(1j * r @ q)
        for tp in range(nat):
            sel = TAU == tp
            M[3 * t:3 * t + 3, 3 * tp:3 * tp + 3] = np.einsum('ajb,j->ab', F[3 * t:3 * t + 3][:, sel], ph[sel]) / np.sqrt(m[t] * m[tp])
    return M


def taylor(s, qh, kind, n):
    F = s['Phi'] if kind == 'D' else s['W']
    nat, X, TAU, cen, m = s['nat'], s['X'], s['TAU'], s['cen'], s['m']
    M = np.zeros((3 * nat, 3 * nat), complex)
    fact = [1, 1, 2][n]
    for t in range(nat):
        r = (X - X[cen[t]]) @ qh
        w = (1j * r) ** n / fact
        for tp in range(nat):
            sel = TAU == tp
            M[3 * t:3 * t + 3, 3 * tp:3 * tp + 3] = np.einsum('ajb,j->ab', F[3 * t:3 * t + 3][:, sel], w[sel]) / np.sqrt(m[t] * m[tp])
    return M


def degenerate_sets(w, tol):
    sets, cur = [], [0]
    for i in range(1, len(w)):
        if abs(w[i] - w[cur[-1]]) <= tol * abs(w[i]) + 1e-14:     # relative: acoustic omega^2 ~ q^2 near Gamma
            cur.append(i)
        else:
            sets.append(cur); cur = [i]
    sets.append(cur)
    return sets


def gamma_direct(s, q):
    D = blocks(s, q, 'D'); W = blocks(s, q, 'W')
    D = (D + D.conj().T) / 2; W = (W + W.conj().T) / 2
    w2, e = np.linalg.eigh(D)
    num = np.zeros(len(w2));
    for st in degenerate_sets(w2, 1e-9):
        sub = e[:, st].conj().T @ W @ e[:, st]
        num[st] = np.linalg.eigvalsh((sub + sub.conj().T) / 2)
    return w2, num, -num / (6 * w2)


def gamma_limit(s, qh):
    nat, m = s['nat'], s['m']
    D0, D1, D2 = (taylor(s, qh, 'D', n) for n in range(3))
    W0, W1, W2 = (taylor(s, qh, 'W', n) for n in range(3))
    T = np.zeros((3 * nat, 3))
    for t in range(nat):
        T[3 * t:3 * t + 3] = np.sqrt(m[t] / m.sum()) * np.eye(3)
    P = np.eye(3 * nat) - T @ T.T
    G = np.linalg.pinv(P @ D0 @ P, rcond=1e-10, hermitian=True) if nat > 1 else np.zeros((3, 3))
    checks = dict(D0T=np.abs(D0 @ T).max(), W0T=np.abs(W0 @ T).max(),
                  TD1T=np.abs(T.T @ D1 @ T).max(), TW1T=np.abs(T.T @ W1 @ T).max(),
                  scaleD=np.abs(D0).max(), scaleW=np.abs(W0).max())
    A = T.T @ D2 @ T - T.T @ D1 @ G @ D1 @ T
    B = T.T @ W2 @ T - T.T @ W1 @ G @ D1 @ T - T.T @ D1 @ G @ W1 @ T + T.T @ D1 @ G @ W0 @ G @ D1 @ T
    A = (A + A.conj().T) / 2; B = (B + B.conj().T) / 2
    a, U = np.linalg.eigh(A)
    b = np.zeros(3)
    for st in degenerate_sets(a, 1e-9):
        sub = U[:, st].conj().T @ B @ U[:, st]
        b[st] = np.linalg.eigvalsh((sub + sub.conj().T) / 2)
    return a, b, -b / (6 * a), checks


if __name__ == '__main__':
    mat = sys.argv[1]; relax = len(sys.argv) > 2 and sys.argv[2] == 'relax'
    s = setup(mat, relax)
    print(f'{mat}: {s["nat"]} atoms, cluster {len(s["X"])} sites, W with {"FOCEX relaxation xi" if relax else "no relaxation"}'
          f'; max |xi_trace| {np.abs(s["xi"]).max():.3e} A')
    conv = read_structure(os.path.join(os.path.dirname(os.path.abspath(__file__)), mat, 'structure.params'))[0]
    dirs = {'a1': conv[0], 'a3(c)': conv[2], 'a1+a2+a3': conv.sum(0), 'uinc': UINC}
    for name, v in dirs.items():
        qh = v / np.linalg.norm(v)
        a, b, gl, ck = gamma_limit(s, qh)
        print(f'\n direction {name} {np.round(qh, 4)}: checks |D0 T| {ck["D0T"]:.1e} |W0 T| {ck["W0T"]:.1e} '
              f'|T D1 T| {ck["TD1T"]:.1e} |T W1 T| {ck["TW1T"]:.1e} (scale D {ck["scaleD"]:.2f}, W {ck["scaleW"]:.2f})')
        print(f'   limit: a = omega^2/q^2 {np.array2string(a, precision=8)}  b {np.array2string(b, precision=8)}')
        print(f'   gamma limit {np.array2string(gl, precision=10)}')
        print('    q (1/A)   gamma direct (3 acoustic)                       |gam-lim| max   (w2/q2-a)/q max   (num/q2-b)/q max')
        for q in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 1e-5, 1e-6):
            w2, num, gd = gamma_direct(s, q * qh)
            dev = np.abs(gd[:3] - gl).max()
            print(f'   {q:8.0e}  {np.array2string(gd[:3], precision=10):48s} {dev:12.3e}   '
                  f'{np.abs(w2[:3] / q**2 - a).max() / q:12.3e}   {np.abs(num[:3] / q**2 - b).max() / q:12.3e}')
