"""gruneisen_focex_check.py MAT RUNDIR: FOCEX's bs_grun.dat (run on the cubic-truncated exact model, mkcubic.py)
against the same quantity from the exact force constants (gruneisen_q0.setup, W with FOCEX's internal
relaxation), at FOCEX's own band-path points (make_kp_bs, kpbs.params) and with FOCEX's rules (freq_vel):
omega^2 at k, eigenvectors at k + 3e-4 gsize uinc (at Gamma: omega^2 at 1e-6 dir, eigenvectors at
3e-4 gsize dir, dir = the path segment, W at q = 0), gamma = -<e|W(k)|e>/(6 omega^2), |gamma| > 1e3 -> 0.
At Gamma the 3 rigid translations are compared with the long-wave limit (gruneisen_q0.gamma_limit).
bs_grun.dat carries 4 decimals; its band columns follow FOCEX's band connection, so bands are matched by
frequency."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gruneisen_q0 import setup, blocks, gamma_limit, UINC
from model import read_structure

mat, rd = sys.argv[1], sys.argv[2]
s = setup(mat, relax=True)
conv, prim, atoms = read_structure(os.path.join(rd, 'structure.params'))
G = 2 * np.pi * np.linalg.inv(prim).T                    # rows g01, g02, g03 (Cartesian, 1/A)
gsize = abs(np.linalg.det(G)) ** (1 / 3)
L = [l.split('#')[0].split() for l in open(os.path.join(rd, 'kpbs.params'))]
units, nkdir, ndir = int(L[0][0]), int(L[1][0]), int(L[2][0])
assert units == 1
ext = np.array([[float(v) for v in L[3 + i][1:4]] for i in range(ndir + 1)]) @ G
kp = np.array([ext[i] + j * (ext[i + 1] - ext[i]) / (nkdir - 1 + 1e-8) for i in range(ndir) for j in range(nkdir)])
ndn = 3 * s['nat']

tok = [float(x) for l in open(os.path.join(rd, 'bs_grun.dat')).read().split('\n')[1:] for x in l.split()]
rec = 1 + 4 + 2 * ndn
fx = np.array(tok).reshape(-1, rec)
assert len(fx) == len(kp), (len(fx), len(kp))
print(f'{mat}: {len(kp)} band-path points; max |k(FOCEX, 4 decimals) - k(here)| = {np.abs(fx[:, 2:5] - kp).max():.1e} 1/A')
gfx = fx[:, 5 + ndn:]


def herm(M):
    return (M + M.conj().T) / 2


def blocks_R(q, kind):
    """D or W with the phase exp(iq.R) of FOCEX's set_dynamical_matrix_std (blocks uses x_j - x_t)"""
    F = s['Phi'] if kind == 'D' else s['W']
    nat, X, TAU, cen, m = s['nat'], s['X'], s['TAU'], s['cen'], s['m']
    M = np.zeros((3 * nat, 3 * nat), complex)
    ph = np.exp(1j * (X - X[np.array(cen)[TAU]]) @ q)
    for t in range(nat):
        for tp in range(nat):
            sel = TAU == tp
            M[3 * t:3 * t + 3, 3 * tp:3 * tp + 3] = np.einsum('ajb,j->ab', F[3 * t:3 * t + 3][:, sel], ph[sel]) / np.sqrt(m[t] * m[tp])
    return M


PHASE = 'x'


def focex_gamma(i):
    bl = blocks if PHASE == 'x' else (lambda s_, q, kind: blocks_R(q, kind))
    k = kp[i]
    if np.linalg.norm(k) < 1e-12:
        d = kp[i + 1] - kp[i] if i + 1 < len(kp) else kp[i] - kp[i - 1]
        d = d / np.linalg.norm(d)
        q, shift = 1e-6 * d, d
    else:
        q, shift = k, UINC
    w2 = np.linalg.eigvalsh(herm(bl(s, q, 'D')))
    _, e = np.linalg.eigh(herm(bl(s, q + 3e-4 * gsize * shift / np.linalg.norm(shift), 'D')))
    W = herm(bl(s, k, 'W'))
    num = np.real(np.einsum('ia,ij,ja->a', e.conj(), W, e))
    g = -num / (6 * w2)
    if np.linalg.norm(k) < 1e-12:
        g[:3] = np.nan      # translations: limit below
    g[np.abs(g) > 1e3] = 0
    return g, w2


# bs_grun.dat lists the bands in FOCEX's band-connection order (band_sort_bs), omega and gamma columns
# together: match them to the energy-sorted exact bands by frequency; inside a group of equal frequencies
# (|d omega| < 0.01 1/cm) gamma is compared as a sorted set (basis-dependent there)
wfx = fx[:, 5:5 + ndn]
cnst = np.median(wfx[1] / np.sqrt(np.abs(np.linalg.eigvalsh(herm(blocks(s, kp[1], 'D'))))))
def compare(label):
    dev, worst, maxw = [], (0, None), 0
    for i in range(len(kp)):
        g, w2 = focex_gamma(i)
        w = np.sqrt(np.abs(w2)) * cnst
        o = np.argsort(wfx[i], kind='stable')
        wf, gf = wfx[i, o], gfx[i, o]
        maxw = max(maxw, np.abs(wf - w).max())
        d = np.zeros(ndn)
        j = 0
        while j < ndn:
            k = j
            while k + 1 < ndn and w[k + 1] - w[j] < 1e-2:
                k += 1
            ge, gx = g[j:k + 1], gf[j:k + 1]
            ok = ~np.isnan(ge)
            if ok.any():
                d[j:k + 1][ok] = np.abs(np.sort(ge[ok]) - np.sort(gx[ok]))
            j = k + 1
        dev.append(d.max())
        if d.max() > worst[0]:
            worst = (d.max(), i, int(np.argmax(d)))
    dev = np.array(dev)
    if PHASE == 'x': print(f'omega (bands matched by frequency): max |FOCEX - exact| = {maxw:.4f} 1/cm (conversion {cnst:.4f})')
    print(f'{label}: gamma, all modes except the Gamma translations: max |FOCEX - exact| = {dev.max():.2e} '
          f'(k-point {worst[1] + 1}, sorted band {worst[2] + 1}); points with max dev > 2e-4: {(dev > 2e-4).sum()} of {len(dev)}')



compare('exact formula (phase x_j - x_t)')
PHASE = 'R'
compare("FOCEX's phase exp(ik.R) (with eigenvectors at k + 3e-4|G|u the two differ near Gamma)")
i0 = int(np.where(np.linalg.norm(kp, axis=1) < 1e-12)[0][0])
d = (kp[i0 + 1] - kp[i0]); d /= np.linalg.norm(d)
a, b, gl, ck = gamma_limit(s, d)
print(f'Gamma (point {i0 + 1}), direction {np.round(d, 4)}: long-wave limit of the 3 translations {np.round(gl, 6)}; '
      f'FOCEX writes {gfx[i0, :3]}')
print(f'  next path point (|k| = {np.linalg.norm(kp[i0 + 1]):.4f} 1/A): FOCEX acoustic gamma {gfx[i0 + 1, :3]}')
