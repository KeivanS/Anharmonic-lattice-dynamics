"""Exact harmonic phonons of the fake_vasp model (nearest-neighbour LJ fcc) and comparison with
FOCEX's bs_freq.dat:   python exact_phonons.py 05_focex/bs_freq.dat"""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
import fake_vasp as fv
r0 = 2**(1/6) * fv.SIG; a = r0 * np.sqrt(2); m = 39.948          # amu (focex-agent MASS table)
h = 1e-5
d2v = (fv.pair(r0 + h)[1] - fv.pair(r0 - h)[1]) / (2 * h)          # V''(r0), eV/A^2
nn = np.array([[s1, s2, 0] for s1 in (1, -1) for s2 in (1, -1)] +
              [[s1, 0, s2] for s1 in (1, -1) for s2 in (1, -1)] +
              [[0, s1, s2] for s1 in (1, -1) for s2 in (1, -1)]) * a / 2
def freqs(k):
    D = np.zeros((3, 3))
    for R in nn:
        e = R / np.linalg.norm(R)
        D += d2v * (1 - np.cos(k @ R)) * np.outer(e, e)
    w2 = np.linalg.eigvalsh(D / m)
    return np.sort(np.sign(w2) * np.sqrt(np.abs(w2)) * 521.1098918)   # 1/cm, as FOCEX's cnst
rows = {}
for l in open(sys.argv[1]):
    if l.lstrip().startswith('#') or not l.strip(): continue
    x = l.split(); rows.setdefault(int(x[0]), (np.array([float(v) for v in x[3:6]]), []))[1].append(float(x[6]))
err = []; worst = None
for nk, (k, w) in rows.items():
    ex = freqs(k); fx = np.sort(w)
    d = np.abs(fx - ex).max()
    err.append(d)
    if worst is None or d > worst[0]: worst = (d, k, fx, ex)
wmax = max(freqs(k).max() for k, w in rows.values())
print(f'{len(rows)} k-points; max |w_FOCEX - w_exact| = {worst[0]:.3f} 1/cm (at k = {np.round(worst[1], 3)}), '
      f'= {worst[0] / wmax:.2%} of the highest frequency; median {np.median(err):.3f} 1/cm')
print('   FOCEX', np.round(worst[2], 2), '\n   exact', np.round(worst[3], 2))
print(f'highest frequency: FOCEX {max(max(w) for k, w in rows.values()):.2f}, exact {max(freqs(k).max() for k, w in rows.values()):.2f} 1/cm')
