import sys, numpy as np
sys.argv = ['', 'run_ZnO/e3']
d = 'run_ZnO/e3'
xe = np.load(f'{d}/x_exact.npy'); cr = np.load(f'{d}/colrank.npy')
with open(f'{d}/amat_dump.bin', 'rb') as f:
    nh, ntr, nro, nhu, nf, nc = np.frombuffer(f.read(24), np.int32)
    f.seek(24 + 8 * nh * nc + 8 * nh)
    A = np.frombuffer(f.read(8 * nf * nc), np.float64).reshape(nc, nf).T
w = np.sqrt(np.mean(A ** 2, 0)); w = np.maximum(w, 1e-8 * w.max())
def rd(run):
    x = []
    for r in (1, 2, 3, 4):
        for l in open(f'{run}/fc{r}_irr.dat'):
            s = l.split()
            if s and not s[0].startswith('#'): x.append(float(s[2 + 2 * r]))
    return np.array(x)
true = np.abs(xe) > 1e-9 * np.abs(xe).max()
for run in ('run_ZnO/e2', 'run_ZnO/e3'):
    x = rd(run); z = w * x; ze = w * xe
    sel = np.abs(z) > 1e-8 * np.abs(z).max()
    lost = true & ~sel
    print(f'{run}: selected (|z|>1e-8 max) {sel.sum()}, true {true.sum()}, true lost {lost.sum()}, by rank',
          [int(np.sum(lost & (cr == r))) for r in (1, 2, 3, 4)])
    if lost.any():
        i = np.argsort(-np.abs(ze) * lost)[:8]
        print('   largest lost: |z_exact|/max', np.round(np.abs(ze[i]) / np.abs(ze).max(), 5), 'rank', cr[i], 'x_exact', np.round(xe[i], 4))
