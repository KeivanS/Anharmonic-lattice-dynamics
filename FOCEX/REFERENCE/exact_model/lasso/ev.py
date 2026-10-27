"""ev.py RUN: lasso_out.bin against x_exact"""
import sys, numpy as np
d = sys.argv[1]
xe = np.load(f'{d}/x_exact.npy'); cr = np.load(f'{d}/colrank.npy')
x = np.fromfile('lasso_out.bin', np.float64)
true = np.abs(xe) > 1e-9 * np.abs(xe).max()
with open(f'{d}/amat_dump.bin', 'rb') as f:
    nh, ntr, nro, nhu, nf, nc = np.frombuffer(f.read(24), np.int32); f.seek(24 + 8 * nh * nc + 8 * nh)
    A = np.frombuffer(f.read(8 * nf * nc), np.float64).reshape(nc, nf).T
w = np.sqrt(np.mean(A ** 2, 0)); w = np.maximum(w, 1e-8 * w.max()); z = w * x
sel = np.abs(z) > 1e-8 * np.abs(z).max()
print('  true lost', int(np.sum(true & ~sel)), 'kept', int(sel.sum()), 'true', int(true.sum()),
      ' max err by rank', ['%.2e' % np.abs(x - xe)[cr == r].max() for r in (1, 2, 3, 4)])
