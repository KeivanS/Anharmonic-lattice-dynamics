"""mkin.py RUN: lasso_in.bin (A, b, K, x0 from the dumps) for drv.x"""
import sys, numpy as np
d = sys.argv[1]
with open(f'{d}/amat_dump.bin', 'rb') as f:
    nh, ntr, nro, nhu, nf, nc = np.frombuffer(f.read(24), np.int32)
    C = np.frombuffer(f.read(8 * nh * nc), np.float64).reshape(nc, nh).T
    dh = np.frombuffer(f.read(8 * nh), np.float64)
    A = np.frombuffer(f.read(8 * nf * nc), np.float64).reshape(nc, nf).T
    b = np.frombuffer(f.read(8 * nf), np.float64)
U, S, Vt = np.linalg.svd(C, full_matrices=True)
r = np.sum(S > 1e-10 * S[0]); K = Vt[r:].T
x0 = Vt[:r].T @ ((U[:, :r].T @ dh) / S[:r])
with open('lasso_in.bin', 'wb') as f:
    np.array([nf, nc, K.shape[1]], np.int32).tofile(f)
    A.ravel(order='F').tofile(f)
    b.tofile(f); K.ravel(order='F').tofile(f); x0.tofile(f)
print('rows', nf, 'FCs', nc, 'kernel', K.shape[1])
