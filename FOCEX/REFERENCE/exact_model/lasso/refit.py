"""refit.py RUN: FOCEX's refit (Gram matrix, eigenvalue cut) vs a QR/SVD refit on the rows themselves,
on the dumped invariance and force rows, supports = true support + random extras."""
import sys, os, numpy as np
d = sys.argv[1]
xe = np.load(f'{d}/x_exact.npy')
with open(f'{d}/amat_dump.bin', 'rb') as f:
    nh, ntr, nro, nhu, nf, nc = np.frombuffer(f.read(24), np.int32)
    C = np.frombuffer(f.read(8 * nh * nc), np.float64).reshape(nc, nh).T
    dh = np.frombuffer(f.read(8 * nh), np.float64)
    A = np.frombuffer(f.read(8 * nf * nc), np.float64).reshape(nc, nf).T
    b = np.frombuffer(f.read(8 * nf), np.float64)
U, S, Vt = np.linalg.svd(C, full_matrices=True)
r = np.sum(S > 1e-10 * S[0]); K = Vt[r:].T
x0 = Vt[:r].T @ ((U[:, :r].T @ dh) / S[:r])
print(f'rows {nf}, FCs {nc}, kernel {K.shape[1]}, |C x0 - d| {np.abs(C @ x0 - dh).max():.1e}, |C xe - d| {np.abs(C @ xe - dh).max():.1e}')
w = np.sqrt(np.mean(A ** 2, 0)); w = np.maximum(w, 1e-8 * w.max())
Q = np.linalg.svd(w[:, None] * K, full_matrices=False)[0]
N = len(b); At = (A / w) @ Q / np.sqrt(N); q = (b - A @ x0) / np.sqrt(N)
G = At.T @ At; c = At.T @ q; z0 = w * x0
te = Q.T @ (w * (xe - x0))          # exact t
print(f'exact t reproduces z: {np.abs(Q @ te + z0 - w * xe).max():.1e}; |At te - q| {np.linalg.norm(At @ te - q):.1e}')
sv = np.linalg.svd(At, compute_uv=False); print(f'sing. values of At: max {sv[0]:.3e} min {sv[-1]:.3e}  cond {sv[0]/sv[-1]:.2e}')
true = np.abs(xe) > 1e-9 * np.abs(xe).max()
rng = np.random.default_rng(1)
def nullsp(B):
    if B.shape[0] == 0: return np.eye(B.shape[1])
    u, s, vt = np.linalg.svd(B); rk = np.sum(s > 1e-9 * max(s[0], 1e-300)); return vt[rk:].T
for extra in (0, 300, 600, 900):
    sel = true.copy(); off = np.where(~sel)[0]; sel[rng.choice(off, min(extra, len(off)), replace=False)] = True
    B = Q[~sel]; Nb = nullsp(B)
    tp = np.linalg.lstsq(B, -z0[~sel], rcond=None)[0] if (~sel).any() else np.zeros(Q.shape[1])
    # FOCEX: eigen-cut 1e-13 on N^T G N
    H = Nb.T @ G @ Nb; hv = Nb.T @ (c - G @ tp); ev, V = np.linalg.eigh(H)
    keep = ev > 1e-13 * ev.max(); s_g = V[:, keep] @ ((V[:, keep].T @ hv) / ev[keep])
    # rows: least squares on At N directly
    s_r = np.linalg.lstsq(At @ Nb, q - At @ tp, rcond=1e-14)[0]
    for name, s in (('Gram 1e-13', s_g), ('rows lstsq', s_r)):
        t = tp + Nb @ s; x = x0 + (Q @ t) / w
        print(f'support {sel.sum():5d} (true {true.sum()}), null dim {Nb.shape[1]}, H eig dropped {np.sum(~keep)}: {name}: '
              f'max|x - exact| {np.abs(x - xe).max():.2e}, resid {np.linalg.norm(At @ t - q)/np.linalg.norm(q):.1e}')
