"""huang.py RUN: test S_abcd - S_cdab = -2 Omega (d_ab sig_cd - d_cd sig_ab) on the exact FC2,
S_abcd = sum_{tau in cell} sum_j Phi(tau a, j b) R_cd, R = X_j - X_tau (unique image: range < inradius)"""
import sys, json, itertools, numpy as np
sys.path.insert(0, '.')
from model3 import Model
run = sys.argv[1]; info = json.load(open(f'{run}/setup.json'))
Xi, latx = np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy')
m = Model(info['rL'], info['Ra']).build(Xi, latx)
N = len(Xi); H = m.hessian(Xi).reshape(N, 3, N, 3)
V = abs(np.linalg.det(latx))
nprim = {'CaF2': 3, 'ZnO': 4, 'HfO2': 12}[info['mat']]
ncell = N // nprim; Om = V / ncell
sig = m.stress(Xi, V)
# one representative per sublattice: sum over all atoms / ncell
inv = np.linalg.inv(latx)
S = np.zeros((3, 3, 3, 3))
for t in range(N):
    R = Xi - Xi[t]; f = R @ inv; R = (f - np.rint(f)) @ latx
    S += np.einsum('ajb,jc,jd->abcd', H[t], R, R)
S /= ncell
d = np.eye(3)
lhs = S - S.transpose(2, 3, 0, 1)
rhs = -2 * Om * (np.einsum('ab,cd->abcd', d, sig) - np.einsum('cd,ab->abcd', d, sig))
print(f"{run}: max|S| {np.abs(S).max():.3f}, max|S_abcd - S_cdab| {np.abs(lhs).max():.4e}, "
      f"max|rhs| {np.abs(rhs).max():.4e}, max|lhs - rhs| {np.abs(lhs - rhs).max():.3e}, max|lhs + rhs| {np.abs(lhs + rhs).max():.3e}, max FC1 {np.abs(m.forces(Xi)).max():.3f}")
