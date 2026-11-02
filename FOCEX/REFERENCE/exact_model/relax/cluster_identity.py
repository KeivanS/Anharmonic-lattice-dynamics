"""finite cluster of the exact model: E2 of u = eps x from energies vs 1/2 u Phi u vs -1/4 sum (du) Phi (du)"""
import sys, os, json, numpy as np
EM = os.path.expanduser('~/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/REFERENCE/exact_model')
sys.path.insert(0, EM)
from model3 import Model, PAR
run = f'{EM}/run_ZnO'; info = json.load(open(f'{run}/setup.json'))
Xi = np.load(f'{run}/Xi.npy')
X = Xi[:40]                                    # a finite cluster
m = Model(info['rL'], info['Ra'], dict(PAR)).build(X, None)
N = len(X); H = m.hessian(X).reshape(N, 3, N, 3); f = m.forces(X)
for ev in (np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.]]), np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0.]])):
    u = X @ ev.T; h = 1e-4
    E2 = (m.energy(X + h * u) - 2 * m.energy(X) + m.energy(X - h * u)) / h**2 / 2
    a = 0.5 * np.einsum('ia,iajb,jb', u, H, u)
    du = u[None, :, :] - u[:, None, :]
    b = -0.25 * np.einsum('ija,iajb,ijb', du, H, du)
    print('E2 energies', E2, ' 1/2 uPhiu', a, ' -1/4 du Phi du', b)
print('ASR: max|sum_b Phi_ab|', np.abs(H.sum(2)).max(), ' max|Phi|', np.abs(H).max())
print('symmetry: max|Phi_ab^ag - Phi_ba^ga|', np.abs(H - H.transpose(2, 3, 0, 1)).max())
print('max |force| in the cluster', np.abs(f).max())
Hsym = 0.5 * (H + H.transpose(0, 3, 2, 1))       # each 3x3 block Phi_ij symmetrized in its Cartesian indices
print('max |antisymmetric part of the 3x3 blocks Phi_ij^{ab}-Phi_ij^{ba}|/2', np.abs(H - Hsym).max())
for ev in (np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.]]), np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0.]])):
    u = X @ ev.T; du = u[None, :, :] - u[:, None, :]
    print('1/2 u Phi_sym u', 0.5 * np.einsum('ia,iajb,jb', u, Hsym, u), ' -1/4 du Phi du', -0.25 * np.einsum('ija,iajb,ijb', du, H, du))
