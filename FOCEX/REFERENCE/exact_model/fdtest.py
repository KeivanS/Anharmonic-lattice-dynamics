"""fdtest.py RUN: analytic derivatives of model3 against finite differences on the supercell.
FC1 (forces) vs -dE/dx and FC2 (Hessian) vs -dF/dx at a displaced snapshot (central, h=1e-4: error
h^2/6 FC3/FC4 ~ 1e-9); FC3, FC4 at the reference vs central differences (step 0.01 A) of the analytic Hessian, which
is quadratic in x, so those differences are exact up to round-off."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from model3 import Model
run = sys.argv[1]
info = json.load(open(f'{run}/setup.json'))
Xi, latx, U = np.load(f'{run}/Xi.npy'), np.load(f'{run}/latx.npy'), np.load(f'{run}/U.npy')
m = Model(info['rL'], info['Ra']).build(Xi, latx)
N = len(Xi); rng = np.random.default_rng(0)
x = Xi + U[0]
h = 1e-4
# FC1: forces vs energy differences, all 3N
F = m.forces(x); Ffd = np.zeros_like(F)
for i in range(N):
    for a in range(3):
        xp, xm = x.copy(), x.copy(); xp[i, a] += h; xm[i, a] -= h
        Ffd[i, a] = -(m.energy(xp) - m.energy(xm)) / (2 * h)
print(f'{run}: forces  analytic vs -dE/dx (FD): max diff {np.abs(F - Ffd).max():.2e}, max |F| {np.abs(F).max():.3f} eV/A')
# FC2: Hessian vs force differences, all (3N)^2
H = m.hessian(x); Hfd = np.zeros_like(H)
for i in range(N):
    for a in range(3):
        xp, xm = x.copy(), x.copy(); xp[i, a] += h; xm[i, a] -= h
        Hfd[:, 3 * i + a] = -(m.forces(xp) - m.forces(xm)).ravel() / (2 * h)
print(f'{run}: FC2     analytic vs -dF/dx (FD): max diff {np.abs(H - Hfd).max():.2e}, max |FC2| {np.abs(H).max():.3f} eV/A^2, '
      f'asymmetry {np.abs(H - H.T).max():.1e}')
# FC2, FC3, FC4 at the reference from model.fc against the Hessian and its differences
H0 = m.hessian(Xi)
pairs = np.array([(m.pi[p], m.pj[p]) for p in range(len(m.pi))])
def rnd_tuple(r):
    # atoms of one term (triplet/quadruplet/pair) so that the FC is non-zero in general
    if r == 4 and len(m.quad) and rng.random() < 0.7:
        at = np.unique(np.concatenate([pairs[p] for p in m.quad[rng.integers(len(m.quad))]]))
    elif r >= 3 and len(m.trip) and rng.random() < 0.7:
        at = np.unique(np.concatenate([pairs[p] for p in m.trip[rng.integers(len(m.trip))]]))
    else:
        at = pairs[rng.integers(len(pairs))]
    return tuple((int(rng.choice(at)), int(rng.integers(3))) for _ in range(r))
h3 = 0.01   # the Hessian is quadratic in x: central differences exact for any step
def hd(tup, sh):
    h = h3
    xx = Xi.copy()
    for (i, a), c in zip(tup, sh):
        xx[i, a] += c * h
    return m.hessian(xx)
for r in (2, 3, 4):
    err, big = 0.0, 0.0
    for _ in range(60):
        t = rnd_tuple(r)
        while len({a for a, _ in t}) < min(r, 2) and r > 2:
            t = rnd_tuple(r)
        an = m.fc(Xi, t)
        (i, a), (j, b) = t[0], t[1]
        if r == 2:
            fd = H0[3 * i + a, 3 * j + b]
        elif r == 3:
            fd = (hd(t[2:], [1]) - hd(t[2:], [-1]))[3 * i + a, 3 * j + b] / (2 * h3)
        else:
            fd = (hd(t[2:], [1, 1]) - hd(t[2:], [1, -1]) - hd(t[2:], [-1, 1]) + hd(t[2:], [-1, -1]))[3 * i + a, 3 * j + b] / (4 * h3 * h3)
        err = max(err, abs(an - fd)); big = max(big, abs(an))
    print(f'{run}: FC{r} at the reference, 60 random terms: model.fc vs differences of the Hessian: max diff {err:.2e}, max |FC| {big:.3f}')
