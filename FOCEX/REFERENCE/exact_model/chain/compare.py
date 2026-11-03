"""compare.py (in the run_chains.sh work directory): FOCEX conductance.dat against ref.py, channels per cell M*a^2."""
import numpy as np
for n, a, lo, hi in (('B', 8.4, 3, 51.0), ('C', 3.5, 3, 46), ('D', 3.5, 3, 27), ('C123', 3.5, 3, 46)):
    ref = np.loadtxt(f'ref_{n}.txt', skiprows=2); f = np.loadtxt(f'run_{n}/conductance.dat', comments='#')
    e = np.interp(f[:, 0], ref[:, 0], ref[:, 1]); m = (f[:, 0] > lo) & (f[:, 0] < hi)
    print(f'{n}: rms|FOCEX-exact|/mean(exact) = {np.sqrt(((f[m,1]-e[m])**2).mean())/e[m].mean():.4f} for {lo} < om < {hi} 1/cm')
    for w in np.linspace(lo, hi, 12):
        i = np.argmin(abs(f[:, 0]-w)); print(f'   om {f[i,0]:7.3f}  M*a^2 FOCEX {f[i,1]*a*a:7.4f}  exact {e[i]*a*a:7.4f}')
# integer plateaus (channels per cell; exact 3 and 6): mean, min, max of FOCEX's M*a^2
for n, a, lo, hi in (('B', 8.4, 8, 48.5), ('B', 8.4, 49.8, 50.6), ('D', 3.5, 5, 26)):
    f = np.loadtxt(f'run_{n}/conductance.dat', comments='#'); m = (f[:, 0] > lo) & (f[:, 0] < hi); v = f[m, 1]*a*a
    print(f'plateau {n} {lo}-{hi} 1/cm: M*a^2 mean {v.mean():.4f}, min {v.min():.4f}, max {v.max():.4f}')
