"""compare.py BASE RUN...: irreducible FCs (fcN_irr.dat, N = 1..4) of each RUN against BASE, term by term.
Rank-N lines: index, group, N pairs (atom, direction), value, ... (the value is column 2 + 2N).
Prints per rank: max |x_run - x_base| and that over max |x_base|. fcN_irr.dat carry ~7 significant digits."""
import sys, os
import numpy as np


def read(d, n):
    f = os.path.join(d, f'fc{n}_irr.dat')
    if not os.path.exists(f):
        return None
    v = []
    for l in open(f):
        w = l.split()
        if len(w) > 2 + 2 * n and w[0].isdigit():
            v.append(float(w[2 + 2 * n]))
    return np.array(v)


base = sys.argv[1]
print(f'{"run":28s}' + ''.join(f'   FC{n}: max|dx|   rel  ' for n in range(1, 5)))
for d in sys.argv[2:]:
    row = f'{os.path.basename(d.rstrip("/")):28s}'
    for n in range(1, 5):
        a, b = read(base, n), read(d, n)
        if a is None or b is None or len(a) != len(b):
            row += f'   {"n/a":>22s}'
            continue
        dx = np.abs(a - b).max()
        row += f'   {dx:10.3e} {dx / max(np.abs(a).max(), 1e-30):9.2e}'
    print(row)
