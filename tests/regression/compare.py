"""Compare every output file of two regression runs.
    python3 compare.py REF/runs/LABEL_A REF/runs/LABEL_B
Lines with times or dates are ignored; numbers are compared exactly (byte-identical after masking),
and for files that differ the largest relative difference of their numbers is printed.
Two FOCEX outputs print uninitialised memory and are masked: `gr` on the log line
"Initial Pressure(GPa), gr,xagne" and the atomposindx column of supercell_N.xyz."""
import os, sys, re, filecmp
A, B = sys.argv[1], sys.argv[2]
mask = re.compile(r'(launched at|ended at|TIME|time|Time|date|CPU|seconds|  [0-9]+\.[0-9]+ s$)')
skip = {'times.dat', 'timeanal.dat', 'focex.out', 'out', 'build.log'}
num = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?')
bad = 0
for case in sorted(os.listdir(A)):
    pa, pb = os.path.join(A, case), os.path.join(B, case)
    if not os.path.isdir(pb): print(f'{case}: missing in B'); bad += 1; continue
    diffs = []
    for f in sorted(os.listdir(pa)):
        fa, fb = os.path.join(pa, f), os.path.join(pb, f)
        if f in skip or os.path.isdir(fa): continue
        if not os.path.exists(fb): diffs.append(f'{f} (missing)'); continue
        if filecmp.cmp(fa, fb, shallow=False): continue
        def clean(path):
            out = []
            for l in open(path, errors='replace'):
                if mask.search(l) or 'gr,xagne' in l: continue
                if os.path.basename(path).startswith('supercell_'):
                    l = re.sub(r'\)\s*(-?\d+|\*+)', ')', l)          # atomposindx, never set
                out.append(l)
            return out
        la, lb = clean(fa), clean(fb)
        if la == lb: continue
        xa = [float(x.replace('D','E')) for x in num.findall(''.join(la))]
        xb = [float(x.replace('D','E')) for x in num.findall(''.join(lb))]
        if len(xa) != len(xb): diffs.append(f'{f} (layout)'); continue
        sc = max([abs(x) for x in xa] + [1e-300])
        m = max(abs(x - y) / max(abs(x), 1e-6 * sc) for x, y in zip(xa, xb)) if xa else 0
        diffs.append(f'{f} ({m:.1e})')
    extra = [f for f in os.listdir(pb) if f not in os.listdir(pa) and f not in skip]
    if extra: diffs.append('only in B: ' + ' '.join(sorted(extra)))
    print(f'{case:24s} ' + ('identical' if not diffs else 'DIFFER: ' + ', '.join(diffs)))
    bad += bool(diffs)
print('ALL IDENTICAL' if not bad else f'{bad} case(s) differ')
