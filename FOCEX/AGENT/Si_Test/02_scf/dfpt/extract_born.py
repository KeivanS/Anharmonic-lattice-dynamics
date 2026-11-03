#!/usr/bin/env python3
"""extract_born.py — Extract Born charges, dielectric and piezoelectric tensors
from a VASP DFPT (IBRION=8 + LEPSILON) or finite-field (LCALCEPS) OUTCAR.
Writes BORN (phonopy format, uses electronic eps_inf) and born_charges.txt.
Usage: python3 extract_born.py [OUTCAR]
"""
import sys, os, re

outcar = sys.argv[1] if len(sys.argv) > 1 else 'OUTCAR'
try:
    text = open(outcar).read()
except FileNotFoundError:
    print(f"ERROR: {outcar} not found"); sys.exit(1)

def tensor3(segment):
    """First 3x3 float block after a dashed separator line."""
    m = re.search(
        r'\n\s*-+\s*\n'
        r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\n'
        r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\n'
        r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', segment)
    if not m:
        return None
    v = [float(x) for x in m.groups()]
    return [v[0:3], v[3:6], v[6:9]]

def find_eps(ionic):
    """Dielectric tensor blocks; ionic=False -> electronic eps_inf (prefer the
    variant including local-field effects), ionic=True -> ionic contribution."""
    best = None
    for m in re.finditer(r'MACROSCOPIC STATIC DIELECTRIC TENSOR[^\n]*', text):
        hdr = m.group()
        if ('IONIC' in hdr) != ionic:
            continue
        t = tensor3(text[m.end():m.end() + 400])
        if t and (best is None or 'local field' in hdr):
            best = t
    return best

def find_piezo(ionic):
    """Piezoelectric stress tensor e (C/m^2): 3 field rows x 6 strain columns."""
    best = None
    for m in re.finditer(r'PIEZOELECTRIC TENSOR[^\n]*', text):
        if ('IONIC' in m.group()) != ionic:
            continue
        rows = re.findall(
            r'\n\s*[xyz]\s+' + r'([-\d.]+)\s+' * 5 + r'([-\d.]+)',
            text[m.end():m.end() + 700])
        if len(rows) >= 3:
            best = [[float(x) for x in r] for r in rows[:3]]
    return best

def add3(a, b):
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]

eps_el = find_eps(ionic=False)
if not eps_el:
    print("ERROR: dielectric tensor not found — did LEPSILON/LCALCEPS run?"); sys.exit(1)
eps_ion = find_eps(ionic=True)
eps_tot = add3(eps_el, eps_ion) if eps_ion else None

# Born effective charges — take the last block (cumulative output)
born_blocks = list(re.finditer(r'BORN EFFECTIVE CHARGES.*?(?=BORN EFFECTIVE CHARGES|\Z)', text, re.DOTALL))
if not born_blocks:
    print("ERROR: Born charges not found in OUTCAR"); sys.exit(1)
ions = re.findall(
    r'ion\s+\d+\s*\n'
    r'\s*\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\n'
    r'\s*\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\n'
    r'\s*\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)',
    born_blocks[-1].group())
if not ions:
    print("ERROR: could not parse Born charge blocks"); sys.exit(1)
born = [[[float(x) for x in m[i*3:(i+1)*3]] for i in range(3)] for m in ions]

# Element label per ion, from the POSCAR next to the OUTCAR
def ion_labels(n):
    try:
        lines = open(os.path.join(os.path.dirname(os.path.abspath(outcar)),
                                  'POSCAR')).read().splitlines()
        syms, counts = lines[5].split(), [int(x) for x in lines[6].split()]
        out = [s for s, c in zip(syms, counts) for _ in range(c)]
        if len(out) == n:
            return out
    except Exception:
        pass
    return [''] * n

labels = ion_labels(len(born))

# Acoustic sum rule: Born charges must sum to zero over all ions
asr = [[sum(z[i][j] for z in born) for j in range(3)] for i in range(3)]
asr_max = max(abs(x) for row in asr for x in row)

pz_el  = find_piezo(ionic=False)
pz_ion = find_piezo(ionic=True)
pz_tot = None
if pz_el:
    pz_tot = ([[pz_el[i][j] + pz_ion[i][j] for j in range(6)] for i in range(3)]
              if pz_ion else pz_el)

# Write BORN (phonopy format: electronic eps_inf + Z*)
with open('BORN', 'w') as f:
    f.write("# Born effective charges and dielectric tensor from VASP OUTCAR\n")
    f.write("14.400\n")
    for row in eps_el:
        f.write("  " + "  ".join(f"{x:12.6f}" for x in row) + "\n")
    for z in born:
        for row in z:
            f.write("  " + "  ".join(f"{x:12.6f}" for x in row) + "\n")
print(f"Written BORN ({len(born)} atoms)")

def w_tensor(f, t, indent='  ', ncol=3):
    for row in t:
        f.write(indent + "  ".join(f"{x:10.5f}" for x in row[:ncol]) + "\n")

with open('born_charges.txt', 'w') as f:
    f.write("Born Effective Charges, Dielectric and Piezoelectric Tensors\n")
    f.write("=" * 62 + "\n\n")
    f.write("Electronic (ion-clamped) dielectric tensor eps_inf:\n")
    w_tensor(f, eps_el)
    if eps_ion:
        f.write("\nIonic contribution to the dielectric tensor:\n")
        w_tensor(f, eps_ion)
        f.write("\nTotal static dielectric tensor eps_0 = eps_inf + ionic:\n")
        w_tensor(f, eps_tot)
        f.write(f"\n  eps_0 (diagonal): {eps_tot[0][0]:.4f}  "
                f"{eps_tot[1][1]:.4f}  {eps_tot[2][2]:.4f}\n")
    else:
        f.write("\n(no ionic contribution found in OUTCAR — eps_0 unavailable)\n")
    f.write("\nBorn effective charge tensors (Z*):\n")
    for i, z in enumerate(born):
        f.write(f"\n  Ion {i+1} {labels[i]}:\n")
        w_tensor(f, z, indent='    ')
    f.write(f"\nAcoustic sum rule: max |sum_ions Z*| component = {asr_max:.5f} e\n")
    f.write("  (should be ~0; large values indicate under-converged k-mesh/EDIFF)\n")
    if pz_el:
        f.write("\nPiezoelectric stress tensor e (C/m^2), rows Ex Ey Ez,\n"
                "columns XX YY ZZ XY YZ ZX:\n")
        f.write("\n  Electronic (clamped-ion):\n")
        w_tensor(f, pz_el, indent='    ', ncol=6)
        if pz_ion:
            f.write("\n  Ionic contribution:\n")
            w_tensor(f, pz_ion, indent='    ', ncol=6)
            f.write("\n  Total (relaxed-ion):\n")
            w_tensor(f, pz_tot, indent='    ', ncol=6)
        f.write("\n  (meaningful only for non-centrosymmetric crystals)\n")
print("Written born_charges.txt")

diag = eps_tot or eps_el
name = 'eps_0' if eps_tot else 'eps_inf'
print(f"\n{name} (diagonal): {diag[0][0]:.4f}  {diag[1][1]:.4f}  {diag[2][2]:.4f}")
print(f"Acoustic sum rule violation (max component): {asr_max:.5f} e")
print("Born charges (diagonal elements):")
for i, z in enumerate(born):
    print(f"  Ion {i+1} {labels[i]}: Zxx={z[0][0]:8.4f}  Zyy={z[1][1]:8.4f}  Zzz={z[2][2]:8.4f}")
