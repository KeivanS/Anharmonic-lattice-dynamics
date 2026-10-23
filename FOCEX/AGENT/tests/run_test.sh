#!/bin/bash
# End-to-end test of focex-agent.py without VASP: a nearest-neighbour Lennard-Jones fcc crystal
# (fake_vasp.py) stands in for VASP; FOCEX's phonons are compared with the exact ones.
#   ./run_test.sh [workdir]
# Needs: ~/dftphonon-env/bin/python (numpy, spglib, pymatgen), VASP-Flow, sc_snaps.x (built with -O2, see
# SC_SNAPS; one built with -fcheck=all aborts), read_outcar.x, a FOCEX binary (FOCEX_X below).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
PY=${PY:-$HOME/dftphonon-env/bin/python}
W=${1:-$(mktemp -d)}
FOCEX_X=${FOCEX_X:-$(ls -t "$HOME"/ALADYN/bin/focex-[0-9]*.x "$HOME"/BIN/v[0-9]* 2>/dev/null | head -1)}
SC_SNAPS_X=${SC_SNAPS_X:-$HOME/BIN/sc_snaps.x}
READ_OUTCAR_X=${READ_OUTCAR_X:-$HOME/ALADYN/bin/read_outcar.x}
mkdir -p "$W"; cd "$W"
printf '#!/bin/bash\nexec %s -I %s "$@"\n' "$PY" "$HERE/fake_vasp.py" > fake_vasp; chmod +x fake_vasp
# a skewed, rotated primitive cell of the equilibrium fcc lattice (tests the cell reduction)
"$PY" -I - <<'EOF'
import numpy as np
a = 2**(1/6) * 3.40 * np.sqrt(2)
P = np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * a
M = np.array([[1, 1, 0], [0, 1, 0], [1, 1, 1]])
rng = np.random.default_rng(5); q, _ = np.linalg.qr(rng.normal(size=(3, 3))); q = q if np.linalg.det(q) > 0 else -q
A = (M @ P) @ q.T
open('POSCAR', 'w').write('Ar fcc (nearest-neighbour LJ), skewed and rotated\n1.0\n' +
    ''.join('  %.12f %.12f %.12f\n' % tuple(v) for v in A) + 'Ar\n1\nDirect\n  0.1 0.2 0.3\n')
EOF
cat > focex-agent-instructions.txt <<EOF
Project: ArLJ
RANKS: 1 1 1 0 0 0 0 0
RANGE: 6 2 0 0 0 0 0
FUNCTIONAL: PBE
MPI: 1
DFPT: off
FOCEX_X: $FOCEX_X
SC_SNAPS_X: $SC_SNAPS_X
READ_OUTCAR_X: $READ_OUTCAR_X
EOF
export VASP_STD="$W/fake_vasp" MPI_NP=1
"$PY" "$HERE/../focex-agent.py" --profile workstation --force   # focex-agent-instructions.txt, POSCAR
printf '   TITEL  = PAW_PBE Ar 07Sep2000 (dummy for the fake-VASP test)\n   ENMAX  =  266.000; ENMIN  =  199.500 eV\n   ZVAL   =    8.000    mass and valenz\n' > ArLJ/POTCAR
"$PY" "$HERE/../focex-agent.py" run "$W/ArLJ" || echo "focex-agent run returned $? (see $W/ArLJ/report.log)"
grep -E "=== .* -- (OK|FAILED)" ArLJ/report.log | sed 's/^\[[^]]*\] //' | cut -c1-160
"$PY" -I "$HERE/exact_phonons.py" ArLJ/05_focex/bs_freq.dat
echo "test folder: $W"
