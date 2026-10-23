#!/bin/bash
# Tests of poscar2focex.x against pymatgen (needs python3 with numpy, pymatgen; e.g. ~/dftphonon-env/bin/python)
#   ./run_tests.sh [python] [poscar2focex.x]
# Works in a temporary folder; prints one line per structure and ALL OK / ALL INVARIANT at the end.
PY=${1:-python3}
X=${2:-$(command -v poscar2focex.x || echo ~/ALADYN/bin/poscar2focex.x)}
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d)
$PY -I "$HERE/make_structures.py" "$W/tests" || exit 1
$PY -I "$HERE/compare_pymatgen.py" "$X" "$W/tests"
$PY -I "$HERE/invariance.py" "$X" "$W/tests"
echo "test folders in $W"
