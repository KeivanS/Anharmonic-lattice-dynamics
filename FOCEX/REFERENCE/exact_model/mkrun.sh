#!/bin/bash
# mkrun.sh MAT E [HUANG] [SUFFIX]: copy run_MAT inputs to run_MAT/eE[SUFFIX] (a new directory: refuses an
# existing one), set enforce_inv=E, Huang flag, stress.params
m=$1; e=$2; h=${3:-1}
d=run_$m/e$e$4
if [ -e $d ]; then echo "mkrun.sh: $d exists" >&2; exit 1; fi
mkdir $d
for f in structure.params dielectric.params latdyn.params kpbs.params default.params POSCAR1 FORCEDISP1; do cp run_$m/$f $d/; done
sed -i '' "5s/^.*$/1 1 $h $e            # transl, rot, Huang (stress.params), enforce_inv/" $d/structure.params
~/dftphonon-env/bin/python -c "
import json,numpy as np; j=json.load(open('run_$m/setup.json')); s=np.array(j['stress_GPa'])
open('$d/stress.params','w').write('\n'.join(' '.join(f'{v:.12f}' for v in r) for r in s)+'\n')"
echo $d
