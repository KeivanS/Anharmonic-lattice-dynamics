#!/bin/bash
# run.sh MAT E [HUANG] [SUFFIX]: make run_MAT/eE[SUFFIX] (mkrun.sh) and run FOCEX there with the dumps
# verify.py needs. FOCEX binary: $FOCEX_X (default ~/ALADYN/bin/focex-9.0.x).
cd "$(dirname "$0")"
x=${FOCEX_X:-~/ALADYN/bin/focex-9.0.x}
d=$(./mkrun.sh $1 $2 ${3:-1} $4) || exit 1
cd $d && FOCEX_DUMP_AMAT=3 FOCEX_DUMP_MAP=1 /usr/bin/time -l $x > focex.out 2> time.err
echo "$d exit $?"
