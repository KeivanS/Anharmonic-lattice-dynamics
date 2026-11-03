#!/bin/bash
# run_chains.sh WORKDIR: chain-conductance test of FOCEX (README.md). FOCEX binary: $FOCEX_X (default ~/ALADYN/bin/focex-9.5.x)
set -e
here=$(cd "$(dirname "$0")" && pwd); x=${FOCEX_X:-~/ALADYN/bin/focex-9.5.x}; py=${PY:-~/myenv/bin/python}
mkdir -p "$1" && cd "$1"
#      name  a(A)  rL  normal  mesh      omega_max  reference mesh
for c in "B 8.4 8.6 0,0,1 8_8_96 55 40_4000" "C 3.5 4.6 0,0,1 12_12_48 52 160_1200" \
         "D 3.5 3.6 0,0,1 8_8_48 30 160_1200" "C123 3.5 4.6 1,2,3 12_12_48 52 160_1200"; do
  set -- $c; n=$1; m=${n:0:1}
  if [ ! -d data_$m ]; then           # C123: the data of C with another normal
    mkdir -p Chain$m && printf "Ar chain\n1.0\n $2 0 0\n 0 $2 0\n 0 0 2.5\nAr\n1\nCartesian\n 0 0 0\n" > Chain$m/POSCAR
    (cd Chain$m && ~/ALADYN/bin/poscar2focex.x POSCAR > p2f.log)
    $py "$here/setup_chain.py" Chain$m data_$m 40 > data_$m.setup.log
  fi
  mkdir -p run_$n && for f in structure.params dielectric.params latdyn.params kpbs.params default.params POSCAR1 FORCEDISP1; do cp data_$m/$f run_$n/; done
  sed -i '' "1s/^[0-9 ]*#/${5//_/ }  #/; 2s/0 0 1 /${4//,/ } /; s/OMMAX/$6/" run_$n/latdyn.params
  (cd run_$n && $x > focex.out 2>&1)
  $py "$here/ref.py" $2 $3 ${7//_/ } $4 > ref_$n.txt
done
$py "$here/compare.py"
