#!/bin/bash
# run.sh MAT E TAG BIN: run BIN with enforce_inv=E on MAT's inputs in MAT/TAG_eE
cd "$(dirname "$0")"
m=$1; e=$2; t=$3; bin=$4
case $m in
  Ge) p=$HOME/PROJECTS/Anharmonic-lattice-dynamics/FOCEX/example/Ge ;;
  BZO) p=$HOME/PROJECTS/MATERIALS/BZO ;;
  MoTe2) p=$HOME/PROJECTS/MATERIALS/MoTe2/TEST ;;
  GeSe) p=$HOME/PROJECTS/MATERIALS/TEST_GeSe ;;
esac
d=$m/${t}_e$e; mkdir -p $d
for f in structure.params dielectric.params latdyn.params kpbs.params default.params POSCAR1 FORCEDISP1; do cp $p/$f $d/; done
awk -v e=$e 'NR==5{$4=e} {print}' $p/structure.params > $d/structure.params
cd $d && /usr/bin/time -l $bin > focex.out 2> time.err
echo "$d exit $? $(grep 'Relative residual' focex.out | tail -1)"
