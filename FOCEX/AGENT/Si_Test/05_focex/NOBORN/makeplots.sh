#!/bin/bash
cd "$(dirname "$0")"
st=0
for f in bvel.plt all.plt qha.plt; do
  gnuplot "$f" > "${f%.plt}.log" 2>&1 && echo "  $f: ok" || { echo "  $f: ERROR, see ${f%.plt}.log"; st=1; }
done
exit $st
