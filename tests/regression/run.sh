#!/bin/bash
# Build FOCEX and THERMACOND from a source tree and run the regression cases.
#   tests/regression/run.sh SRC LABEL      (SRC: a checkout of Anharmonic-lattice-dynamics)
# Inputs in $REF/inputs (REF=${REF:-~/PROJECTS/consolidate-ref}); builds in $REF/build/LABEL,
# runs in $REF/runs/LABEL.  Compare two labels with compare.py.
# STRICT=1: compile with -ffp-contract=off, so that moving code between routines (which changes
# where gfortran fuses multiply-adds) does not change the last bits; for checking refactoring.
set -e
SRC=$(cd "$1" && pwd); LABEL=$2; REF=${REF:-$HOME/PROJECTS/consolidate-ref}
B=$REF/build/$LABEL; RUN=$REF/runs/$LABEL
export SDKROOT=$(xcrun --show-sdk-path) OMPI_FC=/opt/homebrew/bin/gfortran
FX_FLAGS="-O2 -fcheck=all"; TC_OPT="-O2"
if [ "${STRICT:-0}" = 1 ]; then FX_FLAGS="$FX_FLAGS -ffp-contract=off"; TC_OPT="$TC_OPT -ffp-contract=off"; fi
FC=/opt/homebrew/bin/gfortran
rm -rf "$B" "$RUN"; mkdir -p "$B/focex" "$B/thermacond" "$RUN"
cp "$SRC"/FOCEX/*.f90 "$SRC"/FOCEX/*.f "$SRC"/FOCEX/Makefile "$B/focex/" 2>/dev/null || true
[ -d "$SRC/common" ] && cp "$SRC"/common/*.f90 "$B/focex/" 2>/dev/null || true
(cd "$B/focex" && make FF=$FC "FLAGS=$FX_FLAGS" DESTDIR="$B/bin" > build.log 2>&1) || { echo "FOCEX build failed: $B/focex/build.log"; exit 1; }
cp -r "$SRC"/THERMACOND/. "$B/thermacond/"; rm -rf "$B/thermacond/obj_"* "$B/thermacond/"*example
[ -d "$SRC/common" ] && mkdir -p "$B/common" && cp "$SRC"/common/*.f90 "$B/common/"
(cd "$B/thermacond" && make serial FC=$FC "OPT=$TC_OPT" DESTDIR="$B/bin" > build_serial.log 2>&1 && make "OPT=$TC_OPT" DESTDIR="$B/bin" > build_mpi.log 2>&1) || { echo "THERMACOND build failed: $B/thermacond/build_*.log"; exit 1; }
FX=$(command ls "$B"/bin/focex-*.x | head -1)
set +e   # a case that stops (e.g. BZO in THERMACOND, imaginary modes) must not end the run
pids=()
for c in "$REF"/inputs/focex_*; do n=$(basename "$c"); mkdir -p "$RUN/$n"; cp "$c"/* "$RUN/$n/"
  (cd "$RUN/$n" && "$FX" > focex.out 2>&1; echo "$n rc=$?") & pids+=($!); done
for c in "$REF"/inputs/thermacond_*; do n=$(basename "$c"); mkdir -p "$RUN/$n" "$RUN/${n}_mpi3"; cp "$c"/* "$RUN/$n/"; cp "$c"/* "$RUN/${n}_mpi3/"
  (cd "$RUN/$n" && "$B/bin/tcond-2.0-serial.x" > out 2>&1; echo "$n rc=$?") & pids+=($!); done
wait "${pids[@]}"
for c in "$REF"/inputs/thermacond_*; do n=$(basename "$c"); (cd "$RUN/${n}_mpi3" && mpirun -np 3 "$B/bin/tcond-2.0.x" > out 2>&1; echo "${n}_mpi3 rc=$?"); done
echo "runs in $RUN"
