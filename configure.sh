#!/usr/bin/env bash
# configure.sh -- set up the ALADYN suite (FOCEX, SC-SNAPS, THERMACOND, ANFOMOD, the agent and the GUI) on
# this computer: find the hardware and software, write make.local for the Makefiles (read at the end of
# make.inc; not in git, so git pull never conflicts with it) and ~/.aladyn.env for the shell (2026-11-09).
#
#   bash configure.sh              look, test the compilers, write make.local and ~/.aladyn.env (no start-up file
#                                  is changed); also the agent's checks (FOCEX/AGENT/focex-configure.sh)
#   bash configure.sh --install    the same, and add one line to your shell start-up file (~/.bash_profile on
#                                  macOS with bash, ~/.bashrc on Linux, ~/.zshrc with zsh) that reads ~/.aladyn.env
#   bash configure.sh --build      the same, then make (everything; see the Makefile)
# Options:
#   --prefix DIR     where the programs are installed (DESTDIR; default: that of make.local, else make.inc)
#   --fc FILE        the Fortran compiler (default: gfortran, else ifx, else ifort)
#   --mpifc FILE     the MPI Fortran wrapper (default: mpif90)
#   --no-agent       skip the FOCEX agent checks (VASP, POTCARs, Python, VASP-Flow)
#   --makeinc FILE   write make.local there instead (for a test); --env FILE, --agent-env FILE the same for
#                    ~/.aladyn.env and ~/.focex-agent.env
#   (the agent options --potcar, --vasp, --python are passed on to focex-configure.sh)
# Written for bash 3.2 (macOS) and Linux.

set -u
REPO=$(cd "$(dirname "$0")" && pwd -P)
MAKEINC="$REPO/make.local"; ENVFILE="$HOME/.aladyn.env"
AGENTENV="$HOME/.focex-agent.env"; INSTALL=0; BUILD=0; AGENT=1; PREFIX=""; OPT_FC=""; OPT_MPIFC=""; AGENT_OPTS=""
while [ $# -gt 0 ]; do
    case "$1" in
        --install) INSTALL=1 ;;
        --build)   BUILD=1 ;;
        --no-agent) AGENT=0 ;;
        --prefix)  PREFIX="$2"; shift ;;
        --fc)      OPT_FC="$2"; shift ;;
        --mpifc)   OPT_MPIFC="$2"; shift ;;
        --makeinc) MAKEINC="$2"; shift ;;
        --env)     ENVFILE="$2"; shift ;;
        --agent-env) AGENTENV="$2"; shift ;;
        --potcar|--vasp|--python) AGENT_OPTS="$AGENT_OPTS $1 $2"; shift ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "unknown option $1 (see --help)"; exit 2 ;;
    esac
    shift
done
NFAIL=0; NWARN=0; FIXES=""
ok()   { printf '  \033[32mok\033[0m    %s\n' "$1"; }
warn() { printf '  \033[33mWARN\033[0m  %s\n' "$1"; NWARN=$((NWARN+1)); [ -n "${2:-}" ] && FIXES="$FIXES\n  $2"; }
fail() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; NFAIL=$((NFAIL+1)); [ -n "${2:-}" ] && FIXES="$FIXES\n  $2"; }
head_() { printf '\n\033[34m%s\033[0m\n' "$1"; }
expand() { case "$1" in "~"*) echo "$HOME${1#\~}" ;; *) echo "$1" ;; esac; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

# ---------------------------------------------------------------- 1. the computer
head_ "1. This computer"
OS=$(uname -s); ARCH=$(uname -m); HOST=$(hostname -s 2>/dev/null || hostname)
if [ "$OS" = Darwin ]; then
    NCORE=$(sysctl -n hw.physicalcpu 2>/dev/null || echo 1)
    NPERF=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo "$NCORE")
    MEM=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 ))
    ok "macOS $(sw_vers -productVersion 2>/dev/null), $ARCH, $NCORE cores ($NPERF performance), ${MEM} GB"
else
    NCORE=$(lscpu -p=Core,Socket 2>/dev/null | grep -v '^#' | sort -u | wc -l | tr -d ' ')
    { [ -z "$NCORE" ] || [ "$NCORE" = 0 ]; } && NCORE=$(nproc 2>/dev/null || echo 1)
    NPERF=$NCORE
    MEM=$(awk '/MemTotal/ {printf "%d", $2/1048576}' /proc/meminfo 2>/dev/null)
    ok "$OS $ARCH, $NCORE physical cores, ${MEM:-?} GB ($HOST)"
fi
CLUSTER=0
if command -v sbatch >/dev/null 2>&1; then CLUSTER=1; ok "SLURM (sbatch): a cluster -- build on a login node, run in jobs"
else ok "no SLURM"; fi
if command -v module >/dev/null 2>&1 || [ -n "${MODULESHOME:-}" ]; then
    ok "environment modules: loaded now: ${LOADEDMODULES:-none} (load the compiler and MPI modules before configure.sh)"
fi
SHELLNAME=$(basename "${SHELL:-bash}")

# ---------------------------------------------------------------- 2. compilers
head_ "2. Fortran and MPI compilers"
FFC=""
for c in "$OPT_FC" gfortran ifx ifort; do
    [ -n "$c" ] || continue
    p=$(command -v "$c" 2>/dev/null) || continue
    FFC=$p; break
done
if [ -n "$FFC" ]; then
    case "$(basename "$FFC")" in
        gfortran*) FVER=$("$FFC" -dumpversion 2>/dev/null); FAM=gnu ;;
        *)         FVER=$("$FFC" --version 2>/dev/null | head -1); FAM=intel ;;
    esac
    ok "Fortran: $FFC ($FVER)"
    if [ "$FAM" = gnu ] && [ "${FVER%%.*}" -lt 9 ] 2>/dev/null; then
        warn "gfortran $FVER: version 9 or newer is needed" "install a newer gfortran (macOS: brew install gcc)"
    fi
else
    fail "no Fortran compiler (gfortran, ifx, ifort)" "macOS 'brew install gcc'; Linux 'sudo apt install gfortran'; cluster 'module load gcc'"
    FAM=gnu
fi
MPIFC=""
p=$(command -v "${OPT_MPIFC:-mpif90}" 2>/dev/null) && MPIFC=$p
OMPI_FC_SET=""; SDK=""
if [ -n "$MPIFC" ]; then
    UNDER=$("$MPIFC" --showme:command 2>/dev/null | head -1)              # Open MPI
    [ -z "$UNDER" ] && UNDER=$("$MPIFC" -show 2>/dev/null | awk '{print $1}')   # MPICH, Intel MPI
    ok "MPI Fortran: $MPIFC (calls ${UNDER:-?})"
    if [ -n "$FFC" ] && [ -n "$UNDER" ]; then
        U=$(command -v "$UNDER" 2>/dev/null || echo "$UNDER")
        [ "$U" != "$FFC" ] && ok "note: FOCEX and SC-SNAPS use $FFC, MPI's wrapper calls $U (two compilers; fine, THERMACOND is built with the wrapper)"
    fi
    MV=$(mpirun --version 2>&1 | head -1)
    # Homebrew Open MPI on macOS: mpif90 must call the gfortran its mpi.mod was built with, and the linker
    # needs the SDK (install.rst); both go into make.local and ~/.aladyn.env
    if [ "$OS" = Darwin ] && echo "$MV" | grep -q "Open MPI"; then
        BREWGF=$(dirname "$(command -v mpirun)")/gfortran
        if [ -x "$BREWGF" ]; then OMPI_FC_SET=$BREWGF; ok "Open MPI on macOS: OMPI_FC=$BREWGF (its gfortran)"; fi
        SDK=$(xcrun --show-sdk-path 2>/dev/null); [ -n "$SDK" ] && ok "SDKROOT=$SDK"
    fi
else
    warn "no mpif90: THERMACOND is built without MPI (make thermacond-serial)" \
         "macOS 'brew install open-mpi'; Linux 'sudo apt install libopenmpi-dev'; cluster 'module load openmpi'"
fi

# ---------------------------------------------------------------- 3. compile tests
head_ "3. Compiler tests"
[ -n "$OMPI_FC_SET" ] && export OMPI_FC="$OMPI_FC_SET"
[ -n "$SDK" ] && export SDKROOT="$SDK"
OMPFLAG=-fopenmp; [ "$FAM" = intel ] && OMPFLAG=-qopenmp
if [ -n "$FFC" ]; then
    cat > "$T/omp.f90" <<'EOF'
program t
!$ use omp_lib
integer :: n = 1
!$ n = omp_get_max_threads()
print '(a,i0)', 'threads ', n
end program t
EOF
    if (cd "$T" && "$FFC" $OMPFLAG omp.f90 -o omp.x >cc.log 2>&1 && OMP_NUM_THREADS=2 ./omp.x | grep -q 'threads 2'); then
        ok "$(basename "$FFC") $OMPFLAG: compiles and runs on 2 OpenMP threads"
    else fail "$(basename "$FFC") $OMPFLAG does not work:"; sed 's/^/        /' "$T/cc.log" | head -5; fi
fi
if [ -n "$MPIFC" ]; then
    cat > "$T/mpi.f90" <<'EOF'
program t
use mpi
integer ierr, r, n
call MPI_INIT(ierr); call MPI_COMM_RANK(MPI_COMM_WORLD, r, ierr); call MPI_COMM_SIZE(MPI_COMM_WORLD, n, ierr)
if (r == 0) print '(a,i0)', 'ranks ', n
call MPI_FINALIZE(ierr)
end program t
EOF
    if (cd "$T" && "$MPIFC" mpi.f90 -o mpi.x >mc.log 2>&1); then
        if [ $CLUSTER = 1 ]; then ok "$(basename "$MPIFC"): compiles an MPI program (not run on a login node)"
        elif (cd "$T" && mpirun -np 2 ./mpi.x 2>&1 | grep -q 'ranks 2'); then ok "$(basename "$MPIFC"): compiles; mpirun -np 2 runs"
        else warn "MPI program compiled but mpirun -np 2 did not run it" "check mpirun (e.g. mpirun --oversubscribe)"; fi
    else
        fail "$(basename "$MPIFC") cannot compile 'use mpi':"; sed 's/^/        /' "$T/mc.log" | head -5
        FIXES="$FIXES\n  mpif90 and its mpi.mod come from different compilers: set OMPI_FC to the gfortran Open MPI was built with"
    fi
fi

# ---------------------------------------------------------------- 4. make.local
head_ "4. make.local (read at the end of make.inc by every Makefile) -> $MAKEINC"
OLD_DEST=$(sed -n 's/^[[:space:]]*DESTDIR[[:space:]]*=[[:space:]]*\([^#]*\).*/\1/p' "$REPO/make.local" "$REPO/make.inc" 2>/dev/null | head -1 | sed 's/[[:space:]]*$//')
DEST=${PREFIX:-${OLD_DEST:-~/ALADYN/bin}}
if [ "$FAM" = intel ]; then FLAGS="-O2 -check bounds"; else FLAGS="-O2 -fcheck=all -ffree-line-length-0"; fi
[ -f "$MAKEINC" ] && cp "$MAKEINC" "$MAKEINC.bak" && echo "  (the old file is kept as $MAKEINC.bak)"
{
    # GNU make keeps the blanks before a '#' in a value: every comment is on its own line
    echo "# make.local -- written by configure.sh on $(date '+%Y-%m-%d %H:%M') ($HOST, $OS $ARCH); run it again to redo"
    echo "# where make installs the programs (focex-*.x, tcond-*.x, sc_snaps.x, the links of make python-tools)"
    echo "DESTDIR = $DEST"
    echo "# Fortran compiler: FOCEX, FOCEX/utility, ANFOMOD, SC-SNAPS"
    echo "FF = ${FFC:-gfortran}"
    echo "# MPI Fortran wrapper: THERMACOND (its serial build uses the compiler mpif90 calls)"
    echo "MPIFC = ${MPIFC:-mpif90}"
    echo "# SCOP8 (its older Makefile)"
    echo "FC = \$(MPIFC)"
    echo "F90 = \$(MPIFC)"
    echo "# options of FOCEX/utility and ANFOMOD (FOCEX and THERMACOND set their own)"
    echo "FLAGS = $FLAGS"
    if [ -n "$OMPI_FC_SET" ]; then echo "# Open MPI on macOS: the gfortran its mpi.mod was built with"; echo "OMPI_FC = $OMPI_FC_SET"
    else echo "# OMPI_FC: only for an Open MPI whose mpi.mod comes from another gfortran"; fi
    if [ -n "$SDK" ]; then echo "# macOS: the SDK the linker needs"; echo "SDKROOT = $SDK"; fi
} > "$MAKEINC"
sed 's/^/  | /' "$MAKEINC"

# ---------------------------------------------------------------- 5. the FOCEX agent and the GUI
if [ $AGENT = 1 ]; then
    head_ "5. FOCEX agent and GUI (FOCEX/AGENT/focex-configure.sh)"
    # shellcheck disable=SC2086
    ALADYN_CONFIGURE=1 bash "$REPO/FOCEX/AGENT/focex-configure.sh" --env "$AGENTENV" $AGENT_OPTS 2>&1 | sed 's/^/    /'
    AG=${PIPESTATUS[0]}
    [ "$AG" = 0 ] || warn "the agent checks report problems (above); FOCEX itself does not need them"
    PYA=$(sed -n 's/^# Python for the agent: //p' "$AGENTENV" 2>/dev/null)
    if [ -n "$PYA" ]; then
        if "$PYA" -c 'import flask' >/dev/null 2>&1; then ok "Flask in $PYA (aladyn-gui)"
        else warn "no Flask in $PYA: aladyn-gui needs it" "$PYA -m pip install flask"; fi
    fi
fi

# ---------------------------------------------------------------- 6. ~/.aladyn.env
head_ "6. Shell settings -> $ENVFILE"
BIN=$(expand "$DEST")
[ -f "$ENVFILE" ] && cp "$ENVFILE" "$ENVFILE.bak" && echo "  (the old file is kept as $ENVFILE.bak)"
{
    echo "# ALADYN settings, written by configure.sh on $(date '+%Y-%m-%d %H:%M') ($HOST). Read with: . $ENVFILE"
    echo "case \":\$PATH:\" in *\":$BIN:\"*) ;; *) export PATH=\"$BIN:\$PATH\" ;; esac   # the ALADYN programs"
    [ -n "$OMPI_FC_SET" ] && echo "export OMPI_FC=\"$OMPI_FC_SET\"   # Open MPI on macOS: mpif90 calls this gfortran"
    [ -n "$SDK" ] && echo "export SDKROOT=\"\$(xcrun --show-sdk-path 2>/dev/null)\"   # macOS linker SDK (follows Xcode updates)"
    echo "export OPENBLAS_NUM_THREADS=1   # one BLAS thread per MPI rank"
    echo "[ -f \"$AGENTENV\" ] && . \"$AGENTENV\"   # FOCEX agent: VASP, POTCARs, VASP-Flow"
    echo "# FOCEX and THERMACOND choose their threads (FOCEX_NUM_THREADS / TCOND_NUM_THREADS override; $NPERF performance cores here)"
} > "$ENVFILE"
sed 's/^/  | /' "$ENVFILE"
if [ $INSTALL = 1 ]; then
    case "$SHELLNAME" in
        zsh) RC="$HOME/.zshrc" ;;
        *)   if [ "$OS" = Darwin ]; then RC="$HOME/.bash_profile"; else RC="$HOME/.bashrc"; fi ;;
    esac
    LINE="[ -f \"$ENVFILE\" ] && . \"$ENVFILE\"   # ALADYN (configure.sh)"
    if [ -f "$RC" ] && grep -qF "$ENVFILE" "$RC"; then ok "$RC reads $ENVFILE already"
    else
        [ -f "$RC" ] && cp "$RC" "$RC.bak-aladyn"
        echo "$LINE" >> "$RC"; ok "added to $RC (old copy: $RC.bak-aladyn): $LINE"
    fi
fi

# ---------------------------------------------------------------- 7. build
if [ $BUILD = 1 ] && [ $NFAIL = 0 ]; then
    head_ "7. Build (make)"
    if (cd "$REPO" && make 2>&1 | tail -5); then ok "make done; programs in $BIN"; else fail "make failed (see above)"; fi
fi

echo
if [ $NFAIL = 0 ]; then
    printf '\033[32mReady\033[0m (%d warning(s)).' "$NWARN"
    [ $BUILD = 1 ] || printf ' Next: make   (in %s)' "$REPO"
    [ $INSTALL = 1 ] && echo " Open a new terminal (or: . $ENVFILE)." || echo "; bash configure.sh --install adds $ENVFILE to your start-up file."
else
    printf '\033[31m%d problem(s)\033[0m, %d warning(s). To fix:\n' "$NFAIL" "$NWARN"; printf '%b\n' "$FIXES" | awk 'NF && !seen[$0]++'
    echo "Then run configure.sh again."
fi
[ $NFAIL = 0 ]
