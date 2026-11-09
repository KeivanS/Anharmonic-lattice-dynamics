#!/usr/bin/env bash
# focex-configure.sh -- find what the FOCEX agent needs on this computer and write the settings.
#
#   bash focex-configure.sh              look, report, write ~/.focex-agent.env (nothing else is changed)
#   bash focex-configure.sh --install    the same, and add one line to your shell start-up file that
#                                        reads ~/.focex-agent.env in every new terminal and batch job
# Options (when the search finds the wrong thing, or nothing):
#   --potcar DIR     the POTCAR library (PBE): the folder with Si/, O/, Zn/, ... each holding a POTCAR
#   --vasp FILE      the vasp_std executable
#   --python FILE    the Python interpreter of the environment with numpy, spglib, pymatgen
#   --env FILE       where to write the settings (default ~/.focex-agent.env)
#
# What it looks at: the operating system and processor, the cores, SLURM, gfortran and mpif90, the
# FOCEX / SC-SNAPS / read_outcar.x binaries, Python and its packages, VASP-Flow, the POTCAR library,
# vasp_std, MPI and gnuplot. At the end it checks the settings by building a POTCAR with VASP-Flow's
# own functions and by starting focex-agent.py --version. Written for bash 3.2 (macOS) and Linux.

set -u
ENVFILE="$HOME/.focex-agent.env"
INSTALL=0; OPT_POTCAR=""; OPT_VASP=""; OPT_PY=""
while [ $# -gt 0 ]; do
    case "$1" in
        --install) INSTALL=1 ;;
        --potcar)  OPT_POTCAR="$2"; shift ;;
        --vasp)    OPT_VASP="$2"; shift ;;
        --python)  OPT_PY="$2"; shift ;;
        --env)     ENVFILE="$2"; shift ;;
        -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
        *) echo "unknown option $1 (see --help)"; exit 2 ;;
    esac
    shift
done

AGENT_DIR=$(cd "$(dirname "$0")" && pwd -P)          # FOCEX/AGENT
REPO=$(cd "$AGENT_DIR/../.." && pwd -P)               # Anharmonic-lattice-dynamics
NFAIL=0; NWARN=0; FIXES=""
ok()   { printf '  \033[32mok\033[0m    %s\n' "$1"; }
warn() { printf '  \033[33mWARN\033[0m  %s\n' "$1"; NWARN=$((NWARN+1)); [ -n "${2:-}" ] && FIXES="$FIXES\n  $2"; }
fail() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; NFAIL=$((NFAIL+1)); [ -n "${2:-}" ] && FIXES="$FIXES\n  $2"; }
head_() { printf '\n\033[34m%s\033[0m\n' "$1"; }
expand() { case "$1" in "~"*) echo "$HOME${1#\~}" ;; *) echo "$1" ;; esac; }

# ---------------------------------------------------------------- the computer
head_ "1. This computer"
OS=$(uname -s); ARCH=$(uname -m); HOST=$(hostname -s 2>/dev/null || hostname)
if [ "$OS" = Darwin ]; then
    NCORE=$(sysctl -n hw.physicalcpu 2>/dev/null || echo 1)
    NPERF=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo "$NCORE")   # Apple silicon: performance cores
    ok "macOS $(sw_vers -productVersion 2>/dev/null), $ARCH, $NCORE cores ($NPERF performance cores)"
else
    NCORE=$(lscpu -p=Core,Socket 2>/dev/null | grep -v '^#' | sort -u | wc -l | tr -d ' ')
    [ -z "$NCORE" ] || [ "$NCORE" = 0 ] && NCORE=$(nproc 2>/dev/null || echo 1)
    NPERF=$NCORE
    ok "$OS $ARCH, $NCORE physical cores ($HOST)"
fi
CLUSTER=0
if command -v sbatch >/dev/null 2>&1; then
    CLUSTER=1; ok "SLURM found (sbatch): this is a cluster -- VASP runs need a VASP-Flow profile (see the manual)"
else
    ok "no SLURM: VASP runs on this computer (run_workflow.sh is started with bash)"
fi
SHELLNAME=$(basename "${SHELL:-bash}")

# ---------------------------------------------------------------- compilers (only to build)
head_ "2. Compilers (needed only to build FOCEX and SC-SNAPS)"
if command -v gfortran >/dev/null 2>&1; then ok "gfortran: $(command -v gfortran) ($(gfortran -dumpversion 2>/dev/null))"
else warn "gfortran not found" "install gfortran: macOS 'brew install gcc', Linux 'sudo apt install gfortran', cluster 'module load gcc'"; fi
if command -v mpif90 >/dev/null 2>&1; then ok "mpif90: $(command -v mpif90) (THERMACOND with MPI)"
else warn "mpif90 not found (only THERMACOND's MPI build needs it)" "install MPI: macOS 'brew install open-mpi', Linux 'sudo apt install libopenmpi-dev', cluster 'module load openmpi'"; fi

# ---------------------------------------------------------------- ALATDYN binaries
head_ "3. FOCEX, SC-SNAPS, read_outcar.x"
BIN=$(sed -n 's/^[[:space:]]*DESTDIR[[:space:]]*=[[:space:]]*\([^#]*\).*/\1/p' "$REPO/make.local" "$REPO/make.inc" 2>/dev/null | head -1 | sed 's/[[:space:]]*$//')
BIN=$(expand "${BIN:-~/ALADYN/bin}")
FOCEX_X=$(ls -t "$BIN"/focex-[0-9]*.x 2>/dev/null | head -1)
if [ -n "$FOCEX_X" ]; then ok "FOCEX: $FOCEX_X (the newest; the agent takes the same)"
else fail "no focex-*.x in $BIN" "cd $REPO && make focex sc_snaps"; fi
for x in sc_snaps.x read_outcar.x; do
    if [ -x "$BIN/$x" ]; then ok "$x: $BIN/$x"; else fail "$x not in $BIN" "cd $REPO && make focex sc_snaps"; fi
done
if ls "$BIN"/tcond-[0-9]*.x >/dev/null 2>&1; then ok "THERMACOND: $(ls -t "$BIN"/tcond-[0-9]*.x | grep -v serial | head -1)"
else warn "no THERMACOND in $BIN (needed only for THERMACOND_JOB 1-3)" "cd $REPO && make thermacond"; fi
case ":$PATH:" in *":$BIN:"*) INPATH=1 ;; *) INPATH=0 ;; esac

# ---------------------------------------------------------------- Python
head_ "4. Python with numpy, spglib, pymatgen"
import_hint() {   # python module error -> what to do
    case "$3" in
        *GLIBC*|*.so*|*"cannot open shared"*)
            echo "a library is missing: load the compiler/Python modules the environment was built with (e.g. module load gcc), then run this again" ;;
        *"No module named"*)
            echo "$1 -m pip install numpy spglib pymatgen" ;;
        *)
            echo "see the whole message: $1 -c 'import $2'" ;;
    esac
}
PY=""
for p in "$OPT_PY" "${VIRTUAL_ENV:+$VIRTUAL_ENV/bin/python}" "$HOME/myenv/bin/python" "$(command -v python3 2>/dev/null)"; do
    [ -n "$p" ] && [ -x "$p" ] || continue
    if ERR=$("$p" -c 'import numpy, spglib, pymatgen' 2>&1); then PY=$p; break; fi
    # say why (2026-11-09: on a cluster the packages were in site-packages but the import failed, and
    # "lacks" hid the reason): the Python does not start, or which import fails and its last error line
    if ! V=$("$p" -c 'import sys; print(sys.version.split()[0], sys.prefix)' 2>&1); then
        # a venv records the Python it was made from in pyvenv.cfg (home = its bin folder)
        CFG="$(dirname "$(dirname "$p")")/pyvenv.cfg"
        HOMEPY=$(sed -n 's/^home *= *//p' "$CFG" 2>/dev/null | head -1)
        warn "$p does not start: $(printf '%s' "$V" | tail -1)" \
             "$([ -n "$HOMEPY" ] && echo "this environment was made from the Python in $HOMEPY: load its module (module spider python), add that module load to ~/.bashrc and to the profile's modules, then run this again" || echo "load the module this Python comes from (module load ...), or rebuild ~/myenv with a Python that starts")"
        continue
    fi
    for m in numpy spglib pymatgen; do
        if ! E=$("$p" -c "import $m" 2>&1); then
            warn "$p (Python $V): import $m fails: $(printf '%s' "$E" | tail -1)" "$(import_hint "$p" "$m" "$E")"
            break
        fi
    done
done
if [ -n "$PY" ]; then ok "$PY ($("$PY" --version 2>&1)): numpy, spglib, pymatgen import"
else fail "no Python with numpy, spglib and pymatgen" "python3 -m venv ~/myenv && ~/myenv/bin/pip install numpy spglib pymatgen matplotlib"; fi

# ---------------------------------------------------------------- VASP-Flow
head_ "5. VASP-Flow (the agent uses it to write the VASP inputs)"
VF=$(expand "${VASP_FLOW_DIR:-~/VASP-Flow}")
if [ -f "$VF/vasp-agent.py" ] && [ -f "$VF/modules/vasp_input_generator.py" ]; then ok "VASP-Flow: $VF"
else fail "VASP-Flow not found at $VF" "cd ~ && git clone https://github.com/KeivanS/VASP-Flow.git"; VF=""; fi

# ---------------------------------------------------------------- POTCARs
head_ "6. POTCAR library (from_scratch only)"
has_potcars() { [ -d "$1" ] && ls "$1"/*/POTCAR >/dev/null 2>&1; }
POTCAR=""
if [ -n "$OPT_POTCAR" ]; then
    POTCAR=$(expand "$OPT_POTCAR")
    has_potcars "$POTCAR" || { fail "--potcar $POTCAR has no <element>/POTCAR folders"; POTCAR=""; }
elif [ -n "${VASP_POTCAR_DIR:-}" ] && has_potcars "$(expand "$VASP_POTCAR_DIR")"; then
    POTCAR=$(expand "$VASP_POTCAR_DIR"); ok "already set: VASP_POTCAR_DIR=$POTCAR"
else
    # look for folders named like the PBE libraries VASP ships (potpaw_PBE.54, PAW_PBE, potpaw_PBE.64, ...)
    CANDS=""
    for root in "$HOME" /opt /usr/local /share /apps /sw; do
        [ -d "$root" ] || continue
        f=$(find "$root" -maxdepth 4 -type d \( -iname 'potpaw_PBE*' -o -iname 'PAW_PBE*' \) \
               -not -path '*/Library/*' -not -path '*/.Trash/*' 2>/dev/null)
        [ -n "$f" ] && CANDS="$CANDS$f
"
    done
    for c in $(printf '%s' "$CANDS" | sort -r); do       # sort -r: potpaw_PBE.64 before .54
        if has_potcars "$c"; then
            if [ -z "$POTCAR" ]; then POTCAR=$c; ok "found: $c"; else ok "also found: $c (use --potcar to choose it)"; fi
        fi
    done
fi
if [ -n "$POTCAR" ]; then
    LDA=$(echo "$POTCAR" | sed 's/PBE/LDA/')
    if [ "$LDA" != "$POTCAR" ] && has_potcars "$LDA"; then ok "LDA library next to it: $LDA (used for FUNCTIONAL: LDA)"
    else warn "no LDA library next to it (needed only for FUNCTIONAL: LDA)"; fi
else
    fail "no POTCAR library found" "bash $0 --potcar /path/to/potpaw_PBE.54   (the folder with Si/, O/, ...)"
fi

# ---------------------------------------------------------------- VASP and MPI
head_ "7. VASP and MPI (from_scratch only)"
VASP=""
for v in "$OPT_VASP" "${VASP_STD:-}" "$HOME/BIN/vasp_std" "$(command -v vasp_std 2>/dev/null)"; do
    [ -n "$v" ] || continue; v=$(expand "$v")
    [ -x "$v" ] || continue
    # a binary whose libraries are gone (an old build after a system update) cannot run
    if command -v ldd >/dev/null 2>&1 && ldd "$v" 2>/dev/null | grep -q 'not found'; then
        warn "$v cannot run here: $(ldd "$v" 2>/dev/null | grep -c 'not found') of its libraries are missing (ldd $v)" \
             "remove that vasp_std from VASP_STD / PATH; on a cluster the VASP module of the profile is used"
        continue
    fi
    VASP=$v; break
done
VVER=""
if [ -n "$VASP" ]; then
    VVER=$(grep -a -o -E 'vasp\.[0-9]+\.[0-9]+\.[0-9]+' "$VASP" 2>/dev/null | sed 's/^vasp\.//' | sort -u -t. -k1,1n -k2,2n -k3,3n | tail -1)
    ok "vasp_std: $VASP${VVER:+ (VASP $VVER)}"
    if [ -n "$VVER" ] && [ "${VVER%%.*}" -lt 6 ] 2>/dev/null; then
        warn "VASP $VVER has no r2SCAN (the agent's default FUNCTIONAL; VASP 5 runs PBE instead without a word)" \
             "put FUNCTIONAL: PBE (or PBEsol) in the instruction files, or use VASP 6 (--vasp)"
    fi
elif [ $CLUSTER = 1 ]; then ok "no vasp_std set here: on a cluster the VASP module of the SLURM profile is used (./configure finds it)"
else fail "vasp_std not found" "bash $0 --vasp /path/to/vasp_std"; fi
MPI_LAUNCH=""
if [ $CLUSTER = 1 ]; then
    ok "MPI on the cluster: the profile's mpi_cmd (srun) is used"
elif command -v mpirun >/dev/null 2>&1; then
    MV=$(mpirun --version 2>&1 | head -1)
    case "$MV" in
        *"Open MPI"*|*OpenRTE*) MPI_LAUNCH="mpirun -np" ;;
        *) MPI_LAUNCH="mpiexec -n" ;;                    # MPICH, Intel MPI
    esac
    ok "mpirun: $(command -v mpirun) ($MV) -> MPI_LAUNCH=\"$MPI_LAUNCH\""
else
    warn "no mpirun: VASP runs on one core" "macOS 'brew install open-mpi'; Linux 'sudo apt install openmpi-bin'"
fi
if command -v gnuplot >/dev/null 2>&1 && gnuplot -e 'set terminal pdfcairo' 2>&1 | grep -q 'unknown'; then
    warn "$(command -v gnuplot) has no pdfcairo terminal: the agent writes the plot scripts, run makeplots.sh where gnuplot makes PDFs"
elif command -v gnuplot >/dev/null 2>&1; then ok "gnuplot: $(command -v gnuplot) (the PDF plots, pdfcairo)"
else warn "no gnuplot: the agent writes the plot scripts, run makeplots.sh later" "macOS 'brew install gnuplot'; Linux 'sudo apt install gnuplot'"; fi

# ---------------------------------------------------------------- write the settings
head_ "8. Settings written to $ENVFILE"
[ -f "$ENVFILE" ] && cp "$ENVFILE" "$ENVFILE.bak" && echo "  (the old file is kept as $ENVFILE.bak)"
{
    echo "# FOCEX agent settings, written by focex-configure.sh on $(date '+%Y-%m-%d %H:%M') ($HOST, $OS $ARCH)."
    echo "# Read it with:  source $ENVFILE   (focex-configure.sh --install does that in every new shell)."
    if [ -n "$POTCAR" ]; then echo "export VASP_POTCAR_DIR=\"$POTCAR\""
    else echo "# export VASP_POTCAR_DIR=/path/to/potpaw_PBE.54     # NOT FOUND: set it"; fi
    if [ -n "$VASP" ]; then echo "export VASP_STD=\"$VASP\""
    else echo "# export VASP_STD=/path/to/vasp_std                 # not found here"; fi
    [ -n "$VVER" ] && echo "export VASP_VERSION=\"$VVER\"   # the agent warns when a FUNCTIONAL needs a newer VASP"
    [ -n "$MPI_LAUNCH" ] && echo "export MPI_LAUNCH=\"$MPI_LAUNCH\""
    [ -n "$VF" ] && echo "export VASP_FLOW_DIR=\"$VF\""
    echo "case \":\$PATH:\" in *\":$BIN:\"*) ;; *) export PATH=\"$BIN:\$PATH\" ;; esac   # focex-agent, focex-diagnose"
    if [ -n "$PY" ]; then
        echo "# Python for the agent: $PY"
        # a virtual environment: its bin first on PATH, so python3 (and focex-agent, which starts with
        # '#!/usr/bin/env python3') is this one; moved to the front even if already on PATH (a module
        # loaded later, e.g. python, would otherwise come first)
        PYB=$(dirname "$PY")
        if [ -f "$(dirname "$PYB")/pyvenv.cfg" ]; then
            echo "PATH=\"$PYB:\$(printf '%s' \":\$PATH:\" | sed 's#:$PYB:#:#g; s#^:##; s#:\$##')\"; export PATH   # Python environment: python3 = $PY"
        fi
    fi
    echo "# Suggested in the instruction file: MPI: $NPERF   (without an MPI line VASP runs on 1 core)"
} > "$ENVFILE"
sed 's/^/  | /' "$ENVFILE"

# ---------------------------------------------------------------- start-up file
if [ $INSTALL = 1 ]; then
    head_ "9. Start-up file"
    case "$SHELLNAME" in
        zsh) RC="$HOME/.zshrc" ;;
        *)   if [ "$OS" = Darwin ]; then RC="$HOME/.bash_profile"; else RC="$HOME/.bashrc"; fi ;;
    esac
    LINE="[ -f \"$ENVFILE\" ] && source \"$ENVFILE\"   # FOCEX agent (focex-configure.sh)"
    if [ -f "$RC" ] && grep -qF "$ENVFILE" "$RC"; then ok "$RC reads $ENVFILE already"
    else
        [ -f "$RC" ] && cp "$RC" "$RC.bak-focex"
        echo "$LINE" >> "$RC"; ok "added to $RC (old copy: $RC.bak-focex): $LINE"
    fi
fi

# ---------------------------------------------------------------- check the settings in practice
head_ "10. Check"
# shellcheck disable=SC1090
. "$ENVFILE"
if [ -n "$PY" ] && [ -n "$VF" ] && [ -n "$POTCAR" ]; then
    T=$(mktemp -d)
    if (cd "$T" && "$PY" - "$VF" <<'EOF'
import os, sys, importlib.util
vf = sys.argv[1]; sys.path.insert(0, os.path.join(vf, 'modules'))
from vasp_input_generator import potcar_dir_for
d, how = potcar_dir_for('PBE')
s = importlib.util.spec_from_file_location('va', os.path.join(vf, 'vasp-agent.py'))
m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
el = sorted(e for e in os.listdir(d) if os.path.isfile(os.path.join(d, e, 'POTCAR')))[0]
m.build_potcar([el], d, 'POTCAR')
print('  POTCAR test (%s): %s' % (el, open('POTCAR').readline().strip()))
EOF
    ) >"$T/out" 2>&1 && grep -q 'PAW' "$T/out"; then ok "VASP-Flow built a POTCAR from VASP_POTCAR_DIR: $(grep 'POTCAR test' "$T/out" | sed 's/^ *//')"
    else fail "VASP-Flow could not build a POTCAR:"; sed 's/^/        /' "$T/out"; fi
    rm -rf "$T"
fi
if [ -n "$PY" ] && [ -n "$VF" ]; then
    if V=$("$PY" "$AGENT_DIR/focex-agent.py" --version 2>&1); then ok "focex-agent.py starts: $V"
    else fail "focex-agent.py does not start:"; echo "$V" | sed 's/^/        /'; fi
fi

echo
if [ $NFAIL = 0 ]; then
    printf '\033[32mReady\033[0m (%d warning(s)).' "$NWARN"
    if [ $INSTALL = 1 ]; then echo " Open a new terminal (or: source $ENVFILE), then follow the Quick start."
    elif [ -n "${ALADYN_CONFIGURE:-}" ]; then echo     # run by configure, which reads this file from ~/.aladyn.env
    else echo " Next: bash $0 --install   (or add 'source $ENVFILE' to your start-up file yourself)"; fi
else
    printf '\033[31m%d problem(s)\033[0m, %d warning(s). To fix:\n' "$NFAIL" "$NWARN"; printf '%b\n' "$FIXES" | awk 'NF && !seen[$0]++'
    echo "Then run this script again."
fi
[ $NFAIL = 0 ]
