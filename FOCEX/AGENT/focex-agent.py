#!/usr/bin/env python3
"""focex-agent.py 1.10 (FOCEX/AGENT) -- force constants and phonons with FOCEX, from a crystal
structure (MODE: from_scratch, VASP + SC_SNAPS + FOCEX) or from existing force-displacement data
(MODE: extract_only, FOCEX only).  Manual: FOCEX/AGENT/focexagent_manual.pdf

    focex-agent.py [--instr FILE] [--coord FILE] [--profile NAME] [--force]
        sets up the project folder (PROJECT_DIR in the instruction file) with run_workflow.sh;
        FILE defaults: focex-agent-instructions.txt, POSCAR; NAME: a VASP-Flow profile
        (workstation, slurm, or your own = ~/VASP-Flow/profiles/<NAME>.json)
    cd <PROJECT_DIR> && sbatch run_workflow.sh     (bash run_workflow.sh without SLURM)

    focex-agent.py run  <PROJECT_DIR> [--only TASK] [--retry-failed]   what run_workflow.sh calls
    focex-agent.py plots <PROJECT_DIR>                                  remake the plots
    diagnose.py <PROJECT_DIR>         where and why a run stopped; writes diagnosis/fix_it.sh

Tasks (each is recorded in <PROJECT_DIR>/report.log and skipped once done):
  relax      VASP relaxation of the reduced, oriented cell (ACCURACY: Accurate|Normal)
  poscar0    CONTCAR -> Niggli reduction + orientation -> POSCAR0 (frame of every later step)
  scf        SCF on POSCAR0, same INCAR (ISIF=1, NSW=0) and KPOINTS
  dfpt       epsilon_inf and Born charges in 06_dfpt, only if the SCF shows a gap (DFPT: auto|on|off);
             EFIELD_PEAD at half the Zener limit, checks: Zener warning, sum of charges, eps_inf
  snaps      SC_SNAPS: supercell (SUPERCELL: SMALL|LARGE|9 integers, RS_i = sum_j n_ij R0_j) and
             snapshots at TEMPERATURE; if none is given, 30 K tuned to an rms of DISP_TARGET
  supercells one folder per snapshot (INCAR as the SCF + ISYM=0, commensurate KPOINTS)
  vasp_snaps the supercell VASP runs, one after the other
  collect    read_outcar.x -> FORCEDISP1, POSCAR1 = poscar_000; all runs checked
  focex      the five *.params files, FOCEX with enforce_inv 1 (+ NOBORN, and LASSO with refit, when needed)
  plots      PDF plots by PLOTTING/allplots (bands, Gruneisen, velocities, QHA, ...), makeplots.sh
  thermacond THERMACOND in 07_thermacond if THERMACOND_JOB (latdyn.params line 8) is 1, 2 or 3
VASP_RUN: external (THERMACOND_RUN: external): the inputs of each stage are written and the agent stops
as WAITING (exit code 3); it goes on from the outputs once they are in the folders (aladyn-gui, bundle.py).
MODE: extract_only (DATA: <folder> with POSCAR1, FORCEDISP1, [POSCAR2, FORCEDISP2, ...],
optionally POSCAR = primitive cell and dielectric.params, all in the frame of the data) rotates
the data into the FOCEX frame and does focex and plots only.

Every generated input can be replaced by a block in the instruction file, named after the
file:  INCAR ... INCAR_END  (also INCAR_relax, INCAR_scf, INCAR_dfpt, INCAR_snaps),
KPOINTS, KPOINTS_snaps, cell_inp, supercell_inp, snaps_inp, structure_params,
latdyn_params, kpbs_params, dielectric_params, default_params, thermacond_params (the latdyn.params
of 07_thermacond)  (END_<NAME> also accepted,
any case).  INCAR blocks are merged tag by tag; the others are written verbatim.
VASP keywords (FUNCTIONAL, MAGMOM, GGA_U, MPI, ENCUT, ...) are those of VASP-Flow.
"""
import argparse
import datetime
import glob
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True                  # no __pycache__ in utility/ or VASP-Flow
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
# one thread per MPI rank (VASP-Flow's job scripts do the same): Homebrew's OpenBLAS is an OpenMP
# build that ignores OPENBLAS_NUM_THREADS -- ZnO, 96 atoms, 8 ranks: 17 threads per rank and no
# SCF step in 30 minutes
os.environ.setdefault('OMP_NUM_THREADS', os.environ.get('VASP_OMP_NUM_THREADS', '1'))
import numpy as np

__version__ = '1.10'


def gnuplot_pdf():
    """gnuplot if it has the pdfcairo terminal the plot scripts use, else None (2026-11-09: a cluster's gnuplot
    module had no pdfcairo; every plot failed and the run was reported FAILED after a good FOCEX fit)."""
    gp = shutil.which('gnuplot')
    if gp is None:
        return None
    try:
        r = subprocess.run([gp, '-e', 'set terminal pdfcairo'], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    return None if (r.returncode != 0 or 'unknown' in (r.stdout + r.stderr)) else gp
AGENT = os.path.dirname(os.path.realpath(__file__))          # FOCEX/AGENT (realpath: works through a symlink)
HERE = os.path.dirname(AGENT)                                # the FOCEX folder: templates, utility/
VASPFLOW = os.path.expanduser(os.environ.get('VASP_FLOW_DIR', '~/VASP-Flow'))


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


p2f = _load(os.path.join(HERE, 'utility', 'poscar2focex.py'), 'poscar2focex')
sys.path.insert(0, os.path.join(VASPFLOW, 'modules'))
try:
    from instruction_parser import InstructionParser, InstructionError
    from vasp_input_generator import VASPInputGenerator, potcar_dir_for, slurm_account
    vagent = _load(os.path.join(VASPFLOW, 'vasp-agent.py'), 'vasp_agent')
except Exception as e:                                        # pragma: no cover
    sys.exit(f'focex-agent needs VASP-Flow at {VASPFLOW} (set VASP_FLOW_DIR): {e}')

ACCURACY = {'accurate': dict(ENCUT=600, EDIFF='1E-7', EDIFFG=-0.003, PREC='Accurate'),   # EDIFFG -0.001 until 2026-11-09 (user)
            'normal':   dict(ENCUT=450, EDIFF='1E-6', EDIFFG=-0.003, PREC='Normal')}
SUPERCELL = {'small': (10.0, 100), 'large': (14.0, 300)}     # min face distance (A), max atoms
OMEGA2_CM = 521.11 ** 2          # bands.dat: omega^2 in eV/(A^2 amu) -> (1/cm)^2
METALS = set('''Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Tc Ru Rh
 Pd Ag Cd In Sn Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg
 Tl Pb Bi Po Fr Ra Ac Th Pa U Np Pu Am Cm'''.split())
BLOCKS = ('incar', 'incar_relax', 'incar_scf', 'incar_dfpt', 'incar_snaps', 'kpoints',
          'kpoints_snaps', 'cell_inp', 'supercell_inp', 'snaps_inp', 'structure_params',
          'latdyn_params', 'kpbs_params', 'dielectric_params', 'default_params', 'thermacond_params')
INSTR = 'focex-agent-instructions.txt'    # default instruction file; its copy in the project folder
OLD_INSTR = 'instructions.txt'            # name of that copy before 2026-10-23
MODES = ('from_scratch', 'extract_only')
TASKS = ('relax', 'poscar0', 'scf', 'dfpt', 'snaps', 'supercells', 'vasp_snaps', 'collect',
         'focex', 'plots', 'thermacond')
TC_JOB = {0: 'bands, DOS, Gruneisen only (no THERMACOND run)', 1: 'linewidths on the band path',
          2: 'RTA + direct solution', 3: 'RTA + iterative solution'}


# ── instructions ──────────────────────────────────────────────────────────────
def split_blocks(text):
    """(blocks, rest): {name: text} for the FOCEX-style blocks, and the other lines."""
    blocks, rest, cur, body = {}, [], None, []
    for line in text.splitlines():
        key = line.split('#')[0].strip().rstrip(':').strip().lower()
        if cur is None:
            if key in BLOCKS:
                cur, body = key, []
            else:
                rest.append(line)
        elif key in (cur + '_end', 'end_' + cur, cur.replace('_', ' ') + ' end'):
            blocks[cur] = '\n'.join(body) + '\n'
            cur = None
        else:
            body.append(line)
    if cur is not None:
        sys.exit(f'instructions: block {cur.upper()} has no {cur.upper()}_END line')
    return blocks, '\n'.join(rest) + '\n'


def keyword(text, name, default=None):
    m = re.search(r'^[ \t]*' + name + r'[ \t]*[:=][ \t]*([^#\n]*)', text, re.I | re.M)
    return m.group(1).strip() if m and m.group(1).strip() else default


def incar_lines(block):
    return [l.strip() for l in (block or '').splitlines() if '=' in l and not l.strip().startswith(('#', '!'))]


def has_tag(lines, tag):
    return any(l.split('=', 1)[0].strip().upper() == tag for l in lines)


def parse_instructions(path):
    if not os.path.isfile(path):
        sys.exit(f'instruction file {path} not found (default name {INSTR}; another name: --instr FILE)')
    raw = open(path).read()
    blocks, rest = split_blocks(raw)
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as fh:
        fh.write(rest)
    try:
        inst = InstructionParser(fh.name).instructions
    except InstructionError as e:
        sys.exit(str(e))
    finally:
        os.unlink(fh.name)
    fx = {}
    fx['project'] = inst.get('project_name') or 'FOCEX_project'
    acc = (keyword(rest, 'ACCURACY') or keyword(rest, 'PRECISION') or 'Accurate').lower()
    if acc not in ACCURACY:
        sys.exit(f'ACCURACY must be Accurate or Normal, not {acc!r}')
    fx['accuracy'] = acc
    sc = (keyword(rest, 'SUPERCELL') or 'SMALL').split()
    if sc[0].lower() in SUPERCELL:
        fx['supercell'] = sc[0].lower()
    elif len(sc) == 9:
        fx['supercell'] = [[int(x) for x in sc[3 * i:3 * i + 3]] for i in range(3)]
    else:
        sys.exit('SUPERCELL must be SMALL, LARGE or 9 integers (rows, primitive units)')
    fx['nsnap'] = int(keyword(rest, 'SNAPSHOTS', 50))
    # a TEMPERATURE given is kept (e.g. effective phonons at that T); only the default 30 K is
    # retuned so that the rms displacement is DISP_TARGET
    fx['temperature_given'] = keyword(rest, 'TEMPERATURE') is not None
    fx['temperature'] = float(keyword(rest, 'TEMPERATURE', 30))
    fx['avg_freq'] = float(keyword(rest, 'AVG_FREQ', 400))
    fx['disp_target'] = float(keyword(rest, 'DISP_TARGET', 0.02))
    r = keyword(rest, 'RANKS')
    if not r and 'structure_params' not in blocks:
        sys.exit('RANKS must be given, e.g.  RANKS: 1 1 0 0 0 0 0 0  (harmonic) or '
                 '1 1 1 0 0 0 0 0 (+cubic) or 1 1 1 1 0 0 0 0 (+quartic)')
    if r:
        v = [int(x) for x in r.split()]
        fx['ranks'] = v + [0] * (8 - len(v)) if len(v) > 1 else [1] * v[0] + [0] * (8 - v[0])
    shells = [int(x) for x in (keyword(rest, 'RANGE') or '6 2 1 0 0 0 0').split()]
    fx['shells'] = (shells + [0] * 7)[:7]                       # ranks 2..8
    fx['fc2range'] = int(keyword(rest, 'FC2RANGE', 0))
    fx['focex_mesh'] = keyword(rest, 'FOCEX_MESH', 'fine')
    fx['kpath'] = int(keyword(rest, 'KPATH', 0))
    fx['bs_nk'] = int(keyword(rest, 'BS_NKPTS', 40))
    fx['dfpt'] = (keyword(rest, 'DFPT') or 'auto').lower()
    fx['lasso_ratio'] = float(keyword(rest, 'LASSO_RATIO', 2.0))
    fx['born_transpose'] = (keyword(rest, 'BORN_TRANSPOSE') or 'yes').lower() in ('yes', 'true', '1')
    fx['focex_x'] = os.path.expanduser(keyword(rest, 'FOCEX_X') or newest_focex())
    fx['data'] = os.path.expanduser(keyword(rest, 'DATA') or '')
    mode = (keyword(rest, 'MODE') or ('extract_only' if fx['data'] else 'from_scratch')).lower()
    if mode not in MODES:
        sys.exit(f'MODE must be from_scratch or extract_only, not {mode!r}')
    if mode == 'extract_only' and not fx['data']:
        sys.exit('MODE: extract_only needs DATA: <folder with POSCAR1 and FORCEDISP1>')
    if mode == 'from_scratch' and fx['data']:
        sys.exit('DATA is read only with MODE: extract_only')
    fx['mode'] = mode
    # STOP_AFTER: <task> -- the run stops (OK) after that task, e.g. relax to check the structure first;
    # deleting the line and submitting run_workflow.sh again goes on from there
    # STRESS_OFFSET: xx yy zz -- known bias (kB) of VASP's printed stress (e.g. rVV10 in VASP <= 6.2), measured as
    # printed minus -(1/V) dE/d eps; subtracted before stress.params is written. Default: none.
    so = keyword(rest, 'STRESS_OFFSET')
    try:
        fx['stress_offset'] = [float(x) for x in so.split('#')[0].split()] if so else None
    except ValueError:
        sys.exit(f'STRESS_OFFSET must be three numbers (kB: xx yy zz), not {so!r}')
    if fx['stress_offset'] is not None and len(fx['stress_offset']) != 3:
        sys.exit(f'STRESS_OFFSET must be three numbers (kB: xx yy zz), not {so!r}')
    fx['stop_after'] = (keyword(rest, 'STOP_AFTER') or '').split('#')[0].strip().lower()
    if fx['stop_after'] and fx['stop_after'] not in TASKS:
        sys.exit(f'STOP_AFTER must be one of {", ".join(TASKS)}, not {fx["stop_after"]!r}')
    # where the calculations go; relative to the folder of the instruction file
    pdir = os.path.expanduser(keyword(rest, 'PROJECT_DIR') or vagent.slugify(fx['project']))
    fx['project_dir'] = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(path)), pdir))
    fx['species'] = (keyword(rest, 'SPECIES') or '').split()
    fx['sc_snaps_x'] = os.path.expanduser(keyword(rest, 'SC_SNAPS_X') or first_file(
        '~/ALADYN/bin/sc_snaps.x', '~/BIN/sc_snaps.x'))
    fx['read_outcar_x'] = os.path.expanduser(keyword(rest, 'READ_OUTCAR_X') or '~/ALADYN/bin/read_outcar.x')
    # PROFILE: of the instruction file only; the default FOCEX_PROFILE (written to ~/.aladyn.env by ./configure on
    # a cluster) is applied at setup (init) after --profile, and the choice is kept in state.json: at run time
    # it must not replace the --profile given at setup (it did in 1.8-1.9: a --profile workstation test ran srun)
    fx['profile'] = keyword(rest, 'PROFILE') or None
    fx['mpi_given'] = keyword(rest, 'MPI') is not None      # the parser gives MPI = 1 when it is missing
    # THERMACOND: latdyn.params line 8 decides; THERMACOND_JOB, else line 8 of a latdyn_params block
    job = keyword(rest, 'THERMACOND_JOB')
    for blk in ('thermacond_params', 'latdyn_params'):        # line 8 of the THERMACOND latdyn.params
        if job is None and blocks.get(blk):
            lines = blocks[blk].splitlines()
            job = lines[7].split()[0] if len(lines) > 7 and lines[7].split() else None
    try:
        fx['tc_job'] = int(job) if job is not None else 0
    except ValueError:
        sys.exit(f'THERMACOND_JOB must be 0, 1, 2 or 3, not {job!r}')
    if fx['tc_job'] not in TC_JOB:
        sys.exit(f'THERMACOND_JOB must be 0, 1, 2 or 3, not {fx["tc_job"]}')
    fx['tc_job_given'] = job is not None
    fx['tc_mesh'] = (keyword(rest, 'THERMACOND_MESH') or '').split()
    fx['tc_T'] = (keyword(rest, 'THERMACOND_T') or '').split()
    fx['tc_solver'] = keyword(rest, 'THERMACOND_SOLVER')
    if fx['tc_mesh'] and len(fx['tc_mesh']) != 3 or fx['tc_T'] and len(fx['tc_T']) != 3:
        sys.exit('THERMACOND_MESH needs 3 integers (n1 n2 n3), THERMACOND_T 3 numbers (Tmin Tmax nT)')
    fx['thermacond_x'] = os.path.expanduser(keyword(rest, 'THERMACOND_X') or newest_tcond())
    # where VASP and THERMACOND run: 'here' (the agent starts them) or 'external' (the agent writes the
    # inputs and stops as WAITING; aladyn-gui bundles them for a cluster, or submits them there, and the
    # agent continues when the outputs are in the folders)
    for k, name in (('vasp_run', 'VASP_RUN'), ('tc_run', 'THERMACOND_RUN')):
        fx[k] = (keyword(rest, name) or 'here').lower()
        if fx[k] not in ('here', 'external'):
            sys.exit(f'{name} must be here or external, not {fx[k]!r}')
    fx['blocks'] = blocks
    return inst, fx


def slurm_setup(inst, prof, fx):
    """A SLURM profile: NODES:, NTASKS_PER_NODE:, WALLTIME:, PARTITION: of the instruction file replace those of
    the profile (so that the #SBATCH lines and the rank count agree), and the VASP ranks are the job's tasks,
    nodes x ntasks_per_node, unless MPI: N is given. With 'srun' (no rank count) VASP runs on every task, so an
    MPI: N below them is run as 'srun -n N' (KPAR/NCORE are chosen for N); above them is an error.
    Returns (profile, note); note starts with ERROR for the error. Without a SLURM block: (prof, '')."""
    if not prof.get('slurm'):
        return prof, ''
    prof = json.loads(json.dumps(prof))
    s = prof['slurm']
    for k, ik in (('nodes', 'slurm_nodes'), ('ntasks_per_node', 'slurm_ntasks_per_node'),
                  ('time', 'slurm_walltime'), ('partition', 'slurm_partition')):
        if inst.get(ik) is not None:
            s[k] = inst[ik]
    nodes, ntpn = int(s.get('nodes', 1)), int(s.get('ntasks_per_node', 1))
    ntask = nodes * ntpn
    n = int(inst.get('mpi_np') or 1)
    if not fx.get('mpi_given'):
        inst['mpi_np'] = prof['mpi_np'] = ntask
        return prof, f'VASP ranks = the job\'s tasks: {nodes} node(s) x {ntpn} = {ntask} (no MPI line)'
    prof['mpi_np'] = n
    if n == ntask:
        return prof, f'VASP ranks: MPI {n} = {nodes} node(s) x {ntpn}'
    if n > ntask:
        return prof, (f'ERROR MPI: {n} is more than the job\'s {nodes} node(s) x {ntpn} = {ntask} tasks: raise NODES / '
                      'NTASKS_PER_NODE (instruction file or profile), or lower MPI')
    cmd = prof.get('mpi_cmd', 'srun').strip()
    if cmd == 'srun':
        prof['mpi_cmd'] = 'srun -n {np}'
    return prof, (f'WARNING: MPI {n} is fewer than the job\'s {ntask} tasks ({nodes} x {ntpn}): VASP runs on {n} '
                  f'({prof["mpi_cmd"].replace("{np}", str(n))}), the other tasks stay idle')


def first_file(*paths):
    """The first of paths that exists (the first one if none does)."""
    return next((p for p in paths if os.path.isfile(os.path.expanduser(p))), paths[0])


def newest_tcond():
    """The newest MPI build ~/ALADYN/bin/tcond-<version>.x, else (no MPI build, e.g. no mpif90) the newest serial one,
    tcond-<version>-serial.x. THERMACOND_X chooses another; a serial build runs on threads (is_serial_tcond)."""
    c = sorted(glob.glob(os.path.expanduser('~/ALADYN/bin/tcond-[0-9]*.x')), key=os.path.getmtime)
    mpi = [f for f in c if 'serial' not in os.path.basename(f)]
    return (mpi or c or ['~/ALADYN/bin/tcond-2.2.x'])[-1]


def is_serial_tcond(exe):
    """The serial (OpenMP, no MPI) build of THERMACOND: one process with threads, not srun/mpirun."""
    return 'serial' in os.path.basename(exe)


def newest_focex():
    """The newest ~/ALADYN/bin/focex-<version>.x (the Makefile's since 8.19), else the newest
    ~/BIN/v<version> (where 8.18 and earlier were installed)."""
    for pat in ('~/ALADYN/bin/focex-[0-9]*.x', '~/BIN/v[0-9]*'):
        c = sorted(glob.glob(os.path.expanduser(pat)), key=os.path.getmtime)
        if c:
            return c[-1]
    return '~/ALADYN/bin/focex-9.1.x'


# tags that fix the magnetic order and +U: given in any INCAR block they apply to every step (relax, SCF, DFPT,
# supercells), so that all of them describe the same system (2026-11-09, user: MAGMOM of an antiferromagnet in INCAR_scf)
PHYS_TAGS = ('MAGMOM', 'ISPIN', 'NUPDOWN', 'LNONCOLLINEAR', 'LSORBIT', 'SAXIS', 'M_CONSTR', 'I_CONSTRAINED_M', 'LAMBDA',
             'LDAU', 'LDAUTYPE', 'LDAUL', 'LDAUU', 'LDAUJ', 'LMAXMIX')
PER_ATOM_TAGS = {'MAGMOM': (1, 3), 'M_CONSTR': (3,)}          # values per atom (collinear, non-collinear)


def phys_lines(fx):
    """The PHYS_TAGS lines of all INCAR blocks (one per tag); the same tag with two values stops the agent."""
    seen = {}
    for b in ('incar', 'incar_relax', 'incar_scf', 'incar_dfpt', 'incar_snaps'):
        for l in incar_lines(fx['blocks'].get(b)):
            t, v = l.split('=', 1)
            t, v = t.strip().upper(), re.split('[!#]', v)[0].strip()
            if t in PHYS_TAGS:
                if t in seen and expand_values(seen[t][0]) != expand_values(v):
                    sys.exit(f'{t} is {seen[t][0]!r} in the {seen[t][1].upper()} block and {v!r} in {b.upper()}: '
                             'the magnetic order and +U must be the same in every step; give it once (INCAR block)')
                seen.setdefault(t, (v, b))
    return [f'{t} = {v}' for t, (v, b) in seen.items()]


def expand_values(v):
    """'2*5.0 -5 0' -> ['5.0', '5.0', '-5', '0'] (VASP's N*value notation)"""
    out = []
    for tok in v.split():
        if '*' in tok:
            n, x = tok.split('*', 1)
            out += [x] * int(n)
        else:
            out.append(tok)
    return out


def compress_values(vals):
    """['5', '5', '-5'] -> '2*5 -5'"""
    out, i = [], 0
    while i < len(vals):
        j = i
        while j < len(vals) and vals[j] == vals[i]:
            j += 1
        out.append(vals[i] if j - i == 1 else f'{j - i}*{vals[i]}')
        i = j
    return ' '.join(out)


def vasp_inst(inst, fx, step, extra=()):
    """VASP-Flow instructions for one step: accuracy defaults, then the user's blocks."""
    d = json.loads(json.dumps(inst))
    acc = ACCURACY[fx['accuracy']]
    d['encut_val'] = inst.get('encut_val') or acc['ENCUT']
    d['ediffg'] = inst.get('ediffg') or acc['EDIFFG']
    b = fx['blocks']
    user_all, user_step = incar_lines(b.get('incar')), incar_lines(b.get('incar_' + step))
    first = [f"PREC = {acc['PREC']}"] + ([] if step == 'dfpt' else [f"EDIFF = {acc['EDIFF']}"])
    step_list = list(extra) + user_step
    if step == 'dfpt' and not has_tag(user_all + user_step, 'EDIFF'):
        step_list = ['EDIFF = 1E-8'] + step_list             # linear response wants it tight
    if step == 'relax' and not has_tag(user_all + user_step, 'IBRION'):
        # quasi-Newton (forces only): the conjugate-gradient line search (IBRION = 2) stopped with 'ZBRENT: fatal
        # error in bracketing' at 0.003 eV/A near the minimum (FePSe3 r2SCAN+U, 2026-11-09; user's choice)
        step_list = ['IBRION = 1'] + step_list
    # VASP-Flow applies 'all' and then the step's list, later lines winning; scf and snaps
    # INCARs are made from the relaxation INCAR, so their list goes under 'relax'
    # the magnetic and +U tags of any block in every step (PHYS_TAGS); the step's own lines still come last
    d['incar_raw'] = {'all': first + user_all + phys_lines(fx),
                      'relax' if step in ('relax', 'scf', 'snaps') else step: step_list}
    kp = b.get('kpoints_snaps' if step == 'snaps' else 'kpoints')
    d['kpoints_raw'] = {'all': kp} if kp else {}
    return d


# ── small helpers ─────────────────────────────────────────────────────────────
def now():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


class Report:
    def __init__(self, proj):
        self.path = os.path.join(proj, 'report.log')

    def __call__(self, msg=''):
        line = f'[{now()}] {msg}' if msg else ''
        print(line)
        with open(self.path, 'a') as fh:
            fh.write(line + '\n')


def write(path, text, mode=None):
    with open(path, 'w') as fh:
        fh.write(text)
    if mode:
        os.chmod(path, mode)


def set_tags(text, tags):
    """Replace or append INCAR tags (tags: dict, value None removes the tag)."""
    out, seen = [], set()
    for line in text.splitlines():
        t = line.split('=', 1)[0].strip().upper() if '=' in line and not line.strip().startswith('#') else None
        if t in tags:
            seen.add(t)
            if tags[t] is not None:
                out.append(f'{t} = {tags[t]}')
            continue
        out.append(line)
    out += [f'{t} = {v}' for t, v in tags.items() if t not in seen and v is not None]
    return '\n'.join(out) + '\n'


def static_incar(text, header):
    """A relaxation INCAR from VASP-Flow, relabelled for a static run."""
    text = text.replace('# Relaxation calculation', header).replace('# Ionic relaxation', '# Static run')
    return '\n'.join(l for l in text.splitlines()
                     if 'WAVECAR and CHGCAR are written for the SCF' not in l
                     and 'Start from scratch (ignore any stale' not in l) + '\n'


def lattice_params(C):
    n = np.linalg.norm(C, axis=1)
    ang = lambda i, j: math.degrees(math.acos(np.clip(C[i] @ C[j] / n[i] / n[j], -1, 1)))
    return [n[0], n[1], n[2], ang(1, 2), ang(0, 2), ang(0, 1)]


def focex_frame(latp):
    """Conventional vectors (rows) exactly as FOCEX/SC_SNAPS build them from a,b,c,al,be,ga."""
    a, b, c = latp[:3]
    al, be, ga = (math.radians(x) for x in latp[3:])
    v1 = [a, 0, 0]
    v2 = [b * math.cos(ga), b * math.sin(ga), 0]
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(ga) * math.cos(be)) / math.sin(ga)
    return np.array([v1, v2, [cx, cy, math.sqrt(c * c - cx * cx - cy * cy)]])


def read_vasp_cell(path):
    """POSCAR/CONTCAR (VASP-5): dict with A (rows, A), species, counts, F, X (Cartesian), types."""
    p = p2f.read_poscar(path)
    p['X'] = p['F'] @ p['A']
    return p


def write_poscar(path, title, A, species, counts, F):
    out = [title, '1.0'] + ['  %.12f %.12f %.12f' % tuple(p2f.clean(v, 12)) for v in A]
    out += ['  ' + '  '.join(species), '  ' + '  '.join(str(c) for c in counts), 'Direct']
    out += ['  %.12f %.12f %.12f' % tuple(p2f.clean(f % 1.0, 12)) for f in F]
    write(path, '\n'.join(out) + '\n')


def cell_data(V, F, types):
    """Conventional cell of the primitive cell V (rows) in its present frame, which must be
    FOCEX's frame (conventional a along x, b in the xy plane)."""
    ds, _ = p2f.dataset(V, F, types)
    P = np.array(p2f.get(ds, 'transformation_matrix'))      # x_conv = P x_prim
    C = np.linalg.inv(P).T @ V
    latp = lattice_params(C)
    prim_in_conv = P.T                                       # row i: primitive vector i in conv units
    assert np.allclose(prim_in_conv @ C, V, atol=1e-8)
    return dict(V=V, F=F % 1.0, C=C, latp=latp, prim_in_conv=prim_in_conv, Xc=(F % 1.0) @ P.T,
                sg=p2f.get(ds, 'international'), in_frame=bool(np.allclose(C, focex_frame(latp), atol=1e-6)))


def input_moments(inst, fx, n):
    """Initial moments per atom of the input POSCAR given in the instruction file: a per-atom MAGMOM: list, or
    MAGMOM = in an INCAR block, with n (collinear) or 3n (non-collinear) values. A list of n tuples, or None when
    no per-atom moments are given (none, one value, or one per element: they do not tell atoms of a species apart)."""
    lists = []
    spec = inst.get('magmom')
    if spec and spec.get('kind') == 'list':
        lists.append(spec['values'])
    for k, body in fx.get('blocks', {}).items():
        m = re.search(r'^\s*MAGMOM\s*=\s*([^!#\n]*)', body, re.M | re.I) if k.lower().startswith('incar') else None
        if m:
            lists.append([float(v) for v in expand_values(m.group(1))])
    for vals in lists:
        for c in (1, 3):
            if len(vals) == c * n:
                return [tuple(vals[c * i:c * i + c]) for i in range(n)]
    return None


def primitive_input(p, moments=None):
    """The input cell reduced to its primitive cell (same orientation, atoms not moved), or p itself when it is
    primitive already. Atoms whose per-atom initial moments differ count as different kinds, so a magnetic order
    that needs the larger cell (an antiferromagnet) keeps it. Given per-atom moments, the cell is never reduced (the
    MAGMOM list belongs to these atoms): if it could be, the note says so. Returns (p, note)."""
    n = len(p['F'])
    types = list(p['types'])
    if moments is not None:
        kinds = sorted(set((t, tuple(round(x, 2) for x in m)) for t, m in zip(types, moments)))
        types = [kinds.index((t, tuple(round(x, 2) for x in m))) for t, m in zip(types, moments)]
    ds, _ = p2f.dataset(p['A'], p['F'], types)
    cls = list(p2f.get(ds, 'mapping_to_primitive'))
    npr = len(set(cls))
    if npr == n:
        return p, ('' if moments is None else f'; {n} atoms kept: the initial moments (MAGMOM) need this cell')
    if moments is not None:
        return p, (f'; {n} atoms kept because MAGMOM is given per atom, but this magnetic order also fits a cell of '
                   f'{npr} atoms: give that cell and its MAGMOM to save time')
    L = np.array(p2f.get(ds, 'primitive_lattice'), dtype=float)
    M = p['A'] @ np.linalg.inv(L)                    # input cell vectors in primitive units: must be integers
    if not np.allclose(M, np.rint(M), atol=1e-4) or round(abs(np.linalg.det(M))) * npr != n:
        sys.exit(f'could not reduce the {n}-atom input cell to its {npr}-atom primitive cell (spglib lattice not '
                 'commensurate); give the primitive cell as POSCAR')
    if np.linalg.det(L) < 0:
        L = -L
    reps = sorted({c: i for i, c in reversed(list(enumerate(cls)))}.values(), key=lambda i: (p['types'][i], i))
    q = dict(p)
    q['A'], q['types'] = L, [p['types'][i] for i in reps]
    q['F'] = (p['X'][reps] @ np.linalg.inv(L)) % 1.0
    q['X'] = q['F'] @ L
    q['counts'] = [q['types'].count(t) for t in range(len(p['species']))]
    return q, f'; reduced from the {n}-atom input cell to its {npr}-atom primitive cell'


def canonical_cell(p):
    """Niggli-reduced primitive cell (poscar2focex.canonical), rotated so that its conventional
    cell is exactly where FOCEX and SC_SNAPS put it from (a,b,c,alpha,beta,gamma)."""
    V, F, _ = p2f.canonical(p['A'], p['F'], p['types'])
    c = cell_data(V, F, p['types'])
    Q = np.linalg.solve(c['C'], focex_frame(c['latp']))      # C @ Q = C_focex
    if not np.allclose(Q @ Q.T, np.eye(3), atol=1e-8) or np.linalg.det(Q) < 0:
        sys.exit('could not rotate the cell into the FOCEX frame (left-handed conventional cell?)')
    V = V @ Q
    V[np.abs(V) < 1e-12] = 0.0
    c = cell_data(V, F, p['types'])
    assert c['in_frame']
    return c


def frame_of(poscar):
    """cell_data of a POSCAR that is already canonical (POSCAR0); stops if it is not."""
    p = read_vasp_cell(poscar)
    c = cell_data(p['A'], p['F'], p['types'])
    if not c['in_frame']:
        sys.exit(f'{poscar} is not in the FOCEX frame of its conventional cell (was it edited?)')
    return p, c


def read_poscar_any(path, species=()):
    """POSCAR with or without the VASP-5 element line (then species from SPECIES:)."""
    L = open(path).read().splitlines()
    s = float(L[1].split()[0])
    A = np.array([[float(x) for x in L[i].split()[:3]] for i in range(2, 5)])
    if s < 0:
        s = (-s / abs(np.linalg.det(A))) ** (1.0 / 3.0)
    A = A * s
    i = 5
    if L[5].split() and not L[5].split()[0][0].isdigit():
        species, i = L[5].split(), 6
    counts = [int(x) for x in L[i].split()[:len(species) or None]]
    if len(species) != len(counts):
        sys.exit(f'{path}: {len(counts)} atom counts but species {list(species)}; give SPECIES: in the instructions')
    i += 1
    if L[i].strip()[:1] in 'Ss':
        i += 1
    cart = L[i].strip()[:1] in 'CcKk'
    X = np.array([[float(x) for x in L[i + 1 + k].split()[:3]] for k in range(sum(counts))])
    F = X * s @ np.linalg.inv(A) if cart else X
    types = [k for k, n in enumerate(counts) for _ in range(n)]
    return {'A': A, 'species': list(species), 'counts': counts, 'F': F, 'X': F @ A, 'types': types}


def ideal_frame(V, F, types):
    """(V0, R, dev): the primitive cell V (rows, frame of the data) on the symmetrized lattice
    (spglib's idealized lattice parameters) in the FOCEX frame, same basis and atom order, and the
    rotation R (det +1) that takes the data onto it (best fit): vectors (rows) v R, tensors
    R^T T R. dev = max |V R - V0| (A), the size of the idealization. An unsymmetrized lattice
    (DFT data have ~1e-6 A) made FOCEX's force-constant symmetry analysis lose groups."""
    ds, _ = p2f.dataset(V, F, types)
    P = np.array(p2f.get(ds, 'transformation_matrix'))              # x_conv = P x_prim
    latp = lattice_params(np.array(p2f.get(ds, 'std_lattice')))      # symmetrized
    V0 = P.T @ focex_frame(latp)
    U, _, Wt = np.linalg.svd(V.T @ V0)                                # Kabsch: min |V R - V0|
    R = U @ Wt
    if np.linalg.det(R) < 0:
        sys.exit('DATA: the cell is left-handed relative to the FOCEX frame')
    V0[np.abs(V0) < 1e-12] = 0.0
    return V0, R, float(np.abs(V @ R - V0).max())


def rotate_forcedisp(src, dst, R):
    """FORCEDISPi: every line of six (or more) numbers is position (or displacement) and force;
    both rotate as row vectors. Other lines (block headers) are copied."""
    out = []
    for line in open(src).read().splitlines():
        t = line.split()
        try:
            v = [float(x) for x in t[:6]]
            ok = len(t) >= 6 and '=' not in line
        except ValueError:
            ok = False
        if ok:
            a, b = np.array(v[:3]) @ R, np.array(v[3:6]) @ R
            line = '  %15.9f %15.9f %15.9f   %15.9f %15.9f %15.9f' % (*a, *b) + \
                   ('  ' + ' '.join(t[6:]) if len(t) > 6 else '')
        out.append(line)
    write(dst, '\n'.join(out) + '\n')


def snap_files(d):
    return sorted(f for f in glob.glob(os.path.join(d, 'poscar_*')) if re.search(r'poscar_\d+$', f))


# ── supercell choice ──────────────────────────────────────────────────────────
def face_distance(L):
    v = abs(np.linalg.det(L))
    return min(v / np.linalg.norm(np.cross(L[(i + 1) % 3], L[(i + 2) % 3])) for i in range(3))


def choose_supercell(c, nprim, kind):
    """m1 x m2 x m3 multiple of the conventional cell: the fewest atoms whose smallest face
    distance is >= D (SMALL 10 A, LARGE 14 A), within the atom limit; otherwise the largest
    face distance within the limit.  Returns (S rows in primitive units, info)."""
    D, nmax = SUPERCELL[kind]
    conv_in_prim = np.rint(np.linalg.inv(c['prim_in_conv'])).astype(int)   # rows: conv vectors
    nconv = int(round(abs(np.linalg.det(conv_in_prim)))) * nprim
    best = None
    for m in __import__('itertools').product(range(1, 9), repeat=3):
        n = nconv * m[0] * m[1] * m[2]
        if n > nmax:
            continue
        L = np.diag(m) @ c['C']
        d = face_distance(L)
        dr = round(d, 4)               # equal distances must tie (ZnO 4x3x2 beat 3x3x2 by 8e-13 A)
        key = (d < D - 1e-6, n if d >= D - 1e-6 else -dr, -dr, n, max(m) - min(m))
        if best is None or key < best[0]:
            best = (key, m, d, L)
    if best is None:
        sys.exit(f'the conventional cell ({nconv} atoms) is larger than the {kind.upper()} limit of {nmax} atoms;'
                 ' give SUPERCELL explicitly')
    _, m, d, L = best
    S = np.diag(m) @ conv_in_prim
    if np.linalg.det(S) < 0:
        S = -S
    info = dict(m=m, n=int(round(abs(np.linalg.det(S)))), natoms=int(round(abs(np.linalg.det(S)))) * nprim,
                face=d, lmin=min(np.linalg.norm(L, axis=1)), reached=d >= D - 1e-6, D=D, nmax=nmax)
    return S, info


def kmesh_from_kpra(g, kpra):
    """Gamma mesh for kpra k-points per reciprocal atom, as VASP-Flow's _kpoints_from_kpra but with
    N_i = s |b_i*| rounded to the nearest integer (minimum 1): no even or multiple-of-6 rounding
    (user, 2026-10-27). Bound to every generator the agent makes (Agent.gen), so it gives the relax/SCF
    KPOINTS and the FOCEX mesh; an explicit KMESH, a KPOINTS block or FOCEX_MESH n1 n2 n3 still win."""
    try:
        A = read_vasp_cell(g.poscar)['A']
        b = np.linalg.norm(2 * np.pi * np.linalg.inv(A).T, axis=1)
    except Exception:
        return (6, 6, 6)
    flat = bool(g.instructions.get('is_2d', False))
    bb = b[:2] if flat else b
    s = (max(1.0, kpra / g._kpra_atoms()) / float(np.prod(bb))) ** (1.0 / len(bb))
    return tuple(max(1, int(math.floor(s * x + 0.5))) for x in bb) + ((1,) if flat else ())


def sc_kmesh(S, V, Nprim):
    """Commensurate Gamma mesh for supercell S (rows, prim units): N/n_i for a diagonal S,
    otherwise the same k-spacing as the primitive mesh."""
    S = np.asarray(S)
    if all(np.count_nonzero(r) == 1 for r in S):             # diagonal, up to a permutation
        return [max(1, math.ceil(Nprim[int(np.flatnonzero(r)[0])] / abs(r).max())) for r in S]
    bp = np.linalg.norm(np.linalg.inv(V).T, axis=1)
    dk = np.prod(bp / np.array(Nprim)) ** (1 / 3)
    bs = np.linalg.norm(np.linalg.inv(S @ V).T, axis=1)
    return [max(1, int(round(b / dk))) for b in bs]


# ── VASP output checks ────────────────────────────────────────────────────────
def outcar_stress(outcar):
    """Stress of the last ionic step of a VASP OUTCAR in FOCEX's convention (stress.params): GPa,
    sigma = dE/(V d eps), tensile > 0. VASP prints -dE/(V d eps) in kB, positive when the cell wants to
    expand (vasp.at wiki, ISIF), components XX YY ZZ XY YZ ZX, taken as printed. Returns (sigma, pstress
    in kB from the 'Pullay stress' line); sigma None without a full tensor (ISIF < 2: trace only)."""
    if not os.path.isfile(outcar):
        return None, 0.0
    t = open(outcar, errors='replace').read()
    m = re.search(r'ISIF\s*=\s*(\d+)', t)
    kb = re.findall(r'^\s*in kB' + r'\s+([-\d.Ee+*]+)' * 6, t, re.M)
    pul = re.findall(r'Pullay stress\s*=\s*([-\d.Ee+]+)', t)
    p = float(pul[-1]) if pul else 0.0
    if not m or int(m.group(1)) < 2 or not kb:
        return None, p
    try:
        xx, yy, zz, xy, yz, zx = (float(x) for x in kb[-1])
    except ValueError:
        return None, p
    return -np.array([[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]]) / 10.0, p


FUNCTIONAL_IGNORED = 'VASP ignored METAGGA'


def needs_vasp6(functional):
    """True for a functional this agent writes as a meta-GGA that VASP 5.x does not have (r2SCAN came with
    VASP 6; VASP 5.4.4 reads METAGGA = R2SCAN without a word and runs PBE)."""
    return str(functional or '').upper().replace('-', '').startswith('R2SCAN')


def functional_ignored(d, t):
    """Message if the INCAR in d asks for a meta-GGA and the OUTCAR text t says VASP ran without it
    (no "METAGGA = <name>  LMAXTAU" line), else None. Found on Rivanna 2026-11-09: VASP 5.4.4 gave the PBE energy for R2SCAN."""
    try:
        m = re.search(r'^\s*METAGGA\s*=\s*([^\s#!]+)', open(os.path.join(d, 'INCAR')).read(), re.M | re.I)
    except OSError:
        return None
    if not m or m.group(1).upper().strip('.') in ('NONE', 'F', 'FALSE', 'PE'):
        return None
    # an active meta-GGA prints 'METAGGA = <name>  LMAXTAU = ...' (SCAN, PBE, R2SCA with the r2SCAN patch); every
    # OUTCAR also has 'METAGGA=  F  non-selfconsistent MetaGGA calc.', which says nothing about it (2026-11-09: the
    # old test matched that line and stopped a real r2SCAN run). An ignored METAGGA (R2SCAN in unpatched 5.4.4) has
    # only the latter.
    if not re.search(r'^\s*METAGGA\s*=\s*\S+\s+LMAXTAU', t, re.M):
        v = re.search(r'vasp\.(\d+\.\d+\.\d+)', t[:2000])
        return (f'{FUNCTIONAL_IGNORED} = {m.group(1)}: the OUTCAR shows no active meta-GGA (no METAGGA = ... '
                f'LMAXTAU line), so this VASP{" " + v.group(1) if v else ""} ran without it (r2SCAN needs VASP 6, or '
                '5.4.4 with the r2SCAN patch). Use FUNCTIONAL: PBE or PBEsol, or such a build (profile vasp_std / '
                'module), then set up again with --force')
    return None


def outcar_status(d, kind):
    """(ok, message) for the OUTCAR in d; kind: relax | static | dfpt."""
    f = os.path.join(d, 'OUTCAR')
    if not os.path.isfile(f):
        return False, 'no OUTCAR (not run)'
    t = open(f, errors='replace').read()
    ign = functional_ignored(d, t)
    if ign:
        return False, ign
    if 'General timing' not in t:
        tail = ''
        vo = os.path.join(d, 'vasp.out')
        if os.path.isfile(vo):
            tail = ' | vasp.out: ' + ' / '.join(open(vo, errors='replace').read().strip().splitlines()[-3:])
        return False, 'VASP did not finish (crash or time limit)' + tail
    if kind == 'relax' and 'reached required accuracy' not in t:
        return False, 'ionic relaxation did not reach EDIFFG within NSW steps'
    if kind in ('static', 'dfpt') and 'EDIFF is reached' not in t:
        return False, ('electronic SCF did not converge within NELM; '
                       'try INCAR_snaps: ALGO = All or a larger NELM')
    if kind == 'dfpt' and 'BORN EFFECTIVE CHARGES' not in t:
        return False, 'no Born effective charges in OUTCAR'
    m = re.search(r'Elapsed time \(sec\):\s*([\d.]+)', t)
    return True, f'converged, {float(m.group(1)) / 60:.1f} min' if m else 'converged'


def band_gap(d):
    """Gap (eV) from EIGENVAL with occupations; None if unreadable."""
    f = os.path.join(d, 'EIGENVAL')
    try:
        L = open(f).read().splitlines()
        ispin = int(L[0].split()[3])
        nel, nk, nb = (int(x) for x in L[5].split()[:3])
        occ, unocc, i = [], [], 7
        for _ in range(nk):
            i += 1
            for _ in range(nb):
                v = [float(x) for x in L[i].split()]
                for s in range(ispin):
                    e, o = v[1 + s], v[1 + ispin + s]
                    (occ if o > 0.5 else unocc).append(e)
                i += 1
            i += 1
        return max(0.0, min(unocc) - max(occ))
    except Exception:
        return None


def read_born(outcar, natom, transpose):
    """(eps_inf 3x3, [Z 3x3 per atom]) from a LEPSILON / LCALCEPS OUTCAR, in FOCEX's order
    charge(i,j): i = force/displacement, j = field (aux_na.f90 contracts q with j)."""
    t = open(outcar, errors='replace').read()
    num = r'([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)'
    eps, rank = None, -1
    for m in re.finditer(r'MACROSCOPIC STATIC DIELECTRIC TENSOR([^\n]*)\n\s*-+\n'
                         + r'\s*' + num + r'\s*\n\s*' + num + r'\s*\n\s*' + num, t):
        hdr = m.group(1).upper()
        if 'IONIC' in hdr:
            continue                                  # epsilon_inf is the electronic part
        r = 1 if 'INCLUDING LOCAL FIELD' in hdr else 0
        if r >= rank:
            eps, rank = np.array([float(x) for x in m.groups()[1:]]).reshape(3, 3), r
    blocks = list(re.finditer(r'BORN EFFECTIVE CHARGES.*?(?=BORN EFFECTIVE CHARGES|\Z)', t, re.S))
    if eps is None or not blocks:
        return None, None
    rows = re.findall(r'^\s*[123]\s+' + num, blocks[-1].group(), re.M)
    if len(rows) < 3 * natom:
        return None, None
    Z = [np.array(rows[3 * k:3 * k + 3], dtype=float) for k in range(natom)]
    return eps, [z.T if transpose else z for z in Z]


class ExternalRun(Exception):
    """VASP_RUN / THERMACOND_RUN external: the folders dirs are ready to run outside the agent; after maps a
    folder to the one it must wait for (DFPT after the SCF)."""
    def __init__(self, dirs, after=None, kind='vasp'):
        super().__init__(kind)
        self.dirs, self.after, self.kind = list(dirs), dict(after or {}), kind


def copy_if_changed(src, dst_dir):
    """Copy src into dst_dir unless an identical file is there (keeps its time stamp)."""
    dst = os.path.join(dst_dir, os.path.basename(src))
    if not (os.path.isfile(dst) and open(src, 'rb').read() == open(dst, 'rb').read()):
        shutil.copy(src, dst)


def write_if_changed(path, text):
    if not (os.path.isfile(path) and open(path).read() == text):
        write(path, text)


# ── the agent ─────────────────────────────────────────────────────────────────
class Agent:
    def __init__(self, proj):
        self.proj = os.path.abspath(proj)
        self.rep = Report(self.proj)
        st = os.path.join(self.proj, 'state.json')
        self.state = json.load(open(st)) if os.path.isfile(st) else {}
        self.inst, self.fx = parse_instructions(self.instr_file())
        self.profile = vagent.load_profile(self.fx['profile'] or self.state.get('profile') or '')
        self.profile, self.mpi_note = slurm_setup(self.inst, self.profile, self.fx)
        if self.mpi_note.startswith('ERROR'):
            sys.exit(self.mpi_note[6:])

    def instr_file(self):
        f = os.path.join(self.proj, INSTR)
        return f if os.path.isfile(f) or not os.path.isfile(os.path.join(self.proj, OLD_INSTR)) \
            else os.path.join(self.proj, OLD_INSTR)

    def save(self):
        write(os.path.join(self.proj, 'state.json'), json.dumps(self.state, indent=1, default=str))

    def d(self, *a):
        return os.path.join(self.proj, *a)

    def gen(self, poscar, step, extra=()):
        g = VASPInputGenerator(poscar, vasp_inst(self.inst, self.fx, step, extra), profile=self.profile)
        g._kpoints_from_kpra = lambda kpra: kmesh_from_kpra(g, kpra)      # no parity rounding
        return g

    def vasp(self, d, gen):
        """Run VASP in d; an existing OUTCAR is kept as OUTCAR.1, OUTCAR.2, ... (never overwritten).
        VASP_RUN external: the inputs are ready, nothing is run (ExternalRun)."""
        if os.path.isfile(os.path.join(d, 'OUTCAR')):
            k = 1
            while os.path.exists(os.path.join(d, f'OUTCAR.{k}')):
                k += 1
            shutil.move(os.path.join(d, 'OUTCAR'), os.path.join(d, f'OUTCAR.{k}'))
        if self.fx['vasp_run'] == 'external':          # also with the old OUTCAR set aside, so that it is
            raise ExternalRun([d])                     # not read again as the result
        cmd = gen._get_mpi_cmd(gen._get_vasp_exec())
        t0 = time.time()
        subprocess.run(f'{cmd} > vasp.out 2>&1', shell=True, cwd=d)
        return time.time() - t0

    def link_potcar(self, d):
        vagent.link_potcar(d, self.d('POTCAR'))

    @staticmethod
    def clean(d, keep=()):
        for f in ('WAVECAR', 'CHG', 'CHGCAR', 'WFULL', 'TMPCAR', 'PROCAR'):
            if f not in keep and os.path.exists(os.path.join(d, f)):
                os.remove(os.path.join(d, f))

    def is_metal(self, species):
        return all(s.split('_')[0] in METALS for s in species)

    # ---- task 1: relaxation ------------------------------------------------
    def t_relax(self):
        d = self.d('01_relax')
        if not os.path.isfile(os.path.join(d, 'INCAR')):
            self.write_relax()
        for cycle in range(1, 4):
            ok, msg = outcar_status(d, 'relax')
            if ok:
                return True, msg
            if msg.startswith(FUNCTIONAL_IGNORED):              # a restart cannot help
                return False, msg
            if os.path.isfile(os.path.join(d, 'OUTCAR')):        # restart from CONTCAR
                if os.path.isfile(os.path.join(d, 'CONTCAR')) and os.path.getsize(os.path.join(d, 'CONTCAR')) > 0:
                    shutil.copy(os.path.join(d, 'CONTCAR'), os.path.join(d, 'POSCAR'))
                self.rep(f'  relax: {msg}; restart {cycle} from CONTCAR')
            dt = self.vasp(d, self.gen(self.d('POSCAR.init'), 'relax'))
            self.rep(f'  relax run {cycle}: {dt / 60:.1f} min')
        return outcar_status(d, 'relax')

    def write_relax(self):
        d = self.d('01_relax')
        self.gen(self.d('POSCAR.init'), 'relax').generate_relax_input(d)
        self.link_potcar(d)

    # ---- task 2: POSCAR0, SCF, DFPT -----------------------------------------
    def t_poscar0(self):
        cont = self.d('01_relax', 'CONTCAR')
        p = read_vasp_cell(cont)
        c = canonical_cell(p)
        write_poscar(self.d('POSCAR0'), f"{self.fx['project']} relaxed, Niggli-reduced, FOCEX frame",
                     c['V'], p['species'], p['counts'], c['F'])
        G0, G1 = p['A'] @ p['A'].T, c['V'] @ c['V'].T
        same = np.allclose(G0, G1, atol=1e-6) and np.allclose((c['F'] - p['F'] + 0.5) % 1 - 0.5, 0, atol=1e-6)
        self.state['reuse_relax'] = bool(same)
        self.state['sg'] = c['sg']
        return True, (f"space group {c['sg']}; a b c = {' '.join('%.5f' % x for x in c['latp'][:3])} A, "
                      f"angles {' '.join('%.3f' % x for x in c['latp'][3:])}; basis "
                      + ('unchanged (CHGCAR/WAVECAR reused)' if same else 'changed (SCF starts fresh)'))

    def t_scf(self):
        d = self.d('02_scf')
        os.makedirs(d, exist_ok=True)
        ok, msg = outcar_status(d, 'static')
        if ok:
            return True, self.scf_done(d, msg)
        # the relaxation's CHGCAR/WAVECAR only when they are there (a relaxed POSCAR given, or a
        # relaxation run elsewhere, has none)
        reuse = self.state.get('reuse_relax') and os.path.isfile(self.d('01_relax', 'CHGCAR'))
        g = self.gen(self.d('POSCAR0'), 'scf', ['IBRION = -1', 'NSW = 0', 'ISIF = 1'])
        incar = g._generate_incar_relax()
        incar = set_tags(incar, {'ISTART': 1 if reuse else 0, 'ICHARG': 1 if reuse else 2,
                                 'LWAVE': '.TRUE.', 'LCHARG': '.TRUE.'})
        write(os.path.join(d, 'INCAR'), static_incar(incar, '# SCF on POSCAR0 (focex-agent)'))
        shutil.copy(self.d('01_relax', 'KPOINTS'), os.path.join(d, 'KPOINTS'))
        shutil.copy(self.d('POSCAR0'), os.path.join(d, 'POSCAR'))
        self.link_potcar(d)
        if reuse:
            for f in ('CHGCAR', 'WAVECAR'):
                if os.path.isfile(self.d('01_relax', f)):
                    shutil.copy(self.d('01_relax', f), d)
        if self.fx['vasp_run'] == 'external' and self.fx['dfpt'] != 'off':
            # the DFPT inputs go with the SCF, to run after it (its CHGCAR stays where it ran); with
            # DFPT: auto the job skips itself when the SCF shows no gap (bundle.py)
            dd = self.dfpt_dir()
            self.write_dfpt(dd)
            raise ExternalRun([d, dd], after={dd: d})
        self.vasp(d, g)
        ok, msg = outcar_status(d, 'static')
        if ok:
            self.clean(self.d('01_relax'))
            msg = self.scf_done(d, msg)
        return ok, msg

    def scf_done(self, d, msg):
        """Band gap of the converged SCF in d, into state['gap'] and the message."""
        gap = band_gap(d)
        self.state['gap'] = gap
        msg += f'; band gap {gap:.3f} eV' if gap is not None else '; band gap unreadable'
        species = read_vasp_cell(self.d('POSCAR0'))['species']
        if gap is not None and gap < 0.05 and not self.is_metal(species) and not self.fx['blocks'].get('kpoints'):
            msg += ' -- WARNING: metallic, but the coarse k-mesh was used (set KMESH_DENSITY: fine)'
        return msg

    def want_dfpt(self):
        m = self.fx['dfpt']
        return m == 'on' or (m == 'auto' and (self.state.get('gap') or 0) > 0.05)

    def legacy_dir(self, new, old):
        """Folder new; old when only old exists (projects made before 2026-11-07)."""
        n, o = self.d(new), self.d(old)
        return o if os.path.isdir(o) and not os.path.isdir(n) else n

    def dfpt_dir(self):
        return self.legacy_dir('06_dfpt', os.path.join('02_scf', 'dfpt'))      # 06_dfpt as in VASP-Flow

    def tc_dir(self):
        return self.legacy_dir('07_thermacond', '06_thermacond')

    def safe_efield(self, d):
        """EFIELD_PEAD (eV/A) for LCALCEPS: half of VASP's Zener limit e|E||a_i| < E_g/(10 N_i), with the SCF
        gap, the k-mesh of d/KPOINTS and the lattice of POSCAR0; at most VASP's default 0.01. None when the
        gap or the mesh cannot be read."""
        gap = self.state.get('gap')
        try:
            L = open(os.path.join(d, 'KPOINTS')).read().splitlines()
            N = [int(x) for x in L[3].split()[:3]]
        except (OSError, IndexError, ValueError):
            return None
        if not gap or gap <= 0 or L[2].strip()[:1].upper() not in 'GM':
            return None
        A = read_vasp_cell(self.d('POSCAR0'))['A']
        e = 0.5 * min(gap / (10 * n * np.linalg.norm(a)) for n, a in zip(N, A))
        return min(0.01, float('%.2g' % e))

    @staticmethod
    def incar_value(d, tag):
        m = None
        for l in open(os.path.join(d, 'INCAR')).read().splitlines():
            l = re.split('[!#]', l)[0]
            if re.match(r'\s*' + tag + r'\s*=', l, re.I):
                m = l.split('=', 1)[1].split()
        return m

    def set_efield(self, d, e, why):
        """EFIELD_PEAD = e e e in d/INCAR (replacing an earlier line)."""
        f = os.path.join(d, 'INCAR')
        L = [l for l in open(f).read().splitlines() if not re.match(r'\s*EFIELD_PEAD\s*=', l, re.I)]
        write(f, '\n'.join(L + [f'EFIELD_PEAD = {e:g} {e:g} {e:g}   ! focex-agent: {why}']) + '\n')

    def t_dfpt(self):
        d = self.dfpt_dir()
        rel = os.path.relpath(d, self.proj)
        if not self.want_dfpt():
            self.clean(self.d('02_scf'))
            return True, 'skipped (no gap)' if self.fx['dfpt'] == 'auto' else 'skipped (DFPT: off)'
        ok, msg = outcar_status(d, 'dfpt')
        user_field = has_tag(incar_lines(self.fx['blocks'].get('incar')) +
                             incar_lines(self.fx['blocks'].get('incar_dfpt')), 'EFIELD_PEAD')
        if not ok:
            g = self.write_dfpt(d)
            subprocess.run('bash copy_from_scf.sh', shell=True, cwd=d, capture_output=True)
            self.vasp(d, g)
            ok, msg = outcar_status(d, 'dfpt')
        if not ok:
            return ok, msg
        n = sum(read_vasp_cell(self.d('POSCAR0'))['counts'])
        eps, Z = read_born(os.path.join(d, 'OUTCAR'), n, self.fx['born_transpose'])
        if eps is None:
            return False, 'converged, but epsilon/Born charges could not be read from OUTCAR'
        bad = self.dfpt_checks(d, eps, Z)
        e = self.safe_efield(d)
        field = self.incar_value(d, 'EFIELD_PEAD')
        if any('Zener' in b for b in bad) and e and not user_field and \
                (field is None or float(field[0]) > 1.01 * e):
            # VASP's own remedy: a smaller field; the earlier OUTCAR is kept as OUTCAR.<k>
            self.rep(f'  dfpt: {"; ".join(bad)}; rerun with EFIELD_PEAD = {e:g} eV/A')
            self.set_efield(d, e, 'half of the Zener limit, after a run with Zener warnings')
            subprocess.run('bash copy_from_scf.sh', shell=True, cwd=d)
            self.vasp(d, self.gen(self.d('POSCAR0'), 'dfpt'))
            ok, msg = outcar_status(d, 'dfpt')
            if not ok:
                return ok, msg + ' (rerun with a smaller EFIELD_PEAD)'
            eps, Z = read_born(os.path.join(d, 'OUTCAR'), n, self.fx['born_transpose'])
            if eps is None:
                return False, 'converged, but epsilon/Born charges could not be read from OUTCAR'
            bad = self.dfpt_checks(d, eps, Z)
        msg += '; eps_inf diag ' + ' '.join('%.3f' % x for x in np.diag(eps))
        msg += self.born_trace(d, rel)
        self.clean(d)
        self.clean(self.d('02_scf'))
        if bad:
            # not passed on: FOCEX then runs with born_flag 0 (unless a dielectric_params block is given)
            self.state.pop('eps', None)
            self.state.pop('Z', None)
            return True, (msg + ' -- WARNING: ' + '; '.join(bad) + '. Born charges and epsilon NOT used '
                          '(born_flag 0); give them in a dielectric_params block, or fix the DFPT run in '
                          + rel + ' and run --only dfpt')
        self.state['eps'] = eps.tolist()
        self.state['Z'] = [z.tolist() for z in Z]
        return True, msg + '; checks passed'

    def write_dfpt(self, d):
        """The DFPT inputs in d (VASP-Flow's), with EFIELD_PEAD from the SCF gap when it is known; without it
        (SCF not run yet) VASP-Flow's set_efield.py sets it when the job starts."""
        user_field = has_tag(incar_lines(self.fx['blocks'].get('incar')) +
                             incar_lines(self.fx['blocks'].get('incar_dfpt')), 'EFIELD_PEAD')
        g = self.gen(self.d('POSCAR0'), 'dfpt')
        g.generate_dfpt_input(d, self.d('02_scf'))
        self.link_potcar(d)
        e = self.safe_efield(d)
        if e and not user_field and (self.incar_value(d, 'LCALCEPS') or [''])[0].upper().startswith('.T'):
            self.set_efield(d, e, 'half of the Zener limit E_g/(10 N_i |a_i|)')
        return g

    def dfpt_checks(self, d, eps, Z):
        """What is wrong with a DFPT result: VASP's Zener warning (EFIELD_PEAD too large for the gap and
        the k-mesh), charges that do not sum to zero (> 0.1 e in a component; acoustic sum rule), eps_inf
        not between 1 and 100. Also counts VASP's PEAD warnings (Berry-phase strings not well clustered)."""
        bad = []
        t = ''.join(open(os.path.join(d, f), errors='replace').read()
                    for f in ('vasp.out', 'OUTCAR') if os.path.isfile(os.path.join(d, f)))
        if 'Zener' in t:
            bad.append('VASP warns of Zener tunneling (EFIELD_PEAD too large for the gap and the k-mesh)')
        asr = float(np.abs(np.sum(np.array(Z), axis=0)).max())
        if asr > 0.1:
            bad.append('the Born charges sum to %.3f e, not 0' % asr)
        w = np.linalg.eigvalsh(0.5 * (eps + eps.T))
        if w.min() <= 1 or w.max() > 100:
            bad.append('eps_inf eigenvalues %.3g..%.3g, outside 1..100' % (w.min(), w.max()))
        npead = t.count('not well clustered')
        if bad and npead:
            bad.append(f'{npead} PEAD warnings (Berry-phase strings not well clustered)')
        return bad

    def born_trace(self, d, rel):
        """The agent runs VASP itself, not run.sh, whose last line calls extract_born.py: run it here
        (-> born_charges.txt, BORN, extract_born.out in d)."""
        if not os.path.isfile(os.path.join(d, 'extract_born.py')):
            return ''
        r = subprocess.run([sys.executable, 'extract_born.py', 'OUTCAR'], cwd=d, capture_output=True, text=True)
        write(os.path.join(d, 'extract_born.out'), r.stdout + r.stderr)
        if r.returncode == 0 and os.path.isfile(os.path.join(d, 'born_charges.txt')):
            return f'; born_charges.txt written in {rel}'
        return f'; extract_born.py failed (see {rel}/extract_born.out)'

    # ---- task 3: snapshots ---------------------------------------------------
    def t_snaps(self):
        d = self.d('03_snaps')
        os.makedirs(d, exist_ok=True)
        p, c = frame_of(self.d('POSCAR0'))
        if not os.path.isfile(self.fx['sc_snaps_x']):
            return False, f"sc_snaps.x not found at {self.fx['sc_snaps_x']} (set SC_SNAPS_X)"
        b = self.fx['blocks']
        n = len(p['F'])
        if isinstance(self.fx['supercell'], list):
            S, info = np.array(self.fx['supercell']), None
            if np.linalg.det(S) < 0.5:
                return False, 'SUPERCELL matrix must have a positive determinant'
        else:
            S, info = choose_supercell(c, n, self.fx['supercell'])
        self.state['S'] = S.tolist()
        write(os.path.join(d, 'cell.inp'), b.get('cell_inp') or self.cell_inp(c, p))
        write(os.path.join(d, 'supercell.inp'), b.get('supercell_inp')
              or '\n'.join(' '.join(str(x) for x in r) for r in S) + '   # rows: supercell vectors, primitive units\n')
        T = self.fx['temperature']
        lines = []
        if info is None:
            lines.append(f"supercell from SUPERCELL (rows, primitive units): {S.tolist()}, "
                         f"{int(round(np.linalg.det(S))) * n} atoms")
        if info:
            lines.append(f"supercell {self.fx['supercell'].upper()}: {'x'.join(map(str, info['m']))} conventional, "
                         f"{info['n']} primitive cells, {info['natoms']} atoms, "
                         f"min face distance {info['face']:.2f} A (target {info['D']}), shortest vector "
                         f"{info['lmin']:.2f} A" + ('' if info['reached'] else
                         f" -- target not reached within {info['nmax']} atoms"))
        for attempt in range(3):
            write(os.path.join(d, 'snaps.inp'), b.get('snaps_inp') or
                  f"{self.fx['avg_freq']:g}   # average frequency (1/cm)\n{T:g}   # temperature (K)\n"
                  f"{self.fx['nsnap']}   # number of snapshots\n")
            for f in snap_files(d):
                os.remove(f)
            r = subprocess.run(self.fx['sc_snaps_x'], cwd=d, capture_output=True, text=True)
            if not os.path.isfile(os.path.join(d, 'poscar_000')):
                return False, f"sc_snaps.x produced no poscar_000: {(r.stdout + r.stderr)[-400:]}"
            warn = [l.strip() for l in r.stdout.splitlines() if 'WARNING' in l]
            if warn and warn[0] not in lines:
                lines.append(warn[0])
            bad = [os.path.basename(f) for f in snap_files(d) if '*' in open(f).read()]
            if bad:
                return False, (f"sc_snaps.x wrote overflowing positions (***) in {', '.join(bad[:3])}: "
                               'displacements too large for its format (zero-frequency modes in its model; SC-SNAPS before '
                               '2026-10-23 missed periodic images in small supercells)')
            st = self.disp_stats(d)
            lines.append(f"T = {T:g} K: <u> = {' '.join('%+.4f' % x for x in st['mean'])} A, "
                         f"sigma = {' '.join('%.4f' % x for x in st['std'])} A (x y z), rms {st['rms']:.4f} A, "
                         f"mean |u| = {st['norm']:.4f} A; per element sigma: "
                         + ', '.join(f'{k} {v:.4f}' for k, v in st['elem'].items()))
            if b.get('snaps_inp') or self.fx['temperature_given']:
                lines.append('  temperature ' + ('of the snaps_inp block' if b.get('snaps_inp') else
                             'given in TEMPERATURE') + ' kept (no retuning to DISP_TARGET)')
                break
            if abs(st['rms'] / self.fx['disp_target'] - 1) < 0.15:
                break
            T = T * (self.fx['disp_target'] / st['rms']) ** 2
            lines.append(f'  rms differs from DISP_TARGET {self.fx["disp_target"]} A: rerun at T = {T:.1f} K')
        self.state['T_snaps'] = T
        ok, msg = self.check_commensurate(d, c, S)
        lines.append(msg)
        return ok, '\n  '.join(lines)

    def add_snapshots(self):
        """SNAPSHOTS raised after the snapshots were made: the extra ones are made by sc_snaps.x in
        03_snaps/added_NNN-MMM (same cell.inp, supercell.inp, snaps.inp but the count) and copied to
        03_snaps as poscar_NNN..MMM; the existing ones are kept. (ok, message), message '' if nothing to do."""
        d = self.d('03_snaps')
        have = len(snap_files(d)) - 1                       # poscar_000 is the undistorted cell
        want = self.fx['nsnap']
        if want < have:
            return True, (f'SNAPSHOTS {want} < the {have} in 03_snaps: all kept '
                          '(move folders to 04_supercells/unused/ to leave some out)')
        if want == have:
            return True, ''
        if want > 999:
            return False, 'SNAPSHOTS > 999: sc_snaps.x numbers the snapshots with 3 digits'
        new = os.path.join(d, f'added_{have + 1:03d}-{want:03d}')
        os.makedirs(new, exist_ok=True)
        for f in ('cell.inp', 'supercell.inp'):
            shutil.copy(os.path.join(d, f), new)
        inp = open(os.path.join(d, 'snaps.inp')).read().splitlines()
        inp[2] = f'{want - have}   # number of snapshots (added to the {have} in 03_snaps)'
        write(os.path.join(new, 'snaps.inp'), '\n'.join(inp) + '\n')
        for f in snap_files(new):
            os.remove(f)
        r = subprocess.run(self.fx['sc_snaps_x'], cwd=new, capture_output=True, text=True)
        made = snap_files(new)[1:]
        if len(made) != want - have:
            return False, f'sc_snaps.x made {len(made)} snapshots in {new}, not {want - have}: {(r.stdout + r.stderr)[-400:]}'
        if any('*' in open(f).read() for f in made):
            return False, f'sc_snaps.x wrote overflowing positions (***) in {new}'
        for k, f in enumerate(made):
            shutil.copy(f, os.path.join(d, f'poscar_{have + 1 + k:03d}'))
        T = inp[1].split()[0]
        return True, (f'SNAPSHOTS {want} > the {have} in 03_snaps: {want - have} more made at T = {T} K '
                      f'(poscar_{have + 1:03d}..{want:03d}, from {os.path.basename(new)}); the {have} kept')

    def cell_inp(self, c, p):
        latp, pc = c['latp'], c['prim_in_conv']
        masses = [p2f.MASS[s.split('_')[0]] for s in p['species']]
        out = [' '.join('%.10f' % x for x in latp) + '   # a b c alpha beta gamma (conventional)',
               ' '.join('%.12g' % x for x in p2f.clean(pc, 12).ravel()) + '   # primitive vectors in conventional units',
               '1.0 1   # scale, convcoord=1: atoms and supercell.inp in primitive units',
               f"{len(p['species'])}", ' '.join(str(x) for x in p['counts']),
               ' '.join('%.4f' % m for m in masses) + ' ' + ' '.join('0' for _ in masses) + '   # masses, charges',
               ' '.join(p['species'])]
        out += ['%.12f %.12f %.12f' % tuple(p2f.clean(f, 12)) for f in c['F']]
        return '\n'.join(out) + '\n'

    def disp_stats(self, d):
        p0 = read_vasp_cell(os.path.join(d, 'poscar_000'))
        A, inv = p0['A'], np.linalg.inv(p0['A'])
        us = []
        for f in snap_files(d)[1:]:
            x = read_vasp_cell(f)['X']
            df = (x - p0['X']) @ inv
            us.append((df - np.rint(df)) @ A)
        U = np.concatenate(us)
        el = np.array([p0['species'][t] for t in p0['types']] * len(us))
        return dict(mean=U.mean(0), std=U.std(0), rms=float(np.sqrt((U ** 2).mean())),
                    norm=float(np.linalg.norm(U, axis=1).mean()),
                    elem={s: float(np.sqrt((U[el == s] ** 2).mean())) for s in p0['species']})

    def check_commensurate(self, d, c, S):
        """FOCEX's own test, before any supercell VASP job: rs_i integer in r0_i (same frame),
        and every atom of poscar_000 on a POSCAR0 site."""
        p0 = read_vasp_cell(os.path.join(d, 'poscar_000'))
        V = c['V']
        n = p0['A'] @ np.linalg.inv(V)
        if not np.allclose(n, np.rint(n), atol=1e-5):
            return False, f'poscar_000 is NOT commensurate with POSCAR0 in this frame: rs in r0 units =\n{n}'
        X0 = c['F'] @ V
        for x in p0['X']:
            f = (x - X0) @ np.linalg.inv(V)
            if not np.any(np.all(np.abs(f - np.rint(f)) < 1e-4, axis=1)):
                return False, f'an atom of poscar_000 at {x} is not on a POSCAR0 site'
        return True, (f"commensurability checked: supercell = {np.rint(n).astype(int).tolist()} "
                      f"(rows, primitive units), {len(p0['X'])} atoms all on POSCAR0 sites")

    # ---- task 4: supercell folders ------------------------------------------
    SNAP_INPUTS = ('POSCAR', 'INCAR', 'KPOINTS', 'POTCAR')

    def snap_check(self):
        """(expected, missing, extra) snapshot ids (e.g. '007'). Expected: every 03_snaps/poscar_NNN except
        those whose folder was set aside in 04_supercells/unused/ (a deliberate limit, as limit_to_20.sh
        did); missing: expected without a complete 04_supercells/snap_NNN (POSCAR, INCAR, KPOINTS, POTCAR);
        extra: snap_NNN folders that are not expected."""
        top = self.d('04_supercells')
        ids = [f.rsplit('_', 1)[1] for f in snap_files(self.d('03_snaps'))]
        aside = {os.path.basename(d)[5:] for d in glob.glob(os.path.join(top, 'unused', 'snap_*'))}
        expected = [i for i in ids if i not in aside]
        have = {os.path.basename(d)[5:] for d in glob.glob(os.path.join(top, 'snap_*')) if os.path.isdir(d)}
        missing = [i for i in expected
                   if not all(os.path.exists(os.path.join(top, 'snap_' + i, f)) for f in self.SNAP_INPUTS)]
        extra = sorted(have - set(expected))
        return expected, missing, extra

    def snap_dirs(self):
        """The expected snapshot folders, or (None, message) if some are missing or unexpected ones exist."""
        expected, missing, extra = self.snap_check()
        if missing or extra:
            msg = []
            if missing:
                msg.append(f'{len(missing)} snapshot folder(s) missing or incomplete in 04_supercells: '
                           f'{" ".join("snap_" + i for i in missing[:10])}{" ..." if len(missing) > 10 else ""}'
                           ' (the supercells task recreates them)')
            if extra:
                msg.append(f'folder(s) without a snapshot in 03_snaps: {" ".join("snap_" + i for i in extra)} '
                           '(remove them, or move them to 04_supercells/unused/)')
            return None, '; '.join(msg)
        return [self.d('04_supercells', 'snap_' + i) for i in expected], ''

    def t_supercells(self):
        src, top = self.d('03_snaps'), self.d('04_supercells')
        os.makedirs(top, exist_ok=True)
        expected, missing, extra = self.snap_check()
        if extra:
            return False, f'folder(s) without a snapshot in 03_snaps: {" ".join("snap_" + i for i in extra)}'
        snaps = [os.path.join(src, 'poscar_' + i) for i in expected]
        S = np.array(self.state['S'])
        V = read_vasp_cell(self.d('POSCAR0'))['A']
        kind, Np = VASPInputGenerator._parse_kpoints_text(open(self.d('02_scf', 'KPOINTS')).read())
        if kind != 'mesh' and not self.fx['blocks'].get('kpoints_snaps'):
            return False, '02_scf/KPOINTS is not an automatic mesh: give a KPOINTS_snaps block'
        mesh = sc_kmesh(S, V, Np)
        # ISYM=0, time reversal only: k and -k are one point, except the points with k = -k (2 per
        # even mesh direction of a Gamma mesh); VASP 5.4.4 IBZKPT counts on Si for 1,2,3,4,5 and mixed
        # meshes agree (2x2x2: 8, not (8+1)//2 = 4)
        ntrim = int(np.prod([2 if m_ % 2 == 0 else 1 for m_ in mesh]))
        nk = max(1, (int(np.prod(mesh)) + ntrim) // 2)
        base = read_vasp_cell(snaps[0])
        np_ = self.inst.get('mpi_np') or 1
        extra = ['IBRION = -1', 'NSW = 0', 'ISIF = 1', 'ISYM = 0', 'LWAVE = .FALSE.', 'LCHARG = .FALSE.']
        inst_k = {}
        if np_ > 1 and not (self.inst.get('kpar') or self.inst.get('ncore')):
            kpar, ncore = VASPInputGenerator._auto_kpar_ncore(np_, nk)
            inst_k = {'relax_kpar': kpar, 'relax_ncore': ncore}
        kp = self.fx['blocks'].get('kpoints_snaps') or (
            f'Gamma mesh commensurate with the primitive {Np} (focex-agent)\n0\nGamma\n'
            f'{mesh[0]} {mesh[1]} {mesh[2]}\n0 0 0\n')
        for f in snaps:
            i = f.rsplit('_', 1)[1]
            d = os.path.join(top, 'snap_' + i)
            os.makedirs(d, exist_ok=True)
            p = read_vasp_cell(f)
            write_poscar(os.path.join(d, 'POSCAR'), f'snapshot {i} ({self.fx["project"]})',
                         p['A'], p['species'], p['counts'], p['F'])
            write(os.path.join(d, 'KPOINTS'), kp)
            self.link_potcar(d)
        g = self.gen(os.path.join(top, 'snap_000', 'POSCAR'), 'snaps', extra)
        g.instructions.update(inst_k)
        self.magmom_for_supercell(g, base)
        incar = static_incar(g._generate_incar_relax(), '# supercell snapshot, static (focex-agent)')
        if not has_tag(incar_lines(self.fx['blocks'].get('incar_snaps')), 'LORBIT'):
            incar = set_tags(incar, {'LORBIT': None, 'EDIFFG': None, 'ISTART': 0, 'ICHARG': 2})
        incar, note = self.per_atom_to_supercell(incar, base)
        if note.startswith('ERROR'):
            return False, note[6:]
        note += self.magnetic_xyz(incar, [os.path.join(top, 'snap_' + f.rsplit('_', 1)[1], 'POSCAR') for f in snaps])
        for f in snaps:
            i = f.rsplit('_', 1)[1]
            # the undistorted supercell gives the stress of the reference: full tensor (ISIF=1 is the trace only)
            write(os.path.join(top, 'snap_' + i, 'INCAR'), set_tags(incar, {'ISIF': 2}) if i == '000' else incar)
        made = (f' ({len(missing)} missing or incomplete, written now)' if 0 < len(missing) < len(snaps) else '')
        aside = len(snap_files(src)) - len(snaps)
        return True, (f'{len(snaps)} folders{made}' + (f', {aside} set aside in unused/' if aside else '') +
                      f'; KPOINTS {mesh} from the primitive {Np}; '
                      f'{sum(base["counts"])} atoms' + (f'; KPAR {inst_k["relax_kpar"]} NCORE {inst_k["relax_ncore"]}'
                                                       if inst_k else '') + note)

    def per_atom_to_supercell(self, incar, sc):
        """MAGMOM (and M_CONSTR) of the supercell INCAR: values given per atom of POSCAR0 (1 or 3 per atom) are
        carried to every supercell atom through its POSCAR0 site, as magmom_for_supercell does for the MAGMOM
        keyword; a list of the supercell's length is kept; any other length is an error. LDAUL/U/J must have one
        value per species of the supercell (its species are those of POSCAR0). Returns (incar, note)."""
        p = read_vasp_cell(self.d('POSCAR0'))
        npr, nsc = len(p['F']), len(sc['X'])
        inv = np.linalg.inv(p['A'])
        site = []
        for x in sc['X']:
            f = (x - p['F'] @ p['A']) @ inv
            site.append(int(np.argmin(np.abs(f - np.rint(f)).sum(1))))
        notes, tags = [], {}
        for l in incar.splitlines():
            m = re.match(r'\s*([A-Za-z_]+)\s*=\s*([^!#]*)', l)
            if m:
                tags[m.group(1).upper()] = m.group(2).strip()
        for t, comps in PER_ATOM_TAGS.items():
            if t not in tags:
                continue
            vals = expand_values(tags[t])
            c = next((c for c in comps if len(vals) == c * npr), None)
            if c is not None and npr != nsc:
                new = [v for k in site for v in vals[c * k:c * k + c]]
                incar = set_tags(incar, {t: compress_values(new)})
                notes.append(f'{t} of the {npr} POSCAR0 atoms carried to the {nsc} supercell atoms')
            elif not any(len(vals) == c * nsc for c in comps):
                return incar, (f'ERROR {t} has {len(vals)} values: neither {" or ".join(str(c * npr) for c in comps)} '
                               f'(POSCAR0, {npr} atoms) nor {" or ".join(str(c * nsc) for c in comps)} (supercell)')
        nsp = len(sc['species'])
        for t in ('LDAUL', 'LDAUU', 'LDAUJ'):
            if t in tags and len(expand_values(tags[t])) != nsp:
                return incar, f'ERROR {t} has {len(expand_values(tags[t]))} values for the {nsp} species of the supercell'
        if 'LDAUU' in tags:
            notes.append(f"LDAUU {tags['LDAUU']} for {' '.join(sc['species'])}")
        return incar, ''.join('; ' + n for n in notes)

    def magnetic_xyz(self, incar, poscars):
        """04_supercells/snapshots_mag.xyz: one frame per supercell folder, Jmol's xyz format 'El x y z 0 mx my mz'
        with the vector = the initial moment of the INCAR (collinear: along z), for viewing the magnetic order
        (Jmol: vectors on). snapshots.xyz of SC-SNAPS keeps its velocities in these columns. '' without MAGMOM."""
        m = re.search(r'^\s*MAGMOM\s*=\s*([^!#\n]*)', incar, re.M | re.I)
        if not m:
            return ''
        vals = [float(v) for v in expand_values(m.group(1))]
        out = []
        for f in poscars:
            p = read_vasp_cell(f)
            n = len(p['X'])
            if len(vals) == 3 * n:
                mv = np.array(vals).reshape(n, 3)
            elif len(vals) == n:
                mv = np.zeros((n, 3)); mv[:, 2] = vals
            else:
                return f'; snapshots_mag.xyz not written (MAGMOM has {len(vals)} values for {n} atoms)'
            el = [p['species'][t] for t in p['types']]
            out.append(f'{n}\n {os.path.basename(os.path.dirname(f))}: initial moments (muB) as vectors'
                       + (' (collinear: along z)' if len(vals) == n else '') + '\n')
            out += [f'{e:3s} {x[0]:14.8f} {x[1]:14.8f} {x[2]:14.8f}  0.0  {v[0]:8.3f} {v[1]:8.3f} {v[2]:8.3f}\n'
                    for e, x, v in zip(el, p['X'], mv)]
        write(os.path.join(os.path.dirname(os.path.dirname(poscars[0])), 'snapshots_mag.xyz'), ''.join(out))
        return '; magnetic moments in 04_supercells/snapshots_mag.xyz (Jmol: vectors on)'

    def magmom_for_supercell(self, g, sc):
        spec = self.inst.get('magmom')
        if not spec or spec['kind'] != 'list':
            return
        p = read_vasp_cell(self.d('POSCAR0'))
        if len(spec['values']) != len(p['F']):
            return
        inv = np.linalg.inv(p['A'])
        vals = []
        for x in sc['X']:
            f = (x - p['F'] @ p['A']) @ inv
            k = int(np.argmin(np.abs(f - np.rint(f)).sum(1)))
            vals.append(spec['values'][k])
        g.instructions['magmom'] = {'kind': 'list', 'values': vals, 'raw': 'expanded to the supercell'}

    def t_vasp_snaps(self, retry=False):
        dirs, why = self.snap_dirs()
        if dirs is None:
            return False, why
        g = self.gen(os.path.join(dirs[0], 'POSCAR'), 'snaps')
        failed, ran = [], 0
        if self.fx['vasp_run'] == 'external':
            todo = [d for d in dirs if not outcar_status(d, 'static')[0] and
                    (retry or not os.path.isfile(os.path.join(d, 'OUTCAR')) or
                     outcar_status(d, 'static')[1].startswith('VASP did not finish'))]
            if todo:
                raise ExternalRun(todo)                 # all at once: independent runs
        for d in dirs:
            ok, msg = outcar_status(d, 'static')
            if ok:
                continue
            if os.path.isfile(os.path.join(d, 'OUTCAR')):
                if not retry and not msg.startswith('VASP did not finish'):   # finished, not converged
                    failed.append((d, msg))
                    continue
            dt = self.vasp(d, g)
            ran += 1
            ok, msg = outcar_status(d, 'static')
            self.clean(d)
            self.rep(f'  {os.path.basename(d)}: {msg} ({dt / 60:.1f} min)')
            if not ok:
                failed.append((d, msg))
        self.write_resubmit(failed)
        if failed:
            return False, f'{len(failed)} of {len(dirs)} failed:\n  ' + '\n  '.join(
                f'{os.path.basename(d)}: {m}' for d, m in failed) + '\n  fix INCAR_snaps if needed, then ./resubmit_failed.sh'
        return True, f'all {len(dirs)} converged ({ran} run now)'

    def write_workflow(self, poscar):
        """run_workflow.sh: the batch preamble of the profile (or none), then 'focex-agent.py run'."""
        g = self.gen(poscar, 'relax')
        note = ('# The #SBATCH and module load lines come from the profile (~/VASP-Flow/profiles/<name>.json: slurm, modules).\n'
                '# Module names on any cluster: module spider <name> (e.g. module spider vasp), and\n'
                '# module spider <name>/<version> for what has to be loaded before it.\n') if self.profile.get('slurm') else ''
        write(self.d('run_workflow.sh'), '#!/bin/bash\n' + note + g._run_sh_preamble('focex_agent') +
              f'cd "{self.proj}"\n"{sys.executable}" "{os.path.abspath(__file__)}" run "{self.proj}"\n', 0o755)

    def next_step(self):
        sub = 'sbatch' if self.profile.get('slurm') else 'bash'
        self.rep(f'next:  cd {self.proj} && {sub} run_workflow.sh')

    def write_resubmit(self, failed):
        g = self.gen(self.d('POSCAR.init'), 'relax')
        body = (g._run_sh_preamble('focex_retry') +
                f'cd "{self.proj}"\n# failed: {" ".join(os.path.basename(d) for d, _ in failed) or "none"}\n'
                f'{sys.executable} "{os.path.abspath(__file__)}" run "{self.proj}" --retry-failed\n')
        write(self.d('resubmit_failed.sh'), '#!/bin/bash\n' + body, 0o755)

    # ---- task 5: collect -----------------------------------------------------
    def t_collect(self):
        top = self.d('04_supercells')
        dirs, why = self.snap_dirs()
        if dirs is None:
            return False, why
        bad = [os.path.basename(d) for d in dirs if not outcar_status(d, 'static')[0]]
        if bad:
            return False, 'not converged: ' + ' '.join(bad)
        if not os.path.isfile(self.fx['read_outcar_x']):
            return False, f"read_outcar.x not found at {self.fx['read_outcar_x']} (set READ_OUTCAR_X)"
        out, worst_sum, worst_pos = [], 0.0, 0.0
        for d in dirs:
            subprocess.run(self.fx['read_outcar_x'], cwd=d, capture_output=True)
            txt = open(os.path.join(d, 'pos-forc.dat')).read()
            p = read_vasp_cell(os.path.join(d, 'POSCAR'))
            rows = [l.split() for l in txt.splitlines()[2:] if l.strip()]
            a = np.array(rows, dtype=float)
            if len(a) != len(p['X']):
                return False, f'{os.path.basename(d)}: {len(a)} atoms in pos-forc.dat, {len(p["X"])} in POSCAR'
            worst_sum = max(worst_sum, np.linalg.norm(a[:, 3:].sum(0)))
            df = (a[:, :3] - p['X']) @ np.linalg.inv(p['A'])
            worst_pos = max(worst_pos, np.abs((df - np.rint(df)) @ p['A']).max())
            out.append(txt)
        write(os.path.join(top, 'FORCEDISP1'), ''.join(out))
        shutil.copy(self.d('03_snaps', 'poscar_000'), os.path.join(top, 'POSCAR1'))
        msg = (f'FORCEDISP1 from {len(dirs)} OUTCARs; POSCAR1 = poscar_000; max |sum F| = {worst_sum:.2e} eV/A; '
               f'max |r_OUTCAR - r_POSCAR| = {worst_pos:.1e} A')
        sig, pstress = outcar_stress(os.path.join(top, 'snap_000', 'OUTCAR'))
        if sig is None:
            msg += '; no full stress tensor in snap_000/OUTCAR (ISIF < 2): stress.params not written'
        else:
            off = self.fx.get('stress_offset')
            if off:                       # printed - offset = the true printed stress; sigma = -(printed)/10 GPa
                sig = sig + np.diag(off) / 10.0
                msg += f'; STRESS_OFFSET {" ".join(f"{x:g}" for x in off)} kB removed'
            write(os.path.join(top, 'stress.params'),
                  ''.join(' '.join(f'{x:12.6f}' for x in r) + '\n' for r in sig) +
                  '# stress of the undistorted supercell (snap_000/OUTCAR), GPa, tensile > 0 (= -VASP in kB/10)\n')
            msg += f'; stress.params from snap_000: diagonal {np.diag(sig).round(3).tolist()} GPa'
            if pstress != 0:
                msg += f' -- WARNING: PSTRESS = {pstress:g} kB in snap_000, check stress.params'
        if worst_sum > 0.05 or worst_pos > 1e-3:
            msg += ' -- WARNING: check these'
        return True, msg

    # ---- task 6: FOCEX -------------------------------------------------------
    def structure_params(self, line5='1 1 1 1'):
        b = self.fx['blocks'].get('structure_params')
        if b:
            return b
        p, c = frame_of(self.d('POSCAR0'))
        n, fx = len(p['F']), self.fx
        masses = [p2f.MASS[s.split('_')[0]] for s in p['species']]
        L = [' '.join('%.10f' % x for x in c['latp']) + '   # a b c alpha beta gamma of the conventional cell',
             ' '.join('%.12g' % x for x in p2f.clean(c['prim_in_conv'], 12).ravel())
             + '   # primitive vectors in units of the conventional vectors',
             '1.0 0   # scale; 0 = atomic coordinates in units of the conventional vectors',
             ' '.join(str(x) for x in fx['ranks']) + '   # ranks 1..8 included',
             line5 + '   # translational, rotational, Huang invariances; enforce_inv (1 = exact, 3 = LASSO + refit)',
             '0 300   # itemp, tempk',
             f"{self.state.get('npairs', 1)} .false.   # number of POSCARi/FORCEDISPi pairs, verbosity",
             f"{len(p['species'])}   # number of atom types",
             ' '.join('%.4f' % m for m in masses) + '   # masses (amu)',
             ' '.join(p['species']) + '   # names',
             f'{n}   # atoms in the primitive cell',
             f"{fx['fc2range']}   # fc2range (0 = range from the supercell)"]
        for r, s in enumerate(fx['shells'], 2):
            L.append(' '.join([str(s)] * n) + f'   # shells for rank {r}')
        for k in range(n):
            L.append('%d %d %s' % (k + 1, p['types'][k] + 1, ' '.join('%.10f' % v for v in p2f.clean(c['Xc'][k]))))
        return '\n'.join(L) + '\n'

    def latdyn_params(self, wmax=None):
        b = self.fx['blocks'].get('latdyn_params')
        if b:
            return b
        L = open(os.path.join(HERE, 'latdyn.params')).read().splitlines()
        m = self.fx['focex_mesh'].split()
        if len(m) == 3:
            mesh = [int(x) for x in m]
        else:
            g = self.gen(self.d('POSCAR0'), 'scf')
            mesh = list(g._kpoints_from_kpra(g._kpra_of(m[0])))
        L[0] = f'{mesh[0]} {mesh[1]} {mesh[2]}        # kmesh ({m[0]})'
        if wmax:
            L[2] = f'{L[2].split()[0]} {wmax:.1f}    # number of frequency points for the DOS, omega_max (1/cm)'
        return '\n'.join(L) + '\n'

    def kpbs_params(self):
        b = self.fx['blocks'].get('kpbs_params')
        if b:
            return b
        if self.fx['kpath']:
            sys.exit('KPATH is non-zero: give the path in a kpbs_params ... kpbs_params_END block')
        g = self.gen(self.d('POSCAR0'), 'scf')
        sk = g._spglib_kpath()
        if sk is None:
            sys.exit('automatic k-path failed (spglib/pymatgen); give a kpbs_params block')
        segs, coords, name = sk
        pts = [segs[0][0]]
        for s in segs:
            for lbl in s:
                if lbl != pts[-1]:
                    pts.append(lbl)
        L = ['1   # primitive reduced coordinates', f"{self.fx['bs_nk']}   # k-points per direction",
             f'{len(pts) - 1}   # number of directions ({name}: {"|".join("-".join(s) for s in segs)})']
        L += ['%-3s %.6f %.6f %.6f' % ((lbl[:3],) + tuple(coords[lbl])) for lbl in pts]
        return '\n'.join(L) + '\n'

    def dielectric_params(self, flag):
        b = self.fx['blocks'].get('dielectric_params')
        if b:
            return b if flag else '0' + b[b.index('\n'):]
        n = sum(read_vasp_cell(self.d('POSCAR0'))['counts'])
        eps = np.array(self.state.get('eps') or np.eye(3))
        Z = self.state.get('Z') or [np.zeros((3, 3))] * n
        L = [f'{flag}   # born_flag (0: no Born charges / LO-TO)']
        L += [' '.join('%.6f' % x for x in r) + ('   # epsilon_inf' if i == 0 else '') for i, r in enumerate(eps)]
        for k, z in enumerate(Z):
            L += [' '.join('%.6f' % x for x in r) + (f'   # Born charge of atom {k + 1}' if i == 0 else '')
                  for i, r in enumerate(np.array(z))]
        return '\n'.join(L) + '\n'

    def run_focex(self, d, born, line5, label):
        os.makedirs(d, exist_ok=True)
        for i in range(1, self.state.get('npairs', 1) + 1):
            for f in (f'POSCAR{i}', f'FORCEDISP{i}'):
                shutil.copy(self.d('04_supercells', f), d)
        if os.path.isfile(self.d('04_supercells', 'stress.params')):
            shutil.copy(self.d('04_supercells', 'stress.params'), d)
        write(os.path.join(d, 'structure.params'), self.structure_params(line5))
        write(os.path.join(d, 'latdyn.params'), self.latdyn_params())
        write(os.path.join(d, 'kpbs.params'), self.kpbs_params())
        write(os.path.join(d, 'dielectric.params'), self.dielectric_params(born))
        write(os.path.join(d, 'default.params'),
              self.fx['blocks'].get('default_params') or open(os.path.join(HERE, 'default.params')).read())
        res = self.focex_once(d)
        if res['ok'] and res['wmax']:
            wdos = float(open(os.path.join(d, 'latdyn.params')).read().splitlines()[2].split()[1])
            if not (res['wmax'] <= wdos <= 2 * res['wmax']) and not self.fx['blocks'].get('latdyn_params'):
                w = math.ceil(1.1 * res['wmax'] / 10) * 10
                write(os.path.join(d, 'latdyn.params'), self.latdyn_params(w))
                self.rep(f'  {label}: DOS omega_max {wdos:g} -> {w} 1/cm (max frequency {res["wmax"]:.0f}), rerun')
                res = self.focex_once(d)
        self.rep(f'  {label}: ' + self.describe(res))
        return res

    def focex_once(self, d):
        for f in glob.glob(os.path.join(d, 'log*.dat')):
            os.remove(f)
        t0 = time.time()
        with open(os.path.join(d, 'focex.out'), 'w') as fh:
            rc = subprocess.run(self.fx['focex_x'], cwd=d, stdout=fh, stderr=subprocess.STDOUT).returncode
        res = dict(ok=False, minutes=(time.time() - t0) / 60, err=None, eq=None, unk=None, w2min=None, wmax=None)
        logs = sorted(glob.glob(os.path.join(d, 'log*.dat')), key=os.path.getmtime)
        if not logs:
            res['why'] = 'no log file; see focex.out: ' + ' / '.join(open(os.path.join(d, 'focex.out')).read().splitlines()[-3:])
            return res
        t = open(logs[-1], errors='replace').read()
        # enforce_inv 0 prints 'Percent error, || F_dft-F_fit || / || F_dft || =', 1-3 'Relative residual
        # ||F_dft-F_fit||/||F_dft|| =' (focex-agent 1.4 runs enforce_inv 1 and 3)
        m = re.findall(r'(?:Percent error, \|\| F_dft-F_fit \|\| / \|\| F_dft \|\||Relative residual \|\|F_dft-F_fit\|\|/\|\|F_dft\|\|)'
                       r'\s*=\s*([-\d.Ee+]+)', t)
        res['err'] = float(m[-1]) if m else None
        m = re.search(r'Force_constraints, # of independent FCs\s*=\s*(\d+)\s+(\d+)', t)
        if m:
            res['eq'], res['unk'] = int(m.group(1)), int(m.group(2))
        w2 = self.bands_omega2(d)
        if w2 is not None:
            res['w2min'] = float(w2.min() * OMEGA2_CM)
            res['wmax'] = float(math.sqrt(max(w2.max(), 0) * OMEGA2_CM))
        res['ok'] = rc == 0 and res['err'] is not None and w2 is not None
        if not res['ok']:
            # a run can crash after it wrote the fit error and bands.dat (e.g. in the elastic constants)
            res['why'] = (f'FOCEX stopped with exit code {rc}' if rc else 'FOCEX did not finish') + \
                '; see focex.out and ' + os.path.basename(logs[-1])
        return res

    @staticmethod
    def bands_omega2(d):
        """omega^2 (eV/(A^2 amu)) of every band at every point of bands.dat, or None."""
        f = os.path.join(d, 'bands.dat')
        if not os.path.isfile(f):
            return None
        nb = 3 * int(open(os.path.join(d, 'structure.params')).read().splitlines()[10].split()[0])
        toks = [float(x) for l in open(f).read().splitlines() if not l.lstrip().startswith('#') for x in l.split()]
        rec = 5 + 2 * nb
        if not toks or len(toks) % rec:
            return None
        a = np.array(toks).reshape(-1, rec)
        return a[:, 5:5 + nb]

    @staticmethod
    def describe(r):
        if not r['ok']:
            return 'FAILED: ' + r.get('why', '')
        s = (f"{r['err']:.3f} % error; {r['eq']} equations / {r['unk']} unknowns"
             + (f" = {r['eq'] / r['unk']:.2f}" if r['eq'] and r['unk'] else '')
             + f"; min omega^2 = {r['w2min']:.1f} 1/cm^2, max omega = {r['wmax']:.0f} 1/cm; {r['minutes']:.1f} min")
        if r['w2min'] < -25:
            s += ' -- WARNING: imaginary modes (omega^2 < -25 1/cm^2)'
        if r['wmax'] > 2500:
            s += ' -- WARNING: highest frequency above 2500 1/cm'
        return s

    def t_focex(self):
        top = self.d('05_focex')
        if not os.path.isfile(self.fx['focex_x']):
            return False, f"FOCEX binary not found at {self.fx['focex_x']} (set FOCEX_X)"
        born = 1 if self.state.get('Z') or self.fx['blocks'].get('dielectric_params') else 0
        runs = {}
        if born:
            runs['NOBORN'] = self.run_focex(os.path.join(top, 'NOBORN'), 0, '1 1 1 1', 'NOBORN (born_flag=0)')
        # enforce_inv = 1 (invariances exact, error bars) for the main fit; LASSO with the refit (3)
        # (focex-agent 1.4, 2026-11-07, user: was 0 and 2; see the FOCEX FAQ)
        runs['main'] = main = self.run_focex(top, born, '1 1 1 1', f'main (born_flag={born})')
        if main['ok'] and main['eq'] and main['unk'] and main['eq'] / main['unk'] < self.fx['lasso_ratio']:
            self.rep(f"  equations/unknowns = {main['eq'] / main['unk']:.2f} < {self.fx['lasso_ratio']}: LASSO run")
            runs['LASSO'] = self.run_focex(os.path.join(top, 'LASSO'), born, '1 1 1 3', 'LASSO (enforce_inv=3)')
        self.state['focex'] = runs
        if 'LASSO' in runs:
            self.rep('  comparison (positive = min omega^2 > -10 1/cm^2):')
            for k in ('main', 'LASSO'):
                r = runs[k]
                if r['ok']:
                    self.rep(f"    {k:6s} error {r['err']:8.3f} %   min omega^2 {r['w2min']:10.1f}   "
                             f"{'positive' if r['w2min'] > -10 else 'NOT positive'}")
        return main['ok'], self.describe(main)

    # ---- task 7: plots -------------------------------------------------------
    def t_plots(self):
        top = self.d('05_focex')
        msgs, ok = [], True
        for d in [top] + [os.path.join(top, s) for s in ('NOBORN', 'LASSO') if os.path.isdir(os.path.join(top, s))]:
            good, m = self.plots_in(d)
            ok &= good
            msgs.append(f'{os.path.relpath(d, self.proj)}: {m}')
        return ok, '\n  '.join(msgs)

    def plots_in(self, d):
        """makeplots.sh = PLOTTING/allplots on d (2026-11-02, user: the same script outside the agent):
        fills the .plt templates (KTICS.BS, capped bounds, NAME = the folder above 05_focex) and makes
        the PDFs"""
        if not os.path.isfile(os.path.join(d, 'KTICS.BS')):
            return False, 'no KTICS.BS (FOCEX did not finish)'
        allplots = os.path.join(HERE, 'PLOTTING', 'allplots')
        write(os.path.join(d, 'makeplots.sh'),
              f'#!/bin/bash\n# remakes the plots of this FOCEX folder\ncd "$(dirname "$0")"\n'
              f'exec python3 "{allplots}" .\n', 0o755)
        if gnuplot_pdf() is None:
            subprocess.run([sys.executable, allplots, '.', '--no-run'], cwd=d, capture_output=True)
            return True, ('scripts written; no gnuplot with the pdfcairo terminal here, run makeplots.sh where '
                          'there is one')
        r = subprocess.run(['bash', 'makeplots.sh'], cwd=d, capture_output=True, text=True)
        lines = [l.strip() for l in (r.stdout + r.stderr).splitlines()[1:] if l.strip()]
        return r.returncode == 0, '; '.join(lines)

    # ---- task 8: THERMACOND ----------------------------------------------------
    def tc_latdyn(self, base=None):
        """latdyn.params of 07_thermacond: a thermacond_params block as it is, else that of 05_focex with
        lines 1, 6, 8, 9 from THERMACOND_*"""
        fx, job = self.fx, self.fx['tc_job']
        if fx['blocks'].get('thermacond_params'):
            return fx['blocks']['thermacond_params'].rstrip('\n').splitlines()
        L = (base or open(self.d('05_focex', 'latdyn.params')).read()).splitlines()
        if fx['tc_mesh']:
            L[0] = ' '.join(fx['tc_mesh']) + '        # 1 kmesh (THERMACOND_MESH)'
        if fx['tc_T']:
            L[5] = ' '.join(fx['tc_T']) + '        # 6 Tmin(K), Tmax(K), number of temperatures (THERMACOND_T)'
        L[7] = f'{job}            # 8 THERMACOND job ({TC_JOB[job]})'
        if fx['tc_solver'] is not None:
            L[8] = f"{fx['tc_solver']}            # 9 THERMACOND direct solver: 0 conjugate gradient, 1 SVD"
        return L

    def t_thermacond(self):
        """THERMACOND on the force constants of 05_focex, in 07_thermacond, when THERMACOND_JOB
        (latdyn.params line 8) is 1, 2 or 3."""
        job, fx = self.fx['tc_job'], self.fx
        if job == 0:
            return True, ('skipped: THERMACOND_JOB 0' + ('' if fx['tc_job_given'] else ' (not given)')
                          + ' -- set THERMACOND_JOB: 2 or 3 for the thermal conductivity')
        src, d = self.d('05_focex'), self.tc_dir()
        fc3 = os.path.join(src, 'fc3.dat')
        if not os.path.isfile(fc3) or not any(l.split() and not l.lstrip().startswith('#')
                                               for l in open(fc3).read().splitlines()[1:]):
            return False, 'THERMACOND needs cubic force constants: fit rank 3 (RANKS: 1 1 1 0 0 0 0 0)'
        if not os.path.isfile(fx['thermacond_x']):
            return False, (f"THERMACOND binary not found at {fx['thermacond_x']} (make thermacond, or "
                           'set THERMACOND_X)')
        os.makedirs(d, exist_ok=True)
        # the files are copied over those of an earlier run; V3DELTA/ stays: THERMACOND reuses it when
        # lat_fc/fc2/fc3/structure/dielectric and the mesh are those it was made with
        # copied only when changed, so that their time stamps say whether thermacond.out is newer
        ins = ('lat_fc.dat', 'fc2.dat', 'fc3.dat', 'structure.params', 'dielectric.params', 'kpbs.params')
        for f in ins:
            if os.path.isfile(os.path.join(src, f)):
                copy_if_changed(os.path.join(src, f), d)
        L = self.tc_latdyn()
        write_if_changed(os.path.join(d, 'latdyn.params'), '\n'.join(L) + '\n')
        t0 = time.time()
        if fx['tc_run'] == 'external':
            out = os.path.join(d, 'thermacond.out')
            newest = max(os.path.getmtime(os.path.join(d, f)) for f in ins + ('latdyn.params',)
                         if os.path.isfile(os.path.join(d, f)))
            if not (os.path.isfile(out) and os.path.getmtime(out) > newest):
                raise ExternalRun([d], kind='thermacond')
            rc = 0                                      # run elsewhere; its outputs are judged below
        else:
            if is_serial_tcond(fx['thermacond_x']):
                # serial build: one process on threads (MPI ranks of the job/MPI line as the thread count);
                # srun/mpirun would start that many independent copies
                n = int(self.inst.get('mpi_np') or 1)
                cmd = (f'TCOND_NUM_THREADS={n} ' if n > 1 else '') + f'"{fx["thermacond_x"]}"'
            else:
                cmd = self.gen(self.d('POSCAR0'), 'scf')._get_mpi_cmd(fx['thermacond_x'])
                # configure's vasp-env (mpi_cmd 'vasp-env srun') swaps the job's modules for the VASP module;
                # THERMACOND is built with the job's MPI, and under vasp-env srun every rank started as a
                # singleton (40 copies of the whole job; Rivanna 2026-11-09), so it is started without it
                cmd = re.sub(r'^\s*("?)\S*/vasp-env\1\s+', '', cmd)
                # OpenMPI 5 under srun: PMIx's shared-memory store failed before MPI_Init on Rivanna
                # (PMIX_ERR_FILE_OPEN_FAILURE in gds_shmem2.c, exit code 14; plain srun and --mpi=pmix alike);
                # the hash store works (tested with 4 and 40 ranks, 2026-11-09). A PMIX_MCA_gds already set wins.
                if cmd.lstrip().startswith('srun'):
                    cmd = 'PMIX_MCA_gds=${PMIX_MCA_gds:-hash} ' + cmd
            rc = subprocess.run(f'{cmd} > thermacond.out 2>&1', shell=True, cwd=d).returncode
        head = (f"job {job} ({TC_JOB[job]}), mesh {L[0].split('#')[0].strip()}, T {L[5].split('#')[0].strip()}, "
                f'{(time.time() - t0) / 60:.1f} min')
        if job == 1:
            ok = rc == 0 and os.path.isfile(os.path.join(d, 'lw-selfenergy.dat'))
            return ok, head + ('; linewidths in lw-selfenergy.dat' if ok else self.tc_failure(d, rc))
        rows = self.kappa_rows(os.path.join(d, 'kappa_vs_T.dat'))
        if rc != 0 or not rows:
            return False, head + self.tc_failure(d, rc)
        kc = {round(r[0], 2): r for r in self.kappa_rows(os.path.join(d, 'kappa_wigner.dat'))}
        full = 'direct' if job == 2 else 'iterative'
        msg = [head, f'T(K)   kappa_RTA   kappa_{full}   kappa_C   (W/mK, average of xx yy zz)']
        for r in [rows[i] for i in sorted({round(k * (len(rows) - 1) / 6) for k in range(7)})]:   # <= 7, ends included
            c = kc.get(round(r[0], 2))
            msg.append(f'{r[0]:7.1f} {sum(r[1:4]) / 3:10.3f} {sum(r[7:10]) / 3:12.3f} '
                       + (f'{sum(c[13:16]) / 3:9.4f}' if c and len(c) > 15 else '      n/a'))
        low = [r[0] for r in rows if sum(r[7:10]) < 0.99 * sum(r[1:4])]
        if low:
            msg.append(f'WARNING: the {full} solution is below the RTA at {len(low)} temperature(s); the full '
                       'solution should not be lower -- the direct solver has an open issue (compare with job 3)')
        msg.append(self.tc_plots(d))
        return True, '\n  '.join(msg)

    def tc_plots(self, d):
        """THERMACOND 2.2 plots in d: kappa(T), Wigner by bands, cumulative kappa versus mean free path and frequency
        per temperature -- gnuplot (allplots, kappa.plt) when it is there, else matplotlib (kappaplots)."""
        if not os.path.isfile(os.path.join(d, 'kappa_summary.dat')):
            return 'plots: no kappa_summary.dat (THERMACOND before 2.2)'
        name = os.path.basename(self.proj)
        if gnuplot_pdf():
            r = subprocess.run([sys.executable, os.path.join(HERE, 'PLOTTING', 'allplots'), '.', '--name', name],
                               cwd=d, capture_output=True, text=True)
            ok = [l.strip() for l in r.stdout.splitlines() if 'kappa.plt' in l]
            return 'plots (gnuplot): ' + (ok[0] if ok else (r.stdout + r.stderr)[-200:])
        r = subprocess.run([sys.executable, os.path.join(HERE, 'PLOTTING', 'kappaplots'), '.', '--name', name],
                           cwd=d, capture_output=True, text=True)
        return 'plots (matplotlib): ' + (r.stdout + r.stderr).strip()[-300:]

    @staticmethod
    def kappa_rows(path):
        rows = []
        if os.path.isfile(path):
            for l in open(path).read().splitlines():
                if l.strip() and not l.lstrip().startswith('#'):
                    try:
                        rows.append([float(x) for x in l.split()])
                    except ValueError:
                        pass
        return rows

    @staticmethod
    def tc_failure(d, rc):
        tail = [l.strip() for l in open(os.path.join(d, 'thermacond.out'), errors='replace').read().splitlines()
                if l.strip()][-3:] if os.path.isfile(os.path.join(d, 'thermacond.out')) else []
        return f'; THERMACOND failed (exit code {rc}); thermacond.out: ' + ' / '.join(tail)

    # ---- driver --------------------------------------------------------------
    def output_gone(self, t):
        """Why the output of a task marked done is no longer there or out of date ('' if it is fine).
        Checked for the supercell tasks, whose folders users remove or set aside, and for thermacond,
        whose settings users change; tasks provided by DATA are not checked."""
        if str(self.state.get('done', {}).get(t, '')).startswith('provided'):
            return ''
        if t in ('supercells', 'vasp_snaps', 'collect'):
            dirs, why = self.snap_dirs()
            if dirs is None:
                return why
            if t == 'vasp_snaps':
                bad = [os.path.basename(d) for d in dirs if not outcar_status(d, 'static')[0]]
                if bad:
                    return f'no converged OUTCAR in {len(bad)} folder(s): {" ".join(bad[:10])}'
            if t == 'collect':
                f = self.d('04_supercells', 'FORCEDISP1')
                n = open(f).read().count('POSITION') if os.path.isfile(f) else 0
                if n != len(dirs):
                    return f'FORCEDISP1 has {n} snapshots, 04_supercells {len(dirs)}'
        if t == 'thermacond' and self.fx['tc_job'] != 0:
            f = os.path.join(self.tc_dir(), 'latdyn.params')
            if not os.path.isfile(f):
                return f'THERMACOND_JOB is {self.fx["tc_job"]} and THERMACOND has not run'
            if os.path.isfile(self.d('05_focex', 'latdyn.params')) and \
                    open(f).read().splitlines() != self.tc_latdyn():
                return 'the THERMACOND settings (THERMACOND_JOB, _MESH, _T, _SOLVER) changed'
        return ''

    def run(self, only=None, retry=False):
        todo = [only] if only else TASKS
        for t in todo:
            done = self.state.setdefault('done', {})
            if t == 'snaps' and t in done and not only and not str(done[t]).startswith('provided'):
                ok, msg = self.add_snapshots()
                if msg:
                    self.rep(f'=== snaps -- {msg}')
                if not ok:
                    self.state['failed'] = t
                    self.save()
                    return False
            if t in done and not only:
                why = self.output_gone(t)
                if why:
                    # this task and every later one depend on the output that changed: do them again
                    redo = TASKS[TASKS.index(t):]
                    self.rep(f'=== {t} -- marked done, but {why}: redoing {", ".join(x for x in redo if x in done)}')
                    for x in redo:
                        done.pop(x, None)
                    self.save()
            if t in done and not only and not (retry and t in ('vasp_snaps',)):
                continue
            self.rep(f'=== {t} -- started')
            t0 = time.time()
            fn = getattr(self, 't_' + t)
            try:
                ok, msg = fn(retry) if t == 'vasp_snaps' else fn()
            except ExternalRun as x:
                rel = [os.path.relpath(d, self.proj) for d in x.dirs]
                self.state['waiting'] = {'task': t, 'kind': x.kind, 'dirs': rel, 'since': now(),
                                         'after': {os.path.relpath(a, self.proj): os.path.relpath(b, self.proj)
                                                   for a, b in x.after.items()}}
                self.save()
                self.rep(f'=== {t} -- WAITING: {"VASP" if x.kind == "vasp" else "THERMACOND"} inputs ready in '
                         f'{len(rel)} folder(s): {" ".join(rel[:4])}{" ..." if len(rel) > 4 else ""}; run them '
                         'outside the agent (aladyn-gui: bundle or submit), then continue')
                return 'waiting'
            self.state.pop('waiting', None)
            self.rep(f'=== {t} -- {"OK" if ok else "FAILED"} ({(time.time() - t0) / 60:.1f} min): {msg}')
            if not ok:
                self.state['failed'] = t
                self.save()
                self.rep(f'stopped at task {t}; to find out why and get a fix script:  '
                         f'"{sys.executable}" "{os.path.join(AGENT, "diagnose.py")}" "{self.proj}"')
                return False
            done[t] = now()
            self.state.pop('failed', None)
            self.save()
            if t == self.fx.get('stop_after') and not only:
                self.rep(f'=== stopped after {t} (STOP_AFTER: {t}); to go on, delete that line in '
                         f'{os.path.join(self.proj, INSTR)} and submit run_workflow.sh again')
                return True
        return True


# ── init ──────────────────────────────────────────────────────────────────────
def init(args):
    inst, fx = parse_instructions(args.instr)
    proj = fx['project_dir']
    if os.path.exists(os.path.join(proj, 'state.json')) and not args.force:
        sys.exit(f'{proj} exists; use --force to start it again')
    if fx['mode'] == 'from_scratch' and not os.path.isfile(args.coord):
        sys.exit(f'coordinate file {args.coord} not found (default name POSCAR; another name: --coord FILE)')
    args.profile = args.profile or fx['profile'] or os.environ.get('FOCEX_PROFILE') or None
    prof = vagent.load_profile(args.profile or '')
    if prof.get('slurm'):
        try:
            slurm_account(inst, prof)            # a SLURM run needs an account (profile or ACCOUNT:)
        except InstructionError as e:
            sys.exit(f'ERROR: {e}')
    os.makedirs(proj, exist_ok=True)
    shutil.copy(args.instr, os.path.join(proj, INSTR))
    if fx['mode'] == 'extract_only':
        return init_data(args, fx, proj)
    shutil.copy(args.coord, os.path.join(proj, 'POSCAR.input'))
    p = read_vasp_cell(args.coord)
    p, cell_note = primitive_input(p, input_moments(inst, fx, len(p['F'])))
    c = canonical_cell(p)
    write_poscar(os.path.join(proj, 'POSCAR.init'), f"{fx['project']} (Niggli-reduced, FOCEX frame)",
                 c['V'], p['species'], p['counts'], c['F'])
    state = {'profile': args.profile or fx['profile'] or '', 'created': now(), 'mode': 'from_scratch'}
    write(os.path.join(proj, 'state.json'), json.dumps(state, indent=1))
    metal = all(s.split('_')[0] in METALS for s in p['species'])
    # relax and SCF share the mesh: coarse (1000 kpra), fine for metals/intermetallics
    dens = 'fine' if metal else 'coarse'
    with open(os.path.join(proj, INSTR), 'a') as fh:
        if not re.search(r'^\s*(RELAX_)?KMESH(_DENSITY)?\s*:', open(args.instr).read(), re.M | re.I):
            fh.write(f'\n# added by focex-agent init: relax and SCF k-mesh ({"metal" if metal else "not a metal"})\n'
                     f'RELAX_KMESH_DENSITY: {dens}\nKMESH_DENSITY: {dens}\n')
    a = Agent(proj)
    rep = a.rep
    rep(f"focex-agent init: {fx['project']}, {len(p['F'])} atoms, space group {c['sg']}, accuracy "
        f"{fx['accuracy']}, relax/SCF k-mesh {dens}{cell_note}")
    if a.mpi_note:
        rep('  ' + a.mpi_note)
    potcar_dir, _ = potcar_dir_for(inst.get('functional', 'R2SCAN'))
    if potcar_dir and os.path.isdir(potcar_dir):
        vagent.build_potcar(p['species'], potcar_dir, os.path.join(proj, 'POTCAR'))
    else:
        rep('WARNING: no POTCAR library (VASP_POTCAR_DIR); put a POTCAR in the project folder before running')
    for x in ('focex_x', 'sc_snaps_x', 'read_outcar_x') + (('thermacond_x',) if fx['tc_job'] else ()):
        if not os.path.isfile(fx[x]):
            rep(f'WARNING: {x.upper()} = {fx[x]} not found')
    # r2SCAN (the default functional) needs VASP 6; VASP 5.4.4 ignores METAGGA = R2SCAN without a word and
    # runs PBE. VASP_VERSION comes from ./configure (~/.aladyn.env); the relax task checks the OUTCAR as well.
    fun = inst.get('functional', 'R2SCAN')
    vv = os.environ.get('VASP_VERSION', '')
    if needs_vasp6(fun) and re.match(r'\d+', vv) and int(re.match(r'\d+', vv).group()) < 6:
        rep(f'NOTE: FUNCTIONAL {fun} with VASP {vv} (VASP_VERSION, from ./configure): VASP 5.4.4 has r2SCAN only with '
            'the r2SCAN patch (metagga544.diff of gitlab.com/dhamil/r2scan-subroutines); without it VASP runs PBE '
            'without saying so. The relaxation stops at its first OUTCAR if METAGGA was ignored')
    a.write_relax()
    a.write_workflow(os.path.join(proj, 'POSCAR.init'))
    a.write_resubmit([])
    a.rep(f'wrote {proj}: POSCAR.init, 01_relax/, run_workflow.sh')
    a.next_step()


def init_data(args, fx, proj):
    """FOCEX on existing POSCARi/FORCEDISPi (DATA: folder): rotate into the FOCEX frame, write
    POSCAR0 and 04_supercells/, mark the DFT tasks as provided."""
    src = fx['data']
    n = 0
    while os.path.isfile(os.path.join(src, f'POSCAR{n + 1}')) and os.path.isfile(os.path.join(src, f'FORCEDISP{n + 1}')):
        n += 1
    if n == 0:
        sys.exit(f'DATA: no POSCAR1 + FORCEDISP1 in {src}')
    sc = read_poscar_any(os.path.join(src, 'POSCAR1'), fx['species'])
    prim_file = os.path.join(src, 'POSCAR')
    if os.path.isfile(prim_file):
        p = read_poscar_any(prim_file, fx['species'] or sc['species'])
        how = 'POSCAR of the DATA folder'
    else:
        import spglib
        # no_idealize: keep the orientation of the data (find_primitive would rotate the cell)
        V, F, t = spglib.standardize_cell((sc['A'], sc['F'], sc['types']), to_primitive=True,
                                          no_idealize=True, symprec=1e-3)
        o = np.argsort(t, kind='stable')
        V, F, t = np.array(V), np.array(F)[o], np.array(t)[o]
        counts = [int((t == k).sum()) for k in range(len(sc['species']))]
        p = {'A': V, 'F': F, 'types': list(t), 'species': sc['species'], 'counts': counts}
        how = 'primitive cell of POSCAR1 (spglib)'
    V0, R, dev = ideal_frame(p['A'], p['F'], p['types'])
    if dev > 1e-3:
        sys.exit(f'DATA: the primitive cell is {dev:.1e} A from its symmetrized lattice; check the data')
    write_poscar(os.path.join(proj, 'POSCAR0'), f"{fx['project']} ({how}, FOCEX frame)",
                 V0, p['species'], p['counts'], p['F'] % 1.0)
    p0, c = frame_of(os.path.join(proj, 'POSCAR0'))
    top = os.path.join(proj, '04_supercells')
    os.makedirs(top, exist_ok=True)
    for i in range(1, n + 1):
        q = read_poscar_any(os.path.join(src, f'POSCAR{i}'), fx['species'] or sc['species'])
        write_poscar(os.path.join(top, f'POSCAR{i}'), f'POSCAR{i} of {src}, rotated into the FOCEX frame',
                     q['A'] @ R, q['species'], q['counts'], q['F'])
        rotate_forcedisp(os.path.join(src, f'FORCEDISP{i}'), os.path.join(top, f'FORCEDISP{i}'), R)
    state = {'profile': args.profile or fx['profile'] or '', 'created': now(), 'mode': 'extract_only', 'npairs': n,
             'data': os.path.abspath(src), 'rotation': R.tolist(), 'sg': c['sg']}
    diel = os.path.join(src, 'dielectric.params')
    if os.path.isfile(diel) and not fx['blocks'].get('dielectric_params'):
        rows = []
        for l in open(diel).read().splitlines()[1:]:          # line 1: born_flag
            try:
                rows.append([float(x) for x in l.split()[:3]])
            except ValueError:
                continue
        rows = [r for r in rows if len(r) == 3]
        na = len(p['F'])
        if not os.path.isfile(prim_file) and max(p['counts']) > 1:
            sys.exit(f'DATA: {diel} lists Born charges per atom of a primitive cell: put that cell, as POSCAR '
                     f'in the frame of the data, in {src} (the order of the atoms is not defined otherwise)')
        if len(rows) < 3 + 3 * na:
            sys.exit(f'DATA: {diel} has {len(rows)} rows of 3 numbers, {3 + 3 * na} needed')
        eps = np.array(rows[:3])
        Z = [np.array(rows[3 + 3 * k:6 + 3 * k]) for k in range(na)]
        if any(np.abs(z).max() > 0 for z in Z):
            state['eps'] = (R.T @ eps @ R).tolist()
            state['Z'] = [(R.T @ z @ R).tolist() for z in Z]
    state['done'] = {t: 'provided (DATA)' for t in TASKS[:TASKS.index('focex')]}
    write(os.path.join(proj, 'state.json'), json.dumps(state, indent=1))
    a = Agent(proj)
    a.rep(f"focex-agent init (DATA): {fx['project']}, {n} POSCAR/FORCEDISP pair(s) from {src}; {how}: "
          f"{len(p['F'])} atoms, space group {c['sg']}; rotation into the FOCEX frame "
          + ('none' if np.allclose(R, np.eye(3), atol=1e-9) else np.array2string(R, precision=6).replace('\n', ''))
          + f'; lattice symmetrized by {dev:.1e} A'
          + ('; Born charges and epsilon from dielectric.params' if state.get('Z') else ''))
    a.write_workflow(os.path.join(proj, 'POSCAR0'))
    a.rep(f'wrote {proj}: POSCAR0, 04_supercells/ (the rotated data), run_workflow.sh (tasks focex, plots)')
    a.next_step()


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd in ('run', 'plots'):
        ap = argparse.ArgumentParser(prog=f'focex-agent.py {cmd}')
        ap.add_argument('project', help='the project folder (PROJECT_DIR)')
        if cmd == 'run':
            ap.add_argument('--only', choices=TASKS, help='do this one task again')
            ap.add_argument('--retry-failed', action='store_true', help='rerun the failed supercell VASP runs')
        args = ap.parse_args(sys.argv[2:])
        if cmd == 'run':
            r = Agent(args.project).run(args.only, args.retry_failed)
            sys.exit(0 if r is True else 3 if r == 'waiting' else 1)      # 3: waiting for an external run
        sys.exit(0 if Agent(args.project).run('plots') else 1)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--instr', default=INSTR, help=f'instruction file (default {INSTR})')
    ap.add_argument('--coord', default='POSCAR', help='crystal structure, VASP format (default POSCAR; '
                    'MODE: from_scratch only)')
    ap.add_argument('--profile', default=None, help='VASP-Flow profile: workstation, slurm, or the '
                    'name of your own ~/VASP-Flow/profiles/<name>.json (default: PROFILE in the instruction '
                    'file, else no batch system)')
    ap.add_argument('--force', action='store_true', help='set up an existing project folder again')
    ap.add_argument('--version', action='version', version=f'focex-agent {__version__}')
    ap.add_argument('-i', '--instructions', dest='instr', help=argparse.SUPPRESS)   # names before 2026-10-23
    ap.add_argument('-s', '--structure', dest='coord', help=argparse.SUPPRESS)
    argv = sys.argv[2:] if cmd in ('init', 'setup') else sys.argv[1:]
    init(ap.parse_args(argv))


if __name__ == '__main__':
    main()
