#!/usr/bin/env python3
"""diagnose.py 1.0 -- where and why a focex-agent project stopped, and a script that fixes it.

    diagnose.py <PROJECT_DIR> [<PROJECT_DIR> ...]

Diagnosis is read-only.  For every task of the workflow (relax ... plots) it decides done /
FAILED / not reached / running from state.json, report.log and the outputs (OUTCAR, vasp.out,
focex.out, the SLURM logs), finds the cause and writes into <PROJECT_DIR>/diagnosis/:

  report.txt   what failed, where, why (the evidence), what fix_it.sh will change
  fix_it.sh    checks that nothing runs and nothing changed since the diagnosis, backs up the
               files it edits (diagnosis/backup_<date>/), applies the fixes, resubmits
               (sbatch if run_workflow.sh has #SBATCH lines, else bash in the background).
               When the cause is unknown or the fix needs your judgement it is HELD: it only
               prints what to do; bash fix_it.sh --anyway runs its suggested commands.

VASP failures are classified with VASP-Flow's vasp_diagnose.py (same causes and messages):
  timeout     time limit                     walltime x2 in run_workflow.sh, resubmit
  oom         out of memory                  KPAR = 1 (held if it is 1 already)
  scf_nelm    SCF did not converge in NELM   ALGO = All, NELM = 200 (held if ALGO = All already)
  edddav      Davidson eigensolver failed    ALGO = All (held if already)
  zbrent      relaxation line search failed  IBRION = 1 (relaxation)
  relax_unconv relaxation not converged      resubmit (the agent continues from CONTCAR)
  mpi_init    MPI could not start            resubmit unchanged
  not_run     never ran                      resubmit
  others (quota, symmetry, wavecar, chgcar, unknown): held, with the error lines
New INCAR tags go into the existing INCARs of the failed runs and into the matching block of
the project's instruction file (INCAR_relax, INCAR_scf, INCAR_dfpt, INCAR_snaps), so that later
runs use them too.  FOCEX, SC-SNAPS and read_outcar.x stops are recognised by their messages.
"""
import argparse
import datetime
import glob
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
__version__ = '1.0'
AGENT = os.path.dirname(os.path.realpath(__file__))
VASPFLOW = os.path.expanduser(os.environ.get('VASP_FLOW_DIR', '~/VASP-Flow'))
TASKS = ('relax', 'poscar0', 'scf', 'dfpt', 'snaps', 'supercells', 'vasp_snaps', 'collect',
         'focex', 'plots')
INSTR, OLD_INSTR = 'focex-agent-instructions.txt', 'instructions.txt'
VASP_DIR = {'relax': '01_relax', 'scf': '02_scf', 'dfpt': '02_scf/dfpt'}
BLOCK = {'relax': 'INCAR_relax', 'scf': 'INCAR_scf', 'dfpt': 'INCAR_dfpt', 'vasp_snaps': 'INCAR_snaps'}


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


try:
    vd = _load(os.path.join(VASPFLOW, 'vasp_diagnose.py'), 'vasp_diagnose')
except Exception as e:                                         # pragma: no cover
    sys.exit(f'diagnose.py needs VASP-Flow at {VASPFLOW} (set VASP_FLOW_DIR): {e}')

# FOCEX / SC-SNAPS / read_outcar.x messages: (cause, regex, what it means, what to do)
FOCEX_PAT = [
    ('ws_outside', r'R\+tp-t not inside|R\+taup-tau not inside',
     'an FC2 pair of a kept symmetry group lies outside the Wigner-Seitz cell of the supercell: '
     'the supercell does not have the symmetry of the lattice (e.g. 4x3x2 of a hexagonal cell)',
     'use a supercell with the symmetry of the lattice (SUPERCELL: SMALL/LARGE pick one since '
     '2026-10-23, or give 9 integers), then redo from the snapshots'),
    ('sc_map', r'SC ATOM can not be mapped|NO SHIFT VECTOR COULD BE DEFINED|is not integer',
     'POSCAR1 does not fit the primitive cell of structure.params (frame, lattice or positions)',
     'check that POSCAR1 is a supercell of POSCAR0 in the same Cartesian frame'),
    ('sizing', r'sizing disagreement|size of fc2 exceeded',
     'an array-size inconsistency in the force-constant set-up (a FOCEX bug)',
     'report it with the folder; a different RANGE may avoid it'),
    ('forcedisp', r'ERROR in reading the force file|READ_SUPERCELL: could not read',
     'FORCEDISP1 or POSCAR1 could not be read', 'check the files in 04_supercells'),
    ('read_fc', r'READ_FCS_3: the file exists but is empty', 'an fc*.dat to be read is empty',
     'remove it or set include_fc'),
    ('volume', r'[Vv]olume is zero', 'a cell with zero volume (structure.params lines 1-2)',
     'check structure.params'),
]


def read(path, tail=None):
    return vd.read(path, tail=tail)


def now():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


class Finding:
    """One problem: task, where, cause, what it means, evidence, error lines, fix."""
    def __init__(self, task, where, cause, meaning, evidence='', errs=(), fix=(), held=False, todo=''):
        self.task, self.where, self.cause, self.meaning = task, where, cause, meaning
        self.evidence, self.errs, self.fix, self.held, self.todo = evidence, list(errs), list(fix), held, todo


# ── reading a project ─────────────────────────────────────────────────────────
class Project:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        st = os.path.join(self.path, 'state.json')
        if not os.path.isfile(st):
            sys.exit(f'{self.path}: no state.json -- not a focex-agent project (set up with focex-agent.py)')
        self.state = json.load(open(st))
        self.log = read(self.f('report.log'))
        self.instr = self.f(INSTR) if os.path.isfile(self.f(INSTR)) else self.f(OLD_INSTR)
        self.run_sh = read(self.f('run_workflow.sh'))
        self.slurm = '#SBATCH' in self.run_sh

    def f(self, *a):
        return os.path.join(self.path, *a)

    def failed_msg(self, task):
        """The last '=== task -- FAILED' entry of report.log with its continuation lines."""
        m = list(re.finditer(r'=== ' + task + r' -- FAILED \([^)]*\): (.*?)(?=\n\[|\Z)', self.log, re.S))
        return m[-1].group(1).strip() if m else ''

    def running(self):
        """Why the project seems to be running ('' if not): a SLURM job in this folder, or a
        local 'focex-agent.py run' of it."""
        if shutil.which('squeue'):
            out = vd.run(['squeue', '-h', '-u', os.environ.get('USER', ''), '-o', '%i %T %Z'])
            for line in out.splitlines():
                p = line.split()
                if len(p) >= 3 and os.path.realpath(p[2]).startswith(os.path.realpath(self.path)):
                    return f'SLURM job {p[0]} {p[1]} in {p[2]}'
        out = vd.run(['ps', '-axo', 'pid,command'])
        for line in out.splitlines():
            if 'focex-agent.py' in line and ' run ' in line and self.path in line:
                return f'process {line.split()[0]}: {line.split(None, 1)[1][:120]}'
        return ''

    def job_log(self):
        """Text of the newest SLURM output/error files of the project (tails)."""
        logs = sorted(glob.glob(self.f('slurm-*.out')) + glob.glob(self.f('slurm-*.err')) +
                      glob.glob(self.f('*.o[0-9]*')) + glob.glob(self.f('run_workflow.log')), key=os.path.getmtime)
        return '\n'.join(read(x, tail=200000) for x in logs[-2:])


# ── VASP runs ─────────────────────────────────────────────────────────────────
def vasp_cause(d, kind, joblog=''):
    """(cause, evidence, error lines) of the VASP run in folder d, or ('ok', ...).
    kind: relax | static | dfpt."""
    outcar = os.path.join(d, 'OUTCAR')
    text = read(os.path.join(d, 'vasp.out'), tail=200000) + '\n' + joblog + '\n' + read(outcar, tail=20000)
    errs = [e for e in vd.error_lines(text, 12) if 'EDIFF is reached' not in e][:6]
    if not os.path.isfile(outcar):
        if not read(os.path.join(d, 'vasp.out')).strip():
            return 'not_run', 'no OUTCAR and no vasp.out', errs
    elif vd.contains(outcar, 'General timing and accounting'):
        if kind == 'relax' and not vd.contains(outcar, 'reached required accuracy'):
            return 'relax_unconv', 'OUTCAR finished without "reached required accuracy"', errs
        if kind in ('static', 'dfpt') and not vd.contains(outcar, 'aborting loop because EDIFF is reached'):
            return 'scf_nelm', 'OUTCAR finished without "EDIFF is reached"', errs
        if kind == 'dfpt' and not vd.contains(outcar, 'BORN EFFECTIVE CHARGES'):
            return 'no_born', 'DFPT finished without BORN EFFECTIVE CHARGES in OUTCAR', errs
        return 'ok', 'finished and converged', errs
    hit = vd.match_pattern(text)
    if hit:
        return hit[0], f'log: "{hit[2]}"', errs
    lim = re.search(r'DUE TO TIME LIMIT|time limit|' + vd.STOPPED, text, re.I)
    if lim:
        return 'timeout', f'log: "{lim.group(0)}"', errs
    return 'unknown', 'VASP stopped without a known error message', errs


CAUSE_TEXT = dict(vd.CAUSE_TEXT, no_born='DFPT finished but printed no Born effective charges '
                  '(no gap, or LEPSILON was not set)', crash='the program stopped')


def incar_tags(path):
    return vd.incar(path)


def vasp_fix(task, cause, incars):
    """(tags to set or None, held, what to do) for a VASP cause; incars: the INCARs concerned."""
    tags = {}
    for p in incars:
        tags.update(incar_tags(p))
    algo = tags.get('ALGO', '').lower()
    if cause in ('not_run', 'mpi_init', 'relax_unconv'):
        return {}, False, 'resubmit' + (' (the relaxation continues from CONTCAR)' if cause == 'relax_unconv' else '')
    if cause == 'timeout':
        return {}, False, 'walltime x2 in run_workflow.sh and resubmit_failed.sh, resubmit'
    if cause == 'scf_nelm':
        if algo == 'all' and int(tags.get('NELM', '60') or 60) >= 200:
            return None, True, 'ALGO = All and NELM >= 200 already: look at OSZICAR (charge sloshing? magnetism?)'
        return {'ALGO': 'All', 'NELM': '200'}, False, 'ALGO = All, NELM = 200'
    if cause == 'edddav':
        if algo == 'all':
            return None, True, 'ALGO = All already: check the structure (atoms too close?)'
        return {'ALGO': 'All'}, False, 'ALGO = All'
    if cause == 'zbrent' and task == 'relax':
        return {'IBRION': '1'}, False, 'IBRION = 1 (quasi-Newton), continue from CONTCAR'
    if cause == 'oom':
        if tags.get('KPAR', '1') == '1':
            return None, True, 'KPAR is 1 already: use fewer MPI ranks per node or more memory per rank'
        return {'KPAR': '1'}, False, 'KPAR = 1'
    if cause == 'no_born':
        return None, True, 'check the band gap of 02_scf; with no gap set DFPT: off'
    return None, True, 'read the error lines below'


# ── the diagnosis ─────────────────────────────────────────────────────────────
def diagnose(pr):
    done = pr.state.get('done', {})
    failed = pr.state.get('failed')
    status = {t: ('done' if t in done else 'FAILED' if t == failed else 'not reached') for t in TASKS}
    run = pr.running()
    findings = []
    if run:
        return status, run, findings
    joblog = pr.job_log()
    if not failed:
        nxt = next((t for t in TASKS if t not in done), None)
        if nxt is None:
            return status, '', findings
        if 'started' not in pr.log.rsplit('=== ', 1)[-1][:80]:
            return status, f'not started yet: cd {pr.path} && ' + ('sbatch' if pr.slurm else 'bash') + \
                ' run_workflow.sh', findings
        if True:
            # stopped inside a task, without the agent recording a failure: killed from outside
            status[nxt] = 'INTERRUPTED'
            lim = re.search(r'DUE TO TIME LIMIT|CANCELLED|time limit', joblog, re.I)
            oom = re.search(r'oom[-_ ]kill|out of memory|OUT_OF_MEMORY', joblog, re.I)
            if lim:
                findings.append(Finding(nxt, 'the job', 'timeout', vd.CAUSE_TEXT['timeout'], f'job log: "{lim.group(0)}"',
                                        fix=[('walltime', 2), ('resubmit', 'run_workflow.sh')], todo='walltime x2, resubmit'))
            elif oom:
                findings.append(Finding(nxt, 'the job', 'oom', vd.CAUSE_TEXT['oom'], f'job log: "{oom.group(0)}"',
                                        held=True, todo='more memory per rank (fewer ranks per node), then resubmit'))
            else:
                findings.append(Finding(nxt, 'the job', 'interrupted', 'the workflow stopped inside this task '
                                        'without an error (killed, or the computer restarted)',
                                        'report.log ends with "=== ' + nxt + ' -- started"',
                                        fix=[('resubmit', 'run_workflow.sh')], todo='resubmit (finished tasks are skipped)'))
        return status, '', findings
    msg = pr.failed_msg(failed)
    t = failed
    if t in VASP_DIR:
        d = pr.f(VASP_DIR[t])
        cause, ev, errs = vasp_cause(d, 'relax' if t == 'relax' else 'dfpt' if t == 'dfpt' else 'static', joblog)
        incar = os.path.join(d, 'INCAR')
        tags, held, todo = vasp_fix(t, cause, [incar] if os.path.isfile(incar) else [])
        fix = []
        if tags:
            fix += [('incar', incar, tags), ('block', BLOCK[t], tags)]
        if cause == 'timeout':
            fix.append(('walltime', 2))
        fix.append(('resubmit', 'run_workflow.sh'))
        findings.append(Finding(t, VASP_DIR[t], cause, CAUSE_TEXT.get(cause, ''), ev + ('; report.log: ' + msg if msg else ''),
                                errs, fix, held, todo))
    elif t == 'vasp_snaps':
        groups = {}
        for d in sorted(glob.glob(pr.f('04_supercells', 'snap_*'))):
            cause, ev, errs = vasp_cause(d, 'static', '')
            if cause != 'ok':
                groups.setdefault(cause, []).append((d, ev, errs))
        if not groups:
            findings.append(Finding(t, '04_supercells', 'none_now', 'every snapshot run is converged now',
                                    'report.log: ' + msg, fix=[('resubmit', 'run_workflow.sh')], todo='resubmit'))
        for cause, items in groups.items():
            incars = [os.path.join(d, 'INCAR') for d, _, _ in items]
            tags, held, todo = vasp_fix(t, cause, incars)
            fix = [('incar', p, tags) for p in incars] + ([('block', 'INCAR_snaps', tags)] if tags else [])
            if cause == 'timeout' or re.search(r'DUE TO TIME LIMIT', joblog):
                fix.append(('walltime', 2))
            fix.append(('resubmit', 'resubmit_failed.sh'))
            names = ' '.join(os.path.basename(d) for d, _, _ in items)
            findings.append(Finding(t, f'04_supercells: {names}', cause, CAUSE_TEXT.get(cause, ''), items[0][1],
                                    items[0][2], fix, held, todo))
    elif t == 'focex':
        findings += focex_findings(pr, msg)
    else:
        findings.append(other_finding(pr, t, msg))
    return status, '', findings


def focex_findings(pr, msg):
    out = []
    if 'binary not found' in msg:
        return [Finding('focex', '05_focex', 'program', 'the FOCEX binary is missing', msg, held=True,
                        todo='build it (cd FOCEX && make) or set FOCEX_X in the instruction file, then resubmit')]
    for d in [pr.f('05_focex')] + [pr.f('05_focex', s) for s in ('NOBORN', 'LASSO')]:
        fo = os.path.join(d, 'focex.out')
        if not os.path.isfile(fo):
            continue
        text = read(fo, tail=400000)
        logs = sorted(glob.glob(os.path.join(d, 'log*.dat')), key=os.path.getmtime)
        if logs:
            text += '\n' + read(logs[-1], tail=200000)
        if re.search(r'Percent error', text) and os.path.isfile(os.path.join(d, 'bands.dat')) \
                and os.path.isfile(os.path.join(d, 'KTICS.BS')):
            continue
        where = os.path.relpath(d, pr.path)
        for cause, rx, meaning, todo in FOCEX_PAT:
            hit = re.search(rx, text)
            if hit:
                line = next((l.strip() for l in text.splitlines() if hit.group(0) in l), hit.group(0))
                fix = [('reset', 'snaps'), ('resubmit', 'run_workflow.sh')] if cause == 'ws_outside' else []
                out.append(Finding('focex', where, cause, meaning, f'focex.out: "{line[:150]}"', held=True,
                                   todo=todo, fix=fix))
                break
        else:
            tail = [l for l in read(fo).splitlines() if l.strip()][-4:]
            out.append(Finding('focex', where, 'unknown', 'FOCEX stopped without a known message',
                               'report.log: ' + msg, tail, held=True, todo='read the end of focex.out and the log*.dat'))
    if not out:
        out.append(Finding('focex', '05_focex', 'unknown', 'FOCEX failed', 'report.log: ' + msg, held=True,
                           todo='read 05_focex/focex.out'))
    return out


def other_finding(pr, t, msg):
    rules = [('sc_snaps.x not found', 'program', 'the SC-SNAPS program is missing',
              'build it (make sc_snaps in Anharmonic-lattice-dynamics) or set SC_SNAPS_X'),
             ('read_outcar.x not found', 'program', 'read_outcar.x is missing',
              'build it (cd FOCEX/utility && make) or set READ_OUTCAR_X'),
             ('produced no poscar_000', 'sc_snaps', 'SC-SNAPS made no snapshots', 'read the SC-SNAPS output quoted above'),
             ('NOT commensurate', 'frame', 'the SC-SNAPS supercell does not fit POSCAR0', 'check cell.inp/supercell.inp blocks'),
             ('not on a POSCAR0 site', 'frame', 'an atom of the undistorted supercell is not on a POSCAR0 site',
              'check cell.inp/supercell.inp blocks'),
             ('not an automatic mesh', 'kpoints', 'the SCF KPOINTS is not an automatic mesh',
              'give a KPOINTS_snaps block in the instruction file'),
             ('not converged:', 'snaps_unconv', 'some supercell runs are not converged', 'diagnose vasp_snaps: '
              'rerun with ./resubmit_failed.sh after fixing INCAR_snaps'),
             ('pos-forc.dat', 'collect', 'read_outcar.x read a different number of atoms than the POSCAR has',
              'check the OUTCAR of that folder'),
             ('no KTICS.BS', 'plots', 'FOCEX did not finish, so there is nothing to plot', 'diagnose the focex task'),
             ('ERROR, see', 'plots', 'a gnuplot script failed', 'read the *.log next to the .plt in 05_focex')]
    for key, cause, meaning, todo in rules:
        if key in msg:
            fix = [('resubmit', 'resubmit_failed.sh')] if cause == 'snaps_unconv' else []
            return Finding(t, t, cause, meaning, 'report.log: ' + msg[:600], held=not fix, todo=todo, fix=fix)
    return Finding(t, t, 'unknown', 'no known cause', 'report.log: ' + msg[:600], held=True,
                   todo='read report.log and the files of this task')


# ── the report and fix_it.sh ──────────────────────────────────────────────────
def fingerprint(pr):
    out = []
    for f in ('state.json', 'report.log'):
        p = pr.f(f)
        out.append(f'{f}:{os.path.getsize(p)}:{int(os.path.getmtime(p))}' if os.path.isfile(p) else f'{f}:-')
    return ','.join(out)


def write_outputs(pr, status, run, findings):
    out = pr.f('diagnosis')
    os.makedirs(out, exist_ok=True)
    R = [f'diagnose.py {now()}  --  {pr.path}', f"mode {pr.state.get('mode', 'from_scratch')}, "
         f"{'SLURM' if pr.slurm else 'no batch system'}", '']
    R += ['task          status'] + [f'  {t:11s} {s}' for t, s in status.items()] + ['']
    if run.startswith('not started'):
        R.append(run)
    elif run:
        R.append(f'RUNNING: {run} -- nothing to fix now; diagnose again when it has stopped')
    elif not findings:
        warn = [l.split('] ', 1)[-1] for l in pr.log.splitlines() if 'WARNING' in l]
        R.append('all tasks done' if all(s == 'done' for s in status.values()) else 'nothing failed')
        if warn:
            R += ['', 'warnings in report.log (check them):'] + ['  ' + w[:200] for w in warn[-10:]]
    for k, f in enumerate(findings, 1):
        R += [f'[{k}] task {f.task}, {f.where}: {f.cause} -- {f.meaning}', f'    evidence: {f.evidence}']
        R += [f'    error line: {e}' for e in f.errs[:6]]
        R += [f"    {'HELD for review' if f.held else 'fix_it.sh'}: {f.todo}"]
        R += [f'    will set {", ".join(f"{k2} = {v}" for k2, v in x[2].items())} in {os.path.relpath(x[1], pr.path)}'
              for x in f.fix if x[0] == 'incar' and x[2]]
        R += [f'    will add {", ".join(f"{k2} = {v}" for k2, v in x[2].items())} to the {x[1]} block of '
              f'{os.path.basename(pr.instr)}' for x in f.fix if x[0] == 'block' and x[2]]
        R.append('')
    text = '\n'.join(R) + '\n'
    vd_write(os.path.join(out, 'report.txt'), text)
    fix = os.path.join(out, 'fix_it.sh')
    if os.path.isfile(fix):
        os.remove(fix)
    if findings and not run.startswith(('SLURM', 'process')):
        vd_write(fix, fix_script(pr, findings), 0o755)
    return text, (fix if os.path.isfile(fix) else None)


def vd_write(path, text, mode=None):
    with open(path, 'w') as fh:
        fh.write(text)
    if mode:
        os.chmod(path, mode)


def fix_script(pr, findings):
    py, me = sys.executable, os.path.abspath(__file__)
    held = [f for f in findings if f.held]
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    L = ['#!/bin/bash', f'# fix_it.sh -- written by diagnose.py {now()}; read diagnosis/report.txt first',
         'set -e', f'PROJ="{pr.path}"', f'PY="{py}"', f'DIAG="{me}"', 'cd "$PROJ"']
    if held:
        L += ['if [ "$1" != "--anyway" ]; then', '  echo "HELD for review -- nothing changed:"']
        L += [f'  echo "  [{f.task}, {f.where}] {f.cause}: {f.todo}"'.replace('$', '\\$') for f in held]
        L += ['  echo "see diagnosis/report.txt; after your changes run diagnose.py again"',
              '  echo "(bash diagnosis/fix_it.sh --anyway runs the commands below as they are)"', '  exit 1', 'fi']
    L += [f'"$PY" "$DIAG" _guard "$PROJ" "{fingerprint(pr)}"', f'B="$PROJ/diagnosis/backup_{stamp}"', 'mkdir -p "$B"',
          f'cp -p state.json "{os.path.basename(pr.instr)}" run_workflow.sh "$B"/',
          '[ -f resubmit_failed.sh ] && cp -p resubmit_failed.sh "$B"/']
    sub, seen = None, set()
    for f in findings:
        L.append(f'# [{f.task}, {f.where}] {f.cause}')
        for x in f.fix:
            if x[0] == 'incar' and x[2]:
                rel = os.path.relpath(x[1], pr.path)
                L += [f'mkdir -p "$B/{os.path.dirname(rel)}" && cp -p "{rel}" "$B/{rel}"',
                      f'"$PY" "$DIAG" _incar "{rel}" ' + ' '.join(f'"{k}={v}"' for k, v in x[2].items())]
            elif x[0] == 'block' and x[2]:
                L.append(f'"$PY" "$DIAG" _block "{os.path.basename(pr.instr)}" {x[1]} '
                         + ' '.join(f'"{k}={v}"' for k, v in x[2].items()))
            elif x[0] == 'walltime' and 'walltime' not in seen:
                seen.add('walltime')
                L.append(f'"$PY" "$DIAG" _walltime {x[1]} run_workflow.sh resubmit_failed.sh')
            elif x[0] == 'reset':
                L.append(f'"$PY" "$DIAG" _reset "$PROJ" {x[1]} "$B"   # tasks from {x[1]} on are redone; '
                         'their old folders go to the backup')
            elif x[0] == 'resubmit':
                sub = x[1] if sub != 'run_workflow.sh' else sub
    sub = sub or 'run_workflow.sh'
    L = [l for k, l in enumerate(L) if l.startswith('#') or l not in L[:k]]     # each command once
    L.append(f'"$PY" "$DIAG" _unfail "$PROJ"')
    if pr.slurm:
        L.append(f'sbatch {sub}')
    else:
        L += [f'nohup bash {sub} > {sub[:-3]}.log 2>&1 &', f'echo "started {sub} in the background (log: {sub[:-3]}.log)"']
    L.append('echo "done; diagnose again when the job has stopped"')
    return '\n'.join(L) + '\n'


# ── helpers called by fix_it.sh ───────────────────────────────────────────────
def h_guard(proj, fp):
    pr = Project(proj)
    if pr.running():
        sys.exit(f'fix_it.sh: {pr.running()} -- not changing anything while it runs')
    if fingerprint(pr) != fp:
        sys.exit('fix_it.sh: state.json or report.log changed since the diagnosis; run diagnose.py again')


def set_tags(text, tags):
    out, seen = [], set()
    for line in text.splitlines():
        t = line.split('=', 1)[0].strip().upper() if '=' in line and not line.strip().startswith(('#', '!')) else None
        if t in tags:
            if t not in seen:
                out.append(f'{t} = {tags[t]}   # diagnose.py {now()[:10]}')
            seen.add(t)
            continue
        out.append(line)
    out += [f'{t} = {v}   # diagnose.py {now()[:10]}' for t, v in tags.items() if t not in seen]
    return '\n'.join(out) + '\n'


def pairs(args):
    return dict(a.split('=', 1) for a in args)


def h_incar(path, args):
    vd_write(path, set_tags(open(path).read(), pairs(args)))


def h_block(path, block, args):
    """Set tags inside the block of the instruction file (made if missing); other lines kept."""
    tags, L = pairs(args), open(path).read().splitlines()
    key = lambda l: l.split('#')[0].strip().rstrip(':').strip().lower()
    i = next((k for k, l in enumerate(L) if key(l) == block.lower()), None)
    if i is None:
        L += ['', block] + [f'  {t} = {v}   # diagnose.py {now()[:10]}' for t, v in tags.items()] + [block + '_END']
    else:
        j = next(k for k in range(i + 1, len(L)) if key(L[k]) in (block.lower() + '_end', 'end_' + block.lower()))
        body = set_tags('\n'.join(L[i + 1:j]), tags).splitlines()
        L = L[:i + 1] + ['  ' + b.strip() for b in body if b.strip()] + L[j:]
    vd_write(path, '\n'.join(L) + '\n')


def h_walltime(factor, files):
    for f in files:
        if not os.path.isfile(f):
            continue
        t = open(f).read()
        m = re.search(r'(#SBATCH\s+(?:--time=|-t\s*))(\S+)', t)
        if m:
            new = vd.from_sec(vd.to_sec(m.group(2)) * float(factor))
            vd_write(f, t[:m.start(2)] + new + t[m.end(2):])
            print(f'{f}: walltime {m.group(2)} -> {new}')


def h_reset(proj, task, backup):
    """Forget the tasks from task on; their folders are moved into backup (an old supercell's
    converged runs would otherwise be taken as done)."""
    k = TASKS.index(task)
    for d, first in (('03_snaps', 'snaps'), ('04_supercells', 'supercells'), ('05_focex', 'focex')):
        if TASKS.index(first) >= k and os.path.isdir(os.path.join(proj, d)):
            shutil.move(os.path.join(proj, d), os.path.join(backup, d))
            print(f'{d} -> {backup}/{d}')
    st = os.path.join(proj, 'state.json')
    s = json.load(open(st))
    for t in TASKS[TASKS.index(task):]:
        s.get('done', {}).pop(t, None)
    s.pop('failed', None)
    vd_write(st, json.dumps(s, indent=1))


def h_unfail(proj):
    st = os.path.join(proj, 'state.json')
    s = json.load(open(st))
    s.pop('failed', None)
    vd_write(st, json.dumps(s, indent=1))


def main():
    helpers = {'_guard': lambda a: h_guard(a[0], a[1]), '_incar': lambda a: h_incar(a[0], a[1:]),
               '_block': lambda a: h_block(a[0], a[1], a[2:]), '_walltime': lambda a: h_walltime(a[0], a[1:]),
               '_reset': lambda a: h_reset(a[0], a[1], a[2]), '_unfail': lambda a: h_unfail(a[0])}
    if len(sys.argv) > 1 and sys.argv[1] in helpers:
        return helpers[sys.argv[1]](sys.argv[2:])
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0], epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('projects', nargs='+', help='project folders (PROJECT_DIR of the instruction file)')
    ap.add_argument('--version', action='version', version=f'diagnose.py {__version__}')
    args = ap.parse_args()
    for p in args.projects:
        pr = Project(p)
        status, run, findings = diagnose(pr)
        text, fix = write_outputs(pr, status, run, findings)
        print(text)
        if fix:
            print(f'fix script: {fix}' + ('  (HELD: prints what to do)' if any(f.held for f in findings) else ''))
        print(f'report: {pr.f("diagnosis", "report.txt")}\n')


if __name__ == '__main__':
    main()
