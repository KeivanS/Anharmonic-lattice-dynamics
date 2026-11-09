#!/usr/bin/env python3
"""bundle.py -- runs of a focex-agent project outside the agent (used by aladyn-gui; 2026-11-07).

With VASP_RUN: external (THERMACOND_RUN: external) the agent writes the inputs of a stage and stops
as WAITING (state.json 'waiting': task, folders, which folder waits for which). This module then

  prepare(proj)        job.sh in every waiting folder (SLURM header from hpc.json), submit_all.sh
                       (all independent runs submitted at once; DFPT after the SCF; collect after the
                       supercells), pack_results.sh, and for the supercells 04_supercells/collect.py;
  make_bundle(proj)    the same, packed as <project>_<task>_<date>.tar.gz (POTCAR links followed),
                       to copy to the cluster, unpack and run ./submit_all.sh;
  submit_here(proj)    the same, submitted with sbatch where the GUI runs (GUI on the cluster);
  import_results(proj, tar)   unpacks the results_*.tar.gz made by pack_results.sh into the project;
  absorb(proj)         marks the supercell tasks done when collect.py has written FORCEDISP1.

After that, 'focex-agent.py run <proj>' continues from the outputs. Only VASP for now; the per-code
parts are the job lines and the result files (DFT_CODES), so QE or SIESTA can be added there.
Everything run on the cluster is bash or python3 with the standard library only.
"""
import datetime
import glob
import json
import os
import re
import shlex
import shutil
import subprocess
import tarfile

HPC_DEFAULTS = {
    'partition': 'standard', 'account': '', 'nodes': 1, 'ntasks_per_node': 40, 'time': '12:00:00',
    'time_thermacond': '24:00:00', 'mem': '', 'mpi_cmd': 'srun', 'vasp_exec': 'vasp_std',
    'tcond_exec': '~/ALADYN/bin/tcond-2.2.x', 'python': 'python3',
    'modules': '',            # lines run before VASP, e.g. "module load vasp/6.4.2"
    'extra_sbatch': '',       # further #SBATCH lines
}

# per DFT code: the files a finished run leaves that the agent reads (pack_results.sh brings them back)
DFT_CODES = {
    'vasp': {'results': ['OUTCAR', 'CONTCAR', 'OSZICAR', 'EIGENVAL', 'IBZKPT', 'vasp.out', 'born_charges.txt',
                         'BORN', 'extract_born.out', 'slurm-*.out', 'slurm-*.err'],
             'skip': ['vasprun.xml', 'DOSCAR', 'PROCAR', 'XDATCAR', 'WAVECAR', 'CHG', 'PCDAT', 'REPORT']},
}
TC_SKIP = ['V3DELTA']                         # THERMACOND scratch (large), not packed


def profile_hpc():
    """The cluster settings of the VASP-Flow profile FOCEX_PROFILE (written to ~/.aladyn.env by ./configure on a
    cluster), in hpc.json's names; {} without one."""
    name = os.environ.get('FOCEX_PROFILE')
    vf = os.path.expanduser(os.environ.get('VASP_FLOW_DIR', '~/VASP-Flow'))
    f = os.path.join(vf, 'profiles', f'{name}.json') if name else ''
    if not (f and os.path.isfile(f)):
        return {}
    p = json.load(open(f))
    s = p.get('slurm') or {}
    h = {k: s[k] for k in ('partition', 'account', 'nodes', 'ntasks_per_node', 'time') if k in s}
    if p.get('mpi_cmd'):
        h['mpi_cmd'] = p['mpi_cmd']
    if p.get('vasp_std'):
        h['vasp_exec'] = p['vasp_std']
    if p.get('modules'):
        h['modules'] = '\n'.join(f'module load {m}' for m in p['modules'])
    return h


def load_hpc(proj):
    """hpc.json of the project; without it, the defaults and those of the FOCEX_PROFILE profile."""
    f = os.path.join(proj, 'hpc.json')
    h = dict(HPC_DEFAULTS)
    if os.path.isfile(f):
        h.update(json.load(open(f)))
    else:
        h.update(profile_hpc())
    return h


def save_hpc(proj, h):
    with open(os.path.join(proj, 'hpc.json'), 'w') as fh:
        json.dump({k: h.get(k, v) for k, v in HPC_DEFAULTS.items()}, fh, indent=1)


def ranks(h):
    return int(h['nodes']) * int(h['ntasks_per_node'])


def fit_kpar(d, n):
    """KPAR and NCORE of d/INCAR made to fit n MPI ranks (KPAR divides n, NCORE divides n/KPAR), at most the
    values there; the inputs may have been written for another number of ranks. Returns a note or ''."""
    f = os.path.join(d, 'INCAR')
    if not os.path.isfile(f):
        return ''
    L = open(f).read().splitlines()
    val = {}
    for l in L:
        m = re.match(r'\s*(KPAR|NCORE)\s*=\s*(\d+)', l, re.I)
        if m:
            val[m.group(1).upper()] = int(m.group(2))
    if not val:
        return ''
    k0, c0 = val.get('KPAR', 1), val.get('NCORE', 1)
    k = max(x for x in range(1, k0 + 1) if n % x == 0)
    c = max(x for x in range(1, c0 + 1) if (n // k) % x == 0)
    if (k, c) == (k0, c0):
        return ''
    out = []
    for l in L:
        m = re.match(r'\s*(KPAR|NCORE)\s*=', l, re.I)
        out.append(f'{m.group(1).upper()} = {k if m.group(1).upper() == "KPAR" else c}   ! aladyn-gui: fitted to {n} ranks'
                   if m else l)
    with open(f, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return f'{os.path.basename(d)}: KPAR {k0}->{k}, NCORE {c0}->{c} for {n} ranks'


def state(proj):
    return json.load(open(os.path.join(proj, 'state.json')))


def sbatch_header(h, name, time=None, ntasks=None, nodes=None):
    L = ['#!/bin/bash', f'#SBATCH --job-name={name}', f"#SBATCH --partition={h['partition']}",
         f"#SBATCH --nodes={nodes or h['nodes']}",
         f"#SBATCH --ntasks-per-node={ntasks or h['ntasks_per_node']}", f"#SBATCH --time={time or h['time']}",
         '#SBATCH --output=slurm-%j.out', '#SBATCH --error=slurm-%j.err']
    if h.get('account'):
        L.append(f"#SBATCH --account={h['account']}")
    if h.get('mem'):
        L.append(f"#SBATCH --mem={h['mem']}")
    L += [l for l in h.get('extra_sbatch', '').splitlines() if l.strip()]
    L += ['', 'cd "$SLURM_SUBMIT_DIR"', 'export OMP_NUM_THREADS=1', 'rm -f STOPCAR']
    L += [l for l in h.get('modules', '').splitlines() if l.strip()]
    return '\n'.join(L) + '\n'


def mpi_line(h, exe, out):
    if h['mpi_cmd'].strip() == 'srun':
        return f'srun {exe} > {out} 2>&1'
    return f"{h['mpi_cmd']} -np {ranks(h)} {exe} > {out} 2>&1"


GAP_CHECK = r'''#!/usr/bin/env python3
"""gap_check.py (aladyn-gui): exit 1 when the SCF of this DFPT run (copy_from_scf.sh) shows no band gap
(DFPT: auto skips the run then), 0 otherwise or when it cannot tell."""
import os, re, sys
here = os.path.dirname(os.path.abspath(__file__))
scf = os.path.join(here, '..', '02_scf')
cp = os.path.join(here, 'copy_from_scf.sh')
if os.path.isfile(cp):
    m = re.search(r'^SCF_DIR="\$HERE/([^"]+)"', open(cp).read(), re.M)
    if m:
        scf = os.path.join(here, m.group(1))
try:
    L = open(os.path.join(scf, 'EIGENVAL')).read().split('\n')
    nk, nb = [int(float(x)) for x in L[5].split()[1:3]]
    occ, emp, i = [], [], 6
    for _ in range(nk):
        while not L[i].strip():
            i += 1
        i += 1
        for _ in range(nb):
            v = [float(x) for x in L[i].split()[1:]]
            i += 1
            for e, o in ([(v[0], v[1])] if len(v) == 2 else [(v[0], v[2]), (v[1], v[3])]):
                (occ if o > 0.5 else emp).append(e)
    gap = min(emp) - max(occ)
except Exception as err:
    print(f'gap_check.py: gap not readable ({err}); DFPT runs')
    sys.exit(0)
print(f'gap_check.py: SCF gap {gap:.3f} eV')
sys.exit(0 if gap > 0.05 else 1)
'''

COLLECT = r'''#!/usr/bin/env python3
"""collect.py (aladyn-gui) -- run in 04_supercells after the supercell VASP runs (submit_all.sh does it
as a job after them). Does what focex-agent's collect task does with read_outcar.x, in python3 with the
standard library: for every snap_NNN (not unused/), in order, the ionic step(s) of OUTCAR -- t-th
'TOTAL-FORCE' block with the t-th 'free  energy   TOTEN  =' line, as read_outcar.x -- into FORCEDISP1
(' # POSITION      TOTAL FORCE', 't TOTEN', then x y z Fx Fy Fz per atom, Angstrom and eV/A);
POSCAR1 = 03_snaps/poscar_000 (else snap_000/POSCAR); stress.params from snap_000 (GPa, tensile > 0 =
-VASP 'in kB'/10, only with ISIF >= 2); collect.log with the checks (converged, atoms, net force,
OUTCAR positions vs POSCAR). Nothing is written if a run did not converge."""
import glob, math, os, re, shutil, sys

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)
log = []


def say(s):
    print(s)
    log.append(s)


def poscar(path):
    L = open(path).read().splitlines()
    s = float(L[1].split()[0])
    A = [[s * float(x) for x in L[k].split()[:3]] for k in (2, 3, 4)]
    n = sum(int(x) for x in L[6].split())
    k = 7 + (L[7].strip()[:1] in ('s', 'S'))
    direct = L[k].strip()[:1] in 'dD'
    X = []
    for l in L[k + 1:k + 1 + n]:
        f = [float(x) for x in l.split()[:3]]
        X.append([sum(f[j] * A[j][i] for j in range(3)) for i in range(3)] if direct else [s * x for x in f])
    return A, X


def inv3(A):
    a, b, c = A
    det = (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))
    cof = [[b[1] * c[2] - b[2] * c[1], a[2] * c[1] - a[1] * c[2], a[1] * b[2] - a[2] * b[1]],
           [b[2] * c[0] - b[0] * c[2], a[0] * c[2] - a[2] * c[0], a[2] * b[0] - a[0] * b[2]],
           [b[0] * c[1] - b[1] * c[0], a[1] * c[0] - a[0] * c[1], a[0] * b[1] - a[1] * b[0]]]
    return [[cof[i][j] / det for j in range(3)] for i in range(3)]


def steps(outcar):
    L = open(outcar, errors='replace').read().splitlines()
    n = next(int(l[l.index('NIONS =') + 7:].split()[0]) for l in L if 'NIONS =' in l[:99])
    tf = [i for i, l in enumerate(L) if 'TOTAL-FORCE' in l[:99]]
    te = [l for l in L if 'free  energy   TOTEN  =' in l[:99]]
    out = []
    for t in range(min(len(tf), len(te))):
        rows = [[float(x) for x in L[tf[t] + 2 + i].split()[:6]] for i in range(n)]
        e = float(te[t].split('=')[1].split()[0])
        out.append((e, rows))
    return n, out


dirs = sorted(d for d in glob.glob('snap_*') if os.path.isdir(d))
bad = []
for d in dirs:
    f = os.path.join(d, 'OUTCAR')
    t = open(f, errors='replace').read() if os.path.isfile(f) else ''
    if 'General timing' not in t:
        bad.append(f'{d}: VASP did not finish' if t else f'{d}: no OUTCAR')
    elif 'EDIFF is reached' not in t:
        bad.append(f'{d}: electronic SCF not converged (NELM)')
if not dirs or bad:
    say('collect.py: ' + ('no snap_* folders' if not dirs else f'{len(bad)} of {len(dirs)} runs not usable:'))
    for b in bad:
        say('  ' + b)
    say('  FORCEDISP1 not written; rerun those folders (cd snap_NNN && sbatch job.sh), then python3 collect.py')
    open('collect.log', 'w').write('\n'.join(log) + '\n')
    sys.exit(1)
blocks, worst_sum, worst_pos = [], 0.0, 0.0
for d in dirs:
    n, st = steps(os.path.join(d, 'OUTCAR'))
    A, X = poscar(os.path.join(d, 'POSCAR'))
    B = inv3(A)
    if not st or len(X) != n:
        say(f'collect.py: {d}: {n} atoms in OUTCAR, {len(X)} in POSCAR, {len(st)} steps -- stop')
        open('collect.log', 'w').write('\n'.join(log) + '\n')
        sys.exit(1)
    for t, (e, rows) in enumerate(st, 1):
        blocks.append(' # POSITION      TOTAL FORCE \n' + f' {t} {e:.14f}\n' +
                      ''.join('  ' + ' '.join(f'{v:18.10f}' for v in r[:3]) + '  ' +
                              ' '.join(f'{v:18.10E}' for v in r[3:6]) + '\n' for r in rows))
    rows = st[-1][1]
    worst_sum = max(worst_sum, math.sqrt(sum(sum(r[3 + k] for r in rows) ** 2 for k in range(3))))
    for r, x in zip(rows, X):
        dc = [r[k] - x[k] for k in range(3)]
        df = [sum(dc[i] * B[i][j] for i in range(3)) for j in range(3)]
        df = [v - round(v) for v in df]
        worst_pos = max(worst_pos, max(abs(sum(df[j] * A[j][i] for j in range(3))) for i in range(3)))
open('FORCEDISP1', 'w').write(''.join(blocks))
p0 = os.path.join('..', '03_snaps', 'poscar_000')
shutil.copy(p0 if os.path.isfile(p0) else os.path.join('snap_000', 'POSCAR'), 'POSCAR1')
say(f'collect.py: FORCEDISP1 from {len(dirs)} OUTCARs; POSCAR1; max |sum F| = {worst_sum:.2e} eV/A; '
    f'max |r_OUTCAR - r_POSCAR| = {worst_pos:.1e} A' + (' -- WARNING: check these' if worst_sum > 0.05 or worst_pos > 1e-3 else ''))
t = open(os.path.join('snap_000', 'OUTCAR'), errors='replace').read() if os.path.isdir('snap_000') else ''
m = re.search(r'ISIF\s*=\s*(\d+)', t)
kb = re.findall(r'^\s*in kB' + r'\s+([-\d.Ee+*]+)' * 6, t, re.M)
if m and int(m.group(1)) >= 2 and kb:
    xx, yy, zz, xy, yz, zx = (float(x) for x in kb[-1])
    S = [[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]]
    open('stress.params', 'w').write(''.join(' '.join(f'{-v / 10:12.6f}' for v in r) + '\n' for r in S) +
        '# stress of the undistorted supercell (snap_000/OUTCAR), GPa, tensile > 0 (= -VASP in kB/10)\n')
    say(f'collect.py: stress.params from snap_000: diagonal {[round(-S[i][i] / 10, 3) for i in range(3)]} GPa')
else:
    say('collect.py: no full stress tensor in snap_000/OUTCAR (ISIF < 2): stress.params not written')
open('collect.log', 'w').write('\n'.join(log) + '\n')
'''


def waiting(proj):
    w = state(proj).get('waiting')
    if not w:
        raise RuntimeError('nothing is waiting to run: continue the project first (it stops at the next stage to run)')
    return w


def prepare(proj, h=None):
    """job.sh in the waiting folders, submit_all.sh and pack_results.sh at the project top.
    Returns (task, folders relative to proj)."""
    h = h or load_hpc(proj)
    w = waiting(proj)
    name = os.path.basename(os.path.normpath(proj))
    dirs, after, kind = w['dirs'], w.get('after', {}), w.get('kind', 'vasp')
    st = state(proj)
    instr = open(os.path.join(proj, 'focex-agent-instructions.txt')).read()
    dfpt_auto = not re.search(r'^\s*DFPT\s*:\s*on', instr, re.I | re.M)
    py = h['python']
    for rel in dirs:
        d = os.path.join(proj, rel)
        if kind != 'thermacond':
            fit_kpar(d, ranks(h))
        if kind == 'thermacond':
            body = sbatch_header(h, f'{name}_tcond', h['time_thermacond']) + \
                mpi_line(h, h['tcond_exec'], 'thermacond.out') + '\n'
        else:
            pre, post = [], []
            if os.path.basename(rel).endswith('dfpt'):
                pre.append('[ -f copy_from_scf.sh ] && bash copy_from_scf.sh')
                if dfpt_auto:
                    with open(os.path.join(d, 'gap_check.py'), 'w') as fh:
                        fh.write(GAP_CHECK)
                    pre.append(f'{py} gap_check.py || {{ echo "DFPT skipped: no gap in the SCF"; exit 0; }}')
                pre.append(f'[ -f set_efield.py ] && {py} set_efield.py')
                post.append(f'[ -f OUTCAR ] && {py} extract_born.py OUTCAR > extract_born.out 2>&1')
            tag = rel.replace('04_supercells/', '').replace('/', '_')
            body = (sbatch_header(h, f'{name}_{tag}') + ''.join(l + '\n' for l in pre) +
                    mpi_line(h, h['vasp_exec'], 'vasp.out') + '\n' + ''.join(l + '\n' for l in post))
        with open(os.path.join(d, 'job.sh'), 'w') as fh:
            fh.write(body)
        os.chmod(os.path.join(d, 'job.sh'), 0o755)
    snaps = w['task'] == 'vasp_snaps'
    if snaps:
        with open(os.path.join(proj, '04_supercells', 'collect.py'), 'w') as fh:
            fh.write(COLLECT)
        with open(os.path.join(proj, '04_supercells', 'collect_job.sh'), 'w') as fh:
            fh.write(sbatch_header(h, f'{name}_collect', '00:30:00', ntasks=1, nodes=1) + f'{py} collect.py\n')
    # submit_all.sh: in the order given; a folder in 'after' waits for its folder (afterok)
    S = ['#!/bin/bash', f'# aladyn-gui: submit the {w["task"]} runs of {name} ({len(dirs)} job(s))',
         'cd "$(dirname "$0")"', 'set -e', 'ids=""']
    var = {}
    for k, rel in enumerate(dirs):
        dep = f' --dependency=afterok:${var[after[rel]]}' if rel in after and after[rel] in var else ''
        var[rel] = f'j{k}'
        S.append(f'j{k}=$(cd {shlex.quote(rel)} && sbatch --parsable{dep} job.sh); echo "{rel}: job $j{k}"; ids="$ids:$j{k}"')
    if snaps:
        S.append('jc=$(cd 04_supercells && sbatch --parsable --dependency=afterany$ids collect_job.sh); '
                 'echo "collect (FORCEDISP1, after all runs): job $jc"')
    S.append('echo "when the jobs are done: ./pack_results.sh, and import the results_*.tar.gz in aladyn-gui"')
    _write_x(os.path.join(proj, 'submit_all.sh'), '\n'.join(S) + '\n')
    # pack_results.sh: the files the agent reads, per folder; for the supercells only collect.py's output
    # and the small logs (--outcars adds the OUTCARs)
    keep = DFT_CODES['vasp']['results']
    small = ['OSZICAR', 'vasp.out', 'slurm-*.out', 'slurm-*.err']
    P = ['#!/bin/bash', '# aladyn-gui: pack the results of the waiting runs (small files only; --outcars: also',
         '# the supercell OUTCARs)', 'cd "$(dirname "$0")/.."', f'p={shlex.quote(name)}',
         'out="$p/results_${p}_$(date +%Y%m%d-%H%M%S).tar.gz"', 'files=()',
         'dirs="' + ' '.join(dirs) + '"']
    if kind == 'thermacond':
        P.append('for d in $dirs; do for f in "$p/$d"/*; do case "$(basename "$f")" in '
                 + '|'.join(TC_SKIP) + ') ;; *) [ -f "$f" ] && files+=("$f");; esac; done; done')
    elif snaps:
        P.append('names="' + ' '.join(small) + '"; [ "$1" = --outcars ] && names="$names OUTCAR"')
        P.append('for d in $dirs; do for k in $names; do for f in "$p/$d"/$k; do [ -f "$f" ] && files+=("$f"); done; done; done')
        P.append('for f in $p/04_supercells/{FORCEDISP1,POSCAR1,stress.params,collect.log}; do [ -f "$f" ] && files+=("$f"); done')
        P.append('[ -f "$p/04_supercells/FORCEDISP1" ] || echo "WARNING: no 04_supercells/FORCEDISP1 yet (see collect.log; '
                 'python3 collect.py in 04_supercells)"')
    else:
        P.append('names="' + ' '.join(keep) + '"')
        P.append('for d in $dirs; do for k in $names; do for f in "$p/$d"/$k; do [ -f "$f" ] && files+=("$f"); done; done; done')
    P += ['tar czf "$out" "${files[@]}"', 'echo "wrote $out (${#files[@]} files): copy it back and import it in aladyn-gui"']
    _write_x(os.path.join(proj, 'pack_results.sh'), '\n'.join(P) + '\n')
    st.setdefault('bundles', []).append({'task': w['task'], 'dirs': dirs, 'prepared': _now()})
    _save_state(proj, st)
    return w['task'], dirs


def make_bundle(proj, h=None):
    """prepare() and pack: <proj>/bundles/<name>_<task>_<date>.tar.gz with the project's top folder name,
    the waiting folders (inputs; POTCAR links followed), submit_all.sh, pack_results.sh, a README."""
    task, dirs = prepare(proj, h)
    name = os.path.basename(os.path.normpath(proj))
    os.makedirs(os.path.join(proj, 'bundles'), exist_ok=True)
    out = os.path.join(proj, 'bundles', f'{name}_{task}_{datetime.datetime.now():%Y%m%d-%H%M%S}.tar.gz')
    skip = DFT_CODES['vasp']['skip'] + TC_SKIP + ['bundles']
    readme = (f'aladyn-gui bundle: {name}, task {task}, {len(dirs)} folder(s), {_now()}\n\n'
              f'  tar xzf {os.path.basename(out)}        # unpack where the project should live (e.g. scratch)\n'
              f'  cd {name} && ./submit_all.sh\n'
              '  squeue -u $USER                     # wait until the jobs are gone\n'
              '  ./pack_results.sh                   # -> results_*.tar.gz: copy it back, import it in aladyn-gui\n'
              + ('\nThe supercell runs are submitted together; a last job runs 04_supercells/collect.py, which writes\n'
                 'FORCEDISP1, POSCAR1 and stress.params (see collect.log).\n' if task == 'vasp_snaps' else '') +
              ('\nDFPT starts after the SCF (job dependency) and reads its CHGCAR: keep 02_scf and 06_dfpt together.\n'
               if any(d.endswith('dfpt') for d in dirs) else ''))
    with open(os.path.join(proj, 'README_HPC.txt'), 'w') as fh:
        fh.write(readme)
    with tarfile.open(out, 'w:gz', dereference=True) as tar:
        for f in ('submit_all.sh', 'pack_results.sh', 'README_HPC.txt', 'focex-agent-instructions.txt'):
            tar.add(os.path.join(proj, f), arcname=f'{name}/{f}')
        if task == 'vasp_snaps' and os.path.isfile(os.path.join(proj, '03_snaps', 'poscar_000')):
            tar.add(os.path.join(proj, '03_snaps', 'poscar_000'), arcname=f'{name}/03_snaps/poscar_000')
        if task == 'vasp_snaps':
            for f in ('collect.py', 'collect_job.sh'):
                tar.add(os.path.join(proj, '04_supercells', f), arcname=f'{name}/04_supercells/{f}')
        for rel in dirs:
            d = os.path.join(proj, rel)
            for f in sorted(os.listdir(d)):
                p = os.path.join(d, f)
                if f in skip or re.match(r'OUTCAR\.\d+$', f) or f == 'OUTCAR' or (os.path.isdir(p) and f in TC_SKIP):
                    continue
                if os.path.islink(p) and not os.path.exists(p):
                    continue
                tar.add(p, arcname=f'{name}/{rel}/{f}')
    return out


def submit_here(proj, h=None):
    """prepare() and run submit_all.sh with sbatch on this machine (the GUI runs on the cluster)."""
    if shutil.which('sbatch') is None:
        raise RuntimeError('sbatch not found here: make a bundle and submit it on the cluster')
    prepare(proj, h)
    r = subprocess.run(['bash', 'submit_all.sh'], cwd=proj, capture_output=True, text=True)
    return r.returncode == 0, r.stdout + r.stderr


def import_results(proj, tar_path):
    """Unpack a results_*.tar.gz (paths <project name>/...) into proj; returns the files written."""
    name = os.path.basename(os.path.normpath(proj))
    written = []
    with tarfile.open(tar_path, 'r:gz') as tar:
        members = [m for m in tar.getmembers() if m.isfile()]
        # the top folder is the project's name on the cluster (it may differ: renamed, cloned): dropped when
        # every member is under one folder that is not a step folder (NN_name)
        tops = {m.name.split('/')[0] for m in members}
        strip = len(tops) == 1 and all('/' in m.name for m in members) and not re.match(r'\d\d_', next(iter(tops)))
        for m in members:
            parts = m.name.split('/')
            if strip:
                parts = parts[1:]
            rel = os.path.normpath('/'.join(parts))
            if rel.startswith('..') or os.path.isabs(rel):
                continue
            dst = os.path.join(proj, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with tar.extractfile(m) as src, open(dst, 'wb') as fh:
                fh.write(src.read())
            written.append(rel)
    return written


def absorb(proj):
    """The supercell runs are waiting and 04_supercells/FORCEDISP1 (collect.py) is newer than their inputs:
    mark vasp_snaps and collect done, so that the agent goes on to FOCEX. Returns a message or ''."""
    st = state(proj)
    w = st.get('waiting') or {}
    fd = os.path.join(proj, '04_supercells', 'FORCEDISP1')
    if w.get('task') != 'vasp_snaps' or not os.path.isfile(fd):
        return ''
    newest = max(os.path.getmtime(os.path.join(proj, rel, 'POSCAR')) for rel in w['dirs'])
    if os.path.getmtime(fd) < newest:
        return ''
    n = open(fd).read().count('POSITION')
    for t in ('vasp_snaps', 'collect'):
        st.setdefault('done', {})[t] = 'provided (collect.py on the cluster)'
    st.pop('waiting', None)
    _save_state(proj, st)
    return f'FORCEDISP1 with {n} snapshots from the cluster: supercell runs and collect marked done'


def _write_x(path, text):
    with open(path, 'w') as fh:
        fh.write(text)
    os.chmod(path, 0o755)


def _now():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def _save_state(proj, st):
    with open(os.path.join(proj, 'state.json'), 'w') as fh:
        json.dump(st, fh, indent=1, default=str)
