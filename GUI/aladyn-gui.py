#!/usr/bin/env python3
"""aladyn-gui -- browser interface for the ALADYN suite: from a POSCAR to force constants (FOCEX), phonons
and the thermal conductivity (THERMACOND). 2026-11-07.

    aladyn-gui [--workdir DIR] [--port 5070] [--host 127.0.0.1]
    then open http://localhost:5070  (on a cluster: ssh -L 5070:localhost:5070 user@cluster, same URL)

A project is a focex-agent project folder (FOCEX/AGENT/focex-agent.py; manual FOCEX/AGENT/
focexagent_manual.pdf): the GUI writes its instruction file, runs its tasks one by one, and stops before
FOCEX and THERMACOND so that the input files can be reviewed. Each DFT stage (relaxation, SCF + DFPT,
supercells) and THERMACOND can run
  * here       -- on this computer (the agent starts VASP / THERMACOND),
  * as a bundle -- a .tar.gz with job scripts for a SLURM cluster (unpack, ./submit_all.sh,
                  ./pack_results.sh, import the results here), or
  * submitted  -- with sbatch from here, when the GUI itself runs on the cluster.
The bundles are made by FOCEX/AGENT/bundle.py. The SC-SNAPS displacement histograms come from
SC_SNAPS/sc-snaps-gui.py.
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
import threading
import time

from flask import Flask, Response, jsonify, request, send_file

sys.dont_write_bytecode = True
GUI_DIR = os.path.dirname(os.path.realpath(__file__))
REPO = os.path.dirname(GUI_DIR)
AGENT_PY = os.path.join(REPO, 'FOCEX', 'AGENT', 'focex-agent.py')
VERSION = '0.1'


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fa = _load(AGENT_PY, 'focex_agent')                                  # needs VASP-Flow, as the agent
bd = _load(os.path.join(REPO, 'FOCEX', 'AGENT', 'bundle.py'), 'aladyn_bundle')
try:
    snapgui = _load(os.path.join(REPO, 'SC_SNAPS', 'sc-snaps-gui.py'), 'sc_snaps_gui')
except Exception:                                                     # histograms are optional
    snapgui = None

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 ** 3
WORKDIR = os.getcwd()
INSTR = fa.INSTR
TASKS = fa.TASKS
CHECKPOINTS = ('focex', 'thermacond')     # the GUI stops before these until the user starts them
PARAMS = ('structure_params', 'latdyn_params', 'kpbs_params', 'dielectric_params', 'default_params')
JOBS = {}                                 # project path -> job dict
LOCK = threading.Lock()

STEPS = [  # (title, tasks, what it is)
    ('Project', [], 'POSCAR, settings, instruction file'),
    ('Relaxation', ['relax', 'poscar0'], 'VASP relaxation of the primitive cell; POSCAR0 = the cell of every later step'),
    ('SCF and DFPT', ['scf', 'dfpt'], 'band gap; Born charges and epsilon_inf for an insulator'),
    ('Snapshots', ['snaps', 'supercells'], 'SC-SNAPS supercell snapshots and one VASP folder per snapshot'),
    ('Supercell DFT', ['vasp_snaps', 'collect'], 'the supercell VASP runs; FORCEDISP1 and POSCAR1'),
    ('FOCEX', ['focex', 'plots'], 'force constants, phonons, Gruneisen, elastic constants; plots'),
    ('THERMACOND', ['thermacond'], 'linewidths and lattice thermal conductivity'),
]


# ── instruction file ───────────────────────────────────────────────────────────
def instr_path(proj):
    return os.path.join(proj, INSTR)


def read_instr(proj):
    return open(instr_path(proj)).read()


def write_instr(proj, text):
    with open(instr_path(proj), 'w') as fh:
        fh.write(text)


def set_keyword(text, key, value):
    """KEY: value (replaced where it is, else appended); value None removes it."""
    rx = re.compile(r'^[ \t]*' + re.escape(key) + r'[ \t]*[:=].*$', re.I | re.M)
    if value is None:
        return rx.sub('', text)
    line = f'{key}: {value}'
    if rx.search(text):
        return rx.sub(line, text, count=1)
    return text.rstrip('\n') + '\n' + line + '\n'


def get_keyword(text, key, default=''):
    return fa.keyword(fa.split_blocks(text)[1], key, default)


def set_block(text, name, body):
    """Block NAME ... NAME_END (replaced, appended, or removed with body None)."""
    rx = re.compile(r'^[ \t]*' + name + r'[ \t]*:?[ \t]*\n.*?^[ \t]*(?:' + name + r'_END|END_' + name + r')[ \t]*\n?',
                    re.I | re.M | re.S)
    text = rx.sub('', text)
    if body is None:
        return text
    return text.rstrip('\n') + f'\n{name}\n{body.rstrip(chr(10))}\n{name}_END\n'


def get_blocks(proj):
    return fa.split_blocks(read_instr(proj))[0]


# ── projects ──────────────────────────────────────────────────────────────────
def proj_path(name_or_path):
    p = os.path.expanduser(name_or_path or '')
    p = p if os.path.isabs(p) else os.path.join(WORKDIR, p)
    p = os.path.realpath(p)
    if not os.path.isfile(os.path.join(p, 'state.json')):
        raise ValueError(f'{p} is not a project folder (no state.json)')
    return p


def load_state(proj):
    return json.load(open(os.path.join(proj, 'state.json')))


def save_state(proj, st):
    with open(os.path.join(proj, 'state.json'), 'w') as fh:
        json.dump(st, fh, indent=1, default=str)


def gui_cfg(proj):
    f = os.path.join(proj, 'gui.json')
    c = {'local_mpi': max(1, (os.cpu_count() or 2) - 2), 'dft_code': 'vasp'}
    if os.path.isfile(f):
        c.update(json.load(open(f)))
    return c


def save_gui_cfg(proj, c):
    with open(os.path.join(proj, 'gui.json'), 'w') as fh:
        json.dump(c, fh, indent=1)


def list_projects():
    out = []
    for f in sorted(glob.glob(os.path.join(WORKDIR, '*', 'state.json')) + glob.glob(os.path.join(WORKDIR, '*', '*', 'state.json'))):
        d = os.path.dirname(f)
        out.append({'path': d, 'name': os.path.relpath(d, WORKDIR), 'mtime': os.path.getmtime(f)})
    return sorted(out, key=lambda x: -x['mtime'])


def instruction_text(s):
    """The instruction file from the setup form (dict s). Each comment is on its own line above its keyword:
    VASP-Flow reads a keyword's whole line, so a comment after FUNCTIONAL: R2SCAN naming other functionals
    chose one of those (seen 2026-11-07: PBEsol)."""
    def kv(comment, key, value):
        return [f'# {comment}', f'{key}: {value}']
    L = [f"Project: {s['name']}", 'MODE: from_scratch', f"PROJECT_DIR: {s['slug']}",
         '# written by aladyn-gui; every keyword is explained in FOCEX/AGENT/focexagent_manual.pdf']
    L += kv('exchange-correlation functional: PBE, PBEsol, LDA, R2SCAN or HSE06', 'FUNCTIONAL', s.get('functional', 'R2SCAN'))
    L += kv('Accurate (ENCUT 600 eV, EDIFF 1E-7) or Normal', 'ACCURACY', s.get('accuracy', 'Accurate'))
    L += kv('Born charges: auto = when the SCF shows a band gap; on; off', 'DFPT', s.get('dfpt', 'auto'))
    L += kv('MPI ranks of the VASP runs (the GUI sets it for each stage: here or on the cluster)', 'MPI', s.get('mpi', 1))
    L += kv('force-constant ranks 1..8 fitted', 'RANKS', s.get('ranks', '1 1 1 0 0 0 0 0'))
    L += kv('neighbour shells for ranks 2..8', 'RANGE', s.get('range', '6 2 1 0 0 0 0'))
    L += kv('0: FC2 range from the supercell', 'FC2RANGE', s.get('fc2range', 0))
    L += kv('SMALL, LARGE or 9 integers (supercell vectors as rows, primitive units)', 'SUPERCELL', s.get('supercell', 'SMALL'))
    L += kv('number of supercell snapshots', 'SNAPSHOTS', s.get('snapshots', 50))
    if str(s.get('temperature', '')).strip():
        L += kv('temperature of the snapshots (K); without it 30 K tuned to DISP_TARGET', 'TEMPERATURE', s['temperature'])
    L += kv('rms displacement (A) aimed at when TEMPERATURE is not given', 'DISP_TARGET', s.get('disp_target', 0.02))
    L += kv('THERMACOND: 0 none, 1 linewidths, 2 kappa direct, 3 kappa iterative', 'THERMACOND_JOB', s.get('tc_job', 0))
    if str(s.get('tc_mesh', '')).strip():
        L += kv('q mesh of THERMACOND', 'THERMACOND_MESH', s['tc_mesh'])
    if str(s.get('tc_T', '')).strip():
        L += kv('THERMACOND temperatures: Tmin Tmax number', 'THERMACOND_T', s['tc_T'])
    L += kv('where VASP / THERMACOND run: here, or external (bundle or submitted); set by the GUI',
            'VASP_RUN', 'here') + ['THERMACOND_RUN: here']
    extra = s.get('extra', '').strip()
    if extra:
        L += ['# your own lines (VASP-Flow keywords, INCAR blocks, ...)', extra]
    return '\n'.join(L) + '\n'


def agent_cmd(*args):
    return [sys.executable, AGENT_PY] + list(args)


def env():
    e = dict(os.environ)
    e.setdefault('OPENBLAS_NUM_THREADS', '1')
    return e


# ── params drafts ─────────────────────────────────────────────────────────────
def agent_for(proj):
    """An Agent on proj whose POSCAR0 is POSCAR.init until the relaxation made POSCAR0."""
    a = fa.Agent(proj)
    if not os.path.isfile(os.path.join(proj, 'POSCAR0')):
        real = a.d
        a.d = lambda *p: real(*(('POSCAR.init',) if p == ('POSCAR0',) else p))
    return a


def generated_params(proj):
    """The five files as the agent would write them now (blocks ignored), with comments that mark what is
    a placeholder. Returns ({name: text}, notes)."""
    a = agent_for(proj)
    a.fx['blocks'] = {k: v for k, v in a.fx['blocks'].items() if k not in PARAMS}
    notes = []
    unrelaxed = not os.path.isfile(os.path.join(proj, 'POSCAR0'))
    out = {}
    s = a.structure_params('1 1 1 1').splitlines()
    if unrelaxed:
        s[0] = s[0].split('#')[0].rstrip() + '   # a b c alpha beta gamma -- PLACEHOLDER: unrelaxed POSCAR, refreshed after the relaxation'
        notes.append('structure.params: lattice and positions from the unrelaxed POSCAR; refreshed from POSCAR0 after the relaxation (your other lines are kept)')
    s[6] = s[6].split('#')[0].rstrip() + '   # number of POSCARi/FORCEDISPi pairs (1: one supercell data set), verbosity'
    out['structure_params'] = '\n'.join(s) + '\n'
    out['latdyn_params'] = a.latdyn_params()
    k = a.kpbs_params().splitlines()
    out['kpbs_params'] = '\n'.join(k[:3] + [l + '   # high-symmetry point: label, reduced coordinates' for l in k[3:]]) + '\n'
    st = load_state(proj)
    born = 1 if st.get('Z') else 0
    d = a.dielectric_params(born).splitlines()
    if not st.get('Z'):
        d[0] = '0   # born_flag: PLACEHOLDER -- no Born charges yet (DFPT in step 3, or type eps_inf and Z* here and set 1)'
        d = [d[0]] + [l if '#' in l else l + '   # PLACEHOLDER' for l in d[1:]]
        notes.append('dielectric.params: placeholders (zeros, born_flag 0) until DFPT gives the Born charges, or you type them')
    out['dielectric_params'] = '\n'.join(d) + '\n'
    out['default_params'] = open(os.path.join(REPO, 'FOCEX', 'default.params')).read()
    return out, notes


def merge_structure(user, gen):
    """structure.params: the lattice/atom lines of gen (POSCAR0 now), the control lines of user (ranks,
    invariances, itemp, pairs, fc2range, shells) when the number of atoms is the same; else gen."""
    U, G = user.rstrip('\n').splitlines(), gen.rstrip('\n').splitlines()
    try:
        if len(U) != len(G) or U[10].split()[0] != G[10].split()[0]:
            return gen, False
    except IndexError:
        return gen, False
    for i in (3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 18):
        G[i] = U[i]
    return '\n'.join(G) + '\n', True


def params_now(proj):
    """{name: {'text', 'source', 'generated'}} and notes; a structure_params block older than POSCAR0 is
    refreshed (geometry from POSCAR0, the user's control lines kept) and saved again."""
    gen, notes = generated_params(proj)
    blocks = get_blocks(proj)
    cfg = gui_cfg(proj)
    out = {}
    for n in PARAMS:
        if n in blocks:
            text, src = blocks[n], 'yours (saved)'
            if n == 'structure_params' and os.path.isfile(os.path.join(proj, 'POSCAR0')) and \
                    cfg.get('params_saved', 0) < os.path.getmtime(os.path.join(proj, 'POSCAR0')):
                text, ok = merge_structure(text, gen[n])
                write_instr(proj, set_block(read_instr(proj), n, text))
                cfg['params_saved'] = time.time()
                save_gui_cfg(proj, cfg)
                notes.append('structure.params: geometry refreshed from the relaxed POSCAR0' +
                             (', your other lines kept' if ok else ' (atom count changed: regenerated)'))
            if n == 'dielectric_params' and load_state(proj).get('Z'):
                notes.append('dielectric.params: your saved file is used, not the DFPT result')
            out[n] = {'text': text, 'source': src, 'generated': gen[n]}
        else:
            out[n] = {'text': gen[n], 'source': 'generated', 'generated': gen[n]}
    return out, notes


def tc_params_now(proj):
    blocks = get_blocks(proj)
    a = fa.Agent(proj)
    if os.path.isfile(os.path.join(proj, '05_focex', 'latdyn.params')):
        a.fx['blocks'].pop('thermacond_params', None)
        gen = '\n'.join(a.tc_latdyn()) + '\n'
    else:                                   # before FOCEX: its latdyn.params draft with the THERMACOND_* lines
        b = agent_for(proj)
        b.fx['blocks'].pop('thermacond_params', None)
        gen = '\n'.join(b.tc_latdyn(b.latdyn_params())) + '\n'
    if 'thermacond_params' in blocks:
        return {'text': blocks['thermacond_params'], 'source': 'yours (saved)', 'generated': gen}
    return {'text': gen, 'source': 'generated' + ('' if os.path.isfile(os.path.join(proj, '05_focex', 'latdyn.params'))
                                                  else ' (template; FOCEX not run yet)'), 'generated': gen}


# ── guide ─────────────────────────────────────────────────────────────────────
def guide(proj):
    st = load_state(proj)
    done, w, failed = st.get('done', {}), st.get('waiting') or {}, st.get('failed')
    steps, nxt = [], None
    for title, tasks, what in STEPS:
        if not tasks:
            steps.append({'title': title, 'status': 'done', 'what': what, 'detail': ''})
            continue
        if all(t in done for t in tasks):
            prov = [t for t in tasks if str(done[t]).startswith('provided')]
            steps.append({'title': title, 'status': 'done', 'what': what,
                          'detail': ('provided: ' + str(done[prov[0]])[9:].strip('() ')) if prov else ''})
            continue
        if w.get('task') in tasks:
            steps.append({'title': title, 'status': 'waiting', 'what': what,
                          'detail': f"inputs ready in {len(w['dirs'])} folder(s) since {w.get('since', '')}"})
        elif failed in tasks:
            steps.append({'title': title, 'status': 'failed', 'what': what, 'detail': f'task {failed} failed'})
        elif nxt is None:
            steps.append({'title': title, 'status': 'next', 'what': what, 'detail': ''})
        else:
            steps.append({'title': title, 'status': 'later', 'what': what, 'detail': ''})
        if nxt is None:
            nxt = (title, tasks)
    tc_job = fa.parse_instructions(instr_path(proj))[1]['tc_job']
    advice = next_advice(proj, st, steps, tc_job)
    return {'steps': steps, 'advice': advice}


def next_advice(proj, st, steps, tc_job):
    done, w, failed = st.get('done', {}), st.get('waiting') or {}, st.get('failed')
    sbatch = shutil.which('sbatch') is not None
    if JOBS.get(proj, {}).get('running'):
        return ['A job is running (see the Log tab). The guide updates when it ends.']
    if failed:
        last = [l for l in tail(os.path.join(proj, 'report.log'), 60) if f'=== {failed} -- FAILED' in l]
        return [f'Task {failed} failed: ' + (last[-1].split('FAILED', 1)[1][:400] if last else 'see the Log tab'),
                'Fix the cause (Setup: instruction file, or the folder\'s INCAR), then press "Continue".',
                f'diagnose.py writes an explanation and a fix script: {os.path.join(REPO, "FOCEX", "AGENT", "diagnose.py")} {proj}']
    if w:
        kind = 'THERMACOND' if w.get('kind') == 'thermacond' else 'VASP'
        L = [f"{kind} inputs of task {w['task']} are ready in {len(w['dirs'])} folder(s).",
             '1. Run them: "Make bundle" (download, copy to the cluster, tar xzf, ./submit_all.sh)' +
             (', or "Submit here" (sbatch on this machine)' if sbatch else '') + ', or "Run here".',
             '2. When the jobs are done: on the cluster ./pack_results.sh, copy results_*.tar.gz back and "Import results"' +
             (' -- or, when the GUI runs on the cluster, just "Continue"' if sbatch else '') + '.']
        if w['task'] == 'vasp_snaps':
            L.append('The last cluster job runs collect.py: FORCEDISP1, POSCAR1 and stress.params come back, not the OUTCARs.')
        return L
    for s in steps:
        if s['status'] == 'next':
            t = s['title']
            if t == 'Relaxation' and 'relax' in done:
                return ['Next: press "Continue": the relaxed cell becomes POSCAR0, the cell of every later step.']
            if t == 'Relaxation':
                return ['Next: relax the structure. Choose where VASP runs (Run here / Make bundle / Submit here) on the Run tab.',
                        'If your POSCAR is already relaxed with these same settings, use "Skip: POSCAR is relaxed" instead.']
            if t == 'SCF and DFPT':
                return ['Next: SCF (band gap) and, for an insulator, DFPT (Born charges, epsilon_inf). Choose where they run on the Run tab.']
            if t == 'Snapshots':
                return ['Next: make the snapshots and the supercell VASP folders (fast, runs here): press "Continue".',
                        'Check SUPERCELL, SNAPSHOTS and TEMPERATURE on the Setup tab first.']
            if t == 'Supercell DFT':
                return ['Next: the supercell VASP runs -- the expensive step. Usually "Make bundle" (or "Submit here" on the cluster).']
            if t == 'FOCEX':
                return ['Next: FOCEX. Review the five input files on the Params tab (placeholders are marked PLACEHOLDER), save,',
                        'then press "Run FOCEX" (runs here, minutes). The plots follow; see the Results tab.']
            if t == 'THERMACOND':
                if tc_job == 0:
                    return ['FOCEX is done. For the thermal conductivity set line 8 of the THERMACOND latdyn.params to 2 or 3',
                            '(THERMACOND tab), save, then "Run THERMACOND here", or make a bundle for the cluster.']
                return ['Next: THERMACOND. Review its latdyn.params (THERMACOND tab), then run it here (a workstation is enough',
                        'for moderate meshes) or make a bundle / submit it on the cluster.']
    return ['All steps are done. Results: Results tab (plots) and the Log tab (report.log).',
            'To redo with another temperature or supercell: "New project from this one" on the Projects tab.']


def tail(path, n):
    if not os.path.isfile(path):
        return []
    return open(path, errors='replace').read().splitlines()[-n:]


# ── jobs ──────────────────────────────────────────────────────────────────────
def start_job(proj, kind, target):
    with LOCK:
        if JOBS.get(proj, {}).get('running'):
            raise RuntimeError('a job is already running for this project')
        job = {'running': True, 'kind': kind, 'lines': [], 'result': '', 'started': time.time()}
        JOBS[proj] = job
    threading.Thread(target=target, args=(proj, job), daemon=True).start()
    return job


def run_logged(job, cmd, cwd):
    job['lines'].append('$ ' + ' '.join(cmd))
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env())
    for line in p.stdout:
        job['lines'].append(line.rstrip())
        del job['lines'][:-2000]
    return p.wait()


def drive(proj, job, allow=(), post=None):
    """Run the next tasks one by one (focex-agent run --only), stopping at a waiting stage (then post:
    'bundle' or 'submit'), at a failure, or before a checkpoint task not in allow."""
    try:
        while True:
            msg = bd.absorb(proj)
            if msg:
                job['lines'].append(msg)
            st = load_state(proj)
            done = st.get('done', {})
            nxt = next((t for t in TASKS if t not in done), None)
            if nxt is None:
                job['result'] = 'all tasks done'
                return
            # a checkpoint waits for the user, unless it is the waiting stage (its results may be back)
            if nxt in CHECKPOINTS and nxt not in allow and (st.get('waiting') or {}).get('task') != nxt:
                job['result'] = f'stopped before {nxt}: review its input files, then start it'
                return
            if nxt == 'focex':
                params_now(proj)                     # refresh a saved structure.params from POSCAR0
            rc = run_logged(job, agent_cmd('run', proj, '--only', nxt), proj)
            if rc == 3:
                job['result'] = f'{nxt}: inputs ready, waiting for the runs'
                if post == 'bundle':
                    out = bd.make_bundle(proj)
                    job['bundle'] = os.path.relpath(out, proj)
                    job['result'] += f'; bundle {os.path.basename(out)}'
                elif post == 'submit':
                    ok, out = bd.submit_here(proj)
                    job['lines'] += out.splitlines()
                    job['result'] += '; submitted' if ok else '; sbatch FAILED (see the log)'
                return
            if rc != 0:
                job['result'] = f'{nxt} failed (see the log and the guide)'
                return
            if nxt in allow and nxt == 'thermacond':
                job['result'] = 'THERMACOND done'
                return
    except Exception as e:                       # shown in the GUI
        job['result'] = f'error: {e}'
        job['lines'].append(f'ERROR: {e}')
    finally:
        job['running'] = False


def set_mode(proj, vasp=None, tc=None):
    """VASP_RUN / THERMACOND_RUN and the MPI ranks for the next generated inputs (KPAR follows MPI)."""
    t = read_instr(proj)
    if vasp:
        t = set_keyword(t, 'VASP_RUN', vasp)
        n = bd.ranks(bd.load_hpc(proj)) if vasp == 'external' else gui_cfg(proj)['local_mpi']
        t = set_keyword(t, 'MPI', n)
    if tc:
        t = set_keyword(t, 'THERMACOND_RUN', tc)
        if tc == 'here':                         # the agent starts THERMACOND with MPI ranks
            t = set_keyword(t, 'MPI', gui_cfg(proj)['local_mpi'])
    write_instr(proj, t)


# ── API ───────────────────────────────────────────────────────────────────────
def ok(**k):
    return jsonify(ok=True, **k)


def err(msg, code=400):
    return jsonify(ok=False, error=str(msg)), code


@app.route('/')
def index():
    return Response(HTML.replace('__WORKDIR__', WORKDIR).replace('__VERSION__', VERSION), mimetype='text/html')


@app.route('/api/projects')
def api_projects():
    return ok(projects=list_projects(), workdir=WORKDIR, sbatch=shutil.which('sbatch') is not None)


@app.route('/api/new_project', methods=['POST'])
def api_new_project():
    s = request.json or {}
    try:
        name = (s.get('name') or '').strip()
        if not name or not (s.get('poscar') or '').strip():
            return err('a project name and a POSCAR are needed')
        s['name'], s['slug'] = name, re.sub(r'[^A-Za-z0-9_.-]+', '_', name)
        proj = os.path.join(WORKDIR, s['slug'])
        if os.path.exists(os.path.join(proj, 'state.json')):
            return err(f'{proj} exists: choose another name, or open it')
        setup = os.path.join(WORKDIR, '.aladyn_setup_' + s['slug'])
        os.makedirs(setup, exist_ok=True)
        s['mpi'] = s.get('local_mpi') or 1
        with open(os.path.join(setup, INSTR), 'w') as fh:
            fh.write(instruction_text(s))
        with open(os.path.join(setup, 'POSCAR'), 'w') as fh:
            fh.write(s['poscar'].rstrip('\n') + '\n')
        # PROJECT_DIR is relative to the instruction file: ../<slug>
        t = open(os.path.join(setup, INSTR)).read().replace(f"PROJECT_DIR: {s['slug']}", f"PROJECT_DIR: ../{s['slug']}")
        open(os.path.join(setup, INSTR), 'w').write(t)
        r = subprocess.run(agent_cmd('--instr', INSTR, '--coord', 'POSCAR'), cwd=setup, capture_output=True,
                           text=True, env=env())
        shutil.rmtree(setup, ignore_errors=True)
        if r.returncode != 0 or not os.path.isfile(os.path.join(proj, 'state.json')):
            return err('focex-agent could not set up the project:\n' + (r.stdout + r.stderr)[-2000:])
        write_instr(proj, read_instr(proj).replace(f"PROJECT_DIR: ../{s['slug']}", 'PROJECT_DIR: .'))
        cfg = gui_cfg(proj)
        cfg['local_mpi'] = int(s.get('local_mpi') or cfg['local_mpi'])
        save_gui_cfg(proj, cfg)
        if s.get('relaxed'):
            skip_relax(proj)
        return ok(project=proj, log=(r.stdout + r.stderr)[-3000:])
    except Exception as e:
        return err(e, 500)


def skip_relax(proj):
    """The POSCAR is relaxed: 01_relax/CONTCAR = POSCAR.init, task relax provided."""
    shutil.copy(os.path.join(proj, 'POSCAR.init'), os.path.join(proj, '01_relax', 'CONTCAR'))
    st = load_state(proj)
    st.setdefault('done', {})['relax'] = 'provided (relaxed POSCAR given)'
    save_state(proj, st)


@app.route('/api/skip_relax', methods=['POST'])
def api_skip_relax():
    try:
        proj = proj_path((request.json or {}).get('project'))
        skip_relax(proj)
        return ok()
    except Exception as e:
        return err(e)


@app.route('/api/clone', methods=['POST'])
def api_clone():
    """A new project from an existing one: the relaxed cell, SCF and DFPT are kept (tasks up to dfpt
    provided), the snapshots are made again with the new TEMPERATURE / SUPERCELL / SNAPSHOTS."""
    d = request.json or {}
    try:
        src = proj_path(d.get('project'))
        name = re.sub(r'[^A-Za-z0-9_.-]+', '_', (d.get('name') or '').strip())
        if not name:
            return err('a name for the new project is needed')
        dst = os.path.join(WORKDIR, name)
        if os.path.exists(dst):
            return err(f'{dst} exists')
        os.makedirs(dst)
        for f in (INSTR, 'POSCAR.input', 'POSCAR.init', 'POSCAR0', 'hpc.json', 'gui.json'):
            if os.path.isfile(os.path.join(src, f)):
                shutil.copy(os.path.join(src, f), dst)
        if os.path.exists(os.path.join(src, 'POTCAR')):
            shutil.copy(os.path.realpath(os.path.join(src, 'POTCAR')), os.path.join(dst, 'POTCAR'))
        big = {'WAVECAR', 'CHG', 'CHGCAR', 'vasprun.xml', 'DOSCAR', 'PROCAR', 'XDATCAR', 'WFULL', 'TMPCAR'}
        for sub in ('01_relax', '02_scf', '06_dfpt', os.path.join('02_scf', 'dfpt')):
            s_ = os.path.join(src, sub)
            if os.path.isdir(s_):
                os.makedirs(os.path.join(dst, sub), exist_ok=True)
                for f in os.listdir(s_):
                    p = os.path.join(s_, f)
                    if os.path.isfile(p) and f not in big:
                        shutil.copy(os.path.realpath(p), os.path.join(dst, sub, f))
        st0 = load_state(src)
        st = {k: st0[k] for k in ('profile', 'mode', 'sg', 'gap', 'eps', 'Z', 'reuse_relax') if k in st0}
        st['created'] = fa.now()
        st['cloned_from'] = src
        st['done'] = {t: f'provided (from {os.path.basename(src)})' for t in ('relax', 'poscar0', 'scf', 'dfpt')
                      if t in st0.get('done', {})}
        save_state(dst, st)
        t = read_instr(dst)
        for key in ('TEMPERATURE', 'SUPERCELL', 'SNAPSHOTS'):
            v = str(d.get(key.lower(), '')).strip()
            if v:
                t = set_keyword(t, key, v)
        t = set_keyword(t, 'Project', name)
        write_instr(dst, t)
        open(os.path.join(dst, 'report.log'), 'a').write(f'[{fa.now()}] aladyn-gui: new project from {src}\n')
        return ok(project=dst)
    except Exception as e:
        return err(e, 500)


@app.route('/api/status')
def api_status():
    try:
        proj = proj_path(request.args.get('project'))
    except Exception as e:
        return err(e)
    st = load_state(proj)
    job = JOBS.get(proj, {})
    pdfs = sorted(os.path.relpath(f, proj) for f in glob.glob(os.path.join(proj, '05_focex', '**', '*.pdf'), recursive=True) +
                  glob.glob(os.path.join(proj, '07_thermacond', '*.pdf')) + glob.glob(os.path.join(proj, '06_thermacond', '*.pdf')))
    bundles = sorted(os.path.relpath(f, proj) for f in glob.glob(os.path.join(proj, 'bundles', '*.tar.gz')))
    try:
        g = guide(proj)
    except SystemExit as e:                      # the instruction file does not parse
        g = {'steps': [], 'advice': [f'instruction file: {e}']}
    return ok(project=proj, name=os.path.basename(proj), state=st, guide=g, report=tail(os.path.join(proj, 'report.log'), 150),
              job={k: job.get(k) for k in ('running', 'kind', 'result', 'bundle')} | {'lines': job.get('lines', [])[-300:]},
              pdfs=pdfs, bundles=bundles, hpc=bd.load_hpc(proj), cfg=gui_cfg(proj), sbatch=shutil.which('sbatch') is not None,
              instr=read_instr(proj), potcar=os.path.exists(os.path.join(proj, 'POTCAR')),
              snaps=os.path.isfile(os.path.join(proj, '03_snaps', 'snapshots.xyz')))


@app.route('/api/action', methods=['POST'])
def api_action():
    """continue | run_here | bundle | submit | focex | tc_here | tc_bundle | tc_submit | make_bundle_now"""
    d = request.json or {}
    try:
        proj = proj_path(d.get('project'))
        a = d.get('action')
        st = load_state(proj)
        if a in ('bundle', 'submit') and st.get('waiting'):          # already prepared: just pack / submit
            if a == 'bundle':
                return ok(bundle=os.path.relpath(bd.make_bundle(proj), proj))
            r, out = bd.submit_here(proj)
            return ok(output=out) if r else err(out)
        if a == 'run_here':
            set_mode(proj, vasp='here')
            n = gui_cfg(proj)['local_mpi']
            for rel in (st.get('waiting') or {}).get('dirs', []):     # inputs written for the cluster
                bd.fit_kpar(os.path.join(proj, rel), n)
            st.pop('waiting', None)
            save_state(proj, st)
            start_job(proj, a, lambda p, j: drive(p, j))
        elif a in ('bundle', 'submit'):
            set_mode(proj, vasp='external')
            start_job(proj, a, lambda p, j: drive(p, j, post=a))
        elif a == 'continue':
            # never starts VASP here: external mode only writes inputs or reads results that are back
            set_mode(proj, vasp='external')
            start_job(proj, a, lambda p, j: drive(p, j))
        elif a == 'focex':
            start_job(proj, a, lambda p, j: drive(p, j, allow=('focex',)))
        elif a in ('tc_here', 'tc_bundle', 'tc_submit'):
            set_mode(proj, tc='here' if a == 'tc_here' else 'external')
            st.get('done', {}).pop('thermacond', None)
            if a == 'tc_here':
                st.pop('waiting', None)
            else:                                # outputs of an earlier run would be taken as the result
                tc = fa.Agent(proj).tc_dir()
                if os.path.isfile(os.path.join(tc, 'thermacond.out')):
                    os.replace(os.path.join(tc, 'thermacond.out'), os.path.join(tc, 'thermacond.out.prev'))
            save_state(proj, st)
            post = {'tc_bundle': 'bundle', 'tc_submit': 'submit'}.get(a)
            start_job(proj, a, lambda p, j: drive(p, j, allow=('focex', 'thermacond'), post=post))
        else:
            return err(f'unknown action {a}')
        return ok()
    except Exception as e:
        return err(e)


@app.route('/api/import', methods=['POST'])
def api_import():
    try:
        proj = proj_path(request.form.get('project'))
        f = request.files['file']
        os.makedirs(os.path.join(proj, 'bundles'), exist_ok=True)
        path = os.path.join(proj, 'bundles', 'imported_' + os.path.basename(f.filename))
        f.save(path)
        files = bd.import_results(proj, path)
        msg = bd.absorb(proj)
        open(os.path.join(proj, 'report.log'), 'a').write(
            f'[{fa.now()}] aladyn-gui: imported {len(files)} files from {os.path.basename(f.filename)}' + (f'; {msg}' if msg else '') + '\n')
        return ok(files=files, message=msg)
    except Exception as e:
        return err(e)


@app.route('/api/hpc', methods=['POST'])
def api_hpc():
    d = request.json or {}
    try:
        proj = proj_path(d.get('project'))
        h = bd.load_hpc(proj)
        h.update({k: v for k, v in (d.get('hpc') or {}).items() if k in bd.HPC_DEFAULTS})
        bd.save_hpc(proj, h)
        cfg = gui_cfg(proj)
        if d.get('local_mpi'):
            cfg['local_mpi'] = int(d['local_mpi'])
        save_gui_cfg(proj, cfg)
        return ok()
    except Exception as e:
        return err(e)


@app.route('/api/instr', methods=['POST'])
def api_instr():
    d = request.json or {}
    try:
        proj = proj_path(d.get('project'))
        new = d.get('text', '')
        tmp = instr_path(proj) + '.check'
        open(tmp, 'w').write(new)
        try:
            fa.parse_instructions(tmp)
        except SystemExit as e:
            return err(f'not saved: {e}')
        finally:
            os.remove(tmp)
        write_instr(proj, new)
        return ok()
    except Exception as e:
        return err(e)


@app.route('/api/params')
def api_params():
    try:
        proj = proj_path(request.args.get('project'))
        p, notes = params_now(proj)
        return ok(params=p, notes=notes, tc=tc_params_now(proj))
    except SystemExit as e:
        return err(f'{e}')
    except Exception as e:
        return err(e, 500)


@app.route('/api/save_params', methods=['POST'])
def api_save_params():
    """Files changed from the generated ones are kept as blocks of the instruction file (the agent writes
    them as they are); a file equal to the generated one removes its block."""
    d = request.json or {}
    try:
        proj = proj_path(d.get('project'))
        t = read_instr(proj)
        saved = []
        if 'params' in d:
            gen, _ = generated_params(proj)
            for n, text in d['params'].items():
                if n not in PARAMS:
                    continue
                same = text.strip() == gen[n].strip()
                t = set_block(t, n, None if same else text)
                if not same:
                    saved.append(n)
        if 'tc' in d:
            gen = tc_params_now(proj)['generated']
            text = d['tc']
            same = text.strip() == gen.strip()
            t = set_block(t, 'thermacond_params', None if same else text)
            job = text.splitlines()[7].split()[0] if len(text.splitlines()) > 7 else '0'
            t = set_keyword(t, 'THERMACOND_JOB', job)       # line 8 decides; the keyword would win otherwise
            if not same:
                saved.append('thermacond_params')
        tmp = instr_path(proj) + '.check'
        open(tmp, 'w').write(t)
        try:
            fa.parse_instructions(tmp)
        except SystemExit as e:
            return err(f'not saved: {e}')
        finally:
            os.remove(tmp)
        write_instr(proj, t)
        cfg = gui_cfg(proj)
        cfg['params_saved'] = time.time()
        save_gui_cfg(proj, cfg)
        return ok(saved=saved)
    except Exception as e:
        return err(e)


@app.route('/api/file')
def api_file():
    try:
        proj = proj_path(request.args.get('project'))
        p = os.path.realpath(os.path.join(proj, request.args.get('path', '')))
        if not p.startswith(proj + os.sep) or not os.path.isfile(p):
            return err('no such file in the project', 404)
        return send_file(p, as_attachment=request.args.get('download') == '1')
    except Exception as e:
        return err(e)


@app.route('/api/histogram/<kind>')
def api_histogram(kind):
    if snapgui is None:
        return err('SC_SNAPS/sc-snaps-gui.py could not be loaded', 500)
    try:
        proj = proj_path(request.args.get('project'))
        d = os.path.join(proj, '03_snaps')
        fn = {'disp': snapgui._make_displacement_histogram_png,
              'comp': snapgui._make_displacement_component_histogram_png}[kind]
        return Response(fn(os.path.join(d, 'poscar_000'), os.path.join(d, 'snapshots.xyz')), mimetype='image/png')
    except Exception as e:
        return err(e, 500)


# ── page ──────────────────────────────────────────────────────────────────────
HTML = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>ALADYN GUI</title>
<style>
:root{--bg:#f0f4f8;--card:#fff;--border:#e2e8f0;--accent:#0369a1;--accent-light:#e0f2fe;--text:#1e293b;--sub:#64748b;
--term-bg:#0f172a;--term-fg:#e2e8f0;--ok:#15803d;--warn:#b45309;--bad:#b91c1c;}
*{box-sizing:border-box;margin:0;padding:0}
body{font:14px/1.5 system-ui,sans-serif;background:var(--bg);color:var(--text)}
header{background:var(--accent);color:#fff;padding:12px 24px;display:flex;align-items:baseline;gap:14px;flex-wrap:wrap}
header h1{font-size:20px}header span{font-size:13px;opacity:.8}
.wrap{display:flex;gap:18px;padding:18px;max-width:1500px;margin:0 auto;align-items:flex-start}
.guide{width:330px;flex:none;position:sticky;top:12px}
.mainp{flex:1;min-width:0}
.card{background:var(--card);border:1px solid var(--border);border-radius:10px;padding:16px;margin-bottom:16px}
.ct{font-size:11.5px;font-weight:700;text-transform:uppercase;letter-spacing:.06em;color:var(--sub);margin-bottom:10px}
.step{display:flex;gap:10px;padding:7px 0;border-bottom:1px solid var(--border)}
.step:last-child{border:none}.dot{width:22px;height:22px;border-radius:50%;flex:none;display:flex;align-items:center;justify-content:center;color:#fff;font-size:12px;font-weight:700}
.done .dot{background:var(--ok)}.waiting .dot{background:var(--warn)}.failed .dot{background:var(--bad)}.next .dot{background:var(--accent)}.later .dot{background:#cbd5e1}
.step b{font-size:13.5px}.step small{display:block;color:var(--sub);font-size:12px}
.advice{background:var(--accent-light);border-radius:8px;padding:10px;font-size:13px}.advice p{margin:4px 0}
.tabs{display:flex;gap:4px;flex-wrap:wrap;margin-bottom:12px}
.tab{padding:7px 14px;border:1px solid var(--border);border-radius:8px;background:#fff;cursor:pointer;font-weight:600;color:var(--sub)}
.tab.on{background:var(--accent);color:#fff;border-color:var(--accent)}
.pane{display:none}.pane.on{display:block}
label{font-size:12.5px;color:var(--sub);display:block;margin-bottom:2px}
input,select,textarea{font:13px system-ui;border:1px solid var(--border);border-radius:6px;padding:6px 8px;width:100%;background:#fff}
textarea{font-family:ui-monospace,Menlo,monospace;font-size:12.5px;line-height:1.45}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(210px,1fr));gap:10px 14px}
.btn{padding:7px 16px;border:1px solid var(--border);border-radius:7px;background:#fff;cursor:pointer;font-weight:600;margin:4px 6px 4px 0}
.btn.p{background:var(--accent);color:#fff;border-color:var(--accent)}.btn:disabled{opacity:.45;cursor:default}
.term{background:var(--term-bg);color:var(--term-fg);font:12px ui-monospace,Menlo,monospace;padding:10px;border-radius:8px;max-height:460px;overflow:auto;white-space:pre-wrap}
.note{font-size:12.5px;color:var(--warn);margin:6px 0}.muted{color:var(--sub);font-size:12.5px}
.ptabs .tab{font-size:12.5px;padding:5px 10px}
iframe{width:100%;height:640px;border:1px solid var(--border);border-radius:8px}
.row{display:flex;gap:10px;flex-wrap:wrap;align-items:flex-end}
.pill{display:inline-block;padding:1px 8px;border-radius:10px;background:var(--accent-light);color:var(--accent);font-size:12px;font-weight:600}
@media(max-width:900px){.wrap{flex-direction:column}.guide{width:100%;position:static}}
</style></head><body>
<header><h1>ALADYN GUI</h1><span>FOCEX · SC-SNAPS · THERMACOND — v__VERSION__</span><span id="hdrproj"></span></header>
<div class="wrap">
<aside class="guide">
 <div class="card"><div class="ct">Guide: steps in order</div><div id="steps"><p class="muted">Open or create a project.</p></div></div>
 <div class="card"><div class="ct">What to do next</div><div class="advice" id="advice"><p>Create a project from a POSCAR (Projects tab), or open one.</p></div></div>
</aside>
<main class="mainp">
<div class="tabs" id="tabs"></div>

<section class="pane" id="p-projects">
 <div class="card"><div class="ct">Open a project</div>
  <p class="muted">Working folder: __WORKDIR__</p>
  <div class="row"><div style="flex:1"><select id="plist"></select></div><button class="btn p" onclick="openProj()">Open</button></div>
 </div>
 <div class="card"><div class="ct">New project from a POSCAR</div>
  <div class="grid">
   <div><label>Project name</label><input id="n_name" placeholder="e.g. Cu2O"></div>
   <div><label>Functional</label><select id="n_functional"><option>R2SCAN</option><option>PBE</option><option>PBEsol</option><option>LDA</option><option>HSE06</option></select></div>
   <div><label>Accuracy</label><select id="n_accuracy"><option>Accurate</option><option>Normal</option></select></div>
   <div><label>DFPT (Born charges)</label><select id="n_dfpt"><option>auto</option><option>on</option><option>off</option></select></div>
   <div><label>DFT code</label><select id="n_code"><option value="vasp">VASP</option><option disabled>Quantum ESPRESSO (later)</option><option disabled>SIESTA (later)</option></select></div>
   <div><label>Force-constant ranks</label><select id="n_ranks"><option value="1 1 1 0 0 0 0 0">harmonic + cubic (THERMACOND)</option><option value="1 1 0 0 0 0 0 0">harmonic only</option><option value="1 1 1 1 0 0 0 0">harmonic + cubic + quartic</option></select></div>
   <div><label>Shells for ranks 2..8 (RANGE)</label><input id="n_range" value="6 2 1 0 0 0 0"></div>
   <div><label>Supercell</label><input id="n_supercell" value="SMALL" title="SMALL, LARGE or 9 integers"></div>
   <div><label>Snapshots</label><input id="n_snapshots" value="50"></div>
   <div><label>Temperature (K, empty: tuned)</label><input id="n_temperature"></div>
   <div><label>THERMACOND job</label><select id="n_tc"><option value="0">0 none</option><option value="1">1 linewidths</option><option value="2">2 kappa, direct</option><option value="3">3 kappa, iterative</option></select></div>
   <div><label>MPI ranks here (this computer)</label><input id="n_mpi" value="8"></div>
  </div>
  <div style="margin-top:10px"><label>POSCAR (paste, or choose a file)</label><input type="file" id="n_file" onchange="readPos()"><textarea id="n_poscar" rows="10" placeholder="POSCAR"></textarea></div>
  <div style="margin-top:8px"><label><input type="checkbox" id="n_relaxed" style="width:auto"> This POSCAR is already relaxed with these same DFT settings (skip the relaxation)</label></div>
  <div style="margin-top:8px"><label>More instruction lines (VASP-Flow keywords such as MAGMOM, GGA_U, ENCUT; INCAR blocks) — optional</label><textarea id="n_extra" rows="3"></textarea></div>
  <button class="btn p" onclick="newProj()">Create project</button><span id="n_msg" class="muted"></span>
 </div>
 <div class="card"><div class="ct">New project from the open one</div>
  <p class="muted">Keeps the relaxed cell, SCF and DFPT; the snapshots are made again. Another supercell (a FORCEDISP2) is a separate project.</p>
  <div class="grid"><div><label>New name</label><input id="c_name"></div><div><label>Temperature</label><input id="c_temperature"></div>
  <div><label>Supercell</label><input id="c_supercell"></div><div><label>Snapshots</label><input id="c_snapshots"></div></div>
  <button class="btn" onclick="cloneProj()">Create</button><span id="c_msg" class="muted"></span>
 </div>
</section>

<section class="pane" id="p-run">
 <div class="card"><div class="ct">Run the next step</div>
  <p id="runinfo" class="muted"></p>
  <button class="btn p" onclick="act('continue')" title="runs the next steps that are quick here (POSCAR0, snapshots, folders, collect after import), stops at DFT runs and before FOCEX">Continue</button>
  <button class="btn" onclick="act('run_here')" title="the next DFT stage runs on this computer">Run here</button>
  <button class="btn" onclick="act('bundle')" title="the next DFT stage as a .tar.gz with SLURM scripts">Make bundle</button>
  <button class="btn" id="b_submit" onclick="act('submit')" title="sbatch on this machine (the GUI runs on the cluster)">Submit here</button>
  <button class="btn" id="b_skip" onclick="skipRelax()">Skip: POSCAR is relaxed</button>
  <div style="margin-top:8px"><label>Import results (results_*.tar.gz made by pack_results.sh on the cluster)</label>
   <div class="row"><input type="file" id="imp" style="flex:1"><button class="btn" onclick="imp()">Import results</button></div></div>
  <div id="bundles" style="margin-top:8px"></div>
 </div>
 <div class="card"><div class="ct">Cluster (SLURM) settings for bundles and "Submit here"</div>
  <div class="grid" id="hpcgrid"></div>
  <div style="margin-top:8px"><label>Lines run before VASP / THERMACOND in every job (module load ...)</label><textarea id="h_modules" rows="2"></textarea></div>
  <div style="margin-top:8px"><label>More #SBATCH lines</label><textarea id="h_extra_sbatch" rows="2"></textarea></div>
  <button class="btn" onclick="saveHpc()">Save settings</button><span id="h_msg" class="muted"></span>
  <p class="muted">The POTCAR is built here (VASP_POTCAR_DIR) and goes into the bundle. VASP ranks for KPAR: nodes × tasks per node for bundles, "MPI ranks here" for runs here.</p>
 </div>
</section>

<section class="pane" id="p-snaps">
 <div class="card"><div class="ct">Snapshots (SC-SNAPS)</div><div id="snapinfo" class="muted">Not made yet.</div>
  <div class="row" id="hists"></div></div>
</section>

<section class="pane" id="p-params">
 <div class="card"><div class="ct">FOCEX input files</div>
  <p class="muted">Generated from the structure; every line has a comment after the values. Lines marked PLACEHOLDER need your input or are filled in by a later step. Changed files are kept in the instruction file and used as they are.</p>
  <div id="pnotes"></div>
  <div class="tabs ptabs" id="ptabs"></div>
  <textarea id="ptext" rows="24" spellcheck="false"></textarea>
  <div><span class="pill" id="psrc"></span>
  <button class="btn" onclick="resetParam()">Back to the generated file</button>
  <button class="btn p" onclick="saveParams()">Save all</button>
  <button class="btn p" onclick="saveParams(true)">Save and run FOCEX</button><span id="p_msg" class="muted"></span></div>
 </div>
</section>

<section class="pane" id="p-results">
 <div class="card"><div class="ct">FOCEX and THERMACOND results</div><div id="pdflist" class="muted">No plots yet.</div><iframe id="pdfview" style="display:none"></iframe></div>
</section>

<section class="pane" id="p-tc">
 <div class="card"><div class="ct">THERMACOND: latdyn.params of 07_thermacond</div>
  <p class="muted">Line 8 chooses the job (0 none, 1 linewidths, 2 kappa direct, 3 kappa iterative); line 1 the q mesh; line 6 the temperatures.</p>
  <textarea id="tctext" rows="16" spellcheck="false"></textarea><span class="pill" id="tcsrc"></span>
  <div><button class="btn" onclick="saveTc()">Save</button>
  <button class="btn p" onclick="saveTc('tc_here')">Save and run here</button>
  <button class="btn" onclick="saveTc('tc_bundle')">Save and make bundle</button>
  <button class="btn" id="b_tcsub" onclick="saveTc('tc_submit')">Save and submit here</button><span id="tc_msg" class="muted"></span></div>
 </div>
</section>

<section class="pane" id="p-setup">
 <div class="card"><div class="ct">Instruction file (advanced)</div>
  <p class="muted">All settings of the project; the keywords are explained in the focex-agent manual.</p>
  <textarea id="instr" rows="26" spellcheck="false"></textarea>
  <button class="btn p" onclick="saveInstr()">Save</button><span id="i_msg" class="muted"></span></div>
</section>

<section class="pane" id="p-log">
 <div class="card"><div class="ct">Job output</div><div class="term" id="joblog">—</div></div>
 <div class="card"><div class="ct">report.log (one entry per task)</div><div class="term" id="report">—</div></div>
</section>
</main></div>
<script>
const TABS=[['projects','Projects'],['run','Run'],['snaps','Snapshots'],['params','Params'],['results','Results'],['tc','THERMACOND'],['setup','Setup'],['log','Log']];
const HPCF=[['partition','Partition'],['account','Account'],['nodes','Nodes'],['ntasks_per_node','Tasks per node'],['time','Walltime (VASP)'],
 ['time_thermacond','Walltime (THERMACOND)'],['mem','Memory (optional)'],['mpi_cmd','MPI launcher (srun or mpirun)'],['vasp_exec','VASP binary on the cluster'],
 ['tcond_exec','THERMACOND binary on the cluster'],['python','python3 on the cluster']];
let P=null,S=null,PAR=null,cur='structure_params',timer=null;
const $=id=>document.getElementById(id);
function tabs(){ $('tabs').innerHTML=TABS.map(([k,t])=>`<div class="tab" id="t-${k}" onclick="show('${k}')">${t}</div>`).join('');
 $('hpcgrid').innerHTML=HPCF.map(([k,t])=>`<div><label>${t}</label><input id="h_${k}"></div>`).join('')+`<div><label>MPI ranks here (this computer)</label><input id="h_local_mpi"></div>`; }
function show(k){document.querySelectorAll('.pane').forEach(e=>e.classList.remove('on'));document.querySelectorAll('.tab').forEach(e=>e.classList.remove('on'));
 $('p-'+k).classList.add('on');$('t-'+k).classList.add('on'); if(k==='params'||k==='tc') loadParams();}
async function api(u,o){const r=await fetch(u,o);let j;try{j=await r.json()}catch(e){j={ok:false,error:'HTTP '+r.status}} return j;}
const post=(u,b)=>api(u,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(b)});
async function projects(){const j=await api('/api/projects');$('plist').innerHTML=j.projects.map(p=>`<option value="${p.path}">${p.name}</option>`).join('');}
function readPos(){const f=$('n_file').files[0];if(!f)return;f.text().then(t=>{$('n_poscar').value=t;if(!$('n_name').value)$('n_name').value=f.name.replace(/\.(vasp|poscar)$/i,'');});}
async function newProj(){$('n_msg').textContent=' setting up (POTCAR, k-mesh, relaxation inputs) ...';
 const b={name:$('n_name').value,poscar:$('n_poscar').value,functional:$('n_functional').value,accuracy:$('n_accuracy').value,dfpt:$('n_dfpt').value,
 ranks:$('n_ranks').value,range:$('n_range').value,supercell:$('n_supercell').value,snapshots:$('n_snapshots').value,temperature:$('n_temperature').value,
 tc_job:$('n_tc').value,local_mpi:$('n_mpi').value,relaxed:$('n_relaxed').checked,extra:$('n_extra').value};
 const j=await post('/api/new_project',b); if(!j.ok){$('n_msg').textContent=' '+j.error;return;}
 $('n_msg').textContent=' created '+j.project; await projects(); P=j.project; refresh(); show('run');}
async function cloneProj(){if(!P)return;const j=await post('/api/clone',{project:P,name:$('c_name').value,temperature:$('c_temperature').value,supercell:$('c_supercell').value,snapshots:$('c_snapshots').value});
 $('c_msg').textContent=j.ok?' created '+j.project:' '+j.error; if(j.ok){await projects();P=j.project;refresh();}}
function openProj(){P=$('plist').value;PAR=null;refresh();show('run');}
async function refresh(){ if(!P) return; const j=await api('/api/status?project='+encodeURIComponent(P)); if(!j.ok){$('advice').innerHTML='<p>'+j.error+'</p>';return;}
 S=j; $('hdrproj').textContent='project: '+j.project;
 $('steps').innerHTML=j.guide.steps.map((s,i)=>`<div class="step ${s.status}"><div class="dot">${s.status==='done'?'✓':i+1}</div><div><b>${s.title}</b> <span class="muted">${s.status}</span><small>${s.what}${s.detail?' — '+s.detail:''}</small></div></div>`).join('');
 $('advice').innerHTML=j.guide.advice.map(a=>`<p>${esc(a)}</p>`).join('');
 $('joblog').textContent=(j.job.lines||[]).join('\n')+(j.job.result?'\n=> '+j.job.result:'')||'—'; $('report').textContent=j.report.join('\n');
 const t=$('report');t.scrollTop=t.scrollHeight; const jl=$('joblog'); jl.scrollTop=jl.scrollHeight;
 $('b_submit').disabled=!j.sbatch; $('b_tcsub').disabled=!j.sbatch; $('b_skip').disabled=!!(j.state.done||{}).relax;
 const w=j.state.waiting; $('runinfo').textContent=(j.job.running?'Running: '+j.job.kind+' ... ':'')+(w?`Waiting: ${w.task}, ${w.dirs.length} folder(s) ready (${w.dirs.slice(0,3).join(', ')}${w.dirs.length>3?' ...':''}).`:'')+(j.potcar?'':' No POTCAR in the project: put one there.');
 $('bundles').innerHTML=j.bundles.length?'Bundles: '+j.bundles.map(b=>`<a href="/api/file?download=1&project=${encodeURIComponent(P)}&path=${encodeURIComponent(b)}">${b.split('/').pop()}</a>`).join(' · '):'';
 for(const [k] of HPCF){ if(document.activeElement!==$('h_'+k)) $('h_'+k).value=j.hpc[k]??''; }
 if(document.activeElement!==$('h_modules')) $('h_modules').value=j.hpc.modules||''; if(document.activeElement!==$('h_extra_sbatch')) $('h_extra_sbatch').value=j.hpc.extra_sbatch||'';
 if(document.activeElement!==$('h_local_mpi')) $('h_local_mpi').value=j.cfg.local_mpi;
 if(document.activeElement!==$('instr')) $('instr').value=j.instr;
 $('pdflist').innerHTML=j.pdfs.length?j.pdfs.map(f=>`<a href="#" onclick="pdf('${f}');return false">${f}</a>`).join(' · '):'No plots yet.';
 const sn=j.report.filter(l=>l.includes('=== snaps -- OK')); $('snapinfo').textContent=sn.length?sn[sn.length-1].split('OK')[1]:'Not made yet.';
 $('hists').innerHTML=j.snaps?['disp','comp'].map(k=>`<img style="max-width:48%" src="/api/histogram/${k}?project=${encodeURIComponent(P)}&t=${Date.now()}">`).join(''):'';
 clearTimeout(timer); timer=setTimeout(refresh, j.job.running?2500:15000); }
function esc(s){return s.replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]));}
function pdf(f){const v=$('pdfview');v.style.display='block';v.src='/api/file?project='+encodeURIComponent(P)+'&path='+encodeURIComponent(f);}
async function act(a){if(!P)return;const j=await post('/api/action',{project:P,action:a});if(!j.ok)alert(j.error);else if(j.bundle)alert('bundle: '+j.bundle);else if(j.output)alert(j.output);show('log');refresh();}
async function skipRelax(){if(!confirm('Use the POSCAR as the relaxed structure (no relaxation)?'))return;const j=await post('/api/skip_relax',{project:P});if(!j.ok)alert(j.error);refresh();}
async function imp(){const f=$('imp').files[0];if(!f||!P)return;const fd=new FormData();fd.append('project',P);fd.append('file',f);
 const j=await api('/api/import',{method:'POST',body:fd}); alert(j.ok?`imported ${j.files.length} files. ${j.message||''}\nNow press Continue.`:j.error); refresh();}
async function saveHpc(){const h={};for(const [k] of HPCF)h[k]=$('h_'+k).value;h.modules=$('h_modules').value;h.extra_sbatch=$('h_extra_sbatch').value;
 const j=await post('/api/hpc',{project:P,hpc:h,local_mpi:$('h_local_mpi').value});$('h_msg').textContent=j.ok?' saved':' '+j.error;}
async function saveInstr(){const j=await post('/api/instr',{project:P,text:$('instr').value});$('i_msg').textContent=j.ok?' saved':' '+j.error;}
async function loadParams(){if(!P)return;const j=await api('/api/params?project='+encodeURIComponent(P));if(!j.ok){$('pnotes').innerHTML='<p class="note">'+esc(j.error)+'</p>';return;}
 PAR=j.params; for(const k in PAR) PAR[k].edit=PAR[k].text; $('pnotes').innerHTML=j.notes.map(n=>`<p class="note">${esc(n)}</p>`).join('');
 $('ptabs').innerHTML=Object.keys(PAR).map(k=>`<div class="tab" id="pt-${k}" onclick="pick('${k}')">${k.replace('_params','.params')}</div>`).join(''); pick(cur);
 $('tctext').value=j.tc.text; $('tcsrc').textContent=j.tc.source;}
function pick(k){if(PAR&&PAR[cur])PAR[cur].edit=$('ptext').value;cur=k;document.querySelectorAll('.ptabs .tab').forEach(e=>e.classList.remove('on'));$('pt-'+k).classList.add('on');$('ptext').value=PAR[k].edit;$('psrc').textContent=PAR[k].source;}
function resetParam(){$('ptext').value=PAR[cur].generated;}
async function saveParams(run){PAR[cur].edit=$('ptext').value;const p={};for(const k in PAR)p[k]=PAR[k].edit;const j=await post('/api/save_params',{project:P,params:p});
 $('p_msg').textContent=j.ok?(' saved'+(j.saved.length?' (yours: '+j.saved.join(', ')+')':' (all as generated)')):' '+j.error; if(j.ok&&run){act('focex');} else loadParams();}
async function saveTc(a){const j=await post('/api/save_params',{project:P,tc:$('tctext').value});$('tc_msg').textContent=j.ok?' saved':' '+j.error;if(j.ok&&a)act(a);}
tabs();projects();show('projects');
</script></body></html>
"""


def main():
    global WORKDIR
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--workdir', default=os.getcwd(), help='folder of the projects (default: here)')
    ap.add_argument('--port', type=int, default=5070)
    ap.add_argument('--host', default='127.0.0.1', help='127.0.0.1: this computer only (on a cluster use an ssh tunnel)')
    a = ap.parse_args()
    WORKDIR = os.path.realpath(os.path.expanduser(a.workdir))
    print(f'aladyn-gui {VERSION}: projects in {WORKDIR}; open http://localhost:{a.port}')
    app.run(host=a.host, port=a.port, debug=False, threaded=True)


if __name__ == '__main__':
    main()
