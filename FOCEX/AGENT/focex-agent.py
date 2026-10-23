#!/usr/bin/env python3
"""focex-agent.py (FOCEX/AGENT) -- from a primitive POSCAR to FOCEX phonons, with VASP and SC_SNAPS.
Manual: FOCEX/AGENT/focexagent_manual.pdf

    focex-agent.py init [-i instructions.txt] [-s POSCAR] [--profile rivanna]
    focex-agent.py init -i instructions.txt      # with DATA: <folder> -- FOCEX on existing data
    focex-agent.py run  <project>              # every task, in order, restartable
    focex-agent.py run  <project> --retry-failed
    focex-agent.py run  <project> --only snaps  # one task
    focex-agent.py plots <project>             # remake the plot scripts and plots

Tasks (each is recorded in <project>/report.log and skipped once done):
  1 relax      Niggli-reduced, oriented cell; VASP relaxation (ACCURACY: Accurate|Normal)
  2 poscar0    CONTCAR -> Niggli + orientation -> POSCAR0 (the frame of every later step)
    scf        SCF on POSCAR0, same INCAR (ISIF=1, NSW=0) and KPOINTS
    dfpt       epsilon_inf and Born charges, only if the SCF shows a gap (DFPT: auto|on|off)
  3 snaps      SC_SNAPS: nearly cubic supercell (SUPERCELL: SMALL|LARGE|9 integers),
               snapshots; the temperature is tuned to an rms displacement of DISP_TARGET
  4 supercells one folder per snapshot (INCAR as the SCF + ISYM=0, commensurate KPOINTS)
    vasp_snaps the supercell VASP runs, one after the other
  5 collect    read_outcar.x -> FORCEDISP1, POSCAR1 = poscar_000; all runs checked
  6 focex      the five *.params files, FOCEX (+ NOBORN and LASSO runs when needed)
  7 plots      bvel.plt / all.plt / qha.plt filled in, makeplots.sh

Existing data (no VASP, no SC_SNAPS): DATA: <folder> with POSCAR1, FORCEDISP1 (POSCAR2,
FORCEDISP2, ...), optionally POSCAR (the primitive cell, same frame as the data; otherwise the
primitive cell of POSCAR1 found by spglib) and dielectric.params (epsilon and the Born charges
in the order of the atoms of that POSCAR). SPECIES: Mg O names the species when the POSCARs have
no element line. Everything is rotated into the FOCEX frame (POSCARi lattice and positions,
FORCEDISPi positions and forces, epsilon and Z as R^T T R); tasks relax ... collect are marked
"provided" and run goes straight to focex and plots.

Every generated input can be replaced by a block in the instructions file, named after the
file:  INCAR ... INCAR_END  (also INCAR_relax, INCAR_scf, INCAR_dfpt, INCAR_snaps),
KPOINTS, KPOINTS_snaps, cell_inp, supercell_inp, snaps_inp, structure_params,
latdyn_params, kpbs_params, dielectric_params, default_params  (END_<NAME> also accepted,
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

AGENT = os.path.dirname(os.path.abspath(__file__))          # FOCEX/AGENT: this script, its manual
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
    from vasp_input_generator import VASPInputGenerator, potcar_dir_for
    vagent = _load(os.path.join(VASPFLOW, 'vasp-agent.py'), 'vasp_agent')
except Exception as e:                                        # pragma: no cover
    sys.exit(f'focex-agent needs VASP-Flow at {VASPFLOW} (set VASP_FLOW_DIR): {e}')

ACCURACY = {'accurate': dict(ENCUT=600, EDIFF='1E-7', EDIFFG=-0.001, PREC='Accurate'),
            'normal':   dict(ENCUT=450, EDIFF='1E-6', EDIFFG=-0.003, PREC='Normal')}
SUPERCELL = {'small': (10.0, 100), 'large': (14.0, 300)}     # min face distance (A), max atoms
OMEGA2_CM = 521.11 ** 2          # bands.dat: omega^2 in eV/(A^2 amu) -> (1/cm)^2
METALS = set('''Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Tc Ru Rh
 Pd Ag Cd In Sn Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg
 Tl Pb Bi Po Fr Ra Ac Th Pa U Np Pu Am Cm'''.split())
BLOCKS = ('incar', 'incar_relax', 'incar_scf', 'incar_dfpt', 'incar_snaps', 'kpoints',
          'kpoints_snaps', 'cell_inp', 'supercell_inp', 'snaps_inp', 'structure_params',
          'latdyn_params', 'kpbs_params', 'dielectric_params', 'default_params')
TASKS = ('relax', 'poscar0', 'scf', 'dfpt', 'snaps', 'supercells', 'vasp_snaps', 'collect',
         'focex', 'plots')


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
    fx['species'] = (keyword(rest, 'SPECIES') or '').split()
    fx['sc_snaps_x'] = os.path.expanduser(keyword(rest, 'SC_SNAPS_X') or '~/BIN/sc_snaps.x')
    fx['read_outcar_x'] = os.path.expanduser(keyword(rest, 'READ_OUTCAR_X') or '~/ALADYN/bin/read_outcar.x')
    fx['profile'] = keyword(rest, 'PROFILE')
    fx['blocks'] = blocks
    return inst, fx


def newest_focex():
    c = sorted(glob.glob(os.path.expanduser('~/BIN/v[0-9]*')), key=os.path.getmtime)
    return c[-1] if c else '~/BIN/v8.18'


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
    # VASP-Flow applies 'all' and then the step's list, later lines winning; scf and snaps
    # INCARs are made from the relaxation INCAR, so their list goes under 'relax'
    d['incar_raw'] = {'all': first + user_all,
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
    face distance within the limit.  (Skewed Hermite-normal-form supercells of the same size
    can put FC2 pairs outside the supercell's Wigner-Seitz cell, which stops FOCEX in
    set_dynmat_std.)  Returns (S rows in primitive units, info)."""
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
def outcar_status(d, kind):
    """(ok, message) for the OUTCAR in d; kind: relax | static | dfpt."""
    f = os.path.join(d, 'OUTCAR')
    if not os.path.isfile(f):
        return False, 'no OUTCAR (not run)'
    t = open(f, errors='replace').read()
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


def per_band(text, nb):
    """The plot statement after a 'THIS HAS AS MANY LINES AS THERE ARE BANDS' comment gets one
    entry per band: in its first multi-entry file, column 6 becomes 6+b and the other column
    6+nb+b (b = 0..nb-1), e.g. 'ibz_grun.dat' u 12:6 -> u (6+nb+b):(6+b)."""
    L = text.splitlines()
    out, i = [], 0
    while i < len(L):
        out.append(L[i])
        if 'AS MANY LINES AS THERE ARE BANDS' not in L[i]:
            i += 1
            continue
        i += 1
        while i < len(L) and not re.match(r'^\s*plot\b', L[i]):
            out.append(L[i])
            i += 1
        stmt = []
        while i < len(L):
            stmt.append(L[i].rstrip().rstrip('\\').rstrip().rstrip(','))
            i += 1
            if not L[i - 1].rstrip().endswith('\\'):
                break
        stmt[0] = re.sub(r'^\s*plot\s+', '', stmt[0])
        files = [re.match(r"\s*'([^']+)'", e).group(1) if re.match(r"\s*'([^']+)'", e) else None for e in stmt]
        rep = next(f for f in files if f and files.count(f) > 1)
        tpl = stmt[files.index(rep)]
        head = re.match(r"(\s*'[^']+'\s+u(?:sing)?\s+)(.*)", tpl)
        new = [e for e, f in zip(stmt, files) if f != rep]
        for b in range(nb):
            cols = re.sub(r'(?<![\w.])\$?(\d+)(?![\w.])',
                          lambda m: m.group(0).replace(m.group(1), str(6 + b if m.group(1) == '6' else 6 + nb + b)),
                          head.group(2).split()[0])
            new.append(head.group(1) + cols + head.group(2)[len(head.group(2).split()[0]):])
        out.append('plot ' + ' ,\\\n     '.join(e.strip() for e in new))
    return '\n'.join(out) + '\n'


# ── the agent ─────────────────────────────────────────────────────────────────
class Agent:
    def __init__(self, proj):
        self.proj = os.path.abspath(proj)
        self.rep = Report(self.proj)
        st = os.path.join(self.proj, 'state.json')
        self.state = json.load(open(st)) if os.path.isfile(st) else {}
        self.inst, self.fx = parse_instructions(os.path.join(self.proj, 'instructions.txt'))
        self.profile = vagent.load_profile(self.fx['profile'] or self.state.get('profile') or '')

    def save(self):
        write(os.path.join(self.proj, 'state.json'), json.dumps(self.state, indent=1, default=str))

    def d(self, *a):
        return os.path.join(self.proj, *a)

    def gen(self, poscar, step, extra=()):
        return VASPInputGenerator(poscar, vasp_inst(self.inst, self.fx, step, extra), profile=self.profile)

    def vasp(self, d, gen):
        """Run VASP in d; an existing OUTCAR is kept as OUTCAR.1, OUTCAR.2, ... (never overwritten)."""
        if os.path.isfile(os.path.join(d, 'OUTCAR')):
            k = 1
            while os.path.exists(os.path.join(d, f'OUTCAR.{k}')):
                k += 1
            shutil.move(os.path.join(d, 'OUTCAR'), os.path.join(d, f'OUTCAR.{k}'))
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
            return True, msg
        reuse = self.state.get('reuse_relax')
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
        dt = self.vasp(d, g)
        ok, msg = outcar_status(d, 'static')
        if ok:
            self.clean(self.d('01_relax'))
            gap = band_gap(d)
            self.state['gap'] = gap
            msg += f'; band gap {gap:.3f} eV' if gap is not None else '; band gap unreadable'
            species = read_vasp_cell(self.d('POSCAR0'))['species']
            if gap is not None and gap < 0.05 and not self.is_metal(species) and not self.fx['blocks'].get('kpoints'):
                msg += ' -- WARNING: metallic, but the coarse k-mesh was used (set KMESH_DENSITY: fine)'
        return ok, msg

    def want_dfpt(self):
        m = self.fx['dfpt']
        return m == 'on' or (m == 'auto' and (self.state.get('gap') or 0) > 0.05)

    def t_dfpt(self):
        d = self.d('02_scf', 'dfpt')
        if not self.want_dfpt():
            self.clean(self.d('02_scf'))
            return True, 'skipped (no gap)' if self.fx['dfpt'] == 'auto' else 'skipped (DFPT: off)'
        ok, msg = outcar_status(d, 'dfpt')
        if not ok:
            g = self.gen(self.d('POSCAR0'), 'dfpt')
            g.generate_dfpt_input(d, self.d('02_scf'))
            self.link_potcar(d)
            subprocess.run('bash copy_from_scf.sh', shell=True, cwd=d)
            self.vasp(d, g)
            ok, msg = outcar_status(d, 'dfpt')
        if ok:
            n = sum(read_vasp_cell(self.d('POSCAR0'))['counts'])
            eps, Z = read_born(os.path.join(d, 'OUTCAR'), n, self.fx['born_transpose'])
            if eps is None:
                return False, 'converged, but epsilon/Born charges could not be read from OUTCAR'
            self.state['eps'] = eps.tolist()
            self.state['Z'] = [z.tolist() for z in Z]
            self.clean(d)
            self.clean(self.d('02_scf'))
            msg += '; eps_inf diag ' + ' '.join('%.3f' % x for x in np.diag(eps))
        return ok, msg

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
            st = self.disp_stats(d)
            lines.append(f"T = {T:g} K: <u> = {' '.join('%+.4f' % x for x in st['mean'])} A, "
                         f"sigma = {' '.join('%.4f' % x for x in st['std'])} A (x y z), rms {st['rms']:.4f} A, "
                         f"mean |u| = {st['norm']:.4f} A; per element sigma: "
                         + ', '.join(f'{k} {v:.4f}' for k, v in st['elem'].items()))
            if b.get('snaps_inp') or abs(st['rms'] / self.fx['disp_target'] - 1) < 0.15:
                break
            T = T * (self.fx['disp_target'] / st['rms']) ** 2
            lines.append(f'  rms differs from DISP_TARGET {self.fx["disp_target"]} A: rerun at T = {T:.1f} K')
        self.state['T_snaps'] = T
        ok, msg = self.check_commensurate(d, c, S)
        lines.append(msg)
        return ok, '\n  '.join(lines)

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
    def t_supercells(self):
        src, top = self.d('03_snaps'), self.d('04_supercells')
        os.makedirs(top, exist_ok=True)
        snaps = snap_files(src)
        S = np.array(self.state['S'])
        V = read_vasp_cell(self.d('POSCAR0'))['A']
        kind, Np = VASPInputGenerator._parse_kpoints_text(open(self.d('02_scf', 'KPOINTS')).read())
        if kind != 'mesh' and not self.fx['blocks'].get('kpoints_snaps'):
            return False, '02_scf/KPOINTS is not an automatic mesh: give a KPOINTS_snaps block'
        mesh = sc_kmesh(S, V, Np)
        nk = max(1, (int(np.prod(mesh)) + 1) // 2)              # ISYM=0, time reversal only
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
        for f in snaps:
            write(os.path.join(top, 'snap_' + f.rsplit('_', 1)[1], 'INCAR'), incar)
        return True, (f'{len(snaps)} folders; KPOINTS {mesh} from the primitive {Np}; '
                      f'{sum(base["counts"])} atoms' + (f'; KPAR {inst_k["relax_kpar"]} NCORE {inst_k["relax_ncore"]}'
                                                       if inst_k else ''))

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
        dirs = sorted(glob.glob(self.d('04_supercells', 'snap_*')))
        g = self.gen(os.path.join(dirs[0], 'POSCAR'), 'snaps')
        failed, ran = [], 0
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

    def write_resubmit(self, failed):
        g = self.gen(self.d('POSCAR.init'), 'relax')
        body = (g._run_sh_preamble('focex_retry') +
                f'cd "{self.proj}"\n# failed: {" ".join(os.path.basename(d) for d, _ in failed) or "none"}\n'
                f'{sys.executable} "{os.path.abspath(__file__)}" run "{self.proj}" --retry-failed\n')
        write(self.d('resubmit_failed.sh'), '#!/bin/bash\n' + body, 0o755)

    # ---- task 5: collect -----------------------------------------------------
    def t_collect(self):
        top = self.d('04_supercells')
        dirs = sorted(glob.glob(os.path.join(top, 'snap_*')))
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
        if worst_sum > 0.05 or worst_pos > 1e-3:
            msg += ' -- WARNING: check these'
        return True, msg

    # ---- task 6: FOCEX -------------------------------------------------------
    def structure_params(self, line5='1 1 1 0'):
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
             line5 + '   # translational, rotational, Huang invariances; enforce_inv (2 = LASSO)',
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
        m = re.findall(r'Percent error, \|\| F_dft-F_fit \|\| / \|\| F_dft \|\| =\s*([-\d.Ee+]+)', t)
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
            runs['NOBORN'] = self.run_focex(os.path.join(top, 'NOBORN'), 0, '1 1 1 0', 'NOBORN (born_flag=0)')
        runs['main'] = main = self.run_focex(top, born, '1 1 1 0', f'main (born_flag={born})')
        if main['ok'] and main['eq'] and main['unk'] and main['eq'] / main['unk'] < self.fx['lasso_ratio']:
            self.rep(f"  equations/unknowns = {main['eq'] / main['unk']:.2f} < {self.fx['lasso_ratio']}: LASSO run")
            runs['LASSO'] = self.run_focex(os.path.join(top, 'LASSO'), born, '1 1 1 2', 'LASSO (enforce_inv=2)')
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
        if not os.path.isfile(os.path.join(d, 'KTICS.BS')):
            return False, 'no KTICS.BS (FOCEX did not finish)'
        ktics = open(os.path.join(d, 'KTICS.BS')).read().strip()
        kmax = float(re.findall(r'([-\d.]+)\s*\)', ktics)[-1])
        bs = np.loadtxt(os.path.join(d, 'bs_freq.dat'), comments='#', ndmin=2)
        sub = {'MAX_K_VALUE': f'{kmax:.3f}',
               'MAX_FREQUENCY_IN_INVERSE_CM': str(math.ceil(1.05 * bs[:, 6].max() / 10) * 10),
               'MAX_GROUP_VELOCITY': f'{math.ceil(bs[:, 10].max() * 10) / 10:.1f}',
               'NAME': self.fx['project']}
        unknown = set()
        nb = 3 * int(open(os.path.join(d, 'structure.params')).read().splitlines()[10].split()[0])
        for name in ('bvel.plt', 'all.plt', 'qha.plt'):
            out, replace_next = [], False
            for line in per_band(open(os.path.join(HERE, name)).read(), nb).splitlines():
                if replace_next and 'xtics' in line:
                    out.append(ktics)
                    replace_next = False
                    continue
                if 'REPLACE THE FOLLOWING LINE BY THE CONTENTS OF "KTICS.BS"' in line:
                    replace_next = True
                if re.match(r'^\s*pause\b', line):
                    continue
                for k, v in sub.items():
                    line = line.replace(k, v)
                unknown |= set(re.findall(r'\b[A-Z][A-Z0-9]*(?:_[A-Z0-9]+){2,}\b', line.split('#')[0]))
                out.append(line)
            text = '\n'.join(out) + '\n'
            if re.search(r'^\s*set\s+term', text, re.M) and not re.match(r'\s*set\s+term', text):
                text = 'set term unknown\n' + text                  # plots before the pdf terminal
            write(os.path.join(d, name), text)
        write(os.path.join(d, 'makeplots.sh'),
              '#!/bin/bash\ncd "$(dirname "$0")"\nst=0\nfor f in bvel.plt all.plt qha.plt; do\n'
              '  gnuplot "$f" > "${f%.plt}.log" 2>&1 && echo "  $f: ok" || { echo "  $f: ERROR, see ${f%.plt}.log"; st=1; }\n'
              'done\nexit $st\n', 0o755)
        if shutil.which('gnuplot') is None:
            return True, 'scripts written; gnuplot not found, run makeplots.sh where it is'
        r = subprocess.run(['bash', 'makeplots.sh'], cwd=d, capture_output=True, text=True)
        m = r.stdout.strip().replace('\n', ';')
        if unknown:
            m += f'; unfilled placeholders: {" ".join(sorted(unknown))}'
        return r.returncode == 0, m

    # ---- driver --------------------------------------------------------------
    def run(self, only=None, retry=False):
        todo = [only] if only else TASKS
        for t in todo:
            done = self.state.setdefault('done', {})
            if t in done and not only and not (retry and t in ('vasp_snaps',)):
                continue
            self.rep(f'=== {t} -- started')
            t0 = time.time()
            fn = getattr(self, 't_' + t)
            ok, msg = fn(retry) if t == 'vasp_snaps' else fn()
            self.rep(f'=== {t} -- {"OK" if ok else "FAILED"} ({(time.time() - t0) / 60:.1f} min): {msg}')
            if not ok:
                self.state['failed'] = t
                self.save()
                return False
            done[t] = now()
            self.state.pop('failed', None)
            self.save()
        return True


# ── init ──────────────────────────────────────────────────────────────────────
def init(args):
    inst, fx = parse_instructions(args.instructions)
    proj = os.path.abspath(vagent.slugify(fx['project']))
    if os.path.exists(os.path.join(proj, 'state.json')) and not args.force:
        sys.exit(f'{proj} exists; use --force to start it again')
    os.makedirs(proj, exist_ok=True)
    shutil.copy(args.instructions, os.path.join(proj, 'instructions.txt'))
    if fx['data']:
        return init_data(args, fx, proj)
    shutil.copy(args.structure, os.path.join(proj, 'POSCAR.input'))
    p = read_vasp_cell(args.structure)
    c = canonical_cell(p)
    write_poscar(os.path.join(proj, 'POSCAR.init'), f"{fx['project']} (Niggli-reduced, FOCEX frame)",
                 c['V'], p['species'], p['counts'], c['F'])
    state = {'profile': args.profile or fx['profile'] or '', 'created': now()}
    write(os.path.join(proj, 'state.json'), json.dumps(state, indent=1))
    metal = all(s.split('_')[0] in METALS for s in p['species'])
    # relax and SCF share the mesh: coarse (1000 kpra), fine for metals/intermetallics
    dens = 'fine' if metal else 'coarse'
    with open(os.path.join(proj, 'instructions.txt'), 'a') as fh:
        if not re.search(r'^\s*(RELAX_)?KMESH(_DENSITY)?\s*:', open(args.instructions).read(), re.M | re.I):
            fh.write(f'\n# added by focex-agent init: relax and SCF k-mesh ({"metal" if metal else "not a metal"})\n'
                     f'RELAX_KMESH_DENSITY: {dens}\nKMESH_DENSITY: {dens}\n')
    a = Agent(proj)
    rep = a.rep
    rep(f"focex-agent init: {fx['project']}, {len(p['F'])} atoms, space group {c['sg']}, accuracy "
        f"{fx['accuracy']}, relax/SCF k-mesh {dens}")
    potcar_dir, _ = potcar_dir_for(inst.get('functional', 'R2SCAN'))
    if potcar_dir and os.path.isdir(potcar_dir):
        vagent.build_potcar(p['species'], potcar_dir, os.path.join(proj, 'POTCAR'))
    else:
        rep('WARNING: no POTCAR library (VASP_POTCAR_DIR); put a POTCAR in the project folder before running')
    for x in ('focex_x', 'sc_snaps_x', 'read_outcar_x'):
        if not os.path.isfile(fx[x]):
            rep(f'WARNING: {x.upper()} = {fx[x]} not found')
    a.write_relax()
    g = a.gen(os.path.join(proj, 'POSCAR.init'), 'relax')
    write(os.path.join(proj, 'run_workflow.sh'), '#!/bin/bash\n' + g._run_sh_preamble('focex_agent') +
          f'cd "{proj}"\n{sys.executable} "{os.path.abspath(__file__)}" run "{proj}"\n', 0o755)
    a.write_resubmit([])
    rep(f'wrote {proj}: POSCAR.init, 01_relax/, run_workflow.sh  ->  cd {os.path.basename(proj)} && '
        + ('sbatch' if a.profile.get('slurm') else 'bash') + ' run_workflow.sh')


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
    state = {'profile': args.profile or fx['profile'] or '', 'created': now(), 'npairs': n,
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
    a.rep(f"then:  {sys.executable} {os.path.abspath(__file__)} run {proj}   (tasks focex, plots)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    a = sp.add_parser('init')
    a.add_argument('-i', '--instructions', default='instructions.txt')
    a.add_argument('-s', '--structure', default='POSCAR')
    a.add_argument('--profile', default=None, help='VASP-Flow profile (workstation, rivanna, slurm)')
    a.add_argument('--force', action='store_true')
    r = sp.add_parser('run')
    r.add_argument('project')
    r.add_argument('--only', choices=TASKS)
    r.add_argument('--retry-failed', action='store_true')
    q = sp.add_parser('plots')
    q.add_argument('project')
    args = ap.parse_args()
    if args.cmd == 'init':
        init(args)
    elif args.cmd == 'run':
        sys.exit(0 if Agent(args.project).run(args.only, args.retry_failed) else 1)
    else:
        sys.exit(0 if Agent(args.project).run('plots') else 1)


if __name__ == '__main__':
    main()
