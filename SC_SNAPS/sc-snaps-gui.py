#!/usr/bin/env python3
"""
SC-Snaps GUI — browser interface for sc_snaps.x supercell snapshot generator
Run with:  python3 sc-snaps-gui.py   (or:  make run)
Opens http://localhost:5050

The working directory is wherever you invoke this script from.
"""
import os, sys, glob, json, queue, threading, re
from pathlib import Path
from flask import Flask, request, jsonify, Response

app = Flask(__name__)

# ── fixed working directory (set at launch time) ──────────────────────────────
WORKDIR = os.getcwd()


def _request_workdir(payload=None):
    """Return the user-selected working directory or the launch directory."""
    payload = payload or {}
    wd = payload.get('workdir') if isinstance(payload, dict) else None
    if wd is None:
        wd = request.args.get('workdir', '')
    wd = os.path.expanduser(str(wd).strip()) if wd else WORKDIR
    return os.path.abspath(wd)

# ── defaults ─────────────────────────────────────────────────────────────────
DEFAULT_EXEC = '~/BIN/sc_snaps.x'

DEFAULT_CELL = """\
1 1 1   90 90 90
 0 0.5 0.5, 0.5 0 0.5, 0.5 0.5 0
4.247  9  # scale, 0 for conventional reduced units
2         # of atoms types
1 1       # of atoms of each type
24.31  16.00
Mg O
  0 0 0
  0.5 0.5 0.5"""

DEFAULT_SNAPS = """\
400  Avg frequency (1/cm)
300    # temperature in K
10 1   #  of snaps needed """

DEFAULT_SUPERCELL = """\
3 0 0
0 3 0
0 0 3"""

DEFAULT_VISUALIZER = '~/bin/jmol'
DEFAULT_CONVERTER = '~/BIN/poscar2xyz.py'


def _load_defaults_json(workdir=None):
    """Load optional defaults.json from WORKDIR.
    Returns {} if the file is absent or invalid.
    """
    workdir = workdir or WORKDIR
    path = os.path.join(workdir, 'defaults.json')
    if not os.path.isfile(path):
        return {}
    try:
        data = json.loads(Path(path).read_text())
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _effective_defaults(workdir=None):
    """Built-in defaults overridden by defaults.json when present."""
    dj = _load_defaults_json(workdir)
    return {
        'execpath': dj.get('execpath', DEFAULT_EXEC),
        'visualizer': dj.get('visualizer', DEFAULT_VISUALIZER),
        'converter': dj.get('converter', DEFAULT_CONVERTER),
        'cell': dj.get('cell', DEFAULT_CELL),
        'snaps': dj.get('snaps', DEFAULT_SNAPS),
        'supercell': dj.get('supercell', DEFAULT_SUPERCELL),
    }

# ── job management ────────────────────────────────────────────────────────────
_jobs: dict[str, queue.Queue] = {}
_jobs_lock = threading.Lock()
_job_counter = 0

def _next_job_id():
    global _job_counter
    _job_counter += 1
    return f'job_{_job_counter}'

def _run_job(job_id: str, cmd: list, cwd: str):
    q = _jobs[job_id]
    try:
        import subprocess
        proc = subprocess.Popen(
            cmd, cwd=cwd,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1
        )
        for line in iter(proc.stdout.readline, ''):
            q.put(('out', line.rstrip()))
        proc.wait()
        q.put(('done', str(proc.returncode)))
    except Exception as exc:
        q.put(('err', str(exc)))
        q.put(('done', '1'))

# ── cell/supercell parsing helpers ────────────────────────────────────────────

def _parse_cell_inp():
    """Return (element_names, primitive_counts, masses) from WORKDIR/cell.inp."""
    path = os.path.join(WORKDIR, 'cell.inp')
    lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    ntypes = int(lines[3].split()[0])
    counts = [int(x) for x in lines[4].split()[:ntypes]]
    masses = [float(x) for x in lines[5].split()[:ntypes]]
    names  = lines[6].split()[:ntypes]
    return names, counts, masses

def _parse_supercell_inp():
    """Return abs(det(supercell_matrix)) — number of primitive cells in supercell."""
    path = os.path.join(WORKDIR, 'supercell.inp')
    lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    m = [[int(x) for x in lines[i].split()[:3]] for i in range(3)]
    det = (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
         - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
         + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    return abs(det)

# ── API ───────────────────────────────────────────────────────────────────────

@app.route('/api/run', methods=['POST'])
def api_run():
    d = request.json or {}
    workdir = _request_workdir(d)
    defaults = _effective_defaults(workdir)
    execpath = os.path.expanduser(d.get('execpath', defaults['execpath']).strip())

    if not os.path.isfile(execpath):
        return jsonify(error=f'Executable not found: {execpath}'), 404
    if not os.access(execpath, os.X_OK):
        return jsonify(error=f'Executable not executable (check permissions): {execpath}'), 400

    Path(os.path.join(workdir, 'cell.inp')).write_text(d.get('cell', defaults['cell']))
    Path(os.path.join(workdir, 'snaps.inp')).write_text(d.get('snaps', defaults['snaps']))
    Path(os.path.join(workdir, 'supercell.inp')).write_text(d.get('supercell', defaults['supercell']))

    job_id = _next_job_id()
    with _jobs_lock:
        _jobs[job_id] = queue.Queue()

    t = threading.Thread(target=_run_job, args=(job_id, [execpath], workdir), daemon=True)
    t.start()
    return jsonify(job_id=job_id)


@app.route('/api/stream/<job_id>')
def api_stream(job_id):
    def generate():
        with _jobs_lock:
            q = _jobs.get(job_id)
        if q is None:
            yield f'data: {json.dumps({"type":"err","line":"Job not found"})}\n\n'
            return
        while True:
            try:
                typ, line = q.get(timeout=60)
                yield f'data: {json.dumps({"type":typ,"line":line})}\n\n'
                if typ == 'done':
                    break
            except queue.Empty:
                yield 'data: {"type":"ping"}\n\n'
    return Response(generate(), mimetype='text/event-stream',
                    headers={'Cache-Control': 'no-cache', 'X-Accel-Buffering': 'no'})


@app.route('/api/load_inputs', methods=['POST'])
def api_load_inputs():
    d = request.json or {}
    workdir = _request_workdir(d)
    defaults = _effective_defaults(workdir)
    result = {}
    for fname, key in [('cell.inp','cell'), ('snaps.inp','snaps'), ('supercell.inp','supercell')]:
        p = os.path.join(workdir, fname)
        result[key] = Path(p).read_text(errors='replace').strip() if os.path.isfile(p) else defaults[key]
    result['execpath'] = defaults['execpath']
    result['visualizer'] = defaults['visualizer']
    result['converter'] = defaults['converter']
    result['defaults_source'] = 'defaults.json' if os.path.isfile(os.path.join(workdir, 'defaults.json')) else 'built-in'
    result['workdir'] = workdir
    return jsonify(result)


@app.route('/api/save_inputs', methods=['POST'])
def api_save_inputs():
    """Write the three input files to disk without running sc_snaps.x."""
    d = request.json or {}
    workdir = _request_workdir(d)
    defaults = _effective_defaults(workdir)
    Path(os.path.join(workdir, 'cell.inp')).write_text(d.get('cell', defaults['cell']))
    Path(os.path.join(workdir, 'snaps.inp')).write_text(d.get('snaps', defaults['snaps']))
    Path(os.path.join(workdir, 'supercell.inp')).write_text(d.get('supercell', defaults['supercell']))
    return jsonify(ok=True)


@app.route('/api/files', methods=['POST'])
def api_files():
    d = request.json or {}
    workdir = _request_workdir(d)
    poscars = sorted(glob.glob(os.path.join(workdir, 'poscar_*')))
    entries = [{'name': os.path.basename(p), 'type': 'poscar'} for p in poscars]

    snapshots_path = os.path.join(workdir, 'snapshots.xyz')
    if os.path.isfile(snapshots_path):
        entries.append({'name': 'snapshots.xyz', 'type': 'xyz'})

    return jsonify(files=entries)


def _is_shell_script(path: str) -> bool:
    """Return True if the file starts with a #! shebang line."""
    try:
        with open(path, 'rb') as f:
            return f.read(2) == b'#!'
    except Exception:
        return False

def _launch_jmol(visualizer: str, filepath: str, jmolhome: str = ''):
    """Open a file in Jmol. Handles .app bundles, plain binaries, shell scripts.
    If jmolhome is provided, invokes java -jar Jmol.jar directly (bypasses the
    shell script's JMOL_HOME lookup which may fail outside a login shell).
    """
    import subprocess
    env = os.environ.copy()
    if jmolhome:
        jar = os.path.join(os.path.expanduser(jmolhome), 'Jmol.jar')
        subprocess.Popen(['java', '-Xmx512m', '-jar', jar, filepath], env=env)
        return
    if visualizer.endswith('.app'):
        subprocess.Popen(['open', '-a', visualizer, filepath], env=env)
    elif visualizer.endswith('.sh') or _is_shell_script(visualizer):
        subprocess.Popen(['sh', visualizer, filepath], env=env)
    else:
        subprocess.Popen([visualizer, filepath], env=env)


@app.route('/api/open_visual', methods=['POST'])
def api_open_visual():
    """Convert POSCAR → XYZ if needed, then open in Jmol."""
    d          = request.json or {}
    workdir    = _request_workdir(d)
    visualizer = os.path.expanduser(d.get('visualizer', '').strip())
    converter  = os.path.expanduser(d.get('converter',  '').strip())
    jmolhome   = os.path.expanduser(d.get('jmolhome',   '').strip())
    fname      = os.path.basename(d.get('filename', ''))
    filepath   = os.path.join(workdir, fname)

    if not visualizer:
        return jsonify(error='No visualizer path set.'), 400
    if not os.path.isfile(filepath):
        return jsonify(error=f'File not found: {filepath}'), 404

    try:
        if fname.endswith('.xyz'):
            _launch_jmol(visualizer, filepath, jmolhome)
            return jsonify(ok=True, converted=False)
        else:
            if not converter:
                return jsonify(error='poscar2xyz.py path not set.'), 400
            if not os.path.isfile(converter):
                return jsonify(error=f'Converter not found: {converter}'), 404
            import subprocess
            result = subprocess.run(
                [sys.executable, converter, fname],
                cwd=workdir, capture_output=True, text=True
            )
            if result.returncode != 0:
                msg = result.stderr.strip() or result.stdout.strip() or 'conversion failed'
                return jsonify(error=f'poscar2xyz.py error: {msg}'), 500
            xyz_path = filepath + '.xyz'
            if not os.path.isfile(xyz_path):
                return jsonify(error=f'Converter ran but {fname}.xyz was not created.'), 500
            _launch_jmol(visualizer, xyz_path, jmolhome)
            return jsonify(ok=True, converted=True, xyz=fname + '.xyz')
    except Exception as exc:
        return jsonify(error=str(exc)), 500


@app.route('/api/read_file', methods=['POST'])
def api_read_file():
    d     = request.json or {}
    workdir = _request_workdir(d)
    fname = os.path.basename(d.get('filename', ''))
    path  = os.path.join(workdir, fname)
    if not os.path.isfile(path):
        return jsonify(error='File not found'), 404
    return jsonify(content=Path(path).read_text(errors='replace'))


# ── analysis plots ────────────────────────────────────────────────────────────



def _read_poscar_equilibrium(poscar_path):
    """
    Read equilibrium structure from POSCAR/CONTCAR-like file.

    Supports the common VASP 5 style used by poscar_000:
      line 1  comment
      line 2  scale
      line 3-5 lattice vectors
      line 6  species names
      line 7  species counts
      optional line 8 = Selective dynamics
      next line = Direct/Cartesian
      then natoms lines of positions (extra columns ignored)

    Returns
    -------
    lattice : (3, 3) ndarray
        Lattice matrix in Cartesian Å as row vectors.
    species_order : list[str]
        One symbol per atom, in POSCAR order.
    eq_cart : (natoms, 3) ndarray
        Equilibrium positions in Cartesian Å.
    """
    import numpy as np

    lines = [ln.strip() for ln in Path(poscar_path).read_text().splitlines() if ln.strip()]
    if len(lines) < 8:
        raise ValueError(f'{os.path.basename(poscar_path)} is too short to be a valid POSCAR')

    scale = float(lines[1].split()[0])
    lattice = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)], dtype=float)
    lattice *= scale

    names = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    if len(names) != len(counts):
        raise ValueError(f'{os.path.basename(poscar_path)}: mismatch between species names and counts')

    idx = 7
    if lines[idx].lower().startswith('s'):
        idx += 1

    mode = lines[idx].lower()
    idx += 1

    n_total = sum(counts)
    coord_lines = lines[idx:idx + n_total]
    if len(coord_lines) < n_total:
        raise ValueError(f'{os.path.basename(poscar_path)}: not enough coordinate lines')

    coords = np.array([[float(x) for x in ln.split()[:3]] for ln in coord_lines], dtype=float)

    if mode.startswith('d'):
        eq_cart = coords @ lattice
    elif mode.startswith('c') or mode.startswith('k'):
        eq_cart = coords * scale if abs(scale - 1.0) > 1e-14 else coords
    else:
        raise ValueError(f'{os.path.basename(poscar_path)}: unknown coordinate mode "{lines[idx-1]}"')

    species_order = []
    for name, count in zip(names, counts):
        species_order.extend([name] * count)

    return lattice, species_order, eq_cart


def _read_snapshots_xyz(xyz_path):
    """
    Read snapshots.xyz with lines of the form:
        symbol  x  y  z  charge  vx  vy  vz

    Returns
    -------
    species_order : list[str]
        Species order from the first snapshot.
    positions : (n_snaps, n_atoms, 3) ndarray
        Cartesian coordinates in Å from columns 2,3,4.
    velocities : (n_snaps, n_atoms, 3) ndarray
        Cartesian velocities from columns 6,7,8.
    """
    import numpy as np

    with open(xyz_path, 'r') as f:
        raw = [ln.rstrip() for ln in f if ln.strip()]

    if not raw:
        raise ValueError(f'{os.path.basename(xyz_path)} is empty')

    p = 0
    species_ref = None
    pos_blocks = []
    vel_blocks = []
    snap_idx = 0

    while p < len(raw):
        try:
            n_atoms = int(raw[p].split()[0])
        except Exception as exc:
            raise ValueError(f'{os.path.basename(xyz_path)}: expected atom count at line {p+1}') from exc

        if p + 2 + n_atoms > len(raw):
            raise ValueError(f'{os.path.basename(xyz_path)}: truncated snapshot starting at line {p+1}')

        atom_lines = raw[p + 2:p + 2 + n_atoms]
        species = []
        pos = np.zeros((n_atoms, 3), dtype=float)
        vel = np.zeros((n_atoms, 3), dtype=float)

        for i, ln in enumerate(atom_lines):
            cols = ln.split()
            if len(cols) < 8:
                raise ValueError(
                    f'{os.path.basename(xyz_path)}: line "{ln}" has fewer than 8 columns; '
                    'expected: symbol x y z charge vx vy vz'
                )
            species.append(cols[0])
            pos[i] = [float(cols[1]), float(cols[2]), float(cols[3])]
            vel[i] = [float(cols[5]), float(cols[6]), float(cols[7])]

        if species_ref is None:
            species_ref = species
        elif species != species_ref:
            raise ValueError(
                f'{os.path.basename(xyz_path)}: species order changes at snapshot {snap_idx + 1}; '
                'histogram splitting assumes fixed atom ordering'
            )

        pos_blocks.append(pos)
        vel_blocks.append(vel)
        snap_idx += 1
        p += n_atoms + 2

    return species_ref, np.stack(pos_blocks, axis=0), np.stack(vel_blocks, axis=0)


def _species_blocks(species_order):
    """
    Turn a per-atom species list like [Ba, Ba, ..., O, O, ...]
    into consecutive blocks [(Ba, 32), (O, 32)].

    Raises if the atoms are not grouped by species in contiguous blocks.
    """
    blocks = []
    seen = set()
    i = 0
    n = len(species_order)

    while i < n:
        name = species_order[i]
        if name in seen:
            raise ValueError(
                f'Species "{name}" appears in more than one block. '
                'This code assumes atoms are sorted by type in contiguous groups.'
            )
        j = i
        while j < n and species_order[j] == name:
            j += 1
        blocks.append((name, j - i))
        seen.add(name)
        i = j

    return blocks


def _minimum_image_displacements(cart, eq_cart, lattice):
    """
    Compute wrapped displacements cart - eq_cart using the minimum-image convention
    in fractional coordinates for a general 3D cell.

    Parameters
    ----------
    cart : (n_snaps, n_atoms, 3)
    eq_cart : (n_atoms, 3)
    lattice : (3, 3) row-vector lattice in Å

    Returns
    -------
    disp_cart : (n_snaps, n_atoms, 3)
        Wrapped Cartesian displacement vectors in Å.
    """
    import numpy as np

    inv_lat = np.linalg.inv(lattice)
    dcart = cart - eq_cart[None, :, :]
    dfrac = dcart @ inv_lat
    dfrac -= np.round(dfrac)
    return dfrac @ lattice


def _load_histogram_data_from_xyz(poscar_path, xyz_path):
    """
    Single source of truth for the histogram data.

    Uses:
      - poscar_000 for equilibrium positions and species order
      - snapshots.xyz for actual positions and velocities

    Returns
    -------
    disp_by_element : dict[name -> 1D ndarray]
        Displacement magnitudes |r-r0| in Å for each element, pooled over atoms and snapshots.
    velcomp_by_element : dict[name -> 1D ndarray]
        Velocity Cartesian components vx, vy, vz pooled together for each element.
    names : list[str]
    counts : list[int]
    n_snaps : int
    """
    import numpy as np

    lattice, poscar_species, eq_cart = _read_poscar_equilibrium(poscar_path)
    xyz_species, positions, velocities = _read_snapshots_xyz(xyz_path)

    n_atoms = len(poscar_species)
    if len(xyz_species) != n_atoms:
        raise ValueError(
            f'Atom-count mismatch: poscar_000 has {n_atoms} atoms, '
            f'but snapshots.xyz has {len(xyz_species)} atoms per snapshot'
        )

    if poscar_species != xyz_species:
        raise ValueError(
            'Species order mismatch between poscar_000 and snapshots.xyz. '
            'Histogram splitting by atom type would be ambiguous.'
        )

    disp_cart = _minimum_image_displacements(positions, eq_cart, lattice)
    disp_mag = np.linalg.norm(disp_cart, axis=2)

    blocks = _species_blocks(poscar_species)
    names = [b[0] for b in blocks]
    counts = [b[1] for b in blocks]
    n_snaps = positions.shape[0]

    disp_by_element = {}
    velcomp_by_element = {}

    start = 0
    for name, count in blocks:
        sl = slice(start, start + count)
        disp_by_element[name] = disp_mag[:, sl].reshape(-1)
        velcomp_by_element[name] = velocities[:, sl, :].reshape(-1)
        start += count

    return disp_by_element, velcomp_by_element, names, counts, n_snaps


def _make_displacement_histogram_png(poscar_path, xyz_path):
    """
    Overlaid displacement-magnitude histograms by element.
    Uses |r-r0| with minimum-image wrapping, so this is a true histogram of
    displacement from equilibrium rather than raw coordinates.
    """
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from io import BytesIO

    disp_by_element, _, names, counts, n_snaps = _load_histogram_data_from_xyz(poscar_path, xyz_path)
    colors = plt.cm.tab10.colors

    all_vals = np.concatenate(list(disp_by_element.values()))
    if np.allclose(all_vals.std(), 0.0):
        x_max = max(1e-6, float(all_vals.max()) * 1.2 + 1e-6)
        bins = np.linspace(0.0, x_max, 40)
    else:
        bin_width = 3.49 * all_vals.std() * len(all_vals) ** (-1/3)
        if not np.isfinite(bin_width) or bin_width <= 0:
            bin_width = max((all_vals.max() - all_vals.min()) / 40.0, 1e-6)
        x_lo = max(0.0, float(all_vals.min()) - 0.02 * (all_vals.max() - all_vals.min() + 1e-12))
        x_hi = float(all_vals.max()) + 0.04 * (all_vals.max() - all_vals.min() + 1e-12)
        bins = np.arange(x_lo, x_hi + bin_width, bin_width)
        if len(bins) < 10:
            bins = np.linspace(x_lo, max(x_hi, x_lo + 1e-6), 40)

    fig, ax = plt.subplots(figsize=(9, 4.5))

    for i, name in enumerate(names):
        d = disp_by_element[name]
        color = colors[i % len(colors)]
        mu = float(np.mean(d))
        sigma = float(np.std(d, ddof=1)) if len(d) > 1 else 0.0
        n = len(d)
        se_mu = sigma / np.sqrt(n) if n > 0 else 0.0
        lbl = f'{name} ({counts[i]} atoms)   mean = {mu:.4f}±{se_mu:.4f} Å   std = {sigma:.4f} Å'
        ax.hist(d, bins=bins, density=True, alpha=0.35, color=color, edgecolor='none', label=lbl)

    ax.set_xlabel('Displacement from equilibrium, |r - r₀| (Å)', fontsize=12)
    ax.set_ylabel('Probability density', fontsize=12)
    ax.set_title(f'Thermal displacement histogram   ({n_snaps} snapshots)', fontsize=13)
    ax.legend(fontsize=9, framealpha=0.9)
    fig.tight_layout()

    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=120)
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def _make_displacement_component_histogram_png(poscar_path, xyz_path):
    """
    Overlaid histograms of Cartesian displacement components by element.
    For each element, dx, dy, dz are pooled together after minimum-image wrapping,
    so a harmonic thermal distribution should appear as a zero-centered Gaussian.
    """
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.stats import norm
    from io import BytesIO

    lattice, poscar_species, eq_cart = _read_poscar_equilibrium(poscar_path)
    xyz_species, positions, _ = _read_snapshots_xyz(xyz_path)

    if poscar_species != xyz_species:
        raise ValueError(
            'Species order mismatch between poscar_000 and snapshots.xyz. '
            'Histogram splitting by atom type would be ambiguous.'
        )

    disp_cart = _minimum_image_displacements(positions, eq_cart, lattice)
    blocks = _species_blocks(poscar_species)
    names = [b[0] for b in blocks]
    counts = [b[1] for b in blocks]
    n_snaps = positions.shape[0]

    dispcomp_by_element = {}
    start = 0
    for name, count in blocks:
        sl = slice(start, start + count)
        dispcomp_by_element[name] = disp_cart[:, sl, :].reshape(-1)
        start += count

    colors = plt.cm.tab10.colors
    fig, ax = plt.subplots(figsize=(9, 4.5))

    all_vals = np.concatenate(list(dispcomp_by_element.values()))
    if np.allclose(all_vals.std(), 0.0):
        x_max = max(1e-6, float(np.max(np.abs(all_vals))) * 1.2 + 1e-6)
        bins = np.linspace(-x_max, x_max, 40)
    else:
        bin_width = 3.49 * all_vals.std() * len(all_vals) ** (-1/3)
        if not np.isfinite(bin_width) or bin_width <= 0:
            bin_width = max((all_vals.max() - all_vals.min()) / 40.0, 1e-6)
        x_lo, x_hi = float(all_vals.min()), float(all_vals.max())
        pad = 0.04 * (x_hi - x_lo + 1e-12)
        bins = np.arange(x_lo - pad, x_hi + pad + bin_width, bin_width)
        if len(bins) < 10:
            bins = np.linspace(x_lo - pad, x_hi + pad, 40)

    x_fine = np.linspace(bins[0], bins[-1], 600)

    for i, name in enumerate(names):
        d = dispcomp_by_element[name]
        color = colors[i % len(colors)]
        mu, sigma = norm.fit(d)
        n = len(d)
        se_mu = sigma / np.sqrt(n)
        se_sigma = sigma / np.sqrt(2 * n)
        lbl = (f'{name} ({counts[i]} atoms)   '
               f'μ = {mu:+.4f}±{se_mu:.4f} Å   '
               f'σ = {sigma:.4f}±{se_sigma:.4f} Å')
        ax.hist(d, bins=bins, density=True, alpha=0.35, color=color, edgecolor='none')
        ax.plot(x_fine, norm.pdf(x_fine, mu, sigma), color=color, lw=2, label=lbl)

    ax.set_xlabel('Displacement component from equilibrium (Å)', fontsize=12)
    ax.set_ylabel('Probability density', fontsize=12)
    ax.set_title(f'Cartesian displacement components   ({n_snaps} snapshots)', fontsize=13)
    ax.legend(fontsize=9, framealpha=0.9)
    fig.tight_layout()

    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=120)
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def _make_velocity_histogram_png(poscar_path, xyz_path):
    """
    Overlaid histograms of velocity Cartesian components by element.
    For each element, vx, vy, vz are pooled together exactly as in the original
    vel.dat intent, but the data are read from snapshots.xyz columns 6,7,8.
    """
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.stats import norm
    from io import BytesIO

    _, vel_by_element, names, counts, n_snaps = _load_histogram_data_from_xyz(poscar_path, xyz_path)
    colors = plt.cm.tab10.colors

    fig, ax = plt.subplots(figsize=(9, 4.5))

    all_vals = np.concatenate(list(vel_by_element.values()))
    if np.allclose(all_vals.std(), 0.0):
        x_max = max(1e-6, float(np.max(np.abs(all_vals))) * 1.2 + 1e-6)
        bins = np.linspace(-x_max, x_max, 40)
    else:
        bin_width = 3.49 * all_vals.std() * len(all_vals) ** (-1/3)
        if not np.isfinite(bin_width) or bin_width <= 0:
            bin_width = max((all_vals.max() - all_vals.min()) / 40.0, 1e-6)
        x_lo, x_hi = float(all_vals.min()), float(all_vals.max())
        pad = 0.04 * (x_hi - x_lo + 1e-12)
        bins = np.arange(x_lo - pad, x_hi + pad + bin_width, bin_width)
        if len(bins) < 10:
            bins = np.linspace(x_lo - pad, x_hi + pad, 40)

    x_fine = np.linspace(bins[0], bins[-1], 600)
    sigmas = {}

    for i, name in enumerate(names):
        d = vel_by_element[name]
        color = colors[i % len(colors)]
        mu, sigma = norm.fit(d)
        n = len(d)
        se_mu = sigma / np.sqrt(n)
        se_sigma = sigma / np.sqrt(2 * n)
        sigmas[name] = sigma
        lbl = (f'{name} ({counts[i]} atoms)   '
               f'μ = {mu:+.2f}±{se_mu:.2f} m/s   '
               f'σ = {sigma:.2f}±{se_sigma:.2f} m/s')
        ax.hist(d, bins=bins, density=True, alpha=0.35, color=color, edgecolor='none')
        ax.plot(x_fine, norm.pdf(x_fine, mu, sigma), color=color, lw=2, label=lbl)

    ax.set_xlabel('Velocity component (m/s)', fontsize=12)
    ax.set_ylabel('Probability density', fontsize=12)
    ax.set_title(f'Velocity distribution   ({n_snaps} snapshots)', fontsize=13)
    ax.legend(fontsize=9, framealpha=0.9)

    if len(sigmas) > 1:
        sorted_names = sorted(sigmas, key=lambda n: -sigmas[n])
        ratio = sigmas[sorted_names[0]] / sigmas[sorted_names[-1]]
        ax.set_title(
            f'σ({sorted_names[0]})/σ({sorted_names[-1]}) = {ratio:.3f}',
            fontsize=9, loc='right', color='grey', style='italic'
        )

    fig.tight_layout()

    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=120)
    plt.close(fig)
    buf.seek(0)
    return buf.read()

@app.route('/api/plot_histogram')
def api_plot_histogram():
    try:
        workdir = _request_workdir()
        png = _make_displacement_histogram_png(
            os.path.join(workdir, 'poscar_000'),
            os.path.join(workdir, 'snapshots.xyz'))
        return Response(png, mimetype='image/png')
    except FileNotFoundError as exc:
        return jsonify(error=f'Missing required file: {exc.filename}'), 404
    except ImportError as exc:
        return jsonify(error=f'Missing library: {exc}  —  pip install matplotlib scipy numpy'), 500
    except Exception as exc:
        return Response(_error_png(str(exc)), mimetype='image/png')


@app.route('/api/plot_disp_component_histogram')
def api_plot_disp_component_histogram():
    try:
        workdir = _request_workdir()
        png = _make_displacement_component_histogram_png(
            os.path.join(workdir, 'poscar_000'),
            os.path.join(workdir, 'snapshots.xyz'))
        return Response(png, mimetype='image/png')
    except FileNotFoundError as exc:
        return jsonify(error=f'Missing required file: {exc.filename}'), 404
    except ImportError as exc:
        return jsonify(error=f'Missing library: {exc}  —  pip install matplotlib scipy numpy'), 500
    except Exception as exc:
        return Response(_error_png(str(exc)), mimetype='image/png')


@app.route('/api/plot_vel_histogram')
def api_plot_vel_histogram():
    try:
        workdir = _request_workdir()
        png = _make_velocity_histogram_png(
            os.path.join(workdir, 'poscar_000'),
            os.path.join(workdir, 'snapshots.xyz'))
        return Response(png, mimetype='image/png')
    except FileNotFoundError as exc:
        return jsonify(error=f'Missing required file: {exc.filename}'), 404
    except ImportError as exc:
        return jsonify(error=f'Missing library: {exc}  —  pip install matplotlib scipy numpy'), 500
    except Exception as exc:
        return Response(_error_png(str(exc)), mimetype='image/png')


@app.route('/api/plot_freqs')
def api_plot_freqs():
    """Return a PNG spike plot of phonon frequencies from freqs.dat."""
    try:
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from io import BytesIO

        workdir = _request_workdir()
        freqs_path = os.path.join(workdir, 'freqs.dat')
        if not os.path.isfile(freqs_path):
            return jsonify(error='freqs.dat not found — run sc_snaps.x first'), 404

        # Skip any header lines (non-numeric first token)
        skiprows = 0
        with open(freqs_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    skiprows += 1
                    continue
                try:
                    float(line.split()[0])
                    break
                except (ValueError, IndexError):
                    skiprows += 1

        data = np.loadtxt(freqs_path, skiprows=skiprows)
        if data.ndim == 1:
            data = data.reshape(1, -1)

        freqs_thz = data[:, 2]                  # column 3: frequency in THz
        freqs_cm  = freqs_thz * 33.3564         # THz → cm⁻¹

        # Group equal frequencies within tolerance to compute degeneracy
        tol = 0.5   # cm⁻¹
        freqs_sorted = np.sort(freqs_cm)
        groups = []   # (centre_freq, degeneracy)
        i = 0
        while i < len(freqs_sorted):
            f0 = freqs_sorted[i]
            j  = i + 1
            while j < len(freqs_sorted) and (freqs_sorted[j] - f0) < tol:
                j += 1
            groups.append((float(np.mean(freqs_sorted[i:j])), j - i))
            i = j

        centres = np.array([g[0] for g in groups])
        degen   = np.array([g[1] for g in groups])

        fig, ax = plt.subplots(figsize=(10, 4))
        pos_mask = centres >= 0
        neg_mask = ~pos_mask

        if pos_mask.any():
            ax.vlines(centres[pos_mask], 0, degen[pos_mask],
                      colors='#0369a1', lw=1.5, alpha=0.85, label='Real')
        if neg_mask.any():
            ax.vlines(centres[neg_mask], 0, degen[neg_mask],
                      colors='#dc2626', lw=1.5, alpha=0.85, label='Imaginary (neg)')

        ax.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
        ax.set_ylabel('Degeneracy', fontsize=12)
        ax.set_title('Phonon frequency spectrum', fontsize=13)
        ax.set_ylim(bottom=0)
        if neg_mask.any():
            ax.legend(fontsize=10)
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        fig.tight_layout()

        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=120)
        plt.close(fig)
        buf.seek(0)
        return Response(buf.read(), mimetype='image/png')

    except ImportError as exc:
        return jsonify(error=f'Missing library: {exc}  —  pip install matplotlib numpy'), 500
    except Exception as exc:
        return Response(_error_png(str(exc)), mimetype='image/png')


@app.route('/api/debug_dat')
def api_debug_dat():
    """Diagnostic: show basic stats for vel.dat."""
    import numpy as np
    try:
        vel_path = os.path.join(WORKDIR, 'vel.dat')
        flat = np.loadtxt(vel_path).flatten()
        return jsonify({
            'n_values':    int(len(flat)),
            'mean':        float(np.mean(flat)),
            'std':         float(np.std(flat)),
            'min':         float(flat.min()),
            'max':         float(flat.max()),
            'first_10':    flat[:10].tolist(),
        })
    except Exception as exc:
        return jsonify(error=str(exc)), 500


@app.route('/api/read_log')
def api_read_log():
    """Return contents of log.dat."""
    workdir = _request_workdir()
    path = os.path.join(workdir, 'log.dat')
    if not os.path.isfile(path):
        return jsonify(error='log.dat not found'), 404
    return jsonify(content=Path(path).read_text(errors='replace'))


def _error_png(msg: str) -> bytes:
    """Return a tiny PNG that shows an error message (used when plot fails)."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from io import BytesIO
        fig, ax = plt.subplots(figsize=(7, 2))
        ax.text(0.5, 0.5, f'Plot error:\n{msg}', ha='center', va='center',
                transform=ax.transAxes, color='red', fontsize=10, wrap=True)
        ax.axis('off')
        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=80)
        plt.close(fig)
        buf.seek(0)
        return buf.read()
    except Exception:
        # absolute fallback: 1×1 white PNG
        return (b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01'
                b'\x00\x00\x00\x01\x08\x02\x00\x00\x00\x90wS\xde\x00\x00'
                b'\x00\x0cIDATx\x9cc\xf8\x0f\x00\x00\x01\x01\x00\x05\x18'
                b'\xd8N\x00\x00\x00\x00IEND\xaeB`\x82')


@app.route('/')
def index():
    # Embed WORKDIR into the page
    html = _HTML.replace('__WORKDIR__', WORKDIR)
    return html


# ── HTML / CSS / JS ──────────────────────────────────────────────────────────
_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SC-Snaps GUI</title>
<style>
:root {
  --bg:#f0f4f8; --card:#fff; --border:#e2e8f0;
  --accent:#0369a1; --accent-light:#e0f2fe;
  --text:#1e293b; --sub:#64748b;
  --code-bg:#f8fafc; --term-bg:#0f172a; --term-fg:#e2e8f0;
}
*{box-sizing:border-box;margin:0;padding:0;}
body{font:14px/1.5 system-ui,sans-serif;background:var(--bg);color:var(--text);}
header{background:var(--accent);color:#fff;padding:14px 28px;display:flex;align-items:baseline;gap:14px;}
header h1{font-size:20px;font-weight:700;letter-spacing:-.3px;}
header span{font-size:13px;opacity:.75;}
.main{max-width:1280px;margin:0 auto;padding:24px 20px;}
.card{background:var(--card);border:1px solid var(--border);border-radius:10px;padding:20px;margin-bottom:20px;}
.card-title{font-size:11.5px;font-weight:700;text-transform:uppercase;letter-spacing:.06em;color:var(--sub);margin-bottom:14px;}
.row{display:flex;gap:12px;align-items:flex-end;flex-wrap:wrap;}
.f{display:flex;flex-direction:column;gap:4px;flex:1;min-width:180px;}
label{font-size:12px;font-weight:600;color:var(--sub);}
.workdir-display{font:13px/1.5 'Courier New',monospace;background:var(--code-bg);
  border:1px solid var(--border);border-radius:6px;padding:7px 10px;color:var(--text);word-break:break-all;}
input[type=text]{border:1px solid var(--border);border-radius:6px;padding:7px 10px;font:inherit;background:var(--card);color:var(--text);width:100%;}
input:focus{outline:none;border-color:var(--accent);box-shadow:0 0 0 3px var(--accent-light);}
textarea{border:1px solid var(--border);border-radius:6px;padding:10px 12px;
  font:13px/1.7 'Courier New',monospace;background:var(--code-bg);color:var(--text);
  resize:vertical;width:100%;tab-size:4;}
textarea:focus{outline:none;border-color:var(--accent);box-shadow:0 0 0 3px var(--accent-light);}
.files-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin-bottom:20px;}
@media(max-width:860px){.files-grid{grid-template-columns:1fr;}}
.btn{border:1px solid var(--border);border-radius:6px;padding:7px 16px;font:inherit;
     cursor:pointer;background:var(--card);color:var(--text);transition:background .15s;white-space:nowrap;}
.btn:hover{background:#f1f5f9;}
.btn:disabled{opacity:.45;cursor:not-allowed;}
.btn-primary{background:var(--accent);color:#fff;border-color:var(--accent);font-weight:600;padding:9px 22px;font-size:14.5px;}
.btn-primary:hover{background:#0284c7;}
.btn-sm{padding:4px 10px;font-size:12px;}
.btn-ghost{background:transparent;border-color:var(--border);color:var(--sub);}
.btn-ghost:hover{background:#f1f5f9;color:var(--text);}
.run-row{display:flex;gap:14px;align-items:center;margin-bottom:20px;flex-wrap:wrap;}
.run-status{font-size:13px;color:var(--sub);}
.run-status.ok{color:#15803d;font-weight:600;}
.run-status.err{color:#dc2626;font-weight:600;}
/* terminal log */
.term{background:var(--term-bg);color:var(--term-fg);font:13px/1.6 'Courier New',monospace;
      padding:16px;border-radius:8px;height:280px;overflow-y:auto;white-space:pre-wrap;word-break:break-all;}
.term .line-err{color:#f87171;}
.term .line-done-ok{color:#86efac;font-weight:600;}
.term .line-done-err{color:#f87171;font-weight:600;}
/* output files */
.chip-grid{display:flex;flex-wrap:wrap;gap:8px;}
.chip{background:var(--code-bg);border:1px solid var(--border);border-radius:6px;
      padding:5px 13px;font:12.5px 'Courier New',monospace;cursor:pointer;transition:all .15s;}
.chip:hover{background:var(--accent-light);border-color:var(--accent);color:var(--accent);}
.chip.active{background:var(--accent-light);border-color:var(--accent);color:var(--accent);font-weight:700;}
.viewer{background:var(--code-bg);border:1px solid var(--border);border-radius:8px;
        padding:14px 16px;font:12.5px/1.7 'Courier New',monospace;white-space:pre;
        max-height:420px;overflow:auto;margin-top:14px;display:none;}
.badge{display:inline-block;background:var(--accent);color:#fff;border-radius:12px;
       font-size:11px;font-weight:700;padding:1px 9px;margin-left:8px;vertical-align:middle;}
.hint{font-size:11px;color:var(--sub);margin-top:4px;}
.alert{padding:10px 14px;border-radius:6px;font-size:13px;margin-bottom:14px;}
.alert-err{background:#fee2e2;border:1px solid #fca5a5;color:#991b1b;}
.alert-ok {background:#dcfce7;border:1px solid #86efac;color:#166534;}
.sep{border:none;border-top:1px solid var(--border);margin:4px 0 16px;}
/* analysis plots */
.plot-grid{display:grid;grid-template-columns:1fr 1fr;gap:20px;}
@media(max-width:900px){.plot-grid{grid-template-columns:1fr;}}
.plot-box{border:1px solid var(--border);border-radius:8px;overflow:hidden;background:var(--code-bg);}
.plot-box img{width:100%;display:block;}
.plot-label{font-size:11.5px;font-weight:600;color:var(--sub);padding:8px 12px 6px;
            text-transform:uppercase;letter-spacing:.04em;}
.plot-missing{padding:32px;text-align:center;color:var(--sub);font-size:13px;}
</style>
</head>
<body>
<header>
  <h1>SC-Snaps GUI</h1>
  <span>Supercell snapshot generator interface</span>
</header>

<div class="main">
  <div id="alert-top"></div>

  <!-- ── Configuration ─────────────────────────────────────────────────── -->
  <div class="card">
    <div class="card-title">Configuration</div>

    <!-- Working directory — read-only, fixed at launch -->
    <div class="f" style="margin-bottom:16px;">
      <label>Working Directory — files are read and written here (set at launch time)</label>
      <input type="text" id="workdir" value="__WORKDIR__">
      <div class="hint">Starts from the directory where this script was launched. Edit it if you want the GUI to use a different working directory.</div>
    </div>

    <hr class="sep">

    <div class="row">
      <div class="f" style="flex:2;">
        <label>sc_snaps.x — full path to executable</label>
        <input type="text" id="execpath" value="~/BIN/sc_snaps.x">
      </div>
      <div class="f" style="flex:2;">
        <label>Jmol — full path to executable or .app bundle</label>
        <input type="text" id="visualizer" value="~/bin/jmol">
      </div>
 <!--       <div class="f" style="flex:2;">
        <label>JMOL_HOME — directory containing Jmol.jar</label>
        <input type="text" id="jmolhome" placeholder="/path/to/jmol-16.x.x">
        <div class="hint">Set if the jmol script can't find Jmol.jar (e.g. ~/Downloads/jmol-16.1.47).</div>
      </div>  -->
      <div class="f" style="flex:2;">
        <label>POSCAR→XYZ converter — path to poscar2xyz.py</label>
        <div class="hint">Run automatically before opening snapshots in Jmol.</div>
        <input type="text" id="converter" value="~/BIN/poscar2xyz.py">
      </div>
    </div>
  </div>

  <!-- ── Save bar ───────────────────────────────────────────────────────── -->
  <div style="display:flex;gap:12px;align-items:center;margin-bottom:16px;">
    <button class="btn" id="btn-save" onclick="saveInputs()">💾 Save input files</button>
    <button class="btn btn-ghost" onclick="loadInputs()">⬆ Load existing files</button>
    <button class="btn btn-ghost" onclick="resetDefaults()">↺ Reset defaults</button>
    <span id="save-status" style="font-size:12px;color:var(--sub);"></span>
  </div>

  <!-- ── Input files ────────────────────────────────────────────────────── -->
  <div class="files-grid">

    <div class="card">
      <div class="card-title">cell.inp — primitive cell</div>
      <textarea id="cell" rows="12" spellcheck="false">1 1 1   90 90 90
 0 0.5 0.5, 0.5 0 0.5, 0.5 0.5 0
4.247  9
2  
1 1
24.31  16.00
Mg O
  0 0 0
  0.5 0.5 0.5</textarea>
      <div class="hint" style="margin-top:8px;">
        Line 1: conventional cell a b c α β γ<br>
        Line 2: primitive vectors (in terms of conventional)<br>
        Line 3: lattice parameter scale (Å), 0 for conventional units, otherwise primitive units<br>
        Line 4: number of atom types<br>
        Line 5: number of atoms of each type<br>
        Line 6: atomic masses; charges <br>
        Line 7: element names<br>
        Lines 8+: reduced coordinates (conventional lattice if 0 on Line 3, primitive otherwise)
      </div>
    </div>

    <div class="card">
      <div class="card-title">snaps.inp — snapshot parameters</div>
      <textarea id="snaps" rows="5" spellcheck="false">400  Avg frequency (1/cm)
300    # temperature in K
10     #  of snaps needed </textarea>
      <div class="hint" style="margin-top:8px;">
        Line 1: average desired phonon frequency (cm⁻¹)<br>
        Line 2: temperature (K) to adjust displacements magnitudes<br>
        Line 3: number of desired snapshots <br>
        <br>
        Output: <code>poscar_000</code> (unshifted) through <code>poscar_N</code>
      </div>
    </div>

    <div class="card">
      <div class="card-title">supercell.inp — supercell dimensions</div>
      <textarea id="supercell" rows="5" spellcheck="false">3 0 0
0 3 0
0 0 3</textarea>
      <div class="hint" style="margin-top:8px;">
        3×3 integer matrix: supercell vectors expressed<br>
        in terms of the conventional cell vectors.<br>
        A diagonal matrix <em>n n n</em> gives an n×n×n multiple of the conventional cell.
      </div>
    </div>

  </div>

  <!-- ── Run controls ──────────────────────────────────────────────────── -->
  <div class="run-row">
    <button class="btn btn-primary" id="btn-run" onclick="runSnaps()">▶ Run sc_snaps.x</button>
    <span id="run-status" class="run-status"></span>
  </div>

  <!-- ── Log ───────────────────────────────────────────────────────────── -->
  <div class="card">
    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:10px;">
      <div class="card-title" style="margin-bottom:0;">Output log</div>
      <button class="btn btn-ghost btn-sm" onclick="clearLog()">Clear</button>
    </div>
    <div class="term" id="log">Ready — configure inputs above and click Run.</div>
  </div>

  <!-- ── Generated snapshots ───────────────────────────────────────────── -->
  <div class="card" id="output-card" style="display:none;">
    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:12px;">
      <div class="card-title" style="margin-bottom:0;">
        Generated snapshots <span class="badge" id="snap-count">0</span>
      </div>
      <button class="btn btn-ghost btn-sm" onclick="refreshFiles()">↻ Refresh</button>
    </div>
    <div class="chip-grid" id="chip-grid"></div>
    <div id="viewer-toolbar" style="display:none;margin-top:12px;">
      <button class="btn btn-sm" id="btn-open-viewer" onclick="openInViewer()">👁 Open in viewer</button>
      <span id="viewer-fname" style="font-size:12px;font-family:monospace;margin-left:8px;color:var(--sub);"></span>
    </div>
    <div class="viewer" id="file-viewer"></div>
  </div>

  <!-- ── Analysis plots ────────────────────────────────────────────────── -->
  <div class="card" id="analysis-card">
    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:14px;">
      <div class="card-title" style="margin-bottom:0;">Analysis</div>
      <button class="btn btn-ghost btn-sm" onclick="refreshPlots()">↻ Refresh plots</button>
    </div>
    <div class="plot-grid">
      <div class="plot-box">
        <div class="plot-label">Displacement histogram (snapshots.xyz vs poscar_000)</div>
        <div id="hist-wrap">
          <div class="plot-missing">Need poscar_000 and snapshots.xyz</div>
        </div>
      </div>
      <div class="plot-box">
        <div class="plot-label">Cartesian displacement components (snapshots.xyz vs poscar_000)</div>
        <div id="dispcomp-wrap">
          <div class="plot-missing">Need poscar_000 and snapshots.xyz</div>
        </div>
      </div>
    </div>
    <div class="plot-grid" style="margin-top:16px;">
      <div class="plot-box">
        <div class="plot-label">Velocity distribution (snapshots.xyz)</div>
        <div id="vel-wrap">
          <div class="plot-missing">Need poscar_000 and snapshots.xyz</div>
        </div>
      </div>
    </div>
    <div class="plot-grid" style="margin-top:16px;">
      <div class="plot-box">
        <div class="plot-label">Phonon frequency spectrum (freqs.dat)</div>
        <div id="freqs-wrap">
          <div class="plot-missing">Run sc_snaps.x to generate freqs.dat</div>
        </div>
      </div>
    </div>
    <div style="margin-top:16px;">
      <div class="plot-label" style="padding:0 0 6px;">Run log (log.dat)</div>
      <div id="log-dat" class="term" style="height:200px;">log.dat not loaded</div>
    </div>
  </div>

</div><!-- .main -->

<script>
// ── defaults ────────────────────────────────────────────────────────────────
const DEFAULTS = {
  execpath: `~/BIN/sc_snaps.x`,
  visualizer: `~/bin/jmol`,
  converter: `~/BIN/poscar2xyz.py`,
  cell:
`1 1 1   90 90 90
 0 0.5 0.5, 0.5 0 0.5, 0.5 0.5 0
4.247
2
1 1
24.31  16.00
Mg O
  0 0 0
  0.5 0.5 0.5`,
  snaps:
`400  Avg frequency (1/cm)
300    # temperature in K
10 1   #  of snaps needed and supercell type (default=1)`,
  supercell:
`3 0 0
0 3 0
0 0 3`
};

// ── helpers ─────────────────────────────────────────────────────────────────
const v   = id => document.getElementById(id).value;
const el  = id => document.getElementById(id);
const currentWorkdir = () => (el('workdir') ? v('workdir').trim() : '');

function alertTop(msg, type){
  el('alert-top').innerHTML = msg
    ? `<div class="alert alert-${type}">${msg}</div>` : '';
}

async function post(url, body){
  return fetch(url, {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(body)
  });
}

// ── dirty tracking ───────────────────────────────────────────────────────────
let _savedContent = {cell: null, snaps: null, supercell: null};

function _isDirty(){
  return ['cell','snaps','supercell'].some(k =>
    _savedContent[k] !== null && v(k) !== _savedContent[k]);
}

function _markClean(){
  _savedContent.cell      = v('cell');
  _savedContent.snaps     = v('snaps');
  _savedContent.supercell = v('supercell');
}

function resetDefaults(){
  if(_isDirty() && !confirm('Unsaved edits will be lost. Reset to defaults anyway?')) return;
  el('execpath').value   = DEFAULTS.execpath || el('execpath').value;
  el('visualizer').value = DEFAULTS.visualizer || el('visualizer').value;
  el('converter').value  = DEFAULTS.converter || el('converter').value;
  el('cell').value       = DEFAULTS.cell;
  el('snaps').value      = DEFAULTS.snaps;
  el('supercell').value  = DEFAULTS.supercell;
  _markClean();
  alertTop('Inputs reset to defaults.', 'ok');
}

function clearLog(){
  el('log').innerHTML = '';
}

// ── save input files to disk (without running) ───────────────────────────────
async function saveInputs(){
  const btn = el('btn-save');
  const ss  = el('save-status');
  btn.disabled = true;
  ss.textContent = 'Saving…';
  const r = await post('/api/save_inputs', {
    workdir: currentWorkdir(),
    cell: v('cell'), snaps: v('snaps'), supercell: v('supercell')
  });
  const d = await r.json();
  btn.disabled = false;
  if(!r.ok){ alertTop(d.error || 'Save failed.', 'err'); ss.textContent=''; return; }
  _markClean();
  ss.textContent = '✓ Saved';
  setTimeout(()=>{ ss.textContent=''; }, 3000);
  alertTop('Saved cell.inp, snaps.inp, supercell.inp.', 'ok');
}

// ── load existing input files from working directory ─────────────────────────
async function loadInputs(){
  if(_isDirty() && !confirm(
      'You have unsaved edits.\nLoading will overwrite them. Continue?'))
    return;

  alertTop('', '');
  const r   = await post('/api/load_inputs', {workdir: currentWorkdir()});
  const d   = await r.json();
  let loaded = 0;
  if(d.execpath)   { el('execpath').value = d.execpath; }
  if(d.visualizer) { el('visualizer').value = d.visualizer; }
  if(d.converter)  { el('converter').value = d.converter; }
  if(d.workdir)    { el('workdir').value = d.workdir; }
  if(d.cell)      { el('cell').value = d.cell;           loaded++; }
  if(d.snaps)     { el('snaps').value = d.snaps;         loaded++; }
  if(d.supercell) { el('supercell').value = d.supercell; loaded++; }
  _markClean();
  if(loaded === 0)
    alertTop('No existing input files found — showing defaults.', 'ok');
  else
    alertTop(`Loaded ${loaded}/3 input file(s) from ${currentWorkdir() || d.workdir}.`, 'ok');
  await refreshFiles();
}

// ── output file list ─────────────────────────────────────────────────────────
let _activeFile = null;

async function refreshFiles(){
  const r = await post('/api/files', {workdir: currentWorkdir()});
  const d = await r.json();
  renderFileChips(d.files || []);
}

function renderFileChips(files){
  const card  = el('output-card');
  const grid  = el('chip-grid');
  const count = el('snap-count');
  if(files.length === 0){ card.style.display = 'none'; return; }
  card.style.display = 'block';
  count.textContent  = files.length;
  grid.innerHTML = files.map(f => {
    const name   = f.name || f;
    const active = name === _activeFile ? ' active' : '';
    return `<div class="chip${active}" onclick="viewFile('${name}')">${name}</div>`;
  }).join('');
}

async function viewFile(fname){
  if(_activeFile === fname){
    _activeFile = null;
    el('file-viewer').style.display = 'none';
    el('viewer-toolbar').style.display = 'none';
    document.querySelectorAll('.chip').forEach(c => c.classList.remove('active'));
    return;
  }
  _activeFile = fname;
  document.querySelectorAll('.chip').forEach(c =>
    c.classList.toggle('active', c.textContent.trim() === fname));

  const r = await post('/api/read_file', {workdir: currentWorkdir(), filename: fname});
  const d = await r.json();
  const viewer = el('file-viewer');
  viewer.textContent   = d.error ? `Error: ${d.error}` : d.content;
  viewer.style.display = 'block';
  viewer.scrollTop     = 0;

  if(fname.endsWith('.xyz')){
    el('viewer-fname').textContent = `${fname}  →  Jmol`;
    el('btn-open-viewer').textContent = '👁 Open in Jmol';
  } else {
    el('viewer-fname').textContent = `${fname}  →  ${fname}.xyz  →  Jmol`;
    el('btn-open-viewer').textContent = '👁 Convert & open in Jmol';
  }
  el('viewer-toolbar').style.display = v('visualizer').trim() ? 'block' : 'none';
}

async function openInViewer(){
  const visualizer = v('visualizer').trim();
  const converter  = v('converter').trim();
  if(!_activeFile || !visualizer){ return; }

  const btn = el('btn-open-viewer');
  btn.disabled = true;
  btn.textContent = _activeFile.endsWith('.xyz') ? '⏳ Opening…' : '⏳ Converting…';

  const jmolhomeEl = document.getElementById('jmolhome');
  const jmolhome = jmolhomeEl ? jmolhomeEl.value.trim() : '';
  const r = await post('/api/open_visual', {workdir: currentWorkdir(), filename: _activeFile, visualizer, converter, jmolhome});
  const d = await r.json();
  btn.disabled = false;
  btn.textContent = _activeFile.endsWith('.xyz') ? '👁 Open in Jmol' : '👁 Convert & open in Jmol';

  if(!r.ok){
    alertTop(d.error || 'Could not open visualizer.', 'err');
  } else if(d.converted){
    alertTop(`Converted → ${d.xyz}  |  Opened in Jmol.`, 'ok');
  } else {
    alertTop(`Opened ${_activeFile} in Jmol.`, 'ok');
  }
}

// ── analysis plots ────────────────────────────────────────────────────────────
let _plotsLoaded = false;

function _imgTag(src){
  const wd = encodeURIComponent(currentWorkdir());
  return `<img src="${src}?workdir=${wd}&t=${Date.now()}" alt="plot" loading="lazy">`;
}

async function _loadPlot(apiUrl, wrapId){
  el(wrapId).innerHTML = '<div class="plot-missing">Loading…</div>';
  try {
    const wd = encodeURIComponent(currentWorkdir());
    const r = await fetch(`${apiUrl}?workdir=${wd}&t=${Date.now()}`);
    if(r.ok && (r.headers.get('content-type')||'').includes('image/png')){
      el(wrapId).innerHTML = _imgTag(apiUrl);
    } else {
      const d = await r.json().catch(()=>({error:'Plot failed'}));
      el(wrapId).innerHTML = `<div class="plot-missing" style="color:#dc2626;">${d.error}</div>`;
    }
  } catch(e){
    el(wrapId).innerHTML = `<div class="plot-missing" style="color:#dc2626;">${e}</div>`;
  }
}

async function refreshPlots(){
  await Promise.all([
    _loadPlot('/api/plot_histogram',                'hist-wrap'),
    _loadPlot('/api/plot_disp_component_histogram', 'dispcomp-wrap'),
    _loadPlot('/api/plot_vel_histogram',            'vel-wrap'),
    _loadPlot('/api/plot_freqs',                    'freqs-wrap'),
    (async () => {
      const wd = encodeURIComponent(currentWorkdir());
      const r = await fetch(`/api/read_log?workdir=${wd}`);
      const d = await r.json();
      el('log-dat').textContent = d.error ? d.error : d.content;
    })(),
  ]);
  _plotsLoaded = true;
}

// ── run ──────────────────────────────────────────────────────────────────────
async function runSnaps(){
  alertTop('', '');

  const btn    = el('btn-run');
  const status = el('run-status');
  btn.disabled = true;
  status.textContent = 'Running…';
  status.className   = 'run-status';

  const log = el('log');
  log.innerHTML = '';
  el('output-card').style.display = 'none';
  el('file-viewer').style.display = 'none';
  _activeFile = null;

  const logLine = txt => { const s=document.createElement('span'); s.textContent=txt+'\n'; log.appendChild(s); };
  logLine('  cell.inp      : ' + v('cell').split('\n')[0].trim());
  logLine('  snaps.inp     : ' + v('snaps').split('\n')[0].trim());
  logLine('  supercell.inp : ' + v('supercell').split('\n')[0].trim());
  logLine('');
  _markClean();

  const r = await post('/api/run', {
    workdir:   currentWorkdir(),
    execpath:  v('execpath'),
    cell:      v('cell'),
    snaps:     v('snaps'),
    supercell: v('supercell'),
  });
  const d = await r.json();

  if(!r.ok){
    alertTop(d.error || 'Error starting job.', 'err');
    btn.disabled = false;
    status.textContent = '';
    return;
  }

  const es = new EventSource(`/api/stream/${d.job_id}`);

  es.onmessage = async (ev) => {
    const msg = JSON.parse(ev.data);
    if(msg.type === 'ping') return;

    if(msg.type === 'out' || msg.type === 'err'){
      const span = document.createElement('span');
      span.textContent = msg.line + '\n';
      if(msg.type === 'err') span.className = 'line-err';
      log.appendChild(span);
      log.scrollTop = log.scrollHeight;
    }

    if(msg.type === 'done'){
      es.close();
      const ok = msg.line === '0';
      const span = document.createElement('span');
      span.className   = ok ? 'line-done-ok' : 'line-done-err';
      span.textContent = ok
        ? '\n✓ sc_snaps.x completed successfully.\n'
        : `\n✗ sc_snaps.x exited with code ${msg.line}.\n`;
      log.appendChild(span);
      log.scrollTop = log.scrollHeight;

      btn.disabled       = false;
      status.textContent = ok ? '✓ Done' : `✗ Exit code ${msg.line}`;
      status.className   = 'run-status ' + (ok ? 'ok' : 'err');

      if(ok){
        await refreshFiles();
        await refreshPlots();
      }
    }
  };

  es.onerror = () => {
    es.close();
    const span = document.createElement('span');
    span.className   = 'line-err';
    span.textContent = '\nSSE connection lost.\n';
    log.appendChild(span);
    btn.disabled       = false;
    status.textContent = 'Connection error';
    status.className   = 'run-status err';
  };
}

// ── init ─────────────────────────────────────────────────────────────────────
(async function init(){
  // Auto-load existing input files or defaults.json (silently, no alert)
  const r = await post('/api/load_inputs', {workdir: currentWorkdir()});
  const d = await r.json();
  if(d.workdir)    el('workdir').value    = d.workdir;
  if(d.execpath)   el('execpath').value   = d.execpath;
  if(d.visualizer) el('visualizer').value = d.visualizer;
  if(d.converter)  el('converter').value  = d.converter;
  if(d.cell)       el('cell').value       = d.cell;
  if(d.snaps)      el('snaps').value      = d.snaps;
  if(d.supercell)  el('supercell').value  = d.supercell;
  DEFAULTS.execpath = d.execpath || DEFAULTS.execpath;
  DEFAULTS.visualizer = d.visualizer || DEFAULTS.visualizer;
  DEFAULTS.converter = d.converter || DEFAULTS.converter;
  DEFAULTS.cell = d.cell || DEFAULTS.cell;
  DEFAULTS.snaps = d.snaps || DEFAULTS.snaps;
  DEFAULTS.supercell = d.supercell || DEFAULTS.supercell;
  _markClean();
  // Refresh file list and plots
  await refreshFiles();
  await refreshPlots();
})();
</script>
</body>
</html>
"""

# ── entry point ───────────────────────────────────────────────────────────────
if __name__ == '__main__':
    import webbrowser
    port = 5050
    threading.Timer(0.9, lambda: webbrowser.open(f'http://localhost:{port}')).start()
    print(f'SC-Snaps GUI  →  http://localhost:{port}')
    print(f'Working directory: {WORKDIR}')
    app.run(port=port, debug=False)
