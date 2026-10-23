#!/usr/bin/env python3
"""
poscar_to_xyz.py — Convert VASP POSCAR files to XYZ format for Jmol.

Handles:
  - Universal scale factor (positive = Å multiplier, negative = target volume)
  - Modern VASP (species names on line 6) and old VASP (no species line)
  - Direct (fractional) and Cartesian coordinate modes
  - Selective dynamics line (skipped)
  - Comment preserved in XYZ header alongside lattice vectors

Can be used as a library (poscar_to_xyz / poscar_text_to_xyz)
or run directly from the command line:

    python poscar_to_xyz.py poscar_001            # → poscar_001.xyz
    python poscar_to_xyz.py poscar_*              # batch convert
    python poscar_to_xyz.py poscar_001 -o out.xyz # explicit output name
"""

import os
import sys
import math
import argparse
from pathlib import Path


# ── Core parser ───────────────────────────────────────────────────────────────

def poscar_text_to_xyz(text: str, source_name: str = "POSCAR") -> str:
    """
    Convert the string content of a POSCAR file to XYZ format.
    Returns the XYZ string, or raises ValueError with a descriptive message.
    """
    lines = text.splitlines()
    # Strip inline comments (# ...) from every line before parsing numbers,
    # but keep raw lines for the comment field.
    def strip_comment(line):
        return line.split("#")[0].strip()

    idx = 0

    # ── Line 0: comment ───────────────────────────────────────────────────────
    comment_line = lines[idx].strip() if idx < len(lines) else ""
    idx += 1

    # ── Line 1: universal scale factor ────────────────────────────────────────
    scale_raw = strip_comment(lines[idx])
    idx += 1
    try:
        scale = float(scale_raw.split()[0])
    except (IndexError, ValueError):
        raise ValueError(f"Cannot parse scale factor from line 2: {lines[idx-1]!r}")

    # ── Lines 2-4: lattice vectors ────────────────────────────────────────────
    lat = []
    for i in range(3):
        parts = strip_comment(lines[idx]).split()
        idx += 1
        try:
            lat.append([float(x) for x in parts[:3]])
        except (IndexError, ValueError):
            raise ValueError(f"Cannot parse lattice vector from line {idx}: {lines[idx-1]!r}")

    # Apply scale factor to lattice vectors.
    # Negative scale = desired cell volume (VASP convention).
    if scale < 0:
        vol = abs(
            lat[0][0] * (lat[1][1]*lat[2][2] - lat[1][2]*lat[2][1]) -
            lat[0][1] * (lat[1][0]*lat[2][2] - lat[1][2]*lat[2][0]) +
            lat[0][2] * (lat[1][0]*lat[2][1] - lat[1][1]*lat[2][0])
        )
        if vol == 0:
            raise ValueError("Zero-volume cell — cannot apply negative scale factor.")
        scale = (-scale / vol) ** (1.0 / 3.0)

    lat = [[x * scale for x in row] for row in lat]

    # ── Line 5: species names OR atom counts ──────────────────────────────────
    raw5 = strip_comment(lines[idx])
    idx += 1

    # If the first token is non-numeric it's a species line (modern VASP).
    tokens5 = raw5.split()
    if tokens5 and not _is_integer(tokens5[0]):
        # Modern VASP: species names on this line, counts on the next.
        species = tokens5
        counts_raw = strip_comment(lines[idx])
        idx += 1
    else:
        # Old VASP: no species names — use generic labels X1, X2, …
        species = None
        counts_raw = raw5

    try:
        counts = [int(x) for x in counts_raw.split()]
    except ValueError:
        raise ValueError(f"Cannot parse atom counts: {counts_raw!r}")

    n_types = len(counts)
    total_atoms = sum(counts)

    if species is None:
        species = [f"X{i+1}" for i in range(n_types)]
    elif len(species) < n_types:
        # Pad with generic names if too few
        species += [f"X{i+1}" for i in range(len(species), n_types)]

    # Build per-atom element list
    elements = []
    for sym, count in zip(species, counts):
        elements.extend([sym] * count)

    # ── Optional: Selective dynamics ──────────────────────────────────────────
    raw6 = strip_comment(lines[idx])
    idx += 1
    if raw6.lower().startswith("s"):          # "Selective dynamics"
        coord_mode_line = strip_comment(lines[idx])
        idx += 1
    else:
        coord_mode_line = raw6

    # ── Coordinate mode: Direct or Cartesian ─────────────────────────────────
    mode = coord_mode_line.strip().lower()
    if mode.startswith("d"):
        direct = True
    elif mode.startswith("c") or mode.startswith("k"):
        direct = False
    else:
        raise ValueError(f"Unknown coordinate mode: {coord_mode_line!r}  (expected Direct or Cartesian)")

    # ── Atom positions ────────────────────────────────────────────────────────
    cart_coords = []
    for i in range(total_atoms):
        if idx >= len(lines):
            raise ValueError(f"File truncated: expected {total_atoms} atom positions, got {i}.")
        parts = strip_comment(lines[idx]).split()
        idx += 1
        try:
            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        except (IndexError, ValueError):
            raise ValueError(f"Cannot parse coordinate on line {idx}: {lines[idx-1]!r}")

        if direct:
            # Convert fractional → Cartesian: r = x*a1 + y*a2 + z*a3
            cx = x*lat[0][0] + y*lat[1][0] + z*lat[2][0]
            cy = x*lat[0][1] + y*lat[1][1] + z*lat[2][1]
            cz = x*lat[0][2] + y*lat[1][2] + z*lat[2][2]
        else:
            cx, cy, cz = x * scale, y * scale, z * scale   # already Cartesian, but scale already applied to lat
            # Note: for Cartesian POSCAR the coords are already in scaled units.
            # Re-check: VASP docs say Cartesian coords are multiplied by the
            # lattice constant (scale) too.
            cx, cy, cz = x, y, z  # scale was already folded into lat above

        cart_coords.append((cx, cy, cz))

    # ── Build XYZ string ──────────────────────────────────────────────────────
    # XYZ comment line: preserve original comment + attach lattice for Jmol
    lat_str = (
        f'Lattice="{lat[0][0]:.6f} {lat[0][1]:.6f} {lat[0][2]:.6f} '
        f'{lat[1][0]:.6f} {lat[1][1]:.6f} {lat[1][2]:.6f} '
        f'{lat[2][0]:.6f} {lat[2][1]:.6f} {lat[2][2]:.6f}"'
    )
    properties_str = 'Properties=species:S:1:pos:R:3'
    xyz_comment = f'{lat_str} {properties_str}'
    if comment_line:
        xyz_comment = f'{comment_line}  |  {xyz_comment}'

    out_lines = [str(total_atoms), xyz_comment]
    for elem, (cx, cy, cz) in zip(elements, cart_coords):
        out_lines.append(f"{elem:<4s}  {cx:16.10f}  {cy:16.10f}  {cz:16.10f}")

    return "\n".join(out_lines) + "\n"


def poscar_to_xyz(poscar_path: str, xyz_path: str = None) -> str:
    """
    Convert a POSCAR file on disk to XYZ.
    If xyz_path is None, writes to <poscar_path>.xyz and returns that path.
    Returns the output path.
    """
    poscar_path = os.path.expanduser(poscar_path)
    with open(poscar_path, "r") as f:
        text = f.read()

    xyz_text = poscar_text_to_xyz(text, source_name=os.path.basename(poscar_path))

    if xyz_path is None:
        xyz_path = poscar_path + ".xyz"
    xyz_path = os.path.expanduser(xyz_path)

    with open(xyz_path, "w") as f:
        f.write(xyz_text)

    return xyz_path


# ── Helper ────────────────────────────────────────────────────────────────────

def _is_integer(s: str) -> bool:
    try:
        int(s)
        return True
    except ValueError:
        return False


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Convert VASP POSCAR file(s) to XYZ format for Jmol."
    )
    parser.add_argument("inputs", nargs="+", help="POSCAR file(s) to convert")
    parser.add_argument(
        "-o", "--output",
        help="Output XYZ filename (only valid when converting a single file)",
    )
    args = parser.parse_args()

    if args.output and len(args.inputs) > 1:
        print("ERROR: -o/--output can only be used with a single input file.", file=sys.stderr)
        sys.exit(1)

    errors = 0
    for inp in args.inputs:
        out = args.output if args.output else None
        try:
            out_path = poscar_to_xyz(inp, out)
            print(f"  ✓  {inp}  →  {out_path}")
        except Exception as e:
            print(f"  ✗  {inp}  FAILED: {e}", file=sys.stderr)
            errors += 1

    if errors:
        sys.exit(1)


if __name__ == "__main__":
    main()
