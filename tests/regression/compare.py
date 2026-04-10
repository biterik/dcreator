#!/usr/bin/env python3
# dcreator 2.0 - Parallel Volterra dislocation inserter
# Copyright (C) 2000-2026 Erik Bitzek
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE and COPYING.HEADER for full license text and citation info.
"""Compare two IMD-format checkpoint files atom-by-atom.

Exit status 0 iff every atom position agrees within the given tolerance.
Usage: compare.py <a.fcc> <b.fcc> [atol]
"""
import sys
from pathlib import Path


def load(path):
    atoms = {}
    in_body = False
    with open(path) as f:
        for line in f:
            if not in_body:
                if line.startswith('#E'):
                    in_body = True
                continue
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) >= 6:
                atoms[int(parts[0])] = (
                    int(parts[1]), float(parts[2]),
                    float(parts[3]), float(parts[4]), float(parts[5]),
                )
    return atoms


def main():
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    a_path, b_path = sys.argv[1], sys.argv[2]
    atol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10

    a = load(a_path)
    b = load(b_path)
    if set(a.keys()) != set(b.keys()):
        only_a = set(a) - set(b)
        only_b = set(b) - set(a)
        print(f"atom id sets differ: {len(only_a)} only in a, {len(only_b)} only in b")
        return 1

    max_d = 0.0
    worst_id = None
    n_bad = 0
    for k in a:
        _, _, ax, ay, az = a[k]
        _, _, bx, by, bz = b[k]
        dx, dy, dz = ax - bx, ay - by, az - bz
        d = (dx * dx + dy * dy + dz * dz) ** 0.5
        if d > max_d:
            max_d = d
            worst_id = k
        if d > atol:
            n_bad += 1

    print(f"compared {len(a)} atoms")
    print(f"  max ||delta|| = {max_d:.3e}  (worst id {worst_id})")
    print(f"  atoms exceeding tol {atol:.1e}: {n_bad}")
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
