# dcreator 2.0 - Parallel Volterra dislocation inserter
# Copyright (C) 2000-2026 Erik Bitzek
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE and COPYING.HEADER for full license text and citation info.
"""End-to-end test for the Python wrapper.

Skipped automatically if the dcreator binary cannot be found or the
legacy example fixture is missing.
"""
from pathlib import Path

import pytest

from pydcreator import (
    DislocationJob,
    DislocationSpec,
    DisplacementMode,
    DcreatorError,
    find_dcreator_binary,
)

REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "tests" / "regression" / "input.fcc"


def _maybe_skip():
    try:
        find_dcreator_binary()
    except DcreatorError as e:
        pytest.skip(f"dcreator binary not available: {e}")
    if not FIXTURE.is_file():
        pytest.skip(f"fixture missing: {FIXTURE}")


def test_spec_serialization(tmp_path):
    spec = DislocationSpec(
        new_x=(1.0, -1.0, 0.0),
        new_z=(1.0,  1.0, 1.0),
        burgers=(1.76, -1.76, 0.0),
        line=(1.0, 1.0, -2.0),
        glide_plane_normal=(-1.0, -1.0, -1.0),
        origin=(49.78, 0.0, 30.48),
        stripe=50.0,
    )
    text = spec.to_dissparam("in.fcc", "out.data")
    assert "configfile        in.fcc" in text
    assert "displacement_mode upper_half" in text
    assert "burgersv          1.76" in text


def test_run_edge_example(tmp_path):
    _maybe_skip()
    spec = DislocationSpec(
        new_x=(1.0, -1.0, 0.0),
        new_z=(1.0,  1.0, 1.0),
        burgers=(1.76, -1.76, 0.0),
        line=(1.0, 1.0, -2.0),
        glide_plane_normal=(-1.0, -1.0, -1.0),
        origin=(49.7803173955329457, 0.0, 30.4840942132122404),
        stripe=50.0,
        mode=DisplacementMode.UPPER_HALF,
    )
    job = DislocationJob(spec=spec, working_dir=tmp_path)
    out = tmp_path / "out.data"
    cp = job.run(FIXTURE, out)
    assert cp.returncode == 0
    assert out.is_file()
    # Header smoke-check
    head = out.read_text().splitlines()[:6]
    assert any("atoms" in line for line in head)


def test_symmetric_mode_matches_legacy(tmp_path):
    """The legacy `symmetric` mode should match the 1.x reference file
    within floating-point noise (1e-10 Å per atom)."""
    _maybe_skip()
    ref = REPO / "tests" / "regression" / "reference_legacy.fcc"
    if not ref.is_file():
        pytest.skip(f"reference missing: {ref}")

    spec = DislocationSpec(
        new_x=(1.0, -1.0, 0.0),
        new_z=(1.0,  1.0, 1.0),
        burgers=(1.76, -1.76, 0.0),
        line=(1.0, 1.0, -2.0),
        glide_plane_normal=(-1.0, -1.0, -1.0),
        origin=(49.7803173955329457, 0.0, 30.4840942132122404),
        stripe=50.0,
        mode=DisplacementMode.SYMMETRIC,
        input_format="imd",
        output_format="imd",
    )
    out = tmp_path / "out.fcc"
    DislocationJob(spec=spec, working_dir=tmp_path).run(FIXTURE, out)

    def load(path):
        atoms = {}
        in_body = False
        with open(path) as f:
            for line in f:
                if not in_body:
                    if line.startswith("#E"):
                        in_body = True
                    continue
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if len(parts) >= 6:
                    atoms[int(parts[0])] = (
                        float(parts[3]), float(parts[4]), float(parts[5])
                    )
        return atoms

    a = load(out)
    b = load(ref)
    assert set(a) == set(b)
    max_d = max(
        ((a[i][0] - b[i][0]) ** 2
         + (a[i][1] - b[i][1]) ** 2
         + (a[i][2] - b[i][2]) ** 2) ** 0.5
        for i in a
    )
    assert max_d < 1e-10, f"max per-atom delta = {max_d:.3e} > 1e-10"
