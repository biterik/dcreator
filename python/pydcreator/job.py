# dcreator 2.0 - Parallel Volterra dislocation inserter
# Copyright (C) 2000-2026 Erik Bitzek
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE and COPYING.HEADER for full license text and citation info.
"""DislocationJob — runs the dcreator CLI on a given input file."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .spec import DislocationSpec


class DcreatorError(RuntimeError):
    """Raised when the dcreator binary fails or cannot be located."""


def find_dcreator_binary() -> Path:
    """Locate the dcreator binary in priority order.

    1. ``$DCREATOR_BIN`` environment variable (explicit override).
    2. ``dcreator`` on ``$PATH``.
    3. The sibling ``build/`` directory relative to this module (for
       in-tree development runs).
    """
    env = os.environ.get("DCREATOR_BIN")
    if env:
        p = Path(env).expanduser()
        if p.is_file() and os.access(p, os.X_OK):
            return p
        raise DcreatorError(f"$DCREATOR_BIN={env!r} is not an executable file")

    which = shutil.which("dcreator")
    if which:
        return Path(which)

    here = Path(__file__).resolve()
    # repo/python/pydcreator/job.py -> repo/build/dcreator
    candidate = here.parents[2] / "build" / "dcreator"
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return candidate

    raise DcreatorError(
        "dcreator binary not found. Set $DCREATOR_BIN, add dcreator to "
        "$PATH, or build the CLI in <repo>/build."
    )


@dataclass
class DislocationJob:
    """Run the dcreator CLI against an input file.

    Parameters
    ----------
    spec
        A :class:`DislocationSpec` describing the dislocation and I/O
        modes. May be omitted if you plan to call :meth:`run_paramfile`
        directly with a hand-written ``.dissparam`` file.
    binary
        Path to the dcreator executable. If ``None``, uses
        :func:`find_dcreator_binary`.
    working_dir
        Directory in which to create temporary files. If ``None``,
        a fresh temp directory is used per call.
    """
    spec: Optional[DislocationSpec] = None
    binary: Optional[Path] = None
    working_dir: Optional[Path] = None

    def _bin(self) -> Path:
        return Path(self.binary) if self.binary else find_dcreator_binary()

    def run(self, input_file: str | Path, output_file: str | Path,
            capture_output: bool = True) -> subprocess.CompletedProcess:
        """Run dcreator on ``input_file`` writing ``output_file``.

        The parameter file is generated from :attr:`spec` and dropped
        into a temporary directory. Returns the CompletedProcess from
        subprocess; stdout/stderr are captured by default.

        Raises
        ------
        DcreatorError
            If :attr:`spec` is None or the dcreator exit status is nonzero.
        """
        if self.spec is None:
            raise DcreatorError(
                "DislocationJob.run requires a spec; use run_paramfile "
                "for an explicit .dissparam file"
            )
        input_file = Path(input_file).resolve()
        output_file = Path(output_file).resolve()
        if not input_file.is_file():
            raise DcreatorError(f"input file not found: {input_file}")

        if self.working_dir is not None:
            wd = Path(self.working_dir)
            wd.mkdir(parents=True, exist_ok=True)
            param_path = wd / "dcreator.dissparam"
        else:
            tmpdir = tempfile.mkdtemp(prefix="pydcreator-")
            wd = Path(tmpdir)
            param_path = wd / "dcreator.dissparam"

        param_path.write_text(
            self.spec.to_dissparam(str(input_file), str(output_file))
        )
        return self.run_paramfile(param_path, capture_output=capture_output)

    def run_paramfile(self, param_file: str | Path,
                      capture_output: bool = True) -> subprocess.CompletedProcess:
        """Run dcreator on an existing ``.dissparam`` file."""
        param_file = Path(param_file).resolve()
        if not param_file.is_file():
            raise DcreatorError(f"parameter file not found: {param_file}")

        cp = subprocess.run(
            [str(self._bin()), str(param_file)],
            cwd=param_file.parent,
            capture_output=capture_output,
            text=True,
        )
        if cp.returncode != 0:
            msg = cp.stderr.strip() if capture_output and cp.stderr else ""
            raise DcreatorError(
                f"dcreator exited with status {cp.returncode}"
                + (f": {msg}" if msg else "")
            )
        return cp
