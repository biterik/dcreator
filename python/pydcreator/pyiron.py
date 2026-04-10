# dcreator 2.0 - Parallel Volterra dislocation inserter
# Copyright (C) 2000-2026 Erik Bitzek
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE and COPYING.HEADER for full license text and citation info.
"""Optional pyiron integration.

Import this module from :mod:`pydcreator.pyiron`. It adds a
:class:`DislocationInserter` job class that wraps :class:`DislocationJob`
inside pyiron's job-management machinery.

pyiron is not a hard dependency of pydcreator — if you don't import this
module, you don't pay the import cost or need pyiron installed.
"""

from __future__ import annotations

from pathlib import Path

try:
    from pyiron_base import TemplateJob  # type: ignore
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "pyiron_base is required for pydcreator.pyiron; install the "
        "'pyiron' extra: pip install 'pydcreator[pyiron]'"
    ) from e

from .job import DislocationJob, DcreatorError, find_dcreator_binary
from .spec import DislocationSpec, DisplacementMode


class DislocationInserter(TemplateJob):
    """pyiron job: insert a dislocation into an existing structure.

    Example
    -------
    >>> job = pr.create.job.DislocationInserter("edge_dipole")
    >>> job.input.spec = DislocationSpec(
    ...     new_x=(1, -1, 0), new_z=(1, 1, 1),
    ...     burgers=(1.76, -1.76, 0), line=(1, 1, -2),
    ...     glide_plane_normal=(-1, -1, -1),
    ...     origin=(49.78, 0, 30.48), stripe=50.0,
    ... )
    >>> job.input.input_file = "Ni_perfect.data"
    >>> job.run()

    The output file is written to ``job.working_directory / "out.data"``
    by default. Override ``job.input.output_file`` to change that.
    """

    def __init__(self, project, job_name):
        super().__init__(project, job_name)
        self.input.spec = None                 # DislocationSpec
        self.input.input_file = ""             # path to perfect crystal
        self.input.output_file = "out.data"    # relative to working dir
        self.input.binary = None               # optional override

    def validate_ready_to_run(self):
        if not self.input.input_file:
            raise ValueError("input.input_file must be set")
        if self.input.spec is None:
            raise ValueError("input.spec must be set to a DislocationSpec")

    def run_static(self):
        self.validate_ready_to_run()
        wd = Path(self.working_directory)
        wd.mkdir(parents=True, exist_ok=True)
        out_path = wd / self.input.output_file

        job = DislocationJob(
            spec=self.input.spec,
            binary=self.input.binary,
            working_dir=wd,
        )
        cp = job.run(self.input.input_file, out_path, capture_output=True)

        # Publish outputs to pyiron's hdf5 store.
        with self.project_hdf5.open("output") as h:
            h["stdout"] = cp.stdout or ""
            h["stderr"] = cp.stderr or ""
            h["output_file"] = str(out_path)

        self.status.finished = True
