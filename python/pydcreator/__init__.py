# dcreator 2.0 - Parallel Volterra dislocation inserter
# Copyright (C) 2000-2026 Erik Bitzek
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE and COPYING.HEADER for full license text and citation info.
"""pydcreator — Python wrapper around the dcreator CLI.

The high-level entry point is :class:`DislocationJob`, which wraps the
compiled `dcreator` binary:

    >>> from pydcreator import DislocationJob, DislocationSpec
    >>> spec = DislocationSpec(
    ...     new_x=(1.0, -1.0, 0.0),
    ...     new_z=(1.0,  1.0, 1.0),
    ...     burgers=(1.76, -1.76, 0.0),
    ...     line=(1.0, 1.0, -2.0),
    ...     glide_plane_normal=(-1.0, -1.0, -1.0),
    ...     origin=(49.78, 0.0, 30.48),
    ...     stripe=50.0,
    ... )
    >>> job = DislocationJob(spec=spec)
    >>> out = job.run("Ni_perfect.fcc", "Ni_edge.data")

The spec and job can also be serialized to a ``.dissparam`` file by
``spec.write()`` so scripts can stay close to the legacy workflow.

If you don't have pyiron installed the module still imports — the
pyiron integration lives in :mod:`pydcreator.pyiron` and is only loaded
on demand.
"""

from .spec import (
    DislocationSpec,
    DisplacementMode,
    RigidSide,
    InputFormat,
    OutputFormat,
    LammpsAtomStyle,
)
from .job import DislocationJob, DcreatorError, find_dcreator_binary

__version__ = "2.0.0"

__all__ = [
    "DislocationSpec",
    "DislocationJob",
    "DcreatorError",
    "DisplacementMode",
    "RigidSide",
    "InputFormat",
    "OutputFormat",
    "LammpsAtomStyle",
    "find_dcreator_binary",
    "__version__",
]
