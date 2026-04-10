# dcreator

**Parallel Volterra dislocation inserter for atomistic simulations.**

Author: **Erik Bitzek** ([ORCID 0000-0001-7430-3694](https://orcid.org/0000-0001-7430-3694))
Affiliation: Max Planck Institute for Sustainable Materials, Düsseldorf, Germany
Contact: <e.bitzek@mpi-susmat.de>
License: **GPL-3.0-or-later** (see [LICENSE](LICENSE))
Repository: <https://github.com/biterik/dcreator>

`dcreator` takes an atomistic configuration of a crystal, rotates it
into the frame of a prescribed dislocation, applies a *simplified*
piecewise displacement field that is compatible with periodic
boundaries on a rectangular simulation box, and writes the modified
configuration back out. It is a fast, parallel, format-flexible rewrite
of the original dcreator from ~2000, designed to handle very large
systems (hundreds of millions of atoms) and to be driven from Python
and pyiron workflows as well as from the command line.

---

## Why a simplified displacement field?

Unlike the elastic (Volterra) displacement field of an isotropic or
anisotropic continuum, `dcreator` inserts a dislocation by splitting
the crystal into three regions along the glide direction:

1. an **untouched** region far from the dislocation line,
2. a **linear-ramp** strip of configurable width *s* in which one
   half-space is progressively stretched (or sheared) toward the final
   displacement,
3. a **rigid-shift** region in which the half-space is shifted by the
   full Burgers vector.

The jump of the displacement across the glide plane is exactly the
Burgers vector, so the result is a genuine Volterra dislocation.
However — and this is the key property — the displacement field is
**zero outside the ramp+shift strip**, which means a rectangular
simulation box with periodic boundaries in the line direction and
glide-plane-normal direction remains periodic without any box
distortion. No dipole is needed, no box tilt is required, and no atoms
are added or removed: the only thing that changes is the set of
positions.

This makes `dcreator` particularly well-suited for:

- **flat, rectangular simulation cells** with straightforward periodic
  boundary conditions,
- **studies of concentrated solid solutions**, high-entropy alloys, and
  other compositionally disordered crystals where the interaction of a
  dissociated dislocation with the local chemical environment matters.
  In these systems the stacking-fault width between the two Shockley
  partials of an extended dislocation is a physically relevant length
  scale that varies with local composition, and it is often necessary
  to start from a configuration with a chosen initial splitting width.
  `dcreator` makes this trivial: the splitting between two partials is
  a direct consequence of the strip width *s* used when inserting each
  partial, so the user can scan the initial partial separation
  independently of the equilibrium stacking-fault energy.

The approach is adapted from the construction used by D. Rodney to
introduce screw dislocations in fcc simulation cells; see the
[Acknowledgments](#acknowledgments) section.

---

## Physical model

Given a Burgers vector **b**, line direction **l**, glide-plane normal
**n**, origin **O** on the dislocation line, and ramp width *s*, every
atom at sample-frame position *p* is mapped via

```
p'  = R_act2dis · (p − O)                  # rotate into dislocation frame
t   = (rigid_side == right) ? +p'.x : −p'.x

            t < −s/2       →   u = 0                         (untouched)
   |t|  ≤  s/2             →   u = (t + s/2)/s · U(p'.y)     (linear ramp)
            t >  s/2       →   u = U(p'.y)                   (rigid shift)

p_new = O + R_dis2act · (p' + u)
```

where the dislocation frame uses *z_d = l̂*, *y_d = n̂*, *x_d = ŷ_d × ẑ_d*,
and the half-space selector `U(y)` depends on the `displacement_mode`:

| mode         | U(y > 0) | U(y ≤ 0) | physical meaning                               |
|--------------|----------|----------|------------------------------------------------|
| `upper_half` *(default)* | +**b** |    0     | upper half shifted by the full Burgers vector  |
| `lower_half` |    0     | −**b**   | mirror image                                   |
| `symmetric`  | −**b**/2 | +**b**/2 | legacy dcreator 1.x; reproduces old output     |

`rigid_side` (`right` by default, or `left`) chooses which side of the
dislocation line carries the rigid shift — i.e. which end of the cut
extends to the simulation-cell boundary.

---

## Building from source

Dependencies:

- **CMake ≥ 3.20**
- a **C11 compiler with OpenMP** (see platform notes below)
- **libm**
- (optional) **Python ≥ 3.9** for the `pydcreator` wrapper

### macOS (Apple Silicon or Intel)

Apple's bundled `clang` does not ship with the OpenMP runtime, so we
use Homebrew LLVM:

```bash
brew install llvm libomp cmake ninja
git clone https://github.com/biterik/dcreator.git
cd dcreator
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

`CMakeLists.txt` auto-detects Homebrew LLVM at `/opt/homebrew/opt/llvm`
(arm64) or `/usr/local/opt/llvm` (Intel) and points the compiler at it.

### Linux

Any recent `gcc` or `clang` works:

```bash
sudo apt-get install cmake ninja-build build-essential
git clone https://github.com/biterik/dcreator.git
cd dcreator
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

### CMake options

| option                  | default | effect                                      |
|-------------------------|---------|---------------------------------------------|
| `DCREATOR_ENABLE_OPENMP`| `ON`    | link OpenMP; disables parallelism if off    |
| `DCREATOR_NATIVE_ARCH`  | `OFF`   | compile with `-march=native` (fastest binary, non-portable) |
| `DCREATOR_BUILD_TESTS`  | `ON`    | build and register unit + regression tests  |

Install with:

```bash
cmake --install build --prefix /usr/local
```

---

## Usage — command line

```bash
./build/dcreator path/to/my_dislocation.dissparam
```

A parameter file is a plain-text key/value file. The following example
inserts an edge dislocation into an fcc Ni crystal oriented with
*x* ‖ [1 −1 0], *z* ‖ [1 1 1]:

```
# ---- I/O ----
configfile          Ni_perfect.data
outfile             Ni_edge.data

# ---- Sample orientation (any two axes; third is computed) ----
new_x               1.0 -1.0 0.0
new_z               1.0  1.0 1.0

# ---- Dislocation specification (crystal frame) ----
burgersv            1.76 -1.76 0.0        # a/2 <1-10>, magnitude encoded
linev               1.0  1.0 -2.0         # <11-2>
glideplane         -1.0 -1.0 -1.0         # {111}
point               49.78 0.0 30.48       # point on the dislocation line
stripe              50                    # width of the linear-ramp strip

# ---- Behavior knobs (optional) ----
displacement_mode   upper_half            # default
rigid_side          right                 # default

# ---- Formats (auto by default) ----
input_format        auto                  # or lammps_data / lammps_dump / imd
output_format       lammps_data           # or imd
lammps_atom_style   atomic                # atomic / charge / molecular / full / bond
```

### All recognized parameter keys

| key                    | meaning                                                           |
|------------------------|-------------------------------------------------------------------|
| `configfile`           | input atomistic configuration (path)                              |
| `outfile`              | output file (path)                                                |
| `new_x`, `new_y`, `new_z` | sample orientation basis; supply at least two                 |
| `burgersv`             | full Burgers vector in the crystal frame (direction × magnitude)  |
| `linev`                | dislocation line direction (unnormalized accepted)                |
| `glideplane`           | glide-plane normal (unnormalized accepted)                        |
| `point`                | point on the dislocation line in **sample-frame** coordinates     |
| `stripe`               | width of the linear-ramp strip in the glide direction             |
| `displacement_mode`    | `upper_half` (default), `lower_half`, `symmetric`                 |
| `rigid_side`           | `right` (default), `left`                                         |
| `input_format`         | `auto` (default), `lammps_data`, `lammps_dump`, `imd`             |
| `output_format`        | `lammps_data` (default), `imd`                                    |
| `lammps_atom_style`    | `atomic` (default), `charge`, `molecular`, `full`, `bond`         |
| `num_threads`          | OpenMP thread count; `0` means use `OMP_NUM_THREADS`              |

### Inserting extended (dissociated) dislocations

For an extended fcc dislocation dissociated into two Shockley partials,
run `dcreator` **twice** on the same input file, once per partial,
using a different Burgers vector and origin each time. The distance
between the two `point` parameters sets the initial partial separation;
for a given stacking-fault energy the equilibrium separation can be
approached or exceeded independently. This is the mechanism referred
to in the introduction for studying concentrated solid solutions.

---

## Usage — Python

```bash
pip install -e python/
export DCREATOR_BIN=$PWD/build/dcreator     # or add dcreator to $PATH
```

```python
from pydcreator import DislocationJob, DislocationSpec, DisplacementMode

spec = DislocationSpec(
    new_x=(1.0, -1.0, 0.0),
    new_z=(1.0,  1.0, 1.0),
    burgers=(1.76, -1.76, 0.0),
    line=(1.0, 1.0, -2.0),
    glide_plane_normal=(-1.0, -1.0, -1.0),
    origin=(49.78, 0.0, 30.48),
    stripe=50.0,
    mode=DisplacementMode.UPPER_HALF,
)

job = DislocationJob(spec=spec)
job.run("Ni_perfect.data", "Ni_edge.data")
```

The binary is located via `$DCREATOR_BIN`, then `$PATH`, then the
in-tree `build/dcreator` — so in-development runs don't need any
installation.

### pyiron integration

The optional `pydcreator.pyiron` module exposes a `DislocationInserter`
`TemplateJob` subclass:

```python
from pydcreator import DislocationSpec
from pydcreator.pyiron import DislocationInserter

job = pr.create.job.DislocationInserter("edge_dipole")
job.input.spec = DislocationSpec(...)
job.input.input_file = "Ni_perfect.data"
job.run()
```

`pyiron_base` is only imported when you import `pydcreator.pyiron`, so
users without pyiron pay no cost.

---

## File format support

### Reading

- **LAMMPS data**: `atom_style` `atomic`, `charge`, `molecular`,
  `full`, or `bond`. The header (counts, box, `Masses`, coeffs, …) and
  trailing sections (`Velocities`, `Bonds`, …) are preserved verbatim;
  only the `Atoms` section is modified. Per-atom trailing columns
  (image flags, arbitrary extras) round-trip untouched.
- **LAMMPS dump** (text): the first snapshot only. Column names are
  parsed from the `ITEM: ATOMS` header line; `x/y/z`, `xu/yu/zu`, and
  `xs/ys/zs` (scaled — automatically converted to unscaled) are
  recognized. Any extra columns are kept as trailing text.
- **IMD** (ASCII): the legacy format of the IMD molecular-dynamics
  code, identified by `#F A …` on the first line. `number type [mass]
  x y z [vx vy vz]` layouts are supported. Binary IMD is **not**
  supported in v1.
- **Auto-detection**: the default `input_format auto` sniffs the first
  non-blank lines — `ITEM:` → dump, `#F` → IMD, else data.

### Writing

- **LAMMPS data** (default). On a LAMMPS → LAMMPS round-trip the
  original header and trailing sections are passed through byte-for-byte
  except for the `Atoms` block. When converting from IMD, a minimal
  header is synthesized with box bounds copied from the source and
  0-based atom types shifted up to 1-based so they are valid LAMMPS types.
- **IMD** (ASCII). Round-trips IMD exactly. Synthesized from LAMMPS on
  conversion.

### Box handling

The simulation box is **never** modified by `dcreator`. If the rigid-
shift region moves atoms beyond the original box bounds (which happens
for dislocations whose Burgers vector has a component along the glide
direction), either use a slightly oversized box, or run a
`boundary s s p`-style LAMMPS setup that shrink-wraps the relevant
dimension during relaxation.

---

## Compatibility with the legacy dcreator

Old `.dissparam` files from the 2000-era tool run unchanged. They
produce output in the new default mode (`upper_half`), which differs
from the old ±b/2 split by a rigid translation of one half-space.
To reproduce the old output exactly, add

```
displacement_mode symmetric
```

to the parameter file. The bundled regression test confirms that in
`symmetric` mode, the new tool reproduces the original output on the
`DCREATOR_Example` edge dislocation to within 1×10⁻¹⁰ Å per atom
(floating-point noise) for all 48 000 atoms.

---

## Citation

If you use `dcreator` in published work, **please cite** the paper in
which the method was used for the physically relevant study of extended
dislocations in disordered alloys:

> A. Vaid, D. Wei, **E. Bitzek**, S. Nasiri, M. Zaiser,
> *Pinning of extended dislocations in atomically disordered crystals*,
> **Acta Materialia 236 (2022) 118095**.
> https://doi.org/10.1016/j.actamat.2022.118095

BibTeX:

```bibtex
@article{Vaid2022,
    author  = {Vaid, Aviral and Wei, Ding and Bitzek, Erik and
               Nasiri, Samaneh and Zaiser, Michael},
    title   = {Pinning of extended dislocations in atomically disordered crystals},
    journal = {Acta Materialia},
    volume  = {236},
    pages   = {118095},
    year    = {2022},
    doi     = {10.1016/j.actamat.2022.118095},
}
```

GitHub also exposes a "Cite this repository" button — see
[`CITATION.cff`](CITATION.cff).

---

## Acknowledgments

The piecewise displacement construction implemented in `dcreator` was
inspired by the approach used in:

> D. Rodney, *Molecular dynamics simulation of screw dislocations
> interacting with interstitial frank loops in a model FCC crystal*,
> **Acta Materialia 52 (3) (2004) 607–614**.

The original (2000-era) `dcreator` code from which this rewrite derives
its file-format conventions and test fixtures was used in the author's
group at FAU Erlangen-Nürnberg; the IMD checkpoint examples shipped in
the regression test were generated with the LEGO crystal builder by
Arun Prakash and members of the group of Erik Bitzek.

---

## Funding

This work is associated with the NFDI-MatWerk consortium.

> Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research
> Foundation) under the National Research Data Infrastructure –
> NFDI 38/1 – project number 460247524.

The corresponding record for this software is deposited in the
[NFDI-MatWerk Zenodo community](https://zenodo.org/communities/nfdi-matwerk).

## License and warranty

`dcreator` is free software: you can redistribute it and/or modify it
under the terms of the **GNU General Public License v3 or later**
(GPL-3.0-or-later) as published by the Free Software Foundation. See
[`LICENSE`](LICENSE) for the full text and
[`COPYING.HEADER`](COPYING.HEADER) for the short per-file notice.

**This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.** See the GNU
General Public License for more details.

The GPL requires that if you distribute a modified version of
`dcreator` you also make the source code of your modifications
available under the same license. Please send fixes and improvements
back to the upstream repository as pull requests so that the
community benefits; see [`CONTRIBUTING.md`](CONTRIBUTING.md).

---

## Repository layout

```
dcreator/
├── CMakeLists.txt            CMake build (auto-detects Homebrew LLVM on macOS)
├── LICENSE                   GNU GPL v3 (full text)
├── COPYING.HEADER            short per-file copyright/citation notice
├── CITATION.cff              machine-readable citation (GitHub "Cite this repository")
├── CONTRIBUTING.md           how to report bugs and submit patches
├── README.md                 you are here
├── include/dcreator.h        public C API
├── src/
│   ├── core/                 geometry, kernel, parameters, atoms
│   └── io/                   LAMMPS data / dump / IMD readers and writers
├── src/cli/main.c            the `dcreator` command-line tool
├── python/
│   ├── pyproject.toml
│   ├── pydcreator/           subprocess wrapper, pyiron Job class
│   └── tests/                pytest end-to-end tests
└── tests/
    ├── unit/                 C unit tests (geometry, kernel)
    └── regression/           legacy fixture + 1×10⁻¹⁰ Å reference diff
```

---

## Contact

**Erik Bitzek**
Max Planck Institute for Sustainable Materials, Düsseldorf, Germany
ORCID: <https://orcid.org/0000-0001-7430-3694>
Email: <e.bitzek@mpi-susmat.de>

Please report bugs and request features via the
[GitHub issue tracker](https://github.com/biterik/dcreator/issues).
For scientific questions about the method, feel free to email me directly.
