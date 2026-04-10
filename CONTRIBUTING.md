# Contributing to dcreator

Thanks for your interest in improving dcreator. This document explains
how to report bugs, suggest features, and submit code changes.

## Reporting bugs

Please open a GitHub issue with:

- the exact dcreator version (`./dcreator -h` or the tag/commit you built from),
- your OS and compiler version,
- the parameter file, a minimal input file that reproduces the problem,
  and the exact stdout/stderr from the failing run,
- what you expected to happen.

Small, self-contained reproducers make bugs very cheap to fix.

## Suggesting features

If you want a new feature — a new file format, a new displacement mode,
a different parallelization backend — please open an issue first to
discuss the design before writing a lot of code. This saves you time if
the feature overlaps with something already in progress or doesn't fit
the overall direction.

## Submitting changes

dcreator is licensed under the **GNU General Public License v3 or later
(GPL-3.0-or-later)**. By contributing code you agree that your
contribution is distributed under the same license.

1. Fork the repository at <https://github.com/biterik/dcreator> and
   create a feature branch from `main`.
2. Keep the existing code style:
   - C: C11, 4-space indent, braces on the same line, `snake_case` for
     functions and variables, doxygen-free comment style. No tabs.
   - Python: PEP 8, type hints where they help, `black`-compatible
     formatting.
3. Add tests for anything non-trivial you add or fix. New C units go
   under `tests/unit/`; Python tests under `python/tests/`. The
   regression test in `tests/regression/` must continue to pass.
4. Build clean (no new warnings) and run the full test suite locally:
   ```
   cmake -S . -B build -G Ninja
   cmake --build build
   ctest --test-dir build --output-on-failure
   PYTHONPATH=python DCREATOR_BIN=$PWD/build/dcreator pytest python/tests
   ```
5. Open a pull request with a clear description of the problem and the
   approach. Reference the issue you are addressing, if any.

## Reporting improvements back (license requirement)

Because dcreator is GPL-licensed, if you distribute a modified version
of dcreator — inside your group, a collaborating group, a commercial
product, or as part of a paper's supplementary material — you **must**
also distribute the modified source code under the same GPL-3.0-or-later
license. Opening a pull request against this repository is the easiest
way to satisfy this requirement and lets everyone benefit from your fix
or improvement. If upstreaming is not possible for organizational
reasons, please at least send the patch to the author (see README) so it
can be included in a future release.

## Coding conventions cheat sheet

- No hidden globals in new code; the public API uses opaque struct
  pointers.
- Error reporting goes through `dc_errf(errbuf, sz, ...)` so every error
  reaches the caller as a readable string.
- OpenMP directives go on the innermost parallel loop possible to keep
  the serial code path obvious.
- Keep I/O helpers out of the kernel and kernel math out of I/O.

## Questions

If you are unsure about anything, open a draft PR or an issue with the
`question` label. It's always OK to ask before writing code.
