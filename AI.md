# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

IMProToo (Improved MRR Processing Tool) is a Python library that processes raw data from Metek Micro Rain Radars (MRR). It reads MRR "raw data" files (not the vendor's Averaged/Processed products), performs noise removal, peak detection, and Doppler spectrum dealiasing, and derives radar moments (effective reflectivity Ze, Doppler velocity, spectral width) suited for low-SNR targets like snow and light rain. It is distributed on PyPI and versioned with `setuptools_scm` from git tags — there is no hardcoded version string. Requires Python 3.10+.

## Commands

There is no linter config or Makefile in this repo — don't invent flake8/tox invocations that don't exist.

Install for local development (editable, with test dependencies):
```bash
pip install -e ".[test]"
```

Run the test suite:
```bash
pytest tests/
```

Build a distribution:
```bash
python -m build
```

Sanity-check the package imports and version resolve correctly:
```bash
python -c "import IMProToo; print(IMProToo.__version__)"
```

Runtime dependencies: `numpy`, `scipy`, `matplotlib` (plotting only), `netCDF4` (only needed for writing output files). CI runs the test matrix across Python 3.10-3.14 via `.github/workflows/tests.yml`; PyPI releases publish via `.github/workflows/python-publish.yml` on GitHub release, using PyPI Trusted Publishing (OIDC, no token secret).

## Architecture

The package (`IMProToo/`) has two modules, re-exported together via `IMProToo/__init__.py`:

- **`core.py`** — the three main classes, forming a pipeline:
  - `mrrRawData` — reads an MRR raw-data file (optionally gzip-compressed) and parses it into arrays (spectra, heights, transfer function, timestamps, calibration constant).
  - `MrrZe` — takes an `mrrRawData` instance and does the actual signal processing: noise/peak detection (Hildebrand-Sekhon or descending-average method), spectrum dealiasing across the Nyquist range, and computation of Ze/W/spectral width. All tunable behavior lives in the `self.co` dict set up in `__init__` (peak-finding thresholds, dealiasing options, height ranges, netCDF metadata) — this dict is the primary "configuration surface" of the library and must be edited *before* calling `rawToSnow()`. The core call sequence is: construct `MrrZe(rawData)` → optionally `averageSpectra(seconds)` → `rawToSnow()` → `writeNetCDF(path)`.
  - `mrrProcessedData` — separately reads Metek's own "Averaged"/"Processed" MRR output format (not the raw format), for comparison/interoperability purposes.
  - Internal processing methods on `MrrZe` are prefixed with `_` (e.g. `_getPeakHildebrand`, `_dealiaseSpectrum`, `_deAlCoherence`) and are called in sequence from `rawToSnow()`; read that method first to understand the overall processing order before touching any single step.

- **`tools.py`** — stateless helper functions used by `core.py`: time conversion (`date2unix`/`unix2date`), quality-flag file parsing (`getManualQualityArray`, `readQualityFile`), array reshaping (`oneD2twoD`, `limitMaInidces`), and netCDF backend selection (`_get_netCDF_module`; historically supported alternative netCDF3 backends, now always resolves to netCDF4-python since those alternatives are unmaintained and uninstallable on modern Python — `NETCDF3` is kept only as a format-string alias for `NETCDF3_CLASSIC`).

`examples/` contains real-world driver scripts (`batch_convert_rawData.py`, `batch_makeQuicklooks.py`) showing the intended end-to-end usage pattern (raw file → `mrrRawData` → `MrrZe` → averaged → netCDF); these are the best reference for how the classes above are meant to be composed, including site-specific `co` overrides (e.g. `dealiaseSpectrum_heightsWithInterference`). These example scripts are not covered by the test suite or CI and still contain some legacy/unmaintained code paths (e.g. `matplotlib.mlab`, a `Scientific.IO.NetCDF` fallback) — don't assume they run as-is on modern dependency versions.

- **`tests/`** — pytest suite. `test_tools.py` covers the pure helper functions in `tools.py`. `test_core.py` exercises the full `mrrRawData` → `MrrZe.rawToSnow()` → `writeNetCDF()` pipeline against synthetic MRR raw-data fixtures built in `conftest.py` (there is no real MRR raw-data sample file in the repo). One fixture is pure noise (checks the pipeline doesn't crash and produces finite output for real-world messy/no-signal data), the other injects a synthetic Gaussian peak at a known Doppler velocity bin and asserts the recovered `W` matches it — this is a genuine physics-level regression check, not just an import/smoke test. When adding raw-data fixtures, the byte-column layout (9-char fixed-width fields, `startIndex=3` for "new" format) is derived from `mrrRawData.__init__`/`splitMrrRawData` in `core.py`; don't guess at the format without cross-checking the parser.

## Notes

- Requires Python 3.10+; there is no Python 2 support and no `from __future__` imports — don't reintroduce py2 compatibility shims.
- The library is licensed GPLv3; keep the existing copyright/license header style when adding new source files.
- This algorithm is documented in Maahn & Kollias (2012, AMT); non-obvious constants in `MrrZe.co` (e.g. Nyquist velocity bounds, Ze-v dealiasing relations) trace back to that paper and to specific MRR hardware characteristics — don't "simplify" them without understanding the physical justification.
- `.git_archival.txt` (with `export-subst` in `.gitattributes`) lets `setuptools_scm` compute a version when installing from a git archive/tarball that has no `.git` directory (e.g. a GitHub branch tarball). If you ever see version-resolution errors on archive installs, check that file has the `node`/`node-date`/`describe-name`/`ref-names` `$Format:...$` fields intact — a previous, now-removed `setuptools_scm_git_archive` dependency does not provide this on modern setuptools_scm.
- Running `pytest` will print a `DeprecationWarning: Setting the shape on a NumPy array...` from inside netCDF4-python's own `Variable.__setitem__` on NumPy ≥2.5 — this is upstream (netCDF4-python), not something to "fix" in this codebase.
