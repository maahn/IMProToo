import datetime

import numpy as np
import pytest

N_HEIGHTS = 32
N_VELOCITY = 64


def _header_line(dt, cc=1090000, num_spectra=58):
    """One "new" format MRR raw-data header line.

    Field layout is dictated by IMProToo.core.mrrRawData.parseHeaderLine:
    tokens[2] must start with "UTC", the line must end with "RAW", and
    "CC <value>" / "MDQ <x> <numSpectra>" must be present as tokens.
    """
    timestamp = dt.strftime("%y%m%d%H%M%S")
    return f"MRR {timestamp} UTC CC {cc} MDQ 100 {num_spectra} RAW\n"


def _field9_int(value):
    return f"{int(value):9d}"


def _field9_float(value):
    return f"{float(value):9.4f}"


def _height_line(heights):
    # data starts at column 3 for the "new" file format (see splitMrrRawData)
    return "H  " + "".join(_field9_int(h) for h in heights) + "\n"


def _tf_line(transfer_function):
    return "TF " + "".join(_field9_float(v) for v in transfer_function) + "\n"


def _spectrum_bin_line(bin_no, values_over_height):
    return f"F{bin_no:02d}" + "".join(_field9_int(v) for v in values_over_height) + "\n"


def raw_mrr_text(spectra_by_timestep, start=None):
    """Builds "new" format MRR raw-data file content.

    @parameter spectra_by_timestep: sequence of (height, velocity) shaped
        arrays, one per 10s timestep.
    @parameter start: datetime of the first timestep (UTC, naive).
    """
    if start is None:
        start = datetime.datetime(2021, 6, 15, 12, 0, 0)

    heights = [100 * (i + 1) for i in range(N_HEIGHTS)]
    transfer_function = [1.0] * N_HEIGHTS

    lines = []
    for t, spectra in enumerate(spectra_by_timestep):
        assert spectra.shape == (N_HEIGHTS, N_VELOCITY)
        dt = start + datetime.timedelta(seconds=10 * t)
        lines.append(_header_line(dt))
        lines.append(_height_line(heights))
        lines.append(_tf_line(transfer_function))
        for velocity_bin in range(N_VELOCITY):
            lines.append(_spectrum_bin_line(velocity_bin, spectra[:, velocity_bin]))

    return "".join(lines)


@pytest.fixture
def noise_only_raw_file(tmp_path):
    """A few timesteps of low-amplitude noise, no real signal anywhere."""
    rng = np.random.default_rng(42)
    spectra = [
        rng.integers(0, 5, size=(N_HEIGHTS, N_VELOCITY)).astype(float)
        for _ in range(5)
    ]
    content = raw_mrr_text(spectra)
    path = tmp_path / "noise_only.mrr"
    path.write_text(content)
    return path


@pytest.fixture
def peak_raw_file(tmp_path):
    """Timesteps with a synthetic Gaussian peak injected at a known Doppler
    velocity bin, for a band of heights. Used to check the processing
    pipeline actually recovers a signal at the expected velocity.
    """
    rng = np.random.default_rng(7)
    peak_center_bin = 32
    peak_heights = range(5, 16)

    def make_spectra():
        spectra = rng.integers(0, 5, size=(N_HEIGHTS, N_VELOCITY)).astype(float)
        velocity_bins = np.arange(N_VELOCITY)
        peak = 4000 * np.exp(-0.5 * ((velocity_bins - peak_center_bin) / 4) ** 2)
        for h in peak_heights:
            spectra[h] += peak
        return spectra

    spectra = [make_spectra() for _ in range(5)]
    content = raw_mrr_text(spectra)
    path = tmp_path / "peak.mrr"
    path.write_text(content)
    return path, peak_center_bin, peak_heights
