import numpy as np
import netCDF4
import pytest

import IMProToo


def test_mrrRawData_parses_new_format_raw_file(noise_only_raw_file):
    raw = IMProToo.mrrRawData(str(noise_only_raw_file))

    assert raw.shape2D == (5, 32)
    assert raw.shape3D == (5, 32, 64)
    assert raw.mrrRawCC == 1090000
    assert raw.timezone == "UTC"
    assert not raw.mrrRawHeight.mask.any()
    assert not raw.mrrRawTF.mask.any()
    np.testing.assert_array_equal(raw.mrrRawHeight[0], 100 * np.arange(1, 33))


def test_pipeline_runs_end_to_end_on_noise_only_data(noise_only_raw_file, tmp_path):
    raw = IMProToo.mrrRawData(str(noise_only_raw_file))
    processedSpec = IMProToo.MrrZe(raw)
    processedSpec.rawToSnow()

    # co["minH":maxH] -> 31 processed heights out of the 32 raw ones
    assert processedSpec.Ze.shape == (5, 31)
    assert processedSpec.W.shape == (5, 31)
    # noise-only input must never produce NaN/Inf in the output moments,
    # only the missingNumber (-9999) fill value where no peak was found
    assert np.all(np.isfinite(processedSpec.Ze))
    assert np.all(np.isfinite(processedSpec.W))

    outFile = tmp_path / "noise_only.nc"
    processedSpec.writeNetCDF(str(outFile))

    with netCDF4.Dataset(outFile) as ds:
        assert ds.variables["Ze"].shape == (5, 31)
        for var in ("W", "spectralWidth", "quality", "height", "time"):
            assert var in ds.variables


def test_pipeline_recovers_injected_peak_velocity(peak_raw_file):
    path, peak_center_bin, peak_heights = peak_raw_file
    raw = IMProToo.mrrRawData(str(path))
    processedSpec = IMProToo.MrrZe(raw)
    processedSpec.rawToSnow()

    expected_velocity = peak_center_bin * processedSpec.co["nyqVdelta"]

    # Ze/W are indexed by height - co["minH"], so shift the raw height index
    minH = processedSpec.co["minH"]
    for h in peak_heights:
        col = h - minH
        if col < 0 or col >= processedSpec.W.shape[1]:
            continue
        w = processedSpec.W[:, col]
        ze = processedSpec.Ze[:, col]
        # at least some timesteps at each injected height should show a
        # recovered peak close to the known Doppler velocity
        assert np.any(np.abs(w - expected_velocity) < 0.5)
        assert np.any(ze > -9999)
