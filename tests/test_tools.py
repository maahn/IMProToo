import datetime

import numpy as np
import pytest

from IMProToo.tools import (
    date2unix,
    unix2date,
    quantile,
    oneD2twoD,
    limitMaInidces,
    readQualityFile,
    getManualQualityArray,
)


def test_date2unix_unix2date_roundtrip():
    date = datetime.datetime(2021, 6, 15, 12, 30, 45)
    unix = date2unix(date)
    assert unix2date(unix) == date


def test_unix2date_is_naive_utc():
    # utcfromtimestamp-equivalent: naive datetime, no tzinfo attached
    assert unix2date(0).tzinfo is None
    assert unix2date(0) == datetime.datetime(1970, 1, 1, 0, 0, 0)


def test_quantile_matches_numpy_percentile_type7():
    x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    for q in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
        assert quantile(x, q, qtype=7) == pytest.approx(np.percentile(x, q * 100))


def test_quantile_invalid_qtype_returns_none():
    assert quantile([1, 2, 3], 0.5, qtype=10) is None


def test_oneD2twoD_axis0():
    vector = np.array([1, 2, 3])
    result = oneD2twoD(vector, 4, axis=0)
    assert result.shape == (4, 3)
    for row in result:
        assert np.array_equal(row, vector)


def test_oneD2twoD_axis1():
    vector = np.array([1, 2, 3])
    result = oneD2twoD(vector, 4, axis=1)
    assert result.shape == (3, 4)
    for col in result.T:
        assert np.array_equal(col, vector)


def test_oneD2twoD_invalid_axis():
    with pytest.raises(ValueError):
        oneD2twoD(np.array([1, 2, 3]), 4, axis=2)


def test_limitMaInidces_wraps_out_of_range_indices():
    iArray = np.array([-2, -1, 0, 3, 4, 5])
    result = limitMaInidces(iArray, max=4)
    # -2 -> 2, -1 -> 3, 0 -> 0, 3 -> 3, 4 -> 0, 5 -> 1
    assert np.array_equal(result, [2, 3, 0, 3, 0, 1])


def test_limitMaInidces_preserves_shape():
    iArray = np.array([[0, 1], [2, 3]])
    result = limitMaInidces(iArray, max=4)
    assert result.shape == (2, 2)


def test_readQualityFile_parses_multi_entry_dates(tmp_path):
    qualityFile = tmp_path / "quality.txt"
    qualityFile.write_text(
        "120105   3 11.00 17.00 snow on dish\n"
        "            19.00 20.00 snow on dish\n"
        "            22.00 24.00 interference\n"
        "120106   1 00.00 21.00 snow on dish\n"
        "#this is a comment\n"
        "120107   1 00.00 24.00 maintenance\n"
    )

    startTimes, endTimes, comments = readQualityFile(str(qualityFile))

    assert len(startTimes) == len(endTimes) == len(comments) == 5
    assert startTimes[0] == datetime.datetime(2012, 1, 5, 11, 0)
    assert endTimes[0] == datetime.datetime(2012, 1, 5, 17, 0)
    # a "24.00" end time rolls over into 00.00 of the following day
    assert endTimes[2] == datetime.datetime(2012, 1, 6, 0, 0)
    assert endTimes[-1] == datetime.datetime(2012, 1, 8, 0, 0)


def test_getManualQualityArray_flags_time_ranges(tmp_path):
    qualityFile = tmp_path / "quality.txt"
    qualityFile.write_text("120105   1 12.00 13.00 snow on dish\n")

    inside = date2unix(datetime.datetime(2012, 1, 5, 12, 30))
    outside = date2unix(datetime.datetime(2012, 1, 5, 14, 0))
    timeVector = np.array([outside, inside])

    quality = getManualQualityArray(str(qualityFile), timeVector)

    assert list(quality) == [False, True]
