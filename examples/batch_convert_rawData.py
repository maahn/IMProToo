# -*- coding: utf-8 -*-

'''
Copyright (C) 2011-2021 Maximilian Maahn, U Leipzig
maximilian.maahn_AT_uni-leipzig.de

example script for converting mrrRaw data to netcdf using IMProToos
'''
from __future__ import print_function

import sys
import numpy as np
import glob
import os
import datetime
import IMProToo
import gzip


version = IMProToo.__version__


if len(sys.argv) < 4:
    sys.exit('use: python batch_convert_rawData.py pathIn pathOut site')

pathIn = sys.argv[1]
pathOut = sys.argv[2]
site = sys.argv[3]

skipExisting = True

print(pathIn)

try:
    os.mkdir(pathOut)
except OSError:
    pass

# go through all gz compressed files in pathIn/year/month/
for nfile in np.sort(glob.glob(pathIn+"/*raw*")):
    # get the timestamp
    timestamp = None
    if nfile.split('.')[-1] == 'gz':
        f = gzip.open(nfile, 'rt')
    else:
        f = open(nfile, 'r')
    # Sometimes the first MRR timestamps are from the day before, so we cannot take the first date we found. get list of line breaks
    line_offset = []
    offset = 0
    for line in f:
        line_offset.append(offset)
        offset += len(line)
    f.seek(0)

    # Now, to skip 20% of the file
    f.seek(line_offset[len(line_offset)//5])

    # now find the date
    try:
        while True:
            string = str(f.readline())
            if not string:
                break
            if string[:2] == "T:":
                timestamp = datetime.datetime.strptime(
                    string[2:14], "%y%m%d%H%M%S").strftime("%Y%m%d")
                break
            elif string[:4] == "MRR ":
                timestamp = datetime.datetime.strptime(
                    string[4:16], "%y%m%d%H%M%S").strftime("%Y%m%d")
                break

    finally:
        f.close()

    if timestamp is None:
        print("did not find MRR timesamp in %s, Skipping" % nfile)
        continue

    fileOut = pathOut+"/mrr_improtoo_"+version+"_"+site+"_"+timestamp+".nc"

    if skipExisting and (os.path.isfile(fileOut) or os.path.isfile(fileOut+".gz")):
        print("NetCDF file aready exists, skipping: ", timestamp, nfile, fileOut)
        continue

    print(timestamp, nfile, fileOut)

    # load raw data from file
    print("reading...", nfile)
    try:
        rawData = IMProToo.mrrRawData(nfile)
    except:
        print("could not read data")
        continue

    try:
        # convert rawData object
        processedSpec = IMProToo.MrrZe(rawData)
        # average rawData to 60s
        processedSpec.averageSpectra(60)
        # the MRR at 'lyr' was affected by interference for some days, dealiasing routine needs to know about that:
        if site == "lyr" and timestamp in ['20100620', '20100621', '20100622', '20100623', '20100624', '20100625', '20100626', '20100627', '20100628', '20100629', '20100630', '20100701', '20100702', '20100703', '20100704', '20100705', '20100706', '20100707']:
            processedSpec.co['dealiaseSpectrum_heightsWithInterference'] = processedSpec.co[
                'dealiaseSpectrum_heightsWithInterference'] + [25, 26, 27, 28, 29, 30]
        # creator attribute of netCDF file
        processedSpec.co["ncCreator"] = "M.Maahn, IGM University of Cologne"

        # calculate Ze and other moments
        processedSpec.rawToSnow()

        # write all variables to a netCDF file.
        print("writing...", fileOut)
        processedSpec.writeNetCDF(fileOut, ncForm="NETCDF3_CLASSIC")
    except Exception as error:
        print(str(error))
        print("could not process data")
        continue
