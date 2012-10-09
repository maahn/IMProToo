# -*- coding: utf-8 -*-
#Copyright (C) 2011,2012 Maximilian Maahn, IGMK (mmaahn@meteo.uni-koeln.de)

#example script for converting mrrRaw data to netcdf using IMProToos

import sys
import numpy as np
import glob
import os
import IMProToo


version = IMProToo.__version__


if len(sys.argv) < 4:
  sys.exit('use: python batch_convert_rawData.py pathIn pathOut site')

pathIn = sys.argv[1]
pathOut = sys.argv[2]
site =  sys.argv[3]


print pathIn

try: os.mkdir(pathOut)
except OSError: pass

#go through all gz compressed files in pathIn/year/month/
for nfile in np.sort(glob.glob(pathIn+"/*/*/*mrr*raw*")):
  date = nfile[-15:-7] #expects filename as xxxxxxxxxxx_date.mrr.gz
  
  fileOut = pathOut+"/mrr_improtoo_"+version+"_"+site+"_"+date+".nc"
  print date, nfile, fileOut
  
  #load raw data from file
  print "reading...",nfile
  try: rawData = IMProToo.mrrRawData(nfile)
  except: 
    print "print could not read data"
    continue
  
  #convert rawData object
  processedSpec = IMProToo.MrrZe(rawData)
  #average rawData to 60s
  processedSpec.averageSpectra(60)
  #the MRR at 'lyr' was affected by interference for some days, dealiasing routine needs to know about that:
  if site == "lyr" and date in ['20100620','20100621','20100622', '20100623', '20100624', '20100625', '20100626', '20100627','20100628','20100629','20100630','20100701','20100702','20100703','20100704','20100705','20100706', '20100707']:
    processedSpec.co['dealiaseSpectrum_heightsWithInterference'] =  processedSpec.co['dealiaseSpectrum_heightsWithInterference'] + [25,26,27,28,29,30]
  #creator attribute of netCDF file
  processedSpec.co["ncCreator"] = "M.Maahn, IGM University of Cologne"    
    
    
  #calculate Ze and other moments
  processedSpec.rawToSnow()
  
  #write all variables to a netCDF file.
  print "writing...",fileOut
  processedSpec.writeNetCDF(fileOut,ncForm="NETCDF3_CLASSIC")
