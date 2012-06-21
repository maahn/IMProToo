# -*- coding: utf-8 -*-
'''
Copyright (C) 2011,2012 Maximilian Maahn, IGMK (mmaahn@meteo.uni-koeln.de)
make quicklooks from IMProToo NetCDF files.


use: python batch_makeQuicklooks.py pathIn pathOut site

requires:

numpy, matplotlib, netcdf4-python or python-netcdf

'''




import sys
import numpy as np
import glob
import calendar
import datetime
import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc,ticker
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager 
import IMProToo
from IMProTooTools import *
try:
  import netCDF4 as nc
  pyNc = True
except:
  import Scientific.IO.NetCDF as nc
  pyNc = False



def unix2timestamp(unix):
  return datetime.datetime.utcfromtimestamp(unix).strftime("%Y%m%d")
def timestamp2unix(timestamp):
  return calendar.timegm(datetime.datetime(year = int(timestamp[0:4]), month = int(timestamp[4:6]), day = int(timestamp[6:8]), hour = 0, minute = 0, second = 0).timetuple())

def quicklook(site,ncFile,imgFile,imgTitle): 
  """
  Makes Quicklooks of MRR data

  
  @parameter site (str): code for the site where the data was recorded (usually 3 letter)
  @parameter ncFile (str): netcdf file name incl. path, usually "path/mrr_site_yyyymmdd.nc" 
  @parameter imgFile (str): image file name, incl. path, extensions determines file format (e.g. png, eps, pdf ...)
  @parameter imgTitle (str): plot title
  """
  print "##### " + imgTitle + "######"
  if pyNc: ncData = nc.Dataset(ncFile,'r')
  else: ncData  = nc.NetCDFFile(ncFile,'r')

  timestampsNew = ncData.variables["time"][:]
  HNew = ncData.variables["height"][:]
  ZeNew = ncData.variables["Ze"][:]
  noiseAveNew = ncData.variables["etaNoiseAve"][:]
  noiseStdNew = ncData.variables["etaNoiseStd"][:]
  spectralWidthNew = ncData.variables["spectralWidth"][:]
  WNew = ncData.variables["W"][:]
  qualityNew = ncData.variables["quality"][:]
  
  ncData.close()
  #
  date = unix2timestamp(timestampsNew[0])
  starttime = timestamp2unix(date)
  endtime = starttime+60*60*24

  HNew[np.isnan(HNew)] = -9999
  ylim = [np.min(HNew[HNew!=-9999]),np.max(HNew)]
  xlim = [starttime,endtime]
  timestampsNew = oneD2twoD(timestampsNew,ZeNew.shape[1],1)

  fig=plt.figure(figsize=(10, 13))
  
  sp1 = fig.add_subplot(511)
  sp1.set_title(imgTitle)
  levels = np.arange(-15,40,0.1)
  plotCF = sp1.contourf(timestampsNew,HNew, ZeNew, levels,cmap=plt.get_cmap("spectral"), extend="both")#
  cbZe=plt.colorbar(plotCF)
  cbZe.set_label('MRR Ze [dBz]')
  sp1.set_ylim(ylim)
  sp1.set_xlim(xlim)

  sp1.axhline(HNew[-1,2])
  sp1.axhline(HNew[-1,29])

  sp2 = fig.add_subplot(512)
  levels = np.arange(-10,18,0.1)
  plotCF = sp2.contourf(timestampsNew,HNew, WNew, levels,cmap=plt.get_cmap("spectral"), extend="both")#
  cbZe=plt.colorbar(plotCF)
  cbZe.set_label('MRR W [m/s]')
  sp2.set_ylim(ylim)
  sp2.set_xlim(xlim)

  sp2.axhline(HNew[-1,2])
  sp2.axhline(HNew[-1,29])
  
  sp3 = fig.add_subplot(513)
  levels = np.arange(0,1.5,0.1)
  plotCF = sp3.contourf(timestampsNew,HNew, spectralWidthNew, levels,cmap=plt.get_cmap("spectral"), extend="both")#
  cbZe=plt.colorbar(plotCF)
  cbZe.set_label('spectralWidth [m/s]')
  sp3.set_ylim(ylim)
  sp3.set_xlim(xlim)

  sp3.axhline(HNew[-1,2])
  sp3.axhline(HNew[-1,29])

  sp4 = fig.add_subplot(514)
  levels = np.arange(1e-10,1e-8,2e-10)
  plotCF = sp4.contourf(timestampsNew,HNew, noiseAveNew, levels,cmap=plt.get_cmap("spectral"), extend="both")#
  cbZe=plt.colorbar(plotCF)
  cbZe.set_label('mean spectral noise [1/m]')
  sp4.set_ylim(ylim)
  sp4.set_xlim(xlim)
  sp4.axhline(HNew[-1,2])
  sp4.axhline(HNew[-1,29])  
  #import pdb;pdb.set_trace()
    

  sp5 = fig.add_subplot(515)
  levels = np.arange(20)
  for i in levels:
    levels[i] = 2**i
  plotCF = sp5.contourf(timestampsNew,HNew, qualityNew, levels,cmap=plt.get_cmap("spectral"), norm = matplotlib.colors.LogNorm())#
  cbZe=plt.colorbar(plotCF)
  cbZe.set_label('quality array')
  sp5.set_ylim(ylim)
  sp5.set_xlim(xlim)
  sp5.axhline(HNew[-1,2])
  sp5.axhline(HNew[-1,29])  
    
    
  #sp1.set_xlim(np.min(timestampsNew),np.max(timestampsNew))
  sp1.set_xticks(np.arange(sp1.get_xlim()[0],sp1.get_xlim()[1],7200))
  sp1.set_xticklabels([])
  
  #sp2.set_xlim(np.min(timestampsNew),np.max(timestampsNew))
  sp2.set_xticks(np.arange(sp1.get_xlim()[0],sp1.get_xlim()[1],7200))
  sp2.set_xticklabels([])

  #sp3.set_xlim(np.min(timestampsNew),np.max(timestampsNew))
  sp3.set_xticks(np.arange(sp1.get_xlim()[0],sp1.get_xlim()[1],7200))
  sp3.set_xticklabels([])

  #sp4.set_xlim(np.min(timestampsNew),np.max(timestampsNew))
  sp4.set_xticks(np.arange(sp1.get_xlim()[0],sp1.get_xlim()[1],7200))
  sp4.set_xticklabels([])

    #pdb.set_trace()
  #sp5.set_xlim(np.min(timestampsNew)-60,np.max(timestampsNew))
  sp5.set_xticks(np.arange(sp5.get_xlim()[0],sp5.get_xlim()[1]+7200,7200))
  niceDates = list()
  for timestamp in np.arange(sp5.get_xlim()[0],sp5.get_xlim()[1]+7200,7200):
    niceDates.append(str(datetime.datetime.utcfromtimestamp(timestamp).strftime("%H:%M")))
  sp5.set_xticklabels(niceDates)

    



  plt.subplots_adjust(hspace=0.02,left=0.085,right=0.78)
    
  plt.savefig(imgFile)
  print(imgFile)

  plt.close()
  return 

  

if len(sys.argv) < 4:
  print 'use: batch_makeQuicklooks.py pathIn pathOut site'
  sys.exit()

pathIn = sys.argv[1]
pathOut = sys.argv[2]  
site = sys.argv[3]


try: os.mkdir(pathOut)
except OSError: pass


for ncFile in np.sort(glob.glob(pathIn+"/*")):
  date = ncFile[-11:-3]
  print date, ncFile
  imgFile = pathOut + "/mrr_improtoo_"+IMProToo.__version__+'_'+site+"_"+date+".png"
  imgTitle = site + " " + date + " IMProToo " + IMProToo.__version__
  quicklook(site,ncFile,imgFile,imgTitle)

