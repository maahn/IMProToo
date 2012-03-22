# -*- coding: utf-8 -*-
'''
ImProToo
Innominate MRR Processing Tool

(c)2011,2012 Max Maahn, IGMK (mmaahn@meteo.uni-koeln.de)

some helper functions

'''

import numpy as np
import time
import datetime
import calendar
import warnings

warnings.filterwarnings('always', '.*', UserWarning,)


def getManualQualityArray(inputFile,timeVector):
  startTimes, endTimes, comments = _readQualityFile(inputFile)
  quality = np.zeros(timeVector.shape,dtype=bool)
  for startTime, endTime in zip(startTimes,endTimes):
    startTime = date2unix(startTime)
    endTime = date2unix(endTime)
    quality = quality + ((timeVector >=startTime)* (timeVector <endTime))
  return quality
  

  
def _readQualityFile(inputFile):
  '''
  reads manual quality control files, format is in ancient 'Hatpro" format:
  
  file format is:
  date, no of entries for this day, starttime (00.00 - 23:59), endtime (00:01 - 24:00), comment
  120105   3 11.00 17.00 snow on dish
	      19.00 20.00 snow on dish
	      22.00 24.00 interference
  120106   1 00.00 21.00 snow on dish
  #this is a comment
  120107   1 00.00 24.00 maintenance
  
  '''
  startTimes = list()
  endTimes = list()
  comments = list()

  belongsToOldDate = 1
  ll = -1

  f = open(inputFile,"r")
  for line in f.readlines():
    ll += 1  
    fields = line.split()
    try:
      #check for comments
      if fields[0][0:1] == '#':
	print 'comment', ' '.join(fields)
	continue
      # does the line hav its own timestamp?
      if belongsToOldDate == 1:
	startDate = fields[0]
	belongsToOldDate = int(fields[1])
	startTime = fields[2]
	endTime = fields[3]
	comment = ' '.join(fields[4:])
      else:
	startTime = fields[0]
	endTime = fields[1]   
	comment = ' '.join(fields[2:])
	belongsToOldDate -= 1
      #yes, 24.00 is a weired timestamp, correct this to 00.00 next day
      if endTime == '24.00':
	#import pdb;pdb.set_trace()
	endDate = datetime.datetime.strftime((datetime.datetime.strptime(startDate,'%y%m%d')+datetime.timedelta(1)),'%y%m%d')
	endTime = '00.00'
      else:
	endDate = startDate
      comment = ' '.join(fields[4:])
      startTimes.append(datetime.datetime.strptime(startDate+' '+startTime,'%y%m%d %H.%M'))
      endTimes.append(datetime.datetime.strptime(endDate+' '+endTime,'%y%m%d %H.%M'))
      comments.append(comment)
    except:
      belongsToOldDate = 1
      warnings.warn('Could not read line no. '+str(ll) + str(fields))

  f.close()
  
  return startTimes, endTimes, comments

    
def date2unix(date):
  return calendar.timegm(date.timetuple())

def unix2date(unix):
  return datetime.datetime.utcfromtimestamp(unix)   
  
  
def quantile(x, q,  qtype = 7, issorted = False):
  """
  Args:
      x - input data
      q - quantile
      qtype - algorithm
      issorted- True if x already sorted.

  Compute quantiles from input array x given q.For median,
  specify q=0.5.

  References:
      http://reference.wolfram.com/mathematica/ref/Quantile.html
      http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile

  Author:
      Ernesto P.Adorio Ph.D.
      UP Extension Program in Pampanga, Clark Field.
  """
  if not issorted:
      y = sorted(x)
  else:
      y = x
  if not (1 <= qtype <= 9): 
      return None  # error!

  # Parameters for the Hyndman and Fan algorithm
  abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
	  (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
	  (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3

	  (0,   0, 0, 1), # California linear interpolation, R type 4
	  (0.5, 0, 0, 1), # hydrologists method, R type 5
	  (0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6 
	  (1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
	  (1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
	  (3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
	  ]

  a, b, c, d = abcd[qtype-1]
  n = len(x)
  g, j = np.modf( a + (n+b) * q -1)
  if j < 0:
      return y[0]
  elif j >= n:           
      return y[n-1]   # oct. 8, 2010 y[n]???!! uncaught  off by 1 error!!!

  j = int(np.floor(j))
  if g ==  0:
      return y[j]
  else:
      return y[j] + (y[j+1]- y[j])* (c + d * g)    
