    self.co["check4BorderPeaks"] = False # do not use with newer hildebrands!!
    self.co["simpleDealiasingMinPeak"] = 64-25
    self.co["simpleDealiasingMaxPeak"] = 128-25    
    self.co["useHildebrandExtraLimit"] = True  #use a two step hildebrand: loop is stooped if hildebrandExtraLimit is reached. one more bin is added if above hildebrand limit.

    
    #old hildebrand
      def _getPeakHildebrand(self,dataFlat,iMax,noSpecs,h):
    """
    get the peak of the spectrum using Hildebrand algorithm.
    
    @parameter dataFlat (numpy float64): flat spectrum from MRR Raw data
    @parameter iMax (numpy float64): vector containing indices of the maxima
    @parameter Nspec (numpy float64): number of spectra of each averaged spectrum

    @return - iPeakMin, iMax (int float64): edges of each spectrum
    """
    
    # first get the limit reflectivity
    limits = self._noiseHildebrand(dataFlat,noSpecs,h)
    maskHildebrand = np.ones(np.shape(dataFlat),dtype=bool)
    iPeakMax = deepcopy(iMax)
    iPeakMin = deepcopy(iMax)
    
    #import pdb; pdb.set_trace()
    
    if self.co["useHildebrandExtraLimit"] == True:
      #not only uses extra limit, but also starts at the peak!, thus specturm is refolded around peak!
      
      #then get the edges of the peak as index of the spectrum
      for k in np.arange(iMax.shape[0]):
	#unmask the peak
	maskHildebrand[k,iMax[k]] = False
	
	spectrum = np.roll(dataFlat[k],-iMax[k])
	mask = np.roll(maskHildebrand[k],-iMax[k])
	#to the right
	for i in np.arange(1,dataFlat.shape[-1],1):
	  #unmask if above limit (=peak)
	  if spectrum[i]>limits[k]*self.co["getPeak_hildeExtraLimit"]:
	    mask[i] = False
	  #else stop
	  else:
	    #unmask on last bin if between limits[k]*self.co["getPeak_hildeExtraLimit"] and limits[k], but stop in any case!
	    if spectrum[i]>limits[k]: 
	      mask[i] = False
	    break
	#to the left
	for i in np.arange(dataFlat.shape[-1]-1,0-1,-1):
	  if spectrum[i]>limits[k]*self.co["getPeak_hildeExtraLimit"]:
	    mask[i] = False
	  else:
	    if spectrum[i]>limits[k]: 
	      mask[i] = False
	    break

	dataFlat[k] = np.roll(spectrum,iMax[k])
	maskHildebrand[k] = np.roll(mask,iMax[k])

    else:
      #then get the edges of the peak as index of the spectrum
      for k in np.arange(iMax.shape[0]):
	#unmask the peak
	maskHildebrand[k,iMax[k]] = False
	#to the right
	for i in np.arange(iMax[k]+1,dataFlat.shape[-1],1):
	  if dataFlat[k,i]>limits[k]:
	    maskHildebrand[k,i] = False
	  else:
	    break
	#to the left
	for i in np.arange(iMax[k]-1,0-1,-1):
	  if dataFlat[k,i]>limits[k]:
	    maskHildebrand[k,i] = False
	  else:
	    break
	    
	    
    return maskHildebrand
 
     
      
    #if self.co["dealiaseSpectrum"] == 'advanced':
      #shiftedSpecIndex = np.zeros(self._shape3D,dtype=int)
      #shiftedSpecIndex[:] = np.arange(self.co["widthSpectrum"])+self.co["widthSpectrum"]
      #self.specIndex = np.arange(3*self.no_v)
      
      #self.no_v = self.no_v * 3
      ##who applies for dealiasing?
      #extendedRawSpectrum = deepcopy(self.rawSpectrum.data)
      #extendedRawSpectrum = np.concatenate((np.roll(extendedRawSpectrum,1,axis=1),extendedRawSpectrum,np.roll(extendedRawSpectrum,-1,axis=1)),axis=2)
      ##do not apply fo first range gates
      #extendedRawSpectrum[:,0,:self.co["widthSpectrum"]] = 0
      ##and not to the last one
      #extendedRawSpectrum[:,self.co["noH"]-1,2*self.co["widthSpectrum"]:]=0
      #extendedRawSpectrum = np.ma.masked_array(extendedRawSpectrum,True)
      #self.specVel = np.array(list(self.co["nyqVel"] - self.co["widthSpectrum"]*self.co["nyqVdelta"] )+list(self.co["nyqVel"])+list(self.co["nyqVel"] + self.co["widthSpectrum"]*self.co["nyqVdelta"] ))
      #self.specVel3D = np.zeros(np.shape(extendedRawSpectrum))
      #self.specVel3D[:] = self.specVel
      
      ##find folding candidates and fill the gap around 0 m/s
      #foldingCandidates = np.zeros((np.shape(self.rawSpectrum)[0],np.shape(self.rawSpectrum)[1],),dtype=bool)
      #for h in xrange(2,self.co["noH"]):
	#foldingCandidates[:,h] = (self.rawSpectrum.mask[:,h-1,self.co["spectrumBorderMax"][h-1]-1] == False) * (self.rawSpectrum.mask[:,h,self.co["spectrumBorderMin"][h]] == False)
	##+ (spectrum.mask[:,h-1,spectrumMaxs[h-1]-2] == False) * (spectrum.mask[:,h,spectrumMins[h]] == False) + (spectrum.mask[:,h-1,spectrumMaxs[h-1]-1] == False) * (spectrum.mask[:,h,spectrumMins[h]+1] == False)
	#self.rawSpectrum.mask[:,h,0:self.co["spectrumBorderMin"][h]][foldingCandidates[:,h]] = False
	#self.rawSpectrum.mask[:,h-1,self.co["spectrumBorderMax"][h-1]:][foldingCandidates[:,h]] = False
      #self.quality[foldingCandidates] = self.quality[foldingCandidates] + 0b0000000010000000 #9) suspicous of dealiasing
      #deAlSpectrum,shiftedSpecIndex,self.quality,self._foldDirection = self._dealiaseSpectrumAdvanced_OLD_TBD(self.rawSpectrum,shiftedSpecIndex,self.quality,foldingCandidates)

      #if self.co["advancedDealiasingCheck4Coherence"]:
	  #self.rawSpectrum,shiftedSpecIndex,self.quality = self._deAlCoherence(deAlSpectrum,shiftedSpecIndex,self.rawSpectrum,self.quality,self._foldDirection,foldingCandidates)
      #else:
	#self.rawSpectrum = deAlSpectrum
    ##reshape the raw spectrum to 192 bins with the most likely peak marked
    ##import pdb; pdb.set_trace()
      #for t in np.arange(self._shape2D[0]):
	#for h in np.arange(self._shape2D[1]):
	  #peakIndices = shiftedSpecIndex[t,h][~self.rawSpectrum.mask[t,h]]
	  #extendedRawSpectrum.mask[t,h][peakIndices] = False
      #self.rawSpectrum = extendedRawSpectrum
      
    ##import pdb; pdb.set_trace()
    
    ###mask heights 0,1,30 (actual numbers 1,2,31, but zero heigth is not processed)
    ##self.rawSpectrum.mask[:,0:2] = True
    ##self.rawSpectrum.mask[:,self.co["noH"]-1:] = True

  def _getBorderPeaks(self,spectrumFlat,joinedMask):
    '''
    find all spectra with signal at the borders
    '''
    foldingCandidatesLeft = (joinedMask[...,0] == False)
    foldingCandidatesRight = (joinedMask[...,-1] == False)
    
    
    if np.any(foldingCandidatesLeft):
      
      #the getPeakDescendingAve function needs a list of maxima, so give it the right border
      second_iMaxLeft= np.array([spectrumFlat.shape[1]-1]*np.sum(foldingCandidatesLeft))
      
      #get the additional mask, the already found peak is masked out and not treated!
      additionalMaskLeft = self._getPeakDescendingAve(np.ma.masked_array(spectrumFlat[foldingCandidatesLeft],~joinedMask[foldingCandidatesLeft]),second_iMaxLeft)
      #cancel additional peak which are only one bin wide!
      additionalMaskLeft[np.sum(~additionalMaskLeft,axis=-1)<=1] = True
      #add the additional mask to exisiting one
      joinedMask[foldingCandidatesLeft] = joinedMask[foldingCandidatesLeft] * additionalMaskLeft
      
    if np.any(foldingCandidatesRight):  
      #same for the otehr border
      second_iMaxRight= np.array([0]*np.sum(foldingCandidatesRight))

      
      additionalMaskRight =self._getPeakDescendingAve(np.ma.masked_array(spectrumFlat[foldingCandidatesRight],~joinedMask[foldingCandidatesRight]),second_iMaxRight)
      additionalMaskRight[np.sum(~additionalMaskRight,axis=-1)<=1] = True
      
      joinedMask[foldingCandidatesRight] = joinedMask[foldingCandidatesRight] * additionalMaskRight
   
    return joinedMask    
  def _dealiaseSpectrumAdvanced_OLD_TBD(self,newSpectrum,shiftedSpecIndex,qualityArray,foldingCandidates):
  	'''
  	#newSpectrum = deepcopy(rawSpectrum)
  	#newVelocities = deepcopy(specVel)
  
  	#get according velocities
  	'''
  	shiftedIndices = np.arange(0,self.co["widthSpectrum"]*3)
  	checkForDealaiasing = np.where(np.sum(foldingCandidates,axis=-1)>=1)[0]
  
  	foldDirection = np.zeros(np.shape(foldingCandidates),dtype=int)
  
  	noSpecHalf = self.co["widthSpectrum"] / 2
  
  	for t in checkForDealaiasing:
  	  #get all spectra from one time step
  	  singleSpectrum = deepcopy(newSpectrum[t])
  	  #don't do anything if everything is masked!'
  	  if np.all(singleSpectrum.mask == True): 
  	if self.co["debug"] > 1: print t, "all masked"
  	continue
  	  if self.co["debug"] > 0: print t
  	  #find all maxima in the spectrum
  	  iMax = np.ma.masked_array(np.argmax(singleSpectrum,axis=1),np.all(singleSpectrum.mask,axis=-1))
  	  #decide which is the most trustfull one, the one closest to rangegate x
  	  trustfullHeight = np.ma.argmin(np.abs(iMax-20))
  	  #get the shifted index of the max of the trusted height
  	  lastMax = iMax[trustfullHeight]+self.co["widthSpectrum"]
  	  for h in range(trustfullHeight,self.co["noH"]-1):
  	#add lower and upper spectra, according velocities can be found in threeVelocities
  	threeSpectra = np.ma.concatenate((singleSpectrum[h-1],singleSpectrum[h],singleSpectrum[h+1],))
  	# find the most likely peak of the 3 spectra on the basis of the height before
  	# don't do anything if there is no peak
  	if np.all(threeSpectra.mask[lastMax-noSpecHalf:lastMax+noSpecHalf] == True):
  	  if self.co["debug"] > 1: print t,h,"np.all(threeSpectra.mask[lastMax-noSpecHalf:lastMax+noSpecHalf] == True)"
  	  newSpectrum[t,h].mask = True
  	  #self.co["debug"] obnly: newSpectrum[t,h].mask = np.logical_not(newSpectrum[t,h].mask)
  	  continue
  	#get borders and width of the most likely peak
  	peakMinIndex = np.ma.min(np.where(threeSpectra[lastMax-noSpecHalf:lastMax+noSpecHalf])) + lastMax - noSpecHalf
  	peakMaxIndex = np.ma.max(np.where(threeSpectra[lastMax-noSpecHalf:lastMax+noSpecHalf])) + lastMax - noSpecHalf
  	peakWidth = peakMaxIndex - peakMinIndex
  	#how much spectrum do we have to add to the left and to the right?
  	addLeft = (self.co["widthSpectrum"]-peakWidth)/2
  	addRight = (self.co["widthSpectrum"]-peakWidth)/2
  	if (peakWidth % 2 == 1) : addRight +=1
  	#find the new borders, equaly around the peak
  	leftBorder = peakMinIndex-addLeft
  	rightBorder = peakMaxIndex+addRight
  	#what if borders larger than spectrum?
  	if leftBorder < 0:
  	  leftBorder = 0
  	  rightBorder = self.co["widthSpectrum"]
  	elif rightBorder > 3*self.co["widthSpectrum"]:
  	  leftBorder = 2*self.co["widthSpectrum"]
  	  rightBorder = 3*self.co["widthSpectrum"]
  	#finally get the new spectrum and the according velocities!
  	newSpectrum[t,h] = threeSpectra[leftBorder:rightBorder]
  	shiftedSpecIndex[t,h] = shiftedIndices[leftBorder:rightBorder]
  
  	#for statistics, in which direction did we fold
  	if leftBorder< self.co["widthSpectrum"]:
  	  foldDirection[t,h] = -1
  	elif rightBorder >=128:
  	  foldDirection[t,h] = 1    
  	#and remember the max for the next height
  	lastMax = np.argmax(newSpectrum[t,h])+leftBorder
  	if self.co["debug"] > 1: print t,h,leftBorder,lastMax
  
  	  lastMax = iMax[trustfullHeight]+self.co["widthSpectrum"]
  	  for h in range(trustfullHeight-1,1,-1):
  	#if iMax.mask[h] == True:# or foldingCandidates[t,h] == False:
  	  #if self.co["debug"] > 1: print t,h,"iMax.mask[h] == True"#" or foldingCandidates[t,h] == False"
  	  #continue
  	#print t,h
  	threeSpectra = np.ma.concatenate((singleSpectrum[h-1],singleSpectrum[h],singleSpectrum[h+1],))
  
  	#we don't want a shrinked mask, it has to be full size!
  	if threeSpectra.mask.shape == ():
  	  threeSpectra.mask = np.ones(threeSpectra.shape,dtype=bool)*threeSpectra.mask
  
  	if np.all(threeSpectra.mask[lastMax-noSpecHalf:lastMax+noSpecHalf] == True):
  	  if self.co["debug"] > 1: print t,h,"np.all(threeSpectra.mask[lastMax-noSpecHalf:lastMax+noSpecHalf] == True)"
  	  #self.co["debug"] only newSpectrum[t,h].mask = np.logical_not(newSpectrum[t,h].mask)
  	  newSpectrum[t,h].mask = True
  	  continue
  	peakMinIndex = np.ma.min(np.where(threeSpectra[lastMax-noSpecHalf:lastMax+noSpecHalf])) + lastMax - noSpecHalf
  	peakMaxIndex = np.ma.max(np.where(threeSpectra[lastMax-noSpecHalf:lastMax+noSpecHalf])) + lastMax - noSpecHalf
  	peakWidth = peakMaxIndex - peakMinIndex
  	addLeft = (self.co["widthSpectrum"]-peakWidth)/2
  	addRight = (self.co["widthSpectrum"]-peakWidth)/2
  	if (peakWidth % 2 == 1) : addRight +=1
  	leftBorder = peakMinIndex-addLeft
  	rightBorder = peakMaxIndex+addRight
  	if leftBorder < 0:
  	  leftBorder = 0
  	  rightBorder = self.co["widthSpectrum"]
  	elif rightBorder > 3*self.co["widthSpectrum"]:
  	  leftBorder = 2*self.co["widthSpectrum"]
  	  rightBorder = 3*self.co["widthSpectrum"]
  	newSpectrum[t,h] = threeSpectra[leftBorder:rightBorder]
  	shiftedSpecIndex[t,h] = shiftedIndices[leftBorder:rightBorder]
  	#for statistics, in which direction did we fold
  	if leftBorder< self.co["widthSpectrum"]:
  	  foldDirection[t,h] = -1
  	elif rightBorder >=128:
  	  foldDirection[t,h] = 1    
  	lastMax = np.argmax(newSpectrum[t,h])+leftBorder
  	if self.co["debug"] > 1: print t,h,leftBorder,lastMax
  
  	qualityArray[foldDirection == 1] = qualityArray[foldDirection == 1] +   0b0000000000100000 #11) spectrum dealiased to the right 
  	qualityArray[foldDirection == -1] = qualityArray[foldDirection == -1] + 0b0000000001000000 #10) spectrum dealiased to the left 
  
  	return newSpectrum,shiftedSpecIndex,qualityArray,foldDirection
  
  def _deAlCoherence_OLD_TBD(self,newSpectrum,newShiftedIndex,oldSpectrum,qualityArray,foldDirection,foldingCandidates):
  	'''
  	make sure no weired foldings happend
  	'''
  	#checkForWeiredFolding = (np.all(foldDirection>=0,axis=-1) + np.all(foldDirection<=0,axis=-1)) * (np.sum(foldingCandidates,axis=-1)>1)
  	checkForWeiredFolding = np.sum(foldingCandidates,axis=-1)>1
  	foldDirectionPerTimeStep = [cmp(i,0) for i in np.sum(foldDirection,axis=-1)]
  	lastFolding = 0
  	#import pdb; pdb.set_trace()
  	for t in np.arange(np.shape(oldSpectrum)[0]):
  	  if self.co["debug"] > 4: print t,checkForWeiredFolding[t],foldDirectionPerTimeStep[t],lastFolding
  	  if checkForWeiredFolding[t] == False or foldDirectionPerTimeStep[t] == 0:
  	lastFolding = 0
  	continue
  
  	  if lastFolding != foldDirectionPerTimeStep[t] and lastFolding!= 0:
  	if self.co["debug"] > 1: print "ALERT",t,"crazy fold!",foldDirectionPerTimeStep[t]
  
  	#repair here
  	#
  	if foldDirectionPerTimeStep[t] == 1:
  	  iToMove = np.where(foldDirection[t] == 1)[0]
  	  newSpectrum[t,iToMove + 1] = newSpectrum[t,iToMove]
  	  newSpectrum[t,iToMove[0]] = oldSpectrum[t,iToMove[0]]
  	  newShiftedIndex[t,iToMove + 1] = newShiftedIndex[t,iToMove] -  self.co["widthSpectrum"]
  	  newShiftedIndex[t,iToMove[0]] = np.arange(self.co["widthSpectrum"]) + self.co["widthSpectrum"]
  	  qualityArray[t] = qualityArray[t] + 0b0000000000010000 #12) spectrum refolded due to coherence test!
  
  
  	elif foldDirectionPerTimeStep[t] == -1:
  	  iToMove = np.where(foldDirection[t] == -1)[0]
  	  newSpectrum[t,iToMove -1 ] = newSpectrum[t,iToMove]
  	  newSpectrum[t,iToMove[-1]] = oldSpectrum[t,iToMove[-1]]
  	  newShiftedIndex[t,iToMove -1 ] = newShiftedIndex[t,iToMove] +  self.co["widthSpectrum"]
  	  newShiftedIndex[t,iToMove[-1]] = np.arange(self.co["widthSpectrum"]) + self.co["widthSpectrum"]
  	  self.quality[t] = qualityArray[t] + 0b0000000000010000 #12) spectrum refolded due to coherence test!
  
  
  	lastFolding = lastFolding
  	  else:
  	lastFolding = foldDirectionPerTimeStep[t]
  
  	return newSpectrum,newShiftedIndex,qualityArray