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
