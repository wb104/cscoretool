import numpy

def readCscoreFile(cscoreFile, windowSize):
  # assumes that windowSize is correct (could deduce from high-low for first record for each chromosome ...)
  # could do extrapolation / interpolation if file windowSize was different
  
  cscoreDict = {}
  
  with open(cscoreFile, 'rU') as fp:
    for line in fp:
      line = line.strip()
      fields = line.split()
      if len(fields) == 5:
        (chromo, low, high, label, cscore) = fields
      else:
        (chromo, low, high, cscore) = fields
        
      low = int(low)
      high = int(high)
      cscore = float(cscore)
      
      if chromo not in cscoreDict:
        cscoreDict[chromo] = []
        
      cscoreList = cscoreDict[chromo]
      ind = low // windowSize
      n = ind - len(cscoreList)
      if n > 0:
        cscoreList.extend(n*[0.0])
      cscoreList.append(cscore)
        
  for chromo in cscoreDict:
    cscoreDict[chromo] = numpy.array(cscoreDict[chromo])
    
  return cscoreDict
  
def writeCscoreBedFiles(outputFilePrefix, cscoreDict, windowSize, threshold):
                
  outputCscoreFile = '%s_cscore.bed' % outputFilePrefix
  outputAFile = '%s_A.bed' % outputFilePrefix
  outputBFile = '%s_B.bed' % outputFilePrefix

  with open(outputCscoreFile, 'w') as fpc, open(outputAFile, 'w') as fpa, open(outputBFile, 'w') as fpb:
    for chromo in sorted(cscoreDict):
      cscoreArray = cscoreDict[chromo]
      for n, cscore in enumerate(cscoreArray):
        start = n * windowSize
        end = start + windowSize
        fpc.write('{}\t{}\t{}\t{}\t{}\n'.format(chromo, start, end, n, cscore))
        if cscore > threshold:
          fpa.write('{}\t{}\t{}\t{}\t{}\n'.format(chromo, start, end, n, 1))
        elif cscore < -threshold:
          fpb.write('{}\t{}\t{}\t{}\t{}\n'.format(chromo, start, end, n, 1))
  
  