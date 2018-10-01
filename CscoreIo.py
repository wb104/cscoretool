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

def _shouldFlipCscoreSign(cscoreArray, baseCscoreArray):

  n = min(len(cscoreArray), len(baseCscoreArray))
  return numpy.sum(cscoreArray[:n] * baseCscoreArray[:n]) < 0
  
def writeCscoreFile(cscoreFile, cscoreDict, windowSize, baseCscoreDict=None):
  
  if baseCscoreDict is None:
    baseCscoreDict = {}
    
  with open(cscoreFile, 'w') as fp:
    for chromo in sorted(cscoreDict):
      cscoreArray = cscoreDict[chromo]
      if chromo in baseCscoreDict:
        baseCscoreArray = baseCscoreDict[chromo]
        if _shouldFlipCscoreSign(cscoreArray, baseCscoreArray):
          print('Flipping cscore sign for chromo %s' % chromo)
          cscoreDict[chromo] = cscoreArray = -cscoreArray
      for n, cscore in enumerate(cscoreArray):
        fp.write('{}\t{}\t{}\t{}\t{}\n'.format(chromo, n*windowSize, (n+1)*windowSize, n, cscore))
  
def main(cscoreFiles, baseCscoreFile, windowSize):
  
  baseCscoreDict = readCscoreFile(baseCscoreFile, windowSize)
  for cscoreFile in cscoreFiles:
    cscoreDict = readCscoreFile(cscoreFile, windowSize)
    outCscoreFile = '%s_swap.bed' % cscoreFile[:-4]
    writeCscoreFile(outCscoreFile, cscoreDict, windowSize, baseCscoreDict)

if __name__ == '__main__':

  import sys

  if len(sys.argv) < 4:
    print('Arguments: <baseCscore.bed> <windowSize> [ <cscore.bed> ]+')
    sys.exit()

  baseCscoreFile = sys.argv[1]
  windowSize = int(sys.argv[2])
  cscoreFiles = sys.argv[3:]

  main(cscoreFiles, baseCscoreFile, windowSize)
  