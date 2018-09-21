import numpy

import time

import CscoreUtil

MAXC = 0.9999  # not sure really need this

def calcBinFromD(d):
  
  steplength=0.04
  return int((numpy.log10(d+0.01)-3) / steplength)

def readChromoSizes(chromoSizeFile):
  
  chromoSizeDict = {}  # chromo --> size
  with open(chromoSizeFile, 'rU') as fp:
    for line in fp:
      line = line.strip()
      fields = line.split()
      if len(fields) == 2:
        chromo, size = fields
        size = int(size)
        chromoSizeDict[chromo] = size
        print('Chromo {} has size {}'.format(chromo, size))

  return chromoSizeDict
  
def readInteractions(inputSummaryFile, chromoAnalysisDict):
  
  chromos = set(chromoAnalysisDict.keys())
  with open(inputSummaryFile, 'rU') as fp:
    for line in fp:
      line = line.strip()
      fields = line.split()
      if len(fields) >= 6:
        junk, chromo1, position1, junk, chromo2, position2 = fields[:6]
        if chromo1 == chromo2 and chromo1 in chromos:
          chromoAnalysis = chromoAnalysisDict[chromo1]
          position1 = int(position1)
          position2 = int(position2)
          CscoreUtil.updateInteractions(interactions, position1, position2, chromoAnalysis.chromoSize,
                                        chromoAnalysis.windowSize, chromoAnalysis.dBinMin)
  
def readNccFile(nccFile, chromoAnalysisDict):
  
  chromos = set(chromoAnalysisDict.keys())
  with open(nccFile, 'rU') as fp:
    for line in fp:
      chr_a, f_start_a, f_end_a, start_a, end_a, strand_a, chr_b, f_start_b, f_end_b, start_b, end_b, strand_b, ambig_group, pair_id, swap_pair = line.split()
 
      if strand_a == '+':
        pos_a = int(f_start_a)
      else:
        pos_a = int(f_end_a)
 
      if strand_b == '+':
        pos_b = int(f_start_b)
      else:
        pos_b = int(f_end_b)
 
      if chr_a == chr_b and chr_a in chromos:
        chromoAnalysis = chromoAnalysisDict[chr_a]
        pos_a = int(pos_a)
        pos_b = int(pos_b)
        CscoreUtil.updateInteractions(chromoAnalysis.interactions, pos_a, pos_b, chromoAnalysis.chromoSize,
                                      chromoAnalysis.windowSize, chromoAnalysis.dBinMin)
  
class ChromoAnalysis:
  
  def __init__(self, chromo, chromoSize, windowSize, dBinMin):
    # chromo = chromosome name
    # chromoSize = number of base pairs in chromosome
    # windowSize = window size (in base pairs) being considered
    # interactions = matrix of interactions
    
    # procedure:
    #  create this object
    #  fill in self.interactions
    #  call run()
    
    self.chromo = chromo
    self.chromoSize = chromoSize
    self.windowSize = windowSize
    self.dBinMin = dBinMin
    
    self.numWindows = 1 + (chromoSize - 1) // windowSize
    self.interactions = numpy.zeros((self.numWindows, self.numWindows), dtype='int32')
    
  def calcBinFromIndex(self, i, j):
  
    d = self.windowSize * abs(j - i)
    dBin = calcBinFromD(d) - self.dBinMin
  
    return dBin
      
  def initData(self):
  
    numpy.random.seed(123456)
    
    self.dBinMax = self.calcBinFromIndex(0, self.numWindows) + self.dBinMin
    
    self.interactionCount = self.interactions.sum(axis=0).astype('int32')  # Nr in paper
    #print('self.interactionCount', self.interactionCount, numpy.sum(self.interactionCount), self.interactionCount.shape)
    
    self.bias = (self.interactionCount > 0).astype('float')  # 0 or 1 depending on whether any interactions
    
    self.cscore = MAXC * numpy.random.uniform(-1, 1, self.numWindows) * self.bias
    
    self.h = numpy.zeros(self.dBinMax - self.dBinMin + 1, dtype='float')
    
    self.sumbh = numpy.zeros(self.numWindows, dtype='float')
    self.sumbch = numpy.zeros(self.numWindows, dtype='float')
    
    self.interactionDCount = numpy.zeros(self.dBinMax - self.dBinMin + 1, dtype='int32')  # Nd in paper
    
    CscoreUtil.initDCount(self.interactionDCount, self.interactions, self.windowSize, self.dBinMin)
        
  def run(self):

    self.initData()
    
    print('Chromo {} has {} windows and {} interactions'.format(self.chromo, self.numWindows, numpy.sum(self.interactionCount)))
    
    logLikelihood = -1.0e308
   
    hTime = biasTime = cscoreTime = likelihoodTime = 0
    
    while True:
      print('min bias = {}'.format(numpy.min(self.bias[numpy.nonzero(self.bias)])))
      t0 = time.time()
      CscoreUtil.updateH(self.cscore, self.bias, self.h, self.sumbh, self.sumbch,
                         self.interactionDCount, self.windowSize, self.dBinMin)
      t1 = time.time()
      hTime += t1 - t0
      t0 = time.time()
      newLogLikelihood = CscoreUtil.calcLogLikelihood(self.cscore, self.bias, self.h,
                               self.interactionCount, self.interactionDCount, self.interactions)
      t1 = time.time()
      likelihoodTime += t1 - t0
      print('logLikelihood = {}'.format(newLogLikelihood))
      if newLogLikelihood < logLikelihood + 1:
        break
      logLikelihood = newLogLikelihood
      t0 = time.time()
      CscoreUtil.updateBias(self.cscore, self.bias, self.h, self.sumbh, self.sumbch,
                            self.interactionCount, self.windowSize, self.dBinMin)
      t1 = time.time()
      biasTime += t1 - t0
      t0 = time.time()
      CscoreUtil.updateCscore(self.cscore, self.bias, self.h, self.sumbch,
                              self.interactions, self.windowSize, self.dBinMin)
      t1 = time.time()
      cscoreTime += t1 - t0
      
    print('hTime = {}'.format(hTime))
    print('biasTime = {}'.format(biasTime))
    print('cscoreTime = {}'.format(cscoreTime))
    print('likelihoodTime = {}'.format(likelihoodTime))
        
def main(chromoSizeFile, inputSummaryFile, outputPrefix, windowSize, minDis):
  
  dBinMin = calcBinFromD(minDis)
  print('dBinMin = {}'.format(dBinMin))
  
  chromoSizeDict = readChromoSizes(chromoSizeFile)  # could skip this step and use largest position found...
    
  chromoAnalysisDict = {}
  for chromo in chromoSizeDict:
    chromoAnalysisDict[chromo] = ChromoAnalysis(chromo, chromoSizeDict[chromo], windowSize, dBinMin)
  
  if inputSummaryFile.endswith('.ncc'):
    readNccFile(inputSummaryFile, chromoAnalysisDict)
  else:
    readInteractions(inputSummaryFile, chromoAnalysisDict)
    
  for chromo in sorted(chromoSizeDict):
    chromoAnalysisDict[chromo].run()
    
    outputFile = '%s_%s_out.txt' % (outputPrefix, chromo)
    with open(outputFile, 'w') as fp:
      for score in chromoAnalysisDict[chromo].cscore:
        fp.write('{}\n'.format(score))
  
if __name__ == '__main__':
  
  import sys

  if len(sys.argv) != 6:
    print('Arguments: <chrSizes.txt> <input.summary> <OutputPrefix> <windowSize> <minDis>')
    sys.exit()
  
  chromoSizeFile, inputSummaryFile, outputPrefix, windowSize, minDis = sys.argv[1:]
  windowSize = int(windowSize)
  minDis = int(minDis)
  
  main(chromoSizeFile, inputSummaryFile, outputPrefix, windowSize, minDis)
  
  