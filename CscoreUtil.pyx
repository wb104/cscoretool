from numpy cimport ndarray
from numpy import log, log10, zeros

cdef int calcBinFromD(double d):
  
  cdef double steplength=0.04
  
  return int((log10(d+0.01)-3) / steplength)

cdef int calcBinFromIndex(int i, int j, int windowSize, int dBinMin):

  cdef int d = windowSize * abs(j - i)
  cdef int dBin = calcBinFromD(d) - dBinMin

  return dBin
  
def updateInteractions(ndarray[int, ndim=2] interactions,
                       int position1, int position2, int chromoSize, int windowSize, int dBinMin):
    
  if position1 > chromoSize or position2 > chromoSize:
    return
      
  cdef int ind1, ind2, dBin
    
  ind1 = position1 // windowSize
  ind2 = position2 // windowSize
  
  dBin = calcBinFromIndex(ind1, ind2, windowSize, dBinMin)
  if dBin >= 0:
    interactions[ind1][ind2] += 1
    if ind1 != ind2:
      interactions[ind2][ind1] += 1

def initDCount(ndarray[int, ndim=1] interactionDCount,
               ndarray[int, ndim=2] interactions,
               int windowSize, int dBinMin):

  cdef int numWindows = len(interactions)
  cdef int i, j, dBin

  for i in range(numWindows):
    for j in range(i, numWindows):
      dBin = calcBinFromIndex(i, j, windowSize, dBinMin)
      if dBin < 0:
        continue
      interactionDCount[dBin] += interactions[i][j]

def updateBias(ndarray[double, ndim=1] cscore,
               ndarray[double, ndim=1] bias,
               ndarray[double, ndim=1] h,
               ndarray[double, ndim=1] sumbh,
               ndarray[double, ndim=1] sumbch,
               ndarray[int, ndim=1] interactionCount,
	       int windowSize, int dBinMin):
       
  cdef int numWindows = len(cscore)
  cdef int i, j, dBin
  
  for i in range(numWindows):
    if interactionCount[i] > 0:
      x = interactionCount[i] / (sumbh[i] + cscore[i]*sumbch[i])
      bb = x - bias[i]
        
      for j in range(numWindows):
        dBin = calcBinFromIndex(i, j, windowSize, dBinMin)
        if dBin >= 0:
          u = bb * h[dBin]
          sumbh[j] += u
          sumbch[j] += u * cscore[i]
       
      bias[i] = x

def updateH(ndarray[double, ndim=1] cscore,
            ndarray[double, ndim=1] bias,
            ndarray[double, ndim=1] h,
            ndarray[double, ndim=1] sumbh,
            ndarray[double, ndim=1] sumbch,
            ndarray[int, ndim=1] interactionDCount,
            int windowSize, int dBinMin):
  
  cdef int numWindows = len(cscore)
  cdef int dBinMax = calcBinFromIndex(0, numWindows, windowSize, dBinMin) + dBinMin
  cdef int i, j, dBin
  cdef double u, F, G
  
  cdef ndarray[double, ndim=1] sumcc = zeros(dBinMax - dBinMin + 1)
  
  for i in range(numWindows):
    for j in range(i, numWindows):
      dBin = calcBinFromIndex(i, j, windowSize, dBinMin)
      if dBin >= 0:
        sumcc[dBin] += bias[i] * bias[j] * (1 + cscore[i]*cscore[j])

  for i in range(len(interactionDCount)):
    h[i] = interactionDCount[i] / max(sumcc[i], 1.0e-10)

  for i in range(numWindows):
    F = G = 0
    for j in range(numWindows):
      dBin = calcBinFromIndex(i, j, windowSize, dBinMin)
      if dBin >= 0:
        u = bias[j] * h[dBin]
        F += u
        G += u * cscore[j]
    sumbh[i] = F
    sumbch[i] = G  

def updateCscore(ndarray[double, ndim=1] cscore,
                 ndarray[double, ndim=1] bias,
                 ndarray[double, ndim=1] h,
                 ndarray[double, ndim=1] sumbch,
                 ndarray[int, ndim=2] interactions,
                 int windowSize, int dBinMin):

  cdef double MAXC = 0.9999  # not sure really need this
  cdef double maxdevc = 0
  cdef int i, j, nij
  cdef double x, y, z, fl, gl, u, v, cc

  for i, x in enumerate(cscore):
    F = bias[i] * sumbch[i]

    for n in range(100):  # iterations of Newton Raphson
      fl = gl = 0
      for j, nij in enumerate(interactions[i]):
        if nij > 0:
          u = cscore[j] / (1 + x*cscore[j])
          y = nij * u
          z = y * u
          fl += y
          gl += z
    
      v = fl - F
      if x == MAXC and v > 0:
        break
      if x == -MAXC and v < 0:
        break
      if abs(v) < 1.0e-8:
        break
      x += v / gl
      x = min(x, MAXC)
      x = max(x, -MAXC)

    cc = x - cscore[i]
    if abs(cc) > maxdevc:
      maxdevc = abs(bias[i] * cc)

    for j in range(len(sumbch)):
      dBin = calcBinFromIndex(i, j, windowSize, dBinMin)
      if dBin >= 0:
        sumbch[j] += cc * bias[i] * h[dBin]

    cscore[i] = x

def calcLogLikelihood(ndarray[double, ndim=1] cscore,
                 ndarray[double, ndim=1] bias,
                 ndarray[double, ndim=1] h,
                 ndarray[int, ndim=1] interactionCount,
                 ndarray[int, ndim=1] interactionDCount,
                 ndarray[int, ndim=2] interactions):
   
  cdef int numWindows = len(cscore)
  cdef int i, j, dBin
  cdef double logLikelihood = 0
  cdef double z
  
  logLikelihood = 0
  for i in range(numWindows):
    if bias[i] > 0:
      logLikelihood += interactionCount[i] * log(bias[i])
      for j in range(i, numWindows):
        z = 1 + cscore[i]*cscore[j]
        if z > 0:
          logLikelihood += interactions[i][j] * log(z)
         
  for dBin in range(len(interactionDCount)):
    if h[dBin] > 0:
      logLikelihood += interactionDCount[dBin] * log(h[dBin])
  
  return logLikelihood
