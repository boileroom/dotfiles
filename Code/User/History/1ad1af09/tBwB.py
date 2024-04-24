# /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import qetpy as qp
from qetpy import Noise
from qetpy.sim import TESnoise
from qetpy.plotting import compare_noise, plot_noise_sim
from pprint import pprint
import math
import os
import glob
import sys
import linecache
from array import array
import ROOT
from scipy.signal import butter, lfilter
from scipy import interpolate, optimize
import random

def butter_bandpass(lowcut, fs, order=5):
    low = lowcut / (0.5 * fs)
    b, a = butter(order, low, btype='lowpass')
    return b, a


def butter_bandpass_filter(data, lowcut, fs, order=5):
    b, a = butter_bandpass(lowcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

files = glob.glob("Run43_0718_Sapphire902-ABCD_slot5_DG50_QET50pct_PTon_Pulse/*.txt")

tF = ROOT.TFile("Run43_0721_slot5_Sapphire902-ABCD_DG50_QET50pct_NoProbeTune_PTon_Pulse.root")

nF = ROOT.TFile("Run43_0721_slot5_Sapphire902-ABCD_DG50_QET50pct_NoProbeTune_PTon_Noise.root")

saveF = ROOT.TFile("Run43_0718_Sapphire902-ABCD_slot5_DG50_QET50pct_PTon_Pulse_fitPhonon.root","RECREATE")

#volt = 5
#voltD = 10
templateAa = tF.Get("pulse3A_template")
templateBa = tF.Get("pulse3B_template")
templateCa = tF.Get("pulse3C_template")
templateDa = tF.Get("pulse3D_template")

templateAb = templateAa.GetY()
templateBb = templateBa.GetY()
templateCb = templateCa.GetY()
templateDb = templateDa.GetY()

templateA = np.array(np.frombuffer(templateAb, dtype=np.double))
templateB = np.array(np.frombuffer(templateBb, dtype=np.double))
templateC = np.array(np.frombuffer(templateCb, dtype=np.double))
templateD = np.array(np.frombuffer(templateDb, dtype=np.double))

chANoiseA = nF.Get("noiseA_template")
chBNoiseA = nF.Get("noiseB_template")
chCNoiseA = nF.Get("noiseC_template")
chDNoiseA = nF.Get("noiseD_template")

chANoiseB = chANoiseA.GetY()
chBNoiseB = chBNoiseA.GetY()
chCNoiseB = chCNoiseA.GetY()
chDNoiseB = chDNoiseA.GetY()

chANoise = np.array(np.frombuffer(chANoiseB, dtype=np.double))
chBNoise = np.array(np.frombuffer(chBNoiseB, dtype=np.double))
chCNoise = np.array(np.frombuffer(chCNoiseB, dtype=np.double))
chDNoise = np.array(np.frombuffer(chDNoiseB, dtype=np.double))


gain = 50.
gainD = 50.
nSamples = 1900

time = np.arange(nSamples,dtype=np.float64)
samplelength = 1000.
sl_sec = samplelength*1e-9
time2 = time*samplelength

ADCconv = (1)/(gain*1200.*10.*2)  # conversion from Vout to Iin, the first 1000 is a factor from the CAEN data writing
#ADCconvN = 1000*((2.25/16384.))/(gainN*1200.*10.)  # conversion from Vout to Iin, the first 1000 is a factor from the CAEN data writing
#ADCconvc = 1000*(2.25/16384.)
ADCconvD = (1)/(gainD*1200.*10.*2)


nE = 0
nN = 0

time = 0
tt = 0
trigE = 0
trig = 0
tbin = 1

sf = (1/(samplelength*1e-9))

#templateDF =  butter_bandpass_filter(templateD, 50000, sf, order=2)
#chDNoiseF = butter_bandpass_filter(chDNoise, 50000, sf, order=2)

t = ROOT.TTree('data', 'tree with events')

ampOFA = np.zeros(1, dtype=np.float64)
ampOFB = np.zeros(1, dtype=np.float64)
ampOFC = np.zeros(1, dtype=np.float64)
ampOFD = np.zeros(1, dtype=np.float64)

chi2OFA = np.zeros(1, dtype=np.float64)
chi2OFB = np.zeros(1, dtype=np.float64)
chi2OFC = np.zeros(1, dtype=np.float64)
chi2OFD = np.zeros(1, dtype=np.float64)

chi2OFAnop = np.zeros(1, dtype=np.float64)
chi2OFBnop = np.zeros(1, dtype=np.float64)
chi2OFCnop = np.zeros(1, dtype=np.float64)
chi2OFDnop = np.zeros(1, dtype=np.float64)

t0OFA = np.zeros(1, dtype=np.float64)
t0OFB = np.zeros(1, dtype=np.float64)
t0OFC = np.zeros(1, dtype=np.float64)
t0OFD = np.zeros(1, dtype=np.float64)

rtft20A = np.zeros(1, dtype=np.float64)
rtft20B = np.zeros(1, dtype=np.float64)
rtft20C = np.zeros(1, dtype=np.float64)
rtft20D = np.zeros(1, dtype=np.float64)

partXOF = np.zeros(1, dtype=np.float64)
partYOF = np.zeros(1, dtype=np.float64)
delayXOF = np.zeros(1, dtype=np.float64)
delayYOF = np.zeros(1, dtype=np.float64)
delayX = np.zeros(1, dtype=np.float64)
delayY = np.zeros(1, dtype=np.float64)
eventT = np.zeros(1, dtype=np.float64)
meanA = np.zeros(1, dtype=np.float64)
meanB = np.zeros(1, dtype=np.float64)
meanC = np.zeros(1, dtype=np.float64)
meanD = np.zeros(1, dtype=np.float64)

IntTA = np.zeros(1, dtype=np.float64)
IntTB = np.zeros(1, dtype=np.float64)
IntTC = np.zeros(1, dtype=np.float64)
IntTD = np.zeros(1, dtype=np.float64)

stdApre = np.zeros(1, dtype=np.float64)
stdBpre = np.zeros(1, dtype=np.float64)
stdCpre = np.zeros(1, dtype=np.float64)
stdDpre = np.zeros(1, dtype=np.float64)

stdAtail = np.zeros(1, dtype=np.float64)
stdBtail = np.zeros(1, dtype=np.float64)
stdCtail = np.zeros(1, dtype=np.float64)
stdDtail = np.zeros(1, dtype=np.float64)

minA = np.zeros(1, dtype=np.float64)
minB = np.zeros(1, dtype=np.float64)
minC = np.zeros(1, dtype=np.float64)
minD = np.zeros(1, dtype=np.float64)

maxbinA = np.zeros(1, dtype=np.float64)
maxbinB = np.zeros(1, dtype=np.float64)
maxbinC = np.zeros(1, dtype=np.float64)
maxbinD = np.zeros(1, dtype=np.float64)

maxtailA = np.zeros(1, dtype=np.float64)
maxtailB = np.zeros(1, dtype=np.float64)
maxtailC = np.zeros(1, dtype=np.float64)
maxtailD = np.zeros(1, dtype=np.float64)

riseTA = np.zeros(1, dtype=np.float64)
riseTB = np.zeros(1, dtype=np.float64)
riseTC = np.zeros(1, dtype=np.float64)
riseTD = np.zeros(1, dtype=np.float64)

delayTA = np.zeros(1, dtype=np.float64)
delayTB = np.zeros(1, dtype=np.float64)
delayTC = np.zeros(1, dtype=np.float64)
delayTD = np.zeros(1, dtype=np.float64)

riseOFA = np.zeros(1, dtype=np.float64)
riseOFB = np.zeros(1, dtype=np.float64)
riseOFC = np.zeros(1, dtype=np.float64)
riseOFD = np.zeros(1, dtype=np.float64)

riseOFA90 = np.zeros(1, dtype=np.float64)
riseOFB90 = np.zeros(1, dtype=np.float64)
riseOFC90 = np.zeros(1, dtype=np.float64)
riseOFD90 = np.zeros(1, dtype=np.float64)
pminrt90 = np.zeros(1, dtype=np.float64)
pmaxrt90 = np.zeros(1, dtype=np.float64)
pamp = np.zeros(1, dtype=np.float64)
trigger = np.zeros(1, dtype=np.int64)
time = np.zeros(1, dtype=np.int64)

t.Branch( 'ampOFA', ampOFA, 'ampOFA/D')
t.Branch( 'ampOFB', ampOFB, 'ampOFB/D')
t.Branch( 'ampOFC', ampOFC, 'ampOFC/D')
t.Branch( 'ampOFD', ampOFD, 'ampOFD/D')

t.Branch( 'chi2OFA', chi2OFA, 'chi2OFA/D')
t.Branch( 'chi2OFB', chi2OFB, 'chi2OFB/D')
t.Branch( 'chi2OFC', chi2OFC, 'chi2OFC/D')
t.Branch( 'chi2OFD', chi2OFD, 'chi2OFD/D')

t.Branch( 'chi2OFAnop', chi2OFAnop, 'chi2OFAnop/D')
t.Branch( 'chi2OFBnop', chi2OFBnop, 'chi2OFBnop/D')
t.Branch( 'chi2OFCnop', chi2OFCnop, 'chi2OFCnop/D')
t.Branch( 'chi2OFDnop', chi2OFDnop, 'chi2OFDnop/D')

t.Branch( 't0OFA', t0OFA, 't0OFA/D')
t.Branch( 't0OFB', t0OFB, 't0OFB/D')
t.Branch( 't0OFC', t0OFC, 't0OFC/D')
t.Branch( 't0OFD', t0OFD, 't0OFD/D')

t.Branch( 'rtft20A', rtft20A, 'rtft20A/D')
t.Branch( 'rtft20B', rtft20B, 'rtft20A/D')
t.Branch( 'rtft20C', rtft20C, 'rtft20A/D')
t.Branch( 'rtft20D', rtft20D, 'rtft20A/D')

t.Branch( 'partXOF', partXOF, 'partXOF/D')
t.Branch( 'partYOF', partYOF, 'partYOF/D')
t.Branch( 'delayXOF', delayXOF, 'delayXOF/D')
t.Branch( 'delayYOF', delayYOF, 'delayYOF/D')
t.Branch( 'delayX', delayX, 'delayX/D')
t.Branch( 'delayY', delayY, 'delayY/D')
t.Branch( 'eventT', eventT, 'eventT/D')
t.Branch( 'meanA',meanA,'meanA/D') 
t.Branch( 'meanB',meanB,'meanB/D')
t.Branch( 'meanC',meanC,'meanC/D')
t.Branch( 'meanD',meanD,'meanD/D')

t.Branch( 'IntTA',IntTA,'IntTA/D')
t.Branch( 'IntTB',IntTB,'IntTB/D')
t.Branch( 'IntTC',IntTC,'IntTC/D')
t.Branch( 'IntTD',IntTD,'IntTD/D')

t.Branch( 'stdApre',stdApre,'stdApre/D')
t.Branch( 'stdBpre',stdBpre,'stdBpre/D')
t.Branch( 'stdCpre',stdCpre,'stdCpre/D')
t.Branch( 'stdDpre',stdDpre,'stdDpre/D')

t.Branch( 'stdAtail',stdAtail,'stdAtail/D')
t.Branch( 'stdBtail',stdBtail,'stdBtail/D')
t.Branch( 'stdCtail',stdCtail,'stdCtail/D')
t.Branch( 'stdDtail',stdDtail,'stdDtail/D')

t.Branch( 'minA',minA,'minA/D')
t.Branch( 'minB',minB,'minB/D')
t.Branch( 'minC',minC,'minC/D')
t.Branch( 'minD',minD,'minD/D')

t.Branch( 'maxbinA',maxbinA,'maxbinA/D')
t.Branch( 'maxbinB',maxbinB,'maxbinB/D')
t.Branch( 'maxbinC',maxbinC,'maxbinC/D')
t.Branch( 'maxbinD',maxbinD,'maxbinD/D')

t.Branch( 'maxtailA',maxtailA,'maxtailA/D')
t.Branch( 'maxtailB',maxtailB,'maxtailB/D')
t.Branch( 'maxtailC',maxtailC,'maxtailC/D')
t.Branch( 'maxtailD',maxtailD,'maxtailD/D')

t.Branch( 'riseTA', riseTA, 'riseTA/D')
t.Branch( 'riseTB', riseTB, 'riseTB/D')
t.Branch( 'riseTC', riseTC, 'riseTC/D')
t.Branch( 'riseTD', riseTD, 'riseTD/D')

t.Branch( 'delayTA', delayTA, 'delayTA/D')
t.Branch( 'delayTB', delayTB, 'delayTB/D')
t.Branch( 'delayTC', delayTC, 'delayTC/D')
t.Branch( 'delayTD', delayTD, 'delayTD/D')

t.Branch( 'riseOFA', riseOFA, 'riseOFA/D')
t.Branch( 'riseOFB', riseOFB, 'riseOFB/D')
t.Branch( 'riseOFC', riseOFC, 'riseOFC/D')
t.Branch( 'riseOFD', riseOFD, 'riseOFD/D')

t.Branch( 'riseOFA90', riseOFA90, 'riseOFA90/D')
t.Branch( 'riseOFB90', riseOFB90, 'riseOFB90/D')
t.Branch( 'riseOFC90', riseOFC90, 'riseOFC90/D')
t.Branch( 'riseOFD90', riseOFD90, 'riseOFD90/D')
t.Branch( 'trigger',trigger,'trigger/I')
t.Branch( 'pminrt90', pminrt90, 'pminrt90/D')
t.Branch( 'pmaxrt90', pmaxrt90, 'pmaxrt90/D')
t.Branch( 'time',time,'time/I')

pA = np.zeros(nSamples)
pB = np.zeros(nSamples)
pC = np.zeros(nSamples)
pD = np.zeros(nSamples)

pT1A = np.zeros(nSamples)
pT1B = np.zeros(nSamples)
pT1C = np.zeros(nSamples)
pT1D = np.zeros(nSamples)

pT2A = np.zeros(nSamples)
pT2B = np.zeros(nSamples)
pT2C = np.zeros(nSamples)
pT2D = np.zeros(nSamples)

pT3A = np.zeros(nSamples)
pT3B = np.zeros(nSamples)
pT3C = np.zeros(nSamples)
pT3D = np.zeros(nSamples)

pT4A = np.zeros(nSamples)
pT4B = np.zeros(nSamples)
pT4C = np.zeros(nSamples)
pT4D = np.zeros(nSamples)

pT5A = np.zeros(nSamples)
pT5B = np.zeros(nSamples)
pT5C = np.zeros(nSamples)
pT5D = np.zeros(nSamples)


nE = 0
nEf = 0
upperbinA = 0
lowerbinA = 0
upperbinB = 0
lowerbinB = 0
upperbinC = 0
lowerbinC = 0
upperbinD = 0
lowerbinD = 0


maxbinA1 = 0
maxbinB1 = 0
maxbinC1 = 0
maxbinD1 = 0


for filen in files:
  nEf = 0
  f = open(filen, 'r')
  linenmbr = 0
  for line in f:
      linenmbr = linenmbr + 1
      if (linenmbr > 3 and linenmbr < nSamples + 4):
          stuff = line.split()
          try:
              pA[linenmbr - 4] = float(stuff[1])
          except:
              pA[linenmbr - 4] = 10e6 * ADCconv
          try:
              pB[linenmbr - 4] = float(stuff[2])
          except:
              pB[linenmbr - 4] = 10e6 * ADCconv
          try:
              pC[linenmbr - 4] = float(stuff[3])
          except:
              pC[linenmbr - 4] = 10e6 * ADCconv
          try:
              pD[linenmbr - 4] = float(stuff[4])
          except:
              pD[linenmbr - 4] = 10e6 * ADCconv
  nE += 1
  nEf += 1
  triggerTime = .2*time2[nSamples-1]
  eventTime=nE


  if nE % 100 == 0: print(nE)
  meanA[0] = pA[0:400].mean()
  meanB[0] = pB[0:400].mean()
  meanC[0] = pC[0:400].mean()
  meanD[0] = pD[0:400].mean()

  # meanE[0] = pE[0:400].mean()

  pA -= meanA[0]
  pB -= meanB[0]
  pC -= meanC[0]
  pD -= meanD[0]

  # pE -= meanE[0]

  pA *= ADCconv
  pB *= ADCconv
  pC *= ADCconv
  pD *= ADCconv

  # pE *= ADCconv

  stdApre[0] = pA[0:400].std()
  stdBpre[0] = pB[0:400].std()
  stdCpre[0] = pC[0:400].std()
  stdDpre[0] = pD[0:400].std()

  stdAtail[0] = pA[nSamples - 200:nSamples].std()
  stdBtail[0] = pB[nSamples - 200:nSamples].std()
  stdCtail[0] = pC[nSamples - 200:nSamples].std()
  stdDtail[0] = pD[nSamples - 200:nSamples].std()

  minA[0] = pA.min()
  minB[0] = pB.min()
  minC[0] = pC.min()
  minD[0] = pD.min()

  maxtailA[0] = pA[nSamples - 200:nSamples].max()
  maxtailB[0] = pB[nSamples - 200:nSamples].max()
  maxtailC[0] = pC[nSamples - 200:nSamples].max()
  maxtailD[0] = pD[nSamples - 200:nSamples].max()

  if trigger[0] != '10000000000':
      OFA = qp.OptimumFilter(pA, templateA, chANoise, sf)  # initialize the OptimumFilter class
      ampA, t0A, chi2A = OFA.ofamp_withdelay()
      chi2OFAnop[0] = OFA.chi2_nopulse()
      ampOFA[0] = ampA
      Integ = ampA * templateA
      IntTA[0] = Integ.sum()
      chi2OFA[0] = chi2A
      t0OFA[0] = t0A
      # pfilt = butter_bandpass_filter(pA, 50000, sf, order=2)
      pfilt = pA
      maxbin = np.argmax(pfilt)
      maxbinA1 = np.argmax(pfilt)
      # maxbin += 500
      # maxbinA1 += 500
      pmax = pfilt.max()
      lowerbinA = 0
      for i in range(nSamples):
          if pfilt[maxbin - i] < .2 * pmax:
              # m = (pfilt[maxbin-i+1])/sl_sec
              # rtft20A[0] = (.2*pmax)/m + (maxbin-i)*sl_sec
              lowerbinA = maxbin - i
              break

      for i in range(nSamples):
          if maxbin + i > nSamples - 1:
              upperbinA = nSamples - 1
              break
          elif pfilt[maxbin + i] < .2 * pmax:
              upperbinA = maxbin + i
              break

      rtft20A[0] = pfilt[lowerbinA:upperbinA].sum()
      maxbinTA = np.argmax(pfilt[200:(nSamples - 200)])
      maxbinTA += 200
      pmaxTA = pfilt.max()
      riseTA[0] = 0
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .2 * pmaxTA:
              m = (pfilt[maxbinTA - i + 1]) / sl_sec
              riseTA[0] = (.2 * pmaxTA) / m + (maxbinTA - i) * sl_sec
              break
      # Calculate 10% and 40% rise timestamps to find risetime of the pulse
      rise10 = 0
      rise20 = 0
      rise40 = 0
      riseOFA[0] = 0
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .4 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA)).roots()
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .4 * pmaxTA
              # plt.plot(time2,pfilt,'r.',x,f(x),'g--')
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise40 = root
              except (ValueError, RuntimeError):
                  rise40 = 0
              break
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .2 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA)).roots()
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .2 * pmaxTA
              # plt.plot(time2,pfilt,'r.')
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise20 = root
                  delayTA[0] = (root - triggerTime) * 1e-3
              except (ValueError, RuntimeError):
                  rise20 = 0
                  delayTA[0] = 0
              break

      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .1 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .1 * pmaxTA
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise10 = root
              except (ValueError, RuntimeError):
                  rise10 = 0
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.1*pmaxTA),k=1).roots()
              # if (len(spl)==0):
              #    rise10= 0
              # else:
              #    rise10 = spl[0]
              break
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .9 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA)).roots()
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .9 * pmaxTA
              # plt.plot(time2,pfilt,'r.')
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise90 = root
                  delayTA[0] = (root - triggerTime) * 1e-3
              except (ValueError, RuntimeError):
                  rise90 = 0
                  delayTA[0] = 0
              break
      if rise10 != 0 and rise40 != 0:
          riseOFA[0] = (rise40 - rise10) * 1e-3
      if rise10 != 0 and rise90 != 0:
          riseOFA90[0] = (rise90 - rise10) * 1e-3

      OFB = qp.OptimumFilter(pB, templateB, chBNoise, sf)  # initialize the OptimumFilter class
      # ampB, t0B, chi2B = OFB.ofamp_withdelay(nconstrain=600, windowcenter=-280)
      ampB, t0B, chi2B = OFB.ofamp_withdelay()
      # try:
      #    NOFB = qp.OFnonlin(chBNoise,sf,template=None)
      #    risetime = 0
      #    thrA, thrB, thrC, risetime, fall1, fall2, fall3, thrt0 = NOFB.fit_falltimes(pB, npolefit=4)
      # except ValueError:
      #    print("ValueError")
      # riseOFB[0] = risetime
      Integ = ampB * templateB
      IntTB[0] = Integ.sum()
      chi2OFBnop[0] = OFB.chi2_nopulse()
      ampOFB[0] = ampB
      chi2OFB[0] = chi2B
      t0OFB[0] = t0B
      # pfilt = butter_bandpass_filter(pB, 50000, sf, order=2)
      pfilt = pB
      maxbin = np.argmax(pfilt)
      maxbinB1 = np.argmax(pfilt)
      # maxbin += 500
      # maxbinB1 += 500
      pmax = pfilt.max()
      lowerbinB = 0
      for i in range(nSamples):
          if pfilt[maxbin - i] < .2 * pmax:
              lowerbinB = maxbin - i
              break
      for i in range(nSamples):
          if maxbin + i > nSamples - 1:
              upperbinB = nSamples - 1
              break
          elif pfilt[maxbin + i] < .2 * pmax:
              upperbinB = maxbin + i
              break
      rtft20B[0] = pfilt[lowerbinB:upperbinB].sum()
      maxbinTA = np.argmax(pfilt[200:(nSamples - 200)])
      maxbinTA += 200
      pmaxTA = pfilt.max()
      riseTB[0] = 0
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .2 * pmaxTA:
              m = (pfilt[maxbinTA - i + 1]) / sl_sec
              riseTB[0] = (.2 * pmaxTA) / m + (maxbinTA - i) * sl_sec
              break
      # Calculate 10% and 40% rise timestamps to find risetime of the pulse
      rise10 = 0
      rise20 = 0
      rise40 = 0
      rise90 = 0
      riseOFB[0] = 0
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .4 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .4 * pmaxTA
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise40 = root
              except (ValueError, RuntimeError):
                  rise40 = 0
              # print("riseB")
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA),k=1).roots()
              # if (len(spl)==0):
              #    rise40= 0
              # else:
              #    rise40 = spl[0]
              break

      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .2 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA)).roots()
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .2 * pmaxTA
              # plt.plot(time2,pfilt,'r.')
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise20 = root
                  delayTB[0] = (root - triggerTime) * 1e-3
              except (ValueError, RuntimeError):
                  rise20 = 0
                  delayTB[0] = 0
              break
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .1 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .1 * pmaxTA
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise10 = root
              except (ValueError, RuntimeError):
                  rise10 = 0
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.1*pmaxTA),k=1).roots()
              # if (len(spl)==0):
              #    rise10= 0
              # else:
              #    rise10 = spl[0]
              break
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .9 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA)).roots()
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .9 * pmaxTA
              # plt.plot(time2,pfilt,'r.')
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise90 = root
                  delayTB[0] = (root - triggerTime) * 1e-3
              except (ValueError, RuntimeError):
                  rise90 = 0
                  delayTB[0] = 0
              break
      if rise10 != 0 and rise40 != 0:
          riseOFB[0] = (rise40 - rise10) * 1e-3
      if rise10 != 0 and rise90 != 0:
          riseOFB90[0] = (rise90 - rise10) * 1e-3

      OFC = qp.OptimumFilter(pC, templateC, chCNoise, sf)  # initialize the OptimumFilter class
      ampC, t0C, chi2C = OFC.ofamp_withdelay()
      # try:
      #    NOFC = qp.OFnonlin(chCNoise,sf,template=None)
      #    risetime = 0
      #    thrA, thrB, thrC, risetime, fall1, fall2, fall3, thrt0 = NOFC.fit_falltimes(pC, npolefit=4)
      # except ValueError:
      #    print("ValueError")
      # riseOFC[0] = risetime
      Integ = ampC * templateC
      IntTC[0] = Integ.sum()
      chi2OFCnop[0] = OFC.chi2_nopulse()
      ampOFC[0] = ampC
      chi2OFC[0] = chi2C
      t0OFC[0] = t0C
      # pfilt = butter_bandpass_filter(pC, 50000, sf, order=2)
      pfilt = pC
      maxbin = np.argmax(pfilt)
      maxbinC1 = np.argmax(pfilt)
      # maxbin += 500
      # maxbinC1 += 500
      pmax = pfilt.max()
      lowerbinC = 0
      for i in range(nSamples):
          if pfilt[maxbin - i] < .2 * pmax:
              lowerbinC = maxbin - i
              break
      for i in range(nSamples):
          if maxbin + i > nSamples - 1:
              upperbinC = nSamples - 1
              break
          elif pfilt[maxbin + i] < .2 * pmax:
              upperbinC = maxbin + i
              break
      rtft20C[0] = pfilt[lowerbinC:upperbinC].sum()
      maxbinTA = np.argmax(pfilt[200:(nSamples - 200)])
      maxbinTA += 200
      pmaxTA = pfilt.max()
      riseTC[0] = 0
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .2 * pmaxTA:
              m = (pfilt[maxbinTA - i + 1]) / sl_sec
              riseTC[0] = (.2 * pmaxTA) / m + (maxbinTA - i) * sl_sec
              break
      # Calculate 10% and 40% rise timestamps to find risetime of the pulse
      rise10 = 0
      rise20 = 0
      rise40 = 0
      rise90 = 0
      riseOFC[0] = 0
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .4 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .4 * pmaxTA
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise40 = root
              except (ValueError, RuntimeError):
                  rise40 = 0
              # print("riseC")
              # print(y)
              # print(x)
              # print(.4*pmaxTA)
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA),k=1).roots()
              # if (len(spl)==0):
              #    rise40= 0
              # else:
              #    rise40 = spl[0]
              break
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .2 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA)).roots()
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .2 * pmaxTA
              # plt.plot(time2,pfilt,'r.')
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise20 = root
                  delayTC[0] = (root - triggerTime) * 1e-3
              except (ValueError, RuntimeError):
                  rise20 = 0
                  delayTC[0] = 0
              break
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .1 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .1 * pmaxTA
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise10 = root
              except (ValueError, RuntimeError):
                  rise10 = 0
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.1*pmaxTA),k=1).roots()
              # if (len(spl)==0):
              #    rise10= 0
              # else:
              #    rise10 = spl[0]
              break
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .9 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA)).roots()
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .9 * pmaxTA
              # plt.plot(time2,pfilt,'r.')
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise90 = root
                  # delayTC[0] = (root - triggerTime)*1e-3
              except (ValueError, RuntimeError):
                  rise90 = 0
                  delayTC[0] = 0
              break
      if rise10 != 0 and rise40 != 0:
          riseOFC[0] = (rise40 - rise10) * 1e-3
      if rise10 != 0 and rise90 != 0:
          riseOFC90[0] = (rise90 - rise10) * 1e-3

      OFD = qp.OptimumFilter(pD, templateD, chDNoise, sf)
      # initialize the OptimumFilter class
      ampD, t0D, chi2D = OFD.ofamp_withdelay()
      # try:
      #    NOFD = qp.OFnonlin(chDNoise,sf,template=None)
      #    risetime = 0
      #    thrA, thrB, thrC, risetime, fall1, fall2, fall3, thrt0 = NOFD.fit_falltimes(pD, npolefit=4)
      # except ValueError:
      #    print("ValueError")
      # riseOFD[0] = risetime
      Integ = ampD * templateD
      IntTD[0] = Integ.sum()
      # ampD, chi2D = OFD.ofamp_nodelay()
      chi2OFDnop[0] = OFD.chi2_nopulse()
      ampOFD[0] = ampD
      chi2OFD[0] = chi2D
      t0OFD[0] = t0D
      # pfilt = butter_bandpass_filter(pD, 50000, sf, order=2)
      pfilt = pD
      maxbin = np.argmax(pfilt)
      maxbinD1 = np.argmax(pfilt)
      # maxbin += 500
      # maxbinD1 += 500
      pmax = pfilt.max()
      lowerbinD = 0
      for i in range(nSamples):
          if pfilt[maxbin - i] < .2 * pmax:
              lowerbinD = maxbin - i
              break
      for i in range(nSamples):
          if maxbin + i > nSamples - 1:
              upperbinD = nSamples - 1
              break
          elif pfilt[maxbin + i] < .2 * pmax:
              upperbinD = maxbin + i
              break
      rtft20D[0] = pfilt[lowerbinD:upperbinD].sum()
      maxbinTA = np.argmax(pfilt[200:(nSamples - 200)])
      maxbinTA += 200
      pmaxTA = pfilt.max()
      riseTD[0] = 0
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .2 * pmaxTA:
              m = (pfilt[maxbinTA - i + 1]) / sl_sec
              riseTD[0] = (.2 * pmaxTA) / m + (maxbinTA - i) * sl_sec
              break
      time10 = 0
      time90 = 0
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .1 * pmaxTA:
              m = (pfilt[maxbinTA - i + 1]) / sl_sec
              time10 = (.1 * pmaxTA) / m + (maxbinTA - i) * sl_sec
              break

      # Calculate 10% and 40% rise timestamps to find risetime of the pulse
      rise10 = 0
      rise20 = 0
      rise40 = 0
      rise90 = 0
      riseOFD[0] = 0
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .4 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .4 * pmaxTA
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise40 = root
              except (ValueError, RuntimeError):
                  rise40 = 0
              # print("riseD")
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA),k=1).roots()
              # if (len(spl)==0):
              #    rise40= 0
              # else:
              #    rise40 = spl[0]
              break
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .2 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA)).roots()
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .2 * pmaxTA
              # plt.plot(time2,pfilt,'r.')
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise20 = root
                  delayTD[0] = (root - triggerTime) * 1e-3
              except (ValueError, RuntimeError):
                  rise20 = 0
                  delayTD[0] = 0
              break
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .1 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .1 * pmaxTA
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise10 = root
              except (ValueError, RuntimeError):
                  rise10 = 0
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.1*pmaxTA),k=1).roots()
              # if (len(spl)==0):
              #    rise10= 0
              # else:
              #   rise10 = spl[0]
              break
      for i in range(int(nSamples / 4)):
          if pfilt[maxbinTA - i] < .9 * pmaxTA and maxbinTA - i > 0:
              y = np.array(
                  [pfilt[maxbinTA - i - 1], pfilt[maxbinTA - i], pfilt[maxbinTA - i + 1], pfilt[maxbinTA - i + 2]])
              x = np.array(
                  [time2[maxbinTA - i - 1], time2[maxbinTA - i], time2[maxbinTA - i + 1], time2[maxbinTA - i + 2]])
              f1 = interpolate.interp1d(x, y, kind='linear')
              g = lambda x: f1(x) - .9 * pmaxTA
              try:
                  root = optimize.newton(g, time2[maxbinTA - i])
                  rise90 = root
              except (ValueError, RuntimeError):
                  rise90 = 0
              # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.1*pmaxTA),k=1).roots()
              # if (len(spl)==0):
              #    rise10= 0
              # else:
              #   rise10 = spl[0]
              break

      if rise10 != 0 and rise40 != 0:
          riseOFD[0] = (rise40 - rise10) * 1e-3
      if rise10 != 0 and rise90 != 0:
          riseOFD90[0] = (rise90 - rise10) * 1e-3

      # -------------------------------------------------------------------------------------------------Channel E ---------------------------------------------------------------

      # ------------------------------------------------------------------------------------------------------------------------------------------------------------------
      pminrt90[0] = min(riseOFA90[0], riseOFB90[0], riseOFC90[0], riseOFD90[0])
      pamp[0] = max(ampOFA[0], ampOFB[0], ampOFC[0], ampOFD[0])
      if pamp[0] == ampOFA[0]:
          pmaxrt90[0] = riseOFA90[0]
      if pamp[0] == ampOFB[0]:
          pmaxrt90[0] = riseOFB90[0]
      if pamp[0] == ampOFC[0]:
          pmaxrt90[0] = riseOFC90[0]
      if pamp[0] == ampOFD[0]:
          pmaxrt90[0] = riseOFD90[0]

      delayXOF[0] = -1 * (.866 * riseOFD[0] - .866 * riseOFB[0]) * 1e-3
      delayYOF[0] = -1 * (.5 * riseOFD[0] + .5 * riseOFB[0] - riseOFC[0]) * 1e-3

      delayX[0] = -1 * (.866 * delayTD[0] - .866 * delayTB[0]) * 1e-3
      delayY[0] = -1 * (.5 * delayTD[0] + .5 * delayTB[0] - delayTC[0]) * 1e-3
      partXOF[0] = (.866 * ampOFD[0] - .866 * ampOFB[0]) / (ampOFD[0] + ampOFB[0] + ampOFC[0])
      partYOF[0] = (.5 * ampOFD[0] + .5 * ampOFB[0] - ampOFC[0]) / (ampOFD[0] + ampOFB[0] + ampOFC[0])
      # delayXOF[0] = -1*(.866*rtft20D[0] - .866*rtft20A[0])*1e6
      # delayYOF[0] = -1*(.5*rtft20D[0] + .5*rtft20A[0] - rtft20C[0])*1e6
      eventT[0] = eventTime
      maxbinA[0] = maxbinA1
      maxbinB[0] = maxbinB1
      maxbinC[0] = maxbinC1
      maxbinD[0] = maxbinD1

      # if (ampOFA[0]+ampOFC[0]+ampOFD[0]) > .06 and (ampOFA[0]+ampOFC[0]+ampOFD[0]) < .1:
     
      
      if ampOFA[0]+ampOFB[0]+ampOFC[0]+ampOFD[0]<4e-6 and ampOFA[0]+ampOFB[0]+ampOFC[0]+ampOFD[0]>2.5e-6:
          traceA = ROOT.TGraph(nSamples, time2, pA)
          traceA.SetName("traceA" + str(nE))
          traceA.Write()
      
      '''
      if ampOFA[0] <1.1e-6 and ampOFA[0] > 0.9e-6:
          traceA = ROOT.TGraph(nSamples, time2-1000000, pA)
          traceA.SetName("traceA_8KeV" + str(nE))
          traceA.Write()
      '''
      t.Fill()
      #if nE >5000: break

      #f.close()


saveF.Write()
print(nE)
