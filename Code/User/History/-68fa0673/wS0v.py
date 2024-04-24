#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import qetpy as qp
from qetpy import Noise
from qetpy.sim import TESnoise
from qetpy.plotting import compare_noise, plot_noise_sim
from pprint import pprint
import linecache
import os
import glob
import sys
from array import array
import ROOT

#files = glob.glob("")
files = glob.glob("Run43_0721_slot5_Sapphire902-ABCD_DG50_QET50pct_NoProbeTune_PTon_Pulse/*.txt")

saveF = ROOT.TFile("Run43_0721_slot5_Sapphire902-ABCD_DG50_QET50pct_NoProbeTune_PTon_Pulse.root","RECREATE")

gain = 50.
#gainD = 50.
nSamples = 1900

nbin_pre = 400
nbin_post = 200

#nbin_pre = 3000
#nbin_post = 3000

std_A = ROOT.TH1D("std_A_pre"," ; ",1000,0,5e-6)
std_B = ROOT.TH1D("std_B_pre"," ; ",1000,0,5e-6)
std_C = ROOT.TH1D("std_C_pre"," ; ",1000,0,5e-6)
std_D = ROOT.TH1D("std_D_pre"," ; ",1000,0,5e-6)

max_A = ROOT.TH1D("max_A"," ; ",1000,0,40e-6)
max_B = ROOT.TH1D("max_B"," ; ",1000,0,40e-6)
max_C = ROOT.TH1D("max_C"," ; ",1000,0,40e-6)
max_D = ROOT.TH1D("max_D"," ; ",1000,0,40e-6)

max_A_tail = ROOT.TH1D("max_A_tail"," ; ",1000,0,40e-6)
max_B_tail = ROOT.TH1D("max_B_tail"," ; ",1000,0,40e-6)
max_C_tail = ROOT.TH1D("max_C_tail"," ; ",1000,0,40e-6)
max_D_tail = ROOT.TH1D("max_D_tail"," ; ",1000,0,40e-6)

min_A = ROOT.TH1D("min_A"," ; ",1000,-10e-6,0)
min_B = ROOT.TH1D("min_B"," ; ",1000,-10e-6,0)
min_C = ROOT.TH1D("min_C"," ; ",1000,-10e-6,0)
min_D = ROOT.TH1D("min_D"," ; ",1000,-10e-6,0)

stdp_A = ROOT.TH1D("std_A_post"," ; ",1000,0,5e-6)
stdp_B = ROOT.TH1D("std_B_post"," ; ",1000,0,5e-6)
stdp_C = ROOT.TH1D("std_C_post"," ; ",1000,0,5e-6)
stdp_D = ROOT.TH1D("std_D_post"," ; ",1000,0,5e-6)

time = np.arange(nSamples,dtype=np.float64)
samplelength = 1000 #in nanoseconds
time2 = time*samplelength

#ADCconvV = 1000*((2.25/16384.))/(gain*10.*2)  # conversion from Vout to Iin, the first 1000 is a factor from the CAEN data writing
#ADCconvN = 1000*((2.25/16384.))/(gainN*1200.*10.)  # conversion from Vout to Iin, the first 1000 is a factor from the CAEN data writing
#ADCconvQ = 1000*(2.25/16384.)/(gainq*11*2)
ADCconv = (1)/(gain*1200.*10.*2)
#ADCconvD = (1)/(gainD*1200.*10.*2)
scale = 1e-6

time = 0
tt = 0
trigE = 0
trig = 0
tbin = 1

sf = (1/(samplelength*1e-9))


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

for filen in files:
  nEf = 0
  f =  open(filen,'r')
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
  if nE%100 == 0: print(nE)
  if nE==20000: break
  meanA = pA[0:400].mean()
  meanB = pB[0:400].mean()
  meanC = pC[0:400].mean()
  meanD = pD[0:400].mean()

  pA -= meanA
  pB -= meanB
  pC -= meanC
  pD -= meanD

  pA *= ADCconv
  pB *= ADCconv
  pC *= ADCconv
  pD *= ADCconv

  stdA = pA[0:nbin_pre].std()
  stdB = pB[0:nbin_pre].std()
  stdC = pC[0:nbin_pre].std()
  stdD = pD[0:nbin_pre].std()

  stdA2 = pA[nSamples-nbin_post:nSamples].std()
  stdB2 = pB[nSamples-nbin_post:nSamples].std()
  stdC2 = pC[nSamples-nbin_post:nSamples].std()
  stdD2 = pD[nSamples-nbin_post:nSamples].std()

  minA = pA.min()
  minB = pB.min()
  minC = pC.min()
  minD = pD.min()

  maxtailA = pA[nSamples-nbin_post:nSamples].max()
  maxtailB = pB[nSamples-nbin_post:nSamples].max()
  maxtailC = pC[nSamples-nbin_post:nSamples].max()
  maxtailD = pD[nSamples-nbin_post:nSamples].max()

  fitA = False
  fitB = False
  fitC = False
  fitD = False

  #if trigger != '10000000000':
  if True:
      std_A.Fill(stdA)
      std_B.Fill(stdB)
      std_C.Fill(stdC)
      std_D.Fill(stdD)

      stdp_A.Fill(stdA2)
      stdp_B.Fill(stdB2)
      stdp_C.Fill(stdC2)
      stdp_D.Fill(stdD2)

      max_A.Fill(pA.max())
      max_B.Fill(pB.max())
      max_C.Fill(pC.max())
      max_D.Fill(pD.max())

      min_A.Fill(pA.min())
      min_B.Fill(pB.min())
      min_C.Fill(pC.min())
      min_D.Fill(pD.min())

      max_A_tail.Fill(maxtailA)
      max_B_tail.Fill(maxtailB)
      max_C_tail.Fill(maxtailC)
      max_D_tail.Fill(maxtailD)

      ra1 = [0.01 * scale, 0.1 * scale]
      ra2 = [.1 * scale, .2 * scale]
      ra3 = [.2 * scale, .3 * scale]
      ra4 = [.3 * scale, .4 * scale]
      ra5 = [.4 * scale, .5 * scale]

      raa1 = [0 * scale, 1 * scale]
      raa2 = [1 * scale, 1.5 * scale]
      raa3 = [1.5 * scale, 4 * scale]
      raa4 = [4 * scale, 6 * scale]
      raa5 = [6 * scale, 9 * scale]

      re1 = [0.1 * scale, .6 * scale]
      re2 = [0.7 * scale, 1. * scale]
      re3 = [1. * scale, 1.3 * scale]
      re4 = [1.3 * scale, 1.6 * scale]
      re5 = [1.6 * scale, 2 * scale]

      '''
      traceA = ROOT.TGraph(nSamples,time2,pA)
      traceA.SetName("traceA_"+str(nE))
      traceA.Write()
      traceB = ROOT.TGraph(nSamples,time2,pB)
      traceB.SetName("traceB_"+str(nE))
      traceB.Write()
      traceC = ROOT.TGraph(nSamples,time2,pC)
      traceC.SetName("traceC_"+str(nE))
      traceC.Write()
      traceD = ROOT.TGraph(nSamples,time2,pD)
      traceD.SetName("traceD_"+str(nE))
      traceD.Write()
      '''

      if stdA < 0.025 * scale and stdA > 0 * scale and stdA2 < .025 * scale and stdA2 > 0 * scale and minA > -0.03 * scale and maxtailA < 0.1 * scale:
          if pA.max() > ra1[0] and pA.max() < ra1[1]:
              pT1A += pA
          if pA.max() > ra2[0] and pA.max() < ra2[1]:
              pT2A += pA
          if pA.max() > ra3[0] and pA.max() < ra3[1]:
              pT3A += pA
          if pA.max() > ra4[0] and pA.max() < ra4[1]:
              pT4A += pA
          if pA.max() > ra5[0] and pA.max() < ra5[1]:
              pT5A += pA

      if stdB < .01 * scale and stdB > 0 and stdB2 < .01 * scale and stdB2 > 0 and minB > -0.03 * scale and maxtailB < 0.2 * scale:
          if pB.max() > ra1[0] and pB.max() < ra1[1]:
              pT1B += pB
          if pB.max() > ra2[0] and pB.max() < ra2[1]:
              pT2B += pB
          if pB.max() > ra3[0] and pB.max() < ra3[1]:
              pT3B += pB
          if pB.max() > ra4[0] and pB.max() < ra4[1]:
              pT4B += pB
          if pB.max() > ra5[0] and pB.max() < ra5[1]:
              pT5B += pB

      if stdC < .01 * scale and stdC > 0 and stdC2 < .01 * scale and stdC2 > 0 and minC > -0.03 * scale and maxtailC < 0.2 * scale:
          if pC.max() > ra1[0] and pC.max() < ra1[1]:
              pT1C += pC
          if pC.max() > ra2[0] and pC.max() < ra2[1]:
              pT2C += pC
          if pC.max() > ra3[0] and pC.max() < ra3[1]:
              pT3C += pC
          if pC.max() > ra4[0] and pC.max() < ra4[1]:
              pT4C += pC
          if pC.max() > ra5[0] and pC.max() < ra5[1]:
              pT5C += pC
      if stdD < 0.01 * scale and stdD > 0 and stdD2 < 0.01 * scale and stdD2 > 0 and minD > -0.03 * scale and maxtailD < 0.1 * scale:
          if pD.max() > ra1[0] and pD.max() < ra1[1]:
              pT1D += pD
          if pD.max() > ra2[0] and pD.max() < ra2[1]:
              pT2D += pD
          if pD.max() > ra3[0] and pD.max() < ra3[1]:
              pT3D += pD
          if pD.max() > ra4[0] and pD.max() < ra4[1]:
              pT4D += pD
          if pD.max() > ra5[0] and pD.max() < ra5[1]:
              pT5D += pD
  f.close()roo

pT1A /= pT1A.max()
pT1B /= pT1B.max()
pT1C /= pT1C.max()
pT1D /= pT1D.max()

pT2A /= pT2A.max()
pT2B /= pT2B.max()
pT2C /= pT2C.max()
pT2D /= pT2D.max()

pT3A /= pT3A.max()
pT3B /= pT3B.max()
pT3C /= pT3C.max()
pT3D /= pT3D.max()

pT4A /= pT4A.max()
pT4B /= pT4B.max()
pT4C /= pT4C.max()
pT4D /= pT4D.max()

pT5A /= pT5A.max()
pT5B /= pT5B.max()
pT5C /= pT5C.max()
pT5D /= pT5D.max()

pulse1A = ROOT.TGraph(nSamples,time2,pT1A)
pulse1A.SetName("pulse1A_template")
pulse1A.Write()
pulse1B = ROOT.TGraph(nSamples,time2,pT1B)
pulse1B.SetName("pulse1B_template")
pulse1B.Write()
pulse1C = ROOT.TGraph(nSamples,time2,pT1C)
pulse1C.SetName("pulse1C_template")
pulse1C.Write()
pulse1D = ROOT.TGraph(nSamples,time2,pT1D)
pulse1D.SetName("pulse1D_template")
pulse1D.Write()

pulse2A = ROOT.TGraph(nSamples,time2,pT2A)
pulse2A.SetName("pulse2A_template")
pulse2A.Write()
pulse2B = ROOT.TGraph(nSamples,time2,pT2B)
pulse2B.SetName("pulse2B_template")
pulse2B.Write()
pulse2C = ROOT.TGraph(nSamples,time2,pT2C)
pulse2C.SetName("pulse2C_template")
pulse2C.Write()
pulse2D = ROOT.TGraph(nSamples,time2,pT2D)
pulse2D.SetName("pulse2D_template")
pulse2D.Write()

pulse3A = ROOT.TGraph(nSamples,time2,pT3A)
pulse3A.SetName("pulse3A_template")
pulse3A.Write()
pulse3B = ROOT.TGraph(nSamples,time2,pT3B)
pulse3B.SetName("pulse3B_template")
pulse3B.Write()
pulse3C = ROOT.TGraph(nSamples,time2,pT3C)
pulse3C.SetName("pulse3C_template")
pulse3C.Write()
pulse3D = ROOT.TGraph(nSamples,time2,pT3D)
pulse3D.SetName("pulse3D_template")
pulse3D.Write()

pulse4A = ROOT.TGraph(nSamples,time2,pT4A)
pulse4A.SetName("pulse4A_template")
pulse4A.Write()
pulse4B = ROOT.TGraph(nSamples,time2,pT4B)
pulse4B.SetName("pulse4B_template")
pulse4B.Write()
pulse4C = ROOT.TGraph(nSamples,time2,pT4C)
pulse4C.SetName("pulse4C_template")
pulse4C.Write()
pulse4D = ROOT.TGraph(nSamples,time2,pT4D)
pulse4D.SetName("pulse4D_template")
pulse4D.Write()

pulse5A = ROOT.TGraph(nSamples,time2,pT5A)
pulse5A.SetName("pulse5A_template")
pulse5A.Write()
pulse5B = ROOT.TGraph(nSamples,time2,pT5B)
pulse5B.SetName("pulse5B_template")
pulse5B.Write()
pulse5C = ROOT.TGraph(nSamples,time2,pT5C)
pulse5C.SetName("pulse5C_template")
pulse5C.Write()
pulse5D = ROOT.TGraph(nSamples,time2,pT5D)
pulse5D.SetName("pulse5D_template")
pulse5D.Write()


saveF.Write()

print(nE)
