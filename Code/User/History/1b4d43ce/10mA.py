#! /usr/bin/env python3
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

def butter_bandpass(lowcut, fs, order=5):
    low = lowcut / (0.5 * fs)
    b, a = butter(order, low, btype='lowpass')
    return b, a1


def butter_bandpass_filter(data, lowcut, fs, order=5):
    b, a = butter_bandpass(lowcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

#files = glob.glob("Run37_4mmAl2O3_NoTrig_DG50_PTon_pm1V_noise.lvm")
#files = glob.glob("/data/sandro/Run39_Al2O3coins1and2_trigNA_DG50_PTon_QET40_Noise.lvm")
files = glob.glob("/data2/MPHY_Run43/902/Run43_0721_slot5_Sapphire902-ABCD_DG50_QET50pct_NoProbeTune_PTon_Noise*.txt")

saveF = ROOT.TFile("MPHY_Run43/902/Run43_0721_slot5_Sapphire902-ABCD_DG50_QET50pct_NoProbeTune_PTon_Noise.root","RECREATE")
cutoff=8000

std_A = ROOT.TH1D("std_A"," ; ",1000,0,50e-8)
std_B = ROOT.TH1D("std_B"," ; ",1000,0,50e-8)
std_C = ROOT.TH1D("std_C"," ; ",1000,0,50e-8)
std_D = ROOT.TH1D("std_D"," ; ",1000,0,50e-8)

max_A = ROOT.TH1D("max_A"," ; ",1000,0,1e-6)
max_B = ROOT.TH1D("max_B"," ; ",1000,0,1e-6)
max_C = ROOT.TH1D("max_C"," ; ",1000,0,1e-6)
max_D = ROOT.TH1D("max_D"," ; ",1000,0,1e-6)

gain = 50.
#gainD = 50.

nSamples = 9500
nSamples2 = 1900
#volt = 0.5 #in Volts

time = np.arange(nSamples,dtype=np.float64)
samplelength = 1000
time2 = time*samplelength
sf = (1/(samplelength*1e-9))

timeSub = np.arange(nSamples2,dtype=np.float64)
time2Sub = timeSub*samplelength

#sinw = np.linspace(0, ((samplelength*nSamples2*1e-9)/(1./100000.))*2*np.pi, nSamples2)
#sinwave = np.sin(sinw)

#freq = np.arange((nSamples/2)+1,dtype=np.float64)
#freq2 = freq*sf

#input ChC in mVs other Chs in mV
#ADCconvV = 1000*((2.25/16384.))/(gain*10.*2)  # conversion from Vout to Iin, the first 1000 is a factor from the CAEN data writing
ADCconv = (1)/(gain*1200.*10.*2)  # conversion from Vout to Iin, the first 1000 is a factor from the CAEN data writing and x3 from Matt's finding that CAEN is changing the signal
#ADCconvQ = 1000*(2.25/16384.)/(gainq*11*2)
#ADCconvA = (1)/(gainD*1200.*10*2)

time = 0
tt = 0
trigE = 0
trig = 0
tbin = 1

pNA = np.zeros(nSamples)
pNB = np.zeros(nSamples)
pNC = np.zeros(nSamples)
pND = np.zeros(nSamples)


nNoise = len(files)
countr=0
linenmbr=0
print(nNoise)
nNoise = nNoise*math.floor(nSamples/nSamples2)#Multiply the number of noise traces to the chuncks
s = (nNoise,4,nSamples2)
traces_Si17 = np.zeros(s)
noise_Si17 = np.zeros(s)

nE = 0
nN = 0
nEf = 0
countr=0
for filen in files:
  if countr>cutoff:
     break
  f =  open(filen,'r')
  linenmbr = 0

  for line in f:
      linenmbr = linenmbr + 1
      if (linenmbr > 3 and linenmbr < nSamples + 4):
        stuff = line.split()
        try:
            pNA[linenmbr-4] = float(stuff[1])
        except:
            pNA[linenmbr-4] = 10e6*ADCconv
        try:
            pNB[linenmbr-4] = float(stuff[2])
        except:
            pNB[linenmbr-4] = 10e6*ADCconv
            print(1)
        try:
            pNC[linenmbr-4] = float(stuff[3])
        except:
            pNC[linenmbr-4] = 10e6*ADCconv
        try:
            pND[linenmbr-4] = float(stuff[4])
        except:
            pND[linenmbr-4] = 10e6*ADCconv

  nE += 1
  nEf += 1
  countr = countr + 1
  if nE%100 == 0: print(nE)
  meanA = pNA[0:400].mean()
  meanB = pNB[0:400].mean()
  meanC = pNC[0:400].mean()
  meanD = pND[0:400].mean()


  pNA -= meanA
  pNB -= meanB
  pNC -= meanC
  pND -= meanD


  pNA *= ADCconv
  pNB *= ADCconv
  pNC *= ADCconv
  pND *= ADCconv


  for j in range(math.floor(len(pNA)/nSamples2)):
    nN += 1
    noise_Si17[nN-1][0] = pNA[nSamples2*j:nSamples2*(j+1)]
    noise_Si17[nN-1][1] = pNB[nSamples2*j:nSamples2*(j+1)]
    noise_Si17[nN-1][2] = pNC[nSamples2*j:nSamples2*(j+1)]
    noise_Si17[nN-1][3] = pND[nSamples2*j:nSamples2*(j+1)]
    #print(pNA[nSamples2*j:nSamples2*(j+1)])

  f.close()
#plot = noise.plot_psd(lgcoverlay=True)
chANoise = np.zeros(nSamples2)
chBNoise = np.zeros(nSamples2)
chCNoise = np.zeros(nSamples2)
chDNoise = np.zeros(nSamples2)


chAPSD = np.zeros(int(nSamples2/2))
chBPSD = np.zeros(int(nSamples2/2))
chCPSD = np.zeros(int(nSamples2/2))
chDPSD = np.zeros(int(nSamples2/2))


freq1, throw1 = qp.calc_psd(chANoise,sf, folded_over=False)
freq2 = freq1[0:nSamples2-1]
#print(freq2)

nGoodA = 0
nGoodB = 0
nGoodC = 0
nGoodD = 0


for x in range(nNoise-1):
    passA = False
    passB = False
    passC = False
    passD = False

    tA = noise_Si17[x,0]
    tB = noise_Si17[x,1]
    tC = noise_Si17[x,2]
    tD = noise_Si17[x,3]

    #tA = butter_bandpass_filter(tAf, 150000, sf, order=2)
    #tB = butter_bandpass_filter(tBf, 150000, sf, order=2)
    
    std_A.Fill(tA.std())
    std_B.Fill(tB.std())
    std_C.Fill(tC.std())
    std_D.Fill(tD.std())

    max_A.Fill(tA.max())
    max_B.Fill(tB.max())
    max_C.Fill(tC.max())
    max_D.Fill(tD.max())

    if 0.0 not in tA and tA.std() < 8e-9 and tA.std() > 3e-9 and tA.max() < 0.02e-6: passA = True
    if 0.0 not in tB and tB.std() < 7.5e-9 and tB.std() > 2e-9 and tB.max() < 0.02e-6: passB = True
    if 0.0 not in tC and tC.std() < 8.5e-9 and tC.std() > 4e-9 and tC.max() < 0.025e-6: passC = True
    if 0.0 not in tD and tD.std() < 7e-9 and tD.std() > 2e-9 and tD.max() < 0.02e-6: passD = True


    '''
    traceA = ROOT.TGraph(nSamples2,time2Sub,tA)
    traceA.SetName("traceA_"+str(x))
    if passA: traceA.Write()
    traceB = ROOT.TGraph(nSamples2,time2Sub,tB)
    traceB.SetName("traceB_"+str(x))
    if passB: traceB.Write()
    traceC = ROOT.TGraph(nSamples2,time2Sub,tC)
    traceC.SetName("traceC_"+str(x))
    if passC: traceC.Write()
    traceD = ROOT.TGraph(nSamples2,time2Sub,tD)
    traceD.SetName("traceD_"+str(x))
    if passD: traceD.Write()
    '''
 
    if passA:
       chAThrowaway,chANoiseC = qp.calc_psd(tA,sf, folded_over=False)
       chAThrow,chAPSDC = qp.calc_psd(tA,sf, folded_over=True)
       #chAThrow,chAPSDC = qp.calc_psd(sinwave,sf, folded_over=False)
       chAPSD += chAPSDC[0:int(nSamples2/2)]
       chANoise += chANoiseC
       nGoodA += 1
    if passB:
       chBThrowaway,chBNoiseC = qp.calc_psd(tB,sf, folded_over=False)
       chBThrow,chBPSDC = qp.calc_psd(tB,sf, folded_over=True)
       chBPSD += chBPSDC[0:int(nSamples2/2)]
       chBNoise += chBNoiseC
       nGoodB += 1
    if passC:
       chCThrowaway,chCNoiseC = qp.calc_psd(tC,sf, folded_over=False)
       chCThrow,chCPSDC = qp.calc_psd(tC,sf, folded_over=True)
       chCPSD += chCPSDC[0:int(nSamples2/2)]
       chCNoise += chCNoiseC
       nGoodC += 1
    if passD:
       chDThrowaway,chDNoiseC = qp.calc_psd(tD,sf, folded_over=False)
       chDThrow,chDPSDC = qp.calc_psd(tD,sf, folded_over=True)
       chDPSD += chDPSDC[0:int(nSamples2/2)]
       chDNoise += chDNoiseC
       nGoodD += 1


print(float(nGoodA))
print(float(nGoodB))
print(float(nGoodC))
print(float(nGoodD))

chANoise = chANoise/float(nGoodA)
chBNoise = chBNoise/float(nGoodB)
chCNoise = chCNoise/float(nGoodC)
chDNoise = chDNoise/float(nGoodD)

chAPSD = chAPSD/float(nGoodA)
chAPSD = np.sqrt(chAPSD)
chBPSD = chBPSD/float(nGoodB)
chBPSD = np.sqrt(chBPSD)
chCPSD = chCPSD/float(nGoodC)
chCPSD = np.sqrt(chCPSD)
chDPSD = chDPSD/float(nGoodD)
chDPSD = np.sqrt(chDPSD)


noiseA = ROOT.TGraph(nSamples2,freq1,chANoise)
noiseA.SetName("noiseA_template")
noiseA.Write()

noiseB = ROOT.TGraph(nSamples2,freq1,chBNoise)
noiseB.SetName("noiseB_template")
noiseB.Write()

noiseC = ROOT.TGraph(nSamples2,freq1,chCNoise)
noiseC.SetName("noiseC_template")
noiseC.Write()

noiseD = ROOT.TGraph(nSamples2,freq1,chDNoise)
noiseD.SetName("noiseD_template")
noiseD.Write()



PSDA = ROOT.TGraph(int(nSamples2/2),freq2,chAPSD)
PSDA.SetName("PSDA")
PSDA.Write()
PSDB = ROOT.TGraph(int(nSamples2/2),freq2,chBPSD)
PSDB.SetName("PSDB")
PSDB.Write()
PSDC = ROOT.TGraph(int(nSamples2/2),freq2,chCPSD)
PSDC.SetName("PSDC")
PSDC.Write()
PSDD = ROOT.TGraph(int(nSamples2/2),freq2,chDPSD)
PSDD.SetName("PSDD")
PSDD.Write()

saveF.Write()
print(nE)

'''
c1 = ROOT.TCanvas("c1","",800,800)
c1.cd()
PSDA.GetXaxis().SetTitle("Frequency (Hz)")
PSDA.GetYaxis().SetTitle("A/#sqrt{Hz}")
PSDA.SetLineColor(4)
PSDB.SetLineColor(2)
PSDC.SetLineColor(6)
PSDD.SetLineColor(7)
PSDA.SetMaximum(1E-8)
PSDA.SetMinimum(1E-14)
PSDA.GetXaxis().SetLimits(600, 5E5);
ROOT.gPad.SetLogy();
ROOT.gPad.SetLogx();
PSDA.Draw()
PSDB.Draw("same")
PSDC.Draw("same")
PSDD.Draw("same")
leg1 = ROOT.TLegend(0.25,0.2,0.45,0.48)
leg1.SetTextSize(0.06)
leg1.SetBorderSize(0)
leg1.AddEntry(PSDA,"Channel A","l")
leg1.AddEntry(PSDB,"Channel B","l")
leg1.AddEntry(PSDC,"Channel C","l")
leg1.AddEntry(PSDD,"Channel D","l")
leg1.Draw("same")
#c1.SaveAs("noisepsd.png","png")
'''

