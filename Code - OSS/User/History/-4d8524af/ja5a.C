#include "TH1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TFile.h"
#include "TH2.h"

void readData(){
TFile * f = new TFile("Run43_0721_Sapphire902-ABCD_slot5_DG50_QET50pct_PTon_Pulse_fitPhonon.root");
TTree * t=(TTree*)f->Get("data");


Double_t amp1,amp2,amp3,amp4,chi21,chi22,chi23,chi24,eventT;			

TH1F* sapphire = new TH1F("sapphire","",150,0e-6, 5e-06);
TH1F* sapphire_calibrated = new TH1F("sapphire_calibrated","",150,0, 100);

TH1F* sapphireA = new TH1F("sapphireA","ChA",100, 0e-6, 1e-06);
TH1F* sapphireB = new TH1F("sapphireB","ChB",100, 0e-6, 1e-06);
TH1F* sapphireC= new TH1F("sapphireC","ChC",200, 0e-6, 2.5e-06);
TH1F* sapphireD= new TH1F("sapphireD","ChD",100, 0e-6, 1.5e-06);

TH2F* chiA_vs_amp = new TH2F("chiA_vs_amp","ChA;a.u.;#chi^2",1000,0e-6, 3e-06, 1000,0, 16000);
TH2F* chiB_vs_amp = new TH2F("chiB_vs_amp","ChB;a.u.;#chi^2",1000,0e-6, 2e-06, 1000,0, 16000);
TH2F* chiC_vs_amp = new TH2F("chiC_vs_amp","ChC;a.u.;#chi^2",1000,0e-6, 1.2e-06, 1000,0, 16000);
TH2F* chiD_vs_amp = new TH2F("chiD_vs_amp","ChD;a.u.;#chi^2",1000,0e-6, 1.2e-06, 1000,0, 16000);

int nEvents = t->GetEntries();

for(Int_t i=0; i<t->GetEntries(); i++){
  t->SetBranchAddress("ampOFA",&amp1);
  t->SetBranchAddress("ampOFB",&amp2);
  t->SetBranchAddress("ampOFC",&amp3);
  t->SetBranchAddress("ampOFD",&amp4);

  t->SetBranchAddress("chi2OFA",&chi21);
  t->SetBranchAddress("chi2OFB",&chi22);
  t->SetBranchAddress("chi2OFC",&chi23);
  t->SetBranchAddress("chi2OFD",&chi24);

  sapphire->Fill((amp1+amp2+amp3+amp4));
  sapphire_calibrated->Fill(60*(amp1+amp2+amp3+amp4)/3.408e-06);
 
  sapphireA->Fill(amp1);
  sapphireB->Fill(amp2);
  sapphireC->Fill(amp3);
  sapphireD->Fill(amp4);

  chiA_vs_amp->Fill(amp1,chi21);
  chiB_vs_amp->Fill(amp2,chi22);
  chiC_vs_amp->Fill(amp3,chi23);
  chiD_vs_amp->Fill(amp4,chi24);
 
  t->GetEntry(i);
  }

TCanvas *c1 = new TCanvas("c1", "", 500, 500);
c1->cd();

sapphire->GetXaxis()->SetTitle("(a.u.)");
sapphire->GetYaxis()->SetTitle("Counts");
sapphire->Draw();

TCanvas *c11 = new TCanvas("c11", "", 500, 500);
c11->cd();

sapphire_calibrated->GetXaxis()->SetTitle("(KeV)");
sapphire_calibrated->GetYaxis()->SetTitle("Counts");
sapphire_calibrated->Draw();

TCanvas *c22 = new TCanvas("c22", "", 900, 900);
c22->Divide(2,2);
c22->cd(1);
sapphireA->GetXaxis()->SetTitle("a.u.");
sapphireA->GetYaxis()->SetTitle("Counts");
sapphireA->GetXaxis()->SetNdivisions(10);
sapphireA->Draw();

c22->cd(2);
sapphireB->GetXaxis()->SetTitle("a.u.");
sapphireB->GetYaxis()->SetTitle("Counts");
sapphireB->GetXaxis()->SetNdivisions(8);
sapphireB->Draw();

c22->cd(3);
sapphireC->GetXaxis()->SetTitle("a.u.");
sapphireC->GetYaxis()->SetTitle("Counts");
sapphireC->Draw();

c22->cd(4);
sapphireD->GetXaxis()->SetTitle("a.u.");
sapphireD->GetYaxis()->SetTitle("Counts");
sapphireD->Draw();

TCanvas *c555 = new TCanvas("c555", "", 1000, 900);
c555->Divide(2,2);
c555->cd(1);

chiA_vs_amp->GetXaxis()->SetTitleSize(0.03);
chiA_vs_amp->Draw("COLZ");
TF1 *par1 = new TF1("par1","y=2200+2e15*x*x",0,0.3e-6);
par1->SetLineColor(kRed);
par1->Draw("SAME");
c555->cd(2);
chiB_vs_amp->GetXaxis()->SetTitleSize(0.03);
chiB_vs_amp->Draw("COLZ");
TF1 *par2 = new TF1("par2","y=2200+8.555e14*x*x",0,0.8e-6);
par2->SetLineColor(kRed);
par2->Draw("SAME");
c555->cd(3);
chiC_vs_amp->GetXaxis()->SetTitleSize(0.03);
chiC_vs_amp->Draw("COLZ");
TF1 *par3 = new TF1("par3","y=2100+5.5e15*x*x",0,0.6e-6);
par3->SetLineColor(kRed);
par3->Draw("SAME");
c555->cd(4);
chiD_vs_amp->GetXaxis()->SetTitleSize(0.03);
chiD_vs_amp->Draw("COLZ");
TF1 *par4 = new TF1("par3","y=2100+5.55e15*x*x",0,0.6e-6);
par4->SetLineColor(kRed);
par4->Draw("SAME");

}




    