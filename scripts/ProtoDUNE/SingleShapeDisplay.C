////////////////////////////////////////////////////
// CERN ROOT Script to extract some individual    //
// radiological candidate signals.                //
// AFlesher                                       //
////////////////////////////////////////////////////

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"


void SingleShapeDisplay(){
  TFile *cfile = TFile::Open("/path/to/merged/output.root");
  TFile *efile = TFile::Open("/path/to/merged/output_e.root");

  //TFile f_output("Y.root", "RECREATE");
  //TFile f_output("UPlane_high.root", "RECREATE");
  //TFile f_output("VPlane_high.root", "RECREATE");

  TFile *recombFile = TFile::Open("./recombfile.root"); // use RecombPlot.C to generate this
  TGraph *recombGraph_vsEfree = (TGraph*)recombFile->Get("recombGraph_vsEfree");

  TLegend *l = new TLegend(0.6, 0.75, 0.87, 0.87);

  TTree *ctree = (TTree*)cfile->Get("ctree");
  TTree *etree = (TTree*)efile->Get("etree");

  TCanvas *c1 = new TCanvas("c1", "c1");

  //c1->cd();

  const double F = 138.0;  // Gain
  const double I = 0.0236; // Ionization Energy for Ar
  const double R = 0.66;   // Recombination 
  const double K = 5.52;   // Electronics Response Area/Amplitude


  double xWidth, yWidth, channel, time, charge, energy, maxAmp;
  unsigned int ev = etree->GetEntries(), candidates = 0, candidate, t;

  efile->Close();

  gStyle->SetOptStat(0);
 
  TTreeReader cr(ctree);
  
  TTreeReaderValue< vector<int> > ChargeVec(cr, "ChargeVec");
  TTreeReaderValue<double> x(cr, "xWind");
  TTreeReaderValue<double> y(cr, "yWind");
  TTreeReaderValue<double> ch(cr, "channel");
  TTreeReaderValue<double> ti(cr, "time");
  TTreeReaderValue<double> c(cr, "charge");
  TTreeReaderValue<double> e(cr, "energy");
  TTreeReaderValue<double> m(cr, "maxAmp");
  TTreeReaderValue<int> co(cr, "candidateIndex");

  TH2I *cHist = new TH2I("cHist", "Aggregated Integration Window", 
			 10, -5.5, 4.5, 
			 100, -49.5, 50.5);
  TH1F *m2 = new TH1F("m2", "-2", 100, -49.5, 50.5);
  TH1F *m1 = new TH1F("m1", "-1", 100, -49.5, 50.5);
  TH1F *central = new TH1F("central", "0", 100, -49.5, 50.5);
  TH1F *p1 = new TH1F("p1", "+1", 100, -49.5, 50.5);
  TH1F *p2 = new TH1F("p2", "+2", 100, -49.5, 50.5);


  central->SetLineWidth(2);
  central->SetLineColor(kBlue);

  m1->SetLineWidth(2);
  m1->SetLineStyle(2);
  m1->SetLineColor(kRed);

  p1->SetLineWidth(2);
  p1->SetLineStyle(2);
  p1->SetLineColor(kGreen + 3);
  

  l->AddEntry(central, "Central Signal", "L");
  l->AddEntry(m1, "-1 Signal", "L");
  l->AddEntry(p1, "+1 Signal", "L");


  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetShadowColor(0);
  //l->SetTextFont(62);
  l->SetTextSize(0.03);

  while (cr.Next()){
    // ctree->GetEntry(i);
    xWidth = *x, yWidth = *y, channel = *ch, time = *ti, charge = *c, energy = *e, maxAmp = *m, candidate = *co;
    
    // charge = charge * F * I / K;
    // energy = charge / recombGraph_vsEfree->Eval(charge);
   
    // if (t % 10 != 0) {continue;}
      
    t++;
      
    // if (candidates > 100) {continue;}

    if (energy > 700){
      int ci = 0;

      for (int x = -xWidth; x <= xWidth; x++){
	for (int y = -yWidth; y <= yWidth; y++){
	  if (ci < ChargeVec->size()){ //(1 + (2 * xWidth)) * (1 + (2 * yWidth))){
	    cHist->Fill(x, y, ChargeVec->at(ci));
	  
	    if (x == -2)
	      m2->Fill(y, ChargeVec->at(ci));
	    else if (x == -1)
	      m1->Fill(y, ChargeVec->at(ci));
	    else if (x == 0)
	      central->Fill(y, ChargeVec->at(ci));
	    else if (x == 1)
	      p1->Fill(y, ChargeVec->at(ci));
	    else if (x == 2)
	      p2->Fill(y, ChargeVec->at(ci));

	    ci++;
	  }
	}
      }

      central->SetTitle(TString::Format("Candidate %i", candidate)); 
      central->GetXaxis()->SetTitle("Time Tick");
      central->GetYaxis()->SetTitle("Charge [e^{-}]");

      central->Draw("hist");
      m1->Draw("histsame");
      p1->Draw("histsame");
      l->Draw("same");

      c1->SaveAs(TString::Format("./Single_Signal_Shape/%i_%i_%i_Y.png", candidate, (int)channel, (int)time));

      central->Reset();
      m1->Reset();
      p1->Reset();
    }
  }
}
