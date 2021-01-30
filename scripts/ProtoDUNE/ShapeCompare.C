////////////////////////////////////////////////////
// CERN ROOT Script to extract some radiological  //
// candidate signals and plot their average for   //
// various energy values.                         //
// AFlesher                                       //
////////////////////////////////////////////////////

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

double bline(double x1, double y1, double x2, double y2, double x){
  double m, b;

  m = (y1 - y2) / (x1 - x2);
  b = -(m * x1) + y1;

  return ((m * x) + b);
}

double GetMax(TH1F*& h){
  return h->GetBinContent(h->GetMaximumBin());
}

void ShapeCompareMC(){
  TLegend *l = new TLegend(0.6, 0.67, 0.87, 0.87);
  TCanvas *c1 = new TCanvas("c1", "c1");
  //c1->cd();

  TFile *f1 = TFile::Open("/path/to/merged/output.root");

  TTree *ctree = (TTree*)f1->Get("ctree");

  const double F = 138.0;  // Gain
  // const double F = 200.0;  // Gain
  const double I = 0.0236; // Ionization Energy for Ar
  const double R = 0.66;   // Recombination 
  const double K = 5.52;   // Electronics Response Area/Amplitude
  // const double K = 5.02;   // Electronics Response Area/Amplitude
  const float interval = 50.0;
  vector<float> intervals;

  double xWidth, yWidth, channel, charge, energy, maxAmp, confirmed;
  unsigned int candidates = 0, candidates2 = 0; // ev = etree->GetEntries(),

  // ctree->SetBranchAddress("ChargeVec", &ChargeVec);
  // ctree->SetBranchAddress("xWind", &xWidth);
  // ctree->SetBranchAddress("yWind", &yWidth);
  // ctree->SetBranchAddress("channel", &channel);
  // ctree->SetBranchAddress("charge", &charge);
  // ctree->SetBranchAddress("maxAmp", &maxAmp);
  // ctree->SetBranchAddress("confirmed", &confirmed);

  //TTreeReader cr("ctree", cfile);

  gStyle->SetOptStat(0);
 
  TTreeReader cr(ctree);

  
  TTreeReaderValue< vector<int> > ChargeVec(cr, "ChargeVec");
  TTreeReaderValue<double> x(cr, "xWind");
  TTreeReaderValue<double> y(cr, "yWind");
  TTreeReaderValue<double> ch(cr, "channel");
  TTreeReaderValue<double> c(cr, "charge");
  TTreeReaderValue<double> m(cr, "maxAmp");
  TTreeReaderValue<double> co(cr, "confirmed");


  TH2I *cHist = new TH2I("cHist", "Aggregated Integration Window", 
			 10, -5.5, 4.5, 
			 100, -49.5, 50.5);
  TH1F *m2 = new TH1F("m2", "-2", 100, -49.5, 50.5);
  TH1F *m1 = new TH1F("m1", "-1", 100, -49.5, 50.5);
  TH1F *central = new TH1F("central", "0", 100, -49.5, 50.5);
  TH1F *p1 = new TH1F("p1", "+1", 100, -49.5, 50.5);
  TH1F *p2 = new TH1F("p2", "+2", 100, -49.5, 50.5);
  // TH1F *chn = new TH1F("chn", "", 15360, 0, 15359);

  TH1F *m1_2 = new TH1F("m1_2", "-1", 100, -49.5, 50.5);
  TH1F *central_2 = new TH1F("central_2", "0", 100, -49.5, 50.5);
  TH1F *p1_2 = new TH1F("p1_2", "+1", 100, -49.5, 50.5);

  for (int i = 0; i <= 36; i++){
    // for (int i = 0; i <= 1; i++){
    intervals.push_back((float)(i * interval));
  }

  double median, offset_central_p, offset_central_m, offset_m_p, offset_m_m, offset_p_p, offset_p_m;
  vector<double> medvec;
  vector<unsigned int> chhist(15360, 0);

  vector<unsigned short> cskip;


  ////// Skip some noisy channels based on having much higher candidate number than other channels
  // while (cr.Next()){
  //   // chn->Fill(*ch);
  //   chhist[(unsigned int)*ch]++;
  // }

  // cr.Restart();

  // for (unsigned short i = 0; i < 15360; i++){
  //   //medvec.push_back(chn->GetBinContent(i));
  //   if (chhist[i] != 0)    
  //     medvec.push_back(chhist.at(i));
    
  //   // chn->Fill(i, chhist[i]);
  // }

  // sort(medvec.begin(), medvec.end());
  
  // median = 0.5 * (medvec[(unsigned int)((medvec.size() - 1) / 2)] + medvec[(unsigned int)(medvec.size() / 2)]);

  // cout << "Max: " << medvec.back() << ", Min: " << medvec[0] << ", Median: " << median << endl;

  // for (unsigned short i = 1; i <= chn->GetXaxis()->GetNbins(); i++){
  //   if (chn->GetBinContent(i) > 1.5 * median || chn->GetBinContent(i) < 0.5 * median){
  //     cskip.push_back(i - 1);
  //   }
  // }

  // sort(cskip.begin(), cskip.end());



  central->SetLineWidth(2);
  central->SetLineColor(kBlue);

  m1->SetLineWidth(2);
  m1->SetLineColor(kRed);

  p1->SetLineWidth(2);
  p1->SetLineColor(kGreen + 3);

  central_2->SetLineWidth(2);
  m1_2->SetLineWidth(2);
  p1_2->SetLineWidth(2);

  central_2->SetLineStyle(2);
  m1_2->SetLineStyle(2);
  p1_2->SetLineStyle(2);

  central_2->SetLineColor(kBlack);
  m1_2->SetLineColor(kMagenta);
  p1_2->SetLineColor(kCyan + 2);
  

  l->AddEntry(central, "Central Signal new LArG4", "L");
  l->AddEntry(m1, "-1 Signal new LArG4", "L");
  l->AddEntry(p1, "+1 Signal new LArG4", "L");
  l->AddEntry(central_2, "Central Signal Lower Range Cut", "L");
  l->AddEntry(m1_2, "-1 Signal Lower Range Cut", "L");
  l->AddEntry(p1_2, "+1 Signal Lower Range Cut", "L");


  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetShadowColor(0);
  //l->SetTextFont(62);
  l->SetTextSize(0.03);

  for (unsigned short i = 0; i < intervals.size() - 1; i++){

    while (cr.Next()){
      // ctree->GetEntry(i);
      xWidth = *x, yWidth = *y, channel = *ch, charge = *c, maxAmp = *m, confirmed = *c;
    

      // charge = charge * F * I / K;
      // energy = charge / recombGraph_vsEfree->Eval(charge);
      // energy = charge / R;
      energy = charge * 0.0234 * 23.6 / R;

      // if (maxAmp <= 20){continue;}
    

      if (energy >= intervals[i] && energy < intervals[i + 1]){// &&
	  // (double)chhist.at((unsigned int)channel) <= 1.5 * median && (double)chhist.at((unsigned int)channel) >= 0.5 * median){
	// cout << "test" << endl;
	candidates++;
	int ci = 0;

	for (int x = -xWidth; x <= xWidth; x++){
	  for (int y = -yWidth; y <= yWidth; y++){
	    if (ci < ChargeVec->size()){ //(1 + (2 * xWidth)) * (1 + (2 * yWidth))){
	      cHist->Fill(x, y, ChargeVec->at(ci));
	  
	      if (x == -2)
		m2->Fill(y, ChargeVec->at(ci));
	      else if (x == -1){
		m1->Fill(y, ChargeVec->at(ci));		
	      }
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
      }	
    }


    central->SetTitle(TString::Format("%.1f <= Energy < %.1f", intervals[i], intervals[i + 1])); 
    central->GetXaxis()->SetTitle("Time Tick");
    central->GetYaxis()->SetTitle("Average Charge [pdune ADC]");
    
    cout << intervals[i] << ", " << intervals[i + 1] << endl;

    // candidates = ctree->GetEntries(TString::Format("%f >= %f && %f < %f", energy, intervals[i], energy, intervals[i + 1]));
    // // candidates = 10;
    // cout << candidates << endl;
    // central->Scale(1.0 / central->Integral());
    // p1->Scale(1.0 / central->Integral());
    // m1->Scale(1.0 / central->Integral());

    // central->Scale(1.0 / (ev * interval));
    // m1->Scale(1.0 / (ev * interval));
    // p1->Scale(1.0 / (ev * interval));
    central->Scale(1.0 / candidates);
    m1->Scale(1.0 / candidates);
    p1->Scale(1.0 / candidates);

    central->GetYaxis()->SetRangeUser(-0.1 * GetMax(central_2), 1.15 * GetMax(central_2));

    central->Draw("hist");
    m1->Draw("histsame");
    p1->Draw("histsame");

    c1->SaveAs(TString::Format("./Signal_Shape/%i_Y.png", i));

    central->Reset();
    m1->Reset();
    p1->Reset();
    cr.Restart();

    candidates = 0;
    candidates2 = 0;
  }
}

