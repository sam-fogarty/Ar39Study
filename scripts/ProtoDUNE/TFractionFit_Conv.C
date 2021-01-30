#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TObject.h"

void WeierstrassTransform(TH1D*& tmp, TH1D*& h2, double s, unsigned short offset);
void RecombinationSmear(TH1D*& tmp, TH1D*& h2, double s, unsigned short offset);
vector<double> vGain(const string& filename);
double GetMax(TH1D*& h){
  return h->GetBinContent(h->GetMaximumBin());
}

void TFractionFit_Conv(){
  // const int bn = 5100; //bin number
  const int bn = 5000; //bin number
  const double maxE = 4183.1363; //max energy
  // const double maxE = 4417.42; //max energy
  const int linewidth = 3;
  const int linewidthsub = 2;
  const double R = 0.66;

  bool converged = false;

  vector<unsigned short> skip;
  unsigned short filenum;

  // vector<string> Nuclide = {"Ar39", "Co60", "Th232", "K40", "Kr85", "Rn222", "cosmic"};
  vector<string> Nuclide = {"Ar39", "Co60", "Th232", "K40", "Kr85", "U238", "cosmic"};
  // vector<string> Nuclide = {"Ar39", "Co60", "Th232", "K40", "Kr85", "U238", "cosmic"};
  vector<unsigned short> HColors = {kBlue + 4, kGreen + 2, kMagenta + 2, kCyan + 2, kBlue, kOrange + 3, kBlue - 10, kRed + 1, kBlack};

  // TFile *data = TFile::Open("Data/output.root");
  // TFile *cosmics = TFile::Open("MC/output.root");
  // // TFile *radsim = TFile::Open("MC/Ar39-5msLifetime.root");
  // TFile *radsim = TFile::Open("MC/Ar39.root");
  // TFile *radsimkr = TFile::Open("MC/Kr85.root");
  // TFile *radsimu = TFile::Open("MC/U238.root");
  // TFile *radsimk = TFile::Open("MC/K40.root");
  // TFile *radsimco = TFile::Open("MC/Co60.root");
  // TFile *radsimth = TFile::Open("MC/Th232.root");
  // TFile *radsimrn = TFile::Open("MC/Rn222.root");

  // //TFile fout("testf.root", "RECREATE");

  // TTree *d = (TTree*)data->Get("c_tree");
  // TTree *c = (TTree*)cosmics->Get("ctree");
  // TTree *r = (TTree*)radsim->Get("ctree");
  // TTree *kr = (TTree*)radsimkr->Get("ctree");
  // TTree *u = (TTree*)radsimu->Get("ctree");
  // TTree *th = (TTree*)radsimth->Get("ctree");
  // TTree *k40 = (TTree*)radsimk->Get("ctree");  
  // TTree *co = (TTree*)radsimco->Get("ctree");
  // TTree *rn = (TTree*)radsimrn->Get("ctree");

  // TFile *data = TFile::Open("FitFiles/data_5841.root");
  // TFile *cosmics = TFile::Open("FitFiles/cosmic.root");
  // TFile *radsim = TFile::Open("FitFiles/LRC_Radiological/Ar39.root");
  // TFile *radsimkr = TFile::Open("FitFiles/LRC_Radiological/Kr85.root");
  // TFile *radsimu = TFile::Open("FitFiles/LRC_Radiological/U238.root");
  // TFile *radsimk = TFile::Open("FitFiles/LRC_Radiological/K40.root");
  // TFile *radsimco = TFile::Open("FitFiles/LRC_Radiological/Co60.root");
  // TFile *radsimth = TFile::Open("FitFiles/LRC_Radiological/Th232.root");

  TFile *data = TFile::Open("FitFiles/data_5841.root");
  TFile *cosmics = TFile::Open("FitFiles/cosmic.root");
  TFile *radsim = TFile::Open("FitFiles/Ar39.root");
  TFile *radsimkr = TFile::Open("FitFiles/Kr85.root");
  TFile *radsimu = TFile::Open("FitFiles/U238.root");
  TFile *radsimk = TFile::Open("FitFiles/K40.root");
  TFile *radsimco = TFile::Open("FitFiles/Co60.root");
  TFile *radsimth = TFile::Open("FitFiles/Th232.root");
  // TFile *radsimrn = TFile::Open("FitFiles/Rn222.root");

  TTree *d = (TTree*)data->Get("ctree");
  TTree *c = (TTree*)cosmics->Get("ctree");
  TTree *ar = (TTree*)radsim->Get("ctree");
  TTree *kr = (TTree*)radsimkr->Get("ctree");
  TTree *u = (TTree*)radsimu->Get("ctree");
  TTree *k40 = (TTree*)radsimk->Get("ctree");
  TTree *co = (TTree*)radsimco->Get("ctree");
  TTree *th = (TTree*)radsimth->Get("ctree");
  // TTree *rn = (TTree*)radsimrn->Get("ctree");

  vector<double> data_gain = vGain("gain_by_channel.txt");
  vector<double> sim_gain = vGain("null");

  // TChain *ctree = new TChain("ctree");
  // TChain *etree = new TChain("etree");
  
  TH1D *reg = new TH1D("reg","",bn,-0.5,maxE);
  TH1D *mod = new TH1D("mod","",bn,-0.5,maxE);
  TH1D *sub = new TH1D("sub","",bn,-0.5,maxE);

  TH1D *cosm = new TH1D("cosm","",bn,-0.5,maxE);
  TH1D *radar = new TH1D("radar","",bn,-0.5,maxE);
  TH1D *radkr = new TH1D("radkr","",bn,-0.5,maxE);
  TH1D *radu = new TH1D("radu","",bn,-0.5,maxE);
  TH1D *radth = new TH1D("radth","",bn,-0.5,maxE);
  TH1D *radk = new TH1D("radk","",bn,-0.5,maxE);
  TH1D *radco = new TH1D("radco","",bn,-0.5,maxE);
  // TH1D *radrn = new TH1D("radrn","",bn,-0.5,maxE);

  TLegend* leg = new TLegend(0.53,0.4,0.85,0.8);

  TObjArray *MC = new TObjArray(4);

  TLatex ct;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //TCanvas *c1 = new TCanvas("c1","c1");
  //TCanvas *c2 = new TCanvas("c2","c2");

  // float dThresh = 18.0;
  // float dThresh = 8.0;

  double chisquare = 0, minchi = 9999999,  arx, ary = 1, cy = 1, kx, ky = 1, krx, kry = 1, ux, uy = 1, thx, thy = 1, cox, coy = 1, rnx, rny = 1, errors;
  vector<double> EvNum(Nuclide.size() + 1);

  double EThresh = 0.0;
  double median;
  vector<double> medvec;
  vector<unsigned int> chhist(15360, 0);

  TTreeReader RA(d);
  TTreeReaderValue<double> ChargeVal(RA, "charge");
  TTreeReaderValue<double> MaxAmp(RA, "maxAmp");
  TTreeReaderValue<double> nearTrack(RA, "nearTrack");
  TTreeReaderValue<double> channel(RA, "channel");

  while (RA.Next()){
    chhist[(unsigned int)*channel]++;
  }
  RA.Restart();

  for (unsigned short i = 0; i < 15360; i++){
    if (chhist[i] != 0)
      medvec.push_back(chhist.at(i));
  }
  sort(medvec.begin(), medvec.end());
  median = 0.5 * (medvec[(unsigned int)((medvec.size() - 1) / 2)] + medvec[(unsigned int)(medvec.size() / 2)]);

  while (RA.Next()){
    if (chhist[(int)*channel] < 1.5 * median && chhist[(int)*channel] > 0.5 * median) // data
      reg->Fill(1.047 * (*ChargeVal) * data_gain[(*channel)] * 23.6 / R);
  }
  RA.Restart();

  RA.SetTree(ar);
  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radar->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(c);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      cosm->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(k40);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radk->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(kr);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radkr->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(u);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radu->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(co);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radco->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(th);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radth->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  // RA.SetTree(rn);

  // while (RA.Next()){
  //   if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
  //     radrn->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  // }

  TFractionFitter *f1;
  unsigned char loopn = 0;


  float step = 0.01;
  for (unsigned char a = 1; a <= 5; a++){

    if (a == loopn){
      WeierstrassTransform(radar, radar, step, 1000);
      WeierstrassTransform(cosm, cosm, step, 1000);
      WeierstrassTransform(radk, radk, step, 1000);
      WeierstrassTransform(radkr, radkr, step, 1000);
      WeierstrassTransform(radu, radu, step, 1000);
      WeierstrassTransform(radth, radth, step, 1000);
      WeierstrassTransform(radco, radco, step, 1000);
      // WeierstrassTransform(radrn, radrn, step, 1000);

      // RecombinationSmear(radar, radar, step, 1000);
      // RecombinationSmear(cosm, cosm, step, 1000);
      // RecombinationSmear(radk, radk, step, 1000);
      // RecombinationSmear(radkr, radkr, step, 1000);
      // RecombinationSmear(radu, radu, step, 1000);
      // RecombinationSmear(radth, radth, step, 1000);
      // RecombinationSmear(radco, radco, step, 1000);
      // RecombinationSmear(radrn, radrn, step, 1000);
    }else {
      loopn = a;
      WeierstrassTransform(radar, radar, a * step, 1000);
      WeierstrassTransform(cosm, cosm, a * step, 1000);
      WeierstrassTransform(radk, radk, a * step, 1000);
      WeierstrassTransform(radkr, radkr, a * step, 1000);
      WeierstrassTransform(radu, radu, a * step, 1000);
      WeierstrassTransform(radth, radth, a * step, 1000);
      WeierstrassTransform(radco, radco, a * step, 1000);
      // WeierstrassTransform(radrn, radrn, a * step, 1000);

      // RecombinationSmear(radar, radar, a * step, 1000);
      // RecombinationSmear(cosm, cosm, a * step, 1000);
      // RecombinationSmear(radk, radk, a * step, 1000);
      // RecombinationSmear(radkr, radkr, a * step, 1000);
      // RecombinationSmear(radu, radu, a * step, 1000);
      // RecombinationSmear(radth, radth, a * step, 1000);
      // RecombinationSmear(radco, radco, a * step, 1000);
      // RecombinationSmear(radrn, radrn, a * step, 1000);

    }

    loopn++;

    MC->Clear();
    MC->Add(cosm);
    MC->Add(radar);
    MC->Add(radkr);
    MC->Add(radu);
    MC->Add(radth);
    MC->Add(radk);
    MC->Add(radco);
    // MC->Add(radrn);

    f1 = new TFractionFitter(reg, MC);

    f1->Constrain(0, 0.0, 1.0);
    f1->Constrain(1, 0.1, 1.0);
    f1->Constrain(2, 0.0, 1.0);
    f1->Constrain(3, 0.0, 1.0);
    // f1->Constrain(4, 0.1, 0.4);
    f1->Constrain(4, 0.0, 1.0);
    f1->Constrain(5, 0.0, 1.0);
    f1->Constrain(6, 0.0, 1.0);
    // f1->Constrain(7, 0.0, 0.01);

    // f1->SetRangeX(radar->GetXaxis()->FindBin(400),bn); //Set bins to be used in fitting
    // f1->SetRangeX(radar->GetXaxis()->FindBin(100),radar->GetXaxis()->FindBin(600)); //Set bins to be used in fitting

    Int_t status = f1->Fit();

    cout << "fit status: " << status << endl;

    if (status == 0) {
      converged = true;
      
      // mod = (TH1D*) f1->GetPlot();

      mod = (TH1D*) f1->GetPlot();

      // cout << "Degrees of Freedom: " << f1->GetNDF() << endl;

      // (*mod) = (*cosm) + (*radar) + (*radu) + (*radk) + (*radkr) + (*radth) + (*radco) + (*radrn);

      ////Manually calc Chi square                                                                                                                     
      for (int i = 1; i <= bn; i++){
	if (reg->GetBinContent(i) != 0)
	  chisquare += pow(reg->GetBinContent(i) - mod->GetBinContent(i), 2.0)/pow(reg->GetBinError(i), 2.0);
      }

      if (chisquare < minchi){
	minchi = chisquare;
	arx = a * step;
	f1->GetResult(0, cy, errors);
	f1->GetResult(1, ary, errors);
	f1->GetResult(2, kry, errors);
	f1->GetResult(3, uy, errors);
	f1->GetResult(4, thy, errors);
	f1->GetResult(5, ky, errors);
	f1->GetResult(6, coy, errors);
	// f1->GetResult(7, rny, errors);
      }
    }
  }
   

  radar->Reset();
  cosm->Reset();
  radk->Reset();
  radkr->Reset();
  radco->Reset();
  radth->Reset();
  radu->Reset();
  // radrn->Reset();

  RA.Restart();

  RA.SetTree(ar);
  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radar->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(c);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      cosm->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(k40);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radk->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(kr);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radkr->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(u);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radu->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(co);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radco->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  RA.SetTree(th);

  while (RA.Next()){
    if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
      radth->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  }
  RA.Restart();
  // RA.SetTree(rn);

  // while (RA.Next()){
  //   if (*ChargeVal * sim_gain[*channel] * 23.6 / R > EThresh)
  //     radrn->Fill(((*ChargeVal) * sim_gain[(*channel)] * 23.6 / R));
  // }

  WeierstrassTransform(radar, radar, arx, 1000);
  WeierstrassTransform(cosm, cosm, arx, 1000);
  WeierstrassTransform(radk, radk, arx, 1000);
  WeierstrassTransform(radkr, radkr, arx, 1000);
  WeierstrassTransform(radu, radu, arx, 1000);
  WeierstrassTransform(radth, radth, arx, 1000);
  WeierstrassTransform(radco, radco, arx, 1000);
  // WeierstrassTransform(radrn, radrn, arx, 1000);


  mod->Rebin(5);
  reg->Rebin(5);
  sub->Rebin(5);
  cosm->Rebin(5);
  radar->Rebin(5);
  radkr->Rebin(5);
  radk->Rebin(5);
  radth->Rebin(5);
  radu->Rebin(5);
  // radrn->Rebin(5);
  radco->Rebin(5);

  reg->GetXaxis()->SetTitle("Energy [keV]");
  reg->GetYaxis()->SetTitle("Arb. Units");
  
  reg->Scale(1.0/reg->Integral());
  //reg->GetYaxis()->SetRangeUser(0,sub->GetMaximum() + sub->GetMaximum()*0.1);
  // reg->GetYaxis()->SetRangeUser(0.1,reg->GetMaximum() + reg->GetMaximum()*0.1);
  
  reg->GetXaxis()->SetTitleFont(62);
  reg->GetYaxis()->SetTitleFont(62);
  
  reg->GetXaxis()->SetLabelFont(62);
  reg->GetYaxis()->SetLabelFont(62);

  reg->GetXaxis()->SetTitleSize(0.045);
  reg->GetYaxis()->SetTitleSize(0.045);

  reg->GetXaxis()->SetLabelSize(0.045);
  reg->GetYaxis()->SetLabelSize(0.045);

  reg->GetXaxis()->SetRangeUser(0,1000);

  reg->SetLineWidth(linewidth);
  reg->SetLineColor(kBlack);

  //   sub->SetLineColor(kBlack);
  //   sub->SetLineStyle(2);
  //   sub->SetLineWidth(linewidth);

  //   //mod->SetLineColor(kRed);
  //   //mod->Scale(100000.0/mod->Integral());


  //   f1->GetResult(0, cy, errors);
  //   f1->GetResult(1, ary, errors);
  //   f1->GetResult(2, kry, errors);
  //   f1->GetResult(3, uy, errors);
  //   f1->GetResult(4, thy, errors);
  //   f1->GetResult(5, ky, errors);
  //   f1->GetResult(6, coy, errors);
  //   f1->GetResult(7, rny, errors);
    

  cosm->SetLineColor(kGreen + 2);
  cosm->SetLineWidth(linewidthsub);
  cosm->Scale(cy/cosm->Integral());

    
  radar->SetLineColor(kMagenta + 1);
  radar->SetLineWidth(linewidthsub);
  radar->Scale(ary/radar->Integral());

        
  radk->SetLineColor(kCyan + 2);
  radk->SetLineWidth(linewidthsub);
  radk->Scale(ky/radk->Integral());

    
  radkr->SetLineColor(kAzure - 1);
  radkr->SetLineWidth(linewidthsub);
  radkr->Scale(kry/radkr->Integral());

    
  radu->SetLineColor(kOrange + 3);
  radu->SetLineWidth(linewidthsub);
  radu->Scale(uy/radu->Integral());

        
  radth->SetLineColor(kOrange + 10);
  radth->SetLineWidth(linewidthsub);
  radth->Scale(thy/radth->Integral());


  radco->SetLineColor(kBlue + 4);
  radco->SetLineWidth(linewidthsub);
  radco->Scale(coy/radco->Integral());


  // radrn->SetLineColor(kViolet + 4);
  // radrn->SetLineWidth(linewidthsub);
  // radrn->Scale(rny/radrn->Integral());


  // mod = (TH1D*) f1->GetPlot();
  // (*mod) = (*cosm) + (*radar) + (*radu) + (*radk) + (*radkr) + (*radth) + (*radco) + (*radrn);
  (*mod) = (*cosm) + (*radar) + (*radu) + (*radk) + (*radkr) + (*radth) + (*radco);
  //   //(*mod) = (*cosm) + (*radar) + (*radk) + (*radkr);
  mod->SetLineColor(kRed);    
  mod->SetLineWidth(linewidth);

  //   reg->SetTitle("Data");
  //   mod->SetTitle("MC Fit");

  //   cosm->SetTitle("Cosmic");
  //   radar->SetTitle("39Ar");
  //   radu->SetTitle("238U");
  //   radk->SetTitle("40K");
  //   radkr->SetTitle("85Kr");
  //   radth->SetTitle("232Th");
  //   radth->SetTitle("60Co");
  //   radth->SetTitle("222Rn");

  //   //c1->BuildLegend();

  if (converged){
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->SetShadowColor(0);
    leg->SetTextFont(62);
    leg->SetTextSize(0.035);
    // leg->AddEntry(sub,"Data (8 ADC Threshold)","L");
    leg->AddEntry(reg,"Data (Run 5841)","L");
    // leg->AddEntry(reg,"Data (8 ADC Threshold)","L");
    leg->AddEntry(mod,"MC Fit","L");
    leg->AddEntry((TObject*)NULL,"- - - - - - - - - - - - - - - - - - - -","");
    leg->AddEntry(cosm,TString::Format("Cosmic (%.2f%)", (cosm->Integral()/mod->Integral())*100),"L");
    leg->AddEntry(radar,TString::Format("^{39}Ar (%.2f%)", (radar->Integral()/mod->Integral())*100),"L");
    leg->AddEntry(radkr,TString::Format("^{85}Kr (%.2f%)", (radkr->Integral()/mod->Integral())*100),"L");
    leg->AddEntry(radco,TString::Format("^{60}Co (%.2f%)", (radco->Integral()/mod->Integral())*100),"L");
    leg->AddEntry(radk,TString::Format("^{40}K (%.2f%)", (radk->Integral()/mod->Integral())*100),"L");
    leg->AddEntry(radu,TString::Format("^{238}U (%.2f%)", (radu->Integral()/mod->Integral())*100),"L");
    leg->AddEntry(radth,TString::Format("^{232}Th (%.2f%)", (radth->Integral()/mod->Integral())*100),"L");
    // leg->AddEntry(radrn,TString::Format("^{222}Rn (%.2f%)", (radrn->Integral()/mod->Integral())*100),"L");


    // leg->AddEntry(reg,"Data (Run 5841)","L");
    // //leg->AddEntry(reg,"Data (8 ADC Threshold)","L");
    // // leg->AddEntry(mod,"MC Fit","L");
    // // leg->AddEntry((TObject*)NULL,"- - - - - - - - - - - - - - - - - - - -","");
    // leg->AddEntry(cosm,"Simulated Cosmic","L");
    // leg->AddEntry(radar,"Simulated ^{39}Ar","L");
    // leg->AddEntry(radkr,"Simulated ^{85}Kr","L");
    // leg->AddEntry(radco,"Simulated ^{60}Co","L");
    // leg->AddEntry(radk,"Simulated ^{40}K","L");
    // // leg->AddEntry(radu,"Simulated ^{238}U","L");
    // // leg->AddEntry(radth,"Simulated ^{232}Th","L");
    // // leg->AddEntry(radrn,"Simulated ^{222}Rn","L");

    // //TTree reader stuff goes here
    

    reg->GetYaxis()->SetRangeUser(GetMax(reg) * 1E-4, GetMax(reg) * 1.15);
    

    //fout.Write();
    
    reg->Draw("HIST");
    mod->Draw("HISTSAME");
    // //sub->Draw("HISTSAME");
    cosm->Draw("HISTSAME");
    radar->Draw("HISTSAME");
    radk->Draw("HISTSAME");
    radkr->Draw("HISTSAME");
    radu->Draw("HISTSAME");
    radth->Draw("HISTSAME");
    radco->Draw("HISTSAME");
    // radrn->Draw("HISTSAME");
    leg->Draw("SAME");
    
    /*
      reg->Draw("");
      mod->Draw("SAME");
      sub->Draw("SAME");
      cosm->Draw("SAME");
      radar->Draw("SAME");
      radk->Draw("SAME");
      radkr->Draw("SAME");
      radu->Draw("SAME");
      radth->Draw("SAME");
      leg->Draw("SAME");
    */
    //ct.DrawLatex(1100,7300.0,TString::Format("#chi_{r}^{2} =  %.2f", f1->GetChisquare()/f1->GetNDF()));
    ct.SetNDC();
    //ct.DrawLatex(0.56,0.83,TString::Format("#chi_{r}^{2} =  %.2f", f1->GetChisquare()));
    ct.DrawLatex(0.56,0.83,TString::Format("#chi_{r}^{2} =  %.2f", chisquare/f1->GetNDF()));
    ct.DrawLatex(0.56,0.34,TString::Format("Smearing %.2f %%", arx * 100.0));
    // }
  }else{
    cout << "///////////////////////////////////////////" << endl;
    cout << "Nothing converged properly; nothing to draw" << endl;
    cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
  }
  //std::cout << std::endl << arx << ", " << ary << ", " << cy << ", " << ky << ", " << kry << ", " << uy <<std::endl;
  //ct.DrawLatex(2000,12000,TString::Format("#chi_{r}^{2} =  %.2f", chitemp/95.0));

}

vector<double> vGain(const string& filename){
  vector<double> gain(15360, 0.0234);
  unsigned int tmp_chn;

  if (filename == "null"){return gain;}


  ifstream gain_file(filename);

  for (unsigned int i = 0; i < gain.size(); i++){
    gain_file >> tmp_chn >> gain[i];
  }

  gain_file.close();


  return gain;
}

void WeierstrassTransform(TH1D*& tmp, TH1D*& h2, double s, unsigned short offset){
  double val = 0;
  short nbins = h2->GetXaxis()->GetNbins();
  double xmin = h2->GetXaxis()->GetXmin();
  double xmax = h2->GetXaxis()->GetXmax();
  double binwidth = (xmax - xmin) / nbins;

  TH1D *h1 = new TH1D("h1", "", nbins + (2*offset), xmin - (offset * binwidth), xmax + (offset * binwidth));
  TH1D *gh = new TH1D("gh", "", nbins + (2*offset), xmin - (offset * binwidth), xmax + (offset * binwidth));

  TF1 *gf = new TF1("gf", "gausn", gh->GetXaxis()->GetXmin(), gh->GetXaxis()->GetXmax());

  if (s != 0.0){
    for (unsigned short i = 1; i < nbins; i++){
      h1->Fill(tmp->GetXaxis()->GetBinCenter(i), tmp->GetBinContent(i));
    }

    for (unsigned short j = 1; j <= h1->GetXaxis()->GetNbins(); j++){
      gh->Reset();
      if (gh->GetXaxis()->GetBinCenter(j) >= xmin && gh->GetXaxis()->GetBinCenter(j) <= xmax){
	gf->SetParameters(1, 0, s * h1->GetXaxis()->GetBinCenter(j));
	for (unsigned short i = 1; i <= gh->GetXaxis()->GetNbins(); i++){
	  gh->Fill(gh->GetXaxis()->GetBinCenter(i), gf->Eval(gh->GetXaxis()->GetBinCenter(i)));
	}
      }else{
	gh->SetBinContent(1, 1);
      }

      for (unsigned short k = 1; k <= h1->GetXaxis()->GetNbins(); k++){
	val += h1->GetBinContent(j - k) * gh->GetBinContent(k);
      }
    
      h2->SetBinContent(j - (2 * offset), val);

      val = 0;
    }
  }

  delete gf;
  delete h1;
  delete gh;
}

void RecombinationSmear(TH1D*& tmp, TH1D*& h2, double s, unsigned short offset){
  double val = 0;
  short nbins = h2->GetXaxis()->GetNbins();
  double xmin = h2->GetXaxis()->GetXmin();
  double xmax = h2->GetXaxis()->GetXmax();
  double binwidth = (xmax - xmin) / nbins;

  TH1D *h1 = new TH1D("h1", "", nbins + (2*offset), xmin - (offset * binwidth), xmax + (offset * binwidth));
  TH1D *gh = new TH1D("gh", "", nbins + (2*offset), xmin - (offset * binwidth), xmax + (offset * binwidth));

  TF1 *gf = new TF1("gf", "gausn", gh->GetXaxis()->GetXmin(), gh->GetXaxis()->GetXmax());
  TF1 *sq = new TF1("sq", "1/sqrt(x)", gh->GetXaxis()->GetXmin(), gh->GetXaxis()->GetXmax());

  for (unsigned short i = 1; i < nbins; i++){
    h1->Fill(tmp->GetXaxis()->GetBinCenter(i), tmp->GetBinContent(i));
  }

  for (unsigned short j = 1; j <= h1->GetXaxis()->GetNbins(); j++){
    gh->Reset();
    if (gh->GetXaxis()->GetBinCenter(j) >= xmin && gh->GetXaxis()->GetBinCenter(j) <= xmax){
      gf->SetParameters(1, 0, s * sq->Eval(h1->GetXaxis()->GetBinCenter(j)));
      for (unsigned short i = 1; i <= gh->GetXaxis()->GetNbins(); i++){
	gh->Fill(gh->GetXaxis()->GetBinCenter(i), gf->Eval(gh->GetXaxis()->GetBinCenter(i)));
      }
    }else{
      gh->SetBinContent(1, 1);
    }

    for (unsigned short k = 1; k <= h1->GetXaxis()->GetNbins(); k++){
      val += h1->GetBinContent(j - k) * gh->GetBinContent(k);
    }
    
    h2->SetBinContent(j - (2 * offset), val);

    val = 0;
  }

  delete sq;
  delete gf;
  delete h1;
  delete gh;
}
