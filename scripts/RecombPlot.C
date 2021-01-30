/////////////////////////////////////////////////////
// CERN ROOT Script to plot the energy dependent   //
// recombination factor according to mod box model //
// This script was written by Mike Mooney          //
///////////////////////////////////////////////////// 

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

double getRecombVal(double E_start, double step_size = 0.001);

void RecombPlot() {

  double n = 1001;
  double x[1001];
  double x2[1001];
  double y[1001];

  for(int i = 0; i < n; i++) {
    x[i] = 4.0*i/(n-1.0);
    y[i] = getRecombVal(x[i],0.000001);
    x[i] *= 1000.0;
    x2[i] = x[i]*y[i];
    cout << x[i] << " " << y[i] << endl;
  }

  TGraph *recombGraph_vsEtot = new TGraph(n,x,y);
  recombGraph_vsEtot->SetName("recombGraph_vsEtot");
  recombGraph_vsEtot->SetTitle("");

  TGraph *recombGraph_vsEfree = new TGraph(n,x2,y);
  recombGraph_vsEfree->SetName("recombGraph_vsEfree");
  recombGraph_vsEfree->SetTitle("");

  TFile *outFile = new TFile("recombfile.root","RECREATE");
  outFile->cd();
  recombGraph_vsEtot->Write();
  recombGraph_vsEfree->Write();
  outFile->Close();

  TCanvas *recombCanv = new TCanvas();
  recombCanv->cd();
  recombGraph_vsEtot->SetLineWidth(3.0);
  recombGraph_vsEtot->SetLineColor(kBlue);
  recombGraph_vsEtot->Draw("");
  recombGraph_vsEtot->SetTitle("");
  recombGraph_vsEtot->GetXaxis()->SetTitle("Initial Electron Kinetic Energy [keV]");
  recombGraph_vsEtot->GetYaxis()->SetTitle("Effective Recombination Value");
  recombGraph_vsEtot->GetXaxis()->SetTitleSize(0.055);
  recombGraph_vsEtot->GetXaxis()->SetTitleOffset(0.8);
  recombGraph_vsEtot->GetYaxis()->SetTitleSize(0.055);
  recombGraph_vsEtot->GetYaxis()->SetTitleOffset(0.8);
  recombGraph_vsEtot->GetXaxis()->SetRangeUser(0.0,700.0);
  recombGraph_vsEtot->GetYaxis()->SetRangeUser(0.0,0.65);
  recombCanv->SaveAs("recombPlot_vsEtot.png");

  TCanvas *recombCanv2 = new TCanvas();
  recombCanv2->cd();
  recombGraph_vsEfree->SetLineWidth(3.0);
  recombGraph_vsEfree->SetLineColor(kBlue);
  recombGraph_vsEfree->Draw("");
  recombGraph_vsEfree->SetTitle("");
  recombGraph_vsEfree->GetXaxis()->SetTitle("Initial Electron Free Energy [keV]");
  recombGraph_vsEfree->GetYaxis()->SetTitle("Effective Recombination Value");
  recombGraph_vsEfree->GetXaxis()->SetTitleSize(0.055);
  recombGraph_vsEfree->GetXaxis()->SetTitleOffset(0.8);
  recombGraph_vsEfree->GetYaxis()->SetTitleSize(0.055);
  recombGraph_vsEfree->GetYaxis()->SetTitleOffset(0.8);
  recombGraph_vsEfree->GetXaxis()->SetRangeUser(0.0,400.0);
  recombGraph_vsEfree->GetYaxis()->SetRangeUser(0.0,0.65);
  recombCanv2->SaveAs("recombPlot_vsEfree.png");
}

double getRecombVal(double E_start, double step_size = 0.001) {

  // const double Edrift = 0.2739; // MicroBooNE E field value [kv/cm]
  // const double Edrift = 0.5; // DUNE exp E filed value [kv/cm]
  const double Edrift = 0.4867; // ProtoDUNE E field value [kV/cm]
  const double ModBoxA = 0.930; // recomb constant A used in LArG4
  const double ModBoxB = 0.212; // recomb constant B used in LArG4
  const double piVal = 3.14159265;
  const double elec_mass = 0.511;
  const double LAr_density = 1.3973;
  
  double T = E_start;
  double E_dep = 0.0;
  while(T > 0.0) {
    double gamma = (elec_mass+T)/elec_mass;
    double beta = pow(1.0-pow(gamma,-2.0),0.5);
    double dEdx = (2.0*piVal*pow(197.0,2.0)*pow(137.0,-2.0)/(elec_mass*pow(beta,2.0)))*((LAr_density*6.022*pow(10.0,-3.0))/40.0)*18.0*(log((elec_mass*pow(beta,2.0)*T)/(2.0*pow(0.000188,2.0)*(1-pow(beta,2.0))))-log(2.0)*(2.0*pow(1.0-pow(beta,2.0),0.5)-1.0+pow(beta,2.0))+1.0-pow(beta,2.0)+pow(1-pow(1.0-pow(beta,2.0),0.5),2.0)/8.0);

    if ((T < 0.0005) || (dEdx < 0.0)) break;
    
    double Xi = (ModBoxB * dEdx) / (Edrift * LAr_density);
    double recomb = log(ModBoxA + Xi) / Xi;

    E_dep += recomb*dEdx*step_size;

    //cout << T << " " << dEdx << " " << recomb << " " << E_dep << endl;

    T -= dEdx*step_size;
  }

  if (E_start <= 0.0005) {
    return 0.0;
  }
  else {
    return E_dep/E_start;
  }
}
