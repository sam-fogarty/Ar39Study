#ifndef POINTSIGNALFUNCTIONS_H
#define POINTSIGNALFUNCTIONS_H

#include <vector>

#include "TInterpreter.h"
//#include "TString.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"

#include "gallery/Event.h"

#include "Waveforms.h"

// Calculate line from two 2D points
double bline(double x1, double y1, double x2, double y2, double x);

namespace rad_analysis{

  // Select point-like signal and record their waveforms
  std::vector< std::vector<char> > Signal_Select(rad_analysis::Waveforms& nfrw,
                                                std::vector<float>& ICharge,
                                                std::vector<double>& Par,
                                                std::vector<int>& IntWindow,
                                                std::vector<double>& c_info,
                                                std::vector<float>& TEinfo,
                                                TTree& ctree,
                                                int& cind,
                                                const gallery::Event& ev);

  // Confirm candidates based on 3D info on multiple planes; needs improvement; MicroBooNE only
  void Confirm_Candidates(std::vector<short>& fcoor, rad_analysis::Waveforms& nfrw);

}
#endif
