//////////////////////////////////////////
// 39Ar analysis gallery project        //
// Alex Flesher and Mike Mooney         //
// Presently only works in MicroBooNE   //
// and ProtoDUNE.                       //
//////////////////////////////////////////

// Include headers for this project
#include "../inc/Waveforms.h"
#include "../inc/PointSignalFunctions.h"

//some standard C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <numeric>

//some ROOT includes
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

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/Handle.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
//#include "larevt/Filters/ChannelFilter.h"


#define PI 3.141592


using namespace art;
using namespace std;
using namespace std::chrono;


int compare(const void * a, const void * b){  
  return ( *(double*)a - *(double*)b );
}

bool vcompare(vector<short>& a, vector<short>& b){ 
  return (a[1] < b[1]);
}


int main(int argc, char* argv[]) {

  //Begin setting all parameters

  //const int maxEvents = 1000000;
  int maxEvents = 250;
  int maxEventsPerRun = 1000000;
  
  bool nfilter;


  string detector;
  string intag;


  rad_analysis::Waveforms nfrw;


  vector<double> Par(12);

  ifstream pfile("param.txt");

  pfile >> Par.at(0) >> Par.at(1) >> Par.at(2) >> Par.at(3) >> Par.at(4)
	>> Par.at(5) >> Par.at(6) >> Par.at(7) >> Par.at(8) >> Par.at(9)
	>> Par.at(10) >> maxEvents >> maxEventsPerRun >> Par.at(11) 
	>> nfrw.data_properties.MCCX >> detector >> intag >> nfrw.data_properties.SCALE; 
  
  pfile.close();


  InputTag product_tag  { intag };

  
  cout << "Detector: " << detector << endl;

  //Set parameters for the detector
  if (detector == "MicroBooNE"){
    nfrw.detector_properties.DET = 0;    //microboone is detector 0 here
    nfrw.detector_properties.CHN = 8256; //it has 8256 channels
    nfrw.detector_properties.NCW = 3456; //3456 of them are collection
    nfrw.detector_properties.NTT = 6400; //look at 6400 time ticks
    nfrw.detector_properties.ntpcs = 1;  //only 1 tpc
    
    //set noise filter option and indicate if reco
    if (intag == "nfspl1:raw" || intag == "butcher" || intag == "wcNoiseFilter" || intag == "driftWC:orig"){
      nfilter = 0;
      cout << intag << " selected. Noise filter off." << endl;      
    }

    else if (intag == "nfspl1:wiener" || intag == "nfspl1:gauss" || intag == "caldata" ||
	     intag == "simnfspl1:wiener" || intag == "simnfspl1:gauss"){
      nfilter = 0;
      //Par[3] /= 6.05724; //Improve: determine this dynamically
      nfrw.data_properties.recoused = true;

      if (intag == "simnfspl1:wiener" || intag == "simnfspl1:gauss"){
	nfrw.data_properties.simulated = true;
      }

      cout << intag << "selected. Noise filter off. recob used." << endl;
    }

    else if (intag == "daq"){
      nfilter = 1;
      cout << intag << " selected. Noise filter on." << endl;
    }
  }

  else if (detector == "ProtoDUNE"){
    nfrw.detector_properties.DET = 1;
    nfrw.detector_properties.CHN = 15360;
    nfrw.detector_properties.NCW = 5760;
    nfrw.detector_properties.NTT = 6000;
    nfrw.detector_properties.ntpcs = 12;
    nfilter = 1;
  }

  vector<string> filenames;


  nfrw.Fill_Wire_Maps();


  ifstream myfile;


  //Ar39Study has different command syntax
  //"Ar39Study " simply runs Ar39Study with the parameters set in param.txt
  //"Ar39Study file.txt" overrides the default filelist.txt
  //"Ar39Study num1 num2" will override the event numbers and only look between event num1 and num2
  //"Ar39Study 1 2" will only look at event 1 and event 2

  if (argc == 1 || argc == 3){
    myfile.open ("filelist.txt");
    copy(istream_iterator<string>(myfile),istream_iterator<string>(),back_inserter(filenames));

    if (argc == 3){ //can tell range of events from command line
      maxEvents = boost::lexical_cast<int>(argv[2]);
    }
  }else if (argc == 2){ //tell which file to run, not from filelist (grid running)
    filenames.push_back(argv[1]);
  }


  TFile f_output("output.root","RECREATE");
  TTree ctree("ctree", "Candidate info");
  TTree etree("etree", "Event info");


  gStyle->SetNumberContours(255);
  gStyle->SetOptStat(0);


  TH1F *wire8049 = new TH1F("wire8049", "", nfrw.detector_properties.NTT, 0, nfrw.detector_properties.NTT);

  TH1F *ModeHist = new TH1F("ModeHist", "ADC baseline", nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN);
  TH1F *SDHist = new TH1F("SDHist","Standard Deviation per Wire", nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN);  
  //TH1F *SDHist = new TH1F("SDHist","Standard Deviation per Wire", nfrw.detector_properties.NCW, 0, nfrw.detector_properties.NCW);  

  TH1F *ChargeHist = new TH1F("ChargeHist","", 10000, -5000.5, 4999.5);
  TH1F *EnergyHist = new TH1F("EnergyHist","Ar39 Decay Energy", 10000, -5000.5 * 1.556, 4999.5 * 1.556);
  //TH1F *EnergyHist2 = new TH1F("EnergyHist2","Ar39 Decay Energy", 10000, -5000.5 * 1.556, 4999.5 * 1.556);

  TH1I *DetPerEvent = new TH1I("DetPerEvent","Detections of Ar39 Per Event", 1000, -0.5, 999.5);

  TH2F *DetectionMap = new TH2F("DetectionMap", "Map of Points over Threshold", nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN, nfrw.detector_properties.NTT, 0, nfrw.detector_properties.NTT); 
  TH2F *ChargeMap = new TH2F("ChargeMap", "", nfrw.detector_properties.NCW, 0, nfrw.detector_properties.NCW, nfrw.detector_properties.NTT, 0, nfrw.detector_properties.NTT); 
  TH2F *ChargeMap2 = new TH2F("ChargeMap2", "", nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN, nfrw.detector_properties.NTT, 0, nfrw.detector_properties.NTT); 
  TH2F *ZDetectionMap = new TH2F("ZDetectionMap", "Zoomed In Detection Map", 300, 500, 800, 1000, 0, 1000);

  TH2F *EvViewI = new TH2F("EvViewI", "", ((double)nfrw.detector_properties.NCW/nfrw.detector_properties.ntpcs)*3, 0, ((double)nfrw.detector_properties.NCW/nfrw.detector_properties.ntpcs)*3, nfrw.detector_properties.NTT*2, -nfrw.detector_properties.NTT, nfrw.detector_properties.NTT); 
  TH1F *insd = new TH1F("insd","", nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN);  

  TH2F *EvViewO = new TH2F("EvViewO", "", ((double)nfrw.detector_properties.NCW/nfrw.detector_properties.ntpcs)*3, 0, ((double)nfrw.detector_properties.NCW/nfrw.detector_properties.ntpcs)*3, nfrw.detector_properties.NTT*2, -nfrw.detector_properties.NTT, nfrw.detector_properties.NTT); 
  TH1F *osd = new TH1F("osd","", nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN);  

  //variables to control file list loop
  int totalEventNum = 0;
  int writecount = 0;
  int prevRun = -99999;
  int runEventCount = 0;

  int thisRunNum;
  int thisSubrunNum;
  int thisEventNum;

  int cind;

  //Vectors that hold properties of the event and candidates
  vector<float> TEinfo(3);

  vector<int> IntWindow;

  vector<float> ICharge;
  vector< vector<char> > temp;
  //vector< vector<float> > temp;

  vector<double> c_info(15);

  //vector<bool> deadChannel;
  //vector<bool> tExclude;
  //vector<float> wireMd;

  //vector<int> modeVec;


  //define TTree branches with these vectors
  etree.Branch("tep", &TEinfo[0], "tep/F");
  etree.Branch("total", &TEinfo[1], "total/F");
  etree.Branch("remfrac", &TEinfo[2], "remfrac/F");
  etree.Branch("evIndex", &totalEventNum, "evIndex/I");
  etree.Branch("evNum", &thisEventNum, "evNum/I");
  etree.Branch("runNum", &thisRunNum, "runNum/I");
  etree.Branch("subrun", &thisSubrunNum, "subrun/I");

  ctree.Branch("ChargeVec", &IntWindow);
  ctree.Branch("charge", &c_info[0], "charge/D");
  ctree.Branch("energy", &c_info[1], "energy/D"); 
  ctree.Branch("maxAmp", &c_info[2], "maxAmp/D");
  ctree.Branch("wireWidth", &c_info[3], "wireWidth/D");
  ctree.Branch("timeLength", &c_info[4], "timeLength/D");
  ctree.Branch("channel", &c_info[5], "channel/D");
  ctree.Branch("time", &c_info[6], "time/D");
  ctree.Branch("tOverThresh", &c_info[7], "tOverThresh/D");
  ctree.Branch("nearTrack", &c_info[8], "nearTrack/D");
  ctree.Branch("evNum", &thisEventNum, "evNum/I");
  ctree.Branch("runNum", &thisRunNum, "runNum/I");
  ctree.Branch("subrun", &thisSubrunNum, "subrun/I");
  ctree.Branch("evIndex", &totalEventNum, "evIndex/I");
  ctree.Branch("xWind", &Par[5], "xWind/D");
  ctree.Branch("yWind", &Par[6], "yWind/D");
  ctree.Branch("xCheck", &Par[7], "xCheck/D");
  ctree.Branch("yCheck", &Par[8], "yCheck/D");
  ctree.Branch("xCheckWidth", &Par[9], "xCheckWidth/D");
  ctree.Branch("yCheckWidth", &Par[10], "yCheckWidth/D");
  ctree.Branch("avglead", &c_info[9], "avglead/D");
  ctree.Branch("avgtail", &c_info[10], "avgtail/D");
  ctree.Branch("tpc", &c_info[11], "tpc/D");
  ctree.Branch("candidateIndex", &cind, "candidateIndex/I");
  ctree.Branch("confirmed", &c_info[12], "confirmed/D");
  ctree.Branch("y", &c_info[13], "y/D");
  ctree.Branch("z", &c_info[14], "z/D");

  EnergyHist->GetXaxis()->SetTitle("keV");
  EnergyHist->GetYaxis()->SetTitle("Number of Detections");
  EnergyHist->GetYaxis()->SetTitleOffset(1.2);

  srand(time(NULL));


  //cout << "Seems to be working..." << endl; //debug all above this

  //start looping through data files and events
  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
    auto t_begin = high_resolution_clock::now();

    cind = 0;
    double dcounter = 0;    //counts dead channels (also gets divided into percent of detector area)
    int counter = 0;        //used in the track exclusion algorithm


    thisRunNum = ev.eventAuxiliary().run();
    thisSubrunNum = ev.eventAuxiliary().subRun();
    thisEventNum = ev.eventAuxiliary().event();    

    //deadChannel = vector<bool>(nfrw.detector_properties.CHN, 0);
    //tExclude = vector<bool>(nfrw.detector_properties.CHN, 0);

    //wireMd = vector<float>(nfrw.detector_properties.CHN, 0);

    TEinfo[0] = 0.0;
    TEinfo[1] = 0.0;
    TEinfo[2] = 0.0;


    //control the sequence/number of events
    if ((thisRunNum != prevRun) && (runEventCount > 0)){
      prevRun = thisRunNum;
      runEventCount = 1;
    }
    else if ((thisRunNum == prevRun) && (runEventCount >= maxEventsPerRun)){
      continue;
    }
    else {
      prevRun = thisRunNum;
      runEventCount++;
    }

    totalEventNum++;
    writecount++;

    if (totalEventNum > maxEvents) {
      cout << "Event number reached. Closing." << endl;  
      
      f_output.Write();
      f_output.Close();
      myfile.close();
      
      return 0;
    }

    //Keeps progress through long runs if program crashes
    if (writecount == 100){
      f_output.Write();
      writecount = 0;
    }

    //skip to this event
    if (argc > 2 && totalEventNum < boost::lexical_cast<int>(argv[1]))
      continue;

    //debug
    //if (totalEventNum < 400)
    //continue;

    cout << "Processing "
  	 << "Run " << thisRunNum << ", "
  	 << "Subrun " << thisSubrunNum << ", "
  	 << "Event " << thisEventNum << endl;



    ////Declare the raw and reco vectors to fill.
    vector<recob::Wire> wire_vec;
    vector<recob::Wire> wire_vec2;
    vector<raw::RawDigit> allrawdigits_vec;
    vector<int> badchan;

    if (!nfrw.data_properties.recoused){ //rawdigits
      allrawdigits_vec = vector<raw::RawDigit>(*ev.getValidHandle< vector<raw::RawDigit> >(product_tag));
    }
    
    else if (nfrw.data_properties.recoused){ //reco (PDUNE-SP Unsupported)
      if (nfrw.data_properties.MCCX == 8){
    	//keep raw digits for dead channel search until I can figure out how channelservice works.
    	//Also need raw digits for dynamic scale

    	allrawdigits_vec = vector<raw::RawDigit>(*ev.getValidHandle< vector<raw::RawDigit> >("wcNoiseFilter"));
	wire_vec = vector<recob::Wire>(*ev.getValidHandle< vector<recob::Wire> >(product_tag));
      }

      else if (nfrw.data_properties.MCCX == 9){
	if (!nfrw.data_properties.simulated){
	  wire_vec2 = vector<recob::Wire>(*ev.getValidHandle< vector<recob::Wire> >("nfspl1:gauss"));
	  wire_vec = vector<recob::Wire>(*ev.getValidHandle< vector<recob::Wire> >(product_tag));
	  //allrawdigits_vec = vector<raw::RawDigit>(*ev.getValidHandle< vector<raw::RawDigit> >("daq"));
	}else{
	  wire_vec2 = vector<recob::Wire>(*ev.getValidHandle< vector<recob::Wire> >("simnfspl1:gauss"));
	  wire_vec = vector<recob::Wire>(*ev.getValidHandle< vector<recob::Wire> >(product_tag));
	}
      }

      else {
    	cout << "Param nfrw.data_properties.MCCX not set correctly." << endl;
    	return 0;
      }
    }



    DetectionMap->Reset();
    ZDetectionMap->Reset();
    ModeHist->Reset();
    SDHist->Reset();

    
    nfrw.initialize_data (ev, allrawdigits_vec, wire_vec, wire_vec2, nfilter, Par);

    if (nfrw.data_properties.recoused){
      // allrawdigits_vec.clear();
      vector<raw::RawDigit>().swap(allrawdigits_vec); //Clear and free memory from ard_vec

      //nfrw.NoRawDC(badchan); //assigns badchannels to the format used here

      //badchan.clear();
      vector<int>().swap(badchan); //clear badchan after use
    }

    Par[3] /= nfrw.f_reco_scale; //changes the scale of the threshold from adcs to the reconstructed units
 

    // // if (!nfrw.data_properties.recoused){	     
    // for (int t = 0; t < nfrw.detector_properties.NTT; t++){
    //   wire8049->Fill(t, nfrw.ADCvalVec.at(8049).at(t));
    //   // wire8049->Fill(t, allrawdigits_vec.at(7150).ADC(t));
    // }
    // // }

    
    for (int i = 0; i < nfrw.detector_properties.CHN; i++){
      ModeHist->SetBinContent(i, nfrw.baseline[i]);
      SDHist->SetBinContent(i, nfrw.sd_vector[i]);
    }

    
    cout << "try rad_analysis::Signal_Select()" << endl;
    temp = rad_analysis::Signal_Select(nfrw, ICharge, Par, IntWindow, c_info, TEinfo, ctree, cind, ev);


    cout << "rad_analysis::Signal_Select() ended; fill the things with the stuff." << endl;

    // Debug charge maps for looking at events
    // for (unsigned short p = 0; p < nfrw.collection_channel.size(); p++){
    //   for (unsigned short s = 0; s < nfrw.collection_channel[p].size(); s++){
    // 	for (unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
    // 	  ChargeMap->Fill((p * nfrw.collection_channel[p].size()) + s, t, nfrw.ADCvalVec[nfrw.collection_channel[p][s]][t]);
    // 	}
    //   }
    // }

    // for (unsigned short s = 0; s < nfrw.detector_properties.CHN; s++){
    //   for (unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
    // 	ChargeMap2->Fill(s, t, nfrw.ADCvalVec[s][t]);
    //   }
    // }
    

    TEinfo[0] = TEinfo[0] / (6400.0*3456.0);
    etree.Fill();


    //if instructed, draw the event view 
    if (Par[0] == 1 && nfrw.detector_properties.DET == 0){

      //// For induction planes and collection plane
      for (unsigned short p = 0; p < nfrw.collection_channel.size(); p++){
	for (unsigned short s = 0; s < nfrw.collection_channel[p].size(); s++){
	  for(unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
	    DetectionMap->Fill(nfrw.collection_channel[p][s], t, static_cast<int>(temp.at(nfrw.collection_channel[p][s]).at(t)));
	  }
	}
      }      

      for (unsigned short p = 0; p < nfrw.induction_channel_1.size(); p++){
	for (unsigned short s = 0; s < nfrw.induction_channel_1[p].size(); s++){
	  for(unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
	    DetectionMap->Fill(nfrw.induction_channel_1[p][s], t, static_cast<int>(temp.at(nfrw.induction_channel_1[p][s]).at(t)));
	  }
	}
      }      

      for (unsigned short p = 0; p < nfrw.induction_channel_2.size(); p++){
	for (unsigned short s = 0; s < nfrw.induction_channel_2[p].size(); s++){
	  for(unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
	    DetectionMap->Fill(nfrw.induction_channel_2[p][s], t, static_cast<int>(temp.at(nfrw.induction_channel_2[p][s]).at(t)));
	  }
	}
      }  
    }

    
    //only for pdune right now
    else if (Par[0] == 1 && nfrw.detector_properties.DET == 1){
      //Doesn't work for Induction Planes
      for (unsigned short p = 0; p < nfrw.collection_channel.size(); p++){
      	for (unsigned short s = 0, f = 0; s < nfrw.collection_channel[p].size(); s++){
      	  for(unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
      	    f = (p) * nfrw.detector_properties.NCW/12 + s;

      	    DetectionMap->Fill(f, t, static_cast<int>(temp.at(nfrw.collection_channel[p][s]).at(t)));
	    

    	    ////inner planes
    	    if (p == 1)
    	      EvViewI->Fill(s, -t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    else if (p == 5)
    	      EvViewI->Fill(480 + s, -t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    else if (p == 9)
    	      EvViewI->Fill(480*2 + s, -t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    else if (p == 2)
    	      EvViewI->Fill(s, t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    else if (p == 6)
    	      EvViewI->Fill(480 + s, t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    else if (p == 10)
    	      EvViewI->Fill(480*2 + s, t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));


    	    ////outer planes
    	    else if (p == 0)
    	      EvViewO->Fill(s, -t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    else if (p == 4)
    	      EvViewO->Fill(480 + s, -t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    else if (p == 8)
    	      EvViewO->Fill(480*2 + s, -t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    else if (p == 3)
    	      EvViewO->Fill(s, t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    else if (p == 7)
    	      EvViewO->Fill(480 + s, t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    else if (p == 11)
    	      EvViewO->Fill(480*2 + s, t, DetectionMap->GetBinContent(nfrw.collection_channel[p].size()*(p) + s, t));

    	    // if (p == 1 || p == 5 || p == 9){
    	    //   EvViewI->Fill(floor(log2(p)) * nfrw.collection_channel[p].size() + s, -t, static_cast<int>(temp.at(nfrw.collection_channel[p][s]).at(t)));
    	    // }
    	    // else if (p == 2 || p == 6 || p == 10){
    	    //   EvViewI->Fill(floor(log2(p - 1)) * nfrw.collection_channel[p].size() + s, t, static_cast<int>(temp.at(nfrw.collection_channel[p][s]).at(t)));
    	    // }
	    	  
      	  }
      	} 
      }

      // for (unsigned short p = 0; p < 3; p++){
      // 	for (unsigned short s = 0, f = 0; s < nfrw.collection_channel[p].size(); s++){
      // 	  for(unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
      // 	    f = (p + 1) * nfrw.detector_properties.NCW/12 + s;

      // 	    DetectionMap->Fill(f, t, static_cast<int>(temp.at(nfrw.collection_channel[p][s]).at(t)));
	  
      // 	  }
      // 	} 
      // }

      // //inner planes
      // for (unsigned short p = 1, f = 1; (unsigned)(p + 3) < nfrw.collection_channel.size(); p += 4, f++){
      // 	for (unsigned short s = 0, c = 0; s < nfrw.collection_channel[p].size(); s++){
      // 	  for(short t = 0; t < nfrw.detector_properties.NTT; t++){
      // 	    c = f * nfrw.collection_channel[p].size() + s;

      // 	    EvViewI->Fill(c, -t, static_cast<int>(temp.at(nfrw.collection_channel[p][s]).at(t)));
      // 	    EvViewI->Fill(c, t, static_cast<int>(temp.at(nfrw.collection_channel[p + 1][s]).at(t)));
	  
      // 	  }
      // 	} 
      // }

      // //outer planes
      // for (unsigned short p = 0, f = 1; (unsigned)(p + 4) < nfrw.collection_channel.size(); p += 4, f++){
      // 	for (unsigned short s = 0, c = 0; s < nfrw.collection_channel[p].size(); s++){
      // 	  for(short t = -nfrw.detector_properties.NTT + 1; t < nfrw.detector_properties.NTT; t++){
      // 	    c = f * nfrw.collection_channel[p].size() + s;

      // 	    if (t < 0)
      // 	      EvViewO->Fill(c, t, static_cast<int>(temp.at(nfrw.collection_channel[p][s]).at(-t)));

      // 	    else
      // 	      EvViewO->Fill(c, t, static_cast<int>(temp.at(nfrw.collection_channel[p + 3][s]).at(t)));
	  
      // 	  }
      // 	} 
      // }

      // for (int s = 5300; s < 5600; s++){
      // 	for (int t = 0; t < 1000; t++){
      // 	  ZDetectionMap->Fill(s - 4800, t, static_cast<int>(temp.at(s).at(t)));
      // 	}
      // }
    }


    //asemble data acquired from rad_analysis::Signal_Select()
    if (Par[2] == 2){
      
      sort (ICharge.begin(), ICharge.end());

      for (size_t i = 0; i < ICharge.size(); i++){
      	ChargeHist->Fill(ICharge.at(i));
      	EnergyHist->Fill(ICharge.at(i) * 1.556);
      }
           
      DetPerEvent->AddBinContent(totalEventNum, ICharge.size());
      ICharge.erase(ICharge.begin(), ICharge.end());
    }

    cout << "things filled with stuff. How long did it take? Write files." << endl;
    
    auto t_end = high_resolution_clock::now();
    duration<double,std::milli> time_total_ms(t_end-t_begin);
    cout << "\tEvent " << totalEventNum << " took " << time_total_ms.count() << " ms to process." << endl;
  }
  
  
  //testfile.close();
  f_output.Write();
  f_output.Close();  
  myfile.close();

  return 0; 
}
