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

#include <boost/lexical_cast.hpp>

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
  // declare important variables
  //------------------------------

  // Waveform class holds all data, performs noise filtering, baselining etc..
  rad_analysis::Waveforms nfrw;
  
  // detector/data varaiables
  int i_max_events = 250;
  int i_max_events_per_run = 1000000;
  bool b_noise_filter;
  string detector;
  string intag;
  vector<double> parameter_list(12);
  vector<string> filenames;
  ifstream myfile; // file of filenames by their full paths
  ifstream pfile("param.txt");

  // variables to control file list loop
  int i_total_event_num = 0;
  int i_write_count = 0;
  int i_previous_run = -99999;
  int i_run_event_count = 0;
  int i_this_run_num;
  int i_this_subrun_num;
  int i_this_event_num;
  int i_candidate_index;

  // vectors that hold properties of the event and candidates
  vector<float> TEinfo(3);
  vector<int> integration_window;    // this is used to record individual canidate total signal
  vector<float> total_charge;        // this lists total charge for all candidates in an event
  vector< vector<char> > temp_data;  // this vector holds the result from Signal_Select()
  vector<double> candidate_info(15); // this vector holds the candidate stats from Signal_Select()


  // ROOT variables
  TFile f_output("output.root","RECREATE");
  TTree ctree("ctree", "Candidate info");
  TTree etree("etree", "Event info");


  // set important parameters regarding data and detector settings
  //--------------------------------------------------------------

  //Begin setting all parameters from param.txt file
  pfile >> parameter_list.at(0) >> parameter_list.at(1) >> parameter_list.at(2) >> parameter_list.at(3) >> parameter_list.at(4)
	>> parameter_list.at(5) >> parameter_list.at(6) >> parameter_list.at(7) >> parameter_list.at(8) >> parameter_list.at(9)
	>> parameter_list.at(10) >> i_max_events >> i_max_events_per_run >> parameter_list.at(11) 
	>> nfrw.data_properties.MCCX >> detector >> intag >> nfrw.data_properties.SCALE;   
  pfile.close();

  InputTag product_tag  { intag };

  cout << "Detector: " << detector << endl;
  // cout << "Threshold: " << parameter_list[3] << " ADCs" << endl;

  //Set parameters for the detector
  if (detector == "MicroBooNE"){
    nfrw.detector_properties.DET = 0;    //microboone is detector 0 here
    nfrw.detector_properties.CHN = 8256; //it has 8256 channels
    nfrw.detector_properties.NCW = 3456; //3456 of them are collection
    nfrw.detector_properties.NTT = 6400; //look at 6400 time ticks
    nfrw.detector_properties.ntpcs = 1;  //only 1 tpc
    
    //set noise filter option and indicate if reco
    if (intag == "nfspl1:raw" || intag == "butcher" || intag == "wcNoiseFilter" || intag == "driftWC:orig"){
      b_noise_filter = 0;
      cout << intag << " selected. Noise filter off." << endl;      
    }

    else if (intag == "nfspl1:wiener" || intag == "nfspl1:gauss" || intag == "caldata" ||
	     intag == "simnfspl1:wiener" || intag == "simnfspl1:gauss"){
      b_noise_filter = 0;
      //parameter_list[3] /= 6.05724; //Improve: determine this dynamically
      nfrw.data_properties.recoused = true;

      if (intag == "simnfspl1:wiener" || intag == "simnfspl1:gauss"){
	nfrw.data_properties.simulated = true;
      }

      cout << intag << "selected. Noise filter off. recob used." << endl;
    }

    else if (intag == "daq"){
      b_noise_filter = 1;
      cout << intag << " selected. Noise filter on." << endl;
    }
  }

  else if (detector == "ProtoDUNE"){
    nfrw.detector_properties.DET = 1;
    nfrw.detector_properties.CHN = 15360; // total channels
    nfrw.detector_properties.NCW = 5760; // 2*480 ch * 6 APAs collection plane wires. 2*800*6 for induction
    nfrw.detector_properties.NTT = 6000; // time ticks to look at. Is this the total?
    nfrw.detector_properties.ntpcs = 12;
    b_noise_filter = 1;
  }
  else if (detector == "ProtoDUNE-HD"){
    nfrw.detector_properties.DET = 2;
    nfrw.detector_properties.CHN = 10240; // total channels
    nfrw.detector_properties.NCW = 3840; // 2*480 ch * 4 APAs collection plane wires. 2*800*4 for induction
    nfrw.detector_properties.NTT = 6000; // time ticks to look at. Is this the total?
    nfrw.detector_properties.ntpcs = 8; // 4 real TPCs corresponding to the 4 drift volumes, 4 "dummy" TPCs
    b_noise_filter = 1;
    // 0.512 microseconds per tick; 1.953125 MHz sampling rate
  }
  // now that the detector settings are set, fill the wire maps
  nfrw.Fill_Wire_Maps();
  // set TTree branches with vectors to be changed in Signal_Select()
  etree.Branch("tep", &TEinfo[0], "tep/F");
  etree.Branch("total", &TEinfo[1], "total/F");
  etree.Branch("remfrac", &TEinfo[2], "remfrac/F");
  etree.Branch("evIndex", &i_total_event_num, "evIndex/I");
  etree.Branch("evNum", &i_this_event_num, "evNum/I");
  etree.Branch("runNum", &i_this_run_num, "runNum/I");
  etree.Branch("subrun", &i_this_subrun_num, "subrun/I");

  ctree.Branch("ChargeVec", &integration_window);
  ctree.Branch("charge", &candidate_info[0], "charge/D");
  ctree.Branch("energy", &candidate_info[1], "energy/D"); 
  ctree.Branch("maxAmp", &candidate_info[2], "maxAmp/D");
  ctree.Branch("wireWidth", &candidate_info[3], "wireWidth/D");
  ctree.Branch("timeLength", &candidate_info[4], "timeLength/D");
  ctree.Branch("channel", &candidate_info[5], "channel/D");
  ctree.Branch("time", &candidate_info[6], "time/D");
  ctree.Branch("tOverThresh", &candidate_info[7], "tOverThresh/D");
  ctree.Branch("nearTrack", &candidate_info[8], "nearTrack/D");
  ctree.Branch("evNum", &i_this_event_num, "evNum/I");
  ctree.Branch("runNum", &i_this_run_num, "runNum/I");
  ctree.Branch("subrun", &i_this_subrun_num, "subrun/I");
  ctree.Branch("evIndex", &i_total_event_num, "evIndex/I");
  ctree.Branch("xWind", &parameter_list[5], "xWind/D");
  ctree.Branch("yWind", &parameter_list[6], "yWind/D");
  ctree.Branch("xCheck", &parameter_list[7], "xCheck/D");
  ctree.Branch("yCheck", &parameter_list[8], "yCheck/D");
  ctree.Branch("xCheckWidth", &parameter_list[9], "xCheckWidth/D");
  ctree.Branch("yCheckWidth", &parameter_list[10], "yCheckWidth/D");
  ctree.Branch("avglead", &candidate_info[9], "avglead/D");
  ctree.Branch("avgtail", &candidate_info[10], "avgtail/D");
  ctree.Branch("tpc", &candidate_info[11], "tpc/D");
  ctree.Branch("candidateIndex", &i_candidate_index, "candidateIndex/I");
  ctree.Branch("confirmed", &candidate_info[12], "confirmed/D");
  ctree.Branch("y", &candidate_info[13], "y/D");
  ctree.Branch("z", &candidate_info[14], "z/D");
  cout << "Set TTree branches" << endl;
  //// these histograms are for testing and can probably be removed soon
  TH1F *wire8049 = new TH1F("wire8049", "", nfrw.detector_properties.NTT, 0, nfrw.detector_properties.NTT);
  TH1F *ModeHist = new TH1F("ModeHist", "ADC baseline", nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN);
  TH1F *SDHist = new TH1F("SDHist","Standard Deviation per Wire", nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN);  
  TH1F *ChargeHist = new TH1F("ChargeHist","", 10000, -5000.5, 4999.5);
  TH1F *EnergyHist = new TH1F("EnergyHist","Ar39 Decay Energy", 10000, -5000.5 * 1.556, 4999.5 * 1.556); // This histogram is outdated
  TH1I *DetPerEvent = new TH1I("DetPerEvent","Detections of Ar39 Per Event", 1000, -0.5, 999.5);
  TH2F *DetectionMap = new TH2F("DetectionMap", "Map of Points over Threshold", 
                                nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN, 
                                nfrw.detector_properties.NTT, 0, nfrw.detector_properties.NTT); 
  TH2F *ChargeMap = new TH2F("ChargeMap", "", nfrw.detector_properties.NCW, 0, 
                             nfrw.detector_properties.NCW, nfrw.detector_properties.NTT, 0, nfrw.detector_properties.NTT); 
  TH2F *ChargeMap2 = new TH2F("ChargeMap2", "", nfrw.detector_properties.CHN, 0, 
                              nfrw.detector_properties.CHN, nfrw.detector_properties.NTT, 0, nfrw.detector_properties.NTT); 
  TH2F *ZDetectionMap = new TH2F("ZDetectionMap", "Zoomed In Detection Map", 300, 500, 800, 1000, 0, 1000);
  TH2F *EvViewI = new TH2F("EvViewI", "", ((double)nfrw.detector_properties.NCW/nfrw.detector_properties.ntpcs)*3, 0, 
                           ((double)nfrw.detector_properties.NCW/nfrw.detector_properties.ntpcs)*3, 
                           nfrw.detector_properties.NTT*2, -nfrw.detector_properties.NTT, nfrw.detector_properties.NTT); 
  TH1F *insd = new TH1F("insd","", nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN);  
  TH2F *EvViewO = new TH2F("EvViewO", "", ((double)nfrw.detector_properties.NCW/nfrw.detector_properties.ntpcs)*3, 0, 
                           ((double)nfrw.detector_properties.NCW/nfrw.detector_properties.ntpcs)*3, 
                           nfrw.detector_properties.NTT*2, -nfrw.detector_properties.NTT, nfrw.detector_properties.NTT); 
  TH1F *osd = new TH1F("osd","", nfrw.detector_properties.CHN, 0, nfrw.detector_properties.CHN);  
  // set some outdated histogram parameters
  EnergyHist->GetXaxis()->SetTitle("keV");
  EnergyHist->GetYaxis()->SetTitle("Number of Detections");
  EnergyHist->GetYaxis()->SetTitleOffset(1.2);



  //Ar39Study command syntax
  //"Ar39Study " simply runs Ar39Study with the parameters set in param.txt
  //"Ar39Study file.txt" overrides the default filelist.txt
  //"Ar39Study num1 num2" will override the event numbers and only look between event num1 and num2
  //"Ar39Study 1 2" will only look at event 1 and event 2
  if (argc == 1 || argc == 3){
    myfile.open ("filelist.txt");
    copy(istream_iterator<string>(myfile),istream_iterator<string>(),back_inserter(filenames));

    if (argc == 3){ //can tell range of events from command line
      i_max_events = boost::lexical_cast<int>(argv[2]);
    }
  }else if (argc == 2){ //tell which file to run, not from filelist (grid running)
    filenames.push_back(argv[1]);
  }

  //cout << "Seems to be working..." << endl; //debug all above this


  // start looping through data files and events
  //--------------------------------------------

  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
    auto t_begin = high_resolution_clock::now();

    i_candidate_index = 0;
    double dcounter = 0;    //counts dead channels (also gets divided into percent of detector area)
    int counter = 0;        //used in the track exclusion algorithm

    i_this_run_num = ev.eventAuxiliary().run();
    i_this_subrun_num = ev.eventAuxiliary().subRun();
    i_this_event_num = ev.eventAuxiliary().event();    

    TEinfo[0] = 0.0;
    TEinfo[1] = 0.0;
    TEinfo[2] = 0.0;

    // control the sequence/number of events
    if ((i_this_run_num != i_previous_run) && (i_run_event_count > 0)){
      i_previous_run = i_this_run_num;
      i_run_event_count = 1;
    }
    else if ((i_this_run_num == i_previous_run) && (i_run_event_count >= i_max_events_per_run)){
      continue;
    }
    else {
      i_previous_run = i_this_run_num;
      i_run_event_count++;
    }

    i_total_event_num++;
    i_write_count++;

    if (i_total_event_num > i_max_events) {
      cout << "Event number reached. Closing." << endl;  
      
      f_output.Write();
      f_output.Close();
      myfile.close();
      
      return 0;
    }

    //Keeps progress through long runs if program crashes
    if (i_write_count == 100){
      f_output.Write();
      i_write_count = 0;
    }

    //skip to this event
    if (argc > 2 && i_total_event_num < boost::lexical_cast<int>(argv[1]))
      continue;

    //// debug
    //if (i_total_event_num < 400)
    //continue;

    cout << "Processing "
  	 << "Run " << i_this_run_num << ", "
  	 << "Subrun " << i_this_subrun_num << ", "
  	 << "Event " << i_this_event_num << endl;



    //// Declare the raw and reco vectors to fill.
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
    	cout << "parameter_listam nfrw.data_properties.MCCX not set correctly." << endl;
    	return 0;
      }
    }

    DetectionMap->Reset();
    ZDetectionMap->Reset();
    ModeHist->Reset();
    SDHist->Reset();

    // now initialize data for this event
    nfrw.initialize_data (ev, allrawdigits_vec, wire_vec, wire_vec2, b_noise_filter, parameter_list);

    if (nfrw.data_properties.recoused){
      vector<raw::RawDigit>().swap(allrawdigits_vec); //Clear and free memory from ard_vec
      vector<int>().swap(badchan);                    //clear badchan after use
    }

    if (nfrw.f_reco_scale != 1.0) {parameter_list[3] /= nfrw.f_reco_scale;} //changes the scale of the threshold from adcs to the reconstructed units
 
    // // debug 
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
    temp_data = rad_analysis::Signal_Select(nfrw, total_charge, parameter_list, integration_window, candidate_info, TEinfo, ctree, i_candidate_index);

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
    if (parameter_list[0] == 1 && nfrw.detector_properties.DET == 0){

      //// For induction planes and collection plane
      for (unsigned short p = 0; p < nfrw.collection_channel.size(); p++){
	for (unsigned short s = 0; s < nfrw.collection_channel[p].size(); s++){
	  for(unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
	    DetectionMap->Fill(nfrw.collection_channel[p][s], t, static_cast<int>(temp_data.at(nfrw.collection_channel[p][s]).at(t)));
            // DetectionMap->Fill(nfrw.collection_channel[p][s], t, static_cast<int>(nfrw.adc_value.at(nfrw.collection_channel[p][s]).at(t)));
	  }
	}
      }      

      for (unsigned short p = 0; p < nfrw.induction_channel_1.size(); p++){
	for (unsigned short s = 0; s < nfrw.induction_channel_1[p].size(); s++){
	  for(unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
	    DetectionMap->Fill(nfrw.induction_channel_1[p][s], t, static_cast<int>(temp_data.at(nfrw.induction_channel_1[p][s]).at(t)));
	  }
	}
      }      

      for (unsigned short p = 0; p < nfrw.induction_channel_2.size(); p++){
	for (unsigned short s = 0; s < nfrw.induction_channel_2[p].size(); s++){
	  for(unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
	    DetectionMap->Fill(nfrw.induction_channel_2[p][s], t, static_cast<int>(temp_data.at(nfrw.induction_channel_2[p][s]).at(t)));
	  }
	}
      }  
    }

    //only for pdune right now
    else if (parameter_list[0] == 1 && (nfrw.detector_properties.DET == 1 || nfrw.detector_properties.DET == 2)){
      //Doesn't work for Induction Planes
      for (unsigned short p = 0; p < nfrw.collection_channel.size(); p++){
      	for (unsigned short s = 0, f = 0; s < nfrw.collection_channel[p].size(); s++){
      	  for(unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
      	    f = (p) * nfrw.detector_properties.NCW/12 + s;

      	    DetectionMap->Fill(f, t, static_cast<int>(temp_data.at(nfrw.collection_channel[p][s]).at(t)));
            // DetectionMap->Fill(f, t, static_cast<int>(nfrw.adc_value.at(nfrw.collection_channel[p][s]).at(t)));	    

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
    	    //   EvViewI->Fill(floor(log2(p)) * nfrw.collection_channel[p].size() + s, -t, static_cast<int>(temp_data.at(nfrw.collection_channel[p][s]).at(t)));
    	    // }
    	    // else if (p == 2 || p == 6 || p == 10){
    	    //   EvViewI->Fill(floor(log2(p - 1)) * nfrw.collection_channel[p].size() + s, t, static_cast<int>(temp_data.at(nfrw.collection_channel[p][s]).at(t)));
    	    // }
	    	  
      	  }
      	} 
      }

      // for (unsigned short p = 0; p < 3; p++){
      // 	for (unsigned short s = 0, f = 0; s < nfrw.collection_channel[p].size(); s++){
      // 	  for(unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
      // 	    f = (p + 1) * nfrw.detector_properties.NCW/12 + s;

      // 	    DetectionMap->Fill(f, t, static_cast<int>(temp_data.at(nfrw.collection_channel[p][s]).at(t)));
	  
      // 	  }
      // 	} 
      // }

      // //inner planes
      // for (unsigned short p = 1, f = 1; (unsigned)(p + 3) < nfrw.collection_channel.size(); p += 4, f++){
      // 	for (unsigned short s = 0, c = 0; s < nfrw.collection_channel[p].size(); s++){
      // 	  for(short t = 0; t < nfrw.detector_properties.NTT; t++){
      // 	    c = f * nfrw.collection_channel[p].size() + s;

      // 	    EvViewI->Fill(c, -t, static_cast<int>(temp_data.at(nfrw.collection_channel[p][s]).at(t)));
      // 	    EvViewI->Fill(c, t, static_cast<int>(temp_data.at(nfrw.collection_channel[p + 1][s]).at(t)));
	  
      // 	  }
      // 	} 
      // }

      // //outer planes
      // for (unsigned short p = 0, f = 1; (unsigned)(p + 4) < nfrw.collection_channel.size(); p += 4, f++){
      // 	for (unsigned short s = 0, c = 0; s < nfrw.collection_channel[p].size(); s++){
      // 	  for(short t = -nfrw.detector_properties.NTT + 1; t < nfrw.detector_properties.NTT; t++){
      // 	    c = f * nfrw.collection_channel[p].size() + s;

      // 	    if (t < 0)
      // 	      EvViewO->Fill(c, t, static_cast<int>(temp_data.at(nfrw.collection_channel[p][s]).at(-t)));

      // 	    else
      // 	      EvViewO->Fill(c, t, static_cast<int>(temp_data.at(nfrw.collection_channel[p + 3][s]).at(t)));
	  
      // 	  }
      // 	} 
      // }

      // for (int s = 5300; s < 5600; s++){
      // 	for (int t = 0; t < 1000; t++){
      // 	  ZDetectionMap->Fill(s - 4800, t, static_cast<int>(temp_data.at(s).at(t)));
      // 	}
      // }
    }

    //asemble data acquired from rad_analysis::Signal_Select()
    if (parameter_list[2] == 2){
      
      sort (total_charge.begin(), total_charge.end());

      for (size_t i = 0; i < total_charge.size(); i++){
      	ChargeHist->Fill(total_charge.at(i));
      	EnergyHist->Fill(total_charge.at(i) * 1.556);
      }
           
      DetPerEvent->AddBinContent(i_total_event_num, total_charge.size());
      total_charge.erase(total_charge.begin(), total_charge.end());
    }

    cout << "things filled with stuff. How long did it take? Write files." << endl;
    
    auto t_end = high_resolution_clock::now();
    duration<double,std::milli> time_total_ms(t_end-t_begin);
    cout << "\tEvent " << i_total_event_num << " took " << time_total_ms.count() << " ms to process." << endl;
  }
  
  
  //testfile.close();
  f_output.Write();
  f_output.Close();  
  myfile.close();

  return 0; 
}
