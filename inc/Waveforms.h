#ifndef WAVEFORMS_H
#define WAVEFORMS_H

#include <vector>

#include "gallery/Event.h"

#include "lardataobj/RawData/RawDigit.h" 
#include "lardataobj/RecoBase/Wire.h"


namespace rad_analysis {

  struct Detector_P {      // This struct lists the detector properties
    unsigned short NCW;   // Number of Collection Wires
    unsigned short CHN;   // Number of Channels
    unsigned short NTT;   // Number of Time Ticks
    char DET;             // Detector number for implementing in other detectors
    unsigned short ntpcs; // number of tpc modules in the detector
  };                      // DET = 0 is MicroBooNE; DET = 1 is ProtoDUNE


  struct Data_P {          // This struct keeps the data properties -- all set in parameter file
    unsigned short MCCX;  // MCC number because each MCC uses different data products
    bool recoused;        // Is this data reconstructed or raw?
    bool simulated;       // Is this data simulated?
    float SCALE;          // This must be used for a rounding bug in fourier space 
  };                      //  noise-filtered raw data products (ADC = 12 bit integer, but
  //  fourier transform outputs a real number, so rounds to ADC

  // Main class keeping info on the waveforms on each channel
  class Waveforms {
   protected:

    std::vector<int> m_mode_vector; // used in calculating the mode of a waveform
    std::vector<int> m_median_vector;  // used in calculating the median of a waveform
    std::vector<int> m_bad_chan; // used to keep track of bad channels
    std::vector< std::vector<float> > m_region_sd; // Standard deviation of subsets of waveform

    unsigned short m_i_sd_counter;                 // used in calculating m_region_sd
    unsigned short m_i_counter;                    // general short integer counter

   public:

    Detector_P detector_properties;
    Data_P     data_properties;


    std::vector<unsigned short> inner_tpc; // Numbering scheme of inner tpcs
    std::vector<unsigned short> outer_tpc; // Numbering scheme of outer tpcs (pdune geometry)

    //maps wire number to channel number on collection wires
    std::vector< std::vector<short> > collection_channel; //(ntpcs);

    //induction wires
    std::vector< std::vector<short> > induction_channel_1; //(ntpcs);
    std::vector< std::vector<short> > induction_channel_2; //(ntpcs);

    //maps channel number to tpc, wire, plane
    std::vector< std::vector<short> > channel_to_wire; //(nfrw.detector_properties.CHN, vector<short>(3));


    std::vector<float> sd_vector;    // Stores standard deviation on each channel
    std::vector<float> baseline;     // Stores baseline on each channel
    std::vector<bool> dead_channel;  // Boolean whether each channel is dead or not
    std::vector<bool> track_channel; // Does a track cross this channel?


    std::vector< std::vector<float> > adc_value; // Stores ADC waveform 
    std::vector< std::vector<float> > reco_product_value;  // Stores reco waveform (if available and instructed to & usually, hopefully gauss filter product)


    std::vector< std::vector<float> > geometry_info;  // stores geometry info read from detector geom file
    std::vector< std::vector<short> > plane_wire_intersections; // all induction intersections with collection wires


    std::vector< std::vector<double> > candidate_info; // 39Ar candidate info stored here so that it can be easily passed around (total energy etc...)


    float f_reco_scale; // converts ADC threshold to reco units


    // Fill std::vectors as maps 
    // Assumes the detector_properties and data_properties are already filled!  
    // Basic guide for how these vectors work as maps:
    // collection_channel  : [tpc][wire] --> [channel] | channel is a collection wire
    // induction_channel_1 : [tpc][wire] --> [channel] | channel is the first induction wire
    // induction_channel_2 : [tpc][wire] --> [channel] | channel is the second induction wire
    //                                       [tpc  ]
    // channel_to_wire     :   [channel] --> [wire ]
    //                                       [plane]
    void Fill_Wire_Maps(); 
  

    void initialize_data (const gallery::Event& ev,
                          const std::vector<raw::RawDigit>& ard_vec,
                          const std::vector<recob::Wire>& wire_vec,
                          const std::vector<recob::Wire>& wire_vec2,
                          const bool& b_noise_filter,
                          const std::vector<double>& Par);

  };
}

#endif
