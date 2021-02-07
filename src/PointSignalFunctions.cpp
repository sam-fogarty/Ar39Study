#include "../inc/PointSignalFunctions.h"

// #include "Waveforms.h"
#include <iostream>

double bline(double x1, double y1, double x2, double y2, double x){
  double m, b;

  m = (y1 - y2) / (x1 - x2);
  b = -(m * x1) + y1;

  return ((m * x) + b);
} 

rad_analysis::Subarea::Subarea(short a, short b, short c, short d) // Subarea constructor
  : wire1(a)                                                       // wire1 = a
  , wire2(c)                                                       // time1 = b
  , time1(b)                                                       // wire2 = c
  , time2(d)                                                       // time2 = d
{ }

void rad_analysis::Subarea::Draw_Subarea(std::vector< std::vector<char> >& t_data, char val){
  for (int s = wire1; s <= wire2; s++){
    for (int t = time1; t <= time2; t++){
      if (s < (int)t_data.size() && s >= 0 && t < (int)t_data[s].size() && t >= 0){      
        t_data[s][t] = val;
      }
    }
  }
}

bool rad_analysis::Subarea::Check_Track(const std::vector< std::vector<char> >& t_data){
  for (int s = wire1; s <= wire2; s++){
    for (int t = time1; t <= time2; t++){
      if (s < (int)t_data.size() && s >= 0 && t < (int)t_data[s].size() && t >= 0){      

        if (t_data[s][t] == 1) 
          return true;

      }
    }
  }

  return false;
}

bool rad_analysis::Subarea::Check(const std::vector< std::vector<char> >& t_data){
  for (int s = wire1; s <= wire2; s++){
    for (int t = time1; t <= time2; t++){
      if (s < (int)t_data.size() && s >= 0 && t < (int)t_data[s].size() && t >= 0){      

        if (t_data[s][t] != 0)
          return false;

      }
      
      else return false;
    }
  }

  return true;
}

std::vector< std::vector<char> > rad_analysis::Signal_Select(rad_analysis::Waveforms& nfrw,
                                                             std::vector<float>& ICharge,
                                                             std::vector<double>& Par,
                                                             std::vector<int>& IntWindow,
                                                             std::vector<double>& c_info,
                                                             std::vector<float>& TEinfo,
                                                             TTree& ctree,
                                                             int& cind,
                                                             const gallery::Event& ev){

  std::cout << "Signal_Select called" << std::endl;

  // These are used for the asymmetric windowing; these values are hard-coded for now, but we'll want it to depend on the detector/data product
  const short c_center_range[2] = {-16, 12};     // signal center wire
  const short c_side_range[2] = {-14, 10};       // signal side wires

  std::vector<short> coord;  //coordinates to points over threshold
  std::vector<short> fcoor;  //coordinates to points over threshold
  std::vector<short> ecoor;  //coordinates of the initial edges of track veto

  std::vector< std::vector<char> > threshold_area (nfrw.detector_properties.CHN, std::vector<char>(nfrw.detector_properties.NTT, 0));      //map of points over threshold and dead channels
  std::vector< std::vector<char> > track_exclusion_map (nfrw.detector_properties.CHN, std::vector<char>(nfrw.detector_properties.NTT,0));  //map of excluded area due to muon track
  //std::vector< std::vector<float> > adc_value (nfrw.detector_properties.CHN, std::vector<float>(nfrw.detector_properties.NTT, 0)); //actual std::vector of waveforms for event

  // Parameter variables from parameter file (Par is read and filled in main() )
  float f_parameter_threshold   = Par[3];  // threshold for radiological candidates
  bool  b_parameter_noise_study = Par[1];  // 0 is operating normally; 1 randomly selects points in the detector to sample noise
  int   i_parameter_x_window    = Par[5];  // horizontal width parameter for the actual integration window
  int   i_parameter_y_window    = Par[6];  // vertical width parameter for the integration window 
  int   i_parameter_x_buffer    = Par[7];  // buffer between integration window and track check boxes
  int   i_parameter_y_buffer    = Par[8]; 
  int   i_parameter_x_width     = Par[9];  // width of the track check boxes
  int   i_parameter_y_width     = Par[10];

  unsigned short j = 0;
  unsigned short k = 0;
  unsigned short l = 0;

  unsigned short dt = 7;

  srand(time(NULL));
  
  std::cout << "test 1" << std::endl;

  //blank out the edges of each plane so that there isn't accidental integration across two planes
  for (unsigned short a = 0; a < nfrw.collection_channel.size(); a++){ //assumes nfrw.collection_channel and iwx are all the same size
    for (unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
      threshold_area[nfrw.collection_channel[a][0]][t] = -1;
      threshold_area[nfrw.collection_channel[a][nfrw.collection_channel[a].size() - 1]][t] = -1;  // -1 indicates dead or otherwise forbidden
                                                                                                  // signals too close to these are tossed out
      if (nfrw.detector_properties.DET != 1){
	threshold_area[nfrw.induction_channel_1[a][0]][t] = -1;
	threshold_area[nfrw.induction_channel_1[a][nfrw.induction_channel_1[a].size() - 1]][t] = -1;
	threshold_area[nfrw.induction_channel_2[a][0]][t] = -1;
	threshold_area[nfrw.induction_channel_2[a][nfrw.induction_channel_2[a].size() - 1]][t] = -1;
      }
    }
  }

  for (unsigned short a = 0; a < nfrw.collection_channel.size(); a++){
    for (unsigned short s = 0; s < nfrw.collection_channel[a].size(); s++){
      for (unsigned short t = 0; t < i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width; t++){
	threshold_area[nfrw.collection_channel[a][s]][nfrw.detector_properties.NTT - t - 1] = -1;
	threshold_area[nfrw.collection_channel[a][s]][t] = -1;
      }
    }
  }

  // std::cout << "test 2" << std::endl;

  for (unsigned short s = 0; s < nfrw.adc_value.size(); s++){
    if (nfrw.detector_properties.DET != 1 && nfrw.channel_to_wire[s][2] != 2){ //select induction planes
    //if (nfrw.channel_to_wire[s][2] == 0){ //select only U plane
      if (!nfrw.data_properties.recoused){
	for (unsigned short t = 0; t < nfrw.adc_value[s].size() - dt; t++){
	  //if ((nfrw.adc_value[s][t] - nfrw.adc_value[s][t + dt]) > f_parameter_threshold * (4.0/3.0)){
	  if ((nfrw.adc_value[s][t] - nfrw.adc_value[s][t + dt]) > f_parameter_threshold){
	    if (b_parameter_noise_study != 1 && nfrw.adc_value[s][t + dt] < 0){
	      coord.push_back(s);
	      coord.push_back(t + dt);
	      threshold_area[s][t + dt] = 1;
	    }
	  }

	  if (nfrw.dead_channel[s] == 1)
	    threshold_area[s][t] = -1;
	}
      }

      else if (nfrw.data_properties.recoused){
	for (unsigned short t = 0; t < nfrw.adc_value[s].size() - dt; t++){
	  //if ((nfrw.adc_value[s][t] - nfrw.adc_value[s][t + dt]) > f_parameter_threshold * (4.0/3.0)){

	  if (nfrw.adc_value[s][t] > f_parameter_threshold){
	    coord.push_back(s);
	    coord.push_back(t);
	    threshold_area[s][t] = 1;
	  }

	  if (nfrw.dead_channel[s] == 1)
	    threshold_area[s][t] = -1;
	}
      }
    }

    if (nfrw.channel_to_wire[s][2] == 2){ //select collection planes
      for (unsigned short t = 0; t < nfrw.adc_value[s].size(); t++){
	if (abs(nfrw.adc_value[s][t]) > f_parameter_threshold){
	  threshold_area[s][t] = 1;

	  if (b_parameter_noise_study != 1 && nfrw.adc_value[s][t] > 0){
	    coord.push_back(s);
	    coord.push_back(t);
	  }
	}

	if (nfrw.dead_channel[s] == 1)
	  threshold_area[s][t] = -1;

	// else if (nfrw.adc_value[s][t] < -f_parameter_threshold && b_parameter_noise_study == 1){
	//   threshold_area[s][t] = -1;
	// }
	// else if (nfrw.adc_value[s][t] < -20)
	// 	threshold_area[s][t] = -1;
      }

      //if enabled, select random coordinates (noise study, collection plane only)
      if (b_parameter_noise_study == 1 && l < 1000){
	if (j < 1000 / nfrw.collection_channel.size()){
	  // std::cout << "test fill coord" << std::endl;

	  coord.push_back(nfrw.collection_channel[k][rand() % nfrw.collection_channel[k].size()]);
	  coord.push_back(rand() % nfrw.detector_properties.NTT);
	  j++;
	}

	else if (j == (1000 / nfrw.collection_channel.size())){
	  j = 0;
	  k++;
	}
      }

      l++;
    }
  }

  std::cout << "Initial Coord: " << coord.size() << std::endl;

  if (Par[2] == 2){ 

    std::cout << "Track exclusion started: " << std::endl;

    int dim[4];
    bool check[4];
    bool ch;
    
    //track exclusion continued
    for (size_t i_coord = 0; i_coord < coord.size(); i_coord += 2){
      short s = coord.at(i_coord);
      short t = coord.at(i_coord + 1);
      int dums, dumt, tMax = 0;

      for (int i = 0; i < 3; i++){
    	check[i] = false;
	dim[i] = 0;
      }
      
      ch = 0;

      //center on max or minimum
      for (int x = s - i_parameter_x_window; x <= s + i_parameter_x_window; x++){
	for (int y = t - i_parameter_y_window; y <= t + i_parameter_y_window; y++){
	  if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){	    

	    //two different selections; one for induction, and one for collection
	    if (nfrw.adc_value.at(x).at(y) > tMax && nfrw.channel_to_wire[x][2] == 2){ //collection centering
	      tMax = nfrw.adc_value.at(x).at(y);
	      dums = x;
	      dumt = y;
	    }

            //induction centering
	    else if (nfrw.adc_value.at(x).at(y) < tMax && nfrw.channel_to_wire[x][2] != 2 && !nfrw.data_properties.recoused){ 
	      tMax = nfrw.adc_value.at(x).at(y);
	      dums = x;
	      dumt = y;
	    }

            //reco induction centering
            else if (nfrw.adc_value.at(x).at(y) > tMax && nfrw.channel_to_wire[x][2] != 2 && nfrw.data_properties.recoused){ 
	      tMax = nfrw.adc_value.at(x).at(y);
	      dums = x;
	      dumt = y;
	    }
	  }
	}	
      }

      if (nfrw.detector_properties.DET == 1){ //PDUNE-SP Only
      	//Stuck Bit Filter First-Pass
      	if (dums > 0 && (unsigned int)dums < threshold_area.size() && dumt > 0 && dumt < nfrw.detector_properties.NTT - 1){
      	  // if (nfrw.adc_value.at(dums).at(dumt) - nfrw.adc_value.at(dums).at(dumt - 1) > nfrw.adc_value.at(dums).at(dumt) / 2.0 || 
      	  //     nfrw.adc_value.at(dums).at(dumt) - nfrw.adc_value.at(dums).at(dumt + 1) > nfrw.adc_value.at(dums).at(dumt) / 2.0){
	  if (nfrw.adc_value.at(dums).at(dumt) - nfrw.adc_value.at(dums).at(dumt - 1) > 8 && 
      	      nfrw.adc_value.at(dums).at(dumt) - nfrw.adc_value.at(dums).at(dumt + 1) > 8){

	    if (i_coord + 2 < coord.size()){
	      coord.erase(coord.begin() + i_coord, coord.begin() + i_coord + 2);
	    }

            else{
	      coord.erase(coord.begin() + i_coord, coord.begin() + i_coord + 1);
	      coord.erase(coord.begin() + i_coord, coord.begin() + i_coord + 1);
	    }
	    
      	    i_coord -= 2;
      	    continue; 
      	  }
      	}
      }

      if (nfrw.track_channel.at(s) == 1){  //only do this for wires with a muon track (also only collection plane)

	//left box
        rad_analysis::Subarea left_box(s - (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width),
                                       t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width),
                                       s - (i_parameter_x_window + i_parameter_x_buffer),
                                       t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width) );
	
	//top box
        rad_analysis::Subarea top_box(s - (i_parameter_x_window + i_parameter_x_buffer),
                                      t + (i_parameter_y_window + i_parameter_y_buffer),
                                      s + (i_parameter_x_window + i_parameter_x_buffer),
                                      t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width) );
	
	//bottom box
        rad_analysis::Subarea bottom_box(s - (i_parameter_x_window + i_parameter_x_buffer),
                                         t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width),
                                         s + (i_parameter_x_window + i_parameter_x_buffer),
                                         t - (i_parameter_y_window + i_parameter_y_buffer) );
	
	//right box
        rad_analysis::Subarea right_box(s + (i_parameter_x_window + i_parameter_x_buffer),
                                        t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width),
                                        s + (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width),
                                        t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width) );


	//left region to check around candidate
        check[0] = left_box.Check_Track(threshold_area);
        if (check[0])
          ch = true;

	// for (int x = s - (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width - 1); 
        //      x <= s - (i_parameter_x_window + i_parameter_x_buffer); x++){
	//   for (int y = t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); 
        //        y <= t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y++){ 
	                                                                                           
    	//     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){
    	//       //threshold_area.at(x).at(y) = -1;                                                            
    	//       if (nfrw.adc_value.at(x).at(y) > f_parameter_threshold && track_exclusion_map[x][y] != 1){
    	// 	check[0] = true;
    	// 	threshold_area.at(x).at(y) = 1;
    	// 	ch = 1;
    	//       }
    	//       //threshold_area.at(x).at(y) = -1;
    	//     }
    	//   }
    	// }
	
	//top region to check around candidate
        check[1] = top_box.Check_Track(threshold_area);
        if (check[1])
          ch = true;

	// for (int x = s - (i_parameter_x_window + i_parameter_x_buffer); x <= s + (i_parameter_x_window + i_parameter_x_buffer); x++){
	//   for (int y = t + (i_parameter_y_window + i_parameter_y_buffer); y <= t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y++){
	
    	//     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){
    	//       //threshold_area.at(x).at(y) = -1;                                                            
    	//       if (nfrw.adc_value.at(x).at(y) > f_parameter_threshold && track_exclusion_map[x][y] != 1){
    	// 	check[1] = true;
    	// 	threshold_area.at(x).at(y) = 1;
    	// 	ch = 1;
    	//       }
    	//       //threshold_area.at(x).at(y) = -1;
    	//     }
    	//   }
    	// }
	
	//bottom region to check around candidate
        check[2] = bottom_box.Check_Track(threshold_area);
        if (check[2])
          ch = true;


	// for (int x = s - (i_parameter_x_window + i_parameter_x_buffer); x <= s + (i_parameter_x_window + i_parameter_x_buffer); x++){
	//   for (int y = t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y <= t - (i_parameter_y_window + i_parameter_y_buffer); y++){

    	//     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){
    	//       //threshold_area.at(x).at(y) = -1;                                                            
    	//       if (nfrw.adc_value.at(x).at(y) > f_parameter_threshold && track_exclusion_map[x][y] != 1){
    	// 	check[2] = true;
    	// 	threshold_area.at(x).at(y) = 1;
    	// 	ch = 1;
    	//       }
    	//       //threshold_area.at(x).at(y) = -1;
    	//     }
    	//   }
    	// }
	
	//right region to check around candidate
        check[3] = right_box.Check_Track(threshold_area);
        if (check[3])
          ch = true;


	// for (int x = s + (i_parameter_x_window + i_parameter_x_buffer); 
        //      x <= s + (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width - 1); x++){
	//   for (int y = t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); 
        //        y <= t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y++){

    	//     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){
    	//       //threshold_area.at(x).at(y) = -1;                                                            
    	//       if (nfrw.adc_value.at(x).at(y) > f_parameter_threshold && track_exclusion_map[x][y] != 1){
    	// 	check[3] = true;
    	// 	threshold_area.at(x).at(y) = 1;
    	// 	ch = 1;
    	//       }
    	//     }
    	//   }		  
    	// }

    	//rudimentary direction detection/constructing veto bubble
    	if (ch == 1){

    	  if (check[0] || check[3]){ // If charge crosses left or right region
    	    dim[1] = 60;             // draw lines up
	    dim[2] = 60;             // and down
	  }

	  if (check[1] || check[2]){ // If charge crosses top or bottom region
    	    dim[0] = 10;             // draw lines left
	    dim[3] = 10;             // and right
	  }

    	  //refine the track exclusion map so that it's more efficient to draw
    	  for (int x = s - dim[0]; x <= s + dim[3]; x++){
    	    for (int y = t - dim[2]; y <= t + dim[1]; y++){
    	      if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){

		if (track_exclusion_map[x][y] != 1){
		  if (x == s - dim[0] || x == s + dim[3] || y == t - dim[2] || y == t + dim[1]){
		    ecoor.push_back(x);
		    ecoor.push_back(y);
		  }

		  track_exclusion_map[x][y] = 1;
		  TEinfo[0] = TEinfo[0] + 1.0;
		}

    	        //if (threshold_area.at(x).at(y) == 0){
    		  // threshold_area.at(x).at(y) = 0.5;
    		  //threshold_area.at(x).at(y) = -1;

		  // track_exclusion_map.at(x).at(y) = 1;
	          // TEinfo[0] = TEinfo[0] + 1.0;
    		  //threshold_area.at(x).at(y) = -1;
		  //}
              }
    	    }
    	  }
    	}
      }
    }


    // vertical line loop
    for (size_t i_coord = 0; i_coord < ecoor.size(); i_coord += 2){
      std::vector<bool> normal_direc(2, false);

      short s = ecoor.at(i_coord);
      short t = ecoor.at(i_coord + 1);
      unsigned int test = 0;

      if (track_exclusion_map[s][t + 1] == 0){
    	normal_direc[0] = true;
      }

      else if (track_exclusion_map[s][t - 1] == 0){
    	normal_direc[1] = true;
      }

      for (short y = -(Par[11] - 10) * 6 * normal_direc[1]; y < (Par[11] - 10) * 6 * normal_direc[0]; y++){
    	for (short x = 0; x <= 6; x++){
    	  if (x + s < (int)track_exclusion_map.size() && x + s >= 0 && y + t < (int)track_exclusion_map[s].size() && y + t >= 0){
    	    track_exclusion_map[s + x][t + y] = 2;
    	  }
    	}
      }

      std::vector<bool>().swap(normal_direc); //clear memory used by normal
    }

    // horizontal line loop
    for (size_t i_coord = 0; i_coord < ecoor.size(); i_coord += 2){
      std::vector<bool> normal_direc(2, false);

      short s = ecoor.at(i_coord);
      short t = ecoor.at(i_coord + 1);
      unsigned int test = 0;

      if (track_exclusion_map[s + 1][t] == 0 && track_exclusion_map[s - 1][t] != 0){
      	normal_direc[0] = true;
      }

      else if (track_exclusion_map[s - 1][t] == 0 && track_exclusion_map[s + 1][t] != 0){
      	normal_direc[1] = true;
      } 

      for (short x = -(Par[11] - 10) * normal_direc[1] / 2; x < (Par[11] - 10) * normal_direc[0] / 2; x++){
    	for (short y = -(Par[11] - 10) * 6; y <= (Par[11] - 10) * 6; y++){
    	  if (x + s < (int)track_exclusion_map.size() && x + s >= 0 && y + t < (int)track_exclusion_map[s].size() && y + t >= 0){
    	    track_exclusion_map[s + x][t + y] = 1;
    	  }
    	}
      }

      std::vector<bool>().swap(normal_direc); //clear memory used by normal
    }

    std::cout << "Track Exclusion finished. Starting candidate selection from coordinates. " << std::endl;
    std::cout << "Coordinates found: " << coord.size()/2 << std::endl;

    // candidate selection from coordinates 
    for (size_t i_coord = 0; i_coord < coord.size(); i_coord += 2){

      bool check[2];
      bool wireWidthCheck = false;

      short s = coord.at(i_coord);
      short t = coord.at(i_coord + 1);

      int dums, temps, dumt, tempt;

      int tMax = 0, tOverThresh = 0, tWidth = 0, wireWidth = 0;
      int tWidthCounter = 0;

      double avgf = 0, avgt = 0, charge = 0;
      short cf = 0, ct = 0;

      // std::vector<double> sidebandf((i_parameter_x_window * 2) + 1);
      // std::vector<double> sidebandt((i_parameter_x_window * 2) + 1);

      c_info[8] = 0;

      for (int i = 0; i < 2; i++){
	check[i] = true;
      }

      //center on max or minimum
      for (int x = s - i_parameter_x_window; x <= s + i_parameter_x_window; x++){
	for (int y = t - i_parameter_y_window; y <= t + i_parameter_y_window; y++){
	  if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){	    

	    //two different selections; one for induction, and one for collection
	    if (nfrw.adc_value.at(x).at(y) > tMax && nfrw.channel_to_wire[x][2] == 2){ //collection centering
	      tMax = nfrw.adc_value.at(x).at(y);
	      dums = x;
	      dumt = y;
	    }

	    else if (nfrw.adc_value.at(x).at(y) < tMax && nfrw.channel_to_wire[x][2] != 2 && !nfrw.data_properties.recoused){ //induction centering
	      tMax = nfrw.adc_value.at(x).at(y);
	      dums = x;
	      dumt = y;
	    }

            else if (nfrw.adc_value.at(x).at(y) > tMax && nfrw.channel_to_wire[x][2] != 2 && nfrw.data_properties.recoused){ //reco induction centering
	      tMax = nfrw.adc_value.at(x).at(y);
	      dums = x;
	      dumt = y;
	    }

	    if (threshold_area.at(x).at(y) == 1){
	      tOverThresh++;
	      tWidthCounter++;
	      wireWidthCheck = true;
	    }	    
	  }
	}
	
	if (wireWidthCheck == true){
	  wireWidth++;
	  wireWidthCheck = false;
	}

	if (tWidthCounter > tWidth){
	  tWidth = tWidthCounter;
	  tWidthCounter = 0;
	}
      }

      temps = s;
      tempt = t;

      if (b_parameter_noise_study == 0){ //but don't center for the noise study 
	s = dums;
	t = dumt;
      }

      if (track_exclusion_map[temps][tempt] != 0) {c_info[8] = 1;} // flag if candidate is in track region


      //left box
      rad_analysis::Subarea left_box(s - (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width),
                                     t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width),
                                     s - (i_parameter_x_window + i_parameter_x_buffer),
                                     t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width) );
	
      //top box
      rad_analysis::Subarea top_box(s - (i_parameter_x_window + i_parameter_x_buffer),
                                    t + (i_parameter_y_window + i_parameter_y_buffer),
                                    s + (i_parameter_x_window + i_parameter_x_buffer),
                                    t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width) );
	
      //bottom box
      rad_analysis::Subarea bottom_box(s - (i_parameter_x_window + i_parameter_x_buffer),
                                       t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width),
                                       s + (i_parameter_x_window + i_parameter_x_buffer),
                                       t - (i_parameter_y_window + i_parameter_y_buffer) );
	
      //right box
      rad_analysis::Subarea right_box(s + (i_parameter_x_window + i_parameter_x_buffer),
                                      t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width),
                                      s + (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width),
                                      t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width) );
      
      //left region to check around candidate
      check[0] = left_box.Check(threshold_area);
      if (!check[0])
        threshold_area[temps][tempt] = -1;

      // for (int x = s - (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width - 1); x <= s - (i_parameter_x_window + i_parameter_x_buffer); x++){
      //   for (int y = t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y <= t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y++){
	  
      //     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){
      //       if (threshold_area.at(x).at(y) != 0){
      //         threshold_area.at(temps).at(tempt) = -1;
      //         check[0] = false;
      //       }
	    
      //       if (track_exclusion_map.at(x).at(y) != 0)
      //         c_info[8] = 1; 
      //     }

      //     else
      //       check[0] = false;
      //   }
      // }
      
      //top region to check around candidate
      check[0] = top_box.Check(threshold_area);
      if (!check[0])
        threshold_area[temps][tempt] = -1;

      // for (int x = s - (i_parameter_x_window + i_parameter_x_buffer); x <= s + (i_parameter_x_window + i_parameter_x_buffer); x++){
      //   for (int y = t + (i_parameter_y_window + i_parameter_y_buffer); y <= t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y++){
	  
      //     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){
      //       if (threshold_area.at(x).at(y) != 0){
      //         threshold_area.at(temps).at(tempt) = -1;
      //         check[0] = false;
      //       }
	    
      //       if (track_exclusion_map.at(x).at(y) != 0) // If this candidate is in a track exclusion zone,
      //         c_info[8] = 1;                          // flag it, but still record everything; it can be cut from the ttree.
      //     }

      //     else
      //       check[0] = false;
      //   }
      // }
      
      //bottom region to check around candidate
      check[0] = bottom_box.Check(threshold_area);
      if (!check[0])
        threshold_area[temps][tempt] = -1;

      // for (int x = s - (i_parameter_x_window + i_parameter_x_buffer); x <= s + (i_parameter_x_window + i_parameter_x_buffer); x++){
      //   for (int y = t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y <= t - (i_parameter_y_window + i_parameter_y_buffer); y++){
	  
      //     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){
      //       if (threshold_area.at(x).at(y) != 0){
      //         threshold_area.at(temps).at(tempt) = -1;
      //         check[0] = false;
      //       }
            
      //       if (track_exclusion_map.at(x).at(y) != 0) 
      //         c_info[8] = 1;                          
      //     }

      //     else
      //       check[0] = false;
      //   }
      // }
      
      //right region to check around candidate
      check[0] = right_box.Check(threshold_area);
      if (!check[0])
        threshold_area[temps][tempt] = -1;

      // for (int x = s + (i_parameter_x_window + i_parameter_x_buffer); x <= s + (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width - 1); x++){
      //   for (int y = t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y <= t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y++){

      //     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){
      //       if (threshold_area.at(x).at(y) != 0){
      //         threshold_area.at(temps).at(tempt) = -1;
      //         check[0] = false;
      //       }
	    
      //       if (track_exclusion_map.at(x).at(y) != 0)
      //         c_info[8] = 1; 
      //     }

      //     else
      //       check[0] = false;
      //   }		  
      // }
      
      if (check[0] == true){
	//left box
        left_box.Draw_Subarea(threshold_area, -1);

        // rad_analysis::Subarea left_box(s - (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width),
        //                                t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width),
        //                                s - (i_parameter_x_window + i_parameter_x_buffer),
        //                                t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width) );
	// for (int x = s - (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width); x <= s - (i_parameter_x_window + i_parameter_x_buffer); x++){
	//   for (int y = t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y <= t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y++){
	//     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0)
	//       threshold_area.at(x).at(y) = -1;
	//   }
	// }
	
	//top box
        top_box.Draw_Subarea(threshold_area, -1);

        // rad_analysis::Subarea top_box(s - (i_parameter_x_window + i_parameter_x_buffer),
        //                               t + (i_parameter_y_window + i_parameter_y_buffer),
        //                               s + (i_parameter_x_window + i_parameter_x_buffer),
        //                               t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width) );
	// for (int x = s - (i_parameter_x_window + i_parameter_x_buffer); x <= s + (i_parameter_x_window + i_parameter_x_buffer); x++){
	//   for (int y = t + (i_parameter_y_window + i_parameter_y_buffer); y <= t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y++){
	//     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0)
	//       threshold_area.at(x).at(y) = -1;
	//   }
	// }
	
	//bottom box
        bottom_box.Draw_Subarea(threshold_area, -1);

        // rad_analysis::Subarea bottom_box(s - (i_parameter_x_window + i_parameter_x_buffer),
        //                                  t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width),
        //                                  s + (i_parameter_x_window + i_parameter_x_buffer),
        //                                  t - (i_parameter_y_window + i_parameter_y_buffer) );
        // for (int x = s - (i_parameter_x_window + i_parameter_x_buffer); x <= s + (i_parameter_x_window + i_parameter_x_buffer); x++){
	//   for (int y = t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y <= t - (i_parameter_y_window + i_parameter_y_buffer); y++){
	//     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0)
	//       threshold_area.at(x).at(y) = -1;
	//   }
	// }
	
	//right box
        right_box.Draw_Subarea(threshold_area, -1);

        // rad_analysis::Subarea right_box(s + (i_parameter_x_window + i_parameter_x_buffer),
        //                                 t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width),
        //                                 s + (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width),
        //                                 t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width) );
        // for (int x = s + (i_parameter_x_window + i_parameter_x_buffer); x <= s + (i_parameter_x_window + i_parameter_x_buffer + i_parameter_x_width); x++){
	//   for (int y = t - (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y <= t + (i_parameter_y_window + i_parameter_y_buffer + i_parameter_y_width); y++){
	//     if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0)
	//       threshold_area.at(x).at(y) = -1;
	//   }
	// }

	//fill an std::vector of the final coordinates to real candidates
	fcoor.push_back(s);
	fcoor.push_back(t);
      
	for (int x = s - i_parameter_x_window; x <= s + i_parameter_x_window; x++){
	  for (int y = t - i_parameter_y_window; y <= t + i_parameter_y_window; y++){
	    if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){
	      //average front and tail before signal
	      if (y > t && y + i_parameter_y_window < nfrw.detector_properties.NTT){
		avgf += nfrw.adc_value.at(x).at(y + i_parameter_y_window);
		cf++;
	      }
	      
	      if (y < t && y - i_parameter_y_window > 0){
		avgt += nfrw.adc_value.at(x).at(y - i_parameter_y_window);
		ct++;
	      }
	    }
	  }
	}

	nfrw.candidate_info.push_back(std::vector<double> (13, -1));
	nfrw.candidate_info.back()[0] = tMax;
	nfrw.candidate_info.back()[1] = wireWidth;
	nfrw.candidate_info.back()[2] = tWidth;
	nfrw.candidate_info.back()[3] = s;
	nfrw.candidate_info.back()[4] = t; 
	nfrw.candidate_info.back()[5] = tOverThresh;
	nfrw.candidate_info.back()[6] = avgf / cf;
	nfrw.candidate_info.back()[7] = avgt / ct;
	nfrw.candidate_info.back()[8] = c_info[8];
	nfrw.candidate_info.back()[9] = c_info[11];
      }
    }

    std::cout << "candidates looped through" << std::endl;
    
    //3D position reco
    if (nfrw.detector_properties.DET != 1){
      rad_analysis::Confirm_Candidates(fcoor, nfrw);
    }
    
    //Confirm_Candidates() will check which fcoor match with each other based on y position
    std::cout << "candidates confirmed by induction plane " << std::endl;

    for (size_t i_coord = 0; i_coord < fcoor.size(); i_coord += 2){
      unsigned char wire = 0;

      short range[2];

      int s = fcoor[i_coord];
      int t = fcoor[i_coord + 1];
    
      double lineval;
      
      std::vector<double> sidebandf;
      std::vector<double> sidebandt;

      std::vector<double> sf;
      std::vector<double> st; 

      bool stuck_bit = false;
      bool asymmetric;

      //Integrate Charge 
      //if (nfrw.channel_to_wire[s][2] == 2){ //Collection ONLY
      if (nfrw.channel_to_wire[s][2] >= 0){ //Integrate on any plane (for reco waveform or just signal shape on induction plane)
	for (int x = s - i_parameter_x_window; x <= s + i_parameter_x_window; x++){
	  sf.clear();
	  st.clear();

	  sidebandf.push_back(0);
	  sidebandt.push_back(0);


	  if (wire == i_parameter_x_window){ //if not central wire
	    range[0] = c_center_range[0];
	    range[1] = c_center_range[1];
	  }else{
	    range[0] = c_side_range[0];
	    range[1] = c_side_range[1];
	  }

	  if (x < (int)threshold_area.size() && x >= 0 && t + range[1] + 25 < (int)threshold_area.at(s).size() && t + range[1] + 25 >= 0){
	    for (int p = t + range[1] + 1; p <= t + range[1] + 25; p++){
	      sidebandf.at(sidebandf.size() - 1) += nfrw.adc_value.at(x).at(p);
	      sf.push_back(nfrw.adc_value.at(x).at(p));
	    }
	  }

	  if (x < (int)threshold_area.size() && x >= 0 && t + range[0] - 25 < (int)threshold_area.at(s).size() && t + range[0] - 25 >= 0){
	    for (int p = t + range[0] - 25; p < t + range[0]; p++){
	      sidebandt.at(sidebandt.size() - 1) += nfrw.adc_value.at(x).at(p);
	      st.push_back(nfrw.adc_value.at(x).at(p));
	    }
	  }
	  // if (x < (int)threshold_area.size() && x >= 0 && t + range[1] + 13 < (int)threshold_area.at(s).size() && t + range[1] + 13 >= 0){
	  //   for (int p = t + range[1] + 1; p <= t + range[1] + 13; p++){
	  //     sidebandf.at(sidebandf.size() - 1) += nfrw.adc_value.at(x).at(p);
	  //     sf.push_back(nfrw.adc_value.at(x).at(p));
	  //   }
	  // }

	  // if (x < (int)threshold_area.size() && x >= 0 && t - range[0] - 13 < (int)threshold_area.at(s).size() && t - range[0] - 13 >= 0){
	  //   for (int p = t - range[0] - 13; p < t - range[0]; p++){
	  //     sidebandt.at(sidebandt.size() - 1) += nfrw.adc_value.at(x).at(p);
	  //     st.push_back(nfrw.adc_value.at(x).at(p));
	  //   }
	  // }

	  // if (x < (int)threshold_area.size() && x >= 0 && t + i_parameter_y_window + 15 < (int)threshold_area.at(s).size() && t + i_parameter_y_window + 15 >= 0){
	  //   for (int p = t + i_parameter_y_window + 1; p <= t + i_parameter_y_window + 15; p++){
	  //     sidebandf.at(sidebandf.size() - 1) += nfrw.adc_value.at(x).at(p);
	  //     sf.push_back(nfrw.adc_value.at(x).at(p));
	  //   }
	  // }

	  // if (x < (int)threshold_area.size() && x >= 0 && t - i_parameter_y_window - 15 < (int)threshold_area.at(s).size() && t - i_parameter_y_window - 15 >= 0){
	  //   for (int p = t - i_parameter_y_window - 15; p < t - i_parameter_y_window; p++){
	  //     sidebandt.at(sidebandt.size() - 1) += nfrw.adc_value.at(x).at(p);
	  //     st.push_back(nfrw.adc_value.at(x).at(p));
	  //   }
	  // }
	  
	  sort(sf.begin(), sf.end());
	  sort(st.begin(), st.end());

	  sidebandf.at(sidebandf.size() - 1) /= 25.0;
	  sidebandt.at(sidebandt.size() - 1) /= 25.0;
	  // sidebandf.at(sidebandf.size() - 1) /= 12.0;
	  // sidebandt.at(sidebandt.size() - 1) /= 12.0;

	  for (int y = t - i_parameter_y_window; y <= t + i_parameter_y_window; y++){
	    if (x < (int)threshold_area.size() && x >= 0 && y < (int)threshold_area.at(s).size() && y >= 0){
	      if (y - t >= range[0] && y - t <= range[1]){
		asymmetric = 1;
	      }else{
		asymmetric = 0;
	      }

              // Find adjustment for the mode for sagging waveforms
	      lineval = bline(t + range[0] - 13, st[12], t + range[1] + 13, sf[12], y); //using median and +/- 25 bins

	      if (nfrw.adc_value.at(s).at(t) - lineval <= f_parameter_threshold){
		stuck_bit = true; // Not actually a stuck bit, but if this is tripped, the candidate is thrown out.
	      }                   // AKA, after mode adjustment, throw out candidates smaller than threshold

	      if (nfrw.detector_properties.DET == 1) //E-response Sag Fix for PDUNE
		IntWindow.push_back(asymmetric * (nfrw.adc_value.at(x).at(y) - lineval));

	      else if (nfrw.detector_properties.DET == 0 && (nfrw.data_properties.MCCX == 8 || !nfrw.data_properties.recoused))
		IntWindow.push_back(nfrw.adc_value.at(x).at(y));

	      else if (nfrw.detector_properties.DET == 0 && nfrw.data_properties.MCCX == 9 && nfrw.data_properties.recoused){
		IntWindow.push_back(nfrw.reco_product_value.at(x).at(y));
	      }

	      threshold_area.at(x).at(y) = 1;	    

              
              // stuck bit second pass
	      if (y > t - i_parameter_y_window && y < t + i_parameter_y_window){
		if (nfrw.adc_value[x][y] - nfrw.adc_value[x][y - 1] > 8 && 
		    nfrw.adc_value[x][y] - nfrw.adc_value[x][y + 1] > 8){
		  stuck_bit = true;
		}
	      }
	    }
	  }
	  wire++;
	}
      }

      // Clear stuck bit candidates
      if (stuck_bit && b_parameter_noise_study != 1) {
      	sidebandf.clear();
      	sidebandt.clear();
      	IntWindow.clear();	
      	continue;
      }
      
      ICharge.push_back(accumulate(IntWindow.begin(), IntWindow.end(), 0));

      c_info[8] = nfrw.candidate_info[i_coord / 2][8];
      TEinfo[1] = TEinfo[1] + 1.0;

      if (c_info[8] == 1)
	TEinfo[2] = TEinfo[2] + 1.0;
	
      c_info[0] = accumulate(IntWindow.begin(), IntWindow.end(), 0);
      //c_info[0] = charge;
      c_info[1] = c_info[0]; // Energy same as charge for now -- this can probably be removed. Use ROOT scripts for charge --> energy
      c_info[2] = nfrw.candidate_info[i_coord / 2][0];
      c_info[3] = nfrw.candidate_info[i_coord / 2][1];
      c_info[4] = nfrw.candidate_info[i_coord / 2][2];
      c_info[5] = nfrw.candidate_info[i_coord / 2][3];
      c_info[6] = nfrw.candidate_info[i_coord / 2][4];
      c_info[7] = nfrw.candidate_info[i_coord / 2][5];
      c_info[9] = nfrw.candidate_info[i_coord / 2][6];
      c_info[10] = nfrw.candidate_info[i_coord / 2][7];
      cind++;

      for (unsigned short l = 0; l < nfrw.collection_channel.size(); l++){
      	for (unsigned short p = 0; p < nfrw.collection_channel[l].size(); p++){
	 
      	  if ((double)nfrw.collection_channel[l][p] == s)
      	    c_info[11] = l;
 
      	}
      }

      c_info[11] = nfrw.candidate_info[i_coord / 2][9];
      c_info[12] = (double)stuck_bit;
      c_info[13] = nfrw.candidate_info[i_coord / 2][11];
      c_info[14] = nfrw.candidate_info[i_coord / 2][12];

      ctree.Fill();

      // // debug message
      // if(accumulate(IntWindow.begin(), IntWindow.end(), 0) > 3000)
      //   std::cout << "Abnormal charge near time " << t << " and wire " << s - 4800 << std::endl;
	
      IntWindow.clear();
      sidebandf.clear();
      sidebandt.clear();
    }

    std::cout << "finished sideband about to be done" << std::endl;
  }

  //return adc_value;
  return threshold_area;
  // return track_exclusion_map;
}


void rad_analysis::Confirm_Candidates(std::vector<short>& fcoor, rad_analysis::Waveforms& nfrw){

  //ctree already filled at this point
  //fcoor already filled here too

  //make confirm_candidates return a std::vector of candidate numbers and set a confirmed flag if sigselect integrates one.

  unsigned short umin, umax, vmin, vmax; 

  std::vector< std::vector<unsigned short> > matchedset; //collection of matched sets of wires (y, u, v)
  std::vector<double> listy;

  std::vector< std::vector< std::vector<unsigned int> > > tmp(6400, std::vector< std::vector<unsigned int> >()); //tmp is indexed by time tick
                                                                                        //time tick --> (channel, candidate)
  for (size_t i_coord = 0; i_coord < fcoor.size(); i_coord += 2){ 
    tmp[fcoor[i_coord + 1]].push_back(std::vector<unsigned int> ());
    tmp[fcoor[i_coord + 1]][tmp[fcoor[i_coord + 1]].size() - 1].push_back(fcoor[i_coord]); //channel
    tmp[fcoor[i_coord + 1]][tmp[fcoor[i_coord + 1]].size() - 1].push_back(i_coord / 2);    //candidate number
  } //filled tmp with all the candidate coordinates, turning 1d matrix into 3d matrix. 

  //loop through time
  for (unsigned short t = 0; t < nfrw.detector_properties.NTT; t++){
    for (unsigned int dc = 0; dc < tmp[t].size(); dc++){

      if (tmp[t][dc][0] > 4800){
	matchedset.push_back(std::vector<unsigned short>());
	matchedset[matchedset.size() - 1].push_back(tmp[t][dc][1]); //fill collection wire candiate number

  	umin = nfrw.plane_wire_intersections[tmp[t][dc][0]][1];
  	umax = nfrw.plane_wire_intersections[tmp[t][dc][0]][2];
  	vmin = nfrw.plane_wire_intersections[tmp[t][dc][0]][3];
  	vmax = nfrw.plane_wire_intersections[tmp[t][dc][0]][4];

	for (unsigned short t1 = t - 10; t1 < t + 10; t1++){ //look through some tolerance of time ticks for coincident signal

	  if (t1 >= 0 && t1 < nfrw.detector_properties.NTT){ //but make sure it doesn't exceed std::vector sizes

	    for (unsigned int dc1 = 0; dc1 < tmp[t1].size(); dc1++){

	      if ((tmp[t1][dc1][0] >= umin && tmp[t1][dc1][0] <= umax) || (tmp[t1][dc1][0] >= vmin && tmp[t1][dc1][0] <= vmax)){
	      //if (tmp[t1][dc1][0] >= umin && tmp[t1][dc1][0] <= umax){ //ignore V
		matchedset[matchedset.size() - 1].push_back(tmp[t1][dc1][1]); //fill induction wire candidate numbers
	      }

	     // else
	      // 	std::cout << tmp[t1][dc1][0] << "\t" << umin << "\t" << umax << "\t" << vmin << "\t" << vmax << std::endl;

	    }
	  }
	}
      }
    }
  }

  //std::cout << "filled std::vector of matching sets of wires" << std::endl;
  //std::cout << matchedset.size() << " matches" << std::endl;

  {
    short indx;

    int candex = 0;

    double z1, z2, y1, y2, z, test, y;

    for (unsigned int i = 0; i < matchedset.size(); i++){ //loop through matched candidates
      indx = 0;
      y = -999;

      if (matchedset[i].size() != (listy.size() + 1)){
	listy.clear();
	listy = std::vector<double> (matchedset[i].size() - 1, 0);
      }

      if (listy.size() > 1){

	z = nfrw.geometry_info[fcoor[matchedset[i][0] * 2]][5];

	for (unsigned int j = 1; j < listy.size(); j++){
	  z1 = nfrw.geometry_info[fcoor[matchedset[i][j] * 2]][5];
	  z2 = nfrw.geometry_info[fcoor[matchedset[i][j] * 2]][8];
	  y1 = nfrw.geometry_info[fcoor[matchedset[i][j] * 2]][4];
	  y2 = nfrw.geometry_info[fcoor[matchedset[i][j] * 2]][7];

	  listy[j - 1] = bline(z1, y1, z2, y2, z);
	}


	sort(listy.begin(), listy.end());

	for (unsigned int p = 0; p < listy.size() - 1; p++){

	  if (abs(listy[p] - listy[p + 1]) < 3){ //y tolerance
 	    indx++; 
	    y = (listy[p] + listy[p + 1]) / 2;
	  }

	  else {
	    indx = 0;
	  }	  
	}

	//delete elements that don't match with y
	{
	  unsigned int j = 1;
	  
	  while (j < matchedset[i].size()){
	    z1 = nfrw.geometry_info[fcoor[matchedset[i][j] * 2]][5];
	    z2 = nfrw.geometry_info[fcoor[matchedset[i][j] * 2]][8];
	    y1 = nfrw.geometry_info[fcoor[matchedset[i][j] * 2]][4];
	    y2 = nfrw.geometry_info[fcoor[matchedset[i][j] * 2]][7];

	    if (abs(bline(z1, y1, z2, y2, z) - y) > 3){ //y tolerance [cm]
	      matchedset[i].erase(matchedset[i].begin() + j);
	      j = 1;
	    }
	    
	    else {j++;}
	  }
	}
      }


      // //output listy
      // if (y != -1){std::cout << "accepted: " << std::endl;}
      // else {std::cout << "no y pos: " << std::endl;}
      // for (unsigned int p = 0; p < listy.size(); p++){
      // 	std::cout << listy[p];
	
      // 	if (p == listy.size() - 1)
      // 	  std::cout << std::endl;
      // 	else
      // 	  std::cout << ", ";
      // }


      if (y != -999 && matchedset[i].size() == 3){
	candex++;

	for (unsigned int j = 0; j < matchedset[i].size(); j++){
	  nfrw.candidate_info[matchedset[i][j]][10] = candex;
	  nfrw.candidate_info[matchedset[i][j]][11] = y;
	  nfrw.candidate_info[matchedset[i][j]][12] = z;
	}
      }
    
      //if (matchedset[i].size() == 3 && nfrw.channel_to_wire[fcoor[matchedset[i][1] * 2]][2] != nfrw.channel_to_wire[fcoor[matchedset[i][2] * 2]][2]){
      //if (matchedset[i].size() == 2 && nfrw.channel_to_wire[fcoor[matchedset[i][1] * 2]][2] == 0){ //use this when ignoring V plane

	// z1 = nfrw.geometry_info[fcoor[matchedset[i][2] * 2]][5];
	// z2 = nfrw.geometry_info[fcoor[matchedset[i][2] * 2]][8];
	// y1 = nfrw.geometry_info[fcoor[matchedset[i][2] * 2]][4];
	// y2 = nfrw.geometry_info[fcoor[matchedset[i][2] * 2]][7];
	// z = nfrw.geometry_info[fcoor[matchedset[i][0] * 2]][5];

	// test = bline(z1, y1, z2, y2, z);

	// z1 = nfrw.geometry_info[fcoor[matchedset[i][1] * 2]][5];
	// z2 = nfrw.geometry_info[fcoor[matchedset[i][1] * 2]][8];
	// y1 = nfrw.geometry_info[fcoor[matchedset[i][1] * 2]][4];
	// y2 = nfrw.geometry_info[fcoor[matchedset[i][1] * 2]][7];
	// z = nfrw.geometry_info[fcoor[matchedset[i][0] * 2]][5];

	// y = bline(z1, y1, z2, y2, z);	

	// if (abs(y - test) < 0.4){
	  
	//   candex++;

	//   for (unsigned int j = 0; j < matchedset[i].size(); j++){
	//     nfrw.candidate_info[matchedset[i][j]][10] = candex;
	//     nfrw.candidate_info[matchedset[i][j]][11] = y;
	//     nfrw.candidate_info[matchedset[i][j]][12] = z;
	//   }
	// }
        // }
    }
  }
}
