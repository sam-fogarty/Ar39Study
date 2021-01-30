////////////////////////////////////////////////////
// CERN ROOT Script to merge trees of grid files  //
// ROOT's hadd does this, but it also attempts to //
// merge histograms, which causes crashes.        //
// This is much more stable.                      //
// AFlesher                                       //
////////////////////////////////////////////////////

void GenerateRootSpectra(){
  const unsigned short data_run = 5841;
  vector<unsigned short> skip;
  unsigned short filenum;

  // vector<string> Nuclide = {"Ar39", "Co60", "Th232", "K40", "Kr85", "U238", "Rn222", "cosmic"};
  // vector<string> Nuclide = {"Ar39", "Co60", "Th232", "K40", "Kr85", "U238"};
  vector<string> Nuclide = {"39Ar", "60Co", "40K", "85Kr", "232Th", "238U"};      // name your directories in such a way that this 
                                                                                  // script can easily access them with the nuclide
  TChain *ctree = new TChain("ctree");
  TChain *etree = new TChain("etree");

  for (unsigned char r = 0; r < Nuclide.size(); r++){ // let r <= Nuclide.size() to include data, otherwise only look at simulation
    skip.clear();

    // setting up for differing file numbers and usuccessful runs -- This will need to be configured based on your data!
    if (r == Nuclide.size()){ // data
      filenum = 50;
      skip = {8, 43};
    }
    else if (Nuclide[r] == "Ar39" || Nuclide[r] == "39Ar"){
      filenum = 80;
      // skip = {25};
    }
    else if (Nuclide[r] == "K40" || Nuclide[r] == "40K"){filenum = 78;}
    else if (Nuclide[r] == "Co60" || Nuclide[r] == "60Co"){filenum = 62;}
    else if (Nuclide[r] == "U238" || Nuclide[r] == "238U"){filenum = 75;}
    else if (Nuclide[r] == "Rn222"){filenum = 79;}
    else if (Nuclide[r] == "Th232" || Nuclide[r] == "232Th"){
      // skip = {63};
      filenum = 29;
    }
    else if (Nuclide[r] == "cosmic"){filenum = 156;}
    else {filenum = 80;}


    for (int i = 1; i <= filenum; i++){
      if (binary_search(skip.begin(), skip.end(), i)) {continue;}

      if (r < Nuclide.size()){
	ctree->Add(TString::Format("/[path to grid analysis files]/%s/%i/output.root", Nuclide[r].c_str(), i));
	etree->Add(TString::Format("/[path to grid analysis files]/%s/%i/output.root", Nuclide[r].c_str(), i)); // list the file twice to get both event tree
                                                                                                                // and candidate tree
      }
      else{ // data
	ctree->Add(TString::Format("/[path to grid analysis files]/%i/%i/output.root", data_run, i));
	etree->Add(TString::Format("/[path to grid analysis files]/%i/%i/output.root", data_run, i));
      }       
    }


    if (Nuclide[r] == "39Ar"){Nuclide[r] = "Ar39";}
    else if (Nuclide[r] == "40K"){Nuclide[r] = "K40";}
    else if (Nuclide[r] == "60Co"){Nuclide[r] = "Co60";}
    else if (Nuclide[r] == "85Kr"){Nuclide[r] = "Kr85";}
    else if (Nuclide[r] == "232Th"){Nuclide[r] = "Th232";}
    else if (Nuclide[r] == "238U"){Nuclide[r] = "U238";}


    if (r == Nuclide.size()){
      ctree->Merge(TString::Format("./FitFiles/data_%i.root", data_run));
      etree->Merge(TString::Format("./FitFiles/data_%i_e.root", data_run));
    } // data
    else if (r < Nuclide.size()){
      ctree->Merge(TString::Format("./FitFiles/LRC_Radiological/%s.root", Nuclide[r].c_str()));
      etree->Merge(TString::Format("./FitFiles/LRC_Radiological/%s_e.root", Nuclide[r].c_str()));
    }

    ctree->Reset();
    etree->Reset();
  }
}
