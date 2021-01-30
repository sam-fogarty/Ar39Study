#include <iostream>
#include <fstream>
#include <vector> 

using namespace std;

int main(){
  ifstream inf ("ChannelWireGeometry_v2.txt");

  bool outp;
  float dummy;
  double z, zi, zf;

  vector< vector<float> > geominfo (8256, vector<float>(9, 0));
  
  for (int i = 0; i < 8256; i++){
    for (char q = 0; q < 9; q++){
      inf >> geominfo[i][q];
    }
  }

  inf.close();


  for (short c = 0; c < 3456; c++){

    outp = false;

    z = geominfo[c + 4800][5];

    cout << (c + 4800) << "\t";
    //cout << z << "\t";    //debug


    for (short uv = 0; uv < 4800; uv++){

      zi = geominfo[uv][5];
      zf = geominfo[uv][8];

      if (!((zi < z) != (zf < z)) && outp == false){
	//ignore
      }
      else if (((zi < z) != (zf < z)) && outp == false){
	cout << uv << "\t";
	//cout << zi << "\t" << zf << "\t";   //debug
	outp = true;
      }
      else if (((zi < z) != (zf < z)) && outp == true){
	//ignore
      }
      else if (!((zi < z) != (zf < z)) && outp == true){
	cout << (uv - 1) << "\t"; 
	//cout << zi << "\t" << zf << "\t";   //debug
	outp = false;
      }
    }
    
    if (c == 3455)
      cout << 4799;

    cout << endl;
    //cout << "\t" << "in " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - st).count() << "ms" << endl;   //debug 

    
  }
  
  return 0; 
}
