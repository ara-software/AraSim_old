#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include "TTreeIndex.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h" 
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"

//#include <fftw3.h>

using namespace std;

#include "Tools.h"
#include "Constants.h"
#include "Vector.h"
#include "Position.h"
#include "EarthModel.h"
#include "IceModel.h"
#include "Efficiencies.h"
#include "Spectra.h"
#include "Event.h"
#include "Trigger.h"
#include "Detector.h"
#include "Settings.h"

#include "Ray.h"

class EarthModel; //class


string outputdir="outputs";

int main() {


//  Settings *settings1 = new Settings();

  //  Detector *detector=new Detector(settings1->DETECTOR); // builds antenna array, 0 for testbed
//  Detector *detector=0; // builds antenna array, 0 for testbed
  Detector *detector = 0; 
  cout<<"construct detector"<<endl;

  
  TFile *AraFile=new TFile((outputdir+"/AraOut.root").c_str());
  cout<<"AraFile"<<endl;
  TTree *AraTree=(TTree*)AraFile->Get("AraTree");
  cout<<"AraTree"<<endl;
  AraTree->SetBranchAddress("detector",&detector);
  cout<<"branch detector"<<endl;
  
  AraTree->GetEvent(0);
  cout<<"getevent"<<endl;
  cout << "I'm here.\n";
  cout << "freq_step_max is " << detector->GetFreqStepMax() << "\n";
  cout << "Vgain is " << detector->Vgain[0][10] << "\n";
  cout << "GetGain(700,5,0,0) is " << detector->GetGain(700., 5., 0., 0) << "\n";
  cout << "GetGain(10,5,0,0) is " << detector->GetGain(10., 5., 0., 0) << "\n";

  cout<<"station x is "<<detector->stations[0].GetX()<<endl;
  cout<<"string x is "<<detector->stations[0].strings[0].GetX()<<endl;
  cout<<"antenna x is "<<detector->stations[0].strings[0].antennas[0].GetX()<<endl;
  cout<<"antenna Gain(700,5,0) is "<<detector->stations[0].strings[0].antennas[0].GetG(detector,700.,5.,0.)<<endl;

  cout<<"params.number_of_stations : "<<detector->params.number_of_stations<<endl;
  cout<<"params.station_spacing : "<<detector->params.station_spacing<<endl;


}


