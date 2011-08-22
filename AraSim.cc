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


#include "Ray.h"

class EarthModel; //class

void test();

string outputdir="outputs";

int main() {
  
  const int NNU=100;

  // NEED TO FIGURE OUT A GOOD WAY TO READ THIS IN AND STORE THEM.
  // INPUT FILE AGAIN?  SOMETHING ELSE?
  //These were moved here from IceModel under the new compilation scheme
  int ICE_MODEL=0; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP.
  int NOFZ=1; // 1=depth dependent index of refraction,0=off
  int CONSTANTCRUST=0; // set crust density and thickness to constant values.
  int CONSTANTICETHICKNESS=0; // set ice thickness to constant value
  int FIXEDELEVATION=0; // fix the elevation to the thickness of ice.
  int MOOREBAY=0; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data
  double EXPONENT=19.; // 10^19 eV neutrinos only
  

  IceModel *icemodel=new IceModel(ICE_MODEL + NOFZ*10,CONSTANTICETHICKNESS * 1000 + CONSTANTCRUST * 100 + FIXEDELEVATION * 10 + 0,MOOREBAY);// creates Antarctica ice model
  //IceModel inherits from EarthModel  
  Detector *detector=new Detector(); // builds antenna array
  //  Trigger *trigger=new Trigger(detector); // builds the trigger  
  Efficiencies *efficiencies=new Efficiencies(detector->getnRx(),outputdir); // keeps track of efficiencies at each stage of the simulation
  Spectra *spectra=new Spectra(EXPONENT); // gets library (or whatever) of neutrino spectra

  
  TFile *AraFile=new TFile((outputdir+"/AraOut.root").c_str(),"RECREATE","ara");
  TTree *AraTree=new TTree("AraTree","AraTree");
  AraTree->Branch("detector",&detector);


  for (int inu=0;inu<NNU;inu++) { // loop over neutrinos

    
    detector->resetDetector(); // set all signals on antennas to zero


    Event *event=new Event(spectra->GetNuEnergy()); // creates a neutrino interaction
    // contains info about neutrino itself and the interaction
    // input neutrino energy
   


    detector->simulateDetector(event,efficiencies);// simulate detector response
    // including waveforms for analysis
    // this also includes a trigger simulation



    AraTree->Fill();

    //delete event;
  }
  AraTree->Write();
  AraFile->Close();
 efficiencies->summarize(); // summarize the results in an output file  


 delete icemodel;
 delete efficiencies;
 
 delete detector;

 

 return 0;
  
} //end main



void test() {

  cout << "test is " << 0 << "\n";
}
