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
//#include "TObject.h"

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
#include "counting.hh"
#include "Primaries.h"

#include "Ray.h"

class EarthModel; //class

void test();

string outputdir="outputs";

int main() {


    // below is replace by Settings class Initialize() member.
/*    
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
*/  

  Settings *settings1 = new Settings();


  cout<<"\n\tDefault values!"<<endl;
  cout<<"NNU : "<<settings1->NNU<<endl;
  cout<<"ICE_MODEL : "<<settings1->ICE_MODEL<<endl;
  cout<<"NOFZ : "<<settings1->NOFZ<<endl;
  cout<<"CONSTANTICETHICKNESS : "<<settings1->CONSTANTICETHICKNESS<<endl;
  cout<<"FIXEDELEVATION : "<<settings1->FIXEDELEVATION<<endl;
  cout<<"MOOREBAY : "<<settings1->MOOREBAY<<endl;
  cout<<"EXPONENT : "<<settings1->EXPONENT<<endl;
  cout<<"DETECTOR : "<<settings1->DETECTOR<<endl;

  string setupfile = "setup.txt";

  settings1->ReadFile(setupfile);

  cout<<"\n\tNew values!"<<endl;
  cout<<"NNU : "<<settings1->NNU<<endl;
  cout<<"ICE_MODEL : "<<settings1->ICE_MODEL<<endl;
  cout<<"NOFZ : "<<settings1->NOFZ<<endl;
  cout<<"CONSTANTICETHICKNESS : "<<settings1->CONSTANTICETHICKNESS<<endl;
  cout<<"FIXEDELEVATION : "<<settings1->FIXEDELEVATION<<endl;
  cout<<"MOOREBAY : "<<settings1->MOOREBAY<<endl;
  cout<<"EXPONENT : "<<settings1->EXPONENT<<endl;
  cout<<"DETECTOR : "<<settings1->DETECTOR<<endl;


//  IceModel *icemodel=new IceModel(ICE_MODEL + NOFZ*10,CONSTANTICETHICKNESS * 1000 + CONSTANTCRUST * 100 + FIXEDELEVATION * 10 + 0,MOOREBAY);// creates Antarctica ice model
  IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
  //IceModel inherits from EarthModel  

  cout<<endl;
  cout<<"Surface at (log:0, lat:0) : "<<icemodel->Surface(0., 0.)<<endl;
  cout<<"SurfaceAboveGeoid at (log:0, lat:0) : "<<icemodel->SurfaceAboveGeoid(0., 0.)<<endl;
  
  Detector *detector=new Detector(settings1->DETECTOR, icemodel); // builds antenna array, 0 for testbed
//  Detector *detector=new Detector(settings1->DETECTOR); // builds antenna array, 0 for testbed

  //  Trigger *trigger=new Trigger(detector); // builds the trigger  
//  Efficiencies *efficiencies=new Efficiencies(detector->getnRx(),outputdir); // keeps track of efficiencies at each stage of the simulation
  Efficiencies *efficiencies=new Efficiencies(100,outputdir); // keeps track of efficiencies at each stage of the simulation
  
  Spectra *spectra=new Spectra(settings1->EXPONENT); // gets library (or whatever) of neutrino spectra

  

    // 
    // test PickUnbiased in IceModel.
  Counting *count1 = new Counting();
  Primaries *primary1 = new Primaries();
  int whichray = 0; // for test
  Interaction *interaction1=new Interaction("nu",primary1,settings1,whichray,count1);




  TFile *AraFile=new TFile((outputdir+"/AraOut.root").c_str(),"RECREATE","ara");
  TTree *AraTree=new TTree("AraTree","AraTree");    // for single entry
  TTree *AraTree2=new TTree("AraTree2","AraTree2"); //for many entries

  AraTree->Branch("detector",&detector);
  AraTree->Branch("icemodel",&icemodel);
  AraTree->Branch("settings",&settings1);
  AraTree->Branch("spectra",&spectra);
  AraTree2->Branch("interaction",&interaction1);







//--------------------------------------------------
//   for (int inu = 0; inu<30; inu++) {
//       cout<<"\n";
//       int pickunbiased_out = icemodel->PickUnbiased(0, interaction1, icemodel);
//       cout<<"pickunbiased output : "<<pickunbiased_out<<endl;
//       cout<<"nnu : ";
//       interaction1->nnu.Print();
//       cout<<"posnu : ";
//       interaction1->posnu.Print();
//       cout<<"test random # :"<<gRandom->Rndm()<<endl;
//       cout<<"\n";
//       cout<<"\n";
//   }
//-------------------------------------------------- 






   for (int inu=0;inu<settings1->NNU;inu++) { // loop over neutrinos



       icemodel->PickUnbiased(inu, interaction1, icemodel); //pick posnu and nnu unbiased




    
//    detector->resetDetector(); // set all signals on antennas to zero


//    Event *event=new Event(spectra->GetNuEnergy()); // creates a neutrino interaction
    // contains info about neutrino itself and the interaction
    // input neutrino energy
   


//    detector->simulateDetector(event,efficiencies);// simulate detector response
    // including waveforms for analysis
    // this also includes a trigger simulation



    AraTree2->Fill();   //fill interaction every events


    // for test, print first two nnu X
    if (inu<2) {
        cout<<"interaction nnu : "<<interaction1->nnu.GetX()<<endl;
    }

  } // end loop over neutrinos
//--------------------------------------------------
//   cout<<"Total NNU : "<<settings1->NNU<<", PickUnbiased passed NNU : "<<nnu_pass<<endl;
//-------------------------------------------------- 
    
  
  AraTree->Fill();  // fill tree for one entry
  AraFile->Write();


 efficiencies->summarize(); // summarize the results in an output file  



// test
cout<<"station[0] x : "<<detector->stations[0].GetX()<<endl;
cout<<"string[0] x : "<<detector->stations[0].strings[0].GetX()<<endl;
cout<<"antenna[0] x : "<<detector->stations[0].strings[0].antennas[0].GetX()<<endl;
cout<<"antenna[0] Gain(700,5,0) : "<<detector->stations[0].strings[0].antennas[0].GetG(detector,700.,5.,0.)<<endl;
cout<<"GetGain(700,5,0,0) : "<<detector->GetGain(700.,5.,0.,0)<<endl;
cout<<"GetGain(10,5,0,0) : "<<detector->GetGain(10.,5.,0.,0)<<endl;

cout<<"params.number_of_stations : "<<detector->params.number_of_stations<<endl;
cout<<"params.station_spacing : "<<detector->params.station_spacing<<endl;

cout<<"Spectra energy values : "<<spectra->GetE_bin()<<endl;
for (int i=0;i<spectra->GetE_bin();i++) {
    cout<<"energy bin "<<i<<" : "<<spectra->energy[i]<<endl;
}

cout<<"IceModel R_EARTH = "<<icemodel->R_EARTH<<endl;


cout<<"Detector static const double freq_init : "<<detector->Getfreq_init()<<endl;






 delete icemodel;
 delete efficiencies;
 
 delete detector;
 delete settings1;
 delete count1;
 delete primary1;
 delete interaction1;

 delete spectra;

 test();


 return 0;
  
} //end main



void test() {

  cout << "test is " << 0 << "\n";
}
