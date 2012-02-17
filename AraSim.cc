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
#include "signal.hh"
#include "secondaries.hh"

#include "Ray.h"
#include "RaySolver.h"
#include "Report.h"


class EarthModel; //class

void test();

string outputdir="outputs";

//--------------------------------------------------
// extern"C" {
//     void model_(int *ii);
// }
//-------------------------------------------------- 


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
  cout<<"end calling detector"<<endl;
//  Detector *detector=new Detector(settings1->DETECTOR); // builds antenna array, 0 for testbed

  //  Trigger *trigger=new Trigger(detector); // builds the trigger  
//  Efficiencies *efficiencies=new Efficiencies(detector->getnRx(),outputdir); // keeps track of efficiencies at each stage of the simulation
  Efficiencies *efficiencies=new Efficiencies(100,outputdir); // keeps track of efficiencies at each stage of the simulation
  cout<<"called Efficiencies"<<endl;
  
  Spectra *spectra=new Spectra(settings1->EXPONENT); // gets library (or whatever) of neutrino spectra
  cout<<"called Spectra"<<endl;

  Ray *ray = new Ray(); // construct Ray class
  cout<<"called Ray"<<endl;
  

    // 
    // test PickUnbiased in IceModel.
  Counting *count1 = new Counting();
  cout<<"called Counting"<<endl;

  Primaries *primary1 = new Primaries();
  cout<<"called Primaries"<<endl;

  int whichray = 0; // for test
//--------------------------------------------------
//   Interaction *interaction1=new Interaction("nu",primary1,settings1,whichray,count1);
//   cout<<"called Interaction1"<<endl;
//-------------------------------------------------- 

  Event *event = new Event();
  cout<<"called Event"<<endl;

  Report *report = new Report(detector);
  cout<<"called Evt"<<endl;



  TFile *AraFile=new TFile((outputdir+"/AraOut.root").c_str(),"RECREATE","ara");
  TTree *AraTree=new TTree("AraTree","AraTree");    // for single entry
  TTree *AraTree2=new TTree("AraTree2","AraTree2"); //for many entries
  cout<<"assing AraFile, AraTrees"<<endl;

  AraTree->Branch("detector",&detector);
  cout<<"branch detector"<<endl;
  AraTree->Branch("icemodel",&icemodel);
  cout<<"branch icemodel"<<endl;
  AraTree->Branch("settings",&settings1);
  cout<<"branch settings"<<endl;
  AraTree->Branch("spectra",&spectra);
  cout<<"branch spectra"<<endl;
  AraTree2->Branch("event",&event);
  cout<<"branch Evt"<<endl;
  AraTree2->Branch("report",&report);
  cout<<"branch report"<<endl;

  cout<<"finish tree assign"<<endl;


RaySolver *raysolver = new RaySolver;
cout<<"called RaySolver"<<endl;



cout<<"will call secondaries"<<endl;
Secondaries *sec1 = new Secondaries (settings1);
cout<<"will call signal"<<endl;
//Signal *signal = new Signal;
Signal *signal = new Signal (settings1);
signal->SetMedium(0);   // set medium as ice
cout<<"finish calling secondaries and signal"<<endl;




//--------------------------------------------------
// TH1F *hy=new TH1F("hy","hy",100,0.,1.); // histogram for inelasticity
//-------------------------------------------------- 




cout<<"begain looping events!!"<<endl;
   for (int inu=0;inu<settings1->NNU;inu++) { // loop over neutrinos


       event = new Event ( settings1, spectra, primary1, icemodel, detector, signal, sec1 );

       cout<<"inu : "<<inu<<endl;
       cout<<"event->pnu : "<<event->pnu<<endl;
       cout<<"posnu : ";
       event->Nu_Interaction[0].posnu.Print();
       cout<<"nnu : ";
       event->Nu_Interaction[0].nnu.Print();
       cout<<"event->n_interactions : "<<event->n_interactions<<endl;
       cout<<"nu_flavor : "<<event->nuflavor<<endl;
       cout<<"event->Nu_Interaction[0].vmmhz1m[0] : "<<event->Nu_Interaction[0].vmmhz1m[0]<<endl;
       cout<<"pickposnu : "<<event->Nu_Interaction[0].pickposnu<<endl;


       report->Connect_Interaction_Detector (event, detector, raysolver, signal, icemodel, settings1);

       
       for (int i=0;i<1;i++) {  // there are only 1 station for the test!!!
           for (int j=0; j<4;j++) { // 4 strings per station
               for (int k=0;k<4;k++) {  // 4 antennas per string

                   if ( event->Nu_Interaction[0].pickposnu && report->stations[i].strings[j].antennas[k].ray_sol_cnt ) {
                       //cout<<"Evt->pickposnu : "<<event->Nu_Interaction[0].pickposnu<<"\t report->...ray_sol_cnt : "<<report->stations[i].strings[j].antennas[k].ray_sol_cnt<<endl;
                       //event->Nu_Interaction[0].posnu.Print();
                       for (int l=0;l<report->stations[i].strings[j].antennas[k].ray_sol_cnt; l++) {    // loop for number of RaySolver solutions
                           cout<<"reflection : "<<report->stations[i].strings[j].antennas[k].reflection[l]<<endl;
                           cout<<"reflect_ang : "<<report->stations[i].strings[j].antennas[k].reflect_ang[l]<<endl;
                           cout<<"Pol_vector : ";
                           report->stations[i].strings[j].antennas[k].Pol_vector[l].Print();
                           for (int m=0;m<detector->GetFreqBin();m++) {
                               cout<<"evt "<<inu<<"; vmmhz for station["<<i<<"].string["<<j<<"].antenna["<<k<<"].vmmhz["<<l<<"]["<<m<<"] : "<<report->stations[i].strings[j].antennas[k].vmmhz[l][m]<<endl;
                           }
                    
                       }
                   }

               }
           }
       }



    AraTree2->Fill();   //fill interaction every events



  } // end loop over neutrinos


   cout<<" end loop"<<endl;
//--------------------------------------------------
//   cout<<"Total NNU : "<<settings1->NNU<<", PickUnbiased passed NNU : "<<nnu_pass<<endl;
//-------------------------------------------------- 
    
  
  AraTree->Fill();  // fill tree for one entry
  AraFile->Write();


 efficiencies->summarize(); // summarize the results in an output file  




//--------------------------------------------------
// int ii = 1;
// model_(&ii);
//-------------------------------------------------- 


//--------------------------------------------------
//  delete raysolver;
//-------------------------------------------------- 
 delete icemodel;
 delete efficiencies;
 delete ray;
 
 delete detector;
 delete settings1;
 delete count1;
 delete primary1;
 delete event;
 delete report;
//--------------------------------------------------
//  delete interaction1;
//-------------------------------------------------- 

 delete spectra;


 test();


 return 0;
  
} //end main



void test() {

  cout << "test is " << 0 << "\n";
}

