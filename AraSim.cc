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

class EarthModel; //class

void test();

string outputdir="outputs";

extern"C" {
    void model_(int *ii);
}


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
  Interaction *interaction1=new Interaction("nu",primary1,settings1,whichray,count1);
  cout<<"called Interaction1"<<endl;

  Interaction *Evt = new Interaction();
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
  AraTree2->Branch("interaction",&Evt);
  cout<<"branch Evt"<<endl;
//--------------------------------------------------
//   AraTree2->Branch("interaction",&interaction1);
//-------------------------------------------------- 

  cout<<"finish tree assign"<<endl;


RaySolver *raysolver = new RaySolver;
cout<<"called RaySolver"<<endl;
//--------------------------------------------------
// raysolver->Solve_Ray(interaction1->posnu,detector->stations[0].strings[0].antennas[0]);
//-------------------------------------------------- 



cout<<"will call secondaries"<<endl;
Secondaries *sec1 = new Secondaries;
cout<<"will call signal"<<endl;
Signal *signal = new Signal;
signal->SetMedium(0);   // set medium as ice
cout<<"finish calling secondaries and signal"<<endl;

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




int currentint = 0;
int nu_nubar;   // test. output from Gety
string taudecay = "test_taudecay";
TH1F *hy=new TH1F("hy","hy",100,0.,1.); // histogram for inelasticity
double emfrac, hadfrac; // fraction of em shower and had shower
int n_interactions=1; // number of interactions for each events




// below variables for freq domain simulation (SIMULATION_MODE = 0)
           vector <double> vmmhz1m; // signal strength at 1m, certain freq.
           vector <double> vmmhz1m_em;  // same above, em shower part
           vector <double> d_theta_em;  // max spread angle from cherenkov angle (em)
           vector <double> d_theta_had; // same above (had)
           double vmmhz1m_tmp;  // as TaperVmMhz has input vmmhz1m as an output
           double viewangle;   // view angle
           Position launch_vector;  // launch vector from ray source
           Position R1; // vector for detector
           Position R2; // vector for ray source

       vector < vector <double> > ray_output;
// end of freq domain simulation variables







cout<<"begain looping events!!"<<endl;
   for (int inu=0;inu<settings1->NNU;inu++) { // loop over neutrinos



//       Evt = new Interaction(interaction_mode, icemodel, detector);
       Evt = new Interaction(icemodel, detector, settings1, whichray, count1, primary1, spectra); // it will set posnu, nnu, flavor, current, and pnu and store in Evt (Interaction class).
       cout<<"set new Interaction as Evt"<<endl;
       cout<<"pnu : "<<Evt->pnu<<endl;

       double elast_y = primary1->Gety(settings1, Evt->pnu, nu_nubar, currentint);  // set inelasticity
       cout<<"set inelasticity : "<<elast_y<<endl;

       sec1->GetEMFrac( settings1, Evt, taudecay, elast_y, hy, inu, emfrac, hadfrac, n_interactions);   // set em, had frac values.
       cout<<"set emfrac : "<<emfrac<<" hadfrac : "<<hadfrac<<endl;


       if (settings1->SIMULATION_MODE == 0) { // freq domain simulation (old mode)

           // set vmmhz1m (which is generally used for all detector antennas)
           // vmmhz1m is calculated for 1m, cherenkov angle
           //
           for (int i=0; i<detector->GetFreqBin(); i++) {   // for detector freq bin numbers

               if (inu==0) {
                   d_theta_em.push_back(0); // prepare d_theta_em and d_theta_had for GetSpread
                   d_theta_had.push_back(0);
                   vmmhz1m.push_back(0);
                   vmmhz1m_em.push_back(0);
               }

               signal->GetSpread(Evt->pnu, emfrac, hadfrac, detector->GetFreq(i), d_theta_em[i], d_theta_had[i]);   // get max spread angle and save at d_theta_em[i] and d_theta_had[i]
               cout<<"GetSpread, theta_em : "<<d_theta_em[i]<<" theta_had : "<<d_theta_had[i]<<endl;

               vmmhz1m[i] = signal->GetVmMHz1m( Evt->pnu, detector->GetFreq(i) );   // get VmMHz at 1m at cherenkov angle at GetFreq(i)
               cout<<"GetVmMHZ1m : "<<vmmhz1m[i]<<endl;
           }




       int ray_sol_cnt; // counting number of solutions from Solve_Ray

           
           for (int i = 0; i< detector->params.number_of_stations; i++) {

               for (int j=0; j< detector->params.number_of_strings_station; j++) {

                   for (int k=0; k< detector->params.number_of_antennas_string; k++) {

                       // run ray solver, see if solution exist
                       // if not, skip (set something like Sol_No = 0;
                       // if solution exist, calculate view angle and calculate TaperVmMHz

                       if (Evt->pickposnu) {    // if posnu is selected inside the antarctic ic:"<<viewangle<<" th_em:"<<d_theta_em[l]<<" th_had:"<<d_theta_had[l]<<" emfrac:"<<emfrac<<" hadfrac:"<<hadfrac<<" vmmhz1m:"<<vmmhz1m[l]<<endl;e
                           
                           raysolver->Solve_Ray(Evt, detector->stations[i].strings[j].antennas[k], icemodel, ray_output);   // solve ray between source and antenna
                           
//                           cout<<"solution_toggle : "<<raysolver->solution_toggle<<endl;   
                           
                           if (raysolver->solution_toggle && Evt->ray_solver_toggle) {  // if there are solution from raysolver
                               ray_sol_cnt = 0;
//                               cout<<"ray_output size : "<<ray_output[0].size()<<endl;
                               
                               while ( ray_sol_cnt < ray_output[0].size() ) {   // for number of soultions (could be 1 or 2)
//                                   cout<<"Path length : "<<ray_output[0][ray_sol_cnt]<<"\tView angle : "<<ray_output[1][ray_sol_cnt]<<"\tReceipt angle : "<<ray_output[2][ray_sol_cnt]<<endl;
                                   ray_sol_cnt++;

                                   R1 = detector->stations[i].strings[j].antennas[k];
                                   R2 = Evt->posnu;
                                   viewangle = PI/2. - ray_output[1][ray_sol_cnt];

                                   launch_vector = (R1.Cross( R1.Cross(R2) )).Rotate(viewangle, R1.Cross(R2));
                                   viewangle = launch_vector.Angle(Evt->nnu);

                                   for (int l=0; l<detector->GetFreqBin(); l++) {   // for detector freq bin numbers


//                                       cout<<"TaperVmMHz inputs VA:"<<viewangle<<" th_em:"<<d_theta_em[l]<<" th_had:"<<d_theta_had[l]<<" emfrac:"<<emfrac<<" hadfrac:"<<hadfrac<<" vmmhz1m:"<<vmmhz1m[l]<<endl;

                                       vmmhz1m_tmp = vmmhz1m[l];

                                       signal->TaperVmMHz( viewangle, d_theta_em[l], d_theta_had[l], emfrac, hadfrac, vmmhz1m_tmp, vmmhz1m_em[l]);
//                                       cout<<"TaperVmMHz (1m at view angle) at "<<l<<"th bin : "<<vmmhz1m_tmp<<endl;

                                   }// end for freq bin

                               
                               }// end while number of solutions
                           
                           }// end if solution exist
                       
                       }// end if posnu selected

                       

                   }// for number_of_antennas_string

               }// for number_of_strings_station

           }// for number_of_stations




       }// if SIMULATION_MODE = 0 (freq domain old method)








//--------------------------------------------------
//        vector < vector <double> > ray_output;
// 
//        int ray_sol_cnt;
// 
//        if (Evt->pickposnu) {
//            cout<<"posnu lat : "<<Evt->posnu.Lat()<<endl;
//            Evt->posnu.Print();
// //           raysolver->Solve_Ray(Evt->posnu, detector->stations[0].strings[0].antennas[0], icemodel, ray_output);
//            raysolver->Solve_Ray(Evt, detector->stations[0].strings[0].antennas[0], icemodel, ray_output);
//            cout<<"solution_toggle : "<<raysolver->solution_toggle<<endl;
//            if (raysolver->solution_toggle && Evt->ray_solver_toggle) {
//                ray_sol_cnt = 0;
//                cout<<"ray_output size : "<<ray_output[0].size()<<endl;
//                while ( ray_sol_cnt < ray_output[0].size() ) {
//                    cout<<"Path length : "<<ray_output[0][ray_sol_cnt]<<"\tView angle : "<<ray_output[1][ray_sol_cnt]<<"\tReceipt angle : "<<ray_output[2][ray_sol_cnt]<<endl;
//                    ray_sol_cnt++;
//                }
//            }
//        }
//-------------------------------------------------- 









    
//    detector->resetDetector(); // set all signals on antennas to zero


//    detector->simulateDetector(event,efficiencies);// simulate detector response
    // including waveforms for analysis
    // this also includes a trigger simulation



    AraTree2->Fill();   //fill interaction every events



  } // end loop over neutrinos
//--------------------------------------------------
//   cout<<"Total NNU : "<<settings1->NNU<<", PickUnbiased passed NNU : "<<nnu_pass<<endl;
//-------------------------------------------------- 
    
  
  AraTree->Fill();  // fill tree for one entry
  AraFile->Write();


 efficiencies->summarize(); // summarize the results in an output file  




//--------------------------------------------------
// raysolver->test1();
//-------------------------------------------------- 

//--------------------------------------------------
// int test1;
// int ari = 3;
// char* arr[ari];
// arr[1] = "0";
// test1 = sprintf(arr[2], "%f", 100.01);
// 
// cout<<"arr[1] : "<<arr[1]<<endl;
// cout<<"arr[2] : "<<arr[2]<<endl;
// 
//-------------------------------------------------- 

//--------------------------------------------------
// int ii = 1;
// model_(&ii);
//-------------------------------------------------- 


 delete raysolver;
 delete icemodel;
 delete efficiencies;
 delete ray;
 
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

