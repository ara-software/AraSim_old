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

//Include output format to enable reading by analysis software AraRoot
#include "AraRootFormat/UsefulIcrrStationEvent.h"

class EarthModel; //class

void test();

string outputdir="outputs";

//--------------------------------------------------
// extern"C" {
//     void model_(int *ii);
// }
//-------------------------------------------------- 


//int main() {
int main(int argc, char **argv) {   // read setup.txt file



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

  //string setupfile = "setup.txt";
  string setupfile;
  string run_no;
  if (argc<2) { // no setup file input, use default
      setupfile = "setup.txt";
  }
  else if (argc == 2) { // read file!!
      setupfile = string( argv[1] );
  }
  else if (argc == 3) { // read file!!
      setupfile = string( argv[1] );
      run_no = string( argv[2] );
  }
  else { // no mode for argc > 2!
      cout<<"too many info! just use default setup.txt file!"<<endl;
      setupfile = "setup.txt";
  }

  settings1->ReadFile(setupfile);
  cout<<"Read "<<setupfile<<" file!"<<endl;


  cout<<"\n\tNew values!"<<endl;
  cout<<"NNU : "<<settings1->NNU<<endl;
  cout<<"ICE_MODEL : "<<settings1->ICE_MODEL<<endl;
  cout<<"NOFZ : "<<settings1->NOFZ<<endl;
  cout<<"CONSTANTICETHICKNESS : "<<settings1->CONSTANTICETHICKNESS<<endl;
  cout<<"FIXEDELEVATION : "<<settings1->FIXEDELEVATION<<endl;
  cout<<"MOOREBAY : "<<settings1->MOOREBAY<<endl;
  cout<<"EXPONENT : "<<settings1->EXPONENT<<endl;
  cout<<"DETECTOR : "<<settings1->DETECTOR<<endl;
  cout<<"POSNU_RADIUS : "<<settings1->POSNU_RADIUS<<endl;



  // set gRandom as TRandom3 when settings1->RANDOM_MODE = 1
  if (settings1->RANDOM_MODE == 1) {

    // test TRandom3
    TRandom3 *test_randm3 = new TRandom3 (0);
    gRandom = test_randm3;
  }
    //cout<<"first random from TRandom3 : "<<test_randm3->Rndm()<<"\n";
    cout<<"first random : "<<gRandom->Rndm()<<"\n";







//  IceModel *icemodel=new IceModel(ICE_MODEL + NOFZ*10,CONSTANTICETHICKNESS * 1000 + CONSTANTCRUST * 100 + FIXEDELEVATION * 10 + 0,MOOREBAY);// creates Antarctica ice model
  IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
  //IceModel inherits from EarthModel  

  cout<<endl;
  cout<<"Surface at (log:0, lat:0) : "<<icemodel->Surface(0., 0.)<<endl;
  cout<<"SurfaceAboveGeoid at (log:0, lat:0) : "<<icemodel->SurfaceAboveGeoid(0., 0.)<<endl;
  
  Detector *detector=new Detector(settings1, icemodel); // builds antenna array, 0 for testbed
  //Detector *detector=new Detector(settings1->DETECTOR, icemodel); // builds antenna array, 0 for testbed
  cout<<"end calling detector"<<endl;
//  Detector *detector=new Detector(settings1->DETECTOR); // builds antenna array, 0 for testbed

  Trigger *trigger=new Trigger(detector, settings1); // builds the trigger  
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

  //Report *report = new Report(detector, settings1);
  Report *report = new Report();
  cout<<"called Evt"<<endl;


  TFile *AraFile;
   if (argc == 3) {
        AraFile=new TFile((outputdir+"/AraOut."+setupfile+".run"+run_no+".root").c_str(),"RECREATE","ara");
   }
   else {
        AraFile=new TFile((outputdir+"/AraOut.root").c_str(),"RECREATE","ara");
   }

  TTree *AraTree=new TTree("AraTree","AraTree");    // for single entry
  TTree *AraTree2=new TTree("AraTree2","AraTree2"); //for many entries
  cout<<"assing AraFile, AraTrees"<<endl;

  AraTree->Branch("detector",&detector);
  cout<<"branch detector"<<endl;
  AraTree->Branch("icemodel",&icemodel);
  cout<<"branch icemodel"<<endl;
  AraTree->Branch("trigger",&trigger);
  cout<<"branch trigger"<<endl;
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

    cout << "Make output file that is readable by AraRoot" << endl;
    UsefulIcrrStationEvent *theEvent = 0;
    TTree *eventTree;
    eventTree = new TTree("eventTree","Tree of ARA Events");
    eventTree->Branch("UsefulIcrrStationEvent","UsefulIcrrStationEvent",&theEvent);


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

// before start looping events set noise values (this case, thermal)
trigger->SetMeanRmsDiode(settings1, detector, report);
// now in Trigger class, there will be meandiode, rmsdiode values for noise (we need this for trigger later)


double max_dt = 0.; // max arrival time difference

int Total_Global_Pass = 0;  // total global trigger passed number 
double Total_Weight = 0.;

/*

 TCanvas *cFull_window = new TCanvas("cFull_window","A Simple Graph Example",200,10,10000,11200);
 cFull_window->Divide(1,16);

 TGraph *g_Full_window;

TGraph *G_V_threshold_diode;
G_V_threshold_diode = new TGraph(2, threshold_x, threshold_y);

 TCanvas *cFull_window_V = new TCanvas("cFull_window_V","A Simple Graph Example",200,10,3200,2400);
 cFull_window_V->Divide(4,4);

 TGraph *g_Full_window_V;
 */







 double x_V[settings1->NFOUR/2];
 double y_V[settings1->NFOUR/2];



 double xbin[settings1->DATA_BIN_SIZE];
 for (int i=0; i<settings1->DATA_BIN_SIZE; i++) {
     xbin[i] = i;
 }

double threshold_y[2];
double threshold_x[2];

threshold_x[0] = 0.;
threshold_x[1] = (double)settings1->DATA_BIN_SIZE-1.;
threshold_y[0] = (trigger->rmsdiode) * (trigger->powerthreshold);
threshold_y[1] = (trigger->rmsdiode) * (trigger->powerthreshold);





cout<<"powerthreshold : "<<trigger->powerthreshold<<endl;


int check_station_DC;

cout<<"begain looping events!!"<<endl;
   for (int inu=0;inu<settings1->NNU;inu++) { // loop over neutrinos

       check_station_DC = 0;

       std::cerr<<"*";

       event = new Event ( settings1, spectra, primary1, icemodel, detector, signal, sec1 );


  
       report = new Report(detector, settings1);

//--------------------------------------------------
//        cout<<"inu : "<<inu<<endl;
//        cout<<"event->pnu : "<<event->pnu<<endl;
//        cout<<"posnu : ";
//        event->Nu_Interaction[0].posnu.Print();
//        cout<<"nnu : ";
//        event->Nu_Interaction[0].nnu.Print();
//        cout<<"event->n_interactions : "<<event->n_interactions<<endl;
//        cout<<"nu_flavor : "<<event->nuflavor<<endl;
//        cout<<"event->Nu_Interaction[0].vmmhz1m[0] : "<<event->Nu_Interaction[0].vmmhz1m[0]<<endl;
//        cout<<"pickposnu : "<<event->Nu_Interaction[0].pickposnu<<endl;
//-------------------------------------------------- 

       // connect Interaction class (nu interaction with ice) and Detector class (detector properties and layout)
       // save signal, noise at each antennas to Report class
       report->Connect_Interaction_Detector (event, detector, raysolver, signal, icemodel, settings1, trigger);

//       theEvent = new UsefulIcrrStationEvent();
       theEvent = &report->theUsefulEvent;
       eventTree->Fill();
       theEvent = NULL;

       /*

       
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








       cout<<"r_in : ";
       event->Nu_Interaction[0].r_in.Print();
       cout<<"nuexit : ";
       event->Nu_Interaction[0].nuexit.Print();
       cout<<"r_enterice : ";
       event->Nu_Interaction[0].r_enterice.Print();
       cout<<"nuexitice : ";
       event->Nu_Interaction[0].nuexitice.Print();
       cout<<"weight : "<<event->Nu_Interaction[0].weight<<endl;


       cout<<"In trigger, NFOUR : "<<trigger->NFOUR<<"  TIMESTEP : "<<trigger->TIMESTEP<<"  maxt_diode : "<<trigger->maxt_diode<<endl;
       cout<<"meandiode : "<<trigger->meandiode<<"  rmsdiode : "<<trigger->rmsdiode<<endl;
       cout<<"powerthreshold : "<<trigger->powerthreshold<<endl;

       cout<<"\nmaxt_diode : "<<detector->maxt_diode<<"  idelaybefore_peak : "<<detector->idelaybeforepeak<<"  iwindow : "<<detector->iwindow<<endl;


       */




       for (int i=0; i<detector->params.number_of_stations; i++) {
           if (max_dt < report->stations[i].max_arrival_time - report->stations[i].min_arrival_time) max_dt = report->stations[i].max_arrival_time - report->stations[i].min_arrival_time;
           // check the total global trigger passed
           if (report->stations[i].Global_Pass) {
               cout<<"\nGlobal_Pass : "<<report->stations[i].Global_Pass<<" evt : "<<inu<<" added weight : "<<event->Nu_Interaction[0].weight<<"\n"<<endl;

               if ( check_station_DC == 0) { // count trigger pass only once per event

               Total_Global_Pass ++;
               Total_Weight += event->Nu_Interaction[0].weight;

               // test increment weight
               count1->incrementEventsFound( event->Nu_Interaction[0].weight, event );


               }

               check_station_DC++;

               
               /*
               // make plots for all channels
               for (int string=0; string<detector->params.number_of_strings_station; string++) {
                   for (int antenna=0; antenna<detector->params.number_of_antennas_string; antenna++) {

                       //cout<<"plot cd "<<4*string + antenna<<endl;
                       cFull_window->cd( 4*string + antenna + 1 );
                       //g_Full_window = new TGraph( settings1->DATA_BIN_SIZE , xbin, report->Full_window[4*string + antenna] );
                       g_Full_window = new TGraph( settings1->DATA_BIN_SIZE , xbin, trigger->Full_window[4*string + antenna] );
                       g_Full_window->Draw("al");

                       G_V_threshold_diode -> SetLineColor(kRed);
                       G_V_threshold_diode -> Draw("l");




                       cFull_window_V->cd( 4*string + antenna + 1 );
                       for (int p=0; p<settings1->NFOUR/2; p++) {
                           x_V[p] = report->stations[i].strings[string].antennas[antenna].time[p];
                           y_V[p] = report->stations[i].strings[string].antennas[antenna].V_mimic[p];
                       }
                       //cout<<"x_V[0] : "<<x_V[0]<<endl;
                       g_Full_window_V = new TGraph( settings1->NFOUR/2 , x_V, y_V );
                       //g_Full_window_V->GetXaxis()->SetLimits(x_V[0],x_V[settings1->NFOUR/2-1]);
                       g_Full_window_V->Draw("al");




                   }
               }
               */


           }

       }


       // test memory
       //report->stations.clear();



           
       AraTree2->Fill();   //fill interaction every events



//--------------------------------------------------
// 
//  for (int i=0; i<4; i++) { // for strings
//      for (int j=0; j<4; j++) { // for antennas
//          if (report->stations[0].strings[i].antennas[j].ray_sol_cnt == 0 && report->stations[0].Total_ray_sol) { // if there's no raysol (should be only noise)
//              cout<<"evt : "<<inu<<" ch : "<<i*4+j<<endl;
//              //g_Full_window = new TGraph( settings1->DATA_BIN_SIZE , xbin, report->Full_window[i*4+j] );
//              g_Full_window = new TGraph( settings1->DATA_BIN_SIZE , xbin, report->Full_window[0] );
//              cout<<"ray_sol_cnt for 0.0 (plotted?) : "<<report->stations[0].strings[0].antennas[0].ray_sol_cnt<<endl;
//              i = 10;
//              j = 10;
//          }
//      }
//  }
//-------------------------------------------------- 

 //cout<<"evt "<<inu<<endl;


 delete event;
 delete report;
       delete theEvent;


  } // end loop over neutrinos



                       
   //cFull_window_V->Print("test_V_mimic.pdf");





   ofstream weight_file;
   //weight_file.open(("./weight_output/weight_"+setupfile).c_str());
   if (argc == 3) {
        weight_file.open(("./weight_output/weight_"+setupfile+".run"+run_no).c_str());
   }
   else {
        weight_file.open(("./weight_output/weight_"+setupfile).c_str());
   }


   cout<<" end loop"<<endl;
   cout<<"Total_Global_Pass : "<<Total_Global_Pass<<endl;
   cout<<"Total_Weight : "<<Total_Weight<<endl;
   weight_file << "Total_Weight="<<Total_Weight<<endl;

   cout<<"weight bin values : ";
   for (int i=0; i<count1->NBINS-1; i++) {
       cout<<count1->eventsfound_binned[i]<<", ";
       weight_file << count1->eventsfound_binned[i]<<" ";
   }
       cout<<count1->eventsfound_binned[count1->NBINS-1];
       weight_file << count1->eventsfound_binned[count1->NBINS-1]<<"\n";
       weight_file.close();
   cout<<"\n\n";

   double IceVolume;
   IceVolume = PI * (settings1->POSNU_RADIUS) * (settings1->POSNU_RADIUS) * icemodel->IceThickness( detector->stations[0] );
   cout<<"IceVolume : "<<IceVolume<<endl;

   double Veff_test;

   // error bar for weight
   double error_plus = 0;
   double error_minus = 0;
   Counting::findErrorOnSumWeights( count1->eventsfound_binned, error_plus, error_minus );

   Veff_test = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * Total_Weight / (double)(settings1->NNU);

   // account all factors to error
   error_plus = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * error_plus / (double)(settings1->NNU);
   error_minus = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * error_minus / (double)(settings1->NNU);

   cout<<"test Veff : "<<Veff_test<<" m3sr, "<<Veff_test*1.E-9<<" km3sr"<<endl;
   cout<<"And Veff error plus : "<<error_plus*1.E-9<<" and error minus : "<<error_minus*1.E-9<<endl;
//--------------------------------------------------
//   cout<<"Total NNU : "<<settings1->NNU<<", PickUnbiased passed NNU : "<<nnu_pass<<endl;
//-------------------------------------------------- 
    


   // remove noisewaveform info if DATA_SAVE_MODE == 2
   if (settings1->DATA_SAVE_MODE == 2) {
       trigger->v_noise_timedomain.clear();
       trigger->v_noise_timedomain_diode.clear();
   }
  
  AraTree->Fill();  // fill tree for one entry
  AraFile->Write();
  AraFile->Close();


 efficiencies->summarize(); // summarize the results in an output file  


 double freq[detector->GetFreqBin()], Filter[detector->GetFreqBin()];
 double Filter_E[detector->GetFreqBin()];

 for (int i=0; i<detector->GetFreqBin(); i++) {
     freq[i] = detector->GetFreq(i);    // in Hz
     Filter[i] = detector->GetFilterGain(i);    // in dB
     Filter_E[i] = pow(10., (detector->GetFilterGain(i))/20.);
 }




 cout<<"max_dt : "<<max_dt<<endl;


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
 delete trigger;
 delete spectra;
 delete sec1;
 delete signal;


 test();


 return 0;
  
} //end main



void test() {

  cout << "test is " << 0 << "\n";
}

