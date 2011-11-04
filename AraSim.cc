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
  
  cout<<"Surface at (log:0, lat:90) : "<<icemodel->Surface(0., 90.)<<endl;
  cout<<"SurfaceAboveGeoid at (log:0, lat:90) : "<<icemodel->SurfaceAboveGeoid(0., 90.)<<endl;

//  Detector *detector=new Detector(); // builds antenna array
  Detector *detector=new Detector(settings1->DETECTOR, icemodel); // builds antenna array, 0 for testbed
//  Detector *detector=new Detector(settings1->DETECTOR); // builds antenna array, 0 for testbed
//  Detector *detector=new Detector(1); // builds antenna array, 0 for testbed


  cout<<"test3"<<endl;

    double gain_return;

    //this is for Detector 1 or 2!!
//  gain_return = detector->stations[0].strings[0].antennas[0].GetG(detector,1067.0, 20, 50);
//  cout<<"testG_new : "<<*testG_new<<endl;
//    gain_return = *test2;




    if (settings1->DETECTOR == 0) {
        cout<<"Test strings[0] postion x : "<<detector->strings[0].GetX()<<endl;
        cout<<"Test strings[0].antennas[0] postion x : "<<detector->strings[0].antennas[0].GetX()<<endl;
    }
    else {
        cout<<"Test stations[0] postion x : "<<detector->stations[0].GetX()<<endl;
        cout<<"Test stations[0].strings[0] postion x : "<<detector->stations[0].strings[0].GetX()<<endl;
        cout<<"Test stations[0].strings[0].antennas[0] postion x : "<<detector->stations[0].strings[0].antennas[0].GetX()<<endl;
        cout<<"Test stations[0].strings[0].antennas[0] r, theta, phi : "<<detector->stations[0].strings[0].antennas[0].R()<<" "<<detector->stations[0].strings[0].antennas[0].Theta()<<" "<<detector->stations[0].strings[0].antennas[0].Phi()<<endl;
        cout<<"Test stations[0].strings[0].antennas[1] r, theta, phi : "<<detector->stations[0].strings[0].antennas[1].R()<<" "<<detector->stations[0].strings[0].antennas[1].Theta()<<" "<<detector->stations[0].strings[0].antennas[1].Phi()<<endl;
        cout<<"Test stations[0].strings[0].antennas[2] r, theta, phi : "<<detector->stations[0].strings[0].antennas[2].R()<<" "<<detector->stations[0].strings[0].antennas[2].Theta()<<" "<<detector->stations[0].strings[0].antennas[2].Phi()<<endl;
        cout<<"Test stations[0].strings[0].antennas[3] r, theta, phi : "<<detector->stations[0].strings[0].antennas[3].R()<<" "<<detector->stations[0].strings[0].antennas[3].Theta()<<" "<<detector->stations[0].strings[0].antennas[3].Phi()<<endl;
    }


    


    cout<<"\n";
    cout<<"test gain_return : "<<gain_return<<endl;
    cout<<"\n";


//    gain_return = detector->GetGain(711, 27, 357, 0);

//    cout<<"Gain for vpol @ 711MHz, theta 27, phi 357 is : "<<gain_return<<endl;

//    gain_return = detector->GetGain(231, 151, 73, 1);

//    cout<<"Gain for hpol @ 231MHz, theta 151, phi 73 is : "<<gain_return<<endl;












  //  Trigger *trigger=new Trigger(detector); // builds the trigger  
//  Efficiencies *efficiencies=new Efficiencies(detector->getnRx(),outputdir); // keeps track of efficiencies at each stage of the simulation
  Efficiencies *efficiencies=new Efficiencies(100,outputdir); // keeps track of efficiencies at each stage of the simulation
//  Spectra *spectra=new Spectra(EXPONENT); // gets library (or whatever) of neutrino spectra
  
  Spectra *spectra=new Spectra(settings1->EXPONENT); // gets library (or whatever) of neutrino spectra

  


//  for (int inu=0;inu<NNU;inu++) { // loop over neutrinos




    // 
    // test PickUnbiased in IceModel.
  Counting *count1 = new Counting();
  Primaries *primary1 = new Primaries();
  int whichray = 0; // for test
  Interaction *interaction1=new Interaction("nu",primary1,settings1,whichray,count1);




  TFile *AraFile=new TFile((outputdir+"/AraOut.root").c_str(),"RECREATE","ara");
  TTree *AraTree=new TTree("AraTree","AraTree");
//  TTree *AraTree2 = new TTree("AraTree2", "AraTree2");

//  Detector *d_test = detector;

//  AraTree->Branch("detector_test",&d_test);
  AraTree->Branch("detector",&detector);
//  AraTree->Branch("detector",&detector,6400,3);
//  AraTree->Branch("icemodel",&icemodel);
//  AraTree->Branch("settings1",&settings1);
//  AraTree2->Branch("detector",&detector);
//  AraTree2->Branch("spectra",&spectra);
//  AraTree->Branch("interaction1",&interaction1);







  for (int inu = 0; inu<30; inu++) {
      cout<<"\n";
      int pickunbiased_out = icemodel->PickUnbiased(0, interaction1, icemodel);
      cout<<"pickunbiased output : "<<pickunbiased_out<<endl;
      cout<<"nnu : ";
      interaction1->nnu.Print();
      cout<<"posnu : ";
      interaction1->posnu.Print();
      cout<<"test random # :"<<gRandom->Rndm()<<endl;
      cout<<"\n";
      cout<<"\n";
  }
  //
  //







  int nnu_pass = 0; // number of nu events which passed PickUnbiased.
  double posnuX[settings1->NNU];
  double posnuY[settings1->NNU];

  for (int inu=0;inu<settings1->NNU;inu++) { // loop over neutrinos

      // save X, Y of posnus which passed PickUnbiased
      if ( icemodel->PickUnbiased(0, interaction1, icemodel) ) {
          posnuX[nnu_pass] = interaction1->posnu.GetX();
          posnuY[nnu_pass] = interaction1->posnu.GetY();
          nnu_pass++;
      }







    
//    detector->resetDetector(); // set all signals on antennas to zero


//    Event *event=new Event(spectra->GetNuEnergy()); // creates a neutrino interaction
    // contains info about neutrino itself and the interaction
    // input neutrino energy
   


//    detector->simulateDetector(event,efficiencies);// simulate detector response
    // including waveforms for analysis
    // this also includes a trigger simulation



//    AraTree->Fill();

    //delete event;

  } // end loop over neutrinos
  cout<<"Total NNU : "<<settings1->NNU<<", PickUnbiased passed NNU : "<<nnu_pass<<endl;
    
  
  // test AraTree
  AraTree->Fill();
//  AraTree2->Fill();

//  AraTree->Write();
//  AraTree2->Write();
  AraFile->Write();
 efficiencies->summarize(); // summarize the results in an output file  












///////////////////////////////////////////
//  test Detector class
///////////////////////////////////////////




if ( settings1->DETECTOR == 0 ) {

cout<<"\n\t Test reading antenna array infomation !"<<endl;
//cout<<"total number of strings : "<<(int)detector->Detector.params.number_of_strings<<endl;
cout<<"total number of strings : "<<(int)detector->params.number_of_strings<<endl;
cout<<"total number of antennas : "<<(int)detector->params.number_of_antennas<<endl;
cout<<"\nantenna0 position is"<<endl;
cout<<"x : "<<(double)detector->strings[0].GetX()<<endl;
cout<<"y : "<<(double)detector->strings[0].GetY()<<endl;
cout<<"z1 : "<<(double)detector->strings[0].antennas[0].GetZ()<<" type : "<<(int)detector->strings[0].antennas[0].type<<endl;
cout<<"z2 : "<<(double)detector->strings[0].antennas[1].GetZ()<<" type : "<<(int)detector->strings[0].antennas[1].type<<endl;

cout<<"\nantenna1 position is"<<endl;
cout<<"x : "<<(double)detector->strings[1].GetX()<<endl;
cout<<"y : "<<(double)detector->strings[1].GetY()<<endl;
cout<<"z1 : "<<(double)detector->strings[1].antennas[0].GetZ()<<" type : "<<(int)detector->strings[1].antennas[0].type<<endl;
cout<<"z2 : "<<(double)detector->strings[1].antennas[1].GetZ()<<" type : "<<(int)detector->strings[1].antennas[1].type<<endl;

cout<<"\nantenna2 position is"<<endl;
cout<<"x : "<<(double)detector->strings[2].GetX()<<endl;
cout<<"y : "<<(double)detector->strings[2].GetY()<<endl;
cout<<"z1 : "<<(double)detector->strings[2].antennas[0].GetZ()<<" type : "<<(int)detector->strings[2].antennas[0].type<<endl;
cout<<"z2 : "<<(double)detector->strings[2].antennas[1].GetZ()<<" type : "<<(int)detector->strings[2].antennas[1].type<<endl;

cout<<"\nantenna3 position is"<<endl;
cout<<"x : "<<(double)detector->strings[3].GetX()<<endl;
cout<<"y : "<<(double)detector->strings[3].GetY()<<endl;
cout<<"z1 : "<<(double)detector->strings[3].antennas[0].GetZ()<<" type : "<<(int)detector->strings[3].antennas[0].type<<endl;
cout<<"z2 : "<<(double)detector->strings[3].antennas[1].GetZ()<<" type : "<<(int)detector->strings[3].antennas[1].type<<endl;

}






else if ( settings1->DETECTOR == 1 ) {

cout<<"\n\t Test ARA-N array setting !"<<endl;
//cout<<"total number of strings : "<<(int)detector->Detector.params.number_of_strings<<endl;
cout<<"total number of strings : "<<(int)detector->params.number_of_strings<<endl;
cout<<"total number of antennas : "<<(int)detector->params.number_of_antennas<<endl;


TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,2400,700);
c1->Divide(3,1);

double x[(int)detector->params.number_of_stations], y[(int)detector->params.number_of_stations];

for (int i=0;i<(int)detector->params.number_of_stations;i++) {
    x[i] = (double)detector->stations[i].GetX();
    y[i] = (double)detector->stations[i].GetY();
}

TGraph *gr;
gr = new TGraph((int)detector->params.number_of_stations,x,y);

c1->cd(1);
gr->SetTitle("ARA station layout");
gr->GetHistogram()->SetMaximum(10000);
gr->GetHistogram()->SetMinimum(-10000);
gr->GetXaxis()->SetLimits(-10000,10000);
gr->GetHistogram()->SetXTitle("X (m)");
gr->GetHistogram()->SetYTitle("Y (m)");
gr->Draw("a*");


c1->cd(2);

int station_choice = 0;
//int station_choice = (int)detector->params.number_of_stations - 1;

double string_x[4], string_y[4];
double surface_x[4], surface_y[4];

for (int i=0;i<4;i++) {
    string_x[i] = (double)detector->stations[station_choice].strings[i].GetX();
    string_y[i] = (double)detector->stations[station_choice].strings[i].GetY();

    surface_x[i] = (double)detector->stations[station_choice].surfaces[i].GetX();
    surface_y[i] = (double)detector->stations[station_choice].surfaces[i].GetY();
}

TGraph *gr_string;
gr_string = new TGraph(4,string_x,string_y);

gr_string->SetTitle("Strings and surface antennas layout for each station");
//gr_string->GetHistogram()->SetMaximum(5300);
//gr_string->GetHistogram()->SetMinimum(5100);
//gr_string->GetXaxis()->SetLimits(-3100,-2900);
gr_string->GetHistogram()->SetMaximum( (int)detector->stations[station_choice].GetY() + 100);
gr_string->GetHistogram()->SetMinimum( (int)detector->stations[station_choice].GetY() - 100);
gr_string->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
gr_string->SetMarkerColor(4);
gr_string->SetMarkerSize(2);
gr_string->SetMarkerStyle(20);
gr_string->GetHistogram()->SetXTitle("X (m)");
gr_string->GetHistogram()->SetYTitle("Y (m)");
gr_string->Draw("ap");

TGraph *gr_surface;
gr_surface = new TGraph(4,surface_x,surface_y);
gr_surface->SetMarkerColor(2);
gr_surface->SetMarkerSize(2);
gr_surface->SetMarkerStyle(21);
gr_surface->Draw("p");


TLegend *Leg_string_surface = new TLegend(1., 0.95, 0.5,0.8);
Leg_string_surface -> AddEntry(gr_string, "Strings");
Leg_string_surface -> AddEntry(gr_surface, "Surface antennas");
Leg_string_surface -> Draw();


c1->cd(3);

int string_choice = 0;
//int station_choice = (int)detector->params.number_of_stations - 1;

double antenna_x[4], antenna_y[4];   // use x as x, y as z to see the depth layout

for (int i=0;i<4;i++) {
    antenna_x[i] = (double)detector->stations[station_choice].strings[string_choice].GetX();
    antenna_y[i] = (double)detector->stations[station_choice].strings[string_choice].antennas[i].GetZ();
}

TGraph *gr_antenna;
gr_antenna = new TGraph(4,antenna_x,antenna_y);

gr_antenna->SetTitle("Borehole antenna layout for each string");
//gr_string->GetHistogram()->SetMaximum(5300);
//gr_string->GetHistogram()->SetMinimum(5100);
//gr_string->GetXaxis()->SetLimits(-3100,-2900);
gr_antenna->GetHistogram()->SetMaximum( (int)detector->stations[station_choice].strings[string_choice].GetZ() );
gr_antenna->GetHistogram()->SetMinimum( (int)detector->stations[station_choice].strings[string_choice].antennas[0].GetZ() - 20);
gr_antenna->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
gr_antenna->GetYaxis()->SetTitle("Z (depth, m)");
gr_antenna->SetMarkerColor(4);
gr_antenna->SetMarkerSize(2);
gr_antenna->SetMarkerStyle(20);
gr_antenna->GetHistogram()->SetXTitle("X (m)");
gr_antenna->GetHistogram()->SetYTitle("Z (m)");
gr_antenna->Draw("ap");


c1->Print("ARA-37_station_layout.pdf");






}












else if ( settings1->DETECTOR == 2 ) {

cout<<"\n\t Test ARA-37 array setting !"<<endl;
//cout<<"total number of strings : "<<(int)detector->Detector.params.number_of_strings<<endl;
cout<<"total number of strings : "<<(int)detector->params.number_of_strings<<endl;
cout<<"total number of antennas : "<<(int)detector->params.number_of_antennas<<endl;


TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,4000,700);
c1->Divide(5,1);

double x[(int)detector->params.number_of_stations], y[(int)detector->params.number_of_stations];

for (int i=0;i<(int)detector->params.number_of_stations;i++) {
    x[i] = (double)detector->stations[i].GetX();
    y[i] = (double)detector->stations[i].GetY();
}

TGraph *gr;
gr = new TGraph((int)detector->params.number_of_stations,x,y);

c1->cd(1);
gr->SetTitle("ARA station layout");
gr->GetHistogram()->SetMaximum(10000);
gr->GetHistogram()->SetMinimum(-10000);
gr->GetHistogram()->SetXTitle("X (m)");
gr->GetHistogram()->SetYTitle("Y (m)");
gr->GetYaxis()->SetTitleOffset(1.2);
gr->GetHistogram()->SetYTitle("Y (m)");
gr->GetXaxis()->SetLimits(-10000,10000);
gr->Draw("a*");


c1->cd(2);

int station_choice = 0;

double string_x[4], string_y[4];
double surface_x[4], surface_y[4];

for (int i=0;i<4;i++) {
    string_x[i] = (double)detector->stations[station_choice].strings[i].GetX();
    string_y[i] = (double)detector->stations[station_choice].strings[i].GetY();

    surface_x[i] = (double)detector->stations[station_choice].surfaces[i].GetX();
    surface_y[i] = (double)detector->stations[station_choice].surfaces[i].GetY();
}

TGraph *gr_string;
gr_string = new TGraph(4,string_x,string_y);

gr_string->SetTitle("Strings and surface antennas layout for each station");
//gr_string->GetHistogram()->SetMaximum(5300);
//gr_string->GetHistogram()->SetMinimum(5100);
//gr_string->GetXaxis()->SetLimits(-3100,-2900);
gr_string->GetHistogram()->SetMaximum( (int)detector->stations[station_choice].GetY() + 100);
gr_string->GetHistogram()->SetMinimum( (int)detector->stations[station_choice].GetY() - 100);
gr_string->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
gr_string->SetMarkerColor(4);
gr_string->SetMarkerSize(2);
gr_string->SetMarkerStyle(20);
gr_string->GetHistogram()->SetXTitle("X (m)");
gr_string->GetHistogram()->SetYTitle("Y (m)");
gr_string->GetYaxis()->SetTitleOffset(1.2);
gr_string->Draw("ap");

TGraph *gr_surface;
gr_surface = new TGraph(4,surface_x,surface_y);
gr_surface->SetMarkerColor(2);
gr_surface->SetMarkerSize(2);
gr_surface->SetMarkerStyle(21);
gr_surface->Draw("p");


TLegend *Leg_string_surface = new TLegend(1., 0.95, 0.5,0.8);
Leg_string_surface -> AddEntry(gr_string, "Strings");
Leg_string_surface -> AddEntry(gr_surface, "Surface antennas");
Leg_string_surface -> Draw();



c1->cd(3);

int string_choice = 0;
//int station_choice = (int)detector->params.number_of_stations - 1;

double antenna_x[4], antenna_y[4];   // use x as x, y as z to see the depth layout

for (int i=0;i<4;i++) {
    antenna_x[i] = (double)detector->stations[station_choice].strings[string_choice].GetX();
    antenna_y[i] = (double)detector->stations[station_choice].strings[string_choice].antennas[i].GetZ();
}

TGraph *gr_antenna;
gr_antenna = new TGraph(4,antenna_x,antenna_y);

gr_antenna->SetTitle("Borehole antenna layout for each string");
//gr_string->GetHistogram()->SetMaximum(5300);
//gr_string->GetHistogram()->SetMinimum(5100);
//gr_string->GetXaxis()->SetLimits(-3100,-2900);
gr_antenna->GetHistogram()->SetMaximum( (int)detector->stations[station_choice].strings[string_choice].GetZ() );
gr_antenna->GetHistogram()->SetMinimum( (int)detector->stations[station_choice].strings[string_choice].antennas[0].GetZ() - 20);
gr_antenna->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
gr_antenna->GetYaxis()->SetTitle("Z (depth, m)");
gr_antenna->SetMarkerColor(4);
gr_antenna->SetMarkerSize(2);
gr_antenna->SetMarkerStyle(20);
gr_antenna->GetHistogram()->SetXTitle("X (m)");
gr_antenna->GetYaxis()->SetTitleOffset(1.2);
gr_antenna->Draw("ap");



c1->cd(4);


double station_x[(int)detector->params.number_of_stations], station_z[(int)detector->params.number_of_stations];

for (int i=0;i<(int)detector->params.number_of_stations;i++) {
    station_x[i] = (double)detector->stations[i].GetX();
    station_z[i] = (double)detector->stations[i].GetZ();
//    station_z[i] = (double)detector->stations[0].GetZ();
}

TGraph *gr_crosssection;
gr_crosssection = new TGraph((int)detector->params.number_of_stations,station_x,station_z);

gr_crosssection->SetTitle("Station layout CrossSection");
gr_crosssection->GetHistogram()->SetMaximum( (int)detector->stations[0].GetZ() + 10 );
gr_crosssection->GetHistogram()->SetMinimum( (int)detector->stations[0].GetZ() - 10 );
gr_crosssection->GetXaxis()->SetLimits(-10000,10000);
gr_crosssection->GetHistogram()->SetXTitle("X (m)");
gr_crosssection->GetHistogram()->SetYTitle("Z (m)");
gr_crosssection->GetYaxis()->SetTitleOffset(1.2);
gr_crosssection->Draw("a*");


c1->cd(5);

TGraph *gr_posnu;
gr_posnu = new TGraph(nnu_pass,posnuX,posnuY);

gr_posnu->SetTitle("posnu");
//gr_posnu->GetHistogram()->SetMaximum( icemodel->Surface(0.,0.)*sin(icemodel->COASTLINE*RADDEG) );
gr_posnu->GetHistogram()->SetMaximum( icemodel->Surface(0.,0.)*sin(30.*RADDEG) );
gr_posnu->GetHistogram()->SetMinimum( -icemodel->Surface(0.,0.)*sin(30.*RADDEG)  );
gr_posnu->GetHistogram()->SetXTitle("X (m)");
gr_posnu->GetHistogram()->SetYTitle("Y (m)");
gr_posnu->GetYaxis()->SetTitleOffset(1.2);
gr_posnu->GetHistogram()->SetYTitle("Y (m)");
gr_posnu->GetXaxis()->SetLimits(-icemodel->Surface(0.,0.)*sin(30.*RADDEG),icemodel->Surface(0.,0.)*sin(30.*RADDEG));
gr_posnu->Draw("a*");




c1->Print("ARA-37_station_layout.pdf");



}



// test
cout<<"station[0] x : "<<detector->stations[0].GetX()<<endl;
cout<<"string[0] x : "<<detector->stations[0].strings[0].GetX()<<endl;
cout<<"antenna[0] x : "<<detector->stations[0].strings[0].antennas[0].GetX()<<endl;
cout<<"antenna[0] Gain(700,5,0) : "<<detector->stations[0].strings[0].antennas[0].GetG(detector,700.,5.,0.)<<endl;
cout<<"GetGain(700,5,0,0) : "<<detector->GetGain(700.,5.,0.,0)<<endl;
cout<<"GetGain(10,5,0,0) : "<<detector->GetGain(10.,5.,0.,0)<<endl;

cout<<"params.number_of_stations : "<<detector->params.number_of_stations<<endl;
cout<<"params.station_spacing : "<<detector->params.station_spacing<<endl;

/////////////////////////////////////////////



double *energy = spectra->Getenergy();
int Ebin = spectra->GetE_bin();

cout<<"\n";
for (int i=0;i<Ebin;i++) {
    cout<<"energy["<<i<<"] : "<<energy[i]<<endl;
}

///////////////////////////////////////////////////

double E = 19.0;
cout<<"\n";
cout<<"Flux at 19 is : "<<spectra->GetEdNdEdAdt(E)<<endl;

///////////////////////////////////////////////////


TSpline3 *sp1;
sp1 = spectra->GetSEdNdEdAdt();

///////////////////////////////////////////////////





TGraph *GEdN;
GEdN = spectra->GetGEdNdEdAdt();

TCanvas *c2 = new TCanvas("c2","A Simple Graph Example",200,10,1000,700);
c2 -> cd();
GEdN->Draw("al");

sp1->SetLineColor(2);
sp1->Draw("c same");
c2 -> Print("GEdN.pdf");

//////////////////////////////////////////////////






 delete icemodel;
 delete efficiencies;
 
 delete detector;

 delete spectra;

 test();


 return 0;
  
} //end main



void test() {

  cout << "test is " << 0 << "\n";
}
