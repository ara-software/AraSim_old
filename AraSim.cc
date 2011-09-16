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


//  IceModel *icemodel=new IceModel(ICE_MODEL + NOFZ*10,CONSTANTICETHICKNESS * 1000 + CONSTANTCRUST * 100 + FIXEDELEVATION * 10 + 0,MOOREBAY);// creates Antarctica ice model
  IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
  //IceModel inherits from EarthModel  

    int mode = 2;

//  Detector *detector=new Detector(); // builds antenna array
  Detector *detector=new Detector(mode); // builds antenna array, 0 for testbed
//  Detector *detector=new Detector(1); // builds antenna array, 0 for testbed

  //  Trigger *trigger=new Trigger(detector); // builds the trigger  
//  Efficiencies *efficiencies=new Efficiencies(detector->getnRx(),outputdir); // keeps track of efficiencies at each stage of the simulation
  Efficiencies *efficiencies=new Efficiencies(100,outputdir); // keeps track of efficiencies at each stage of the simulation
//  Spectra *spectra=new Spectra(EXPONENT); // gets library (or whatever) of neutrino spectra
  
//  Spectra *spectra=new Spectra(settings1->EXPONENT); // gets library (or whatever) of neutrino spectra

  
  TFile *AraFile=new TFile((outputdir+"/AraOut.root").c_str(),"RECREATE","ara");
  TTree *AraTree=new TTree("AraTree","AraTree");
//  AraTree->Branch("detector",&detector);


//  for (int inu=0;inu<NNU;inu++) { // loop over neutrinos
  for (int inu=0;inu<settings1->NNU;inu++) { // loop over neutrinos

    
//    detector->resetDetector(); // set all signals on antennas to zero


//    Event *event=new Event(spectra->GetNuEnergy()); // creates a neutrino interaction
    // contains info about neutrino itself and the interaction
    // input neutrino energy
   


//    detector->simulateDetector(event,efficiencies);// simulate detector response
    // including waveforms for analysis
    // this also includes a trigger simulation



    AraTree->Fill();

    //delete event;
  }
  AraTree->Write();
  AraFile->Close();
 efficiencies->summarize(); // summarize the results in an output file  


///////////////////////////////////////////
//  test Detector class
///////////////////////////////////////////

if ( mode == 0 ) {

cout<<"\n\t Test reading antenna array infomation !"<<endl;
//cout<<"total number of strings : "<<(int)detector->Detector.params.number_of_strings<<endl;
cout<<"total number of strings : "<<(int)detector->params.number_of_strings<<endl;
cout<<"total number of antennas : "<<(int)detector->params.number_of_antennas<<endl;
cout<<"\nantenna0 position is"<<endl;
cout<<"x : "<<(float)detector->strings[0].x<<endl;
cout<<"y : "<<(float)detector->strings[0].y<<endl;
cout<<"z1 : "<<(float)detector->strings[0].antennas[0].z<<" model : "<<(int)detector->strings[0].antennas[0].model<<endl;
cout<<"z2 : "<<(float)detector->strings[0].antennas[1].z<<" model : "<<(int)detector->strings[0].antennas[1].model<<endl;

cout<<"\nantenna1 position is"<<endl;
cout<<"x : "<<(float)detector->strings[1].x<<endl;
cout<<"y : "<<(float)detector->strings[1].y<<endl;
cout<<"z1 : "<<(float)detector->strings[1].antennas[0].z<<" model : "<<(int)detector->strings[1].antennas[0].model<<endl;
cout<<"z2 : "<<(float)detector->strings[1].antennas[1].z<<" model : "<<(int)detector->strings[1].antennas[1].model<<endl;

cout<<"\nantenna2 position is"<<endl;
cout<<"x : "<<(float)detector->strings[2].x<<endl;
cout<<"y : "<<(float)detector->strings[2].y<<endl;
cout<<"z1 : "<<(float)detector->strings[2].antennas[0].z<<" model : "<<(int)detector->strings[2].antennas[0].model<<endl;
cout<<"z2 : "<<(float)detector->strings[2].antennas[1].z<<" model : "<<(int)detector->strings[2].antennas[1].model<<endl;

cout<<"\nantenna3 position is"<<endl;
cout<<"x : "<<(float)detector->strings[3].x<<endl;
cout<<"y : "<<(float)detector->strings[3].y<<endl;
cout<<"z1 : "<<(float)detector->strings[3].antennas[0].z<<" model : "<<(int)detector->strings[3].antennas[0].model<<endl;
cout<<"z2 : "<<(float)detector->strings[3].antennas[1].z<<" model : "<<(int)detector->strings[3].antennas[1].model<<endl;

}






else if ( mode == 1 ) {

cout<<"\n\t Test ARA-N array setting !"<<endl;
//cout<<"total number of strings : "<<(int)detector->Detector.params.number_of_strings<<endl;
cout<<"total number of strings : "<<(int)detector->params.number_of_strings<<endl;
cout<<"total number of antennas : "<<(int)detector->params.number_of_antennas<<endl;


TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,1400,700);
c1->Divide(2,1);

float x[(int)detector->params.number_of_stations], y[(int)detector->params.number_of_stations];

for (int i=0;i<(int)detector->params.number_of_stations;i++) {
    x[i] = (float)detector->stations[i].x;
    y[i] = (float)detector->stations[i].y;
}

TGraph *gr;
gr = new TGraph((int)detector->params.number_of_stations,x,y);

c1->cd(1);
gr->SetTitle("ARA station layout");
gr->GetHistogram()->SetMaximum(10000);
gr->GetHistogram()->SetMinimum(-10000);
gr->GetXaxis()->SetLimits(-10000,10000);
gr->Draw("a*");


c1->cd(2);

int station_choice = 0;
//int station_choice = (int)detector->params.number_of_stations - 1;

float string_x[4], string_y[4];
float surface_x[4], surface_y[4];

for (int i=0;i<4;i++) {
    string_x[i] = (float)detector->stations[station_choice].ARA_strings[i].x;
    string_y[i] = (float)detector->stations[station_choice].ARA_strings[i].y;

    surface_x[i] = (float)detector->stations[station_choice].surfs[i].x;
    surface_y[i] = (float)detector->stations[station_choice].surfs[i].y;
}

TGraph *gr_string;
gr_string = new TGraph(4,string_x,string_y);

gr_string->SetTitle("Strings and surface antennas layout for each station");
//gr_string->GetHistogram()->SetMaximum(5300);
//gr_string->GetHistogram()->SetMinimum(5100);
//gr_string->GetXaxis()->SetLimits(-3100,-2900);
gr_string->GetHistogram()->SetMaximum( (int)detector->stations[station_choice].y + 100);
gr_string->GetHistogram()->SetMinimum( (int)detector->stations[station_choice].y - 100);
gr_string->GetXaxis()->SetLimits( (int)detector->stations[station_choice].x - 100, (int)detector->stations[station_choice].x + 100 );
gr_string->SetMarkerColor(4);
gr_string->SetMarkerSize(2);
gr_string->SetMarkerStyle(20);
gr_string->Draw("ap");

TGraph *gr_surface;
gr_surface = new TGraph(4,surface_x,surface_y);
gr_surface->SetMarkerColor(2);
gr_surface->SetMarkerSize(2);
gr_surface->SetMarkerStyle(21);
gr_surface->Draw("p");


TLegend *Leg_string_surf = new TLegend(1., 0.95, 0.5,0.8);
Leg_string_surf -> AddEntry(gr_string, "Strings");
Leg_string_surf -> AddEntry(gr_surface, "Surface antennas");
Leg_string_surf -> Draw();


c1->Print("ARA-37_station_layout.pdf");






}












else if ( mode == 2 ) {

cout<<"\n\t Test ARA-37 array setting !"<<endl;
//cout<<"total number of strings : "<<(int)detector->Detector.params.number_of_strings<<endl;
cout<<"total number of strings : "<<(int)detector->params.number_of_strings<<endl;
cout<<"total number of antennas : "<<(int)detector->params.number_of_antennas<<endl;


TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,1400,700);
c1->Divide(2,1);

float x[(int)detector->params.number_of_stations], y[(int)detector->params.number_of_stations];

for (int i=0;i<(int)detector->params.number_of_stations;i++) {
    x[i] = (float)detector->stations[i].x;
    y[i] = (float)detector->stations[i].y;
}

TGraph *gr;
gr = new TGraph((int)detector->params.number_of_stations,x,y);

c1->cd(1);
gr->SetTitle("ARA station layout");
gr->GetHistogram()->SetMaximum(10000);
gr->GetHistogram()->SetMinimum(-10000);
gr->GetXaxis()->SetLimits(-10000,10000);
gr->Draw("a*");


c1->cd(2);

int station_choice = 0;

float string_x[4], string_y[4];
float surface_x[4], surface_y[4];

for (int i=0;i<4;i++) {
    string_x[i] = (float)detector->stations[station_choice].ARA_strings[i].x;
    string_y[i] = (float)detector->stations[station_choice].ARA_strings[i].y;

    surface_x[i] = (float)detector->stations[station_choice].surfs[i].x;
    surface_y[i] = (float)detector->stations[station_choice].surfs[i].y;
}

TGraph *gr_string;
gr_string = new TGraph(4,string_x,string_y);

gr_string->SetTitle("Strings and surface antennas layout for each station");
//gr_string->GetHistogram()->SetMaximum(5300);
//gr_string->GetHistogram()->SetMinimum(5100);
//gr_string->GetXaxis()->SetLimits(-3100,-2900);
gr_string->GetHistogram()->SetMaximum( (int)detector->stations[station_choice].y + 100);
gr_string->GetHistogram()->SetMinimum( (int)detector->stations[station_choice].y - 100);
gr_string->GetXaxis()->SetLimits( (int)detector->stations[station_choice].x - 100, (int)detector->stations[station_choice].x + 100 );
gr_string->SetMarkerColor(4);
gr_string->SetMarkerSize(2);
gr_string->SetMarkerStyle(20);
gr_string->Draw("ap");

TGraph *gr_surface;
gr_surface = new TGraph(4,surface_x,surface_y);
gr_surface->SetMarkerColor(2);
gr_surface->SetMarkerSize(2);
gr_surface->SetMarkerStyle(21);
gr_surface->Draw("p");


TLegend *Leg_string_surf = new TLegend(1., 0.95, 0.5,0.8);
Leg_string_surf -> AddEntry(gr_string, "Strings");
Leg_string_surf -> AddEntry(gr_surface, "Surface antennas");
Leg_string_surf -> Draw();


c1->Print("ARA-37_station_layout.pdf");






}



/////////////////////////////////////////////







 delete icemodel;
 delete efficiencies;
 
 delete detector;

 test();


 return 0;
  
} //end main



void test() {

  cout << "test is " << 0 << "\n";
}
