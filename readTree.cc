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
#include "counting.hh"
#include "Primaries.h"

#include "Ray.h"

class EarthModel; //class


string outputdir="outputs";

int main() {


//  Settings *settings = new Settings();

  //  Detector *detector=new Detector(settings->DETECTOR); // builds antenna array, 0 for testbed
//  Detector *detector=0; // builds antenna array, 0 for testbed
  Detector *detector = 0; 
  Settings *settings = 0;
  Spectra *spectra = 0;
  IceModel *icemodel = 0;
  Interaction *interaction = 0;
  cout<<"construct detector"<<endl;

  
  TFile *AraFile=new TFile((outputdir+"/AraOut.root").c_str());
  cout<<"AraFile"<<endl;
  TTree *AraTree=(TTree*)AraFile->Get("AraTree");
  TTree *AraTree2=(TTree*)AraFile->Get("AraTree2");
  cout<<"AraTree"<<endl;
  AraTree->SetBranchAddress("detector",&detector);
  AraTree->SetBranchAddress("settings",&settings);
  AraTree->SetBranchAddress("spectra",&spectra);
  AraTree->SetBranchAddress("icemodel",&icemodel);
  AraTree2->SetBranchAddress("interaction",&interaction);
  cout<<"branch detector"<<endl;
  
  AraTree->GetEvent(0);
  cout<<"getevent"<<endl;
  cout << "I'm here.\n";
  cout << "Vgain is " << detector->Vgain[0][10] << "\n";
  cout << "GetGain(700,5,0,0) is " << detector->GetGain(700., 5., 0., 0) << "\n";
  cout << "GetGain(10,5,0,0) is " << detector->GetGain(10., 5., 0., 0) << "\n";

  cout<<"station x is "<<detector->stations[0].GetX()<<endl;
  cout<<"string x is "<<detector->stations[0].strings[0].GetX()<<endl;
  cout<<"antenna x is "<<detector->stations[0].strings[0].antennas[0].GetX()<<endl;
  cout<<"antenna Gain(700,5,0) is "<<detector->stations[0].strings[0].antennas[0].GetG(detector,700.,5.,0.)<<endl;

  cout<<"params.number_of_stations : "<<detector->params.number_of_stations<<endl;
  cout<<"params.station_spacing : "<<detector->params.station_spacing<<endl;

  cout<<"\n"<<endl;
  cout<<"Settings->NNU : "<<settings->NNU<<endl;
  cout<<"Settings->DETECTOR : "<<settings->DETECTOR<<endl;


cout<<"random energy from Spectra : "<<spectra->GetNuEnergy()<<endl;

cout<<"Detector static const double freq_init : "<<detector->Getfreq_init()<<endl;

cout<<"icemodel surface : "<<icemodel->Surface(0.,0.)<<endl;

AraTree2->GetEvent(0);
cout<<"interaction nnu : "<<interaction->nnu.GetX()<<endl;
AraTree2->GetEvent(1);
cout<<"interaction nnu : "<<interaction->nnu.GetX()<<endl;

  int nnu_pass = 0; // number of nu events which passed PickUnbiased.
  double posnuX[settings->NNU];
  double posnuY[settings->NNU];

  for (int inu=0;inu<settings->NNU;inu++) { // loop over neutrinos


      AraTree2->GetEvent(inu);

      // save X, Y of posnus which passed PickUnbiased
      if ( interaction->pickunbiased ) {
          posnuX[nnu_pass] = interaction->posnu.GetX();
          posnuY[nnu_pass] = interaction->posnu.GetY();
          nnu_pass++;
      }


  } // end loop over neutrinos



///////////////////////////////////////////
//  test Detector class
///////////////////////////////////////////




if ( settings->DETECTOR == 0 ) {

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






else if ( settings->DETECTOR == 1 ) {

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












else if ( settings->DETECTOR == 2 ) {

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

cout<<"Spectra energy values : "<<spectra->GetE_bin()<<endl;
for (int i=0;i<spectra->GetE_bin();i++) {
    cout<<"energy bin "<<i<<" : "<<spectra->energy[i]<<endl;
}

cout<<"IceModel R_EARTH : "<<icemodel->R_EARTH<<endl;

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















}


