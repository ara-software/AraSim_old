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
#include "Report.h"

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
  Event *event = 0;
  Report *report = 0;
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
  AraTree2->SetBranchAddress("event",&event);
  AraTree2->SetBranchAddress("report",&report);
  cout<<"branch detector"<<endl;
  
  AraTree->GetEvent(0);
  cout<<"getevent"<<endl;
  cout << "I'm here.\n";
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


cout<<"Detector -> freq_forfft[0] : "<<detector->freq_forfft[0]<<endl;
cout<<"Detector -> freq_forfft[9] : "<<detector->freq_forfft[9]<<endl;
cout<<"Detector -> freq_forfft[100] : "<<detector->freq_forfft[100]<<endl;

cout<<"icemodel surface : "<<icemodel->Surface(0.,0.)<<endl;

AraTree2->GetEvent(0);
cout<<"nnu x : "<<event->nnu.GetX()<<endl;
AraTree2->GetEvent(1);
cout<<"nnu x : "<<event->nnu.GetX()<<endl;

  int nnu_pass = 0; // number of nu events which passed PickUnbiased.
  double posnuX[settings->NNU];
  double posnuY[settings->NNU];
  double posnuR[settings->NNU];

  for (int inu=0;inu<settings->NNU;inu++) { // loop over neutrinos


      AraTree2->GetEvent(inu);

      // save X, Y of posnus which passed PickUnbiased
      if ( event->Nu_Interaction[0].pickposnu ) {
          posnuX[nnu_pass] = event->Nu_Interaction[0].posnu.GetX();
          posnuY[nnu_pass] = event->Nu_Interaction[0].posnu.GetY();
          posnuR[nnu_pass] = event->Nu_Interaction[0].posnu.R();
          nnu_pass++;

          cout<<"evt no "<<inu<<"stations[0].strings[1].antennas[2].ray_sol_cnt : "<<report->stations[0].strings[1].antennas[2].ray_sol_cnt<<endl;

/*
          if ( interaction->ray_solver_toggle ) {   // if ray_solver succeeded to get soutions
              cout<<"pass evt : "<<nnu_pass<<"\t";
              for (int i=0; i<interaction->ray_output[0].size(); i++) {
                  for (int j=0; j<interaction->ray_output.size(); j++) {
                      cout<<j<<"th ray_output : "<<interaction->ray_output[j][i]<<"\t";
                  }
                  cout<<"\n";
              }
          }// end if ray_solver_toggle
*/
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




/*

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





*/






//else if ( settings->DETECTOR == 2 ) {
else {

cout<<"\n\t Test ARA-37 array setting !"<<endl;
//cout<<"total number of strings : "<<(int)detector->Detector.params.number_of_strings<<endl;
cout<<"total number of strings : "<<(int)detector->params.number_of_strings<<endl;
cout<<"total number of antennas : "<<(int)detector->params.number_of_antennas<<endl;


TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,4800,700);
c1->Divide(6,1);

double x[(int)detector->params.number_of_stations], y[(int)detector->params.number_of_stations];

for (int i=0;i<(int)detector->params.number_of_stations;i++) {
    x[i] = (double)detector->stations[i].GetX();
    y[i] = (double)detector->stations[i].GetY();
}

TGraph *gr;
gr = new TGraph((int)detector->params.number_of_stations,x,y);

c1->cd(1);
gr->SetTitle("ARA station layout");
gr->GetHistogram()->SetMaximum(detector->params.core_y + 10000);
gr->GetHistogram()->SetMinimum(detector->params.core_y - 10000);
gr->GetHistogram()->SetXTitle("X (m)");
gr->GetHistogram()->SetYTitle("Y (m)");
gr->GetYaxis()->SetTitleOffset(1.2);
gr->GetHistogram()->SetYTitle("Y (m)");
gr->GetXaxis()->SetLimits(detector->params.core_x-10000, detector->params.core_x+10000);
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
//--------------------------------------------------
// gr_crosssection->GetHistogram()->SetMaximum( (int)detector->stations[0].GetZ() + 10 );
// gr_crosssection->GetHistogram()->SetMinimum( (int)detector->stations[0].GetZ() - 10 );
//-------------------------------------------------- 
gr_crosssection->GetHistogram()->SetMaximum( detector->stations[detector->params.number_of_stations-1].GetZ() + 1000 );
gr_crosssection->GetHistogram()->SetMinimum( detector->stations[detector->params.number_of_stations-1].GetZ() - 1000 );
gr_crosssection->GetXaxis()->SetLimits(detector->params.core_x-10000, detector->params.core_x+10000);
gr_crosssection->GetHistogram()->SetXTitle("X (m)");
gr_crosssection->GetHistogram()->SetYTitle("Z (m)");
gr_crosssection->GetYaxis()->SetTitleOffset(1.2);
gr_crosssection->Draw("a*");


c1->cd(5);

TGraph *gr_posnu;
gr_posnu = new TGraph(nnu_pass,posnuX,posnuY);

gr_posnu->SetTitle("posnu");
//--------------------------------------------------
// gr_posnu->GetHistogram()->SetMaximum( icemodel->Surface(0.,0.)*sin(30.*RADDEG) );
// gr_posnu->GetHistogram()->SetMinimum( -icemodel->Surface(0.,0.)*sin(30.*RADDEG)  );
//-------------------------------------------------- 
gr_posnu->GetHistogram()->SetMaximum( detector->params.core_y+3000. );
gr_posnu->GetHistogram()->SetMinimum( detector->params.core_y-3000.  );
gr_posnu->GetHistogram()->SetXTitle("X (m)");
gr_posnu->GetHistogram()->SetYTitle("Y (m)");
gr_posnu->GetYaxis()->SetTitleOffset(1.2);
gr_posnu->GetHistogram()->SetYTitle("Y (m)");
gr_posnu->GetXaxis()->SetLimits( detector->params.core_x-3000., detector->params.core_x+3000.);
//--------------------------------------------------
// gr_posnu->GetXaxis()->SetLimits(-icemodel->Surface(0.,0.)*sin(30.*RADDEG),icemodel->Surface(0.,0.)*sin(30.*RADDEG));
//-------------------------------------------------- 
gr_posnu->Draw("a*");


c1->cd(6);

double cd6_x[settings->NNU];
double cd6_y_top[settings->NNU];
double cd6_y_bot[settings->NNU];

for(int i =0; i<settings->NNU; i++) {
      
    AraTree2->GetEvent(i);

    cd6_x[i] = i;
    cd6_y_top[i] = icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() );
    cd6_y_bot[i] = icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() ) - icemodel->IceThickness( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() );
    
    if (icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() ) - posnuR[i] < 0) { 
        cout<<"Surface - posnuR : "<<icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() ) - posnuR[i]<<endl;
        cout<<"!offsurface"<<endl;
    }

    if (posnuR[i] - ( icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat()) - icemodel->IceThickness(event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat()) ) < 0)  {
        cout<<"posnuR - Icebottom : "<<posnuR[i] - ( icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() ) - icemodel->IceThickness(event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat()) ) <<endl;
        cout<<"!offbottomice"<<endl;
    }

}


TGraph *gr_depth;
gr_depth = new TGraph(settings->NNU, cd6_x, posnuR);

TGraph *gr_top;
gr_top = new TGraph(settings->NNU, cd6_x, cd6_y_top);

TGraph *gr_bot;
gr_bot = new TGraph(settings->NNU, cd6_x, cd6_y_bot);


gr_depth->SetTitle("posnu_position");
gr_depth->GetHistogram()->SetMaximum( cd6_y_top[0] + 1000. );
gr_depth->GetHistogram()->SetMinimum( cd6_y_bot[0] - 1000.  );
gr_depth->GetHistogram()->SetXTitle("posnu evt number");
gr_depth->GetHistogram()->SetYTitle("depth (m)");
gr_depth->GetYaxis()->SetTitleOffset(1.2);
gr_depth->SetMarkerColor(2);
gr_depth->SetMarkerSize(1);
gr_depth->SetMarkerStyle(21);
gr_depth->Draw("ap");

gr_top->SetMarkerColor(3);
gr_top->SetMarkerSize(1);
gr_top->SetMarkerStyle(21);
gr_top->Draw("p");

gr_bot->SetMarkerColor(4);
gr_bot->SetMarkerSize(1);
gr_bot->SetMarkerStyle(21);
gr_bot->Draw("p");


TLegend *Leg_top_bot = new TLegend(1., 0.95, 0.5,0.8);
Leg_top_bot -> AddEntry(gr_top, "Ice Surface");
Leg_top_bot -> AddEntry(gr_depth, "posnu depth");
Leg_top_bot -> AddEntry(gr_bot, "Bedrock");
Leg_top_bot -> Draw();

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


//--------------------------------------------------
// TSpline3 *sp1;
// sp1 = spectra->GetSEdNdEdAdt();
//-------------------------------------------------- 

///////////////////////////////////////////////////





TGraph *GEdN;
//--------------------------------------------------
// GEdN = spectra->GetGEdNdEdAdt();
//-------------------------------------------------- 
GEdN = new TGraph( spectra->GetE_bin(), spectra->Getenergy(), spectra->GetEdNdEdAdt() );

TCanvas *c2 = new TCanvas("c2","A Simple Graph Example",200,10,1000,700);
c2 -> cd();
c2->SetLogy();
GEdN->SetTitle("Neutrino flux (ESS)");
GEdN->GetHistogram()->SetXTitle("log E");
GEdN->GetHistogram()->SetYTitle("Flux EdNdEdAdt (cm^2*str*s)");
GEdN->Draw("al");

//--------------------------------------------------
// sp1->SetLineColor(2);
// sp1->Draw("c same");
//-------------------------------------------------- 
c2 -> Print("GEdN.pdf");

//////////////////////////////////////////////////

// roughly 1 deg from the south pole, (approx 100km from south pole)
// ang resolution 0.1 deg.

int ang_step = 100;
double max_ang = 10.; //1 deg lat
//--------------------------------------------------
// double max_ang = 1.; //1 deg lat
//-------------------------------------------------- 

double lat[ang_step];
double Surf[ang_step];
double surf_abv_geo[ang_step];
double ice_bot[ang_step];
double ice_bot_ex[ang_step];

double ice_bot_min = 0.;
double ice_bot_min_ex = icemodel->Surface(0.,0.);;

for (int i=0;i<ang_step;i++) {
    lat[i] = (max_ang/(double)ang_step) * (double)i;
    Surf[i] = icemodel->Surface(0.,lat[i]);
    surf_abv_geo[i] = icemodel->SurfaceAboveGeoid(0.,lat[i]);
    ice_bot[i] = surf_abv_geo[i] - icemodel->IceThickness(0.,lat[i]);
    ice_bot_ex[i] = Surf[i] - icemodel->IceThickness(0.,lat[i]);
    if (ice_bot[i] < ice_bot_min) ice_bot_min = ice_bot[i];
    if (ice_bot_ex[i] < ice_bot_min_ex) ice_bot_min_ex = ice_bot_ex[i];
}

TGraph *G_surf_abv_geo;
G_surf_abv_geo = new TGraph(ang_step, lat, surf_abv_geo);
    
TGraph *G_ice_bot;
G_ice_bot = new TGraph(ang_step, lat, ice_bot);

TCanvas *cGeo = new TCanvas("cGeo","A Simple Graph Example",200,10,2000,700);
cGeo -> Divide(2,1);
cGeo -> cd(1);
G_surf_abv_geo->SetTitle("Ice surface and Ice thickness wrt Geoid");
G_surf_abv_geo->GetHistogram()->SetMaximum( surf_abv_geo[0] + 100);
G_surf_abv_geo->GetHistogram()->SetMinimum( ice_bot_min - 100);
//--------------------------------------------------
// G_surf_abv_geo->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
//-------------------------------------------------- 
G_surf_abv_geo->SetMarkerColor(4);
G_surf_abv_geo->SetMarkerSize(2);
G_surf_abv_geo->SetMarkerStyle(20);
G_surf_abv_geo->GetHistogram()->SetXTitle("Lat (deg)");
G_surf_abv_geo->GetHistogram()->SetYTitle("Z (m)");
G_surf_abv_geo->Draw("ap");

G_ice_bot->SetMarkerColor(2);
G_ice_bot->SetMarkerSize(2);
G_ice_bot->SetMarkerStyle(21);
G_ice_bot->Draw("p");

TLegend *Leg_Geo = new TLegend(1., 0.95, 0.5,0.8);
Leg_Geo -> AddEntry(G_surf_abv_geo, "Ice surface");
Leg_Geo -> AddEntry(G_ice_bot, "Ice bottom");
Leg_Geo -> Draw();


    
cGeo -> cd(2);

TGraph *G_surf_abv_geo2;
G_surf_abv_geo2 = new TGraph(ang_step, lat, Surf);
    
TGraph *G_ice_bot2;
G_ice_bot2 = new TGraph(ang_step, lat, ice_bot_ex);

G_surf_abv_geo2->SetTitle("Ice surface and Ice thickness");
G_surf_abv_geo2->GetHistogram()->SetMaximum( Surf[0] + 100);
G_surf_abv_geo2->GetHistogram()->SetMinimum( ice_bot_min_ex - 100);
//--------------------------------------------------
// G_surf_abv_geo->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
//-------------------------------------------------- 
G_surf_abv_geo2->SetMarkerColor(4);
G_surf_abv_geo2->SetMarkerSize(2);
G_surf_abv_geo2->SetMarkerStyle(20);
G_surf_abv_geo2->GetHistogram()->SetXTitle("Lat (deg)");
G_surf_abv_geo2->GetHistogram()->SetYTitle("Z (m)");
G_surf_abv_geo2->Draw("ap");

G_ice_bot2->SetMarkerColor(2);
G_ice_bot2->SetMarkerSize(2);
G_ice_bot2->SetMarkerStyle(21);
G_ice_bot2->Draw("p");

//--------------------------------------------------
// TLegend *Leg_Geo2 = new TLegend(1., 0.95, 0.5,0.8);
// Leg_Geo2 -> AddEntry(G_surf_abv_geo2, "Ice surface");
// Leg_Geo2 -> AddEntry(G_ice_bot2, "Ice bottom");
// Leg_Geo2 -> Draw();
//-------------------------------------------------- 




cGeo -> Print("GEOID1.pdf");



//////////////////////////////////////////////////


int evt_n, station_n, string_n, antenna_n, ray_sol_n;

evt_n = 0;
station_n = 0;
string_n = 1;
//string_n = 0;
antenna_n = 1;
//antenna_n = 0;
ray_sol_n = 0;

AraTree2->GetEvent(evt_n);

double time[settings->NFOUR/2];
double V[settings->NFOUR/2];
double V_org[settings->NFOUR/2];

cout<<"view angle : "<<report->stations[station_n].strings[string_n].antennas[antenna_n].view_ang[ray_sol_n]*DEGRAD<<endl;

for (int i=0;i<settings->NFOUR/2;i++) {
    time[i] = report->stations[station_n].strings[string_n].antennas[antenna_n].time[ray_sol_n][i];
    V_org[i] = report->stations[station_n].strings[string_n].antennas[antenna_n].V[ray_sol_n][i];
    if (i<settings->NFOUR/4) {
        V[i+settings->NFOUR/4] = report->stations[station_n].strings[string_n].antennas[antenna_n].V[ray_sol_n][i];
    }
    else {
        V[i-settings->NFOUR/4] = report->stations[station_n].strings[string_n].antennas[antenna_n].V[ray_sol_n][i];
    }
}


TGraph *G_V_time;
G_V_time = new TGraph(settings->NFOUR/2, time, V);

TGraph *G_V_time_org;
G_V_time_org = new TGraph(settings->NFOUR/2, time, V_org);

TCanvas *cV_time = new TCanvas("cV_time","V(t)", 200,10,1400,700);
cV_time->Divide(2,1);
cV_time -> cd(1);
G_V_time_org->SetTitle("V(t) for evt %d, station[%d].string[%d].antenna[%d] (org)");
G_V_time_org->GetHistogram()->SetXTitle("time (s)");
G_V_time_org->GetHistogram()->SetYTitle("Voltage (V)");
G_V_time_org->Draw("al");

cV_time -> cd(2);
G_V_time->SetTitle("V(t) for evt 0, station[0].string[0].antenna[0] (fixed?)");
G_V_time->GetHistogram()->SetXTitle("time (s)");
G_V_time->GetHistogram()->SetYTitle("Voltage (V)");
G_V_time->Draw("al");

cV_time -> Print("V_time_evt0.pdf");


TGraph *G_V_time_zoom;
G_V_time_zoom = new TGraph(settings->NFOUR/2, time, V);

TCanvas *cV_time_zoom = new TCanvas("cV_time_zoom","V(t) for evt 0, station0.string0.antenna0",200,10,1400,700);
cV_time_zoom->Divide(2,1);
cV_time_zoom -> cd(1);
G_V_time->Draw("al");

cV_time_zoom -> cd(2);
G_V_time_zoom->SetTitle("V(t) for evt 0, station[0].string[0].antenna[0] (Zoomed)");
G_V_time_zoom->GetHistogram()->SetXTitle("time (s)");
G_V_time_zoom->GetHistogram()->SetYTitle("Voltage (V)");
G_V_time_zoom->GetXaxis()->SetLimits(8E-8,12E-8);
G_V_time_zoom->Draw("al");

cV_time_zoom -> Print("V_time_evt0_zoom.pdf");



int N_freq;
N_freq = detector->GetFreqBin();

double freq[N_freq];
double vmmhz_freq[N_freq];



for (int i=0;i<N_freq;i++) {
    freq[i] = detector->GetFreq(i); // freq in Hz
    vmmhz_freq[i] = report->stations[station_n].strings[string_n].antennas[antenna_n].vmmhz[ray_sol_n][i];
}


TGraph *G_vmmhz_freq;
G_vmmhz_freq = new TGraph(N_freq, freq, vmmhz_freq);


TCanvas *cVmMHz = new TCanvas("cVmMHz","V/m/MHz", 200,10,1000,700);

cVmMHz -> cd();
cVmMHz -> SetLogy();

G_vmmhz_freq->SetTitle("VmMHz for evt 0, station[0].string[0].antenna[0] (before antenna)");
G_vmmhz_freq->GetHistogram()->SetXTitle("freq (Hz)");
G_vmmhz_freq->GetHistogram()->SetYTitle("Signal (V/m/MHz)");
G_vmmhz_freq->Draw("al");

cVmMHz -> Print("VmMHz_evt0.pdf");


AraTree2->GetEvent(3);
for (int i=0; i<detector->params.number_of_stations; i++) {
    for (int j=0; j<detector->params.number_of_strings_station; j++) {
        for (int k=0; k<detector->params.number_of_antennas_string; k++) {
            cout<<"Rank of stations["<<i<<"].strings["<<j<<"].antennas["<<k<<"] = "<<report->stations[i].strings[j].antennas[k].Rank[0]<<endl;
            cout<<"PeakV of stations["<<i<<"].strings["<<j<<"].antennas["<<k<<"] = "<<report->stations[i].strings[j].antennas[k].PeakV[0]<<endl;
        }
    }
}




}


