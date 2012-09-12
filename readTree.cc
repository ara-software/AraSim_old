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
  Trigger *trigger = 0;
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
  AraTree->SetBranchAddress("trigger",&trigger);
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

      if ( report->stations[0].Global_Pass ) {
          for ( int j=0; j<detector->params.number_of_strings_station; j++) {
              for (int k=0; k<detector->params.number_of_antennas_string; k++) {
                  cout<<"noise_ID : "<<report->stations[0].strings[j].antennas[k].noise_ID[0]<<endl;
              }
          }
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

evt_n = 34;
station_n = 0;
string_n = 0;
antenna_n = 1;
ray_sol_n = 0;

AraTree2->GetEvent(evt_n);




int N_freq;
double NFOUR;
N_freq = detector->GetFreqBin();
NFOUR = (double)settings->NFOUR;

double dfreq_detector;  // dfreq at detector related freq bin
double dfreq_fft;       // dfreq at fft bin
double dfreq_databin;       // dfreq at data bin
dfreq_detector = (detector->GetFreq(1) - detector->GetFreq(0));
dfreq_fft = 1./ ((double)settings->NFOUR/2 * settings->TIMESTEP);
dfreq_databin = 1./ ((double)settings->DATA_BIN_SIZE * settings->TIMESTEP);

double freq[N_freq];
double vmmhz_freq[N_freq];
double vmmhz_antfactors[N_freq];
double vmmhz_filter[N_freq];

double freq_fft[settings->NFOUR/4];
double freq_databin[settings->DATA_BIN_SIZE/2];
double vmmhz_fft[settings->NFOUR/4];

//double Vfft_noise_before[settings->NFOUR/4];
//double Vfft_noise_after[settings->NFOUR/4];
double Vfft_noise_before[settings->DATA_BIN_SIZE/2];
double Vfft_noise_after[settings->DATA_BIN_SIZE/2];

double total_mean_power_noise_before = 0.;
double total_mean_power_noise_after = 0.;
double total_mean_power_signal = 0.;

//double Filter[settings->NFOUR/4];
double Filter[settings->DATA_BIN_SIZE/2];


cout<<"define"<<endl;

if ( report->stations[station_n].strings[string_n].antennas[antenna_n].ray_sol_cnt >= ray_sol_n && settings->DATA_SAVE_MODE==0) { // if there is ray sol && output file have all information



for (int i=0;i<N_freq;i++) {
    freq[i] = detector->GetFreq(i); // freq in Hz
    vmmhz_freq[i] = report->stations[station_n].strings[string_n].antennas[antenna_n].vmmhz[ray_sol_n][i] / (1.E6); // from V/m/MHz to V/m/Hz
    vmmhz_freq[i] = vmmhz_freq[i]*vmmhz_freq[i] * dfreq_detector * 2. / 50.;

    vmmhz_antfactors[i] = report->stations[station_n].strings[string_n].antennas[antenna_n].VHz_antfactor[ray_sol_n][i];
    vmmhz_antfactors[i] = vmmhz_antfactors[i]*vmmhz_antfactors[i] * dfreq_detector * 2. / 50.;

    vmmhz_filter[i] = report->stations[station_n].strings[string_n].antennas[antenna_n].VHz_filter[ray_sol_n][i];
    vmmhz_filter[i] = vmmhz_filter[i]*vmmhz_filter[i] * dfreq_detector * 2. / 50.;

}

cout<<"first loop"<<endl;

for (int i=0;i<settings->NFOUR/4;i++) {
    freq_fft[i] = (double)i * dfreq_fft;

    vmmhz_fft[i] = pow(report->stations[station_n].strings[string_n].antennas[antenna_n].Vfft[ray_sol_n][i*2],2.) + pow(report->stations[station_n].strings[string_n].antennas[antenna_n].Vfft[ray_sol_n][i*2 + 1],2.);

    vmmhz_fft[i] = vmmhz_fft[i] * pow( NFOUR/4.,2.) * 2. / ( ( NFOUR * NFOUR / 4.) * 50. * dfreq_fft );
    // first pow(NFOUR/4, 2) is to cancel out the factor which is meaningful after the inverse FFT (MakeArraysforFFT apply it before inverse FFT)
    // second 2/(NFOUR*NFOUR/4) is for 2/(N*N) in equation (for mean power)

    // dont know why they are different!!!
    //vmmhz_fft[i] *= 3.;
    //



    total_mean_power_signal += vmmhz_fft[i];

}

cout<<"second loop"<<endl;

for (int i=0; i<settings->DATA_BIN_SIZE/2; i++) {

    freq_databin[i] = (double)i * dfreq_databin;


    Vfft_noise_after[i] = pow(trigger->Vfft_noise_after[i*2],2.) + pow(trigger->Vfft_noise_after[i*2 + 1],2.);    // after have R, I part
    Vfft_noise_before[i] = pow(trigger->Vfft_noise_before[i],2.);

    Vfft_noise_after[i] = Vfft_noise_after[i] * 2. / ( ( settings->DATA_BIN_SIZE * settings->DATA_BIN_SIZE) * 50. * dfreq_databin );
    Vfft_noise_before[i] = Vfft_noise_before[i] * 2. / ( ( settings->DATA_BIN_SIZE * settings->DATA_BIN_SIZE) * 50. * dfreq_databin );


    //Filter[i] = pow(10., ( detector->GetFilterGain_fft(i) )/20.);   // from dB to unitless gain for voltage
    //cout<<"filter gain in db : "<<detector->GetFilterGain_databin(i)<<endl;
    Filter[i] = pow(10., ( detector->GetFilterGain_databin(i) )/20.);   // from dB to unitless gain for voltage

    total_mean_power_noise_before += Vfft_noise_before[i];
    total_mean_power_noise_after += Vfft_noise_after[i];

}



TGraph *G_vmmhz_freq;
G_vmmhz_freq = new TGraph(N_freq, freq, vmmhz_freq);

TGraph *G_vmmhz_factor;
G_vmmhz_factor = new TGraph(N_freq, freq, vmmhz_antfactors);

TGraph *G_vmmhz_filter;
G_vmmhz_filter = new TGraph(N_freq, freq, vmmhz_filter);

TGraph *G_vmmhz_fft;
G_vmmhz_fft = new TGraph(settings->NFOUR/4, freq_fft, vmmhz_fft);

TGraph *G_Vfft_noise_before;
//G_Vfft_noise_before = new TGraph(settings->NFOUR/4, freq_fft, Vfft_noise_before);
G_Vfft_noise_before = new TGraph(settings->DATA_BIN_SIZE/2, freq_databin, Vfft_noise_before);

TGraph *G_Vfft_noise_after;
//G_Vfft_noise_after = new TGraph(settings->NFOUR/4, freq_fft, Vfft_noise_after);
G_Vfft_noise_after = new TGraph(settings->DATA_BIN_SIZE/2, freq_databin, Vfft_noise_after);


TGraph *G_Filter;
//G_Filter = new TGraph(settings->NFOUR/4, freq_fft, Filter);
G_Filter = new TGraph(settings->DATA_BIN_SIZE/2, freq_databin, Filter);


TCanvas *cVmMHz = new TCanvas("cVmMHz","V/m/MHz", 200,10,4000,1400);
cVmMHz->Divide(4,2);

cVmMHz -> cd(1);
cVmMHz -> SetLogy();

G_vmmhz_freq->SetTitle("spectrum (before antenna W/m^2/Hz)");
G_vmmhz_freq->GetHistogram()->SetXTitle("freq (Hz)");
G_vmmhz_freq->GetHistogram()->SetYTitle("mean power spectrum (W/m^2/Hz)");
G_vmmhz_freq->Draw("al");



cVmMHz -> cd(2);
cVmMHz -> SetLogy();

G_vmmhz_factor->SetTitle("spectrum (after antenna W/Hz)");
G_vmmhz_factor->GetHistogram()->SetXTitle("freq (Hz)");
G_vmmhz_factor->GetHistogram()->SetYTitle("mean power spectrum (W/Hz)");
G_vmmhz_factor->Draw("al");

cVmMHz -> cd(3);
cVmMHz -> SetLogy();

G_vmmhz_filter->SetTitle("spectrum (after antenna, filter W/Hz)");
G_vmmhz_filter->GetHistogram()->SetXTitle("freq (Hz)");
G_vmmhz_filter->GetHistogram()->SetYTitle("mean power spectrum (W/Hz)");
G_vmmhz_filter->Draw("al");


cVmMHz -> cd(4);
cVmMHz -> SetLogy();

G_vmmhz_fft->SetTitle("spectrum (after antenna, filter, changed for fft W/Hz)");
G_vmmhz_fft->GetHistogram()->SetXTitle("freq (Hz)");
G_vmmhz_fft->GetHistogram()->SetYTitle("mean power spectrum (W/Hz)");
G_vmmhz_fft->Draw("al");


//cVmMHz -> cd(5)->SetLogy();


cVmMHz -> cd(5);
cVmMHz -> SetLogy();

G_Filter->SetTitle("Filter gain (unitless)");
G_Filter->GetHistogram()->SetXTitle("freq (Hz)");
G_Filter->GetHistogram()->SetYTitle("gain (Vo/Vi)");
G_Filter->Draw("al");


cVmMHz -> cd(6);

G_Vfft_noise_before->SetTitle("Noise spectrum (after filter W/Hz)");
G_Vfft_noise_before->GetHistogram()->SetXTitle("freq (Hz)");
G_Vfft_noise_before->GetHistogram()->SetYTitle("mean power spectrum (W/Hz)");
G_Vfft_noise_before->Draw("al");

cVmMHz -> cd(7);

G_Vfft_noise_after->SetTitle("Noise spectrum (after filter & random_rician W/Hz)");
G_Vfft_noise_after->GetHistogram()->SetXTitle("freq (Hz)");
G_Vfft_noise_after->GetHistogram()->SetYTitle("mean power spectrum (W/Hz)");
G_Vfft_noise_after->Draw("al");




cVmMHz -> Print("VmMHz_evt0.pdf");




} // end if there is ray sol




AraTree2->GetEvent(0);
for (int i=0; i<detector->params.number_of_stations; i++) {
    for (int j=0; j<detector->params.number_of_strings_station; j++) {
        for (int k=0; k<detector->params.number_of_antennas_string; k++) {
            cout<<"Rank of stations["<<i<<"].strings["<<j<<"].antennas["<<k<<"] = "<<report->stations[i].strings[j].antennas[k].Rank[0]<<endl;
            cout<<"PeakV of stations["<<i<<"].strings["<<j<<"].antennas["<<k<<"] = "<<report->stations[i].strings[j].antennas[k].PeakV[0]<<endl;
            cout<<"pathtime of stations["<<i<<"].strings["<<j<<"].antennas["<<k<<"] = "<<report->stations[i].strings[j].antennas[k].arrival_time[0]<<endl;
        }
    }
}


cout<<"noise Vfft_org : "<<report->Vfft_noise_org<<endl;

cout<<"total_mean_power_noise_before : "<<total_mean_power_noise_before<<endl;
cout<<"total_mean_power_noise_after : "<<total_mean_power_noise_after<<endl;
cout<<"total_mean_power_signal : "<<total_mean_power_signal<<endl;
//cout<<"Peak at evt0, 0110 : "<<report->stations[station_n].strings[string_n].antennas[antenna_n].PeakV[0]<<endl;


}


