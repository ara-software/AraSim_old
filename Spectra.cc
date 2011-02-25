#include "Spectra.h"
#include "Tools.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

Spectra::Spectra(int EXPONENT) {



  double Emuons[12]; // E^2 dN/dE/dA/dt for neutrinos that are produced as muon neutrinos or muon antineutrinos.
  double Eelectrons[12];// E^2 dN/dE/dA/dt for neutrinos that are produced as electron neutrinos or muon antineutrinos.
  
  for (int i=0;i<12;i++) {
    energy[i]=16.+((double)i)/2.;
    Emuons[i]=-30.;
    Eelectrons[i]=-30.;
  } //for
  if(EXPONENT>100 && EXPONENT<200)//using Fenfang-digitized fluxes
    {
      switch(EXPONENT)
	{
	case 101:
	  GetFlux("ahlers.dat");
	  break;
	case 102:
	  GetFlux("allard.dat");
	  break;
	case 103:
	  GetFlux("ave_max.dat");
	  break;
	case 104:
	  GetFlux("ave_min.dat");
	  break;
	case 105:
	  GetFlux("essfig9.dat");
	  break;
	case 106:
	  GetFlux("essbaseline.dat");
	  break;
	case 107:
	  GetFlux("ess_n0.dat");
	  break;
	case 108:
	  GetFlux("kkss_envo.dat");
	  break;
	case 109:
	  GetFlux("gzk_peter.dat");
	  break;
	case 110:
	  GetFlux("waxgzk.dat");
	  break;
	case 111:
	  GetFlux("e-2.dat");
	  break;
	case 112:
	  GetFlux("yuksel_grb.dat");
	  break;
	case 113:
	  GetFlux("yuksel_qso.dat");
	  break;
	case 114:
	  GetFlux("yuksel_sfh.dat");
	  break;
	case 115:
	  GetFlux("berezinsky_saturate.dat");
	  break;
	default:
	  cout<<"Error: Wrong input of EXPONENT"<< endl;
	  
	}
    }
  
  
  
  else if (EXPONENT==0) {
    // what I was using previously, from upper curve of ANITA proposal
    //    E2dNdEdAdt[0]=-9.6; // 16.
    //    E2dNdEdAdt[1]=-8.9; // 16.5
    //    E2dNdEdAdt[2]=-8.1; // 17.
    //    E2dNdEdAdt[3]=-7.5; // 17.5
    //    E2dNdEdAdt[4]=-7.2; // 18.
    //    E2dNdEdAdt[5]=-6.8; // 18.5
    //    E2dNdEdAdt[6]=-6.7; // 19
    //    E2dNdEdAdt[7]=-6.8; // 19.5
    //    E2dNdEdAdt[8]=-7.2; // 20.
    //    E2dNdEdAdt[9]=-7.5; // 20.5
    //    E2dNdEdAdt[10]=-8.2; // 21.0
    //    E2dNdEdAdt[11]=-9.1; // 21.5
    
    // electron component of Figure 4 of ES&S
    // astro-ph/0101216
    Eelectrons[0]=-17.2; // 16.
    Eelectrons[1]=-17.35; // 16.5
    Eelectrons[2]=-17.2; // 17.
    Eelectrons[3]=-17.1; // 17.5
    Eelectrons[4]=-17.2; // 18.
    Eelectrons[5]=-17.5; // 18.5
    Eelectrons[6]=-18.0; // 19
    Eelectrons[7]=-18.5; // 19.5
    Eelectrons[8]=-19.4; // 20.
    Eelectrons[9]=-30.; // 20.5 punt
    Eelectrons[10]=-30.; // 21.0 punt
    Eelectrons[11]=-30.; // 21.5 punt
    
    // muon component of Figure 4 of ES&S
    // astro-ph/0101216
    //    Emuons[0]=-17.8; // 16.
    //    Emuons[1]=-17.4; // 16.5
    //    Emuons[2]=-17.; // 17.
    //    Emuons[3]=-16.75; // 17.5
    //    Emuons[4]=-16.9; // 18.
    //    Emuons[5]=-17.2; // 18.5
    //    Emuons[6]=-17.7; // 19
    //    Emuons[7]=-18.3; // 19.5
    //    Emuons[8]=-19.1; // 20.
    //    Emuons[9]=-30.; // 20.5 punt
    //    Emuons[10]=-30.; // 21.0 punt
    //    Emuons[11]=-30.; // 21.5 punt
    
    // lower curve of Figure 9 of ES&S
    // astro-ph/0101216
       Emuons[0]=-17.1;  //16.
       Emuons[1]=-16.6;  //16.5
       Emuons[2]=-16.3;  //17.
       Emuons[3]=-16.2; // 17.5
       Emuons[4]=-16.4; // 18.
       Emuons[5]=-16.7; // 18.5
       Emuons[6]=-17.3; // 19
       Emuons[7]=-17.95; // 19.5
       Emuons[8]=-18.85; // 20.
       Emuons[9]=-19.9; // 20.5 punt
       Emuons[10]=-30.; // 21.0 punt
       Emuons[11]=-30.; // 21.5 punt
    
    
//     // upper curve in Figure 9 of ES&S
//     // astro-ph/0101216
//     Emuons[0]=-16.85;  //16.
//     Emuons[1]=-16.4;  //16.5
//     Emuons[2]=-16.05;  //17.
//     Emuons[3]=-16.; // 17.5
//     Emuons[4]=-16.15; // 18.
//     Emuons[5]=-16.5; // 18.5
//     Emuons[6]=-17.1; // 19
//     Emuons[7]=-17.7; // 19.5
//     Emuons[8]=-18.65; // 20.
//     Emuons[9]=-19.75; // 20.5 punt
//     Emuons[10]=-30.; // 21.0 punt
//     Emuons[11]=-30.; // 21.5 punt
    
  
    for (int i=0;i<12;i++) {
      EdNdEdAdt[i]=pow(10.,Eelectrons[i])+pow(10.,Emuons[i]);
      //      cout << "EdNdEdAdt is " << EdNdEdAdt[i] << "\n";
      //      cout << "energy is " << energy[i] << "\n";
      //      cout << "log(EdNdEdAdt) is " << log10(EdNdEdAdt[i]) << "\n";
    }
     

  } // end EXPONENT=0, baseline ES&S

  else if (EXPONENT==1) { // E^-1 spectrum
    for (int i=0;i<12;i++) {
      energy[i]=17.5+((double)i)/2.;
    } //for
    
    for (int i=0;i<12;i++) {
      EdNdEdAdt[i]=1.;
    } //for
  }
  else if (EXPONENT==2) { // E^-2 spectrum

    double E2dNdEdAdt[12]; //log(brightness)

    for (int i=0;i<12;i++) {
      energy[i]=17.5+((double)i)/2.;
    } //for
  
    for (int i=0;i<12;i++) {
      E2dNdEdAdt[i]=1.; 
    } //for

    for (int i=0;i<12;i++) {
      EdNdEdAdt[i]=pow(10,E2dNdEdAdt[i]-(energy[i]-9.)); // convert E^2 to E, where E is in GeV and the energy array is the exponent of the energy in eV
    }  //for
 
  } // end E^-2
  else if (EXPONENT==3) {
    
    double E2dNdEdAdt[12]; //log(brightness)
    
    for (int i=0;i<12;i++) {
      energy[i]=17.5+((double)i)/2.;
    } //for
    
    for (int i=0;i<12;i++) {
      E2dNdEdAdt[i]=-(energy[i]-9); 
    } //for
    
    for (int i=0;i<12;i++) {
      EdNdEdAdt[i]=pow(10,E2dNdEdAdt[i]-(energy[i]-9.));
    }  //for
    
  } //E^-3
  else if (EXPONENT==4) {

    double E2dNdEdAdt[12]; //log(brightness)
    
    for (int i=0;i<12;i++) {
      energy[i]=17.5+((double)i)/2.;
    } //for
    
    for (int i=0;i<12;i++) {
    E2dNdEdAdt[i]=-2.*(energy[i]-9); 
    } //for
    
    for (int i=0;i<12;i++) {
      EdNdEdAdt[i]=pow(10,E2dNdEdAdt[i]-(energy[i]-9.));
    }  //for
    
  } //E^-4

  // fenfang's getgzk
  else if (EXPONENT==5) {
    
    //double E2dNdEdAdt[12]; //log(brightness)
    
    double Emuons[12]; // E^2 dN/dE/dA/dt for neutrinos that are produced as muon neutrinos or muon antineutrinos.
    double Eelectrons[12];// E^2 dN/dE/dA/dt for neutrinos that are produced as electron neutrinos or muon antineutrinos.
    double emuratio[12];
    for (int i=0;i<12;i++) {
      energy[i]=16.+((double)i)/2.;
      Emuons[i]=-30.;
      Eelectrons[i]=-30.;
 } //for
    // what I was using previously, from upper curve of ANITA proposal
    //    E2dNdEdAdt[0]=-9.6; // 16.
    //    E2dNdEdAdt[1]=-8.9; // 16.5
    //    E2dNdEdAdt[2]=-8.1; // 17.
    //    E2dNdEdAdt[3]=-7.5; // 17.5
    //    E2dNdEdAdt[4]=-7.2; // 18.
    //    E2dNdEdAdt[5]=-6.8; // 18.5
    //    E2dNdEdAdt[6]=-6.7; // 19
    //    E2dNdEdAdt[7]=-6.8; // 19.5
    //    E2dNdEdAdt[8]=-7.2; // 20.
    //    E2dNdEdAdt[9]=-7.5; // 20.5
    //    E2dNdEdAdt[10]=-8.2; // 21.0
    //    E2dNdEdAdt[11]=-9.1; // 21.5
    
    // electron component of Figure 4 of ES&S
    // astro-ph/0101216
    Eelectrons[0]=-17.2; // 16.
    Eelectrons[1]=-17.35; // 16.5
    Eelectrons[2]=-17.2; // 17.
    Eelectrons[3]=-17.1; // 17.5
    Eelectrons[4]=-17.2; // 18.
    Eelectrons[5]=-17.5; // 18.5
    Eelectrons[6]=-18.0; // 19
    Eelectrons[7]=-18.5; // 19.5
    Eelectrons[8]=-19.4; // 20.
    Eelectrons[9]=-30.; // 20.5 punt
    Eelectrons[10]=-30.; // 21.0 punt
    Eelectrons[11]=-30.; // 21.5 punt
    
    // muon component of Figure 4 of ES&S
    // astro-ph/0101216
    Emuons[0]=-17.8; // 16.
    Emuons[1]=-17.4; // 16.5
    Emuons[2]=-17.; // 17.
    Emuons[3]=-16.75; // 17.5
    Emuons[4]=-16.9; // 18.
    Emuons[5]=-17.2; // 18.5
    Emuons[6]=-17.7; // 19
    Emuons[7]=-18.3; // 19.5
    Emuons[8]=-19.1; // 20.
    Emuons[9]=-30.; // 20.5 punt
    Emuons[10]=-30.; // 21.0 punt
    Emuons[11]=-30.; // 21.5 punt
    
    for(int i=0;i<12;i++)
      emuratio[i]=Eelectrons[i]/Emuons[i];
    // lower curve of Figure 9 of ES&S
    // astro-ph/0101216
    //    Emuons[0]=-17.1;  //16.
    //    Emuons[1]=-16.6;  //16.5
    //    Emuons[2]=-16.3;  //17.
    //    Emuons[3]=-16.2; // 17.5
    //    Emuons[4]=-16.4; // 18.
    //    Emuons[5]=-16.7; // 18.5
    //    Emuons[6]=-17.3; // 19
    //    Emuons[7]=-17.95; // 19.5
    //    Emuons[8]=-18.85; // 20.
    //    Emuons[9]=-19.9; // 20.5 punt
    //    Emuons[10]=-30.; // 21.0 punt
    //    Emuons[11]=-30.; // 21.5 punt
    
    
    // upper curve in Figure 9 of ES&S
    // astro-ph/0101216
    Emuons[0]=-16.85;  //16.
    Emuons[1]=-16.4;  //16.5
    Emuons[2]=-16.05;  //17.
    Emuons[3]=-16.; // 17.5
    Emuons[4]=-16.15; // 18.
    Emuons[5]=-16.5; // 18.5
    Emuons[6]=-17.1; // 19
    Emuons[7]=-17.7; // 19.5
    Emuons[8]=-18.65; // 20.
    Emuons[9]=-19.75; // 20.5 punt
    Emuons[10]=-30.; // 21.0 punt
    Emuons[11]=-30.; // 21.5 punt
    
    
    for(int i=0;i<12;i++)
      Eelectrons[i]=Emuons[i]*emuratio[i];
  
    for (int i=0;i<12;i++) {
      EdNdEdAdt[i]=pow(10.,Eelectrons[i])+pow(10.,Emuons[i]); 
    }
    
  } // end if cosmological constant
// Incident neutrinos
  else if (EXPONENT==6) {

    for (int i=0;i<12;i++) {
      energy[i]=12.+((double)i);
    } //for
    
    EdNdEdAdt[0]=-17.1; // 12.
    EdNdEdAdt[1]=-16.5; // 13.
    EdNdEdAdt[2]=-16.1; // 14.
    EdNdEdAdt[3]=-16.5; // 15.
    EdNdEdAdt[4]=-17.3; // 16.
    EdNdEdAdt[5]=-16.8; // 17.
    EdNdEdAdt[6]=-17.2; // 18.
    EdNdEdAdt[7]=-19.; // 19.
    EdNdEdAdt[8]=-0.; // 20.
    EdNdEdAdt[9]=-0.; // 21.
    EdNdEdAdt[10]=-0.; // 22.
    EdNdEdAdt[11]=-0.; // 23.
    
    for (int i=0;i<12;i++) 
      {
	EdNdEdAdt[i]=pow(10,EdNdEdAdt[i]);
	
      } //for
        
    // now throw at a dartboard.
    
    //   double thisenergy=16.;
    //   double thisflux=2.;
    //   double max=1.;
    //   int energybin=0;
    //   while(thisflux>max) {
    //     // pick an energy  
    //     thisenergy=Rand3.Rndm()*(Tools::dMax(energy,8)-GetMin(energy,8));
    //     energybin=(int)(thisenergy);
    //     max=EdNdEdAdt[energybin];
    //     thisflux=Rand3.Rndm(); 
    //   } //while
    
    //   for (int i=0;i<8;i++) {
    //     EdNdEdAdt[i]=EdNdEdAdt[i]*maxflux;
    //   } //for 
    
    //   return pow(10.,thisenergy+GetMin(energy,8));
  } //GetGZKIron
  else if (EXPONENT==200) { // flat spectrum - for generating a very hard spectrum
    cout << "I'm in exponent==200.\n";
    for (int i=0;i<12;i++) {
      energy[i]=17.5+((double)i)/2.;
    } //for
    double E2dNdEdAdt[12];
    for (int i=0;i<12;i++) {
      E2dNdEdAdt[i]=1.; 
    } //for
    // start with E^-2
    for (int i=0;i<12;i++) {
      EdNdEdAdt[i]=pow(10,E2dNdEdAdt[i]-(energy[i]-9.)); // convert E^2 to E, where E is in GeV and the energy array is the exponent of the energy in eV
      
      EdNdEdAdt[i]*=pow(2,i);
    }  //for


  }
  // find the max so we can normalise it
  maxflux=Tools::dMax(EdNdEdAdt,12);

}
double  Spectra::GetNuEnergy() {
  
  TRandom3 Rand3;

  double thisenergy=16.; // arbitrary initialisation
  double thisflux=2.; // initialise higher than max
  double max=1.;
  int energybin=0; // arbitrary initialisation
  for (int i=0;i<12;i++) {
    //cout << "energy, EdNdEdAdt, flux/max are " << energy[i] << " " << EdNdEdAdt[i] << " " << EdNdEdAdt[i]/maxflux << "\n";
  }
  while(thisflux>max) {
    // pick an energy  
    thisenergy=Rand3.Rndm()*(Tools::dMax(energy,12)-Tools::dMin(energy,12)); // pick energy at random between the highest and lowest
    // the energy array is actually filled with energy exponents 
    // and thisenergy starts from 0 so it has an offset
    energybin=(int)(thisenergy/0.5); // this assumes the bins are in increments of half an order of magnitude
    max=EdNdEdAdt[energybin]/maxflux; // this is the maximum the normalized flux can be in this bin, always less than 1
    thisflux=Rand3.Rndm(); // pick the flux at random between 0 and 1, if it's less than max it's a winner
  } //while
  //cout << "energy is " << pow(10.,thisenergy+GetMin(energy,12)) << "\n";
  return pow(10.,thisenergy+Tools::dMin(energy,12));
	
} //Pick Neutrino Energy

void Spectra::GetFlux(string filename)
{ifstream influx(("./fluxes/"+filename).c_str());
 const int NLINES=11;
 double flux[NLINES];//E2dNdE GeV cm^-2 s^-1 sr^-1
 double junk;
 cout<<"We are using "<<filename.c_str()<<" as the flux data."<<endl;
 for(int i=0;i<NLINES;i++)
   {influx>>junk>>flux[i];
   }
 flux[11]=-15.;//I don't have the flux value for 10^21.5 eV, so just make it small
 for(int i=0;i<12;i++)
   {
     EdNdEdAdt[i]=pow(10.,(flux[i]+9.-energy[i]));
     
   }
 // maxflux=Tools::dMax(EdNdEdAdt,12);
 // cout<<maxflux<<endl;

}
