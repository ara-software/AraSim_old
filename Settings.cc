#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "Settings.h"

ClassImp(Settings);

using namespace std;

Settings::Settings() {
    Initialize();
//    ReadFile();

}

void Settings::Initialize() {


// below : values from icemc Settings class
  NDISCONES_PASS=3;

  DEBUG=false;                   // debugging option
outputdir="outputs"; // directory where outputs go
 FREQ_LOW_SEAVEYS=200.E6;
 FREQ_HIGH_SEAVEYS=1200.E6;
 BW_SEAVEYS=FREQ_HIGH_SEAVEYS-FREQ_LOW_SEAVEYS;
 SIGMAPARAM=1;  // Connolly et al. 2011 default cross section parametrization
 SIGMA_FACTOR=1.;   // default sigma factor : 1
 YPARAM=1;  // Connolly et al. 2011 default y parametrization
 UNBIASED_SELECTION=1.; // (0) pick neutrino interaction in the ice and neutrino from any direction or (1) choose neutrino interaction point in the horizon on the balloon in the ice and neutrino direction on the cerenkov cone

// end of values from icemc



  NNU=100;

  // NEED TO FIGURE OUT A GOOD WAY TO READ THIS IN AND STORE THEM.
  // INPUT FILE AGAIN?  SOMETHING ELSE?
  //These were moved here from IceModel under the new compilation scheme
  ICE_MODEL=0; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP.
  NOFZ=1; // 1=depth dependent index of refraction,0=off
  CONSTANTCRUST=0; // set crust density and thickness to constant values.
  CONSTANTICETHICKNESS=0; // set ice thickness to constant value
  FIXEDELEVATION=0; // fix the elevation to the thickness of ice.
  MOOREBAY=0; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data
  
  EXPONENT=19.; // 10^19 eV neutrinos only

  DETECTOR=1;   //ARA layout with small number of stations

  INTERACTION_MODE=1;   //PickNear mode

  POSNU_RADIUS=2000;    //radius for PickNear method

  WHICHPARAMETERIZATION=0;  //

  SIMULATION_MODE=0;    // default freq domain simulation

  EVENT_TYPE=0;         // default neutrino only events

  WAVE_TYPE=0;          // default wave type : plane wave (inside the ice)

  LPM=1;                //default : enable LPM effect

  SECONDARIES=1;        //default : enable secondary interactions

  TAUDECAY=1;           //default : let taudecay as secondary interactions

  TIMESTEP=(1./2.6)*1.E-9;  // default, same with icemc

  PHASE=90.;            // default : 90 deg phase (it means all imaginary values)

  NFOUR=1024;           // default : 1024, same as in icemc
    
  NOISE=0;              // degault : 0, flat thermal noise

  ATMOSPHERE=1;         // default : 1, include atmosphere

  POWERTHRESHOLD=-4.41; // default : -4.41 (same as icemc).

  MAXT_DIODE=70.E-9;    // default : 70 ns

  IDELAYBEFOREPEAK_DIODE=(int)(13.E-9 / TIMESTEP);    // default : 13.e-9/TIMESTEP = 33

  IWINDOW_DIODE=(int)(4.E-9 / TIMESTEP);           // default : 4.e-9 / TIMESTEP = 10

  DATA_BIN_SIZE=16384;   // default : 16384

  NOISE_TEMP=325.;      // default : 325 K

  PURE_NOISE_ANALYSIS=0;    // default : 0, no pure noise analysis

  TRIG_TIMEOUT=1.E-6;       // default : 1us

  TRIG_WINDOW=2.5E-7;       // default : 250 ns

  NOISE_EVENTS=1000;        // default : 1000 events

  DATA_SAVE_MODE=0;         // default : 0 (full mode)

}

void Settings::ReadFile(string setupfile) {

  ifstream setFile (setupfile.c_str());
  
  string line, label;

  if ( setFile.is_open() ) {
      while (setFile.good() ) {
          getline (setFile,line);

          if (line[0] != "/"[0]) {
              label = line.substr(0, line.find_first_of("="));
              
              if (label == "NNU") {
                  NNU = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "ICE_MODEL") {
                  ICE_MODEL = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NOFZ") {
                  NOFZ = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "CONSTANTCRUST") {
                  CONSTANTCRUST = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "CONSTANTICETHICKNESS") {
                  CONSTANTICETHICKNESS = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "FIXEDELEVATION") {
                  FIXEDELEVATION = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "MOOREBAY") {
                  MOOREBAY = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "EXPONENT") {
                  EXPONENT = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "DETECTOR") {
                  DETECTOR = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "INTERACTION_MODE") {
                  INTERACTION_MODE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "POSNU_RADIUS") {
                  POSNU_RADIUS = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "WHICHPARAMETERIZATION") {
                  WHICHPARAMETERIZATION = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "SIMULATION_MODE") {
                  SIMULATION_MODE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "EVENT_TYPE") {
                  EVENT_TYPE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "WAVE_TYPE") {
                  WAVE_TYPE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "LPM") {
                  LPM = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "SECONDARIES") {
                  SECONDARIES = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TAUDECAY") {
                  TAUDECAY = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TIMESTEP") {
                  TIMESTEP = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "PHASE") {
                  PHASE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NFOUR") {
                  NFOUR = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NOISE") {
                  NOISE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "ATMOSPHERE") {
                  ATMOSPHERE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "POWERTHRESHOLD") {
                  POWERTHRESHOLD = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "DATA_BIN_SIZE") {
                  DATA_BIN_SIZE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NOISE_TEMP") {
                  NOISE_TEMP = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "PURE_NOISE_ANALYSIS") {
                  PURE_NOISE_ANALYSIS = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TRIG_TIMEOUT") {
                  TRIG_TIMEOUT = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TRIG_WINDOW") {
                  TRIG_WINDOW = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NOISE_EVENTS") {
                  NOISE_EVENTS = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "DATA_SAVE_MODE") {
                  DATA_SAVE_MODE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
          }
      }
      setFile.close();
  }
  else cout<<"Unable to open setup.txt file!"<<endl;
  return;
}




