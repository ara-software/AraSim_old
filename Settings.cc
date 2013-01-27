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

Settings::~Settings() {
    //default destructor
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

  TIMESTEP=(1./2.6)*1.E-9;  // default, same with icemc (in sec)

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

  TRIG_ANALYSIS_MODE=0;    // default : 0, signal + noise

  TRIG_TIMEOUT=1.E-6;       // default : 1us

  TRIG_WINDOW=2.5E-7;       // default : 250 ns

  NOISE_EVENTS=1000;        // default : 1000 events

  DATA_SAVE_MODE=0;         // default : 0 (full mode)

  N_TRIG=3;                 // default : 3 (3 out of all channels in a station)

  RANDOM_MODE=1;            // default : 1 (seed is unique in time/space)

  BORE_HOLE_ANTENNA_LAYOUT=0;   // default : 0 (VHVH)

  WRITE_ALL_EVENTS=0; //default : 0 (writes only globally triggered events)

  RAYSOL_RANGE=3000; // default : 3000 m

  PICK_POSNU_DEPTH=0;     //default : 0 pick posnu depth from 0 to ice depth

  MAX_POSNU_DEPTH=200.;     // default : 200m depth max

  NNU_THIS_THETA=0;         // default : nnu angle pure random

  NNU_THETA=0.785;          // default : nnu theta : 45 deg

  NNU_D_THETA=0.0873;       // default : nnu d_theta : 5 deg

    
    CALPULSER_ON=0; // default : calpulsers off
    
    TESTBED_ON=0; // default : 0 stations[0] is ARA1 not Testbed
    
    READGEOM=0; // default : 0 : use idealized geometry and do not read in from sqlite database
    
    V_MIMIC_MODE = 1; // default : 1 - write out window for non-triggered events that begins with the last bin in the trigger window
    
    USE_INSTALLED_TRIGGER_SETTINGS = 0; // default : 0 - use idealized settings for the trigger
    
    NUM_INSTALLED_STATIONS = 2;
    
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
              else if (label == "TRIG_ANALYSIS_MODE") {
                  TRIG_ANALYSIS_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
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
              else if (label == "N_TRIG") {
                  N_TRIG = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "RANDOM_MODE") {
                  RANDOM_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "BORE_HOLE_ANTENNA_LAYOUT") {
                  BORE_HOLE_ANTENNA_LAYOUT = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "RAYSOL_RANGE") {
                  RAYSOL_RANGE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "CALPULSER_ON") {
                  CALPULSER_ON = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TESTBED_ON") {
                  TESTBED_ON = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "READGEOM") {
                  cout << "Read in READGEOM" << endl;
                  READGEOM = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "PICK_POSNU_DEPTH") {
                  PICK_POSNU_DEPTH = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "MAX_POSNU_DEPTH") {
                  MAX_POSNU_DEPTH = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NNU_THIS_THETA") {
                  NNU_THIS_THETA = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NNU_THETA") {
                  NNU_THETA = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NNU_D_THETA") {
                  NNU_D_THETA = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "WRITE_ALL_EVENTS") {
                  WRITE_ALL_EVENTS = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "V_MIMIC_MODE") {
                  V_MIMIC_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "USE_INSTALLED_TRIGGER_SETTINGS") {
                  USE_INSTALLED_TRIGGER_SETTINGS = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NUM_INSTALLED_STATIONS") {
                  NUM_INSTALLED_STATIONS = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }              
          }
      }
      setFile.close();
  }
  else cout<<"Unable to open "<<setupfile<<" file!"<<endl;
  return;
}




