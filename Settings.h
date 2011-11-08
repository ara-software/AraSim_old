#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <string>
#include "TObject.h"

using namespace std;
//using std::string;

class Settings 
{
    protected:

    public:

        Settings();
        void Initialize();
        void ReadFile(string setupfile);

        int NNU;

  // NEED TO FIGURE OUT A GOOD WAY TO READ THIS IN AND STORE THEM.
  // INPUT FILE AGAIN?  SOMETHING ELSE?
  //These were moved here from IceModel under the new compilation scheme
        int ICE_MODEL; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP.
        int NOFZ; // 1=depth dependent index of refraction,0=off
        int CONSTANTCRUST; // set crust density and thickness to constant values.
        int CONSTANTICETHICKNESS; // set ice thickness to constant value
        int FIXEDELEVATION; // fix the elevation to the thickness of ice.
        int MOOREBAY; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data
        double EXPONENT; // 10^19 eV neutrinos only

        int DETECTOR;   // choose detector layout


 // below : values from icemc


  int UNBIASED_SELECTION;
 int WHICH; // which payload to use 0=Anita-lite,1=Ross,2=Smex,3=make your own
  int CYLINDRICALSYMMETRY; // is it cylindrically symmetric =1 if which=1,2, =0 if which=0
  // if which=3 then 0 or 1  
  double SIGMA_FACTOR; // factor to multiply cross section by for error analysis
  int SIGMAPARAM; // 0=Reno, 1=Connolly et al. 2011 for cross section parametrization
int YPARAM; // 0=Reno, 1=Connolly et al. 2011 for cross section parametrization
 int SIGNAL_FLUCT;  // 1=add noise fluctuation to signal or 0=do not
  int TRIGGERSCHEME;  // frequency domain voltage, frequency domain energy, time domain diode integration
  int ZEROSIGNAL;  // zero the signal to see how many of our hits are noise hits
    int REMOVEPOLARIZATION; //Disable polarizations
  
  int EVENTSMAP;//whether draw the events distribution map
  
  int WHICHRAYS;  // how many rays to look at (1) direct only (2) direct and down-going.

// trigger
int LCPRCP; // 1 for circular polarization trigger, 0 for V and H
int JUSTVPOL; // 0 for both polarizations, 1 for just V polarization
// doesn't allow for both LCPRCP=1 and JUSTVPOL=1
//int FIFTHBAND; // 1 to include 0.2-1.2 GHz as a frequency band if JUSTVPOL==1
//int NFOLD=3;  // how many channels must pass the trigger - in old mechanism - only used for anita-lite
int NFOLD;  // how many channels must pass the trigger - in old mechanism - only used for anita-lite


//int CHMASKING=1; // whether or not to include channel masking
//int PHIMASKING=1; // whether or not to include phi masking  
int CHMASKING; // whether or not to include channel masking
int PHIMASKING; // whether or not to include phi masking  

//int NLAYERS=0;
//int NANTENNAS=0;

int NLAYERS;
int NANTENNAS;

/* int ONLYFINAL=1; // only write to final histogram */
/* int HIST_MAX_ENTRIES=10000; //maximum number of events to put in histograms */
/* int HIST=1;          //write to histograms  */

int ONLYFINAL; // only write to final histogram
int HIST_MAX_ENTRIES; //maximum number of events to put in histograms
int HIST;          //write to histograms 
double BW; // BANDWIDTH
//int DISCONES=1; // whether or not to use discones
int DISCONES; // whether or not to use discones

//double NDISCONES_PASS=3; // number of discones needed to pass
double NDISCONES_PASS; // number of discones needed to pass

int BORESIGHTS; // whether to loop over boresights
int SLAC; // whether or not we are simulating the slac run
double SLACSLOPE; // slope of the ice
double SLACICELENGTH;  // length of the block of ice
double SLAC_HORIZDIST; // horizontal distance from interaction to center of payload at slac beam test
double SLAC_DEPTH; // vertical depth of interaction at slac beam test
double SLAC_HORIZ_DEPTH; // horizontal depth of interaction at slac

int ROUGHNESS; // include effects of surface roughness
int FIRN; // whether or not to include the firn

//int SLOPEY=1; // 1=slopeyness on, 0=slopeyness off
//double SLOPEYSIZE=0.012; // This determines size of the slopeyness (0.10=5.4, 0.20=7.4 deg mean)

int SLOPEY; // 1=slopeyness on, 0=slopeyness off
double SLOPEYSIZE; // This determines size of the slopeyness (0.10=5.4, 0.20=7.4 deg mean)

 bool DEBUG;
string outputdir; // directory where outputs go

//double THERMALNOISE_FACTOR=1.0; // factor to multiply thermal noise for error analysis
double THERMALNOISE_FACTOR; // factor to multiply thermal noise for error analysis

//double FREQ_LOW_SEAVEYS=200.E6; // min frequency for seaveys
//const double FREQ_HIGH_SEAVEYS=1200.E6; // max frequency for seaveys

double FREQ_LOW_SEAVEYS; // min frequency for seaveys
double FREQ_HIGH_SEAVEYS; // max frequency for seaveys
 double BW_SEAVEYS;
 //int FORSECKEL=1; // Make array of strength of signal across frequencies for different viewing angles.
int FORSECKEL; // Make array of strength of signal across frequencies for different viewing angles.

double ROUGHSIZE; // roughness size


/* int ICE_MODEL=0; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP. */
/* int NOFZ=1; // 1=depth dependent index of refraction,0=off */
/* int CONSTANTCRUST=0; // set crust density and thickness to constant values. */
/* int CONSTANTICETHICKNESS=0; // set ice thickness to constant value */
/* int FIXEDELEVATION=0; // fix the elevation to the thickness of ice. */
/* int MOOREBAY=0; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data */
/* int USEPOSITIONWEIGHTS=1;// whether or not to restrict the neutrino position so it is within the horizon of the balloon */
/* int WRITE_FILE=0; //Select whether or not to write a new input file for CreateHorizons */

int USEPOSITIONWEIGHTS;// whether or not to restrict the neutrino position so it is within the horizon of the balloon
int WRITE_FILE; //Select whether or not to write a new input file for CreateHorizons

int MINRAY;
int MAXRAY;

int horizontal_banana_points;
  int vertical_banana_points;



 // end of values from icemc


  ClassDef(Settings,1);


};
#endif

