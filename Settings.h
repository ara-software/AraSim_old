#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <string>
#include "TObject.h"


class Detector;

using namespace std;
//using std::string;

class Settings 
{
    protected:

    public:

        Settings();
        ~Settings();
        void Initialize();
        void ReadFile(string setupfile);


        int CheckCompatibilities(Detector *detector);// check if settings are not compatible to each other


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

        int INTERACTION_MODE;   // method to choose interaction point posnu. 0 : PickUnbiased, 1 : PickNear

        double POSNU_RADIUS;    //PickNear radius in meter

        int WHICHPARAMETERIZATION;  //

        int SIMULATION_MODE;    // 0 : old freq domain mode, 1: new time domain mode

        int EVENT_TYPE;         // 0 : neutrino only events,  1 : blackhole evnet? ... etc

        int WAVE_TYPE;          // 0 : plane wave,  1 : spherical wave

        int LPM;                // 0 : no LPM effect, 1 : with LPM effect(default)

        int SECONDARIES;        // 0 : no secondary interaction, 1 : with secondary interactions

        int TAUDECAY;           // 0 : don't account tau decay as secondaries, 1 : let tau decay as secondaries

        double TIMESTEP;        // time step after fft. t domain bin width

        double PHASE;           // phase factor for fft. default 90 deg (in deg !!)

        int NFOUR;              // number of total bins for FFT. has to be power of 2 values

        int NOISE;              // noise condition settings degault 0 ( : thermal flat noise), 1 : Rayleigh dist. fit for installed TestBed geom

        int ATMOSPHERE;         // include atmosphere 1, no 0

        double POWERTHRESHOLD;  // power threshold value. default -4.41 (same with icemc powerthreshold)

        double MAXT_DIODE;      // diode model max time, default : 70.e-9s

        int IDELAYBEFOREPEAK_DIODE; // diode model bins before the peak, default : (int) 13.e-9 / TIMESTEP = 33

        int IWINDOW_DIODE;            // diode model interesting bins after the peak, default : (int) 4.e-9 / TIMESTEP = 10

        int DATA_BIN_SIZE;          // bin size which mimic the data (time delay between antennas), default : 16384

        double NOISE_TEMP;          // noise temperature (default : 325 K = 230K (ice) + 95K (receiver), from Peter's spreadsheet)

        int TRIG_ANALYSIS_MODE; // trigger mode 0 : signal + noise, 1 : only pure signal, 2 : only pure noise, default : 0

        double TRIG_TIMEOUT;    // time out after the trigger (we have to wait this amount of time to trig next event), default : 1us = 1.E-6
        double TRIG_WINDOW;     // coincidence time window for trigger default : 250 ns

        int NOISE_EVENTS;       // number of noise events which will be stored in Trigger class for later use. This will also used to calculate mean, rms noise (with diode convlv). default : 1000

        int DATA_SAVE_MODE;     // 0 : save all information which are generated during the processes. 1 : light mode, remove most of data before the final global trigger (most of data except geometric info, final data V_mimic will be removed). 2 : light mode 2, except for physics information (energy, geometric info), all waveform information are removed (include noise waveforms in trigger class)

        int N_TRIG;         // number of coincidence channels triggered in TRIG_WINDOW default : 3

        int RANDOM_MODE;    // 0 : same random number generate, 1 : TRandom3(0) used. purely randomly generated (seed is guaranteed to be unique in space and time)

        int BORE_HOLE_ANTENNA_LAYOUT;   // 0 = (V-H-V-H), 1 = (V-H-V), 2 = (V-H-V-V), 3 = (V-H-H-H), 4 = (V-H-H) default : 0
    
        int WRITE_ALL_EVENTS; // 0 only write globally triggered events, 1 Write all event events including events that are not globally triggered
        // Note: NNU is the number of neutrinos that have been thrown in total, not just globally triggered events
        // When writing all events, the waveform stored in UsefulAraSimEvent->VoltsRF[] is just the untriggered, noiseless waveform of the initial signal, and it has not propagated to the antenna yet.
        // Hong added : 2 don't write any info to UsefulAraSimEvent (to reduce output root file)

        double RAYSOL_RANGE;    // direct distance limit to do raysolver. If distance between posnu and antenna is bigger than RAYSOL_RANGE, AraSim will not try RaySol for that event. Default : 3000 m

        int CALPULSER_ON; // 0: no calpulsers in event list, 1: only throws calpulser 1 events, 2: only throws calpulser 2 events, 3: throws both calpulser 1 and 2 events alternating between them, 4: throws calpulser 1 and 2 events integrated with the simulated data (not yet implemented)
    
        int TESTBED_ON;
    
        int READGEOM;
    
        int V_MIMIC_MODE; // default : 0 - write out all chs where global triggered bin is center of the window
                        // 1 - same as above 0 mode but apply TestBed ch delay - average BH ch delay
                        // 2 - same as above 0 mode but apply TestBed ch delay - average BH ch delay + additional delay to match with actual TestBed data waveforms
    
        int USE_INSTALLED_TRIGGER_SETTINGS; // default : 0 - use idealized settings for the trigger
        //other options:  1 - use trigger settings for installed stations, i.e. trigger window, etc.
    
        int NUM_INSTALLED_STATIONS; // the number of stations including the testbed that are in reality installed in the ice and have position and electronics information
    
        int PICK_POSNU_DEPTH;  // whether use MAX_POSNU_DEPTH or not. 0 : pick posnu depth full ice depth, 1 : pick posnu depth only MAX_POSNU_DEPTH

        double MAX_POSNU_DEPTH;  // maximum posnu depth when above PICK_POSNU_DEPTH=1

        int NNU_THIS_THETA;     // if nnu theta angle will be selected randomly from [0, PI] (default=0) or set nnu theta to near some angle (=1)

        double NNU_THETA;       // nnu theta when NNU_THIS_THETA=1

        double NNU_D_THETA;     // nnu theta variation from NNU_THETA, when NNU_THIS_THETA=1 case

        double CALPUL_OFFCONE_ANGLE;    // for calpulser events, what's the offcone angle value?

        double CALPUL_AMP;    // for calpulser events, how strong the calpulser waveforms?

        int TRIG_ONLY_BH_ON;    // if trigger will occur with all chs (0, default) or only borehole chs (1)

        int TRIG_THRES_MODE;    // if trigger threshold will use no specific offset value (0, default). or use data/threshold_offset.csv file values as a threshold offset for each chs

        int NOISE_TEMP_MODE;    // if using same noise temp for all chs (0, default), using diff temp for all chs (1), using diff temp for first 8 chs and share same temp for other chs (2)

        double CONST_MEANDIODE;
    
        double CONST_RMSDIODE;  // in case NOISE_TEMP_MODE = 1, just using this CONST_RMSDIODE value for threshold


        int USE_TESTBED_RFCM_ON;    // use RFCM measurement for testbed or not (default 0)
    
        double RFCM_OFFSET;  // if above USE_TESTBED_RFCM_ON = 1, we need RFCM attenuator factor cancel (default 80)

        int USE_MANUAL_GAINOFFSET; //if use gain offset file to read values or just use constant gain offset from setup file (default 0 : use file)
    
        double MANUAL_GAINOFFSET_VALUE; // gain offset value

        int NOISE_WAVEFORM_GENERATE_MODE; // default 0 : generate new noise waveforms for each evts, if you set other values, noise waveforms will be generated once and use them for all evts


        int USE_CH_GAINOFFSET; // default 0 : don't apply any individul channels' gain offset. 1 : apply ch gain offset by using data/preamp_gain_offset.csv file (only installed TestBed mode available)

        int GETCHORD_MODE; // which Getchord function to use. default 0 : old Getchord function (not correct tau weight, weight don't have ice inside interaction probability in it). 1 : new Getchord from icemc. This has new tau weight calculation and ice interaction probability applied to weight factor.

        int taumodes; // taumodes 1 : tau created in the rock

        int BH_ANT_SEP_DIST_ON; // if we are going to use separate bore hole antenna distance or not. By default it's 0 (don't use separate dist)

        int TRIG_MODE; // default 0 : if any antennas got passed N_TRIG or more, global trig. 1 : either Vpol or Hpol antennas got passed N_TRIG_V or N_TRIG_H respectively, global trig.

        int N_TRIG_V; // default : 3 (3 out of Vpolchannels in a station)
    
        int N_TRIG_H; // default : 3 (3 out of Hpol channels in a station)

        int FILL_TREE_MODE; // default 0 : fill tree for all events, 1 : fill tree only usable posnu events, 2 : fill tree only trigger passed events




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

