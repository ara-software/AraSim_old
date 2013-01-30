////////////////////////////////////////////////////////////////////////////////////////////////
//class Trigger:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef TRIGGER_H
#define TRIGGER_H

#include <vector>
#include "TObject.h"
#include <cstdlib>

using namespace std;

class Efficiencies;
class Detector;
class Settings;
class Report;

class Trigger {


 private:



 public:

           double Full_window[16][16384];  // test with array, not vector, diode response
           double Full_window_V[16][16384];  // test with array, not vector, voltage waveform


     double TIMESTEP;   // will copy from Detector class
     double maxt_diode; // will copy from Detector class
     int maxt_diode_bin; // will copy from Detector class
     int NFOUR;         // will copy from Detector class
     int DATA_BIN_SIZE; // will copy from settings class

     double V_noise_freqbin;    // thermal noise freq bin value


     double meandiode;
     double rmsdiode;
     double rmsvoltage;// rms voltage value without diode response

     vector < vector <double> > v_noise_timedomain;   // time domain noise waveform examples
     vector < vector <double> > v_noise_timedomain_diode; // time domain diode convlved noise waveforms examples

     double powerthreshold; // threshold for the trigger

     int iminbin;   // same with icemc trigger
     int imaxbin;
     
     Trigger();
     Trigger(Detector *detector, Settings *settings1);
     ~Trigger();

     void SetMeanRmsDiode(Settings *settings1, Detector *detector, Report *report);

     int CheckChannelsPass( vector <double> &V_total_diode);
     
     
     void myconvlv(vector <double> &data, const int NFOUR, vector <double> &fdiode, vector <double> &diodeconv);
     void myconvlv(double *data, const int NFOUR, vector <double> &fdiode, vector <double> &diodeconv);

     

     void myconvlv_half(vector <double> &data, const int NFOUR, vector <double> &fdiode, vector <double> &diodeconv);
     void myconvlv_half(double *data, const int NFOUR, vector <double> &fdiode, vector <double> &diodeconv);



           //double Full_window[16][16384];  // test with array, not vector
           //vector < vector <double> >  Full_window;  // test with array, not vector

     
     ClassDef(Trigger,1);

};

#endif //TRIGGER_H
