////////////////////////////////////////////////////////////////////////////////////////////////
//class Detector:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef DETECTOR_H
#define DETECTOR_H


//#include "TObject.h"
#include <vector>
#include "Trigger.h"
#include "Position.h"
#include "Vector.h"
//#include "IceModel.h"


using namespace std;

//#include "Trigger.h"
//class Event;
//class Efficiencies;
//class Position;

class Detector;
class IceModel;
class Settings;
class TF1;

//struct Parameters {
//class Parameters : public TObject {
class Parameters {

    private:

    public:
    int stations_per_side;           //number of stations on one side of pentagon.
    double station_spacing;          // space between stations.

    int antenna_orientation;        //antenna orientation setting.

    int number_of_stations;         //total stations
    int number_of_strings_station;  // strings per station
    int number_of_antennas_string;  //antennas per string
    int number_of_surfaces_station;    //surface antennas per station
    int number_of_channels; // number of channels for each regular (non-TestBed) station
    
    int number_of_strings_station_TB; //number of strings in the TestBed
    int number_of_antennas_string_TB; //number of antennas in each string in the Testbed
    int number_of_surfaces_station_TB; //number of surface stations for the TestBed
    int number_of_channels_TB; // number of channels for the TestBed
    
    int num_of_channels[2];
    
    int bore_hole_antenna_layout;   // bore hole antenna layout, 0 : VHVH, 1 : VHV, 2 : VHVV

    int number_of_strings;
    int number_of_antennas;

    double core_x;
    double core_y;
    
    int DeployedStations;
    
        
//    static const int freq_step = 60;
//    static const int ang_step = 2664;
//    static const double freq_width = 16.667;  // this value could be changed when Nec2 condition changes!
//    static const double freq_init = 83.333;  // this value could be changed when Nec2 condition changes!
//    static const int freq_width = 16.667;  // this value could be changed when Nec2 condition changes!
//    static const int freq_init = 83.333;  // this value could be changed when Nec2 condition changes!
//    static const float freq_width = 16.667;  // this value could be changed when Nec2 condition changes!
//    static const float freq_init = 83.333;  // this value could be changed when Nec2 condition changes!

//    int freq_step = 60;
//    int ang_step = 2664;
//    double freq_width = 16.667;  // this value could be changed when Nec2 condition changes!
//    double freq_init = 83.333;  // this value could be changed when Nec2 condition changes!

    int freq_step;  // this value will be obtained through settings class and copy to Detector->freq_step
    int ang_step;   // also copy to Detector->ang_step
    double freq_width;  // this value could be changed when Nec2 condition changes! copy to Detector->freq_width
    double freq_init;  // this value could be changed when Nec2 condition changes! copy to Detector->freq_init


    double TestBed_Ch_delay[16];
    int TestBed_Ch_delay_bin[16]; // in bin
    double TestBed_BH_Mean_delay;
    int TestBed_BH_Mean_delay_bin; // in bin

    double TestBed_WFtime_offset_ns; // waveform time offset to match TestBed data waveform


    ClassDef(Parameters,1);


};





//struct Surface_antenna : public Position, public TObject {
//class Surface_antenna : public Position, public TObject {
class Surface_antenna : public Position {

    public:

//    double x, y;
    int type;   // need to be defined
    int orient; //0 : facing x, 1 : y, 2 : -x, 3 : -y.

    double GetG(Detector *D, double freq, double theta, double phi);    // read gain value from Detector class

    ClassDef(Surface_antenna,1);

};


    
//struct Antenna : public Position, public TObject {
//class Antenna : public Position, public TObject {
class Antenna : public Position {
    public:

//    double z;
    int type;  // type 0 : v-pol (bicone), type 1 : h-pol (bowtie for testbed, QSC for ARA)
    int orient; // 0 : facing x, 1 : y, 2 : -x, 3 : -y.

    int DAQchan;    // DAQ channel type (from AraGeomTools). 0 : discone (BH chs), 1 : BAT (shallow, surf; not first 8 chs)

    int manual_delay;   // to fit the waveform to actual TestBed waveform, added manual delay time

    int manual_delay_bin;   // to fit the waveform to actual TestBed waveform, added manual delay bin
    
    double GetG(Detector *D, double freq, double theta, double phi);    // read gain value from Detector class

    //ClassDef(Antenna,1);
    ClassDef(Antenna,3);
    
};



//struct Antenna_string : public Position, public TObject {
//class Antenna_string : public Position, public TObject {
class Antenna_string : public Position {
//    double x, y;
//    int number_of_antennas;
    public:
    vector <Antenna> antennas;

    ClassDef(Antenna_string,1);
};

//struct ARA_station : public Position, public TObject {
//class ARA_station : public Position, public TObject {
class ARA_station : public Position {
//    double x, y;
    public:
    vector <Antenna_string> strings;
    vector <Surface_antenna> surfaces;
    int StationID;
    double TRIG_WINDOW; // in ns, the size of the trigger window used for the antennas
    int NFOUR; // 2 X nbins for readout waveform - for fourier tranform
    double TIMESTEP; // trigger and readout timestep
    double DATA_BIN_SIZE;
    
    ClassDef(ARA_station,1);
};



class InstalledStation {

    public:
    int nSurfaces;
    int nStrings;
    vector < int > surfaceChannels;
    vector < vector < int > > VHChannel;
    int nChannels;
    int nChannelsVH;
    vector < vector < int > > VHID;
    vector < int > surfaceID;

    ClassDef(InstalledStation,1);

};



//class Detector : public TObject {
class Detector {
    private:
        static const int freq_step_max = 60;
        static const int ang_step_max = 2664;
        void ReadVgain(string filename);
        void ReadHgain(string filename);
        double Vgain[freq_step_max][ang_step_max];
        double Hgain[freq_step_max][ang_step_max];
        double Freq[freq_step_max];

        void ReadFilter(string filename, Settings *settings1);
        double FilterGain[freq_step_max];   // Filter gain (dB) for Detector freq bin array
        vector <double> FilterGain_databin;   // Filter gain (dB) for DATA_BIN_SIZE bin array

        void ReadPreamp(string filename, Settings *settings1);
        double PreampGain[freq_step_max];   // Filter gain (dB) for Detector freq bin array
        vector <double> PreampGain_databin;   // Filter gain (dB) for DATA_BIN_SIZE bin array


        void ReadFOAM(string filename, Settings *settings1);
        double FOAMGain[freq_step_max];   // Filter gain (dB) for Detector freq bin array
        vector <double> FOAMGain_databin;   // Filter gain (dB) for DATA_BIN_SIZE bin array


        void FlattoEarth_ARA(IceModel *icesurface);
        void FlattoEarth_ARA_sharesurface(IceModel *icesurface);  // each station share the lowest surface

        int freq_step;
        int ang_step;
        double freq_width;
        double freq_init;
        int Detector_mode;

    public:
        Parameters params;
        Detector ();    //default constructor
        Detector (Settings *settings1, IceModel *icesurface);
        //Detector (int mode, IceModel *icesurface);
        vector <ARA_station> stations;
        vector <Antenna_string> strings;

        vector <double> freq_forfft;

        double GetGain(double freq, double theta, double phi, int ant_m, int ant_o);    //read antenna gain at certain angle, certain type, and certain orientation
        double GetGain(double freq, double theta, double phi, int ant_m);   //read antenna gain at certain angle, certain type. (orientation : default)
        double GetFilterGain(int bin) { return FilterGain[bin]; }   // same bin with Vgain, Hgain
        double GetFilterGain_databin(int bin) { return FilterGain_databin[bin]; }   // bin for FFT

        double GetPreampGain(int bin) { return PreampGain[bin]; }   // same bin with Vgain, Hgain
        double GetPreampGain_databin(int bin) { return PreampGain_databin[bin]; }   // bin for FFT

        double GetFOAMGain(int bin) { return FOAMGain[bin]; }   // same bin with Vgain, Hgain
        double GetFOAMGain_databin(int bin) { return FOAMGain_databin[bin]; }   // bin for FFT

        
        double Getfreq_init() {return freq_init;}

        int Get_mode() {return Detector_mode;}

        int GetFreqBin() {return freq_step;}
        double GetFreq(int bin) {return Freq[bin]*1.e6;} //from MHz to Hz

        vector <double> diode_real; // NFOUR/2 array of t domain tunnel diode response. same with icemc -> anita -> diode_real  but only full bandwidth array 4
        vector <double> fdiode_real_databin;    // NFOUR array of f domain tunnel diode response (FFT of diode_real). also same with icemc -> anita -> fdiode_real  but only full bandwidth array 4
        vector <double> fdiode_real;    // NFOUR/2 array of f domain tunnel diode response (FFT of diode_real). also same with icemc -> anita -> fdiode_real  but only full bandwidth array 4
        vector <double> fdiode_real_double;    // NFOUR array of f domain tunnel diode response (FFT of diode_real). also same with icemc -> anita -> fdiode_real  but only full bandwidth array 4
        
        double TIMESTEP;    // will copy TIMESTEP from Settings
        int NFOUR;          // also will copy NFOUR from Settings

        //TF1 fdiode;   // removed. so not exactly same as icemc, but I believe it doesn't matter
        double maxt_diode;
        int maxt_diode_bin; // above maxt_diode in bin
        int idelaybeforepeak;
        int iwindow;
        int ibinshift;

        void getDiodeModel(Settings *settings1);   // similar with icemc -> anita -> getDiodeModel().  set diode_real and fdiode_real values.
    
//    vector < vector < vector < int > > > ChannelfromStringAntenna;
//    void SetChannelStringAntennaMap();
    int GetChannelfromStringAntenna (int stationNum, int stringnum, int antennanum );
    void GetSSAfromChannel ( int stationNum, int channelNum, int * antennaNum, int * stringNum );
    
    void UseAntennaInfo (int stationNum, Settings *settings1);
   
    
    /*
    struct InstalledStation {
        int nSurfaces;
        int nStrings;
        vector < int > surfaceChannels;
        vector < vector < int > > VHChannel;
        int nChannels;
        int nChannelsVH;
        vector < vector < int > > VHID;
        vector < int > surfaceID;

    };
    */
    
    vector < InstalledStation > InstalledStations;
    
    struct IdealStation{
        int nSurfaces;
        int nStrings;
        vector < int > surfaceChannels;
        vector < vector < int > > VHChannel;
        int nChannels;
        int nChannelsVH;
        vector < vector < int > > VHID;
        vector < int > surfaceID;
        vector < int > IDSurface;
        vector < int > IDAntenna;
        vector < int > IDString;
        
    };
    
    vector < IdealStation > IdealStations;
    

    void SetupInstalledStations();
    void PrepareVectorsInstalled();

    void SetupIdealStations();

    int getAntennafromArbAntID( int stationID, int ant_ID);
    int getStringfromArbAntID( int stationID, int ant_ID);
    
    
    
        ~Detector();    //destructor

        ClassDef(Detector,1);
        
    
    
    
};



#endif //DETECTOR_H
