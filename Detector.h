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

    int number_of_strings;
    int number_of_antennas;
        
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

    ClassDef(ARA_station,1);
    
};





//class Detector : public TObject {
class Detector {
    private:
        static const int freq_step_max = 60;
        static const int ang_step_max = 2664;
        void ReadVgain(string filename);
        void ReadHgain(string filename);
//        double Vgain[freq_step_max][ang_step_max];
        double Hgain[freq_step_max][ang_step_max];
        double Freq[freq_step_max];

        void FlattoEarth_ARA(IceModel *icesurface);

        int freq_step;
        int ang_step;
        double freq_width;
        double freq_init;

    public:
        double Vgain[freq_step_max][ang_step_max];
        Parameters params;
        Detector ();    //default constructor
        Detector (int mode, IceModel *icesurface);
        vector <ARA_station> stations;
        vector <Antenna_string> strings;
        double GetGain(double freq, double theta, double phi, int ant_m, int ant_o);    //read antenna gain at certain angle, certain type, and certain orientation
        double GetGain(double freq, double theta, double phi, int ant_m);   //read antenna gain at certain angle, certain type. (orientation : default)
        
        double Getfreq_init() {return freq_init;}

        ~Detector();    //destructor


        ClassDef(Detector,1);
        

};



#endif //DETECTOR_H
