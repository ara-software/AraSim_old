////////////////////////////////////////////////////////////////////////////////////////////////
//class Detector:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef DETECTOR_H
#define DETECTOR_H


#include <TObject.h>
#include <vector>
#include "Trigger.h"
#include "Position.h"

using namespace std;

//#include "Trigger.h"
//class Event;
//class Efficiencies;
//class Position;

class Detector;

struct Parameters {

    int stations_per_side;           //number of stations on one side of pentagon.
    double station_spacing;          // space between stations.

    int antenna_orientation;        //antenna orientation setting.

    int number_of_stations;         //total stations
    int number_of_strings_station;  // strings per station
    int number_of_antennas_string;  //antennas per string
    int number_of_surfaces_station;    //surface antennas per station

    int number_of_strings;
    int number_of_antennas;
        
    static const int freq_step = 60;
    static const int ang_step = 2664;
    static const double freq_width = 16.667;  // this value could be changed when Nec2 condition changes!
    static const double freq_init = 83.333;  // this value could be changed when Nec2 condition changes!


};





struct Surface_antenna {

    double x, y;
    int type;   // need to be defined
    int orient; //0 : facing x, 1 : y, 2 : -x, 3 : -y.

    double GetG(Detector *D, double freq, double theta, double phi);    // read gain value from Detector class

};


    
struct Antenna {

    double z;
    int type;  // type 0 : v-pol (bicone), type 1 : h-pol (bowtie for testbed, QSC for ARA)
    int orient; // 0 : facing x, 1 : y, 2 : -x, 3 : -y.

    double GetG(Detector *D, double freq, double theta, double phi);    // read gain value from Detector class
    
};



struct Antenna_string {
    double x, y;
//    int number_of_antennas;
    vector <Antenna> antennas;
};

struct ARA_station {
    double x, y;
    vector <Antenna_string> strings;
    vector <Surface_antenna> surfaces;

    
};





class Detector {
    private:
        static const int freq_step_max = 60;
        static const int ang_step_max = 2664;
        void ReadVgain(string filename);
        void ReadHgain(string filename);
        double Vgain[freq_step_max][ang_step_max];
        double Hgain[freq_step_max][ang_step_max];
        double Freq[freq_step_max];

    public:
        Parameters params;
        Detector (int mode);
        vector <ARA_station> stations;
        vector <Antenna_string> strings;
        double GetGain(double freq, double theta, double phi, int ant_m, int ant_o);    //read antenna gain at certain angle, certain type, and certain orientation
        double GetGain(double freq, double theta, double phi, int ant_m);   //read antenna gain at certain angle, certain type. (orientation : default)

        ~Detector();
        

};



#endif //DETECTOR_H
