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

struct Surf_antenna {
    float x, y;
    int model;
};

struct Antenna {
    float z;
    int model;  // model 0 : v-pol (bicone), model 1 : h-pol (bowtie)
    //add more later
};

struct Antenna_string {
    float x, y;
//    int number_of_antennas;
    vector <Antenna> antennas;
};

struct ARA_station {
    float x, y;
    vector <Antenna_string> ARA_strings;
    vector <Surf_antenna> surfs;
};

/*
struct Parameters {
    int number_of_strings;
    int number_of_antennas;
} params;
*/

struct Parameters {

    int stations_per_side;           //number of stations on one side of pentagon.
    float station_spacing;          // space between stations.

    int number_of_stations;         //total stations
    int number_of_strings_station;  // strings per station
    int number_of_antennas_string;  //antennas per string
    int number_of_surfs_station;    //surface antennas per station

    int number_of_strings;
    int number_of_antennas;
};
/*
struct Event {
    // add more later
};
*/
//        vector <Antenna_string> strings;

class Detector {
    private:


    public:
        Detector (int mode);
        vector <ARA_station> stations;
        vector <Antenna_string> strings;
        Parameters params;
        ~Detector();
        

};



#endif //DETECTOR_H
