//--------------------------------------------------
// class Antenna_Response
//-------------------------------------------------- 
//

#ifndef REPORT_H
#define REPORT_H

#include <vector>
#include "Position.h"
#include "Vector.h"

class Position;
class Detector;

using namespace std;

class Surface_antenna_r {
    public:

        ClassDef(Surface_antenna_r,1);
};

class Antenna_r {
    public:
        // one dimention for number of solutions, (if there,) another dimention for array of information
        //

        int ray_sol_cnt;    // number of RaySolver solutions

        //vector <int> trg;    // if antenna recieved any signal or not. 0 : no signal,  1 : yes signal

        vector <double> view_ang;    //viewing angle
        vector <double> rec_ang;     //receiving angle
        vector <double> reflect_ang; // surface reflection angle (if 100 : non reflected case)
        vector <double> Dist;        //Distance between posnu and antenna
        vector <int> reflection;     // non-reflected : 0,  reflected : 1

        // below freq domain simulation output
        vector < vector <double> > vmmhz;  // signal V/m/MHz for each freq bin
        //


        // below time domain simulation output
        vector < vector <double> > time;   // time of time domain Askaryan radiation
        vector < vector <double> > Ax;     // vector potential x component
        vector < vector <double> > Ay;
        vector < vector <double> > Az;
        //

        ClassDef(Antenna_r,1);
};

class String_r {
    public:
        //int trg;    // if any antenna trigg in the event. 0 : no antenna in the string trg
                    //                                    1 : 1 or more antenna trg

        vector <Antenna_r> antennas;

        ClassDef(String_r,1);
};

class Station_r {
    public:
        //int trg;    // if any antenna trigg in the event. 0 : no antenna trg
                    //                                    1: 1 or more antenna trg 
        vector <String_r> strings;
        vector <Surface_antenna_r> surfaces;

        ClassDef(Station_r,1);
};



class Report {
    private:

    public:
        //int trg;    // if any antenna in entire detectors trg. 0 : no antenna trg
                    //                                         1 : 1 or more antenna trg

        Report ();
        Report (Detector *detector);

        void Initialize (Detector *detector);

        vector <Station_r> stations;
        vector <String_r> strings;

        ClassDef(Report,1);

};

#endif  //REPORT_H
