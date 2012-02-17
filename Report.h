//--------------------------------------------------
// class Antenna_Response
//-------------------------------------------------- 
//

#ifndef REPORT_H
#define REPORT_H

#include <vector>

class Detector;
class Event;
class RaySolver;
class Signal;
class IceModel;
class Settings;
class Vector;

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
        vector <double> reflect_ang; // surface reflection angle (if 100 : no reflection case)
        vector <double> Dist;        //Distance between posnu and antenna
        vector <int> reflection;     // non-reflected : 0,  reflected : 1
        vector <Position> Pol_vector;   // polarization vector at the antenna
        //vector <Position> n_H;  // normalized vector for H pol
        //vector <Position> n_V;  // normalized vector for V pol

        // below freq domain simulation output
        vector < vector <double> > vmmhz;  // signal V/m/MHz for each freq bin
        //


        // below time domain simulation output
        vector < vector <double> > time;   // time of time domain Askaryan radiation
        vector < vector <double> > Ax;     // vector potential x component
        vector < vector <double> > Ay;
        vector < vector <double> > Az;
        //
        
        void clear ();  // clear all vector format information for next event

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


        void Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1);

        Vector GetPolarization (Vector &nnu, Vector &launch_vector);
        void GetParameters (Position &src, Position &trg, Vector &nnu, double &viewangle, double receive_angle, Vector &launch_vector, Vector &receive_vector, Vector &n_trg_slappy, Vector &n_trg_pokey );    // get viewangle, launch, receive vectors  (it reads launch angle as a viewangle and returns actual viewangle)

        vector <Station_r> stations;
        vector <String_r> strings;

        ClassDef(Report,1);

};

#endif  //REPORT_H
