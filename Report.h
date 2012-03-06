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
        vector <double> arrival_time;        //time from posnu to antenna (t:0 at posnu)
        vector <int> reflection;     // non-reflected : 0,  reflected : 1
        vector <Position> Pol_vector;   // polarization vector at the antenna
        //vector <Position> n_H;  // normalized vector for H pol
        //vector <Position> n_V;  // normalized vector for V pol

        // below freq domain simulation output
        vector < vector <double> > vmmhz;  // signal V/m/MHz for each freq bin
        //
        vector < vector <double> > Vfft;  // signal V preparing for FFT


        // below time domain simulation output
        vector < vector <double> > time;   // time of time domain Askaryan radiation
        vector < vector <double> > Ax;     // vector potential x component
        vector < vector <double> > Ay;
        vector < vector <double> > Az;
        vector < vector <double> > V;   // volt signal at the backend of antenna (from fft)
        //
        //
        vector <double> PeakV;  // peak voltage in time domain
        vector <int> Rank;      // rank of peak voltage between antennas (Rank = 0 for 0 signal)

        
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
        Report (Detector *detector, Settings *settings1);

        void Initialize (Detector *detector, Settings *settings1);


        void Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1);

        Vector GetPolarization (Vector &nnu, Vector &launch_vector);

        void GetParameters (Position &src, Position &trg, Vector &nnu, double &viewangle, double receive_angle, Vector &launch_vector, Vector &receive_vector, Vector &n_trg_slappy, Vector &n_trg_pokey );    // get viewangle, launch, receive vectors  (it reads launch angle as a viewangle and returns actual viewangle)

        double GaintoHeight(double gain, double freq, double n_medium);

        void ApplyAntFactors(Settings *settings1, double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &vmmhz);

        void GetAngleAnt(Vector &rec_vector, Position &antenna, double &ant_theta, double &ant_phi);

        void MakeArraysforFFT(Settings *settings1, Detector *detector, vector <double> &vsignal_array, double *vsignal_forfft);

        void MakeNoiseArraysforFFT(Settings *settings1, double vnoise, double *vnoisesignal_forfft);


        void ReadFilter (string filename, int &N, vector <double> &xfreq, vector <double> &ygain);

        double FindPeak (double *waveform, int n);  // same with icemc; trigger->AntTrigger::FindPeak

        void SetRank(Detector *detector); // set rank (rank of strength of signal at each antenna)



        vector <Station_r> stations;
        vector <String_r> strings;

        ClassDef(Report,1);

};

#endif  //REPORT_H
