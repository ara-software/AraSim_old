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
        vector < vector <double> > VHz_antfactor;  // after applying ApplyAntFactors to vmmhz above ( 1/sqrt2 * 1/dt * 0.5 * heff * pol_factor )
        vector < vector <double> > VHz_filter;  // after applying ApplyAntFactors above and then apply filter gain from detector->GetFilterGain
        vector < vector <double> > Vfft;  // signal V preparing for FFT
        vector < vector <double> > Vfft_noise;  // noise V preparing for FFT


        // below time domain simulation output
        vector < vector <double> > time;   // time of time domain Askaryan radiation
        vector < vector <double> > Ax;     // vector potential x component
        vector < vector <double> > Ay;
        vector < vector <double> > Az;
        vector < vector <double> > V;   // volt signal with all factors applied (as far as we can) (from fft)

        vector < vector <double> > V_noise; // volt noise signal (with all factors applied as far as we can) (from thermal noise + fft)
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
        vector <double> noise_phase;    // random noise phase generated in GetNoisePhase()

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

        void ApplyAntFactors(double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &vmmhz);

        void ApplyFilter(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFilter_fft(int bin_n, Detector *detector, double &vmmhz);

        void GetAngleAnt(Vector &rec_vector, Position &antenna, double &ant_theta, double &ant_phi);

        void GetNoiseWaveforms(Settings *settings1, Detector *detector, double vhz_noise, double *vnoise);
        void GetNoisePhase(Settings *settings1);

        void MakeArraysforFFT(Settings *settings1, Detector *detector, vector <double> &vsignal_array, double *vsignal_forfft);
        void MakeArraysforFFT_noise(Settings *settings1, Detector *detector, vector <double> &vsignal_array, double *vsignal_forfft);


        double FindPeak (double *waveform, int n);  // same with icemc; trigger->AntTrigger::FindPeak

        void SetRank(Detector *detector); // set rank (rank of strength of signal at each antenna)

        vector <double> Vfft_noise_after;   // noise Vfft after get_random_rician
        vector <double> Vfft_noise_before;   // noise Vfft before get_random_rician
        double Vfft_noise_org;              // V/Hz for thermal noise from Johnson-Nyquist


        vector <Station_r> stations;
        vector <String_r> strings;

        ClassDef(Report,1);

};

#endif  //REPORT_H
