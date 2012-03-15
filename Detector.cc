#include "Detector.h"
#include "Tools.h"
#include "Event.h"
#include "IceModel.h"
#include "Settings.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "Constants.h"

ClassImp(Detector);
ClassImp(Parameters);
ClassImp(Surface_antenna);
ClassImp(Antenna);
ClassImp(Antenna_string);
ClassImp(ARA_station);


Detector::Detector() {
    //Default constructor
}


Detector::Detector(Settings *settings1, IceModel *icesurface) {
//Detector::Detector(int mode, IceModel *icesurface) {

    // set freq_forfft for later use
    //

    // set freq_forfft array
    // same with icemc anita class initialization function

    double freqstep=1./(double)(settings1->NFOUR/2)/(settings1->TIMESTEP);
    
    //for (int i=0;i<HALFNFOUR/2;i++) {
    for (int i=0;i<settings1->NFOUR/4;i++) {
	//--------------------------------------------------
	// freq_forfft[2*i]=(double)i*freqstep;
	// freq_forfft[2*i+1]=(double)i*freqstep;
	//-------------------------------------------------- 
	freq_forfft.push_back( (double)i*freqstep );    // even numbers
	freq_forfft.push_back( (double)i*freqstep );    // odd numbers
	
    }
    for (int i=settings1->NFOUR/4;i<settings1->NFOUR/2;i++) {
	//--------------------------------------------------
	// freq_forfft[2*i]=(double)i*freqstep;
	// freq_forfft[2*i+1]=(double)i*freqstep;
	//-------------------------------------------------- 
	freq_forfft.push_back( (double)i*freqstep );    // even numbers
	freq_forfft.push_back( (double)i*freqstep );    // odd numbers
	
    }
    // end of settings freq_forfft
    



    //set mode ex) mode 0 = testbed,
    // mode 1 = ARA_1
    // mode 2 = ARA_2
    // mode 3 = ARA_37
    int mode = settings1->DETECTOR;
    Detector_mode = mode;

    int string_id = -1;
//    int string_id = 0;
    int antenna_id = 0;


//    Parameters params;
//    vector <Antenna_string> strings;
    ARA_station temp_station;
    Antenna_string temp;
    Antenna_string temp_string;
    Antenna temp_antenna;
    Surface_antenna temp_surface;

    params.number_of_strings = 0;
    params.number_of_antennas = 0;

    
    //initialize few params values.
    params.freq_step = 60;
    params.ang_step = 2664;
    params.freq_width = 16.667;
    params.freq_init = 83.333;
    //end initialize


    //copy freq_width, freq_init in params to Detector freq_width, freq_init
    freq_step = params.freq_step;
    ang_step = params.ang_step;
    freq_width = params.freq_width;
    freq_init = params.freq_init;
    //end copy


    string testbed_file = "testbed_info.txt";
    string ARA_N_file = "ARA_N_info.txt";
    string ARA37_file = "ARA37_info.txt";

    string line, label;

//    IceModel *icesurface = new IceModel;
    cout<<"Ice surface at 0,0 : "<<icesurface->Geoid(0.)<<endl;


////////////////////////////////////////////////////////////////////////////////////    

    if (mode == 0) {
        cout<<"\n\tDector mode 0 : testbed"<<endl;
        ifstream testbed( testbed_file.c_str() );
        cout<<"We use "<<testbed_file.c_str()<<" as antenna info."<<endl;


        if ( testbed.is_open() ) {
            while (testbed.good() ) {
                getline (testbed, line);

                if (line[0] != "/"[0]) {
                    label = line.substr(0, line.find_first_of("=") );

                    if (label == "number_of_strings") {
                        params.number_of_strings = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                        for (int i=0; i<(int) params.number_of_strings; i++) {
                            strings.push_back(temp);
                        }
                        cout<<"read numner_of_strings"<<endl;
//                        Parameters.number_of_strings = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                    }
                    else if (label == "antenna_string") {
                        string_id++;
                        antenna_id = 0;
                        cout<<"read antenna_string"<<endl;
//                        if (string_id + 1 >= params.number_of_strings) {
//                            cout<<"Error! Too many strings!"<<endl;
//                        }
                    }
                    else if (label == "x") {
                        //strings[string_id].x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        strings[string_id].SetX( atof( line.substr( line.find_first_of("=") + 1).c_str() ) );
                        cout<<"read x : "<<(double)strings[string_id].GetX()<<" string_id : "<<string_id<<endl;
                    }
                    else if (label == "y") {
                        //strings[string_id].y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        strings[string_id].SetY( atof( line.substr( line.find_first_of("=") + 1).c_str() ) );
                        cout<<"read y : "<<(double)strings[string_id].GetY()<<" string_id : "<<string_id<<endl;
                    }
                    else if (label == "z") {
                        strings[string_id].antennas.push_back(temp_antenna);
                        //strings[string_id].antennas[antenna_id].z = atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() );
                        strings[string_id].antennas[antenna_id].SetZ( atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() ) );
                        strings[string_id].antennas[antenna_id].type = atoi( line.substr( line.find_first_of(",") + 1).c_str() );
                        cout<<"read z : "<<(double)strings[string_id].antennas[antenna_id].GetZ()<<" string_id : "<<string_id<<" antenna_id : "<<antenna_id<<" type : "<<(int)strings[string_id].antennas[antenna_id].type<<endl;
                        antenna_id++;
                        params.number_of_antennas++;
//                        Parameters.number_of_antennas++;
                    }


                }
            }
            testbed.close();
        }
        
        else {
            cout<<"Unable to open antenna array file !"<<endl;
//            return 1;
        }


        // testbed version of FlattoEarth_ARA 
        // strings and antennas on the strings use geoid surface!
        double Dist = 0.;   //for sqrt(x^2 + y^2)
        double R1 = icesurface->Surface(0.,0.); // from core of earth to surface at theta, phi = 0.
        double theta_tmp;
        double phi_tmp;
                    
        // set same theta, phi to all antennas in same string
        for (int i=0; i<params.number_of_strings; i++) {

            Dist = sqrt( pow(strings[i].GetX(),2) + pow(strings[i].GetY(),2) );
            theta_tmp = Dist/R1;    // assume R1 is constant (which is not)
            phi_tmp = atan2(strings[i].GetY(),strings[i].GetX());

            if (phi_tmp<0.) phi_tmp += 2.*PI;

            // set theta, phi for strings.
            strings[i].SetThetaPhi(theta_tmp, phi_tmp);
            //set R for strings.
            strings[i].SetR( icesurface->Surface( strings[i].Lon(), strings[i].Lat()) );

            cout<<"R, Theta, Phi : "<<strings[i].R()<<" "<<strings[i].Theta()<<" "<<strings[i].Phi()<<endl;

            // set antennas r, theta, phi
            for (int j=0; j<antenna_id; j++) {
                strings[i].antennas[j].SetRThetaPhi( strings[i].R() + strings[i].antennas[j].GetZ() , strings[i].Theta(), strings[i].Phi() );
            }
        }






    }

/////////////////////////////////////////////////////////////////////////////




    else if (mode == 1) {
        cout<<"\n\tDector mode 1 : Specific number of stations (less than 7 stations) !"<<endl;
        ifstream ARA_N( ARA_N_file.c_str() );
        cout<<"We use "<<ARA_N_file.c_str()<<" as antenna info."<<endl;


        // initialize info
        params.number_of_stations = 1;
        params.number_of_strings_station = 4;   // ARA-1 has 4 strings
        params.number_of_antennas_string = 4; // 4 antennas on each strings
        params.number_of_surfaces_station = 4;

        //double core_x = 0.; 
        //double core_y = 0.;
        params.core_x = 0.; 
        params.core_y = 0.;
        double R_string = 10.;  // all units are in meter
        double R_surface = 60.;
        double z_max = 200.;
        double z_btw = 20.;
        params.stations_per_side = 4;       // total 37 stations
        params.station_spacing = 2000.;     // 2km spacing
        params.antenna_orientation = 0;     // all antenna facing x
        // finish initialization
        //






        // Read new parameters if there are...
        if ( ARA_N.is_open() ) {
            while (ARA_N.good() ) {
                getline (ARA_N, line);

                if (line[0] != "/"[0]) {
                    label = line.substr(0, line.find_first_of("=") );

                    if (label == "core_x") {
                        params.core_x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_x"<<endl;
                    }
                    else if (label == "core_y") {
                        params.core_y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_y"<<endl;
                    }
                    else if (label == "R_string") {
                        R_string = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_string"<<endl;
                    }
                    else if (label == "R_surface") {
                        R_surface = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_surface"<<endl;
                    }
                    else if (label == "z_max") {
                        z_max = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_max"<<endl;
                    }
                    else if (label == "z_btw") {
                        z_btw = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw"<<endl;
                    }
                    else if (label == "number_of_stations") {
                        params.number_of_stations = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read stations_per_side"<<endl;
                    }
                    else if (label == "station_spacing") {
                        params.station_spacing = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read station_spacting"<<endl;
                    }
                    else if (label == "antenna_orientation") {
                        params.antenna_orientation = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read antenna_orientation"<<endl;
                    }
                }
            }
            ARA_N.close();
        }
        // finished reading new parameters





        //
        // caculate number of stations, strings, antennas 
        params.number_of_strings = params.number_of_stations * params.number_of_strings_station;
        params.number_of_antennas = params.number_of_strings * params.number_of_antennas_string;
        //
        //



        //
        // prepare vectors
        for (int i=0; i<params.number_of_stations; i++) {
            stations.push_back(temp_station);

            for (int j=0; j<params.number_of_surfaces_station; j++) {
                stations[i].surfaces.push_back(temp_surface);
            }

            for (int k=0; k<params.number_of_strings_station; k++) {
                stations[i].strings.push_back(temp_string);

                for (int l=0; l<params.number_of_antennas_string; l++) {
                    stations[i].strings[k].antennas.push_back(temp_antenna);
                }

            }


        }
        // end prepare vectors
        //





        //
        // for ARA-37 (or more than 1 station case), need code for setting position for all 37 stations here!
        //
        int station_count = 0;

        for (int istation = 0; istation < (int)params.number_of_stations; istation++) {
            if (station_count < (int)params.number_of_stations - 1) {
                    //stations[station_count].x = core_x + (double)params.station_spacing * cos( (PI/3.) * (double)station_count );
                    //stations[station_count].y = core_y + (double)params.station_spacing * sin( (PI/3.) * (double)station_count );
                    stations[station_count].SetX( params.core_x + (double)params.station_spacing * cos( (PI/3.) * (double)station_count ) );
                    stations[station_count].SetY( params.core_y + (double)params.station_spacing * sin( (PI/3.) * (double)station_count ) );
                    station_count++;
            }
            else if (station_count < (int)params.number_of_stations) {
                    //stations[station_count].x = core_x;
                    //stations[station_count].y = core_y;
                    stations[station_count].SetX( params.core_x );
                    stations[station_count].SetY( params.core_y );
                    station_count++;
            }
            else {
                    cout<<"\n\tError, too many stations !"<<endl;
            }
        }
        // finished setting all stations' position

                
        cout<<"total station_count : "<<station_count<<endl;
        if (station_count != (int)params.number_of_stations) cout<<"\n\tError, station number not match !"<<endl;


        //
        // set antenna values from parameters
        // set station positions
        for (int i=0; i<params.number_of_stations; i++) {

            //
            // set string postions based on station position
            stations[i].strings[0].SetX( stations[i].GetX() - (R_string * cos(PI/4.)) );
            stations[i].strings[0].SetY( stations[i].GetY() + (R_string * sin(PI/4.)) );

            stations[i].strings[1].SetX( stations[i].GetX() + (R_string * cos(PI/4.)) );
            stations[i].strings[1].SetY( stations[i].GetY() + (R_string * sin(PI/4.)) );
            
            stations[i].strings[2].SetX( stations[i].GetX() - (R_string * cos(PI/4.)) );
            stations[i].strings[2].SetY( stations[i].GetY() - (R_string * sin(PI/4.)) );
            
            stations[i].strings[3].SetX( stations[i].GetX() + (R_string * cos(PI/4.)) );
            stations[i].strings[3].SetY( stations[i].GetY() - (R_string * sin(PI/4.)) );

            //
            // set antenna postions in borehole
            // and set type (h or v pol antenna) and set orientation (facing x or y)
            for (int j=0; j<params.number_of_strings_station; j++) {
                for (int k=0; k<params.number_of_antennas_string; k++) {
                    stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw*k );

                    if (k%2 == 0) {
                        stations[i].strings[j].antennas[k].type = 0;   // v-pol
                    }
                    else {
                        stations[i].strings[j].antennas[k].type = 1;   // h-pol
                    }

                    if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                        stations[i].strings[j].antennas[k].orient = 0;
                    }
                    else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                        if ( j==0||j==3 ) {
                            if ( k==0||k==1 ) {
                                stations[i].strings[j].antennas[k].orient = 0;
                            }
                            else {
                                stations[i].strings[j].antennas[k].orient = 1;
                            }
                        }
                        else {
                            if ( k==0||k==1 ) {
                                stations[i].strings[j].antennas[k].orient = 1;
                            }
                            else {
                                stations[i].strings[j].antennas[k].orient = 0;
                            }
                        }

                    }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                }
            }

            //
            // set surface antenna postions
            stations[i].surfaces[0].SetX( stations[i].GetX() + (R_surface * cos(PI/3.)) );
            stations[i].surfaces[0].SetY( stations[i].GetY() + (R_surface * sin(PI/3.)) );

            stations[i].surfaces[1].SetX( stations[i].GetX() + (R_surface * cos(-PI/3.)) );
            stations[i].surfaces[1].SetY( stations[i].GetY() + (R_surface * sin(-PI/3.)) );

            stations[i].surfaces[2].SetX( stations[i].GetX() + (R_surface * cos(PI)) );
            stations[i].surfaces[2].SetY( stations[i].GetY() );

            stations[i].surfaces[3].SetX( stations[i].GetX() );
            stations[i].surfaces[3].SetY( stations[i].GetY() );

        }



        // change coordinate from flat surface to curved Earth surface
        FlattoEarth_ARA(icesurface);





        // test read V-pol gain file!!
        ReadVgain("ARA_bicone6in_output.txt");
        // test read H-pol gain file!!
        ReadHgain("ARA_dipoletest1_output.txt");
        // read filter file!!
        ReadFilter("./data/filter.csv", settings1);


    }







/////////////////////////////////////////////////////////////////////////////////    



    else if (mode == 2) {
        cout<<"\n\tDector mode 2 : Pentagon"<<endl;
        cout<<"\n\tBy default, ARA-37 is set"<<endl;
        ifstream ARA37( ARA37_file.c_str() );
        cout<<"We use "<<ARA37_file.c_str()<<" as antenna info."<<endl;


        //
        // initialize info
        params.number_of_stations = 37;
        params.number_of_strings_station = 4;   // ARA-1 has 4 strings
        params.number_of_antennas_string = 4; // 4 antennas on each strings
        params.number_of_surfaces_station = 4;

        //double core_x = 0.;  // all units are in meter
        //double core_y = 0.;
        params.core_x = 0.;  // all units are in meter
        params.core_y = 0.;
        double R_string = 10.;
        double R_surface = 60.;
        double z_max = 200.;
        double z_btw = 20.;
        params.stations_per_side = 4;       // total 37 stations
        params.station_spacing = 2000.;     // 2km spacing
        params.antenna_orientation = 0;     // all antenna facing x
        // finish initialization
        //






        // Read new parameters if there are...
        if ( ARA37.is_open() ) {
            while (ARA37.good() ) {
                getline (ARA37, line);

                if (line[0] != "/"[0]) {
                    label = line.substr(0, line.find_first_of("=") );

                    if (label == "core_x") {
                        params.core_x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_x"<<endl;
                    }
                    else if (label == "core_y") {
                        params.core_y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_y"<<endl;
                    }
                    else if (label == "R_string") {
                        R_string = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_string"<<endl;
                    }
                    else if (label == "R_surface") {
                        R_surface = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_surface"<<endl;
                    }
                    else if (label == "z_max") {
                        z_max = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_max"<<endl;
                    }
                    else if (label == "z_btw") {
                        z_btw = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw"<<endl;
                    }
                    else if (label == "stations_per_side") {
                        params.stations_per_side = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read stations_per_side"<<endl;
                    }
                    else if (label == "station_spacing") {
                        params.station_spacing = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read station_spacting"<<endl;
                    }
                    else if (label == "antenna_orientation") {
                        params.antenna_orientation = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read antenna_orientation"<<endl;
                    }
                }
            }
            ARA37.close();
        }
        // finished reading new parameters




        //
        // caculate number of stations, strings, antennas 
        params.number_of_stations = 1 + (3 * params.stations_per_side) * (params.stations_per_side - 1);

        params.number_of_strings = params.number_of_stations * params.number_of_strings_station;
        params.number_of_antennas = params.number_of_strings * params.number_of_antennas_string;
        // 



        //
        // prepare vectors
        for (int i=0; i<params.number_of_stations; i++) {
            stations.push_back(temp_station);

            for (int j=0; j<params.number_of_surfaces_station; j++) {
                stations[i].surfaces.push_back(temp_surface);
            }

            for (int k=0; k<params.number_of_strings_station; k++) {
                stations[i].strings.push_back(temp_string);

                for (int l=0; l<params.number_of_antennas_string; l++) {
                    stations[i].strings[k].antennas.push_back(temp_antenna);
                }

            }


        }
        // end perpare vectors
        //








        //
        // for ARA-37 (or more than 1 station case), need code for setting position for all 37 stations here!
        //
        //
        // here, this only works for pentagon shape!
        //
        double y_offset = (double)params.station_spacing * sqrt(3) / 2.;

        int station_count = 0;

        for (int irow = 0; irow < ((int)params.stations_per_side * 2)-1; irow++) {
            double current_y = y_offset * ( (double)params.stations_per_side - 1 - irow) + params.core_y;
            int stations_this_row = (2 * (int)params.stations_per_side - 1) - abs((int)params.stations_per_side - 1 - irow);

            for (int istation = 0; istation < stations_this_row; istation++) {
                if (station_count < (int)params.number_of_stations) {
                    stations[station_count].SetY( current_y );
                    stations[station_count].SetX( (double)params.station_spacing * ((double)istation - ((double)stations_this_row - 1.) / 2.) + params.core_x );
                    station_count++;
                }
                else {
                    cout<<"\n\tError, too many stations !"<<endl;
                }
            }
        }
        // finished setting all stations' position

                
        cout<<"total station_count : "<<station_count<<endl;
        if (station_count != (int)params.number_of_stations) cout<<"\n\tError, station number not match !"<<endl;


        //
        // set antenna values from parameters
        // set station positions
        for (int i=0; i<params.number_of_stations; i++) {

            //
            // set string postions based on station position
//            for (int j=0; j<params.number_of_strings_station; j++) {
//            stations[i].string[0].x = stations[i].x - (R_string / 1.414);
            stations[i].strings[0].SetX( stations[i].GetX() - (R_string * cos(PI/4.)) );
            stations[i].strings[0].SetY( stations[i].GetY() + (R_string * sin(PI/4.)) );

            stations[i].strings[1].SetX( stations[i].GetX() + (R_string * cos(PI/4.)) );
            stations[i].strings[1].SetY( stations[i].GetY() + (R_string * sin(PI/4.)) );
            
            stations[i].strings[2].SetX( stations[i].GetX() - (R_string * cos(PI/4.)) );
            stations[i].strings[2].SetY( stations[i].GetY() - (R_string * sin(PI/4.)) );
            
            stations[i].strings[3].SetX( stations[i].GetX() + (R_string * cos(PI/4.)) );
            stations[i].strings[3].SetY( stations[i].GetY() - (R_string * sin(PI/4.)) );



            //
            // set antenna postions in borehole
            // and set type (h or v pol antenna) and set orientation (facing x or y)
            for (int j=0; j<params.number_of_strings_station; j++) {
                for (int k=0; k<params.number_of_antennas_string; k++) {
                    stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw*k );

                    if (k%2 == 0) {
                        stations[i].strings[j].antennas[k].type = 0;   // v-pol
                    }
                    else {
                        stations[i].strings[j].antennas[k].type = 1;   // h-pol
                    }

                    if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                        stations[i].strings[j].antennas[k].orient = 0;
                    }
                    else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                        if ( j==0||j==3 ) {
                            if ( k==0||k==1 ) {
                                stations[i].strings[j].antennas[k].orient = 0;
                            }
                            else {
                                stations[i].strings[j].antennas[k].orient = 1;
                            }
                        }
                        else {
                            if ( k==0||k==1 ) {
                                stations[i].strings[j].antennas[k].orient = 1;
                            }
                            else {
                                stations[i].strings[j].antennas[k].orient = 0;
                            }
                        }

                    }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                }
            }





            //
            // set surface antenna postions
            stations[i].surfaces[0].SetX( stations[i].GetX() + (R_surface * cos(PI/3.)) );
            stations[i].surfaces[0].SetY( stations[i].GetY() + (R_surface * sin(PI/3.)) );

            stations[i].surfaces[1].SetX( stations[i].GetX() + (R_surface * cos(-PI/3.)) );
            stations[i].surfaces[1].SetY( stations[i].GetY() + (R_surface * sin(-PI/3.)) );

            stations[i].surfaces[2].SetX( stations[i].GetX() + (R_surface * cos(PI)) );
//            stations[i].surfaces[2].y = stations[i].y + (R_surface * sin(PI));
            stations[i].surfaces[2].SetY( stations[i].GetY() );

            stations[i].surfaces[3].SetX( stations[i].GetX() );
            stations[i].surfaces[3].SetY( stations[i].GetY() );

        }



        // change coordinate from flat surface to curved Earth surface
        FlattoEarth_ARA(icesurface);

        
        
        // test read V-pol gain file!!
        ReadVgain("ARA_bicone6in_output.txt");
        // test read H-pol gain file!!
        ReadHgain("ARA_dipoletest1_output.txt");
        // read filter file!!
        ReadFilter("./data/filter.csv", settings1);


    }


/////////////////////////////////////////////////////////////////////////////////    




//    return 0;

    cout<<"test2"<<endl;
}



inline void Detector::ReadVgain(string filename) {
    ifstream NecOut( filename.c_str() );

    string line;

    if ( NecOut.is_open() ) {
        while (NecOut.good() ) {

            for (int i=0; i<freq_step; i++) {
                getline (NecOut, line);
                if ( line.substr(0, line.find_first_of(":")) == "freq ") {
                    Freq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
//                    cout<<"freq["<<i<<"] = "<<Freq[i]<<" MHz"<<endl;
                    getline (NecOut, line); //read SWR
                    getline (NecOut, line); //read names

                    for (int j=0; j<ang_step; j++) {
                        getline (NecOut, line); //read data line
                        Vgain[i][j] = atof( line.substr( 18 ).c_str() );  // read gain (not dB)
//                        cout<<"Gain : "<<Vgain[i][j]<<endl;

                    }// end ang_step

                }// end check freq label

            }// end freq_step

        }// end while NecOut.good
        NecOut.close();
    }// end if file open

}// end ReadVgain



inline void Detector::ReadHgain(string filename) {
    ifstream NecOut( filename.c_str() );

    string line;

    if ( NecOut.is_open() ) {
        while (NecOut.good() ) {

            for (int i=0; i<freq_step; i++) {
                getline (NecOut, line);
                if ( line.substr(0, line.find_first_of(":")) == "freq ") {
                    Freq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
//                    cout<<"freq["<<i<<"] = "<<Freq[i]<<" MHz"<<endl;
                    getline (NecOut, line); //read SWR
                    getline (NecOut, line); //read names

                    for (int j=0; j<ang_step; j++) {
                        getline (NecOut, line); //read data line
                        Hgain[i][j] = atof( line.substr( 20 ).c_str() );  // read gain (not dB)
//                        cout<<"Gain : "<<Hgain[i][j]<<endl;

                    }// end ang_step

                }// end check freq label

            }// end freq_step

        }// end while NecOut.good
        NecOut.close();
    }// end if file open

}// end ReadHgain




double Detector::GetGain(double freq, double theta, double phi, int ant_m, int ant_o) { // using Interpolation on multidimentions!
//double GetGain(double freq, double theta, double phi, int ant_m, int ant_o) { // using Interpolation on multidimentions!

    Parameters params;

    // change antenna facing orientation
    if (ant_o == 0) {
        // no change...
    }
    else if (ant_o == 1) {
        if (phi - 90. >= 0.) {
            phi = phi - 90.;
        }
        else {
            phi = 360. + phi - 90.;
        }
    }
    else if (ant_o == 2) {
        if (phi - 180. >= 0.) {
            phi = phi - 180.;
        }
        else {
            phi = 360. + phi - 180.;
        }
    }
    else if (ant_o == 3) {
        if (phi - 270. >= 0.) {
            phi = phi - 270.;
        }
        else {
            phi = 360. + phi - 270.;
        }
    }
    else {
        cout<<"Wrong option selected for antenna orientation "<<ant_o<<" !!"<<endl;
        cout<<"ant_o will be replaced from "<<ant_o<<" to 0"<<endl;
    }
    // end changing antenna orientation


    int i = (int)(theta/5.);
    int j = (int)(phi/5.);

    double thetai = 5.*( (int)(theta/5.) );
    double thetai1 = 5.*( (int)(theta/5.) + 1.);
    double phij = 5.*( (int)(phi/5.) );
    double phij1 = 5.*( (int)(phi/5.) + 1.);

    double t = (theta - thetai)/(thetai1 - thetai);
    double u = (phi - phij)/(phij1 - phij);

    // in case when freq is out of nec2 freq range. use nearest min/max freq bin value. 
    if ( freq < freq_init ) {
        cout<<"Frequency value is smaller than frequency range with Gain."<<endl;
        cout<<"Frequency value "<<freq<<" will be replaced to minimum frequency value "<<freq_init<<endl;
        freq = freq_init;
    }
    else if ( freq > (freq_init + freq_width*((double)freq_step-1.) ) ) {
        cout<<"Frequency value is bigger than frequency range with Gain."<<endl;
        cout<<"Frequency value "<<freq<<" will be replaced to maximum frequency value "<< freq_init + freq_width*((double)freq_step-1.) - 0.01 <<endl;
        freq = freq_init + freq_width*((double)freq_step-1.) - 0.01;
    }


//    int fx1 = (int)( (freq + (freq_width/2.) - freq_init)/freq_width );
    int fx1 = (int)( (freq - freq_init)/freq_width );
    int fx2 = fx1 + 1;
//    cout<<"fx1 : "<<fx1<<endl;
//    cout<<"fx2 : "<<fx2<<endl;

    double Gij, Gi1j, Gij1, Gi1j1, Gout1, Gout2, Gout;

    if (ant_m == 0) {   // for V pol antenna!!
        Gij = Vgain[fx1][(int)(37*j+i)];
        Gi1j = Vgain[fx1][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx1][(int)(i)];
            Gi1j1 = Vgain[fx1][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx1][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx1][(int)(37*(j+1)+i+1)];
        }

        Gout1 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest smaller freq bin

        Gij = Vgain[fx2][(int)(37*j+i)];
        Gi1j = Vgain[fx2][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx2][(int)(i)];
            Gi1j1 = Vgain[fx2][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx2][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx2][(int)(37*(j+1)+i+1)];
        }

        Gout2 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest higher freq bin
    
    }

    else if (ant_m == 1) {   // for H pol antenna!!
        Gij = Hgain[fx1][(int)(37*j+i)];
        Gi1j = Hgain[fx1][(int)(37*j+i+1)];
        if ( j == 71 ) {
            Gij1 = Hgain[fx1][(int)(i)];
            Gi1j1 = Hgain[fx1][(int)(i+1)];
        }
        else {
            Gij1 = Hgain[fx1][(int)(37*(j+1)+i)];
            Gi1j1 = Hgain[fx1][(int)(37*(j+1)+i+1)];
        }

        Gout1 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest smaller freq bin

        Gij = Vgain[fx2][(int)(37*j+i)];
        Gi1j = Vgain[fx2][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx2][(int)(i)];
            Gi1j1 = Vgain[fx2][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx2][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx2][(int)(37*(j+1)+i+1)];
        }

        Gout2 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest higher freq bin
    }

    else {
        cout<<"There is no antenna type : "<<ant_m<<" !!"<<endl;
        cout<<"Will return Gain = 0 !!"<<endl;
        Gout1 = 0.;
        Gout2 = 0.;
    }

    Gout = ((Gout2 - Gout1)/freq_width) * ( freq - (freq_init + fx1*freq_width) ) + Gout1; // get linear interpolation between two nearest freq bin.


    return Gout;

// ant_o face x = 0, y = 1, -x = 2, -y = 3

}


double Detector::GetGain(double freq, double theta, double phi, int ant_m) {
//double GetGain(double freq, double theta, double phi, int ant_m) {

    Parameters params;

    int i = (int)(theta/5.);
    int j = (int)(phi/5.);

    double thetai = 5.*( (int)(theta/5.) );
    double thetai1 = 5.*( (int)(theta/5.) + 1.);
    double phij = 5.*( (int)(phi/5.) );
    double phij1 = 5.*( (int)(phi/5.) + 1.);

    double t = (theta - thetai)/(thetai1 - thetai);
    double u = (phi - phij)/(phij1 - phij);


    // in case when freq is out of nec2 freq range. use nearest min/max freq bin value. 
    if ( freq < freq_init ) {
        cout<<"Frequency value is smaller than frequency range with Gain."<<endl;
        cout<<"Frequency value "<<freq<<" will be replaced to minimum frequency value "<<freq_init<<endl;
        freq = freq_init;
    }
    else if ( freq > (freq_init + freq_width*((double)freq_step - 1.) ) ) {
        cout<<"Frequency value is bigger than frequency range with Gain."<<endl;
        cout<<"Frequency value "<<freq<<" will be replaced to maximum frequency value "<< freq_init + freq_width*((double)freq_step-1.) - 0.01 <<endl;
        freq = freq_init + freq_width*((double)freq_step-1.) - 0.01;
    }


//    int fx1 = (int)( (freq + (freq_width/2.) - freq_init)/freq_width );
    int fx1 = (int)( (freq - freq_init)/freq_width );
    int fx2 = fx1 + 1;

    double Gij, Gi1j, Gij1, Gi1j1, Gout1, Gout2, Gout;

    if (ant_m == 0) {   // for V pol antenna!!
        Gij = Vgain[fx1][(int)(37*j+i)];
        Gi1j = Vgain[fx1][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx1][(int)(i)];
            Gi1j1 = Vgain[fx1][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx1][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx1][(int)(37*(j+1)+i+1)];
        }

        Gout1 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest smaller freq bin

        Gij = Vgain[fx2][(int)(37*j+i)];
        Gi1j = Vgain[fx2][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx2][(int)(i)];
            Gi1j1 = Vgain[fx2][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx2][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx2][(int)(37*(j+1)+i+1)];
        }

        Gout2 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest higher freq bin
    
    }

    else if (ant_m == 1) {   // for H pol antenna!!
        Gij = Hgain[fx1][(int)(37*j+i)];
        Gi1j = Hgain[fx1][(int)(37*j+i+1)];
        if ( j == 71 ) {
            Gij1 = Hgain[fx1][(int)(i)];
            Gi1j1 = Hgain[fx1][(int)(i+1)];
        }
        else {
            Gij1 = Hgain[fx1][(int)(37*(j+1)+i)];
            Gi1j1 = Hgain[fx1][(int)(37*(j+1)+i+1)];
        }

        Gout1 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest smaller freq bin

        Gij = Vgain[fx2][(int)(37*j+i)];
        Gi1j = Vgain[fx2][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx2][(int)(i)];
            Gi1j1 = Vgain[fx2][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx2][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx2][(int)(37*(j+1)+i+1)];
        }

        Gout2 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest higher freq bin
    }

    else {
        cout<<"There is no antenna type : "<<ant_m<<" !!"<<endl;
        cout<<"Will return Gain = 0 !!"<<endl;
        Gout1 = 0.;
        Gout2 = 0.;
    }

    Gout = ((Gout2 - Gout1)/freq_width) * ( freq - (freq_init + fx1*freq_width) ) + Gout1; // get linear interpolation between two nearest freq bin.


    return Gout;


}


double Antenna::GetG(Detector *D, double freq, double theta, double phi) {

    return D->GetGain(freq, theta, phi, type, orient);
}



double Surface_antenna::GetG(Detector *D, double freq, double theta, double phi) {

    return D->GetGain(freq, theta, phi, type, orient);
}


inline void Detector::FlattoEarth_ARA(IceModel *icesurface) {
    
    double Dist = 0.;   //for sqrt(x^2 + y^2)
    double R1 = icesurface->Surface(0.,0.); // from core of earth to surface at theta, phi = 0.
//--------------------------------------------------
//     double R1 = icesurface->Geoid(0.); // from core of earth to surface at theta, phi = 0.
//-------------------------------------------------- 
    double theta_tmp;
    double phi_tmp;

        // stations
        // stations, strings, and borehole antennas use geoid surface !!
        for (int i=0; i<params.number_of_stations; i++) {

            Dist = sqrt( pow(stations[i].GetX(),2) + pow(stations[i].GetY(),2) );
            theta_tmp = Dist/R1;
            phi_tmp = atan2(stations[i].GetY(),stations[i].GetX());

            if (phi_tmp<0.) phi_tmp += 2.*PI;

            // set theta, phi for stations.
            stations[i].SetThetaPhi(theta_tmp, phi_tmp);
            //set R for stations.
            stations[i].SetR( icesurface->Surface( stations[i].Lon(), stations[i].Lat()) );


            // strings
            for (int j=0; j<params.number_of_strings_station; j++) {
                Dist = sqrt( pow(stations[i].strings[j].GetX(),2) + pow(stations[i].strings[j].GetY(),2) );
                theta_tmp = Dist/R1;
                phi_tmp = atan2(stations[i].strings[j].GetY(),stations[i].strings[j].GetX());

                if (phi_tmp<0.) phi_tmp += 2.*PI;

                stations[i].strings[j].SetThetaPhi(theta_tmp, phi_tmp);
                // string Vector points the position where string meets the ice surface!
                stations[i].strings[j].SetR( icesurface->Surface( stations[i].strings[j].Lon(), stations[i].strings[j].Lat()) );
                
                

                // borehole antennas
                for (int k=0; k<params.number_of_antennas_string; k++) {
                    stations[i].strings[j].antennas[k].SetRThetaPhi( stations[i].strings[j].R() + stations[i].strings[j].antennas[k].GetZ() , stations[i].strings[j].Theta(), stations[i].strings[j].Phi() );
                }


            }
            
            // surface antennas
            // surface antennas are on actual ice surface (not geoid surface)
            for (int l=0; l<params.number_of_surfaces_station; l++) {
                Dist = sqrt( pow(stations[i].surfaces[l].GetX(),2) + pow(stations[i].surfaces[l].GetY(),2) );
                theta_tmp = Dist/R1;
                phi_tmp = atan2(stations[i].surfaces[l].GetY(),stations[i].surfaces[l].GetX());

                if (phi_tmp<0.) phi_tmp += 2.*PI;

                stations[i].surfaces[l].SetThetaPhi(theta_tmp, phi_tmp);
                stations[i].surfaces[l].SetR( icesurface->Surface( stations[i].surfaces[l].Lon(), stations[i].surfaces[l].Lat()) );
            }


        }


}



inline void Detector::ReadFilter(string filename, Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain

    ifstream Filter( filename.c_str() );

    string line;

    int N=-1;

    vector <double> xfreq_tmp;
    vector <double> ygain_tmp;

    if ( Filter.is_open() ) {
        while (Filter.good() ) {

            getline (Filter, line);
            //xfreq[N] = atof( line.substr(0, line.find_first_of(",")).c_str() );
            xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) );
            //xfreq.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) * 1.E6 );  // from MHz to Hz

            //xfreq[N] = xfreq[N] * 1.E6; // from MHz to Hz

            //ygain[N] = atof( line.substr(line.find_first_of(",") + 1).c_str() );
            ygain_tmp.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) );
            //ygain.push_back( pow(pow(10,atof( line.substr(line.find_first_of(",") + 1).c_str() ) /10.),0.5) );  // from dB to unitless gain for voltage

            //ygain[N] = pow(pow(10,yFilter[i]/10.0),0.5);    // from dB to field strength unitless gain
            
            N++;

        }
        Filter.close();
    }

    else cout<<"Filter file can not opened!!"<<endl;

    double xfreq[N], ygain[N];  // need array for Tools::SimpleLinearInterpolation
    double xfreq_fft[settings1->NFOUR/4];   // array for FFT freq bin
    double ygain_fft[settings1->NFOUR/4];   // array for gain in FFT bin
    double df_fft;

    df_fft = 1./ ( (double)(settings1->NFOUR/2) * settings1->TIMESTEP );

    for (int i=0;i<N;i++) { // copy values
        xfreq[i] = xfreq_tmp[i];
        ygain[i] = ygain_tmp[i];
    }
    for (int i=0;i<settings1->NFOUR/4;i++) {    // this one is for FFT freq bin
        xfreq_fft[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }


    // Tools::SimpleLinearInterpolation will return Filter array (in dB)
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, freq_step, Freq, FilterGain );

    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->NFOUR/4, xfreq_fft, ygain_fft );

    for (int i=0;i<settings1->NFOUR/4;i++) {
        FilterGain_fft.push_back( ygain_fft[i] );
    }




}





Detector::~Detector() {
    cout<<"Destruct class Detector"<<endl;
}



