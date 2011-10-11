#include "Detector.h"
#include "Tools.h"
#include "Event.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "Constants.h"


Detector::Detector(int mode) {

    //set mode ex) mode 0 = testbed,
    // mode 1 = ARA_1
    // mode 2 = ARA_2
    // mode 3 = ARA_37

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

    string testbed_file = "testbed_info.txt";
    string ARA_N_file = "ARA_N_info.txt";
    string ARA37_file = "ARA37_info.txt";

    string line, label;


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
                        strings[string_id].x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read x : "<<(double)strings[string_id].x<<" string_id : "<<string_id<<endl;
                    }
                    else if (label == "y") {
                        strings[string_id].y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read y : "<<(double)strings[string_id].y<<" string_id : "<<string_id<<endl;
                    }
                    else if (label == "z") {
                        strings[string_id].antennas.push_back(temp_antenna);
                        strings[string_id].antennas[antenna_id].z = atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() );
                        strings[string_id].antennas[antenna_id].type = atoi( line.substr( line.find_first_of(",") + 1).c_str() );
                        cout<<"read z : "<<(double)strings[string_id].antennas[antenna_id].z<<" string_id : "<<string_id<<" antenna_id : "<<antenna_id<<" type : "<<(int)strings[string_id].antennas[antenna_id].type<<endl;
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

        double core_x = 0.;  // all units are in meter
        double core_y = 0.;
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
        if ( ARA_N.is_open() ) {
            while (ARA_N.good() ) {
                getline (ARA_N, line);

                if (line[0] != "/"[0]) {
                    label = line.substr(0, line.find_first_of("=") );

                    if (label == "core_x") {
                        core_x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_x"<<endl;
                    }
                    else if (label == "core_y") {
                        core_y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
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
                    stations[station_count].x = core_x + (double)params.station_spacing * cos( (PI/3.) * (double)station_count );
                    stations[station_count].y = core_y + (double)params.station_spacing * sin( (PI/3.) * (double)station_count );
                    station_count++;
            }
            else if (station_count < (int)params.number_of_stations) {
                    stations[station_count].x = core_x;
                    stations[station_count].y = core_y;
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
//            for (int j=0; j<params.number_of_strings_station; j++) {
//            stations[i].string[0].x = stations[i].x - (R_string / 1.414);
            stations[i].strings[0].x = stations[i].x - (R_string * cos(PI/4.));
            stations[i].strings[0].y = stations[i].y + (R_string * sin(PI/4.));

            stations[i].strings[1].x = stations[i].x + (R_string * cos(PI/4.));
            stations[i].strings[1].y = stations[i].y + (R_string * sin(PI/4.));
            
            stations[i].strings[2].x = stations[i].x - (R_string * cos(PI/4.));
            stations[i].strings[2].y = stations[i].y - (R_string * sin(PI/4.));
            
            stations[i].strings[3].x = stations[i].x + (R_string * cos(PI/4.));
            stations[i].strings[3].y = stations[i].y - (R_string * sin(PI/4.));

            //
            // set antenna postions in borehole
            // and set type (h or v pol antenna) and set orientation (facing x or y)
            for (int j=0; j<params.number_of_strings_station; j++) {
                for (int k=0; k<params.number_of_antennas_string; k++) {
                    stations[i].strings[j].antennas[k].z = -z_max + z_btw*k;

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
            stations[i].surfaces[0].x = stations[i].x + (R_surface * cos(PI/3.));
            stations[i].surfaces[0].y = stations[i].y + (R_surface * sin(PI/3.));

            stations[i].surfaces[1].x = stations[i].x + (R_surface * cos(-PI/3.));
            stations[i].surfaces[1].y = stations[i].y + (R_surface * sin(-PI/3.));

            stations[i].surfaces[2].x = stations[i].x + (R_surface * cos(PI));
//            stations[i].surfaces[2].y = stations[i].y + (R_surface * sin(PI));
            stations[i].surfaces[2].y = stations[i].y;

            stations[i].surfaces[3].x = stations[i].x;
            stations[i].surfaces[3].y = stations[i].y;

        }



        // test read V-pol gain file!!
        ReadVgain("ARA_bicone6in_output.txt");
        // test read H-pol gain file!!
        ReadHgain("ARA_dipoletest1_output.txt");


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

        double core_x = 0.;  // all units are in meter
        double core_y = 0.;
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
                        core_x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_x"<<endl;
                    }
                    else if (label == "core_y") {
                        core_y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
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
            double current_y = y_offset * ( (double)params.stations_per_side - 1 - irow) + core_y;
            int stations_this_row = (2 * (int)params.stations_per_side - 1) - abs((int)params.stations_per_side - 1 - irow);

            for (int istation = 0; istation < stations_this_row; istation++) {
                if (station_count < (int)params.number_of_stations) {
                    stations[station_count].y = current_y;
                    stations[station_count].x = (double)params.station_spacing * ((double)istation - ((double)stations_this_row - 1.) / 2.) + core_x;
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
            stations[i].strings[0].x = stations[i].x - (R_string * cos(PI/4.));
            stations[i].strings[0].y = stations[i].y + (R_string * sin(PI/4.));

            stations[i].strings[1].x = stations[i].x + (R_string * cos(PI/4.));
            stations[i].strings[1].y = stations[i].y + (R_string * sin(PI/4.));
            
            stations[i].strings[2].x = stations[i].x - (R_string * cos(PI/4.));
            stations[i].strings[2].y = stations[i].y - (R_string * sin(PI/4.));
            
            stations[i].strings[3].x = stations[i].x + (R_string * cos(PI/4.));
            stations[i].strings[3].y = stations[i].y - (R_string * sin(PI/4.));



            //
            // set antenna postions in borehole
            // and set type (h or v pol antenna) and set orientation (facing x or y)
            for (int j=0; j<params.number_of_strings_station; j++) {
                for (int k=0; k<params.number_of_antennas_string; k++) {
                    stations[i].strings[j].antennas[k].z = -z_max + z_btw*k;

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
            stations[i].surfaces[0].x = stations[i].x + (R_surface * cos(PI/3.));
            stations[i].surfaces[0].y = stations[i].y + (R_surface * sin(PI/3.));

            stations[i].surfaces[1].x = stations[i].x + (R_surface * cos(-PI/3.));
            stations[i].surfaces[1].y = stations[i].y + (R_surface * sin(-PI/3.));

            stations[i].surfaces[2].x = stations[i].x + (R_surface * cos(PI));
//            stations[i].surfaces[2].y = stations[i].y + (R_surface * sin(PI));
            stations[i].surfaces[2].y = stations[i].y;

            stations[i].surfaces[3].x = stations[i].x;
            stations[i].surfaces[3].y = stations[i].y;

        }

        
        // test read V-pol gain file!!
        ReadVgain("ARA_bicone6in_output.txt");
        // test read H-pol gain file!!
        ReadHgain("ARA_dipoletest1_output.txt");


    }


/////////////////////////////////////////////////////////////////////////////////    



//    return 0;

}



inline void Detector::ReadVgain(string filename) {
    ifstream NecOut( filename.c_str() );

    string line;

    if ( NecOut.is_open() ) {
        while (NecOut.good() ) {

            for (int i=0; i<params.freq_step; i++) {
                getline (NecOut, line);
                if ( line.substr(0, line.find_first_of(":")) == "freq ") {
                    Freq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
                    cout<<"freq["<<i<<"] = "<<Freq[i]<<" MHz"<<endl;
                    getline (NecOut, line); //read SWR
                    getline (NecOut, line); //read names

                    for (int j=0; j<params.ang_step; j++) {
                        getline (NecOut, line); //read data line
                        Vgain[i][j] = atof( line.substr( 18 ).c_str() );  // read gain (not dB)
                        cout<<"Gain : "<<Vgain[i][j]<<endl;

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

            for (int i=0; i<params.freq_step; i++) {
                getline (NecOut, line);
                if ( line.substr(0, line.find_first_of(":")) == "freq ") {
                    Freq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
                    cout<<"freq["<<i<<"] = "<<Freq[i]<<" MHz"<<endl;
                    getline (NecOut, line); //read SWR
                    getline (NecOut, line); //read names

                    for (int j=0; j<params.ang_step; j++) {
                        getline (NecOut, line); //read data line
                        Hgain[i][j] = atof( line.substr( 20 ).c_str() );  // read gain (not dB)
                        cout<<"Gain : "<<Hgain[i][j]<<endl;

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
    if ( freq < params.freq_init ) {
        cout<<"Frequency value is smaller than frequency range with Gain."<<endl;
        cout<<"Frequency value "<<freq<<" will be replaced to minimum frequency value "<<params.freq_init<<endl;
        freq = params.freq_init;
    }
    else if ( freq > (params.freq_init + params.freq_width*((double)params.freq_step-1.) ) ) {
        cout<<"Frequency value is bigger than frequency range with Gain."<<endl;
        cout<<"Frequency value "<<freq<<" will be replaced to maximum frequency value "<< params.freq_init + params.freq_width*((double)params.freq_step-1.) - 0.01 <<endl;
        freq = params.freq_init + params.freq_width*((double)params.freq_step-1.) - 0.01;
    }


//    int fx1 = (int)( (freq + (params.freq_width/2.) - params.freq_init)/params.freq_width );
    int fx1 = (int)( (freq - params.freq_init)/params.freq_width );
    int fx2 = fx1 + 1;
    cout<<"fx1 : "<<fx1<<endl;
    cout<<"fx2 : "<<fx2<<endl;

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

    Gout = ((Gout2 - Gout1)/params.freq_width) * ( freq - (params.freq_init + fx1*params.freq_width) ) + Gout1; // get linear interpolation between two nearest freq bin.


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
    if ( freq < params.freq_init ) {
        cout<<"Frequency value is smaller than frequency range with Gain."<<endl;
        cout<<"Frequency value "<<freq<<" will be replaced to minimum frequency value "<<params.freq_init<<endl;
        freq = params.freq_init;
    }
    else if ( freq > (params.freq_init + params.freq_width*((double)params.freq_step - 1.) ) ) {
        cout<<"Frequency value is bigger than frequency range with Gain."<<endl;
        cout<<"Frequency value "<<freq<<" will be replaced to maximum frequency value "<< params.freq_init + params.freq_width*((double)params.freq_step-1.) - 0.01 <<endl;
        freq = params.freq_init + params.freq_width*((double)params.freq_step-1.) - 0.01;
    }


//    int fx1 = (int)( (freq + (params.freq_width/2.) - params.freq_init)/params.freq_width );
    int fx1 = (int)( (freq - params.freq_init)/params.freq_width );
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

    Gout = ((Gout2 - Gout1)/params.freq_width) * ( freq - (params.freq_init + fx1*params.freq_width) ) + Gout1; // get linear interpolation between two nearest freq bin.


    return Gout;


}


double Antenna::GetG(Detector *D, double freq, double theta, double phi) {

    return D->GetGain(freq, theta, phi, type, orient);
}



double Surface_antenna::GetG(Detector *D, double freq, double theta, double phi) {

    return D->GetGain(freq, theta, phi, type, orient);
}


Detector::~Detector() {
    cout<<"Destruct class Detector"<<endl;
}



