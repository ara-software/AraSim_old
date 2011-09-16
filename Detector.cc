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
    Surf_antenna temp_surf;

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
                        cout<<"read x : "<<(float)strings[string_id].x<<" string_id : "<<string_id<<endl;
                    }
                    else if (label == "y") {
                        strings[string_id].y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read y : "<<(float)strings[string_id].y<<" string_id : "<<string_id<<endl;
                    }
                    else if (label == "z") {
                        strings[string_id].antennas.push_back(temp_antenna);
                        strings[string_id].antennas[antenna_id].z = atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() );
                        strings[string_id].antennas[antenna_id].model = atoi( line.substr( line.find_first_of(",") + 1).c_str() );
                        cout<<"read z : "<<(float)strings[string_id].antennas[antenna_id].z<<" string_id : "<<string_id<<" antenna_id : "<<antenna_id<<" model : "<<(int)strings[string_id].antennas[antenna_id].model<<endl;
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
        params.number_of_surfs_station = 4;

        float core_x = 0.;  // all units are in meter
        float core_y = 0.;
        float R_string = 10.;
        float R_surf = 60.;
        float z_max = 200.;
        float z_btw = 20.;
        params.stations_per_side = 4;       // total 37 stations
        params.station_spacing = 2000.;     // 2km spacing
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
                    else if (label == "R_surf") {
                        R_surf = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_surf"<<endl;
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

            for (int j=0; j<params.number_of_surfs_station; j++) {
                stations[i].surfs.push_back(temp_surf);
            }

            for (int k=0; k<params.number_of_strings_station; k++) {
                stations[i].ARA_strings.push_back(temp_string);

                for (int l=0; l<params.number_of_antennas_string; l++) {
                    stations[i].ARA_strings[k].antennas.push_back(temp_antenna);
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
                    stations[station_count].x = core_x + (float)params.station_spacing * cos( (PI/3.) * (float)station_count );
                    stations[station_count].y = core_y + (float)params.station_spacing * sin( (PI/3.) * (float)station_count );
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
            stations[i].ARA_strings[0].x = stations[i].x - (R_string * cos(PI/4.));
            stations[i].ARA_strings[0].y = stations[i].y + (R_string * sin(PI/4.));

            stations[i].ARA_strings[1].x = stations[i].x + (R_string * cos(PI/4.));
            stations[i].ARA_strings[1].y = stations[i].y + (R_string * sin(PI/4.));
            
            stations[i].ARA_strings[2].x = stations[i].x - (R_string * cos(PI/4.));
            stations[i].ARA_strings[2].y = stations[i].y - (R_string * sin(PI/4.));
            
            stations[i].ARA_strings[3].x = stations[i].x + (R_string * cos(PI/4.));
            stations[i].ARA_strings[3].y = stations[i].y - (R_string * sin(PI/4.));

            //
            // set antenna postions in borehole
            for (int j=0; j<params.number_of_strings_station; j++) {
                for (int k=0; k<params.number_of_antennas_string; k++) {
                    stations[i].ARA_strings[j].antennas[k].z = -z_max + z_btw*k;
                }
            }

            //
            // set surface antenna postions
            stations[i].surfs[0].x = stations[i].x + (R_surf * cos(PI/3.));
            stations[i].surfs[0].y = stations[i].y + (R_surf * sin(PI/3.));

            stations[i].surfs[1].x = stations[i].x + (R_surf * cos(-PI/3.));
            stations[i].surfs[1].y = stations[i].y + (R_surf * sin(-PI/3.));

            stations[i].surfs[2].x = stations[i].x + (R_surf * cos(PI));
//            stations[i].surfs[2].y = stations[i].y + (R_surf * sin(PI));
            stations[i].surfs[2].y = stations[i].y;

            stations[i].surfs[3].x = stations[i].x;
            stations[i].surfs[3].y = stations[i].y;

        }



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
        params.number_of_surfs_station = 4;

        float core_x = 0.;  // all units are in meter
        float core_y = 0.;
        float R_string = 10.;
        float R_surf = 60.;
        float z_max = 200.;
        float z_btw = 20.;
        params.stations_per_side = 4;       // total 37 stations
        params.station_spacing = 2000.;     // 2km spacing
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
                    else if (label == "R_surf") {
                        R_surf = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_surf"<<endl;
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

            for (int j=0; j<params.number_of_surfs_station; j++) {
                stations[i].surfs.push_back(temp_surf);
            }

            for (int k=0; k<params.number_of_strings_station; k++) {
                stations[i].ARA_strings.push_back(temp_string);

                for (int l=0; l<params.number_of_antennas_string; l++) {
                    stations[i].ARA_strings[k].antennas.push_back(temp_antenna);
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
                    stations[station_count].x = (float)params.station_spacing * ((float)istation - ((float)stations_this_row - 1.) / 2.) + core_x;
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
            stations[i].ARA_strings[0].x = stations[i].x - (R_string * cos(PI/4.));
            stations[i].ARA_strings[0].y = stations[i].y + (R_string * sin(PI/4.));

            stations[i].ARA_strings[1].x = stations[i].x + (R_string * cos(PI/4.));
            stations[i].ARA_strings[1].y = stations[i].y + (R_string * sin(PI/4.));
            
            stations[i].ARA_strings[2].x = stations[i].x - (R_string * cos(PI/4.));
            stations[i].ARA_strings[2].y = stations[i].y - (R_string * sin(PI/4.));
            
            stations[i].ARA_strings[3].x = stations[i].x + (R_string * cos(PI/4.));
            stations[i].ARA_strings[3].y = stations[i].y - (R_string * sin(PI/4.));

            //
            // set antenna postions in borehole
            for (int j=0; j<params.number_of_strings_station; j++) {
                for (int k=0; k<params.number_of_antennas_string; k++) {
                    stations[i].ARA_strings[j].antennas[k].z = -z_max + z_btw*k;
                }
            }

            //
            // set surface antenna postions
            stations[i].surfs[0].x = stations[i].x + (R_surf * cos(PI/3.));
            stations[i].surfs[0].y = stations[i].y + (R_surf * sin(PI/3.));

            stations[i].surfs[1].x = stations[i].x + (R_surf * cos(-PI/3.));
            stations[i].surfs[1].y = stations[i].y + (R_surf * sin(-PI/3.));

            stations[i].surfs[2].x = stations[i].x + (R_surf * cos(PI));
//            stations[i].surfs[2].y = stations[i].y + (R_surf * sin(PI));
            stations[i].surfs[2].y = stations[i].y;

            stations[i].surfs[3].x = stations[i].x;
            stations[i].surfs[3].y = stations[i].y;

        }



    }


/////////////////////////////////////////////////////////////////////////////////    


//    return 0;

}





Detector::~Detector() {
    cout<<"Destruct class Detector"<<endl;
}



