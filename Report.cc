#include "Detector.h"
#include "Report.h"
#include "Event.h"
#include "RaySolver.h"
#include "signal.hh"
#include "IceModel.h"
#include "Settings.h"
#include "Vector.h"

#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "Constants.h"

ClassImp(Report);
ClassImp(Antenna_r);
ClassImp(Surface_antenna_r);
ClassImp(String_r);
ClassImp(Station_r);

Report::Report() {
}

Report::Report(Detector *detector) {
    // Default constructor

    Initialize(detector);
    
}


void Report::Initialize(Detector *detector) {
    
    // clear information stored in previous event
    //
    stations.clear();


    // tmp for push_back vector structure
    Antenna_r tmp_antenna;
    String_r tmp_string;
    Station_r tmp_station;
    Surface_antenna_r tmp_surface;

    for (int i=0; i<detector->params.number_of_stations; i++) {
        // vector stations
        stations.push_back(tmp_station);
        for (int j=0; j<detector->params.number_of_strings_station; j++) {
            // vector strings
            stations[i].strings.push_back(tmp_string);
            for (int k=0; k<detector->params.number_of_antennas_string; k++) {
                // vector antennas
                stations[i].strings[j].antennas.push_back(tmp_antenna);
            }
        }
        for (int j=0; j<detector->params.number_of_surfaces_station; j++) {
            // vector surface antennas
            stations[i].surfaces.push_back(tmp_surface);
        }
    }


}

void Antenna_r::clear() {   // if any vector variable added in Antenna_r, need to be added here!
    view_ang.clear();
    rec_ang.clear();
    reflect_ang.clear();
    Dist.clear();
    reflection.clear();
    Pol_vector.clear();
    vmmhz.clear();

    time.clear();
    Ax.clear();
    Ay.clear();
    Az.clear();

}


void Report::Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1) {

    int ray_sol_cnt;
    double viewangle;
    Position launch_vector; // direction of ray at the source
    Position receive_vector;    // direction of ray at the target antenna
    vector < vector <double> > ray_output;

    double vmmhz1m_tmp, vmmhz1m_sum, vmmhz1m_em;    // currently not using vmmhz1m_em
    Position Pol_vector;                            // polarization vector at the source
    double mag;                                     // magnification factor. it can vary in case of plane / spherical wave
    double fresnel;                                 // fresnel factor
    double tmp; // for non use return values


           for (int i = 0; i< detector->params.number_of_stations; i++) {

               for (int j=0; j< detector->params.number_of_strings_station; j++) {

                   for (int k=0; k< detector->params.number_of_antennas_string; k++) {


                       stations[i].strings[j].antennas[k].clear();  // clear data in antenna which stored in previous event

                       // run ray solver, see if solution exist
                       // if not, skip (set something like Sol_No = 0;
                       // if solution exist, calculate view angle and calculate TaperVmMHz

                       if (event->Nu_Interaction[0].pickposnu) {    // if posnu is selected inside the antarctic ic:"<<viewangle<<" th_em:"<<d_theta_em[l]<<" th_had:"<<d_theta_had[l]<<" emfrac:"<<emfrac<<" hadfrac:"<<hadfrac<<" vmmhz1m:"<<vmmhz1m[l]<<endl;e
                           
                           raysolver->Solve_Ray(event->Nu_Interaction[0].posnu, detector->stations[i].strings[j].antennas[k], icemodel, ray_output);   // solve ray between source and antenna
                           
                           ray_sol_cnt = 0;
                           
                           if (raysolver->solution_toggle) {  // if there are solution from raysolver
                              cout<<"ray_output size : "<<ray_output[0].size()<<endl;
                               
                               while ( ray_sol_cnt < ray_output[0].size() ) {   // for number of soultions (could be 1 or 2)
                                  cout<<"Path length : "<<ray_output[0][ray_sol_cnt]<<"\tView angle : "<<ray_output[1][ray_sol_cnt]<<"\tReceipt angle : "<<ray_output[2][ray_sol_cnt]<<endl;

                                  // set viewangle, launch_vector, receive vectors
                                  viewangle = ray_output[1][ray_sol_cnt];
                                  GetParameters ( event->Nu_Interaction[0].posnu,   // posnu
                                                detector->stations[i].strings[j].antennas[k],   // trg antenna
                                                event->nnu,                         // nnu
                                                viewangle,         // inputs launch_angle, returns viewangle
                                                ray_output[2][ray_sol_cnt],         // receive_angle
                                                launch_vector, receive_vector );


                                   // store information to report
                                   stations[i].strings[j].antennas[k].view_ang.push_back(viewangle);
                                   stations[i].strings[j].antennas[k].rec_ang.push_back(ray_output[2][ray_sol_cnt]);
                                   stations[i].strings[j].antennas[k].Dist.push_back(ray_output[0][ray_sol_cnt]);
                                   stations[i].strings[j].antennas[k].reflect_ang.push_back(ray_output[3][ray_sol_cnt]);
                                   
                                   stations[i].strings[j].antennas[k].vmmhz.resize(ray_sol_cnt+1);

                                   // calculate the polarization vector at the source
                                   Pol_vector = GetPolarization (event->nnu, launch_vector);

                                   icemodel->GetFresnel( ray_output[1][ray_sol_cnt],    // launch_angle
                                                        ray_output[2][ray_sol_cnt],     // rec_angle
                                                        ray_output[3][ray_sol_cnt],     // reflect_angle
                                                        event->Nu_Interaction[0].posnu,
                                                        launch_vector,
                                                        receive_vector,
                                                        settings1,
                                                        fresnel,
                                                        mag,
                                                        Pol_vector);                    // input src Pol and return Pol at trg

                                   

                                   if ( ray_output[3][ray_sol_cnt] < PI/2. ) {  // when not reflected at the surface, angle = 100
                                       stations[i].strings[j].antennas[k].reflection.push_back(1);  // say this is reflected ray
                                   }
                                   else {
                                       stations[i].strings[j].antennas[k].reflection.push_back(0);  // say this is not reflected ray
                                   }


                                   stations[i].strings[j].antennas[k].Pol_vector.push_back(Pol_vector); // this Pol_vector is for the target antenna



                                   vmmhz1m_sum = 0;

                                   for (int l=0; l<detector->GetFreqBin(); l++) {   // for detector freq bin numbers


                                      cout<<"TaperVmMHz inputs VA:"<<viewangle<<" th_em:"<<event->Nu_Interaction[0].d_theta_em[l]<<" th_had:"<<event->Nu_Interaction[0].d_theta_had[l]<<" emfrac:"<<event->Nu_Interaction[0].emfrac<<" hadfrac:"<<event->Nu_Interaction[0].hadfrac<<" vmmhz1m:"<<event->Nu_Interaction[0].vmmhz1m[l]<<endl;

                                       vmmhz1m_tmp = event->Nu_Interaction[0].vmmhz1m[l];

                                       signal->TaperVmMHz( viewangle, event->Nu_Interaction[0].d_theta_em[l], event->Nu_Interaction[0].d_theta_had[l], event->Nu_Interaction[0].emfrac, event->Nu_Interaction[0].hadfrac, vmmhz1m_tmp, vmmhz1m_em);
                                      cout<<"TaperVmMHz (1m at view angle) at "<<l<<"th bin : "<<vmmhz1m_tmp<<endl;


                                      // multiply all factors for traveling ice
                                      vmmhz1m_tmp = vmmhz1m_tmp / ray_output[0][ray_sol_cnt] * exp(-0.5*ray_output[0][ray_sol_cnt]/icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0)) * mag * fresnel;  // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
                                      cout<<"AttenLength : "<<icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0)<<endl;


                                      vmmhz1m_sum += vmmhz1m_tmp;


                                       stations[i].strings[j].antennas[k].vmmhz[ray_sol_cnt].push_back( vmmhz1m_tmp );


                                   }// end for freq bin
                                   
                                   cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"].vmmhz1m["<<ray_sol_cnt<<"][0] : "<<stations[i].strings[j].antennas[k].vmmhz[ray_sol_cnt][0]<<endl;
                                       

                                   ray_sol_cnt++;
                               
                               }// end while number of solutions
                           
                           }// end if solution exist

                           else {
                               
                               //cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"].trg = "<<stations[i].strings[j].antennas[k].trg[ray_sol_cnt]<<"  No vmmhz1m data!"<<endl;

                           }
                       
                       }// end if posnu selected

                           else {
                               cout<<" No posnu!!!!!! No signals calculated at all!!"<<endl;
                           }


                           stations[i].strings[j].antennas[k].ray_sol_cnt = ray_sol_cnt;    // save number of RaySolver solutions
                       

                   }// for number_of_antennas_string

               }// for number_of_strings_station

           }// for number_of_stations





}


Vector Report::GetPolarization (Vector &nnu, Vector &launch_vector) {
    // copy from icemc GetPolarization

    // Want to find a unit vector in the same plane as
    // nnu and launch_vector, but perpendicular to launch_vector, pointing away
    // from nnu.
    
    // cross nnu with launch_vector to get the direction of the B field.
    Vector n_bfield = nnu.Cross(launch_vector);
    
    // cross b-field with nrf2_iceside to get the polarization vector.
    Vector n_pol = n_bfield.Cross(launch_vector);
    
    n_pol = n_pol.Unit();
    
    // check and make sure E-field is pointing in the right direction.
    if (nnu*launch_vector>0 && n_pol*nnu>0)
	cout << "error in GetPolarization.  \n";
    
    
    
    return n_pol;
} //GetPolarization


void Report::GetParameters( Position &src, Position &trg, Vector &nnu, double &viewangle, double receive_angle, Vector &launch_vector, Vector &receive_vector) {

    viewangle = PI/2. - viewangle;  // viewangle was actually launch angle

    launch_vector = (trg.Cross( trg.Cross(src) )).Rotate(viewangle, trg.Cross(src));
    launch_vector = launch_vector.Unit();
    viewangle = launch_vector.Angle(nnu);

                                   cout<<"launch_vector angle between R1 (trg) : "<<launch_vector.Angle(trg)<<"\n";

    receive_vector = trg.Rotate( receive_angle, src.Cross(trg) );
    receive_vector = receive_vector.Unit();


}



