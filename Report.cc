#include "Detector.h"
#include "Report.h"
#include "Event.h"
#include "RaySolver.h"
#include "signal.hh"
#include "IceModel.h"
#include "Settings.h"
#include "Vector.h"
#include "Tools.h"

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

Report::Report(Detector *detector, Settings *settings1) {
    // Default constructor

    Initialize(detector, settings1);
    
}


void Report::Initialize(Detector *detector, Settings *settings1) {
    
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
    arrival_time.clear();
    reflection.clear();
    Pol_vector.clear();
    vmmhz.clear();
    VHz_antfactor.clear();
    VHz_filter.clear();
    Vfft.clear();
    Vfft_noise.clear();

    time.clear();
    Ax.clear();
    Ay.clear();
    Az.clear();

    V.clear();
    V_noise.clear();

    PeakV.clear();
    Rank.clear();

}


void Report::Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1) {

    int ray_sol_cnt;
    double viewangle;
    Position launch_vector; // direction of ray at the source
    Position receive_vector;    // direction of ray at the target antenna
    Vector n_trg_pokey;         // unit pokey vector at the target
    Vector n_trg_slappy;        // unit slappy vector at the target
    vector < vector <double> > ray_output;

    double vmmhz1m_tmp, vmmhz1m_sum, vmmhz1m_em;    // currently not using vmmhz1m_em
    Position Pol_vector;                            // polarization vector at the source
    double mag;                                     // magnification factor. it can vary in case of plane / spherical wave
    double fresnel;                                 // fresnel factor
    double tmp; // for non use return values

    double freq_tmp, heff, antenna_theta, antenna_phi;  // values needed for apply antenna gain factor and prepare fft, trigger
    double volts_forfft[settings1->NFOUR/2];       // array for fft
    double volts_noise_forfft[settings1->NFOUR/2]; // t domain noise voltage array


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
                                  cout<<"Path length : "<<ray_output[0][ray_sol_cnt]<<"\tView angle : "<<ray_output[1][ray_sol_cnt]<<"\tReceipt angle : "<<ray_output[2][ray_sol_cnt]<<"\tpath time : "<<ray_output[4][ray_sol_cnt]<<endl;

                                  // set viewangle, launch_vector, receive vectors
                                  viewangle = ray_output[1][ray_sol_cnt];
                                  GetParameters ( event->Nu_Interaction[0].posnu,   // posnu
                                                detector->stations[i].strings[j].antennas[k],   // trg antenna
                                                event->nnu,                         // nnu
                                                viewangle,         // inputs launch_angle, returns viewangle
                                                ray_output[2][ray_sol_cnt],         // receive_angle
                                                launch_vector, receive_vector,
                                                n_trg_slappy, n_trg_pokey);


                                  // check viewangle that if ray in near Cherenkov cone
                                  //
                                  if (viewangle*DEGRAD >55. && viewangle*DEGRAD <57.) { // if viewangle is 56 deg +- 1 deg
                                      cout<<"near cone! view angle : "<<viewangle*DEGRAD<<"station["<<i<<"].string["<<j<<"].antenna["<<k<<"] with  ray_sol_cnt : "<<ray_sol_cnt<<endl;
                                  }

                                   // store information to report
                                   stations[i].strings[j].antennas[k].view_ang.push_back(viewangle);
                                   stations[i].strings[j].antennas[k].rec_ang.push_back(ray_output[2][ray_sol_cnt]);
                                   stations[i].strings[j].antennas[k].Dist.push_back(ray_output[0][ray_sol_cnt]);
                                   stations[i].strings[j].antennas[k].arrival_time.push_back(ray_output[4][ray_sol_cnt]);
                                   stations[i].strings[j].antennas[k].reflect_ang.push_back(ray_output[3][ray_sol_cnt]);
                                   
                                   stations[i].strings[j].antennas[k].vmmhz.resize(ray_sol_cnt+1);
                                   
                                   stations[i].strings[j].antennas[k].VHz_antfactor.resize(ray_sol_cnt+1);
                                   stations[i].strings[j].antennas[k].VHz_filter.resize(ray_sol_cnt+1);
                                   stations[i].strings[j].antennas[k].Vfft.resize(ray_sol_cnt+1);
                                   stations[i].strings[j].antennas[k].Vfft_noise.resize(ray_sol_cnt+1);

                                   stations[i].strings[j].antennas[k].V.resize(ray_sol_cnt+1);
                                   stations[i].strings[j].antennas[k].V_noise.resize(ray_sol_cnt+1);
                                   stations[i].strings[j].antennas[k].time.resize(ray_sol_cnt+1);

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


                                   GetAngleAnt(receive_vector, detector->stations[i].strings[j].antennas[k], antenna_theta, antenna_phi);   // get theta, phi for signal ray arrived at antenna
                                   cout<<"antenna theta : "<<antenna_theta<<"  phi : "<<antenna_phi<<endl;  




                                   for (int l=0; l<detector->GetFreqBin(); l++) {   // for detector freq bin numbers


                                      cout<<"TaperVmMHz inputs VA:"<<viewangle<<" th_em:"<<event->Nu_Interaction[0].d_theta_em[l]<<" th_had:"<<event->Nu_Interaction[0].d_theta_had[l]<<" emfrac:"<<event->Nu_Interaction[0].emfrac<<" hadfrac:"<<event->Nu_Interaction[0].hadfrac<<" vmmhz1m:"<<event->Nu_Interaction[0].vmmhz1m[l]<<endl;

                                       vmmhz1m_tmp = event->Nu_Interaction[0].vmmhz1m[l];

                                       signal->TaperVmMHz( viewangle, event->Nu_Interaction[0].d_theta_em[l], event->Nu_Interaction[0].d_theta_had[l], event->Nu_Interaction[0].emfrac, event->Nu_Interaction[0].hadfrac, vmmhz1m_tmp, vmmhz1m_em);
                                      cout<<"TaperVmMHz (1m at view angle) at "<<l<<"th bin : "<<vmmhz1m_tmp<<endl;


                                      // multiply all factors for traveling ice
                                      vmmhz1m_tmp = vmmhz1m_tmp / ray_output[0][ray_sol_cnt] * exp(-ray_output[0][ray_sol_cnt]/icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0)) * mag * fresnel;  // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
                                      cout<<"AttenLength : "<<icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0)<<endl;


                                      vmmhz1m_sum += vmmhz1m_tmp;


                                       stations[i].strings[j].antennas[k].vmmhz[ray_sol_cnt].push_back( vmmhz1m_tmp );


                                       freq_tmp = detector->GetFreq(l); // freq in Hz

                                       heff = GaintoHeight(detector->stations[i].strings[j].antennas[k].GetG(detector, freq_tmp*1.E-6, // to MHz
                                                   antenna_theta, antenna_phi), 
                                                   freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]) );
                                       cout<<"n_medium : "<<icemodel->GetN(detector->stations[i].strings[j].antennas[k])<<endl;
                                       cout<<"gain : "<<detector->stations[i].strings[j].antennas[k].GetG(detector, freq_tmp*1.E-6, antenna_theta, antenna_phi)<<endl;
                                       cout<<"heff : "<<heff<<endl;

                                       // apply pol factor, heff
                                       ApplyAntFactors(heff, n_trg_pokey, n_trg_slappy, Pol_vector, detector->stations[i].strings[j].antennas[k].type, vmmhz1m_tmp);

                                       //stations[i].strings[j].antennas[k].Vfft[ray_sol_cnt].push_back( vmmhz1m_tmp );
                                       stations[i].strings[j].antennas[k].VHz_antfactor[ray_sol_cnt].push_back( vmmhz1m_tmp );

                                       // apply filter
                                       ApplyFilter(l, detector, vmmhz1m_tmp);

                                       //stations[i].strings[j].antennas[k].Vfft[ray_sol_cnt].push_back( vmmhz1m_tmp );
                                       stations[i].strings[j].antennas[k].VHz_filter[ray_sol_cnt].push_back( vmmhz1m_tmp );




                                   }// end for freq bin
                                   
                                   //cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"].vmmhz1m["<<ray_sol_cnt<<"][0] : "<<stations[i].strings[j].antennas[k].vmmhz[ray_sol_cnt][0]<<endl;

                                   MakeArraysforFFT(settings1, detector, stations[i].strings[j].antennas[k].VHz_filter[ray_sol_cnt], volts_forfft);


                                   // save freq domain array which is prepaired for realft
                                   for (int n=0; n<settings1->NFOUR/2; n++) {
                                       stations[i].strings[j].antennas[k].Vfft[ray_sol_cnt].push_back( volts_forfft[n] );
                                   }




                                   // now, after realft, volts_forfft is time domain signal at backend of antenna
                                   Tools::realft(volts_forfft,-1,settings1->NFOUR/2);
                                   //Tools::realft(volts_forfft,1,settings1->NFOUR/2);


                                   cout<<"Finished getting V signal part!!"<<endl;
                                   
                                   stations[i].strings[j].antennas[k].PeakV.push_back( FindPeak(volts_forfft, settings1->NFOUR/2) );

                                   // test; noise from flat spectrum Kb * T
                                   // T_ice = 240K
                                   double T_noise = 240.;
                                   Vfft_noise_org = sqrt( (double)(settings1->NFOUR/2) * 50. * KBOLTZ * T_noise / (settings1->TIMESTEP * 2.) );
                                   // Vfft_noise_org is in fft freq bin!!
                                   // same unit with Vfft [V] but filter not applied

                                   

                                   cout<<"start GetNoiseWaveforms!!"<<endl;
                                   // from Vfft_noise_org (flat thermal noise), get random noise waveform
                                   // GetNoiseWaveforms will apply filter in it.
                                   GetNoiseWaveforms(settings1, detector, Vfft_noise_org, volts_noise_forfft);
                                   cout<<"Finished GetNoiseWaveforms!!"<<endl;

                                   Tools::NormalTimeOrdering(settings1->NFOUR/2, volts_forfft);
                                   Tools::NormalTimeOrdering(settings1->NFOUR/2, volts_noise_forfft);
                                   cout<<"finished NormalTimeOrdering!!"<<endl;


                                   for (int n=0; n<settings1->NFOUR/2; n++) {
                                       stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back( volts_forfft[n] );
                                       stations[i].strings[j].antennas[k].V_noise[ray_sol_cnt].push_back( volts_noise_forfft[n] );
                                       stations[i].strings[j].antennas[k].time[ray_sol_cnt].push_back( stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt] + (double)(n - settings1->NFOUR/4)* settings1->TIMESTEP );   // time at 0 s is when ray started at the posnu
                                   }
                                   cout<<"finished push_back V, V_noise, and time!!"<<endl;
                                   
                                       

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


           // after all values are stored in Report, set ranking of signal between antennas
           SetRank(detector);



}   // end Connect_Interaction_Detector


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


void Report::GetParameters( Position &src, Position &trg, Vector &nnu, double &viewangle, double receive_angle, Vector &launch_vector, Vector &receive_vector, Vector &n_trg_slappy, Vector &n_trg_pokey) {

    viewangle = PI/2. - viewangle;  // viewangle was actually launch angle

    launch_vector = (trg.Cross( trg.Cross(src) )).Rotate(viewangle, trg.Cross(src));
    launch_vector = launch_vector.Unit();
    viewangle = launch_vector.Angle(nnu);

                                   cout<<"launch_vector angle between R1 (trg) : "<<launch_vector.Angle(trg)<<"\n";

    receive_vector = trg.Rotate( receive_angle, src.Cross(trg) );
    receive_vector = receive_vector.Unit();

    n_trg_pokey = trg.Unit();
    n_trg_slappy = (trg.Cross(src)).Unit();


}


double Report::GaintoHeight(double gain, double freq, double n_medium) {


    // from gain=4*pi*A_eff/lambda^2
    // and h_eff=2*sqrt(A_eff*Z_rx/Z_air)
    // gain is unitless value
    
    return 2*sqrt(gain/4/PI*CLIGHT*CLIGHT/(freq*freq*n_medium*n_medium)*Zr/(Z0/n_medium));  // n_medium parts are changed from icemc(I believe this is correct one; E. Hong)
}


void Report::ApplyAntFactors(double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &vmmhz) {  // vmmhz is input and output. output will have some antenna factors on it

    double pol_factor;

    if (ant_type == 0) {    // if v pol
        pol_factor = n_trg_pokey * Pol_vector;
    }
    else if (ant_type == 1) {   // if h pol
        pol_factor = n_trg_slappy * Pol_vector;
    }

    // apply 3dB spliter, d nu to prepare FFT
    // now actually vmmhz is not V/m/MHz but V/m/Hz unit
    //vmmhz = vmmhz/sqrt(2.)/(settings1->TIMESTEP*1.E6); //sqrt(2) for 3dB spliter for TURF, SURF
    vmmhz = vmmhz/sqrt(2.)/(1.E6); //sqrt(2) for 3dB spliter for TURF, SURF. 1/TIMESTEP moved to MakeArraysforFFT
    // 1/(1.E6) for V/MHz to V/Hz


    // apply antenna effective height and 0.5 (to calculate power with heff), and polarization factor
    // not vmmhz is actually V/Hz unit
    vmmhz = vmmhz * 0.5 * heff * abs(pol_factor);

    // now we have to use MakeArraysforFFT to change signal arrays for FFT
}


void Report::ApplyFilter(int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetFilterGain(bin_n) )/20.);   // from dB to unitless gain for voltage

}


void Report::ApplyFilter_fft(int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetFilterGain_fft(bin_n) )/20.);   // from dB to unitless gain for voltage

}



void Report::GetAngleAnt(Vector &rec_vector, Position &antenna, double &ant_theta, double &ant_phi) {   //ant_theta and ant_phi is in degree 

    // need to fix some parts.
    // currently phi is not correct.
    // I think we have to actually rotate (0,-1,0) to antenna location and get phi from rotated (0,-1,0)
    ant_theta = rec_vector.Angle(antenna)*DEGRAD;

    Vector unit0_10;
    Vector unit100;
    unit0_10.SetXYZ(0,-1,0);
    unit100.SetXYZ(1,0,0);
    Vector proj;

    proj = antenna.Cross( rec_vector.Cross(antenna) );

    if (proj * unit100 > 0.) {  // rec_vector is pointing to positive x direction
        ant_phi = 360. - ( proj.Angle(unit0_10) )*DEGRAD;
    }
    else {
        ant_phi = proj.Angle(unit0_10) * DEGRAD;
    }

}


void Report::GetNoiseWaveforms(Settings *settings1, Detector *detector, double v_noise, double *vnoise) {

    if (settings1->NOISE == 0) {    // NOISE == 0 : flat thermal noise with Johnson-Nyquist noise
        //V_noise_fft_bin = sqrt( (double)(NFOUR/2) * 50. * KBOLTZ * T_noise / (2. * TIMESTEP) ); 

        Vfft_noise_after.clear();  // remove previous Vfft_noise values
        Vfft_noise_before.clear();  // remove previous Vfft_noise values


        double V_tmp; // copy original flat H_n [V] value
        double current_amplitude, current_phase;

        GetNoisePhase(settings1); // get random phase for noise

        //MakeArraysforFFT_noise(settings1, detector, vhz_noise_tmp , vnoise);
        // returned array vnoise currently have real value = imag value as no phase term applied

        for (int k=0; k<settings1->NFOUR/4; k++) {

            V_tmp = v_noise; // copy original flat H_n [V] value

            ApplyFilter_fft(k, detector, V_tmp);
            Vfft_noise_before.push_back( V_tmp );


            current_phase = noise_phase[k];

            Tools::get_random_rician( 0., 0., sqrt(2./ M_PI) * V_tmp, current_amplitude, current_phase);    // use real value array value

            // vnoise is currently noise spectrum (before fft, unit : V)
           //vnoise[2 * k] = sqrt(current_amplitude) * cos(noise_phase[k]);
           //vnoise[2 * k + 1] = sqrt(current_amplitude) * sin(noise_phase[k]);
           vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
           vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);



           //vnoise[2 * k] = (V_tmp) * cos(noise_phase[k]);
           //vnoise[2 * k + 1] = (V_tmp) * sin(noise_phase[k]);
           

            Vfft_noise_after.push_back( vnoise[2*k] );
            Vfft_noise_after.push_back( vnoise[2*k+1] );

            // inverse FFT normalization factor!
            vnoise[2 * k] *= 2./((double)settings1->NFOUR/2);
            vnoise[2 * k + 1] *= 2./((double)settings1->NFOUR/2);


        }


        // now vnoise is time domain waveform
        Tools::realft( vnoise, -1, settings1->NFOUR/2);

    }
    else {  // currently there are no more options!
        cout<<"no noise option for NOISE = "<<settings1->NOISE<<endl;
    }


}


void Report::GetNoisePhase(Settings *settings1) {

    noise_phase.clear();    // remove all previous noise_phase vector values

    for (int i=0; i<settings1->NFOUR/4; i++) {
        noise_phase.push_back(2*PI*gRandom->Rndm());  // get random phase for flat thermal noise
    }
}


    

void Report::MakeArraysforFFT(Settings *settings1, Detector *detector, vector <double> &vsignal_array, double *vsignal_forfft) {



    // from icemc, anita class MakeArraysforFFT

    int NFOUR = settings1->NFOUR;
    Tools::Zero(vsignal_forfft,NFOUR/2);
    
    double previous_value_e_even=0.;
    double previous_value_e_odd=0.;
    int count_nonzero=0;
    int iprevious=0;
    int ifirstnonzero=-1;
    int ilastnonzero=2000;
    //for (int i=0;i<NFREQ;i++) {
    for (int i=0;i<detector->GetFreqBin();i++) {
	
	// freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
	// but there are only NFOUR/4 different values
	// it's the index among the NFOUR/4 that we're interested in
	int ifour=Tools::Getifreq(detector->GetFreq(i),detector->freq_forfft[0],detector->freq_forfft[NFOUR/2-1],NFOUR/4);
	
	if (ifour!=-1 && 2*ifour+1<NFOUR/2) {
	    count_nonzero++;
	    if (ifirstnonzero==-1)
		ifirstnonzero=ifour;
	    
	    vsignal_forfft[2*ifour]=vsignal_array[i]*2/((double)NFOUR/2)/(settings1->TIMESTEP); // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft
	    
	    //      cout << "ifour, vsignal is " << ifour << " " << vsignal_e_forfft[2*ifour] << "\n";
	    
	    vsignal_forfft[2*ifour+1]=vsignal_array[i]*2/((double)NFOUR/2)/(settings1->TIMESTEP); // phase is 90 deg.
	    // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting
	    
	    // how about we interpolate instead of doing a box average
	    
	    for (int j=iprevious+1;j<ifour;j++) {
		vsignal_forfft[2*j]=previous_value_e_even+(vsignal_forfft[2*ifour]-previous_value_e_even)*(double)(j-iprevious)/(double)(ifour-iprevious);
		//	cout << "j, vsignal is " << j << " " << vsignal_e_forfft[2*j] << "\n";
		
		vsignal_forfft[2*j+1]=previous_value_e_odd+(vsignal_forfft[2*ifour+1]-previous_value_e_odd)*(double)(j-iprevious)/(double)(ifour-iprevious);
	    }
	    
	    ilastnonzero=ifour;
	    iprevious=ifour;
	    previous_value_e_even=vsignal_forfft[2*ifour];
	    previous_value_e_odd=vsignal_forfft[2*ifour+1];
	}
	
    } // end loop over nfreq
    


    // don't applying any extra factor for the change in array (change of bin size)
    // as change in the bin size doesn't matter for the total energy
    // total energy is just same as integral over frequency range and this frequency range will not change
    //
      //cout << "ifirstnonzero, ilastnonzero are " << ifirstnonzero << " " << ilastnonzero << "\n";
      //cout << "non zero count is " << count_nonzero << "\n";
      //cout << "ratio is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
    for (int j=0;j<NFOUR/4;j++) {
	vsignal_forfft[2*j]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
	vsignal_forfft[2*j+1]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
    }



    
    //  Tools::InterpolateComplex(vsignal_e_forfft,NFOUR/4);
    //Tools::InterpolateComplex(vsignal_h_forfft,NFOUR/4);
    for (int ifour=0;ifour<NFOUR/4;ifour++) {



	vsignal_forfft[2*ifour]*=cos(settings1->PHASE*PI/180.);
	vsignal_forfft[2*ifour+1]*=sin(settings1->PHASE*PI/180.);
	
	//--------------------------------------------------
	// if (!PULSER) {
	//     
	//     vsignal_forfft[2*ifour]*=cos(phase*PI/180.);
	//     vsignal_forfft[2*ifour+1]*=sin(phase*PI/180.);
	//     
	//     
	// }
	// else {
	//     vsignal_forfft[2*ifour]*=cos(v_phases[ifour]*PI/180.);
	//     vsignal_forfft[2*ifour+1]*=sin(v_phases[ifour]*PI/180.);
	//     
	// }	  	  	  
	//-------------------------------------------------- 
	
	
    }
    
    
}




void Report::MakeArraysforFFT_noise(Settings *settings1, Detector *detector, vector <double> &vsignal_array, double *vsignal_forfft) {



    // from icemc, anita class MakeArraysforFFT

    int NFOUR = settings1->NFOUR;
    Tools::Zero(vsignal_forfft,NFOUR/2);
    
    double previous_value_e_even=0.;
    double previous_value_e_odd=0.;
    int count_nonzero=0;
    int iprevious=0;
    int ifirstnonzero=-1;
    int ilastnonzero=2000;
    //for (int i=0;i<NFREQ;i++) {
    for (int i=0;i<detector->GetFreqBin();i++) {
	
	// freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
	// but there are only NFOUR/4 different values
	// it's the index among the NFOUR/4 that we're interested in
	int ifour=Tools::Getifreq(detector->GetFreq(i),detector->freq_forfft[0],detector->freq_forfft[NFOUR/2-1],NFOUR/4);
	
	if (ifour!=-1 && 2*ifour+1<NFOUR/2) {
	    count_nonzero++;
	    if (ifirstnonzero==-1)
		ifirstnonzero=ifour;
	    
	    vsignal_forfft[2*ifour]=vsignal_array[i]*2/((double)NFOUR/2)/(settings1->TIMESTEP); // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft
	    
	    //      cout << "ifour, vsignal is " << ifour << " " << vsignal_e_forfft[2*ifour] << "\n";
	    
	    vsignal_forfft[2*ifour+1]=vsignal_array[i]*2/((double)NFOUR/2)/(settings1->TIMESTEP); // phase is 90 deg.
	    // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting
	    
	    // how about we interpolate instead of doing a box average
	    
	    for (int j=iprevious+1;j<ifour;j++) {
		vsignal_forfft[2*j]=previous_value_e_even+(vsignal_forfft[2*ifour]-previous_value_e_even)*(double)(j-iprevious)/(double)(ifour-iprevious);
		//	cout << "j, vsignal is " << j << " " << vsignal_e_forfft[2*j] << "\n";
		
		vsignal_forfft[2*j+1]=previous_value_e_odd+(vsignal_forfft[2*ifour+1]-previous_value_e_odd)*(double)(j-iprevious)/(double)(ifour-iprevious);
	    }
	    
	    ilastnonzero=ifour;
	    iprevious=ifour;
	    previous_value_e_even=vsignal_forfft[2*ifour];
	    previous_value_e_odd=vsignal_forfft[2*ifour+1];
	}
	
    } // end loop over nfreq
    


    // don't applying any extra factor for the change in array (change of bin size)
    // as change in the bin size doesn't matter for the total energy
    // total energy is just same as integral over frequency range and this frequency range will not change
    //
      //cout << "ifirstnonzero, ilastnonzero are " << ifirstnonzero << " " << ilastnonzero << "\n";
      //cout << "non zero count is " << count_nonzero << "\n";
      //cout << "ratio is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
    for (int j=0;j<NFOUR/4;j++) {
	vsignal_forfft[2*j]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
	vsignal_forfft[2*j+1]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
    }


    // phase for noise will be applied in GetNoiseWaveforms

    
}



double Report::FindPeak(double *waveform, int n) {  // same with icemc, trigger-> AntTrigger::FindPeak

    double peak=0.;

    for (int i=0;i<n;i++) {
        if (fabs(waveform[i])>peak)
            peak = fabs(waveform[i]);
    }
    return peak;
}

void Report::SetRank(Detector *detector) {
    // think about this!
    // test 1

    int current = 0;    // current checking rank
    int check;  // check is any change in rank
    double maxpeak; // maxpeak value
    double pre_maxpeak; // maxpeak at previous round


        for (int i = 0; i< detector->params.number_of_stations; i++) {

            for (int j=0; j< detector->params.number_of_strings_station; j++) {

                for (int k=0; k< detector->params.number_of_antennas_string; k++) {

                    if (stations[i].strings[j].antennas[k].PeakV[0] == 0.) {
                        stations[i].strings[j].antennas[k].Rank.push_back(0);  // rank 0, PeakV is 0, non-countable rank.
                    }
                    else {
                        stations[i].strings[j].antennas[k].Rank.push_back(current+1);  // elses, if PeakV is not 0, set Rank as 1 (for now)
                    }
                }
            }
        }


    cout<<"Start while loop for Ranking!!"<<endl;

    check = 1;
    pre_maxpeak = 1.E5; // unreasonably big value which real PeakV will never reach
    while (check!=0) {
        check=0;
        maxpeak = 0.;
        for (int i = 0; i< detector->params.number_of_stations; i++) {

            for (int j=0; j< detector->params.number_of_strings_station; j++) {

                for (int k=0; k< detector->params.number_of_antennas_string; k++) {

                    //for (int l=0; l< detector->stations[i].strings[j].antennas[k].ray_sol_cnt; l++) {
                    //
                    //lets start with 1st ray_sol only
                    //

                    if (stations[i].strings[j].antennas[k].Rank[0] != 0) {  // there is non-zero value!
                    if (stations[i].strings[j].antennas[k].PeakV[0] < pre_maxpeak) {

                        if (maxpeak < stations[i].strings[j].antennas[k].PeakV[0] ) {
                            maxpeak = stations[i].strings[j].antennas[k].PeakV[0];

                            stations[i].strings[j].antennas[k].Rank[0] = current+1;
                            check++;
                        } // is maxpeak < PeakV
                        else {
                            stations[i].strings[j].antennas[k].Rank[0] = current + 2;
                        } // else

                    } // if rank > current
                    }

                } // for antennas
            } // for strings
        } // for stations

        current++;
        pre_maxpeak = maxpeak;


    }   // while


    cout<<"END Ranking!!"<<endl;




}






