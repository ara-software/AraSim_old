#include "Event.h"
#include <string>
#include "Settings.h"
#include "Spectra.h"
#include "Primaries.h"
#include "Report.h"

ClassImp(Event);

using namespace std;


Event::Event () {
    //default constructor
    Initialize ();
}

void Event::Initialize () {

    Nu_Interaction.clear();
    //test_report.clear();
    n_interactions = 1;  // total number of interaction (including secondaries)

}


Event::Event (Settings *settings1, Spectra *spectra1, Primaries *primary1, IceModel *icemodel, Detector *detector, Signal *signal, Secondaries *sec1 ) {

    Initialize ();

    Choose_Evt_Type (settings1);

    if (Event_type == 0) { // if only neutrino events exist

        pnu = spectra1->GetNuEnergy();
        nuflavor = primary1->GetNuFlavor();

        if (settings1->NNU_THIS_THETA==1) {    // set specific theta angle for nnu
            nnu = primary1->GetThatDirection(settings1->NNU_THETA, settings1->NNU_D_THETA);
        }
        else { // nnu angle random
            nnu = primary1->GetAnyDirection();
        }
        
        if (nuflavor=="nue")
            nuflavorint=1;
        else if (nuflavor=="numu")
            nuflavorint=2;
        else if (nuflavor=="nutau")
            nuflavorint=3;

        Interaction *Nu_temp;
        //Report *report_tmp;

        Nu_temp = new Interaction (pnu, nnu, nuflavor, n_interactions, icemodel, detector, settings1, primary1, signal, sec1 );
        //report_tmp = new Report(detector ,settings1);

        Nu_Interaction.push_back(*Nu_temp);  // for the first interaction
        //test_report.push_back(*report_tmp);

        delete Nu_temp;

        // for multiple interactions...
/*
        while (interaction_count < n_interactions) {    // not sure if this will work???

            Nu_tmp = new Interaction (...., n_interactions );

            Nu_Interaction.push_back( Nu_tmp );
        }
*/

        
    }





}

Event::~Event() {
    Nu_Interaction.clear();
}



void Event::Choose_Evt_Type (Settings *settings1) {

    if (settings1->EVENT_TYPE==0) {
        //cout<<"Only Neutrino Evnets!"<<endl;
        Event_type = 0;
    }
    else {
        cout<<"Currently, only neutrino events possible!"<<endl;
        cout<<"Change Evt_type from "<<settings1->EVENT_TYPE<<" to 0"<<endl;
        Event_type = 0;
    }

}


