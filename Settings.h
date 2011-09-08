#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <string>

using namespace std;
//using std::string;

class Settings 
{

    public:

        Settings();
        void Initialize();
        void ReadFile(string setupfile);

        int NNU;

  // NEED TO FIGURE OUT A GOOD WAY TO READ THIS IN AND STORE THEM.
  // INPUT FILE AGAIN?  SOMETHING ELSE?
  //These were moved here from IceModel under the new compilation scheme
        int ICE_MODEL; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP.
        int NOFZ; // 1=depth dependent index of refraction,0=off
        int CONSTANTCRUST; // set crust density and thickness to constant values.
        int CONSTANTICETHICKNESS; // set ice thickness to constant value
        int FIXEDELEVATION; // fix the elevation to the thickness of ice.
        int MOOREBAY; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data
        double EXPONENT; // 10^19 eV neutrinos only

};
#endif

