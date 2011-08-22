////////////////////////////////////////////////////////////////////////////////////////////////
//class Event:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EVENT_H
#define EVENT_H

#include <vector>
#include <string>
//#include "Position.h"


using std::vector;
using std::string;

class Position;

// Treat this like a struct and make data members public
class Event {

 private:
  
  protected:
  
  public:

  //  Event(int);
  //  Event(double);
 Event(double) : inu(0), nnu(Position()), pnu(1.E18), posnu(1,Position()), int_type(0), nu_flav(0), y(0.2), eField1m(3) {}
  ~Event();
  int inu; // index for this event
  Position nnu;  // unit vector pointing in direction of neutrino trajectory
  double pnu;  // neutrino energy


// direction and charge density of a charged lepton
  Vector sh_dir;

  //there should be same number of posnu's as charge_per_lengths
  //0th posnu is the hadronic shower
  //starting with 1st posnu are leptonic showers
  //there is one less charge_per_lengths than posnu's
  //0th charge_per_length is charge/length just after primary leptonic interaction
  //1st charge_per_length is just after first secondary leptonic interaction

  vector<Position> posnu;  // position of neutrino interaction
  vector<double> charge_per_length;
  vector<string> sh_type; // same length as posnu.
  // 0th element always indicates hadronic
  //choices are "em", "had"
  vector<double> sh_energy; // same length as posnu
  



  int int_type; // type of interaction - cc (0) or nc(1)
  int nu_flav;// neutrino flavor - e (0), mu (1), tau (2)
  double y; // inelasticity

  // something that describes the radiation that emerges from the shower
  // this may be a root file or text file that is generated from jaime's program
  vector<double> eField1m; // E field - 3d vector as a function of frequency at 1 m from the interaction
  
 


};



#endif //EVENT_H
