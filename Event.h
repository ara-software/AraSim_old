////////////////////////////////////////////////////////////////////////////////////////////////
//class Event:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EVENT_H
#define EVENT_H

#include <vector>
#include "Position.h"


using namespace std;

class Position;

// Treat this like a struct and make data members public
class Event {

 private:
  
  protected:
  
  public:

  //  Event(int);
  //  Event(double);
  ~Event();
  int inu; // index for this event
  Position nnu;  // unit vector pointing in direction of neutrino trajectory
  double pnu;  // neutrino energy
  Position posnu;  // position of neutrino interaction

  int int_type; // type of interaction - cc (0) or nc(1)
  int nu_flav;// neutrino flavor - e (0), mu (1), tau (2)
  double y; // inelasticity

  // something that describes the radiation that emerges from the shower
  // this may be a root file or text file that is generated from jaime's program
  vector<double> eField1m; // E field - 3d vector as a function of frequency at 1 m from the interaction
  
 
 Event(double) : inu(0), nnu(Position()), pnu(1.E18), posnu(Position()), int_type(0), nu_flav(0), y(0.2), eField1m(3) {}


};

#endif //EVENT_H
