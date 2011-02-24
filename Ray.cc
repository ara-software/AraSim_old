#include "Ray.h"
#include "Tools.h"

Ray::Ray() {

}
vector<double> Ray::getEField(Event *event,Position pos) {

  double d=event->posnu.Distance(pos);
  
  vector<double> efield;
  efield.push_back(event->eField1m[0]/d);
  efield.push_back(event->eField1m[1]/d);
  efield.push_back(event->eField1m[2]/d);
  
  return efield;
  
}
