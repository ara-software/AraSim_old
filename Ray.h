////////////////////////////////////////////////////////////////////////////////////////////////
//class Ray:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef RAY_H
#define RAY_H

#include "Event.h"
#include <vector>
using namespace std;

class Ray {


 public:
  Ray();
  vector<double> getEField(Event *event,vector<double> pos); // Get E field from the shower at a particular location pos
  
  protected:
  
  private:

};

#endif //RAY_H

