////////////////////////////////////////////////////////////////////////////////////////////////
//class Ray:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef RAY_H
#define RAY_H

#include <cmath>
#include <iostream>

class Event;
class Position;
class Vector;
class IceModel;

using std::vector;
using std::cout;


class Ray {


 public:
  Ray();
  vector<double> getEField(Event *event,Position pos); // Get E field from the shower at a particular location pos


  //
  // below WhereDoesItLeave member is copied from ray.hh in icemc.
  static int WhereDoesItLeave(int inu,const Position &posnu,
			      const Vector &ntemp,IceModel *antarctica,
			      Position &r_out) {
    
    double distance=0;
    double posnu_length=posnu.Mag(); // distance from center of earth to interaction
    
    double lon,lat,lon_old,lat_old; //latitude, longitude indices for 1st and 2nd iteration
    lon = posnu.Lon(); // what latitude, longitude does interaction occur at
    lat = posnu.Lat();
    lon_old=lon; // save this longitude and latitude so we can refer to it later
    lat_old=lat;
    
    // use law of cosines to get distance from interaction to exit point for the ray
    // need to solve for that distance using the quadratic formula
    
    // angle between posnu and ntemp vector for law of cosines.
    double costheta=-1*(posnu*ntemp)/posnu_length;
    
    // a,b,c for quadratic formula, needed to solve for 
    double a=1;
    double b=-1*2*posnu_length*costheta;
    double c=posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2);

    
    if (b*b-4*a*c<0.) {
   
      return 0;
    }
    // else
//       cout << "positive.  c is " << c << "\n";
   

    // use the "+" solution because the other one is where the ray is headed downward toward the rock
    distance=(-1*b+sqrt(b*b-4*a*c))/2;
    
    
    // now here is the exit point for the ray
    r_out = posnu + distance*ntemp;
    
    lon = r_out.Lon(); // latitude and longitude of exit point
    lat = r_out.Lat();
    
    
    c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
    // sometimes though the new surface is lower than the one over posnu which causes a problem.
    if (b*b-4*a*c<0.) {
      //cout << "inu is " << inu << "\n";  
      // try halving the distance
      distance=distance/2.;
      //cout << "bad.  distance 1/2 is " << distance << "\n";
      r_out = posnu + distance*ntemp;
      lon = r_out.Lon(); // latitude and longitude of exit point
      lat = r_out.Lat();
      c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
      if (b*b-4*a*c<0.) { // if we still have the problem back up more
	distance=distance/2.; // now we are at 1/4 the distance
	//cout << "bad.  distance 1/4 is " << distance << "\n";
	r_out = posnu + distance*ntemp;
	lon = r_out.Lon(); // latitude and longitude of exit point
	lat = r_out.Lat();
	c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	if (b*b-4*a*c<0.) { // the problem is less then 1/4 of the way in
	  
	  distance=distance/2.; // now we are at 1/8 the distance
	  //cout << "bad.  distance 1/8 is " << distance << "\n";
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines

	  if (b*b-4*a*c<0.) {
	    // still have the problem so just make the distance 0
	    distance=0.; // now we are at 1/8 the distance
	    //cout << "bad.  distance is " << distance << "\n";
	    lon = posnu.Lon(); // latitude and longitude of exit point
	    lat = posnu.Lat();
	    r_out=antarctica->Surface(lon,lat)/posnu.Mag()*posnu;
	  }
// 	  else
// 	    cout << "good.\n";





	} // now we are at 1/8 the distance
	else {// if this surface is ok problem is between 1/4 and 1/2
	  distance=distance*1.5; // now we are at 3/8 the distance
	  //	cout << "good.  distance 3/8 is " << distance << "\n";
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines

	  if (b*b-4.*a*c<0.) {
	    distance=distance*2./3.; // go back to 1/4
	    r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	  //cout << "good at distance 1/4 is " << distance << "\n";
	  }
// 	  else
// 	    cout << "good.\n";

	} // now we are at 3/8 the distance
      

      } // now we are at 1/4 the distance
      else { // if this surface at 1/2 distance is ok see if we can go a little further
	distance=distance*1.5; // now we are at 3/4 the distance
	//cout << "good.  distance 3/4 is " << distance << "\n";
	r_out = posnu + distance*ntemp;
	lon = r_out.Lon(); // latitude and longitude of exit point
	lat = r_out.Lat();
	c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	if (b*b-4*a*c<0.) { // the problem is between 1/2 and 3/4 of the way in
	  
	  distance=distance*5./6.; // now we are at 5/8 the distance
	  //cout << "bad.  distance 5/8 is " << distance << "\n";
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines

	  if (b*b-4*a*c<0.) {
	    distance=distance*4./5.;
 	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	  //cout << "good at distance 3/4 is " << distance << "\n";
	  }
// 	  else 
// 	    cout << "good at distance 5/8.\n";



	} // now we are at 1/8 the distance
	else {// if this surface is ok problem is between 1/4 and 1/2
	  distance=distance*7./6.; // now we are at 7/8 the distance
	  //cout << "good.  distance 7/8 is " << distance << "\n";
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines

	  if (b*b-4*a*c<0) {
	    // now found the problem so go back to 3/4 distance
	    distance=distance*6./7.;
	    //cout << "good at  distance 3/4 is " << distance << "\n";
	    r_out = posnu + distance*ntemp;
	    lon = r_out.Lon(); // latitude and longitude of exit point
	    lat = r_out.Lat();
	    c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines


	  }
// 	  else
// 	    cout << "good.\n";

	} // now we are at 3/8 the distance







      } // now we are at 3/4 distance
    
    
    } // if exit point we initially found was not ok
    else {
      distance=(-1*b+sqrt(b*b-4*a*c))/2; // and quadratic formula
      r_out = posnu + distance*ntemp;
    }
    if (inu==1159426) {
      cout << "a,b,c,b*b-4*a*c are " << a << " " << b << " " << c << " " << b*b-4*a*c << "\n";
      cout << "posnu is ";posnu.Print();
      cout << "distance is " << distance << "\n";
      cout << "ntemp is ";ntemp.Print();
    }
    

    
    
    
    return 1;
  }

  
  protected:
  
  private:

};

#endif //RAY_H

