#include "IceModel.h"
#include "EarthModel.h"
#include "Constants.h"
#include "Vector.h"

#include <iostream>
#include <fstream>
using namespace std;




//Parameters of the BEDMAP ice model. (See http://www.antarctica.ac.uk/aedc/bedmap/download/)
int nCols_ice=1200; //number of columns in data, set by header file (should be 1200)
int nRows_ice=1000; //number of rows in data, set by header file (should be 1000)
int cellSize=5000; //in meters, set by header file (should be 5000) - same for both files
int xLowerLeft_ice=-3000000; 
int yLowerLeft_ice=-2500000;
int nCols_ground=1068;
int nRows_ground=869;
int xLowerLeft_ground=-2661600;
int yLowerLeft_ground=-2149967;
int nCols_water=1200;
int nRows_water=1000;
int xLowerLeft_water=-3000000;
int yLowerLeft_water=-2500000;
int NODATA=-9999;

//Variables for conversion between BEDMAP polar stereographic coordinates and lat/lon.  Conversion equations from ftp://164.214.2.65/pub/gig/tm8358.2/TM8358_2.pdf
const double scale_factor=0.97276901289;  //scale factor at pole corresponding to 71 deg S latitude of true scale (used in BEDMAP)
const double ellipsoid_inv_f = 298.257223563; //of Earth
const double ellipsoid_b = R_EARTH*(1-(1/ellipsoid_inv_f));
const double eccentricity = sqrt((1/ellipsoid_inv_f)*(2-(1/ellipsoid_inv_f)));
const double bedmap_a_bar = pow(eccentricity,2)/2 + 5*pow(eccentricity,4)/24 + pow(eccentricity,6)/12 + 13*pow(eccentricity,8)/360;
const double bedmap_b_bar = 7*pow(eccentricity,4)/48 + 29*pow(eccentricity,6)/240 + 811*pow(eccentricity,8)/11520;
const double bedmap_c_bar = 7*pow(eccentricity,6)/120 + 81*pow(eccentricity,8)/1120;
const double bedmap_d_bar = 4279*pow(eccentricity,8)/161280;
const double bedmap_c_0 = (2*R_EARTH / sqrt(1-pow(eccentricity,2))) * pow(( (1-eccentricity) / (1+eccentricity) ),eccentricity/2);
double bedmap_R = scale_factor*bedmap_c_0 * pow(( (1 + eccentricity*sin(71*RADDEG)) / (1 - eccentricity*sin(71*RADDEG)) ),eccentricity/2) * tan((PI/4) - (71*RADDEG)/2); //varies with latitude, defined here for 71 deg S latitude
const double bedmap_nu = bedmap_R / cos(71*RADDEG);


IceModel::IceModel(int model,int earth_model,int moorebay) : EarthModel(earth_model),mooreBayFlag(moorebay) {

  
  setUpIceModel(model);

 


 }

void IceModel::setUpIceModel(int model) {
  
  DEPTH_DEPENDENT_N = (int) (model / 10);
  model -= DEPTH_DEPENDENT_N * 10;
  ice_model=model;

  if (ice_model != 0 && ice_model != 1) {
    cout<<"Error!  Unknown ice model requested!  Defaulting to Crust 2.0.\n";
    ice_model = 0;
  } //if
  else if (ice_model==1) {
    ReadIceThickness();
    ReadWaterDepth();
    ReadGroundBed();
  } //else if (BEDMAP)
 //read in attenuation length data for direct signals
  int i=0;
  ifstream sheetup("data/icesheet_attenlength_up.txt");
  if(sheetup.fail())
    {
      cerr << "Failed to open icesheet_attenlength_up.txt" << endl;
      exit(1);
    }
  
  i=0;
  while(sheetup>>d_sheetup[i]>>l_sheetup[i])
    {
      i++;
    }
  sheetup.close();
   
  ifstream shelfup("data/iceshelf_attenlength_up.txt");
  if(shelfup.fail())
    {
      cerr << "Failed to open iceshelf_attenlength_up.txt" << endl;
      exit(1);
    }
  
  i=0;
  while(shelfup>>d_shelfup[i]>>l_shelfup[i])
    {
      i++;
    }
  shelfup.close();
   
  ifstream westlandup("data/westland_attenlength_up.txt");
  if(westlandup.fail())
    {cerr << "Failed to open westland_attenlength_up.txt";
      exit(1);
    }
  i=0;
  while(westlandup>>d_westlandup[i]>>l_westlandup[i])
    {
      i++;
    }
  westlandup.close();
   
  //read in attenuation length for downgoing signals
  ifstream sheetdown("data/icesheet_attenlength_down.txt");
  if(sheetdown.fail())
    {
      cerr << "Failed to open icesheet_attenlength_down.txt" << endl;
      exit(1);
    }
  
  i=0;
  while(sheetdown>>d_sheetdown[i]>>l_sheetdown[i])
    {
      i++;
    }
  sheetdown.close();

  
  ifstream shelfdown("data/iceshelf_attenlength_down.txt");
  if(shelfdown.fail())
    {
      cerr << "Failed to open iceshelf_attenlength_down.txt" << endl;
      exit(1);
    }
  
  i=0;
  while(shelfdown>>d_shelfdown[i]>>l_shelfdown[i])
    {
      i++;
    }
  shelfdown.close();
   
  ifstream westlanddown("data/westland_attenlength_down.txt");
  if(westlanddown.fail())
    {cerr << "Failed to open westland_attenlength_down.txt";
      exit(1);
    }
  i=0;
  while(westlanddown>>d_westlanddown[i]>>l_westlanddown[i])
    {
      i++;
    }
  westlanddown.close();

 }


 //constructor IceModel(int model)
inline Position IceModel::PickBalloonPosition() const {
  Vector temp;
  return temp;

}

inline Vector IceModel::GetSurfaceNormal(const Position &r_out) const {
  Vector n_surf = r_out.Unit();
  if (FLATSURFACE) 
    return n_surf;

  if (ice_model==0) {
    double theta=r_out.Theta();

    int ilon,ilat;
    GetILonILat(r_out,ilon,ilat);
    
    int ilon_previous=ilon-1;
    if (ilon_previous<0)
      ilon_previous=NLON-1;
    
    int ilon_next=ilon+1;
    if (ilon_next==NLON)
      ilon_next=0;
    
    double r=(geoid[ilat]+surfacer[ilon][ilat])*sin(theta);
    
    double slope_phi=(surfacer[ilon_next][ilat]-surfacer[ilon_previous][ilat])/(r*2*phistep);
    
    int ilat_previous=ilat-1;
    if (ilat_previous<0)
      ilat_previous=0;
    
    int ilat_next=ilat+1;
    if (ilat_next==NLAT)
      ilat_next=NLAT-1;
     
    double slope_costheta=(surfacer[ilon][ilat_next]-surfacer[ilon][ilat_previous])/((geoid[ilat]+surfacer[ilon][ilat])*2*thetastep);
    
    // first rotate n_surf according to tilt in costheta and position on continent - rotate around the y axis.
    double angle=atan(slope_costheta);

    n_surf = n_surf.RotateY(angle);

    // now rotate n_surf according to tilt in phi - rotate around the z axis.
    angle=atan(slope_phi);
  
    n_surf = n_surf.RotateZ(angle);
  } //end if(Crust 2.0)
  else if (ice_model==1) {
    double dist_to_check = 7500; //check elevations at this distance north, south, east and west of event
    double lon,lat;
    double lon_prev,lon_next;
    double lat_prev,lat_next;
    lon = r_out.Lon();
    lat = r_out.Lat(); //longitude and latitude of interaction
    double local_surface_elevation = Surface(lon,lat);

    lat_next = lat + dist_to_check * (180 / (local_surface_elevation * PI)); //the latitude 7.5 km south of the interaction
    lat_prev = lat - dist_to_check * (180 / (local_surface_elevation * PI)); //the latitude 7.5 km north of the interaction

    lon_next = lon + dist_to_check * (180 / (sin(lat*RADDEG) * local_surface_elevation * PI)); 
    lon_prev = lon - dist_to_check * (180 / (sin(lat*RADDEG) * local_surface_elevation * PI)); 

    if (lat_next > 90) {
      //cout<<"lat_next is > 90"<<endl;
      lat_next = 90 - (lat_next - 90);  //if we went past the pole, set coordinates for the other side
      lon_next += 180;
      lon_prev += 180;
    } //end if
    //cout<<"lon, lat: "<<lon<<" , "<<lat<<endl;
    //correct any out of range longitudes
    if (lon_next > 360) {
      //cout<<"lon_next > 360\n";
      lon_next -= 360;
    }
    else if (lon_next < 0) {
      //cout<<"lon_next < 0\n";
      lon_next += 360;
    }
    if (lon_prev > 360) {
      //cout<<"lon_prev > 360\n";
      lon_prev -= 360;
    }
    else if (lon_prev < 0) {
      //cout << "lon_prev < 0";
      lon_prev += 360;
    }
   
    double slope_phi=(SurfaceAboveGeoid(lon_next,lat)-SurfaceAboveGeoid(lon_prev,lat))/(2*dist_to_check);

    double slope_costheta=(SurfaceAboveGeoid(lon,lat_next)-SurfaceAboveGeoid(lon,lat_prev))/(2*dist_to_check);
    
    // first rotate n_surf according to tilt in costheta - rotate around the y axis.
    double angle=atan(slope_costheta);

    n_surf = n_surf.RotateY(angle);
   
    // now rotate n_surf according to tilt in phi - rotate around the z axis.
    angle=atan(slope_phi);

    n_surf = n_surf.RotateZ(angle);
  } //end if(BEDMAP)
    
  return n_surf;
    
} //method GetSurfaceNormal

inline Position IceModel::WhereDoesItEnterIce(const Position &posnu,
				       const Vector &nnu,
				       double stepsize) const {
  // now get exit point...
  //   see my geometry notes.
  // parameterize the neutrino trajectory and just see where it
  // crosses the earth radius.

  Position r_enterice;
  double distance=0;
  int left_edge=0;
  Position x = posnu;
  double x2;
  
  Position x_previous = posnu;

  double x_previous2= x_previous * x_previous;
  x2=x_previous2;
  
  double lon = x.Lon(),lat = x.Lat();
  double lon_old = lon,lat_old = lat;
  double local_surface = Surface(lon,lat);
  double rock_previous2= pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
  double surface_previous2=pow(local_surface,2);

  double rock2=rock_previous2;
  double surface2=surface_previous2;


  while (distance<2*local_surface+1000) {

    distance+=stepsize;

    x -= stepsize*nnu;
    x2=x*x;

    lon = x.Lon();
    lat = x.Lat();

    if (lon!=lon_old || lat!=lat_old) {
      local_surface = Surface(lon,lat);

      if (lat>COASTLINE) 
	left_edge=1;
      
      rock2=pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
      surface2=pow(local_surface,2);    

      if (ice_model==0) {
	if ((int)(lat)==COASTLINE && rock_previous2 < x2 && surface2 > x2)
	  left_edge=1;
      } //if (Crust 2.0)
    } //if (neutrino has stepped into new lon/lat bin)

    if ((x_previous2>rock_previous2 && x2<rock2)
	|| (x_previous2<surface_previous2 && x2>surface2)
	|| left_edge) {
       
      r_enterice = x;
      // this gets you out of the loop.
      distance=3*Geoid(lat);
    } //if

    x_previous = x;
    x_previous2 = x2;

    if (lon!=lon_old || lat!=lat_old) {
      rock_previous2 = rock2;
      surface_previous2 = surface2;
      lat_old = lat;
      lon_old = lon;
    } //if

  } //while

  return r_enterice;
}//WhereDoesItEnterIce


inline double IceModel::IceThickness(double lon, double lat) const {
  //This method returns the thickness of the ice in meters at a location specified by a latitude and longitude (in degrees).  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  double ice_thickness=0;
  //cout << "ice_model is " << ice_model << "\n";
  //cout << "icethkarray is " << icethkarray[(int)(lon/2)][(int)(lat/2)]*1000. << "\n";
  if (ice_model==1) {
    int e_coord=0;
    int n_coord=0;
    IceLonLattoEN(lon,lat,e_coord,n_coord);
    if (e_coord <= 1200 && e_coord >= 0 && n_coord <= 1000 && n_coord > 0)
      ice_thickness = ice_thickness_array[e_coord][n_coord]; //if this region has BEDMAP data, use it.
    else
      ice_thickness = icethkarray[(int)(lon/2)][(int)(lat/2)]*1000.; //if the location given is not covered by BEDMAP, use Crust 2.0 data
  } //BEDMAP ice thickness
  else if (ice_model==0) {
    ice_thickness = icethkarray[(int)(lon/2)][(int)(lat/2)]*1000.;
    //cout << "ilon, ilat are " << (int)(lon/2) << " " << (int)(lat/2) << "\n";
  } //Crust 2.0 ice thickness

  return ice_thickness;
} //method IceThickness
inline double IceModel::IceThickness(const Position &pos) const {
  //This method returns the thickness of the ice in meters at a location under a given position vector.  Code by Stephen Hoover.

  return IceThickness(pos.Lon(),pos.Lat());
} //method IceThickness(position)

inline double IceModel::SurfaceAboveGeoid(double lon, double lat) const {
  //This method returns the elevation above the geoid of the surface of the ice (or bare ground, if no ice is present) in meters, at a location specified by a latitude and longitude (in degrees).  In areas covered by water where no ice present, the method returns 0.  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  // lon must be 0 to 360
  double surface=0;

  if (ice_model==1) {
    int e_coord_ice=0;
    int n_coord_ice=0;
    int e_coord_ground=0;
    int n_coord_ground=0;
    IceLonLattoEN(lon,lat,e_coord_ice,n_coord_ice);
    GroundLonLattoEN(lon,lat,e_coord_ground,n_coord_ground);
    if (e_coord_ground <= 1068 && e_coord_ground >= 0 && n_coord_ground <= 869 && n_coord_ground >= 0 && e_coord_ice <= 1200 && e_coord_ice >= 0 && n_coord_ice <= 1000 && n_coord_ice >= 0)
      surface = ground_elevation[e_coord_ground][n_coord_ground] + ice_thickness_array[e_coord_ice][n_coord_ice] + water_depth[e_coord_ice][n_coord_ice];
    else
      surface = surfacer[(int)(lon/2)][(int)(lat/2)]; //If the position requested is outside the bounds of the BEDMAP data, use the Crust 2.0 data, regardless of the ice_model flag.
  } //Elevation of surface above geoid according to BEDMAP
  else if (ice_model==0) {
    surface = surfacer[(int)(lon/2)][(int)(lat/2)];
  } //Elevation of surface above geoid according to Crust 2.0

  return surface;
} //method SurfaceAboveGeoid

inline double IceModel::SurfaceAboveGeoid(const Position &pos) const {
  //This method returns the elevation above the geoid of the surface of the ice (or bare ground, if no ice is present) in meters, at a location specified by a position vector.  Code by Stephen Hoover.

  return SurfaceAboveGeoid(pos.Lon(),pos.Lat());
} //method SurfaceAboveGeoid(position)

inline double IceModel::Surface(double lon,double lat) const {
  return (SurfaceAboveGeoid(lon,lat) + Geoid(lat)); // distance from center of the earth to surface
} //Surface

inline double IceModel::Surface(const Position& pos) const {
  return Surface(pos.Lon(),pos.Lat());
} //Surface

inline double IceModel::WaterDepth(double lon, double lat) const {
  //This method returns the depth of water beneath ice shelves in meters, at a location specified by a latitude and longitude (in degrees).  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  double water_depth_value=0;

  if (ice_model==0) {
      water_depth_value = waterthkarray[(int)(lon/2)][(int)(lat/2)]*1000;
  } //if(Crust 2.0)
  else if (ice_model==1) {
    int e_coord=0;
    int n_coord=0;
    WaterLonLattoEN(lon,lat,e_coord,n_coord);
    if (e_coord <= 1200 && e_coord >= 0 && n_coord <= 1000 && n_coord >= 0)
      water_depth_value = water_depth[e_coord][n_coord];
    else
      water_depth_value = waterthkarray[(int)(lon/2)][(int)(lat/2)]*1000;
  } //else if(BEDMAP)

  return water_depth_value;
} //method WaterDepth(longitude, latitude)
inline double IceModel::WaterDepth(const Position &pos) const {
  //This method returns the depth of water beneath ice shelves in meters, at a location specified by a position vector.  Code by Stephen Hoover.

  return WaterDepth(pos.Lon(),pos.Lat());
} //method WaterDepth(position)

inline int IceModel::IceOnWater(const Position &pos) const
{
  if(IceThickness(pos)>0.&&WaterDepth(pos)>0.)
    return 1;
  else return 0;

}
inline int IceModel::RossIceShelf(const Position &pos) const {
  int ilon,ilat;

  GetILonILat(pos,ilon,ilat);

  if ((ilat==2 && ilon>=5 && ilon<=14) ||
      (ilat==3 && (ilon>=168 || ilon<=14)) ||
      (ilat==4 && (ilon>=168 || ilon<=13)) ||
      (ilat==5 && (ilon>=168 || ilon<=14)))
    return 1;
  else
    return 0;
}//RossIceShelf

inline int IceModel::RossExcept(const Position &pos) const{
  int ilon,ilat;
  GetILonILat(pos,ilon,ilat);
if(ilon<=178&&ilon>=174&&ilat>=4&&ilat<=5)
    return 1;
  else 
    return 0;
}


inline int IceModel::RonneIceShelf(const Position &pos) const {
  int ilon,ilat;

  GetILonILat(pos,ilon,ilat);

  if ((ilat==4 && ilon>=52 && ilon<=74) ||
      (ilat==5 && ilon>=50 && ilon<=71) ||
      (ilat==6 && ilon>=55 && ilon<=64))
    return 1;
  else
    return 0;

}//RonneIceShelf

inline int IceModel::WestLand(const Position &pos) const {
  double lon = pos.Lon() , lat = pos.Lat();

  if((lat>=4&&lat<=26)&&((lon>=0&&lon<=180)||lon>=336))
    return 1;
  else return 0;

}//WestLand


inline int IceModel::OutsideAntarctica(const Position &pos) const {
  return (pos.Lat() >= COASTLINE);
} //OutsideAntarctica(Position)

inline int IceModel::OutsideAntarctica(double lat) const {
  return (lat >= COASTLINE);
} //OutsideAntarctica(double lat)

inline int IceModel::AcceptableRfexit(const Vector &nsurf_rfexit,const Position &rfexit,const Vector &n_exit2rx) const {

  //Make sure there's actually ice where the ray leaves
  if (rfexit.Lat()>COASTLINE || IceThickness(rfexit)<0.0001) {
    cout << "latitude is " << rfexit.Lat() << " compared to COASTLINE at " << COASTLINE << "\n";
    cout << "ice thickness is " << IceThickness(rfexit) << "\n";
    return 0;

  } //if

  if (nsurf_rfexit*n_exit2rx<0) {
    cout << "dot product is " << nsurf_rfexit*n_exit2rx << "\n";
    return 0;
  } //if

  return 1;
} //AcceptableRfexit

inline double IceModel::GetN(double altitude) const {
  // these are Peter's fit parameters
  double a1=0.463251;
  double b1=0.0140157;
  double n=0;

  if (altitude < FIRNDEPTH) 
    n=NICE;
  else if (altitude >= FIRNDEPTH && altitude <=0 && DEPTH_DEPENDENT_N) 
    //    N_DEPTH=NFIRN-(4.6198+13.62*(altitude_int/1000.))*
    //(altitude_int/1000.);   // Besson's equation for n(z)
    n=NFIRN+a1*(1.0-exp(b1*altitude));   // Peter's equation for n(z)
  else if (altitude > 0)
    cout<<"Error!  N requested for position in air!\n";
  else if (!DEPTH_DEPENDENT_N)
    n = NFIRN;

  return n;
} //GetN(altitude)

inline double IceModel::GetN(const Position &pos) const{
  return GetN(pos.Mag() - Surface(pos.Lon(),pos.Lat()));
} //GetN(Position)

inline double IceModel::EffectiveAttenuationLength(const Position &pos,const int &whichray) const {
  double localmaxdepth = IceThickness(pos);
  double depth = Surface(pos) - pos.Mag();
  
  int depth_index=0;
  double attenuation_length=0.0;
//   if (inu<10) {
//     cout << "pos is ";pos.Print();
//     cout << "surface is " << Surface(pos) << "\n";
//   }
  if(WestLand(pos) && !CONSTANTICETHICKNESS) 
    {
      depth_index=int(depth*419.9/localmaxdepth);//use 420 m ice shelf attenuation length data as the standard, squeeze or stretch if localmaxdepth is longer or shorter than 420m.
      if(RossIceShelf(pos) || RonneIceShelf(pos)) 
	{	  
	  if(whichray==0)
	    attenuation_length=l_shelfup[depth_index];
	  else if(whichray==1)
	    attenuation_length=l_shelfdown[depth_index];
	  else
	    cerr << " wrong attenuation length " <<endl;
	  
	  //for sanity check
	  if((depth_index+0.5)!=d_shelfup[depth_index])
	    {
	      cerr << "the index of the array l_iceshelfup is wrong!" << endl;
	      exit(1);
	    }
	}
      else //in ice sheet of westland
	{
	  if(whichray==0)
	    attenuation_length=l_westlandup[depth_index]; 
	  else if(whichray==1)
	    attenuation_length=l_westlanddown[depth_index];
	  else
	    cerr << " wrong attenuation length " <<endl;
      	}
       
      if(mooreBayFlag)//if use Moore's Bay measured data for the west land
	attenuation_length*=1.717557; //about 450 m (field attenuation length) for one whole way when assuming -3dB for the power loss at the bottom
    }
  else //in east antarctica or constant ice thickness
     { 
//        if (inu<10) {
//        cout << "localmaxdepth is " << localmaxdepth << "\n";
//        cout << "depth is " << depth << "\n";
//        }
       depth_index =int(depth*(2809.9/localmaxdepth));
       //if (inu<10)
	 //       cout << "depth_index is " << depth_index << "\n";

       if(whichray==0)
	 attenuation_length =l_sheetup[depth_index];
       else if(whichray==1)
	 attenuation_length =l_sheetdown[depth_index];
       else
	 cerr << " wrong attenuation length " <<endl;
     } //else

  return attenuation_length;
} //EffectiveAttenuationLengthUp

inline double IceModel::Area(double latitude) const {
  //Returns the area of one square of the BEDMAP data at a given latitude. 
  double lat_rad = (90 - latitude) * RADDEG;

  return (pow(cellSize* ((1 + sin(71*RADDEG)) / (1 + sin(lat_rad))),2));
} //method Area

inline void IceModel::LonLattoEN(double lon, double lat, double xLowerLeft, double yLowerLeft, int& e_coord, int& n_coord) const {
  //takes as input a latitude and longitude (in degrees) and converts to indicies for BEDMAP matricies. Needs a location for the corner of the matrix, as not all the BEDMAP files cover the same area.  Code by Stephen Hoover.

  double easting=0;
  double northing=0;

  double lon_rad = (lon - 180) * RADDEG; //convert to radians, and shift origin to conventional spot
  double lat_rad = (90 - lat) * RADDEG;

  bedmap_R = scale_factor*bedmap_c_0 * pow(( (1 + eccentricity*sin(lat_rad)) / (1 - eccentricity*sin(lat_rad)) ),eccentricity/2) * tan((PI/4) - lat_rad/2);

  easting = bedmap_R * sin(lon_rad);
  northing = bedmap_R * cos(lon_rad);

  //  cout << "bedmap_R is " << bedmap_R << "\n";
  //cout << "easting, northing are " << easting << " " << northing << "\n";

  e_coord = (int)((easting - xLowerLeft) / cellSize);
  n_coord = (int)((-1*northing - yLowerLeft) / cellSize);

  return;
} //method LonLattoEN

inline void IceModel::IceLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
  //Converts a latitude and longitude (in degrees) to indicies for BEDMAP ice thickness data.  Code by Stephen Hoover.
  LonLattoEN(lon, lat, xLowerLeft_ice, yLowerLeft_ice, e_coord, n_coord);
}//IceLonLattoEN
inline void IceModel::GroundLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
  //Converts a latitude and longitude (in degrees) to indicies for BEDMAP ground elevation data.  Code by Stephen Hoover.
  LonLattoEN(lon, lat, xLowerLeft_ground, yLowerLeft_ground, e_coord, n_coord);
}//GroundLonLattoEN
inline void IceModel::WaterLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
  //Converts a latitude and longitude (in degrees) to indicies for BEDMAP water depth data.  Code by Stephen Hoover.
  LonLattoEN(lon, lat, xLowerLeft_water, yLowerLeft_water, e_coord, n_coord);
}//WaterLonLattoEN

inline void IceModel::ENtoLonLat(int e_coord, int n_coord, double xLowerLeft, double yLowerLeft, double& lon, double& lat) const {
  //Takes as input the indicies from a BEDMAP data set, and turns them into latitude and longitude coordinates.  Information on which data set (surface data, ice depth, water depth) is necessary, in the form of coordinates of a corner of the map.  Code by Stephen Hoover.

  double isometric_lat=0;
  double easting = xLowerLeft+(cellSize*(e_coord+0.5)); //Add offset of 0.5 to get coordinates of middle of cell instead of edges.
  double northing = -1*(yLowerLeft+(cellSize*(n_coord+0.5)));

  //  cout << "easting, northing are " << easting << " " << northing << "\n";

  //first set longitude

  if (northing!=0)
    lon = atan(easting/northing);
  else
    lon = 90*RADDEG;

  // this puts lon between -pi and pi
  if (easting > 0 && lon < 0) //adjust sign of longitude
    lon += PI;
  else if (easting < 0 && lon > 0)
    lon -= PI;
  else if (easting == 0 && northing < 0)
    lon += PI;

  //  now find latitude

  if (easting != 0)
    bedmap_R = fabs(easting/sin(lon));
  else if (easting == 0 && northing != 0)
    bedmap_R = fabs(northing);
  else {
    lat = 0; //at the pole, set lat=0 degrees
    lon = lon*DEGRAD; // now put lon between 180 and 180 (only at pol)
    return;
  } //else

  isometric_lat = (PI/2) - 2*atan(bedmap_R/(scale_factor*bedmap_c_0));

  lat = isometric_lat + bedmap_a_bar*sin(2*isometric_lat) + bedmap_b_bar*sin(4*isometric_lat) + bedmap_c_bar*sin(6*isometric_lat) + bedmap_d_bar*sin(8*isometric_lat);

  lon = lon * DEGRAD + 180;  //convert to degrees, shift 0 to line up with bin 0 of Crust 2.0
  lat = 90 - lat*DEGRAD; //convert to degrees, with 0 degrees at the south pole

  //  if (lon>160 && lon<165)
  //cout << "e_coord, n_coord, easting, northing, lon are " << e_coord << " " << n_coord << " " << easting << " " << northing << " " << lon << "\n";
  return;
  
} //method ENtoLonLat

inline void IceModel::IceENtoLonLat(int e, int n, double& lon, double& lat) const {
  //Converts indicies of the BEDMAP ice thickness matrix into longitude and latitude.  Code by Stephen Hoover.
  // cout << "I'm inside IceENtoLonLat.\n";
  ENtoLonLat(e,n,xLowerLeft_ice,yLowerLeft_ice,lon,lat);
}//IceENtoLonLat
inline void IceModel::GroundENtoLonLat(int e, int n, double& lon, double& lat) const {
  //Converts indicies of the BEDMAP ground elevation matrix into longitude and latitude.  Code by Stephen Hoover.
  ENtoLonLat(e,n,xLowerLeft_ground,yLowerLeft_ground,lon,lat);
}//GroundENtoLonLat
inline void IceModel::WaterENtoLonLat(int e, int n, double& lon, double& lat) const {
  //Converts indicies of the BEDMAP water depth matrix into longitude and latitude.  Code by Stephen Hoover.
  ENtoLonLat(e,n,xLowerLeft_water,yLowerLeft_water,lon,lat);
}//WaterENtoLonLat



void IceModel::ReadIceThickness() {
  //Reads the BEDMAP ice thickness data.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover

  ifstream IceThicknessFile("data/icethic.asc");
  if(!IceThicknessFile) {
    cerr << "Couldn't open: data/icethic.asc" << endl;
    exit(1);
  }

  cout<<"Reading in BEDMAP data on ice thickness.\n";

  string tempBuf1;
  string tempBuf2;
  string tempBuf3;
  string tempBuf4;
  string tempBuf5;
  string tempBuf6;
  int temp1,temp2,temp3,temp4,temp5,temp6;
  
  IceThicknessFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
		   >> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
		   >> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
  
  if(tempBuf1 == string("ncols")) {
    nCols_ice=temp1;
  }
  if(tempBuf2 == string("nrows")) {
    nRows_ice=temp2;
  }
  if(tempBuf3 == string("xllcorner")) {
    xLowerLeft_ice=temp3;
  }
  if(tempBuf4 == string("yllcorner")) {
    yLowerLeft_ice=temp4;
  }
  if(tempBuf5 == string("cellsize")) {
    cellSize=temp5;
  }
  if(tempBuf6 == string("NODATA_value")) {
    NODATA=temp6;
  }
  //cout<<"nCols_ice, nRows_ice "<<nCols_ice<<" , "<<nRows_ice<<endl;
  //cout<<"xLL_ice, yLL_ice, cellsize "<<xLowerLeft_ice<<" , "<<yLowerLeft_ice<<" , "<<cellSize<<endl<<endl;
  
  double theValue;
  for(int rowNum=0;rowNum<nRows_ice;rowNum++) {
    for(int colNum=0;colNum<nCols_ice;colNum++) {
      IceThicknessFile >> theValue;
      if(theValue==NODATA)
	theValue=0; //Set ice depth to 0 where we have no data.
      ice_thickness_array[colNum][rowNum] = double(theValue); //This stores data as ice_thickness_array[easting][northing]
    }//for
  }//for
  
  IceThicknessFile.close();
  return;
} //method ReadIceThickness

void IceModel::ReadGroundBed() {
  //Reads the BEDMAP data on the elevation of the ground beneath the ice.  If there is water beneath the ice, the ground elevation is given the value 0.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
  ifstream GroundBedFile("data/groundbed.asc");
  if(!GroundBedFile) {
    cerr << "Couldn't open: data/groundbed.asc" << endl;
    exit(1);
  }

  cout<<"Reading in BEDMAP data on elevation of ground.\n";

  string tempBuf1;
  string tempBuf2;
  string tempBuf3;
  string tempBuf4;
  string tempBuf5;
  string tempBuf6;
  int temp1,temp2,temp3,temp4,temp5,temp6;
  
  GroundBedFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
		>> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
		>> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
  
  if(tempBuf1 == string("ncols")) {
    nCols_ground=temp1;
  }
  if(tempBuf2 == string("nrows")) {
    nRows_ground=temp2;
  }
  if(tempBuf3 == string("xllcorner")) {
    xLowerLeft_ground=temp3;
  }
  if(tempBuf4 == string("yllcorner")) {
    yLowerLeft_ground=temp4;
  }
  if(tempBuf5 == string("cellsize")) {
    cellSize=temp5;
  }
  if(tempBuf6 == string("NODATA_value")) {
    NODATA=temp6;
  }

  //cout<<"nCols_ground, nRows_ground "<<nCols_ground<<" , "<<nRows_ground<<endl;
  //cout<<"xLL_ground, yLL_ground, cellsize "<<xLowerLeft_ground<<" , "<<yLowerLeft_ground<<" , "<<cellSize<<endl<<endl;
  
  double theValue;
  for(int rowNum=0;rowNum<nRows_ground;rowNum++) {
    for(int colNum=0;colNum<nCols_ground;colNum++) {
      GroundBedFile >> theValue;
      
      if(theValue==NODATA)
	theValue=0; //Set elevation to 0 where we have no data.
      ground_elevation[colNum][rowNum] = double(theValue);
      //if (theValue != -96 && theValue != 0)
      //cout<<"ground_elevation: "<<theValue<<endl;
    }//for
  }//for
  
  GroundBedFile.close();
  return;
} //method ReadGroundBed

void IceModel::ReadWaterDepth() {
  //Reads BEDMAP data on the depth of water beneath the ice.  Where no water is present, the value 0 is entered.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
  ifstream WaterDepthFile("data/water.asc");
  if(!WaterDepthFile) {
    cerr << "Couldn't open: data/water.asc" << endl;
    exit(1);
  }

  cout<<"Reading in BEDMAP data on water depth.\n";

  string tempBuf1;
  string tempBuf2;
  string tempBuf3;
  string tempBuf4;
  string tempBuf5;
  string tempBuf6;
  int temp1,temp2,temp3,temp4,temp5,temp6;
  
  WaterDepthFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
		 >> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
		 >> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
  
  if(tempBuf1 == string("ncols")) {
    nCols_water=temp1;
  }
  if(tempBuf2 == string("nrows")) {
    nRows_water=temp2;
  }
  if(tempBuf3 == string("xllcorner")) {
    xLowerLeft_water=temp3;
  }
  if(tempBuf4 == string("yllcorner")) {
    yLowerLeft_water=temp4;
  }
  if(tempBuf5 == string("cellsize")) {
    cellSize=temp5;
  }
  if(tempBuf6 == string("NODATA_value")) {
    NODATA=temp6;
  }

  //cout<<"nCols_water, nRows_water "<<nCols_water<<" , "<<nRows_water<<endl;
  //cout<<"xLL_water, yLL_water, cellsize "<<xLowerLeft_water<<" , "<<yLowerLeft_water<<" , "<<cellSize<<endl<<endl;
  
  double theValue;
  for(int rowNum=0;rowNum<nRows_water;rowNum++) {
    for(int colNum=0;colNum<nCols_water;colNum++) {
      WaterDepthFile >> theValue;
      
      if(theValue==NODATA)
	theValue=0; //Set depth to 0 where we have no data.
      water_depth[colNum][rowNum] = double(theValue);
    }//for
  }//for
  
  WaterDepthFile.close();
  return;
} //method ReadWaterDepth

//void IceModel::FillArraysforTree(double icethck[1200][1000],double elev[1068][869],double lon_ground[1068][869],double lat_ground[1068][869],double lon_ice[1200][1000],double lat_ice[1200][1000],double h20_depth[1200][1000],double lon_water[1200][1000],double lat_water[1200][1000]) {
// void IceModel::FillArraysforTree(double lon_ground[1068][869],double lat_ground[1068][869],double lon_ice[1200][1000],double lat_ice[1200][1000],double lon_water[1200][1000],double lat_water[1200][1000]) {
 
  //  for (int rowNum=0;rowNum<nRows_ice;rowNum++) {
  //for (int colNum=0;colNum<nCols_ice;colNum++) {
//   for (int i=0;i<nRows_ice;i++) {
//     for (int j=0;j<nCols_ice;j++) {
//       //     this->IceENtoLonLat(colNum,rowNum,lon_ice[colNum][rowNum],lat_ice[colNum][rowNum]); //Recall that the e / n coordinates in horizon were picked from the ground bed array.
//       double test1,test2;
//       cout << "rowNum, colNum are " << i << " " << j << "\n";
//       //      cout << "lon_ice is " << rowNum << " " << colNum << " " << lon_ice[rowNum][colNum] << "\n";
//       //this->IceENtoLonLat(i,j,test1,test2); //Recall that the e / n coordinates in horizon were picked from the ground bed array.
//       //icethck[colNum][rowNum]=IceThickness(lon_ice[colNum][rowNum],lat_ice[colNum][rowNum]);
     
//       //      WaterENtoLonLat(colNum,rowNum,lon_water[colNum][rowNum],lat_water[colNum][rowNum]);
//       //h20_depth[colNum][rowNum]=WaterDepth(lon_water[colNum][rowNum],lat_water[colNum][rowNum]);

//     }
//   }h
//   for (int n_coord=0;n_coord<nRows_ground;n_coord++) {
//     for (int e_coord=0;e_coord<nCols_ground;e_coord++) {
//       GroundENtoLonLat(e_coord,n_coord,lon_ground[e_coord][n_coord],lat_ground[e_coord][n_coord]);
//       //     elev[e_coord][n_coord] = SurfaceAboveGeoid(lon_ice[e_coord][n_coord],lat_ice[e_coord][n_coord]);

//     }
//   }





//}
