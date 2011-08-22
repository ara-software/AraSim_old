#include "TSpline.h"
#include <string>
#include "TRandom3.h"

using namespace std;

class Spectra {
private:
  double maxflux;      


public:  
  static const int NSPECTRA_MAX=300;  
  double energy[12]; // energies that correspond to the fluxes in the previous array  
  double EdNdEdAdt[12]; //flux of incident neutrinos vs. energy E*dN/dE/dA/dt
  double E2dNdEdAdt[12]; //flux of incident neutrinos vs. energy E^2*dN/dE/dA/dt
  TSpline5 *spectrum;
  Spectra(int EXPONENT); // constructor  
  double GetNuEnergy(); // get the neutrino energy  
  void GetFlux(string filename);

  TGraph *gspectrum[NSPECTRA_MAX];      
  

}; //class Spectra
