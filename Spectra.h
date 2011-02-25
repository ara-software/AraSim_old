////////////////////////////////////////////////////////////////////////////////////////////////
//class Spectra:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SPECTRA_H
#define SPECTRA_H

#include <string>


using namespace std;

class Spectra {

private:
  double maxflux;

public:
  double energy[12]; // energies that correspond to the fluxes in the previous array
  double EdNdEdAdt[12]; //flux of incident neutrinos vs. energy
  Spectra(int EXPONENT); // constructor
  double GetNuEnergy(); // get the neutrino energy
  void GetFlux(string filename);
}; //class Spectra

#endif //SPECTRA_H
