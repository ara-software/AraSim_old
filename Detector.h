////////////////////////////////////////////////////////////////////////////////////////////////
//class Detector:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef DETECTOR_H
#define DETECTOR_H

#include <TObject.h>
#include <vector>
#include "Trigger.h"
using namespace std;

//#include "Trigger.h"
class Event;
class Efficiencies;

class Detector : public TObject {

  //const std::auto_ptr<Trigger> trigger;

  int nRx; // number of antennas
  Trigger *trigger;
  int nSamples; // number of samples in waveforms for analysis
  vector< vector<double> > rxpos; // position of antennas 
  vector< vector<double> > waveforms;// waveforms for analysis 
  vector<int> rxtype; // antenna types (0=dipole, 1=slot) 
  
  




 public:
 Detector(): nRx(100), trigger(new Trigger(nRx)),  nSamples(128), rxpos(nRx,vector<double>(3)), waveforms(nRx,vector<double>(nSamples)), rxtype(nRx)  {} 
  void simulateDetector(Event *event,Efficiencies *efficiencies); // simulates the detector response, including the waveforms at each antenna for analysis
  void resetDetector();
  int getnRx();
  Trigger* getTrigger();
  ~Detector();

 protected:
  

			   
    ClassDef(Detector,1);



};

#endif //DETECTOR_H
