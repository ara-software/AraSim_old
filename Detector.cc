#include "Detector.h"
#include "Tools.h"
#include "Event.h"

ClassImp(Detector);

void Detector::simulateDetector(Event *event,Efficiencies *efficiencies) { // simulates the detector response, including the waveforms for analysis
  double r;
  for (int i=0;i<nRx;i++) {
    
    r=Tools::distance(event->posnu,rxpos[i]);
  }
  trigger->simulateTrigger(nRx,waveforms,efficiencies);
}
void Detector::resetDetector() {
  //trigger->resetTrigger();
  for (int i=0;i<nRx;i++) {
    for (int j=0;j<nSamples;j++) {
      waveforms[i][j]=0.;
    }
  }

}
int Detector::getnRx() {
  return nRx;
}
Trigger* Detector::getTrigger() {
  return trigger;
}
Detector::~Detector() {

}
