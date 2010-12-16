#include "Trigger.h"
#include "Efficiencies.h"
#include "Tools.h"
#include <vector>
ClassImp(Trigger);
Trigger::Trigger() {

}

    


//Trigger::~Trigger() {

//}
void Trigger::simulateTrigger(int nRx,vector <vector <double> > waveforms,Efficiencies *efficiencies) {
  l1Hits=l1Trigger(nRx,waveforms); // loop over waveforms and decide if they pass the l1 triggers
  //l1Hits=1; // loop over waveforms and decide if they pass the l1 triggers
  efficiencies->incrementL1Counter(l1Hits);

}
std::vector<int> Trigger::l1Trigger(int nRx, std::vector< vector<double> > waveforms) {
  //std::vector<int> l1Hits;
  for (int i=0;i<nRx;i++) {
    if (Tools::getMaxMagnitude(waveforms[i])>threshold)
      l1Hits.push_back(1);
    else
      l1Hits.push_back(0);
  }
  


  return l1Hits;
}
std::vector<int> Trigger::l2Trigger(std::vector<int> l1Hits) {
  std::vector<int> l2Hits;
  for (int i=0;i<(int)l1Hits.size();i++) {
    l2Hits.push_back(l1Hits[i]);
  }

  return l2Hits;

}
std::vector<int> Trigger::l3Trigger(std::vector<int> l2Hits) {
  std::vector<int> l3Hits;
  for (int i=0;i<(int)l2Hits.size();i++) {
    l3Hits.push_back(l2Hits[i]);
  }


  return l3Hits;
}
void Trigger::resetTrigger() {

}
