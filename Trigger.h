////////////////////////////////////////////////////////////////////////////////////////////////
//class Trigger:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef TRIGGER_H
#define TRIGGER_H

#include <TObject.h>
#include <vector>
using namespace std;

class Efficiencies;
class Trigger : public TObject{


 private:

  vector<int> l1Hits;
  vector<int> l2Hits;
  vector<int> l3Hits;
  int nSamples;
  vector< vector<double> > waveforms; // these are the waveforms in the trigger path, different from the waveforms in the signal path used for analysis
  double threshold;

  vector<int> l1Trigger(int nRx,vector< vector<double> > waveforms);
  vector<int> l2Trigger(vector<int> l1Hits);
  vector<int> l3Trigger(vector<int> l2Hits);
  

 public:
  Trigger();
 Trigger(int nRx) : l1Hits(nRx), l2Hits(nRx), l3Hits(nRx), nSamples(128), waveforms(nRx,vector<double>(nSamples)), threshold(3.) {} ;

  //  ~Trigger();
  void resetTrigger();
  void simulateTrigger(int nRx,vector <vector <double> > waveforms,Efficiencies *efficiencies);

  
 protected:



  ClassDef(Trigger,1);

};

#endif //TRIGGER_H
