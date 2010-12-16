////////////////////////////////////////////////////////////////////////////////////////////////
//class Tools:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef TOOLS_H
#define TOOLS_H

#include "TSpline.h"
#include "TRandom3.h"
#include <vector>
#include <string>

using namespace std;

class Tools {


public:
  static double distance(vector<double> vec1,vector<double> vec2);
  static double dSquare(double*);
  static double dMax(double,double);
  static double dMax(const double*,int);
  static double dvMax(const vector<double>);
  static double dsMax(TSpline5 *sp);
  static double dMin(const double*,int);
  static double dMinNotZero(const double*,int);
  static double dMin(double,double);
  static double getMaxMagnitude(vector<double> v);
  static int iMin(int,int);
  static int Getifreq(double freq,double freq_low,double freq_high,int n);
  static void InterpolateComplex(double *array, const int n);

  static void four1(double *data, const int isign,int nsize);
  static void realft(double *data, const int isign, int nsize);

  static void SWAP(double &a, double &b) // swaps two numbers
  {double dum=a; a=b; b=dum;}
  static void NormalTimeOrdering(const int n,double *volts);
  static void ShiftLeft(double *x,const int n,int ishift);  
  static void ShiftRight(double *x,const int n,int ishift);  
  static void Zero(double *anarray,int n);
  static void Zero(int *anarray,int n);
protected:

};


#endif //TOOLS_H
