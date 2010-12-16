#include "Tools.h"
#include <iostream>
#include <cmath>

 double Tools::getMaxMagnitude(vector<double> v) {
  double mag=0.;
  for (int i=0;i<(int)v.size();i++) {
    if (v[i]>mag)
      mag=v[i];

  }
  return mag;

}
 double Tools::distance(vector<double> vec1,vector<double> vec2) {
  int size;
  if (vec1.size()!=vec2.size()) {
    cout << "Warning!!!! Finding distance between two vectors of different sizes!!\n";
    size=Tools::iMin((int)vec1.size(),(int)vec2.size());
  }
  else
    size=vec1.size();
  double length=0.;
  for (int i=0;i<size;i++) {
    length+=(vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
  }
  return sqrt(length);
}

 double Tools::dSquare(double *p) {
  return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
} //dSquare

 void Tools::Zero(int *anarray,int n) {
  for (int i=0;i<n;i++) {
    anarray[i]=0;
  } //for
} //Zero (int*,int)

 void Tools::Zero(double *anarray,int n) {
  for (int i=0;i<n;i++) {
    anarray[i]=0.;
  } //for
} //Zero (int*,int)
void Tools::ShiftLeft(double *x,const int n,int ishift) {

  double x_temp[n];
  // shift the x array to the left by ishift bins and fill the gap with zeroes
  for (int i=0;i<n;i++) {
    x_temp[i]=x[i];
  }
  for (int i=0;i<n-ishift;i++) {
    x[i]=x_temp[i+ishift];
  }
  for (int i=n-ishift;i<n;i++) {
    x[i]=0.;
  }

}
void Tools::ShiftRight(double *x,const int n,int ishift) {

  double x_temp[n];
  // shift the x array to the right by ishift bins and fill the gap with zeroes
  for (int i=0;i<n;i++) {
    x_temp[i]=x[i];
  }
  for (int i=ishift;i<n;i++) {
    x[i]=x_temp[i-ishift];
  }
  for (int i=0;i<ishift;i++) {
    x[i]=0.;
  }

}

void Tools::realft(double *data, const int isign, int nsize){
  int i, i1, i2, i3, i4;
  double c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;
  theta=3.141592653589793238/(nsize>>1);
  if (isign == 1) {
                c2 = -0.5;
                four1(data,1,nsize);
        } else {
                c2=0.5;
                theta = -theta;
        }
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0+wpr;
        wi=wpi;
        for (i=1;i<(nsize>>2);i++) {
                i2=1+(i1=i+i);
                i4=1+(i3=nsize-i1);
                h1r=c1*(data[i1]+data[i3]);
                h1i=c1*(data[i2]-data[i4]);
                h2r= -c2*(data[i2]+data[i4]);
                h2i=c2*(data[i1]-data[i3]);
                data[i1]=h1r+wr*h2r-wi*h2i;
                data[i2]=h1i+wr*h2i+wi*h2r;
                data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4]= -h1i+wr*h2i+wi*h2r;
                wr=(wtemp=wr)*wpr-wi*wpi+wr;
                wi=wi*wpr+wtemp*wpi+wi;
        }
        if (isign == 1) {
                data[0] = (h1r=data[0])+data[1];
                data[1] = h1r-data[1];
        } else {
                data[0]=c1*((h1r=data[0])+data[1]);
                data[1]=c1*(h1r-data[1]);
                four1(data,-1,nsize);
        }
}

void Tools::four1(double *data, const int isign,int nsize) {
	int n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

	int nn=nsize/2;
	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j-1],data[i-1]);
			SWAP(data[j],data[i]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

 double Tools::dMinNotZero(const double *x,int n) {
  double min=dMax(x,n);
  if (min==0)
    cout << "max is 0.\n";
  for (int k=1;k<n;k++) {
    if (x[k]<min && x[k]!=0)
      min=x[k];
  }
  return min;
} //dMinNotZero(double*, int)

 double Tools::dMin(const double *x,int n) {
  double min=x[0];
  for (int k=1;k<n;k++) {
    if (x[k]<min)
      min=x[k];
  }
  return min;
} 
 int Tools::iMin(int x,int y) {

  int min;
  if (x<y)
    min=x;
  else
    min=y;
  
  return min;

}


 double Tools::dMin(double x,double y) {
  double min=1.E22;
  if (x<y)
    min=x;
  else
    min=y;
  
  return min;
} //dMin(double,double)


 double Tools::dMax(const double *x,int n) {
  
  double max=x[0];
  for (int k=1;k<n;k++) {
    if (x[k]>max)
      max=x[k];
  }
  return max;
} //dMax(double*, int)


 double Tools::dvMax(const vector<double> x) {
  
  double max=x[0];
  for (int k=1;k<(int)x.size();k++) {
    if (x[k]>max)
      max=x[k];
  }
  return max;
} //dMax(double*, int)
 double Tools::dsMax(TSpline5 *sp) {
  vector<double> y;
  double maxn;
 double blah1,blah2;
  for (int i=0;i<sp->GetNp();i++) {
    sp->GetKnot(i,blah1,blah2);
    y.push_back(blah2);
  }
  maxn=Tools::dvMax(y);
  return maxn;
}

 double Tools::dMax(double a,double b) {
  if (a>b)
    return a;
  else if (a<b)
    return b;
  else if (a==b)
    return a;
  return 0;
} //dMax(double,double
 int Tools::Getifreq(double freq,double freq_low,double freq_high,int n) {

  if (freq>=freq_high)
    return -1;
  if (freq<freq_low)
    return -1;

  return (int)((freq-freq_low)/(freq_high-freq_low)*(double)n);
} //Getifreq
void Tools::InterpolateComplex(double *array, const int n) {
  // to get rid of the zero bins
  double previous_nonzero=0.;
  double next_nonzero=0.;
  double check;
  int ifirstnonzero=0;
  int ilastnonzero=0;
  int k;
  int m=0;
  int count_nonzero=0;

  // find the first nonzero even element
  while (array[2*m]==0) {
    m++;
  }
  ifirstnonzero=m;
  
  
  // count the nonzero elements
  for (int i=0;i<n;i++) {
    if (array[2*i]!=0)
      count_nonzero++;
  }
  if (count_nonzero!=0) {

    // loop through the elements of the array and replace the zeros with interpolated values
  for (int i=ifirstnonzero;i<n;i++) {
    
    if (array[2*i]!=0.) {
      // set the lower nonzero value that we are interpolating from 
      previous_nonzero=array[2*i];
    
    }
    else {
      check=0.;
      k=i;
      while (check==0. && k<n) {
	check=array[2*k];
	k++;
      }
      if (k<n) {
      next_nonzero=check;
    
      for (int j=i;j<k;j++) {
	array[2*j]=previous_nonzero+(next_nonzero-previous_nonzero)*(double)(j-(i-1))/(double)(k-(i-1));
	array[2*j+1]=array[2*j];
	
      }
      i=k-1;
      
      previous_nonzero=next_nonzero;
      }
      else {
	ilastnonzero=i-1;
	i=n;
      }
    } // end if array=0
  } // end loop over i
  
  //if (inu==49416)
  //cout << "inu, count_nonzero, diff are " << inu << " " << count_nonzero << " " << ilastnonzero << " " << ifirstnonzero << "\n";
  //cout << "factor is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
  for (int j=0;j<n;j++) {

    array[2*j]*=sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
    array[2*j+1]*=sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
  }
  }

}
void Tools::NormalTimeOrdering(const int n,double *volts) {
  double volts_temp[n];
  for (int i=0;i<n/2;i++) {
    volts_temp[i]=volts[i+n/2];
    volts_temp[i+n/2]=volts[i];
  }
  for (int i=0;i<n;i++) {
    volts[i]=volts_temp[i];
  }

}
