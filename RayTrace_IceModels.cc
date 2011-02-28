#include "RayTrace_IceModels.h"

#include <limits>
#include <stdexcept>

exponentialRefractiveIndex::exponentialRefractiveIndex(double n_surface, double n_deep, double transition):
A(n_deep),B(n_surface-n_deep),C(transition){}

double exponentialRefractiveIndex::indexOfRefraction(double z) const{
	if(z<0.0)
		return(1.0);
	return(A+B*exp(C*z));
}

double exponentialRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z<0.0)
		return(0.0);
	return(B*C*exp(C*z));
}

void exponentialRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z<0.0){
		n=1.0;
		dndz=0.0;
	}
	else{
		n=B*exp(C*z);
		dndz=C*n;
		n+=A;
	}
}

RayTrace::indexOfRefractionModel::RayEstimate exponentialRefractiveIndex::estimateRayAngle(double sourceDepth, double receiverDepth, double distance) const{
	if(B==0.0){
		//in this degenerate case, n(z)==A for all z
		//which causes problems when we divide by (n-A) and (n0-A) below
		//however, this case is really simple, so we can handle it directly
		double theta=atan(distance/(receiverDepth-sourceDepth));
		if(theta<0.0)
			theta+=RayTrace::pi;
		return(RayEstimate(SOLUTION,theta));
	}
	
	bool swap=(sourceDepth>receiverDepth);
	if(swap)
		std::swap(sourceDepth,receiverDepth);
	double n0=A+B*exp(C*sourceDepth);
	double n=A+B*exp(C*receiverDepth);
	double sDiff=C*distance;
	
	double eps=1e-4;
	double t1=eps, t2=(n<n0?asin(n/n0):RayTrace::pi/2);
	double a,b,c,s;
	
	double f;
	s=sin(t1);
	a=s*s*n0*n0;
	b=sqrt(A*A-a);
	c=A/b;
	f=log((((sqrt(n*n-a)+b)/(n-A))+c)/(((sqrt(n0*n0-a)+b)/(n0-A))+c))+((b*sDiff)/(s*n0));
	
	double fmid;
	s=sin(t2);
	a=s*s*n0*n0;
	b=sqrt(A*A-a);
	c=A/b;
	fmid=log((((sqrt(n*n-a)+b)/(n-A))+c)/(((sqrt(n0*n0-a)+b)/(n0-A))+c))+((b*sDiff)/(s*n0));
	
	if(f*fmid>=0.0 || std::isnan(fmid)){
		//std::cout << "Did not bracket a root" << std::endl;
		if(swap)
			return(RayEstimate(LOWER_LIMIT,RayTrace::pi-asin((n0/n)*sin(t2))));
		return(RayEstimate(LOWER_LIMIT,t2));
	}
	
	double dt,tmid;
	double rtb=(f<0.0?(dt=t2-t1,t1):(dt=t1-t2,t2));
	const unsigned int maxIter=40;
	unsigned int i;
	for(i=0; i<maxIter; i++){
		dt*=0.5;
		tmid=rtb+dt;
		s=sin(tmid);
		a=s*s*n0*n0;
		b=sqrt(A*A-a);
		c=A/b;
		fmid=log((((sqrt(n*n-a)+b)/(n-A))+c)/(((sqrt(n0*n0-a)+b)/(n0-A))+c))+((b*sDiff)/(s*n0));
		if(fmid<=0.0)
			rtb=tmid;
		if(std::abs(fmid)<1.0e-4)
			break;
	}
	if(i==maxIter)
		return(RayEstimate());
	if(swap)
		return(RayEstimate(SOLUTION,RayTrace::pi-asin((n0/n)*sin(tmid))));
	return(RayEstimate(SOLUTION,tmid));
}


inverseExponentialRefractiveIndex::inverseExponentialRefractiveIndex(double n_surface, double n_deep, double transition):
A(2*n_surface-n_deep),B(2.*(n_deep-n_surface)),C(transition){}

double inverseExponentialRefractiveIndex::indexOfRefraction(double z) const{
	if(z<0.0)
		return(1.0);
	return(A+B/(1.+exp(C*z)));
}
double inverseExponentialRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z<0.0)
		return(0.0);
	return(-B*C*exp(C*z)/((1.+exp(C*z))*(1.+exp(C*z))));
}
void inverseExponentialRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z<0.0){
		n=1.0;
		dndz=0.0;
	}
	else{
		double e=exp(C*z);
		n=A+B/(1.+e);
		dndz=-B*C*e/((1.+e)*(1.+e));
	}
}


simpleExponentialRefractiveIndex::simpleExponentialRefractiveIndex(double a, double b):
A(a),B(b){}

double simpleExponentialRefractiveIndex::indexOfRefraction(double z) const{
	if(z<0.0)
		return(1.0);
	return(A*exp(B*z));
}
double simpleExponentialRefractiveIndex::indexOfRefractionDerivative(double z) const{
	return(A*B*exp(B*z));
}
void simpleExponentialRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z<0.0){
		n=1.0;
		dndz=0.0;
	}
	else{
		n=A*exp(B*z);
		dndz=B*n;
	}
}

quadraticRefractiveIndex::quadraticRefractiveIndex(double a, double b, double c):
A(a),B(b),C(c),maxPoint(B/(2.*C)),maxVal(A+(B*B)/(4.*C)){}

double quadraticRefractiveIndex::indexOfRefraction(double z) const{
	if(z<0.0)
		return(1.0);
	if(z>maxPoint)
		return(maxVal);
	return(A+(B-C*z)*z);
}
double quadraticRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z<0.0 || z>maxPoint)
		return(0.0);
	return(B-2.*C*z);
}
void quadraticRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z<0.0){
		n=1.0;
		dndz=0.0;
	}
	else if(z>maxPoint){
		n=maxVal;
		dndz=0.0;
	}
	else{
		n=(A+(B-C*z)*z);
		dndz=B-2*C*z;
	}
}


double todorDensity(double z){
	if (z<0)
		throw std::domain_error("todorDensity is defined only for positive depths");
	float rho=1.;
	static const float par[7]={0.9283,2.375,0.0249,-3.095,0.0386,1.354,0.0635};
	float f1a=par[2]*z; 
	if(f1a<75.)
		rho-=par[1]/exp(f1a);
	float f2a=par[4]*z; 
	if(f2a<75.)
		rho-=par[3]/exp(f2a);
	float f3a=par[6]*z; 
	if(f3a<75.)
		rho-=par[5]/exp(f3a);
	rho*=par[0];
	return(rho);
}

double approxTodorDensity(double z){
	const double a=0.0;
	const double b=500.0;
	const double rho_b=0.928291380405426;
	const unsigned int nCoeffs=16;
	const static double coeffs[nCoeffs] = 
	{1.630184163164813,0.1974340731867368,-0.1301940665484578,0.063661880435938,
		-0.02284342020414062,0.008295837629225287,-0.007815359409374134,0.0102804719047748,
		-0.01079129490788108,0.009077891732174149,-0.006456928164827926,0.004035689213371977,
		-0.002274127986029174,0.001175871855042831,-0.0005648884951338422,0.0002544329712929343};
	
	if(z>b)
		return(rho_b);
	double d=0.0, dd=0.0, sv, y, y2;
	y2=2.*(y=(2.*z-a-b)/(b-a));
	for(int i=nCoeffs-1; i>0; i--){
		sv=d;
		d=y2*d-dd+coeffs[i];
		dd=sv;
	}
	return(y*d-dd+.5*coeffs[0]);
}

double approxTodorDensityDerivative(double z){
	const double a=0.0;
	const double b=500.0;
	const unsigned int nCoeffs=16;
	const static double coeffs[nCoeffs] = 
	{0.005183430938414656,-0.004694588496049445,0.003603958352920761,-0.00261148343127412,
		0.002076073222458249,-0.00188049398474162,0.001744239717289238,-0.001505356733091662,
		0.001168533290621849,-0.0008147138589872724,0.0005149250859053101,-0.0002981596058010382,
		0.0001597844351285761,-7.984331914223742e-05,3.749376220412173e-05,-1.657580768724709e-05};
	
	if(z>b)
		return(0.0);
	double d=0.0, dd=0.0, sv, y, y2;
	y2=2.*(y=(2.*z-a-b)/(b-a));
	for(int i=nCoeffs-1; i>0; i--){
		sv=d;
		d=y2*d-dd+coeffs[i];
		dd=sv;
	}
	return(y*d-dd+.5*coeffs[0]);
}


todorLinearRefractiveIndex::todorLinearRefractiveIndex(double a, double b):
nref0(a),rho0(b){}

double todorLinearRefractiveIndex::indexOfRefraction(double z) const{
	if(z<0.0)
		return(1.0);
	return(1.+((nref0-1.)/rho0)*approxTodorDensity(z));
}
double todorLinearRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z<0.0 || z>500.0)
		return(0.0);
	return(approxTodorDensityDerivative(z));
}
void todorLinearRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z<0.0){
		n=1.0;
		dndz=0.0;
	}
	else if(z>500.0){
		n=1.+((nref0-1.)/rho0)*todorDensity(z);
		dndz=0.0;
	}
	else{
		n=1.+((nref0-1.)/rho0)*todorDensity(z);
		dndz=((nref0-1.)/rho0)*approxTodorDensityDerivative(z);
	}
}

todorChiRefractiveIndex::todorChiRefractiveIndex(double a, double b):
nref0(a),rho0(b){}

double todorChiRefractiveIndex::indexOfRefraction(double z) const{
	if(z<0.0)
		return(1.0);
	double a=(nref0*nref0-1)/rho0;
	return(sqrt(1.+a*approxTodorDensity(z)));
}
double todorChiRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z<0.0 || z>500.0)
		return(0.0);
	double a=(nref0*nref0-1)/rho0;
	return(0.5*(a*approxTodorDensityDerivative(z))/sqrt(1.+a*approxTodorDensity(z)));
}
void todorChiRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z<0.0){
		n=1.0;
		dndz=0.0;
	}
	double a=(nref0*nref0-1)/rho0;
	if(z>500.0){
		n=sqrt(1.+a*approxTodorDensity(z));
		dndz=0.0;
	}
	else{
		n=sqrt(1.+a*approxTodorDensity(z));
		dndz=0.5*a*approxTodorDensityDerivative(z)/n;
	}
}


todorLLRefractiveIndex::todorLLRefractiveIndex(double a, double b):
nref0(a),rho0(b){}

double todorLLRefractiveIndex::indexOfRefraction(double z) const{
	if(z<0.0)
		return(1.0);
	double a=(nref0*nref0-1)/(nref0*nref0+2)/rho0;
	double a_rho=a*approxTodorDensity(z);
	return(sqrt((1.+2.*a_rho)/(1.-a_rho)));
}
double todorLLRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z<0.0 || z>500.0)
		return(0.0);
	double a=(nref0*nref0-1)/(nref0*nref0+2)/rho0;
	double a_rho=a*approxTodorDensity(z);
	return(3*a/(sqrt((1.+2.*a_rho)/(1.-a_rho))*(2.*a_rho*a_rho-4.*a_rho+2))*approxTodorDensityDerivative(z));
}
void todorLLRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z<0.0){
		n=1.0;
		dndz=0.0;
	}
	double a=(nref0*nref0-1)/(nref0*nref0+2)/rho0;
	double a_rho=a*approxTodorDensity(z);
	if(z>500.0){
		n=sqrt((1.+2.*a_rho)/(1.-a_rho));
		dndz=0.0;
	}
	else{
		n=sqrt((1.+2.*a_rho)/(1.-a_rho));
		dndz=(3*a/(n*(2.*a_rho*a_rho-4.*a_rho+2))*approxTodorDensityDerivative(z));
	}
}


double negligibleAttenuationModel::attenuationLength(double z, double frequency) const{
	return(std::numeric_limits<double>::infinity());
}


double basicAttenuationModel::temperature(double z) const{
	return(-51.5 + z*(4.5319e-3 + 5.822e-6*z));
}

double basicAttenuationModel::attenuationLength(double z, double frequency) const{
	if(z<0.0)
		return(std::numeric_limits<double>::infinity());
	double t = temperature(z);
	const double f0=0.0001, f2=3.16;
	const double w0=log(f0), w1=0.0, w2=log(f2), w=log(frequency);
	const double b0=-6.74890+t*(0.026709-t*0.000884);
	const double b1=-6.22121-t*(0.070927+t*0.001773);
	const double b2=-4.09468-t*(0.002213+t*0.000332);
	double a,bb;
	if(frequency<1.){
		a=(b1*w0-b0*w1)/(w0-w1);
		bb=(b1-b0)/(w1-w0);
	}
	else{
		a=(b2*w1-b1*w2)/(w1-w2);
		bb=(b2-b1)/(w2-w1);
	}
	return 1./exp(a+bb*w);
}
