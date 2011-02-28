#include "RayTrace.h"
#include <cmath>
#include <iterator>
#include <algorithm>
#include <limits>
#include <queue>

namespace RayTrace{
	
	const double TraceRecord::noReflection(100.0);
	
	//========== minimalRayPosition ==========//
	
	minimalRayPosition::minimalRayPosition():
	x(0.0),z(0.0),theta(0.0){}
	
	minimalRayPosition::minimalRayPosition(double x_, double z_, double theta_):
	x(x_),z(z_),theta(theta_){}
	
	minimalRayPosition& minimalRayPosition::operator +=(const minimalRayPosition& p){
		x+=p.x;
		z+=p.z;
		theta+=p.theta;
		return(*this);
	}
	
	minimalRayPosition& minimalRayPosition::operator -=(const minimalRayPosition& p){
		x-=p.x;
		z-=p.z;
		theta-=p.theta;
		return(*this);
	}
	
	minimalRayPosition& minimalRayPosition::operator *=(double m){
		x*=m;
		z*=m;
		theta*=m;
		return(*this);
	}
	
	const minimalRayPosition minimalRayPosition::operator +(const minimalRayPosition& p) const{
		minimalRayPosition res=*this;
		return(res+=p);
	}
	
	const minimalRayPosition minimalRayPosition::operator -(const minimalRayPosition& p) const{
		minimalRayPosition res=*this;
		return(res-=p);
	}
	
	const minimalRayPosition minimalRayPosition::operator *(double m) const{
		minimalRayPosition res=*this;
		return(res*=m);
	}
	
	const minimalRayPosition operator *(double m, const minimalRayPosition& p){
		return(p*m);
	}
	
	minimalRayPosition abs(const minimalRayPosition& p){
		return(minimalRayPosition(std::abs(p.x), std::abs(p.z), std::abs(p.theta)));
	}
	
	void minimalRayPosition::makeTiny(double s){
		x=s;
		z=s;
		theta=s;
	}
	
	void minimalRayPosition::giveData(TraceRecord& trace) const{
		trace.pathTime=-1.0;
		trace.receiptAngle=theta;
		trace.attenuation=-1.0;
	}
	
	//========== rayPosition ==========//
	
	rayPosition::rayPosition():
	x(0.0),z(0.0),theta(0.0),time(0.0){}
	
	rayPosition::rayPosition(double x_, double z_, double theta_, double time_):
	x(x_),z(z_),theta(theta_),time(time_){}
	
	rayPosition& rayPosition::operator +=(const rayPosition& p){
		x+=p.x;
		z+=p.z;
		theta+=p.theta;
		time+=p.time;
		return(*this);
	}
	
	rayPosition& rayPosition::operator -=(const rayPosition& p){
		x-=p.x;
		z-=p.z;
		theta-=p.theta;
		time-=p.time;
		return(*this);
	}
	
	rayPosition& rayPosition::operator *=(double m){
		x*=m;
		z*=m;
		theta*=m;
		time*=m;
		return(*this);
	}
	
	const rayPosition rayPosition::operator +(const rayPosition& p) const{
		rayPosition res=*this;
		return(res+=p);
	}
	
	const rayPosition rayPosition::operator -(const rayPosition& p) const{
		rayPosition res=*this;
		return(res-=p);
	}
	
	const rayPosition rayPosition::operator *(double m) const{
		rayPosition res=*this;
		return(res*=m);
	}
	
	const rayPosition operator *(double m, const rayPosition& p){
		return(p*m);
	}
	
	rayPosition abs(const rayPosition& p){
		return(rayPosition(std::abs(p.x), std::abs(p.z), std::abs(p.theta), std::abs(p.time)));
	}
	
	void rayPosition::makeTiny(double s){
		x=s;
		z=s;
		theta=s;
		time=s;
	}
	
	void rayPosition::giveData(TraceRecord& trace) const{
		trace.pathTime=time;
		trace.receiptAngle=theta;
		trace.attenuation=-1.0;
	}

	//========== fullRayPosition ==========//
	
	fullRayPosition::fullRayPosition():
	x(0.0),z(0.0),theta(0.0),time(0.0),attenuation(1.0){}
	
	fullRayPosition::fullRayPosition(double x_, double z_, double theta_, double time_, double attenuation_):
	x(x_),z(z_),theta(theta_),time(time_),attenuation(attenuation_){}
	
	fullRayPosition& fullRayPosition::operator +=(const fullRayPosition& p){
		x+=p.x;
		z+=p.z;
		theta+=p.theta;
		time+=p.time;
		attenuation+=p.attenuation;
		return(*this);
	}
	
	fullRayPosition& fullRayPosition::operator -=(const fullRayPosition& p){
		x-=p.x;
		z-=p.z;
		theta-=p.theta;
		time-=p.time;
		attenuation-=p.attenuation;
		return(*this);
	}
	
	fullRayPosition& fullRayPosition::operator *=(double m){
		x*=m;
		z*=m;
		theta*=m;
		time*=m;
		attenuation*=m;
		return(*this);
	}
	
	const fullRayPosition fullRayPosition::operator +(const fullRayPosition& p) const{
		fullRayPosition res=*this;
		return(res+=p);
	}
	
	const fullRayPosition fullRayPosition::operator -(const fullRayPosition& p) const{
		fullRayPosition res=*this;
		return(res-=p);
	}
	
	const fullRayPosition fullRayPosition::operator *(double m) const{
		fullRayPosition res=*this;
		return(res*=m);
	}
	
	const fullRayPosition operator *(double m, const fullRayPosition& p){
		return(p*m);
	}
	
	fullRayPosition abs(const fullRayPosition& p){
		return(fullRayPosition(std::abs(p.x), std::abs(p.z), std::abs(p.theta), std::abs(p.time), std::abs(p.attenuation)));
	}
	
	void fullRayPosition::makeTiny(double s){
		x=s;
		z=s;
		theta=s;
		time=s;
		attenuation=s;
	}
	
	void fullRayPosition::giveData(TraceRecord& trace) const{
		trace.pathTime=time;
		trace.receiptAngle=theta;
		trace.attenuation=attenuation;
	}

	//========== TraceFinder ==========//

	const double TraceFinder::maximum_ice_depth = 2850.0;
	
	///\brief Computes the derivatives of the position coordinates with respect to path length
	///
	///This specialization computes the minimum usable set of coordinates: the radial and vertical 
	///position and the angle. 
	///
	///\param pos The position at which the derivatives are to be computed
	///\param der The object into which to record the calculated derivatives
	///\param frequency The frequency of the signal being propogated, in GHz
	template<>
	void TraceFinder::computeRayDerivatives<minimalRayPosition>(const minimalRayPosition& pos, minimalRayPosition& der, double frequency) const{
		//put dn/dz into der.theta
		double n;
		rModel->indexOfRefractionWithDerivative(pos.z,n,der.theta);
		der.x = sin(pos.theta);
		der.z = cos(pos.theta);
		der.theta = -(der.x * der.theta)/n;
	}
	
	///\brief Computes the derivatives of the position coordinates with respect to path length
	///
	///This specialization computes an intermediate set of coordinates: the radial and 
	///vertical position, the angle, and the time of flight. 
	///
	///\param pos The position at which the derivatives are to be computed
	///\param der The object into which to record the calculated derivatives
	///\param frequency The frequency of the signal being propogated, in GHz
	template<>
	void TraceFinder::computeRayDerivatives<rayPosition>(const rayPosition& pos, rayPosition& der, double frequency) const{
		//put n int der.time, put dn/dz into der.theta
		rModel->indexOfRefractionWithDerivative(pos.z,der.time,der.theta);
		der.x = sin(pos.theta);
		der.z = cos(pos.theta);
		der.theta = -(der.x * der.theta)/der.time;
		der.time /= speedOfLight;
	}
	
	///\brief Computes the derivatives of the position coordinates with respect to path length
	///
	///This specialization computes a complete set of coordinates: the radial and vertical 
	///position, the angle, the time of flight, and the amplitude correction. 
	///
	///\param pos The position at which the derivatives are to be computed
	///\param der The object into which to record the calculated derivatives
	///\param frequency The frequency of the signal being propogated, in GHz
	template<>
	void TraceFinder::computeRayDerivatives<fullRayPosition>(const fullRayPosition& pos, fullRayPosition& der, double frequency) const{
		//put n int der.time, put dn/dz into der.theta
		rModel->indexOfRefractionWithDerivative(pos.z,der.time,der.theta);
		der.x = sin(pos.theta);
		der.z = cos(pos.theta);
		der.theta = -(der.x * der.theta)/der.time;
		der.time /= speedOfLight;
		der.attenuation=-pos.attenuation/aModel->attenuationLength(pos.z,frequency);
		//std::cout << "  eval with z=" << pos.z << " att=" << pos.attenuation << " yields " << der.attenuation << std::endl;
	}

	template<typename positionType>
	void TraceFinder::rkStep(const positionType& pos, const typename positionType::derivativeType& der, const double h, positionType& newPos, positionType& errors, const double frequency) const{
		//const double a[5] = {0.2, 0.3, 0.6, 1.0, 0.875};
		const double b[15] = {0.2, 0.075, 0.3, -11.0/54.0, 1631.0/55296.0, 
			0.225, -0.9, 2.5, 175.0/512.0, 
			1.2, -70.0/27.0, 575.0/13824.0, 
			35.0/27.0, 44275.0/110592.0, 
			253.0/4096.0};
		const double c[6] = {37.0/378.0, 0, 250.0/621.0, 125.0/594.0, 0, 512.0/1771.0};
		const double cp[6] = {c[0]-2825.0/27648.0, 0, c[2]-18575.0/48384.0, c[3]-13525.0/55296.0, c[4]-277.0/14336.0, c[5]-0.25};
		
		positionType temp;
		typename positionType::derivativeType ak2, ak3, ak4, ak5, ak6; //derivatives
		//first step
		temp = pos + h*b[0]*der;
		//second step
		computeRayDerivatives(temp,ak2,frequency);
		temp = pos + h*(b[1]*der + b[5]*ak2);
		//third step
		computeRayDerivatives(temp,ak3,frequency);
		temp = pos + h*(b[2]*der + b[6]*ak2 + b[9]*ak3);
		//fourth step
		computeRayDerivatives(temp,ak4,frequency);
		temp = pos + h*(b[3]*der + b[7]*ak2 + b[10]*ak3 + b[12]*ak4);
		//fifth step
		computeRayDerivatives(temp,ak5,frequency);
		temp = pos + h*(b[4]*der + b[8]*ak2 + b[11]*ak3 + b[13]*ak4 + b[14]*ak5);
		//sixth step
		computeRayDerivatives(temp,ak6,frequency);
		newPos = pos + h*(c[0]*der + c[2]*ak3 + c[3]*ak4 + c[5]*ak6);
		//estimate error
		errors = h*(cp[0]*der + cp[2]*ak3 + cp[3]*ak4 + cp[4]*ak5 + cp[5]*ak6);
	}
	
	template<typename positionType>
	void TraceFinder::rkStep(const positionRecordingWrapper<positionType>& pos, const typename positionRecordingWrapper<positionType>::derivativeType& der, const double h, positionRecordingWrapper<positionType>& newPos, positionRecordingWrapper<positionType>& errors, const double frequency) const{
		//const double a[5] = {0.2, 0.3, 0.6, 1.0, 0.875};
		const double b[15] = {0.2, 0.075, 0.3, -11.0/54.0, 1631.0/55296.0, 
			0.225, -0.9, 2.5, 175.0/512.0, 
			1.2, -70.0/27.0, 575.0/13824.0, 
			35.0/27.0, 44275.0/110592.0, 
			253.0/4096.0};
		const double c[6] = {37.0/378.0, 0, 250.0/621.0, 125.0/594.0, 0, 512.0/1771.0};
		const double cp[6] = {c[0]-2825.0/27648.0, 0, c[2]-18575.0/48384.0, c[3]-13525.0/55296.0, c[4]-277.0/14336.0, c[5]-0.25};
		
		positionRecordingWrapper<positionType> temp;
		typename positionRecordingWrapper<positionType>::derivativeType ak2, ak3, ak4, ak5, ak6; //derivatives
		//first step
		temp = pos + h*b[0]*der;
		temp.recordStep(1);
		//second step
		computeRayDerivatives(temp,ak2,frequency);
		temp = pos + h*(b[1]*der + b[5]*ak2);
		temp.recordStep(2);
		//third step
		computeRayDerivatives(temp,ak3,frequency);
		temp = pos + h*(b[2]*der + b[6]*ak2 + b[9]*ak3);
		temp.recordStep(3);
		//fourth step
		computeRayDerivatives(temp,ak4,frequency);
		temp = pos + h*(b[3]*der + b[7]*ak2 + b[10]*ak3 + b[12]*ak4);
		temp.recordStep(4);
		//fifth step
		computeRayDerivatives(temp,ak5,frequency);
		temp = pos + h*(b[4]*der + b[8]*ak2 + b[11]*ak3 + b[13]*ak4 + b[14]*ak5);
		temp.recordStep(5);
		//sixth step
		computeRayDerivatives(temp,ak6,frequency);
		newPos = pos + h*(c[0]*der + c[2]*ak3 + c[3]*ak4 + c[5]*ak6);
		newPos.takeStepData(temp);
		//estimate error
		errors = h*(cp[0]*der + cp[2]*ak3 + cp[3]*ak4 + cp[4]*ak5 + cp[5]*ak6);
	}
	
	void TraceFinder::replayRkStep(const rkStepRecord& step, const double atten, const double& attenDer, double& newAtten, const double frequency) const{
		const double b[15] = {0.2, 0.075, 0.3, -11.0/54.0, 1631.0/55296.0, 
			0.225, -0.9, 2.5, 175.0/512.0, 
			1.2, -70.0/27.0, 575.0/13824.0, 
			35.0/27.0, 44275.0/110592.0, 
			253.0/4096.0};
		const double c[6] = {37.0/378.0, 0, 250.0/621.0, 125.0/594.0, 0, 512.0/1771.0};
		
		double temp;
		double ak2, ak3, ak4, ak5, ak6; //derivatives
		
		//first step
		temp = atten + step.length*b[0]*attenDer;
		//second step
		ak2=-temp/aModel->attenuationLength(step.z[1],frequency);
		//std::cout << "  eval with z=" << step.z[1] << " att=" << temp << " yields " << ak2 << std::endl;
		temp = atten + step.length*(b[1]*attenDer + b[5]*ak2);
		//third step
		ak3=-temp/aModel->attenuationLength(step.z[2],frequency);
		//std::cout << "  eval with z=" << step.z[2] << " att=" << temp << " yields " << ak3 << std::endl;
		temp = atten + step.length*(b[2]*attenDer + b[6]*ak2 + b[9]*ak3);
		//fourth step
		ak4=-temp/aModel->attenuationLength(step.z[3],frequency);
		//std::cout << "  eval with z=" << step.z[3] << " att=" << temp << " yields " << ak4 << std::endl;
		temp = atten + step.length*(b[3]*attenDer + b[7]*ak2 + b[10]*ak3 + b[12]*ak4);
		//fifth step
		ak5=-temp/aModel->attenuationLength(step.z[4],frequency);
		//std::cout << "  eval with z=" << step.z[4] << " att=" << temp << " yields " << ak5 << std::endl;
		temp = atten + step.length*(b[4]*attenDer + b[8]*ak2 + b[11]*ak3 + b[13]*ak4 + b[14]*ak5);
		//sixth step
		ak6=-temp/aModel->attenuationLength(step.z[5],frequency);
		//std::cout << "  eval with z=" << step.z[5] << " att=" << temp << " yields " << ak6 << std::endl;
		newAtten = atten + step.length*(c[0]*attenDer + c[2]*ak3 + c[3]*ak4 + c[5]*ak6);
	}
	
	///\brief Computes the maximum of a set of coordinate errors, scaled by given factors
	///
	///This specialization computes the minimum usable set of coordinates: the radial and 
	///vertical position and the angle. 
	///
	///\param errors The coordinate errors
	///\param scale The scaling factors for the errors
	///\return The maximum error divided by its associated scale
	template<>
	double TraceFinder::maxError<minimalRayPosition>(const minimalRayPosition& errors, const minimalRayPosition& scale) const{
		double result=std::abs(errors.x/scale.x);
		result=std::max(result,std::abs(errors.z/scale.z));
		result=std::max(result,std::abs(errors.theta/scale.theta));
		return(result);
	}
	
	///\brief Computes the maximum of a set of coordinate errors, scaled by given factors
	///
	///This specialization computes an intermediate set of coordinates: the radial and 
	///vertical position, the angle, and the time of flight. 
	///
	///\param errors The coordinate errors
	///\param scale The scaling factors for the errors
	///\return The maximum error divided by its associated scale
	template<>
	double TraceFinder::maxError<rayPosition>(const rayPosition& errors, const rayPosition& scale) const{
		double result=std::abs(errors.x/scale.x);
		result=std::max(result,std::abs(errors.z/scale.z));
		result=std::max(result,std::abs(errors.theta/scale.theta));
		result=std::max(result,std::abs(errors.time/scale.time));
		return(result);
	}
	
	///\brief Computes the maximum of a set of coordinate errors, scaled by given factors
	///
	///This specialization computes a complete set of coordinates: the radial and vertical 
	///position, the angle, the time of flight, and the amplitude correction. 
	///
	///\param errors The coordinate errors
	///\param scale The scaling factors for the errors
	///\return The maximum error divided by its associated scale
	template<>
	double TraceFinder::maxError<fullRayPosition>(const fullRayPosition& errors, const fullRayPosition& scale) const{
		double result=std::abs(errors.x/scale.x);
		result=std::max(result,std::abs(errors.z/scale.z));
		result=std::max(result,std::abs(errors.theta/scale.theta));
		result=std::max(result,std::abs(errors.time/scale.time));
		result=std::max(result,std::abs(errors.attenuation/scale.attenuation));
		return(result);
	}
	
	template<typename positionType>
	void confirmFinalStepPosition(positionType& pos, const positionType& temp, const double& length){
		pos=temp;
	}
	
	template<typename positionType>
	void confirmFinalStepPosition(positionRecordingWrapper<positionType>& pos, const positionRecordingWrapper<positionType>& temp, const double& length){
		double startZ=pos.z;
		pos=temp;
		pos.setFirstStep(startZ);
		pos.setStepLength(length);
		//pos.dumpData();
	}

	//requires htry>0 !!!
	template<typename positionType>
	void TraceFinder::rkStepControl(double& length, double frequency, positionType& pos, typename positionType::derivativeType& der, const positionType& scale, const double htry, const double eps, double& hdid, double& hnext) const{
		const double SAFETY=0.9;
		//const double growthPower=-0.2;
		//const double shrinkPower=-0.25;
		const double growLimit=1.89e-4; //=(5.0/SAFETY)**(1/growthPower)
		const double shrinkLimit=6561.0; //=(.1/SAFETY)**(1/shrinkPower)
		
		//std::cout << "\tTrying RK(CK) step of size " << htry << std::endl;
		
		double h=htry, errMax;
		positionType errors, temp;
		while(true){
			//take step
			rkStep(pos,der,h,temp,errors,frequency);
			errMax=maxError(errors,scale)/eps;
			if(errMax<=1.0)
				break; //error is tolerable
			//need to decrease step size
			if(errMax < shrinkLimit)
				h=SAFETY*h/sqrt(sqrt(errMax)); // 1/sqrt(sqrt(errMax)) == errMax**shrinkPower
			else
				h*=0.1;
			//std::cout << "\tError too big, reducing step size to " << h << std::endl;
			if((length+h)==length)
				throw std::runtime_error("TraceFinder::rkStepControl: stepsize underflow");
		}
		//increase step size
		if(errMax > growLimit){
			//Here, we use a taylor expansion (about 1) for errMax**growthPower. This 
			//expansion always underestimates the true function (in the domain [0,1])
			//so it will give conservative advice for the next step's size
			errMax-=1.0;
			hnext = SAFETY*h* (1.+errMax*(-.2+errMax*(.12+errMax*(-.088+errMax*.0704))));
		}
		else
			hnext = 5.0*h;
		//std::cout << "\tError acceptable, plan next step size to be " << hnext << std::endl;
		length+=(hdid=h);
		confirmFinalStepPosition(pos,temp,hdid);//pos=temp;
	}
	
	double fresnelReflect(double theta, double& polarization, double n1, double n2){
		double s=sin(theta);
		if(n1*s > n2){ //total internal reflection!
			//std::cout << "Total internal refraction" << std::endl;
			return(1.0);
		}
		double c=cos(theta);
		double d=sqrt(n2*n2 - n1*n1*s*s);
		double rPerpedicular=(n1*c-d)/(n1*c+d);
		double rParallel=(n2*n2*c-n1*d)/(n2*n2*c+n1*d);
		double cp=cos(polarization);
		double sp=sin(polarization);
		double R=sqrt(rPerpedicular*rPerpedicular*cp*cp + rParallel*rParallel*sp*sp);
		polarization=atan((rParallel*sp)/(rPerpedicular*cp));
		return(R);
	}
	
	double fresnelTransmit(double theta, double& polarization, double n1, double n2){
		double s=sin(theta);
		if(n1*s > n2){ //total internal reflection!
			//std::cout << "Total internal refraction" << std::endl;
			return(0.0);
		}
		double c=cos(theta);
		double d=sqrt(n2*n2 - n1*n1*s*s);
		double tPerpedicular=(2.*n1*c)/(n1*c+d);
		double tParallel=(2*n1*n2*c)/(n2*n2*c+n1*d);
		double cp=cos(polarization);
		double sp=sin(polarization);
		double T=sqrt(tPerpedicular*tPerpedicular*cp*cp + tParallel*tParallel*sp*sp);
		polarization=atan((tParallel*sp)/(tPerpedicular*cp));
		return(T);
	}
	
	template <>
	void correctAmplitudeReflect<fullRayPosition>(fullRayPosition& pos, double polarization, double n1, double n2){
		pos.attenuation*=fresnelReflect(pos.theta, polarization, n1, n2);
	}
	
	template <>
	void correctAmplitudeTransmit<fullRayPosition>(fullRayPosition& pos, double polarization, double n1, double n2){
		pos.attenuation*=fresnelTransmit(pos.theta, polarization, n1, n2);
	}
	
	double TraceFinder::recalculateAmplitude(const traceReplayRecord& trace, double frequency, double polarization){
		double attenuation=1.0, attenuationDerivative,nextAttenuation;
		//std::cout << " att: " << attenuation << std::endl;
		for(std::vector<stepRecord>::const_iterator step=trace.steps.begin(), end=trace.steps.end(); step!=end; ++step){
			switch(step->stepType){
				case RK_FIRST_STEP:
					//do nothing; should never occur
					break;
				case RK_AIR_STEP:
					//TODO: implement this!
					break;
				case RK_STEP:
					attenuationDerivative=-attenuation/aModel->attenuationLength(step->rkData.z[0],frequency);
					//std::cout << "  eval with z=" << step->rkData.z[0] << " att=" << attenuation << " yields " << attenuationDerivative << std::endl;
					replayRkStep(step->rkData, attenuation, attenuationDerivative, nextAttenuation, frequency);
					attenuation=nextAttenuation;
					break;
				case RK_REFLECT_STEP:
					if(step->angle<(pi/2.)){ //ray was reflected down form the ice surface
						//use the angle _after_ reflection
						attenuation*=fresnelReflect(step->angle, polarization, rModel->indexOfRefraction(0.0),1.0);
					}
					else{ //ray was reflected up from bedrock
						//use the angle _before_ reflection
						attenuation*=fresnelReflect(pi-step->angle, polarization, rModel->indexOfRefraction(maximum_ice_depth),1.0);
					}
					break;
			}
			//std::cout << " att: " << attenuation << std::endl;
		}
		return(attenuation);
	}

	std::pair<bool, double> TraceFinder::traceMax(double emit_depth, const rayTargetRecord& target, double left, double right) const{
		//unsigned int steps=1;
		std::pair<bool,double> result;
		double a=left; //the left boundary of the search interval
		double b=right; //the right boundary of the search interval
		const double w=(3.0-sqrt(5.0))/2.0;
		double c=a+w*(b-a); //the best estimate point
		//frequency and polarization are unimportant
		double cy=doTrace<minimalRayPosition>(emit_depth,c,target,SurfaceReflection,0.0,0.0).miss;
		double d,dy;
		const double tol=1e-3;
		while(cy<0.0 && fabs(b-a)>tol*c){
			//std::cout << "\tsearch domain is now [" << a << ',' << b << ']' << std::endl;
			if((c-a)>=(b-c)){ //the left subinterval is larger, so search within it
				d=c-w*(c-a);
				dy=doTrace<minimalRayPosition>(emit_depth,d,target,SurfaceReflection,0.0,0.0).miss; //evaluate the miss distance
				//std::cout << "\t " << d << ' ' << dy << std::endl;
				if(dy>cy){
					b=c;
					c=d;
					cy=dy;
				}
				else
					a=d;
			}
			else{ //otherwise seach within the right subinterval
				d=c+w*(b-c);
				dy=doTrace<minimalRayPosition>(emit_depth,d,target,SurfaceReflection,0.0,0.0).miss; //evaluate the miss distance
				//std::cout << "\t " << d << ' ' << dy << std::endl;
				if(dy>cy){
					a=c;
					c=d;
					cy=dy;
				}
				else
					b=d;
			}
			//steps++;
		}
		//std::cout << "traceMax stopped at c=" << c << " f(c)=" << cy << std::endl;
		result.first = cy>0.0; //whether the best point indicates the existance of roots
		result.second=c;
		//std::cout << "traceMax took " << steps << " steps" << std::endl;
		return(result);
	}

	std::pair<bool, double> TraceFinder::traceMin(double emit_depth, const rayTargetRecord& target, double left, double right) const{
		//unsigned int steps=1;
		std::pair<bool,double> result;
		double a=left; //the left boundary of the search interval
		double b=right; //the right boundary of the search interval
		const double w=(3.0-sqrt(5.0))/2.0;
		double c=a+w*(b-a); //the best estimate point
		//frequency and polarization are unimportant
		double cy=doTrace<minimalRayPosition>(emit_depth,c,target,BedrockReflection,0.0,0.0).miss;
		double d,dy;
		const double tol=1e-3;
		while(cy>0.0 && fabs(b-a)>tol*c){
			if((c-a)>=(b-c)){ //the left subinterval is larger, so search within it
				d=c-w*(c-a);
				dy=doTrace<minimalRayPosition>(emit_depth,d,target,BedrockReflection,0.0,0.0).miss; //evaluate the miss distance
				if(dy>cy){
					b=c;
					c=d;
					cy=dy;
				}
				else
					a=d;
			}
			else{ //otherwise seach within the right subinterval
				d=c+w*(b-c);
				dy=doTrace<minimalRayPosition>(emit_depth,d,target,BedrockReflection,0.0,0.0).miss; //evaluate the miss distance
				if(dy>cy){
					a=c;
					c=d;
					cy=dy;
				}
				else
					b=d;
			}
			//steps++;
		}
		//std::cout << "traceMin stopped at c=" << c << " f(c)=" << cy << std::endl;
		result.first = cy<0.0; //whether the best point indicates the existance of roots
		result.second=c;
		//std::cout << "traceMin took " << steps << " steps" << std::endl;
		return(result);
	}

	double TraceFinder::traceRootImpl(double emit_depth, const rayTargetRecord& target, bool rising, unsigned short allowedReflections, double requiredAccuracy, 
									  double a, TraceRecord& aTrace, double c, TraceRecord& cTrace, double angle) const{
		const double miss_eps=requiredAccuracy/100.0;
		const double angle_eps=1e-10;
		double lastMiss;
		TraceRecord trace;
		
		double e=0.0,ep=0.0,p,q,r,s,t;
		
		//initialize these so that the miss change check won't fail accidentally
		trace.miss=1.0e3*miss_eps;
		lastMiss=-1.0e3*miss_eps;
		
		while((c-a)>angle_eps && fabs(trace.miss-lastMiss)>miss_eps){
			//std::cout << "\tAngular range=" << (c-a) << ", miss change=" << fabs(trace.miss-lastMiss) << std::endl;
			lastMiss=trace.miss;
			//std::cout << " Angular range is now [" << a << ',' << c << "]\n";
			//std::cout << " Trying angle=" << angle << std::endl;
			trace=doTrace<minimalRayPosition>(emit_depth,angle,target,allowedReflections,0.0,0.0);
			//steps++;
			//std::cout << " miss distance was " << trace.miss << '\n';
			if(std::abs(trace.miss) < requiredAccuracy){
				//std::cout << "traceRoot took " << steps << " steps" << std::endl;
				//doTrace<fullRayPosition>(emit_depth,angle,target,allowedReflections,frequency,polarization)
				return(angle);
			}
			
			//calculate the quadratic interpolation
			r=trace.miss/cTrace.miss;
			s=trace.miss/aTrace.miss;
			t=aTrace.miss/cTrace.miss;
			p=s*(t*(r-t)*(c-angle)+(r-1)*(angle-a));
			q=(r-1)*(s-1)*(t-1);
			
			//collapse the interval
			if((trace.miss>0.0) != rising){
				a=angle;
				aTrace=trace;
			}
			else{
				c=angle;
				cTrace=trace;
			}
			
			//accept the interpolation only if it falls within the current boundaries
			if((angle+(p/q))>=a && (angle+(p/q))<=c){
				//would like to interpolate, but only do so if recent convergence has been suitably rapid
				if(std::abs(p/q) >= 0.5*std::abs(e)){ //it has not; bisect instead
					//std::cout << " will bisect" << std::endl;
					e=ep;
					ep=0.5*(a+c)-angle;
					angle=0.5*(a+c);
				}
				else{ //it has; use the interpolation
					//std::cout << " will interpolate" << std::endl;
					angle+=p/q;
					e=ep;
					ep=p/q;
				}
			}
			else{ //otherwise, bisect
				//std::cout << " will bisect" << std::endl;
				e=ep;
				ep=0.5*(a+c)-angle;
				angle=0.5*(a+c);
			}
		}
		//std::cout << "traceRoot took " << steps << " steps and bailed out" << std::endl;
		//doTrace<fullRayPosition>(emit_depth,angle,target,allowedReflections,frequency,polarization)
		return(angle);
	}
	
	//rough implementation of Brent's method
	//This implementation neglects some of the checks made by Brent to ensure 
	//rapidity of convergence, and does not bother to use secant interpolation. 
	//In practice it seems not to matter, as it is still capable of the same 
	//best-case behavior as the full Brent's method. 
	//This implementation is more conveniently controllable for this purpose 
	//than NR's zbrent, with regard to the stopping conditions
	double TraceFinder::traceRoot(double emit_depth, const rayTargetRecord& target, double minAngle, double maxAngle, bool rising, unsigned short allowedReflections, double requiredAccuracy) const{
		double a=minAngle,c=maxAngle;
		TraceRecord aTrace, cTrace;
		//frequency and polarization are irrelevant on all traces except the final one
		aTrace = doTrace<minimalRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0);
		cTrace = doTrace<minimalRayPosition>(emit_depth,c,target,allowedReflections,0.0,0.0);
		double angle=0.5*(minAngle+maxAngle);
		//std::cout << " end point miss values are " << aTrace.miss << ',' << cTrace.miss << std::endl;
		return(traceRootImpl(emit_depth, target, rising, allowedReflections, requiredAccuracy, a, aTrace, c, cTrace, angle));
	}
	
	double TraceFinder::refineRoot(double emit_depth, const rayTargetRecord& target, const TraceRecord& seed, bool rising, unsigned short allowedReflections, double requiredAccuracy) const{
		//bracket the root
		double a, c;
		TraceRecord aTrace, cTrace;
		double testDisp=.01;
		const unsigned int maxTests=(unsigned int)std::ceil(0.5*(sqrt(1.0+(8*pi/testDisp))-1.0));
		if((seed.miss>0.0)==rising)
			testDisp*=-1;
		unsigned int i;
		a=seed.launchAngle;
		for(i=0; i<maxTests; i++){
			a+=testDisp;
			testDisp*=2.0;
			if(a<0.0)
				a=0.0;
			else if(a>pi)
				a=pi;
			aTrace=doTrace<minimalRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0);
			if(std::abs(aTrace.miss) < requiredAccuracy)
				return(a);
				//return(doTrace<fullRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0));
			if((aTrace.miss*seed.miss)<0.0)
				break;
		}
		if(i==maxTests){
			//std::cerr << "Last attempt at bracketing was: [" << seed.launchAngle << "->" << seed.miss << ',' << a << "->" << aTrace.miss << ']' << std::endl;
			//std::cerr << "\trising edge = " << std::boolalpha << rising << std::endl;
			throw std::runtime_error("TraceFinder::refineRoot: exceeded maximum allowed number of steps for bracketing.");
		}
		if((seed.miss>0.0)==rising){
			c=seed.launchAngle;
			cTrace=seed;
		}
		else{
			c=a;
			cTrace=aTrace;
			a=seed.launchAngle;
			aTrace=seed;
		}
		
		double angle=0.5*(a+c);
		return(traceRootImpl(emit_depth, target, rising, allowedReflections, requiredAccuracy, a, aTrace, c, cTrace, angle));
	}
	
	bool shorterPath(const TraceRecord& a, const TraceRecord& b){
		return(a.pathLen < b.pathLen);
	}

	std::vector<TraceRecord> TraceFinder::findPaths(Vector sourcePos, Vector targetPos, double frequency, double polarization, unsigned short allowedReflections, double requiredAccuracy, std::vector<traceReplayRecord>* replayBuffer) const{
		std::vector<TraceRecord> results;
		pathRecorder<fullRayPosition>* recorder=NULL;
		if(replayBuffer!=NULL)
			recorder=new pathRecorder<fullRayPosition>();
		//std::cout << "Finding paths from (" << sourcePos.x << ',' << sourcePos.y << ',' << sourcePos.z << ") to (" << targetPos.x << ',' << targetPos.y << ',' << targetPos.z << ')' << std::endl;
		if(sourcePos.GetZ()>maximum_ice_depth || targetPos.GetZ()>maximum_ice_depth)
			return(results);
		bool fullyContained=(sourcePos.GetZ()>=0.0 && targetPos.GetZ()>=0.0); //whether the ray is entirely inside the ice
		double dist = sqrt((targetPos.GetX()-sourcePos.GetX())*(targetPos.GetX()-sourcePos.GetX())+(targetPos.GetY()-sourcePos.GetY())*(targetPos.GetY()-sourcePos.GetY()));
		rayTargetRecord target(targetPos.GetZ(),dist);
		
		//special case for pesky near-vertical rays
		if(dist<=requiredAccuracy){
			//std::cout << "Computing direct ray" << std::endl;
			if(replayBuffer==NULL)
				results.push_back(doVerticalTrace<fullRayPosition>(sourcePos.GetZ(), (sourcePos.GetZ()<targetPos.GetZ()?0.0:pi), target, NoReflection, frequency, polarization));
			else{
				results.push_back(doVerticalTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(), (sourcePos.GetZ()<targetPos.GetZ()?0.0:pi), target, NoReflection, frequency, polarization, recorder));
				replayBuffer->push_back(recorder->getData());
			}
			if((allowedReflections & SurfaceReflection) && fullyContained){
				//std::cout << "Computing surface-reflected ray" << std::endl;
				if(replayBuffer==NULL)
					results.push_back(doVerticalTrace<fullRayPosition>(sourcePos.GetZ(), pi, target, SurfaceReflection, frequency, polarization));
				else{
					results.push_back(doVerticalTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(), pi, target, SurfaceReflection, frequency, polarization, recorder));
					replayBuffer->push_back(recorder->getData());
				}
			}
			if(allowedReflections & BedrockReflection){
				//std::cout << "Computing bedrock-reflected ray" << std::endl;
				if(replayBuffer==NULL)
					results.push_back(doVerticalTrace<fullRayPosition>(sourcePos.GetZ(), 0.0, target, BedrockReflection, frequency, polarization));
				else{
					results.push_back(doVerticalTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(), 0.0, target, BedrockReflection, frequency, polarization, recorder));
					replayBuffer->push_back(recorder->getData());
				}
			}
			std::sort(results.begin(),results.end(),&shorterPath);
			return(results);
		}
		
		indexOfRefractionModel::RayEstimate est;
		if(!fullyContained){
			//TODO: Hande recording!
			results.push_back(findUncontainedFast(sourcePos, targetPos, frequency, polarization, requiredAccuracy));
			return(results);
		}
		if(fullyContained){
			//std::cout << "Looking for estimate" << std::endl;
			est=rModel->estimateRayAngle(sourcePos.GetZ(), targetPos.GetZ(), dist);
		}
		if(est.status==indexOfRefractionModel::SOLUTION){
			//std::cout << "Got fast solution" << std::endl;
			if(replayBuffer==NULL)
				results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(), est.angle, target, NoReflection, frequency, polarization)); //this is supposed to be a direct ray, so no reflections should be needed
			else{
				results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(), est.angle, target, NoReflection, frequency, polarization, recorder)); //this is supposed to be a direct ray, so no reflections should be needed
				replayBuffer->push_back(recorder->getData());
			}
			if(std::abs(results.front().miss) > requiredAccuracy){
				//std::cout << "Fast Solution not close enough; refining" << std::endl;
				est.angle=refineRoot(sourcePos.GetZ(), target, results.front(), true, NoReflection, requiredAccuracy);
				if(replayBuffer==NULL)
					results.front()=doTrace<fullRayPosition>(sourcePos.GetZ(), est.angle, target, NoReflection, frequency, polarization);
				else{
					results.front()=doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(), est.angle, target, NoReflection, frequency, polarization, recorder);
					replayBuffer->front()=recorder->getData();
				}
			}
			if(allowedReflections & SurfaceReflection){
				//std::cout << "Looking for surface reflected solution" << std::endl;
				double angle=traceRoot(sourcePos.GetZ(),target,est.angle+0.01,pi,false,SurfaceReflection,requiredAccuracy);
				if(replayBuffer==NULL)
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization));
				else{
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization,recorder));
					replayBuffer->push_back(recorder->getData());
				}
			}
			if(allowedReflections & BedrockReflection){
				//std::cout << "Looking for bedrock reflected solution" << std::endl;
				double angle=traceRoot(sourcePos.GetZ(),target,0.0,est.angle-0.01,false,BedrockReflection,requiredAccuracy);
				if(replayBuffer==NULL)
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization));
				else{
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization,recorder));
					replayBuffer->push_back(recorder->getData());
				}
			}
		}
		else{
			double a=0.0,b=pi; //the endpoints of a range bracketing a root
			std::pair<bool, double> minRes, maxRes;
			switch(est.status){
				case indexOfRefractionModel::SOLUTION:
					//impossible
					break;
				case indexOfRefractionModel::UPPER_LIMIT:
					b=est.angle;//std::min(pi,est.angle+.1); //add a small amount incase the solution is right at the boundary given by the estimate
					maxRes = traceMax(sourcePos.GetZ(),target,b,pi);
					minRes = traceMin(sourcePos.GetZ(),target,0.0,b);
					break;
				case indexOfRefractionModel::LOWER_LIMIT:
					a=est.angle;//std::max(0.0,est.angle-.1);
					maxRes = traceMax(sourcePos.GetZ(),target,a,pi);
					minRes = traceMin(sourcePos.GetZ(),target,0.0,a);
					break;
				case indexOfRefractionModel::NO_SOLUTION:
					//Note that getting indexOfRefractionModel::NO_SOLUTION implies no _direct_ solution exists
					//if we were instructed to find reflected rays, we still need to do the work to look for them
					
					//fall through
				case indexOfRefractionModel::UNKNOWN:
					maxRes = traceMax(sourcePos.GetZ(),target,0.0,pi);
					minRes = traceMin(sourcePos.GetZ(),target,0.0,pi);
					break;
			}
			
			if(minRes.first && (allowedReflections & BedrockReflection)){
				//std::cout << "Looking for bedrock reflected solution" << std::endl;
				double angle=traceRoot(sourcePos.GetZ(),target,0.0,minRes.second,false,BedrockReflection,requiredAccuracy);
				if(replayBuffer==NULL)
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization));
				else{
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization,recorder));
					replayBuffer->push_back(recorder->getData());
				}
				a=std::max(minRes.second,a);
			}
			if(maxRes.first && (allowedReflections & SurfaceReflection) && fullyContained){
				//std::cout << "Looking for surface reflected solution" << std::endl;
				double angle=traceRoot(sourcePos.GetZ(),target,maxRes.second,pi,false,SurfaceReflection,requiredAccuracy);
				if(replayBuffer==NULL)
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization));
				else{
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization,recorder));
					replayBuffer->push_back(recorder->getData());
				}
				b=std::min(maxRes.second,b);
			}
			if((minRes.first || maxRes.first) && est.status!=indexOfRefractionModel::NO_SOLUTION){
				//std::cout << "Looking for direct solution" << std::endl;
				double angle=traceRoot(sourcePos.GetZ(),target,a,b,true,NoReflection,requiredAccuracy);
				if(replayBuffer==NULL)
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization));
				else{
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization,recorder));
					replayBuffer->push_back(recorder->getData());
				}
				if(std::abs(results.back().miss) > 10.0*requiredAccuracy){ //if it wasn't really a solution, throw it away
					//std::cout << "Supposed solution missed by " << results.back().miss << " meters, rejecting" << std::endl;
					results.pop_back();
					if(replayBuffer!=NULL)
						replayBuffer->pop_back();
				}
			}
		}
		if(results.size() > 1)
			std::sort(results.begin(),results.end(),&shorterPath);
		return(results);
	}
	
	TraceRecord TraceFinder::findUncontainedFast(Vector sourcePos, Vector targetPos, double frequency, double polarization, double requiredAccuracy) const{
		TraceRecord trace;
		//TODO: if both points are above surface, give obvious answers
		if(sourcePos.GetZ()>maximum_ice_depth || targetPos.GetZ()>maximum_ice_depth)
			return(trace);
		//the index of refraction of the ice at the surface
		const double surfaceN = rModel->indexOfRefraction(0.0);
		double wholeDist = sqrt((targetPos.GetX()-sourcePos.GetX())*(targetPos.GetX()-sourcePos.GetX())+(targetPos.GetY()-sourcePos.GetY())*(targetPos.GetY()-sourcePos.GetY()));
		double upper=std::min(sourcePos.GetZ(),targetPos.GetZ()); //upper position in our coordinates has the smaller z coordinate!
		double lower=std::max(sourcePos.GetZ(),targetPos.GetZ());
		rayTargetRecord target(lower,wholeDist);
		double minDist=0.0, maxDist=wholeDist;
		double dist=(minDist+maxDist)/2.;
		double thetaA=0.0,thetaB=0.0;
		while((maxDist-minDist)>requiredAccuracy/10.){
			//std::cout << "Trying dist=" << dist << std::endl;
			//thetaA is the angle of the straight line ray above the ice
			thetaA=atan((wholeDist-dist)/-upper);
			//thetaB is the angle the straight line ray would have below the ice
			thetaB=asin(sin(thetaA)/surfaceN);
			
			indexOfRefractionModel::RayEstimate est=rModel->estimateRayAngle(0.0, lower, dist);
			if(est.status == indexOfRefractionModel::SOLUTION){
				if(est.angle<thetaB)
					minDist=dist;
				else if(est.angle>thetaB)
					maxDist=dist;
				else //we got lucky
					break;
			}
			else{
				target.distance = dist;
				trace = doTrace<minimalRayPosition>(0.0,thetaB,target,NoReflection,0.0,0.0);
				if(std::abs(trace.miss) <= requiredAccuracy)
					break;
				if(trace.miss<0.0)
					maxDist=dist;
				else if(trace.miss>0.0)
					minDist=dist;
			}
			
			dist=(minDist+maxDist)/2.;
		}
		//do the final version trace
		target.distance = dist;
		trace = doTrace<fullRayPosition>(0.0,thetaB,target,NoReflection,frequency,polarization);
		//add on the above ice segment
		double extraLength=sqrt(upper*upper+(wholeDist-dist)*(wholeDist-dist));
		//double oldLenth=trace.pathLen;
		trace.pathLen+=extraLength;
		trace.pathTime+=extraLength/speedOfLight;
		//trace.amplitude=(trace.amplitude*oldLenth)/trace.pathLen;
		//fix up all the angles, and take care of the refraction
		if(sourcePos.GetZ() < targetPos.GetZ()){ //upper is the source depth
			trace.launchAngle=thetaA;
			//receipt angle is as calculated by doTrace
			trace.reflectionAngle=-thetaA;
			trace.attenuation*=fresnelTransmit(thetaA, polarization, 1.0, surfaceN);
		}
		else{ //lower is the source dpeth
			trace.launchAngle=pi-trace.receiptAngle;
			trace.receiptAngle=pi-thetaA;
			trace.reflectionAngle=thetaB-(pi/2.);
			trace.attenuation*=fresnelTransmit(-trace.reflectionAngle, polarization, 1.0, surfaceN);
		}
		return(trace);
	}
	
	double TraceFinder::signalStrength(const TraceRecord& ray, const Vector& src, const Vector& trg, unsigned short allowedReflections) const{
		const double changeThresh=1e-3; //10^(-n) should give about n digits of accuracy
		const double startDelta=.002;
		//std::cout << "Considering focusing from " << src << " to " << trg << " at angle " << ray.launchAngle << std::endl;
		double dist = sqrt((trg.GetX()-src.GetX())*(trg.GetX()-src.GetX())+(trg.GetY()-src.GetY())*(trg.GetY()-src.GetY()));
		if(ray.launchAngle<.01 || ray.launchAngle>(pi-.01)){
			//std::cout << "ray too close to vertical" << std::endl;
			return(-1.0);
		}
		
		//std::cout << "Original ray miss distance: " << ray.miss << std::endl;
		unsigned int steps=1;
		double scaledChange=1e10;
		double lastChange=0.0;
		unsigned int scaleIncrease=0;
		const double multInc=2.0;
		double ratio,accel;
		std::vector<double> terms;
		std::queue<double> past;
		for(double deltaTheta=startDelta; scaledChange>changeThresh; deltaTheta/=multInc){
			//std::cout << "\tdeltaTheta = " << deltaTheta << std::endl;
			//std::cout << "\tlaunch angle will be " << ray.launchAngle-deltaTheta << std::endl;
			TraceRecord testRay = doTrace<minimalRayPosition>(src.GetZ(), ray.launchAngle-deltaTheta, rayTargetRecord(trg.GetZ(),dist), allowedReflections, 0.0, 0.0);
			//std::cout << "\t\ttest ray miss distance: " << testRay.miss << std::endl;
			ratio=std::abs(deltaTheta/(ray.miss-testRay.miss));
			//std::cout << "\t\tdeltaTheta/deltaZ =  " << ratio << std::endl;
			terms.push_back(ratio);
			if(!past.empty()){
				double old=past.back();
				double mult=multInc;
				accel=ratio;
				past.push(ratio);
				for(unsigned int i=0; i<steps-1; i++){
					mult*=multInc;
					accel=(mult*accel-past.front())/(mult-1.0);
					past.pop();
					past.push(accel);
				}
				//std::cout << "\t\textrapolated value = " << accel << std::endl;
				scaledChange=std::abs(old-accel)/accel;
				//std::cout << "\t\tscaled change = " << scaledChange << std::endl;
				if(scaledChange>lastChange){
					//std::cout << "\t\tWARNING: SCALE INCREASE!" << std::endl;
					scaleIncrease++;
					if(scaleIncrease>=3){
						//std::cout << "\t\tWARNING: scale increased over " << scaleIncrease <<  " consecutive steps" << std::endl;
						//scaleIncrease>=n means terms.size()>=n+1
						//accel=(8**(terms.end()-2)-6**(terms.end()-1)+terms.back())/3.0; //n=2 case
						accel=(64.**(terms.end()-3)-56.**(terms.end()-2)+14.**(terms.end()-1)-terms.back())/21.0; //n=3 case
						break;
					}
					//break;
				}
				else
					scaleIncrease=0;
				lastChange=scaledChange;
			}
			else{
				past.push(ratio);
				accel=ratio;
			}
			steps++;
		}
		
		double ir=rModel->indexOfRefraction(trg.GetZ())/rModel->indexOfRefraction(src.GetZ());
		return(sqrt(ir*accel/dist));
	}
	
} //namespace RayTrace
