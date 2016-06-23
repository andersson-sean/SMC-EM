#include "ConfinedDiffusionModel.h"
#include <boost/random/uniform_01.hpp>
#include "../Utilities/Utilities.h"
#include <limits>

// Evaluates the transition function given the current position and the previous position.
// For a confined diffusion, this is given by the solution of the diffusion equation
// for a finite slab; see (Carslaw and Jaeger, 1959, pg.  101 and 114-11) for a 
// derivation.
double ConfinedDiffusionModel::Evaluate
	(const double& currentPosition, const double& previousPosition) const{
	return Evaluate(currentPosition, previousPosition, _diffCoeff, _length);
}

double ConfinedDiffusionModel::Evaluate
	(const double& currentPosition, const double& previousPosition, const double& diffCoeff, const double& length) const{
	const double tol = 1e-10; // Just hard code this for now.
	
	// Input positions are assumed to be between [Center - L/2, Center + L/2],
	// but the pdf assumes they are between [0,L].  We therefore need to subtract
	// "Center - L/2" to recenter.
	const double shiftCurrPos = currentPosition - (_center - length / 2.0);
	const double shiftPrevPos = previousPosition - (_center - length / 2.0);
	
	double pdf = 1.0 / length; // initial condition
	unsigned int m = 1; // initial series term number
	bool addTerms = true; // keep adding terms until tolerance is achieved
	while(addTerms) {
		double expTerm = exp(-diffCoeff * _dt * (m*M_PI / length)*(m*M_PI / length));
		if(fabs(expTerm) < tol)
			addTerms = false; // tolerance has been met
		pdf += (2.0 / length) * expTerm * 
			cos(m * M_PI * shiftCurrPos / length) * 
        	cos(m * M_PI * shiftPrevPos / length);
    	m = m+1;
	}
	return pdf;	
}

// Equivalent to ConfinedDiffusionModel::Evaluate but with derivative wrt "diffCoeff".
double ConfinedDiffusionModel::EvaluateDerivativeWRTDiffCoeff
	(const double& currentPosition, const double& previousPosition, const double& diffCoeff, const double& length) const{
	const double tol = 1e-10; // Just hard code this for now.
	
	// Input positions are assumed to be between [Center - L/2, Center + L/2],
	// but the pdf assumes they are between [0,L].  We therefore need to subtract
	// "Center - L/2" to recenter.
	const double shiftCurrPos = currentPosition - (_center - length / 2.0);
	const double shiftPrevPos = previousPosition - (_center - length / 2.0);
	
	double pdf = 0.0; // initial condition
	unsigned int m = 1; // initial series term number
	bool addTerms = true; // keep adding terms until tolerance is achieved
	while(addTerms) {
		double expTerm = (-_dt * (m*M_PI / length)*(m*M_PI / length)) * 
			exp(-diffCoeff * _dt * (m*M_PI / length)*(m*M_PI / length));
		if(fabs(expTerm) < tol)
			addTerms = false; // tolerance has been met
		pdf = pdf + (2.0 / length) * expTerm * 
			cos(m * M_PI * shiftCurrPos / length) * 
        	cos(m * M_PI * shiftPrevPos / length);
    	m = m+1;
	}
	return pdf;	
}


// Equivalent to ConfinedDiffusionModel::Evaluate but with derivative wrt "currentPosition".
double ConfinedDiffusionModel::EvaluateDerivativeWRTPos
	(const double& currentPosition, const double& previousPosition) const{
	const double tol = 1e-10; // Just hard code this for now.
	
	// Input positions are assumed to be between [Center - L/2, Center + L/2],
	// but the pdf assumes they are between [0,L].  We therefore need to subtract
	// "Center - L/2" to recenter.
	const double shiftCurrPos = currentPosition - (_center - _length / 2.0);
	const double shiftPrevPos = previousPosition - (_center - _length / 2.0);
	
	double pdf = 0.0;
	bool addTerms = true;
	unsigned int m = 1;
	while(addTerms) {
		double expTerm = (m * M_PI / _length) * 
			exp(-_diffCoeff * _dt * (m*M_PI / _length)*(m*M_PI / _length));
		if(fabs(expTerm) < tol)
			addTerms = false;
		pdf = pdf + (2.0 / _length) * expTerm * 
			-sin(m * M_PI * shiftCurrPos / _length) * 
        	cos(m * M_PI * shiftPrevPos / _length);
    	m = m+1;
	}
	return pdf;	
}

// Equivalent to ConfinedDiffusionModel::Evaluate but with second derivative wrt 
// "currentPosition".
double ConfinedDiffusionModel::EvaluateSecondDerivativeWRTPos
	(const double& currentPosition, const double& previousPosition) const{
	const double tol = 1e-10; // Just hard code this for now.
	
	// Input positions are assumed to be between [Center - L/2, Center + L/2],
	// but the pdf assumes they are between [0,L].  We therefore need to subtract
	// "Center - L/2" to recenter.
	const double shiftCurrPos = currentPosition - (_center - _length / 2.0);
	const double shiftPrevPos = previousPosition - (_center - _length / 2.0);
	
	double pdf = 0.0;
	bool addTerms = true;
	unsigned int m = 1;
	while(addTerms) {
		double expTerm = (m * M_PI / _length) * (m * M_PI / _length) *
			exp(-_diffCoeff * _dt * (m*M_PI / _length)*(m*M_PI / _length));
		if(fabs(expTerm) < tol)
			addTerms = false;
		pdf = pdf + (2.0 / _length) * expTerm * 
			-cos(m * M_PI * shiftCurrPos / _length) * 
        	cos(m * M_PI * shiftPrevPos / _length);
    	m = m+1;
	}
	return pdf;	
}

// Equivalent to ConfinedDiffusionModel::Evaluate but with third derivative wrt 
// "currentPosition".
double ConfinedDiffusionModel::EvaluateThirdDerivativeWRTPos
	(const double& currentPosition, const double& previousPosition) const{
	const double tol = 1e-10; // Just hard code this for now.
	
	// Input positions are assumed to be between [Center - L/2, Center + L/2],
	// but the pdf assumes they are between [0,L].  We therefore need to subtract
	// "Center - L/2" to recenter.
	const double shiftCurrPos = currentPosition - (_center - _length / 2.0);
	const double shiftPrevPos = previousPosition - (_center - _length / 2.0);
	
	double pdf = 0.0;
	bool addTerms = true;
	unsigned int m = 1;
	while(addTerms) {
		double expTerm = (m * M_PI / _length) * (m * M_PI / _length) * (m * M_PI / _length) * 
			exp(-_diffCoeff * _dt * (m*M_PI / _length)*(m*M_PI / _length));
		if(fabs(expTerm) < tol)
			addTerms = false;
		pdf = pdf + (2.0 / _length) * expTerm * 
			sin(m * M_PI * shiftCurrPos / _length) * 
        	cos(m * M_PI * shiftPrevPos / _length);
    	m = m+1;
	}
	return pdf;	
}

// Generates a sample from a confined diffusion given the previous position.
// This is done via importance sampling.
double ConfinedDiffusionModel::Sample(const double& previousPosition){
	
	// Rejection sampling (RS) works by drawing from a pdf g(x) over the 
	// domain of interest. The domain of interest here is [0,L].  It 
	// is "easy" to uniformly sample over this domain, so g(x) = 1/L.
    // RS also requires that f(x) < k*g(x), where f(.) is the pdf from which 
    // we want to sample.  Hence, we choose k = L*max(f(.)).
    
    // Calculate the maximum of the pdf given the previous sample using 
    // the tangent hyperbolas method (aka "Halley's Method"... see the 
    // Wolfram page on it.)
    // This method is slower than bisection, but it tends to perform better
    // near the boundaries.
	
	// Evaluate the maximum of the confined pdf given the previous position.
	unsigned int numIter = 1; 
	unsigned int maxIter = 1000;
    double thresh = 1e-16;
    double tmp = previousPosition; // initial condition is the previous location
    while(numIter < maxIter) {
        double d0 = EvaluateDerivativeWRTPos(tmp,previousPosition);
        double d1 = EvaluateSecondDerivativeWRTPos(tmp,previousPosition);
        double d2 = EvaluateThirdDerivativeWRTPos(tmp,previousPosition);
        double inc = -2*d0*d1 / (2*d1*d1 - d0*d2);
        tmp += inc;
        //std::cout << d0 << " " << d1 << " " << d2 << std::endl;
        if(fabs(inc) < thresh)
            break;
        numIter = numIter + 1;
    }
    //std::cout << Evaluate(tmp-0.2,previousPosition) << " " << Evaluate(tmp,previousPosition) << " " << Evaluate(tmp+0.2,previousPosition) << std::endl;
    
    double boundVal = Evaluate(tmp,previousPosition) / (1.0 / _length);
    
    // Rejection sampling
    boost::random::uniform_01<> unif; // continuous (0,1) uniform distribution
    bool generate = true;
    double candidate = 0;
    while(generate) {
        // Generate a u ~ Uniform(0,1)
        double u = unif(gen); 
        
        // Generate a candidate ~ Uniform(center - L/2 , center + L/2)
        candidate = _center + _length*(unif(gen) - 0.5);
        
        // Calculate f(candidate)
        double fx = Evaluate(candidate,previousPosition);
        
        // If u < f(x)/(k*g(x)), keep the sample; otherwise, reject and try again.
        if(u < fx / (boundVal / _length))
            generate = 0;
	}
	return candidate;
}

// Given data from the particle filter/smoother, evaluate the next best parameters
void ConfinedDiffusionModel::UpdateParameters(){
	switch(_axis){
		case DataManager::X:
			_length = CalculateLength(_dm->xParticle);
			_diffCoeff = CalculateDiffusionCoefficient(_dm->xParticle);
			std::cout << "X Diffusion Coefficient: " << _diffCoeff << std::endl;
			std::cout << "X Confinement Length: " << _length << std::endl;
			break;
		case DataManager::Y:
			_length = CalculateLength(_dm->yParticle);
			_diffCoeff = CalculateDiffusionCoefficient(_dm->yParticle);
			std::cout << "Y Diffusion Coefficient: " << _diffCoeff << std::endl;
			std::cout << "Y Confinement Length: " << _length << std::endl;
			break;
		case DataManager::Z:
			_length = CalculateLength(_dm->zParticle);
			_diffCoeff = CalculateDiffusionCoefficient(_dm->zParticle);
			std::cout << "Z Diffusion Coefficient: " << _diffCoeff << std::endl;
			std::cout << "Z Confinement Length: " << _length << std::endl;
			break;
		default:
			break;
	}

}

// Utility function used by UpdateParameters that updates the length parameter.
double ConfinedDiffusionModel::CalculateLength
	(const std::vector<std::vector<double> >& particle) {
	// Initial conditions always start at the largest and smallest values for a double:
	double minVal = std::numeric_limits<double>::max();
	double maxVal = std::numeric_limits<double>::min();
	// Iterate through every position and determine whether is the largest or the
	// smallest we've seen so far.
	for (size_t t = 0; t < _dm->NumTimesteps(); ++t) {
		for (size_t i = 0; i < _dm->NumParticles(); ++i) {
			minVal = fmin(particle[t][i],minVal);
			maxVal = fmax(particle[t][i],maxVal);
		}
	}
	// The optimal length parameter is the deviation between the maximum particle location
	// and the minimum particle location.
	return maxVal - minVal;
}

// Utility function used by UpdateParameters that updates the diffCoeff parameter.
double ConfinedDiffusionModel::CalculateDiffusionCoefficient
	(const std::vector<std::vector<double> >& particle) {
	// Find the maximizing diffusion coefficient.  There is no analytical solution,
  	// so a numerical method must be used.  
  	// Note that this value depends on L, so that value should be updated first.
  	// Here, we use the bisection method to calculate the maximizing diffusion coefficient.
  	const double paramThresh = 1e-10; // if length of interval is smaller, then quit
  	const double derivThresh = 1e-10; // if abs(derivative) is smaller, then quit
  	const int maxIter = 100; // maximum number of iterations
  	double left = _diffCoeff * 0.01; // initial value for left boundary
  	double right = _diffCoeff * 10; // initial value for right boundary
  	// Note that the choice of the left and right boundaries don't really make much of 
  	// a difference.  If the interval [left,right] does not contain the maximizer, then
  	// bisection should get close to either "left" or "right"; although the result
  	// may not be the maximizer, it is still an improvement over the previous value and
  	// EM, in most cases, is robust to this and should continue to improve the estimate
  	// in future iterations.
  	double center = 0.0;
	for (int iter = 0; iter < maxIter; ++iter) {
		center = (left + right) / 2.0;
		std::cout << left << " " << center << " " << right << std::endl;
    	double centerVal = 0; 
    	double leftVal = 0;
    	double rightVal = 0;
    	// Evaluate the derivative at both the left and center locations
    	for (size_t t = 0; t < _dm->NumTimesteps()-1; t++) {
	  		for (size_t i = 0; i < _dm->NumParticles(); i++) {
	      		for(size_t j = 0; j < _dm->NumParticles(); j++) {
		  			leftVal += _dm->wJoint[t][i][j] * 
		  				EvaluateDerivativeWRTDiffCoeff(particle[t+1][j], particle[t][i], left, _length) / 
		  				Evaluate(particle[t+1][j], particle[t][i], left, _length);
		  				
		  			centerVal += _dm->wJoint[t][i][j] * 
		  				EvaluateDerivativeWRTDiffCoeff(particle[t+1][j], particle[t][i], center, _length) / 
		  				Evaluate(particle[t+1][j], particle[t][i], center, _length);
		  				
		  			rightVal += _dm->wJoint[t][i][j] * 
		  				EvaluateDerivativeWRTDiffCoeff(particle[t+1][j], particle[t][i], right, _length) / 
		  				Evaluate(particle[t+1][j], particle[t][i], right, _length);
				}
	    	}
		}
		std::cout << leftVal << " " << centerVal << " " << rightVal << std::endl << std::endl;
		// Determine whether to terminate the algorithm
		if ( fabs(centerVal) < derivThresh || (right - left) / 2.0 < paramThresh )
			break;
		else { // Continue updating the parameter:		
			if ( emlib::sgn(centerVal) == emlib::sgn(leftVal) )
				left = center;
			else
				right = center;
		}
	}
	return center; 
}

// Output stream implementation
void ConfinedDiffusionModel::print(std::ostream& os) const {
	os << _diffCoeff << " " << _length << " ";
}
