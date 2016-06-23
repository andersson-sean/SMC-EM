#include "NonparaxialDebyeWidefieldPSF.h"
#include "cubature.h"
#include <boost/math/special_functions/bessel.hpp>

// An adapter method to calculate the real part of the Debye integral
int NonparaxialDebyeWidefieldPSF_cubature::AdapterMethod_realPart(
	unsigned inputDim, 
	const double *input, 
	void *param, 
	unsigned outputDim, 
	double *output) {
	
	// Parameters are assumed to be of the form of an AdapterStruct. 
	// If they are not of this form, then bad things will certainly happen.
	AdapterStruct* adaptedParam = static_cast<AdapterStruct*>(param);
	
	// Evaluate the real part of the PSF for a given input specified by cubature
	(*output) = (*(adaptedParam->ptr)).realPart(adaptedParam->input,*input);
	
	return 0;
}

// The class method to calculate the real part of the Debye integral
double NonparaxialDebyeWidefieldPSF::realPart
	(const std::vector<double>& input, const double& theta) {
	using namespace boost::math;
	
	const double x = input[0];
	const double y = input[1];
	const double z = input[2];
	return sin(theta) * sqrt(cos(theta)) * cos(_kappa * z * cos(theta)) *
		cyl_bessel_j(0,_kappa * sqrt(x*x + y*y) * sin(theta));
}

// An adapter method to calculate the imaginary part of the Debye integral
int NonparaxialDebyeWidefieldPSF_cubature::AdapterMethod_imagPart(
	unsigned inputDim, 
	const double *input, 
	void *param, 
	unsigned outputDim, 
	double *output) {
	
	// Parameters are assumed to be of the form of an AdapterStruct. 
	// If they are not of this form, then bad things will certainly happen.
	AdapterStruct* adaptedParam = static_cast<AdapterStruct*>(param);
	
	// Evaluate the real part of the PSF for a given input specified by cubature
	(*output) = (*(adaptedParam->ptr)).imagPart(adaptedParam->input,*input);
	
	return 0;
}

// The class method to calculate the real part of the Debye integral
double NonparaxialDebyeWidefieldPSF::imagPart
	(const std::vector<double>& input, const double& theta) {
	using namespace boost::math;
	
	const double x = input[0];
	const double y = input[1];
	const double z = input[2];
	return sin(theta) * sqrt(cos(theta)) * sin(_kappa * z * cos(theta)) *
		cyl_bessel_j(0,_kappa * sqrt(x*x + y*y) * sin(theta));
}

// Evaluate the PSF
double NonparaxialDebyeWidefieldPSF::operator()
	(const std::vector<double>& origin, const std::vector<double>& sensor) {

	// Integration parameters:
	const double absError = 0; // absolute error
	const double relError = 1e-2; // relative error
	const size_t maxEval = 0; // maximum number of evaluations
	const unsigned outputDim = 1; // output dimension of integrand
	
	// Integration boundaries:
	const unsigned boundDim = 1; // dimension of boundary
	const double boundMin[1] = {0}; // lower bound on boundary
	const double boundMax[1] = {_alpha}; // upper bound on boundary

	// Create an adapter structure so we can call a C++ member function within
  	// cubature (which is a C library that does not support member functions).
	NonparaxialDebyeWidefieldPSF_cubature::AdapterStruct adaptedParam;
	adaptedParam.input.resize(sensor.size());
	adaptedParam.input[0] = sensor[0]-origin[0];
	adaptedParam.input[1] = sensor[1]-origin[1];
	adaptedParam.input[2] = sensor[2]-origin[2];
	adaptedParam.ptr = this;

	// Evaluate real part:
	double realVal = 0.0;
	double err = 0.0;
	pcubature(outputDim, 
		NonparaxialDebyeWidefieldPSF_cubature::AdapterMethod_realPart, 
		static_cast<void*>(&adaptedParam), 
		boundDim, boundMin, boundMax, 
		maxEval, absError, relError, ERROR_INDIVIDUAL, &realVal, &err);
  
  	// Evaluate imaginary part:
	double imagVal = 0.0;
	err = 0.0;
	pcubature(outputDim, 
		NonparaxialDebyeWidefieldPSF_cubature::AdapterMethod_imagPart, 
		static_cast<void*>(&adaptedParam), 
		boundDim, boundMin, boundMax, 
		maxEval, absError, relError, ERROR_INDIVIDUAL, &imagVal, &err);
  
  	// Take the norm and then scale it:
  	double val = _scale*(realVal*realVal + imagVal*imagVal);
  	
  	// Never return a negative value:
  	return (val >= 0)? val : 0;
}