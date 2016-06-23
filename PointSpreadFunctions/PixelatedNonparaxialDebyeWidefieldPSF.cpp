#include "PixelatedGaussianPSF.h"

#include "cubature.h"

// A method to be called by Cubature for use with the PixelatedGaussianPSF class.
int PixelatedGaussianPSF_cubature::AdapterMethod(
	unsigned inputDim, // dimension of input
	const double *input, // input values 
	void *param, // parameters
	unsigned outputDim, // dimension of output
	double *output) {

	// Parameters are assumed to be of the form of an AdapterStruct. 
	// If they are not of this form, then bad things will certainly happen.
	AdapterStruct* adaptedParam = static_cast<AdapterStruct*>(param);
	
	// Sensor locations are varied over the integration interval as specified
	// by the cubature subroutine.
	std::vector<double> sensorCopy(3,0);
	std::copy(
		adaptedParam->sensor.begin(),adaptedParam->sensor.end(),
		sensorCopy.begin()
		);
	sensorCopy[0] += input[0];
	sensorCopy[1] += input[1];
	
	// Output value is equal to the value of the Gaussian PSF given the
	// particle location and the sensor location (as adjusted by Cubature).
	// This syntax is messy; here's how you should parse it:
	// 1) The GaussianPSF object implements the (vector<double>,vector<double>) method
	// 2) Take adaptedParam and grab its psf
	// 3) The psf acquired in step (2) needs to be dereferenced since it is a pointer.
	// 4) Given the dereferenced pointer, use it to call its own 
	// 		(vector<double>,vector<double>) method on the particle and sensor locations.
	(*output) = (*(adaptedParam->psf))(adaptedParam->particle,sensorCopy);
	
	return 0;	
}

// Evaluate 
double PixelatedGaussianPSF::operator()
	(const std::vector<double>& particle, const std::vector<double>& sensor) {

	// Integration parameters:
	const double absError = 0; // absolute error
	const double relError = 1e-2; // relative error
	const size_t maxEval = 0; // max number of evaluations
	const unsigned outputDim = 1; // output dimension of integrand
	
	// Integration boundaries:
	const unsigned boundDim = 2; // dimension of boundaries
	const double boundMin[2] = {-_pixelSize/2.0, -_pixelSize/2.0};
	const double boundMax[2] = { _pixelSize/2.0,  _pixelSize/2.0};
  
  	// Create an adapter structure so we can call a C++ member function within
  	// cubature (which is a C library that does not support member functions).
  	PixelatedGaussianPSF_cubature::AdapterStruct adaptedParam;
  	adaptedParam.particle.resize(particle.size());
  	std::copy( // copy particle location to adapter
  		particle.begin(),particle.end(),
  		adaptedParam.particle.begin()
  		);
  	adaptedParam.sensor.resize(sensor.size());
  	std::copy( // copy particle location to adapter
  		sensor.begin(),sensor.end(),
  		adaptedParam.sensor.begin()
  		);
  	adaptedParam.psf = &_gaussPSF; // provide adapter with a reference to the PSF

	// Allocate space:
	double val = 0.0; // output value
	double err = 0.0; // output error

	// Integrate:
	pcubature(outputDim, 
		PixelatedGaussianPSF_cubature::AdapterMethod, 
		static_cast<void*>(&adaptedParam), 
		boundDim, boundMin, boundMax, 
	    maxEval, absError, relError, ERROR_INDIVIDUAL, &val, &err);

	return _scale * val;
	
}




