#ifndef NonparaxialDebyeWidefieldPSF_H
#define NonparaxialDebyeWidefieldPSF_H

#include <vector>
#include <cmath>
#include <iostream>

class NonparaxialDebyeWidefieldPSF;

// These are adapter methods and structures that allow for use with cubature:
namespace NonparaxialDebyeWidefieldPSF_cubature{
	struct AdapterStruct {
		NonparaxialDebyeWidefieldPSF* ptr;
		std::vector<double> input;
	};
	int AdapterMethod_realPart(unsigned, const double*, void*, unsigned, double*);
	int AdapterMethod_imagPart(unsigned, const double*, void*, unsigned, double*);
}

class NonparaxialDebyeWidefieldPSF {
private:
	// Parameters describing the PSF
	double _NA; // numerical aperture
	double _n; // index of refraction
	double _lambda; // wavelength of emission
	double _kappa; // wavenumber of emission
	double _alpha; // maximum semiangle of light collected by objective
	double _scale; // a value set so that the origin is unity
	
	// Hide the default constructor:
	NonparaxialDebyeWidefieldPSF() {}
	
	// Methods for calculating the real and imaginary parts:
	double realPart(const std::vector<double>&, const double&); 
	double imagPart(const std::vector<double>&, const double&);
	
	// Friend declarations for real and imaginary parts for cubature adapters: 
	friend int NonparaxialDebyeWidefieldPSF_cubature::AdapterMethod_realPart(unsigned, const double*, void*, unsigned, double*);
	friend int NonparaxialDebyeWidefieldPSF_cubature::AdapterMethod_imagPart(unsigned, const double*, void*, unsigned, double*);
	
public:
	NonparaxialDebyeWidefieldPSF(
		const double& NA,
		const double& n,
		const double& lambda) :
			_NA(NA),
			_n(n),
			_lambda(lambda) {
		_kappa = n * (2 * M_PI / lambda); // wavenumber of light
		_alpha = asin(NA / n); // maximum semiangle of light into objective
		_scale = 1; 
		
		// set scale so that origin is unity
		std::vector<double> tmp(3,0);
		_scale = 1.0 / (*this)(tmp,tmp);
	}
		
	double operator()(const std::vector<double>&, const std::vector<double>&);
};

#endif