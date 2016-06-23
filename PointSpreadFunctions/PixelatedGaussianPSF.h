#ifndef PixelatedGaussianPSF_H
#define PixelatedGaussianPSF_H

#include <vector>
#include <algorithm>
#include "GaussianPSF.h"

namespace PixelatedGaussianPSF_cubature {
	
	// An adapter structure to mesh the PixelatedGaussianPSF object
	// with the cubature library.
	struct AdapterStruct {
		GaussianPSF* psf;
		std::vector<double> particle;
		std::vector<double> sensor;
	};
	
	// A method to be called by Cubature for use with the PixelatedGaussianPSF class.
	int AdapterMethod(unsigned, const double*, void*, unsigned, double*);	
}

class PixelatedGaussianPSF {
private:
	double _pixelSize;
	GaussianPSF _gaussPSF;
	double _scale;
	
	PixelatedGaussianPSF() {} // = delete;
	
public:
	
	PixelatedGaussianPSF(const std::vector<double>& widths, double pixelSize) :
		_pixelSize(pixelSize),
		_gaussPSF(widths) { 
			// Set scale to unity, calculate the PSF value at the origin,
			// and then set the scale as the inverse of that;
			// This will make the PSF value at the origin be unity.
			_scale = 1.0;
			std::vector<double> tmp(3,0);
			_scale = 1.0 / (*this)(tmp,tmp);
		}
		
	double operator()(const std::vector<double>&, const std::vector<double>&);
};


#endif