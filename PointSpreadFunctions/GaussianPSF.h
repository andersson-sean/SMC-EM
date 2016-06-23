// A functor representing a Gaussian field of the form:
// F(x,y) = exp[-0.5*( (x0-y0)^2 / (2*s0^2) + ... )]

// Note that the vector<double> constructor takes in (s0,s1,s2,...) values and that
// they are automatically squared, reciprocal-ized, and scaled to allow for more efficient
// computation.  

// Standard usage:
// // Specify the widths:
// std::vector<double> widths(3);
// widths[0] = 0.1; widths[1] = 0.1; widths[2] = 0.1;
// // Construct the functor with the chosen widths:
// GaussianPSF psf(widths);
// // Specify the origin of the PSF:
// std::vector<double> origin(3);
// origin[0] = 1; origin[2] = 2; origin[3] = 3;
// // Specify the evaluation point:
// std::vector<double> evalPoint(3);
// evalPoint[0] = 1; evalPoint[2] = 2; evalPoint[3] = 3;
// // Call the functor on the evaluation point:
// double eval = psf(origin,evalPoint);

#ifndef GaussianPSF_H
#define GaussianPSF_H

#include <vector>
#include <algorithm>
#include <functional>

class GaussianPSF {
private:
	std::vector<double> _widthsScaled;

public:
	GaussianPSF() {
		_widthsScaled.push_back(1.0); // default width is one-dimensional
	}
	
	GaussianPSF(const std::vector<double>& widths) {
		// Allocate space:
		_widthsScaled.resize(widths.size());
		// Copy used-defined widths:
		std::copy(
			widths.begin(),widths.end(), // input
			_widthsScaled.begin() // output
		);
		// Squared the widths:
		std::transform( // element-wise square
			_widthsScaled.begin(),_widthsScaled.end(), // input #1
			_widthsScaled.begin(), // input #2
			_widthsScaled.begin(), // output
			std::multiplies<double>() // method
		);
		// Invert the squared widths and multiply by -1/2:
		std::transform( // element-wise reciprocal and multiply by -1/2
			_widthsScaled.begin(),_widthsScaled.end(), // input
			_widthsScaled.begin(), // output
			std::bind1st(std::divides<double>(),-0.5) // method
		);
	}
		
	double operator()(const std::vector<double>&, const std::vector<double>&);
};

#endif