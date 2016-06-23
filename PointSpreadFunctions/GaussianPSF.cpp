#include "GaussianPSF.h"
#include <cmath>
#include <numeric>


// Evaluate 
double GaussianPSF::operator()
	(const std::vector<double>& origin, const std::vector<double>& sensor) {
	std::vector<double> tmp(origin.size()); // scratch
	std::transform( // element-wise subtract
		origin.begin(), origin.end(), // input #1
		sensor.begin(), // input #2
		tmp.begin(), // output
		std::minus<double>() // method
	);
	std::transform( // element-wise square
		tmp.begin(), tmp.end(), // input #1
		tmp.begin(), // input#2
		tmp.begin(), // output
		std::multiplies<double>() // method
	);
	std::transform( // element-wise multiply by _widthsScaled
		tmp.begin(), tmp.end(), // input #1
		_widthsScaled.begin(), // input #2
		tmp.begin(), // output
		std::multiplies<double>() // method
	);
	double val = std::accumulate(
		tmp.begin(), tmp.end(), // input
		0.0 // initial condition
	);
	val = exp(val);
	return (val > 0)? val : 0.0;
}
