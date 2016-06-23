#include "PoissonModel.h"
#include <boost/math/distributions/poisson.hpp>

// Evaluates the transition function given the current position and the previous position.
// For free diffusion, this is a normal distribution with mean set to "previousPosition"
// and standard deviation set to sqrt(2D*dt).
double PoissonModel::Evaluate
	(const double& measurement,
		const std::vector<double>& origin,
		const std::vector<double>& evalPoint){
	double mean = _scale * _mean(origin,evalPoint) + _bias;
	boost::math::poisson_distribution<> dist(mean);
	return boost::math::pdf(dist,measurement);
}

// Given data from the particle filter/smoother, evaluate the next best value for 
// the scale.
void PoissonModel::UpdateParameters(){

	// Calculate the peak intensity value
    // This must be done numerically, so we use the Tangent Hyperbolas method,
    // which is a variant of Newton's method.
    const double numIter = 100; // Maximum number of iterations
    const double tol = 0.001; // Tolerance to achieve before terminating
    double tmp = _scale;
    for (unsigned int iter = 0; iter < numIter; iter++)
    {
    	double func0 = 0; // first derivative
    	double func1 = 0; // second derivative
    	double func2 = 0; // third derivative
    	for (size_t t = 0; t < _dm->NumTimesteps(); ++t) {
    		for (size_t part = 0; part < _dm->NumParticles(); ++part) {
    			std::vector<double> origin(3);
    			origin[0] = _dm->xParticle[t][part];
    			origin[1] = _dm->yParticle[t][part];
    			origin[2] = _dm->zParticle[t][part];
    			for(size_t pix = 0; pix < _dm->NumPixels(); ++pix) {
    				
    				std::vector<double> evalPoint(3);
					evalPoint[0] = _dm->xSensor[t][pix];
					evalPoint[1] = _dm->ySensor[t][pix];
					evalPoint[2] = _dm->zSensor[t][pix];					
    				const double lambda = _mean(origin,evalPoint);
					
					//std::cout << evalPoint[0] << " " << evalPoint[1] << " " << evalPoint[2] << std::endl;
					//std::cout << origin[0] << " " << origin[1] << " " << origin[2] << std::endl; 
					//std::cout << lambda << std::endl;
					//int i = 0;
					//std::cin >> i;
					
    				double den = (tmp * lambda + _bias);
					func0 += _dm->wSmooth[t][part] * lambda * (_dm->Data[t][pix] / den - 1);
					func1 -= _dm->wSmooth[t][part] * lambda * lambda * (_dm->Data[t][pix] / (den * den));
					func2 += _dm->wSmooth[t][part] * lambda * lambda * lambda * (2 * _dm->Data[t][pix] / (den * den * den));
				}
    		}
    	}
    	// double decrement = func1; // Gradient descent
    	//double decrement = func1 / func2; // Newton's method
    	double decrement = 2 * func0 * func1 / ((2 * func1 * func1) - (func0 * func2)); // Tangent hyperbolas method
    	tmp -= decrement;
    	//std::cout << tmp << std::endl;
    	if(abs(decrement) < tol)
    		break;
    }
    _scale = tmp;
    std::cout << "PSF Gain: " << _scale << std::endl;
}

// Output stream implementation
void PoissonModel::print(std::ostream& os) const {
	os << _scale << " ";
}