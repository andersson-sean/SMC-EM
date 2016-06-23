#include "NormalInitialCondition.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/random/normal_distribution.hpp>
#include <cmath>
#include <iostream>
// Generates a sample from a normally distributed initial condition model.
double NormalInitialCondition::Sample(){
	boost::random::normal_distribution<> dist(_mean,_stdDev);
	return dist(gen);
}

// Given data from the particle filter/smoother, evaluate the next mean and standard
// deviation.  This method just evaluates what axis it is and then passes the relevant
// data to the other methods.
// Important Note: The order in which you calculate the mean and standard deviations
// is critical.  The mean _must_ be updated first as it is used in the evaluation of the
// standard deviation.

void NormalInitialCondition::UpdateParameters(){
	switch(_axis){
		case DataManager::X:
			_mean = CalculateMean(_dm->xParticle);
			_stdDev = CalculateStandardDeviation(_dm->xParticle);
			std::cout << "X Init Cond ~ N(" << _mean << " , " << _stdDev << ")" << std::endl;
			break;
		case DataManager::Y:
			_mean = CalculateMean(_dm->yParticle);
			_stdDev = CalculateStandardDeviation(_dm->yParticle);
			std::cout << "Y Init Cond ~ N(" << _mean << " , " << _stdDev << ")" << std::endl;
			break;
		case DataManager::Z:
			_mean = CalculateMean(_dm->zParticle);
			_stdDev = CalculateStandardDeviation(_dm->zParticle);
			std::cout << "Z Init Cond ~ N(" << _mean << " , " << _stdDev << ")" << std::endl;
			break;
		default:
			break;
	}	
}

// Calculates the mean for the initial condition for a given axis
// given the data collected with the Data Manager.
double NormalInitialCondition::CalculateMean
	(const std::vector<std::vector<double> >& particle) {
	double tmp = 0;
	for(size_t i = 0; i < _dm->NumParticles(); ++i) {
		//std::cout << _dm->wSmooth[0][i] << " " << particle[0][i] << std::endl;
		tmp += (_dm->wSmooth[0][i] * particle[0][i]);
	}
	return tmp;
	
}

// Calculates the standard deviation for the initial condition for a given axis
// given the data collected with the Data Manager.
double NormalInitialCondition::CalculateStandardDeviation
	(const std::vector<std::vector<double> >& particle) {
	
	double tmp = 0;
	for(size_t i = 0; i < _dm->NumParticles(); ++i)
		tmp += _dm->wSmooth[0][i]*
			(particle[0][i]-_mean)*
				(particle[0][i]-_mean);
	return sqrt(tmp);
}

// Output stream implementation
void NormalInitialCondition::print(std::ostream& os) const {
	os << _mean << " " << _stdDev << " ";
}

