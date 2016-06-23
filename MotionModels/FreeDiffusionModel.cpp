#include "FreeDiffusionModel.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/random/normal_distribution.hpp>

// Evaluates the transition function given the current position and the previous position.
// For free diffusion, this is a normal distribution with mean set to "previousPosition"
// and standard deviation set to sqrt(2D*dt).
double FreeDiffusionModel::Evaluate
	(const double& currentPosition, const double& previousPosition) const{
	const double stdDev = sqrt(2*_diffCoeff*_dt);
	const double mean = previousPosition;
	boost::math::normal_distribution<> dist(mean,stdDev);
	return boost::math::pdf(dist,currentPosition);	
}

// Generates a sample from a free diffusion given the previous position.
double FreeDiffusionModel::Sample(const double& previousPosition){
	const double stdDev = sqrt(2*_diffCoeff*_dt);
	const double mean = previousPosition;
	boost::random::normal_distribution<> dist(mean,stdDev);
	return dist(gen);
}

// Given data from the particle filter/smoother, evaluate the next best parameters
void FreeDiffusionModel::UpdateParameters(){
	switch(_axis){
		case DataManager::X:
			_diffCoeff = CalculateDiffusionCoefficient(_dm->xParticle);
			std::cout << "X Diffusion Coefficient: " << _diffCoeff << std::endl;
			break;
		case DataManager::Y:
			_diffCoeff = CalculateDiffusionCoefficient(_dm->yParticle);
			std::cout << "Y Diffusion Coefficient: " << _diffCoeff << std::endl;
			break;
		case DataManager::Z:
			_diffCoeff = CalculateDiffusionCoefficient(_dm->zParticle);
			std::cout << "Z Diffusion Coefficient: " << _diffCoeff << std::endl;
			break;
		default:
			break;
	}

}

// Utility function used by UpdateParameters
double FreeDiffusionModel::CalculateDiffusionCoefficient
	(const std::vector<std::vector<double> >& particle) {
	double tmp = 0;
	for(size_t t = 0; t < _dm->NumTimesteps() - 1; ++t)
		for(size_t j = 0; j < _dm->NumParticles(); ++j)
			for(size_t i = 0; i < _dm->NumParticles(); ++i)
				tmp += _dm->wJoint[t][i][j] * 
					(particle[t+1][j] - particle[t][i])*
						(particle[t+1][j] - particle[t][i]);
	return tmp / ((double)(_dm->NumTimesteps()-1) * 2 * _dt );
}

// Output stream implementation
void FreeDiffusionModel::print(std::ostream& os) const {
	os << _diffCoeff << " ";
}