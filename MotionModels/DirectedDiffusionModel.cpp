#include "DirectedDiffusionModel.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/random/normal_distribution.hpp>

// Evaluates the transition function given the current position and the previous position.
// For directed diffusion according to the Ito SDE
// dX(t) = V*dt + sqrt(2*D)*dW,
// this is a normal distribution with parameters
// mean = prevPos + V*dt,
// variance = 2*D*dt,
// where V is the speed, D is the diffusion coefficient, and dt is the sampling period.
double DirectedDiffusionModel::Evaluate
	(const double& currentPosition, const double& previousPosition) const{
	const double stdDev = sqrt(2*_diffCoeff*_dt);
	const double mean = previousPosition + (_speed * _dt);
	boost::math::normal_distribution<> dist(mean,stdDev);
	return boost::math::pdf(dist,currentPosition);	
}

// Generates a sample from a directed diffusion given the previous position.
double DirectedDiffusionModel::Sample(const double& previousPosition){
	const double stdDev = sqrt(2*_diffCoeff*_dt);
	const double mean = previousPosition + (_speed * _dt);
	boost::random::normal_distribution<> dist(mean,stdDev);
	return dist(gen);
}

// Given data from the particle filter/smoother, evaluate the next parameters
void DirectedDiffusionModel::UpdateParameters(){
	switch(_axis){
		case DataManager::X:
			_speed = CalculateSpeed(_dm->xParticle);
			_diffCoeff = CalculateDiffusionCoefficient(_dm->xParticle);
			break;
		case DataManager::Y:
			_speed = CalculateSpeed(_dm->yParticle);
			_diffCoeff = CalculateDiffusionCoefficient(_dm->yParticle);
			break;
		case DataManager::Z:
			_speed = CalculateSpeed(_dm->zParticle);
			_diffCoeff = CalculateDiffusionCoefficient(_dm->zParticle);
			break;
		default:
			break;
	}

}

// Utility function used by UpdateParameters
double DirectedDiffusionModel::CalculateSpeed
	(const std::vector<std::vector<double> >& particle) {
	double tmp = 0;
	for(size_t t = 0; t < _dm->NumTimesteps() - 1; ++t)
		for(size_t j = 0; j < _dm->NumParticles(); ++j)
			for(size_t i = 0; i < _dm->NumParticles(); ++i)
				tmp += _dm->wJoint[t][i][j] * (particle[t+1][j] - particle[t][i]);
	return tmp / ((double)(_dm->NumTimesteps()-1) * _dt );
}

// Utility function used by UpdateParameters
double DirectedDiffusionModel::CalculateDiffusionCoefficient
	(const std::vector<std::vector<double> >& particle) {
	double bias = _speed * _dt;
	double tmp = 0;
	for(size_t t = 0; t < _dm->NumTimesteps() - 1; ++t)
		for(size_t j = 0; j < _dm->NumParticles(); ++j)
			for(size_t i = 0; i < _dm->NumParticles(); ++i)
				tmp += _dm->wJoint[t][i][j] * 
					(particle[t+1][j] - particle[t][i] - bias)*
						(particle[t+1][j] - particle[t][i] - bias);
	return tmp / ((double)(_dm->NumTimesteps()-1) * 2 * _dt );
}

// Output stream implementation
void DirectedDiffusionModel::print(std::ostream& os) const {
	os << _diffCoeff << " " << _speed << " ";
}