#include "TetheredDiffusionModel.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/random/normal_distribution.hpp>

// Evaluates the transition function given the current position and the previous position.
// For tethered diffusion (Ito SDE)
// dX(t) = -A*X(t)dt + sqrt(2D)*W(t),
// this is a normal distribution with and variance specified as
// mean = prevPos * exp(-A*dt)
// variance = (D/A) * (1-exp(-2*dt*A))
// where A is the stiffness parameter, D is the diffusion coefficient, and dt is the
// sampling period.  Note that it is critical for A>0, otherwise the generated position
// will diverge from zero.
double TetheredDiffusionModel::Evaluate
	(const double& currentPosition, const double& previousPosition) const{
	const double stdDev = sqrt( (_diffCoeff / _stiffCoeff) * (1-exp(-2*_dt*_stiffCoeff)) );
	const double mean = previousPosition * exp(-_dt * _stiffCoeff);
	boost::math::normal_distribution<> dist(mean,stdDev);
	return boost::math::pdf(dist,currentPosition);	
}

// Generates a sample from a tethered diffusion given the previous position.
double TetheredDiffusionModel::Sample(const double& previousPosition){
	const double stdDev = sqrt( (_diffCoeff / _stiffCoeff) * (1-exp(-2*_dt*_stiffCoeff)) );
	const double mean = previousPosition * exp(-_dt * _stiffCoeff);
	boost::random::normal_distribution<> dist(mean,stdDev);
	return dist(gen);
}

// Given data from the particle filter/smoother, evaluate the next best diffusion
// coefficient.
void TetheredDiffusionModel::UpdateParameters(){
	switch(_axis){
		case DataManager::X:
			_stiffCoeff = CalculateStiffnessCoefficient(_dm->xParticle);
			_diffCoeff = CalculateDiffusionCoefficient(_dm->xParticle);
			break;
		case DataManager::Y:
			_stiffCoeff = CalculateStiffnessCoefficient(_dm->yParticle);
			_diffCoeff = CalculateDiffusionCoefficient(_dm->yParticle);
			break;
		case DataManager::Z:
			_stiffCoeff = CalculateStiffnessCoefficient(_dm->zParticle);
			_diffCoeff = CalculateDiffusionCoefficient(_dm->zParticle);
			break;
		default:
			break;
	}
}

// Utility function used by UpdateParameters to update the stiffCoeff parameter.
double TetheredDiffusionModel::CalculateStiffnessCoefficient 
	(const std::vector<std::vector<double> >& particle){
  	double num = 0;
  	double den = 0;
  	for (size_t t = 0; t < _dm->NumTimesteps()-1; ++t) {
    	for (size_t j = 0; j < _dm->NumParticles(); ++j) {
	  		for (size_t i = 0; i < _dm->NumParticles(); ++i) {
	      		num += _dm->wJoint[t][i][j] * particle[t+1][j] * particle[t][i];
	      		den += _dm->wJoint[t][i][j] * particle[t][i] * particle[t][i];
	    	}
		}
	}
	return -log(num / den) / _dt;
}

// Utility function used by UpdateParameters to update the diffCoeff parameter.
double TetheredDiffusionModel::CalculateDiffusionCoefficient
	(const std::vector<std::vector<double> >& particle) {
	const double scale = exp(-_dt * _stiffCoeff);
	double tmp = 0;
	for(size_t t = 0; t < _dm->NumTimesteps() - 1; ++t)
		for(size_t j = 0; j < _dm->NumParticles(); ++j)
			for(size_t i = 0; i < _dm->NumParticles(); ++i)
				tmp += _dm->wJoint[t][i][j] * 
					(particle[t+1][j] - scale*particle[t][i])*
						(particle[t+1][j] - scale*particle[t][i]);
	return tmp / ((double)(_dm->NumTimesteps()-1) * 2 * _dt );
}

// Output stream implementation
void TetheredDiffusionModel::print(std::ostream& os) const {
	os << _diffCoeff << " " << _stiffCoeff << " ";
}