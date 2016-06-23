#ifndef ConfinedDiffusionModel_H
#define ConfinedDiffusionModel_H

#include "IMotionModel.h"
#include "../DataManager/DataManager.h"
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <iostream>

// Implements the confined diffusion motion model for a Brownian motion with diffusion
// coefficient "diffCoeff" confined to the interval 
// ["center" - "length"/2.0, "center" + "length"/2.0]
class ConfinedDiffusionModel : public IMotionModel {

private:
	double _diffCoeff; // diffusion coefficient
	double _length; // length of confinement channel
	double _center; // center position of confinement channel
	double _dt; // sampling interval
	
	DataManager::AxisEnum _axis; // axis to which this model is associated
	DataManager* _dm; // reference to the data manager
	boost::random::mt19937 gen; // Mersenne Twister random number generator
	double CalculateLength(const std::vector<std::vector<double> >&);
	double CalculateDiffusionCoefficient(const std::vector<std::vector<double> >&);
	
	void print(std::ostream&) const;
public:
	ConfinedDiffusionModel () :
		IMotionModel(),
		_diffCoeff(0),
		_length(1.0),
		_center(0.0),
		_dt(0),
		_axis(DataManager::X) // default to the X axis
		 {
			gen.seed(0); // by default set seed to zero
		};
		
	ConfinedDiffusionModel
		(const double& dt, 
			const double& diffCoeff, 
			const double& length,
			const double& center,
			const DataManager::AxisEnum& axis,
			const unsigned int& seed = 0) : 
		IMotionModel(),
		_diffCoeff(diffCoeff),
		_length(length),
		_center(center),
		_dt(dt),
		_axis(axis) {
			gen.seed(seed); // set the RNG seed
		};
	
	double Evaluate(const double&, const double&) const;
	double Evaluate(const double&, const double&, const double&, const double&) const;
	double EvaluateDerivativeWRTDiffCoeff(const double&, const double&, const double&, const double&) const;
	double EvaluateDerivativeWRTPos(const double&, const double&) const;
	double EvaluateSecondDerivativeWRTPos(const double&, const double&) const;
	double EvaluateThirdDerivativeWRTPos(const double&, const double&) const;
	double Sample(const double&);
	void UpdateParameters();
	
	void SetDataManager(DataManager* dm) { _dm = dm; }

};

#endif