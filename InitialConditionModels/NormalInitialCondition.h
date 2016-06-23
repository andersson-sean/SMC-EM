#ifndef NormalInitialCondition_H
#define NormalInitialCondition_H

#include "IInitialConditionModel.h"
#include "../DataManager/DataManager.h"
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <iostream>

class NormalInitialCondition : public IInitialConditionModel {

private:
	double _mean;
	double _stdDev;
	DataManager::AxisEnum _axis; // axis to which this model is associated
	DataManager* _dm; // reference to the data manager
	boost::random::mt19937 gen; // Mersenne Twister random number generator
	
	double CalculateMean (const std::vector<std::vector<double> >&);
	double CalculateStandardDeviation (const std::vector<std::vector<double> >&);
	
	void print(std::ostream&) const;
public:
	NormalInitialCondition () :
		IInitialConditionModel(),
		_mean(0),
		_stdDev(1),
		_axis(DataManager::X) // default to the X axis
		 {
			gen.seed(0); // by default set seed to zero
		};
		
	NormalInitialCondition
		(const double& mean, 
			const double& stdDev, 
			const DataManager::AxisEnum& axis,
			const unsigned int& seed = 0) : 
		IInitialConditionModel(),
		_mean(mean),
		_stdDev(stdDev),
		_axis(axis) {
			gen.seed(seed); // set the RNG seed
		};
	
	double Sample();
	void UpdateParameters();
	
	void SetDataManager(DataManager* dm) { _dm = dm; }
};

#endif