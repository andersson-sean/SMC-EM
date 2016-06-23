#ifndef PoissonModel_H
#define PoissonModel_H

#include "IObservationModel.h"
#include "../DataManager/DataManager.h"
#include <vector>
#include <functional>
#include <iostream>
#include <boost/function.hpp>

class PoissonModel : public IObservationModel {

private:
	boost::function<double(const std::vector<double>, const std::vector<double>)> _mean;
	double _bias;
	double _scale;
	DataManager* _dm; // reference to the data manager
	
	void print(std::ostream&) const;
	
public:
	PoissonModel () :
		IObservationModel(),
		_bias(0),
		_scale(0) {};
	
	PoissonModel (
		boost::function<double(const std::vector<double>, const std::vector<double>)> mean,
		const double& bias,
		const double& scale) :
		IObservationModel(),
		_mean(mean),
		_bias(bias),
		_scale(scale) {};
	
	double Evaluate(const double&, const std::vector<double>&, const std::vector<double>&);
	void UpdateParameters();
	
	void SetDataManager(DataManager* dm) { _dm = dm; }
};

#endif