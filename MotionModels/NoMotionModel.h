#ifndef NoMotionModel_H
#define NoMotionModel_H

#include "IMotionModel.h"
#include "../DataManager/DataManager.h"
#include <iostream>

class NoMotionModel : public IMotionModel {

private:
	double _bias;
	DataManager::AxisEnum _axis; // axis to which this model is associated
	DataManager* _dm; // reference to the data manager
	
	void print(std::ostream&) const;
public:
	NoMotionModel (double bias) :
		IMotionModel(),
		_bias(bias),
		_axis(DataManager::X) {} // default to the X axis
	
	double Evaluate(const double&, const double&) const;
	double Sample(const double&);
	void UpdateParameters();

	void SetDataManager(DataManager* dm) { _dm = dm; }

};

#endif