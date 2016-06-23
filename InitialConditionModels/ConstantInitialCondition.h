#ifndef ConstantInitialCondition_H
#define ConstantInitialCondition_H

#include "IInitialConditionModel.h"
#include "DataManager/DataManager.h"
#include <iostream>

class ConstantInitialCondition : public IInitialConditionModel {

private:
	double _bias;
	DataManager* _dm;
	
	void print(std::ostream&) const;
public:
	ConstantInitialCondition (double bias) :
		IInitialConditionModel(),
		_bias(bias) {}
		
	double Sample();
	void UpdateParameters();
	
	void SetDataManager(DataManager* dm) { _dm = dm; }
};

#endif