#ifndef IInitialConditionModel_H
#define IInitialConditionModel_H

#include <iostream>

// This is the interface for a motion model:
// 1) It must provide a method for sampling from its transition density.
// 2) It must provide a method for updating its parameters.
class IInitialConditionModel{

public:
	IInitialConditionModel() {}
	virtual ~IInitialConditionModel() {}
	virtual double Sample() = 0;
	virtual void UpdateParameters() = 0;
	
	// Overload the output stream operator to allow easy access to model parameters
	friend std::ostream& operator<<(std::ostream& os, const IInitialConditionModel& ic) {
		ic.print(os);
		return os;
	}
	friend std::ostream& operator<<(std::ostream& os, const IInitialConditionModel* ic){
		ic->print(os);
		return os;
	}
	
private:
	virtual void print(std::ostream&) const = 0;
};

#endif