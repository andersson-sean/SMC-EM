#ifndef IMotionModel_H
#define IMotionModel_H

#include <iostream>

// This is the interface for a motion model:
// 1) It must provide a method for evaluating its transition density.
// 2) It must provide a method for sampling from its transition density.
// 3) It must provide a method for updating its parameters.
class IMotionModel{

public:
	IMotionModel() {}
	virtual ~IMotionModel() {}
	virtual double Evaluate(const double&,const double&) const = 0;
	virtual double Sample(const double&) = 0;
	virtual void UpdateParameters() = 0;
	
	// Overload the output stream operator to allow easy access to model parameters
	friend std::ostream& operator<<(std::ostream& os, const IMotionModel& mm) {
		mm.print(os);
		return os;
	}
	friend std::ostream& operator<<(std::ostream& os, const IMotionModel* mm) {
		mm->print(os);
		return os;
	}
	
private:
	virtual void print(std::ostream&) const = 0;
};

#endif