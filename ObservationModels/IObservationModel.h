#ifndef IObservationModel_H
#define IObservationModel_H

#include <vector>
#include <iostream>

// This is the interface for a observation model:
// 1) It must provide a method for evaluating its density.
// 2) It must provide a method for updating its parameters.
class IObservationModel{

public:
	IObservationModel() {}
	virtual ~IObservationModel() {}
	virtual double Evaluate(const double&, const std::vector<double>&,const std::vector<double>&) = 0;
	virtual void UpdateParameters() = 0;

	// Overload the output stream operator to allow easy access to model parameters
	friend std::ostream& operator<<(std::ostream& os, const IObservationModel& obs) {
		obs.print(os);
		return os;
	}
	friend std::ostream& operator<<(std::ostream& os, const IObservationModel* obs) {
		obs->print(os);
		return os;
	}
private:
	virtual void print(std::ostream&) const = 0;

};

#endif