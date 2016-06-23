#include "NoMotionModel.h"

// Evaluates the transition function given the current position and the previous position.
// For no motion, the probability density is just a delta function.
double NoMotionModel::Evaluate
	(const double& currentPosition, const double& previousPosition) const{
	if(previousPosition == currentPosition)
		return 1.0;	
	else
		return 0.0;
}

// Just returns the bias.
double NoMotionModel::Sample(const double& previousPosition){
	return _bias;
}

// Update the parameters.
void NoMotionModel::UpdateParameters(){
	// No parameters to update.
}

// Output stream implementation
void NoMotionModel::print(std::ostream& os) const {
	// No parameters to print out
}