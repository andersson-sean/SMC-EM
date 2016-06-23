#include "ConstantInitialCondition.h"


double ConstantInitialCondition::Sample(){
	return _bias;
}


void ConstantInitialCondition::UpdateParameters(){
}

// Output stream implementation
void ConstantInitialCondition::print(std::ostream& os) const {
}

