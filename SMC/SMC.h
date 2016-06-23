#ifndef SMC_H
#define SMC_H

#include "../InitialConditionModels/IInitialConditionModel.h"
#include "../MotionModels/IMotionModel.h"
#include "../ObservationModels/IObservationModel.h"
#include "../DataManager/DataManager.h"
#include <vector>

class SMC {
private:
	
	std::vector<IInitialConditionModel*> _ic;
	std::vector<IMotionModel*> _mm;
	std::vector<IObservationModel*> _obsv;
	DataManager* _dm;

	bool NormalizeLog(std::vector<double>&);
	void Equalize(std::vector<double>&);
	void GenerateResamplingIndices(const std::vector<double>&, std::vector<double>&);
	void Resample(std::vector<double>&, const std::vector<double>&);
	
public:

	void Filter();
	void Smooth();
	SMC() {}

	// Mutators
	void SetDataManager(DataManager* dm) { _dm = dm; }
	void SetInitialConditionModel(const std::vector<IInitialConditionModel*>& ic) {
		_ic.resize(ic.size());
		std::copy(ic.begin(), ic.end(), _ic.begin());
	}
	void SetMotionModel(const std::vector<IMotionModel*>& mm) {
		_mm.resize(mm.size());
		std::copy(mm.begin(), mm.end(), _mm.begin());
	}
	void SetObservationModel(const std::vector<IObservationModel*>& obsv) {
		_obsv.resize(obsv.size());
		std::copy(obsv.begin(), obsv.end(), _obsv.begin());
	}

};

#endif