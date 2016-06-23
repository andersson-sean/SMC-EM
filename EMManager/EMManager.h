#ifndef EMManager_H
#define EMManager_H

#include <vector>
#include <algorithm>

#include "../SMC/SMC.h"
#include "../InitialConditionModels/IInitialConditionModel.h"
#include "../MotionModels/IMotionModel.h"
#include "../ObservationModels/IObservationModel.h"
#include "../FileManager/FileManager.h"

class EMManager {

private:
	unsigned int _numIterations;
	
	SMC* _smc;
	FileManager* _fm;
	std::vector<IInitialConditionModel*> _ic;
	std::vector<IMotionModel*> _mm;
	std::vector<IObservationModel*> _obsv;

	EMManager();
	EMManager(const EMManager&);
	EMManager& operator=(const EMManager&);
	
public:
	EMManager(unsigned int numIterations) :
		_numIterations(numIterations) {}
		
	void Execute();
	
	void SetFileManager(FileManager* fm){ _fm = fm; }
	void SetSMC(SMC* smc){ _smc = smc; }
	void SetMotionModel(const std::vector<IMotionModel*> mm) { 
		_mm.resize(mm.size());
		std::copy(mm.begin(),mm.end(),_mm.begin());
	}
	void SetInitialConditionModel(const std::vector<IInitialConditionModel*> ic) { 
		_ic.resize(ic.size());
		std::copy(ic.begin(),ic.end(),_ic.begin());
	}
	void SetObservationModel(const std::vector<IObservationModel*> obsv) { 
		_obsv.resize(obsv.size());
		std::copy(obsv.begin(),obsv.end(),_obsv.begin());
	}
};

#endif