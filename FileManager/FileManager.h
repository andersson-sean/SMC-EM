#ifndef FileManager_H
#define FileManager_H


#include <string>
#include "../DataManager/DataManager.h"
#include "../MotionModels/IMotionModel.h"
#include "../ObservationModels/IObservationModel.h"
#include "../InitialConditionModels/IInitialConditionModel.h"

class FileManager {

private:	

	DataManager* _dm;
	
	std::vector<IInitialConditionModel*> _ic;
	std::vector<IMotionModel*> _mm;
	std::vector<IObservationModel*> _obsv;
	
	std::string _sensorFilename;
	std::string _imgDir;
	std::string _resultsDir;

	FileManager();
	FileManager(const FileManager&);
	FileManager& operator=(const FileManager&);
	
	void VerifyPaths();
	
public:
	FileManager(
		const std::string& sensorFilename,
		const std::string& imgDir,
		const std::string& resultsDir) : 
			_sensorFilename(sensorFilename),
			_imgDir(imgDir),
			_resultsDir(resultsDir) {
			
			VerifyPaths();
		};
			
	void LoadData();
	void SaveData(unsigned int);
	
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