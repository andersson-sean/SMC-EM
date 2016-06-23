#include <string>
#include <iostream>
#include <vector> 

#include <boost/function.hpp>

#include "DataManager/DataManager.h"

#include "InitialConditionModels/NormalInitialCondition.h"
#include "InitialConditionModels/ConstantInitialCondition.h"

#include "MotionModels/FreeDiffusionModel.h"
#include "MotionModels/DirectedDiffusionModel.h"
#include "MotionModels/NoMotionModel.h"

#include "ObservationModels/PoissonModel.h"
#include "PointSpreadFunctions/GaussianPSF.h"
#include "PointSpreadFunctions/PixelatedGaussianPSF.h"

#include "SMC/SMC.h"
#include "EMManager/EMManager.h"
#include "FileManager/FileManager.h"


int main(int argc, char* argv[]) {

	///////////////////////////////////////////////////////////////
	// BEGIN User-modifiable code
	
	
	///////////////////////////
	// Inference parameters:
	const size_t numEMIterations = 100; // number of EM iterations
	const size_t numParticles = 512; // number of SMC estimates
	const size_t numTimesteps = 200; // number of images
	const size_t numPixels = 49; // number of pixels per image
	const double imagingPeriod =  0.033; // image acquisition period [s]
	const double pixelSize = 0.204; // size (e.g., length or width) of a square pixel [um]

	
	/////////////////////////////////////
	// Create initial condition models:
	const double xInitCondMean = 1.66; // First guess at X's mean initial condition
	const double xInitCondStdDev = 0.1; // First guess at X's initial condition std. dev.
	const unsigned int xInitCondSeed = 1; // X's initial condition seed
	
	const double yInitCondMean = 1.93; // First guess at Y's mean initial condition
	const double yInitCondStdDev = 0.1; // First guess at Y's initial condition std. dev.
	const unsigned int yInitCondSeed = 2; // Y's initial condition seed
	
	const double zInitCondMean = 0; // First guess at Z's mean initial condition
	const double zInitCondStdDev = 0.1; // First guess at Z's initial condition std. dev.
	const unsigned int zInitCondSeed = 3; // Z's initial condition seed
	
	NormalInitialCondition* icX = 
		new NormalInitialCondition(xInitCondMean,xInitCondStdDev,DataManager::X,xInitCondSeed);
	NormalInitialCondition* icY = 
		new NormalInitialCondition(yInitCondMean,yInitCondStdDev,DataManager::Y,yInitCondSeed);
	NormalInitialCondition* icZ =
		new NormalInitialCondition(zInitCondMean,zInitCondStdDev,DataManager::Z,zInitCondSeed);


	///////////////////////////
	// Create motion models:
	const double xDiffCoeff = 0.001; // First guess at X's diffusion coefficient
	const double xSpeed = 0.1; // First guess at X's speed
	const unsigned int xDiffCoeffSeed = 1; // X's motion model seed
	
	const double yDiffCoeff = 0.001; // First guess at Y's diffusion coefficient
	const double ySpeed = 0.1; // First guess at Y's speed
	const unsigned int yDiffCoeffSeed = 2; // Y's motion model seed
	
	const double zDiffCoeff = 0.0001; // First guess at Z's diffusion coefficient
	const unsigned int zDiffCoeffSeed = 3; // Z's motion model seed
	
	DirectedDiffusionModel* mmX = 
		new DirectedDiffusionModel(imagingPeriod,xDiffCoeff,xSpeed,DataManager::X,xDiffCoeffSeed);
	DirectedDiffusionModel* mmY = 
		new DirectedDiffusionModel(imagingPeriod,yDiffCoeff,ySpeed,DataManager::Y,yDiffCoeffSeed);
	FreeDiffusionModel* mmZ = 
		new FreeDiffusionModel(imagingPeriod,zDiffCoeff,DataManager::Z,zDiffCoeffSeed);
	

	///////////////////////////
	// Construct the PSF
	std::vector<double> widths(3); // allocate space for Gaussian widths
	widths[0] = 0.2625; widths[1] = 0.2625; widths[2] = 0.394; // specify the widths- x,y as 0.6*lambda/NA, z as 1.5 of that.
	PixelatedGaussianPSF psf(widths,pixelSize); // create Gaussian functor
	boost::function<
		double(const std::vector<double>,const std::vector<double>)> 
			psfFunctor = psf; // create functor wrapper
	
	// Create observation model:
	const double bias = 1646;
	const double gain = 1489;
	PoissonModel* obsvModel = new PoissonModel(psfFunctor,bias,gain);

	// END User-modifiable code
	////////////////////////////////////////////////////////////////////

	// Concatenate initial condition, motion, and observation models
	std::vector<IInitialConditionModel*> icVec(3);
	std::vector<IMotionModel*> mmVec(3);
	std::vector<IObservationModel*> obsvVec(1);
	icVec[0] = icX; icVec[1] = icY; icVec[2] = icZ;
	mmVec[0] = mmX; mmVec[1] = mmY; mmVec[2] = mmZ;
	obsvVec[0] = obsvModel;

	// Create the data manager:
	DataManager* dm = new DataManager(numParticles,numTimesteps,numPixels);
	
	// Inform the models about the data manager
	icX->SetDataManager(dm);
	icY->SetDataManager(dm);
	icZ->SetDataManager(dm);
	mmX->SetDataManager(dm);
	mmY->SetDataManager(dm);
	mmZ->SetDataManager(dm);
	obsvModel->SetDataManager(dm);
	
	
	// Create the SMC object:
	SMC* smc = new SMC();
	smc->SetDataManager(dm);
	smc->SetInitialConditionModel(icVec);
	smc->SetMotionModel(mmVec);
	smc->SetObservationModel(obsvVec);

	// Parse the input from the user.  Should be in the format
	// "./program_name 'data_path' 'results_path' 'sensorpath/name"
	std::string DataPath;
	std::string ResultsPath;
	std::string SensorFilePathAndName;
	if(argc != 4) {
		std::cout << "Error: program must be called as" << std::endl;
		std::cout << "./program_name 'data_path' 'output_path' 'sensor_file_path/name" << std::endl;
		std::cout << "Example:" << std::endl;
		std::cout << "./main '/project/smcfilt/ConfocalData/1/' '/project/smcfilt/ConfocalData/1/Results/' /project/smcfilt/ConfocalData/1/Sensor/sensor" << std::endl; 
		std::cout << "Aborting operation." << std::endl;
		return 0;
	} else {
		DataPath.assign(argv[1]);
		ResultsPath.assign(argv[2]);
		SensorFilePathAndName.assign(argv[3]);
	}
	
	// Create the FileManager:
	FileManager* fm = new FileManager(SensorFilePathAndName,DataPath,ResultsPath);
	fm->SetDataManager(dm);
	fm->SetInitialConditionModel(icVec);
	fm->SetMotionModel(mmVec);
	fm->SetObservationModel(obsvVec);
	fm->LoadData();

	// Create the EM manager:
	EMManager* em = new EMManager(numEMIterations);
	em->SetSMC(smc);
	em->SetFileManager(fm);
	em->SetInitialConditionModel(icVec);
	em->SetMotionModel(mmVec);
	em->SetObservationModel(obsvVec);

	// Run the main algorithm
	em->Execute();

	// Clean up the program
	delete dm;
	delete icX;
	delete icY;
	delete icZ;
	delete mmX;
	delete mmY;
	delete mmZ;
	delete obsvModel;
	delete smc;
	delete fm;
	delete em;

}