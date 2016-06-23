#include "EMManager.h"
#include <algorithm>
#include <functional>
#include <iostream>

void EMManager::Execute() {
	for(unsigned int iter = 0; iter < _numIterations; ++iter){
		
		std::cout << "---------------------------------------------" << std::endl;
		std::cout << "Starting EM iteration #" << iter << std::endl;
		
		// Expectation step:
		std::cout << "Filtering........";
		_smc->Filter();
		std::cout << "complete!" << std::endl;
		
		std::cout << "Smoothing........";
		_smc->Smooth();
		std::cout << "complete!" << std::endl;
		
		// Maximization step:
		std::cout << "Maximizing ICs:" << std::endl;
		std::for_each(_ic.begin(),_ic.end(),std::mem_fun(&IInitialConditionModel::UpdateParameters));
		std::cout << "Maximizing MMs:" << std::endl;
		std::for_each(_mm.begin(),_mm.end(),std::mem_fun(&IMotionModel::UpdateParameters));
		std::cout << "Maximizing OBs:" << std::endl;
		std::for_each(_obsv.begin(),_obsv.end(),std::mem_fun(&IObservationModel::UpdateParameters));
		
		std::cout << "Saving data........";
		// Save results to disk:
		_fm->SaveData(iter);
		std::cout << "complete!" << std::endl;
		std::cout << "---------------------------------------------" << std::endl;
	}
}