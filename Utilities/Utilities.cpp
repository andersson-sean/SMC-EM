#ifndef DataManager_H
#define DataManager_H

#include <vector>

// The DataManager is the gatekeeper to all the important data within the program.
// Despite this, it's actually a dumb class; you can think of it as a POD structure.
// It is assumed that once the DataManager is constructed, the size of the data it contains
// will not change; thus, it is a semi-static structure.
class DataManager {
private:
	// These should never change after construction:
	double _numParticles;
	double _numTimesteps;
	double _numPixels;
	
	// Make the default constructor private to avoid misuse:
	DataManager();
	
	// Make the copy constructor and copy-assign operators invisible so they're not
	// accidentally called.  You definitely _do not_ want to accidentally copy all of
	// this data during execution.
	DataManager(const DataManager&);
	DataManager& operator=(const DataManager&);
	
public:
	// Some useful typedefs
	enum AxisEnum {X,Y,Z};
  	typedef std::vector<double> Array1D;
  	typedef std::vector<std::vector<double> > Array2D;
  	typedef std::vector<std::vector<std::vector<double> > > Array3D;
	
	// Public Data:
	Array2D xSensor, ySensor, zSensor; // Sensor Positions [time][pixel]
	Array2D xParticle, yParticle, zParticle; // Particle Positions [time][particle]
	Array2D wFilter, wSmooth; // Weights [time][particle]
	Array2D Data; // Measurements [time][pixel]
  	Array3D wJoint; // Joint Weights [time][particle][particle]
  	
  	// One-and-done constructor:
	DataManager(size_t numParticles, size_t numTimesteps, size_t numPixels) :
		_numParticles(numParticles),
		_numTimesteps(numTimesteps),
		_numPixels(numPixels) {	
			xSensor.resize(numTimesteps, std::vector<double>(numPixels));;
			ySensor.resize(numTimesteps, std::vector<double>(numPixels));;
			zSensor.resize(numTimesteps, std::vector<double>(numPixels));;
			xParticle.resize(numTimesteps, std::vector<double>(numParticles));
			yParticle.resize(numTimesteps, std::vector<double>(numParticles));
			zParticle.resize(numTimesteps, std::vector<double>(numParticles));
			wFilter.resize(numTimesteps, std::vector<double>(numParticles));
			wSmooth.resize(numTimesteps, std::vector<double>(numParticles));
			Data.resize(numTimesteps, std::vector<double>(numPixels));
			wJoint.resize(numTimesteps-1,std::vector<std::vector<double> >(numParticles,std::vector<double>(numParticles)));
	}
	
	// Public accessors:
	size_t NumTimesteps() { return _numTimesteps;}
	size_t NumParticles() { return _numParticles;}
	size_t NumPixels() {return _numPixels;}

};

#endif