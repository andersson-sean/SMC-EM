#include "FileManager.h"
#include <fstream>
#include <cstring>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <exception>
#include <boost/filesystem.hpp>
#include <tiffio.h>
#include <algorithm>


// Loads the sensor data and the corresponding images
void FileManager::LoadData(){

	// Load the sensor data
  	std::ifstream file(_sensorFilename.c_str());
  	if(file.is_open()) {
      	for (size_t t = 0; t < _dm->NumTimesteps(); ++t) {
      		for(size_t pix = 0; pix < _dm->NumPixels(); ++pix) {
      			if(!file.eof()){
      				file >> _dm->xSensor[t][pix];
					file >> _dm->ySensor[t][pix];
					file >> _dm->zSensor[t][pix];
					//std::cout << "t = " << t << ": " << _dm->xSensor[t][pix] << " " << _dm->ySensor[t][pix] << " " << _dm->zSensor[t][pix] << std::endl;
      			}
      			else {
      				std::cout << "Sensor file length mismatch." << std::endl;
      				throw std::exception();
				}
			}	
		}
    }
    //std::cin.get();
    file.close();  
  
	// Load the images
	for (size_t t = 0; t < _dm->NumTimesteps(); ++t) {
		// Format the image number with leading zeros.
		// Total number of digits should be five.
		std::string imgNumber = boost::lexical_cast<std::string>(t+1);
		imgNumber = std::string(5-imgNumber.length(),'0').append(imgNumber);

		// Open up the TIFF
		std::string tmpString = _imgDir + imgNumber + ".tif";
		char fileName[1024];
		strncpy(fileName,tmpString.c_str(),sizeof(fileName));
		fileName[sizeof(fileName)-1] = 0;
		TIFF* tif = TIFFOpen(fileName,"r");
      
		// If the TIFF exists...
		if(tif) {
			// Query the image size
			uint32 imageLength, imageWidth;
			TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
			TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
	  
	  		// Check to see if the image is the correct size
	  		if(imageLength*imageWidth != _dm->NumPixels()) {
	  			std::cout << "Size of image at time t = " << t << " is not the correct size." << std::endl;
	  			throw std::exception();
	  		}
	  
			// Allocate space for the scanline
			tdata_t buffer = _TIFFmalloc(TIFFScanlineSize(tif));
			for (uint32 row = 0; row < imageLength; row++) {
				// Grab the scanline and store it in buffer
				TIFFReadScanline(tif, buffer, row);
				for (uint32 col = 0; col < imageWidth; col++) {
					// Grab the pixel value and store it into memory
					unsigned int pixelIndex = col + row * imageWidth;
					//unsigned int pixelIndex = row + col * imageWidth;
					_dm->Data[t][pixelIndex] = ((unsigned short*)(buffer))[col];
					//std::cout << _dm->Data[t][pixelIndex] << std::endl;
				}
			}
			_TIFFfree(buffer);
			TIFFClose(tif);
		}
	}
}

// Saves the particle filter and smoother data as well as the estimated model parameters.
void FileManager::SaveData(unsigned int iterNumber){
  
	// Step 1: 
	// Write the particle data to "results_directory/loc#" where # is iteration number
	char fileName[1024];
	std::string tmpString1 = _resultsDir +  "loc" + boost::lexical_cast<std::string>(iterNumber);
	strncpy(fileName,tmpString1.c_str(),sizeof(fileName));
	fileName[sizeof(fileName)-1] = 0;
	std::ofstream file(fileName, std::fstream::trunc);
	for (size_t t = 0; t < _dm->NumTimesteps(); t++) { 
		//Particle indices 0 through M-1:
		for (size_t part = 0; part < _dm->NumParticles()-1; part++) {
			file << _dm->xParticle[t][part] << " " 
				<< _dm->yParticle[t][part] << " " 
				<< _dm->zParticle[t][part] << " "
				<< _dm->wSmooth[t][part] << " ";
		}
		//Final particle index: 
		file << _dm->xParticle[t][_dm->NumParticles()-1] << " " 
			<< _dm->yParticle[t][_dm->NumParticles()-1] << " " 
			<< _dm->zParticle[t][_dm->NumParticles()-1] << " "
			<< _dm->wSmooth[t][_dm->NumParticles()-1] << std::endl;
	}
	file.close();
  
  	// Step 2: 
  	//Write the parameter data
  	std::string tmpString2 = _resultsDir + "param" + boost::lexical_cast<std::string>(iterNumber);
  	strncpy(fileName,tmpString2.c_str(),sizeof(fileName));
  	fileName[sizeof(fileName)-1] = 0;
  	file.open(fileName, std::fstream::trunc);
  	std::copy(_ic.begin(),_ic.end(),std::ostream_iterator<IInitialConditionModel*>(file));
  	std::copy(_mm.begin(),_mm.end(),std::ostream_iterator<IMotionModel*>(file));
  	std::copy(_obsv.begin(),_obsv.end(),std::ostream_iterator<IObservationModel*>(file));
  	file << std::endl;
  	file.close();
  	
  	// Step 3:
  	// Write the joint density to a separate file
  	std::string tmpString3 = _resultsDir + "joint" + boost::lexical_cast<std::string>(iterNumber);
	strncpy(fileName,tmpString3.c_str(),sizeof(fileName));
  	fileName[sizeof(fileName)-1] = 0;
  	file.open(fileName, std::fstream::trunc);
  	for (size_t i = 0; i < _dm->NumParticles(); i++) {
      	for (size_t j = 0; j < _dm->NumParticles(); j++) {
	  		for (size_t t = 0; t < _dm->NumTimesteps()-1; t++)	{
	      		file << _dm->wJoint[t][i][j] << " ";
	    	}
	  		file << std::endl;
		}
    }
  	file.close();
  
  	// Step 4: 
  	// To save space, remove the previous joint density file
  	if ( iterNumber > 0) {
      	std::string tmpString4 = _resultsDir + "joint" + boost::lexical_cast<std::string>(iterNumber-1);
      	strncpy(fileName,tmpString4.c_str(),sizeof(fileName));
      	fileName[sizeof(fileName)-1] = 0;
      	boost::filesystem::remove(fileName);
    }
}

// Verifies all paths are correctly set up.
void FileManager::VerifyPaths() {
	
	// Make sure the image directory is formatted correctly:
	if(*(_imgDir.end()-1) != '/') {
		_imgDir.push_back('/');
		std::cout << "Warning: Appended '/' to image directory." << std::endl;
	}
	
	// Make sure the image path exists:
	if(!boost::filesystem::exists(_imgDir)) {
		std::cout << "Error: Image directory does not exist." << std::endl;
		throw std::exception();
	}
	
	// Make sure the sensor file exists:
	if(!boost::filesystem::exists(_sensorFilename)) {
		std::cout << "Error: Sensor file does not exist." << std::endl;
		throw std::exception();
	}
	
	// Make sure the results directory is formatted correctly:
	if(*(_resultsDir.end()-1) != '/') {
		_resultsDir.push_back('/');
		std::cout << "Warning: Appended '/' to results directory." << std::endl;
	}
		
	// If the results directory exists, destroy its contents:
	if(boost::filesystem::exists(_resultsDir)) {
		std::cout << "Warning: Results directory exists. Destroying its contents." << std::endl;
		boost::filesystem::remove_all(_resultsDir);
	} 
	
	// Create the results directory:
	boost::filesystem::create_directories(_resultsDir);
	
}
