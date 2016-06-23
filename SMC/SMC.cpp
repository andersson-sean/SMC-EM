#include "SMC.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <iterator>


// Employs a Sampling-Importance-Resampling (SIR) particle filter to the models specified.
// Resampling is performed on every time index.
void SMC::Filter() {
	//using std::cout; using std::endl;
	
	// Initial condition:
	for(size_t part = 0; part < _dm->NumParticles(); ++part) {

		// Sample from initial condition model:
		_dm->xParticle[0][part] = _ic[0]->Sample();
		_dm->yParticle[0][part] = _ic[1]->Sample();
		_dm->zParticle[0][part] = _ic[2]->Sample();
		
		//cout << _dm->xParticle[0][part] << " " << _dm->yParticle[0][part] << " " << _dm->zParticle[0][part] << endl;
		
		// Evaluate (log-)likelihood of each particle:
		std::vector<double> partPos(3);
		partPos[0] = _dm->xParticle[0][part];
		partPos[1] = _dm->yParticle[0][part];
		partPos[2] = _dm->zParticle[0][part];
		double logLikelihood = 0;
		for(size_t pix = 0; pix < _dm->NumPixels(); ++pix) {
			std::vector<double> sensorPos(3);
			sensorPos[0] = _dm->xSensor[0][pix];
			sensorPos[1] = _dm->ySensor[0][pix];
			sensorPos[2] = _dm->zSensor[0][pix];
			double pdfVal = _obsv[0]->Evaluate(
				_dm->Data[0][pix],
				partPos,
				sensorPos);
			//cout << _dm->Data[0][pix] << " ";
			//cout << sensorPos[0] << " " << sensorPos[1] << " " << sensorPos[2] << endl;
			//cout << partPos[0] << " " << partPos[1] << " " << partPos[2] << endl << endl;
			if(pdfVal != 0) logLikelihood += log(pdfVal);
			else logLikelihood += log(1e-300);
		}
		_dm->wFilter[0][part] = logLikelihood;	
	}
	
	// Normalize the weights so they sum to unity:
	NormalizeLog(_dm->wFilter[0]);

	// Generate resampling indices:
	std::vector<double> resamplingIndices;
	GenerateResamplingIndices(_dm->wFilter[0], resamplingIndices);
	
	//std::copy(_dm->wFilter[0].begin(),_dm->wFilter[0].end(),std::ostream_iterator<double>(cout," "));
	//cout << endl;
	//std::copy(resamplingIndices.begin(),resamplingIndices.end(),std::ostream_iterator<double>(cout," "));
	//cout << endl << endl;
	
	// Resample the pertinent vectors:
	Resample(_dm->xParticle[0], resamplingIndices);
	Resample(_dm->yParticle[0], resamplingIndices);
	Resample(_dm->zParticle[0], resamplingIndices);
	
	// Reset the weights to 1.0 / (Number of Particles)
	Equalize(_dm->wFilter[0]);
	
	
	// Remaining time steps:
	for(size_t t = 1; t < _dm->NumTimesteps(); ++t){
		//cout << t << endl;
		for(size_t part = 0; part < _dm->NumParticles(); ++part){
			// Sample from the motion model given the previous particle value:
			_dm->xParticle[t][part] = _mm[0]->Sample(_dm->xParticle[t-1][part]);
			_dm->yParticle[t][part] = _mm[1]->Sample(_dm->yParticle[t-1][part]);
			_dm->zParticle[t][part] = _mm[2]->Sample(_dm->zParticle[t-1][part]);
		
		
			//cout << _dm->xParticle[t][part] << " " << _dm->yParticle[t][part] << " " << _dm->zParticle[t][part] << endl;
			//std::cin.get();
			// Evaluate (log-)likelihood of each particle:
			std::vector<double> partPos(3);
			partPos[0] = _dm->xParticle[t][part];
			partPos[1] = _dm->yParticle[t][part];
			partPos[2] = _dm->zParticle[t][part];
			double logLikelihood = 0;
			for(size_t pix = 0; pix < _dm->NumPixels(); ++pix) {
				std::vector<double> sensorPos(3);
				sensorPos[0] = _dm->xSensor[t][pix];
				sensorPos[1] = _dm->ySensor[t][pix];
				sensorPos[2] = _dm->zSensor[t][pix];
				double pdfVal = _obsv[0]->Evaluate(
					_dm->Data[t][pix],
					partPos,
					sensorPos);
			if(pdfVal != 0) logLikelihood += log(pdfVal);
			else logLikelihood += log(1e-300);
			}
			_dm->wFilter[t][part] = logLikelihood;	
		}
		
		// Normalize the weights so they sum to unity:
		NormalizeLog(_dm->wFilter[t]);
	
		// Generate resampling indices:
		std::vector<double> resamplingIndices;
		GenerateResamplingIndices(_dm->wFilter[t], resamplingIndices);
	
		//std::copy(_dm->wFilter[t].begin(),_dm->wFilter[t].end(),std::ostream_iterator<double>(cout," "));
		//cout << endl;
		//std::copy(resamplingIndices.begin(),resamplingIndices.end(),std::ostream_iterator<double>(cout," "));
		//cout << endl << endl;
		// Resample the pertinent vectors:
		Resample(_dm->xParticle[t], resamplingIndices);
		Resample(_dm->yParticle[t], resamplingIndices);
		Resample(_dm->zParticle[t], resamplingIndices);
	
		// Reset the weights to 1.0 / (Number of Particles)
		Equalize(_dm->wFilter[t]);	
	}
}

// Smooths the resulting weights using the Backward Smoothing algorithm.
void SMC::Smooth() {

	// Pre-calculate transition densities for combined diffusion and directed motion.
	// These are pre-calculated and stored to speed up the backward-smoothing process below. 
	std::vector<std::vector<std::vector<double> > > transDensity;
	transDensity.resize(
		_dm->NumTimesteps()-1,std::vector<std::vector<double> >(
		_dm->NumParticles(),std::vector<double>(
		_dm->NumParticles())));
	for (size_t t = 0; t < _dm->NumTimesteps()-1; ++t)
		for (size_t i = 0; i < _dm->NumParticles(); ++i)
			for (size_t j = 0; j < _dm->NumParticles(); ++j)
				transDensity[t][j][i] = 
					_mm[0]->Evaluate(_dm->xParticle[t+1][j],_dm->xParticle[t][i])*
					_mm[1]->Evaluate(_dm->yParticle[t+1][j],_dm->yParticle[t][i])*
					_mm[2]->Evaluate(_dm->zParticle[t+1][j],_dm->zParticle[t][i]);

  
	// Calculate the smoothed weights by going evaluating them backward in time.
	// The initial condition for the smoothed weights is the same as the filtered weights:
	std::copy(
		_dm->wFilter[_dm->NumTimesteps()-1].begin(),_dm->wFilter[_dm->NumTimesteps()-1].end(), // input
		_dm->wSmooth[_dm->NumTimesteps()-1].begin() // output
	);
  	// Now propagate it through the remaining time samples.
	for (int t = _dm->NumTimesteps()-2; t>=0; t--) {
		for (size_t i = 0; i < _dm->NumParticles(); i++) {
			double v1 = 0;
			for (size_t k = 0; k < _dm->NumParticles(); k++) {
				double v2 = 0;
				for (size_t l = 0; l < _dm->NumParticles(); l++)
					v2 += _dm->wFilter[t][l] * transDensity[t][k][l];
				
				v1 += _dm->wSmooth[t+1][k] * transDensity[t][k][i]/v2;
			}
			_dm->wSmooth[t][i] = _dm->wFilter[t][i] * v1;
		}
    }

	// Calculate joint weights.
	for (size_t t = 0; t < _dm->NumTimesteps()-1; ++t) {
		for (size_t j = 0; j < _dm->NumParticles(); ++j) {
			double v = 0;
			for (size_t l = 0; l < _dm->NumParticles(); ++l)
				v += _dm->wFilter[t][l] * transDensity[t][j][l]; 
			
			for (size_t i = 0; i < _dm->NumParticles(); ++i)
				_dm->wJoint[t][i][j] = 
					_dm->wFilter[t][i] * 
					_dm->wSmooth[t+1][j] * 
					transDensity[t][j][i] / v;
		}
    }
}

// Normalizes a vector of numbers that have been scaled by their natural logarithm.
// For example, if
// 		input = [log(0.1), log(0.01), log(0.001)]
// then the method will scale each element as if it were [0.1, 0.01, 0.001] such that
// 		output = [0.1, 0.01, 0.001] / sum([0.1, 0.01, 0.001])
// This is method is more numerically stable for very large and small values.
bool SMC::NormalizeLog(std::vector<double>& input){
	std::vector<double> tmp(input.size());
	for(size_t i = 0; i < tmp.size(); ++i){
		tmp[i] = 0.0;		
		for(size_t j = 0; j < tmp.size(); ++j)
			tmp[i] += exp(input[j] - input[i]);
	}
	for(size_t i = 0; i < tmp.size(); ++i)
		if(tmp[i] != 0) input[i] = 1.0 / tmp[i];
		else return true;
		
	return false;
}

// Sets all elements of a vector to the inverse of the number of elements
void SMC::Equalize(std::vector<double>& input) {
	size_t numEl = input.size();
	double invVal = 1.0 / static_cast<double>(numEl);
	std::fill(
		input.begin(), input.end(),
		invVal
	);
}

// Generates a vector of resampling indices given an input weight vector 
void SMC::GenerateResamplingIndices(
	const std::vector<double>& input, // input array
	std::vector<double>& output) { // output array
	
	// Make sure the output is of the same size as the input:
	output.resize(input.size());

	//Generate a vector of Uniform(0,1) random variates 
	//by first calling rand() and then by multiplying each element
	//by 1.0/RAND_MAX.  Recall that rand() generates discrete unsigned 
	//uniform [0,RAND_MAX).
	std::vector<double> UnifRandVar(input.size());;
	std::generate(
		UnifRandVar.begin(), UnifRandVar.end(), // output
		rand // method
	);
	const double RAND_MAX_RECIPROCAL = 1.0 / static_cast<double>(RAND_MAX);
	std::transform(
		UnifRandVar.begin(), UnifRandVar.end(), // input
		UnifRandVar.begin(), // output
			std::bind2nd(std::multiplies<double>(), RAND_MAX_RECIPROCAL) // method
	);
   
	// Calculate the cumulative sum of the likelihoods stored in input array.
	std::vector<double> CumSum(input.size());
	std::partial_sum(
		input.begin(), input.end(), // input
		CumSum.begin() // output
	);

	// Generate resampled indices.
	// Iterators make this a bit difficult to understand.  Here, we:
	// 1) Iterate through each of the uniform random variables (URVs) in UnifRandVar.
	// 2) Find the smallest value in Cumulative_Sum that is larger than the current URV.
	// 3) Store the resulting index in the output vector.
	std::vector<double>::iterator uv_it;
	for (uv_it = UnifRandVar.begin() ; uv_it != UnifRandVar.end() ; ++uv_it) {
		const size_t URV_index = uv_it - UnifRandVar.begin();
		std::vector<double>::iterator cs_it = 
			std::find_if(
				CumSum.begin(), CumSum.end(), // input
				std::bind1st(std::less<double>(),*uv_it) // method
			);
			output[URV_index] = cs_it - CumSum.begin();
    }
}


// Resamples a vector according to specified indices.
void SMC::Resample(
	std::vector<double>& input, // input array
	const std::vector<double>& indices) { // array of indices by which to resample
	
	// Store a temporary copy of the row
	std::vector<double> tmp(input.size());
	std::copy(
		input.begin(), input.end(), // input
		tmp.begin() // output
	);
  
	// Replace the values according to the resampling indices.
	std::vector<double>::iterator it;
	for (it = input.begin() ; it != input.end() ; ++it)
		*it = tmp[indices[it - input.begin()]];
}

