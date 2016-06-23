#ifndef Utilities_H
#define Utilities_H

namespace emlib {

	// Template function for calculating the signum of a number
	template <typename T> int sgn(T val) {
    	return (T(0) < val) - (val < T(0));
	}
}

#endif
