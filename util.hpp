#ifndef STRANG_UTIL_LIB
#define STRANG_UTIL_LIB

#include <vector>
#include <cmath>
#include <complex>	// requires c++14 (gnu++14)
#include <algorithm>
#include <cassert>

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include "types.hpp"
#include "strang.hpp"


namespace util {

	std::vector<double> CreateLaplacian1D(const unsigned int N);
	std::vector<double> CreateGrid1D(const unsigned int N);

	Eigen::MatrixXcd CreateLaplacian(const unsigned int N);

	Eigen::MatrixXcd IntializePotential(const unsigned int N, HarmonicPotential& pot);
}

#endif /* STRANG_UTIL_LIB */