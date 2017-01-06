#ifndef STRANG_SPLITTING_LIB
#define STRANG_SPLITTING_LIB

#include <iostream>
#include <cmath>
#include <complex>	// requires c++14
#include <algorithm>
#include <vector>

#include <Eigen/Core>

//#include <fftw3.h>
#include <unsupported/Eigen/FFT>

//#include "../kiss_fft130/kiss_fft.h"

#include "types.hpp"


const double eps = 0.01;
const double ieps = 1.0/eps;

//double harmonic(const double x, const double v0=8, const double beta=0.25) {

// PRE: x is a discrete coordinate
double harmonic(const double x);

// PRE: N is the number of discrete grid points
// POST: return the laplacian (FIXME pi)
std::vector<double> create_laplacian_1d(const unsigned int N);

std::vector<double> create_grid_1d(const unsigned int N);

// pass function...
std::vector<double> initialvalue(const std::vector<double>& grid);

void split();


#endif /* STRANG_SPLITTNG_LIB */