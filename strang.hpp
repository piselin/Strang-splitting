#ifndef STRANG_SPLITTING_LIB
#define STRANG_SPLITTING_LIB

#include <iostream>
#include <cmath>
#include <complex>	// requires c++14
#include <algorithm>
#include <vector>

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/FFT>

//#include <fftw3.h>
//#include "../kiss_fft130/kiss_fft.h"

#include "types.hpp"

// Parameter
// fixme pi: move to config file?
const double kEps = 0.01;
const double kIEps = 1.0/kEps;
const auto kScaleGrid = M_PI;


//fixme: should be some OO based desing!

// class StrangSplitter {
// public:
// 	StrangSplitter(const unsigned int dimension, const double epsilon)
// 		:
// 		dim_(dimension),
// 		eps_(epsilon)
// 		{}
	

// private:
// 	const unsigned int dim_;
// 	const double eps_ = 0.01;
// 	const double ieps_ = 1.0/eps_;

// };

// PRE: x is a discrete coordinate
double harmonic(const double x);

// pass function...
std::vector<std::complex<double>> Initialvalue(const std::vector<double>& grid);

// function object (functor)
template<typename T>
class HarmonicPotential {
public:
	real_t operator()(const T& x) const {
		return 0.5*x*x;
	}
};


std::vector<double> CreateLaplacian1D(const unsigned int N);
std::vector<double> CreateGrid1D(const unsigned int N);

Eigen::MatrixXcd CreateLaplacian(const unsigned int N);

template<typename T>
std::vector<double> IntializePotential(const unsigned int n, HarmonicPotential<T>& pot);
//Eigen::MatrixXcd IntializePotential(const unsigned int N, HarmonicPotential<T>& pot);
double norm(const std::vector<std::complex<double>>& v);

void InitializeExponentialA(
		const size_t N, 
		const std::complex<double> k, 
		const std::vector<double>& laplacian, 
		std::vector<std::complex<double>>& ea);

void InitializeExponentialB(
		const size_t N, 
		const std::complex<double> k, 
		const std::vector<double>& V, 
		std::vector<std::complex<double>>& eb);


void split(
	const unsigned int N,	
	const int n_timesteps,
	const std::vector<std::complex<double>>& ea,
	const std::vector<std::complex<double>>& eb,
	std::vector<std::complex<double>>& v,
	std::vector<std::complex<double>>& v_freq);

void run();

#endif /* STRANG_SPLITTNG_LIB */