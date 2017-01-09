#include "strang.hpp"


//double harmonic(const double x, const double v0=8, const double beta=0.25) {

// PRE: x is a discrete coordinate
double harmonic(const double x) {
	return 0.5*(x*x);
}


// pass function...
std::vector<double> initialvalue(const std::vector<double>& grid) {
	std::vector<double> result(grid.size(), 0.0);

	//////////// g1 ////////////////////////////////
	double a = std::sqrt(std::sqrt(ieps/M_PI));
	for(size_t i = 0; i < grid.size(); i++) {
		double x = grid[i];
		x*=x; // x^2
		result[i] = a*std::exp(-1*0.5*ieps*x);
	}
	////////////////////////////////////////////////

	return result;
}


void split() {
	//auto f = [=] (auto x) {std::cout << x << std::endl;};

	// 1. Set up environment
	const auto scale = M_PI;
	const real_t Dim = 2;

	// time
	const double tend = 5.0*std::sqrt(ieps);
	const int n_timesteps = tend*10;
	std::cout << "Doing " << n_timesteps << " timesteps" << std::endl;
	
	const double delta_t = tend/n_timesteps; // delta_t is h in python
	complex_t delta_t_c = delta_t*1i;				// this is ih in python


	// space discretisation
	const unsigned int l = 4;
	const unsigned int N = 1<<l; // calculate 2^l	
	std::cout << N << std::endl;

	// Laplacian
	std::vector<double> laplacian = util::CreateLaplacian1D(N);
	//CMatrix<N,N> A = laplacian.asDiagonal();
	CMatrix<N,N> A = util::CreateLaplacian(N);
	std::cout << A << std::endl;
	
					// std::cout << laplacian.back() << std::endl;
					// std::cout << laplacian[laplacian.size()/2-1] << std::endl;
					// std::cout << laplacian[laplacian.size()/2] << std::endl;
					// std::cout << laplacian[laplacian.size()/2+1] << std::endl;

	// Potential
	//std::vector<double> x = util::CreateGrid1D(N);
	// for(auto i : x)
	// 	std::cout << i << std::endl;

	// functor
	HarmonicPotential<double> potential;	// this is to V in python
	CMatrix<N,N> V = util::IntializePotential(N, potential);


	// std::vector<double> v = initialvalue(x);
	// for(auto i: v)
	// 	std::cout << i << std::endl;

	Eigen::FFT<double> fft;
	std::vector<std::complex<double>> u;

	// u => freqvec
	// v => timevec
	fft.fwd(u, v);
	// for(auto x : u)
	// 	std::cout << x << std::endl;




	// Exponentials
	A *= delta_t_c;
	CMatrix<N,N> ea = A.exp();
	std::cout << "norm = " << ea.norm() << std::endl;

	CMatrix<N,N> V = 


	/*
	=====================================================
	EXAMPLE eigen fft useage
	*/

	// Eigen::FFT<float> fft;

	// std::vector<float> timevec  = {0,1,2,3,4,5,6,7};
	// std::vector<std::complex<float>> freqvec;

	// std::cout << "forward" << std::endl;

	// fft.fwd(freqvec, timevec);
	// for(auto x : freqvec)
	// 	std::cout << x << std::endl;

	// std::cout << "backward" << std::endl;

	// fft.inv(timevec, freqvec);
	// for(auto x : timevec)
	// 	std::cout << x << std::endl;

	/*
	=====================================================
	*/

	// HarmonicPotential<double> hp;
	// std::cout << hp(2) << std::endl;

	const int D = 2;
	Eigen::Matrix<double, D,D> m1;
	std::cout << m1.size() << std::endl;
	
}
