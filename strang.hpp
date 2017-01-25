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


#include "types.hpp"		// from libwaveblocks
#include "hdf5writer.hpp"


#define SIMULATION_1D ;
//#define SIMULATION_2D;

template <dim_t N>
class StrangSplitter {
public:
	/**
	Create the Simulator
	with the most important parameters
	*/
	StrangSplitter(const double epsilon, const double T, const double dt)
		:
		n_(N),
		eps_(epsilon),
		tend_(T),
		dt_(dt)
		{
			n_timesteps_ = tend_ / dt_;
			dt_complex_ = dt_*1i;
			ieps_ = 1/eps_; 

			CreateLaplacian();

			CreateGrid();

			InitializePotential();

			InitialValue();

			InitializeExponentialA();

			InitializeExponentialB();
		}

		void Split() {
			// Hdf5Writer<4> writer("strang_1D_cpp.hdf5");
			// writer.tester();
			Eigen::FFT<double> fft;

			CMatrix<1,N> u_freq;

			// Strang Splitting Algorithm
			for(size_t step = 1; step < n_timesteps_; ++step) {
			
				// 1. u = eb*u
				for(size_t i = 0; i < n_; ++i)
					u_[i] = eb_[i]*u_[i];

				// 2. u = fft(u)
				fft.fwd(u_freq,u_);

				// 3. u = ea*u
				for(size_t i = 0; i < n_; ++i)
					u_freq[i] = ea_[i]*u_freq[i];

				// 4. u = ifft(u)
				fft.inv(u_,u_freq);
				
				// 5. u = eb*u
				for(size_t i = 0; i < n_; ++i)
					u_[i] = eb_[i]*u_[i];
			}


				//writer.StoreResult();
				//std::cout << "Norm in " << step << " = " << norm(v) << std::endl;

			std::cout << "after" << std::endl;
			std::cout << u_ << std::endl;

			std::cout << "Norm = " << u_.norm() << std::endl;
		}
	
private: /* Member Functions */

	/**
	Pre: 	Number of gridpoints and eps_ initialized.
	Post: 	In 1D the Laplacian has the form [0,1,4,9,....,9,4,1] * k
			where k depends on the eps_.
	*/
	void CreateLaplacian() {
		const int n = n_;

		std::vector<double> a(n/2, 0);
		std::vector<double> b(n/2, 0);

		std::iota(a.begin(), a.end(), 0);		// fill range evenly spaced starting at 0
		std::iota(b.begin(), b.end(), -1*n/2);	// fill range evenly spaced starting at -n/2

		for(auto x : b)
			a.push_back(x);

		// calculate x^2*1/2*eps for each element
		std::for_each(a.begin(), a.end(), [&](double& x) {x*=x*0.5*eps_;});
		for(size_t i = 0; i < n_; i++) {
			laplacian_[i] = a[i];
		}
	}

	/**
	Post:	The 1D grid is of the form [-N/2, -N/2+1,...,-1,0,1,2,...N/2-1]
			where N = n_ ,i.e. the number of gridpoints.
	*/
	void CreateGrid() {
		grid_[0] = -1;
		double increment = 2.0/n_;
		for(size_t i = 1; i < n_; i++) {
			grid_[i] = grid_[i-1] + increment;
		}

		// scale 
		for(size_t i = 0; i < n_; i++) {
			grid_[i] = grid_scale_*grid_[i];
		}

	}

	/**
	Fixme pi: The potential function should be exchangeable
	Pre:	Gridpoints have to be defined
	Post:	v_pot_ contains the potential function evaluated
			at each grid point.
	*/
	void InitializePotential() {
		HarmonicPotential<double> potential;
		for(size_t i = 0; i < n_; ++i)
			v_pot_[i] = ieps_*potential(grid_[i]);
	}

	/**
	Fixme pi: the initial value should be exchangeable
	Post: 	creates v(t=0, x) for each grid point. The initial value
			is currently of the form g(x) = a * e^(-1/2 x^2)
	*/
	void InitialValue() {
		auto a = std::sqrt(std::sqrt(ieps_/grid_scale_));
		for(size_t i = 0; i < n_; i++) {
			auto x = grid_[i];
			x*=x; // x^2
			u_[i] = a*std::exp(-1*0.5*ieps_*x);
		}
	}

	/**
	Post:	Define exponentials according to 
	*/
	void InitializeExponentialA() {
		for(size_t i = 0; i < n_; ++i)
			ea_[i] = std::exp(-dt_complex_*laplacian_[i]);
	}

	void InitializeExponentialB() {
		for(size_t i = 0; i < n_; ++i)
			eb_[i] = std::exp(-0.5*dt_complex_*v_pot_[i]);
	}

	/**
	The Potential for 
	*/
	template<typename T>
	class HarmonicPotential {
	public:
		real_t operator()(const T& x) const {
			return 0.5*x*x;
		}
	};



public: /*Data Members*/
	/**
	Fixme pi: 	this should not be public. It is currently public for easy
				GTesting. Remember, unit testing should just be a black-box
				test. So change the gtest unit tests?
	*/
	#ifdef SIMULATION_1D
		const dim_t dim_ = 1;
		CMatrix<1,N> u_;	// the solution of the schrödinger equation
		RMatrix<1,N> laplacian_;
		RMatrix<1,N> grid_;
		RMatrix<1,N> v_pot_;
		CMatrix<1,N> ea_;
		CMatrix<1,N> eb_;
	#endif

private:
	#ifdef SIMULATION_2D
		// not yet supported
	#endif

	//const unsigned int dim_;
	const unsigned int n_; // number of grid points
	const double eps_ = 0.01;
	double ieps_ = 1.0/eps_;
	const double grid_scale_ = M_PI;
	const double tend_;
	const double dt_;
	unsigned int n_timesteps_;
	complex_t dt_complex_;
};

#endif /* STRANG_SPLITTNG_LIB */