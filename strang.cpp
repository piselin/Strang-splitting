#include "strang.hpp"


//double harmonic(const double x, const double v0=8, const double beta=0.25) {

// PRE: x is a discrete coordinate
double harmonic(const double x) {
	return 0.5*(x*x);
}

// PRE: N is the number of discrete grid points
// POST: return the laplacian (FIXME pi)
std::vector<double> create_laplacian_1d(const unsigned int N) {

	const int n = N;

	std::vector<double> a(n/2, 0);
	std::vector<double> b(n/2, 0);

	std::iota(a.begin(), a.end(), 0);		// fill range evenly spaced starting at 0
	std::iota(b.begin(), b.end(), -n/2);	// fill range evenly spaced starting at -n/2

	for(auto x : b)
		a.push_back(x);

	// calculate x^2*1/2*eps for each element
	std::for_each(a.begin(), a.end(), [](double& x) {x*=x*0.5*eps;});

	return a; 
}

std::vector<double> create_grid_1d(const unsigned int N) {
	std::vector<double> a(N, 0.0);
	a[0] = -1;
	double increment = 2.0/N;
	for(size_t i = 1; i < a.size(); i++) {
		a[i] = a[i-1] + increment;
	}

	// scale 
	for(size_t i = 0; i < a.size(); i++) {
		a[i] = M_PI*a[i];
	}

	return a;
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

std::vector<float> MakeMyData() {
	return {0,1,2,3,4,5,6,7};
}

void split() {
	//auto f = [=] (auto x) {std::cout << x << std::endl;};

	// 1. Set up environment
	auto scale = M_PI;

	const double tend = 5.0*std::sqrt(ieps);
	const int n_timesteps = tend*10;
	std::cout << "Doing " << n_timesteps << " timesteps" << std::endl;

	const double h = tend/n_timesteps; // h is delta_t
	complex_t ih = h*1i;


	// space discretisation
	const unsigned int l = 4;
	const unsigned int N = 1<<l; // calculate 2^l	
	std::cout << N << std::endl;

	// Laplacian
	std::vector<double> laplacian = create_laplacian_1d(N);
	// std::cout << laplacian.back() << std::endl;
	// std::cout << laplacian[laplacian.size()/2-1] << std::endl;
	// std::cout << laplacian[laplacian.size()/2] << std::endl;
	// std::cout << laplacian[laplacian.size()/2+1] << std::endl;

	// Potential
	std::vector<double> x = create_grid_1d(N);
	// for(auto i : x)
	// 	std::cout << i << std::endl;

	std::vector<double> v = initialvalue(x);
	// for(auto i: v)
	// 	std::cout << i << std::endl;

	// /* FFTW */
	// fftw_complex in[v.size()];
	// for (size_t i = 0; i < v.size(); ++i) {
	// 	in[i][0] = v[i];	// real 
	// 	in[i][1] = 0.0;		// imaginary
 //    }

	// fftw_plan p_inplace = fftw_plan_dft_1d(N, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
	// fftw_execute(p_inplace);

	// for(unsigned int i = 0; i < N; i++)
	// 	std::cout << in[i][0] << " " << in[i][1] << std::endl;
	

	/*

	just testing

	*/

	Eigen::FFT<float> fft;

	std::vector<float> timevec = MakeMyData();
	std::vector<std::complex<float>> freqvec;

	std::cout << "forward" << std::endl;

	fft.fwd(freqvec, timevec);
	for(auto x : freqvec)
		std::cout << x << std::endl;

	std::cout << "backward" << std::endl;

	fft.inv(timevec, freqvec);
	for(auto x : timevec)
		std::cout << x << std::endl;

	/*

	just testing

	*/

	

	// 2. strang splitting

	// 3. print solution
}
