#include "strang.hpp"


//double harmonic(const double x, const double v0=8, const double beta=0.25) {

// PRE: x is a discrete coordinate
double harmonic(const double x) {
	return 0.5*(x*x);
}


// pass function...
std::vector<std::complex<double>> Initialvalue(const std::vector<double>& grid) {
	std::vector<std::complex<double>> result(grid.size(), 0.0);

	//////////// g1 initial value  /////////////////
	double a = std::sqrt(std::sqrt(kIEps/M_PI));
	for(size_t i = 0; i < grid.size(); i++) {
		double x = grid[i];
		x*=x; // x^2
		result[i] = a*std::exp(-1*0.5*kIEps*x);
	}
	////////////////////////////////////////////////

	return result;
}


// Eigen::MatrixXcd IntializePotential(const unsigned int n, HarmonicPotential<double>& pot) {
// 	Eigen::MatrixXcd m(n,n);
// 	return m;
// }
std::vector<double> IntializePotential(const unsigned int n, HarmonicPotential<double>& pot) {
	std::vector<double> V(n, 0.0);

	return V;
}

Eigen::MatrixXcd CreateLaplacian(const unsigned int n) {
	assert(n>0);
	std::vector<double> l = CreateLaplacian1D(n);
	Eigen::MatrixXcd m(n,n);

	for(unsigned int i = 0; i < n; ++i)
		m(i,i) = l[i];
	return m;
}


// PRE: N is the number of discrete grid points
// POST: return the laplacian (FIXME pi)
std::vector<double> CreateLaplacian1D(const unsigned int N) {
	double kEps = 0.01;
	const int n = N;

	std::vector<double> a(n/2, 0);
	std::vector<double> b(n/2, 0);

	std::iota(a.begin(), a.end(), 0);		// fill range evenly spaced starting at 0
	std::iota(b.begin(), b.end(), -n/2);	// fill range evenly spaced starting at -n/2

	for(auto x : b)
		a.push_back(x);

	// calculate x^2*1/2*kEps for each element
	std::for_each(a.begin(), a.end(), [&](double& x) {x*=x*0.5*kEps;});

	return a; 
}


std::vector<double> CreateGrid1D(const unsigned int N) {
	std::vector<double> a(N, 0.0);
	a[0] = -1;
	double increment = 2.0/N;
	for(size_t i = 1; i < a.size(); i++) {
		a[i] = a[i-1] + increment;
	}

	// scale 
	for(size_t i = 0; i < a.size(); i++) {
		a[i] = kScaleGrid*a[i];
	}

	return a;
}

double norm(const std::vector<std::complex<double>>& v) {
	double r = 0.;
	for(size_t i = 0; i < v.size(); ++i) {
		r += std::norm(v[i]);
	}
	
	return std::sqrt(r);
}

void InitializeExponentialA(
		const size_t N, 
		const std::complex<double> k, 
		const std::vector<double>& laplacian, 
		std::vector<std::complex<double>>& ea) 
{

	for(size_t i = 0; i < N; ++i) {
		ea[i] = std::exp(-k*laplacian[i]);
	}
}

void InitializeExponentialB(
		const size_t N, 
		const std::complex<double> k, 
		const std::vector<double>& V, 
		std::vector<std::complex<double>>& eb)
{

	for(size_t i = 0; i < N; ++i) {
		eb[i] = std::exp(-0.5*k*V[i]);
	}
}


void split() {
	//auto f = [=] (auto x) {std::cout << x << std::endl;};

	// 1. Set up environment
	//const auto scale = M_PI;
	const real_t Dim = 2;

	// time
	const double tend = 5.0*std::sqrt(kIEps);
	const int n_timesteps = tend*10;
	//std::cout << "Doing " << n_timesteps << " timesteps" << std::endl;
	
	const double delta_t = tend/n_timesteps; // delta_t is h in python
	complex_t delta_t_c = delta_t*1i;				// this is ih in python


	// space discretisation
	const unsigned int l = 4;
	const unsigned int N = 1<<l; // calculate 2^l	
	//std::cout << N << std::endl;

	// Laplacian
	std::vector<double> laplacian = CreateLaplacian1D(N);
	//CMatrix<N,N> A = laplacian.asDiagonal();
	//CMatrix<N,N> A = CreateLaplacian(N);
	//std::cout << A << std::endl;
	
					// std::cout << laplacian.back() << std::endl;
					// std::cout << laplacian[laplacian.size()/2-1] << std::endl;
					// std::cout << laplacian[laplacian.size()/2] << std::endl;
					// std::cout << laplacian[laplacian.size()/2+1] << std::endl;

	// Potential
	std::vector<double> x = CreateGrid1D(N);
	// for(auto i : x)
	// 	std::cout << i << std::endl;

	// functor
	HarmonicPotential<double> potential;
	//CMatrix<N,N> V = IntializePotential(N, potential);
	//std::vector<double> pot = IntializePotential(N, potential);
	std::vector<double> V(N, 0.0);
	for(size_t i = 0; i < N; ++i) {
		V[i] = kIEps*potential(x[i]);
	}
	// std::cout << "pot" << std::endl;
	// for(auto x : V)
	// 	std::cout << x << std::endl;

	std::vector<std::complex<double>> v = Initialvalue(x);
	// for(auto i: v)
	// 	std::cout << i << std::endl;

	Eigen::FFT<double> fft;
	std::vector<std::complex<double>> u;

	// u => freqvec
	// v => timevec
	// std::cout << "before" << std::endl;
	// for(auto x : v)
	// 	std::cout << x << std::endl;

	fft.fwd(u, v);

	// std::cout << "u" << std::endl;
	// for(auto x : u)
	// 	std::cout << x << std::endl;

	// std::cout << "after" << std::endl;
	// for(auto x : v)
	// 	std::cout << x << std::endl;



	// Exponentials
	std::vector<std::complex<double>> ea(N, 0.0);
	InitializeExponentialA(N, delta_t_c, laplacian, ea);

	std::vector<std::complex<double>> eb(N, 0.0);
	InitializeExponentialB(N, delta_t_c, V, eb);

	// for(auto x : eb)
	// 	std::cout << x << std::endl;


	// Strang Splitting Algorithm
	for(size_t step = 1; step <= n_timesteps; ++step) {
		//std::cout << "running" << std::endl;
		// 1. v = eb*v
		for(size_t i = 0; i < N; ++i)
			v[i] = eb[i]*v[i];

		// 2. v = fft(v)
		fft.fwd(v,v);

		// 3. v = ea*v
		for(size_t i = 0; i < N; ++i)
			v[i] = ea[i]*v[i];

		// 4. v = ifft(v)
		fft.inv(v,v);

		// 5. v = eb*v
		for(size_t i = 0; i < N; ++i)
			v[i] = eb[i]*v[i];

		std::cout << "Norm in " << step << " = " << norm(v) << std::endl;

	}

	std::cout << "after" << std::endl;
	for(auto x : v)
		std::cout << x << std::endl;

	std::cout << "Norm = " << norm(v) << std::endl;


	//std::vector<std::complex<double>> eb(N, 0.0);


	// A *= delta_t_c;
	// CMatrix<N,N> ea = A.exp();
	// std::cout << "norm = " << ea.norm() << std::endl;

	//CMatrix<N,N> V = 


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

	// const int D = 2;
	// Eigen::Matrix<double, D,D> m1;
	// std::cout << m1.size() << std::endl;
	
}
