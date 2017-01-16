#include "util.hpp"

namespace util  {

	Eigen::MatrixXcd IntializePotential(const unsigned int n, HarmonicPotential<double>& pot) {
		Eigen::MatrixXcd m(n,n);
		return m;
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
		double eps = 0.01;
		const int n = N;

		std::vector<double> a(n/2, 0);
		std::vector<double> b(n/2, 0);

		std::iota(a.begin(), a.end(), 0);		// fill range evenly spaced starting at 0
		std::iota(b.begin(), b.end(), -n/2);	// fill range evenly spaced starting at -n/2

		for(auto x : b)
			a.push_back(x);

		// calculate x^2*1/2*eps for each element
		std::for_each(a.begin(), a.end(), [&](double& x) {x*=x*0.5*eps;});

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
			a[i] = M_PI*a[i];
		}

		return a;
	}
}