#include <iostream>
#include "gtest/gtest.h"

// my own libraries
#include "strang.hpp"

class FFtransformTest : public ::testing::Test {
	// setup
	virtual void SetUp() {}
	// if we would write
	//virtual void SetUp() = 0;
	// this would be a pure virtual function, making the class abstract.
	// not sure if this would work for GTest
};


class SetupTester : public ::testing::Test {
protected:
	virtual void SetUp() {
		tend = 5.0*std::sqrt(ieps);
		n_timesteps = tend*10;	
		delta_t = tend/n_timesteps; // delta_t is h in python;
		delta_t_c = delta_t*1i;		
		N = 16;
	}

	double tend;
	int n_timesteps;
	double delta_t;
	complex_t delta_t_c;

	int N;

};

TEST(Util, LaplacianContent) {
	const unsigned size = 16;

	auto laplacian = CreateLaplacian1D(size);
	std::vector<double> ref = {0, 0.005, 0.02, 0.045, 0.08, 0.125, 0.18, 0.245, 0.32, 0.245, 0.18, 0.125, 0.08, 0.045, 0.02, 0.005};

	// sum up the elemets and compare them
	double s1 = std::accumulate(laplacian.begin(), laplacian.end(), 0);
	double s2 = std::accumulate(ref.begin(), ref.end(), 0);
	EXPECT_EQ(s1, s2);

	for(unsigned int i = 0; i < size; ++i)
		EXPECT_FLOAT_EQ(ref[i], laplacian[i]);
}


TEST(Util, GridContent) {
	const unsigned grid_size = 16;
	auto grid = CreateGrid1D(grid_size);
	std::vector<double> ref = {-3.14159265, -2.74889357, -2.35619449, -1.96349541, -1.57079633,
       -1.17809725, -0.78539816, -0.39269908,  0.        ,  0.39269908,
        0.78539816,  1.17809725,  1.57079633,  1.96349541,  2.35619449,
        2.74889357};

    for(unsigned int i = 0; i < grid_size; ++i)
		EXPECT_FLOAT_EQ(ref[i], grid[i]);
}

TEST(Util, Norm) {
	std::vector<std::complex<double>> v;

	double n = 0;
	v.push_back(0);
	EXPECT_FLOAT_EQ(0, n*n);

	v.push_back(1+1i);
	n = norm(v);
	EXPECT_FLOAT_EQ(2.0, n*n);
}

TEST_F(SetupTester, ExponentialA) {

	std::vector<double> laplacian = CreateLaplacian1D(N);
	std::vector<std::complex<double>> ea(N, 0.0);
	InitializeExponentialA(N, delta_t_c, laplacian, ea);

	double n = norm(ea);
	EXPECT_FLOAT_EQ(4.0, n);	// from python
}

TEST_F(SetupTester, ExponentialB) {

	std::vector<double> x = CreateGrid1D(N);
	// functor
	HarmonicPotential<double> potential;
	//CMatrix<N,N> V = IntializePotential(N, potential);
	//std::vector<double> pot = IntializePotential(N, potential);
	std::vector<double> V(N, 0.0);
	for(size_t i = 0; i < N; ++i) {
		V[i] = ieps*potential(x[i]);
	}

	std::vector<std::complex<double>> eb(N, 0.0);
	InitializeExponentialB(N, delta_t_c, V, eb);
	
	double n = norm(eb);
	EXPECT_FLOAT_EQ(4.0, n);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}