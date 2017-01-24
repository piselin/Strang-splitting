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

class StrangSplitterTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		eps = 0.01;
		ieps = 1./eps;
		tend = 5.0*std::sqrt(ieps);
		n_timesteps = tend*10;
		delta_t = tend/n_timesteps; // delta_t is h in python
		N = 16;

	}
	double eps;
	double ieps;
	double tend;
	double n_timesteps;
	double delta_t; 
	unsigned int N;
};


TEST_F(StrangSplitterTest, CheckLaplacian) {
	StrangSplitter<16> system(eps, tend, delta_t);


	auto laplacian = system.laplacian_;
	unsigned size = laplacian.size();
	std::vector<double> ref = {0, 0.005, 0.02, 0.045, 0.08, 0.125, 0.18, 0.245, 0.32, 0.245, 0.18, 0.125, 0.08, 0.045, 0.02, 0.005};

	EXPECT_EQ(size, ref.size());

	for(unsigned int i = 0; i < size; ++i)
		EXPECT_FLOAT_EQ(ref[i], laplacian[i]);

}

TEST_F(StrangSplitterTest, GridContent) {
	StrangSplitter<16> system(eps, tend, delta_t);
	// //auto grid = CreateGrid1D(grid_size);
	std::vector<double> ref = {-3.14159265, -2.74889357, -2.35619449, -1.96349541, -1.57079633,
       -1.17809725, -0.78539816, -0.39269908,  0.        ,  0.39269908,
        0.78539816,  1.17809725,  1.57079633,  1.96349541,  2.35619449,
        2.74889357};

    for(unsigned int i = 0; i < ref.size(); ++i)
		EXPECT_FLOAT_EQ(system.grid_[i], ref[i]);
}

TEST_F(StrangSplitterTest, InitialValue) {
	StrangSplitter<16> system(eps, tend, delta_t);

	EXPECT_FLOAT_EQ(2.375268, system.v_.norm());
}

TEST_F(StrangSplitterTest, ExponentialA) {
	StrangSplitter<16> system(eps, tend, delta_t);

	EXPECT_FLOAT_EQ(4.0, system.ea_.norm());
}

TEST_F(StrangSplitterTest, ExponentialB) {
	StrangSplitter<16> system(eps, tend, delta_t);

	EXPECT_FLOAT_EQ(4.0, system.eb_.norm());
}



// class SetupTester : public ::testing::Test {
// protected:
// 	virtual void SetUp() {
// 		tend = 5.0*std::sqrt(kIEps);
// 		n_timesteps = tend*10;	
// 		delta_t = tend/n_timesteps; // delta_t is h in python;
// 		delta_t_c = delta_t*1i;		
// 		N = 16;
// 	}

// 	double tend;
// 	int n_timesteps;
// 	double delta_t;
// 	complex_t delta_t_c;

// 	int N;

// };

// TEST(Util, LaplacianContent) {
// 	const unsigned size = 16;

// 	auto laplacian = CreateLaplacian1D(size);
// 	std::vector<double> ref = {0, 0.005, 0.02, 0.045, 0.08, 0.125, 0.18, 0.245, 0.32, 0.245, 0.18, 0.125, 0.08, 0.045, 0.02, 0.005};

// 	// sum up the elemets and compare them
// 	double s1 = std::accumulate(laplacian.begin(), laplacian.end(), 0);
// 	double s2 = std::accumulate(ref.begin(), ref.end(), 0);
// 	EXPECT_EQ(s1, s2);

// 	for(unsigned int i = 0; i < size; ++i)
// 		EXPECT_FLOAT_EQ(ref[i], laplacian[i]);
// }


// TEST(Util, GridContent) {
// 	const unsigned grid_size = 16;
// 	auto grid = CreateGrid1D(grid_size);
// 	std::vector<double> ref = {-3.14159265, -2.74889357, -2.35619449, -1.96349541, -1.57079633,
//        -1.17809725, -0.78539816, -0.39269908,  0.        ,  0.39269908,
//         0.78539816,  1.17809725,  1.57079633,  1.96349541,  2.35619449,
//         2.74889357};

//     for(unsigned int i = 0; i < grid_size; ++i)
// 		EXPECT_FLOAT_EQ(ref[i], grid[i]);
// }

// TEST(Util, Norm) {
// 	std::vector<std::complex<double>> v;

// 	double n = 0;
// 	v.push_back(0);
// 	EXPECT_FLOAT_EQ(0, n*n);

// 	v.push_back(1+1i);
// 	n = norm(v);
// 	EXPECT_FLOAT_EQ(2.0, n*n);
// }

// TEST_F(SetupTester, ExponentialA) {

// 	std::vector<double> laplacian = CreateLaplacian1D(N);
// 	std::vector<std::complex<double>> ea(N, 0.0);
// 	InitializeExponentialA(N, delta_t_c, laplacian, ea);

// 	double n = norm(ea);
// 	EXPECT_FLOAT_EQ(4.0, n);	// from python
// }

// TEST_F(SetupTester, ExponentialB) {

// 	std::vector<double> x = CreateGrid1D(N);
// 	// functor
// 	HarmonicPotential<double> potential;
// 	//CMatrix<N,N> V = IntializePotential(N, potential);
// 	//std::vector<double> pot = IntializePotential(N, potential);
// 	std::vector<double> V(N, 0.0);
// 	for(size_t i = 0; i < N; ++i) {
// 		V[i] = kIEps*potential(x[i]);
// 	}

// 	std::vector<std::complex<double>> eb(N, 0.0);
// 	InitializeExponentialB(N, delta_t_c, V, eb);

// 	double n = norm(eb);
// 	EXPECT_FLOAT_EQ(4.0, n);
// }

// TEST_F(SetupTester, InitialSolution) {
// 	std::vector<double> x = CreateGrid1D(N);
// 	std::vector<std::complex<double>> v = Initialvalue(x);

// 	EXPECT_FLOAT_EQ(2.375268, norm(v));
// }

// TEST_F(SetupTester, TestSplitting) {
// 	// store initial configuration
// 	//writer.StoreResult();

// 	/*
// 	for(timesteps){
// 		strang_splitting();
// 		writer.StoreResult();
// 	}
// 	writer.Finalize();
// 	*/
// 	// -4.18727594e-04 +4.30123296e-04j   1.02119717e-02 +2.36600989e-03j
// 	// -4.65856888e-04 -4.66465895e-04j  -5.72318819e-04 +1.67151603e-04j
// 	// 5.77668735e-03 +5.22103975e-05j   9.06869911e-04 -2.91216590e-03j
// 	// -7.85651928e-06 +1.84632275e-04j   1.00783377e-02 +3.66927688e-02j
// 	// 1.33689647e+00 +1.96249916e+00j   1.00783377e-02 +3.66927688e-02j
// 	// -7.85651928e-06 +1.84632275e-04j   9.06869911e-04 -2.91216590e-03j
// 	// 5.77668735e-03 +5.22103975e-05j  -5.72318819e-04 +1.67151603e-04j
// 	// -4.65856888e-04 -4.66465895e-04j   1.02119717e-02 +2.36600989e-03j]
// 	EXPECT_FLOAT_EQ(Split(), 2.375268);
// }

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}