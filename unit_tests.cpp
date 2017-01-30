#include <iostream>
#include <vector>
#include "gtest/gtest.h"

// my own libraries
#include "strang.hpp"

class StrangSplitterTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		eps_ = 0.01;
		ieps_ = 1./eps_;
		tend_ = 5.0*std::sqrt(ieps_);
		n_timesteps_ = tend_*10;
		dt_ = tend_/n_timesteps_; // delta_t is h in python
		
	}
	double eps_;
	double ieps_;
	double tend_;
	unsigned int n_timesteps_;
	double dt_; 
};

TEST_F(StrangSplitterTest, CheckTimeReversal) {
	const auto exponent = 4;
	const auto size = 2 << exponent;

	StrangSplitter<size> system(eps_, tend_, dt_);

	// Advance by 1 step to get the Solution at t=1
	system.Advance();

	// Store results
	auto u1 = system.GetSolution();
	auto norm_initial = system.Norm();

	
	// Run the System until t = T
	for(unsigned int step = 0; step < n_timesteps_; ++step){
		system.Advance();
	}

	system.Reverse();
	
	// Run the System "backwards" until t=1 again
	for(unsigned int step = n_timesteps_; step > 0; --step) {
		system.Advance();
	}

	// Store results of the reversed run
	auto u2 = system.GetSolution();
	auto norm_end = system.Norm();

	// compare
	EXPECT_FLOAT_EQ(norm_initial, norm_end);
	EXPECT_EQ(u2.size(), u1.size());
	for(unsigned int i = 0; i < u1.size(); ++i) {
		EXPECT_FLOAT_EQ(u1[i].real(), u2[i].real());
		EXPECT_FLOAT_EQ(u1[i].imag(), u2[i].imag());
	}
}

TEST_F(StrangSplitterTest, CheckDifferentSystemSizes) {
	const auto max_exponent = 13;
	const auto exponent = 8;
	assert(exponent <= max_exponent);

	const auto size = 2 << exponent;
	
	// we check the sizes
	// 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096
	std::map<unsigned int, double> reference;
	reference[4] = 2.37526752924;
	reference[8] =  2.37526752924;
	reference[16] = 2.37526800605;
	reference[32] = 2.42502403458;
	reference[64] = 3.1916522201;
	reference[128] = 4.51351666838;
	reference[256] = 6.38307648642;
	reference[512] = 9.02703333676;
	reference[1024] = 12.7661529728;
	reference[2048] = 18.0540666735;
	reference[4096] = 25.5323059457;
	reference[8192] = 36.1081333471;

	StrangSplitter<size> system(eps_, tend_, dt_);

	for(real_t t = dt_; t < tend_; t+=dt_) {
		system.Advance();
	}

	EXPECT_FLOAT_EQ(system.Norm(), reference[size]);

}

/* 

These following Tests are no longer usefull because they violate the
principle of data encapsulation. They were usefull during developement of
The StrangSplitter class. But now the class should be tested only
with true unit tests. These are Blackbox tests thus we should only test
against the final result.

*/

////////////////////////////////////////////////////////////////////////
// TEST_F(StrangSplitterTest, CheckLaplacian) {
// 	StrangSplitter<16> system(eps, tend, delta_t);


// 	auto laplacian = system.laplacian_;
// 	unsigned size = laplacian.size();
// 	std::vector<double> ref = {0, 0.005, 0.02, 0.045, 0.08, 0.125, 0.18, 0.245, 0.32, 0.245, 0.18, 0.125, 0.08, 0.045, 0.02, 0.005};

// 	EXPECT_EQ(size, ref.size());

// 	for(unsigned int i = 0; i < size; ++i)
// 		EXPECT_FLOAT_EQ(ref[i], laplacian[i]);

// }

// TEST_F(StrangSplitterTest, GridContent) {
// 	StrangSplitter<16> system(eps, tend, delta_t);
// 	// //auto grid = CreateGrid1D(grid_size);
// 	std::vector<double> ref = {-3.14159265, -2.74889357, -2.35619449, -1.96349541, -1.57079633,
//        -1.17809725, -0.78539816, -0.39269908,  0.        ,  0.39269908,
//         0.78539816,  1.17809725,  1.57079633,  1.96349541,  2.35619449,
//         2.74889357};

//     for(unsigned int i = 0; i < ref.size(); ++i)
// 		EXPECT_FLOAT_EQ(system.grid_[i], ref[i]);
// }

// TEST_F(StrangSplitterTest, InitialValue) {
// 	StrangSplitter<16> system(eps, tend, delta_t);

// 	EXPECT_FLOAT_EQ(2.375268, system.u_.norm());
// }

// TEST_F(StrangSplitterTest, ExponentialA) {
// 	StrangSplitter<16> system(eps, tend, delta_t);

// 	EXPECT_FLOAT_EQ(4.0, system.ea_.norm());
// }

// TEST_F(StrangSplitterTest, ExponentialB) {
// 	StrangSplitter<16> system(eps, tend, delta_t);

// 	EXPECT_FLOAT_EQ(4.0, system.eb_.norm());
// }
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}