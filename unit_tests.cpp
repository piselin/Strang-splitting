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

	EXPECT_FLOAT_EQ(2.375268, system.u_.norm());
}

TEST_F(StrangSplitterTest, ExponentialA) {
	StrangSplitter<16> system(eps, tend, delta_t);

	EXPECT_FLOAT_EQ(4.0, system.ea_.norm());
}

TEST_F(StrangSplitterTest, ExponentialB) {
	StrangSplitter<16> system(eps, tend, delta_t);

	EXPECT_FLOAT_EQ(4.0, system.eb_.norm());
}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}