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

TEST_F(FFtransformTest, JustADummyTest) {
	EXPECT_EQ(1,1.0);
}

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

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}