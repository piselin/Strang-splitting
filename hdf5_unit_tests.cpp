#include <iostream>
#include <string>
#include <vector>
#include "gtest/gtest.h"

#include "H5Cpp.h"

// my own libraries
#include "hdf5writer.hpp"

using namespace H5;

class TestHDF : public ::testing::Test {
protected: // according to gtest primer
	TestHDF() {
		std::cout << "ctor" << std::endl;
	}

	virtual ~TestHDF() {
		std::cout << "dtor" << std::endl;
		//testfile_.close();
	}

	void SetUp() {

	}
	H5File testfile_;
	H5std_string filename_;

};


TEST_F(TestHDF, SomeTest) {
	EXPECT_EQ(0, 0.);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}