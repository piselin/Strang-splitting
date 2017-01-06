#include <iostream>
#include "gtest/gtest.h"
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

// class DisassemblyTest : public ::testing::Test {
// protected:
// 	virtual void SetUp() {}

// 	Disassembler d;
// };

// // TEST_F(DisassemblyTest, CountDigits) {
	
// // 	EXPECT_EQ(1, d.CountDigits(0x1));
// // 	EXPECT_EQ(3, d.CountDigits(0x321));
// // 	EXPECT_EQ(8, d.CountDigits(0x3200000f));
// // 	EXPECT_EQ(0, d.CountDigits(0));
// // }


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}