#include "H5Cpp.h" // CPP header for HDF interface

#include <string>
#include <vector>

/* remove me */
#include <iostream>

#include <Eigen/Core>


using namespace H5;

/**
Main functionality is to store the result of the Strang Splitting
in order to compare the output here to the results of WaveBlocks

The class is templated on the dimension D of the simulation
*/

template<unsigned int D>
class Hdf5Writer
{
public:

	/**
	Construct a new writer by passing the filename
	*/
	Hdf5Writer(std::string filename) : filename_(filename), file_(filename_,H5F_ACC_TRUNC) {
		//H5File file(FILE_NAME, H5F_ACC_TRUNC);
	}

	void tester() {
		std::cout << "Class is working, dimension = " << D << std::endl;
		//std::cout << filename_ << std::endl;
	}

private:
	H5std_string filename_;
	H5File file_;
};
