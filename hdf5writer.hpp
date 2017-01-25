#ifndef STRANG_SPLITTING_HDF5_WRITER
#define STRANG_SPLITTING_HDF5_WRITER
#include <string>
#include <vector>
#include <memory>
#include <iostream> /* Fixme remove me */

#include <Eigen/Core>
#include "H5Cpp.h" // CPP header for HDF interface

#include "strang.hpp"

using namespace H5;

/**
Main functionality is to store the result of the Strang Splitting
in order to compare the output here to the results of WaveBlocks

Fixme (this is not so correct): The class is templated on the dimension N of the simulation
*/
template<dim_t N>
class Hdf5Writer
{
public:

	/**
	Construct a new writer by passing the filename
	*/
	Hdf5Writer(std::string filename, const StrangSplitter<N>& system) 
		: 
		filename_(filename), 
		file_(filename_,H5F_ACC_TRUNC),
		hdf_complex_type_(sizeof(ctype)) 
	{
		// set up complex datatype for hdf
		hdf_complex_type_.insertMember("r", HOFFSET(ctype, real), PredType::NATIVE_DOUBLE);
		hdf_complex_type_.insertMember("i", HOFFSET(ctype, imag), PredType::NATIVE_DOUBLE);

		Setup(system);
	}

	void tester() {
		std::cout << "Class is working, dimension = " << N << std::endl;
		//std::cout << filename_ << std::endl;
	}

	/**
	Post:	Extends the HDF file with the results of the current timestep
	*/
	void StoreResult(const StrangSplitter<N>& system){

	}

	/**
	Post:	After the simulation is done and all results are collected
			this function will clean up the created hdf5 file and make it
			as small as possible.
	*/
	void Finalize() {

	}

	/**
	Post:	Changing the timestep size to another value than 1 (which is the default)
			will store results only at every t'th timestep
	*/
	void SetTimestepSize(const unsigned int t) {
		timestep_size_ = t;
	}

	/**
	The HDF5 library needs access to the members of a type
	so we cannot use std::complex
	*/
	struct ctype {
		double real=0.;
		double imag=0.;
	};
private:
	/**
	Pre: 	Constructor was already called. This function should not be
			called somewhere else
	Post: 	HDF5 is created with all the correct metadata in place
	*/
	void Setup(const StrangSplitter<N>& system) {
		// see prestructuring
		// define group structure
		// define dimension
		// 

		DefineGroupStructure();
		WriteSimulationParameters(system);
		WriteGrid(system); 
	}

	void DefineGroupStructure() {
		g_parameter_ = std::make_shared<Group>(file_.createGroup(path_parameters_));
		g_grid_ = std::make_shared<Group>(file_.createGroup(path_grid_));
	}

	/**
	Post:	Write the simulation parameters as attributes into the file
	*/
	void WriteSimulationParameters(const StrangSplitter<N>& system) {
		constexpr dim_t nx = N;
		constexpr dim_t ny = 1;

		hsize_t size_parameters[] = {1};	// Parameters are always just 1 element
		DataSpace s_parameters(RANK1_, size_parameters);

		// Create Attributes
		a_dt_ = g_parameter_->createAttribute("dt", PredType::NATIVE_DOUBLE, s_parameters);
		a_eps_ = g_parameter_->createAttribute("epsilon", PredType::NATIVE_DOUBLE, s_parameters);
		a_scale_ = g_parameter_->createAttribute("scale", PredType::NATIVE_DOUBLE, s_parameters);
		a_n_timesteps_ = g_parameter_->createAttribute("n_timesteps", PredType::NATIVE_UINT, s_parameters);
		// dimensions of the grid, they go into the grid group
		a_nx_ = g_grid_->createAttribute("Nx", PredType::NATIVE_UINT, s_parameters); 
		a_ny_ = g_grid_->createAttribute("Ny", PredType::NATIVE_UINT, s_parameters);

		auto dt = system.GetDt();
		auto eps = system.GetEpsilon();
		auto scale = system.GetGridScale();
		auto nts = system.GetNumberOfTimeSteps();
		
		// write() requires the pointer to the input data
		a_dt_.write(PredType::NATIVE_DOUBLE, &dt);
		a_eps_.write(PredType::NATIVE_DOUBLE, &eps);
		a_scale_.write(PredType::NATIVE_DOUBLE, &scale);
		a_n_timesteps_.write(PredType::NATIVE_UINT, &nts);
		a_nx_.write(PredType::NATIVE_UINT, &nx);
		a_ny_.write(PredType::NATIVE_UINT, &ny);
	}

	/**
	Post: 	Write the spacial grid into the file
	*/
	void WriteGrid(const StrangSplitter<N>& system) {
		
		hsize_t size_parameters[] = {1};	// Parameters are always just 1 element
		DataSpace s_parameters(RANK1_, size_parameters);

		hsize_t size_grid[] = {N};
		DataSpace s_grid(RANK1_, size_grid);
		const auto& grid = system.GetGrid();

		double transformed_grid[N];
		transform(transformed_grid, grid);

		std::shared_ptr<DataSet> d_grid;
		d_grid = std::make_shared<DataSet>(file_.createDataSet(path_grid_+"/Computational_Grid", PredType::NATIVE_DOUBLE, s_grid));	

		d_grid->write(transformed_grid, PredType::NATIVE_DOUBLE);
	}

	/**
	Post:	Transforms the Eigen::Matrix into a native array of the same type.
			(Obviously this does not work for std::complex, there a different approach
			is necessary.)
	*/
	template <class T>
	void transform(T transformed_data[], const Eigen::Matrix<T,1,N>& mat) const {
		//transformed_data.resize(N);
		for(size_t i = 0; i < N; ++i)
			transformed_data[i] = mat[i];
	}

private:
	H5std_string filename_;
	H5File file_;
	CompType hdf_complex_type_;
	unsigned int timestep_size_ = 1; //how often do we print results. 1 means every timestep

	/* The Groups of the HDF file*/
	std::shared_ptr<Group> g_parameter_;// group that holds the attributes which again hold the simulation parameters
	std::shared_ptr<Group> g_grid_;		// space grid
	std::shared_ptr<Group> g_wave_;		// solution of the simulation

	/* Attributes contain the simulation parameters*/
	Attribute a_dt_;
	Attribute a_n_timesteps_;
	Attribute a_eps_;
	Attribute a_scale_;
	Attribute a_nx_;
	Attribute a_ny_;


	/* DataSpace */
	const int RANK1_=1;
    const int RANK2_=2;
    const int RANK3_=3;

	/* These strings describe the paths inside the HDF file */
	const H5std_string path_parameters_="Simulation_Parameters";
	const H5std_string path_grid_ = "Grid";
};

#endif /* STRANG_SPLITTING_HDF5_WRITER */