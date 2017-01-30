#ifndef STRANG_SPLITTING_HDF5_WRITER
#define STRANG_SPLITTING_HDF5_WRITER
#include <string>
#include <vector>
#include <memory>

#include <Eigen/Core>
#include "H5Cpp.h" // CPP header for HDF5 interface

#include "strang.hpp"

using namespace H5;

/*If we want N*N matrices, remove this and replace each occurence 
 with the template parameter N*/
const int DIMENSION_Y = 1;

/**
Main functionality is to store the result of the Strang Splitting
in order to compare the output here to the results of WaveBlocks
*/
template<dim_t N>
class StrangWriter
{
public:

	/**
	Construct a new writer by passing the filename and the system we want to store

	Post:	hdf5 compatible complex datatype and file structure ready to go.
	*/
	StrangWriter(std::string filename, const StrangSplitter<N>& system) 
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
	~StrangWriter() {
		Finalize();
	}

	/*
	Post:	Extends the HDF file with the results of the current timestep
	*/
	void StoreResult(const StrangSplitter<N>& system){

		store_counter_++;	// count how often this function was called

		// If timestep_size is 1, we write every timestep to the file
		// otherwise we do it only every n'th step.
		if (store_counter_ % timestep_size_ == 0) {
			std::cout << "Storing step " << store_counter_<< std::endl;

			SelectWritingSpace();

			// Eigen Objects can't be written straight to HDF
			// so we need to transform it first. 
			std::vector<ctype> transformed_data;
			const auto& u = system.GetSolution();
			transform(transformed_data, u);

			u_ds_->write(transformed_data.data(), hdf_complex_type_, u_basic_space_, u_space_);

			// Set up the Writing space for the next iteration
			AdvanceWriter();
		}
		
	}

	

	/**
	Post:	Changing the timestep size to another value than 1 (which is the default)
			will store results only at every t'th timestep
	*/
	void StoreAfterNSteps(const unsigned int t) {
		timestep_size_ = t;
	}

	
private:
	/**
	The HDF5 library needs access to the members of a type
	so we cannot use std::complex
	*/
	struct ctype {
		double real=0.;
		double imag=0.;
	} instanceof;

	/**
	Pre: 	Constructor was already called. This function should not be
			called somewhere else
	Post: 	HDF5 is created with all the correct metadata in place
	*/
	void Setup(const StrangSplitter<N>& system) {
		DefineGroupStructure();
		// For the simple static data, ie. parameters and grid
		WriteSimulationParameters(system);
		WriteGrid(system); 
		// For the more complicated, growing data, ie. the solution.
		SetChunks();
		SetDataSpaces();
		AllocateDataSets();	
		SetBasicSpace();
		SelectBasicHyperslabs();
	}

	void DefineGroupStructure() {
		g_parameter_ = std::make_shared<Group>(file_.createGroup(path_parameters_));
		g_grid_ = std::make_shared<Group>(file_.createGroup(path_grid_));
		g_solution_ = std::make_shared<Group>(file_.createGroup(path_solution_));

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
		a_dt_ = g_parameter_->createAttribute(attr_name_dt_, PredType::NATIVE_DOUBLE, s_parameters);
		a_eps_ = g_parameter_->createAttribute(attr_name_eps_, PredType::NATIVE_DOUBLE, s_parameters);
		a_scale_ = g_parameter_->createAttribute(attr_name_scale_, PredType::NATIVE_DOUBLE, s_parameters);
		a_n_timesteps_ = g_parameter_->createAttribute(attr_name_n_timesteps_, PredType::NATIVE_UINT, s_parameters);
		// dimensions of the grid, they go into the grid group
		a_nx_ = g_grid_->createAttribute(attr_name_nx_, PredType::NATIVE_UINT, s_parameters); 
		a_ny_ = g_grid_->createAttribute(attr_name_ny_, PredType::NATIVE_UINT, s_parameters);

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
		
		hsize_t size_grid[] = {N};
		DataSpace s_grid(RANK1_, size_grid);
		const auto& grid = system.GetGrid();

		real_t transformed_grid[N];
		transform(transformed_grid, grid);

		std::shared_ptr<DataSet> d_grid;
		d_grid = std::make_shared<DataSet>(file_.createDataSet(path_grid_+data_name_grid_, PredType::NATIVE_DOUBLE, s_grid));	

		d_grid->write(transformed_grid, PredType::NATIVE_DOUBLE);
	}

	/**
	Post:	Property list describing dimensionality of chunked data
	*/
	void SetChunks() {
		hsize_t chunk_dim[] = {1,N,DIMENSION_Y};
		propertylist_u_.setChunk(RANK3_, chunk_dim);
		propertylist_u_.setFillValue(hdf_complex_type_, &instanceof);

	}
	/**
	Post: 	Initial DataSpace is set with the unlimited tag
	*/
	void SetDataSpaces() {
		const hsize_t dim[] = {1,N,DIMENSION_Y};
		const hsize_t maxdim[] = {H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};
		DataSpace s_u(RANK3_, dim, maxdim); 
		u_space_ = s_u;
	}
	/**
	Post:	Initial DataSet is created within the file
	*/
	void AllocateDataSets() {
		u_ds_ = std::make_shared<DataSet>(
			file_.createDataSet(path_solution_+data_name_u_, hdf_complex_type_, u_space_, propertylist_u_));
	}
	/**
	This only needs to be called once because the source space never changes
	Post:	Define the source space, this layout never changes
	*/
	void SetBasicSpace() {
		u_basic_dim_[0] = 1;
		u_basic_dim_[1] = N;
		u_basic_dim_[2] = DIMENSION_Y;
		DataSpace s_temp(RANK3_, u_basic_dim_);
		u_basic_space_ = s_temp;
	}
	/**
	This only needs to be called once because the source space never changes
	Post:	Source space is selected
	*/
	void SelectBasicHyperslabs() {
		hsize_t count[]={1,N,DIMENSION_Y};
		hsize_t start[]={0,0,0};
		hsize_t stride[]={1,1,1};
		hsize_t block[]={1,1,1};

		u_basic_space_.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	}
	/**
	Post:	Dynamic writing space is selected
	*/
	void SelectWritingSpace() {
		hsize_t count[]={1,N,DIMENSION_Y};
		
		hsize_t start[3];
		start[0] = timestep_index_-1;
		start[1] = 0;
		start[2] = 0;
		
		hsize_t stride[]={1,1,1};
		hsize_t block[]={1,1,1};

		u_space_.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	}
	/**
	Post:	HDF5 file is now ready write the next timestep
	*/
	void AdvanceWriter() {
		timestep_index_+=1;
		ExtendWritingSpace();
		ExtendDataSet();
		UpdateFileSpace();

	}
	/**
	Post:	Dimensions for the new DataSpace are set
	*/
	void ExtendWritingSpace() {
		extension_dim[0] = timestep_index_;
		extension_dim[1] = N;
		extension_dim[2] = DIMENSION_Y;
	}
	/**
	Post:	Dimension for the new DataSpace are now extended
	*/
	void ExtendDataSet() {
		u_ds_->extend(extension_dim);
	}
	/**
	Post:	Recover the newly created DataSpace for the current DataSet	
	*/
	void UpdateFileSpace() {
		u_space_=u_ds_->getSpace();
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

	/**
	Post:	content of complex Eigen::Matrix is stored in a vector of the ctype 
			which in turn is used to fill the compound datatype hdf_complex_type_
	*/
	void transform(std::vector<ctype>& transformed_data, const CMatrix<1,N>& mat) const {
		transformed_data.resize(N);
		for(int i = 0; i < N; ++i) {
			transformed_data[i].real = mat[i].real();
			transformed_data[i].imag = mat[i].imag();
		}
	}
	/**
	Post:	After the simulation is done and all results are collected
			this function will clean up the created hdf5 file and make it
			as small as possible.
	*/
	void Finalize() {
		timestep_index_ -=1;
		ExtendWritingSpace();
		ExtendDataSet();
	}

private:
	H5std_string filename_;
	H5File file_;
	CompType hdf_complex_type_;
	unsigned int store_counter_ = 0; //counts how often the store function has been called
	
	unsigned int timestep_size_ = 1; //how often do we print results. 1 means every timestep

	hsize_t timestep_index_ = 1; // index of the 
	

	/* The Groups of the HDF file*/
	std::shared_ptr<Group> g_parameter_;// group that holds the attributes which again hold the simulation parameters
	std::shared_ptr<Group> g_grid_;		// space grid
	std::shared_ptr<Group> g_solution_;		// solution of the simulation

	/* DataSpace */
	DataSpace u_space_;
	DataSpace u_basic_space_;
	hsize_t u_basic_dim_[3];

	hsize_t extension_dim[3];

	/* DataSet */
	std::shared_ptr<DataSet> u_ds_;

	/* Properties for DataSets*/
	DSetCreatPropList propertylist_u_;

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
	const H5std_string path_parameters_			= "Simulation_Parameters";
	const H5std_string path_grid_ 				= "Grid";
	const H5std_string path_solution_ 			= "Solution";

	const H5std_string attr_name_dt_ 			= "dt";
	const H5std_string attr_name_eps_ 			= "epsilon";
	const H5std_string attr_name_scale_ 		= "scale";
	const H5std_string attr_name_n_timesteps_ 	= "n_timesteps";
	const H5std_string attr_name_nx_			= "Nx";
	const H5std_string attr_name_ny_			= "Ny";

	const H5std_string data_name_grid_			= "/Computational_Grid";
	const H5std_string data_name_u_				= "/u";


	
};

#endif /* STRANG_SPLITTING_HDF5_WRITER */