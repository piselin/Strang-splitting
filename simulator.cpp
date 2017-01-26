#include "strang.hpp"
#include "hdf5writer.hpp"

/* 
* This is an example on how to use this
* Schrödinger Equation Solver. 
*/

int main() {

	const auto eps = 0.01;					// Problemdependent constant...
	const auto ieps = 1./eps;				// ... and its inverse
	const auto tend = 5.0*std::sqrt(ieps);	// T final
	const auto n_timesteps = tend*10;		// Total number of timesteps
	const auto dt = tend/n_timesteps;		// Size of a timestep
		
	const auto N = 16;

	// Setup the Solver
	StrangSplitter<N> system(eps, tend, dt);

	Hdf5Writer<N> writer("strang_1D_cpp.hdf5", system);
	writer.StoreResult(system); // store the initial configuration
	writer.StoreAfterNSteps(100);
	
	std::cout << "Solver will run through " << n_timesteps << " time steps" << std::endl;

	for(real_t t = dt; t < tend; t+=dt) {
		system.Advance();
		writer.StoreResult(system);
	}

	writer.Finalize();
	std::cout << "... and we are done." << std::endl;
	std::cout << "Norm(u) = " << system.Norm() << std::endl;

	return 0;
}