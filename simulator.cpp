#include "strang.hpp"
#include "hdf5writer.hpp"

/* Schrödinger Equation Solver */

int main() {

	const auto eps = 0.01;
	const auto ieps = 1./eps;
	const auto tend = 5.0*std::sqrt(ieps);
	const auto n_timesteps = tend*10;
	const auto dt = tend/n_timesteps; // dt is h in python
		

	const auto N = 16;

	StrangSplitter<N> system(eps, tend, dt);
	//system.Split();


	// store initial configuration
	//writer.StoreResult();
	Hdf5Writer<N> writer("strang_1D_cpp.hdf5", system);
	writer.StoreResult(system);
	//writer.tester();

	for(real_t t = dt; t < tend; t+=dt) {
		system.Advance();
		writer.StoreResult(system);
	}

	writer.Finalize();
	std::cout << system.Norm() << std::endl;

	/*
	for(timesteps){
		strang_splitting();
		writer.StoreResult();
	}
	writer.Finalize();
	*/


	//run();

	return 0;
}