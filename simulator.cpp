#include "strang.hpp"

/* Schrödinger Equation Solver */

int main() {

	const auto eps = 0.01;
	const auto ieps = 1./eps;
	const auto tend = 5.0*std::sqrt(ieps);
	const auto n_timesteps = tend*10;
	const auto delta_t = tend/n_timesteps; // delta_t is h in python
		

	const auto N = 16;

	StrangSplitter<N> system(eps, tend, delta_t);
	system.Split();

	// store initial configuration
	//writer.StoreResult();

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