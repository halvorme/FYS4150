#include "algo.hpp"
#include "utils.hpp"

int main()
{
	std::string folder = "data/";
	double x_min = 0.;
	double x_max = 1.;

	int prec = 15;
	std::string filename;

	// Problem 2
	int n_steps = 1000;
	filename = folder + "u_" + std::to_string(n_steps) + ".txt";

	arma::mat xu = uExact(n_steps, x_min, x_max);

	matToFile(xu, filename, prec);


	// Problem 7 and 8

	// Set maximal number of steps (10^max_steps)
	int max_steps = 5;

	// Creates file to save maximum relative error
	std::ofstream max_error;
	max_error.open(folder + "max_error.txt");

	// Solves the Poisson equation numerically for variyng numbers of steps. 
	// Prints both results and errors to file
	for (int i = 1; i <= max_steps; i++)
	{
		n_steps = std::pow(10,i);
		filename = folder + "x_v_" + std::to_string(i) + ".txt";

		// Solves Poisson equation
		arma::mat xv = solvePoissonGen(n_steps, x_min, x_max);
		matToFile(xv, filename, prec);

		filename = folder + "err_" + std::to_string(i) + ".txt";

		// Calculates error
		arma::mat err = numError(xv.col(0), xv.col(1));
		matToFile(err, filename, prec);

		// Finds maximal error
		max_error << std::setprecision(3) << std::scientific << std::exp(err.col(2).max()) << std::endl;
	}

	max_error.close();


	// Problem 10

	// Set maximal number of steps (10^max_steps) and how many repetitions for each choich of n_steps.
	max_steps = 5;
	int reps = 100;
	
	// Output options
	filename = folder + "timing.txt";
	prec = 3;

	// Results of timing are printed to file "timing.txt"
	timeAlgo(max_steps, reps, filename, prec, x_min, x_max);

	return 0;
}
