#include "algo.hpp"
#include "utils.hpp"

int main()
{
	// Set maximal number of steps (10^max_steps) for problem 7-8 and 10
	int max_steps78 = 5;
	int max_steps10 = 5;


	// Sets the boundaries of the x-axis
	double x_min = 0.;
	double x_max = 1.;

	// Formats the output
	int prec = 14;
	
	std::string folder = "data/";
	std::string filename;

	// Problem 2
	int n_steps = 1000;
	filename = folder + "u_" + std::to_string(n_steps) + ".txt";

	arma::mat xu = uExact(n_steps, x_min, x_max);

	matToFile(xu, filename, prec);


	// Problem 7 and 8

	// Creates file to save maximum relative error
	std::ofstream max_error;
	max_error.open(folder + "max_error.txt");

	// Solves the Poisson equation numerically for variyng numbers of steps. 
	// Prints both results and errors to file
	for (int i = 1; i <= max_steps78; i++)
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

	// Set the number of repetitions of the algorithm for each choice of n_steps.
	int reps = 100;

	// Output options
	filename = folder + "timing.txt";
	prec = 3;

	// Results of timing are printed to file "timing.txt"
	timeAlgo(max_steps10, reps, filename, prec, x_min, x_max);

	return 0;
}
