#include "algo.hpp"
#include "utils.hpp"

int main()
{
	std::string folder = "data/";
	int prec = 15;
	int n_steps = 1000;
	double x_min = 0.;
	double x_max = 1.;

	std::string filename;


	filename = folder + "u_" + std::to_string(n_steps) + ".txt";
	arma::mat xu = uExact(n_steps, x_min, x_max, prec);
	toFile(xu, filename, prec);

	n_steps = 10;
	for (int i = 1; i < 5; i++)
	{
		n_steps = std::pow(10,i);
		filename = folder + "x_v_" + std::to_string(n_steps) + ".txt";

		solvePoissonGen(n_steps, x_min, x_max, filename, prec);
	}
	return 0;
}
