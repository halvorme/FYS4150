#include <string>
#include <cmath>
#include <armadillo> // compile with flag -larmadillo

#include "utils.hpp"

double func(double x);
arma::vec ux(int n, double h, arma::vec x);
int uExact(int n_steps, double x_min, double x_max, std::string filename, int prec);

int main()
{
    std::string filename;
    int n_steps = 10000;
    int prec = 8;

    double x_min = 0.;
    double x_max = 1.;

	filename = "u_" + std::to_string(n_steps) + ".txt";
	uExact(n_steps, x_min, x_max, filename, prec);

    return 0;
}

int uExact(int n_steps, double x_min, double x_max, std::string filename, int prec)
{
    int n = n_steps-1;
    double h = (x_max - x_min) / n_steps;

	// Initialise the x-axis (does not include endpoints)
	arma::vec x = xAxis(n_steps, x_min, h);

    arma::vec u = ux(n,h,x);

    toFile(x, u, filename, prec);

    return 0;
}

arma::vec ux(int n, double h, arma::vec x)
{
	arma::vec u = arma::vec(n);
	for (int i = 0; i < n; i++) {
		u(i) = func(x(i));
	}
	return u;
}


double func(double x){
    return 1. - (1.-std::exp(-10))*x - std::exp(-10.*x);
}
