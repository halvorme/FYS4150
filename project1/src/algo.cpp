#include <armadillo>  // compile with flag -larmadillo
#include <cmath>
#include <string>

#include "utils.hpp"
#include "algo.hpp"

// Support functions for uExact
double uFunc(double x);
arma::vec ux(int n, arma::vec x);


// The following are support functions for solvePoissonGen
double f(double x);
arma::vec gx(int n, double h, arma::vec x);
// The implementation of the actual general algorithm
arma::vec genAlgo(arma::vec a, arma::vec b, arma::vec c, arma::vec g);


arma::mat uExact(int n_steps, double x_min, double x_max, int prec)
{
    int n = n_steps-1;
    double h = (x_max - x_min) / n_steps;

	// Initialise Armadillo matrix for x and u axis
	arma::mat xu = arma::mat(n,2);

	// Initialise the x-axis (does not include endpoints)
	xu.col(0) = xAxis(n_steps, x_min, h);

    xu.col(1) = ux(n, xu.col(0));

    return xu;
}


// Solves the Poisson equation and prints the solution to file
int solvePoissonGen(int n_steps, double x_min, double x_max, std::string filename, int prec)
{
	int n = n_steps - 1;
	double h = (x_max - x_min)/n_steps;

	arma::mat xv = arma::mat(n,2);

	// Initialise the x-axis (does not include endpoints)
	xv.col(0) = xAxis(n_steps, x_min, h);

	// Initialise tridiagonal matrix 
	arma::vec a = arma::vec(n-1).fill(-1.);
	arma::vec c = arma::vec(n-1).fill(-1.);
	arma::vec b = arma::vec(n).fill(2.);

	// Initialise g(x)
	arma::vec g = gx(n, h, xv.col(0));

	// Solve the matrix equation for v
	xv.col(1) = genAlgo(a, b, c, g);

	// Write to file
	toFile(xv, filename, prec);

	return 0;
}


arma::vec ux(int n, arma::vec x)
{
	arma::vec u = arma::vec(n);
	for (int i = 0; i < n; i++) {
		u(i) = uFunc(x(i));
	}
	return u;
}


double uFunc(double x){
    return 1. - (1.-std::exp(-10))*x - std::exp(-10.*x);
}


// Initialise g(x)
arma::vec gx(int n, double h, arma::vec x)
{
	arma::vec g = arma::vec(n);
	for (int i = 0; i < n; i++) {
		g(i) = std::pow(h,2)*f(x(i));
	}
	return g;
}

// Performs the general algorithm for solving the matrix equation Av=g, 
// where A is tridiagonal
arma::vec genAlgo(arma::vec a, arma::vec b, arma::vec c, arma::vec g)
{
	int n = arma::size(g)[0];

	arma::vec v = arma::vec(n);

	// Forward substitution
	double d;
	for (int i = 0; i < n-1; i++) {
		d = a(i)/b(i);
		b(i+1) -= d*c(i);
		g(i+1) -= d*g(i);
	}

	// Backward substitution
	v(n-1) = g(n-1)/b(n-1);

	for (int i = 2; i <= n; i++) {
		v(n-i) = (g(n-i)-c(n-i)*v(n-i+1))/b(n-i);
	}

	return v;
}


double f(double x)
{
    return 100*std::exp(-10.*x);
}
