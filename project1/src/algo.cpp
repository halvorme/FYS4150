#include <armadillo>  // compile with flag -larmadillo
#include <cmath>
#include <string>

#include "utils.hpp"
#include "algo.hpp"

// Support functions for uExact
double uFunc(double x);
arma::vec ux(arma::vec x);


// Support functions for solvePoissonGen
double f(double x);
arma::vec gx(double h, arma::vec x);

// Implementation of the general algorithm
arma::vec genAlgo(arma::vec a, arma::vec b, arma::vec c, arma::vec g);
// Implementation of the special algorithm
arma::vec specAlgo(arma::vec g);


// Gives the exact solution of the Poisson equation
// Inputs are the number of points we want to calculate(+1) and the endpoints of the x-axis
// Returns a matrix, xu, with coloumns x and u(x)
arma::mat uExact(int n_steps, double x_min, double x_max)
{
	int n = n_steps-1;
	double h = (x_max - x_min) / n_steps;

	// Initialise Armadillo matrix for x and u axis
	arma::mat xu = arma::mat(n,2);

	// Initialise the x-axis (does not include endpoints)
	xu.col(0) = xAxis(n_steps, x_min, h);

	xu.col(1) = ux(xu.col(0));

	return xu;
}


// Solves the Poisson equation numerically, using the general algorithm
// Inputs are the precision of the discretisation and the endpoints of the x-axis
// Returns a matrix, 'xv', with coloumns x and v(x) (The endpoints are not included in xv)
arma::mat solvePoissonGen(int n_steps, double x_min, double x_max)
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
	arma::vec g = gx(h, xv.col(0));

	// Solve the matrix equation for v
	xv.col(1) = genAlgo(a, b, c, g);

	// Write to file
	return xv;
}


// Solves the Poisson equation numerically, using the special algorithm
// Inputs are the precision of the discretisation and the endpoints of the x-axis
// Returns a matrix, 'xv', with coloumns x and v(x) (The endpoints are not included in xv)
arma::mat solvePoissonSpec(int n_steps, double x_min, double x_max)
{
	int n = n_steps - 1;
	double h = (x_max - x_min)/n_steps;

	arma::mat xv = arma::mat(n,2);

	// Initialise the x-axis (does not include endpoints)
	xv.col(0) = xAxis(n_steps, x_min, h);

	// Initialise g(x)
	arma::vec g = gx(h, xv.col(0));

	// Solve the matrix equation for v
	xv.col(1) = specAlgo(g);

	// Write to file
	return xv;
}


// Performs the general algorithm for solving the matrix equation Av=g, where A is tridiagonal
// Inputs are the three diagonals, 'a', 'b' and 'c', defining A and the modified source term 'g'
// Returns the solution 'v'
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


// Performs the special algorithm to solve the Poisson equation
// Input is the modified source term 'g'
// Returns the numerical solution 'v'
arma::vec specAlgo(arma::vec g)
{
	int n = arma::size(g)[0];

	arma::vec v = arma::vec(n);
	arma::vec b = arma::vec(n);

	// Forward substitution
	b(0) = 2.;
	for (int i = 0; i < n-1; i++) {
		b(i+1) = 2. - 1./b(i);
		g(i+1) += g(i)/b(i);
	}

	// Backward substitution
	v(n-1) = g(n-1)/b(n-1);

	for (int i = 2; i <= n; i++) {
		v(n-i) = (g(n-i) + v(n-i+1))/b(n-i);
	}

	return v;
}


// Calculates the logarithm of the absolute and relative error of v(x)
// Inputs are an x-axis, the exact solution u(x) and the approximate solution v(x)
// Return a matrix with the coloumns x, absolute error and relative error. 
arma::mat numError(arma::vec x, arma::vec v)
{
	int n = arma::size(x)[0];
	arma::mat err(n,3);
	err.col(0) = x;

	arma::vec u = ux(x);

	for (int i = 0; i < n; i++)
	{
		err(i,1) = std::log10(std::abs(u(i)-v(i)));
		err(i,2) = std::log10(std::abs((u(i)-v(i))/u(i)));
	}
	return err;
}


// Tests the running time of the general and special algorithm
// Inputs are the number of steps of the most fine-grained x-axis (10^max_steps steps), number of 
// repetitions for each choice of n_steps.
// Results are printed to file 'filename', with precision 'prec'
int timeAlgo(int max_steps, int reps, std::string filename, int prec, double x_min, double x_max)
{
	int n_steps;
	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();
	double gen_time;
	double spec_time;
	double frac;
	// Open file
	std::ofstream timer;
	timer.open("data/timing.txt");

	for (int i = 1; i <= max_steps; i++)
	{
		n_steps = std::pow(10,i);
		// Time general algorithm
		t1 = std::chrono::high_resolution_clock::now();
		for (int j = 0; j < reps; j++)
		{
			solvePoissonGen(n_steps,x_min,x_max);
		}
		t2 = std::chrono::high_resolution_clock::now();
		gen_time = std::chrono::duration<double>(t2 - t1).count();

		// Time special algorithm
		t1 = std::chrono::high_resolution_clock::now();
		for (int j = 0; j < reps; j++)
		{
			solvePoissonSpec(n_steps,x_min,x_max);
		}
		t2 = std::chrono::high_resolution_clock::now();
		spec_time = std::chrono::duration<double>(t2 - t1).count();

		frac = gen_time/spec_time;

		// Write results to file "timing.txt"
		timer << std::setw(prec + 8) << std::setprecision(prec) << std::scientific << gen_time
			  << std::setw(prec + 8) << std::setprecision(prec) << std::scientific << spec_time
			  << std::setw(prec + 8) << std::setprecision(prec) << std::scientific << frac << std::endl;
	}

	timer.close();

	return 0;
}


// Returns the exact solution to the equation 'uFunc' on the axis 'x' 
arma::vec ux(arma::vec x)
{
	int n = arma::size(x)[0];
	arma::vec u = arma::vec(n);
	for (int i = 0; i < n; i++) {
		u(i) = uFunc(x(i));
	}
	return u;
}


// Returns the value of 'u' at the point 'x'
double uFunc(double x){
	return 1. - (1.-std::exp(-10))*x - std::exp(-10.*x);
}


// Initialises the modified source term 'g'
// Depends on the source term and stepsize 'h'
arma::vec gx(double h, arma::vec x)
{
	int n = arma::size(x)[0];
	arma::vec g = arma::vec(n);
	for (int i = 0; i < n; i++) {
		g(i) = std::pow(h,2)*f(x(i));
	}
	return g;
}


// Defines the source term
double f(double x)
{
	return 100*std::exp(-10.*x);
}
