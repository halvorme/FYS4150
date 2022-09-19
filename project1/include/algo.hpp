#ifndef __algo_hpp__
#define __algo_hpp__

#include <armadillo>	// compile with flag -larmadillo
#include <cmath>

// Gives the exact solution of the Poisson equation
// Inputs are the number of points we want to calculate(+1) and the endpoints of the x-axis
// Returns a matrix, xu, with coloumns x and u(x)
arma::mat uExact(int n_steps, double x_min, double x_max);


// Solves the Poisson equation numerically
// Inputs are the precision of the discretisation and the endpoints of the x axis
// Returns a matrix, xv, with coloumns x and v(x) (The endpoints are not included in xv)
arma::mat solvePoissonGen(int n_steps, double x_min, double x_max);


// Calculates the logarithm of the absolute and relative error of v(x)
// Inputs are an x-axis, the exact solution u(x) and the approximate solution v(x)
// Return a matrix with the coloumns x, absolute error and relative error. 
arma::mat numError(arma::vec x, arma::vec v);


// Tests the running time of the general and special algorithm
// Inputs are the number of steps of the most fine-grained x-axis (10^max_steps steps), number of 
// repetitions for each choice of n_steps.
// Results are printed to file 'filename', with precision 'prec'.
int timeAlgo(int max_steps, int reps, std::string filename,int prec, double x_min, double x_max);

#endif