#ifndef __algo_hpp__
#define __algo_hpp__

#include <armadillo>  // compile with flag -larmadillo
#include <cmath>

// Calculates the exact solution u(x)
int uExact(int n_steps, double x_min, double x_max, std::string filename, int prec);


// Solves the Poisson equation and prints the solution to file
int solvePoissonGen(int n_steps, double x_min, double x_max, std::string filename, int prec);

#endif