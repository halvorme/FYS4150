#ifndef __utils_hpp__
#define __utils_hpp__

#include <armadillo> // Compile with flag -larmadillo
#include <iomanip>
#include <iostream>
#include <fstream>


// Write two vectors to file
int matToFile(arma::mat xy, std::string filename, int prec); 

// Initialise an x-axis (not including endpoints)
arma::vec xAxis(int n_steps, double x_min, double h);

#endif