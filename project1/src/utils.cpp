#include <armadillo> // Compile with flag -larmadillo
#include <iomanip>
#include <iostream>
#include <fstream>

#include "utils.hpp"

// Write two vectors to file
int toFile(arma::vec x, arma::vec y, std::string filename, int prec) 
{
	int n = arma::size(x)[0];

	// Open file
	std::ofstream ofile;
	ofile.open(filename);

	// Set formatting
	int width = prec + 8;

	//Writes initial values to file
	ofile << std::setw(width) << std::setprecision(prec) << std::scientific << 0.
			  << std::setw(width) << std::setprecision(prec) << std::scientific << 0.
			  << std::endl;
	// Write vectors to file in scientific format
	for (int i = 0; i < n; i++){
		ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x(i)
			  << std::setw(width) << std::setprecision(prec) << std::scientific << y(i)
			  << std::endl;
	}
	// Write final values to file
	ofile << std::setw(width) << std::setprecision(prec) << std::scientific << 1.
			  << std::setw(width) << std::setprecision(prec) << std::scientific << 0.
			  << std::endl;


	// Close file
	ofile.close();
	return 0;
}

// Initialise an x-axis (not including endpoints)
arma::vec xAxis(int n_steps, double x_min, double h)
{
	int n = n_steps - 1;

	arma::vec x = arma::vec(n);
	x(0) = x_min + h;
	for (int i = 1; i < n; i++) {
		x(i) = x(i-1) + h;
	}
	return x;
}
