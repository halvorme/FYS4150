#include <armadillo> // Compile with flag -larmadillo
#include <iomanip>
#include <iostream>
#include <fstream>

#include "utils.hpp"

// Write Armadillo matrix to file
int matToFile(arma::mat xy, std::string filename, int prec) 
{
	int n = arma::size(xy)[0];
	int n_col = arma::size(xy)[1];

	// Open file
	std::ofstream ofile;
	ofile.open(filename);

	// Set formatting
	int width = prec + 10;

	//Writes initial values to file
	ofile << std::setw(width) << std::setprecision(prec) << std::scientific << 0.;
	for (int j = 0; j < n_col-1; j++)
	{
		ofile << std::setw(width) << std::setprecision(prec) << std::scientific << 0.;
	}
	ofile << std::endl;

	// Write vectors to file in scientific format
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n_col; j++)
		{
			ofile << std::setw(width) << std::setprecision(prec) << std::scientific << xy(i,j);
		}
		ofile << std::endl;
	}
	// Write final values to file
	ofile << std::setw(width) << std::setprecision(prec) << std::scientific << 1.;
	for (int j = 0; j < n_col-1; j++)
	{
		ofile << std::setw(width) << std::setprecision(prec) << std::scientific << 0.;
	}
	ofile << std::endl;

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
