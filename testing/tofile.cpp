#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

int main(){
	 // Set a filename
	 std::string filename = "x_y.txt";
	
	// Create and open the output file. Or, technically, create 
	// an "output file stream" (type std::ofstream) and connect it to our filename.
	std::ofstream ofile;
	ofile.open(filename);
	
	// Set some parameters for our computation
	double x_min = 0.0;
	double x_max = 1.0;
	int n_steps = 100;
	double h = (x_max - x_min) / n_steps;
	
	// Some width and precision parameters we will use to format the output
	int width = 12;
	int prec  = 4;
	
	// Set the first x and y values
	double x = x_min;
	double y = x*x;
	
	// Loop over steps
	for (int i = 0; i <= n_steps; i++)
	{
		// Write a line with the current x and y values (nicely formatted) to file
		ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x
			<< std::setw(width) << std::setprecision(prec) << std::scientific << y
			<< std::endl;
		
		// Update x and y values
		x += h;
		y = x*x;
	}  
	
	// Close the output file
	ofile.close();
	
	// All is well. Exit program with return code 0.
	return 0;
}
