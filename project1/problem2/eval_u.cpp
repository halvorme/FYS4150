#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>

#include <typeinfo>

double u(double x);

int main(){
    std::string filename = "x_u_exact.txt";

    std::ofstream ofile;
    ofile.open(filename);

    double x_min = 0.0;
    double x_max = 1.0;
    int n_steps = 100;
    double h = (x_max - x_min) / n_steps;

    int width = 11;
    int prec = 3;

    double x = x_min;
    double y = u(x);

    for (int i = 0; i <= n_steps; i++){
        // Write a line with the current x and y values (nicely formatted) to file
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x
              << std::setw(width) << std::setprecision(prec) << std::scientific << y
              << std::endl;

        // Update x and y values
        x += h;
        y = u(x);
    }

    ofile.close();

    return 0;
}

double u(double x){
    return 1. - (1.-std::exp(-10))*x - std::exp(-10.*x);
}
