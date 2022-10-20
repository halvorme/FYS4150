#include "simulations.hpp"

#include <armadillo>
#include <fstream>
#include <iomanip>

#include "PenningTrap.hpp"
#include "Particle.hpp"

int single_part(int n, double t, arma::vec3 init_pos, arma::vec3 init_vel)
{
    PenningTrap trap;
    Particle p(init_pos, init_vel);

    trap.add_particle(p);

    double dt = t/n;

    std::ofstream ofile;
	ofile.open("data/singlepart.txt");

    int prec = 8;
    int width = prec + 10;

    ofile << std::setw(width) << std::setprecision(prec) 
                << std::scientific << 0.;
    for (int j = 0; j < 3; j++)
    {
        ofile << std::setw(width) << std::setprecision(prec) 
                << std::scientific << trap.parts[0].pos()(j);
        ofile << std::setw(width) << std::setprecision(prec) 
                << std::scientific << trap.parts[0].vel()(j);
    }
    ofile << std::endl;

    for (int i = 0; i < n; i++)
    {
        trap.evolve_forward_Euler(dt);

        ofile << std::setw(width) << std::setprecision(prec) 
                << std::scientific << (i+1)*dt;
        for (int j = 0; j < 3; j++)
        {
            ofile << std::setw(width) << std::setprecision(prec) 
                    << std::scientific << trap.parts[0].pos()(j);
            ofile << std::setw(width) << std::setprecision(prec) 
                    << std::scientific << trap.parts[0].vel()(j);
        }
        ofile << std::endl;
    }
    
    return 0;
}