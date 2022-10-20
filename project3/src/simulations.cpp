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
        trap.evolve_RK4(dt);

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

    ofile.close();

    return 0;
}

int two_parts(int n, double t, arma::vec3 ipos1, arma::vec3 ivel1, 
                arma::vec3 ipos2, arma::vec3 ivel2, bool interaction)
{
    PenningTrap trap;
    Particle p1(ipos1, ivel1), p2(ipos2, ivel2);

    trap.add_particle(p1);
    trap.add_particle(p2);

    double dt = t/n;

    std::vector<std::ofstream> ofile(2);

    for (int i = 0; i < 2; i++)
    {
        ofile[i].open("data/twoparts" + std::to_string(i+1) + ".txt");
    }

    int prec = 8;
    int width = prec + 10;

    for (int i = 0; i < 2; i++)
    {
        ofile[i] << std::setw(width) << std::setprecision(prec) 
                    << std::scientific << 0.;
        for (int j = 0; j < 3; j++)
        {
            ofile[i] << std::setw(width) << std::setprecision(prec) 
                        << std::scientific << trap.parts[i].pos()(j);
            ofile[i] << std::setw(width) << std::setprecision(prec) 
                        << std::scientific << trap.parts[i].vel()(j);
        }
    }
    for (int i = 0; i < 2; i++)
    {
        ofile[i] << std::endl;
    }


    for (int m = 0; m < n; m++)
    {
        trap.evolve_RK4(dt, interaction);
        // trap.evolve_forward_Euler(dt, interaction);

        for (int i = 0; i < 2; i++)
        {
            ofile[i] << std::setw(width) << std::setprecision(prec) 
                        << std::scientific << (i+1)*dt;
            for (int j = 0; j < 3; j++)
            {
                ofile[i] << std::setw(width) << std::setprecision(prec) 
                            << std::scientific << trap.parts[i].pos()(j);
                ofile[i] << std::setw(width) << std::setprecision(prec) 
                            << std::scientific << trap.parts[i].vel()(j);
            }
            ofile[i] << std::endl;
        }
    }

    for (int i = 0; i < 2; i++)
    {
        ofile[i].close();
    }

    return 0;
}

