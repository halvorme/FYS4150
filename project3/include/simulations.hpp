#ifndef __simulations_hpp__
#define __simulations_hpp__

#include <armadillo>

#include "PenningTrap.hpp"
#include "Particle.hpp"

int single_part(int n, double t);

int two_parts(int n, double t);

int error_sim(std::vector<int> n, double t);

int single_exact(std::vector<int> n, double t);

int trapped_broad(int n, double t);

int resonance(int n, double t, double omega_low, double omega_up, bool interction);

int search(int n, double t, int n_parts, double omega_low, double omega_up, 
            double d_omega, std::vector<double> f, bool interaction, 
            std::ofstream& ofile);

int trapped_broad_new(int n, double t);

#endif