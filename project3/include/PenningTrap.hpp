#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <vector>
#include <armadillo>
#include <string>

#include "Particle.hpp"

class PenningTrap
{
private:
    double B0_, V0_, d_;
    double Vd_;

    double k = 1.3893533e5;

public:
    std::vector<Particle> parts;

    // Constructors
    PenningTrap();
    PenningTrap(double B0, double V0, double d);
    PenningTrap(double B0, double V0, double d, std::vector<Particle> parts);
    PenningTrap(std::vector<Particle> parts);

    // Electric and magnetic field at position 'r'
    const arma::vec3 E_field(arma::vec3 r);
    const arma::vec3 B_field(arma::vec3 r);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // Force on particle_i from particle_j
    const arma::vec3 force_particle(int i, int j);

    // The total force on particle_i from the external fields
    const arma::vec3 total_force_external(int i);

    // The total force on particle_i from the other particles
    const arma::vec3 total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    const arma::vec3 total_force(int i, bool interaction);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt, bool interaction = true);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt, bool interaction = true);

    // Run the system in 'trap' for time 't'
    void runExperiment(int n, double t, std::string filename, 
                        bool interaction = true);

};


#endif