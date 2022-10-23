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
    // Vd = V0/d^2
    double Vd_;

    const double k = 1.3893533e5;

public:
    std::vector<Particle> parts;

    // Constructors
    PenningTrap();
    PenningTrap(double B0, double V0, double d);
    PenningTrap(double B0, double V0, double d, std::vector<Particle> parts);
    PenningTrap(std::vector<Particle> parts);

    PenningTrap(int n);

    // Electric and magnetic field at position 'r'
    const arma::vec3 E_field(arma::vec3 r);
    const arma::vec3 B_field(arma::vec3 r);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // Counts number of particles in the trap
    int num_parts_in_trap();

    // Force on particle_i from particle_j
    const arma::vec3 force_particle(int i, int j);

    // The total force on particle_i from the external fields
    const arma::vec3 total_force_external(int i);

    // The total force on particle_i from the other particles
    const arma::vec3 total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    const arma::vec3 total_force(int i, bool interaction);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_Euler(double dt, bool interaction = true);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt, bool interaction = true);

    // Run the system in 'trap' for time 't'
    int run_experiment(int n, double t, std::string filename, 
                        bool interaction = true, 
                        std::string method = "RK4");
    
    double omega_0(Particle p);
    double omega_z(Particle p);

    // The phases omega_+ and omega_-
    double omega_1(Particle p);
    double omega_2(Particle p);

    int exact_sol(int n, double t, Particle p, std::string filename);

};


#endif