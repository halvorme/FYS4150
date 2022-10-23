#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <vector>
#include <armadillo>
#include <string>

#include "Particle.hpp"


// Defines the state of a Penning trap and interface to simulate 
// particle dynamics in the trap
class PenningTrap
{
public:
    // List of particles in the trap
    std::vector<Particle> parts;

    // Initialises an empty trap with default parameters
    PenningTrap();
    // Initialises an empty trap with specified parameters
    PenningTrap(double B0, double V0, double d);
    // Initialises a trap with specified parameters and 
    // the particles given in 'parts'
    PenningTrap(double B0, double V0, double d, std::vector<Particle> parts);
    // Initialises a trap with default parameters and 
    // the particles given in 'parts'
    PenningTrap(std::vector<Particle> parts);

    PenningTrap(bool interaction);

    // Initialises a trap with default parameters and 'n' particles with 
    // random positions and velocities
    PenningTrap(int n, bool interaction);


    // Add a particle 'p' to the trap
    void add_particle(Particle p);

    // Counts number of particles within the electric and 
    // magnetic field of the trap
    const int num_parts_in_trap();

    // Run the system in 'trap' for time 't'
    int run_experiment(int n, double t, std::string filename,
                        std::string method = "RK4");

    // Writes the exact solution for a single particle 'p' to file
    const int exact_sol(int n, double t, Particle p, 
                        std::string filename);

private:
    // Coulomb constant
    const double k = 1.3893533e5;

    // External magnetic field strength
    double B0_ = 96.5;
    // Electric potential field strength
    double V0_ = 2.41e6;
    // Characteristic dimension of trap
    double d_ = 500.;

    // Switch Coulomb interaction between particles 
    // on ('true') or off ('false')
    bool interaction_ = true;


    // Returns the external electric and magnetic field 
    // at position 'r'
    const arma::vec3 E_field(arma::vec3 r);
    const arma::vec3 B_field(arma::vec3 r);

    // Returns frequencies for the exact solution of a single particle in the trap
    const double omega_0(Particle p);
    const double omega_z(Particle p);
    // The frequencies omega_+ and omega_-
    const double omega_1(Particle p);
    const double omega_2(Particle p);


    // Force on particle 'i' from particle 'j'
    const arma::vec3 interaction_force(int i, int j);
    // The total force on particle 'i' from the external fields
    const arma::vec3 external_force(int i);
    // The total force on particle 'i' from the other particles
    const arma::vec3 total_interaction_force(int i);
    // The total force on particle 'i' from both external fields 
    // and other particles
    const arma::vec3 force(int i);


    // Evolve the system one time step (dt) using Euler method
    void evolve_Euler(double dt);
    // Evolve the system one time step (dt) using 
    // Runge-Kutta 4th order
    void evolve_RK4(double dt);

};


#endif