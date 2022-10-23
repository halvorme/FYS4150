#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>


// Particle, described by its charge, mass, position and velocity.
// Defaults to a ion with mass 40 u and charge +1 e.
class Particle
{
public:
    // Constructors
    Particle();
    Particle(arma::vec3 r, arma::vec3 v);
    Particle(double q, double m, arma::vec3 r, arma::vec3 v);

    // Read particle values
    const double charge();
    const double mass();
    const arma::vec3 pos();
    const arma::vec3 vel();

    // Change particle position and velocity
    void set_pos(arma::vec3 r);
    void set_vel(arma::vec3 v);

private:
    // Charge and mass of particle
    double q_ = 1.;
    double m_ = 40.;

    // Position and velocity of particle
    arma::vec3 r_ = arma::vec3(arma::fill::zeros);
    arma::vec3 v_ = arma::vec3(arma::fill::zeros);
};

#endif