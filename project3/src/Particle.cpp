#include "Particle.hpp"

#include <armadillo>

// Constructors
Particle::Particle() {}

Particle::Particle(arma::vec3 r, arma::vec3 v)
    : r_(r), v_(v)
{}

Particle::Particle(double q, double m, arma::vec3 r, arma::vec3 v)
    : q_(q), m_(m), r_(r), v_(v)
{}

// Read particle charge
const double Particle::charge()
{
    return q_;
}

// Read particle mass
const double Particle::mass()
{
    return m_;
}

// Read particle position
const arma::vec3 Particle::pos()
{
    return r_;
}

// Read particle velocity
const arma::vec3 Particle::vel()
{
    return v_;
}

// Change position
void Particle::set_pos(arma::vec3 r)
{
    r_ = r;
}

// Change velocity
void Particle::set_vel(arma::vec3 v)
{
    v_ = v;
}
