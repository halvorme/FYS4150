#include "Particle.hpp"

#include <armadillo>

// Constructors
Particle::Particle()
{
    q_ = 1.;
    m_ = 40.;
    qm_ = q_/m_;

    r_ = arma::vec3(arma::fill::zeros);
    v_ = arma::vec3(arma::fill::zeros);
}

Particle::Particle(arma::vec3 r, arma::vec3 v)
{
    q_ = 1.;
    m_ = 40.;
    qm_ = q_/m_;

    r_ = r;
    v_ = v;
}

Particle::Particle(double q, double m, arma::vec3 r, arma::vec3 v)
{
    q_ = q;
    m_ = m;
    qm_ = q/m;

    r_ = r;
    v_ = v;
}

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
