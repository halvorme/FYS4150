#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle
{
private:
    double q_, m_;
    double qm_;
    arma::vec3 r_, v_;

    double k = 1.3893533e5;

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

    // Coloumb force from particle 'p2'
    // const arma::vec3 F_Coloumb(Particle p2);

};

#endif