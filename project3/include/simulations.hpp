#ifndef __simulations_hpp__
#define __simulations_hpp__

#include <armadillo>

#include "PenningTrap.hpp"
#include "Particle.hpp"

int single_part(int n, double t, arma::vec3 init_pos, arma::vec3 init_vel);

int two_parts(int n, double t, arma::vec3 ipos1, arma::vec3 ivel1, 
                arma::vec3 ipos2, arma::vec3 ivel2, bool interaction);

#endif