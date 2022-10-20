#include <string>
#include <iostream>

#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "simulations.hpp"

int main()
{
    int n = 4000;
    double t = 20;

    arma::vec3 init_pos = {20, 0, 20};
    arma::vec3 init_vel = {0, 25, 0};

    single_part(n, t, init_pos, init_vel);

    return 0;
}
