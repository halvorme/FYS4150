#include <string>
#include <iostream>

#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "simulations.hpp"

int main()
{
    int n = 4000;
    double t = 50;

    arma::vec3 ipos1 = {20, 0, 20};
    arma::vec3 ivel1 = {0, 25, 0};

    arma::vec3 ipos2 = {25, 25, 0};
    arma::vec3 ivel2 = {0, 40, 5};

    // arma::vec3 init_vel = {0, 0, 0};

    single_part(n, t, ipos1, ivel1);

    std::cout << "two" << std::endl;
    two_parts(n, t, ipos1, ivel1, ipos2, ivel2, true);

    return 0;
}
