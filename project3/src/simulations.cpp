#include "simulations.hpp"

#include <armadillo>
#include <fstream>
#include <iomanip>
#include <string>

#include "PenningTrap.hpp"
#include "Particle.hpp"


int single_part()
{
    int n = 4000;
    double t = 50;

    arma::vec3 init_pos = {20, 0, 20};
    arma::vec3 init_vel = {0, 25, 0};

    PenningTrap trap;

    Particle p(init_pos, init_vel);

    trap.add_particle(p);

    trap.runExperiment(n, t, "data/singlepart");

    return 0;
}

int two_parts()
{
    int n = 4000;
    double t = 50;

    arma::vec3 ipos1 = {20, 0, 20};
    arma::vec3 ivel1 = {0, 25, 0};

    arma::vec3 ipos2 = {25, 25, 0};
    arma::vec3 ivel2 = {0, 40, 5};

    PenningTrap trap_int, trap_noint;

    Particle p1(ipos1, ivel1), p2(ipos2, ivel2);

    trap_int.add_particle(p1);
    trap_int.add_particle(p2);

    trap_noint.add_particle(p1);
    trap_noint.add_particle(p2);

    trap_int.runExperiment(n, t, "data/twoparts_int");
    trap_noint.runExperiment(n, t, "data/twoparts_noint", false);

    return 0;
}
