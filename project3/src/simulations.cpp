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

    trap.run_experiment(n, t, "data/singlepart");

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

    trap_int.run_experiment(n, t, "data/twoparts_int");
    trap_noint.run_experiment(n, t, "data/twoparts_noint", false);

    return 0;
}


int error_sim()
{
    double t = 50.;

    std::vector<int> n(4);
    n[0] = 4000;
    n[1] = 8000;
    n[2] = 16000;
    n[3] = 32000;

    arma::vec3 r0 = {20 ,0, 20};
    arma::vec3 v0 = {0, 25, 0};

    Particle p0(r0, v0);

    std::vector<PenningTrap> traps(8);

    for (int i = 0; i < 4; i++)
    {
        traps[i].add_particle(p0);

        traps[i].run_experiment(n[i], t, "data/errorRK4_" 
                                + std::to_string(n[i]) + "_");

        traps[i+4].add_particle(p0);
        traps[i+4].run_experiment(n[i], t, "data/errorEuler_" + std::to_string(n[i]) 
                                     + "_", true, "Euler");
    }

    return 0;
}


int single_exact()
{
    double t = 50.;

    std::vector<int> n(4);
    n[0] = 4000;
    n[1] = 8000;
    n[2] = 16000;
    n[3] = 32000;

    arma::vec3 r0 = {20 ,0, 20};
    arma::vec3 v0 = {0, 25, 0};

    PenningTrap init_trap;
    Particle p(r0, v0);

    init_trap.add_particle(p);

    std::string filename = "data/exact_";

    for (int i = 0; i < 4; i++)
    {
        init_trap.exact_sol(n[i], t, init_trap.parts[0], 
                            filename + std::to_string(n[i]) + ".txt");
    }

    return 0;
}