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

    std::cout << trap.num_parts_in_trap() << std::endl;

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

    PenningTrap trap_int, trap_noint(false);

    Particle p1(ipos1, ivel1), p2(ipos2, ivel2);

    trap_int.add_particle(p1);
    trap_int.add_particle(p2);

    trap_noint.add_particle(p1);
    trap_noint.add_particle(p2);

    trap_int.run_experiment(n, t, "data/twoparts_int");
    trap_noint.run_experiment(n, t, "data/twoparts_noint");

    std::cout << trap_int.num_parts_in_trap() << std::endl;
    std::cout << trap_noint.num_parts_in_trap() << std::endl;

    return 0;
}


int error_sim()
{
    double t = 50.;

    std::vector<int> n{4000, 8000, 16000, 32000};

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
        traps[i+4].run_experiment(n[i], t, "data/errorEuler_" 
                                    + std::to_string(n[i]) + "_", "Euler");
    }

    return 0;
}


int single_exact()
{
    double t = 50.;

    std::vector<int> n{4000, 8000, 16000, 32000};
    // n[0] = 2000;
    // n[1] = 8000;
    // n[2] = 16000;
    // n[3] = 32000;

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


int trapped_broad()
{
    double t = 500.;

    int n = 4000;

    std::vector<double> f{.1, .4, .7};
    double d_omega = .02;
    bool interaction = false;
    int n_parts = 100;

    double omega_V = 0;

    std::ofstream ofile;
    ofile.open("data/resonance_broad.txt");

    int prec = 14;
    int width = prec + 10;

    while (omega_V <= 2.5)
    {
        std::vector<int> remaining(3);
        omega_V += d_omega;

        for (int i = 0; i < 3; i++)
        {
            arma::arma_rng::set_seed(1);
            PenningTrap trap(n_parts, f[i], omega_V, interaction);

            trap.run_experiment_noprint_test(n, t);

            remaining[i] = trap.num_parts_in_trap();
        }

        ofile << std::setw(width) << std::setprecision(prec)
                << std::scientific << omega_V;
        for (int i = 0; i < 3; i++)
        {
            ofile << std::setw(width) << std::setprecision(prec)
                << std::scientific << remaining[i];
        }
        ofile << std::endl;

        std::cout << "Yo man" << std::endl;
    }

    ofile.close();

    return 0;

}