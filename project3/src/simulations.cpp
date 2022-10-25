#include "simulations.hpp"

#include <armadillo>
#include <fstream>
#include <iomanip>
#include <string>

#include "PenningTrap.hpp"
#include "Particle.hpp"


int single_part(int n, double t)
{
    arma::vec3 init_pos = {20, 0, 20};
    arma::vec3 init_vel = {0, 25, 0};

    PenningTrap trap;

    Particle p(init_pos, init_vel);

    trap.add_particle(p);

    trap.run_experiment(n, t, "data/singlepart");

    return 0;
}

int two_parts(int n, double t)
{
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

    return 0;
}


int error_sim(std::vector<int> n, double t)
{
    arma::vec3 r0 = {20 ,0, 20};
    arma::vec3 v0 = {0, 25, 0};

    Particle p0(r0, v0);

    std::vector<PenningTrap> traps(2*n.size());

    for (int i = 0; i < n.size(); i++)
    {
        traps[i].add_particle(p0);

        traps[i].run_experiment(n[i], t, "data/errorRK4_" 
                                + std::to_string(n[i]) + "_");

        traps[i+n.size()].add_particle(p0);

        traps[i+n.size()].run_experiment(n[i], t, "data/errorEuler_" 
                                    + std::to_string(n[i]) + "_", "Euler");
    }

    return 0;
}


int single_exact(std::vector<int> n, double t)
{
    arma::vec3 r0 = {20 ,0, 20};
    arma::vec3 v0 = {0, 25, 0};

    PenningTrap init_trap;
    Particle p(r0, v0);

    init_trap.add_particle(p);

    std::string filename = "data/exact_";

    for (int i = 0; i < n.size(); i++)
    {
        init_trap.exact_sol(n[i], t, init_trap.parts[0], 
                            filename + std::to_string(n[i]) + ".txt");
    }

    return 0;
}


int search(int n, double t, int n_parts, double omega_low, double omega_up, 
            double d_omega, std::vector<double> f, bool interaction, 
            std::ofstream& ofile)
{
    double omega_V = omega_low;

    int prec = 4;
    int width = prec + 10;

    while (omega_V <= omega_up)
    {
        std::vector<double> frac_remaining(3);
        omega_V += d_omega;

        for (int i = 0; i < f.size(); i++)
        {
            arma::arma_rng::set_seed(1);
            PenningTrap trap(n_parts, f[i], omega_V, interaction);

            trap.run_experiment_noprint_test(n, t);

            frac_remaining[i] = trap.num_parts_in_trap()/double(n_parts);
        }

        ofile << std::setw(width) << std::setprecision(prec)
                << std::scientific << omega_V;
        for (int i = 0; i < f.size(); i++)
        {
            ofile << std::setw(width) << frac_remaining[i];
        }
        ofile << std::endl;

        std::cout << omega_V << std::endl;
    }

    return 0;
}


int trapped_broad_new(int n, double t)
{
    std::vector<double> f{.1, .4, .7};
    double d_omega = .02;
    bool interaction = false;
    int n_parts = 100;

    double omega_low = 0.;
    double omega_up = 2.5;

    std::ofstream ofile;
    ofile.open("data/resonance_broad_" + std::to_string(int(t)) + "_" 
                + std::to_string(n) +".txt");

    search(n, t, n_parts, omega_low, omega_up, d_omega, f, interaction, ofile);

    ofile.close();

    return 0;
}


int resonance(int n, double t, double omega_low, double omega_up, bool interction)
{
    return 0;
}