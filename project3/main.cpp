#include <string>
#include <iostream>

#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "simulations.hpp"

int main()
{
    // 4000 to generate smooth plots
    // 400 gives rel error ~10^-3 for t = 50
    int n = 6000;
    double t = 500.;


    // single_part(n, t);

    // two_parts(n, t);

    // std::vector<int> n_list{4000, 8000, 16000, 32000};
    // std::vector<int> n_list{2000, 4000, 6000, 8000, 16000};
    // t = 50;

    // error_sim(n_list, t);

    // single_exact(n_list, t);

    n = 4000.;
    t = 50.;

    trapped_broad(n, t);

    return 0;
}
