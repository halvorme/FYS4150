#include "PenningTrap.hpp"

#include <vector>

#include "Particle.hpp"

// Constructors
PenningTrap::PenningTrap()
{
    B0_ = 96.5;
    V0_ = 2.41e6;
    d_ = 500.;
    Vd_ = V0_/(d_*d_);
}

PenningTrap::PenningTrap(double B0, double V0, double d)
{
    B0_ = B0;
    V0_ = V0;
    d_ = d;
    Vd_ = V0/(d*d);
}

PenningTrap::PenningTrap(double B0, double V0, double d, 
                            std::vector<Particle> particles)
{
    B0_ = B0;
    V0_ = V0;
    d_ = d;
    Vd_ = V0/(d*d);
    parts = particles;
}

PenningTrap::PenningTrap(std::vector<Particle> particles)
{
    B0_ = 96.5;
    V0_ = 2.41e6;
    d_ = 500.;
    Vd_ = V0_/(d_*d_);
    parts = particles;
}

// Electric and magnetic field
const arma::vec3 PenningTrap::E_field(arma::vec3 r)
{
    arma::vec3 E;
    E(0) = Vd_*r(0);
    E(1) = Vd_*r(1);
    E(2) = -2.*Vd_*r(2);
    return E;
}

const arma::vec3 PenningTrap::B_field(arma::vec3 r)
{
    arma::vec3 B(arma::fill::zeros);
    B(2) = B0_;
    return B;
}

void PenningTrap::add_particle(Particle p_in)
{
    parts.push_back(p_in);
}


const arma::vec3 PenningTrap::force_particle(int i, int j)
{
    Particle p1 = parts[i];
    Particle p2 = parts[j];

    arma::vec3 x = p1.pos() - p2.pos();
    arma::vec3 F = k * p1.charge() * p2.charge() * x / pow(arma::norm(x), 3);

    return F;
}

// The total force on particle_i from the external fields
const arma::vec3 PenningTrap::total_force_external(int i)
{
    Particle p = parts[i];

    arma::vec3 F_ext = p.charge() * (E_field(p.pos()) 
                            + arma::cross(p.vel(), B_field(p.pos())));
    
    return F_ext;
}

// The total force on particle_i from the other particles
const arma::vec3 PenningTrap::total_force_particles(int i)
{
    int n = parts.size();

    arma::vec3 F(arma::fill::zeros);

    for (int j  = 0; j < n; j++)
    {
        if (j != i)
        {
            F += force_particle(i,j);
        }
    }

    return F;
}

// The total force on particle_i from both external fields and other particles
const arma::vec3 PenningTrap::total_force(int i)
{
    arma::vec3 F_tot = total_force_external(i) + total_force_particles(i);

    return F_tot;
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
    int n = parts.size();
    std::vector<Particle> old = parts;
    arma::vec3 newpos, newvel;

    std::vector<arma::vec3> forces(n);
    for (int i = 0; i < n; i++)
    {
        forces[i] = total_force(i);
    }

    for (int i = 0; i < n; i++)
    {
        newpos = old[i].pos() + dt * old[i].vel();
        newvel = old[i].vel() + dt / old[i].mass() * forces[i];

        parts[i].set_pos(newpos);
        parts[i].set_vel(newvel);
    }
}
