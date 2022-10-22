#include "PenningTrap.hpp"

#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>

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
    arma::vec3 F = k * p1.charge() * p2.charge() * x / std::pow(arma::norm(x), 3);

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
const arma::vec3 PenningTrap::total_force(int i, bool interaction)
{
    arma::vec3 F_tot = total_force_external(i);

    if (interaction)
    {
        F_tot += total_force_particles(i);
    }

    return F_tot;
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_Euler(double dt, bool interaction)
{
    int n = parts.size();

    std::vector<Particle> old = parts;
    arma::vec3 newpos, newvel;

    std::vector<arma::vec3> forces(n);
    for (int i = 0; i < n; i++)
    {
        forces[i] = total_force(i, interaction);
    }

    for (int i = 0; i < n; i++)
    {
        newpos = old[i].pos() + dt * old[i].vel();
        newvel = old[i].vel() + dt / old[i].mass() * forces[i];

        parts[i].set_pos(newpos);
        parts[i].set_vel(newvel);
    }
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, bool interaction)
{
    int n = parts.size();
    std::vector<Particle> copy = parts;
    std::vector<arma::vec3> kx1(n), kx2(n), kx3(n), kx4(n);
    std::vector<arma::vec3> kv1(n), kv2(n), kv3(n), kv4(n);

    // Find k1
    for (int i = 0; i < n; i++)
    {
        kx1[i] = dt * parts[i].vel();
        kv1[i] = dt / parts[i].mass() * total_force(i, interaction);
    }

    // Set first intermediate position and velocity of all particles
    for (int i = 0; i < n; i++)
    {
        parts[i].set_pos(copy[i].pos() + kx1[i]/2.);
        parts[i].set_vel(copy[i].vel() + kv1[i]/2.);
    }

    // Find k2
    for (int i = 0; i < n; i++)
    {
        kx2[i] = dt * parts[i].vel();
        kv2[i] = dt / parts[i].mass() * total_force(i, interaction);
    }

    // Set second intermediate position and velocity of all particles
    for (int i = 0; i < n; i++)
    {
        parts[i].set_pos(copy[i].pos() + kx2[i]/2.);
        parts[i].set_vel(copy[i].vel() + kv2[i]/2.);
    }

    // Find k3
    for (int i = 0; i < n; i++)
    {
        kx3[i] = dt * parts[i].vel();
        kv3[i] = dt / parts[i].mass() * total_force(i, interaction);
    }

    // Set third intermediate position and velocity of all particles
    for (int i = 0; i < n; i++)
    {
        parts[i].set_pos(copy[i].pos() + kx3[i]);
        parts[i].set_vel(copy[i].vel() + kv3[i]);
    }

    // Find k4
    for (int i = 0; i < n; i++)
    {
        kx4[i] = dt * parts[i].vel();
        kv4[i] = dt / parts[i].mass() * total_force(i, interaction);
    }

    // Set final position and velocity
    for (int i = 0; i < n; i++)
    {
        parts[i].set_pos(copy[i].pos() 
                            + (kx1[i] + 2.*(kx2[i]+kx3[i]) + kx4[i]) / 6.);
        parts[i].set_vel(copy[i].vel() 
                            + (kv1[i] + 2.*(kv2[i]+kv3[i]) + kv4[i]) / 6.);
    }
}

// Run the system in 'trap' for time 't'
int PenningTrap::runExperiment(int n, double t, std::string filename,
                                bool interaction, std::string method)
{
    double dt = t/n;
    int n_parts = parts.size();

    std::vector<std::ofstream> ofile(n_parts);

    for (int i = 0; i < n_parts; i++)
    {
        ofile[i].open(filename + std::to_string(i+1) + ".txt");
    }

    int prec = 10;
    int width = prec + 10;

    for (int i = 0; i < n_parts; i++)
    {
        ofile[i] << std::setw(width) << std::setprecision(prec)
                    << std::scientific << 0.;
        for (int j = 0; j < 3; j++)
        {
            ofile[i] << std::setw(width) << std::setprecision(prec)
                        << std::scientific << parts[i].pos()(j);
            ofile[i] << std::setw(width) << std::setprecision(prec)
                        << std::scientific << parts[i].vel()(j);
        }
    }
    for (int i = 0; i < n_parts; i++)
    {
        ofile[i] << std::endl;
    }

    for (int m = 0; m < n; m++)
    {
        if (method == "RK4")
        {
            evolve_RK4(dt, interaction);
        }
        else if (method == "Euler")
        {
            evolve_Euler(dt, interaction);
        }
        else
        {
            std::cout << "Error: No method called " << method << std::endl;
            return 1;
        }

        for (int i = 0; i < n_parts; i++)
        {
            ofile[i] << std::setw(width) << std::setprecision(prec)
                        << std::scientific << (m+1)*dt;
            for (int j = 0; j < 3; j++)
            {
                ofile[i] << std::setw(width) << std::setprecision(prec)
                            << std::scientific << parts[i].pos()(j);
                ofile[i] << std::setw(width) << std::setprecision(prec)
                            << std::scientific << parts[i].vel()(j);
            }
            ofile[i] << std::endl;
        }
    }

    for (int i = 0; i < n_parts; i++)
    {
        ofile[i].close();
    }

    return 0;
}


double PenningTrap::omega_0(Particle p)
{
    return p.charge() * B0_ / p.mass();
}


double PenningTrap::omega_z(Particle p)
{
    return std::sqrt(2.*p.charge() * Vd_ / p.mass());
}


double PenningTrap::omega_1(Particle p)
{
    return (omega_0(p) + std::sqrt(std::pow(omega_0(p),2)
            - 2.*std::pow(omega_z(p),2))) / 2.;
}


double PenningTrap::omega_2(Particle p)
{
    return (omega_0(p) - std::sqrt(std::pow(omega_0(p),2)
            - 2.*std::pow(omega_z(p),2))) / 2.;
}


int PenningTrap::exact_sol(int n, double t_tot, Particle p,
                            std::string filename)
{
    double dt = t_tot/n;

    double wz = omega_z(p);
    double w1 = omega_1(p);
    double w2 = omega_2(p);

    double x0 = p.pos()(0);
    double v0 = p.vel()(1);
    double z0 = p.pos()(2);

    double A1 = (v0 + w2 * x0)/(w1 - w2);
    double A2 = (v0 + w1 * x0)/(w1 - w2);

    // t, x, vx, y, vy, z, vz
    std::vector<double> exactval(7);
    double t;

    std::ofstream ofile;
    ofile.open(filename);
    int prec = 20;
    int width = prec + 10;

    for (int i = 0; i < n+1; i++)
    {
        t = i * dt;
        exactval[0] = t;
        exactval[1] = - A1 * std::cos(w1 * t) + A2 * std::cos(w2 * t);
        exactval[2] = A1 * w1 * std::sin(w1 * t)
                        - A2 * w2 * std::sin(w2 * t);
        exactval[3] = A1 * std::sin(w1 * t) - A2 * std::sin(w2 * t);
        exactval[4] = A1 * w1 * std::cos(w1 * t)
                        - A2 * w2 * std::cos(w2 * t);
        exactval[5] = z0 * std::cos(wz * t);
        exactval[6] = - z0 * wz * std::sin(wz * t);

        for (int j = 0; j < 7; j++)
        {
            ofile << std::setw(width) << std::setprecision(prec)
                << std::scientific << exactval[j];
        }
        ofile << std::endl;
    }
    return 0;
}