#include "config.h"
#include "sim_ini.h"
#include "voronoi_particle.h"
#include <cmath>
#include "sim.h"
#include <vector>
#include "r2.h"
#include "my_random.h"

// color the fluid with two different colors
void sim_ini::two_colors(double c_top, double border, double c_bottom, Sim* sim)
{
    VoronoiParticle* particle;
    double y = 0;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];
        y = particle->m_location.y;

        if (y > border)
        {
            particle->m_color = c_top;
        }
        else
        {
            particle->m_color = c_bottom;
        }
    }
}

// perturbation function for Kelvin-Helmholtz
void sim_ini::perturbing_v_y_arepo(double y1, double y2, double w_0, double sigma, Sim* sim)
{
    VoronoiParticle* particle = NULL;
    double x = 0;
    double y = 0;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];
        x = particle->m_location.x;
        y = particle->m_location.y;

        particle->m_velocity.y = w_0 * sin(4 * M_PI * x) * (exp(-pow(y - y1, 2) / (2 * sigma * sigma)) +
                                                            exp(-pow(y - y2, 2) / (2 * sigma * sigma)));
    }
}

// initialize the primitive variables for the Rayleigh-Taylor instability (Arepo)
void sim_ini::primVarIniRT_arepo(Sim* sim)
{
    // box: x: [0,0.5]	y: [0,1.5]
    double g = -0.1;
    double P_0 = 2.5;
    double w_0 = 0.0025;

    VoronoiParticle* particle = NULL;
    double x = 0;
    double y = 0;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];
        x = particle->m_location.x;
        y = particle->m_location.y;

        if (y > 0.75)
        {
            particle->m_density = 2;
        }
        else
        {
            particle->m_density = 1;
        }

        particle->m_pressure = P_0 + g * (y - 0.75) * particle->m_density;

        particle->m_velocity.x = 0;
        particle->m_velocity.y = w_0 * (1 - cos(4. * M_PI * x)) * (1 - cos(4 * M_PI * y / 3.));
    }
}

// initialize the primitive variables for the Rayleigh-Taylor instability (Arepo, square box)
void sim_ini::primVarIniRT_arepo_square(Sim* sim)
{
    // box: x: [0,1]	y: [0,1]
    double g = -0.1;
    double P_0 = 2.5;
    double w_0 = 0.0015;

    VoronoiParticle* particle = NULL;
    double x = 0;
    double y = 0;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];
        x = particle->m_location.x;
        y = particle->m_location.y;

        if (y > 0.5)
        {
            particle->m_density = 2;
        }
        else
        {
            particle->m_density = 1;
        }

        particle->m_pressure = P_0 + g * (y - 0.5) * particle->m_density;

        particle->m_velocity.x = 0;
        particle->m_velocity.y = w_0 * (1 - cos(2. * M_PI * x)) * (1 - cos(2 * M_PI * y));
    }
}

// initialize the primitive variables for the Rayleigh-Taylor instability
void sim_ini::primVarIniRT(Sim* sim, double border)
{
    double g = Constants::getACCELERATION_Y();
    double P_0 = 2.5;

    VoronoiParticle* particle = NULL;

    double y = 0;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];

        y = particle->m_location.y;

        if (y > border)
        {
            particle->m_density = 2;
        }
        else
        {
            particle->m_density = 1;
        }

        particle->m_pressure = P_0 + g * (y - border) * particle->m_density;

        particle->m_velocity.x = myRandom::randomDouble(0, 0.1);
        particle->m_velocity.y = myRandom::randomDouble(0, 0.1);
    }
}

// the pressure is set according to the hydrostatic equilibrium
void sim_ini::set_pressure_hydrostatic_equilibrium(Sim* sim, double P_0, double height_P_0)
{
    // double middle = (Constants::getBOX_YMAX() + Constants::getBOX_YMIN()) * 0.5;
    double middle = height_P_0;

    double g = Constants::getACCELERATION_Y();
    double y = 0;
    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];

        y = particle->m_location.y;

        particle->m_pressure = P_0 + g * (y - middle) * particle->m_density;
    }
}

// initialize the primitive variables for a quad shock simulation
void sim_ini::prim_var_ini_quad_shock(Sim* sim)
{
    VoronoiParticle* vp = NULL;
    Point location = Point(0, 0);

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        vp = &sim->m_particles[i];
        location = vp->m_location;

        if (location.x >= 0 && location.y >= 0) // PP
        {
            vp->m_density = 1;
            vp->m_velocity = Vector(0.75, -0.5);
            vp->m_pressure = 1;
        }
        else if (location.x < 0 && location.y > 0) // MP
        {
            vp->m_density = 2;
            vp->m_velocity = Vector(0.75, 0.5);
            vp->m_pressure = 1;
        }
        else if (location.x <= 0 && location.y <= 0) // MM
        {
            vp->m_density = 1;
            vp->m_velocity = Vector(-0.75, 0.5);
            vp->m_pressure = 1;
        }
        else if (location.x > 0 && location.y < 0) // PM
        {
            vp->m_density = 3;
            vp->m_velocity = Vector(-0.75, -0.5);
            vp->m_pressure = 1;
        }
        else
        {
            assert(false);
        }
    }
}

// add a boost velocity in x-direction
void sim_ini::x_boost(Sim* sim, double v_x)
{
    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].m_velocity.x += v_x;
    }
}

// add a boost velocity in y-direction
void sim_ini::y_boost(Sim* sim, double v_y)
{
    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].m_velocity.y += v_y;
    }
}