#include "config.h"
#include "obstacle.h"
#include "sim.h"
#include "sim_graph.h"
#include "sim_core.h"

Obstacle::Obstacle(Sim* sim)
        :
        m_sim(sim),
        m_obstacle_particles(nof_obstacle_particles), // global value, s.a. sim.cpp
        m_velocity(0, 0)
{

}

Obstacle::~Obstacle()
{

}

// binary input
void Obstacle::bin_input(unsigned int nof_particles_in_use)
{
    m_nof_obstacle_particles = nof_particles_in_use;

    unsigned int nof_particles_half = (nof_particles_in_use - 1) / 2;

    for (unsigned int i = 0; i < nof_particles_half; i++)
    {
        // inner
        m_sim->m_particles.at(i).m_obstacle_particle = true;
        m_sim->m_particles.at(i).m_inner_obstacle_particle = true;

        // outer
        m_sim->m_particles.at(i + nof_particles_half).m_obstacle_particle = true;

        // inner
        m_obstacle_particles.current() = &m_sim->m_particles.at(i);
        m_obstacle_particles.counterInc();

        // outer
        m_obstacle_particles.current() = &m_sim->m_particles.at(i + nof_particles_half);
        m_obstacle_particles.counterInc();
    }

    m_sim->m_particles.at(2 * nof_particles_half).m_obstacle_particle = true;
    m_sim->m_particles.at(2 * nof_particles_half).m_inner_obstacle_particle = true;

    m_obstacle_particles.current() = &m_sim->m_particles.at(2 * nof_particles_half);
    m_obstacle_particles.counterInc();
}

// minimum distance to the obstacle
double Obstacle::minimum_distance_square()
{
    // same value as dr
    return (pow(delta_r() * 0.5, 2) + pow(delta_phi() * 0.5, 2)) * 1.05;
}

double Obstacle::delta_r()
{
    return m_sim->m_box->m_dh * Constants::getCIRCLE_DR();
}

double Obstacle::delta_phi()
{
    return (m_sim->m_box->m_dh * Constants::getCIRCLE_DPHI());
}

// position of the obstacle
Point Obstacle::obstacle_center()
{
    return m_obstacle_particles.at(m_nof_obstacle_particles - 1)->m_location;
}

// handle a particle inside the obstacle
// assumption: obstacle is at (0.0)
void Obstacle::handle_inside_obstacle(VoronoiParticle* vp)
{
    Vector normal_vector = vp->m_location - obstacle_center();
    normal_vector.normalize();

    // set new location to the rim of the circle
    vp->m_location = normal_vector * (Constants::getRADIUS() + sqrt(minimum_distance_square()));

    // retransform
    vp->m_location += obstacle_center();
}

// create a disk obstacle
void Obstacle::disk()
{
    // center
    Point center = Point(Constants::getCIRCLE_X(), Constants::getCIRCLE_Y());

    double radius = Constants::getRADIUS();
    double dr = delta_r();

    unsigned int particles_phi_direction = (int) 2 * M_PI * radius / delta_phi();

    // number of obstacle particles
    m_nof_obstacle_particles = particles_phi_direction * 2 + 1;

    double radius_inner = radius - dr * 0.5;
    double radius_outer = radius + dr * 0.5;

    // coordinates of the temporary point
    double x = 0;
    double y = 0;
    double phi = 0;
    double distance = 0;

    // remove particles inside the obstacle

    for (unsigned int k = 0; k < m_sim->m_particles.counter(); k++)
    {
        distance = sqrt((center - m_sim->m_particles[k].m_location).length_square());

        if (distance <= Constants::getRADIUS() + sqrt(minimum_distance_square()))
        {
            m_sim->m_particles.del(&m_sim->m_particles[k]);
            k--;
        }
    }

    // free space for obstacle particles
    for (unsigned int l = 0; l < m_nof_obstacle_particles; l++)
    {
        m_sim->m_particles.current() = m_sim->m_particles.at(l);
        m_sim->m_particles.at(l) = VoronoiParticle();
        m_sim->m_particles.counterInc();
    }

    // insert obstacle particles
    for (unsigned int i = 0; i < particles_phi_direction; i++)
    {

        phi = i / (double) particles_phi_direction * (2 * M_PI);

        // inner circle
        x = radius_inner * cos(phi);
        y = radius_inner * sin(phi);

        x += center.x;
        y += center.y;

        m_sim->m_particles.at(i).m_location.x = x;
        m_sim->m_particles.at(i).m_location.y = y;
        m_sim->m_particles.at(i).m_obstacle_particle = true;
        m_sim->m_particles.at(i).m_inner_obstacle_particle = true;

        m_obstacle_particles.current() = &m_sim->m_particles.at(i);

        m_obstacle_particles.counterInc();

        // outer circle
        x = radius_outer * cos(phi);
        y = radius_outer * sin(phi);

        x += center.x;
        y += center.y;

        m_sim->m_particles.at(i + particles_phi_direction).m_location.x = x;
        m_sim->m_particles.at(i + particles_phi_direction).m_location.y = y;
        m_sim->m_particles.at(i + particles_phi_direction).m_obstacle_particle = true;

        m_obstacle_particles.current() = &m_sim->m_particles.at(i + particles_phi_direction);

        m_obstacle_particles.counterInc();
    }

    // add one obstacle particle into the middle
    m_sim->m_particles.at(2 * particles_phi_direction).m_location.x = center.x;
    m_sim->m_particles.at(2 * particles_phi_direction).m_location.y = center.y;
    m_sim->m_particles.at(2 * particles_phi_direction).m_obstacle_particle = true;
    m_sim->m_particles.at(2 * particles_phi_direction).m_inner_obstacle_particle = true;

    m_obstacle_particles.current() = &m_sim->m_particles.at(2 * particles_phi_direction);

    m_obstacle_particles.counterInc();
}

// is a particle inside the obstacle
// assumption: obstacle is at (0,0)
bool Obstacle::is_inside(Point* location)
{
    return (sqrt((*location - obstacle_center()).length_square()) <
            Constants::getRADIUS() + sqrt(minimum_distance_square()));
}

// assign velocities to the obstacle particles
void Obstacle::assign_velocities()
{
    Vector surface_force = Vector(0, 0);

    if (Constants::getOBS_FLOATING())
    {
        surface_force = calculate_force();
    }

    // acceleration: external plus acceleration due to the fluid
    m_velocity.x = m_velocity.x +
                   sim_core::simulation_time * (Constants::getOBS_AX() + surface_force.x / Constants::getOBS_MASS());
    m_velocity.y = m_velocity.y +
                   sim_core::simulation_time * (Constants::getOBS_AY() + surface_force.y / Constants::getOBS_MASS());
    printf("surface force: %f	%f\n", surface_force.x, surface_force.y);

    for (unsigned int i = 0; i < m_nof_obstacle_particles; i++)
    {
        m_obstacle_particles[i]->m_particleVelocity.x = m_velocity.x;
        m_obstacle_particles[i]->m_particleVelocity.y = m_velocity.y;
    }
}

void Obstacle::delete_zero_length_edges()
{
    for (unsigned int i = 0; i < m_nof_obstacle_particles - 2; i += 2)
    {
        HalfEdge* start = m_obstacle_particles[i]->m_outerComponent;
        HalfEdge* currentEdge = start->m_next;
        HalfEdge* next = NULL;

        do
        {
            next = currentEdge->m_next;

            currentEdge->updateLengthandCenter();

            if (currentEdge->is_zero())
            {
                assert(currentEdge != start);
                sim_graph::delete_edge(currentEdge);
            }

            currentEdge = next;

        } while (currentEdge != start);

        // start half edge
        currentEdge->updateLengthandCenter();

        if (currentEdge->is_zero())
        {
            sim_graph::delete_edge(currentEdge);
        }
    }
}

// copy primitive variables and cell centers to the inner obstacle cells
void Obstacle::copyPrimitiveVariablesCellCentersAndVolume()
{
    Vector n = Vector(0, 0);
    Vector s = Vector(0, 0);

    // inner obstacle particles only
    for (unsigned int i = 0; i < m_nof_obstacle_particles - 2; i += 2)
    {
        n = m_obstacle_particles[i + 1]->m_location -
            m_obstacle_particles[i]->m_location; // outer particle - inner particle
        n.normalize();
        s = n;
        s.rotate(M_PI * 0.5);

        // primitive variables
        m_obstacle_particles[i]->m_density = m_obstacle_particles[i + 1]->m_density;
        m_obstacle_particles[i]->m_pressure = m_obstacle_particles[i + 1]->m_pressure;

        m_obstacle_particles[i]->m_velocity = n * (-1.) * m_obstacle_particles[i + 1]->m_velocity.scalarProduct(n) +
                                              s * m_obstacle_particles[i + 1]->m_velocity.scalarProduct(s);
        // relevant for moving obstacles
        m_obstacle_particles[i]->m_velocity += m_obstacle_particles[i]->m_particleVelocity * 2;

        // cell centers
        Point particles_center = (m_obstacle_particles[i]->m_location + m_obstacle_particles[i + 1]->m_location) * 0.5;
        Vector relative = m_obstacle_particles[i + 1]->m_cellCenterOfVolume - particles_center;
        Vector relative_mirrored = n * (-1.) * relative.scalarProduct(n) + s * relative.scalarProduct(s);
        m_obstacle_particles[i]->m_cellCenterOfVolume = particles_center + relative_mirrored;

        // volume
        m_obstacle_particles[i]->m_cellVolume = m_obstacle_particles[i + 1]->m_cellVolume;
    }
}

// copy gradients to the inner obstacle cells
void Obstacle::copyGradients()
{
    Vector n = Vector(0, 0);
    Vector s = Vector(0, 0);
    Vector grad_vx_mirrored = Vector(0, 0);
    Vector grad_vy_mirrored = Vector(0, 0);

    for (unsigned int i = 0; i < m_nof_obstacle_particles - 2; i += 2)
    {
        n = m_obstacle_particles[i + 1]->m_location -
            m_obstacle_particles[i]->m_location; // outer particle - inner particle
        n.normalize();
        s = n;
        s.rotate(M_PI * 0.5);

        // density and pressure gradient
        m_obstacle_particles[i]->m_gradientDensity =
                n * (-1.) * m_obstacle_particles[i + 1]->m_gradientDensity.scalarProduct(n) +
                s * m_obstacle_particles[i + 1]->m_gradientDensity.scalarProduct(s);
        m_obstacle_particles[i]->m_gradientPressure =
                n * (-1.) * m_obstacle_particles[i + 1]->m_gradientPressure.scalarProduct(n) +
                s * m_obstacle_particles[i + 1]->m_gradientPressure.scalarProduct(s);

        // velocity_x and velocity_y gradient
        grad_vx_mirrored = n * (-1.) * m_obstacle_particles[i + 1]->m_gradientVelocityX.scalarProduct(n) +
                           s * m_obstacle_particles[i + 1]->m_gradientVelocityX.scalarProduct(s);
        grad_vy_mirrored = n * (-1.) * m_obstacle_particles[i + 1]->m_gradientVelocityY.scalarProduct(n) +
                           s * m_obstacle_particles[i + 1]->m_gradientVelocityY.scalarProduct(s);

        m_obstacle_particles[i]->m_gradientVelocityX =
                grad_vx_mirrored * (-n.x * n.x + s.x * s.x) + grad_vy_mirrored * (-n.x * n.y + s.x * s.y);
        m_obstacle_particles[i]->m_gradientVelocityY =
                grad_vx_mirrored * (-n.x * n.y + s.x * s.y) + grad_vy_mirrored * (-n.y * n.y + s.y * s.y);

    }
}

// calculate the force of the fluid
Vector Obstacle::calculate_force()
{
    HalfEdge* half_edge = NULL;
    VoronoiParticle* outer = NULL;
    Vector force = Vector(0, 0);
    Vector direction = Vector(0, 0); // direction of the force on an edge

    // loop over inner obstacle particles
    for (unsigned int i = 0; i < m_nof_obstacle_particles - 2; i += 2)
    {
        half_edge = m_obstacle_particles[i]->m_outerComponent;

        while (!(half_edge->m_twin->m_vp == m_obstacle_particles[i + 1]))
        {
            half_edge = half_edge->m_next;
        }

        outer = m_obstacle_particles[i + 1];

        direction = m_obstacle_particles[i]->m_location - outer->m_location;
        direction.normalize();

        force += direction * half_edge->m_length * outer->m_pressure;
    }

    return force;
}