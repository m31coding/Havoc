#include "config.h"
#include "sim_core.h"
#include "debug.h"
#include "sweepline.h"
#include "my_random.h"
#include "binary_input_output.h"
#include "sim_graph.h"
#include <limits>
#include "point_location.h"

// global time step
double sim_core::global_timestep = 0; // calculated every time step

// simulation time
double sim_core::simulation_time = 0; // calculated every time step

// time for the next binary output
double sim_core::next_output_time = 0;

// is it time for a binary output?
bool sim_core::binary_output = false;

// output counter
unsigned int sim_core::counter = 0; // output file, e.g.: 5 => 00005.hc

// calculate the fluxes
void sim_core::calculateFluxes(Sim* sim)
{
    const unsigned int counter = sim->m_halfEdges.counter();

    CO_PRINTF("\n---calculating the fluxes across the edges---\n\n");

    for (unsigned int i = 0; i < counter; i++)
    {
        if (sim->m_halfEdges[i].m_inUse)
        {
            sim->m_halfEdges[i].updateFlux();
        }
    }
}

// relax the initial data, the particles move to their centers, hexagons are formed
void sim_core::relax(Sim* sim, unsigned int steps)
{
    VoronoiParticle* particle = NULL;

    for (unsigned int k = 0; k < steps; k++)
    {
        printf("relax %u / %u\n", k, steps);

        sim_core::updateCellVolumesAndCenters(sim);

        // move particles
        for (unsigned int j = 0; j < sim->m_particles.counter(); j++)
        {
            particle = &sim->m_particles[j];
            particle->m_location += (particle->m_cellCenterOfVolume - particle->m_location) * 0.5;
            assert(sim->m_box->isInside(&particle->m_location));
        }

        sim->m_box->determineGhostParticles();
        sim->freeGraph();
        sim->m_box->createGhostParticles();
        sim->m_box->reset();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
    }
}

// shift the outer components
void sim_core::shiftOuterComponents(Sim* sim)
{
    VoronoiParticle* particle = NULL;
    int rdm_int;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];
        rdm_int = myRandom::randomInteger(0, 5);

        for (int k = 0; k < rdm_int; k++)
        {
            particle->m_outerComponent = particle->m_outerComponent->m_next;
        }
    }
}

// helper function to adapt the time step for a binary output
static void adapt_timestep()
{
    if (sim_core::simulation_time + sim_core::global_timestep >= sim_core::next_output_time)
    {
        sim_core::global_timestep = sim_core::next_output_time - sim_core::simulation_time;

        assert(sim_core::global_timestep > 0);

        sim_core::binary_output = true;
    }
}

// determine the global time step
void sim_core::constant_global_timestep(Sim* sim)
{
    CO_PRINTF("\n---determining global time step (constant) ---\n\n");

    global_timestep = Constants::getDELTA_T();

    if (Constants::getEVERY_NTH_STEPS_OUTPUT() == -1)
    {
        adapt_timestep();
    }

    CO_PRINTF("global time step: %f\n", global_timestep);
}

// determine the global time step zero_grad
void sim_core::variable_global_timestep(Sim* sim)
{
    CO_PRINTF("\n---determining global time step (zero_grad) ---\n\n");

    double minimum_timestep = std::numeric_limits<double>::max();
    double localTimestep = 0;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        localTimestep = sim->m_particles[i].local_timestep();

        if (localTimestep < minimum_timestep)
        {
            minimum_timestep = localTimestep;
        }
    }

    global_timestep = minimum_timestep;

    adapt_timestep();

    CO_PRINTF("global time step: %f\n", global_timestep);
}

// check whether it's time for a binary output
void sim_core::maybe_binary_output(Sim* sim, int step)
{
    if (Constants::getEVERY_NTH_STEPS_OUTPUT() != -1)
    {
        if (step % Constants::getEVERY_NTH_STEPS_OUTPUT() == 0)
        {
            binary_output = true;
        }
    }

    if (binary_output)
    {
        bin_io::binaryOutput(NULL, counter, sim);
        counter++;
        binary_output = false;
        next_output_time += Constants::getDELTA_T_OUTPUT();
    }
}

// calculate and set the volumes and the centers of volume of all cells
void sim_core::updateCellVolumesAndCenters(Sim* sim)
{
    CO_PRINTF("\n---updating cell volumes and centers of volume---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].updateCellVolumeAndCenter();
    }
}

// calculate and set the length and the center of every half edge
void sim_core::updateHalfEdgeLengthsAndCenters(Sim* sim)
{
    CO_PRINTF("\n---updating lengths and centers of the half edges---\n\n");

    for (unsigned int i = 0; i < sim->m_halfEdges.counter(); i++)
    {
        if (sim->m_halfEdges[i].m_inUse)
        {
            sim->m_halfEdges[i].updateLengthandCenter();
        }
    }
}

// calculate and set the length of every half edge
void sim_core::updateHalfEdgeLengths(Sim* sim)
{
    CO_PRINTF("\n---updating lengths of the half edges---\n\n");

    for (unsigned int i = 0; i < sim->m_halfEdges.counter(); i++)
    {
        if (sim->m_halfEdges[i].m_inUse)
        {
            sim->m_halfEdges[i].updateLength();
        }
    }
}

// delete the edges with zero length
void sim_core::delete_zero_length_edges(Sim* sim)
{
    updateHalfEdgeLengths(sim);

    CO_PRINTF("\n---deleting zero length edges---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].delete_zero_length_edges();
    }
}

// calculate the primitive fluid variables
void sim_core::calculatePrimitiveVariables(Sim* sim)
{
    CO_PRINTF("\n---calculating primitive fluid variables---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].calculatePrimitiveVariables();
    }
}

// calculate the conserved fluid variables
void sim_core::calculateConservedVariables(Sim* sim)
{
    CO_PRINTF("\n---calculating conserved fluid variables---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].calculateConservedVariables();
    }
}

// update the conserved variables
void sim_core::updateConservedVariables(Sim* sim)
{
    CO_PRINTF("\n---updating conserved fluid variables with the fluxes---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].updateConservedVariables();
    }
}

void sim_core::calculate_gradients(Sim* sim)
{
    CO_PRINTF("\n---calculating gradients---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].calculate_gradient(DENSITY);
        sim->m_particles[i].calculate_gradient(PRESSURE);
        sim->m_particles[i].calculate_gradient(VELOCITY_X);
        sim->m_particles[i].calculate_gradient(VELOCITY_Y);
    }
}

// slope limit the gradients
void sim_core::slope_limit_gradients_arepo(Sim* sim)
{
    CO_PRINTF("\n---slope limiting gradients (Arepo method)---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].slope_limit_gradient_arepo(DENSITY);
        sim->m_particles[i].slope_limit_gradient_arepo(PRESSURE);
        sim->m_particles[i].slope_limit_gradient_arepo(VELOCITY_X);
        sim->m_particles[i].slope_limit_gradient_arepo(VELOCITY_Y);
    }
}

void sim_core::slope_limit_gradients_tess(Sim* sim)
{
    CO_PRINTF("\n---slope limiting gradients (Tess method)---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].slope_limit_gradient_tess(DENSITY);
        sim->m_particles[i].slope_limit_gradient_tess(PRESSURE);
        sim->m_particles[i].slope_limit_gradient_tess(VELOCITY_X);
        sim->m_particles[i].slope_limit_gradient_tess(VELOCITY_Y);
    }
}

// update the velocities of the particles
void sim_core::updateParticleVelocitiesLagrangian(Sim* sim)
{
    CO_PRINTF("\n---updating velocities (Lagrangian)---\n\n");

    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];

        particle->m_particleVelocity = particle->m_velocity;

        // idea: reconstruct:
        //particle->m_particleVelocity.x = particle->m_velocity.x +
        //particle->m_gradientVelocityX.scalarProduct(particle->m_location - particle->m_cellCenterOfVolume);
        //particle->m_particleVelocity.y = particle->m_velocity.y +
        //particle->m_gradientVelocityY.scalarProduct(particle->m_location - particle->m_cellCenterOfVolume);
    }
}

// update the velocities of the particles in order to achieve mesh regularity
void sim_core::do_mesh_regulation_arepo_1(Sim* sim)
{
    CO_PRINTF("\n---mesh regulation (Arepo 1)---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].mesh_regulation_arepo_1();
    }
}

void sim_core::do_mesh_regulation_arepo_2(Sim* sim)
{
    CO_PRINTF("\n---mesh regulation (Arepo 2)---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].mesh_regulation_arepo_2();
    }
}

void sim_core::roundness_correction(Sim* sim)
{
    CO_PRINTF("\n---mesh regulation (havoc)---\n\n");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].roundness_correction();
    }
}

// move all Voronoi particles
void sim_core::moveVoronoiParticles(Sim* sim)
{
    if (sim->m_obstacle == NULL)
    {
        for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
        {
            sim->m_particles[i].move();

            if (!sim->m_box->isInside(&sim->m_particles[i].m_location))
            {
                sim->m_box->handle_outside_boundary(&sim->m_particles[i], i);
            }
        }
    }
    else
    {
        for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
        {
            sim->m_particles[i].move();

            if (!sim->m_box->isInside(&sim->m_particles[i].m_location))
            {
                sim->m_box->handle_outside_boundary(&sim->m_particles[i], i);
            }

            if (!sim->m_particles[i].m_obstacle_particle)
            {
                if (sim->m_obstacle->is_inside(&sim->m_particles[i].m_location))
                {
                    sim->m_obstacle->handle_inside_obstacle(&sim->m_particles[i]);
                }
            }
        }
    }
}

// reset the flag m_fluxUpdated of every half edge
void sim_core::reset_flux_update_flag(Sim* sim)
{
    for (unsigned int i = 0; i < sim->m_halfEdges.counter(); i++)
    {
        sim->m_halfEdges[i].m_fluxUpdated = false;
    }
}

// adaptive mesh refinement
void sim_core::amr_distribute_variables(Sim* sim)
{
    VoronoiParticle* particle = NULL;
    double c = 0; // ratio of the volumes
    double v1 = 0; // normalized volume of particle
    double v2 = 0; // normalized volume of twin particle

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];

        if (particle->m_amr_twin_particle != NULL)
        {
            c = particle->m_cellVolume / particle->m_amr_twin_particle->m_cellVolume;
            v1 = c / (c + 1);
            v2 = 1 / (c + 1);

            particle->m_amr_twin_particle->m_mass = v2 * particle->m_mass;
            particle->m_amr_twin_particle->m_energy = v2 * particle->m_energy;
            particle->m_amr_twin_particle->m_momentum = particle->m_momentum * v2;

            particle->m_mass *= v1;
            particle->m_energy *= v1;
            particle->m_momentum *= v1;

            particle->m_amr_twin_particle = NULL;
        }
    }
}

void sim_core::amr_refine(Sim* sim)
{
    VoronoiParticle* particle = NULL; // particle with big cell volume
    VoronoiParticle* vp = NULL; // new particle
    Vector particle_velocity_normalized = Vector(0, 0);

    const unsigned int nof_particles = sim->m_particles.counter();

    for (unsigned int i = 0; i < nof_particles; i++)
    {
        particle = &sim->m_particles[i];

        if (!particle->m_obstacle_particle &&
            particle->m_cellVolume > Constants::getMAX_CELL_VOLUME() * sim->m_box->m_DV)
        {
            vp = &sim->m_particles.current();
            *vp = VoronoiParticle(); //reset

            particle_velocity_normalized = particle->m_particleVelocity.return_normalized();

            if (particle_velocity_normalized.x != 0 || particle_velocity_normalized.y != 0)
            {
                vp->m_location = particle->m_location + particle->m_particleVelocity.return_normalized() * 0.000001;
            }
            else
            {
                vp->m_location = particle->m_location + Vector(myRandom::randomDouble(0, 1),
                                                               myRandom::randomDouble(0, 1)).return_normalized() *
                                                        0.000001;
            }

            sim->m_particles.counterInc();
            particle->m_amr_twin_particle = vp;
        }
    }
}

// calculate the values of primitive variables with point location and the gradient
//

double sim_core::point_density(Point p, Sim* sim, VoronoiParticle** start_and_found)
{
    *start_and_found = point_location::findCell(&p, sim, *start_and_found);

    if ((*start_and_found)->m_inner_obstacle_particle)
    {
        return 0;
    }

    return ((*start_and_found)->m_density +
            (*start_and_found)->m_gradientDensity.scalarProduct(p - (*start_and_found)->m_cellCenterOfVolume));
}

double sim_core::point_pressure(Point p, Sim* sim, VoronoiParticle** start_and_found)
{
    *start_and_found = point_location::findCell(&p, sim, *start_and_found);

    if ((*start_and_found)->m_inner_obstacle_particle)
    {
        return 0;
    }

    return ((*start_and_found)->m_pressure +
            (*start_and_found)->m_gradientPressure.scalarProduct(p - (*start_and_found)->m_cellCenterOfVolume));

}

double sim_core::point_v_x(Point p, Sim* sim, VoronoiParticle** start_and_found)
{
    *start_and_found = point_location::findCell(&p, sim, *start_and_found);

    if ((*start_and_found)->m_inner_obstacle_particle)
    {
        return 0;
    }

    return ((*start_and_found)->m_velocity.x +
            (*start_and_found)->m_gradientVelocityX.scalarProduct(p - (*start_and_found)->m_cellCenterOfVolume));
}

double sim_core::point_v_y(Point p, Sim* sim, VoronoiParticle** start_and_found)
{
    *start_and_found = point_location::findCell(&p, sim, *start_and_found);

    if ((*start_and_found)->m_inner_obstacle_particle)
    {
        return 0;
    }

    return ((*start_and_found)->m_velocity.y +
            (*start_and_found)->m_gradientVelocityY.scalarProduct(p - (*start_and_found)->m_cellCenterOfVolume));
}

// find outer half edges (debug)
void sim_core::printOuterHalfEdges(Sim* sim)
{
    for (unsigned int i = 0; i < sim->m_halfEdges.counter(); i++)
    {
        if (sim->m_halfEdges[i].m_inUse)
        {
            if (sim->m_halfEdges[i].m_twin == NULL)
            {
                printf("outer half edge: origin: %f %f\n", sim->m_halfEdges[i].m_origin->m_position.x,
                       sim->m_halfEdges[i].m_origin->m_position.y);
            }
        }
    }
}

// find particle indice (debug)
unsigned int sim_core::findParticleIndice(Sim* sim, double x, double y)
{
    unsigned int index_result = 0;
    double distance_square = std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        if (pow(sim->m_particles[i].m_location.x - x, 2) + pow(sim->m_particles[i].m_location.y - y, 2) <
            distance_square)
        {
            index_result = i;
            distance_square =
                    pow(sim->m_particles[i].m_location.x - x, 2) + pow(sim->m_particles[i].m_location.y - y, 2);
        }
    }

    printf("particle [%d]: %f %f\n", index_result, sim->m_particles[index_result].m_location.x,
           sim->m_particles[index_result].m_location.y);

    return index_result;
}

// find half edge indice (debug)
unsigned int sim_core::findHalfEdgeIndice(Sim* sim, double x, double y)
{
    unsigned int index_result = 0;
    double distance_square = std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < sim->m_halfEdges.counter(); i++)
    {
        if (pow(sim->m_halfEdges[i].m_center.x - x, 2) + pow(sim->m_halfEdges[i].m_center.y - y, 2) < distance_square)
        {
            index_result = i;
            distance_square = pow(sim->m_halfEdges[i].m_center.x - x, 2) + pow(sim->m_halfEdges[i].m_center.y - y, 2);
        }
    }

    unsigned int index_twin = 0;

    for (unsigned int j = 0; j < sim->m_halfEdges.counter(); j++)
    {
        if (&sim->m_halfEdges[j] == sim->m_halfEdges[index_result].m_twin)
        {
            index_twin = j;
            break;
        }
    }

    printf("half edge [%d]: %f %f\n", index_result, sim->m_halfEdges[index_result].m_center.x,
           sim->m_halfEdges[index_result].m_center.y);
    printf("half edge twin [%d]: %f %f\n", index_twin, sim->m_halfEdges[index_twin].m_center.x,
           sim->m_halfEdges[index_twin].m_center.y);

    return index_result;
}

// find ghost indice (debug)
unsigned int sim_core::findGhostIndice(Sim* sim, double x, double y)
{
    unsigned int index_result = 0;
    double distance_square = std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < sim->m_ghosts.counter(); i++)
    {
        if (pow(sim->m_ghosts[i].m_location.x - x, 2) + pow(sim->m_ghosts[i].m_location.y - y, 2) < distance_square)
        {
            index_result = i;
            distance_square = pow(sim->m_ghosts[i].m_location.x - x, 2) + pow(sim->m_ghosts[i].m_location.y - y, 2);
        }
    }

    printf("ghost [%d]: %f %f\n", index_result, sim->m_ghosts[index_result].m_location.x,
           sim->m_ghosts[index_result].m_location.y);

    return index_result;
}