#include "config.h"
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include "graph.h"
#include "binary_search_tree.h"
#include "event.h"
#include "sim.h"
#include "r2.h"
#include "voronoi_particle.h"
#include "sim_graph.h"
#include "my_random.h"
#include "debug.h"
#include "constants.h"
#include "sweepline.h"
#include "sim_core.h"
#include "riemann_solver.h"
#include "time.h"
#include "info.h"
#include "options.h"
#include "binary_input_output.h"
#include "tools.h"
#include "sim_ini.h"

// analyze memory leaks
//#include <crtdbg.h>

/**
	main function	
	\param argc number of arguments
	\param argv arguments 
	\return error code
*/
int main(int argc, char* argv[])
{
    // parse options
    options::parse_options(argc, argv);

    // check options
    options::check_options(argc, argv);

    // use a havoc tool if the option is set
    if (tools::use_a_tool())
    {
        return EXIT_SUCCESS;
    }

    // read constants from file
    Constants::readConstants();

    // write into the info file
    info::print_header();
    info::print_version();
    info::print_args(argc, argv);
    Constants::fprint(Constants::getP_INFO_FILE());

    // check and adapt constants
    Constants::check_and_adapt();

    // initialize the Riemann solver
    exact_riemann_solver::ini();

    // initialize the random number generator via a seed
    myRandom::random_ini();

    // create a new simulation
    Sim* sim_1 = new Sim((Constants::getNOF_PARTICLES_X() * Constants::getNOF_PARTICLES_Y()) +
                         nof_obstacle_particles + nof_amr_particles);

    // create simulation boundary
    sim_1->create_box();

    // initialize sites
    sim_1->m_box->sitesIniGrid(Constants::getNOF_PARTICLES_X(), Constants::getNOF_PARTICLES_Y());
    //sim_1->sitesIniFile("../sites");
    //sim_1->sitesIniRandomDouble(-0.5, -0.5, 0.5, 0.5, 25*25);
    //sim_1->sitesIniPolarGrid(50, 200, 1. / sqrt(2));
    //ReflectingRectangle* box = static_cast<ReflectingRectangle*>(sim_1->m_box);
    //box->sitesIniGrid_diagonal(100, 10);

    // debug
    //sim_graph::printSites(sim_1,"../graph/sites");

    // create ghost particles
    sim_1->m_box->create_ghost_particles(0.1);

    // create an obstacle (circle, not moving)
    if (Constants::getOBSTACLE())
    {
        sim_1->create_obstacle();
    }

    // print general info
    info::print_general(sim_1);

    // initialize the event queue and calculate the tesselation
    sim_1->eventQueueIni();
    sim_1->eventQueueIniGhosts();
    sweepline::makeTesselation(sim_1);

    // debug
    //sim_graph::print_graph_all(sim_1);

    sim_graph::assertRealCellsClosed(sim_1);
    (Constants::get()->function_delete_zero_length_edges())(sim_1);
    sim_core::updateHalfEdgeLengthsAndCenters(sim_1);
    sim_graph::assertRealCellsClosed(sim_1);

    // relax the first graph
    if (Constants::getRELAX_FIRST_GRAPH())
    {
        sim_core::relax(sim_1, 1000);
    }

    // set the initial conditions
    sim_ini::initialize(sim_1);

    // binary input
    if (options::input_filename != NULL)
    {
        sim_1->freeGraph();
        delete sim_1;
        sim_1 = bin_io::binaryInput(options::input_filename);
        sim_1->create_box();

        // initialize the event queue and calculate the tesselation
        sim_1->eventQueueIni();
        sim_1->eventQueueIniGhosts();
        sweepline::makeTesselation(sim_1);
        sim_graph::assertRealCellsClosed(sim_1);
        (Constants::get()->function_delete_zero_length_edges())(sim_1);
        sim_core::updateHalfEdgeLengthsAndCenters(sim_1);
        sim_graph::assertRealCellsClosed(sim_1);
        sim_core::updateCellVolumesAndCenters(sim_1); // only real particles
        sim_core::calculateConservedVariables(sim_1); // only real particles
    }

    // continuation of a simulation
    if (options::input_filename != NULL)
    {
        sim_core::next_output_time = sim_core::simulation_time + Constants::getDELTA_T_OUTPUT();
        sim_core::binary_output = false;
    }
    else // new simulation
    {
        sim_core::next_output_time = 0;
        sim_core::binary_output = true;
    }

    // start time measurement
    time_t start_time = time(NULL);

    // note: must generate a tesselation out of particles and ghost particles and
    // calculate conserved variables before the main loop
    unsigned int max_steps = Constants::getMAX_STEPS();
    double max_time = Constants::getMAX_T() + 0.001;
    bool euler = (Constants::getVELOCITY_UPDATE_FUNCTION() == 0);

    // loop begin
    unsigned int k = 0;
    for (; k < max_steps && (sim_core::simulation_time <= max_time); k++)
    {
        // print step number and starting time of the time step
        printf("step: %u   time: %f", k, sim_core::simulation_time);

        if (!euler || k == 0)
        {
            // determine new ghost particles (old ghost particles needed)
            sim_1->m_box->determineGhostParticles();

            // reset data
            sim_1->freeGraph();

            // create new ghost particles
            sim_1->m_box->createGhostParticles();
            sim_1->m_box->reset();

            // initialize the event queue and calculate the tesselation
            sim_1->eventQueueIni();
            sim_1->eventQueueIniGhosts();
            sweepline::makeTesselation(sim_1);
            // note: for the ghost particles the attributes are copied, not calculated

            // optionally delete the edges with zero length
            (Constants::get()->function_delete_zero_length_edges())(sim_1);

            // delete zero length edges of the obstacle
            if (sim_1->m_obstacle != NULL)
            {
                sim_1->m_obstacle->delete_zero_length_edges();
            }

            // update the half edges (only real half edges and their twins)
            sim_core::updateHalfEdgeLengthsAndCenters(sim_1);

            // update the Voronoi particles (only real particles)
            sim_core::updateCellVolumesAndCenters(sim_1);
        }

        if (euler)
        {
            // reset the flag m_fluxUpdated of every half edge (half edges are reused)
            sim_core::reset_flux_update_flag(sim_1);
        }

        // distribute the conserved variables in refined cells
        if (Constants::getADAPTIVE_MESH_REFINEMENT())
        {
            sim_core::amr_distribute_variables(sim_1);
        }

        // calculate the primitive fluid variables and the sound speed in every cell (only real cells)
        sim_core::calculatePrimitiveVariables(sim_1);

        // check whether it's time for a binary output
        sim_core::maybe_binary_output(sim_1, k);

        // copy primitive variables, cell centers and volume to the ghost cells
        sim_1->m_box->copyPrimitiveVariablesCellCentersAndVolume();

        // copy primitive variables, cell centers and volume to the inner obstacle cells
        if (sim_1->m_obstacle != NULL)
        {
            sim_1->m_obstacle->copyPrimitiveVariablesCellCentersAndVolume();
        }

        // calculate gradients
        sim_core::calculate_gradients(sim_1);

        // slope limiter
        (Constants::get()->function_slopeLimitGradients())(sim_1);

        // copy gradients to ghost cells
        sim_1->m_box->copyGradients();

        // copy gradients to the inner obstacle cells
        if (sim_1->m_obstacle != NULL)
        {
            sim_1->m_obstacle->copyGradients();
        }

        // assign velocities to the Voronoi particles (only real particles)
        (Constants::get()->function_updateParticleVelocities())(sim_1);

        // add velocities in order to relax the Voronoi cells (only real particles)
        (Constants::get()->function_regulateMesh())(sim_1);

        // add velocities in order to do a correction of roundness
        (Constants::get()->function_roundness_correction())(sim_1);

        // assign velocities to the ghost particles
        sim_1->m_box->copyParticleVelocities();

        // assign velocities to the obstacle particles
        if (sim_1->m_obstacle != NULL)
        {
            sim_1->m_obstacle->assign_velocities();
        }

        // determine the global time step (only real cells are considered)
        (Constants::get()->function_determine_global_timestep())(sim_1);

        // adjust the simulation time
        sim_core::simulation_time += sim_core::global_timestep;

        // print end time of time step
        printf(" => %f\n", sim_core::simulation_time);

        // calculate the fluxes of all half edges (only real half edges and their twins)
        sim_core::calculateFluxes(sim_1);

        // update the conserved variables (only real particles)
        sim_core::updateConservedVariables(sim_1);

        // move the Voronoi particles (only real particles)
        sim_core::moveVoronoiParticles(sim_1);

        // move the ghost particles (needed for zero gradient boundaries)
        sim_1->m_box->move_ghost_particles();

        // refine the mesh if necessary
        if (Constants::getADAPTIVE_MESH_REFINEMENT())
        {
            sim_core::amr_refine(sim_1);
        }
    }

    // free memory of the graph
    sim_1->freeGraph();

    // free memory
    delete sim_1;

    // duration of the simulation
    unsigned long int duration = time(NULL) - start_time;

    // print duration into the info file
    info::print_duration(duration, k - 1);

    // free memory
    Constants::del();
    options::clear();

    printf("\ndone\n");

    return EXIT_SUCCESS;
}