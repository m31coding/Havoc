#ifndef SIM_CORE
#define SIM_CORE

#include "sim.h"

namespace sim_core
{
    /// calculate the fluxes
    void calculateFluxes(Sim* sim);

    /// relax the initial data, the particles move to their centers, hexagons are formed
    void relax(Sim* sim, unsigned int steps);

    /// shift the outer components
    void shiftOuterComponents(Sim* sim);

    /// global time step
    extern double global_timestep;

    /// simulation time
    extern double simulation_time;

    /// time for the next binary output
    extern double next_output_time;

    /// is it time for a binary output?
    extern bool binary_output;

    /// output counter
    extern unsigned int counter;

    /// determine the constant global time step
    void constant_global_timestep(Sim* sim);

    /// determine the global time step (variable)
    void variable_global_timestep(Sim* sim);

    /// check whether it's time for a binary output
    void maybe_binary_output(Sim* sim, int step);

    /// calculate and set the volumes and the centers of volume of all cells
    void updateCellVolumesAndCenters(Sim* sim);

    /// calculate and set the length and the center of every half edge
    void updateHalfEdgeLengthsAndCenters(Sim* sim);

    /// calculate and set the length of every half edge
    void updateHalfEdgeLengths(Sim* sim);

    /// delete the edges with zero length
    void delete_zero_length_edges(Sim* sim);

    /// calculate the primitive fluid variables
    void calculatePrimitiveVariables(Sim* sim);

    /// calculate the conserved fluid variables (update cell volumes first!)
    void calculateConservedVariables(Sim* sim);

    /// update the conserved variables
    void updateConservedVariables(Sim* sim);

    /// calculate the gradients
    void calculate_gradients(Sim* sim);

    /// slope limit the gradients (Arepo version)
    void slope_limit_gradients_arepo(Sim* sim);

    /// slope limit the gradients (Tess version)
    void slope_limit_gradients_tess(Sim* sim);

    /// update the velocities of the particles
    void updateParticleVelocitiesLagrangian(Sim* sim);

    /// update the velocities of the particles in order to achieve mesh regularity
    void do_mesh_regulation_arepo_1(Sim* sim);

    /// update the velocities of the particles in order to achieve mesh regularity
    void do_mesh_regulation_arepo_2(Sim* sim);

    /// fix roundness
    void roundness_correction(Sim* sim);

    /// move all Voronoi particles
    void moveVoronoiParticles(Sim* sim);

    /// reset the flag m_fluxUpdated of every half edge
    void reset_flux_update_flag(Sim* sim);

    /// adaptive mesh refinement: distribute variables
    void amr_distribute_variables(Sim* sim);

    /// adaptive mesh refinement: refine cells
    void amr_refine(Sim* sim);

    /// calculate the density with point location and the gradient
    double point_density(Point p, Sim* sim, VoronoiParticle** start_and_found);

    /// calculate the pressure with point location and the gradient
    double point_pressure(Point p, Sim* sim, VoronoiParticle** start_and_found);

    /// calculate the x-velocity with point location and the gradient
    double point_v_x(Point p, Sim* sim, VoronoiParticle** start_and_found);

    /// calculate the y-velocity with point location and the gradient
    double point_v_y(Point p, Sim* sim, VoronoiParticle** start_and_found);

    /// find outer half edges (debug)
    void printOuterHalfEdges(Sim* sim);

    /// find particle indice (debug)
    unsigned int findParticleIndice(Sim* sim, double x, double y);

    /// find half edge indice (debug)
    unsigned int findHalfEdgeIndice(Sim* sim, double x, double y);

    /// find ghost indice (debug)
    unsigned int findGhostIndice(Sim* sim, double x, double y);

    /// dummy function, function pointer is set to this function if nothing has to be done
    inline void doNothing(Sim* sim)
    {}
};

#endif