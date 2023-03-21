#ifndef VORONOI_PARTICLE_H
#define VORONOI_PARTICLE_H

#include <vector>
#include "r2.h"
#include <cstdio>
#include "my_array.h"
#include <cmath>
#include "constants.h"

class HalfEdge;

/// primitive variables
enum prim_var
{
    DENSITY, PRESSURE, VELOCITY_X, VELOCITY_Y
};

/// Voronoi particles are generators of the Voronoi tesselation and move according to the velocity field
class VoronoiParticle
{
public:

    Point m_location; ///< the location in R^2 of the Voronoi particle
    Vector m_particleVelocity; ///< the velocity of the Voronoi particle

    HalfEdge* m_outerComponent; ///< a pointer to an outer half edge of the face corresponding to the Voronoi particle

    Point m_cellCenterOfVolume; ///< the center of volume of the corresponding Voronoi cell
    double m_cellVolume; ///< the volume of the cell

    // primitive variables
    Vector m_velocity; ///< the velocity of the Voronoi particle
    double m_density; ///< the density of the Voronoi particle
    double m_pressure; ///< the pressure of the Voronoi particle

    // conserved fluid variables
    double m_mass; ///< the total mass contained in the cell
    Vector m_momentum; ///< the momentum contained in the cell
    double m_energy; ///< the energy contained in the cell

    // gradients
    Point m_gradientVelocityX;
    Point m_gradientVelocityY;
    Point m_gradientDensity;
    Point m_gradientPressure;

    VoronoiParticle* m_amr_twin_particle;

    double m_color; ///< color the particle in order to track it
    int m_flag_1; ///< flags

    bool m_obstacle_particle; ///< is the particle an obstacle particle
    bool m_inner_obstacle_particle; ///< is it an inner obstacle particle

    /// m_copied != NULL means, that the particle has been copied over a part of the boundary represented by a pointer
    myArray<VoronoiParticle*>* m_copiedOver;

    double& var(prim_var pv); ///< get a primitive variable
    Vector& grad(prim_var pv); ///< get a gradient

    /// calculate the speed of sound of the cell
    inline double soundSpeed()
    {
        assert(m_density != 0);
        return sqrt(Constants::getGAMMA() * m_pressure / m_density);
    }

    /// calculate the effective radius of the cell
    inline double radius()
    {
        return sqrt(m_cellVolume / M_PI);
    }

    double roundness(); ///< calculate the roundness of the cell
    double vorticity(); ///< calculate the vorticity of the cell
    double turbulence_kinetic_energy(); ///< calculate the turbulent kinetic energy of the cell

    /// calculate the divergence of the velocity
    double div_v()
    {
        return (m_gradientVelocityX.x + m_gradientVelocityY.y);
    }

    VoronoiParticle(); ///< default constructor
    VoronoiParticle(const Point* location); ///< constructor, sets m_location
    VoronoiParticle(double x, double y); ///< constructor, sets m_location

    virtual ~VoronoiParticle()
    {};

    double local_timestep(); ///< calculate the local time step of the cell
    bool virtual isGhost(); ///< is the particle a ghost particle?
    void
    updateCellVolumeAndCenter(); ///< calculate and set the volume and the center of volume of the corresponding cell
    void calculate_gradient(prim_var pv); /// calculate the gradient
    void slope_limit_gradient_arepo(prim_var pv); ///< slope limiter
    void slope_limit_gradient_tess(prim_var pv); ///< slope limiter
    void updateConservedVariables(); ///< update the conserved variables with the fluxes

    //update the velocities of the particles in order to achieve mesh regularity
    //

    void mesh_regulation_arepo_1(); ///< normal flows
    void mesh_regulation_arepo_2(); ///< cold flows
    void roundness_correction(); ///< additional correction

    void calculatePrimitiveVariables(); ///< calculate the primitive variables of the cell from the conserved variables
    void calculateConservedVariables(); ///< calculate the conserved variables of the cell from the primitive variables

    void delete_zero_length_edges(); ///< delete the half edges which have zero length

    void move(); ///< move the particle according to its velocity
    void move_backwards(); ///move the particle backward

    void info(FILE* pFile); ///< write stored data into a file (debug)
    void info_vertices(FILE* pFile); ///< write vertex positions into a file (debug)

    inline virtual void reset() ///< resets a particle, needed between two time steps
    {
        m_outerComponent = NULL;
        m_copiedOver = NULL;
    }

    int number_of_neighbours();
};

#endif