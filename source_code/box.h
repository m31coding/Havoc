#ifndef BOX_H
#define BOX_H

#include "r2.h"
#include "graph.h"
#include "my_array.h"

/// abstract base class
class Box
{
public:

    virtual ~Box()
    {}

    /// corresponding simulation
    Sim* m_sim;

    /// distance of particles
    double m_dh;
    double m_dv;
    double m_DV;

    /// check whether a point is inside of the box
    virtual bool isInside(Point* point) = 0;

    /// determine ghost particles
    virtual void determineGhostParticles() = 0;

    /// create ghost particles out of particles close to the border
    virtual void create_ghost_particles(double distance) = 0;

    /// create ghost particles, previous determined via determineGhostParticles
    virtual void createGhostParticles() = 0;

    /// copy the primitive variables to the ghost particles
    virtual void copyPrimitiveVariablesCellCentersAndVolume() = 0;

    /// copy the gradients to the ghost particles
    virtual void copyGradients() = 0;

    /// assign velocities to the ghost particles
    virtual void copyParticleVelocities() = 0;

    /// handle a particle which is outside of the boundary
    virtual void handle_outside_boundary(VoronoiParticle* vp, unsigned int& iterator) = 0;

    /// move the ghost particles
    virtual void move_ghost_particles() = 0;

    /// reset arrays
    virtual void reset() = 0;

    /// ini functions
    virtual void sitesIniGrid(unsigned int particles_x_dimension, unsigned int particles_y_dimension) = 0;
};

#endif