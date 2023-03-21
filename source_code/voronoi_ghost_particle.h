#ifndef VORONOI_GHOST_PARTICLE_H
#define VORONOI_GHOST_PARTICLE_H

#include "voronoi_particle.h"

class GhostParticle : public VoronoiParticle
{
public:

    /// the corresponding real particle
    VoronoiParticle* m_realFriend;

    GhostParticle();

    /// is the particle a ghost particle?
    bool isGhost();

    /// resets a ghost particle, needed between two time steps
    void reset()
    {
        m_outerComponent = NULL;
        m_copiedOver = NULL;
        m_realFriend = NULL;
    }

    /// check whether the ghost cell is closed and has a real neighbour
    bool closed_and_inner();

    /// check whether the ghost cell is closed
    bool closed();

    /// update lengths and centers of the half edges
    void update_half_edges();
};

#endif