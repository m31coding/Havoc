#ifndef POINT_LOCATION_H
#define POINT_LOCATION_H

#include "r2.h"
#include "voronoi_particle.h"
#include "sim.h"

namespace point_location
{
    /**
     check whether a point is inside a cell
     note: on the rim is also inside, this means that a point can be inside of multiple cells
     */
    bool isInside(Point* point, VoronoiParticle* vp);

    /// find the (ghost) cell which contains the point, start searching at cell start
    VoronoiParticle* findCell(Point* point, Sim* sim, VoronoiParticle* search_start);
};

#endif