#ifndef SWEEPLINE_H
#define SWEEPLINE_H

#include "sim.h"

namespace sweepline
{
    /// calculate the Voronoi tesselation via sweepline algorithm
    void makeTesselation(Sim* sim);

    /// handle first site and horizontal line degeneracy
    void handleHorizontalDegeneracy(Sim* sim);

    /**
    convert remaining breakpoints to vertices
    note: they are not inserted into m_vertices of the simulation but in sim->m_breakpointVertices
    */
    void breakpointsToVertices(Node* root, Sim* sim); // recursive
};

#endif