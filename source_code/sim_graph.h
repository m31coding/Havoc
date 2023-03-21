#ifndef SIMGRAPH_H
#define SIMGRAPH_H

#include "sim.h"

/// functions operating on the graph of a simulation
namespace sim_graph
{
    /// delete a half edge
    void delete_edge(HalfEdge* half_edge);

    /// paint the cells with random colors
    void randomColors(Sim* sim);

    /// check whether all real cells are closed
    void assertRealCellsClosed(Sim* sim);

    /// print the sites
    void printSites(Sim* sim, const char* filename);

    /// print the positions of the ghost particles
    void printGhosts(Sim* sim, const char* filename);

    /// print the completed faces in gnuplot readable format
    void printFaces(Sim* sim, const char* filename);

    /// print the ghost faces of the graph
    void printGhostFaces(Sim* sim, const char* filename);

    /// print the completed faces with half edges in gnuplot readable format
    void printHalfFaces(Sim* sim, const char* filename);

    /// print data of the half edges
    void printHalfEdgesData(Sim* sim, const char* filename);

    /// print the vertices
    void printVertices(Sim* sim, const char* filename);

    /// print the breakpoint vertices
    void printBreakPointVertices(Sim* sim, const char* filename);

    /// print the edges adjacent to the breakpoint vertices
    void printBreakPointEdges(Sim* sim, const char* filename);

    /// print the face centers
    void printFaceCenters(Sim* sim, const char* filename);

    /// print the graph for gnuplot pm3d map splot
    void printGraphColor(Sim* sim, const char* filename);

    /// print the whole graph
    void print_graph_all(Sim* sim);
};

#endif