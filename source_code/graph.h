#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <cstdio>
#include "r2.h"
#include "voronoi_particle.h"

// Computational Geometry, p.32: representing the graph as a doubly connected half edge list
// adapted to Voronoi, instead of creating a class face, a half edge stores a Voronoi particle

class Vertex;

class HalfEdge
{
public:

    bool m_inUse; ///< is the half edge in the array in use or not
    Vertex* m_origin;
    HalfEdge* m_twin; ///< the twin of the half edge
    VoronoiParticle* m_vp;
    HalfEdge* m_next;
    HalfEdge* m_prev;

    // geometric data of the edge
    double m_length;
    Point m_center;
    bool m_geoDataSet; ///< flag, tells whether m_length and m_center has already been calculated via the twin

    // flux
    bool m_fluxUpdated; ///< flag, tells whether the flux members have been set in the current time step.
    double m_massFlux;
    Vector m_momentumFlux;
    double m_energyFlux;

    HalfEdge(Vertex* origin, HalfEdge* twin, VoronoiParticle* vp, HalfEdge* next, HalfEdge* prev);

    HalfEdge();

    void updateLength(); ///< calculates and sets the length of the half edge
    void updateLengthandCenter(); ///< calculates and sets the length and the center of the half edge and its twin
    // warning: updates only real half edges and their twins

    void update_length_and_center_force(); /// update length and center without constrains

    void updateFlux(); ///< update the flux across the half edge
    // warning: updates only real half edges and their twins

    inline void set() ///< set member variables when using a new half edge
    {
        m_inUse = true;
        m_origin = NULL;
        m_twin = NULL;
        m_vp = NULL;
        m_next = NULL;
        m_prev = NULL;

        m_geoDataSet = false;
        m_fluxUpdated = false;
    }

    Vector velocity(); ///< calculate the velocity of the face (Springel, p. 805)

    bool is_zero(); ///< does the half edge have zero length?

    void info(FILE* pFile); ///< write data into a file
};

class Vertex
{
public:

    bool m_inUse; ///< is the half edge in the array in use or not
    Point m_position; ///< coordinates
    HalfEdge* m_edge; ///< incident Edge

    Vertex(const Point* position, HalfEdge* IncidentEdge);

    Vertex(const Point* position);

    Vertex();

    ~Vertex();

    inline void set(Point* position, HalfEdge* edge) ///< set all member variables when using a new vertex
    {
        m_inUse = true;
        m_position = *position;
        m_edge = edge;
    }

    void info(FILE* pFile); ///< write data into a file
};

namespace graph
{
    /// check whether the half face starting at start is closed
    bool isClosed(HalfEdge* start);

    /// print a face in gnuplot readable format
    void printFace(FILE* pFile, HalfEdge* start);

    /**
    print a half edge arrow, helper function for printHalfFace
    \sa printHalfFace
    */
    void printHalfEdge(FILE* pFile, const HalfEdge* halfEdge);

    /// print a half Face
    void printHalfFace(FILE* pFile, HalfEdge* start);
};
#endif