#ifndef EVENT_H
#define EVENT_H

#include <queue>
#include "graph.h"
#include "binary_search_tree.h"
#include "r2.h"
#include "voronoi_particle.h"
#include "black_red_tree.h"

// different event types
enum e_event
{
    SITE_EVENT, CIRCLE_EVENT, PRINT_EVENT
};

class CompareEvents;

class Sim;

// events in PQ

/**
Events are stored in the priority queue and handled if the sweep line meets the y coordinate of m_point.
*/
class Event : public BR_Node
{
public:

    Event()
    {};

    virtual ~Event()
    {};

    virtual void handleEvent(Sim* sim) = 0;

    virtual e_event whatisit() = 0;
};

class SiteEvent : public Event
{
public:

    VoronoiParticle* m_vp; ///< stores a site in m_vp->location

    //constructor/destructor
    SiteEvent(VoronoiParticle* vp);

    SiteEvent();

    ~SiteEvent()
    {}

    // handle site event
    void handleEvent(Sim* sim);

    // check type
    e_event whatisit();
};

class CircleEvent : public Event
{
public:

    // stores the lowest point of the circle in m_vp->location
    Arc* m_arc;
    double m_intersection_x;
    double m_intersection_y;

    CircleEvent(double intersection_x, double intersection_y, double x, double y, Arc* arc);

    CircleEvent();

    ~CircleEvent();

    // handle circle event
    void handleEvent(Sim* sim);

    // check type
    e_event whatisit();
};

/// An Event for debugging and creating plots for Gnuplot
class PrintEvent : public Event
{
public:

    PrintEvent(double x, double y);

    ~PrintEvent();

    //handle print event
    void handleEvent(Sim* sim);

    //check type
    e_event whatisit();
};

#endif