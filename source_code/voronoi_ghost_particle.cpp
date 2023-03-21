#include "config.h"
#include "voronoi_ghost_particle.h"
#include "graph.h"


GhostParticle::GhostParticle()
        : m_realFriend(NULL)
{

}

// is the particle a ghost particle?
bool GhostParticle::isGhost()
{
    return true;
}

// check whether the ghost cell is closed and has a real neighbour
bool GhostParticle::closed_and_inner()
{
    bool real_neighbour = false;

    HalfEdge* start = m_outerComponent;
    HalfEdge* current = m_outerComponent;
    HalfEdge* next = NULL;

    do
    {
        next = current->m_next;
        if (next == NULL)
        { return false; } // ghost face not closed

        if (!real_neighbour)
        {
            real_neighbour = (!current->m_twin->m_vp->isGhost());
        }

        current = next;

    } while (current != start);

    return real_neighbour;
}

// check whether the ghost cell is closed
bool GhostParticle::closed()
{
    HalfEdge* start = m_outerComponent;
    HalfEdge* current = m_outerComponent;
    HalfEdge* next = NULL;

    do
    {
        next = current->m_next;
        if (next == NULL)
        { return false; } // ghost face not closed

        current = next;

    } while (current != start);

    return true;
}

// update lengths and centers of the half edges
void GhostParticle::update_half_edges()
{
    HalfEdge* start = m_outerComponent;
    HalfEdge* current = m_outerComponent;

    do
    {
        assert(current != NULL);

        current->update_length_and_center_force();
        current = current->m_next;

    } while (current != start);
}