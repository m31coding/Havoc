#include "config.h"
#include "point_location.h"
#include "graph.h"
#include "debug.h"

// check whether a point is inside a cell
bool point_location::isInside(Point* point, VoronoiParticle* vp)
{
    assert(vp->m_outerComponent != NULL);

    HalfEdge* start = vp->m_outerComponent;
    HalfEdge* current = vp->m_outerComponent;

    Point p = Point(0, 0);
    Point q = Point(0, 0);

    do
    {
        p = current->m_origin->m_position;
        q = current->m_twin->m_origin->m_position;

        // right turn but with numerical criteria (see also r2_geometry::rightTurn)
        double result = (q.y - p.y) * (point->x - p.x) + (p.x - q.x) * (point->y - p.y);

        if (result > 0.000001)
        {
            return false;
        }

        current = current->m_next;

    } while (current != start);

    return true;
}

// find the cell which contains the point, starts searching at cell start
VoronoiParticle* point_location::findCell(Point* point, Sim* sim, VoronoiParticle* search_start)
{
    PL_PRINTF("point: ( %f , %f )\n", point->x, point->y);

    if (!sim->m_box->isInside(point))
    {
        fprintf(stderr, "ERROR in point_location::findCell: point is not inside of the box!\n");
        fprintf(stderr, "\tpoint: ( %f , %f )\n", point->x, point->y);
        exit(1);
    }

    VoronoiParticle* vp = NULL;

    if (search_start != NULL)
    {
        // start with Voronoi particle start
        vp = search_start;
    }
    else
    {
        // start with a Voronoi particle near the "middle"
        vp = &sim->m_particles[((int) (sim->m_particles.counter() * 0.5))];
    }

    double distance_square = (vp->m_location - *point).length_square();
    double new_distance_square = 0;

    // walk around the cell
    HalfEdge* start = NULL;
    HalfEdge* current = NULL;

    bool better = true;

    while (better)
    {
        better = false;

        start = vp->m_outerComponent;
        current = start;

        do
        {
            PL_PRINTF("\tcurrent cell: ( %f , %f )\n", current->m_vp->m_location.x, current->m_vp->m_location.y);

            new_distance_square = (current->m_twin->m_vp->m_location - *point).length_square();

            if (new_distance_square < distance_square)
            {
                distance_square = new_distance_square;
                vp = current->m_twin->m_vp;
                better = true;
            }

            current = current->m_next;

        } while (current != start);
    }

    PL_PRINTF("final cell: ( %f , %f )\n", current->m_vp->m_location.x, current->m_vp->m_location.y);

    if (!isInside(point, current->m_vp))
    {
        fprintf(stderr, "ERROR in point_location::findCell: point location failed for point %f	%f!\n", point->x,
                point->y);
        fprintf(stderr, "found %f	%f\n", current->m_vp->m_location.x, current->m_vp->m_location.y);
        exit(1);
    }

    if (current->m_vp->isGhost())
    {
        fprintf(stderr, "ERROR in point_location::findCell: cell is a ghost\n");
        exit(1);
    }

    return current->m_vp;
}