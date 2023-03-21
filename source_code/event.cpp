#include "config.h"
#include "event.h"
#include <cstdio>
#include "sim.h"
#include "debug.h"
#include <assert.h>
#include "constants.h"
#include "black_red_tree.h"
#include <cfloat>

// check event type
e_event PrintEvent::whatisit()
{
    return PRINT_EVENT;
}

e_event SiteEvent::whatisit()
{
    return SITE_EVENT;
}

e_event CircleEvent::whatisit()
{
    return CIRCLE_EVENT;
}

CircleEvent::CircleEvent(double intersection_x, double intersection_y, double x, double y, Arc* arc)
{
    m_intersection_x = intersection_x;
    m_intersection_y = intersection_y;
    m_x = x;
    m_y = y;
    m_arc = arc;
}

CircleEvent::CircleEvent()
{
    m_arc = NULL;
}

CircleEvent::~CircleEvent()
{

}

PrintEvent::PrintEvent(double x, double y)
{
    m_x = x;
    m_y = y;
}

PrintEvent::~PrintEvent()
{

}

SiteEvent::SiteEvent(VoronoiParticle* vp)
{
    m_vp = vp;
    m_x = vp->m_location.x;
    m_y = vp->m_location.y;
}

SiteEvent::SiteEvent()
{
    m_vp = NULL;
}

// handle a print event
void PrintEvent::handleEvent(Sim* sim)
{
    // debug
    SD_PRINTF("x: %f y: %f ", m_x, m_y);
    SD_PRINTF("(print event)\n");

    // print vector sites
    FILE* pFile;
    char filename[50];
    sprintf(filename, "../PrintEvent/%.2f_sites", Breakpoint::sweepline_y);

    pFile = fopen(filename, "w");
    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        fprintf(pFile, "%f %f\n", sim->m_particles[i].m_location.x, sim->m_particles[i].m_location.y);
    }
    fclose(pFile);

    // print the tree
    filename[0] = '\0';
    sprintf(filename, "../PrintEvent/%.2f_tree", Breakpoint::sweepline_y);
    //sprintf(filename,"PrintEvent/tree");
    pFile = fopen(filename, "w");
    tree::fPrintTreeKeys(pFile, sim->m_root, 0, 0, 1);
    fclose(pFile);

    filename[0] = '\0';
    sprintf(filename, "../PrintEvent/%.2f_tree.plt", Breakpoint::sweepline_y);
    pFile = fopen(filename, "w");
    fprintf(pFile, "unset border\nunset grid\nunset xtics\nunset ytics\n");
    fprintf(pFile, "plot \"%.2f_tree\" u 1:2:3 w lp ls 1 @orchid, \"%.2f_tree\" u 1:2:3 w labels offset 2",
            Breakpoint::sweepline_y, Breakpoint::sweepline_y);
    fclose(pFile);

    // print parabolas: use sites above the sweep line only, printed in Gnuplot readable format
    filename[0] = '\0';
    sprintf(filename, "../PrintEvent/%.2f_parabolas.plt", Breakpoint::sweepline_y);
    //sprintf(filename,"PrintEvent/parabolas.plt");
    pFile = fopen(filename, "w");

    fprintf(pFile, "sweep_y=%f\n", Breakpoint::sweepline_y);

    int unsigned n = 1;
    double offset = 5.0;
    double xrangemin = 0;
    double xrangemax = 0;
    double yrangemin = Breakpoint::sweepline_y - offset;
    double yrangemax = 0;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        if (sim->m_particles[i].m_location.y <= Breakpoint::sweepline_y)
        {
            continue;
        }
        else
        {
            if (sim->m_particles[i].m_location.x < xrangemin + offset)
            { xrangemin = sim->m_particles[i].m_location.x - offset; }
            if (sim->m_particles[i].m_location.x > xrangemax - offset)
            { xrangemax = sim->m_particles[i].m_location.x + offset; }
            if (sim->m_particles[i].m_location.y > yrangemax - offset)
            { yrangemax = sim->m_particles[i].m_location.y + offset; }

            fprintf(pFile, "x%d=%f\n", n, sim->m_particles[i].m_location.x);
            fprintf(pFile, "y%d=%f\n", n, sim->m_particles[i].m_location.y);
            fprintf(pFile, "f%d(x)=(y%d**2-sweep_y**2+x**2-2*x*x%d+x%d**2)/(2*y%d-2*sweep_y)\n", n, n, n, n, n);
        }

        n++;
    }

    fprintf(pFile, "set xrange[%f:%f]\n", xrangemin, xrangemax);
    fprintf(pFile, "set yrange[%f:%f]\n", yrangemin, yrangemax);

    fprintf(pFile, "set grid\n");

    fprintf(pFile, "plot sweep_y\\\n");

    for (unsigned int i = 1; i < n; i++)
    {
        fprintf(pFile, ",f%d(x)\\\n", i);

    }

    fprintf(pFile, ",\"%.2f_sites\"\\\n", Breakpoint::sweepline_y);
    //fprintf(pFile,",\"sites\"\\\n");
    fclose(pFile);
}

// helper function for Site Event and Circle Event
static inline bool handleTriple(const Point* p, const Point* q, const Point* r, Arc* arc, Point* intersection, Sim* sim)
{
    if (r2_geometry::breakpointsConverge(p, q, r, intersection)) // breakpoints converge and an intersection exists
    {
        double lowestPoint_y = r2_geometry::lowestPointY(intersection, q);

        if (lowestPoint_y <= Breakpoint::sweepline_y) // this would always be true with arbitrary precision arithmetic
        {
            CircleEvent* cEvent = new CircleEvent(intersection->x, intersection->y, 0, lowestPoint_y, arc);
            BR_Node::insert(cEvent, sim);

            arc->m_circleEvent = cEvent;

            SD_PRINTF("\t\tlowest point: %f\n", lowestPoint_y);

            return true;
        }
        else
        {
            fprintf(Constants::getP_WARNINGS_FILE(),
                    "WARNING: numerical inaccuracy in Event.cpp:handleTriple:, future circle event is above the sweepline\n");

            // set the circle event on the r2_geometry
            CircleEvent* cEvent = new CircleEvent(intersection->x, intersection->y, 0, Breakpoint::sweepline_y, arc);
            BR_Node::insert(cEvent, sim);

            arc->m_circleEvent = cEvent;

            SD_PRINTF("\t\tlowest point: %f\n", Breakpoint::sweepline_y);

            return true;
        }
    }
    else
    {
        return false;
    }
}

// handle site event
void SiteEvent::handleEvent(Sim* sim)
{
    // debug
    SD_PRINTF("x: %f y: %f ", m_vp->m_location.x, m_vp->m_location.y);
    SD_PRINTF("(site event)\n");

    // implementation follows computational geometry p.158

    // m_vp corresponds to p_i (new,lower site/Voronoi particle)
    // alpha->m_vp corresponds to p_j (upper site/Voronoi particle)

    // step 1
    //

    if (sim->m_root == NULL) // tree is empty
    {
        Node* arc = new Arc(m_vp, NULL);
        tree::TreeInsert(&(sim->m_root), arc);
        return;
    }

    // step 2
    //
    Arc* alpha = tree::searchArc(sim->m_root, m_vp); // alpha is the arc vertically above p_i

    // circle event is false alarm
    if (alpha->m_circleEvent != NULL)
    {
        delete BR_Node::remove(alpha->m_circleEvent, sim);
        alpha->m_circleEvent = NULL;

    }

    // step 3
    //

    // create a subtree
    Breakpoint* breakPointLeft = new Breakpoint();
    Breakpoint* breakPointRight = new Breakpoint();
    Arc* arcLeft = new Arc();
    Arc* arcMiddle = new Arc(); // new arc
    Arc* arcRight = new Arc();

    breakPointLeft->left = arcLeft;
    arcLeft->parent = breakPointLeft;

    breakPointLeft->right = breakPointRight;
    breakPointRight->parent = breakPointLeft;

    breakPointRight->left = arcMiddle;
    arcMiddle->parent = breakPointRight;

    breakPointRight->right = arcRight;
    arcRight->parent = breakPointRight;

    arcMiddle->m_vp = m_vp;
    arcLeft->m_vp = alpha->m_vp;
    arcRight->m_vp = alpha->m_vp;

    breakPointLeft->m_vpLeft = alpha->m_vp;
    breakPointLeft->m_vpRight = m_vp;

    breakPointRight->m_vpLeft = m_vp;
    breakPointRight->m_vpRight = alpha->m_vp;

    // substitute alpha with subtree

    if (alpha->parent == NULL) // alpha is the root
    {
        sim->m_root = breakPointLeft;
    }
    else // alpha has a parent
    {
        if (alpha->parent->right == alpha) // alpha is right child
        {
            alpha->parent->right = breakPointLeft;
            breakPointLeft->parent = alpha->parent;
        }
        else // alpha is left child
        {
            alpha->parent->left = breakPointLeft;
            breakPointLeft->parent = alpha->parent;
        }
    }

    // step 4
    //

    // create two half edges (one per breakpoint) and store them in the breakpoints

    HalfEdge* halfEdgeBPleft = &(sim->m_halfEdges.current());
    halfEdgeBPleft->set();
    sim->m_halfEdges.counterInc();

    HalfEdge* halfEdgeBPright = &(sim->m_halfEdges.current());
    halfEdgeBPright->set();
    sim->m_halfEdges.counterInc();

    breakPointLeft->m_halfEdge = halfEdgeBPleft;
    breakPointRight->m_halfEdge = halfEdgeBPright;

    // set twin and incident Voronoi particle (face)

    // the two breakpoints trace the same edge
    halfEdgeBPleft->m_twin = halfEdgeBPright;
    halfEdgeBPright->m_twin = halfEdgeBPleft;

    // the right Breakpoint moves to the right, therefore we store the upper site (Voronoi particle) in the half edge,
    // the lower site is stored in the half edge of the left Breakpoint, moving to the left

    halfEdgeBPright->m_vp = alpha->m_vp;
    halfEdgeBPleft->m_vp = m_vp;

    // store outer half edges in the Voronoi particles
    m_vp->m_outerComponent = halfEdgeBPleft;

    if (alpha->m_vp->m_outerComponent == NULL)
    {
        alpha->m_vp->m_outerComponent = halfEdgeBPright;
    }

    // step 5
    //

    // check the triple of consecutive arcs where the new arc for pi is the left/right arc to see if the breakpoints converge
    //
    Point intersection = Point(0, 0);

    // p_i is left arc
    // check whether the new arc has two successor arcs in the tree and whether the breakpoints converge

    Node* breakPointRightRight = tree::treeSuccessor(arcRight);

    if (breakPointRightRight != NULL) // a tree-successor exists (it is a breakpoint)
    {
        Arc* arcRightRight = static_cast<Arc*>(tree::treeSuccessor(breakPointRightRight));

        SD_PRINTF("\tnew arc is left arc:\n");
        SD_PRINTF("\tx1=%f y1=%f\n\tx2=%f y2=%f\n\tx3=%f y3=%f\n", arcMiddle->m_vp->m_location.x,
                  arcMiddle->m_vp->m_location.y, arcRight->m_vp->m_location.x, arcRight->m_vp->m_location.y,
                  arcRightRight->m_vp->m_location.x, arcRightRight->m_vp->m_location.y);

        if (handleTriple(&(arcMiddle->m_vp->m_location), &(arcRight->m_vp->m_location),
                         &(arcRightRight->m_vp->m_location), arcRight, &intersection, sim))
        {
            SD_PRINTF("\t\tintersection: x=%f y=%f\n", intersection.x, intersection.y);
        }
        else
        {
            SD_PRINTF("\t\tno intersection\n");
        }
    }

    Point intersection2 = Point(0, 0);

    // p_i is right arc
    // check whether the new arc has two predecessor arcs in the tree and whether the breakpoints converge

    Node* breakPointLeftLeft = tree::treePredecessor(arcLeft);

    if (breakPointLeftLeft != NULL)// a tree predecessor exists (it is a breakpoint)
    {
        Arc* arcLeftLeft = static_cast<Arc*>(tree::treePredecessor(breakPointLeftLeft));

        SD_PRINTF("\tnew arc is right arc:\n");
        SD_PRINTF("\tx1=%f y1=%f\n\tx2=%f y2=%f\n\tx3=%f y3=%f\n", arcLeftLeft->m_vp->m_location.x,
                  arcLeftLeft->m_vp->m_location.y, arcLeft->m_vp->m_location.x, arcLeft->m_vp->m_location.y,
                  arcMiddle->m_vp->m_location.x, arcMiddle->m_vp->m_location.y);

        if (handleTriple(&(arcLeftLeft->m_vp->m_location), &(arcLeft->m_vp->m_location), &(arcMiddle->m_vp->m_location),
                         arcLeft, &intersection2, sim))
        {
            SD_PRINTF("\t\tintersection: x=%f y=%f\n", intersection2.x, intersection2.y);
        }
        else
        {
            SD_PRINTF("\t\tno intersection\n");
        }
    }

    delete alpha;
}

// handle CircleEvent
void CircleEvent::handleEvent(Sim* sim)
{
    // debug
    SD_PRINTF("x: %f y: %f ", m_x, m_y);
    SD_PRINTF("(circle event)\n");

    // step 0
    //

    // skip circle events with false alarm
    // events with false alarm are deleted in the red black tree

    // step 1
    //

    // m_arc = alpha, the disappearing arc

    // store the adjacent breakpoints and arcs
    Breakpoint* prevBP = static_cast<Breakpoint*>(tree::treePredecessor(m_arc));
    Breakpoint* nextBP = static_cast<Breakpoint*>(tree::treeSuccessor(m_arc));

    Arc* nextArc = static_cast<Arc*>(tree::treeSuccessor(tree::treeSuccessor(m_arc)));
    Arc* prevArc = static_cast<Arc*>(tree::treePredecessor(tree::treePredecessor(m_arc)));

    // the center of the circle
    Point center = Point(m_intersection_x, m_intersection_y);
    SD_PRINTF("\tcenter: ( %f , %f )\n", center.x, center.y);

    // delete the circle events involving alpha
    if (nextArc->m_circleEvent != NULL)
    {
        delete BR_Node::remove(nextArc->m_circleEvent, sim);
        nextArc->m_circleEvent = NULL;
    }

    if (prevArc->m_circleEvent != NULL)
    {
        delete BR_Node::remove(prevArc->m_circleEvent, sim);
        prevArc->m_circleEvent = NULL;
    }

    // identify the remaining breakpoint (m_arc and father will be deleted) and update it
    Node* father = m_arc->parent;
    Breakpoint* remainingBP = NULL;

    if (father == prevBP)
    {
        remainingBP = nextBP;
    }
    else //father == nextBP
    {
        remainingBP = prevBP;
    }

    // step 2
    //

    // create half edge list entries

    HalfEdge* halfEdgePrevBP = prevBP->m_halfEdge;
    HalfEdge* halfEdgePrevBPtwin = halfEdgePrevBP->m_twin;

    HalfEdge* halfEdgeNextBP = nextBP->m_halfEdge;
    HalfEdge* halfEdgeNextBPtwin = halfEdgeNextBP->m_twin;

    // delete m_arc and father in this order
    delete tree::TreeDeleteOneChild(&(sim->m_root), m_arc);
    delete tree::TreeDeleteOneChild(&(sim->m_root), father);

    // use a new vertex
    sim->m_vertices.current().set(&center, halfEdgePrevBPtwin);

    // new edge
    HalfEdge* halfEdgeRemainingBP = &(sim->m_halfEdges.current());
    halfEdgeRemainingBP->set();
    sim->m_halfEdges.counterInc();

    HalfEdge* halfEdgeRemainingBPtwin = &(sim->m_halfEdges.current());
    halfEdgeRemainingBPtwin->set();
    sim->m_halfEdges.counterInc();

    // set pointers appropriately
    halfEdgeRemainingBP->m_twin = halfEdgeRemainingBPtwin;
    halfEdgeRemainingBPtwin->m_twin = halfEdgeRemainingBP;

    // m_origin
    halfEdgePrevBPtwin->m_origin = &(sim->m_vertices.current());
    halfEdgeNextBPtwin->m_origin = &(sim->m_vertices.current());
    halfEdgeRemainingBP->m_origin = &(sim->m_vertices.current());

    // jump to next unused vertex
    sim->m_vertices.counterInc();

    // m_next/m_prev
    halfEdgePrevBP->m_next = halfEdgeNextBPtwin;
    halfEdgeNextBPtwin->m_prev = halfEdgePrevBP;

    halfEdgeRemainingBPtwin->m_next = halfEdgePrevBPtwin;
    halfEdgePrevBPtwin->m_prev = halfEdgeRemainingBPtwin;

    halfEdgeNextBP->m_next = halfEdgeRemainingBP;
    halfEdgeRemainingBP->m_prev = halfEdgeNextBP;

    // set m_vp of the half edges
    halfEdgeRemainingBP->m_vp = nextArc->m_vp;
    halfEdgeRemainingBPtwin->m_vp = prevArc->m_vp;

    // update the remaining breakpoint
    remainingBP->m_vpLeft = prevArc->m_vp;
    remainingBP->m_vpRight = nextArc->m_vp;
    remainingBP->m_halfEdge = halfEdgeRemainingBP;

    // step 3
    //

    // check triples of consecutive arcs
    //

    Point intersection = Point(0, 0);

    // nextArc is middle arc

    Node* breakPointNextNext = tree::treeSuccessor(nextArc);

    if (breakPointNextNext != NULL) // a tree-successor exists (it is a breakpoint)
    {
        Arc* arcNextNext = static_cast<Arc*>(tree::treeSuccessor(breakPointNextNext));

        SD_PRINTF("\tnext arc to the vanishing arc is the middle arc:\n");
        SD_PRINTF("\tx1=%f y1=%f\n\tx2=%f y2=%f\n\tx3=%f y3=%f\n", prevArc->m_vp->m_location.x,
                  prevArc->m_vp->m_location.y, nextArc->m_vp->m_location.x, nextArc->m_vp->m_location.y,
                  arcNextNext->m_vp->m_location.x, arcNextNext->m_vp->m_location.y);

        if (handleTriple(&(prevArc->m_vp->m_location), &(nextArc->m_vp->m_location), &(arcNextNext->m_vp->m_location),
                         nextArc, &intersection, sim))
        {
            SD_PRINTF("\t\tintersection: x=%f y=%f\n", intersection.x, intersection.y);
        }
        else
        {
            SD_PRINTF("\t\tno intersection\n");
        }
    }

    // prevArc is middle arc

    Node* breakPointPrevPrev = tree::treePredecessor(prevArc);

    if (breakPointPrevPrev != NULL) // a tree predecessor exists (it is a breakpoint)
    {
        Arc* arcPrevPrev = static_cast<Arc*>(tree::treePredecessor(breakPointPrevPrev));

        SD_PRINTF("\tpreview arc to the vanishing arc is the middle arc:\n");
        SD_PRINTF("\tx1=%f y1=%f\n\tx2=%f y2=%f\n\tx3=%f y3=%f\n", arcPrevPrev->m_vp->m_location.x,
                  arcPrevPrev->m_vp->m_location.y, prevArc->m_vp->m_location.x, prevArc->m_vp->m_location.y,
                  nextArc->m_vp->m_location.x, nextArc->m_vp->m_location.y);

        if (handleTriple(&(arcPrevPrev->m_vp->m_location), &(prevArc->m_vp->m_location), &(nextArc->m_vp->m_location),
                         prevArc, &intersection, sim))
        {
            SD_PRINTF("\t\tintersection: x=%f y=%f\n", intersection.x, intersection.y);
        }
        else
        {
            SD_PRINTF("\t\tno intersection\n");
        }
    }
}