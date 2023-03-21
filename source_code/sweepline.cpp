#include "config.h"
#include "sweepline.h"
#include "binary_search_tree.h"
#include "event.h"
#include "debug.h"
#include "assert.h"
#include "constants.h"

static GhostParticle top_ghost;

// calculate the Voronoi tesselation with the sweep line algorithm
void sweepline::makeTesselation(Sim* sim)
{
    SD_PRINTF("\n---sweep line algorithm---\n\n");

    // stop the algorithm when all parabolas are outside the box
    double finalSweeplinePosition = sim->getSweepY();

    // handle horizontal line degeneracy
    //handleHorizontalDegeneracy(sim);

    // add a single site in order to avoid horizontal degeneracy
    top_ghost.m_location.x = 0.5 * (sim->getCDXmin() + sim->getCDXmax());
    top_ghost.m_location.y = 2 * sim->getCDYmax() - sim->getSweepY() + 0.1;
    Event* p_siteEvent = new SiteEvent(&top_ghost);
    BR_Node::insertSite(p_siteEvent, sim);

    BR_Node* delMe;

    // main loop
    Event* p_Event = NULL;
    while (sim->m_BR_root != BR_Node::p_nil)
    {
        // take the first Event
        p_Event = static_cast<Event*>(BR_Node::treeMaximum(sim->m_BR_root));
        delMe = BR_Node::remove(p_Event, sim);

        Breakpoint::sweepline_y = p_Event->m_y;

        // stop the algorithm when all parabolas are outside the box
        if (Breakpoint::sweepline_y < finalSweeplinePosition)
        {
            delete delMe;
            break;
        }

        // handle the Event (polymorphic)
        p_Event->handleEvent(sim);

        // free memory
        delete delMe;
    }

    // set the sweepline position
    Breakpoint::sweepline_y = finalSweeplinePosition;
    breakpointsToVertices(sim->m_root, sim);

    GN_PRINTF("\nm_particles counter (faces): %u\n", sim->m_particles.counter());
    GN_PRINTF("m_ghosts counter: %u\n", sim->m_ghosts.counter());
    GN_PRINTF("ratio: %f\n", (double) sim->m_particles.counter() / sim->m_ghosts.counter());
    GN_PRINTF("m_vertices counter: %u\n", sim->m_vertices.counter());
    GN_PRINTF("number of active vertices: %u\n", sim->m_vertices.NOFactive());
    GN_PRINTF("ratio: %f\n", (((double) sim->m_vertices.counter()) / sim->m_vertices.NOFactive()));
    GN_PRINTF("m_halfEdges counter: %u\n", sim->m_halfEdges.counter());
    GN_PRINTF("number of active half edges: %u\n", sim->m_halfEdges.NOFactive());
    GN_PRINTF("ratio: %f\n", ((double) sim->m_halfEdges.counter()) / sim->m_halfEdges.NOFactive());
}

// handle first site and horizontal line degeneracy
void sweepline::handleHorizontalDegeneracy(Sim* sim)
{
    // there should be at least two Voronoi particles
    assert(sim->m_BR_root != BR_Node::p_nil &&
           (sim->m_BR_root->m_left != BR_Node::p_nil || sim->m_BR_root->m_right != BR_Node::p_nil));

    double y_coord = 0;
    SiteEvent* p_Event = NULL;
    BR_Node* max = NULL;
    BR_Node* delMe = NULL;

    // handle the first site
    //

    max = BR_Node::treeMaximum(sim->m_BR_root);

    if (static_cast<Event*>(max)->whatisit() != SITE_EVENT)
    {
        // first event is print event, change print events
        fprintf(stderr, "ERROR in sweepline::handleHorizontalDegeneracy: first event is no site event\n");
        return;
    }

    p_Event = static_cast<SiteEvent*>(max);
    delMe = BR_Node::remove(max, sim);

    // save the y coordinate of the first site
    y_coord = p_Event->m_y;

    // debug
    SD_PRINTF("x: %f y: %f ", p_Event->m_x, p_Event->m_y);
    SD_PRINTF("(first site event)\n");

    // insert first site into the tree
    assert(sim->m_root == NULL); //tree should be empty

    Arc* arc = new Arc(p_Event->m_vp, NULL);
    tree::TreeInsert(&(sim->m_root), arc);

    // handle all next sites with the same y coordinate as the first site
    //

    // store the previous Voronoi particle
    VoronoiParticle* oldParticle = p_Event->m_vp;

    delete delMe;

    max = BR_Node::treeMaximum(sim->m_BR_root);

    if (static_cast<Event*>(max)->whatisit() != SITE_EVENT)
    {
        // second event is print event and therefore the next site event has not the same y-coord as the first site event
        return;
    }

    assert(static_cast<Event*>(max)->whatisit() == SITE_EVENT);
    p_Event = static_cast<SiteEvent*>(max);

    while ((sim->m_BR_root != BR_Node::p_nil) && (p_Event->m_vp->m_location.y == y_coord))
    {
        // debug
        SD_PRINTF("x: %f y: %f ", p_Event->m_x, p_Event->m_y);
        SD_PRINTF("(site event horizontally degenerated)\n");

        delMe = BR_Node::remove(max, sim);

        // create a new breakpoint
        Breakpoint* newBreakpoint = new Breakpoint;

        // create a new arc
        Arc* newArc = new Arc;

        // create half edges
        HalfEdge* newHe = &(sim->m_halfEdges.current());
        newHe->m_inUse = true;
        sim->m_halfEdges.counterInc();

        HalfEdge* newHeTwin = &(sim->m_halfEdges.current());
        newHeTwin->m_inUse = true;
        sim->m_halfEdges.counterInc();

        // set the member variables of the new nodes
        newBreakpoint->m_vpLeft = oldParticle;
        newBreakpoint->m_vpRight = p_Event->m_vp;
        newArc->m_vp = p_Event->m_vp;
        newBreakpoint->m_halfEdge = newHe;

        // insert the new nodes into the tree
        newBreakpoint->left = sim->m_root;
        sim->m_root->parent = newBreakpoint;
        sim->m_root = newBreakpoint;

        newBreakpoint->right = newArc;
        newArc->parent = newBreakpoint;

        // set the member variables of the half edges
        newHe->m_twin = newHeTwin;
        newHeTwin->m_twin = newHe;

        newHe->m_vp = p_Event->m_vp;
        newHeTwin->m_vp = oldParticle;

        // set outer components of the Voronoi particle
        p_Event->m_vp->m_outerComponent = newHe;
        oldParticle->m_outerComponent = newHeTwin;

        // create a vertex outside the box
        Vertex* upperVertex = new Vertex();
        sim->m_breakpointVertices.push_back(upperVertex);

        double ymax = 0;

        ymax = sim->getCDYmax();

        upperVertex->m_position.x = (oldParticle->m_location.x + p_Event->m_vp->m_location.x) / 2;
        upperVertex->m_position.y = ymax * 2;

        upperVertex->m_edge = newHe;
        newHe->m_origin = upperVertex;

        // update and free memory
        oldParticle = p_Event->m_vp;
        delete delMe;

        // return if only two sites exist
        if (sim->m_BR_root == BR_Node::p_nil)
        {
            return;
        }

        max = BR_Node::treeMaximum(sim->m_BR_root);

        // skip possible print event
        if (static_cast<Event*>(max)->whatisit() == PRINT_EVENT)
        {
            return;
        }

        assert(static_cast<Event*>(max)->whatisit() == SITE_EVENT);
        p_Event = static_cast<SiteEvent*>(max);
    }
}


// convert a remaining breakpoint to a vertex (helper function for breakpointsToVertices)
static void breakpointToVertex(const Breakpoint* breakpoint, Sim* sim)
{
    double breakpoint_x = breakpoint->xcoord();
    Point P_Breakpoint(breakpoint_x, r2_geometry::breakp_y(&(breakpoint->m_vpLeft->m_location),
                                                           breakpoint_x, Breakpoint::sweepline_y));

    // the breakpoint vertices have to be outside the computational domain
    assert(!sim->isInsideComputationalDomain(&P_Breakpoint));

    // create a new vertex
    Vertex* vertex = new Vertex(&P_Breakpoint, breakpoint->m_halfEdge->m_twin);
    breakpoint->m_halfEdge->m_twin->m_origin = vertex;
    sim->m_breakpointVertices.push_back(vertex);
}

// convert remaining breakpoints to vertices
void sweepline::breakpointsToVertices(Node* root, Sim* sim)
{
    if (root == NULL)
    { return; }

    if (!root->isLeaf()) // the node is a breakpoint
    {
        const Breakpoint* breakpoint = static_cast<const Breakpoint*>(root);

        breakpointToVertex(breakpoint, sim);

        breakpointsToVertices(root->left, sim);
        breakpointsToVertices(root->right, sim);
    }

    return;
} // warning: recursive