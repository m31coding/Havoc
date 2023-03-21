#include "config.h"
#include "graph.h"
#include "debug.h"
#include "assert.h"
#include <cmath>
#include <stdlib.h>
#include "sim_core.h"
#include "riemann_solver.h"

HalfEdge::HalfEdge(Vertex* origin, HalfEdge* twin, VoronoiParticle* vp, HalfEdge* next, HalfEdge* prev)
        :
        m_inUse(false),
        m_origin(origin),
        m_twin(twin),
        m_vp(vp),
        m_next(next),
        m_prev(prev),
        m_length(0),
        m_center(0, 0),
        m_geoDataSet(false),
        m_fluxUpdated(false),
        m_massFlux(0),
        m_momentumFlux(0, 0),
        m_energyFlux(0)
{

}

HalfEdge::HalfEdge()
        :
        m_inUse(false),
        m_origin(NULL),
        m_twin(NULL),
        m_vp(NULL),
        m_next(NULL),
        m_prev(NULL),
        m_length(0),
        m_center(0, 0),
        m_geoDataSet(false),
        m_fluxUpdated(false),
        m_massFlux(0),
        m_momentumFlux(0, 0),
        m_energyFlux(0)
{

}

void HalfEdge::updateLengthandCenter()
{
    if (m_geoDataSet || m_vp->isGhost())
    {
        return;
    }

    Point a = m_origin->m_position;
    Point b = m_next->m_origin->m_position;
    // m_twin->m_origin if you don't want to ignore ghost particles

    // calculate the length
    double length = sqrt(pow(b.x - a.x, 2) + pow(b.y - a.y, 2));

    // calculate the center
    Point center = (a + b) * 0.5;

    // update
    m_length = length;
    m_center = center;
    m_geoDataSet = true;

    HE_PRINTF("\nhalf edge center: ( %f %f )\n", m_center.x, m_center.y);
    HE_PRINTF("\tlength: %f\n", length);
    HE_PRINTF("\tcenter: ( %f , %f )\n", center.x, center.y);

    assert(m_twin != NULL);

    m_twin->m_length = length;
    m_twin->m_center = center;
    m_twin->m_geoDataSet = true;

    return;
}

// update length and center without constraints
void HalfEdge::update_length_and_center_force()
{
    Point a = m_origin->m_position;
    Point b = m_next->m_origin->m_position;
    // m_twin->m_origin if you don't want to ignore ghost particles

    // calculate the length
    m_length = sqrt(pow(b.x - a.x, 2) + pow(b.y - a.y, 2));

    // calculate the center
    m_center = (a + b) * 0.5;

    m_twin->m_length = m_length;
    m_twin->m_center = m_center;
}

void HalfEdge::updateLength()
{
    if (m_vp->isGhost())
    {
        return;
    }

    Point a = m_origin->m_position;
    Point b = m_next->m_origin->m_position;
    // m_twin->m_origin if you don't want to ignore ghost particles

    // calculate the length
    double length = sqrt(pow(b.x - a.x, 2.) + pow(b.y - a.y, 2.));

    // update
    m_length = length;
    m_twin->m_length = length;

    HE_PRINTF("\tlength: %f\n", length);

    return;
}

// update the flux across the half edge
void HalfEdge::updateFlux()
{
    if (m_fluxUpdated || m_vp->isGhost())
    {
        return;
    }

    // velocity of the half edge
    Vector w = velocity();

    VoronoiParticle* vp_left = m_vp;
    VoronoiParticle* vp_right = m_twin->m_vp;

    if (vp_left->m_inner_obstacle_particle && vp_right->m_inner_obstacle_particle)
    {
        return;
    }

    // left state (on the left side of the half edge)
    double rho_left = vp_left->m_density;
    Vector v_left = vp_left->m_velocity;
    double p_left = vp_left->m_pressure;

    // right state
    double rho_right = vp_right->m_density;
    Vector v_right = vp_right->m_velocity;
    double p_right = vp_right->m_pressure;

    HE_PRINTF("states in the lab frame:\n");

    HE_PRINTF("\tleft state velocity: ( %f %f )\n", v_left.x, v_left.y);
    HE_PRINTF("\tright state velocity: ( %f %f )\n", v_right.x, v_right.y);

    // transform to the rest frame
    v_left -= w;
    v_right -= w;

    HE_PRINTF("states in the rest frame:\n");

    HE_PRINTF("left state: (%f %f )\n", vp_left->m_location.x, vp_left->m_location.y);
    HE_PRINTF("\tdensity: %f\n", rho_left);
    HE_PRINTF("\tvelocity: ( %f %f )\n", v_left.x, v_left.y);
    HE_PRINTF("\tpressure: %f\n", p_left);
    HE_PRINTF("\t\tgradient density: ( %f %f )\n", vp_left->m_gradientDensity.x, vp_left->m_gradientDensity.y);
    HE_PRINTF("\t\tgradient v_x: ( %f %f )\n", vp_left->m_gradientVelocityX.x, vp_left->m_gradientVelocityX.y);
    HE_PRINTF("\t\tgradient v_y: ( %f %f )\n", vp_left->m_gradientVelocityY.x, vp_left->m_gradientVelocityY.y);
    HE_PRINTF("\t\tgradient pressure: ( %f %f )\n", vp_left->m_gradientPressure.x, vp_left->m_gradientPressure.y);

    HE_PRINTF("right state: (%f %f )\n", vp_right->m_location.x, vp_right->m_location.y);
    HE_PRINTF("\tdensity: %f\n", rho_right);
    HE_PRINTF("\tvelocity: ( %f %f )\n", v_right.x, v_right.y);
    HE_PRINTF("\tpressure: %f\n", p_right);
    HE_PRINTF("\t\tgradient density: ( %f %f )\n", vp_right->m_gradientDensity.x, vp_right->m_gradientDensity.y);
    HE_PRINTF("\t\tgradient v_x: ( %f %f )\n", vp_right->m_gradientVelocityX.x, vp_right->m_gradientVelocityX.y);
    HE_PRINTF("\t\tgradient v_y: ( %f %f )\n", vp_right->m_gradientVelocityY.x, vp_right->m_gradientVelocityY.y);
    HE_PRINTF("\t\tgradient pressure: ( %f %f )\n", vp_right->m_gradientPressure.x, vp_right->m_gradientPressure.y);

    // predict the states at the centroid of the face and forward in time half a time step
    double halfTimeStep = sim_core::global_timestep * 0.5;

    Vector l_left = m_center - vp_left->m_cellCenterOfVolume;
    Vector l_right = m_center - vp_right->m_cellCenterOfVolume;

    HE_PRINTF("predicted states in space:\n");

    HE_PRINTF("left state:\n");
    HE_PRINTF("\tdensity: %f\n", rho_left + vp_left->m_gradientDensity.scalarProduct(l_left));
    HE_PRINTF("\tvelocity: ( %f %f )\n", v_left.x + vp_left->m_gradientVelocityX.scalarProduct(l_left),
              v_left.y + vp_left->m_gradientVelocityY.scalarProduct(l_left));
    HE_PRINTF("\tpressure: %f\n", p_left + vp_left->m_gradientPressure.scalarProduct(l_left));

    HE_PRINTF("right state:\n");
    HE_PRINTF("\tdensity: %f\n", rho_right + vp_right->m_gradientDensity.scalarProduct(l_right));
    HE_PRINTF("\tvelocity: ( %f %f )\n", v_right.x + vp_right->m_gradientVelocityX.scalarProduct(l_right),
              v_right.y + vp_right->m_gradientVelocityY.scalarProduct(l_right));
    HE_PRINTF("\tpressure: %f\n", p_right + vp_right->m_gradientPressure.scalarProduct(l_right));

    double rho_left_predicted = rho_left + vp_left->m_gradientDensity.scalarProduct(l_left) +
                                halfTimeStep * (-v_left.scalarProduct(vp_left->m_gradientDensity) -
                                                rho_left * vp_left->div_v());

    Vector v_left_predicted = Vector(0, 0);

    v_left_predicted.x = v_left.x + vp_left->m_gradientVelocityX.scalarProduct(l_left) +
                         halfTimeStep * (-v_left.scalarProduct(vp_left->m_gradientVelocityX) -
                                         (1. / rho_left * vp_left->m_gradientPressure.x));

    v_left_predicted.y = v_left.y + vp_left->m_gradientVelocityY.scalarProduct(l_left) +
                         halfTimeStep * (-v_left.scalarProduct(vp_left->m_gradientVelocityY) -
                                         (1. / rho_left * vp_left->m_gradientPressure.y));

    double p_left_predicted = p_left + vp_left->m_gradientPressure.scalarProduct(l_left) +
                              halfTimeStep * (-v_left.scalarProduct(vp_left->m_gradientPressure) -
                                              Constants::getGAMMA() * p_left * vp_left->div_v());

    double rho_right_predicted = rho_right + vp_right->m_gradientDensity.scalarProduct(l_right) +
                                 halfTimeStep * (-v_right.scalarProduct(vp_right->m_gradientDensity) -
                                                 rho_right * vp_right->div_v());

    Vector v_right_predicted = Vector(0, 0);

    v_right_predicted.x = v_right.x + vp_right->m_gradientVelocityX.scalarProduct(l_right) +
                          halfTimeStep * (-v_right.scalarProduct(vp_right->m_gradientVelocityX) -
                                          (1. / rho_right * vp_right->m_gradientPressure.x));

    v_right_predicted.y = v_right.y + vp_right->m_gradientVelocityY.scalarProduct(l_right) +
                          halfTimeStep * (-v_right.scalarProduct(vp_right->m_gradientVelocityY) -
                                          (1. / rho_right * vp_right->m_gradientPressure.y));

    double p_right_predicted = p_right + vp_right->m_gradientPressure.scalarProduct(l_right) +
                               halfTimeStep * (-v_right.scalarProduct(vp_right->m_gradientPressure) -
                                               Constants::getGAMMA() * p_right * vp_right->div_v());

    if (sqrt(v_left_predicted.length_square()) > 1000 || sqrt(v_right_predicted.length_square()) > 1000)
    {
        fprintf(stderr, "WARNING: predicted velocity is very high in HalfEdge::updateFlux\n");
    }

    HE_PRINTF("predicted states in space and time:\n");

    HE_PRINTF("left state:\n");
    HE_PRINTF("\tdensity: %f\n", rho_left_predicted);
    HE_PRINTF("\tvelocity: ( %f %f )\n", v_left_predicted.x, v_left_predicted.y);
    HE_PRINTF("\tpressure: %f\n", p_left_predicted);

    HE_PRINTF("right state:\n");
    HE_PRINTF("\tdensity: %f\n", rho_right_predicted);
    HE_PRINTF("\tvelocity: ( %f %f )\n", v_right_predicted.x, v_right_predicted.y);
    HE_PRINTF("\tpressure: %f\n", p_right_predicted);

    // rotate the states such that the edge is parallel to the y-axis

    double phi = 0;

    if (m_length != 0) // only rotate if the length of the half edge is greater than zero
    {
        // the angle between the half edge and the x-axis
        phi = (vp_right->m_location - vp_left->m_location).phi();

        HE_PRINTF("angle of the half edge: %f\n", phi);

        v_left_predicted.rotate(-phi);
        v_right_predicted.rotate(-phi);

        HE_PRINTF("rotated states:\n");
        HE_PRINTF("\tleft state velocity: ( %f %f )\n", v_left_predicted.x, v_left_predicted.y);
        HE_PRINTF("\tright state velocity: ( %f %f )\n", v_right_predicted.x, v_right_predicted.y);
    }

    // solve the Riemann problem
    double rho_middle = 0;
    Vector u_middle = Vector(0, 0);
    double p_middle = 0;

    exact_riemann_solver::solveAndSampleZeroSpeed(&rho_middle, &u_middle.x, &p_middle, rho_left_predicted,
                                                  v_left_predicted.x, p_left_predicted, rho_right_predicted,
                                                  v_right_predicted.x, p_right_predicted);

    HE_PRINTF("state on the face:\n");
    HE_PRINTF("\tdensity: %f\n", rho_middle);
    HE_PRINTF("\tvelocity: ( %f %f )\n", u_middle.x, u_middle.y);
    HE_PRINTF("\tpressure: %f\n", p_middle);

    v_left.rotate(-phi);
    v_right.rotate(-phi);

    if (v_left.x + v_right.x > 0)
    {
        u_middle.y = v_left_predicted.y;
    }
    else
    {
        u_middle.y = v_right_predicted.y;
    }

    // retransform the state to the lab frame

    u_middle.rotate(phi);
    u_middle += w;

    HE_PRINTF("lab state:\n");
    HE_PRINTF("\tdensity: %f\n", rho_middle);
    HE_PRINTF("\tvelocity: ( %f %f )\n", u_middle.x, u_middle.y);
    HE_PRINTF("\tpressure: %f\n", p_middle);

    // calculate the net flux
    Vector r = u_middle - w; //(v_lab - w)
    Vector n = vp_right->m_location - vp_left->m_location; // the normal vector of the half edge
    n.normalize();

    double mass_net_flux = r.scalarProduct(n) * rho_middle;
    double momentum_net_flux_x = rho_middle * u_middle.x * r.scalarProduct(n) + p_middle * n.x;
    double momentum_net_flux_y = rho_middle * u_middle.y * r.scalarProduct(n) + p_middle * n.y;

    assert(rho_middle != 0);
    double e_lab = 0.5 * u_middle.length_square() + p_middle / ((Constants::getGAMMA() - 1) * rho_middle);

    double energy_net_flux = rho_middle * e_lab * r.scalarProduct(n) + p_middle * u_middle.scalarProduct(n);

    HE_PRINTF("net fluxes:\n");
    HE_PRINTF("\tmass: %f\n", mass_net_flux);
    HE_PRINTF("\tmomentum_x: %f\n", momentum_net_flux_x);
    HE_PRINTF("\tmomentum_y: %f\n", momentum_net_flux_y);
    HE_PRINTF("\tenergy: %f\n", energy_net_flux);

    m_massFlux = mass_net_flux;
    m_momentumFlux.x = momentum_net_flux_x;
    m_momentumFlux.y = momentum_net_flux_y;
    m_energyFlux = energy_net_flux;
    m_fluxUpdated = true;

    m_twin->m_massFlux = -mass_net_flux;
    m_twin->m_momentumFlux.x = -momentum_net_flux_x;
    m_twin->m_momentumFlux.y = -momentum_net_flux_y;
    m_twin->m_energyFlux = -energy_net_flux;
    m_twin->m_fluxUpdated = true;
}

Vertex::Vertex(const Point* position, HalfEdge* IncidentEdge)
{
    m_inUse = false;
    m_position = *position;
    m_edge = IncidentEdge;
}

Vertex::Vertex(const Point* position)
{
    m_inUse = false;
    m_position = *position;
    m_edge = NULL;
}

Vertex::Vertex()
{
    m_inUse = false;
    m_position = Point(0, 0);
    m_edge = NULL;
}

Vertex::~Vertex()
{

}

// calculate the velocity of the face (Springel, p. 805)
Vector HalfEdge::velocity()
{
    HE_PRINTF("\n\n\nhalf edge center: ( %f %f )\n", m_center.x, m_center.y);

    assert(m_vp->isGhost() ==
           false); // don't calculate the velocity of a ghost half edge, because m_center has not been calculated

    Point r_L = m_vp->m_location;
    Point r_R = m_twin->m_vp->m_location;

    Vector w_L = m_vp->m_particleVelocity;
    Vector w_R = m_twin->m_vp->m_particleVelocity;

    Vector primary = (w_R + w_L) * 0.5;
    Vector secondary;

    if ((r_R - r_L).length_square() < pow(10., -6))
    {
        secondary = Vector(0, 0);
    }
    else
    {
        secondary =
                (r_R - r_L) * (w_L - w_R).scalarProduct(m_center - ((r_R + r_L) * 0.5)) / (r_R - r_L).length_square();
    }

    if (sqrt((primary + secondary).length_square()) > 50)
    {
        fprintf(stderr, "WARNING: half edge velocity is very high\n");
    }

    return (primary + secondary);
}

// does the half edge have zero length?
bool HalfEdge::is_zero()
{
    return (m_length < pow(10., -12.));
}

// write data into in a file
void HalfEdge::info(FILE* pFile)
{
    fprintf(pFile, "address: %p\n", this);
    fprintf(pFile, "origin: %p\n", m_origin);

    if (m_origin != NULL)
    {
        fprintf(pFile, "\tx=%f y=%f\n", m_origin->m_position.x, m_origin->m_position.y);
    }

    if (m_next != NULL)
    {
        if (m_next->m_origin != NULL)
        {
            fprintf(pFile, "end: x=%f y=%f\n", m_next->m_origin->m_position.x, m_next->m_origin->m_position.y);
        }
    }

    if (m_vp != NULL)
    {
        fprintf(pFile, "voronoi particle: x=%f y=%f\n", m_vp->m_location.x, m_vp->m_location.y);
    }

    fprintf(pFile, "twin: %p\n", m_twin);
    fprintf(pFile, "next: %p\n", m_next);
    fprintf(pFile, "prev: %p\n", m_prev);
    fprintf(pFile, "length: %f\n", m_length);
}

// write data into a file
void Vertex::info(FILE* pFile)
{
    fprintf(pFile, "address: %p\n", this);
    fprintf(pFile, "position: x=%f, y=%f\n", m_position.x, m_position.y);
    fprintf(pFile, "incident edge: %p\n", m_edge);
}

// check whether the half face starting at start is closed
bool graph::isClosed(HalfEdge* start)
{
    HalfEdge* startEdge = start;
    HalfEdge* currentEdge = start;
    HalfEdge* next = NULL;

    do
    {
        next = currentEdge->m_next;
        if (next == NULL)
        { return false; }
        currentEdge = next;

    } while (currentEdge != startEdge);

    return true;
}

// print a face in gnuplot readable format
void graph::printFace(FILE* pFile, HalfEdge* start)
{
    if (!isClosed(start))
    {
        fprintf(stderr, "ERROR in graph::printFace: face not closed around the voronoi particle: x=%f, y=%f\n",
                start->m_vp->m_location.x, start->m_vp->m_location.y);
        return;

    } // face is not closed, return

    HalfEdge* startEdge = start;
    HalfEdge* currentEdge = start;
    HalfEdge* next = NULL;

    do
    {
        if (currentEdge->m_origin == NULL)
        {
            fprintf(stderr, "ERROR in graph::printFace: half edge has no origin\n");

            return;
        }
        fprintf(pFile, "%f\t%f\n", currentEdge->m_origin->m_position.x, currentEdge->m_origin->m_position.y);

        next = currentEdge->m_next;
        currentEdge = next;

    } while (currentEdge != startEdge);

    fprintf(pFile, "%f\t%f\n\n", start->m_origin->m_position.x,
            start->m_origin->m_position.y); //print first point twice for Gnuplot (with linespoints)
}

// print a half edge arrow
void graph::printHalfEdge(FILE* pFile, const HalfEdge* halfEdge)
{
    double shift = 0.0015;
    double arrowlength = shift / 2;

    if (halfEdge->m_origin == NULL)
    {
        fprintf(stderr, "ERROR in printHalfEdge: half edge has no origin\n");
        return;
    }

    if (halfEdge->m_vp == NULL)
    {
        fprintf(stderr, "ERROR in printHalfEdge: half edge has no voronoi particle\n");
        // this should not happen at all!
        return;
    }

    // vector from vertex to center of volume
    Point v(halfEdge->m_vp->m_cellCenterOfVolume.x - halfEdge->m_origin->m_position.x,
            halfEdge->m_vp->m_cellCenterOfVolume.y - halfEdge->m_origin->m_position.y);

    if (!(v.x == 0 && v.y == 0))
    {
        v.normalize();
    }

    v.x *= shift;
    v.y *= shift;

    Point origin = halfEdge->m_origin->m_position;
    origin = origin + v;

    fprintf(pFile, "%f\t%f\n", origin.x, origin.y);

    if (halfEdge->m_prev->m_origin == NULL)
    {
        fprintf(stderr, "ERROR in printHalfEdge: previous half edge has no origin\n");
        return;
    }

    // vector from vertex to vertex
    Point w(halfEdge->m_origin->m_position.x - halfEdge->m_prev->m_origin->m_position.x,
            halfEdge->m_origin->m_position.y - halfEdge->m_prev->m_origin->m_position.y);

    if (!(w.x == 0 && w.y == 0))
    {
        w.normalize();
    }

    // print arrowhead
    Point ortho = w.orthogonal();

    ortho.x *= arrowlength;
    ortho.y *= arrowlength;

    w.x *= arrowlength * 2;
    w.y *= arrowlength * 2;

    Point head = (origin + ortho) - w;

    fprintf(pFile, "%f\t%f\n", head.x, head.y);
    fprintf(pFile, "%f\t%f\n", origin.x, origin.y);

    head = (origin - ortho) - w;
    fprintf(pFile, "%f\t%f\n", head.x, head.y);
    fprintf(pFile, "%f\t%f\n", origin.x, origin.y);
}

// print a half Face
void graph::printHalfFace(FILE* pFile, HalfEdge* start)
{
    if (!isClosed(start))
    { return; } // the half face is not closed, return

    HalfEdge* startEdge = start;
    HalfEdge* currentEdge = start;
    HalfEdge* next = NULL;

    do
    {
        printHalfEdge(pFile, currentEdge);

        next = currentEdge->m_next;
        currentEdge = next;

    } while (currentEdge != startEdge);

    // print fist point twice for Gnuplot (with linespoints)
    printHalfEdge(pFile, start);
    fprintf(pFile, "\n");
}