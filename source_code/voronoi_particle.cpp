#include "config.h"
#include "voronoi_particle.h"
#include "assert.h"
#include "constants.h"
#include "graph.h"
#include "debug.h"
#include <cmath>
#include "sim_core.h"
#include <limits>
#include "sim_graph.h"

VoronoiParticle::VoronoiParticle(const Point* location)
        :
        m_location(*location),
        m_particleVelocity(0, 0),
        m_outerComponent(NULL),
        m_cellCenterOfVolume(0, 0),
        m_velocity(0, 0),
        m_density(0),
        m_pressure(0),
        m_mass(0),
        m_momentum(0, 0),
        m_energy(0),
        m_gradientVelocityX(0, 0),
        m_gradientVelocityY(0, 0),
        m_gradientDensity(0, 0),
        m_gradientPressure(0, 0),
        m_amr_twin_particle(NULL),
        m_color(0),
        m_flag_1(0),
        m_obstacle_particle(false),
        m_inner_obstacle_particle(false),
        m_copiedOver(NULL)
{

}

VoronoiParticle::VoronoiParticle()
        :
        m_location(0, 0),
        m_particleVelocity(0, 0),
        m_outerComponent(NULL),
        m_cellCenterOfVolume(0, 0),
        m_velocity(0, 0),
        m_density(0),
        m_pressure(0),
        m_mass(0),
        m_momentum(0, 0),
        m_energy(0),
        m_gradientVelocityX(0, 0),
        m_gradientVelocityY(0, 0),
        m_gradientDensity(0, 0),
        m_gradientPressure(0, 0),
        m_amr_twin_particle(NULL),
        m_color(0),
        m_flag_1(0),
        m_obstacle_particle(false),
        m_inner_obstacle_particle(false),
        m_copiedOver(NULL)
{

}

VoronoiParticle::VoronoiParticle(double x, double y)
        :
        m_location(x, y),
        m_particleVelocity(0, 0),
        m_outerComponent(NULL),
        m_cellCenterOfVolume(0, 0),
        m_velocity(0, 0),
        m_density(0),
        m_pressure(0),
        m_mass(0),
        m_momentum(0, 0),
        m_energy(0),
        m_gradientVelocityX(0, 0),
        m_gradientVelocityY(0, 0),
        m_gradientDensity(0, 0),
        m_gradientPressure(0, 0),
        m_amr_twin_particle(NULL),
        m_color(0),
        m_flag_1(0),
        m_obstacle_particle(false),
        m_inner_obstacle_particle(false),
        m_copiedOver(NULL)
{

}

// get a primitive variable
double& VoronoiParticle::var(prim_var pv)
{
    if (pv == DENSITY)
    {
        return m_density;
    }

    if (pv == PRESSURE)
    {
        return m_pressure;
    }

    if (pv == VELOCITY_X)
    {
        return m_velocity.x;
    }

    if (pv == VELOCITY_Y)
    {
        return m_velocity.y;
    }

    fprintf(stderr, "ERROR in cell::var: invalid primitive variable\n");
    exit(1);
}

// get a gradient
Vector& VoronoiParticle::grad(prim_var pv)
{
    if (pv == DENSITY)
    {
        return m_gradientDensity;
    }

    if (pv == PRESSURE)
    {
        return m_gradientPressure;
    }

    if (pv == VELOCITY_X)
    {
        return m_gradientVelocityX;
    }

    if (pv == VELOCITY_Y)
    {
        return m_gradientVelocityY;
    }

    fprintf(stderr, "ERROR in cell::grad: invalid primitive variable\n");
    exit(1);
}

// calculate the local time step of the cell
double VoronoiParticle::local_timestep()
{
    if (m_inner_obstacle_particle)
    {
        return std::numeric_limits<double>::max();
    }

    double denominator = soundSpeed() + sqrt((m_velocity - m_particleVelocity).length_square());

    if (denominator == 0)
    {
        return std::numeric_limits<double>::max();
    }
    else
    {
        return (Constants::getCFL() * radius() / denominator);
    }
}

// is the particle a ghost particle?
bool VoronoiParticle::isGhost()
{
    return false;
}

// write stored data into a file
void VoronoiParticle::info(FILE* pFile)
{
    fprintf(pFile, "address: %p\n", this);
    fprintf(pFile, "location: x=%f, y=%f\n", m_location.x, m_location.y);
    fprintf(pFile, "outer half edge: %p\n", m_outerComponent);
}

// write vertex positions into a file (debug)
void VoronoiParticle::info_vertices(FILE* pFile)
{
    HalfEdge* start = m_outerComponent;
    HalfEdge* current = m_outerComponent;

    do
    {
        fprintf(pFile, "vertex: (%f, %f)\n", current->m_origin->m_position.x, current->m_origin->m_position.y);

        current = current->m_next;

    } while (current != start);
}

// calculate the center of volume of a triangle
static Point triangleCenterOfVolume(Point* a, Point* b, Point* c)
{
    Point center;
    center.x = (a->x + b->x + c->x) / 3;
    center.y = (a->y + b->y + c->y) / 3;

    return center;
}

// calculate the volume of a triangle
static double triangleArea(Point* a, Point* b, Point* c)
{
    return (0.5 * (a->x * (b->y - c->y) + b->x * (c->y - a->y) + c->x * (a->y - b->y)));
}

// calculate and set the volume and the center of volume of the corresponding cell
void VoronoiParticle::updateCellVolumeAndCenter()
{
    VP_PRINTF("Voronoi particle: ( %f , %f )\n", m_location.x, m_location.y);

    assert(m_outerComponent != NULL);
    assert(graph::isClosed(m_outerComponent));

    // let's walk around the cell
    HalfEdge* start = m_outerComponent;
    HalfEdge* currentHalfEdge = start;

    // startpoint and endpoint of the half edge
    Point endPoint = start->m_next->m_origin->m_position;
    Point startPoint = start->m_origin->m_position;

    // entire area of the cell
    double entireArea = 0;

    // area of a triangle
    double triArea = 0;

    Point cellCenter(0, 0);

    do
    {
        triArea = triangleArea(&m_location, &startPoint, &endPoint);

        entireArea += triArea;
        cellCenter += (triangleCenterOfVolume(&m_location, &startPoint, &endPoint) * triArea);

        currentHalfEdge = currentHalfEdge->m_next;
        startPoint = endPoint;
        endPoint = currentHalfEdge->m_next->m_origin->m_position;


    } while (currentHalfEdge != start);

    // set the member variables
    m_cellVolume = entireArea;

    if (entireArea != 0)
    {
        m_cellCenterOfVolume = cellCenter / entireArea;

    }
    else
    {
        m_cellCenterOfVolume = m_location;
    }

    VP_PRINTF("\tcell center of volume: ( %f , %f )\n", m_cellCenterOfVolume.x, m_cellCenterOfVolume.y);
    VP_PRINTF("\tcell volume: %f\n", m_cellVolume);
}

// calculate the gradient
void VoronoiParticle::calculate_gradient(prim_var pv)
{
    VP_PRINTF("particle: ( %f , %f )\n", m_location.x, m_location.y);
    assert(graph::isClosed(m_outerComponent));

    grad(pv) = Vector(0, 0);
    Vector c_ij = Vector(0, 0);
    Vector r_ij = Vector(0, 0);
    double r_ij_length = 0;

    VoronoiParticle* vp_j = NULL;
    HalfEdge* start = m_outerComponent;
    HalfEdge* current = m_outerComponent;

    do
    {
        vp_j = current->m_twin->m_vp;
        c_ij = current->m_center - ((m_location + vp_j->m_location) * 0.5);
        r_ij = m_location - vp_j->m_location;
        r_ij_length = sqrt(r_ij.length_square());

        grad(pv) += ((c_ij / r_ij_length * (vp_j->var(pv) - var(pv))) -
                     (r_ij / r_ij_length * (var(pv) + vp_j->var(pv)) * 0.5)) * current->m_length;

        current = current->m_next;

    } while (current != start);

    grad(pv) /= m_cellVolume;
}

// helper functions for the slope limiter
//

static double neighbour_maximum(VoronoiParticle* particle, prim_var pv)
{
    double maximum = particle->var(pv);

    HalfEdge* start = particle->m_outerComponent;
    HalfEdge* current = particle->m_outerComponent;

    do
    {
        if (current->m_twin->m_vp->var(pv) > maximum)
        {
            maximum = current->m_twin->m_vp->var(pv);
        }

        current = current->m_next;

    } while (current != start);

    return maximum;
}

static double neighbour_minimum(VoronoiParticle* particle, prim_var pv)
{
    double minimum = particle->var(pv);

    HalfEdge* start = particle->m_outerComponent;
    HalfEdge* current = particle->m_outerComponent;

    do
    {
        if (current->m_twin->m_vp->var(pv) < minimum)
        {
            minimum = current->m_twin->m_vp->var(pv);
        }

        current = current->m_next;

    } while (current != start);

    return minimum;
}

// slope limiter
void VoronoiParticle::slope_limit_gradient_arepo(prim_var pv)
{
    double alpha = 1;
    double psi_j = 0;
    double phi_j = 0;

    double maximum = neighbour_maximum(this, pv);
    double minimum = neighbour_minimum(this, pv);

    HalfEdge* start = m_outerComponent;
    HalfEdge* current = m_outerComponent;

    do
    {
        phi_j = this->grad(pv).scalarProduct(current->m_center - m_cellCenterOfVolume);

        if (phi_j > 0)
        {
            psi_j = (maximum - var(pv)) / phi_j;
        }

        else if (phi_j < 0)
        {
            psi_j = (minimum - var(pv)) / phi_j;
        }

        else
        {
            psi_j = 1;
        }

        if (psi_j < alpha)
        {
            alpha = psi_j;
        }

        current = current->m_next;

    } while (current != start);

    grad(pv) *= alpha;
}

static double maximum(double a, double b)
{
    if (a > b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

// slope limiter
void VoronoiParticle::slope_limit_gradient_tess(prim_var pv)
{
    double theta = Constants::getSL_THETA();
    double alpha = 1;
    double psi_j = 0;
    double phi_j = 0;

    HalfEdge* start = m_outerComponent;
    HalfEdge* current = m_outerComponent;

    do
    {
        phi_j = this->grad(pv).scalarProduct(current->m_center - m_cellCenterOfVolume);

        if (phi_j > 0)
        {
            psi_j = maximum(theta * (current->m_twin->m_vp->var(pv) - var(pv)) / phi_j, 0);
        }
        else if (phi_j < 0)
        {
            psi_j = maximum(theta * (current->m_twin->m_vp->var(pv) - var(pv)) / phi_j, 0);
        }
        else
        {
            psi_j = 1;
        }

        if (psi_j < alpha)
        {
            alpha = psi_j;
        }

        current = current->m_next;

    } while (current != start);

    grad(pv) *= alpha;
}

// update the conserved variables with the fluxes
void VoronoiParticle::updateConservedVariables()
{
    if (m_inner_obstacle_particle)
    {
        m_mass = 1;
        m_momentum = Vector(0, 0);
        m_energy = 1;

        return;
    }

    VP_PRINTF("particle: ( %f , %f )\n", m_location.x, m_location.y);

    VP_PRINTF("\n");
    VP_PRINTF("old mass: %f\n", m_mass);
    VP_PRINTF("old momentum: ( %f %f )\n", m_momentum.x, m_momentum.y);
    VP_PRINTF("old energy: %f\n", m_energy);
    VP_PRINTF("\n");

    assert(graph::isClosed(m_outerComponent));

    HalfEdge* start = m_outerComponent;
    HalfEdge* current = m_outerComponent;

    do
    {
        m_mass -= sim_core::global_timestep * current->m_length * current->m_massFlux;
        m_momentum -= (current->m_momentumFlux * sim_core::global_timestep * current->m_length);
        m_energy -= sim_core::global_timestep * current->m_length * current->m_energyFlux;

        current = current->m_next;

    } while (current != start);

    // external acceleration
    if (Constants::getENABLE_ACCELERATION())
    {
        Vector acceleration = Vector(Constants::getACCELERATION_X(), Constants::getACCELERATION_Y());

        m_momentum += acceleration * sim_core::global_timestep * m_cellVolume * m_density;
        m_energy += acceleration.scalarProduct(m_velocity) * m_density * sim_core::global_timestep * m_cellVolume;
    }

    if (Constants::getMASS_DENSITY_FLOOR() != -1)
    {
        if (m_mass / m_cellVolume < Constants::getMASS_DENSITY_FLOOR())
        {
            m_mass = Constants::getMASS_DENSITY_FLOOR() * m_cellVolume;
        }
    }

    if (Constants::getENERGY_DENSITY_FLOOR() != -1)
    {
        if (m_energy / m_cellVolume < Constants::getENERGY_DENSITY_FLOOR())
        {
            m_energy = Constants::getENERGY_DENSITY_FLOOR() * m_cellVolume;
        }

    }

    if (m_mass <= 0 || m_energy <= 0)
    {
        fprintf(stderr, "ERROR in VoronoiParticle::updateConservedVariables: invalid variables: mass: %f, energy: %f\n",
                m_mass, m_energy);
        fprintf(stderr, "location: ( %f , %f )\n", m_location.x, m_location.y);
        if (m_inner_obstacle_particle)
        {
            fprintf(stderr, "(inner obstacle particle)\n");
        }
        exit(1);
    }

    VP_PRINTF("new mass: %f\n", m_mass);
    VP_PRINTF("new momentum: ( %f %f )\n", m_momentum.x, m_momentum.y);
    VP_PRINTF("new energy: %f\n", m_energy);
}

// update the velocities of the particles in order to achieve mesh regularity
void VoronoiParticle::mesh_regulation_arepo_1()
{
    // distance between the center of mass and the Voronoi particle
    double d = sqrt((m_cellCenterOfVolume - m_location).length_square());

    // cell radius
    double cellRadius = radius();
    assert(cellRadius != 0);

    // mesh regularity parameters
    double eta = Constants::getMESH_ETA();
    double chi = Constants::getMESH_CHI();

    Vector v = Vector(0, 0);
    double a = d / (eta * cellRadius);

    if (a < 0.9)
    {
        // v stays zero
    }
    else if (a >= 1.1)
    {
        v = (m_cellCenterOfVolume - m_location) * soundSpeed() / d;
    }
    else // a is between 0.9 and 1.1
    {
        v = (m_cellCenterOfVolume - m_location) * soundSpeed() / d * (d - 0.9 * eta * cellRadius) /
            (0.2 * eta * cellRadius);
    }

    m_particleVelocity += v * chi;
}

void VoronoiParticle::mesh_regulation_arepo_2()
{
    fprintf(stderr, "ERROR in VoronoiParticle::mesh_regulation_arepo_2: this mesh regulation is not recommended\n");
    exit(1);

    // distance between the center of mass and the Voronoi particle
    double d = sqrt((m_cellCenterOfVolume - m_location).length_square());

    // cell radius
    double cellRadius = radius();
    assert(cellRadius != 0);

    // mesh regularity parameters
    double eta = Constants::getMESH_ETA();
    double chi = Constants::getMESH_CHI();

    Vector v = Vector(0, 0);
    double a = d / (eta * cellRadius);

    if (a < 0.9)
    {
        // v stays zero
    }
    else if (a >= 1.1)
    {
        v = (m_cellCenterOfVolume - m_location) / sim_core::global_timestep;
    }
    else // a is between 0.9 and 1.1
    {
        v = (m_cellCenterOfVolume - m_location) / sim_core::global_timestep * (d - 0.9 * eta * cellRadius) /
            (0.2 * eta * cellRadius);
    }

    if (sqrt(v.length_square()) > 1000)
    {
        fprintf(stderr, "WARNING: particle velocity is very high because of mesh_regulation_arepo2\n");
    }

    m_particleVelocity += v * chi;
}

void VoronoiParticle::roundness_correction()
{
    double rs = roundness();

    if (rs > Constants::getMESH_OMEGA())
    {
        m_particleVelocity -= m_velocity * rs * Constants::getMESH_PSI();
    }
}

// calculate the primitive variables of the cell
void VoronoiParticle::calculatePrimitiveVariables()
{
    if (m_inner_obstacle_particle)
    {
        m_density = 1;
        m_pressure = 1;
        m_velocity = Vector(0, 0);
        return;
    }

    VP_PRINTF("particle: ( %f , %f )\n", m_location.x, m_location.y);

    assert(m_cellVolume != 0);
    assert(m_mass != 0);

    m_density = m_mass / m_cellVolume;

    m_velocity = m_momentum / m_mass;

    m_pressure = (m_energy / m_cellVolume - 0.5 * m_density * m_velocity.length_square()) *
                 (Constants::getGAMMA() - 1.);

    if (Constants::getPRESSURE_FLOOR() != -1)
    {
        if (m_pressure < Constants::getPRESSURE_FLOOR())
        {
            m_pressure = Constants::getPRESSURE_FLOOR();
        }
    }

    if (m_density <= 0 || m_pressure < 0)
    {
        fprintf(stderr,
                "ERROR in VoronoiParticle::calculatePrimitiveVariables: invalid variables: density: %f, pressure: %f\n",
                m_density, m_pressure);
        fprintf(stderr, "location: ( %f , %f )\n", m_location.x, m_location.y);

        if (m_inner_obstacle_particle)
        {
            fprintf(stderr, "(inner obstacle particle)\n");
        }
        exit(1);
    }

    VP_PRINTF("\t density: %f\n", m_density);
    VP_PRINTF("\t velocity: ( %f , %f )\n", m_velocity.x, m_velocity.y);
    VP_PRINTF("\t pressure: %f\n", m_pressure);
}

// calculate the conserved variables of the cell
void VoronoiParticle::calculateConservedVariables()
{
    VP_PRINTF("Voronoi particle: ( %f , %f )\n", m_location.x, m_location.y);

    m_mass = m_density * m_cellVolume;
    m_momentum = m_velocity * m_density * m_cellVolume;
    m_energy = (0.5 * m_density * m_velocity.length_square() + m_pressure / (Constants::getGAMMA() - 1.)) *
               m_cellVolume;

    if (m_mass < 0 || m_energy < 0)
    {
        fprintf(stderr,
                "ERROR in VoronoiParticle::calculateConservedVariables(): invalid variables: mass: %f, energy: %f\n",
                m_mass, m_energy);
        exit(1);
    }

    VP_PRINTF("\t mass: %f\n", m_mass);
    VP_PRINTF("\t momentum: ( %f , %f )\n", m_momentum.x, m_momentum.y);
    VP_PRINTF("\t energy: %f\n", m_energy);
}

// delete the half edges which have zero length
void VoronoiParticle::delete_zero_length_edges()
{
    HalfEdge* start = m_outerComponent;
    HalfEdge* currentEdge = m_outerComponent;
    HalfEdge* next = NULL;

    do
    {
        next = currentEdge->m_next;

        if (currentEdge->is_zero())
        {
            assert(currentEdge != start);
            sim_graph::delete_edge(currentEdge);
        }

        currentEdge = next;

    } while (currentEdge != start);
}

// move the particle according to its velocity
void VoronoiParticle::move()
{
    m_location += (m_particleVelocity * sim_core::global_timestep);
}

void VoronoiParticle::move_backwards()
{
    m_location -= (m_particleVelocity * sim_core::global_timestep);
}

int VoronoiParticle::number_of_neighbours()
{
    HalfEdge* start = m_outerComponent;
    HalfEdge* current = m_outerComponent;

    int counter = 0;

    do
    {
        current = current->m_next;
        counter++;

    } while (start != current);

    return counter;
}

// calculate the roundness of the cell
double VoronoiParticle::roundness()
{
    HalfEdge* start = m_outerComponent;
    HalfEdge* current = m_outerComponent;

    double max_distance = 0;

    do
    {
        if (sqrt((current->m_origin->m_position - m_cellCenterOfVolume).length_square()) > max_distance)
        {
            max_distance = sqrt((current->m_origin->m_position - m_cellCenterOfVolume).length_square());
        }

        current = current->m_next;

    } while (start != current);

    return ((max_distance - radius()) / radius());
}

// calculate the vorticity of the cell
double VoronoiParticle::vorticity()
{
    return (m_gradientVelocityY.x - m_gradientVelocityX.y);
}