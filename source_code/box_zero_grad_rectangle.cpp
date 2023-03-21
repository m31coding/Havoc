#include "config.h"
#include "debug.h"
#include "box_zero_grad_rectangle.h"
#include <cstdlib>
#include "sim.h"
#include <cmath>
#include "point_location.h"

/// flags for zero gradient boundaries ghost particles
enum e_zgb
{
    INNER = 0,
    FIXED = 1,
    FIXED_TOP = 2,
    FIXED_BOTTOM = 3,
    FIXED_LEFT = 4,
    FIXED_RIGHT = 5,
    VP_MOVED_OUT = 10,
    GHOST_MOVED_IN = 11,
    GHOST_MOVED_OUT = 12,
    VERY_INNER = 13,
    INNER_UPDATED = 14,
    INNER_IN_PROGRESS = 15
};

#define IS_FIXED(n) (n <= 5 && n >= 1)
#define IS_FIXED_INNER(n) (n <= 5 && n >= 2)

const unsigned int LAYERS = 8;

// constructor
ZeroGradRectangle::ZeroGradRectangle(Sim* sim, double Xmin, double Xmax, double Ymin, double Ymax)
        :
        m_Xmin(Xmin),
        m_Xmax(Xmax),
        m_Ymin(Ymin),
        m_Ymax(Ymax),
        m_deltaX(Xmax - Xmin),
        m_deltaY(Ymax - Ymin),
        m_topRight(Xmax, Ymax),
        m_topLeft(Xmin, Ymax),
        m_bottomLeft(Xmin, Ymin),
        m_bottomRight(Xmax, Ymin),
        m_middle(Xmin + (Xmax - Xmin) * 0.5, Ymin + (Ymax - Ymin) * 0.5),
        m_ghost_box_x_min(0),
        m_ghost_box_x_max(0),
        m_ghost_box_y_min(0),
        m_ghost_box_y_max(0),
        m_ghost_points(NULL),
        m_nof_ghost_points(0)
{
    m_sim = sim;

    BOX_PRINTF("\nbottom right: ( %f , %f )\n", m_bottomRight.x, m_bottomRight.y);
    BOX_PRINTF("top right: ( %f , %f )\n", m_topRight.x, m_topRight.y);
    BOX_PRINTF("top left: ( %f , %f )\n", m_topLeft.x, m_topLeft.y);
    BOX_PRINTF("bottom left: ( %f , %f )\n", m_bottomLeft.x, m_bottomLeft.y);
    BOX_PRINTF("delta x: %f\ndelta y: %f\n", m_deltaX, m_deltaY);

    // number of ghost points = number of maximum ghost particles (s.a. sim.cpp)
    m_ghost_points = new zgr_ghost[(int) (factorGhostParticles * sim->getMaxNumberOfParticles() + 10000)];

    // horizontal distance of the particles
    m_dh = getDeltaX() / Constants::getNOF_PARTICLES_X();

    // vertical distance of the particles
    m_dv = getDeltaY() / Constants::getNOF_PARTICLES_Y();

    m_DV = m_dh * m_dv;

    m_ghost_box_x_max = m_Xmax + m_dh * (LAYERS - 2);
    m_ghost_box_x_min = m_Xmin - m_dh * (LAYERS - 2);
    m_ghost_box_y_max = m_Ymax + m_dv * (LAYERS - 2);
    m_ghost_box_y_min = m_Ymin - m_dv * (LAYERS - 2);

    BOX_PRINTF("m_ghost_box_x_max = %f\n", m_ghost_box_x_max);
    BOX_PRINTF("m_ghost_box_x_min = %f\n", m_ghost_box_x_min);
    BOX_PRINTF("m_ghost_box_y_max = %f\n", m_ghost_box_y_max);
    BOX_PRINTF("m_ghost_box_y_min = %f\n", m_ghost_box_y_min);
}

// destructor
ZeroGradRectangle::~ZeroGradRectangle()
{
    delete[] m_ghost_points;
}

// reset function
void ZeroGradRectangle::reset()
{
    BOX_PRINTF("\n--- reset ---\n\n");

    m_nof_ghost_points = 0;
}

// check whether a point is inside the box
bool ZeroGradRectangle::isInside(Point* point)
{
    return (point->x <= m_Xmax && point->x >= m_Xmin && point->y <= m_Ymax && point->y >= m_Ymin);
}

// check whether a point is inside the ghost box
bool ZeroGradRectangle::is_inside_ghost_box(Point* point)
{
    return (point->x <= m_ghost_box_x_max && point->x >= m_ghost_box_x_min && point->y <= m_ghost_box_y_max &&
            point->y >= m_ghost_box_y_min);
}

// check the fixed ghosts
bool ZeroGradRectangle::is_fixed_left(Point* point)
{
    return (point->x <= m_ghost_box_x_min && point->x >= m_ghost_box_x_min - m_dh &&
            point->y <= m_ghost_box_y_max - m_dv && point->y >= m_ghost_box_y_min);
}

bool ZeroGradRectangle::is_fixed_right(Point* point)
{
    return (point->x <= m_ghost_box_x_max + m_dh && point->x >= m_ghost_box_x_max && point->y <= m_ghost_box_y_max &&
            point->y >= m_ghost_box_y_min + m_dv);
}

bool ZeroGradRectangle::is_fixed_top(Point* point)
{
    return (point->x >= m_ghost_box_x_min && point->x <= m_ghost_box_x_max - m_dh && point->y >= m_ghost_box_y_max &&
            point->y <= m_ghost_box_y_max + m_dv);
}

bool ZeroGradRectangle::is_fixed_bottom(Point* point)
{
    return (point->x <= m_ghost_box_x_max && point->x >= m_ghost_box_x_min + m_dh && point->y <= m_ghost_box_y_min &&
            point->y >= m_ghost_box_y_min - m_dv);
}

// check the volume of a fixed cell and insert a new particle if needed
void ZeroGradRectangle::potential_particle_insertion(GhostParticle* ghost)
{
    if (IS_FIXED_INNER(ghost->m_flag_1))
    {
        ghost->updateCellVolumeAndCenter();

        if (ghost->m_cellVolume >= 1.4999999 * m_dh * m_dv)
        {
            if (ghost->m_flag_1 == FIXED_LEFT)
            {
                m_ghost_points[m_nof_ghost_points].location = Point(ghost->m_location.x + m_dh, ghost->m_location.y);
            }
            else if (ghost->m_flag_1 == FIXED_RIGHT)
            {
                m_ghost_points[m_nof_ghost_points].location = Point(ghost->m_location.x - m_dh, ghost->m_location.y);
            }
            else if (ghost->m_flag_1 == FIXED_TOP)
            {
                m_ghost_points[m_nof_ghost_points].location = Point(ghost->m_location.x, ghost->m_location.y - m_dv);
            }
            else if (ghost->m_flag_1 == FIXED_BOTTOM)
            {
                m_ghost_points[m_nof_ghost_points].location = Point(ghost->m_location.x, ghost->m_location.y + m_dv);
            }
            else
            {
                assert(false);
            }

            m_ghost_points[m_nof_ghost_points].flag = INNER;
            m_nof_ghost_points++;
        }
    }
}

// determine ghost particles
void ZeroGradRectangle::determineGhostParticles()
{
    BOX_PRINTF("\n--- determining ghost particles ---\n\n");

    GhostParticle* ghost = NULL;

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        ghost = &m_sim->m_ghosts[i];

        if (IS_FIXED(ghost->m_flag_1))
        {
            // keep the ghost
            m_ghost_points[m_nof_ghost_points].location = ghost->m_location;
            m_ghost_points[m_nof_ghost_points].flag = ghost->m_flag_1;

            m_nof_ghost_points++;

            potential_particle_insertion(ghost);
        }
        else if (ghost->m_flag_1 == INNER_UPDATED)
        {
            m_ghost_points[m_nof_ghost_points].location = ghost->m_location;
            m_ghost_points[m_nof_ghost_points].flag = INNER;

            m_nof_ghost_points++;
        }
        else if (ghost->m_flag_1 == VP_MOVED_OUT)
        {
            // keep the ghost
            m_ghost_points[m_nof_ghost_points].location = ghost->m_location;
            m_ghost_points[m_nof_ghost_points].flag = INNER;

            m_nof_ghost_points++;
        }
        else if (ghost->m_flag_1 == GHOST_MOVED_IN || ghost->m_flag_1 == GHOST_MOVED_OUT)
        {
            //don't keep the ghost
        }
        else if (ghost->m_flag_1 == INNER) // first simulation step
        {
            m_ghost_points[m_nof_ghost_points].location = ghost->m_location;
            m_ghost_points[m_nof_ghost_points].flag = INNER;

            m_nof_ghost_points++;
        }
    }
}

// create ghost particles out of particles close to the border
void ZeroGradRectangle::create_ghost_particles(double distance)
{
    BOX_PRINTF("\n--- creating ghost particles (distance) ---\n\n");

    BOX_PRINTF("number of ghost points: %d\n", m_nof_ghost_points);

    // create ghost particles
    for (unsigned int i = 0; i < m_nof_ghost_points; i++)
    {
        m_sim->m_ghosts.current() = GhostParticle(); //reset
        m_sim->m_ghosts.current().m_location = m_ghost_points[i].location;
        m_sim->m_ghosts.current().m_flag_1 = m_ghost_points[i].flag;

        m_sim->m_ghosts.counterInc();
    }

    reset();
}

// create ghost particles
void ZeroGradRectangle::createGhostParticles()
{
    create_ghost_particles_zero_grad();
}

// create zero gradient ghost particles
void ZeroGradRectangle::create_ghost_particles_zero_grad()
{
    BOX_PRINTF("\n--- creating ghost particles ---\n\n");

    BOX_PRINTF("number of ghost points: %d\n", m_nof_ghost_points);

    assert(m_nof_ghost_points > 0);

    // create ghost particles
    for (unsigned int i = 0; i < m_nof_ghost_points; i++)
    {
        m_sim->m_ghosts.current() = GhostParticle(); //reset
        m_sim->m_ghosts.current().m_location = m_ghost_points[i].location;
        m_sim->m_ghosts.current().m_flag_1 = m_ghost_points[i].flag;

        BOX_PRINTF("ghost: (%f, %f)\n", m_sim->m_ghosts.current().m_location.x, m_sim->m_ghosts.current().m_location.y);

        m_sim->m_ghosts.counterInc();
    }
}

// helper function
static void assign_values(GhostParticle* ghost)
{
    HalfEdge* start = ghost->m_outerComponent;
    HalfEdge* current = ghost->m_outerComponent;
    VoronoiParticle* twin_vp = NULL;

    double density = 0;
    double pressure = 0;
    double v_x = 0;
    double v_y = 0;

    double length = 0;
    double current_length = 0;

    do
    {
        twin_vp = current->m_twin->m_vp;
        current_length = current->m_length;

        if (!twin_vp->isGhost())
        {
            density += twin_vp->m_density * current_length;
            pressure += twin_vp->m_pressure * current_length;

            v_x += twin_vp->m_velocity.x * current_length;
            v_y += twin_vp->m_velocity.y * current_length;

            length += current_length;
        }

        current = current->m_next;

    } while (current != start);

    if (length != 0)
    {
        ghost->m_density = density / length;
        ghost->m_pressure = pressure / length;
        ghost->m_velocity.x = v_x / length;
        ghost->m_velocity.y = v_y / length;

    }
    else
    {
        ghost->m_density = 1;
        ghost->m_pressure = 1;
        ghost->m_velocity.x = 0;
        ghost->m_velocity.y = 0;
    }

    // calculate cell centers and volume
    ghost->updateCellVolumeAndCenter();

    BOX_PRINTF("\tdensity: %f\n\tpressure: %f\n\tv_x: %f\n\tv_y:%f\n", ghost->m_density, ghost->m_pressure,
               ghost->m_velocity.x, ghost->m_velocity.y);
    BOX_PRINTF("\t\tcell volume: %f\n\t\tcell center: ( %f , %f )\n\n", ghost->m_cellVolume,
               ghost->m_cellCenterOfVolume.x, ghost->m_cellCenterOfVolume.y);
}

// copy the primitive variables to the ghost particles
void ZeroGradRectangle::copyPrimitiveVariablesCellCentersAndVolume()
{
    // primitive variables: weighted arithmetic mean

    BOX_PRINTF("\n---copying primitive variables---\n\n");

    GhostParticle* ghost = NULL;

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        ghost = &m_sim->m_ghosts[i];

        if (!IS_FIXED(ghost->m_flag_1))
        {
            ghost->update_half_edges();
        }

        if (ghost->closed_and_inner())
        {
            assert(!IS_FIXED(ghost->m_flag_1));

            ghost->m_flag_1 = VERY_INNER;
            assign_values(ghost);
        }
        else
        {
            BOX_PRINTF("\touter ghost particle; values stay zero\n\n");
        }
    }
}

// copy the gradients to the ghost particles
void ZeroGradRectangle::copyGradients()
{
    // zero gradient

#ifndef NDEBUG

    GhostParticle* ghost = NULL;

    for(unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        ghost = &m_sim->m_ghosts[i];
        assert(ghost->m_gradientVelocityX.x == 0);
        assert(ghost->m_gradientVelocityX.y == 0);
        assert(ghost->m_gradientVelocityY.x == 0);
        assert(ghost->m_gradientVelocityY.y == 0);
        assert(ghost->m_gradientPressure.x == 0);
        assert(ghost->m_gradientPressure.y == 0);
        assert(ghost->m_gradientDensity.x == 0);
        assert(ghost->m_gradientDensity.y == 0);
    }

#endif

}

static void assign_particle_velocity_inner(GhostParticle* ghost)
{
    HalfEdge* start = ghost->m_outerComponent;
    HalfEdge* current = ghost->m_outerComponent;
    VoronoiParticle* twin_vp = NULL;

    double p_v_x = 0;
    double p_v_y = 0;

    double length = 0;
    double current_length = 0;

    do
    {
        twin_vp = current->m_twin->m_vp;
        current_length = current->m_length;

        if (!twin_vp->isGhost() && !current->is_zero())
        {
            p_v_x += twin_vp->m_particleVelocity.x * current_length;
            p_v_y += twin_vp->m_particleVelocity.y * current_length;

            length += current_length;
        }

        current = current->m_next;

    } while (current != start);

    if (length != 0)
    {
        ghost->m_particleVelocity.x = p_v_x / length;
        ghost->m_particleVelocity.y = p_v_y / length;

        ghost->m_flag_1 = INNER_UPDATED;
    }
    else
    {
        ghost->m_flag_1 = INNER;
    }

    BOX_PRINTF("\tparticle_velocity_x: %f\n\tparticle_velocity_y: %f\n", ghost->m_particleVelocity.x,
               ghost->m_particleVelocity.y);
}

static void assign_particle_velocity(GhostParticle* ghost)
{
    HalfEdge* start = ghost->m_outerComponent;
    HalfEdge* current = ghost->m_outerComponent;
    VoronoiParticle* twin_vp = NULL;

    double p_v_x = 0;
    double p_v_y = 0;

    double length = 0;
    double current_length = 0;

    do
    {
        twin_vp = current->m_twin->m_vp;
        current_length = current->m_length;

        if (twin_vp->m_flag_1 == INNER_UPDATED && !current->is_zero())
        {
            p_v_x += twin_vp->m_particleVelocity.x * current_length;
            p_v_y += twin_vp->m_particleVelocity.y * current_length;

            length += current_length;
        }

        current = current->m_next;

    } while (current != start);

    if (length != 0)
    {
        ghost->m_particleVelocity.x = p_v_x / length;
        ghost->m_particleVelocity.y = p_v_y / length;

        ghost->m_flag_1 = INNER_IN_PROGRESS;
    }

    BOX_PRINTF("\tparticle_velocity_x: %f\n\tparticle_velocity_y: %f\n", ghost->m_particleVelocity.x,
               ghost->m_particleVelocity.y);
}

// assign velocities to the ghost particles
void ZeroGradRectangle::copyParticleVelocities()
{
    BOX_PRINTF("\n--- copying particle velocities ---\n\n");

    GhostParticle* ghost;

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        ghost = &m_sim->m_ghosts[i];

        BOX_PRINTF("ghost: ( %f , %f )\n", ghost->m_location.x, ghost->m_location.y);

        // outer ghost particle
        if (!(ghost->m_flag_1 == VERY_INNER))
        {
            BOX_PRINTF("\touter ghost particle; values stay zero\n\n");
        }
        else //inner ghost particle
        {
            assign_particle_velocity_inner(ghost); //flag is set to INNER or INNER_UPDATED
        }
    }

    bool repeat = true;
    int count = 0;

    while (repeat)
    {
        count++;
        assert(count < 20);

        for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
        {
            ghost = &m_sim->m_ghosts[i];

            if (!(ghost->m_flag_1 == INNER_UPDATED) && !IS_FIXED(ghost->m_flag_1))
            {
                assign_particle_velocity(ghost);
            }
        }

        repeat = false;

        for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
        {
            ghost = &m_sim->m_ghosts[i];

            if (ghost->m_flag_1 == INNER_IN_PROGRESS)
            {
                ghost->m_flag_1 = INNER_UPDATED;
                repeat = true;
            }
        }
    }
}

// move a particle which is outside the boundary
void ZeroGradRectangle::handle_outside_boundary(VoronoiParticle* vp, unsigned int& iterator)
{
    BOX_PRINTF("\n--- handling outside boundary ---\n\n");

    m_sim->m_ghosts.current().m_location = vp->m_location;
    m_sim->m_ghosts.current().m_flag_1 = VP_MOVED_OUT;

    BOX_PRINTF("VP => ghost (%f, %f)\n", m_sim->m_ghosts.current().m_location.x,
               m_sim->m_ghosts.current().m_location.y);

    m_sim->m_ghosts.counterInc();

    m_sim->m_particles.del(vp);
    iterator--;
}

// move the ghost particles
void ZeroGradRectangle::move_ghost_particles()
{
    BOX_PRINTF("\n--- moving ghost particles ---\n\n");

    GhostParticle* ghost = NULL;
    VoronoiParticle* vp = NULL;

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        ghost = &m_sim->m_ghosts[i];

        if (!(ghost->m_flag_1 == VP_MOVED_OUT))
        {
            ghost->move();

            if (isInside(&ghost->m_location)) //ghost moved inside
            {
                ghost->m_flag_1 = GHOST_MOVED_IN;

                vp = &m_sim->m_particles.current();
                *vp = VoronoiParticle(); //reset

                // assign members
                //

                vp->m_location = ghost->m_location;

                //n o flux, primitive variables are the same as in the previous step
                vp->m_density = ghost->m_density;
                vp->m_pressure = ghost->m_pressure;
                vp->m_velocity = ghost->m_velocity;

                vp->m_particleVelocity = ghost->m_particleVelocity;
                vp->m_cellCenterOfVolume = ghost->m_cellCenterOfVolume;
                vp->m_cellVolume = ghost->m_cellVolume;

                // calculate the conserved variables
                vp->calculateConservedVariables();

                BOX_PRINTF("ghost => VP (%f, %f)\n", vp->m_location.x, vp->m_location.y);

                m_sim->m_particles.counterInc();
            }

            if (!IS_FIXED(ghost->m_flag_1) && !is_inside_ghost_box(&ghost->m_location)) //ghost moved outside
            {
                ghost->m_flag_1 = GHOST_MOVED_OUT;
            }
        }
    }
}

// ini function
void ZeroGradRectangle::sitesIniGrid(unsigned int particles_x_dimension, unsigned int particles_y_dimension)
{
    if (particles_x_dimension * particles_y_dimension > m_sim->getMaxNumberOfParticles())
    {
        fprintf(stderr,
                "ERROR in HybridRectangle::sitesIniGrid: too many particles, create a simulation with more maximum number of particles\n");
        exit(1);
    }

    int j_min = -(int) LAYERS;
    int j_max = (int) particles_y_dimension + (int) LAYERS - 1;

    int i_min = -(int) LAYERS;
    int i_max = (int) particles_x_dimension + (int) LAYERS - 1;

    for (int j = j_min; j <= j_max; j++)
    {
        for (int i = i_min; i <= i_max; i++)
        {
            // Voronoi particles
            if (j >= 0 && j < (int) particles_y_dimension && i >= 0 && i < (int) particles_x_dimension)
            {
                m_sim->m_particles.current().m_location.x = getXmin() + (i + 0.5) * m_dh;
                m_sim->m_particles.current().m_location.y = getYmin() + (j + 0.5) * m_dv;
                m_sim->m_particles.counterInc();
            }
            else //ghost particles
            {
                m_ghost_points[m_nof_ghost_points].location.x = getXmin() + (i + 0.5) * m_dh;
                m_ghost_points[m_nof_ghost_points].location.y = getYmin() + (j + 0.5) * m_dv;

                if (is_inside_ghost_box(&m_ghost_points[m_nof_ghost_points].location))
                {
                    m_ghost_points[m_nof_ghost_points].flag = INNER;
                }
                else
                {
                    if (is_fixed_bottom(&m_ghost_points[m_nof_ghost_points].location))
                    {
                        m_ghost_points[m_nof_ghost_points].flag = FIXED_BOTTOM;
                    }
                    else if (is_fixed_top(&m_ghost_points[m_nof_ghost_points].location))
                    {
                        m_ghost_points[m_nof_ghost_points].flag = FIXED_TOP;
                    }
                    else if (is_fixed_left(&m_ghost_points[m_nof_ghost_points].location))
                    {
                        m_ghost_points[m_nof_ghost_points].flag = FIXED_LEFT;
                    }
                    else if (is_fixed_right(&m_ghost_points[m_nof_ghost_points].location))
                    {
                        m_ghost_points[m_nof_ghost_points].flag = FIXED_RIGHT;
                    }
                    else
                    {
                        m_ghost_points[m_nof_ghost_points].flag = FIXED;
                    }
                }

                m_nof_ghost_points++;
            }
        }
    }
}