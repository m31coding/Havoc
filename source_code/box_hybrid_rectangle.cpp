#include "config.h"
#include "debug.h"
#include "box_hybrid_rectangle.h"

// constructor
HybridRectangle::HybridRectangle(Sim* sim, double Xmin, double Xmax, double Ymin, double Ymax)
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
        m_rightBorder(sim->getMaxNumberOfParticles() * factorGhostParticles),
        m_topBorder(sim->getMaxNumberOfParticles() * factorGhostParticles),
        m_leftBorder(sim->getMaxNumberOfParticles() * factorGhostParticles),
        m_bottomBorder(sim->getMaxNumberOfParticles() * factorGhostParticles),
        m_top_right_border(1),
        m_top_left_border(1),
        m_bottom_right_border(1),
        m_bottom_left_border(1)
{
    m_sim = sim;

    // Layer of ten particles...
    m_distance_x = 10 * m_deltaX / sqrt(sim->getMaxNumberOfParticles() * m_deltaX / m_deltaY);
    m_distance_y = 10 * m_deltaY / sqrt(sim->getMaxNumberOfParticles() * m_deltaY / m_deltaX);

    // ... or minimum distance of 10% of the box length
    if (m_distance_x < 0.1 * m_deltaX)
    {
        m_distance_x = 0.1 * m_deltaX;
    }

    if (m_distance_y < 0.1 * m_deltaY)
    {
        m_distance_y = 0.1 * m_deltaY;
    }

    BOX_PRINTF("\nbottom right: ( %f , %f )\n", m_bottomRight.x, m_bottomRight.y);
    BOX_PRINTF("top right: ( %f , %f )\n", m_topRight.x, m_topRight.y);
    BOX_PRINTF("top left: ( %f , %f )\n", m_topLeft.x, m_topLeft.y);
    BOX_PRINTF("bottom left: ( %f , %f )\n", m_bottomLeft.x, m_bottomLeft.y);
    BOX_PRINTF("delta x: %f\ndelta y: %f\n", m_deltaX, m_deltaY);

    // horizontal distance of the particles
    m_dh = m_deltaX / Constants::getNOF_PARTICLES_X();

    // vertical distance of the particles
    m_dv = m_deltaY / Constants::getNOF_PARTICLES_Y();

    m_DV = m_dh * m_dv;
}

// destructor
HybridRectangle::~HybridRectangle()
{

}

// reset the arrays
void HybridRectangle::reset()
{
    m_rightBorder.reset();
    m_topBorder.reset();
    m_leftBorder.reset();
    m_bottomBorder.reset();
}

// check whether a point is inside the box
bool HybridRectangle::isInside(Point* point)
{
    return (point->x <= m_Xmax && point->x >= m_Xmin && point->y <= m_Ymax && point->y >= m_Ymin);
}

// handle a point
Point HybridRectangle::handle_point(Point* point, myArray<VoronoiParticle*>* segment)
{
    Point result;

    if (segment == &m_topBorder)
    {
        result.x = point->x;
        result.y = 2 * getYmax() - point->y;
    }
    else if (segment == &m_bottomBorder)
    {
        result.x = point->x;
        result.y = 2 * getYmin() - point->y;
    }
    else if (segment == &m_leftBorder)
    {
        result.x = point->x + getDeltaX();
        result.y = point->y;
    }
    else if (segment == &m_rightBorder)
    {
        result.x = point->x - getDeltaX();
        result.y = point->y;
    }
    else if (segment == &m_top_right_border)
    {
        result.x = point->x - getDeltaX();
        result.y = 2 * getYmax() - point->y;
    }
    else if (segment == &m_top_left_border)
    {
        result.x = point->x + getDeltaX();
        result.y = 2 * getYmax() - point->y;
    }
    else if (segment == &m_bottom_right_border)
    {
        result.x = point->x - getDeltaX();
        result.y = 2 * getYmin() - point->y;
    }
    else if (segment == &m_bottom_left_border)
    {
        result.x = point->x + getDeltaX();
        result.y = 2 * getYmin() - point->y;
    }
    else
    {
        fprintf(stderr, "ERROR in HybridRectangle::mirrorPoint: segment is not part of the box\n");
        exit(1);
    }

    return result;
}

Vector HybridRectangle::handle_vector(Vector* vector, myArray<VoronoiParticle*>* segment)
{
    Vector result;

    if (segment == &m_topBorder || segment == &m_bottomBorder)
    {
        result.x = vector->x;
        result.y = -vector->y;
    }
    else if (segment == &m_leftBorder || segment == &m_rightBorder)
    {
        result = *vector;
    }
    else if (segment == &m_top_right_border || segment == &m_top_left_border || segment == &m_bottom_right_border ||
             segment == &m_bottom_left_border)
    {
        result.x = vector->x;
        result.y = -vector->y;
    }
    else
    {
        fprintf(stderr, "ERROR in HybridRectangle::mirrorPoint: segment is not part of the box\n");
        exit(1);
    }

    return result;
}

// create ghost particles out of particles close to the border
void HybridRectangle::create_ghost_particles(double distance_x, double distance_y)
{
    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < m_sim->m_particles.counter(); i++)
    {
        particle = &m_sim->m_particles[i];

        // top
        if (m_Ymax - particle->m_location.y <= distance_y)
        {
            m_sim->m_ghosts.current().m_location = handle_point(&particle->m_location, &m_topBorder);
            m_sim->m_ghosts.current().m_realFriend = particle;
            m_sim->m_ghosts.current().m_copiedOver = &m_topBorder;

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in HybridRectangle::create_ghost_particles:\n particle outside of the computational domain after mirroring: ( %f , %f )\n",
                        particle->m_location.x, particle->m_location.y);
                exit(1);
            }

            m_sim->m_ghosts.counterInc();
        }

        // bottom
        if (particle->m_location.y - m_Ymin <= distance_y)
        {
            m_sim->m_ghosts.current().m_location = handle_point(&particle->m_location, &m_bottomBorder);
            m_sim->m_ghosts.current().m_realFriend = particle;
            m_sim->m_ghosts.current().m_copiedOver = &m_bottomBorder;

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in HybridRectangle::create_ghost_particles:\n particle outside of the computational domain after mirroring: ( %f , %f )\n",
                        particle->m_location.x, particle->m_location.y);
                exit(1);
            }

            m_sim->m_ghosts.counterInc();
        }

        // left
        if (particle->m_location.x - m_Xmin <= distance_x)
        {
            m_sim->m_ghosts.current().m_location = handle_point(&particle->m_location, &m_leftBorder);
            m_sim->m_ghosts.current().m_realFriend = particle;
            m_sim->m_ghosts.current().m_copiedOver = &m_leftBorder;

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in HybridRectangle::create_ghost_particles:\n particle outside of the computational domain after mirroring: ( %f , %f )\n",
                        particle->m_location.x, particle->m_location.y);
                exit(1);
            }

            m_sim->m_ghosts.counterInc();
        }

        // right
        if (m_Xmax - particle->m_location.x <= distance_x)
        {
            m_sim->m_ghosts.current().m_location = handle_point(&particle->m_location, &m_rightBorder);
            m_sim->m_ghosts.current().m_realFriend = particle;
            m_sim->m_ghosts.current().m_copiedOver = &m_rightBorder;

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in HybridRectangle::create_ghost_particles:\n particle outside of the computational domain after mirroring: ( %f , %f )\n",
                        particle->m_location.x, particle->m_location.y);
                exit(1);
            }

            m_sim->m_ghosts.counterInc();
        }

        // top right corner
        if (m_Ymax - particle->m_location.y <= distance_y && m_Xmax - particle->m_location.x <= distance_x)
        {
            Point temp = handle_point(&particle->m_location, &m_rightBorder);
            m_sim->m_ghosts.current().m_location = handle_point(&temp, &m_topBorder);
            m_sim->m_ghosts.current().m_realFriend = particle;
            m_sim->m_ghosts.current().m_copiedOver = &m_top_right_border;

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in HybridRectangle::create_ghost_particles:\n particle outside of the computational domain after mirroring: ( %f , %f )\n",
                        particle->m_location.x, particle->m_location.y);
                exit(1);
            }

            m_sim->m_ghosts.counterInc();
        }

        // bottom right corner
        if (particle->m_location.y - m_Ymin <= distance_y && m_Xmax - particle->m_location.x <= distance_x)
        {
            Point temp = handle_point(&particle->m_location, &m_rightBorder);
            m_sim->m_ghosts.current().m_location = handle_point(&temp, &m_bottomBorder);
            m_sim->m_ghosts.current().m_realFriend = particle;
            m_sim->m_ghosts.current().m_copiedOver = &m_bottom_right_border;

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in HybridRectangle::create_ghost_particles:\n particle outside of the computational domain after mirroring: ( %f , %f )\n",
                        particle->m_location.x, particle->m_location.y);
                exit(1);
            }

            m_sim->m_ghosts.counterInc();
        }

        // bottom left corner
        if (particle->m_location.y - m_Ymin <= distance_y && particle->m_location.x - m_Xmin <= distance_x)
        {
            Point temp = handle_point(&particle->m_location, &m_leftBorder);
            m_sim->m_ghosts.current().m_location = handle_point(&temp, &m_bottomBorder);
            m_sim->m_ghosts.current().m_realFriend = particle;
            m_sim->m_ghosts.current().m_copiedOver = &m_bottom_left_border;

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in HybridRectangle::create_ghost_particles:\n particle outside of the computational domain after mirroring: ( %f , %f )\n",
                        particle->m_location.x, particle->m_location.y);
                exit(1);
            }

            m_sim->m_ghosts.counterInc();
        }

        // top left corner
        if (m_Ymax - particle->m_location.y <= distance_y && particle->m_location.x - m_Xmin <= distance_x)
        {
            Point temp = handle_point(&particle->m_location, &m_leftBorder);
            m_sim->m_ghosts.current().m_location = handle_point(&temp, &m_topBorder);
            m_sim->m_ghosts.current().m_realFriend = particle;
            m_sim->m_ghosts.current().m_copiedOver = &m_top_left_border;

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in HybridRectangle::create_ghost_particles:\n particle outside of the computational domain after mirroring: ( %f , %f )\n",
                        particle->m_location.x, particle->m_location.y);
                exit(1);
            }

            m_sim->m_ghosts.counterInc();
        }
    }
}

void HybridRectangle::create_ghost_particles(double distance)
{
    create_ghost_particles(distance, distance);
}

// determine ghost particles
void HybridRectangle::determineGhostParticles()
{
    return;
}

// create ghost particles, previous determined via determineGhostParticles
void HybridRectangle::createGhostParticles()
{
    create_ghost_particles(m_distance_x, m_distance_y);
}

// copy the primitive variables to the ghost particles
void HybridRectangle::copyPrimitiveVariablesCellCentersAndVolume()
{
    VAR_PRINTF("\n---copying primitive variables---\n\n");

    GhostParticle* ghost = NULL; // temporary ghost

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        VAR_PRINTF("ghost particle: ( %f , %f )\n", m_sim->m_ghosts[i].m_location.x, m_sim->m_ghosts[i].m_location.y);

        ghost = &m_sim->m_ghosts[i];

        // primitive variables
        ghost->m_density = ghost->m_realFriend->m_density;
        ghost->m_pressure = ghost->m_realFriend->m_pressure;
        ghost->m_velocity = handle_vector(&ghost->m_realFriend->m_velocity, ghost->m_copiedOver);

        // cell volume
        ghost->m_cellVolume = ghost->m_realFriend->m_cellVolume;

        // cell center
        ghost->m_cellCenterOfVolume = handle_point(&ghost->m_realFriend->m_cellCenterOfVolume, ghost->m_copiedOver);

        VAR_PRINTF("\t density: %f\n", ghost->m_density);
        VAR_PRINTF("\t velocity: ( %f , %f )\n", ghost->m_velocity.x, ghost->m_velocity.y);
        VAR_PRINTF("\t pressure: %f\n", ghost->m_pressure);

        VAR_PRINTF("\t\treal friend: ( %f , %f )\n", ghost->m_realFriend->m_location.x,
                   ghost->m_realFriend->m_location.y);
        VAR_PRINTF("\t\t\t density: %f\n", ghost->m_realFriend->m_density);
        VAR_PRINTF("\t\t\t velocity: ( %f , %f )\n", ghost->m_realFriend->m_velocity.x,
                   ghost->m_realFriend->m_velocity.y);
        VAR_PRINTF("\t\t\t pressure: %f\n", ghost->m_realFriend->m_pressure);
    }
}

// calculate the ghost gradients
Vector HybridRectangle::ghost_gradient_x(Vector* real_gradient, myArray<VoronoiParticle*>* segment)
{
    Vector result;

    if (segment == &m_topBorder || segment == &m_bottomBorder || segment == &m_top_right_border ||
        segment == &m_top_left_border || segment == &m_bottom_right_border || segment == &m_bottom_left_border)
    {
        result.x = real_gradient->x;
        result.y = -real_gradient->y;
    }
    else if (segment == &m_leftBorder || segment == &m_rightBorder)
    {
        result.x = real_gradient->x;
        result.y = real_gradient->y;
    }
    else
    {
        fprintf(stderr, "ERROR in box_refelcting_rectangle::ghost_gradient_y: segment is not part of the box\n");
        exit(1);
    }

    return result;
}

Vector HybridRectangle::ghost_gradient_y(Vector* real_gradient, myArray<VoronoiParticle*>* segment)
{
    Vector result;

    if (segment == &m_topBorder || segment == &m_bottomBorder || segment == &m_top_right_border ||
        segment == &m_top_left_border || segment == &m_bottom_right_border || segment == &m_bottom_left_border)
    {
        result.x = real_gradient->x;
        result.y = -real_gradient->y;
        result *= -1;
    }

    else if (segment == &m_leftBorder || segment == &m_rightBorder)
    {
        result.x = real_gradient->x;
        result.y = real_gradient->y;
    }
    else
    {
        fprintf(stderr, "ERROR in box_refelcting_rectangle::ghost_gradient_y: segment is not part of the box\n");
        exit(1);
    }

    return result;
}

// copy the gradients to the ghost particles
void HybridRectangle::copyGradients()
{
    GRAD_PRINTF("\n---copying gradients---\n\n");

    GhostParticle* ghost = NULL; //temporary ghost

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        ghost = &m_sim->m_ghosts[i];

        GRAD_PRINTF("ghost particle: ( %f , %f )\n", ghost->m_location.x, ghost->m_location.y);

        ghost->m_gradientDensity = handle_vector(&ghost->m_realFriend->m_gradientDensity, ghost->m_copiedOver);
        ghost->m_gradientPressure = handle_vector(&ghost->m_realFriend->m_gradientPressure, ghost->m_copiedOver);
        ghost->m_gradientVelocityX = ghost_gradient_x(&ghost->m_realFriend->m_gradientVelocityX, ghost->m_copiedOver);
        ghost->m_gradientVelocityY = ghost_gradient_y(&ghost->m_realFriend->m_gradientVelocityY, ghost->m_copiedOver);

        GRAD_PRINTF("\tgrad(v_x): ( %f , %f)\n", ghost->m_gradientVelocityX.x, ghost->m_gradientVelocityX.y);
        GRAD_PRINTF("\tgrad(v_y): ( %f , %f)\n", ghost->m_gradientVelocityY.x, ghost->m_gradientVelocityY.y);
        GRAD_PRINTF("\tgrad(density): ( %f , %f)\n", ghost->m_gradientDensity.x, ghost->m_gradientDensity.y);
        GRAD_PRINTF("\tgrad(pressure): ( %f , %f)\n", ghost->m_gradientPressure.x, ghost->m_gradientPressure.y);

        GRAD_PRINTF("\t\treal friend: ( %f , %f )\n", ghost->m_realFriend->m_location.x,
                    ghost->m_realFriend->m_location.y);
        GRAD_PRINTF("\t\t\tgrad(v_x): ( %f , %f)\n", ghost->m_realFriend->m_gradientVelocityX.x,
                    ghost->m_realFriend->m_gradientVelocityX.y);
        GRAD_PRINTF("\t\t\tgrad(v_y): ( %f , %f)\n", ghost->m_realFriend->m_gradientVelocityY.x,
                    ghost->m_realFriend->m_gradientVelocityY.y);
        GRAD_PRINTF("\t\t\tgrad(density): ( %f , %f)\n", ghost->m_realFriend->m_gradientDensity.x,
                    ghost->m_realFriend->m_gradientDensity.y);
        GRAD_PRINTF("\t\t\tgrad(pressure): ( %f , %f)\n", ghost->m_realFriend->m_gradientPressure.x,
                    ghost->m_realFriend->m_gradientPressure.y);
    }
}

// assign velocities to the ghost particles
void HybridRectangle::copyParticleVelocities()
{
    GhostParticle* ghost = NULL;

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        ghost = &m_sim->m_ghosts[i];

        ghost->m_particleVelocity = handle_vector(&ghost->m_realFriend->m_particleVelocity, ghost->m_copiedOver);
    }
}

// initialization of m_particles, the graph is a cartesian grid
void HybridRectangle::sitesIniGrid(unsigned int particles_x_dimension, unsigned int particles_y_dimension)
{
    if (particles_x_dimension * particles_y_dimension > m_sim->getMaxNumberOfParticles())
    {
        fprintf(stderr,
                "ERROR in HybridRectangle::sitesIniGrid: too many particles, create a simulation with more maximum number of particles\n");
        exit(1);
    }

    // horizontal distance of the particles
    double dh = getDeltaX() / (particles_x_dimension);

    // vertical distance of the particles
    double dv = getDeltaY() / (particles_y_dimension);

    for (unsigned int j = 0; j < particles_y_dimension; j++)
    {
        for (unsigned int i = 0; i < particles_x_dimension; i++)
        {

            m_sim->m_particles.current().m_location.x = getXmin() + (i + 0.5) * dh;
            m_sim->m_particles.current().m_location.y = getYmin() + (j + 0.5) * dv;

            m_sim->m_particles.counterInc();
        }
    }
}

// move a particle which is outside the boundary
void HybridRectangle::handle_outside_boundary(VoronoiParticle* vp, unsigned int& iterator)
{
    if (vp->m_location.x > getXmax())
    {
        vp->m_location.x -= getDeltaX();
    }

    else if (vp->m_location.x < getXmin())
    {
        vp->m_location.x += getDeltaX();
    }
}