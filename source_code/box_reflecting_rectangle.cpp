#include "config.h"
#include "box_reflecting_rectangle.h"
#include "debug.h"
#include "box.h"
#include "r2.h"
#include "sim.h"

// constructor
ReflectingRectangle::ReflectingRectangle(Sim* sim, double Xmin, double Xmax, double Ymin, double Ymax)
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
        m_bottomBorder(sim->getMaxNumberOfParticles() * factorGhostParticles)
{
    m_sim = sim;
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
ReflectingRectangle::~ReflectingRectangle()
{

}

// reset the arrays
void ReflectingRectangle::reset()
{
    m_rightBorder.reset();
    m_topBorder.reset();
    m_leftBorder.reset();
    m_bottomBorder.reset();
}

// check whether a point is inside the box
bool ReflectingRectangle::isInside(Point* point)
{
    return (point->x <= m_Xmax && point->x >= m_Xmin && point->y <= m_Ymax && point->y >= m_Ymin);
}

// search for all neighbours of the particle which are inside the box
void ReflectingRectangle::searchNeighboursInside(VoronoiParticle* particle)
{
    NE_PRINTF("searching for neighbours: ( %f , %f )\n", particle->m_location.x, particle->m_location.y);

    HalfEdge* start = particle->m_outerComponent;
    HalfEdge* current = start;
    assert(start != NULL);

    do
    {
        if (!(current->m_twin->m_vp->isGhost()))
        {
            NE_PRINTF("\tneighbour: ( %f , %f)\n", current->m_twin->m_vp->m_location.x,
                      current->m_twin->m_vp->m_location.y);
            handleNeighbours(particle, current->m_twin->m_vp);
        }
        current = current->m_next;

    } while (current != start && current != NULL);

    if (current == NULL)
    {
        current = start->m_prev;

        while (current != NULL)
        {
            if (!(current->m_twin->m_vp->isGhost()))
            {
                NE_PRINTF("\tneighbour: ( %f , %f)\n", current->m_twin->m_vp->m_location.x,
                          current->m_twin->m_vp->m_location.y);
                handleNeighbours(particle, current->m_twin->m_vp);
            }
            current = current->m_prev;
        }
    }
}

// handle two neighbours
void ReflectingRectangle::handleNeighbours(VoronoiParticle* particle_outside, VoronoiParticle* particle_inside)
{
    myArray<VoronoiParticle*>* intersection = intersectOnce(particle_outside->m_location, particle_inside->m_location);
    assert(intersection != NULL); // there should be an intersection

    intersection->current() = particle_inside;
    intersection->counterInc();

    handleNeighboursOfNeighbours(particle_inside, intersection);
}

// handle the neighbours of the neighbours
void ReflectingRectangle::handleNeighboursOfNeighbours(VoronoiParticle* particle_inside,
                                                       myArray<VoronoiParticle*>* intersection)
{
    NE_PRINTF("searching for neighbours of neighbour: ( %f , %f )\n", particle_inside->m_location.x,
              particle_inside->m_location.y);

    HalfEdge* start = particle_inside->m_outerComponent;
    HalfEdge* current = start;

    assert(start != NULL);
    assert(graph::isClosed(start));

    do
    {
        if (!(current->m_twin->m_vp->isGhost()))
        {
            NE_PRINTF("\tneighbour of neighbour: ( %f , %f)\n", current->m_twin->m_vp->m_location.x,
                      current->m_twin->m_vp->m_location.y);

            intersection->current() = current->m_twin->m_vp;
            intersection->counterInc();
        }
        current = current->m_next;

    } while (current != start);
}

// determine ghost particles
void ReflectingRectangle::determineGhostParticles()
{
    BOX_PRINTF("\n---determining ghost particles---\n\n");

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        searchNeighboursInside(&m_sim->m_ghosts[i]);
    }
}

// create ghost particles
void ReflectingRectangle::createGhostParticles()
{
    BOX_PRINTF("\n---creating ghost particles---\n\n");

    assert(m_sim->m_ghosts.counter() == 0); // the ghost particle array should be empty

    // right border
    //
    for (unsigned int i = 0; i < m_rightBorder.counter(); i++)
    {
        // copy all particles only once over each border segment
        if (m_rightBorder[i]->m_copiedOver == &m_rightBorder)
        {
            continue;
        }

        // copy over right border
        m_sim->m_ghosts.current().m_location.x = -m_rightBorder[i]->m_location.x + 2 * m_Xmax;
        m_sim->m_ghosts.current().m_location.y = m_rightBorder[i]->m_location.y;
        m_sim->m_ghosts.current().m_realFriend = m_rightBorder[i];
        m_sim->m_ghosts.current().m_copiedOver = &m_rightBorder;

        BOX_PRINTF("right border: ghost particle: ( %f , %f ) <= ( %f, %f)\n", m_sim->m_ghosts.current().m_location.x,
                   m_sim->m_ghosts.current().m_location.y,
                   m_rightBorder[i]->m_location.x, m_rightBorder[i]->m_location.y);

        m_sim->m_ghosts.counterInc();
        m_rightBorder[i]->m_copiedOver = &m_rightBorder;
    }

    // top border
    //
    for (unsigned int j = 0; j < m_topBorder.counter(); j++)
    {
        // copy all particles only once over each border segment
        if (m_topBorder[j]->m_copiedOver == &m_topBorder)
        {
            continue;
        }

        // copy over top border
        m_sim->m_ghosts.current().m_location.x = m_topBorder[j]->m_location.x;
        m_sim->m_ghosts.current().m_location.y = -m_topBorder[j]->m_location.y + 2 * m_Ymax;
        m_sim->m_ghosts.current().m_realFriend = m_topBorder[j];
        m_sim->m_ghosts.current().m_copiedOver = &m_topBorder;

        BOX_PRINTF("top border: ghost particle: ( %f , %f ) <= ( %f, %f)\n", m_sim->m_ghosts.current().m_location.x,
                   m_sim->m_ghosts.current().m_location.y,
                   m_topBorder[j]->m_location.x, m_topBorder[j]->m_location.y);

        m_sim->m_ghosts.counterInc();
        m_topBorder[j]->m_copiedOver = &m_topBorder;
    }

    // left border
    //
    for (unsigned int k = 0; k < m_leftBorder.counter(); k++)
    {
        // copy all particles only once over each border segment
        if (m_leftBorder[k]->m_copiedOver == &m_leftBorder)
        {
            continue;
        }

        // copy over left border
        m_sim->m_ghosts.current().m_location.x = -m_leftBorder[k]->m_location.x + 2 * m_Xmin;
        m_sim->m_ghosts.current().m_location.y = m_leftBorder[k]->m_location.y;
        m_sim->m_ghosts.current().m_realFriend = m_leftBorder[k];
        m_sim->m_ghosts.current().m_copiedOver = &m_leftBorder;

        BOX_PRINTF("left border: ghost particle: ( %f , %f ) <= ( %f, %f)\n", m_sim->m_ghosts.current().m_location.x,
                   m_sim->m_ghosts.current().m_location.y,
                   m_leftBorder[k]->m_location.x, m_leftBorder[k]->m_location.y);

        m_sim->m_ghosts.counterInc();
        m_leftBorder[k]->m_copiedOver = &m_leftBorder;
    }

    // bottom border
    //
    for (unsigned int l = 0; l < m_bottomBorder.counter(); l++)
    {
        // copy all particles only once over each border segment
        if (m_bottomBorder[l]->m_copiedOver == &m_bottomBorder)
        {
            continue;
        }

        // copy over bottom border
        m_sim->m_ghosts.current().m_location.x = m_bottomBorder[l]->m_location.x;
        m_sim->m_ghosts.current().m_location.y = -m_bottomBorder[l]->m_location.y + 2 * m_Ymin;
        m_sim->m_ghosts.current().m_realFriend = m_bottomBorder[l];
        m_sim->m_ghosts.current().m_copiedOver = &m_bottomBorder;

        BOX_PRINTF("bottom border: ghost particle: ( %f , %f ) <= ( %f, %f)\n", m_sim->m_ghosts.current().m_location.x,
                   m_sim->m_ghosts.current().m_location.y,
                   m_bottomBorder[l]->m_location.x, m_bottomBorder[l]->m_location.y);

        m_sim->m_ghosts.counterInc();
        m_bottomBorder[l]->m_copiedOver = &m_bottomBorder;
    }
}

// copy the primitive variables to the ghost particles
void ReflectingRectangle::copyPrimitiveVariablesCellCentersAndVolume()
{
    VAR_PRINTF("\n---copying primitive variables---\n\n");

    GhostParticle* ghost = NULL; // temporary ghost

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        VAR_PRINTF("ghost particle: ( %f , %f )\n", m_sim->m_ghosts[i].m_location.x, m_sim->m_ghosts[i].m_location.y);

        ghost = &m_sim->m_ghosts[i];

        ghost->m_density = ghost->m_realFriend->m_density;
        ghost->m_pressure = ghost->m_realFriend->m_pressure;

        // cell volume
        ghost->m_cellVolume = ghost->m_realFriend->m_cellVolume;

        if (ghost->m_copiedOver == &m_topBorder ||
            ghost->m_copiedOver == &m_bottomBorder) // copied over top or bottom border
        {
            ghost->m_velocity.x = ghost->m_realFriend->m_velocity.x;
            ghost->m_velocity.y = -ghost->m_realFriend->m_velocity.y;

            if (ghost->m_copiedOver == &m_topBorder) // top border
            {
                ghost->m_cellCenterOfVolume.x = ghost->m_realFriend->m_cellCenterOfVolume.x;
                ghost->m_cellCenterOfVolume.y = 2 * getYmax() - ghost->m_realFriend->m_cellCenterOfVolume.y;
            }

            else // bottom border
            {
                ghost->m_cellCenterOfVolume.x = ghost->m_realFriend->m_cellCenterOfVolume.x;
                ghost->m_cellCenterOfVolume.y = 2 * getYmin() - ghost->m_realFriend->m_cellCenterOfVolume.y;
            }
        }

        else // copied over left or right border
        {
            ghost->m_velocity.x = -ghost->m_realFriend->m_velocity.x;
            ghost->m_velocity.y = ghost->m_realFriend->m_velocity.y;

            if (ghost->m_copiedOver == &m_leftBorder) // left border
            {
                ghost->m_cellCenterOfVolume.x = 2 * getXmin() - ghost->m_realFriend->m_cellCenterOfVolume.x;
                ghost->m_cellCenterOfVolume.y = ghost->m_realFriend->m_cellCenterOfVolume.y;
            }
            else // right border
            {
                ghost->m_cellCenterOfVolume.x = 2 * getXmax() - ghost->m_realFriend->m_cellCenterOfVolume.x;
                ghost->m_cellCenterOfVolume.y = ghost->m_realFriend->m_cellCenterOfVolume.y;
            }
        }

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

// copy the gradients to the ghost particles
void ReflectingRectangle::copyGradients()
{
    GRAD_PRINTF("\n---copying gradients---\n\n");

    GhostParticle* ghost = NULL; // temporary ghost

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        ghost = &m_sim->m_ghosts[i];

        GRAD_PRINTF("ghost particle: ( %f , %f )\n", ghost->m_location.x, ghost->m_location.y);

        ghost->m_gradientDensity = mirrorVector(&ghost->m_realFriend->m_gradientDensity, ghost->m_copiedOver);
        ghost->m_gradientPressure = mirrorVector(&ghost->m_realFriend->m_gradientPressure, ghost->m_copiedOver);
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
void ReflectingRectangle::copyParticleVelocities()
{
    GhostParticle* ghost = NULL;

    for (unsigned int i = 0; i < m_sim->m_ghosts.counter(); i++)
    {
        ghost = &m_sim->m_ghosts[i];

        ghost->m_particleVelocity = mirrorVector(&ghost->m_realFriend->m_particleVelocity,
                                                 ghost->m_realFriend->m_copiedOver);
    }
}

// calculate the ghost gradients
Vector ReflectingRectangle::ghost_gradient_x(Vector* real_gradient, myArray<VoronoiParticle*>* segment)
{
    Vector result;

    if (segment == &m_topBorder || segment == &m_bottomBorder)
    {
        result.x = real_gradient->x;
        result.y = -real_gradient->y;
    }
    else if (segment == &m_leftBorder || segment == &m_rightBorder)
    {
        result.x = -real_gradient->x;
        result.y = real_gradient->y;
        result *= -1;
    }
    else
    {
        fprintf(stderr, "ERROR in box_refelcting_rectangle::ghost_gradient_y: segment is not part of the box\n");
        exit(1);
    }

    return result;
}


Vector ReflectingRectangle::ghost_gradient_y(Vector* real_gradient, myArray<VoronoiParticle*>* segment)
{
    Vector result;

    if (segment == &m_topBorder || segment == &m_bottomBorder)
    {
        result.x = real_gradient->x;
        result.y = -real_gradient->y;
        result *= -1;
    }
    else if (segment == &m_leftBorder || segment == &m_rightBorder)
    {
        result.x = -real_gradient->x;
        result.y = real_gradient->y;
    }
    else
    {
        fprintf(stderr, "ERROR in box_refelcting_rectangle::ghost_gradient_y: segment is not part of the box\n");
        exit(1);
    }

    return result;
}

// mirror a vector on a box segment
Vector ReflectingRectangle::mirrorVector(Vector* vector, myArray<VoronoiParticle*>* segment)
{
    Vector result;

    if (segment == &m_topBorder || segment == &m_bottomBorder)
    {
        result.x = vector->x;
        result.y = -vector->y;
    }
    else if (segment == &m_leftBorder || segment == &m_rightBorder)
    {
        result.x = -vector->x;
        result.y = vector->y;
    }
    else
    {
        fprintf(stderr, "ERROR in ReflectingRectangle::mirrorVector: segment is not part of the box\n");
        exit(1);
    }

    return result;
}

// mirror a point on a box segment
Point ReflectingRectangle::mirrorPoint(Point* point, myArray<VoronoiParticle*>* segment)
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
        result.x = 2 * getXmin() - point->x;
        result.y = point->y;
    }
    else if (segment == &m_rightBorder)
    {
        result.x = 2 * getXmax() - point->x;
        result.y = point->y;
    }
    else
    {
        fprintf(stderr, "ERROR in ReflectingRectangle::mirrorPoint: segment is not part of the box\n");
        exit(1);
    }

    return result;
}

// mirror particles with a distance to a wall that is less than or equal to "distance".
void ReflectingRectangle::create_ghost_particles(double distance)
{
    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < m_sim->m_particles.counter(); i++)
    {
        particle = &m_sim->m_particles[i];

        // mirror on ymax
        if (m_Ymax - particle->m_location.y <= distance)
        {
            m_sim->m_ghosts.current().m_location = mirrorPoint(&particle->m_location, &m_topBorder);
            m_sim->m_ghosts.counterInc();

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in ReflectingRectangle::mirrorParticles:\n particle outside of the computational domain after mirroring: ( %f , %f ) -> ( %f , %f ) \n",
                        particle->m_location.x, particle->m_location.y, m_sim->m_ghosts.current().m_location.x,
                        m_sim->m_ghosts.current().m_location.y);
                exit(1);
            }
        }

        // mirror on ymin
        if (particle->m_location.y - m_Ymin <= distance)
        {
            m_sim->m_ghosts.current().m_location = mirrorPoint(&particle->m_location, &m_bottomBorder);
            m_sim->m_ghosts.counterInc();

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in ReflectingRectangle::mirrorParticles:\n particle outside of the computational domain after mirroring: ( %f , %f ) -> ( %f , %f ) \n",
                        particle->m_location.x, particle->m_location.y, m_sim->m_ghosts.current().m_location.x,
                        m_sim->m_ghosts.current().m_location.y);
                exit(1);
            }
        }

        // mirror on xmin
        if (particle->m_location.x - m_Xmin <= distance)
        {
            m_sim->m_ghosts.current().m_location = mirrorPoint(&particle->m_location, &m_leftBorder);
            m_sim->m_ghosts.counterInc();

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in ReflectingRectangle::mirrorParticles:\n particle outside of the computational domain after mirroring: ( %f , %f ) -> ( %f , %f ) \n",
                        particle->m_location.x, particle->m_location.y, m_sim->m_ghosts.current().m_location.x,
                        m_sim->m_ghosts.current().m_location.y);
                exit(1);
            }
        }

        // mirror on xmax
        if (m_Xmax - particle->m_location.x <= distance)
        {
            m_sim->m_ghosts.current().m_location = mirrorPoint(&particle->m_location, &m_rightBorder);
            m_sim->m_ghosts.counterInc();

            if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
            {
                fprintf(stderr,
                        "ERROR in ReflectingRectangle::mirrorParticles:\n particle outside of the computational domain after mirroring: ( %f , %f ) -> ( %f , %f ) \n",
                        particle->m_location.x, particle->m_location.y, m_sim->m_ghosts.current().m_location.x,
                        m_sim->m_ghosts.current().m_location.y);
                exit(1);
            }
        }
    }
}

// mirror all particles at the walls
void ReflectingRectangle::mirrorAllParticles()
{
    unsigned int nofParticles = m_sim->m_particles.counter();

    if (nofParticles * 4 > m_sim->m_ghosts.max_size())
    {
        fprintf(stderr,
                "ERROR in ReflectingRectangle::mirrorAllParticles: too many particles, adapt global double factorGhostParticles in sim.cpp\n");
        exit(1);
    }

    for (unsigned int i = 0; i < nofParticles; i++)
    {
        // mirror on ymax
        m_sim->m_ghosts.current().m_location.x = m_sim->m_particles[i].m_location.x;
        m_sim->m_ghosts.current().m_location.y = -m_sim->m_particles[i].m_location.y + 2 * m_Ymax;
        if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
        {
            fprintf(stderr,
                    "ERROR in ReflectingRectangle::mirrorAllParticles:\n particle outside of the computational domain after mirroring: ( %f , %f ) -> ( %f , %f ) \n",
                    m_sim->m_particles[i].m_location.x, m_sim->m_particles[i].m_location.y,
                    m_sim->m_ghosts.current().m_location.x, m_sim->m_ghosts.current().m_location.y);
            exit(1);
        }
        m_sim->m_ghosts.counterInc();

        // mirror on ymin
        m_sim->m_ghosts.current().m_location.x = m_sim->m_particles[i].m_location.x;
        m_sim->m_ghosts.current().m_location.y = -m_sim->m_particles[i].m_location.y + 2 * m_Ymin;
        if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
        {
            fprintf(stderr,
                    "ERROR in ReflectingRectangle::mirrorAllParticles:\n particle outside of the computational domain after mirroring: ( %f , %f ) -> ( %f , %f ) \n",
                    m_sim->m_particles[i].m_location.x, m_sim->m_particles[i].m_location.y,
                    m_sim->m_ghosts.current().m_location.x, m_sim->m_ghosts.current().m_location.y);
            exit(1);
        }
        m_sim->m_ghosts.counterInc();

        // mirror on xmin
        m_sim->m_ghosts.current().m_location.x = -m_sim->m_particles[i].m_location.x + 2 * m_Xmin;
        m_sim->m_ghosts.current().m_location.y = m_sim->m_particles[i].m_location.y;
        if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
        {
            fprintf(stderr,
                    "ERROR in ReflectingRectangle::mirrorAllParticles:\n particle outside of the computational domain after mirroring: ( %f , %f ) -> ( %f , %f ) \n",
                    m_sim->m_particles[i].m_location.x, m_sim->m_particles[i].m_location.y,
                    m_sim->m_ghosts.current().m_location.x, m_sim->m_ghosts.current().m_location.y);
            exit(1);
        }
        m_sim->m_ghosts.counterInc();

        // mirror on xmax
        m_sim->m_ghosts.current().m_location.x = -m_sim->m_particles[i].m_location.x + 2 * m_Xmax;
        m_sim->m_ghosts.current().m_location.y = m_sim->m_particles[i].m_location.y;
        if (!m_sim->isInsideComputationalDomain(&(m_sim->m_ghosts.current().m_location)))
        {
            fprintf(stderr,
                    "ERROR in ReflectingRectangle::mirrorAllParticles:\n particle outside of the computational domain after mirroring: ( %f , %f ) -> ( %f , %f ) \n",
                    m_sim->m_particles[i].m_location.x, m_sim->m_particles[i].m_location.y,
                    m_sim->m_ghosts.current().m_location.x, m_sim->m_ghosts.current().m_location.y);
            exit(1);
        }
        m_sim->m_ghosts.counterInc();
    }
}

// intersect a line segment with the box
myArray<VoronoiParticle*>* ReflectingRectangle::intersectOnce(const Point& outerPoint, const Point& innerPoint)
{
    BOX_PRINTF("intersecting once:\n");
    BOX_PRINTF("\touter Point: ( %f , %f )\n", outerPoint.x, outerPoint.y);
    BOX_PRINTF("\tinnerPoint: ( %f , %f )\n", innerPoint.x, innerPoint.y);

    Point TRrelative = (*getTopRight()) - outerPoint;
    Point TLrelative = (*getTopLeft()) - outerPoint;
    Point BLrelative = (*getBottomLeft()) - outerPoint;
    Point BRrelative = (*getBottomRight()) - outerPoint;

    Point relative = innerPoint - outerPoint;

    double t = 0; // temporary t

    if (r2_geometry::isBetween(&TRrelative, &TLrelative, &relative)) // intersection with y = YMAX
    {
        if (relative.y == 0)
        {
            fprintf(stderr, "ERROR in ReflectingRectangle::intersectOnce: denominator is zero\n");
        }

        t = (m_Ymax - outerPoint.y) / relative.y;

        if (t <= 1) // line segment intersects
        {
            BOX_PRINTF("\t\tintersection with top border\n");
            return &m_topBorder;
        }
    }

    if (r2_geometry::isBetween(&TLrelative, &BLrelative, &relative)) // intersection with x = XMIN
    {
        if (relative.x == 0)
        {
            fprintf(stderr, "ERROR in ReflectingRectangle::intersectOnce: denominator is zero\n");
        }

        t = (m_Xmin - outerPoint.x) / relative.x;

        if (t <= 1) // line segment intersects
        {
            BOX_PRINTF("\t\tintersection with left border\n");
            return &m_leftBorder;
        }
    }

    if (r2_geometry::isBetween(&BLrelative, &BRrelative, &relative)) // intersection with y = YMIN
    {
        if (relative.y == 0)
        {
            fprintf(stderr, "ERROR in ReflectingRectangle::intersectOnce: denominator is zero\n");
        }

        t = (m_Ymin - outerPoint.y) / relative.y;

        if (t <= 1) // linesegment intersects
        {
            BOX_PRINTF("\t\tintersection with bottom border\n");
            return &m_bottomBorder;
        }
    }

    if (r2_geometry::isBetween(&BRrelative, &TRrelative, &relative)) // intersection with x = XMAX
    {
        if (relative.x == 0)
        {
            fprintf(stderr, "ERROR in ReflectingRectangle::intersectOnce: denominator is zero\n");
        }

        t = (m_Xmax - outerPoint.x) / relative.x;

        if (t <= 1) // line segment intersects
        {
            BOX_PRINTF("\t\tintersection with right border\n");
            return &m_rightBorder;
        }
    }

    fprintf(stderr, "ERROR in ReflectingRectangle::intersectOnce: no intersection\n");
    fprintf(stderr, "\touter Point: ( %f , %f )\n", outerPoint.x, outerPoint.y);
    fprintf(stderr, "\tinnerPoint: ( %f , %f )\n", innerPoint.x, innerPoint.y);
    exit(1);
    return NULL;
}

// initialization of m_particles, the graph is a cartesian grid
void ReflectingRectangle::sitesIniGrid(unsigned int particlesPerDimension)
{
    if (particlesPerDimension * particlesPerDimension > m_sim->getMaxNumberOfParticles())
    {
        fprintf(stderr,
                "ERROR in ReflectingRectangle::sitesIniGrid: too many particles, create a simulation with more maximum number of particles\n");
        exit(1);
    }

    // horizontal distance of the particles
    double dh = getDeltaX() / (particlesPerDimension);

    // vertical distance of the particles
    double dv = getDeltaY() / (particlesPerDimension);

    for (unsigned int j = 0; j < particlesPerDimension; j++)
    {
        for (unsigned int i = 0; i < particlesPerDimension; i++)
        {
            m_sim->m_particles.current().m_location.x = getXmin() + (i + 0.5) * dh;
            m_sim->m_particles.current().m_location.y = getYmin() + (j + 0.5) * dv;
            m_sim->m_particles.counterInc();
        }
    }
}

// initialization of m_particles, the graph is a cartesian grid
void ReflectingRectangle::sitesIniGrid(unsigned int particles_x_dimension, unsigned int particles_y_dimension)
{
    if (particles_x_dimension * particles_y_dimension > m_sim->getMaxNumberOfParticles())
    {
        fprintf(stderr,
                "ERROR in ReflectingRectangle::sitesIniGrid: too many particles, create a simulation with more maximum number of particles\n");
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

// initialization of m_particles, the graph is a cartesian grid but rotated by 45 degrees
void ReflectingRectangle::sitesIniGrid_diagonal(unsigned int particles_x_dimension, unsigned int particles_y_dimension)
{
    if (particles_x_dimension * particles_y_dimension > m_sim->getMaxNumberOfParticles())
    {
        fprintf(stderr,
                "ERROR in ReflectingRectangle::sitesIniGrid_diagonal: too many particles, create a simulation with more maximum number of particles\n");
        exit(1);
    }

    // horizontal distance of the particles
    double dh = getDeltaX() / (particles_x_dimension);

    // vertical distance of the particles
    double dv = getDeltaY() / (particles_y_dimension);

    Point location = Point(0, 0);
    Point pivot = Point(getXmin() + 0.5 * getDeltaX(), getYmin() + 0.5 * getDeltaY());

    for (unsigned int i = 0; i < 30 * particles_x_dimension; i++)
    {
        for (unsigned int j = 0; j < 30 * particles_y_dimension; j++)
        {
            location.x = getXmin() - 10 * getDeltaX() + (i + 0.5) * dh;
            location.y = getYmin() - 10 * getDeltaY() + (j + 0.5) * dv;

            location.rotate(M_PI / 4., &pivot);

            if (isInside(&location))
            {
                m_sim->m_particles.current().m_location.x = location.x;
                m_sim->m_particles.current().m_location.y = location.y;
                m_sim->m_particles.counterInc();
            }
        }
    }
}