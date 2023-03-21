#ifndef BOX_REFLECTING_RECTANGLE_H
#define BOX_REFLECTING_RECTANGLE_H

#include "r2.h"
#include "box.h"
#include "my_array.h"

class Sim;

/// a reflecting rectangle
class ReflectingRectangle : public Box
{

private:

    // box parameters
    double m_Xmin;
    double m_Xmax;
    double m_Ymin;
    double m_Ymax;

    // additional box parameters
    double m_deltaX;
    double m_deltaY;
    Point m_topRight;
    Point m_topLeft;
    Point m_bottomLeft;
    Point m_bottomRight;
    Point m_middle;

    // arrays storing border particles to copy
    myArray<VoronoiParticle*> m_rightBorder;
    myArray<VoronoiParticle*> m_topBorder;
    myArray<VoronoiParticle*> m_leftBorder;
    myArray<VoronoiParticle*> m_bottomBorder;

public:

    // getters and setters
    inline Point* getTopRight()
    {
        return &m_topRight;
    }

    inline Point* getTopLeft()
    {
        return &m_topLeft;
    }

    inline Point* getBottomLeft()
    {
        return &m_bottomLeft;
    }

    inline Point* getBottomRight()
    {
        return &m_bottomRight;
    }

    inline Point* getMiddle()
    {
        return &m_middle;
    }

    inline double getXmin()
    {
        return m_Xmin;
    }

    inline double getXmax()
    {
        return m_Xmax;
    }

    inline double getYmin()
    {
        return m_Ymin;
    }

    inline double getYmax()
    {
        return m_Ymax;
    }

    inline double getDeltaX()
    {
        return m_deltaX;
    }

    inline double getDeltaY()
    {
        return m_deltaY;
    }

    /// constructor
    ReflectingRectangle(Sim* sim, double Xmin, double Xmax, double Ymin, double Ymax);

    /// destructor
    ~ReflectingRectangle();

    /// reset the arrays
    void reset();

    /// check whether a point is inside of the box
    bool isInside(Point* point);

    /// search for all neighbours of the particle, which are inside the box
    void searchNeighboursInside(VoronoiParticle* particle);

    /// handle two neighbours
    void handleNeighbours(VoronoiParticle* particle_outside, VoronoiParticle* particle_inside);

    /// handle the neighbours of the neighbours
    void handleNeighboursOfNeighbours(VoronoiParticle* particle_inside, myArray<VoronoiParticle*>* intersection);

    /// determine ghost particles
    void determineGhostParticles();

    /// create ghost particles, previous determined via determineGhostParticles
    void createGhostParticles();

    /// copy the primitive variables to the ghost particles
    void copyPrimitiveVariablesCellCentersAndVolume();

    /// copy the gradients to the ghost particles
    void copyGradients();

    /// assign velocities to the ghost particles
    void copyParticleVelocities();

    /// calculate the ghost gradient in x-direction
    Vector ghost_gradient_x(Vector* real_gradient, myArray<VoronoiParticle*>* segment);

    /// calculate the ghost gradient in y-direction
    Vector ghost_gradient_y(Vector* real_gradient, myArray<VoronoiParticle*>* segment);

    /// mirror a vector on a box segment
    Vector mirrorVector(Vector* vector, myArray<VoronoiParticle*>* segment);

    /// mirror a point on a box segment
    Point mirrorPoint(Point* point, myArray<VoronoiParticle*>* segment);

    /// mirror particles with distance less equal "distance" to a wall
    void create_ghost_particles(double distance);

    /// mirror all particles at the walls
    void mirrorAllParticles();

    /// move a particle which is outside of the boundary
    void handle_outside_boundary(VoronoiParticle* vp, unsigned int& iterator)
    {}; // not needed here

    /// move the ghost particles
    void move_ghost_particles()
    {}; // not needed here

    /**
    intersect a linesegment with the box
    \param outerPoint the location of the point outside of the box
    \param innerPoint the location of the point inside of the box
    \return address of the intersected border or NULL
    */
    myArray<VoronoiParticle*>* intersectOnce(const Point& outerPoint, const Point& innerPoint);

    /**
    intersect a linesegment with the box
    \param intersection if an intersection exists its stored here
    \param outerPoint the location of the point outside of the box
    \param secondPoint the other point, it can be outside or inside or on the rim of the box
    \return does an intersection exist?
    */

    /// initialization of m_particles, the graph is a cartesian grid
    void sitesIniGrid(unsigned int particlesPerDimension);

    /// initialization of m_particles, the graph is a cartesian grid
    void sitesIniGrid(unsigned int particles_x_dimension, unsigned int particles_y_dimension);

    /// initialization of m_particles, the graph is a cartesian grid but rotated by 45 degrees
    void sitesIniGrid_diagonal(unsigned int particles_x_dimension, unsigned int particles_y_dimension);
};

#endif