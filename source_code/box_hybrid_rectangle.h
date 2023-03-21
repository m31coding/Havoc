#ifndef BOX_HYBRID_RECTANGLE_H
#define BOX_HYBRID_RECTANGLE_H

#include "r2.h"
#include "box.h"
#include "voronoi_particle.h"
#include "sim.h"

/// a box with periodic boundaries left-right and reflecting boundaries top and bottom
class HybridRectangle : public Box
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

    // flags
    myArray<VoronoiParticle*> m_top_right_border;
    myArray<VoronoiParticle*> m_top_left_border;
    myArray<VoronoiParticle*> m_bottom_right_border;
    myArray<VoronoiParticle*> m_bottom_left_border;

    // distance for ghost particles
    double m_distance_x;
    double m_distance_y;

public:

    // getters and setters
    inline double get_distance_x()
    {
        return m_distance_x;
    }

    inline double get_distance_y()
    {
        return m_distance_y;
    }

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
    HybridRectangle(Sim* sim, double Xmin, double Xmax, double Ymin, double Ymax);

    /// destructor
    ~HybridRectangle();

    /// reset the arrays
    void reset();

    /// check whether a point is inside of the box
    bool isInside(Point* point);

    /// handle a point
    Point handle_point(Point* point, myArray<VoronoiParticle*>* segment);

    /// handle a vector
    Vector handle_vector(Vector* vector, myArray<VoronoiParticle*>* segment);

    /// create ghost particles out of particles close to the border
    void create_ghost_particles(double distance_x, double distance_y);

    void create_ghost_particles(double distance);

    /// determine ghost particles
    void determineGhostParticles();

    /// create ghost particles, previous determined via determineGhostParticles
    void createGhostParticles();

    /// copy the primitive variables to the ghost particles
    void copyPrimitiveVariablesCellCentersAndVolume();

    /// copy the gradients to the ghost particles
    void copyGradients();

    /// calculate the ghost gradient in x-direction
    Vector ghost_gradient_x(Vector* real_gradient, myArray<VoronoiParticle*>* segment);

    /// calculate the ghost gradient in x-direction
    Vector ghost_gradient_y(Vector* real_gradient, myArray<VoronoiParticle*>* segment);

    /// assign velocities to the ghost particles
    void copyParticleVelocities();

    /// move a particle which is outside of the boundary
    void handle_outside_boundary(VoronoiParticle* vp, unsigned int& iterator);

    /// move the ghost particles
    void move_ghost_particles()
    {}; // not needed here

    /// ini functions
    void sitesIniGrid(unsigned int particles_x_dimension, unsigned int particles_y_dimension);
};

#endif