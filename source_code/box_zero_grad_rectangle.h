#ifndef BOX_ZGRADIENT_RECTANGLE_H
#define BOX_ZGRADIENT_RECTANGLE_H

#include "box.h"
#include "r2.h"

class GhostParticle;

struct zgr_ghost
{
    Point location;
    int flag;
};

class ZeroGradRectangle : public Box
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

    // ghost box
    double m_ghost_box_x_min;
    double m_ghost_box_x_max;
    double m_ghost_box_y_min;
    double m_ghost_box_y_max;

    // locations of the ghost particles
    zgr_ghost* m_ghost_points;
    unsigned int m_nof_ghost_points;

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

    // constructor
    ZeroGradRectangle(Sim* sim, double Xmin, double Xmax, double Ymin, double Ymax);

    // destructor
    ~ZeroGradRectangle();

    /// reset function
    void reset();

    /// check whether a point is inside of the box
    bool isInside(Point* point);

    /// check whether a point is inside of the ghost box
    bool is_inside_ghost_box(Point* point);

    // check the fixed ghosts
    bool is_fixed_left(Point* point);

    bool is_fixed_right(Point* point);

    bool is_fixed_top(Point* point);

    bool is_fixed_bottom(Point* point);

    /// determine ghost particles
    void determineGhostParticles();

    /// check the volume of a fixed cell and insert a new particle if needed
    void potential_particle_insertion(GhostParticle* ghost);

    /// create ghost particles out of particles close to the border
    void create_ghost_particles(double distance);

    /// create ghost particles
    void createGhostParticles();

    /// create zero gradient ghost particles
    void create_ghost_particles_zero_grad();

    /// copy the primitive variables to the ghost particles
    void copyPrimitiveVariablesCellCentersAndVolume();

    /// copy the gradients to the ghost particles
    void copyGradients();

    /// assign velocities to the ghost particles
    void copyParticleVelocities();

    /// move a particle which is outside of the boundary
    void handle_outside_boundary(VoronoiParticle* vp, unsigned int& iterator);

    /// move the ghost particles
    void move_ghost_particles();

    /// ini function
    void sitesIniGrid(unsigned int particles_x_dimension, unsigned int particles_y_dimension);
};

#endif