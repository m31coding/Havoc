#ifndef OBSTACLE_H
#define OBSTACLE_H

#include "my_array.h"
#include "r2.h"

class VoronoiParticle;

class Sim;

class Obstacle
{
public:

    Sim* m_sim;

    /**
    includes all the obstacle particles
    order: inner particle, corresponding outer particle, inner particle, corresponding outer particle, a.s.o
    */
    myArray<VoronoiParticle*> m_obstacle_particles;

    unsigned int m_nof_obstacle_particles;

    Obstacle(Sim* sim);

    ~Obstacle();

    /// the current velocity of the obstacle
    Point m_velocity;

    /// binary input
    void bin_input(unsigned int nof_particles_in_use);

    /// create a disk obstacle
    void disk();

    /// delta r
    double delta_r();

    /// delta phi
    double delta_phi();

    /// position
    Point obstacle_center();

    /// is point inside of the obstacle?
    bool is_inside(Point* location);

    /// assign velocities to the obstacle particles
    void assign_velocities();

    /// delete zero length edges
    void delete_zero_length_edges();

    /// minimum distance to the obstacle
    double minimum_distance_square();

    /// handle a particle inside the obstacle
    void handle_inside_obstacle(VoronoiParticle* vp);

    /// copy primitive variables and cell centers to the inner obstacle cells
    void copyPrimitiveVariablesCellCentersAndVolume();

    /// copy gradients to the inner obstacle cells
    void copyGradients();

    /// calculate the force of the fluid
    Vector calculate_force();
};

#endif