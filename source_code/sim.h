#ifndef SIM_H
#define SIM_H

#include <vector>
#include "voronoi_particle.h"
#include "event.h"
#include "my_array.h"
#include "box.h"
#include "box_reflecting_rectangle.h"
#include "voronoi_ghost_particle.h"
#include "obstacle.h"

class Node;

class Vertex;

// maxNumberOfParticles * factorGhostParticles + 10000 = number of ghost particles
extern double factorGhostParticles;

extern unsigned int nof_obstacle_particles;
extern unsigned int nof_amr_particles;

/**
Simulation class, stores all important data and provides initialization and output functions
*/
class Sim
{
private:
    unsigned int m_maxNumberOfParticles; ///< maximum number of particles

    // computational domain
    double m_cd_Xmin;
    double m_cd_Ymin;
    double m_cd_Xmax;
    double m_cd_Ymax;

    double m_sweepY; ///< algorithm ends at this sweepline position

public:

    // getters and setters
    inline unsigned int getMaxNumberOfParticles()
    {
        return m_maxNumberOfParticles;
    }

    inline void setMaxNumberOfParticles(unsigned int maxNumberOfParticles)
    {
        m_maxNumberOfParticles = maxNumberOfParticles;
    }

    inline double getSweepY()
    {
        return m_sweepY;
    }

    inline double getCDXmin()
    {
        return m_cd_Xmin;
    }

    inline double getCDYmin()
    {
        return m_cd_Ymin;
    }

    inline double getCDXmax()
    {
        return m_cd_Xmax;
    }

    inline double getCDYmax()
    {
        return m_cd_Ymax;
    }

    /// the boundary
    Box* m_box;

    /// an obstacle
    Obstacle* m_obstacle;

    /// create a box
    void create_box();

    /// create an obstacle
    void create_obstacle();

    /// check whether a point is inside the computational domain
    bool isInsideComputationalDomain(Point* point);

    // general data
    myArray<VoronoiParticle> m_particles; ///< an array with all particles
    myArray<GhostParticle> m_ghosts; ///< an array with ghost particles
    myArray<Vertex> m_vertices; ///< an array with all vertices of the Voronoi graph
    myArray<HalfEdge> m_halfEdges; ///< an array with all half edges of the Voronoi graph

    // data for the sweepline algorithm
    BR_Node* m_BR_root; ///< a pointer to the root of the red black tree (event queue)
    Node* m_root;    ///< a pointer to a node representing a binary search tree
    std::vector<Vertex*> m_breakpointVertices; ///< a vector with all vertices of the graph which have less then three edges

    Sim(unsigned int maxNumberOfParticles); ///< constructor
    ~Sim(); ///< destructor

    void deleteVertex(Vertex* vertex); ///< deactivate a vertex
    void deleteHalfEdge(HalfEdge* halfEdge); ///< deactivate a half edge

    void freeGraph(); ///< clears the graph and sweepline data between two timesteps

    // initialization of the primitive variables
    //

    /// initialization of the primitive variables, one constant state
    void primVarIni1State(double rho, double v_x, double v_y, double pressure);

    /// initialization of the primitive variables, two constants states left and right of border_x
    void primVarIni2States(double rho_left, double v_x_left, double p_left, double border_x, double rho_right,
                           double v_x_right, double p_right);

    /// initialization of the primitive variables, two constant states, top and bottom of border_y
    void primVarIni2States_vertical(double rho_top, double v_y_top, double p_top, double border_y, double rho_bottom,
                                    double v_y_bottom, double p_bottom);

    /// initialization of the primitive variables, three states left middle and right,
    /// borders: border_lm between left state and middle and border_mr between middle and right state
    void primVarIni3States(double rho_left, double v_x_left, double p_left, double border_lm, double rho_middle,
                           double v_x_middle, double p_middle, double border_mr, double rho_right, double v_x_right,
                           double p_right);

    /// initialization of the primitive variables
    void primVarIni3States_vertical(double rho_top, double v_x_top, double p_top, double border_top, double rho_middle,
                                    double v_x_middle, double p_middle, double border_bottom, double rho_bottom,
                                    double v_x_bottom, double p_bottom);

    /// inject energy
    void injectEnergy(double x_min, double y_min, double x_max, double y_max, double energy);

    void injectEnergyParticle(int i, double energy);

    void injectEnergyMiddleParticle(double energy);

    void injectEnergyGauss(double energy, double radius, double x, double y);

    // initialization of the sites
    //

    /// initializes the member m_particles and sets m_location of every particle
    void sitesIniFile(const char* filename);

    /// initializes the member m_particles and sets m_location of every particle
    void sitesIniRandomDouble(double lowerBound_x, double lowerBound_y, double upperBound_x, double upperBound_y,
                              unsigned int nofSites);

    /// initializes the member m_particles and sets m_location of every particle
    void sitesIniRandomInteger(double lowerBound_x, double lowerBound_y, double upperBound_x, double upperBound_y,
                               unsigned int nofSites);

    /// random initialization of m_particles in a Circle
    void sitesIniRandomDoubleCircle(unsigned int nofSites, double radius);

    /// initialization of m_particles in polar arrangement
    void sitesIniPolarGrid(unsigned int particlesInRhoDimension, unsigned int particlesInPhiDimension, double radius);

    /// initialization of m_particles in cartesian arrangement inside of a circle
    void sitesIniGridinCircle(unsigned int particlesPerDimension, double radius);

    // ini methods for the sweep line algorithm
    //

    void eventQueueIni(); ///< initializes the member m_eventQueue with particles
    void eventQueueIniGhosts(); ///< initializes the member m_eventQueue with ghost particles

    /**
    Inserts a print event into the priority queue m_eventQueue.
    \param yCoords an array with doubles representing the y coordinate of the print event
    \param nofElements the number of elements of the array (first parameter)
    \sa breakp_x()
    */
    void insertPrintEvents(double yCoords[], unsigned int nofElements); ///< Inserts print Events into the queue

    // m_particles
    void printCellCentersOfVolume(const char* filename); ///< print the cell center of volume of each cell
    void printGhostCellCentersOfVolume(const char* filename); ///< print the cell center of volume of each ghost cell
};

#endif