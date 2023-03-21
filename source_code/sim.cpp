#include "config.h"
#include "debug.h"
#include "constants.h"
#include "sim.h"
#include "binary_search_tree.h"
#include "my_random.h"
#include "assert.h"
#include "black_red_tree.h"
#include <cmath>
#include "box_reflecting_rectangle.h"
#include "box_hybrid_rectangle.h"
#include "box_zero_grad_rectangle.h"

// create a box
void Sim::create_box()
{
    unsigned int type = Constants::getBOX_TYPE();
    double Xmin = Constants::getBOX_XMIN();
    double Ymin = Constants::getBOX_YMIN();
    double Xmax = Constants::getBOX_XMAX();
    double Ymax = Constants::getBOX_YMAX();

    BOX_PRINTF("\n---box---\n\n");

    BOX_PRINTF("x_min: %f\ny_min: %f\nx_max: %f\ny_max: %f\n\n", Xmin, Ymin, Xmax, Ymax);

    // computational domain (determines sweepline end position)
    m_cd_Xmax = 2 * Xmax - Xmin;
    m_cd_Ymax = 2 * Ymax - Ymin;

    m_cd_Xmin = 2 * Xmin - Xmax;
    m_cd_Ymin = 2 * Ymin - Ymax;

    BOX_PRINTF("computational domain:\n\nx_min: %f\ny_min: %f\nx_max: %f\ny_max: %f\n\n", m_cd_Xmin, m_cd_Ymin,
               m_cd_Xmax, m_cd_Ymax);

    double cd_deltaY = m_cd_Ymax - m_cd_Ymin;
    double cd_deltaX = m_cd_Xmax - m_cd_Xmin;

    m_sweepY = m_cd_Ymin - cd_deltaY - cd_deltaX;

    BOX_PRINTF("final sweepline position: %f\n", m_sweepY);

    switch (type)
    {
        case 0: // reflecting rectangle
        {
            m_box = new ReflectingRectangle(this, Xmin, Xmax, Ymin, Ymax);
        }
            break;

        case 1: // hybrid rectangle
        {
            m_box = new HybridRectangle(this, Xmin, Xmax, Ymin, Ymax);
        }
            break;

        case 2: // zero_grad rectangle
        {
            m_box = new ZeroGradRectangle(this, Xmin, Xmax, Ymin, Ymax);
        }
            break;

        default:
        {
            fprintf(stderr, "ERROR in Sim::createBox(): invalid box type\n");
            exit(1);
        }
    }
}

// create an obstacle
void Sim::create_obstacle()
{
    // disk
    m_obstacle = new Obstacle(this);
    m_obstacle->disk();

    // assign initial velocities
    m_obstacle->m_velocity.x = Constants::getOBS_VX();
    m_obstacle->m_velocity.y = Constants::getOBS_VY();
}

// check whether a point is inside the computational domain
bool Sim::isInsideComputationalDomain(Point* point)
{
    return (point->x <= m_cd_Xmax && point->x >= m_cd_Xmin && point->y <= m_cd_Ymax && point->y >= m_cd_Ymin);
}

// storage management parameters
//

// more half edges and vertices are needed for the construction than for the final graph
static double ratio = 2.;

// maxNumberOfParticles * factorGhostParticles = number of ghost particles
double factorGhostParticles = 1;

// number of obstacle particles + amr particles
unsigned int nof_amr_particles = 100000;
unsigned int nof_obstacle_particles = 10000;

Sim::Sim(unsigned int maxNumberOfParticles)
        :
        m_maxNumberOfParticles(maxNumberOfParticles),
        m_cd_Xmin(0),
        m_cd_Ymin(0),
        m_cd_Xmax(0),
        m_cd_Ymax(0),
        m_sweepY(0),
        m_box(NULL),
        m_obstacle(NULL),
        m_particles(maxNumberOfParticles),
        m_ghosts(factorGhostParticles * maxNumberOfParticles + 10000),
        m_vertices((int) ratio * (2 * (maxNumberOfParticles + factorGhostParticles * maxNumberOfParticles) + 2)),
        m_halfEdges((int) ratio * (6 * (maxNumberOfParticles + factorGhostParticles * maxNumberOfParticles) + 2)),
        m_BR_root(BR_Node::p_nil),
        m_root(NULL) // max: 2n-1 arcs, 2n-2 breakpoints -> maximum number of elements < 4n
{
    m_breakpointVertices = std::vector<Vertex*>();

    // reserve memory
    m_breakpointVertices.reserve(maxNumberOfParticles + factorGhostParticles * maxNumberOfParticles);
}

Sim::~Sim()
{
    // delete the content of m_vertices
    // delete the content of m_particles
    // delete the content of m_halfEdges
    // delete the content of m_boxSegments
    // => not needed

    // delete the box pointer
    delete m_box;

    // delete the obstacle pointer
    if (m_obstacle != NULL)
    {
        delete m_obstacle;
    }
}

// deactivate a vertex 
void Sim::deleteVertex(Vertex* vertex)
{
    vertex->m_inUse = false;
    m_vertices.deactivatedInc();
}

// deactivate a half edge
void Sim::deleteHalfEdge(HalfEdge* halfEdge)
{
    halfEdge->m_inUse = false;
    m_halfEdges.deactivatedInc();
}

// clears the graph and sweepline data between two time steps
void Sim::freeGraph()
{
    // reset the Voronoi particles
    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        m_particles[i].reset();
    }

    for (unsigned int k = 0; k < m_ghosts.counter(); k++)
    {
        m_ghosts[k].reset();
    }

    // reset counters
    m_ghosts.reset();
    m_vertices.reset();
    m_halfEdges.reset();

    // delete the content of the event queue
    BR_Node::freeTree(m_BR_root);
    m_BR_root = BR_Node::p_nil;

    // delete the tree. note: the whole tree is deleted recursively
    tree::freeTree(m_root);
    m_root = NULL;

    // delete content of m_breakpointVertices
    for (unsigned int j = 0; j < m_breakpointVertices.size(); j++)
    {
        delete m_breakpointVertices[j];
    }
    m_breakpointVertices.clear();

    // the following queues should be empty
    assert(m_breakpointVertices.empty());

    // don't reset the box here! data from previous time step is needed!
}

// initialization of the primitive variables, one constant state
void Sim::primVarIni1State(double rho, double v_x, double v_y, double pressure)
{
    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        particle = &m_particles[i];

        // set primitive variables
        particle->m_velocity = Vector(v_x, v_y);
        particle->m_density = rho;
        particle->m_pressure = pressure;
    }
}

// initialization of the primitive variables, two constants states left and right of border_x
void Sim::primVarIni2States(double rho_left, double v_x_left, double p_left, double border_x, double rho_right,
                            double v_x_right, double p_right)
{
    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        particle = &m_particles[i];

        if (particle->m_location.x < border_x)
        {
            // set primitive variables
            particle->m_velocity = Vector(v_x_left, 0);
            particle->m_density = rho_left;
            particle->m_pressure = p_left;
        }
        else
        {
            // set primitive variables
            particle->m_velocity = Vector(v_x_right, 0);
            particle->m_density = rho_right;
            particle->m_pressure = p_right;
        }
    }
}

// initialization of the primitive variables, two constant states, top and bottom of border_y
void Sim::primVarIni2States_vertical(double rho_top, double v_y_top, double p_top, double border_y, double rho_bottom,
                                     double v_y_bottom, double p_bottom)
{
    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        particle = &m_particles[i];

        if (particle->m_location.y > border_y)
        {
            // set primitive variables
            particle->m_velocity = Vector(0, v_y_top);
            particle->m_density = rho_top;
            particle->m_pressure = p_top;
        }
        else
        {
            // set primitive variables
            particle->m_velocity = Vector(0, v_y_bottom);
            particle->m_density = rho_bottom;
            particle->m_pressure = p_bottom;
        }
    }
}

// initialization of the primitive variables, three states left middle and right,
// borders: border_lm between left state and middle and border_mr between middle and right state
void Sim::primVarIni3States(double rho_left, double v_x_left, double p_left, double border_lm, double rho_middle,
                            double v_x_middle, double p_middle, double border_mr, double rho_right, double v_x_right,
                            double p_right)
{
    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        particle = &m_particles[i];

        if (particle->m_location.x < border_lm)
        {
            // set primitive variables
            particle->m_velocity = Vector(v_x_left, 0);
            particle->m_density = rho_left;
            particle->m_pressure = p_left;
        }
        else if (particle->m_location.x >= border_lm && particle->m_location.x <= border_mr)
        {
            // set primitive variables
            particle->m_velocity = Vector(v_x_middle, 0);
            particle->m_density = rho_middle;
            particle->m_pressure = p_middle;
        }
        else
        {
            // set primitive variables
            particle->m_velocity = Vector(v_x_right, 0);
            particle->m_density = rho_right;
            particle->m_pressure = p_right;
        }
    }
}

// initialization of the primitive variables
void Sim::primVarIni3States_vertical(double rho_top, double v_x_top, double p_top, double border_top, double rho_middle,
                                     double v_x_middle, double p_middle, double border_bottom, double rho_bottom,
                                     double v_x_bottom, double p_bottom)
{
    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        particle = &m_particles[i];

        if (particle->m_location.y > border_top)
        {
            // set primitive variables
            particle->m_velocity = Vector(v_x_top, 0);
            particle->m_density = rho_top;
            particle->m_pressure = p_top;
        }
        else if (particle->m_location.y < border_bottom)
        {
            // set primitive variables
            particle->m_velocity = Vector(v_x_bottom, 0);
            particle->m_density = rho_bottom;
            particle->m_pressure = p_bottom;
        }
        else
        {
            // set primitive variables
            particle->m_velocity = Vector(v_x_middle, 0);
            particle->m_density = rho_middle;
            particle->m_pressure = p_middle;
        }
    }
}

// inject energy
void Sim::injectEnergy(double x_min, double y_min, double x_max, double y_max, double energy)
{
    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        particle = &m_particles[i];

        if (particle->m_location.x >= x_min && particle->m_location.x <= x_max && particle->m_location.y >= y_min &&
            particle->m_location.y <= y_max)
        {
            particle->m_energy = energy;

            return; // return after one cell gets the energy
        }
    }

    fprintf(stderr, "ERROR in Sim::injectEnergy: no particle found inside the given domain, no energy injected\n");
    exit(1);
}

// inject energy into one particle
void Sim::injectEnergyParticle(int i, double energy)
{
    m_particles.at(i).m_energy = energy;
}

// inject energy
void Sim::injectEnergyMiddleParticle(double energy)
{
    int nof_particles = m_particles.counter();

    if (nof_particles % 2 == 0)
    {
        fprintf(stderr,
                "ERROR in Sim::injectEnergyMiddleParticle: number of particles is even, odd number is required\n");
        exit(1);
    }

    m_particles[(int) (nof_particles * 0.5)].m_energy = energy;
}

// inject energy
void Sim::injectEnergyGauss(double energy, double radius, double x, double y)
{
    VoronoiParticle* particle = NULL;
    double r = 0;
    double sigma = radius / 4.;

    Point relative_location = Point(0, 0);
    Point center = Point(x, y);

    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        particle = &m_particles[i];

        relative_location = center * (-1.) + particle->m_location;

        r = sqrt(relative_location.length_square()); // distance to center

        if (r < radius)
        {
            particle->m_energy = 1. / (sigma * sqrt(2 * M_PI)) * exp(-0.5 * pow(r / sigma, 2));
        }
    }

    // sum up energy
    double energy_sum = 0;

    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        particle = &m_particles[i];

        relative_location = center * (-1.) + particle->m_location;

        r = sqrt(relative_location.length_square());

        if (r < radius)
        {
            energy_sum += particle->m_energy;
        }
    }

    // normalize energy to 1
    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        particle = &m_particles[i];

        relative_location = center * (-1.) + particle->m_location;

        r = sqrt(relative_location.length_square());

        if (r < radius)
        {
            particle->m_energy /= energy_sum;
        }
    }
}

// initialization of m_particles
void Sim::sitesIniFile(const char* filename)
{
    FILE* pFile;

    pFile = fopen(filename, "r");

    if (pFile == NULL)
    {
        fprintf(stderr, "ERROR in sitesIniFile: can't open File %s\n", filename);
        exit(1);
    }

    double x = 0;
    double y = 0;

    while (true)
    {
        if (fscanf(pFile, "%lf %lf", &x, &y) == EOF)
        {
            break;
        }

        m_particles.current().m_location.x = x;
        m_particles.current().m_location.y = y;

        // check whether the particles are inside the box
        if (!(m_box->isInside(&(m_particles.current().m_location))))
        {
            fprintf(stderr, "ERROR in Sim::sitesIniFile: site is outside of the box\n");

            fclose(pFile);

            exit(1);
        }

        m_particles.counterInc();
    }

    fclose(pFile);
}

// initialization of m_particles
void Sim::sitesIniRandomDouble(double lowerBound_x, double lowerBound_y, double upperBound_x, double upperBound_y,
                               unsigned int nofSites)
{
    // check whether the particles are inside the box
    Point BL(lowerBound_x, lowerBound_y); // bottom left point
    Point TR(upperBound_x, upperBound_y); // top right point

    if (!(m_box->isInside(&BL) && m_box->isInside(&TR)))
    {
        fprintf(stderr, "ERROR in Sim::sitesIniRandom: sites could be outside of the box\n");
        exit(1);
    }

    for (unsigned int i = 0; i < nofSites; i++)
    {
        m_particles.current().m_location.x = myRandom::randomDouble(lowerBound_x, upperBound_x);
        m_particles.current().m_location.y = myRandom::randomDouble(lowerBound_y, upperBound_y);
        m_particles.counterInc();
    }
}

// initialization of m_particles
void Sim::sitesIniRandomInteger(double lowerBound_x, double lowerBound_y, double upperBound_x, double upperBound_y,
                                unsigned int nofSites)
{
    // check whether the particles are inside the box
    Point BL(lowerBound_x, lowerBound_y); // bottom left point
    Point TR(upperBound_x, upperBound_y); // top right point

    if (!(m_box->isInside(&BL) && m_box->isInside(&TR)))
    {
        fprintf(stderr, "ERROR in Sim::sitesIniRandom: sites could be outside of the box\n");
        exit(1);
    }

    for (unsigned int i = 0; i < nofSites; i++)
    {
        m_particles.current().m_location.x = myRandom::randomInteger(lowerBound_x, upperBound_x);
        m_particles.current().m_location.y = myRandom::randomInteger(lowerBound_y, upperBound_y);
        m_particles.counterInc();
    }
}


// initialization of m_particles in a circle
void Sim::sitesIniRandomDoubleCircle(unsigned int nofSites, double radius)
{
    double x = 0;
    double y = 0;
    Point p = Point(0, 0);

    for (unsigned int i = 0; i < nofSites; i++)
    {
        x = myRandom::randomDouble(-radius, radius);
        y = myRandom::randomDouble(-radius, radius);

        p.x = x;
        p.y = y;

        if (!(x * x + y * y <= radius * radius)) // the point has to be inside the circle
        {
            i--;
            continue;
        }

        m_particles.current().m_location = p;
        m_particles.counterInc();
    }
}

void Sim::sitesIniPolarGrid(unsigned int particlesInRhoDimension, unsigned int particlesInPhiDimension, double radius)
{
    if (particlesInRhoDimension * particlesInPhiDimension + 1 > m_maxNumberOfParticles)
    {
        fprintf(stderr, "ERROR in Sim::sitesIniCircleGrid: too many particles\n");
        exit(1);
    }

    // set one particle in the center
    m_particles.current().m_location.x = 0;
    m_particles.current().m_location.y = 0;
    m_particles.counterInc();

    // radial distance of the particles
    double dr = radius / (particlesInRhoDimension + 1);

    // angular distance of the particles in circular measure
    double dphi = 2 * M_PI / particlesInPhiDimension;

    for (unsigned int i = 0; i < particlesInRhoDimension; i++)
    {
        for (unsigned int j = 0; j < particlesInPhiDimension; j++)
        {
            m_particles.current().m_location.x = (i + 1) * dr * cos(j * dphi);
            m_particles.current().m_location.y = (i + 1) * dr * sin(j * dphi);

            if (this->m_box->isInside(&m_particles.current().m_location))
            {
                m_particles.counterInc();
            }
        }
    }
}

void Sim::sitesIniGridinCircle(unsigned int particlesPerDiameter, double radius)
{
    if (particlesPerDiameter * particlesPerDiameter > m_maxNumberOfParticles)
    {
        fprintf(stderr, "ERROR in Sim::sitesIniGrid: too many particles, adapt the constant MAX_NOF_PARTICLES\n");
        exit(1);
    }

    // horizontal distance of the particles
    double dh = 2 * radius / (particlesPerDiameter);

    // vertical distance of the particles
    double dv = 2 * radius / (particlesPerDiameter);

    // position of a particle
    Point p = Point(0, 0);

    for (unsigned int i = 0; i < particlesPerDiameter; i++)
    {
        for (unsigned int j = 0; j < particlesPerDiameter; j++)
        {
            p.x = -radius + (i + 0.5) * dh;
            p.y = -radius + (j + 0.5) * dv;

            if (p.x * p.x + p.y * p.y <= radius * radius) // the point has to be inside the circle
            {
                m_particles.current().m_location.x = -radius + (i + 0.5) * dh;
                m_particles.current().m_location.y = -radius + (j + 0.5) * dv;
                m_particles.counterInc();
            }
        }
    }
}

// initialization of the event queue
void Sim::eventQueueIni()
{
    if (m_particles.counter() < 0)
    {
        fprintf(stderr, "ERROR in Sim::eventQueueIni: less than two Voronoi particles\n");
        exit(1); // initialize the Voronoi particles first!
    }

    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        assert(m_box->isInside(&(m_particles[i].m_location))); // the particles have to be inside the box

        Event* p_siteEvent = new SiteEvent(&(m_particles[i]));
        BR_Node::insertSite(p_siteEvent, this);
    }
}

// initializes the member m_eventQueue with ghost particles
void Sim::eventQueueIniGhosts()
{
    for (unsigned int i = 0; i < m_ghosts.counter(); i++)
    {
        assert(isInsideComputationalDomain(
                &(m_ghosts[i].m_location))); // the particles have to be inside the computational domain

        Event* p_siteEvent = new SiteEvent(&(m_ghosts[i]));
        BR_Node::insertSite(p_siteEvent, this);
    }
}

// print events
void Sim::insertPrintEvents(double yCoords[], unsigned int nofElements)
{
    for (unsigned int i = 0; i < nofElements; i++)
    {
        Event* event = new PrintEvent(0, yCoords[i]);
        BR_Node::insert(event, this);
    }
}

// print the cell center of volume of each cell
void Sim::printCellCentersOfVolume(const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int i = 0; i < m_particles.counter(); i++)
    {
        fprintf(pFile, "%f %f\n", m_particles[i].m_cellCenterOfVolume.x, m_particles[i].m_cellCenterOfVolume.y);
    }

    fclose(pFile);
}

// print the cell center of volume of each ghost cell
void Sim::printGhostCellCentersOfVolume(const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int i = 0; i < m_ghosts.counter(); i++)
    {
        fprintf(pFile, "%f %f\n", m_ghosts[i].m_cellCenterOfVolume.x, m_ghosts[i].m_cellCenterOfVolume.y);
    }

    fclose(pFile);
}