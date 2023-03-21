#include "config.h"
#include "binary_input_output.h"
#include "debug.h"
#include "sim_core.h"
#include "string.h"
#include "options.h"
#include <limits>

// 100 bytes header in output file
static unsigned int headersize = 200;

// output version:
static unsigned int version = 8;

#define _3_(...) if(outputVersion >= 3)__VA_ARGS__
#define _2_(...) if(outputVersion >= 2)__VA_ARGS__
#define _4_(...) if(outputVersion >= 4)__VA_ARGS__
#define _5_(...) if(outputVersion >= 5)__VA_ARGS__
#define _6_(...) if(outputVersion >= 6)__VA_ARGS__
#define _7_(...) if(outputVersion >= 7)__VA_ARGS__
#define _8_(...) if(outputVersion >= 7)__VA_ARGS__

// write most important data into a binary file
void bin_io::binaryOutput(const char* file, unsigned int filenumber, Sim* sim)
{
    size_t pointSize = sizeof(Point);
    size_t doubleSize = sizeof(double);
    size_t int_size = sizeof(int);

    char filename[9];
    sprintf(filename, "%05u.hc", filenumber);

    //open file
    FILE* pFile;

    if (file == NULL)
    {
        printf(" (binary output %u)", filenumber);
        IO_PRINTF("\n\n---binary output of file %u---\n", filenumber);

        pFile = fopen(filename, "w");
    }
    else
    {
        pFile = fopen(file, "w");
    }

    if (pFile == NULL)
    {
        fprintf(stderr, "ERROR in Sim::binaryOutput: cant't open file %s\n", filename);
        exit(1);
    }

    // header
    //

    char fileEnding[] = "hc"; // 3 bytes (hc + terminating character)
    unsigned int outputVersion = version; // 4 bytes
    double time = sim_core::simulation_time; // 8 bytes
    unsigned int output_counter = sim_core::counter; // 8 bytes
    unsigned int particlesCounter = sim->m_particles.counter(); // 4 bytes
    unsigned int ghostsCounter = sim->m_ghosts.counter(); // 4 bytes
    unsigned int nof_obstacle_particles = 0; // 4 bytes
    double obstacle_v_x = 0; // 8 bytes
    double obstacle_v_y = 0; // 8 bytes
    unsigned int numberOfActiveVertices = sim->m_vertices.NOFactive(); // 4 bytes
    unsigned int numberOfActiveHalfEdges = sim->m_halfEdges.NOFactive(); // 4 bytes

    if (Constants::getOBSTACLE())
    {
        nof_obstacle_particles = sim->m_obstacle->m_nof_obstacle_particles;
        obstacle_v_x = sim->m_obstacle->m_velocity.x;
        obstacle_v_y = sim->m_obstacle->m_velocity.y;
    }

    fwrite(fileEnding, sizeof(fileEnding), 1, pFile);
    fwrite(&outputVersion, sizeof(unsigned int), 1, pFile);
    fwrite(&time, sizeof(double), 1, pFile);
    fwrite(&output_counter, sizeof(unsigned int), 1, pFile);
    fwrite(&particlesCounter, sizeof(unsigned int), 1, pFile);
    fwrite(&ghostsCounter, sizeof(unsigned int), 1, pFile);
    _7_(fwrite(&nof_obstacle_particles, sizeof(unsigned int), 1, pFile);)
    _8_(fwrite(&obstacle_v_x, sizeof(double), 1, pFile);)
    _8_(fwrite(&obstacle_v_y, sizeof(double), 1, pFile);)

    fwrite(&numberOfActiveVertices, sizeof(unsigned int), 1, pFile);
    fwrite(&numberOfActiveHalfEdges, sizeof(unsigned int), 1, pFile);

    IO_PRINTF("output version: %u\n", outputVersion);
    IO_PRINTF("m_particles counter: %u\n", particlesCounter);

    // set position indicator
    fseek(pFile, headersize, SEEK_SET);

    unsigned int particleNofVertices = 0;
    HalfEdge* first = NULL;
    HalfEdge* current = NULL;
    fpos_t posAfterParticle;
    fpos_t posNofVertices;

    VoronoiParticle* particle = NULL;

    // output of the particles
    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];

        // reset
        particleNofVertices = 0;
        first = particle->m_outerComponent;
        current = particle->m_outerComponent;

        fwrite(&(particle->m_location), pointSize, 1, pFile);

        IO_PRINTF("particle %u: %f %f %f %f %f %f %f\n", i + 1, particle->m_location.x, particle->m_location.y,
                  particle->m_velocity.x, particle->m_velocity.y,
                  particle->m_density,
                  particle->m_pressure,
                  particle->m_color);

        fwrite(&(particle->m_velocity), pointSize, 1, pFile);
        fwrite(&(particle->m_density), doubleSize, 1, pFile);
        fwrite(&(particle->m_pressure), doubleSize, 1, pFile);
        fwrite(&(particle->m_color), doubleSize, 1, pFile);

        // number of vertices of the particle (not yet known)
        fgetpos(pFile, &posNofVertices);
        fwrite(&particleNofVertices, sizeof(unsigned int), 1, pFile); //write value zero

        // walk around the particle
        //

        if (first != NULL)
        {
            assert(graph::isClosed(first)); // the cell has to be closed

            // walk
            do
            {
                fwrite(&(current->m_origin->m_position), pointSize, 1, pFile);
                particleNofVertices++;

                current = current->m_next;

            } while (first != current);

            // overwrite number of vertices (its currently zero)
            fgetpos(pFile, &posAfterParticle);

            fsetpos(pFile, &posNofVertices);
            fwrite(&particleNofVertices, sizeof(unsigned int), 1, pFile);

            fsetpos(pFile, &posAfterParticle);
        }
        else
        {
            fprintf(Constants::getP_WARNINGS_FILE(), "WARNING in Sim::binaryOutput: particle has no outer component\n");

        }
    }

    GhostParticle* ghost = NULL;

    // output of the ghost particles
    for (unsigned int j = 0; j < sim->m_ghosts.counter(); j++)
    {
        ghost = &sim->m_ghosts[j];
        fwrite(&(ghost->m_location), pointSize, 1, pFile);

        _6_(fwrite(&ghost->m_flag_1, int_size, 1, pFile));

        IO_PRINTF("ghost particle %u: %f %f\n", j, ghost->m_location.x, ghost->m_location.y);
    }

    // close file
    fclose(pFile);
}

// read in data from binary file
Sim* bin_io::binaryInput(const char* filename)
{
    size_t pointSize = sizeof(Point);
    size_t doubleSize = sizeof(double);
    size_t int_size = sizeof(int);

    IO_PRINTF("\n\n---binary input of file %s---\n", filename);

    FILE* pFile;

    pFile = fopen(filename, "r");

    if (pFile == NULL)
    {
        fprintf(stderr, "ERROR in Sim::binaryInput: cant't open file %s\n", filename);
        exit(1);
    }

    char fileEnding[] = "xx"; //3 bytes
    unsigned int outputVersion = 0; //4 bytes
    unsigned int particlesCounter = 0; //4 bytes
    unsigned int ghostsCounter = 0; //4 bytes

    unsigned int nof_obstacle_particles = 0; // 4bytes
    double obstacle_v_x = 0;
    double obstacle_v_y = 0;

    fread(fileEnding, sizeof(fileEnding), 1, pFile);
    fread(&outputVersion, sizeof(unsigned int), 1, pFile);
    _3_(fread(&sim_core::simulation_time, sizeof(double), 1, pFile));
    _5_(fread(&sim_core::counter, sizeof(unsigned int), 1, pFile));
    _5_(sim_core::counter++);
    fread(&particlesCounter, sizeof(unsigned int), 1, pFile);
    _4_(fread(&ghostsCounter, sizeof(unsigned int), 1, pFile));

    _7_(fread(&nof_obstacle_particles, sizeof(unsigned int), 1, pFile);)
    _8_(fread(&obstacle_v_x, sizeof(unsigned int), 1, pFile);)
    _8_(fread(&obstacle_v_y, sizeof(unsigned int), 1, pFile);)

    if (!(fileEnding[0] == 'h' && fileEnding[1] == 'c'))
    {
        fprintf(stderr, "ERROR in Sim::binaryInput: input file is no .hc file\n");
        exit(1);
    }

    IO_PRINTF("output version: %u\n", outputVersion);
    IO_PRINTF("simulation time: %f\n", sim_core::simulation_time);
    IO_PRINTF("counter: %f\n", sim_core::counter);
    IO_PRINTF("m_particles counter: %u\n", particlesCounter);
    IO_PRINTF("m_ghosts counter: %u\n", ghostsCounter);

    // create a new simulation
    Sim* sim = new Sim(particlesCounter + 10000);

    // set position indicator
    fseek(pFile, headersize, SEEK_SET);

    unsigned int particleNofVertices = 0;
    Point vertex = Point(0, 0);

    // read in particle locations and primitive variables
    for (unsigned int i = 0; i < particlesCounter; i++)
    {
        particleNofVertices = 0;

        fread(&(sim->m_particles.current().m_location), pointSize, 1, pFile);
        fread(&(sim->m_particles.current().m_velocity), pointSize, 1, pFile);
        fread(&(sim->m_particles.current().m_density), doubleSize, 1, pFile);
        fread(&(sim->m_particles.current().m_pressure), doubleSize, 1, pFile);
        _2_(fread(&(sim->m_particles.current().m_color), doubleSize, 1, pFile));

        IO_PRINTF("particle %u: %f %f %f %f %f %f %f\n", i + 1, sim->m_particles.current().m_location.x,
                  sim->m_particles.current().m_location.y,
                  sim->m_particles.current().m_velocity.x, sim->m_particles.current().m_velocity.y,
                  sim->m_particles.current().m_density,
                  sim->m_particles.current().m_pressure,
                  sim->m_particles.current().m_color);

        sim->m_particles.counterInc();

        fread(&particleNofVertices, sizeof(unsigned int), 1, pFile);

        IO_PRINTF("vertices: %u\n", particleNofVertices);

        // skip the vertices
        for (unsigned int j = 0; j < particleNofVertices; j++)
        {
            fread(&vertex, pointSize, 1, pFile);
        }
    }

    // read in ghost particles
    _4_(for (unsigned int j = 0; j < ghostsCounter; j++)
        {
            fread(&(sim->m_ghosts.current().m_location), pointSize, 1, pFile);
            IO_PRINTF("ghost particle %u: %f %f\n", j, sim->m_ghosts.current().m_location.x,
                      sim->m_ghosts.current().m_location.y);

            _6_(fread(&(sim->m_ghosts.current().m_flag_1), int_size, 1, pFile);)

            sim->m_ghosts.counterInc();
        })

    fclose(pFile);

    // obstacle
    if (nof_obstacle_particles != 0)
    {
        sim->m_obstacle = new Obstacle(sim);
        sim->m_obstacle->bin_input(nof_obstacle_particles);
        sim->m_obstacle->m_velocity.x = obstacle_v_x;
        sim->m_obstacle->m_velocity.y = obstacle_v_y;
    }

    return sim;
}

// helper functions
static Point intersection_x(Point& a, Point& b, double x_0)
{
    IO_PRINTF("a: %f %f, b: %f %f\n", a.x, a.y, b.x, b.y);

    double t = 0;

    Point relative = a - b;

    if (relative.x == 0)
    {
        fprintf(stderr, "ERROR in intersection_x: denominator is zero\n");
    }

    t = (x_0 - b.x) / relative.x;

    IO_PRINTF("intersection: %f %f\n\n", (b + relative * t).x, (b + relative * t).y);

    Point intersection = Point(0, 0);
    intersection.y = b.y + relative.y * t;
    intersection.x = x_0;

    return intersection;
}

static Point intersection_y(Point& a, Point& b, double y_0)
{
    double t = 0;

    Point relative = a - b;

    if (relative.y == 0)
    {
        fprintf(stderr, "ERROR in intersection_y: denominator is zero\n");
    }

    t = (y_0 - b.y) / relative.y;

    Point intersection = Point(0, 0);
    intersection.x = b.x + relative.x * t;
    intersection.y = y_0;

    return intersection;
}

static bool is_inside(Point& p, double x_min, double y_min, double x_max, double y_max)
{
    return (p.x >= x_min && p.x <= x_max && p.y >= y_min && p.y <= y_max);
}

static bool is_border_cell(std::vector <Point>& vectorList, double x_min, double y_min, double x_max, double y_max)
{
    for (unsigned i = 0; i < vectorList.size(); i++)
    {
        if (!is_inside(vectorList[i], x_min, y_min, x_max, y_max))
        {
            return true;
        }
    }

    return false;
}

static void shift(std::vector <Point>& vectorList, double dx, double dy)
{
    for (unsigned i = 0; i < vectorList.size(); i++)
    {
        vectorList[i].x += dx;
        vectorList[i].y += dy;
    }
}

static void insert_intersections_x(std::vector <Point>& vectorList, double x_0)
{
    vectorList.push_back(vectorList[0]);

    std::vector <Point> vectorListTmp = vectorList;
    vectorList.clear();

    Point previous = Point(0, 0);
    bool previous_inside;
    Point intersection = Point(0, 0);

    vectorList.push_back(vectorListTmp[0]);
    previous = vectorListTmp[0];

    if (vectorListTmp[0].x > x_0) // first point is outside
    {
        previous_inside = false;
    }
    else // first point is inside
    {
        previous_inside = true;
    }

    for (unsigned i = 1; i < vectorListTmp.size(); i++)
    {
        if (vectorListTmp[i].x > x_0) // point is outside
        {
            if (previous_inside) // intersection
            {
                intersection = intersection_x(previous, vectorListTmp[i], x_0);
                vectorList.push_back(intersection);
            }

            vectorList.push_back(vectorListTmp[i]);
            previous_inside = false;
        }
        else // point is inside
        {
            if (!previous_inside) // intersection
            {
                intersection = intersection_x(vectorListTmp[i], previous, x_0);
                vectorList.push_back(intersection);
            }

            vectorList.push_back(vectorListTmp[i]);
            previous_inside = true;
        }

        previous = vectorListTmp[i];
    }
}

static void insert_intersections_y(std::vector <Point>& vectorList, double y_0)
{
    vectorList.push_back(vectorList[0]);

    std::vector <Point> vectorListTmp = vectorList;
    vectorList.clear();

    Point previous = Point(0, 0);
    bool previous_inside;
    Point intersection = Point(0, 0);

    vectorList.push_back(vectorListTmp[0]);
    previous = vectorListTmp[0];

    if (vectorListTmp[0].y > y_0) // first point is outside
    {
        previous_inside = false;
    }
    else // first point is inside
    {
        previous_inside = true;
    }

    for (unsigned i = 1; i < vectorListTmp.size(); i++)
    {
        if (vectorListTmp[i].y > y_0) // point is outside
        {
            if (previous_inside) // intersection
            {
                intersection = intersection_y(previous, vectorListTmp[i], y_0);
                vectorList.push_back(intersection);
            }

            vectorList.push_back(vectorListTmp[i]);
            previous_inside = false;
        }
        else // point is inside
        {
            if (!previous_inside) // intersection
            {
                intersection = intersection_y(vectorListTmp[i], previous, y_0);
                vectorList.push_back(intersection);
            }

            vectorList.push_back(vectorListTmp[i]);
            previous_inside = true;
        }

        previous = vectorListTmp[i];
    }
}

static void cut_y_top(std::vector <Point>& vectorList, double y_0)
{
    insert_intersections_y(vectorList, y_0);

    std::vector <Point> vectorListTmp;

    // remove vertices outside
    for (unsigned int k = 0; k < vectorList.size(); k++)
    {
        if (vectorList[k].y <= y_0)
        {
            vectorListTmp.push_back(vectorList[k]);
        }
    }

    vectorList = vectorListTmp;
}

static void cut_y_bottom(std::vector <Point>& vectorList, double y_0)
{
    insert_intersections_y(vectorList, y_0);

    std::vector <Point> vectorListTmp;

    // remove vertices outside
    for (unsigned int k = 0; k < vectorList.size(); k++)
    {
        if (vectorList[k].y >= y_0)
        {
            vectorListTmp.push_back(vectorList[k]);
        }
    }

    vectorList = vectorListTmp;
}

static void cut_x_right(std::vector <Point>& vectorList, double x_0)
{
    insert_intersections_x(vectorList, x_0);

    std::vector <Point> vectorListTmp;

    // remove vertices outside
    for (unsigned int k = 0; k < vectorList.size(); k++)
    {
        if (vectorList[k].x <= x_0)
        {
            vectorListTmp.push_back(vectorList[k]);
        }
    }

    vectorList = vectorListTmp;
}

static void cut_x_left(std::vector <Point>& vectorList, double x_0)
{
    insert_intersections_x(vectorList, x_0);

    std::vector <Point> vectorListTmp;

    // remove vertices outside
    for (unsigned int k = 0; k < vectorList.size(); k++)
    {
        if (vectorList[k].x >= x_0)
        {
            vectorListTmp.push_back(vectorList[k]);
        }
    }

    vectorList = vectorListTmp;
}

// convert the ghost data to ASCII data
void bin_io::ghost_convert(const char* filename, Sim* sim)
{
    unsigned int name_length = strlen(filename) - 3;

    char name[name_length + 1];
    name[0] = '\0';
    strncat(name, filename, name_length);

    char ghostsFilename[64];
    strcpy(ghostsFilename, name);
    const char ending[] = ".ghosts";
    strcat(ghostsFilename, ending);

    char facesFilename[64];
    strcpy(facesFilename, name);
    const char endingFaces[] = ".ghost_faces";
    strcat(facesFilename, endingFaces);

    FILE* pGhosts = fopen(ghostsFilename, "w");
    FILE* pFaces = fopen(facesFilename, "w");

    GhostParticle* ghost = NULL;

    HalfEdge* first = NULL;
    HalfEdge* current = NULL;

    for (unsigned int i = 0; i < sim->m_ghosts.counter(); i++)
    {
        ghost = &sim->m_ghosts[i];

        // location
        fprintf(pGhosts, "%f	%f\n", ghost->m_location.x, ghost->m_location.y);

        // face
        first = ghost->m_outerComponent;
        current = ghost->m_outerComponent;

        if (graph::isClosed(first))//only closed faces
        {
            do
            {
                fprintf(pFaces, "%f	%f\n", current->m_origin->m_position.x, current->m_origin->m_position.y);

                current = current->m_next;

            } while (first != current);

            // print the last point twice
            fprintf(pFaces, "%f	%f\n\n\n", current->m_origin->m_position.x, current->m_origin->m_position.y);
        }
    }

    fclose(pGhosts);
    fclose(pFaces);
}

// convert a binary file to ASCII
void bin_io::convert(const char* filename)
{
    FILE* pFile;

    unsigned int name_length = strlen(filename) - 3;

    // check whether the filename length is > 3
    if (name_length <= 0)
    {
        fprintf(stderr, "ERROR in bin_io::convert: invalid filename %s\n", filename);
        exit(1);
    }

    // check whether the file ending is .hc
    if (!(filename[name_length] == '.' && filename[name_length + 1] == 'h' && filename[name_length + 2] == 'c'))
    {
        fprintf(stderr, "ERROR in bin_io::convert: filename %s has not the ending .hc\n", filename);
        exit(1);
    }

    pFile = fopen(filename, "r"); // open the binary file

    if (pFile == NULL)
    {
        fprintf(stderr, "ERROR in bin_io::convert: can't open file %s\n", filename);
        exit(1);
    }

    // initialization of the output files
    //

    char name[name_length + 1];
    name[0] = '\0';
    strncat(name, filename, name_length);

    char dataFilename[64];
    strcpy(dataFilename, name);
    const char ending[] = ".data";
    strcat(dataFilename, ending);

    char densityFilename[64];
    strcpy(densityFilename, name);
    const char endingDensity[] = ".density2D";
    strcat(densityFilename, endingDensity);

    char pressureFilename[64];
    strcpy(pressureFilename, name);
    const char endingPressure[] = ".pressure2D";
    strcat(pressureFilename, endingPressure);

    char colorFilename[64];
    strcpy(colorFilename, name);
    const char endingColor[] = ".color2D";
    strcat(colorFilename, endingColor);

    char facesFilename[64];
    strcpy(facesFilename, name);
    const char endingFaces[] = ".faces";
    strcat(facesFilename, endingFaces);

    char obstacleFilename[64];
    strcpy(obstacleFilename, name);
    const char endingObstacle[] = ".obstacle";
    strcat(obstacleFilename, endingObstacle);

    // read in data
    //

    size_t pointSize = sizeof(Point);
    size_t doubleSize = sizeof(double);

    char fileEnding[] = "xx"; //3 bytes
    unsigned int outputVersion = 0; //4 bytes
    double time = 0; //8 bytes
    unsigned int output_counter = 0; //4 bytes
    unsigned int particlesCounter = 0; //4 bytes
    unsigned int ghostsCounter = 0; //4 bytes
    unsigned int nof_obstacle_particles = 0; //4 bytes

    fread(fileEnding, sizeof(fileEnding), 1, pFile);
    fread(&outputVersion, sizeof(unsigned int), 1, pFile);
    _3_(fread(&time, sizeof(double), 1, pFile));
    _5_(fread(&output_counter, sizeof(unsigned int), 1, pFile));
    fread(&particlesCounter, sizeof(unsigned int), 1, pFile);
    _4_(fread(&ghostsCounter, sizeof(unsigned int), 1, pFile));

    if (options::obstacle)
    {
        fread(&nof_obstacle_particles, sizeof(unsigned int), 1, pFile);
    }

    if (!(fileEnding[0] == 'h' && fileEnding[1] == 'c'))
    {
        fprintf(stderr, "ERROR in bin_io::convert: input file is no .hc file\n");
        exit(1);
    }

    FILE* pData;
    FILE* pDensity;
    FILE* pPressure;
    FILE* pColor;
    FILE* pFaces;
    FILE* pObstacle;

    if (options::obstacle)
    {
        pObstacle = fopen(obstacleFilename, "w");
    }
    else
    {
        pObstacle = fopen("/dev/null", "w");
    }

    if (options::density_only)
    {
        pData = fopen("/dev/null", "w");
        pDensity = fopen(densityFilename, "w");
        pPressure = fopen("/dev/null", "w");
        pColor = fopen("/dev/null", "w");
        pFaces = fopen("/dev/null", "w");
    }
    else if (options::color_only)
    {
        pData = fopen("/dev/null", "w");
        pDensity = fopen("/dev/null", "w");
        pPressure = fopen("/dev/null", "w");
        pColor = fopen(colorFilename, "w");
        pFaces = fopen("/dev/null", "w");
    }
    else
    {
        pData = fopen(dataFilename, "w");
        pDensity = fopen(densityFilename, "w");
        pPressure = fopen(pressureFilename, "w");
        pColor = fopen(colorFilename, "w");
        pFaces = fopen(facesFilename, "w");
    }

    printf("output version: %u\n", outputVersion);
    printf("m_particles counter: %u\n", particlesCounter);

    // set position indicator
    fseek(pFile, headersize, SEEK_SET);

    Point location = Point(0, 0);
    Vector velocity = Vector(0, 0);
    double density = 0;
    double pressure = 0;
    double color = 0;
    Point vertex = Point(0, 0);
    std::vector <Point> vectorList;
    vectorList.reserve(16);

    unsigned int particleNofVertices = 0;

    fprintf(pData, "#time: %f\n", time);
    fprintf(pData, "#x		y		v_x		v_y		rho		p		color\n");

    if (options::cut)
    {
        // box parameters
        double x_min = options::x_min;
        double y_min = options::y_min;
        double x_max = options::x_max;
        double y_max = options::y_max;

        double delta_x_cutbox = x_max - x_min;
        double delta_y_cutbox = y_max - y_min;

        double x_middle = 0.5 * (x_min + x_max);
        double y_middle = 0.5 * (y_min + y_max);


        for (unsigned int n = 0; n < particlesCounter; n++)
        {
            // reset
            //

            vectorList.clear();

            // read in
            //

            fread(&location, pointSize, 1, pFile);
            fread(&velocity, pointSize, 1, pFile);
            fread(&density, doubleSize, 1, pFile);
            fread(&pressure, doubleSize, 1, pFile);
            _2_(fread(&color, doubleSize, 1, pFile));


            fread(&particleNofVertices, sizeof(unsigned int), 1, pFile);

            fprintf(pData, "%f	%f	%f	%f	%f	%f	%f\n", location.x, location.y, velocity.x, velocity.y,
                    density, pressure, color);

            for (unsigned int j = 0; j < particleNofVertices; j++)
            {
                fread(&vertex, pointSize, 1, pFile);
                vectorList.push_back(vertex);
            }

            if (is_border_cell(vectorList, x_min, y_min, x_max, y_max))
            {

                // color = 1;

                if (location.x > x_middle)
                {
                    cut_x_right(vectorList, x_max);
                }

                if (location.x < x_middle)
                {
                    cut_x_left(vectorList, x_min);
                }

                if (location.y > y_middle)
                {
                    cut_y_top(vectorList, y_max);
                }

                if (location.y < y_middle)
                {
                    cut_y_bottom(vectorList, y_min);
                }
            }

            // create 2D plotting data
            //

            int size = vectorList.size();

            int i(0);
            int j(size - 1);

            while (j - i > 1)
            {
                // density
                //

                fprintf(pDensity, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, density);
                fprintf(pDensity, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, density);

                fprintf(pDensity, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, density);
                fprintf(pDensity, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, density);

                // pressure
                //

                fprintf(pPressure, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, pressure);
                fprintf(pPressure, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, pressure);

                fprintf(pPressure, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, pressure);
                fprintf(pPressure, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, pressure);

                // color
                //

                fprintf(pColor, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, color);
                fprintf(pColor, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, color);

                fprintf(pColor, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, color);
                fprintf(pColor, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, color);

                if (options::obstacle)
                {
                    if (n < (nof_obstacle_particles - 1) / 2 || n == nof_obstacle_particles - 1)
                    {
                        fprintf(pObstacle, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, 0.);
                        fprintf(pObstacle, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, 0.);

                        fprintf(pObstacle, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, 0.);
                        fprintf(pObstacle, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, 0.);
                    }
                }

                i++;
                j--;
            }

            if (location.x > x_min - delta_x_cutbox && location.x < x_max + delta_x_cutbox &&
                location.y < y_max + delta_y_cutbox && location.y > y_min - delta_y_cutbox)
            {
                // create faces
                for (unsigned k = 0; k < vectorList.size(); k++)
                {
                    fprintf(pFaces, "%f	%f\n", vectorList[k].x, vectorList[k].y);
                }

                fprintf(pFaces, "%f	%f\n\n\n", vectorList[0].x, vectorList[0].y); //close face
            }
        }
    }

    else if (options::hybrid)
    {
        // read constants from file
        Constants::readConstants();

        std::vector <Point> vectorList2;
        vectorList2.reserve(16);
        bool border_cell = false;

        // box parameters
        double x_middle = (Constants::getBOX_XMAX() + Constants::getBOX_XMIN()) * 0.5;
        double delta_x = Constants::getBOX_XMAX() - Constants::getBOX_XMIN();

        for (unsigned int n = 0; n < particlesCounter; n++)
        {
            // reset
            //

            vectorList.clear();
            vectorList2.clear();

            // read in
            //

            fread(&location, pointSize, 1, pFile);
            fread(&velocity, pointSize, 1, pFile);
            fread(&density, doubleSize, 1, pFile);
            fread(&pressure, doubleSize, 1, pFile);
            _2_(fread(&color, doubleSize, 1, pFile));

            fread(&particleNofVertices, sizeof(unsigned int), 1, pFile);

            fprintf(pData, "%f	%f	%f	%f	%f	%f	%f\n", location.x, location.y, velocity.x, velocity.y,
                    density, pressure, color);

            for (unsigned int j = 0; j < particleNofVertices; j++)
            {
                fread(&vertex, pointSize, 1, pFile);
                vectorList.push_back(vertex);
            }

            if (is_border_cell(vectorList, Constants::getBOX_XMIN(), -std::numeric_limits<double>::max(),
                               Constants::getBOX_XMAX(), std::numeric_limits<double>::max()))
            {
                border_cell = true;
                // color = 1;
                vectorList2 = vectorList;


                if (location.x > x_middle)
                {
                    shift(vectorList2, -delta_x, 0);
                    cut_x_left(vectorList2, Constants::getBOX_XMIN());
                    cut_x_right(vectorList, Constants::getBOX_XMAX());
                }
                else
                {
                    shift(vectorList2, delta_x, 0);
                    cut_x_right(vectorList2, Constants::getBOX_XMAX());
                    cut_x_left(vectorList, Constants::getBOX_XMIN());
                }
            }
            else
            {
                border_cell = false;
            }

            // create 2D plotting data
            //

            int size = vectorList.size();

            if (size < 3)
            {
                fprintf(stderr, "ERROR in bin_io::convert: Voronoi particle ( %f , %f ) has less then three vertices\n",
                        location.x, location.y);
                exit(1);
            }

            int i(0);
            int j(size - 1);

            while (j - i > 1)
            {
                // density
                //

                fprintf(pDensity, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, density);
                fprintf(pDensity, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, density);

                fprintf(pDensity, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, density);
                fprintf(pDensity, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, density);

                // pressure
                //

                fprintf(pPressure, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, pressure);
                fprintf(pPressure, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, pressure);

                fprintf(pPressure, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, pressure);
                fprintf(pPressure, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, pressure);

                //color
                //

                fprintf(pColor, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, color);
                fprintf(pColor, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, color);

                fprintf(pColor, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, color);
                fprintf(pColor, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, color);

                if (options::obstacle)
                {
                    if (n < (nof_obstacle_particles - 1) / 2 || n == nof_obstacle_particles - 1)
                    {
                        fprintf(pObstacle, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, 0.);
                        fprintf(pObstacle, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, 0.);

                        fprintf(pObstacle, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, 0.);
                        fprintf(pObstacle, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, 0.);
                    }
                }

                i++;
                j--;
            }

            // create faces
            for (unsigned k = 0; k < vectorList.size(); k++)
            {
                fprintf(pFaces, "%f	%f\n", vectorList[k].x, vectorList[k].y);
            }

            fprintf(pFaces, "%f	%f\n\n\n", vectorList[0].x, vectorList[0].y); //close face

            // border cell
            if (border_cell)
            {
                // create 2D plotting data
                size = vectorList2.size();

                if (size < 3)
                {
                    fprintf(stderr,
                            "ERROR in bin_io::convert: Voronoi particle ( %f , %f ) has less then three vertices\n",
                            location.x, location.y);
                    exit(1);
                }

                int k(0);
                int l(size - 1);

                while (l - k > 1)
                {
                    // density
                    //

                    fprintf(pDensity, "%f	%f	%f\n", vectorList2[k].x, vectorList2[k].y, density);
                    fprintf(pDensity, "%f	%f	%f\n\n", vectorList2[k + 1].x, vectorList2[k + 1].y, density);

                    fprintf(pDensity, "%f	%f	%f\n", vectorList2[l].x, vectorList2[l].y, density);
                    fprintf(pDensity, "%f	%f	%f\n\n\n", vectorList2[l - 1].x, vectorList2[l - 1].y, density);

                    // pressure
                    //

                    fprintf(pPressure, "%f	%f	%f\n", vectorList2[k].x, vectorList2[k].y, pressure);
                    fprintf(pPressure, "%f	%f	%f\n\n", vectorList2[k + 1].x, vectorList2[k + 1].y, pressure);

                    fprintf(pPressure, "%f	%f	%f\n", vectorList2[l].x, vectorList2[l].y, pressure);
                    fprintf(pPressure, "%f	%f	%f\n\n\n", vectorList2[l - 1].x, vectorList2[l - 1].y, pressure);

                    // color
                    //

                    fprintf(pColor, "%f	%f	%f\n", vectorList2[k].x, vectorList2[k].y, color);
                    fprintf(pColor, "%f	%f	%f\n\n", vectorList2[k + 1].x, vectorList2[k + 1].y, color);

                    fprintf(pColor, "%f	%f	%f\n", vectorList2[l].x, vectorList2[l].y, color);
                    fprintf(pColor, "%f	%f	%f\n\n\n", vectorList2[l - 1].x, vectorList2[l - 1].y, color);

                    if (options::obstacle)
                    {
                        if (n < (nof_obstacle_particles - 1) / 2 || n == nof_obstacle_particles - 1)
                        {
                            fprintf(pObstacle, "%f	%f	%f\n", vectorList[k].x, vectorList[i].y, 0.);
                            fprintf(pObstacle, "%f	%f	%f\n\n", vectorList[k + 1].x, vectorList[k + 1].y, 0.);

                            fprintf(pObstacle, "%f	%f	%f\n", vectorList[l].x, vectorList[l].y, 0.);
                            fprintf(pObstacle, "%f	%f	%f\n\n\n", vectorList[l - 1].x, vectorList[l - 1].y, 0.);
                        }
                    }

                    k++;
                    l--;
                }

                // create faces
                for (unsigned k = 0; k < vectorList2.size(); k++)
                {
                    fprintf(pFaces, "%f	%f\n", vectorList2[k].x, vectorList2[k].y);
                }
                fprintf(pFaces, "%f	%f\n\n\n", vectorList2[0].x, vectorList2[0].y); //close face
            }
        }
    }

    else // convert default
    {
        for (unsigned int n = 0; n < particlesCounter; n++)
        {
            // reset
            //

            vectorList.clear();

            // read in
            //

            fread(&location, pointSize, 1, pFile);
            fread(&velocity, pointSize, 1, pFile);
            fread(&density, doubleSize, 1, pFile);
            fread(&pressure, doubleSize, 1, pFile);
            _2_(fread(&color, doubleSize, 1, pFile));

            fread(&particleNofVertices, sizeof(unsigned int), 1, pFile);

            fprintf(pData, "%f	%f	%f	%f	%f	%f	%f\n", location.x, location.y, velocity.x, velocity.y,
                    density, pressure, color);

            for (unsigned int j = 0; j < particleNofVertices; j++)
            {
                fread(&vertex, pointSize, 1, pFile);
                vectorList.push_back(vertex);
            }

            // create 2D plotting data
            //

            int size = vectorList.size();

            if (size < 3)
            {
                fprintf(stderr, "ERROR in bin_io::convert: Voronoi particle ( %f , %f ) has less then three vertices\n",
                        location.x, location.y);
                exit(1);
            }

            int i(0);
            int j(size - 1);

            while (j - i > 1)
            {
                // density
                //

                fprintf(pDensity, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, density);
                fprintf(pDensity, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, density);

                fprintf(pDensity, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, density);
                fprintf(pDensity, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, density);

                // pressure
                //

                fprintf(pPressure, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, pressure);
                fprintf(pPressure, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, pressure);

                fprintf(pPressure, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, pressure);
                fprintf(pPressure, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, pressure);

                // color
                //

                fprintf(pColor, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, color);
                fprintf(pColor, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, color);

                fprintf(pColor, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, color);
                fprintf(pColor, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, color);

                if (options::obstacle)
                {
                    if (n < (nof_obstacle_particles - 1) / 2 || n == nof_obstacle_particles - 1)
                    {
                        fprintf(pObstacle, "%f	%f	%f\n", vectorList[i].x, vectorList[i].y, 0.);
                        fprintf(pObstacle, "%f	%f	%f\n\n", vectorList[i + 1].x, vectorList[i + 1].y, 0.);

                        fprintf(pObstacle, "%f	%f	%f\n", vectorList[j].x, vectorList[j].y, 0.);
                        fprintf(pObstacle, "%f	%f	%f\n\n\n", vectorList[j - 1].x, vectorList[j - 1].y, 0.);
                    }
                }

                i++;
                j--;
            }

            // create faces
            for (unsigned k = 0; k < vectorList.size(); k++)
            {
                fprintf(pFaces, "%f	%f\n", vectorList[k].x, vectorList[k].y);
            }

            fprintf(pFaces, "%f	%f\n\n\n", vectorList[0].x, vectorList[0].y); //close face
        }
    }

    fclose(pData);
    fclose(pDensity);
    fclose(pPressure);
    fclose(pColor);
    fclose(pFaces);
    fclose(pFile);
    fclose(pObstacle);
}

// convert an integer to a filename, e.g. 25 => 00025.hc
void bin_io::int_to_filename(char* filename, int filenumber)
{
    sprintf(filename, "%05u.hc", filenumber);
}

// convert a filename to an integer, e.g. 00004.hc => 4
int bin_io::filename_to_int(const char* filename)
{
    unsigned int name_length = strlen(filename) - 3;
    char name[name_length + 1];
    name[0] = '\0';
    strncat(name, filename, name_length);
    return (atoi(name));
}

// does the file exist?
bool bin_io::file_exists(const char* filename)
{
    FILE* pFile;

    pFile = fopen(filename, "r");
    bool result = (pFile != NULL);

    if (result)
    {
        fclose(pFile);
    }

    return result;
}