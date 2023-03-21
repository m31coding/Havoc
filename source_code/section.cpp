#include "config.h"
#include "section.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "options.h"
#include "r2.h"
#include "sim.h"
#include "binary_input_output.h"
#include "sim_graph.h"
#include "sweepline.h"
#include "point_location.h"

// do a cross section and write ASCII data
void section::cross_section(const char* filename)
{
    unsigned int name_length = strlen(filename) - 3;

    // check whether the filename length is > 3
    if (name_length <= 0)
    {
        fprintf(stderr, "ERROR in bin_io::cross_section: invalid filename %s\n", filename);
        exit(1);
    }

    // check whether the file ending is .hc
    if (!(filename[name_length] == '.' && filename[name_length + 1] == 'h' && filename[name_length + 2] == 'c'))
    {
        fprintf(stderr, "ERROR in bin_io::cross_section: filename %s has not the ending .hc\n", filename);
        exit(1);
    }

    // binary input and calculation of the tesselation
    Sim* sim = bin_io::binaryInput(filename);
    sim->create_box();
    sim->eventQueueIni();
    sim->eventQueueIniGhosts();
    sweepline::makeTesselation(sim);
    sim_graph::assertRealCellsClosed(sim);

    //initialization of the output file
    char name[name_length + 1];
    name[0] = '\0';
    strncat(name, filename, name_length);

    char sectionFilename[64];
    strcpy(sectionFilename, name);
    const char ending[] = ".section";
    strcat(sectionFilename, ending);

    FILE* pSection;

    pSection = fopen(sectionFilename, "w");

    fprintf(pSection, "#section: px = %f, py = %f, qx = %f, qy = %f, nofp = %d\n", options::p.x, options::p.y,
            options::q.x, options::q.y, options::nofp);

    assert(options::var != 0);

    if (options::var == 1)
    {
        fprintf(pSection, "#x		y		l		rho\n");
    }
    else if (options::var == 2)
    {
        fprintf(pSection, "#x		y		l		p\n");
    }

    Vector relative = options::q - options::p;
    relative /= (options::nofp - 1); // nofp >= 2, see options::check_options

    double dl = sqrt(relative.length_square()); // distance between two points
    double l = 0;

    Point start = options::p;
    Point current = Point(0, 0);

    VoronoiParticle* cell = NULL;

    for (unsigned int i = 0; i < options::nofp; i++)
    {
        printf("%d\n", i);

        current = start + relative * i;
        l = i * dl;

        printf("current: %f %f \n", current.x, current.y);
        cell = point_location::findCell(&current, sim, cell);

        if (options::var == 1)
        {
            fprintf(pSection, "%f	%f	%f	%f\n", current.x, current.y, l, cell->m_density);
        }
        else if (options::var == 2)
        {
            fprintf(pSection, "%f	%f	%f	%f\n", current.x, current.y, l, cell->m_pressure);
        }
    }

    sim->freeGraph();
    delete sim;
    fclose(pSection);
}