#include "config.h"
#include "info.h"
#include "constants.h"

// print the header
void info::print_header()
{
    fprintf(Constants::getP_INFO_FILE(), "---info file---\n\n");
}

// print version information into the info file
void info::print_version()
{
    fprintf(Constants::getP_INFO_FILE(), "\n-version-\n\n");
    fprintf(Constants::getP_INFO_FILE(), "this file was compiled on %s, %s.\n", __DATE__, __TIME__);
    fprintf(Constants::getP_INFO_FILE(), "last git commit: %s\n", GIT_COMMIT);
}

// print the program arguments
void info::print_args(int argc, char** argv)
{
    fprintf(Constants::getP_INFO_FILE(), "\n-call-\n\n");

    for (int i = 0; i < argc; i++)
    {
        fprintf(Constants::getP_INFO_FILE(), "%s ", argv[i]);
    }

    fprintf(Constants::getP_INFO_FILE(), "\n");
}

// print general info
void info::print_general(Sim* sim)
{
    fprintf(Constants::getP_INFO_FILE(), "\n-general-\n\n");
    fprintf(Constants::getP_INFO_FILE(), "random seed: %u\n", Constants::getSEED());
    fprintf(Constants::getP_INFO_FILE(), "number of voronoi particles: %d\n", sim->m_particles.counter());
}

// print duration
void info::print_duration(unsigned long int duration, unsigned int steps)
{
    unsigned int hours = duration / 3600;
    unsigned int minutes = (duration % 3600) / 60;
    unsigned int seconds = (duration % 60);

    fprintf(Constants::getP_INFO_FILE(), "\n-done-\n\n");
    fprintf(Constants::getP_INFO_FILE(), "simulation steps: %d\n", steps);
    fprintf(Constants::getP_INFO_FILE(), "duration: %d hours, %d minutes, %d seconds\n", hours, minutes, seconds);

    fprintf(Constants::getP_INFO_FILE(), "\n\n\n\n");
}