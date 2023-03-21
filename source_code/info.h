#ifndef INFO_H
#define INFO_H

#include "sim.h"

namespace info
{
    /// print the header
    void print_header();

    /// print version information into the info file
    void print_version();

    /// print the program arguments
    void print_args(int argc, char** argv);

    /// print general info
    void print_general(Sim* sim);

    /// print duration
    void print_duration(unsigned long int duration, unsigned int steps);
};

#endif