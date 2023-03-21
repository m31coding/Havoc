#ifndef SECTION_H
#define SECTION_H

#include <cstdio>

class Sim;

namespace tools
{
    /// use a havoc tool
    bool use_a_tool();

    /// do a cross section and write ASCII data
    void cross_section(const char* filename);

    /// calculate gradients and write ASCII data
    void calculate_gradients(const char* filename);

    /// find minimum / maximum value of a primitive variable
    void find_minimum_maximum(Sim* sim);

    /// print minimum and maximum pressure
    void print_minimum_maximum(FILE* pFile);

    /// color a particle
    void color(double x, double y, double color, Sim* sim);

    /// color the particles according to their volume
    void color_volume(Sim* sim);

    /// color the particles according to their roundness
    void color_roundness(Sim* sim);

    /// color the particles according to their number of neighbours
    void color_nof_neighbours(Sim* sim);

    /// color the particles according to their vorticity
    void color_vorticity(Sim* sim);

    /// downsampling
    void downsampling(unsigned int pixel_per_length, char* filename, Sim* sim);
};

#endif