#ifndef OPTIONS_H
#define OPTIONS_H

class Point;

namespace options
{
    extern char* input_filename;

    extern char* convert_filename;
    extern bool zero_grad;
    extern double cut_percent;
    extern bool hybrid;
    extern bool cut;
    extern bool color_only;
    extern bool density_only;
    extern double x_min;
    extern double y_min;
    extern double x_max;
    extern double y_max;

    extern char* graph_filename;

    extern char* section_filename;
    extern Point p;
    extern Point q;
    extern unsigned int nofp;

    extern char* gradient_filename;

    extern unsigned int var; // 1 for density, 2 for pressure, 3 for velocity_x, 4 for velocity_y

    extern char* minmax_filename;
    extern bool subsequent;

    extern double color;
    extern bool color_set;
    extern char* color_filename;

    extern char* color_volume_filename;
    extern char* color_roundness_filename;
    extern char* color_nof_neighbours_filename;
    extern char* color_vorticity_filename;

    extern char* indice_filename;
    extern unsigned int object;
    extern bool object_set;

    extern char* convert_ghost_filename;

    extern bool obstacle;

    extern bool whatever;
    extern char* whatever_filename;

    extern char* add_up_filename;

    extern bool abs_vorticity; // absolute value of the vorticity
    extern bool kinetic_energy;
    extern bool turbulence_kinetic_energy;
    extern bool energy;
    extern bool kinetic_energy_y;
    extern bool mass_std_dev_square;
    extern bool color_on_the_fly;

    extern char* downsampling_filename;
    extern unsigned int pixel_per_length;

    /// parse program options
    void parse_options(int argc, char** argv);

    /// check the options
    void check_options(int argc, char** argv);

    /// print usage of the program
    void usage(int argc, char** argv);

    /// free memory
    void clear();
};

#endif