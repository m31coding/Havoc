#ifndef BINARY_INPUT_OUTPUT_H
#define BINARY_INPUT_OUTPUT_H

#include "sim.h"

namespace bin_io
{
    /// write most important data into a binary file
    void binaryOutput(const char* file, unsigned int filenumber, Sim* sim);

    /// read in data from binary file (particle locations and primitive variables), returns a pointer to a new simulation
    Sim* binaryInput(const char* filename);

    /// convert a binary file to ASCII
    void convert(const char* filename);

    /// convert the ghost data to ASCII data
    void ghost_convert(const char* filename, Sim* sim);

    /// convert a integer to a filename, e.g. 25 => 00025.hc
    void int_to_filename(char* filename, int filenumber);

    /// convert a filename to a integer, e.g. 00004.hc => 4
    int filename_to_int(const char* filename);

    /// does the file exist?
    bool file_exists(const char* filename);
};

#endif