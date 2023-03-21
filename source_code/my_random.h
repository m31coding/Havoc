#ifndef MY_RANDOM_H
#define MY_RANDOM_H

#include <cstdlib>

namespace myRandom
{
    /// initialize rand()
    void random_ini();

    /**
    create a random double in a given interval
    */
    double randomDouble(double lowerBound, double upperBound);

    /**
    create a random integer in a given interval
    */
    int randomInteger(int lowerBound, int upperBound);
};

#endif