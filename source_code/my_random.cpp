#include "my_random.h"
#include "config.h"
#include "debug.h"
#include "constants.h"

// initialize rand()
void myRandom::random_ini()
{
    SD_PRINTF("initializing seed for srand()\n");

    if (Constants::getSEED() != 0)
    {
        srand(Constants::getSEED());
        SD_PRINTF("\tseed: %u\n", Constants::getSEED());
    }
    else
    {
        FILE* urandomFile;
        unsigned int seed;
        urandomFile = fopen("/dev/urandom", "r");
        if (urandomFile != NULL)
        {
            if (fread(&seed, sizeof(seed), 1, urandomFile) < 1)
            {
                fprintf(stderr, "ERROR: 0 bytes read from urandom\n");
                fclose(urandomFile);
                exit(1);
            }
            fclose(urandomFile);
        }
        else
        {
            fprintf(stderr, "ERROR: can't open urandom\n");
            exit(1);
        }

        srand(seed);
        SD_PRINTF("\trandom seed: %u\n", seed);
        Constants::setSEED(seed);
    }
}

// create a random double between two given doubles
double myRandom::randomDouble(double lowerBound, double upperBound)
{
    return ((upperBound - lowerBound) / RAND_MAX * rand() + lowerBound);
}

// create a random integer in a given interval
int myRandom::randomInteger(int lowerBound, int upperBound)
{
    return ((rand() % (upperBound - lowerBound + 1)) + lowerBound);
}