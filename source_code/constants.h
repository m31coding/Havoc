#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cstdio>
#include "r2.h"

class Sim;

/**
a singleton class storing all constants
**/
class Constants
{
public:

    static Constants* m_constants; ///< static member, stores a single instance

    /// get the instance
    inline static Constants* get()
    {
        return m_constants;
    }

    /// delete the singleton
    inline static void del()
    {
        fclose(m_constants->P_INFO_FILE);
        fclose(m_constants->P_WARNINGS_FILE);
        delete[] m_constants->INFO_FILENAME;
        delete[] m_constants->WARNINGS_FILENAME;
        delete m_constants;
    }

    /// read the constants from file "constants"
    static void readConstants();

    /// print the constants in a file
    static void fprint(FILE* pFile);

    /// create further constants
    static void createFurtherConstants();

    /// check and adapt constants
    static void check_and_adapt();

    // getters
    inline static bool getADAPTIVE_MESH_REFINEMENT()
    {
        return m_constants->ADAPTIVE_MESH_REFINEMENT;
    }

    inline static double getMAX_CELL_VOLUME()
    {
        return m_constants->MAX_CELL_VOLUME;
    }

    inline static char* getINFO_FILENAME()
    {
        return m_constants->INFO_FILENAME;
    }

    inline static char* getWARNINGS_FILENAME()
    {
        return m_constants->WARNINGS_FILENAME;
    }

    inline static unsigned int getSEED()
    {
        return m_constants->SEED;
    }

    inline static unsigned int getBOX_TYPE()
    {
        return m_constants->BOX_TYPE;
    }

    inline static double getBOX_XMIN()
    {
        return m_constants->BOX_XMIN;
    }

    inline static double getBOX_XMAX()
    {
        return m_constants->BOX_XMAX;
    }

    inline static double getBOX_YMIN()
    {
        return m_constants->BOX_YMIN;
    }

    inline static double getBOX_YMAX()
    {
        return m_constants->BOX_YMAX;
    }

    inline static unsigned int getNOF_PARTICLES_X()
    {
        return m_constants->NOF_PARTICLES_X;
    }

    inline static unsigned int getNOF_PARTICLES_Y()
    {
        return m_constants->NOF_PARTICLES_Y;
    }

    inline static bool getRELAX_FIRST_GRAPH()
    {
        return m_constants->RELAX_FIRST_GRAPH;
    }

    inline static bool getDEL_ZERO_EDGES()
    {
        return m_constants->DEL_ZERO_EDGES;
    }

    inline static double getDELTA_T_OUTPUT()
    {
        return m_constants->DELTA_T_OUTPUT;
    }

    inline static int getEVERY_NTH_STEPS_OUTPUT()
    {
        return m_constants->EVERY_NTH_STEPS_OUTPUT;
    }

    inline static double getDELTA_T()
    {
        return m_constants->DELTA_T;
    }

    inline static double getMAX_T()
    {
        return m_constants->MAX_T;
    }

    inline static int getMAX_STEPS()
    {
        return m_constants->MAX_STEPS;
    }

    inline static FILE* getP_INFO_FILE()
    {
        return m_constants->P_INFO_FILE;
    }

    inline static FILE* getP_WARNINGS_FILE()
    {
        return m_constants->P_WARNINGS_FILE;
    }

    inline static double getGAMMA()
    {
        return m_constants->GAMMA;
    }

    inline static double getCFL()
    {
        return m_constants->CFL;
    }

    inline static unsigned int getVELOCITY_UPDATE_FUNCTION()
    {
        return m_constants->VELOCITY_UPDATE_FUNCTION;
    }

    inline static void (* function_updateParticleVelocities())(Sim*)
    {
        return m_constants->updateParticleVelocities;
    }

    inline static unsigned int getSLOPE_LIMITER()
    {
        return m_constants->SLOPE_LIMITER;
    }

    inline static void (* function_slopeLimitGradients())(Sim*)
    {
        return m_constants->slopeLimitGradients;
    }

    inline static double getSL_THETA()
    {
        return m_constants->SL_THETA;
    }

    inline static double getMESH_ETA()
    {
        return m_constants->MESH_ETA;
    }

    inline static double getMESH_CHI()
    {
        return m_constants->MESH_CHI;
    }

    inline static bool getMESH_ROUNDNESS_CORRECTION()
    {
        return m_constants->MESH_ROUNDNESS_CORRECTION;
    }

    inline static double getMESH_OMEGA()
    {
        return m_constants->MESH_OMEGA;
    }

    inline static double getMESH_PSI()
    {
        return m_constants->MESH_PSI;
    }

    inline static double getMASS_DENSITY_FLOOR()
    {
        return m_constants->MASS_DENSITY_FLOOR;
    }

    inline static double getENERGY_DENSITY_FLOOR()
    {
        return m_constants->ENERGY_DENSITY_FLOOR;
    }

    inline static double getPRESSURE_FLOOR()
    {
        return m_constants->PRESSURE_FLOOR;
    }

    inline static double getACCELERATION_X()
    {
        return m_constants->ACCELERATION_X;
    }

    inline static double getACCELERATION_Y()
    {
        return m_constants->ACCELERATION_Y;
    }

    inline static bool getENABLE_ACCELERATION()
    {
        return m_constants->ENABLE_ACCELERATION;
    }

    inline static unsigned int getMESH_REGULATION_FUNCTION()
    {
        return m_constants->MESH_REGULATION_FUNCTION;
    }

    inline static void (* function_regulateMesh())(Sim*)
    {
        return m_constants->regulateMesh;
    }

    inline static void (* function_determine_global_timestep())(Sim*)
    {
        return m_constants->determine_global_timestep;
    }

    inline static void (* function_delete_zero_length_edges())(Sim*)
    {
        return m_constants->delete_zero_length_edges;
    }

    inline static void (* function_roundness_correction())(Sim*)
    {
        return m_constants->roundness_correction;
    }

    inline static bool getOBSTACLE()
    {
        return m_constants->OBSTACLE;
    }

    inline static double getRADIUS()
    {
        return m_constants->RADIUS;
    }

    inline static double getCIRCLE_DR()
    {
        return m_constants->CIRCLE_DR;
    }

    inline static double getCIRCLE_DPHI()
    {
        return m_constants->CIRCLE_DPHI;
    }

    inline static double getCIRCLE_X()
    {
        return m_constants->CIRCLE_X;
    }

    inline static double getCIRCLE_Y()
    {
        return m_constants->CIRCLE_Y;
    }

    inline static double getOBS_VX()
    {
        return m_constants->OBS_VX;
    }

    inline static double getOBS_VY()
    {
        return m_constants->OBS_VY;
    }

    inline static double getOBS_AX()
    {
        return m_constants->OBS_AX;
    }

    inline static double getOBS_AY()
    {
        return m_constants->OBS_AY;
    }

    inline static double getOBS_MASS()
    {
        return m_constants->OBS_MASS;
    }

    inline static bool getOBS_FLOATING()
    {
        return m_constants->OBS_FLOATING;
    }

    // setters
    inline static void setSEED(unsigned int seed)
    {
        m_constants->SEED = seed;
    }

private:

    // constants from file

    // filename for infofile
    char* INFO_FILENAME;

    // filename for warnings
    char* WARNINGS_FILENAME;

    // random seed
    unsigned int SEED;

    // box type
    unsigned int BOX_TYPE;

    // box dimensions
    double BOX_XMIN;
    double BOX_YMIN;
    double BOX_XMAX;
    double BOX_YMAX;

    // number of particles
    unsigned int NOF_PARTICLES_X;
    unsigned int NOF_PARTICLES_Y;

    // relax the first graph?
    bool RELAX_FIRST_GRAPH;

    // delete zero length edges?
    bool DEL_ZERO_EDGES;

    // time step size
    double DELTA_T;

    // time interval between binary outputs
    double DELTA_T_OUTPUT;

    // binary output every nth step
    int EVERY_NTH_STEPS_OUTPUT;

    // maximum simulation time
    double MAX_T;

    // maximum simulation time steps
    int MAX_STEPS;

    // adiabatic exponent
    double GAMMA;

    // CFL coefficient
    double CFL;

    // velocity update method
    unsigned int VELOCITY_UPDATE_FUNCTION;

    // slope limiter function
    unsigned int SLOPE_LIMITER;

    // parameter for the tess slope limiter
    double SL_THETA;

    // parameters for Springel1 mesh regulation
    bool MESH_ROUNDNESS_CORRECTION;
    double MESH_ETA;
    double MESH_CHI;

    // parameters for havoc mesh regulation
    double MESH_OMEGA;
    double MESH_PSI;

    // mesh regulation function
    unsigned int MESH_REGULATION_FUNCTION;

    // adaptive mesh refinement
    bool ADAPTIVE_MESH_REFINEMENT;
    double MAX_CELL_VOLUME;

    // floors
    double MASS_DENSITY_FLOOR;
    double ENERGY_DENSITY_FLOOR;
    double PRESSURE_FLOOR;

    // acceleration
    bool ENABLE_ACCELERATION;
    double ACCELERATION_X;
    double ACCELERATION_Y;

    // obstacle
    bool OBSTACLE;
    double RADIUS;

    // circle
    double CIRCLE_DR;
    double CIRCLE_DPHI;

    // position
    double CIRCLE_X;
    double CIRCLE_Y;

    // movement of the obstacle
    double OBS_VX;
    double OBS_VY;
    double OBS_AX;
    double OBS_AY;

    // mass
    double OBS_MASS;

    // floating
    bool OBS_FLOATING;

    // further constants

    // stream for the info file
    FILE* P_INFO_FILE;

    // stream for the warnings file
    FILE* P_WARNINGS_FILE;

    // function pointer: update velocities of the particles
    void (* updateParticleVelocities)(Sim*);

    // function pointer: slope limit the gradients
    void (* slopeLimitGradients)(Sim*);

    // function pointer: mesh regulation
    void (* regulateMesh)(Sim*);

    // function pointer: global time step
    void (* determine_global_timestep)(Sim*);

    // function pointer: delete zero length edges
    void (* delete_zero_length_edges)(Sim*);

    // function pointer: roundness correction
    void (* roundness_correction)(Sim*);
};

#endif