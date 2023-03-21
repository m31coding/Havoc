#include "config.h"
#include "constants.h"
#include "debug.h"
#include "read_constants.h"
#include <cstdlib>
#include "sim_core.h"
#include <limits>
#include <iostream>

// initialize m_constants
Constants* Constants::m_constants = NULL;

// read the constants from file "constants"
void Constants::readConstants()
{
    // create constants
    m_constants = new Constants;

    ReadConstants::readConstants("constants");

    if (!ReadConstants::keyExists("INFO_FILENAME"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key INFO_FILENAME does not exist\n");
        exit(1);
    }

    Constants::get()->INFO_FILENAME = new char[ReadConstants::GetString("INFO_FILENAME").length() + 1];
    strcpy(Constants::get()->INFO_FILENAME, ReadConstants::GetString("INFO_FILENAME").c_str());

    if (!ReadConstants::keyExists("WARNINGS_FILENAME"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key WARNINGS_FILENAME does not exist\n");
        exit(1);
    }

    Constants::get()->WARNINGS_FILENAME = new char[ReadConstants::GetString("WARNINGS_FILENAME").length() + 1];
    strcpy(Constants::get()->WARNINGS_FILENAME, ReadConstants::GetString("WARNINGS_FILENAME").c_str());

    if (!ReadConstants::keyExists("SEED"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key SEED does not exist\n");
        exit(1);
    }
    Constants::get()->SEED = ReadConstants::GetUnsignedInt("SEED");

    if (!ReadConstants::keyExists("BOX_TYPE"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key BOX_TYPE does not exist\n");
        exit(1);
    }
    Constants::get()->BOX_TYPE = ReadConstants::GetUnsignedInt("BOX_TYPE");

    if (!ReadConstants::keyExists("BOX_XMIN"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key BOX_XMIN does not exist\n");
        exit(1);
    }
    Constants::get()->BOX_XMIN = ReadConstants::GetDouble("BOX_XMIN");

    if (!ReadConstants::keyExists("BOX_XMAX"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key BOX_XMAX does not exist\n");
        exit(1);
    }
    Constants::get()->BOX_XMAX = ReadConstants::GetDouble("BOX_XMAX");

    if (!ReadConstants::keyExists("BOX_YMIN"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key BOX_YMIN does not exist\n");
        exit(1);
    }
    Constants::get()->BOX_YMIN = ReadConstants::GetDouble("BOX_YMIN");

    if (!ReadConstants::keyExists("BOX_YMAX"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key BOX_YMAX does not exist\n");
        exit(1);
    }
    Constants::get()->BOX_YMAX = ReadConstants::GetDouble("BOX_YMAX");

    if (!ReadConstants::keyExists("NOF_PARTICLES_X"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key NOF_PARTICLES_X does not exist\n");
        exit(1);
    }
    Constants::get()->NOF_PARTICLES_X = ReadConstants::GetUnsignedInt("NOF_PARTICLES_X");

    if (!ReadConstants::keyExists("NOF_PARTICLES_Y"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key NOF_PARTICLES_Y does not exist\n");
        exit(1);
    }
    Constants::get()->NOF_PARTICLES_Y = ReadConstants::GetUnsignedInt("NOF_PARTICLES_Y");

    if (!ReadConstants::keyExists("RELAX_FIRST_GRAPH"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key RELAX_FIRST_GRAPH does not exist\n");
        exit(1);
    }
    Constants::get()->RELAX_FIRST_GRAPH = ReadConstants::GetBool("RELAX_FIRST_GRAPH");

    if (!ReadConstants::keyExists("DEL_ZERO_EDGES"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key DEL_ZERO_EDGES does not exist\n");
        exit(1);
    }
    Constants::get()->DEL_ZERO_EDGES = ReadConstants::GetBool("DEL_ZERO_EDGES");

    if (!ReadConstants::keyExists("DELTA_T_OUTPUT"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key DELTA_T_OUTPUT does not exist\n");
        exit(1);
    }
    Constants::get()->DELTA_T_OUTPUT = ReadConstants::GetDouble("DELTA_T_OUTPUT");

    if (!ReadConstants::keyExists("EVERY_NTH_STEPS_OUTPUT"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key EVERY_NTH_STEPS_OUTPUT does not exist\n");
        exit(1);
    }
    Constants::get()->EVERY_NTH_STEPS_OUTPUT = ReadConstants::GetInt("EVERY_NTH_STEPS_OUTPUT");

    if (!ReadConstants::keyExists("DELTA_T"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key DELTA_T does not exist\n");
        exit(1);
    }
    Constants::get()->DELTA_T = ReadConstants::GetDouble("DELTA_T");

    if (!ReadConstants::keyExists("MAX_T"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key MAX_T does not exist\n");
        exit(1);
    }
    Constants::get()->MAX_T = ReadConstants::GetDouble("MAX_T");

    if (!ReadConstants::keyExists("MAX_STEPS"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key MAX_STEPS does not exist\n");
        exit(1);
    }
    Constants::get()->MAX_STEPS = ReadConstants::GetInt("MAX_STEPS");

    if (!ReadConstants::keyExists("GAMMA"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key GAMMA does not exist\n");
        exit(1);
    }
    Constants::get()->GAMMA = ReadConstants::GetDouble("GAMMA");

    if (!ReadConstants::keyExists("CFL"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key CFL does not exist\n");
        exit(1);
    }
    Constants::get()->CFL = ReadConstants::GetDouble("CFL");

    if (!ReadConstants::keyExists("VELOCITY_UPDATE_FUNCTION"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key VELOCITY_UPDATE_FUNCTION does not exist\n");
        exit(1);
    }
    Constants::get()->VELOCITY_UPDATE_FUNCTION = ReadConstants::GetUnsignedInt("VELOCITY_UPDATE_FUNCTION");

    if (!ReadConstants::keyExists("SLOPE_LIMITER"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key SLOPE_LIMITER does not exist\n");
        exit(1);
    }
    Constants::get()->SLOPE_LIMITER = ReadConstants::GetUnsignedInt("SLOPE_LIMITER");

    if (!ReadConstants::keyExists("SL_THETA"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key SL_THETA does not exist\n");
        exit(1);
    }
    Constants::get()->SL_THETA = ReadConstants::GetDouble("SL_THETA");

    if (!ReadConstants::keyExists("MESH_ETA"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key MESH_ETA does not exist\n");
        exit(1);
    }
    Constants::get()->MESH_ETA = ReadConstants::GetDouble("MESH_ETA");


    if (!ReadConstants::keyExists("MESH_CHI"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key MESH_CHI does not exist\n");
        exit(1);
    }
    Constants::get()->MESH_CHI = ReadConstants::GetDouble("MESH_CHI");


    if (!ReadConstants::keyExists("MESH_ROUNDNESS_CORRECTION"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key MESH_ROUNDNESS_CORRECTION does not exist\n");
        exit(1);
    }
    Constants::get()->MESH_ROUNDNESS_CORRECTION = ReadConstants::GetBool("MESH_ROUNDNESS_CORRECTION");

    if (!ReadConstants::keyExists("MESH_OMEGA"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key MESH_OMEGA does not exist\n");
        exit(1);
    }
    Constants::get()->MESH_OMEGA = ReadConstants::GetDouble("MESH_OMEGA");

    if (!ReadConstants::keyExists("MESH_PSI"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key MESH_PSI does not exist\n");
        exit(1);
    }
    Constants::get()->MESH_PSI = ReadConstants::GetDouble("MESH_PSI");

    if (!ReadConstants::keyExists("MESH_REGULATION_FUNCTION"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key MESH_REGULATION_FUNCTION does not exist\n");
        exit(1);
    }
    Constants::get()->MESH_REGULATION_FUNCTION = ReadConstants::GetUnsignedInt("MESH_REGULATION_FUNCTION");

    if (!ReadConstants::keyExists("MASS_DENSITY_FLOOR"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key MASS_DENSITY_FLOOR does not exist\n");
        exit(1);
    }
    Constants::get()->MASS_DENSITY_FLOOR = ReadConstants::GetDouble("MASS_DENSITY_FLOOR");

    if (!ReadConstants::keyExists("ENERGY_DENSITY_FLOOR"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key ENERGY_DENSITY_FLOOR does not exist\n");
        exit(1);
    }
    Constants::get()->ENERGY_DENSITY_FLOOR = ReadConstants::GetDouble("ENERGY_DENSITY_FLOOR");

    if (!ReadConstants::keyExists("PRESSURE_FLOOR"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key PRESSURE_FLOOR does not exist\n");
        exit(1);
    }
    Constants::get()->PRESSURE_FLOOR = ReadConstants::GetDouble("PRESSURE_FLOOR");

    if (!ReadConstants::keyExists("ENABLE_ACCELERATION"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key ENABLE_ACCELERATION does not exist\n");
        exit(1);
    }
    Constants::get()->ENABLE_ACCELERATION = ReadConstants::GetBool("ENABLE_ACCELERATION");

    if (!ReadConstants::keyExists("ACCELERATION_X"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key ACCELERATION_X does not exist\n");
        exit(1);
    }
    Constants::get()->ACCELERATION_X = ReadConstants::GetDouble("ACCELERATION_X");

    if (!ReadConstants::keyExists("ACCELERATION_Y"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key ACCELERATION_Y does not exist\n");
        exit(1);
    }
    Constants::get()->ACCELERATION_Y = ReadConstants::GetDouble("ACCELERATION_Y");

    if (!ReadConstants::keyExists("OBSTACLE"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key OBSTACLE does not exist\n");
        exit(1);
    }
    Constants::get()->OBSTACLE = ReadConstants::GetBool("OBSTACLE");

    if (!ReadConstants::keyExists("CIRCLE_DR"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key CIRCLE_DR does not exist\n");
        Constants::get()->CIRCLE_DR = 0.1;
    }
    else
    {
        Constants::get()->CIRCLE_DR = ReadConstants::GetDouble("CIRCLE_DR");
    }

    if (!ReadConstants::keyExists("CIRCLE_DPHI"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key CIRCLE_DPHI does not exist\n");
        Constants::get()->CIRCLE_DPHI = 0.5;
    }
    else
    {
        Constants::get()->CIRCLE_DPHI = ReadConstants::GetDouble("CIRCLE_DPHI");
    }

    if (!ReadConstants::keyExists("CIRCLE_X"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key CIRCLE_X does not exist\n");
        Constants::get()->CIRCLE_X = 0.;
    }
    else
    {
        Constants::get()->CIRCLE_X = ReadConstants::GetDouble("CIRCLE_X");
    }

    if (!ReadConstants::keyExists("CIRCLE_Y"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key CIRCLE_Y does not exist\n");
        Constants::get()->CIRCLE_Y = 0.;
    }
    else
    {
        Constants::get()->CIRCLE_Y = ReadConstants::GetDouble("CIRCLE_Y");
    }

    if (!ReadConstants::keyExists("OBS_VX"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key OBS_VX does not exist\n");
        Constants::get()->OBS_VX = 0;
    }
    else
    {
        Constants::get()->OBS_VX = ReadConstants::GetDouble("OBS_VX");
    }

    if (!ReadConstants::keyExists("OBS_VY"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key OBS_VY does not exist\n");
        Constants::get()->OBS_VY = 0;
    }
    else
    {
        Constants::get()->OBS_VY = ReadConstants::GetDouble("OBS_VY");
    }

    if (!ReadConstants::keyExists("OBS_AX"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key OBS_AX does not exist\n");
        Constants::get()->OBS_AX = 0;
    }
    else
    {
        Constants::get()->OBS_AX = ReadConstants::GetDouble("OBS_AX");
    }

    if (!ReadConstants::keyExists("OBS_AY"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key OBS_AY does not exist\n");
        Constants::get()->OBS_AY = 0;
    }
    else
    {
        Constants::get()->OBS_AY = ReadConstants::GetDouble("OBS_AY");
    }

    if (!ReadConstants::keyExists("RADIUS"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key RADIUS does not exist\n");
        exit(1);
    }
    Constants::get()->RADIUS = ReadConstants::GetDouble("RADIUS");

    if (!ReadConstants::keyExists("OBS_MASS"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key OBS_MASS does not exist\n");
        Constants::get()->OBS_MASS = 1;
    }
    else
    {
        Constants::get()->OBS_MASS = ReadConstants::GetDouble("OBS_MASS");
    }

    if (!ReadConstants::keyExists("OBS_FLOATING"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key OBS_FLOATING does not exist\n");
        Constants::get()->OBS_FLOATING = false;
    }
    else
    {
        Constants::get()->OBS_FLOATING = ReadConstants::GetBool("OBS_FLOATING");
    }

    if (!ReadConstants::keyExists("ADAPTIVE_MESH_REFINEMENT"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key ADAPTIVE_MESH_REFINEMENT does not exist\n");
        Constants::get()->ADAPTIVE_MESH_REFINEMENT = false;
    }
    else
    {
        Constants::get()->ADAPTIVE_MESH_REFINEMENT = ReadConstants::GetBool("ADAPTIVE_MESH_REFINEMENT");
    }

    if (!ReadConstants::keyExists("MAX_CELL_VOLUME"))
    {
        fprintf(stderr, "ERROR in Constants::readConstants(): key MAX_CELL_VOLUME does not exist\n");
        Constants::get()->MAX_CELL_VOLUME = 1;
    }
    else
    {
        Constants::get()->MAX_CELL_VOLUME = ReadConstants::GetDouble("MAX_CELL_VOLUME");
    }

    // free memory
    ReadConstants::del();

    // create further constants
    createFurtherConstants();

    // print constants
#ifdef CONSTANTS_DEBUG
    fprint(stdout);
#endif

}

// print the constants to a file
void Constants::fprint(FILE* pFile)
{
    fprintf(pFile, "\n-constants-\n\n");

    fprintf(pFile, "INFO_FILENAME\t%s\n", Constants::getINFO_FILENAME());
    fprintf(pFile, "WARNINGS_FILENAME\t%s\n", Constants::getWARNINGS_FILENAME());
    fprintf(pFile, "SEED\t%u\n", Constants::getSEED());
    fprintf(pFile, "BOX_TYPE\t%u\n", Constants::getBOX_TYPE());
    fprintf(pFile, "BOX_XMIN\t%f\n", Constants::getBOX_XMIN());
    fprintf(pFile, "BOX_YMIN\t%f\n", Constants::getBOX_YMIN());
    fprintf(pFile, "BOX_XMAX\t%f\n", Constants::getBOX_XMAX());
    fprintf(pFile, "BOX_YMAX\t%f\n", Constants::getBOX_YMAX());
    fprintf(pFile, "NOF_PARTICLES_X\t%u\n", Constants::getNOF_PARTICLES_X());
    fprintf(pFile, "NOF_PARTICLES_Y\t%u\n", Constants::getNOF_PARTICLES_Y());
    fprintf(pFile, "RELAX_FIRST_GRAPH\t%s\n", (Constants::getRELAX_FIRST_GRAPH()) ? "yes" : "no");
    fprintf(pFile, "DELTA_T\t%f\n", Constants::getDELTA_T());
    fprintf(pFile, "DELTA_T_OUTPUT\t%f\n", Constants::getDELTA_T_OUTPUT());
    fprintf(pFile, "EVERY_STEPS_OUTPUT\t%d\n", Constants::getEVERY_NTH_STEPS_OUTPUT());
    fprintf(pFile, "MAX_T\t%f\n", Constants::getMAX_T());
    fprintf(pFile, "MAX_STEPS\t%d\n", Constants::getMAX_STEPS());
    fprintf(pFile, "GAMMA\t%f\n", Constants::getGAMMA());
    fprintf(pFile, "CFL\t%f\n", Constants::getCFL());
    fprintf(pFile, "VELOCITY_UPDATE_FUNCTION\t%u\n", Constants::getVELOCITY_UPDATE_FUNCTION());
    fprintf(pFile, "SLOPE_LIMITER\t%u\n", Constants::getSLOPE_LIMITER());
    fprintf(pFile, "SL_THETA\t%f\n", Constants::getSL_THETA());
    fprintf(pFile, "MESH_ETA\t%f\n", Constants::getMESH_ETA());
    fprintf(pFile, "MESH_CHI\t%f\n", Constants::getMESH_CHI());
    fprintf(pFile, "MESH_ROUNDNESS_CORRECTION\t%s\n", (Constants::getMESH_ROUNDNESS_CORRECTION()) ? "yes" : "no");
    fprintf(pFile, "MESH_OMEGA\t%f\n", Constants::getMESH_OMEGA());
    fprintf(pFile, "MESH_PSI\t%f\n", Constants::getMESH_PSI());
    fprintf(pFile, "MESH_REGULATION_FUNCTION\t%u\n", Constants::getMESH_REGULATION_FUNCTION());
    fprintf(pFile, "ADAPTIVE_MESH_REFINEMENT\t%s\n", Constants::getADAPTIVE_MESH_REFINEMENT() ? "yes" : "no");
    fprintf(pFile, "MAX_CELL_VOLUME\t%f\n", Constants::getMAX_CELL_VOLUME());
    fprintf(pFile, "DEL_ZERO_EDGES\t%s\n", (Constants::getDEL_ZERO_EDGES()) ? "yes" : "no");
    fprintf(pFile, "MASS_DENSITY_FLOOR\t%f\n", Constants::getMASS_DENSITY_FLOOR());
    fprintf(pFile, "ENERGY_DENSITY_FLOOR\t%f\n", Constants::getENERGY_DENSITY_FLOOR());
    fprintf(pFile, "PRESSURE_FLOOR\t%f\n", Constants::getPRESSURE_FLOOR());
    fprintf(pFile, "ENABLE_ACCELERATION\t%s\n", (Constants::getENABLE_ACCELERATION()) ? "yes" : "no");
    fprintf(pFile, "ACCELERATION_X\t%f\n", Constants::getACCELERATION_X());
    fprintf(pFile, "ACCELERATION_Y\t%f\n", Constants::getACCELERATION_Y());
    fprintf(pFile, "OBSTACLE\t%s\n", (Constants::getDEL_ZERO_EDGES()) ? "yes" : "no");
    fprintf(pFile, "RADIUS\t%f\n", Constants::getRADIUS());
    fprintf(pFile, "CIRCLE_X\t%f\n", Constants::getCIRCLE_X());
    fprintf(pFile, "CIRCLE_Y\t%f\n", Constants::getCIRCLE_Y());
    fprintf(pFile, "CIRCLE_DR\t%f\n", Constants::getCIRCLE_DR());
    fprintf(pFile, "CIRCLE_DPHI\t%f\n", Constants::getCIRCLE_DPHI());
    fprintf(pFile, "OBS_VX\t%f\n", Constants::getOBS_VX());
    fprintf(pFile, "OBS_VY\t%f\n", Constants::getOBS_VY());
    fprintf(pFile, "OBS_AX\t%f\n", Constants::getOBS_AX());
    fprintf(pFile, "OBS_AY\t%f\n", Constants::getOBS_AY());
    fprintf(pFile, "OBS_MASS\t%f\n", Constants::getOBS_MASS());
    fprintf(pFile, "OBS_FLOATING\t%s\n", (Constants::getOBS_FLOATING()) ? "yes" : "no");
}

// create further constants
void Constants::createFurtherConstants()
{
    SD_PRINTF("creating further constants\n");

    m_constants->P_INFO_FILE = fopen(m_constants->INFO_FILENAME, "a");
    if (m_constants->P_INFO_FILE == NULL)
    {
        fprintf(stderr, "ERROR in Constants::createFurtherConstants(): can't open file %s\n",
                m_constants->INFO_FILENAME);
        exit(1);
    }

    m_constants->P_WARNINGS_FILE = fopen(m_constants->WARNINGS_FILENAME, "w");
    if (m_constants->P_WARNINGS_FILE == NULL)
    {
        fprintf(stderr, "ERROR in Constants::createFurtherConstants(): can't open file %s\n",
                m_constants->WARNINGS_FILENAME);
        exit(1);
    }

    // set function pointer of updateParticleVelocities
    switch (m_constants->VELOCITY_UPDATE_FUNCTION)
    {
        case 0:
        {
            m_constants->updateParticleVelocities = &sim_core::doNothing;
            FP_PRINTF("updateParticleVelocities = doNothing\n");
        }
            break;

        case 1:
        {
            m_constants->updateParticleVelocities = &sim_core::updateParticleVelocitiesLagrangian;
            FP_PRINTF("updateParticleVelocities = updateParticleVelocitiesLagrangian\n");
        }
            break;

        default:
        {
            fprintf(stderr, "ERROR in Constants::createFurtherConstants: invalid constant VELOCITY_UPDATE_FUNCTION\n");
            exit(1);
        }
    }

    // set function pointer of the slope limiter
    switch (m_constants->SLOPE_LIMITER)
    {
        case 0:
        {
            m_constants->slopeLimitGradients = &sim_core::doNothing;
            FP_PRINTF("slopeLimitGradients = doNothing\n");
        }
            break;

        case 1:
        {
            m_constants->slopeLimitGradients = &sim_core::slope_limit_gradients_arepo;
            FP_PRINTF("slopeLimitGradients = slope_limit_gradients_arepo\n");
        }
            break;

        case 2:
        {
            m_constants->slopeLimitGradients = &sim_core::slope_limit_gradients_tess;
            FP_PRINTF("slopeLimitGradients = slope_limit_gradients_tess\n");
        }
            break;

        default:
        {
            fprintf(stderr, "ERROR in Constants::createFurtherConstants: invalid constant SLOPE_LIMITER\n");
            exit(1);
        }
    }

    // set function pointer of the mesh regulation
    switch (m_constants->MESH_REGULATION_FUNCTION)
    {
        case 0:
        {
            m_constants->regulateMesh = &sim_core::doNothing;
            FP_PRINTF("regulateMesh = doNothing\n");
        }
            break;

        case 1:
        {
            m_constants->regulateMesh = &sim_core::do_mesh_regulation_arepo_1;
            FP_PRINTF("regulateMesh = calculateMeshRegularitySpringel1\n");
        }
            break;

        case 2:
        {
            m_constants->regulateMesh = &sim_core::do_mesh_regulation_arepo_2;
            FP_PRINTF("regulateMesh = calculateMeshRegularitySpringel1\n");
        }
            break;

        default:
        {
            fprintf(stderr, "ERROR in Constants::createFurtherConstants: invalid constant MESH_REGULATION_FUNCTION\n");
            exit(1);
        }
    }

    // set function pointer of the global time step calculation
    if (m_constants->DELTA_T == -1)
    {
        m_constants->determine_global_timestep = &sim_core::variable_global_timestep;
    }
    else
    {
        m_constants->determine_global_timestep = &sim_core::constant_global_timestep;
    }

    // set function pointer for havoc mesh correction
    if (m_constants->MESH_ROUNDNESS_CORRECTION == true)
    {
        m_constants->roundness_correction = &sim_core::roundness_correction;
    }
    else
    {
        m_constants->roundness_correction = &sim_core::doNothing;
    }

    // set function pointer for deleting zero length edges
    if (m_constants->DEL_ZERO_EDGES == true)
    {
        m_constants->delete_zero_length_edges = &sim_core::delete_zero_length_edges;
    }
    else
    {
        m_constants->delete_zero_length_edges = &sim_core::doNothing;
    }
}

// check and adapt constants
void Constants::check_and_adapt()
{
    if (m_constants->MAX_T == -1)
    {
        m_constants->MAX_T = std::numeric_limits<double>::max();
    }
    else if (m_constants->MAX_T < 0)
    {
        fprintf(stderr, "ERROR in Constants::check_and_adapt(): invalid constant value: MAX_T\t%f\n",
                m_constants->MAX_T);
        exit(1);
    }

    if (m_constants->MAX_STEPS == -1)
    {
        m_constants->MAX_STEPS = std::numeric_limits<unsigned int>::max();
    }
    else if (m_constants->MAX_STEPS < 0)
    {
        fprintf(stderr, "ERROR in Constants::check_and_adapt(): invalid constant value: MAX_STEPS\t%d\n",
                m_constants->MAX_STEPS);
        exit(1);
    }

    if (m_constants->DELTA_T_OUTPUT <= 0)
    {
        fprintf(stderr, "ERROR in Constants::check_and_adapt(): invalid constant value: DELTA_T_OUTPUT\t%f\n",
                m_constants->DELTA_T_OUTPUT);
        exit(1);
    }

    if (m_constants->VELOCITY_UPDATE_FUNCTION == 0 &&
        (m_constants->MESH_REGULATION_FUNCTION == 1 || m_constants->MESH_ROUNDNESS_CORRECTION == true))
    {
        fprintf(stderr,
                "ERROR in Constants::check_and_adapt(): velocity update function is 0 but mesh regulation function is 1\n");
        exit(1);
    }

    if (m_constants->DELTA_T != -1 && m_constants->DELTA_T <= 0)
    {
        fprintf(stderr, "ERROR in Constants::check_and_adapt(): invalid constant value: DELTA_T\t%f\n",
                m_constants->DELTA_T);
        exit(1);
    }

    if (m_constants->DELTA_T != -1 && m_constants->DELTA_T_OUTPUT / m_constants->DELTA_T -
                                      (int) (m_constants->DELTA_T_OUTPUT / m_constants->DELTA_T) != 0.)
    {
        fprintf(stderr, "ERROR in Constants::check_and_adapt(): DELTA_T_OUTPUT is no multiple of DELTA_T\n");
        exit(1);
    }

    if (m_constants->MASS_DENSITY_FLOOR != -1 && m_constants->MASS_DENSITY_FLOOR <= 0)
    {
        fprintf(stderr, "ERROR in Constants::check_and_adapt(): MASS_DENSITY_FLOOR has to be positive\n");
    }

    if (m_constants->ENERGY_DENSITY_FLOOR != -1 && m_constants->ENERGY_DENSITY_FLOOR <= 0)
    {
        fprintf(stderr, "ERROR in Constants::check_and_adapt(): ENERGY_DENSITY_FLOOR has to be positive\n");
    }

    if (m_constants->PRESSURE_FLOOR != -1 && m_constants->PRESSURE_FLOOR < 0)
    {
        fprintf(stderr, "ERROR in Constants::check_and_adapt(): PRESSURE_FLOOR has to be non-negative\n");
    }

    if (m_constants->RADIUS <= 0)
    {
        fprintf(stderr, "ERROR in Constants::check_and_adapt(): RADIUS has to be positive\n");
    }
}