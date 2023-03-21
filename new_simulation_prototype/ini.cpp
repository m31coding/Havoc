#include "sim_ini.h"
#include "sim.h"
#include "sim_core.h"

void sim_ini::initialize(Sim* sim)
{
    //sim->primVarIni3States(1, 0, 3./5, -0.4, 1+0.000001, 1, 3./5, -0.35, 1, 0, 3./5);
    //sim->primVarIni1State(1,2*sqrt(Constants::getGAMMA()),0,1); // density, velocity x, velocity y, pressure
    //sim_ini::two_colors(1,0,2, sim); // color top, broder, color bottom
    //sim->primVarIni1State(1,0,0,0.1); // density, velocity x, velocity y, pressure // SEDOV
    //sim->primVarIni1State(1,0,0,1); // density, velocity x, velocity y, pressure
    //sim->primVarIni2States(1.,0,1.,0,2,0,1); // density, velocity, pressure, border,...
    //sim->primVarIni2States_vertical(1.,0,1.,0,2,0,1);
    //sim->primVarIni3States(1,0.1,3./5,-1,1.5,0.1,3./5,1,1,0.1,3./5);

    // hydrostatic equilibrium
    //sim_ini::set_pressure_hydrostatic_equilibrium(sim, 2, 1);

    //sim_ini::prim_var_ini_quad_shock(sim); // QUAD SHOCK
    //sim_ini::primVarIniRT(sim,0); // RAYLEIGH TAYLOR
    //sim_ini::primVarIniRT_arepo_square(sim); // RAYLEIGH TAYLOR

    // KELVIN HELMHOLTZ
    // d v_x p border d v_x p border d v_x p
    sim->primVarIni3States_vertical(1, -0.5, 2.5, 0.25, 2, 0.5, 2.5, -0.25, 1, -0.5, 2.5);

    // perturb the y-velocities in order to trigger a single mode
    sim_ini::perturbing_v_y_arepo(0.25, -0.25, 0.1, 0.05 / sqrt(2), sim);

    // calculate the conserved fluid variables
    sim_core::updateCellVolumesAndCenters(sim); // only real particles
    sim_core::calculateConservedVariables(sim); // only real particles

    // SEDOV inject energy
    //sim->injectEnergyGauss(1, 0.03, 0., 0.);
    //sim->injectEnergyMiddleParticle(1);

    // add a constant boost velocity
    //sim_core::calculatePrimitiveVariables(sim);
    //sim_ini::x_boost(sim, 10);
    //sim_core::calculateConservedVariables(sim);
}