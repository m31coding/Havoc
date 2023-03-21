#ifndef SIM_INI_H
#define SIM_INI_H

class Sim;

namespace sim_ini
{
    /// sets the initial conditions
    void initialize(Sim* sim);

    /// color the fluid with two different colors
    void two_colors(double c_top, double border, double c_bottom, Sim* sim);

    /// perturbation function for Kelvin-Helmholtz
    void perturbing_v_y_arepo(double y1, double y2, double w_0, double sigma, Sim* sim);

    /// initialize the primitive variables for the Rayleigh-Taylor instability (Arepo)
    void primVarIniRT_arepo(Sim* sim); //Rayleigh-Taylor instability

    /// initialize the primitive variables for the Rayleigh-Taylor instability (Arepo, square box)
    void primVarIniRT_arepo_square(Sim* sim); //Rayleigh-Taylor instability

    /// initialize the primitive variables for the Rayleigh-Taylor instability
    void primVarIniRT(Sim* sim, double border); //Rayleigh-Taylor instability

    /// initialize the primitive variables for a quad shock simulation
    void prim_var_ini_quad_shock(Sim* sim);

    /// add a boost velocity in x-direction
    void x_boost(Sim* sim, double v_x);

    /// add a boost velocity in y-direction
    void y_boost(Sim* sim, double v_y);

    /// the pressure is set according to the hydrostatic equilibrium
    void set_pressure_hydrostatic_equilibrium(Sim* sim, double P_0, double height_P_0);
};

#endif