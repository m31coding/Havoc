#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H

/**
 * Implementation of an exact Riemann solver following the book
 * Riemann Solvers and Numerical Methods for Fluid Dynamics (Toro) chapter 4 and 6
 * Many thanks to Prof. Eleuterio Toro for his kind permission to include this implementation in the repository.
 */
namespace exact_riemann_solver
{
    /// initialize constants (call this function once)
    void ini();

    /// calculate variables
    void calculateVariables(double rho_left, double u_left, double p_left,
                            double rho_right, double u_right, double p_right);

    /// clear variables (for purpose of debug)
    void resetVariables();

    /// calculate pressure and velocity in the star region
    void starPressureAndVelocity(double* starPressure, double* starVelocity);

    /// calculate start pressure for Newton-Raphson (linearized solution)
    double guessPressure();

    /**
	sample the solution at speed S	
	\param D the sampled density (output)
	\param U the sampled velocity (output)
	\param P the sampled pressure (output)
	\param PM the pressure in the star region
	\param UM the velocity in the star region
	\param S the speed s = x / t
    */
    void sample(double* D, double* U, double* P, double PM, double UM, double S);

    /**
	solve the Riemann problem with an exact solver and sample the solution along the speed s = x / t = 0
	\param rho_middle the sampled density (output)
	\param u_middle the sampled velocity (output)
	\param p_middle the sampled pressure (output)
	\param rho_left the density of the left side of the Riemann problem
	\param u_left the velocity of the left side
	\param p_left the pressure of the left side
	\param rho_right the density of the right side of the Riemann problem
	\param u_right the velocity of the right side
	\param p_right the pressure of the right side
    */
    void solveAndSampleZeroSpeed(double* rho_middle, double* u_middle, double* p_middle, double rho_left, double u_left,
                                 double p_left, double rho_right, double u_right, double p_right);

};

#endif