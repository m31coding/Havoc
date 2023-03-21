#include "config.h"
#include "riemann_solver.h"
#include "constants.h"
#include "debug.h"
#include <stdlib.h>
#include <cmath>

/*
 * Implementation of an exact Riemann solver following the book
 * Riemann Solvers and Numerical Methods for Fluid Dynamics (Toro) chapter 4 and 6
 * Many thanks to Prof. Eleuterio Toro for his kind permission to include this implementation in the repository.
 */

static double TOL = pow(10, -6); // tolerance
unsigned int MAX_ITERATIONS = 20; // maximum number of iterations

// gamma related constants
static double GAMMA = 0;
static double G1 = 0;
static double G2 = 0;
static double G3 = 0;
static double G4 = 0;
static double G5 = 0;
static double G6 = 0;
static double G7 = 0;
static double G8 = 0;

// variables
static double CL = 0; // speed of sound of the left state
static double CR = 0; // speed of sound of the right state
static double DL = 0; // initial density of the left state
static double UL = 0; // initial velocity of the left state
static double PL = 0; // initial pressure of the left state
static double DR = 0; // initial density of the right state
static double UR = 0; // initial velocity of the right state
static double PR = 0; // initial pressure of the right state

// initialize constants
void exact_riemann_solver::ini()
{
    RS_PRINTF("\n---initializing Riemann solver---\n\n");

    if (Constants::getGAMMA() == 0 || Constants::getGAMMA() == 1)
    {
        fprintf(stderr, "ERROR in exact_riemann_solver::ini: invalid constant GAMMA, denominator is zero\n");
        exit(1);
    }

    GAMMA = Constants::getGAMMA();
    G1 = (Constants::getGAMMA() - 1.) * 0.5 / Constants::getGAMMA();
    G2 = (Constants::getGAMMA() + 1.) * 0.5 / Constants::getGAMMA();
    G3 = 2. * Constants::getGAMMA() / (Constants::getGAMMA() - 1.);
    G4 = 2. / (Constants::getGAMMA() - 1.);
    G5 = 2. / (Constants::getGAMMA() + 1.);
    G6 = (Constants::getGAMMA() - 1.) / (Constants::getGAMMA() + 1.);
    G7 = (Constants::getGAMMA() - 1.) / 2.;
    G8 = Constants::getGAMMA() - 1.;

    RS_PRINTF("TOL = %f\n", TOL);
    RS_PRINTF("MAX_ITERATIONS = %u\n", MAX_ITERATIONS);
    RS_PRINTF("GAMMA = %f\n", GAMMA);
    RS_PRINTF("G1 = %f\n", G1);
    RS_PRINTF("G2 = %f\n", G2);
    RS_PRINTF("G3 = %f\n", G3);
    RS_PRINTF("G4 = %f\n", G4);
    RS_PRINTF("G5 = %f\n", G5);
    RS_PRINTF("G6 = %f\n", G6);
    RS_PRINTF("G7 = %f\n", G7);
    RS_PRINTF("G8 = %f\n", G8);
}

// calculate variables
void exact_riemann_solver::calculateVariables(double rho_left, double u_left, double p_left, double rho_right,
                                              double u_right, double p_right)
{
    RS_PRINTF("calculating Riemann solver variables\n\n");

    // don't allow zero or negative densities or negative pressures
    if (rho_left <= 0 || rho_right <= 0 || p_left < 0 || p_right < 0)
    {
        fprintf(stderr, "ERROR in exact_riemann_solver::calculateVariables: invalid states!\n");
        fprintf(stderr, "\tdensity left = %f, velocity left = %f, pressure left = %f\n", rho_left, u_left, p_left);
        fprintf(stderr, "\tdensity right = %f, velocity right = %f, pressure right = %f\n", rho_right, u_right,
                p_right);
        exit(1);
    }

    DL = rho_left;
    UL = u_left;
    PL = p_left;
    DR = rho_right;
    UR = u_right;
    PR = p_right;

    // calculate speeds of sound
    CL = sqrt(GAMMA * PL / DL);
    CR = sqrt(GAMMA * PR / DR);

    RS_PRINTF("DL = %f\n", DL);
    RS_PRINTF("UL = %f\n", UL);
    RS_PRINTF("PL = %f\n", PL);
    RS_PRINTF("DR = %f\n", DR);
    RS_PRINTF("UR = %f\n", UR);
    RS_PRINTF("PR = %f\n", PR);

    RS_PRINTF("\n");

    RS_PRINTF("CL = %f\n", CL);
    RS_PRINTF("CR = %f\n", CR);

    // only allow velocities near the speed of sound
}

// clear variables (for debug purpose)
void exact_riemann_solver::resetVariables()
{
    DL = 0;
    UL = 0;
    PL = 0;
    DR = 0;
    UR = 0;
    PR = 0;
    CL = 0;
    CR = 0;
}


// evaluate pressure functions (helper for starPressureAndVelocity)
static void pref(double* F, double* FD, double P, double DK, double PK, double CK)
{
    if (P <= PK) // rarefaction wave
    {
        double PRAT = P / PK;
        *F = G4 * CK * (pow(PRAT, G1) - 1);
        *FD = (1. / (DK * CK)) * pow(PRAT, -G2);
    }

    else // shock wave
    {
        double AK = G5 / DK;
        double BK = G6 * PK;
        double QRT = sqrt(AK / (BK + P));
        *F = (P - PK) * QRT;
        *FD = (1. - 0.5 * (P - PK) / (BK + P)) * QRT;
    }
}

// calculate pressure and velocity in the star region
void exact_riemann_solver::starPressureAndVelocity(double* starPressure, double* starVelocity)
{
    RS_PRINTF("\ncalculating the pressure in the star region\n");

    // test pressure positivity condition
    if (G4 * (CL + CR) <= (UR - UL))
    {
        fprintf(stderr, "ERROR in exact_riemann_solver::starPressure: initial data generated vacuum\n");
        fprintf(stderr, "\tdensity left = %f, velocity left = %f, pressure left = %f\n", DL, UL, PL);
        fprintf(stderr, "\tdensity right = %f, velocity right = %f, pressure right = %f\n", DR, UR, PR);
        fprintf(stderr, "\tsound speed left = %f, sound speed right = %f\n", CL, CR);
        exit(1);
    }

    double P_OLD = guessPressure();
    RS_PRINTF("\tguess pressure: %f\n\n", P_OLD);

    double P = 0;
    double CHANGE = 0;
    double U_DIFF = UR - UL;
    double FL = 0;
    double FR = 0;
    double FLD = 0;
    double FRD = 0;
    bool converged = false;

    unsigned int i = 0;

    for (; i < MAX_ITERATIONS; i++)
    {

        pref(&FL, &FLD, P_OLD, DL, PL, CL);
        pref(&FR, &FRD, P_OLD, DR, PR, CR);

        if (FLD + FRD == 0)
        {
            fprintf(stderr, "ERROR in exact_riemann_solver::starPressure(): denominator is zero\n");
            exit(1);
        }

        P = P_OLD - (FL + FR + U_DIFF) / (FLD + FRD);

        CHANGE = 2 * fabs((P - P_OLD) / (P + P_OLD));

        if (CHANGE <= TOL)
        {
            converged = true;
            break;
        }

        if (P < 0)
        {
            P = TOL;
        }

        P_OLD = P;
    }

    if (!converged)
    {
        RS_PRINTF("WARNING in exact_riemann_solver::starPressure: doesn't reach tolerance within %u iterations.\n",
                  MAX_ITERATIONS);
        fprintf(Constants::getP_WARNINGS_FILE(),
                "WARNING in exact_riemann_solver::starPressure: tolerance not reached within %u iterations.\n",
                MAX_ITERATIONS);

        fprintf(stderr, "ERROR in exact_riemann_solver::starPressure: tolerance not reached within %u iterations.\n",
                MAX_ITERATIONS);
        fprintf(stderr, "change: %f\n", CHANGE);
        exit(1);
    }

    *starPressure = P;
    *starVelocity = 0.5 * (UL + UR + FR - FL);

    RS_PRINTF("number of iterations: %u\n", i + 1);
    RS_PRINTF("pressure in star region: %f\n", *starPressure);
    RS_PRINTF("velocity in star region: %f\n", *starVelocity);
}

// calculate start pressure for Newton-Raphson (linearized solution)
double exact_riemann_solver::guessPressure()
{
    double p_pv = 0.5 * (PL + PR) - 0.125 * (UR - UL) * (DL + DR) * (CL + CR);

    if (p_pv > TOL)
    {
        return p_pv;
    }
    else
    {
        return TOL;
    }
}

// sample the solution at speed S
// warning: no denominator check
void exact_riemann_solver::sample(double* D, double* U, double* P, double PM, double UM, double S)
{
    RS_PRINTF("\nsampling the solution at speed %f\n", S);

    if (S <= UM) // sampling point lies to the left of the contact discontinuity
    {
        RS_PRINTF("\tsampling point lies to the left of the contact discontinuity\n");

        if (PM <= PL) // left rarefaction
        {
            RS_PRINTF("\t\tleft wave is a rarefaction\n");

            double SHL = UL - CL; // speed of the head of the rarefaction

            if (S <= SHL) // sampled point is left data state
            {
                RS_PRINTF("\t\t\tsampled point is left data state\n");

                *D = DL;
                *U = UL;
                *P = PL;

                return;
            }
            else
            {
                double CML = CL * pow((PM / PL), G1);
                double STL = UM - CML; // speed of the tail of the rarefaction

                if (S > STL) // sampled point is star left state
                {
                    RS_PRINTF("\t\t\tsampled point is star left state\n");

                    *D = DL * pow((PM / PL), 1. / GAMMA);
                    *U = UM;
                    *P = PM;

                    return;
                }
                else // sampled point is inside left fan
                {
                    RS_PRINTF("\t\t\tsampled point is inside left fan\n");

                    *U = G5 * (CL + G7 * UL + S);
                    double C = G5 * (CL + G7 * (UL - S));
                    *D = DL * pow((C / CL), G4);
                    *P = PL * pow((C / CL), G3);

                    return;
                }
            }
        }
        else // left shock
        {
            RS_PRINTF("\t\tleft wave is a shock\n");

            double PML = PM / PL;
            double SL = UL - CL * sqrt(G2 * PML + G1); // shock speeds

            if (S <= SL) // sampled point is left data state
            {
                RS_PRINTF("\t\t\tsampled point is left data state\n");

                *D = DL;
                *U = UL;
                *P = PL;

                return;
            }
            else // sampled point is star left state
            {
                RS_PRINTF("\t\t\tsampled point is star left state\n");

                *D = DL * (PML + G6) / (PML * G6 + 1.);
                *U = UM;
                *P = PM;

                return;
            }
        }
    }
    else // sampling point lies to the right of the contact discontinuity
    {
        RS_PRINTF("\tsampling point lies to the right of the contact discontinuity\n");

        if (PM > PR) // right shock
        {
            RS_PRINTF("\t\tright wave is a shock\n");

            double PMR = PM / PR;
            double SR = UR + CR * sqrt(G2 * PMR + G1); // shock speed

            if (S >= SR) // sampled point is right data state
            {
                RS_PRINTF("\t\t\tsampled point is right data state\n");

                *D = DR;
                *U = UR;
                *P = PR;

                return;
            }
            else // sampled point is star right state
            {
                RS_PRINTF("\t\t\tsampled point is star right state\n");

                *D = DR * (PMR + G6) / (PMR * G6 + 1.);
                *U = UM;
                *P = PM;

                return;
            }
        }
        else // right rarefaction
        {
            RS_PRINTF("\t\tright wave is a rarefaction\n");

            double SHR = UR + CR; //speed of the rarefaction

            if (S >= SHR) //sample point is right data state
            {
                RS_PRINTF("\t\t\tsampled point is right data state\n");

                *D = DR;
                *U = UR;
                *P = PR;

                return;
            }
            else
            {
                double CMR = CR * pow((PM / PR), G1);
                double STR = UM + CMR; // speed of the tail of the rarefaction

                if (S <= STR) // sampled point is star right state
                {
                    RS_PRINTF("\t\t\tsampled point is star right state\n");

                    *D = DR * pow((PM / PR), 1. / GAMMA);
                    *U = UM;
                    *P = PM;

                    return;
                }
                else // sampled point is inside left fan
                {
                    RS_PRINTF("\t\t\tsampled point is inside left fan\n");

                    *U = G5 * (-CR + G7 * UR + S);
                    double C = G5 * (CR - G7 * (UR - S));
                    *D = DR * pow((C / CR), G4);
                    *P = PR * pow((C / CR), G3);

                    return;
                }
            }
        }
    }
}

// solve the Riemann problem with an exact solver and sample the solution along the speed s = x / t = 0
void exact_riemann_solver::solveAndSampleZeroSpeed(double* rho_middle, double* u_middle, double* p_middle,
                                                   double rho_left, double u_left, double p_left,
                                                   double rho_right, double u_right, double p_right)
{
    RS_PRINTF("\n---solving the Riemann problem and sampling it along zero speed---\n\n");

    calculateVariables(rho_left, u_left, p_left, rho_right, u_right, p_right);

    double starPressure = 0;
    double starVelocity = 0;

    // calculate star pressure and star velocity
    starPressureAndVelocity(&starPressure, &starVelocity);

    // sample the solution
    sample(rho_middle, u_middle, p_middle, starPressure, starVelocity, 0);

    RS_PRINTF("solution:\n");
    RS_PRINTF("density = %f\nvelocity = %f\npressure = %f\n", *rho_middle, *u_middle, *p_middle);

    resetVariables();
}