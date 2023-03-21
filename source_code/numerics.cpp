#include "numerics.h"
#include <cmath>
#include <cstdio>

// http://en.wikipedia.org/wiki/Loss_of_significance
// a stable way of solving a quadratic equation.
void numerics::solve_quadratic_equation(double a, double b, double c, double* solution_plus, double* solution_minus)
{
    // solution plus/minus corresponds to the abc-formula solution

    double p = -b / a;
    double q = c / a;

    if (p > 0)
    {
        if (a > 0)
        {
            *solution_plus = p * 0.5 + sqrt(p * p * 0.25 - q);
            *solution_minus = q / (*solution_plus);
        }
        else
        {
            *solution_minus = p * 0.5 + sqrt(p * p * 0.25 - q);
            *solution_plus = q / (*solution_minus);
        }
    }
    else // p<=0
    {
        if (a > 0)
        {
            *solution_minus = p * 0.5 - sqrt(p * p * 0.25 - q);
            *solution_plus = q / (*solution_minus);
        }
        else
        {
            *solution_plus = p * 0.5 - sqrt(p * p * 0.25 - q);
            *solution_minus = q / (*solution_plus);
        }
    }
}