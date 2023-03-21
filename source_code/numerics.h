#ifndef NUMERICS_H
#define NUMERICS_H

namespace numerics
{
    /// solve the quadratic equation ax^2+bx+c=0
    void solve_quadratic_equation(double a, double b, double c, double* solution_plus, double* solution_minus);
};
#endif