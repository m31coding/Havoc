#include "config.h"
#include "r2.h"
#include <cmath>
#include <cstdio>
#include "constants.h"
#include "assert.h"
#include <cstdlib>
#include "numerics.h"

Point::Point()
        :
        x(0),
        y(0)
{

}

Point::Point(double d_x, double d_y)
        : x(d_x),
          y(d_y)
{

}

// sets the point to its orthogonal point via cross product
Point Point::orthogonal()
{
    Point temp;
    temp.x = -y;
    temp.y = x;

    return temp;
}

// normalizes the vector to length 1
void Point::normalize()
{
    double length = sqrt(x * x + y * y);
    if (length == 0)
    {
        fprintf(stderr, "ERROR in Point::normalize: vector has zero length\n");
        exit(1);
        return;
    }

    x /= length;
    y /= length;

}

Vector Vector::return_normalized()
{
    double length = sqrt(x * x + y * y);
    if (length == 0)
    {
        return Vector(0, 0);
    }

    return Vector(x / length, y / length);
}

// rotate the vector, angle = phi
void Point::rotate(double phi)
{
    double x_rotated = cos(phi) * x - sin(phi) * y;
    y = sin(phi) * x + cos(phi) * y;

    x = x_rotated;
}

// rotate the point, angle = phi, center of rotation: pivot_point
void Point::rotate(double phi, Point* pivot_point)
{
    x -= pivot_point->x;
    y -= pivot_point->y;

    this->rotate(phi);

    x += pivot_point->x;
    y += pivot_point->y;
}

// calculate the phi coordinate of the point
double Point::phi() const
{
    // if both arguments are zero an error occurs in atan2
    if ((x == 0) && (y == 0))
    {
        fprintf(stderr, "ERROR in Point::phi: both arguments are zero\n");
        return 0;
    }

    double p = atan2(y, x); // in [-pi,pi]

    if (p < 0)
    { return p + 2 * M_PI; }
    return p;
}

// calculate the x-coordinate of the breakpoint of the parabolas, given two sites and the position of the sweep line
// source: analytic calculation
double r2_geometry::breakp_x(const Point* left_site, const Point* right_site, double sweep_y)
{
    // transformation
    //double x1=0; // left_site->x-left_site->x, var x1 not needed
    double y1 = left_site->y - sweep_y;

    double x2 = right_site->x - left_site->x;
    double y2 = right_site->y - sweep_y;

    double x = 0;
    double y12 = y1 * y2;

    // calculation of x
    if (y12 != 0 && y1 != y2) // both points are not on the sweepline and have not the same y coordinate.
    {
        // stable calculation of the root
        double solution_plus = 0;
        double solution_minus = 0;

        numerics::solve_quadratic_equation(y2 - y1, 2 * x2 * y1, (y1 - y2) * y12 - x2 * x2 * y1, &solution_plus,
                                           &solution_minus);
        x = solution_plus;

        // skip newton iterations
        goto no_instability_fix;

        if (fabs(y1 - y2) < 0.000001) // numerical instability
        {
            // newton iterations
            double change = 0;
            double x_old = 0;
            bool converged = false;

            // at most four iterations because there is the danger that the wrong root is found
            for (int i = 0; i < 4; i++)
            {
                x_old = x;
                x = x -
                    0.5 * (y2 * (y1 * y1 + x * x) - y1 * (y2 * y2 + (x - x2) * (x - x2))) / (x * (y2 - y1) + y1 * x2);
                change = 2 * fabs((x - x_old) / (x + x_old));
                if (change <= 0.000001)
                {
                    converged = true;
                    break;
                }
            }

            if (!converged)
            {
                fprintf(stderr, "WARNING in r2_geometry::breakp_x: tolerance not reached within maximum iterations.\n");
                fprintf(stderr, "sites: x1=%f y1=%f, x2=%f, y2=%f\n", left_site->x, left_site->y, right_site->x,
                        right_site->y);
            }
        }

        no_instability_fix:;

        // retransform
        x += left_site->x;

        return x;
    }

    if (y1 == y2 && x2 != 0 && y2 != 0) // both points are not on the sweepline  but have the same y coordinate
    {
        x = x2 * 0.5;

        // retransform
        x += left_site->x;

        return x;
    }

    if (y2 == 0 && y1 != 0) // point right site on sweep line
    {
        x = x2;

        // retransform
        x += left_site->x;

        return x;
    }

    if (y1 == 0 && y2 != 0) // point left site on sweep line
    {
        x = 0;

        // retransform
        x += left_site->x;

        return x;
    }

#ifdef SPD_FIX_IN_PROCESS
    if(0 == x2 && y1 == y2) // point left_site == point right_site
    {
        return left_site->x;
    }
#endif

    // point left_site == point right_site, or y1=y2=0.

    fprintf(stderr, "ERROR in r2_geometry::breakp_x: can't calculate x-coordinate of the breakpoint\n");
    fprintf(stderr, "x1=%.8f y1=%.8f, x2=%.8f, y2=%.8f\n", 0.0, y1, x2, y2);
    fprintf(stderr, "sites: x1=%.8f y1=%.8f, x2=%.8f, y2=%.8f\n", left_site->x, left_site->y, right_site->x,
            right_site->y);

    if (left_site == right_site)
    {
        fprintf(stderr, "same sites!\n");
    }
    else
    {
        fprintf(stderr, "NOT the same sites!\n");
    }

    exit(1);
    return 0;
}

// calculate the y-coordinate of the breakpoint of the parabola
double r2_geometry::breakp_y(const Point* site, double xCoord, double sweepY)
{
    double denominator = 2 * (site->y - sweepY);

    if (denominator == 0)
    {
        fprintf(stderr, "ERROR in breakp_y: can't calculate y-coordinate of the breakpoint\n");
        return 0;
    }

    return ((pow(site->y, 2) - pow(sweepY, 2) + pow(site->x - xCoord, 2)) / (2 * (site->y - sweepY)));
}

#ifdef SPD_FIX_IN_PROCESS
// r1, p == q
static bool breakpointsConvergeHelper1(const Point* q, const Point* r, Point* intersection)
{
    double q1 = q->x;
    double q2 = q->y;

    double r1 = r->x;
    double r2 = r->y;

    if(r2 >= q2) return false;

    intersection->x = q1;

    // analytical calculation
    assert(q2 - r2 != 0);

    intersection->y = 0.5 * (q2*q2 - r2*r2 - pow(r1-q1, 2)) / (q2-r2);

    return true;
}

// helper2, q == r
static bool breakpointsConvergeHelper2(const Point* p, const Point* q, Point* intersection)
{
    double p1 = p->x;
    double p2 = p->y;

    double q1 = q->x;
    double q2 = q->y;

    if(p2 >= q2) return false;

    intersection->x = q1;

    // analytical calculation
    assert(p2 - q2 != 0);

    intersection->y = 0.5 * (p2*p2 - q2*q2 + pow(q1-p1, 2)) / (p2-q2);

    return true;
}

// helper3, p == r
static bool breakpointsConvergeHelper3(const Point* p, const Point* q, Point* intersection)
{
    double p1 = p->x;
    double p2 = p->y;

    double q1 = q->x;
    double q2 = q->y;

    if(p2 >= q2) return false;

    intersection->x = p1;

    // analytical calculation
    assert(p2 - q2 != 0);

    intersection->y = 0.5 * (p2*p2 - q2*q2 - pow(q1-p1, 2)) / (p2-q2);

    return true;
}
#endif

// in a triangle the three perpendicular bisectors intersect in one point, calculate this intersection
// (check a triple of consecutive arcs in the sweep line algorithm)
// not needed, contained in breakpointsConverge
bool r2_geometry::checkTriple(const Point* p, const Point* q, const Point* r, Point* intersection)
{
#ifdef SPD_FIX_IN_PROCESS
    // handle same point degeneracy
    if(*p == *q) return (breakpointsConvergeHelper1(q,r,intersection));
    if(*q == *r) return (breakpointsConvergeHelper2(p,q,intersection));
    if(*p == *r) return (breakpointsConvergeHelper3(p,q,intersection));
#endif

    double p1 = p->x;
    double p2 = p->y;

    double q1 = q->x;
    double q2 = q->y;

    double r1 = r->x;
    double r2 = r->y;

    // denominator of y
    double denominator = 2 * ((p2 - q2) * (q1 - r1) - (q2 - r2) * (p1 - q1));

    if (denominator == 0)
    {
        return false;
    }

    double temp1 = (r1 * r1 + r2 * r2 - q1 * q1 - q2 * q2);
    double temp2 = (q1 * q1 + q2 * q2 - p1 * p1 - p2 * p2);

    // numerator of y
    double nominator = temp1 * (p1 - q1) - temp2 * (q1 - r1);

    intersection->y = nominator / denominator;

    // numerator of x
    nominator = temp1 * (p2 - q2) - temp2 * (q2 - r2);

    intersection->x = nominator / denominator * (-1);

    return true;
}


// draw a straight line through p and q (in this order)
// the function checks whether the point r is in the right half-plane
// Rolf Klein: Algorithmische Geometrie, page 27
// not needed, contained in breakpointsConverge
bool r2_geometry::rightTurn(const Point* p, const Point* q, const Point* r)
{
    double result = (q->y - p->y) * (r->x - p->x) + (p->x - q->x) * (r->y - p->y);

    return (result > 0);
}


// check whether the breakpoints of a consecutive triple of arcs converge and if so store the intersection
bool r2_geometry::breakpointsConverge(const Point* p, const Point* q, const Point* r, Point* intersection)
{
#ifdef SPD_FIX_IN_PROCESS
    // handle same point degeneracy
    if(*p == *q) return (breakpointsConvergeHelper1(q,r,intersection));
    if(*q == *r) return (breakpointsConvergeHelper2(p,q,intersection));
    if(*p == *r) return (breakpointsConvergeHelper3(p,q,intersection));
#endif

    double p1 = p->x;
    double p2 = p->y;

    double q1 = q->x;
    double q2 = q->y;

    double r1 = r->x;
    double r2 = r->y;

    // denominator of y
    double denominator = 2 * ((p2 - q2) * (q1 - r1) - (q2 - r2) * (p1 - q1));

    if (denominator <= 0) // no right turn
    {
        return false;
    }

    double temp1 = (r1 * r1 + r2 * r2 - q1 * q1 - q2 * q2);
    double temp2 = (q1 * q1 + q2 * q2 - p1 * p1 - p2 * p2);

    // numerator of y
    double numerator = temp1 * (p1 - q1) - temp2 * (q1 - r1);

    intersection->y = numerator / denominator;

    // numerator of x
    numerator = temp1 * (p2 - q2) - temp2 * (q2 - r2);

    intersection->x = numerator / denominator * (-1);

    return true;
}

// checks whether c is in between a and b (a,b,c are in a plane, in between means inside the smaller angle)
bool r2_geometry::isBetween(const Point* a, const Point* b, const Point* c)
{
    double denominator = (a->x) * (b->y) - (b->x) * (a->y);

    if (denominator == 0)
    { return true; } // the angle between a and b is 180°

    double t2 = ((a->x) * (c->y) - (a->y) * (c->x)) / denominator;
    double t1 = ((c->y) * (b->x) - (b->y) * (c->x)) / denominator * (-1);

    return (t1 >= 0 && t2 >= 0);
}

// calculate the distance between two points
double r2_geometry::distance(const Point* a, const Point* b)
{
    return sqrt(pow((a->x) - (b->x), 2) + pow((a->y) - (b->y), 2));
}

// calculate the y coordinate of the lowest point of the circle in the circle event
double r2_geometry::lowestPointY(const Point* intersection, const Point* site)
{
    return (intersection->y - distance(intersection, site));
}