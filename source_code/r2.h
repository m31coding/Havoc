#ifndef R2_H
#define R2_H

/// a Point in R^2 (vector)
class Point
{

public:

    double x; ///< x coordinate
    double y; ///< y coordinate

    Point();

    Point(double d_x, double d_y);

    // standard operators
    inline Point operator+(const Point& rhs) const ///< add two vectors in R^2 (addition of each component)
    {
        Point temp;
        temp.x = x + rhs.x;
        temp.y = y + rhs.y;
        return temp;
    }

    inline Point operator-(const Point& rhs) const ///< subtract two vectors in R^2 (subtraction of each component)
    {
        Point temp;
        temp.x = x - rhs.x;
        temp.y = y - rhs.y;
        return temp;
    }

    inline void operator+=(const Point& rhs)
    {
        x += rhs.x;
        y += rhs.y;
    }

    inline void operator-=(const Point& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
    }

    inline Point operator*(const double& rhs) ///< multiply the vector with a scalar
    {
        Point temp;
        temp.x = x * rhs;
        temp.y = y * rhs;
        return temp;
    }

    inline void operator*=(const double& rhs) ///< multiply the vector with a scalar
    {
        x *= rhs;
        y *= rhs;
    }

    inline Point operator/(const double& rhs) ///< divide the vector by a scalar
    {
        Point temp;
        temp.x = x / rhs;
        temp.y = y / rhs;
        return temp;
    }

    inline void operator/=(const double& rhs) ///< divide the vector by a scalar
    {
        x /= rhs;
        y /= rhs;
    }

    inline double scalarProduct(const Point& rhs) ///< scalar product
    {
        return (x * rhs.x + y * rhs.y);
    }

    /// check whether two points are equal
    inline bool operator==(const Point& rhs) const
    {
        return (x == rhs.x && y == rhs.y);
    }

    inline double length_square() const
    {
        return (x * x + y * y);
    }

    Point orthogonal(); ///< returns the orthogonal vector via cross product: (0,0,1) x (x,y,0)

    void normalize(); ///< normalizes the vector to length 1
    Point return_normalized();

    void rotate(double phi); /// rotate the vector, angle = phi
    void rotate(double phi, Point* pivot_point); /// rotate the point, angle = phi, center of rotation: pivot_point

    /// calculate the phi coordinate of the point
    double phi() const;
};

typedef Point Vector;


namespace r2_geometry
{

// geometry, requires an explicit metric (distance function)
// here: Euclidean metric

    /**
    calculate the x-coordinate of the breakpoint of the parabolas, given two sites and the position of the sweep line
    source: analytic calculation
    */
    double breakp_x(const Point* leftSite, const Point* rightSite, double sweepY);

    /**
    calculate the y-coordinate of the breakpoint of the parabolas
    source: analytic calculation
    */
    double breakp_y(const Point* site, double xCoord, double sweepY);

    /**
    in a triangle the three perpendicular bisectors intersect in one point, calculate this intersection
    (check a triple of consecutive arcs in the sweepline algorithm)
    Source: analytic calculation
    \param p a point (corresponding to the first arc)
    \param q a point (corresponding to the second arc)
    \param r a point (corresponding to the third arc)
    \param intersection the intersection will be stored in this pointer
    \return returns true if there is an intersection (of the two bisectors). False otherwise.
    */
    bool checkTriple(const Point* p, const Point* q, const Point* r, Point* intersection);

    /**
    draw a straight line through p and q (in this order). The function checks whether the point r is in the right half-plane
    */
    bool rightTurn(const Point* p, const Point* q, const Point* r);

    /**
    check whether the breakpoints of a consecutive triple of arcs converge and if so store the intersection
    \param p the point representing the left arc
    \param q the point representing the middle arc
    \param r the point representing the right arc
    \param intersection the intersection is stored here
    \return do the breakpoints converge and does an intersection exist
    */
    bool breakpointsConverge(const Point* p, const Point* q, const Point* r, Point* intersection);

    /**
    checks whether c is in between a and b (a,b,c are in a plane, in between means inside the smaller angle)
    returns false if the angle between a and b is 180°
    */
    bool isBetween(const Point* a, const Point* b, const Point* c);

    /// calculate the distance between two points
    double distance(const Point* a, const Point* b);

    /**
    calculate the y coordinate of the lowest point of the circle in the circle event
    \param intersection a pointer to the intersection point
    \param site the pointer to a site on the circle
    */
    double lowestPointY(const Point* intersection, const Point* site);
};

#endif