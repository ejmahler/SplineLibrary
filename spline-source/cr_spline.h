#ifndef CATMULLROMSPLINE_H
#define CATMULLROMSPLINE_H

#include <vector>
#include <unordered_map>

#include "cubic_hermite_spline.h"

class CRSpline final: public CubicHermiteSpline
{
public:

    CRSpline(const std::vector<Vector3D> &points, double alpha = 0.0);
};

#endif // CATMULLROMSPLINE_H
