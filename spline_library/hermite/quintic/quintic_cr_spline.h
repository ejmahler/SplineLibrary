#ifndef QUINTICCATMULLROMSPLINE_H
#define QUINTICCATMULLROMSPLINE_H

#include <vector>
#include <unordered_map>

#include "quintic_hermite_spline.h"

class QuinticCRSpline final: public QuinticHermiteSpline
{
public:
    QuinticCRSpline(const std::vector<Vector3D> &points, double alpha = 0.0);
};

#endif // QUINTICCATMULLROMSPLINE_H
