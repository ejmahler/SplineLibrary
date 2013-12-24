#ifndef QUINTICCATMULLROMSPLINE_H
#define QUINTICCATMULLROMSPLINE_H

#include <vector>
#include <unordered_map>

#include "quintic_hermite_spline.h"

class QuinticCRSpline : public QuinticHermiteSpline
{
public:
    QuinticCRSpline(const std::vector<Vector3D> &points);
};

#endif // QUINTICCATMULLROMSPLINE_H
