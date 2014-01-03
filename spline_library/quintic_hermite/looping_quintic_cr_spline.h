#ifndef LOOPINGQUINTICCRSPLINE_H
#define LOOPINGQUINTICCRSPLINE_H

#include <vector>
#include <unordered_map>

#include "looping_quintic_hermite_spline.h"

class LoopingQuinticCRSpline final: public LoopingQuinticHermiteSpline
{
public:
    LoopingQuinticCRSpline(const std::vector<Vector3D> &points, double alpha = 0.0);
};

#endif // LOOPINGQUINTICCRSPLINE_H
