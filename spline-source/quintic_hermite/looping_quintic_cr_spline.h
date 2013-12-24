#ifndef LOOPINGQUINTICCRSPLINE_H
#define LOOPINGQUINTICCRSPLINE_H

#include <vector>
#include <unordered_map>

#include "looping_quintic_hermite_spline.h"

class LoopingQuinticCRSpline : public LoopingQuinticHermiteSpline
{
public:
    LoopingQuinticCRSpline(const std::vector<Vector3D> &points);
};

#endif // LOOPINGQUINTICCRSPLINE_H
