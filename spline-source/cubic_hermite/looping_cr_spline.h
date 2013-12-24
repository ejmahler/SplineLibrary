#ifndef LOOPING_CR_SPLINE_H
#define LOOPING_CR_SPLINE_H

#include <vector>
#include <unordered_map>

#include "looping_cubic_hermite_spline.h"

class LoopingCRSpline final: public LoopingCubicHermiteSpline
{
public:
    LoopingCRSpline(const std::vector<Vector3D> &points, double alpha = 0.0);
};

#endif // LOOPING_CR_SPLINE_H
