#ifndef LOOPING_B_SPLINE_H
#define LOOPING_B_SPLINE_H

#include "cubic_b_spline.h"

class LoopingCubicBSpline : public CubicBSpline
{
public:
    LoopingCubicBSpline(const std::vector<Vector3D> &points);

    virtual Vector3D getPosition(double x) const;
    virtual InterpolatedPT getTangent(double x) const;
    virtual InterpolatedPTC getCurvature(double x) const;
    virtual InterpolatedPTCW getWiggle(double x) const;

    virtual bool isLooping(void) const;
};

#endif // LOOPING_B_SPLINE_H
