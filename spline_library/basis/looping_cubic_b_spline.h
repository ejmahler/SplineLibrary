#ifndef LOOPING_B_SPLINE_H
#define LOOPING_B_SPLINE_H

#include "cubic_b_spline.h"

class LoopingCubicBSpline : public CubicBSpline
{
public:
    LoopingCubicBSpline(const std::vector<Vector3D> &points);

    virtual Vector3D getPosition(double x) const override;
    virtual InterpolatedPT getTangent(double x) const override;
    virtual InterpolatedPTC getCurvature(double x) const override;
    virtual InterpolatedPTCW getWiggle(double x) const override;

    virtual bool isLooping(void) const override;
};

#endif // LOOPING_B_SPLINE_H
