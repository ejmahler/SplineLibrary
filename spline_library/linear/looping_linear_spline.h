#ifndef LOOPING_LINEAR_SPLINE_H
#define LOOPING_LINEAR_SPLINE_H

#include "linear_spline.h"

class LoopingLinearSpline : public LinearSpline
{
//constructors
public:
    LoopingLinearSpline(const std::vector<Vector3D> &points, double alpha = 0.0);
protected:
    //you're only allowed to create one of these without point data if a subclass is providing the point data
    LoopingLinearSpline();

//methods
protected:
    virtual Vector3D getPosition(double x) const override;
    virtual InterpolatedPT getTangent(double x) const override;
    virtual InterpolatedPTC getCurvature(double x) const override;
    virtual InterpolatedPTCW getWiggle(double x) const override;

    virtual bool isLooping(void) const override;
};

#endif // LOOPING_LINEAR_SPLINE_H
