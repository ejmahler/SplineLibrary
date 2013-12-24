#ifndef LOOPING_QUINTIC_HERMITE_SPLINE_H
#define LOOPING_QUINTIC_HERMITE_SPLINE_H

#include "quintic_hermite_spline.h"

class LoopingQuinticHermiteSpline : public QuinticHermiteSpline
{
public:
    virtual Vector3D getPosition(double x) const;
    virtual InterpolatedPT getTangent(double x) const;
    virtual InterpolatedPTC getCurvature(double x) const;

    virtual bool isLooping(void) const;

protected:
    LoopingQuinticHermiteSpline();
};

#endif // LOOPING_QUINTIC_HERMITE_SPLINE_H
