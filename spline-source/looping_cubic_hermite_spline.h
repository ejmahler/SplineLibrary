#ifndef LOOPING_CUBIC_HERMITE_SPLINE_H
#define LOOPING_CUBIC_HERMITE_SPLINE_H

#include "cubic_hermite_spline.h"

class LoopingCubicHermiteSpline : public CubicHermiteSpline
{
public:
    virtual Vector3D getPosition(double x) const;
    virtual InterpolatedPT getTangent(double x) const;
    virtual InterpolatedPTC getCurvature(double x) const;

    virtual bool isLooping(void) const;

protected:
    LoopingCubicHermiteSpline();
};

#endif // LOOPING_CUBIC_HERMITE_SPLINE_H
