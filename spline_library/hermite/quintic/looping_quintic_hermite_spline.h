#ifndef LOOPING_QUINTIC_HERMITE_SPLINE_H
#define LOOPING_QUINTIC_HERMITE_SPLINE_H

#include "quintic_hermite_spline.h"

class LoopingQuinticHermiteSpline : public QuinticHermiteSpline
{
//constructors
public:
    LoopingQuinticHermiteSpline(const std::vector<Vector3D> &points,
                                const std::vector<Vector3D> &tangents,
                                const std::vector<Vector3D> &curvatures,
                                float alpha
                                );
protected:
    //you're only allowed to create one of these without point data if a subclass is providing the point data
    LoopingQuinticHermiteSpline();


//methods
protected:
    virtual Vector3D getPosition(double x) const;
    virtual InterpolatedPT getTangent(double x) const;
    virtual InterpolatedPTC getCurvature(double x) const;
    virtual InterpolatedPTCW getWiggle(double x) const;

    virtual bool isLooping(void) const;
};

#endif // LOOPING_QUINTIC_HERMITE_SPLINE_H
