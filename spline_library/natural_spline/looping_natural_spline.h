#ifndef LOOPINGNATURALSPLINE_H
#define LOOPINGNATURALSPLINE_H

#include "natural_spline.h"

class LoopingNaturalSpline : public NaturalSpline
{
public:
    //constructors
    public:
        LoopingNaturalSpline(const std::vector<Vector3D> &points, double alpha = 0.0);

    //methods
    public:
        virtual Vector3D getPosition(double x) const;
        virtual InterpolatedPT getTangent(double x) const;
        virtual InterpolatedPTC getCurvature(double x) const;
        virtual InterpolatedPTCW getWiggle(double x) const;

        virtual bool isLooping(void) const;
};

#endif // LOOPINGNATURALSPLINE_H
