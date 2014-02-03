#ifndef SPLINE_H
#define SPLINE_H

#include <map>
#include <vector>
#include "vector3d.h"

class Spline
{
public:
    struct InterpolatedPT;
    struct InterpolatedPTC;
    struct InterpolatedPTCW;

	virtual Vector3D getPosition(double x) const = 0;
    virtual InterpolatedPT getTangent(double x) const = 0;
    virtual InterpolatedPTC getCurvature(double x) const = 0;
    virtual InterpolatedPTCW getWiggle(double x) const = 0;

	virtual double getT(int index) const = 0;
    virtual double getMaxT(void) const = 0;

	virtual const std::vector<Vector3D> &getPoints(void) const = 0;
    virtual bool isLooping(void) const = 0;
};

struct Spline::InterpolatedPT
{
    Vector3D position;
    Vector3D tangent;

    InterpolatedPT(const Vector3D &p, const Vector3D &t)
        :position(p),tangent(t)
    {}
};

struct Spline::InterpolatedPTC
{
    Vector3D position;
    Vector3D tangent;
    Vector3D curvature;

    InterpolatedPTC(const Vector3D &p, const Vector3D &t, const Vector3D &c)
        :position(p),tangent(t),curvature(c)
    {}
};

struct Spline::InterpolatedPTCW
{
    Vector3D position;
    Vector3D tangent;
    Vector3D curvature;
    Vector3D wiggle;

    InterpolatedPTCW(const Vector3D &p, const Vector3D &t, const Vector3D &c, const Vector3D &w)
        :position(p),tangent(t),curvature(c), wiggle(w)
    {}
};

#endif // SPLINE_H
