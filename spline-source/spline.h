#ifndef SPLINE_H
#define SPLINE_H

#include <map>
#include <vector>
#include "vector3d.h"

struct InterpolatedPT
{
	Vector3D position;
    Vector3D tangent;

    InterpolatedPT(const Vector3D &p, const Vector3D &t)
        :position(p),tangent(t)
	{}
};

struct InterpolatedPTC
{
	Vector3D position;
    Vector3D tangent;
    Vector3D curvature;

    InterpolatedPTC(const Vector3D &p, const Vector3D &t, const Vector3D &c)
        :position(p),tangent(t),curvature(c)
	{}
};

class Spline
{
public:
	virtual Vector3D getPosition(double x) const = 0;
    virtual InterpolatedPT getTangent(double x) const = 0;
    virtual InterpolatedPTC getCurvature(double x) const = 0;

	virtual double getT(int index) const = 0;
	virtual double getMaxT(void) const = 0;
	virtual int getNumSegments(void) const = 0;

	virtual const std::vector<Vector3D> &getPoints(void) const = 0;
	virtual bool isLoop(void) const = 0;
};
#endif // SPLINE_H
