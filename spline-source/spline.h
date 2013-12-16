#ifndef SPLINE_H
#define SPLINE_H

#include <map>
#include <vector>
#include "vector3d.h"

struct InterpolatedPV
{
	Vector3D position;
	Vector3D velocity;

	InterpolatedPV(const Vector3D &p, const Vector3D &v)
		:position(p),velocity(v)
	{}
};

struct InterpolatedPVA
{
	Vector3D position;
	Vector3D velocity;
	Vector3D acceleration;

	InterpolatedPVA(const Vector3D &p, const Vector3D &v, const Vector3D &a)
		:position(p),velocity(v),acceleration(a)
	{}
};

class Spline
{
public:
	virtual Vector3D getPosition(double x) const = 0;
	virtual InterpolatedPV getPositionVelocity(double x) const = 0;
	virtual InterpolatedPVA getPositionVelocityAcceleration(double x) const = 0;

	virtual double getT(int index) const = 0;
	virtual double getMaxT(void) const = 0;
	virtual int getNumSegments(void) const = 0;

	virtual const std::vector<Vector3D> &getPoints(void) const = 0;
	virtual bool isLoop(void) const = 0;
};
#endif // SPLINE_H
