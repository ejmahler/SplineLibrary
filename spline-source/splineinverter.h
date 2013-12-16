#ifndef SplineInverter_H
#define SplineInverter_H

#include <vector>
#include <memory>

#include "vector3d.h"
class Spline;

class SplineInverter
{
public:
	SplineInverter(const std::shared_ptr<Spline> &spline, int samplesPerT);
	~SplineInverter();

	double findClosestFast(const Vector3D &queryPoint) const;
	double findClosestPrecise(const Vector3D &queryPoint) const;

private:
	struct DistanceSample;
	struct SplineSample;

	double findClosestSample(const Vector3D &queryPoint) const;
	

	double circleProjectionMethod(const Vector3D &queryPoint, double lowerBound, double upperBound) const;
	double bisectionMethod(const Vector3D &queryPoint, double lowerBound, double upperBound) const;

	int findSampleIndex(double xValue) const;
	double getDistanceSlope(const Vector3D &queryPoint, double t) const;

	std::shared_ptr<Spline> spline;

	//how far apart samples are taken (in units of T), should be 10+ samples per point
	double sampleStep;

	//when we find a T whose abs(distance slope) is less than this tolerance, we return
	double slopeTolerance;

	//a list of all spline samples, sorted by x coordinate
	std::vector<SplineSample> splineSamples;
};

struct SplineInverter::DistanceSample
{
	//we could theoretically store the actual distance here as well as the derivatives,
	//but nothing needs it so we're just leaving it out
	double velocity;
	double acceleration;

	DistanceSample(double v, double a)
		:velocity(v),acceleration(a)
	{}
};

struct SplineInverter::SplineSample
{
	double t;
	Vector3D position;

	bool operator<(const SplineSample& other) const
	{
		return position.x() < other.position.x();
	}
};

#endif // SplineInverter_H
