#ifndef SplineInverter_H
#define SplineInverter_H

#include <vector>
#include <memory>

#include "vector3d.h"
class Spline;

class SplineInverter
{
public:
	SplineInverter(const std::shared_ptr<Spline> &spline, int samplesPerT = 10);
	~SplineInverter();

	double findClosestFast(const Vector3D &queryPoint) const;
	double findClosestPrecise(const Vector3D &queryPoint) const;

private: //methods
	double circleProjectionMethod(const Vector3D &queryPoint, double lowerBound, double upperBound) const;
    double bisectionMethod(const Vector3D &queryPoint, double lowerBound, double upperBound) const;

	double getDistanceSlope(const Vector3D &queryPoint, double t) const;

private: //data
    std::shared_ptr<Spline> spline;

	//when we find a T whose abs(distance slope) is less than this tolerance, we return
    double slopeTolerance;

    //distance in t between samples
    double sampleStep;

    //inner class used to provide an abstraction between the spline inverter and nanoflann
    class SampleTree;
    std::unique_ptr<SampleTree> sampleTree;
};

#endif // SplineInverter_H
