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

    double findClosestT(const Vector3D &queryPoint) const;

private: //methods
    double brentsMethod(const Vector3D &queryPoint, double a, double aValue, double b, double bValue) const;

	double getDistanceSlope(const Vector3D &queryPoint, double t) const;

private: //data
    std::shared_ptr<Spline> spline;

    //distance in t between samples
    double sampleStep;

	//when we find a T whose abs(distance slope) is less than this tolerance, we return
    double slopeTolerance;

    //inner class used to provide an abstraction between the spline inverter and nanoflann
    class SampleTree;
    std::unique_ptr<SampleTree> sampleTree;
};

#endif // SplineInverter_H
