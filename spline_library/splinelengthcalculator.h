#ifndef SPLINELENGTHCALCULATOR_H
#define SPLINELENGTHCALCULATOR_H

#include <memory>
#include <atomic>

class Spline;
class Vector3D;

class SplineLengthCalculator
{
public:
    SplineLengthCalculator(const std::shared_ptr<Spline> &spline);

    //precise but possibly slower method to find an approximation of the spline length from the begin T to the end T
    //if the spline is a looping spline and useShortestPath is true, this will try to find a shorter path around the spline by going "backwards" if possible
    //useShortestPath has no effect on non-looping splines
    double findLength(double beginT, double endT, bool useShortestPath=false) const;

private: //methods

    double computeLength(double beginT, double endT) const;

    //recursively compute the length of the given segment - if the angle between the tangents is too large, sut the segment in half and recursively compute the length of each half
    double computeLengthHelper(double beginT, const Vector3D &beginPosition, const Vector3D &beginTangentNormalized,
                                      double endT, const Vector3D &endPosition, const Vector3D &endTangentNormalized) const;

private: //data

    std::shared_ptr<Spline> spline;
    double maxT;
    mutable std::atomic<double> splineLength;
};

#endif // SPLINELENGTHCALCULATOR_H
