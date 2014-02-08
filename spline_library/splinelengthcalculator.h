#ifndef SPLINELENGTHCALCULATOR_H
#define SPLINELENGTHCALCULATOR_H

#include <memory>

class Spline;
class Vector3D;

class SplineLengthCalculator
{
public:
    SplineLengthCalculator(const std::shared_ptr<Spline> &spline);

    //fast but possibly imprecise method to find an approximation of the spline length from the begin T to the end T
    //if the spline is a looping spline and useShortestPath is true, this will try to find a shorter path around the spline by going "backwards" if possible
    //useShortestPath has no effect on non-looping splines
    double findLengthFast(double beginT, double endT, bool useShortestPath=false) const;

    //precise but possibly slower method to find an approximation of the spline length from the begin T to the end T
    //if the spline is a looping spline and useShortestPath is true, this will try to find a shorter path around the spline by going "backwards" if possible
    //useShortestPath has no effect on non-looping splines
    double findLengthPrecise(double beginT, double endT, bool useShortestPath=false) const;

private: //methods

    double computeLengthFast(double beginT, double endT) const;
    double computeLengthPrecise(double beginT, double endT) const;

    //compute the length of a circle arc passing through the specified points, and roughly passing though the specified tangents
    double computeArcLength(const Vector3D &beginPosition, const Vector3D &beginTangent,
                            const Vector3D &endPosition, const Vector3D &endTangent) const;

    //recursively compute the length of the given segment - if the angle between the tangents is too large, sut the segment in half and recursively compute the length of each half
    double computeLengthPreciseHelper(double beginT, const Vector3D &beginPosition, const Vector3D &beginTangent,
                                      double endT, const Vector3D &endPosition, const Vector3D &endTangent) const;

private: //data

    std::shared_ptr<Spline> spline;
    double maxT;
    double splineLength;
};

#endif // SPLINELENGTHCALCULATOR_H
