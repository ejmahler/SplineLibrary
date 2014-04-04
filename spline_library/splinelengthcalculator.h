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
    double findLength(double beginT, double endT, bool useShortestPath=false, double eps = 0.0025) const;

private: //methods
    double computeLength(double beginT, double endT, double eps) const;

    //recursively compute the length of the given segment
    double computeLengthHelper(double beginT, const Vector3D &beginPosition,
                                double endT, const Vector3D &endPosition, double eps) const;

private: //data

    std::shared_ptr<Spline> spline;
    double maxT;
    mutable std::atomic<double> atomic_splineLength;
};

#endif // SPLINELENGTHCALCULATOR_H
