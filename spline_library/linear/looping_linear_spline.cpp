#include "looping_linear_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cmath>
#include <cassert>

LoopingLinearSpline::LoopingLinearSpline()
{
}

LoopingLinearSpline::LoopingLinearSpline(const std::vector<Vector3D> &points, double alpha)
{
    assert(points.size() >= 2);

    this->points = points;

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 0;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        InterpolationData segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at((i + 1)%size);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D LoopingLinearSpline::getPosition(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return computePosition(t, segment);
}

Spline::InterpolatedPT LoopingLinearSpline::getTangent(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPT(
        computePosition(t, segment),
        computeTangent(segment)
        );
}

Spline::InterpolatedPTC LoopingLinearSpline::getCurvature(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTC(
        computePosition(t, segment),
        computeTangent(segment),
        Vector3D() //curvature is always 0 for linear spline
        );
}

Spline::InterpolatedPTCW LoopingLinearSpline::getWiggle(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(segment),
                Vector3D(), //curvature and wiggle are always 0 for linear spline
                Vector3D()
                );
}

bool LoopingLinearSpline::isLooping(void) const
{
    return true;
}
