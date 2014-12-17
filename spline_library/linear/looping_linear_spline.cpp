#include "looping_linear_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cmath>
#include <cassert>

LoopingLinearSpline::LoopingLinearSpline(const std::vector<Vector3D> &points, double alpha)
    :Spline(points)
{
    assert(points.size() >= 2);

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 0;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        LinearSplineKernel::InterpolationData<Vector3D> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at((i + 1)%size);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D LoopingLinearSpline::getPosition(double globalT) const
{
    //use modular arithmetic to bring x into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

Spline::InterpolatedPT LoopingLinearSpline::getTangent(double globalT) const
{
    //use modular arithmetic to bring x into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPT(
                segment.computePosition(localT),
                segment.computeTangent()
                );
}

Spline::InterpolatedPTC LoopingLinearSpline::getCurvature(double globalT) const
{
    //use modular arithmetic to bring x into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTC(
                segment.computePosition(localT),
                segment.computeTangent(),
                Vector3D() //curvature is always 0 for linear spline
                );
}

Spline::InterpolatedPTCW LoopingLinearSpline::getWiggle(double globalT) const
{
    //use modular arithmetic to bring x into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTCW(
                segment.computePosition(localT),
                segment.computeTangent(),
                Vector3D(), //curvature and wiggle are always 0 for linear spline
                Vector3D()
                );
}

double LoopingLinearSpline::getT(int index) const
{
    return indexToT.at(index);
}

double LoopingLinearSpline::getMaxT(void) const
{
    return maxT;
}

bool LoopingLinearSpline::isLooping(void) const
{
    return true;
}
