#include "looping_cubic_b_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cassert>
#include <cmath>

LoopingCubicBSpline::LoopingCubicBSpline(const std::vector<Vector3D> &points)
    :Spline(points)
{
    assert(points.size() >= 3);
    double alpha = 0.0;

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 1;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < size + 1; i++)
    {
        CubicBSplineKernel::InterpolationData<Vector3D> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.beforePoint = points.at((i - 1 + size)%size);
        segment.p0 = points.at((i + size)%size);
        segment.p1 = points.at((i + 1 + size)%size);
        segment.afterPoint = points.at((i + 2 + size)%size);

        segmentData.push_back(segment);
    }
}

Vector3D LoopingCubicBSpline::getPosition(double globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

Spline::InterpolatedPT LoopingCubicBSpline::getTangent(double globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPT(
                segment.computePosition(localT),
                segment.computeTangent(localT)
                );
}

Spline::InterpolatedPTC LoopingCubicBSpline::getCurvature(double globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTC(
                segment.computePosition(localT),
                segment.computeTangent(localT),
                segment.computeCurvature(localT)
                );
}

Spline::InterpolatedPTCW LoopingCubicBSpline::getWiggle(double globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTCW(
                segment.computePosition(localT),
                segment.computeTangent(localT),
                segment.computeCurvature(localT),
                segment.computeWiggle()
                );
}

double LoopingCubicBSpline::getT(int index) const
{
    return indexToT.at(index);
}

double LoopingCubicBSpline::getMaxT(void) const
{
    return maxT;
}

bool LoopingCubicBSpline::isLooping(void) const
{
    return true;
}
