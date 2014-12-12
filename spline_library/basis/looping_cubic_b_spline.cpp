#include "looping_cubic_b_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cassert>
#include <cmath>

LoopingCubicBSpline::LoopingCubicBSpline(const std::vector<Vector3D> &points)
{
    assert(points.size() >= 3);
    double alpha = 0.0;

    this->points = points;

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 1;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < size + 1; i++)
    {
        InterpolationData segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.beforePoint = points.at((i - 1 + size)%size);
        segment.p0 = points.at((i + size)%size);
        segment.p1 = points.at((i + 1 + size)%size);
        segment.afterPoint = points.at((i + 2 + size)%size);

        segmentData.push_back(segment);
    }
}

Vector3D LoopingCubicBSpline::getPosition(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    double t = x - segment.t0;

    return computePosition(t, segment);
}

Spline::InterpolatedPT LoopingCubicBSpline::getTangent(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    double t = x - segment.t0;

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(t, segment)
                );
}

Spline::InterpolatedPTC LoopingCubicBSpline::getCurvature(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    double t = x - segment.t0;

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment)
                );
}

Spline::InterpolatedPTCW LoopingCubicBSpline::getWiggle(double x) const
{
    //use modular arithmetic to bring x into an acceptable range
    x = fmod(x, numSegments);
    if(x < 0)
        x += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    double t = x - segment.t0;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment),
                computeWiggle(segment)
                );
}

bool LoopingCubicBSpline::isLooping(void) const
{
    return true;
}
