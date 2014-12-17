#include "cubic_b_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cmath>
#include <cassert>

CubicBSpline::CubicBSpline(const std::vector<Vector3D> &points)
    :Spline(points)
{
    assert(points.size() >= 4);
    double alpha = 0.0;

    int size = points.size();
    int padding = 1;
    numSegments = size - 3;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValues(points, alpha, padding);
    maxT = indexToT.at(padding + numSegments);

    //pre-arrange the data needed for interpolation
    for(int i = padding; i < padding + numSegments; i++)
    {
        CubicBSplineKernel::InterpolationData<Vector3D> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.beforePoint = points.at(i - 1);
        segment.p0 = points.at(i);
        segment.p1 = points.at(i + 1);
        segment.afterPoint = points.at(i + 2);

        segmentData.push_back(segment);
    }
}

Vector3D CubicBSpline::getPosition(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

Spline::InterpolatedPT CubicBSpline::getTangent(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPT(
                segment.computePosition(localT),
                segment.computeTangent(localT)
                );
}

Spline::InterpolatedPTC CubicBSpline::getCurvature(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTC(
                segment.computePosition(localT),
                segment.computeTangent(localT),
                segment.computeCurvature(localT)
                );
}

Spline::InterpolatedPTCW CubicBSpline::getWiggle(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTCW(
                segment.computePosition(localT),
                segment.computeTangent(localT),
                segment.computeCurvature(localT),
                segment.computeWiggle()
                );
}

double CubicBSpline::getT(int index) const
{
    return indexToT.at(index);
}

double CubicBSpline::getMaxT(void) const
{
    return maxT;
}

bool CubicBSpline::isLooping(void) const
{
    return false;
}
