#include "cubic_b_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cmath>
#include <cassert>

CubicBSpline::CubicBSpline()
{

}

CubicBSpline::CubicBSpline(const std::vector<Vector3D> &points)
{
    assert(points.size() >= 4);
    double alpha = 0.0;

    this->points = points;

    int size = points.size();
    int padding = 1;
    numSegments = size - 3;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValues(points, alpha, padding);
    maxT = indexToT.at(padding + numSegments);

    //pre-arrange the data needed for interpolation
    for(int i = padding; i < padding + numSegments; i++)
    {
        InterpolationData segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.beforePoint = points.at(i - 1);
        segment.p0 = points.at(i);
        segment.p1 = points.at(i + 1);
        segment.afterPoint = points.at(i + 2);

        segmentData.push_back(segment);
    }
}

Vector3D CubicBSpline::getPosition(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    double t = x - segment.t0;

    return computePosition(t, segment);
}

Spline::InterpolatedPT CubicBSpline::getTangent(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    double t = x - segment.t0;

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(t, segment)
                );
}

Spline::InterpolatedPTC CubicBSpline::getCurvature(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    double t = x - segment.t0;

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment)
                );
}

Spline::InterpolatedPTCW CubicBSpline::getWiggle(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    double t = x - segment.t0;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment),
                computeWiggle(segment)
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

const std::vector<Vector3D> &CubicBSpline::getPoints(void) const
{
    return points;
}

bool CubicBSpline::isLooping(void) const
{
    return false;
}
