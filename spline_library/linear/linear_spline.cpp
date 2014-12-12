#include "linear_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cmath>
#include <cassert>

LinearSpline::LinearSpline()
{

}

LinearSpline::LinearSpline(const std::vector<Vector3D> &points, double alpha)
{
    assert(points.size() >= 2);

    this->points = points;

    int size = points.size();
    numSegments = size - 1;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValues(points, alpha, 0);
    maxT = indexToT.at(numSegments);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        InterpolationData segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at(i + 1);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D LinearSpline::getPosition(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return computePosition(t, segment);
}

Spline::InterpolatedPT LinearSpline::getTangent(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(segment)
                );
}

Spline::InterpolatedPTC LinearSpline::getCurvature(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(segment),
                Vector3D() //curvature is always 0 for linear spline
                );
}

Spline::InterpolatedPTCW LinearSpline::getWiggle(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(segment),
                Vector3D(), //curvature and wiggle are always 0 for linear spline
                Vector3D()
                );
}

double LinearSpline::getT(int index) const
{
    return indexToT.at(index);
}

double LinearSpline::getMaxT(void) const
{
    return maxT;
}

const std::vector<Vector3D> &LinearSpline::getPoints(void) const
{
    return points;
}

bool LinearSpline::isLooping(void) const
{
    return false;
}
