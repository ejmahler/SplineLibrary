#include "cubic_hermite_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cmath>
#include <cassert>

CubicHermiteSpline::CubicHermiteSpline()
{

}

CubicHermiteSpline::CubicHermiteSpline(const std::vector<Vector3D> &points, const std::vector<Vector3D> &tangents, double alpha)
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());

    this->points = points;

    int size = points.size();
    int firstTangent = 0;
    numSegments = size - 1;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValues(points, alpha, firstTangent);
    maxT = indexToT.at(firstTangent + numSegments);

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

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangents.at(i) * tDistance;
        segment.m1 = tangents.at(i + 1) * tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D CubicHermiteSpline::getPosition(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return computePosition(t, segment);
}

Spline::InterpolatedPT CubicHermiteSpline::getTangent(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPT(
                computePosition(t, segment),
                computeTangent(t, segment)
                );
}

Spline::InterpolatedPTC CubicHermiteSpline::getCurvature(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTC(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment)
                );
}

Spline::InterpolatedPTCW CubicHermiteSpline::getWiggle(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment),
                computeWiggle(segment)
                );
}

double CubicHermiteSpline::getT(int index) const
{
    return indexToT.at(index);
}

double CubicHermiteSpline::getMaxT(void) const
{
    return maxT;
}

const std::vector<Vector3D> &CubicHermiteSpline::getPoints(void) const
{
    return points;
}

bool CubicHermiteSpline::isLooping(void) const
{
    return false;
}
