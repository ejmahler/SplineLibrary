#include "quintic_hermite_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cassert>
#include <cmath>

/*
Vector3D QuinticHermiteSpline::computePosition2(double t, const InterpolationData &segment) const
{
    double a = segment.x0;
    double b = segment.x1;
    double c = segment.x2;
    double d = segment.x3;
    double e = segment.x4;
    double f = segment.x5;

    Vector3D W01 = segment.P0 * (b - t) / (b - a) + segment.P1 * (t - a) / (b - a);
    Vector3D W12 = segment.P1 * (c - t) / (c - b) + segment.P2 * (t - b) / (c - b);
    Vector3D W23 = segment.P2 * (d - t) / (d - c) + segment.P3 * (t - c) / (d - c);
    Vector3D W34 = segment.P3 * (e - t) / (e - d) + segment.P4 * (t - d) / (e - d);
    Vector3D W45 = segment.P4 * (f - t) / (f - e) + segment.P5 * (t - e) / (f - e);

    Vector3D X1 = W01 * (c - t) / (c - a) + W12 * (t - a) / (c - a);
    Vector3D X2 = W12 * (d - t) / (d - b) + W23 * (t - b) / (d - b);
    Vector3D X3 = W23 * (e - t) / (e - c) + W34 * (t - c) / (e - c);
    Vector3D X4 = W34 * (f - t) / (f - d) + W45 * (t - d) / (f - d);

    Vector3D Y12 = X1 * (c - t) / (c - b) + X2 * (t - b) / (c - b);
    Vector3D Y23 = X2 * (d - t) / (d - c) + X3 * (t - c) / (d - c);
    Vector3D Y34 = X3 * (e - t) / (e - d) + X4 * (t - d) / (e - d);

    Vector3D Z2 = Y12 * (d - t) / (d - b) + Y23 * (t - b) / (d - b);
    Vector3D Z3 = Y23 * (e - t) / (e - c) + Y34 * (t - c) / (e - c);

    return Z2 * (d - t) / (d - c) + Z3 * (t - c) / (d - c);
}*/

QuinticHermiteSpline::QuinticHermiteSpline()
{

}

QuinticHermiteSpline::QuinticHermiteSpline(
        const std::vector<Vector3D> &points,
        const std::vector<Vector3D> &tangents,
        const std::vector<Vector3D> &curvatures,
        float alpha
       )
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());
    assert(points.size() == curvatures.size());

    this->points = points;

    int size = points.size();
    int padding = 0;
    numSegments = size - 1;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValues(points, alpha, padding);
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

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangents.at(i) * tDistance;
        segment.m1 = tangents.at(i + 1) * tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.c0 = curvatures.at(i) * tDistance;
        segment.c1 = curvatures.at(i + 1) * tDistance * tDistance * tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D QuinticHermiteSpline::getPosition(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return computePosition(t, segment);
}

Spline::InterpolatedPT QuinticHermiteSpline::getTangent(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPT(
        computePosition(t, segment),
        computeTangent(t, segment)
        );
}

Spline::InterpolatedPTC QuinticHermiteSpline::getCurvature(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTC(
        computePosition(t, segment),
        computeTangent(t, segment),
        computeCurvature(t, segment)
        );
}

Spline::InterpolatedPTCW QuinticHermiteSpline::getWiggle(double x) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, x);
    auto t = (x - segment.t0) * segment.tDistanceInverse;

    return InterpolatedPTCW(
                computePosition(t, segment),
                computeTangent(t, segment),
                computeCurvature(t, segment),
                computeWiggle(t, segment)
                );
}

double QuinticHermiteSpline::getT(int index) const
{
    return indexToT.at(index);
}

double QuinticHermiteSpline::getMaxT(void) const
{
    return maxT;
}

const std::vector<Vector3D> &QuinticHermiteSpline::getPoints(void) const
{
    return points;
}


bool QuinticHermiteSpline::isLooping(void) const
{
    return false;
}
