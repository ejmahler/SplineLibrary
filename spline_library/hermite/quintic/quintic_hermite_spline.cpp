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

QuinticHermiteSpline::QuinticHermiteSpline(
        const std::vector<Vector3D> &points,
        const std::vector<Vector3D> &tangents,
        const std::vector<Vector3D> &curvatures,
        float alpha
       )
    :Spline(points)
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());
    assert(points.size() == curvatures.size());

    int size = points.size();
    int padding = 0;
    numSegments = size - 1;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValues(points, alpha, padding);
    maxT = indexToT.at(numSegments);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        QuinticHermiteSplineKernel::InterpolationData<Vector3D> segment;

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

QuinticHermiteSpline::QuinticHermiteSpline(const std::vector<Vector3D> &points, double alpha)
    :Spline(points)
{
    assert(points.size() >= 6);

    int size = points.size();
    int firstTangent = 1;
    int lastTangent = size - 2;
    int firstCurvature = 2;
    numSegments = size - 5;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValues(points, alpha, firstCurvature);
    maxT = indexToT.at(firstCurvature + numSegments);


    //compute the tangents
    std::map<int, Vector3D> tangentMap;
    for(int i = firstTangent; i <= lastTangent; i++)
    {
        double tPrev = indexToT.at(i - 1);
        double tCurrent = indexToT.at(i);
        double tNext = indexToT.at(i + 1);

        Vector3D pPrev = points.at(i - 1);
        Vector3D pCurrent = points.at(i);
        Vector3D pNext = points.at(i + 1);

        //the tangent is the standard catmull-rom spline tangent calculation
        tangentMap[i] =
                  pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
                + pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

             //plus a little something extra - this is derived from the pyramid contruction
             //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0,
             //yielding the standard catmull-rom formula
                - pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
    }

    //compute the curvatures
    std::map<int, Vector3D> curveMap;
    int lastCurvature = firstCurvature + numSegments;
    for(int i = firstCurvature; i <= lastCurvature; i++)
    {
        double tPrev = indexToT.at(i - 1);
        double tCurrent = indexToT.at(i);
        double tNext = indexToT.at(i + 1);

        Vector3D pPrev = tangentMap.at(i - 1);
        Vector3D pCurrent = tangentMap.at(i);
        Vector3D pNext = tangentMap.at(i + 1);

        //the tangent is the standard catmull-rom spline tangent calculation
        curveMap[i] =
                  pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
                + pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))

             //plus a little something extra - this is derived from the pyramid contruction
             //when the t values are evenly spaced (ie when alpha is 0), this whole line collapses to 0,
             //yielding the standard catmull-rom formula
                - pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
    }


    //pre-arrange the data needed for interpolation
    int lastSegment = firstCurvature + numSegments;
    for(int i = firstCurvature; i < lastSegment; i++)
    {
        QuinticHermiteSplineKernel::InterpolationData<Vector3D> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at(i + 1);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangentMap.at(i) * tDistance;
        segment.m1 = tangentMap.at(i + 1) * tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.c0 = curveMap.at(i) * tDistance * tDistance;
        segment.c1 = curveMap.at(i + 1) * tDistance * tDistance;

        segmentData.push_back(segment);
    }
}


Vector3D QuinticHermiteSpline::getPosition(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

Spline::InterpolatedPT QuinticHermiteSpline::getTangent(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPT(
        segment.computePosition(localT),
        segment.computeTangent(localT)
        );
}

Spline::InterpolatedPTC QuinticHermiteSpline::getCurvature(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTC(
        segment.computePosition(localT),
        segment.computeTangent(localT),
        segment.computeCurvature(localT)
        );
}

Spline::InterpolatedPTCW QuinticHermiteSpline::getWiggle(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTCW(
                segment.computePosition(localT),
                segment.computeTangent(localT),
                segment.computeCurvature(localT),
                segment.computeWiggle(localT)
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


bool QuinticHermiteSpline::isLooping(void) const
{
    return false;
}
