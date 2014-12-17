#include "cubic_hermite_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cmath>
#include <cassert>

CubicHermiteSpline::CubicHermiteSpline(const std::vector<Vector3D> &points, const std::vector<Vector3D> &tangents, double alpha)
    :Spline(points)
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());

    int size = points.size();
    int firstTangent = 0;
    numSegments = size - 1;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValues(points, alpha, firstTangent);
    maxT = indexToT.at(firstTangent + numSegments);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        CubicHermiteSplineKernel::InterpolationData<Vector3D> segment;

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

CubicHermiteSpline::CubicHermiteSpline(const std::vector<Vector3D> &points, double alpha)
    :Spline(points)
{
    assert(points.size() >= 4);

    int size = points.size();
    int firstTangent = 1;
    numSegments = size - 3;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValues(points, alpha, firstTangent);
    maxT = indexToT.at(firstTangent + numSegments);

    //compute the tangents
    std::map<int, Vector3D> tangentMap;
    for(int i = firstTangent; i < firstTangent + numSegments + 1; i++)
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

    //pre-arrange the data needed for interpolation
    for(int i = firstTangent; i < firstTangent + numSegments; i++)
    {
        CubicHermiteSplineKernel::InterpolationData<Vector3D> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at(i + 1);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangentMap.at(i) * tDistance;
        segment.m1 = tangentMap.at(i + 1) * tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D CubicHermiteSpline::getPosition(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

Spline::InterpolatedPT CubicHermiteSpline::getTangent(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPT(
                segment.computePosition(localT),
                segment.computeTangent(localT)
                );
}

Spline::InterpolatedPTC CubicHermiteSpline::getCurvature(double globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return InterpolatedPTC(
                segment.computePosition(localT),
                segment.computeTangent(localT),
                segment.computeCurvature(localT)
                );
}

Spline::InterpolatedPTCW CubicHermiteSpline::getWiggle(double globalT) const
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

double CubicHermiteSpline::getT(int index) const
{
    return indexToT.at(index);
}

double CubicHermiteSpline::getMaxT(void) const
{
    return maxT;
}

bool CubicHermiteSpline::isLooping(void) const
{
    return false;
}
