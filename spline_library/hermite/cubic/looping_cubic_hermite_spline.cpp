#include "looping_cubic_hermite_spline.h"

#include "spline_library/utils/spline_setup.h"

#include <cmath>
#include <cassert>

LoopingCubicHermiteSpline::LoopingCubicHermiteSpline(
        const std::vector<Vector3D> &points,
        const std::vector<Vector3D> &tangents,
        double alpha)
    :Spline(points)
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 0;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        CubicHermiteSplineKernel::InterpolationData<Vector3D> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at((i + 1)%size);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangents.at(i) * tDistance;
        segment.m1 = tangents.at((i + 1)%size) * tDistance;

        segmentData.push_back(segment);
    }
}

LoopingCubicHermiteSpline::LoopingCubicHermiteSpline(const std::vector<Vector3D> &points, double alpha)
    :Spline(points)
{
    assert(points.size() >= 4);

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 1;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //compute the tangents
    std::map<int, Vector3D> tangentMap;
    for(int i = 0; i < size + 1; i++)
    {
        double tPrev = indexToT.at(i - 1);
        double tCurrent = indexToT.at(i);
        double tNext = indexToT.at(i + 1);

        Vector3D pPrev = points.at((i - 1 + size)%size);
        Vector3D pCurrent = points.at((i + size)%size);
        Vector3D pNext = points.at((i + 1 + size)%size);

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
    for(int i = 0; i < numSegments; i++)
    {
        CubicHermiteSplineKernel::InterpolationData<Vector3D> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at((i + 1)%size);

        double tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangentMap.at(i) * tDistance;
        segment.m1 = tangentMap.at(i + 1) * tDistance;

        segmentData.push_back(segment);
    }
}

Vector3D LoopingCubicHermiteSpline::getPosition(double globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

Spline::InterpolatedPT LoopingCubicHermiteSpline::getTangent(double globalT) const
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

Spline::InterpolatedPTC LoopingCubicHermiteSpline::getCurvature(double globalT) const
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

Spline::InterpolatedPTCW LoopingCubicHermiteSpline::getWiggle(double globalT) const
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

double LoopingCubicHermiteSpline::getT(int index) const
{
    return indexToT.at(index);
}

double LoopingCubicHermiteSpline::getMaxT(void) const
{
    return maxT;
}

bool LoopingCubicHermiteSpline::isLooping(void) const
{
    return true;
}
