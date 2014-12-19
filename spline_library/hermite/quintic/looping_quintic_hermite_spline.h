#ifndef LOOPING_QUINTIC_HERMITE_SPLINE_H
#define LOOPING_QUINTIC_HERMITE_SPLINE_H

#include <unordered_map>

#include "spline_library/spline.h"
#include "spline_library/hermite/quintic/quintic_hermite_spline_kernel.h"

#include "spline_library/utils/spline_setup.h"

template<class InterpolationType, typename floating_t=float>
class LoopingQuinticHermiteSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    LoopingQuinticHermiteSpline(const std::vector<InterpolationType> &points,
                                const std::vector<InterpolationType> &tangents,
                                const std::vector<InterpolationType> &curvatures,
                                floating_t alpha = 0.0
                                );

    LoopingQuinticHermiteSpline(const std::vector<InterpolationType> &points, floating_t alpha = 0.0);

//methods
public:
    InterpolationType getPosition(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t x) const override;

    floating_t getT(int index) const override;
    floating_t getMaxT(void) const override;

    bool isLooping(void) const override;

//data
private:
    //a vector containing pre-computed datasets, one per segment
    //there will be lots of duplication of data here,
    //but precomputing this really speeds up the interpolation
    int numSegments;
    std::vector<QuinticHermiteSplineKernel::InterpolationData<InterpolationType, floating_t>> segmentData;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
LoopingQuinticHermiteSpline<InterpolationType,floating_t>::LoopingQuinticHermiteSpline(
        const std::vector<InterpolationType> &points,
        const std::vector<InterpolationType> &tangents,
        const std::vector<InterpolationType> &curvatures,
        floating_t alpha
        )
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() >= 2);
    assert(points.size() == tangents.size());
    assert(points.size() == curvatures.size());

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 0;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < numSegments; i++)
    {
        QuinticHermiteSplineKernel::InterpolationData<InterpolationType, floating_t> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at((i + 1)%size);

        floating_t tDistance = segment.t1 - segment.t0;
        segment.tDistanceInverse = 1 / tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.m0 = tangents.at(i) * tDistance;
        segment.m1 = tangents.at((i + 1)%size) * tDistance;

        //we scale the tangents by this segment's t distance, because wikipedia says so
        segment.c0 = curvatures.at(i) * tDistance * tDistance * tDistance;
        segment.c1 = curvatures.at((i + 1)%size) * tDistance * tDistance * tDistance;

        segmentData.push_back(segment);
    }
}

template<class InterpolationType, typename floating_t>
LoopingQuinticHermiteSpline<InterpolationType,floating_t>::LoopingQuinticHermiteSpline(const std::vector<InterpolationType> &points, floating_t alpha)
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() >= 3);

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 2;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //compute the tangents
    std::unordered_map<int, InterpolationType> tangentMap;
    for(int i = -1; i < size + 2; i++)
    {
        floating_t tPrev = indexToT.at(i - 1);
        floating_t tCurrent = indexToT.at(i);
        floating_t tNext = indexToT.at(i + 1);

        InterpolationType pPrev = points.at((i - 1 + size)%size);
        InterpolationType pCurrent = points.at((i + size)%size);
        InterpolationType pNext = points.at((i + 1 + size)%size);

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
    std::unordered_map<int, InterpolationType> curveMap;
    for(int i = 0; i < size + 1; i++)
    {
        floating_t tPrev = indexToT.at(i - 1);
        floating_t tCurrent = indexToT.at(i);
        floating_t tNext = indexToT.at(i + 1);

        InterpolationType pPrev = tangentMap.at(i - 1);
        InterpolationType pCurrent = tangentMap.at(i);
        InterpolationType pNext = tangentMap.at(i + 1);

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
    for(int i = 0; i < numSegments; i++)
    {
        QuinticHermiteSplineKernel::InterpolationData<InterpolationType, floating_t> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.p0 = points.at(i);
        segment.p1 = points.at((i + 1)%size);

        floating_t tDistance = segment.t1 - segment.t0;
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

template<class InterpolationType, typename floating_t>
InterpolationType LoopingQuinticHermiteSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    LoopingQuinticHermiteSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
        segment.computePosition(localT),
        segment.computeTangent(localT)
        );
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    LoopingQuinticHermiteSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
        segment.computePosition(localT),
        segment.computeTangent(localT),
        segment.computeCurvature(localT)
        );
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTCW
    LoopingQuinticHermiteSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
    //use modular arithmetic to bring globalT into an acceptable range
    globalT = fmod(globalT, numSegments);
    if(globalT < 0)
        globalT += numSegments;

    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                segment.computePosition(localT),
                segment.computeTangent(localT),
                segment.computeCurvature(localT),
                segment.computeWiggle(localT)
                );
}

template<class InterpolationType, typename floating_t>
floating_t LoopingQuinticHermiteSpline<InterpolationType,floating_t>::getT(int index) const
{
    return indexToT.at(index);
}

template<class InterpolationType, typename floating_t>
floating_t LoopingQuinticHermiteSpline<InterpolationType,floating_t>::getMaxT(void) const
{
    return maxT;
}

template<class InterpolationType, typename floating_t>
bool LoopingQuinticHermiteSpline<InterpolationType,floating_t>::isLooping(void) const
{
    return true;
}

#endif // LOOPING_QUINTIC_HERMITE_SPLINE_H
