#ifndef LOOPING_B_SPLINE_H
#define LOOPING_B_SPLINE_H

#include <unordered_map>

#include "spline_library/spline.h"
#include "spline_library/basis/cubic_b_spline_kernel.h"

#include "spline_library/utils/spline_setup.h"

template<class InterpolationType, typename floating_t=float>
class LoopingCubicBSpline final : public Spline<InterpolationType, floating_t>
{
public:
    LoopingCubicBSpline(const std::vector<InterpolationType> &points);

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
protected:
    //a vector containing pre-computed datasets, one per segment
    //there will be lots of duplication of data here,
    //but precomputing this really speeds up the interpolation
    int numSegments;
    std::vector<CubicBSplineKernel::InterpolationData<InterpolationType>> segmentData;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
LoopingCubicBSpline<InterpolationType,floating_t>::LoopingCubicBSpline(const std::vector<InterpolationType> &points)
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() >= 3);
    floating_t alpha = 0.0;

    int size = points.size();
    numSegments = size;

    //compute the T values for each point
    int padding = 1;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //pre-arrange the data needed for interpolation
    for(int i = 0; i < size + 1; i++)
    {
        CubicBSplineKernel::InterpolationData<InterpolationType> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.beforePoint = points.at((i - 1 + size)%size);
        segment.p0 = points.at((i + size)%size);
        segment.p1 = points.at((i + 1 + size)%size);
        segment.afterPoint = points.at((i + 2 + size)%size);

        segmentData.push_back(segment);
    }
}

template<class InterpolationType, typename floating_t>
InterpolationType LoopingCubicBSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
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
    LoopingCubicBSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
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
    LoopingCubicBSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
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
    LoopingCubicBSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
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
                segment.computeWiggle()
                );
}

template<class InterpolationType, typename floating_t>
floating_t LoopingCubicBSpline<InterpolationType,floating_t>::getT(int index) const
{
    return indexToT.at(index);
}

template<class InterpolationType, typename floating_t>
floating_t LoopingCubicBSpline<InterpolationType,floating_t>::getMaxT(void) const
{
    return maxT;
}

template<class InterpolationType, typename floating_t>
bool LoopingCubicBSpline<InterpolationType,floating_t>::isLooping(void) const
{
    return true;
}

#endif // LOOPING_B_SPLINE_H
