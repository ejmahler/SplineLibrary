#ifndef B_SPLINE_H
#define B_SPLINE_H

#include "spline_library/spline.h"
#include "spline_library/basis/cubic_b_spline_kernel.h"
#include "spline_library/utils/spline_setup.h"

#include <cassert>
#include <unordered_map>

template<class InterpolationType, typename floating_t=float>
class CubicBSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    CubicBSpline(const std::vector<InterpolationType> &points);

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
    std::vector<CubicBSplineKernel::InterpolationData<InterpolationType, floating_t>> segmentData;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
CubicBSpline<InterpolationType,floating_t>::CubicBSpline(const std::vector<InterpolationType> &points)
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() >= 4);
    floating_t alpha = 0.0;

    int size = points.size();
    int padding = 1;
    int numSegments = size - 3;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValuesWithInnerPadding(points, alpha, padding);
    maxT = indexToT.at(padding + numSegments);

    //pre-arrange the data needed for interpolation
    for(int i = padding; i < padding + numSegments; i++)
    {
        CubicBSplineKernel::InterpolationData<InterpolationType, floating_t> segment;

        segment.t0 = indexToT.at(i);
        segment.t1 = indexToT.at(i + 1);

        segment.beforePoint = points.at(i - 1);
        segment.p0 = points.at(i);
        segment.p1 = points.at(i + 1);
        segment.afterPoint = points.at(i + 2);

        segmentData.push_back(segment);
    }
}

template<class InterpolationType, typename floating_t>
InterpolationType CubicBSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return segment.computePosition(localT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    CubicBSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    auto segment = SplineSetup::getSegmentForT(segmentData, globalT);
    auto localT = segment.computeLocalT(globalT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                segment.computePosition(localT),
                segment.computeTangent(localT)
                );
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    CubicBSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
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
    CubicBSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
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
floating_t CubicBSpline<InterpolationType,floating_t>::getT(int index) const
{
    return indexToT.at(index);
}

template<class InterpolationType, typename floating_t>
floating_t CubicBSpline<InterpolationType,floating_t>::getMaxT(void) const
{
    return maxT;
}

template<class InterpolationType, typename floating_t>
bool CubicBSpline<InterpolationType,floating_t>::isLooping(void) const
{
    return false;
}

#endif // B_SPLINE_H
