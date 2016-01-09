#ifndef UNIFORM_CUBIC_BSPLINE_H
#define UNIFORM_CUBIC_BSPLINE_H

#include "spline_library/spline.h"
#include "spline_library/basis/uniform_cubic_bspline_common.h"
#include "spline_library/utils/spline_setup.h"

#include <cassert>
#include <unordered_map>

template<class InterpolationType, typename floating_t=float>
class UniformCubicBSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    UniformCubicBSpline(const std::vector<InterpolationType> &points);

//methods
public:
    InterpolationType getPosition(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t x) const override;

    floating_t getT(int index) const override { return indexToT.at(index); }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return false; }

//data
protected:
    UniformCubicBSplineCommon<InterpolationType, floating_t> common;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
UniformCubicBSpline<InterpolationType,floating_t>::UniformCubicBSpline(const std::vector<InterpolationType> &points)
    :Spline<InterpolationType,floating_t>(points), common(points)
{
    assert(points.size() >= 4);
    floating_t alpha = 0.0;

    int size = points.size();
    int padding = 1;
    int numSegments = size - 3;

    //compute the T values for each point
    indexToT = SplineSetup::computeTValuesWithInnerPadding(points, alpha, padding);
    maxT = indexToT.at(padding + numSegments);
}

template<class InterpolationType, typename floating_t>
InterpolationType UniformCubicBSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    return common.getPosition(globalT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    UniformCubicBSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    return common.getTangent(globalT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    UniformCubicBSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
    return common.getCurvature(globalT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTCW
    UniformCubicBSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
    return common.getWiggle(globalT);
}

#endif // UNIFORM_CUBIC_BSPLINE_H
