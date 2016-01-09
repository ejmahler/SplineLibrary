#ifndef LOOPING_B_SPLINE_H
#define LOOPING_B_SPLINE_H

#include <unordered_map>

#include "spline_library/spline.h"
#include "spline_library/basis/uniform_cubic_bspline_common.h"

#include "spline_library/utils/spline_setup.h"

template<class InterpolationType, typename floating_t=float>
class LoopingUniformCubicBSpline final : public Spline<InterpolationType, floating_t>
{
public:
    LoopingUniformCubicBSpline(const std::vector<InterpolationType> &points);

//methods
public:
    InterpolationType getPosition(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t x) const override;

    floating_t getT(int index) const override { return indexToT.at(index); }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return true; }

//data
protected:
    UniformCubicBSplineCommon<InterpolationType, floating_t> common;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int, floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
LoopingUniformCubicBSpline<InterpolationType,floating_t>::LoopingUniformCubicBSpline(const std::vector<InterpolationType> &points)
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() >= 3);
    floating_t alpha = 0.0;

    int size = points.size();
    int degree = 3;

    //compute the T values for each point
    int padding = 1;
    indexToT = SplineSetup::computeLoopingTValues(points, alpha, padding);
    maxT = indexToT.at(size);

    //we need enough space to repeat the last 'degree' elements
    std::vector<InterpolationType> positions(points.size() + degree);

    //it would be easiest to just copy the points vector to the position vector, then copy the first 'degree' elements again
    //this DOES work, but interpolation begins in the wrong place (ie getPosition(0) occurs at the wrong place on the spline)
    //to fix this, we effectively "rotate" the position vector backwards, by copying point[size-1] to the beginning
    //then copying the points vector in after, then copying degree-1 elements from the beginning
    positions[0] = points[size - 1];
    std::copy(points.begin(), points.end(), positions.begin() + 1);
    std::copy_n(points.begin(), degree - 1, positions.end() - (degree - 1));

    common = UniformCubicBSplineCommon<InterpolationType, floating_t>(positions);
}

template<class InterpolationType, typename floating_t>
InterpolationType LoopingUniformCubicBSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getPosition(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    LoopingUniformCubicBSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getTangent(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    LoopingUniformCubicBSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getCurvature(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTCW
    LoopingUniformCubicBSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getWiggle(wrappedT);
}

#endif // LOOPING_B_SPLINE_H
