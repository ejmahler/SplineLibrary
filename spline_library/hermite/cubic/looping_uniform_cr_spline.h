#pragma once

#include "../../spline.h"
#include "uniform_cr_spline_common.h"

#include "../../utils/spline_common.h"

template<class InterpolationType, typename floating_t=float>
class LoopingUniformCRSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    LoopingUniformCRSpline(const std::vector<InterpolationType> &points);

//methods
public:
    InterpolationType getPosition(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t x) const override;

    floating_t arcLength(floating_t a, floating_t b) const override;
    floating_t totalLength(void) const override { return ArcLength::arcLength::totalLength(*this); }

    floating_t getT(int index) const override { return index; }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return true; }

    size_t segmentCount(void) const override { return common.segmentCount(); }
    size_t segmentForT(floating_t t) const override { return common.segmentForT(t); }
    floating_t segmentT(size_t segmentIndex) const override { return segmentIndex; }
    floating_t segmentArcLength(size_t segmentIndex, floating_t a, floating_t b) const override { return common.segmentLength(segmentIndex, a, b); }

//data
private:
    UniformCRSplineCommon<InterpolationType, floating_t> common;

    floating_t maxT;
};

template<class InterpolationType, typename floating_t>
LoopingUniformCRSpline<InterpolationType,floating_t>::LoopingUniformCRSpline(const std::vector<InterpolationType> &points)
    :Spline<InterpolationType,floating_t>(points), maxT(points.size())
{
    assert(points.size() >= 4);

    //we need enough space to repeat the last 'degree' elements
    std::vector<InterpolationType> positions(points.size() + 3);

    //it would be easiest to just copy the points vector to the position vector, then copy the first 'degree' elements again
    //this DOES work, but interpolation begins in the wrong place (ie getPosition(0) occurs at the wrong place on the spline)
    //to fix this, we effectively "rotate" the position vector backwards, by copying point[size-1] to the beginning
    //then copying the points vector in after, then copying degree-1 elements from the beginning
    positions[0] = points.back();
    std::copy(points.begin(), points.end(), positions.begin() + 1);
    std::copy_n(points.begin(), 2, positions.end() - 2);

    common = UniformCRSplineCommon<InterpolationType, floating_t>(std::move(positions));
}



template<class InterpolationType, typename floating_t>
InterpolationType LoopingUniformCRSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    floating_t wrappedT = SplineCommon::wrapGlobalT(globalT, maxT);
    return common.getPosition(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    LoopingUniformCRSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    floating_t wrappedT = SplineCommon::wrapGlobalT(globalT, maxT);
    return common.getTangent(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    LoopingUniformCRSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
    floating_t wrappedT = SplineCommon::wrapGlobalT(globalT, maxT);
    return common.getCurvature(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTCW
    LoopingUniformCRSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
    floating_t wrappedT = SplineCommon::wrapGlobalT(globalT, maxT);
    return common.getWiggle(wrappedT);
}

template<class InterpolationType, typename floating_t>
floating_t LoopingUniformCRSpline<InterpolationType,floating_t>::arcLength(floating_t a, floating_t b) const
{
    floating_t wrappedA =  SplineCommon::wrapGlobalT(a, maxT);
    floating_t wrappedB =  SplineCommon::wrapGlobalT(b, maxT);

    return ArcLength::arcLength(*this, wrappedA, wrappedB);
}
