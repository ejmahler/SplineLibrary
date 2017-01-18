#pragma once

#include <cassert>

#include "../../spline.h"
#include "uniform_cr_spline_common.h"
#include "../../arclength.h"

template<class InterpolationType, typename floating_t=float>
class UniformCRSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    UniformCRSpline(const std::vector<InterpolationType> &points);

//methods
public:
    InterpolationType getPosition(floating_t t) const override { return common.getPosition(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t t) const override { return common.getTangent(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t t) const override { return common.getCurvature(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t t) const override { return common.getWiggle(t); }

    floating_t arcLength(floating_t a, floating_t b) const override { if(a > b) std::swap(a,b); return ArcLength::arcLength(*this,a,b); }
    floating_t totalLength(void) const override { return ArcLength::totalLength(*this); }

    floating_t getT(int index) const override { return index - 1; }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return false; }

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
UniformCRSpline<InterpolationType,floating_t>::UniformCRSpline(const std::vector<InterpolationType> &points)
    :Spline<InterpolationType,floating_t>(points), common(points), maxT(points.size() - 3)
{
    assert(points.size() >= 4);
}
