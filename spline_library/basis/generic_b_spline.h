#pragma once

#include "../spline.h"
#include "../utils/spline_common.h"
#include "../arclength.h"

#include "generic_b_spline_common.h"

#include <unordered_map>

template<class InterpolationType, typename floating_t=float>
class GenericBSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    GenericBSpline(const std::vector<InterpolationType> &points, int degree);

//methods
public:
    InterpolationType getPosition(floating_t t) const override { return common.getPosition(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t t) const override { return common.getTangent(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t t) const override { return common.getCurvature(t); }
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t t) const override { return common.getWiggle(t); }

    floating_t arcLength(floating_t a, floating_t b) const override { if(a > b) std::swap(a,b); return ArcLength::arcLength(*this,a,b); }
    floating_t totalLength(void) const override { return ArcLength::totalLength(*this); }

    floating_t getT(int index) const override { return indexToT.at(index); }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return false; }

    size_t segmentCount(void) const override { return common.segmentCount(); }
    size_t segmentForT(floating_t t) const override { return common.segmentForT(t); }
    floating_t segmentT(size_t segmentIndex) const override { return common.segmentT(segmentIndex); }
    floating_t segmentArcLength(size_t segmentIndex, floating_t a, floating_t b) const override { return common.segmentLength(segmentIndex, a, b); }

//methods
private:
    InterpolationType computeDeboor(size_t knotIndex, int degree, float globalT) const;
    InterpolationType computeDeboorDerivative(size_t knotIndex, int degree, float globalT, int derivativeLevel) const;

//data
private:
    GenericBSplineCommon<InterpolationType, floating_t> common;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
GenericBSpline<InterpolationType,floating_t>::GenericBSpline(const std::vector<InterpolationType> &points, int degree)
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() > size_t(degree));

    int size = points.size();
    int padding = degree - 1;

    //compute the T values for each point
    indexToT = SplineCommon::computeTValuesWithOuterPadding(points, 0.0f, padding);
    maxT = indexToT[size - degree];

    //for purposes of actual interpolation, we don't need the negative indexes found in indexToT
    //so we're going to add the minimum possible value to every entry and stick them in a vector
    std::vector<floating_t> knots = std::vector<floating_t>(indexToT.size());
    for(int i = -padding; i < size + padding; i++)
    {
        knots[i + padding] = indexToT[i];
    }

    common = GenericBSplineCommon<InterpolationType, floating_t>(points, std::move(knots), degree);
}
