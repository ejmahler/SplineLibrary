#ifndef LOOPING_GENERIC_B_SPLINE
#define LOOPING_GENERIC_B_SPLINE


#include "spline_library/spline.h"
#include "spline_library/utils/spline_setup.h"
#include "generic_b_spline_common.h"

#include <unordered_map>

template<class InterpolationType, typename floating_t=float>
class LoopingGenericBSpline final : public Spline<InterpolationType, floating_t>
{
//constructors
public:
    LoopingGenericBSpline(const std::vector<InterpolationType> &points, int degree);

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
private:
    GenericBSplineCommon<InterpolationType, floating_t> common;

    floating_t maxT;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
LoopingGenericBSpline<InterpolationType,floating_t>::LoopingGenericBSpline(const std::vector<InterpolationType> &points, int degree)
    :Spline<InterpolationType,floating_t>(points)
{
    assert(points.size() > size_t(degree));

    int size = points.size();
    int padding = degree - 1;

    //compute the T values for each point
    indexToT = SplineSetup::computeLoopingTValues(points, 0.0f, padding);
    maxT = indexToT[size];

    //we need enough space to repeat the last 'degree' elements
    std::vector<InterpolationType> positions(points.size() + degree);

    //it would be easiest to just copy the points vector to the position vector, then copy the first 'degree' elements again
    //this DOES work, but interpolation begins in the wrong place (ie getPosition(0) occurs at the wrong place on the spline)
    //to fix this, we effectively "rotate" the position vector backwards, by copying point[size-1] to the beginning
    //then copying the points vector in after, then copying degree-1 elements from the beginning
    positions[0] = points[size - 1];
    std::copy(points.begin(), points.end(), positions.begin() + 1);
    std::copy_n(points.begin(), padding, positions.end() - padding);

    //the index to t calculation computes an extra T value that we don't need, we just discard it here when creating the knot vector
    std::vector<floating_t> knots(indexToT.size() - 1);

    //for purposes of actual interpolation, we don't need the negative indexes found in indexToT
    //so we're going to add the minimum possible value to every entry and stick them in a vector
    for(int i = -padding; i < size + padding + 1; i++)
    {
        knots[i + padding] = indexToT[i];
    }

    common = GenericBSplineCommon<InterpolationType, floating_t>(std::move(positions), std::move(knots), degree);
}

template<class InterpolationType, typename floating_t>
InterpolationType LoopingGenericBSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getPosition(wrappedT);
}



template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    LoopingGenericBSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getTangent(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    LoopingGenericBSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getCurvature(wrappedT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTCW
    LoopingGenericBSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    return common.getWiggle(wrappedT);
}


#endif // LOOPING_GENERIC_B_SPLINE

