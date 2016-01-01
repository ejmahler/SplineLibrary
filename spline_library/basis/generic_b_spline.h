#ifndef GENERIC_B_SPLINE
#define GENERIC_B_SPLINE


#include "spline_library/spline.h"
#include "spline_library/utils/spline_setup.h"

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
    InterpolationType getPosition(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t x) const override;
    typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t x) const override;

    floating_t getT(int index) const override { return indexToT.at(index); }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return false; }

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
    indexToT = SplineSetup::computeTValuesWithOuterPadding(points, 0.0f, padding);
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

template<class InterpolationType, typename floating_t>
InterpolationType GenericBSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    return common.getPosition(globalT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    GenericBSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    return common.getTangent(globalT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    GenericBSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
    return common.getCurvature(globalT);
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTCW
    GenericBSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
    return common.getWiggle(globalT);
}

#endif // GENERIC_B_SPLINE

