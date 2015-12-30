#ifndef LOOPING_GENERIC_B_SPLINE
#define LOOPING_GENERIC_B_SPLINE


#include "spline_library/spline.h"
#include "spline_library/utils/spline_setup.h"

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

    floating_t getT(int index) const override { indexToT.at(index); }
    floating_t getMaxT(void) const override { return maxT; }

    bool isLooping(void) const override { return true; }

//methods
private:
    InterpolationType computeDeboor(size_t knotIndex, int degree, float globalT) const;
    InterpolationType computeDeboorDerivative(size_t knotIndex, int degree, float globalT, int derivativeLevel) const;
//data
private:
    std::vector<InterpolationType> positions;
    std::vector<floating_t> knots;

    floating_t maxT;

    int splineDegree;

    //map from index to t value. it's a map and not an array so we can store negative indexes
    std::unordered_map<int,floating_t> indexToT;
};

template<class InterpolationType, typename floating_t>
LoopingGenericBSpline<InterpolationType,floating_t>::LoopingGenericBSpline(const std::vector<InterpolationType> &points, int degree)
    :Spline<InterpolationType,floating_t>(points), splineDegree(degree)
{
    assert(points.size() > splineDegree);

    int size = points.size();
    int padding = degree - 1;

    //compute the T values for each point
    indexToT = SplineSetup::computeLoopingBSplineKnots(points, 0.0f, padding);
    maxT = indexToT[size];

    //we need enough space to repeat the last 'degree' elements
    positions = std::vector<InterpolationType>(points.size() + degree);
    positions[0] = points[size - 1];
    std::copy(points.begin(), points.end(), positions.begin() + 1);
    std::copy_n(points.begin(), padding, positions.end() -padding);

    knots = std::vector<floating_t>(indexToT.size() - 1);

    for(int i = -padding; i < size + padding + 1; i++)
    {
        knots[i + padding] = indexToT[i];
    }
}

template<class InterpolationType, typename floating_t>
InterpolationType LoopingGenericBSpline<InterpolationType,floating_t>::getPosition(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    size_t startIndex = SplineSetup::getIndexForT(knots, wrappedT);

    return computeDeboor(startIndex + 1, splineDegree, wrappedT);
}



template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPT
    LoopingGenericBSpline<InterpolationType,floating_t>::getTangent(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    size_t startIndex = SplineSetup::getIndexForT(knots, wrappedT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                computeDeboor(startIndex + 1, splineDegree, wrappedT),
                computeDeboorDerivative(startIndex + 1, splineDegree, wrappedT, 1)
                );
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTC
    LoopingGenericBSpline<InterpolationType,floating_t>::getCurvature(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    size_t startIndex = SplineSetup::getIndexForT(knots, wrappedT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                computeDeboor(startIndex + 1, splineDegree, wrappedT),
                computeDeboorDerivative(startIndex + 1, splineDegree, wrappedT, 1),
                computeDeboorDerivative(startIndex + 1, splineDegree, wrappedT, 2)
                );
}

template<class InterpolationType, typename floating_t>
typename Spline<InterpolationType,floating_t>::InterpolatedPTCW
    LoopingGenericBSpline<InterpolationType,floating_t>::getWiggle(floating_t globalT) const
{
    floating_t wrappedT = SplineSetup::wrapGlobalT(globalT, maxT);
    size_t startIndex = SplineSetup::getIndexForT(knots, wrappedT);

    return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                computeDeboor(startIndex + 1, splineDegree, wrappedT),
                computeDeboorDerivative(startIndex + 1, splineDegree, wrappedT, 1),
                computeDeboorDerivative(startIndex + 1, splineDegree, wrappedT, 2),
                computeDeboorDerivative(startIndex + 1, splineDegree, wrappedT, 3)
                );
}

template<class InterpolationType, typename floating_t>
InterpolationType LoopingGenericBSpline<InterpolationType,floating_t>::computeDeboor(size_t knotIndex, int degree, float globalT) const
{
    if(degree == 0)
    {
        return positions[knotIndex];
    }
    else
    {
        floating_t alpha = (globalT - knots[knotIndex - 1]) / (knots[knotIndex + splineDegree - degree] - knots[knotIndex - 1]);

        InterpolationType leftRecursive = computeDeboor(knotIndex - 1, degree - 1, globalT);
        InterpolationType rightRecursive = computeDeboor(knotIndex, degree - 1, globalT);

        InterpolationType blended = leftRecursive * (1 - alpha) + rightRecursive * alpha;

        return blended;
    }
}

template<class InterpolationType, typename floating_t>
InterpolationType LoopingGenericBSpline<InterpolationType,floating_t>::computeDeboorDerivative(size_t knotIndex, int degree, float globalT, int derivativeLevel) const
{
    if(degree == 0)
    {
        return InterpolationType();
    }
    else
    {
        float multiplier = degree / (knots[knotIndex + splineDegree - degree] - knots[knotIndex - 1]);

        if(derivativeLevel <= 1)
        {
            //once we reach this point we don't want to compute the derivative anymore
            return multiplier *
                    (computeDeboor(knotIndex, degree - 1, globalT)
                   - computeDeboor(knotIndex - 1, degree - 1, globalT)
                     );
        }
        else
        {
            //recursively call the derivative function to compute a higher derivative
            return multiplier *
                    (computeDeboorDerivative(knotIndex, degree - 1, globalT, derivativeLevel - 1)
                   - computeDeboorDerivative(knotIndex - 1, degree - 1, globalT, derivativeLevel - 1)
                     );
        }
    }
}

#endif // LOOPING_GENERIC_B_SPLINE

