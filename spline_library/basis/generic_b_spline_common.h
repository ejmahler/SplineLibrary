#pragma once

#include <vector>

#include "../utils/spline_common.h"

template<class InterpolationType, typename floating_t>
class GenericBSplineCommon
{
public:
    inline GenericBSplineCommon(void) = default;
    inline GenericBSplineCommon(std::vector<InterpolationType> positions, std::vector<floating_t> knots, int splineDegree)
        :positions(std::move(positions)), knots(std::move(knots)), splineDegree(splineDegree)
    {}

    inline size_t segmentCount(void) const
    {
        return positions.size() - splineDegree;
    }

    inline size_t segmentForT(floating_t t) const
    {
        if(t < 0) {
            return 0;
        }

        size_t segmentIndex = SplineCommon::getIndexForT(knots, t) - (splineDegree - 1);
        if(segmentIndex > segmentCount() - 1)
        {
            return segmentCount() - 1;
        }
        else
        {
            return segmentIndex;
        }
    }

    inline floating_t segmentT(size_t segmentIndex) const
    {
        return knots[segmentIndex + splineDegree - 1];
    }

    inline InterpolationType getPosition(floating_t globalT) const
    {
        size_t startIndex = SplineCommon::getIndexForT(knots, globalT);
        if(startIndex >= positions.size() - 2)
        {
            startIndex = positions.size() - 2;
        }

        return computeDeboor(startIndex + 1, splineDegree, globalT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        size_t startIndex = SplineCommon::getIndexForT(knots, globalT);
        if(startIndex >= positions.size() - 2)
        {
            startIndex = positions.size() - 2;
        }

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computeDeboor(startIndex + 1, splineDegree, globalT),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 1)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        size_t startIndex = SplineCommon::getIndexForT(knots, globalT);
        if(startIndex >= positions.size() - 2)
        {
            startIndex = positions.size() - 2;
        }

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    computeDeboor(startIndex + 1, splineDegree, globalT),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 1),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 2)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        size_t startIndex = SplineCommon::getIndexForT(knots, globalT);
        if(startIndex >= positions.size() - 2)
        {
            startIndex = positions.size() - 2;
        }

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    computeDeboor(startIndex + 1, splineDegree, globalT),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 1),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 2),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 3)
                    );
    }

    inline floating_t segmentLength(size_t segmentIndex, floating_t a, floating_t b) const {

        auto innerIndex = segmentIndex + splineDegree - 1;
        auto tDistance = knots[innerIndex + 1] - knots[innerIndex];

        //it's perfectly legal for Bspline segments to have a T distance of 0, in which case the arc length is 0
        if(tDistance > 0)
        {
            return computeSegmentLength(innerIndex, knots[innerIndex] + a * tDistance, knots[innerIndex] + b * tDistance);
        }
        else
        {
            return 0;
        }
    }

private: //methods
    InterpolationType computeDeboor(size_t knotIndex, int degree, float globalT) const;
    InterpolationType computeDeboorDerivative(size_t knotIndex, int degree, float globalT, int derivativeLevel) const;

    inline floating_t computeSegmentLength(size_t index, floating_t from, floating_t to) const
    {
        auto segmentFunction = [this, index](floating_t t) -> floating_t {
            auto tangent = computeDeboorDerivative(index + 1, splineDegree, t, 1);
            return tangent.length();
        };

        return SplineLibraryCalculus::gaussLegendreQuadratureIntegral(segmentFunction, from, to);
    }

private: //data
    std::vector<InterpolationType> positions;
    std::vector<floating_t> knots;
    int splineDegree;
};

template<class InterpolationType, typename floating_t>
InterpolationType GenericBSplineCommon<InterpolationType,floating_t>::computeDeboor(size_t knotIndex, int degree, float globalT) const
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
InterpolationType GenericBSplineCommon<InterpolationType,floating_t>::computeDeboorDerivative(size_t knotIndex, int degree, float globalT, int derivativeLevel) const
{
    if(degree == 0)
    {
        //if we hit degree 0 before derivative level 0, then this spline's
        //degree isn't high enough to support whatever derivative level was requested
        return InterpolationType();
    }
    else
    {
        floating_t multiplier = degree / (knots[knotIndex + splineDegree - degree] - knots[knotIndex - 1]);

        if(derivativeLevel <= 1)
        {
            //once we reach this point, the derivative calculation is "complete"
            //in that from here, we go back to the normal deboor calculation deeper in the recursive tree
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
