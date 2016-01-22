#pragma once

#include <vector>

#include "../utils/spline_setup.h"

template<class InterpolationType, typename floating_t>
class GenericBSplineCommon
{
public:
    inline GenericBSplineCommon(void) = default;
    inline GenericBSplineCommon(std::vector<InterpolationType> positions, std::vector<floating_t> knots, int splineDegree)
        :positions(std::move(positions)), knots(std::move(knots)), splineDegree(splineDegree)
    {}

    inline InterpolationType getPosition(floating_t globalT) const
    {
        size_t startIndex = SplineSetup::getIndexForT(knots, globalT);
        return computeDeboor(startIndex + 1, splineDegree, globalT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        size_t startIndex = SplineSetup::getIndexForT(knots, globalT);

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computeDeboor(startIndex + 1, splineDegree, globalT),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 1)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        size_t startIndex = SplineSetup::getIndexForT(knots, globalT);

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    computeDeboor(startIndex + 1, splineDegree, globalT),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 1),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 2)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        size_t startIndex = SplineSetup::getIndexForT(knots, globalT);

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    computeDeboor(startIndex + 1, splineDegree, globalT),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 1),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 2),
                    computeDeboorDerivative(startIndex + 1, splineDegree, globalT, 3)
                    );
    }

    inline floating_t getLength(floating_t a, floating_t b) const
    {
        //get the knot indices for the beginning and end
        size_t aIndex = SplineSetup::getIndexForT(knots, a);
        size_t bIndex = SplineSetup::getIndexForT(knots, b);

        size_t numSegments = positions.size() - (splineDegree - 1);

        //if a and b occur inside the same segment, compute the length within that segment
        //but excude cases where a > b, because that means we need to wrap around
        if(aIndex == bIndex && a <= b) {
            return computeSegmentLength(aIndex, a - aIndex, b - aIndex);
        }
        else {
            //a and b occur in different segments, so compute one length for every segment
            floating_t result{0};

            //first segment
            result += computeSegmentLength(aIndex, a, knots[aIndex + 1]);

            //last segment
            result += computeSegmentLength(bIndex, knots[bIndex], b);

            //if b index is less than a index, that means the user wants to wrap around the end of the spline and back to the beginning
            //if so, add the number of points in the spline to bIndex, and we'll use mod to make sure it stays in range
            if(bIndex <= aIndex)
                bIndex += numSegments;

            //middle segments
            auto padding = splineDegree - 1;
            for(size_t i = aIndex + 1; i < bIndex; i++) {
                size_t wrappedIndex = (i - padding)%numSegments + padding;
                result += computeSegmentLength(wrappedIndex, knots[wrappedIndex], knots[wrappedIndex + 1]);
            }

            return result;
        }
    }

    inline floating_t getTotalLength(void) const
    {
        floating_t result{0};
        size_t numSegments = positions.size() - splineDegree;
        for(size_t segment = 0; segment < numSegments; segment++) {
            auto index = segment + splineDegree - 1;
            result += computeSegmentLength(index, knots[index], knots[index + 1]);
        }
        return result;
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

        return SplineLibraryCalculus::adaptiveSimpsonsIntegral(segmentFunction, from, to);
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
        float multiplier = degree / (knots[knotIndex + splineDegree - degree] - knots[knotIndex - 1]);

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
