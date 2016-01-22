#pragma once

#include <vector>

#include "../utils/calculus.h"

template<class InterpolationType, typename floating_t>
class UniformCubicBSplineCommon
{
public:

    inline UniformCubicBSplineCommon(void) = default;
    inline UniformCubicBSplineCommon(std::vector<InterpolationType> points)
        :points(std::move(points))
    {}

    inline InterpolationType getPosition(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        floating_t localT = globalT - knotIndex;

        return computePosition(knotIndex, localT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        floating_t localT = globalT - knotIndex;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computePosition(knotIndex, localT),
                    computeTangent(knotIndex, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        floating_t localT = globalT - knotIndex;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    computePosition(knotIndex, localT),
                    computeTangent(knotIndex, localT),
                    computeCurvature(knotIndex, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        floating_t localT = globalT - knotIndex;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    computePosition(knotIndex, localT),
                    computeTangent(knotIndex, localT),
                    computeCurvature(knotIndex, localT),
                    computeWiggle(knotIndex)
                    );
    }

    inline floating_t getLength(floating_t a, floating_t b) const
    {
        //get the knot indices for the beginning and end
        size_t aIndex = size_t(a);
        size_t bIndex = size_t(b);

        size_t numSegments = points.size() - 4;

        if(aIndex > numSegments)
            aIndex = numSegments;
        if(bIndex > numSegments)
            bIndex = numSegments;

        //if a and b occur inside the same segment, compute the length within that segment
        //but excude cases where a > b, because that means we need to wrap around
        if(aIndex == bIndex && a <= b) {
            return computeSegmentLength(aIndex, a - aIndex, b - aIndex);
        }
        else {
            //a and b occur in different segments, so compute one length for every segment
            floating_t result{0};

            //first segment
            result += computeSegmentLength(aIndex, a - aIndex, 1);

            //last segment
            result += computeSegmentLength(bIndex, 0, b - bIndex);

            //if b index is less than a index, that means the user wants to wrap around the end of the spline and back to the beginning
            //if so, add the number of points in the spline to bIndex, and we'll use mod to make sure it stays in range
            if(bIndex <= aIndex)
                bIndex += numSegments;

            //middle segments
            for(size_t i = aIndex + 1; i < bIndex; i++) {
                result += computeSegmentLength(i%numSegments, 0, 1);
            }

            return result;
        }
    }

    inline floating_t getTotalLength(void) const
    {
        floating_t result{0};
        for(size_t i = 0; i <= points.size() - 4; i++) {
            result += computeSegmentLength(i, 0, 1);
        }
        return result;
    }

private: //methods
    inline InterpolationType computePosition(size_t index, floating_t t) const
    {
        return (
                    points[index] * ((1 - t) * (1 - t) * (1 - t)) +
                    points[index + 1] * (t * t * 3 * (t - 2) + 4) +
                    points[index + 2] * (t * (t * (-3 * t + 3) + 3) + 1) +
                    points[index + 3] * (t * t * t)
                ) / 6;
    }

    inline InterpolationType computeTangent(size_t index, floating_t t) const
    {
        return (
                    points[index] * (-(1 - t) * (1 - t)) +
                    points[index + 1] * (t * (3 * t - 4)) +
                    points[index + 2] * ((3 * t + 1) * (1 - t)) +
                    points[index + 3] * (t * t)
                ) / 2;
    }

    inline InterpolationType computeCurvature(size_t index, floating_t t) const
    {
        return (
                    points[index] * (1 - t) +
                    points[index + 1] * (3 * t - 2) +
                    points[index + 2] * (1 - 3 * t) +
                    points[index + 3] * (t)
                );
    }

    inline InterpolationType computeWiggle(size_t index) const
    {
        return 3 * (points[index + 1] - points[index + 2]) + (points[index + 3] - points[index]);
    }

    inline floating_t computeSegmentLength(size_t index, floating_t from, floating_t to) const
    {
        auto segmentFunction = [this, index](floating_t t) -> floating_t {
            auto tangent = computeTangent(index, t);
            return tangent.length();
        };

        return SplineLibraryCalculus::adaptiveSimpsonsIntegral(segmentFunction, from, to);
    }

private: //data
    std::vector<InterpolationType> points;
};
