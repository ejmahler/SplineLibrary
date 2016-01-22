#pragma once

#include <vector>

#include "../../utils/calculus.h"

template<class InterpolationType, typename floating_t>
class UniformCRSplineCommon
{
public:

    inline UniformCRSplineCommon(void) = default;
    inline UniformCRSplineCommon(std::vector<InterpolationType> points)
        :points(std::move(points))
    {}

    inline InterpolationType getPosition(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        floating_t localT = globalT - knotIndex;

        return computePosition(knotIndex + 1, localT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        size_t knotIndex = size_t(globalT);

        //make sure the knot index stays in-bounds
        if(knotIndex > points.size() - 4)
            knotIndex = points.size() - 4;

        floating_t localT = globalT - knotIndex;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computePosition(knotIndex + 1, localT),
                    computeTangent(knotIndex + 1, localT)
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
                    computePosition(knotIndex + 1, localT),
                    computeTangent(knotIndex + 1, localT),
                    computeCurvature(knotIndex + 1, localT)
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
                    computePosition(knotIndex + 1, localT),
                    computeTangent(knotIndex + 1, localT),
                    computeCurvature(knotIndex + 1, localT),
                    computeWiggle(knotIndex + 1)
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
        auto beforeTangent = computeTangentAtIndex(index);
        auto afterTangent = computeTangentAtIndex(index + 1);

        auto oneMinusT = 1 - t;

        auto basis00 = (1 + 2*t) * oneMinusT * oneMinusT;
        auto basis10 = t * oneMinusT * oneMinusT;

        auto basis11 = t * t * -oneMinusT;
        auto basis01 = t * t * (3 - 2*t);

        return
                basis00 * points[index] +
                basis10 * beforeTangent +

                basis11 * afterTangent +
                basis01 * points[index + 1];
    }

    inline InterpolationType computeTangent(size_t index, floating_t t) const
    {
        auto beforeTangent = computeTangentAtIndex(index);
        auto afterTangent = computeTangentAtIndex(index + 1);

        auto oneMinusT = 1 - t;

        auto d_basis00 = 6 * t * (t - 1);
        auto d_basis10 = (1 - 3*t) * oneMinusT;

        auto d_basis11 = t * (3 * t - 2);
        auto d_basis01 = -d_basis00;

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the derivative of the position function and nothing else
        //if you know why please let me know
        return
                d_basis00 * points[index] +
                d_basis10 * beforeTangent +

                d_basis11 * afterTangent +
                d_basis01 * points[index + 1];
    }

    inline InterpolationType computeCurvature(size_t index, floating_t t) const
    {
        auto beforeTangent = computeTangentAtIndex(index);
        auto afterTangent = computeTangentAtIndex(index + 1);

        auto d2_basis00 = 6 * (2 * t - 1);
        auto d2_basis10 = 2 * (3 * t - 2);

        auto d2_basis11 = 2 * (3 * t - 1);
        auto d2_basis01 = -d2_basis00;

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the 2nd derivative of the position function and nothing else
        //if you know why please let me know
        return
                d2_basis00 * points[index] +
                d2_basis10 * beforeTangent +

                d2_basis11 * afterTangent +
                d2_basis01 * points[index + 1];
    }

    inline InterpolationType computeWiggle(size_t index) const
    {
        auto beforeTangent = computeTangentAtIndex(index);
        auto afterTangent = computeTangentAtIndex(index + 1);

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the 2nd derivative of the position function and nothing else
        //if you know why please let me know
        return 12 * (points[index] - points[index + 1]) + 6 * (beforeTangent + afterTangent);
    }

    inline InterpolationType computeTangentAtIndex(size_t i) const
    {
        return (points[i + 1] - points[i - 1]) / 2;
    }


    inline floating_t computeSegmentLength(size_t index, floating_t from, floating_t to) const
    {
        auto segmentFunction = [this, index](floating_t t) -> floating_t {
            auto tangent = computeTangent(index + 1, t);
            return tangent.length();
        };

        return SplineLibraryCalculus::adaptiveSimpsonsIntegral(segmentFunction, from, to);
    }

private: //data
    std::vector<InterpolationType> points;
};
