#pragma once

#include <vector>

#include "../utils/spline_setup.h"
#include "../utils/calculus.h"

template<class InterpolationType, typename floating_t>
class NaturalSplineCommon
{
public:
    struct alignas(16) NaturalSplineSegment
    {
        InterpolationType a, c;
    };

    inline NaturalSplineCommon(void) = default;
    inline NaturalSplineCommon(std::vector<NaturalSplineSegment> segments, std::vector<floating_t> knots)
        :segments(std::move(segments)), knots(std::move(knots))
    {}

    inline InterpolationType getPosition(floating_t globalT) const
    {
        size_t segmentIndex = SplineSetup::getIndexForT(knots, globalT);
        if(segmentIndex >= knots.size() - 1)
            segmentIndex--;

        floating_t localT = globalT - knots[segmentIndex];
        floating_t tDiff = knots[segmentIndex + 1] - knots[segmentIndex];

        return computePosition(segmentIndex, tDiff, localT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        size_t segmentIndex = SplineSetup::getIndexForT(knots, globalT);
        if(segmentIndex >= knots.size() - 1)
            segmentIndex--;

        floating_t localT = globalT - knots[segmentIndex];
        floating_t tDiff = knots[segmentIndex + 1] - knots[segmentIndex];

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computePosition(segmentIndex, tDiff, localT),
                    computeTangent(segmentIndex, tDiff, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        size_t segmentIndex = SplineSetup::getIndexForT(knots, globalT);
        if(segmentIndex >= knots.size() - 1)
            segmentIndex--;

        floating_t localT = globalT - knots[segmentIndex];
        floating_t tDiff = knots[segmentIndex + 1] - knots[segmentIndex];

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    computePosition(segmentIndex, tDiff, localT),
                    computeTangent(segmentIndex, tDiff, localT),
                    computeCurvature(segmentIndex, tDiff, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        size_t segmentIndex = SplineSetup::getIndexForT(knots, globalT);
        if(segmentIndex >= knots.size() - 1)
            segmentIndex--;

        floating_t localT = globalT - knots[segmentIndex];
        floating_t tDiff = knots[segmentIndex + 1] - knots[segmentIndex];

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    computePosition(segmentIndex, tDiff, localT),
                    computeTangent(segmentIndex, tDiff, localT),
                    computeCurvature(segmentIndex, tDiff, localT),
                    computeWiggle(segmentIndex, tDiff)
                    );
    }

    inline floating_t getLength(floating_t a, floating_t b) const
    {
        //get the knot indices for the beginning and end
        size_t numSegments = knots.size() - 1;
        size_t aIndex = SplineSetup::getIndexForT(knots, a);
        if(aIndex >= numSegments)
            aIndex--;
        size_t bIndex = SplineSetup::getIndexForT(knots, b);
        if(bIndex >= numSegments)
            bIndex--;

        //if a and b occur inside the same segment, compute the length within that segment
        //but excude cases where a > b, because that means we need to wrap around
        if(aIndex == bIndex && a <= b) {
            return computeSegmentLength(aIndex, a - knots[aIndex], b - knots[aIndex]);
        }
        else {
            //a and b occur in different segments, so compute one length for every segment
            floating_t result{0};

            //first segment
            result += computeSegmentLength(aIndex, a - knots[aIndex], knots[aIndex + 1] - knots[aIndex]);

            //last segment
            result += computeSegmentLength(bIndex, 0, b - knots[bIndex]);

            //if b index is less than a index, that means the user wants to wrap around the end of the spline and back to the beginning
            //if so, add the number of points in the spline to bIndex, and we'll use mod to make sure it stays in range
            if(bIndex <= aIndex)
                bIndex += numSegments;

            //middle segments
            for(size_t i = aIndex + 1; i < bIndex; i++) {
                floating_t segmentEnd = knots[(i + 1)%knots.size()] - knots[i%knots.size()];
                result += computeSegmentLength(i%numSegments, 0, segmentEnd);
            }

            return result;
        }
    }

    inline floating_t getTotalLength(void) const
    {
        floating_t result{0};
        for(size_t i = 0; i < knots.size() - 1; i++) {
            floating_t segmentEnd = knots[i + 1] - knots[i];
            result += computeSegmentLength(i, 0, segmentEnd);
        }
        return result;
    }

private: //methods
    inline InterpolationType computePosition(size_t index, floating_t tDiff, floating_t t) const
    {
        auto b = computeB(index, tDiff);
        auto d = computeD(index, tDiff);

        return segments[index].a + t * (b + t * (segments[index].c + t * d));
    }

    inline InterpolationType computeTangent(size_t index, floating_t tDiff, floating_t t) const
    {
        auto b = computeB(index, tDiff);
        auto d = computeD(index, tDiff);

        //compute the derivative of the position function
        return b + t * (floating_t(2) * segments[index].c + (3 * t) * d);
    }

    inline InterpolationType computeCurvature(size_t index, floating_t tDiff, floating_t t) const
    {
        auto d = computeD(index, tDiff);

        //compute the 2nd derivative of the position function
        return floating_t(2) * segments[index].c + (6 * t) * d;
    }

    inline InterpolationType computeWiggle(size_t index, floating_t tDiff) const
    {
        auto d = computeD(index, tDiff);

        //compute the 3rd derivative of the position function
        return floating_t(6) * d;
    }


    //B is the tangent at t=0 for a segment, and D is effectively the wiggle for a segment
    //we COULD precompute these and store them alongside a and c in the segment
    //testing shows that it's faster (because of cache, and pipelining, etc) to just recompute them every time
    inline InterpolationType computeB(size_t index, floating_t tDiff) const
    {
        return (segments[index+1].a - segments[index].a) / tDiff - (tDiff / 3) * (segments[index+1].c + floating_t(2)*segments[index].c);
    }
    inline InterpolationType computeD(size_t index, floating_t tDiff) const
    {
        return (segments[index+1].c - segments[index].c) / (3 * tDiff);
    }



    inline floating_t computeSegmentLength(size_t index, floating_t from, floating_t to) const
    {
        floating_t tDiff = knots[index + 1] - knots[index];
        auto segmentFunction = [=](floating_t t) -> floating_t {
            auto tangent = computeTangent(index, tDiff, t);
            return tangent.length();
        };

        return SplineLibraryCalculus::gaussLegendreQuadratureIntegral(segmentFunction, from, to);
    }

private: //data
    std::vector<NaturalSplineSegment> segments;
    std::vector<floating_t> knots;
};
