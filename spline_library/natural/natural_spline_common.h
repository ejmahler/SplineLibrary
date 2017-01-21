#pragma once

#include <vector>

#include "../utils/spline_common.h"
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

    inline size_t segmentCount(void) const
    {
        return segments.size() - 1;
    }

    inline size_t segmentForT(floating_t t) const
    {
        size_t segmentIndex = SplineCommon::getIndexForT(knots, t);
        if(segmentIndex >= segmentCount())
            return segmentCount() - 1;
        else
            return segmentIndex;
    }

    inline floating_t segmentT(size_t segmentIndex) const
    {
        return knots[segmentIndex];
    }

    inline InterpolationType getPosition(floating_t globalT) const
    {
        size_t segmentIndex = SplineCommon::getIndexForT(knots, globalT);
        if(segmentIndex >= knots.size() - 1)
            segmentIndex--;

        floating_t localT = globalT - knots[segmentIndex];
        floating_t tDiff = knots[segmentIndex + 1] - knots[segmentIndex];

        return computePosition(segmentIndex, tDiff, localT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        size_t segmentIndex = SplineCommon::getIndexForT(knots, globalT);
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
        size_t segmentIndex = SplineCommon::getIndexForT(knots, globalT);
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
        size_t segmentIndex = SplineCommon::getIndexForT(knots, globalT);
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

    inline floating_t segmentLength(size_t segmentIndex, floating_t a, floating_t b) const {

        floating_t tDiff = knots[segmentIndex + 1] - knots[segmentIndex];
        auto segmentFunction = [=](floating_t t) -> floating_t {
            auto tangent = computeTangent(segmentIndex, tDiff, t);
            return tangent.length();
        };

        floating_t localA = a - knots[segmentIndex];
        floating_t localB = b - knots[segmentIndex];

        return SplineLibraryCalculus::gaussLegendreQuadratureIntegral(segmentFunction, localA, localB);
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

    }

private: //data
    std::vector<NaturalSplineSegment> segments;
    std::vector<floating_t> knots;
};
