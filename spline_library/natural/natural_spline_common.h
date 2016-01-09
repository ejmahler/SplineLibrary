#ifndef NATURAL_SPLINE_COMMON_H
#define NATURAL_SPLINE_COMMON_H

#include <vector>

#include "../utils/spline_setup.h"

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
        return b + t * (2 * segments[index].c + (3 * t) * d);
    }

    inline InterpolationType computeCurvature(size_t index, floating_t tDiff, floating_t t) const
    {
        auto d = computeD(index, tDiff);

        //compute the 2nd derivative of the position function
        return 2 * segments[index].c + (6 * t) * d;
    }

    inline InterpolationType computeWiggle(size_t index, floating_t tDiff) const
    {
        auto d = computeD(index, tDiff);

        //compute the 3rd derivative of the position function
        return 6 * d;
    }


    //B is the tangent at t=0 for a segment, and D is effectively the wiggle for a segment
    //we COULD precompute these and store them alongside a and c in the segment
    //testing shows that it's faster (because of cache, and pipelining, etc) to just recompute them every time
    inline InterpolationType computeB(size_t index, floating_t tDiff) const
    {
        return (segments[index+1].a - segments[index].a) / tDiff - (tDiff / 3) * (segments[index+1].c + 2*segments[index].c);
    }
    inline InterpolationType computeD(size_t index, floating_t tDiff) const
    {
        return (segments[index+1].c - segments[index].c) / (3 * tDiff);
    }

private: //data
    std::vector<NaturalSplineSegment> segments;
    std::vector<floating_t> knots;
};




#endif // NATURAL_SPLINE_COMMON_H

