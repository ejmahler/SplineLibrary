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
        InterpolationType a, b, c, d;

        inline InterpolationType computePosition(floating_t t) const
        {
            return a + t * (b + t * (c + t * d));
        }

        inline InterpolationType computeTangent(floating_t t) const
        {
            //compute the derivative of the position function
            return b + t * (2 * c + (3 * t) * d);
        }

        inline InterpolationType computeCurvature(floating_t t) const
        {
            //compute the 2nd derivative of the position function
            return 2 * c + (6 * t) * d;
        }

        inline InterpolationType computeWiggle(void) const
        {
            //compute the 3rd derivative of the position function
            return 6 * d;
        }
    };

    inline NaturalSplineCommon(void) = default;
    inline NaturalSplineCommon(std::vector<NaturalSplineSegment> segments, std::vector<floating_t> knots)
        :segments(std::move(segments)), knots(std::move(knots))
    {}

    inline InterpolationType getPosition(floating_t globalT) const
    {
        size_t segmentIndex = SplineSetup::getIndexForT(knots, globalT);
        floating_t localT = globalT - knots[segmentIndex];
        return segments[segmentIndex].computePosition(localT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        size_t segmentIndex = SplineSetup::getIndexForT(knots, globalT);
        floating_t localT = globalT - knots[segmentIndex];

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    segments[segmentIndex].computePosition(localT),
                    segments[segmentIndex].computeTangent(localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        size_t segmentIndex = SplineSetup::getIndexForT(knots, globalT);
        floating_t localT = globalT - knots[segmentIndex];

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    segments[segmentIndex].computePosition(localT),
                    segments[segmentIndex].computeTangent(localT),
                    segments[segmentIndex].computeCurvature(localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        size_t segmentIndex = SplineSetup::getIndexForT(knots, globalT);
        floating_t localT = globalT - knots[segmentIndex];

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    segments[segmentIndex].computePosition(localT),
                    segments[segmentIndex].computeTangent(localT),
                    segments[segmentIndex].computeCurvature(localT),
                    segments[segmentIndex].computeWiggle()
                    );
    }

private: //data
    std::vector<NaturalSplineSegment> segments;
    std::vector<floating_t> knots;
};




#endif // NATURAL_SPLINE_COMMON_H

