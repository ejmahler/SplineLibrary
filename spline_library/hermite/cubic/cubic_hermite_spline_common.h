#pragma once

#include <vector>

#include "../../utils/calculus.h"
#include "../../utils/spline_setup.h"

template<class InterpolationType, typename floating_t>
class CubicHermiteSplineCommon
{
public:
    struct alignas(8) CubicHermiteSplinePoint
    {
        InterpolationType position, tangent;
    };

    inline CubicHermiteSplineCommon(void) = default;
    inline CubicHermiteSplineCommon(std::vector<CubicHermiteSplinePoint> points, std::vector<floating_t> knots)
        :points(std::move(points)), knots(std::move(knots))
    {}

    inline InterpolationType getPosition(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return computePosition(knotIndex, tDiff, localT);
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPT getTangent(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPT(
                    computePosition(knotIndex, tDiff, localT),
                    computeTangent(knotIndex, tDiff, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTC getCurvature(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTC(
                    computePosition(knotIndex, tDiff, localT),
                    computeTangent(knotIndex, tDiff, localT),
                    computeCurvature(knotIndex, tDiff, localT)
                    );
    }

    inline typename Spline<InterpolationType,floating_t>::InterpolatedPTCW getWiggle(floating_t globalT) const
    {
        //get the knot index. if it's the final knot, back it up by one
        size_t knotIndex = SplineSetup::getIndexForT(knots, globalT);
        if(knotIndex >= knots.size() - 1)
            knotIndex--;

        floating_t tDiff = (knots[knotIndex + 1] - knots[knotIndex]);
        floating_t localT = (globalT - knots[knotIndex]) / tDiff;

        return typename Spline<InterpolationType,floating_t>::InterpolatedPTCW(
                    computePosition(knotIndex, tDiff, localT),
                    computeTangent(knotIndex, tDiff, localT),
                    computeCurvature(knotIndex, tDiff, localT),
                    computeWiggle(knotIndex, tDiff)
                    );
    }


    inline floating_t getLength(floating_t a, floating_t b) const
    {
        //get the knot indices for the beginning and end
        size_t aIndex = SplineSetup::getIndexForT(knots, a);
        size_t bIndex = SplineSetup::getIndexForT(knots, b);

        if(aIndex > knots.size() - 2)
            aIndex = knots.size() - 2;
        if(bIndex > knots.size() - 2)
            bIndex = knots.size() - 2;

        floating_t aPercent = (a - knots[aIndex]) / (knots[aIndex + 1] - knots[aIndex]);
        floating_t bPercent = (b - knots[bIndex]) / (knots[bIndex + 1] - knots[bIndex]);

        //if a and b occur inside the same segment, compute the length within that segment
        //but excude cases where a > b, because that means we need to wrap around
        if(aIndex == bIndex && a <= b) {
            return computeSegmentLength(aIndex, aPercent, bPercent);
        }
        else {
            //a and b occur in different segments, so compute one length for every segment
            floating_t result{0};

            //first segment
            result += computeSegmentLength(aIndex, aPercent, 1);;

            //last segment
            result += computeSegmentLength(bIndex, 0, bPercent);

            //if b index is less than a index, that means the user wants to wrap around the end of the spline and back to the beginning
            //if so, add the number of points in the spline to bIndex, and we'll use mod to make sure it stays in range
            size_t numSegments = knots.size() - 1;
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
        for(size_t i = 0; i < knots.size() - 1; i++) {
            result += computeSegmentLength(i, 0, 1);
        }
        return result;
    }


private: //methods
    inline InterpolationType computePosition(size_t index, floating_t tDiff, floating_t t) const
    {
        auto oneMinusT = 1 - t;

        auto basis00 = (1 + 2*t) * oneMinusT * oneMinusT;
        auto basis10 = t * oneMinusT * oneMinusT;

        auto basis11 = t * t * -oneMinusT;
        auto basis01 = t * t * (3 - 2*t);

        return
                basis00 * points[index].position +
                basis10 * tDiff * points[index].tangent +

                basis11 * tDiff * points[index + 1].tangent +
                basis01 * points[index + 1].position;
    }

    inline InterpolationType computeTangent(size_t index, floating_t tDiff, floating_t t) const
    {
        auto oneMinusT = 1 - t;

        auto d_basis00 = 6 * t * (t - 1);
        auto d_basis10 = (1 - 3*t) * oneMinusT;

        auto d_basis11 = t * (3 * t - 2);
        auto d_basis01 = -d_basis00;

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the derivative of the position function and nothing else
        //if you know why please let me know
        return (
                d_basis00 * points[index].position +
                d_basis10 * tDiff * points[index].tangent +

                d_basis11 * tDiff * points[index + 1].tangent +
                d_basis01 * points[index + 1].position
                ) / tDiff;
    }

    inline InterpolationType computeCurvature(size_t index, floating_t tDiff, floating_t t) const
    {
        auto d2_basis00 = 6 * (2 * t - 1);
        auto d2_basis10 = 2 * (3 * t - 2);

        auto d2_basis11 = 2 * (3 * t - 1);
        auto d2_basis01 = -d2_basis00;

        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the 2nd derivative of the position function and nothing else
        //if you know why please let me know
        return (
                d2_basis00 * points[index].position +
                d2_basis10 * tDiff * points[index].tangent +

                d2_basis11 * tDiff * points[index + 1].tangent +
                d2_basis01 * points[index + 1].position
                ) / (tDiff * tDiff);
    }

    inline InterpolationType computeWiggle(size_t index, floating_t tDiff) const
    {
        //tests and such have shown that we have to scale this by the inverse of the t distance, and i'm not sure why
        //intuitively it would just be the 2nd derivative of the position function and nothing else
        //if you know why please let me know
        return (
                    floating_t(12) * (points[index].position - points[index + 1].position) + floating_t(6) * tDiff * (points[index].tangent + points[index + 1].tangent)
                ) / (tDiff * tDiff * tDiff);
    }


    inline floating_t computeSegmentLength(size_t index, floating_t from, floating_t to) const
    {
        floating_t tDiff = knots[index + 1] - knots[index];
        auto segmentFunction = [this, index, tDiff](floating_t t) -> floating_t {
            auto tangent = computeTangent(index, tDiff, t);
            return tangent.length();
        };

        return tDiff * SplineLibraryCalculus::gaussLegendreQuadratureIntegral(segmentFunction, from, to);
    }

private: //data
    std::vector<CubicHermiteSplinePoint> points;
    std::vector<floating_t> knots;
};
